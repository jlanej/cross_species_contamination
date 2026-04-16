#!/usr/bin/env bash
# =============================================================================
# classify_array.sh
#
# SLURM array job: run Kraken2 taxonomic classification on per-sample FASTQ
# files previously extracted by extract_unmapped_array.sh.
#
# Runs entirely inside an Apptainer (Singularity) container so that no
# local Kraken2 or Python installation is required on the HPC node.
# The container image is pulled automatically on first use.
#
# Each array task processes one sample listed in the classify manifest TSV
# produced by submit_classify.sh.
#
# Output per sample:
#   <CLASSIFY_OUTDIR>/<SAMPLE_ID>/<SAMPLE_ID>.kraken2.report.txt
#   <CLASSIFY_OUTDIR>/<SAMPLE_ID>/<SAMPLE_ID>.kraken2.output.txt
#
# Usage (see submit_classify.sh for the recommended submission wrapper):
#
#   sbatch --array=1-N classify_array.sh
#
# Required environment variables (exported by submit_classify.sh or set manually):
#   CLASSIFY_MANIFEST – TSV with columns SAMPLE_ID, R1, R2
#                       (R2 empty → single-end; non-empty → paired-end)
#   CLASSIFY_OUTDIR   – output directory for per-sample classify results
#   EXTRACT_OUTDIR    – directory containing the extracted FASTQ files
#                       (needed for Apptainer bind mount)
#   DB                – path to the Kraken2 database directory
#
# Optional environment variables:
#   THREADS           – Kraken2 threads (default: SLURM_CPUS_PER_TASK or 4)
#   CONFIDENCE        – Kraken2 confidence threshold (default: 0.0)
#   MEMORY_MAPPING    – set to "1" to use memory-mapped DB (lower RAM, slower)
#   CONTAINER_SIF     – path to the Apptainer SIF image.
#                       Defaults to <CLASSIFY_OUTDIR>/csc.sif; auto-pulled from
#                       CONTAINER_IMAGE if absent.
#   CONTAINER_IMAGE   – Docker URI for auto-pull
#                       (default: ghcr.io/jlanej/cross_species_contamination:latest)
#
# AI assistance acknowledgment: developed with AI assistance.
# =============================================================================

#SBATCH --job-name=csc_classify
#SBATCH --output=logs/classify_%A_%a.out
#SBATCH --error=logs/classify_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=normal

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --------------------------------------------------------------------------- #
# Configurable defaults                                                         #
# --------------------------------------------------------------------------- #
CLASSIFY_MANIFEST="${CLASSIFY_MANIFEST:-}"
CLASSIFY_OUTDIR="${CLASSIFY_OUTDIR:-}"
EXTRACT_OUTDIR="${EXTRACT_OUTDIR:-}"
DB="${DB:-}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
CONFIDENCE="${CONFIDENCE:-0.0}"
MEMORY_MAPPING="${MEMORY_MAPPING:-0}"

CONTAINER_IMAGE="${CONTAINER_IMAGE:-ghcr.io/jlanej/cross_species_contamination:latest}"
# Default to CLASSIFY_OUTDIR so that direct sbatch invocations use a writable
# location.  submit_classify.sh always passes this via --export.
CONTAINER_SIF="${CONTAINER_SIF:-${CLASSIFY_OUTDIR}/csc.sif}"

# --------------------------------------------------------------------------- #
# Apptainer bootstrap (identical pattern to extract_unmapped_array.sh)         #
# --------------------------------------------------------------------------- #

if command -v apptainer &>/dev/null; then
    APPTAINER_CMD="apptainer"
elif command -v singularity &>/dev/null; then
    APPTAINER_CMD="singularity"
else
    echo "ERROR: Neither 'apptainer' nor 'singularity' found in PATH." >&2
    echo "       Load the apptainer module (e.g. 'module load apptainer') or" >&2
    echo "       install it: https://apptainer.org/docs/user/latest/quick_start.html" >&2
    exit 1
fi

pull_container() {
    [[ -f "${CONTAINER_SIF}" ]] && return 0

    echo "Container SIF not found: ${CONTAINER_SIF}"
    echo "Pulling ${CONTAINER_IMAGE} – this only happens once per cluster..."

    mkdir -p "$(dirname "${CONTAINER_SIF}")"

    local lock="${CONTAINER_SIF}.lock"
    local max_wait=300
    local waited=0

    if [[ -f "${lock}" ]]; then
        echo "Another task is pulling the image; waiting up to ${max_wait}s..."
        while [[ -f "${lock}" && ${waited} -lt ${max_wait} ]]; do
            sleep 5; waited=$(( waited + 5 ))
        done
        [[ -f "${CONTAINER_SIF}" ]] && { echo "Container now available."; return 0; }
        echo "WARNING: Timed out waiting for pull; retrying." >&2
    fi

    touch "${lock}"
    if "${APPTAINER_CMD}" pull --force "${CONTAINER_SIF}" "docker://${CONTAINER_IMAGE}"; then
        rm -f "${lock}"
        echo "Container pulled successfully: ${CONTAINER_SIF}"
    else
        rm -f "${lock}"
        echo "ERROR: Failed to pull container image '${CONTAINER_IMAGE}'." >&2
        exit 1
    fi
}

# Run a command inside the container, binding every path that the tool will
# need to read or write.
container_run() {
    local -a bind_args=()
    bind_args+=("--bind" "${CLASSIFY_OUTDIR}:${CLASSIFY_OUTDIR}")
    bind_args+=("--bind" "${EXTRACT_OUTDIR}:${EXTRACT_OUTDIR}")
    bind_args+=("--bind" "${DB}:${DB}")
    "${APPTAINER_CMD}" exec "${bind_args[@]}" "${CONTAINER_SIF}" "$@"
}

# --------------------------------------------------------------------------- #
# Validate environment                                                          #
# --------------------------------------------------------------------------- #
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Submit with --array." >&2
    exit 1
fi

if [[ -z "${CLASSIFY_MANIFEST}" ]]; then
    echo "ERROR: CLASSIFY_MANIFEST is not set." >&2
    exit 1
fi

if [[ ! -f "${CLASSIFY_MANIFEST}" ]]; then
    echo "ERROR: Classify manifest not found: ${CLASSIFY_MANIFEST}" >&2
    exit 1
fi

if [[ -z "${DB}" ]]; then
    echo "ERROR: DB is not set. Export DB=/path/to/kraken2_database" >&2
    exit 1
fi

if [[ ! -d "${DB}" ]]; then
    echo "ERROR: Kraken2 database directory not found: ${DB}" >&2
    exit 1
fi

if [[ -z "${CLASSIFY_OUTDIR}" ]]; then
    echo "ERROR: CLASSIFY_OUTDIR is not set." >&2
    exit 1
fi

if [[ -z "${EXTRACT_OUTDIR}" ]]; then
    echo "ERROR: EXTRACT_OUTDIR is not set." >&2
    exit 1
fi

# Pull/verify the container (no-op if SIF already exists)
pull_container

# --------------------------------------------------------------------------- #
# Load the sample line from the classify manifest                               #
# Header is line 1; sample data starts at line 2.                              #
# SLURM_ARRAY_TASK_ID=1 → manifest line 2, etc.                               #
# --------------------------------------------------------------------------- #
LINE_NUM=$(( SLURM_ARRAY_TASK_ID + 1 ))
MANIFEST_LINE="$(sed -n "${LINE_NUM}p" "${CLASSIFY_MANIFEST}")"

if [[ -z "${MANIFEST_LINE}" ]]; then
    echo "ERROR: No data at line ${LINE_NUM} in ${CLASSIFY_MANIFEST} (TASK_ID=${SLURM_ARRAY_TASK_ID})" >&2
    exit 1
fi

# Parse TSV columns
SAMPLE_ID="$(printf '%s' "${MANIFEST_LINE}" | cut -f1)"
R1="$(printf '%s' "${MANIFEST_LINE}"        | cut -f2)"
R2="$(printf '%s' "${MANIFEST_LINE}"        | cut -f3)"

if [[ -z "${SAMPLE_ID}" || -z "${R1}" ]]; then
    echo "ERROR: Could not parse manifest line ${LINE_NUM}: ${MANIFEST_LINE}" >&2
    exit 1
fi

echo "======================================================"
echo "  Sample        : ${SAMPLE_ID}"
echo "  R1             : ${R1}"
echo "  R2             : ${R2:-<single-end>}"
echo "  DB             : ${DB}"
echo "  Classify outdir: ${CLASSIFY_OUTDIR}"
echo "  Threads        : ${THREADS}"
echo "  Confidence     : ${CONFIDENCE}"
echo "  Memory mapping : ${MEMORY_MAPPING}"
echo "  Container      : ${CONTAINER_SIF}"
echo "======================================================"

# --------------------------------------------------------------------------- #
# Prepare per-sample output directory                                           #
# --------------------------------------------------------------------------- #
SAMPLE_CLASSIFY_DIR="${CLASSIFY_OUTDIR}/${SAMPLE_ID}"
mkdir -p "${SAMPLE_CLASSIFY_DIR}"

# Idempotence: skip if classification report already exists and is non-empty
REPORT="${SAMPLE_CLASSIFY_DIR}/${SAMPLE_ID}.kraken2.report.txt"
if [[ -s "${REPORT}" ]]; then
    echo "Classification already complete for ${SAMPLE_ID}. Skipping."
    exit 0
fi

# Validate input FASTQ files
if [[ ! -s "${R1}" ]]; then
    echo "ERROR: R1 FASTQ not found or empty: ${R1}" >&2
    exit 1
fi
if [[ -n "${R2}" && ! -s "${R2}" ]]; then
    echo "ERROR: R2 FASTQ not found or empty: ${R2}" >&2
    exit 1
fi

# --------------------------------------------------------------------------- #
# Build csc-classify command                                                    #
# --------------------------------------------------------------------------- #
CLASSIFY_ARGS=(
    csc-classify
    "--db"         "${DB}"
    "-o"           "${SAMPLE_CLASSIFY_DIR}"
    "--sample-id"  "${SAMPLE_ID}"
    "--confidence" "${CONFIDENCE}"
    "--threads"    "${THREADS}"
)

if [[ -n "${R2}" ]]; then
    # Paired-end: prepend both FASTQs and add --paired
    CLASSIFY_ARGS=("${CLASSIFY_ARGS[0]}" "${R1}" "${R2}" "${CLASSIFY_ARGS[@]:1}" "--paired")
else
    # Single-end
    CLASSIFY_ARGS=("${CLASSIFY_ARGS[0]}" "${R1}" "${CLASSIFY_ARGS[@]:1}")
fi

if [[ "${MEMORY_MAPPING}" == "1" ]]; then
    CLASSIFY_ARGS+=("--memory-mapping")
fi

# --------------------------------------------------------------------------- #
# Run classification inside the container                                       #
# --------------------------------------------------------------------------- #
echo "Running classification for ${SAMPLE_ID}..."
EXIT_CODE=0
container_run "${CLASSIFY_ARGS[@]}" || EXIT_CODE=$?

if [[ ${EXIT_CODE} -ne 0 ]]; then
    echo "ERROR: Classification failed for ${SAMPLE_ID} (exit ${EXIT_CODE})" >&2
    # Remove partial outputs so the job is not considered done on restart
    rm -f "${SAMPLE_CLASSIFY_DIR}/${SAMPLE_ID}.kraken2"*.txt
    exit ${EXIT_CODE}
fi

echo "Done: ${SAMPLE_ID}"
echo "Output files:"
ls -lh "${SAMPLE_CLASSIFY_DIR}/${SAMPLE_ID}.kraken2"*.txt 2>/dev/null \
    || echo "  (no output files found)"
