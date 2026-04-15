#!/usr/bin/env bash
# =============================================================================
# extract_unmapped_array.sh
#
# SLURM array job: extract unmapped reads from 1000 Genomes CRAMs.
#
# Runs entirely inside an Apptainer (Singularity) container so that no
# local tool installation is required – just Apptainer on the HPC node.
# The container image is pulled automatically on first run.
#
# Each array task processes one CRAM listed in the manifest. Only the
# unmapped-read section of the remote CRAM is fetched (via the CRAI index),
# so the full file is never downloaded.
#
# Output per sample:
#   <OUTDIR>/<SAMPLE_ID>/<SAMPLE_ID>_unmapped_R1.fastq.gz
#   <OUTDIR>/<SAMPLE_ID>/<SAMPLE_ID>_unmapped_R2.fastq.gz   (if paired)
#   <OUTDIR>/<SAMPLE_ID>/<SAMPLE_ID>_unmapped_singleton.fastq.gz (unpaired)
#
# Usage (see submit_extract.sh for a convenience wrapper):
#
#   # All 3202 samples
#   sbatch --array=1-3202 extract_unmapped_array.sh
#
#   # First 10 samples
#   sbatch --array=1-10 extract_unmapped_array.sh
#
#   # Specific samples by 1-based line number (comma/hyphen list)
#   sbatch --array=1,5,42-50 extract_unmapped_array.sh
#
# Required environment variables (set defaults below or export before sbatch):
#   MANIFEST      – path to manifest.tsv (default: this script's directory)
#   OUTDIR        – output directory     (default: ./output)
#
# Optional environment variables:
#   CONTAINER_SIF – path to the Apptainer SIF image.
#                   Defaults to <SCRIPT_DIR>/csc.sif; auto-pulled from
#                   CONTAINER_IMAGE if absent.
#   CONTAINER_IMAGE – Docker URI for auto-pull
#                   (default: ghcr.io/jlanej/cross_species_contamination:latest)
#   REFERENCE     – path/URL to GRCh38 reference FASTA for CRAM decoding.
#                   If absent, samtools uses lossy decoding for unmapped reads,
#                   which is usually fine but may emit a warning.
#   THREADS       – samtools threads (default: matches --cpus-per-task)
#   KEEP_CRAM     – if set to "1", also save the raw unmapped CRAM before FASTQ
#                   conversion (useful for downstream csc-extract / re-processing)
#
# AI assistance acknowledgment: developed with AI assistance.
# =============================================================================

#SBATCH --job-name=1kg_extract_unmapped
#SBATCH --output=logs/extract_%A_%a.out
#SBATCH --error=logs/extract_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=normal

set -euo pipefail

# --------------------------------------------------------------------------- #
# Configurable defaults                                                         #
# --------------------------------------------------------------------------- #
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${MANIFEST:-${SCRIPT_DIR}/manifest.tsv}"
OUTDIR="${OUTDIR:-${SCRIPT_DIR}/output}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
KEEP_CRAM="${KEEP_CRAM:-0}"

# Container settings – no pre-setup required; image is pulled automatically
CONTAINER_IMAGE="${CONTAINER_IMAGE:-ghcr.io/jlanej/cross_species_contamination:latest}"
CONTAINER_SIF="${CONTAINER_SIF:-${SCRIPT_DIR}/csc.sif}"

# --------------------------------------------------------------------------- #
# Apptainer bootstrap                                                           #
# --------------------------------------------------------------------------- #

# Detect the apptainer/singularity executable
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

# Pull the SIF image on first use (batteries-included – no pre-setup required).
# A lock file prevents concurrent pulls from racing across array tasks.
pull_container() {
    [[ -f "${CONTAINER_SIF}" ]] && return 0

    echo "Container SIF not found: ${CONTAINER_SIF}"
    echo "Pulling ${CONTAINER_IMAGE} – this only happens once per cluster..."

    mkdir -p "$(dirname "${CONTAINER_SIF}")"

    local lock="${CONTAINER_SIF}.lock"
    local max_wait=300
    local waited=0

    # If another task is already pulling, wait for it to finish
    if [[ -f "${lock}" ]]; then
        echo "Another task is pulling the image; waiting up to ${max_wait}s..."
        while [[ -f "${lock}" && ${waited} -lt ${max_wait} ]]; do
            sleep 5; waited=$(( waited + 5 ))
        done
        [[ -f "${CONTAINER_SIF}" ]] && { echo "Container now available."; return 0; }
        echo "WARNING: Timed out waiting for pull; retrying." >&2
    fi

    # Claim the lock and pull
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

# Run a command inside the container, binding the output directory and (if
# provided) the directory containing the reference FASTA.
container_run() {
    local -a bind_args=("--bind" "${OUTDIR}:${OUTDIR}")
    if [[ -n "${REFERENCE:-}" ]]; then
        local ref_dir
        ref_dir="$(dirname "${REFERENCE}")"
        bind_args+=("--bind" "${ref_dir}:${ref_dir}")
    fi
    "${APPTAINER_CMD}" exec "${bind_args[@]}" "${CONTAINER_SIF}" "$@"
}

# --------------------------------------------------------------------------- #
# Validate environment                                                          #
# --------------------------------------------------------------------------- #
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Submit with --array." >&2
    exit 1
fi

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}" >&2
    exit 1
fi

# Pull/verify the container (no-op if SIF already exists)
pull_container

# --------------------------------------------------------------------------- #
# Load the sample line from the manifest                                        #
# The manifest has a header on line 1; sample data starts at line 2.           #
# SLURM_ARRAY_TASK_ID=1 → first sample (manifest line 2), etc.                #
# --------------------------------------------------------------------------- #
LINE_NUM=$(( SLURM_ARRAY_TASK_ID + 1 ))   # +1 to skip header
MANIFEST_LINE="$(sed -n "${LINE_NUM}p" "${MANIFEST}")"

if [[ -z "${MANIFEST_LINE}" ]]; then
    echo "ERROR: No data at line ${LINE_NUM} in ${MANIFEST} (TASK_ID=${SLURM_ARRAY_TASK_ID})" >&2
    exit 1
fi

# Parse TSV columns
SAMPLE_ID="$(echo "${MANIFEST_LINE}" | cut -f1)"
CRAM_URL="$(echo "${MANIFEST_LINE}"  | cut -f2)"
CRAI_URL="$(echo "${MANIFEST_LINE}"  | cut -f3)"

if [[ -z "${SAMPLE_ID}" || -z "${CRAM_URL}" || -z "${CRAI_URL}" ]]; then
    echo "ERROR: Could not parse manifest line ${LINE_NUM}: ${MANIFEST_LINE}" >&2
    exit 1
fi

echo "======================================================"
echo "  Sample     : ${SAMPLE_ID}"
echo "  CRAM URL   : ${CRAM_URL}"
echo "  CRAI URL   : ${CRAI_URL}"
echo "  Output dir : ${OUTDIR}"
echo "  Threads    : ${THREADS}"
echo "  Container  : ${CONTAINER_SIF}"
echo "======================================================"

# --------------------------------------------------------------------------- #
# Prepare output directory                                                      #
# --------------------------------------------------------------------------- #
SAMPLE_DIR="${OUTDIR}/${SAMPLE_ID}"
mkdir -p "${SAMPLE_DIR}"

# Skip if outputs already exist (restart-safe)
R1="${SAMPLE_DIR}/${SAMPLE_ID}_unmapped_R1.fastq.gz"
if [[ -s "${R1}" ]]; then
    echo "Outputs already present for ${SAMPLE_ID}. Skipping."
    exit 0
fi

# --------------------------------------------------------------------------- #
# Build samtools view arguments                                                 #
# We fetch ONLY the '*' (unmapped) virtual contig from the remote CRAM,        #
# using the remote CRAI so samtools can seek directly to that section.         #
# --------------------------------------------------------------------------- #
VIEW_ARGS=(
    samtools view
    --threads "${THREADS}"
    -u          # uncompressed BAM on stdout (piped to fastq)
    -f 4        # FLAG: read unmapped
    -X "${CRAI_URL}"
    "${CRAM_URL}"
    '*'
)

if [[ -n "${REFERENCE:-}" ]]; then
    VIEW_ARGS=(
        samtools view
        -T "${REFERENCE}"
        --threads "${THREADS}"
        -u -f 4
        -X "${CRAI_URL}"
        "${CRAM_URL}"
        '*'
    )
fi

# --------------------------------------------------------------------------- #
# Optionally save intermediate unmapped CRAM                                   #
# --------------------------------------------------------------------------- #
COLLATE_TMP="${SAMPLE_DIR}/.collate_tmp_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${COLLATE_TMP}"

if [[ "${KEEP_CRAM}" == "1" ]]; then
    UNMAPPED_CRAM="${SAMPLE_DIR}/${SAMPLE_ID}_unmapped.cram"
    echo "Saving intermediate unmapped CRAM: ${UNMAPPED_CRAM}"
    container_run "${VIEW_ARGS[@]}" -C -o "${UNMAPPED_CRAM}"
fi

# --------------------------------------------------------------------------- #
# Stream unmapped reads → collate by name → convert to FASTQ                  #
# --------------------------------------------------------------------------- #
echo "Extracting unmapped reads and converting to FASTQ..."

if [[ "${KEEP_CRAM}" == "1" ]]; then
    container_run \
        samtools collate \
            --threads "${THREADS}" \
            -u -O \
            "${SAMPLE_DIR}/${SAMPLE_ID}_unmapped.cram" \
            "${COLLATE_TMP}/tmp" \
    | container_run \
        samtools fastq \
            --threads "${THREADS}" \
            -1 "${R1}" \
            -2 "${SAMPLE_DIR}/${SAMPLE_ID}_unmapped_R2.fastq.gz" \
            -s "${SAMPLE_DIR}/${SAMPLE_ID}_unmapped_singleton.fastq.gz" \
            -0 "${SAMPLE_DIR}/${SAMPLE_ID}_unmapped_other.fastq.gz" \
            -
else
    # Single pipeline: view → collate → fastq (all tools run in the container)
    container_run "${VIEW_ARGS[@]}" \
    | container_run \
        samtools collate \
            --threads "${THREADS}" \
            -u -O \
            - \
            "${COLLATE_TMP}/tmp" \
    | container_run \
        samtools fastq \
            --threads "${THREADS}" \
            -1 "${R1}" \
            -2 "${SAMPLE_DIR}/${SAMPLE_ID}_unmapped_R2.fastq.gz" \
            -s "${SAMPLE_DIR}/${SAMPLE_ID}_unmapped_singleton.fastq.gz" \
            -0 "${SAMPLE_DIR}/${SAMPLE_ID}_unmapped_other.fastq.gz" \
            -
fi

EXIT_CODE=$?

# Clean up temporary collate files
rm -rf "${COLLATE_TMP}"

if [[ ${EXIT_CODE} -ne 0 ]]; then
    echo "ERROR: samtools pipeline failed for ${SAMPLE_ID} (exit ${EXIT_CODE})" >&2
    # Remove partial outputs so the job is not considered done on restart
    rm -f "${SAMPLE_DIR}/${SAMPLE_ID}_unmapped"*.fastq.gz
    exit ${EXIT_CODE}
fi

# Remove empty output files (e.g. no singletons)
find "${SAMPLE_DIR}" -name "${SAMPLE_ID}_unmapped*.fastq.gz" -empty -delete

echo "Done: ${SAMPLE_ID}"
echo "Output files:"
ls -lh "${SAMPLE_DIR}/${SAMPLE_ID}_unmapped"*.fastq.gz 2>/dev/null \
    || echo "  (no output files – sample may have no unmapped reads)"
