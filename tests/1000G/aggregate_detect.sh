#!/usr/bin/env bash
# =============================================================================
# aggregate_detect.sh
#
# SLURM job: aggregate all per-sample Kraken2 reports into a sample-by-taxon
# matrix (csc-aggregate), then run statistical outlier detection on that matrix
# (csc-detect) to flag cross-species contamination.
#
# This script is designed to run after classify_array.sh completes.
# submit_classify.sh submits it with --dependency=afterok:<classify_job_id>
# so it starts automatically when all classification array tasks have finished.
#
# Output:
#   <AGG_OUTDIR>/taxa_matrix.tsv          – primary compatibility matrix (CPM by default)
#   <AGG_OUTDIR>/taxa_matrix_cpm.tsv      – all-sample CPM matrix
#   <AGG_OUTDIR>/taxa_matrix_raw.tsv      – all-sample raw-count matrix
#   <AGG_OUTDIR>/taxa_matrix_S.tsv        – species-level matrix
#   <AGG_OUTDIR>/taxa_matrix_G.tsv        – genus-level matrix
#   <AGG_OUTDIR>/taxa_matrix_F.tsv        – family-level matrix
#   <AGG_OUTDIR>/aggregation_metadata.json
#   <DETECT_OUTDIR>/flagged_samples.tsv
#   <DETECT_OUTDIR>/qc_summary.json
#   <DETECT_OUTDIR>/quarantine_list.txt
#   <DETECT_OUTDIR>/{S,G,F}/…             – per-rank detect outputs
#
# Usage (see submit_classify.sh for the recommended submission wrapper):
#
#   sbatch aggregate_detect.sh
#
# Required environment variables:
#   CLASSIFY_OUTDIR – directory containing per-sample classify outputs
#   AGG_OUTDIR      – output directory for aggregate results
#
# Optional environment variables:
#   DETECT_OUTDIR   – output directory for detect results
#                     (default: <AGG_OUTDIR>/../detect)
#   THREADS         – CPUs for aggregate/detect (default: SLURM_CPUS_PER_TASK or 4)
#   MIN_READS       – minimum direct-read count per taxon (default: 0)
#   NO_NORMALIZE    – set to "1" to output raw counts instead of CPM (default: 0)
#   RANK_FILTER_CODES – colon-separated rank codes for per-rank matrices
#                       (default: S:G:F)
#   DETECT_MATRIX   – matrix type for csc-detect input: cpm or raw (default: cpm)
#   DETECT_METHOD   – outlier detection method: mad or iqr (default: mad)
#   MAD_THRESHOLD   – MAD threshold (default: 3.5)
#   IQR_MULTIPLIER  – IQR multiplier (default: 1.5)
#   SKIP_DETECT     – set to "1" to skip outlier detection (default: 0)
#   CONTAINER_SIF   – path to the Apptainer SIF image
#   CONTAINER_IMAGE – Docker URI for auto-pull
#                     (default: ghcr.io/jlanej/cross_species_contamination:latest)
#
# AI assistance acknowledgment: developed with AI assistance.
# =============================================================================

#SBATCH --job-name=csc_aggregate_detect
#SBATCH --output=logs/aggregate_detect_%j.out
#SBATCH --error=logs/aggregate_detect_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=normal

set -euo pipefail

# --------------------------------------------------------------------------- #
# Configurable defaults                                                         #
# --------------------------------------------------------------------------- #
CLASSIFY_OUTDIR="${CLASSIFY_OUTDIR:-}"
AGG_OUTDIR="${AGG_OUTDIR:-}"
# DETECT_OUTDIR defaults to a sibling directory of AGG_OUTDIR; computed below
# after we know AGG_OUTDIR.
DETECT_OUTDIR="${DETECT_OUTDIR:-}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
MIN_READS="${MIN_READS:-0}"
NO_NORMALIZE="${NO_NORMALIZE:-0}"
# Colon-separated rank codes to avoid spaces in the --export string
RANK_FILTER_CODES="${RANK_FILTER_CODES:-S:G:F}"
DETECT_MATRIX="${DETECT_MATRIX:-cpm}"
DETECT_METHOD="${DETECT_METHOD:-mad}"
MAD_THRESHOLD="${MAD_THRESHOLD:-3.5}"
IQR_MULTIPLIER="${IQR_MULTIPLIER:-1.5}"
SKIP_DETECT="${SKIP_DETECT:-0}"

CONTAINER_IMAGE="${CONTAINER_IMAGE:-ghcr.io/jlanej/cross_species_contamination:latest}"
CONTAINER_SIF="${CONTAINER_SIF:-${AGG_OUTDIR}/csc.sif}"

# --------------------------------------------------------------------------- #
# Validate required variables                                                   #
# --------------------------------------------------------------------------- #
if [[ -z "${CLASSIFY_OUTDIR}" ]]; then
    echo "ERROR: CLASSIFY_OUTDIR is not set." >&2
    exit 1
fi

if [[ ! -d "${CLASSIFY_OUTDIR}" ]]; then
    echo "ERROR: CLASSIFY_OUTDIR not found: ${CLASSIFY_OUTDIR}" >&2
    exit 1
fi

if [[ -z "${AGG_OUTDIR}" ]]; then
    echo "ERROR: AGG_OUTDIR is not set." >&2
    exit 1
fi
if [[ "${DETECT_MATRIX}" != "cpm" && "${DETECT_MATRIX}" != "raw" ]]; then
    echo "ERROR: DETECT_MATRIX must be 'cpm' or 'raw'." >&2
    exit 1
fi

# Default DETECT_OUTDIR now that AGG_OUTDIR is validated
if [[ -z "${DETECT_OUTDIR}" ]]; then
    DETECT_OUTDIR="$(dirname "${AGG_OUTDIR}")/detect"
fi

# --------------------------------------------------------------------------- #
# Apptainer bootstrap                                                           #
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
        echo "Another process is pulling the image; waiting up to ${max_wait}s..."
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

container_run() {
    local -a bind_args=()
    bind_args+=("--bind" "${CLASSIFY_OUTDIR}:${CLASSIFY_OUTDIR}")
    bind_args+=("--bind" "${AGG_OUTDIR}:${AGG_OUTDIR}")
    if [[ "${SKIP_DETECT}" != "1" ]]; then
        bind_args+=("--bind" "${DETECT_OUTDIR}:${DETECT_OUTDIR}")
    fi
    "${APPTAINER_CMD}" exec "${bind_args[@]}" "${CONTAINER_SIF}" "$@"
}

# Pull/verify the container
pull_container

# --------------------------------------------------------------------------- #
# Discover classification reports                                               #
# --------------------------------------------------------------------------- #
mapfile -t REPORTS < <(find "${CLASSIFY_OUTDIR}" -name "*.kraken2.report.txt" | sort)

if [[ ${#REPORTS[@]} -eq 0 ]]; then
    echo "ERROR: No .kraken2.report.txt files found under ${CLASSIFY_OUTDIR}" >&2
    exit 1
fi

echo "Found ${#REPORTS[@]} classification report(s) to aggregate."

# Convert colon-separated rank codes to space-separated for display/logging
RANK_FILTER_DISPLAY="${RANK_FILTER_CODES//:/ }"
IFS=':' read -ra RANK_CODES <<< "${RANK_FILTER_CODES}"

echo "======================================================"
echo "  Classify outdir : ${CLASSIFY_OUTDIR}"
echo "  Aggregate outdir: ${AGG_OUTDIR}"
echo "  Detect outdir   : ${DETECT_OUTDIR}"
echo "  Reports found   : ${#REPORTS[@]}"
echo "  Min reads       : ${MIN_READS}"
echo "  Normalize       : $([[ "${NO_NORMALIZE}" == "1" ]] && echo "no (raw counts)" || echo "yes (CPM)")"
echo "  Detect matrix   : ${DETECT_MATRIX}"
echo "  Rank filter     : ${RANK_FILTER_DISPLAY}"
echo "  Skip detect     : ${SKIP_DETECT}"
if [[ "${SKIP_DETECT}" != "1" ]]; then
    echo "  Detect method   : ${DETECT_METHOD}"
    echo "  MAD threshold   : ${MAD_THRESHOLD}"
    echo "  IQR multiplier  : ${IQR_MULTIPLIER}"
fi
echo "  Container       : ${CONTAINER_SIF}"
echo "======================================================"

# --------------------------------------------------------------------------- #
# Step 1: Aggregate reports → taxa matrix                                       #
# --------------------------------------------------------------------------- #
mkdir -p "${AGG_OUTDIR}"

AGGREGATE_ARGS=(
    csc-aggregate
    "${REPORTS[@]}"
    -o "${AGG_OUTDIR}"
)
if [[ "${MIN_READS}" -gt 0 ]]; then
    AGGREGATE_ARGS+=("--min-reads" "${MIN_READS}")
fi
if [[ "${NO_NORMALIZE}" == "1" ]]; then
    AGGREGATE_ARGS+=("--no-normalize")
fi
AGGREGATE_ARGS+=("--rank-filter" "${RANK_CODES[@]}")

echo ""
echo "=== Step 1: Aggregating ${#REPORTS[@]} reports ==="
EXIT_CODE=0
container_run "${AGGREGATE_ARGS[@]}" || EXIT_CODE=$?

if [[ ${EXIT_CODE} -ne 0 ]]; then
    echo "ERROR: csc-aggregate failed (exit ${EXIT_CODE})" >&2
    exit ${EXIT_CODE}
fi

echo "Aggregation complete. Matrix: ${AGG_OUTDIR}/taxa_matrix.tsv"

# --------------------------------------------------------------------------- #
# Step 2: Detect outliers (optional)                                            #
# --------------------------------------------------------------------------- #
if [[ "${SKIP_DETECT}" == "1" ]]; then
    echo "Skipping outlier detection (SKIP_DETECT=1)."
    exit 0
fi

MATRIX="${AGG_OUTDIR}/taxa_matrix_${DETECT_MATRIX}.tsv"
if [[ ! -s "${MATRIX}" ]]; then
    echo "ERROR: Expected aggregate matrix not found: ${MATRIX}" >&2
    exit 1
fi

mkdir -p "${DETECT_OUTDIR}"

DETECT_ARGS=(
    csc-detect
    "${MATRIX}"
    -o "${DETECT_OUTDIR}"
    --method       "${DETECT_METHOD}"
    --mad-threshold "${MAD_THRESHOLD}"
    --iqr-multiplier "${IQR_MULTIPLIER}"
    --rank-filter  "${RANK_CODES[@]}"
)

echo ""
echo "=== Step 2: Detecting outliers ==="
EXIT_CODE=0
container_run "${DETECT_ARGS[@]}" || EXIT_CODE=$?

if [[ ${EXIT_CODE} -ne 0 ]]; then
    echo "ERROR: csc-detect failed (exit ${EXIT_CODE})" >&2
    exit ${EXIT_CODE}
fi

echo "Detection complete. Results: ${DETECT_OUTDIR}"
echo ""
echo "=== Pipeline complete ==="
echo "  Aggregate : ${AGG_OUTDIR}"
echo "  Detect    : ${DETECT_OUTDIR}"
