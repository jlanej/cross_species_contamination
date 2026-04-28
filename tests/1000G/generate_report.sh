#!/usr/bin/env bash
# =============================================================================
# generate_report.sh
#
# SLURM job: generate a static, self-contained HTML contamination report from
# the outputs of aggregate_detect.sh (csc-aggregate + csc-detect).
#
# This script is designed to run after aggregate_detect.sh completes.
# submit_classify.sh submits it with --dependency=afterok:<agg_job_id>
# so it starts automatically when aggregation and detection have finished.
#
# Output:
#   <REPORT_OUTDIR>/contamination_report.html  – self-contained HTML report
#   <REPORT_OUTDIR>/report_manifest.json       – machine-readable sidecar
#
# Usage (see submit_classify.sh for the recommended submission wrapper):
#
#   sbatch generate_report.sh
#
# Required environment variables:
#   AGG_OUTDIR    – directory containing csc-aggregate outputs
#                   (taxa_matrix_raw.tsv, taxa_matrix_cpm.tsv, etc.)
#
# Optional environment variables:
#   DETECT_OUTDIR – directory containing csc-detect outputs
#                   (flagged_samples.tsv, qc_summary.json); when set, outlier
#                   detection parameters and flagged samples are included in the
#                   report (default: <AGG_OUTDIR>/../detect if that directory
#                   exists, otherwise omitted)
#   REPORT_OUTDIR – output directory for the HTML report
#                   (default: <AGG_OUTDIR>/../report)
#   REPORT_FILE   – full path for the output HTML file
#                   (default: <REPORT_OUTDIR>/contamination_report.html)
#   REPORT_TITLE  – title string shown in the HTML <h1> and <title>
#                   (default: "Cross-Species Contamination Summary Report")
#   TOP_N         – number of taxa to show in the cohort-wide top table
#                   (default: 10)
#   VARIANT_IMPACT_THRESHOLD_PPM
#                 – absolute-burden threshold (reads per million total
#                   sequenced reads) above which samples are flagged in the
#                   Variant-Calling Impact section (default: 1000 ppm = 0.1%)
#   SKIP_DETECT_IN_REPORT
#                 – set to "1" to omit detect outputs from the report even if
#                   DETECT_OUTDIR is set (default: 0)
#   CONTAINER_SIF   – path to the Apptainer SIF image
#   CONTAINER_IMAGE – Docker URI for auto-pull
#                     (default: ghcr.io/jlanej/cross_species_contamination:latest)
#
# AI assistance acknowledgment: developed with AI assistance.
# =============================================================================

#SBATCH --job-name=csc_report
#SBATCH --output=logs/generate_report_%j.out
#SBATCH --error=logs/generate_report_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --time=04:30:00
#SBATCH --partition=normal

set -euo pipefail

# --------------------------------------------------------------------------- #
# Configurable defaults                                                         #
# --------------------------------------------------------------------------- #
AGG_OUTDIR="${AGG_OUTDIR:-}"
DETECT_OUTDIR="${DETECT_OUTDIR:-}"
REPORT_OUTDIR="${REPORT_OUTDIR:-}"
REPORT_FILE="${REPORT_FILE:-}"
REPORT_TITLE="${REPORT_TITLE:-Cross-Species Contamination Summary Report}"
TOP_N="${TOP_N:-10}"
VARIANT_IMPACT_THRESHOLD_PPM="${VARIANT_IMPACT_THRESHOLD_PPM:-}"
SKIP_DETECT_IN_REPORT="${SKIP_DETECT_IN_REPORT:-0}"

CONTAINER_IMAGE="${CONTAINER_IMAGE:-ghcr.io/jlanej/cross_species_contamination:latest}"
CONTAINER_SIF="${CONTAINER_SIF:-}"

# --------------------------------------------------------------------------- #
# Validate required variables                                                   #
# --------------------------------------------------------------------------- #
if [[ -z "${AGG_OUTDIR}" ]]; then
    echo "ERROR: AGG_OUTDIR is not set." >&2
    exit 1
fi

if [[ ! -d "${AGG_OUTDIR}" ]]; then
    echo "ERROR: AGG_OUTDIR not found: ${AGG_OUTDIR}" >&2
    exit 1
fi

# Resolve defaults that depend on AGG_OUTDIR
if [[ -z "${REPORT_OUTDIR}" ]]; then
    REPORT_OUTDIR="$(dirname "${AGG_OUTDIR}")/report"
fi
if [[ -z "${REPORT_FILE}" ]]; then
    REPORT_FILE="${REPORT_OUTDIR}/contamination_report.html"
fi
if [[ -z "${CONTAINER_SIF}" ]]; then
    CONTAINER_SIF="$(dirname "${AGG_OUTDIR}")/csc.sif"
fi

# Auto-discover DETECT_OUTDIR when not explicitly set
if [[ -z "${DETECT_OUTDIR}" ]]; then
    candidate="$(dirname "${AGG_OUTDIR}")/detect"
    if [[ -d "${candidate}" ]]; then
        DETECT_OUTDIR="${candidate}"
    fi
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
    bind_args+=("--bind" "${AGG_OUTDIR}:${AGG_OUTDIR}")
    bind_args+=("--bind" "${REPORT_OUTDIR}:${REPORT_OUTDIR}")
    if [[ "${SKIP_DETECT_IN_REPORT}" != "1" && -n "${DETECT_OUTDIR}" && -d "${DETECT_OUTDIR}" ]]; then
        bind_args+=("--bind" "${DETECT_OUTDIR}:${DETECT_OUTDIR}")
    fi
    "${APPTAINER_CMD}" exec "${bind_args[@]}" "${CONTAINER_SIF}" "$@"
}

# Pull/verify the container
pull_container

# --------------------------------------------------------------------------- #
# Print configuration summary                                                   #
# --------------------------------------------------------------------------- #
echo "======================================================"
echo "  Aggregate outdir : ${AGG_OUTDIR}"
if [[ "${SKIP_DETECT_IN_REPORT}" != "1" && -n "${DETECT_OUTDIR}" && -d "${DETECT_OUTDIR}" ]]; then
    echo "  Detect outdir    : ${DETECT_OUTDIR}"
else
    echo "  Detect outdir    : (not included in report)"
fi
echo "  Report outdir    : ${REPORT_OUTDIR}"
echo "  Report file      : ${REPORT_FILE}"
echo "  Report title     : ${REPORT_TITLE}"
echo "  Top-N taxa       : ${TOP_N}"
if [[ -n "${VARIANT_IMPACT_THRESHOLD_PPM}" ]]; then
    echo "  VI threshold ppm : ${VARIANT_IMPACT_THRESHOLD_PPM}"
else
    echo "  VI threshold ppm : (default: 1000 ppm)"
fi
echo "  Container        : ${CONTAINER_SIF}"
echo "======================================================"

# --------------------------------------------------------------------------- #
# Generate HTML report                                                          #
# --------------------------------------------------------------------------- #
mkdir -p "${REPORT_OUTDIR}"

REPORT_ARGS=(
    csc-report
    "${AGG_OUTDIR}"
    -o "${REPORT_FILE}"
    --title "${REPORT_TITLE}"
    --top-n "${TOP_N}"
)

if [[ "${SKIP_DETECT_IN_REPORT}" != "1" && -n "${DETECT_OUTDIR}" && -d "${DETECT_OUTDIR}" ]]; then
    REPORT_ARGS+=("--detect-dir" "${DETECT_OUTDIR}")
fi

if [[ -n "${VARIANT_IMPACT_THRESHOLD_PPM}" ]]; then
    REPORT_ARGS+=("--variant-impact-threshold-ppm" "${VARIANT_IMPACT_THRESHOLD_PPM}")
fi

echo ""
echo "=== Generating HTML report ==="
EXIT_CODE=0
container_run "${REPORT_ARGS[@]}" || EXIT_CODE=$?

if [[ ${EXIT_CODE} -ne 0 ]]; then
    echo "ERROR: csc-report failed (exit ${EXIT_CODE})" >&2
    exit ${EXIT_CODE}
fi

echo ""
echo "=== Report complete ==="
echo "  HTML report : ${REPORT_FILE}"
echo "  Manifest    : ${REPORT_OUTDIR}/report_manifest.json"
