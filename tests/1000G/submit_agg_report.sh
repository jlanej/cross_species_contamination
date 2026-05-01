#!/usr/bin/env bash
# =============================================================================
# submit_agg_report.sh
#
# Submission wrapper for the CSC aggregation, detection, and report pipeline.
#
# This script submits:
#   1. aggregate_detect.sh – aggregates per-sample Kraken2 reports, runs
#      statistical outlier detection (csc-aggregate + csc-detect).
#   2. generate_report.sh  – produces a self-contained HTML contamination
#      report from the aggregate and detect outputs (csc-report).
#
# Designed as a standalone submission that runs after classification is
# complete.  Use --dependency <JOB_ID> to chain it automatically after a
# classify array job (e.g. the job submitted by submit_classify.sh).
#
# Quick start:
#
# Usage:
#   ./submit_agg_report.sh [options]
#
#   # Run agg/detect/report on existing classify outputs
#   ./submit_agg_report.sh \
#       --outdir         /scratch/me/1kg_classify \
#       --extract-outdir /scratch/me/1kg_out \
#       --db             /data/kraken2/PrackenDB
#
#   # Chain after a classify job (depends on classify array)
#   ./submit_agg_report.sh \
#       --outdir         /scratch/me/1kg_classify \
#       --extract-outdir /scratch/me/1kg_out \
#       --db             /data/kraken2/PrackenDB \
#       --dependency     <CLASSIFY_JOB_ID>
#
#   # Dry-run: print sbatch commands without submitting
#   ./submit_agg_report.sh \
#       --outdir         /scratch/me/1kg_classify \
#       --extract-outdir /scratch/me/1kg_out \
#       --db             /data/kraken2/PrackenDB \
#       --dry-run
#
#   # Skip detect (aggregate only)
#   ./submit_agg_report.sh \
#       --outdir         /scratch/me/1kg_classify \
#       --extract-outdir /scratch/me/1kg_out \
#       --db             /data/kraken2/PrackenDB \
#       --skip-detect
#
# Options:
#   --outdir        DIR    Output base directory (same as submit_classify.sh
#                          --outdir); classify/, aggregate/, detect/, and
#                          report/ subdirectories are derived from here
#                          [default: ./classify_output]
#   --extract-outdir DIR   Extraction output directory containing per-sample
#                          `{sample}.reads_summary.json` idxstats sidecars.
#                          Required unless --skip-idxstats-metrics.
#                          [default: ./output]
#   --db            DIR    Path to Kraken2 database directory (contains
#                          taxonomy/nodes.dmp).  Used for lineage-aware domain
#                          annotation and confidence-threshold aggregation.
#                          Defaults to --db-path if set; may be omitted when
#                          neither feature is needed.
#   --partition     STR    SLURM partition              [default: normal]
#   --min-reads     N      Min direct reads per taxon in aggregate [default: 0]
#   --rank-filter   STR    Space-separated rank codes for per-rank matrices
#                          [default: "S G F"]
#   --db-path       DIR    Kraken2 DB dir for lineage-aware domain annotation;
#                          defaults to the value of --db (set to "" to disable)
#   --confidence-thresholds STR  Colon-separated Kraken2 confidence cutoffs
#                          (e.g. "0.1:0.5") for the high-confidence tier.
#                          Each value > 0 produces a parallel matrix set with
#                          suffix _conf{T} using per-read kraken2 outputs
#                          already produced by classify.  Requires --db-path.
#                          [default: "0.1" – dual-tier (sensitive + 0.1);
#                          set to "" to disable]
#   --detect-matrix STR    Matrix type for detect input: cpm or raw [default: cpm]
#   --detect-method STR    Outlier detection method: all, mad, iqr, or gmm
#                          [default: all]
#   --mad-threshold FLOAT  MAD threshold for outlier detection [default: 3.5]
#   --iqr-multiplier FLOAT IQR multiplier for outlier detection [default: 1.5]
#   --gmm-threshold FLOAT  GMM posterior probability threshold [default: 0.5]
#   --skip-detect          Run aggregate but skip the detect step
#   --no-abs-detection     Disable the absolute-burden side pass that
#                          csc-detect runs by default when an
#                          absolute-burden sibling matrix is available
#                          (see docs/detect.md)
#   --skip-idxstats-metrics  Skip idxstats-based absolute burden metrics in
#                          aggregate/reporting (default: off)
#   --agg-cpus      N      CPUs for aggregate/detect job [default: 4]
#   --agg-mem       STR    Memory for aggregate/detect job [default: 16G]
#   --agg-time      STR    Wall-clock time for aggregate/detect [default: 02:00:00]
#   --skip-report          Skip HTML report generation (default: off)
#   --report-title  STR    Title shown in the HTML report
#                          [default: "Cross-Species Contamination Summary Report"]
#   --report-top-n  N      Number of top taxa shown in the report [default: 10]
#   --variant-impact-threshold-ppm FLOAT
#                          Absolute-burden threshold (ppm) for the Variant-Calling
#                          Impact section of the report [default: 1000]
#   --report-cpus   N      CPUs for the report job [default: 2]
#   --report-mem    STR    Memory for the report job [default: 8G]
#   --report-time   STR    Wall-clock time for the report job [default: 00:30:00]
#   --dependency    JOB_ID Optional SLURM job dependency (afterok:<JOB_ID>);
#                          use to chain this script after a classify array job
#   --container     FILE   Path to the Apptainer SIF image
#                          [default: <outdir>/csc.sif; auto-pulled if absent]
#   --image         URI    Docker URI to pull the container from
#                          [default: ghcr.io/jlanej/cross_species_contamination:latest]
#   --dry-run              Print sbatch commands; do not submit
#   -h, --help             Show this help message
#
# AI assistance acknowledgment: developed with AI assistance.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Defaults ─────────────────────────────────────────────────────────────────
OUTDIR="${SCRIPT_DIR}/classify_output"
EXTRACT_OUTDIR="${SCRIPT_DIR}/output"
DB=""
PARTITION="normal"
MIN_READS=0
RANK_FILTER="S G F"
DB_PATH=""
CONFIDENCE_THRESHOLDS="0.1"  # colon-separated; "" = sensitive tier only.
                             # Default 0.1 enables dual-tier reporting
                             # (sensitive 0.0 + high-confidence 0.1) per
                             # Wood et al. 2019 / Marcelino et al. 2020.
DETECT_MATRIX="cpm"
DETECT_METHOD="all"
MAD_THRESHOLD=3.5
IQR_MULTIPLIER=1.5
GMM_THRESHOLD=0.5
SKIP_DETECT=0
NO_ABS_DETECTION=0
SKIP_IDXSTATS_METRICS=0
AGG_CPUS=4
AGG_MEM="16G"
AGG_WALLTIME="02:00:00"
SKIP_REPORT=0
REPORT_TITLE="Cross-Species Contamination Summary Report"
REPORT_TOP_N=10
VARIANT_IMPACT_THRESHOLD_PPM=""
REPORT_CPUS=2
REPORT_MEM="8G"
REPORT_WALLTIME="00:30:00"
DEPENDENCY=""
CONTAINER_SIF=""          # resolved later to an absolute path under OUTDIR
CONTAINER_IMAGE="ghcr.io/jlanej/cross_species_contamination:latest"
DRY_RUN=0
# Keep usage output focused on the documented header/options block.
USAGE_LINES=95

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | sed 's/^# \{0,2\}//' | head -"${USAGE_LINES}"
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)         OUTDIR="$2";         shift 2 ;;
        --extract-outdir) EXTRACT_OUTDIR="$2"; shift 2 ;;
        --db)             DB="$2";             shift 2 ;;
        --partition)      PARTITION="$2";      shift 2 ;;
        --min-reads)      MIN_READS="$2";      shift 2 ;;
        --rank-filter)    RANK_FILTER="$2";    shift 2 ;;
        --db-path)        DB_PATH="$2";        shift 2 ;;
        --confidence-thresholds) CONFIDENCE_THRESHOLDS="$2"; shift 2 ;;
        --detect-matrix)  DETECT_MATRIX="$2";  shift 2 ;;
        --detect-method)  DETECT_METHOD="$2";  shift 2 ;;
        --mad-threshold)  MAD_THRESHOLD="$2";  shift 2 ;;
        --iqr-multiplier) IQR_MULTIPLIER="$2"; shift 2 ;;
        --gmm-threshold)  GMM_THRESHOLD="$2";  shift 2 ;;
        --skip-detect)    SKIP_DETECT=1;       shift ;;
        --no-abs-detection) NO_ABS_DETECTION=1; shift ;;
        --skip-idxstats-metrics) SKIP_IDXSTATS_METRICS=1; shift ;;
        --agg-cpus)       AGG_CPUS="$2";       shift 2 ;;
        --agg-mem)        AGG_MEM="$2";        shift 2 ;;
        --agg-time)       AGG_WALLTIME="$2";   shift 2 ;;
        --skip-report)    SKIP_REPORT=1;       shift ;;
        --report-title)   REPORT_TITLE="$2";   shift 2 ;;
        --report-top-n)   REPORT_TOP_N="$2";   shift 2 ;;
        --variant-impact-threshold-ppm) VARIANT_IMPACT_THRESHOLD_PPM="$2"; shift 2 ;;
        --report-cpus)    REPORT_CPUS="$2";    shift 2 ;;
        --report-mem)     REPORT_MEM="$2";     shift 2 ;;
        --report-time)    REPORT_WALLTIME="$2"; shift 2 ;;
        --dependency)     DEPENDENCY="$2";     shift 2 ;;
        --container)      CONTAINER_SIF="$2";  shift 2 ;;
        --image)          CONTAINER_IMAGE="$2"; shift 2 ;;
        --dry-run)        DRY_RUN=1;           shift ;;
        -h|--help)        usage ;;
        *) echo "ERROR: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ── Validate required arguments ───────────────────────────────────────────────
if [[ "${DETECT_METHOD}" != "all" && "${DETECT_METHOD}" != "mad" && "${DETECT_METHOD}" != "iqr" && "${DETECT_METHOD}" != "gmm" ]]; then
    echo "ERROR: --detect-method must be 'all', 'mad', 'iqr', or 'gmm'." >&2
    exit 1
fi
if [[ "${DETECT_MATRIX}" != "cpm" && "${DETECT_MATRIX}" != "raw" ]]; then
    echo "ERROR: --detect-matrix must be 'cpm' or 'raw'." >&2
    exit 1
fi

if [[ "${DRY_RUN}" -eq 0 ]] && [[ ! -d "${OUTDIR}" ]]; then
    echo "ERROR: Output directory not found: ${OUTDIR}" >&2
    exit 1
fi

if [[ "${SKIP_IDXSTATS_METRICS}" -eq 0 ]]; then
    if [[ ! -d "${EXTRACT_OUTDIR}" ]]; then
        echo "ERROR: Extract output directory not found: ${EXTRACT_OUTDIR}" >&2
        echo "       Set --skip-idxstats-metrics to bypass idxstats-based metrics." >&2
        exit 1
    fi
fi

# ── Derived paths ─────────────────────────────────────────────────────────────
CLASSIFY_OUTDIR="${OUTDIR}/classify"
AGG_OUTDIR="${OUTDIR}/aggregate"
DETECT_OUTDIR="${OUTDIR}/detect"
REPORT_OUTDIR="${OUTDIR}/report"
LOG_DIR="${OUTDIR}/logs"

# Resolve container SIF to an absolute, user-writable path so SLURM worker
# nodes (which copy the script to /var/spool/slurmd/…) never need to pull.
if [[ -z "${CONTAINER_SIF}" ]]; then
    CONTAINER_SIF="${OUTDIR}/csc.sif"
fi
CONTAINER_SIF="$(realpath -m "${CONTAINER_SIF}")"

# Encode rank filter as colon-separated (no spaces) for safe export to jobs
RANK_FILTER_CODES="${RANK_FILTER// /:}"

# Default DB_PATH to the Kraken2 DB (which contains taxonomy/nodes.dmp) so
# that csc-aggregate performs lineage-aware domain annotation automatically.
# Pass --db-path "" explicitly to disable domain annotation.
if [[ -z "${DB_PATH}" ]] && [[ -n "${DB}" ]]; then
    DB_PATH="${DB}"
fi

# ── Print configuration summary ───────────────────────────────────────────────
echo "======================================================"
echo "  Outdir          : ${OUTDIR}"
echo "  Classify outdir : ${CLASSIFY_OUTDIR}"
echo "  Aggregate outdir: ${AGG_OUTDIR}"
echo "  Detect outdir   : ${DETECT_OUTDIR}"
echo "  Report outdir   : ${REPORT_OUTDIR}"
if [[ "${SKIP_IDXSTATS_METRICS}" -eq 0 ]]; then
    echo "  Extract outdir  : ${EXTRACT_OUTDIR}"
else
    echo "  Extract outdir  : (idxstats metrics skipped)"
fi
if [[ -n "${DB_PATH}" ]]; then
    echo "  DB path         : ${DB_PATH} (domain annotation enabled)"
else
    echo "  DB path         : (not set – domain annotation disabled)"
fi
echo "  Min reads       : ${MIN_READS}"
echo "  Detect matrix   : ${DETECT_MATRIX}"
echo "  Rank filter     : ${RANK_FILTER}"
echo "  Skip detect     : ${SKIP_DETECT}"
if [[ "${SKIP_DETECT}" -ne 1 ]]; then
    echo "  Detect method   : ${DETECT_METHOD}"
    echo "  MAD threshold   : ${MAD_THRESHOLD}"
    echo "  IQR multiplier  : ${IQR_MULTIPLIER}"
    echo "  GMM threshold   : ${GMM_THRESHOLD}"
fi
if [[ -n "${CONFIDENCE_THRESHOLDS}" ]]; then
    echo "  Confidence tiers: ${CONFIDENCE_THRESHOLDS//:/ } (high-confidence aggregation enabled)"
else
    echo "  Confidence tiers: (sensitive tier only)"
fi
echo "  Agg CPUs        : ${AGG_CPUS}"
echo "  Agg mem         : ${AGG_MEM}"
echo "  Agg time        : ${AGG_WALLTIME}"
if [[ -n "${DEPENDENCY}" ]]; then
    echo "  Dependency      : afterok:${DEPENDENCY}"
fi
echo "  Container       : ${CONTAINER_SIF}"
echo "======================================================"

# ── Pull container once before any jobs start ─────────────────────────────────
if [[ "${DRY_RUN}" -eq 0 ]] && [[ ! -f "${CONTAINER_SIF}" ]]; then
    if command -v apptainer &>/dev/null; then
        APPTAINER_CMD="apptainer"
    elif command -v singularity &>/dev/null; then
        APPTAINER_CMD="singularity"
    else
        echo "ERROR: Neither 'apptainer' nor 'singularity' found in PATH." >&2
        exit 1
    fi
    echo "Pulling container image: ${CONTAINER_IMAGE}"
    echo "Saving to: ${CONTAINER_SIF}"
    mkdir -p "$(dirname "${CONTAINER_SIF}")"
    if "${APPTAINER_CMD}" pull --force "${CONTAINER_SIF}" "docker://${CONTAINER_IMAGE}"; then
        echo "Container pulled successfully: ${CONTAINER_SIF}"
    else
        echo "ERROR: Failed to pull container image '${CONTAINER_IMAGE}'." >&2
        exit 1
    fi
elif [[ "${DRY_RUN}" -eq 0 ]]; then
    echo "Container already present: ${CONTAINER_SIF}"
fi

# ── Submit aggregate/detect job ───────────────────────────────────────────────
mkdir -p "${AGG_OUTDIR}" "${DETECT_OUTDIR}" "${LOG_DIR}"

AGG_EXPORTS=(
    "CLASSIFY_OUTDIR=${CLASSIFY_OUTDIR}"
    "EXTRACT_OUTDIR=${EXTRACT_OUTDIR}"
    "AGG_OUTDIR=${AGG_OUTDIR}"
    "DETECT_OUTDIR=${DETECT_OUTDIR}"
    "THREADS=${AGG_CPUS}"
    "MIN_READS=${MIN_READS}"
    "RANK_FILTER_CODES=${RANK_FILTER_CODES}"
    "DETECT_MATRIX=${DETECT_MATRIX}"
    "DETECT_METHOD=${DETECT_METHOD}"
    "MAD_THRESHOLD=${MAD_THRESHOLD}"
    "IQR_MULTIPLIER=${IQR_MULTIPLIER}"
    "GMM_THRESHOLD=${GMM_THRESHOLD}"
    "SKIP_DETECT=${SKIP_DETECT}"
    "NO_ABS_DETECTION=${NO_ABS_DETECTION}"
    "SKIP_IDXSTATS_METRICS=${SKIP_IDXSTATS_METRICS}"
    "DB_PATH=${DB_PATH}"
    "CONFIDENCE_THRESHOLDS=${CONFIDENCE_THRESHOLDS}"
    "CONTAINER_SIF=${CONTAINER_SIF}"
    "CONTAINER_IMAGE=${CONTAINER_IMAGE}"
)
AGG_EXPORT_STR="ALL,$(IFS=,; echo "${AGG_EXPORTS[*]}")"

SBATCH_AGG=(
    sbatch
    --job-name=csc_aggregate_detect
    --cpus-per-task="${AGG_CPUS}"
    --mem="${AGG_MEM}"
    --time="${AGG_WALLTIME}"
    --partition="${PARTITION}"
    --output="${LOG_DIR}/aggregate_detect_%j.out"
    --error="${LOG_DIR}/aggregate_detect_%j.err"
    --export="${AGG_EXPORT_STR}"
    "${SCRIPT_DIR}/aggregate_detect.sh"
)

# Add dependency when requested (e.g. afterok:<classify_job_id>)
if [[ -n "${DEPENDENCY}" ]]; then
    SBATCH_AGG=("${SBATCH_AGG[0]}" "--dependency=afterok:${DEPENDENCY}" "${SBATCH_AGG[@]:1}")
fi

echo ""
echo "Aggregate/detect sbatch command:"
printf '  %s \\\n' "${SBATCH_AGG[@]}"
echo ""

AGG_JOB_ID=""
if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "(Dry-run: aggregate/detect job not submitted)"
    AGG_JOB_ID="<dry-run>"
else
    AGG_JOB_OUTPUT="$("${SBATCH_AGG[@]}")"
    echo "${AGG_JOB_OUTPUT}"
    AGG_JOB_ID="$(echo "${AGG_JOB_OUTPUT}" | grep -oE '[0-9]+$')"
    echo "Aggregate/detect job submitted: ${AGG_JOB_ID}"
fi

# ── Submit report job (with dependency on aggregate/detect) ──────────────────
if [[ "${SKIP_REPORT}" -eq 1 ]]; then
    echo "Skipping report generation (--skip-report)."
    if [[ "${DRY_RUN}" -eq 0 ]] && [[ -n "${AGG_JOB_ID}" ]]; then
        echo "Monitor with: squeue -j ${AGG_JOB_ID}"
    fi
    exit 0
fi

mkdir -p "${REPORT_OUTDIR}"

REPORT_EXPORTS=(
    "AGG_OUTDIR=${AGG_OUTDIR}"
    "DETECT_OUTDIR=${DETECT_OUTDIR}"
    "REPORT_OUTDIR=${REPORT_OUTDIR}"
    "REPORT_TITLE=${REPORT_TITLE}"
    "TOP_N=${REPORT_TOP_N}"
    "SKIP_DETECT_IN_REPORT=${SKIP_DETECT}"
    "CONTAINER_SIF=${CONTAINER_SIF}"
    "CONTAINER_IMAGE=${CONTAINER_IMAGE}"
)
if [[ -n "${VARIANT_IMPACT_THRESHOLD_PPM}" ]]; then
    REPORT_EXPORTS+=("VARIANT_IMPACT_THRESHOLD_PPM=${VARIANT_IMPACT_THRESHOLD_PPM}")
fi
REPORT_EXPORT_STR="ALL,$(IFS=,; echo "${REPORT_EXPORTS[*]}")"

SBATCH_REPORT=(
    sbatch
    --job-name=csc_report
    --cpus-per-task="${REPORT_CPUS}"
    --mem="${REPORT_MEM}"
    --time="${REPORT_WALLTIME}"
    --partition="${PARTITION}"
    --output="${LOG_DIR}/generate_report_%j.out"
    --error="${LOG_DIR}/generate_report_%j.err"
    --export="${REPORT_EXPORT_STR}"
    "${SCRIPT_DIR}/generate_report.sh"
)

# Add dependency on aggregate/detect job
if [[ -n "${AGG_JOB_ID}" && "${AGG_JOB_ID}" != "<dry-run>" ]]; then
    SBATCH_REPORT=("${SBATCH_REPORT[0]}" "--dependency=afterok:${AGG_JOB_ID}" "${SBATCH_REPORT[@]:1}")
elif [[ "${DRY_RUN}" -eq 1 && -n "${AGG_JOB_ID}" ]]; then
    SBATCH_REPORT=("${SBATCH_REPORT[0]}" "--dependency=afterok:${AGG_JOB_ID}" "${SBATCH_REPORT[@]:1}")
fi

echo ""
echo "Report sbatch command:"
printf '  %s \\\n' "${SBATCH_REPORT[@]}"
echo ""

if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "(Dry-run: report job not submitted)"
else
    REPORT_JOB_OUTPUT="$("${SBATCH_REPORT[@]}")"
    echo "${REPORT_JOB_OUTPUT}"
    REPORT_JOB_ID="$(echo "${REPORT_JOB_OUTPUT}" | grep -oE '[0-9]+$')"
    echo "Report job submitted: ${REPORT_JOB_ID}"
    echo ""
    echo "Pipeline submitted:"
    echo "  Aggregate/detect:  job ${AGG_JOB_ID}"
    echo "  Report:            job ${REPORT_JOB_ID}  [depends on ${AGG_JOB_ID}]"
    echo ""
    echo "Monitor with: squeue -j ${AGG_JOB_ID},${REPORT_JOB_ID}"
fi
