#!/usr/bin/env bash
# =============================================================================
# submit_classify.sh
#
# Submission wrapper for the CSC classification pipeline.
#
# This script:
#   1. Scans the extraction output directory (from extract_unmapped_array.sh)
#      and consolidates the list of samples that have extracted FASTQ files.
#   2. Submits classify_array.sh as a SLURM array job (one task per sample).
#
# After classification completes, run submit_agg_report.sh to aggregate,
# detect, and generate the HTML report.  Use its --dependency option to chain
# it automatically after this job (pass the job ID printed here).
#
# Quick start:
#   # Classify all extracted samples
#   ./submit_classify.sh \
#       --extract-outdir /scratch/me/1kg_out \
#       --outdir         /scratch/me/1kg_classify \
#       --db             /data/kraken2/PrackenDB
#
#   # Dry-run: print sbatch commands without submitting
#   ./submit_classify.sh \
#       --extract-outdir /scratch/me/1kg_out \
#       --outdir         /scratch/me/1kg_classify \
#       --db             /data/kraken2/PrackenDB \
#       --dry-run
#
#   # Only a subset of samples
#   ./submit_classify.sh \
#       --extract-outdir /scratch/me/1kg_out \
#       --outdir         /scratch/me/1kg_classify \
#       --db             /data/kraken2/PrackenDB \
#       --samples        my_subset.txt
#
# Options:
#   --extract-outdir DIR   Directory with FASTQ files from extract_unmapped_array.sh
#                          [default: ./output]
#   --outdir        DIR    Output base directory for classify outputs
#                          [default: ./classify_output]
#   --db            DIR    Path to Kraken2 database directory (required)
#   --samples       FILE   Text file with SAMPLE_IDs to process (one per line)
#   --limit         N      Process only the first N extracted samples
#   --range         STR    Explicit SLURM --array range string (e.g. "1-50,60")
#   --partition     STR    SLURM partition              [default: normal]
#   --cpus          N      CPUs per classify task       [default: 4]
#   --mem           STR    Memory per classify task     [default: 16G]
#   --time          STR    Wall-clock time limit        [default: 04:00:00]
#   --max-concurrent-jobs N  Max concurrent classify array tasks [default: 200]
#   --confidence    FLOAT  Kraken2 confidence threshold [default: 0.0]
#   --memory-mapping       Use memory-mapped Kraken2 DB (lower RAM, slower)
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
EXTRACT_OUTDIR="${SCRIPT_DIR}/output"
OUTDIR="${SCRIPT_DIR}/classify_output"
DB=""
SAMPLES_FILE=""
LIMIT=""
RANGE=""
PARTITION="normal"
CPUS=4
MEM="16G"
WALLTIME="04:00:00"
MAX_CONCURRENT_JOBS=200
MAX_CONCURRENT_JOBS_SET=0
CONFIDENCE=0.0
MEMORY_MAPPING=0
CONTAINER_SIF=""          # resolved later to an absolute path under OUTDIR
CONTAINER_IMAGE="ghcr.io/jlanej/cross_species_contamination:latest"
DRY_RUN=0
# Keep usage output focused on the documented header/options block.
USAGE_LINES=65

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | sed 's/^# \{0,2\}//' | head -"${USAGE_LINES}"
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --extract-outdir) EXTRACT_OUTDIR="$2"; shift 2 ;;
        --outdir)         OUTDIR="$2";         shift 2 ;;
        --db)             DB="$2";             shift 2 ;;
        --samples)        SAMPLES_FILE="$2";   shift 2 ;;
        --limit)          LIMIT="$2";          shift 2 ;;
        --range)          RANGE="$2";          shift 2 ;;
        --partition)      PARTITION="$2";      shift 2 ;;
        --cpus)           CPUS="$2";           shift 2 ;;
        --mem)            MEM="$2";            shift 2 ;;
        --time)           WALLTIME="$2";       shift 2 ;;
        --max-concurrent-jobs) MAX_CONCURRENT_JOBS="$2"; MAX_CONCURRENT_JOBS_SET=1; shift 2 ;;
        --confidence)     CONFIDENCE="$2";     shift 2 ;;
        --memory-mapping) MEMORY_MAPPING=1;    shift ;;
        --container)      CONTAINER_SIF="$2";  shift 2 ;;
        --image)          CONTAINER_IMAGE="$2"; shift 2 ;;
        --dry-run)        DRY_RUN=1;           shift ;;
        -h|--help)        usage ;;
        *) echo "ERROR: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ── Validate required arguments ───────────────────────────────────────────────
if [[ -z "${DB}" ]]; then
    echo "ERROR: --db is required. Supply the path to a Kraken2 database directory." >&2
    exit 1
fi

if [[ "${DRY_RUN}" -eq 0 ]] && [[ ! -d "${DB}" ]]; then
    echo "ERROR: Kraken2 database not found: ${DB}" >&2
    exit 1
fi

if ! [[ "${MAX_CONCURRENT_JOBS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --max-concurrent-jobs must be a positive integer." >&2
    exit 1
fi

if [[ ! -d "${EXTRACT_OUTDIR}" ]]; then
    echo "ERROR: Extract output directory not found: ${EXTRACT_OUTDIR}" >&2
    exit 1
fi

# ── Derived paths ─────────────────────────────────────────────────────────────
CLASSIFY_OUTDIR="${OUTDIR}/classify"
CLASSIFY_MANIFEST="${OUTDIR}/classify_manifest.tsv"
LOG_DIR="${OUTDIR}/logs"

# Resolve container SIF to an absolute, user-writable path so SLURM worker
# nodes (which copy the script to /var/spool/slurmd/…) never need to pull.
if [[ -z "${CONTAINER_SIF}" ]]; then
    CONTAINER_SIF="${OUTDIR}/csc.sif"
fi
CONTAINER_SIF="$(realpath -m "${CONTAINER_SIF}")"

# ── Scan extraction output for samples with FASTQ files ──────────────────────
echo "Scanning ${EXTRACT_OUTDIR} for extracted samples..."

declare -a ALL_SAMPLE_IDS=()
declare -a ALL_R1_PATHS=()
declare -a ALL_R2_PATHS=()

while IFS= read -r sample_dir; do
    [[ -d "${sample_dir}" ]] || continue
    sid="$(basename "${sample_dir}")"
    r1="${sample_dir}/${sid}_unmapped_R1.fastq.gz"
    r2="${sample_dir}/${sid}_unmapped_R2.fastq.gz"
    if [[ -s "${r1}" ]]; then
        ALL_SAMPLE_IDS+=("${sid}")
        ALL_R1_PATHS+=("${r1}")
        r2_val=""
        [[ -s "${r2}" ]] && r2_val="${r2}"
        ALL_R2_PATHS+=("${r2_val}")
    fi
done < <(find "${EXTRACT_OUTDIR}" -mindepth 1 -maxdepth 1 -type d | sort)

if [[ ${#ALL_SAMPLE_IDS[@]} -eq 0 ]]; then
    echo "ERROR: No extracted FASTQ files found under ${EXTRACT_OUTDIR}" >&2
    echo "       Run extract_unmapped_array.sh (via submit_extract.sh) first." >&2
    exit 1
fi

echo "Found ${#ALL_SAMPLE_IDS[@]} extracted sample(s) in ${EXTRACT_OUTDIR}"

# ── Apply --samples / --limit filter ─────────────────────────────────────────
declare -a SELECTED_IDS=()
declare -a SELECTED_R1=()
declare -a SELECTED_R2=()

if [[ -n "${SAMPLES_FILE}" ]]; then
    if [[ ! -f "${SAMPLES_FILE}" ]]; then
        echo "ERROR: Samples file not found: ${SAMPLES_FILE}" >&2
        exit 1
    fi

    echo "Filtering to samples listed in ${SAMPLES_FILE}..."
    # Build a quick lookup of ALL_SAMPLE_IDS → index
    declare -A sid_to_idx
    for i in "${!ALL_SAMPLE_IDS[@]}"; do
        sid_to_idx["${ALL_SAMPLE_IDS[$i]}"]="${i}"
    done

    while IFS= read -r sid; do
        [[ -z "${sid}" ]] && continue
        if [[ -v "sid_to_idx[${sid}]" ]]; then
            idx="${sid_to_idx[${sid}]}"
            SELECTED_IDS+=("${ALL_SAMPLE_IDS[$idx]}")
            SELECTED_R1+=("${ALL_R1_PATHS[$idx]}")
            SELECTED_R2+=("${ALL_R2_PATHS[$idx]}")
        else
            echo "WARNING: Sample '${sid}' has no extracted FASTQs in ${EXTRACT_OUTDIR}; skipping." >&2
        fi
    done < "${SAMPLES_FILE}"

    unset sid_to_idx

    if [[ ${#SELECTED_IDS[@]} -eq 0 ]]; then
        echo "ERROR: No samples from ${SAMPLES_FILE} had extracted FASTQ files." >&2
        exit 1
    fi
    echo "Selected ${#SELECTED_IDS[@]} sample(s) from samples file."

else
    SELECTED_IDS=("${ALL_SAMPLE_IDS[@]}")
    SELECTED_R1=("${ALL_R1_PATHS[@]}")
    SELECTED_R2=("${ALL_R2_PATHS[@]}")
fi

# Apply --limit (take first N)
if [[ -n "${LIMIT}" ]]; then
    if ! [[ "${LIMIT}" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: --limit must be a positive integer." >&2
        exit 1
    fi
    if (( LIMIT < ${#SELECTED_IDS[@]} )); then
        SELECTED_IDS=("${SELECTED_IDS[@]:0:${LIMIT}}")
        SELECTED_R1=("${SELECTED_R1[@]:0:${LIMIT}}")
        SELECTED_R2=("${SELECTED_R2[@]:0:${LIMIT}}")
        echo "Limiting to first ${LIMIT} sample(s)."
    fi
fi

TOTAL_SELECTED=${#SELECTED_IDS[@]}
echo "Processing ${TOTAL_SELECTED} sample(s)."

# ── Write classify manifest ───────────────────────────────────────────────────
mkdir -p "${OUTDIR}" "${LOG_DIR}"

{
    printf 'SAMPLE_ID\tR1\tR2\n'
    for i in "${!SELECTED_IDS[@]}"; do
        printf '%s\t%s\t%s\n' \
            "${SELECTED_IDS[$i]}" \
            "${SELECTED_R1[$i]}" \
            "${SELECTED_R2[$i]}"
    done
} > "${CLASSIFY_MANIFEST}"

echo "Classify manifest written: ${CLASSIFY_MANIFEST} (${TOTAL_SELECTED} samples)"

# ── Helper: expand a SLURM array spec into individual indices ─────────────────
expand_array_spec() {
    local spec="$1"
    local token start end i
    local -a tokens=()
    IFS=',' read -ra tokens <<< "${spec}"
    for token in "${tokens[@]}"; do
        [[ -z "${token}" ]] && continue
        if [[ "${token}" =~ ^([0-9]+)-([0-9]+)$ ]]; then
            start="${BASH_REMATCH[1]}"
            end="${BASH_REMATCH[2]}"
            if (( start > end )); then
                echo "ERROR: Invalid range '${token}' (start > end)." >&2
                return 1
            fi
            for (( i=start; i<=end; i++ )); do
                echo "${i}"
            done
        elif [[ "${token}" =~ ^[0-9]+$ ]]; then
            echo "${token}"
        else
            echo "ERROR: Invalid array token '${token}'." >&2
            return 1
        fi
    done
}

# Convert a sorted list of integers (one per line on stdin) to compact SLURM
# range notation, merging consecutive runs into N-M ranges.
# Example: 1 2 3 5 6 10  →  1-3,5-6,10
compress_indices_to_range() {
    awk '
    BEGIN { start = -1; prev = -1 }
    {
        n = $1 + 0
        if (start < 0) { start = n; prev = n; next }
        if (n == prev + 1) { prev = n; next }
        if (prev > start) printf "%d-%d,", start, prev
        else printf "%d,", start
        start = n; prev = n
    }
    END {
        if (start >= 0) {
            if (prev > start) printf "%d-%d", start, prev
            else printf "%d", start
        }
    }'
}

# Check if a sample (by 1-based manifest index) is already classified
sample_classified_for_index() {
    local idx="$1"
    local sid report
    sid="${SELECTED_IDS[$((idx - 1))]:-}"
    [[ -z "${sid}" ]] && return 1
    report="${CLASSIFY_OUTDIR}/${sid}/${sid}.kraken2.report.txt"
    [[ -s "${report}" ]]
}

# ── Resolve array range ───────────────────────────────────────────────────────
if [[ -n "${RANGE}" ]]; then
    ARRAY_SPEC="${RANGE}"
else
    ARRAY_SPEC="1-${TOTAL_SELECTED}"
fi

# ── Idempotence: skip already-classified samples ──────────────────────────────
COMPLETED_COUNT=0
PENDING_INDICES=()

while IFS= read -r idx; do
    [[ -z "${idx}" ]] && continue
    # Validate index is within the manifest bounds
    if (( idx < 1 || idx > TOTAL_SELECTED )); then
        echo "WARNING: Index ${idx} is out of manifest bounds (1-${TOTAL_SELECTED}); skipping." >&2
        continue
    fi
    if sample_classified_for_index "${idx}"; then
        COMPLETED_COUNT=$(( COMPLETED_COUNT + 1 ))
    else
        PENDING_INDICES+=("${idx}")
    fi
done < <(expand_array_spec "${ARRAY_SPEC}")

if [[ ${#PENDING_INDICES[@]} -eq 0 ]]; then
    echo "All selected samples already have classification output. Nothing to classify."
    exit 0
else
    if [[ "${COMPLETED_COUNT}" -gt 0 ]]; then
        echo "Skipping ${COMPLETED_COUNT} already-classified sample(s); submitting ${#PENDING_INDICES[@]} pending sample(s)."
    fi

    # Rebuild the array spec from pending indices, compressing consecutive runs
    PENDING_ARRAY_SPEC="$(printf '%s\n' "${PENDING_INDICES[@]}" | sort -n | uniq | compress_indices_to_range)"

    # Apply concurrency throttle
    if [[ "${PENDING_ARRAY_SPEC}" == *%* ]]; then
        if [[ "${MAX_CONCURRENT_JOBS_SET}" -eq 1 ]]; then
            echo "ERROR: Cannot use --max-concurrent-jobs when --range already specifies '%' concurrency." >&2
            exit 1
        fi
    else
        PENDING_ARRAY_SPEC="${PENDING_ARRAY_SPEC}%${MAX_CONCURRENT_JOBS}"
    fi

    CLASSIFY_JOB_ID=""   # filled in after submission
fi

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

# ── Build classify exports string ─────────────────────────────────────────────
CLASSIFY_EXPORTS=(
    "CLASSIFY_MANIFEST=${CLASSIFY_MANIFEST}"
    "CLASSIFY_OUTDIR=${CLASSIFY_OUTDIR}"
    "EXTRACT_OUTDIR=${EXTRACT_OUTDIR}"
    "DB=${DB}"
    "THREADS=${CPUS}"
    "CONFIDENCE=${CONFIDENCE}"
    "MEMORY_MAPPING=${MEMORY_MAPPING}"
    "CONTAINER_SIF=${CONTAINER_SIF}"
    "CONTAINER_IMAGE=${CONTAINER_IMAGE}"
)
CLASSIFY_EXPORT_STR="ALL,$(IFS=,; echo "${CLASSIFY_EXPORTS[*]}")"

# ── Submit classify array job (if there are pending samples) ──────────────────
mkdir -p "${CLASSIFY_OUTDIR}"

if [[ ${#PENDING_INDICES[@]} -gt 0 ]]; then
    SBATCH_CLASSIFY=(
        sbatch
        --job-name=csc_classify
        --array="${PENDING_ARRAY_SPEC}"
        --cpus-per-task="${CPUS}"
        --mem="${MEM}"
        --time="${WALLTIME}"
        --partition="${PARTITION}"
        --output="${LOG_DIR}/classify_%A_%a.out"
        --error="${LOG_DIR}/classify_%A_%a.err"
        --export="${CLASSIFY_EXPORT_STR}"
        "${SCRIPT_DIR}/classify_array.sh"
    )

    echo ""
    echo "Classify sbatch command:"
    printf '  %s \\\n' "${SBATCH_CLASSIFY[@]}"
    echo ""

    if [[ "${DRY_RUN}" -eq 0 ]]; then
        CLASSIFY_JOB_OUTPUT="$("${SBATCH_CLASSIFY[@]}")"
        echo "${CLASSIFY_JOB_OUTPUT}"
        # Extract job ID from "Submitted batch job 12345"
        CLASSIFY_JOB_ID="$(echo "${CLASSIFY_JOB_OUTPUT}" | grep -oE '[0-9]+$')"
        echo "Classify job submitted: ${CLASSIFY_JOB_ID}"
    else
        echo "(Dry-run: classify job not submitted)"
        CLASSIFY_JOB_ID="<dry-run>"
    fi
fi

if [[ "${DRY_RUN}" -eq 0 ]] && [[ -n "${CLASSIFY_JOB_ID:-}" ]]; then
    echo ""
    echo "Classification job submitted: ${CLASSIFY_JOB_ID}"
    echo "  Classify (array):  job ${CLASSIFY_JOB_ID}  [${#PENDING_INDICES[@]} task(s)]"
    echo ""
    echo "To aggregate, detect, and generate a report, run:"
    echo "  ./submit_agg_report.sh \\"
    echo "      --outdir ${OUTDIR} \\"
    echo "      --extract-outdir ${EXTRACT_OUTDIR} \\"
    echo "      --db ${DB} \\"
    echo "      --dependency ${CLASSIFY_JOB_ID}"
    echo ""
    echo "Monitor classify with: squeue -j ${CLASSIFY_JOB_ID}"
fi
