#!/usr/bin/env bash
# =============================================================================
# submit_extract.sh
#
# Convenience wrapper to submit extract_unmapped_array.sh as a SLURM array
# job for any subset of the 1000 Genomes manifest.
#
# Usage:
#   # All samples
#   ./submit_extract.sh
#
#   # Only the first 100 samples
#   ./submit_extract.sh --limit 100
#
#   # A named subset (one SAMPLE_ID per line)
#   ./submit_extract.sh --samples my_subset.txt
#
#   # Explicit 1-based array range
#   ./submit_extract.sh --range 1-50
#
#   # Custom output directory and SLURM partition
#   ./submit_extract.sh --outdir /scratch/me/1kg_out --partition gpu
#
#   # Dry-run: print the sbatch command without submitting
#   ./submit_extract.sh --dry-run --limit 10
#
# Options:
#   --manifest FILE   Path to manifest.tsv  [default: same dir as this script]
#   --outdir   DIR    Output base directory [default: ./output]
#   --samples  FILE   Text file with SAMPLE_IDs to process (one per line)
#   --limit    N      Process only the first N samples
#   --range    STR    Explicit SLURM --array range string (e.g. "1-50,60,100")
#   --partition STR   SLURM partition       [default: normal]
#   --cpus     N      CPUs per task         [default: 4]
#   --mem      STR    Memory per task       [default: 8G]
#   --time     STR    Wall-clock time limit [default: 02:00:00]
#   --max-concurrent-jobs N  Max concurrently running array tasks [default: 300]
#   --reference FILE  Reference FASTA for CRAM decoding (optional)
#   --container FILE  Path to the Apptainer SIF image
#                     [default: <script_dir>/csc.sif; auto-pulled if absent]
#   --image     URI   Docker URI to pull the container from
#                     [default: ghcr.io/jlanej/cross_species_contamination:latest]
#   --keep-cram       Also save intermediate unmapped CRAM per sample
#   --dry-run         Print sbatch command; do not submit
#   -h, --help        Show this help message
#
# AI assistance acknowledgment: developed with AI assistance.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Defaults ────────────────────────────────────────────────────────────────
MANIFEST="${SCRIPT_DIR}/manifest.tsv"
OUTDIR="${SCRIPT_DIR}/output"
SAMPLES_FILE=""
LIMIT=""
RANGE=""
PARTITION="normal"
CPUS=4
MEM="8G"
WALLTIME="02:00:00"
REFERENCE=""
CONTAINER_SIF="${SCRIPT_DIR}/csc.sif"
CONTAINER_IMAGE="ghcr.io/jlanej/cross_species_contamination:latest"
KEEP_CRAM="0"
DRY_RUN=0
MAX_CONCURRENT_JOBS=300
MAX_CONCURRENT_JOBS_SET=0
# Keep usage output focused on the documented header/options block and long
# enough to include all currently documented options. If usage text grows,
# increase this so it still covers the full top comment usage/options block.
USAGE_LINES=80

# ── Argument parsing ─────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | sed 's/^# \{0,2\}//' | head -"${USAGE_LINES}"
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)  MANIFEST="$2";    shift 2 ;;
        --outdir)    OUTDIR="$2";      shift 2 ;;
        --samples)   SAMPLES_FILE="$2"; shift 2 ;;
        --limit)     LIMIT="$2";       shift 2 ;;
        --range)     RANGE="$2";       shift 2 ;;
        --partition) PARTITION="$2";   shift 2 ;;
        --cpus)      CPUS="$2";        shift 2 ;;
        --mem)       MEM="$2";         shift 2 ;;
        --time)      WALLTIME="$2";    shift 2 ;;
        --max-concurrent-jobs) MAX_CONCURRENT_JOBS="$2"; MAX_CONCURRENT_JOBS_SET=1; shift 2 ;;
        --reference) REFERENCE="$2";  shift 2 ;;
        --container) CONTAINER_SIF="$2"; shift 2 ;;
        --image)     CONTAINER_IMAGE="$2"; shift 2 ;;
        --keep-cram) KEEP_CRAM="1";   shift ;;
        --dry-run)   DRY_RUN=1;       shift ;;
        -h|--help)   usage ;;
        *) echo "ERROR: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ── Validate manifest ────────────────────────────────────────────────────────
if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}" >&2
    exit 1
fi

if ! [[ "${MAX_CONCURRENT_JOBS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --max-concurrent-jobs must be a positive integer." >&2
    exit 1
fi

# Total number of samples in the manifest (header = line 1, samples start at 2)
TOTAL_SAMPLES=$(( $(wc -l < "${MANIFEST}") - 1 ))
echo "Manifest: ${MANIFEST} (${TOTAL_SAMPLES} samples)"

# ── Resolve array range ──────────────────────────────────────────────────────
if [[ -n "${RANGE}" ]]; then
    ARRAY_SPEC="${RANGE}"

elif [[ -n "${SAMPLES_FILE}" ]]; then
    # Convert SAMPLE_IDs → 1-based line numbers in the manifest
    if [[ ! -f "${SAMPLES_FILE}" ]]; then
        echo "ERROR: Samples file not found: ${SAMPLES_FILE}" >&2
        exit 1
    fi

    # Build a lookup: SAMPLE_ID → manifest line number (1-based, skipping header)
    echo "Resolving sample IDs from ${SAMPLES_FILE}..."
    INDICES=()
    while IFS= read -r sid; do
        [[ -z "${sid}" ]] && continue
        # grep for the sample id in column 1 and get its line number relative to data
        IDX=$(awk -v s="${sid}" 'NR>1 && $1==s {print NR-1; exit}' FS='\t' "${MANIFEST}")
        if [[ -z "${IDX}" ]]; then
            echo "WARNING: Sample '${sid}' not found in manifest; skipping." >&2
        else
            INDICES+=("${IDX}")
        fi
    done < "${SAMPLES_FILE}"

    if [[ ${#INDICES[@]} -eq 0 ]]; then
        echo "ERROR: No matching samples found in manifest." >&2
        exit 1
    fi

    # Sort and de-duplicate indices, then format as SLURM array spec
    ARRAY_SPEC="$(printf '%s\n' "${INDICES[@]}" | sort -n | uniq | paste -sd',')"
    echo "Resolved ${#INDICES[@]} samples → array spec: ${ARRAY_SPEC}"

elif [[ -n "${LIMIT}" ]]; then
    ARRAY_SPEC="1-${LIMIT}"
else
    ARRAY_SPEC="1-${TOTAL_SAMPLES}"
fi

# Ensure SLURM array throttling is set as --array=<spec>%<max_concurrent_jobs>:
# preserve explicit '%' in --range, otherwise apply the configurable default.
if [[ "${ARRAY_SPEC}" == *%* ]]; then
    if [[ "${MAX_CONCURRENT_JOBS_SET}" -eq 1 ]]; then
        echo "ERROR: Cannot use --max-concurrent-jobs when --range already specifies '%' concurrency. Use only one concurrency method." >&2
        exit 1
    fi
else
    ARRAY_SPEC="${ARRAY_SPEC}%${MAX_CONCURRENT_JOBS}"
fi

# ── Create directories ───────────────────────────────────────────────────────
mkdir -p "${OUTDIR}" "${SCRIPT_DIR}/logs"

# ── Resolve container image – pull once before any array tasks start ─────────
# Ensures the SIF exists at an absolute, user-writable path so that SLURM
# worker nodes (which copy the script to /var/spool/slurmd/…) never need to
# pull or write the lock file themselves.
CONTAINER_SIF="$(realpath -m "${CONTAINER_SIF}")"

if [[ "${DRY_RUN}" -eq 0 ]] && [[ ! -f "${CONTAINER_SIF}" ]]; then
    # Detect apptainer / singularity
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

# ── Build sbatch command ─────────────────────────────────────────────────────
# CONTAINER_SIF is always passed explicitly so array tasks never fall back to
# ${SCRIPT_DIR}/csc.sif (which resolves to SLURM's /var/spool/slurmd/… temp dir).
ARRAY_JOB_EXPORTS="CONTAINER_SIF=${CONTAINER_SIF},CONTAINER_IMAGE=${CONTAINER_IMAGE}"
[[ -n "${REFERENCE}" ]] && ARRAY_JOB_EXPORTS+=",REFERENCE=${REFERENCE}"

SBATCH_CMD=(
    sbatch
    --job-name=1kg_extract
    --array="${ARRAY_SPEC}"
    --cpus-per-task="${CPUS}"
    --mem="${MEM}"
    --time="${WALLTIME}"
    --partition="${PARTITION}"
    --output="${SCRIPT_DIR}/logs/extract_%A_%a.out"
    --error="${SCRIPT_DIR}/logs/extract_%A_%a.err"
    --export="ALL,MANIFEST=${MANIFEST},OUTDIR=${OUTDIR},THREADS=${CPUS},KEEP_CRAM=${KEEP_CRAM},${ARRAY_JOB_EXPORTS}"
    "${SCRIPT_DIR}/extract_unmapped_array.sh"
)

# ── Submit or print ──────────────────────────────────────────────────────────
echo ""
echo "sbatch command:"
printf '  %s \\\n' "${SBATCH_CMD[@]}"
echo ""

if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "(Dry-run: not submitting)"
else
    "${SBATCH_CMD[@]}"
fi
