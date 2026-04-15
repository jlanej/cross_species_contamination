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
#   ASPERA_SSH_KEY – optional explicit path to the Aspera SSH key file
#   ASPERA_BANDWIDTH – Aspera transfer cap (default: 300m)
#   ASPERA_PORT   – Aspera transfer port (default: 33001)
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
# Default to OUTDIR (a user-writable path) so that if this script is submitted
# directly (without submit_extract.sh) it still uses a writable location.
# When submitted via submit_extract.sh, CONTAINER_SIF is always passed
# explicitly as an absolute path via --export.
CONTAINER_SIF="${CONTAINER_SIF:-${OUTDIR}/csc.sif}"

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

find_aspera_key() {
    local candidates=(
        "${ASPERA_SSH_KEY:-}"
        "/etc/asperaweb_id_dsa.openssh"
        "/usr/etc/asperaweb_id_dsa.openssh"
        "${HOME}/.aspera/connect/etc/asperaweb_id_dsa.openssh"
        "/opt/aspera/etc/asperaweb_id_dsa.openssh"
    )
    local key
    for key in "${candidates[@]}"; do
        [[ -n "${key}" && -f "${key}" ]] && { echo "${key}"; return 0; }
    done
    return 1
}

download_crai_with_aspera() {
    local ftp_url="$1"
    local dest="$2"
    local aspera_path=""
    local aspera_user=""
    local aspera_host=""
    local aspera_key=""
    local aspera_port="${ASPERA_PORT:-33001}"
    local aspera_bandwidth="${ASPERA_BANDWIDTH:-300m}"
    local dest_dir remote_basename downloaded

    command -v ascp &>/dev/null || return 1
    aspera_key="$(find_aspera_key)" || return 1

    if [[ "${ftp_url}" == *"ftp.sra.ebi.ac.uk"* ]]; then
        aspera_path="${ftp_url#ftp://ftp.sra.ebi.ac.uk/}"
        aspera_user="era-fasp"
        aspera_host="fasp.sra.ebi.ac.uk"
    elif [[ "${ftp_url}" == *"ftp.1000genomes.ebi.ac.uk"* ]]; then
        aspera_path="${ftp_url#ftp://ftp.1000genomes.ebi.ac.uk/}"
        aspera_user="fasp-g1k"
        aspera_host="fasp.1000genomes.ebi.ac.uk"
    else
        return 1
    fi

    dest_dir="$(dirname "${dest}")"
    remote_basename="$(basename "${aspera_path}")"
    downloaded="${dest_dir}/${remote_basename}"

    # -T: disable encryption for speed on public data; -r: resume interrupted
    # transfer support; -Q: adaptive flow control; -L-: no local log dir.
    if ! ascp -i "${aspera_key}" \
        -Tr -Q -l "${aspera_bandwidth}" -P"${aspera_port}" -L- \
        "${aspera_user}@${aspera_host}:${aspera_path}" \
        "${dest_dir}/"; then
        rm -f "${downloaded}" "${dest}"
        return 1
    fi

    if [[ "${downloaded}" != "${dest}" && -f "${downloaded}" && ! -f "${dest}" ]]; then
        mv -f "${downloaded}" "${dest}"
    elif [[ "${downloaded}" != "${dest}" && -f "${downloaded}" && -f "${dest}" ]]; then
        echo "WARNING: Aspera produced both '${downloaded}' and '${dest}'; keeping '${dest}'." >&2
        rm -f "${downloaded}"
    fi

    [[ -s "${dest}" ]] || { rm -f "${dest}"; return 1; }
    return 0
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
# Download CRAI index to a local file                                           #
#                                                                               #
# samtools view -X does not reliably support FTP URLs for the index file       #
# (some htslib builds report "Exec format error" when the URL scheme is ftp:// #
# even though the same CRAM data stream works fine).  Downloading the small    #
# index file locally first avoids this limitation while preserving the         #
# efficient seek-to-unmapped behaviour.                                         #
# --------------------------------------------------------------------------- #
CRAI_LOCAL="${SAMPLE_DIR}/.${SAMPLE_ID}.crai.tmp"
echo "Downloading CRAI index for ${SAMPLE_ID}..."
if download_crai_with_aspera "${CRAI_URL}" "${CRAI_LOCAL}"; then
    echo "CRAI downloaded via Aspera ($(wc -c < "${CRAI_LOCAL}") bytes)"
elif curl -fsSL --retry 3 --retry-delay 5 -o "${CRAI_LOCAL}" "${CRAI_URL}"; then
    echo "CRAI downloaded via curl ($(wc -c < "${CRAI_LOCAL}") bytes)"
else
    echo "ERROR: Failed to download CRAI index: ${CRAI_URL}" >&2
    rm -f "${CRAI_LOCAL}"
    exit 1
fi

# --------------------------------------------------------------------------- #
# Build samtools view arguments                                                 #
# We fetch ONLY the '*' (unmapped) virtual contig from the remote CRAM,        #
# using the local CRAI so samtools can seek directly to that section without   #
# reading the entire file.                                                      #
# --------------------------------------------------------------------------- #
VIEW_ARGS=(
    samtools view
    --threads "${THREADS}"
    -u          # uncompressed BAM on stdout (piped to fastq)
    -f 4        # FLAG: read unmapped
    -X "${CRAM_URL}" "${CRAI_LOCAL}"
    '*'
)

if [[ -n "${REFERENCE:-}" ]]; then
    VIEW_ARGS=(
        samtools view
        -T "${REFERENCE}"
        --threads "${THREADS}"
        -u -f 4
        -X "${CRAM_URL}" "${CRAI_LOCAL}"
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
    container_run "${VIEW_ARGS[@]}" -C -o "${UNMAPPED_CRAM}" || {
        echo "ERROR: Failed to save intermediate CRAM for ${SAMPLE_ID}" >&2
        exit 1
    }
fi

# --------------------------------------------------------------------------- #
# Stream unmapped reads → collate by name → convert to FASTQ                  #
# Run the entire samtools pipeline in a single container invocation so that    #
# the inter-process pipes stay within one container context.                   #
# --------------------------------------------------------------------------- #
echo "Extracting unmapped reads and converting to FASTQ..."

R2="${SAMPLE_DIR}/${SAMPLE_ID}_unmapped_R2.fastq.gz"
SINGLETON="${SAMPLE_DIR}/${SAMPLE_ID}_unmapped_singleton.fastq.gz"
OTHER="${SAMPLE_DIR}/${SAMPLE_ID}_unmapped_other.fastq.gz"

# Write the pipeline as a temporary shell script so variables with special
# characters (paths, URLs) are handled safely without shell re-expansion.
# All user-controlled values are escaped with printf '%q' before being
# written into the script to prevent any shell injection.
PIPELINE_SCRIPT="${SAMPLE_DIR}/.pipeline_${SLURM_ARRAY_TASK_ID}.sh"

# Shell-escape every value that comes from external input
q_threads=$(printf '%q' "${THREADS}")
q_crai_local=$(printf '%q' "${CRAI_LOCAL}")
q_cram_url=$(printf '%q' "${CRAM_URL}")
q_collate_tmp=$(printf '%q' "${COLLATE_TMP}/tmp")
q_r1=$(printf '%q' "${R1}")
q_r2=$(printf '%q' "${R2}")
q_singleton=$(printf '%q' "${SINGLETON}")
q_other=$(printf '%q' "${OTHER}")

{
    printf '#!/bin/bash\nset -euo pipefail\n'

    if [[ "${KEEP_CRAM}" == "1" ]]; then
        q_unmapped_cram=$(printf '%q' "${SAMPLE_DIR}/${SAMPLE_ID}_unmapped.cram")
        printf 'samtools collate --threads %s -u -O %s %s \\\n' \
            "${q_threads}" "${q_unmapped_cram}" "${q_collate_tmp}"
    else
        # samtools view: fetch only the unmapped virtual contig ('*') via local CRAI
        # The CRAI has been downloaded locally to avoid FTP URL limitations in
        # some htslib builds (see download step above).
        if [[ -n "${REFERENCE:-}" ]]; then
            q_ref=$(printf '%q' "${REFERENCE}")
            printf 'samtools view -T %s --threads %s -u -f 4 -X %s %s %s \\\n' \
                "${q_ref}" "${q_threads}" "${q_cram_url}" "${q_crai_local}" "'*'"
        else
            printf 'samtools view --threads %s -u -f 4 -X %s %s %s \\\n' \
                "${q_threads}" "${q_cram_url}" "${q_crai_local}" "'*'"
        fi
        printf '| samtools collate --threads %s -u -O - %s \\\n' \
            "${q_threads}" "${q_collate_tmp}"
    fi

    printf '| samtools fastq --threads %s \\\n    -1 %s \\\n    -2 %s \\\n    -s %s \\\n    -0 %s \\\n    -\n' \
        "${q_threads}" "${q_r1}" "${q_r2}" "${q_singleton}" "${q_other}"
} > "${PIPELINE_SCRIPT}"

chmod +x "${PIPELINE_SCRIPT}"
EXIT_CODE=0
container_run bash "${PIPELINE_SCRIPT}" || EXIT_CODE=$?

# Clean up temporary files
rm -rf "${COLLATE_TMP}"
rm -f "${PIPELINE_SCRIPT}"
rm -f "${CRAI_LOCAL}"

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
