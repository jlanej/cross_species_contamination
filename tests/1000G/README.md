# 1000 Genomes – Unmapped-Read Extraction

Scripts in this directory extract unmapped reads from 1000 Genomes
high-coverage CRAMs hosted on EBI endpoints.

Only the unmapped section of each remote CRAM is fetched (using the CRAI index
as a seek pointer), so the full ~30 GB file is never downloaded.

The CRAI index is downloaded to a local temporary file before calling
`samtools view -X` (Aspera first when available, `curl` fallback). This avoids
`htslib`/`samtools` failures seen when passing FTP URLs directly as the `-X`
index argument ("Exec format error"). The data CRAM stream is still read
directly from EBI.

All tools run inside an **Apptainer** (Singularity) container that is pulled
automatically on first use — no local software installation is required beyond
Apptainer itself.

---

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| **Apptainer ≥ 1.0** (or Singularity ≥ 3.x) | `module load apptainer` on most HPCs |
| **SLURM** | For job scheduling |
| Internet access to `ghcr.io` | Container is pulled on first run |
| Internet access to EBI CRAM/CRAI endpoints | `ftp.sra.ebi.ac.uk` and/or `ftp.1000genomes.ebi.ac.uk` |

The container (`ghcr.io/jlanej/cross_species_contamination:latest`) bundles
**samtools**, **Kraken2**, and all CSC Python tools.  It is built and published
automatically via GitHub Actions on every push to `main`.

---

## Files

| File | Purpose |
|------|---------|
| `manifest.tsv` | 3202-sample manifest with FTP URLs for every 1KG CRAM/CRAI |
| `extract_unmapped_array.sh` | SLURM array job – one task per sample |
| `submit_extract.sh` | Submission helper with subsetting options |

---

## Quick start

```bash
# Navigate to this directory on your HPC
cd tests/1000G

# Make scripts executable
chmod +x submit_extract.sh extract_unmapped_array.sh

# Load Apptainer (if not already in PATH)
module load apptainer   # command varies by HPC; skip if already available

# --- Option A: all 3202 samples (container is pulled automatically) ---
./submit_extract.sh

# --- Option B: first 50 samples (good for testing) ---
./submit_extract.sh --limit 50

# --- Option C: named subset ---
echo -e "NA12718\nNA12748\nNA18488" > my_subset.txt
./submit_extract.sh --samples my_subset.txt

# --- Option D: dry-run to preview the sbatch command ---
./submit_extract.sh --limit 10 --dry-run
```

The container SIF is saved as `csc.sif` next to `submit_extract.sh` and reused by all
subsequent jobs.  The pull happens **once**, in `submit_extract.sh` before the
SLURM array is submitted, so individual worker nodes never need write access to
pull the image themselves.  To use a pre-downloaded image:

```bash
./submit_extract.sh --container /shared/containers/csc.sif
```

---

## Output structure

```
output/
├── NA12718/
│   ├── NA12718_unmapped_R1.fastq.gz
│   ├── NA12718_unmapped_R2.fastq.gz
│   └── NA12718_unmapped_singleton.fastq.gz   # if any
├── NA12748/
│   └── ...
└── ...
```

Each sample directory contains gzip-compressed FASTQ files for the
unmapped reads.  Samples with zero unmapped reads will have an empty
directory.

By default, extraction also writes per-sample idxstats sidecars:

- `{sample}.idxstats.tsv`
- `{sample}.reads_summary.json`

These are used downstream to compute idxstats-based total-sequenced-read
metrics (absolute burden matrices and reporting).

---

## Options (`submit_extract.sh`)

| Option | Default | Description |
|--------|---------|-------------|
| `--manifest FILE` | `manifest.tsv` (same dir) | Path to the manifest |
| `--outdir DIR` | `./output` | Output base directory |
| `--samples FILE` | – | Text file with SAMPLE\_IDs (one per line) |
| `--limit N` | – | Process only the first N samples |
| `--range STR` | – | Raw SLURM `--array` spec, e.g. `1-50,60,100` |
| `--partition STR` | `normal` | SLURM partition name |
| `--cpus N` | `4` | CPUs per task |
| `--mem STR` | `8G` | Memory per task |
| `--time STR` | `02:00:00` | Wall-clock time limit |
| `--max-concurrent-jobs N` | `300` | Max concurrently running array tasks (`--array=<spec>%N`) |
| `--reference FILE` | – | GRCh38 FASTA for CRAM decoding (recommended) |
| `--container FILE` | `./csc.sif` | Path to pre-pulled Apptainer SIF image |
| `--image URI` | `ghcr.io/jlanej/cross_species_contamination:latest` | Docker URI for auto-pull |
| `--keep-cram` | off | Also save an intermediate unmapped CRAM |
| `--skip-idxstats` | off | Skip idxstats sidecars (default computes and requires them) |
| `--dry-run` | off | Print sbatch command without submitting |

---

## Container image

The default image is
`ghcr.io/jlanej/cross_species_contamination:latest`. You can pin a version:

```bash
# Pin to a specific release tag
./submit_extract.sh --image ghcr.io/jlanej/cross_species_contamination:v0.2.0
```

To build and use a local image:

```bash
docker build -t csc:local .
apptainer build csc_local.sif docker-daemon://csc:local
./submit_extract.sh --container csc_local.sif
```

---

## Resource guidance

| Samples | Recommended `--array` strategy |
|---------|-------------------------------|
| Test (≤10) | `--limit 10` |
| Small batch (≤100) | `--limit 100` |
| Full cohort | default (no flags) |

Each task typically completes in **10–30 minutes** for 30× coverage samples,
and uses **< 2 GB RAM**.  Most of the time is network I/O from the EBI FTP
server.  If your HPC has an egress firewall, ensure outbound FTP (port 21/20)
to `ftp.sra.ebi.ac.uk` and HTTPS (port 443) to `ghcr.io` are allowed.

---

## Restarting failed jobs

`submit_extract.sh` automatically excludes samples whose
`_unmapped_R1.fastq.gz` already exists and is non-empty, and exits cleanly if
everything selected is already complete. `extract_unmapped_array.sh` also
checks per-sample output and skips completed work at task runtime.

If you intentionally skipped idxstats sidecars during extraction, run
`submit_classify.sh --skip-idxstats-metrics` to bypass idxstats-based
absolute metrics in aggregation/reporting.

---

## Using a reference genome (recommended)

Providing a GRCh38 reference FASTA allows samtools to decode CRAM blocks
correctly and avoids warning messages:

```bash
./submit_extract.sh \
    --reference /path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

The EBI FTP also hosts the reference used for the 1KG CRAMs:
```
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

---

## Troubleshooting

### `Exec format error` when opening the CRAI

Some htslib builds (in particular older samtools packages distributed via
`apt`) cannot open FTP URLs as the index argument to `samtools view -X`.
The symptom is:

```
[E::hts_hopen] Failed to open file ftp://...NA12718.final.cram.crai
samtools view: failed to open "ftp://..." for reading: Exec format error
```

**This is fixed in `extract_unmapped_array.sh`**: the CRAI is downloaded
locally (Aspera first when available, `curl` fallback) before the samtools
pipeline runs. The CRAM data stream itself is still read directly from EBI
(only the unmapped section is transferred).

If you encounter this on a custom samtools invocation, apply the same pattern:

```bash
# 1. Download the small CRAI index locally
curl -fsSL --retry 3 -o /tmp/sample.crai "${CRAI_FTP_URL}"

# 2. Use the local CRAI; stream only the unmapped contig from the remote CRAM
samtools view -u -f 4 -X "${CRAM_FTP_URL}" /tmp/sample.crai '*' \
  | samtools collate -u -O - /tmp/collate_tmp \
  | samtools fastq -1 R1.fastq.gz -2 R2.fastq.gz -s singleton.fastq.gz -0 other.fastq.gz -
```

### CI / Automated testing of remote access

A dedicated GitHub Actions workflow
(`.github/workflows/remote-integration.yml`) exercises the full remote
pipeline against the EBI FTP server.  It runs automatically every Monday and
can also be triggered manually via the **Actions → Remote Integration** tab.
This means regressions in FTP connectivity or the samtools pipeline are
caught in CI rather than discovered during production runs.

---

## Downstream processing

The extracted FASTQ files can be fed directly into the full CSC pipeline:

```bash
# Build a samples CSV pointing at the extracted FASTQs
find output -name '*_R1.fastq.gz' | \
  awk -F'/' '{gsub(/_unmapped_R1.fastq.gz/,"",$NF); print $NF","$0}' \
  > samples_for_classify.csv

# Run the CSC Nextflow pipeline (classify → aggregate → detect)
nextflow run nextflow/main.nf \
    --input_csv  samples_for_classify.csv \
    --kraken2_db /path/to/prackendb \
    --outdir     results_1kg/ \
    -profile     slurm,apptainer
```

---

## SLURM array indexing

The manifest has a one-line header; sample data begins on line 2.
SLURM task IDs map to manifest lines as:

```
TASK_ID 1  → manifest line 2  (first sample)
TASK_ID 2  → manifest line 3  (second sample)
...
TASK_ID N  → manifest line N+1
```

So `--array=1-3202` covers all 3202 samples.

---

*AI assistance acknowledgment: scripts developed with AI assistance.*
