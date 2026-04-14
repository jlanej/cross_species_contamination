# 1000 Genomes – Unmapped-Read Extraction

Scripts in this directory extract unmapped reads from the 1000 Genomes Project
high-coverage CRAMs hosted on the EBI FTP server.

Only the unmapped section of each remote CRAM is fetched (using the CRAI index
as a seek pointer), so the full ~30 GB file is never downloaded.

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
# Clone / navigate to this directory on your HPC
cd tests/1000G

# Make scripts executable
chmod +x submit_extract.sh extract_unmapped_array.sh

# --- Option A: all 3 202 samples ---
./submit_extract.sh

# --- Option B: first 50 samples (good for testing) ---
./submit_extract.sh --limit 50

# --- Option C: named subset ---
echo -e "NA12718\nNA12748\nNA18488" > my_subset.txt
./submit_extract.sh --samples my_subset.txt

# --- Option D: dry-run to preview the sbatch command ---
./submit_extract.sh --limit 10 --dry-run
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
| `--reference FILE` | – | GRCh38 FASTA for CRAM decoding (recommended) |
| `--keep-cram` | off | Also save an intermediate unmapped CRAM |
| `--dry-run` | off | Print sbatch command without submitting |

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
to `ftp.sra.ebi.ac.uk` is allowed.

---

## Restarting failed jobs

The job script skips samples whose `_unmapped_R1.fastq.gz` already exists and
is non-empty.  Simply resubmit the same `--array` spec and only failed samples
will be reprocessed.

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
    -profile     slurm,singularity
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
