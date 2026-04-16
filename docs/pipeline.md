# End-to-End Pipeline Usage

The CSC pipeline runs four stages in sequence:

```
BAM/CRAM ──► extract ──► classify ──► aggregate ──► detect ──► summary report
```

## Quick Start

```bash
nextflow run nextflow/main.nf \
    --input_csv  samples.csv \
    --kraken2_db /data/kraken2/PlusPF \
    --outdir     results/
```

## Input CSV Format

The pipeline reads a CSV file with required columns `sample_id` and `file`.
An optional `reference` column can supply per-sample reference FASTA paths
(needed for CRAM inputs).

```csv
sample_id,file,reference
SAMPLE_A,/data/bams/sample_a.bam,
SAMPLE_B,/data/crams/sample_b.cram,/data/refs/GRCh38.fa
```

## Parameters

### General

| Parameter     | Default    | Description                              |
|---------------|------------|------------------------------------------|
| `--input_csv` | *required* | Sample sheet CSV (columns: sample_id, file) |
| `--kraken2_db`| *required* | Path to Kraken2 database directory       |
| `--outdir`    | `results`  | Output directory                         |

### Extract Stage

| Parameter           | Default | Description                              |
|---------------------|---------|------------------------------------------|
| `--mapq`            | `null`  | MAPQ threshold (null = unmapped only)    |
| `--extract_cpus`    | `4`     | CPUs per extraction task                 |
| `--extract_memory`  | `4 GB`  | Memory per extraction task               |
| `--extract_time`    | `2h`    | Wall-time limit per extraction task      |

### Classify Stage

| Parameter             | Default | Description                            |
|-----------------------|---------|----------------------------------------|
| `--confidence`        | `0.0`   | Kraken2 confidence threshold           |
| `--memory_mapping`    | `false` | Use memory mapping for Kraken2 DB      |
| `--classify_cpus`     | `4`     | CPUs per classification task           |
| `--classify_memory`   | `8 GB`  | Memory per classification task         |
| `--classify_time`     | `4h`    | Wall-time limit per classification task|

### Aggregate Stage

| Parameter              | Default | Description                           |
|------------------------|---------|---------------------------------------|
| `--min_reads`          | `0`     | Min reads per taxon for inclusion     |
| `--no_normalize`       | `false` | Use raw counts for legacy primary matrix names |
| `--detect_matrix`      | `cpm`   | Matrix type used by detect stage (`cpm` or `raw`) |
| `--aggregate_cpus`     | `2`     | CPUs for aggregation                  |
| `--aggregate_memory`   | `4 GB`  | Memory for aggregation               |
| `--aggregate_time`     | `1h`    | Wall-time limit for aggregation      |

### Detect Stage

| Parameter                   | Default | Description                     |
|-----------------------------|---------|---------------------------------|
| `--detect_method`           | `mad`   | Outlier method (`mad` or `iqr`) |
| `--mad_threshold`           | `3.5`   | MAD threshold                   |
| `--iqr_multiplier`          | `1.5`   | IQR fence multiplier            |
| `--kitome_taxa`             | `null`  | Space-separated NCBI tax IDs    |
| `--no_subtract_background`  | `false` | Skip background subtraction     |
| `--detect_cpus`             | `2`     | CPUs for detection              |
| `--detect_memory`           | `4 GB`  | Memory for detection            |
| `--detect_time`             | `1h`    | Wall-time limit for detection   |

## Execution Profiles

Run with `-profile <name>` to select an execution environment:

```bash
# Local execution (default)
nextflow run nextflow/main.nf -profile standard ...

# Docker
nextflow run nextflow/main.nf -profile docker ...

# Singularity
nextflow run nextflow/main.nf -profile singularity ...

# SLURM cluster
nextflow run nextflow/main.nf -profile slurm ...

# Minimal resources for testing
nextflow run nextflow/main.nf -profile test,docker ...
```

### Profile Details

| Profile       | Executor  | Container Engine | Notes                    |
|---------------|-----------|------------------|--------------------------|
| `standard`    | local     | none             | Default profile          |
| `docker`      | local     | Docker           | Requires Docker          |
| `singularity` | local     | Singularity      | Auto-mounts enabled      |
| `slurm`       | SLURM     | none             | Set `--queue` as needed  |
| `local`       | local     | none             | Alias for standard       |
| `test`        | local     | none             | Reduced resource limits  |

Combine profiles with a comma: `-profile test,docker`.

## Output Structure

```
results/
├── extract/
│   └── SAMPLE_A/
│       ├── SAMPLE_A.unmapped.R1.fastq.gz
│       └── SAMPLE_A.unmapped.R2.fastq.gz
├── classify/
│   └── SAMPLE_A/
│       ├── SAMPLE_A.kraken2.report.txt
│       └── SAMPLE_A.kraken2.output.txt
├── aggregate/
│   ├── taxa_matrix.tsv
│   ├── taxa_matrix_cpm.tsv
│   ├── taxa_matrix_raw.tsv
│   └── aggregation_metadata.json
├── detect/
│   ├── flagged_samples.tsv
│   ├── qc_summary.json
│   └── quarantine_list.txt
├── pipeline_report.html
├── csc_pipeline_mqc.yaml
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    ├── trace.txt
    └── dag.html
```

## MultiQC Integration

The pipeline produces a `csc_pipeline_mqc.yaml` file that is compatible
with [MultiQC](https://multiqc.info/) custom content. Run MultiQC on the
output directory to include CSC results in a combined QC report:

```bash
multiqc results/
```

## Resume Capability

Nextflow tracks task state in its `work/` directory. If a run fails
partway, simply add `-resume` to restart from the last successful step:

```bash
nextflow run nextflow/main.nf -resume \
    --input_csv samples.csv \
    --kraken2_db /data/kraken2/PlusPF \
    --outdir results/
```

## Error Handling

- Each process uses `errorStrategy 'finish'` — other independent samples
  will continue even if one sample fails.
- Processes can be retried once by default (`maxRetries 1`).
- On pipeline failure the completion handler prints the error, duration,
  and work directory for diagnosis.
- All logs are preserved in the Nextflow work directory for post-mortem
  inspection.

## Docker

Build the container locally:

```bash
docker build -t csc-pipeline .
```

Then run with the Docker profile:

```bash
nextflow run nextflow/main.nf -profile docker \
    --input_csv samples.csv \
    --kraken2_db /data/kraken2/PlusPF \
    --outdir results/
```
