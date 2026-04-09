# Cross-Species Contamination Detection Pipeline

A multi-module Python pipeline for detecting cross-species contamination in
whole-genome sequencing BAM/CRAM files.

> **AI Acknowledgment:** This project was developed with AI assistance.
> Best practices in bioinformatics should always take precedence over specific
> implementation details.

## Overview

The pipeline is organised into four modules that run sequentially:

```
BAM/CRAM ──► extract ──► classify ──► aggregate ──► detect
                │            │            │            │
            FASTQ files   taxa labels   summary     report
```

| Module | Package | Status |
|--------|---------|--------|
| **Extract** | `csc.extract` | ✅ Implemented |
| **Classify** | `csc.classify` | ✅ Implemented |
| **Aggregate** | `csc.aggregate` | ✅ Implemented |
| **Detect** | `csc.detect` | ✅ Implemented |

Shared helpers live in `csc.utils` and pipeline-wide settings are managed by
`csc.config` (see [docs/configuration.md](docs/configuration.md)).

## Quick Start

### Installation

```bash
# From source
pip install .

# Or with test dependencies
pip install ".[test]"
```

Requires **samtools ≥ 1.12** on `PATH`.  
For classification, requires **Kraken2** on `PATH` and a Kraken2 database.  
**Recommended:** [PrackenDB](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads) —
a curated database with one genome per species for robust species-level detection
(see [Kraken2 Database Management](#kraken2-database-management)).

### Basic Usage

```bash
# Extract unmapped reads from a BAM file
csc-extract sample.bam -o output_dir/

# Extract unmapped reads from a CRAM file (requires reference)
csc-extract sample.cram -o output_dir/ --reference GRCh38.fa

# Also include poorly mapped reads (MAPQ < 10)
csc-extract sample.bam -o output_dir/ --mapq 10

# Use 4 decompression threads
csc-extract sample.bam -o output_dir/ --threads 4

# Interleaved output (single FASTQ file)
csc-extract sample.bam -o output_dir/ --interleaved
```

### Classify Reads

```bash
# Classify extracted reads with Kraken2
csc-classify output_dir/sample.unmapped.R1.fastq.gz \
             output_dir/sample.unmapped.R2.fastq.gz \
    --db /data/kraken2/PlusPF -o classified/ --paired

# Single-end classification with confidence threshold
csc-classify reads.fastq.gz --db /data/kraken2/PlusPF -o classified/ --confidence 0.2

# Memory-mapped mode (lower RAM usage)
csc-classify reads.fastq.gz --db /data/kraken2/PlusPF -o classified/ --memory-mapping
```

### Kraken2 Database Management

The `csc-db` CLI manages Kraken2 database downloads, caching, and verification.
Databases are cached in `~/.csc/db` by default (override with `--cache-dir` or
the `CSC_DB_CACHE` environment variable).

> **Recommended:** Use **PrackenDB** — a curated Kraken2 database containing one
> NCBI reference genome per species (bacteria, archaea, fungi, protists, viruses,
> human, and UniVec Core).  Because each species contributes exactly one reference,
> k-mer counts are unambiguous and LCA assignments are not inflated by redundant
> genomes.  This is essential for robust species-level contamination detection.
>
> See the [KDF PrackenDB documentation](https://github.com/jlanej/kmer_denovo_filter/blob/main/docs/kraken2_bacterial_detection.md#the-prackendb-reference-database)
> for details on why single-genome-per-species matters.

```bash
# Fetch the recommended PrackenDB database (one genome per species)
csc-db fetch prackendb

# Fetch a database from a public URL
csc-db fetch https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240605.tar.gz

# Fetch with SHA-256 verification
csc-db fetch https://example.com/db.tar.gz --sha256 abc123...

# Fetch from S3 (requires AWS CLI on PATH)
csc-db fetch s3://my-bucket/kraken2/PlusPF.tar.gz --name PlusPF

# Use a local directory directly (no download)
csc-db fetch /data/kraken2/PlusPF

# List cached databases
csc-db list

# Show database metadata (file sizes, SHA-256 hashes)
csc-db info /data/kraken2/PlusPF
csc-db info /data/kraken2/PlusPF --json   # JSON output

# Verify a database directory
csc-db verify /data/kraken2/PlusPF

# Clean the cache (all or specific)
csc-db clean
csc-db clean --name old_db
```

#### Docker / CI / Cloud Workflows

```bash
# Docker: download a small test DB at build time
docker build --build-arg DB_URL=https://example.com/k2_mini.tar.gz -t csc .

# CI (GitHub Actions): cache the DB between runs
- uses: actions/cache@v4
  with:
    path: ~/.csc/db
    key: kraken2-db-${{ hashFiles('db-version.txt') }}
- run: csc-db fetch https://example.com/k2_standard_08gb.tar.gz

# HPC shared filesystem: point to a pre-built DB
export CSC_DB_CACHE=/shared/kraken2_dbs
csc-db fetch /shared/kraken2/PlusPF --name PlusPF
```

#### Python API

```python
from csc.classify import (
    fetch_database, fetch_prackendb, database_info,
    list_databases, is_prackendb, validate_taxonomy,
)

# Download and cache PrackenDB (recommended)
db_path = fetch_prackendb()

# Check taxonomy files
tax = validate_taxonomy(db_path)
print(tax)  # {'taxonomy/nodes.dmp': True, 'taxonomy/names.dmp': True}

# Check if a database is PrackenDB-compatible
print(is_prackendb(db_path))  # True

# Download and cache a database from URL
db_path = fetch_database(
    "https://example.com/k2_standard.tar.gz",
    name="standard",
    expected_hash="abc123...",
    hash_algorithm="sha256",
)

# Inspect a database
info = database_info(db_path)
print(info["valid"], info["total_size_bytes"], info["sha256"])

# List cached databases
for db in list_databases():
    print(db["name"], db["size_bytes"], db["valid"])
```

### Configuration

All modules read a shared YAML config.  Override defaults by setting the
`CSC_CONFIG` environment variable or passing `--config`:

```bash
export CSC_CONFIG=my_config.yaml
```

See [docs/configuration.md](docs/configuration.md) for the full reference.

## End-to-End Nextflow Pipeline

For biobank-scale cohorts, use the included Nextflow pipeline that runs all
four stages in sequence:

```bash
nextflow run nextflow/main.nf \
    --input_csv  samples.csv \
    --kraken2_db /data/kraken2/PlusPF \
    --outdir     results/
```

The CSV file must have columns `sample_id` and `file`, with an optional
`reference` column:

```csv
sample_id,file,reference
SAMPLE_001,/data/SAMPLE_001.bam,
SAMPLE_002,/data/SAMPLE_002.cram,/ref/GRCh38.fa
```

The pipeline produces a MultiQC-compatible summary (`csc_pipeline_mqc.yaml`)
and an HTML report.  Run `multiqc results/` to combine with other QC tools.

### Individual Stage Pipelines

The individual extraction and classification pipelines are also available:

```bash
# Extraction only
nextflow run nextflow/extract_unmapped.nf --input_csv samples.csv --outdir results/

# Classification only
nextflow run nextflow/classify.nf --input_csv samples.csv --kraken2_db /data/kraken2/PlusPF
```

### Execution Profiles

```bash
# Local Docker execution
nextflow run nextflow/main.nf -profile docker --input_csv samples.csv --kraken2_db /data/k2

# Singularity on HPC
nextflow run nextflow/main.nf -profile singularity --input_csv samples.csv --kraken2_db /data/k2

# SLURM cluster
nextflow run nextflow/main.nf -profile slurm --input_csv samples.csv --kraken2_db /data/k2

# Resume after failure
nextflow run nextflow/main.nf -resume --input_csv samples.csv --kraken2_db /data/k2
```

See [docs/pipeline.md](docs/pipeline.md) for the full parameter reference,
output structure, and resource configuration.

## Docker

```bash
# Build (includes samtools and Kraken2)
docker build -t csc-pipeline .

# Run extraction
docker run --rm -v /data:/data csc-pipeline \
    /data/sample.bam -o /data/output/

# Run classification (mount Kraken2 database)
docker run --rm \
    -v /data/kraken2/PlusPF:/data/kraken2_db:ro \
    -v /data:/data \
    --entrypoint csc-classify \
    csc-pipeline \
    /data/output/sample.unmapped.fastq.gz \
    --db /data/kraken2_db -o /data/classified/
```

## Testing

Tests use synthetic BAM data generated by `tests/generate_test_data.py`—no
external datasets are required.

```bash
# Install test dependencies
pip install pysam ".[test]"

# Run all tests
pytest -v tests/

# Generate test data manually (optional)
python tests/generate_test_data.py /tmp/test_data
```

### What the Tests Cover

| Test | Description |
|------|-------------|
| `test_extract.py::TestValidation` | Input validation (missing files, bad extensions) |
| `test_extract.py::TestBuildCommand` | Command construction for samtools |
| `test_extract.py::TestExtractUnmappedOnly` | Unmapped-read extraction count verification |
| `test_extract.py::TestExtractWithMAPQ` | MAPQ-filtered extraction correctness |
| `test_extract.py::TestCLI` | CLI entry-point and error handling |
| `test_classify.py::TestValidateDatabase` | Kraken2 database validation |
| `test_classify.py::TestBuildClassifyCommand` | Kraken2 command construction |
| `test_classify.py::TestClassifyReads` | Classification with mocked Kraken2 |
| `test_classify.py::TestCLI` | Classify CLI entry-point and error handling |
| `test_db.py::TestComputeHash` | Hash computation (MD5, SHA-256) |
| `test_db.py::TestVerifyHash` | Hash verification and mismatch detection |
| `test_db.py::TestFetchDatabase` | DB download from local/HTTP/S3 (mocked) |
| `test_db.py::TestListDatabases` | Cache listing and validation status |
| `test_db.py::TestCleanCache` | Cache cleanup (specific and all) |
| `test_db.py::TestDatabaseInfo` | DB metadata and SHA-256 reporting |
| `test_db.py::TestDBCLI` | csc-db CLI entry-point and subcommands |
| `test_db.py::TestValidateTaxonomy` | Taxonomy file validation (nodes.dmp, names.dmp) |
| `test_db.py::TestIsPrackenDB` | PrackenDB detection heuristic |
| `test_db.py::TestPrackenDBConstants` | PrackenDB constants and URL verification |
| `test_db.py::TestFetchPrackenDB` | PrackenDB fetch with taxonomy validation |
| `test_db.py::TestConfigRecommendedDB` | Default config recommends PrackenDB |
| `test_config.py` | Config loading, merging, env-var override, error handling |

## Project Structure

```
├── csc/                        # Main package
│   ├── __init__.py             # Version & top-level docstring
│   ├── config.py               # Central YAML config loader
│   ├── default_config.yaml     # Built-in default settings
│   ├── extract/                # Extraction module
│   │   ├── __init__.py
│   │   ├── extract.py          # Core extraction logic
│   │   └── cli.py              # CLI entry point
│   ├── classify/               # Classification module
│   │   ├── __init__.py
│   │   ├── classify.py         # Core Kraken2 classification logic
│   │   ├── cli.py              # Classify CLI entry point
│   │   ├── db.py               # Database download, cache & hash verification
│   │   └── db_cli.py           # csc-db CLI entry point
│   ├── aggregate/              # Aggregation module
│   │   ├── __init__.py
│   │   ├── aggregate.py
│   │   └── cli.py
│   ├── detect/                 # Detection module
│   │   ├── __init__.py
│   │   ├── detect.py
│   │   ├── report.py
│   │   └── cli.py
│   └── utils/                  # Shared utility helpers
│       └── __init__.py
├── docs/                       # Module & pipeline documentation
│   ├── README.md
│   ├── configuration.md
│   ├── extract.md
│   ├── classify.md
│   ├── aggregate.md
│   └── detect.md
├── nextflow/
│   ├── main.nf                 # End-to-end pipeline (extract→classify→aggregate→detect)
│   ├── modules/
│   │   ├── extract.nf          # Extract process module
│   │   ├── classify.nf         # Classify process module
│   │   ├── aggregate.nf        # Aggregate process module
│   │   ├── detect.nf           # Detect process module
│   │   └── summary.nf          # Pipeline summary report module
│   ├── extract_unmapped.nf     # Standalone extraction pipeline
│   ├── classify.nf             # Standalone classification pipeline
│   └── nextflow.config         # Nextflow profiles & resource settings
├── tests/
│   ├── conftest.py             # Pytest fixtures
│   ├── generate_test_data.py   # Synthetic BAM generation
│   ├── test_extract.py         # Extraction test suite
│   ├── test_classify.py        # Classification test suite
│   ├── test_db.py              # DB management test suite
│   ├── test_config.py          # Config loader tests
│   ├── test_aggregate.py       # Aggregation test suite
│   ├── test_detect.py          # Detection test suite
│   └── test_nextflow_pipeline.py  # Pipeline structure tests
├── Dockerfile
├── .github/workflows/ci.yml   # GitHub Actions CI
├── pyproject.toml
└── README.md
```

## Documentation

See the [docs/](docs/) directory for per-module documentation:

- [Pipeline usage](docs/pipeline.md) — end-to-end Nextflow workflow
- [Extract](docs/extract.md) — BAM/CRAM read extraction
- [Classify](docs/classify.md) — Kraken2 taxonomic classification
- [Aggregate](docs/aggregate.md) — sample-by-taxon matrix building
- [Detect](docs/detect.md) — statistical outlier detection
- [Configuration](docs/configuration.md) — YAML config reference

## License

MIT – see [LICENSE](LICENSE).
