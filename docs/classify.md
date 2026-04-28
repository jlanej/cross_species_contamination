# Classify Module

> **Status:** ✅ Implemented

## Overview

The **classify** module (`csc.classify`) provides taxonomic classification of
FASTQ reads extracted by the `csc.extract` module using
[Kraken2](https://ccb.jhu.edu/software/kraken2/).

It wraps the Kraken2 binary, validates databases, and produces standardized
report and per-read classification output files suitable for downstream
aggregation by `csc.aggregate`.

## Python API

```python
from csc.classify import classify_reads, validate_database

# Validate a Kraken2 database
db_path = validate_database("/data/kraken2/PlusPF")

# Classify single-end reads
result = classify_reads(
    ["sample.unmapped.fastq.gz"],
    output_dir="results/",
    db="/data/kraken2/PlusPF",
    confidence=0.2,
    threads=4,
)

# Classify paired-end reads
result = classify_reads(
    ["sample.R1.fastq.gz", "sample.R2.fastq.gz"],
    output_dir="results/",
    db="/data/kraken2/PlusPF",
    paired=True,
)

# Result dict
print(result["report"])       # Path to Kraken2 report file
print(result["output"])       # Path to per-read classification file
print(result["sample_id"])    # Sample identifier
print(result["input_files"])  # List of input file paths
```

### `classify_reads()` Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_files` | `list[Path]` | required | FASTQ files to classify |
| `output_dir` | `Path` | required | Output directory |
| `db` | `Path` | required | Kraken2 database directory |
| `sample_id` | `str` | auto | Override output file basename |
| `confidence` | `float` | `0.0` | Minimum confidence threshold (0.0–1.0) |
| `threads` | `int` | `1` | Number of Kraken2 threads |
| `memory_mapping` | `bool` | `False` | Use memory mapping (less RAM, slower) |
| `paired` | `bool` | `False` | Paired-end mode (requires 2 input files) |

### `validate_database()` Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `db_path` | `Path` | Path to Kraken2 database directory |

Raises `FileNotFoundError` if the directory does not exist, or `ValueError`
if required database files (`hash.k2d`, `opts.k2d`, `taxo.k2d`) are missing.

## CLI Usage

```bash
# Basic single-end classification
csc-classify reads.fastq.gz --db /data/kraken2/PlusPF -o results/

# Paired-end classification
csc-classify R1.fastq.gz R2.fastq.gz --db /data/kraken2/PlusPF -o results/ --paired

# With confidence threshold and multiple threads
csc-classify reads.fastq.gz --db /data/kraken2/PlusPF -o results/ \
    --confidence 0.2 --threads 8

# Memory-mapped mode (lower RAM footprint)
csc-classify reads.fastq.gz --db /data/kraken2/PlusPF -o results/ --memory-mapping

# Custom sample ID
csc-classify reads.fastq.gz --db /data/kraken2/PlusPF -o results/ --sample-id SAMPLE_001

# Structured JSON logging
csc-classify reads.fastq.gz --db /data/kraken2/PlusPF -o results/ --json-log

# Verbose debug logging
csc-classify reads.fastq.gz --db /data/kraken2/PlusPF -o results/ -v
```

### CLI Options

| Flag | Description |
|------|-------------|
| `input` | Input FASTQ file(s) (positional) |
| `-o, --output-dir` | Output directory (required) |
| `--db` | Kraken2 database path (required) |
| `--confidence` | Confidence threshold (default: 0.0) |
| `--threads` | Thread count (default: 1) |
| `--memory-mapping` | Use memory-mapped DB |
| `--paired` | Paired-end mode |
| `--sample-id` | Override sample ID |
| `--json-log` | JSON structured logging |
| `-v, --verbose` | Debug logging |
| `--version` | Show version |

> **Note – sensitive vs. high-confidence reporting.**  We recommend
> running `csc-classify` with the default `--confidence 0.0` (maximum
> sensitivity) so the per-read kraken2 output captures every k-mer
> match, then producing a high-confidence tier downstream in
> `csc-aggregate` via `--confidence-threshold`.  This avoids re-running
> Kraken2: the per-read confidence is recomputed from the existing
> `*.kraken2.output.txt` file using the same definition Kraken2 itself
> uses (in-clade k-mers ÷ non-ambiguous k-mers).  See
> [Aggregate – High-Confidence Tier](aggregate.md#high-confidence-tier-dual-tier-reporting).

## Chaining with Extract

The classify module is designed to accept output directly from `csc-extract`:

```bash
# Step 1: Extract unmapped reads
csc-extract sample.bam -o extracted/

# Step 2: Classify the extracted reads
csc-classify extracted/sample.unmapped.R1.fastq.gz \
             extracted/sample.unmapped.R2.fastq.gz \
    --db /data/kraken2/PlusPF \
    -o classified/ \
    --paired
```

## Output Files

For a sample named `SAMPLE_001`, the module produces:

| File | Description |
|------|-------------|
| `SAMPLE_001.kraken2.report.txt` | Kraken2 summary report (tab-separated) |
| `SAMPLE_001.kraken2.output.txt` | Per-read classification assignments |

### Report Format

The Kraken2 report is a tab-separated file with columns:

1. Percentage of reads rooted at this taxon
2. Number of reads rooted at this taxon
3. Number of reads directly assigned to this taxon
4. Rank code (U, R, D, P, C, O, F, G, S)
5. NCBI taxonomy ID
6. Scientific name

### Output Format

The per-read output file has columns:

1. C (classified) or U (unclassified)
2. Read name
3. Taxonomy ID (0 if unclassified)
4. Read length
5. LCA k-mer mapping

## Memory Management

Kraken2 databases are large (PrackenDB is ~70 GB on disk).  By default,
Kraken2 loads the entire ``hash.k2d`` file into RAM before classifying reads.
On systems where available RAM is limited, this can cause out-of-memory (OOM)
errors or excessive swapping.

### Estimating RAM Before Classification

Use `csc-db estimate-memory` to check how much RAM will be required to load a
database and whether memory-mapping is recommended on the current host:

```bash
# Quick check — text output
csc-db estimate-memory /data/kraken2/PlusPF

# Machine-readable JSON output
csc-db estimate-memory /data/kraken2/PlusPF --json
```

Example output:

```
Database: /data/kraken2/PlusPF
  hash.k2d size (primary RAM consumer): 55.3 GiB
  Total database size:                  56.1 GiB
  Estimated RAM without --memory-mapping: 55.3 GiB
  Available system RAM:                  62.0 GiB
  Recommendation: --memory-mapping not required  (database fits in available RAM)
```

The Python API is also available:

```python
from csc.classify import estimate_db_memory

est = estimate_db_memory("/data/kraken2/PlusPF")
print(est["human"]["estimated_ram"])   # e.g. "55.3 GiB"
print(est["recommend_memory_mapping"]) # True / False
```

### Memory Mapping (`--memory-mapping`)

Memory mapping instructs Kraken2 to access the database via OS-level page
mapping instead of loading the full hash table into RAM upfront.  Individual
database pages are loaded on demand as reads are classified, which means:

- **Lower peak RAM usage:** Only the pages actually accessed are loaded.
- **Slower classification:** On-demand page faults add latency — expect
  roughly 2–5× slower throughput compared to a fully in-RAM database.
- **Shared page cache:** Multiple concurrent Kraken2 processes pointing at the
  same database file share OS page-cache pages, making memory-mapping
  especially efficient on multi-sample HPC nodes.

Enable memory mapping via the CLI or Python API:

```bash
# CLI
csc-classify reads.fastq.gz --db /data/kraken2/PlusPF -o results/ --memory-mapping
```

```python
result = classify_reads(
    ["reads.fastq.gz"],
    output_dir="results/",
    db="/data/kraken2/PlusPF",
    memory_mapping=True,
)
```

**Guideline:** use `csc-db estimate-memory` to check the RAM requirement; if
the recommended action is `--memory-mapping`, pass that flag to `csc-classify`
(or set `memory_mapping: true` in your config file).

### Troubleshooting: Limited or Variable Resource Environments

| Symptom | Likely cause | Solution |
|---------|-------------|----------|
| `Killed` or OOM during classification | Database does not fit in RAM | Use `--memory-mapping` |
| Very slow classification | Memory mapping with slow spinning disk | Pre-copy DB to local SSD; or increase RAM |
| `kraken2 not found on PATH` | Kraken2 not installed | Install Kraken2 (`conda install -c bioconda kraken2`) |
| Repeated DB downloads | `CSC_DB_CACHE` not set persistently | `export CSC_DB_CACHE=/shared/kraken2_dbs` in your shell profile |
| OOM on HPC with many parallel jobs | Each job loads the DB independently | Use `--memory-mapping`; jobs share OS page cache for the same file |
| Docker: DB not found | DB not mounted | `docker run -v /host/db:/db ... --db /db` |

#### Nextflow / HPC: memory-safe classification

Set memory-mapping in your Nextflow run to avoid OOM on nodes with limited RAM:

```bash
nextflow run nextflow/main.nf \
    --kraken2_db /data/kraken2/PlusPF \
    --memory_mapping \
    --classify_memory "16 GB"
```

Alternatively, set the default in `nextflow/nextflow.config`:

```groovy
params {
    memory_mapping = true
}
```

#### CI / Cloud: caching the database

```bash
# GitHub Actions: cache DB between workflow runs
- uses: actions/cache@v4
  with:
    path: ~/.csc/db
    key: kraken2-db-${{ hashFiles('db-version.txt') }}
- run: csc-db fetch prackendb
- run: csc-db estimate-memory ~/.csc/db/prackendb
```



### Recommended: PrackenDB

**PrackenDB** is the recommended Kraken2 database for the CSC pipeline.  It is
a curated, pre-built database published by the Kraken2 project containing one
NCBI reference genome per species — bacteria, archaea, fungi, protists,
viruses, human, and UniVec Core.

**Why PrackenDB?**

- **Unambiguous species-level detection:** Because each species contributes
  exactly one reference genome, k-mer counts are unambiguous and LCA
  assignments are not inflated by redundant genomes.
- **Robust contamination detection:** Avoids false positives from shared k-mers
  between redundant reference genomes of the same species.
- **Taxonomy files included:** PrackenDB includes `taxonomy/nodes.dmp` and
  `taxonomy/names.dmp`, enabling lineage-aware classification (e.g. tracing a
  species taxid to its domain) and human-readable taxon names.
- **Comparable with KDF:** Results are directly comparable with the
  [kmer_denovo_filter](https://github.com/jlanej/kmer_denovo_filter) pipeline.

See the [KDF PrackenDB documentation](https://github.com/jlanej/kmer_denovo_filter/blob/main/docs/kraken2_bacterial_detection.md#the-prackendb-reference-database)
for full details on why single-genome-per-species databases matter.

```bash
# Fetch PrackenDB (recommended)
csc-db fetch prackendb

# Verify taxonomy files
csc-db verify /path/to/prackendb
```

```python
from csc.classify import fetch_prackendb, is_prackendb, validate_taxonomy

# Fetch PrackenDB
db_path = fetch_prackendb()

# Check PrackenDB compatibility
assert is_prackendb(db_path)

# Verify taxonomy files
tax = validate_taxonomy(db_path)
assert tax["taxonomy/nodes.dmp"]
assert tax["taxonomy/names.dmp"]
```

> **Warning:** If you use a non-PrackenDB database, the pipeline will emit a
> warning about possible species-level ambiguity.  This does not prevent
> execution, but results may be less reliable for species-level detection.

### Required Database Files

A valid Kraken2 database directory must contain:

- `hash.k2d` — hash table
- `opts.k2d` — options/parameters
- `taxo.k2d` — taxonomy data

### Public Databases

Pre-built databases from
[Ben Langmead's index zone](https://benlangmead.github.io/aws-indexes/k2):

| Database | Size | Description |
|----------|------|-------------|
| **PrackenDB** | ~70 GB | **Recommended.** One genome per species (bacteria, archaea, fungi, protists, viruses, human, UniVec Core) |
| **PlusPF** | ~70 GB | Standard + protozoa + fungi (multiple genomes per species) |
| **Standard** | ~50 GB | RefSeq archaea, bacteria, viral, plasmid, human, UniVec_Core |
| **PlusPFP** | ~150 GB | PlusPF + plant |
| **Standard-8** | ~8 GB | Standard capped at 8 GB |

```bash
# Download PrackenDB (recommended)
csc-db fetch prackendb
```

### Custom Databases

Build a custom database using `kraken2-build`:

```bash
kraken2-build --download-taxonomy --db /data/kraken2/custom_db
kraken2-build --download-library bacteria --db /data/kraken2/custom_db
kraken2-build --build --db /data/kraken2/custom_db --threads 16
```

See the [Kraken2 manual](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)
for full details on custom database construction.

## Configuration

Default settings in `csc/default_config.yaml`:

```yaml
classify:
  tool: "kraken2"
  db: null                   # path to classifier database
  confidence: 0.0            # minimum confidence threshold
  threads: 1
  memory_mapping: false      # use memory mapping instead of loading DB into RAM
  recommended_db: "prackendb"  # recommended: PrackenDB (one genome per species)
```

Override via a user config file:

```yaml
# my_config.yaml
classify:
  db: "/data/kraken2/PlusPF"
  confidence: 0.2
  threads: 8
  memory_mapping: true
```

```bash
export CSC_CONFIG=my_config.yaml
```

## Nextflow Integration

```bash
nextflow run nextflow/classify.nf \
    --input_csv samples.csv \
    --kraken2_db /data/kraken2/PlusPF \
    --outdir results/ \
    --confidence 0.2 \
    --threads 8
```

Input CSV format:

```csv
sample_id,fastq
SAMPLE_001,/data/extracted/SAMPLE_001.unmapped.fastq.gz
SAMPLE_002,/data/extracted/SAMPLE_002.unmapped.fastq.gz
```

For paired-end data, add a `fastq2` column:

```csv
sample_id,fastq,fastq2
SAMPLE_001,/data/extracted/SAMPLE_001.R1.fastq.gz,/data/extracted/SAMPLE_001.R2.fastq.gz
```

## Docker

```bash
# Build image (includes Kraken2)
docker build -t csc-pipeline .

# Run classification with mounted database
docker run --rm \
    -v /data/kraken2/PlusPF:/data/kraken2_db:ro \
    -v /data/fastqs:/data/input:ro \
    -v /data/output:/data/output \
    --entrypoint csc-classify \
    csc-pipeline \
    /data/input/sample.fastq.gz \
    --db /data/kraken2_db \
    -o /data/output/
```
