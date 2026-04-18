# Aggregate Module

> **Status:** Implemented.

## Overview

The **aggregate** module (`csc.aggregate`) collects per-sample Kraken2
classification reports and produces sample-by-taxon matrices in **both**
raw direct-read counts and counts-per-million (CPM).

The implementation is designed to handle very large cohorts (100 K+ samples)
efficiently by processing reports in configurable chunks.

## Python API

```python
from csc.aggregate import parse_kraken2_report, aggregate_reports

# Parse a single report
records = parse_kraken2_report("sample.kraken2.report.txt")

# Build both raw and CPM matrices from multiple reports
result = aggregate_reports(
    ["sampleA.kraken2.report.txt", "sampleB.kraken2.report.txt"],
    output_dir="results/",
    min_reads=10,       # exclude taxa with fewer than 10 reads per sample
    chunk_size=500,     # reports per processing chunk
    rank_filter=("S", "G", "F"),  # per-rank matrices (default)
)
print(result["matrix_raw_path"])  # results/taxa_matrix_raw.tsv
print(result["matrix_cpm_path"])  # results/taxa_matrix_cpm.tsv
print(result["metadata_path"])    # results/aggregation_metadata.json
print(result["rank_matrices_raw"])  # {'S': Path('results/taxa_matrix_raw_S.tsv'), ...}
print(result["rank_matrices_cpm"])  # {'S': Path('results/taxa_matrix_cpm_S.tsv'), ...}
```

### Key Functions

| Function | Description |
|---|---|
| `parse_kraken2_report(path)` | Parse a Kraken2 report into a list of `TaxonRecord` dicts |
| `aggregate_reports(paths, output_dir, ...)` | Build the sample×taxon matrix; returns `AggregationResult` |
| `sample_id_from_report(path)` | Derive sample ID by stripping `.kraken2.report.txt` suffix |

### Return Types

**`AggregationResult`** (TypedDict):
- `matrix_raw_path` – Path to the unfiltered raw-count TSV matrix
- `matrix_cpm_path` – Path to the unfiltered CPM-normalised TSV matrix
- `metadata_path` – Path to the JSON metadata file
- `sample_count` – Number of samples in the matrix
- `taxon_count` – Number of taxa (rows) in the unfiltered matrix
- `rank_matrices_raw` – Dict mapping rank code to per-rank raw matrix path
- `rank_matrices_cpm` – Dict mapping rank code to per-rank CPM matrix path
- `rank_metadata_path` – Path to the rank-filter sidecar JSON

**`TaxonRecord`** (TypedDict):
- `tax_id` – NCBI taxonomy ID
- `name` – Scientific name
- `rank` – Rank code (U, R, D, P, C, O, F, G, S, …)
- `clade_reads` – Reads rooted at this taxon
- `direct_reads` – Reads directly assigned to this taxon
- `percentage` – Percentage of total reads

## CLI Usage

```bash
# Basic usage – writes both CPM and raw matrices
csc-aggregate reports/*.kraken2.report.txt -o results/

# Filter low-count taxa
csc-aggregate reports/*.kraken2.report.txt -o results/ --min-reads 50

# Custom rank filter (species only)
csc-aggregate reports/*.kraken2.report.txt -o results/ --rank-filter S

# Multiple rank filters
csc-aggregate reports/*.kraken2.report.txt -o results/ --rank-filter S G

# Verbose + JSON logging
csc-aggregate reports/*.kraken2.report.txt -o results/ -v --json-log
```

### CLI Options

| Flag | Description | Default |
|---|---|---|
| `input` | Kraken2 report files (positional, one or more) | required |
| `-o, --output-dir` | Output directory | required |
| `--min-reads` | Minimum direct reads to include a taxon | `0` |
| `--chunk-size` | Reports per processing chunk | `500` |
| `--rank-filter` | Taxonomy rank codes for per-rank matrices | `S G F` |
| `--json-log` | Structured JSON logging | `False` |
| `-v, --verbose` | DEBUG-level logging | `False` |

## Output Files

### `taxa_matrix_raw.tsv` and `taxa_matrix_cpm.tsv`

Tab-separated matrix.  Rows are taxa, columns are samples.

```
tax_id  name                   sampleA   sampleB   sampleC
562     Escherichia coli       500       12500     0
1280    Staphylococcus aureus  2500      0         3000
```

First two columns are `tax_id` and `name`.  Remaining columns are
sample IDs derived from report file names.  Missing taxa are filled
with `0`.

`taxa_matrix_raw.tsv` stores integer direct-read counts.  
`taxa_matrix_cpm.tsv` stores CPM values where each sample column sums to
approximately 1,000,000.

> **CPM denominator**: The CPM denominator uses the sum of all classified
> direct reads *before* any `min_reads` filtering.  This ensures CPM
> values are comparable across samples even when different `min_reads`
> thresholds are used, avoiding compositional bias.

> **CPM limitations**: CPM normalisation accounts for sequencing depth
> but **not** genome size or gene length.  Species with larger genomes
> will have higher CPM at equal true abundance.  For applications where
> genome size matters, consider additional RPKM or TPM normalisation
> outside this pipeline.  CPM values should be interpreted as relative
> read-count proportions, not true relative abundances.

> **Confidence intervals**: The current output does not include
> confidence intervals or error bars for count/CPM values.  For
> low-read-count taxa, statistical uncertainty may be substantial.
> Users requiring uncertainty quantification should consider Poisson
> confidence intervals (e.g. ±1.96 × √N) as a post-processing step.

### `taxa_matrix_raw_S.tsv`, `taxa_matrix_cpm_S.tsv`, etc.

Per-rank filtered matrices containing only taxa of the specified rank.
Both a raw-count version and a CPM-normalised version are written for
each rank code requested via `--rank-filter` (e.g. `S`, `G`, `F`).

For species-rank (`S`) matrices, values are `direct_reads` (reads
assigned directly to that species).  For higher-rank matrices (`G`, `F`,
etc.), values are `clade_reads` (reads rooted at that taxon, including
all descendant species).  This ensures genus/family matrices capture
total abundance rather than being misleadingly sparse.

### `rank_filter_metadata.json`

```json
{
  "rank_filter": ["S", "G", "F"],
  "ranks": {
    "S": {
      "matrix_raw_path": "results/taxa_matrix_raw_S.tsv",
      "matrix_cpm_path": "results/taxa_matrix_cpm_S.tsv",
      "taxon_count": 15,
      "taxa": [
        {"tax_id": 562, "name": "Escherichia coli"},
        {"tax_id": 1280, "name": "Staphylococcus aureus"}
      ]
    }
  }
}
```

### `aggregation_metadata.json`

```json
{
  "sample_count": 3,
  "taxon_count": 25,
  "min_reads": 10,
  "matrix_paths": {
    "raw": "taxa_matrix_raw.tsv",
    "cpm": "taxa_matrix_cpm.tsv"
  },
  "samples": ["sampleA", "sampleB", "sampleC"],
  "errors": []
}
```

## Configuration

Default thresholds are set in `csc/default_config.yaml`:

```yaml
aggregate:
  min_reads: 10
  rank_filter:
    - "S"
    - "G"
    - "F"
```

These can be overridden via a user config file or the CLI `--min-reads`
flag.

### Configuration vs CLI Precedence

When using the Python CLI (`csc-aggregate`), the `--min-reads` default
is **0** (include all taxa).  The `default_config.yaml` sets
`min_reads: 10`, which is used by the Nextflow pipeline and when
loading config programmatically via `load_config()`.

**Precedence order** (highest to lowest):
1. Explicit CLI argument (e.g. `--min-reads 50`)
2. CLI default value (`0`)

To use the config file value of `10`, pass it explicitly:
`--min-reads 10`.

### Duplicate Sample IDs

If two report files produce the same sample ID (e.g. same filename in
different directories), a warning is logged and the later report
overwrites the earlier one.  Rename report files to ensure unique
sample IDs.

### Taxonomy Consistency

If a `tax_id` is seen with different names or rank codes across input
reports (e.g. due to different Kraken2 database versions), a warning
is logged.  The last-seen name/rank is used.  Users should ensure all
input reports were generated with the same Kraken2 database.
