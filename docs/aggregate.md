# Aggregate Module

> **Status:** Implemented.

## Overview

The **aggregate** module (`csc.aggregate`) collects per-sample Kraken2
classification reports and produces a sample-by-taxon matrix.  Values can
be raw direct-read counts or counts-per-million (CPM) normalised.

The implementation is designed to handle very large cohorts (100 K+ samples)
efficiently by processing reports in configurable chunks.

## Python API

```python
from csc.aggregate import parse_kraken2_report, aggregate_reports

# Parse a single report
records = parse_kraken2_report("sample.kraken2.report.txt")

# Build a matrix from multiple reports
result = aggregate_reports(
    ["sampleA.kraken2.report.txt", "sampleB.kraken2.report.txt"],
    output_dir="results/",
    min_reads=10,       # exclude taxa with fewer than 10 reads per sample
    normalize=True,     # CPM normalisation (default)
    chunk_size=500,     # reports per processing chunk
    rank_filter=("S", "G", "F"),  # per-rank matrices (default)
)
print(result["matrix_path"])    # results/taxa_matrix.tsv
print(result["metadata_path"])  # results/aggregation_metadata.json
print(result["rank_matrices"])  # {'S': Path('results/taxa_matrix_S.tsv'), ...}
```

### Key Functions

| Function | Description |
|---|---|
| `parse_kraken2_report(path)` | Parse a Kraken2 report into a list of `TaxonRecord` dicts |
| `aggregate_reports(paths, output_dir, ...)` | Build the sample×taxon matrix; returns `AggregationResult` |
| `sample_id_from_report(path)` | Derive sample ID by stripping `.kraken2.report.txt` suffix |

### Return Types

**`AggregationResult`** (TypedDict):
- `matrix_path` – Path to the unfiltered output TSV matrix
- `metadata_path` – Path to the JSON metadata file
- `sample_count` – Number of samples in the matrix
- `taxon_count` – Number of taxa (rows) in the unfiltered matrix
- `rank_matrices` – Dict mapping rank code to per-rank filtered matrix path
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
# Basic usage – CPM-normalised matrix
csc-aggregate reports/*.kraken2.report.txt -o results/

# Filter low-count taxa
csc-aggregate reports/*.kraken2.report.txt -o results/ --min-reads 50

# Raw counts (no normalisation)
csc-aggregate reports/*.kraken2.report.txt -o results/ --no-normalize

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
| `--no-normalize` | Output raw counts instead of CPM | `False` |
| `--chunk-size` | Reports per processing chunk | `500` |
| `--rank-filter` | Taxonomy rank codes for per-rank matrices | `S G F` |
| `--json-log` | Structured JSON logging | `False` |
| `-v, --verbose` | DEBUG-level logging | `False` |

## Output Files

### `taxa_matrix.tsv`

Tab-separated matrix.  Rows are taxa, columns are samples.

```
tax_id  name                   sampleA   sampleB   sampleC
562     Escherichia coli       500.0000  12500.00  0.0000
1280    Staphylococcus aureus  2500.00   0.0000    3000.00
```

First two columns are `tax_id` and `name`.  Remaining columns are
sample IDs derived from report file names.  Missing taxa are filled
with `0`.

When `normalize=True` (default), values are CPM (counts per million):
each sample column sums to 1,000,000.

### `taxa_matrix_S.tsv`, `taxa_matrix_G.tsv`, `taxa_matrix_F.tsv`

Per-rank filtered matrices containing only taxa of the specified
rank.  Produced for each rank in `--rank-filter` (default: S, G, F).

### `rank_filter_metadata.json`

```json
{
  "rank_filter": ["S", "G", "F"],
  "ranks": {
    "S": {
      "matrix_path": "results/taxa_matrix_S.tsv",
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
  "normalized": true,
  "normalization_method": "CPM",
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
