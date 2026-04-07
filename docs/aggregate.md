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
)
print(result["matrix_path"])    # results/taxa_matrix.tsv
print(result["metadata_path"])  # results/aggregation_metadata.json
```

### Key Functions

| Function | Description |
|---|---|
| `parse_kraken2_report(path)` | Parse a Kraken2 report into a list of `TaxonRecord` dicts |
| `aggregate_reports(paths, output_dir, ...)` | Build the sample×taxon matrix; returns `AggregationResult` |
| `sample_id_from_report(path)` | Derive sample ID by stripping `.kraken2.report.txt` suffix |

### Return Types

**`AggregationResult`** (TypedDict):
- `matrix_path` – Path to the output TSV matrix
- `metadata_path` – Path to the JSON metadata file
- `sample_count` – Number of samples in the matrix
- `taxon_count` – Number of taxa (rows) in the matrix

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
```

These can be overridden via a user config file or the CLI `--min-reads`
flag.
