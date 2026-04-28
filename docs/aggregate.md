# Aggregate Module

> **Status:** Implemented.

## Overview

The **aggregate** module (`csc.aggregate`) collects per-sample Kraken2
classification reports and produces sample-by-taxon matrices in **both**
raw direct-read counts and counts-per-million (CPM).

The implementation is designed to handle very large cohorts (100 K+ samples)
efficiently by processing reports in configurable chunks.

## High-Confidence Tier (Dual-Tier Reporting)

By default Kraken2 is run with `--confidence 0.0` (maximum sensitivity).
This is great for catching trace contamination but can produce false
positives caused by low-complexity regions or sparse k-mer matches.
The aggregate module supports an optional **high-confidence tier**
that recomputes per-read confidence from the existing per-read
output files (`*.kraken2.output.txt`) and demotes weak classifications
to "unclassified" before building the matrices.  This avoids re-running
Kraken2 entirely.

The confidence definition matches Kraken2's own algorithm
(`classify.cc`):

```
confidence(taxon) = (k-mers whose taxon ∈ clade rooted at taxon)
                    -------------------------------------------------
                    (total k-mers excluding ambiguous A:N runs)
```

Multiple thresholds may be supplied; each non-zero threshold produces
a parallel matrix set with suffix `_conf{T}` (e.g.
`taxa_matrix_raw_conf0p50.tsv`, `taxa_matrix_cpm_S_conf0p50.tsv`).
The canonical sensitive tier (no suffix) is always written.

### CLI

```bash
csc-aggregate \
    sampleA.kraken2.report.txt sampleB.kraken2.report.txt \
    -o results/ \
    --db-path /path/to/kraken2_db \
    --kraken2-output sampleA.kraken2.output.txt sampleB.kraken2.output.txt \
    --confidence-threshold 0.1 0.5
```

Files supplied to `--kraken2-output` are matched to samples by
filename (the `.kraken2.output.txt` suffix is stripped).  The flag
requires `--db-path` (we need `taxonomy/nodes.dmp` to walk the
clade).  `--confidence-threshold 0.0` is treated as a no-op.

### Python API

```python
result = aggregate_reports(
    [...],
    output_dir="results/",
    db_path="/path/to/kraken2_db",
    kraken2_output_paths=[...],
    confidence_thresholds=[0.1, 0.5],
)
for tier_suffix, tier in result.get("confidence_tiers", {}).items():
    print(tier_suffix, tier["threshold"], tier["matrix_raw_path"])
```

`aggregation_metadata.json` records every tier under
`confidence_tiers.<suffix>` with the threshold, output filenames,
and per-sample `reads_demoted_to_unclassified` counters for
auditability.

### Downstream propagation

`csc-detect` automatically discovers sibling confidence-tier matrices
next to the input matrix and runs detection on each, writing results
into `<output>/<tier_suffix>/` (e.g. `detect/conf0p50/`).  Disable with
`--no-confidence-tiers`.  Per-rank tier matrices are picked up
automatically as well (`<output>/<tier_suffix>/<rank>/`).

The Nextflow workflow exposes `--confidence_thresholds '0.1 0.5'`
and the 1000G driver exposes `CONFIDENCE_THRESHOLDS=0.1:0.5`
(colon-separated to avoid spaces in `sbatch --export`).



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
| `--idxstats` | Per-sample `reads_summary.json` paths (from `csc-extract`); enables `taxa_matrix_abs.tsv` | _none_ |
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

### `taxa_matrix_abs.tsv` (absolute contaminant burden)

Produced **only** when `--idxstats PATH [PATH ...]` supplies the
per-sample `reads_summary.json` sidecars emitted by `csc-extract`
(see [docs/extract.md](extract.md)).

Each cell is the contaminant read count expressed as
`raw_count / total_reads * 1_000_000` – i.e. **per-million total
sequenced reads**, not per-million classified reads.  This is the
absolute contamination burden metric recommended by
[Natarajan et al., Nat Biotechnol 2023](https://www.nature.com/articles/s41587-023-01732-6)
and is the correct denominator when judging downstream impact of
contamination on variant calling, assembly, or expression quantitation.

Samples whose `sample_id` does not match any supplied idxstats sidecar
are written as `NA` in every row of the absolute matrix, and listed
under `samples_without_idxstats` in `aggregation_metadata.json` so
downstream reports can flag them clearly.

> **Denominator cheat-sheet**
>
> | Matrix | Cell value | Denominator | Interpretation |
> |---|---|---|---|
> | `taxa_matrix_raw.tsv` | `direct_reads` (or `clade_reads` for G/F) | – | Integer counts |
> | `taxa_matrix_cpm.tsv` | `raw / classified_total * 1e6` | Sum of classified direct reads (pre-filter) | Composition **among contaminants** |
> | `taxa_matrix_abs.tsv` | `raw / total_reads * 1e6` | `total_mapped + total_unmapped` at extract | **Absolute burden** per million sequenced |

Per-rank absolute matrices (`taxa_matrix_abs_S.tsv`,
`taxa_matrix_abs_G.tsv`, …) are written alongside the raw/CPM rank
matrices for the same rank codes requested via `--rank-filter`.

### `rank_filter_metadata.json`

```json
{
  "rank_filter": ["S", "G", "F"],
  "ranks": {
    "S": {
      "matrix_raw_path": "results/taxa_matrix_raw_S.tsv",
      "matrix_cpm_path": "results/taxa_matrix_cpm_S.tsv",
      "matrix_abs_path": "results/taxa_matrix_abs_S.tsv",
      "taxon_count": 15,
      "taxa": [
        {"tax_id": 562, "name": "Escherichia coli"},
        {"tax_id": 1280, "name": "Staphylococcus aureus"}
      ]
    }
  }
}
```

The `matrix_abs_path` entry is present only when idxstats sidecars were
supplied.

### `aggregation_metadata.json`

```json
{
  "schema_version": "1.1",
  "sample_count": 3,
  "taxon_count": 25,
  "min_reads": 10,
  "matrix_paths": {
    "raw": "taxa_matrix_raw.tsv",
    "cpm": "taxa_matrix_cpm.tsv",
    "abs": "taxa_matrix_abs.tsv"
  },
  "absolute_burden_enabled": true,
  "samples": ["sampleA", "sampleB", "sampleC"],
  "samples_without_idxstats": [],
  "sample_provenance": {
    "sampleA": {
      "total_reads": 912469134,
      "total_mapped": 900123456,
      "total_unmapped": 12345678,
      "source_input": "/abs/path/sampleA.bam",
      "extraction_time": "2024-01-01T12:00:00+00:00"
    }
  },
  "errors": []
}
```

`sample_provenance` carries the exact denominator used for each
sample's absolute-burden column, plus the source BAM path and
extraction timestamp.  This is the canonical join key for reporting
modules (e.g. the upcoming **PR #55 non-human contamination summary
report**) – it lets figure captions cite "denominator = `total_reads`
from extract, captured at `extraction_time`, for input `source_input`".

## Downstream integration — breadcrumbs for PR #55

The PR #55 summary-report agent should:

1. Load `taxa_matrix_abs.tsv` alongside `taxa_matrix_cpm.tsv` so every
   heatmap / scatterplot carries both a **relative** (composition among
   contaminants) and an **absolute** (per-million sequenced) panel.
2. Join sample columns against
   `aggregation_metadata.json → sample_provenance[sample_id]` for
   figure-caption provenance (input BAM, `total_reads`,
   `extraction_time`).
3. Honour `samples_without_idxstats` by rendering such samples greyed
   out / "NA" in absolute-burden plots rather than silently omitting
   them.
4. Check `schema_version` in the metadata and refuse to render if the
   schema is newer than what the report knows about.

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
