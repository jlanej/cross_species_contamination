# Report Module

The `csc-report` CLI generates a static, self-contained HTML summary of
non-human content in a cohort of human whole-genome sequencing (WGS)
samples directly from the outputs of `csc-aggregate` (and optionally
`csc-detect`).

> **AI Acknowledgment:** This module was developed with AI assistance.
> Best practices in bioinformatics should always take precedence over
> specific implementation details.

## Why a dedicated report module?

The three-tier normalisation emitted by `csc-aggregate`
(`taxa_matrix_raw.tsv`, `taxa_matrix_cpm.tsv`, `taxa_matrix_abs.tsv`)
answers three scientifically distinct questions and is easy to
mis-interpret if reported without labelled denominators.  This module
ensures that:

1. **CPM** (compositional, per million *classified* reads) and
   **absolute burden** (per million *total sequenced* reads) are shown
   side-by-side in every summary table, with denominators in column
   headers.
2. Every figure caption states its denominator explicitly.
3. A dedicated **Variant-Calling Impact** section uses the absolute
   denominator to flag samples whose non-human content exceeds a
   configurable threshold (default 0.1% of total sequencing).
4. The Methods section distinguishes the three read counts that the
   pipeline produces: *total sequenced*, *extracted for classification*,
   and *classified by Kraken2*.
5. A **Methods Transparency Checklist** records schema versions,
   parameters, and relative links to every input matrix so the report
   is fully reproducible and review-proof.

The renderer uses only the Python standard library (no Jinja2, pandas
or matplotlib dependencies) and emits one HTML file with inline CSS and
inline SVG figures.

## Usage

```bash
# Minimum: aggregate outputs only
csc-report aggregate_out/ -o report/contamination_report.html

# With detect outputs (adds §2.4 and flagged-sample counts)
csc-report aggregate_out/ -o report/contamination_report.html \
    --detect-dir detect_out/

# Tighter variant-calling impact threshold (100 ppm = 0.01% of sequencing)
csc-report aggregate_out/ -o report/contamination_report.html \
    --variant-impact-threshold-ppm 100

# Show more than the default top 10 taxa in the cohort-wide table
csc-report aggregate_out/ -o report/contamination_report.html --top-n 25
```

## Required inputs

`<aggregate_dir>/` must contain:

| File | Purpose |
|---|---|
| `taxa_matrix_raw.tsv` | Raw integer read counts |
| `taxa_matrix_cpm.tsv` | CPM (compositional) |
| `aggregation_metadata.json` | Schema version, parameters, provenance |

Optional but **strongly recommended** (enables §4 Variant-Calling
Impact and the absolute-burden figure):

| File | Purpose |
|---|---|
| `taxa_matrix_abs.tsv` | Absolute burden per million total sequenced reads |

Run `csc-aggregate --idxstats <path>/reads_summary.json ...` with the
per-sample sidecars emitted by `csc-extract` to produce it.

Optional detect-module inputs (from `csc-detect -o <detect_dir>`):

| File | Purpose |
|---|---|
| `flagged_samples.tsv` | Per (sample, taxon) outlier flags |
| `qc_summary.json` | Detect parameters, flagged samples |
| `quarantine_list.txt` | Plain-text quarantine list |

## Outputs

| File | Purpose |
|---|---|
| `<output>.html` | Self-contained HTML report |
| `<output-parent>/report_manifest.json` | Machine-readable manifest with the report schema version, generated timestamp, chosen threshold, and the list of samples flagged by §4 |

## Report sections

1. **Executive Summary** — cohort headline metrics (sample/taxon
   counts, median non-human composition and absolute burden, samples
   exceeding the variant-calling threshold).
2. **Methods** — extract / classify / aggregate / detect procedures and
   a clear statement distinguishing the three read counts.
3. **Results** — per-sample summary table, cohort-wide top-N
   contaminating taxa, and stacked-bar composition figures for both
   CPM and absolute burden.
4. **Variant-Calling Impact** — samples exceeding the configured
   absolute-burden threshold.
5. **Discussion & Caveats** — interpretation caveats, kitome /
   reagent-contamination warnings, spike-in / negative-control guidance.
6. **Methods Transparency Checklist** — versioned references and
   download links for every input matrix.

## Defensible denominator choices

| Matrix | Denominator | Report section where used |
|---|---|---|
| `taxa_matrix_raw.tsv` | — (integer counts) | Reproducibility table, transparency checklist |
| `taxa_matrix_cpm.tsv` | Classified direct reads | Composition stacked bar, % non-human of classified, diversity indices |
| `taxa_matrix_abs.tsv` | Total sequenced reads (idxstats) | Absolute-burden figure, §4 Variant-Calling Impact, manuscript QC tables |

See [docs/aggregate.md](aggregate.md) for the full denominator
cheat-sheet and [docs/extract.md](extract.md) for how the total
sequenced read count is captured at extract time.

## Python API

```python
from csc.report import load_inputs, generate_html_report

inputs = load_inputs("aggregate_out/", detect_dir="detect_out/")
generate_html_report(
    inputs,
    "report/contamination_report.html",
    top_n=25,
    threshold_ppm=1000.0,  # 0.1% of total sequencing
)
```
