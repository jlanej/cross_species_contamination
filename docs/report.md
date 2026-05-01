# Report Module

The `csc-report` CLI generates a static, self-contained HTML summary of
non-human content in a cohort of human whole-genome sequencing (WGS)
samples directly from the outputs of `csc-aggregate` (and optionally
`csc-detect`).

> **AI Acknowledgment:** This module was developed with AI assistance.
> Best practices in bioinformatics should always take precedence over
> specific implementation details.

## What changed in schema 2.1 (dual-tier confidence integration)

`REPORT_SCHEMA_VERSION` is now **`2.1`** (additive only).  When
`csc-aggregate` emits sibling high-confidence matrices
(`taxa_matrix_*_conf*.tsv`, the new default — see
[aggregate.md](aggregate.md#high-confidence-tier-dual-tier-reporting)),
`csc-report` discovers them automatically and renders **both tiers in a
single self-contained HTML** under a top-level "Confidence tier"
selector:

* The §1 executive summary, §3 cohort landscape, §4 variant-calling
  impact, §5 detection summary and §6 per-sample appendix are all
  rendered once per tier.  Both tier sections are present in the DOM;
  the selector toggles their `display` style client-side, so flipping
  between them is instantaneous and does not lose scroll state.
* A new always-visible **§5.1 Sensitive vs High-Confidence
  concordance** subsection summarises the flag-set agreement /
  disagreement between the two tiers, including a 3-cell tile diagram
  (both / sensitive-only / high-confidence-only) and the per-sample
  `reads_demoted_to_unclassified` counts from
  `aggregation_metadata.json`.
* The Methods §2 gains a new **§2.6 Sensitive vs high-confidence
  reporting** subsection that explains the Kraken2 confidence
  recomputation and cites Wood *et al.* 2019, Marcelino *et al.* 2020
  and Lu &amp; Salzberg 2020.
* The Discussion §7 gains an auto-generated bullet list comparing the
  two tiers' top species and flag overlap.

`report_manifest.json` exposes a new `confidence_tiers[]` array with
one entry per tier (`tier_suffix`, `threshold`, source matrix
filenames, per-tier `samples_flagged_primary_detect` /
`samples_flagged_abs_detect`, sidecar TSV filename).  When only the
sensitive tier is present (legacy outputs or
`confidence_thresholds: []`), the array is empty, the toggle is
omitted, and the report is byte-compatible with the schema-2.0 single-
tier rendering — fully backward-compatible.

## What changed in schema 2.0 (cohort layout)

`REPORT_SCHEMA_VERSION` was **`2.0`** before this release (breaking,
additive only).  The default report layout (`--layout cohort`) is
**species-centric** rather than per-sample-centric so cohorts of 3K+
samples remain interpretable in seconds rather than tens of pages of
scrollable tables.  The previous per-sample layout is preserved for
one release window behind `--layout legacy` for byte-level diffing.

| What | Old layout (`legacy`) | New layout (`cohort`, default) |
|---|---|---|
| Unit of analysis | per sample | per **species** |
| §3 contents | per-sample summary table + cohort top-N table + per-sample stacked bars | §3.1 species summary table (paginated, sortable, filterable, with sparklines), §3.2 prevalence × abundance map, §3.3 rank-abundance, §3.4 core/accessory/rare partition, §3.5 cohort-wide distribution figures, §3.6 hierarchical-cluster heatmap, §3.7 PCoA β-diversity ordination, §3.8 per-species drill-down |
| §4 Variant-Calling Impact | linear table of flagged samples | histogram of cohort burden + species attribution stacked bar + paginated flagged-samples table |
| Detection summary | folded inside §3 | promoted to its own §5 with an UpSet-style 2-set tile diagram |
| Per-sample table | full 3K-row literal table | demoted to §6, paginated (default 25/page) + sortable + filterable + a sibling `per_sample_summary.tsv` for offline analysis |
| Discussion | static text | static text + auto-generated species-keyed paragraphs from the §3.4 partition |

### `report_manifest.json` migration

All schema-1.0 keys remain.  New schema-2.0 keys are *additive*:

- `layout` (`"cohort"` or `"legacy"`)
- `page_size`, `top_species`, `cluster_method`, `cluster_distance`
- `prevalence_core`, `prevalence_rare`, `max_samples_cluster`
- `min_reads_for_prevalence`
- `species_summary[]` (top-200 compact JSON of the §3.1 table)
- `partition_counts` (`{core, accessory, rare, core_threshold, rare_threshold}`)
- `top_species_by_burden`, `top_species_by_prevalence`
- `per_sample_tsv` (filename of the §6 sidecar)

If a downstream consumer was already pinning to schema `1.0`, either
ask `csc-report --layout legacy` (manifest reports `"layout":
"legacy"`) or upgrade to ignore the new keys.

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
3. The cohort layout uses **robust statistics** (median, IQR, P95,
   robust CV = MAD / median) over positives only, never bare means –
   contamination distributions are heavy-tailed.
4. A dedicated **Variant-Calling Impact** section uses the absolute
   denominator to flag samples whose non-human content exceeds a
   configurable threshold (default 0.1 % of total sequencing).
5. The Methods section distinguishes the three read counts that the
   pipeline produces: *total sequenced*, *extracted for classification*,
   and *classified by Kraken2*.
6. A **Methods Transparency Checklist** records schema versions,
   parameters, and relative links to every input matrix so the report
   is fully reproducible and review-proof.

The renderer uses only the Python standard library (no Jinja2, pandas,
matplotlib, scipy, or numpy dependencies).  Hierarchical clustering
(Lance–Williams over a precomputed distance matrix) and PCoA (classical
multidimensional scaling via power iteration with deflation) are
implemented in `csc/report/cohort.py`; SVG primitives in
`csc/report/svg.py`; and inline JS / CSS in `csc/report/interactive.py`.
The HTML stays a single self-contained artefact.

## Usage

```bash
# Cohort layout (default).  Single-pass everything.
csc-report aggregate_out/ -o report/contamination_report.html

# With detect outputs (populates §2.4 and §5)
csc-report aggregate_out/ -o report/contamination_report.html \
    --detect-dir detect_out/

# Tighter variant-calling impact threshold (100 ppm = 0.01 % of sequencing)
csc-report aggregate_out/ -o report/contamination_report.html \
    --variant-impact-threshold-ppm 100

# Smaller pages, fewer drill-downs, ward linkage, Jaccard distance.
csc-report aggregate_out/ -o report/contamination_report.html \
    --page-size 50 --drilldown-top 10 \
    --cluster-method ward --cluster-distance jaccard

# Larger cohort: cap pairwise distance computation
csc-report aggregate_out/ -o report/contamination_report.html \
    --max-samples-cluster 1000

# Reproduce the schema-1.0 per-sample layout (one release window only)
csc-report aggregate_out/ -o report/legacy_report.html --layout legacy
```

## How to read the contamination quadrants in 30 seconds (§3.2)

The §3.2 *Prevalence × Abundance* scatter plots one dot per species.

- **x = prevalence** (fraction of samples with ≥1 read for that taxon).
- **y = log10 median absolute burden** (ppm of total sequenced reads,
  positives only) when `taxa_matrix_abs.tsv` is available, or
  log10 median CPM otherwise.
- **size** = log10 cohort raw reads (so dominant species are visually
  larger).
- **colour** = taxonomic domain.

Quadrant guides at `prevalence = 50 %` and the median y-value give:

| Quadrant | Meaning | Typical action |
|---|---|---|
| **Core kitome** (high prev, low burden) | Reagent / library contamination | Add to `csc-detect --kitome-taxa`; exclude from biological interpretation |
| **Systemic burden** (high prev, high burden) | Cohort-wide contamination episode (or a real biological signal) | Investigate batches, reagent lots; cross-reference with §3.6 heatmap |
| **Outburst** (low prev, high burden) | Per-sample handling errors | Inspect listed samples in §3.8 / §6 |
| **Trace** (low prev, low burden) | Ignorable noise | None |

## Required inputs

`<aggregate_dir>/` must contain:

| File | Purpose |
|---|---|
| `taxa_matrix_raw.tsv` | Raw integer read counts |
| `taxa_matrix_cpm.tsv` | CPM (compositional) |
| `aggregation_metadata.json` | Schema version, parameters, provenance |

Optional but **strongly recommended** (enables §4 Variant-Calling
Impact and the absolute-burden figures):

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
| `<output-parent>/report_manifest.json` | Machine-readable manifest with the report schema version, generated timestamp, chosen threshold, the list of samples flagged by §4, the compact species summary, and the partition counts |
| `<output-parent>/per_sample_summary.tsv` | (Cohort layout only) full per-sample table written alongside §6 for offline analysis |

## Report sections (cohort layout)

1. **Executive Summary** – cohort headline metrics: sample / taxon
   counts, median non-human composition and absolute burden, top
   species by cohort burden and by prevalence, partition counts, and a
   mini horizontal stacked bar of cohort burden by domain.
2. **Methods** – extract / classify / aggregate / detect procedures
   and a clear statement distinguishing the three read counts.
3. **Cohort species landscape** – the new primary section.  Eight
   sub-sections (§3.1 species summary table, §3.2 prevalence ×
   abundance map, §3.3 rank-abundance, §3.4 core/accessory/rare
   partition, §3.5 cohort-wide distribution figures, §3.6
   hierarchical-cluster heatmap, §3.7 PCoA, §3.8 per-species
   drill-down).
4. **Variant-Calling Impact** – histogram of cohort burden, species
   attribution stacked bar, and a paginated table of flagged samples.
5. **Detection summary** – relocated from the legacy §3.3.  Adds an
   UpSet-style 2-set tile diagram (primary ∩ abs).
6. **Per-sample appendix** – paginated, sortable, filterable.  Default
   25 rows / page.  A sidecar `per_sample_summary.tsv` is written for
   offline use.
7. **Discussion & Caveats** – kitome / outburst guidance auto-generated
   from §3.4.
8. **Methods Transparency Checklist** – versioned references and
   download links for every input matrix.

## Defensible denominator choices

| Matrix | Denominator | Report section where used |
|---|---|---|
| `taxa_matrix_raw.tsv` | — (integer counts) | Reproducibility table, transparency checklist |
| `taxa_matrix_cpm.tsv` | Classified direct reads | Composition figures, % non-human of classified, diversity indices, species summary CPM columns |
| `taxa_matrix_abs.tsv` | Total sequenced reads (idxstats) | §3.2 prevalence × abundance, §3.3 rank-abundance, §4 Variant-Calling Impact, manuscript QC tables, species summary ppm columns |

See [docs/aggregate.md](aggregate.md) for the full denominator
cheat-sheet and [docs/extract.md](extract.md) for how the total
sequenced read count is captured at extract time.

## When to use `--layout legacy`

The legacy per-sample layout is retained for one release window for:

- **Byte-level diffing** of historical reports against the new layout.
- **Tiny cohorts** (n < ~30) where a literal per-sample stacked-bar
  figure is genuinely more readable than aggregated distributions.
- **Downstream consumers** still pinned to schema 1.0.  In schema 2.0
  manifests the `layout` key reports `"legacy"` so consumers can
  detect the older shape.

For cohorts of 200+ samples the cohort layout is the recommended
default.

## Performance budget

The cohort layout's expensive steps are the heatmap and PCoA, which
need pairwise sample distances – O(n²) memory and time.  The
`--max-samples-cluster` flag (default 2000) caps these computations and
deterministically sub-samples (every-k-th sample) above the cap; the
caption of the resulting figure documents the sub-sampling factor.
For very large cohorts (n > 5000) start with
`--max-samples-cluster 1500` and increase only if needed.

## Python API

```python
from csc.report import load_inputs, generate_html_report

inputs = load_inputs("aggregate_out/", detect_dir="detect_out/")
generate_html_report(
    inputs,
    "report/contamination_report.html",
    threshold_ppm=1000.0,        # 0.1 % of total sequencing
    layout="cohort",              # default
    page_size=25,
    top_species=50,
    cluster_method="average",     # or "ward" / "single"
    cluster_distance="bray",      # or "jaccard"
    prevalence_core=0.5,
    prevalence_rare=0.1,
    max_samples_cluster=2000,
)
```
