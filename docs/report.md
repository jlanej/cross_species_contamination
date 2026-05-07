# Report Module

The `csc-report` CLI generates a static, self-contained HTML summary of
non-human content in a cohort of human whole-genome sequencing (WGS)
samples directly from the outputs of `csc-aggregate` (and optionally
`csc-detect`).

> **AI Acknowledgment:** This module was developed with AI assistance.
> Best practices in bioinformatics should always take precedence over
> specific implementation details.

## What changed in schema 2.2 (renderable for thousand-sample cohorts)

`REPORT_SCHEMA_VERSION` is now **`2.2`** (additive only).  For cohorts
with thousands of samples the previous layout could produce >80 MB
non-responsive HTML files because every sample and every taxon was
materialised as a DOM row / SVG element in the HTML.  Schema 2.2
refocuses the static report on cohort-level summaries and curated
examples, while always writing the full per-sample / per-taxon detail
to TSV sidecars for offline auditability:

* **§6 is now *Notable samples & per-sample sidecar*.**  The literal
  per-sample HTML table has been removed; in its place are short
  highlight tables for the most informative samples (top by absolute
  burden, top by composition, most diverse, and any sample(s) flagged
  by `csc-detect`).  The complete per-sample summary remains written
  to `per_sample_summary.tsv` (linked prominently from the section).
* **§3.1 species summary table** is capped at the top
  `--species-table-top` taxa (default 200) sorted by cohort burden.
  The complete table is always written to
  `species_summary.tsv` (and tier-suffixed siblings) regardless.
* **§4 Variant-Calling Impact flagged-samples table** is capped at the
  top `--variant-flagged-top` rows (default 100) sorted by burden.
  The complete list is always written to
  `variant_impact_flagged.tsv`.
* **§3.6 heatmap** has its sample cap split from the clustering cap
  (`--max-samples-heatmap`, default `min(--max-samples-cluster, 500)`)
  because each cell adds inline-SVG markup; per-cell `<title>`
  tooltips are also automatically suppressed when the grid exceeds
  ~5 000 cells (sub-pixel cells cannot be hovered individually anyway).
* `report_manifest.json` records the new caps
  (`species_table_top`, `variant_flagged_top`, `notable_top`,
  `max_samples_heatmap`) and the new sidecar filenames
  (`species_summary_tsv`, `variant_impact_flagged_tsv`).

These changes are backward-compatible with previous schema-2.x manifest
consumers (additive only) and the legacy per-sample layout remains
available behind `--layout legacy` for one release window.

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
| §3 contents | per-sample summary table + cohort top-N table + per-sample stacked bars | §3.1 species summary table (paginated, sortable, filterable, with sparklines), §3.2 prevalence × abundance map, §3.3 rank-abundance, §3.4 core/accessory/rare partition, §3.5 cohort-wide distribution figures, §3.6 hierarchical-cluster heatmap (with vertical colour-bar legend), §3.6.1 per-sample domain composition (ordered by §3.6 column order), §3.6.2 prevalence × mean-abundance scatter, §3.6.3 cohort burden distribution by domain, §3.6.4 per-sample diversity overview (Shannon + richness histograms), §3.6.5 absolute-burden parallel of §3.6 (when `taxa_matrix_abs.tsv` is available), §3.7 PCoA β-diversity ordination, §3.8 per-species drill-down |
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
| `<detect_dir>/<tier>/abs/qc_summary.json` | Absolute-burden side-pass detect summary; presence of this file is what `csc-report` uses to decide whether to render the §2.5 "absolute-burden side pass" section |

> **Note** — the report decides the absolute-burden side pass is
> "available" by looking for `<detect_dir>/<tier>/abs/qc_summary.json`,
> not by inspecting `aggregate/`.  If you re-use a `detect/` directory
> that pre-dates the absolute-burden matrix, the §2.5 callout will
> warn that the side pass was not run even when
> `taxa_matrix_abs.tsv` is present in `aggregate/`.  Re-running
> `csc-detect` against the current `aggregate/` directory (no
> re-aggregation needed) populates `detect/<tier>/abs/` and resolves
> the warning.

## Outputs

| File | Purpose |
|---|---|
| `<output>.html` | Self-contained HTML report |
| `<output-parent>/report_manifest.json` | Machine-readable manifest with the report schema version, generated timestamp, chosen threshold, the list of samples flagged by §4, the compact species summary, the partition counts, and the names of every TSV sidecar |
| `<output-parent>/per_sample_summary.tsv` | (Cohort layout) **Complete** per-sample table (every sample × every metric).  §6 in the HTML only renders curated highlight tables; this TSV has the full data. |
| `<output-parent>/species_summary.tsv` | (Cohort layout) **Complete** §3.1 species summary (every non-human taxon × every metric).  The HTML §3.1 table is capped at `--species-table-top` for renderability. |
| `<output-parent>/variant_impact_flagged.tsv` | (Cohort layout) **Complete** list of samples flagged by §4 Variant-Calling Impact.  The HTML §4 table is capped at `--variant-flagged-top`. |

## Report sections (cohort layout)

1. **Executive Summary** – cohort headline metrics: sample / taxon
   counts, median non-human composition and absolute burden, top
   species by cohort burden and by prevalence, partition counts, and a
   mini horizontal stacked bar of cohort burden by domain.
2. **Methods** – extract / classify / aggregate / detect procedures
   and a clear statement distinguishing the three read counts.
3. **Cohort species landscape** – the new primary section.  Sub-sections:
   §3.1 species summary table, §3.2 prevalence × abundance map, §3.3
   rank-abundance, §3.4 core/accessory/rare partition, §3.5 cohort-wide
   distribution figures, §3.6 hierarchical-cluster heatmap (with a
   vertical colour-bar legend on the right that maps cell shade to the
   underlying log1p-CPM value, anchored at 0 / 10× / 100× / max), and
   five §3.6.x companion summary plots driven from the same clustered
   grid as §3.6:
   * **§3.6 Sample × species heatmap** — the kraken2 `unclassified`
     (tax_id 0) and `root` (tax_id 1) pseudo-taxa are excluded
     up-front from the top-K species selection.  Their prevalence and
     magnitude would otherwise dominate every column and produce
     uninformative dendrograms.  When confidence-tier sibling matrices
     are present (e.g. `taxa_matrix_*_conf0p10.tsv`), the species set,
     Bray–Curtis distance matrix, and both (sample, taxa) dendrograms
     are recomputed independently for each tier — flipping the tier
     picker can therefore change the row/column order, which is the
     point of the comparison.
   * **§3.6.1 Per-sample domain composition** — a 100%-stacked column
     bar where each column is one sample (ordered identically to §3.6
     columns).  Solid blocks of one colour across a §3.6 cluster
     indicate host-depletion or batch artefacts; a mosaic suggests
     genuine per-sample contamination heterogeneity.  Detect-flagged
     samples are marked with an orange dot.  Unlike §3.6 / §3.6.5 /
     §3.7, the `Unclassified` domain is **kept** here because the size
     of the unclassified bucket is itself a useful QC signal.
   * **§3.6.2 Prevalence × abundance scatter** — top-200 species
     plotted as prevalence (fraction of cohort) × median CPM among
     positives; dot size scales with p95 CPM.  Surfaces rare-but-spiky
     contaminants (large dots in the bottom-left) that the heatmap
     dilutes.  Points are **interactive on hover**: a JS-driven
     tooltip appears immediately with the species name, tax_id,
     domain, prevalence, median CPM and p95 CPM (the same content is
     also available as a native SVG `<title>` for the no-JS
     fallback).  The same hover affordance is wired into §3.2 and
     §3.7.
   * **§3.6.3 Cohort burden distribution by domain** — boxplot of
     per-sample CPM totals by taxonomic domain; quickly answers "is
     this a host-depletion-rate cohort or a contamination cohort?"
     without reading every cell of the heatmap.
   * **§3.6.4 Per-sample diversity overview** — Shannon entropy and
     observed-species histograms, with the count of detect-flagged
     samples annotated.  Often pulls out PCR-duplication /
     undersequenced samples.
   * **§3.6.5 Absolute-burden parallel of §3.6** — same row/column
     order as §3.6 (and therefore the same `unclassified` / `root`
     exclusion) but cells encode log1p(ppm) instead of log1p(CPM)
     (rendered only when `taxa_matrix_abs.tsv` is available).  A
     cluster that is dark in §3.6 but pale here is a host-depletion
     artefact; a cluster that is dark in both is a genuine
     contamination episode.

   Followed by §3.7 PCoA — which applies the *same* `unclassified` +
   `root` exclusion as §3.6 and reports each sample's dominant
   *informative* domain on hover (Human, Unclassified and root are
   skipped) — and §3.8 per-species drill-down.
4. **Variant-Calling Impact** – histogram of cohort burden, species
   attribution stacked bar, and a paginated table of flagged samples.
5. **Detection summary** – relocated from the legacy §3.3.  Adds an
   UpSet-style 2-set tile diagram (primary ∩ abs).
6. **Notable samples & per-sample sidecar** – curated highlight
   tables (top by absolute burden, top by composition, most diverse,
   `csc-detect`-flagged exemplars).  The literal full per-sample table
   is **not** rendered in the HTML for cohorts of thousands of samples
   — the complete per-sample summary is always written to the sidecar
   `per_sample_summary.tsv` and prominently linked from the section.
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
`--max-samples-cluster` flag (default 2000) caps the pairwise distance
work and deterministically sub-samples (every-k-th sample) above the
cap; the caption of the resulting figure documents the sub-sampling
factor.  The `--max-samples-heatmap` flag (default
`min(--max-samples-cluster, 500)`) further caps the heatmap because
each cell is an inline SVG element — for cohorts of 1K+ samples a
tighter heatmap cap keeps the static HTML small without impacting the
PCoA which renders one circle per sample.  When the heatmap grid would
exceed ~5 000 cells the renderer also suppresses per-cell `<title>`
tooltips automatically (sub-pixel cells cannot be hovered
individually anyway).

To control HTML size for thousand-sample cohorts the renderer caps the
*HTML* projections of the long tables at user-controlled top-N values
and writes the complete tables to TSV sidecars:

| Flag | Defaults | Bounds |
|---|---|---|
| `--species-table-top` | 200 | §3.1 HTML rows |
| `--variant-flagged-top` | 100 | §4 flagged-samples HTML rows |
| `--notable-top` | 15 | rows per §6 highlight table |
| `--max-samples-heatmap` | min(2000, 500) | §3.6 heatmap columns |

For very large cohorts (n > 5000) start with
`--max-samples-cluster 1500 --max-samples-heatmap 400` and increase
only if needed.  Empirically, a synthetic 1500-sample × 400-taxon
cohort renders to ~2.4 MB HTML under the schema-2.2 defaults compared
to 80+ MB under the unconstrained schema-2.1 layout.

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
    max_samples_heatmap=500,
    species_table_top=200,        # §3.1 HTML rows; full table → species_summary.tsv
    variant_flagged_top=100,      # §4 HTML rows; full list → variant_impact_flagged.tsv
    notable_top=15,               # rows per §6 highlight table
)
```
