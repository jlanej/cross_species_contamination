# Detect Module

> **Status:** Implemented

## Overview

The **detect** module (`csc.detect`) applies statistical outlier detection
to aggregated classification results (the sample-by-taxon matrix produced by
`csc-aggregate`) to flag samples with significant cross-species contamination.

Three complementary detection methods are available and are run together
by default (`--method all`):

| Method | Strengths | Assumption |
|--------|-----------|------------|
| **MAD** | Highly robust to extreme outliers; works with small cohorts | Majority of samples are clean |
| **IQR** | Simple fence-based detection; well-characterised thresholds | Majority of samples are clean |
| **GMM** | Handles widespread contamination and batch effects | *None* — explicitly models both clean and contaminated populations |

Running all three together (the default) ensures comprehensive coverage:
MAD and IQR catch sparse outliers quickly, while the GMM detects
contamination even when a majority of samples are affected.

## Methods

### MAD (Median Absolute Deviation)

For each taxon, the population median and MAD are computed.  A sample
is flagged if its modified Z-score exceeds the configured threshold
(default: 3.5) **and** its value lies above the population median:

    modified_z = |value − median| / (MAD × 1.4826)

The constant 1.4826 rescales MAD to be consistent with the standard
deviation for normally distributed data.  MAD is highly robust to
extreme outliers—up to 50 % of the data can be contaminated without
affecting the estimate.

Only the **upper tail** is flagged because contamination produces
*excess* abundance for a foreign taxon; a sample with an unusually low
abundance for a given taxon is not informative for cross-species
contamination detection.

**Default threshold: 3.5** — corresponds to approximately 3.5 standard
deviations for normally distributed data, offering a good balance between
sensitivity and specificity for metagenomic contamination detection
(Iglewicz and Hoaglin, 1993).

### IQR (Inter-Quartile Range)

A fence-based method.  A value is flagged when it exceeds the upper
fence `Q3 + k × IQR`, where *k* is the configurable multiplier
(default: 1.5).  As with MAD, only the upper tail is evaluated; the
lower fence ``Q1 − k × IQR`` is intentionally not applied because
samples falling below it indicate undersampling rather than
contamination.  The IQR method works well with larger sample sizes
(≥ 10 recommended).

**Default multiplier: 1.5** — the classic Tukey fence threshold that
identifies "mild" outliers.  Increase to 3.0 to restrict flagging to
only "extreme" outliers.

### GMM (Gaussian Mixture Model)

A two-component Gaussian mixture model fitted via Expectation-Maximization
(EM).  Unlike MAD and IQR, the GMM does **not** assume that the majority
of samples are clean.  It explicitly models a *background* component and
a *contamination* component, then flags samples whose posterior
probability of belonging to the high-abundance component exceeds the
threshold.

Three safeguards prevent false positives:

1. **BIC model selection** — a two-component model is only accepted when
   the Bayesian Information Criterion (BIC) prefers it over a single
   Gaussian.  This naturally penalises unnecessary complexity and prevents
   over-fitting when the data is unimodal.
2. **Minimum effective component size** — the contamination component must
   have ≥ 2 effective members (`sum(posteriors) ≥ 2`).  This prevents the
   GMM from treating a single outlier as its own component.
3. **Separation criterion** — the component means must be at least 3× the
   larger component sigma apart (Cohen's d ≥ 3), confirming genuine
   bimodality rather than noise.

**Default threshold: 0.5** — the natural decision boundary of a
two-component mixture.  A sample is flagged when it is more likely to
belong to the contaminated component than the clean component.

**When GMM adds the most value:**  When contamination or batch effects
are widespread (affecting ≥ 50 % of samples), MAD and IQR will treat
the contaminated majority as "normal" and miss the contamination
entirely.  The GMM correctly identifies both populations regardless
of their relative sizes.

### Combined Mode (`all`)

The default `--method all` runs MAD, IQR, and GMM together on each
taxon and merges their flags.  Each flagged (sample, taxon) entry
records which method produced it in the `method` column.  Duplicate
flags from different methods are preserved so that downstream review
can assess concordance.

### Population-Mean Subtraction

Before applying any method the per-taxon **arithmetic mean** across
samples can be subtracted (enabled by default).  This centres the data
around zero so that detection focuses on *excess* signal rather than
absolute abundance.  Negative values after subtraction are preserved
to keep the distribution shape intact.

> **Note:** MAD and IQR are translation-invariant, so subtracting a
> per-taxon constant does **not** change which samples are flagged by
> those methods — it only re-centres the displayed values for
> readability.  The GMM is *not* translation-invariant; for GMM the
> centring re-anchors the background component near zero, which can
> stabilise the EM initialisation.

### Kitome / Environmental Background Exclusion

Taxa known to be laboratory or environmental contaminants (so-called
"kitome" taxa) can be excluded by providing their NCBI taxonomy IDs.
This prevents reagent-derived sequences from triggering false positives.

## Input matrix types

`csc-aggregate` writes three taxa-by-sample matrices that differ only
in their **denominator**, and the choice of which one to feed to
`csc-detect` determines what kind of contamination it can see.  All
three are written to the aggregate output directory side-by-side:

| Matrix | Filename | Denominator | What an outlier means |
|--------|----------|-------------|-----------------------|
| **Raw counts** | `taxa_matrix_raw.tsv` | None (integer reads) | A sample has unusually many reads assigned to a taxon — but library size effects (sequencing depth, host-depletion efficiency) confound interpretation. |
| **CPM** (compositional) | `taxa_matrix_cpm.tsv` | Sum of *classified direct reads* (pre-`min_reads`-filter) per sample | A sample has an unusually large *fraction* of its non-human content from a taxon.  Robust to differences in sequencing depth and to differences in *how much* host material was depleted; sensitive to changes in *composition*. |
| **Absolute burden** | `taxa_matrix_abs.tsv` | *Total sequenced reads* per sample (from `samtools idxstats`, captured at extract time) | A sample has an unusually large *fraction of its whole sequencing run* from a taxon.  Robust to host-depletion-rate differences; sensitive to absolute contamination magnitude. |

> The `abs` matrix is only emitted when `csc-aggregate` is invoked
> with `--idxstats` and per-sample `reads_summary.json` sidecars.
> Without those sidecars, only `raw` and `cpm` are written.

### Why the default is CPM

The CPM matrix is **compositional** and therefore invariant to
sequencing depth and host-depletion efficiency.  It directly answers
the question "what fraction of the non-human content in this sample
is *Salmonella*?" and is the appropriate denominator when a reviewer
asks "which samples have an unusual contaminant *profile*?".  This is
the most common QC question in cohort studies and is the reason
`csc-detect` ships with `taxa_matrix_cpm.tsv` as the documented
input.

### Why CPM alone is not enough — the abs side pass

CPM has a critical blind spot: because it is compositional, a sample
that retains *much* more host than the cohort norm (e.g. due to a
failed depletion step) can carry orders of magnitude more
contaminating reads than its peers and still look perfectly normal in
CPM space, because the *shape* of its non-human content is unchanged.

To plug this gap, `csc-detect` automatically runs a parallel
**absolute-burden side pass** on `taxa_matrix_abs.tsv` whenever:

1. the input matrix is the canonical CPM or raw matrix
   (`taxa_matrix_cpm.tsv` / `taxa_matrix_raw.tsv`); and
2. a sibling `taxa_matrix_abs.tsv` exists in the same directory.

Outputs are written to `<output>/abs/` (and to
`<output>/abs/<rank>/`, `<output>/abs/<tier>/`, etc. when
`--rank-filter` and confidence tiers are also in play), mirroring the
primary directory structure exactly.  Each abs `qc_summary.json`
records `"matrix_type": "abs"` so the HTML report can render the two
flag sets side-by-side and highlight samples that only one pass
catches.

The cost is a single extra detection pass over the same cohort —
trivially cheap relative to extraction or classification — and the
benefit is a contamination call that is independent of host-depletion
variation.  Disable with `--no-abs-detection` if you need the
compositional view alone.



```yaml
detect:
  method: "all"              # "all", "mad", "iqr", or "gmm"
  mad_threshold: 3.5         # scaled-MAD units above median to flag
  iqr_multiplier: 1.5        # IQR fence multiplier
  gmm_threshold: 0.5         # posterior probability threshold for GMM
  subtract_background: true  # subtract population-mean per taxon
  kitome_taxa: []            # NCBI taxonomy IDs to exclude
```

## CLI Usage

```bash
# Default: run all three methods together (recommended).  When a
# sibling taxa_matrix_abs.tsv exists in the same directory the
# absolute-burden side pass is run automatically; results land in
# results/detect/abs/ alongside the primary outputs.
csc-detect results/taxa_matrix_cpm.tsv -o results/detect/

# Disable the absolute-burden side pass (compositional CPM only)
csc-detect results/taxa_matrix_cpm.tsv -o results/detect/ --no-abs-detection

# Run detection directly on the absolute-burden matrix instead of
# CPM.  No further side pass is triggered (the input is already abs).
csc-detect results/taxa_matrix_abs.tsv -o results/detect_abs/

# MAD only
csc-detect results/taxa_matrix_cpm.tsv -o results/detect/ --method mad

# IQR method with custom multiplier
csc-detect results/taxa_matrix_cpm.tsv -o results/detect/ --method iqr --iqr-multiplier 2.0

# GMM only with custom posterior threshold
csc-detect results/taxa_matrix_cpm.tsv -o results/detect/ --method gmm --gmm-threshold 0.8

# Exclude known kitome taxa
csc-detect results/taxa_matrix_cpm.tsv -o results/detect/ --kitome-taxa 9606 562

# Skip population-background subtraction
csc-detect results/taxa_matrix_cpm.tsv -o results/detect/ --no-subtract-background

# Run detection on species-only matrix (auto-discovers matching rank matrix names)
csc-detect results/taxa_matrix_cpm.tsv -o results/detect/ --rank-filter S

# Run detection on species and genus matrices
csc-detect results/taxa_matrix_raw.tsv -o results/detect/ --rank-filter S G F

# Verbose + JSON logging
csc-detect results/taxa_matrix_cpm.tsv -o results/detect/ -v --json-log
```

When `--rank-filter` is specified (default: S G F), the tool looks for
per-rank matrices in the same directory as the input matrix. For example,
`taxa_matrix_cpm.tsv` maps to `taxa_matrix_cpm_S.tsv`;
`taxa_matrix_raw.tsv` maps to `taxa_matrix_raw_S.tsv`.
Results for each rank are written to subdirectories
(e.g. `results/detect/S/`, `results/detect/G/`).

## Output Files

| File                    | Format | Description                                      |
|-------------------------|--------|--------------------------------------------------|
| `flagged_samples.tsv`   | TSV    | One row per (sample, taxon) outlier flag          |
| `qc_summary.json`       | JSON   | Aggregate QC statistics and detection parameters  |
| `quarantine_list.txt`   | Text   | Plain list of sample IDs recommended for review   |
| `<rank>/…`              | dir    | Per-rank-filtered results (when `--rank-filter` is active) |
| `conf<T>/…`             | dir    | Per-confidence-tier results (auto-discovered).  Produced by **default** since v0.3 — `csc-aggregate` ships with `confidence_thresholds: [0.1]`, so a `conf0p10/` subdirectory is normally present alongside the sensitive-tier outputs.  The directory contains the same `flagged_samples.tsv` / `qc_summary.json` triple as the primary pass, plus any `<rank>/` or `abs/` sub-directories.  See [aggregate.md](aggregate.md#high-confidence-tier-dual-tier-reporting) for the rationale. |
| `abs/…`                 | dir    | Absolute-burden side-pass results, with the same `flagged_samples.tsv` / `qc_summary.json` / `quarantine_list.txt` triple plus any `abs/<rank>/`, `abs/conf<T>/` sub-directories.  `qc_summary.json` records `"matrix_type": "abs"` to disambiguate from the primary pass. |

### flagged_samples.tsv columns

| Column             | Description                                           |
|--------------------|-------------------------------------------------------|
| `sample_id`        | Sample identifier                                     |
| `tax_id`           | NCBI taxonomy ID of the flagged taxon                 |
| `taxon_name`       | Scientific name                                       |
| `value`            | Observed value (after background subtraction if used) |
| `population_median`| Median of the taxon across all samples                |
| `deviation`        | Modified Z-score (MAD), distance above fence (IQR), or posterior probability (GMM) |
| `method`           | Detection method that produced this flag (`mad`, `iqr`, or `gmm`) |

### qc_summary.json fields

In addition to overall counts (`total_samples`, `total_taxa_analysed`,
`flagged_count`, `flagged_samples`, `kitome_taxa_excluded`,
`background_subtracted`), `qc_summary.json` also reports diagnostic
counters when the relevant method is run:

| Field                     | Description                                                              |
|---------------------------|--------------------------------------------------------------------------|
| `matrix_type`             | Which input matrix this summary corresponds to: `raw`, `cpm`, or `abs` (omitted if the input filename does not match the canonical typed pattern). |
| `taxa_skipped_mad_zero`   | Taxa skipped by the MAD method because MAD = 0 (≥ 50 % identical values) |
| `taxa_skipped_iqr_zero`   | Taxa skipped by the IQR method because IQR = 0                            |

These counters help users interpret apparently low flag counts: a
detection run on a sparse matrix will frequently see most taxa
skipped, and that information is now visible in the QC output rather
than silently absorbed.

## Python API

```python
from csc.detect import detect_outliers, generate_report

# Run all methods (default, recommended)
result = detect_outliers(
    "results/taxa_matrix_cpm.tsv",
    kitome_taxa=[9606],
)

# Run a specific method
result = detect_outliers(
    "results/taxa_matrix_cpm.tsv",
    method="gmm",
    gmm_threshold=0.5,
    kitome_taxa=[9606],
)

reports = generate_report(result, "results/detect/")
```

## Default Parameter Rationale

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `method` | `all` | Running MAD + IQR + GMM together provides complementary coverage.  MAD/IQR catch sparse outliers; GMM handles widespread contamination and batch effects. |
| `mad_threshold` | `3.5` | Approximately 3.5 standard deviations for normal data.  Standard robust-statistics threshold balancing sensitivity and specificity (Iglewicz & Hoaglin, 1993). |
| `iqr_multiplier` | `1.5` | Tukey's classic "mild outlier" fence.  Well-characterised and widely used in exploratory data analysis. |
| `gmm_threshold` | `0.5` | Natural decision boundary of a two-component mixture: a sample is flagged when it is more likely contaminated than clean. |
| `subtract_background` | `true` | Centres each taxon's distribution around zero so detection focuses on excess signal rather than absolute baseline abundance. |

## Caveats and Recommendations

- **Minimum sample size:** At least 10 samples are recommended for
  robust outlier detection.  With fewer samples (e.g., 5), a single
  extreme outlier can corrupt IQR estimates.  MAD is more resilient but
  still benefits from larger cohorts.  The GMM requires at least 2
  contaminated samples to form a reliable contamination component.

- **MAD / IQR = 0 handling:** When MAD equals zero (more than half the
  values are identical) the modified Z-score is undefined; similarly
  when IQR equals zero the upper fence collapses onto Q3.  In both
  cases the affected taxon is skipped (no flags raised) and counted in
  `qc_summary.json` under `taxa_skipped_mad_zero` /
  `taxa_skipped_iqr_zero`.  This is statistically conservative but
  avoids false positives on taxa with no measurable variance.

- **GMM model selection:** The GMM uses BIC (Bayesian Information
  Criterion) to decide whether two components are warranted.  On
  unimodal data or very small samples, BIC will prefer a single
  Gaussian and no GMM flags will be raised for that taxon—MAD/IQR
  then serve as the primary detectors.

- **Widespread contamination:** If you suspect that contamination or
  batch effects affect a majority of samples, the GMM is the
  appropriate method.  MAD and IQR assume the majority are clean and
  may fail in this scenario.

- **Background subtraction:** Enabled by default.  This helps
  discriminate true contamination from population-level baseline
  abundance.  Disable it (`--no-subtract-background`) when the
  population itself is heterogeneous or when working with very few
  samples.

- **Kitome lists:** Maintain a curated list of known environmental /
  reagent-derived taxa for your laboratory.  Common candidates include
  *Homo sapiens* (9606) and reagent-associated organisms.

- **Downstream review:** Flagged samples should be reviewed manually
  or with additional targeted analyses before any data is discarded.
  The quarantine list is an *alert*, not an automatic exclusion.
