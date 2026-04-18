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
(default: 3.5):

    modified_z = |value − median| / (MAD × 1.4826)

The constant 1.4826 rescales MAD to be consistent with the standard
deviation for normally distributed data.  MAD is highly robust to
extreme outliers—up to 50 % of the data can be contaminated without
affecting the estimate.

**Default threshold: 3.5** — corresponds to approximately 3.5 standard
deviations for normally distributed data, offering a good balance between
sensitivity and specificity for metagenomic contamination detection
(Iglewicz and Hoaglin, 1993).

### IQR (Inter-Quartile Range)

A fence-based method.  A value is flagged when it exceeds
`Q3 + k × IQR`, where *k* is the configurable multiplier (default: 1.5).
The IQR method works well with larger sample sizes (≥ 10 recommended).

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

Before applying any method the per-taxon population mean can be
subtracted (enabled by default).  This centres the data around zero so
that detection focuses on *excess* signal rather than absolute abundance.
Negative values after subtraction are preserved to maintain the
distribution shape, which keeps MAD and IQR reliable.

### Kitome / Environmental Background Exclusion

Taxa known to be laboratory or environmental contaminants (so-called
"kitome" taxa) can be excluded by providing their NCBI taxonomy IDs.
This prevents reagent-derived sequences from triggering false positives.

## Configuration

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
# Default: run all three methods together (recommended)
csc-detect results/taxa_matrix_cpm.tsv -o results/detect/

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

- **MAD = 0 handling:** When MAD equals zero (more than half the values
  are identical), no outliers are flagged for that taxon.  This is
  statistically conservative but avoids false positives on taxa with
  no measurable variance.

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
