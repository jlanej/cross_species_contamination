# Detect Module

> **Status:** Implemented

## Overview

The **detect** module (`csc.detect`) applies statistical outlier detection
to aggregated classification results (the sample-by-taxon matrix produced by
`csc-aggregate`) to flag samples with significant cross-species contamination.

## Methods

### MAD (Median Absolute Deviation)

The default method.  For each taxon, the population median and MAD are
computed.  A sample is flagged if its modified Z-score exceeds the
configured threshold (default: 3.5):

    modified_z = |value − median| / (MAD × 1.4826)

The constant 1.4826 rescales MAD to be consistent with the standard
deviation for normally distributed data.  MAD is highly robust to
extreme outliers—up to 50 % of the data can be contaminated without
affecting the estimate.

### IQR (Inter-Quartile Range)

An alternative fence-based method.  A value is flagged when it exceeds
`Q3 + k × IQR`, where *k* is the configurable multiplier (default: 1.5).
The IQR method works well with larger sample sizes (≥ 10 recommended).

### Population-Mean Subtraction

Before applying either method the per-taxon population mean can be
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
  method: "mad"              # "mad" or "iqr"
  mad_threshold: 3.5         # scaled-MAD units above median to flag
  iqr_multiplier: 1.5        # IQR fence multiplier
  subtract_background: true  # subtract population-mean per taxon
  kitome_taxa: []            # NCBI taxonomy IDs to exclude
```

## CLI Usage

```bash
# MAD-based detection on CPM matrix (default recommendation)
csc-detect results/taxa_matrix_cpm.tsv -o results/detect/

# IQR method with custom multiplier
csc-detect results/taxa_matrix_cpm.tsv -o results/detect/ --method iqr --iqr-multiplier 2.0

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
| `deviation`        | Modified Z-score (MAD) or distance above fence (IQR)  |
| `method`           | Detection method used                                 |

## Python API

```python
from csc.detect import detect_outliers, generate_report

result = detect_outliers(
    "results/taxa_matrix_cpm.tsv",
    method="mad",
    mad_threshold=3.5,
    kitome_taxa=[9606],
)

reports = generate_report(result, "results/detect/")
```

## Caveats and Recommendations

- **Minimum sample size:** At least 10 samples are recommended for
  robust outlier detection.  With fewer samples (e.g., 5), a single
  extreme outlier can corrupt IQR estimates.  MAD is more resilient but
  still benefits from larger cohorts.

- **MAD = 0 handling:** When MAD equals zero (more than half the values
  are identical), no outliers are flagged for that taxon.  This is
  statistically conservative but avoids false positives on taxa with
  no measurable variance.

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
