# Testing Guide

This document describes the test strategy, how to run tests, how to
generate and update test fixtures, and how to interpret the outputs.

## Quick Start

```bash
# Install test dependencies
pip install pysam pytest-cov
pip install ".[test]"

# Run the full suite
pytest -v tests/

# Run with coverage
pytest -v --cov=csc --cov-report=term-missing tests/
```

## Test Organisation

| Test file | Module(s) | Description |
|-----------|-----------|-------------|
| `test_extract.py` | `csc.extract` | BAM/CRAM extraction — input validation, command building, integration |
| `test_classify.py` | `csc.classify` | Kraken2 classification — DB validation, command building, mocked runs |
| `test_db.py` | `csc.classify.db` | DB management — hash, fetch, cache, clean |
| `test_config.py` | `csc.config` | YAML configuration loading, merging, env-var override |
| `test_aggregate.py` | `csc.aggregate` | Report parsing, matrix building, CPM normalisation, edge cases |
| `test_detect.py` | `csc.detect` | Statistical helpers, MAD/IQR outlier detection, precision/recall |
| `test_golden.py` | `csc.aggregate`, `csc.detect` | Golden-output regression tests (exact output matching) |
| `test_end_to_end.py` | `csc.aggregate` → `csc.detect` | End-to-end pipeline with spike-in controls, accuracy metrics |
| `test_nextflow_pipeline.py` | Nextflow | Pipeline file structure, DSL2 syntax, process definitions |

## Synthetic Test Data

### BAM/CRAM Fixtures

Synthetic BAM and CRAM files are generated at test time by
`tests/generate_test_data.py`.  They contain a deterministic mixture of:

| Read type | Count (default) | MAPQ | Purpose |
|-----------|----------------|------|---------|
| Well-mapped pairs | 20 | 60 | Reads that should **not** be extracted |
| Unmapped pairs | 10 | 0 | Simulated contaminant reads |
| Low-MAPQ pairs | 5 | 3 | Borderline reads for MAPQ-threshold testing |

**Seed:** `random.seed(42)` for reproducibility.

Reference FASTA (`reference.fa`) contains two 500 bp chromosomes (`chr1`,
`chr2`) and is automatically indexed.

### Kraken2 Report Fixtures

Each test creates Kraken2 report files in `tmp_path`.  These are plain-text
files with the standard Kraken2 six-column format and deterministic read
counts.

### Generating Test Data Manually

```bash
# BAM/CRAM/reference
python tests/generate_test_data.py /tmp/test_data

# Golden output files
python tests/generate_golden.py
```

## Golden-Output Regression Testing

Golden reference files are stored in `tests/golden/`:

| File | Module | What it captures |
|------|--------|------------------|
| `aggregate_matrix.tsv` | aggregate | Expected raw-count matrix for a two-sample fixture |
| `aggregate_metadata.json` | aggregate | Expected aggregation metadata |
| `detect_flagged_samples.tsv` | detect | Expected flagged-sample table (MAD method) |
| `detect_qc_summary.json` | detect | Expected QC summary JSON |
| `detect_quarantine_list.txt` | detect | Expected quarantine list |

If a golden test fails it means pipeline behaviour has changed.  Review
the diff carefully:

1. **Intentional change** — regenerate golden files:
   ```bash
   python tests/generate_golden.py
   ```
   Commit the updated golden files with a description of what changed.

2. **Unintentional regression** — investigate the code change that caused
   the output to diverge and fix the bug.

## Negative and Positive Controls

### Negative Controls (Clean Cohort)

`TestNegativeControlCleanCohort` in `test_end_to_end.py` creates a cohort
of 10 samples with uniform read counts across all taxa.  The expected
result is **zero** flags from both MAD and IQR methods.

### Positive Controls (Contaminated Cohort)

`TestPositiveControlContaminatedCohort` in `test_end_to_end.py` creates a
cohort of 15 clean + 3 contaminated samples.  Contaminated samples have a
spike-in of 5 000 reads for *E. coli* (taxon 562).  The tests validate:

* **Precision ≥ 0.8** — most flagged samples are truly contaminated
* **Recall ≥ 0.8** — most contaminated samples are detected
* **FDR ≤ 0.2** — false discovery rate is controlled

## Accuracy Metrics

The test suite continuously monitors:

| Metric | Threshold | Test |
|--------|-----------|------|
| Precision | ≥ 0.8 | `test_precision_above_threshold` |
| Recall | ≥ 0.8 | `test_recall_above_threshold` |
| False discovery rate | ≤ 0.2 | `test_fdr_below_threshold` |

These thresholds apply to the MAD method on well-separated synthetic data.
If detection parameters change, review the impact on these metrics.

## CI Configuration

The GitHub Actions workflow (`.github/workflows/ci.yml`) runs:

* **Python matrix:** 3.9, 3.10, 3.11, 3.12
* **Coverage check:** minimum 60 % line coverage via `pytest-cov`
* **samtools:** installed via `apt-get`
* **pysam:** installed via `pip`

## Interpreting Outputs

### Flagged Samples TSV

The `flagged_samples.tsv` output has columns:

| Column | Description |
|--------|-------------|
| `sample_id` | Sample that triggered the flag |
| `tax_id` | NCBI taxonomy ID of the anomalous taxon |
| `taxon_name` | Human-readable taxon name |
| `value` | Observed value (after background subtraction if enabled) |
| `population_median` | Median value across all samples for this taxon |
| `deviation` | Modified Z-score (MAD) or distance above fence (IQR) |
| `method` | Detection method used |

### QC Summary JSON

Key fields in `qc_summary.json`:

| Field | Description |
|-------|-------------|
| `total_samples` | Number of samples in the matrix |
| `total_taxa_analysed` | Taxa remaining after kitome filtering |
| `method` | Detection method |
| `flagged_count` | Number of (sample, taxon) flags |
| `flagged_samples` | Unique sample IDs that were flagged |
| `kitome_taxa_excluded` | Count of excluded kitome taxa |
| `background_subtracted` | Whether population-mean subtraction was applied |

### Quarantine List

`quarantine_list.txt` contains one sample ID per line, sorted
alphabetically.  These samples are recommended for further investigation
or exclusion from downstream analyses.
