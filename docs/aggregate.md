# Aggregate Module

> **Status:** Stub – not yet implemented.

## Planned Functionality

The **aggregate** module (`csc.aggregate`) will collect per-sample
classification outputs and produce summary tables across a cohort.

### Planned Features

- Merge per-sample Kraken2 reports into a single matrix
- Filter taxa by minimum read count
- Normalize counts for cross-sample comparison

## Configuration

```yaml
aggregate:
  min_reads: 10
```
