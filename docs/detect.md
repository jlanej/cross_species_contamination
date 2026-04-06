# Detect Module

> **Status:** Stub – not yet implemented.

## Planned Functionality

The **detect** module (`csc.detect`) will apply statistical methods to
aggregated classification results to flag samples with significant
cross-species contamination.

### Planned Features

- Statistical test for each taxon vs. background
- FDR correction for multiple testing
- Generate final contamination report

## Configuration

```yaml
detect:
  method: "statistical"
  fdr: 0.05
```
