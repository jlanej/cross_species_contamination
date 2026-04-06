# Classify Module

> **Status:** Stub – not yet implemented.

## Planned Functionality

The **classify** module (`csc.classify`) will provide taxonomic
classification of reads extracted by the extract module.

### Planned Features

- Wrapper around Kraken2 for taxonomic classification
- Support for custom databases
- Configurable confidence thresholds
- Output in standardized format for downstream aggregation

## Configuration

```yaml
classify:
  tool: "kraken2"
  db: null
  confidence: 0.0
  threads: 1
```
