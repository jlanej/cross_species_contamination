# CSC Pipeline Documentation

This directory contains documentation for the Cross-Species Contamination
(CSC) detection pipeline.

## Modules

| Module | Status | Description |
|--------|--------|-------------|
| [extract](extract.md) | ✅ Implemented | Streaming extraction of unmapped/poorly-mapped reads |
| [classify](classify.md) | ✅ Implemented | Taxonomic classification of extracted reads |
| [aggregate](aggregate.md) | ✅ Implemented | Aggregation of classification results |
| [detect](detect.md) | 🔲 Stub | Statistical contamination detection |

## Configuration

All modules share a central YAML configuration file.  See [configuration.md](configuration.md)
for details on default values and how to override them.

## Testing

See [testing.md](testing.md) for the full testing guide, including how to
generate synthetic fixtures, run golden-output regression tests, and
interpret accuracy metrics (precision, recall, FDR).

## Pipeline Overview

```
BAM/CRAM ──► extract ──► classify ──► aggregate ──► detect
                │            │            │            │
            FASTQ files   taxa labels   summary     report
```
