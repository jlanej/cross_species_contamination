# Extract Module

## Overview

The **extract** module (`csc.extract`) provides streaming extraction of
unmapped and optionally poorly mapped reads from BAM/CRAM files using
samtools.  No intermediate disk files are created.

## CLI Usage

```bash
# Extract unmapped reads
csc-extract sample.bam -o output_dir/

# Include poorly mapped reads (MAPQ < 10)
csc-extract sample.bam -o output_dir/ --mapq 10

# CRAM input with reference
csc-extract sample.cram -o output_dir/ --reference GRCh38.fa

# Interleaved output
csc-extract sample.bam -o output_dir/ --interleaved

# Multiple threads
csc-extract sample.bam -o output_dir/ --threads 4
```

## Python API

```python
from csc.extract import extract_reads, build_extract_command

# Extract reads programmatically
outputs = extract_reads("sample.bam", "output_dir/", mapq_threshold=10)

# Build samtools commands without executing
cmds = build_extract_command(Path("sample.bam"), mapq_threshold=10)
```

## Output Files

| File | Description |
|------|-------------|
| `{sample}.unmapped.R1.fastq.gz` | Read 1 of unmapped pairs |
| `{sample}.unmapped.R2.fastq.gz` | Read 2 of unmapped pairs |
| `{sample}.unmapped.singleton.fastq.gz` | Singleton reads |
| `{sample}.unmapped.other.fastq.gz` | Other reads |
| `{sample}.unmapped.fastq.gz` | Interleaved output (if `--interleaved`) |

## Configuration

Relevant keys in the config YAML:

```yaml
extract:
  mapq_threshold: null  # null = unmapped only; integer = also low-MAPQ
  threads: 1
  interleaved: false
```
