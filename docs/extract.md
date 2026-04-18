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
| `{sample}.idxstats.tsv` | Raw `samtools idxstats` per-chromosome counts |
| `{sample}.reads_summary.json` | Aggregate mapped / unmapped / total reads (see below) |

### `{sample}.idxstats.tsv`

Produced by running `samtools idxstats` on the input BAM/CRAM.  The
command reads only the BAM index (milliseconds per sample) and produces
one TSV row per reference contig plus a final `*` row capturing reads
unmapped **and** without coordinates (mate also unmapped).  Columns:

```
chrom<TAB>length<TAB>mapped<TAB>unmapped
```

### `{sample}.reads_summary.json`

Structured summary parsed from the idxstats TSV.  Schema (versioned via
`schema_version`):

```json
{
  "schema_version": "1.0",
  "sample_id": "NA12878",
  "input": "/abs/path/NA12878.bam",
  "extraction_time": "2024-01-01T12:00:00+00:00",
  "total_mapped": 900123456,
  "total_unmapped": 12345678,
  "total_reads": 912469134,
  "per_chromosome": [
    {"chrom": "chr1", "length": 248956422, "mapped": 70123456, "unmapped": 0},
    {"chrom": "*",    "length": 0,         "mapped": 0,        "unmapped": 12345678}
  ]
}
```

`total_reads` is used by `csc-aggregate` (see
[docs/aggregate.md](aggregate.md)) as the denominator for the
**absolute contamination burden** matrix.

> **Interpretation note – denominator provenance.**  `total_unmapped`
> counts reads in the source BAM with the unmapped flag set.  The
> number of reads that ultimately reach the taxonomic classifier may be
> smaller (adapter trimming, host-read removal, `--mapq` filter) or
> larger (if poorly-mapped reads are routed through classification).
> Always use `reads_summary.json` (captured **at extract time**, before
> any downstream filtering) as the authoritative denominator for
> absolute-burden calculations and manuscript reporting.

### Batch summary TSV (`--summary`)

When `--summary summary.tsv` is passed, the CLI writes a cohort-level
TSV containing per-sample `read_count`, `total_mapped`,
`total_unmapped` and `total_reads` columns – enough to compute
cohort-level absolute burden without re-reading every
`reads_summary.json`.

## Configuration

Relevant keys in the config YAML:

```yaml
extract:
  mapq_threshold: null  # null = unmapped only; integer = also low-MAPQ
  threads: 1
  interleaved: false
```
