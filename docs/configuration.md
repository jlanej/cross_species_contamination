# Configuration

All CSC pipeline modules share a central YAML configuration file located at
`csc/default_config.yaml`.

## Loading Configuration

```python
from csc.config import load_config

# Load defaults only
cfg = load_config()

# Load with user overrides
cfg = load_config("my_config.yaml")
```

The config can also be set via the `CSC_CONFIG` environment variable:

```bash
export CSC_CONFIG=/path/to/my_config.yaml
```

If both an explicit path and the environment variable are set, the explicit
path takes precedence.

## Default Values

```yaml
extract:
  mapq_threshold: null
  threads: 1
  interleaved: false
  output_suffix: ".unmapped"

classify:
  tool: "kraken2"
  db: null
  confidence: 0.0
  threads: 1
  recommended_db: "prackendb"

aggregate:
  min_reads: 10
  rank_filter:
    - "S"
    - "G"
    - "F"

detect:
  method: "mad"
  mad_threshold: 3.5
  iqr_multiplier: 1.5
  subtract_background: true
  kitome_taxa: []

logging:
  level: "INFO"
```

`csc-aggregate` always writes both `taxa_matrix_raw.tsv` (integer direct-read
counts) and `taxa_matrix_cpm.tsv` (counts-per-million) plus typed per-rank
matrices such as `taxa_matrix_raw_S.tsv` and `taxa_matrix_cpm_S.tsv`.
Use `--detect_matrix cpm` (default) or `--detect_matrix raw` in the pipeline
to select which matrix is passed to `csc-detect`.

## Overriding Defaults

Create a YAML file with only the keys you want to change.  Unspecified keys
retain their default values.  The merge is recursive so nested keys are
preserved:

```yaml
# my_config.yaml – only override what you need
extract:
  threads: 8
  mapq_threshold: 10

classify:
  db: /data/kraken2_db

logging:
  level: DEBUG
```

### CLI vs Configuration File Precedence

When using the CSC CLI tools directly (e.g. `csc-aggregate`, `csc-classify`),
**CLI arguments take precedence** over configuration file values.  CLI defaults
may differ from the config file defaults:

| Parameter | CLI Default | Config Default | Notes |
|-----------|-------------|----------------|-------|
| `min_reads` | `0` | `10` | CLI default includes all taxa; config applies a filter |
| `confidence` | `0.0` | `0.0` | Consistent between CLI and config |

The Nextflow pipeline reads values from the configuration file.  When running
CLI tools standalone, pass parameters explicitly to match config behavior
(e.g. `--min-reads 10`).
