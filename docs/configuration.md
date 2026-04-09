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

detect:
  method: "statistical"
  fdr: 0.05

logging:
  level: "INFO"
```

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
