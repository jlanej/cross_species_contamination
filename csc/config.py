"""Central configuration loader for the CSC pipeline.

Reads the bundled ``default_config.yaml`` and optionally merges in a
user-provided YAML file.  The user config path can be supplied via:

1. The ``CSC_CONFIG`` environment variable, or
2. An explicit *path* argument to :func:`load_config`.

The merge is a shallow-per-section update: top-level keys in the user
config replace (not deep-merge) the corresponding default sections.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import yaml

_DEFAULT_CONFIG_PATH = Path(__file__).parent / "default_config.yaml"


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge *override* into *base*, returning a new dict."""
    merged = dict(base)
    for key, value in override.items():
        if (
            key in merged
            and isinstance(merged[key], dict)
            and isinstance(value, dict)
        ):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """Load and return the pipeline configuration dictionary.

    Parameters
    ----------
    path:
        Optional path to a user YAML config file.  If *None*, the
        ``CSC_CONFIG`` environment variable is checked.  If neither is
        set, only the built-in defaults are returned.

    Returns
    -------
    dict
        Merged configuration dictionary.
    """
    with open(_DEFAULT_CONFIG_PATH) as fh:
        config: dict[str, Any] = yaml.safe_load(fh)

    user_path = path or os.environ.get("CSC_CONFIG")
    if user_path is not None:
        user_path = Path(user_path)
        if not user_path.is_file():
            raise FileNotFoundError(f"Config file not found: {user_path}")
        with open(user_path) as fh:
            user_config: dict[str, Any] = yaml.safe_load(fh) or {}
        config = _deep_merge(config, user_config)

    return config
