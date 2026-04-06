"""Shared utility helpers for the CSC pipeline.

Common functions used across multiple CSC modules live here so they
can be imported as ``from csc.utils import ...``.
"""

from __future__ import annotations

import logging


def setup_logging(level: str = "INFO") -> None:
    """Configure root logging with a consistent format.

    Parameters
    ----------
    level:
        Log level name (e.g. ``"DEBUG"``, ``"INFO"``).
    """
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )
