"""Classify module – taxonomic classification of extracted reads.

This module wraps Kraken2 to assign taxonomic labels to reads extracted
by the ``csc.extract`` module.

Public API
----------
.. autofunction:: classify_reads
.. autofunction:: build_classify_command
.. autofunction:: validate_database
"""

from csc.classify.classify import (
    ClassificationResult,
    build_classify_command,
    classify_reads,
    validate_database,
)

__all__ = [
    "ClassificationResult",
    "build_classify_command",
    "classify_reads",
    "validate_database",
]
