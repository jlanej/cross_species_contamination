"""Classify module – taxonomic classification of extracted reads.

This module wraps Kraken2 to assign taxonomic labels to reads extracted
by the ``csc.extract`` module.

Public API
----------
.. autofunction:: classify_reads
.. autofunction:: build_classify_command
.. autofunction:: validate_database
.. autofunction:: fetch_database
.. autofunction:: database_info
.. autofunction:: list_databases
"""

from csc.classify.classify import (
    ClassificationResult,
    build_classify_command,
    classify_reads,
    validate_database,
)
from csc.classify.db import (
    clean_cache,
    compute_hash,
    database_info,
    fetch_database,
    get_cache_dir,
    list_databases,
    verify_hash,
)

__all__ = [
    "ClassificationResult",
    "build_classify_command",
    "classify_reads",
    "clean_cache",
    "compute_hash",
    "database_info",
    "fetch_database",
    "get_cache_dir",
    "list_databases",
    "validate_database",
    "verify_hash",
]
