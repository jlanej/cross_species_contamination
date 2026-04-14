"""Classify module – taxonomic classification of extracted reads.

This module wraps Kraken2 to assign taxonomic labels to reads extracted
by the ``csc.extract`` module.

PrackenDB is the recommended Kraken2 database — it contains one genome
per species, enabling unambiguous species-level detection.  Use
:func:`fetch_prackendb` to download and validate PrackenDB.

Public API
----------
.. autofunction:: classify_reads
.. autofunction:: build_classify_command
.. autofunction:: validate_database
.. autofunction:: fetch_database
.. autofunction:: fetch_prackendb
.. autofunction:: validate_taxonomy
.. autofunction:: is_prackendb
.. autofunction:: database_info
.. autofunction:: list_databases
.. autofunction:: estimate_db_memory
"""

from csc.classify.classify import (
    ClassificationResult,
    build_classify_command,
    classify_reads,
    validate_database,
)
from csc.classify.db import (
    PRACKENDB_NAME,
    PRACKENDB_URL,
    TAXONOMY_FILES,
    clean_cache,
    compute_hash,
    database_info,
    estimate_db_memory,
    fetch_database,
    fetch_prackendb,
    get_cache_dir,
    is_prackendb,
    list_databases,
    validate_taxonomy,
    verify_hash,
)

__all__ = [
    "ClassificationResult",
    "PRACKENDB_NAME",
    "PRACKENDB_URL",
    "TAXONOMY_FILES",
    "build_classify_command",
    "classify_reads",
    "clean_cache",
    "compute_hash",
    "database_info",
    "estimate_db_memory",
    "fetch_database",
    "fetch_prackendb",
    "get_cache_dir",
    "is_prackendb",
    "list_databases",
    "validate_database",
    "validate_taxonomy",
    "verify_hash",
]
