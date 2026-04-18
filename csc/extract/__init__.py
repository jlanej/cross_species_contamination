"""Extraction module – streaming unmapped/poorly-mapped read extraction.

Public API
----------
.. autofunction:: extract_reads
.. autofunction:: build_extract_command
.. autofunction:: run_idxstats
"""

from csc.extract.extract import (
    READS_SUMMARY_SCHEMA_VERSION,
    ExtractionResult,
    ReadsSummary,
    build_extract_command,
    extract_reads,
    run_idxstats,
)

__all__ = [
    "ExtractionResult",
    "ReadsSummary",
    "READS_SUMMARY_SCHEMA_VERSION",
    "build_extract_command",
    "extract_reads",
    "run_idxstats",
]
