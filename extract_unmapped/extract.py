"""Backward-compatibility shim – use ``csc.extract.extract`` instead."""

from csc.extract.extract import (  # noqa: F401
    SUPPORTED_EXTENSIONS,
    _count_reads,
    _find_samtools,
    _validate_input,
    build_extract_command,
    extract_reads,
)
