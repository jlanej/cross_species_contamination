"""Extraction module – streaming unmapped/poorly-mapped read extraction.

Public API
----------
.. autofunction:: extract_reads
.. autofunction:: build_extract_command
"""

from csc.extract.extract import build_extract_command, extract_reads

__all__ = ["build_extract_command", "extract_reads"]
