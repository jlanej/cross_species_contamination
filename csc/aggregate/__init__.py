"""Aggregate module – aggregation of Kraken2 classification results.

This module collects per-sample Kraken2 report files and produces a
sample-by-taxon count matrix in both raw counts and CPM-normalised form.
It is designed to handle very large cohorts (100K+ samples) efficiently.

Public API
----------
parse_kraken2_report
    Parse a single Kraken2 report into structured records.
aggregate_reports
    Build sample×taxon matrices (raw + CPM) from a list of report paths.
sample_id_from_report
    Derive a sample identifier from a report file path.
AggregationResult
    TypedDict describing the return value of :func:`aggregate_reports`.
TaxonRecord
    TypedDict for a single parsed taxon entry.
"""

from csc.aggregate.aggregate import (
    AggregationResult,
    DEFAULT_RANK_FILTER,
    TaxonRecord,
    VALID_RANK_CODES,
    aggregate_reports,
    parse_kraken2_report,
    sample_id_from_report,
    typed_matrix_filename,
    typed_rank_matrix_filename,
)
from csc.aggregate.taxonomy import (
    assign_domains,
    load_taxonomy_tree,
)

__all__ = [
    "AggregationResult",
    "DEFAULT_RANK_FILTER",
    "TaxonRecord",
    "VALID_RANK_CODES",
    "aggregate_reports",
    "assign_domains",
    "load_taxonomy_tree",
    "parse_kraken2_report",
    "sample_id_from_report",
    "typed_matrix_filename",
    "typed_rank_matrix_filename",
]
