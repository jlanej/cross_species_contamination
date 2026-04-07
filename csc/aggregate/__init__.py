"""Aggregate module – aggregation of Kraken2 classification results.

This module collects per-sample Kraken2 report files and produces a
sample-by-taxon count matrix with optional CPM normalisation.  It is
designed to handle very large cohorts (100K+ samples) efficiently.

Public API
----------
parse_kraken2_report
    Parse a single Kraken2 report into structured records.
aggregate_reports
    Build a sample×taxon matrix from a list of report paths.
sample_id_from_report
    Derive a sample identifier from a report file path.
AggregationResult
    TypedDict describing the return value of :func:`aggregate_reports`.
TaxonRecord
    TypedDict for a single parsed taxon entry.
"""

from csc.aggregate.aggregate import (
    AggregationResult,
    TaxonRecord,
    aggregate_reports,
    parse_kraken2_report,
    sample_id_from_report,
)

__all__ = [
    "AggregationResult",
    "TaxonRecord",
    "aggregate_reports",
    "parse_kraken2_report",
    "sample_id_from_report",
]
