"""Detect module – statistical detection of cross-species contamination.

This module applies statistical methods (MAD, IQR, GMM) to aggregated
classification results to flag samples with significant cross-species
contamination.  It also supports kitome/environmental-background
exclusion and generates reviewable reports and quarantine lists.

Public API
----------
detect_outliers
    Run outlier detection on a sample-by-taxon matrix.
generate_report
    Write flagged-sample table, QC summary, and quarantine list.
load_matrix
    Read a taxa-matrix TSV produced by :mod:`csc.aggregate`.
filter_kitome
    Remove known kitome/environmental taxa before analysis.
subtract_population_background
    Subtract per-taxon population mean from every sample.
DetectionResult
    TypedDict describing the return value of :func:`detect_outliers`.
FlaggedSample
    TypedDict for a single flagged outlier observation.
"""

from csc.detect.detect import (
    DetectionResult,
    FlaggedSample,
    detect_outliers,
    filter_kitome,
    load_matrix,
    subtract_population_background,
)
from csc.detect.report import generate_report

__all__ = [
    "DetectionResult",
    "FlaggedSample",
    "detect_outliers",
    "filter_kitome",
    "generate_report",
    "load_matrix",
    "subtract_population_background",
]
