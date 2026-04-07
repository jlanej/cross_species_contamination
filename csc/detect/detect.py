"""Core outlier-detection logic for the CSC pipeline.

Applies statistical methods to a sample-by-taxon count matrix
(as produced by :mod:`csc.aggregate`) to identify samples whose
taxonomic profile deviates significantly from the population,
suggesting cross-species contamination.

Two complementary strategies are implemented:

* **MAD (Median Absolute Deviation)** – robust measure of dispersion
  that is insensitive to extreme outliers.
* **IQR (Inter-Quartile Range)** – flags values beyond
  ``Q1 - k·IQR`` or ``Q3 + k·IQR``.

An optional *kitome* / environmental-background exclusion list removes
taxa known to originate from laboratory reagents or the environment so
that they do not trigger false positives.

AI assistance acknowledgment: This module was developed with AI assistance.
Best practices in the bioinformatics field should always take precedence over
specific implementation details.
"""

from __future__ import annotations

import csv
import logging
import math
from pathlib import Path
from typing import Any, TypedDict

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public types
# ---------------------------------------------------------------------------


class FlaggedSample(TypedDict):
    """A single sample flagged as an outlier for a specific taxon."""

    sample_id: str
    tax_id: int
    taxon_name: str
    value: float
    population_median: float
    deviation: float
    method: str


class DetectionResult(TypedDict):
    """Return type for :func:`detect_outliers`."""

    flagged: list[FlaggedSample]
    summary: dict[str, Any]


# ---------------------------------------------------------------------------
# Matrix I/O
# ---------------------------------------------------------------------------


def load_matrix(path: str | Path) -> tuple[list[str], list[dict[str, Any]], dict[int, str]]:
    """Read a taxa-matrix TSV produced by :mod:`csc.aggregate`.

    Returns
    -------
    sample_ids : list[str]
        Column headers (sample identifiers).
    rows : list[dict]
        One dict per taxon row with keys ``tax_id`` (int) and one float
        entry per sample.
    tax_names : dict[int, str]
        Mapping from tax_id to taxon name.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Matrix file not found: {path}")

    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)  # tax_id, name, sample1, sample2, …

    if len(header) < 3:
        raise ValueError(
            f"Matrix must have at least 3 columns (tax_id, name, and one "
            f"sample), got {len(header)}"
        )

    sample_ids = header[2:]
    rows: list[dict[str, Any]] = []
    tax_names: dict[int, str] = {}

    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)  # skip header
        for lineno, cols in enumerate(reader, 2):
            if len(cols) < 3:
                logger.warning("Skipping line %d: too few columns", lineno)
                continue
            try:
                tid = int(cols[0])
            except ValueError:
                logger.warning("Skipping line %d: non-integer tax_id", lineno)
                continue
            name = cols[1]
            tax_names[tid] = name
            entry: dict[str, Any] = {"tax_id": tid}
            for i, sid in enumerate(sample_ids):
                try:
                    entry[sid] = float(cols[2 + i])
                except (ValueError, IndexError):
                    entry[sid] = 0.0
            rows.append(entry)

    if not rows:
        raise ValueError(f"No valid taxon rows in matrix: {path}")

    return sample_ids, rows, tax_names


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------


def _median(values: list[float]) -> float:
    """Return the median of *values*."""
    n = len(values)
    if n == 0:
        return 0.0
    s = sorted(values)
    mid = n // 2
    if n % 2 == 1:
        return s[mid]
    return (s[mid - 1] + s[mid]) / 2.0


def _mad(values: list[float]) -> float:
    """Return the Median Absolute Deviation of *values*."""
    med = _median(values)
    deviations = [abs(v - med) for v in values]
    return _median(deviations)


def _quartiles(values: list[float]) -> tuple[float, float, float]:
    """Return (Q1, median, Q3) using the exclusive method."""
    s = sorted(values)
    n = len(s)
    if n == 0:
        return (0.0, 0.0, 0.0)
    med = _median(s)
    mid = n // 2
    lower = s[:mid]
    upper = s[mid:] if n % 2 == 0 else s[mid + 1:]
    if not lower:
        lower = s[:1]
    if not upper:
        upper = s[-1:]
    q1 = _median(lower)
    q3 = _median(upper)
    return q1, med, q3


# ---------------------------------------------------------------------------
# Background / kitome filtering
# ---------------------------------------------------------------------------


def filter_kitome(
    rows: list[dict[str, Any]],
    kitome_taxa: list[int] | None,
) -> list[dict[str, Any]]:
    """Remove rows whose ``tax_id`` appears in the kitome list.

    Parameters
    ----------
    rows:
        Taxon rows from :func:`load_matrix`.
    kitome_taxa:
        NCBI taxonomy IDs known to be laboratory / environmental
        contaminants.  If *None* or empty, no filtering is performed.

    Returns
    -------
    list[dict]
        Filtered rows.
    """
    if not kitome_taxa:
        return rows

    kitome_set = set(kitome_taxa)
    before = len(rows)
    filtered = [r for r in rows if r["tax_id"] not in kitome_set]
    removed = before - len(filtered)
    if removed:
        logger.info(
            "Kitome filter removed %d / %d taxa", removed, before
        )
    return filtered


# ---------------------------------------------------------------------------
# Population-mean subtraction
# ---------------------------------------------------------------------------


def subtract_population_background(
    rows: list[dict[str, Any]],
    sample_ids: list[str],
) -> list[dict[str, Any]]:
    """Subtract the per-taxon population mean from every sample.

    After subtraction a sample sitting at the population mean will have
    a value of zero.  Negative values are preserved so that the
    distribution shape (and therefore MAD / IQR) remains meaningful.

    Parameters
    ----------
    rows:
        Taxon rows from :func:`load_matrix`.
    sample_ids:
        Sample identifiers (column names).

    Returns
    -------
    list[dict]
        Adjusted rows (new list; originals are not mutated).
    """
    adjusted: list[dict[str, Any]] = []
    for row in rows:
        values = [row.get(sid, 0.0) for sid in sample_ids]
        mean_val = sum(values) / len(values) if values else 0.0
        new_row: dict[str, Any] = {"tax_id": row["tax_id"]}
        for sid in sample_ids:
            new_row[sid] = row.get(sid, 0.0) - mean_val
        adjusted.append(new_row)
    return adjusted


# ---------------------------------------------------------------------------
# Outlier detection
# ---------------------------------------------------------------------------

# Scaling constant so that MAD ≈ standard deviation for normal data.
_MAD_SCALE = 1.4826


def detect_outliers(
    matrix_path: str | Path,
    *,
    method: str = "mad",
    mad_threshold: float = 3.5,
    iqr_multiplier: float = 1.5,
    kitome_taxa: list[int] | None = None,
    subtract_background: bool = True,
) -> DetectionResult:
    """Run outlier detection on a sample-by-taxon matrix.

    Parameters
    ----------
    matrix_path:
        Path to the TSV matrix produced by ``csc-aggregate``.
    method:
        ``"mad"`` for MAD-based detection, ``"iqr"`` for IQR-based
        detection.
    mad_threshold:
        Number of scaled-MAD units above the median to consider a
        value an outlier (used when *method* is ``"mad"``).
    iqr_multiplier:
        Multiplier for the IQR fence (used when *method* is
        ``"iqr"``).
    kitome_taxa:
        NCBI taxonomy IDs to exclude before analysis.
    subtract_background:
        If ``True``, subtract the population mean per taxon before
        testing.

    Returns
    -------
    DetectionResult
        A dict with ``flagged`` (list of :class:`FlaggedSample`) and
        ``summary`` (aggregate statistics).

    Raises
    ------
    ValueError
        If *method* is not one of ``"mad"`` or ``"iqr"``.
    """
    method = method.lower()
    if method not in ("mad", "iqr"):
        raise ValueError(f"Unknown detection method: {method!r}")

    sample_ids, rows, tax_names = load_matrix(matrix_path)
    rows = filter_kitome(rows, kitome_taxa)

    if subtract_background:
        rows = subtract_population_background(rows, sample_ids)

    flagged: list[FlaggedSample] = []

    for row in rows:
        tid = row["tax_id"]
        values = [row.get(sid, 0.0) for sid in sample_ids]

        if method == "mad":
            med = _median(values)
            mad_val = _mad(values)
            if mad_val == 0.0:
                # MAD=0 means more than half the values are identical.
                # No robust dispersion is detectable, so skip this taxon.
                continue

            scaled_mad = mad_val * _MAD_SCALE
            for sid, v in zip(sample_ids, values):
                modified_z = abs(v - med) / scaled_mad
                if modified_z > mad_threshold and v > med:
                    flagged.append(
                        FlaggedSample(
                            sample_id=sid,
                            tax_id=tid,
                            taxon_name=tax_names.get(tid, ""),
                            value=v,
                            population_median=med,
                            deviation=modified_z,
                            method="mad",
                        )
                    )

        else:  # iqr
            q1, med, q3 = _quartiles(values)
            iqr = q3 - q1
            upper_fence = q3 + iqr_multiplier * iqr

            for sid, v in zip(sample_ids, values):
                if v > upper_fence and v > 0:
                    flagged.append(
                        FlaggedSample(
                            sample_id=sid,
                            tax_id=tid,
                            taxon_name=tax_names.get(tid, ""),
                            value=v,
                            population_median=med,
                            deviation=v - upper_fence,
                            method="iqr",
                        )
                    )

    summary: dict[str, Any] = {
        "total_samples": len(sample_ids),
        "total_taxa_analysed": len(rows),
        "method": method,
        "flagged_count": len(flagged),
        "flagged_samples": sorted({f["sample_id"] for f in flagged}),
        "kitome_taxa_excluded": len(kitome_taxa) if kitome_taxa else 0,
        "background_subtracted": subtract_background,
    }

    logger.info(
        "Detection complete: %d flags across %d samples",
        len(flagged),
        len(summary["flagged_samples"]),
    )

    return DetectionResult(flagged=flagged, summary=summary)
