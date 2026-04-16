"""Core aggregation logic for Kraken2 classification reports.

This module provides functions to parse Kraken2 report files, build
sample-by-taxon count matrices, and normalize counts (e.g. CPM).
It is designed to handle large cohorts (100K+ samples) efficiently
by processing reports in configurable chunks and using a streaming
dictionary-of-counters approach that avoids loading all data into
a single in-memory structure at once.

AI assistance acknowledgment: This module was developed with AI assistance.
Best practices in the bioinformatics field should always take precedence over
specific implementation details.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any, TypedDict

logger = logging.getLogger(__name__)

# Kraken2 report columns (tab-separated):
#   0: percentage of reads rooted at this taxon
#   1: number of reads rooted at this taxon (clade)
#   2: number of reads directly assigned to this taxon
#   3: rank code (U, R, D, P, C, O, F, G, S, …)
#   4: NCBI taxonomy ID
#   5: scientific name (may have leading spaces)
_COL_PCT = 0
_COL_CLADE_READS = 1
_COL_DIRECT_READS = 2
_COL_RANK = 3
_COL_TAX_ID = 4
_COL_NAME = 5

# Number of expected columns in a well-formed Kraken2 report line.
_MIN_COLUMNS = 6


class TaxonRecord(TypedDict):
    """A single taxon entry parsed from a Kraken2 report."""

    tax_id: int
    name: str
    rank: str
    clade_reads: int
    direct_reads: int
    percentage: float


class AggregationResult(TypedDict):
    """Return type for :func:`aggregate_reports`."""

    matrix_path: Path
    matrix_raw_path: Path
    matrix_cpm_path: Path
    metadata_path: Path
    sample_count: int
    taxon_count: int
    rank_matrices: dict[str, Path]
    rank_matrices_raw: dict[str, Path]
    rank_matrices_cpm: dict[str, Path]
    rank_metadata_path: Path


# Valid Kraken2 rank codes.
VALID_RANK_CODES = ("U", "R", "D", "P", "C", "O", "F", "G", "S")

# Default ranks for which per-rank filtered matrices are produced.
DEFAULT_RANK_FILTER: tuple[str, ...] = ("S", "G", "F")


def rank_matrix_filename(rank: str) -> str:
    """Return the canonical filename for a rank-filtered matrix."""
    return f"taxa_matrix_{rank}.tsv"


def typed_matrix_filename(matrix_type: str) -> str:
    """Return canonical filename for a typed unfiltered matrix."""
    return f"taxa_matrix_{matrix_type}.tsv"


def typed_rank_matrix_filename(rank: str, matrix_type: str) -> str:
    """Return canonical filename for a typed rank-filtered matrix."""
    return f"taxa_matrix_{matrix_type}_{rank}.tsv"


# ---- Parsing -----------------------------------------------------------------


def parse_kraken2_report(path: str | Path) -> list[TaxonRecord]:
    """Parse a single Kraken2 report file.

    Parameters
    ----------
    path:
        Path to a Kraken2-format report (tab-separated).

    Returns
    -------
    list[TaxonRecord]
        One entry per taxon line in the report.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If the file is empty or contains no parseable taxon lines.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Report file not found: {path}")

    records: list[TaxonRecord] = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < _MIN_COLUMNS:
                logger.warning(
                    "%s:%d: skipping malformed line (%d columns, need %d)",
                    path,
                    lineno,
                    len(cols),
                    _MIN_COLUMNS,
                )
                continue
            try:
                records.append(
                    TaxonRecord(
                        tax_id=int(cols[_COL_TAX_ID].strip()),
                        name=cols[_COL_NAME].strip(),
                        rank=cols[_COL_RANK].strip(),
                        clade_reads=int(cols[_COL_CLADE_READS].strip()),
                        direct_reads=int(cols[_COL_DIRECT_READS].strip()),
                        percentage=float(cols[_COL_PCT].strip()),
                    )
                )
            except (ValueError, IndexError) as exc:
                logger.warning(
                    "%s:%d: skipping unparseable line: %s", path, lineno, exc
                )
    if not records:
        raise ValueError(f"No valid taxon lines found in {path}")
    return records


def sample_id_from_report(path: str | Path) -> str:
    """Derive a sample ID from a Kraken2 report file path.

    Strips the ``.kraken2.report.txt`` suffix (if present) and returns
    the remaining file stem.
    """
    name = Path(path).name
    for suffix in (".kraken2.report.txt",):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(path).stem


# ---- Matrix construction -----------------------------------------------------


def _collect_sample_counts(
    records: list[TaxonRecord],
    min_reads: int = 0,
) -> dict[int, int]:
    """Return {tax_id: direct_reads} for records passing *min_reads*."""
    return {
        r["tax_id"]: r["direct_reads"]
        for r in records
        if r["direct_reads"] >= min_reads
    }


def aggregate_reports(
    report_paths: list[str | Path],
    output_dir: str | Path,
    *,
    min_reads: int = 0,
    normalize: bool = True,
    chunk_size: int = 500,
    rank_filter: tuple[str, ...] | list[str] = DEFAULT_RANK_FILTER,
) -> AggregationResult:
    """Build a sample-by-taxon matrix from Kraken2 reports.

    The matrix rows are taxa (identified by NCBI taxonomy ID) and the
    columns are samples.  Values are either raw direct-read counts or,
    when *normalize* is ``True``, counts-per-million (CPM) based on
    each sample's total classified reads.

    Processing is chunked so that memory stays bounded even for very
    large cohorts.

    In addition to the unfiltered matrix, per-rank filtered matrices are
    written for each rank code in *rank_filter* (e.g. ``taxa_matrix_S.tsv``
    for species).  A sidecar ``rank_filter_metadata.json`` records which
    taxa were retained in each rank.

    Parameters
    ----------
    report_paths:
        Paths to Kraken2 report files (one per sample).
    output_dir:
        Directory for output files.  Created if it does not exist.
    min_reads:
        Minimum direct-read count for a taxon to be included per sample.
    normalize:
        If ``True``, values in the matrix are reads-per-million (CPM).
    chunk_size:
        Number of reports to process before flushing intermediate state.
        Helps keep memory bounded for very large cohorts.
    rank_filter:
        Taxonomy rank codes for which per-rank matrices are produced.
        Defaults to ``("S", "G", "F")`` (species, genus, family).

    Returns
    -------
    AggregationResult
        Paths to the output matrix TSV and metadata JSON, plus summary
        counts.

    Raises
    ------
    ValueError
        If *report_paths* is empty.
    """
    if not report_paths:
        raise ValueError("At least one report path is required.")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rank_filter = tuple(rank_filter)
    invalid = [r for r in rank_filter if r not in VALID_RANK_CODES]
    if invalid:
        raise ValueError(
            f"Invalid rank code(s): {invalid}. "
            f"Valid codes are: {VALID_RANK_CODES}"
        )

    # Accumulate per-sample data: {sample_id: {tax_id: count}}
    sample_data: dict[str, dict[int, int]] = {}
    # Map tax_id -> scientific name (last one wins; should be consistent)
    tax_names: dict[int, str] = {}
    # Map tax_id -> rank code (last one wins; should be consistent)
    tax_ranks: dict[int, str] = {}
    # Track errors
    errors: list[dict[str, str]] = []

    resolved = [Path(p) for p in report_paths]

    for i, rpath in enumerate(resolved):
        sid = sample_id_from_report(rpath)
        try:
            records = parse_kraken2_report(rpath)
        except (FileNotFoundError, ValueError) as exc:
            logger.warning("Skipping %s: %s", rpath, exc)
            errors.append({"file": str(rpath), "error": str(exc)})
            continue

        counts = _collect_sample_counts(records, min_reads=min_reads)
        sample_data[sid] = counts

        for rec in records:
            if rec["direct_reads"] >= min_reads:
                tax_names[rec["tax_id"]] = rec["name"]
                tax_ranks[rec["tax_id"]] = rec["rank"]

        if (i + 1) % chunk_size == 0:
            logger.info("Processed %d / %d reports", i + 1, len(resolved))

    logger.info(
        "Processed %d reports (%d skipped)",
        len(sample_data),
        len(errors),
    )

    # Build the sorted union of all taxa across samples
    all_taxa = sorted(tax_names.keys())
    sample_ids = sorted(sample_data.keys())

    # Compute per-sample totals (for CPM normalisation)
    sample_totals: dict[str, int] = {}
    for sid in sample_ids:
        sample_totals[sid] = sum(sample_data[sid].values())

    # Always write both unfiltered matrices (raw + CPM)
    matrix_raw_path = output_dir / typed_matrix_filename("raw")
    matrix_cpm_path = output_dir / typed_matrix_filename("cpm")
    _write_matrix(
        matrix_raw_path,
        sample_ids=sample_ids,
        all_taxa=all_taxa,
        tax_names=tax_names,
        sample_data=sample_data,
        sample_totals=sample_totals,
        normalize=False,
    )
    _write_matrix(
        matrix_cpm_path,
        sample_ids=sample_ids,
        all_taxa=all_taxa,
        tax_names=tax_names,
        sample_data=sample_data,
        sample_totals=sample_totals,
        normalize=True,
    )

    # Backward-compatible primary matrix name
    matrix_path = output_dir / "taxa_matrix.tsv"
    _write_matrix(
        matrix_path,
        sample_ids=sample_ids,
        all_taxa=all_taxa,
        tax_names=tax_names,
        sample_data=sample_data,
        sample_totals=sample_totals,
        normalize=normalize,
    )

    # Write per-rank filtered matrices (always raw + CPM), plus primary legacy names
    rank_matrices: dict[str, Path] = {}
    rank_matrices_raw: dict[str, Path] = {}
    rank_matrices_cpm: dict[str, Path] = {}
    rank_sidecar: dict[str, Any] = {}
    for rank in rank_filter:
        rank_taxa = [t for t in all_taxa if tax_ranks.get(t) == rank]
        if not rank_taxa:
            logger.info("Rank '%s': no taxa found, skipping matrix", rank)
            continue

        rank_raw_path = output_dir / typed_rank_matrix_filename(rank, "raw")
        _write_matrix(
            rank_raw_path,
            sample_ids=sample_ids,
            all_taxa=rank_taxa,
            tax_names=tax_names,
            sample_data=sample_data,
            sample_totals=sample_totals,
            normalize=False,
        )
        rank_matrices_raw[rank] = rank_raw_path

        rank_cpm_path = output_dir / typed_rank_matrix_filename(rank, "cpm")
        _write_matrix(
            rank_cpm_path,
            sample_ids=sample_ids,
            all_taxa=rank_taxa,
            tax_names=tax_names,
            sample_data=sample_data,
            sample_totals=sample_totals,
            normalize=True,
        )
        rank_matrices_cpm[rank] = rank_cpm_path

        rank_path = output_dir / rank_matrix_filename(rank)
        _write_matrix(
            rank_path,
            sample_ids=sample_ids,
            all_taxa=rank_taxa,
            tax_names=tax_names,
            sample_data=sample_data,
            sample_totals=sample_totals,
            normalize=normalize,
        )
        rank_matrices[rank] = rank_path
        rank_sidecar[rank] = {
            "matrix_path": str(rank_path),
            "matrix_raw_path": str(rank_raw_path),
            "matrix_cpm_path": str(rank_cpm_path),
            "taxon_count": len(rank_taxa),
            "taxa": [
                {"tax_id": t, "name": tax_names.get(t, "")}
                for t in rank_taxa
            ],
        }
        logger.info(
            "Rank '%s': wrote %d taxa to %s", rank, len(rank_taxa), rank_path
        )

    # Write rank-filter metadata sidecar
    rank_metadata_path = output_dir / "rank_filter_metadata.json"
    rank_meta_doc: dict[str, Any] = {
        "rank_filter": list(rank_filter),
        "ranks": rank_sidecar,
    }
    with open(rank_metadata_path, "w") as fh:
        json.dump(rank_meta_doc, fh, indent=2)

    # Write metadata
    metadata_path = output_dir / "aggregation_metadata.json"
    meta: dict[str, Any] = {
        "sample_count": len(sample_ids),
        "taxon_count": len(all_taxa),
        "min_reads": min_reads,
        "normalized": normalize,
        "normalization_method": "CPM" if normalize else "raw",
        "matrix_paths": {
            "primary": matrix_path.name,
            "raw": matrix_raw_path.name,
            "cpm": matrix_cpm_path.name,
        },
        "rank_filter": list(rank_filter),
        "samples": sample_ids,
        "errors": errors,
    }
    with open(metadata_path, "w") as fh:
        json.dump(meta, fh, indent=2)

    return AggregationResult(
        matrix_path=matrix_path,
        matrix_raw_path=matrix_raw_path,
        matrix_cpm_path=matrix_cpm_path,
        metadata_path=metadata_path,
        sample_count=len(sample_ids),
        taxon_count=len(all_taxa),
        rank_matrices=rank_matrices,
        rank_matrices_raw=rank_matrices_raw,
        rank_matrices_cpm=rank_matrices_cpm,
        rank_metadata_path=rank_metadata_path,
    )


def _write_matrix(
    path: Path,
    *,
    sample_ids: list[str],
    all_taxa: list[int],
    tax_names: dict[int, str],
    sample_data: dict[str, dict[int, int]],
    sample_totals: dict[str, int],
    normalize: bool,
) -> None:
    """Write the taxa-by-sample matrix as a TSV file.

    Rows are taxa; columns are samples.  First two columns are
    ``tax_id`` and ``name``.
    """
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        # Header
        writer.writerow(["tax_id", "name"] + sample_ids)
        # Data rows
        for tid in all_taxa:
            row: list[str] = [str(tid), tax_names.get(tid, "")]
            for sid in sample_ids:
                raw = sample_data[sid].get(tid, 0)
                if normalize:
                    total = sample_totals[sid]
                    value = (raw / total * 1_000_000) if total > 0 else 0.0
                    row.append(f"{value:.4f}")
                else:
                    row.append(str(raw))
            writer.writerow(row)

    logger.info("Wrote matrix to %s", path)
