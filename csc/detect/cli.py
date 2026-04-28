"""Command-line interface for the CSC outlier-detection module.

Usage examples::

    # Run MAD-based outlier detection on the CPM aggregate matrix
    csc-detect results/taxa_matrix_cpm.tsv -o results/detect/

    # Use IQR method with custom multiplier
    csc-detect results/taxa_matrix_cpm.tsv -o results/detect/ --method iqr --iqr-multiplier 2.0

    # Exclude known kitome taxa
    csc-detect results/taxa_matrix_cpm.tsv -o results/detect/ --kitome-taxa 9606 562

    # Skip population-background subtraction
    csc-detect results/taxa_matrix_cpm.tsv -o results/detect/ --no-subtract-background
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys
from pathlib import Path

from csc import __version__
from csc.aggregate.aggregate import (
    DEFAULT_RANK_FILTER,
    typed_rank_matrix_filename,
)
from csc.detect.detect import detect_outliers
from csc.detect.report import generate_report
from csc.utils import setup_logging

TYPED_BASE_MATRIX_PATTERN = re.compile(r"taxa_matrix_(raw|cpm)(?:_(conf\d+p\d+))?\.tsv")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="csc-detect",
        description=(
            "Statistical outlier detection on a sample-by-taxon matrix.  "
            "Identifies samples with abnormally high taxonomic counts that "
            "may indicate cross-species contamination."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  csc-detect results/taxa_matrix_cpm.tsv -o detect_out/\n"
            "  csc-detect results/taxa_matrix_cpm.tsv -o detect_out/ --method iqr\n"
            "  csc-detect results/taxa_matrix_cpm.tsv -o detect_out/ --kitome-taxa 9606 562\n"
        ),
    )
    parser.add_argument(
        "matrix",
        type=Path,
        help="Path to the sample-by-taxon matrix TSV (from csc-aggregate).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for reports and alert files.",
    )
    parser.add_argument(
        "--method",
        choices=["all", "mad", "iqr", "gmm"],
        default="all",
        help=(
            "Outlier detection method (default: all).  'all' runs MAD, "
            "IQR, and GMM together — recommended because the methods "
            "are complementary."
        ),
    )
    parser.add_argument(
        "--mad-threshold",
        type=float,
        default=3.5,
        help=(
            "Threshold for MAD-based detection: number of scaled-MAD "
            "units above the median (default: 3.5)."
        ),
    )
    parser.add_argument(
        "--iqr-multiplier",
        type=float,
        default=1.5,
        help=(
            "Multiplier for the IQR fence in IQR-based detection "
            "(default: 1.5)."
        ),
    )
    parser.add_argument(
        "--gmm-threshold",
        type=float,
        default=0.5,
        help=(
            "Posterior-probability threshold for the contamination "
            "component in GMM-based detection (default: 0.5)."
        ),
    )
    parser.add_argument(
        "--kitome-taxa",
        nargs="+",
        type=int,
        default=None,
        help=(
            "NCBI taxonomy IDs of known kitome / environmental taxa to "
            "exclude before analysis."
        ),
    )
    parser.add_argument(
        "--no-subtract-background",
        action="store_true",
        help="Skip population-mean subtraction before outlier detection.",
    )
    parser.add_argument(
        "--no-confidence-tiers",
        action="store_true",
        help=(
            "Disable automatic discovery of sibling confidence-tier "
            "matrices (taxa_matrix_<type>_conf<T>.tsv).  By default, "
            "each tier is detected on independently and written to its "
            "own subdirectory (e.g. <output>/conf0p50/)."
        ),
    )
    parser.add_argument(
        "--rank-filter",
        nargs="+",
        default=list(DEFAULT_RANK_FILTER),
        metavar="RANK",
        help=(
            "Taxonomy rank codes to run detection on (default: S G F).  "
            "For each rank, the tool looks for a rank-filtered matrix "
            "(e.g. taxa_matrix_cpm_S.tsv) next to the input matrix.  Results "
            "are written to per-rank subdirectories."
        ),
    )
    parser.add_argument(
        "--json-log",
        action="store_true",
        help="Emit structured JSON log lines instead of human-readable text.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose (DEBUG) logging.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def _rank_matrix_candidates(matrix: Path, rank: str) -> list[Path]:
    """Return candidate rank-matrix paths in highest-to-lowest priority order.

    When the input is a typed base matrix (e.g. ``taxa_matrix_cpm.tsv``
    or ``taxa_matrix_cpm_conf0p50.tsv``), the matching typed rank
    matrix is tried first (e.g. ``taxa_matrix_cpm_S.tsv`` /
    ``taxa_matrix_cpm_S_conf0p50.tsv``).  When the input filename is
    not a recognised typed base matrix, both typed variants are
    returned as candidates.
    """
    matrix_dir = matrix.parent
    matrix_name = matrix.name
    typed_match = TYPED_BASE_MATRIX_PATTERN.fullmatch(matrix_name)
    if typed_match:
        matrix_type = typed_match.group(1)
        tier_suffix = typed_match.group(2) or ""
        return [
            matrix_dir / typed_rank_matrix_filename(rank, matrix_type, tier_suffix),
        ]
    return [
        matrix_dir / typed_rank_matrix_filename(rank, "cpm"),
        matrix_dir / typed_rank_matrix_filename(rank, "raw"),
    ]


CONFIDENCE_TIER_PATTERN = re.compile(r"^conf\d+p\d+$")


def _discover_confidence_tier_matrices(matrix: Path) -> list[Path]:
    """Find sibling confidence-tier matrices for the same matrix type.

    When *matrix* is ``taxa_matrix_cpm.tsv``, returns any sibling
    files matching ``taxa_matrix_cpm_conf*p*.tsv``.  An empty list is
    returned if *matrix* is itself a confidence-tier matrix or its
    name does not match the canonical pattern.
    """
    typed_match = TYPED_BASE_MATRIX_PATTERN.fullmatch(matrix.name)
    if not typed_match:
        return []
    if typed_match.group(2):  # already a tier matrix
        return []
    matrix_type = typed_match.group(1)
    pattern = re.compile(
        rf"^taxa_matrix_{matrix_type}_(conf\d+p\d+)\.tsv$"
    )
    return sorted(
        p for p in matrix.parent.iterdir()
        if p.is_file() and pattern.fullmatch(p.name)
    )


def main(argv: list[str] | None = None) -> int:
    """CLI entry point.  Returns 0 on success, 1 on failure."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    setup_logging(
        level="DEBUG" if args.verbose else "INFO",
        json_format=args.json_log,
    )
    log = logging.getLogger(__name__)

    # --- Fail-early validation ---
    if not args.matrix.exists():
        log.error("Matrix file not found: %s", args.matrix)
        return 1
    if not os.access(args.matrix, os.R_OK):
        log.error("Matrix file not readable: %s", args.matrix)
        return 1

    try:
        _run_detection_for_matrix(args.matrix, args.output_dir, args, log)

        # Discover sibling confidence-tier matrices and run detection for each.
        if not args.no_confidence_tiers:
            tier_matrices = _discover_confidence_tier_matrices(args.matrix)
            for tier_matrix in tier_matrices:
                tier_match = TYPED_BASE_MATRIX_PATTERN.fullmatch(tier_matrix.name)
                tier_suffix = tier_match.group(2) if tier_match else tier_matrix.stem
                tier_out = args.output_dir / tier_suffix
                log.info(
                    "Confidence tier '%s': running detection on %s",
                    tier_suffix, tier_matrix,
                )
                print(f"  confidence tier {tier_suffix}:")
                _run_detection_for_matrix(tier_matrix, tier_out, args, log)

        return 0

    except Exception as exc:
        log.error("Detection failed: %s", exc)
        return 1


def _run_detection_for_matrix(
    matrix: Path,
    output_dir: Path,
    args: argparse.Namespace,
    log: logging.Logger,
) -> None:
    """Run detection on a single matrix (and its per-rank siblings)."""
    result = detect_outliers(
        matrix,
        method=args.method,
        mad_threshold=args.mad_threshold,
        iqr_multiplier=args.iqr_multiplier,
        gmm_threshold=args.gmm_threshold,
        kitome_taxa=args.kitome_taxa,
        subtract_background=not args.no_subtract_background,
    )

    reports = generate_report(result, output_dir)

    summary = result["summary"]
    print(f"  method: {summary['method']}")
    print(f"  samples analysed: {summary['total_samples']}")
    print(f"  taxa analysed: {summary['total_taxa_analysed']}")
    print(f"  flags raised: {summary['flagged_count']}")
    print(f"  flagged samples: {len(summary['flagged_samples'])}")
    print(f"  reports: {output_dir}")
    for name, p in reports.items():
        print(f"    {name}: {p}")

    # Run on per-rank filtered matrices when available
    for rank in args.rank_filter:
        rank_matrix: Path | None = None
        for candidate in _rank_matrix_candidates(matrix, rank):
            if candidate.exists():
                rank_matrix = candidate
                break
        if rank_matrix is None:
            log.info(
                "No rank-%s matrix found for %s, skipping", rank, matrix
            )
            continue

        rank_result = detect_outliers(
            rank_matrix,
            method=args.method,
            mad_threshold=args.mad_threshold,
            iqr_multiplier=args.iqr_multiplier,
            gmm_threshold=args.gmm_threshold,
            kitome_taxa=args.kitome_taxa,
            subtract_background=not args.no_subtract_background,
        )

        rank_out = output_dir / rank
        rank_reports = generate_report(rank_result, rank_out)

        rank_summary = rank_result["summary"]
        print(f"  rank {rank}:")
        print(f"    taxa analysed: {rank_summary['total_taxa_analysed']}")
        print(f"    flags raised: {rank_summary['flagged_count']}")
        for name, p in rank_reports.items():
            print(f"    {name}: {p}")


if __name__ == "__main__":
    sys.exit(main())
