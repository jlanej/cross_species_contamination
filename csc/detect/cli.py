"""Command-line interface for the CSC outlier-detection module.

Usage examples::

    # Run MAD-based outlier detection on an aggregated matrix
    csc-detect results/taxa_matrix.tsv -o results/detect/

    # Use IQR method with custom multiplier
    csc-detect results/taxa_matrix.tsv -o results/detect/ --method iqr --iqr-multiplier 2.0

    # Exclude known kitome taxa
    csc-detect results/taxa_matrix.tsv -o results/detect/ --kitome-taxa 9606 562

    # Skip population-background subtraction
    csc-detect results/taxa_matrix.tsv -o results/detect/ --no-subtract-background
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

from csc import __version__
from csc.aggregate.aggregate import DEFAULT_RANK_FILTER
from csc.detect.detect import detect_outliers
from csc.detect.report import generate_report
from csc.utils import setup_logging


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
            "  csc-detect results/taxa_matrix.tsv -o detect_out/\n"
            "  csc-detect results/taxa_matrix.tsv -o detect_out/ --method iqr\n"
            "  csc-detect results/taxa_matrix.tsv -o detect_out/ --kitome-taxa 9606 562\n"
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
        choices=["mad", "iqr"],
        default="mad",
        help="Outlier detection method (default: mad).",
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
        "--rank-filter",
        nargs="+",
        default=list(DEFAULT_RANK_FILTER),
        metavar="RANK",
        help=(
            "Taxonomy rank codes to run detection on (default: S G F).  "
            "For each rank, the tool looks for a rank-filtered matrix "
            "(e.g. taxa_matrix_S.tsv) next to the input matrix.  Results "
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
        # Always run on the primary (unfiltered) matrix
        result = detect_outliers(
            args.matrix,
            method=args.method,
            mad_threshold=args.mad_threshold,
            iqr_multiplier=args.iqr_multiplier,
            kitome_taxa=args.kitome_taxa,
            subtract_background=not args.no_subtract_background,
        )

        reports = generate_report(result, args.output_dir)

        summary = result["summary"]
        print(f"  method: {summary['method']}")
        print(f"  samples analysed: {summary['total_samples']}")
        print(f"  taxa analysed: {summary['total_taxa_analysed']}")
        print(f"  flags raised: {summary['flagged_count']}")
        print(f"  flagged samples: {len(summary['flagged_samples'])}")
        print(f"  reports: {args.output_dir}")
        for name, p in reports.items():
            print(f"    {name}: {p}")

        # Run on per-rank filtered matrices when available
        matrix_dir = args.matrix.parent
        for rank in args.rank_filter:
            rank_matrix = matrix_dir / f"taxa_matrix_{rank}.tsv"
            if not rank_matrix.exists():
                log.info(
                    "No rank-%s matrix found at %s, skipping", rank, rank_matrix
                )
                continue

            rank_result = detect_outliers(
                rank_matrix,
                method=args.method,
                mad_threshold=args.mad_threshold,
                iqr_multiplier=args.iqr_multiplier,
                kitome_taxa=args.kitome_taxa,
                subtract_background=not args.no_subtract_background,
            )

            rank_out = args.output_dir / rank
            rank_reports = generate_report(rank_result, rank_out)

            rank_summary = rank_result["summary"]
            print(f"  rank {rank}:")
            print(f"    taxa analysed: {rank_summary['total_taxa_analysed']}")
            print(f"    flags raised: {rank_summary['flagged_count']}")
            for name, p in rank_reports.items():
                print(f"    {name}: {p}")

        return 0

    except Exception as exc:
        log.error("Detection failed: %s", exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
