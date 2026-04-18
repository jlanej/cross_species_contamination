"""Command-line interface for the CSC aggregation module.

Usage examples::

    # Build sample-by-taxon matrices (raw + CPM) from Kraken2 reports
    csc-aggregate reports/*.kraken2.report.txt -o results/

    # Only include taxa with ≥ 50 direct reads per sample
    csc-aggregate reports/*.kraken2.report.txt -o results/ --min-reads 50

    # Structured JSON logging
    csc-aggregate reports/*.kraken2.report.txt -o results/ --json-log
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

from csc import __version__
from csc.aggregate.aggregate import DEFAULT_RANK_FILTER, aggregate_reports
from csc.utils import setup_logging


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="csc-aggregate",
        description=(
            "Aggregate Kraken2 classification reports into sample-by-taxon "
            "matrices.  Accepts the report files produced by csc-classify and "
            "always writes both raw-count (taxa_matrix_raw.tsv) and CPM "
            "(taxa_matrix_cpm.tsv) TSV matrices."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  csc-aggregate reports/*.kraken2.report.txt -o results/\n"
            "  csc-aggregate reports/*.kraken2.report.txt -o results/ --min-reads 50\n"
        ),
    )
    parser.add_argument(
        "input",
        nargs="+",
        type=Path,
        help="Kraken2 report files (one per sample).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for the taxa matrix and metadata files.",
    )
    parser.add_argument(
        "--min-reads",
        type=int,
        default=0,
        help=(
            "Minimum direct-read count for a taxon to be included per "
            "sample (default: 0).  Note: the default_config.yaml sets "
            "min_reads: 10 which takes effect when using the Nextflow "
            "pipeline; CLI invocations use 0 unless overridden."
        ),
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=500,
        help=(
            "Number of reports to process in each chunk for memory "
            "efficiency (default: 500)."
        ),
    )
    parser.add_argument(
        "--rank-filter",
        nargs="+",
        default=list(DEFAULT_RANK_FILTER),
        metavar="RANK",
        help=(
            "Taxonomy rank codes for which per-rank filtered matrices "
            "are produced (default: S G F).  Common codes: S (species), "
            "G (genus), F (family), D (domain)."
        ),
    )
    parser.add_argument(
        "--db-path",
        type=Path,
        default=None,
        help=(
            "Path to the Kraken2 database directory containing "
            "taxonomy/nodes.dmp.  When provided, a lineage-aware "
            "'domain' column is added to every output matrix."
        ),
    )
    parser.add_argument(
        "--idxstats",
        nargs="+",
        type=Path,
        default=None,
        metavar="JSON",
        help=(
            "Paths to per-sample 'reads_summary.json' sidecars produced "
            "by csc-extract.  When supplied, an absolute-burden matrix "
            "(taxa_matrix_abs.tsv) is written using total_reads as the "
            "denominator.  Samples without a matching sidecar are "
            "written as NA in the absolute matrix."
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

    # --- Fail-early: validate all inputs exist and are readable ---
    missing = [p for p in args.input if not p.exists()]
    if missing:
        for p in missing:
            log.error("Input file not found: %s", p)
        return 1

    unreadable = [p for p in args.input if p.exists() and not os.access(p, os.R_OK)]
    if unreadable:
        for p in unreadable:
            log.error("Input file not readable: %s", p)
        return 1

    try:
        result = aggregate_reports(
            args.input,
            args.output_dir,
            min_reads=args.min_reads,
            chunk_size=args.chunk_size,
            rank_filter=tuple(args.rank_filter),
            db_path=args.db_path,
            idxstats_paths=args.idxstats,
        )
        print(f"  matrix (raw): {result['matrix_raw_path']}")
        print(f"  matrix (cpm): {result['matrix_cpm_path']}")
        if result.get("matrix_abs_path"):
            print(f"  matrix (abs): {result['matrix_abs_path']}")
        print(f"  metadata: {result['metadata_path']}")
        print(
            f"  samples: {result['sample_count']}, "
            f"taxa: {result['taxon_count']}"
        )
        if result["rank_matrices_raw"]:
            for rank, rpath in sorted(result["rank_matrices_raw"].items()):
                print(f"  rank {rank} matrix (raw): {rpath}")
        if result["rank_matrices_cpm"]:
            for rank, rpath in sorted(result["rank_matrices_cpm"].items()):
                print(f"  rank {rank} matrix (cpm): {rpath}")
        for rank, rpath in sorted((result.get("rank_matrices_abs") or {}).items()):
            print(f"  rank {rank} matrix (abs): {rpath}")
        return 0
    except Exception as exc:
        log.error("Aggregation failed: %s", exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
