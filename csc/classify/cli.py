"""Command-line interface for the CSC classification module.

Usage examples::

    # Classify extracted FASTQ reads with Kraken2
    csc-classify reads.fastq.gz --db /data/kraken2/PlusPF -o output_dir/

    # Paired-end classification
    csc-classify R1.fastq.gz R2.fastq.gz --db /data/kraken2/PlusPF -o output_dir/ --paired

    # Use memory mapping to reduce RAM usage
    csc-classify reads.fastq.gz --db /data/kraken2/PlusPF -o output_dir/ --memory-mapping

    # Set confidence threshold
    csc-classify reads.fastq.gz --db /data/kraken2/PlusPF -o output_dir/ --confidence 0.2

    # Structured JSON logging
    csc-classify reads.fastq.gz --db /data/kraken2/PlusPF -o output_dir/ --json-log
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Any

from csc import __version__
from csc.classify.classify import classify_reads
from csc.utils import setup_logging


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="csc-classify",
        description=(
            "Classify FASTQ reads taxonomically using Kraken2.  Accepts "
            "output from the csc-extract module and produces standardized "
            "Kraken2 report and per-read classification files."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  csc-classify reads.fastq.gz --db /data/kraken2db -o results/\n"
            "  csc-classify R1.fq.gz R2.fq.gz --db /data/kraken2db -o results/ --paired\n"
            "  csc-classify reads.fastq.gz --db /data/kraken2db -o results/ --memory-mapping\n"
        ),
    )
    parser.add_argument(
        "input",
        nargs="+",
        type=Path,
        help="Input FASTQ file(s).  For paired-end, supply R1 and R2.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for Kraken2 report and classification files.",
    )
    parser.add_argument(
        "--db",
        type=Path,
        required=True,
        help="Path to the Kraken2 database directory.",
    )
    parser.add_argument(
        "--confidence",
        type=float,
        default=0.0,
        help="Minimum confidence score for Kraken2 classification (default: 0.0).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for Kraken2 (default: 1).",
    )
    parser.add_argument(
        "--memory-mapping",
        action="store_true",
        help=(
            "Use memory mapping for the Kraken2 database instead of loading "
            "it into RAM.  Reduces memory usage at the cost of speed."
        ),
    )
    parser.add_argument(
        "--paired",
        action="store_true",
        help="Treat input files as paired-end reads (requires exactly 2 inputs).",
    )
    parser.add_argument(
        "--sample-id",
        default=None,
        help="Override sample ID used for output file names.",
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
        result = classify_reads(
            args.input,
            args.output_dir,
            db=args.db,
            sample_id=args.sample_id,
            confidence=args.confidence,
            threads=args.threads,
            memory_mapping=args.memory_mapping,
            paired=args.paired,
        )
        print(f"  report: {result['report']}")
        print(f"  output: {result['output']}")
        return 0
    except Exception as exc:
        log.error("Classification failed: %s", exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
