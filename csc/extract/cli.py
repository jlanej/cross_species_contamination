"""Command-line interface for the CSC extraction module.

Usage examples::

    # Extract unmapped reads from a BAM file
    csc-extract input.bam -o output_dir/

    # Extract unmapped + poorly mapped (MAPQ < 10) reads from a CRAM file
    csc-extract input.cram -o output_dir/ --mapq 10 --reference ref.fa

    # Batch extraction from a list of files
    csc-extract *.bam -o output_dir/ --threads 4

    # Structured JSON logging
    csc-extract input.bam -o output_dir/ --json-log

    # Write a batch summary TSV
    csc-extract *.bam -o output_dir/ --summary summary.tsv
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
from pathlib import Path
from typing import Any

from csc import __version__
from csc.extract.extract import extract_reads
from csc.utils import setup_logging


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="csc-extract",
        description=(
            "Extract unmapped (and optionally poorly mapped) reads from "
            "BAM/CRAM files using samtools.  All output is streamed without "
            "disk intermediates."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  csc-extract sample.bam -o results/\n"
            "  csc-extract sample.cram -o results/ --reference ref.fa\n"
            "  csc-extract *.bam -o results/ --mapq 10 --threads 4\n"
        ),
    )
    parser.add_argument(
        "input",
        nargs="+",
        type=Path,
        help="Input BAM or CRAM file(s).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for FASTQ files.",
    )
    parser.add_argument(
        "--mapq",
        type=int,
        default=None,
        metavar="THRESHOLD",
        help=(
            "Also extract mapped reads with MAPQ below this threshold. "
            "If not set, only truly unmapped reads (flag 4) are extracted."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Additional decompression threads for samtools (default: 1).",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=None,
        help="Reference FASTA file (required for CRAM input).",
    )
    parser.add_argument(
        "--interleaved",
        action="store_true",
        help="Write a single interleaved FASTQ instead of split R1/R2 files.",
    )
    parser.add_argument(
        "--sample-id",
        default=None,
        help="Override sample ID used for output file names.",
    )
    parser.add_argument(
        "--skip-idxstats",
        action="store_true",
        help=(
            "Skip idxstats sidecars and idxstats-derived metrics. "
            "By default idxstats is required and extraction fails if it "
            "cannot be computed."
        ),
    )
    parser.add_argument(
        "--json-log",
        action="store_true",
        help="Emit structured JSON log lines instead of human-readable text.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=None,
        metavar="TSV",
        help="Write a batch summary TSV to this path after all extractions.",
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


def _write_summary_tsv(
    path: Path,
    results: list[dict[str, Any]],
    errors: list[dict[str, str]],
) -> None:
    """Write a batch summary TSV.

    Columns include per-sample idxstats totals (``total_mapped``,
    ``total_unmapped``, ``total_reads``) when available, so cohort-level
    absolute burden can be computed directly from the summary without
    re-reading every ``reads_summary.json``.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "sample_id",
                "input",
                "status",
                "read_count",
                "total_mapped",
                "total_unmapped",
                "total_reads",
                "output_files",
            ]
        )
        for r in results:
            files = ";".join(str(p) for p in r["files"].values())
            writer.writerow([
                r["sample_id"],
                r["input"],
                "OK",
                r["read_count"],
                r.get("total_mapped", ""),
                r.get("total_unmapped", ""),
                r.get("total_reads", ""),
                files,
            ])
        for e in errors:
            writer.writerow([
                Path(e["input"]).stem, e["input"], "FAILED", 0, "", "", "", e["error"],
            ])


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

    results: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []

    for input_path in args.input:
        try:
            result = extract_reads(
                input_path,
                args.output_dir,
                sample_id=args.sample_id,
                mapq_threshold=args.mapq,
                threads=args.threads,
                reference=args.reference,
                interleaved=args.interleaved,
                skip_idxstats=args.skip_idxstats,
            )
            results.append(result)
            for kind, path in sorted(result["files"].items()):
                print(f"  {kind}: {path}")
        except Exception as exc:
            log.error(
                "Failed to process %s: %s", input_path, exc
            )
            errors.append({"input": str(input_path), "error": str(exc)})

    if args.summary is not None:
        _write_summary_tsv(args.summary, results, errors)
        log.info("Batch summary written to %s", args.summary)

    if errors:
        log.error(
            "Failed to process %d file(s): %s",
            len(errors),
            ", ".join(e["input"] for e in errors),
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
