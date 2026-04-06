"""Command-line interface for extract_unmapped.

Usage examples::

    # Extract unmapped reads from a BAM file
    extract-unmapped input.bam -o output_dir/

    # Extract unmapped + poorly mapped (MAPQ < 10) reads from a CRAM file
    extract-unmapped input.cram -o output_dir/ --mapq 10 --reference ref.fa

    # Batch extraction from a list of files
    extract-unmapped *.bam -o output_dir/ --threads 4
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from extract_unmapped import __version__
from extract_unmapped.extract import extract_reads


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="extract-unmapped",
        description=(
            "Extract unmapped (and optionally poorly mapped) reads from "
            "BAM/CRAM files using samtools.  All output is streamed without "
            "disk intermediates."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  extract-unmapped sample.bam -o results/\n"
            "  extract-unmapped sample.cram -o results/ --reference ref.fa\n"
            "  extract-unmapped *.bam -o results/ --mapq 10 --threads 4\n"
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

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    errors: list[str] = []
    for input_path in args.input:
        try:
            outputs = extract_reads(
                input_path,
                args.output_dir,
                sample_id=args.sample_id,
                mapq_threshold=args.mapq,
                threads=args.threads,
                reference=args.reference,
                interleaved=args.interleaved,
            )
            for kind, path in sorted(outputs.items()):
                print(f"  {kind}: {path}")
        except Exception as exc:
            logging.getLogger(__name__).error(
                "Failed to process %s: %s", input_path, exc
            )
            errors.append(str(input_path))

    if errors:
        logging.getLogger(__name__).error(
            "Failed to process %d file(s): %s", len(errors), ", ".join(errors)
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
