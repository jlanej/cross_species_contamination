"""Command-line interface for the CSC HTML report module.

Usage examples::

    # Generate the report from an aggregation output directory
    csc-report aggregate_out/ -o report/contamination_report.html

    # Include detect-module outputs so §2.4 and the flagged-sample
    # counts appear
    csc-report aggregate_out/ -o report/contamination_report.html \\
        --detect-dir detect_out/

    # Use a stricter variant-calling impact threshold (0.01 % = 100 ppm)
    csc-report aggregate_out/ -o report/contamination_report.html \\
        --variant-impact-threshold-ppm 100
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from csc import __version__
from csc.report.report import (
    DEFAULT_VARIANT_IMPACT_THRESHOLD_PPM,
    generate_html_report,
    load_inputs,
)
from csc.utils import setup_logging


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="csc-report",
        description=(
            "Generate a static, self-contained HTML summary report of "
            "non-human content in a cohort of WGS samples.  Consumes "
            "outputs of csc-aggregate (and optionally csc-detect) and "
            "produces a single HTML file plus a report_manifest.json "
            "sidecar.  The report labels CPM (compositional) and "
            "absolute burden (per million total sequenced reads) "
            "explicitly in every table and figure caption, and includes "
            "a Variant-Calling Impact section that uses the absolute "
            "denominator to flag samples of concern."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "aggregate_dir",
        type=Path,
        help=(
            "Directory containing csc-aggregate outputs "
            "(taxa_matrix_raw.tsv, taxa_matrix_cpm.tsv, optional "
            "taxa_matrix_abs.tsv, and aggregation_metadata.json)."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the output HTML file.",
    )
    parser.add_argument(
        "--detect-dir",
        type=Path,
        default=None,
        help=(
            "Optional directory containing csc-detect outputs "
            "(flagged_samples.tsv, qc_summary.json).  When supplied, "
            "outlier-detection parameters and flagged samples are "
            "included in the report."
        ),
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Number of taxa to show in the cohort-wide top table (default: 10).",
    )
    parser.add_argument(
        "--variant-impact-threshold-ppm",
        type=float,
        default=DEFAULT_VARIANT_IMPACT_THRESHOLD_PPM,
        help=(
            "Absolute-burden threshold (reads per million total "
            f"sequenced reads) above which samples are flagged in §4 "
            f"Variant-Calling Impact (default: "
            f"{DEFAULT_VARIANT_IMPACT_THRESHOLD_PPM:.0f} ppm, = 0.1%% "
            "of total sequencing).  Document the chosen value in your "
            "manuscript methods."
        ),
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Cross-Species Contamination Summary Report",
        help="Report title (appears in <h1> and <title>).",
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

    if args.top_n < 1:
        log.error("--top-n must be >= 1, got %d", args.top_n)
        return 1
    if args.variant_impact_threshold_ppm < 0:
        log.error(
            "--variant-impact-threshold-ppm must be >= 0, got %s",
            args.variant_impact_threshold_ppm,
        )
        return 1

    try:
        inputs = load_inputs(args.aggregate_dir, detect_dir=args.detect_dir)
    except (FileNotFoundError, ValueError) as exc:
        log.error("Failed to load inputs: %s", exc)
        return 1

    try:
        html_path = generate_html_report(
            inputs,
            args.output,
            top_n=args.top_n,
            threshold_ppm=args.variant_impact_threshold_ppm,
            title=args.title,
        )
    except (ValueError, OSError) as exc:
        log.error("Failed to generate report: %s", exc)
        return 1

    print(f"  report: {html_path}")
    print(f"  manifest: {html_path.with_name('report_manifest.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
