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
        "--layout",
        choices=("cohort", "legacy"),
        default="cohort",
        help=(
            "Report layout.  'cohort' (default) is the species-centric "
            "layout suitable for cohorts of 3K+ samples.  'legacy' "
            "preserves the per-sample manuscript layout for one release "
            "window for byte-level diffing."
        ),
    )
    parser.add_argument(
        "--page-size",
        type=int,
        default=25,
        help=(
            "Default page size for paginated tables (species summary, "
            "variant-impact flagged samples, per-sample appendix).  Cohort "
            "layout only.  Default: 25."
        ),
    )
    parser.add_argument(
        "--top-species",
        type=int,
        default=50,
        help=(
            "Top-K species (by cohort burden) shown in the heatmap "
            "and used for PCoA distance computation.  Cohort layout "
            "only.  Default: 50."
        ),
    )
    parser.add_argument(
        "--drilldown-top",
        type=int,
        default=25,
        help=(
            "Number of species expanded as drill-down details in §3.8.  "
            "Cohort layout only.  Default: 25."
        ),
    )
    parser.add_argument(
        "--cluster-method",
        choices=("average", "single", "ward"),
        default="average",
        help=(
            "Hierarchical-clustering linkage method for §3.6 / §3.7.  "
            "Default: average (UPGMA)."
        ),
    )
    parser.add_argument(
        "--cluster-distance",
        choices=("bray", "jaccard"),
        default="bray",
        help=(
            "Pairwise distance for hierarchical clustering and PCoA.  "
            "'bray' (Bray–Curtis on log1p-CPM, default) is the standard "
            "metagenomic choice; 'jaccard' is presence/absence based."
        ),
    )
    parser.add_argument(
        "--prevalence-core",
        type=float,
        default=0.5,
        help=(
            "Minimum prevalence (fraction of samples) for a species to "
            "be classified as 'core' in §3.4.  Default: 0.5."
        ),
    )
    parser.add_argument(
        "--prevalence-rare",
        type=float,
        default=0.1,
        help=(
            "Maximum prevalence (fraction of samples) for a species to "
            "be classified as 'rare' in §3.4 (species below this and "
            "above 0 are 'rare', between this and --prevalence-core are "
            "'accessory').  Default: 0.1."
        ),
    )
    parser.add_argument(
        "--max-samples-cluster",
        type=int,
        default=2000,
        help=(
            "Cap on samples included in O(n²) pairwise-distance "
            "calculations (heatmap + PCoA).  When the cohort exceeds "
            "this cap a deterministic every-k-th sub-sample is taken.  "
            "Default: 2000."
        ),
    )
    parser.add_argument(
        "--min-reads-for-prevalence",
        type=int,
        default=5,
        help=(
            "Per-taxon minimum read count for the second prevalence "
            "column in §3.1.  Default: 5."
        ),
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
    if args.page_size < 1:
        log.error("--page-size must be >= 1, got %d", args.page_size)
        return 1
    if args.top_species < 1:
        log.error("--top-species must be >= 1, got %d", args.top_species)
        return 1
    if args.max_samples_cluster < 3:
        log.error(
            "--max-samples-cluster must be >= 3, got %d",
            args.max_samples_cluster,
        )
        return 1
    if not (0 < args.prevalence_rare <= args.prevalence_core <= 1):
        log.error(
            "Need 0 < --prevalence-rare <= --prevalence-core <= 1, got "
            "%s / %s", args.prevalence_rare, args.prevalence_core,
        )
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
            layout=args.layout,
            page_size=args.page_size,
            top_species=args.top_species,
            drilldown_top=args.drilldown_top,
            cluster_method=args.cluster_method,
            cluster_distance=args.cluster_distance,
            prevalence_core=args.prevalence_core,
            prevalence_rare=args.prevalence_rare,
            max_samples_cluster=args.max_samples_cluster,
            min_reads_for_prevalence=args.min_reads_for_prevalence,
        )
    except (ValueError, OSError) as exc:
        log.error("Failed to generate report: %s", exc)
        return 1

    print(f"  report: {html_path}")
    print(f"  manifest: {html_path.with_name('report_manifest.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
