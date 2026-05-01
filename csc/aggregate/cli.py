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
        "--kraken2-output",
        nargs="+",
        type=Path,
        default=None,
        metavar="OUT",
        help=(
            "Paths to per-read Kraken2 output files "
            "(*.kraken2.output.txt).  Required when --confidence-threshold "
            "is used with any value > 0.0; ignored otherwise.  Files are "
            "matched to samples by filename (e.g. <sample_id>.kraken2.output.txt)."
        ),
    )
    parser.add_argument(
        "--confidence-threshold",
        nargs="+",
        type=float,
        default=None,
        metavar="T",
        help=(
            "Kraken2 confidence cutoff(s) in [0.0, 1.0].  For each "
            "threshold T > 0.0, a parallel 'high-confidence' matrix "
            "set is written with filenames suffixed _conf{T} (e.g. "
            "taxa_matrix_raw_conf0p10.tsv).  Requires --kraken2-output "
            "and --db-path.  Multiple values yield multiple tiers (e.g. "
            "--confidence-threshold 0.1 0.5).  When omitted the value "
            "from default_config.yaml ([0.1] by default) is used; pass "
            "--no-confidence-tiers to force sensitive-only output."
        ),
    )
    parser.add_argument(
        "--no-confidence-tiers",
        action="store_true",
        help=(
            "Disable the high-confidence tier(s) configured in "
            "default_config.yaml (aggregate.confidence_thresholds).  "
            "Equivalent to --confidence-threshold (no value)."
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
        # ── Resolve confidence-tier defaults from default_config.yaml ──
        # If the user neither supplied --confidence-threshold nor
        # --no-confidence-tiers, fall back to the YAML default (0.1
        # ships as the new default).  We only honour the default when
        # the dependencies (kraken2 outputs + db-path) are available;
        # otherwise we log a one-line warning and proceed sensitive-only
        # so that ad-hoc CLI invocations without those inputs continue
        # to work without surprising failures.
        if args.confidence_threshold is None and not args.no_confidence_tiers:
            try:
                from csc.config import load_config
                cfg = load_config()
                cfg_thresholds = (
                    cfg.get("aggregate", {}).get("confidence_thresholds") or []
                )
            except Exception as cfg_exc:  # pragma: no cover - defensive
                log.debug("Could not load config defaults: %s", cfg_exc)
                cfg_thresholds = []
            if cfg_thresholds:
                if args.kraken2_output and args.db_path is not None:
                    args.confidence_threshold = [float(t) for t in cfg_thresholds]
                    log.info(
                        "Defaulting --confidence-threshold to %s from "
                        "default_config.yaml (pass --no-confidence-tiers "
                        "to disable).", args.confidence_threshold,
                    )
                else:
                    log.warning(
                        "Skipping default high-confidence tier %s: requires "
                        "--db-path and --kraken2-output.  Pass "
                        "--no-confidence-tiers to silence this warning.",
                        cfg_thresholds,
                    )

        # Validate confidence-tier dependencies fail-early.
        if args.confidence_threshold:
            high_thresholds = [t for t in args.confidence_threshold if t > 0.0]
            if high_thresholds and not args.kraken2_output:
                log.error(
                    "--confidence-threshold > 0.0 requires --kraken2-output"
                )
                return 1
            if high_thresholds and args.db_path is None:
                log.error(
                    "--confidence-threshold > 0.0 requires --db-path "
                    "(needs taxonomy/nodes.dmp for confidence calculation)"
                )
                return 1
            missing_outputs = [
                p for p in (args.kraken2_output or []) if not p.exists()
            ]
            if missing_outputs:
                for p in missing_outputs:
                    log.error("Kraken2 output file not found: %s", p)
                return 1

        result = aggregate_reports(
            args.input,
            args.output_dir,
            min_reads=args.min_reads,
            chunk_size=args.chunk_size,
            rank_filter=tuple(args.rank_filter),
            db_path=args.db_path,
            idxstats_paths=args.idxstats,
            confidence_thresholds=args.confidence_threshold,
            kraken2_output_paths=args.kraken2_output,
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
        for tier_suffix, tier in sorted(
            (result.get("confidence_tiers") or {}).items()
        ):
            print(
                f"  confidence tier {tier_suffix} "
                f"(threshold={tier['threshold']}):"
            )
            print(f"    matrix (raw): {tier['matrix_raw_path']}")
            print(f"    matrix (cpm): {tier['matrix_cpm_path']}")
            if tier.get("matrix_abs_path"):
                print(f"    matrix (abs): {tier['matrix_abs_path']}")
        return 0
    except Exception as exc:
        log.error("Aggregation failed: %s", exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
