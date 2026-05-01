#!/usr/bin/env python3
"""Generate comprehensive synthetic data for the HTML report pipeline.

This module creates a realistic, varied cohort of 20 synthetic samples
spanning multiple contamination scenarios:

* **clean** – predominantly human; negligible non-human burden
* **moderate** – bacterial background at levels commonly seen in WGS labs
* **heavy** – high bacterial/fungal/viral contamination
* **mixed** – multiple co-contaminating species

All domain categories surfaced by the report (Bacteria, Archaea, Fungi,
Viruses, Metazoa_other, Viridiplantae, Human) are represented so that
every section of the cohort-layout report is exercised.

The function :func:`generate_report_test_data` drives the full pipeline:

1. Writes per-sample Kraken2-format ``.kraken2.report.txt`` files.
2. Writes ``reads_summary.json`` sidecars so the absolute-burden matrix is
   populated.
3. Calls :func:`csc.aggregate.aggregate.aggregate_reports` to build the
   aggregate directory.
4. Calls the ``csc-detect`` CLI to produce detection outputs.
5. Calls the ``csc-report`` CLI to render ``docs/index.html``.

It can also be run as a standalone script::

    python tests/generate_report_test_data.py

which writes output to ``/tmp/csc_demo_report/`` and renders the final
HTML to ``docs/index.html`` (relative to the repository root).
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Taxonomy catalogue
# ---------------------------------------------------------------------------

# Each entry: (tax_id, name, rank_code, domain)
# The catalogue covers all domain categories the report can display.
_TAXA: list[tuple[int, str, str, str]] = [
    # Human
    (9606, "Homo sapiens", "S", "Human"),
    # Bacteria – common lab contaminants / environmental
    (1280, "Staphylococcus aureus", "S", "Bacteria"),
    (562, "Escherichia coli", "S", "Bacteria"),
    (1423, "Bacillus subtilis", "S", "Bacteria"),
    (1352, "Enterococcus faecium", "S", "Bacteria"),
    (287, "Pseudomonas aeruginosa", "S", "Bacteria"),
    (1314, "Streptococcus pyogenes", "S", "Bacteria"),
    (210, "Helicobacter pylori", "S", "Bacteria"),
    (485, "Neisseria gonorrhoeae", "S", "Bacteria"),
    (727, "Haemophilus influenzae", "S", "Bacteria"),
    (1639, "Listeria monocytogenes", "S", "Bacteria"),
    (1773, "Mycobacterium tuberculosis", "S", "Bacteria"),
    (194, "Campylobacter jejuni", "S", "Bacteria"),
    (1301, "Streptococcus pneumoniae", "S", "Bacteria"),
    (620, "Klebsiella pneumoniae", "S", "Bacteria"),
    # Genera (G)
    (1279, "Staphylococcus", "G", "Bacteria"),
    (561, "Escherichia", "G", "Bacteria"),
    (1386, "Bacillus", "G", "Bacteria"),
    # Archaea
    (2190, "Methanobrevibacter smithii", "S", "Archaea"),
    (187420, "Nitrosopumilus maritimus", "S", "Archaea"),
    # Fungi
    (5207, "Cryptococcus neoformans", "S", "Fungi"),
    (4932, "Saccharomyces cerevisiae", "S", "Fungi"),
    (5476, "Candida albicans", "S", "Fungi"),
    (5234, "Aspergillus fumigatus", "S", "Fungi"),
    # Viruses
    (11234, "Measles morbillivirus", "S", "Viruses"),
    (12227, "Human adenovirus 5", "S", "Viruses"),
    (10359, "Human cytomegalovirus", "S", "Viruses"),
    (10376, "Epstein-Barr virus", "S", "Viruses"),
    # Metazoa (non-human)
    (10090, "Mus musculus", "S", "Metazoa_other"),
    (10116, "Rattus norvegicus", "S", "Metazoa_other"),
    (9986, "Oryctolagus cuniculus", "S", "Metazoa_other"),
    # Plants
    (3702, "Arabidopsis thaliana", "S", "Viridiplantae"),
    (4577, "Zea mays", "S", "Viridiplantae"),
    # UniVec / synthetic
    (81077, "artificial sequences", "S", "UniVec_Core"),
    # Unclassified (tax_id 0 – always present)
    (0, "unclassified", "U", "Unclassified"),
    # Root / higher-rank nodes
    (1, "root", "R", "Unclassified"),
    (2, "Bacteria", "D", "Bacteria"),
    (2157, "Archaea", "D", "Archaea"),
    (4751, "Fungi", "D", "Fungi"),
    (10239, "Viruses", "D", "Viruses"),
    (33208, "Metazoa", "D", "Metazoa_other"),
]

# Index by tax_id for quick lookup
_TAXA_BY_ID: dict[int, tuple[int, str, str, str]] = {t[0]: t for t in _TAXA}


# ---------------------------------------------------------------------------
# Scenario definitions
# ---------------------------------------------------------------------------

# Each scenario is a dict of {tax_id: direct_reads}.  The _build_report
# function fills clade_reads from parent→child relationships (simplified:
# clade_reads = direct_reads for leaves; higher ranks accumulate).

_SCENARIOS: dict[str, dict[int, int]] = {
    # --- Clean samples (low non-human) ---
    "clean_A": {9606: 950_000, 1280: 30, 562: 10, 5476: 5, 0: 5_000, 1: 50_000},
    "clean_B": {9606: 970_000, 562: 15, 1423: 5, 0: 3_000, 1: 30_000},
    "clean_C": {9606: 990_000, 1280: 5, 0: 2_000, 1: 10_000},
    "clean_D": {9606: 945_000, 1301: 20, 727: 10, 0: 4_500, 1: 55_000},
    # --- Moderate bacterial contamination ---
    "mod_staph": {9606: 800_000, 1280: 5_000, 562: 1_000, 1279: 6_000,
                  0: 20_000, 1: 200_000, 2: 12_000},
    "mod_ecoli": {9606: 820_000, 562: 8_000, 561: 9_000, 1423: 500,
                  0: 15_000, 1: 180_000, 2: 10_000},
    "mod_mixed_bact": {9606: 790_000, 1280: 3_000, 562: 3_000, 620: 2_000,
                       287: 1_000, 1314: 500, 0: 25_000, 1: 210_000, 2: 10_000},
    "mod_mycobact": {9606: 850_000, 1773: 4_000, 1352: 800, 0: 10_000,
                     1: 140_000, 2: 5_000},
    # --- Heavy bacterial contamination ---
    "heavy_staph": {9606: 500_000, 1280: 60_000, 1279: 65_000, 562: 5_000,
                    0: 50_000, 1: 500_000, 2: 70_000},
    "heavy_ecoli": {9606: 480_000, 562: 80_000, 561: 85_000, 1423: 3_000,
                    194: 1_000, 0: 40_000, 1: 490_000, 2: 90_000},
    "heavy_klebsiella": {9606: 460_000, 620: 70_000, 562: 20_000, 287: 10_000,
                         0: 60_000, 1: 470_000, 2: 100_000},
    # --- Multi-domain contamination ---
    "fungal_contam": {9606: 700_000, 5476: 15_000, 4932: 5_000, 5207: 3_000,
                      1280: 2_000, 0: 30_000, 1: 290_000, 4751: 23_000},
    "viral_contam": {9606: 750_000, 10376: 8_000, 10359: 5_000, 11234: 2_000,
                     1280: 1_000, 0: 20_000, 1: 240_000, 10239: 15_000},
    "mouse_contam": {9606: 680_000, 10090: 40_000, 10116: 5_000,
                     1280: 2_000, 562: 1_000, 0: 30_000, 1: 310_000,
                     33208: 45_000},
    "archaea_contam": {9606: 780_000, 2190: 12_000, 187420: 5_000,
                       1280: 1_000, 0: 25_000, 1: 220_000, 2157: 17_000},
    "plant_contam": {9606: 720_000, 3702: 20_000, 4577: 8_000,
                     562: 500, 0: 35_000, 1: 280_000},
    "synthetic_contam": {9606: 810_000, 81077: 10_000, 1280: 3_000,
                         0: 20_000, 1: 190_000},
    # --- Heavy multi-domain ---
    "zoo_mix": {9606: 400_000, 10090: 50_000, 1280: 30_000, 562: 20_000,
                5476: 10_000, 10376: 5_000, 3702: 5_000, 0: 80_000,
                1: 600_000, 2: 60_000, 33208: 55_000, 4751: 10_000,
                10239: 5_000},
    "super_contaminated": {9606: 200_000, 1280: 100_000, 562: 80_000,
                           620: 50_000, 287: 30_000, 5476: 20_000,
                           10090: 15_000, 10376: 10_000, 0: 100_000,
                           1: 800_000, 2: 260_000, 4751: 20_000,
                           10239: 10_000, 33208: 15_000},
    "outlier_single_species": {9606: 600_000, 1280: 120_000, 1279: 125_000,
                                0: 40_000, 1: 400_000, 2: 130_000},
}

# Assign deterministic total_reads per scenario (mapped + unmapped ≈ 1 M each)
_TOTAL_READS: dict[str, tuple[int, int]] = {
    "clean_A":               (950_100, 50_000),
    "clean_B":               (970_030, 30_000),
    "clean_C":               (990_020, 10_000),
    "clean_D":               (945_080, 55_000),
    "mod_staph":             (830_000, 170_000),
    "mod_ecoli":             (840_000, 160_000),
    "mod_mixed_bact":        (825_000, 175_000),
    "mod_mycobact":          (860_000, 140_000),
    "heavy_staph":           (620_000, 380_000),
    "heavy_ecoli":           (600_000, 400_000),
    "heavy_klebsiella":      (590_000, 410_000),
    "fungal_contam":         (735_000, 265_000),
    "viral_contam":          (766_000, 234_000),
    "mouse_contam":          (728_000, 272_000),
    "archaea_contam":        (799_000, 201_000),
    "plant_contam":          (749_000, 251_000),
    "synthetic_contam":      (824_000, 176_000),
    "zoo_mix":               (620_000, 380_000),
    "super_contaminated":    (605_000, 395_000),
    "outlier_single_species": (745_000, 255_000),
}

# Higher-rank taxa that aggregate clade_reads from child species
_PARENT_CHILDREN: dict[int, list[int]] = {
    1: list({t[0] for t in _TAXA if t[0] not in (0, 1)}),
    2: [1280, 562, 1423, 1352, 287, 1314, 210, 485, 727, 1639, 1773, 194,
        1301, 620, 1279, 561, 1386],
    4751: [5207, 4932, 5476, 5234],
    10239: [11234, 12227, 10359, 10376],
    33208: [10090, 10116, 9986, 9606],
    2157: [2190, 187420],
    1279: [1280],
    561: [562],
    1386: [1423],
}


# ---------------------------------------------------------------------------
# Kraken2 report builder
# ---------------------------------------------------------------------------

def _build_report_lines(
    direct_reads: dict[int, int],
    total_reads: int,
) -> str:
    """Produce a Kraken2-format report from a {tax_id: direct_reads} map.

    Clade reads for parent nodes are accumulated from child direct reads.
    The Kraken2 report format is::

        pct  clade_reads  direct_reads  rank_code  tax_id  name
    """
    # Compute clade reads (parent = sum of own + all child clade)
    clade: dict[int, int] = dict(direct_reads)
    # BFS from children up to root, accumulate
    for parent_id, children in _PARENT_CHILDREN.items():
        child_sum = sum(clade.get(c, 0) for c in children)
        if parent_id in clade:
            clade[parent_id] = clade[parent_id] + child_sum
        elif child_sum > 0:
            clade[parent_id] = child_sum

    lines: list[str] = []
    for tid, entry in _TAXA_BY_ID.items():
        dr = direct_reads.get(tid, 0)
        cr = clade.get(tid, 0)
        if dr == 0 and cr == 0:
            continue
        _, name, rank_code, _ = entry
        pct = cr / total_reads * 100 if total_reads > 0 else 0.0
        indent = "  " if rank_code not in ("U", "R", "D") else ""
        lines.append(
            f"{pct:.2f}\t{cr}\t{dr}\t{rank_code}\t{tid}\t{indent}{name}"
        )
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_report_test_data(
    output_dir: str | Path | None = None,
    report_html: str | Path | None = None,
    *,
    seed: int = 42,
    run_detect: bool = True,
) -> dict[str, Path]:
    """Generate comprehensive mock data and render ``docs/index.html``.

    Parameters
    ----------
    output_dir:
        Directory for intermediate pipeline artefacts (Kraken2 reports,
        aggregate outputs, detect outputs).  Defaults to
        ``/tmp/csc_demo_report``.
    report_html:
        Destination path for the rendered HTML.  Defaults to
        ``<repo_root>/docs/index.html``.
    seed:
        Random seed for reproducibility (currently unused, data is fully
        deterministic).
    run_detect:
        If ``True`` (default), also run ``csc-detect`` so the report
        includes detection-summary sections.

    Returns
    -------
    dict
        Keys: ``aggregate_dir``, ``detect_dir`` (or ``None``),
        ``report_html``.
    """
    random.seed(seed)

    # --- Resolve paths -------------------------------------------------------
    if output_dir is None:
        output_dir = Path("/tmp/csc_demo_report")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if report_html is None:
        # Locate the repo root relative to this file
        repo_root = Path(__file__).resolve().parent.parent
        report_html = repo_root / "docs" / "index.html"
    report_html = Path(report_html)
    report_html.parent.mkdir(parents=True, exist_ok=True)

    # --- Write Kraken2 reports + reads_summary sidecars ----------------------
    kraken_dir = output_dir / "kraken_reports"
    kraken_dir.mkdir(exist_ok=True)
    sidecar_paths: list[Path] = []

    for scenario_name, direct_reads in _SCENARIOS.items():
        total_reads_tuple = _TOTAL_READS[scenario_name]
        mapped, unmapped = total_reads_tuple
        total = mapped + unmapped

        report_text = _build_report_lines(direct_reads, total)
        report_path = kraken_dir / f"{scenario_name}.kraken2.report.txt"
        report_path.write_text(report_text)

        sidecar = {
            "sample_id": scenario_name,
            "total_reads": total,
            "total_mapped": mapped,
            "total_unmapped": unmapped,
            "input": f"/fake/{scenario_name}.bam",
            "extraction_time": "2025-01-15T12:00:00Z",
        }
        sidecar_path = output_dir / f"{scenario_name}.reads_summary.json"
        sidecar_path.write_text(json.dumps(sidecar, indent=2))
        sidecar_paths.append(sidecar_path)

    # --- Aggregate -----------------------------------------------------------
    from csc.aggregate.aggregate import aggregate_reports

    agg_dir = output_dir / "aggregate_out"
    aggregate_reports(
        sorted(kraken_dir.glob("*.kraken2.report.txt")),
        agg_dir,
        idxstats_paths=sidecar_paths,
    )

    # --- Detect (optional) ---------------------------------------------------
    detect_dir: Path | None = None
    if run_detect:
        from csc.detect.cli import main as detect_main

        detect_dir = output_dir / "detect_out"
        rc = detect_main([
            str(agg_dir / "taxa_matrix_cpm.tsv"),
            "-o", str(detect_dir),
        ])
        if rc != 0:
            import warnings
            warnings.warn(
                f"csc-detect exited with code {rc}; detection section "
                "will be omitted from the report."
            )
            detect_dir = None

    # --- Report --------------------------------------------------------------
    from csc.report.cli import main as report_main

    argv = [
        str(agg_dir),
        "-o", str(report_html),
        "--title", "Cross-Species Contamination – Demo Cohort Report",
        "--layout", "cohort",
    ]
    if detect_dir is not None:
        argv += ["--detect-dir", str(detect_dir)]

    rc = report_main(argv)
    if rc != 0:
        raise RuntimeError(f"csc-report exited with code {rc}")

    return {
        "aggregate_dir": agg_dir,
        "detect_dir": detect_dir,
        "report_html": report_html,
    }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate comprehensive mock data and render docs/index.html."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Intermediate artefacts directory (default: /tmp/csc_demo_report).",
    )
    parser.add_argument(
        "--report-html",
        type=Path,
        default=None,
        help="Destination HTML path (default: docs/index.html in repo root).",
    )
    parser.add_argument(
        "--no-detect",
        action="store_true",
        help="Skip csc-detect (report will lack the detection-summary section).",
    )
    args = parser.parse_args()

    paths = generate_report_test_data(
        output_dir=args.output_dir,
        report_html=args.report_html,
        run_detect=not args.no_detect,
    )
    print(f"aggregate_dir : {paths['aggregate_dir']}")
    print(f"detect_dir    : {paths['detect_dir']}")
    print(f"report_html   : {paths['report_html']}")
    sys.exit(0)
