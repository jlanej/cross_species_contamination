#!/usr/bin/env python3
"""Regenerate golden output files used by tests/test_golden.py.

Run this script from the repository root whenever pipeline output
intentionally changes::

    python tests/generate_golden.py

Always review the diff before committing updated golden files.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

GOLDEN_DIR = Path(__file__).resolve().parent / "golden"

# Fields that vary between runs and must be excluded from golden comparison.
_VOLATILE_METADATA_KEYS = {"timestamp", "date", "time"}


def _write_report(path: Path, content: str) -> Path:
    path.write_text(content)
    return path


def generate_aggregate_golden() -> None:
    """Produce golden aggregate matrix and metadata."""
    import tempfile

    from csc.aggregate.aggregate import aggregate_reports

    tmpdir = Path(tempfile.mkdtemp())
    d = tmpdir / "reports"
    d.mkdir()

    _write_report(
        d / "sampleA.kraken2.report.txt",
        "50.00\t500\t200\tU\t0\tunclassified\n"
        "50.00\t500\t50\tR\t1\troot\n"
        "20.00\t200\t100\tD\t2\tBacteria\n"
        "10.00\t100\t80\tG\t1279\tStaphylococcus\n"
        "5.00\t50\t40\tS\t1280\tStaphylococcus aureus\n"
        "5.00\t50\t30\tS\t562\tEscherichia coli\n",
    )

    _write_report(
        d / "sampleB.kraken2.report.txt",
        "40.00\t400\t150\tU\t0\tunclassified\n"
        "60.00\t600\t50\tR\t1\troot\n"
        "50.00\t500\t200\tD\t2\tBacteria\n"
        "5.00\t50\t30\tG\t1279\tStaphylococcus\n"
        "3.00\t30\t20\tS\t1280\tStaphylococcus aureus\n"
        "30.00\t300\t250\tS\t562\tEscherichia coli\n",
    )

    out = tmpdir / "agg_out"
    aggregate_reports(
        sorted(d.glob("*.kraken2.report.txt")),
        out,
        normalize=False,
        min_reads=0,
    )

    # Copy matrix
    (GOLDEN_DIR / "aggregate_matrix.tsv").write_text(
        (out / "taxa_matrix.tsv").read_text()
    )

    # Copy metadata (strip volatile fields)
    meta = json.loads((out / "aggregation_metadata.json").read_text())
    meta = {k: v for k, v in meta.items() if k not in _VOLATILE_METADATA_KEYS}
    (GOLDEN_DIR / "aggregate_metadata.json").write_text(
        json.dumps(meta, indent=2) + "\n"
    )

    print("✓ Aggregate golden files updated")


def generate_detect_golden() -> None:
    """Produce golden detect outputs."""
    import tempfile

    from csc.detect.detect import detect_outliers
    from csc.detect.report import generate_report

    tmpdir = Path(tempfile.mkdtemp())

    header = [
        "tax_id", "name",
        "clean_01", "clean_02", "clean_03", "clean_04", "clean_05",
        "contam_01",
    ]
    rows = [
        ["1279", "Staphylococcus", "100", "102", "99", "101", "100", "105"],
        ["562", "Escherichia coli", "50", "51", "49", "50", "52", "5000"],
        ["0", "unclassified", "800", "810", "795", "805", "798", "803"],
    ]

    matrix_path = tmpdir / "matrix.tsv"
    with open(matrix_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for r in rows:
            writer.writerow(r)

    result = detect_outliers(matrix_path, method="mad")
    out = tmpdir / "detect_out"
    generate_report(result, out)

    (GOLDEN_DIR / "detect_flagged_samples.tsv").write_text(
        (out / "flagged_samples.tsv").read_text()
    )
    (GOLDEN_DIR / "detect_qc_summary.json").write_text(
        (out / "qc_summary.json").read_text()
    )
    (GOLDEN_DIR / "detect_quarantine_list.txt").write_text(
        (out / "quarantine_list.txt").read_text()
    )

    print("✓ Detect golden files updated")


if __name__ == "__main__":
    GOLDEN_DIR.mkdir(exist_ok=True)
    generate_aggregate_golden()
    generate_detect_golden()
    print("\nAll golden files regenerated in:", GOLDEN_DIR)
