"""Golden-output regression tests.

Compares pipeline module outputs against versioned golden files in
``tests/golden/``.  If any of these tests fail it means the pipeline
behaviour has changed—review the diff and update golden files only when
the change is intentional.

To regenerate golden files run::

    python tests/generate_golden.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from csc.aggregate.aggregate import aggregate_reports
from csc.detect.detect import detect_outliers
from csc.detect.report import generate_report

GOLDEN_DIR = Path(__file__).resolve().parent / "golden"

# ── Helpers ──────────────────────────────────────────────────────────────────


def _write_report(path: Path, content: str) -> Path:
    path.write_text(content)
    return path


def _make_aggregate_reports(tmp_path: Path) -> Path:
    """Create the two deterministic Kraken2 reports used for the golden."""
    d = tmp_path / "reports"
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
    return d


def _make_detect_matrix(tmp_path: Path) -> Path:
    """Create the deterministic contaminated matrix used for the golden."""
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
    p = tmp_path / "matrix.tsv"
    with open(p, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for r in rows:
            writer.writerow(r)
    return p


# ── Aggregate golden tests ──────────────────────────────────────────────────


class TestAggregateGolden:
    """Verify aggregate output matches the versioned golden matrix."""

    def test_matrix_matches_golden(self, tmp_path: Path) -> None:
        report_dir = _make_aggregate_reports(tmp_path)
        out = tmp_path / "agg_out"
        aggregate_reports(
            sorted(report_dir.glob("*.kraken2.report.txt")),
            out,
            normalize=False,
            min_reads=0,
        )

        actual = (out / "taxa_matrix.tsv").read_text()
        golden = (GOLDEN_DIR / "aggregate_matrix.tsv").read_text()
        assert actual == golden, (
            "Aggregate matrix differs from golden file.  "
            "If the change is intentional, regenerate with: "
            "python tests/generate_golden.py"
        )

    def test_metadata_matches_golden(self, tmp_path: Path) -> None:
        report_dir = _make_aggregate_reports(tmp_path)
        out = tmp_path / "agg_out"
        aggregate_reports(
            sorted(report_dir.glob("*.kraken2.report.txt")),
            out,
            normalize=False,
            min_reads=0,
        )

        actual = json.loads((out / "aggregation_metadata.json").read_text())
        golden = json.loads(
            (GOLDEN_DIR / "aggregate_metadata.json").read_text()
        )

        # Remove volatile fields that change across runs
        for key in ("timestamp", "date", "time"):
            actual.pop(key, None)
            golden.pop(key, None)

        assert actual == golden, (
            "Aggregate metadata differs from golden file."
        )


# ── Detect golden tests ─────────────────────────────────────────────────────


class TestDetectGolden:
    """Verify detect output matches the versioned golden files."""

    def test_qc_summary_matches_golden(self, tmp_path: Path) -> None:
        matrix = _make_detect_matrix(tmp_path)
        result = detect_outliers(matrix, method="mad")
        out = tmp_path / "detect_out"
        generate_report(result, out)

        actual = json.loads((out / "qc_summary.json").read_text())
        golden = json.loads(
            (GOLDEN_DIR / "detect_qc_summary.json").read_text()
        )
        assert actual == golden, (
            "Detect QC summary differs from golden file."
        )

    def test_flagged_samples_matches_golden(self, tmp_path: Path) -> None:
        matrix = _make_detect_matrix(tmp_path)
        result = detect_outliers(matrix, method="mad")
        out = tmp_path / "detect_out"
        generate_report(result, out)

        actual = (out / "flagged_samples.tsv").read_text()
        golden = (GOLDEN_DIR / "detect_flagged_samples.tsv").read_text()
        assert actual == golden, (
            "Detect flagged samples TSV differs from golden file."
        )

    def test_quarantine_list_matches_golden(self, tmp_path: Path) -> None:
        matrix = _make_detect_matrix(tmp_path)
        result = detect_outliers(matrix, method="mad")
        out = tmp_path / "detect_out"
        generate_report(result, out)

        actual = (out / "quarantine_list.txt").read_text()
        golden = (GOLDEN_DIR / "detect_quarantine_list.txt").read_text()
        assert actual == golden, (
            "Detect quarantine list differs from golden file."
        )
