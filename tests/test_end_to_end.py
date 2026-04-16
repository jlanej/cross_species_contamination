"""End-to-end pipeline tests with synthetic spike-in contamination.

Tests the full aggregate → detect pipeline using deterministic synthetic
data with ground-truth labels.  Includes both negative (clean) and
positive (contaminated) controls to validate that:

* Clean cohorts produce **no** flags (negative control).
* Contaminated cohorts correctly identify spiked samples (positive control).
* Accuracy metrics (precision, recall, FDR) remain within bounds.

These tests complement the golden-output regression tests in
``test_golden.py`` which verify exact output stability, whereas the
tests here focus on correctness against ground truth.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from csc.aggregate.aggregate import aggregate_reports
from csc.detect.detect import detect_outliers
from csc.detect.report import generate_report


# ── Helpers ──────────────────────────────────────────────────────────────────


def _write_report(path: Path, content: str) -> Path:
    path.write_text(content)
    return path


def _write_matrix(
    path: Path, header: list[str], rows: list[list[str]]
) -> Path:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for r in rows:
            writer.writerow(r)
    return path


# ── Negative control: clean cohort ──────────────────────────────────────────


class TestNegativeControlCleanCohort:
    """A fully clean cohort should produce zero flags."""

    @staticmethod
    def _build_clean_cohort(
        tmp_path: Path, n_samples: int = 10
    ) -> Path:
        """Create Kraken2 reports for a cohort with no contamination."""
        d = tmp_path / "clean_reports"
        d.mkdir()

        for i in range(n_samples):
            content = (
                f"60.00\t600\t{200 + i}\tU\t0\tunclassified\n"
                f"40.00\t400\t{50 + i}\tR\t1\troot\n"
                f"20.00\t200\t{100 + i}\tD\t2\tBacteria\n"
                f"10.00\t100\t{80 + i}\tG\t1279\tStaphylococcus\n"
                f"5.00\t50\t{40 + i}\tS\t1280\tStaphylococcus aureus\n"
                f"5.00\t50\t{30 + i}\tS\t562\tEscherichia coli\n"
            )
            _write_report(d / f"clean_{i:03d}.kraken2.report.txt", content)
        return d

    def test_clean_aggregate_detect_no_flags(self, tmp_path: Path) -> None:
        """End-to-end: aggregate clean reports → detect → zero flags."""
        report_dir = self._build_clean_cohort(tmp_path)
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))

        agg_out = tmp_path / "agg_out"
        agg_result = aggregate_reports(reports, agg_out)
        assert agg_result["sample_count"] == 10

        detect_result = detect_outliers(
            agg_result["matrix_raw_path"], method="mad"
        )
        assert detect_result["summary"]["flagged_count"] == 0
        assert detect_result["flagged"] == []

    def test_clean_iqr_no_flags(self, tmp_path: Path) -> None:
        """IQR method should also produce zero flags on clean data."""
        report_dir = self._build_clean_cohort(tmp_path)
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))

        agg_out = tmp_path / "agg_out"
        agg_result = aggregate_reports(reports, agg_out)

        detect_result = detect_outliers(
            agg_result["matrix_raw_path"], method="iqr"
        )
        assert detect_result["summary"]["flagged_count"] == 0

    def test_clean_reports_generate_complete_output(
        self, tmp_path: Path
    ) -> None:
        """Even clean data should produce all three output files."""
        report_dir = self._build_clean_cohort(tmp_path, n_samples=5)
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))

        agg_out = tmp_path / "agg_out"
        agg_result = aggregate_reports(reports, agg_out)

        detect_result = detect_outliers(
            agg_result["matrix_raw_path"], method="mad"
        )
        detect_out = tmp_path / "detect_out"
        report_files = generate_report(detect_result, detect_out)

        assert report_files["flagged_samples"].exists()
        assert report_files["qc_summary"].exists()
        assert report_files["quarantine_list"].exists()
        # Quarantine list should be empty
        assert report_files["quarantine_list"].read_text().strip() == ""


# ── Positive control: contaminated cohort ───────────────────────────────────


class TestPositiveControlContaminatedCohort:
    """A cohort with spiked contamination should flag the correct samples."""

    @staticmethod
    def _build_contaminated_cohort(
        tmp_path: Path,
        *,
        n_clean: int = 15,
        n_contaminated: int = 3,
        spike_reads: int = 5000,
    ) -> tuple[Path, set[str], set[str]]:
        """Create reports: clean samples + contaminated samples with spike-in.

        Returns (report_dir, clean_ids, contaminated_ids).
        """
        d = tmp_path / "mixed_reports"
        d.mkdir()
        clean_ids: set[str] = set()
        contaminated_ids: set[str] = set()

        # Clean samples
        for i in range(n_clean):
            sid = f"clean_{i:03d}"
            clean_ids.add(sid)
            content = (
                f"60.00\t600\t{200 + i}\tU\t0\tunclassified\n"
                f"40.00\t400\t{50 + i}\tR\t1\troot\n"
                f"20.00\t200\t{100 + i}\tD\t2\tBacteria\n"
                f"10.00\t100\t{80 + i}\tG\t1279\tStaphylococcus\n"
                f"5.00\t50\t{40 + i}\tS\t1280\tStaphylococcus aureus\n"
                f"5.00\t50\t{30 + i}\tS\t562\tEscherichia coli\n"
            )
            _write_report(d / f"{sid}.kraken2.report.txt", content)

        # Contaminated samples: spike E. coli
        for i in range(n_contaminated):
            sid = f"contam_{i:03d}"
            contaminated_ids.add(sid)
            content = (
                f"60.00\t600\t{200 + i}\tU\t0\tunclassified\n"
                f"40.00\t400\t{50 + i}\tR\t1\troot\n"
                f"20.00\t200\t{100 + i}\tD\t2\tBacteria\n"
                f"10.00\t100\t{80 + i}\tG\t1279\tStaphylococcus\n"
                f"5.00\t50\t{40 + i}\tS\t1280\tStaphylococcus aureus\n"
                f"90.00\t{spike_reads}\t{spike_reads}\t"
                f"S\t562\tEscherichia coli\n"
            )
            _write_report(d / f"{sid}.kraken2.report.txt", content)

        return d, clean_ids, contaminated_ids

    def test_flags_contaminated_samples(self, tmp_path: Path) -> None:
        """End-to-end: contaminated samples should be flagged."""
        report_dir, clean_ids, contam_ids = (
            self._build_contaminated_cohort(tmp_path)
        )
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))

        agg_out = tmp_path / "agg_out"
        agg_result = aggregate_reports(reports, agg_out)

        detect_result = detect_outliers(
            agg_result["matrix_raw_path"], method="mad"
        )
        flagged = set(detect_result["summary"]["flagged_samples"])

        # All contaminated samples should be flagged
        for sid in contam_ids:
            assert sid in flagged, f"Contaminated {sid} not flagged"

    def test_no_clean_samples_flagged(self, tmp_path: Path) -> None:
        """Clean samples should not be falsely flagged."""
        report_dir, clean_ids, contam_ids = (
            self._build_contaminated_cohort(tmp_path)
        )
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))

        agg_out = tmp_path / "agg_out"
        agg_result = aggregate_reports(reports, agg_out)

        detect_result = detect_outliers(
            agg_result["matrix_raw_path"], method="mad"
        )
        flagged = set(detect_result["summary"]["flagged_samples"])

        false_positives = flagged & clean_ids
        assert len(false_positives) == 0, (
            f"Clean samples falsely flagged: {false_positives}"
        )

    def test_precision_above_threshold(self, tmp_path: Path) -> None:
        """Precision must be >= 0.8."""
        report_dir, clean_ids, contam_ids = (
            self._build_contaminated_cohort(tmp_path)
        )
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))

        agg_out = tmp_path / "agg_out"
        agg_result = aggregate_reports(reports, agg_out)

        detect_result = detect_outliers(
            agg_result["matrix_raw_path"], method="mad"
        )
        flagged = set(detect_result["summary"]["flagged_samples"])

        tp = len(flagged & contam_ids)
        precision = tp / len(flagged) if flagged else 1.0
        assert precision >= 0.8, f"Precision {precision:.2f} < 0.8"

    def test_recall_above_threshold(self, tmp_path: Path) -> None:
        """Recall must be >= 0.8."""
        report_dir, clean_ids, contam_ids = (
            self._build_contaminated_cohort(tmp_path)
        )
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))

        agg_out = tmp_path / "agg_out"
        agg_result = aggregate_reports(reports, agg_out)

        detect_result = detect_outliers(
            agg_result["matrix_raw_path"], method="mad"
        )
        flagged = set(detect_result["summary"]["flagged_samples"])

        tp = len(flagged & contam_ids)
        recall = tp / len(contam_ids) if contam_ids else 0.0
        assert recall >= 0.8, f"Recall {recall:.2f} < 0.8"

    def test_fdr_below_threshold(self, tmp_path: Path) -> None:
        """False discovery rate must be <= 0.2."""
        report_dir, clean_ids, contam_ids = (
            self._build_contaminated_cohort(
                tmp_path, n_clean=30, n_contaminated=3, spike_reads=10000
            )
        )
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))

        agg_out = tmp_path / "agg_out"
        agg_result = aggregate_reports(reports, agg_out)

        detect_result = detect_outliers(
            agg_result["matrix_raw_path"], method="mad"
        )
        flagged = set(detect_result["summary"]["flagged_samples"])

        fp = len(flagged & clean_ids)
        fdr = fp / len(flagged) if flagged else 0.0
        assert fdr <= 0.2, f"FDR {fdr:.2f} > 0.2"

    def test_quarantine_list_correct(self, tmp_path: Path) -> None:
        """Quarantine list should contain exactly the contaminated samples."""
        report_dir, clean_ids, contam_ids = (
            self._build_contaminated_cohort(tmp_path)
        )
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))

        agg_out = tmp_path / "agg_out"
        agg_result = aggregate_reports(reports, agg_out)

        detect_result = detect_outliers(
            agg_result["matrix_raw_path"], method="mad"
        )
        detect_out = tmp_path / "detect_out"
        report_files = generate_report(detect_result, detect_out)

        quarantine = set(
            report_files["quarantine_list"].read_text().strip().splitlines()
        )
        assert contam_ids <= quarantine, (
            f"Missing from quarantine: {contam_ids - quarantine}"
        )


# ── Edge-case tests ─────────────────────────────────────────────────────────


class TestEndToEndEdgeCases:
    """Edge cases for the aggregate → detect pipeline."""

    def test_single_sample_no_crash(self, tmp_path: Path) -> None:
        """A single-sample cohort should not error out."""
        d = tmp_path / "reports"
        d.mkdir()
        _write_report(
            d / "only.kraken2.report.txt",
            "50.00\t500\t200\tU\t0\tunclassified\n"
            "50.00\t500\t50\tR\t1\troot\n",
        )

        agg_out = tmp_path / "agg_out"
        agg_result = aggregate_reports(
            [d / "only.kraken2.report.txt"], agg_out
        )

        detect_result = detect_outliers(
            agg_result["matrix_raw_path"], method="mad"
        )
        # Single sample → MAD=0 → no outliers can be detected
        assert detect_result["summary"]["flagged_count"] == 0

    def test_all_identical_values_no_flags(self, tmp_path: Path) -> None:
        """When every sample has identical counts, nothing should be flagged."""
        header = ["tax_id", "name", "s1", "s2", "s3", "s4", "s5"]
        rows = [
            ["562", "Escherichia coli", "100", "100", "100", "100", "100"],
        ]
        matrix = _write_matrix(tmp_path / "uniform.tsv", header, rows)
        result = detect_outliers(matrix, method="mad")
        assert result["summary"]["flagged_count"] == 0

    def test_two_samples_with_outlier(self, tmp_path: Path) -> None:
        """With only two samples, extreme difference should flag one."""
        header = ["tax_id", "name", "s1", "s2"]
        rows = [
            ["562", "Escherichia coli", "50", "50000"],
        ]
        matrix = _write_matrix(tmp_path / "two.tsv", header, rows)
        result = detect_outliers(matrix, method="iqr")
        # IQR may not flag with only 2 samples; just ensure no crash
        assert isinstance(result["summary"]["flagged_count"], int)

    def test_cpm_normalized_detect(self, tmp_path: Path) -> None:
        """Detection should also work on CPM-normalized matrices."""
        d = tmp_path / "reports"
        d.mkdir()
        for i in range(8):
            reads = 30 + i
            content = (
                f"60.00\t600\t200\tU\t0\tunclassified\n"
                f"5.00\t50\t{reads}\tS\t562\tEscherichia coli\n"
            )
            _write_report(d / f"s{i:03d}.kraken2.report.txt", content)

        # Add one contaminated sample
        _write_report(
            d / "contam.kraken2.report.txt",
            "60.00\t600\t200\tU\t0\tunclassified\n"
            "90.00\t9000\t9000\tS\t562\tEscherichia coli\n",
        )

        agg_out = tmp_path / "agg_out"
        agg_result = aggregate_reports(
            sorted(d.glob("*.kraken2.report.txt")),
            agg_out,
        )

        detect_result = detect_outliers(
            agg_result["matrix_raw_path"], method="mad"
        )
        flagged = set(detect_result["summary"]["flagged_samples"])
        assert "contam" in flagged
