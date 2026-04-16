"""Tests for the outlier-detection module (csc.detect).

Covers:
- Matrix loading and validation
- Statistical helpers (median, MAD, quartiles)
- Kitome/background filtering and population-mean subtraction
- MAD and IQR outlier detection on synthetic data
- Report generation (flagged table, QC summary, quarantine list)
- CLI entry point
- Precision/recall on simulated contaminated vs. clean datasets
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from unittest import mock

import pytest

from csc.detect.detect import (
    DetectionResult,
    FlaggedSample,
    _mad,
    _median,
    _quartiles,
    detect_outliers,
    filter_kitome,
    load_matrix,
    subtract_population_background,
)
from csc.detect.report import (
    generate_report,
    write_flagged_samples,
    write_qc_summary,
    write_quarantine_list,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_matrix(path: Path, header: list[str], rows: list[list[str]]) -> Path:
    """Write a TSV matrix file."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)
    return path


@pytest.fixture()
def clean_matrix(tmp_path: Path) -> Path:
    """Matrix where all samples have similar values (no outliers)."""
    header = ["tax_id", "name", "sA", "sB", "sC", "sD", "sE"]
    rows = [
        ["1279", "Staphylococcus", "100", "105", "98", "102", "101"],
        ["562", "Escherichia coli", "50", "48", "52", "49", "51"],
        ["0", "unclassified", "800", "810", "795", "805", "798"],
    ]
    return _write_matrix(tmp_path / "clean.tsv", header, rows)


@pytest.fixture()
def contaminated_matrix(tmp_path: Path) -> Path:
    """Matrix with one sample (sC) spiked in one taxon (E. coli)."""
    header = ["tax_id", "name", "sA", "sB", "sC", "sD", "sE"]
    rows = [
        ["1279", "Staphylococcus", "100", "105", "98", "102", "101"],
        ["562", "Escherichia coli", "50", "48", "5000", "49", "51"],
        ["0", "unclassified", "800", "810", "795", "805", "798"],
    ]
    return _write_matrix(tmp_path / "contaminated.tsv", header, rows)


@pytest.fixture()
def contaminated_iqr_matrix(tmp_path: Path) -> Path:
    """Matrix with enough samples for IQR to work (sC is spiked in E. coli)."""
    samples = [f"s{i}" for i in range(10)]
    header = ["tax_id", "name"] + samples
    # Background taxon – uniform
    staph_vals = ["100", "105", "98", "102", "101", "103", "99", "104", "97", "100"]
    # E. coli: sC (index 2) is spiked
    ecoli_vals = ["50", "48", "5000", "49", "51", "52", "47", "50", "53", "48"]
    rows = [
        ["1279", "Staphylococcus"] + staph_vals,
        ["562", "Escherichia coli"] + ecoli_vals,
    ]
    return _write_matrix(tmp_path / "contaminated_iqr.tsv", header, rows)


@pytest.fixture()
def kitome_matrix(tmp_path: Path) -> Path:
    """Matrix with a kitome taxon that should be excluded."""
    header = ["tax_id", "name", "sA", "sB", "sC"]
    rows = [
        ["1279", "Staphylococcus", "100", "105", "98"],
        ["562", "Escherichia coli", "50", "48", "5000"],
        ["9606", "Homo sapiens", "900", "910", "920"],
    ]
    return _write_matrix(tmp_path / "kitome.tsv", header, rows)


# ---------------------------------------------------------------------------
# Unit tests — statistical helpers
# ---------------------------------------------------------------------------

class TestStatisticalHelpers:
    def test_median_odd(self) -> None:
        assert _median([3, 1, 2]) == 2.0

    def test_median_even(self) -> None:
        assert _median([4, 1, 3, 2]) == 2.5

    def test_median_single(self) -> None:
        assert _median([42]) == 42.0

    def test_median_empty(self) -> None:
        assert _median([]) == 0.0

    def test_mad_constant(self) -> None:
        assert _mad([5, 5, 5, 5]) == 0.0

    def test_mad_symmetric(self) -> None:
        # values: 1,2,3,4,5 → median=3, deviations=2,1,0,1,2 → MAD=1
        assert _mad([1, 2, 3, 4, 5]) == 1.0

    def test_quartiles_odd(self) -> None:
        q1, med, q3 = _quartiles([1, 2, 3, 4, 5])
        assert med == 3.0
        assert q1 == 1.5
        assert q3 == 4.5

    def test_quartiles_even(self) -> None:
        q1, med, q3 = _quartiles([1, 2, 3, 4])
        assert med == 2.5
        assert q1 == 1.5
        assert q3 == 3.5

    def test_quartiles_empty(self) -> None:
        assert _quartiles([]) == (0.0, 0.0, 0.0)


# ---------------------------------------------------------------------------
# Unit tests — matrix loading
# ---------------------------------------------------------------------------

class TestLoadMatrix:
    def test_basic_load(self, clean_matrix: Path) -> None:
        sids, rows, names = load_matrix(clean_matrix)
        assert sids == ["sA", "sB", "sC", "sD", "sE"]
        assert len(rows) == 3
        assert names[1279] == "Staphylococcus"
        assert names[562] == "Escherichia coli"

    def test_values_parsed(self, clean_matrix: Path) -> None:
        sids, rows, _ = load_matrix(clean_matrix)
        staph = rows[0]
        assert staph["tax_id"] == 1279
        assert staph["sA"] == 100.0
        assert staph["sB"] == 105.0

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="not found"):
            load_matrix(tmp_path / "nope.tsv")

    def test_empty_matrix_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "empty.tsv"
        p.write_text("tax_id\tname\tsA\n")
        with pytest.raises(ValueError, match="No valid taxon rows"):
            load_matrix(p)

    def test_too_few_columns_raises(self, tmp_path: Path) -> None:
        p = tmp_path / "bad.tsv"
        p.write_text("tax_id\tname\n1\tFoo\n")
        with pytest.raises(ValueError, match="at least 3 columns"):
            load_matrix(p)


# ---------------------------------------------------------------------------
# Unit tests — kitome filtering
# ---------------------------------------------------------------------------

class TestKitomeFilter:
    def test_no_kitome_returns_all(self) -> None:
        rows = [{"tax_id": 1}, {"tax_id": 2}]
        assert filter_kitome(rows, None) is rows
        assert filter_kitome(rows, []) is rows

    def test_removes_kitome_taxa(self) -> None:
        rows = [{"tax_id": 1}, {"tax_id": 9606}, {"tax_id": 2}]
        filtered = filter_kitome(rows, [9606])
        assert len(filtered) == 2
        assert all(r["tax_id"] != 9606 for r in filtered)

    def test_removes_multiple(self) -> None:
        rows = [{"tax_id": 1}, {"tax_id": 9606}, {"tax_id": 562}]
        filtered = filter_kitome(rows, [9606, 562])
        assert len(filtered) == 1


# ---------------------------------------------------------------------------
# Unit tests — background subtraction
# ---------------------------------------------------------------------------

class TestBackgroundSubtraction:
    def test_subtracts_mean(self) -> None:
        rows = [{"tax_id": 1, "sA": 10.0, "sB": 20.0, "sC": 30.0}]
        result = subtract_population_background(rows, ["sA", "sB", "sC"])
        # mean = 20, so: -10, 0, 10
        assert result[0]["sA"] == pytest.approx(-10.0)
        assert result[0]["sB"] == pytest.approx(0.0)
        assert result[0]["sC"] == pytest.approx(10.0)

    def test_does_not_mutate_original(self) -> None:
        rows = [{"tax_id": 1, "sA": 10.0, "sB": 20.0}]
        subtract_population_background(rows, ["sA", "sB"])
        assert rows[0]["sA"] == 10.0  # unchanged

    def test_preserves_negative_values(self) -> None:
        rows = [{"tax_id": 1, "sA": 1.0, "sB": 100.0}]
        result = subtract_population_background(rows, ["sA", "sB"])
        # mean = 50.5, sA = 1-50.5 = -49.5
        assert result[0]["sA"] == pytest.approx(-49.5)


# ---------------------------------------------------------------------------
# Unit tests — outlier detection (MAD)
# ---------------------------------------------------------------------------

class TestDetectOutliersMAD:
    def test_clean_no_flags(self, clean_matrix: Path) -> None:
        result = detect_outliers(clean_matrix, method="mad")
        assert result["summary"]["flagged_count"] == 0
        assert result["flagged"] == []

    def test_contaminated_flags_spike(self, contaminated_matrix: Path) -> None:
        result = detect_outliers(contaminated_matrix, method="mad")
        flagged_ids = {f["sample_id"] for f in result["flagged"]}
        assert "sC" in flagged_ids
        # Only sC should be flagged for E. coli
        ecoli_flags = [f for f in result["flagged"] if f["tax_id"] == 562]
        assert len(ecoli_flags) >= 1
        assert ecoli_flags[0]["sample_id"] == "sC"

    def test_result_types(self, contaminated_matrix: Path) -> None:
        result = detect_outliers(contaminated_matrix, method="mad")
        assert "flagged" in result
        assert "summary" in result
        assert isinstance(result["flagged"], list)
        assert isinstance(result["summary"], dict)

    def test_invalid_method_raises(self, clean_matrix: Path) -> None:
        with pytest.raises(ValueError, match="Unknown detection method"):
            detect_outliers(clean_matrix, method="banana")


# ---------------------------------------------------------------------------
# Unit tests — outlier detection (IQR)
# ---------------------------------------------------------------------------

class TestDetectOutliersIQR:
    def test_clean_no_flags(self, clean_matrix: Path) -> None:
        result = detect_outliers(clean_matrix, method="iqr")
        assert result["summary"]["flagged_count"] == 0

    def test_contaminated_flags_spike(self, contaminated_iqr_matrix: Path) -> None:
        result = detect_outliers(contaminated_iqr_matrix, method="iqr")
        flagged_ids = {f["sample_id"] for f in result["flagged"]}
        assert "s2" in flagged_ids  # s2 is the spiked sample


# ---------------------------------------------------------------------------
# Unit tests — kitome exclusion in detection
# ---------------------------------------------------------------------------

class TestKitomeExclusion:
    def test_kitome_excludes_taxon(self, kitome_matrix: Path) -> None:
        # Without kitome exclusion, Homo sapiens is present
        result_no_kit = detect_outliers(kitome_matrix, kitome_taxa=None)
        taxa_analysed_no_kit = result_no_kit["summary"]["total_taxa_analysed"]

        # With kitome exclusion
        result_kit = detect_outliers(kitome_matrix, kitome_taxa=[9606])
        taxa_analysed_kit = result_kit["summary"]["total_taxa_analysed"]

        assert taxa_analysed_kit == taxa_analysed_no_kit - 1
        assert result_kit["summary"]["kitome_taxa_excluded"] == 1


# ---------------------------------------------------------------------------
# Unit tests — report generation
# ---------------------------------------------------------------------------

class TestReportGeneration:
    def test_generate_report_files(self, contaminated_matrix: Path, tmp_path: Path) -> None:
        result = detect_outliers(contaminated_matrix, method="mad")
        out = tmp_path / "report_out"
        reports = generate_report(result, out)

        assert "flagged_samples" in reports
        assert "qc_summary" in reports
        assert "quarantine_list" in reports

        assert reports["flagged_samples"].exists()
        assert reports["qc_summary"].exists()
        assert reports["quarantine_list"].exists()

    def test_flagged_samples_tsv(self, contaminated_matrix: Path, tmp_path: Path) -> None:
        result = detect_outliers(contaminated_matrix, method="mad")
        out = tmp_path / "rpt"
        reports = generate_report(result, out)

        with open(reports["flagged_samples"]) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
        assert len(rows) > 0
        assert "sample_id" in rows[0]
        assert "tax_id" in rows[0]

    def test_qc_summary_json(self, contaminated_matrix: Path, tmp_path: Path) -> None:
        result = detect_outliers(contaminated_matrix, method="mad")
        out = tmp_path / "rpt"
        reports = generate_report(result, out)

        with open(reports["qc_summary"]) as fh:
            data = json.load(fh)
        assert "total_samples" in data
        assert "method" in data
        assert "flagged_count" in data

    def test_quarantine_list(self, contaminated_matrix: Path, tmp_path: Path) -> None:
        result = detect_outliers(contaminated_matrix, method="mad")
        out = tmp_path / "rpt"
        reports = generate_report(result, out)

        lines = reports["quarantine_list"].read_text().strip().splitlines()
        assert "sC" in lines

    def test_empty_flagged_produces_files(self, clean_matrix: Path, tmp_path: Path) -> None:
        result = detect_outliers(clean_matrix, method="mad")
        out = tmp_path / "rpt"
        reports = generate_report(result, out)

        assert reports["flagged_samples"].exists()
        # Quarantine list should be empty
        assert reports["quarantine_list"].read_text().strip() == ""


# ---------------------------------------------------------------------------
# Unit tests — CLI
# ---------------------------------------------------------------------------

class TestCLI:
    def test_version(self) -> None:
        from csc.detect.cli import main
        with pytest.raises(SystemExit, match="0"):
            main(["--version"])

    def test_missing_matrix(self, tmp_path: Path) -> None:
        from csc.detect.cli import main
        rc = main([
            str(tmp_path / "nonexistent.tsv"),
            "-o", str(tmp_path / "out"),
        ])
        assert rc == 1

    def test_successful_run(self, contaminated_matrix: Path, tmp_path: Path) -> None:
        from csc.detect.cli import main
        out = tmp_path / "cli_out"
        rc = main([
            str(contaminated_matrix),
            "-o", str(out),
            "--method", "mad",
        ])
        assert rc == 0
        assert (out / "flagged_samples.tsv").exists()
        assert (out / "qc_summary.json").exists()
        assert (out / "quarantine_list.txt").exists()

    def test_iqr_method(self, contaminated_matrix: Path, tmp_path: Path) -> None:
        from csc.detect.cli import main
        out = tmp_path / "cli_iqr"
        rc = main([
            str(contaminated_matrix),
            "-o", str(out),
            "--method", "iqr",
            "--iqr-multiplier", "1.5",
        ])
        assert rc == 0

    def test_kitome_cli_flag(self, kitome_matrix: Path, tmp_path: Path) -> None:
        from csc.detect.cli import main
        out = tmp_path / "cli_kit"
        rc = main([
            str(kitome_matrix),
            "-o", str(out),
            "--kitome-taxa", "9606",
        ])
        assert rc == 0

    def test_no_subtract_background_flag(self, contaminated_matrix: Path, tmp_path: Path) -> None:
        from csc.detect.cli import main
        out = tmp_path / "cli_no_bg"
        rc = main([
            str(contaminated_matrix),
            "-o", str(out),
            "--no-subtract-background",
        ])
        assert rc == 0

    def test_verbose_flag(self, contaminated_matrix: Path, tmp_path: Path) -> None:
        from csc.detect.cli import main
        out = tmp_path / "cli_verbose"
        rc = main([
            str(contaminated_matrix),
            "-o", str(out),
            "-v",
        ])
        assert rc == 0


# ---------------------------------------------------------------------------
# Precision / recall on synthetic datasets
# ---------------------------------------------------------------------------

class TestPrecisionRecall:
    """Simulate a cohort with known contaminated and clean samples,
    then measure precision and recall of the detector."""

    @staticmethod
    def _build_cohort(
        tmp_path: Path,
        n_clean: int = 20,
        n_contaminated: int = 5,
        spike_value: float = 5000.0,
    ) -> tuple[Path, set[str], set[str]]:
        """Build a matrix with clean + contaminated samples.

        Returns (matrix_path, clean_ids, contaminated_ids).
        """
        samples = []
        clean_ids: set[str] = set()
        contaminated_ids: set[str] = set()
        for i in range(n_clean):
            sid = f"clean_{i:03d}"
            samples.append(sid)
            clean_ids.add(sid)
        for i in range(n_contaminated):
            sid = f"contam_{i:03d}"
            samples.append(sid)
            contaminated_ids.add(sid)

        header = ["tax_id", "name"] + samples
        # Taxon 1: background (~100 CPM ± small noise)
        row_bg = ["1279", "Staphylococcus"]
        for sid in samples:
            row_bg.append(str(100 + (hash(sid) % 10)))
        # Taxon 2: E. coli – spiked in contaminated samples
        row_ecoli = ["562", "Escherichia coli"]
        for sid in samples:
            if sid in contaminated_ids:
                row_ecoli.append(str(spike_value))
            else:
                row_ecoli.append(str(50 + (hash(sid) % 5)))

        p = tmp_path / "cohort.tsv"
        _write_matrix(p, header, [row_bg, row_ecoli])
        return p, clean_ids, contaminated_ids

    def test_precision_recall_mad(self, tmp_path: Path) -> None:
        matrix, clean, contam = self._build_cohort(tmp_path)
        result = detect_outliers(matrix, method="mad")
        flagged_ids = set(result["summary"]["flagged_samples"])

        true_positives = flagged_ids & contam

        # Recall: should detect all contaminated samples
        recall = len(true_positives) / len(contam) if contam else 0.0
        assert recall >= 0.8, f"Recall too low: {recall}"

        # Precision: should not flag clean samples
        assert not (flagged_ids & clean), "Clean samples should not be flagged"
        precision = (
            len(true_positives) / len(flagged_ids)
            if flagged_ids
            else 1.0
        )
        assert precision >= 0.8, f"Precision too low: {precision}"

    def test_precision_recall_iqr(self, tmp_path: Path) -> None:
        matrix, clean, contam = self._build_cohort(tmp_path)
        result = detect_outliers(matrix, method="iqr")
        flagged_ids = set(result["summary"]["flagged_samples"])

        true_positives = flagged_ids & contam

        recall = len(true_positives) / len(contam) if contam else 0.0
        assert recall >= 0.8, f"Recall too low: {recall}"

        assert not (flagged_ids & clean), "Clean samples should not be flagged"
        precision = (
            len(true_positives) / len(flagged_ids)
            if flagged_ids
            else 1.0
        )
        assert precision >= 0.8, f"Precision too low: {precision}"

    def test_false_discovery_rate(self, tmp_path: Path) -> None:
        """FDR should be below 0.2 on a well-separated cohort."""
        matrix, clean, contam = self._build_cohort(
            tmp_path, n_clean=50, n_contaminated=5, spike_value=10000.0
        )
        result = detect_outliers(matrix, method="mad")
        flagged_ids = set(result["summary"]["flagged_samples"])

        false_positives = flagged_ids & clean
        fdr = len(false_positives) / len(flagged_ids) if flagged_ids else 0.0
        assert fdr <= 0.2, f"FDR too high: {fdr}"


# ---------------------------------------------------------------------------
# Integration: aggregate → detect pipeline
# ---------------------------------------------------------------------------

class TestPipelineIntegration:
    """End-to-end test: write a matrix, run detection, verify reports."""

    def test_end_to_end(self, tmp_path: Path) -> None:
        header = ["tax_id", "name", "s1", "s2", "s3", "s4", "s5"]
        rows = [
            ["1279", "Staphylococcus", "100", "102", "99", "101", "100"],
            ["562", "Escherichia coli", "50", "51", "49", "50", "8000"],
        ]
        matrix = _write_matrix(tmp_path / "matrix.tsv", header, rows)

        result = detect_outliers(matrix, method="mad")
        reports = generate_report(result, tmp_path / "detect_out")

        # Verify s5 was flagged
        with open(reports["flagged_samples"]) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            flagged_rows = list(reader)

        flagged_sample_ids = {r["sample_id"] for r in flagged_rows}
        assert "s5" in flagged_sample_ids

        # Verify quarantine
        quarantine = (tmp_path / "detect_out" / "quarantine_list.txt").read_text()
        assert "s5" in quarantine

        # Verify QC summary
        with open(reports["qc_summary"]) as fh:
            qc = json.load(fh)
        assert qc["total_samples"] == 5
        assert qc["flagged_count"] >= 1


# ---------------------------------------------------------------------------
# Rank-filter tests for detect CLI
# ---------------------------------------------------------------------------

class TestRankFilterDetect:
    """Tests for --rank-filter in csc-detect CLI."""

    def test_cli_rank_filter_with_rank_matrices(self, tmp_path: Path) -> None:
        """When rank-filtered matrices exist, detect should process them."""
        from csc.detect.cli import main

        # Write a typed CPM main matrix
        header = ["tax_id", "name", "sA", "sB", "sC", "sD", "sE"]
        rows = [
            ["1279", "Staphylococcus", "100", "105", "98", "102", "101"],
            ["562", "Escherichia coli", "50", "48", "5000", "49", "51"],
            ["2", "Bacteria", "800", "810", "795", "805", "798"],
        ]
        matrix_dir = tmp_path / "agg_out"
        matrix_dir.mkdir()
        _write_matrix(matrix_dir / "taxa_matrix_cpm.tsv", header, rows)

        # Write a species-only typed CPM matrix
        species_rows = [
            ["562", "Escherichia coli", "50", "48", "5000", "49", "51"],
        ]
        _write_matrix(matrix_dir / "taxa_matrix_cpm_S.tsv", header, species_rows)

        out = tmp_path / "detect_out"
        rc = main([
            str(matrix_dir / "taxa_matrix_cpm.tsv"),
            "-o", str(out),
            "--rank-filter", "S",
        ])
        assert rc == 0

        # Per-rank outputs should exist
        assert (out / "S" / "flagged_samples.tsv").exists()
        assert (out / "S" / "qc_summary.json").exists()
        assert (out / "S" / "quarantine_list.txt").exists()

        # Primary (unfiltered) outputs should also exist
        assert (out / "flagged_samples.tsv").exists()

    def test_cli_rank_filter_with_typed_matrix_names(self, tmp_path: Path) -> None:
        """Rank discovery should follow typed matrix names (taxa_matrix_cpm/raw)."""
        from csc.detect.cli import main

        header = ["tax_id", "name", "sA", "sB", "sC", "sD", "sE"]
        rows = [
            ["1279", "Staphylococcus", "100", "105", "98", "102", "101"],
            ["562", "Escherichia coli", "50", "48", "5000", "49", "51"],
        ]
        matrix_dir = tmp_path / "agg_out"
        matrix_dir.mkdir()
        _write_matrix(matrix_dir / "taxa_matrix_cpm.tsv", header, rows)
        species_rows = [
            ["562", "Escherichia coli", "50", "48", "5000", "49", "51"],
        ]
        _write_matrix(matrix_dir / "taxa_matrix_cpm_S.tsv", header, species_rows)

        out = tmp_path / "detect_out"
        rc = main([
            str(matrix_dir / "taxa_matrix_cpm.tsv"),
            "-o", str(out),
            "--rank-filter", "S",
        ])
        assert rc == 0
        assert (out / "S" / "qc_summary.json").exists()

    def test_cli_rank_filter_no_rank_matrices(self, tmp_path: Path) -> None:
        """When no rank matrices exist, detection still runs on main matrix."""
        from csc.detect.cli import main

        header = ["tax_id", "name", "sA", "sB", "sC", "sD", "sE"]
        rows = [
            ["562", "Escherichia coli", "50", "48", "5000", "49", "51"],
        ]
        matrix_dir = tmp_path / "agg_out"
        matrix_dir.mkdir()
        _write_matrix(matrix_dir / "taxa_matrix.tsv", header, rows)

        out = tmp_path / "detect_out"
        rc = main([
            str(matrix_dir / "taxa_matrix.tsv"),
            "-o", str(out),
            "--rank-filter", "S", "G",
        ])
        assert rc == 0

        # Main output should exist
        assert (out / "flagged_samples.tsv").exists()
        # No per-rank subdirs since no rank matrices exist
        assert not (out / "S").exists()
        assert not (out / "G").exists()

    def test_aggregate_detect_rank_pipeline(self, tmp_path: Path) -> None:
        """End-to-end: aggregate with rank filter → detect on rank matrices."""
        from csc.aggregate.aggregate import aggregate_reports
        from csc.detect.cli import main as detect_main

        # Create reports with multiple rank levels
        d = tmp_path / "reports"
        d.mkdir()
        for sid in ("s1", "s2", "s3", "s4", "s5"):
            # s5 gets a spike in E. coli (species level)
            ecoli_reads = "5000" if sid == "s5" else "50"
            (d / f"{sid}.kraken2.report.txt").write_text(
                "60.00\t600\t200\tU\t0\tunclassified\n"
                "40.00\t400\t50\tR\t1\troot\n"
                "20.00\t200\t100\tD\t2\tBacteria\n"
                "10.00\t100\t80\tG\t1279\tStaphylococcus\n"
                f"5.00\t50\t{ecoli_reads}\tS\t562\tEscherichia coli\n"
                "5.00\t50\t50\tS\t1280\tStaphylococcus aureus\n"
            )

        agg_out = tmp_path / "agg_out"
        agg_result = aggregate_reports(
            sorted(d.glob("*.kraken2.report.txt")),
            agg_out,
            rank_filter=("S", "G"),
        )

        # Species typed matrices should exist
        assert "S" in agg_result["rank_matrices_raw"]
        assert agg_result["rank_matrices_raw"]["S"].exists()

        # Run detect with rank filter on the CPM matrix
        detect_out = tmp_path / "detect_out"
        rc = detect_main([
            str(agg_result["matrix_cpm_path"]),
            "-o", str(detect_out),
            "--rank-filter", "S", "G",
        ])
        assert rc == 0

        # Species-level detection should flag s5
        assert (detect_out / "S" / "flagged_samples.tsv").exists()
        with open(detect_out / "S" / "qc_summary.json") as fh:
            qc = json.load(fh)
        # Species matrix should only have S-rank taxa
        assert qc["total_taxa_analysed"] == 2  # E. coli + S. aureus
