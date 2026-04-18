"""Tests for the Kraken2 aggregation module (csc.aggregate).

Covers report parsing, matrix construction, CPM normalisation, edge
cases (missing taxa, empty reports, malformed lines), the CLI entry
point, and a stress test simulating a large cohort.
"""

from __future__ import annotations

import csv
import json
import logging
import math
from pathlib import Path
from unittest import mock

import pytest

from csc.aggregate.aggregate import (
    AggregationResult,
    DEFAULT_RANK_FILTER,
    TaxonRecord,
    VALID_RANK_CODES,
    aggregate_reports,
    parse_kraken2_report,
    sample_id_from_report,
    _collect_sample_counts,
    _collect_sample_clade_counts,
    _compute_pre_filter_total,
)


# ---------------------------------------------------------------------------
# Helpers for generating mock Kraken2 report files
# ---------------------------------------------------------------------------

_BASIC_REPORT = (
    "50.00\t500\t100\tU\t0\tunclassified\n"
    "50.00\t500\t10\tR\t1\troot\n"
    "30.00\t300\t200\tD\t2\tBacteria\n"
    "10.00\t100\t80\tG\t1279\tStaphylococcus\n"
    "5.00\t50\t50\tS\t1280\tStaphylococcus aureus\n"
    "5.00\t50\t50\tS\t562\tEscherichia coli\n"
)


def _write_report(path: Path, content: str = _BASIC_REPORT) -> Path:
    """Write a mock Kraken2 report and return its path."""
    path.write_text(content)
    return path


@pytest.fixture()
def basic_report(tmp_path: Path) -> Path:
    """A well-formed Kraken2 report file."""
    return _write_report(tmp_path / "sample_A.kraken2.report.txt")


@pytest.fixture()
def report_dir(tmp_path: Path) -> Path:
    """Directory containing three distinct report files."""
    d = tmp_path / "reports"
    d.mkdir()

    # Sample A – Staphylococcus-heavy
    _write_report(
        d / "sampleA.kraken2.report.txt",
        "60.00\t600\t200\tU\t0\tunclassified\n"
        "40.00\t400\t50\tR\t1\troot\n"
        "20.00\t200\t150\tD\t2\tBacteria\n"
        "10.00\t100\t80\tG\t1279\tStaphylococcus\n"
        "5.00\t50\t50\tS\t1280\tStaphylococcus aureus\n",
    )

    # Sample B – E. coli heavy
    _write_report(
        d / "sampleB.kraken2.report.txt",
        "40.00\t400\t100\tU\t0\tunclassified\n"
        "60.00\t600\t20\tR\t1\troot\n"
        "50.00\t500\t300\tD\t2\tBacteria\n"
        "30.00\t300\t250\tS\t562\tEscherichia coli\n",
    )

    # Sample C – very few reads
    _write_report(
        d / "sampleC.kraken2.report.txt",
        "90.00\t90\t80\tU\t0\tunclassified\n"
        "10.00\t10\t5\tR\t1\troot\n"
        "5.00\t5\t3\tS\t1280\tStaphylococcus aureus\n",
    )

    return d


# ---------------------------------------------------------------------------
# Unit tests — parsing
# ---------------------------------------------------------------------------

class TestParseKraken2Report:
    def test_basic_parsing(self, basic_report: Path) -> None:
        records = parse_kraken2_report(basic_report)
        assert len(records) == 6
        assert all(isinstance(r, dict) for r in records)

    def test_fields(self, basic_report: Path) -> None:
        records = parse_kraken2_report(basic_report)
        first = records[0]
        assert first["tax_id"] == 0
        assert first["name"] == "unclassified"
        assert first["rank"] == "U"
        assert first["clade_reads"] == 500
        assert first["direct_reads"] == 100
        assert first["percentage"] == pytest.approx(50.0)

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="not found"):
            parse_kraken2_report(tmp_path / "nonexistent.txt")

    def test_empty_file_raises(self, tmp_path: Path) -> None:
        empty = tmp_path / "empty.txt"
        empty.write_text("")
        with pytest.raises(ValueError, match="No valid taxon lines"):
            parse_kraken2_report(empty)

    def test_comment_only_file_raises(self, tmp_path: Path) -> None:
        f = tmp_path / "comments.txt"
        f.write_text("# This is a comment\n# Another comment\n")
        with pytest.raises(ValueError, match="No valid taxon lines"):
            parse_kraken2_report(f)

    def test_malformed_lines_skipped(self, tmp_path: Path) -> None:
        f = tmp_path / "partial.txt"
        f.write_text(
            "not enough\tcolumns\n"
            "50.00\t500\t100\tU\t0\tunclassified\n"
        )
        records = parse_kraken2_report(f)
        assert len(records) == 1

    def test_bad_numeric_skipped(self, tmp_path: Path) -> None:
        f = tmp_path / "bad_num.txt"
        f.write_text(
            "abc\t500\t100\tU\t0\tunclassified\n"
            "50.00\t500\t100\tU\t0\tunclassified\n"
        )
        records = parse_kraken2_report(f)
        assert len(records) == 1


# ---------------------------------------------------------------------------
# Unit tests — sample ID derivation
# ---------------------------------------------------------------------------

class TestSampleIdFromReport:
    def test_standard_suffix(self) -> None:
        assert sample_id_from_report("sample_A.kraken2.report.txt") == "sample_A"

    def test_full_path(self) -> None:
        assert sample_id_from_report("/data/reports/SAMPLE.kraken2.report.txt") == "SAMPLE"

    def test_non_standard_extension(self) -> None:
        # Falls back to Path.stem
        assert sample_id_from_report("report.txt") == "report"

    def test_path_object(self) -> None:
        assert sample_id_from_report(Path("/a/b/S1.kraken2.report.txt")) == "S1"


# ---------------------------------------------------------------------------
# Unit tests — _collect_sample_counts
# ---------------------------------------------------------------------------

class TestCollectSampleCounts:
    def test_basic(self) -> None:
        records: list[TaxonRecord] = [
            TaxonRecord(tax_id=1, name="root", rank="R", clade_reads=100, direct_reads=10, percentage=10.0),
            TaxonRecord(tax_id=2, name="Bacteria", rank="D", clade_reads=90, direct_reads=50, percentage=9.0),
        ]
        counts = _collect_sample_counts(records)
        assert counts == {1: 10, 2: 50}

    def test_min_reads_filter(self) -> None:
        records: list[TaxonRecord] = [
            TaxonRecord(tax_id=1, name="root", rank="R", clade_reads=100, direct_reads=5, percentage=10.0),
            TaxonRecord(tax_id=2, name="Bacteria", rank="D", clade_reads=90, direct_reads=50, percentage=9.0),
        ]
        counts = _collect_sample_counts(records, min_reads=10)
        assert 1 not in counts
        assert counts[2] == 50


# ---------------------------------------------------------------------------
# Unit tests — aggregate_reports
# ---------------------------------------------------------------------------

class TestAggregateReports:
    def test_empty_paths_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="At least one report"):
            aggregate_reports([], tmp_path / "out")

    def test_single_report(self, basic_report: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        result = aggregate_reports([basic_report], out)

        assert result["sample_count"] == 1
        assert result["taxon_count"] == 6
        assert result["matrix_raw_path"].exists()
        assert result["matrix_cpm_path"].exists()
        assert result["metadata_path"].exists()

    def test_multiple_reports(self, report_dir: Path, tmp_path: Path) -> None:
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))
        out = tmp_path / "out"
        result = aggregate_reports(reports, out)

        assert result["sample_count"] == 3
        # All three samples share some taxa + have some unique ones
        assert result["taxon_count"] > 0
        assert result["matrix_raw_path"].exists()
        assert result["matrix_cpm_path"].exists()

    def test_output_dir_created(self, basic_report: Path, tmp_path: Path) -> None:
        out = tmp_path / "nested" / "deep" / "out"
        result = aggregate_reports([basic_report], out)
        assert out.is_dir()
        assert result["matrix_raw_path"].parent == out

    def test_min_reads_filter(self, basic_report: Path, tmp_path: Path) -> None:
        out_all = tmp_path / "all"
        result_all = aggregate_reports([basic_report], out_all, min_reads=0)

        out_filt = tmp_path / "filt"
        result_filt = aggregate_reports([basic_report], out_filt, min_reads=100)

        # Filtering should reduce the number of taxa
        assert result_filt["taxon_count"] <= result_all["taxon_count"]

    def test_raw_counts(self, basic_report: Path, tmp_path: Path) -> None:
        out = tmp_path / "raw"
        result = aggregate_reports([basic_report], out)

        assert result["matrix_raw_path"].exists()
        assert result["matrix_cpm_path"].exists()

        rows = _read_matrix(out / "taxa_matrix_raw.tsv")
        # All values should be integers (no decimal points)
        for row in rows:
            for val in row["values"]:
                assert "." not in val

        # CPM matrix should still be produced
        cpm_rows = _read_matrix(out / "taxa_matrix_cpm.tsv")
        sample0_cpm_sum = sum(float(row["values"][0]) for row in cpm_rows)
        assert sample0_cpm_sum == pytest.approx(1_000_000, rel=1e-4)

    def test_cpm_normalisation(self, basic_report: Path, tmp_path: Path) -> None:
        out = tmp_path / "cpm"
        aggregate_reports([basic_report], out)

        rows = _read_matrix(out / "taxa_matrix_cpm.tsv")
        # Sum of CPM values for the single sample should be ≈ 1,000,000
        total_cpm = sum(float(row["values"][0]) for row in rows)
        assert total_cpm == pytest.approx(1_000_000, rel=1e-4)

    def test_metadata_json(self, basic_report: Path, tmp_path: Path) -> None:
        out = tmp_path / "meta"
        result = aggregate_reports([basic_report], out)

        with open(result["metadata_path"]) as fh:
            meta = json.load(fh)

        assert meta["sample_count"] == 1
        assert "matrix_paths" in meta
        assert meta["matrix_paths"]["raw"].endswith("taxa_matrix_raw.tsv")
        assert meta["matrix_paths"]["cpm"].endswith("taxa_matrix_cpm.tsv")
        assert "primary" not in meta["matrix_paths"]
        assert "normalized" not in meta
        assert isinstance(meta["samples"], list)

    def test_missing_report_skipped(self, basic_report: Path, tmp_path: Path) -> None:
        out = tmp_path / "skip"
        missing = tmp_path / "no_such.txt"
        result = aggregate_reports([basic_report, missing], out)

        # The valid report should still be processed
        assert result["sample_count"] == 1

        with open(result["metadata_path"]) as fh:
            meta = json.load(fh)
        assert len(meta["errors"]) == 1

    def test_matrix_has_all_taxa_across_samples(
        self, report_dir: Path, tmp_path: Path
    ) -> None:
        """Taxa absent in one sample should appear as 0 in the matrix."""
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))
        out = tmp_path / "union"
        aggregate_reports(reports, out)

        rows = _read_matrix(out / "taxa_matrix_raw.tsv")
        # Each row should have a value for every sample (even if 0)
        for row in rows:
            assert len(row["values"]) == 3  # sampleA, sampleB, sampleC

    def test_chunk_size_does_not_affect_result(
        self, report_dir: Path, tmp_path: Path
    ) -> None:
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))

        out1 = tmp_path / "chunk1"
        r1 = aggregate_reports(reports, out1, chunk_size=1)

        out2 = tmp_path / "chunk500"
        r2 = aggregate_reports(reports, out2, chunk_size=500)

        assert r1["sample_count"] == r2["sample_count"]
        assert r1["taxon_count"] == r2["taxon_count"]

        # Matrices should be identical
        m1 = (out1 / "taxa_matrix_raw.tsv").read_text()
        m2 = (out2 / "taxa_matrix_raw.tsv").read_text()
        assert m1 == m2


# ---------------------------------------------------------------------------
# Edge-case tests
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_zero_total_reads_no_division_error(self, tmp_path: Path) -> None:
        """A report where all direct_reads are 0 should not cause ZeroDivisionError."""
        f = tmp_path / "zero.kraken2.report.txt"
        f.write_text(
            "100.00\t100\t0\tU\t0\tunclassified\n"
            "0.00\t0\t0\tR\t1\troot\n"
        )
        out = tmp_path / "out"
        result = aggregate_reports([f], out)
        assert result["sample_count"] == 1

        rows = _read_matrix(out / "taxa_matrix_cpm.tsv")
        for row in rows:
            # All CPM values should be 0.0 (not NaN or Inf)
            for v in row["values"]:
                assert math.isfinite(float(v))

    def test_single_taxon_cpm(self, tmp_path: Path) -> None:
        """A report with only one taxon getting all reads → CPM = 1,000,000."""
        f = tmp_path / "one.kraken2.report.txt"
        f.write_text("100.00\t500\t500\tS\t562\tEscherichia coli\n")
        out = tmp_path / "out"
        aggregate_reports([f], out)

        rows = _read_matrix(out / "taxa_matrix_cpm.tsv")
        assert len(rows) == 1
        assert float(rows[0]["values"][0]) == pytest.approx(1_000_000)

    def test_duplicate_tax_ids_across_reports(self, tmp_path: Path) -> None:
        """Same taxon in multiple reports should merge correctly."""
        d = tmp_path / "dup"
        d.mkdir()
        _write_report(
            d / "s1.kraken2.report.txt",
            "50.00\t50\t30\tS\t562\tEscherichia coli\n"
            "50.00\t50\t20\tS\t1280\tStaphylococcus aureus\n",
        )
        _write_report(
            d / "s2.kraken2.report.txt",
            "80.00\t80\t70\tS\t562\tEscherichia coli\n"
            "20.00\t20\t10\tS\t9606\tHomo sapiens\n",
        )

        out = tmp_path / "out"
        result = aggregate_reports(
            sorted(d.glob("*.kraken2.report.txt")), out
        )

        # s1 has taxa 562 and 1280; s2 has taxa 562 and 9606
        # Union should be 3 taxa
        assert result["taxon_count"] == 3

        rows = _read_matrix(out / "taxa_matrix_raw.tsv")
        by_tid = {int(r["tax_id"]): r for r in rows}

        # s1 column (index 0), s2 column (index 1)
        assert by_tid[562]["values"][0] == "30"  # s1: E. coli
        assert by_tid[562]["values"][1] == "70"  # s2: E. coli
        assert by_tid[1280]["values"][0] == "20"  # s1: S. aureus
        assert by_tid[1280]["values"][1] == "0"  # s2: not present
        assert by_tid[9606]["values"][0] == "0"  # s1: not present
        assert by_tid[9606]["values"][1] == "10"  # s2: H. sapiens


# ---------------------------------------------------------------------------
# Stress test — large cohort simulation
# ---------------------------------------------------------------------------

class TestStress:
    def test_large_cohort_simulation(self, tmp_path: Path) -> None:
        """Simulate 1000 samples to verify efficiency and correctness.

        This is a reduced version of the 100K+ simulation; CI can scale
        up via environment variables if desired.
        """
        import os
        n_samples = int(os.environ.get("CSC_STRESS_SAMPLES", "1000"))
        n_taxa = 50

        d = tmp_path / "reports"
        d.mkdir()

        # Generate n_samples report files
        for i in range(n_samples):
            lines: list[str] = []
            for t in range(n_taxa):
                tid = 1000 + t
                reads = (i * 7 + t * 3) % 100  # deterministic, varied
                lines.append(
                    f"1.00\t{reads}\t{reads}\tS\t{tid}\ttaxon_{tid}\n"
                )
            (d / f"sample_{i:06d}.kraken2.report.txt").write_text(
                "".join(lines)
            )

        out = tmp_path / "out"
        result = aggregate_reports(
            sorted(d.glob("*.kraken2.report.txt")),
            out,
            chunk_size=200,
        )

        assert result["sample_count"] == n_samples
        assert result["taxon_count"] == n_taxa
        assert result["matrix_cpm_path"].exists()

        # Spot-check: CPM per sample should sum to ≈ 1,000,000
        rows = _read_matrix(out / "taxa_matrix_cpm.tsv")
        # Check first and last sample columns
        for col_idx in [0, n_samples - 1]:
            total = sum(float(r["values"][col_idx]) for r in rows)
            if total > 0:
                assert total == pytest.approx(1_000_000, rel=1e-3)


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------

class TestCLI:
    def test_version(self) -> None:
        from csc.aggregate.cli import main

        with pytest.raises(SystemExit) as exc_info:
            main(["--version"])
        assert exc_info.value.code == 0

    def test_missing_required_args(self) -> None:
        from csc.aggregate.cli import main

        with pytest.raises(SystemExit) as exc_info:
            main([])
        assert exc_info.value.code != 0

    def test_missing_input_returns_error(self, tmp_path: Path) -> None:
        from csc.aggregate.cli import main

        rc = main([
            str(tmp_path / "no_such.txt"),
            "-o", str(tmp_path / "out"),
        ])
        assert rc == 1

    def test_successful_cli_run(
        self, basic_report: Path, tmp_path: Path
    ) -> None:
        from csc.aggregate.cli import main

        rc = main([
            str(basic_report),
            "-o", str(tmp_path / "cli_out"),
        ])
        assert rc == 0
        assert not (tmp_path / "cli_out" / "taxa_matrix.tsv").exists()
        assert (tmp_path / "cli_out" / "taxa_matrix_raw.tsv").exists()
        assert (tmp_path / "cli_out" / "taxa_matrix_cpm.tsv").exists()
        assert (tmp_path / "cli_out" / "aggregation_metadata.json").exists()

    def test_cli_passes_params(
        self, basic_report: Path, tmp_path: Path
    ) -> None:
        from csc.aggregate.cli import main

        rc = main([
            str(basic_report),
            "-o", str(tmp_path / "cli_out"),
            "--min-reads", "50",
            "--chunk-size", "100",
        ])
        assert rc == 0

        with open(tmp_path / "cli_out" / "aggregation_metadata.json") as fh:
            meta = json.load(fh)
        assert meta["min_reads"] == 50
        assert "normalized" not in meta

    @mock.patch(
        "csc.aggregate.cli.aggregate_reports",
        side_effect=RuntimeError("test failure"),
    )
    def test_cli_handles_error(
        self, mock_agg: mock.Mock, basic_report: Path, tmp_path: Path
    ) -> None:
        from csc.aggregate.cli import main

        rc = main([
            str(basic_report),
            "-o", str(tmp_path / "out"),
        ])
        assert rc == 1

    def test_verbose_flag(
        self, basic_report: Path, tmp_path: Path
    ) -> None:
        from csc.aggregate.cli import main

        rc = main([
            str(basic_report),
            "-o", str(tmp_path / "out"),
            "--verbose",
        ])
        assert rc == 0


# ---------------------------------------------------------------------------
# Integration: end-to-end pipeline mock
# ---------------------------------------------------------------------------

class TestPipelineIntegration:
    """Simulate the classify → aggregate pipeline using mocked data."""

    def test_classify_output_feeds_aggregate(self, tmp_path: Path) -> None:
        """
        Mock the classify module output and verify that the aggregation
        module correctly consumes the generated report files.
        """
        classify_out = tmp_path / "classify_out"
        classify_out.mkdir()

        # Simulate two ClassificationResult outputs
        for sid in ("SAMPLE_001", "SAMPLE_002"):
            report = classify_out / f"{sid}.kraken2.report.txt"
            _write_report(
                report,
                "60.00\t600\t200\tU\t0\tunclassified\n"
                "40.00\t400\t50\tR\t1\troot\n"
                "20.00\t200\t100\tD\t2\tBacteria\n"
                "10.00\t100\t80\tG\t1279\tStaphylococcus\n"
                f"5.00\t50\t{30 if sid.endswith('1') else 40}\t"
                f"S\t1280\tStaphylococcus aureus\n",
            )

        agg_out = tmp_path / "agg_out"
        reports = sorted(classify_out.glob("*.kraken2.report.txt"))

        result = aggregate_reports(reports, agg_out)

        assert result["sample_count"] == 2
        assert result["matrix_raw_path"].exists()
        assert result["matrix_cpm_path"].exists()

        # Verify both sample IDs are in the matrix header
        with open(result["matrix_raw_path"]) as fh:
            header = fh.readline().strip().split("\t")
        assert "SAMPLE_001" in header
        assert "SAMPLE_002" in header


# ---------------------------------------------------------------------------
# Rank filter tests
# ---------------------------------------------------------------------------

class TestRankFilter:
    """Tests for the rank_filter parameter of aggregate_reports."""

    def test_default_rank_filter(self) -> None:
        assert DEFAULT_RANK_FILTER == ("S", "G", "F")

    def test_valid_rank_codes_contains_defaults(self) -> None:
        for code in DEFAULT_RANK_FILTER:
            assert code in VALID_RANK_CODES

    def test_invalid_rank_code_raises(self, basic_report: Path, tmp_path: Path) -> None:
        """An invalid rank code should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid rank code"):
            aggregate_reports(
                [basic_report], tmp_path / "out", rank_filter=("X",)
            )

    def test_rank_matrices_produced(self, report_dir: Path, tmp_path: Path) -> None:
        """Per-rank matrices are created for ranks present in the data."""
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))
        out = tmp_path / "out"
        result = aggregate_reports(reports, out)

        # The report_dir fixture has S, G, D ranks present
        assert "S" in result["rank_matrices_raw"]
        assert "S" in result["rank_matrices_cpm"]
        assert result["rank_matrices_raw"]["S"].exists()
        assert result["rank_matrices_cpm"]["S"].exists()

        # No legacy compat matrices should exist
        assert not (out / "taxa_matrix.tsv").exists()
        assert not (out / "taxa_matrix_S.tsv").exists()

    def test_rank_matrix_contains_only_matching_rank(
        self, report_dir: Path, tmp_path: Path
    ) -> None:
        """The species-only matrix should contain only S-rank taxa."""
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))
        out = tmp_path / "out"
        result = aggregate_reports(reports, out)

        if "S" not in result["rank_matrices_raw"]:
            pytest.skip("No species-rank taxa in fixture")

        # Read the species matrix
        rows = _read_matrix(result["rank_matrices_raw"]["S"])
        # All rows should correspond to taxa that were S in the reports
        assert len(rows) > 0
        # Verify that the species matrix has fewer rows than the full matrix
        full_rows = _read_matrix(result["matrix_raw_path"])
        assert len(rows) < len(full_rows)

    def test_rank_filter_metadata_sidecar(
        self, report_dir: Path, tmp_path: Path
    ) -> None:
        """The rank_filter_metadata.json sidecar should be written."""
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))
        out = tmp_path / "out"
        result = aggregate_reports(reports, out)

        assert result["rank_metadata_path"].exists()

        with open(result["rank_metadata_path"]) as fh:
            meta = json.load(fh)

        assert "rank_filter" in meta
        assert meta["rank_filter"] == list(DEFAULT_RANK_FILTER)
        assert "ranks" in meta

        # Ranks that have matrices should be present
        for rank, info in meta["ranks"].items():
            assert "matrix_raw_path" in info
            assert "matrix_cpm_path" in info
            assert "taxon_count" in info
            assert "taxa" in info
            assert info["taxon_count"] == len(info["taxa"])

    def test_rank_filter_in_aggregation_metadata(
        self, report_dir: Path, tmp_path: Path
    ) -> None:
        """The aggregation_metadata.json should include rank_filter."""
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))
        out = tmp_path / "out"
        result = aggregate_reports(reports, out)

        with open(result["metadata_path"]) as fh:
            meta = json.load(fh)

        assert "rank_filter" in meta
        assert meta["rank_filter"] == list(DEFAULT_RANK_FILTER)

    def test_custom_rank_filter(self, report_dir: Path, tmp_path: Path) -> None:
        """A custom rank_filter should produce only the requested ranks."""
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))
        out = tmp_path / "out"
        result = aggregate_reports(
            reports, out, rank_filter=("D",)
        )

        # D (domain) should be in rank_matrices_raw if Bacteria is present
        assert "D" in result["rank_matrices_raw"]
        assert "S" not in result["rank_matrices_raw"]
        assert "G" not in result["rank_matrices_raw"]

    def test_empty_rank_filter(self, report_dir: Path, tmp_path: Path) -> None:
        """An empty rank_filter should produce no rank matrices."""
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))
        out = tmp_path / "out"
        result = aggregate_reports(
            reports, out, rank_filter=()
        )

        assert result["rank_matrices_raw"] == {}
        assert result["rank_matrices_cpm"] == {}
        # Typed unfiltered matrices should still exist
        assert result["matrix_raw_path"].exists()
        assert result["matrix_cpm_path"].exists()

    def test_rank_filter_with_no_matching_taxa(
        self, tmp_path: Path
    ) -> None:
        """A rank filter for a rank with no taxa should produce no matrix."""
        f = tmp_path / "only_species.kraken2.report.txt"
        f.write_text("100.00\t500\t500\tS\t562\tEscherichia coli\n")
        out = tmp_path / "out"
        result = aggregate_reports(
            [f], out, rank_filter=("G",)
        )

        # No genus-rank taxa exist, so no G matrix
        assert "G" not in result["rank_matrices_raw"]
        assert "G" not in result["rank_matrices_cpm"]

    def test_rank_matrix_values_match_full_matrix(
        self, report_dir: Path, tmp_path: Path
    ) -> None:
        """Values in species rank-filtered matrix should match the unfiltered matrix.

        For higher-rank matrices (G, F), clade_reads are used instead of
        direct_reads, so they are expected to differ from the unfiltered matrix.
        """
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))
        out = tmp_path / "out"
        result = aggregate_reports(reports, out)

        full_rows = _read_matrix(result["matrix_raw_path"])
        full_by_tid = {r["tax_id"]: r for r in full_rows}

        # Only species-rank matrices should exactly match the full matrix
        if "S" in result["rank_matrices_raw"]:
            rank_rows = _read_matrix(result["rank_matrices_raw"]["S"])
            for rrow in rank_rows:
                tid = rrow["tax_id"]
                assert tid in full_by_tid
                assert rrow["values"] == full_by_tid[tid]["values"]

    def test_cli_rank_filter(self, basic_report: Path, tmp_path: Path) -> None:
        from csc.aggregate.cli import main

        rc = main([
            str(basic_report),
            "-o", str(tmp_path / "cli_out"),
            "--rank-filter", "S",
        ])
        assert rc == 0
        assert not (tmp_path / "cli_out" / "taxa_matrix.tsv").exists()
        assert not (tmp_path / "cli_out" / "taxa_matrix_S.tsv").exists()
        assert (tmp_path / "cli_out" / "taxa_matrix_raw_S.tsv").exists()
        assert (tmp_path / "cli_out" / "taxa_matrix_cpm_S.tsv").exists()
        assert not (tmp_path / "cli_out" / "taxa_matrix_raw_G.tsv").exists()

    def test_cli_rank_filter_multiple(
        self, basic_report: Path, tmp_path: Path
    ) -> None:
        from csc.aggregate.cli import main

        rc = main([
            str(basic_report),
            "-o", str(tmp_path / "cli_out"),
            "--rank-filter", "S", "G",
        ])
        assert rc == 0
        assert (tmp_path / "cli_out" / "taxa_matrix_raw_S.tsv").exists()
        assert (tmp_path / "cli_out" / "taxa_matrix_cpm_S.tsv").exists()
        assert (tmp_path / "cli_out" / "taxa_matrix_raw_G.tsv").exists()
        assert (tmp_path / "cli_out" / "taxa_matrix_cpm_G.tsv").exists()


# ---------------------------------------------------------------------------
# Issue-specific tests
# ---------------------------------------------------------------------------

class TestCPMDenominator:
    """Issue 1: CPM denominator should use pre-filter totals."""

    def test_cpm_uses_pre_filter_denominator(self, tmp_path: Path) -> None:
        """CPM values should be comparable across different min_reads settings.

        Two runs with identical raw data but different min_reads should
        produce CPM values for shared taxa that are proportional to the
        same denominator (all classified reads).
        """
        report = _write_report(
            tmp_path / "sample.kraken2.report.txt",
            "50.00\t50\t30\tS\t562\tEscherichia coli\n"
            "30.00\t30\t20\tS\t1280\tStaphylococcus aureus\n"
            "20.00\t20\t5\tS\t9606\tHomo sapiens\n",
        )

        # Run with min_reads=0 (all taxa pass)
        out_all = tmp_path / "out_all"
        aggregate_reports([report], out_all, min_reads=0)
        rows_all = _read_matrix(out_all / "taxa_matrix_cpm.tsv")
        cpm_ecoli_all = float(
            next(r for r in rows_all if r["tax_id"] == "562")["values"][0]
        )

        # Run with min_reads=10 (only E. coli and S. aureus pass)
        out_filt = tmp_path / "out_filt"
        aggregate_reports([report], out_filt, min_reads=10)
        rows_filt = _read_matrix(out_filt / "taxa_matrix_cpm.tsv")
        cpm_ecoli_filt = float(
            next(r for r in rows_filt if r["tax_id"] == "562")["values"][0]
        )

        # The CPM for E. coli should be the same regardless of min_reads
        # because the denominator uses *all* direct reads (30+20+5=55)
        assert cpm_ecoli_all == pytest.approx(cpm_ecoli_filt, rel=1e-4)

    def test_pre_filter_total_computation(self) -> None:
        records: list[TaxonRecord] = [
            TaxonRecord(tax_id=1, name="a", rank="S", clade_reads=100, direct_reads=50, percentage=50.0),
            TaxonRecord(tax_id=2, name="b", rank="S", clade_reads=60, direct_reads=30, percentage=30.0),
            TaxonRecord(tax_id=3, name="c", rank="S", clade_reads=20, direct_reads=5, percentage=5.0),
        ]
        assert _compute_pre_filter_total(records) == 85  # 50+30+5


class TestCladeReadsForHigherRanks:
    """Issue 2: Genus/family matrices should use clade_reads."""

    def test_genus_matrix_uses_clade_reads(self, tmp_path: Path) -> None:
        """A genus with many species should show clade_reads in the G matrix."""
        report = _write_report(
            tmp_path / "sample.kraken2.report.txt",
            # Genus Staphylococcus: clade=100, direct=5
            # Species S. aureus: clade=50, direct=50
            # Species S. epidermidis: clade=45, direct=45
            "50.00\t100\t5\tG\t1279\tStaphylococcus\n"
            "25.00\t50\t50\tS\t1280\tStaphylococcus aureus\n"
            "22.50\t45\t45\tS\t1282\tStaphylococcus epidermidis\n",
        )
        out = tmp_path / "out"
        result = aggregate_reports([report], out, rank_filter=("G", "S"))

        # Genus matrix should have clade_reads (100), not direct_reads (5)
        g_rows = _read_matrix(result["rank_matrices_raw"]["G"])
        staph = next(r for r in g_rows if r["tax_id"] == "1279")
        assert staph["values"][0] == "100"

        # Species matrix should still use direct_reads
        s_rows = _read_matrix(result["rank_matrices_raw"]["S"])
        aureus = next(r for r in s_rows if r["tax_id"] == "1280")
        assert aureus["values"][0] == "50"

    def test_collect_sample_clade_counts(self) -> None:
        records: list[TaxonRecord] = [
            TaxonRecord(tax_id=1, name="root", rank="R", clade_reads=100, direct_reads=10, percentage=10.0),
            TaxonRecord(tax_id=2, name="Bacteria", rank="D", clade_reads=90, direct_reads=50, percentage=9.0),
        ]
        counts = _collect_sample_clade_counts(records)
        assert counts == {1: 100, 2: 90}

    def test_collect_sample_clade_counts_with_min_reads(self) -> None:
        records: list[TaxonRecord] = [
            TaxonRecord(tax_id=1, name="root", rank="R", clade_reads=8, direct_reads=5, percentage=10.0),
            TaxonRecord(tax_id=2, name="Bacteria", rank="D", clade_reads=90, direct_reads=50, percentage=9.0),
        ]
        counts = _collect_sample_clade_counts(records, min_reads=10)
        # Taxon 1 is excluded because clade_reads (8) < 10
        assert 1 not in counts
        assert counts[2] == 90

    def test_clade_counts_retains_high_clade_low_direct(self) -> None:
        """A genus with low direct_reads but high clade_reads is retained."""
        records: list[TaxonRecord] = [
            TaxonRecord(tax_id=1279, name="Staphylococcus", rank="G",
                        clade_reads=10000, direct_reads=2, percentage=50.0),
        ]
        counts = _collect_sample_clade_counts(records, min_reads=10)
        # clade_reads=10000 >= 10, so taxon is kept despite direct_reads=2
        assert counts[1279] == 10000


class TestDuplicateSampleIds:
    """Issue 6: Duplicate sample IDs should produce a warning."""

    def test_duplicate_sample_ids_warn(self, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
        """Two reports yielding the same sample ID should log a warning."""
        d = tmp_path / "reports"
        d.mkdir()
        sub1 = d / "dir1"
        sub1.mkdir()
        sub2 = d / "dir2"
        sub2.mkdir()

        _write_report(
            sub1 / "sample.kraken2.report.txt",
            "100.00\t100\t100\tS\t562\tEscherichia coli\n",
        )
        _write_report(
            sub2 / "sample.kraken2.report.txt",
            "100.00\t100\t50\tS\t562\tEscherichia coli\n",
        )

        out = tmp_path / "out"
        with caplog.at_level(logging.WARNING, logger="csc.aggregate.aggregate"):
            aggregate_reports(
                [sub1 / "sample.kraken2.report.txt",
                 sub2 / "sample.kraken2.report.txt"],
                out,
            )

        assert any("Duplicate sample ID" in msg for msg in caplog.messages)


class TestTaxNameRankInconsistency:
    """Issue 7: Inconsistent tax_id → name/rank across samples."""

    def test_inconsistent_name_logs_warning(self, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
        d = tmp_path / "reports"
        d.mkdir()

        _write_report(
            d / "s1.kraken2.report.txt",
            "100.00\t100\t100\tS\t562\tEscherichia coli\n",
        )
        _write_report(
            d / "s2.kraken2.report.txt",
            "100.00\t100\t100\tS\t562\tE. coli (renamed)\n",
        )

        out = tmp_path / "out"
        with caplog.at_level(logging.WARNING, logger="csc.aggregate.aggregate"):
            aggregate_reports(sorted(d.glob("*.kraken2.report.txt")), out)

        assert any("Inconsistent name for tax_id 562" in msg for msg in caplog.messages)

    def test_inconsistent_rank_logs_warning(self, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
        d = tmp_path / "reports"
        d.mkdir()

        _write_report(
            d / "s1.kraken2.report.txt",
            "100.00\t100\t100\tS\t562\tEscherichia coli\n",
        )
        _write_report(
            d / "s2.kraken2.report.txt",
            "100.00\t100\t100\tG\t562\tEscherichia coli\n",
        )

        out = tmp_path / "out"
        with caplog.at_level(logging.WARNING, logger="csc.aggregate.aggregate"):
            aggregate_reports(sorted(d.glob("*.kraken2.report.txt")), out)

        assert any("Inconsistent rank for tax_id 562" in msg for msg in caplog.messages)


class TestConfidenceValidation:
    """Issue 10: Bounds-checking on --confidence argument."""

    def test_cli_confidence_too_high(self, tmp_path: Path) -> None:
        from csc.classify.cli import main

        dummy = tmp_path / "reads.fastq.gz"
        dummy.touch()
        rc = main([
            str(dummy),
            "--db", str(tmp_path),
            "-o", str(tmp_path / "out"),
            "--confidence", "1.5",
        ])
        assert rc == 1

    def test_cli_confidence_negative(self, tmp_path: Path) -> None:
        from csc.classify.cli import main

        dummy = tmp_path / "reads.fastq.gz"
        dummy.touch()
        rc = main([
            str(dummy),
            "--db", str(tmp_path),
            "-o", str(tmp_path / "out"),
            "--confidence", "-0.1",
        ])
        assert rc == 1

    def test_cli_confidence_valid_boundary(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Valid boundary confidence values should pass validation."""
        from csc.classify.cli import main

        dummy = tmp_path / "reads.fastq.gz"
        dummy.touch()

        for val in ("0.0", "1.0", "0.5"):
            caplog.clear()
            with caplog.at_level(logging.ERROR):
                main([
                    str(dummy),
                    "--db", str(tmp_path),
                    "-o", str(tmp_path / f"out_{val}"),
                    "--confidence", val,
                ])
            # The call may fail (no kraken2 installed) but the error
            # must NOT be about confidence validation.
            assert not any(
                "Invalid --confidence" in msg for msg in caplog.messages
            ), f"confidence={val} incorrectly rejected"

    def test_api_confidence_validation(self, tmp_path: Path) -> None:
        from csc.classify.classify import classify_reads

        dummy = tmp_path / "reads.fastq.gz"
        dummy.touch()
        with pytest.raises(ValueError, match="confidence must be between"):
            classify_reads(
                [dummy],
                tmp_path / "out",
                db=tmp_path,
                confidence=2.0,
            )


# ---------------------------------------------------------------------------
# Matrix reading helper (for test assertions)
# ---------------------------------------------------------------------------

def _read_matrix(path: Path) -> list[dict[str, Any]]:
    """Read a taxa_matrix.tsv and return rows as dicts.

    Each dict has keys ``tax_id``, ``name``, and ``values`` (a list of
    string cell values, one per sample).
    """
    from typing import Any

    rows: list[dict[str, Any]] = []
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        sample_cols = header[2:]  # skip tax_id, name
        for line in reader:
            rows.append(
                {
                    "tax_id": line[0],
                    "name": line[1],
                    "values": line[2:],
                }
            )
    return rows


# ---------------------------------------------------------------------------
# Absolute-burden matrix (idxstats ingestion)
# ---------------------------------------------------------------------------


def _write_reads_summary(
    path: Path, sample_id: str, total_mapped: int, total_unmapped: int
) -> Path:
    """Write a minimal reads_summary.json sidecar."""
    doc = {
        "schema_version": "1.0",
        "sample_id": sample_id,
        "input": f"/fake/{sample_id}.bam",
        "extraction_time": "2024-01-01T00:00:00+00:00",
        "total_mapped": total_mapped,
        "total_unmapped": total_unmapped,
        "total_reads": total_mapped + total_unmapped,
        "per_chromosome": [
            {"chrom": "chr1", "length": 1000, "mapped": total_mapped,
             "unmapped": 0},
            {"chrom": "*", "length": 0, "mapped": 0,
             "unmapped": total_unmapped},
        ],
    }
    path.write_text(json.dumps(doc))
    return path


class TestAbsoluteBurden:
    """Tests for the absolute-burden matrix produced from idxstats sidecars."""

    def test_idxstats_loader(self, tmp_path: Path) -> None:
        from csc.aggregate.aggregate import load_idxstats_map, load_reads_summary

        s1 = _write_reads_summary(tmp_path / "s1.reads_summary.json", "s1", 90, 10)
        s2 = _write_reads_summary(tmp_path / "s2.reads_summary.json", "s2", 1000, 0)
        m = load_idxstats_map([s1, s2])
        assert set(m.keys()) == {"s1", "s2"}
        assert m["s1"]["total_reads"] == 100

        # Individual loader
        doc = load_reads_summary(s1)
        assert doc["total_mapped"] == 90

    def test_load_reads_summary_missing(self, tmp_path: Path) -> None:
        from csc.aggregate.aggregate import load_reads_summary
        with pytest.raises(FileNotFoundError):
            load_reads_summary(tmp_path / "nope.json")

    def test_load_reads_summary_invalid(self, tmp_path: Path) -> None:
        from csc.aggregate.aggregate import load_reads_summary
        bad = tmp_path / "bad.json"
        bad.write_text(json.dumps({"not": "valid"}))
        with pytest.raises(ValueError, match="Invalid reads_summary"):
            load_reads_summary(bad)

    def test_load_idxstats_map_skips_invalid(self, tmp_path: Path) -> None:
        from csc.aggregate.aggregate import load_idxstats_map
        good = _write_reads_summary(tmp_path / "good.reads_summary.json", "good", 1, 1)
        bad = tmp_path / "bad.reads_summary.json"
        bad.write_text("{}")
        missing = tmp_path / "missing.reads_summary.json"
        m = load_idxstats_map([good, bad, missing])
        assert list(m.keys()) == ["good"]

    def test_abs_matrix_written_when_idxstats_supplied(
        self, report_dir: Path, tmp_path: Path
    ) -> None:
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))
        sidecars = [
            _write_reads_summary(
                tmp_path / f"{sid}.reads_summary.json", sid, mapped, unmapped
            )
            for sid, mapped, unmapped in [
                ("sampleA", 900, 100),
                ("sampleB", 500_000, 500_000),
                ("sampleC", 50, 50),
            ]
        ]

        result = aggregate_reports(
            reports, tmp_path / "out", idxstats_paths=sidecars,
        )
        assert "matrix_abs_path" in result
        assert result["matrix_abs_path"].exists()
        assert "S" in result["rank_matrices_abs"]

        # Metadata exposes schema version, provenance, and the abs matrix path
        with open(result["metadata_path"]) as fh:
            meta = json.load(fh)
        assert meta["absolute_burden_enabled"] is True
        assert "abs" in meta["matrix_paths"]
        assert meta["schema_version"]
        assert meta["sample_provenance"]["sampleA"]["total_reads"] == 1000
        assert meta["samples_without_idxstats"] == []

        # Verify absolute-burden math: sampleA has direct_reads=50 for
        # Staphylococcus aureus (tax_id 1280) and total_reads=1000
        # → 50 / 1000 * 1e6 = 50000 CPM-of-total.
        with open(result["matrix_abs_path"]) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = {int(r["tax_id"]): r for r in reader}
        assert rows[1280]["sampleA"] == f"{50_000.0:.4f}"

    def test_abs_matrix_na_for_missing_sidecar(
        self, report_dir: Path, tmp_path: Path
    ) -> None:
        """Samples lacking a sidecar must appear as NA in the abs matrix."""
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))
        # Only provide a sidecar for sampleA
        sidecars = [
            _write_reads_summary(
                tmp_path / "sampleA.reads_summary.json", "sampleA", 900, 100
            ),
        ]
        result = aggregate_reports(
            reports, tmp_path / "out", idxstats_paths=sidecars,
        )
        with open(result["matrix_abs_path"]) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
        # Non-sampleA columns should be NA
        for r in rows:
            assert r["sampleB"] == "NA"
            assert r["sampleC"] == "NA"
            # sampleA has a numeric value (may be '0.0000')
            float(r["sampleA"])  # must parse
        with open(result["metadata_path"]) as fh:
            meta = json.load(fh)
        assert set(meta["samples_without_idxstats"]) == {"sampleB", "sampleC"}

    def test_abs_matrix_all_mapped(self, tmp_path: Path) -> None:
        """Edge case: sample is entirely mapped (no unmapped reads)."""
        report = _write_report(
            tmp_path / "allmapped.kraken2.report.txt",
            "100.00\t1000\t1000\tU\t0\tunclassified\n"
            "0.00\t0\t0\tR\t1\troot\n",
        )
        sc = _write_reads_summary(
            tmp_path / "allmapped.reads_summary.json", "allmapped", 1000, 0
        )
        result = aggregate_reports(
            [report], tmp_path / "out", idxstats_paths=[sc]
        )
        assert result["matrix_abs_path"].exists()

    def test_abs_matrix_all_unmapped(self, tmp_path: Path) -> None:
        """Edge case: >99% unmapped (common for heavy contamination)."""
        report = _write_report(
            tmp_path / "heavy.kraken2.report.txt",
            "10.00\t10\t10\tU\t0\tunclassified\n"
            "90.00\t990\t990\tS\t562\tEscherichia coli\n",
        )
        sc = _write_reads_summary(
            tmp_path / "heavy.reads_summary.json", "heavy", 5, 995
        )
        result = aggregate_reports(
            [report], tmp_path / "out", idxstats_paths=[sc]
        )
        with open(result["matrix_abs_path"]) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = {int(r["tax_id"]): r for r in reader}
        # 990 E. coli reads / 1000 total sequenced * 1e6 = 990_000
        assert rows[562]["heavy"] == f"{990_000.0:.4f}"

    def test_abs_matrix_zero_total_reads(self, basic_report: Path, tmp_path: Path) -> None:
        """Zero total_reads should not divide-by-zero; produces NA."""
        sc = _write_reads_summary(
            tmp_path / "sample_A.reads_summary.json", "sample_A", 0, 0
        )
        result = aggregate_reports(
            [basic_report], tmp_path / "out", idxstats_paths=[sc]
        )
        with open(result["matrix_abs_path"]) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for r in reader:
                assert r["sample_A"] == "NA"

    def test_abs_matrix_not_written_without_idxstats(
        self, report_dir: Path, tmp_path: Path
    ) -> None:
        """No ``matrix_abs_path`` when no sidecars are supplied."""
        reports = sorted(report_dir.glob("*.kraken2.report.txt"))
        result = aggregate_reports(reports, tmp_path / "out")
        assert "matrix_abs_path" not in result
        assert not (tmp_path / "out" / "taxa_matrix_abs.tsv").exists()
        with open(result["metadata_path"]) as fh:
            meta = json.load(fh)
        assert meta["absolute_burden_enabled"] is False

    def test_cli_idxstats_flag(
        self, report_dir: Path, tmp_path: Path
    ) -> None:
        """``--idxstats`` CLI flag must trigger abs-matrix emission."""
        from csc.aggregate.cli import main

        reports = sorted(report_dir.glob("*.kraken2.report.txt"))
        sc = _write_reads_summary(
            tmp_path / "sampleA.reads_summary.json", "sampleA", 900, 100
        )
        out = tmp_path / "cli_abs"
        rc = main(
            [str(reports[0]), str(reports[1]), str(reports[2]),
             "-o", str(out), "--idxstats", str(sc)]
        )
        assert rc == 0
        assert (out / "taxa_matrix_abs.tsv").exists()


class TestTypedFilenames:
    """Smoke tests for the typed filename helpers supporting 'abs'."""

    def test_abs_unfiltered_name(self) -> None:
        from csc.aggregate.aggregate import typed_matrix_filename
        assert typed_matrix_filename("abs") == "taxa_matrix_abs.tsv"

    def test_abs_rank_name(self) -> None:
        from csc.aggregate.aggregate import typed_rank_matrix_filename
        assert typed_rank_matrix_filename("S", "abs") == "taxa_matrix_abs_S.tsv"

    def test_invalid_matrix_type(self) -> None:
        from csc.aggregate.aggregate import typed_matrix_filename
        with pytest.raises(ValueError):
            typed_matrix_filename("bogus")
