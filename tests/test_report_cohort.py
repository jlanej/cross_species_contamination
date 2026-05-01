"""Unit tests for the cohort-oriented report layout (§3 species centre).

These tests cover the new pure-stdlib statistics in
:mod:`csc.report.cohort` and verify that the rendered report contains
the new sections / sidecar artefacts.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from csc.aggregate.aggregate import aggregate_reports
from csc.report import (
    REPORT_SCHEMA_VERSION,
    generate_html_report,
    load_inputs,
)
from csc.report import cohort as ch
from csc.report.report import _parse_matrix
from csc.report.cli import main as report_cli


_REPORT_A = (
    "50.00\t500\t100\tU\t0\tunclassified\n"
    "50.00\t500\t10\tR\t1\troot\n"
    "40.00\t400\t0\tD\t2\tBacteria\n"
    "30.00\t300\t200\tG\t1279\tStaphylococcus\n"
    "20.00\t200\t200\tS\t1280\tStaphylococcus aureus\n"
    "10.00\t100\t100\tS\t562\tEscherichia coli\n"
    "5.00\t50\t50\tS\t9606\tHomo sapiens\n"
)
_REPORT_B = (
    "30.00\t300\t50\tU\t0\tunclassified\n"
    "70.00\t700\t20\tR\t1\troot\n"
    "60.00\t600\t0\tD\t2\tBacteria\n"
    "40.00\t400\t400\tS\t1280\tStaphylococcus aureus\n"
    "10.00\t100\t100\tS\t562\tEscherichia coli\n"
    "10.00\t100\t100\tS\t1423\tBacillus subtilis\n"
    "2.00\t20\t20\tS\t9606\tHomo sapiens\n"
)
_REPORT_C = (
    "95.00\t950\t950\tS\t9606\tHomo sapiens\n"
    "5.00\t50\t10\tD\t2\tBacteria\n"
    "5.00\t50\t40\tS\t1280\tStaphylococcus aureus\n"
)


def _write_reads_summary(
    path: Path, sample_id: str, mapped: int, unmapped: int
) -> Path:
    path.write_text(json.dumps({
        "sample_id": sample_id,
        "total_reads": mapped + unmapped,
        "total_mapped": mapped,
        "total_unmapped": unmapped,
        "input": f"/fake/{sample_id}.bam",
        "extraction_time": "2024-01-01T00:00:00Z",
    }))
    return path


@pytest.fixture
def aggregate_outputs(tmp_path: Path) -> dict[str, Path]:
    report_dir = tmp_path / "reports"
    report_dir.mkdir()
    (report_dir / "sampleA.kraken2.report.txt").write_text(_REPORT_A)
    (report_dir / "sampleB.kraken2.report.txt").write_text(_REPORT_B)
    (report_dir / "sampleC.kraken2.report.txt").write_text(_REPORT_C)
    sidecars = [
        _write_reads_summary(
            tmp_path / "sampleA.reads_summary.json", "sampleA", 900_000, 100_000
        ),
        _write_reads_summary(
            tmp_path / "sampleB.reads_summary.json", "sampleB", 500_000, 500_000
        ),
        _write_reads_summary(
            tmp_path / "sampleC.reads_summary.json", "sampleC", 999_000, 1_000
        ),
    ]
    out = tmp_path / "aggregate_out"
    aggregate_reports(
        sorted(report_dir.glob("*.kraken2.report.txt")),
        out,
        idxstats_paths=sidecars,
    )
    return {"aggregate_dir": out, "tmp_path": tmp_path}


# ---------------------------------------------------------------------------
# Stats primitives
# ---------------------------------------------------------------------------


class TestPercentilesAndQuantiles:
    def test_percentile_linear_interp(self) -> None:
        vals = [1.0, 2.0, 3.0, 4.0, 5.0]
        assert ch._percentile(vals, 0.0) == 1.0
        assert ch._percentile(vals, 1.0) == 5.0
        assert ch._percentile(vals, 0.5) == 3.0
        # Linear interpolation between 2 and 3
        assert ch._percentile(vals, 0.3) == pytest.approx(2.2)

    def test_cohort_quantiles_drops_none_and_nan(self) -> None:
        q = ch.cohort_quantiles([1.0, None, 2.0, float("nan"), 3.0])
        assert q["n"] == 3
        assert q["median"] == 2.0
        assert q["min"] == 1.0
        assert q["max"] == 3.0


class TestHistogram:
    def test_linear_histogram_assigns_all_values(self) -> None:
        edges, counts = ch.histogram(list(range(10)), n_bins=5)
        assert len(counts) == 5
        assert sum(counts) == 10

    def test_log_histogram_drops_non_positive(self) -> None:
        edges, counts = ch.histogram([0, 1, 10, 100, 1000], n_bins=4, log=True)
        assert sum(counts) == 4  # 0 dropped


# ---------------------------------------------------------------------------
# Species summary rows
# ---------------------------------------------------------------------------


class TestSpeciesSummaryRows:
    def test_excludes_human(self, aggregate_outputs: dict) -> None:
        inputs = load_inputs(aggregate_outputs["aggregate_dir"])
        rows = ch.species_summary_rows(
            inputs.matrix_raw, inputs.matrix_cpm, inputs.matrix_abs,
        )
        assert all(r["tax_id"] != 9606 for r in rows)

    def test_prevalence_and_burden(self, aggregate_outputs: dict) -> None:
        inputs = load_inputs(aggregate_outputs["aggregate_dir"])
        rows = ch.species_summary_rows(
            inputs.matrix_raw, inputs.matrix_cpm, inputs.matrix_abs,
        )
        by_tid = {r["tax_id"]: r for r in rows}
        # S. aureus appears in all 3 samples → prevalence 1.0
        assert by_tid[1280]["prevalence"] == pytest.approx(1.0)
        # E. coli in A and B but not C
        assert by_tid[562]["prevalence"] == pytest.approx(2 / 3)
        # B. subtilis only in B
        assert by_tid[1423]["prevalence"] == pytest.approx(1 / 3)
        # Sorted by burden desc: top row should have largest cohort burden ppm
        burdens = [
            r["cohort_burden_ppm"]
            for r in rows
            if not (isinstance(r["cohort_burden_ppm"], float)
                    and math.isnan(r["cohort_burden_ppm"]))
        ]
        assert burdens == sorted(burdens, reverse=True)

    def test_sparkline_is_non_negative(self, aggregate_outputs: dict) -> None:
        inputs = load_inputs(aggregate_outputs["aggregate_dir"])
        rows = ch.species_summary_rows(
            inputs.matrix_raw, inputs.matrix_cpm, inputs.matrix_abs,
        )
        for r in rows:
            assert all(c >= 0 for c in r["sparkline"])
            assert len(r["sparkline"]) == 12


class TestPrevalencePartition:
    def test_thresholds_partition_correctly(self) -> None:
        rows = [
            {"tax_id": 1, "name": "core", "domain": "Bacteria",
             "prevalence": 0.9, "cohort_burden_ppm": 100.0,
             "cohort_raw_total": 999},
            {"tax_id": 2, "name": "accessory", "domain": "Bacteria",
             "prevalence": 0.3, "cohort_burden_ppm": 50.0,
             "cohort_raw_total": 100},
            {"tax_id": 3, "name": "rare", "domain": "Bacteria",
             "prevalence": 0.05, "cohort_burden_ppm": 1.0,
             "cohort_raw_total": 1},
        ]
        p = ch.prevalence_partition(rows, core_threshold=0.5, rare_threshold=0.1)
        assert p["core_count"] == 1
        assert p["accessory_count"] == 1
        assert p["rare_count"] == 1
        assert p["core_top"][0]["name"] == "core"
        assert p["accessory_top"][0]["name"] == "accessory"
        assert p["rare_top"][0]["name"] == "rare"


class TestRankAbundance:
    def test_descending_by_metric(self) -> None:
        rows = [
            {"name": "a", "tax_id": 1, "domain": "X",
             "cohort_burden_ppm": 100.0, "cohort_raw_total": 100},
            {"name": "b", "tax_id": 2, "domain": "X",
             "cohort_burden_ppm": 1.0, "cohort_raw_total": 10},
            {"name": "c", "tax_id": 3, "domain": "X",
             "cohort_burden_ppm": 50.0, "cohort_raw_total": 50},
        ]
        out = ch.rank_abundance(rows)
        assert [r["name"] for r in out] == ["a", "c", "b"]
        assert [r["rank"] for r in out] == [1, 2, 3]


# ---------------------------------------------------------------------------
# Bray–Curtis, hclust, PCoA
# ---------------------------------------------------------------------------


class TestBrayCurtis:
    def test_identical_samples_distance_zero(self, tmp_path: Path) -> None:
        p = tmp_path / "m.tsv"
        p.write_text(
            "tax_id\tname\tA\tB\n"
            "1\tx\t10\t10\n"
            "2\ty\t5\t5\n"
        )
        m = _parse_matrix(p)
        D = ch.bray_curtis_matrix(m, ["A", "B"])
        assert D[0][1] == pytest.approx(0.0)

    def test_disjoint_samples_distance_one(self, tmp_path: Path) -> None:
        p = tmp_path / "m.tsv"
        p.write_text(
            "tax_id\tname\tA\tB\n"
            "1\tx\t10\t0\n"
            "2\ty\t0\t10\n"
        )
        m = _parse_matrix(p)
        D = ch.bray_curtis_matrix(m, ["A", "B"])
        assert D[0][1] == pytest.approx(1.0)

    def test_symmetric_and_zero_diagonal(self, tmp_path: Path) -> None:
        p = tmp_path / "m.tsv"
        p.write_text(
            "tax_id\tname\tA\tB\tC\n"
            "1\tx\t10\t5\t1\n"
            "2\ty\t1\t5\t10\n"
        )
        m = _parse_matrix(p)
        D = ch.bray_curtis_matrix(m, ["A", "B", "C"])
        for i in range(3):
            assert D[i][i] == 0.0
            for j in range(3):
                assert D[i][j] == pytest.approx(D[j][i])


class TestHClust:
    def test_average_linkage_known_tree(self) -> None:
        # 4 points in 1D at positions 0, 1, 10, 11.
        D = [
            [0.0, 1.0, 10.0, 11.0],
            [1.0, 0.0, 9.0, 10.0],
            [10.0, 9.0, 0.0, 1.0],
            [11.0, 10.0, 1.0, 0.0],
        ]
        out = ch.hclust(D, method="average")
        # Order should keep {0,1} together and {2,3} together.
        order = out["order"]
        assert set(order[:2]) == {0, 1}
        assert set(order[2:]) == {2, 3}
        # First two merges are the closest pairs.
        first = out["merges"][0]
        second = out["merges"][1]
        assert first[2] == pytest.approx(1.0)
        assert second[2] == pytest.approx(1.0)
        # Final merge happens at the UPGMA average distance over all
        # 4 leaf pairs across the two clusters: (10+11+9+10)/4 = 10.
        assert out["merges"][2][2] == pytest.approx(10.0)

    def test_single_linkage_min_distance(self) -> None:
        D = [
            [0.0, 1.0, 10.0],
            [1.0, 0.0, 5.0],
            [10.0, 5.0, 0.0],
        ]
        out = ch.hclust(D, method="single")
        # First merge {0,1} at d=1; merge with {2} at min(d_{02}, d_{12}) = 5.
        assert out["merges"][0][2] == pytest.approx(1.0)
        assert out["merges"][1][2] == pytest.approx(5.0)


class TestPCoA:
    def test_eigenvalue_ordering_and_explained(self) -> None:
        # 4 points placed on a square in 2D → exact 2D embedding.
        coords_true = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]
        n = len(coords_true)
        D = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                D[i][j] = math.dist(coords_true[i], coords_true[j])
        out = ch.pcoa_2d(D)
        assert len(out["coords"]) == n
        # Top eigenvalue must be >= second.
        assert out["eigvals"][0] >= out["eigvals"][1] - 1e-6
        # Trace explained should sum to <= 1 + small tolerance.
        assert sum(out["explained_var"]) <= 1.0 + 1e-9


# ---------------------------------------------------------------------------
# Rendered report assertions
# ---------------------------------------------------------------------------


class TestRenderedCohortReport:
    def test_pagination_and_sidecar_emitted(
        self, aggregate_outputs: dict
    ) -> None:
        inputs = load_inputs(aggregate_outputs["aggregate_dir"])
        out = aggregate_outputs["tmp_path"] / "report.html"
        generate_html_report(inputs, out, threshold_ppm=100.0)
        text = out.read_text(encoding="utf-8")
        # Inline JS is present so tables become interactive offline.
        assert "table.paginated" in text
        assert "querySelectorAll('table.paginated')" in text
        # All cohort sub-sections present.
        for marker in (
            "3.1 Species summary",
            "3.2 Prevalence",
            "3.3 Rank-abundance",
            "3.4 Core",
            "3.5 Cohort-wide distribution figures",
            "3.6 Sample × species heatmap",
            "3.8 Per-species drill-down",
        ):
            assert marker in text, f"missing section marker: {marker}"
        # Per-sample TSV sidecar exists and is non-empty.
        tsv = out.with_name("per_sample_summary.tsv")
        assert tsv.exists()
        lines = tsv.read_text().strip().splitlines()
        assert len(lines) == 1 + 3  # header + 3 samples
        assert lines[0].split("\t")[0] == "sample_id"
        # Per-sample appendix paginated, not a literal full-cohort dump
        # in the body of the HTML when n is large – here it's tiny but
        # the data-page-size attribute must be present.
        assert 'data-page-size="25"' in text

    def test_manifest_has_cohort_keys(self, aggregate_outputs: dict) -> None:
        inputs = load_inputs(aggregate_outputs["aggregate_dir"])
        out = aggregate_outputs["tmp_path"] / "report.html"
        generate_html_report(inputs, out, threshold_ppm=100.0)
        m = json.loads(out.with_name("report_manifest.json").read_text())
        assert m["schema_version"] == "2.0"
        assert m["layout"] == "cohort"
        assert "species_summary" in m
        assert "partition_counts" in m
        assert {"core", "accessory", "rare"} <= set(m["partition_counts"])
        assert m["cluster_method"] == "average"
        assert m["cluster_distance"] == "bray"
        # Top species by burden is present (one of S. aureus / E. coli /
        # B. subtilis / Bacteria / Staphylococcus / root in this fixture).
        assert m["top_species_by_burden"] is not None

    def test_legacy_layout_still_works(
        self, aggregate_outputs: dict
    ) -> None:
        inputs = load_inputs(aggregate_outputs["aggregate_dir"])
        out = aggregate_outputs["tmp_path"] / "legacy_report.html"
        generate_html_report(inputs, out, layout="legacy")
        text = out.read_text(encoding="utf-8")
        # Old §3 Results heading is restored under the legacy flag.
        assert "3. Results" in text
        assert "5. Discussion" in text
        m = json.loads(out.with_name("report_manifest.json").read_text())
        assert m["layout"] == "legacy"


class TestCohortCLI:
    def test_cli_layout_flag(
        self, aggregate_outputs: dict, capsys: pytest.CaptureFixture
    ) -> None:
        out = aggregate_outputs["tmp_path"] / "cli_report.html"
        rc = report_cli([
            str(aggregate_outputs["aggregate_dir"]),
            "-o", str(out),
            "--layout", "legacy",
            "--page-size", "10",
        ])
        assert rc == 0
        m = json.loads(out.with_name("report_manifest.json").read_text())
        assert m["layout"] == "legacy"

    def test_cli_invalid_prevalence(self, aggregate_outputs: dict) -> None:
        out = aggregate_outputs["tmp_path"] / "r.html"
        rc = report_cli([
            str(aggregate_outputs["aggregate_dir"]),
            "-o", str(out),
            "--prevalence-core", "0.1",
            "--prevalence-rare", "0.5",
        ])
        assert rc == 1

    def test_cli_invalid_page_size(self, aggregate_outputs: dict) -> None:
        out = aggregate_outputs["tmp_path"] / "r.html"
        rc = report_cli([
            str(aggregate_outputs["aggregate_dir"]),
            "-o", str(out),
            "--page-size", "0",
        ])
        assert rc == 1


class TestSchemaVersion:
    def test_schema_version_bumped(self) -> None:
        assert REPORT_SCHEMA_VERSION == "2.0"
