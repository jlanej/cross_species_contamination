"""Tests for the static HTML contamination report module (``csc.report``)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from csc.aggregate.aggregate import aggregate_reports
from csc.report import (
    REPORT_SCHEMA_VERSION,
    generate_html_report,
    load_inputs,
)
from csc.report.cli import main as report_cli
from csc.report.report import (
    DEFAULT_VARIANT_IMPACT_THRESHOLD_PPM,
    _domain_composition,
    _parse_matrix,
    _per_sample_stats,
    _shannon_index,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


_REPORT_A = (
    "50.00\t500\t100\tU\t0\tunclassified\n"
    "50.00\t500\t10\tR\t1\troot\n"
    "40.00\t400\t0\tD\t2\tBacteria\n"
    "30.00\t300\t200\tG\t1279\tStaphylococcus\n"
    "20.00\t200\t200\tS\t1280\tStaphylococcus aureus\n"
    "10.00\t100\t100\tS\t562\tEscherichia coli\n"
    "5.00\t50\t50\tS\t9606\tHomo sapiens\n"
)

# Sample B: more diverse, more contamination
_REPORT_B = (
    "30.00\t300\t50\tU\t0\tunclassified\n"
    "70.00\t700\t20\tR\t1\troot\n"
    "60.00\t600\t0\tD\t2\tBacteria\n"
    "40.00\t400\t400\tS\t1280\tStaphylococcus aureus\n"
    "10.00\t100\t100\tS\t562\tEscherichia coli\n"
    "10.00\t100\t100\tS\t1423\tBacillus subtilis\n"
    "2.00\t20\t20\tS\t9606\tHomo sapiens\n"
)

# Sample C: mostly human (clean)
_REPORT_C = (
    "95.00\t950\t950\tS\t9606\tHomo sapiens\n"
    "5.00\t50\t10\tD\t2\tBacteria\n"
    "5.00\t50\t40\tS\t1280\tStaphylococcus aureus\n"
)


def _write_reads_summary(
    path: Path, sample_id: str, mapped: int, unmapped: int
) -> Path:
    doc = {
        "sample_id": sample_id,
        "total_reads": mapped + unmapped,
        "total_mapped": mapped,
        "total_unmapped": unmapped,
        "input": f"/fake/{sample_id}.bam",
        "extraction_time": "2024-01-01T00:00:00Z",
    }
    path.write_text(json.dumps(doc))
    return path


@pytest.fixture
def aggregate_outputs(tmp_path: Path) -> dict[str, Path]:
    """Produce a realistic aggregate output dir with raw+cpm+abs matrices."""
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
# Core unit tests
# ---------------------------------------------------------------------------


class TestParseMatrix:
    """``_parse_matrix`` must preserve ``NA`` as ``None``."""

    def test_na_preserved_as_none(self, tmp_path: Path) -> None:
        p = tmp_path / "m.tsv"
        p.write_text(
            "tax_id\tname\tsampleA\tsampleB\n"
            "1280\tStaphylococcus aureus\t42.5\tNA\n"
            "562\tEscherichia coli\t0\t10\n"
        )
        m = _parse_matrix(p)
        assert m.sample_ids == ["sampleA", "sampleB"]
        assert m.values[1280]["sampleA"] == pytest.approx(42.5)
        assert m.values[1280]["sampleB"] is None
        assert m.values[562]["sampleB"] == pytest.approx(10.0)

    def test_domain_column_parsed(self, tmp_path: Path) -> None:
        p = tmp_path / "m.tsv"
        p.write_text(
            "tax_id\tname\tdomain\tsampleA\n"
            "1280\tStaphylococcus aureus\tBacteria\t5\n"
        )
        m = _parse_matrix(p)
        assert m.tax_domains[1280] == "Bacteria"


class TestShannonAndComposition:
    def test_shannon_uniform_vs_concentrated(self, tmp_path: Path) -> None:
        p = tmp_path / "m.tsv"
        p.write_text(
            "tax_id\tname\tsA\tsB\n"
            "1\tfoo\t50\t99\n"
            "2\tbar\t50\t1\n"
        )
        m = _parse_matrix(p)
        # Uniform should have higher diversity than concentrated
        assert _shannon_index(m, "sA") > _shannon_index(m, "sB")

    def test_domain_composition_buckets_by_annotation(
        self, tmp_path: Path
    ) -> None:
        p = tmp_path / "m.tsv"
        p.write_text(
            "tax_id\tname\tdomain\tsA\n"
            "1\ta\tBacteria\t10\n"
            "2\tb\tBacteria\t20\n"
            "3\tc\tViruses\t5\n"
        )
        m = _parse_matrix(p)
        comp = _domain_composition(m, "sA")
        assert comp == {"Bacteria": 30.0, "Viruses": 5.0}


# ---------------------------------------------------------------------------
# End-to-end report generation tests
# ---------------------------------------------------------------------------


class TestGenerateReport:
    def test_report_generated_with_abs(self, aggregate_outputs: dict) -> None:
        inputs = load_inputs(aggregate_outputs["aggregate_dir"])
        out = aggregate_outputs["tmp_path"] / "report.html"
        html = generate_html_report(inputs, out, threshold_ppm=100.0)
        assert html.exists()
        manifest_path = html.with_name("report_manifest.json")
        assert manifest_path.exists()

        text = html.read_text(encoding="utf-8")
        # Structural sanity
        assert "<!DOCTYPE html>" in text
        assert "<title>Cross-Species Contamination Summary Report</title>" in text
        # All eight required sections (cohort layout, default).
        assert "1. Executive Summary" in text
        assert "2. Methods" in text
        assert "3. Cohort species landscape" in text
        assert "4. Variant-Calling Impact" in text
        assert "5. Detection summary" in text
        assert "6. Per-sample appendix" in text
        assert "7. Discussion" in text
        assert "8. Methods Transparency Checklist" in text
        # Denominator labelling – CPM denominator on tables/figures
        assert "per million classified reads" in text
        assert "ppm of total sequenced reads" in text
        # The three-reads table must appear in Methods
        assert "Total sequenced reads" in text
        assert "Reads extracted for classification" in text
        assert "Reads classified by Kraken2" in text
        # Diversity indices appear in the per-sample appendix
        assert "Shannon" in text
        assert "Simpson" in text
        # Cohort-layout artefacts: paginated tables + sidecar TSV
        assert 'class="paginated"' in text
        assert (out.with_name("per_sample_summary.tsv")).exists()
        # Manifest schema sanity
        m = json.loads(manifest_path.read_text())
        assert m["schema_version"] == REPORT_SCHEMA_VERSION
        assert m["absolute_burden_enabled"] is True
        assert m["sample_count"] == 3
        assert m["layout"] == "cohort"
        # Cohort-layout manifest keys
        assert "species_summary" in m
        assert "partition_counts" in m
        assert "top_species_by_burden" in m
        assert "top_species_by_prevalence" in m

    def test_variant_impact_flags_expected_samples(
        self, aggregate_outputs: dict
    ) -> None:
        """Samples with heavy non-human burden must be flagged in §4.

        sampleC is dominated by human (small non-human burden); sampleA
        and sampleB both have substantial non-human direct-read counts.
        A threshold between sampleC and sampleA/B should flag exactly
        sampleA and sampleB.
        """
        inputs = load_inputs(aggregate_outputs["aggregate_dir"])
        out = aggregate_outputs["tmp_path"] / "report.html"
        generate_html_report(inputs, out, threshold_ppm=200.0)
        manifest = json.loads(out.with_name("report_manifest.json").read_text())
        flagged = set(manifest["samples_flagged_variant_impact"])
        assert "sampleA" in flagged and "sampleB" in flagged
        assert "sampleC" not in flagged
        assert manifest["variant_impact_threshold_ppm"] == 200.0

        # With a very high threshold, no samples flagged.
        generate_html_report(inputs, out, threshold_ppm=1_000_000.0)
        manifest = json.loads(out.with_name("report_manifest.json").read_text())
        assert manifest["samples_flagged_variant_impact"] == []

        # With a zero threshold, every sample with any non-human burden
        # must be flagged.
        generate_html_report(inputs, out, threshold_ppm=0.0)
        manifest = json.loads(out.with_name("report_manifest.json").read_text())
        assert set(manifest["samples_flagged_variant_impact"]) == {
            "sampleA", "sampleB", "sampleC"
        }

    def test_report_without_abs_disables_variant_impact(
        self, tmp_path: Path
    ) -> None:
        """When no idxstats sidecars were supplied, §4 is disabled."""
        report_dir = tmp_path / "reports"
        report_dir.mkdir()
        (report_dir / "sampleA.kraken2.report.txt").write_text(_REPORT_A)
        (report_dir / "sampleB.kraken2.report.txt").write_text(_REPORT_B)

        out_dir = tmp_path / "aggregate_out"
        aggregate_reports(
            sorted(report_dir.glob("*.kraken2.report.txt")), out_dir
        )

        inputs = load_inputs(out_dir)
        assert inputs.matrix_abs is None

        html_path = tmp_path / "report.html"
        generate_html_report(inputs, html_path)
        text = html_path.read_text()
        # §4 renders a callout instead of a table
        assert "4. Variant-Calling Impact" in text
        assert "absolute-burden matrix" in text.lower()
        # Manifest reflects disabled state
        manifest = json.loads(
            html_path.with_name("report_manifest.json").read_text()
        )
        assert manifest["absolute_burden_enabled"] is False
        assert manifest["samples_flagged_variant_impact"] == []

    def test_per_sample_stats_nonhuman_computation(
        self, aggregate_outputs: dict
    ) -> None:
        """Non-human = classified - Homo sapiens (9606)."""
        inputs = load_inputs(aggregate_outputs["aggregate_dir"])
        stats = _per_sample_stats(inputs)
        by_sample = {s["sample_id"]: s for s in stats}

        # sampleA: direct reads for non-human species = 200 + 100 = 300
        # (9606 has 50 direct → human; classified total per sample excludes
        # unclassified, and includes only taxa present as direct reads.)
        # The exact value depends on which rows the matrix includes, but
        # the non-human percentage must be strictly > 0 and < 100.
        for sid in ("sampleA", "sampleB"):
            assert by_sample[sid]["nonhuman_reads"] > 0
            assert by_sample[sid]["nonhuman_pct_of_classified"] is not None
            assert 0 < by_sample[sid]["nonhuman_pct_of_classified"] < 100

        # sampleC: dominated by human (950 direct), small non-human
        assert (
            by_sample["sampleC"]["nonhuman_pct_of_classified"]
            < by_sample["sampleA"]["nonhuman_pct_of_classified"]
        )

        # Absolute-burden columns populated when idxstats supplied.
        for sid in ("sampleA", "sampleB", "sampleC"):
            assert by_sample[sid]["nonhuman_abs_ppm_total"] is not None

    def test_cpm_and_abs_denominators_differ(
        self, aggregate_outputs: dict
    ) -> None:
        """CPM and absolute columns must not produce the same numbers.

        Because total sequenced reads ≠ classified reads, the two
        denominators must yield different values even for the same
        taxon.  This guards against a refactor that accidentally swaps
        the matrices.
        """
        inputs = load_inputs(aggregate_outputs["aggregate_dir"])
        cpm = inputs.matrix_cpm
        abs_m = inputs.matrix_abs
        assert abs_m is not None
        # For S. aureus in sampleA, cpm must be orders of magnitude
        # larger than abs (a few hundred reads out of ~1M sequenced vs
        # out of a few hundred classified).
        cpm_val = cpm.values[1280]["sampleA"]
        abs_val = abs_m.values[1280]["sampleA"]
        assert cpm_val is not None and abs_val is not None
        assert cpm_val > abs_val


class TestCLI:
    def test_cli_success(
        self, aggregate_outputs: dict, capsys: pytest.CaptureFixture
    ) -> None:
        out = aggregate_outputs["tmp_path"] / "cli_report.html"
        rc = report_cli([
            str(aggregate_outputs["aggregate_dir"]),
            "-o", str(out),
        ])
        assert rc == 0
        assert out.exists()
        assert out.with_name("report_manifest.json").exists()
        captured = capsys.readouterr().out
        assert "report:" in captured

    def test_cli_missing_aggregate_dir(
        self, tmp_path: Path
    ) -> None:
        rc = report_cli([
            str(tmp_path / "does-not-exist"),
            "-o", str(tmp_path / "r.html"),
        ])
        assert rc == 1

    def test_cli_invalid_top_n(
        self, aggregate_outputs: dict
    ) -> None:
        out = aggregate_outputs["tmp_path"] / "r.html"
        rc = report_cli([
            str(aggregate_outputs["aggregate_dir"]),
            "-o", str(out),
            "--top-n", "0",
        ])
        assert rc == 1

    def test_default_threshold_constant(self) -> None:
        # Guard against silent changes to the default.
        assert DEFAULT_VARIANT_IMPACT_THRESHOLD_PPM == 1000.0


# ---------------------------------------------------------------------------
# Absolute-burden detection side-pass tests
# ---------------------------------------------------------------------------


class TestAbsDetectInReport:
    """Reports must surface the abs side-pass results when present."""

    def _run_detect_and_load(
        self, aggregate_outputs: dict, no_abs: bool = False
    ) -> tuple:
        from csc.detect.cli import main as detect_main

        agg_dir = aggregate_outputs["aggregate_dir"]
        detect_out = aggregate_outputs["tmp_path"] / "detect_out"
        argv = [
            str(agg_dir / "taxa_matrix_cpm.tsv"),
            "-o", str(detect_out),
        ]
        if no_abs:
            argv.append("--no-abs-detection")
        rc = detect_main(argv)
        assert rc == 0
        inputs = load_inputs(agg_dir, detect_dir=detect_out)
        return detect_out, inputs

    def test_load_inputs_picks_up_abs_detect(
        self, aggregate_outputs: dict
    ) -> None:
        detect_out, inputs = self._run_detect_and_load(aggregate_outputs)
        assert inputs.detect_summary is not None
        assert inputs.abs_detect_summary is not None
        # qc_summary records the matrix_type for both passes
        assert inputs.detect_summary.get("matrix_type") == "cpm"
        assert inputs.abs_detect_summary["matrix_type"] == "abs"
        assert inputs.abs_detect_dir == detect_out / "abs"

    def test_report_renders_abs_detect_section(
        self, aggregate_outputs: dict
    ) -> None:
        _, inputs = self._run_detect_and_load(aggregate_outputs)
        out = aggregate_outputs["tmp_path"] / "report.html"
        generate_html_report(inputs, out)
        text = out.read_text()
        # §2.5 + §5 should be present and self-explanatory
        assert "Absolute-burden side pass" in text
        assert "5. Detection summary" in text
        assert "absolute-burden" in text.lower()
        # The comparison table headings appear (carried over from legacy
        # detection comparison renderer)
        assert "Flagged by both passes" in text
        assert "Flagged only by absolute-burden pass" in text

        manifest = json.loads(
            out.with_name("report_manifest.json").read_text()
        )
        assert manifest["abs_detection_enabled"] is True
        assert manifest["primary_matrix_type"] == "cpm"
        assert "samples_flagged_abs_detect" in manifest
        assert "samples_flagged_primary_detect" in manifest

    def test_report_no_abs_detect_callout(
        self, aggregate_outputs: dict
    ) -> None:
        """When --no-abs-detection was used, §2.5 explains the gap."""
        _, inputs = self._run_detect_and_load(
            aggregate_outputs, no_abs=True
        )
        assert inputs.abs_detect_summary is None
        out = aggregate_outputs["tmp_path"] / "report.html"
        generate_html_report(inputs, out)
        text = out.read_text()
        # Methods §2.5 still appears with an explanatory callout
        assert "Absolute-burden side pass" in text
        assert "was <strong>not run</strong>" in text
        manifest = json.loads(
            out.with_name("report_manifest.json").read_text()
        )
        assert manifest["abs_detection_enabled"] is False
        assert manifest["samples_flagged_abs_detect"] == []
