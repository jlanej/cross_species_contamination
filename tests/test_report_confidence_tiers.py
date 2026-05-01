"""Integration tests for the report's dual-tier confidence integration.

Covers:
 - tier discovery in :func:`csc.report.report.load_inputs`
 - rendered HTML contains the tier picker, tier sections and the
   §5.1 concordance subsection
 - manifest exposes ``confidence_tiers`` with thresholds + suffixes
 - per-tier per-sample appendix sidecars are written
 - single-tier (legacy) fallback is preserved when no high-confidence
   sibling matrices are present.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from csc.aggregate.aggregate import aggregate_reports
from csc.report import generate_html_report, load_inputs


# Synthetic taxonomy used by the confidence-tier tests in test_confidence.py.
# Tree: 1 → 2 (Bacteria) → 469 → 1210098;  1 → 2759 → 33208 → 9606 (Human)
_TREE: dict[int, int] = {
    1: 1, 2: 1, 469: 2, 1210098: 469, 2759: 1, 33208: 2759, 9606: 33208,
}


def _write_nodes_dmp(db: Path, tree: dict[int, int]) -> None:
    db.mkdir(parents=True, exist_ok=True)
    tax = db / "taxonomy"
    tax.mkdir()
    lines = []
    for child, parent in tree.items():
        # tax_id | parent | rank | ... (NCBI dump uses '\t|\t' field sep)
        lines.append(f"{child}\t|\t{parent}\t|\tno rank\t|\t\t|\n")
    (tax / "nodes.dmp").write_text("".join(lines))


def _make_report(text: str) -> str:
    return text


@pytest.fixture
def dual_tier_aggregate_dir(tmp_path: Path) -> Path:
    """Run csc-aggregate with one high-confidence tier + return the output dir."""
    db = tmp_path / "db"
    _write_nodes_dmp(db, _TREE)

    reports = tmp_path / "reports"
    reports.mkdir()
    # Two samples with similar, but not identical, profiles so the
    # detect outputs (when present) can differ between tiers.
    (reports / "S1.kraken2.report.txt").write_text(
        "20.00\t100\t100\tU\t0\tunclassified\n"
        "80.00\t400\t10\tR\t1\troot\n"
        "60.00\t300\t10\tD\t2\tBacteria\n"
        "50.00\t250\t250\tS\t1210098\tAcinetobacter sp.\n"
        "10.00\t50\t50\tS\t9606\tHomo sapiens\n"
    )
    (reports / "S2.kraken2.report.txt").write_text(
        "30.00\t150\t150\tU\t0\tunclassified\n"
        "70.00\t350\t10\tR\t1\troot\n"
        "40.00\t200\t5\tD\t2\tBacteria\n"
        "30.00\t150\t150\tS\t1210098\tAcinetobacter sp.\n"
        "30.00\t150\t150\tS\t9606\tHomo sapiens\n"
    )

    outs = tmp_path / "outs"
    outs.mkdir()
    # Build per-read output files where roughly half the bacterial
    # reads pass a 0.5-confidence threshold and half do not.  Use the
    # same patterns as test_confidence.py to keep deterministic.
    def _write_kraken_output(path: Path, n_high: int, n_low: int) -> None:
        lines = []
        for i in range(n_high):
            # 100 confident k-mers vs 10 ambiguous
            lines.append(f"C\trh{i}\t9606\t150\t9606:100 0:10")
        for i in range(n_low):
            # 5 confident vs 71 noise → fails 0.5 threshold
            lines.append(f"C\trb{i}\t1210098\t150\t0:71 1210098:5 0:40")
        # Some unclassified reads
        for i in range(20):
            lines.append(f"U\tru{i}\t0\t150\t0:116")
        path.write_text("\n".join(lines) + "\n")

    _write_kraken_output(outs / "S1.kraken2.output.txt", n_high=50, n_low=250)
    _write_kraken_output(outs / "S2.kraken2.output.txt", n_high=150, n_low=150)

    out_dir = tmp_path / "aggregate_out"
    aggregate_reports(
        sorted(reports.glob("*.kraken2.report.txt")),
        out_dir,
        db_path=db,
        confidence_thresholds=[0.5],
        kraken2_output_paths=sorted(outs.glob("*.kraken2.output.txt")),
    )
    return out_dir


# ---------------------------------------------------------------------------
# Tier discovery
# ---------------------------------------------------------------------------


class TestTierDiscovery:
    def test_discover_tier_suffixes(self, dual_tier_aggregate_dir: Path) -> None:
        from csc.report.report import _discover_tier_suffixes

        tiers = _discover_tier_suffixes(dual_tier_aggregate_dir)
        assert len(tiers) == 1
        suffix, threshold = tiers[0]
        assert suffix == "conf0p50"
        assert threshold == pytest.approx(0.5)

    def test_load_inputs_attaches_tier_bundle(
        self, dual_tier_aggregate_dir: Path
    ) -> None:
        inputs = load_inputs(dual_tier_aggregate_dir)
        assert "conf0p50" in inputs.confidence_tiers
        tier = inputs.confidence_tiers["conf0p50"]
        assert tier.tier_suffix == "conf0p50"
        assert tier.tier_threshold == pytest.approx(0.5)
        # Sensitive bundle's matrices have at least as many tax ids as
        # the high-confidence tier (high-confidence demotes some to 0).
        assert (
            len(tier.matrix_raw.sample_ids) == len(inputs.matrix_raw.sample_ids)
        )

    def test_parse_tier_suffix_handles_floats(self) -> None:
        from csc.report.report import _parse_tier_suffix

        assert _parse_tier_suffix("taxa_matrix_raw_conf0p10") == (
            "conf0p10", pytest.approx(0.10)
        )
        assert _parse_tier_suffix("taxa_matrix_cpm_conf0p50") == (
            "conf0p50", pytest.approx(0.50)
        )
        assert _parse_tier_suffix("foo_bar") is None


# ---------------------------------------------------------------------------
# Rendered HTML
# ---------------------------------------------------------------------------


class TestDualTierHtml:
    def test_html_contains_tier_picker_and_two_sections(
        self, dual_tier_aggregate_dir: Path, tmp_path: Path
    ) -> None:
        inputs = load_inputs(dual_tier_aggregate_dir)
        out = tmp_path / "report.html"
        generate_html_report(inputs, out, threshold_ppm=100.0)
        html = out.read_text()

        # Tier picker control + two tier sections.
        assert "csc-tier-select" in html
        assert "data-tier='_sensitive'" in html
        assert "data-tier='conf0p50'" in html
        # Both <section class='csc-tier'> blocks present.
        assert html.count("class='csc-tier'") == 2
        # JS that powers the toggle.
        assert "csc-tier-active" in html

    def test_html_contains_concordance_section(
        self, dual_tier_aggregate_dir: Path, tmp_path: Path
    ) -> None:
        inputs = load_inputs(dual_tier_aggregate_dir)
        out = tmp_path / "report.html"
        generate_html_report(inputs, out, threshold_ppm=100.0)
        html = out.read_text()
        assert "Sensitive vs High-Confidence concordance" in html
        # Citation hint should be present.
        assert "Wood" in html and "Marcelino" in html

    def test_html_contains_methods_subsection(
        self, dual_tier_aggregate_dir: Path, tmp_path: Path
    ) -> None:
        inputs = load_inputs(dual_tier_aggregate_dir)
        out = tmp_path / "report.html"
        generate_html_report(inputs, out, threshold_ppm=100.0)
        html = out.read_text()
        assert "2.6 Sensitive vs high-confidence reporting" in html
        # Confidence formula present
        assert "confidence(taxon)" in html

    def test_manifest_records_tier(
        self, dual_tier_aggregate_dir: Path, tmp_path: Path
    ) -> None:
        inputs = load_inputs(dual_tier_aggregate_dir)
        out = tmp_path / "report.html"
        generate_html_report(inputs, out, threshold_ppm=100.0)
        manifest = json.loads(out.with_name("report_manifest.json").read_text())
        assert manifest["schema_version"] == "2.1"
        tiers = manifest.get("confidence_tiers")
        assert isinstance(tiers, list)
        assert len(tiers) == 1
        assert tiers[0]["tier_suffix"] == "conf0p50"
        assert tiers[0]["threshold"] == pytest.approx(0.5)
        assert tiers[0]["matrix_raw"].endswith("conf0p50.tsv")

    def test_per_tier_appendix_tsv_written(
        self, dual_tier_aggregate_dir: Path, tmp_path: Path
    ) -> None:
        inputs = load_inputs(dual_tier_aggregate_dir)
        out = tmp_path / "report.html"
        generate_html_report(inputs, out, threshold_ppm=100.0)
        # Sensitive sidecar always written
        assert (out.with_name("per_sample_summary.tsv")).exists()
        # High-confidence sidecar uses the tier-suffixed name
        assert (out.with_name("per_sample_summary_conf0p50.tsv")).exists()

    def test_demoted_reads_table_present(
        self, dual_tier_aggregate_dir: Path, tmp_path: Path
    ) -> None:
        inputs = load_inputs(dual_tier_aggregate_dir)
        out = tmp_path / "report.html"
        generate_html_report(inputs, out, threshold_ppm=100.0)
        html = out.read_text()
        # The reads_demoted_to_unclassified table should appear.
        assert "Reads demoted to unclassified" in html


# ---------------------------------------------------------------------------
# Single-tier (legacy) backward compatibility
# ---------------------------------------------------------------------------


class TestSingleTierFallback:
    @pytest.fixture
    def sensitive_only_dir(self, tmp_path: Path) -> Path:
        reports = tmp_path / "r"
        reports.mkdir()
        (reports / "S1.kraken2.report.txt").write_text(
            "20.00\t100\t100\tU\t0\tunclassified\n"
            "80.00\t400\t10\tR\t1\troot\n"
            "60.00\t300\t10\tD\t2\tBacteria\n"
            "50.00\t250\t250\tS\t1210098\tAcinetobacter sp.\n"
        )
        out = tmp_path / "agg"
        aggregate_reports(
            sorted(reports.glob("*.kraken2.report.txt")),
            out,
        )
        return out

    def test_no_tier_picker_when_no_sibling_tier(
        self, sensitive_only_dir: Path, tmp_path: Path
    ) -> None:
        inputs = load_inputs(sensitive_only_dir)
        assert inputs.confidence_tiers == {}
        out = tmp_path / "report.html"
        generate_html_report(inputs, out, threshold_ppm=100.0)
        html = out.read_text()
        # No picker, no extra tier section.
        assert "csc-tier-select" not in html
        assert "Sensitive vs High-Confidence concordance" not in html
        # Sensitive-only manifest still includes the (empty) tier list
        manifest = json.loads(out.with_name("report_manifest.json").read_text())
        assert manifest["confidence_tiers"] == []
