"""Tests for the per-read Kraken2 confidence recomputation module.

Covers:
 - kmer-string parsing (paired-end, ambiguous tokens, malformed input)
 - per-read confidence calculation against a small synthetic taxonomy
 - lineage `is_descendant` semantics
 - end-to-end report-record reconstruction with thresholding
 - aggregate_reports dual-tier integration
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from csc.aggregate.aggregate import aggregate_reports
from csc.aggregate.confidence import (
    compute_read_confidence,
    filter_records_by_confidence,
    format_threshold_suffix,
    is_descendant,
    iter_kraken2_output,
    map_outputs_to_samples,
    parse_kmer_string,
    sample_id_from_output,
)


# ── Synthetic taxonomy ───────────────────────────────────────────────────────
# Tree: root(1) -> Bacteria(2) -> Acinetobacter(469) -> 1210098
#                        -> Eukaryota(2759) -> Metazoa(33208) -> Human(9606)
TREE: dict[int, int] = {
    1: 1,
    2: 1,
    469: 2,
    1210098: 469,
    2759: 1,
    33208: 2759,
    9606: 33208,
}


# ── parse_kmer_string ────────────────────────────────────────────────────────


class TestParseKmerString:
    def test_basic(self) -> None:
        toks = parse_kmer_string("0:71 1210098:5 0:40 |:| 9606:3 0:113")
        assert toks == [
            ("0", 71), ("1210098", 5), ("0", 40),
            ("9606", 3), ("0", 113),
        ]

    def test_empty(self) -> None:
        assert parse_kmer_string("") == []
        assert parse_kmer_string("   ") == []

    def test_ambiguous_token_preserved(self) -> None:
        # 'A:5' represents 5 kmers spanning N runs; preserved with label "A"
        toks = parse_kmer_string("A:5 0:10 1210098:3")
        labels = [t[0] for t in toks]
        assert "A" in labels

    def test_malformed_tokens_skipped(self) -> None:
        toks = parse_kmer_string("0:71 garbage 1210098:notanumber 9606:3")
        assert toks == [("0", 71), ("9606", 3)]

    def test_zero_count_skipped(self) -> None:
        toks = parse_kmer_string("0:0 1210098:5")
        assert toks == [("1210098", 5)]


# ── is_descendant ────────────────────────────────────────────────────────────


class TestIsDescendant:
    def test_self(self) -> None:
        assert is_descendant(1210098, 1210098, TREE)

    def test_direct_parent(self) -> None:
        assert is_descendant(1210098, 469, TREE)

    def test_grandparent(self) -> None:
        assert is_descendant(1210098, 2, TREE)

    def test_unrelated(self) -> None:
        assert not is_descendant(9606, 1210098, TREE)
        assert not is_descendant(1210098, 9606, TREE)

    def test_root(self) -> None:
        assert is_descendant(1210098, 1, TREE)

    def test_zero(self) -> None:
        assert not is_descendant(0, 1210098, TREE)
        assert not is_descendant(1210098, 0, TREE)


# ── compute_read_confidence ──────────────────────────────────────────────────


class TestComputeReadConfidence:
    def test_issue_example(self) -> None:
        # From the issue: assigned 1210098, kmers "0:71 1210098:5 0:40 |:| 9606:3 0:113"
        # in-clade(1210098): 1210098:5 -> 5
        # total non-ambiguous: 71 + 5 + 40 + 3 + 113 = 232
        # confidence = 5/232
        c = compute_read_confidence(
            1210098, "0:71 1210098:5 0:40 |:| 9606:3 0:113", TREE,
        )
        assert c == pytest.approx(5 / 232)

    def test_high_confidence_human(self) -> None:
        # Human read with most kmers landing on 9606
        c = compute_read_confidence(9606, "9606:100 0:10 33208:5", TREE)
        # in-clade(9606) includes 9606 (=100) and 33208 (Metazoa) is an
        # *ancestor*, NOT a descendant, so it doesn't count.
        # in-clade = 100; total = 115
        assert c == pytest.approx(100 / 115)

    def test_clade_includes_descendants(self) -> None:
        # Assigned to genus Acinetobacter (469); kmers at species 1210098
        # SHOULD count toward in-clade.
        c = compute_read_confidence(469, "469:5 1210098:5 0:10", TREE)
        assert c == pytest.approx(10 / 20)

    def test_ambiguous_excluded_from_denominator(self) -> None:
        c = compute_read_confidence(1210098, "1210098:5 A:100", TREE)
        assert c == pytest.approx(1.0)  # A:100 is excluded entirely

    def test_unclassified(self) -> None:
        assert compute_read_confidence(0, "0:50", TREE) == 0.0

    def test_no_kmers(self) -> None:
        assert compute_read_confidence(9606, "", TREE) == 0.0
        assert compute_read_confidence(9606, "A:50", TREE) == 0.0


# ── filter_records_by_confidence end-to-end ──────────────────────────────────


def _write_output(path: Path, lines: list[str]) -> Path:
    path.write_text("\n".join(lines) + "\n")
    return path


def test_filter_records_by_confidence(tmp_path: Path) -> None:
    out = _write_output(
        tmp_path / "S1.kraken2.output.txt",
        [
            # Confident human read (100/115 ≈ 0.87)
            "C\tr1\t9606\t150\t9606:100 0:10 33208:5",
            # Low-confidence Acinetobacter read (5/232 ≈ 0.02)
            "C\tr2\t1210098\t150|150\t0:71 1210098:5 0:40 |:| 9606:3 0:113",
            # Unclassified read
            "U\tr3\t0\t150\t0:116",
        ],
    )
    # Threshold 0.0 keeps all classified reads.
    recs = filter_records_by_confidence(
        out,
        threshold=0.0,
        tree=TREE,
        names={9606: "Homo sapiens", 1210098: "Acinetobacter sp."},
        ranks={9606: "S", 1210098: "S"},
    )
    by_tid = {r["tax_id"]: r for r in recs}
    assert by_tid[9606]["direct_reads"] == 1
    assert by_tid[1210098]["direct_reads"] == 1
    assert by_tid[0]["direct_reads"] == 1  # unclassified preserved

    # Threshold 0.5 should drop the low-confidence Acinetobacter read.
    recs = filter_records_by_confidence(
        out, threshold=0.5, tree=TREE,
        names={9606: "Homo sapiens"}, ranks={9606: "S"},
    )
    by_tid = {r["tax_id"]: r for r in recs}
    assert by_tid[9606]["direct_reads"] == 1
    assert 1210098 not in by_tid  # demoted to unclassified
    assert by_tid[0]["direct_reads"] == 2  # original 1 + demoted 1


def test_filter_propagates_clade_reads(tmp_path: Path) -> None:
    out = _write_output(
        tmp_path / "S1.kraken2.output.txt",
        [
            "C\tr1\t1210098\t150\t1210098:100",  # confident species hit
            "C\tr2\t1210098\t150\t1210098:100",
        ],
    )
    recs = filter_records_by_confidence(out, threshold=0.5, tree=TREE)
    by_tid = {r["tax_id"]: r for r in recs}
    # Direct reads at species
    assert by_tid[1210098]["direct_reads"] == 2
    assert by_tid[1210098]["clade_reads"] == 2
    # Ancestor clade reads should propagate up through 469 -> 2 -> 1
    assert by_tid[469]["clade_reads"] == 2
    assert by_tid[469]["direct_reads"] == 0
    assert by_tid[2]["clade_reads"] == 2
    assert by_tid[1]["clade_reads"] == 2


def test_invalid_threshold(tmp_path: Path) -> None:
    out = _write_output(tmp_path / "S1.kraken2.output.txt", ["U\tr\t0\t150\t0:50"])
    with pytest.raises(ValueError):
        filter_records_by_confidence(out, threshold=1.5, tree=TREE)


def test_iter_kraken2_output_skips_malformed(tmp_path: Path) -> None:
    p = _write_output(
        tmp_path / "x.kraken2.output.txt",
        [
            "C\tr1\t9606\t150\t9606:100",
            "garbage_line",                           # too few cols
            "C\tr2\tNOTANINT\t150\t0:50",             # bad taxid
            "U\tr3\t0\t150\t0:50",
        ],
    )
    rows = list(iter_kraken2_output(p))
    assert [r[1] for r in rows] == [9606, 0]


def test_sample_id_from_output(tmp_path: Path) -> None:
    assert sample_id_from_output("/data/HG01121.kraken2.output.txt") == "HG01121"
    assert sample_id_from_output("/data/foo.txt") == "foo"


def test_map_outputs_to_samples(tmp_path: Path) -> None:
    p1 = tmp_path / "A.kraken2.output.txt"
    p2 = tmp_path / "B.kraken2.output.txt"
    p1.write_text("")
    p2.write_text("")
    m = map_outputs_to_samples([p1, p2])
    assert m == {"A": p1, "B": p2}


def test_format_threshold_suffix() -> None:
    assert format_threshold_suffix(0.5) == "conf0p50"
    assert format_threshold_suffix(0.10) == "conf0p10"
    assert format_threshold_suffix(1.0) == "conf1p00"


def test_detect_cli_tier_discovery(tmp_path: Path) -> None:
    """csc-detect's confidence-tier discovery picks up sibling matrices."""
    from csc.detect.cli import (
        _discover_confidence_tier_matrices,
        _rank_matrix_candidates,
    )

    for name in (
        "taxa_matrix_cpm.tsv", "taxa_matrix_cpm_conf0p50.tsv",
        "taxa_matrix_cpm_conf0p10.tsv", "taxa_matrix_raw.tsv",
        "taxa_matrix_cpm_S.tsv", "taxa_matrix_cpm_S_conf0p50.tsv",
    ):
        (tmp_path / name).write_text("")

    base = tmp_path / "taxa_matrix_cpm.tsv"
    tiers = sorted(p.name for p in _discover_confidence_tier_matrices(base))
    assert tiers == [
        "taxa_matrix_cpm_conf0p10.tsv",
        "taxa_matrix_cpm_conf0p50.tsv",
    ]
    # When called on a tier matrix directly, no further tiers are discovered.
    tier = tmp_path / "taxa_matrix_cpm_conf0p50.tsv"
    assert _discover_confidence_tier_matrices(tier) == []
    # Rank candidates resolve to the tier-suffixed rank matrix.
    assert [c.name for c in _rank_matrix_candidates(tier, "S")] == [
        "taxa_matrix_cpm_S_conf0p50.tsv"
    ]


# ── End-to-end aggregate_reports dual-tier integration ───────────────────────


def _write_nodes_dmp(db_root: Path, tree: dict[int, int]) -> None:
    """Write a minimal Kraken2-style nodes.dmp under db_root/taxonomy."""
    tax = db_root / "taxonomy"
    tax.mkdir(parents=True, exist_ok=True)
    with open(tax / "nodes.dmp", "w") as fh:
        for child, parent in tree.items():
            rank = {
                1: "no rank", 2: "domain", 469: "genus",
                1210098: "species", 2759: "domain", 33208: "kingdom",
                9606: "species",
            }.get(child, "no rank")
            fh.write(f"{child}\t|\t{parent}\t|\t{rank}\t|\n")
    with open(tax / "names.dmp", "w") as fh:
        for tid, name in {
            1: "root", 2: "Bacteria", 469: "Acinetobacter",
            1210098: "Acinetobacter sp.", 9606: "Homo sapiens",
            33208: "Metazoa", 2759: "Eukaryota",
        }.items():
            fh.write(f"{tid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")


def _write_report(path: Path, content: str) -> Path:
    path.write_text(content)
    return path


def test_aggregate_reports_with_confidence_tiers(tmp_path: Path) -> None:
    """End-to-end: aggregate_reports produces both sensitive and high-conf tiers."""
    # Build a minimal Kraken2 DB with nodes.dmp
    db = tmp_path / "db"
    _write_nodes_dmp(db, TREE)

    # Write a sensitive-tier report (used for the canonical aggregation)
    reports = tmp_path / "reports"
    reports.mkdir()
    _write_report(
        reports / "S1.kraken2.report.txt",
        "33.33\t1\t1\tU\t0\tunclassified\n"
        "66.67\t2\t0\tR\t1\troot\n"
        "33.33\t1\t0\tD\t2\tBacteria\n"
        "33.33\t1\t1\tS\t1210098\tAcinetobacter sp.\n"
        "33.33\t1\t1\tS\t9606\tHomo sapiens\n",
    )

    # Write the matching kraken2 per-read output
    outs = tmp_path / "outs"
    outs.mkdir()
    _write_output(
        outs / "S1.kraken2.output.txt",
        [
            "C\tr_human\t9606\t150\t9606:100 0:10",       # high-conf human
            "C\tr_bact\t1210098\t150\t0:71 1210098:5 0:40",  # low-conf bact
            "U\tr_un\t0\t150\t0:116",
        ],
    )

    out_dir = tmp_path / "agg"
    result = aggregate_reports(
        [reports / "S1.kraken2.report.txt"],
        out_dir,
        db_path=db,
        confidence_thresholds=[0.5],
        kraken2_output_paths=[outs / "S1.kraken2.output.txt"],
    )

    # Sensitive tier present
    assert result["matrix_raw_path"].exists()
    assert result["matrix_cpm_path"].exists()

    # Confidence tier present
    tiers = result.get("confidence_tiers")
    assert tiers is not None and "conf0p50" in tiers
    tier = tiers["conf0p50"]
    assert tier["threshold"] == pytest.approx(0.5)
    assert tier["matrix_raw_path"].exists()
    # Filename contains the tier suffix
    assert "conf0p50" in tier["matrix_raw_path"].name

    # High-confidence raw matrix should have 1 read for 9606 and 0 for 1210098
    raw = tier["matrix_raw_path"].read_text().strip().splitlines()
    header = raw[0].split("\t")
    sid_idx = header.index("S1")
    rows = {int(line.split("\t")[0]): line.split("\t") for line in raw[1:]}
    assert int(rows[9606][sid_idx]) == 1
    # 1210098 should be 0 (demoted)
    if 1210098 in rows:
        assert int(rows[1210098][sid_idx]) == 0

    # Sensitive raw matrix retains both
    raw_sens = result["matrix_raw_path"].read_text().strip().splitlines()
    header_s = raw_sens[0].split("\t")
    sid_idx_s = header_s.index("S1")
    rows_s = {int(line.split("\t")[0]): line.split("\t") for line in raw_sens[1:]}
    assert int(rows_s[9606][sid_idx_s]) == 1
    assert int(rows_s[1210098][sid_idx_s]) == 1

    # Metadata records the tier
    meta = json.loads(result["metadata_path"].read_text())
    assert "confidence_tiers" in meta
    assert "conf0p50" in meta["confidence_tiers"]
    tier_meta = meta["confidence_tiers"]["conf0p50"]
    assert tier_meta["threshold"] == pytest.approx(0.5)
    # Reads demoted: r_bact is the only confident-tier-failing classified read
    assert tier_meta["reads_demoted_to_unclassified"]["S1"] == 1


def test_aggregate_reports_threshold_zero_no_tier(tmp_path: Path) -> None:
    """Threshold 0.0 alone should not create a confidence tier."""
    db = tmp_path / "db"
    _write_nodes_dmp(db, TREE)
    reports = tmp_path / "reports"
    reports.mkdir()
    _write_report(
        reports / "S1.kraken2.report.txt",
        "100.00\t1\t1\tU\t0\tunclassified\n",
    )
    out_dir = tmp_path / "agg"
    result = aggregate_reports(
        [reports / "S1.kraken2.report.txt"],
        out_dir,
        db_path=db,
        confidence_thresholds=[0.0],
        kraken2_output_paths=None,
    )
    assert not result.get("confidence_tiers")
