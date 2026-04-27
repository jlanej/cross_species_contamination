"""Tests for the lineage-aware taxonomy module (csc.aggregate.taxonomy).

Covers taxonomy tree loading, lineage traversal, and domain assignment
for the canonical buckets (Bacteria, Archaea, Fungi, Protists, Viruses,
UniVec_Core, Human, Unclassified), plus integration with aggregate_reports
and the detect load_matrix reader.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path

import pytest

from csc.aggregate.taxonomy import (
    DOMAIN_ARCHAEA,
    DOMAIN_BACTERIA,
    DOMAIN_FUNGI,
    DOMAIN_HUMAN,
    DOMAIN_METAZOA_OTHER,
    DOMAIN_PROTISTS,
    DOMAIN_UNCLASSIFIED,
    DOMAIN_UNIVEC_CORE,
    DOMAIN_VIRIDIPLANTAE,
    DOMAIN_VIRUSES,
    TAXID_ARCHAEA,
    TAXID_BACTERIA,
    TAXID_EUKARYOTA,
    TAXID_FUNGI,
    TAXID_HUMAN,
    TAXID_METAZOA,
    TAXID_ROOT,
    TAXID_UNIVEC_CORE,
    TAXID_VIRIDIPLANTAE,
    TAXID_VIRUSES,
    _assign_single_domain,
    _get_lineage,
    assign_domains,
    load_taxonomy_tree,
)


# ---------------------------------------------------------------------------
# Helpers – write a minimal nodes.dmp
# ---------------------------------------------------------------------------

def _write_nodes_dmp(path: Path, entries: list[tuple[int, int]]) -> Path:
    """Write a minimal ``taxonomy/nodes.dmp`` file.

    Each *entry* is a (child_taxid, parent_taxid) tuple.
    """
    tax_dir = path / "taxonomy"
    tax_dir.mkdir(parents=True, exist_ok=True)
    nodes = tax_dir / "nodes.dmp"
    with open(nodes, "w") as fh:
        for child, parent in entries:
            # nodes.dmp uses "\t|\t" separated fields
            fh.write(f"{child}\t|\t{parent}\t|\tno rank\t|\n")
    return path


# A minimal tree that covers all canonical domains:
#   1 (root)
#   ├── 2 (Bacteria)
#   │   └── 1279 (Staphylococcus)
#   │       └── 1280 (S. aureus)
#   ├── 2157 (Archaea)
#   │   └── 2158 (Euryarchaeota)
#   ├── 2759 (Eukaryota)
#   │   ├── 4751 (Fungi)
#   │   │   └── 4890 (Ascomycota)
#   │   ├── 33208 (Metazoa)
#   │   │   ├── 9606 (Homo sapiens)
#   │   │   └── 10090 (Mus musculus)
#   │   ├── 33090 (Viridiplantae)
#   │   │   └── 3700 (Brassicales)
#   │   └── 5794 (Apicomplexa – a protist)
#   ├── 10239 (Viruses)
#   │   └── 11118 (Coronaviridae)
#   └── 81077 (UniVec Core)
#       └── 81078 (synthetic construct)
_MINI_TREE = [
    (1, 1),
    (2, 1),
    (1279, 2),
    (1280, 1279),
    (562, 2),  # E. coli under Bacteria
    (2157, 1),
    (2158, 2157),
    (2759, 1),
    (4751, 2759),
    (4890, 4751),
    (33208, 2759),
    (9606, 33208),
    (10090, 33208),  # Mus musculus under Metazoa
    (33090, 2759),
    (3700, 33090),
    (5794, 2759),
    (10239, 1),
    (11118, 10239),
    (81077, 1),
    (81078, 81077),
]


@pytest.fixture()
def mini_db(tmp_path: Path) -> Path:
    """A mock Kraken2 DB directory with a minimal nodes.dmp."""
    return _write_nodes_dmp(tmp_path / "db", _MINI_TREE)


# ---------------------------------------------------------------------------
# Tests – load_taxonomy_tree
# ---------------------------------------------------------------------------

class TestLoadTaxonomyTree:
    def test_basic_load(self, mini_db: Path) -> None:
        tree = load_taxonomy_tree(mini_db)
        assert tree[1] == 1  # root points to itself
        assert tree[2] == 1
        assert tree[1279] == 2
        assert tree[1280] == 1279

    def test_missing_nodes_dmp_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="nodes.dmp not found"):
            load_taxonomy_tree(tmp_path)

    def test_tree_size(self, mini_db: Path) -> None:
        tree = load_taxonomy_tree(mini_db)
        assert len(tree) == len(_MINI_TREE)

    def test_fallback_to_db_root(self, tmp_path: Path) -> None:
        """nodes.dmp in db root (not taxonomy/) should be found via fallback."""
        db = tmp_path / "db"
        db.mkdir()
        nodes = db / "nodes.dmp"
        with open(nodes, "w") as fh:
            for child, parent in _MINI_TREE:
                fh.write(f"{child}\t|\t{parent}\t|\tno rank\t|\n")
        tree = load_taxonomy_tree(db)
        assert tree[1] == 1
        assert len(tree) == len(_MINI_TREE)

    def test_fallback_to_db_root_emits_warning(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A warning should be logged when falling back to DB root location."""
        db = tmp_path / "db"
        db.mkdir()
        nodes = db / "nodes.dmp"
        with open(nodes, "w") as fh:
            for child, parent in _MINI_TREE:
                fh.write(f"{child}\t|\t{parent}\t|\tno rank\t|\n")
        with caplog.at_level(logging.WARNING, logger="csc.aggregate.taxonomy"):
            load_taxonomy_tree(db)
        assert any("falling back to DB root" in msg for msg in caplog.messages)


# ---------------------------------------------------------------------------
# Tests – _get_lineage
# ---------------------------------------------------------------------------

class TestGetLineage:
    def test_root_lineage(self, mini_db: Path) -> None:
        tree = load_taxonomy_tree(mini_db)
        lineage = _get_lineage(1, tree)
        assert lineage == [1]

    def test_species_lineage(self, mini_db: Path) -> None:
        tree = load_taxonomy_tree(mini_db)
        lineage = _get_lineage(1280, tree)
        assert lineage == [1280, 1279, 2, 1]

    def test_unknown_taxid(self, mini_db: Path) -> None:
        tree = load_taxonomy_tree(mini_db)
        # Taxid not in tree – should return single-element list
        lineage = _get_lineage(99999, tree)
        assert lineage == [99999]


# ---------------------------------------------------------------------------
# Tests – _assign_single_domain
# ---------------------------------------------------------------------------

class TestAssignSingleDomain:
    @pytest.fixture(autouse=True)
    def _load_tree(self, mini_db: Path) -> None:
        self.tree = load_taxonomy_tree(mini_db)

    def test_bacteria(self) -> None:
        assert _assign_single_domain(2, self.tree) == DOMAIN_BACTERIA
        assert _assign_single_domain(1280, self.tree) == DOMAIN_BACTERIA
        assert _assign_single_domain(562, self.tree) == DOMAIN_BACTERIA

    def test_archaea(self) -> None:
        assert _assign_single_domain(2157, self.tree) == DOMAIN_ARCHAEA
        assert _assign_single_domain(2158, self.tree) == DOMAIN_ARCHAEA

    def test_fungi(self) -> None:
        assert _assign_single_domain(4751, self.tree) == DOMAIN_FUNGI
        assert _assign_single_domain(4890, self.tree) == DOMAIN_FUNGI

    def test_viruses(self) -> None:
        assert _assign_single_domain(10239, self.tree) == DOMAIN_VIRUSES
        assert _assign_single_domain(11118, self.tree) == DOMAIN_VIRUSES

    def test_univec_core(self) -> None:
        assert _assign_single_domain(81077, self.tree) == DOMAIN_UNIVEC_CORE
        assert _assign_single_domain(81078, self.tree) == DOMAIN_UNIVEC_CORE

    def test_human(self) -> None:
        assert _assign_single_domain(9606, self.tree) == DOMAIN_HUMAN

    def test_protist(self) -> None:
        # Apicomplexa – eukaryote not in Metazoa/Fungi/Viridiplantae
        assert _assign_single_domain(5794, self.tree) == DOMAIN_PROTISTS

    def test_eukaryota_metazoa_other(self) -> None:
        # Metazoa (non-human) → Metazoa_other
        assert _assign_single_domain(33208, self.tree) == DOMAIN_METAZOA_OTHER

    def test_mouse_metazoa_other(self) -> None:
        # Mouse (taxid 10090) → Metazoa_other
        assert _assign_single_domain(10090, self.tree) == DOMAIN_METAZOA_OTHER

    def test_eukaryota_viridiplantae(self) -> None:
        assert _assign_single_domain(3700, self.tree) == DOMAIN_VIRIDIPLANTAE

    def test_unclassified_taxid_zero(self) -> None:
        assert _assign_single_domain(0, self.tree) == DOMAIN_UNCLASSIFIED

    def test_root(self) -> None:
        # Root itself doesn't belong to any domain
        assert _assign_single_domain(1, self.tree) == DOMAIN_UNCLASSIFIED

    def test_unknown_taxid(self) -> None:
        assert _assign_single_domain(99999, self.tree) == DOMAIN_UNCLASSIFIED


# ---------------------------------------------------------------------------
# Tests – assign_domains (batch)
# ---------------------------------------------------------------------------

class TestAssignDomains:
    def test_batch_assignment(self, mini_db: Path) -> None:
        tree = load_taxonomy_tree(mini_db)
        taxids = {2, 562, 2157, 4890, 9606, 5794, 11118, 81078, 0, 10090, 3700}
        result = assign_domains(taxids, tree)
        assert result[2] == DOMAIN_BACTERIA
        assert result[562] == DOMAIN_BACTERIA
        assert result[2157] == DOMAIN_ARCHAEA
        assert result[4890] == DOMAIN_FUNGI
        assert result[9606] == DOMAIN_HUMAN
        assert result[5794] == DOMAIN_PROTISTS
        assert result[11118] == DOMAIN_VIRUSES
        assert result[81078] == DOMAIN_UNIVEC_CORE
        assert result[0] == DOMAIN_UNCLASSIFIED
        assert result[10090] == DOMAIN_METAZOA_OTHER
        assert result[3700] == DOMAIN_VIRIDIPLANTAE


# ---------------------------------------------------------------------------
# Integration – aggregate_reports with db_path
# ---------------------------------------------------------------------------

def _write_report(path: Path, content: str) -> Path:
    path.write_text(content)
    return path


class TestAggregateWithDomain:
    """Integration tests for domain-annotated matrix output."""

    def test_domain_column_present(self, mini_db: Path, tmp_path: Path) -> None:
        """When db_path is provided, matrices have a 'domain' column."""
        report = _write_report(
            tmp_path / "s1.kraken2.report.txt",
            "50.00\t50\t30\tS\t562\tEscherichia coli\n"
            "50.00\t50\t20\tS\t1280\tStaphylococcus aureus\n",
        )
        out = tmp_path / "out"
        from csc.aggregate.aggregate import aggregate_reports

        result = aggregate_reports([report], out, db_path=mini_db)

        # Read the raw matrix and check for 'domain' header
        with open(result["matrix_raw_path"]) as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
        assert header[2] == "domain"

        # Read data rows and verify domains
        rows = {}
        with open(result["matrix_raw_path"]) as fh:
            reader = csv.reader(fh, delimiter="\t")
            next(reader)  # skip header
            for cols in reader:
                rows[int(cols[0])] = cols[2]  # domain column

        assert rows[562] == DOMAIN_BACTERIA
        assert rows[1280] == DOMAIN_BACTERIA

    def test_no_domain_without_db_path(self, tmp_path: Path) -> None:
        """When db_path is not provided, no domain column."""
        report = _write_report(
            tmp_path / "s1.kraken2.report.txt",
            "100.00\t100\t100\tS\t562\tEscherichia coli\n",
        )
        out = tmp_path / "out"
        from csc.aggregate.aggregate import aggregate_reports

        result = aggregate_reports([report], out)

        with open(result["matrix_raw_path"]) as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
        assert header[2] != "domain"

    def test_domain_in_cpm_matrix(self, mini_db: Path, tmp_path: Path) -> None:
        """CPM matrices should also have the domain column."""
        report = _write_report(
            tmp_path / "s1.kraken2.report.txt",
            "100.00\t100\t100\tS\t562\tEscherichia coli\n",
        )
        out = tmp_path / "out"
        from csc.aggregate.aggregate import aggregate_reports

        result = aggregate_reports([report], out, db_path=mini_db)

        with open(result["matrix_cpm_path"]) as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
        assert header[2] == "domain"

    def test_domain_in_rank_matrices(self, mini_db: Path, tmp_path: Path) -> None:
        """Per-rank matrices should also have the domain column."""
        report = _write_report(
            tmp_path / "s1.kraken2.report.txt",
            "50.00\t50\t30\tS\t562\tEscherichia coli\n"
            "50.00\t50\t20\tG\t1279\tStaphylococcus\n",
        )
        out = tmp_path / "out"
        from csc.aggregate.aggregate import aggregate_reports

        result = aggregate_reports(
            [report], out, db_path=mini_db, rank_filter=("S",)
        )

        if "S" in result["rank_matrices_raw"]:
            with open(result["rank_matrices_raw"]["S"]) as fh:
                reader = csv.reader(fh, delimiter="\t")
                header = next(reader)
            assert header[2] == "domain"

    def test_metadata_records_domain_annotated(
        self, mini_db: Path, tmp_path: Path
    ) -> None:
        """Metadata JSON should include domain_annotated flag."""
        report = _write_report(
            tmp_path / "s1.kraken2.report.txt",
            "100.00\t100\t100\tS\t562\tEscherichia coli\n",
        )
        out = tmp_path / "out"
        from csc.aggregate.aggregate import aggregate_reports

        result = aggregate_reports([report], out, db_path=mini_db)

        with open(result["metadata_path"]) as fh:
            meta = json.load(fh)
        assert meta["domain_annotated"] is True

    def test_metadata_no_domain_without_db(self, tmp_path: Path) -> None:
        """Metadata should record domain_annotated=False without db_path."""
        report = _write_report(
            tmp_path / "s1.kraken2.report.txt",
            "100.00\t100\t100\tS\t562\tEscherichia coli\n",
        )
        out = tmp_path / "out"
        from csc.aggregate.aggregate import aggregate_reports

        result = aggregate_reports([report], out)

        with open(result["metadata_path"]) as fh:
            meta = json.load(fh)
        assert meta["domain_annotated"] is False

    def test_mixed_domains(self, mini_db: Path, tmp_path: Path) -> None:
        """Taxa from different domains should be correctly labelled."""
        report = _write_report(
            tmp_path / "s1.kraken2.report.txt",
            "20.00\t20\t20\tS\t562\tEscherichia coli\n"
            "20.00\t20\t20\tS\t2158\tEuryarchaeota\n"
            "20.00\t20\t20\tS\t4890\tAscomycota\n"
            "20.00\t20\t20\tS\t11118\tCoronaviridae\n"
            "20.00\t20\t20\tS\t9606\tHomo sapiens\n",
        )
        out = tmp_path / "out"
        from csc.aggregate.aggregate import aggregate_reports

        aggregate_reports([report], out, db_path=mini_db)

        # Read domain assignments from the raw matrix
        domains = {}
        with open(out / "taxa_matrix_raw.tsv") as fh:
            reader = csv.reader(fh, delimiter="\t")
            next(reader)  # skip header
            for cols in reader:
                domains[int(cols[0])] = cols[2]

        assert domains[562] == DOMAIN_BACTERIA
        assert domains[2158] == DOMAIN_ARCHAEA
        assert domains[4890] == DOMAIN_FUNGI
        assert domains[11118] == DOMAIN_VIRUSES
        assert domains[9606] == DOMAIN_HUMAN


# ---------------------------------------------------------------------------
# Integration – detect load_matrix with domain column
# ---------------------------------------------------------------------------

class TestDetectLoadMatrixWithDomain:
    """Verify the detect module can read domain-annotated matrices."""

    def test_load_matrix_with_domain(self, tmp_path: Path) -> None:
        """load_matrix should parse domain column correctly."""
        mat = tmp_path / "matrix.tsv"
        mat.write_text(
            "tax_id\tname\tdomain\tsample1\tsample2\n"
            "562\tEscherichia coli\tBacteria\t100\t200\n"
            "9606\tHomo sapiens\tHuman\t5\t10\n"
        )
        from csc.detect.detect import load_matrix

        sample_ids, rows, tax_names = load_matrix(mat)
        assert sample_ids == ["sample1", "sample2"]
        assert len(rows) == 2
        assert rows[0]["tax_id"] == 562
        assert rows[0]["domain"] == "Bacteria"
        assert rows[0]["sample1"] == 100.0
        assert rows[1]["domain"] == "Human"

    def test_load_matrix_without_domain(self, tmp_path: Path) -> None:
        """load_matrix should still work without domain column."""
        mat = tmp_path / "matrix.tsv"
        mat.write_text(
            "tax_id\tname\tsample1\n"
            "562\tEscherichia coli\t100\n"
        )
        from csc.detect.detect import load_matrix

        sample_ids, rows, tax_names = load_matrix(mat)
        assert sample_ids == ["sample1"]
        assert rows[0]["tax_id"] == 562
        assert "domain" not in rows[0]

    def test_end_to_end_aggregate_detect(
        self, mini_db: Path, tmp_path: Path
    ) -> None:
        """Aggregate with domain → detect load_matrix round-trip."""
        report = _write_report(
            tmp_path / "s1.kraken2.report.txt",
            "50.00\t50\t30\tS\t562\tEscherichia coli\n"
            "50.00\t50\t20\tS\t1280\tStaphylococcus aureus\n",
        )
        out = tmp_path / "out"
        from csc.aggregate.aggregate import aggregate_reports
        from csc.detect.detect import load_matrix

        result = aggregate_reports([report], out, db_path=mini_db)
        sample_ids, rows, tax_names = load_matrix(result["matrix_raw_path"])

        assert len(sample_ids) == 1
        assert all("domain" in r for r in rows)
        # All taxa should be Bacteria
        for r in rows:
            assert r["domain"] == DOMAIN_BACTERIA


# ---------------------------------------------------------------------------
# CLI – db-path flag
# ---------------------------------------------------------------------------

class TestCLIDbPath:
    def test_cli_with_db_path(self, mini_db: Path, tmp_path: Path) -> None:
        report = _write_report(
            tmp_path / "s1.kraken2.report.txt",
            "100.00\t100\t100\tS\t562\tEscherichia coli\n",
        )
        from csc.aggregate.cli import main

        rc = main([
            str(report),
            "-o", str(tmp_path / "cli_out"),
            "--db-path", str(mini_db),
        ])
        assert rc == 0

        with open(tmp_path / "cli_out" / "taxa_matrix_raw.tsv") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
        assert header[2] == "domain"

    def test_cli_without_db_path(self, tmp_path: Path) -> None:
        report = _write_report(
            tmp_path / "s1.kraken2.report.txt",
            "100.00\t100\t100\tS\t562\tEscherichia coli\n",
        )
        from csc.aggregate.cli import main

        rc = main([
            str(report),
            "-o", str(tmp_path / "cli_out"),
        ])
        assert rc == 0

        with open(tmp_path / "cli_out" / "taxa_matrix_raw.tsv") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
        assert header[2] != "domain"

    def test_cli_with_flat_db_path(self, tmp_path: Path) -> None:
        """--db-path with nodes.dmp at DB root (no taxonomy/ subdir) should work.

        This is the k2_NCBI_reference_20251007 layout that triggered issue #64.
        """
        # Build a flat DB: nodes.dmp directly under db/, no taxonomy/ dir.
        flat_db = tmp_path / "flat_db"
        flat_db.mkdir()
        with open(flat_db / "nodes.dmp", "w") as fh:
            for child, parent in _MINI_TREE:
                fh.write(f"{child}\t|\t{parent}\t|\tno rank\t|\n")

        report = _write_report(
            tmp_path / "s1.kraken2.report.txt",
            "100.00\t100\t100\tS\t562\tEscherichia coli\n",
        )
        from csc.aggregate.cli import main

        rc = main([
            str(report),
            "-o", str(tmp_path / "cli_out"),
            "--db-path", str(flat_db),
        ])
        assert rc == 0

        # Read matrix once: verify domain column header and data rows.
        rows = {}
        with open(tmp_path / "cli_out" / "taxa_matrix_raw.tsv") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            assert header[2] == "domain"
            for cols in reader:
                rows[int(cols[0])] = cols[2]
        assert rows[562] == DOMAIN_BACTERIA
