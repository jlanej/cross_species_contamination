"""Lineage-aware taxonomy tree loading and domain assignment.

Loads ``taxonomy/nodes.dmp`` (and optionally ``taxonomy/names.dmp``) from a
Kraken2 database directory and assigns each taxon ID to a canonical domain
by traversing the tree up to the root.

Canonical domain taxon IDs
--------------------------
==========  ======  =========================================
Domain      TaxID   Rule
==========  ======  =========================================
Bacteria    2       taxid 2 and all descendants
Archaea     2157    taxid 2157 and all descendants
Fungi       4751    taxid 4751 and all descendants (under Eukaryota)
Viruses     10239   taxid 10239 and all descendants
UniVec_Core 81077   taxid 81077 and all descendants
Human       9606    taxid 9606 and all descendants
Protists    —       Eukaryota (2759) minus Metazoa (33208),
                    Fungi (4751), and Viridiplantae (33090)
==========  ======  =========================================

AI assistance acknowledgment: This module was developed with AI assistance.
Best practices in the bioinformatics field should always take precedence over
specific implementation details.
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# ── Canonical domain anchor taxids ────────────────────────────────────────────

TAXID_BACTERIA = 2
TAXID_ARCHAEA = 2157
TAXID_FUNGI = 4751
TAXID_EUKARYOTA = 2759
TAXID_METAZOA = 33208
TAXID_VIRIDIPLANTAE = 33090
TAXID_VIRUSES = 10239
TAXID_UNIVEC_CORE = 81077
TAXID_HUMAN = 9606
TAXID_ROOT = 1

# Eukaryotic subtrees excluded from the "Protists" bucket.
_EUKARYOTA_EXCLUSIONS = frozenset({TAXID_METAZOA, TAXID_FUNGI, TAXID_VIRIDIPLANTAE})

# Domain label constants.
DOMAIN_BACTERIA = "Bacteria"
DOMAIN_ARCHAEA = "Archaea"
DOMAIN_FUNGI = "Fungi"
DOMAIN_VIRUSES = "Viruses"
DOMAIN_UNIVEC_CORE = "UniVec_Core"
DOMAIN_HUMAN = "Human"
DOMAIN_PROTISTS = "Protists"
DOMAIN_UNCLASSIFIED = "Unclassified"


# ── Tree loading ──────────────────────────────────────────────────────────────


def load_taxonomy_tree(db_path: str | Path) -> dict[int, int]:
    """Parse ``taxonomy/nodes.dmp`` and return a child→parent mapping.

    Parameters
    ----------
    db_path:
        Path to the Kraken2 database directory containing
        ``taxonomy/nodes.dmp``.

    Returns
    -------
    dict[int, int]
        Mapping from each taxon ID to its parent taxon ID.
        The root node (taxid 1) maps to itself.

    Raises
    ------
    FileNotFoundError
        If ``taxonomy/nodes.dmp`` does not exist under *db_path*.
    """
    nodes_path = Path(db_path) / "taxonomy" / "nodes.dmp"
    if not nodes_path.exists():
        raise FileNotFoundError(
            f"taxonomy/nodes.dmp not found in {db_path}"
        )

    tree: dict[int, int] = {}
    with open(nodes_path) as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if len(parts) < 2:
                continue
            try:
                child = int(parts[0].strip())
                parent = int(parts[1].strip())
            except (ValueError, IndexError):
                continue
            tree[child] = parent

    logger.info(
        "Loaded taxonomy tree with %d nodes from %s", len(tree), nodes_path
    )
    return tree


# ── Lineage helpers ───────────────────────────────────────────────────────────


def _get_lineage(taxid: int, tree: dict[int, int]) -> list[int]:
    """Return the lineage from *taxid* up to (and including) the root.

    The returned list starts with *taxid* and ends with taxid 1 (root).
    If a cycle is detected the traversal stops early.
    """
    lineage: list[int] = []
    seen: set[int] = set()
    current = taxid
    while current not in seen:
        lineage.append(current)
        seen.add(current)
        parent = tree.get(current)
        if parent is None or parent == current:
            break
        current = parent
    return lineage


def _assign_single_domain(taxid: int, tree: dict[int, int]) -> str:
    """Assign a canonical domain label to a single taxon ID.

    Traverses the lineage from *taxid* to root and checks membership
    against the canonical domain anchor taxids.
    """
    if taxid == 0:
        return DOMAIN_UNCLASSIFIED

    lineage = _get_lineage(taxid, tree)
    lineage_set = frozenset(lineage)

    # Check non-eukaryotic domains first (unambiguous).
    if TAXID_BACTERIA in lineage_set:
        return DOMAIN_BACTERIA
    if TAXID_ARCHAEA in lineage_set:
        return DOMAIN_ARCHAEA
    if TAXID_VIRUSES in lineage_set:
        return DOMAIN_VIRUSES
    if TAXID_UNIVEC_CORE in lineage_set:
        return DOMAIN_UNIVEC_CORE

    # Eukaryotic sub-domains: Human > Fungi > Protists (residual).
    if TAXID_EUKARYOTA in lineage_set:
        # Human is under Metazoa; check specifically for 9606 lineage.
        if TAXID_HUMAN in lineage_set:
            return DOMAIN_HUMAN
        if TAXID_FUNGI in lineage_set:
            return DOMAIN_FUNGI
        # Metazoa (non-human) and Viridiplantae are excluded from Protists
        # but don't have their own domain bucket in the CSC scheme.
        if lineage_set & _EUKARYOTA_EXCLUSIONS:
            return DOMAIN_UNCLASSIFIED
        return DOMAIN_PROTISTS

    return DOMAIN_UNCLASSIFIED


def assign_domains(
    taxids: set[int] | list[int],
    tree: dict[int, int],
) -> dict[int, str]:
    """Assign canonical domain labels to a collection of taxon IDs.

    Parameters
    ----------
    taxids:
        Taxon IDs to classify.
    tree:
        Child→parent mapping as returned by :func:`load_taxonomy_tree`.

    Returns
    -------
    dict[int, str]
        Mapping from each taxon ID to its domain label.
    """
    return {tid: _assign_single_domain(tid, tree) for tid in taxids}
