"""Per-read confidence recomputation for Kraken2 classifications.

This module provides utilities to compute Kraken2-style **confidence
scores** from existing per-read output files (``*.kraken2.output.txt``),
*without* re-running Kraken2.  It supports building a "high-confidence"
tier of taxonomic assignments alongside the original "sensitive" tier.

Background
----------
Kraken2 classifies a read by looking up each k-mer's lowest common
ancestor (LCA) in the database and walking the path from each k-mer's
LCA up to the root, summing weighted counts to find the *root-to-leaf*
path with the most weight.  When the ``--confidence T`` flag is passed
to Kraken2, it walks back **down** the chosen path and returns the
deepest taxon at which the fraction of in-clade k-mers meets ``T``.

Definition (matches Kraken2's ``classify.cc``)::

    confidence(taxon) = (k-mers whose taxon ∈ clade rooted at taxon)
                        ----------------------------------------------
                        (total k-mers excluding ambiguous ``A`` runs)

Note that k-mers labelled ``0`` (no match) are *included* in the
denominator because they represent kmers that the database failed to
place in any clade — they are evidence *against* the assigned clade.

For a read assigned taxon ``X`` with confidence ``c``, we treat the
classification as "confidently assigned at ``X``" if ``c >= threshold``.
Reads failing the threshold are reported as **unclassified**
(taxid ``0``) in the high-confidence tier.

This implementation works on the per-read kmer LCA string emitted by
Kraken2 in column 5 of ``--output``::

    0:71 1210098:5 0:40 |:| 9606:3 0:113

Each ``taxid:count`` token contributes ``count`` k-mers labelled with
``taxid``.  The ``|:|`` token separates the two reads of a pair and is
treated as a no-op (the kmers from both mates are summed, matching
Kraken2's paired-end semantics).  The ``A:N`` token (ambiguous run of
``N`` ``N`` characters) is excluded from the denominator.

AI assistance acknowledgment: This module was developed with AI
assistance.  Best practices in the bioinformatics field should always
take precedence over specific implementation details.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Iterator

from csc.aggregate.aggregate import TaxonRecord

logger = logging.getLogger(__name__)


# ── Kraken2 per-read output columns ──────────────────────────────────────────
# Format (tab-separated):
#   0: status ('C' classified / 'U' unclassified)
#   1: read ID
#   2: assigned taxon ID
#   3: read length(s) (e.g. "150" or "150|150")
#   4: kmer LCA string (e.g. "0:71 1210098:5 |:| 9606:3 0:113")
_OUT_COL_STATUS = 0
_OUT_COL_READ_ID = 1
_OUT_COL_TAXID = 2
_OUT_COL_LENGTH = 3
_OUT_COL_KMERS = 4
_OUT_MIN_COLUMNS = 5

# Token separating the two reads of a pair in the kmer string.
PAIR_SEP_TOKEN = "|:|"

# Ambiguous-kmer label (kmers spanning ``N`` characters).
AMBIGUOUS_TOKEN = "A"


# ── Tree helpers ─────────────────────────────────────────────────────────────


def is_descendant(taxid: int, ancestor: int, tree: dict[int, int]) -> bool:
    """Return ``True`` if *taxid* lies in the clade rooted at *ancestor*.

    A taxon is considered to be in its own clade.  Walks parent pointers
    until reaching the root or a cycle.  Unknown taxids (not present in
    *tree*) are treated as descendants of themselves only.
    """
    if taxid == ancestor:
        return True
    if taxid <= 0 or ancestor <= 0:
        return False
    seen: set[int] = set()
    current = taxid
    while current not in seen:
        seen.add(current)
        parent = tree.get(current)
        if parent is None or parent == current:
            return False
        if parent == ancestor:
            return True
        current = parent
    return False


# ── Per-read parsing ─────────────────────────────────────────────────────────


def parse_kmer_string(kmers: str) -> list[tuple[str, int]]:
    """Parse a Kraken2 kmer LCA string into ``(label, count)`` tokens.

    The pair separator ``|:|`` is dropped; ambiguous runs are returned
    with label ``"A"`` so callers can choose to exclude them from the
    denominator.  Returns an empty list for empty / whitespace-only
    input.
    """
    out: list[tuple[str, int]] = []
    if not kmers:
        return out
    for tok in kmers.split():
        if tok == PAIR_SEP_TOKEN:
            continue
        if ":" not in tok:
            # Malformed token; skip silently (Kraken2 itself won't emit these).
            continue
        label, _, count_str = tok.rpartition(":")
        try:
            count = int(count_str)
        except ValueError:
            continue
        if count <= 0:
            continue
        out.append((label, count))
    return out


def compute_read_confidence(
    assigned_taxid: int,
    kmers: str,
    tree: dict[int, int],
) -> float:
    """Compute Kraken2-style confidence for a single read.

    Parameters
    ----------
    assigned_taxid:
        The taxon ID Kraken2 assigned to this read (column 3 of the
        per-read output).  ``0`` means unclassified.
    kmers:
        The kmer LCA string (column 5 of the per-read output).
    tree:
        Child→parent mapping from :func:`csc.aggregate.taxonomy.load_taxonomy_tree`.

    Returns
    -------
    float
        Confidence in the range ``[0.0, 1.0]``.  Returns ``0.0`` when
        the read is unclassified or has no non-ambiguous k-mers.
    """
    if assigned_taxid <= 0:
        return 0.0
    in_clade = 0
    total = 0
    for label, count in parse_kmer_string(kmers):
        if label == AMBIGUOUS_TOKEN:
            # Excluded from the denominator (Kraken2 semantics).
            continue
        total += count
        try:
            label_taxid = int(label)
        except ValueError:
            continue
        if is_descendant(label_taxid, assigned_taxid, tree):
            in_clade += count
    if total == 0:
        return 0.0
    return in_clade / total


# ── Per-sample iteration ─────────────────────────────────────────────────────


def iter_kraken2_output(path: str | Path) -> Iterator[tuple[str, int, str]]:
    """Yield ``(status, assigned_taxid, kmer_string)`` for each read.

    Malformed lines are logged at WARNING level and skipped so a single
    corrupt line does not abort the whole sample.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Kraken2 output file not found: {p}")
    with open(p) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n\r")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < _OUT_MIN_COLUMNS:
                logger.warning(
                    "%s:%d: skipping malformed line (%d columns, need %d)",
                    p, lineno, len(cols), _OUT_MIN_COLUMNS,
                )
                continue
            try:
                taxid = int(cols[_OUT_COL_TAXID].strip())
            except ValueError:
                logger.warning(
                    "%s:%d: unparseable taxid %r", p, lineno, cols[_OUT_COL_TAXID],
                )
                continue
            yield cols[_OUT_COL_STATUS].strip(), taxid, cols[_OUT_COL_KMERS]


# ── Filtered TaxonRecord construction ────────────────────────────────────────


# Kraken2 rank codes recognised in standard taxonomy dumps (subset that
# we map directly to single-letter rank codes used elsewhere in csc).
_KRAKEN_RANK_CODES = {
    "no rank": "-",
    "superkingdom": "D",
    "domain": "D",
    "kingdom": "K",
    "phylum": "P",
    "class": "C",
    "order": "O",
    "family": "F",
    "genus": "G",
    "species": "S",
}


def filter_records_by_confidence(
    output_path: str | Path,
    *,
    threshold: float,
    tree: dict[int, int],
    names: dict[int, str] | None = None,
    ranks: dict[int, str] | None = None,
    stats: dict[str, int] | None = None,
) -> list[TaxonRecord]:
    """Build :class:`TaxonRecord` list from per-read output, post-filter.

    Reads whose recomputed confidence is below *threshold* are dropped
    (treated as unclassified).  The remaining reads are tallied per
    taxon ID (``direct_reads``) and propagated up the lineage to give
    ``clade_reads`` for every ancestor.

    The returned records have the same shape as those produced by
    :func:`csc.aggregate.aggregate.parse_kraken2_report` and can be fed
    into the matrix-construction pipeline.

    Parameters
    ----------
    output_path:
        Path to a ``*.kraken2.output.txt`` file.
    threshold:
        Confidence cutoff in ``[0.0, 1.0]``.  Reads scoring strictly
        below this value are excluded.  ``0.0`` retains all classified
        reads (this matches Kraken2's default behaviour).
    tree:
        Child→parent mapping from :func:`csc.aggregate.taxonomy.load_taxonomy_tree`.
    names:
        Optional ``{tax_id: scientific name}`` map.  Missing names are
        emitted as the empty string (matching Kraken2's behaviour for
        unrecognised taxa).
    ranks:
        Optional ``{tax_id: rank code}`` map (single-letter codes such
        as ``"S"`` / ``"G"`` / ``"F"``).  Missing ranks default to
        ``"-"`` (no rank).
    stats:
        Optional dict that, if supplied, is populated with provenance
        counters: ``total_reads``, ``originally_unclassified``,
        ``demoted_to_unclassified`` (reads that were classified by
        Kraken2 but failed the confidence threshold), and
        ``classified_after_filter``.
    """
    if not (0.0 <= threshold <= 1.0):
        raise ValueError(f"threshold must be in [0.0, 1.0], got {threshold}")
    names = names or {}
    ranks = ranks or {}

    direct: dict[int, int] = {}
    total_reads = 0
    originally_unclassified = 0
    demoted = 0

    for _status, taxid, kmers in iter_kraken2_output(output_path):
        total_reads += 1
        if taxid <= 0:
            originally_unclassified += 1
            continue
        if threshold > 0.0:
            conf = compute_read_confidence(taxid, kmers, tree)
            if conf < threshold:
                demoted += 1
                continue
        direct[taxid] = direct.get(taxid, 0) + 1

    unclassified = originally_unclassified + demoted
    if stats is not None:
        stats["total_reads"] = total_reads
        stats["originally_unclassified"] = originally_unclassified
        stats["demoted_to_unclassified"] = demoted
        stats["classified_after_filter"] = total_reads - unclassified

    # Propagate direct counts up to ancestors to produce clade_reads.
    clade: dict[int, int] = dict(direct)
    for tid, count in direct.items():
        seen: set[int] = {tid}
        parent = tree.get(tid)
        while parent is not None and parent not in seen:
            seen.add(parent)
            clade[parent] = clade.get(parent, 0) + count
            grand = tree.get(parent)
            if grand is None or grand == parent:
                break
            parent = grand

    records: list[TaxonRecord] = []
    # Unclassified row (taxid 0) — kept to mirror the standard report layout.
    if unclassified or total_reads == 0:
        pct = (unclassified / total_reads * 100.0) if total_reads else 0.0
        records.append(
            TaxonRecord(
                tax_id=0,
                name=names.get(0, "unclassified"),
                rank="U",
                clade_reads=unclassified,
                direct_reads=unclassified,
                percentage=pct,
            )
        )
    for tid in sorted(clade.keys()):
        clade_n = clade[tid]
        direct_n = direct.get(tid, 0)
        pct = (clade_n / total_reads * 100.0) if total_reads else 0.0
        records.append(
            TaxonRecord(
                tax_id=tid,
                name=names.get(tid, ""),
                rank=ranks.get(tid, "-"),
                clade_reads=clade_n,
                direct_reads=direct_n,
                percentage=pct,
            )
        )
    return records


# ── Optional names.dmp / rank loaders ────────────────────────────────────────


def load_names_dmp(db_path: str | Path) -> dict[int, str]:
    """Parse ``taxonomy/names.dmp`` and return ``{tax_id: scientific name}``.

    Falls back to ``names.dmp`` in the database root (some Kraken2
    builds do not use a ``taxonomy/`` subdirectory).  Returns an empty
    dict when no names file is found, which is acceptable: downstream
    matrices simply use empty name strings.
    """
    db = Path(db_path)
    candidates = [db / "taxonomy" / "names.dmp", db / "names.dmp"]
    names_path = next((p for p in candidates if p.exists()), None)
    if names_path is None:
        logger.warning("names.dmp not found under %s; names will be empty", db)
        return {}
    out: dict[int, str] = {}
    with open(names_path) as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("\t|")]
            if len(parts) < 4:
                continue
            try:
                tid = int(parts[0])
            except ValueError:
                continue
            name = parts[1]
            name_class = parts[3].rstrip("\t|").strip()
            if name_class == "scientific name":
                out[tid] = name
    return out


def load_ranks_dmp(db_path: str | Path) -> dict[int, str]:
    """Parse ``taxonomy/nodes.dmp`` and return ``{tax_id: rank code}``.

    Maps Kraken2's full rank strings to the single-letter codes used in
    this codebase (``"D" / "P" / "C" / "O" / "F" / "G" / "S"``); other
    ranks map to ``"-"``.
    """
    db = Path(db_path)
    candidates = [db / "taxonomy" / "nodes.dmp", db / "nodes.dmp"]
    nodes_path = next((p for p in candidates if p.exists()), None)
    if nodes_path is None:
        return {}
    out: dict[int, str] = {}
    with open(nodes_path) as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("\t|")]
            if len(parts) < 3:
                continue
            try:
                tid = int(parts[0])
            except ValueError:
                continue
            rank_str = parts[2].lower()
            out[tid] = _KRAKEN_RANK_CODES.get(rank_str, "-")
    return out


# ── Helpers for sample / file pairing ────────────────────────────────────────


def sample_id_from_output(path: str | Path) -> str:
    """Derive a sample ID from a Kraken2 per-read output filename.

    Strips the canonical ``.kraken2.output.txt`` suffix; falls back to
    :py:meth:`pathlib.Path.stem` for unrecognised names.
    """
    name = Path(path).name
    suffix = ".kraken2.output.txt"
    if name.endswith(suffix):
        return name[: -len(suffix)]
    return Path(path).stem


def format_threshold_suffix(threshold: float) -> str:
    """Render a threshold as a filesystem-safe suffix (e.g. ``conf0p50``).

    The dot in the threshold is replaced with ``p`` so that filenames
    like ``taxa_matrix_raw_conf0p50.tsv`` remain easy to glob.  Two
    decimal places of precision are kept (Kraken2 confidences below
    that resolution are not statistically meaningful given typical
    150 bp k-mer counts).
    """
    return f"conf{threshold:0.2f}".replace(".", "p")


def map_outputs_to_samples(
    output_paths: Iterable[str | Path],
) -> dict[str, Path]:
    """Build a ``{sample_id: path}`` map from output filenames.

    Duplicate sample IDs emit a warning; the latter wins.
    """
    out: dict[str, Path] = {}
    for p in output_paths:
        sid = sample_id_from_output(p)
        if sid in out:
            logger.warning(
                "Duplicate sample_id '%s' from kraken2 outputs (%s replaces %s)",
                sid, p, out[sid],
            )
        out[sid] = Path(p)
    return out
