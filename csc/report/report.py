"""Core rendering logic for the CSC static HTML contamination report.

The report is organised into the following sections, in the order
requested by the upstream issue:

1. **Executive Summary** – cohort-level headline metrics.
2. **Methods** – sample processing, classification, aggregation and
   outlier-detection procedures, plus a clear statement distinguishing
   (a) reads extracted, (b) reads classified by Kraken2, and
   (c) total sequenced reads.
3. **Results** – summary tables (both CPM *and* absolute burden, with
   explicit denominators in column headers and figure captions) plus
   inline SVG visualisations.
4. **Variant-Calling Impact** – flags samples whose *absolute* non-human
   burden exceeds a configurable threshold.  Absolute burden (per
   million total sequenced reads) is the scientifically appropriate
   denominator when judging downstream impact on variant calling,
   assembly, or expression quantitation.  See Natarajan et al.,
   *Nat. Biotechnol.* 41, 1520–1530 (2023).
5. **Discussion** – interpretation caveats.
6. **Methods Transparency Checklist** – versioned references, filter
   settings, and download links for every input matrix so the report
   is fully reproducible and review-proof.

The rendering uses only the Python standard library and emits one
self-contained HTML file with inline CSS and inline SVG plots.

AI assistance acknowledgment: This module was developed with AI
assistance.  Best practices in the bioinformatics field should always
take precedence over specific implementation details.
"""

from __future__ import annotations

import csv
import html
import json
import logging
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# Schema version for the HTML report output.  Bump on breaking changes
# to the structure of ``report_manifest.json``.
REPORT_SCHEMA_VERSION = "1.0"

# Default absolute-burden threshold (reads per million total sequenced
# reads) above which a sample is flagged as potentially impacting
# downstream variant calling.  0.1% = 1,000 per-million is a commonly
# cited heuristic; configurable via ``generate_html_report``.
DEFAULT_VARIANT_IMPACT_THRESHOLD_PPM = 1000.0


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


@dataclass
class Matrix:
    """A parsed taxa-by-sample matrix.

    ``values`` is ``{tax_id: {sample_id: float | None}}`` where ``None``
    represents ``NA`` (no denominator available, e.g. a sample missing
    from the idxstats sidecars in the absolute-burden matrix).
    """

    sample_ids: list[str]
    tax_ids: list[int]
    tax_names: dict[int, str]
    tax_domains: dict[int, str]
    values: dict[int, dict[str, float | None]]


@dataclass
class ReportInputs:
    """All parsed inputs required to render the HTML report."""

    aggregate_metadata: dict[str, Any]
    matrix_raw: Matrix
    matrix_cpm: Matrix
    matrix_abs: Matrix | None = None
    detect_summary: dict[str, Any] | None = None
    flagged_samples: list[dict[str, Any]] = field(default_factory=list)
    # Optional parallel detection pass on the absolute-burden matrix.
    # Populated when ``<detect_dir>/abs/qc_summary.json`` exists,
    # i.e. csc-detect ran its abs side pass (the default when an
    # absolute-burden sibling matrix is available).  This pass uses
    # 'reads per million total sequenced reads' as its denominator and
    # therefore highlights contamination episodes that compositional
    # CPM detection can mask when host-depletion rates differ between
    # samples.
    abs_detect_summary: dict[str, Any] | None = None
    abs_flagged_samples: list[dict[str, Any]] = field(default_factory=list)
    # Source paths (preserved for the transparency checklist and for
    # the relative download links that appear in the rendered HTML).
    aggregate_dir: Path | None = None
    detect_dir: Path | None = None
    abs_detect_dir: Path | None = None


def _parse_matrix(path: Path) -> Matrix:
    """Parse a CSC taxa matrix TSV, preserving ``NA`` as ``None``.

    Handles both the two-metadata-column layout
    (``tax_id, name, sample1, …``) and the domain-annotated
    three-column layout (``tax_id, name, domain, sample1, …``).
    """
    if not path.exists():
        raise FileNotFoundError(f"Matrix file not found: {path}")

    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"Empty matrix file: {path}")

        has_domain = len(header) > 2 and header[2] == "domain"
        meta_cols = 3 if has_domain else 2
        if len(header) < meta_cols + 1:
            raise ValueError(
                f"Matrix must have at least {meta_cols + 1} columns "
                f"(tax_id, name,{' domain,' if has_domain else ''} "
                f"and one sample), got {len(header)}: {path}"
            )

        sample_ids = header[meta_cols:]
        tax_ids: list[int] = []
        tax_names: dict[int, str] = {}
        tax_domains: dict[int, str] = {}
        values: dict[int, dict[str, float | None]] = {}

        for lineno, cols in enumerate(reader, 2):
            if len(cols) < meta_cols + 1:
                logger.warning(
                    "Skipping line %d in %s: too few columns", lineno, path
                )
                continue
            try:
                tid = int(cols[0])
            except ValueError:
                logger.warning(
                    "Skipping line %d in %s: non-integer tax_id", lineno, path
                )
                continue
            tax_ids.append(tid)
            tax_names[tid] = cols[1]
            if has_domain:
                tax_domains[tid] = cols[2]
            row: dict[str, float | None] = {}
            for i, sid in enumerate(sample_ids):
                raw = cols[meta_cols + i] if meta_cols + i < len(cols) else ""
                if raw == "" or raw == "NA":
                    row[sid] = None
                else:
                    try:
                        row[sid] = float(raw)
                    except ValueError:
                        row[sid] = None
            values[tid] = row

    return Matrix(
        sample_ids=sample_ids,
        tax_ids=tax_ids,
        tax_names=tax_names,
        tax_domains=tax_domains,
        values=values,
    )


def _parse_flagged(path: Path) -> list[dict[str, Any]]:
    """Parse a ``flagged_samples.tsv`` emitted by :mod:`csc.detect`."""
    flagged: list[dict[str, Any]] = []
    if not path.exists():
        return flagged
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            flagged.append(dict(row))
    return flagged


def load_inputs(
    aggregate_dir: str | Path,
    detect_dir: str | Path | None = None,
) -> ReportInputs:
    """Load all files the report needs from existing pipeline outputs.

    Parameters
    ----------
    aggregate_dir:
        Directory containing ``taxa_matrix_raw.tsv``,
        ``taxa_matrix_cpm.tsv`` and ``aggregation_metadata.json``
        (and optionally ``taxa_matrix_abs.tsv``).
    detect_dir:
        Optional directory containing ``flagged_samples.tsv`` and
        ``qc_summary.json`` from :mod:`csc.detect`.  When the directory
        also contains an ``abs/`` subdirectory (written by
        ``csc-detect``'s absolute-burden side pass), its
        ``qc_summary.json`` and ``flagged_samples.tsv`` are loaded as
        well so the report can summarise both passes.

    Raises
    ------
    FileNotFoundError
        If any required aggregate output is missing.
    """
    aggregate_dir = Path(aggregate_dir)
    if not aggregate_dir.is_dir():
        raise FileNotFoundError(
            f"aggregate_dir does not exist or is not a directory: {aggregate_dir}"
        )

    meta_path = aggregate_dir / "aggregation_metadata.json"
    if not meta_path.exists():
        raise FileNotFoundError(
            f"Missing aggregation_metadata.json in {aggregate_dir}"
        )
    with open(meta_path) as fh:
        aggregate_metadata = json.load(fh)

    matrix_raw = _parse_matrix(aggregate_dir / "taxa_matrix_raw.tsv")
    matrix_cpm = _parse_matrix(aggregate_dir / "taxa_matrix_cpm.tsv")

    abs_path = aggregate_dir / "taxa_matrix_abs.tsv"
    matrix_abs: Matrix | None = None
    if abs_path.exists():
        matrix_abs = _parse_matrix(abs_path)

    detect_summary: dict[str, Any] | None = None
    flagged: list[dict[str, Any]] = []
    detect_dir_path: Path | None = None
    abs_detect_summary: dict[str, Any] | None = None
    abs_flagged: list[dict[str, Any]] = []
    abs_detect_dir_path: Path | None = None
    if detect_dir is not None:
        detect_dir_path = Path(detect_dir)
        qc_path = detect_dir_path / "qc_summary.json"
        if qc_path.exists():
            with open(qc_path) as fh:
                detect_summary = json.load(fh)
        flagged = _parse_flagged(detect_dir_path / "flagged_samples.tsv")

        # csc-detect writes its abs side-pass results to <detect_dir>/abs/.
        # Pick them up automatically so the report can compare CPM-based
        # and absolute-burden-based flag sets side-by-side.
        abs_dir = detect_dir_path / "abs"
        abs_qc_path = abs_dir / "qc_summary.json"
        if abs_qc_path.exists():
            abs_detect_dir_path = abs_dir
            with open(abs_qc_path) as fh:
                abs_detect_summary = json.load(fh)
            abs_flagged = _parse_flagged(abs_dir / "flagged_samples.tsv")

    return ReportInputs(
        aggregate_metadata=aggregate_metadata,
        matrix_raw=matrix_raw,
        matrix_cpm=matrix_cpm,
        matrix_abs=matrix_abs,
        detect_summary=detect_summary,
        flagged_samples=flagged,
        abs_detect_summary=abs_detect_summary,
        abs_flagged_samples=abs_flagged,
        aggregate_dir=aggregate_dir,
        detect_dir=detect_dir_path,
        abs_detect_dir=abs_detect_dir_path,
    )


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------


def _column_sum(matrix: Matrix, sample_id: str) -> float:
    total = 0.0
    for tid in matrix.tax_ids:
        v = matrix.values[tid].get(sample_id)
        if v is not None and not math.isnan(v):
            total += v
    return total


def _shannon_index(matrix: Matrix, sample_id: str) -> float:
    """Shannon-Wiener diversity index (natural log) for one sample.

    Computed on raw counts so the value is invariant under the choice
    of CPM or absolute normalisation.  Returns 0 for samples with zero
    classified reads (or when every classified read hits a single
    taxon).
    """
    counts = []
    for tid in matrix.tax_ids:
        v = matrix.values[tid].get(sample_id)
        if v is not None and v > 0:
            counts.append(v)
    total = sum(counts)
    if total <= 0:
        return 0.0
    h = 0.0
    for c in counts:
        p = c / total
        h -= p * math.log(p)
    return h


def _simpson_index(matrix: Matrix, sample_id: str) -> float:
    """Gini–Simpson diversity index (1 − Σ p²) on raw counts."""
    counts = []
    for tid in matrix.tax_ids:
        v = matrix.values[tid].get(sample_id)
        if v is not None and v > 0:
            counts.append(v)
    total = sum(counts)
    if total <= 0:
        return 0.0
    return 1.0 - sum((c / total) ** 2 for c in counts)


def _richness(matrix: Matrix, sample_id: str) -> int:
    """Number of taxa with at least one classified read."""
    n = 0
    for tid in matrix.tax_ids:
        v = matrix.values[tid].get(sample_id)
        if v is not None and v > 0:
            n += 1
    return n


def _top_taxa(
    matrix: Matrix, sample_id: str, top_n: int = 10
) -> list[tuple[int, float]]:
    """Return the top-N (tax_id, value) pairs for a sample, desc by value.

    Entries with ``None`` (NA) values are ignored.
    """
    pairs: list[tuple[int, float]] = []
    for tid in matrix.tax_ids:
        v = matrix.values[tid].get(sample_id)
        if v is not None and v > 0:
            pairs.append((tid, v))
    pairs.sort(key=lambda kv: kv[1], reverse=True)
    return pairs[:top_n]


def _domain_composition(
    matrix: Matrix, sample_id: str
) -> dict[str, float]:
    """Aggregate a sample's values by taxonomic domain.

    When the matrix lacks a ``domain`` annotation, all taxa are pooled
    under the key ``"Unannotated"`` – this keeps the stacked-bar plot
    correct even if the aggregation was run without ``--db-path``.
    """
    comp: dict[str, float] = {}
    for tid in matrix.tax_ids:
        v = matrix.values[tid].get(sample_id)
        if v is None or v <= 0:
            continue
        dom = matrix.tax_domains.get(tid, "Unannotated")
        comp[dom] = comp.get(dom, 0.0) + v
    return comp


# ---------------------------------------------------------------------------
# Rendering helpers
# ---------------------------------------------------------------------------


# A small color palette chosen to be colour-blind friendly (Okabe–Ito)
# with a deterministic ordering.  If a sample has more domains than
# palette entries, we cycle — acceptable for the small number of
# domains CSC emits (Human, Bacteria, Archaea, Fungi, Protists,
# Viruses, UniVec_Core, Metazoa_other, Viridiplantae, Unclassified).
_DOMAIN_PALETTE = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#999999",  # neutral grey
    "#000000",  # black
    "#66C2A5",  # teal
]


def _domain_colour(domain: str, idx: int) -> str:
    return _DOMAIN_PALETTE[idx % len(_DOMAIN_PALETTE)]


def _fmt_num(v: float | None, *, ndigits: int = 2) -> str:
    """Render a float for a human-readable table cell, preserving NA."""
    if v is None:
        return "NA"
    if math.isnan(v):
        return "NA"
    if v == 0:
        return "0"
    if abs(v) >= 1_000_000:
        return f"{v:,.0f}"
    if abs(v) >= 1:
        return f"{v:,.{ndigits}f}"
    return f"{v:.{max(ndigits, 4)}g}"


def _fmt_int(v: int | None) -> str:
    if v is None:
        return "NA"
    return f"{v:,}"


def _stacked_bar_svg(
    compositions: list[tuple[str, dict[str, float]]],
    title: str,
    caption: str,
    *,
    width: int = 720,
    bar_height: int = 22,
    label_width: int = 180,
    normalize: bool = True,
) -> str:
    """Render a stacked-bar SVG figure (one row per sample).

    ``compositions`` is a list of ``(sample_id, {domain: value})``.
    When ``normalize=True`` each sample's bar sums to 100 %; the
    caption should then describe the composition.  When
    ``normalize=False`` the bars share an absolute x-axis and the
    caption should state the denominator explicitly (e.g. "per
    million total sequenced reads").
    """
    n = len(compositions)
    if n == 0:
        return (
            '<div class="figure"><p class="fig-title">'
            + html.escape(title)
            + "</p><p><em>No samples to plot.</em></p></div>"
        )
    # Collect all domains in a deterministic order (by total abundance desc).
    domain_totals: dict[str, float] = {}
    for _, comp in compositions:
        for d, v in comp.items():
            domain_totals[d] = domain_totals.get(d, 0.0) + v
    domains = sorted(domain_totals, key=lambda d: -domain_totals[d])

    plot_width = max(100, width - label_width - 30)
    height = n * (bar_height + 6) + 60

    x_max = 1.0 if normalize else max(
        (sum(c.values()) for _, c in compositions),
        default=1.0,
    ) or 1.0

    parts: list[str] = []
    parts.append(f'<div class="figure">')
    parts.append(f'<p class="fig-title">{html.escape(title)}</p>')
    parts.append(
        f'<svg viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img" '
        f'aria-label="{html.escape(title)}">'
    )

    # Bars
    for i, (sid, comp) in enumerate(compositions):
        y = 20 + i * (bar_height + 6)
        parts.append(
            f'<text x="{label_width - 5}" y="{y + bar_height * 0.7:.1f}" '
            f'text-anchor="end" font-size="12" '
            f'font-family="sans-serif">{html.escape(sid)}</text>'
        )
        total = sum(comp.values())
        if normalize and total > 0:
            scale = plot_width / 1.0
            running = 0.0
            for j, d in enumerate(domains):
                frac = comp.get(d, 0.0) / total
                if frac <= 0:
                    continue
                x = label_width + running * scale
                w = frac * scale
                parts.append(
                    f'<rect x="{x:.2f}" y="{y}" width="{w:.2f}" '
                    f'height="{bar_height}" '
                    f'fill="{_domain_colour(d, j)}">'
                    f'<title>{html.escape(sid)} / {html.escape(d)}: '
                    f'{frac * 100:.2f}%</title></rect>'
                )
                running += frac
        else:
            scale = plot_width / x_max if x_max > 0 else 0
            running = 0.0
            for j, d in enumerate(domains):
                val = comp.get(d, 0.0)
                if val <= 0:
                    continue
                x = label_width + running * scale
                w = val * scale
                parts.append(
                    f'<rect x="{x:.2f}" y="{y}" width="{w:.2f}" '
                    f'height="{bar_height}" '
                    f'fill="{_domain_colour(d, j)}">'
                    f'<title>{html.escape(sid)} / {html.escape(d)}: '
                    f'{val:.2f}</title></rect>'
                )
                running += val

    # X-axis label
    axis_y = height - 30
    parts.append(
        f'<line x1="{label_width}" y1="{axis_y}" '
        f'x2="{label_width + plot_width}" y2="{axis_y}" '
        f'stroke="#444" stroke-width="1" />'
    )
    if normalize:
        for frac, lbl in [(0.0, "0%"), (0.25, "25%"), (0.5, "50%"),
                          (0.75, "75%"), (1.0, "100%")]:
            x = label_width + frac * plot_width
            parts.append(
                f'<line x1="{x:.2f}" y1="{axis_y}" x2="{x:.2f}" '
                f'y2="{axis_y + 4}" stroke="#444" stroke-width="1" />'
                f'<text x="{x:.2f}" y="{axis_y + 16}" '
                f'text-anchor="middle" font-size="11" '
                f'font-family="sans-serif">{lbl}</text>'
            )

    # Legend
    legend_y = height - 12
    legend_x = label_width
    for j, d in enumerate(domains):
        parts.append(
            f'<rect x="{legend_x}" y="{legend_y - 9}" width="10" height="10" '
            f'fill="{_domain_colour(d, j)}" />'
        )
        parts.append(
            f'<text x="{legend_x + 14}" y="{legend_y - 1}" font-size="11" '
            f'font-family="sans-serif">{html.escape(d)}</text>'
        )
        legend_x += 14 + 7 * len(d) + 18

    parts.append("</svg>")
    parts.append(f'<p class="fig-caption">{html.escape(caption)}</p>')
    parts.append("</div>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Section renderers
# ---------------------------------------------------------------------------


_PAGE_STYLE = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
       "Helvetica Neue", Arial, sans-serif; color: #222; max-width: 1080px;
       margin: 2em auto; padding: 0 1em; line-height: 1.55; }
h1 { border-bottom: 2px solid #333; padding-bottom: 0.2em; }
h2 { margin-top: 2em; border-bottom: 1px solid #ccc; padding-bottom: 0.15em; }
h3 { margin-top: 1.4em; }
table { border-collapse: collapse; margin: 0.6em 0 1.2em; width: 100%; }
th, td { border: 1px solid #ccc; padding: 4px 8px; font-size: 13px;
         text-align: right; }
th:first-child, td:first-child { text-align: left; }
th { background: #f4f4f4; }
tr:nth-child(even) td { background: #fafafa; }
.meta-key { color: #555; font-weight: 500; }
.figure { margin: 1.2em 0; }
.fig-title { font-weight: 600; margin-bottom: 0.3em; }
.fig-caption { font-size: 12.5px; color: #555; margin-top: 0.3em;
               font-style: italic; }
.callout { background: #fff8e1; border-left: 4px solid #ffb300;
           padding: 0.6em 0.9em; margin: 1em 0; border-radius: 3px; }
.callout.danger { background: #ffebee; border-left-color: #d32f2f; }
.callout.ok { background: #e8f5e9; border-left-color: #2e7d32; }
code { background: #f1f1f1; padding: 1px 4px; border-radius: 3px; }
.download-list li { margin: 0.2em 0; }
.denom { color: #555; font-weight: normal; font-size: 11.5px; display: block; }
.checklist { list-style: none; padding-left: 0; }
.checklist li::before { content: "\\2611  "; color: #2e7d32;
                        margin-right: 0.4em; }
.footer { margin-top: 3em; font-size: 12px; color: #666;
          border-top: 1px solid #ddd; padding-top: 0.6em; }
""".strip()


def _render_executive_summary(
    inputs: ReportInputs,
    per_sample_stats: list[dict[str, Any]],
    variant_flagged: list[dict[str, Any]],
    threshold_ppm: float,
) -> str:
    meta = inputs.aggregate_metadata
    n_samples = meta.get("sample_count", len(inputs.matrix_raw.sample_ids))
    n_taxa = meta.get("taxon_count", len(inputs.matrix_raw.tax_ids))
    abs_enabled = inputs.matrix_abs is not None

    # Cohort-level aggregate numbers
    total_classified = sum(s["classified_reads"] for s in per_sample_stats)
    median_richness = _median([s["richness"] for s in per_sample_stats])

    nonhuman_pct_values = [
        s["nonhuman_pct_of_classified"]
        for s in per_sample_stats
        if s["nonhuman_pct_of_classified"] is not None
    ]
    median_nonhuman_pct_classified = (
        _median(nonhuman_pct_values) if nonhuman_pct_values else None
    )

    if abs_enabled:
        nonhuman_burden_values = [
            s["nonhuman_abs_ppm_total"]
            for s in per_sample_stats
            if s["nonhuman_abs_ppm_total"] is not None
        ]
        median_burden = (
            _median(nonhuman_burden_values) if nonhuman_burden_values else None
        )
    else:
        median_burden = None

    items: list[str] = []
    items.append(
        f"<li><span class='meta-key'>Samples analysed:</span> {n_samples}</li>"
    )
    items.append(
        f"<li><span class='meta-key'>Distinct taxa observed:</span> {n_taxa}</li>"
    )
    items.append(
        f"<li><span class='meta-key'>Total classified reads (cohort):</span> "
        f"{_fmt_int(int(total_classified))}</li>"
    )
    if median_nonhuman_pct_classified is not None:
        items.append(
            f"<li><span class='meta-key'>Median % non-human "
            f"<em>among classified reads</em> (CPM denominator):</span> "
            f"{median_nonhuman_pct_classified:.3f}%</li>"
        )
    if median_burden is not None:
        items.append(
            f"<li><span class='meta-key'>Median non-human absolute burden "
            f"<em>per million total sequenced reads</em>:</span> "
            f"{median_burden:,.2f} ppm "
            f"({median_burden / 1e4:.4f}%)</li>"
        )
    items.append(
        f"<li><span class='meta-key'>Median species richness per sample:</span> "
        f"{median_richness}</li>"
    )
    if inputs.detect_summary is not None:
        n_flag = len(inputs.detect_summary.get("flagged_samples", []))
        cpm_label = (
            inputs.detect_summary.get("matrix_type") or "primary"
        )
        items.append(
            f"<li><span class='meta-key'>Samples flagged by "
            f"<code>csc-detect</code> (<em>{html.escape(str(cpm_label))}"
            f"</em> matrix):</span> {n_flag}</li>"
        )
    if inputs.abs_detect_summary is not None:
        abs_flag = len(inputs.abs_detect_summary.get("flagged_samples", []))
        items.append(
            "<li><span class='meta-key'>Samples flagged by "
            "<code>csc-detect</code> (<em>absolute-burden</em> "
            "side pass &mdash; ppm of total sequenced reads):</span> "
            f"{abs_flag}</li>"
        )
    items.append(
        f"<li><span class='meta-key'>Samples exceeding variant-calling "
        f"impact threshold "
        f"(&gt; {threshold_ppm:,.0f} ppm of total reads, i.e. "
        f"{threshold_ppm / 1e4:.3f}% of sequencing):</span> "
        f"{len(variant_flagged) if abs_enabled else 'N/A (no idxstats)'}</li>"
    )

    callout = ""
    if not abs_enabled:
        callout = (
            "<div class='callout'>⚠ <strong>Absolute-burden matrix "
            "(<code>taxa_matrix_abs.tsv</code>) not available.</strong> "
            "Re-run <code>csc-aggregate --idxstats …</code> with the "
            "per-sample <code>reads_summary.json</code> sidecars emitted "
            "by <code>csc-extract</code> to enable absolute-burden "
            "reporting and the variant-calling impact section. "
            "Without it, this report can only describe contaminant "
            "<em>composition</em> – not the absolute fraction of "
            "sequencing affected.</div>"
        )

    return (
        "<h2>1. Executive Summary</h2>"
        "<p>This report summarises non-human content detected in a cohort "
        "of human whole-genome sequencing (WGS) samples.  It was "
        "generated directly from <code>csc-aggregate</code> outputs and "
        "is intended to be suitable both for pipeline QC and for "
        "inclusion in manuscripts.  All tables and figures label their "
        "denominators explicitly (composition within contaminants vs. "
        "fraction of total sequencing); see §2 for the distinction "
        "between <em>reads extracted</em>, <em>reads classified</em>, "
        "and <em>total reads sequenced</em>.</p>"
        + callout
        + "<ul>" + "\n".join(items) + "</ul>"
    )


def _median(values: list[float]) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    n = len(s)
    m = n // 2
    if n % 2 == 1:
        return float(s[m])
    return (s[m - 1] + s[m]) / 2.0


def _render_methods(inputs: ReportInputs) -> str:
    meta = inputs.aggregate_metadata
    schema = meta.get("schema_version", "unknown")
    min_reads = meta.get("min_reads", "unknown")
    rank_filter = meta.get("rank_filter", [])
    domain_annotated = meta.get("domain_annotated", False)
    abs_enabled = meta.get("absolute_burden_enabled", False)

    parts = [
        "<h2>2. Methods</h2>",
        "<h3>2.1 Sample processing and classification</h3>",
        "<ol>"
        "<li><strong>Extract</strong> (<code>csc-extract</code>): "
        "streaming extraction of unmapped reads (and optionally "
        "low-MAPQ reads) from aligned BAM/CRAM files using "
        "<code>samtools</code>.  Per-sample totals from "
        "<code>samtools idxstats</code> are captured in a "
        "<code>reads_summary.json</code> sidecar <em>before</em> any "
        "downstream filtering so they can serve as the authoritative "
        "absolute-burden denominator.</li>"
        "<li><strong>Classify</strong> (<code>csc-classify</code>): "
        "taxonomic assignment of the extracted reads with Kraken2 "
        "against a curated reference database.  The recommended "
        "database is <strong>PrackenDB</strong> (one NCBI reference "
        "genome per species for bacteria, archaea, fungi, protists, "
        "viruses, human, and UniVec Core).</li>"
        "<li><strong>Aggregate</strong> (<code>csc-aggregate</code>): "
        "per-taxon read counts are collated across samples into three "
        "matrices with distinct, explicitly labelled denominators "
        "(see §2.2).  Rank-filtered matrices "
        f"(<code>{', '.join(rank_filter) or 'none'}</code>) are emitted "
        "alongside the unfiltered matrices.</li>"
        "<li><strong>Detect</strong> (<code>csc-detect</code>, optional): "
        "statistical outlier detection using MAD, IQR, and/or a "
        "two-component GMM to flag samples whose profile deviates "
        "from the cohort.</li>"
        "</ol>",
        "<h3>2.2 Denominator provenance – three reads, three numbers</h3>",
        "<p>Three distinct read counts appear in this pipeline and in "
        "this report; they are <strong>not</strong> interchangeable:</p>"
        "<table><thead><tr>"
        "<th>Quantity</th><th>Where counted</th>"
        "<th>Why it differs from the others</th></tr></thead><tbody>"
        "<tr><td><strong>Total sequenced reads</strong></td>"
        "<td><code>samtools idxstats</code> on the input BAM/CRAM "
        "(captured at extract time into "
        "<code>reads_summary.json</code>)</td>"
        "<td>Largest number; includes mapped + unmapped.  "
        "Denominator for <strong>absolute contamination burden</strong> – "
        "the scientifically appropriate metric for judging impact on "
        "variant calling, assembly, or expression quantitation.</td></tr>"
        "<tr><td><strong>Reads extracted for classification</strong></td>"
        "<td>Output FASTQs from <code>csc-extract</code></td>"
        "<td>Subset of total reads: unmapped (+ optional low-MAPQ).  "
        "Adapter trimming / host filtering may further reduce this "
        "count before classification.</td></tr>"
        "<tr><td><strong>Reads classified by Kraken2</strong></td>"
        "<td>Sum of direct-read counts across all taxa in each "
        "sample's Kraken2 report</td>"
        "<td>Subset of extracted reads: only those to which Kraken2 "
        "could assign a taxonomic label.  "
        "Denominator for <strong>compositional CPM</strong> – i.e., "
        "&ldquo;what fraction of the <em>non-human</em> content is "
        "this taxon?&rdquo;</td></tr>"
        "</tbody></table>",
        "<p>Accordingly the report presents three normalisation tiers:</p>"
        "<table><thead><tr><th>Matrix</th><th>Denominator</th>"
        "<th>Answers</th></tr></thead><tbody>"
        "<tr><td><code>taxa_matrix_raw.tsv</code></td><td>None (integer "
        "counts)</td><td>How many reads hit this taxon?</td></tr>"
        "<tr><td><code>taxa_matrix_cpm.tsv</code></td>"
        "<td>Sum of classified direct reads (pre-filter)</td>"
        "<td>What fraction of the non-human content is this taxon? "
        "(compositional; use for stacked bars and diversity indices.)</td></tr>"
        "<tr><td><code>taxa_matrix_abs.tsv</code></td>"
        "<td>Total sequenced reads (from idxstats)</td>"
        "<td>What fraction of <em>all sequencing</em> is this taxon? "
        "(absolute burden; use for variant-calling impact and "
        "manuscript QC tables.)</td></tr>"
        "</tbody></table>",
        "<h3>2.3 Aggregation parameters</h3>"
        "<table><tbody>"
        f"<tr><td>Aggregation schema version</td><td><code>{html.escape(str(schema))}</code></td></tr>"
        f"<tr><td>Minimum reads per taxon</td><td><code>{html.escape(str(min_reads))}</code></td></tr>"
        f"<tr><td>Rank filter</td><td><code>{html.escape(', '.join(rank_filter) or 'none')}</code></td></tr>"
        f"<tr><td>Lineage-aware domain annotation</td>"
        f"<td><code>{'enabled' if domain_annotated else 'disabled'}</code></td></tr>"
        f"<tr><td>Absolute-burden matrix available</td>"
        f"<td><code>{'yes' if abs_enabled else 'no'}</code></td></tr>"
        "</tbody></table>",
    ]

    if inputs.detect_summary is not None:
        ds = inputs.detect_summary
        parts.append("<h3>2.4 Outlier-detection parameters</h3>")
        parts.append(
            "<p><code>csc-detect</code> applies statistical outlier "
            "detection (MAD, IQR, and/or a two-component GMM) per "
            "taxon, flagging samples whose abundance exceeds the "
            "cohort distribution.  The detector identifies "
            "<em>upper-tail</em> outliers only; under-representation "
            "is not informative for contamination.  The same "
            "parameters are applied to every matrix variant analysed "
            "(rank-filtered matrices, confidence-tier matrices, and "
            "the absolute-burden side pass).</p>"
            "<table><tbody>"
        )
        for k in (
            "method",
            "matrix_type",
            "min_reads",
            "mad_multiplier",
            "iqr_multiplier",
            "gmm_threshold",
            "subtract_background",
            "kitome_taxa",
        ):
            if k in ds:
                parts.append(
                    f"<tr><td>{html.escape(k)}</td>"
                    f"<td><code>{html.escape(json.dumps(ds[k]))}</code></td></tr>"
                )
        parts.append("</tbody></table>")

        if inputs.abs_detect_summary is not None:
            parts.append(
                "<h3>2.5 Absolute-burden side pass</h3>"
                "<p>In addition to the primary detection pass, "
                "<code>csc-detect</code> automatically ran a parallel "
                "pass on <code>taxa_matrix_abs.tsv</code> "
                "(reads per million <em>total sequenced reads</em>, "
                "from <code>samtools idxstats</code>).  The "
                "compositional CPM denominator is invariant to how "
                "many reads were classified, so two libraries with "
                "very different host-depletion efficiencies can "
                "produce indistinguishable CPM profiles even when "
                "their absolute non-human burdens differ by orders of "
                "magnitude.  The absolute-burden pass uses a "
                "denominator that is independent of host depletion and "
                "therefore catches contamination episodes the CPM "
                "pass can mask.  Results are written to "
                "<code>detect/abs/</code> and are summarised "
                "side-by-side with the primary pass in §3.3.</p>"
            )
        else:
            parts.append(
                "<h3>2.5 Absolute-burden side pass</h3>"
                "<div class='callout'>The absolute-burden side pass "
                "was <strong>not run</strong> for this report (either "
                "no <code>taxa_matrix_abs.tsv</code> was available, "
                "or detection was invoked with "
                "<code>--no-abs-detection</code>).  Re-run "
                "<code>csc-aggregate</code> with the per-sample "
                "<code>reads_summary.json</code> sidecars and re-run "
                "<code>csc-detect</code> to enable.  Without it, "
                "outlier detection is purely compositional and may "
                "miss contamination episodes that are masked by "
                "host-depletion-rate differences.</div>"
            )

    return "\n".join(parts)


def _per_sample_stats(
    inputs: ReportInputs,
    human_tax_ids: tuple[int, ...] = (9606,),
) -> list[dict[str, Any]]:
    """Compute one stats row per sample."""
    rows: list[dict[str, Any]] = []
    raw = inputs.matrix_raw
    cpm = inputs.matrix_cpm
    abs_m = inputs.matrix_abs
    provenance = inputs.aggregate_metadata.get("sample_provenance", {}) or {}

    for sid in raw.sample_ids:
        classified = _column_sum(raw, sid)

        # Non-human raw reads = classified - human counts (sum of 9606 etc.)
        human_reads = 0.0
        for h in human_tax_ids:
            v = raw.values.get(h, {}).get(sid)
            if v is not None and v > 0:
                human_reads += v
        nonhuman_reads = max(0.0, classified - human_reads)

        nonhuman_pct_classified = (
            (nonhuman_reads / classified) * 100.0 if classified > 0 else None
        )

        richness = _richness(raw, sid)
        shannon = _shannon_index(raw, sid)
        simpson = _simpson_index(raw, sid)

        prov = provenance.get(sid, {}) if isinstance(provenance, dict) else {}
        total_reads = prov.get("total_reads")
        total_unmapped = prov.get("total_unmapped")

        # Absolute non-human burden (reads per million total sequenced)
        nonhuman_abs_ppm_total: float | None = None
        if abs_m is not None:
            total_abs = 0.0
            any_value = False
            for tid in abs_m.tax_ids:
                if tid in human_tax_ids:
                    continue
                v = abs_m.values.get(tid, {}).get(sid)
                if v is None:
                    continue
                any_value = True
                if v > 0:
                    total_abs += v
            nonhuman_abs_ppm_total = total_abs if any_value else None

        rows.append(
            {
                "sample_id": sid,
                "total_reads": total_reads,
                "total_unmapped": total_unmapped,
                "classified_reads": classified,
                "human_reads": human_reads,
                "nonhuman_reads": nonhuman_reads,
                "nonhuman_pct_of_classified": nonhuman_pct_classified,
                "nonhuman_abs_ppm_total": nonhuman_abs_ppm_total,
                "richness": richness,
                "shannon": shannon,
                "simpson": simpson,
            }
        )

    return rows


def _render_results(
    inputs: ReportInputs,
    per_sample_stats: list[dict[str, Any]],
    top_n: int,
) -> str:
    abs_enabled = inputs.matrix_abs is not None

    # --- Per-sample summary table ---
    header_cells = [
        "Sample",
        "Total sequenced reads<span class='denom'>idxstats</span>",
        "Classified reads<span class='denom'>Kraken2 direct</span>",
        "Non-human reads<span class='denom'>classified − Homo sapiens</span>",
        "% Non-human<span class='denom'>of classified (CPM denom.)</span>",
    ]
    if abs_enabled:
        header_cells.append(
            "Non-human burden "
            "<span class='denom'>ppm of total sequenced reads</span>"
        )
    header_cells.extend([
        "Richness<span class='denom'>taxa with ≥1 read</span>",
        "Shannon H<span class='denom'>natural log</span>",
        "Simpson (1−D)<span class='denom'>Gini–Simpson</span>",
    ])

    rows_html: list[str] = []
    for s in per_sample_stats:
        cells = [
            html.escape(s["sample_id"]),
            _fmt_int(
                int(s["total_reads"]) if s["total_reads"] is not None else None
            ),
            _fmt_int(int(s["classified_reads"])),
            _fmt_int(int(s["nonhuman_reads"])),
            (
                "NA"
                if s["nonhuman_pct_of_classified"] is None
                else f"{s['nonhuman_pct_of_classified']:.3f}"
            ),
        ]
        if abs_enabled:
            cells.append(
                "NA"
                if s["nonhuman_abs_ppm_total"] is None
                else f"{s['nonhuman_abs_ppm_total']:,.2f}"
            )
        cells.extend([
            _fmt_int(s["richness"]),
            f"{s['shannon']:.3f}",
            f"{s['simpson']:.3f}",
        ])
        rows_html.append(
            "<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>"
        )

    summary_table = (
        "<h3>3.1 Per-sample summary</h3>"
        "<p>Every quantitative column labels its denominator.  "
        "<strong>% Non-human</strong> is compositional (within classified "
        "reads).  <strong>Non-human burden (ppm)</strong> is the "
        "absolute fraction of <em>total sequencing</em> attributed to "
        "non-human taxa.  The two metrics answer different scientific "
        "questions and must not be conflated.</p>"
        "<table><thead><tr>"
        + "".join(f"<th>{c}</th>" for c in header_cells)
        + "</tr></thead><tbody>" + "\n".join(rows_html) + "</tbody></table>"
    )

    # --- Cohort-wide top taxa ---
    cohort_totals: dict[int, float] = {}
    cohort_cpm_sum: dict[int, float] = {}
    cohort_abs_sum: dict[int, float] = {}
    for tid in inputs.matrix_raw.tax_ids:
        total = 0.0
        for sid in inputs.matrix_raw.sample_ids:
            v = inputs.matrix_raw.values[tid].get(sid)
            if v:
                total += v
        if total > 0:
            cohort_totals[tid] = total
            cohort_cpm_sum[tid] = sum(
                (inputs.matrix_cpm.values.get(tid, {}).get(sid) or 0.0)
                for sid in inputs.matrix_cpm.sample_ids
            )
            if inputs.matrix_abs is not None:
                cohort_abs_sum[tid] = sum(
                    (inputs.matrix_abs.values.get(tid, {}).get(sid) or 0.0)
                    for sid in inputs.matrix_abs.sample_ids
                )

    top_cohort = sorted(cohort_totals.items(), key=lambda kv: -kv[1])[:top_n]
    n_samples = len(inputs.matrix_raw.sample_ids) or 1
    cohort_header = [
        "Rank", "Taxon", "Tax ID", "Domain",
        "Cohort raw reads",
        "Mean CPM<span class='denom'>per million classified reads</span>",
    ]
    if abs_enabled:
        cohort_header.append(
            "Mean burden<span class='denom'>ppm of total sequenced reads</span>"
        )
    cohort_rows: list[str] = []
    for rank_i, (tid, total) in enumerate(top_cohort, start=1):
        cells = [
            str(rank_i),
            html.escape(inputs.matrix_raw.tax_names.get(tid, "")),
            str(tid),
            html.escape(inputs.matrix_raw.tax_domains.get(tid, "—")),
            _fmt_int(int(total)),
            f"{cohort_cpm_sum.get(tid, 0.0) / n_samples:,.2f}",
        ]
        if abs_enabled:
            cells.append(
                f"{cohort_abs_sum.get(tid, 0.0) / n_samples:,.2f}"
            )
        cohort_rows.append(
            "<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>"
        )

    cohort_table = (
        f"<h3>3.2 Top {len(top_cohort)} contaminating taxa (cohort-wide)</h3>"
        "<p>Ranked by total raw read count across the cohort.  Both CPM "
        "(composition within classified reads) and absolute burden "
        "(when available) are shown side-by-side with explicit "
        "denominators; use the absolute column when assessing "
        "downstream variant-calling impact.</p>"
        "<table><thead><tr>"
        + "".join(f"<th>{c}</th>" for c in cohort_header)
        + "</tr></thead><tbody>"
        + ("\n".join(cohort_rows) or
           "<tr><td colspan='7'><em>No non-zero taxa.</em></td></tr>")
        + "</tbody></table>"
    )

    # --- Stacked-bar composition figure (CPM) ---
    comp_cpm = [
        (sid, _domain_composition(inputs.matrix_cpm, sid))
        for sid in inputs.matrix_cpm.sample_ids
    ]
    fig_cpm = _stacked_bar_svg(
        comp_cpm,
        title="Figure 1 — Composition of classified reads by taxonomic domain",
        caption=(
            "Stacked bars are normalised to 100% of <strong>classified "
            "reads</strong> per sample (CPM denominator = sum of classified "
            "direct reads, pre-filter).  This figure describes "
            "<em>composition within the non-human fraction</em>; it does "
            "not indicate how much sequencing was contaminated.  For "
            "absolute burden see Figure 2."
        ),
        normalize=True,
    )

    fig_abs = ""
    if abs_enabled:
        comp_abs = [
            (sid, _domain_composition(inputs.matrix_abs, sid))  # type: ignore[arg-type]
            for sid in inputs.matrix_abs.sample_ids  # type: ignore[union-attr]
        ]
        fig_abs = _stacked_bar_svg(
            comp_abs,
            title=(
                "Figure 2 — Absolute non-human burden by taxonomic domain"
            ),
            caption=(
                "Bars share a common x-axis in reads per million "
                "<strong>total sequenced reads</strong> (denominator from "
                "<code>samtools idxstats</code>).  This is the metric to "
                "inspect when asking whether non-human content is "
                "extensive enough to impact variant calling, assembly, "
                "or expression quantitation."
            ),
            normalize=False,
        )

    detect_block = _render_detection_comparison(inputs)

    return summary_table + cohort_table + fig_cpm + fig_abs + detect_block


def _flagged_taxon_breakdown(
    flagged: list[dict[str, Any]], top_n: int = 5
) -> list[tuple[str, int]]:
    """Return the *top_n* most frequently flagged taxa.

    Each entry is ``(taxon_name, count)`` where ``count`` is the
    number of distinct ``(sample_id, taxon)`` flag rows for that taxon.
    """
    counts: dict[str, int] = {}
    for row in flagged:
        name = row.get("taxon_name") or row.get("tax_id") or "(unknown)"
        counts[name] = counts.get(name, 0) + 1
    return sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]


def _render_detection_comparison(inputs: ReportInputs) -> str:
    """Render §3.3 outlier-detection summary (CPM + abs side-by-side).

    Returns an empty string when no detection results are available.
    """
    primary = inputs.detect_summary
    abs_d = inputs.abs_detect_summary
    if primary is None and abs_d is None:
        return ""

    sections: list[str] = ["<h3>3.3 Outlier-detection results</h3>"]
    sections.append(
        "<p>The table below summarises every <code>csc-detect</code> "
        "pass that contributed to this report.  When both the primary "
        "(usually CPM, compositional) and the absolute-burden side "
        "passes are available, comparing their flagged-sample sets "
        "is highly informative: samples flagged by <em>both</em> are "
        "high-confidence contamination calls, while samples flagged "
        "by only one expose how the choice of denominator changes "
        "the answer.</p>"
    )

    rows: list[str] = []

    def _row(label: str, summary: dict[str, Any] | None) -> None:
        if summary is None:
            return
        flagged_samples = summary.get("flagged_samples") or []
        total_taxa = summary.get("total_taxa_analysed", "—")
        method = summary.get("method", "—")
        rows.append(
            "<tr>"
            f"<td>{html.escape(label)}</td>"
            f"<td><code>{html.escape(str(method))}</code></td>"
            f"<td>{_fmt_int(total_taxa) if isinstance(total_taxa, int) else html.escape(str(total_taxa))}</td>"
            f"<td>{_fmt_int(int(summary.get('flagged_count', 0)))}</td>"
            f"<td>{_fmt_int(len(flagged_samples))}</td>"
            "</tr>"
        )

    if primary is not None:
        primary_label = (
            "Primary "
            f"({primary.get('matrix_type', 'compositional')})"
        )
        _row(primary_label, primary)
    if abs_d is not None:
        _row("Absolute-burden side pass (abs)", abs_d)

    sections.append(
        "<table><thead><tr>"
        "<th>Pass</th>"
        "<th>Method</th>"
        "<th>Taxa analysed</th>"
        "<th>Flag rows<span class='denom'>(sample, taxon) pairs</span></th>"
        "<th>Distinct flagged samples</th>"
        "</tr></thead><tbody>"
        + "\n".join(rows)
        + "</tbody></table>"
    )

    # Set comparison – samples flagged by primary, by abs, by both.
    if primary is not None and abs_d is not None:
        primary_set = set(primary.get("flagged_samples") or [])
        abs_set = set(abs_d.get("flagged_samples") or [])
        both = primary_set & abs_set
        only_primary = primary_set - abs_set
        only_abs = abs_set - primary_set
        sections.append(
            "<table><thead><tr><th>Set</th><th>Count</th>"
            "<th>Sample IDs</th></tr></thead><tbody>"
            f"<tr><td>Flagged by both passes</td><td>{len(both)}</td>"
            f"<td><code>{html.escape(', '.join(sorted(both)) or '—')}</code></td></tr>"
            f"<tr><td>Flagged only by primary "
            f"({html.escape(str(primary.get('matrix_type', 'compositional')))})"
            f"</td><td>{len(only_primary)}</td>"
            f"<td><code>{html.escape(', '.join(sorted(only_primary)) or '—')}</code></td></tr>"
            f"<tr><td>Flagged only by absolute-burden pass</td>"
            f"<td>{len(only_abs)}</td>"
            f"<td><code>{html.escape(', '.join(sorted(only_abs)) or '—')}</code></td></tr>"
            "</tbody></table>"
            "<p class='fig-caption'>Samples appearing in both sets "
            "are robust contamination calls (independent of denominator). "
            "A sample flagged only by the absolute-burden pass typically "
            "indicates <em>more</em> total contaminating reads than the "
            "cohort but whose composition does not deviate &mdash; "
            "common when host-depletion efficiency varies.  A sample "
            "flagged only by the primary pass has unusual "
            "<em>composition</em> within its non-human reads.</p>"
        )

    # Top flagged taxa per pass (helps reviewers spot whether one pass
    # is dominated by a different organism than the other).
    primary_top = _flagged_taxon_breakdown(inputs.flagged_samples)
    abs_top = _flagged_taxon_breakdown(inputs.abs_flagged_samples)
    if primary_top or abs_top:

        def _fmt_top(items: list[tuple[str, int]]) -> str:
            if not items:
                return "—"
            return "; ".join(f"{html.escape(n)} ({c})" for n, c in items)

        sections.append(
            "<table><thead><tr><th>Pass</th>"
            "<th>Top flagged taxa<span class='denom'>"
            "(name → flag count)</span></th></tr></thead><tbody>"
            f"<tr><td>Primary</td><td>{_fmt_top(primary_top)}</td></tr>"
            f"<tr><td>Absolute-burden</td><td>{_fmt_top(abs_top)}</td></tr>"
            "</tbody></table>"
        )

    return "\n".join(sections)


def _render_variant_impact(
    per_sample_stats: list[dict[str, Any]],
    threshold_ppm: float,
    abs_enabled: bool,
) -> tuple[str, list[dict[str, Any]]]:
    """Render §4 and return the list of samples exceeding the threshold."""
    if not abs_enabled:
        body = (
            "<h2>4. Variant-Calling Impact</h2>"
            "<div class='callout'>This section requires the "
            "absolute-burden matrix (<code>taxa_matrix_abs.tsv</code>).  "
            "Re-run <code>csc-aggregate</code> with <code>--idxstats</code> "
            "and the per-sample <code>reads_summary.json</code> sidecars "
            "to enable.</div>"
        )
        return body, []

    flagged = [
        s for s in per_sample_stats
        if s["nonhuman_abs_ppm_total"] is not None
        and s["nonhuman_abs_ppm_total"] > threshold_ppm
    ]

    rows = []
    for s in sorted(
        flagged, key=lambda s: -(s["nonhuman_abs_ppm_total"] or 0)
    ):
        rows.append(
            "<tr>"
            f"<td>{html.escape(s['sample_id'])}</td>"
            f"<td>{_fmt_int(int(s['total_reads']) if s['total_reads'] else None)}</td>"
            f"<td>{s['nonhuman_abs_ppm_total']:,.2f}</td>"
            f"<td>{(s['nonhuman_abs_ppm_total'] or 0) / 1e4:.4f}%</td>"
            "</tr>"
        )
    table = (
        "<table><thead><tr><th>Sample</th>"
        "<th>Total sequenced reads</th>"
        "<th>Non-human burden (ppm)</th>"
        "<th>Non-human burden (% of sequencing)</th>"
        "</tr></thead><tbody>"
        + (
            "\n".join(rows)
            if rows
            else "<tr><td colspan='4'><em>No samples exceed the threshold.</em>"
                 "</td></tr>"
        )
        + "</tbody></table>"
    )

    callout_cls = "danger" if flagged else "ok"
    msg = (
        f"{len(flagged)} sample(s) exceed the configured threshold of "
        f"{threshold_ppm:,.0f} ppm of total sequenced reads "
        f"({threshold_ppm / 1e4:.3f}% of sequencing)."
    )
    body = (
        "<h2>4. Variant-Calling Impact</h2>"
        "<p>The <strong>absolute</strong> fraction of sequencing assigned "
        "to non-human taxa is the denominator reviewers typically "
        "scrutinise when judging whether contamination can affect "
        "variant calls, assemblies, or expression estimates.  The "
        "threshold below (configurable via "
        "<code>--variant-impact-threshold-ppm</code>) defaults to "
        f"{DEFAULT_VARIANT_IMPACT_THRESHOLD_PPM:,.0f} ppm "
        f"(= {DEFAULT_VARIANT_IMPACT_THRESHOLD_PPM / 1e4:.3f}% of total "
        "reads); set it to match the sensitivity of your downstream "
        "application and document the chosen value in methods.</p>"
        f"<div class='callout {callout_cls}'>{msg}</div>"
        + table
    )
    return body, flagged


def _render_discussion(
    inputs: ReportInputs,
    variant_flagged: list[dict[str, Any]],
) -> str:
    abs_enabled = inputs.matrix_abs is not None
    parts = [
        "<h2>5. Discussion &amp; Caveats</h2>",
        "<ul>",
        "<li>Kraken2 assignments are only as reliable as the reference "
        "database.  Using a redundant database (multiple genomes per "
        "species) inflates low-common-ancestor assignments and should be "
        "cross-checked against PrackenDB-style single-genome-per-species "
        "references.</li>",
        "<li>Compositional metrics (CPM, stacked bars, diversity "
        "indices) describe the <em>shape</em> of the non-human fraction "
        "but are invariant to the size of that fraction.  Two samples "
        "with identical stacked bars can have vastly different "
        "downstream impact.  Always pair CPM with absolute burden "
        "(§4) and, when available, with the absolute-burden detection "
        "side pass (§3.3).</li>",
        "<li>Laboratory and reagent contaminants ('kitome') can "
        "dominate low-input samples; consider excluding known kitome "
        "taxa in <code>csc-detect</code> via <code>--kitome-taxa</code>.</li>",
        "<li>Include negative controls and, where possible, known "
        "spike-ins to benchmark sensitivity and specificity before "
        "drawing biological conclusions from contaminant signal.</li>",
    ]
    if inputs.abs_detect_summary is not None and inputs.detect_summary is not None:
        primary_set = set(inputs.detect_summary.get("flagged_samples") or [])
        abs_set = set(inputs.abs_detect_summary.get("flagged_samples") or [])
        only_abs = abs_set - primary_set
        if only_abs:
            parts.append(
                f"<li>{len(only_abs)} sample(s) were flagged "
                "<strong>only</strong> by the absolute-burden side "
                "pass (§3.3).  Such samples typically have unusually "
                "high <em>total</em> non-human read counts but a "
                "composition that resembles the cohort &mdash; a "
                "pattern consistent with poor host depletion rather "
                "than a single-organism contaminant.  Investigate the "
                "library prep / extraction conditions.</li>"
            )
    if not abs_enabled:
        parts.append(
            "<li><strong>No absolute-burden denominator was available "
            "for this cohort.</strong>  Any claim about the magnitude of "
            "contamination (as opposed to its composition) should be "
            "deferred until <code>reads_summary.json</code> sidecars are "
            "supplied to <code>csc-aggregate</code>.  This also "
            "disables the absolute-burden detection side pass; "
            "outlier flags reported here are compositional only.</li>"
        )
    if variant_flagged:
        parts.append(
            f"<li>{len(variant_flagged)} sample(s) were flagged in §4 as "
            "potentially impacting variant calling.  These should be "
            "inspected manually before downstream use; consider "
            "re-extracting with a stricter MAPQ filter or excluding them "
            "from analyses where non-human content could bias results.</li>"
        )
    parts.append("</ul>")
    return "\n".join(parts)


def _render_transparency_checklist(
    inputs: ReportInputs,
    args_echo: dict[str, Any],
) -> str:
    agg_meta = inputs.aggregate_metadata
    aggregate_dir = inputs.aggregate_dir or Path(".")
    downloads: list[tuple[str, str]] = []

    # Primary matrices
    downloads.append(
        (
            f"{aggregate_dir / 'taxa_matrix_raw.tsv'}",
            "Raw taxa-by-sample matrix (integer direct reads)",
        )
    )
    downloads.append(
        (
            f"{aggregate_dir / 'taxa_matrix_cpm.tsv'}",
            "CPM matrix (per million classified reads – compositional)",
        )
    )
    if inputs.matrix_abs is not None:
        downloads.append(
            (
                f"{aggregate_dir / 'taxa_matrix_abs.tsv'}",
                "Absolute-burden matrix (per million total sequenced reads)",
            )
        )
    downloads.append(
        (
            f"{aggregate_dir / 'aggregation_metadata.json'}",
            "Aggregation metadata (schema version, parameters, provenance)",
        )
    )
    rank_meta = aggregate_dir / "rank_filter_metadata.json"
    if rank_meta.exists():
        downloads.append(
            (str(rank_meta), "Per-rank matrix metadata")
        )
    if inputs.detect_dir is not None:
        for name, label in [
            ("flagged_samples.tsv", "Detect: flagged (sample, taxon) pairs"),
            ("qc_summary.json", "Detect: QC summary"),
            ("quarantine_list.txt", "Detect: quarantine list"),
        ]:
            p = inputs.detect_dir / name
            if p.exists():
                downloads.append((str(p), label))
    if inputs.abs_detect_dir is not None:
        for name, label in [
            (
                "flagged_samples.tsv",
                "Detect (abs): flagged pairs from absolute-burden pass",
            ),
            ("qc_summary.json", "Detect (abs): QC summary"),
            ("quarantine_list.txt", "Detect (abs): quarantine list"),
        ]:
            p = inputs.abs_detect_dir / name
            if p.exists():
                downloads.append((str(p), label))

    # Rank matrices
    for rank in agg_meta.get("rank_filter", []):
        for tag, label in [
            ("raw", "raw"), ("cpm", "cpm"), ("abs", "abs"),
        ]:
            p = aggregate_dir / f"taxa_matrix_{tag}_{rank}.tsv"
            if p.exists():
                downloads.append(
                    (str(p), f"Rank-filtered {label.upper()} matrix ({rank})")
                )

    dl_html = "\n".join(
        f'<li><a href="{html.escape(str(Path(p).name))}">'
        f'{html.escape(str(Path(p).name))}</a> — {html.escape(label)}</li>'
        for p, label in downloads
    )

    checklist_items = [
        f"CSC aggregate schema version: <code>{html.escape(str(agg_meta.get('schema_version', 'unknown')))}</code>",
        f"Report schema version: <code>{REPORT_SCHEMA_VERSION}</code>",
        f"Minimum reads per taxon filter: <code>{html.escape(str(agg_meta.get('min_reads', 'unknown')))}</code>",
        f"Rank filter: <code>{html.escape(', '.join(agg_meta.get('rank_filter', [])) or 'none')}</code>",
        f"Lineage-aware domain annotation: <code>{'enabled' if agg_meta.get('domain_annotated') else 'disabled'}</code>",
        f"Absolute-burden matrix: <code>{'present' if inputs.matrix_abs is not None else 'absent'}</code>",
        f"Absolute-burden detection side pass: <code>{'present' if inputs.abs_detect_summary is not None else 'absent'}</code>",
        f"Variant-calling impact threshold (ppm): <code>{args_echo.get('threshold_ppm')}</code>",
        f"Samples missing idxstats sidecars: <code>{len(agg_meta.get('samples_without_idxstats', []))}</code>",
        "Every figure caption explicitly names its denominator.",
        "Every table column labels CPM vs absolute burden explicitly.",
        "All matrices are linked below for byte-identical reproduction.",
    ]

    return (
        "<h2>6. Methods Transparency Checklist</h2>"
        "<ul class='checklist'>"
        + "\n".join(f"<li>{item}</li>" for item in checklist_items)
        + "</ul>"
        "<h3>6.1 Downloadable artefacts</h3>"
        "<p>The following files are referenced by this report.  Links "
        "are relative to the report's output directory (copy the files "
        "alongside the HTML to preserve the links).</p>"
        f"<ul class='download-list'>{dl_html}</ul>"
    )


# ---------------------------------------------------------------------------
# Top-level renderer
# ---------------------------------------------------------------------------


def generate_html_report(
    inputs: ReportInputs,
    output_path: str | Path,
    *,
    top_n: int = 10,
    threshold_ppm: float = DEFAULT_VARIANT_IMPACT_THRESHOLD_PPM,
    title: str = "Cross-Species Contamination Summary Report",
) -> Path:
    """Render the static HTML report.

    Parameters
    ----------
    inputs:
        Parsed pipeline outputs, typically from :func:`load_inputs`.
    output_path:
        Destination HTML file.
    top_n:
        Number of taxa to include in the cohort-wide summary table
        (default 10).
    threshold_ppm:
        Absolute-burden threshold, in reads per million total sequenced
        reads, above which a sample is flagged in §4 Variant-Calling
        Impact.  Defaults to
        :data:`DEFAULT_VARIANT_IMPACT_THRESHOLD_PPM` (= 1,000 ppm,
        i.e. 0.1% of total sequencing).  Document your chosen value in
        the manuscript methods.

    Returns
    -------
    pathlib.Path
        Path of the written HTML file.  A sibling
        ``report_manifest.json`` is also written so downstream tooling
        can introspect the inputs without re-parsing the HTML.
    """
    if top_n < 1:
        raise ValueError("top_n must be >= 1")
    if threshold_ppm < 0:
        raise ValueError("threshold_ppm must be >= 0")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    per_sample = _per_sample_stats(inputs)
    variant_section, variant_flagged = _render_variant_impact(
        per_sample,
        threshold_ppm=threshold_ppm,
        abs_enabled=inputs.matrix_abs is not None,
    )
    exec_section = _render_executive_summary(
        inputs, per_sample, variant_flagged, threshold_ppm
    )
    methods_section = _render_methods(inputs)
    results_section = "<h2>3. Results</h2>" + _render_results(
        inputs, per_sample, top_n=top_n
    )
    discussion_section = _render_discussion(inputs, variant_flagged)
    checklist_section = _render_transparency_checklist(
        inputs, {"threshold_ppm": threshold_ppm, "top_n": top_n}
    )

    generated_at = datetime.now(tz=timezone.utc).isoformat(
        timespec="seconds"
    )
    footer = (
        f"<p class='footer'>Generated by <code>csc-report</code> "
        f"(schema {REPORT_SCHEMA_VERSION}) at {generated_at}.  "
        f"This report is static and self-contained – no external "
        f"resources are required to view it.</p>"
    )

    body = (
        f"<h1>{html.escape(title)}</h1>"
        + exec_section
        + methods_section
        + results_section
        + variant_section
        + discussion_section
        + checklist_section
        + footer
    )

    doc = (
        "<!DOCTYPE html>\n"
        "<html lang=\"en\">\n<head>\n"
        "<meta charset=\"utf-8\" />\n"
        f"<title>{html.escape(title)}</title>\n"
        f"<style>{_PAGE_STYLE}</style>\n"
        "</head>\n<body>\n"
        + body
        + "\n</body>\n</html>\n"
    )

    output_path.write_text(doc, encoding="utf-8")
    logger.info("Wrote HTML report to %s", output_path)

    # Sidecar manifest – useful for machine-readable introspection and
    # golden-file diffing.
    manifest = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "generated_at": generated_at,
        "title": title,
        "variant_impact_threshold_ppm": threshold_ppm,
        "top_n": top_n,
        "sample_count": len(inputs.matrix_raw.sample_ids),
        "taxon_count": len(inputs.matrix_raw.tax_ids),
        "absolute_burden_enabled": inputs.matrix_abs is not None,
        "abs_detection_enabled": inputs.abs_detect_summary is not None,
        "samples_flagged_variant_impact": [
            s["sample_id"] for s in variant_flagged
        ],
        "samples_flagged_primary_detect": (
            list(inputs.detect_summary.get("flagged_samples", []))
            if inputs.detect_summary is not None
            else []
        ),
        "samples_flagged_abs_detect": (
            list(inputs.abs_detect_summary.get("flagged_samples", []))
            if inputs.abs_detect_summary is not None
            else []
        ),
        "primary_matrix_type": (
            inputs.detect_summary.get("matrix_type")
            if inputs.detect_summary is not None
            else None
        ),
        "aggregate_metadata_schema_version": inputs.aggregate_metadata.get(
            "schema_version"
        ),
    }
    manifest_path = output_path.with_name("report_manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)
    logger.info("Wrote report manifest to %s", manifest_path)

    return output_path
