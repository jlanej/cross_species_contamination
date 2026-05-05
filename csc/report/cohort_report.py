"""Renderers for the cohort-oriented HTML report layout (§1–§8).

Sectioning:

* §1 Executive summary (cohort)  – :func:`render_executive_summary_v2`
* §2 Methods (unchanged)          – :func:`csc.report.report._render_methods`
* §3 Cohort species landscape     – :func:`render_cohort_landscape`
* §4 Variant-calling impact       – :func:`render_variant_impact_v2`
* §5 Detection summary            – :func:`render_detection_section_v2`
* §6 Per-sample appendix          – :func:`render_per_sample_appendix`
* §7 Discussion & caveats         – :func:`render_discussion_v2`
* §8 Methods transparency         – :func:`csc.report.report._render_transparency_checklist`

This module focuses on **species** as the unit of analysis.  The
older per-sample renderers in :mod:`csc.report.report` remain
available behind ``--layout legacy`` for one release window.
"""

from __future__ import annotations

import csv
import html
import json
import logging
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from csc.report import cohort as _cohort
from csc.report import svg as _svg

logger = logging.getLogger(__name__)


# Maximum number of per-sample rows to emit in the §5.1 concordance
# "Reads demoted to unclassified" table before truncating with a
# follow-up note pointing at aggregation_metadata.json.  Keeps very
# large cohorts (3K+ samples) from blowing up the report's HTML size.
_DEMOTED_TABLE_MAX_ROWS = 200


def _esc(s: object) -> str:
    return html.escape(str(s))


def _fmt_num(v: float | None, *, ndigits: int = 2) -> str:
    if v is None:
        return "NA"
    if isinstance(v, float) and math.isnan(v):
        return "NA"
    if v == 0:
        return "0"
    av = abs(v)
    if av >= 1_000_000:
        return f"{v:,.0f}"
    if av >= 1:
        return f"{v:,.{ndigits}f}"
    return f"{v:.{max(ndigits, 4)}g}"


def _fmt_int(v: int | None) -> str:
    if v is None:
        return "NA"
    return f"{v:,}"


def _fmt_pct(v: float | None, *, ndigits: int = 1) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "NA"
    return f"{v * 100:.{ndigits}f}%"


# ---------------------------------------------------------------------------
# §1 Executive summary
# ---------------------------------------------------------------------------


def render_executive_summary_v2(
    inputs: Any,
    per_sample_stats: Sequence[Mapping[str, Any]],
    species_rows: Sequence[Mapping[str, Any]],
    partition: Mapping[str, Any],
    variant_flagged: Sequence[Mapping[str, Any]],
    threshold_ppm: float,
) -> str:
    """Cohort-headline executive summary.

    Mirrors the field labels of the legacy renderer so that downstream
    consumers (and tests) can still find ``"per million classified
    reads"`` etc.
    """
    meta = inputs.aggregate_metadata
    n_samples = meta.get("sample_count", len(inputs.matrix_raw.sample_ids))
    n_taxa = meta.get("taxon_count", len(inputs.matrix_raw.tax_ids))
    abs_enabled = inputs.matrix_abs is not None

    # Cohort-level numbers reused from per-sample stats.
    total_classified = sum(s["classified_reads"] for s in per_sample_stats)

    pct_classified = [
        s["nonhuman_pct_of_classified"]
        for s in per_sample_stats
        if s["nonhuman_pct_of_classified"] is not None
    ]
    burden_ppm = [
        s["nonhuman_abs_ppm_total"]
        for s in per_sample_stats
        if s["nonhuman_abs_ppm_total"] is not None
    ]
    median_pct = (sorted(pct_classified)[len(pct_classified) // 2]
                  if pct_classified else None)
    median_burden = (sorted(burden_ppm)[len(burden_ppm) // 2]
                     if burden_ppm else None)

    # Species-centric headlines.
    top_burden = species_rows[0] if species_rows else None
    top_prev = max(species_rows, key=lambda r: r["prevalence"]) if species_rows else None
    domain_mix = _cohort.domain_burden_mix(
        species_rows,
        metric="cohort_burden_ppm" if abs_enabled else "cohort_raw_total",
    )

    items: list[str] = []
    items.append(f"<li><span class='meta-key'>Samples analysed:</span> {n_samples}</li>")
    items.append(f"<li><span class='meta-key'>Distinct taxa observed:</span> {n_taxa}</li>")
    items.append(
        f"<li><span class='meta-key'>Total classified reads (cohort):</span> "
        f"{_fmt_int(int(total_classified))}</li>"
    )
    if median_pct is not None:
        items.append(
            f"<li><span class='meta-key'>Median % non-human "
            f"<em>among classified reads</em> (per million classified reads "
            f"denominator):</span> {median_pct:.3f}%</li>"
        )
    if median_burden is not None:
        items.append(
            f"<li><span class='meta-key'>Median non-human absolute burden "
            f"<em>per million total sequenced reads (ppm of total sequenced "
            f"reads)</em>:</span> {median_burden:,.2f} ppm "
            f"({median_burden / 1e4:.4f}%)</li>"
        )
    items.append(
        f"<li><span class='meta-key'>Species partition by prevalence "
        f"(thresholds {partition['core_threshold']:.0%} / "
        f"{partition['rare_threshold']:.0%}):</span> "
        f"<strong>core</strong> {partition['core_count']}, "
        f"<strong>accessory</strong> {partition['accessory_count']}, "
        f"<strong>rare</strong> {partition['rare_count']}</li>"
    )
    if top_burden is not None:
        burden_v = top_burden.get("cohort_burden_ppm")
        burden_str = (
            f"{burden_v:,.1f} ppm cohort burden"
            if burden_v is not None and not (isinstance(burden_v, float) and math.isnan(burden_v))
            else f"{int(top_burden['cohort_raw_total']):,} cohort raw reads"
        )
        items.append(
            f"<li><span class='meta-key'>Top species by cohort burden:</span> "
            f"<em>{_esc(top_burden['name'])}</em> "
            f"(tax_id {top_burden['tax_id']}, {_esc(top_burden['domain'])}; "
            f"{burden_str})</li>"
        )
    if top_prev is not None:
        items.append(
            f"<li><span class='meta-key'>Most prevalent species:</span> "
            f"<em>{_esc(top_prev['name'])}</em> "
            f"({_fmt_pct(top_prev['prevalence'])} of samples)</li>"
        )
    if inputs.detect_summary is not None:
        n_flag = len(inputs.detect_summary.get("flagged_samples", []))
        items.append(
            f"<li><span class='meta-key'>Samples flagged by "
            f"<code>csc-detect</code> (primary pass):</span> {n_flag}</li>"
        )
    if inputs.abs_detect_summary is not None:
        abs_flag = len(inputs.abs_detect_summary.get("flagged_samples", []))
        items.append(
            f"<li><span class='meta-key'>Samples flagged by "
            f"<code>csc-detect</code> (absolute-burden side pass):</span> "
            f"{abs_flag}</li>"
        )
    items.append(
        f"<li><span class='meta-key'>Samples exceeding variant-calling "
        f"impact threshold (&gt; {threshold_ppm:,.0f} ppm of total reads, "
        f"i.e. {threshold_ppm / 1e4:.3f}% of sequencing):</span> "
        f"{len(variant_flagged) if abs_enabled else 'N/A (no idxstats)'}</li>"
    )

    callout = ""
    if not abs_enabled:
        callout = (
            "<div class='callout'>⚠ <strong>Absolute-burden matrix "
            "(<code>taxa_matrix_abs.tsv</code>) not available.</strong> "
            "Re-run <code>csc-aggregate --idxstats …</code> with the "
            "per-sample <code>reads_summary.json</code> sidecars to enable "
            "absolute-burden reporting and the variant-calling impact "
            "section.  Without it this report can only describe "
            "<em>composition</em> – not the absolute fraction of "
            "sequencing affected.</div>"
        )

    mini_bar = ""
    if domain_mix:
        mini_bar = (
            "<p class='meta-key'>Cohort burden by taxonomic domain "
            "(metric: " +
            ("absolute ppm" if abs_enabled else "cohort raw reads") +
            "):</p>"
            + _svg.mini_stacked_bar_svg(domain_mix)
            + _svg.domain_legend_html(sorted(domain_mix.keys()))
        )

    return (
        "<h2>1. Executive Summary</h2>"
        "<p>This report summarises non-human content detected in a cohort of "
        "human whole-genome sequencing (WGS) samples.  It treats "
        "<strong>species</strong> as the unit of analysis: see §3 for the "
        "cohort-wide species landscape, §4 for variant-calling impact, "
        "§5 for outlier-detection summary, and §6 for the (paginated) "
        "per-sample appendix.  All tables and figures label their "
        "denominators explicitly (composition within contaminants vs. "
        "fraction of total sequencing).</p>"
        + callout
        + "<ul>" + "\n".join(items) + "</ul>"
        + mini_bar
    )


# ---------------------------------------------------------------------------
# §3 Cohort species landscape
# ---------------------------------------------------------------------------


def render_cohort_landscape(
    inputs: Any,
    species_rows: Sequence[Mapping[str, Any]],
    partition: Mapping[str, Any],
    *,
    page_size: int = 25,
    cluster_method: str = "average",
    cluster_distance: str = "bray",
    max_samples_cluster: int = 2000,
    max_samples_heatmap: int | None = None,
    top_species_heatmap: int = 50,
    drilldown_top: int = 25,
    species_table_top: int = 200,
    species_table_tsv: Path | None = None,
) -> str:
    """Render the new §3 — *Cohort species landscape*.

    Subsections:

    * §3.1 Species summary table (paginated, sortable, filterable).
    * §3.2 Prevalence-vs-abundance map ("contamination quadrants").
    * §3.3 Rank-abundance / Whittaker curve.
    * §3.4 Core / accessory / rare partition.
    * §3.5 Cohort-wide composition figures (boxplots, mini-stacked bar).
    * §3.6 Sample × species heatmap with hierarchical-cluster ordering.
    * §3.7 Sample β-diversity ordination (PCoA).
    * §3.8 Per-species drill-down for the top-N species.

    The §3.1 species summary table is capped at ``species_table_top``
    rows in the HTML (sorted by cohort burden); when truncated the
    full table is written to ``species_table_tsv`` so reviewers can
    audit every taxon offline.  The heatmap is capped independently
    via ``max_samples_heatmap`` (defaults to ``max_samples_cluster``)
    because every cell adds inline-SVG markup.
    """
    abs_enabled = inputs.matrix_abs is not None
    n_species = len(species_rows)
    n_samples = len(inputs.matrix_raw.sample_ids)
    parts: list[str] = ["<h2>3. Cohort species landscape</h2>"]

    if max_samples_heatmap is None:
        max_samples_heatmap = max_samples_cluster

    logger.info("§3.1 species summary table (%d species)", n_species)
    parts.append(_render_species_table(species_rows, page_size=page_size,
                                        abs_enabled=abs_enabled,
                                        top_n=species_table_top,
                                        tsv_path=species_table_tsv))
    logger.info("§3.2 prevalence–abundance map")
    parts.append(_render_prevalence_abundance_map(species_rows,
                                                   abs_enabled=abs_enabled))
    logger.info("§3.3 rank-abundance curve")
    parts.append(_render_rank_abundance(species_rows, abs_enabled=abs_enabled))
    logger.info("§3.4 core/accessory/rare partition")
    parts.append(_render_partition(partition, abs_enabled=abs_enabled))
    logger.info("§3.5 cohort distribution figures")
    parts.append(_render_distribution_figures(inputs, species_rows))
    logger.info(
        "§3.6 heatmap (top %d species, %d samples cap, method=%s, dist=%s)",
        top_species_heatmap, max_samples_heatmap, cluster_method, cluster_distance,
    )
    parts.append(_render_heatmap(
        inputs,
        species_rows,
        cluster_method=cluster_method,
        cluster_distance=cluster_distance,
        max_samples_cluster=max_samples_heatmap,
        top_species=top_species_heatmap,
    ))
    logger.info(
        "§3.7 PCoA β-diversity ordination (%d samples, dist=%s)",
        n_samples, cluster_distance,
    )
    parts.append(_render_pcoa(
        inputs,
        species_rows,
        cluster_distance=cluster_distance,
        max_samples=max_samples_cluster,
        top_species=top_species_heatmap,
    ))
    logger.info("§3.8 per-species drill-down (top %d)", drilldown_top)
    parts.append(_render_species_drilldown(
        inputs, species_rows, top_n=drilldown_top, abs_enabled=abs_enabled,
    ))
    logger.info("§3 cohort landscape complete (§3.1–§3.8)")

    return "\n".join(parts)


def _render_species_table(
    species_rows: Sequence[Mapping[str, Any]],
    *,
    page_size: int,
    abs_enabled: bool,
    top_n: int = 200,
    tsv_path: Path | None = None,
) -> str:
    n_total = len(species_rows)
    truncated = n_total > top_n
    display_rows = species_rows[:top_n] if truncated else species_rows

    # Always write the full species summary as a TSV sidecar so
    # reviewers can audit every taxon offline even when the HTML
    # table is truncated for renderability.
    sidecar_link = ""
    if tsv_path is not None:
        try:
            _write_species_summary_tsv(species_rows, tsv_path,
                                       abs_enabled=abs_enabled)
            sidecar_link = (
                f"  The full table ({_fmt_int(n_total)} taxa) is also "
                f"written to <a class='dl-link' href='{_esc(tsv_path.name)}'>"
                f"{_esc(tsv_path.name)}</a> for offline analysis."
            )
        except OSError as exc:  # pragma: no cover - defensive
            logger.warning("Could not write species summary TSV: %s", exc)

    parts: list[str] = ["<h3>3.1 Species summary table</h3>"]
    truncation_note = ""
    if truncated:
        truncation_note = (
            f"<div class='callout'>For renderability the table below "
            f"shows the <strong>top {top_n}</strong> non-human taxa "
            f"sorted by cohort burden; "
            f"<strong>{_fmt_int(n_total - top_n)} additional taxa</strong> "
            f"are omitted from the HTML and only present in the sidecar "
            f"TSV.</div>"
        )
    parts.append(
        "<p>One row per non-human taxon (excluding tax_id 9606).  Robust "
        "summaries: <strong>median (positives only)</strong>, <strong>IQR</strong> "
        "and <strong>95th percentile</strong> over samples that detected the "
        "taxon.  Compositional values are <em>per million classified reads</em> "
        "(CPM); absolute values are <em>ppm of total sequenced reads</em>.  "
        "<code>CV (robust)</code> is MAD&nbsp;/&nbsp;median over the same "
        "positives — high values flag ‘outburst’ species, low values flag "
        "baseline / kitome-like species.  Click any column header to sort; "
        "use the search box and the domain dropdown to filter; bottom controls "
        f"paginate at {page_size} rows per page." + sidecar_link + "</p>"
    )
    parts.append(truncation_note)

    headers = [
        ("Taxon", None, None),
        ("Tax ID", None, None),
        ("Domain", "domain", None),
        ("Prevalence<span class='denom'>any read</span>", None, None),
        ("Prevalence<span class='denom'>≥ min reads</span>", None, None),
        ("Median CPM<span class='denom'>positives</span>", None, None),
        ("IQR CPM<span class='denom'>positives</span>", None, None),
        ("P95 CPM<span class='denom'>positives</span>", None, None),
    ]
    if abs_enabled:
        headers.extend([
            ("Median ppm<span class='denom'>positives, ppm of total sequenced reads</span>",
             None, None),
            ("IQR ppm<span class='denom'>positives</span>", None, None),
            ("P95 ppm<span class='denom'>positives</span>", None, None),
            ("Cohort burden<span class='denom'>Σ ppm of total sequenced reads</span>",
             None, None),
        ])
    headers.extend([
        ("Cohort raw reads", None, None),
        ("Flagged (primary)", None, None),
        ("Flagged (abs)", None, None),
        ("CV (robust)<span class='denom'>MAD/median, CPM positives</span>", None, None),
        ("Distribution<span class='denom'>log-CPM positives</span>", None, "nosort"),
    ])

    head_html = []
    for label, ftype, opt in headers:
        attrs = []
        if ftype:
            attrs.append(f'data-filter="{ftype}"')
        if opt == "nosort":
            attrs.append("data-nosort")
        head_html.append(f"<th {' '.join(attrs)}>{label}</th>")

    # Default sort: cohort burden desc when abs_enabled (column index 11),
    # else cohort raw reads desc (column index 8).
    default_sort_col = 11 if abs_enabled else 8
    rows_html: list[str] = []
    for r in display_rows:
        cells = [
            f"<td>{_esc(r['name'])}</td>",
            f"<td>{r['tax_id']}</td>",
            f"<td>{_esc(r['domain'])}</td>",
            f"<td data-sort='{r['prevalence']}'>{_fmt_pct(r['prevalence'])}</td>",
            f"<td data-sort='{r['prevalence_at_min']}'>"
            f"{_fmt_pct(r['prevalence_at_min'])}</td>",
            f"<td>{_fmt_num(r['cpm_median_pos'])}</td>",
            f"<td>{_fmt_num(r['cpm_iqr_pos'])}</td>",
            f"<td>{_fmt_num(r['cpm_p95_pos'])}</td>",
        ]
        if abs_enabled:
            cells.extend([
                f"<td>{_fmt_num(r['abs_median_pos'])}</td>",
                f"<td>{_fmt_num(r['abs_iqr_pos'])}</td>",
                f"<td>{_fmt_num(r['abs_p95_pos'])}</td>",
                f"<td>{_fmt_num(r['cohort_burden_ppm'])}</td>",
            ])
        cells.extend([
            f"<td>{_fmt_int(int(r['cohort_raw_total']))}</td>",
            f"<td>{_fmt_int(int(r['n_flagged_primary']))}</td>",
            f"<td>{_fmt_int(int(r['n_flagged_abs']))}</td>",
            f"<td>{_fmt_num(r['cv_robust'])}</td>",
            f"<td>{_svg.sparkline_svg(r['sparkline'])}</td>",
        ])
        rows_html.append("<tr>" + "".join(cells) + "</tr>")

    if not rows_html:
        rows_html.append(
            f"<tr><td colspan='{len(headers)}'><em>No non-human taxa observed.</em>"
            "</td></tr>"
        )

    parts.append(
        f'<table class="paginated" data-page-size="{page_size}" '
        f'data-default-sort="{default_sort_col}|desc">'
        f'<thead><tr>{"".join(head_html)}</tr></thead>'
        f'<tbody>{chr(10).join(rows_html)}</tbody></table>'
    )
    return "\n".join(parts)


def _write_species_summary_tsv(
    species_rows: Sequence[Mapping[str, Any]],
    tsv_path: Path,
    *,
    abs_enabled: bool,
) -> None:
    """Write the full §3.1 species summary as a TSV sidecar."""
    header = [
        "tax_id", "name", "domain", "prevalence", "prevalence_at_min",
        "positives_n", "cpm_median_pos", "cpm_iqr_pos", "cpm_p95_pos",
    ]
    if abs_enabled:
        header += [
            "abs_median_pos", "abs_iqr_pos", "abs_p95_pos",
            "cohort_burden_ppm",
        ]
    header += [
        "cohort_raw_total", "n_flagged_primary", "n_flagged_abs", "cv_robust",
    ]

    def _f(v: Any) -> str:
        if v is None:
            return "NA"
        if isinstance(v, float) and math.isnan(v):
            return "NA"
        if isinstance(v, float):
            return f"{v:.6g}"
        return str(v)

    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for r in species_rows:
            row = [
                str(r["tax_id"]),
                str(r["name"]),
                str(r["domain"]),
                _f(r.get("prevalence")),
                _f(r.get("prevalence_at_min")),
                _f(r.get("positives_n")),
                _f(r.get("cpm_median_pos")),
                _f(r.get("cpm_iqr_pos")),
                _f(r.get("cpm_p95_pos")),
            ]
            if abs_enabled:
                row += [
                    _f(r.get("abs_median_pos")),
                    _f(r.get("abs_iqr_pos")),
                    _f(r.get("abs_p95_pos")),
                    _f(r.get("cohort_burden_ppm")),
                ]
            row += [
                _f(int(r["cohort_raw_total"])),
                _f(int(r["n_flagged_primary"])),
                _f(int(r["n_flagged_abs"])),
                _f(r.get("cv_robust")),
            ]
            w.writerow(row)


def _render_prevalence_abundance_map(
    species_rows: Sequence[Mapping[str, Any]],
    *,
    abs_enabled: bool,
) -> str:
    parts = ["<h3>3.2 Prevalence-vs-abundance map (contamination quadrants)</h3>"]
    parts.append(
        "<p>One dot per species.  <em>x</em> = prevalence (fraction of "
        "samples with ≥1 read).  <em>y</em> = log10 median absolute burden "
        "(ppm of total sequenced reads, positives only) when available, "
        "otherwise log10 median CPM (per million classified reads).  Dot "
        "size is proportional to log10 cohort raw reads, colour by domain.  "
        "Quadrant guidance (50% prevalence × the median y-value of the "
        "scatter) separates <em>core kitome</em>, <em>systemic burden</em>, "
        "<em>outburst</em>, and <em>trace</em> species.</p>"
    )

    metric = "abs_median_pos" if abs_enabled else "cpm_median_pos"
    pts = []
    domains = sorted({r["domain"] for r in species_rows})
    cmap = _svg.domain_colour_map(domains)
    ymetric_vals: list[float] = []
    for r in species_rows:
        v = r.get(metric)
        if v is None or (isinstance(v, float) and math.isnan(v)) or v <= 0:
            continue
        if r["prevalence"] <= 0:
            continue
        ymetric_vals.append(float(v))
        cohort_raw = r.get("cohort_raw_total") or 1
        size = max(2.5, min(10.0, 2.0 + math.log10(max(cohort_raw, 1)) * 0.9))
        pts.append({
            "x": r["prevalence"],
            "y": float(v),
            "size": size,
            "domain": r["domain"],
            "colour": cmap.get(r["domain"]),
            "label": (
                f"{r['name']} (tax_id {r['tax_id']}) — "
                f"prev {r['prevalence'] * 100:.1f}%, "
                f"median {('ppm' if abs_enabled else 'CPM')} "
                f"{float(v):.3g}"
            ),
        })

    qy = math.log10(sorted(ymetric_vals)[len(ymetric_vals) // 2]) if ymetric_vals else 0.0
    quadrant_labels = (
        "core kitome",            # high prev, low burden  (top-left in plot, but x is prev so prev high right)
        "systemic burden",        # high prev, high burden
        "trace",                  # low prev, low burden
        "outburst",               # low prev, high burden — hmm orientation
    )
    # plot uses (top-left, top-right, bottom-left, bottom-right). After
    # log y, "high y" = top.
    quadrant_labels = (
        "outburst (low prev, high burden)",
        "systemic (high prev, high burden)",
        "trace (low prev, low burden)",
        "core kitome (high prev, low burden)",
    )

    fig = _svg.scatter_svg(
        pts,
        title=("Figure 3.2 — prevalence vs absolute burden"
               if abs_enabled else
               "Figure 3.2 — prevalence vs compositional CPM"),
        x_label="prevalence (fraction of samples)",
        y_label=("log10 median ppm of total sequenced reads (positives)"
                 if abs_enabled else
                 "log10 median CPM, per million classified reads (positives)"),
        x_log=False,
        y_log=True,
        quadrant_lines=(0.5, qy),
        quadrant_labels=quadrant_labels,
        legend_domains=domains,
    )
    parts.append(fig)
    parts.append(
        "<p class='fig-caption'>The denominator behind the y-axis is "
        "stated in its label.  Always inspect this map together with §4 — "
        "high prevalence at low absolute burden typically reflects "
        "reagent / library contamination ('kitome') rather than "
        "biological signal.</p>"
    )
    return "\n".join(parts)


def _render_rank_abundance(
    species_rows: Sequence[Mapping[str, Any]],
    *,
    abs_enabled: bool,
) -> str:
    metric = "cohort_burden_ppm" if abs_enabled else "cohort_raw_total"
    series = _cohort.rank_abundance(species_rows, metric=metric)
    parts = [
        "<h3>3.3 Rank-abundance (Whittaker) curve</h3>",
        "<p>Log-rank vs log-cohort-burden.  A steep slope indicates strong "
        "cohort-level dominance by a few species; a shallow slope indicates "
        "broad, even contamination consistent with environmental / reagent "
        "background.</p>",
        _svg.rank_abundance_svg(
            series,
            title=("Figure 3.3 — Rank-abundance (cohort burden, ppm of total "
                   "sequenced reads)"
                   if abs_enabled else
                   "Figure 3.3 — Rank-abundance (cohort raw reads)"),
        ),
    ]
    return "\n".join(parts)


def _render_partition(
    partition: Mapping[str, Any], *, abs_enabled: bool,
) -> str:
    def _band_table(title: str, rows: Sequence[Mapping[str, Any]]) -> str:
        if not rows:
            return (
                f"<div><h4>{_esc(title)} (0)</h4>"
                "<p><em>No species in this band.</em></p></div>"
            )
        body = []
        for r in rows:
            metric = r.get("cohort_burden_ppm")
            metric_str = (
                f"{metric:,.1f} ppm"
                if metric is not None and not (isinstance(metric, float) and math.isnan(metric))
                else f"{int(r['cohort_raw_total']):,} reads"
            )
            body.append(
                f"<tr><td>{_esc(r['name'])}</td>"
                f"<td>{_esc(r['domain'])}</td>"
                f"<td>{_fmt_pct(r['prevalence'])}</td>"
                f"<td>{metric_str}</td></tr>"
            )
        return (
            f"<div><h4>{_esc(title)} ({len(rows)})</h4>"
            f"<table><thead><tr><th>Taxon</th><th>Domain</th>"
            f"<th>Prev.</th><th>Burden</th></tr></thead>"
            f"<tbody>{''.join(body)}</tbody></table></div>"
        )

    return (
        "<h3>3.4 Core / accessory / rare partition</h3>"
        "<p>Species are partitioned by prevalence into three bands: "
        f"<strong>core</strong> (≥{partition['core_threshold'] * 100:.0f}% of "
        f"samples), <strong>accessory</strong> "
        f"({partition['rare_threshold'] * 100:.0f}–"
        f"{partition['core_threshold'] * 100:.0f}%), "
        f"<strong>rare</strong> (&lt;{partition['rare_threshold'] * 100:.0f}%).  "
        "Within each band the top species are ranked by cohort burden "
        "(ppm of total sequenced reads when available, else cohort raw "
        "reads).</p>"
        "<div class='partition-grid'>"
        + _band_table("core", partition["core_top"])
        + _band_table("accessory", partition["accessory_top"])
        + _band_table("rare", partition["rare_top"])
        + "</div>"
    )


def _render_distribution_figures(
    inputs: Any,
    species_rows: Sequence[Mapping[str, Any]],
) -> str:
    """§3.5 — boxplots of % non-human and ppm by domain, plus cohort mini bar."""
    from csc.report.report import _per_sample_stats  # local import: avoid cycle

    per_sample = _per_sample_stats(inputs)
    abs_enabled = inputs.matrix_abs is not None

    # Per-sample → grouped by *dominant* domain (largest contributor).
    by_dom_pct: dict[str, list[float]] = {}
    by_dom_ppm: dict[str, list[float]] = {}
    for s in per_sample:
        # Largest-contributor domain at the per-sample level
        # (computed from the CPM matrix; falls back to "Unannotated").
        dom = _largest_domain(inputs.matrix_cpm, s["sample_id"])
        if s["nonhuman_pct_of_classified"] is not None:
            by_dom_pct.setdefault(dom, []).append(
                float(s["nonhuman_pct_of_classified"])
            )
        if abs_enabled and s["nonhuman_abs_ppm_total"] is not None:
            by_dom_ppm.setdefault(dom, []).append(
                float(s["nonhuman_abs_ppm_total"])
            )

    pct_groups = sorted(by_dom_pct.items(), key=lambda kv: -len(kv[1]))
    fig_pct = _svg.boxplot_svg(
        pct_groups,
        title=(
            "Figure 3.5a — % non-human reads (compositional, denominator: "
            "per million classified reads), grouped by sample's dominant "
            "non-human domain"
        ),
        y_label="% non-human of classified reads",
    )
    fig_ppm = ""
    if abs_enabled and by_dom_ppm:
        ppm_groups = sorted(by_dom_ppm.items(), key=lambda kv: -len(kv[1]))
        fig_ppm = _svg.boxplot_svg(
            ppm_groups,
            title=(
                "Figure 3.5b — Non-human absolute burden (ppm of total "
                "sequenced reads) by sample's dominant non-human domain"
            ),
            y_label="ppm of total sequenced reads (log10)",
            log_y=True,
        )

    # Cohort-aggregate mini stacked bar
    domain_mix = _cohort.domain_burden_mix(
        species_rows,
        metric="cohort_burden_ppm" if abs_enabled else "cohort_raw_total",
    )
    mini_bar = ""
    if domain_mix:
        mini_bar = (
            "<p class='fig-caption'><strong>Cohort-aggregate composition</strong> "
            "(single bar = entire cohort summed; metric: " +
            ("ppm of total sequenced reads" if abs_enabled else "raw cohort reads") +
            "):</p>"
            + _svg.mini_stacked_bar_svg(domain_mix, width=720)
            + _svg.domain_legend_html(sorted(domain_mix.keys()))
        )

    n_samples = len(inputs.matrix_raw.sample_ids)
    per_sample_warning = (
        f"<p class='fig-caption'>The legacy per-sample stacked-bar figures "
        f"(one row per sample) are <strong>not</strong> rendered by default "
        f"because the cohort has {n_samples} sample(s).  Use "
        f"<code>--layout legacy</code> to obtain them when n &lt; 200.</p>"
    )

    return (
        "<h3>3.5 Cohort-wide distribution figures</h3>"
        + fig_pct
        + fig_ppm
        + mini_bar
        + per_sample_warning
    )


def _largest_domain(matrix: Any, sample_id: str) -> str:
    """Return the taxonomic domain that contributes the most to a sample.

    Excludes ``Human`` so the output is informative for contamination
    grouping.
    """
    totals: dict[str, float] = {}
    for tid in matrix.tax_ids:
        if tid == 9606:
            continue
        v = matrix.values.get(tid, {}).get(sample_id)
        if v is None or v <= 0:
            continue
        d = matrix.tax_domains.get(tid, "Unannotated") or "Unannotated"
        totals[d] = totals.get(d, 0.0) + float(v)
    if not totals:
        return "no non-human"
    return max(totals.items(), key=lambda kv: kv[1])[0]


def _render_heatmap(
    inputs: Any,
    species_rows: Sequence[Mapping[str, Any]],
    *,
    cluster_method: str,
    cluster_distance: str,
    max_samples_cluster: int,
    top_species: int,
) -> str:
    matrix = inputs.matrix_cpm
    species_rows = list(species_rows)[:top_species]
    tax_ids = [r["tax_id"] for r in species_rows]
    sample_ids = list(matrix.sample_ids)
    sub_note = ""
    if len(sample_ids) > max_samples_cluster:
        # Stable deterministic sub-sample: every k-th sample.
        step = math.ceil(len(sample_ids) / max_samples_cluster)
        sample_ids = sample_ids[::step][:max_samples_cluster]
        sub_note = (
            f"<p class='fig-caption'>The cohort has &gt; {max_samples_cluster} "
            f"samples; clustering was performed on a deterministic "
            f"every-{step}-th sub-sample (n = {len(sample_ids)}).  Adjust with "
            f"<code>--max-samples-cluster</code>.</p>"
        )
        logger.info(
            "_render_heatmap: sub-sampled to %d samples (step=%d) "
            "for clustering", len(sample_ids), step,
        )

    if not tax_ids or not sample_ids:
        return (
            "<h3>3.6 Sample × species heatmap</h3>"
            "<p><em>Insufficient data to render heatmap.</em></p>"
        )

    logger.info(
        "_render_heatmap: clustering %d samples × %d species "
        "(method=%s, dist=%s)",
        len(sample_ids), len(tax_ids), cluster_method, cluster_distance,
    )

    # Build value matrix [taxa][samples]
    values = [
        [
            (matrix.values.get(tid, {}).get(sid) or 0.0)
            for sid in sample_ids
        ]
        for tid in tax_ids
    ]

    # Cluster samples via Bray–Curtis on the same top-K subset.
    if cluster_distance == "jaccard":
        D_samp = _jaccard_matrix(values)
    else:  # bray (default)
        D_samp = _bray_from_value_matrix(values)
    sample_clust = _cohort.hclust(D_samp, method=cluster_method)

    # Cluster species via Bray–Curtis on rows (species axis).
    D_taxa = _bray_from_rows(values)
    taxa_clust = _cohort.hclust(D_taxa, method=cluster_method)

    row_labels = [str(r["name"]) for r in species_rows]
    col_labels = list(sample_ids)

    fig = _svg.heatmap_with_dendrogram_svg(
        values,
        row_labels=row_labels,
        col_labels=col_labels,
        row_order=taxa_clust["order"] or list(range(len(row_labels))),
        col_order=sample_clust["order"] or list(range(len(col_labels))),
        title=(
            f"Figure 3.6 — Top {len(row_labels)} species × {len(col_labels)} "
            f"samples, log1p(CPM); rows and columns reordered by "
            f"{cluster_method}-linkage hierarchical clustering "
            f"({cluster_distance} distance)"
        ),
        show_col_labels=len(col_labels) <= 50,
    )
    return (
        "<h3>3.6 Sample × species heatmap (hierarchical clustering)</h3>"
        + fig
        + sub_note
        + "<p class='fig-caption'>Distinct sample clusters in this heatmap "
        "frequently correspond to sequencing batch, kit, or lab "
        "differences rather than biological signal.  Cluster method "
        f"= <code>{_esc(cluster_method)}</code>, distance "
        f"= <code>{_esc(cluster_distance)}</code> on log1p-CPM over the "
        f"top-{len(row_labels)} species by cohort burden.</p>"
    )


def _bray_from_value_matrix(values: Sequence[Sequence[float]]) -> list[list[float]]:
    """Bray–Curtis between samples (columns of `values`)."""
    if not values or not values[0]:
        return []
    n_samp = len(values[0])
    sums = [0.0] * n_samp
    for row in values:
        for j, v in enumerate(row):
            if v is None:
                continue
            try:
                fv = float(v)
            except (TypeError, ValueError):
                continue
            if math.isnan(fv) or fv <= 0:
                continue
            sums[j] += fv
    D = [[0.0] * n_samp for _ in range(n_samp)]
    for i in range(n_samp):
        si = sums[i]
        for j in range(i + 1, n_samp):
            sj = sums[j]
            denom = si + sj
            if denom <= 0:
                D[i][j] = D[j][i] = 0.0
                continue
            num = 0.0
            for row in values:
                a = row[i] if row[i] is not None else 0.0
                b = row[j] if row[j] is not None else 0.0
                a = max(0.0, float(a))
                b = max(0.0, float(b))
                num += abs(a - b)
            D[i][j] = D[j][i] = num / denom
    return D


def _bray_from_rows(values: Sequence[Sequence[float]]) -> list[list[float]]:
    """Bray–Curtis between rows of `values` (i.e. species)."""
    n = len(values)
    if n == 0:
        return []
    sums = []
    for row in values:
        s = 0.0
        for v in row:
            if v is None:
                continue
            try:
                fv = float(v)
            except (TypeError, ValueError):
                continue
            if math.isnan(fv) or fv <= 0:
                continue
            s += fv
        sums.append(s)
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            denom = sums[i] + sums[j]
            if denom <= 0:
                D[i][j] = D[j][i] = 0.0
                continue
            num = 0.0
            for k in range(len(values[0])):
                a = values[i][k] if values[i][k] is not None else 0.0
                b = values[j][k] if values[j][k] is not None else 0.0
                a = max(0.0, float(a))
                b = max(0.0, float(b))
                num += abs(a - b)
            D[i][j] = D[j][i] = num / denom
    return D


def _jaccard_matrix(values: Sequence[Sequence[float]]) -> list[list[float]]:
    """Binary Jaccard distance between sample columns of `values`."""
    if not values or not values[0]:
        return []
    n = len(values[0])
    presence = [[bool(v and v > 0) for v in row] for row in values]
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            inter = union = 0
            for row in presence:
                a, b = row[i], row[j]
                if a or b:
                    union += 1
                    if a and b:
                        inter += 1
            d = 1.0 - (inter / union) if union > 0 else 0.0
            D[i][j] = D[j][i] = d
    return D


def _render_pcoa(
    inputs: Any,
    species_rows: Sequence[Mapping[str, Any]],
    *,
    cluster_distance: str,
    max_samples: int,
    top_species: int,
) -> str:
    matrix = inputs.matrix_cpm
    sample_ids = list(matrix.sample_ids)
    if not sample_ids:
        return ""

    sub_note = ""
    if len(sample_ids) > max_samples:
        step = math.ceil(len(sample_ids) / max_samples)
        sample_ids = sample_ids[::step][:max_samples]
        sub_note = (
            f" (deterministic every-{step}-th sub-sample; n = {len(sample_ids)})"
        )
        logger.info(
            "_render_pcoa: sub-sampled to %d samples (step=%d) for PCoA",
            len(sample_ids), step,
        )

    species_rows = list(species_rows)[:top_species]
    tax_ids = [r["tax_id"] for r in species_rows]
    if not tax_ids or len(sample_ids) < 3:
        return ""
    logger.info(
        "_render_pcoa: computing distance matrix (%d samples, %d species, dist=%s)",
        len(sample_ids), len(tax_ids), cluster_distance,
    )
    values = [
        [(matrix.values.get(tid, {}).get(sid) or 0.0) for sid in sample_ids]
        for tid in tax_ids
    ]
    if cluster_distance == "jaccard":
        D = _jaccard_matrix(values)
    else:
        D = _bray_from_value_matrix(values)
    pcoa = _cohort.pcoa_2d(D)
    coords = pcoa["coords"]
    if not coords:
        return ""

    flagged_set = set()
    if inputs.detect_summary is not None:
        flagged_set = set(inputs.detect_summary.get("flagged_samples") or [])
    if inputs.abs_detect_summary is not None:
        flagged_set |= set(inputs.abs_detect_summary.get("flagged_samples") or [])

    domains = sorted({r["domain"] for r in species_rows}) + ["flagged sample"]
    cmap = _svg.domain_colour_map(domains)
    pts = []
    for i, sid in enumerate(sample_ids):
        dom = _largest_domain(matrix, sid)
        flagged = sid in flagged_set
        colour = cmap.get(
            "flagged sample" if flagged else dom, _svg._colour_for(0)
        )
        pts.append({
            "x": coords[i][0],
            "y": coords[i][1],
            "size": 5 if flagged else 3,
            "colour": colour,
            "label": (
                f"{sid} — dominant {dom}" + (" — flagged" if flagged else "")
            ),
        })
    ev = pcoa["explained_var"]
    ev_str = (
        f" (axes explain {ev[0] * 100:.1f}% / {ev[1] * 100:.1f}% of B-trace)"
        if len(ev) >= 2 else ""
    )
    fig = _svg.scatter_svg(
        pts,
        title=(
            f"Figure 3.7 — PCoA (classical MDS) on {cluster_distance} "
            f"distance over top-{len(tax_ids)} species{sub_note}"
        ),
        x_label=f"PCo1{ev_str}",
        y_label="PCo2",
        legend_domains=domains,
    )
    return (
        "<h3>3.7 Sample β-diversity ordination (PCoA)</h3>"
        + fig
        + "<p class='fig-caption'>Each dot is a sample, projected by "
        "classical multidimensional scaling on a "
        f"<code>{_esc(cluster_distance)}</code>-distance matrix over the "
        "top-K species.  Tight clusters / streaks here often co-localise "
        "with batch / kit / lab effects when cross-referenced with the "
        "heatmap (§3.6); samples flagged by <code>csc-detect</code> are "
        "highlighted.</p>"
    )


def _render_species_drilldown(
    inputs: Any,
    species_rows: Sequence[Mapping[str, Any]],
    *,
    top_n: int,
    abs_enabled: bool,
) -> str:
    parts = ["<h3>3.8 Per-species drill-down</h3>"]
    parts.append(
        "<p>Collapsible details for the top "
        f"{min(top_n, len(species_rows))} species by cohort burden.  Each "
        "section lists the most-affected samples for that species (top 10 "
        "by absolute ppm when available, else by raw reads).</p>"
    )
    matrix = inputs.matrix_abs or inputs.matrix_raw
    metric_label = "ppm of total sequenced reads" if inputs.matrix_abs is not None else "raw reads"
    for r in list(species_rows)[:top_n]:
        tid = r["tax_id"]
        per_samp = []
        for sid in matrix.sample_ids:
            v = matrix.values.get(tid, {}).get(sid)
            if v is None or v <= 0:
                continue
            per_samp.append((sid, float(v)))
        per_samp.sort(key=lambda kv: -kv[1])
        body_rows = "".join(
            f"<tr><td>{_esc(s)}</td><td>{_fmt_num(v)}</td></tr>"
            for s, v in per_samp[:10]
        ) or "<tr><td colspan='2'><em>No positive samples.</em></td></tr>"
        parts.append(
            f"<details class='species-details'>"
            f"<summary>{_esc(r['name'])} (tax_id {r['tax_id']}, "
            f"{_esc(r['domain'])}) — prevalence {_fmt_pct(r['prevalence'])}, "
            f"cohort burden { _fmt_num(r['cohort_burden_ppm']) if abs_enabled else _fmt_int(int(r['cohort_raw_total'])) }"
            f"{ ' ppm' if abs_enabled else ' reads' }</summary>"
            f"<table><thead><tr><th>Sample</th>"
            f"<th>{_esc(metric_label)}</th></tr></thead>"
            f"<tbody>{body_rows}</tbody></table>"
            f"</details>"
        )
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# §4 Variant-calling impact (cohort-aware)
# ---------------------------------------------------------------------------


def render_variant_impact_v2(
    inputs: Any,
    per_sample_stats: Sequence[Mapping[str, Any]],
    species_rows: Sequence[Mapping[str, Any]],
    threshold_ppm: float,
    page_size: int = 25,
    *,
    flagged_table_top: int = 100,
    flagged_tsv_path: Path | None = None,
) -> tuple[str, list[Mapping[str, Any]]]:
    abs_enabled = inputs.matrix_abs is not None
    if not abs_enabled:
        body = (
            "<h2>4. Variant-Calling Impact</h2>"
            "<div class='callout'>This section requires the "
            "absolute-burden matrix (<code>taxa_matrix_abs.tsv</code>).  "
            "Re-run <code>csc-aggregate --idxstats …</code> with the "
            "per-sample <code>reads_summary.json</code> sidecars to "
            "enable.</div>"
        )
        return body, []

    flagged = sorted(
        [s for s in per_sample_stats
         if s["nonhuman_abs_ppm_total"] is not None
         and s["nonhuman_abs_ppm_total"] > threshold_ppm],
        key=lambda s: -(s["nonhuman_abs_ppm_total"] or 0),
    )

    burdens = [
        s["nonhuman_abs_ppm_total"]
        for s in per_sample_stats
        if s["nonhuman_abs_ppm_total"] is not None
    ]
    edges, counts = _cohort.histogram(burdens, n_bins=30, log=True) if burdens else ([0, 1], [])
    hist = _svg.histogram_svg(
        edges,
        counts,
        title=("Figure 4 — Distribution of non-human burden across the cohort "
               "(log10 ppm of total sequenced reads); orange dashed line = "
               "configured threshold"),
        x_label="ppm of total sequenced reads (log10)",
        y_label="# samples",
        vline=threshold_ppm if burdens else None,
        log_x=True,
    )

    # Species attribution of cohort burden ABOVE threshold.
    flagged_ids = {s["sample_id"] for s in flagged}
    matrix = inputs.matrix_abs
    attribution: dict[str, float] = {}
    if matrix is not None and flagged_ids:
        for r in species_rows:
            tid = r["tax_id"]
            tot = 0.0
            for sid in flagged_ids:
                v = matrix.values.get(tid, {}).get(sid)
                if v is None or v <= 0:
                    continue
                tot += float(v)
            if tot > 0:
                attribution[r["name"]] = tot
    attribution_bar = ""
    if attribution:
        # Keep top 10 species, group the rest as "other".
        items = sorted(attribution.items(), key=lambda kv: -kv[1])
        if len(items) > 10:
            other = sum(v for _, v in items[10:])
            items = items[:10] + [("other species", other)]
        comp = dict(items)
        attribution_bar = (
            "<p class='fig-caption'>Species attribution of cohort burden "
            "above the threshold (single horizontal bar; ppm of total "
            "sequenced reads, summed across the flagged samples):</p>"
            + _svg.mini_stacked_bar_svg(comp, width=720)
            + _svg.domain_legend_html(list(comp.keys()))
        )

    # Always write the full flagged-samples list to a TSV sidecar so
    # the HTML can stay small while reviewers retain full auditability.
    sidecar_link = ""
    if flagged_tsv_path is not None:
        try:
            _write_variant_impact_tsv(flagged, flagged_tsv_path)
            sidecar_link = (
                f"  Full list ({_fmt_int(len(flagged))} sample(s)) is "
                f"written to <a class='dl-link' href='{_esc(flagged_tsv_path.name)}'>"
                f"{_esc(flagged_tsv_path.name)}</a>."
            )
        except OSError as exc:  # pragma: no cover - defensive
            logger.warning("Could not write variant-impact TSV: %s", exc)

    n_flagged = len(flagged)
    truncated = n_flagged > flagged_table_top
    display_flagged = flagged[:flagged_table_top] if truncated else flagged

    rows_html = []
    for s in display_flagged:
        rows_html.append(
            "<tr>"
            f"<td>{_esc(s['sample_id'])}</td>"
            f"<td>{_fmt_int(int(s['total_reads']) if s['total_reads'] else None)}</td>"
            f"<td data-sort='{s['nonhuman_abs_ppm_total']}'>"
            f"{s['nonhuman_abs_ppm_total']:,.2f}</td>"
            f"<td>{(s['nonhuman_abs_ppm_total'] or 0) / 1e4:.4f}%</td>"
            "</tr>"
        )
    if not rows_html:
        rows_html.append(
            "<tr><td colspan='4'><em>No samples exceed the threshold.</em></td></tr>"
        )
    truncation_note = ""
    if truncated:
        truncation_note = (
            f"<div class='callout'>Showing the top "
            f"<strong>{flagged_table_top}</strong> flagged samples by "
            f"absolute burden; <strong>{_fmt_int(n_flagged - flagged_table_top)}</strong> "
            f"additional flagged sample(s) are present in the sidecar TSV "
            f"only.</div>"
        )
    table = (
        "<h4>Flagged samples (top by burden)</h4>"
        "<p>Sortable and filterable.  Default sort: highest absolute "
        f"burden first.{sidecar_link}</p>"
        + truncation_note
        + f'<table class="paginated" data-page-size="{page_size}" '
        'data-default-sort="2|desc">'
        "<thead><tr><th>Sample</th>"
        "<th>Total sequenced reads</th>"
        "<th>Non-human burden (ppm)<span class='denom'>ppm of total sequenced reads</span></th>"
        "<th>Non-human burden (% of sequencing)</th>"
        "</tr></thead><tbody>" + "\n".join(rows_html) + "</tbody></table>"
    )

    msg = (
        f"{n_flagged} sample(s) exceed the configured threshold of "
        f"{threshold_ppm:,.0f} ppm of total sequenced reads "
        f"({threshold_ppm / 1e4:.3f}% of sequencing)."
    )
    callout_cls = "danger" if flagged else "ok"
    body = (
        "<h2>4. Variant-Calling Impact</h2>"
        "<p>The <strong>absolute</strong> fraction of sequencing assigned "
        "to non-human taxa is the metric reviewers typically scrutinise "
        "when judging whether contamination can affect variant calls, "
        "assemblies, or expression estimates.  The threshold below "
        "(configurable via <code>--variant-impact-threshold-ppm</code>) "
        "should match the sensitivity of your downstream application.</p>"
        f"<div class='callout {callout_cls}'>{msg}</div>"
        + hist
        + attribution_bar
        + table
    )
    return body, list(flagged)


def _write_variant_impact_tsv(
    flagged: Sequence[Mapping[str, Any]],
    tsv_path: Path,
) -> None:
    """Write the full §4 flagged-samples list as a TSV sidecar."""
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "sample_id",
            "total_reads",
            "nonhuman_abs_ppm_total",
            "nonhuman_pct_of_sequencing",
        ])
        for s in flagged:
            ppm = s["nonhuman_abs_ppm_total"]
            w.writerow([
                str(s["sample_id"]),
                str(int(s["total_reads"])) if s["total_reads"] is not None else "NA",
                f"{ppm:.6f}" if ppm is not None else "NA",
                f"{(ppm or 0) / 1e4:.6f}" if ppm is not None else "NA",
            ])


# ---------------------------------------------------------------------------
# §5 Detection summary
# ---------------------------------------------------------------------------


def render_detection_section_v2(inputs: Any) -> str:
    """Promote the legacy 3.3 block to its own §5 with an UpSet-style banner."""
    primary = inputs.detect_summary
    abs_d = inputs.abs_detect_summary
    if primary is None and abs_d is None:
        return (
            "<h2>5. Detection summary</h2>"
            "<p><em>No <code>csc-detect</code> outputs were supplied; "
            "this section is empty.</em></p>"
        )

    # Reuse the existing implementation for the body; it already covers
    # the side-by-side flag table, top-flagged-taxa table, and set
    # comparison.  We just relocate the heading.
    from csc.report.report import _render_detection_comparison

    body = _render_detection_comparison(inputs)
    # Remove the legacy "3.3 Outlier-detection results" subheading and
    # promote to §5.
    body = body.replace("<h3>3.3 Outlier-detection results</h3>", "")

    # Add an UpSet-style 2-set tile diagram when both passes ran.
    upset = ""
    if primary is not None and abs_d is not None:
        primary_set = set(primary.get("flagged_samples") or [])
        abs_set = set(abs_d.get("flagged_samples") or [])
        both = primary_set & abs_set
        only_primary = primary_set - abs_set
        only_abs = abs_set - primary_set
        upset = (
            "<div class='quadrant-grid' style='grid-template-columns: repeat(3, 1fr);'>"
            f"<div class='quadrant'><h4>Both passes</h4>"
            f"<p style='font-size:22px;margin:0.2em 0'>{len(both)}</p>"
            "<p>High-confidence contamination calls (independent of denominator).</p></div>"
            f"<div class='quadrant'><h4>Primary only</h4>"
            f"<p style='font-size:22px;margin:0.2em 0'>{len(only_primary)}</p>"
            "<p>Unusual <em>composition</em> within non-human reads.</p></div>"
            f"<div class='quadrant'><h4>Abs only</h4>"
            f"<p style='font-size:22px;margin:0.2em 0'>{len(only_abs)}</p>"
            "<p>Unusual <em>total</em> non-human read count — typically "
            "host-depletion variability.</p></div>"
            "</div>"
        )

    return "<h2>5. Detection summary</h2>" + upset + body


# ---------------------------------------------------------------------------
# §6 Per-sample appendix (paginated + sidecar TSV)
# ---------------------------------------------------------------------------


def render_per_sample_appendix(
    inputs: Any,
    per_sample_stats: Sequence[Mapping[str, Any]],
    *,
    page_size: int,
    tsv_path: Path,
    notable_top: int = 15,
) -> str:
    """Render §6 — *Notable samples & per-sample sidecar*.

    For large cohorts (thousands of samples) emitting one HTML row per
    sample produces a multi-megabyte, unresponsive appendix.  Instead,
    this section:

    * Always writes the **complete** per-sample table to ``tsv_path``
      (same columns and units as before) so reviewers can audit every
      sample offline.
    * Renders only **curated highlight tables** in the HTML — top
      samples by absolute burden, by composition, by richness/diversity,
      and any sample(s) flagged by ``csc-detect`` — so the static
      document stays focused on actionable cohort summaries and on the
      worst / most informative individual samples.

    The ``notable_top`` argument bounds each highlight table.  The
    legacy ``page_size`` argument is retained for forward-compatibility
    (it no longer governs an HTML table here, but is documented in the
    sidecar caption alongside the link to the TSV).
    """
    abs_enabled = inputs.matrix_abs is not None

    # --- Always-on TSV sidecar (unchanged columns) -------------------
    tsv_header = [
        "sample_id", "total_reads", "classified_reads", "nonhuman_reads",
        "nonhuman_pct_of_classified",
    ]
    if abs_enabled:
        tsv_header.append("nonhuman_abs_ppm_total")
    tsv_header += ["richness", "shannon", "simpson"]

    tsv_rows: list[list[str]] = []
    for s in per_sample_stats:
        row = [
            str(s["sample_id"]),
            str(s["total_reads"]) if s["total_reads"] is not None else "NA",
            str(int(s["classified_reads"])),
            str(int(s["nonhuman_reads"])),
            "NA" if s["nonhuman_pct_of_classified"] is None
            else f"{s['nonhuman_pct_of_classified']:.6f}",
        ]
        if abs_enabled:
            v = s["nonhuman_abs_ppm_total"]
            row.append("NA" if v is None else f"{v:.6f}")
        row += [
            str(s["richness"]),
            f"{s['shannon']:.6f}",
            f"{s['simpson']:.6f}",
        ]
        tsv_rows.append(row)

    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(tsv_header)
        writer.writerows(tsv_rows)

    # --- Curated highlight tables ------------------------------------
    n = len(per_sample_stats)
    # Detect-flagged sample IDs (union over both passes when available)
    flagged_ids: set[str] = set()
    if inputs.detect_summary is not None:
        flagged_ids |= set(inputs.detect_summary.get("flagged_samples") or [])
    if inputs.abs_detect_summary is not None:
        flagged_ids |= set(inputs.abs_detect_summary.get("flagged_samples") or [])

    parts: list[str] = ["<h2>6. Notable samples &amp; per-sample sidecar</h2>"]
    parts.append(
        "<p>To keep the report easily renderable for cohorts of "
        f"thousands of samples, the per-sample HTML table is "
        "intentionally <strong>not</strong> emitted here.  The "
        f"complete per-sample summary ({_fmt_int(n)} sample(s) × "
        f"{len(tsv_header)} columns) is written alongside this report "
        f"as <a class='dl-link' href='{_esc(tsv_path.name)}'>"
        f"{_esc(tsv_path.name)}</a> for offline analysis.  The "
        "highlight tables below surface a small number of samples that "
        "are most useful for interpreting the cohort: outliers by "
        "burden / composition, the most diverse non-human metagenomes, "
        "and any samples that <code>csc-detect</code> flagged as "
        "contamination calls.</p>"
    )

    # Highlight A — top by absolute burden (variant-impact relevant)
    if abs_enabled:
        top_burden = sorted(
            (s for s in per_sample_stats
             if s["nonhuman_abs_ppm_total"] is not None),
            key=lambda s: -float(s["nonhuman_abs_ppm_total"] or 0.0),
        )[:notable_top]
        parts.append(_render_notable_table(
            top_burden,
            heading=(
                f"6.1 Highest absolute non-human burden "
                f"(top {len(top_burden)} of {n})"
            ),
            description=(
                "Samples whose non-human reads make up the largest "
                "fraction of <em>total sequencing</em> — these are the "
                "most likely to perturb variant calls, assemblies and "
                "expression estimates."
            ),
            abs_enabled=abs_enabled,
            flagged_ids=flagged_ids,
            sort_metric="abs",
        ))

    # Highlight B — top by composition (% non-human of classified)
    top_comp = sorted(
        (s for s in per_sample_stats
         if s["nonhuman_pct_of_classified"] is not None),
        key=lambda s: -float(s["nonhuman_pct_of_classified"] or 0.0),
    )[:notable_top]
    parts.append(_render_notable_table(
        top_comp,
        heading=(
            f"6.{2 if abs_enabled else 1} Highest non-human composition "
            f"(top {len(top_comp)} of {n})"
        ),
        description=(
            "Samples whose <em>classified-read</em> fraction is dominated "
            "by non-human taxa.  When the absolute burden is small but "
            "this percentage is high, the sample is typically a "
            "successfully host-depleted or low-yield WGS run rather than "
            "a contamination event — cross-reference with the absolute "
            "burden column."
        ),
        abs_enabled=abs_enabled,
        flagged_ids=flagged_ids,
        sort_metric="pct",
    ))

    # Highlight C — most diverse non-human content (richness × Shannon)
    # Use Shannon (more robust to small richness inflation from kitome).
    top_div = sorted(
        per_sample_stats,
        key=lambda s: -float(s.get("shannon") or 0.0),
    )[:notable_top]
    sec_idx = 3 if abs_enabled else 2
    parts.append(_render_notable_table(
        top_div,
        heading=(
            f"6.{sec_idx} Most diverse non-human metagenomes "
            f"(top {len(top_div)} by Shannon H)"
        ),
        description=(
            "Samples with the most varied non-human taxa.  High diversity "
            "concurrent with high burden often reflects environmental / "
            "skin commensal carry-over; high diversity at low burden is "
            "characteristic of trace kitome backgrounds."
        ),
        abs_enabled=abs_enabled,
        flagged_ids=flagged_ids,
        sort_metric="shannon",
    ))

    # Highlight D — detect-flagged exemplars (when present)
    if flagged_ids:
        flagged_samples = [
            s for s in per_sample_stats if s["sample_id"] in flagged_ids
        ]
        # Sort flagged exemplars by burden (or composition fallback) so
        # the most actionable ones appear first.
        if abs_enabled:
            flagged_samples.sort(
                key=lambda s: -float(s.get("nonhuman_abs_ppm_total") or 0.0)
            )
        else:
            flagged_samples.sort(
                key=lambda s: -float(s.get("nonhuman_pct_of_classified") or 0.0)
            )
        n_flagged_total = len(flagged_samples)
        flagged_samples = flagged_samples[:notable_top]
        sec_idx = 4 if abs_enabled else 3
        parts.append(_render_notable_table(
            flagged_samples,
            heading=(
                f"6.{sec_idx} Samples flagged by csc-detect "
                f"(showing {len(flagged_samples)} of {n_flagged_total})"
            ),
            description=(
                "The union of samples flagged by the primary and "
                "absolute-burden detect passes.  These are the "
                "statistically-defined contamination calls; review the "
                "per-taxon flag table in §5 for the species responsible."
            ),
            abs_enabled=abs_enabled,
            flagged_ids=flagged_ids,
            sort_metric="abs" if abs_enabled else "pct",
        ))

    return "\n".join(parts)


def _render_notable_table(
    samples: Sequence[Mapping[str, Any]],
    *,
    heading: str,
    description: str,
    abs_enabled: bool,
    flagged_ids: set[str],
    sort_metric: str,
) -> str:
    """Render a small, static highlight table of notable samples."""
    headers = [
        "Sample",
        "Total sequenced reads<span class='denom'>idxstats</span>",
        "Classified reads<span class='denom'>Kraken2 direct</span>",
        "Non-human reads<span class='denom'>classified − Homo sapiens</span>",
        "% Non-human<span class='denom'>per million classified reads denominator</span>",
    ]
    if abs_enabled:
        headers.append(
            "Non-human burden <span class='denom'>ppm of total sequenced reads</span>"
        )
    headers.extend([
        "Richness<span class='denom'>taxa with ≥1 read</span>",
        "Shannon H<span class='denom'>natural log</span>",
        "Simpson (1−D)<span class='denom'>Gini–Simpson</span>",
        "csc-detect flag",
    ])
    head_html = "".join(f"<th>{h}</th>" for h in headers)

    rows_html: list[str] = []
    for s in samples:
        pct = s["nonhuman_pct_of_classified"]
        pct_cell = "NA" if pct is None else f"{pct:.3f}"
        cells = [
            f"<td>{_esc(s['sample_id'])}</td>",
            f"<td>{_fmt_int(int(s['total_reads']) if s['total_reads'] is not None else None)}</td>",
            f"<td>{_fmt_int(int(s['classified_reads']))}</td>",
            f"<td>{_fmt_int(int(s['nonhuman_reads']))}</td>",
            f"<td>{pct_cell}</td>",
        ]
        if abs_enabled:
            v = s["nonhuman_abs_ppm_total"]
            v_cell = "NA" if v is None else f"{v:,.2f}"
            cells.append(f"<td>{v_cell}</td>")
        cells.extend([
            f"<td>{_fmt_int(s['richness'])}</td>",
            f"<td>{s['shannon']:.3f}</td>",
            f"<td>{s['simpson']:.3f}</td>",
            "<td>" + ("⚑ flagged" if s["sample_id"] in flagged_ids else "—") + "</td>",
        ])
        rows_html.append("<tr>" + "".join(cells) + "</tr>")
    if not rows_html:
        rows_html.append(
            f"<tr><td colspan='{len(headers)}'><em>No samples to show.</em>"
            "</td></tr>"
        )

    return (
        f"<h3>{_esc(heading)}</h3>"
        f"<p>{description}</p>"
        f"<table><thead><tr>{head_html}</tr></thead>"
        f"<tbody>{chr(10).join(rows_html)}</tbody></table>"
    )


# ---------------------------------------------------------------------------
# §7 Discussion (auto-keyed off cohort partition)
# ---------------------------------------------------------------------------


def render_discussion_v2(
    inputs: Any,
    species_rows: Sequence[Mapping[str, Any]],
    partition: Mapping[str, Any],
    variant_flagged: Sequence[Mapping[str, Any]],
) -> str:
    abs_enabled = inputs.matrix_abs is not None
    parts = [
        "<h2>7. Discussion &amp; Caveats</h2>",
        "<ul>",
        "<li>Kraken2 assignments are only as reliable as the reference "
        "database; use a non-redundant database (e.g. PrackenDB) where "
        "possible to avoid LCA inflation.</li>",
        "<li>Compositional metrics (CPM, stacked bars, diversity) describe "
        "the <em>shape</em> of the non-human fraction but are invariant "
        "to its size.  Always pair CPM with absolute burden (§4) when "
        "available.</li>",
        "<li>Laboratory and reagent contaminants ('kitome') frequently "
        "appear as <strong>core, low-burden</strong> species in §3.2; "
        "consider adding them to <code>csc-detect --kitome-taxa</code>.</li>",
        "<li>Per-sample handling errors typically appear as "
        "<strong>outburst</strong> species (low prevalence, very high "
        "absolute burden in a handful of samples) — cross-reference §3.2 "
        "with §4.</li>",
    ]

    # Auto-generate species-keyed advice from the partition.
    if partition.get("core_top"):
        names = ", ".join(_esc(r["name"]) for r in partition["core_top"][:3])
        parts.append(
            f"<li><strong>Core species in this cohort</strong>: {names}.  "
            "Their high prevalence at low-to-moderate burden is the "
            "fingerprint of reagent / library contamination (kitome).  "
            "Consider whether they should be excluded from any biological "
            "interpretation or added to a kitome ignore list.</li>"
        )
    if partition.get("rare_top"):
        rare_high = [
            r for r in partition["rare_top"]
            if (r.get("cohort_burden_ppm") or 0) > 0
        ][:3]
        if rare_high:
            names = ", ".join(_esc(r["name"]) for r in rare_high)
            parts.append(
                f"<li><strong>Outburst species</strong>: {names}.  Rare "
                "across the cohort but with substantial burden where they "
                "appear — investigate the affected samples (§3.8 / §6) "
                "for handling errors or sample-specific contamination.</li>"
            )

    if not abs_enabled:
        parts.append(
            "<li><strong>No absolute-burden denominator was available.</strong>  "
            "Any claim about the magnitude of contamination should be "
            "deferred until <code>reads_summary.json</code> sidecars are "
            "supplied to <code>csc-aggregate</code>.</li>"
        )
    if variant_flagged:
        parts.append(
            f"<li>{len(variant_flagged)} sample(s) were flagged in §4 as "
            "potentially impacting variant calling.  These should be "
            "inspected manually before downstream use.</li>"
        )

    if inputs.abs_detect_summary is not None and inputs.detect_summary is not None:
        primary_set = set(inputs.detect_summary.get("flagged_samples") or [])
        abs_set = set(inputs.abs_detect_summary.get("flagged_samples") or [])
        only_abs = abs_set - primary_set
        if only_abs:
            parts.append(
                f"<li>{len(only_abs)} sample(s) were flagged "
                "<strong>only</strong> by the absolute-burden side pass — "
                "characteristic of variable host-depletion efficiency rather "
                "than a single-organism contaminant.</li>"
            )

    parts.append("</ul>")

    # Confidence-tier auto-paragraph: when a high-confidence sibling
    # bundle is available, append a comparison of the two flag sets
    # (this paragraph is emitted only when both tiers exist; otherwise
    # the section is omitted).
    if getattr(inputs, "confidence_tiers", None):
        parts.append(_render_discussion_confidence_paragraph(inputs))

    return "\n".join(parts)


def _render_discussion_confidence_paragraph(inputs: Any) -> str:
    """Auto-generated paragraph comparing sensitive vs high-confidence tiers."""
    primary_set = set((inputs.detect_summary or {}).get("flagged_samples") or [])
    abs_set = set((inputs.abs_detect_summary or {}).get("flagged_samples") or [])

    bullets: list[str] = []
    for tier_suffix, tier_inputs in inputs.confidence_tiers.items():
        t_primary = set(
            (tier_inputs.detect_summary or {}).get("flagged_samples") or []
        )
        t_abs = set(
            (tier_inputs.abs_detect_summary or {}).get("flagged_samples") or []
        )
        only_sensitive = primary_set - t_primary
        only_high = t_primary - primary_set
        shared = primary_set & t_primary
        thresh = tier_inputs.tier_threshold
        bullets.append(
            f"<li>At Kraken2 confidence ≥ {thresh:.2f} (tier "
            f"<code>{_esc(tier_suffix)}</code>), the primary detect pass "
            f"shares <strong>{len(shared)}</strong> flagged sample(s) "
            f"with the sensitive tier; <strong>{len(only_sensitive)}</strong> "
            f"sample(s) are flagged only at sensitive confidence (likely "
            f"low-complexity / sparse-kmer false positives) and "
            f"<strong>{len(only_high)}</strong> sample(s) emerge only at "
            f"high confidence (genuine, well-supported hits that the "
            f"unfiltered tier may bury under noise).</li>"
        )
        if abs_set or t_abs:
            shared_abs = abs_set & t_abs
            bullets.append(
                f"<li>The absolute-burden side pass concurs on "
                f"<strong>{len(shared_abs)}</strong> sample(s) across the "
                f"two tiers — a robust indicator of true non-human burden "
                f"that survives both denominator and confidence "
                f"perturbations.</li>"
            )
    if not bullets:
        return ""
    return (
        "<h3>Sensitive vs high-confidence agreement</h3>"
        "<p>Findings that survive a stricter Kraken2 confidence threshold "
        "are the most defensible contamination calls.  Findings unique "
        "to the sensitive tier should be inspected manually before "
        "exclusion — sensitive-only flags caused by low-complexity "
        "repeats are common in shotgun WGS and are exactly what the "
        "high-confidence tier is designed to suppress (Wood <em>et al.</em> "
        "2019; Marcelino <em>et al.</em> 2020).</p>"
        "<ul>" + "\n".join(bullets) + "</ul>"
    )


# ---------------------------------------------------------------------------
# Sensitive vs High-Confidence concordance section
# ---------------------------------------------------------------------------


def render_confidence_concordance(inputs: Any) -> str:
    """Render the §5.x "Sensitive vs High-Confidence concordance" section.

    Always emits both tiers' summary data side-by-side regardless of
    which tier the toggle is currently showing, so the user can compare
    counts at a glance.  Returns an empty string when no high-confidence
    tier sibling is present.
    """
    tiers = getattr(inputs, "confidence_tiers", None) or {}
    if not tiers:
        return ""

    parts: list[str] = [
        "<h2 class='confidence-concordance'>5.1 Sensitive vs "
        "High-Confidence concordance</h2>",
        "<p>Kraken2's sensitive setting (<code>--confidence 0.0</code>) "
        "is required for trace-contamination detection but produces "
        "false positives caused by low-complexity regions or sparse "
        "k-mer matches.  The high-confidence tier(s) below were "
        "produced by recomputing the per-read Kraken2 confidence from "
        "the existing per-read output files and demoting reads scoring "
        "below the threshold to <em>unclassified</em> "
        "(see §2 Methods).  Findings that persist across tiers are the "
        "most defensible.</p>",
    ]

    primary = inputs.detect_summary
    abs_d = inputs.abs_detect_summary
    primary_set = set((primary or {}).get("flagged_samples") or [])
    abs_set = set((abs_d or {}).get("flagged_samples") or [])

    # Per-tier flag-set Venn-style summary.
    rows: list[str] = [
        "<thead><tr>"
        "<th>Tier</th>"
        "<th>Threshold</th>"
        "<th>Taxa retained</th>"
        "<th>Samples flagged (primary)</th>"
        "<th>… ∩ sensitive</th>"
        "<th>… only at this tier</th>"
        "<th>Samples flagged (abs side pass)</th>"
        "</tr></thead><tbody>",
        "<tr>"
        "<td><strong>Sensitive</strong> "
        "(<code>--confidence 0.0</code>)</td>"
        "<td>0.00</td>"
        f"<td>{_fmt_int(len(inputs.matrix_raw.tax_ids))}</td>"
        f"<td>{_fmt_int(len(primary_set))}</td>"
        "<td>—</td><td>—</td>"
        f"<td>{_fmt_int(len(abs_set))}</td>"
        "</tr>",
    ]
    for tier_suffix, tier_inputs in tiers.items():
        t_primary = set(
            (tier_inputs.detect_summary or {}).get("flagged_samples") or []
        )
        t_abs = set(
            (tier_inputs.abs_detect_summary or {}).get("flagged_samples") or []
        )
        rows.append(
            "<tr>"
            f"<td><strong>{_esc(tier_suffix)}</strong></td>"
            f"<td>{tier_inputs.tier_threshold:.2f}</td>"
            f"<td>{_fmt_int(len(tier_inputs.matrix_raw.tax_ids))}</td>"
            f"<td>{_fmt_int(len(t_primary))}</td>"
            f"<td>{_fmt_int(len(t_primary & primary_set))}</td>"
            f"<td>{_fmt_int(len(t_primary - primary_set))}</td>"
            f"<td>{_fmt_int(len(t_abs))}</td>"
            "</tr>"
        )
    rows.append("</tbody>")
    parts.append("<table>" + "\n".join(rows) + "</table>")

    # 2-set quadrant tile per tier (primary detect pass).
    for tier_suffix, tier_inputs in tiers.items():
        t_primary = set(
            (tier_inputs.detect_summary or {}).get("flagged_samples") or []
        )
        both = primary_set & t_primary
        only_sens = primary_set - t_primary
        only_high = t_primary - primary_set
        parts.append(
            "<h3>Primary-pass flag overlap — "
            f"sensitive ↔ <code>{_esc(tier_suffix)}</code> "
            f"(threshold {tier_inputs.tier_threshold:.2f})</h3>"
            "<div class='quadrant-grid' "
            "style='grid-template-columns: repeat(3, 1fr);'>"
            "<div class='quadrant'><h4>Both tiers</h4>"
            f"<p style='font-size:22px;margin:0.2em 0'>{len(both)}</p>"
            "<p>Robust contamination calls — flagged regardless "
            "of confidence stringency.</p></div>"
            "<div class='quadrant'><h4>Sensitive only</h4>"
            f"<p style='font-size:22px;margin:0.2em 0'>{len(only_sens)}</p>"
            "<p>Likely low-complexity / sparse-kmer false positives "
            "the high-confidence tier filters out.</p></div>"
            "<div class='quadrant'><h4>High-confidence only</h4>"
            f"<p style='font-size:22px;margin:0.2em 0'>{len(only_high)}</p>"
            "<p>Genuine signal that the sensitive tier may have "
            "buried under noise.</p></div>"
            "</div>"
        )

    # Per-sample reads_demoted_to_unclassified counts from
    # aggregation_metadata.json (auditability).  Only emit this table
    # when the metadata exposes the counters – older aggregations may
    # lack them.
    meta_tiers = (inputs.aggregate_metadata or {}).get("confidence_tiers") or {}
    if meta_tiers:
        demoted_rows: list[str] = []
        for tier_suffix, tier_meta in meta_tiers.items():
            demoted = tier_meta.get("reads_demoted_to_unclassified") or {}
            if not demoted:
                continue
            sample_ids = sorted(demoted.keys())
            for sid in sample_ids[:_DEMOTED_TABLE_MAX_ROWS]:
                demoted_rows.append(
                    f"<tr><td>{_esc(sid)}</td>"
                    f"<td>{_esc(tier_suffix)}</td>"
                    f"<td>{_fmt_int(int(demoted[sid]))}</td></tr>"
                )
            if len(sample_ids) > _DEMOTED_TABLE_MAX_ROWS:
                demoted_rows.append(
                    "<tr><td colspan='3'><em>… "
                    f"{len(sample_ids) - _DEMOTED_TABLE_MAX_ROWS} additional "
                    "sample(s) omitted; see aggregation_metadata.json for "
                    "the full list.</em></td></tr>"
                )
        if demoted_rows:
            parts.append(
                "<h3>Reads demoted to unclassified per sample</h3>"
                "<p>Counts from <code>aggregation_metadata.json</code> "
                "(<code>confidence_tiers.&lt;suffix&gt;."
                "reads_demoted_to_unclassified</code>).  These reads "
                "passed Kraken2 but failed the high-confidence "
                "threshold — the larger this number, the more "
                "low-complexity content was filtered out.</p>"
                "<table><thead><tr>"
                "<th>Sample</th><th>Tier</th>"
                "<th>Reads demoted</th>"
                "</tr></thead><tbody>"
                + "\n".join(demoted_rows)
                + "</tbody></table>"
            )

    return "\n".join(parts)
