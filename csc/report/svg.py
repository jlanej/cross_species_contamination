"""Inline-SVG figure primitives for the cohort report.

All figures are emitted as standalone ``<svg>`` strings (no external
files, no JS).  Captions and figure containers are added by the
caller – these helpers focus on the geometry only.
"""

from __future__ import annotations

import html
import math
from typing import Iterable, Mapping, Sequence


# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------


# Okabe–Ito + a few extras.  Matched (in order) to the palette used by
# ``csc.report.report._domain_colour`` so colour assignments stay
# consistent across the report.
_DOMAIN_PALETTE = (
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999", "#000000", "#66C2A5",
)


def _colour_for(idx: int) -> str:
    return _DOMAIN_PALETTE[idx % len(_DOMAIN_PALETTE)]


def _esc(s: object) -> str:
    return html.escape(str(s))


def domain_colour_map(domains: Iterable[str]) -> dict[str, str]:
    """Stable colour mapping for the given domain labels."""
    seen: dict[str, str] = {}
    for d in domains:
        if d not in seen:
            seen[d] = _colour_for(len(seen))
    return seen


def _safe_log10(x: float, fallback: float = -6.0) -> float:
    if x is None or x <= 0 or math.isnan(x):
        return fallback
    return math.log10(x)


# ---------------------------------------------------------------------------
# Sparkline
# ---------------------------------------------------------------------------


def sparkline_svg(
    counts: Sequence[int],
    *,
    width: int = 110,
    height: int = 22,
    fill: str = "#56B4E9",
) -> str:
    """Tiny vertical-bar sparkline of histogram counts.

    Used inside the species summary table to convey the shape of the
    per-sample CPM distribution without spending a row of vertical
    real-estate on a full chart.
    """
    n = len(counts)
    if n == 0 or max(counts) <= 0:
        return (
            f'<svg width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" '
            f'xmlns="http://www.w3.org/2000/svg" aria-label="empty">'
            f'<line x1="0" y1="{height - 1}" x2="{width}" y2="{height - 1}" '
            f'stroke="#bbb"/></svg>'
        )
    mx = max(counts)
    bw = width / n
    parts = [
        f'<svg width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" '
        f'aria-label="distribution sparkline">'
    ]
    for i, c in enumerate(counts):
        if c <= 0:
            continue
        h = (c / mx) * (height - 2)
        x = i * bw
        y = height - h - 1
        parts.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" '
            f'width="{max(bw - 0.6, 0.6):.2f}" height="{h:.2f}" fill="{fill}"/>'
        )
    parts.append('</svg>')
    return "".join(parts)


# ---------------------------------------------------------------------------
# Mini horizontal stacked bar (executive summary domain mix)
# ---------------------------------------------------------------------------


def mini_stacked_bar_svg(
    composition: Mapping[str, float],
    *,
    width: int = 360,
    height: int = 18,
) -> str:
    items = [(k, float(v)) for k, v in composition.items() if v and v > 0]
    if not items:
        return ""
    items.sort(key=lambda kv: -kv[1])
    total = sum(v for _, v in items)
    cmap = domain_colour_map(k for k, _ in items)
    parts = [
        f'<svg width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img" '
        f'aria-label="cohort domain composition">'
    ]
    x = 0.0
    for k, v in items:
        w = (v / total) * width
        parts.append(
            f'<rect x="{x:.2f}" y="0" width="{w:.2f}" height="{height}" '
            f'fill="{cmap[k]}"><title>{_esc(k)}: {v / total * 100:.1f}%</title></rect>'
        )
        x += w
    parts.append('</svg>')
    return "".join(parts)


def domain_legend_html(domains: Sequence[str]) -> str:
    """Inline coloured-square legend for domain colours."""
    cmap = domain_colour_map(domains)
    items = "".join(
        f'<span class="lg-item"><span class="lg-sw" '
        f'style="background:{c}"></span>{_esc(d)}</span>'
        for d, c in cmap.items()
    )
    return f'<div class="legend">{items}</div>'


# ---------------------------------------------------------------------------
# Scatter (used for prevalence-vs-abundance, PCoA)
# ---------------------------------------------------------------------------


def scatter_svg(
    points: Sequence[Mapping[str, object]],
    *,
    title: str,
    x_label: str,
    y_label: str,
    width: int = 720,
    height: int = 440,
    x_log: bool = False,
    y_log: bool = False,
    quadrant_lines: tuple[float, float] | None = None,
    quadrant_labels: tuple[str, str, str, str] | None = None,
    legend_domains: Sequence[str] | None = None,
) -> str:
    """Render a colour-coded scatter plot.

    Each point dict accepts: ``x``, ``y``, ``size`` (radius pixels –
    optional, default 4), ``colour`` (hex), ``label`` (tooltip),
    ``domain`` (used for legend colour-fallback when ``colour`` is
    omitted).  Quadrant lines (vertical x, horizontal y in *plot*
    coordinates, after log-transform) are optional.
    """
    pad_left, pad_right = 70, 18
    pad_top, pad_bottom = 30, 56
    plot_w = width - pad_left - pad_right
    plot_h = height - pad_top - pad_bottom
    if not points:
        return (
            f'<div class="figure"><p class="fig-title">{_esc(title)}</p>'
            f'<p><em>No data to plot.</em></p></div>'
        )

    cmap: dict[str, str] = {}
    if legend_domains is not None:
        cmap = domain_colour_map(legend_domains)

    def _xform(x: float) -> float:
        return _safe_log10(x) if x_log else float(x)

    def _yform(y: float) -> float:
        return _safe_log10(y) if y_log else float(y)

    xs = [_xform(float(p["x"])) for p in points]
    ys = [_yform(float(p["y"])) for p in points]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    if x_max == x_min:
        x_max = x_min + 1.0
    if y_max == y_min:
        y_max = y_min + 1.0

    def _sx(x: float) -> float:
        return pad_left + (x - x_min) / (x_max - x_min) * plot_w

    def _sy(y: float) -> float:
        return pad_top + (1 - (y - y_min) / (y_max - y_min)) * plot_h

    parts: list[str] = []
    parts.append(f'<div class="figure">')
    parts.append(f'<p class="fig-title">{_esc(title)}</p>')
    parts.append(
        f'<svg viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img" '
        f'aria-label="{_esc(title)}">'
    )

    # Plot area frame
    parts.append(
        f'<rect x="{pad_left}" y="{pad_top}" width="{plot_w}" '
        f'height="{plot_h}" fill="#fafafa" stroke="#ccc"/>'
    )

    # Quadrant guide lines + labels (in transformed coords).
    if quadrant_lines is not None:
        qx, qy = quadrant_lines
        if x_min <= qx <= x_max:
            x = _sx(qx)
            parts.append(
                f'<line x1="{x:.1f}" y1="{pad_top}" x2="{x:.1f}" '
                f'y2="{pad_top + plot_h}" stroke="#888" stroke-dasharray="4 3"/>'
            )
        if y_min <= qy <= y_max:
            y = _sy(qy)
            parts.append(
                f'<line x1="{pad_left}" y1="{y:.1f}" x2="{pad_left + plot_w}" '
                f'y2="{y:.1f}" stroke="#888" stroke-dasharray="4 3"/>'
            )
        if quadrant_labels is not None:
            tl, tr, bl, br = quadrant_labels
            for txt, gx, gy, anchor in (
                (tl, pad_left + 6, pad_top + 14, "start"),
                (tr, pad_left + plot_w - 6, pad_top + 14, "end"),
                (bl, pad_left + 6, pad_top + plot_h - 6, "start"),
                (br, pad_left + plot_w - 6, pad_top + plot_h - 6, "end"),
            ):
                parts.append(
                    f'<text x="{gx}" y="{gy}" font-size="11" fill="#666" '
                    f'text-anchor="{anchor}" font-style="italic">{_esc(txt)}</text>'
                )

    # Tick marks (5 on each axis)
    def _fmt_tick(v: float, log: bool) -> str:
        if log:
            return f"10^{v:.1f}" if abs(v) > 2 else f"{10 ** v:.3g}"
        if abs(v) >= 1000 or (0 < abs(v) < 0.01):
            return f"{v:.2g}"
        return f"{v:.2f}"

    for k in range(5):
        frac = k / 4
        xv = x_min + frac * (x_max - x_min)
        yv = y_min + frac * (y_max - y_min)
        x_pix = _sx(xv)
        y_pix = _sy(yv)
        parts.append(
            f'<line x1="{x_pix:.1f}" y1="{pad_top + plot_h}" '
            f'x2="{x_pix:.1f}" y2="{pad_top + plot_h + 4}" stroke="#444"/>'
            f'<text x="{x_pix:.1f}" y="{pad_top + plot_h + 16}" '
            f'font-size="10" text-anchor="middle">{_fmt_tick(xv, x_log)}</text>'
        )
        parts.append(
            f'<line x1="{pad_left - 4}" y1="{y_pix:.1f}" '
            f'x2="{pad_left}" y2="{y_pix:.1f}" stroke="#444"/>'
            f'<text x="{pad_left - 6}" y="{y_pix + 3:.1f}" '
            f'font-size="10" text-anchor="end">{_fmt_tick(yv, y_log)}</text>'
        )

    # Axis labels
    parts.append(
        f'<text x="{pad_left + plot_w / 2:.1f}" y="{height - 18}" '
        f'font-size="12" text-anchor="middle">{_esc(x_label)}</text>'
    )
    parts.append(
        f'<text x="14" y="{pad_top + plot_h / 2:.1f}" font-size="12" '
        f'text-anchor="middle" '
        f'transform="rotate(-90 14 {pad_top + plot_h / 2:.1f})">{_esc(y_label)}</text>'
    )

    # Points
    for i, p in enumerate(points):
        x = _sx(_xform(float(p["x"])))
        y = _sy(_yform(float(p["y"])))
        r = float(p.get("size", 4) or 4)
        colour = p.get("colour")
        if not colour:
            domain = str(p.get("domain") or "")
            colour = cmap.get(domain) or _colour_for(i)
        label = p.get("label", "")
        parts.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{colour}" '
            f'fill-opacity="0.7" stroke="#222" stroke-width="0.4">'
            f'<title>{_esc(label)}</title></circle>'
        )

    parts.append('</svg>')
    if legend_domains is not None:
        parts.append(domain_legend_html(legend_domains))
    parts.append('</div>')
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Boxplot
# ---------------------------------------------------------------------------


def boxplot_svg(
    groups: Sequence[tuple[str, Sequence[float]]],
    *,
    title: str,
    y_label: str,
    width: int = 720,
    height: int = 360,
    log_y: bool = False,
) -> str:
    """Horizontal box-and-whisker, one row per group.

    Each box shows q25/q50/q75 + whiskers at q05/q95 + a coloured
    median tick.  Used for cohort-wide distribution figures (§3.5).
    """
    pad_left, pad_right = 140, 24
    pad_top, pad_bottom = 30, 50
    plot_w = width - pad_left - pad_right
    n = max(len(groups), 1)
    plot_h = height - pad_top - pad_bottom
    row_h = plot_h / n

    # Compute global x range
    all_vals: list[float] = []
    for _, vs in groups:
        for v in vs:
            if v is None:
                continue
            if isinstance(v, float) and math.isnan(v):
                continue
            if log_y and v <= 0:
                continue
            all_vals.append(float(v))
    if not all_vals:
        return (
            f'<div class="figure"><p class="fig-title">{_esc(title)}</p>'
            f'<p><em>No data to plot.</em></p></div>'
        )
    if log_y:
        x_min = math.log10(min(all_vals))
        x_max = math.log10(max(all_vals))
    else:
        x_min = min(all_vals)
        x_max = max(all_vals)
    if x_max == x_min:
        x_max = x_min + 1

    def _sx(v: float) -> float:
        if log_y:
            v = math.log10(v) if v > 0 else x_min
        return pad_left + (v - x_min) / (x_max - x_min) * plot_w

    def _q(sorted_vals: Sequence[float], q: float) -> float:
        if not sorted_vals:
            return 0.0
        n = len(sorted_vals)
        if n == 1:
            return float(sorted_vals[0])
        pos = q * (n - 1)
        lo = int(math.floor(pos))
        hi = int(math.ceil(pos))
        if lo == hi:
            return float(sorted_vals[lo])
        f = pos - lo
        return float(sorted_vals[lo]) * (1 - f) + float(sorted_vals[hi]) * f

    parts = [
        f'<div class="figure">',
        f'<p class="fig-title">{_esc(title)}</p>',
        f'<svg viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img" '
        f'aria-label="{_esc(title)}">',
        f'<rect x="{pad_left}" y="{pad_top}" width="{plot_w}" '
        f'height="{plot_h}" fill="#fafafa" stroke="#ccc"/>',
    ]

    cmap = domain_colour_map(g for g, _ in groups)
    for i, (gname, vs) in enumerate(groups):
        clean = sorted(
            float(v) for v in vs
            if v is not None and not (isinstance(v, float) and math.isnan(v))
            and (not log_y or v > 0)
        )
        y = pad_top + i * row_h + row_h / 2
        bh = row_h * 0.4
        parts.append(
            f'<text x="{pad_left - 8}" y="{y + 4:.1f}" font-size="12" '
            f'text-anchor="end">{_esc(gname)}</text>'
        )
        if not clean:
            continue
        q05 = _q(clean, 0.05)
        q25 = _q(clean, 0.25)
        q50 = _q(clean, 0.50)
        q75 = _q(clean, 0.75)
        q95 = _q(clean, 0.95)
        x05, x25, x50, x75, x95 = (_sx(q05), _sx(q25), _sx(q50), _sx(q75), _sx(q95))
        colour = cmap[gname]
        # Whiskers
        parts.append(
            f'<line x1="{x05:.2f}" y1="{y:.2f}" x2="{x95:.2f}" y2="{y:.2f}" '
            f'stroke="#444"/>'
        )
        parts.append(
            f'<line x1="{x05:.2f}" y1="{y - bh / 2:.2f}" x2="{x05:.2f}" '
            f'y2="{y + bh / 2:.2f}" stroke="#444"/>'
            f'<line x1="{x95:.2f}" y1="{y - bh / 2:.2f}" x2="{x95:.2f}" '
            f'y2="{y + bh / 2:.2f}" stroke="#444"/>'
        )
        # Box
        parts.append(
            f'<rect x="{x25:.2f}" y="{y - bh:.2f}" '
            f'width="{max(x75 - x25, 0.5):.2f}" height="{bh * 2:.2f}" '
            f'fill="{colour}" fill-opacity="0.5" stroke="#222"/>'
        )
        # Median
        parts.append(
            f'<line x1="{x50:.2f}" y1="{y - bh:.2f}" x2="{x50:.2f}" '
            f'y2="{y + bh:.2f}" stroke="#000" stroke-width="2"/>'
        )

    # X-axis ticks
    for k in range(5):
        frac = k / 4
        xv = x_min + frac * (x_max - x_min)
        x_pix = pad_left + frac * plot_w
        if log_y:
            tick = f"{10 ** xv:.3g}"
        else:
            tick = f"{xv:.3g}"
        parts.append(
            f'<line x1="{x_pix:.2f}" y1="{pad_top + plot_h:.2f}" '
            f'x2="{x_pix:.2f}" y2="{pad_top + plot_h + 4:.2f}" stroke="#444"/>'
            f'<text x="{x_pix:.2f}" y="{pad_top + plot_h + 16:.2f}" '
            f'font-size="10" text-anchor="middle">{tick}</text>'
        )
    parts.append(
        f'<text x="{pad_left + plot_w / 2:.1f}" y="{height - 14}" '
        f'font-size="12" text-anchor="middle">{_esc(y_label)}</text>'
    )
    parts.append('</svg></div>')
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Histogram
# ---------------------------------------------------------------------------


def histogram_svg(
    edges: Sequence[float],
    counts: Sequence[int],
    *,
    title: str,
    x_label: str,
    y_label: str = "samples",
    width: int = 720,
    height: int = 280,
    vline: float | None = None,
    log_x: bool = False,
) -> str:
    if not counts or max(counts) <= 0:
        return (
            f'<div class="figure"><p class="fig-title">{_esc(title)}</p>'
            f'<p><em>No samples to plot.</em></p></div>'
        )
    pad_left, pad_right = 60, 18
    pad_top, pad_bottom = 24, 50
    plot_w = width - pad_left - pad_right
    plot_h = height - pad_top - pad_bottom
    n = len(counts)
    cmax = max(counts)

    def _sx(x: float) -> float:
        v = math.log10(x) if log_x and x > 0 else x
        lo = math.log10(edges[0]) if log_x and edges[0] > 0 else edges[0]
        hi = math.log10(edges[-1]) if log_x and edges[-1] > 0 else edges[-1]
        if hi == lo:
            return pad_left
        return pad_left + (v - lo) / (hi - lo) * plot_w

    parts = [
        f'<div class="figure">',
        f'<p class="fig-title">{_esc(title)}</p>',
        f'<svg viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img">',
        f'<rect x="{pad_left}" y="{pad_top}" width="{plot_w}" height="{plot_h}" '
        f'fill="#fafafa" stroke="#ccc"/>',
    ]
    for i in range(n):
        c = counts[i]
        if c <= 0:
            continue
        x0 = _sx(edges[i])
        x1 = _sx(edges[i + 1])
        h = (c / cmax) * plot_h
        y = pad_top + plot_h - h
        parts.append(
            f'<rect x="{x0:.2f}" y="{y:.2f}" '
            f'width="{max(x1 - x0 - 0.5, 0.5):.2f}" height="{h:.2f}" '
            f'fill="#56B4E9" stroke="#222" stroke-width="0.4">'
            f'<title>{c} samples</title></rect>'
        )
    if vline is not None:
        x = _sx(vline)
        if pad_left - 1 <= x <= pad_left + plot_w + 1:
            parts.append(
                f'<line x1="{x:.2f}" y1="{pad_top}" x2="{x:.2f}" '
                f'y2="{pad_top + plot_h}" stroke="#D55E00" stroke-width="2" '
                f'stroke-dasharray="6 4"/>'
                f'<text x="{x:.2f}" y="{pad_top - 6}" font-size="11" '
                f'text-anchor="middle" fill="#D55E00">threshold</text>'
            )
    # X axis ticks
    for k in range(5):
        frac = k / 4
        xv = edges[0] + frac * (edges[-1] - edges[0])
        if log_x and edges[0] > 0:
            xv = 10 ** (math.log10(edges[0]) + frac * (math.log10(edges[-1]) - math.log10(edges[0])))
        x_pix = pad_left + frac * plot_w
        parts.append(
            f'<line x1="{x_pix:.2f}" y1="{pad_top + plot_h:.2f}" '
            f'x2="{x_pix:.2f}" y2="{pad_top + plot_h + 4:.2f}" stroke="#444"/>'
            f'<text x="{x_pix:.2f}" y="{pad_top + plot_h + 16:.2f}" '
            f'font-size="10" text-anchor="middle">{xv:.3g}</text>'
        )
    # Y axis label + max tick
    parts.append(
        f'<text x="{pad_left - 6}" y="{pad_top + 10}" font-size="10" '
        f'text-anchor="end">{cmax}</text>'
    )
    parts.append(
        f'<text x="{pad_left - 6}" y="{pad_top + plot_h}" font-size="10" '
        f'text-anchor="end">0</text>'
    )
    parts.append(
        f'<text x="{pad_left + plot_w / 2:.1f}" y="{height - 14}" '
        f'font-size="12" text-anchor="middle">{_esc(x_label)}</text>'
    )
    parts.append(
        f'<text x="14" y="{pad_top + plot_h / 2:.1f}" font-size="12" '
        f'text-anchor="middle" transform="rotate(-90 14 '
        f'{pad_top + plot_h / 2:.1f})">{_esc(y_label)}</text>'
    )
    parts.append('</svg></div>')
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Rank-abundance / Whittaker
# ---------------------------------------------------------------------------


def rank_abundance_svg(
    series: Sequence[Mapping[str, object]],
    *,
    title: str,
    width: int = 720,
    height: int = 320,
    max_rank: int = 200,
) -> str:
    """Log-rank vs log-burden curve, coloured by domain."""
    if not series:
        return (
            f'<div class="figure"><p class="fig-title">{_esc(title)}</p>'
            f'<p><em>No data.</em></p></div>'
        )
    series = list(series)[:max_rank]
    domains = sorted({str(p.get("domain", "Unannotated")) for p in series})
    cmap = domain_colour_map(domains)
    points = [
        {
            "x": float(p["rank"]),
            "y": float(p["value"]),
            "domain": p.get("domain", "Unannotated"),
            "label": f"#{p['rank']} {p['name']} ({p.get('domain')})",
            "size": 3.2,
            "colour": cmap.get(str(p.get("domain", "Unannotated"))),
        }
        for p in series
    ]
    return scatter_svg(
        points,
        title=title,
        x_label="rank (log10)",
        y_label="cohort burden (log10)",
        width=width,
        height=height,
        x_log=True,
        y_log=True,
        legend_domains=domains,
    )


# ---------------------------------------------------------------------------
# Heatmap with dendrogram
# ---------------------------------------------------------------------------


def heatmap_with_dendrogram_svg(
    values: Sequence[Sequence[float]],
    *,
    row_labels: Sequence[str],
    col_labels: Sequence[str],
    row_order: Sequence[int],
    col_order: Sequence[int],
    title: str,
    row_axis_label: str = "species",
    col_axis_label: str = "samples",
    width: int = 980,
    cell_min_h: float = 1.0,
    cell_min_w: float = 1.0,
    show_col_labels: bool = False,
    show_row_labels: bool = True,
) -> str:
    """Render a coloured heatmap with the rows and columns reordered.

    The dendrogram shapes themselves are *not* rendered (a 3K-leaf
    dendrogram is unreadable and dwarfs the heatmap on a typical
    monitor).  Instead we colour-code rows by domain through the
    palette and rely on the precomputed ``row_order`` /
    ``col_order`` to expose cluster structure visually.
    """
    n_rows = len(values)
    n_cols = len(values[0]) if values else 0
    if n_rows == 0 or n_cols == 0:
        return (
            f'<div class="figure"><p class="fig-title">{_esc(title)}</p>'
            f'<p><em>No data.</em></p></div>'
        )

    pad_left = 220 if show_row_labels else 60
    pad_top = 40
    pad_bottom = 80 if show_col_labels else 30
    pad_right = 30
    cell_w = max(cell_min_w, (width - pad_left - pad_right) / n_cols)
    cell_h = max(cell_min_h, 12.0 if n_rows < 30 else (8.0 if n_rows < 80 else 4.0))
    plot_w = cell_w * n_cols
    plot_h = cell_h * n_rows
    height = pad_top + plot_h + pad_bottom

    # Find global max for colour scaling (log1p) – ignore None.
    vmax = 0.0
    for row in values:
        for v in row:
            if v is None:
                continue
            v = float(v)
            if math.isnan(v):
                continue
            if v > vmax:
                vmax = v
    if vmax <= 0:
        vmax = 1.0
    log_vmax = math.log1p(vmax)

    def _cell_colour(v: float | None) -> str:
        if v is None or v <= 0:
            return "#ffffff"
        v = float(v)
        if math.isnan(v):
            return "#ffffff"
        t = math.log1p(v) / log_vmax
        # Cool-warm gradient: white -> deep blue.
        # Linear interp between #f7fbff and #08306b.
        r = int(247 + (8 - 247) * t)
        g = int(251 + (48 - 251) * t)
        b = int(255 + (107 - 255) * t)
        return f"#{r:02x}{g:02x}{b:02x}"

    parts: list[str] = []
    parts.append(f'<div class="figure">')
    parts.append(f'<p class="fig-title">{_esc(title)}</p>')
    parts.append(
        f'<svg viewBox="0 0 {width} {height:.0f}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img" '
        f'aria-label="{_esc(title)}">'
    )

    parts.append(
        f'<rect x="{pad_left}" y="{pad_top}" width="{plot_w:.1f}" '
        f'height="{plot_h:.1f}" fill="#fff" stroke="#ccc"/>'
    )

    for ri, r in enumerate(row_order):
        for ci, c in enumerate(col_order):
            v = values[r][c] if r < len(values) and c < len(values[r]) else None
            colour = _cell_colour(v)
            x = pad_left + ci * cell_w
            y = pad_top + ri * cell_h
            # Skip transparent rectangle generation for empty cells when
            # cells are tiny – keeps the SVG manageably sized.
            if v is None or v <= 0:
                continue
            parts.append(
                f'<rect x="{x:.2f}" y="{y:.2f}" width="{cell_w:.2f}" '
                f'height="{cell_h:.2f}" fill="{colour}">'
                f'<title>{_esc(row_labels[r])} / {_esc(col_labels[c])}: '
                f'{float(v):.3g}</title></rect>'
            )

        if show_row_labels:
            label = row_labels[r] if r < len(row_labels) else ""
            parts.append(
                f'<text x="{pad_left - 6}" y="{pad_top + (ri + 0.5) * cell_h + 3:.2f}" '
                f'font-size="11" text-anchor="end">{_esc(label)}</text>'
            )

    if show_col_labels:
        for ci, c in enumerate(col_order):
            label = col_labels[c] if c < len(col_labels) else ""
            x = pad_left + (ci + 0.5) * cell_w
            y = pad_top + plot_h + 6
            parts.append(
                f'<text x="{x:.2f}" y="{y:.2f}" font-size="9" '
                f'transform="rotate(-60 {x:.2f} {y:.2f})">{_esc(label)}</text>'
            )

    parts.append(
        f'<text x="{pad_left + plot_w / 2:.1f}" y="{height - 8:.0f}" '
        f'font-size="11" text-anchor="middle">{_esc(col_axis_label)} '
        f'(n = {n_cols}; ordered by hierarchical clustering)</text>'
    )
    parts.append('</svg></div>')
    return "\n".join(parts)
