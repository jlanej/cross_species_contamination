"""Cohort / metagenomic-summary statistics for the CSC HTML report.

Pure-stdlib helpers that summarise a cohort of WGS samples around the
**species** axis instead of the per-sample axis.  Output of these
functions feeds the new §3 *Cohort species landscape* section
(species summary table, prevalence-vs-abundance map, rank-abundance,
core/accessory/rare partition, hierarchical-clustering heatmap, and
PCoA ordination).

Design notes
------------
* No numpy / scipy.  Everything operates on ``dict[int, dict[str, float]]``
  matrices as parsed by :func:`csc.report.report._parse_matrix`.
* All summaries use **robust** statistics: median, IQR, percentiles,
  CV-of-positives.  Bare means of long-tailed contamination
  distributions are misleading and deliberately avoided.
* Compositional (CPM) and absolute (ppm of total sequenced reads)
  denominators are summarised side by side, never silently mixed.
* Hierarchical clustering is implemented via the
  Lance–Williams recurrence so any of {average, single, ward} works
  on the same generic driver.

AI assistance acknowledgment: this module was developed with AI
assistance.  Best practices in the bioinformatics field should always
take precedence over specific implementation details.
"""

from __future__ import annotations

import math
import statistics
from typing import Any, Iterable, Mapping, Sequence

# A "matrix" in this module mirrors :class:`csc.report.report.Matrix`
# (we accept the parsed dataclass directly to avoid duplication).


# ---------------------------------------------------------------------------
# Robust descriptive statistics
# ---------------------------------------------------------------------------


def _percentile(sorted_vals: Sequence[float], q: float) -> float:
    """Linear-interpolation percentile on an already-sorted sequence.

    ``q`` is a fraction in ``[0, 1]``.  Returns ``float('nan')`` for
    an empty input.
    """
    n = len(sorted_vals)
    if n == 0:
        return float("nan")
    if n == 1:
        return float(sorted_vals[0])
    pos = q * (n - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return float(sorted_vals[lo])
    frac = pos - lo
    return float(sorted_vals[lo]) * (1 - frac) + float(sorted_vals[hi]) * frac


def cohort_quantiles(values: Iterable[float]) -> dict[str, float]:
    """Return q5/q25/q50/q75/q95 + min/max for a sequence of values.

    NaNs and ``None`` are dropped.  An empty input returns NaNs.
    """
    cleaned: list[float] = []
    for v in values:
        if v is None:
            continue
        if isinstance(v, float) and math.isnan(v):
            continue
        cleaned.append(float(v))
    cleaned.sort()
    return {
        "n": len(cleaned),
        "min": cleaned[0] if cleaned else float("nan"),
        "q05": _percentile(cleaned, 0.05),
        "q25": _percentile(cleaned, 0.25),
        "median": _percentile(cleaned, 0.50),
        "q75": _percentile(cleaned, 0.75),
        "q95": _percentile(cleaned, 0.95),
        "max": cleaned[-1] if cleaned else float("nan"),
    }


def histogram(
    values: Iterable[float],
    n_bins: int = 30,
    *,
    log: bool = False,
) -> tuple[list[float], list[int]]:
    """Build a histogram of ``values``.

    Returns ``(edges, counts)`` where ``len(edges) == n_bins + 1``.
    When ``log=True`` the bins are log10-spaced over positive values
    only (zeros / non-positive values are dropped); the caller is
    responsible for labelling axes appropriately.
    """
    cleaned = [float(v) for v in values if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if log:
        cleaned = [v for v in cleaned if v > 0]
    if not cleaned:
        return ([0.0, 1.0], [0] * n_bins)
    if log:
        lo = math.log10(min(cleaned))
        hi = math.log10(max(cleaned))
    else:
        lo = min(cleaned)
        hi = max(cleaned)
    if hi <= lo:
        hi = lo + 1.0
    width = (hi - lo) / n_bins
    edges = [lo + i * width for i in range(n_bins + 1)]
    counts = [0] * n_bins
    for v in cleaned:
        x = math.log10(v) if log else v
        idx = int((x - lo) / width)
        if idx == n_bins:
            idx = n_bins - 1
        if 0 <= idx < n_bins:
            counts[idx] += 1
    if log:
        edges = [10 ** e for e in edges]
    return (edges, counts)


# ---------------------------------------------------------------------------
# Species summary table (the new §3.1 centrepiece)
# ---------------------------------------------------------------------------


def species_summary_rows(
    matrix_raw: Any,
    matrix_cpm: Any,
    matrix_abs: Any | None,
    *,
    flagged_taxa_primary: Mapping[int, int] | None = None,
    flagged_taxa_abs: Mapping[int, int] | None = None,
    min_reads_for_prevalence: int = 5,
    human_tax_ids: tuple[int, ...] = (9606,),
    sparkline_bins: int = 12,
) -> list[dict[str, Any]]:
    """Compute one summary row per taxon (excluding human).

    Each row contains:

    ``tax_id``, ``name``, ``domain``, ``prevalence`` (0–1),
    ``prevalence_at_min`` (0–1), ``positives_n``,
    ``cpm_median_pos``, ``cpm_iqr_pos``, ``cpm_p95_pos``,
    ``abs_median_pos``, ``abs_iqr_pos``, ``abs_p95_pos`` (NaN when no
    abs matrix), ``cohort_raw_total``, ``cohort_burden_ppm`` (Σ ppm
    when abs available, else NaN), ``n_flagged_primary``,
    ``n_flagged_abs``, ``cv_of_positives`` (MAD/median; robust),
    ``sparkline`` (list of int counts in ``sparkline_bins`` log-spaced
    bins of CPM positives).
    """
    flagged_taxa_primary = flagged_taxa_primary or {}
    flagged_taxa_abs = flagged_taxa_abs or {}

    n_samples = max(len(matrix_raw.sample_ids), 1)
    rows: list[dict[str, Any]] = []

    for tid in matrix_raw.tax_ids:
        if tid in human_tax_ids:
            continue
        # Raw positives
        raw_pos = 0
        raw_pos_min = 0
        raw_total = 0.0
        for sid in matrix_raw.sample_ids:
            v = matrix_raw.values.get(tid, {}).get(sid)
            if v is None or v <= 0:
                continue
            raw_pos += 1
            raw_total += v
            if v >= min_reads_for_prevalence:
                raw_pos_min += 1

        # CPM positives
        cpm_pos: list[float] = []
        for sid in matrix_cpm.sample_ids:
            v = matrix_cpm.values.get(tid, {}).get(sid)
            if v is not None and v > 0:
                cpm_pos.append(float(v))
        cpm_pos.sort()
        cpm_q = _summary_pos(cpm_pos)

        # ABS positives (when available)
        if matrix_abs is not None:
            abs_pos: list[float] = []
            abs_total_ppm = 0.0
            for sid in matrix_abs.sample_ids:
                v = matrix_abs.values.get(tid, {}).get(sid)
                if v is not None and v > 0:
                    abs_pos.append(float(v))
                    abs_total_ppm += float(v)
            abs_pos.sort()
            abs_q = _summary_pos(abs_pos)
            cohort_burden_ppm = abs_total_ppm
        else:
            abs_q = {"n": 0, "median": float("nan"), "iqr": float("nan"),
                     "p95": float("nan")}
            cohort_burden_ppm = float("nan")

        # Robust CV (MAD / median) over CPM positives
        if len(cpm_pos) >= 2 and cpm_q["median"] > 0:
            med = cpm_q["median"]
            mad = statistics.median(abs(v - med) for v in cpm_pos)
            cv_robust = mad / med
        else:
            cv_robust = float("nan")

        # Sparkline = log-binned histogram of CPM positives
        if cpm_pos:
            _, sparkline = histogram(cpm_pos, n_bins=sparkline_bins, log=True)
        else:
            sparkline = [0] * sparkline_bins

        rows.append({
            "tax_id": tid,
            "name": matrix_raw.tax_names.get(tid, str(tid)),
            "domain": matrix_raw.tax_domains.get(tid, "Unannotated"),
            "prevalence": raw_pos / n_samples,
            "prevalence_at_min": raw_pos_min / n_samples,
            "positives_n": raw_pos,
            "min_reads_for_prevalence": min_reads_for_prevalence,
            "cpm_median_pos": cpm_q["median"],
            "cpm_iqr_pos": cpm_q["iqr"],
            "cpm_p95_pos": cpm_q["p95"],
            "abs_median_pos": abs_q["median"],
            "abs_iqr_pos": abs_q["iqr"],
            "abs_p95_pos": abs_q["p95"],
            "cohort_raw_total": raw_total,
            "cohort_burden_ppm": cohort_burden_ppm,
            "n_flagged_primary": int(flagged_taxa_primary.get(tid, 0)),
            "n_flagged_abs": int(flagged_taxa_abs.get(tid, 0)),
            "cv_robust": cv_robust,
            "sparkline": sparkline,
        })

    # Default sort = cohort burden (ppm) desc, falling back to raw_total.
    def _sort_key(r: dict[str, Any]) -> tuple[float, float]:
        bur = r["cohort_burden_ppm"]
        if isinstance(bur, float) and math.isnan(bur):
            bur = -1.0
        return (-float(bur), -float(r["cohort_raw_total"]))

    rows.sort(key=_sort_key)
    return rows


def _summary_pos(sorted_vals: Sequence[float]) -> dict[str, float]:
    if not sorted_vals:
        return {"n": 0, "median": float("nan"), "iqr": float("nan"),
                "p95": float("nan")}
    return {
        "n": len(sorted_vals),
        "median": _percentile(sorted_vals, 0.5),
        "iqr": _percentile(sorted_vals, 0.75) - _percentile(sorted_vals, 0.25),
        "p95": _percentile(sorted_vals, 0.95),
    }


# ---------------------------------------------------------------------------
# Prevalence-band partition (core / accessory / rare)
# ---------------------------------------------------------------------------


def prevalence_partition(
    rows: Sequence[Mapping[str, Any]],
    *,
    core_threshold: float = 0.5,
    rare_threshold: float = 0.1,
    top_n: int = 10,
) -> dict[str, Any]:
    """Split species rows into core / accessory / rare prevalence bands.

    A species is *core* if ``prevalence >= core_threshold``,
    *accessory* if ``rare_threshold <= prevalence < core_threshold``,
    else *rare*.  Each band is sorted by ``cohort_burden_ppm``
    (falling back to ``cohort_raw_total``); only the top-``top_n``
    rows of each band are returned, plus the band counts.
    """
    core: list[Mapping[str, Any]] = []
    accessory: list[Mapping[str, Any]] = []
    rare: list[Mapping[str, Any]] = []
    for r in rows:
        p = float(r.get("prevalence", 0.0) or 0.0)
        if p >= core_threshold:
            core.append(r)
        elif p >= rare_threshold:
            accessory.append(r)
        else:
            rare.append(r)

    def _key(r: Mapping[str, Any]) -> tuple[float, float]:
        bur = r.get("cohort_burden_ppm")
        if bur is None or (isinstance(bur, float) and math.isnan(bur)):
            bur = -1.0
        return (-float(bur), -float(r.get("cohort_raw_total", 0.0) or 0.0))

    return {
        "core_threshold": core_threshold,
        "rare_threshold": rare_threshold,
        "core_count": len(core),
        "accessory_count": len(accessory),
        "rare_count": len(rare),
        "core_top": sorted(core, key=_key)[:top_n],
        "accessory_top": sorted(accessory, key=_key)[:top_n],
        "rare_top": sorted(rare, key=_key)[:top_n],
    }


def rank_abundance(
    rows: Sequence[Mapping[str, Any]],
    *,
    metric: str = "cohort_burden_ppm",
) -> list[dict[str, Any]]:
    """Return a rank-abundance series.

    Each entry is ``{rank, name, tax_id, domain, value}``, sorted
    descending by ``metric`` (defaults to absolute cohort burden;
    falls back to raw read total when the abs matrix is unavailable).
    """
    series: list[tuple[float, Mapping[str, Any]]] = []
    for r in rows:
        v = r.get(metric)
        if v is None or (isinstance(v, float) and math.isnan(v)):
            v = float(r.get("cohort_raw_total") or 0.0)
        v = float(v)
        if v <= 0:
            continue
        series.append((v, r))
    series.sort(key=lambda kv: -kv[0])
    out: list[dict[str, Any]] = []
    for i, (v, r) in enumerate(series, 1):
        out.append({
            "rank": i,
            "name": r.get("name"),
            "tax_id": r.get("tax_id"),
            "domain": r.get("domain", "Unannotated"),
            "value": v,
        })
    return out


# ---------------------------------------------------------------------------
# Pairwise distances
# ---------------------------------------------------------------------------


def bray_curtis_matrix(
    matrix: Any,
    sample_ids: Sequence[str],
    tax_ids: Sequence[int] | None = None,
) -> list[list[float]]:
    """Pairwise Bray–Curtis dissimilarity over the given samples.

    Restricts to ``tax_ids`` when supplied (recommended: top-K species
    by cohort burden, to keep the cost O(K · n²)).  ``None`` cells in
    the matrix are treated as zero.
    """
    if tax_ids is None:
        tax_ids = list(matrix.tax_ids)
    # Pre-extract per-sample vectors as plain lists (positional) for speed
    vecs: list[list[float]] = []
    sums: list[float] = []
    for sid in sample_ids:
        v = []
        s = 0.0
        for tid in tax_ids:
            x = matrix.values.get(tid, {}).get(sid)
            if x is None or (isinstance(x, float) and math.isnan(x)):
                x = 0.0
            else:
                x = float(x)
                if x < 0:
                    x = 0.0
            v.append(x)
            s += x
        vecs.append(v)
        sums.append(s)

    n = len(sample_ids)
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        si = sums[i]
        vi = vecs[i]
        for j in range(i + 1, n):
            sj = sums[j]
            denom = si + sj
            if denom <= 0:
                d = 0.0
            else:
                vj = vecs[j]
                num = 0.0
                for a, b in zip(vi, vj):
                    num += abs(a - b)
                d = num / denom
            D[i][j] = D[j][i] = d
    return D


# ---------------------------------------------------------------------------
# Hierarchical clustering (Lance–Williams)
# ---------------------------------------------------------------------------


def hclust(
    distance_matrix: Sequence[Sequence[float]],
    *,
    method: str = "average",
) -> dict[str, Any]:
    """Agglomerative hierarchical clustering via Lance–Williams.

    Supports ``"average"`` (UPGMA), ``"single"``, and ``"ward"``.
    Returns a dict with:

    * ``order`` – leaf ordering produced by post-order traversal
      (suitable for plotting a heatmap).
    * ``merges`` – list of ``(i, j, distance, size)`` tuples in merge
      order (i, j refer to the cluster labels at merge time;
      ``n + k`` denotes the cluster created by the *k*-th merge).
    * ``method`` – echoed for transparency.

    For a small *n* this is O(n³), which is fine for the cohort-scale
    sub-sample we feed it (default cap = 2000).
    """
    n = len(distance_matrix)
    if n == 0:
        return {"order": [], "merges": [], "method": method}
    if n == 1:
        return {"order": [0], "merges": [], "method": method}

    # Working symmetric distance matrix as dict-of-dicts so we can
    # delete clusters as they are merged.
    D: dict[int, dict[int, float]] = {}
    for i in range(n):
        D[i] = {}
        for j in range(n):
            if i != j:
                D[i][j] = float(distance_matrix[i][j])

    sizes: dict[int, int] = {i: 1 for i in range(n)}
    children: dict[int, tuple[int, int] | None] = {i: None for i in range(n)}
    merges: list[tuple[int, int, float, int]] = []

    next_id = n
    while len(D) > 1:
        # Find the closest pair.
        best = None
        best_d = float("inf")
        for i, row in D.items():
            for j, d in row.items():
                if i < j and d < best_d:
                    best_d = d
                    best = (i, j)
        assert best is not None
        i, j = best
        ni, nj = sizes[i], sizes[j]
        new_size = ni + nj
        new_id = next_id
        next_id += 1

        # Lance–Williams update for all remaining clusters k.
        new_row: dict[int, float] = {}
        for k in list(D.keys()):
            if k == i or k == j:
                continue
            d_ik = D[i][k]
            d_jk = D[j][k]
            if method == "single":
                d_new = min(d_ik, d_jk)
            elif method == "ward":
                nk = sizes[k]
                d_ij = best_d
                # squared-distance variant of Ward (treats input D as
                # already a metric – fine for Bray–Curtis although
                # technically Ward assumes Euclidean).
                num = (
                    (ni + nk) * (d_ik ** 2)
                    + (nj + nk) * (d_jk ** 2)
                    - nk * (d_ij ** 2)
                )
                d_new = math.sqrt(num / (ni + nj + nk))
            else:  # average / UPGMA
                d_new = (ni * d_ik + nj * d_jk) / (ni + nj)
            new_row[k] = d_new
            D[k][i] = D[k][j] = 0.0  # placeholder
            del D[k][i]
            del D[k][j]
            D[k][new_id] = d_new

        del D[i]
        del D[j]
        D[new_id] = new_row

        sizes[new_id] = new_size
        children[new_id] = (i, j)
        merges.append((i, j, best_d, new_size))

    # Compute leaf order via post-order traversal of the binary tree.
    root = next(iter(D.keys()))

    order: list[int] = []
    stack = [root]
    while stack:
        node = stack.pop()
        ch = children.get(node)
        if ch is None:
            order.append(node)
        else:
            # push right then left so left is processed first
            stack.append(ch[1])
            stack.append(ch[0])

    return {"order": order, "merges": merges, "method": method}


# ---------------------------------------------------------------------------
# Classical-MDS / PCoA
# ---------------------------------------------------------------------------


def pcoa_2d(
    distance_matrix: Sequence[Sequence[float]],
    *,
    n_components: int = 2,
    power_iters: int = 200,
    tol: float = 1e-8,
) -> dict[str, Any]:
    """Classical-MDS embedding of a distance matrix into ``n_components``.

    Implements the standard double-centring transform
    ``B = -1/2 J D² J`` followed by power-iteration with deflation to
    extract the top-``n_components`` eigenpairs.  Pure stdlib – no
    numpy.

    Returns a dict with:

    * ``coords`` – list of ``[x, y, …]`` coordinates per row of
      ``distance_matrix``.
    * ``eigvals`` – the extracted eigenvalues in descending order.
    * ``explained_var`` – the same scaled to fractions of the trace
      of B (handy for figure captions).
    """
    n = len(distance_matrix)
    if n == 0:
        return {"coords": [], "eigvals": [], "explained_var": []}
    if n == 1:
        return {"coords": [[0.0] * n_components], "eigvals": [],
                "explained_var": []}

    # Build B = -1/2 (D² - row_mean - col_mean + grand_mean)
    D2 = [[float(distance_matrix[i][j]) ** 2 for j in range(n)] for i in range(n)]
    row_mean = [sum(D2[i]) / n for i in range(n)]
    col_mean = row_mean  # symmetric
    grand = sum(row_mean) / n
    B = [
        [-0.5 * (D2[i][j] - row_mean[i] - col_mean[j] + grand) for j in range(n)]
        for i in range(n)
    ]

    # Trace of B = sum of eigenvalues.
    trace = sum(B[i][i] for i in range(n))

    eigvals: list[float] = []
    eigvecs: list[list[float]] = []
    for _ in range(n_components):
        # Power iteration on the deflated B.
        v = [1.0 / math.sqrt(n)] * n
        prev_lam = 0.0
        lam = 0.0
        for _ in range(power_iters):
            w = [sum(B[i][j] * v[j] for j in range(n)) for i in range(n)]
            norm = math.sqrt(sum(x * x for x in w))
            if norm < tol:
                break
            v = [x / norm for x in w]
            # Rayleigh quotient
            lam = sum(v[i] * sum(B[i][j] * v[j] for j in range(n)) for i in range(n))
            if abs(lam - prev_lam) < tol * max(1.0, abs(lam)):
                break
            prev_lam = lam
        eigvals.append(lam)
        eigvecs.append(v)
        # Deflation: B <- B - lam v vᵀ
        for i in range(n):
            for j in range(n):
                B[i][j] -= lam * v[i] * v[j]

    # Coordinates = sqrt(max(lam, 0)) * v
    coords = [[0.0] * n_components for _ in range(n)]
    for k, (lam, v) in enumerate(zip(eigvals, eigvecs)):
        scale = math.sqrt(max(lam, 0.0))
        for i in range(n):
            coords[i][k] = scale * v[i]

    explained_var: list[float] = []
    if trace > 0:
        explained_var = [lam / trace for lam in eigvals]
    return {"coords": coords, "eigvals": eigvals, "explained_var": explained_var}


# ---------------------------------------------------------------------------
# Convenience
# ---------------------------------------------------------------------------


def domain_burden_mix(
    rows: Sequence[Mapping[str, Any]],
    *,
    metric: str = "cohort_burden_ppm",
) -> dict[str, float]:
    """Aggregate burden-by-domain (for the executive-summary mini bar)."""
    out: dict[str, float] = {}
    for r in rows:
        v = r.get(metric)
        if v is None or (isinstance(v, float) and math.isnan(v)):
            v = float(r.get("cohort_raw_total") or 0.0)
        v = float(v)
        if v <= 0:
            continue
        d = r.get("domain", "Unannotated") or "Unannotated"
        out[d] = out.get(d, 0.0) + v
    return out


def flagged_taxon_counts(
    flagged: Iterable[Mapping[str, Any]],
) -> dict[int, int]:
    """Count distinct samples flagged per taxon from a flagged_samples row stream."""
    by_tax: dict[int, set[str]] = {}
    for row in flagged:
        try:
            tid = int(row.get("tax_id") or row.get("taxid") or 0)
        except (TypeError, ValueError):
            continue
        sid = str(row.get("sample_id") or row.get("sample") or "")
        if not sid:
            continue
        by_tax.setdefault(tid, set()).add(sid)
    return {tid: len(s) for tid, s in by_tax.items()}
