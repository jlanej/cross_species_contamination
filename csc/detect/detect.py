"""Core outlier-detection logic for the CSC pipeline.

Applies statistical methods to a sample-by-taxon count matrix
(as produced by :mod:`csc.aggregate`) to identify samples whose
taxonomic profile deviates significantly from the population,
suggesting cross-species contamination.

Three complementary strategies are implemented:

* **MAD (Median Absolute Deviation)** – robust measure of dispersion
  that is insensitive to extreme outliers.  Flags samples whose
  modified Z-score (``|x − median| / (1.4826·MAD)``) exceeds a
  threshold *and* whose value is **above** the population median.
  Only the upper tail is flagged because contamination manifests as
  excess abundance, not depletion.
* **IQR (Inter-Quartile Range)** – flags values above Tukey's upper
  fence (``Q3 + k·IQR``).  As with MAD, only the upper tail is
  flagged: a sample with an unusually *low* abundance for a given
  taxon is not informative for contamination detection.
* **GMM (Gaussian Mixture Model)** – a two-component mixture model
  fitted via Expectation-Maximization.  Unlike MAD/IQR, the GMM does
  **not** assume that the majority of samples are clean; it explicitly
  models both a *background* component and a *contamination* component
  and flags samples whose posterior probability of belonging to the
  high-abundance component exceeds a configurable threshold.

An optional *kitome* / environmental-background exclusion list removes
taxa known to originate from laboratory reagents or the environment so
that they do not trigger false positives.

AI assistance acknowledgment: This module was developed with AI assistance.
Best practices in the bioinformatics field should always take precedence over
specific implementation details.
"""

from __future__ import annotations

import csv
import logging
import math
from pathlib import Path
from typing import Any, TypedDict

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public types
# ---------------------------------------------------------------------------


class FlaggedSample(TypedDict):
    """A single sample flagged as an outlier for a specific taxon."""

    sample_id: str
    tax_id: int
    taxon_name: str
    value: float
    population_median: float
    deviation: float
    method: str


class DetectionResult(TypedDict):
    """Return type for :func:`detect_outliers`."""

    flagged: list[FlaggedSample]
    summary: dict[str, Any]


# ---------------------------------------------------------------------------
# Matrix I/O
# ---------------------------------------------------------------------------


def load_matrix(path: str | Path) -> tuple[list[str], list[dict[str, Any]], dict[int, str]]:
    """Read a taxa-matrix TSV produced by :mod:`csc.aggregate`.

    Handles both the classic two-metadata-column layout
    (``tax_id, name, sample1, …``) and the domain-annotated three-column
    layout (``tax_id, name, domain, sample1, …``).

    Returns
    -------
    sample_ids : list[str]
        Column headers (sample identifiers).
    rows : list[dict]
        One dict per taxon row with keys ``tax_id`` (int) and one float
        entry per sample.  If a ``domain`` column is present the dict
        also contains a ``"domain"`` key.
    tax_names : dict[int, str]
        Mapping from tax_id to taxon name.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Matrix file not found: {path}")

    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)  # tax_id, name, [domain,] sample1, sample2, …

    # Detect whether a "domain" metadata column is present.
    has_domain = len(header) > 2 and header[2] == "domain"
    meta_cols = 3 if has_domain else 2

    if len(header) < meta_cols + 1:
        raise ValueError(
            f"Matrix must have at least {meta_cols + 1} columns "
            f"(tax_id, name, {'domain, ' if has_domain else ''}and one "
            f"sample), got {len(header)}"
        )

    sample_ids = header[meta_cols:]
    rows: list[dict[str, Any]] = []
    tax_names: dict[int, str] = {}

    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)  # skip header
        for lineno, cols in enumerate(reader, 2):
            if len(cols) < meta_cols + 1:
                logger.warning("Skipping line %d: too few columns", lineno)
                continue
            try:
                tid = int(cols[0])
            except ValueError:
                logger.warning("Skipping line %d: non-integer tax_id", lineno)
                continue
            name = cols[1]
            tax_names[tid] = name
            entry: dict[str, Any] = {"tax_id": tid}
            if has_domain:
                entry["domain"] = cols[2]
            for i, sid in enumerate(sample_ids):
                try:
                    entry[sid] = float(cols[meta_cols + i])
                except (ValueError, IndexError):
                    entry[sid] = 0.0
            rows.append(entry)

    if not rows:
        raise ValueError(f"No valid taxon rows in matrix: {path}")

    return sample_ids, rows, tax_names


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------


def _median(values: list[float]) -> float:
    """Return the median of *values*."""
    n = len(values)
    if n == 0:
        return 0.0
    s = sorted(values)
    mid = n // 2
    if n % 2 == 1:
        return s[mid]
    return (s[mid - 1] + s[mid]) / 2.0


def _mad(values: list[float]) -> float:
    """Return the Median Absolute Deviation of *values*."""
    med = _median(values)
    deviations = [abs(v - med) for v in values]
    return _median(deviations)


def _quartiles(values: list[float]) -> tuple[float, float, float]:
    """Return (Q1, median, Q3) using the exclusive method."""
    s = sorted(values)
    n = len(s)
    if n == 0:
        return (0.0, 0.0, 0.0)
    med = _median(s)
    mid = n // 2
    lower = s[:mid]
    upper = s[mid:] if n % 2 == 0 else s[mid + 1:]
    if not lower:
        lower = s[:1]
    if not upper:
        upper = s[-1:]
    q1 = _median(lower)
    q3 = _median(upper)
    return q1, med, q3


# ---------------------------------------------------------------------------
# GMM (Gaussian Mixture Model) helpers
# ---------------------------------------------------------------------------

# Default number of EM iterations and convergence tolerance.
_GMM_MAX_ITER = 200
_GMM_TOL = 1e-6

# Floor for component standard deviations to prevent numerical collapse.
_GMM_STD_FLOOR = 1e-8

# Minimum Cohen's d (separation / max sigma) required between the two
# GMM components.  A value of 3 demands an "extreme" effect size,
# preventing the model from flagging minor noise as contamination.
_GMM_MIN_COHENS_D = 3.0


def _normal_pdf(x: float, mu: float, sigma: float) -> float:
    """Univariate normal density, returning a small floor when *sigma* ≈ 0."""
    if sigma < _GMM_STD_FLOOR:
        return 1e-300
    return (
        math.exp(-0.5 * ((x - mu) / sigma) ** 2)
        / (sigma * math.sqrt(2.0 * math.pi))
    )


class _GMMResult:
    """Outcome of a two-component GMM fit on a single taxon."""

    __slots__ = (
        "converged",
        "mu_bg",
        "sigma_bg",
        "mu_contam",
        "sigma_contam",
        "weight_contam",
        "posteriors",
        "log_likelihood",
    )

    def __init__(
        self,
        converged: bool,
        mu_bg: float,
        sigma_bg: float,
        mu_contam: float,
        sigma_contam: float,
        weight_contam: float,
        posteriors: list[float],
        log_likelihood: float,
    ) -> None:
        self.converged = converged
        self.mu_bg = mu_bg
        self.sigma_bg = sigma_bg
        self.mu_contam = mu_contam
        self.sigma_contam = sigma_contam
        self.weight_contam = weight_contam
        self.posteriors = posteriors
        self.log_likelihood = log_likelihood


def _sample_variance(values: list[float], mu: float) -> float:
    """Compute sample variance of *values* around *mu*."""
    return sum((v - mu) ** 2 for v in values) / len(values)


def _fit_two_component_gmm(
    values: list[float],
    *,
    max_iter: int = _GMM_MAX_ITER,
    tol: float = _GMM_TOL,
) -> _GMMResult | None:
    """Fit a two-component univariate Gaussian mixture via EM.

    Returns *None* when the data has zero variance (all values
    identical) because no meaningful two-component model exists.

    Component labels are assigned so that *mu_bg ≤ mu_contam* (the
    component with the higher mean is the "contamination" component).

    Parameters
    ----------
    values:
        Per-sample abundances for a single taxon (length ≥ 2).
    max_iter:
        Maximum number of EM iterations.
    tol:
        Log-likelihood convergence tolerance.

    Returns
    -------
    _GMMResult or None
    """
    n = len(values)
    if n < 2:
        return None

    # Check for zero variance (all identical) → no mixture needed
    vmin = min(values)
    vmax = max(values)
    if vmax - vmin < _GMM_STD_FLOOR:
        return None

    # ---- Initialisation (K-means++ inspired) ----
    # Sort and split into lower / upper halves to seed the two components.
    s = sorted(values)
    mid = n // 2
    lower = s[:mid] if mid > 0 else s[:1]
    upper = s[mid:] if mid < n else s[-1:]

    mu_a = sum(lower) / len(lower)
    mu_b = sum(upper) / len(upper)

    var_a = _sample_variance(lower, mu_a) if len(lower) > 1 else (vmax - vmin) ** 2 / 4
    var_b = _sample_variance(upper, mu_b) if len(upper) > 1 else (vmax - vmin) ** 2 / 4

    sigma_a = max(math.sqrt(var_a), _GMM_STD_FLOOR)
    sigma_b = max(math.sqrt(var_b), _GMM_STD_FLOOR)

    pi_a = 0.5
    pi_b = 0.5

    prev_ll = -math.inf
    converged = False

    for _ in range(max_iter):
        # ---- E-step: compute responsibilities ----
        resp_a: list[float] = []
        resp_b: list[float] = []
        log_likelihood = 0.0

        for v in values:
            pa = pi_a * _normal_pdf(v, mu_a, sigma_a)
            pb = pi_b * _normal_pdf(v, mu_b, sigma_b)
            total = pa + pb
            if total < 1e-300:
                total = 1e-300
            resp_a.append(pa / total)
            resp_b.append(pb / total)
            log_likelihood += math.log(max(total, 1e-300))

        # ---- Convergence check ----
        if abs(log_likelihood - prev_ll) < tol:
            converged = True
            break
        prev_ll = log_likelihood

        # ---- M-step: update parameters ----
        n_a = sum(resp_a)
        n_b = sum(resp_b)

        if n_a < 1e-10 or n_b < 1e-10:
            # One component has collapsed – the data is effectively unimodal
            return None

        pi_a = n_a / n
        pi_b = n_b / n

        mu_a = sum(r * v for r, v in zip(resp_a, values)) / n_a
        mu_b = sum(r * v for r, v in zip(resp_b, values)) / n_b

        var_a = sum(r * (v - mu_a) ** 2 for r, v in zip(resp_a, values)) / n_a
        var_b = sum(r * (v - mu_b) ** 2 for r, v in zip(resp_b, values)) / n_b

        sigma_a = max(math.sqrt(var_a), _GMM_STD_FLOOR)
        sigma_b = max(math.sqrt(var_b), _GMM_STD_FLOOR)

    # ---- Label assignment: bg = lower-mean component ----
    if mu_a <= mu_b:
        posteriors = resp_b  # P(contamination | x_i)
        return _GMMResult(
            converged=converged,
            mu_bg=mu_a,
            sigma_bg=sigma_a,
            mu_contam=mu_b,
            sigma_contam=sigma_b,
            weight_contam=pi_b,
            posteriors=posteriors,
            log_likelihood=log_likelihood,
        )
    else:
        posteriors = resp_a
        return _GMMResult(
            converged=converged,
            mu_bg=mu_b,
            sigma_bg=sigma_b,
            mu_contam=mu_a,
            sigma_contam=sigma_a,
            weight_contam=pi_a,
            posteriors=posteriors,
            log_likelihood=log_likelihood,
        )


def _bic_favors_two_components(values: list[float], gmm: _GMMResult) -> bool:
    """Return True when BIC prefers the 2-component GMM over a single Gaussian.

    BIC = −2·LL + k·ln(n)

    * 1-component model: k = 2 (mean, variance)
    * 2-component model: k = 5 (2 means, 2 variances, 1 mixing weight)

    Lower BIC is better.  Using BIC for model selection is a
    standard approach that naturally penalises unnecessary complexity
    and prevents the GMM from over-fitting noise as a second component.
    """
    n = len(values)
    if n < 2:
        return False

    # 1-component log-likelihood
    mu = sum(values) / n
    var = sum((v - mu) ** 2 for v in values) / n
    sigma = max(math.sqrt(var), _GMM_STD_FLOOR)
    ll_one = sum(
        math.log(max(_normal_pdf(v, mu, sigma), 1e-300)) for v in values
    )

    log_n = math.log(n)
    bic_one = -2.0 * ll_one + 2.0 * log_n
    bic_two = -2.0 * gmm.log_likelihood + 5.0 * log_n

    return bic_two < bic_one


# ---------------------------------------------------------------------------
# Background / kitome filtering
# ---------------------------------------------------------------------------


def filter_kitome(
    rows: list[dict[str, Any]],
    kitome_taxa: list[int] | None,
) -> list[dict[str, Any]]:
    """Remove rows whose ``tax_id`` appears in the kitome list.

    Parameters
    ----------
    rows:
        Taxon rows from :func:`load_matrix`.
    kitome_taxa:
        NCBI taxonomy IDs known to be laboratory / environmental
        contaminants.  If *None* or empty, no filtering is performed.

    Returns
    -------
    list[dict]
        Filtered rows.
    """
    if not kitome_taxa:
        return rows

    kitome_set = set(kitome_taxa)
    before = len(rows)
    filtered = [r for r in rows if r["tax_id"] not in kitome_set]
    removed = before - len(filtered)
    if removed:
        logger.info(
            "Kitome filter removed %d / %d taxa", removed, before
        )
    return filtered


# ---------------------------------------------------------------------------
# Population-mean subtraction
# ---------------------------------------------------------------------------


def subtract_population_background(
    rows: list[dict[str, Any]],
    sample_ids: list[str],
) -> list[dict[str, Any]]:
    """Subtract the per-taxon population mean from every sample.

    After subtraction a sample sitting at the population mean will have
    a value of zero.  Negative values are preserved so that the
    distribution shape is not distorted.

    .. note::
       MAD and IQR are **translation-invariant** statistics: subtracting
       a per-taxon constant (the mean) shifts every sample by the same
       amount, so the median, MAD, quartiles and IQR all shift by that
       same amount and the set of flagged samples is unchanged.  This
       function therefore exists for **presentational** clarity (centred
       deviations are easier to interpret) rather than to alter
       detection sensitivity for MAD/IQR.  The GMM operates on the
       centred values too, but it is *not* translation-invariant; for
       GMM the centring re-anchors the background component near zero
       which can stabilise EM initialisation.

    Parameters
    ----------
    rows:
        Taxon rows from :func:`load_matrix`.
    sample_ids:
        Sample identifiers (column names).

    Returns
    -------
    list[dict]
        Adjusted rows (new list; originals are not mutated).
    """
    adjusted: list[dict[str, Any]] = []
    for row in rows:
        values = [row.get(sid, 0.0) for sid in sample_ids]
        mean_val = sum(values) / len(values) if values else 0.0
        new_row: dict[str, Any] = {"tax_id": row["tax_id"]}
        for sid in sample_ids:
            new_row[sid] = row.get(sid, 0.0) - mean_val
        adjusted.append(new_row)
    return adjusted


# ---------------------------------------------------------------------------
# Outlier detection
# ---------------------------------------------------------------------------

# Scaling constant so that MAD ≈ standard deviation for normal data.
_MAD_SCALE = 1.4826

# Valid individual methods and the special "all" keyword.
_VALID_METHODS = frozenset({"mad", "iqr", "gmm", "all"})


def _run_mad(
    sample_ids: list[str],
    rows: list[dict[str, Any]],
    tax_names: dict[int, str],
    *,
    mad_threshold: float,
) -> tuple[list[FlaggedSample], int]:
    """Apply MAD-based outlier detection to every taxon row.

    Only samples *above* the population median are flagged: contamination
    is an additive signal, so a value far below the median is not
    informative for cross-species contamination detection.  When the MAD
    is exactly zero (more than half of the values are identical) the
    modified Z-score is undefined and the taxon is skipped; the count of
    skipped taxa is returned alongside the flag list so callers can
    surface it in QC output.
    """
    flagged: list[FlaggedSample] = []
    zero_mad_taxa = 0
    for row in rows:
        tid = row["tax_id"]
        values = [row.get(sid, 0.0) for sid in sample_ids]
        med = _median(values)
        mad_val = _mad(values)
        if mad_val == 0.0:
            zero_mad_taxa += 1
            continue
        scaled_mad = mad_val * _MAD_SCALE
        for sid, v in zip(sample_ids, values):
            modified_z = abs(v - med) / scaled_mad
            if modified_z > mad_threshold and v > med:
                flagged.append(
                    FlaggedSample(
                        sample_id=sid,
                        tax_id=tid,
                        taxon_name=tax_names.get(tid, ""),
                        value=v,
                        population_median=med,
                        deviation=modified_z,
                        method="mad",
                    )
                )
    return flagged, zero_mad_taxa


def _run_iqr(
    sample_ids: list[str],
    rows: list[dict[str, Any]],
    tax_names: dict[int, str],
    *,
    iqr_multiplier: float,
) -> tuple[list[FlaggedSample], int]:
    """Apply IQR-fence outlier detection to every taxon row.

    Flags samples whose value exceeds Tukey's upper fence
    ``Q3 + k·IQR`` and is strictly positive.  The lower fence
    ``Q1 - k·IQR`` is intentionally **not** evaluated because
    contamination produces excess (upper-tail) abundance; samples
    falling below the lower fence indicate undersampling and are not
    relevant for contamination detection.  Taxa with ``IQR == 0``
    (more than half of the values share both the median and one of the
    quartiles) are skipped and counted; the count is returned to the
    caller.
    """
    flagged: list[FlaggedSample] = []
    zero_iqr_taxa = 0
    for row in rows:
        tid = row["tax_id"]
        values = [row.get(sid, 0.0) for sid in sample_ids]
        q1, med, q3 = _quartiles(values)
        iqr = q3 - q1
        if iqr == 0.0:
            zero_iqr_taxa += 1
            continue
        upper_fence = q3 + iqr_multiplier * iqr
        for sid, v in zip(sample_ids, values):
            if v > upper_fence and v > 0:
                flagged.append(
                    FlaggedSample(
                        sample_id=sid,
                        tax_id=tid,
                        taxon_name=tax_names.get(tid, ""),
                        value=v,
                        population_median=med,
                        deviation=v - upper_fence,
                        method="iqr",
                    )
                )
    return flagged, zero_iqr_taxa


def _run_gmm(
    sample_ids: list[str],
    rows: list[dict[str, Any]],
    tax_names: dict[int, str],
    *,
    gmm_threshold: float,
) -> list[FlaggedSample]:
    """Apply two-component GMM outlier detection to every taxon row.

    Three guards prevent false positives:

    1. **BIC model selection** – a 2-component model is only accepted
       when BIC prefers it over a single Gaussian.
    2. **Minimum effective component size** – the contamination component
       must have ≥ 2 effective members (``sum(posteriors) ≥ 2``).
    3. **Separation criterion** – the component means must be at least
       3× the larger component sigma apart (Cohen's d ≥ 3).
    """
    flagged: list[FlaggedSample] = []
    for row in rows:
        tid = row["tax_id"]
        values = [row.get(sid, 0.0) for sid in sample_ids]
        med = _median(values)
        gmm_result = _fit_two_component_gmm(values)
        if gmm_result is None:
            continue
        if not _bic_favors_two_components(values, gmm_result):
            continue
        eff_contam = sum(gmm_result.posteriors)
        if eff_contam < 2.0:
            continue
        separation = gmm_result.mu_contam - gmm_result.mu_bg
        if separation < _GMM_MIN_COHENS_D * max(gmm_result.sigma_bg,
                                                gmm_result.sigma_contam):
            continue
        for sid, v, posterior in zip(
            sample_ids, values, gmm_result.posteriors
        ):
            if posterior > gmm_threshold:
                flagged.append(
                    FlaggedSample(
                        sample_id=sid,
                        tax_id=tid,
                        taxon_name=tax_names.get(tid, ""),
                        value=v,
                        population_median=med,
                        deviation=posterior,
                        method="gmm",
                    )
                )
    return flagged


def detect_outliers(
    matrix_path: str | Path,
    *,
    method: str = "all",
    mad_threshold: float = 3.5,
    iqr_multiplier: float = 1.5,
    gmm_threshold: float = 0.5,
    kitome_taxa: list[int] | None = None,
    subtract_background: bool = True,
    matrix_type: str | None = None,
) -> DetectionResult:
    """Run outlier detection on a sample-by-taxon matrix.

    Parameters
    ----------
    matrix_path:
        Path to the TSV matrix produced by ``csc-aggregate``.
    method:
        Detection method.  One of:

        * ``"all"`` *(default)* — run MAD, IQR, and GMM together and
          merge their flags.  Each sample–taxon pair is annotated with the
          method that produced the flag.  This is the recommended setting
          because the three methods are complementary: MAD and IQR are
          robust when the majority of samples are clean, while the GMM
          does not make that assumption and can detect contamination even
          when batch effects are widespread.
        * ``"mad"`` — MAD-based detection only.
        * ``"iqr"`` — IQR-based detection only.
        * ``"gmm"`` — Gaussian-mixture-model detection only.

    mad_threshold:
        Modified Z-score threshold above the population median for MAD
        detection (default: 3.5).  The default of 3.5 corresponds to
        approximately 3.5 standard deviations for normally distributed
        data, which provides a good balance between sensitivity and
        specificity for metagenomic contamination detection (Iglewicz
        and Hoaglin, 1993).
    iqr_multiplier:
        Multiplier *k* for the IQR upper fence ``Q3 + k·IQR``
        (default: 1.5).  The classic Tukey fence of 1.5 flags "mild"
        outliers; increase to 3.0 for "extreme" outliers only.
    gmm_threshold:
        Posterior-probability threshold for the contamination component
        in GMM detection (default: 0.5).  A sample is flagged when its
        posterior probability of belonging to the high-abundance
        component exceeds this value.  The default of 0.5 corresponds
        to the natural decision boundary of a two-component mixture
        (the sample is more likely to be contaminated than clean).
    kitome_taxa:
        NCBI taxonomy IDs to exclude before analysis.
    subtract_background:
        If ``True``, subtract the population mean per taxon before
        testing.
    matrix_type:
        Optional label identifying which type of input matrix this is
        (``"raw"``, ``"cpm"``, ``"abs"``, or any caller-defined value).
        Stored verbatim in the ``summary`` dict so downstream tooling
        (e.g. the HTML report) can disambiguate flag sets that were
        produced from different denominators.  When ``None`` the
        ``matrix_type`` key is omitted from the summary.

    Returns
    -------
    DetectionResult
        A dict with ``flagged`` (list of :class:`FlaggedSample`) and
        ``summary`` (aggregate statistics).

    Raises
    ------
    ValueError
        If *method* is not one of ``"all"``, ``"mad"``, ``"iqr"``,
        or ``"gmm"``.
    """
    method = method.lower()
    if method not in _VALID_METHODS:
        raise ValueError(f"Unknown detection method: {method!r}")

    sample_ids, rows, tax_names = load_matrix(matrix_path)
    rows = filter_kitome(rows, kitome_taxa)

    if subtract_background:
        rows = subtract_population_background(rows, sample_ids)

    flagged: list[FlaggedSample] = []
    zero_mad_taxa = 0
    zero_iqr_taxa = 0

    methods_to_run: list[str] = (
        ["mad", "iqr", "gmm"] if method == "all" else [method]
    )

    for m in methods_to_run:
        if m == "mad":
            mad_flags, zero_mad_taxa = _run_mad(
                sample_ids, rows, tax_names,
                mad_threshold=mad_threshold,
            )
            flagged.extend(mad_flags)
        elif m == "iqr":
            iqr_flags, zero_iqr_taxa = _run_iqr(
                sample_ids, rows, tax_names,
                iqr_multiplier=iqr_multiplier,
            )
            flagged.extend(iqr_flags)
        else:  # gmm
            flagged.extend(
                _run_gmm(sample_ids, rows, tax_names,
                         gmm_threshold=gmm_threshold)
            )

    summary: dict[str, Any] = {
        "total_samples": len(sample_ids),
        "total_taxa_analysed": len(rows),
        "method": method,
        "flagged_count": len(flagged),
        "flagged_samples": sorted({f["sample_id"] for f in flagged}),
        "kitome_taxa_excluded": len(kitome_taxa) if kitome_taxa else 0,
        "background_subtracted": subtract_background,
    }
    if matrix_type is not None:
        summary["matrix_type"] = matrix_type
    # Diagnostic counters: surface taxa skipped by MAD/IQR due to zero
    # dispersion so users can interpret apparent low flag counts.  Only
    # included for methods that were actually run.
    if "mad" in methods_to_run:
        summary["taxa_skipped_mad_zero"] = zero_mad_taxa
    if "iqr" in methods_to_run:
        summary["taxa_skipped_iqr_zero"] = zero_iqr_taxa

    logger.info(
        "Detection complete: %d flags across %d samples",
        len(flagged),
        len(summary["flagged_samples"]),
    )

    return DetectionResult(flagged=flagged, summary=summary)
