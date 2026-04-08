"""Report generation for the CSC detect module.

Writes three output files:

* **flagged_samples.tsv** – one row per (sample, taxon) outlier flag.
* **qc_summary.json** – aggregate QC statistics.
* **quarantine_list.txt** – plain list of sample IDs recommended for
  quarantine or further investigation.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

from csc.detect.detect import DetectionResult, FlaggedSample

logger = logging.getLogger(__name__)


def write_flagged_samples(
    flagged: list[FlaggedSample],
    path: Path,
) -> Path:
    """Write the flagged-sample table as TSV.

    Columns: sample_id, tax_id, taxon_name, value, population_median,
    deviation, method.
    """
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "sample_id",
            "tax_id",
            "taxon_name",
            "value",
            "population_median",
            "deviation",
            "method",
        ])
        for f in flagged:
            writer.writerow([
                f["sample_id"],
                f["tax_id"],
                f["taxon_name"],
                f"{f['value']:.4f}",
                f"{f['population_median']:.4f}",
                f"{f['deviation']:.4f}",
                f["method"],
            ])
    logger.info("Wrote flagged samples to %s", path)
    return path


def write_qc_summary(
    summary: dict[str, Any],
    path: Path,
) -> Path:
    """Write the QC summary as a JSON file."""
    with open(path, "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    logger.info("Wrote QC summary to %s", path)
    return path


def write_quarantine_list(
    flagged_samples: list[str],
    path: Path,
) -> Path:
    """Write a plain-text quarantine / alert list (one sample per line)."""
    with open(path, "w") as fh:
        for sid in sorted(flagged_samples):
            fh.write(sid + "\n")
    logger.info("Wrote quarantine list (%d samples) to %s", len(flagged_samples), path)
    return path


def generate_report(
    result: DetectionResult,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Generate all report files from a :class:`DetectionResult`.

    Parameters
    ----------
    result:
        Output of :func:`csc.detect.detect.detect_outliers`.
    output_dir:
        Directory for output files.  Created if it does not exist.

    Returns
    -------
    dict[str, Path]
        Mapping of report name to file path:
        ``"flagged_samples"``, ``"qc_summary"``, ``"quarantine_list"``.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    flagged_path = write_flagged_samples(
        result["flagged"],
        output_dir / "flagged_samples.tsv",
    )
    qc_path = write_qc_summary(
        result["summary"],
        output_dir / "qc_summary.json",
    )
    quarantine_path = write_quarantine_list(
        result["summary"]["flagged_samples"],
        output_dir / "quarantine_list.txt",
    )

    return {
        "flagged_samples": flagged_path,
        "qc_summary": qc_path,
        "quarantine_list": quarantine_path,
    }
