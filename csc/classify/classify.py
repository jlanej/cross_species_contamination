"""Core classification logic using Kraken2.

This module provides a Python API for running Kraken2 taxonomic classification
on FASTQ files extracted by the ``csc.extract`` module.  It validates the
database, builds the Kraken2 command, runs classification, and returns
structured results.

AI assistance acknowledgment: This module was developed with AI assistance.
Best practices in the bioinformatics field should always take precedence over
specific implementation details.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any, TypedDict

logger = logging.getLogger(__name__)

#: Files that must be present in a valid Kraken2 database directory.
REQUIRED_DB_FILES = ["hash.k2d", "opts.k2d", "taxo.k2d"]


class ClassificationResult(TypedDict):
    """Return type for :func:`classify_reads`."""

    report: Path
    output: Path
    sample_id: str
    input_files: list[str]


def _find_kraken2() -> str:
    """Return the path to the kraken2 executable or raise."""
    path = shutil.which("kraken2")
    if path is None:
        raise FileNotFoundError(
            "kraken2 not found on PATH. Install Kraken2 to use this module."
        )
    return path


def _get_kraken2_version(kraken2: str) -> str:
    """Return the Kraken2 version string."""
    try:
        proc = subprocess.run(
            [kraken2, "--version"],
            capture_output=True,
            text=True,
            check=False,
        )
        # First line typically: "Kraken version 2.x.x"
        for line in proc.stdout.splitlines():
            if "version" in line.lower():
                return line.strip()
        return proc.stdout.splitlines()[0].strip() if proc.stdout else "unknown"
    except Exception:
        return "unknown"


def validate_database(db_path: str | Path) -> Path:
    """Validate that *db_path* is a usable Kraken2 database directory.

    Parameters
    ----------
    db_path:
        Path to the Kraken2 database directory.

    Returns
    -------
    Path
        Resolved absolute path to the database.

    Raises
    ------
    FileNotFoundError
        If *db_path* does not exist.
    ValueError
        If required database files are missing.
    """
    db = Path(db_path).resolve()
    if not db.is_dir():
        raise FileNotFoundError(f"Kraken2 database directory not found: {db}")

    missing = [f for f in REQUIRED_DB_FILES if not (db / f).exists()]
    if missing:
        raise ValueError(
            f"Kraken2 database at {db} is missing required files: "
            f"{', '.join(missing)}. "
            f"Expected files: {', '.join(REQUIRED_DB_FILES)}"
        )
    return db


def build_classify_command(
    input_files: list[Path],
    *,
    db: Path,
    output: Path,
    report: Path,
    confidence: float = 0.0,
    threads: int = 1,
    memory_mapping: bool = False,
    paired: bool = False,
) -> list[str]:
    """Build the Kraken2 command for taxonomic classification.

    Parameters
    ----------
    input_files:
        One or two FASTQ input files (single-end or paired-end).
    db:
        Path to the Kraken2 database directory.
    output:
        Path for the per-read classification output file.
    report:
        Path for the Kraken2 summary report file.
    confidence:
        Minimum confidence score threshold (0.0–1.0).
    threads:
        Number of threads for Kraken2 to use.
    memory_mapping:
        If ``True``, use memory mapping instead of loading the full DB into
        RAM (``--memory-mapping``).  Reduces memory at the cost of speed.
    paired:
        If ``True``, treat the two input files as paired-end reads.

    Returns
    -------
    list[str]
        The command list suitable for :func:`subprocess.run`.
    """
    kraken2 = _find_kraken2()

    cmd = [
        kraken2,
        "--db", str(db),
        "--output", str(output),
        "--report", str(report),
        "--confidence", str(confidence),
        "--threads", str(threads),
    ]

    if memory_mapping:
        cmd.append("--memory-mapping")

    if paired and len(input_files) == 2:
        cmd.append("--paired")

    cmd.extend(str(f) for f in input_files)
    return cmd


def classify_reads(
    input_files: list[str | Path],
    output_dir: str | Path,
    *,
    db: str | Path,
    sample_id: str | None = None,
    confidence: float = 0.0,
    threads: int = 1,
    memory_mapping: bool = False,
    paired: bool = False,
) -> ClassificationResult:
    """Classify FASTQ reads using Kraken2.

    Parameters
    ----------
    input_files:
        List of FASTQ file paths to classify.  For paired-end data, supply
        exactly two files (R1, R2) and set ``paired=True``.
    output_dir:
        Directory for Kraken2 output files.
    db:
        Path to the Kraken2 database directory.
    sample_id:
        Basename for output files.  Defaults to the first input file stem.
    confidence:
        Minimum confidence score threshold (0.0–1.0).
    threads:
        Number of threads for Kraken2.
    memory_mapping:
        If ``True``, use ``--memory-mapping`` to reduce RAM usage.
    paired:
        If ``True``, treat input files as paired-end reads.

    Returns
    -------
    ClassificationResult
        A dict with ``report``, ``output``, ``sample_id``, and
        ``input_files`` keys.

    Raises
    ------
    FileNotFoundError
        If kraken2 is not installed, the database is invalid, or input
        files do not exist.
    ValueError
        If input validation fails.
    RuntimeError
        If the kraken2 command exits with a non-zero status.
    """
    # Resolve paths
    resolved_inputs = [Path(f).resolve() for f in input_files]
    output_dir = Path(output_dir).resolve()
    db_path = validate_database(db)

    # Validate input files exist
    for fp in resolved_inputs:
        if not fp.exists():
            raise FileNotFoundError(f"Input file not found: {fp}")

    if not resolved_inputs:
        raise ValueError("At least one input file is required.")

    if paired and len(resolved_inputs) != 2:
        raise ValueError(
            f"Paired-end mode requires exactly 2 input files, "
            f"got {len(resolved_inputs)}."
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    # Derive sample ID
    if sample_id is None:
        stem = resolved_inputs[0].name
        # Strip common FASTQ extensions
        for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if stem.endswith(ext):
                stem = stem[: -len(ext)]
                break
        sample_id = stem

    report_path = output_dir / f"{sample_id}.kraken2.report.txt"
    output_path = output_dir / f"{sample_id}.kraken2.output.txt"

    # Log classification info
    kraken2 = _find_kraken2()
    version = _get_kraken2_version(kraken2)
    logger.info("Kraken2 version: %s", version)
    logger.info("Database: %s", db_path)
    logger.info(
        "Classifying %d file(s) for sample '%s'",
        len(resolved_inputs),
        sample_id,
    )

    cmd = build_classify_command(
        resolved_inputs,
        db=db_path,
        output=output_path,
        report=report_path,
        confidence=confidence,
        threads=threads,
        memory_mapping=memory_mapping,
        paired=paired,
    )
    logger.debug("Running: %s", " ".join(cmd))

    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Kraken2 failed (exit {proc.returncode}): {proc.stderr}"
        )

    if proc.stderr:
        # Kraken2 prints stats to stderr; log them at INFO level
        for line in proc.stderr.strip().splitlines():
            logger.info("kraken2: %s", line)

    logger.info("Classification complete for sample '%s'", sample_id)

    return {
        "report": report_path,
        "output": output_path,
        "sample_id": sample_id,
        "input_files": [str(f) for f in resolved_inputs],
    }
