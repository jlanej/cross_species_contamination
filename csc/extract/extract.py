"""Core extraction logic for unmapped/poorly-mapped reads.

This module provides streaming extraction of unmapped reads (SAM flag 4) and
optionally poorly mapped reads (below a MAPQ threshold) from BAM/CRAM files
using samtools.  No intermediate disk files are created; all data is streamed
through pipes.

AI assistance acknowledgment: This module was developed with AI assistance.
Best practices in the bioinformatics field should always take precedence over
specific implementation details.
"""

from __future__ import annotations

import datetime as _dt
import json
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, TypedDict

logger = logging.getLogger(__name__)

SUPPORTED_EXTENSIONS = {".bam", ".cram"}

# Schema version for reads_summary.json.  Bump when breaking changes are
# introduced so downstream consumers (aggregate, reporting, PR #55
# summary report) can detect and adapt.
READS_SUMMARY_SCHEMA_VERSION = "1.0"


class ReadsSummary(TypedDict):
    """Per-sample mapped/unmapped read counts derived from ``samtools idxstats``.

    Written to ``{sample_id}.reads_summary.json`` alongside the raw
    idxstats TSV.  Downstream modules (e.g. :mod:`csc.aggregate`) use
    ``total_reads`` as the denominator for absolute contamination
    burden.
    """

    schema_version: str
    sample_id: str
    input: str
    extraction_time: str
    total_mapped: int
    total_unmapped: int
    total_reads: int
    per_chromosome: list[dict[str, Any]]


class ExtractionResult(TypedDict, total=False):
    """Return type for :func:`extract_reads`."""

    files: dict[str, Path]
    read_count: int
    sample_id: str
    input: str
    # idxstats-derived fields (always present when idxstats succeeds)
    total_mapped: int
    total_unmapped: int
    total_reads: int
    idxstats_path: Path
    reads_summary_path: Path


def run_idxstats(
    input_path: str | Path,
    output_dir: str | Path,
    *,
    sample_id: str | None = None,
    reference: str | Path | None = None,
    threads: int = 1,
) -> ReadsSummary:
    """Run ``samtools idxstats`` on an indexed BAM/CRAM and write sidecars.

    Two files are written into *output_dir*:

    * ``{sample_id}.idxstats.tsv`` – raw TSV from ``samtools idxstats``
      with columns ``chrom<TAB>length<TAB>mapped<TAB>unmapped``.  The
      final row (``*``) reports reads unmapped and without coordinates
      (i.e. reads whose mate is also unmapped).
    * ``{sample_id}.reads_summary.json`` – structured summary with
      ``total_mapped``, ``total_unmapped``, ``total_reads`` (sum),
      schema version, and a per-chromosome breakdown for transparency
      and reproducibility.

    ``samtools idxstats`` reads only the BAM index and is effectively
    free (milliseconds per sample) compared to a full scan.

    Notes
    -----
    ``total_unmapped`` returned here is the count of reads in the BAM
    file with the unmapped flag set.  Downstream classifier totals may
    differ because post-extraction filtering (e.g. adapter trimming,
    host-read removal, or the ``--mapq`` filter in :func:`extract_reads`)
    can drop or add reads before classification.  Always use
    ``reads_summary.json`` (emitted at extract time) as the authoritative
    denominator for absolute contamination burden.

    Parameters
    ----------
    input_path:
        BAM or CRAM file.  Must be indexed (``.bai`` or ``.crai``).
    output_dir:
        Directory to write sidecar files into.  Created if needed.
    sample_id:
        Output file basename.  Defaults to the input file stem.
    reference:
        Reference FASTA for CRAM (optional for BAM).
    threads:
        Additional samtools decompression threads.

    Returns
    -------
    ReadsSummary
        Parsed counts plus provenance (sample_id, input, extraction_time).
    """
    input_path = Path(input_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if sample_id is None:
        sample_id = input_path.stem
        if sample_id.endswith(".sorted"):
            sample_id = sample_id[: -len(".sorted")]

    samtools = _find_samtools()
    # Note: ``samtools idxstats`` reads only the BAM/CRAM index; it does
    # not decode records and therefore does not accept --reference.  The
    # *reference* argument is kept in the API for signature parity with
    # other extract helpers but is intentionally unused here.
    del reference  # unused on purpose; index-only operation

    cmd = [samtools, "idxstats", "-@", str(threads), str(input_path)]
    logger.debug("Running: %s", " ".join(cmd))
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"samtools idxstats failed (exit {proc.returncode}): "
            f"{proc.stderr.decode(errors='replace')}"
        )
    raw = proc.stdout.decode(errors="replace")

    # Persist raw TSV verbatim.
    idxstats_path = output_dir / f"{sample_id}.idxstats.tsv"
    idxstats_path.write_text(raw)

    # Parse into structured per-chromosome rows.
    per_chrom: list[dict[str, Any]] = []
    total_mapped = 0
    total_unmapped = 0
    for line in raw.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            logger.warning("Skipping malformed idxstats line: %r", line)
            continue
        chrom, length_s, mapped_s, unmapped_s = parts[0], parts[1], parts[2], parts[3]
        try:
            length = int(length_s)
            mapped = int(mapped_s)
            unmapped = int(unmapped_s)
        except ValueError:
            logger.warning("Skipping non-numeric idxstats line: %r", line)
            continue
        per_chrom.append(
            {
                "chrom": chrom,
                "length": length,
                "mapped": mapped,
                "unmapped": unmapped,
            }
        )
        total_mapped += mapped
        total_unmapped += unmapped

    summary: ReadsSummary = {
        "schema_version": READS_SUMMARY_SCHEMA_VERSION,
        "sample_id": sample_id,
        "input": str(input_path),
        "extraction_time": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "total_mapped": total_mapped,
        "total_unmapped": total_unmapped,
        "total_reads": total_mapped + total_unmapped,
        "per_chromosome": per_chrom,
    }

    reads_summary_path = output_dir / f"{sample_id}.reads_summary.json"
    with open(reads_summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)

    logger.info(
        "idxstats for %s: mapped=%d unmapped=%d total=%d",
        sample_id,
        total_mapped,
        total_unmapped,
        total_mapped + total_unmapped,
    )
    return summary


def _find_samtools() -> str:
    """Return the path to the samtools executable or raise."""
    path = shutil.which("samtools")
    if path is None:
        raise FileNotFoundError(
            "samtools not found on PATH. Install samtools >= 1.12 to use this tool."
        )
    return path


def _validate_input(input_path: Path, reference: Path | None) -> None:
    """Validate that the input file exists and has a supported extension."""
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    ext = input_path.suffix.lower()
    if ext not in SUPPORTED_EXTENSIONS:
        raise ValueError(
            f"Unsupported file format '{ext}'. Supported: {SUPPORTED_EXTENSIONS}"
        )
    if ext == ".cram" and reference is not None and not reference.exists():
        raise FileNotFoundError(f"Reference file not found: {reference}")


def _resolve_reference(
    reference: Path | None,
    input_path: Path,
) -> Path | None:
    """Determine the reference FASTA to use for a given input file.

    Resolution order:

    1. Explicit *reference* argument (CLI ``--reference``).
    2. ``REF_PATH`` environment variable.
    3. ``None`` – acceptable for BAM, but will raise for CRAM.
    """
    if reference is not None:
        return reference

    env_ref = os.environ.get("REF_PATH")
    if env_ref:
        ref = Path(env_ref)
        if ref.is_file():
            logger.debug("Using reference from REF_PATH env: %s", ref)
            return ref
        logger.warning("REF_PATH is set but '%s' is not a file; ignoring.", env_ref)

    if input_path.suffix.lower() == ".cram":
        raise FileNotFoundError(
            f"CRAM input requires a reference FASTA. Supply --reference or "
            f"set the REF_PATH environment variable. Input: {input_path}"
        )
    return None


def build_extract_command(
    input_path: Path,
    *,
    output_r1: Path | None = None,
    output_r2: Path | None = None,
    output_single: Path | None = None,
    output_other: Path | None = None,
    output_interleaved: Path | None = None,
    mapq_threshold: int | None = None,
    threads: int = 1,
    reference: Path | None = None,
) -> list[list[str]]:
    """Build the samtools command(s) for extracting unmapped reads.

    Returns a list of command lists.  When *mapq_threshold* is ``None`` a
    single ``samtools fastq -f 4`` command is returned.  When a threshold is
    given, a ``samtools view`` piped into ``samtools fastq`` is returned so
    that both unmapped and poorly-mapped reads are captured.

    Parameters
    ----------
    input_path:
        Path to the input BAM or CRAM file.
    output_r1, output_r2:
        Output paths for paired-end read 1 / read 2 FASTQ files.
    output_single:
        Output path for singleton reads.
    output_other:
        Output path for reads that are neither R1, R2, nor singleton.
    output_interleaved:
        Output path for interleaved FASTQ (all reads in one file).
    mapq_threshold:
        If set, also extract mapped reads with MAPQ < this value.
    threads:
        Number of additional samtools decompression threads.
    reference:
        Reference FASTA required for CRAM input (optional for BAM).
    """
    samtools = _find_samtools()

    # Common reference arg
    ref_args: list[str] = []
    if reference is not None:
        ref_args = ["--reference", str(reference)]

    thread_args = ["-@", str(threads)]

    # Build output args for samtools fastq
    out_args: list[str] = []
    if output_interleaved is not None:
        out_args += ["-o", str(output_interleaved)]
    if output_r1 is not None:
        out_args += ["-1", str(output_r1)]
    if output_r2 is not None:
        out_args += ["-2", str(output_r2)]
    if output_single is not None:
        out_args += ["-s", str(output_single)]
    if output_other is not None:
        out_args += ["-0", str(output_other)]

    if mapq_threshold is None:
        # Simple mode: extract only unmapped reads
        cmd = [
            samtools,
            "fastq",
            "-f",
            "4",
            *thread_args,
            *ref_args,
            *out_args,
            str(input_path),
        ]
        return [cmd]
    else:
        # Advanced mode: unmapped OR MAPQ < threshold
        # Uses samtools expression filter (requires samtools >= 1.12)
        view_cmd = [
            samtools,
            "view",
            "-h",
            "-F",
            "256",  # exclude secondary alignments
            "-e",
            f"flag.unmap || mapq < {mapq_threshold}",
            *thread_args,
            *ref_args,
            str(input_path),
        ]
        fastq_cmd = [
            samtools,
            "fastq",
            *thread_args,
            *out_args,
            "-",
        ]
        return [view_cmd, fastq_cmd]


def extract_reads(
    input_path: str | Path,
    output_dir: str | Path,
    *,
    sample_id: str | None = None,
    mapq_threshold: int | None = None,
    threads: int = 1,
    reference: str | Path | None = None,
    interleaved: bool = False,
) -> ExtractionResult:
    """Extract unmapped (and optionally low-MAPQ) reads from a BAM/CRAM file.

    All output files are gzip-compressed FASTQ.  Returns a dict with keys:

    * ``"files"`` – dict mapping output type to :class:`Path`
    * ``"read_count"`` – total reads written
    * ``"sample_id"`` – resolved sample identifier
    * ``"input"`` – resolved input path

    Parameters
    ----------
    input_path:
        BAM or CRAM input file.
    output_dir:
        Directory for output FASTQ files.
    sample_id:
        Basename for output files.  Defaults to the input file stem.
    mapq_threshold:
        If set, also extract reads with MAPQ below this value.
    threads:
        Additional decompression threads for samtools.
    reference:
        Reference FASTA for CRAM files.
    interleaved:
        If True, write a single interleaved FASTQ instead of split files.
    """
    input_path = Path(input_path).resolve()
    output_dir = Path(output_dir).resolve()
    reference = Path(reference).resolve() if reference else None

    _validate_input(input_path, reference)
    reference = _resolve_reference(reference, input_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    if sample_id is None:
        # Strip .bam / .cram and any .sorted etc.
        sample_id = input_path.stem
        if sample_id.endswith(".sorted"):
            sample_id = sample_id[: -len(".sorted")]

    outputs: dict[str, Path] = {}
    if interleaved:
        outputs["interleaved"] = output_dir / f"{sample_id}.unmapped.fastq.gz"
    else:
        outputs["r1"] = output_dir / f"{sample_id}.unmapped.R1.fastq.gz"
        outputs["r2"] = output_dir / f"{sample_id}.unmapped.R2.fastq.gz"
        outputs["singleton"] = output_dir / f"{sample_id}.unmapped.singleton.fastq.gz"
        outputs["other"] = output_dir / f"{sample_id}.unmapped.other.fastq.gz"

    if interleaved:
        cmds = build_extract_command(
            input_path,
            output_interleaved=outputs["interleaved"],
            mapq_threshold=mapq_threshold,
            threads=threads,
            reference=reference,
        )
    else:
        cmds = build_extract_command(
            input_path,
            output_r1=outputs["r1"],
            output_r2=outputs["r2"],
            output_single=outputs["singleton"],
            output_other=outputs["other"],
            mapq_threshold=mapq_threshold,
            threads=threads,
            reference=reference,
        )

    logger.info("Extracting reads from %s", input_path)
    if mapq_threshold is not None:
        logger.info("Including reads with MAPQ < %d", mapq_threshold)

    if len(cmds) == 1:
        # Single command
        cmd = cmds[0]
        logger.debug("Running: %s", " ".join(cmd))
        proc = subprocess.run(cmd, stderr=subprocess.PIPE, check=False)
        if proc.returncode != 0:
            raise RuntimeError(
                f"samtools failed (exit {proc.returncode}): "
                f"{proc.stderr.decode(errors='replace')}"
            )
    else:
        # Piped commands (view | fastq)
        logger.debug(
            "Running pipeline: %s",
            " | ".join(" ".join(c) for c in cmds),
        )
        view_proc = subprocess.Popen(
            cmds[0],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        fastq_proc = subprocess.Popen(
            cmds[1],
            stdin=view_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # Allow view_proc to receive SIGPIPE if fastq_proc exits
        if view_proc.stdout is None:
            raise RuntimeError("samtools view process has no stdout handle")
        view_proc.stdout.close()
        _, fastq_err = fastq_proc.communicate()
        view_proc.wait()

        if fastq_proc.returncode != 0:
            raise RuntimeError(
                f"samtools fastq failed (exit {fastq_proc.returncode}): "
                f"{fastq_err.decode(errors='replace')}"
            )
        if view_proc.returncode != 0:
            view_err = view_proc.stderr.read() if view_proc.stderr else b""
            raise RuntimeError(
                f"samtools view failed (exit {view_proc.returncode}): "
                f"{view_err.decode(errors='replace')}"
            )

    # Filter out empty output files from the result
    final_outputs: dict[str, Path] = {}
    for key, path in outputs.items():
        if path.exists() and path.stat().st_size > 0:
            final_outputs[key] = path

    read_count = _count_reads(final_outputs)
    logger.info("Extracted %d reads to %s", read_count, output_dir)

    result: ExtractionResult = {
        "files": final_outputs,
        "read_count": read_count,
        "sample_id": sample_id,
        "input": str(input_path),
    }

    # Always emit idxstats sidecars.  Failures are logged but do not
    # abort extraction – idxstats requires an index which may be absent
    # for test fixtures or streamed inputs.
    try:
        summary = run_idxstats(
            input_path,
            output_dir,
            sample_id=sample_id,
            reference=reference,
            threads=threads,
        )
        result["total_mapped"] = summary["total_mapped"]
        result["total_unmapped"] = summary["total_unmapped"]
        result["total_reads"] = summary["total_reads"]
        result["idxstats_path"] = output_dir / f"{sample_id}.idxstats.tsv"
        result["reads_summary_path"] = output_dir / f"{sample_id}.reads_summary.json"
    except Exception as exc:  # noqa: BLE001 - non-fatal sidecar
        logger.warning("idxstats sidecar skipped for %s: %s", sample_id, exc)

    return result


def _count_reads(outputs: dict[str, Path]) -> int:
    """Quick line-count based read estimate from gzipped FASTQ."""
    import gzip

    total = 0
    for path in outputs.values():
        if not path.exists():
            continue
        try:
            with gzip.open(path, "rt") as fh:
                lines = sum(1 for _ in fh)
            total += lines // 4
        except Exception:
            pass
    return total
