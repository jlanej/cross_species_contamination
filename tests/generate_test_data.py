#!/usr/bin/env python3
"""Generate synthetic BAM test data for the extraction pipeline.

Creates small BAM files containing a mixture of mapped, unmapped, and
low-MAPQ reads so the extraction module can be validated in CI without
any external data dependencies.

This script requires ``pysam`` (``pip install pysam``).
"""

from __future__ import annotations

import os
import random
import tempfile
from pathlib import Path

import pysam

# Fixed seed for reproducibility
random.seed(42)

# A minimal reference: two short "chromosomes"
REFERENCE_SEQS = {
    "chr1": "A" * 500,
    "chr2": "C" * 500,
}

CONTAMINANT_SEQ = "GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT"  # 50 bp


def _random_seq(length: int = 100) -> str:
    return "".join(random.choices("ACGT", k=length))


def _qual_string(length: int = 100) -> str:
    return "I" * length  # Phred 40 across the board


def generate_reference_fasta(outdir: Path) -> Path:
    """Write a minimal FASTA reference and its .fai index."""
    ref_path = outdir / "reference.fa"
    with open(ref_path, "w") as fh:
        for name, seq in REFERENCE_SEQS.items():
            fh.write(f">{name}\n{seq}\n")
    pysam.faidx(str(ref_path))
    return ref_path


def generate_test_bam(
    outdir: Path,
    ref_path: Path,
    *,
    n_mapped: int = 20,
    n_unmapped: int = 10,
    n_low_mapq: int = 5,
    low_mapq_value: int = 3,
    paired: bool = True,
    filename: str = "test_sample.bam",
) -> Path:
    """Create a synthetic BAM file with known read composition.

    Parameters
    ----------
    outdir : Path
        Directory to write to.
    ref_path : Path
        Reference FASTA (must have .fai index).
    n_mapped : int
        Number of well-mapped reads/pairs.
    n_unmapped : int
        Number of fully unmapped reads/pairs.
    n_low_mapq : int
        Number of mapped reads with MAPQ below *low_mapq_value*.
    low_mapq_value : int
        MAPQ value to assign to low-quality mapped reads.
    paired : bool
        If True, generate paired-end reads.
    filename : str
        Output BAM file name.

    Returns
    -------
    Path
        Path to the sorted, indexed BAM file.
    """
    bam_path = outdir / filename
    header = pysam.AlignmentHeader.from_references(
        list(REFERENCE_SEQS.keys()),
        [len(s) for s in REFERENCE_SEQS.values()],
    )

    reads: list[pysam.AlignedSegment] = []
    read_idx = 0

    # --- Well-mapped reads ---
    for i in range(n_mapped):
        rname = f"mapped_{read_idx}"
        chrom_idx = i % len(REFERENCE_SEQS)
        pos = random.randint(0, 350)
        seq = _random_seq(100)

        if paired:
            # Read 1
            a = pysam.AlignedSegment(header)
            a.query_name = rname
            a.query_sequence = seq
            a.flag = 0x1 | 0x40  # paired, first in pair
            a.reference_id = chrom_idx
            a.reference_start = pos
            a.mapping_quality = 60
            a.cigar = [(0, 100)]  # 100M
            a.query_qualities = pysam.qualitystring_to_array(_qual_string(100))
            a.mate_is_reverse = True
            a.next_reference_id = chrom_idx
            a.next_reference_start = pos + 150
            a.template_length = 250
            reads.append(a)

            # Read 2
            b = pysam.AlignedSegment(header)
            b.query_name = rname
            b.query_sequence = _random_seq(100)
            b.flag = 0x1 | 0x80 | 0x10  # paired, second, reverse
            b.reference_id = chrom_idx
            b.reference_start = pos + 150
            b.mapping_quality = 60
            b.cigar = [(0, 100)]
            b.query_qualities = pysam.qualitystring_to_array(_qual_string(100))
            b.mate_is_reverse = False
            b.next_reference_id = chrom_idx
            b.next_reference_start = pos
            b.template_length = -250
            reads.append(b)
        else:
            a = pysam.AlignedSegment(header)
            a.query_name = rname
            a.query_sequence = seq
            a.flag = 0
            a.reference_id = chrom_idx
            a.reference_start = pos
            a.mapping_quality = 60
            a.cigar = [(0, 100)]
            a.query_qualities = pysam.qualitystring_to_array(_qual_string(100))
            reads.append(a)
        read_idx += 1

    # --- Unmapped reads (simulated contaminant) ---
    for i in range(n_unmapped):
        rname = f"unmapped_{read_idx}"
        seq = CONTAMINANT_SEQ + _random_seq(50)  # 100 bp contaminant-like

        if paired:
            # Both mates unmapped
            a = pysam.AlignedSegment(header)
            a.query_name = rname
            a.query_sequence = seq
            a.flag = 0x1 | 0x4 | 0x8 | 0x40  # paired, unmapped, mate unmapped, R1
            a.reference_id = -1
            a.reference_start = -1
            a.mapping_quality = 0
            a.query_qualities = pysam.qualitystring_to_array(_qual_string(100))
            reads.append(a)

            b = pysam.AlignedSegment(header)
            b.query_name = rname
            b.query_sequence = _random_seq(100)
            b.flag = 0x1 | 0x4 | 0x8 | 0x80  # paired, unmapped, mate unmapped, R2
            b.reference_id = -1
            b.reference_start = -1
            b.mapping_quality = 0
            b.query_qualities = pysam.qualitystring_to_array(_qual_string(100))
            reads.append(b)
        else:
            a = pysam.AlignedSegment(header)
            a.query_name = rname
            a.query_sequence = seq
            a.flag = 0x4  # unmapped
            a.reference_id = -1
            a.reference_start = -1
            a.mapping_quality = 0
            a.query_qualities = pysam.qualitystring_to_array(_qual_string(100))
            reads.append(a)
        read_idx += 1

    # --- Low-MAPQ reads ---
    for i in range(n_low_mapq):
        rname = f"lowmapq_{read_idx}"
        chrom_idx = i % len(REFERENCE_SEQS)
        pos = random.randint(0, 350)
        seq = _random_seq(100)

        if paired:
            a = pysam.AlignedSegment(header)
            a.query_name = rname
            a.query_sequence = seq
            a.flag = 0x1 | 0x40
            a.reference_id = chrom_idx
            a.reference_start = pos
            a.mapping_quality = low_mapq_value
            a.cigar = [(0, 100)]
            a.query_qualities = pysam.qualitystring_to_array(_qual_string(100))
            a.next_reference_id = chrom_idx
            a.next_reference_start = pos + 150
            a.template_length = 250
            reads.append(a)

            b = pysam.AlignedSegment(header)
            b.query_name = rname
            b.query_sequence = _random_seq(100)
            b.flag = 0x1 | 0x80 | 0x10
            b.reference_id = chrom_idx
            b.reference_start = pos + 150
            b.mapping_quality = low_mapq_value
            b.cigar = [(0, 100)]
            b.query_qualities = pysam.qualitystring_to_array(_qual_string(100))
            b.next_reference_id = chrom_idx
            b.next_reference_start = pos
            b.template_length = -250
            reads.append(b)
        else:
            a = pysam.AlignedSegment(header)
            a.query_name = rname
            a.query_sequence = seq
            a.flag = 0
            a.reference_id = chrom_idx
            a.reference_start = pos
            a.mapping_quality = low_mapq_value
            a.cigar = [(0, 100)]
            a.query_qualities = pysam.qualitystring_to_array(_qual_string(100))
            reads.append(a)
        read_idx += 1

    # Sort and write
    reads.sort(key=lambda r: (r.reference_id, r.reference_start))
    tmp_unsorted = outdir / f"_unsorted_{filename}"
    with pysam.AlignmentFile(str(tmp_unsorted), "wb", header=header) as outf:
        for r in reads:
            outf.write(r)

    pysam.sort("-o", str(bam_path), str(tmp_unsorted))
    pysam.index(str(bam_path))
    tmp_unsorted.unlink()

    return bam_path


def generate_test_data(outdir: str | Path | None = None) -> dict[str, Path]:
    """Generate a complete set of test data and return a dict of paths.

    If *outdir* is ``None``, a temporary directory is used.
    """
    if outdir is None:
        outdir = Path(tempfile.mkdtemp(prefix="extract_unmapped_test_"))
    else:
        outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ref = generate_reference_fasta(outdir)
    bam = generate_test_bam(outdir, ref)

    return {
        "reference": ref,
        "bam": bam,
        "outdir": outdir,
    }


if __name__ == "__main__":
    import sys

    dest = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    paths = generate_test_data(dest)
    for k, v in paths.items():
        print(f"{k}: {v}")
