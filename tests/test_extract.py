"""Tests for the unmapped-read extraction pipeline.

The synthetic test BAM (generated in conftest.py) contains:
  - 20 well-mapped paired-end reads  (MAPQ 60)
  - 10 unmapped paired-end reads     (flag 4, simulated contaminant)
  -  5 low-MAPQ paired-end reads     (MAPQ 3)

These counts are used to verify that extraction is correct.
"""

from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

import pytest

from extract_unmapped.extract import (
    build_extract_command,
    extract_reads,
    _find_samtools,
    _validate_input,
)

# Expected counts from generate_test_data defaults
EXPECTED_UNMAPPED_READS = 10  # pairs
EXPECTED_LOW_MAPQ_READS = 5  # pairs


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def count_fastq_reads(path: Path) -> int:
    """Count reads in a (possibly gzipped) FASTQ file."""
    if not path.exists():
        return 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return sum(1 for line in fh) // 4


# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------

class TestValidation:
    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            _validate_input(tmp_path / "nonexistent.bam", None)

    def test_unsupported_extension_raises(self, tmp_path: Path) -> None:
        bad = tmp_path / "reads.fastq"
        bad.touch()
        with pytest.raises(ValueError, match="Unsupported"):
            _validate_input(bad, None)

    def test_valid_bam_extension(self, test_bam: Path) -> None:
        _validate_input(test_bam, None)  # should not raise


class TestBuildCommand:
    def test_unmapped_only_single_command(self, test_bam: Path) -> None:
        cmds = build_extract_command(test_bam)
        assert len(cmds) == 1
        cmd = cmds[0]
        assert "fastq" in cmd
        assert "-f" in cmd
        assert "4" in cmd

    def test_mapq_filter_produces_pipeline(self, test_bam: Path) -> None:
        cmds = build_extract_command(test_bam, mapq_threshold=10)
        assert len(cmds) == 2
        # First command is samtools view with expression filter
        assert "view" in cmds[0]
        assert any("mapq < 10" in arg for arg in cmds[0])
        # Second command is samtools fastq reading from stdin
        assert "fastq" in cmds[1]
        assert "-" in cmds[1]


# ---------------------------------------------------------------------------
# Integration tests – require samtools
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _require_samtools() -> None:
    """Skip integration tests if samtools is not installed."""
    try:
        _find_samtools()
    except FileNotFoundError:
        pytest.skip("samtools not available")


class TestExtractUnmappedOnly:
    """Extract only unmapped reads (default mode, no MAPQ filter)."""

    def test_produces_output_files(
        self, test_bam: Path, tmp_path: Path
    ) -> None:
        outputs = extract_reads(test_bam, tmp_path / "out")
        # At least one output file should exist
        assert len(outputs) > 0
        for path in outputs.values():
            assert path.exists()

    def test_correct_unmapped_count(
        self, test_bam: Path, tmp_path: Path
    ) -> None:
        outputs = extract_reads(test_bam, tmp_path / "out")
        total = sum(count_fastq_reads(p) for p in outputs.values())
        # Each unmapped pair produces 2 reads in the FASTQ output
        # We expect exactly the unmapped reads (10 pairs = 20 reads)
        assert total == EXPECTED_UNMAPPED_READS * 2

    def test_interleaved_mode(
        self, test_bam: Path, tmp_path: Path
    ) -> None:
        outputs = extract_reads(
            test_bam, tmp_path / "out_interleaved", interleaved=True
        )
        assert "interleaved" in outputs
        total = count_fastq_reads(outputs["interleaved"])
        assert total == EXPECTED_UNMAPPED_READS * 2


class TestExtractWithMAPQ:
    """Extract unmapped + low-MAPQ reads."""

    def test_more_reads_with_mapq_filter(
        self, test_bam: Path, tmp_path: Path
    ) -> None:
        out_unmapped = tmp_path / "unmapped_only"
        out_mapq = tmp_path / "with_mapq"

        outputs_unmapped = extract_reads(test_bam, out_unmapped)
        outputs_mapq = extract_reads(test_bam, out_mapq, mapq_threshold=10)

        total_unmapped = sum(
            count_fastq_reads(p) for p in outputs_unmapped.values()
        )
        total_mapq = sum(
            count_fastq_reads(p) for p in outputs_mapq.values()
        )
        # With MAPQ filter we should get *more* reads
        assert total_mapq > total_unmapped

    def test_mapq_extracts_expected_count(
        self, test_bam: Path, tmp_path: Path
    ) -> None:
        outputs = extract_reads(test_bam, tmp_path / "out", mapq_threshold=10)
        total = sum(count_fastq_reads(p) for p in outputs.values())
        # unmapped (10 pairs) + low-MAPQ (5 pairs) = 15 pairs = 30 reads
        expected = (EXPECTED_UNMAPPED_READS + EXPECTED_LOW_MAPQ_READS) * 2
        assert total == expected


class TestCLI:
    """Test the CLI entry point."""

    def test_version(self) -> None:
        from extract_unmapped.cli import main

        with pytest.raises(SystemExit) as exc_info:
            main(["--version"])
        assert exc_info.value.code == 0

    def test_extract_via_cli(
        self, test_bam: Path, tmp_path: Path
    ) -> None:
        from extract_unmapped.cli import main

        rc = main([str(test_bam), "-o", str(tmp_path / "cli_out")])
        assert rc == 0
        # Verify output was produced
        out_files = list((tmp_path / "cli_out").glob("*.fastq.gz"))
        assert len(out_files) > 0

    def test_missing_input_returns_error(self, tmp_path: Path) -> None:
        from extract_unmapped.cli import main

        rc = main([str(tmp_path / "no_such.bam"), "-o", str(tmp_path / "out")])
        assert rc == 1
