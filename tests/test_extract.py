"""Tests for the unmapped-read extraction pipeline.

The synthetic test BAM (generated in conftest.py) contains:
  - 20 well-mapped paired-end reads  (MAPQ 60)
  - 10 unmapped paired-end reads     (flag 4, simulated contaminant)
  -  5 low-MAPQ paired-end reads     (MAPQ 3)

These counts are used to verify that extraction is correct.
"""

from __future__ import annotations

import csv
import gzip
import json
import logging
import subprocess
from pathlib import Path

import pytest

from csc.extract.extract import (
    build_extract_command,
    extract_reads,
    _find_samtools,
    _resolve_reference,
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

    def test_valid_cram_extension(self, test_cram: Path, test_reference: Path) -> None:
        _validate_input(test_cram, test_reference)  # should not raise

    def test_cram_missing_reference_raises(self, test_cram: Path, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="Reference file not found"):
            _validate_input(test_cram, tmp_path / "no_ref.fa")


class TestResolveReference:
    def test_explicit_reference_returned(self, test_reference: Path, test_bam: Path) -> None:
        assert _resolve_reference(test_reference, test_bam) == test_reference

    def test_env_ref_path(
        self, test_reference: Path, test_bam: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("REF_PATH", str(test_reference))
        result = _resolve_reference(None, test_bam)
        assert result == test_reference

    def test_env_ref_path_nonexistent_ignored(
        self, tmp_path: Path, test_bam: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("REF_PATH", str(tmp_path / "missing.fa"))
        result = _resolve_reference(None, test_bam)
        assert result is None

    def test_cram_without_reference_raises(self, tmp_path: Path) -> None:
        cram = tmp_path / "sample.cram"
        cram.touch()
        with pytest.raises(FileNotFoundError, match="CRAM input requires"):
            _resolve_reference(None, cram)

    def test_explicit_beats_env(
        self, test_reference: Path, test_bam: Path, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        other_ref = tmp_path / "other.fa"
        other_ref.touch()
        monkeypatch.setenv("REF_PATH", str(other_ref))
        result = _resolve_reference(test_reference, test_bam)
        assert result == test_reference


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
        result = extract_reads(test_bam, tmp_path / "out")
        outputs = result["files"]
        # At least one output file should exist
        assert len(outputs) > 0
        for path in outputs.values():
            assert path.exists()

    def test_correct_unmapped_count(
        self, test_bam: Path, tmp_path: Path
    ) -> None:
        result = extract_reads(test_bam, tmp_path / "out")
        outputs = result["files"]
        total = sum(count_fastq_reads(p) for p in outputs.values())
        # Each unmapped pair produces 2 reads in the FASTQ output
        # We expect exactly the unmapped reads (10 pairs = 20 reads)
        assert total == EXPECTED_UNMAPPED_READS * 2

    def test_result_contains_stats(
        self, test_bam: Path, tmp_path: Path
    ) -> None:
        result = extract_reads(test_bam, tmp_path / "out_stats")
        assert "read_count" in result
        assert "sample_id" in result
        assert "input" in result
        assert result["read_count"] == EXPECTED_UNMAPPED_READS * 2

    def test_interleaved_mode(
        self, test_bam: Path, tmp_path: Path
    ) -> None:
        result = extract_reads(
            test_bam, tmp_path / "out_interleaved", interleaved=True
        )
        outputs = result["files"]
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

        result_unmapped = extract_reads(test_bam, out_unmapped)
        result_mapq = extract_reads(test_bam, out_mapq, mapq_threshold=10)

        total_unmapped = sum(
            count_fastq_reads(p) for p in result_unmapped["files"].values()
        )
        total_mapq = sum(
            count_fastq_reads(p) for p in result_mapq["files"].values()
        )
        # With MAPQ filter we should get *more* reads
        assert total_mapq > total_unmapped

    def test_mapq_extracts_expected_count(
        self, test_bam: Path, tmp_path: Path
    ) -> None:
        result = extract_reads(test_bam, tmp_path / "out", mapq_threshold=10)
        total = sum(count_fastq_reads(p) for p in result["files"].values())
        # unmapped (10 pairs) + low-MAPQ (5 pairs) = 15 pairs = 30 reads
        expected = (EXPECTED_UNMAPPED_READS + EXPECTED_LOW_MAPQ_READS) * 2
        assert total == expected


# ---------------------------------------------------------------------------
# CRAM integration tests
# ---------------------------------------------------------------------------

class TestExtractCRAM:
    """Extract reads from CRAM files with reference."""

    def test_cram_unmapped_only(
        self, test_cram: Path, test_reference: Path, tmp_path: Path
    ) -> None:
        result = extract_reads(
            test_cram, tmp_path / "cram_out", reference=test_reference
        )
        outputs = result["files"]
        assert len(outputs) > 0
        total = sum(count_fastq_reads(p) for p in outputs.values())
        assert total == EXPECTED_UNMAPPED_READS * 2

    def test_cram_with_mapq(
        self, test_cram: Path, test_reference: Path, tmp_path: Path
    ) -> None:
        result = extract_reads(
            test_cram, tmp_path / "cram_mapq", reference=test_reference,
            mapq_threshold=10,
        )
        total = sum(count_fastq_reads(p) for p in result["files"].values())
        expected = (EXPECTED_UNMAPPED_READS + EXPECTED_LOW_MAPQ_READS) * 2
        assert total == expected

    def test_cram_interleaved(
        self, test_cram: Path, test_reference: Path, tmp_path: Path
    ) -> None:
        result = extract_reads(
            test_cram, tmp_path / "cram_interleaved",
            reference=test_reference, interleaved=True,
        )
        assert "interleaved" in result["files"]
        total = count_fastq_reads(result["files"]["interleaved"])
        assert total == EXPECTED_UNMAPPED_READS * 2

    def test_cram_via_ref_path_env(
        self, test_cram: Path, test_reference: Path, tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        monkeypatch.setenv("REF_PATH", str(test_reference))
        result = extract_reads(test_cram, tmp_path / "cram_env")
        outputs = result["files"]
        total = sum(count_fastq_reads(p) for p in outputs.values())
        assert total == EXPECTED_UNMAPPED_READS * 2


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------

class TestCLI:
    """Test the CLI entry point."""

    def test_version(self) -> None:
        from csc.extract.cli import main

        with pytest.raises(SystemExit) as exc_info:
            main(["--version"])
        assert exc_info.value.code == 0

    def test_extract_via_cli(
        self, test_bam: Path, tmp_path: Path
    ) -> None:
        from csc.extract.cli import main

        rc = main([str(test_bam), "-o", str(tmp_path / "cli_out")])
        assert rc == 0
        # Verify output was produced
        out_files = list((tmp_path / "cli_out").glob("*.fastq.gz"))
        assert len(out_files) > 0

    def test_missing_input_returns_error(self, tmp_path: Path) -> None:
        from csc.extract.cli import main

        rc = main([str(tmp_path / "no_such.bam"), "-o", str(tmp_path / "out")])
        assert rc == 1

    def test_unreadable_input_returns_error(self, tmp_path: Path) -> None:
        from csc.extract.cli import main
        import os

        bad = tmp_path / "unreadable.bam"
        bad.touch()
        os.chmod(bad, 0o000)
        try:
            rc = main([str(bad), "-o", str(tmp_path / "out")])
            assert rc == 1
        finally:
            os.chmod(bad, 0o644)

    def test_batch_missing_fails_early(self, test_bam: Path, tmp_path: Path) -> None:
        """If any file in a batch is missing, fail before processing any."""
        from csc.extract.cli import main

        rc = main([
            str(test_bam),
            str(tmp_path / "no_such.bam"),
            "-o", str(tmp_path / "out"),
        ])
        assert rc == 1
        # Output dir should not have been created for the valid file
        # because the CLI fails early
        out_files = list((tmp_path / "out").glob("*.fastq.gz"))
        assert len(out_files) == 0


class TestCLIJsonLog:
    """Test JSON structured logging via --json-log."""

    def test_json_log_flag(
        self, test_bam: Path, tmp_path: Path, capfd: pytest.CaptureFixture[str]
    ) -> None:
        from csc.extract.cli import main

        # Reset root logger handlers for clean test
        root = logging.getLogger()
        for h in root.handlers[:]:
            root.removeHandler(h)

        rc = main([str(test_bam), "-o", str(tmp_path / "json_out"), "--json-log", "-v"])
        assert rc == 0
        stderr = capfd.readouterr().err
        # At least one JSON log line should be present
        json_lines = [
            line for line in stderr.strip().splitlines() if line.startswith("{")
        ]
        assert len(json_lines) > 0
        parsed = json.loads(json_lines[0])
        assert "timestamp" in parsed
        assert "level" in parsed
        assert "message" in parsed


class TestCLISummary:
    """Test batch summary TSV via --summary."""

    def test_summary_tsv_produced(
        self, test_bam: Path, tmp_path: Path
    ) -> None:
        from csc.extract.cli import main

        summary_path = tmp_path / "summary.tsv"
        rc = main([
            str(test_bam), "-o", str(tmp_path / "sum_out"),
            "--summary", str(summary_path),
        ])
        assert rc == 0
        assert summary_path.exists()

        with open(summary_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 1
        assert rows[0]["status"] == "OK"
        assert int(rows[0]["read_count"]) == EXPECTED_UNMAPPED_READS * 2

    def test_summary_tsv_includes_errors(self, tmp_path: Path) -> None:
        """Summary should record failed files when processing errors occur."""
        from csc.extract.cli import main

        # Create a valid-looking but corrupt BAM file
        bad = tmp_path / "corrupt.bam"
        bad.write_bytes(b"\x00" * 100)

        summary_path = tmp_path / "summary.tsv"
        rc = main([
            str(bad), "-o", str(tmp_path / "sum_out"),
            "--summary", str(summary_path),
        ])
        assert rc == 1
        assert summary_path.exists()

        with open(summary_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 1
        assert rows[0]["status"] == "FAILED"
