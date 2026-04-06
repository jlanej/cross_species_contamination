"""Tests for the Kraken2 classification module (csc.classify).

Unit tests mock the kraken2 binary so they can run without Kraken2 installed.
Integration tests (marked with _require_kraken2) require kraken2 on PATH and
a mock database.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from unittest import mock

import pytest

from csc.classify.classify import (
    ClassificationResult,
    _find_kraken2,
    _get_kraken2_version,
    build_classify_command,
    classify_reads,
    validate_database,
    REQUIRED_DB_FILES,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def mock_db(tmp_path: Path) -> Path:
    """Create a minimal mock Kraken2 database directory."""
    db_dir = tmp_path / "mock_kraken2_db"
    db_dir.mkdir()
    for fname in REQUIRED_DB_FILES:
        (db_dir / fname).write_bytes(b"\x00" * 64)
    return db_dir


@pytest.fixture()
def mock_fastq(tmp_path: Path) -> Path:
    """Create a minimal FASTQ file for testing."""
    fq = tmp_path / "sample.unmapped.fastq.gz"
    # Write a minimal uncompressed FASTQ (not gzipped, but sufficient for
    # command construction and mocked execution tests).
    fq.write_text(
        "@read1\nACGTACGT\n+\nIIIIIIII\n"
        "@read2\nTGCATGCA\n+\nIIIIIIII\n"
    )
    return fq


@pytest.fixture()
def mock_fastq_pair(tmp_path: Path) -> tuple[Path, Path]:
    """Create a minimal paired-end FASTQ pair for testing."""
    r1 = tmp_path / "sample.unmapped.R1.fastq.gz"
    r2 = tmp_path / "sample.unmapped.R2.fastq.gz"
    r1.write_text("@read1/1\nACGTACGT\n+\nIIIIIIII\n")
    r2.write_text("@read1/2\nTGCATGCA\n+\nIIIIIIII\n")
    return r1, r2


# ---------------------------------------------------------------------------
# Unit tests — database validation
# ---------------------------------------------------------------------------

class TestValidateDatabase:
    def test_valid_database(self, mock_db: Path) -> None:
        result = validate_database(mock_db)
        assert result == mock_db.resolve()

    def test_missing_directory_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="not found"):
            validate_database(tmp_path / "nonexistent")

    def test_missing_files_raises(self, tmp_path: Path) -> None:
        db_dir = tmp_path / "incomplete_db"
        db_dir.mkdir()
        # Only create one of the required files
        (db_dir / "hash.k2d").write_bytes(b"\x00")
        with pytest.raises(ValueError, match="missing required files"):
            validate_database(db_dir)

    def test_file_not_directory_raises(self, tmp_path: Path) -> None:
        not_dir = tmp_path / "not_a_dir"
        not_dir.touch()
        with pytest.raises(FileNotFoundError, match="not found"):
            validate_database(not_dir)


# ---------------------------------------------------------------------------
# Unit tests — command building
# ---------------------------------------------------------------------------

class TestBuildClassifyCommand:
    def test_basic_command(self, mock_fastq: Path, mock_db: Path) -> None:
        with mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2"):
            cmd = build_classify_command(
                [mock_fastq],
                db=mock_db,
                output=Path("/tmp/out.txt"),
                report=Path("/tmp/report.txt"),
            )
        assert cmd[0] == "/usr/bin/kraken2"
        assert "--db" in cmd
        assert str(mock_db) in cmd
        assert "--output" in cmd
        assert "--report" in cmd
        assert "--confidence" in cmd
        assert "0.0" in cmd
        assert "--memory-mapping" not in cmd
        assert "--paired" not in cmd

    def test_memory_mapping_flag(self, mock_fastq: Path, mock_db: Path) -> None:
        with mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2"):
            cmd = build_classify_command(
                [mock_fastq],
                db=mock_db,
                output=Path("/tmp/out.txt"),
                report=Path("/tmp/report.txt"),
                memory_mapping=True,
            )
        assert "--memory-mapping" in cmd

    def test_paired_flag(self, mock_fastq_pair: tuple[Path, Path], mock_db: Path) -> None:
        r1, r2 = mock_fastq_pair
        with mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2"):
            cmd = build_classify_command(
                [r1, r2],
                db=mock_db,
                output=Path("/tmp/out.txt"),
                report=Path("/tmp/report.txt"),
                paired=True,
            )
        assert "--paired" in cmd
        assert str(r1) in cmd
        assert str(r2) in cmd

    def test_custom_confidence_and_threads(self, mock_fastq: Path, mock_db: Path) -> None:
        with mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2"):
            cmd = build_classify_command(
                [mock_fastq],
                db=mock_db,
                output=Path("/tmp/out.txt"),
                report=Path("/tmp/report.txt"),
                confidence=0.5,
                threads=8,
            )
        idx_conf = cmd.index("--confidence")
        assert cmd[idx_conf + 1] == "0.5"
        idx_threads = cmd.index("--threads")
        assert cmd[idx_threads + 1] == "8"

    def test_input_files_at_end(self, mock_fastq: Path, mock_db: Path) -> None:
        with mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2"):
            cmd = build_classify_command(
                [mock_fastq],
                db=mock_db,
                output=Path("/tmp/out.txt"),
                report=Path("/tmp/report.txt"),
            )
        # Input file should be the last argument
        assert cmd[-1] == str(mock_fastq)


# ---------------------------------------------------------------------------
# Unit tests — classify_reads with mocked subprocess
# ---------------------------------------------------------------------------

class TestClassifyReads:
    @mock.patch("csc.classify.classify.subprocess.run")
    @mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2")
    def test_successful_classification(
        self, mock_find: mock.Mock, mock_run: mock.Mock,
        mock_fastq: Path, mock_db: Path, tmp_path: Path,
    ) -> None:
        mock_run.return_value = mock.Mock(
            returncode=0,
            stderr="  2 sequences classified (100.00%)\n  0 sequences unclassified (0.00%)\n",
        )

        out_dir = tmp_path / "classify_out"
        result = classify_reads(
            [mock_fastq], out_dir, db=mock_db,
        )

        assert result["sample_id"] == "sample.unmapped"
        assert result["report"] == out_dir / "sample.unmapped.kraken2.report.txt"
        assert result["output"] == out_dir / "sample.unmapped.kraken2.output.txt"
        assert len(result["input_files"]) == 1
        # subprocess.run is called twice: once for --version, once for classification
        assert mock_run.call_count == 2

    @mock.patch("csc.classify.classify.subprocess.run")
    @mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2")
    def test_custom_sample_id(
        self, mock_find: mock.Mock, mock_run: mock.Mock,
        mock_fastq: Path, mock_db: Path, tmp_path: Path,
    ) -> None:
        mock_run.return_value = mock.Mock(returncode=0, stderr="")
        out_dir = tmp_path / "classify_out"

        result = classify_reads(
            [mock_fastq], out_dir, db=mock_db, sample_id="SAMPLE_001",
        )
        assert result["sample_id"] == "SAMPLE_001"
        assert "SAMPLE_001" in result["report"].name

    @mock.patch("csc.classify.classify.subprocess.run")
    @mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2")
    def test_paired_mode(
        self, mock_find: mock.Mock, mock_run: mock.Mock,
        mock_fastq_pair: tuple[Path, Path], mock_db: Path, tmp_path: Path,
    ) -> None:
        mock_run.return_value = mock.Mock(returncode=0, stderr="")
        r1, r2 = mock_fastq_pair
        out_dir = tmp_path / "classify_out"

        result = classify_reads(
            [r1, r2], out_dir, db=mock_db, paired=True,
        )
        assert len(result["input_files"]) == 2

        # Verify --paired was passed in the command
        call_args = mock_run.call_args
        cmd = call_args[0][0]
        assert "--paired" in cmd

    @mock.patch("csc.classify.classify.subprocess.run")
    @mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2")
    def test_memory_mapping_passed(
        self, mock_find: mock.Mock, mock_run: mock.Mock,
        mock_fastq: Path, mock_db: Path, tmp_path: Path,
    ) -> None:
        mock_run.return_value = mock.Mock(returncode=0, stderr="")
        out_dir = tmp_path / "classify_out"

        classify_reads(
            [mock_fastq], out_dir, db=mock_db, memory_mapping=True,
        )

        call_args = mock_run.call_args
        cmd = call_args[0][0]
        assert "--memory-mapping" in cmd

    @mock.patch("csc.classify.classify.subprocess.run")
    @mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2")
    def test_kraken2_failure_raises(
        self, mock_find: mock.Mock, mock_run: mock.Mock,
        mock_fastq: Path, mock_db: Path, tmp_path: Path,
    ) -> None:
        mock_run.return_value = mock.Mock(
            returncode=1,
            stderr="Error: database load failed",
        )
        out_dir = tmp_path / "classify_out"

        with pytest.raises(RuntimeError, match="Kraken2 failed"):
            classify_reads([mock_fastq], out_dir, db=mock_db)

    def test_missing_input_raises(self, mock_db: Path, tmp_path: Path) -> None:
        with mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2"):
            with pytest.raises(FileNotFoundError, match="Input file not found"):
                classify_reads(
                    [tmp_path / "missing.fastq.gz"],
                    tmp_path / "out",
                    db=mock_db,
                )

    def test_empty_input_raises(self, mock_db: Path, tmp_path: Path) -> None:
        with mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2"):
            with pytest.raises(ValueError, match="At least one input"):
                classify_reads([], tmp_path / "out", db=mock_db)

    def test_paired_wrong_count_raises(
        self, mock_fastq: Path, mock_db: Path, tmp_path: Path,
    ) -> None:
        with mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2"):
            with pytest.raises(ValueError, match="exactly 2 input files"):
                classify_reads(
                    [mock_fastq], tmp_path / "out",
                    db=mock_db, paired=True,
                )

    def test_invalid_db_raises(self, mock_fastq: Path, tmp_path: Path) -> None:
        with mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2"):
            with pytest.raises(FileNotFoundError, match="not found"):
                classify_reads(
                    [mock_fastq], tmp_path / "out",
                    db=tmp_path / "no_db",
                )

    @mock.patch("csc.classify.classify.subprocess.run")
    @mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2")
    def test_output_dir_created(
        self, mock_find: mock.Mock, mock_run: mock.Mock,
        mock_fastq: Path, mock_db: Path, tmp_path: Path,
    ) -> None:
        mock_run.return_value = mock.Mock(returncode=0, stderr="")
        out_dir = tmp_path / "nested" / "classify_out"

        classify_reads([mock_fastq], out_dir, db=mock_db)
        assert out_dir.is_dir()

    @mock.patch("csc.classify.classify.subprocess.run")
    @mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2")
    def test_sample_id_strips_fastq_gz(
        self, mock_find: mock.Mock, mock_run: mock.Mock,
        mock_db: Path, tmp_path: Path,
    ) -> None:
        fq = tmp_path / "reads.fastq.gz"
        fq.write_text("@r1\nACGT\n+\nIIII\n")
        mock_run.return_value = mock.Mock(returncode=0, stderr="")

        result = classify_reads([fq], tmp_path / "out", db=mock_db)
        assert result["sample_id"] == "reads"

    @mock.patch("csc.classify.classify.subprocess.run")
    @mock.patch("csc.classify.classify._find_kraken2", return_value="/usr/bin/kraken2")
    def test_confidence_passed(
        self, mock_find: mock.Mock, mock_run: mock.Mock,
        mock_fastq: Path, mock_db: Path, tmp_path: Path,
    ) -> None:
        mock_run.return_value = mock.Mock(returncode=0, stderr="")

        classify_reads(
            [mock_fastq], tmp_path / "out",
            db=mock_db, confidence=0.5,
        )

        call_args = mock_run.call_args
        cmd = call_args[0][0]
        idx = cmd.index("--confidence")
        assert cmd[idx + 1] == "0.5"


# ---------------------------------------------------------------------------
# Unit tests — helper functions
# ---------------------------------------------------------------------------

class TestFindKraken2:
    def test_not_found_raises(self) -> None:
        with mock.patch("csc.classify.classify.shutil.which", return_value=None):
            with pytest.raises(FileNotFoundError, match="kraken2 not found"):
                _find_kraken2()

    def test_found_returns_path(self) -> None:
        with mock.patch("csc.classify.classify.shutil.which", return_value="/usr/bin/kraken2"):
            assert _find_kraken2() == "/usr/bin/kraken2"


class TestGetKraken2Version:
    def test_parses_version_line(self) -> None:
        with mock.patch("csc.classify.classify.subprocess.run") as mock_run:
            mock_run.return_value = mock.Mock(
                stdout="Kraken version 2.1.3\nSome other line\n",
                stderr="",
            )
            result = _get_kraken2_version("/usr/bin/kraken2")
            assert "2.1.3" in result

    def test_returns_unknown_on_failure(self) -> None:
        with mock.patch("csc.classify.classify.subprocess.run", side_effect=OSError):
            result = _get_kraken2_version("/usr/bin/kraken2")
            assert result == "unknown"


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------

class TestCLI:
    def test_version(self) -> None:
        from csc.classify.cli import main

        with pytest.raises(SystemExit) as exc_info:
            main(["--version"])
        assert exc_info.value.code == 0

    def test_missing_required_args(self) -> None:
        from csc.classify.cli import main

        with pytest.raises(SystemExit) as exc_info:
            main([])
        assert exc_info.value.code != 0

    def test_missing_input_returns_error(self, tmp_path: Path) -> None:
        from csc.classify.cli import main

        rc = main([
            str(tmp_path / "no_such.fastq.gz"),
            "--db", str(tmp_path / "db"),
            "-o", str(tmp_path / "out"),
        ])
        assert rc == 1

    @mock.patch("csc.classify.cli.classify_reads")
    def test_successful_cli_run(
        self, mock_classify: mock.Mock,
        mock_fastq: Path, mock_db: Path, tmp_path: Path,
    ) -> None:
        from csc.classify.cli import main

        mock_classify.return_value = {
            "report": tmp_path / "out" / "sample.kraken2.report.txt",
            "output": tmp_path / "out" / "sample.kraken2.output.txt",
            "sample_id": "sample",
            "input_files": [str(mock_fastq)],
        }

        rc = main([
            str(mock_fastq),
            "--db", str(mock_db),
            "-o", str(tmp_path / "out"),
        ])
        assert rc == 0
        mock_classify.assert_called_once()

    @mock.patch("csc.classify.cli.classify_reads")
    def test_cli_passes_params(
        self, mock_classify: mock.Mock,
        mock_fastq: Path, mock_db: Path, tmp_path: Path,
    ) -> None:
        from csc.classify.cli import main

        mock_classify.return_value = {
            "report": tmp_path / "report.txt",
            "output": tmp_path / "output.txt",
            "sample_id": "sample",
            "input_files": [str(mock_fastq)],
        }

        main([
            str(mock_fastq),
            "--db", str(mock_db),
            "-o", str(tmp_path / "out"),
            "--confidence", "0.3",
            "--threads", "4",
            "--memory-mapping",
            "--sample-id", "MY_SAMPLE",
        ])

        call_kwargs = mock_classify.call_args
        assert call_kwargs[1]["confidence"] == 0.3
        assert call_kwargs[1]["threads"] == 4
        assert call_kwargs[1]["memory_mapping"] is True
        assert call_kwargs[1]["sample_id"] == "MY_SAMPLE"

    @mock.patch("csc.classify.cli.classify_reads", side_effect=RuntimeError("kraken2 fail"))
    def test_cli_handles_error(
        self, mock_classify: mock.Mock,
        mock_fastq: Path, mock_db: Path, tmp_path: Path,
    ) -> None:
        from csc.classify.cli import main

        rc = main([
            str(mock_fastq),
            "--db", str(mock_db),
            "-o", str(tmp_path / "out"),
        ])
        assert rc == 1
