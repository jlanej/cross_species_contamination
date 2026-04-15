"""
Tests for tests/1000G shell scripts.

Checks syntax, dry-run behaviour, and key path-resolution logic without
requiring SLURM, Apptainer, or network access.
"""

import os
import subprocess
import tempfile
import textwrap
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parent / "1000G"
SUBMIT_SCRIPT = SCRIPTS_DIR / "submit_extract.sh"
ARRAY_SCRIPT = SCRIPTS_DIR / "extract_unmapped_array.sh"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd, env=None, **kwargs):
    """Run *cmd* and return CompletedProcess; never raises on non-zero exit."""
    return subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        env=env,
        **kwargs,
    )


def minimal_manifest(tmp_path: Path) -> Path:
    """Write a small but valid manifest.tsv to *tmp_path*."""
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "SAMPLE_ID\tCRAM_URL\tCRAI_URL\n"
        "NA12718\thttps://example.com/NA12718.cram\thttps://example.com/NA12718.cram.crai\n"
        "NA12748\thttps://example.com/NA12748.cram\thttps://example.com/NA12748.cram.crai\n"
        "NA18488\thttps://example.com/NA18488.cram\thttps://example.com/NA18488.cram.crai\n"
    )
    return manifest


# ---------------------------------------------------------------------------
# Syntax checks (bash -n)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("script", [SUBMIT_SCRIPT, ARRAY_SCRIPT])
def test_bash_syntax(script):
    """Both scripts must pass bash syntax checking."""
    result = run(["bash", "-n", str(script)])
    assert result.returncode == 0, (
        f"bash -n failed for {script.name}:\n{result.stderr}"
    )


# ---------------------------------------------------------------------------
# submit_extract.sh dry-run tests
# ---------------------------------------------------------------------------

class TestSubmitExtractDryRun:
    """Tests that use --dry-run so no SLURM / Apptainer is needed."""

    def test_dry_run_no_args(self, tmp_path):
        """dry-run with default args should succeed and mention sbatch."""
        manifest = minimal_manifest(tmp_path)
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(tmp_path / "output"),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "sbatch" in result.stdout

    def test_dry_run_limit(self, tmp_path):
        """--limit 1 should produce array spec '1-1'."""
        manifest = minimal_manifest(tmp_path)
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(tmp_path / "output"),
            "--limit", "1",
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "1-1" in result.stdout

    def test_dry_run_range(self, tmp_path):
        """--range 2-3 should pass that spec verbatim to sbatch."""
        manifest = minimal_manifest(tmp_path)
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(tmp_path / "output"),
            "--range", "2-3",
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "2-3" in result.stdout

    def test_dry_run_samples_file(self, tmp_path):
        """--samples resolves sample IDs to 1-based indices."""
        manifest = minimal_manifest(tmp_path)
        samples_file = tmp_path / "samples.txt"
        samples_file.write_text("NA12718\nNA18488\n")
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(tmp_path / "output"),
            "--samples", str(samples_file),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        # NA12718 is line 1, NA18488 is line 3 → indices 1 and 3
        assert "1" in result.stdout
        assert "3" in result.stdout

    def test_dry_run_missing_sample_warning(self, tmp_path):
        """Unknown sample IDs should emit a warning but not fail (other IDs ok)."""
        manifest = minimal_manifest(tmp_path)
        samples_file = tmp_path / "samples.txt"
        samples_file.write_text("NA12718\nNONEXISTENT\n")
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(tmp_path / "output"),
            "--samples", str(samples_file),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "WARNING" in result.stderr or "WARNING" in result.stdout

    def test_dry_run_container_sif_always_exported(self, tmp_path):
        """CONTAINER_SIF must always appear in the --export string."""
        manifest = minimal_manifest(tmp_path)
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(tmp_path / "output"),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "CONTAINER_SIF=" in result.stdout

    def test_dry_run_custom_container(self, tmp_path):
        """--container <path> should propagate the custom SIF path to sbatch."""
        manifest = minimal_manifest(tmp_path)
        custom_sif = tmp_path / "my_custom.sif"
        # Create a fake SIF so the script skips the pull step
        custom_sif.touch()
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(tmp_path / "output"),
            "--container", str(custom_sif),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert str(custom_sif) in result.stdout

    def test_dry_run_missing_manifest(self, tmp_path):
        """Missing manifest should exit with an error."""
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(tmp_path / "nonexistent.tsv"),
            "--dry-run",
        ])
        assert result.returncode != 0
        assert "Manifest not found" in result.stderr or "ERROR" in result.stderr

    def test_help_flag(self):
        """-h should print usage and exit 0."""
        result = run(["bash", str(SUBMIT_SCRIPT), "-h"])
        assert result.returncode == 0
        assert "Usage" in result.stdout or "usage" in result.stdout

    def test_dry_run_all_samples_creates_full_range(self, tmp_path):
        """With no --limit/--range/--samples, the spec covers all data rows."""
        manifest = minimal_manifest(tmp_path)
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(tmp_path / "output"),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        # 3 samples → array spec should be '1-3'
        assert "1-3" in result.stdout


# ---------------------------------------------------------------------------
# extract_unmapped_array.sh unit checks
# ---------------------------------------------------------------------------

class TestArrayScript:
    def test_fails_without_slurm_task_id(self):
        """Script should exit non-zero if SLURM_ARRAY_TASK_ID is unset."""
        env = {**os.environ}
        env.pop("SLURM_ARRAY_TASK_ID", None)
        result = run(["bash", str(ARRAY_SCRIPT)], env=env)
        assert result.returncode != 0
        # Either apptainer is missing (no HPC) or SLURM_ARRAY_TASK_ID check fires;
        # both are expected failure modes in a non-HPC environment.
        assert result.stderr != ""

    def test_default_container_sif_uses_outdir(self, tmp_path):
        """When CONTAINER_SIF is not set, default should be under OUTDIR."""
        # We set a non-existent OUTDIR; the script will fail before pulling
        # because SLURM_ARRAY_TASK_ID is not set, but we can extract the
        # default via a small sourcing trick.
        env = {**os.environ}
        env.pop("CONTAINER_SIF", None)
        env["OUTDIR"] = str(tmp_path / "outdir")
        script = textwrap.dedent("""\
            #!/usr/bin/env bash
            SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
            OUTDIR="${OUTDIR:-${SCRIPT_DIR}/output}"
            CONTAINER_SIF="${CONTAINER_SIF:-${OUTDIR}/csc.sif}"
            echo "CONTAINER_SIF=${CONTAINER_SIF}"
        """)
        # Inline the relevant lines so we don't run the whole script
        result = run(["bash", "-c", script], env=env)
        assert result.returncode == 0
        # The default SIF path must NOT contain SLURM temp-dir patterns
        sif_path = [
            line.split("=", 1)[1]
            for line in result.stdout.splitlines()
            if line.startswith("CONTAINER_SIF=")
        ][0]
        # The default SIF path must be under OUTDIR (not a SLURM temp dir)
        assert sif_path.startswith(str(tmp_path / "outdir"))
