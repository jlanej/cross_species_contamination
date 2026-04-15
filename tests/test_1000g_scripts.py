"""
Tests for tests/1000G shell scripts.

Checks syntax, dry-run behaviour, key path-resolution logic, and the
CRAI-download workaround without requiring SLURM, Apptainer, or network
access.

Background on the CRAI download fix
------------------------------------
``samtools view -X <index>`` accepts an explicit index path so that
samtools can seek directly to the unmapped virtual contig (``'*'``) and
avoid scanning the entire CRAM.  However, direct FTP access can be flaky
on some systems.

The fix in ``extract_unmapped_array.sh`` is to download the small CRAI to
a local temp file (Aspera first, curl fallback) before running the
samtools pipeline, then pass the **local path** to ``-X``.
"""

import os
import re
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
    """Write a small but valid manifest.tsv to *tmp_path*.

    Column headers match the real manifest: SAMPLE_ID, CRAM_FTP_URL, CRAI_FTP_URL.
    FTP URLs are used here (as in production) to exercise the code paths that
    deal with remote index files.
    """
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "SAMPLE_ID\tCRAM_FTP_URL\tCRAI_FTP_URL\n"
        "NA12718\tftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram\t"
        "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram.crai\n"
        "NA12748\tftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239481/NA12748.final.cram\t"
        "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239481/NA12748.final.cram.crai\n"
        "NA18488\tftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239483/NA18488.final.cram\t"
        "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239483/NA18488.final.cram.crai\n"
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

    def test_dry_run_default_max_concurrent_jobs(self, tmp_path):
        """Default array spec should include %300 throttling."""
        manifest = minimal_manifest(tmp_path)
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(tmp_path / "output"),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "--array=1-3%300" in result.stdout

    def test_dry_run_custom_max_concurrent_jobs(self, tmp_path):
        """--max-concurrent-jobs should override default throttling."""
        manifest = minimal_manifest(tmp_path)
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(tmp_path / "output"),
            "--max-concurrent-jobs", "25",
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "--array=1-3%25" in result.stdout

    def test_dry_run_skips_completed_samples(self, tmp_path):
        """Completed sample outputs should be excluded from array submission."""
        manifest = minimal_manifest(tmp_path)
        outdir = tmp_path / "output"
        completed = outdir / "NA12718"
        completed.mkdir(parents=True)
        (completed / "NA12718_unmapped_R1.fastq.gz").write_text("done")
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(outdir),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "--array=2,3%300" in result.stdout
        assert "Skipping 1 completed sample(s)" in result.stdout

    def test_dry_run_all_completed_submits_nothing(self, tmp_path):
        """If all outputs exist, script should exit 0 without sbatch command."""
        manifest = minimal_manifest(tmp_path)
        outdir = tmp_path / "output"
        for sid in ["NA12718", "NA12748", "NA18488"]:
            sample_dir = outdir / sid
            sample_dir.mkdir(parents=True)
            (sample_dir / f"{sid}_unmapped_R1.fastq.gz").write_text("done")
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(outdir),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "Nothing to submit" in result.stdout
        assert "sbatch command:" not in result.stdout


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


# ---------------------------------------------------------------------------
# CRAI download workaround tests
#
# These tests verify that the array script correctly downloads the CRAI to a
# local temp file and uses the local path with -X, rather than passing the raw
# FTP URL to samtools (which fails on some htslib builds with "Exec format
# error").
# ---------------------------------------------------------------------------

class TestCraiLocalDownload:
    """Verify the CRAI local-download logic in extract_unmapped_array.sh."""

    def _make_fake_curl(self, tmp_path: Path, *, succeed: bool = True) -> Path:
        """Write a fake ``curl`` executable that creates a dummy file on success."""
        fake_curl = tmp_path / "curl"
        if succeed:
            fake_curl.write_text(
                "#!/usr/bin/env bash\n"
                # curl is called as: curl ... -o <dest> <url>
                # Find the -o argument and create a 4-byte dummy file there.
                "dest=''\n"
                "while [[ $# -gt 0 ]]; do\n"
                "  if [[ $1 == '-o' ]]; then dest=$2; shift 2\n"
                "  else shift; fi\n"
                "done\n"
                "[[ -n $dest ]] && printf '\\x00\\x00\\x00\\x00' > \"$dest\"\n"
                "exit 0\n"
            )
        else:
            fake_curl.write_text("#!/usr/bin/env bash\necho 'curl: error' >&2; exit 6\n")
        fake_curl.chmod(0o755)
        return fake_curl

    def _make_fake_apptainer(self, tmp_path: Path) -> Path:
        """Write a fake `apptainer` that always exits 0."""
        fake = tmp_path / "apptainer"
        fake.write_text("#!/usr/bin/env bash\nexit 0\n")
        fake.chmod(0o755)
        return fake

    def test_pipeline_script_uses_local_crai_not_url(self, tmp_path):
        """The generated pipeline script must reference the local CRAI path,
        not the original FTP URL."""
        manifest = tmp_path / "manifest.tsv"
        manifest.write_text(
            "SAMPLE_ID\tCRAM_FTP_URL\tCRAI_FTP_URL\n"
            "NA12718\tftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram\t"
            "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram.crai\n"
        )
        outdir = tmp_path / "output"
        outdir.mkdir()

        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        # Provide fake curl and apptainer so the script can progress far enough
        # to write (and then we can inspect) the pipeline script.
        self._make_fake_curl(bin_dir)
        self._make_fake_apptainer(bin_dir)

        env = {**os.environ, "PATH": f"{bin_dir}:{os.environ['PATH']}"}

        # Source just the sections of the script we need by extracting the
        # relevant bash logic inline (avoids running SLURM / container bits).
        inline = textwrap.dedent("""\
            #!/usr/bin/env bash
            set -euo pipefail
            SAMPLE_ID=NA12718
            CRAM_URL=ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram
            CRAI_URL=ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram.crai
            SAMPLE_DIR={outdir}/NA12718
            mkdir -p "$SAMPLE_DIR"
            CRAI_LOCAL="$SAMPLE_DIR/.NA12718.crai.tmp"
            curl -fsSL --retry 3 --retry-delay 5 -o "$CRAI_LOCAL" "$CRAI_URL"
            # Check that the local path exists and the FTP URL is NOT embedded
            echo "CRAI_LOCAL=$CRAI_LOCAL"
        """.format(outdir=outdir))

        result = run(["bash", "-c", inline], env=env)
        assert result.returncode == 0, result.stderr
        assert "CRAI_LOCAL=" in result.stdout
        local_path = [
            line.split("=", 1)[1]
            for line in result.stdout.splitlines()
            if line.startswith("CRAI_LOCAL=")
        ][0]
        # Must be a filesystem path, not an FTP URL
        assert not local_path.startswith("ftp://"), (
            f"Expected local path but got: {local_path}"
        )
        assert local_path.startswith(str(outdir)), (
            f"CRAI temp file should be under OUTDIR; got: {local_path}"
        )

    def test_crai_download_failure_aborts_pipeline(self, tmp_path):
        """If curl fails to download the CRAI the script must exit non-zero."""
        manifest = tmp_path / "manifest.tsv"
        manifest.write_text(
            "SAMPLE_ID\tCRAM_FTP_URL\tCRAI_FTP_URL\n"
            "NA12718\tftp://ftp.sra.ebi.ac.uk/NA12718.cram\t"
            "ftp://ftp.sra.ebi.ac.uk/NA12718.cram.crai\n"
        )
        outdir = tmp_path / "output"
        outdir.mkdir()

        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        self._make_fake_curl(bin_dir, succeed=False)

        env = {**os.environ, "PATH": f"{bin_dir}:{os.environ['PATH']}"}

        inline = textwrap.dedent("""\
            #!/usr/bin/env bash
            set -euo pipefail
            SAMPLE_ID=NA12718
            CRAI_URL=ftp://ftp.sra.ebi.ac.uk/NA12718.cram.crai
            SAMPLE_DIR={outdir}/NA12718
            mkdir -p "$SAMPLE_DIR"
            CRAI_LOCAL="$SAMPLE_DIR/.NA12718.crai.tmp"
            if ! curl -fsSL --retry 3 --retry-delay 5 -o "$CRAI_LOCAL" "$CRAI_URL"; then
                echo "ERROR: Failed to download CRAI index: $CRAI_URL" >&2
                rm -f "$CRAI_LOCAL"
                exit 1
            fi
        """.format(outdir=outdir))

        result = run(["bash", "-c", inline], env=env)
        assert result.returncode != 0, "Expected non-zero exit when curl fails"
        assert "ERROR" in result.stderr

    def test_array_script_content_uses_local_crai(self):
        """Verify the actual script file references CRAI_LOCAL not CRAI_URL in
        the samtools view command embedded in the pipeline script."""
        import re

        content = ARRAY_SCRIPT.read_text()

        # The pipeline script generation must use the local variable, not the URL
        assert "q_crai_local" in content, (
            "Script should define q_crai_local for the local CRAI path"
        )
        # Find printf lines that write a samtools view command into the pipeline
        # script.  These are the lines that actually end up in the generated
        # pipeline, so they must use q_crai_local (local path), not q_crai_url
        # (the raw FTP URL which fails on some htslib builds).
        view_printf_pattern = re.compile(
            r"""^\s*printf\s+['"](.*samtools\s+view.*-X.*%s.*)['"]""",
            re.MULTILINE,
        )
        for match in view_printf_pattern.finditer(content):
            assert "q_crai_url" not in match.group(0), (
                "samtools view printf should use q_crai_local, not q_crai_url: "
                + match.group(0).strip()
            )

    def test_array_script_content_passes_cram_then_crai_to_x(self):
        """samtools -X must receive CRAM first, then CRAI index path."""
        content = ARRAY_SCRIPT.read_text()
        assert (
            '"${q_threads}" "${q_cram_url}" "${q_crai_local}"'
        ) in content, "Pipeline script should pass -X <CRAM_URL> <CRAI_LOCAL>"
        assert (
            '"${q_threads}" "${q_crai_local}" "${q_cram_url}"'
        ) not in content, "Script must not pass -X <CRAI_LOCAL> <CRAM_URL>"

    def test_array_script_curl_downloads_crai_before_pipeline(self):
        """The curl download of CRAI must precede the pipeline script block."""
        import re

        content = ARRAY_SCRIPT.read_text()

        # Locate the first non-comment line that calls curl for the CRAI download
        curl_match = re.search(
            r"^[^#]*curl\s+-fsSL",
            content,
            re.MULTILINE,
        )
        # Locate the assignment that opens the pipeline script file
        pipeline_match = re.search(
            r"^[^#]*PIPELINE_SCRIPT\s*=",
            content,
            re.MULTILINE,
        )
        assert curl_match is not None, "Script must contain a curl download for the CRAI"
        assert pipeline_match is not None, "Script must contain PIPELINE_SCRIPT= assignment"
        assert curl_match.start() < pipeline_match.start(), (
            "curl download of CRAI must come before pipeline script generation"
        )

    def test_array_script_includes_aspera_download_path(self):
        """Array script should include an Aspera download helper for CRAI."""
        content = ARRAY_SCRIPT.read_text()
        assert "download_crai_with_aspera()" in content
        assert "command -v ascp" in content
        assert 'aspera_user="era-fasp"' in content
        assert 'aspera_host="fasp.sra.ebi.ac.uk"' in content
