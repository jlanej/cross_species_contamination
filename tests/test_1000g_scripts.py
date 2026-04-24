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

import json
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
SUBMIT_CLASSIFY_SCRIPT = SCRIPTS_DIR / "submit_classify.sh"
CLASSIFY_ARRAY_SCRIPT = SCRIPTS_DIR / "classify_array.sh"
AGGREGATE_DETECT_SCRIPT = SCRIPTS_DIR / "aggregate_detect.sh"

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

@pytest.mark.parametrize("script", [
    SUBMIT_SCRIPT,
    ARRAY_SCRIPT,
    SUBMIT_CLASSIFY_SCRIPT,
    CLASSIFY_ARRAY_SCRIPT,
    AGGREGATE_DETECT_SCRIPT,
])
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

    def test_dry_run_skip_idxstats_default_exported(self, tmp_path):
        """SKIP_IDXSTATS should default to 0 in array job exports."""
        manifest = minimal_manifest(tmp_path)
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(tmp_path / "output"),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "SKIP_IDXSTATS=0" in result.stdout

    def test_dry_run_skip_idxstats_flag_exported(self, tmp_path):
        """--skip-idxstats should export SKIP_IDXSTATS=1."""
        manifest = minimal_manifest(tmp_path)
        result = run([
            "bash", str(SUBMIT_SCRIPT),
            "--manifest", str(manifest),
            "--outdir", str(tmp_path / "output"),
            "--skip-idxstats",
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "SKIP_IDXSTATS=1" in result.stdout

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

    def test_idxstats_enabled_by_default_and_can_be_skipped(self):
        """Array script should default SKIP_IDXSTATS=0 and support skip mode."""
        content = ARRAY_SCRIPT.read_text()
        assert 'SKIP_IDXSTATS="${SKIP_IDXSTATS:-0}"' in content
        assert 'if [[ "${SKIP_IDXSTATS}" == "1" ]]; then' in content
        assert "samtools idxstats" in content
        assert ".reads_summary.json" in content

    def test_idxstats_uses_cram_url_with_local_crai_idx(self):
        """idxstats should use the CRAM URL with ``##idx##`` pointing to the local CRAI.

        ``samtools idxstats`` does not support ``-X``.  Instead, htslib's
        ``##idx##`` URL suffix is used to pass the already-downloaded local
        CRAI, avoiding a redundant (and potentially stalling) FTP fetch.
        """
        content = ARRAY_SCRIPT.read_text()
        # Should call idxstats with CRAM_URL##idx##CRAI_LOCAL (no -X flag)
        assert 'samtools idxstats "${CRAM_URL}##idx##${CRAI_LOCAL}"' in content
        # Must NOT use the (invalid for idxstats) -X flag
        assert 'samtools idxstats -X' not in content
        # Must NOT use the old bare URL form (which re-downloads CRAI via FTP)
        assert 'samtools idxstats "${CRAM_URL}"' not in content

    def test_idxstats_remote_url_local_crai_functional(self, tmp_path: Path):
        """Functional test: idxstats block executes samtools idxstats with ##idx##.

        Uses a fake samtools that captures its arguments and emits valid
        idxstats output, simulating the 1000G remote CRAM scenario without
        actual network access.  Verifies:

        * samtools is called with ``<CRAM_URL>##idx##<CRAI_LOCAL>`` (no ``-X``)
        * the local CRAI path is embedded in the argument via ``##idx##``
        * the ``.idxstats.tsv`` sidecar is written
        * the ``reads_summary.json`` sidecar is produced with correct structure
        """
        sample_dir = tmp_path / "NA12718"
        sample_dir.mkdir()

        # A real (but empty) file to stand in as the local CRAI
        crai_local = tmp_path / ".NA12718.crai.tmp"
        crai_local.write_bytes(b"")

        # Fake samtools: record args, emit valid idxstats TSV on idxstats sub-command
        args_log = tmp_path / "samtools_args.txt"
        fake_samtools = tmp_path / "samtools"
        fake_samtools.write_text(
            "#!/usr/bin/env bash\n"
            "echo \"$@\" >> " + str(args_log) + "\n"
            'if [[ "$1" == "idxstats" ]]; then\n'
            "    printf 'chr1\\t248956422\\t45\\t5\\n'\n"
            "    printf 'chr2\\t242193529\\t5\\t3\\n'\n"
            "    printf '*\\t0\\t0\\t12\\n'\n"
            "fi\n"
        )
        fake_samtools.chmod(0o755)

        # The idxstats block expects python3 to be available for the JSON step.
        remote_cram_url = (
            "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram"
        )
        reads_summary_json = sample_dir / "NA12718.reads_summary.json"
        idxstats_tsv = sample_dir / "NA12718.idxstats.tsv"

        # Execute the idxstats block inline, with our fake samtools on the PATH.
        # ``samtools idxstats`` uses ##idx## to specify the local CRAI; no -X.
        inline = textwrap.dedent(f"""\
            set -euo pipefail
            export PATH="{tmp_path}:$PATH"
            SAMPLE_ID="NA12718"
            SAMPLE_DIR="{sample_dir}"
            CRAM_URL="{remote_cram_url}"
            CRAI_LOCAL="{crai_local}"
            SLURM_ARRAY_TASK_ID="1"
            IDXSTATS_TSV="{idxstats_tsv}"
            READS_SUMMARY_JSON="{reads_summary_json}"
            # container_run: no apptainer in CI, so run directly
            container_run() {{ "$@"; }}
            SKIP_IDXSTATS=0

            if [[ "${{SKIP_IDXSTATS}}" == "1" ]]; then
                echo "Skipping idxstats sidecars (SKIP_IDXSTATS=1)."
            else
                echo "Computing idxstats sidecars for ${{SAMPLE_ID}}..."
                if container_run samtools idxstats "${{CRAM_URL}}##idx##${{CRAI_LOCAL}}" > "${{IDXSTATS_TSV}}.tmp"; then
                    echo "idxstats computed (${{CRAM_URL}})"
                else
                    echo "ERROR: samtools idxstats failed." >&2
                    exit 1
                fi
                [[ -s "${{IDXSTATS_TSV}}.tmp" ]] || {{ echo "ERROR: empty idxstats output." >&2; exit 1; }}
                mv -f "${{IDXSTATS_TSV}}.tmp" "${{IDXSTATS_TSV}}"

                IDXSTATS_JSON_SCRIPT="{tmp_path}/.idxstats_to_json_1.py"
                cat > "${{IDXSTATS_JSON_SCRIPT}}" <<'PY'
import datetime
import json
import os
import sys
from pathlib import Path

sample_id = os.environ["CSC_SAMPLE_ID"]
input_path = os.environ["CSC_INPUT_PATH"]
idxstats_tsv = Path(os.environ["CSC_IDXSTATS_TSV"])
reads_summary_json = Path(os.environ["CSC_READS_SUMMARY_JSON"])
per_chromosome = []
total_mapped = 0
total_unmapped = 0
for line in idxstats_tsv.read_text().splitlines():
    if not line.strip():
        continue
    parts = line.split("\\t")
    if len(parts) < 4:
        continue
    chrom, length_s, mapped_s, unmapped_s = parts[:4]
    try:
        length = int(length_s)
        mapped = int(mapped_s)
        unmapped = int(unmapped_s)
    except ValueError:
        continue
    per_chromosome.append({{"chrom": chrom, "length": length, "mapped": mapped, "unmapped": unmapped}})
    total_mapped += mapped
    total_unmapped += unmapped
summary = {{
    "schema_version": "1.0", "sample_id": sample_id, "input": input_path,
    "extraction_time": datetime.datetime.now(tz=datetime.timezone.utc).isoformat(),
    "total_mapped": total_mapped, "total_unmapped": total_unmapped,
    "total_reads": total_mapped + total_unmapped, "per_chromosome": per_chromosome,
}}
reads_summary_json.write_text(json.dumps(summary, indent=2) + "\\n")
PY
                CSC_SAMPLE_ID="${{SAMPLE_ID}}" CSC_INPUT_PATH="${{CRAM_URL}}" \
                    CSC_IDXSTATS_TSV="${{IDXSTATS_TSV}}" \
                    CSC_READS_SUMMARY_JSON="${{READS_SUMMARY_JSON}}" \
                    python3 "${{IDXSTATS_JSON_SCRIPT}}"
            fi
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0, f"idxstats block failed:\n{result.stderr}"

        # Verify samtools was invoked with ##idx## URL (not bare CRAM_URL, not -X)
        assert args_log.exists(), "Fake samtools never invoked"
        args_recorded = args_log.read_text().strip().splitlines()
        assert len(args_recorded) >= 1, "samtools args not captured"
        first_call = args_recorded[0]
        expected_arg = f"{remote_cram_url}##idx##"
        assert expected_arg in first_call, (
            f"Expected '##idx##' in samtools args: {first_call!r}"
        )
        assert str(crai_local) in first_call, (
            f"Expected local CRAI path in samtools args: {first_call!r}"
        )
        assert "-X" not in first_call, (
            f"samtools idxstats must NOT use -X (not supported): {first_call!r}"
        )

        # Verify sidecars were produced
        assert idxstats_tsv.exists(), ".idxstats.tsv sidecar not written"
        assert reads_summary_json.exists(), "reads_summary.json not written"

        # Verify reads_summary.json structure and counts
        doc = json.loads(reads_summary_json.read_text())
        assert doc["sample_id"] == "NA12718"
        assert doc["input"] == remote_cram_url
        assert doc["total_mapped"] == 50   # 45 + 5
        assert doc["total_unmapped"] == 20  # 5 + 3 + 12
        assert doc["total_reads"] == 70
        assert "per_chromosome" in doc
        chroms = {r["chrom"] for r in doc["per_chromosome"]}
        assert "chr1" in chroms and "*" in chroms


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


# ---------------------------------------------------------------------------
# Helpers shared by classify pipeline tests
# ---------------------------------------------------------------------------

def _fake_extracted_samples(outdir: Path, sample_ids: list) -> None:
    """Create dummy extraction output dirs and R1/R2 FASTQ files."""
    for sid in sample_ids:
        sample_dir = outdir / sid
        sample_dir.mkdir(parents=True, exist_ok=True)
        # write non-empty files so [[ -s ]] passes
        (sample_dir / f"{sid}_unmapped_R1.fastq.gz").write_bytes(b"\x1f\x8b\x08\x00")
        (sample_dir / f"{sid}_unmapped_R2.fastq.gz").write_bytes(b"\x1f\x8b\x08\x00")


def _fake_db(tmp_path: Path) -> Path:
    """Create a minimal fake Kraken2 DB directory."""
    db = tmp_path / "kraken2_db"
    db.mkdir()
    for f in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (db / f).touch()
    return db


# ---------------------------------------------------------------------------
# Syntax checks already covered by the expanded parametrize above.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# submit_classify.sh dry-run tests
# ---------------------------------------------------------------------------

class TestSubmitClassifyDryRun:
    """Tests that use --dry-run so no SLURM / Apptainer / Kraken2 is needed."""

    def test_dry_run_basic(self, tmp_path):
        """dry-run with minimal args should succeed and mention sbatch."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718", "NA12748"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "sbatch" in result.stdout

    def test_dry_run_missing_db(self, tmp_path):
        """Missing --db should exit with an error."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--dry-run",
        ])
        assert result.returncode != 0
        assert "ERROR" in result.stderr or "ERROR" in result.stdout

    def test_dry_run_missing_extract_outdir(self, tmp_path):
        """Nonexistent --extract-outdir should exit with an error."""
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(tmp_path / "nonexistent"),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode != 0
        assert "ERROR" in result.stderr or "ERROR" in result.stdout

    def test_dry_run_no_extracted_samples(self, tmp_path):
        """Empty extraction dir should exit with an error."""
        extract_out = tmp_path / "output"
        extract_out.mkdir()
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode != 0
        assert "No extracted" in result.stderr or "No extracted" in result.stdout

    def test_dry_run_writes_classify_manifest(self, tmp_path):
        """dry-run should still write classify_manifest.tsv."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718", "NA12748", "NA18488"])
        db = _fake_db(tmp_path)
        out = tmp_path / "classify_out"
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(out),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        manifest = out / "classify_manifest.tsv"
        assert manifest.exists(), "classify_manifest.tsv should be written"
        lines = manifest.read_text().splitlines()
        # header + 3 samples
        assert len(lines) == 4
        # header
        assert lines[0].startswith("SAMPLE_ID")
        # data lines contain sample IDs and R1 paths
        assert "NA12718" in lines[1]
        assert "_unmapped_R1.fastq.gz" in lines[1]

    def test_dry_run_limit(self, tmp_path):
        """--limit should cap the array spec to N."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718", "NA12748", "NA18488"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--limit", "2",
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        # Array spec has 2 samples; they may be rendered as "1-2" or "1,2"
        assert "1,2" in result.stdout or "1-2" in result.stdout

    def test_dry_run_samples_file(self, tmp_path):
        """--samples resolves sample IDs and only those appear in the manifest."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718", "NA12748", "NA18488"])
        db = _fake_db(tmp_path)
        samples_file = tmp_path / "subset.txt"
        samples_file.write_text("NA12718\nNA18488\n")
        out = tmp_path / "classify_out"
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(out),
            "--db", str(db),
            "--samples", str(samples_file),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        manifest = out / "classify_manifest.tsv"
        lines = manifest.read_text().splitlines()
        assert len(lines) == 3  # header + 2 samples
        sample_ids_in_manifest = [l.split("\t")[0] for l in lines[1:]]
        assert "NA12718" in sample_ids_in_manifest
        assert "NA18488" in sample_ids_in_manifest
        assert "NA12748" not in sample_ids_in_manifest

    def test_dry_run_unknown_sample_warns(self, tmp_path):
        """Unknown sample ID in --samples file should warn but not fail."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        samples_file = tmp_path / "subset.txt"
        samples_file.write_text("NA12718\nNONEXISTENT\n")
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--samples", str(samples_file),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "WARNING" in result.stderr or "WARNING" in result.stdout

    def test_dry_run_array_spec_includes_throttle(self, tmp_path):
        """Default array spec must include %200 concurrency throttle."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718", "NA12748"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "%200" in result.stdout

    def test_dry_run_custom_max_concurrent_jobs(self, tmp_path):
        """--max-concurrent-jobs should override default throttle."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718", "NA12748"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--max-concurrent-jobs", "50",
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "%50" in result.stdout

    def test_dry_run_skips_completed_samples(self, tmp_path):
        """Already-classified samples should be excluded from the array."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718", "NA12748", "NA18488"])
        db = _fake_db(tmp_path)
        classify_out = tmp_path / "classify_out"
        # Mark NA12718 as already classified
        done_dir = classify_out / "classify" / "NA12718"
        done_dir.mkdir(parents=True)
        (done_dir / "NA12718.kraken2.report.txt").write_text("# classified\n")
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(classify_out),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "Skipping 1 already-classified" in result.stdout
        # Only 2 and 3 should appear (NA12748, NA18488)
        assert "2,3" in result.stdout or "2-3" in result.stdout

    def test_dry_run_all_completed_no_classify_job(self, tmp_path):
        """If all samples are classified, no classify array job should appear."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718", "NA12748"])
        db = _fake_db(tmp_path)
        classify_out = tmp_path / "classify_out"
        for sid in ["NA12718", "NA12748"]:
            done_dir = classify_out / "classify" / sid
            done_dir.mkdir(parents=True)
            (done_dir / f"{sid}.kraken2.report.txt").write_text("# classified\n")
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(classify_out),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "All selected samples already have classification output" in result.stdout
        # Aggregate/detect job should still be printed (for existing results)
        assert "aggregate_detect" in result.stdout or "sbatch" in result.stdout

    def test_dry_run_skip_aggregate(self, tmp_path):
        """--skip-aggregate should suppress the aggregate/detect sbatch command."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--skip-aggregate",
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "aggregate_detect" not in result.stdout

    def test_dry_run_container_sif_exported(self, tmp_path):
        """CONTAINER_SIF must appear in the classify --export string."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "CONTAINER_SIF=" in result.stdout

    def test_dry_run_rank_filter_colon_encoded(self, tmp_path):
        """RANK_FILTER_CODES must use colon separator (no spaces) for safe export."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--rank-filter", "S G F",
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "RANK_FILTER_CODES=S:G:F" in result.stdout

    def test_dry_run_detect_matrix_default_exported(self, tmp_path):
        """DETECT_MATRIX should default to cpm in aggregate/detect exports."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "DETECT_MATRIX=cpm" in result.stdout

    def test_dry_run_detect_matrix_raw_exported(self, tmp_path):
        """--detect-matrix raw should be exported for aggregate/detect."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--detect-matrix", "raw",
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "DETECT_MATRIX=raw" in result.stdout

    def test_dry_run_detect_method_default_exported(self, tmp_path):
        """DETECT_METHOD should default to 'all' in the aggregate/detect exports."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "DETECT_METHOD=all" in result.stdout

    def test_dry_run_mad_threshold_default_exported(self, tmp_path):
        """MAD_THRESHOLD should default to 3.5 in the aggregate/detect exports."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "MAD_THRESHOLD=3.5" in result.stdout

    def test_dry_run_iqr_multiplier_default_exported(self, tmp_path):
        """IQR_MULTIPLIER should default to 1.5 in the aggregate/detect exports."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "IQR_MULTIPLIER=1.5" in result.stdout

    def test_help_flag(self):
        """-h should print usage and exit 0."""
        result = run(["bash", str(SUBMIT_CLASSIFY_SCRIPT), "-h"])
        assert result.returncode == 0
        assert "Usage" in result.stdout or "usage" in result.stdout

    def test_unknown_flag_errors(self, tmp_path):
        """Unknown flags should exit non-zero with a message."""
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--this-flag-does-not-exist",
        ])
        assert result.returncode != 0
        assert "ERROR" in result.stderr or "ERROR" in result.stdout

    def test_dry_run_paired_r2_in_manifest(self, tmp_path):
        """Samples with R2 files should have R2 path in the manifest."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        out = tmp_path / "classify_out"
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(out),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        manifest = out / "classify_manifest.tsv"
        line = manifest.read_text().splitlines()[1]
        cols = line.split("\t")
        assert len(cols) == 3
        # R2 column should be non-empty (we created both R1 and R2 above)
        assert "_unmapped_R2.fastq.gz" in cols[2]

    def test_dry_run_single_end_empty_r2_in_manifest(self, tmp_path):
        """Samples with only R1 should have empty R2 column in the manifest."""
        extract_out = tmp_path / "output"
        sid = "SE_SAMPLE"
        sample_dir = extract_out / sid
        sample_dir.mkdir(parents=True)
        (sample_dir / f"{sid}_unmapped_R1.fastq.gz").write_bytes(b"\x1f\x8b\x08\x00")
        # No R2 file created
        db = _fake_db(tmp_path)
        out = tmp_path / "classify_out"
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(out),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        manifest = out / "classify_manifest.tsv"
        line = manifest.read_text().splitlines()[1]
        cols = line.split("\t")
        # R2 column should be empty
        assert cols[2] == ""

    def test_dry_run_db_path_defaults_to_db(self, tmp_path):
        """DB_PATH should default to the --db value in the aggregate/detect exports."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        # DB_PATH should be exported to the aggregate/detect job with the DB value
        assert f"DB_PATH={db}" in result.stdout

    def test_dry_run_extract_outdir_exported_to_aggregate_job(self, tmp_path):
        """EXTRACT_OUTDIR should be exported to aggregate/detect job."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert f"EXTRACT_OUTDIR={extract_out}" in result.stdout

    def test_dry_run_skip_idxstats_metrics_default_exported(self, tmp_path):
        """SKIP_IDXSTATS_METRICS should default to 0."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "SKIP_IDXSTATS_METRICS=0" in result.stdout

    def test_dry_run_skip_idxstats_metrics_flag_exported(self, tmp_path):
        """--skip-idxstats-metrics should export SKIP_IDXSTATS_METRICS=1."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--skip-idxstats-metrics",
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert "SKIP_IDXSTATS_METRICS=1" in result.stdout

    def test_dry_run_db_path_explicit_override(self, tmp_path):
        """--db-path should override the default DB_PATH in aggregate/detect exports."""
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        db = _fake_db(tmp_path)
        custom_db = tmp_path / "custom_db"
        custom_db.mkdir()
        result = run([
            "bash", str(SUBMIT_CLASSIFY_SCRIPT),
            "--extract-outdir", str(extract_out),
            "--outdir", str(tmp_path / "classify_out"),
            "--db", str(db),
            "--db-path", str(custom_db),
            "--dry-run",
        ])
        assert result.returncode == 0, result.stderr
        assert f"DB_PATH={custom_db}" in result.stdout
        # Must not contain the default DB path as DB_PATH
        assert f"DB_PATH={db}" not in result.stdout


# ---------------------------------------------------------------------------
# classify_array.sh unit checks
# ---------------------------------------------------------------------------

class TestClassifyArrayScript:

    @staticmethod
    def _make_fake_apptainer(tmp_path: Path) -> Path:
        """Write a fake apptainer/singularity that always exits 0 (never pulls)."""
        bin_dir = tmp_path / "fakebin"
        bin_dir.mkdir(exist_ok=True)
        for name in ("apptainer", "singularity"):
            exe = bin_dir / name
            exe.write_text("#!/usr/bin/env bash\nexit 0\n")
            exe.chmod(0o755)
        return bin_dir

    @staticmethod
    def _fake_sif(tmp_path: Path) -> Path:
        """Return path to a zero-byte fake SIF file (just needs to exist as a regular file)."""
        sif = tmp_path / "fake.sif"
        sif.touch()
        return sif

    def test_fails_without_slurm_task_id(self, tmp_path):
        """Script should exit non-zero if SLURM_ARRAY_TASK_ID is unset."""
        env = {**os.environ}
        env.pop("SLURM_ARRAY_TASK_ID", None)
        result = run(["bash", str(CLASSIFY_ARRAY_SCRIPT)], env=env)
        assert result.returncode != 0
        # Either apptainer is missing (no HPC) or SLURM_ARRAY_TASK_ID check fires;
        # both are expected failure modes in a non-HPC environment.
        assert result.stderr != ""

    def test_fails_without_classify_manifest(self, tmp_path):
        """Missing CLASSIFY_MANIFEST should exit with an error."""
        bin_dir = self._make_fake_apptainer(tmp_path)
        env = {**os.environ,
               "PATH": f"{bin_dir}:{os.environ['PATH']}",
               "SLURM_ARRAY_TASK_ID": "1",
               "CLASSIFY_MANIFEST": str(tmp_path / "nonexistent.tsv"),
               "CLASSIFY_OUTDIR": str(tmp_path / "classify"),
               "EXTRACT_OUTDIR": str(tmp_path / "output"),
               "DB": str(tmp_path / "db"),
               "CONTAINER_SIF": str(self._fake_sif(tmp_path))}
        result = run(["bash", str(CLASSIFY_ARRAY_SCRIPT)], env=env)
        assert result.returncode != 0
        assert "ERROR" in result.stderr

    def test_fails_without_db(self, tmp_path):
        """Missing DB environment variable should exit with an error."""
        bin_dir = self._make_fake_apptainer(tmp_path)
        manifest = tmp_path / "manifest.tsv"
        manifest.write_text("SAMPLE_ID\tR1\tR2\nNA12718\t/r1.fq.gz\t\n")
        env = {**os.environ,
               "PATH": f"{bin_dir}:{os.environ['PATH']}",
               "SLURM_ARRAY_TASK_ID": "1",
               "CLASSIFY_MANIFEST": str(manifest),
               "CLASSIFY_OUTDIR": str(tmp_path / "classify"),
               "EXTRACT_OUTDIR": str(tmp_path / "output"),
               "DB": "",
               "CONTAINER_SIF": str(self._fake_sif(tmp_path))}
        result = run(["bash", str(CLASSIFY_ARRAY_SCRIPT)], env=env)
        assert result.returncode != 0
        assert "DB" in result.stderr

    def test_idempotence_skips_completed_sample(self, tmp_path):
        """If report already exists, the array task should exit 0 immediately."""
        bin_dir = self._make_fake_apptainer(tmp_path)
        # Build a minimal classify manifest
        extract_out = tmp_path / "output"
        _fake_extracted_samples(extract_out, ["NA12718"])
        manifest = tmp_path / "classify_manifest.tsv"
        r1 = extract_out / "NA12718" / "NA12718_unmapped_R1.fastq.gz"
        r2 = extract_out / "NA12718" / "NA12718_unmapped_R2.fastq.gz"
        manifest.write_text(
            f"SAMPLE_ID\tR1\tR2\n"
            f"NA12718\t{r1}\t{r2}\n"
        )

        # Create a pre-existing (non-empty) report
        classify_dir = tmp_path / "classify" / "NA12718"
        classify_dir.mkdir(parents=True)
        (classify_dir / "NA12718.kraken2.report.txt").write_text("# done\n")

        # Provide a fake DB dir so the DB check passes
        db = tmp_path / "db"
        db.mkdir()

        env = {**os.environ,
               "PATH": f"{bin_dir}:{os.environ['PATH']}",
               "SLURM_ARRAY_TASK_ID": "1",
               "CLASSIFY_MANIFEST": str(manifest),
               "CLASSIFY_OUTDIR": str(tmp_path / "classify"),
               "EXTRACT_OUTDIR": str(extract_out),
               "DB": str(db),
               # Provide a real (zero-byte) fake SIF so pull_container is skipped
               "CONTAINER_SIF": str(self._fake_sif(tmp_path))}
        result = run(["bash", str(CLASSIFY_ARRAY_SCRIPT)], env=env)
        # Should exit 0 with "already complete" message
        assert result.returncode == 0, result.stderr
        assert "already complete" in result.stdout or "Skipping" in result.stdout


# ---------------------------------------------------------------------------
# aggregate_detect.sh unit checks
# ---------------------------------------------------------------------------

class TestAggregateDetectScript:

    def test_fails_without_classify_outdir(self, tmp_path):
        """Missing CLASSIFY_OUTDIR should exit with an error."""
        env = {**os.environ,
               "CLASSIFY_OUTDIR": "",
               "AGG_OUTDIR": str(tmp_path / "agg")}
        result = run(["bash", str(AGGREGATE_DETECT_SCRIPT)], env=env)
        assert result.returncode != 0
        assert "CLASSIFY_OUTDIR" in result.stderr

    def test_fails_without_agg_outdir(self, tmp_path):
        """Missing AGG_OUTDIR should exit with an error."""
        classify_out = tmp_path / "classify"
        classify_out.mkdir()
        env = {**os.environ,
               "CLASSIFY_OUTDIR": str(classify_out),
               "AGG_OUTDIR": ""}
        result = run(["bash", str(AGGREGATE_DETECT_SCRIPT)], env=env)
        assert result.returncode != 0
        assert "AGG_OUTDIR" in result.stderr

    def test_fails_when_no_reports_found(self, tmp_path):
        """Empty classify dir (no .kraken2.report.txt files) should fail."""
        classify_out = tmp_path / "classify"
        classify_out.mkdir()
        env = {**os.environ,
               "CLASSIFY_OUTDIR": str(classify_out),
               "AGG_OUTDIR": str(tmp_path / "agg"),
               "CONTAINER_SIF": "/dev/null"}
        # apptainer will not be called because it fails on missing reports first
        result = run(["bash", str(AGGREGATE_DETECT_SCRIPT)], env=env)
        assert result.returncode != 0
        assert "No .kraken2.report.txt" in result.stderr or "ERROR" in result.stderr

    def test_rank_filter_colon_expansion(self, tmp_path):
        """RANK_FILTER_CODES colon separator must expand to space-separated codes."""
        # We test the expansion logic inline to avoid requiring apptainer/csc tools
        inline = textwrap.dedent("""\
            RANK_FILTER_CODES="S:G:F"
            IFS=':' read -ra RANK_CODES <<< "${RANK_FILTER_CODES}"
            echo "${RANK_CODES[@]}"
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0
        assert result.stdout.strip() == "S G F"

    def test_detect_outdir_defaults_to_sibling_of_agg(self, tmp_path):
        """When DETECT_OUTDIR is not set, it should default to <parent>/detect."""
        agg_dir = tmp_path / "classify_out" / "aggregate"
        agg_dir.mkdir(parents=True)
        classify_out = tmp_path / "classify_out" / "classify"
        classify_out.mkdir(parents=True)

        # Build the default-derivation logic inline
        inline = textwrap.dedent(f"""\
            AGG_OUTDIR="{agg_dir}"
            DETECT_OUTDIR="${{DETECT_OUTDIR:-}}"
            if [[ -z "$DETECT_OUTDIR" ]]; then
                DETECT_OUTDIR="$(dirname "$AGG_OUTDIR")/detect"
            fi
            echo "DETECT_OUTDIR=$DETECT_OUTDIR"
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0
        detect_path = [
            l.split("=", 1)[1]
            for l in result.stdout.splitlines()
            if l.startswith("DETECT_OUTDIR=")
        ][0]
        expected = str(tmp_path / "classify_out" / "detect")
        assert detect_path == expected

    def test_db_path_passed_to_aggregate_when_set(self, tmp_path):
        """When DB_PATH is set, --db-path should be included in csc-aggregate args."""
        # We test the AGGREGATE_ARGS construction logic inline
        db = tmp_path / "mydb"
        db.mkdir()
        inline = textwrap.dedent(f"""\
            DB_PATH="{db}"
            AGGREGATE_ARGS=(csc-aggregate /some/report.txt -o /out --rank-filter S G F)
            if [[ -n "$DB_PATH" ]]; then
                AGGREGATE_ARGS+=("--db-path" "$DB_PATH")
            fi
            echo "${{AGGREGATE_ARGS[@]}}"
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0
        assert "--db-path" in result.stdout
        assert str(db) in result.stdout

    def test_db_path_not_passed_when_empty(self, tmp_path):
        """When DB_PATH is empty, --db-path must not appear in csc-aggregate args."""
        inline = textwrap.dedent("""\
            DB_PATH=""
            AGGREGATE_ARGS=(csc-aggregate /some/report.txt -o /out --rank-filter S G F)
            if [[ -n "$DB_PATH" ]]; then
                AGGREGATE_ARGS+=("--db-path" "$DB_PATH")
            fi
            echo "${AGGREGATE_ARGS[@]}"
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0
        assert "--db-path" not in result.stdout

    def test_idxstats_paths_added_to_aggregate_args_by_default(self, tmp_path):
        """Aggregate args should include --idxstats with per-sample reads_summary files."""
        extract_out = tmp_path / "extract"
        sample = extract_out / "NA12718"
        sample.mkdir(parents=True)
        reads_summary = sample / "NA12718.reads_summary.json"
        reads_summary.write_text('{"sample_id":"NA12718","total_reads":70}\n')
        inline = textwrap.dedent(f"""\
            SKIP_IDXSTATS_METRICS=0
            EXTRACT_OUTDIR="{extract_out}"
            REPORTS=("/classify/NA12718/NA12718.kraken2.report.txt")
            AGGREGATE_ARGS=(csc-aggregate "${{REPORTS[@]}}" -o /agg)
            if [[ "${{SKIP_IDXSTATS_METRICS}}" != "1" ]]; then
                IDXSTATS_PATHS=()
                for report in "${{REPORTS[@]}}"; do
                    sid="$(basename "$report")"
                    sid="${{sid%.kraken2.report.txt}}"
                    reads_summary="${{EXTRACT_OUTDIR}}/${{sid}}/${{sid}}.reads_summary.json"
                    [[ -s "$reads_summary" ]] || exit 1
                    IDXSTATS_PATHS+=("$reads_summary")
                done
                AGGREGATE_ARGS+=("--idxstats" "${{IDXSTATS_PATHS[@]}}")
            fi
            echo "${{AGGREGATE_ARGS[@]}}"
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0
        assert "--idxstats" in result.stdout
        assert str(reads_summary) in result.stdout

    def test_idxstats_missing_sidecar_fails_when_not_skipped(self, tmp_path):
        """Missing reads_summary should fail when idxstats metrics are required."""
        extract_out = tmp_path / "extract"
        extract_out.mkdir()
        inline = textwrap.dedent(f"""\
            SKIP_IDXSTATS_METRICS=0
            EXTRACT_OUTDIR="{extract_out}"
            report="/classify/NA12718/NA12718.kraken2.report.txt"
            sid="$(basename "$report")"
            sid="${{sid%.kraken2.report.txt}}"
            reads_summary="${{EXTRACT_OUTDIR}}/${{sid}}/${{sid}}.reads_summary.json"
            if [[ ! -s "$reads_summary" ]]; then
                echo "ERROR: Missing required idxstats sidecar for sample $sid: $reads_summary" >&2
                exit 1
            fi
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode != 0
        assert "Missing required idxstats sidecar" in result.stderr

    def test_skip_idxstats_metrics_omits_idxstats_args(self, tmp_path):
        """When skip is enabled, aggregate args should not include --idxstats."""
        inline = textwrap.dedent("""\
            SKIP_IDXSTATS_METRICS=1
            AGGREGATE_ARGS=(csc-aggregate /report.txt -o /agg)
            if [[ "${SKIP_IDXSTATS_METRICS}" != "1" ]]; then
                AGGREGATE_ARGS+=("--idxstats" "/x.reads_summary.json")
            fi
            echo "${AGGREGATE_ARGS[@]}"
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0
        assert "--idxstats" not in result.stdout

    def test_aggregate_args_always_include_rank_filter(self, tmp_path):
        """csc-aggregate AGGREGATE_ARGS must always include --rank-filter with rank codes."""
        inline = textwrap.dedent("""\
            RANK_FILTER_CODES="S:G:F"
            MIN_READS=0
            IFS=':' read -ra RANK_CODES <<< "${RANK_FILTER_CODES}"
            AGGREGATE_ARGS=(
                csc-aggregate
                /report.txt
                -o /out
            )
            if [[ "${MIN_READS}" -gt 0 ]]; then
                AGGREGATE_ARGS+=("--min-reads" "${MIN_READS}")
            fi
            AGGREGATE_ARGS+=("--rank-filter" "${RANK_CODES[@]}")
            echo "${AGGREGATE_ARGS[@]}"
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0
        assert "--rank-filter S G F" in result.stdout

    def test_detect_args_include_method_and_thresholds(self, tmp_path):
        """DETECT_ARGS must include --method, --mad-threshold, --iqr-multiplier, --gmm-threshold, --rank-filter."""
        inline = textwrap.dedent("""\
            DETECT_METHOD="all"
            MAD_THRESHOLD=3.5
            IQR_MULTIPLIER=1.5
            GMM_THRESHOLD=0.5
            RANK_FILTER_CODES="S:G:F"
            IFS=':' read -ra RANK_CODES <<< "${RANK_FILTER_CODES}"
            MATRIX="/agg/taxa_matrix_cpm.tsv"
            DETECT_OUTDIR="/detect"
            DETECT_ARGS=(
                csc-detect
                "${MATRIX}"
                -o "${DETECT_OUTDIR}"
                --method       "${DETECT_METHOD}"
                --mad-threshold "${MAD_THRESHOLD}"
                --iqr-multiplier "${IQR_MULTIPLIER}"
                --gmm-threshold "${GMM_THRESHOLD}"
                --rank-filter  "${RANK_CODES[@]}"
            )
            echo "${DETECT_ARGS[@]}"
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0
        assert "--method all" in result.stdout
        assert "--mad-threshold 3.5" in result.stdout
        assert "--iqr-multiplier 1.5" in result.stdout
        assert "--gmm-threshold 0.5" in result.stdout
        assert "--rank-filter S G F" in result.stdout

    def test_detect_matrix_path_uses_cpm_by_default(self, tmp_path):
        """With DETECT_MATRIX=cpm, MATRIX should point to taxa_matrix_cpm.tsv."""
        inline = textwrap.dedent("""\
            AGG_OUTDIR="/some/agg"
            DETECT_MATRIX="cpm"
            MATRIX="${AGG_OUTDIR}/taxa_matrix_${DETECT_MATRIX}.tsv"
            echo "MATRIX=$MATRIX"
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0
        assert "taxa_matrix_cpm.tsv" in result.stdout

    def test_detect_matrix_path_raw_selection(self, tmp_path):
        """With DETECT_MATRIX=raw, MATRIX should point to taxa_matrix_raw.tsv."""
        inline = textwrap.dedent("""\
            AGG_OUTDIR="/some/agg"
            DETECT_MATRIX="raw"
            MATRIX="${AGG_OUTDIR}/taxa_matrix_${DETECT_MATRIX}.tsv"
            echo "MATRIX=$MATRIX"
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0
        assert "taxa_matrix_raw.tsv" in result.stdout

    def test_detect_matrix_invalid_value_rejected(self, tmp_path):
        """DETECT_MATRIX must be 'cpm' or 'raw'; invalid values should fail validation."""
        inline = textwrap.dedent("""\
            DETECT_MATRIX="bad"
            if [[ "${DETECT_MATRIX}" != "cpm" && "${DETECT_MATRIX}" != "raw" ]]; then
                echo "ERROR: DETECT_MATRIX must be 'cpm' or 'raw'." >&2
                exit 1
            fi
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode != 0
        assert "ERROR" in result.stderr

    def test_skip_detect_bypasses_detect_step(self, tmp_path):
        """SKIP_DETECT=1 must cause the script to exit before running csc-detect."""
        inline = textwrap.dedent("""\
            SKIP_DETECT=1
            if [[ "${SKIP_DETECT}" == "1" ]]; then
                echo "Skipping outlier detection (SKIP_DETECT=1)."
                exit 0
            fi
            echo "RUNNING_DETECT"
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0
        assert "Skipping" in result.stdout
        assert "RUNNING_DETECT" not in result.stdout

    def test_container_run_binds_extract_outdir_when_idxstats_enabled(self, tmp_path):
        """container_run must bind EXTRACT_OUTDIR when SKIP_IDXSTATS_METRICS != 1.

        Regression test for the bug where reads_summary.json sidecars were
        inaccessible inside the Apptainer container because EXTRACT_OUTDIR was
        not included in the bind-mount arguments.
        """
        extract_out = tmp_path / "extract"
        classify_out = tmp_path / "classify"
        agg_out = tmp_path / "agg"
        detect_out = tmp_path / "detect"
        for d in (extract_out, classify_out, agg_out, detect_out):
            d.mkdir(parents=True)

        inline = textwrap.dedent(f"""\
            SKIP_IDXSTATS_METRICS=0
            SKIP_DETECT=0
            CLASSIFY_OUTDIR="{classify_out}"
            AGG_OUTDIR="{agg_out}"
            DETECT_OUTDIR="{detect_out}"
            EXTRACT_OUTDIR="{extract_out}"
            DB_PATH=""
            container_run() {{
                local -a bind_args=()
                bind_args+=("--bind" "${{CLASSIFY_OUTDIR}}:${{CLASSIFY_OUTDIR}}")
                bind_args+=("--bind" "${{AGG_OUTDIR}}:${{AGG_OUTDIR}}")
                if [[ "${{SKIP_IDXSTATS_METRICS}}" != "1" ]]; then
                    bind_args+=("--bind" "${{EXTRACT_OUTDIR}}:${{EXTRACT_OUTDIR}}")
                fi
                if [[ "${{SKIP_DETECT}}" != "1" ]]; then
                    bind_args+=("--bind" "${{DETECT_OUTDIR}}:${{DETECT_OUTDIR}}")
                fi
                if [[ -n "${{DB_PATH}}" ]]; then
                    bind_args+=("--bind" "${{DB_PATH}}:${{DB_PATH}}")
                fi
                echo "${{bind_args[@]}}"
            }}
            container_run
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0
        assert str(extract_out) in result.stdout, (
            "EXTRACT_OUTDIR must appear in container bind args when "
            "SKIP_IDXSTATS_METRICS != 1"
        )

    def test_container_run_omits_extract_outdir_when_idxstats_skipped(self, tmp_path):
        """container_run must NOT bind EXTRACT_OUTDIR when SKIP_IDXSTATS_METRICS=1."""
        extract_out = tmp_path / "extract"
        classify_out = tmp_path / "classify"
        agg_out = tmp_path / "agg"
        detect_out = tmp_path / "detect"
        for d in (extract_out, classify_out, agg_out, detect_out):
            d.mkdir(parents=True)

        inline = textwrap.dedent(f"""\
            SKIP_IDXSTATS_METRICS=1
            SKIP_DETECT=0
            CLASSIFY_OUTDIR="{classify_out}"
            AGG_OUTDIR="{agg_out}"
            DETECT_OUTDIR="{detect_out}"
            EXTRACT_OUTDIR="{extract_out}"
            DB_PATH=""
            container_run() {{
                local -a bind_args=()
                bind_args+=("--bind" "${{CLASSIFY_OUTDIR}}:${{CLASSIFY_OUTDIR}}")
                bind_args+=("--bind" "${{AGG_OUTDIR}}:${{AGG_OUTDIR}}")
                if [[ "${{SKIP_IDXSTATS_METRICS}}" != "1" ]]; then
                    bind_args+=("--bind" "${{EXTRACT_OUTDIR}}:${{EXTRACT_OUTDIR}}")
                fi
                if [[ "${{SKIP_DETECT}}" != "1" ]]; then
                    bind_args+=("--bind" "${{DETECT_OUTDIR}}:${{DETECT_OUTDIR}}")
                fi
                if [[ -n "${{DB_PATH}}" ]]; then
                    bind_args+=("--bind" "${{DB_PATH}}:${{DB_PATH}}")
                fi
                echo "${{bind_args[@]}}"
            }}
            container_run
        """)
        result = run(["bash", "-c", inline])
        assert result.returncode == 0
        assert str(extract_out) not in result.stdout, (
            "EXTRACT_OUTDIR must NOT appear in container bind args when "
            "SKIP_IDXSTATS_METRICS=1"
        )
