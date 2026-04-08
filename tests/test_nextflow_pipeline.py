"""Tests for the Nextflow pipeline configuration and structure.

Validates that all Nextflow pipeline files exist, are well-formed,
and reference the expected processes and parameters.  These are
static / structural tests that do **not** require Nextflow to be
installed.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

# ── Helpers ──────────────────────────────────────────────────────────────────

NEXTFLOW_DIR = Path(__file__).resolve().parent.parent / "nextflow"
MODULES_DIR = NEXTFLOW_DIR / "modules"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


# ── File existence ───────────────────────────────────────────────────────────


class TestPipelineFilesExist:
    """All required Nextflow files must be present."""

    def test_main_nf_exists(self):
        assert (NEXTFLOW_DIR / "main.nf").is_file()

    def test_nextflow_config_exists(self):
        assert (NEXTFLOW_DIR / "nextflow.config").is_file()

    @pytest.mark.parametrize(
        "module",
        ["extract.nf", "classify.nf", "aggregate.nf", "detect.nf", "summary.nf"],
    )
    def test_module_exists(self, module: str):
        assert (MODULES_DIR / module).is_file(), f"Module {module} not found"


# ── main.nf structure ────────────────────────────────────────────────────────


class TestMainWorkflow:
    """Validate the main workflow file."""

    @pytest.fixture(autouse=True)
    def _load(self):
        self.content = _read(NEXTFLOW_DIR / "main.nf")

    def test_dsl2_enabled(self):
        assert "nextflow.enable.dsl = 2" in self.content

    def test_includes_all_modules(self):
        for proc in [
            "EXTRACT_UNMAPPED",
            "CLASSIFY_READS",
            "AGGREGATE_REPORTS",
            "DETECT_OUTLIERS",
            "PIPELINE_SUMMARY",
        ]:
            assert proc in self.content, f"Process {proc} not included"

    def test_input_csv_param(self):
        assert "params.input_csv" in self.content

    def test_kraken2_db_param(self):
        assert "params.kraken2_db" in self.content

    def test_outdir_param(self):
        assert "params.outdir" in self.content

    def test_validation_input_csv(self):
        assert "input_csv" in self.content
        # Should contain an error message about missing input_csv
        assert re.search(r"error.*input_csv", self.content, re.IGNORECASE)

    def test_validation_kraken2_db(self):
        assert re.search(r"error.*kraken2_db", self.content, re.IGNORECASE)

    def test_workflow_on_complete(self):
        assert "workflow.onComplete" in self.content

    def test_workflow_on_error(self):
        assert "workflow.onError" in self.content

    def test_stage_ordering(self):
        """Processes should appear in the right order in the workflow."""
        extract_pos = self.content.index("EXTRACT_UNMAPPED(")
        classify_pos = self.content.index("CLASSIFY_READS(")
        aggregate_pos = self.content.index("AGGREGATE_REPORTS(")
        detect_pos = self.content.index("DETECT_OUTLIERS(")
        summary_pos = self.content.index("PIPELINE_SUMMARY(")
        assert extract_pos < classify_pos < aggregate_pos < detect_pos < summary_pos


# ── Module processes ─────────────────────────────────────────────────────────


class TestExtractModule:
    """Validate the extract module."""

    @pytest.fixture(autouse=True)
    def _load(self):
        self.content = _read(MODULES_DIR / "extract.nf")

    def test_process_defined(self):
        assert "process EXTRACT_UNMAPPED" in self.content

    def test_calls_csc_extract(self):
        assert "csc-extract" in self.content

    def test_emits_fastqs(self):
        assert "emit: fastqs" in self.content

    def test_configurable_cpus(self):
        assert "params.extract_cpus" in self.content

    def test_configurable_memory(self):
        assert "params.extract_memory" in self.content

    def test_publishdir(self):
        assert "publishDir" in self.content


class TestClassifyModule:
    """Validate the classify module."""

    @pytest.fixture(autouse=True)
    def _load(self):
        self.content = _read(MODULES_DIR / "classify.nf")

    def test_process_defined(self):
        assert "process CLASSIFY_READS" in self.content

    def test_calls_csc_classify(self):
        assert "csc-classify" in self.content

    def test_emits_reports(self):
        assert "emit: reports" in self.content

    def test_paired_end_support(self):
        assert "paired" in self.content.lower()

    def test_configurable_cpus(self):
        assert "params.classify_cpus" in self.content

    def test_configurable_memory(self):
        assert "params.classify_memory" in self.content


class TestAggregateModule:
    """Validate the aggregate module."""

    @pytest.fixture(autouse=True)
    def _load(self):
        self.content = _read(MODULES_DIR / "aggregate.nf")

    def test_process_defined(self):
        assert "process AGGREGATE_REPORTS" in self.content

    def test_calls_csc_aggregate(self):
        assert "csc-aggregate" in self.content

    def test_emits_matrix(self):
        assert "emit: matrix" in self.content

    def test_emits_metadata(self):
        assert "emit: metadata" in self.content

    def test_configurable_cpus(self):
        assert "params.aggregate_cpus" in self.content


class TestDetectModule:
    """Validate the detect module."""

    @pytest.fixture(autouse=True)
    def _load(self):
        self.content = _read(MODULES_DIR / "detect.nf")

    def test_process_defined(self):
        assert "process DETECT_OUTLIERS" in self.content

    def test_calls_csc_detect(self):
        assert "csc-detect" in self.content

    def test_emits_flagged(self):
        assert "emit: flagged" in self.content

    def test_emits_qc_summary(self):
        assert "emit: qc_summary" in self.content

    def test_emits_quarantine(self):
        assert "emit: quarantine" in self.content

    def test_configurable_cpus(self):
        assert "params.detect_cpus" in self.content

    def test_method_param(self):
        assert "params.detect_method" in self.content


class TestSummaryModule:
    """Validate the summary module."""

    @pytest.fixture(autouse=True)
    def _load(self):
        self.content = _read(MODULES_DIR / "summary.nf")

    def test_process_defined(self):
        assert "process PIPELINE_SUMMARY" in self.content

    def test_emits_report(self):
        assert "emit: report" in self.content

    def test_emits_multiqc_data(self):
        assert "emit: multiqc_data" in self.content

    def test_produces_html(self):
        assert "pipeline_report.html" in self.content

    def test_produces_multiqc_yaml(self):
        assert "csc_pipeline_mqc.yaml" in self.content


# ── nextflow.config ──────────────────────────────────────────────────────────


class TestNextflowConfig:
    """Validate nextflow.config."""

    @pytest.fixture(autouse=True)
    def _load(self):
        self.content = _read(NEXTFLOW_DIR / "nextflow.config")

    def test_manifest_name(self):
        assert "cross_species_contamination" in self.content

    def test_main_script(self):
        assert "main.nf" in self.content

    @pytest.mark.parametrize(
        "profile",
        ["standard", "docker", "singularity", "slurm", "test"],
    )
    def test_profile_exists(self, profile: str):
        assert profile in self.content

    def test_docker_profile(self):
        assert "docker.enabled" in self.content

    def test_singularity_profile(self):
        assert "singularity.enabled" in self.content

    def test_slurm_profile(self):
        assert "'slurm'" in self.content

    def test_error_strategy(self):
        assert "errorStrategy" in self.content

    def test_timeline_enabled(self):
        assert "timeline" in self.content

    def test_trace_enabled(self):
        assert "trace" in self.content

    def test_container_defined(self):
        assert "container" in self.content


# ── Documentation ────────────────────────────────────────────────────────────


class TestPipelineDocs:
    """Validate that pipeline documentation exists."""

    def test_pipeline_docs_exist(self):
        docs = Path(__file__).resolve().parent.parent / "docs" / "pipeline.md"
        assert docs.is_file()

    def test_pipeline_docs_content(self):
        docs = Path(__file__).resolve().parent.parent / "docs" / "pipeline.md"
        content = docs.read_text()
        assert "input_csv" in content
        assert "kraken2_db" in content
        assert "MultiQC" in content
        assert "Docker" in content
        assert "Singularity" in content
        assert "SLURM" in content
        assert "-resume" in content
