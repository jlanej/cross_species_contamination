"""Shared pytest fixtures for extraction pipeline tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.generate_test_data import generate_test_data


@pytest.fixture(scope="session")
def test_data(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """Generate synthetic BAM + CRAM + reference once per test session."""
    outdir = tmp_path_factory.mktemp("test_data")
    return generate_test_data(outdir)


@pytest.fixture(scope="session")
def test_bam(test_data: dict[str, Path]) -> Path:
    return test_data["bam"]


@pytest.fixture(scope="session")
def test_cram(test_data: dict[str, Path]) -> Path:
    return test_data["cram"]


@pytest.fixture(scope="session")
def test_reference(test_data: dict[str, Path]) -> Path:
    return test_data["reference"]
