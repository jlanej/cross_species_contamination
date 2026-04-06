"""Tests for the central configuration loader (csc.config)."""

from __future__ import annotations

from pathlib import Path

import pytest

from csc.config import load_config


class TestLoadDefaults:
    """Loading without any user override returns built-in defaults."""

    def test_returns_dict(self) -> None:
        cfg = load_config()
        assert isinstance(cfg, dict)

    def test_has_extract_section(self) -> None:
        cfg = load_config()
        assert "extract" in cfg

    def test_has_classify_section(self) -> None:
        cfg = load_config()
        assert "classify" in cfg

    def test_has_aggregate_section(self) -> None:
        cfg = load_config()
        assert "aggregate" in cfg

    def test_has_detect_section(self) -> None:
        cfg = load_config()
        assert "detect" in cfg

    def test_has_logging_section(self) -> None:
        cfg = load_config()
        assert "logging" in cfg

    def test_default_extract_threads(self) -> None:
        cfg = load_config()
        assert cfg["extract"]["threads"] == 1

    def test_default_mapq_threshold_is_none(self) -> None:
        cfg = load_config()
        assert cfg["extract"]["mapq_threshold"] is None

    def test_default_log_level(self) -> None:
        cfg = load_config()
        assert cfg["logging"]["level"] == "INFO"


class TestLoadUserOverride:
    """Loading with a user YAML merges values into defaults."""

    def test_override_single_key(self, tmp_path: Path) -> None:
        user_cfg = tmp_path / "user.yaml"
        user_cfg.write_text("extract:\n  threads: 8\n")
        cfg = load_config(user_cfg)
        assert cfg["extract"]["threads"] == 8
        # Other defaults preserved
        assert cfg["extract"]["mapq_threshold"] is None

    def test_add_new_section(self, tmp_path: Path) -> None:
        user_cfg = tmp_path / "user.yaml"
        user_cfg.write_text("custom:\n  key: value\n")
        cfg = load_config(user_cfg)
        assert cfg["custom"]["key"] == "value"
        # Defaults still present
        assert "extract" in cfg

    def test_deep_merge(self, tmp_path: Path) -> None:
        user_cfg = tmp_path / "user.yaml"
        user_cfg.write_text("classify:\n  confidence: 0.5\n")
        cfg = load_config(user_cfg)
        assert cfg["classify"]["confidence"] == 0.5
        # Other classify defaults preserved
        assert cfg["classify"]["tool"] == "kraken2"


class TestLoadFromEnv:
    """CSC_CONFIG environment variable points to user config."""

    def test_env_override(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        user_cfg = tmp_path / "env.yaml"
        user_cfg.write_text("extract:\n  threads: 16\n")
        monkeypatch.setenv("CSC_CONFIG", str(user_cfg))
        cfg = load_config()
        assert cfg["extract"]["threads"] == 16

    def test_explicit_path_beats_env(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        env_cfg = tmp_path / "env.yaml"
        env_cfg.write_text("extract:\n  threads: 16\n")
        explicit_cfg = tmp_path / "explicit.yaml"
        explicit_cfg.write_text("extract:\n  threads: 32\n")
        monkeypatch.setenv("CSC_CONFIG", str(env_cfg))
        cfg = load_config(explicit_cfg)
        assert cfg["extract"]["threads"] == 32


class TestLoadErrors:
    """Config loader raises on bad input."""

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="Config file not found"):
            load_config(tmp_path / "nonexistent.yaml")
