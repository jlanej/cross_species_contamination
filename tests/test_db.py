"""Tests for the Kraken2 database management module (csc.classify.db).

All tests mock external I/O (network, S3, filesystem copies) so they run
without network access or the AWS CLI.
"""

from __future__ import annotations

import hashlib
import io
import json
import shutil
import tarfile
from pathlib import Path
from unittest import mock

import pytest

from csc.classify.classify import REQUIRED_DB_FILES, validate_database
from csc.classify.db import (
    DEFAULT_CACHE_DIR,
    PRACKENDB_NAME,
    PRACKENDB_URL,
    SUPPORTED_HASH_ALGORITHMS,
    TAXONOMY_FILES,
    _NON_PRACKENDB_WARNING,
    _extract_tarball,
    _human_size,
    _is_s3_uri,
    _is_url,
    clean_cache,
    compute_hash,
    database_info,
    estimate_db_memory,
    fetch_database,
    fetch_prackendb,
    get_cache_dir,
    is_prackendb,
    list_databases,
    validate_taxonomy,
    verify_hash,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def mock_db(tmp_path: Path) -> Path:
    """Create a minimal valid Kraken2 database directory."""
    db_dir = tmp_path / "mock_kraken2_db"
    db_dir.mkdir()
    for fname in REQUIRED_DB_FILES:
        (db_dir / fname).write_bytes(b"\x00" * 64)
    return db_dir


@pytest.fixture()
def cache_dir(tmp_path: Path) -> Path:
    """Return a fresh temporary cache directory."""
    d = tmp_path / "db_cache"
    d.mkdir()
    return d


@pytest.fixture()
def sample_file(tmp_path: Path) -> Path:
    """Create a small file for hashing tests."""
    f = tmp_path / "sample.bin"
    f.write_bytes(b"hello world")
    return f


def _make_db_tarball(tmp_path: Path, db_name: str = "testdb") -> Path:
    """Create a .tar.gz containing a valid mock Kraken2 database."""
    db_dir = tmp_path / db_name
    db_dir.mkdir()
    for fname in REQUIRED_DB_FILES:
        (db_dir / fname).write_bytes(b"\x00" * 64)

    archive = tmp_path / f"{db_name}.tar.gz"
    with tarfile.open(archive, "w:gz") as tf:
        tf.add(db_dir, arcname=db_name)
    return archive


@pytest.fixture()
def prackendb_db(tmp_path: Path) -> Path:
    """Create a mock PrackenDB-like database with taxonomy files."""
    db_dir = tmp_path / "prackendb"
    db_dir.mkdir()
    for fname in REQUIRED_DB_FILES:
        (db_dir / fname).write_bytes(b"\x00" * 64)
    # Add taxonomy files
    tax_dir = db_dir / "taxonomy"
    tax_dir.mkdir()
    (tax_dir / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n")
    (tax_dir / "names.dmp").write_text("1\t|\troot\t|\n")
    return db_dir


def _make_prackendb_tarball(tmp_path: Path, db_name: str = "prackendb") -> Path:
    """Create a .tar.gz containing a mock PrackenDB database."""
    db_dir = tmp_path / db_name
    db_dir.mkdir()
    for fname in REQUIRED_DB_FILES:
        (db_dir / fname).write_bytes(b"\x00" * 64)
    tax_dir = db_dir / "taxonomy"
    tax_dir.mkdir()
    (tax_dir / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n")
    (tax_dir / "names.dmp").write_text("1\t|\troot\t|\n")

    archive = tmp_path / f"{db_name}.tar.gz"
    with tarfile.open(archive, "w:gz") as tf:
        tf.add(db_dir, arcname=db_name)
    return archive


# ---------------------------------------------------------------------------
# Hash helpers
# ---------------------------------------------------------------------------


class TestComputeHash:
    def test_sha256(self, sample_file: Path) -> None:
        expected = hashlib.sha256(b"hello world").hexdigest()
        assert compute_hash(sample_file, "sha256") == expected

    def test_md5(self, sample_file: Path) -> None:
        expected = hashlib.md5(b"hello world").hexdigest()
        assert compute_hash(sample_file, "md5") == expected

    def test_unsupported_algorithm_raises(self, sample_file: Path) -> None:
        with pytest.raises(ValueError, match="Unsupported hash algorithm"):
            compute_hash(sample_file, "sha512")

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            compute_hash(tmp_path / "no_such_file", "sha256")


class TestVerifyHash:
    def test_matching_hash(self, sample_file: Path) -> None:
        expected = hashlib.sha256(b"hello world").hexdigest()
        assert verify_hash(sample_file, expected, "sha256") is True

    def test_mismatched_hash(self, sample_file: Path) -> None:
        assert verify_hash(sample_file, "deadbeef", "sha256") is False

    def test_case_insensitive(self, sample_file: Path) -> None:
        expected = hashlib.sha256(b"hello world").hexdigest().upper()
        assert verify_hash(sample_file, expected, "sha256") is True


# ---------------------------------------------------------------------------
# Source detection
# ---------------------------------------------------------------------------


class TestSourceDetection:
    def test_is_s3_uri(self) -> None:
        assert _is_s3_uri("s3://bucket/key.tar.gz") is True
        assert _is_s3_uri("http://example.com") is False
        assert _is_s3_uri("/local/path") is False

    def test_is_url(self) -> None:
        assert _is_url("https://example.com/db.tar.gz") is True
        assert _is_url("http://example.com/db.tar.gz") is True
        assert _is_url("s3://bucket/key") is False
        assert _is_url("/local/path") is False


# ---------------------------------------------------------------------------
# Cache directory management
# ---------------------------------------------------------------------------


class TestGetCacheDir:
    def test_explicit_path(self, tmp_path: Path) -> None:
        d = tmp_path / "my_cache"
        result = get_cache_dir(d)
        assert result == d.resolve()
        assert result.is_dir()

    def test_env_var(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        d = tmp_path / "env_cache"
        monkeypatch.setenv("CSC_DB_CACHE", str(d))
        result = get_cache_dir()
        assert result == d.resolve()
        assert result.is_dir()

    def test_default_when_no_override(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("CSC_DB_CACHE", raising=False)
        # Don't actually create in home dir; just check the path logic
        with mock.patch("csc.classify.db.DEFAULT_CACHE_DIR", Path("/fake/.csc/db")):
            with mock.patch.object(Path, "mkdir"):
                result = get_cache_dir()
                assert str(result).endswith(".csc/db")


class TestListDatabases:
    def test_empty_cache(self, cache_dir: Path) -> None:
        dbs = list_databases(cache_dir)
        assert dbs == []

    def test_lists_valid_db(self, cache_dir: Path) -> None:
        db = cache_dir / "mydb"
        db.mkdir()
        for f in REQUIRED_DB_FILES:
            (db / f).write_bytes(b"\x00" * 64)
        dbs = list_databases(cache_dir)
        assert len(dbs) == 1
        assert dbs[0]["name"] == "mydb"
        assert dbs[0]["valid"] is True
        assert dbs[0]["size_bytes"] == 64 * 3

    def test_lists_invalid_db(self, cache_dir: Path) -> None:
        db = cache_dir / "bad_db"
        db.mkdir()
        (db / "hash.k2d").write_bytes(b"\x00")
        dbs = list_databases(cache_dir)
        assert len(dbs) == 1
        assert dbs[0]["valid"] is False


class TestCleanCache:
    def test_clean_specific(self, cache_dir: Path) -> None:
        (cache_dir / "db1").mkdir()
        (cache_dir / "db2").mkdir()
        removed = clean_cache(cache_dir, name="db1")
        assert removed == ["db1"]
        assert not (cache_dir / "db1").exists()
        assert (cache_dir / "db2").exists()

    def test_clean_all(self, cache_dir: Path) -> None:
        (cache_dir / "db1").mkdir()
        (cache_dir / "db2").mkdir()
        removed = clean_cache(cache_dir)
        assert sorted(removed) == ["db1", "db2"]
        assert list(cache_dir.iterdir()) == []

    def test_clean_nonexistent(self, cache_dir: Path) -> None:
        removed = clean_cache(cache_dir, name="nope")
        assert removed == []


# ---------------------------------------------------------------------------
# Database info
# ---------------------------------------------------------------------------


class TestDatabaseInfo:
    def test_valid_db_info(self, mock_db: Path) -> None:
        info = database_info(mock_db)
        assert info["valid"] is True
        assert info["total_size_bytes"] == 64 * 3
        assert set(info["files"].keys()) == set(REQUIRED_DB_FILES)
        assert len(info["sha256"]) == 3

    def test_invalid_path(self, tmp_path: Path) -> None:
        info = database_info(tmp_path / "nonexistent")
        assert info["valid"] is False
        assert info["total_size_bytes"] == 0


# ---------------------------------------------------------------------------
# Tarball extraction
# ---------------------------------------------------------------------------


class TestExtractTarball:
    def test_extract_single_dir(self, tmp_path: Path) -> None:
        archive = _make_db_tarball(tmp_path, "mydb")
        dest = tmp_path / "extract_dest"
        result = _extract_tarball(archive, dest)
        # Should return the inner directory
        assert result.name == "mydb"
        for f in REQUIRED_DB_FILES:
            assert (result / f).exists()

    def test_unsafe_path_raises(self, tmp_path: Path) -> None:
        # Create a tarball with an unsafe path
        archive = tmp_path / "evil.tar.gz"
        with tarfile.open(archive, "w:gz") as tf:
            data = b"evil content"
            info = tarfile.TarInfo(name="../../../etc/passwd")
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))

        with pytest.raises(ValueError, match="unsafe path"):
            _extract_tarball(archive, tmp_path / "dest")


# ---------------------------------------------------------------------------
# Fetch database
# ---------------------------------------------------------------------------


class TestFetchDatabase:
    def test_fetch_local_directory(self, mock_db: Path, cache_dir: Path) -> None:
        """Fetching a local valid DB directory should return it directly."""
        result = fetch_database(str(mock_db), cache_dir=cache_dir)
        assert result == mock_db.resolve()

    def test_fetch_local_tarball(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        archive = _make_db_tarball(tmp_path, "localdb")

        result = fetch_database(str(archive), name="localdb", cache_dir=cache)
        assert result.name == "localdb"
        assert (result / "hash.k2d").exists()

    def test_fetch_local_tarball_with_hash(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        archive = _make_db_tarball(tmp_path, "hashdb")
        sha = compute_hash(archive, "sha256")

        result = fetch_database(
            str(archive),
            name="hashdb",
            cache_dir=cache,
            expected_hash=sha,
            hash_algorithm="sha256",
        )
        assert (result / "hash.k2d").exists()

    def test_fetch_local_tarball_bad_hash(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        archive = _make_db_tarball(tmp_path, "baddb")

        with pytest.raises(ValueError, match="Hash verification failed"):
            fetch_database(
                str(archive),
                name="baddb",
                cache_dir=cache,
                expected_hash="deadbeef",
                hash_algorithm="sha256",
            )

    def test_fetch_missing_local_raises(self, cache_dir: Path) -> None:
        with pytest.raises(FileNotFoundError, match="Local source not found"):
            fetch_database("/no/such/path", cache_dir=cache_dir)

    @mock.patch("csc.classify.db._download_http")
    def test_fetch_http(
        self, mock_dl: mock.Mock, tmp_path: Path
    ) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        # Create the archive that _download_http would produce
        archive = _make_db_tarball(tmp_path, "httpdb")

        def fake_download(url: str, dest: Path) -> Path:
            shutil.copy2(archive, dest)
            return dest

        mock_dl.side_effect = fake_download

        result = fetch_database(
            "https://example.com/httpdb.tar.gz",
            name="httpdb",
            cache_dir=cache,
        )
        assert result.name == "httpdb"
        assert (result / "hash.k2d").exists()
        mock_dl.assert_called_once()

    @mock.patch("csc.classify.db._copy_s3")
    def test_fetch_s3(
        self, mock_s3: mock.Mock, tmp_path: Path
    ) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        archive = _make_db_tarball(tmp_path, "s3db")

        def fake_s3(uri: str, dest: Path) -> Path:
            shutil.copy2(archive, dest)
            return dest

        mock_s3.side_effect = fake_s3

        result = fetch_database(
            "s3://my-bucket/s3db.tar.gz",
            name="s3db",
            cache_dir=cache,
        )
        assert result.name == "s3db"
        assert (result / "hash.k2d").exists()
        mock_s3.assert_called_once()


# ---------------------------------------------------------------------------
# Human-size helper
# ---------------------------------------------------------------------------


class TestHumanSize:
    def test_bytes(self) -> None:
        assert "B" in _human_size(500)

    def test_kib(self) -> None:
        assert "KiB" in _human_size(2048)

    def test_mib(self) -> None:
        assert "MiB" in _human_size(5 * 1024 * 1024)


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestDBCLI:
    def test_version(self) -> None:
        from csc.classify.db_cli import main

        with pytest.raises(SystemExit) as exc_info:
            main(["--version"])
        assert exc_info.value.code == 0

    def test_no_command_returns_error(self) -> None:
        from csc.classify.db_cli import main

        rc = main([])
        assert rc == 1

    def test_list_empty(self, cache_dir: Path) -> None:
        from csc.classify.db_cli import main

        rc = main(["--cache-dir", str(cache_dir), "list"])
        assert rc == 0

    def test_verify_valid_db(self, mock_db: Path) -> None:
        from csc.classify.db_cli import main

        rc = main(["verify", str(mock_db)])
        assert rc == 0

    def test_verify_invalid_db(self, tmp_path: Path) -> None:
        from csc.classify.db_cli import main

        rc = main(["verify", str(tmp_path / "nope")])
        assert rc == 1

    def test_info_valid_db(self, mock_db: Path) -> None:
        from csc.classify.db_cli import main

        rc = main(["info", str(mock_db)])
        assert rc == 0

    def test_info_json(self, mock_db: Path, capsys: pytest.CaptureFixture[str]) -> None:
        from csc.classify.db_cli import main

        rc = main(["info", str(mock_db), "--json"])
        assert rc == 0
        output = capsys.readouterr().out
        data = json.loads(output)
        assert data["valid"] is True
        assert "taxonomy" in data
        assert "prackendb_compatible" in data

    def test_clean_empty(self, cache_dir: Path) -> None:
        from csc.classify.db_cli import main

        rc = main(["--cache-dir", str(cache_dir), "clean"])
        assert rc == 0

    def test_clean_named(self, cache_dir: Path) -> None:
        from csc.classify.db_cli import main

        (cache_dir / "mydb").mkdir()
        rc = main(["--cache-dir", str(cache_dir), "clean", "--name", "mydb"])
        assert rc == 0
        assert not (cache_dir / "mydb").exists()

    def test_fetch_local_dir(self, mock_db: Path, cache_dir: Path) -> None:
        from csc.classify.db_cli import main

        rc = main([
            "--cache-dir", str(cache_dir),
            "fetch", str(mock_db),
        ])
        assert rc == 0

    def test_fetch_failure(self, cache_dir: Path) -> None:
        from csc.classify.db_cli import main

        rc = main([
            "--cache-dir", str(cache_dir),
            "fetch", "/no/such/path",
        ])
        assert rc == 1

    def test_verify_shows_prackendb_status(
        self, prackendb_db: Path, capsys: pytest.CaptureFixture[str],
    ) -> None:
        from csc.classify.db_cli import main

        rc = main(["verify", str(prackendb_db)])
        assert rc == 0
        output = capsys.readouterr().out
        assert "PrackenDB-compatible: yes" in output
        assert "taxonomy/nodes.dmp [OK]" in output
        assert "taxonomy/names.dmp [OK]" in output

    def test_verify_non_prackendb_warns(
        self, mock_db: Path, capsys: pytest.CaptureFixture[str],
    ) -> None:
        from csc.classify.db_cli import main

        rc = main(["verify", str(mock_db)])
        assert rc == 0
        output = capsys.readouterr().out
        assert "PrackenDB-compatible: no" in output

    def test_info_shows_prackendb_status(
        self, prackendb_db: Path, capsys: pytest.CaptureFixture[str],
    ) -> None:
        from csc.classify.db_cli import main

        rc = main(["info", str(prackendb_db)])
        assert rc == 0
        output = capsys.readouterr().out
        assert "PrackenDB-compatible: yes" in output


# ---------------------------------------------------------------------------
# Taxonomy validation
# ---------------------------------------------------------------------------


class TestValidateTaxonomy:
    def test_all_taxonomy_present(self, prackendb_db: Path) -> None:
        result = validate_taxonomy(prackendb_db)
        assert result["taxonomy/nodes.dmp"] is True
        assert result["taxonomy/names.dmp"] is True

    def test_no_taxonomy(self, mock_db: Path) -> None:
        result = validate_taxonomy(mock_db)
        assert result["taxonomy/nodes.dmp"] is False
        assert result["taxonomy/names.dmp"] is False

    def test_partial_taxonomy(self, tmp_path: Path) -> None:
        db_dir = tmp_path / "partial_db"
        db_dir.mkdir()
        for fname in REQUIRED_DB_FILES:
            (db_dir / fname).write_bytes(b"\x00" * 64)
        tax_dir = db_dir / "taxonomy"
        tax_dir.mkdir()
        (tax_dir / "nodes.dmp").write_text("1\t|\t1\t|\n")
        # names.dmp missing
        result = validate_taxonomy(db_dir)
        assert result["taxonomy/nodes.dmp"] is True
        assert result["taxonomy/names.dmp"] is False


# ---------------------------------------------------------------------------
# PrackenDB detection
# ---------------------------------------------------------------------------


class TestIsPrackenDB:
    def test_prackendb_detected(self, prackendb_db: Path) -> None:
        assert is_prackendb(prackendb_db) is True

    def test_non_prackendb_without_taxonomy(self, mock_db: Path) -> None:
        assert is_prackendb(mock_db) is False

    def test_invalid_db_returns_false(self, tmp_path: Path) -> None:
        assert is_prackendb(tmp_path / "nonexistent") is False

    def test_partial_taxonomy_not_prackendb(self, tmp_path: Path) -> None:
        db_dir = tmp_path / "partial"
        db_dir.mkdir()
        for fname in REQUIRED_DB_FILES:
            (db_dir / fname).write_bytes(b"\x00" * 64)
        (db_dir / "taxonomy").mkdir()
        (db_dir / "taxonomy" / "nodes.dmp").write_text("1\t|\t1\t|\n")
        # names.dmp missing → not PrackenDB
        assert is_prackendb(db_dir) is False


# ---------------------------------------------------------------------------
# PrackenDB constants
# ---------------------------------------------------------------------------


class TestPrackenDBConstants:
    def test_prackendb_url_is_https(self) -> None:
        assert PRACKENDB_URL.startswith("https://")
        assert PRACKENDB_URL.startswith("https://genome-idx.s3.amazonaws.com/")

    def test_prackendb_name(self) -> None:
        assert PRACKENDB_NAME == "prackendb"

    def test_taxonomy_files(self) -> None:
        assert "taxonomy/nodes.dmp" in TAXONOMY_FILES
        assert "taxonomy/names.dmp" in TAXONOMY_FILES

    def test_non_prackendb_warning_content(self) -> None:
        assert "PrackenDB" in _NON_PRACKENDB_WARNING
        assert "one genome per species" in _NON_PRACKENDB_WARNING


# ---------------------------------------------------------------------------
# fetch_prackendb
# ---------------------------------------------------------------------------


class TestFetchPrackenDB:
    @mock.patch("csc.classify.db._download_http")
    def test_fetch_prackendb_with_taxonomy(
        self, mock_dl: mock.Mock, tmp_path: Path,
    ) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        archive = _make_prackendb_tarball(tmp_path, "prackendb")

        def fake_download(url: str, dest: Path) -> Path:
            shutil.copy2(archive, dest)
            return dest

        mock_dl.side_effect = fake_download

        result = fetch_prackendb(cache_dir=cache)
        assert result.name == PRACKENDB_NAME
        assert (result / "hash.k2d").exists()
        assert (result / "taxonomy" / "nodes.dmp").exists()
        assert (result / "taxonomy" / "names.dmp").exists()
        assert is_prackendb(result) is True

    @mock.patch("csc.classify.db._download_http")
    def test_fetch_prackendb_warns_on_missing_taxonomy(
        self, mock_dl: mock.Mock, tmp_path: Path,
    ) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        # Use regular tarball without taxonomy files
        archive = _make_db_tarball(tmp_path, "prackendb")

        def fake_download(url: str, dest: Path) -> Path:
            shutil.copy2(archive, dest)
            return dest

        mock_dl.side_effect = fake_download

        with mock.patch("csc.classify.db.logger") as mock_logger:
            result = fetch_prackendb(cache_dir=cache)
            # Should have warned about missing taxonomy files
            warning_calls = [
                str(c) for c in mock_logger.warning.call_args_list
                if "taxonomy file not found" in str(c)
            ]
            assert len(warning_calls) == 2  # nodes.dmp and names.dmp
            assert any("nodes.dmp" in w for w in warning_calls)
            assert any("names.dmp" in w for w in warning_calls)

    @mock.patch("csc.classify.db._download_http")
    def test_cli_fetch_prackendb(
        self, mock_dl: mock.Mock, tmp_path: Path,
    ) -> None:
        from csc.classify.db_cli import main

        cache = tmp_path / "cache"
        cache.mkdir()
        archive = _make_prackendb_tarball(tmp_path, "prackendb")

        def fake_download(url: str, dest: Path) -> Path:
            shutil.copy2(archive, dest)
            return dest

        mock_dl.side_effect = fake_download

        rc = main(["--cache-dir", str(cache), "fetch", "prackendb"])
        assert rc == 0


# ---------------------------------------------------------------------------
# Config: recommended_db
# ---------------------------------------------------------------------------


class TestConfigRecommendedDB:
    def test_default_config_has_recommended_db(self) -> None:
        from csc.config import load_config

        cfg = load_config()
        assert cfg["classify"]["recommended_db"] == "prackendb"


# ---------------------------------------------------------------------------
# estimate_db_memory
# ---------------------------------------------------------------------------


class TestEstimateDbMemory:
    def test_returns_expected_keys(self, mock_db: Path) -> None:
        est = estimate_db_memory(mock_db)
        assert "db_path" in est
        assert "hash_k2d_bytes" in est
        assert "total_db_bytes" in est
        assert "estimated_ram_bytes" in est
        assert "available_ram_bytes" in est
        assert "recommend_memory_mapping" in est
        assert "human" in est

    def test_hash_bytes_equals_hash_file_size(self, mock_db: Path) -> None:
        est = estimate_db_memory(mock_db)
        expected = (mock_db / "hash.k2d").stat().st_size
        assert est["hash_k2d_bytes"] == expected

    def test_estimated_ram_equals_hash_bytes(self, mock_db: Path) -> None:
        est = estimate_db_memory(mock_db)
        assert est["estimated_ram_bytes"] == est["hash_k2d_bytes"]

    def test_total_db_bytes_is_sum_of_files(self, mock_db: Path) -> None:
        est = estimate_db_memory(mock_db)
        expected_total = sum(
            f.stat().st_size for f in mock_db.rglob("*") if f.is_file()
        )
        assert est["total_db_bytes"] == expected_total

    def test_human_readable_strings_present(self, mock_db: Path) -> None:
        est = estimate_db_memory(mock_db)
        h = est["human"]
        assert isinstance(h["hash_k2d"], str)
        assert isinstance(h["total_db"], str)
        assert isinstance(h["estimated_ram"], str)
        assert isinstance(h["available_ram"], str)

    def test_recommend_memory_mapping_when_db_exceeds_ram(
        self, mock_db: Path
    ) -> None:
        # Force available RAM to be smaller than the DB
        with mock.patch("csc.classify.db._available_ram_bytes", return_value=1):
            est = estimate_db_memory(mock_db)
        assert est["recommend_memory_mapping"] is True

    def test_no_recommend_memory_mapping_when_db_fits_in_ram(
        self, mock_db: Path
    ) -> None:
        # Force available RAM to be much larger than the DB
        with mock.patch(
            "csc.classify.db._available_ram_bytes", return_value=100 * 1024 ** 3
        ):
            est = estimate_db_memory(mock_db)
        assert est["recommend_memory_mapping"] is False

    def test_recommend_memory_mapping_when_available_ram_unknown_large_db(
        self, tmp_path: Path
    ) -> None:
        # Build a DB with a small hash.k2d (< 1 GiB) — should NOT recommend mm
        db_dir = tmp_path / "bigdb"
        db_dir.mkdir()
        for fname in REQUIRED_DB_FILES:
            (db_dir / fname).write_bytes(b"\x00" * 64)
        # Without patching: small test file → no recommendation even with unknown RAM
        with mock.patch("csc.classify.db._available_ram_bytes", return_value=None):
            est = estimate_db_memory(db_dir)
        assert est["recommend_memory_mapping"] is False  # tiny test file

    def test_invalid_db_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            estimate_db_memory(tmp_path / "nonexistent")

    def test_missing_db_files_raises(self, tmp_path: Path) -> None:
        d = tmp_path / "emptydb"
        d.mkdir()
        with pytest.raises(ValueError):
            estimate_db_memory(d)

    def test_available_ram_bytes_field(self, mock_db: Path) -> None:
        # With a known mock value
        with mock.patch("csc.classify.db._available_ram_bytes", return_value=8 * 1024 ** 3):
            est = estimate_db_memory(mock_db)
        assert est["available_ram_bytes"] == 8 * 1024 ** 3
        assert est["human"]["available_ram"] != "unknown"

    def test_available_ram_bytes_unknown(self, mock_db: Path) -> None:
        with mock.patch("csc.classify.db._available_ram_bytes", return_value=None):
            est = estimate_db_memory(mock_db)
        assert est["available_ram_bytes"] is None
        assert est["human"]["available_ram"] == "unknown"


# ---------------------------------------------------------------------------
# CLI: estimate-memory subcommand
# ---------------------------------------------------------------------------


class TestDBCLIEstimateMemory:
    def test_estimate_memory_valid_db(
        self, mock_db: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        from csc.classify.db_cli import main

        with mock.patch(
            "csc.classify.db._available_ram_bytes", return_value=100 * 1024 ** 3
        ):
            rc = main(["estimate-memory", str(mock_db)])
        assert rc == 0
        output = capsys.readouterr().out
        assert "hash.k2d" in output
        assert "Estimated RAM" in output

    def test_estimate_memory_json(
        self, mock_db: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        from csc.classify.db_cli import main

        with mock.patch(
            "csc.classify.db._available_ram_bytes", return_value=100 * 1024 ** 3
        ):
            rc = main(["estimate-memory", str(mock_db), "--json"])
        assert rc == 0
        output = capsys.readouterr().out
        data = json.loads(output)
        assert "hash_k2d_bytes" in data
        assert "estimated_ram_bytes" in data
        assert "recommend_memory_mapping" in data

    def test_estimate_memory_invalid_db(self, tmp_path: Path) -> None:
        from csc.classify.db_cli import main

        rc = main(["estimate-memory", str(tmp_path / "nonexistent")])
        assert rc == 1

    def test_estimate_memory_recommends_memory_mapping(
        self, mock_db: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        from csc.classify.db_cli import main

        with mock.patch("csc.classify.db._available_ram_bytes", return_value=1):
            rc = main(["estimate-memory", str(mock_db)])
        assert rc == 0
        output = capsys.readouterr().out
        assert "--memory-mapping" in output
