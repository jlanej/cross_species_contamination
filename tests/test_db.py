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
    SUPPORTED_HASH_ALGORITHMS,
    _extract_tarball,
    _human_size,
    _is_s3_uri,
    _is_url,
    clean_cache,
    compute_hash,
    database_info,
    fetch_database,
    get_cache_dir,
    list_databases,
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
