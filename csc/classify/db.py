"""Kraken2 database management: download, verify, and cache.

Provides helpers to fetch Kraken2 databases from local paths, HTTP(S) URLs,
or S3 URIs, validate them via MD5/SHA-256 checksums, and manage a local
cache directory.

PrackenDB is the recommended Kraken2 database for the CSC pipeline.  It
contains one NCBI reference genome per species (bacteria, archaea, fungi,
protists, viruses, human, and UniVec Core), which avoids inflated LCA
assignments from redundant genomes and enables robust species-level
detection.  See :func:`fetch_prackendb` for a convenience wrapper that
downloads and validates PrackenDB.

AI assistance acknowledgment: This module was developed with AI assistance.
Best practices in the bioinformatics field should always take precedence over
specific implementation details.
"""

from __future__ import annotations

import hashlib
import logging
import os
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
from pathlib import Path
from typing import Any

from csc.classify.classify import REQUIRED_DB_FILES, validate_database

logger = logging.getLogger(__name__)

#: Default cache directory for downloaded databases.
DEFAULT_CACHE_DIR = Path.home() / ".csc" / "db"

#: Supported hash algorithms for verification.
SUPPORTED_HASH_ALGORITHMS = ("md5", "sha256")

#: Chunk size for streaming downloads / hash computation (8 MiB).
_CHUNK_SIZE = 8 * 1024 * 1024

# ---------------------------------------------------------------------------
# PrackenDB constants
# ---------------------------------------------------------------------------

#: Default download URL for the PrackenDB Kraken2 database.
PRACKENDB_URL = (
    "https://genome-idx.s3.amazonaws.com/kraken/k2_NCBI_reference_20251007.tar.gz"
)

#: Default name for PrackenDB in the cache directory.
PRACKENDB_NAME = "prackendb"

#: Taxonomy files expected in a PrackenDB (or compatible) database.
TAXONOMY_FILES = ("taxonomy/nodes.dmp", "taxonomy/names.dmp")

#: Warning emitted when the user selects a non-PrackenDB database.
_NON_PRACKENDB_WARNING = (
    "The selected database does not appear to be PrackenDB. "
    "PrackenDB (one genome per species) is recommended for the CSC pipeline "
    "because it provides unambiguous per-species k-mer counts and avoids "
    "inflated LCA assignments from redundant genomes. Using other databases "
    "may result in species-level ambiguity or loss of granularity. "
    "See: https://github.com/jlanej/kmer_denovo_filter/blob/main/"
    "docs/kraken2_bacterial_detection.md#the-prackendb-reference-database"
)


# ---------------------------------------------------------------------------
# Hash helpers
# ---------------------------------------------------------------------------


def compute_hash(path: Path, algorithm: str = "sha256") -> str:
    """Compute the hex-digest of a file using *algorithm*.

    Parameters
    ----------
    path:
        Path to the file to hash.
    algorithm:
        Hash algorithm name (``"md5"`` or ``"sha256"``).

    Returns
    -------
    str
        Hex-encoded hash digest.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If *algorithm* is not supported.
    """
    if algorithm not in SUPPORTED_HASH_ALGORITHMS:
        raise ValueError(
            f"Unsupported hash algorithm: {algorithm!r}. "
            f"Supported: {', '.join(SUPPORTED_HASH_ALGORITHMS)}"
        )
    if not path.is_file():
        raise FileNotFoundError(f"File not found: {path}")

    h = hashlib.new(algorithm)
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(_CHUNK_SIZE)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def verify_hash(
    path: Path, expected: str, algorithm: str = "sha256"
) -> bool:
    """Return ``True`` when the hash of *path* matches *expected*.

    Parameters
    ----------
    path:
        Path to the file to verify.
    expected:
        Expected hex-encoded hash digest.
    algorithm:
        Hash algorithm to use.

    Returns
    -------
    bool
    """
    actual = compute_hash(path, algorithm)
    match = actual == expected.lower().strip()
    if not match:
        logger.warning(
            "Hash mismatch for %s (%s): expected %s, got %s",
            path,
            algorithm,
            expected,
            actual,
        )
    return match


# ---------------------------------------------------------------------------
# Source detection
# ---------------------------------------------------------------------------

def _is_s3_uri(source: str) -> bool:
    return source.startswith("s3://")


def _is_url(source: str) -> bool:
    return source.startswith("http://") or source.startswith("https://")


# ---------------------------------------------------------------------------
# Download / copy helpers
# ---------------------------------------------------------------------------


def _download_http(url: str, dest: Path) -> Path:
    """Download a file from *url* to *dest* and return *dest*."""
    logger.info("Downloading %s -> %s", url, dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "csc-db/0.2"})
    with urllib.request.urlopen(req) as resp, open(dest, "wb") as out:  # noqa: S310
        total = 0
        while True:
            chunk = resp.read(_CHUNK_SIZE)
            if not chunk:
                break
            out.write(chunk)
            total += len(chunk)
    logger.info("Downloaded %d bytes to %s", total, dest)
    return dest


def _copy_s3(uri: str, dest: Path) -> Path:
    """Copy a file from an S3 URI using the ``aws`` CLI."""
    aws = shutil.which("aws")
    if aws is None:
        raise FileNotFoundError(
            "AWS CLI (aws) not found on PATH. "
            "Install the AWS CLI to download from S3."
        )
    dest.parent.mkdir(parents=True, exist_ok=True)
    logger.info("S3 copy %s -> %s", uri, dest)
    proc = subprocess.run(
        [aws, "s3", "cp", uri, str(dest)],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"aws s3 cp failed: {proc.stderr}")
    return dest


def _copy_local(source: Path, dest: Path) -> Path:
    """Copy a local file or directory to *dest*."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(source, dest)
    else:
        shutil.copy2(source, dest)
    return dest


# ---------------------------------------------------------------------------
# Tarball extraction
# ---------------------------------------------------------------------------


def _extract_tarball(archive: Path, dest_dir: Path) -> Path:
    """Extract a tarball into *dest_dir* and return the extracted DB path.

    If the tarball contains a single top-level directory, that directory is
    returned.  Otherwise *dest_dir* itself is returned.
    """
    logger.info("Extracting %s -> %s", archive, dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:*") as tf:
        # Security: check for path traversal
        for member in tf.getmembers():
            member_path = os.path.normpath(member.name)
            if member_path.startswith("..") or os.path.isabs(member_path):
                raise ValueError(
                    f"Tarball contains unsafe path: {member.name}"
                )
        tf.extractall(dest_dir)  # noqa: S202

    # If a single directory was extracted, return it
    entries = list(dest_dir.iterdir())
    if len(entries) == 1 and entries[0].is_dir():
        return entries[0]
    return dest_dir


# ---------------------------------------------------------------------------
# Cache directory management
# ---------------------------------------------------------------------------


def get_cache_dir(cache_dir: str | Path | None = None) -> Path:
    """Return the resolved cache directory path, creating it if needed.

    Priority: explicit argument > ``CSC_DB_CACHE`` env var > default.
    """
    if cache_dir is not None:
        d = Path(cache_dir)
    else:
        env = os.environ.get("CSC_DB_CACHE")
        d = Path(env) if env else DEFAULT_CACHE_DIR
    d = d.resolve()
    d.mkdir(parents=True, exist_ok=True)
    return d


def list_databases(cache_dir: str | Path | None = None) -> list[dict[str, Any]]:
    """List Kraken2 databases present in the cache directory.

    Returns
    -------
    list[dict]
        Each dict has keys ``name``, ``path``, ``size_bytes``, ``valid``.
    """
    d = get_cache_dir(cache_dir)
    dbs: list[dict[str, Any]] = []
    if not d.is_dir():
        return dbs
    for entry in sorted(d.iterdir()):
        if not entry.is_dir():
            continue
        size = sum(f.stat().st_size for f in entry.rglob("*") if f.is_file())
        try:
            validate_database(entry)
            valid = True
        except (FileNotFoundError, ValueError):
            valid = False
        dbs.append(
            {"name": entry.name, "path": entry, "size_bytes": size, "valid": valid}
        )
    return dbs


def clean_cache(
    cache_dir: str | Path | None = None, name: str | None = None
) -> list[str]:
    """Remove databases from the cache.

    Parameters
    ----------
    cache_dir:
        Cache directory path.
    name:
        If given, remove only the named database.  Otherwise remove all.

    Returns
    -------
    list[str]
        Names of removed databases.
    """
    d = get_cache_dir(cache_dir)
    removed: list[str] = []
    if name is not None:
        target = d / name
        if target.is_dir():
            shutil.rmtree(target)
            removed.append(name)
            logger.info("Removed database: %s", target)
        else:
            logger.warning("Database not found in cache: %s", name)
    else:
        for entry in d.iterdir():
            if entry.is_dir():
                shutil.rmtree(entry)
                removed.append(entry.name)
        if removed:
            logger.info("Removed %d database(s) from cache", len(removed))
    return removed


# ---------------------------------------------------------------------------
# Database info
# ---------------------------------------------------------------------------


def database_info(db_path: str | Path) -> dict[str, Any]:
    """Return metadata about a Kraken2 database.

    Parameters
    ----------
    db_path:
        Path to the Kraken2 database directory.

    Returns
    -------
    dict
        Keys: ``path``, ``valid``, ``files``, ``total_size_bytes``,
        ``sha256`` (dict mapping file name to hex digest).
    """
    db = Path(db_path).resolve()
    try:
        validate_database(db)
        valid = True
    except (FileNotFoundError, ValueError):
        valid = False

    files: dict[str, int] = {}
    hashes: dict[str, str] = {}
    total_size = 0
    if db.is_dir():
        for fp in sorted(db.iterdir()):
            if fp.is_file():
                sz = fp.stat().st_size
                files[fp.name] = sz
                total_size += sz
                hashes[fp.name] = compute_hash(fp, "sha256")

    return {
        "path": str(db),
        "valid": valid,
        "files": files,
        "total_size_bytes": total_size,
        "sha256": hashes,
    }


# ---------------------------------------------------------------------------
# Fetch / download entry-point
# ---------------------------------------------------------------------------


def fetch_database(
    source: str,
    *,
    name: str | None = None,
    cache_dir: str | Path | None = None,
    expected_hash: str | None = None,
    hash_algorithm: str = "sha256",
) -> Path:
    """Fetch a Kraken2 database from *source* into the cache.

    Parameters
    ----------
    source:
        One of:

        * A local directory path (used directly or copied to cache).
        * An HTTP(S) URL pointing to a ``.tar.gz`` archive.
        * An S3 URI (``s3://bucket/key``) pointing to a ``.tar.gz`` archive.
    name:
        Name for the database inside the cache directory.  Defaults to
        the basename of *source*.
    cache_dir:
        Override cache directory.
    expected_hash:
        If provided, the downloaded archive is verified against this
        digest before extraction.
    hash_algorithm:
        Hash algorithm for verification (``"md5"`` or ``"sha256"``).

    Returns
    -------
    Path
        Path to the ready-to-use Kraken2 database directory.

    Raises
    ------
    ValueError
        If hash verification fails.
    FileNotFoundError
        If the source is not found or required tools are missing.
    """
    cache = get_cache_dir(cache_dir)
    source_str = str(source).rstrip("/")

    # Determine DB name
    if name is None:
        name = Path(source_str.split("?")[0].split("#")[0]).stem
        # Strip common archive suffixes
        for suffix in (".tar.gz", ".tar", ".tgz"):
            if name.endswith(suffix.replace(".", "")):
                name = name[: -len(suffix.replace(".", ""))]
                break
            if source_str.endswith(suffix):
                name = Path(source_str.split("?")[0]).name
                for s in (".tar.gz", ".tar", ".tgz"):
                    if name.endswith(s):
                        name = name[: -len(s)]
                        break
                break

    db_dest = cache / name

    # --- Local directory: validate in-place or copy ---
    local_path = Path(source_str)
    if not _is_s3_uri(source_str) and not _is_url(source_str):
        if local_path.is_dir():
            validated = validate_database(local_path)
            logger.info(
                "Using local database: %s (%s)",
                validated,
                _human_size(_dir_size(validated)),
            )
            return validated
        if local_path.is_file():
            # Treat as tarball
            archive = local_path
        else:
            raise FileNotFoundError(f"Local source not found: {source_str}")
    else:
        # --- Remote: download archive ---
        archive_name = name + ".tar.gz"
        archive = cache / archive_name
        if _is_s3_uri(source_str):
            _copy_s3(source_str, archive)
        else:
            _download_http(source_str, archive)

    # --- Hash verification on archive ---
    if expected_hash is not None:
        if not verify_hash(archive, expected_hash, hash_algorithm):
            archive.unlink(missing_ok=True)
            raise ValueError(
                f"Hash verification failed for {archive} "
                f"(expected {hash_algorithm}:{expected_hash})"
            )
        logger.info("Hash verified (%s): %s", hash_algorithm, expected_hash)

    # --- Extract ---
    with tempfile.TemporaryDirectory(dir=cache) as tmp:
        extracted = _extract_tarball(archive, Path(tmp))
        # Move extracted DB to final location
        if db_dest.exists():
            shutil.rmtree(db_dest)
        shutil.move(str(extracted), str(db_dest))

    # Clean up archive
    if archive.exists() and archive.parent == cache:
        archive.unlink()

    # Validate the result
    validated = validate_database(db_dest)
    info = database_info(validated)
    logger.info(
        "Database ready: %s (%s, %d files)",
        validated,
        _human_size(info["total_size_bytes"]),
        len(info["files"]),
    )
    return validated


# ---------------------------------------------------------------------------
# Taxonomy validation
# ---------------------------------------------------------------------------


def validate_taxonomy(db_path: str | Path) -> dict[str, bool]:
    """Check for expected taxonomy files in a Kraken2 database.

    PrackenDB (and other well-formed databases) include
    ``taxonomy/nodes.dmp`` and ``taxonomy/names.dmp``.  These files enable
    lineage-aware classification (e.g. tracing a species taxid back to its
    domain) and human-readable taxon names.

    Parameters
    ----------
    db_path:
        Path to the Kraken2 database directory.

    Returns
    -------
    dict[str, bool]
        Mapping of taxonomy file relative paths to presence status.
    """
    db = Path(db_path).resolve()
    result: dict[str, bool] = {}
    for relpath in TAXONOMY_FILES:
        result[relpath] = (db / relpath).is_file()
    return result


def is_prackendb(db_path: str | Path) -> bool:
    """Heuristic check for whether *db_path* looks like a PrackenDB database.

    Returns ``True`` when the database directory is valid *and* contains
    the expected taxonomy files (``taxonomy/nodes.dmp`` and
    ``taxonomy/names.dmp``), which are characteristic of PrackenDB.

    This is a lightweight heuristic — it does not verify that the database
    truly contains one genome per species.
    """
    db = Path(db_path).resolve()
    try:
        validate_database(db)
    except (FileNotFoundError, ValueError):
        return False
    tax = validate_taxonomy(db)
    return all(tax.values())


def _emit_non_prackendb_warning(db_path: Path) -> None:
    """Log a warning if *db_path* does not look like PrackenDB."""
    if not is_prackendb(db_path):
        logger.warning(_NON_PRACKENDB_WARNING)


# ---------------------------------------------------------------------------
# PrackenDB convenience wrapper
# ---------------------------------------------------------------------------


def fetch_prackendb(
    *,
    url: str = PRACKENDB_URL,
    name: str = PRACKENDB_NAME,
    cache_dir: str | Path | None = None,
    expected_hash: str | None = None,
    hash_algorithm: str = "sha256",
) -> Path:
    """Download and cache the PrackenDB Kraken2 database.

    This is a convenience wrapper around :func:`fetch_database` with
    PrackenDB-specific defaults.  After fetching, the taxonomy files are
    validated and any missing files are reported as warnings.

    Parameters
    ----------
    url:
        Download URL for the PrackenDB tarball.
    name:
        Name for the database inside the cache directory.
    cache_dir:
        Override cache directory.
    expected_hash:
        Optional expected hash of the archive.
    hash_algorithm:
        Hash algorithm for verification.

    Returns
    -------
    Path
        Path to the ready-to-use PrackenDB database directory.
    """
    db_path = fetch_database(
        url,
        name=name,
        cache_dir=cache_dir,
        expected_hash=expected_hash,
        hash_algorithm=hash_algorithm,
    )

    # Validate taxonomy files
    tax = validate_taxonomy(db_path)
    for relpath, present in tax.items():
        if not present:
            logger.warning(
                "PrackenDB: taxonomy file not found: %s/%s. "
                "Lineage-aware classification may be degraded.",
                db_path,
                relpath,
            )
        else:
            logger.info("PrackenDB: taxonomy file OK: %s/%s", db_path, relpath)

    return db_path


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _dir_size(path: Path) -> int:
    """Total size of all files under *path* in bytes."""
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def _human_size(nbytes: int) -> str:
    """Return a human-friendly file-size string."""
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(nbytes) < 1024:
            return f"{nbytes:.1f} {unit}"
        nbytes /= 1024  # type: ignore[assignment]
    return f"{nbytes:.1f} PiB"
