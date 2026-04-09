"""Command-line interface for Kraken2 database management.

Usage examples::

    # Fetch the recommended PrackenDB database
    csc-db fetch prackendb

    # Fetch a database from a URL
    csc-db fetch https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb.tar.gz

    # Fetch with hash verification
    csc-db fetch https://example.com/db.tar.gz --sha256 abc123...

    # Use a local database directory
    csc-db fetch /data/kraken2/PlusPF --name PlusPF

    # List cached databases
    csc-db list

    # Show database info
    csc-db info /data/kraken2/PlusPF

    # Verify a database
    csc-db verify /data/kraken2/PlusPF

    # Clean the cache
    csc-db clean
    csc-db clean --name old_db
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from csc import __version__
from csc.classify.db import (
    PRACKENDB_NAME,
    PRACKENDB_URL,
    clean_cache,
    database_info,
    fetch_database,
    fetch_prackendb,
    get_cache_dir,
    is_prackendb,
    list_databases,
    validate_taxonomy,
    verify_hash,
)
from csc.classify.classify import validate_database
from csc.utils import setup_logging


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="csc-db",
        description=(
            "Manage Kraken2 databases for the CSC classification module. "
            "Download, verify, cache, and inspect Kraken2 databases."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "--json-log",
        action="store_true",
        help="Emit structured JSON log lines.",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose (DEBUG) logging.",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help=(
            "Override the database cache directory. "
            "Default: CSC_DB_CACHE env var or ~/.csc/db"
        ),
    )

    sub = parser.add_subparsers(dest="command", help="Sub-command to run.")

    # -- fetch --
    fetch = sub.add_parser(
        "fetch",
        help="Fetch a Kraken2 database from a local path, URL, or S3 URI.",
        description=(
            "Fetch a Kraken2 database.  Use 'csc-db fetch prackendb' to "
            "download the recommended PrackenDB database (one genome per "
            "species).  Other sources (local paths, HTTP(S) URLs, S3 URIs) "
            "are also supported."
        ),
    )
    fetch.add_argument(
        "source",
        help=(
            "Database source: 'prackendb' (recommended), local directory, "
            "HTTP(S) URL to a .tar.gz, or S3 URI (s3://bucket/key.tar.gz)."
        ),
    )
    fetch.add_argument(
        "--name",
        default=None,
        help="Name for the database in the cache. Defaults to source basename.",
    )
    fetch.add_argument(
        "--sha256",
        default=None,
        help="Expected SHA-256 hash of the archive for verification.",
    )
    fetch.add_argument(
        "--md5",
        default=None,
        help="Expected MD5 hash of the archive for verification.",
    )

    # -- verify --
    verify = sub.add_parser(
        "verify",
        help="Verify that a database directory is a valid Kraken2 database.",
    )
    verify.add_argument("path", help="Path to the Kraken2 database directory.")

    # -- info --
    info = sub.add_parser(
        "info",
        help="Show metadata about a Kraken2 database.",
    )
    info.add_argument("path", help="Path to the Kraken2 database directory.")
    info.add_argument(
        "--json",
        action="store_true",
        dest="json_output",
        help="Output in JSON format.",
    )

    # -- list --
    sub.add_parser(
        "list",
        help="List databases in the cache directory.",
    )

    # -- clean --
    clean = sub.add_parser(
        "clean",
        help="Remove databases from the cache.",
    )
    clean.add_argument(
        "--name",
        default=None,
        help="Remove only the named database. Omit to remove all.",
    )

    return parser


def _cmd_fetch(args: argparse.Namespace, log: logging.Logger) -> int:
    expected_hash = None
    algorithm = "sha256"
    if args.sha256:
        expected_hash = args.sha256
        algorithm = "sha256"
    elif args.md5:
        expected_hash = args.md5
        algorithm = "md5"

    try:
        # Handle "prackendb" as a special convenience source
        if args.source.lower() == "prackendb":
            db_path = fetch_prackendb(
                cache_dir=args.cache_dir,
                expected_hash=expected_hash,
                hash_algorithm=algorithm,
            )
        else:
            db_path = fetch_database(
                args.source,
                name=args.name,
                cache_dir=args.cache_dir,
                expected_hash=expected_hash,
                hash_algorithm=algorithm,
            )
            # Emit warning for non-PrackenDB databases
            if not is_prackendb(db_path):
                from csc.classify.db import _NON_PRACKENDB_WARNING

                log.warning(_NON_PRACKENDB_WARNING)

        # Report taxonomy status
        tax = validate_taxonomy(db_path)
        for relpath, present in tax.items():
            status = "OK" if present else "MISSING"
            print(f"  Taxonomy: {relpath} [{status}]")

        print(f"Database ready: {db_path}")
        return 0
    except Exception as exc:
        log.error("Fetch failed: %s", exc)
        return 1


def _cmd_verify(args: argparse.Namespace, log: logging.Logger) -> int:
    try:
        db = validate_database(args.path)
        info = database_info(db)
        pracken = is_prackendb(db)
        print(f"Valid Kraken2 database: {db}")
        print(f"  PrackenDB-compatible: {'yes' if pracken else 'no'}")
        print(f"  Files: {len(info['files'])}")
        print(f"  Total size: {info['total_size_bytes']} bytes")
        tax = validate_taxonomy(db)
        for relpath, present in tax.items():
            status = "OK" if present else "MISSING"
            print(f"  Taxonomy: {relpath} [{status}]")
        for fname, digest in info["sha256"].items():
            print(f"  {fname}: sha256:{digest}")
        if not pracken:
            from csc.classify.db import _NON_PRACKENDB_WARNING

            log.warning(_NON_PRACKENDB_WARNING)
        return 0
    except (FileNotFoundError, ValueError) as exc:
        log.error("Verification failed: %s", exc)
        return 1


def _cmd_info(args: argparse.Namespace, log: logging.Logger) -> int:
    try:
        info = database_info(args.path)
    except Exception as exc:
        log.error("Info failed: %s", exc)
        return 1

    # Enrich with taxonomy/PrackenDB data
    db = Path(args.path).resolve()
    tax = validate_taxonomy(db) if db.is_dir() else {}
    pracken = is_prackendb(db) if db.is_dir() else False
    info["taxonomy"] = tax
    info["prackendb_compatible"] = pracken

    if getattr(args, "json_output", False):
        print(json.dumps(info, indent=2, default=str))
    else:
        print(f"Database: {info['path']}")
        print(f"  Valid: {info['valid']}")
        print(f"  PrackenDB-compatible: {'yes' if pracken else 'no'}")
        print(f"  Total size: {info['total_size_bytes']} bytes")
        for relpath, present in tax.items():
            status = "OK" if present else "MISSING"
            print(f"  Taxonomy: {relpath} [{status}]")
        for fname, size in info["files"].items():
            digest = info["sha256"].get(fname, "n/a")
            print(f"  {fname}: {size} bytes (sha256:{digest})")
    return 0


def _cmd_list(args: argparse.Namespace, log: logging.Logger) -> int:
    dbs = list_databases(args.cache_dir)
    if not dbs:
        print(f"No databases in cache: {get_cache_dir(args.cache_dir)}")
        return 0
    print(f"Databases in {get_cache_dir(args.cache_dir)}:")
    for db in dbs:
        status = "valid" if db["valid"] else "INVALID"
        print(f"  {db['name']}: {db['size_bytes']} bytes [{status}]")
    return 0


def _cmd_clean(args: argparse.Namespace, log: logging.Logger) -> int:
    removed = clean_cache(args.cache_dir, name=args.name)
    if removed:
        print(f"Removed: {', '.join(removed)}")
    else:
        print("Nothing to remove.")
    return 0


def main(argv: list[str] | None = None) -> int:
    """CLI entry point.  Returns 0 on success, 1 on failure."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    setup_logging(
        level="DEBUG" if args.verbose else "INFO",
        json_format=args.json_log,
    )
    log = logging.getLogger(__name__)

    if args.command is None:
        parser.print_help()
        return 1

    dispatch = {
        "fetch": _cmd_fetch,
        "verify": _cmd_verify,
        "info": _cmd_info,
        "list": _cmd_list,
        "clean": _cmd_clean,
    }
    handler = dispatch.get(args.command)
    if handler is None:
        parser.print_help()
        return 1

    return handler(args, log)


if __name__ == "__main__":
    sys.exit(main())
