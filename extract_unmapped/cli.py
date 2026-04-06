"""Backward-compatibility shim – delegates to ``csc.extract.cli``."""

from csc.extract.cli import main  # noqa: F401

if __name__ == "__main__":
    import sys

    sys.exit(main())
