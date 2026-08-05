"""Entry point: ``python3 -m nmranalyzer``."""

import sys

from .app import main

if __name__ == "__main__":
    sys.exit(main())
