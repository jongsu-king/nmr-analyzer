"""Entry point for the frozen builds.

PyInstaller starts from a plain script; a package's ``__main__`` cannot be
used as an entry point once frozen.
"""

import multiprocessing

if __name__ == "__main__":
    multiprocessing.freeze_support()
    from nmranalyzer.app import main
    main()
