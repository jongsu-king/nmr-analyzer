# Contributing

Bug reports, spectra that fail to load, and patches are all welcome.

## Reporting a problem

Open an issue at
<https://github.com/jongsu-king/nmr-analyzer/issues> and say what you did,
what happened, and what you expected. For a file that will not open, the most
useful thing you can attach is the file itself, or — if the data is not yours
to share — the vendor parameter files (`acqus`, `procs`) with the raw data
removed, which is usually enough to reproduce a format problem.

Reports of a spectrum that *loads* but processes wrongly are especially
valuable. Two-dimensional support has so far only been checked against
synthetic data built to the format specification, so a real COSY or HSQC that
comes out wrong is a genuinely useful bug report.

## Making a change

```bash
git clone https://github.com/jongsu-king/nmr-analyzer.git
cd nmr-analyzer
python3 run_tests.py          # 98 tests, no third-party packages needed
python3 -m nmranalyzer        # run it
```

Two conventions hold throughout:

- **No third-party dependencies.** This is the point of the project. If a
  change needs NumPy, it belongs somewhere else.
- **Tests generate their own data.** No binary fixtures are committed; a test
  that needs a Bruker dataset writes one. See `tests/fixtures.py`.

Please add a test with any change to processing or format handling, and run
the suite before opening a pull request. Comments should explain why the code
is the way it is, particularly where a format has a trap in it — several of
the existing comments exist because the obvious implementation was wrong.

## Scope

The project aims to be a readable, installable tool for everyday 1D and 2D
work, not a complete replacement for a vendor package. Contributions that keep
it small and dependency-free are easier to accept than ones that broaden it.
