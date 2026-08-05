---
title: 'NMR Analyzer: a dependency-free desktop application for 1D and 2D NMR spectra'
tags:
  - Python
  - NMR
  - spectroscopy
  - analytical chemistry
  - JCAMP-DX
authors:
  - name: Jongsu Lim
    orcid: 0009-0000-3395-1484
    affiliation: 1
affiliations:
  - name: "AFFILIATION TO BE COMPLETED"   # e.g. Department, University, City, Country
    index: 1
date: 6 August 2026
bibliography: paper.bib
---

# Summary

`NMR Analyzer` is a desktop application for inspecting and quantifying nuclear
magnetic resonance spectra. It reads Bruker experiment folders, ACD/Labs `.esp`
exports and JCAMP-DX files; transforms raw free induction decays; and provides
phase and baseline correction, peak picking, integration, multiplet
classification with coupling constants, pseudo-Voigt deconvolution of
overlapping signals, contour display of two-dimensional experiments, and a
consistency check of measured integrals against a structure supplied as SMILES.

Its distinguishing property is that it has no dependencies. Everything —
the fast Fourier transform, the lineshape fitting, the contouring, the
structure layout, the graphical interface — is implemented against the Python
standard library. The program runs from a clone with `python3 -m nmranalyzer`,
and pre-built applications for Windows, macOS and Linux require no Python
installation at all.

# Statement of need

Processing an NMR spectrum is routine in synthetic chemistry: a chemist records
a spectrum, integrates it, and asks whether the result is consistent with the
compound they intended to make. The software for this divides into two groups,
and both leave a gap.

Vendor and commercial packages — TopSpin, MestReNova, ACD/Spectrus — are
capable but licensed, tied to particular machines, and closed. A student
cannot read how an integral is computed, and a laboratory cannot install the
software on an arbitrary computer.

Scriptable open-source libraries fill part of the gap. `nmrglue`
[@Helmus2013] is the reference implementation for reading vendor formats in
Python and is widely used. It is, however, a library rather than an
application: it requires NumPy and SciPy, and it assumes the user writes code.
For a chemist who wants to open a file, drag across a peak and read off an
integral, that is a different task from the one the library solves.

`NMR Analyzer` targets that second task. Because it carries no dependencies it
can be run on a locked-down instrument computer, a teaching laboratory machine,
or a personal laptop, without a package manager and without administrator
rights. Because every algorithm is written out in plain Python rather than
delegated to a compiled array library, the processing chain is readable: the
Fourier transform, the digital-filter correction and the integration are each
a few dozen lines that a student can follow.

# Functionality

**Formats.** Bruker experiment directories and archives are read, including
raw `fid` and `ser` data, processed `1r`/`1i` and `2rr`, and the acquisition
and processing parameter files. JCAMP-DX is supported in its PAC, SQZ, DIF and
DUP compressed forms [@Davies1993], along with `NTUPLES` tables and multi-block
`LINK` documents. The ACD/Labs `.esp` format is undocumented; its structure was
determined from the files themselves and is described in the repository
documentation.

**Processing.** Free induction decays are apodised, zero-filled and
transformed, with the Bruker digital-filter group delay removed by a
first-order phase correction that handles its fractional part. Phase
correction can be fitted automatically by minimising negative spectral area.
Baseline correction offers a polynomial fit and a spline through local anchors;
the latter matters when an intense solvent line leaves a rolling baseline
across the region of interest.

**Analysis.** Peaks are picked using both an absolute threshold and a
prominence test, which suppresses the ripples that sit on the flanks of broad
lines. Multiplets are classified from line spacings and intensity ratios, and
coupling constants reported. Overlapping signals can be resolved by fitting a
sum of pseudo-Voigt lines with Levenberg--Marquardt, giving areas that do not
depend on where an integration boundary was drawn. Two-dimensional data are
displayed as contours computed for the region on screen at the resolution of
the screen, with cross-peak picking and skyline projections.

**Structure.** A structure entered as SMILES [@Weininger1988] is parsed,
depicted, and reduced to its distinct proton and carbon environments.
Diastereotopic protons are separated, so that a CH~2~ adjacent to a
stereocentre counts as two signals and the isopropyl group of valine as two
three-proton signals rather than one six-proton one. Approximate shifts are
estimated by additivity, and the measured integrals are scaled to the
molecular formula and matched to predicted environments on both size and
shift. Regions can also be assigned to specific atoms by clicking.

# Quality control

The package ships 98 tests that generate their own data, so no binary fixtures
are committed. They cover the numerical claims and the format traps that a
synthetic test would not catch. Three are worth naming, because each was a
real defect found by the corresponding test:

- In JCAMP-DX, `E` is simultaneously the SQZ digit 5 and the exponent marker,
  so honouring scientific notation reads `1946E434` as $1946\times10^{434}$.
- Bruker stores processed 2D data in submatrices; reading the file as a plain
  matrix produces a scrambled spectrum.
- Each free induction decay in a Bruker `ser` file begins on a 1024-byte
  boundary, and an early version of the padding test used a row length that
  happened to be an exact multiple of that boundary, so it tested nothing.

Transforming a raw Bruker `fid` and matching the vendor's own zero-order phase
reproduces the processed spectrum TopSpin wrote with zero point offset across
all 65536 points. Deconvolution recovers the areas of two overlapping
Lorentzian lines to better than 0.25 %.

# Limitations

Two-dimensional support has been verified against synthetic data constructed
to the format specification, but no experimentally acquired 2D file has yet
been processed with it. Chemical-shift estimates are additivity rules and are
unreliable for multiply substituted carbons. Proton equivalence is determined
from connectivity and declared stereodescriptors, so hindered rotation is not
modelled. Diastereotopic protons are labelled H~a~/H~b~ rather than pro-*R*
and pro-*S*, which would require Cahn--Ingold--Prelog ranking.

# Acknowledgements

The residual-solvent reference shifts are those tabulated by Fulmer and
co-workers [@Fulmer2010].

# References
