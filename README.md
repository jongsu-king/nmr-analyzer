# NMR Analyzer

A desktop viewer and analyser for **1D and 2D NMR spectra**. Reads Bruker
experiment folders, zipped Bruker data, ACD/Labs `.esp` exports and JCAMP-DX
files; processes raw FIDs; and does peak picking, multiplet analysis,
lineshape fitting, integration with conversion and molar-ratio readout, 2D
contour display with cross-peak picking, and a SMILES structure check against
the integrals.

Pure Python with tkinter — **no third-party packages.**

## Install and run

Nothing to install: clone the repository and run it.

```bash
git clone https://github.com/jongsu-king/nmr-analyzer.git
cd nmr-analyzer
python3 -m nmranalyzer                          # empty, then File > Open
python3 -m nmranalyzer spectrum.esp data.zip    # or open files directly
```

On macOS you can double-click **`run.command`** in Finder; on Windows,
**`run.bat`**; on Linux, `./run.sh`.

To put `nmr-analyzer` on your PATH instead:

```bash
pip install .
nmr-analyzer
```

Or install a built wheel from the
[latest release](https://github.com/jongsu-king/nmr-analyzer/releases/latest):

```bash
pip install nmranalyzer-1.2.0-py3-none-any.whl
```

Requires Python 3.8 or newer with tkinter, which ships with the python.org and
Homebrew builds. On Debian/Ubuntu: `sudo apt install python3-tk`.

## Tests

```bash
python3 run_tests.py
```

81 tests covering the format traps (the `E`/SQZ ambiguity, DIF checkpoints,
Bruker submatrix de-tiling, the `.esp` axis flip) and the numerical claims
(FFT against a direct DFT, deconvolution areas, aromatic hydrogen counts).
They generate their own data, so no fixtures are committed.

## Supported formats

| Input | What is read |
|---|---|
| Bruker folder or `.zip` | every experiment containing `acqus`; raw `fid`, processed `pdata/1/1r` + `1i`, `procs`, `title`; 2D via `acqu2s` + `pdata/1/2rr` |
| ACD/Labs `.esp` | real and imaginary spectrum, referenced ppm axis, acquisition metadata |
| JCAMP-DX `.jdx` / `.dx` | `XYDATA` and NTUPLES `DATA TABLE` in PAC, SQZ, DIF and DUP forms; multi-block LINK files; FIDs |

A zip containing several experiment numbers loads all of them at once. If
`pdata/1` is absent — because the data was reprocessed into `pdata/2` or later
— the lowest processing number that actually holds data is used, and the Info
tab records which one.

## Workflow

1. **Open** the data. Each spectrum gets a colour and a row in the list; click
   the *On* column to toggle it in the plot.
2. **Process** (left panel). Line broadening and zero filling need the raw FID;
   phase correction and baseline correction work on any spectrum. *Auto Phase*
   fits PH0/PH1 by minimising the negative spectral area. *Apply to all
   spectra* propagates the settings across the whole list, which is what you
   want before an overlay comparison.
3. **Navigate.** Drag to zoom, right-drag to pan, wheel to zoom, shift+wheel
   for intensity, double-click to reset. A crosshair follows the pointer with
   a live ppm readout, and dragging shows the width of the selection in both
   ppm and Hz. The view is clamped to the data, so zooming out cannot lose the
   spectrum off-screen.
4. **Pick peaks** — operates on the *displayed* region, so zoom into the range
   you care about first. Raise *Sensitivity* to reject more noise.
5. **Integrate.** Either drag across a signal in Integrate mode, or press
   *Auto Integrate* to group the picked peaks into regions automatically.
6. **Assign.** Select a region, type its proton count and press *Assign*. Every
   other integral is then reported relative to it. Regions you have not
   assigned are reported as relative integrals, never as invented `nH` values.
7. **Report / export.** The Report tab writes the experimental-section string
   plus composition; File > Export writes peak and integral CSVs, an SVG plot
   or a PostScript plot.

### Calibration

Press *Calibrate to solvent peak* and the residual solvent line is moved to its
book value (Fulmer *et al.*, *Organometallics* **2010**, *29*, 2176). The
solvent is taken from the file metadata and recognised in both short and long
spellings (`TFA-d`, `TRIFLUOROACETIC ACID-d`).

For anything else, switch to Reference mode, click a peak (it snaps to the
nearest maximum), type its true shift and press Apply. Either way the axis
moves and existing peaks and integration regions move with it, so integrals
stay attached to their signals.

### Lineshape fitting

Select a region and press *Fit lines*. Overlapping signals are fitted as a sum
of pseudo-Voigt lines and the Fit tab lists each resolved line with its shift,
FWHM and area; the fitted envelope and its components are drawn over the
spectrum.

Use this when two multiplets overlap, because a plain integral then depends on
where you put the boundary. Note that fitted areas **exclude the local
baseline**, so they run below the plain integral — compare fitted areas with
each other, not with the integral column. The Fit tab reports the residual as a
percentage of the tallest line; anything much above ~10 % means the model did
not describe the data well.

### Conversion and mole ratio

Assign proton counts to two or more regions and the Report tab adds a
composition block: integral-per-proton, mol % for each species, and the
conversion computed by treating the **highest-shift assigned region** as
starting material. With three or more, the product ratio is reported too.

### 2D spectra

Open a Bruker folder or zip that contains an `acqu2s` and the 2D viewer opens
in its own window; further ones are listed under `Tools > 2D Spectra`.

Contours are recomputed for whatever region is on screen at the resolution of
the screen, so zooming in reveals detail rather than magnifying a coarse grid.
*Base level %* sets the lowest contour as a fraction of the strongest point,
and *Lower* / *Raise* step it geometrically — the usual way to dig a weak
correlation out of the noise. Negative contours are drawn in red, which is
what separates the phase-sensitive HSQC/HMBC responses from the positive ones.

Skyline projections run along the top (F2) and left (F1), and homonuclear
experiments get a diagonal guide. *Pick Peaks* lists the cross peaks and, for
homonuclear spectra, automatically excludes the diagonal.

### Structure check

`Tools > Structure` takes a SMILES string, draws the structure and reports the
formula, molecular mass, DBE and ring count, plus the ¹H and ¹³C environments
the connectivity predicts and an estimated shift for each.

Press *Check against spectrum* and the integrals of the active spectrum are
scaled so their total equals the proton count of that formula. Each region is
then matched to an environment on **both** its size and its predicted shift.

Three caveats, all deliberate:

* Integral size alone is **not** the test. With a free scale factor almost any
  structure can be made to land near integers — an early version of this scored
  toluene higher than ethanol on a real ethanol spectrum. Requiring the shift
  to agree as well is what makes a wrong structure fail: a methyl-sized
  integral sitting at 7.9 ppm cannot be an aliphatic CH₃ however neatly the
  areas divide.
* Shifts are additivity estimates, roughly ±0.35 ppm for ¹H and ±5–12 ppm for
  ¹³C, and worse for anything strained, charged or unusually conjugated. Every
  estimate is shown as a window rather than a single number. Multiply
  substituted carbons are the weak spot: chloroform comes out at 8.1 against a
  real 7.26, because three chlorine α-effects do not simply add.
* Equivalence is topological, so diastereotopic protons are not split apart
  and hindered rotation is ignored — DMF's two N-methyls count as one
  environment. That needs stereochemistry; this only knows the connectivity.

It is a consistency check, not proof: overlapping signals and regions you
forgot to integrate look the same as a wrong structure, and the window says so.

### Sessions

*File > Save Session* writes a `.nmrs` JSON file holding calibration,
processing settings, integration regions and proton assignments — but not the
spectra themselves, which are re-read from their original paths on open. Move
the raw data and the session will tell you which sources it could not find.

## Files

All modules live in the `nmranalyzer/` package.

| File | Contents |
|---|---|
| `app.py` | tkinter GUI, plotting, interaction |
| `nmrio.py` | format readers and the `Spectrum` model |
| `dsp.py` | FFT, apodisation, phasing, baseline, noise estimation |
| `analysis.py` | peak picking, integration, multiplets, composition |
| `fitting.py` | pseudo-Voigt deconvolution (Levenberg-Marquardt) |
| `solvents.py` | residual-solvent shifts and auto-calibration |
| `export.py` | SVG plot export and session save/load |
| `prefs.py` | recent files, window geometry, remembered defaults |
| `nmr2d.py` | Bruker 2D reading, de-tiling, cross-peak picking |
| `contour.py` | marching-squares contouring |
| `plot2d.py` | the 2D contour window |
| `smiles.py` | SMILES parser, formula/MW/DBE, symmetry classes |
| `depict.py` | 2D structure layout and drawing |
| `shifts.py` | additivity estimates of ¹H and ¹³C shifts |
| `structwin.py` | the structure window and the integral cross-check |

Plus `tests/` (the suite and its data generators), `run_tests.py`, and
`run.command` / `run.sh` / `run.bat` launchers.

## Documents

The app works on a session document. `File > Save` writes to the current
`.nmrs` file and only asks for a name the first time; the title bar carries a
`*` while there are unsaved changes, and New, Open and Quit all offer to save
first. `File > Open Recent` remembers the last ten files.

Analysis actions (picking, integrating, assigning, calibrating) go through a
40-step undo stack, and the Edit menu names what will be undone.

## Implementation notes

**Axis convention.** Every reader normalises to Bruker's layout — index 0 is
the left-hand (high ppm) edge, and `ppm(i) = OFFSET - i * SW / (SF * SI)`.

**FID processing.** `dsp.transform` applies apodisation, zero-fills, scales the
first point by `FCOR`, transforms, then removes the Bruker digital-filter group
delay by multiplying bin *k* by `exp(2πi·k·GRPDLY/N)` — a circular shift that,
unlike an integer roll, handles the fractional part of `GRPDLY`. Output bin *i*
sits at frequency `(N/2 − i)·SW/N`.

Validated against TopSpin: transforming the raw `fid` of a 500 MHz dataset and
matching TopSpin's own `PHC0` reproduces the shipped `1r` with zero point
offset across all 65536 points.

**Baseline.** The default `spline` method interpolates between many local
anchors (the median of the lowest 20 % of each segment, then median-smoothed),
which follows the rolling baseline left by the tail of an intense solvent line.
A global `poly` fit cannot, and is offered only for gently tilted baselines. On
TFA-d spectra where the solvent line is ~1000× the analyte, this lifts the
aromatic region from S/N ≈ 10 to S/N ≈ 115.

**Peak picking.** A candidate must exceed `sensitivity × noise` *and* be
prominent by the same margin, where prominence is the rise above the higher of
its two neighbouring valleys. Without the prominence test, ripples riding on a
broad line are picked as dozens of spurious "peaks". Apex positions are refined
by parabolic interpolation.

**Multiplets.** Equal line spacings plus Pascal-triangle intensities give
s/d/t/q/quint/sext/sept; four equal lines with spacings *J₂, J₁−J₂, J₂* give a
`dd`; six lines in a 2-2 pattern give a `dt`. Anything else stays `m` rather
than being forced into a pattern the data does not support.

**Fitting.** A sum of pseudo-Voigt lines — `η·Lorentzian + (1−η)·Gaussian` with
a shared η — plus a linear baseline, optimised by Levenberg-Marquardt with
analytic derivatives. Areas are computed in closed form, so they include the
intensity in the tails that lies outside the integration window.

Validated on synthetic data: two Lorentzians 4 Hz apart with 3 Hz linewidths
and a 2:1 area ratio are recovered to better than 0.25 % in area and 0.001 Hz
in width, where plain summation over the same window is 5.5 % low and cannot
separate them at all.

### JCAMP-DX

Validated against 10 real files: the Robert J. Lancashire conformance suite
(which exercises every ASDF form) and an ACD-exported 1H spectrum. All ten
decode to exactly their declared `##NPOINTS`, and the NMR file reproduces its
declared `FIRSTY`, `MAXY` and `MINY`.

Four things that a synthetic test would never have caught:

* In SQZ and DIF files the abscissa runs straight into the first ordinate with
  no separator (`5000.03B1399T...`), so the line cannot be split on
  whitespace. The whole line is decoded as ASDF and the leading value dropped.
* Abscissae carry decimal points, which the ordinate alphabets do not.
* **`E` is ambiguous**: it is the SQZ digit 5 *and* the exponent marker, so
  `1946E434` reads as 1946×10⁴³⁴ = infinity if AFFN exponents are honoured.
  JCAMP data tables never use scientific notation, so exponents are not parsed.
* Real files are often multi-block `LINK` documents with structure, assignment
  and peak blocks alongside the data; the spectrum block has to be selected.

`blckpac1.jdx` appears to mismatch but does not: its header declares
`MAXY= -.006136` and `MINY= .19`, i.e. a maximum below its minimum. The first
ordinate in its data really is −0.006136, so the file's own labels are wrong.

### Bruker 2D submatrix layout

`pdata/N/2rr` is **not** stored row by row. It is cut into tiles of `XDIM`
points along F2 (from `procs`) by `XDIM` points along F1 (from `proc2s`), and
the tiles are written one after another, row-major both between and within
tiles. Reading the file as a plain matrix gives a scrambled spectrum.

Verified against a synthetic dataset with peaks at known shifts: the de-tiled
matrix is bit-identical to the source, cross peaks land within 0.008 ppm, and
reading the same file row-major differs by the full peak height — so the test
would have caught a wrong layout.

**SMILES and hydrogen counting.** Implicit hydrogens come from standard
valences, with one subtlety that is easy to get wrong: an aromatic atom
carrying a formal double bond in the Kekule structure uses one more bond than
its sigma count. Carbon always does; nitrogen only does in the pyridine sense,
because once it has three sigma neighbours it must be donating its lone pair
instead. Adding the bonus unconditionally invents a hydrogen on every
N-substituted pyrrole-type nitrogen — it turned caffeine into C8H13N4O2.
Checked against 22 structures including pyridine, N-methylpyrrole, furan,
thiophene, indole and caffeine.

**Depiction.** Rings are regular polygons; a fused ring is reflected across the
bond it shares with one already placed; separate ring systems are rotated so
they grow away from the bond that reached them. The second line of a double or
aromatic bond must be drawn on the *inside* of its ring — offsetting it to a
fixed side of the bond vector makes fused aromatics look like buckled polygons.

**PNG vs SVG.** There is no PNG export. Rasterising the plot would mean
rendering text, and tkinter gives no access to a font engine that can write
into an image — a hand-rolled PNG would lose every label. SVG keeps labels as
real text, scales without blurring, handles non-ASCII spectrum titles, and
imports directly into Word, PowerPoint and Illustrator.

### The ACD `.esp` format

Undocumented, so this was worked out from the files themselves. After the two
length-prefixed strings `(C) ACD 1994` and `.ESP.( V 1.0 )` comes a simple
tag/length/value header — one byte of tag, one byte of length, then that many
bytes — terminated by a zero tag. The tags that matter:

| Tag | Meaning |
|---|---|
| `0x03` | spectral width, Hz (float32) |
| `0x04` | spectrum centre, Hz from the reference (float32) |
| `0x05` | observe frequency, MHz (float32) |
| `0x10` | **half** the number of stored real points (int32) |
| `0x0A` `0x11` `0x19` `0x0B` `0x0E` | title, solvent, pulse program, date, temperature |

Two float32 blocks follow the header — real then imaginary, each of
`2 × tag 0x10` points — ordered **low to high ppm**, the reverse of Bruker. The
axis is `right = tag0x04/SF − SW/(2·SF)`, ascending.

Cross-check: three unrelated `.esp` files from this instrument put the TFA-d
line at 11.4998, 11.5004 and 11.5001 ppm under this reading, i.e. ACD had
referenced all three to the solvent at exactly 11.50.

## Author

Jongsu Lim

## Licence

MIT — see [LICENSE](LICENSE).
