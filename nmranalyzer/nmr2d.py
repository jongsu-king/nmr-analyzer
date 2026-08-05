"""Bruker 2D NMR: reading processed data and picking cross peaks.

The awkward part of the format is that ``pdata/N/2rr`` is not stored row by
row.  It is cut into submatrices ("tiles") of ``XDIM`` points along F2 by
``XDIM`` points along F1 -- the F2 tile width comes from ``procs`` and the F1
tile height from ``proc2s`` -- and those tiles are written one after another.
Reading the file as a plain matrix produces a scrambled spectrum, so
:func:`_untile` is the piece that has to be right.

Axis convention matches the 1D side: index 0 is the high-ppm edge of each
dimension.  F2 is the direct (horizontal) axis, F1 the indirect (vertical) one.
"""

from __future__ import annotations

import os
import struct

from . import dsp
from . import nmrio


class Axis:
    """One frequency axis of a 2D spectrum."""

    def __init__(self, size, sf, sw_hz, offset_ppm, label="", nucleus=""):
        self.size = size
        self.sf = sf                    # MHz
        self.sw_hz = sw_hz
        self.offset_ppm = offset_ppm    # ppm of index 0
        self.label = label
        self.nucleus = nucleus

    @property
    def sw_ppm(self):
        return self.sw_hz / self.sf if self.sf else 0.0

    @property
    def delta(self):
        return -self.sw_ppm / self.size if self.size else 0.0

    def ppm(self, index):
        return self.offset_ppm + index * self.delta

    def index(self, ppm):
        return 0 if not self.delta else int(round((ppm - self.offset_ppm) / self.delta))

    def clamp(self, index):
        return max(0, min(self.size - 1, index))

    @property
    def limits(self):
        return self.ppm(0), self.ppm(self.size - 1)


class Spectrum2D:
    """A processed 2D spectrum: ``data[f1_index][f2_index]``."""

    def __init__(self, name, data, f2, f1, meta=None, source=""):
        self.name = name
        self.data = data
        self.f2 = f2                    # direct dimension, horizontal
        self.f1 = f1                    # indirect dimension, vertical
        self.meta = meta or {}
        self.source = source
        self.peaks = []                 # list of Peak2D

        self._noise = None
        self._max = None

    @property
    def rows(self):
        return self.f1.size

    @property
    def cols(self):
        return self.f2.size

    def value(self, f1_index, f2_index):
        return self.data[self.f1.clamp(f1_index)][self.f2.clamp(f2_index)]

    def max_intensity(self):
        if self._max is None:
            self._max = max(max(row) for row in self.data)
        return self._max

    def noise(self):
        """Robust noise estimate from the quietest corner of the matrix.

        The corners of a 2D map are almost always signal-free, which makes
        them a better noise reference than any global statistic.
        """
        if self._noise is not None:
            return self._noise
        rows, cols = self.rows, self.cols
        rh, cw = max(4, rows // 16), max(4, cols // 16)
        best = None
        for r0 in (0, rows - rh):
            for c0 in (0, cols - cw):
                block = [v for r in range(r0, r0 + rh)
                         for v in self.data[r][c0:c0 + cw]]
                if not block:
                    continue
                mean = sum(block) / len(block)
                var = sum((v - mean) ** 2 for v in block) / len(block)
                if best is None or var < best:
                    best = var
        self._noise = (best or 0.0) ** 0.5
        return self._noise

    def is_homonuclear(self):
        a = (self.f2.nucleus or "").strip()
        b = (self.f1.nucleus or "").strip()
        return bool(a) and a == b


class Peak2D:
    def __init__(self, f2_ppm, f1_ppm, intensity):
        self.f2_ppm = f2_ppm
        self.f1_ppm = f1_ppm
        self.intensity = intensity

    def __repr__(self):
        return "Peak2D(F2 %.3f, F1 %.3f)" % (self.f2_ppm, self.f1_ppm)


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def _untile(flat, rows, cols, tile_rows, tile_cols):
    """Reassemble Bruker's submatrix layout into a normal matrix.

    Tiles are stored in row-major order of tiles, and row-major within each
    tile.  Any partial tile at the right or bottom edge is tolerated.
    """
    matrix = [[0.0] * cols for _ in range(rows)]
    if tile_rows <= 0 or tile_cols <= 0:
        tile_rows, tile_cols = rows, cols

    tiles_down = (rows + tile_rows - 1) // tile_rows
    tiles_across = (cols + tile_cols - 1) // tile_cols

    position = 0
    limit = len(flat)
    for tile_row in range(tiles_down):
        for tile_col in range(tiles_across):
            base_r = tile_row * tile_rows
            base_c = tile_col * tile_cols
            for r in range(tile_rows):
                if position >= limit:
                    return matrix
                row_slice = flat[position:position + tile_cols]
                position += tile_cols
                target_r = base_r + r
                if target_r >= rows:
                    continue
                target = matrix[target_r]
                for c, value in enumerate(row_slice):
                    target_c = base_c + c
                    if target_c < cols:
                        target[target_c] = value
    return matrix


def _axis_from_procs(procs, acqus, label):
    size = int(procs.get("SI", 0))
    sf = float(procs.get("SF", 0.0)) or float(acqus.get("SFO1", 0.0))
    sw_hz = float(procs.get("SW_p", 0.0)) or float(acqus.get("SW_h", 0.0))
    offset = float(procs.get("OFFSET", 0.0))
    return Axis(size, sf, sw_hz, offset, label=label,
                nucleus=str(acqus.get("NUC1", "")))


def find_2d_experiments(store):
    """Folders that hold an ``acqu2s`` -- i.e. a second dimension."""
    roots = []
    for name in store.names:
        if name.endswith("acqu2s") and "/pdata/" not in name:
            roots.append(name[: -len("acqu2s")])
    return sorted(set(roots))


def read_bruker_2d(path, procno=None):
    """Read every 2D experiment found at ``path`` (folder or zip)."""
    store = nmrio._Store(path)
    try:
        out = []
        for root in find_2d_experiments(store):
            spec = _read_one(store, root, path, procno)
            if spec is not None:
                out.append(spec)
        return out
    finally:
        store.close()


def _read_one(store, root, path, procno):
    # Use the requested processing number, else the lowest one that actually
    # holds a 2rr matrix.
    candidates = [procno] if procno else nmrio.available_procnos(store, root) or [1]
    chosen = None
    for number in candidates:
        trial = "%spdata/%s/" % (root, number)
        if (store.exists(trial + "2rr") and store.exists(trial + "procs")
                and store.exists(trial + "proc2s")):
            chosen = number
            break
    if chosen is None:
        return None
    pdata = "%spdata/%s/" % (root, chosen)

    acqus = nmrio.parse_jcamp_params(store.read(root + "acqus").decode("latin-1"))
    acqu2s = nmrio.parse_jcamp_params(store.read(root + "acqu2s").decode("latin-1"))
    procs = nmrio.parse_jcamp_params(store.read(pdata + "procs").decode("latin-1"))
    proc2s = nmrio.parse_jcamp_params(store.read(pdata + "proc2s").decode("latin-1"))

    f2 = _axis_from_procs(procs, acqus, "F2")
    f1 = _axis_from_procs(proc2s, acqu2s, "F1")
    if not (f2.size and f1.size):
        return None

    big = int(procs.get("BYTORDP", 0)) == 1
    dtype = "float64" if int(procs.get("DTYPP", 0)) == 2 else "int32"
    raw = nmrio._unpack(store.read(pdata + "2rr"), dtype, big)
    scale = 2.0 ** float(procs.get("NC_proc", 0))
    if scale != 1.0:
        raw = [v * scale for v in raw]

    tile_cols = int(procs.get("XDIM", f2.size) or f2.size)
    tile_rows = int(proc2s.get("XDIM", f1.size) or f1.size)
    data = _untile(raw, f1.size, f2.size, tile_rows, tile_cols)

    title = ""
    if store.exists(pdata + "title"):
        lines = store.read(pdata + "title").decode("latin-1").strip().splitlines()
        title = lines[0] if lines else ""
    expno = root.rstrip("/").split("/")[-1]
    name = title or "%s [%s]" % (os.path.basename(path), expno)

    meta = {
        "Experiment": expno,
        "Pulse program": acqus.get("PULPROG", ""),
        "F2 nucleus": acqus.get("NUC1", ""),
        "F1 nucleus": acqu2s.get("NUC1", ""),
        "F2 size": f2.size,
        "F1 size": f1.size,
        "F2 frequency (MHz)": round(f2.sf, 4),
        "F1 frequency (MHz)": round(f1.sf, 4),
        "Solvent": acqus.get("SOLVENT", ""),
        "Scans": acqus.get("NS", ""),
        "Format": "Bruker 2D",
        "Processing no.": chosen,
        "Submatrix": "%d x %d" % (tile_rows, tile_cols),
    }
    return Spectrum2D(name, data, f2, f1, meta=meta, source=path)


def load_2d(path, allow_raw=True, progress=None):
    """Load 2D spectra from a folder or zip; empty list if there are none.

    Processed ``2rr`` data is preferred.  If an experiment has a second
    dimension but was never processed in TopSpin, the raw ``ser`` is
    transformed instead so the data is still usable.
    """
    if not (os.path.isdir(path) or path.lower().endswith(".zip")):
        return []
    try:
        processed = read_bruker_2d(path)
    except Exception:
        processed = []
    if processed or not allow_raw:
        return processed
    try:
        return read_bruker_ser(path, progress=progress)
    except Exception:
        return []


# ---------------------------------------------------------------------------
# Cross-peak picking
# ---------------------------------------------------------------------------


def pick_peaks_2d(spec, sensitivity=12.0, max_peaks=400, window=None,
                  skip_diagonal_ppm=0.0):
    """Find local maxima of the 2D map.

    ``window`` limits the search to ``(f2_lo, f2_hi, f1_lo, f1_hi)`` in ppm.
    ``skip_diagonal_ppm`` drops peaks within that distance of the diagonal,
    which is what you want when reading a COSY.
    """
    threshold = spec.noise() * sensitivity
    if threshold <= 0:
        threshold = spec.max_intensity() * 0.02

    if window:
        f2_lo, f2_hi, f1_lo, f1_hi = window
        c0 = spec.f2.clamp(spec.f2.index(f2_hi))
        c1 = spec.f2.clamp(spec.f2.index(f2_lo))
        r0 = spec.f1.clamp(spec.f1.index(f1_hi))
        r1 = spec.f1.clamp(spec.f1.index(f1_lo))
    else:
        c0, c1, r0, r1 = 0, spec.cols - 1, 0, spec.rows - 1
    if c0 > c1:
        c0, c1 = c1, c0
    if r0 > r1:
        r0, r1 = r1, r0

    found = []
    data = spec.data
    for r in range(max(r0, 1), min(r1, spec.rows - 2) + 1):
        row = data[r]
        above = data[r - 1]
        below = data[r + 1]
        for c in range(max(c0, 1), min(c1, spec.cols - 2) + 1):
            v = row[c]
            if v < threshold:
                continue
            if (v >= row[c - 1] and v > row[c + 1]
                    and v >= above[c] and v > below[c]
                    and v >= above[c - 1] and v >= below[c + 1]
                    and v >= above[c + 1] and v >= below[c - 1]):
                found.append((v, r, c))

    found.sort(key=lambda item: -item[0])
    peaks = []
    for v, r, c in found:
        f2_ppm = spec.f2.ppm(c)
        f1_ppm = spec.f1.ppm(r)
        if skip_diagonal_ppm and abs(f2_ppm - f1_ppm) <= skip_diagonal_ppm:
            continue
        peaks.append(Peak2D(f2_ppm, f1_ppm, v))
        if len(peaks) >= max_peaks:
            break
    return peaks


def projection(spec, axis="f2", mode="max"):
    """Collapse the matrix onto one axis (a skyline or summed projection)."""
    if axis == "f2":
        if mode == "sum":
            return [sum(spec.data[r][c] for r in range(spec.rows))
                    for c in range(spec.cols)]
        out = [float("-inf")] * spec.cols
        for row in spec.data:
            for c, v in enumerate(row):
                if v > out[c]:
                    out[c] = v
        return out
    if mode == "sum":
        return [sum(row) for row in spec.data]
    return [max(row) for row in spec.data]


# ---------------------------------------------------------------------------
# Raw 2D data (``ser``)
# ---------------------------------------------------------------------------

# Bruker ##$FnMODE in acqu2s: how the indirect dimension was detected.
FNMODE_NAMES = {
    0: "undefined", 1: "QF", 2: "QSEQ", 3: "TPPI",
    4: "States", 5: "States-TPPI", 6: "Echo-Antiecho",
}

SER_BLOCK = 1024        # every FID in a ser file starts on a 1024-byte boundary


def _read_ser_rows(raw, td1, td2, dtype, big_endian):
    """Split a ``ser`` file into rows of complex points.

    Each row is padded out to a 1024-byte boundary, so the stride is not
    simply the number of points; reading it as one flat array shears the
    spectrum along F1.
    """
    width = 8 if dtype == "float64" else 4
    row_bytes = ((td2 * width + SER_BLOCK - 1) // SER_BLOCK) * SER_BLOCK
    rows = []
    for index in range(td1):
        start = index * row_bytes
        chunk = raw[start:start + td2 * width]
        if len(chunk) < td2 * width:
            break
        values = nmrio._unpack(chunk, dtype, big_endian)
        rows.append([complex(values[i], values[i + 1])
                     for i in range(0, len(values) - 1, 2)])
    return rows


def _combine_f1(columns, fnmode):
    """Turn the detected rows into one complex interferogram per F2 point.

    ``columns`` is indexed ``[t1_row][f2_point]``; the return value is indexed
    ``[f2_point][t1_point]``.
    """
    n_rows = len(columns)
    if not n_rows:
        return []
    width = len(columns[0])

    if fnmode in (FNMODE_STATES := 4, 5, 6):
        pairs = n_rows // 2
        out = [[0j] * pairs for _ in range(width)]
        for k in range(pairs):
            first = columns[2 * k]
            second = columns[2 * k + 1]
            for f in range(width):
                if fnmode == 6:                       # Echo-Antiecho
                    p, n = first[f], second[f]
                    cosine = (p + n) / 2.0
                    sine = (p - n) / 2.0j
                else:
                    cosine, sine = first[f], second[f]
                value = complex(cosine.real, sine.real)
                if fnmode == 5 and k % 2:             # States-TPPI
                    value = -value
                out[f][k] = value
        return out

    if fnmode == 3:                                   # TPPI: real modulation
        out = [[0j] * n_rows for _ in range(width)]
        for k, row in enumerate(columns):
            for f in range(width):
                out[f][k] = complex(row[f].real, 0.0)
        return out

    # QF and anything unrecognised: treat the rows as already complex.
    out = [[0j] * n_rows for _ in range(width)]
    for k, row in enumerate(columns):
        for f in range(width):
            out[f][k] = row[f]
    return out


def process_ser(rows, sw2_hz, sw1_hz, si2, si1, fnmode=4, grpdly=0.0,
                lb2=0.3, lb1=0.3, magnitude_f1=False, progress=None):
    """Transform a raw 2D FID into a ``[f1][f2]`` matrix of real intensities.

    F2 is transformed first, row by row, exactly as for a 1D spectrum.  The
    indirect dimension is then assembled from the detected rows according to
    ``fnmode`` and transformed.
    """
    dwell2 = 1.0 / sw2_hz
    f2_spectra = []
    for index, row in enumerate(rows):
        spectrum = dsp.transform(row, dwell2, si2, lb=lb2, grpdly=grpdly)
        f2_spectra.append(spectrum)
        if progress and index % 16 == 0:
            progress("F2", index, len(rows))

    interferograms = _combine_f1(f2_spectra, fnmode)
    dwell1 = 1.0 / sw1_hz if sw1_hz else 1.0

    n_f1 = si1 if fnmode != 3 else si1 * 2
    columns = []
    for index, series in enumerate(interferograms):
        spectrum = dsp.transform(series, dwell1, n_f1, lb=lb1)
        if fnmode == 3:                       # TPPI keeps half the transform
            spectrum = spectrum[:si1]
        columns.append(spectrum)
        if progress and index % 64 == 0:
            progress("F1", index, len(interferograms))

    rows_out = []
    for r in range(min(si1, len(columns[0]) if columns else 0)):
        if magnitude_f1 or fnmode in (1, 2):
            rows_out.append([abs(column[r]) for column in columns])
        else:
            rows_out.append([column[r].real for column in columns])
    return rows_out


def read_bruker_ser(path, si2=None, si1=None, lb2=0.3, lb1=0.3,
                    progress=None):
    """Read and transform every raw 2D experiment at ``path``.

    Used when a dataset was acquired but never processed in TopSpin, so there
    is no ``2rr`` to read.
    """
    store = nmrio._Store(path)
    try:
        out = []
        for root in find_2d_experiments(store):
            if not store.exists(root + "ser"):
                continue
            spec = _read_ser_one(store, root, path, si2, si1, lb2, lb1,
                                 progress)
            if spec is not None:
                out.append(spec)
        return out
    finally:
        store.close()


def _read_ser_one(store, root, path, si2, si1, lb2, lb1, progress):
    acqus = nmrio.parse_jcamp_params(store.read(root + "acqus").decode("latin-1"))
    acqu2s = nmrio.parse_jcamp_params(store.read(root + "acqu2s").decode("latin-1"))

    td2 = int(acqus.get("TD", 0))
    td1 = int(acqu2s.get("TD", 0))
    sw2 = float(acqus.get("SW_h", 0.0))
    sw1 = float(acqu2s.get("SW_h", 0.0))
    sf2 = float(acqus.get("SFO1", 0.0))
    sf1 = float(acqu2s.get("SFO1", 0.0)) or sf2
    if not (td2 and td1 and sw2 and sf2):
        return None

    fnmode = int(acqu2s.get("FnMODE", 4) or 4)
    if fnmode in (0, 2):
        fnmode = 4                      # undefined in practice means States
    big = int(acqus.get("BYTORDA", 0)) == 1
    dtype = "float64" if int(acqus.get("DTYPA", 0)) == 2 else "int32"
    grpdly = float(acqus.get("GRPDLY", 0.0) or 0.0)

    rows = _read_ser_rows(store.read(root + "ser"), td1, td2, dtype, big)
    if not rows:
        return None

    size2 = si2 or dsp.next_pow2(len(rows[0]))
    pairs = len(rows) // 2 if fnmode in (4, 5, 6) else len(rows)
    size1 = si1 or dsp.next_pow2(max(pairs, 2))

    data = process_ser(rows, sw2, sw1 or sw2, size2, size1, fnmode=fnmode,
                       grpdly=grpdly, lb2=lb2, lb1=lb1, progress=progress)
    if not data:
        return None

    o1 = float(acqus.get("O1", 0.0))
    o2 = float(acqu2s.get("O1", o1))
    f2 = Axis(len(data[0]), sf2, sw2, o1 / sf2 + (sw2 / sf2) / 2.0,
              label="F2", nucleus=str(acqus.get("NUC1", "")))
    f1 = Axis(len(data), sf1, sw1 or sw2,
              o2 / sf1 + ((sw1 or sw2) / sf1) / 2.0,
              label="F1", nucleus=str(acqu2s.get("NUC1", "")))

    expno = root.rstrip("/").split("/")[-1]
    meta = {
        "Experiment": expno,
        "Pulse program": acqus.get("PULPROG", ""),
        "F2 nucleus": acqus.get("NUC1", ""),
        "F1 nucleus": acqu2s.get("NUC1", ""),
        "F2 size": f2.size,
        "F1 size": f1.size,
        "F2 frequency (MHz)": round(sf2, 4),
        "F1 frequency (MHz)": round(sf1, 4),
        "Solvent": acqus.get("SOLVENT", ""),
        "Scans": acqus.get("NS", ""),
        "Detection (FnMODE)": FNMODE_NAMES.get(fnmode, str(fnmode)),
        "Format": "Bruker 2D (from raw ser)",
    }
    name = "%s [%s] raw" % (os.path.basename(path), expno)
    return Spectrum2D(name, data, f2, f1, meta=meta, source=path)
