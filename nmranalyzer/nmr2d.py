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


def read_bruker_2d(path, procno="1"):
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
    pdata = "%spdata/%s/" % (root, procno)
    if not (store.exists(pdata + "2rr") and store.exists(pdata + "procs")
            and store.exists(pdata + "proc2s")):
        return None

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
        "Submatrix": "%d x %d" % (tile_rows, tile_cols),
    }
    return Spectrum2D(name, data, f2, f1, meta=meta, source=path)


def load_2d(path):
    """Load 2D spectra from a folder or zip; empty list if there are none."""
    if not (os.path.isdir(path) or path.lower().endswith(".zip")):
        return []
    try:
        return read_bruker_2d(path)
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
