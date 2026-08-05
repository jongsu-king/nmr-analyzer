"""Synthetic NMR data used by the tests.

Everything the suite needs is generated here rather than committed as binary
fixtures, so the tests stay self-contained and the expected values are visible
in the code that produces them.
"""

from __future__ import annotations

import cmath
import math
import os
import struct
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


# ---------------------------------------------------------------------------
# Bruker
# ---------------------------------------------------------------------------


def _write_params(path, pairs, title="synthetic"):
    with open(path, "w") as fh:
        fh.write("##TITLE= %s\n" % title)
        for key, value in pairs:
            fh.write("##$%s= %s\n" % (key, value))
        fh.write("##END=\n")


def write_bruker_1d(root, peaks=((7.26, 1.0), (2.10, 0.4)), si=8192,
                    sf=500.13, sw=5000.0, offset=12.0):
    """A processed Bruker 1D dataset with Lorentzian peaks at known shifts."""
    os.makedirs(os.path.join(root, "pdata", "1"), exist_ok=True)
    step = (sw / sf) / si
    data = []
    for i in range(si):
        ppm = offset - i * step
        value = 0.0
        for centre, height in peaks:
            u = 2.0 * (ppm - centre) * sf / 2.0        # 2 Hz linewidth
            value += height / (1.0 + u * u)
        data.append(int(value * 1e6))
    with open(os.path.join(root, "pdata", "1", "1r"), "wb") as fh:
        fh.write(struct.pack("<%di" % si, *data))

    _write_params(os.path.join(root, "acqus"),
                  [("SFO1", sf), ("SW_h", sw), ("NUC1", "<1H>"),
                   ("PULPROG", "<zg30>"), ("NS", 16), ("SOLVENT", "<CDCl3>"),
                   ("TD", si), ("BYTORDA", 0), ("DTYPA", 2), ("GRPDLY", 0)])
    _write_params(os.path.join(root, "pdata", "1", "procs"),
                  [("SI", si), ("SF", sf), ("SW_p", sw), ("OFFSET", offset),
                   ("NC_proc", 0), ("BYTORDP", 0), ("DTYPP", 0), ("LB", 0.3),
                   ("FCOR", 0.5)])
    with open(os.path.join(root, "pdata", "1", "title"), "w") as fh:
        fh.write("synthetic 1D\n")
    return root


COSY_PEAKS = [(7.20, 7.20, 1.0), (3.50, 3.50, 0.9),
              (7.20, 3.50, 0.45), (3.50, 7.20, 0.45),
              (1.25, 1.25, 0.8), (1.25, 3.50, 0.30), (3.50, 1.25, 0.30)]


def write_bruker_2d(root, si2=512, si1=256, xdim2=128, xdim1=64,
                    sf=500.13, sw=5000.0, offset=10.0, peaks=COSY_PEAKS):
    """A processed 2D COSY-like dataset, written in Bruker submatrix order.

    Returns ``(path, expected_matrix)`` so a test can check the de-tiling is
    exact rather than merely plausible.
    """
    os.makedirs(os.path.join(root, "pdata", "1"), exist_ok=True)
    step2 = (sw / sf) / si2
    step1 = (sw / sf) / si1

    def col(ppm):
        return int(round((offset - ppm) / step2))

    def row(ppm):
        return int(round((offset - ppm) / step1))

    matrix = [[0.0] * si2 for _ in range(si1)]
    for f2, f1, height in peaks:
        c0, r0 = col(f2), row(f1)
        for r in range(max(0, r0 - 10), min(si1, r0 + 11)):
            for c in range(max(0, c0 - 10), min(si2, c0 + 11)):
                matrix[r][c] += height * math.exp(
                    -(((c - c0) / 2.5) ** 2 + ((r - r0) / 2.5) ** 2))
    expected = [[float(int(v * 1e7)) for v in row_] for row_ in matrix]

    flat = []
    for tile_row in range(si1 // xdim1):
        for tile_col in range(si2 // xdim2):
            for r in range(xdim1):
                source = expected[tile_row * xdim1 + r]
                flat.extend(int(v) for v in
                            source[tile_col * xdim2:(tile_col + 1) * xdim2])
    with open(os.path.join(root, "pdata", "1", "2rr"), "wb") as fh:
        fh.write(struct.pack("<%di" % len(flat), *flat))

    _write_params(os.path.join(root, "acqus"),
                  [("SFO1", sf), ("SW_h", sw), ("NUC1", "<1H>"),
                   ("PULPROG", "<cosygpqf>"), ("NS", 4), ("SOLVENT", "<CDCl3>"),
                   ("TD", 1024), ("BYTORDA", 0), ("DTYPA", 2)])
    _write_params(os.path.join(root, "acqu2s"),
                  [("SFO1", sf), ("SW_h", sw), ("NUC1", "<1H>"), ("TD", 128)])
    _write_params(os.path.join(root, "pdata", "1", "procs"),
                  [("SI", si2), ("SF", sf), ("SW_p", sw), ("OFFSET", offset),
                   ("NC_proc", 0), ("XDIM", xdim2), ("BYTORDP", 0), ("DTYPP", 0)])
    _write_params(os.path.join(root, "pdata", "1", "proc2s"),
                  [("SI", si1), ("SF", sf), ("SW_p", sw), ("OFFSET", offset),
                   ("NC_proc", 0), ("XDIM", xdim1), ("BYTORDP", 0), ("DTYPP", 0)])
    with open(os.path.join(root, "pdata", "1", "title"), "w") as fh:
        fh.write("Synthetic COSY\n")
    return root, expected


# ---------------------------------------------------------------------------
# JCAMP-DX
# ---------------------------------------------------------------------------


JCAMP_HEADER = [
    "##TITLE=synthetic", "##JCAMP-DX=5.01", "##DATA TYPE=NMR SPECTRUM",
    "##XUNITS=PPM", "##YUNITS=ARBITRARY UNITS", "##.OBSERVE FREQUENCY=400.0",
    "##.OBSERVE NUCLEUS=^1H", "##.SOLVENT NAME=CDCl3",
]


def jcamp_file(path, body_lines, firstx=10.0, lastx=3.0, npoints=8,
               yfactor=1.0):
    lines = list(JCAMP_HEADER) + [
        "##FIRSTX=%s" % firstx, "##LASTX=%s" % lastx,
        "##YFACTOR=%s" % yfactor, "##NPOINTS=%d" % npoints,
        "##XYDATA=(X++(Y..Y))",
    ] + list(body_lines) + ["##END="]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    return path


# ---------------------------------------------------------------------------
# ACD .esp
# ---------------------------------------------------------------------------


def write_esp(path, peaks=((7.26, 1.0), (2.10, 0.4)), npts=4096,
              sf=500.203, sw=9999.85, centre_hz=2946.93):
    """A minimal ACD ``.esp`` file: header tags then real and imaginary blocks.

    Data is stored low-to-high ppm, which is the reverse of Bruker; the reader
    has to flip it.
    """
    sw_ppm = sw / sf
    right = centre_hz / sf - sw_ppm / 2.0

    ascending = []
    for i in range(npts):
        ppm = right + i * sw_ppm / npts
        value = 0.0
        for centre, height in peaks:
            u = 2.0 * (ppm - centre) * sf / 3.0
            value += height / (1.0 + u * u)
        ascending.append(value * 1e9)

    def tag(number, payload):
        return bytes([number, len(payload)]) + payload

    header = b""
    header += tag(0x03, struct.pack("<f", sw))
    header += tag(0x04, struct.pack("<f", centre_hz))
    header += tag(0x05, struct.pack("<f", sf))
    header += tag(0x10, struct.pack("<i", npts // 2))
    header += tag(0x0A, b"synthetic esp")
    header += tag(0x11, b"TRIFLUOROACETIC ACID-d")
    header += tag(0x19, b"zg30")
    header += tag(0x0F, struct.pack("<h", 16))
    header += tag(0x0E, struct.pack("<f", 25.0))
    header += b"\x00\x00"

    with open(path, "wb") as fh:
        fh.write(struct.pack("<H", 12) + b"(C) ACD 1994")
        fh.write(struct.pack("<H", 14) + b".ESP.( V 1.0 )")
        fh.write(header)
        fh.write(struct.pack("<%df" % npts, *ascending))
        fh.write(struct.pack("<%df" % npts, *([0.0] * npts)))
    return path


def write_bruker_ser(root, peaks=((7.2, 3.5), (3.5, 7.2), (2.0, 2.0)),
                     td2=512, td1=128, sf=500.0, sw=5000.0, o1=None,
                     fnmode=4, decay=0.06):
    """Write a raw 2D ``ser`` with States detection and known cross peaks.

    ``peaks`` are (F2 ppm, F1 ppm).  Rows alternate cosine- and
    sine-modulated, which is what States detection produces, and each row is
    padded out to a 1024-byte boundary the way TopSpin writes them.
    """
    os.makedirs(root, exist_ok=True)
    if o1 is None:
        o1 = 5.0 * sf                      # carrier at 5 ppm
    dwell2 = 1.0 / sw
    dwell1 = 1.0 / sw
    n_complex = td2 // 2
    pairs = td1 // 2

    def offset_hz(ppm):
        return ppm * sf - o1               # relative to the carrier

    rows = []
    for k in range(pairs):
        t1 = k * dwell1
        for component in ("cos", "sin"):
            row = []
            for n in range(n_complex):
                t2 = n * dwell2
                total = 0j
                for f2_ppm, f1_ppm in peaks:
                    w1 = 2.0 * math.pi * offset_hz(f1_ppm)
                    w2 = 2.0 * math.pi * offset_hz(f2_ppm)
                    modulation = (math.cos(w1 * t1) if component == "cos"
                                  else math.sin(w1 * t1))
                    total += (modulation
                              * cmath.exp(1j * w2 * t2)
                              * math.exp(-t1 / decay) * math.exp(-t2 / decay))
                row.append(total * 1e5)
            rows.append(row)

    width = 4                              # int32
    row_bytes = ((td2 * width + 1023) // 1024) * 1024
    padding = row_bytes - td2 * width
    with open(os.path.join(root, "ser"), "wb") as fh:
        for row in rows:
            flat = []
            for value in row:
                flat.append(int(value.real))
                flat.append(int(value.imag))
            fh.write(struct.pack("<%di" % len(flat), *flat))
            fh.write(b"\x00" * padding)

    _write_params(os.path.join(root, "acqus"),
                  [("SFO1", sf), ("SW_h", sw), ("O1", o1), ("NUC1", "<1H>"),
                   ("PULPROG", "<cosygpqf>"), ("NS", 2), ("SOLVENT", "<CDCl3>"),
                   ("TD", td2), ("BYTORDA", 0), ("DTYPA", 0), ("GRPDLY", 0)])
    _write_params(os.path.join(root, "acqu2s"),
                  [("SFO1", sf), ("SW_h", sw), ("O1", o1), ("NUC1", "<1H>"),
                   ("TD", td1), ("FnMODE", fnmode)])
    return root, list(peaks)
