#!/usr/bin/env python3
"""Draw the application icon and write it out for every platform.

The artwork is a spectrum, and it is a real one: the trace is a sum of
Lorentzian lines evaluated with the same shape function the analyser fits,
rather than a hand-drawn squiggle.

Standard library only.  PNG is written directly (zlib plus four chunks), the
Windows ``.ico`` is assembled from PNG frames, and the macOS ``.icns`` is
produced with the system ``iconutil``.

    python3 packaging/make_icon.py
"""

from __future__ import annotations

import math
import os
import struct
import subprocess
import sys
import zlib

SIZE = 1024                 # master resolution
CORNER = 0.22               # corner radius as a fraction of the side

# Background: a slight vertical gradient reads better than a flat fill at
# large sizes and costs nothing at small ones.
TOP_RGB = (32, 62, 118)
BOTTOM_RGB = (16, 30, 62)
TRACE_RGB = (255, 255, 255)
BASELINE_RGB = (120, 160, 220)

# (centre, height, half-width) in fractions of the plot width/height.
# Chosen to read as a spectrum at a glance: a doublet, a tall singlet, a
# second doublet.  Anything busier turns to mush at 32 px.
PEAKS = [
    (0.20, 0.44, 0.0105),
    (0.29, 0.44, 0.0105),
    (0.51, 0.95, 0.0115),
    (0.72, 0.35, 0.0105),
    (0.81, 0.35, 0.0105),
]


def lorentzian(x, centre, height, half_width):
    u = (x - centre) / half_width
    return height / (1.0 + u * u)


def spectrum(x):
    return sum(lorentzian(x, c, h, w) for c, h, w in PEAKS)


# ---------------------------------------------------------------------------
# Canvas
# ---------------------------------------------------------------------------


class Canvas:
    """A tiny RGBA canvas with the two primitives this icon needs."""

    def __init__(self, size):
        self.size = size
        self.buf = bytearray(size * size * 4)

    def _blend(self, x, y, rgb, alpha):
        if alpha <= 0.0 or not (0 <= x < self.size and 0 <= y < self.size):
            return
        if alpha > 1.0:
            alpha = 1.0
        i = (y * self.size + x) * 4
        buf = self.buf
        dst_a = buf[i + 3] / 255.0
        out_a = alpha + dst_a * (1.0 - alpha)
        if out_a <= 0.0:
            return
        for k in range(3):
            src = rgb[k] / 255.0
            dst = buf[i + k] / 255.0
            buf[i + k] = int(round(((src * alpha + dst * dst_a * (1.0 - alpha))
                                    / out_a) * 255.0))
        buf[i + 3] = int(round(out_a * 255.0))

    def rounded_rect(self, radius, top_rgb, bottom_rgb):
        """Fill a rounded square with a vertical gradient, edges antialiased."""
        n = self.size
        for y in range(n):
            t = y / (n - 1)
            rgb = tuple(int(round(top_rgb[k] + (bottom_rgb[k] - top_rgb[k]) * t))
                        for k in range(3))
            for x in range(n):
                # Distance outside the rounded rectangle, in pixels.
                dx = max(radius - x, x - (n - 1 - radius), 0.0)
                dy = max(radius - y, y - (n - 1 - radius), 0.0)
                if dx > 0.0 and dy > 0.0:
                    d = math.hypot(dx, dy) - radius
                else:
                    d = max(dx, dy) - radius
                coverage = 0.5 - d          # 1 px wide transition band
                if coverage >= 1.0:
                    self._blend(x, y, rgb, 1.0)
                elif coverage > 0.0:
                    self._blend(x, y, rgb, coverage)

    def stamp(self, cx, cy, radius, rgb, alpha=1.0):
        """An antialiased filled disc, used to give the trace its thickness."""
        x0, x1 = int(math.floor(cx - radius - 1)), int(math.ceil(cx + radius + 1))
        y0, y1 = int(math.floor(cy - radius - 1)), int(math.ceil(cy + radius + 1))
        for y in range(y0, y1 + 1):
            for x in range(x0, x1 + 1):
                d = math.hypot(x + 0.5 - cx, y + 0.5 - cy)
                coverage = radius - d + 0.5
                if coverage >= 1.0:
                    self._blend(x, y, rgb, alpha)
                elif coverage > 0.0:
                    self._blend(x, y, rgb, coverage * alpha)

    def stroke(self, points, width, rgb, alpha=1.0):
        radius = width / 2.0
        for x, y in points:
            self.stamp(x, y, radius, rgb, alpha)

    def rows(self):
        n = self.size
        return [bytes(self.buf[y * n * 4:(y + 1) * n * 4]) for y in range(n)]


# ---------------------------------------------------------------------------
# PNG
# ---------------------------------------------------------------------------


def png_bytes(width, height, rows):
    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF))

    raw = b"".join(b"\x00" + row for row in rows)          # filter type 0
    header = struct.pack(">IIBBBBB", width, height, 8, 6, 0, 0, 0)
    return (b"\x89PNG\r\n\x1a\n"
            + chunk(b"IHDR", header)
            + chunk(b"IDAT", zlib.compress(raw, 9))
            + chunk(b"IEND", b""))


def write_png(path, canvas):
    with open(path, "wb") as fh:
        fh.write(png_bytes(canvas.size, canvas.size, canvas.rows()))


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------


def draw_icon(size=SIZE):
    canvas = Canvas(size)
    canvas.rounded_rect(size * CORNER, TOP_RGB, BOTTOM_RGB)

    margin = size * 0.16
    plot_w = size - 2 * margin
    baseline = size * 0.78
    peak_room = baseline - size * 0.20

    # Baseline, kept faint so it reads as a spectrum without competing.
    canvas.stroke([(margin + i * 0.5, baseline)
                   for i in range(int(plot_w * 2) + 1)],
                  size * 0.012, BASELINE_RGB, alpha=0.55)

    # The trace, sampled finely enough that consecutive stamps overlap.
    steps = int(plot_w * 3)
    points = []
    for i in range(steps + 1):
        fx = i / steps
        points.append((margin + fx * plot_w,
                       baseline - spectrum(fx) * peak_room))
    canvas.stroke(points, size * 0.036, TRACE_RGB)
    return canvas


# ---------------------------------------------------------------------------
# Platform packaging
# ---------------------------------------------------------------------------


def resize(src, dst, size):
    subprocess.run(["sips", "-z", str(size), str(size), src, "--out", dst],
                   check=True, capture_output=True)


def build_ico(master_png, out_path, sizes=(16, 24, 32, 48, 64, 128, 256)):
    """Assemble a Windows .ico from PNG frames.

    Vista and later accept PNG-compressed entries, which keeps this to a
    header plus the frames themselves.
    """
    frames = []
    work = os.path.dirname(out_path)
    for size in sizes:
        tmp = os.path.join(work, "_ico_%d.png" % size)
        resize(master_png, tmp, size)
        with open(tmp, "rb") as fh:
            frames.append((size, fh.read()))
        os.remove(tmp)

    header = struct.pack("<HHH", 0, 1, len(frames))
    offset = len(header) + 16 * len(frames)
    entries, payload = b"", b""
    for size, data in frames:
        entries += struct.pack("<BBBBHHII",
                               0 if size >= 256 else size,
                               0 if size >= 256 else size,
                               0, 0, 1, 32, len(data), offset)
        payload += data
        offset += len(data)
    with open(out_path, "wb") as fh:
        fh.write(header + entries + payload)


def build_icns(master_png, out_path):
    """Build a macOS .icns via the system iconutil."""
    work = os.path.splitext(out_path)[0] + ".iconset"
    os.makedirs(work, exist_ok=True)
    for size in (16, 32, 128, 256, 512):
        resize(master_png, os.path.join(work, "icon_%dx%d.png" % (size, size)), size)
        resize(master_png, os.path.join(work, "icon_%dx%d@2x.png" % (size, size)),
               size * 2)
    subprocess.run(["iconutil", "-c", "icns", work, "-o", out_path],
                   check=True, capture_output=True)
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    master = os.path.join(here, "icon.png")

    print("drawing %dx%d ..." % (SIZE, SIZE))
    write_png(master, draw_icon(SIZE))
    print("  %s  %.0f KB" % (os.path.basename(master),
                             os.path.getsize(master) / 1024))

    if sys.platform == "darwin":
        ico = os.path.join(here, "icon.ico")
        build_ico(master, ico)
        print("  icon.ico   %.0f KB" % (os.path.getsize(ico) / 1024))
        icns = os.path.join(here, "icon.icns")
        build_icns(master, icns)
        print("  icon.icns  %.0f KB" % (os.path.getsize(icns) / 1024))
    else:
        print("  (.ico/.icns need macOS tools; run this on a Mac)")


if __name__ == "__main__":
    main()
