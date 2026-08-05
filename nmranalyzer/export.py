"""Vector plot export and session persistence.

SVG rather than PNG: tkinter can only rasterise text through a font engine it
does not expose, so a hand-rolled PNG would lose every label.  SVG keeps the
labels as real text, scales without blurring, and drops straight into Word,
PowerPoint and Illustrator.
"""

from __future__ import annotations

import json
import os

from . import analysis
from . import nmrio

SESSION_FORMAT = "nmr-analyzer-session"
SESSION_VERSION = 1


# ---------------------------------------------------------------------------
# SVG
# ---------------------------------------------------------------------------

_ESCAPE = {"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;"}


def _esc(text):
    return "".join(_ESCAPE.get(ch, ch) for ch in str(text))


def _envelope(spec, left, right, x0, x1, base_y, scale, top, bottom):
    """Screen coordinates for a trace, one min/max pair per output column."""
    a = spec.clamp(spec.index(left))
    b = spec.clamp(spec.index(right))
    if a > b:
        a, b = b, a
    if b <= a:
        return []

    span = left - right
    width = x1 - x0
    columns = max(1, int(width))
    count = b - a + 1
    data = spec.real
    points = []

    def sx(ppm):
        return x0 + (left - ppm) / span * width

    def sy(value):
        return max(top - 400.0, min(bottom + 400.0, base_y - value * scale))

    if count <= columns * 2:
        for i in range(a, b + 1):
            points.append((sx(spec.ppm(i)), sy(data[i])))
    else:
        per = count / columns
        idx = a
        for col in range(columns):
            stop = min(a + int((col + 1) * per), b + 1)
            if stop <= idx:
                stop = idx + 1
            chunk = data[idx:stop]
            if not chunk:
                break
            x = x0 + col
            points.append((x, sy(max(chunk))))
            points.append((x, sy(min(chunk))))
            idx = stop
    return points


def _nice_step(span):
    import math
    if span <= 0:
        return 1.0
    raw = span / 10.0
    mag = 10.0 ** math.floor(math.log10(raw))
    for mult in (1.0, 2.0, 5.0, 10.0):
        if raw <= mult * mag:
            return mult * mag
    return 10.0 * mag


def write_svg(path, spectra, left, right, width=1000, height=520,
              stack=False, normalise_each=True, show_peaks=True,
              show_integrals=True, active=None, title=""):
    """Write the current plot as a standalone SVG file."""
    visible = [s for s in spectra if s.visible]
    if not visible or left <= right:
        raise ValueError("nothing to export")

    margin_l, margin_r, margin_t, axis_h = 14, 14, 22, 40
    x0, x1 = margin_l, width - margin_r
    plot_top = margin_t + (16 if title else 0)
    plot_bottom = height - axis_h

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
        'viewBox="0 0 %d %d" font-family="Helvetica, Arial, sans-serif">'
        % (width, height, width, height),
        '<rect width="%d" height="%d" fill="white"/>' % (width, height),
    ]
    if title:
        parts.append('<text x="%d" y="%d" font-size="13" fill="#333">%s</text>'
                     % (x0, margin_t, _esc(title)))

    lanes = len(visible) if stack else 1
    lane_h = (plot_bottom - plot_top) / lanes

    def trace_max(spec):
        a = spec.clamp(spec.index(left))
        b = spec.clamp(spec.index(right))
        if a > b:
            a, b = b, a
        data = spec.real
        if b - a < 4:
            chunk = data[a:b + 1]
            return max(chunk) if chunk else 1.0
        best = 0.0
        for i in range(a + 1, b):
            smoothed = (data[i - 1] + data[i] + data[i + 1]) / 3.0
            if smoothed > best:
                best = smoothed
        return best or 1.0

    global_max = max((trace_max(s) for s in visible), default=1.0)

    # Clip so a tall line cannot spill over the axis or the title.
    parts.append('<clipPath id="plot"><rect x="%d" y="%d" width="%d" '
                 'height="%d"/></clipPath>'
                 % (x0, plot_top, x1 - x0, plot_bottom - plot_top))

    for i, spec in enumerate(visible):
        base_y = plot_bottom - (i * lane_h if stack else 0)
        usable = lane_h * 0.92 if stack else (plot_bottom - plot_top) * 0.95
        ref = trace_max(spec) if normalise_each else global_max
        scale = usable * spec.scale / (ref or 1.0)

        if show_integrals and spec.regions:
            norm = analysis.normalise(spec.regions)
            for region, value in zip(spec.regions, norm):
                rx1 = x0 + (left - region.hi) / (left - right) * (x1 - x0)
                rx2 = x0 + (left - region.lo) / (left - right) * (x1 - x0)
                if rx2 < x0 or rx1 > x1:
                    continue
                parts.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" '
                             'fill="#eef4fb"/>'
                             % (rx1, plot_top, max(rx2 - rx1, 0.8),
                                base_y - plot_top))
                label = "%gH" % region.protons if region.protons else "%.2f" % value
                parts.append('<text x="%.2f" y="%.2f" font-size="10" fill="#2266aa" '
                             'text-anchor="middle">%s</text>'
                             % ((rx1 + rx2) / 2, plot_top + 10, _esc(label)))

        points = _envelope(spec, left, right, x0, x1, base_y, scale,
                           plot_top, plot_bottom)
        if len(points) >= 2:
            d = " ".join("%.2f,%.2f" % p for p in points)
            parts.append('<polyline points="%s" fill="none" stroke="%s" '
                         'stroke-width="%s" clip-path="url(#plot)"/>'
                         % (d, spec.color, "1.2" if spec is active else "0.9"))

        if show_peaks and spec is active and spec.peaks:
            last_x = None
            for peak in spec.peaks:
                if not (right <= peak.ppm <= left):
                    continue
                px = x0 + (left - peak.ppm) / (left - right) * (x1 - x0)
                py = max(plot_top, base_y - peak.height * scale)
                parts.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" '
                             'stroke="#aa3333" stroke-width="0.8"/>'
                             % (px, py - 4, px, py - 11))
                if last_x is None or abs(px - last_x) > 34:
                    parts.append('<text x="%.2f" y="%.2f" font-size="9" '
                                 'fill="#aa3333" text-anchor="middle">%.3f</text>'
                                 % (px, py - 13, peak.ppm))
                    last_x = px

        if stack:
            parts.append('<text x="%.2f" y="%.2f" font-size="10" fill="%s" '
                         'text-anchor="end">%s</text>'
                         % (x1 - 4, base_y - usable + 8, spec.color,
                            _esc(spec.name)))

    # Axis
    parts.append('<line x1="%d" y1="%.2f" x2="%d" y2="%.2f" stroke="#333" '
                 'stroke-width="1"/>' % (x0, plot_bottom, x1, plot_bottom))
    step = _nice_step(left - right)
    value = int(right / step) * step
    while value <= left + step:
        if right - 1e-9 <= value <= left + 1e-9:
            tx = x0 + (left - value) / (left - right) * (x1 - x0)
            parts.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" '
                         'stroke="#333" stroke-width="1"/>'
                         % (tx, plot_bottom, tx, plot_bottom + 5))
            parts.append('<text x="%.2f" y="%.2f" font-size="10" fill="#333" '
                         'text-anchor="middle">%g</text>'
                         % (tx, plot_bottom + 17, round(value, 4)))
        value += step
    parts.append('<text x="%d" y="%.2f" font-size="10" fill="#666" '
                 'text-anchor="end">ppm</text>' % (x1, plot_bottom + 32))
    parts.append("</svg>")

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(parts))
    return path


# ---------------------------------------------------------------------------
# Sessions
# ---------------------------------------------------------------------------


def save_session(path, spectra, view, options):
    """Write the analysis state (not the data) as JSON."""
    entries = []
    for spec in spectra:
        entries.append({
            "source": spec.source,
            "source_index": getattr(spec, "source_index", 0),
            "name": spec.name,
            "color": spec.color,
            "visible": spec.visible,
            "scale": spec.scale,
            "ref_shift": spec.ref_shift,
            "p0": spec.p0,
            "p1": spec.p1,
            "lb": spec.lb,
            "si": spec.si,
            "baseline_on": spec.baseline_on,
            "baseline_method": spec.baseline_method,
            "regions": [{"lo": r.lo, "hi": r.hi, "protons": r.protons,
                         "label": r.label} for r in spec.regions],
        })
    payload = {
        "format": SESSION_FORMAT,
        "version": SESSION_VERSION,
        "view": {"left": view[0], "right": view[1]},
        "options": dict(options),
        "spectra": entries,
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
    return path


def load_session(path):
    """Re-open the sources named in a session and restore the analysis state.

    Returns ``(spectra, view, options, warnings)``.  Sources that have moved or
    been deleted are reported in ``warnings`` rather than aborting the load.
    """
    with open(path, "r", encoding="utf-8") as fh:
        payload = json.load(fh)
    if payload.get("format") != SESSION_FORMAT:
        raise ValueError("not an NMR Analyzer session file")

    spectra = []
    warnings = []
    cache = {}
    for entry in payload.get("spectra", []):
        source = entry.get("source", "")
        if source not in cache:
            if not os.path.exists(source):
                warnings.append("missing: %s" % source)
                cache[source] = []
            else:
                try:
                    cache[source] = nmrio.load(source)
                except Exception as exc:
                    warnings.append("%s: %s" % (os.path.basename(source), exc))
                    cache[source] = []
        loaded = cache[source]
        index = entry.get("source_index", 0)
        if index >= len(loaded):
            continue
        spec = loaded[index]

        spec.name = entry.get("name", spec.name)
        spec.color = entry.get("color", spec.color)
        spec.visible = entry.get("visible", True)
        spec.scale = entry.get("scale", 1.0)
        spec.ref_shift = entry.get("ref_shift", 0.0)
        spec.p0 = entry.get("p0", 0.0)
        spec.p1 = entry.get("p1", 0.0)
        spec.lb = entry.get("lb", 0.0)
        spec.si = entry.get("si", spec.npoints)
        spec.baseline_on = entry.get("baseline_on", False)
        spec.baseline_method = entry.get("baseline_method", "spline")
        spec.source_index = index
        if spec.p0 or spec.p1 or spec.baseline_on or (spec.fid and spec.lb):
            spec.reprocess()

        spec.regions = []
        for r in entry.get("regions", []):
            region = analysis.Region(r["lo"], r["hi"], label=r.get("label", ""))
            region.protons = r.get("protons")
            analysis.integrate_region(spec, region)
            spec.regions.append(region)
        spectra.append(spec)

    view = payload.get("view", {})
    options = payload.get("options", {})
    return spectra, (view.get("left", 12.0), view.get("right", -1.0)), options, warnings
