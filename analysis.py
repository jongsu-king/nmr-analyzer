"""Peak picking, integration and multiplet interpretation."""

from __future__ import annotations

import math

import dsp

MULTIPLICITY_NAMES = {
    1: ("s", "singlet"),
    2: ("d", "doublet"),
    3: ("t", "triplet"),
    4: ("q", "quartet"),
    5: ("quint", "quintet"),
    6: ("sext", "sextet"),
    7: ("sept", "septet"),
}


class Peak:
    def __init__(self, index, ppm, height, hz=0.0):
        self.index = index
        self.ppm = ppm
        self.height = height
        self.hz = hz

    def __repr__(self):
        return "Peak(%.4f ppm, %.3g)" % (self.ppm, self.height)


class Region:
    """An integration region, in ppm."""

    def __init__(self, lo, hi, value=0.0, label=""):
        self.lo = min(lo, hi)
        self.hi = max(lo, hi)
        self.value = value
        self.label = label
        self.protons = None     # user-assigned H count, drives normalisation
        self.peaks = []
        self.multiplet = None
        self.fit = None         # fitting.FitResult once deconvoluted

    @property
    def center(self):
        return (self.lo + self.hi) / 2.0


# ---------------------------------------------------------------------------
# Peak picking
# ---------------------------------------------------------------------------


def _prominence(data, i, span):
    """How far a local maximum stands above the higher of its two valleys.

    Ripples riding on the flank of a broad line have a large absolute height
    but almost no prominence, so this is what separates real multiplet lines
    from digitisation noise.
    """
    top = data[i]
    left_min = top
    j = i - 1
    limit = max(0, i - span)
    while j >= limit:
        if data[j] > top:
            break
        left_min = min(left_min, data[j])
        j -= 1
    right_min = top
    j = i + 1
    limit = min(len(data) - 1, i + span)
    while j <= limit:
        if data[j] > top:
            break
        right_min = min(right_min, data[j])
        j += 1
    return top - max(left_min, right_min)


def pick_peaks(spec, lo_ppm=None, hi_ppm=None, sensitivity=8.0,
               min_sep_hz=1.0, max_peaks=400):
    """Find local maxima that rise ``sensitivity`` x noise above the baseline.

    A candidate must also be *prominent* by the same margin, which suppresses
    the ripples that sit on top of broad lines.  The apex position is refined
    by fitting a parabola through the maximum and its two neighbours, which
    recovers sub-point accuracy on sharp lines.
    """
    data = spec.real
    n = len(data)
    if n < 3:
        return []

    lo = 0 if hi_ppm is None else spec.clamp(spec.index(hi_ppm))
    hi = n - 1 if lo_ppm is None else spec.clamp(spec.index(lo_ppm))
    if lo > hi:
        lo, hi = hi, lo

    noise = dsp.noise_level(data)
    threshold = noise * sensitivity
    if threshold <= 0.0:
        threshold = max(data[lo:hi + 1] or [0.0]) * 0.01

    per_hz = spec.hz_per_point() or 1.0
    min_sep = max(1, int(round(min_sep_hz / per_hz)))
    prom_span = max(min_sep * 4, int(round(30.0 / per_hz)))

    candidates = []
    for i in range(max(lo, 1), min(hi, n - 2) + 1):
        y = data[i]
        if y < threshold:
            continue
        if y >= data[i - 1] and y > data[i + 1]:
            if _prominence(data, i, prom_span) >= threshold:
                candidates.append(i)

    # Keep the tallest peak within each min_sep window.
    candidates.sort(key=lambda i: -data[i])
    kept = []
    for i in candidates:
        if all(abs(i - j) >= min_sep for j in kept):
            kept.append(i)
        if len(kept) >= max_peaks:
            break
    kept.sort()

    peaks = []
    for i in kept:
        y0, y1, y2 = data[i - 1], data[i], data[i + 1]
        denom = y0 - 2.0 * y1 + y2
        shift = 0.5 * (y0 - y2) / denom if denom else 0.0
        shift = max(-0.5, min(0.5, shift))
        ppm = spec.ppm(0) + (i + shift) * spec.delta_ppm
        peaks.append(Peak(i, ppm, y1, ppm * spec.sf))
    return peaks


# ---------------------------------------------------------------------------
# Integration
# ---------------------------------------------------------------------------


def integrate(spec, lo_ppm, hi_ppm):
    """Area under the spectrum between two ppm bounds."""
    a = spec.clamp(spec.index(hi_ppm))
    b = spec.clamp(spec.index(lo_ppm))
    if a > b:
        a, b = b, a
    step = abs(spec.delta_ppm)
    return sum(spec.real[a:b + 1]) * step


def refresh_peaks(spec):
    """Re-read peak heights from the current data.

    Peaks remember the height they had when they were picked, so after
    re-phasing, re-transforming or a baseline correction their markers would
    otherwise float at heights the spectrum no longer has.
    """
    for peak in spec.peaks:
        index = spec.clamp(spec.index(peak.ppm))
        peak.index = index
        peak.height = spec.real[index]


def integrate_region(spec, region):
    region.value = integrate(spec, region.lo, region.hi)
    region.peaks = pick_peaks(spec, region.lo, region.hi)
    region.multiplet = analyse_multiplet(region.peaks, spec.sf)
    return region


def normalise(regions):
    """Scale every integral so that assigned regions match their H counts.

    The reference factor is the average of ``value / protons`` over all
    regions the user has annotated; if none are annotated the smallest
    integral is taken as 1 H.
    """
    assigned = [r for r in regions if r.protons]
    if assigned:
        factors = [r.value / r.protons for r in assigned if r.protons]
        factor = sum(factors) / len(factors)
    else:
        positives = [r.value for r in regions if r.value > 0]
        factor = min(positives) if positives else 1.0
    if factor == 0.0:
        factor = 1.0
    return [r.value / factor for r in regions]


# ---------------------------------------------------------------------------
# Multiplet interpretation
# ---------------------------------------------------------------------------


def _binomial(n):
    row = [1]
    for k in range(n - 1):
        row = [1] + [row[i] + row[i + 1] for i in range(len(row) - 1)] + [1]
    return row


def _looks_binomial(heights, tol=0.40):
    n = len(heights)
    expected = _binomial(n)
    scale = max(heights) / max(expected)
    for h, e in zip(heights, expected):
        target = e * scale
        if target == 0:
            continue
        if abs(h - target) / target > tol:
            return False
    return True


class Multiplet:
    def __init__(self, center_ppm, pattern, name, couplings, width_hz, n_lines):
        self.center_ppm = center_ppm
        self.pattern = pattern          # "d", "dd", "m", ...
        self.name = name                # human readable
        self.couplings = couplings      # list of J values in Hz
        self.width_hz = width_hz
        self.n_lines = n_lines

    def describe(self):
        if self.couplings:
            js = ", ".join("J = %.1f Hz" % j for j in self.couplings)
            return "%s (%s)" % (self.pattern, js)
        return self.pattern


def analyse_multiplet(peaks, sf, tol_hz=0.7, tol_frac=0.12):
    """Classify a group of lines as s/d/t/q/dd/... and extract J values.

    ``tol_hz`` / ``tol_frac`` set how close two line spacings must be before
    they are treated as the same coupling constant.
    """
    if not peaks:
        return None

    peaks = sorted(peaks, key=lambda p: -p.ppm)
    hz = [p.ppm * sf for p in peaks]
    heights = [p.height for p in peaks]
    weight = sum(heights) or 1.0
    center = sum(p.ppm * p.height for p in peaks) / weight
    width = abs(hz[0] - hz[-1]) if len(hz) > 1 else 0.0
    n = len(peaks)

    if n == 1:
        return Multiplet(center, "s", "singlet", [], 0.0, 1)

    spacings = [abs(hz[i] - hz[i + 1]) for i in range(n - 1)]

    def close(a, b):
        return abs(a - b) <= max(tol_hz, tol_frac * max(a, b))

    # First-order binomial multiplet: equal spacings and Pascal intensities.
    if n <= 7 and all(close(s, spacings[0]) for s in spacings) \
            and _looks_binomial(heights):
        pattern, name = MULTIPLICITY_NAMES[n]
        j = sum(spacings) / len(spacings)
        return Multiplet(center, pattern, name, [j], width, n)

    # Doublet of doublets: spacings J2, J1-J2, J2 with four equal lines.
    if n == 4 and close(spacings[0], spacings[2]):
        j2 = (spacings[0] + spacings[2]) / 2.0
        j1 = spacings[1] + j2
        if j1 > j2 and _roughly_equal(heights, 0.45):
            return Multiplet(center, "dd", "doublet of doublets",
                             [max(j1, j2), min(j1, j2)], width, n)

    # Doublet of triplets / triplet of doublets: six lines, two couplings.
    if n == 6 and close(spacings[0], spacings[1]) and close(spacings[3], spacings[4]):
        small = (spacings[0] + spacings[1] + spacings[3] + spacings[4]) / 4.0
        large = spacings[2] + small
        return Multiplet(center, "dt", "doublet of triplets",
                         [max(large, small), min(large, small)], width, n)

    return Multiplet(center, "m", "multiplet", [], width, n)


def _roughly_equal(values, tol):
    top = max(values)
    if top == 0:
        return False
    return all(abs(v - top) / top <= tol for v in values)


# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------


class Component:
    """One species in a mixture, defined by a region and its proton count."""

    def __init__(self, label, region, protons):
        self.label = label
        self.region = region
        self.protons = protons

    @property
    def moles(self):
        """Integral per proton — proportional to the amount of substance."""
        if not self.protons:
            return 0.0
        return self.region.value / self.protons


def composition(components):
    """Mole fractions and conversion for a set of labelled components.

    Returns ``(rows, conversion)`` where each row is
    ``(label, moles, mole_fraction_percent)``.  ``conversion`` is the fraction
    of everything that is *not* the first component, i.e. treat the first
    component as starting material and the rest as products.
    """
    usable = [c for c in components if c.protons]
    total = sum(c.moles for c in usable)
    rows = []
    for c in usable:
        frac = 100.0 * c.moles / total if total else 0.0
        rows.append((c.label, c.moles, frac))
    if not usable or total <= 0:
        return rows, 0.0
    remaining = usable[0].moles / total
    return rows, 100.0 * (1.0 - remaining)


def format_composition(components):
    rows, conversion = composition(components)
    if not rows:
        return "Assign a proton count to at least two regions."
    width = max(len(r[0]) for r in rows)
    lines = ["%-*s  %12s  %9s" % (width, "Component", "Integral/H", "mol %")]
    for label, moles, frac in rows:
        lines.append("%-*s  %12.4g  %8.2f%%" % (width, label, moles, frac))
    lines.append("")
    lines.append("Treating %r as starting material:" % rows[0][0])
    lines.append("  conversion = %.2f%%" % conversion)
    if len(rows) > 2:
        products = rows[1:]
        product_total = sum(r[1] for r in products)
        if product_total:
            ratio = " : ".join("%.2f" % (r[1] / min(p[1] for p in products))
                               for r in products)
            lines.append("  product ratio (%s) = %s"
                         % (", ".join(r[0] for r in products), ratio))
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def format_report(spec, regions):
    """Render regions in the house style used for SI / paper experimentals.

    Proton counts are only written as ``nH`` where the user has actually
    assigned one; unassigned regions carry their relative integral instead, so
    the report never invents an integration the data does not support.
    """
    if not regions:
        return ""
    ordered = sorted(regions, key=lambda r: -r.center)
    norm = normalise(ordered)
    any_assigned = any(r.protons for r in ordered)
    nucleus = spec.meta.get("Nucleus", "1H") or "1H"
    freq = spec.meta.get("Frequency (MHz)", spec.sf)
    solvent = spec.meta.get("Solvent", "")
    header = "%s NMR (%.0f MHz, %s) delta " % (nucleus, float(freq), solvent)

    parts = []
    for region, value in zip(ordered, norm):
        m = region.multiplet
        if m and m.pattern != "m":
            shift = "%.2f" % m.center_ppm
        else:
            shift = "%.2f-%.2f" % (region.hi, region.lo)
        bits = [shift]
        if m and m.pattern != "m":
            bits.append(m.pattern)
            for j in m.couplings:
                bits.append("J = %.1f Hz" % j)
        elif m:
            bits.append("m")
        if region.protons:
            bits.append("%gH" % region.protons)
        elif any_assigned:
            bits.append("%.2fH" % value)
        else:
            bits.append("rel. %.2f" % value)
        parts.append("(%s)" % ", ".join(bits))
    return header + ", ".join(parts)
