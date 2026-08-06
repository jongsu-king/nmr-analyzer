"""Signal processing primitives for 1D NMR, standard library only.

Everything works on plain Python lists of ``complex`` (time domain) or
``float`` (frequency domain).  The routines are written so that a 64k point
transform stays under a couple of seconds on CPython.
"""

from __future__ import annotations

import cmath
import math

# ---------------------------------------------------------------------------
# FFT
# ---------------------------------------------------------------------------


def next_pow2(n: int) -> int:
    """Smallest power of two greater than or equal to ``n``."""
    p = 1
    while p < n:
        p <<= 1
    return p


def _bit_reverse(a: list) -> None:
    n = len(a)
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            a[i], a[j] = a[j], a[i]


def fft(data: list, inverse: bool = False) -> list:
    """In-place iterative radix-2 Cooley-Tukey FFT.

    ``data`` is padded with zeros up to the next power of two.  The result is
    a new list in standard DFT bin order (DC first).
    """
    n = next_pow2(len(data))
    a = [complex(v) for v in data]
    if len(a) < n:
        a.extend([0j] * (n - len(a)))
    _bit_reverse(a)

    length = 2
    while length <= n:
        # Twiddle factors for this stage, computed once and reused.
        ang = (2.0 if inverse else -2.0) * math.pi / length
        half = length >> 1
        tw = [cmath.exp(1j * ang * k) for k in range(half)]
        for start in range(0, n, length):
            mid = start + half
            for k in range(half):
                u = a[start + k]
                v = a[mid + k] * tw[k]
                a[start + k] = u + v
                a[mid + k] = u - v
        length <<= 1

    if inverse:
        inv = 1.0 / n
        a = [v * inv for v in a]
    return a


def fftshift(a: list) -> list:
    """Move the zero-frequency component to the centre of the spectrum."""
    h = len(a) // 2
    return a[h:] + a[:h]


def hilbert(real: list) -> list:
    """Reconstruct the imaginary counterpart of a real spectrum.

    Used when a dataset ships only the real part (some vendors drop ``1i``)
    but phase correction still has to be possible.
    """
    n = next_pow2(len(real))
    spec = fft(real)
    half = n // 2
    spec[0] = spec[0]
    for k in range(1, half):
        spec[k] *= 2.0
    for k in range(half + 1, n):
        spec[k] = 0j
    analytic = fft(spec, inverse=True)
    return [v.imag for v in analytic[: len(real)]]


# ---------------------------------------------------------------------------
# Time domain
# ---------------------------------------------------------------------------


def apodize(fid: list, dwell: float, lb: float = 0.0, gb: float = 0.0) -> list:
    """Apply exponential (``lb``, Hz) and/or Gaussian (``gb``, Hz) windows."""
    if lb == 0.0 and gb == 0.0:
        return list(fid)
    out = []
    for n, v in enumerate(fid):
        t = n * dwell
        w = 1.0
        if lb:
            w *= math.exp(-math.pi * lb * t)
        if gb:
            w *= math.exp(-((math.pi * gb * t) ** 2) / (4.0 * math.log(2.0)))
        out.append(v * w)
    return out


def zero_fill(fid: list, size: int) -> list:
    """Extend the FID with zeros up to ``size`` points."""
    if size <= len(fid):
        return list(fid[:size])
    return list(fid) + [0j] * (size - len(fid))


def group_delay_phase(spec: list, grpdly: float) -> list:
    """Undo the Bruker digital-filter group delay.

    Multiplying bin ``k`` by ``exp(2*pi*i*k*grpdly/N)`` is exactly a circular
    left shift of the FID by ``grpdly`` samples, and unlike an integer roll it
    handles the fractional part correctly.
    """
    if not grpdly:
        return spec
    n = len(spec)
    f = 2.0 * math.pi * grpdly / n
    return [v * cmath.exp(1j * f * k) for k, v in enumerate(spec)]


def transform(fid: list, dwell: float, si: int, lb: float = 0.0,
              gb: float = 0.0, grpdly: float = 0.0, fcor: float = 0.5) -> list:
    """FID -> complex spectrum ordered from high to low frequency.

    The output follows the Bruker convention: index ``i`` sits at frequency
    ``(N/2 - i) * SW / N``, so index 0 is the left-hand (high ppm) edge.
    """
    work = apodize(fid, dwell, lb, gb)
    work = zero_fill(work, si)
    if work:
        work[0] = work[0] * fcor  # scales the DC baseline offset
    spec = fft(work)
    spec = group_delay_phase(spec, grpdly)
    n = len(spec)
    half = n // 2
    return [spec[(half - i) % n] for i in range(n)]


# ---------------------------------------------------------------------------
# Phase
# ---------------------------------------------------------------------------


def phase(spec: list, p0: float, p1: float, pivot: float = 0.0) -> list:
    """Apply zero- and first-order phase correction, angles in degrees.

    ``pivot`` is the fractional position (0..1) across the spectrum at which
    the first-order correction is zero.
    """
    n = len(spec)
    if n == 0:
        return []
    r0 = math.radians(p0)
    r1 = math.radians(p1)
    out = []
    for k, v in enumerate(spec):
        ang = r0 + r1 * (k / n - pivot)
        out.append(v * cmath.exp(-1j * ang))
    return out


# A first-order correction is only justified when the data demands it, so a
# small penalty is charged for using one.  Without it, a spectrum whose peaks
# sit in a narrow window leaves p1 almost unconstrained: on a synthetic
# spectrum rotated by a known 40 degrees of pure zero-order phase, the
# unpenalised search returns a first-order term of 22 degrees, while with this
# penalty it returns zero and recovers the rotation exactly.  The value is
# small enough that a spectrum genuinely needing a large first-order
# correction still gets one.
P1_PENALTY = 5.0e-4


def _phase_penalty(spec: list, p0: float, p1: float, pivot: float,
                   step: int) -> float:
    """Negative-area objective used by :func:`autophase`.

    A correctly phased absorption spectrum is non-negative, so the summed
    magnitude of everything below the baseline is a good thing to minimise.
    """
    n = len(spec)
    r0 = math.radians(p0)
    r1 = math.radians(p1)
    neg = 0.0
    total = 0.0
    for k in range(0, n, step):
        v = spec[k]
        ang = r0 + r1 * (k / n - pivot)
        c = math.cos(ang)
        s = math.sin(ang)
        re = v.real * c + v.imag * s
        total += abs(re)
        if re < 0.0:
            neg -= re
    if total == 0.0:
        return 0.0
    return neg / total + P1_PENALTY * abs(p1)


def autophase(spec: list, pivot: float = 0.0, fit_p1: bool = True) -> tuple:
    """Grid search for (p0, p1) minimising the negative spectral area."""
    n = len(spec)
    step = max(1, n // 4096)

    best = (0.0, 0.0)
    best_val = float("inf")
    p1_range = [-180.0, -90.0, 0.0, 90.0, 180.0] if fit_p1 else [0.0]
    for p1 in p1_range:
        for p0 in range(-180, 180, 15):
            val = _phase_penalty(spec, p0, p1, pivot, step)
            if val < best_val:
                best_val = val
                best = (float(p0), p1)

    # Successive refinement around the coarse optimum.
    p0, p1 = best
    span0, span1 = 15.0, 90.0
    for _ in range(6):
        improved = False
        for d0 in (-span0, -span0 / 2, 0.0, span0 / 2, span0):
            for d1 in ((-span1, -span1 / 2, 0.0, span1 / 2, span1)
                       if fit_p1 else (0.0,)):
                val = _phase_penalty(spec, p0 + d0, p1 + d1, pivot, step)
                if val < best_val - 1e-9:
                    best_val = val
                    best = (p0 + d0, p1 + d1)
                    improved = True
        p0, p1 = best
        span0 /= 2.0
        span1 /= 2.0
        if not improved and span0 < 0.5:
            break
    return best


# ---------------------------------------------------------------------------
# Baseline
# ---------------------------------------------------------------------------


def solve_linear(matrix: list, rhs: list) -> list:
    """Solve ``matrix . x = rhs`` by Gaussian elimination with partial pivoting."""
    n = len(matrix)
    m = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-12:
            continue
        m[col], m[piv] = m[piv], m[col]
        inv = 1.0 / m[col][col]
        for r in range(col + 1, n):
            factor = m[r][col] * inv
            if factor:
                for c in range(col, n + 1):
                    m[r][c] -= factor * m[col][c]
    sol = [0.0] * n
    for r in range(n - 1, -1, -1):
        if abs(m[r][r]) < 1e-12:
            continue
        acc = m[r][n] - sum(m[r][c] * sol[c] for c in range(r + 1, n))
        sol[r] = acc / m[r][r]
    return sol


def baseline_nodes(data: list, segments: int = 64, quantile: float = 0.2) -> list:
    """Pick ``(index, level)`` anchors that most likely sit on the baseline.

    Each segment contributes the median of its lowest-``quantile`` fraction of
    points, which ignores peaks without needing a peak list.
    """
    n = len(data)
    if n == 0:
        return []
    seg = max(1, n // max(1, segments))
    nodes = []
    for start in range(0, n, seg):
        chunk = sorted(data[start:start + seg])
        if not chunk:
            continue
        take = max(1, int(len(chunk) * quantile))
        low = chunk[:take]
        nodes.append((start + seg // 2, low[len(low) // 2]))
    return nodes


def _median(values: list) -> float:
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def baseline_correct(data: list, order: int = 3, segments: int = 256,
                     method: str = "spline") -> list:
    """Remove the baseline.

    ``spline`` (default) interpolates between many local anchors, so it can
    follow the rolling baseline left by the tail of an intense solvent line.
    ``poly`` fits a single polynomial of the given ``order`` and is the safer
    choice when the baseline is only gently tilted.
    """
    if method == "poly":
        return _baseline_poly(data, order, segments)

    n = len(data)
    nodes = baseline_nodes(data, segments)
    if len(nodes) < 3:
        return list(data)

    # Median-smooth the anchor levels so a node that landed on a peak flank
    # gets pulled back down by its neighbours.
    levels = [v for _, v in nodes]
    smoothed = []
    for i in range(len(levels)):
        lo = max(0, i - 2)
        hi = min(len(levels), i + 3)
        smoothed.append(_median(levels[lo:hi]))

    xs = [i for i, _ in nodes]
    out = list(data)
    for seg in range(len(xs) - 1):
        x0, x1 = xs[seg], xs[seg + 1]
        y0, y1 = smoothed[seg], smoothed[seg + 1]
        span = x1 - x0
        if span <= 0:
            continue
        slope = (y1 - y0) / span
        for i in range(x0, min(x1, n)):
            out[i] = data[i] - (y0 + slope * (i - x0))
    for i in range(0, min(xs[0], n)):
        out[i] = data[i] - smoothed[0]
    for i in range(xs[-1], n):
        out[i] = data[i] - smoothed[-1]
    return out


def _baseline_poly(data: list, order: int, segments: int) -> list:
    """Subtract a polynomial fitted through the baseline anchors."""
    n = len(data)
    nodes = baseline_nodes(data, segments)
    if len(nodes) <= order:
        return list(data)

    # Least squares on a normalised abscissa keeps the normal equations sane.
    deg = order + 1
    xs = [(i / n) * 2.0 - 1.0 for i, _ in nodes]
    ys = [v for _, v in nodes]
    ata = [[0.0] * deg for _ in range(deg)]
    atb = [0.0] * deg
    for x, y in zip(xs, ys):
        powers = [x ** k for k in range(deg)]
        for r in range(deg):
            atb[r] += powers[r] * y
            for c in range(deg):
                ata[r][c] += powers[r] * powers[c]
    coef = solve_linear(ata, atb)

    out = []
    for i, v in enumerate(data):
        x = (i / n) * 2.0 - 1.0
        base = 0.0
        p = 1.0
        for c in coef:
            base += c * p
            p *= x
        out.append(v - base)
    return out


def noise_level(data: list, segments: int = 32) -> float:
    """Robust noise estimate: the standard deviation of the quietest segment."""
    n = len(data)
    if n < 2:
        return 0.0
    seg = max(8, n // segments)
    best = None
    for start in range(0, n - seg + 1, seg):
        chunk = data[start:start + seg]
        mean = sum(chunk) / len(chunk)
        var = sum((v - mean) ** 2 for v in chunk) / len(chunk)
        if best is None or var < best:
            best = var
    return math.sqrt(best) if best else 0.0
