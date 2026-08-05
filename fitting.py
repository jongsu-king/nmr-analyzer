"""Lineshape fitting (deconvolution) for overlapping NMR signals.

Simple summation gets an integral wrong when two multiplets overlap, because
the answer then depends on where you draw the boundary.  Fitting a sum of
pseudo-Voigt lines instead assigns the shared intensity to each component
according to its shape, and the areas come out analytically.

The optimiser is Levenberg-Marquardt with analytic derivatives; no third-party
packages are involved.
"""

from __future__ import annotations

import math

import dsp

LN2 = math.log(2.0)
SQRT_PI_OVER_LN2 = math.sqrt(math.pi / LN2)


class FittedPeak:
    """One resolved component line."""

    def __init__(self, ppm, height, fwhm_hz, area):
        self.ppm = ppm
        self.height = height
        self.fwhm_hz = fwhm_hz
        self.area = area

    def __repr__(self):
        return "FittedPeak(%.4f ppm, fwhm %.2f Hz, area %.4g)" % (
            self.ppm, self.fwhm_hz, self.area)


class FitResult:
    def __init__(self, peaks, eta, baseline, rms, rel_rms, iterations,
                 converged, center_hz, sf):
        self.peaks = peaks
        self.eta = eta                  # 1 = pure Lorentzian, 0 = pure Gaussian
        self.baseline = baseline        # (offset, slope per Hz)
        self.rms = rms
        self.rel_rms = rel_rms          # rms / peak height, a fit-quality gauge
        self.iterations = iterations
        self.converged = converged
        self.center_hz = center_hz
        self.sf = sf

    @property
    def total_area(self):
        return sum(p.area for p in self.peaks)

    def evaluate(self, ppm_values):
        """Model intensity at the given ppm positions, baseline included."""
        params = []
        for peak in self.peaks:
            params.extend((peak.ppm * self.sf - self.center_hz,
                           peak.height, peak.fwhm_hz))
        params.extend((self.eta, self.baseline[0], self.baseline[1]))
        xs = [p * self.sf - self.center_hz for p in ppm_values]
        return [_model(params, len(self.peaks), x) for x in xs]

    def component(self, index, ppm_values):
        """One component line on its own, without the baseline."""
        peak = self.peaks[index]
        x0 = peak.ppm * self.sf - self.center_hz
        out = []
        for p in ppm_values:
            x = p * self.sf - self.center_hz
            out.append(_shape(x - x0, peak.height, peak.fwhm_hz, self.eta))
        return out


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


def _shape(dx, height, fwhm, eta):
    if fwhm <= 0.0:
        return 0.0
    u = 2.0 * dx / fwhm
    lorentz = height / (1.0 + u * u)
    gauss = height * math.exp(-LN2 * u * u)
    return eta * lorentz + (1.0 - eta) * gauss


def _model(params, n_peaks, x):
    total = params[3 * n_peaks + 1] + params[3 * n_peaks + 2] * x
    eta = params[3 * n_peaks]
    for k in range(n_peaks):
        x0, h, w = params[3 * k], params[3 * k + 1], params[3 * k + 2]
        total += _shape(x - x0, h, w, eta)
    return total


def _model_and_jacobian(params, n_peaks, xs):
    """Model values and d(model)/d(parameter) at every sample point."""
    eta = params[3 * n_peaks]
    b0 = params[3 * n_peaks + 1]
    b1 = params[3 * n_peaks + 2]
    n_par = 3 * n_peaks + 3

    values = []
    jac = []
    for x in xs:
        total = b0 + b1 * x
        row = [0.0] * n_par
        d_eta = 0.0
        for k in range(n_peaks):
            x0 = params[3 * k]
            h = params[3 * k + 1]
            w = params[3 * k + 2]
            if w <= 0.0:
                continue
            u = 2.0 * (x - x0) / w
            denom = 1.0 + u * u
            lor_unit = 1.0 / denom
            gau_unit = math.exp(-LN2 * u * u)
            unit = eta * lor_unit + (1.0 - eta) * gau_unit
            total += h * unit

            # d/du of the unit-height shape
            d_unit_du = (eta * (-2.0 * u / (denom * denom))
                         + (1.0 - eta) * (-2.0 * LN2 * u * gau_unit))
            row[3 * k] = h * d_unit_du * (-2.0 / w)       # d/dx0
            row[3 * k + 1] = unit                          # d/dh
            row[3 * k + 2] = h * d_unit_du * (-u / w)      # d/dw
            d_eta += h * (lor_unit - gau_unit)
        row[3 * n_peaks] = d_eta
        row[3 * n_peaks + 1] = 1.0
        row[3 * n_peaks + 2] = x
        values.append(total)
        jac.append(row)
    return values, jac


def _area(height, fwhm, eta):
    """Analytic area of a pseudo-Voigt line."""
    lorentz = height * math.pi * fwhm / 2.0
    gauss = height * fwhm * 0.5 * SQRT_PI_OVER_LN2
    return eta * lorentz + (1.0 - eta) * gauss


# ---------------------------------------------------------------------------
# Initial guesses
# ---------------------------------------------------------------------------


def _estimate_fwhm(xs, ys, index, height, floor):
    """Half-height width around a maximum, in the units of ``xs``."""
    half = floor + (height - floor) / 2.0
    left = index
    while left > 0 and ys[left] > half:
        left -= 1
    right = index
    while right < len(ys) - 1 and ys[right] > half:
        right += 1
    width = abs(xs[right] - xs[left])
    return width if width > 0 else abs(xs[1] - xs[0]) * 2.0


# ---------------------------------------------------------------------------
# Levenberg-Marquardt
# ---------------------------------------------------------------------------


def fit_region(spec, lo_ppm, hi_ppm, seed_peaks=None, max_iterations=60,
               tolerance=1e-6, fit_shape=True):
    """Deconvolute the signals between ``lo_ppm`` and ``hi_ppm``.

    ``seed_peaks`` are :class:`analysis.Peak` objects used as starting
    positions; if omitted the local maxima of the region are used.  Returns a
    :class:`FitResult`, or ``None`` when there is nothing to fit.
    """
    a = spec.clamp(spec.index(hi_ppm))
    b = spec.clamp(spec.index(lo_ppm))
    if a > b:
        a, b = b, a
    if b - a < 6:
        return None

    ppm_values = [spec.ppm(i) for i in range(a, b + 1)]
    ys = list(spec.real[a:b + 1])
    center_hz = ((ppm_values[0] + ppm_values[-1]) / 2.0) * spec.sf
    xs = [p * spec.sf - center_hz for p in ppm_values]

    # Seed positions
    if seed_peaks:
        seeds = [p for p in seed_peaks if lo_ppm <= p.ppm <= hi_ppm]
    else:
        seeds = []
    if not seeds:
        import analysis
        seeds = analysis.pick_peaks(spec, lo_ppm, hi_ppm, sensitivity=4.0)
    if not seeds:
        top = max(range(len(ys)), key=lambda i: ys[i])
        seeds = [type("S", (), {"ppm": ppm_values[top], "height": ys[top]})()]

    floor = min(ys)
    params = []
    for seed in seeds:
        idx = min(range(len(ppm_values)),
                  key=lambda i: abs(ppm_values[i] - seed.ppm))
        height = ys[idx] - floor
        if height <= 0:
            height = max(ys) - floor or 1.0
        fwhm = _estimate_fwhm(xs, ys, idx, ys[idx], floor)
        params.extend((xs[idx], height, max(fwhm, 0.2)))
    n_peaks = len(seeds)
    params.extend((0.8 if fit_shape else 1.0, floor, 0.0))

    n_par = len(params)
    lam = 1e-3
    prev_cost = None
    iterations = 0
    converged = False

    for iterations in range(1, max_iterations + 1):
        values, jac = _model_and_jacobian(params, n_peaks, xs)
        residual = [y - v for y, v in zip(ys, values)]
        cost = sum(r * r for r in residual)

        # Normal equations J^T J dp = J^T r
        jtj = [[0.0] * n_par for _ in range(n_par)]
        jtr = [0.0] * n_par
        for row, r in zip(jac, residual):
            for i in range(n_par):
                ri = row[i]
                if ri == 0.0:
                    continue
                jtr[i] += ri * r
                target = jtj[i]
                for j in range(i, n_par):
                    target[j] += ri * row[j]
        for i in range(n_par):                 # mirror the symmetric half
            for j in range(i):
                jtj[i][j] = jtj[j][i]

        improved = False
        for _ in range(8):                     # damping search
            damped = [row[:] for row in jtj]
            for i in range(n_par):
                damped[i][i] += lam * (jtj[i][i] if jtj[i][i] > 0 else 1.0)
            step = dsp.solve_linear(damped, jtr)
            trial = [p + s for p, s in zip(params, step)]

            # Keep the parameters physical.
            for k in range(n_peaks):
                trial[3 * k + 2] = max(0.05, trial[3 * k + 2])
            trial[3 * n_peaks] = min(1.0, max(0.0, trial[3 * n_peaks]))
            if not fit_shape:
                trial[3 * n_peaks] = params[3 * n_peaks]

            trial_values = [_model(trial, n_peaks, x) for x in xs]
            trial_cost = sum((y - v) ** 2 for y, v in zip(ys, trial_values))
            if trial_cost < cost:
                params = trial
                lam = max(lam * 0.3, 1e-9)
                improved = True
                break
            lam = min(lam * 10.0, 1e9)

        if not improved:
            # The damping search ran all the way to a near-zero step without
            # lowering the cost, which means this is the optimum.  Reporting
            # it as "not converged" would make a perfect fit look suspect.
            converged = True
            break
        if prev_cost is not None and abs(prev_cost - trial_cost) <= \
                tolerance * max(prev_cost, 1e-30):
            converged = True
            prev_cost = trial_cost
            break
        prev_cost = trial_cost

    eta = params[3 * n_peaks]
    peaks = []
    for k in range(n_peaks):
        x0, h, w = params[3 * k], params[3 * k + 1], params[3 * k + 2]
        peaks.append(FittedPeak(ppm=(x0 + center_hz) / spec.sf,
                                height=h, fwhm_hz=w,
                                area=_area(h, w, eta)))
    peaks.sort(key=lambda p: -p.ppm)

    final = [_model(params, n_peaks, x) for x in xs]
    rms = math.sqrt(sum((y - v) ** 2 for y, v in zip(ys, final)) / len(ys))
    tallest = max((p.height for p in peaks), default=0.0)
    rel = rms / tallest if tallest else float("inf")

    return FitResult(peaks, eta, (params[3 * n_peaks + 1],
                                  params[3 * n_peaks + 2]),
                     rms, rel, iterations, converged, center_hz, spec.sf)


def area_to_ppm_scale(area_hz, sf):
    """Convert a fitted area (intensity x Hz) to the ppm-based integral scale.

    :func:`analysis.integrate` sums intensity x ppm, so fitted areas have to be
    divided by the observe frequency before they can be compared with it.
    """
    return area_hz / sf if sf else area_hz
