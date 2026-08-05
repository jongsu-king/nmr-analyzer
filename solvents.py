"""Residual-solvent shifts and one-click calibration.

Reference values are the 1H residual-protio shifts tabulated by Fulmer et al.,
*Organometallics* 2010, 29, 2176, which is the table most groups reference
against.  TFA-d is listed at the shift ACD/Labs uses when it auto-references
these spectra.
"""

from __future__ import annotations


class Solvent:
    def __init__(self, key, label, shifts, aliases=()):
        self.key = key
        self.label = label
        self.shifts = shifts          # residual 1H shifts, ppm, strongest first
        self.aliases = aliases

    @property
    def primary(self):
        return self.shifts[0]

    def __repr__(self):
        return "Solvent(%s, %.2f ppm)" % (self.label, self.primary)


SOLVENTS = [
    Solvent("cdcl3", "CDCl3", [7.26], ("chloroform", "chloroform-d")),
    Solvent("dmso", "DMSO-d6", [2.50], ("dimethyl sulfoxide", "dmso-d6", "d6-dmso")),
    Solvent("d2o", "D2O", [4.79], ("deuterium oxide", "water-d2")),
    Solvent("acetone", "Acetone-d6", [2.05], ("acetone-d6", "d6-acetone")),
    Solvent("ch3cn", "CD3CN", [1.94], ("acetonitrile", "acetonitrile-d3", "cd3cn")),
    Solvent("meod", "CD3OD", [3.31, 4.87], ("methanol", "methanol-d4", "cd3od")),
    Solvent("c6d6", "C6D6", [7.16], ("benzene", "benzene-d6")),
    Solvent("thf", "THF-d8", [3.58, 1.72], ("tetrahydrofuran", "thf-d8")),
    Solvent("toluene", "Toluene-d8", [2.09, 7.09], ("toluene-d8",)),
    Solvent("dmf", "DMF-d7", [8.03, 2.92, 2.75], ("dimethylformamide", "dmf-d7")),
    Solvent("ch2cl2", "CD2Cl2", [5.32], ("dichloromethane", "methylene chloride")),
    Solvent("pyridine", "Pyridine-d5", [8.74, 7.58, 7.22], ("pyridine-d5",)),
    Solvent("tfa", "TFA-d", [11.50], ("trifluoroacetic acid", "trifluoroacetic acid-d",
                                      "tfa", "cf3cooh")),
]

BY_KEY = {s.key: s for s in SOLVENTS}


def _normalise(text):
    return "".join(ch for ch in text.lower() if ch.isalnum())


def identify(name):
    """Match a solvent string from file metadata to a known solvent.

    Handles both the short Bruker form (``TFA-d``) and the long ACD form
    (``TRIFLUOROACETIC ACID-d``).
    """
    if not name:
        return None
    target = _normalise(name)
    if not target:
        return None

    best = None
    best_len = 0
    for solvent in SOLVENTS:
        for candidate in (solvent.label,) + tuple(solvent.aliases) + (solvent.key,):
            token = _normalise(candidate)
            if not token:
                continue
            if token == target:
                return solvent
            # Prefer the longest alias that is contained in the metadata string,
            # so "trifluoroaceticacidd" does not match a short unrelated token.
            if (token in target or target in token) and len(token) > best_len:
                best, best_len = solvent, len(token)
    return best


def expected_shift(spec):
    """The residual-solvent shift expected for this spectrum, or ``None``."""
    solvent = identify(spec.meta.get("Solvent", ""))
    return solvent.primary if solvent else None


def calibrate(spec, target_ppm=None, search_ppm=0.6):
    """Return the ppm offset that puts the solvent line at its book value.

    Looks for the tallest point within ``search_ppm`` of where the residual
    solvent signal is expected and returns ``(delta, found_ppm, solvent)``.
    Returns ``None`` when the solvent is unknown or nothing is found.
    """
    solvent = identify(spec.meta.get("Solvent", ""))
    if target_ppm is None:
        if solvent is None:
            return None
        target_ppm = solvent.primary

    lo = spec.clamp(spec.index(target_ppm + search_ppm))
    hi = spec.clamp(spec.index(target_ppm - search_ppm))
    if lo > hi:
        lo, hi = hi, lo
    if hi <= lo:
        return None

    best = max(range(lo, hi + 1), key=lambda i: spec.real[i])
    found = spec.ppm(best)

    # Refine to sub-point accuracy with a parabola through the apex.
    if 0 < best < spec.npoints - 1:
        y0, y1, y2 = spec.real[best - 1], spec.real[best], spec.real[best + 1]
        denom = y0 - 2.0 * y1 + y2
        if denom:
            shift = max(-0.5, min(0.5, 0.5 * (y0 - y2) / denom))
            found = spec.ppm(0) + (best + shift) * spec.delta_ppm
    return target_ppm - found, found, solvent
