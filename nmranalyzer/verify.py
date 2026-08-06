"""Decide whether a set of measured signals is consistent with a structure.

The measurement is reduced to what an integration actually yields — a list of
(chemical shift, number of protons) — and matched against the environments the
structure predicts, on **both** counts: a signal must be the right size *and*
sit where that environment is expected.

Size alone is not a test. The scale factor relating integrals to protons is
free, so with a handful of signals almost any structure can be made to land on
plausible integers; requiring the shift to agree as well is what lets a wrong
structure fail.

This module holds the logic so that the interface and the benchmark exercise
the same code.
"""

from __future__ import annotations

from . import shifts

# A signal may differ from the predicted proton count by this much before the
# match is rejected: an absolute floor, or a fraction of the count, whichever
# is larger.
SIZE_FLOOR = 0.3
SIZE_FRACTION = 0.15

# How far outside its predicted window a shift may fall, as a multiple of the
# window's own half-width.
SHIFT_TOLERANCE = 1.0

EXCHANGEABLE = ("O", "N", "S")


class Match:
    def __init__(self, signal, environment, estimate):
        self.signal = signal              # (shift, protons)
        self.environment = environment    # ProtonEnvironment, or None
        self.estimate = estimate          # shifts.Estimate, or None

    @property
    def matched(self):
        return self.environment is not None

    def __repr__(self):
        if not self.matched:
            return "Match(%.2f ppm, %.2f H -> none)" % self.signal
        return "Match(%.2f ppm, %.2f H -> %dH %s)" % (
            self.signal[0], self.signal[1],
            self.environment.count, self.environment.label)


class Result:
    def __init__(self, matches, unaccounted, expected_protons, skipped):
        self.matches = matches
        self.unaccounted = unaccounted      # predicted environments not used
        self.expected_protons = expected_protons
        self.skipped = skipped              # exchangeable protons excluded

    @property
    def matched(self):
        return sum(1 for m in self.matches if m.matched)

    @property
    def total(self):
        return len(self.matches)

    @property
    def consistent(self):
        """Every signal matched, and no predicted environment left over."""
        return (self.total > 0
                and self.matched == self.total
                and not self.unaccounted)

    def summary(self):
        if self.consistent:
            return ("all %d signals match a predicted environment"
                    % self.total)
        parts = []
        if self.matched < self.total:
            parts.append("%d of %d signals match no environment"
                         % (self.total - self.matched, self.total))
        if self.unaccounted:
            parts.append("%d environment(s) unaccounted for"
                         % len(self.unaccounted))
        return "; ".join(parts)


def candidates(molecule, ignore_exchangeable=False):
    """Predicted proton environments, with their estimated shifts."""
    out = list(shifts.predict_proton_environments(molecule))
    if ignore_exchangeable:
        out = [(env, est) for env, est in out
               if env.carrier.symbol not in EXCHANGEABLE]
    return out


def expected_protons(molecule, ignore_exchangeable=False):
    """Protons the formula predicts, and how many were set aside."""
    total = molecule.formula_counts().get("H", 0)
    if not ignore_exchangeable:
        return total, 0
    skipped = sum(env.count for env in molecule.proton_environments()
                  if env.carrier.symbol in EXCHANGEABLE)
    return total - skipped, skipped


def _fits(environment, estimate, protons, shift):
    size_ok = abs(environment.count - protons) <= max(
        SIZE_FLOOR, SIZE_FRACTION * environment.count)
    if estimate.contains(shift):
        shift_ok = True
    else:
        gap = min(abs(shift - estimate.low), abs(shift - estimate.high))
        shift_ok = gap <= SHIFT_TOLERANCE * estimate.window
    return size_ok, shift_ok


def check(molecule, signals, ignore_exchangeable=False):
    """Match measured ``(shift, protons)`` signals to predicted environments.

    Signals are taken in order of decreasing shift and greedily assigned to
    the best remaining environment, which is what a chemist does by hand.
    """
    wanted, skipped = expected_protons(molecule, ignore_exchangeable)
    remaining = candidates(molecule, ignore_exchangeable)

    matches = []
    for shift, protons in sorted(signals, key=lambda s: -s[0]):
        best = None
        best_cost = None
        for env, estimate in remaining:
            size_off = abs(env.count - protons) / max(env.count, 1.0)
            if estimate.contains(shift):
                shift_off = 0.0
            else:
                gap = min(abs(shift - estimate.low), abs(shift - estimate.high))
                shift_off = gap / max(estimate.window, 0.2)
            cost = size_off + shift_off
            if best_cost is None or cost < best_cost:
                best, best_cost = (env, estimate), cost

        if best is not None:
            env, estimate = best
            size_ok, shift_ok = _fits(env, estimate, protons, shift)
            if size_ok and shift_ok:
                remaining.remove(best)
                matches.append(Match((shift, protons), env, estimate))
                continue
        matches.append(Match((shift, protons), None, None))

    return Result(matches, [env for env, _e in remaining], wanted, skipped)


def scale_to_formula(values, expected):
    """Scale raw integrals so that they total the expected proton count."""
    total = sum(v for v in values if v > 0)
    if total <= 0 or expected <= 0:
        return None
    return expected / total
