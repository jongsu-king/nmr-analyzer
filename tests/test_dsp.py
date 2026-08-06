"""FFT, transform conventions, phasing and baseline."""

import cmath
import math
import unittest

import fixtures  # noqa: F401  (puts the project root on sys.path)
from nmranalyzer import dsp


def naive_dft(data):
    n = len(data)
    return [sum(data[k] * cmath.exp(-2j * math.pi * k * j / n)
                for k in range(n)) for j in range(n)]


class TestFFT(unittest.TestCase):
    def test_matches_naive_dft(self):
        data = [complex(math.sin(i * 0.7), math.cos(i * 0.3)) for i in range(64)]
        fast = dsp.fft(data)
        slow = naive_dft(data)
        worst = max(abs(a - b) for a, b in zip(fast, slow))
        self.assertLess(worst, 1e-9, "FFT disagrees with the direct DFT")

    def test_round_trip(self):
        data = [complex(i % 7, (i * 3) % 5) for i in range(128)]
        back = dsp.fft(dsp.fft(data), inverse=True)
        worst = max(abs(a - b) for a, b in zip(data, back))
        self.assertLess(worst, 1e-9)

    def test_pads_to_power_of_two(self):
        self.assertEqual(len(dsp.fft([1] * 100)), 128)


class TestTransform(unittest.TestCase):
    """Index 0 must be the high-frequency edge, Bruker style."""

    def test_axis_convention(self):
        n = 512
        sw = 1000.0
        dwell = 1.0 / sw
        # A tone at +sw/4 must land a quarter of the way in from the left.
        fid = [cmath.exp(2j * math.pi * (sw / 4.0) * k * dwell) for k in range(n)]
        spec = dsp.transform(fid, dwell, n, fcor=1.0)
        peak = max(range(n), key=lambda i: abs(spec[i]))
        self.assertEqual(peak, n // 4)

    def test_group_delay_is_a_pure_phase(self):
        spec = [complex(1.0, 0.0)] * 64
        shifted = dsp.group_delay_phase(spec, 7.5)
        for a, b in zip(spec, shifted):
            self.assertAlmostEqual(abs(a), abs(b), places=9)


class TestPhase(unittest.TestCase):
    def test_zero_order_rotates(self):
        spec = [complex(1.0, 0.0)] * 8
        out = dsp.phase(spec, 90.0, 0.0)
        for value in out:
            self.assertAlmostEqual(value.real, 0.0, places=9)
            self.assertAlmostEqual(value.imag, -1.0, places=9)

    def test_autophase_recovers_a_known_rotation(self):
        """A proper complex Lorentzian: real absorption, imaginary dispersion.

        The imaginary part matters.  With it set to zero every rotation under
        90 degrees leaves the real part non-negative, so the negative-area
        criterion has nothing to work with and the problem is ill-posed.
        """
        n = 2048
        spec = []
        for i in range(n):
            u = (i - n / 2) / 4.0
            spec.append(1.0 / complex(1.0, u))
        rotated = dsp.phase(spec, -40.0, 0.0)
        p0, _p1 = dsp.autophase(rotated, fit_p1=False)
        corrected = dsp.phase(rotated, p0, 0.0)
        self.assertGreater(max(v.real for v in corrected), 0.95)
        self.assertGreater(min(v.real for v in corrected), -0.06,
                           "dispersion left in the real part")


class TestBaseline(unittest.TestCase):
    def test_spline_follows_a_curved_baseline(self):
        n = 4096
        data = []
        for i in range(n):
            drift = 50.0 * math.sin(math.pi * i / n)
            peak = 100.0 if abs(i - 2000) < 8 else 0.0
            data.append(drift + peak)
        fixed = dsp.baseline_correct(data, method="spline", segments=128)
        flat = [v for i, v in enumerate(fixed) if abs(i - 2000) >= 40]
        self.assertLess(max(abs(v) for v in flat), 6.0,
                        "spline baseline left a curved residue")
        self.assertGreater(max(fixed), 80.0, "spline baseline ate the peak")

    def test_noise_level_is_positive_and_small(self):
        import random
        random.seed(3)
        data = [random.gauss(0.0, 2.0) for _ in range(4096)]
        estimate = dsp.noise_level(data)
        self.assertGreater(estimate, 0.5)
        self.assertLess(estimate, 4.0)


class TestSolver(unittest.TestCase):
    def test_solves_a_small_system(self):
        matrix = [[2.0, 1.0], [1.0, 3.0]]
        solution = dsp.solve_linear(matrix, [5.0, 10.0])
        self.assertAlmostEqual(solution[0], 1.0, places=9)
        self.assertAlmostEqual(solution[1], 3.0, places=9)


if __name__ == "__main__":
    unittest.main()


class TestPhaseRegularisation(unittest.TestCase):
    """A first-order term should only appear when the data calls for one."""

    def _lines(self, rotation_deg=0.0, n=4096):
        import cmath
        out = []
        for i in range(n):
            value = 0j
            for centre, height, width in ((1200, 1.0, 8), (2000, 0.6, 8),
                                          (2600, 0.9, 8)):
                u = (i - centre) / width
                value += height * (1 / (1 + u * u) + 1j * u / (1 + u * u))
            out.append(value * cmath.exp(-1j * math.radians(rotation_deg)))
        return out

    def test_pure_zero_order_rotation_gives_no_first_order_term(self):
        p0, p1 = dsp.autophase(self._lines(40.0))
        self.assertAlmostEqual(abs(p0), 40.0, delta=3.0)
        self.assertAlmostEqual(p1, 0.0, delta=2.0,
                               msg="a pure PH0 error should not need PH1")

    def test_an_already_phased_spectrum_is_left_nearly_alone(self):
        p0, p1 = dsp.autophase(self._lines(0.0))
        self.assertAlmostEqual(p0, 0.0, delta=5.0)
        self.assertAlmostEqual(p1, 0.0, delta=5.0)

    def test_the_penalty_does_not_block_a_needed_correction(self):
        """A real first-order error must still be found."""
        import cmath
        n = 4096
        spectrum = []
        for i in range(n):
            value = 0j
            for centre, height, width in ((600, 1.0, 8), (3400, 1.0, 8)):
                u = (i - centre) / width
                value += height * (1 / (1 + u * u) + 1j * u / (1 + u * u))
            # phase ramping across the spectrum: a genuine first-order error
            spectrum.append(value * cmath.exp(-1j * math.radians(120.0 * i / n)))
        _p0, p1 = dsp.autophase(spectrum)
        self.assertGreater(abs(p1), 40.0,
                           "a real first-order error was suppressed")
