"""Peak picking, integration, multiplets, composition and lineshape fitting."""

import math
import unittest

import fixtures
from nmranalyzer import analysis
from nmranalyzer import fitting
from nmranalyzer import solvents
from nmranalyzer.nmrio import Spectrum

SF, SW, N = 500.0, 5000.0, 16384      # 10 ppm wide at 500 MHz
OFFSET = 10.0


def lorentzian_spectrum(peaks, width_hz=2.0, noise=0.0, seed=1):
    """``peaks`` is a list of (ppm, height)."""
    import random
    random.seed(seed)
    step = (SW / SF) / N
    data = []
    for i in range(N):
        ppm = OFFSET - i * step
        value = 0.0
        for centre, height in peaks:
            u = 2.0 * (ppm - centre) * SF / width_hz
            value += height / (1.0 + u * u)
        if noise:
            value += random.gauss(0.0, noise)
        data.append(value)
    return Spectrum("synthetic", data, SF, SW, OFFSET)


class TestPeakPicking(unittest.TestCase):
    def test_finds_isolated_peaks(self):
        spec = lorentzian_spectrum([(7.5, 1.0), (3.0, 0.6), (1.2, 0.3)],
                                   noise=0.002)
        peaks = analysis.pick_peaks(spec, sensitivity=8.0)
        found = sorted(round(p.ppm, 2) for p in peaks)
        self.assertEqual(found, [1.2, 3.0, 7.5])

    def test_prominence_rejects_ripples_on_a_broad_line(self):
        """Without a prominence test, noise on a big peak is picked as peaks."""
        spec = lorentzian_spectrum([(5.0, 1.0)], width_hz=40.0, noise=0.004)
        peaks = analysis.pick_peaks(spec, 4.0, 6.0, sensitivity=8.0)
        self.assertLessEqual(len(peaks), 3,
                             "ripples on the broad line were picked as peaks")

    def test_apex_is_interpolated(self):
        spec = lorentzian_spectrum([(7.333, 1.0)])
        peaks = analysis.pick_peaks(spec, 7.0, 7.6, sensitivity=8.0)
        self.assertEqual(len(peaks), 1)
        self.assertAlmostEqual(peaks[0].ppm, 7.333, places=3)


class TestIntegration(unittest.TestCase):
    def test_area_ratio_matches_height_ratio_for_equal_widths(self):
        spec = lorentzian_spectrum([(7.0, 3.0), (2.0, 1.0)])
        a = analysis.integrate(spec, 6.5, 7.5)
        b = analysis.integrate(spec, 1.5, 2.5)
        self.assertAlmostEqual(a / b, 3.0, delta=0.05)

    def test_refresh_peaks_tracks_the_data(self):
        spec = lorentzian_spectrum([(7.0, 1.0)])
        spec.peaks = analysis.pick_peaks(spec, sensitivity=8.0)
        before = spec.peaks[0].height
        spec.real = [v * 0.5 for v in spec.real]
        analysis.refresh_peaks(spec)
        self.assertAlmostEqual(spec.peaks[0].height, before * 0.5, places=6)


class TestMultiplets(unittest.TestCase):
    def make(self, centre, spacing_hz, heights):
        peaks = []
        n = len(heights)
        for k, height in enumerate(heights):
            offset = (k - (n - 1) / 2.0) * spacing_hz / SF
            peaks.append(analysis.Peak(0, centre + offset, height))
        return peaks

    def test_singlet(self):
        m = analysis.analyse_multiplet(self.make(7.0, 0, [1.0]), SF)
        self.assertEqual(m.pattern, "s")

    def test_doublet(self):
        m = analysis.analyse_multiplet(self.make(7.0, 8.0, [1.0, 1.0]), SF)
        self.assertEqual(m.pattern, "d")
        self.assertAlmostEqual(m.couplings[0], 8.0, delta=0.2)

    def test_triplet(self):
        m = analysis.analyse_multiplet(self.make(7.0, 7.0, [1.0, 2.0, 1.0]), SF)
        self.assertEqual(m.pattern, "t")
        self.assertAlmostEqual(m.couplings[0], 7.0, delta=0.2)

    def test_quartet(self):
        m = analysis.analyse_multiplet(
            self.make(7.0, 7.0, [1.0, 3.0, 3.0, 1.0]), SF)
        self.assertEqual(m.pattern, "q")

    def test_doublet_of_doublets(self):
        # lines at -(J1+J2)/2, -(J1-J2)/2, +(J1-J2)/2, +(J1+J2)/2
        j1, j2 = 12.0, 4.0
        offsets = [-(j1 + j2) / 2, -(j1 - j2) / 2, (j1 - j2) / 2, (j1 + j2) / 2]
        peaks = [analysis.Peak(0, 7.0 + o / SF, 1.0) for o in offsets]
        m = analysis.analyse_multiplet(peaks, SF)
        self.assertEqual(m.pattern, "dd")
        self.assertAlmostEqual(max(m.couplings), j1, delta=0.3)
        self.assertAlmostEqual(min(m.couplings), j2, delta=0.3)

    def test_irregular_stays_a_multiplet(self):
        peaks = [analysis.Peak(0, 7.0, 1.0), analysis.Peak(0, 6.98, 0.3),
                 analysis.Peak(0, 6.90, 0.9)]
        m = analysis.analyse_multiplet(peaks, SF)
        self.assertEqual(m.pattern, "m",
                         "an irregular group was forced into a pattern")


class TestComposition(unittest.TestCase):
    def test_conversion_and_ratio(self):
        class R:
            def __init__(self, v):
                self.value = v
        parts = [analysis.Component("SM", R(4.0), 4),
                 analysis.Component("P1", R(3.0), 3),
                 analysis.Component("P2", R(2.0), 2)]
        rows, conversion = analysis.composition(parts)
        self.assertEqual(len(rows), 3)
        for _label, moles, _frac in rows:
            self.assertAlmostEqual(moles, 1.0, places=9)
        self.assertAlmostEqual(conversion, 200.0 / 3.0, places=6)

    def test_report_does_not_invent_proton_counts(self):
        spec = lorentzian_spectrum([(7.0, 1.0)])
        region = analysis.Region(6.5, 7.5)
        analysis.integrate_region(spec, region)
        text = analysis.format_report(spec, [region])
        self.assertIn("rel.", text)
        # the only "1H" allowed is the nucleus in the header, not a made-up
        # integration
        body = text.split("delta", 1)[1]
        self.assertNotIn("H", body)


class TestFitting(unittest.TestCase):
    def test_resolves_two_overlapping_lines(self):
        width = 3.0
        spec = lorentzian_spectrum([(1.000, 1.0), (0.992, 0.5)],
                                   width_hz=width)
        result = fitting.fit_region(spec, 0.96, 1.03)
        self.assertEqual(len(result.peaks), 2)
        expected = [math.pi * 1.0 * width / 2.0, math.pi * 0.5 * width / 2.0]
        for peak, want in zip(result.peaks, expected):
            self.assertAlmostEqual(peak.area / want, 1.0, delta=0.02)
        self.assertAlmostEqual(result.peaks[0].fwhm_hz, width, delta=0.1)

    def test_reports_its_own_quality(self):
        spec = lorentzian_spectrum([(5.0, 1.0)], width_hz=3.0)
        result = fitting.fit_region(spec, 4.9, 5.1)
        self.assertLess(result.rel_rms, 0.05)
        self.assertTrue(result.converged)


class TestSolvents(unittest.TestCase):
    def test_identifies_short_and_long_spellings(self):
        self.assertEqual(solvents.identify("TFA-d").label, "TFA-d")
        self.assertEqual(
            solvents.identify("TRIFLUOROACETIC ACID-d").label, "TFA-d")
        self.assertEqual(solvents.identify("CDCl3").label, "CDCl3")
        self.assertEqual(solvents.identify("Chloroform-D").label, "CDCl3")
        self.assertIsNone(solvents.identify("banana"))

    def test_calibration_offset(self):
        spec = lorentzian_spectrum([(7.10, 1.0)])
        spec.meta["Solvent"] = "CDCl3"
        delta, found, solvent = solvents.calibrate(spec)
        self.assertAlmostEqual(found, 7.10, places=2)
        self.assertAlmostEqual(delta, 7.26 - found, places=6)
        self.assertEqual(solvent.label, "CDCl3")


if __name__ == "__main__":
    unittest.main()


class TestAssignment(unittest.TestCase):
    """Regions can be linked to specific atoms of a structure."""

    def setUp(self):
        self.region = analysis.Region(1.10, 1.30)
        self.region.value = 3.0

    def test_a_region_starts_unassigned(self):
        self.assertIsNone(self.region.assignment)
        self.assertEqual(self.region.assignment_label, "")

    def test_assignment_reaches_the_report(self):
        spec = fixtures.synthetic_spectrum()
        region = analysis.Region(1.10, 1.30)
        analysis.integrate_region(spec, region)
        region.protons = 3
        region.assignment = [0]
        region.assignment_label = "C-H3"
        text = analysis.format_report(spec, [region])
        self.assertIn("C-H3", text)
        self.assertIn("3H", text)

    def test_report_omits_the_label_when_unassigned(self):
        spec = fixtures.synthetic_spectrum()
        region = analysis.Region(1.10, 1.30)
        analysis.integrate_region(spec, region)
        region.protons = 3
        text = analysis.format_report(spec, [region])
        self.assertNotIn("C-H3", text)
