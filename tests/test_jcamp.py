"""JCAMP-DX: ASDF decoding and the traps that real files contain."""

import os
import tempfile
import unittest

import fixtures
from nmranalyzer import nmrio


class TestASDF(unittest.TestCase):
    def decode(self, text):
        return nmrio._decode_asdf(text)[0]

    def test_pac(self):
        self.assertEqual(self.decode("+1 +2 +3"), [1, 2, 3])
        self.assertEqual(self.decode("1 -2 3"), [1, -2, 3])

    def test_sqz(self):
        self.assertEqual(self.decode("A2B2C2"), [12, 22, 32])
        self.assertEqual(self.decode("a2b2"), [-12, -22])
        self.assertEqual(self.decode("@0"), [0])

    def test_dif(self):
        # DIF letters are leading digits, so J0 is +10 rather than +0.
        self.assertEqual(self.decode("A0J0J0"), [10, 20, 30])
        self.assertEqual(self.decode("A0J5j5"), [10, 25, 10])
        self.assertEqual(self.decode("A0%%%"), [10, 10, 10, 10])

    def test_dup(self):
        self.assertEqual(self.decode("A5U"), [15, 15, 15])          # value
        self.assertEqual(self.decode("A0J2T"), [10, 22, 34])        # difference
        self.assertEqual(self.decode("A0J2V"), [10, 22, 34, 46, 58])

    def test_dup_count_is_multi_digit(self):
        # S1 is a count of 11, not 1 then 1.
        self.assertEqual(len(self.decode("A0J1S1")), 12)

    def test_decimal_abscissa(self):
        # SQZ files run the abscissa straight into the first ordinate.
        self.assertEqual(self.decode("5000.03B1399"), [5000.03, 21399])

    def test_E_is_a_sqz_digit_not_an_exponent(self):
        """The single nastiest trap in the format.

        ``E`` means the digit 5 in SQZ.  Reading ``1946E434`` as scientific
        notation yields 1946e434, i.e. infinity, and silently destroys the
        spectrum.
        """
        values = self.decode("1946E434")
        self.assertEqual(values, [1946.0, 5434.0])
        self.assertTrue(all(abs(v) < 1e6 for v in values))


class TestJcampFiles(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def path(self, name):
        return os.path.join(self.dir, name)

    def test_fixed_form(self):
        path = fixtures.jcamp_file(
            self.path("fix.jdx"),
            ["10.0 1 2 3 4", "6.0 5 6 7 8"], npoints=8)
        spec = nmrio.load(path)[0]
        self.assertEqual(spec.npoints, 8)
        # stored high-ppm first
        self.assertEqual(spec.real[0], 1.0)
        self.assertEqual(spec.real[-1], 8.0)

    def test_yfactor_is_applied(self):
        path = fixtures.jcamp_file(
            self.path("yf.jdx"), ["10.0 1 2 3 4", "6.0 5 6 7 8"],
            npoints=8, yfactor=2.5)
        spec = nmrio.load(path)[0]
        self.assertAlmostEqual(spec.real[0], 2.5, places=9)

    def test_dif_checkpoint_is_dropped(self):
        """A DIF row restates the previous row's last ordinate."""
        path = fixtures.jcamp_file(
            self.path("dif.jdx"),
            ["10.0 A00J3J4%%", "6.0  A27j8j9%"], npoints=8)
        spec = nmrio.load(path)[0]
        self.assertEqual(spec.npoints, 8)
        self.assertEqual([int(v) for v in spec.real],
                         [100, 113, 127, 127, 127, 109, 90, 90])

    def test_axis_spans_firstx_to_lastx(self):
        path = fixtures.jcamp_file(
            self.path("ax.jdx"), ["10.0 1 2 3 4", "6.0 5 6 7 8"],
            firstx=10.0, lastx=3.0, npoints=8)
        spec = nmrio.load(path)[0]
        self.assertAlmostEqual(spec.ppm(0), 10.0, places=6)
        self.assertAlmostEqual(spec.ppm(spec.npoints - 1), 3.0, places=6)

    def test_hz_abscissa_converts_to_ppm(self):
        lines = list(fixtures.JCAMP_HEADER)
        lines[3] = "##XUNITS=HZ"
        body = "\n".join(lines + [
            "##FIRSTX=4000.0", "##LASTX=400.0", "##YFACTOR=1.0",
            "##NPOINTS=8", "##XYDATA=(X++(Y..Y))",
            "4000.0 1 2 3 4", "2000.0 5 6 7 8", "##END="])
        path = self.path("hz.jdx")
        with open(path, "w") as fh:
            fh.write(body)
        spec = nmrio.load(path)[0]
        self.assertAlmostEqual(spec.ppm(0), 10.0, places=6)   # 4000 Hz / 400 MHz

    def test_multi_block_link_picks_the_spectrum(self):
        """Real exports wrap the data in a LINK document with other blocks."""
        body = "\n".join([
            "##TITLE=whole thing", "##JCAMP-DX=5.01", "##DATA TYPE=LINK",
            "##BLOCKS=2",
            "##TITLE=peak table", "##DATA TYPE=NMR PEAK ASSIGNMENTS",
            "##NPOINTS=2", "##XUNITS=PPM", "##FIRSTX=1.0", "##LASTX=2.0",
            "##YFACTOR=1.0", "##XYDATA=(X++(Y..Y))", "1.0 9 9", "##END=",
            "##TITLE=the data", "##DATA TYPE=NMR SPECTRUM",
            "##.OBSERVE FREQUENCY=400.0", "##XUNITS=PPM",
            "##FIRSTX=10.0", "##LASTX=3.0", "##YFACTOR=1.0", "##NPOINTS=8",
            "##XYDATA=(X++(Y..Y))", "10.0 1 2 3 4", "6.0 5 6 7 8", "##END=",
            "##END=",
        ])
        path = self.path("link.jdx")
        with open(path, "w") as fh:
            fh.write(body)
        spec = nmrio.load(path)[0]
        self.assertEqual(spec.npoints, 8)
        self.assertEqual(spec.sf, 400.0)


if __name__ == "__main__":
    unittest.main()
