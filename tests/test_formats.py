"""Bruker 1D, ACD .esp and the 2D submatrix layout."""

import os
import struct
import tempfile
import unittest

import fixtures
from nmranalyzer import nmrio
from nmranalyzer import nmr2d


class TestBruker1D(unittest.TestCase):
    def test_reads_processed_data_and_axis(self):
        root = fixtures.write_bruker_1d(
            os.path.join(tempfile.mkdtemp(), "1"))
        spec = nmrio.load(root)[0]
        self.assertEqual(spec.npoints, 8192)
        self.assertAlmostEqual(spec.ppm(0), 12.0, places=6)
        self.assertEqual(spec.meta["Solvent"], "CDCl3")
        self.assertEqual(spec.meta["Nucleus"], "1H")

    def test_peaks_land_where_they_were_put(self):
        root = fixtures.write_bruker_1d(
            os.path.join(tempfile.mkdtemp(), "1"),
            peaks=((7.26, 1.0), (2.10, 0.4)))
        spec = nmrio.load(root)[0]
        top = max(range(spec.npoints), key=lambda i: spec.real[i])
        self.assertAlmostEqual(spec.ppm(top), 7.26, places=2)


class TestEsp(unittest.TestCase):
    def test_reads_header_and_flips_the_axis(self):
        path = fixtures.write_esp(
            os.path.join(tempfile.mkdtemp(), "t.esp"))
        spec = nmrio.load(path)[0]
        self.assertEqual(spec.meta["Format"], "ACD .esp")
        self.assertEqual(spec.meta["Solvent"], "TRIFLUOROACETIC ACID-d")
        # index 0 must be the high-ppm end after the flip
        self.assertGreater(spec.ppm(0), spec.ppm(spec.npoints - 1))

    def test_peak_position_survives_the_flip(self):
        path = fixtures.write_esp(
            os.path.join(tempfile.mkdtemp(), "t.esp"),
            peaks=((7.26, 1.0),))
        spec = nmrio.load(path)[0]
        top = max(range(spec.npoints), key=lambda i: spec.real[i])
        self.assertAlmostEqual(spec.ppm(top), 7.26, places=2)

    def test_point_count_is_twice_the_header_tag(self):
        path = fixtures.write_esp(
            os.path.join(tempfile.mkdtemp(), "t.esp"), npts=2048)
        spec = nmrio.load(path)[0]
        self.assertEqual(spec.npoints, 2048)


class TestBruker2D(unittest.TestCase):
    """The submatrix layout is the part most likely to be wrong."""

    def setUp(self):
        self.root, self.expected = fixtures.write_bruker_2d(
            os.path.join(tempfile.mkdtemp(), "1"))
        self.spec = nmr2d.read_bruker_2d(self.root)[0]

    def test_detiling_is_exact(self):
        worst = max(abs(a - b)
                    for row_a, row_b in zip(self.spec.data, self.expected)
                    for a, b in zip(row_a, row_b))
        self.assertEqual(worst, 0.0, "de-tiled matrix differs from the source")

    def test_row_major_read_would_have_been_wrong(self):
        """Proves the previous test has teeth.

        If the file happened to be row-major, de-tiling would be a no-op and
        the exactness check above would pass for the wrong reason.
        """
        rows, cols = self.spec.rows, self.spec.cols
        with open(os.path.join(self.root, "pdata", "1", "2rr"), "rb") as fh:
            raw = struct.unpack("<%di" % (rows * cols), fh.read())
        naive = [list(raw[r * cols:(r + 1) * cols]) for r in range(rows)]
        worst = max(abs(a - b)
                    for row_a, row_b in zip(naive, self.expected)
                    for a, b in zip(row_a, row_b))
        self.assertGreater(worst, 0.0)

    def test_axes(self):
        self.assertAlmostEqual(self.spec.f2.ppm(0), 10.0, places=6)
        self.assertAlmostEqual(self.spec.f1.ppm(0), 10.0, places=6)
        self.assertTrue(self.spec.is_homonuclear())

    def test_cross_peaks_are_found_at_the_right_shifts(self):
        peaks = nmr2d.pick_peaks_2d(self.spec, sensitivity=6.0)
        self.assertEqual(len(peaks), len(fixtures.COSY_PEAKS))
        for f2, f1, _height in fixtures.COSY_PEAKS:
            best = min(peaks, key=lambda p: abs(p.f2_ppm - f2) + abs(p.f1_ppm - f1))
            self.assertAlmostEqual(best.f2_ppm, f2, delta=0.05)
            self.assertAlmostEqual(best.f1_ppm, f1, delta=0.05)

    def test_diagonal_can_be_excluded(self):
        peaks = nmr2d.pick_peaks_2d(self.spec, sensitivity=6.0,
                                    skip_diagonal_ppm=0.3)
        self.assertEqual(len(peaks), 4)     # 7 total minus 3 on the diagonal


class TestContour(unittest.TestCase):
    def test_a_single_hill_gives_closed_contours(self):
        from nmranalyzer import contour
        size = 40
        grid = [[100.0 * (1.0 - ((r - 20) ** 2 + (c - 20) ** 2) / 400.0)
                 for c in range(size)] for r in range(size)]
        segments = contour.segments(grid, [20.0, 50.0])
        self.assertGreater(len(segments[20.0]), 8)
        self.assertGreater(len(segments[50.0]), 8)
        # the higher level must enclose a smaller area, so fewer segments
        self.assertLess(len(segments[50.0]), len(segments[20.0]))

    def test_flat_grid_has_no_contours(self):
        from nmranalyzer import contour
        grid = [[5.0] * 10 for _ in range(10)]
        segments = contour.segments(grid, [1.0, 9.0])
        self.assertEqual(sum(len(v) for v in segments.values()), 0)


if __name__ == "__main__":
    unittest.main()


class TestProcessingNumbers(unittest.TestCase):
    """TopSpin writes reprocessed data to pdata/2, pdata/3, ..."""

    def test_falls_back_when_pdata_1_is_absent(self):
        base = os.path.join(tempfile.mkdtemp(), "1")
        fixtures.write_bruker_1d(base)
        os.rename(os.path.join(base, "pdata", "1"),
                  os.path.join(base, "pdata", "2"))
        spec = nmrio.load(base)[0]
        self.assertEqual(spec.npoints, 8192)
        self.assertEqual(spec.meta["Processing no."], 2)

    def test_prefers_the_lowest_available(self):
        import shutil
        base = os.path.join(tempfile.mkdtemp(), "1")
        fixtures.write_bruker_1d(base)
        shutil.copytree(os.path.join(base, "pdata", "1"),
                        os.path.join(base, "pdata", "3"))
        spec = nmrio.load(base)[0]
        self.assertEqual(spec.meta["Processing no."], 1)

    def test_2d_falls_back_too(self):
        base = os.path.join(tempfile.mkdtemp(), "1")
        fixtures.write_bruker_2d(base)
        os.rename(os.path.join(base, "pdata", "1"),
                  os.path.join(base, "pdata", "4"))
        specs = nmr2d.read_bruker_2d(base)
        self.assertEqual(len(specs), 1)
        self.assertEqual(specs[0].meta["Processing no."], 4)
