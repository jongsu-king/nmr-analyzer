"""Session round-trip and SVG export."""

import os
import tempfile
import unittest
import xml.etree.ElementTree as ET

import fixtures
from nmranalyzer import analysis
from nmranalyzer import export
from nmranalyzer import nmrio


class TestSession(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        root = fixtures.write_bruker_1d(os.path.join(self.dir, "1"))
        self.spec = nmrio.load(root)[0]
        for lo, hi in ((7.0, 7.5), (1.9, 2.3)):
            region = analysis.Region(lo, hi)
            analysis.integrate_region(self.spec, region)
            self.spec.regions.append(region)
        self.spec.regions[0].protons = 2
        self.spec.ref_shift = -0.15
        self.spec.baseline_on = True
        self.spec.color = "#123456"

    def test_round_trip_preserves_analysis(self):
        path = os.path.join(self.dir, "s.nmrs")
        export.save_session(path, [self.spec], (12.0, -1.0),
                            {"stack": True, "sensitivity": 11.0})
        restored, view, options, warnings = export.load_session(path)

        self.assertEqual(warnings, [])
        self.assertEqual(len(restored), 1)
        back = restored[0]
        self.assertEqual(view, (12.0, -1.0))
        self.assertTrue(options["stack"])
        self.assertEqual(options["sensitivity"], 11.0)
        self.assertEqual(back.color, "#123456")
        self.assertAlmostEqual(back.ref_shift, -0.15, places=9)
        self.assertTrue(back.baseline_on)
        self.assertEqual(len(back.regions), 2)
        self.assertEqual(back.regions[0].protons, 2)
        for before, after in zip(self.spec.regions, back.regions):
            self.assertAlmostEqual(before.lo, after.lo, places=9)
            self.assertAlmostEqual(before.hi, after.hi, places=9)

    def test_missing_source_is_reported_not_fatal(self):
        path = os.path.join(self.dir, "s.nmrs")
        export.save_session(path, [self.spec], (12.0, -1.0), {})
        import json
        with open(path) as fh:
            payload = json.load(fh)
        payload["spectra"][0]["source"] = "/nowhere/at/all"
        with open(path, "w") as fh:
            json.dump(payload, fh)
        restored, _view, _options, warnings = export.load_session(path)
        self.assertEqual(restored, [])
        self.assertTrue(any("missing" in w for w in warnings))

    def test_rejects_a_foreign_file(self):
        path = os.path.join(self.dir, "other.json")
        with open(path, "w") as fh:
            fh.write('{"format": "something else"}')
        self.assertRaises(ValueError, export.load_session, path)


class TestSvgExport(unittest.TestCase):
    def test_writes_well_formed_svg_with_labels(self):
        directory = tempfile.mkdtemp()
        root = fixtures.write_bruker_1d(os.path.join(directory, "1"))
        spec = nmrio.load(root)[0]
        region = analysis.Region(7.0, 7.5)
        analysis.integrate_region(spec, region)
        region.protons = 1
        spec.regions.append(region)

        path = os.path.join(directory, "plot.svg")
        export.write_svg(path, [spec], 9.0, 6.0, active=spec)
        root_element = ET.parse(path).getroot()
        self.assertTrue(root_element.tag.endswith("svg"))
        texts = [e.text for e in root_element.iter()
                 if e.tag.endswith("text") and e.text]
        self.assertIn("ppm", texts)
        self.assertIn("1H", texts)          # the assigned integral label

    def test_refuses_an_empty_plot(self):
        self.assertRaises(ValueError, export.write_svg,
                          "/tmp/none.svg", [], 9.0, 6.0)


if __name__ == "__main__":
    unittest.main()
