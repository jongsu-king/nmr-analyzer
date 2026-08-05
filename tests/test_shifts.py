"""Chemical-shift estimates and carbon environments.

The tolerances here are deliberately loose: these are additivity rules, not a
calculation, and the point of the tests is that the estimates stay in the
right part of the spectrum rather than that they are precise.
"""

import unittest

import fixtures
from nmranalyzer import shifts, smiles


class TestProtonShifts(unittest.TestCase):
    # smiles, environment label fragment, literature shift in CDCl3
    CASES = [
        ("CCO", "C-H2", 3.69), ("CCO", "C-H3", 1.19),
        ("Cc1ccccc1", "aromatic", 7.20),
        ("O=Cc1ccccc1", "C-H", 10.02),
        ("c1ccccc1", "aromatic", 7.26),
        ("COc1ccccc1", "C-H3", 3.80),
        ("CC(=O)O", "O-H", 11.4),
    ]

    def _estimate(self, text, fragment):
        mol = smiles.parse(text)
        for env, estimate in shifts.predict_proton_environments(mol):
            if fragment in env.label:
                return estimate
        self.fail("no environment matching %r in %s" % (fragment, text))

    def test_estimates_are_in_the_right_region(self):
        for text, fragment, literature in self.CASES:
            with self.subTest(smiles=text, environment=fragment):
                estimate = self._estimate(text, fragment)
                self.assertLess(abs(estimate.value - literature), 1.0,
                                "%s %s: predicted %.2f, literature %.2f"
                                % (text, fragment, estimate.value, literature))

    def test_window_usually_brackets_the_literature_value(self):
        inside = 0
        for text, fragment, literature in self.CASES:
            if self._estimate(text, fragment).contains(literature):
                inside += 1
        self.assertGreaterEqual(inside, len(self.CASES) - 2,
                                "the quoted windows are too optimistic")

    def test_aromatic_substituent_effects_have_the_right_sign(self):
        """Nitro deshields the ring, methoxy shields it."""
        def ortho(text):
            mol = smiles.parse(text)
            values = [e.value for env, e in shifts.predict_proton_environments(mol)
                      if "aromatic" in env.label]
            return max(values)

        self.assertGreater(ortho("O=[N+]([O-])c1ccccc1"), 8.0)
        self.assertLess(min(
            e.value for env, e in
            shifts.predict_proton_environments(smiles.parse("COc1ccccc1"))
            if "aromatic" in env.label), 7.26)

    def test_exchangeable_protons_get_wide_windows(self):
        for text, fragment in (("CCO", "O-H"), ("CC(=O)O", "O-H")):
            estimate = self._estimate(text, fragment)
            self.assertGreaterEqual(estimate.window, 1.0,
                                    "OH shifts are concentration dependent and "
                                    "should not be quoted tightly")


class TestCarbonShifts(unittest.TestCase):
    def test_carbonyls_land_downfield(self):
        cases = [("CC(=O)C", 190.0), ("CC(=O)O", 165.0),
                 ("CCOC(C)=O", 160.0), ("O=Cc1ccccc1", 180.0)]
        for text, floor in cases:
            with self.subTest(smiles=text):
                mol = smiles.parse(text)
                top = max(e.value for _env, e
                          in shifts.predict_carbon_environments(mol))
                self.assertGreater(top, floor)

    def test_aromatic_carbons_near_128(self):
        mol = smiles.parse("c1ccccc1")
        env, estimate = shifts.predict_carbon_environments(mol)[0]
        self.assertAlmostEqual(estimate.value, 128.5, delta=3.0)

    def test_oxygen_deshields_the_attached_carbon(self):
        mol = smiles.parse("CCO")
        by_label = {env.label: est for env, est
                    in shifts.predict_carbon_environments(mol)}
        self.assertGreater(by_label["CH2"].value, by_label["CH3"].value + 20)


class TestCarbonEnvironments(unittest.TestCase):
    # smiles, number of distinct 13C signals
    CASES = [("CCO", 2), ("c1ccccc1", 1), ("Cc1ccccc1", 5),
             ("CCOC(C)=O", 4), ("CC(C)C", 2),
             ("c1ccc2[nH]c3ccccc3c2c1", 6)]     # carbazole is C2v: 12 C, 6 lines

    def test_signal_counts(self):
        for text, expected in self.CASES:
            with self.subTest(smiles=text):
                mol = smiles.parse(text)
                self.assertEqual(len(mol.carbon_environments()), expected)

    def test_every_carbon_belongs_to_exactly_one_environment(self):
        for text, _expected in self.CASES:
            mol = smiles.parse(text)
            covered = sum(env.count for env in mol.carbon_environments())
            self.assertEqual(covered, mol.formula_counts().get("C", 0))

    def test_quaternary_carbons_are_included(self):
        """Unlike 1H, a carbon with no hydrogens still gives a line."""
        mol = smiles.parse("CC(C)(C)c1ccccc1")
        labels = [env.label for env in mol.carbon_environments()]
        self.assertIn("quaternary C", labels)
        self.assertIn("aromatic C", labels)


if __name__ == "__main__":
    unittest.main()
