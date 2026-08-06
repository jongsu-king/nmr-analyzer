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


class TestTopicity(unittest.TestCase):
    """Diastereotopic protons give separate signals; equivalent ones do not."""

    def counts(self, text, split=True):
        mol = smiles.parse(text)
        return sorted((e.count for e in mol.proton_environments(split)),
                      reverse=True)

    def test_plain_molecules_are_untouched(self):
        """Nothing splits without a stereocentre, a ring or a prochiral centre."""
        for text, expected in (("CCO", [3, 2, 1]),
                               ("CCOC(C)=O", [3, 3, 2]),
                               ("Cc1ccccc1", [3, 2, 2, 1]),
                               ("CC(O)C", [6, 1, 1])):
            with self.subTest(smiles=text):
                self.assertEqual(self.counts(text), expected)

    def test_ch2_next_to_a_stereocentre_splits(self):
        """Butan-2-ol: the C3 protons are diastereotopic."""
        self.assertEqual(self.counts("CC(O)CC", split=False), [3, 3, 2, 1, 1])
        self.assertEqual(self.counts("CC(O)CC"), [3, 3, 1, 1, 1, 1])

    def test_declared_stereochemistry_is_honoured(self):
        mol = smiles.parse("C[C@H](O)CC")
        self.assertTrue(mol.stereocentres())
        self.assertEqual(self.counts("C[C@H](O)CC"), [3, 3, 1, 1, 1, 1])

    def test_isopropyl_methyls_split_next_to_a_stereocentre(self):
        """Valine's two methyls give two doublets, not one six-proton signal."""
        self.assertEqual(self.counts("CC(C)C(N)C(=O)O", split=False),
                         [6, 2, 1, 1, 1])
        self.assertEqual(self.counts("CC(C)C(N)C(=O)O"), [3, 3, 2, 1, 1, 1])

    def test_ring_ch2_protons_split(self):
        """The two faces of a substituted ring are not the same."""
        self.assertEqual(self.counts("OC1CCCCC1", split=False), [4, 4, 2, 1, 1])
        self.assertEqual(self.counts("OC1CCCCC1"), [2, 2, 2, 2, 2, 1, 1])

    def test_prochiral_centre_without_a_stereocentre(self):
        """Glycerol has no stereocentre but C2 is prochiral."""
        mol = smiles.parse("OCC(O)CO")
        self.assertFalse(mol.stereocentres())
        self.assertTrue(mol.prochiral_centres())
        self.assertEqual(self.counts("OCC(O)CO"), [2, 2, 2, 1, 1])

    def test_splitting_conserves_the_proton_count(self):
        for text in ("CC(O)CC", "CC(C)C(N)C(=O)O", "OC1CCCCC1",
                     "CC(C)CC(N)C(=O)O", "OCC(O)CO"):
            with self.subTest(smiles=text):
                mol = smiles.parse(text)
                total = sum(e.count for e in mol.proton_environments())
                self.assertEqual(total, mol.formula_counts()["H"])

    def test_labels_name_the_right_group(self):
        """A split methyl must not be labelled as a CH2."""
        mol = smiles.parse("CC(C)C(N)C(=O)O")
        split = [e for e in mol.proton_environments() if e.diastereotopic]
        self.assertTrue(split)
        for env in split:
            self.assertIn("H3", env.label)


class TestExplicitHydrogens(unittest.TestCase):
    """Database SMILES write every hydrogen out; both spellings must agree."""

    PAIRS = [
        ("CCO", "C([H])([H])C([H])([H])O[H]"),
        ("Cc1ccccc1", "C([H])([H])([H])c1c([H])c([H])c([H])c([H])c1[H]"),
        ("CC(=O)O", "C([H])([H])([H])C(=O)O[H]"),
    ]

    def test_formula_is_the_same_either_way(self):
        for implicit, explicit in self.PAIRS:
            with self.subTest(smiles=implicit):
                self.assertEqual(smiles.parse(implicit).formula(),
                                 smiles.parse(explicit).formula())

    def test_environments_are_the_same_either_way(self):
        for implicit, explicit in self.PAIRS:
            with self.subTest(smiles=implicit):
                a = [(e.count, e.label)
                     for e in smiles.parse(implicit).proton_environments()]
                b = [(e.count, e.label)
                     for e in smiles.parse(explicit).proton_environments()]
                self.assertEqual(a, b)

    def test_explicit_hydrogens_leave_no_stray_atoms(self):
        mol = smiles.parse("C([H])([H])([H])O[H]")
        self.assertFalse([a for a in mol.atoms if a.symbol == "H"],
                         "explicit H should have been folded into its neighbour")
        self.assertEqual([a.index for a in mol.atoms],
                         list(range(len(mol.atoms))), "indices must stay dense")

    def test_deuterium_is_not_a_proton(self):
        """CD3OH shows one proton signal, not four."""
        mol = smiles.parse("C([2H])([2H])([2H])O[H]")
        environments = mol.proton_environments()
        self.assertEqual(len(environments), 1)
        self.assertEqual(environments[0].count, 1)
        self.assertIn("O-H", environments[0].label)
