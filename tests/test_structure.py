"""SMILES parsing, molecular arithmetic, symmetry and 2D depiction."""

import math
import unittest

import fixtures
import depict
import smiles


class TestParsing(unittest.TestCase):
    CASES = [
        # smiles, formula, MW, DBE, proton environments (counts, descending)
        ("CCO", "C2H6O", 46.07, 0, [3, 2, 1]),
        ("c1ccccc1", "C6H6", 78.11, 4, [6]),
        ("Cc1ccccc1", "C7H8", 92.14, 4, [3, 2, 2, 1]),
        ("CCOC(C)=O", "C4H8O2", 88.11, 1, [3, 3, 2]),
        ("CC(=O)O", "C2H4O2", 60.05, 1, [3, 1]),
        ("ClC(Cl)Cl", "CHCl3", 119.37, 0, [1]),
        ("O", "H2O", 18.02, 0, [2]),
        ("CC(C)C", "C4H10", 58.12, 0, [9, 1]),
        ("OCCO", "C2H6O2", 62.07, 0, [4, 2]),
        ("Ic1ccccc1", "C6H5I", 204.01, 4, [2, 2, 1]),
        ("CC(C)(C)c1ccccc1", "C10H14", 134.22, 4, [9, 2, 2, 1]),
        ("O=Cc1ccccc1", "C7H6O", 106.12, 5, [2, 2, 1, 1]),
    ]

    def test_formula_mass_dbe_and_environments(self):
        for text, formula, mass, dbe, envs in self.CASES:
            with self.subTest(smiles=text):
                mol = smiles.parse(text)
                self.assertEqual(mol.formula(), formula)
                self.assertAlmostEqual(mol.molecular_weight(), mass, delta=0.05)
                self.assertEqual(mol.dbe(), dbe)
                self.assertEqual([e.count for e in mol.proton_environments()],
                                 envs)

    def test_rejects_malformed_input(self):
        for bad in ("", "C(", "C)", "c1ccccc", "C%9"):
            with self.subTest(smiles=bad):
                self.assertRaises(smiles.SmilesError, smiles.parse, bad)


class TestAromaticHydrogens(unittest.TestCase):
    """Aromatic nitrogen is the case that is easy to get wrong.

    Carbon always carries a formal double bond in the Kekule structure, so it
    uses one bond more than its sigma count.  Nitrogen only does so in the
    pyridine sense: once it has three sigma neighbours it is donating its lone
    pair instead, and adding the extra bond invents a hydrogen.
    """

    CASES = [
        ("c1ccncc1", "C5H5N", [2, 2, 1]),              # pyridine, 2-coordinate n
        ("Cn1cccc1", "C5H7N", [3, 2, 2]),              # N-methylpyrrole, 3-coord
        ("[nH]1cccc1", "C4H5N", [2, 2, 1]),            # pyrrole, explicit H
        ("c1ccoc1", "C4H4O", [2, 2]),                  # furan
        ("c1ccsc1", "C4H4S", [2, 2]),                  # thiophene
        ("c1ccc2[nH]c3ccccc3c2c1", "C12H9N", [2, 2, 2, 2, 1]),   # carbazole
        ("c1ccc2[nH]ccc2c1", "C8H7N", [1] * 7),        # indole
        ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", "C8H10N4O2", [3, 3, 3, 1]),  # caffeine
    ]

    def test_counts(self):
        for text, formula, envs in self.CASES:
            with self.subTest(smiles=text):
                mol = smiles.parse(text)
                self.assertEqual(mol.formula(), formula)
                self.assertEqual([e.count for e in mol.proton_environments()],
                                 envs)

    def test_caffeine_regression(self):
        """Adding the aromatic bond unconditionally gave C8H13N4O2."""
        mol = smiles.parse("Cn1cnc2c1c(=O)n(C)c(=O)n2C")
        self.assertEqual(mol.formula_counts()["H"], 10)


class TestSymmetry(unittest.TestCase):
    def test_equivalent_sites_are_merged(self):
        mol = smiles.parse("Cc1ccccc1")             # toluene
        envs = mol.proton_environments()
        self.assertEqual(len(envs), 4)
        methyl = envs[0]
        self.assertEqual(methyl.count, 3)
        self.assertEqual(len(methyl.atoms), 1)
        ortho = [e for e in envs if e.count == 2]
        self.assertEqual(len(ortho), 2)
        self.assertEqual(len(ortho[0].atoms), 2)     # two equivalent CH sites

    def test_tert_butyl_is_one_environment(self):
        mol = smiles.parse("CC(C)(C)c1ccccc1")
        self.assertEqual(mol.proton_environments()[0].count, 9)


class TestRingPerception(unittest.TestCase):
    CASES = [("c1ccccc1", 1, [6]), ("C1CC1", 1, [3]),
             ("c1ccc2ccccc2c1", 2, [6, 6]),
             ("c1ccc2[nH]c3ccccc3c2c1", 3, [5, 6, 6]),
             ("CCO", 0, []), ("c1ccccc1-c1ccccc1", 2, [6, 6])]

    def test_counts_and_sizes(self):
        for text, count, sizes in self.CASES:
            with self.subTest(smiles=text):
                mol = smiles.parse(text)
                self.assertEqual(mol.ring_count(), count)
                found = sorted(len(r) for r in depict.find_rings(mol))
                self.assertEqual(found, sorted(sizes))


class TestLayout(unittest.TestCase):
    STRUCTURES = ["c1ccccc1", "Cc1ccccc1", "c1ccc2ccccc2c1",
                  "c1ccc2[nH]c3ccccc3c2c1", "c1ccccc1-c1ccccc1",
                  "CC(=O)Oc1ccccc1C(=O)O", "c1ccc2cc3ccccc3cc2c1",
                  "CCOC(C)=O", "c1ccccc1CCc1ccccc1"]

    def test_bond_lengths_are_uniform(self):
        for text in self.STRUCTURES:
            with self.subTest(smiles=text):
                mol = smiles.parse(text)
                coords = depict.layout(mol)
                lengths = [
                    math.hypot(coords[b.a.index][0] - coords[b.b.index][0],
                               coords[b.a.index][1] - coords[b.b.index][1])
                    for b in mol.bonds]
                self.assertGreater(min(lengths), 0.85)
                self.assertLess(max(lengths), 1.2)

    def test_atoms_do_not_overlap(self):
        for text in self.STRUCTURES:
            with self.subTest(smiles=text):
                mol = smiles.parse(text)
                coords = depict.layout(mol)
                bonded = {frozenset((b.a.index, b.b.index)) for b in mol.bonds}
                indices = [a.index for a in mol.atoms]
                closest = min(
                    (math.hypot(coords[i][0] - coords[j][0],
                                coords[i][1] - coords[j][1])
                     for n, i in enumerate(indices) for j in indices[n + 1:]
                     if frozenset((i, j)) not in bonded), default=9.0)
                self.assertGreater(closest, 0.62)

    def test_rings_are_regular_polygons(self):
        mol = smiles.parse("c1ccc2[nH]c3ccccc3c2c1")
        coords = depict.layout(mol)
        for ring in depict.find_rings(mol):
            cx = sum(coords[a.index][0] for a in ring) / len(ring)
            cy = sum(coords[a.index][1] for a in ring) / len(ring)
            radii = [math.hypot(coords[a.index][0] - cx,
                                coords[a.index][1] - cy) for a in ring]
            self.assertLess(max(radii) - min(radii), 0.02)

    def test_svg_is_well_formed(self):
        import xml.etree.ElementTree as ET
        mol = smiles.parse("c1ccc2[nH]c3ccccc3c2c1")
        svg = depict.render_svg(mol, 300, 240)
        root = ET.fromstring(svg)
        self.assertTrue(root.tag.endswith("svg"))
        self.assertGreater(
            sum(1 for e in root.iter() if e.tag.endswith("line")), 10)


if __name__ == "__main__":
    unittest.main()
