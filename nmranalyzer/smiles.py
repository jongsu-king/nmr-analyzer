"""A small SMILES parser, plus the molecular arithmetic an NMR user needs.

Enough of the language for ordinary organic structures: the organic subset,
bracket atoms with charge/isotope/explicit hydrogens, aromatic lowercase
atoms, branches, ring closures (including ``%nn``) and disconnected parts.
Stereochemistry is parsed and then ignored, because a flat depiction and a
proton count do not depend on it.

The point of the symmetry pass is that topologically equivalent protons give
one NMR signal, so counting equivalence classes predicts how many signals a
structure should show.
"""

from __future__ import annotations

import re

# Standard atomic weights, IUPAC 2021, for the elements that turn up in
# ordinary organic chemistry.
ATOMIC_WEIGHT = {
    "H": 1.008, "D": 2.014, "B": 10.81, "C": 12.011, "N": 14.007,
    "O": 15.999, "F": 18.998, "Na": 22.990, "Mg": 24.305, "Al": 26.982,
    "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "K": 39.098,
    "Ca": 40.078, "Fe": 55.845, "Cu": 63.546, "Zn": 65.38, "Br": 79.904,
    "I": 126.904, "Sn": 118.71, "Pd": 106.42, "Pb": 207.2,
}

# Valences used for implicit hydrogen counting.
DEFAULT_VALENCE = {
    "B": (3,), "C": (4,), "N": (3, 5), "O": (2,), "P": (3, 5),
    "S": (2, 4, 6), "F": (1,), "Cl": (1,), "Br": (1,), "I": (1,),
}

ORGANIC_SUBSET = ["Cl", "Br", "B", "C", "N", "O", "P", "S", "F", "I",
                  "b", "c", "n", "o", "p", "s"]

_BRACKET = re.compile(
    r"\[(\d*)([A-Za-z][a-z]?)(@{0,2}|@TH[12]|@AL[12])?(H\d*)?([+-]\d*|[+-]+)?(?::\d+)?\]")


class Atom:
    def __init__(self, symbol, aromatic=False, charge=0, isotope=None,
                 explicit_h=None, index=0, chirality=""):
        self.symbol = symbol            # capitalised element symbol
        self.aromatic = aromatic
        self.charge = charge
        self.isotope = isotope
        self.explicit_h = explicit_h    # set only for bracket atoms
        self.chirality = chirality      # "@" / "@@" if the SMILES declared it
        self.index = index
        self.bonds = []                 # list of Bond
        self.n_hydrogens = 0            # filled in by _assign_hydrogens

    @property
    def degree(self):
        return len(self.bonds)

    def neighbours(self):
        return [b.other(self) for b in self.bonds]

    def __repr__(self):
        return "Atom(%d %s%s)" % (self.index, self.symbol,
                                  "*" if self.aromatic else "")


class Bond:
    def __init__(self, a, b, order=1, aromatic=False):
        self.a = a
        self.b = b
        self.order = order
        self.aromatic = aromatic

    def other(self, atom):
        return self.b if atom is self.a else self.a

    def __repr__(self):
        return "Bond(%d-%d, %s)" % (self.a.index, self.b.index,
                                    "ar" if self.aromatic else self.order)


class Molecule:
    def __init__(self, atoms, bonds, smiles=""):
        self.atoms = atoms
        self.bonds = bonds
        self.smiles = smiles

    # -- composition --------------------------------------------------------

    def formula_counts(self):
        counts = {}
        for atom in self.atoms:
            counts[atom.symbol] = counts.get(atom.symbol, 0) + 1
        hydrogens = sum(a.n_hydrogens for a in self.atoms)
        if hydrogens:
            counts["H"] = counts.get("H", 0) + hydrogens
        return counts

    def formula(self):
        """Hill notation: C first, then H, then everything else alphabetically."""
        counts = self.formula_counts()
        parts = []
        for symbol in ("C", "H"):
            if counts.get(symbol):
                parts.append(symbol + (str(counts[symbol])
                                       if counts[symbol] > 1 else ""))
        for symbol in sorted(k for k in counts if k not in ("C", "H")):
            parts.append(symbol + (str(counts[symbol])
                                   if counts[symbol] > 1 else ""))
        return "".join(parts)

    def molecular_weight(self):
        total = 0.0
        for symbol, count in self.formula_counts().items():
            total += ATOMIC_WEIGHT.get(symbol, 0.0) * count
        return total

    def dbe(self):
        """Double bond equivalents: rings plus pi bonds.

        DBE = (2*C + 2 + N - H - X) / 2, with halogens counted as hydrogens.
        """
        counts = self.formula_counts()
        carbon = counts.get("C", 0)
        nitrogen = counts.get("N", 0)
        hydrogen = counts.get("H", 0)
        halogen = sum(counts.get(x, 0) for x in ("F", "Cl", "Br", "I"))
        return (2 * carbon + 2 + nitrogen - hydrogen - halogen) / 2.0

    def ring_count(self):
        """Independent rings, from the cycle rank of each connected part."""
        return len(self.bonds) - len(self.atoms) + self._component_count()

    def _component_count(self):
        seen = set()
        components = 0
        for atom in self.atoms:
            if atom.index in seen:
                continue
            components += 1
            stack = [atom]
            while stack:
                current = stack.pop()
                if current.index in seen:
                    continue
                seen.add(current.index)
                stack.extend(current.neighbours())
        return components

    # -- symmetry -----------------------------------------------------------

    def equivalence_classes(self):
        """Group topologically equivalent atoms (Morgan-style refinement).

        Two atoms end up in the same class when no amount of looking outward
        through the bond graph can tell them apart, which is exactly the
        condition for their protons to give one NMR signal.
        """
        ranks = {}
        for atom in self.atoms:
            ranks[atom.index] = (atom.symbol, atom.aromatic, atom.charge,
                                 atom.degree, atom.n_hydrogens)

        previous = -1
        for _ in range(len(self.atoms) + 2):
            distinct = len(set(ranks.values()))
            if distinct == previous:
                break
            previous = distinct
            refined = {}
            for atom in self.atoms:
                neighbourhood = sorted(
                    (str(ranks[b.other(atom).index]),
                     "ar" if b.aromatic else str(b.order))
                    for b in atom.bonds)
                refined[atom.index] = (str(ranks[atom.index]),
                                       tuple(neighbourhood))
            # Re-key to compact tuples so the strings do not grow without bound.
            ordered = sorted(set(refined.values()))
            lookup = {value: position for position, value in enumerate(ordered)}
            ranks = {index: lookup[value] for index, value in refined.items()}

        groups = {}
        for atom in self.atoms:
            groups.setdefault(ranks[atom.index], []).append(atom)
        return list(groups.values())

    def proton_environments(self, split_diastereotopic=True):
        """Distinct proton environments, most protons first.

        The two protons of a diastereotopic CH2 are reported separately, since
        they genuinely give different signals; pass
        ``split_diastereotopic=False`` for the purely constitutional count.
        """
        splittable = ({a.index for a in self.diastereotopic_carbons()}
                      if split_diastereotopic else set())
        branches = (self.diastereotopic_branches()
                    if split_diastereotopic else [])
        out = []
        for group in self.equivalence_classes():
            per_atom = group[0].n_hydrogens
            if not per_atom:
                continue
            members = frozenset(a.index for a in group)
            if members in branches:
                # Two identical branches that a stereocentre pulls apart.
                for atom, tag in zip(group, ("a", "b")):
                    out.append(ProtonEnvironment([atom], per_atom,
                                                 diastereotopic=tag))
            elif group[0].index in splittable:
                for tag in ("a", "b"):
                    out.append(ProtonEnvironment(group, len(group),
                                                 diastereotopic=tag))
            else:
                out.append(ProtonEnvironment(group, per_atom * len(group)))
        out.sort(key=lambda env: (-env.count, env.label))
        return out


    def carbon_environments(self):
        """Distinct carbon environments, the 13C analogue of the above.

        Quaternary carbons are included: unlike 1H they still give a signal,
        just a weak one.
        """
        out = []
        for group in self.equivalence_classes():
            if group[0].symbol != "C":
                continue
            out.append(CarbonEnvironment(group))
        out.sort(key=lambda env: (-len(env.atoms), env.label))
        return out


    # -- topicity -----------------------------------------------------------

    def _neighbour_classes(self):
        """Map each atom index to its equivalence-class id."""
        lookup = {}
        for cid, group in enumerate(self.equivalence_classes()):
            for atom in group:
                lookup[atom.index] = cid
        return lookup

    def stereocentres(self):
        """Carbons that make the molecule chiral.

        Either declared in the SMILES with ``@``/``@@``, or found by the
        constitutional test: four substituents that are all different, counting
        an implicit hydrogen as one of them.
        """
        classes = self._neighbour_classes()
        found = []
        for atom in self.atoms:
            if atom.symbol != "C" or atom.aromatic:
                continue
            if atom.chirality:
                found.append(atom)
                continue
            if atom.degree + atom.n_hydrogens != 4 or atom.n_hydrogens > 1:
                continue
            seen = [classes[n.index] for n in atom.neighbours()]
            if len(set(seen)) == len(seen) == (4 - atom.n_hydrogens):
                found.append(atom)
        return found

    def prochiral_centres(self):
        """Atoms carrying exactly two constitutionally identical branches.

        Replacing one branch and then the other gives either enantiomers or
        diastereomers, never the same molecule, which is what makes the two
        branches -- and the protons of an attached CH2 -- distinguishable.
        """
        classes = self._neighbour_classes()
        out = []
        for atom in self.atoms:
            if atom.aromatic:
                continue
            neighbours = atom.neighbours()
            if len(neighbours) + atom.n_hydrogens != 4:
                continue
            counts = {}
            for n in neighbours:
                counts[classes[n.index]] = counts.get(classes[n.index], 0) + 1
            paired = [cid for cid, k in counts.items() if k == 2]
            if len(paired) != 1:
                continue
            rest = [str(classes[n.index]) for n in neighbours
                    if classes[n.index] != paired[0]]
            rest += ["H"] * atom.n_hydrogens
            if len(set(rest)) == len(rest):
                out.append((atom, paired[0]))
        return out

    def diastereotopic_branches(self):
        """Pairs of identical branches that a stereocentre makes inequivalent.

        The classic case is an isopropyl group next to a stereocentre: its two
        methyls give two separate doublets rather than one six-proton one.
        """
        if not self.stereocentres():
            return []
        classes = self._neighbour_classes()
        pairs = []
        for atom, cid in self.prochiral_centres():
            branch = [n.index for n in atom.neighbours()
                      if classes[n.index] == cid]
            if len(branch) == 2:
                pairs.append(frozenset(branch))
        return pairs

    def diastereotopic_carbons(self):
        """CH2 groups whose two protons are inequivalent.

        A CH2 is *prochiral* when its other two substituents differ.  Its
        protons then become diastereotopic -- genuinely different shifts --
        once the molecule contains a stereocentre, or when the carbon sits in
        a ring, where the two faces are not the same.

        This catches the cases that matter in practice.  It does not catch
        pseudo-asymmetric centres such as glycerol, where the protons are
        diastereotopic without any stereocentre being present.
        """
        classes = self._neighbour_classes()
        has_stereocentre = bool(self.stereocentres())
        ring_atoms = self._ring_atom_indices()
        prochiral = {a.index for a, _cid in self.prochiral_centres()}

        out = []
        for atom in self.atoms:
            if atom.symbol != "C" or atom.n_hydrogens != 2:
                continue
            heavy = [n for n in atom.neighbours()]
            if len(heavy) != 2:
                continue
            if classes[heavy[0].index] == classes[heavy[1].index]:
                continue                     # not prochiral itself
            attached_prochiral = any(n.index in prochiral for n in heavy)
            if has_stereocentre or atom.index in ring_atoms or attached_prochiral:
                out.append(atom)
        return out

    def _ring_atom_indices(self):
        """Atoms that lie on a cycle, via a bridge-free spanning test."""
        indices = set()
        seen = set()
        order = {}
        low = {}
        counter = [0]
        bridges = set()

        def visit(atom, parent_bond):
            order[atom.index] = low[atom.index] = counter[0]
            counter[0] += 1
            for bond in atom.bonds:
                if bond is parent_bond:
                    continue
                other = bond.other(atom)
                if other.index not in order:
                    visit(other, bond)
                    low[atom.index] = min(low[atom.index], low[other.index])
                    if low[other.index] > order[atom.index]:
                        bridges.add(id(bond))
                else:
                    low[atom.index] = min(low[atom.index], order[other.index])

        import sys as _sys
        limit = _sys.getrecursionlimit()
        _sys.setrecursionlimit(max(limit, len(self.atoms) * 4 + 100))
        try:
            for atom in self.atoms:
                if atom.index not in order:
                    visit(atom, None)
        finally:
            _sys.setrecursionlimit(limit)

        for bond in self.bonds:
            if id(bond) not in bridges:
                indices.add(bond.a.index)
                indices.add(bond.b.index)
        return indices


class CarbonEnvironment:
    def __init__(self, atoms):
        self.atoms = atoms

    @property
    def carrier(self):
        return self.atoms[0]

    @property
    def count(self):
        """Carbons in this environment (they all give one line)."""
        return len(self.atoms)

    @property
    def label(self):
        atom = self.carrier
        hydrogens = atom.n_hydrogens
        kind = {0: "quaternary C", 1: "CH", 2: "CH2", 3: "CH3"}.get(
            hydrogens, "C-H%d" % hydrogens)
        if atom.aromatic:
            kind = "aromatic " + ("C" if hydrogens == 0 else "CH")
        return kind

    def describe(self):
        sites = len(self.atoms)
        where = " (%d equivalent)" % sites if sites > 1 else ""
        return "%s%s" % (self.label, where)


class ProtonEnvironment:
    def __init__(self, atoms, count, diastereotopic=""):
        self.atoms = atoms          # equivalent heavy atoms carrying the H
        self.count = count          # total protons in this environment
        self.diastereotopic = diastereotopic   # "a"/"b" for a split CH2

    @property
    def carrier(self):
        return self.atoms[0]

    @property
    def label(self):
        atom = self.carrier
        kind = atom.symbol
        if atom.aromatic:
            kind = "aromatic " + kind
        if self.diastereotopic:
            return "%s-H%d (H%s, diastereotopic)" % (kind, atom.n_hydrogens,
                                                     self.diastereotopic)
        return "%s-H%s" % (kind, "" if atom.n_hydrogens == 1 else
                           str(atom.n_hydrogens))

    def describe(self):
        sites = len(self.atoms)
        where = " (%d equivalent sites)" % sites if sites > 1 else ""
        return "%dH  %s%s" % (self.count, self.label, where)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


class SmilesError(ValueError):
    pass


def parse(text):
    """Parse a SMILES string into a :class:`Molecule`."""
    text = (text or "").strip()
    if not text:
        raise SmilesError("empty SMILES")

    atoms = []
    bonds = []
    branch_stack = []
    ring_bonds = {}             # ring number -> (atom, pending bond order)
    previous = None
    pending_order = None
    pending_aromatic = False

    i = 0
    n = len(text)
    while i < n:
        ch = text[i]

        if ch == "(":
            if previous is None:
                raise SmilesError("branch opens before any atom")
            branch_stack.append(previous)
            i += 1
            continue
        if ch == ")":
            if not branch_stack:
                raise SmilesError("unbalanced ')'")
            previous = branch_stack.pop()
            i += 1
            continue
        if ch in "-=#$:":
            pending_order = {"-": 1, "=": 2, "#": 3, "$": 4, ":": 1}[ch]
            pending_aromatic = ch == ":"
            i += 1
            continue
        if ch in "/\\":                 # cis/trans marker, no effect here
            pending_order = 1
            i += 1
            continue
        if ch == ".":
            previous = None
            pending_order = None
            i += 1
            continue

        # ring closure
        if ch.isdigit() or ch == "%":
            if ch == "%":
                match = re.match(r"%(\d\d)", text[i:])
                if not match:
                    raise SmilesError("bad %% ring label at position %d" % i)
                number = int(match.group(1))
                i += 3
            else:
                number = int(ch)
                i += 1
            if previous is None:
                raise SmilesError("ring closure before any atom")
            if number in ring_bonds:
                partner, order, aromatic = ring_bonds.pop(number)
                order = pending_order or order or 1
                aromatic = (aromatic or pending_aromatic
                            or (partner.aromatic and previous.aromatic))
                bonds.append(_link(partner, previous, order, aromatic))
            else:
                ring_bonds[number] = (previous, pending_order, pending_aromatic)
            pending_order = None
            pending_aromatic = False
            continue

        # atom
        atom, consumed = _read_atom(text, i, len(atoms))
        i += consumed
        atoms.append(atom)
        if previous is not None:
            order = pending_order or 1
            aromatic = (pending_aromatic
                        or (previous.aromatic and atom.aromatic
                            and pending_order is None))
            bonds.append(_link(previous, atom, order, aromatic))
        previous = atom
        pending_order = None
        pending_aromatic = False

    if branch_stack:
        raise SmilesError("unbalanced '('")
    if ring_bonds:
        raise SmilesError("unclosed ring label(s): %s"
                          % ", ".join(str(k) for k in sorted(ring_bonds)))

    molecule = Molecule(atoms, bonds, smiles=text)
    _assign_hydrogens(molecule)
    return molecule


def _link(a, b, order, aromatic):
    bond = Bond(a, b, order=order, aromatic=aromatic)
    a.bonds.append(bond)
    b.bonds.append(bond)
    return bond


def _read_atom(text, i, index):
    if text[i] == "[":
        end = text.find("]", i)
        if end < 0:
            raise SmilesError("unclosed '[' at position %d" % i)
        match = _BRACKET.match(text[i:end + 1])
        if not match:
            raise SmilesError("cannot read bracket atom %r" % text[i:end + 1])
        isotope, symbol, chirality, hydrogens, charge = match.groups()

        aromatic = symbol[0].islower()
        element = symbol.capitalize()

        explicit_h = 0
        if hydrogens:
            explicit_h = int(hydrogens[1:]) if len(hydrogens) > 1 else 1

        total_charge = 0
        if charge:
            if charge in ("+", "-") or set(charge) == {charge[0]}:
                total_charge = len(charge) * (1 if charge[0] == "+" else -1)
            else:
                total_charge = int(charge)

        atom = Atom(element, aromatic=aromatic, charge=total_charge,
                    isotope=int(isotope) if isotope else None,
                    explicit_h=explicit_h, index=index,
                    chirality=chirality or "")
        return atom, end + 1 - i

    for symbol in ORGANIC_SUBSET:
        if text.startswith(symbol, i):
            # "Cl"/"Br" must not swallow the C of e.g. "C" followed by "l"
            aromatic = symbol[0].islower()
            atom = Atom(symbol.capitalize(), aromatic=aromatic, index=index)
            return atom, len(symbol)

    if text[i] == "*":
        return Atom("*", index=index), 1
    raise SmilesError("unexpected character %r at position %d" % (text[i], i))


def _assign_hydrogens(molecule):
    """Fill in implicit hydrogens from standard valences."""
    for atom in molecule.atoms:
        if atom.explicit_h is not None:
            atom.n_hydrogens = atom.explicit_h
            continue

        used = 0.0
        for bond in atom.bonds:
            used += 1 if bond.aromatic else bond.order

        # An aromatic atom that carries a formal double bond in the Kekule
        # structure uses one more bond than its sigma count.  Carbon always
        # does.  Nitrogen only does in the pyridine sense: once it has three
        # sigma neighbours it must be donating its lone pair instead, as in
        # N-methylpyrrole, and adding the extra bond would invent a hydrogen.
        # Aromatic O and S are always lone-pair donors.
        if atom.aromatic:
            if atom.symbol == "C":
                used += 1
            elif atom.symbol in ("N", "P") and atom.degree < 3:
                used += 1

        valences = DEFAULT_VALENCE.get(atom.symbol)
        if not valences:
            atom.n_hydrogens = 0
            continue

        target = valences[0]
        for candidate in valences:
            if candidate >= used:
                target = candidate
                break
        else:
            target = valences[-1]

        # A charge changes how many bonds the atom wants.
        if atom.symbol == "N" and atom.charge > 0:
            target += atom.charge
        elif atom.charge < 0:
            target += atom.charge
        elif atom.charge > 0 and atom.symbol != "N":
            target -= atom.charge

        atom.n_hydrogens = max(0, int(round(target - used)))
