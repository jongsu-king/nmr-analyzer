"""Rough chemical-shift estimates from connectivity.

These are the classical additivity rules, not a calculation: a base value for
the kind of environment plus increments for what is attached.  Aromatic rings
follow the substituent tables that every textbook prints; aliphatic carbons get
alpha and attenuated beta corrections.

Accuracy is roughly +/- 0.3 ppm for 1H and +/- 5 ppm for 13C on ordinary
molecules, and much worse on anything strained, charged or conjugated in an
unusual way.  Every estimate therefore carries a window, and the callers are
expected to treat it as a hint for assigning peaks rather than a prediction.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# 1H
# ---------------------------------------------------------------------------

BENZENE_H = 7.26

# substituent -> (ortho, meta, para) increment on the ring protons
AROMATIC_INCREMENTS = {
    "alkyl": (-0.17, -0.09, -0.18),
    "OH": (-0.56, -0.12, -0.45),
    "OR": (-0.48, -0.09, -0.44),
    "OC=O": (-0.25, 0.03, -0.13),
    "NH2": (-0.75, -0.25, -0.65),
    "NR2": (-0.66, -0.18, -0.67),
    "NC=O": (0.12, -0.07, -0.28),
    "F": (-0.26, 0.00, -0.20),
    "Cl": (0.02, -0.06, -0.04),
    "Br": (0.22, -0.13, -0.03),
    "I": (0.40, -0.26, -0.03),
    "NO2": (0.95, 0.26, 0.38),
    "CHO": (0.58, 0.21, 0.27),
    "COOH": (0.85, 0.18, 0.27),
    "COR": (0.62, 0.14, 0.21),
    "COOR": (0.74, 0.07, 0.20),
    "CN": (0.36, 0.18, 0.28),
    "aryl": (0.18, 0.00, 0.08),
    "C=C": (0.06, -0.03, -0.10),
}

ALIPHATIC_BASE = {1: 1.50, 2: 1.20, 3: 0.86}     # CH, CH2, CH3

# effect on a proton of a group attached to the same carbon
ALPHA_H = {
    "OH": 2.30, "OR": 2.35, "OC=O": 3.00, "NH2": 1.55, "NR2": 1.55,
    "NC=O": 1.95, "F": 3.20, "Cl": 2.20, "Br": 2.10, "I": 1.80,
    "NO2": 3.00, "CHO": 1.15, "COOH": 1.05, "COR": 1.05, "COOR": 1.05,
    "CN": 1.00, "aryl": 1.40, "C=C": 0.80, "SR": 1.25, "alkyl": 0.00,
}
BETA_FACTOR = 0.22       # a group one bond further away counts for much less

# environments where additivity is meaningless and only a range can be given
SPECIAL_H = {
    "aldehyde": (9.70, 0.4),
    "COOH": (11.50, 1.5),
    "OH": (3.00, 2.0),
    "phenol OH": (5.50, 2.5),
    "NH amine": (1.80, 1.5),
    "NH amide": (7.00, 1.5),
    "NH aromatic": (8.50, 2.0),
    "SH": (1.60, 0.8),
    "alkene": (5.30, 0.7),
    "alkyne": (2.40, 0.4),
}

# ---------------------------------------------------------------------------
# 13C
# ---------------------------------------------------------------------------

CARBON_BASE = {
    "CH3": (14.0, 8.0), "CH2": (25.0, 10.0), "CH": (32.0, 10.0),
    "Cq": (35.0, 12.0),
    "aromatic CH": (128.5, 6.0), "aromatic C": (135.0, 12.0),
    "alkene CH": (125.0, 10.0), "alkene C": (135.0, 12.0),
    "alkyne": (75.0, 10.0),
    "nitrile": (118.0, 5.0),
    "aldehyde": (192.0, 5.0), "ketone": (205.0, 8.0),
    "acid": (178.0, 5.0), "ester": (171.0, 5.0), "amide": (170.0, 6.0),
}
CARBON_ALPHA = {
    "OH": 35.0, "OR": 40.0, "OC=O": 40.0, "NH2": 20.0, "NR2": 22.0,
    "NC=O": 18.0, "F": 60.0, "Cl": 25.0, "Br": 10.0, "I": -10.0,
    "NO2": 55.0, "CN": 2.0, "aryl": 10.0, "C=C": 8.0, "SR": 10.0,
    "CHO": 6.0, "COOH": 8.0, "COR": 8.0, "COOR": 8.0, "alkyl": 0.0,
}


class Estimate:
    """A predicted shift with an honest window around it."""

    def __init__(self, value, window, basis=""):
        self.value = value
        self.window = window
        self.basis = basis

    @property
    def low(self):
        return self.value - self.window

    @property
    def high(self):
        return self.value + self.window

    def contains(self, shift):
        return self.low <= shift <= self.high

    def __repr__(self):
        return "%.2f +/- %.2f" % (self.value, self.window)

    def text(self):
        return "%.2f (%.1f-%.1f)" % (self.value, self.low, self.high)


# ---------------------------------------------------------------------------
# Classifying what is attached
# ---------------------------------------------------------------------------


def _is_carbonyl(atom):
    return atom.symbol == "C" and any(
        b.order == 2 and b.other(atom).symbol == "O" for b in atom.bonds)


def substituent_kind(atom, came_from):
    """Name the group that ``atom`` starts, seen from ``came_from``."""
    symbol = atom.symbol

    if symbol in ("F", "Cl", "Br", "I"):
        return symbol

    if symbol == "O":
        others = [n for n in atom.neighbours() if n is not came_from]
        if atom.n_hydrogens:
            return "OH"
        if any(_is_carbonyl(n) for n in others):
            return "OC=O"
        return "OR"

    if symbol == "N":
        others = [n for n in atom.neighbours() if n is not came_from]
        if any(_is_carbonyl(n) for n in others):
            return "NC=O"
        if any(b.order == 3 for b in atom.bonds):
            return "CN"
        oxygens = sum(1 for n in atom.neighbours() if n.symbol == "O")
        if oxygens >= 2:
            return "NO2"
        return "NH2" if atom.n_hydrogens >= 1 else "NR2"

    if symbol == "S":
        return "SR"

    if symbol == "C":
        if any(b.order == 3 for b in atom.bonds):
            partner = [n for n in atom.neighbours() if n is not came_from]
            if any(n.symbol == "N" for n in partner):
                return "CN"
            return "C=C"
        if _is_carbonyl(atom):
            oxygens = [n for n in atom.neighbours()
                       if n.symbol == "O" and n is not came_from]
            single = [o for o in oxygens
                      if not any(b.order == 2 and b.other(o) is atom
                                 for b in o.bonds)]
            if atom.n_hydrogens:
                return "CHO"
            for oxygen in single:
                return "COOH" if oxygen.n_hydrogens else "COOR"
            if any(n.symbol == "N" for n in atom.neighbours()):
                return "NC=O"
            return "COR"
        if atom.aromatic:
            return "aryl"
        if any(b.order == 2 for b in atom.bonds):
            return "C=C"
        return "alkyl"
    return "alkyl"


# ---------------------------------------------------------------------------
# Aromatic ring geometry
# ---------------------------------------------------------------------------


def _benzene_ring(atom, rings):
    """The six-membered all-aromatic ring this atom sits in, if any."""
    for ring in rings:
        if len(ring) == 6 and atom in ring and all(a.aromatic for a in ring):
            return ring
    return None


def _aromatic_h(atom, rings):
    ring = _benzene_ring(atom, rings)
    if ring is None:
        # Fused or heteroaromatic: additivity does not transfer, so widen.
        return Estimate(7.60, 1.0, "aromatic CH, ring not a simple benzene")

    order = list(ring)
    here = order.index(atom)
    total = BENZENE_H
    described = []
    for position, member in enumerate(order):
        if member is atom:
            continue
        outside = [n for n in member.neighbours() if n not in order]
        if not outside:
            continue
        separation = min((position - here) % 6, (here - position) % 6)
        if separation not in (1, 2, 3):
            continue
        for neighbour in outside:
            kind = substituent_kind(neighbour, member)
            increments = AROMATIC_INCREMENTS.get(kind)
            if increments is None:
                continue
            total += increments[separation - 1]
            described.append("%s %s" % (["ortho", "meta", "para"][separation - 1],
                                        kind))
    basis = "benzene 7.26" + (" + " + ", ".join(described) if described else "")
    return Estimate(total, 0.35, basis)


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------


def proton_shift(atom, rings):
    """Estimate the shift of the protons on ``atom``."""
    if atom.symbol == "O":
        neighbours = [n for n in atom.neighbours()]
        if any(n.aromatic for n in neighbours):
            value, window = SPECIAL_H["phenol OH"]
            return Estimate(value, window, "phenol OH, concentration dependent")
        if any(_is_carbonyl(n) for n in neighbours):
            value, window = SPECIAL_H["COOH"]
            return Estimate(value, window, "carboxylic acid OH")
        value, window = SPECIAL_H["OH"]
        return Estimate(value, window, "alcohol OH, concentration dependent")

    if atom.symbol == "N":
        if atom.aromatic:
            value, window = SPECIAL_H["NH aromatic"]
            return Estimate(value, window, "aromatic NH")
        if any(_is_carbonyl(n) for n in atom.neighbours()):
            value, window = SPECIAL_H["NH amide"]
            return Estimate(value, window, "amide NH")
        value, window = SPECIAL_H["NH amine"]
        return Estimate(value, window, "amine NH, concentration dependent")

    if atom.symbol == "S":
        value, window = SPECIAL_H["SH"]
        return Estimate(value, window, "thiol SH")

    if atom.symbol != "C":
        return Estimate(3.0, 3.0, "no rule for %s-H" % atom.symbol)

    if atom.aromatic:
        return _aromatic_h(atom, rings)

    if _is_carbonyl(atom) and atom.n_hydrogens:
        value, window = SPECIAL_H["aldehyde"]
        return Estimate(value, window, "aldehyde CH")

    if any(b.order == 3 for b in atom.bonds):
        value, window = SPECIAL_H["alkyne"]
        return Estimate(value, window, "alkyne CH")

    if any(b.order == 2 for b in atom.bonds):
        value, window = SPECIAL_H["alkene"]
        return Estimate(value, window, "alkene CH")

    # sp3: base for CH/CH2/CH3 plus alpha and beta group effects
    base = ALIPHATIC_BASE.get(atom.n_hydrogens, 1.50)
    total = base
    described = []
    for neighbour in atom.neighbours():
        kind = substituent_kind(neighbour, atom)
        alpha = ALPHA_H.get(kind, 0.0)
        if alpha:
            total += alpha
            described.append("alpha %s" % kind)
        # one bond further out
        for second in neighbour.neighbours():
            if second is atom:
                continue
            beta_kind = substituent_kind(second, neighbour)
            beta = ALPHA_H.get(beta_kind, 0.0) * BETA_FACTOR
            if beta >= 0.15:
                total += beta
                described.append("beta %s" % beta_kind)
    label = {3: "CH3", 2: "CH2", 1: "CH"}.get(atom.n_hydrogens, "CH")
    basis = "%s %.2f" % (label, base)
    if described:
        basis += " + " + ", ".join(described)
    return Estimate(total, 0.4, basis)


def carbon_shift(atom, rings):
    """Estimate the 13C shift of ``atom``."""
    if atom.symbol != "C":
        return None

    if _is_carbonyl(atom):
        kind = substituent_kind(atom, None)
        mapping = {"CHO": "aldehyde", "COOH": "acid",
                   "COOR": "ester", "NC=O": "amide", "COR": "ketone"}
        base, window = CARBON_BASE[mapping.get(kind, "ketone")]
        return Estimate(base, window, mapping.get(kind, "ketone"))

    if atom.aromatic:
        key = "aromatic CH" if atom.n_hydrogens else "aromatic C"
        base, window = CARBON_BASE[key]
        return Estimate(base, window, key)

    if any(b.order == 3 for b in atom.bonds):
        if any(n.symbol == "N" for n in atom.neighbours()):
            base, window = CARBON_BASE["nitrile"]
            return Estimate(base, window, "nitrile")
        base, window = CARBON_BASE["alkyne"]
        return Estimate(base, window, "alkyne")

    if any(b.order == 2 for b in atom.bonds):
        key = "alkene CH" if atom.n_hydrogens else "alkene C"
        base, window = CARBON_BASE[key]
        return Estimate(base, window, key)

    key = {3: "CH3", 2: "CH2", 1: "CH"}.get(atom.n_hydrogens, "Cq")
    base, window = CARBON_BASE[key]
    total = base
    described = []
    for neighbour in atom.neighbours():
        kind = substituent_kind(neighbour, atom)
        shift = CARBON_ALPHA.get(kind, 0.0)
        if shift:
            total += shift
            described.append(kind)
    basis = "%s %.0f" % (key, base)
    if described:
        basis += " + " + ", ".join(described)
    return Estimate(total, window, basis)


def predict_proton_environments(mol):
    """``[(environment, Estimate), ...]`` ordered by descending shift."""
    from . import depict
    rings = depict.find_rings(mol)
    out = [(env, proton_shift(env.carrier, rings))
           for env in mol.proton_environments()]
    out.sort(key=lambda pair: -pair[1].value)
    return out


def predict_carbon_environments(mol):
    from . import depict
    rings = depict.find_rings(mol)
    out = [(env, carbon_shift(env.carrier, rings))
           for env in mol.carbon_environments()]
    out.sort(key=lambda pair: -pair[1].value)
    return out
