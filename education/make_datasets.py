#!/usr/bin/env python3
"""Generate the teaching datasets, each with an answer the instructor knows.

Spectra are synthesised from Lorentzian lines at chosen shifts and written as
JCAMP-DX, which every NMR program reads.  Because the composition is chosen
rather than measured, the correct answer to every exercise is known exactly,
which is what makes the worksheet markable.

    python3 education/make_datasets.py
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SF = 400.0                 # observe frequency, MHz
SW_PPM = 12.0              # 12 ppm wide
POINTS = 16384
NOISE = 0.0015             # fraction of the tallest line


def lorentzian(x, centre, height, width_hz, phase_deg=0.0):
    """Lorentzian line, optionally rotated out of pure absorption.

    A zero-order phase error mixes in the dispersion lineshape, so building
    the line that way gives a spectrum that is mis-phased in exactly the sense
    a real one is -- and that a phase correction therefore fixes.
    """
    u = 2.0 * (x - centre) * SF / width_hz
    absorption = 1.0 / (1.0 + u * u)
    dispersion = u / (1.0 + u * u)
    if not phase_deg:
        return height * absorption
    angle = math.radians(phase_deg)
    return height * (absorption * math.cos(angle) + dispersion * math.sin(angle))


def multiplet(centre, protons, pattern, j_hz, width_hz=1.2):
    """Lines of a first-order multiplet, with binomial intensities."""
    n = {"s": 1, "d": 2, "t": 3, "q": 4, "quint": 5, "sept": 7}[pattern]
    weights = [1]
    for _ in range(n - 1):
        weights = [1] + [weights[i] + weights[i + 1]
                         for i in range(len(weights) - 1)] + [1]
    total = float(sum(weights))
    spacing = j_hz / SF
    first = centre - spacing * (n - 1) / 2.0
    return [(first + i * spacing, protons * w / total, width_hz)
            for i, w in enumerate(weights)]


def build(lines, seed=1, phase_deg=0.0):
    import random
    rng = random.Random(seed)
    top = 10.0 / SW_PPM * POINTS      # arbitrary intensity scale
    data = []
    for i in range(POINTS):
        ppm = SW_PPM - i * SW_PPM / POINTS
        value = sum(lorentzian(ppm, c, h, w, phase_deg) for c, h, w in lines)
        data.append(value * top)
    peak = max(data) or 1.0
    return [v + rng.gauss(0.0, NOISE * peak) for v in data]


def write_jcamp(path, title, data, solvent="CDCl3"):
    """Write a minimal but valid JCAMP-DX 1H spectrum (AFFN, X++(Y..Y))."""
    first_x, last_x = SW_PPM, SW_PPM - (POINTS - 1) * SW_PPM / POINTS
    scale = max(abs(v) for v in data) / 1e6
    ints = [int(round(v / scale)) for v in data]
    with open(path, "w", encoding="ascii") as fh:
        fh.write("##TITLE=%s\n##JCAMP-DX=4.24\n##DATA TYPE=NMR SPECTRUM\n" % title)
        fh.write("##ORIGIN=NMR Analyzer teaching set\n##OWNER=public domain\n")
        fh.write("##.OBSERVE FREQUENCY=%.4f\n##.OBSERVE NUCLEUS=^1H\n" % SF)
        fh.write("##.SOLVENT NAME=%s\n" % solvent)
        fh.write("##XUNITS=PPM\n##YUNITS=ARBITRARY UNITS\n")
        fh.write("##FIRSTX=%.6f\n##LASTX=%.6f\n##NPOINTS=%d\n" %
                 (first_x, last_x, POINTS))
        fh.write("##XFACTOR=1.0\n##YFACTOR=%.8E\n" % scale)
        fh.write("##FIRSTY=%.6E\n" % data[0])
        fh.write("##XYDATA=(X++(Y..Y))\n")
        per_line = 8
        for i in range(0, POINTS, per_line):
            x = SW_PPM - i * SW_PPM / POINTS
            fh.write("%.4f %s\n" % (x, " ".join(str(v)
                                                for v in ints[i:i + per_line])))
        fh.write("##END=\n")


# --- the datasets ----------------------------------------------------------
# Each entry: filename, title, the lines, and the answer the worksheet expects.

DATASETS = [
    dict(
        name="01-ethyl-acetate",
        title="Unknown A",
        answer="ethyl acetate, CCOC(C)=O.  4.12 q (2H, J 7.1), "
               "2.04 s (3H), 1.26 t (3H, J 7.1)",
        lines=(multiplet(4.12, 2, "q", 7.1)
               + multiplet(2.04, 3, "s", 0.0)
               + multiplet(1.26, 3, "t", 7.1)),
    ),
    dict(
        name="02-ethanol-misphased",
        title="Ethanol (needs phasing)",
        answer="ethanol, CCO.  3.69 q (2H), 2.60 s (1H, OH), 1.19 t (3H). "
               "Deliberately left unphased so the integrals are wrong "
               "until the student corrects it.",
        lines=(multiplet(3.69, 2, "q", 7.0)
               + multiplet(2.60, 1, "s", 0.0)
               + multiplet(1.19, 3, "t", 7.0)),
        phase=35.0,
    ),
    dict(
        name="03-conversion",
        title="Reaction mixture, aliquot",
        answer="60 % conversion.  Starting material CH3 singlet 2.10 (3H) at "
               "40 % of the mixture; product CH3 singlet 2.35 (3H) at 60 %.",
        lines=(multiplet(2.10, 3 * 0.40, "s", 0.0)
               + multiplet(2.35, 3 * 0.60, "s", 0.0)
               + multiplet(7.30, 5 * 1.00, "s", 0.0, width_hz=2.0)),
    ),
    dict(
        name="04-overlapping",
        title="Overlapping signals",
        answer="two singlets that overlap: 2H at 7.41 and 1H at 7.39, each "
               "6 Hz wide.  They merge into one lump with a shoulder, so any "
               "boundary drawn between them is arbitrary; fitting recovers "
               "2:1.",
        lines=(multiplet(7.41, 2, "s", 0.0, width_hz=6.0)
               + multiplet(7.39, 1, "s", 0.0, width_hz=6.0)),
    ),
    dict(
        name="05-unknown",
        title="Unknown B",
        answer="4-methylacetophenone, CC(=O)c1ccc(C)cc1.  "
               "7.86 d (2H), 7.25 d (2H), 2.58 s (3H), 2.41 s (3H)",
        lines=(multiplet(7.86, 2, "d", 8.1)
               + multiplet(7.25, 2, "d", 8.1)
               + multiplet(2.58, 3, "s", 0.0)
               + multiplet(2.41, 3, "s", 0.0)),
    ),
]


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "datasets")
    os.makedirs(out, exist_ok=True)

    answers = ["# Answer key\n",
               "Generated by `make_datasets.py`; the composition is chosen, "
               "so these are exact.\n"]
    for spec in DATASETS:
        data = build(spec["lines"], phase_deg=spec.get("phase", 0.0))
        path = os.path.join(out, spec["name"] + ".jdx")
        write_jcamp(path, spec["title"], data)
        print("  %-28s %6.0f KB" % (spec["name"] + ".jdx",
                                    os.path.getsize(path) / 1024))
        answers.append("\n## %s\n\n%s\n" % (spec["name"], spec["answer"]))

    with open(os.path.join(out, "ANSWERS.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(answers))
    print("  ANSWERS.md")


if __name__ == "__main__":
    main()
