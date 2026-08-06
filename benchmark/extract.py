#!/usr/bin/env python3
"""Turn nmrshiftdb2 into a benchmark of structures and their proton counts.

Each record in the NMReData file carries a structure and a ¹H spectrum whose
lines are assigned to individual hydrogen atoms.  Collapsing that assignment
gives exactly what an experimentalist obtains by integrating: a list of
(chemical shift, number of protons).  That is the input the structure check
consumes, so the benchmark exercises the real code path on real assignments.

    python3 benchmark/extract.py path/to/nmrshiftdb2.nmredata.sd

Writes ``benchmark/dataset.json``.
"""

from __future__ import annotations

import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from nmranalyzer import smiles                                    # noqa: E402

MERGE_PPM = 0.02        # lines closer than this would be one integral
MIN_SIGNALS = 3         # fewer than this carries too little information
MAX_HEAVY = 40          # keep the structures tractable
MAX_RECORDS = 0         # 0 = no limit


def records(path):
    """Yield each SDF record as text."""
    buf = []
    with open(path, encoding="latin-1", errors="replace") as fh:
        for line in fh:
            if line.startswith("$$$$"):
                yield "".join(buf)
                buf = []
            else:
                buf.append(line)


def tag_block(record, name):
    """The lines of one ``> <TAG>`` block."""
    match = re.search(r"> <%s>\s*\n" % re.escape(name), record)
    if not match:
        return []
    out = []
    for line in record[match.end():].splitlines():
        if line.startswith("> <") or not line.strip():
            break
        out.append(line.strip().rstrip("\\").strip())
    return out


def parse_assignment(record):
    """``label -> number of atoms`` from the NMREDATA_ASSIGNMENT block."""
    counts = {}
    for line in tag_block(record, "NMREDATA_ASSIGNMENT"):
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 3:
            continue
        label = parts[0]
        atoms = [p for p in parts[2:] if p]
        counts[label] = len(atoms)
    return counts


def parse_proton_lines(record):
    """``[(shift, label)]`` from the NMREDATA_1D_1H block."""
    out = []
    for line in tag_block(record, "NMREDATA_1D_1H"):
        if "=" in line and not line.startswith(("-", ".")) and "," not in line:
            continue                      # Larmor=, Spectrum_Location= etc.
        parts = [p.strip() for p in line.split(",")]
        try:
            shift = float(parts[0])
        except (ValueError, IndexError):
            continue
        label = None
        for p in parts[1:]:
            if p.upper().startswith("L="):
                label = p[2:].strip()
        out.append((shift, label))
    return out


def collapse(lines, assignment):
    """Merge nearby lines into the signals an integration would report."""
    signals = []
    for shift, label in sorted(lines):
        protons = assignment.get(label, 1) if label else 1
        if signals and abs(shift - signals[-1][0]) <= MERGE_PPM:
            prev_shift, prev_protons = signals[-1]
            total = prev_protons + protons
            # intensity-weighted centre, as an integral's centre would be
            centre = (prev_shift * prev_protons + shift * protons) / total
            signals[-1] = (centre, total)
        else:
            signals.append((shift, protons))
    return signals


def build(path, limit=MAX_RECORDS):
    out = []
    seen = 0
    stats = {"no 1H": 0, "no SMILES": 0, "unparsable": 0,
             "too few signals": 0, "too large": 0,
             "proton count mismatch": 0, "kept": 0}

    for record in records(path):
        seen += 1
        if limit and seen > limit:
            break
        if "NMREDATA_1D_1H" not in record:
            stats["no 1H"] += 1
            continue
        match = re.search(r"> <NMREDATA_SMILES>\s*\n(.+)", record)
        if not match:
            stats["no SMILES"] += 1
            continue
        text = match.group(1).strip().rstrip("\\")
        try:
            mol = smiles.parse(text)
        except Exception:
            stats["unparsable"] += 1
            continue
        if len(mol.atoms) > MAX_HEAVY * 3:
            stats["too large"] += 1
            continue

        signals = collapse(parse_proton_lines(record), parse_assignment(record))
        if len(signals) < MIN_SIGNALS:
            stats["too few signals"] += 1
            continue

        # The assignment must account for the protons the formula says exist,
        # otherwise the "integral" is not comparable with a prediction.
        formula_h = mol.formula_counts().get("H", 0)
        assigned_h = sum(n for _s, n in signals)
        if not formula_h or abs(assigned_h - formula_h) > 0:
            stats["proton count mismatch"] += 1
            continue

        name = tag_block(record, "CHEMNAME")
        solvent = tag_block(record, "NMREDATA_SOLVENT")
        out.append({
            "smiles": text,
            "formula": mol.formula(),
            "name": name[0] if name else "",
            "solvent": solvent[0] if solvent else "",
            "signals": [[round(s, 4), n] for s, n in signals],
            "protons": formula_h,
        })
        stats["kept"] += 1

    return out, stats, seen


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    source = sys.argv[1]
    data, stats, seen = build(source)

    here = os.path.dirname(os.path.abspath(__file__))
    target = os.path.join(here, "dataset.json")
    with open(target, "w", encoding="utf-8") as fh:
        json.dump(data, fh)

    print("records read      : %d" % seen)
    for key, value in stats.items():
        print("  %-22s %6d" % (key, value))
    print("\nwrote %s  (%d compounds, %.1f MB)"
          % (target, len(data), os.path.getsize(target) / 1e6))

    formulas = {}
    for entry in data:
        formulas.setdefault(entry["formula"], []).append(entry)
    shared = {f: v for f, v in formulas.items() if len(v) > 1}
    print("distinct formulae : %d" % len(formulas))
    print("formulae shared by >1 compound (decoy pool): %d, covering %d compounds"
          % (len(shared), sum(len(v) for v in shared.values())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
