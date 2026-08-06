#!/usr/bin/env python3
"""How well does an integral set discriminate between candidate structures?

For each compound the true structure is checked against its own assigned
spectrum (does it pass?), and then against decoys — other compounds of the
same molecular formula (do they fail?).  Sensitivity and specificity follow.

    python3 benchmark/run.py [--limit N] [--exchangeable]

Writes ``benchmark/results.json`` and prints a summary.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from nmranalyzer import smiles, verify                            # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DECOYS_PER_COMPOUND = 5


def load():
    with open(os.path.join(HERE, "dataset.json"), encoding="utf-8") as fh:
        return json.load(fh)


def by_formula(entries):
    groups = {}
    for entry in entries:
        groups.setdefault(entry["formula"], []).append(entry)
    return groups


def evaluate(limit=0, ignore_exchangeable=False, seed=0):
    entries = load()
    groups = by_formula(entries)
    rng = random.Random(seed)

    # Only compounds whose formula is shared can be given a same-formula decoy.
    testable = [e for e in entries if len(groups[e["formula"]]) > 1]
    if limit:
        testable = testable[:limit]

    parsed = {}

    def molecule(text):
        if text not in parsed:
            try:
                parsed[text] = smiles.parse(text)
            except Exception:
                parsed[text] = None
        return parsed[text]

    true_pass = true_total = 0
    decoy_pass = decoy_total = 0
    per_compound = []
    started = time.time()

    for index, entry in enumerate(testable):
        signals = [(s, n) for s, n in entry["signals"]]
        mol = molecule(entry["smiles"])
        if mol is None:
            continue

        result = verify.check(mol, signals, ignore_exchangeable)
        true_total += 1
        accepted = result.consistent
        true_pass += accepted

        pool = [e for e in groups[entry["formula"]]
                if e["smiles"] != entry["smiles"]]
        rng.shuffle(pool)
        decoy_results = []
        for decoy in pool[:DECOYS_PER_COMPOUND]:
            dmol = molecule(decoy["smiles"])
            if dmol is None:
                continue
            dres = verify.check(dmol, signals, ignore_exchangeable)
            decoy_total += 1
            decoy_pass += dres.consistent
            decoy_results.append(dres.consistent)

        per_compound.append({
            "name": entry["name"][:60],
            "formula": entry["formula"],
            "signals": len(signals),
            "protons": entry["protons"],
            "true_accepted": accepted,
            "decoys": len(decoy_results),
            "decoys_accepted": sum(decoy_results),
            "summary": result.summary(),
        })

        if index and index % 100 == 0:
            print("   ... %d/%d" % (index, len(testable)), flush=True)

    elapsed = time.time() - started
    return {
        "compounds": true_total,
        "true_accepted": true_pass,
        "sensitivity": true_pass / true_total if true_total else 0.0,
        "decoys": decoy_total,
        "decoys_accepted": decoy_pass,
        "specificity": 1.0 - (decoy_pass / decoy_total) if decoy_total else 0.0,
        "ignore_exchangeable": ignore_exchangeable,
        "seconds": round(elapsed, 1),
        "per_compound": per_compound,
    }


def report(result):
    print("\n" + "=" * 62)
    print("  compounds tested            %6d" % result["compounds"])
    print("  true structure accepted     %6d   sensitivity %.3f"
          % (result["true_accepted"], result["sensitivity"]))
    print("  same-formula decoys tested  %6d" % result["decoys"])
    print("  decoys wrongly accepted     %6d   specificity %.3f"
          % (result["decoys_accepted"], result["specificity"]))
    print("  elapsed                     %6.1f s" % result["seconds"])
    print("=" * 62)

    rows = result["per_compound"]
    missed = [r for r in rows if not r["true_accepted"]]
    if missed:
        print("\n  why true structures were rejected (top reasons):")
        reasons = {}
        for r in missed:
            reasons[r["summary"]] = reasons.get(r["summary"], 0) + 1
        for reason, count in sorted(reasons.items(), key=lambda x: -x[1])[:6]:
            print("    %5d  %s" % (count, reason))

    by_signals = {}
    for r in rows:
        bucket = min(r["signals"], 10)
        hit, total = by_signals.get(bucket, (0, 0))
        by_signals[bucket] = (hit + r["true_accepted"], total + 1)
    print("\n  sensitivity by number of signals:")
    for bucket in sorted(by_signals):
        hit, total = by_signals[bucket]
        label = "%d" % bucket if bucket < 10 else "10+"
        print("    %-4s n=%-5d %.3f" % (label, total, hit / total))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--exchangeable", action="store_true",
                    help="ignore OH/NH protons, as an experimentalist often must")
    args = ap.parse_args()

    result = evaluate(args.limit, args.exchangeable)
    report(result)
    target = os.path.join(HERE, "results.json")
    with open(target, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=1)
    print("\nwrote %s" % target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
