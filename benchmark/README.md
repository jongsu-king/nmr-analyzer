# Benchmark: how far can proton counts alone verify a structure?

The structure check in this program scales measured integrals to a proposed
molecular formula and asks whether each signal matches a predicted proton
environment in both size and chemical shift. This directory measures how well
that works on real assigned spectra.

## Data

[nmrshiftdb2](https://nmrshiftdb.nmr.uni-koeln.de/) is distributed as an
NMReData SD file in which each ¹H line is assigned to individual hydrogen
atoms. Collapsing that assignment gives exactly what integration yields — a
list of (shift, number of protons) — so the benchmark feeds the real code path
with real assignments.

```bash
python3 benchmark/extract.py nmrshiftdb2.nmredata.sd    # -> dataset.json
python3 benchmark/run.py                                # -> results.json
```

Of 64,723 records, 1,301 survive filtering. The strict requirement is that the
assignment accounts for **every** proton the formula implies; most entries list
only resolved signals, and an incomplete set cannot be scaled to the formula.
581 of the survivors share a molecular formula with another compound and can
therefore be given a same-formula decoy.

## What was measured

**As a binary verdict** — demanding that every signal match and no predicted
environment be left over — the check has specificity 0.988 but sensitivity
0.060. It rejects nearly everything, correct structures included, so it is not
usable as confirmation.

The limiting factor is shift prediction, not proton counting. Per signal, 39 %
fail on shift alone against 13 % on size alone, and where a shift fails the
median miss is 1.4 ppm *outside* the quoted window. Additivity rules are poor
on the natural products that dominate nmrshiftdb2.

**As a ranking** — scoring candidates by total mismatch and asking whether the
true structure outranks a decoy — it does considerably better:

| | pairs | AUC |
|---|---|---|
| all same-formula decoys | 1156 | 0.660 |
| excluding decoys with an identical environment multiset | 519 | **0.799** |

637 of the pairs were discarded on the second row because the decoy predicts
exactly the same set of proton counts as the true structure. No integral set
can separate those, so scoring them measures nothing.

## Reading this honestly

The useful claim is the ranking one: given several candidate structures of the
same formula, integrals plus additivity shifts put the right one first about
four times in five. The binary check should be read as a filter that rules
structures out, not one that confirms them — which is how the interface
words its verdict.

Improving this means better shift prediction, not better matching.
