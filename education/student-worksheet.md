# ¹H NMR: processing and interpretation

Open the software with `python3 -m nmranalyzer`, then `File > Open Data Files`
and choose a spectrum from `education/datasets`.

Useful controls: drag to zoom, right-drag to pan, double-click for the whole
spectrum, and `Integrate` mode to drag across a signal.

---

## 1 — Integrals are ratios (dataset 01)

Open **`01-ethyl-acetate.jdx`**. This is an unknown, *A*.

1. Integrate each of the three signals. Write down the raw integral values.

   | Signal (ppm) | Raw integral |
   |---|---|
   | | |

2. The raw numbers are meaningless on their own. Select the signal near
   4.1 ppm, type `2` into *Set H* and press **Assign**.

   Now record the *Relative* column for all three:

   | Signal (ppm) | Protons |
   |---|---|
   | | |

3. **Q1.** Why was it necessary to declare one signal as 2 H before the others
   could be read as proton counts? What would have changed if you had declared
   it as 4 H instead? Try it.

4. **Q2.** The *Multiplet* and *J (Hz)* columns are filled in automatically.
   Write down the multiplicity and coupling constant of the signal at 4.1 ppm
   and of the one at 1.26 ppm. What does it mean that they are equal?

5. **Q3.** Propose a structure for *A* consistent with your integrals and
   multiplicities. Check it: `Tools > Structure`, enter your SMILES, press
   *Draw*, then *Check against spectrum*.

---

## 2 — Phase errors change the answer (dataset 02)

Open **`02-ethanol-misphased.jdx`**. This spectrum has been left badly phased.

1. Without correcting anything, integrate the signal near 3.7 ppm and the one
   near 1.2 ppm. Assign the 3.7 signal as 2 H. What proton count does the
   1.2 ppm signal come out as?

   Measured: ____________  (it should be 3 H)

2. Look at the baseline on either side of each peak. Describe its shape. A
   correctly phased peak sits on a flat baseline and is positive on both sides.

3. Correct the phase: drag the **PH0** slider until every peak is upright and
   the baseline either side is flat. Re-integrate.

   Measured after phasing: ____________

4. **Q4.** By what percentage did the apparent proton count change? Which is
   larger — the error in the tall signal or in the short one — and why?

5. **Q5.** Press **Auto Phase** and compare with your manual result. Note the
   PH0 and PH1 values it chose.

---

## 3 — How complete is the reaction? (dataset 03)

Open **`03-conversion.jdx`**. A reaction converts a starting material into a
product. The methyl singlet of the starting material is at 2.10 ppm; that of
the product is at 2.35 ppm. Both are 3 H.

1. Integrate both singlets. Assign **3** protons to each.
2. Read the *Report* tab.

   Conversion: ____________ %

3. **Q6.** The two signals you integrated are both 3 H, so why is it valid to
   compare their integrals directly to get a mole ratio? Under what
   circumstance would comparing integrals *not* give a mole ratio?

4. **Q7.** There is a third signal near 7.3 ppm belonging to a part of the
   molecule that does not change in the reaction. Integrate it as well. Does
   including it change the conversion you calculate? Should it?

---

## 4 — When integration is not enough (dataset 04)

Open **`04-overlapping.jdx`**. Two signals overlap.

1. Zoom in on the region near 7.4 ppm. How many signals can you see?
2. Integrate across the whole lump as one region. Record the value.
3. Now try to integrate the two components separately by dragging two adjacent
   regions. Do it twice, choosing the boundary in a different place each time.

   Attempt 1: ______ : ______   Attempt 2: ______ : ______

4. **Q8.** How much did your answer depend on where you put the boundary?
5. Select the whole region and press **Fit lines**. Read the Fit tab.

   Fitted ratio: ______ : ______

6. **Q9.** Explain in one or two sentences why fitting can separate these two
   signals when drawing a boundary cannot.

---

## 5 — Identify the unknown (dataset 05)

Open **`05-unknown.jdx`**. Compound *B* has molecular formula **C₉H₁₀O**.

1. Integrate every signal and determine the proton count of each, using the
   formula as your total.
2. Record multiplicities and coupling constants.
3. Calculate the degrees of unsaturation from the formula. What does the value
   suggest?
4. Propose a structure. Enter it under `Tools > Structure` and use *Check
   against spectrum*.
5. **Q10.** The check reports whether each integral matches a predicted
   environment in both size *and* chemical shift. Deliberately enter a *wrong*
   isomer with the same formula and see what the check says. What does this
   tell you about how much confidence the check deserves?

---

## Reflection

**Q11.** In question 2 you changed a processing parameter and the measured
proton count changed. If a published paper reports an integral as "2.05 H",
what would you now want to know about how the spectrum was processed?
