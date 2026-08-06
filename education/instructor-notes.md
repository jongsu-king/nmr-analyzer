# Instructor notes

## Setup

Nothing needs installing. Copy the repository folder to each machine, or have
students clone it. Check once, on the machine they will actually use:

```bash
python3 run_tests.py            # 101 tests, should end "OK"
python3 education/make_datasets.py
python3 -m nmranalyzer
```

If `python3` is absent, hand out the pre-built application from the
[releases page](https://github.com/jongsu-king/nmr-analyzer/releases) — it
needs no Python at all. On Windows the first launch shows an unsigned-publisher
warning (*More info → Run anyway*); on macOS, right-click → Open. Do this once
before the class rather than in front of thirty students.

## Timing

| | Minutes |
|---|---|
| Demonstration: open a spectrum, zoom, integrate | 15 |
| Part 1 — integrals are ratios | 25 |
| Part 2 — phase errors | 30 |
| Part 3 — conversion | 20 |
| Part 4 — overlapping signals | 25 |
| Part 5 — the unknown | 40 |
| Discussion of Q11 | 15 |

Parts 1–3 stand alone if only ninety minutes are available. Part 4 is the one
to cut if time runs short; part 5 is the one worth protecting.

## Expected answers

**Part 1.** Ethyl acetate. 4.12 q (2H, *J* 7.1), 2.04 s (3H), 1.26 t (3H,
*J* 7.1). Measured relative values come out 1.00 : 1.49 : 1.50 against the
4.12 signal set to 2 H, i.e. 2 : 3 : 3.

- **Q1.** An integral is an area in arbitrary units; only ratios carry
  information. Declaring the 4.1 signal as 4 H doubles every other number,
  which is the point — the scale is arbitrary until one signal is fixed.
- **Q2.** Quartet and triplet, both *J* = 7.1 Hz. Equal couplings mean the two
  groups are coupled *to each other*: an ethyl group.
- **Q3.** `CCOC(C)=O`. The structure check matches all three regions.

**Part 2.** Ethanol, deliberately rotated by 35° of zero-order phase.

- As loaded, the CH₃:CH₂ ratio reads about **1.31** instead of 1.50 — roughly
  **13 % low**.
- After correcting PH0 (about +40 on the slider) it reads **1.50**.
- **Q4.** The error is proportional to the dispersion contribution, which has
  long tails; the narrower, taller signal loses relatively more of its area
  outside the integration limits.
- **Q5.** Auto Phase should land close to the manual answer. Students often
  find it satisfying that the automatic result is not magic — it minimises the
  negative area, which is what they were doing by eye.

**Part 3.** 60 % conversion. Measured: **60.4 %**.

- **Q6.** Comparing integrals gives a mole ratio only when the signals
  represent the same number of protons per molecule — here both are 3 H
  methyls. Otherwise each integral must first be divided by its proton count.
- **Q7.** The 7.3 ppm signal belongs to a fragment common to both species, so
  it carries no information about conversion. Including it as a third component
  changes nothing about the ratio of the two methyls, but if a student assigns
  it a proton count and includes it in the composition it dilutes both
  percentages. This is a good place to discuss what an internal standard is
  for.

**Part 4.** A 2 H and a 1 H singlet, 8 Hz apart, each 6 Hz wide.

- Manual boundaries typically give anything from 1.6:1 to 2.5:1 depending on
  where the student draws the line.
- Fitting recovers **2.00 : 1.00**.
- **Q9.** Integration assigns every point to exactly one region, so intensity
  in the overlap has to be given to one signal or the other. Fitting instead
  models both lineshapes and lets them share the overlap in the proportion
  their shapes require.

**Part 5.** 4-Methylacetophenone, `CC(=O)c1ccc(C)cc1`, C₉H₁₀O.

- 7.86 d (2H, *J* 8.1), 7.25 d (2H, *J* 8.1), 2.58 s (3H), 2.41 s (3H).
  Measured: 2.00, 1.99, 3.01, 3.02 H.
- DBE = 5: a benzene ring (4) plus the carbonyl (1).
- Two doublets with equal *J* and 2 H each is the *para*-disubstituted pattern.
- **Q10.** The correct structure matches all four regions with nothing left
  over. Wrong isomers of the same formula — 2-phenylpropanal, or
  3-phenylpropanal — match only two of four and leave four predicted
  environments unaccounted for. The discussion to draw out: the check tests
  *consistency*, and consistency is not proof. A structure that fails is
  ruled out; a structure that passes is only "not contradicted by this
  spectrum".

**Q11.** The intended answer: how it was phased and baseline-corrected, where
the integration limits were placed, and whether overlapping signals were
integrated or fitted. Students who have just watched a number move by 13 %
because of a slider tend to arrive at this on their own.

## Common difficulties

- **Students integrate before phasing.** This is worth allowing to happen in
  part 1 and then confronting in part 2, rather than warning against.
- **Dragging in the wrong mouse mode.** *Zoom* and *Integrate* look similar.
  Point at the mode selector in the demonstration.
- **Expecting exact integers.** 2.98 H is a good measurement. Discuss what
  size of deviation should worry them.
- **Trusting the structure check.** The point of Q10 is to puncture this.

## Assessment

The worksheet questions are markable as written, but for evidence of learning
the following pairs work as pre- and post-session items. Ask the pre-items
before any software is opened.

1. *A spectrum shows two signals with integrals 41.2 and 62.0. How many
   protons does each represent?* — tests objective 1; the correct response
   notes that the question is unanswerable without a reference.
2. *How many ¹H signals would you expect from p-xylene? From o-xylene?* —
   objective 2.
3. *A colleague reports a 3:2 integral ratio. Name two processing choices
   that could have changed that number, and say roughly how much.* —
   objective 4; the post-session version should draw on the 13 % they measured.
4. *Two signals overlap. Describe how you would obtain the area of each, and
   state the main source of uncertainty in your approach.* — objective 5.

A short confidence survey before and after ("how confident are you that you
could determine the ratio of two components in a mixture from its ¹H
spectrum?", five-point scale) gives a second measure and takes two minutes.
