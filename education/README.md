# Teaching materials

A three-hour session in which students process and interpret ¹H NMR spectra
themselves rather than reading printed ones.

Everything here is self-contained: `make_datasets.py` generates the spectra, so
their composition — and therefore the correct answer to every question — is
known exactly and the worksheet can be marked objectively.

```bash
python3 education/make_datasets.py     # writes education/datasets/
python3 -m nmranalyzer                 # the software
```

## Learning objectives

By the end of the session a student should be able to:

1. **Relate integral ratios to proton counts**, and state why an integral is
   only meaningful relative to another.
2. **Predict the number of ¹H signals from a structure** using molecular
   symmetry, and explain a case where the naive count is wrong.
3. **Extract a coupling constant** from line spacings and classify a multiplet.
4. **Explain how phase and baseline errors distort an integral**, having
   measured the distortion rather than been told about it.
5. **Judge whether a spectrum is consistent with a proposed structure**, and
   say what the evidence does and does not establish.

Objective 4 is the one that is hard to teach from printed spectra: the student
has to be able to *break* the processing and see the number change.

## Why this software for this session

- **Nothing to install.** It runs from a folder on a locked-down teaching
  machine, with no administrator rights and no package manager.
- **The processing is not a black box.** A student who asks how an integral is
  computed can be shown the function that computes it, in about fifteen lines
  of ordinary Python. `dsp.py` and `analysis.py` are written to be read.
- **No licence limits.** Every student can have it, including on their own
  laptop, during and after the course.

## Files

| File | What it is |
|---|---|
| `student-worksheet.md` | The exercises, in order |
| `instructor-notes.md` | Timing, setup, expected answers, common difficulties |
| `make_datasets.py` | Generates the spectra |
| `datasets/ANSWERS.md` | The exact composition of each spectrum |

## A note on the datasets

The spectra are simulated. That is a deliberate choice: it makes every answer
exact, so a student who reports 2.9 H knows the target was 3, and an instructor
can mark a class of thirty consistently. Real spectra should follow — the
software reads Bruker, JCAMP-DX and ACD `.esp` — but the simulated set removes
the ambiguity from the first pass.
