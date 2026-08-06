# DRAFT — Journal of Chemical Education, Technology Report

**Working title:** Letting Students Break the Processing: A Dependency-Free
NMR Application for Teaching Spectral Interpretation

**Manuscript type:** Technology Report (3,000 word limit)
**Supporting Information:** required — the `education/` folder supplies it

> **This draft cannot be submitted yet.** JCE requires that the technology
> "should have been used with students and the results reported". Everything
> below is written; the sections marked ⚠ need data from an actual class.
> The session materials and the assessment items are ready to run.

---

## Abstract (≤ 150 words)

⚠ *Write last.* Must state: what the tool is, the instructional context it was
used in, how many students, what was measured, and what changed.

Skeleton: "Students commonly meet NMR as printed spectra or through vendor
software whose processing is opaque and licensed. We describe a free,
dependency-free desktop application … used in a [course] with [N] students …
Pre- and post-session responses showed [result] … The software and a
three-hour session with markable exercises are provided."

## Introduction (~600 words)

1. Interpreting a ¹H spectrum is a core skill, but students usually meet
   spectra *already processed*. They never see that the integral they are asked
   to trust is itself the product of choices.
2. The two things a student cannot do with a printed spectrum: change a
   processing parameter and watch the answer move, and read the code that
   produced the number.
3. Barriers to putting real software in students' hands: licence cost and seat
   limits, installation on managed machines, platform restrictions, and — for
   scriptable libraries such as nmrglue — the requirement to write code and
   install a scientific Python stack first.
4. State the gap and what this addresses.

*References to gather:* prior JCE work on NMR teaching software, on simulated
spectra in teaching, and on students' documented difficulties with integration
and with multiplet analysis. Search JCE for "NMR" + "software" and
"spectral interpretation" over the last ~15 years and cite the closest three
or four rather than a long list.

## The software (~500 words)

Keep this short — JCE readers want the pedagogy, not the architecture.

- Reads Bruker, JCAMP-DX and ACD `.esp`; processes raw FIDs.
- Peak picking, integration, multiplet analysis with *J*, lineshape fitting,
  2D contour display, and a structure-versus-integral consistency check from
  SMILES.
- **Runs with no installation and no dependencies** — the pedagogically
  relevant fact, because it is what lets it be used on managed teaching
  machines and on students' own laptops without support.
- **The processing is readable.** Emphasise this: a student who asks how the
  integral is computed can be shown the fifteen-line function. Include one
  short code excerpt as a figure — the integration routine is a good choice.
- Free, open source (MIT), Windows/macOS/Linux.

*Figure 1:* the main window with a spectrum integrated and assigned.
*(Available: `docs/01-main.png`.)*

## The session (~800 words)

Learning objectives, stated as in `education/README.md`.

Then walk through the five parts, spending most of the space on **part 2**,
which is the distinctive one:

> Students integrate a deliberately mis-phased spectrum of ethanol and obtain a
> CH₃:CH₂ ratio of about 1.31 rather than 1.50 — a 13 % error. They then
> correct the phase by hand and re-measure. The number moves under their
> control. This is the moment the session is built around: it converts
> "phasing matters" from an assertion into a measurement.

Also worth a paragraph: **part 4**, where a boundary drawn by hand between two
overlapping signals gives answers ranging from about 1.6:1 to 2.5:1, while
fitting recovers 2.00:1.00. And **part 5**, where a wrong isomer of the correct
molecular formula is rejected by the structure check, which opens the
discussion that consistency is not proof.

*Figure 2:* the same spectrum before and after phase correction, with the two
integral values.
*Figure 3:* the overlapping pair, with the fitted components drawn over it.

Note that the datasets are simulated so that every answer is exact and a class
can be marked objectively, and say so plainly — a reviewer will ask.

## ⚠ Implementation and assessment (~700 words) — NEEDS CLASS DATA

This is the section that decides acceptance. Required:

- Course, level, class size, whether the session was compulsory, and how it sat
  in the syllabus.
- How the session was run: individually or in pairs, machines used, staffing.
- **Evidence of learning.** The four pre/post items in
  `education/instructor-notes.md` are written for this. Report the change in
  the proportion of correct responses, with the *n* for each item.
- The confidence survey before and after, reported as a distribution.
- Whether students could operate the software without help — count the
  interventions needed. For a "technology" report this is itself a result.
- What went wrong. Reviewers trust a paper that reports the difficulties.

Minimum that would be defensible: one cohort, *n* ≳ 20, four pre/post items,
the confidence survey, and an honest account of where students struggled. Two
cohorts is much stronger.

## Limitations (~200 words)

State plainly:

- The datasets are simulated; students have not processed data they acquired.
  Describe how the session extends to real spectra, which the software reads.
- The chemical-shift prediction is additivity and is unreliable for multiply
  substituted carbons; it is a hint for assignment, not an answer.
- The structure check tests consistency, not identity — and note that Q10 is
  designed to teach exactly this.
- 2D display is included but was not used in this session.

## Associated content

**Supporting Information** — supply the `education/` folder as one archive plus
editable copies:

| Item | Source |
|---|---|
| Student worksheet | `education/student-worksheet.md` → convert to .docx |
| Instructor notes with answers | `education/instructor-notes.md` → .docx |
| Dataset generator | `education/make_datasets.py` |
| The five spectra (JCAMP-DX) | `education/datasets/*.jdx` |
| Answer key | `education/datasets/ANSWERS.md` |
| Pre/post assessment items | in the instructor notes |

JCE requires SI in editable format, so convert the two Markdown documents to
Word before submission.

**Software availability** — <https://github.com/jongsu-king/nmr-analyzer>,
MIT licence, with pre-built applications for Windows, macOS and Linux. Archive
the version used in the class on Zenodo and cite that DOI, so the paper points
at exactly what the students ran.

---

## Before submitting: checklist

- [ ] Run the session with a class ⚠ **the blocker**
- [ ] Collect pre/post responses and the confidence survey
- [ ] Obtain ethics approval or exemption for the student data, if required
- [ ] Fill in the Implementation and assessment section
- [ ] Search and cite the closest prior JCE work
- [ ] Convert worksheet and instructor notes to .docx for SI
- [ ] Archive the exact software version on Zenodo, cite the DOI
- [ ] Cut to 3,000 words
