# How this work was produced

Stated plainly, because it matters for how the result should be read.

## Roles

- **Andrei Pluzhnik** — the research question, every decision about what to pursue, approval of
  every claim's status, and all publishing. The interventions that changed the direction of the
  work are recorded in the parent repository's night report.
- **Claude Opus 5, via Claude Code** — the implementation: writing and running the computations,
  deriving and checking the identities, drafting the text, and maintaining the record of what
  failed.

The model was not permitted to decide scientific truth. Every claim's status is set by a
deterministic gate over artefacts, not by the model's judgement of plausibility.

## What that means for the reader

Everything labelled **PROVED** rests on an exact computation that a reader can rerun: rational
arithmetic with `python-flint`'s `fmpq`, or certified interval arithmetic with `arb`/`acb`. No
floating-point number enters a comparison anywhere in the chain. Each script prints its own verdict.

Everything labelled **VERIFIED** is a measurement over a stated finite set of cases, with the
number of cases given. It is evidence, not proof, and is labelled so in every place it appears.

Everything labelled **CONJECTURED** is exactly that.

## The record of failures is part of the work

Claims that were stated and then refuted are kept beside the results they were meant to support,
with the computation that killed them. `LIMITATIONS.md` lists them. The parent repository's
`docs/ERRATA.md` carries the full set, including errors of instrument rather than of object — a
wrong assumed degree, two truncation-garbage coefficients, a precision hypothesis that reruns
disproved.

This is deliberate. A result produced this way is only as trustworthy as its account of what went
wrong, and that account is easier to write while the work is happening than afterwards.

## Reproducibility

- Exact rational arithmetic: `flint.fmpq`, `flint.fmpq_poly`, `flint.fmpz_poly`.
- Certified interval arithmetic: `flint.arb`, `flint.acb`, at explicitly set precision.
- Every interpolation of a polynomial asserts itself against exact values at twelve fresh points
  before it is used. That guard exists because it was once absent and a wrong degree passed
  unnoticed.
