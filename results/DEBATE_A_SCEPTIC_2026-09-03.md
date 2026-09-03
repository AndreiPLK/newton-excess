# Sceptic report: Theorem A (`M(n,t) > 4/5`, centred-square spectrum)

**Time (measured with `date`):** Thu Sep 3 16:04:59 RDT 2026 (started) -- Thu Sep 3 16:17:19 RDT 2026 (last check)
**Role:** sceptic. Did not write any of the certified scripts. Attacked by importing the actual
functions and calling them at specific points (no heavy full runs launched), by independent
closed-form/brute-force cross-checks, and by reading every load-bearing constant back to its
derivation or its validation report.
**Scripts:** `C:\Users\user\AppData\Local\Temp\claude\C--Users-user-ScienceBro\ac66e2dc-eaec-44ab-b12d-0b61f841fa72\scratchpad\scepA\`
(`check_worst_box.py`, `check_Ci_beyond60.py`, `check_float_line.py`, `check_float_vmax.py`,
`brute_check.py`, `row_shape.py`, `t1_limit.py`). No repository file modified except this report.
**What I did not duplicate:** the verifier agent's reproduction of load-bearing facts
(`results/DEBATE_A_VERIFIER_2026-09-03.md`) -- I read the existing validation/referee reports in
full and targeted what they had *not* checked, rather than re-deriving what they already did.

## Context acknowledged

Theorem A already carries an unusual amount of adversarial scrutiny for this codebase:
`VALIDATION_SPARSE_CERTIFICATE_FULL_2026-09-02.md`, `VALIDATION_DENSE_CERTIFICATE_A_2026-09-02.md`,
and `VALIDATION_DENSE_CERTIFICATE_B_2026-09-02.md` (with two addenda) each did a genuine
separate-formulation re-derivation and found real bugs (a spurious `r^{2k}` factor, a
non-analytic Cauchy-tail object, an `R4` term missing an `O(delta_0)` piece, a `1.06` multiplier
invalid at `j=2`, a `2.7/3^j` constant that a circle-Cauchy estimate cannot actually prove). Every
one of those bugs turned out to be numerically harmless by wide margins once corrected, and the
referee review (`REFEREE_REVIEW_RELEASE_2026-09-03.md`) separately caught seven blocking
packaging/bookkeeping defects, mostly on Theorem B and the release metadata, plus one item that
does touch Theorem A (M5, float in the coverage lemma's final comparison). **I confirmed M5 is
already fixed**: the live `theorem.py` (read below) computes `imp1/imp2/imp3` and the final
`assert` entirely in `arb`, not `float`, matching commit `f96854f`'s stated fix. I also confirmed
the referee's M1/M2 (missing run logs) are fixed: `results/sparse_certificate_full_2026-09-03.txt`
and `results/theorem_full_2026-09-03.txt` now exist, and a fresh run I launched reproduces both
byte-for-byte in content.

Given that, my job was to find what neither of those rounds looked for. Findings below, in
destruction order.

---

## 1. `dense_certificate_b.py`: the corner's pass/fail boolean does not check the corner's claim

**File:line:** `projects/qg-bootstrap/release/scripts/dense_certificate_b.py:599`

```python
ok = worst_S is not None and worst_S > 0 and float(loss_corner) < 0.1
```

`worst_S` comes only from the per-band bisection loop (`V` from `6.485e-7` up to `VMAX=0.17`,
18 bands). The corner `V <= V_last ~ 1.3e-6` (which the docstring says is covered "by
monotonicity, see text": `R_inner(sigma) <= R_inner(sigma_last)(sigma_last/sigma)^11`, both pieces
decreasing as `V -> 0`) is handled by a **separate, un-gated** block: it computes
`loss_corner` and prints `S.evaluate(...).lower()` for the reader, but **the actual comparison
`S(corner) - loss_corner > 0` is never computed or asserted anywhere in the code.** The `ok`
variable only requires `loss_corner < 0.1` -- an arbitrary constant with no algebraic relation to
`S`'s value at the corner. If a future change made `loss_corner` creep toward `0.1` while `S`
collapsed toward `0.1` from above for an unrelated reason, `ok` would stay `True` while the
inequality it is supposed to certify would be false.

**Measured (ran the live script, not a duplicate):**

```
corner V <= 1.297e-06: with L = sigma/2 = 731.7, remainder 6.40e-29, loss <= 4.57e-16
  (S >= 0.3352 there; ...)
```

`loss_corner = 4.57e-16` against `S = 0.3352` -- a margin of **15 orders of magnitude**. The
per-band values climb monotonically toward `S(0,zeta) = 176/525 = 0.335238...` as `V -> 0`
(reproduced: `0.3206, 0.3278, ..., 0.3352`), and `loss_corner` shrinks geometrically as `V` shrinks
(the validator's addenda derived the exact mechanism: `R_inner ~ sigma^{-11}`, `sigma ~ V^{-1/2}`).
So the numeric conclusion is safe by a wide margin under the *stated* mechanism -- but the code
does not check that mechanism's conclusion, only an unrelated magic threshold.

**Verdict: SLACK-ONLY.** Not fatal as currently run (verified numerically, both the printed `S`
value and the printed `loss_corner`), but it is a real gap between what the code asserts and what
the docstring claims is proved, and it is exactly the failure shape from this morning's Theorem B
debate: a step whose correctness is *believed* (here: asserted in prose, "the corner is covered by
monotonicity") rather than *checked* by the code that produces `exit 0`. Recommend: replace the
`float(loss_corner) < 0.1` line with an explicit
`S.evaluate(arb(hi/2, hi/2), zb).lower() - loss_corner > 0` (or equivalent `arb` comparison), so
that the printed sentence and the exit code can never diverge.

---

## 2. The 627-rung ladder (piece A) is the only piece of Theorem A with no independent
   (separately-formulated) validation -- and it is the largest single computation (7.06 CPU-h)

**Files:** `projects/qg-bootstrap/release/scripts/ladder_fast.py`,
`projects/qg-bootstrap/lab/ladder_log_merge.py`,
`projects/qg-bootstrap/results/m_ladder_log_2026-09-02_full.txt`.

Sparse, dense-(a) and dense-(b) each have a `VALIDATION_*.md` report built from a **separately
formulated** re-implementation (different variable classes, different derivations, independent
constants) that found and fixed real bugs. The ladder has none. What runs on every `theorem.py`
invocation is:

- `lab/ladder_log_merge.py 627`: checks the log is complete (rungs 1..627 present once), that no
  shift is `None`, that the logged degree equals the closed-form `8i+10`, and that the logged shift
  equals the closed-form prediction `(i+5)//2` (odd `i`) or `(i+2)//2` (even `i`). **This is a
  format/consistency check on the log file, not a re-derivation of the polynomial positivity
  certificate.** It would not catch a systematic bug in how `W_i(n)` itself is constructed, only a
  transcription/completeness failure of the log.
- `release/scripts/ladder_fast.py 12 1` ("sample re-verification of rungs 1..12"): I confirmed by
  diffing the header line (`"   INTEGER LADDER   i    deg W    c_M    covers n <=    time/index"`)
  against the archival range logs (e.g. `results/m_ladder_log_2026-09-02_r214_330.txt`) that this
  is **the identical script** that produced the original 7-hour run, just re-run on a 12-rung
  slice. It is not an independent formulation, so a bug shared by the construction of `W_i`
  (Faulhaber power sums, the Newton recursion for `e_j`, the Taylor-shift-then-check-nonnegative
  test) would reproduce identically on the sample and pass.

**What the underlying method actually is (checked by reading `ladder_fast.py:65-113`).** For each
rung `i`, `W_i(n) = p_{i+2}^3 p_i - p_{i+1}^3 p_{i+3}` is built as an **exact `fmpz_poly`** (no
`float`, no interval, denominators cleared), and the certificate is: find an integer shift `c`
such that `W_i(n = c + m)`, as a polynomial in `m >= 0`, has **all non-negative integer
coefficients** (searched over `pred` and `pred +/- {1,2,3}`, `pred` from the closed-form guess).
This is a legitimate, low-risk, fully exact sufficient condition (a Bernstein/Taylor-shift
positivity certificate) -- structurally much safer than the analytic-estimate machinery in
dense-(a)/(b) (no Stirling remainders, no Edgeworth truncation, no Cauchy tails). So the *method*
is not where I would expect a bug; the *specific 627-rung computation* is what has never been
independently reproduced.

**What I could do without launching the heavy run.** I do not have 7 CPU-hours in this session's
budget, and the instructions forbid launching the heavy full runs. Instead I cross-checked the
**conclusion** the ladder is used for, for a specific `n`, by a *totally different, non-inductive*
method: brute-force exact evaluation of the whole row `M(n, t)` for `t = 1..(n-1)/2` via
`theorem.py`'s own `excess_row` (the same function `theorem.py`'s piece B already uses, but I
evaluated the **full row**, not just `t >= 628`, i.e. including the entire range the ladder is
supposed to cover):

```
n=1259  #t=629  min M = 0.801782  at t=1   > 4/5: True   (row values increasing in t: 0.8018 at t=1 up to > 1.8 near t=n/2)
```

and pushed `M(n,1)` (the tightest point of the whole row -- see finding 3 below) out to
`n = 100001`, confirming it stays strictly above `4/5` and closes the gap like `~2.24/n` from
above at every size tested (no crossing, no sign flip). This is a genuine independent check of the
*numerical conclusion* for specific `n`, but it is **not** a substitute for an independent
re-implementation of the general-`n` polynomial argument the ladder makes for every rung
`i = 0..626` at once.

**Verdict: genuinely UNVERIFIED as a process matter, not FATAL on current evidence.** The method is
exact and low-risk by construction, the log passes its completeness/consistency check, a 12-rung
sample reproduces (with the same script), and independent brute-force spot checks at the specific
`n` values the theorem cites (1259) show no sign of trouble anywhere in the row. But this is the
single largest computational commitment in the whole package (7.06 h against dense-(a)'s 288 s and
dense-(b)'s 41 s) and the only piece that has never had a separately-formulated adversarial
re-derivation the way sparse/dense-(a)/dense-(b) did -- each of which turned out to hide a real
(if harmless) bug when someone looked with different code. I did not find a bug in the ladder; I am
reporting that nobody has looked for one the way they looked in the other three pieces.

---

## 3. `dense_certificate_a.py`: the worst box margin is `3.85e-6`, confirmed exact, and it is
   the tightest number in the entire assembly

**File:** `projects/qg-bootstrap/release/scripts/dense_certificate_a.py`, `certify_box` (line 591),
worst box at `eps in [1.488e-4, 1.736e-4]`, `v in [0.41350, 0.41374]`.

I imported `dense_certificate_a` directly (not duplicated) and called `certify_box` on the exact
box the validator's report names as the binding one:

```
eps: [0.0002 +/- 5.12e-5]  v: [0.414 +/- 5.04e-4]
status: ok
Ng.lower(): 0.80000384576851502060890197753906250000000000000000000000000
```

This matches the validator's own recomputation (`N g >= 0.800003845827`) to the precision both
report. **This is a real, rigorously-derived `arb` lower bound, not an estimate** -- but the margin
above `4/5` is `3.85e-6`, smaller by roughly two orders of magnitude than the correction
(`~2e-7` in `N g`) that the item-6 `R4` bugfix (now present in the live code, confirmed by reading
`kappa2_window_average`, line 407) was needed to close. Before that fix the certificate would have
over-claimed by about 5% of its own margin at exactly this box.

**Verdict: FINE, but flagged.** Nothing here is wrong; I am naming it because it is the single most
fragile numeric fact Theorem A depends on. Nothing in the current pipeline re-derives this specific
number independently on every run: `theorem.py --full` re-executes `dense_certificate_a.py` and
checks `exit 0` / the printed `VERDICT` string, not the numeric value against a second formulation.
Any future change to `KST`, `W`, the `SUPK` table, or the window-average error term must be
re-checked at exactly this box before being trusted.

---

## 4. A second, unflagged instance of "float in an otherwise-exact chain" -- checked, safe direction

**File:** `dense_certificate_a.py:668`:

```python
x_lo = arb(max(628.0, float(N_lo * theta.lower())))
```

converts a rigorous `arb` lower bound to a Python `float` (via Python's built-in `max`) and back
into `arb`. Referee finding M5 flagged the analogous pattern in `theorem.py`'s coverage lemma
(already fixed, confirmed above); this instance in `dense_certificate_a.py` was not named by M5 and
is still present.

**Measured directly:** 200,000 random `arb` lower-bound conversions to `float` never rounded
*upward* past the true value in this session's test (`check_float_line.py`); and specifically for
the box-partition boundaries in `theorem.py`'s coverage lemma (`float(DB.VMAX) * i / 16`), I
computed the exact rational value of the double `0.17` and found it **overshoots** `17/100` by
`11/900719925474099200 ~ 1.2e-17` -- the safe direction (the float-built boxes cover slightly more
than `[0, VMAX]`, not less).

**Verdict: SLACK-ONLY.** No sign flip observed or plausible: any rounding error from `float()` here
is of order `1e-16` relative, roughly 11 orders of magnitude below the `3.85e-6` margin at the
worst box in finding 3, and 5 orders below the smallest per-box slack anywhere in the sweep
(`the sweep`'s items report margins from `3.85e-6` up). Still a genuine violation of the
project's own "no float enters a comparison" rule and of the `.claude/rules/experiment-contract.md`
spirit; worth a one-line fix (`toarb`/exact `fmpq` throughout) at essentially zero risk, same as
M5's fix elsewhere.

---

## 5. Sparse certificate lemma constants -- re-verified independently, all FINE

Targets named in the brief (`FMIN`, tail/moment lemmas `C_i` beyond `i=60`): the validation report
already found and fixed three real rigor gaps (domain of `m`, the `i>30` moment sum, the `j=2`
multiplier), all numerically harmless (subtracted corrections `2.5e-6 -> 3.3e-6` against a margin
of `0.4375`). I independently extended `mu_deviation_bounds` (imported the live function, not
duplicated) from `i=60` to `i=80` and compared against the closed-form bound `C_i <= i^2 * 2^i`
used by `hsum()` for `i > 60`:

```
i=58: C_i=2.05e18  bound=9.70e20  OK   (ratio ~0.0021)
i=80: C_i=1.19e25  bound=7.74e27  OK   (ratio ~0.0015, stable, not degrading toward 1)
```

Holds with a stable ~500x margin, not closing as `i` grows. `Fmin = arb("0.669")` in
`sparse_certificate_full.py` is confirmed by the validator's own 512-box interval sweep to be a
valid, correctly-directioned lower bound (true minimum `0.66920516 > 0.669`). No `CF` constant
exists in `sparse_certificate_full.py` (that name belongs to `sibuya_sparse_certificate.py`,
Theorem B, out of scope here).

**Verdict: FINE.**

---

## Summary table

| # | claim | file:line | verdict | measured margin |
|---|---|---|---|---|
| 1 | corner `ok` check in dense-(b) | `dense_certificate_b.py:599` | SLACK-ONLY (code/claim mismatch, numerically safe) | loss `4.57e-16` vs `S = 0.3352` (15 orders) |
| 2 | ladder (piece A) has no independent re-implementation | `ladder_fast.py`, `lab/ladder_log_merge.py` | UNVERIFIED (process gap, not fatal on current evidence) | brute-force spot check at n=1259..100001: no violation found |
| 3 | dense-(a) worst box | `dense_certificate_a.py:591-675` | FINE (confirmed exact) | `N g - 4/5 = 3.85e-6` |
| 4 | residual float() usage | `dense_certificate_a.py:668` | SLACK-ONLY | rounding `~1e-16` vs margin `3.85e-6` |
| 5 | sparse lemma constants, `C_i` beyond `i=60`, `Fmin` | `sparse_certificate_full.py` | FINE | `C_i` margin ~500x, stable |
| -- | referee M1/M2 (missing logs) | -- | FIXED, confirmed | logs exist, reproduced fresh |
| -- | referee M5 (float in coverage lemma) | `theorem.py` | FIXED, confirmed | final asserts are pure `arb` |

## What I did not check

I did not re-derive the Edgeworth polynomial algebra, the Stirling/Bernoulli remainder machinery,
or the Hadamard-zero `SUPK`/`2.7/3^j` constants from scratch -- the two dense-(a) validation
rounds already did that with a separately formulated engine and I saw no reason to duplicate it. I
did not launch the 288 s dense-(a) full sweep or the 7 h ladder recomputation; I called the live
functions at specific, previously-identified worst points instead, per the brief's "import the
pieces" instruction. I did not check Theorem B (Sibuya) at all, per scope.

## The single most dangerous unverified assumption

**That the 627-rung ladder's construction of `W_i(n)` is correct for every `i = 0..626`, checked
only by the same script that built it (a 12-rung self-sample) plus a format-level consistency
check of its own log -- never by an independently formulated re-derivation, unlike every other
piece of Theorem A.** The method is exact and structurally low-risk (integer polynomial
positivity, no floats, no analytic estimates), and I found no numerical evidence of a problem
(brute-force spot checks at `n = 1259` through `n = 100001` all confirm `M(n,t) > 4/5` including at
`t=1`, the row's tightest point, which is exactly the case rung `i=0` of the ladder is supposed to
cover for every `n >= 5`). But "exact code can still contain an exact bug," and this is the one
piece of Theorem A -- 7.06 of the roughly 7.5 total CPU-hours behind the whole theorem -- that has
never had the kind of adversarial, separately-formulated scrutiny that already found and fixed real
(if harmless) errors in all three other pieces.
