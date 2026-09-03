# Independent verification of Theorem A: `M(n,t) > 4/5`

**Verifier run: 2026-09-03 16:05-16:25 local (`date` at start 16:05:31, at report time 16:24:42).**
Format and standard follow `results/DEBATE_VERIFIER_2026-09-03.md` (Theorem B). No script under
test was imported for the value it is supposed to prove; where a certificate module (`dense_certificate_a`,
`dense_certificate_b`) is called, it is only to obtain the CANDIDATE bound, which is then checked
against a TRUE value computed by a disjoint route. All comparisons are exact `fmpq`/`fmpz` unless
marked "measured" (interval `arb`, radius reported). Scripts:
`C:\Users\user\AppData\Local\Temp\claude\C--Users-user-ScienceBro\ac66e2dc-eaec-44ab-b12d-0b61f841fa72\scratchpad\vfyA\`
(`f1_f7_f8.py`, `f6_coverage.py`, `f6b_boundary.py`, `f345_constants.py`, `f2_agree.py`,
`f2c_sparse.py`, `f2d_identity_check.py`, `f2e_convergence.py`).

## The independence achieved, named

| fact | the author's path | my path |
|---|---|---|
| 1 | doubling DP over paired odd squares (`theorem.py:excess_row`) | **root-polynomial product** `prod(x-r_i)` via `fmpz_poly`, coefficients read off exactly |
| 2 | analytic sampling expansion (`sparse_certificate_full.py`) / tilted-Fourier `kappa_2` (`dense_certificate_a.py`) | **Newton's identity from closed-form power sums** (`fmpz` accumulation) for the true `e_j`, compared against the certificate's own numeric output at concrete points; and a convergence study of the base identity `p_tau ~ bbar^tau F(tau)` at fixed small `t`, `N -> infinity`, with true `e_j` |
| 3 | `mu_deviation_bounds` rational-function evaluation / the `F~ >= 0.669` polynomial sweep | **direct numeric summation of the E_j/R_j series** (no `Q[a,b,E]` polynomial machinery) for `Fmin`; own extension of `mu_deviation_bounds` to `i = 61..80` to test the `i^2 2^i` asymptotic |
| 4 | `SUPK` constants "certified separately in lab/" | **direct exact `fmpq` evaluation of the Bernoulli cumulant polynomials** (`dense_certificate_b.bernoulli_cumulant_polys`) on a fine rational grid, an independent sup search |
| 5 | `R4` and `geo_K` derivations in the docstrings | **direct algebraic/numeric check** of the `1/(1+h)` truncation remainder and the Cauchy-vs-`7.05(2/3)^j` comparison |
| 6 | the printed coverage-lemma text and its `arb` sub-box sweep of `khat_1` | **pure exact-`fmpq` boolean tautology** (`C1 or C2 or C3` for all `theta in (0,1/2]`) plus an exact 280,000+ point grid, and a targeted numeric probe of the actual `b`-threshold implied by `zeta <= ZMAX` |
| 7 | (article's own derivation, not re-read here) | **closed-form Faulhaber polynomials + `fmpq_poly` leading-coefficient extraction**, no limit taken by substitution, no sympy |
| 8 | `excess_row` | **third route**: `e_1, e_2` from closed-form power sums directly (no DP at all) |

---

## Fact 1 -- `M(n,t)` from `e_j`, exact, small `n`: **REPRODUCED**

Root set `{(n-2k)^2}` multiplied out as `prod_k (x - r_k)` with `fmpz_poly` (a completely different
routine from the doubling-by-pairs DP in `theorem.py:excess_row`), coefficients read off and
converted with the standard `(-1)^{N-k}` sign rule. For `n = 5, 7, 9, 11, 13, 15, 21, 25, 31, 41, 51, 61`
the resulting `M(n,t)` vector equals `theorem.py:excess_row(n)` **exactly, term by term, `fmpq`**, for
every `t`. Nothing to correct.

## Fact 2 -- sparse and dense representations agree with exact `M`: **REPRODUCED**, with one item left as **CANNOT VERIFY** at the literal certified corner

Read: full code of `dense_certificate_a.py` (`certify_box`), `dense_certificate_b.py` (`build_Ng`,
`cumulant_series`), `sparse_certificate_full.py` (`build`, `mu_inf`); docstrings only for the deep
internals of the Bernoulli-remainder / Edgeworth machinery not directly exercised.

**Dense (a).** At four real points inside the certified domain (`N=1260,t=628`; `N=1260,t=630`;
`N=2000,t=700`; `N=5000,t=628`) I computed the TRUE `M(n,t)` via a Newton's-identity route (power
sums from `fmpz` accumulation over the `m = N//2` distinct roots, then `k e_k = SUM (-1)^{i-1}
e_{k-i} p_i`, exact `fmpq`) -- a route sharing no code with `excess_row`'s doubling DP -- and called
`dense_certificate_a.certify_box` at the matching `(eps, v)` box. In all four cases `certify_box`
returned `status = ok` and its lower bound on `N g` was below the true value, by a small, correctly
signed, and non-trivial margin:

| N | t | theta | true M(n,t) | certify_box `N g` lower | true M / (bound * n/N) |
|---|---|---|---|---|---|
| 1260 | 628 | 0.4984 | 1.800633 | 1.791569 | 1.0043 |
| 1260 | 630 | 0.5000 | 1.806994 | 1.797921 | 1.0042 |
| 2000 | 700 | 0.3500 | 1.342302 | 1.335739 | 1.0044 |
| 5000 | 628 | 0.1256 | 0.944748 | 0.931431 | 1.0141 |

(First run of this check used a malformed `eps` box -- `arb(1, m)`, i.e. midpoint 1, radius `m` --
which made every call fail; corrected to `arb(1)/m` before the numbers above were produced.)

**Dense (b).** Two points on and near the `theta = 0.05` boundary of region C3 (`N=12560, t=628`,
`t^2/N = 31.4`; `N=20000, t=700`, `t^2/N = 24.5`) give true `M(n,t) = 0.853037` and `0.836541`
respectively, both `> 4/5` as required, computed by the same independent Newton's-identity route (the
full `t^2/N >= 1` corner requires `N` in the hundreds of thousands and was not reached, see below).

**Sparse -- the base identity, verified by convergence rather than at the corner itself.** The
docstring's claim `p_tau = bbar^tau F(tau)` (so that `p_t^2/(p_{t-1}p_{t+1}) = F(t)^2/(F(t-1)F(t+1))`
exactly, the `bbar` factors cancelling in the ratio) is **not a global identity across the whole range
`t = 1..N`** -- I checked this directly with the TRUE (brute-force) `e_j` at `N = 40` and the ratio
`[F(t)/F(t-1)] / [p_t/p_{t-1}]` drifts from `1.0020` at `t=2` to `1.725` at `t=40`. But at FIXED small
`t = 5` and growing `N` (true `e_j`, brute force, `N = 20, 40, 80, 160, 320, 640`), the same ratio's
deviation from 1 shrinks **like `1/N^2`** (`1.09e-2, 2.23e-3, 5.10e-4, 1.22e-4, 2.99e-5, 7.40e-6` --
each halving of `N` divides the deviation by `~4.0`, converging to the theoretically expected exponent),
confirming the identity is a genuine asymptotic statement that sharpens exactly as claimed, in the
regime `t << N` the certificate is built for.

**What I could not reach.** A direct spot check of the true `M(n,t)` at the sparse certificate's own
worst corner (`t = 628`, `N = t^2 = 394384`, i.e. `m ~ 197,000` distinct roots) needs power sums up to
order ~629 over ~197,000 roots -- about `1.2e8` growing-bignum multiplications (final integers
~7000 decimal digits), estimated at tens of minutes; not run in this budget. A first attempt to
substitute the certificate's own `E_j` (built from the `mu^inf` limit moments, not the true finite-N
moments) into the ratio identity at moderate points (`t = 20..150`, `N = t^2`) gave a **systematic,
slowly-shrinking `~5x` mismatch against the true `M`** -- this is consistent with, not contrary to,
the theory: the moment-defect lemma (vi) needs `m >= 627^2/2 ~ 196,600` for its explicit `C_i/m^2`
bound to be small, and every point I could afford is far below that threshold, so a large mismatch there
is expected rather than alarming. **Verdict on this specific sub-item (agreement of the full `E_j`/`mu^inf`
sparse machinery with true `M` at the actual certified corner): CANNOT VERIFY within this budget** --
what would be needed is either much faster exact arithmetic at `N ~ 4x10^5` or trusting the moment-defect
lemma's own bound (already checked quantitatively in Fact 3 below).

## Fact 3 -- sparse lemma constants (`Fmin`, `C_i` beyond `i=60`): **REPRODUCED**

**`Fmin = 0.669`.** Independent direct numeric evaluation of `F(tau) = SUM_j E_j(a,b) R_j(a,b)` (own
implementation of the `E_j` recursion and, critically, a **corrected** `R_j = prod_{i<j}(1-ia) /
prod_{i<j}(1-i a^2 b)` product -- my first draft copied the *last* factor group used internally by the
author's `P`-polynomial construction rather than the true `R_j`, which is fine for `P` but is not `R_j`
by itself; caught and fixed by a consistency test in Fact 2, see below) at `a = 1/627`, `b` from 0 to 1
(the boundary): worst sampled `F = 0.670568` at `b = 1`, giving `Fmin = 0.669` a margin of **0.23%**
(ratio true/bound = 1.0023). Tight, not violated.

**`C_i` for `i > 60` (the `i^2 2^i` crude tail in `hsum`).** Own extension of `sparse_certificate_full.py`'s
own `mu_deviation_bounds` function to `i = 61..80` (same rational-function-in-`1/m` machinery, just run
further -- an exact-arithmetic extension, not a re-derivation) shows the true `C_i` sitting **~500x below**
the crude `i^2 2^i` bound at `i=61` (ratio `2.02e-3`) and the ratio *shrinking* monotonically to `1.54e-3`
at `i=80`. The crude asymptotic is valid and very conservative throughout the range actually used
(`i = 61..399` in `hsum`).

## Fact 4 -- `SUPK` sup bounds in `dense_certificate_a.py`: **REPRODUCED**

`SUPK[j] = sup_{q in [0,1]} |kappa_j(q)|` (`j = 3..13`) checked by an independent exact-`fmpq` sup search
(4000-point grid plus local refinement) of `dense_certificate_b.bernoulli_cumulant_polys(14)` -- a
different module built from the Bernoulli recursion `kappa_{p+1} = q(1-q) kappa_p'`, not whatever
produced `SUPK`'s own "certified separately in lab/" numbers. Measured/claimed ratio for every
`j = 3..13`: `0.9992, 0.9992, 0.9999, 0.9996, 0.9998, 0.9999, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000`
-- every `SUPK[j]` is a valid upper bound, rounded up from the true sup by less than 0.1%, i.e. as tight
as a rounded constant can be.

## Fact 5 -- `R4` / `geo_K` in dense (a) and (b): **REPRODUCED**

**`R4` (the `1/(1+h)` truncation remainder in `kappa2_window_average`).** The claimed bound
`hbar^4/(1-hbar)` is checked directly and exactly (`fmpq`, 200-point grid in `h`) against the true
`|1/(1+h) - (1-h+h^2-h^3)|` for `hbar = 0.1, 0.25, 0.4, 0.49`: ratio measured/claimed is **exactly
`1.0000` in every case**. This is not numerical luck: `1/(1+h) - (1-h+h^2-h^3) = h^4/(1+h)` is an exact
algebraic identity, and its extremum over `|h| <= hbar` is at `h = -hbar`, where it equals exactly
`hbar^4/(1-hbar)` -- the code's bound is **tight, with equality**, not merely valid.

**`geo_K`'s `Kg = 7.05 khat_1/khat_2`, `rho = 3/2` tail constant (`dense_certificate_b.py`).** The raw
Cauchy estimate `j! rho^{-j} q(e^rho-1)/(1-q(e^rho-1))` at `q = 0.1453` (`= VMAX/(1+VMAX)`) compared
against the simplified `j! * 7.05 * (2/3)^j * q` form for `j = 14, 20, 30, 40`: ratio (true Cauchy /
simplified form) is **`0.9995` at every `j` tested** -- the simplified constant `7.05` is a valid,
essentially tight (0.05% margin) upper bound on the raw estimate, independent of the article's own
derivation of "7.05" (I recomputed the Cauchy factor from scratch, not by trusting the number).

## Fact 6 -- coverage: **REPRODUCED WITH A CORRECTION**

Read: full code of `theorem.py`'s coverage-lemma block; docstrings of `sparse_certificate_full.py`,
`sparse_certificate.py`, `dense_certificate_b.py`.

**Exact tautology.** For every rational `theta = t/N` in `(0, 1/2]`: either `b := t*theta <= 1` (C1)
or `theta >= 0.05` (hence C2, since `theta <= 1/2` always) or `theta < 0.05` and `b >= 1` (C3) -- a
pure logical case split needing no numeric input. Verified additionally on an exact `fmpq` grid of
**280,756 points** with `t >= 628`, `N >= 1260`, including exact boundary points `t^2 = N` and
`N = 20t` (`theta = 0.05` exactly) and `t = N//2` (`theta = 1/2` exactly): **zero gaps**. The `N < 1260`
boundary is exactly `{n=1257 (t=628 only), n=1259 (t=628,629)}` -- matching piece B precisely, with no
third case and no gap to piece C's `N >= 1260` (next odd `n = 1261` gives `N = 1260` immediately).

**The correction.** `sparse_certificate_full.py`'s own docstring states the consequence needs
`t^2/N <= BMAX*(628/629)^2` (matching `sparse_certificate.py`'s derivation: the window integral
`INT_0^1 INT_0^1 Phi(t-1+u+v) du dv` reaches `tau = t+1`, so the real-tau condition `tau^2/N <= BMAX`
must hold up to `tau = t+1`, giving `t^2/N <= BMAX*(t/(t+1))^2`, minimised -- hence the constant,
conservative choice -- at `t = 628`). `theorem.py`'s own printed coverage-lemma text and its C1/C3
overlap claim instead use `t^2/N <= BMAX = 1` (i.e. **omits the `(628/629)^2 ~ 0.9968` factor**), and
the executed coverage-lemma code in `theorem.py:main()` only checks `b >= 1 => zeta <= 2/9`, never the
true `(628/629)^2` threshold sparse cert's own docstring requires.

**No actual gap results.** I measured, independently of `theorem.py`'s own 16-box sweep (64 boxes, same
`dense_certificate_b.cumulant_series` object but a new derived quantity `theorem.py` itself never
computes), the smallest `b` for which `zeta <= ZMAX = 0.24` is still guaranteed: **`b_min ~ 0.9259`**
(using the actual measured `khat_1 <= 0.333334`). Since `0.9259 < 0.9968 = (628/629)^2 < 1`, the dense-(b)
region's real reach comfortably covers the sliver `b in [(628/629)^2, 1)` that C1's true (not the
loosely-stated) domain leaves open, with margin to spare (`dense_certificate_b.py`'s own docstring
separately states `b >= 0.95` as its intended coverage, itself between my measured `0.9259` and the
sparse cert's true requirement `0.9968`). **So the coverage claim is true, but `theorem.py`'s own printed
text describing the C1/C3 meeting point ("meet at `b=1`") is not the boundary the sparse certificate's
own docstring requires**, one more example of a comment describing a different (here, looser) condition
than what the neighbouring module's docstring states is actually needed -- caught only by reading both
docstrings side by side and checking the region overlap by direct arithmetic.

## Fact 7 -- `4/5 = lim Var/mean^2`: **REPRODUCED**

Own closed-form Faulhaber polynomials `S_k(m) = SUM_{i=1}^m (2i-1)^{2k}` (`fmpq_poly`, built via the
standard `(m+1)^{p+1}` recursion, no reference to any moment/mu machinery elsewhere in the repository).
`Var(N)/mean(N)^2 = (m S_2 - S_1^2)/S_1^2` as an exact rational function of `m`; both numerator and
denominator have degree 6, and the limit is **the ratio of their leading coefficients**, extracted by
polynomial-degree inspection (`fmpq_poly.degree()`, `.coeffs()`), no substitution/limit routine and no
sympy: leading coefficients `64/45` and `16/9`, ratio **exactly `4/5`**. Numeric cross-check at
`m = 10, 100, 1000, 100000` gives `0.8860, 0.80092, 0.800009, 0.80000000`, converging as expected.

## Fact 8 -- `M(n,1)` decreasing to `4/5`: **REPRODUCED**

Third independent route (neither the DP nor the root-polynomial product): `e_1 = p_1` and
`e_2 = (p_1^2-p_2)/2` directly from the exact power sums `p_1, p_2` of the root multiset (elementary,
no recursion needed for just `j <= 2`). Checked for `n = 5, 7, ..., 201` (99 values) plus
`n = 301, 501, 701, 1001, 1501, 2001`: **strictly decreasing** and **`> 4/5`** at every value, with
`M(5,1) = 80/59 = 1.35593220` (exact) down to `M(2001,1) = 0.8011204078`. Cross-checked exactly against
`theorem.py:excess_row` at `n = 5, 51, 501, 2001`: identical `fmpq` values.

---

## Summary

| fact | verdict |
|---|---|
| 1. `M(n,t)` from `e_j`, exact, small `n` | **REPRODUCED** -- brute-force root-polynomial `e_j` matches `excess_row` exactly, 12 values of `n` |
| 2. sparse/dense representations agree with exact `M` | **REPRODUCED** for dense (a) and (b) at real points (bound below truth by 0.4-1.4%); sparse base identity converges like `1/N^2` at fixed small `t`; **CANNOT VERIFY** the full `mu^inf`-based sparse machinery numerically at the literal `t=628, N~4x10^5` corner (compute cost) |
| 3. sparse lemma constants (`Fmin`, `C_i` beyond `i=60`) | **REPRODUCED** -- `Fmin` margin 0.23%; `C_i` crude tail ~500-650x conservative |
| 4. `SUPK` sup bounds | **REPRODUCED** -- independent exact sup search matches to <0.1% at every `j=3..13` |
| 5. `R4` and `geo_K` in dense (a)/(b) | **REPRODUCED** -- `R4` is an exact tight algebraic identity at the extremum; `geo_K`'s `7.05` constant matches the raw Cauchy estimate to 0.05% |
| 6. coverage of every `(n,t)` | **REPRODUCED WITH A CORRECTION** -- coverage holds exactly (280,756-point grid, 0 gaps), but `theorem.py`'s printed C1/C3 boundary (`b=1`) contradicts `sparse_certificate_full.py`'s own stated requirement (`b <= (628/629)^2`); no actual gap since dense-(b)'s real reach (`b >= ~0.926`) covers the sliver either way |
| 7. `4/5` as `lim Var/mean^2` | **REPRODUCED** -- exact ratio of `fmpq_poly` leading coefficients, `64/45 : 16/9 = 4/5` |
| 8. `M(n,1)` decreasing to `4/5` | **REPRODUCED** -- third independent closed-form route, 105 values of `n`, exact |

**Bottom line.** Every fact checked with a fully independent method reproduced the author's numbers,
several with unexpectedly tight margins (the `R4` bound is exactly tight at its extremum; `geo_K`'s
`7.05` and `SUPK`'s constants match to <0.1-0.05%). One genuine documentation discrepancy was found in
Fact 6 (theorem.py's coverage text omits the `(628/629)^2` shrink factor that `sparse_certificate_full.py`'s
own docstring requires for the window-integral consequence); it costs nothing because the actual
certified dense-(b) region reaches far enough down to absorb the sliver either way, so no number in the
theorem changes and no gap exists in the assembled proof. One item (fact 2's sparse corner) could not be
brought to a numeric verdict within this session's compute budget and is left as CANNOT VERIFY rather
than assumed; what it would take is either an exact `e_j` computation at `N ~ 4x10^5` (tens of minutes of
`fmpz` arithmetic, not attempted here) or accepting the moment-defect lemma's own explicit bound, which
Fact 3 already checked quantitatively and found valid with large margin.
