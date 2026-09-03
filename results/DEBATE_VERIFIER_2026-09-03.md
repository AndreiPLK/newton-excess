# Independent verification of the six load-bearing Sibuya facts

**Verifier run: 2026-09-03 13:41 local.** Written by the VERIFIER role in the debate. No
script under test was imported except in fact 6, where the task is precisely to compare two
implementations; everything else is my own code in `python-flint`, exact `fmpq`/`fmpz` or certified
`arb`/`acb`. No float appears in any comparison a verdict rests on (floats appear only in printed
diagnostics, and where they do not I say so).

## The independence achieved, named

| fact | the author's path | my path |
|---|---|---|
| 1 | generating-function product / Gamma ratio | **dynamic programming over the root set**, exact `fmpq`, root by root |
| 2 | row recurrence for `[n m]` | **unsigned Stirling from its own Pascal recursion** `c(n,k) = c(n-1,k-1) + (n-1)c(n-1,k)`, exact `fmpq`, plus the algebraic identity checked term by term |
| 3 | `arb_series` of `lgamma` | **exact `fmpq` DP** at small `N` and a **second `arb` path through Hurwitz-zeta power sums** `exp(SUM_r (-1)^{r-1} P_r x^r/r)` at large `N` -- no flint routine in common with `lgamma` |
| 4 | the bracket as coded | `h_m` and `E_k` rebuilt from Hurwitz-zeta power sums; the bracket defect computed **exactly rather than bounded**, then compared with the bound |
| 5 | the analytic estimate | the ratio `E_k(N)/Ehat_k(H_N) - 1` **measured** and compared with the claimed bound |
| 6 | Stirling series in two forms | **direct `acb.polygamma` at 1200 bits**, no Stirling series, no Bernoulli numbers |

Scripts: `C:\Users\user\AppData\Local\Temp\claude\C--Users-user-ScienceBro\ac66e2dc-eaec-44ab-b12d-0b61f841fa72\scratchpad\vfy\`
(`f12.py`, `f3.py`, `f45.py`, `f45b.py`, `f4scan.py`, `f6.py`, `f6a.py`).

---

## Fact 1 -- the identity behind everything: **REPRODUCED**

For every one of 14 pairs `(N, j)` with `N` up to 30, `e_j(1..N) = N! e_{N-j}(1,1/2,...,1/N)` holds as an
exact `fmpq` equality **for every `j`**, not only the listed one. The Sibuya ratio transcription

```
p_j^2/(p_{j-1} p_{j+1})  =  [E_{j'}^2/(E_{j'+1} E_{j'-1})] * j'(N-j') / ((j'+1)(N-j'+1)),   j' = N-j
```

is an exact `fmpq` equality on all 14 cases (`p_j = e_j/C(N,j)`). And `1/(3n-j) = 1/(2N+3+j')` with
`N = n-1`, `j' = N-j`, checked exactly for every `n < 60` and every `j`. Nothing to correct.

## Fact 2 -- the target transcription: **REPRODUCED WITH A CORRECTION**

With `[n m]` from my own Pascal recursion, `[n m] = e_{n-m}(1..n-1)` exactly for every `2 <= n <= 40` and
every `m`. On all **702** pairs `(n, j)`, `n = 5..40`, `j = 2..n-2`:

* the boolean `s_m > s_{m+1}` (computed from `s_m = ((m-1)(2n+m)/(n-m+1)) [n m]/[n m-1]` in exact `fmpq`)
  agrees with the boolean `p_j^2/(p_{j-1}p_{j+1}) > 1 + 1/(3n-j)` in **702/702** cases;
* the underlying identity
  `[m(2n+m+1)(n-m+1)/((m-1)(2n+m)(n-m))] * [j(N-j)/((j+1)(N-j+1))] = 1 + 1/(3n-j)`
  holds as an exact `fmpq` equality in **702/702** cases;
* the inequality itself is strict everywhere on that range; worst exact margin `119/828360 = 1.437e-4`
  at `n = 40, j = 2`.

**The correction.** The relation is `>=`, **not** `>`. At `j = 1` -- Sibuya's separately named pair
`m = n-1, n` -- the two sides are **exactly equal** for every `n = 5..40`: `R_1 - (1 + 1/(3n-1)) = 0` as an
`fmpq`, i.e. `s_{n-1} = s_n` identically (both equal 6). Any transcription written as a strict inequality
over `j = 1..n-2` is false at `j = 1`. The repository's own literature note
(`results/LITERATURE_VERDICT_STIRLING_NEWTON_GAP_2026-09-02.md`, lines 85-90) already records this
correctly, and `sibuya_harmonic.py`/`sibuya_corner_grid.py` never touch `j = 1` (they run at `j >= 1001`),
so the certificates are not affected; but the task statement's `>=` is the right sign and `>` would not be.

## Fact 3 -- the Gamma-ratio series and its precision: **REPRODUCED**

`SUM_k E_k x^k = Gamma(N+1+x)/(Gamma(N+1)Gamma(1+x))`:

* every exact `fmpq` coefficient from the root DP lies inside the `arb` ball of the Gamma-ratio series,
  for `(N, kmax) = (12,12), (50,40), (200,40), (1000,50)`;
* at `N = 10^6` and `N = 485165195` the series satisfies `E_k(N) = E_k(N-1) + E_{k-1}(N-1)/N` for every
  `1 <= k <= 60`;
* a second, structurally unrelated `arb` path (Hurwitz-zeta power sums into `arb_series.exp`) overlaps the
  Gamma-ratio series for every `k <= 120` at `N = 1000, 10^6, 485165195`.

**Precision.** First, what actually runs: with `JPMAX = 399` the script calls `series_E(N, 401)` and
`u_coeffs(402)`, so `prec_for(kmax) = max(4000, 8*kmax) = 4000` bits -- **the `8*kmax` branch never fires in
the released run.** Measured worst *relative* radius over all coefficients, against a rerun at ten times
the bits (which agrees everywhere):

| N | kmax | prec | worst relative radius |
|---|---|---|---|
| 10^6 | 401 | 4000 | `1e-718` |
| 485165195 | 401 | 4000 | `1e-778` |
| 10^6 | 800 | 6400 | `1e-715` |
| 485165195 | 800 | 6400 | `1e-861` |
| 485165195 | 800 | 4000 | `1e-138` |

So the precision is ample by hundreds of orders of magnitude, and the claim "8*kmax bits is enough at
k = 800" is true with enormous room. **The comment's number is not reproduced**: `prec_for`'s
"measured at k = 800: 6400 bits give a relative radius 2e-6" is off by ~850 orders of magnitude. It errs in
the conservative direction (reality is far better than claimed), so nothing rests on it, but it should be
corrected or deleted rather than left as a "measurement".

## Fact 4 -- the third-order step bracket: **REPRODUCED WITH A CORRECTION**

The four stated identities all hold (`N = 10^6 -> 1003004`, `arb` overlap, 4000 bits):

* `h_1 = H_{N'} - H_N = D` -- ok;
* `h_2 = (D^2 - S_2)/2`, `S_2 = SUM_{i=N+1}^{N'} 1/i^2` -- ok;
* `h_3 = (D^3 - 3 D S_2 + 2 S_3)/6` (Newton), and `h_m <= D^m/m!` decided for every `2 <= m <= 59`
  (at `m = 1` it is an equality, so the ball comparison is undecidable -- not a violation),
  hence `h_3 <= D^3/6` -- ok;
* `E_k(N') = SUM_m h_m E_{k-m}(N)` at `k = 33, 100, 200, 399, 401` -- ok.

**The `1.5` is justified.** `SUM_{m>=4} h_m E_{k-m}(N) <= 1.5 (D^4/24) E_{k-4}(N')` holds at every point
tested -- `N = 10^6, 10^7, 10^8`, `D = 0.003 .. 0.26`, `j' = 33 .. 900` (at 12000 bits) -- with the largest
ratio to the bound **0.664** (at `j' = 33`), falling to 0.53 at `j' = 900`. The constant actually needed is
therefore ~1.0 and 1.5 gives about 50% headroom. It does not degrade with `j'`.

**The correction, and it is load-bearing.** The docstring's list of bounded pieces is incomplete. Writing
out the bracket,

```
E_k(N') - [E_k + D E_{k-1} + (D^2/2) E_{k-2} + (D^3/6) E_{k-3}]
        = -(S_2/2) E_{k-2}  +  ((3 D S_2 - 2 S_3)/6) * (-E_{k-3})  +  SUM_{m>=4} h_m E_{k-m} ,
```

the middle term -- the `h_3` slack -- is **not** among the pieces `rest` claims to bound, and it is of the
**same order** as the `S_2` piece that is (their ratio is `D * E_{k-3}/E_{k-2}`, about 0.006 at the released
step size). The code is nonetheless a valid enclosure, but for a reason the comment does not state: `rest`
evaluates the `S_2` term at `Eb = E(N')`, not at `E(N)`, and `E_{k-2}(N') - E_{k-2}(N) ~ D E_{k-3}(N)`
absorbs exactly the missing piece. Substituting `E(N)` for `E(N')` there breaks it, measured at
`N = 10^6`, `j' = 33`, `k = 32`:

| step `N'-N` | `D` | `|defect|/rest` (code, with `E(N')`) | with `E(N)` instead |
|---|---|---|---|
| 300000 | 2.62e-01 | 0.426066 | 0.750304 |
| 3000 | 3.00e-03 | 0.970802 | 0.977642 |
| 300 | 3.00e-04 | **0.999969** | **1.000673 -- FAILS** |
| 30 | 3.00e-05 | 0.999998 | 1.000069 -- FAILS |
| 1 | 1.00e-06 | 0.999998 | 1.000001 -- FAILS |

So the true margin of the bracket is not the ~50% the `1.5` suggests; it is `O(D)` -- 2.9e-2 at the step the
released run actually uses (1904 steps over `H` from 14.393 to 20.577 gives an average `D = 3.25e-3`,
`results/sibuya_corner_grid_2026-09-03.txt`) and 2e-6 at `D = 1e-6`. It never fails, but it is
asymptotically tight and it survives by an unnamed mechanism. Recommended: state the `h_3` slack explicitly
among the bounded pieces and note that `Eb` (not `E`) is required in the `S_2` term.

**Two smaller notes.** (i) In `grid_ok`, `S2` and `Dm_eff` are carried as Python floats
(`S2 = float(acb(N+1).polygamma(...) - ...)`, `Dm_eff = arb(float(harmonic(N2) - H))`) -- an exact quantity
represented by a float inside a chain a verdict rests on. The induced relative error is ~1e-16, seven orders
below the run's own worst relative margin of 2.083e-9, so it is harmless here; it is not exact.
(ii) **A false alarm I raised and then killed, recorded so nobody re-derives it:** my first pass at 4000 bits
reported the bracket failing at `j' = 800` (ratios 4.95 and 5.3e6). That was my own cancellation, not the
code's: at 12000 bits the same points give 0.543 and 0.484. The instrument, not the object.

## Fact 5 -- the H-model: **REPRODUCED WITH A CORRECTION**

`|eta| = |E_k(N)/Ehat_k(H_N) - 1| <= 0.6 (2j')^2 e^{gamma-H}/H^2` holds at every point tested, and `Ehat_k`
is positive throughout. In the H-model's own domain (`H >= H1 = 20`, i.e. `N >= 4.85e8`), which is where
`tail_ok` uses it:

| N | H | `j'` | worst `|eta| / bound` |
|---|---|---|---|
| 485165195 | 20.577 | 399 | 0.3129 |
| 485165195 | 20.577 | 800 | 0.3415 |
| 5e9 | 22.910 | 800 | 0.3179 |
| 1e11 | 25.906 | 800 | 0.2971 |

The bound is a factor ~3 conservative and the ratio is flat in `H`, so it does not decay at large `H`.
Outside the domain (`N = 10^6`, `H = 14.393`) the ratio rises to 0.4218 -- `|eta| = 7.780e-4` against a bound
of `1.844e-3` -- still valid. Not wildly conservative, not tight: a factor of about 2.4 to 3.

**The correction.** The justifying step written in `tail_ok`'s docstring, "`1/N <= e^{gamma-H}`", is
**backwards**. Because `H_N > log N + gamma` for every `N`, one has `e^{gamma-H_N} < 1/N`. Measured:

```
N = 1000        1/N = 1.000000e-03   e^(gamma-H_N) = 9.995002e-04   ratio 0.99950021
N = 10^6        1/N = 1.000000e-06   e^(gamma-H_N) = 9.999995e-07   ratio 0.99999950
N = 485165195   1/N = 2.061154e-09   e^(gamma-H_N) = 2.061154e-09   ratio 1.00000000 (to 8 places)
```

The substitution is therefore in the unsafe direction. It costs a factor `1/(1 - 1/(2N)) <= 1 + 1.1e-9` in
the regime where the H-model is applied, which the measured factor-3 slack absorbs with ease, so the
conclusion stands and no number in the certificate changes. The derivation as written does not stand and
should be repaired (replace `1/N <= e^{gamma-H}` by `e^{gamma-H} <= 1/N <= e^{gamma-H}/(1 - 1/(2N))`, or
simply bound with `1/N` directly).

## Fact 6 -- the wedge's ratio form: **REPRODUCED**

Analytically `w^n Phi_n(w) = w^n N^n psi^{(n)}(N w) = z^n psi^{(n)}(z)` with `z = w/eps = 1/u`, so for
`n >= 1` the dependence on `u` alone is an **exact identity, not an approximation**. Verified against a
direct `acb.polygamma` reference at 1200 bits: four pairs `(eps, w)` with `u = 1/1500` and `w` spanning
`1.5e-3 .. 150` (a factor `10^5`) give identical values to 12 digits for `n = 1..5`, and
`sibuya_top_w.phi_ratio` contains the true value at every pair and returns the identical interval at all
four.

**Except at `n = 0`, where the docstring over-reaches.** `Phi_0(w) = psi(1/u) + log u + log w` carries a
separate `log w`. Measured: the differences between same-`u` points are exactly `log(w_2/w_1)` --
`3.40119738166` = `log 30`, `-6.90775527898` = `log 1e-3`, `4.60517018599` = `log 100`, all agreeing to 12
digits. **The code is right** -- its `n == 0` branch keeps `w.log()` explicitly -- only the sentence
"`w^n Phi_n(w) = B_n(u)`, the Stirling series with the `w^{-n}` factored out" is inaccurate at `n = 0`.

Against `sibuya_dense_a.phi_stirling_real`: at thin arguments the two agree and both contain the direct
polygamma value in 15/15 tested `(n, eps, w)` (`n = 0..5`, `eps` from `4.5e-4` to `1e-7`), with identical
widths from `2.0e-50` to `6.4e-15`. At 5%-wide `w`-bands the direct form's enclosure is already `+/-inf` or
`nan` at `n = 1`, while the ratio form has relative width `8.4e-7 .. 4.2e-5` and encloses all 9 sampled true
values across the band -- reproducing, and in fact exceeding, the docstring's claim that the direct form
"returns nan from k = 3 on".

---

## Summary

| fact | verdict |
|---|---|
| 1 identity `e_j = N! E_{N-j}` and the ratio transcription | **REPRODUCED** |
| 2 Sibuya (3.4) `<=>` `p_j^2/(p_{j-1}p_{j+1}) >= 1 + 1/(3n-j)` | **REPRODUCED WITH A CORRECTION** -- `>=`, with equality exactly at `j = 1`, never `>` there |
| 3 Gamma-ratio series and `prec_for` | **REPRODUCED** -- the precision claim holds with vast room; the comment's "2e-6 at k=800" is wrong by ~850 orders (conservatively) |
| 4 third-order step bracket, the `1.5` | **REPRODUCED WITH A CORRECTION** -- valid, `1.5` justified (needed ~1.0), but the `h_3` slack is missing from the stated bracket and the enclosure survives only because `E(N')` is used in the `S_2` term; true margin is `O(D)`, 0.99997 of the bound at `D = 3e-4` |
| 5 H-model relative error bound | **REPRODUCED WITH A CORRECTION** -- bound true, factor ~3 loose; its stated justification `1/N <= e^{gamma-H}` is backwards, harmlessly |
| 6 `w^n Phi_n(w)` depends only on `eps/w` | **REPRODUCED** -- exactly, for `n >= 1`; false at `n = 0`, where the code is nevertheless correct |

Nothing found here invalidates a number in the released certificates. Four docstrings state something other
than what the code relies on, and two of those (fact 4's bracket, fact 5's `1/N` substitution) are proof
text, not commentary.
