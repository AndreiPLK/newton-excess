# What the debate broke, and what replaced it

*Two agents attacked the Sibuya package independently on 3 September 2026: one hunting for holes
(`DEBATE_SCEPTIC_2026-09-03.md`), one re-deriving the load-bearing facts by disjoint code paths
(`DEBATE_VERIFIER_2026-09-03.md`).  They converged on the same place from opposite directions.  This file
records the two broken steps, the measurements that broke them, and the replacement.*

## The two steps that were false as written

Both lived in the docstring of `sibuya_corner_grid.tail_ok`, which is proof text: it is the only argument
covering `N` beyond the grid, for every certified index.  The conclusion was never wrong; the derivation was.

**1. A ratio assumed to be Poisson, at `k/H` up to 40.**  The relative defect was bounded by
`|eta_k| <= 0.6 (2j')^2 e^{gamma-H} / H^2`.  Getting there needs `Ehat_{k-2}(H)/Ehat_k(H) <= k^2/H^2` -- true
for `H^k/k!`, and silently assumed for the H-model coefficients.  Measured directly:

    H = 20:  worst  Ehat_{k-2}/Ehat_k * H^2/k^2  =  1.697   at k = 900
    H = 12:                                          2.793   at k = 500
    H = 30:                                          1.344   at k = 900

The assumption is false by up to a factor of 2.8.  The bound it supports survived only on the unstated slack
between `k^2` (with `k <= j'+1`) and the `(2j')^2` that was written -- a factor near 4.  An independent check
from the other side measured the FINAL ratio `|eta| / bound` at 0.31 to 0.42, i.e. the conclusion holds with
about 3x to spare.  A conclusion that holds for a reason its own derivation does not give is not proved.

**2. A substitution in the unsafe direction.**  The same derivation used `1/N <= e^{gamma-H}`.  The truth is the
reverse: `H_N = log N + gamma + d` with `0 < d < 1/(2N)`, so `e^{gamma-H_N} < 1/N`.  Measured ratios
`e^{gamma-H_N} * N`: 0.99950 at `N = 1000`, 0.9999995 at `N = 10^6`.  The same slip made the companion constant
`T(N) <= T_inf (1 + 1.5 e^{gamma-H})` fail on every case tested:

              N       j'     T/T_inf - 1      1.5 e^{gamma-H}
        1000000        1     1.500001e-06     1.499999e-06     fails
        1000000      802     1.500602e-06     1.499999e-06     fails
     1000000000      802     1.500001e-09     1.500000e-09     fails

By a relative `4e-7` to `4e-4` -- far inside the certificate's own slack, and outside its stated derivation.

## The replacement

`tail_ok` no longer estimates a ratio and no longer uses an asymptotic form of `Ehat`.  The defect is carried as
a POLYNOMIAL and the same Taylor shift reads its sign:

    E_k = [x^k] Phi(x) exp(R(x)),   R(x) = SUM_{r>=2} psi^{(r-1)}(N+1) x^r / r! ,
    |psi^{(r-1)}(N+1)| <= (r-2)!/N^{r-1}   =>   |v_r| <= 1/(r(r-1) N^{r-1}) ,
    SUM_{s>=2} |c_s| x^s  <=  exp(g(x)) - 1,    g(x) = SUM_{r>=2} t^{r-1} x^r/(r(r-1)),   t >= 1/N ,

every coefficient of `exp(g)-1` being `t` times a polynomial in `t` with non-negative coefficients.  Evaluating
those at the largest admissible `t` gives scalars `beta_s` with `|c_s| <= t beta_s`, and

    |E_k - Ehat_k|  <=  t * SUM_{s>=2} beta_s Ehat_{k-s}(H) ,

a polynomial in `H` obtained from ONE convolution of `beta` with the H-model coefficients, shared by every
index.  The direction of step 2 is repaired at the source: `1/N = e^{gamma-H_N} e^{d} <= 1.01 e^{gamma-h1}`, and
the companion constant is now 1.6, proved from `T/T_inf - 1 = 1/(2N+3+j') + 1/(N-j') + cross <= 1.51/N` given
`j'/N <= 1e-5` (which holds because `H >= 20` forces `N >= 2.7e8`).

**Cost and effect.**  The new bound is about 2.8x TIGHTER than the one it replaces, and runs in about a second
at `j' = 802`.  Both extreme indices re-certified immediately:

    j' = 399:  True,  shift at h1 = 20, polynomial defect (t1 = 3.71e-09, q = 20)   [0 s]
    j' = 802:  True,  same                                                          [1 s]

## A third finding, of a different kind

The verdict line of the same script PRINTED "exact recursion N <= 10^6 by sibuya_harmonic.py" while that module
was never imported and never ran in the process.  The run existed, elsewhere and earlier -- but the sentence was
an f-string, not evidence, and a referee reading the script would have taken it for one.  The verdict now reads
the log of that run, checks it says `ok = True`, and refuses to print unless the log covers every index the run
asserts.  Related: `exact_log_ok` in the same file, and THE LOG LAW in `CLAUDE.md`.

## Two more repairs, cosmetic but recorded

A comment claimed "measured at k = 800: 6400 bits give a relative radius 2e-6"; the measured value is near
`1e-861`.  Wrong by 850 orders, in the safe direction, and not a measurement.  And the third-order step bracket
omits the slack of `h_3 <= D^3/6` from its stated pieces; the enclosure survives because the tail is taken at
`E(N')` rather than `E(N)`, which absorbs it, so the bracket's true margin is `O(D)` -- 2.9e-2 at the released
step -- and not the ~50% its constant suggests.  Both docstrings now say so.

## What did not break

Every number in the released certificates stands.  The identity `e_j(1..N) = N! e_{N-j}(1,1/2,...,1/N)` was
reproduced exactly by a path sharing no code (dynamic programming over the root set, exact rationals), on 14
cases; the transcription of Sibuya's inequality was reproduced on all 702 pairs with `n <= 40`; the Gamma-ratio
series was reproduced through Hurwitz zeta power sums, sharing no flint routine with the `lgamma` path; and the
ratio form of the cumulants was shown to be an exact identity rather than an approximation for `n >= 1`.
