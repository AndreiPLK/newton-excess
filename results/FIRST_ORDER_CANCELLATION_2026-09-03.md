# The index ceiling was a bookkeeping loss, not a property of the problem

*Found on 3 September 2026 while asking why the H-model certificate stops at one index and not another.
The answer moved the ceiling from about 1300 to at least 3000 at the same starting height, and cost nothing
in rigour.  Measurements below are from `sibuya_corner_grid.py` at 6000 bits.*

## What was being paid for

The finite-`N` defect enters through `E_k = Ehat_k + v_2 Ehat_{k-2} + ...` with `v_2 = psi'(N+1)/2 ~ 1/(2N)`.
Bounding each `E_k` separately pays `|v_2| Ehat_{k-2}` three times over, and

    Ehat_{k-2}/Ehat_k  ~  k^2/H^2 ,

so the bill is a relative `k^2/H^2` -- **1600 at `j' = 802, H = 20`**.  Against a Newton gap of relative size
`~ 1/(k H)`, that is what fixes the largest index a given starting height can carry.

## What is actually there

The first-order defect does not enter `F = E_{j'}^2 - T E_{j'+1} E_{j'-1}` three times; it enters once, in

    B  =  2 Ehat_{j'} Ehat_{j'-2}  -  T_inf ( Ehat_{j'+1} Ehat_{j'-3} + Ehat_{j'-1}^2 ) ,

and the `k^2` cancels inside it.  In the Poisson limit `Ehat_k = H^k/k!` the bracket is **exactly**
`-2 H^{2k-2}/(k!)^2`: expanding the falling factorials, the `k^2` and `k` terms vanish identically and a
constant `-2` is left.  For the true H-model, `B H^2 / Ehat_{j'}^2` measures

        j'        5         20        100        400       1000
    H = 20     -2.034     -0.963     11.01      74.67     235.88
    H = 12     -2.096      0.647     37.21     312.35    1308.10
    H = 40     -2.008     -1.669      2.469     22.19      66.03

against the `k^2` that was being paid: a gain of **12x at `j' = 5` and 4240x at `j' = 1000`** (H = 20), and
15000x at H = 40.  The same cancellation happens at second order: `C = Ehat_{j'-2}^2 - T_inf Ehat_{j'-1}
Ehat_{j'-3}` is positive because its Newton factor `(j'+1)/j'` is the weaker one, so that term helps too.

## Why the sign has to be handled, and how

`B` is positive where it matters and negative in the far tail -- the Poisson limit above says `B -> -2 ...`
as `H -> infinity` at fixed `k`.  Measured signs:

        j'        H=20   H=25   H=30   H=40   H=60   H=100
        40          +      +      -      -      -      -
       200          +      +      +      +      +      +
       802          +      +      +      +      +      +
      1300          +      +      +      +      +      +

The turning point moves out with the index.  So the certificate is written in two pieces meeting at an `h2`
chosen per index: below `h2` the bracket is certified non-negative and the first-order term is DROPPED (it
helps, so dropping it is a valid lower bound); above `h2` the per-term bound suffices, because `t <= e^{gamma-h2}`
has fallen far enough, and one Taylor shift settles `[h2, infinity)` at once.

## The instrument that made it usable

Certifying `B >= 0` on a FINITE interval is not what a Taylor shift does -- a shift at `a` with non-negative
coefficients asserts positivity on all of `[a, infinity)`, which here is false.  The bisection fallback went
exponential at degree 1600.  Replaced by the Bernstein / Moebius certificate: every coefficient of

    (1 + u)^d  p( a + (b - a) u/(1 + u) )

non-negative implies `p >= 0` on `[a, b]`, in one pass of `d` linear steps.  **Rule worth carrying: when a
quantity is positive on a range and negative outside it, the shift is the wrong certificate and the Moebius
substitution is the right one.**

## Effect, measured

    j' = 400    True   [2 s]        j' = 2000   True   [114 s]
    j' = 802    True   [8 s]        j' = 3000   True   [419 s]

The naive version failed at `j' = 1500` at the same starting height.  Cost grows about as `j'^3.5`, so the
remaining limit is arithmetic, not the bound.

## What it does NOT do

It does not close the open region.  The certified window is `33 <= j' <= 802` because that is where the GRID
has run, and the grid, not the tail, is the expensive half.  Extending the window by a factor is grinding: the
open region is `j' >= 803` and stays infinite whatever finite window is closed.  Recorded here so that the
next reader knows the ceiling was never the obstruction.
