# The H-model tail is a classical theorem, not a computation

*Found while trying to make the H-model uniform in the index. Verified numerically the same hour; the machine
certificates it replaces are in `sibuya_corner_grid.py` (770 Taylor shifts, one per index).*

## The observation

The `H`-model of `sibuya_corner_grid.py` works with

    Ehat_k(H) = [x^k] e^{(H - gamma)x} / Gamma(1+x) ,

the `N -> infinity` limit of `E_k(N) = e_k(1, 1/2, ..., 1/N)` at fixed harmonic number `H = H_N`. Using the
Weierstrass product `1/Gamma(1+x) = e^{gamma x} PROD_{n>=1} (1 + x/n) e^{-x/n}`, the generating function is

    Phi(x) = SUM_k Ehat_k(H) x^k = e^{Hx} PROD_{n>=1} (1 + x/n) e^{-x/n} .

That is a **Laguerre-Polya function**: genus one, no Gaussian factor, and its zeros are exactly
`x = -1, -2, -3, ...` -- all real, all negative, all simple. Laguerre-Polya functions are locally uniform limits
of real-rooted polynomials (`e^{Hx} = lim (1 + Hx/m)^m`), so Newton's inequality passes to them by Hurwitz, in
the limit form obtained from `p_k = a_k/C(N,k)` as `N -> infinity`:

    a_k^2  >=  ((k+1)/k) a_{k-1} a_{k+1} ,        a_k = Ehat_k(H) .

**That is exactly the statement `Fhat >= 0` that the 770 Taylor-shift certificates verify** -- `Fhat` is
`Ehat_{j'}^2 - ((j'+1)/j') Ehat_{j'+1} Ehat_{j'-1}`, and `T_infinity = (j'+1)/j'` is precisely the limiting
Newton factor.

## What this explains and what it does not

**Explains:** why every shift certificate passed, at every index and every starting height tried (`j' = 33..802`
at `h1 = 20`; `j' = 1000, 2000, 3000` at `h1 = 20..80`). They were re-deriving a classical inequality, one index
at a time. It also says the shifts will keep passing for every index -- the H-model needs no uniformity proof.

**Does not solve:** the finite-`N` problem. Sibuya's inequality is Newton's inequality TIMES the extra factor
`1 + 1/(2N+3+j')`, and that factor is exactly what vanishes in the limit. So the limit statement is free, while
the finite-`N` statement needs a quantitative lower bound on the Laguerre-Polya gap

    G(k, H) := a_k^2 - ((k+1)/k) a_{k-1} a_{k+1}  >  (the O(1/N) defect) ,

which is what `sibuya_corner_grid.tail_ok` bounds numerically and what the corner of the wedge could not decide
(there the two sides cancel to `+/- 1e-12`).

## Verification

`u_m = [x^m] e^{-gamma x}/Gamma(1+x)` at 4000 bits, `Ehat_k(H) = SUM_m u_m H^{k-m}/(k-m)!`:

    H = 3, 12, 40;  k = 2..59:   a_k^2 - ((k+1)/k) a_{k-1} a_{k+1} > 0 in every case, no violation.

## Why it matters for the write-up

A referee reading 770 machine certificates for a limit statement will ask what they are for. Replacing them with
"the limit function is Laguerre-Polya, so Newton's inequality holds for its coefficients" removes 770 certificates
from the proof and puts a classical citation in their place. The machine work that remains is the finite-`N`
quantitative part, which is where the mathematics actually is.

Related: `results/LITERATURE_VERDICT_STIRLING_NEWTON_GAP_2026-09-02.md` (Sibuya 1988 and what is known),
`release/ARTICLE.md` section 5 (the pieces), `projects/qg-bootstrap/OBSTRUCTION.md` (the open region).
