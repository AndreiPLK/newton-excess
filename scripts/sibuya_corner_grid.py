"""SIBUYA'S CONJECTURE (1988, eq. 3.4), the window 33 <= j' <= JPMAX of missing indices, EVERY N.

`j' = N - j`.  The statement is `F := E_{j'}^2 - T E_{j'+1} E_{j'-1} > 0` with `E_k = e_k(1, 1/2, ..., 1/N)` and
`T = (1 + 1/(2N+3+j')) (j'+1)(N-j'+1)/(j'(N-j'))`  (see sibuya_harmonic.py for the derivation).

Three ranges, one instrument each:

  1. `N <= N_EXACT`: the exact incremental recursion `E_k(N) = E_k(N-1) + E_{k-1}(N-1)/N` (sibuya_harmonic.exact_pass),
     positive terms only, no cancellation.
  2. `N_EXACT < N <= e^{H1}`: a grid in `N` with EXACT values at the grid points from the Gamma-ratio series
     `SUM_k E_k x^k = Gamma(N+1+x)/(Gamma(N+1) Gamma(1+x))` (arb power series: at a single integer `N` the
     cancellation is a precision issue, and precision is free), and a SECOND-ORDER bracket between them:
         E_k(N') = E_k(N) + D E_{k-1}(N) + theta,  0 <= theta <= D^2 E_{k-2}(N_b),  D = H_{N'} - H_N in [0, D_max],
     so that `F` is a quadratic in `D` with exactly-computed coefficients (the first-order terms of the numerator
     and the denominator nearly cancel, which is why the crude monotone bracket is 40x too lossy).
  3. `N > e^{H1}`: `E_k(N) = Ehat_k(H) + O(e^{-H})` with `Ehat_k(H) = [x^k] e^{(H-gamma)x}/Gamma(1+x)` -- a
     polynomial in `H` with EXACT coefficients (no interval power sums).  `Fhat(H)` is shifted to `H = H1 + s`
     and certified by the LADDER TRICK: all coefficients of the shifted polynomial non-negative => positive for
     every `s >= 0`.  The `O(e^{-H})` difference is bounded explicitly and is below `1e-8` at `H1 = 20`.

Run:  uv run python projects/qg-bootstrap/release/scripts/sibuya_corner_grid.py [JPMAX] [--tail-only]
"""

from __future__ import annotations

import os
import re
import sys
import time
from math import factorial

from flint import acb, arb, arb_poly, arb_series, ctx, fmpq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

JPMIN = 33
JPMAX = int(sys.argv[1]) if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else 802
N_EXACT = 10**6  # closed by sibuya_harmonic.exact_pass
H1 = arb(20)  # above this the Gamma-ratio and its H-model differ by less than e^{-H} < 2.1e-9
PREC = 4000  # raised with the index: E_k ~ H^k/k! is around 10^-1150 at k = 800, so the Gamma-ratio series


# needs enough digits to carry both the value and the cancellation inside exp()
def prec_for(kmax):
    # measured (independent check, 2026-09-03): at kmax = 401 the 4000-bit floor already gives a relative
    # radius near 1e-778 on E_k, and 6400 bits at k = 800 give near 1e-861 -- the 8*kmax branch is a large
    # safety margin, not a necessity.  (An earlier comment here claimed "2e-6 at k = 800", which was wrong by
    # some 850 orders of magnitude; it was conservative, so nothing rested on it, but it was not a measurement.)
    return max(4000, 8 * kmax)


EULER = None


def series_E(N, kmax, prec=None):
    """All E_k(N), k <= kmax, at a single integer N: the Taylor coefficients of the Gamma ratio."""
    ctx.prec = prec or prec_for(kmax)
    ctx.cap = kmax + 2
    x = arb_series([0, 1])
    s = (x + (N + 1)).lgamma() - arb(N + 1).lgamma() - (x + 1).lgamma()
    return s.exp().coeffs()


def u_coeffs(mmax, prec=None):
    """u_m = [x^m] e^{-gamma x}/Gamma(1+x): the exact coefficients of the H-model."""
    ctx.prec = prec or prec_for(mmax)
    ctx.cap = mmax + 2
    x = arb_series([0, 1])
    return (-arb.const_euler() * x - (x + 1).lgamma()).exp().coeffs()


def harmonic(N):
    return arb(N + 1).digamma() + arb.const_euler()


def T_of(N, jp):
    N = arb(N)
    return (1 + 1 / (2 * N + 3 + jp)) * arb(fmpq(jp + 1, jp)) * (N - jp + 1) / (N - jp)


def pmulD(a, b, deg=7):
    out = [arb(0)] * deg
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j < deg:
                out[i + j] += x * y
    return out


def step_margin(E, Eb, N, Dm, jp, S2=0.0):
    """Lower bound of F over one step, THIRD order in D.

    E_k(N') = SUM_m h_m E_{k-m}(N), h_m = e_m(1/(N+1), ..., 1/N') >= 0 the elementary symmetric functions of the
    new reciprocals: h_1 = D = H_{N'} - H_N, h_2 = (D^2 - S_2)/2 with S_2 = SUM 1/i^2 <= S2, h_3 <= D^3/6, and
    SUM_{m>=4} h_m E_{k-m} <= 1.5 (D^4/24) E_{k-4}(N').  D stays symbolic in all three brackets, so the
    first-order cancellation between the numerator and the denominator is kept exactly; the pieces that are not
    polynomial in D (the S_2 defect and the m >= 4 tail) enter as small intervals.  The third order buys a step
    ~15x longer than the quadratic bracket at the same margin.

    TWO THINGS THIS BRACKET DOES AND THE LIST ABOVE DOES NOT SAY (independent check, 2026-09-03).
    (i) `h_3 <= D^3/6` is an inequality, and its slack `((3 D S_2 - 2 S_3)/6) E_{k-3}` is the SAME order as the
    `S_2` piece that is listed.  The enclosure is still valid, but only because `rest` is taken at `Eb = E(N')`
    and the difference `E_{k-2}(N') - E_{k-2}(N) ~ D E_{k-3}(N)` absorbs it: substituting `E(N)` there makes the
    check fail (measured |defect|/rest = 1.000673 at D = 3e-4).  So the true margin of the bracket is O(D)
    -- 2.9e-2 at the released step size D ~ 3.25e-3 -- not the ~50% that the constant 1.5 suggests.
    (ii) `S2` and `Dm` arrive as Python floats, so this bracket carries a ~1e-16 relative slop.  That is seven
    orders below the run's own worst margin (2.083e-9) and therefore harmless, but it is not exact arithmetic.
    """
    T = T_of(N, jp)
    D4 = arb(0, float(Dm.upper() ** 4 / 24 * 1.5))

    def br(k):
        pol = [E[k], E[k - 1], E[k - 2] / 2, E[k - 3] / 6]
        rest = arb(0, float(S2 / 2 * Eb[k - 2])) + D4 * Eb[max(k - 4, 0)]
        return pol, rest

    p0, r0 = br(jp)
    pp, rp = br(jp + 1)
    pm, rm = br(jp - 1)
    sq = pmulD(p0, p0)
    cr = pmulD(pp, pm)
    poly = [sq[i] - T * cr[i] for i in range(len(sq))]
    Db = arb(float(Dm.upper()) / 2, float(Dm.upper()) / 2)
    val = arb(0)
    for c in reversed(poly):
        val = val * Db + c
    slack = (
        2 * abs(E[jp]) * r0 + r0 * r0 + T * (abs(E[jp + 1]) * rm + abs(E[jp - 1]) * rp + rp * rm)
    )
    return (val - abs(slack)).lower(), poly[0]


# ------------------------------------------------------------------ range 2: the grid
def grid_ok(N_lo, N_hi, jpmax, verbose=True):
    """Every integer N in [N_lo, N_hi] and every j' in [JPMIN, jpmax]."""
    t0 = time.time()
    N = N_lo
    steps = 0
    worst = None
    Dm_guess = arb("0.02")
    while N < N_hi:
        E = series_E(N, jpmax + 2)
        H = harmonic(N)
        Dm = Dm_guess
        N2 = min(int(float(arb(N) * Dm.exp())), N_hi)
        if N2 <= N:
            N2 = N + 1
        Eb = series_E(
            N2, jpmax + 2
        )  # an upper bound for every smaller step as well (E_k increases in N)
        for _ in range(12):
            Dm_eff = arb(float(harmonic(N2) - H))
            S2 = float(acb(N + 1).polygamma(acb(1)).real - acb(N2 + 1).polygamma(acb(1)).real)
            ok_here = True
            for jp in range(JPMIN, jpmax + 1):
                lo, c0 = step_margin(E, Eb, N, Dm_eff, jp, S2)
                if not (lo > 0):
                    ok_here = False
                    break
                r = float(lo / (E[jp] ** 2))
                if worst is None or r < worst[0]:
                    worst = (r, N, jp)
            if ok_here:
                break
            N2 = N + max(1, (N2 - N) // 3)  # shrink the step; Eb stays a valid upper bound
        else:
            return False, f"grid FAILED at N = {N} (step could not be shrunk enough)"
        steps += 1
        Dm_guess = arb(float(harmonic(N2) - H)) * 2
        if verbose and steps % 500 == 0:
            print(
                f"    grid: N = {N:.4e}, H = {float(H):.3f}, {steps} steps, worst rel {worst[0]:.2e}"
                f"  [{time.time() - t0:.0f} s]",
                flush=True,
            )
        N = N2
    return (
        True,
        f"{steps} grid steps {N_lo:.3e} -> {N_hi:.3e}, worst relative margin {worst[0]:.3e}  [{time.time() - t0:.0f} s]",
    )


# ------------------------------------------------------------------ range 3: the H-model with the ladder trick
def poly_mul(a, b):
    """arb_poly multiplication (C level: the Python double loop was the bottleneck of the tail check)."""
    return list((arb_poly(list(a)) * arb_poly(list(b))).coeffs())


def poly_add(a, b):
    """sum of two coefficient lists of different lengths."""
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else arb(0)) + (b[i] if i < len(b) else arb(0)) for i in range(n)]


def positive_on_interval(c, a, b):
    """p > 0 on [a, b], by ONE pass instead of bisection.

    `p >= 0` on `[a, b]` follows if every coefficient of `(1+u)^d p(a + (b-a) u/(1+u))` is non-negative, since
    `u >= 0` covers `H in [a, b)` and `(1+u)^d > 0`.  This is the Bernstein / Moebius positivity certificate,
    and unlike the Taylor shift at `a` it does NOT assert anything beyond `b` -- which is what is needed for a
    quantity that is positive on a finite range and turns negative outside it.  Cost is one Taylor shift plus
    `d` linear steps; the bisection it replaces went exponential at degree 1600.
    """
    sh = shift_poly(c, a)
    d = len(sh) - 1
    ba = b - a
    acc = [sh[d]]
    ypow = [arb(1)]
    for i in range(d - 1, -1, -1):
        ypow = poly_mul(ypow, [arb(1), arb(1)])
        acc = [arb(0)] + [ba * x for x in acc]
        acc = poly_add(acc, [sh[i] * y for y in ypow])
    return all(x.lower() >= 0 for x in acc)


def shift_poly(c, h):
    """coefficients of p(h + s) in s, by composition with (s + h) in arb_poly."""
    return list(arb_poly(list(c))(arb_poly([h, 1])).coeffs())


def peval(c, H):
    tot = arb(0)
    for x in reversed(c):
        tot = tot * H + x
    return tot


def positive_on(c, lo, hi, depth=0):
    """p > 0 on [lo, hi].  The polynomial is TAYLOR-SHIFTED to the left endpoint before evaluating: a degree-800
    polynomial with alternating coefficients evaluated directly at H ~ 10^3 loses every digit, while the shifted
    form has mostly non-negative coefficients and, when all of them are non-negative, settles [lo, infinity) at once."""
    sh = shift_poly(c, lo)
    if all(x.lower() >= 0 for x in sh):
        return True
    w = arb(0, float(hi - lo))
    if peval(sh, w).lower() > 0:
        return True
    if depth > 60:
        return False
    mid = (lo + hi) / 2
    return positive_on(c, lo, mid, depth + 1) and positive_on(c, mid, hi, depth + 1)


_BETA_CACHE = {}


def beta_conv(kmax, u, h1, smin=2):
    """`C[n] = SUM_{s=2..n} beta_s u_{n-s}`, computed once: the whole defect of the H-model, convolved.

    `E_k(N) = [x^k] Phi(x) exp(R(x))` with `R(x) = SUM_{r>=2} psi^{(r-1)}(N+1) x^r / r!`, so
    `E_k - Ehat_k = SUM_{s>=2} c_s Ehat_{k-s}` with `exp(R) = SUM c_s x^s`.  From
    `|psi^{(r-1)}(N+1)| <= (r-2)!/N^{r-1}` one gets `|v_r| <= 1/(r(r-1)N^{r-1})`, hence
        SUM_{s>=2} |c_s| x^s  <=  exp(g(x)) - 1,      g(x) = SUM_{r>=2} t^{r-1} x^r / (r(r-1)),
    where `t >= 1/N`.  Every coefficient of `exp(g)-1` is `t` times a polynomial in `t` with non-negative
    coefficients, so dividing by `t` and evaluating at the largest admissible `t` gives scalars `beta_s` with
    `|c_s| <= t beta_s` for every `N` in range.  NOTHING here is an asymptotic estimate and no ratio
    `Ehat_{k-s}/Ehat_k` is used: the defect stays a POLYNOMIAL in `H` and the shift trick reads its sign.

    (The earlier version bounded the relative defect by `0.6 (2j')^2 e^{gamma-H}/H^2`, which silently assumed
    `Ehat_{k-2}/Ehat_k <= k^2/H^2` -- the Poisson estimate, applied at `k/H` up to 40.  Measured, that ratio
    reaches 1.70 at `H = 20, k = 900` and 2.79 at `H = 12, k = 500`: the assumption is false, and the bound
    survived only on the slack between `k^2` and `(2j')^2`.)
    """
    key = (kmax, float(h1), smin)
    if key in _BETA_CACHE:
        return _BETA_CACHE[key]
    ctx.prec = prec_for(kmax)
    ctx.cap = kmax + 2
    # t_1 >= 1/N for every N with H_N >= h1:  H_N = log N + gamma + d, 0 < d < 1/(2N), so
    # 1/N = e^{gamma-H_N} e^{d} <= 1.01 e^{gamma-h1} once N >= 100.  (The reverse substitution
    # "1/N <= e^{gamma-H}" used before is FALSE -- e^{gamma-H_N} < 1/N -- and this is its repair.)
    t1 = arb("1.01") * (arb.const_euler() - h1).exp()
    g = arb_series([arb(0)] * 2 + [t1 ** (r - 1) / arb(r * (r - 1)) for r in range(2, kmax + 2)])
    beta = [c / t1 ** (smin - 1) for c in list((g.exp()).coeffs())]
    beta = (beta + [arb(0)] * (kmax + 2))[: kmax + 2]
    for i in range(min(smin, len(beta))):
        beta[i] = arb(0)
    C = poly_mul(beta[: kmax + 1], list(u[: kmax + 1]))[: kmax + 1]
    _BETA_CACHE.clear()
    _BETA_CACHE[key] = (C, t1)
    return C, t1


def tail_bracket(jp, u, h1, h2):
    """F(N) > 0 for every N with H_N >= h1, by ONE polynomial inequality in H -- WITH the first-order
    cancellation kept.

    `E_k = Ehat_k + v_2 Ehat_{k-2} + r_k`, `v_2 = psi'(N+1)/2 in (0, t/2]`, `|r_k| <= t^2 R_k(H)` with
    `R_k = SUM_{s>=3} gamma_s Ehat_{k-s}` (every `c_s` with `s >= 3` carries `t^2` or more).  Substituting,

        F  =  Fhat  +  v_2 B  +  v_2^2 C  -  tau T_inf E_{k+1} E_{k-1}  +  (terms in r) ,
        B  =  2 Ehat_j' Ehat_{j'-2} - T_inf ( Ehat_{j'+1} Ehat_{j'-3} + Ehat_{j'-1}^2 ) ,
        C  =  Ehat_{j'-2}^2 - T_inf Ehat_{j'-1} Ehat_{j'-3} .

    **`B` is where the index ceiling came from.**  Bounding each `E_k` separately pays `|v_2| Ehat_{k-2}` three
    times over, i.e. a relative `k^2/H^2` -- 1600 at `j' = 802, H = 20`.  In the combination the `k^2` cancels:
    in the Poisson limit `Ehat_k = H^k/k!` the bracket is EXACTLY `-2 H^{2k-2}/(k!)^2`, and for the true model
    `B H^2/Ehat_{j'}^2` measures -2.03, -0.96, 11.0, 74.7, 235.9 at `j' = 5, 20, 100, 400, 1000` (H = 20) --
    a gain of 12x to 15000x over the naive bound.  Better still, `B >= 0` wherever it matters, and then the
    first-order term HELPS: it is dropped, not bounded.  `C > Fhat_{j'-2} > 0` likewise (its Newton factor
    `(j'+1)/j' < (j'-1)/(j'-2)` is the weaker one), so it is dropped too.

    What is left to pay for is `tau <= 1.6 t` on `T_inf E_{j'+1} E_{j'-1}` -- relative size 1.6, not 1600 -- and
    the `r` terms, which carry `t^2`.  Multiplying by `H^q`, `q = [h1] <= h1`, and using `e^{-H} H^q <=
    e^{-h1} h1^q` on `H >= h1 >= q` turns the whole statement into one polynomial inequality settled by the
    Taylor shift at `h1`.  If the shift of `B` shows a negative coefficient the routine reports failure rather
    than falling back on the naive bound: at that point the index is genuinely outside the model's range.
    """
    lo = jp - 3
    Ehat = {}
    for k in range(lo, jp + 2):
        c = [arb(0)] * (k + 1)
        for m in range(min(k, len(u) - 1) + 1):
            c[k - m] = u[m] / arb(factorial(k - m))
        Ehat[k] = c
    C3, t1 = beta_conv(jp + 1, u, h1, smin=3)
    R = {k: [C3[k - i] / arb(factorial(i)) for i in range(k - 2)] for k in (jp - 1, jp, jp + 1)}
    Tinf = arb(fmpq(jp + 1, jp))
    G = poly_mul(Ehat[jp], Ehat[jp])
    P = poly_mul(Ehat[jp + 1], Ehat[jp - 1])
    F = [G[i] - (Tinf * P[i] if i < len(P) else arb(0)) for i in range(len(G))]
    # the two top coefficients vanish IDENTICALLY -- that is the statement that T_inf = (j'+1)/j' is exactly
    # Newton's factor for the Laguerre-Polya function e^{(H-gamma)x}/Gamma(1+x); numerically they are balls
    # around 0, and a ball times H^{2j'} would be unbounded, so they are checked and then set to exact zero
    scale = max(float(abs(c).upper()) for c in F)
    for i in (2 * jp, 2 * jp - 1):
        if not (float(abs(F[i]).upper()) <= 1e-30 * scale):
            return False, f"top coefficient at degree {i} did not vanish: {F[i].str(4)}"
        F[i] = arb(0)
    B = poly_add(
        [2 * c for c in poly_mul(Ehat[jp], Ehat[jp - 2])],
        [
            -Tinf * c
            for c in poly_add(poly_mul(Ehat[jp + 1], Ehat[jp - 3]), poly_mul(Ehat[jp - 1], Ehat[jp - 1]))
        ],
    )
    if not positive_on_interval(B, h1, h2):
        return False, f"the bracket B is not >= 0 on [{float(h1):.0f}, {float(h2):.0f}] at j' = {jp}"
    half = t1 / 2
    up = {k: poly_add(Ehat[k], [half * c for c in Ehat[k - 2]]) for k in (jp - 1, jp, jp + 1)}
    up = {k: poly_add(v, [t1 * t1 * c for c in R[k]]) for k, v in up.items()}
    # order t : the T defect on E_{j'+1} E_{j'-1}
    W1 = [arb("1.6") * Tinf * c for c in poly_mul(up[jp + 1], up[jp - 1])]
    # order t^2 : everything carrying r
    W2 = poly_add(
        [2 * c for c in poly_mul(up[jp], R[jp])],
        [
            Tinf * c
            for c in poly_add(poly_mul(up[jp + 1], R[jp - 1]), poly_mul(up[jp - 1], R[jp + 1]))
        ],
    )
    W2 = poly_add(W2, [t1 * c for c in poly_mul(R[jp + 1], R[jp - 1])])
    q = int(float(h1))
    fac = arb("1.01") * (arb.const_euler() - h1).exp() * h1**q
    W = poly_add(W1, [t1 * c for c in W2])
    Z = poly_add([arb(0)] * q + list(F), [-fac * c for c in W])
    if not positive_on_interval(Z, h1, h2):
        return False, f"Fhat H^q - defect is not positive on [{float(h1):.0f}, {float(h2):.0f}]"
    return True, f"bracket B >= 0 dropped on [{float(h1):.0f}, {float(h2):.0f}]"


def tail_naive(jp, u, h1=H1):
    """F(N) > 0 for every N with H_N >= h1, bounding each E_k on its own -- the version WITHOUT the first-order
    cancellation.  It is used on the upper half [h2, infinity), where t is already small enough that the
    cancellation is not needed; `tail_ok` keeps the cancellation on [h1, h2], where it is.

    Write `E_k = Ehat_k + d_k` with `|d_k| <= t dt_k(H)`, `dt_k(H) = SUM_{s>=2} beta_s Ehat_{k-s}(H)` from
    `beta_conv`, and `T(N) <= T_inf (1 + 1.6 t)` (from `T/T_inf - 1 = 1/(2N+3+j') + 1/(N-j') + cross
    <= 1.51/N` once `j'/N <= 1e-5`, which holds because `H >= 20` means `N >= 2.7e8`).  Then, using
    `(a+d)^2 >= a^2 - 2 a dt t` and `t <= t_1`,
        F  >=  Fhat  -  t W,     W = 2 Ehat_j' dt_j' + T_inf Qt + 1.6 T_inf (P + t_1 Qt),
        Qt = Ehat_{j'+1} dt_{j'-1} + Ehat_{j'-1} dt_{j'+1} + t_1 dt_{j'+1} dt_{j'-1},   P = Ehat_{j'+1} Ehat_{j'-1}.
    `W >= 0` pointwise (every factor is a positive quantity), so multiplying by `H^q` with `q = [h1] <= h1` and
    using `e^{-H} H^q <= e^{-h1} h1^q` for `H >= h1 >= q` turns it into the single polynomial statement
        Fhat(H) H^q  -  1.01 e^{gamma-h1} h1^q W(H)  >  0     on [h1, infinity),
    of degree `2j'-2+q` against `2j'`, settled for all `H >= h1` at once by the Taylor shift at `h1`.
    """
    Ehat = {}
    for k in (jp - 1, jp, jp + 1):
        c = [arb(0)] * (k + 1)
        for m in range(min(k, len(u) - 1) + 1):
            c[k - m] = u[m] / arb(factorial(k - m))
        Ehat[k] = c
    C, t1 = beta_conv(jp + 1, u, h1)
    dt = {}
    for k in (jp - 1, jp, jp + 1):
        dt[k] = [C[k - i] / arb(factorial(i)) for i in range(k - 1)]
    G = poly_mul(Ehat[jp], Ehat[jp])
    P = poly_mul(Ehat[jp + 1], Ehat[jp - 1])
    Tinf = arb(fmpq(jp + 1, jp))
    F = [G[i] - (Tinf * P[i] if i < len(P) else arb(0)) for i in range(len(G))]
    # the two top coefficients vanish IDENTICALLY: [H^{2j'}] = (1/j'!)^2 - ((j'+1)/j')/((j'+1)!(j'-1)!) = 0 and
    # likewise at 2j'-1 (both are the statement that T_inf = (j'+1)/j' is exactly Newton's factor, i.e. that the
    # limit function e^{(H-gamma)x}/Gamma(1+x) is Laguerre-Polya); numerically they come out as balls around 0,
    # and a ball times H^{2j'} would be unbounded as H -> infinity, so they are verified against a tolerance and
    # then set to exact zero
    scale = max(float(abs(c).upper()) for c in F)
    for i in (2 * jp, 2 * jp - 1):
        if not (float(abs(F[i]).upper()) <= 1e-30 * scale):
            return False, f"top coefficient at degree {i} did not vanish: {F[i].str(4)}"
        F[i] = arb(0)
    Qt = poly_add(
        poly_add(poly_mul(Ehat[jp + 1], dt[jp - 1]), poly_mul(Ehat[jp - 1], dt[jp + 1])),
        [t1 * c for c in poly_mul(dt[jp + 1], dt[jp - 1])],
    )
    W = poly_add(
        [2 * c for c in poly_mul(Ehat[jp], dt[jp])],
        poly_add(
            [Tinf * c for c in Qt],
            [arb("1.6") * Tinf * c for c in poly_add(P, [t1 * c for c in Qt])],
        ),
    )
    q = int(float(h1))
    fac = arb("1.01") * (arb.const_euler() - h1).exp() * h1**q
    lhs = [arb(0)] * q + list(F)  # Fhat(H) * H^q
    Z = poly_add(lhs, [-fac * c for c in W])
    sh = shift_poly(Z, h1)
    neg = [i for i, c in enumerate(sh) if not (c.lower() >= 0)]
    if neg:
        return False, f"shifted coefficients of Fhat H^q - defect negative at degrees {neg[:5]}"
    return True, f"shift at h1 = {float(h1):.0f}, polynomial defect (t1 = {float(t1):.2e}, q = {q})"


H2_LADDER = (25, 30, 40, 55, 75, 100, 140)


def tail_ok(jp, u, h1=H1):
    """F(N) > 0 for every N with H_N >= h1, in two pieces that meet at an h2 chosen per index.

    The first-order defect enters F through the bracket `B` of `tail_bracket`, in which the `k^2` of the naive
    per-term bound cancels.  `B >= 0` on a range of `H` that GROWS with the index (measured: `B > 0` at every
    `H <= 100` for `j' >= 200`, while at `j' = 40` it turns negative near `H = 30`; the Poisson limit
    `B -> -2 H^{2k-2}/(k!)^2` says it must turn eventually, at `H` of order `k`).  So:

      * on `[h1, h2]` the bracket is verified non-negative and the first-order term is DROPPED, leaving only
        the `1.6 t` defect of `T` and the `t^2` remainder -- relative size 1.6 instead of `k^2/H^2`;
      * on `[h2, infinity)` the naive per-term bound suffices, because `t <= e^{gamma-h2}` has fallen far
        enough; that half is settled for all `H >= h2` at once by the pure Taylor shift.

    `h2` comes from a ladder, smallest first: a larger `h2` makes the lower half harder and the upper half
    easier.  Returning False means no rung worked, not that F <= 0.
    """
    last = ""
    for h2v in H2_LADDER:
        h2 = arb(h2v)
        if h2.lower() <= h1.lower():
            continue
        ok_hi, msg_hi = tail_naive(jp, u, h2)
        if not ok_hi:
            last = f"upper half at h2 = {h2v}: {msg_hi}"
            continue
        ok_lo, msg_lo = tail_bracket(jp, u, h1, h2)
        if not ok_lo:
            last = f"lower half [{float(h1):.0f}, {h2v}]: {msg_lo}"
            continue
        return True, f"{msg_lo}; naive shift from h2 = {h2v}"
    return False, last or "no rung of the h2 ladder worked"


EXACT_LOG = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "results",
    "sibuya_exact_pass_803_2026-09-03.txt",
)


def exact_log_ok(jpmax):
    """The low-N half of the coverage is an artifact, not a claim: check that the log of
    `sibuya_harmonic.exact_pass` exists, says ok = True, and covers every index this run asserts.

    Written after a debate agent found the verdict PRINTING "exact recursion by sibuya_harmonic.py" while the
    module was never imported and never run in this process.  A sentence in an f-string is not evidence.
    """
    if not os.path.exists(EXACT_LOG):
        return False, f"no log at {EXACT_LOG}"
    txt = open(EXACT_LOG, encoding="utf-8").read()
    m = re.search(
        r"every j' < (\d+) and every N in \[\d+ \+ j', 10\^(\d+)\]: ok = (True|False)", txt
    )
    if m is None:
        return False, "log does not carry the coverage line"
    jp_log, e_log, ok = int(m.group(1)), int(m.group(2)), m.group(3) == "True"
    if not ok:
        return False, "log says ok = False"
    if jp_log <= jpmax:
        return False, f"log covers j' < {jp_log}, this run needs j' <= {jpmax}"
    if 10**e_log < N_EXACT:
        return False, f"log covers N <= 10^{e_log}, this run assumes N <= {N_EXACT}"
    w = re.search(r"worst relative margin ([0-9.e+-]+)", txt)
    return True, f"j' < {jp_log}, N <= 10^{e_log}, worst rel {w.group(1) if w else '?'}"


def main():
    t0 = time.time()
    jpmax = JPMAX
    print(f"Sibuya (3.4), the window {JPMIN} <= j' <= {jpmax}, every N >= 1001 + j':", flush=True)
    u = u_coeffs(jpmax + 3)
    bad = []
    for jp in range(JPMIN, jpmax + 1):
        ok, msg = tail_ok(jp, u)
        if not ok:
            bad.append((jp, msg))
    print(
        f"  H-model tail (H >= {float(H1)}, i.e. N >= {float(H1.exp()):.2e}): {jpmax - JPMIN + 1 - len(bad)}/{jpmax - JPMIN + 1} indices certified by the shift"
        f"  [{time.time() - t0:.0f} s]",
        flush=True,
    )
    if bad:
        print(f"  FAILED indices (first 5): {bad[:5]}")
        sys.exit(1)
    if "--tail-only" in sys.argv:
        return
    ok_ex, msg_ex = exact_log_ok(jpmax)
    print(f"  exact recursion N <= {N_EXACT}: {'ok' if ok_ex else 'MISSING'}  {msg_ex}", flush=True)
    if not ok_ex:
        sys.exit(1)
    N_hi = int(float(H1.exp()))
    ok, msg = grid_ok(N_EXACT, N_hi, jpmax)
    print(f"  grid {N_EXACT} -> {N_hi}: {'ok' if ok else 'FAIL'}  {msg}", flush=True)
    if not ok:
        sys.exit(1)
    print()
    print(
        f"VERDICT: Sibuya (3.4) holds for {JPMIN} <= j' <= {jpmax} and every N (exact recursion N <= {N_EXACT}"
        f" from {os.path.basename(EXACT_LOG)}, grid to {N_hi}, H-model beyond): True  [{time.time() - t0:.0f} s]"
    )


if __name__ == "__main__":
    main()
