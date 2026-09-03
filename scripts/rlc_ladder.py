"""The RLC ladder: ratio log-concavity of p_j = e_j/C(N,j), index by index, for every m.

    RLC at index i :   p_{i+1}^3 p_{i-1}  >  p_i^3 p_{i+2}

With  M(n,1) > 4/5  (proved in base.py) this closes the conjecture at every index the ladder
reaches, because RLC makes M(n,t) increasing in t, putting its minimum at t = 1.

The rung at index i is a one-variable polynomial in m of degree exactly 16i + 6. Clearing
denominators to fmpz_poly and substituting m -> m + floor(i/2) + 2 makes every coefficient
non-negative with a positive constant term, which proves the rung for every admissible m.
Both laws -- degree 16i+6 and shift floor(i/2)+2 -- are exact at every index computed.

Unlike the earlier certificate ladder this statement carries NO CONSTANT. **It is NOT cheaper
at comparable indices** -- an earlier version of this docstring said "orders of magnitude
cheaper, 67 rungs in 14 s against ~100 s per rung", which extrapolated from the cheap early
rungs. Measured: rungs 1..144 take about 10 minutes together, and rung 145 alone takes over
49 minutes, while the M-ladder was doing 130 s per rung at index 627. At comparable indices this
ladder is SLOWER. Its value is the constant-free statement, not speed.

    Run:  python scripts/rlc_ladder.py [max index]
"""

import sys
from math import comb, factorial, gcd

from flint import fmpq, fmpq_poly, fmpz_poly

NE = int(sys.argv[1]) if len(sys.argv) > 1 else 24


# Faulhaber recursion: S_p(x) = sum_{k=1}^{x} k^p, one O(p) step per p instead of an
# O(p^2) Lagrange interpolation -- the setup was the bottleneck, not the ladder.
def faulhaber(P):
    S = [fmpq_poly([0, 1])]  # S_0(x) = x
    for p in range(1, P + 1):
        acc = fmpq_poly([1, 1]) ** (p + 1) - fmpq_poly([1])
        for j in range(p):
            acc = acc - S[j] * fmpq(comb(p + 1, j))
        S.append(acc / fmpq(p + 1))
    return S


# Incremental: rung i needs power sums only to order i+2, so nothing is paid up front and
# every rung is printed the moment it is proved.
FA = faulhaber(6)
Pw = [None]
E = [fmpq_poly([1])]
N = fmpq_poly([0, 2])


def grow(order):
    """Extend the Faulhaber table, the power sums and the elementary symmetric polynomials."""
    global FA
    if 2 * order >= len(FA):
        FA = faulhaber(2 * order + 2)
    while len(Pw) <= order:
        r = len(Pw)
        c = list(FA[2 * r].coeffs())
        scaled = fmpq_poly([c[j] * fmpq(2) ** j for j in range(len(c))])
        Pw.append((scaled - FA[2 * r] * fmpq(2) ** (2 * r)) * fmpq(2))
    while len(E) <= order:
        j = len(E)
        acc = fmpq_poly([0])
        for i2 in range(1, j + 1):
            acc = acc + (E[j - i2] * Pw[i2]) * fmpq((-1) ** (i2 - 1))
        E.append(acc / fmpq(j))


def C(k):
    r = fmpq_poly([1])
    for a in range(k):
        r = r * (N - fmpq(a))
    return r / fmpq(factorial(k))


print("  i  degree  deg-(16i+6)  shift  floor(i/2)+2  leading coefficient")
for i in range(1, NE + 1):
    grow(i + 2)
    Dp = E[i + 1] ** 3 * E[i - 1] * (C(i) ** 3 * C(i + 2)) - E[i] ** 3 * E[i + 2] * (
        C(i + 1) ** 3 * C(i - 1)
    )
    den = 1
    for c in Dp.coeffs():
        d = int(c.denom())
        den = den * d // gcd(den, d)
    Z = fmpz_poly([int(c * den) for c in Dp.coeffs()])
    g = 0
    for c in Z.coeffs():
        g = gcd(g, abs(int(c)))
    Z = fmpz_poly([int(c) // g for c in Z.coeffs()])
    sh = None
    for sv in range(
        1, i // 2 + 8
    ):  # the law says floor(i/2)+2; a fixed 60 silently failed at i=116
        cs = [int(c) for c in Z(fmpz_poly([sv, 1])).coeffs()]
        if all(c >= 0 for c in cs) and cs[0] > 0:
            sh = sv
            break
    print(
        f"{i:3d} {Dp.degree():6d} {Dp.degree() - (16 * i + 6):10d} "
        f"{sh!s:>6} {i // 2 + 2:10d}   {len(str(int(Z.coeffs()[-1])))} digits",
        flush=True,
    )
