"""The certificate's leading coefficient in closed form, and where 176/175 comes from.

Reproduces the lemma of 2 September 2026:

    L_i := [n^(8i+10)] ( p_{i+2}^3 p_i - p_{i+1}^3 p_{i+3} ) = (176/175) . 2^(8i+10) / 3^(4i+6)

together with the two ingredients that make it a derivation rather than an observation:

    a_i = (4/3)^i                      (leading coefficients; Newton's EQUALITY case,
                                        which is why the top two orders cancel)
    u_i = -i(i-1)/5                    (first correction, from f_2/f_1^2 and the binomial)

Exact rational arithmetic throughout (python-flint fmpq / fmpq_poly). No floating point.
Run:  python scripts/leading.py
"""

from flint import fmpq, fmpq_poly
from math import comb
import sys

IMAX = 24


def bernoulli(nmax):
    B = [fmpq(0)] * (nmax + 1)
    B[0] = fmpq(1)
    for m in range(1, nmax + 1):
        s = fmpq(0)
        for k in range(m):
            s += fmpq(comb(m + 1, k)) * B[k]
        B[m] = -s / fmpq(m + 1)
    return B


def faulhaber(p, B):
    c = [fmpq(0)] * (p + 2)
    for j in range(p + 1):
        c[p + 1 - j] += fmpq(comb(p + 1, j)) * B[j]
    return (fmpq_poly(c) * (fmpq(1) / fmpq(p + 1)))(fmpq_poly([1, 1]))


def powersums(K, B):
    """S_j(m) = 2 * sum_{k=1..m} (2k-1)^(2j), as polynomials in m."""
    return [
        (faulhaber(2 * k, B)(fmpq_poly([0, 2])) - faulhaber(2 * k, B) * fmpq(4**k)) * 2
        for k in range(1, K + 1)
    ]


def binom_poly(i):
    """C(2m, i) as a polynomial in m."""
    P = fmpq_poly([1])
    for r in range(i):
        P = P * fmpq_poly([-r, 2])
    f = fmpq(1)
    for r in range(1, i + 1):
        f *= r
    return P * (fmpq(1) / f)


K = IMAX + 4
B = bernoulli(2 * K + 2)
ps = powersums(K, B)
e = [fmpq_poly([1])]
_cache = {}


def p_poly(i):
    """p_i = e_i / C(2m, i), exactly, as a polynomial in m of degree 2i."""
    while len(e) <= i:
        mm = len(e)
        s = fmpq_poly([0])
        for k in range(1, mm + 1):
            term = e[mm - k] * ps[k - 1]
            s = s + term if (k % 2 == 1) else s - term
        e.append(s * (fmpq(1) / fmpq(mm)))
    if i not in _cache:
        q, r = divmod(e[i], binom_poly(i))
        assert r == fmpq_poly([0]), "e_i is not divisible by C(2m,i) -- construction is wrong"
        _cache[i] = q
    return _cache[i]


f1, f2, f3, g1 = fmpq(8, 3), fmpq(-16, 5), fmpq(128, 21), fmpq(-2, 3)


def u_pred(i):
    return -fmpq(i * (i - 1), 5)


def w_pred(i):
    a2 = (
        fmpq(i * (i - 1) * (i - 2) * (i - 3), 2) * (f2 / f1**2) ** 2
        + fmpq(i * (i - 1) * (i - 2)) * (f3 / f1**3)
        + fmpq(i) * (g1 / f1)
    )
    S = fmpq(i * (i - 1), 2)
    Q = fmpq((i - 1) * i * (2 * i - 1), 6)
    return a2 + Q / 8 + S * S / 8 + fmpq(i * (i - 1)) * (f2 / f1**2) * (S / 2)


print("the four constants of the expansion:  f1=%s  f2=%s  f3=%s  g1=%s" % (f1, f2, f3, g1))
print("S_1(m) = %s   -- only odd powers of m, which is what makes the expansion work" % ps[0])
print()

bad = 0
for i in range(IMAX + 1):
    P = p_poly(i)
    d = 2 * i
    assert P.degree() == d, (i, P.degree())
    a = P.coeffs()[d]
    if a != fmpq(4, 3) ** i:
        print("  a_%d MISMATCH" % i)
        bad += 1
    if d >= 1 and P.coeffs()[d - 1] / a != u_pred(i):
        print("  u_%d MISMATCH" % i)
        bad += 1
    if d >= 2 and P.coeffs()[d - 2] / a != w_pred(i):
        print("  w_%d MISMATCH" % i)
        bad += 1
print(
    "a_i = (4/3)^i, u_i = -i(i-1)/5, w_i from the four constants:  %s  (i = 0..%d)"
    % ("all exact" if bad == 0 else "%d MISMATCHES" % bad, IMAX)
)


def W1_minus_W2(i):
    u, w = u_pred, w_pred
    W1 = 3 * w(i + 2) + 3 * u(i + 2) ** 2 + 3 * u(i + 2) * u(i) + w(i)
    W2 = 3 * w(i + 1) + 3 * u(i + 1) ** 2 + 3 * u(i + 1) * u(i + 3) + w(i + 3)
    return W1 - W2


vals = {W1_minus_W2(i) for i in range(60)}
print("W1 - W2 over i = 0..59 :  %s   -- the i-dependence cancels" % vals)

bad = 0
for i in range(IMAX + 1):
    Z = (p_poly(i + 2) ** 3) * p_poly(i) - (p_poly(i + 1) ** 3) * p_poly(i + 3)
    assert Z.degree() == 8 * i + 10, (i, Z.degree())
    L = Z.coeffs()[Z.degree()]
    pred = fmpq(176, 175) * fmpq(2) ** (8 * i + 10) / fmpq(3) ** (4 * i + 6)
    if L != pred:
        print("  L_%d MISMATCH: %s vs %s" % (i, L, pred))
        bad += 1
print(
    "L_i = (176/175) 2^(8i+10)/3^(4i+6) against the exact certificate:  %s  (i = 0..%d)"
    % ("all exact" if bad == 0 else "%d MISMATCHES" % bad, IMAX)
)
print()
print("176/175 is also the proved minimum gap of the limiting inequality (see endgame.py),")
print("reached there at theta = 0 by a computation that shares no code with this one.")
sys.exit(1 if bad else 0)
