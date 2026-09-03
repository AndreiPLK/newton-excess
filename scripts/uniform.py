"""The uniform expansion, and the single bound the conjecture now rests on.

    n g(n,t)  =  H(theta)  +  [ H(theta) + (1/2) R''(theta) ] / N  +  r ,

with g = -Delta^2 log p, theta = t/N, N = n-1, H the limit shape and
R(theta) = log[ s(theta)/(theta(1-theta)) ].

Both leading coefficients are proved positive elsewhere in this package:
  H >= 4/5              (variance_limit.py, and Section 8 of the paper)
  (1/2)R'' >= 0.1021555 (majorant over exact rational Taylor coefficients)
so  n g >= 4/5 + 0.9021555/N - |r| , and since M = n(e^g - 1) > n g the conjecture
follows from  |N^2 r| <= 0.9021555 N .

This script checks that bound at every index of every odd n in a range, and prints the
supremum. Exact rational arithmetic for the spectrum; certified interval arithmetic (arb)
for H and R''.   Run:  python scripts/uniform.py [nmax]
"""

from flint import fmpq, arb, ctx
from math import comb
import sys

ctx.prec = 250
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 101
NMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 7


def A(q):
    return arb(q.numer().str()) / arb(q.denom().str())


def theta_of_v(v):
    return arb(1) - v.atan() / v


_vc = {}


def v_of_theta(th):
    k = str(th)
    if k in _vc:
        return _vc[k]
    lo, hi = arb("1e-16"), arb(120)
    for _ in range(220):
        mid = (lo + hi) / 2
        if theta_of_v(mid) < th:
            lo = mid
        else:
            hi = mid
        lo, hi = arb(lo.mid()), arb(hi.mid())
    r = arb(((lo + hi) / 2).mid())
    _vc[k] = r
    return r


def H(th):
    v = v_of_theta(th)
    t2 = theta_of_v(v)
    return arb(2) / ((1 - t2) - 1 / (1 + v * v)) - 1 / (t2 * (1 - t2))


_h = arb(1) / arb(100000)


def Rpp(th):
    def L(x):
        v = v_of_theta(x)
        t2 = theta_of_v(v)
        s = ((1 - t2) - 1 / (1 + v * v)) / 2
        return (s / (x * (1 - x))).log()

    return (L(th - _h) - 2 * L(th) + L(th + _h)) / (_h * _h)


ALLOW = 0.9021555  # a proved lower bound for H + (1/2)R''
worst, wn, wt, cells, fails = 0.0, None, None, 0, 0

for n in range(NMIN, NMAX + 1, 2):
    N = n - 1
    m = N // 2
    e = [fmpq(1)]
    for k in range(1, m + 1):
        c = fmpq((2 * k - 1) ** 2)
        for _ in range(2):
            e = [e[q] + (c * e[q - 1] if q else 0) for q in range(len(e))] + [e[-1] * c]
    p = [e[j] / fmpq(comb(N, j)) for j in range(N + 1)]
    for t in range(1, N // 2 + 1):
        th = arb(t) / arb(N)
        ng = arb(n) * A(p[t] ** 2 / (p[t - 1] * p[t + 1])).log()
        r = ng - H(th) - (H(th) + Rpp(th) / 2) / arb(N)
        v = abs(float(arb(N) * arb(N) * r))
        cells += 1
        if v > worst:
            worst, wn, wt = v, n, t
        if v > ALLOW * N:
            fails += 1
            print("  FAIL n=%d t=%d : |N^2 r| = %.4f > %.4f" % (n, t, v, ALLOW * N))

print("cells checked            : %d   (odd n from %d to %d, every index)" % (cells, NMIN, NMAX))
print("worst |N^2 r|            : %.5f   at n = %d, t = %d" % (worst, wn, wt))
print("requirement |N^2 r| <= %.7f N   failures : %d" % (ALLOW, fails))
print()
print("the allowance grows linearly in N while |N^2 r| decreases, so the margin widens;")
print("n = 5 is the single exception and is proved directly (base.py), n = 3 is where the")
print("conjecture is false.")
sys.exit(1 if fails else 0)
