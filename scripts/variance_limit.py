"""PROVED: the limiting variance inequality behind the conjecture.

With A(v) = arctan(v)/v and B(v) = 1/(1+v^2),

    F(v) := 10 A(1-A) - (A-B)(5 + 4A - 4A^2)   >   0    for every v > 0.

This is the m -> infinity form of the statement

    Var(q) . (5n + 4 N Q)  >  4 N Q^2 ,   Q = qbar(1-qbar) ,

for the saddle probabilities q_k = b_k r/(1+b_k r) of the centred-square spectrum, which together
with  g >= (sigma_bin^2 - sigma^2)/(sigma^2 sigma_bin^2)  gives M(n,t) > 4/5.

Three regions, no gaps:
  (0, 0.2]   exact rational series in w = v^2, first term (352/945) w^3, alternating tail bounded
  [0.2, 3.5] certified interval bisection with arb
  [3.5, inf) analytic: F = 5A - 14A^2 + 4A^3 + B(5+4A-4A^2); the B-term is positive since
             5+4A-4A^2 >= 5 on [0,1], and 5-14A+4A^2 > 0 for A < (14-sqrt(116))/8 = 0.40372,
             while A is decreasing with A(3.5) = 0.36928.

Run:  python scripts/variance_limit.py
"""

from flint import fmpq, arb, ctx
import sys

ctx.prec = 200
K = 40
Z = fmpq(0)


def mul(a, b):
    r = [Z] * K
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if i + j >= K:
                break
            if bj:
                r[i + j] += ai * bj
    return r


def series_F():
    """F as an exact rational series in w = v^2."""
    A = [fmpq((-1) ** k, 2 * k + 1) for k in range(K)]  # arctan(v)/v
    B = [fmpq((-1) ** k) for k in range(K)]  # 1/(1+w)
    one = [Z] * K
    one[0] = fmpq(1)
    Q = mul(A, [one[i] - A[i] for i in range(K)])
    AB = [A[i] - B[i] for i in range(K)]
    T = mul(AB, [5 * one[j] + 4 * Q[j] for j in range(K)])
    return [10 * Q[i] - T[i] for i in range(K)]


def F(v):
    A = v.atan() / v
    B = 1 / (1 + v * v)
    return 10 * A * (1 - A) - (A - B) * (5 + 4 * A - 4 * A * A)


ok = True

# --- region 1: the series, 0 < v <= 0.2  (w <= 0.04) ---------------------------------
Fs = series_F()
s = next(i for i in range(K) if Fs[i] != 0)
G = [Fs[s + i] for i in range(K - s)]
assert s == 3 and G[0] == fmpq(352, 945), (s, G[0])


def fl(q):
    return float(arb(q.numer().str()) / arb(q.denom().str()))


W = fmpq(4, 100)
KK = 26
Cg = max(abs(fl(G[i])) / (i + 1) for i in range(KK))
tot = fmpq(0)
pw = W
for i in range(1, KK):
    tot += abs(G[i]) * pw
    pw *= W
w = arb(W.numer().str()) / arb(W.denom().str())
tail = arb(str(Cg + 0.001)) * w**KK * ((KK + 1) - KK * w) / ((1 - w) ** 2)
low = fl(G[0]) - fl(tot) - float(tail)
print(
    "(0, 0.2]   series: first term (352/945) w^3;  F/w^3 >= %.10f   tail <= %.2e   %s"
    % (low, float(tail), "OK" if low > 0 else "FAILED")
)
ok &= low > 0


# --- region 2: interval bisection on [0.2, 3.5] -------------------------------------
def sweep(lo, hi):
    stack = [(arb(lo), arb(hi))]
    boxes = 0
    while stack:
        a, b = stack.pop()
        boxes += 1
        if boxes > 2000000:
            return False, boxes
        iv = arb((a + b) / 2, ((b - a) / 2).mid())
        if F(iv).lower() > 0:
            continue
        if (b - a) < arb("1e-13"):
            return False, boxes
        mid = (a + b) / 2
        stack.append((a, mid))
        stack.append((mid, b))
    return True, boxes


good, boxes = sweep("0.2", "3.5")
print("[0.2, 3.5] interval bisection: %s   %d boxes" % ("PROVED F>0" if good else "FAILED", boxes))
ok &= good

# --- region 3: analytic, v >= 3.5 ----------------------------------------------------
A35 = arb("3.5").atan() / arb("3.5")
root = (14 - arb(116).sqrt()) / 8
print(
    "[3.5, inf)  analytic: A(3.5) = %.8f  <  (14-sqrt(116))/8 = %.8f   %s"
    % (float(A35), float(root), "OK" if A35 < root else "FAILED")
)
ok &= bool(A35 < root)

print()
print("THE LIMIT VARIANCE INEQUALITY IS %s for every v > 0." % ("PROVED" if ok else "NOT PROVED"))
sys.exit(0 if ok else 1)
