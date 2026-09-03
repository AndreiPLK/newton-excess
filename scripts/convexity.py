"""H'' > 0 on the whole of [0, 1/2] : the limit shape is convex. PROVED, two ranges that overlap.

This is the leading-order form of the last open lemma, `Delta^4 log p < 0`, because
`Delta^4 log p = -Delta^2 g` and `g = H(theta)/n + O(1/n^2)`.

    theta in [0, 0.04]      exact rational Taylor series of H''; 72 coefficients, all
                            non-negative, c_0 = 352/175 exactly, tail below 1.3e-87
    theta in [0.038, 0.5]   certified interval bisection in v, 1 267 758 boxes

The ranges overlap, so the two together cover [0, 1/2].

The series needs both vanishing factors divided out first -- D = A - B = w E with E(0) = 2/3,
and 1 - A = w R with R(0) = 1/3 -- after which H = (1/w)[2/E - 1/(A R)] is regular at w = 0.
Reverting theta = w R and composing gives H(theta) in Q[[theta]]; two derivatives give H''.

Exact rational arithmetic for the series, certified interval arithmetic for the sweep; no
floating point enters a comparison.
"""

import sys
import time

from flint import arb, ctx, fmpq, fmpq_poly

ctx.prec = 200
NT = 76
VSTAR = arb("2.3311224")  # theta(v*) >= 1/2
VLO = arb("0.35")  # theta(VLO) = 0.03807 < 0.04, so the ranges overlap


def tr(p):
    return fmpq_poly(list(p.coeffs())[:NT])


def inv(p):
    r = fmpq_poly([1 / p.coeffs()[0]])
    for _ in range(9):
        r = tr(r * (fmpq_poly([2]) - tr(p * r)))
    return r


def sh(p):
    return fmpq_poly(list(p.coeffs())[1:])


def comp(f, g):
    out = fmpq_poly([0])
    pw = fmpq_poly([1])
    for k in range(NT):
        c = f.coeffs()[k] if k < len(f.coeffs()) else fmpq(0)
        out = out + pw * c
        pw = tr(pw * g)
    return tr(out)


def hpp_series():
    A = fmpq_poly([fmpq((-1) ** k, 2 * k + 1) for k in range(NT)])
    B = fmpq_poly([fmpq((-1) ** k) for k in range(NT)])
    E = sh(A - B)
    R = sh(fmpq_poly([1]) - A)
    Pq = tr(A * R)
    H_w = sh(tr(fmpq_poly([2]) * inv(E)) - inv(Pq))
    th_w = tr(fmpq_poly([0, 1]) * R)
    w_th = fmpq_poly([0, 3])
    for _ in range(10):
        w_th = tr(
            w_th - tr(tr(comp(th_w, w_th) - fmpq_poly([0, 1])) * inv(comp(th_w.derivative(), w_th)))
        )
    return comp(H_w, w_th).derivative().derivative()


def Hpp(v):
    A = v.atan() / v
    B = 1 / (1 + v * v)
    D = A - B
    T = D / v
    Dp = -D / v + 2 * v * B * B
    Tp = -2 * D / (v * v) + 2 * B * B
    Dpp = -Dp / v + D / (v * v) + 2 * B * B - 8 * v * v * B * B * B
    P = A * (1 - A)
    Pp = -T * (1 - 2 * A)
    Ppp = -Tp * (1 - 2 * A) - 2 * T * T
    Hv = -2 * Dp / (D * D) + Pp / (P * P)
    Hvv = -2 * Dpp / (D * D) + 4 * Dp * Dp / (D**3) + Ppp / (P * P) - 2 * Pp * Pp / (P**3)
    return (Hvv * T - Hv * Tp) / (T**3)


def sweep(lo, hi, maxd=70):
    stack = [(lo, hi, 0)]
    boxes = 0
    while stack:
        a, b, d = stack.pop()
        boxes += 1
        if Hpp(arb((a + b) / 2, (b - a) / 2)) > 0:
            continue
        if d > maxd:
            return False, boxes
        mid = (a + b) / 2
        stack.append((a, mid, d + 1))
        stack.append((mid, b, d + 1))
    return True, boxes


# ---- 1. the edge, by exact series -------------------------------------------------
cs = list(hpp_series().coeffs())
clean = 0
for k in range(len(cs)):
    if 0.9 <= float(cs[k]) / ((k + 1) * (k + 2)) <= 1.1:
        clean = k
    else:
        break
allpos = all(cs[k] >= 0 for k in range(clean + 1))
TH = fmpq(4, 100)
tail = sum(fmpq(10) ** 10 * (k + 1) * (k + 2) * TH**k for k in range(clean + 1, clean + 400))
edge_lo = float(cs[0]) - float(tail)

print("edge  theta in [0, 0.04], exact rational Taylor series")
print("  coefficients clean through k = %d   all non-negative: %s" % (clean, allpos))
print("  c_0 = %s  (= 352/175: %s)" % (cs[0], cs[0] == fmpq(352, 175)))
print("  tail with the wild bound |c_k| <= 1e10 (k+1)(k+2) : %.2e" % float(tail))
print("  H'' >= %.10f on the edge" % edge_lo)

# ---- 2. the bulk, by certified bisection -------------------------------------------
print()
print("bulk  theta in [0.03807, 0.5], certified interval bisection in v")
total, allok, lo = 0, True, VLO
t0 = time.time()
for c in ["0.4", "0.45", "0.5", "0.6", "0.8", "1.0"]:
    hi = arb(c)
    ok, n = sweep(lo, hi)
    total += n
    allok = allok and ok
    lo = hi
ok, n = sweep(lo, VSTAR)
total += n
allok = allok and ok
print("  boxes %d, %.1f s, certified: %s" % (total, time.time() - t0, allok))

ok_all = allpos and edge_lo > 0 and allok and float(1 - VLO.atan() / VLO) < 0.04
print()
print("ranges overlap at theta in [0.03807, 0.04] :", float(1 - VLO.atan() / VLO) < 0.04)
print("VERDICT: H'' > 0 on [0, 1/2] :", ok_all)
sys.exit(0 if ok_all else 1)
