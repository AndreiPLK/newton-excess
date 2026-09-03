from flint import arb, ctx
import sys

ctx.prec = 1200
sys.setrecursionlimit(20000)


def F(v):
    A = v.atan()
    a = A / v
    th = 1 - a
    x = v * v
    p = a * (1 + x) - 1
    Q = p / (1 + x)
    Qp = -1 + 2 * x / (p * (1 + x))
    return (1 / (th * th) - 1 / (a * a)) - 2 * Qp / (Q * Q)


def iv(lo, hi):
    return arb(lo).union(arb(hi))


def certify(lo, hi, depth=0):
    f = F(iv(lo, hi))
    if f > 0:
        return True, 1
    if depth > 60:
        return False, 1
    mid = (lo + hi) / 2
    o1, b1 = certify(lo, mid, depth + 1)
    if not o1:
        return False, b1
    o2, b2 = certify(mid, hi, depth + 1)
    return o2, b1 + b2


VSTAR = 2.3311223
for LOW in (0.5, 0.3, 0.2, 0.15):
    ok, boxes = certify(LOW, VSTAR)
    print(
        "   interval arithmetic on v in [%.2f, %.6f] : %s   boxes %d"
        % (LOW, VSTAR, "PROVED F>0" if ok else "failed", boxes)
    )
    sys.stdout.flush()
    if not ok:
        break
