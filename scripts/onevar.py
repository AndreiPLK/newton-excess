from flint import arb, ctx
import sys

ctx.prec = 600


def vof(th):
    lo, hi = arb(1) / arb(10**9), arb(10**7)
    for _ in range(300):
        mid = (lo * hi).sqrt()
        if 1 - mid.atan() / mid < th:
            lo = mid
        else:
            hi = mid
    return (lo * hi).sqrt()


def LHS_RHS(th):
    v = vof(th)
    v2 = v * v
    Q = (1 - th) - 1 / (1 + v2)
    Qp = -1 + 2 * v2 / (Q * (1 + v2) * (1 + v2))
    L = 2 * Qp / (Q * Q)
    R = 1 / (th * th) - 1 / ((1 - th) * (1 - th))
    return L, R, Q, Qp


print("   THE ONE-VARIABLE INEQUALITY :   2 Q'/Q^2  <=  1/theta^2 - 1/(1-theta)^2   on (0, 1/2)")
print("     theta       LHS            RHS           RHS - LHS        Q         Q'")
bad = 0
tot = 0
xs = [arb(k) / arb(1000) for k in list(range(1, 60, 4)) + list(range(60, 501, 20))]
for th in xs:
    L, R, Q, Qp = LHS_RHS(th)
    tot += 1
    ok = R - L > 0
    if not ok:
        bad += 1
    print(
        "     %.3f   %+13.6f  %+13.6f   %+13.6f   %8.5f  %+8.5f  %s"
        % (
            float(th),
            float(L),
            float(R),
            float(R - L),
            float(Q),
            float(Qp),
            "" if ok else "   <-- VIOLATED",
        )
    )
    sys.stdout.flush()
print()
print("   points tested: %d      violations: %d" % (tot, bad))
print()
print("   behaviour as theta -> 0  (both sides blow up like 1/theta^2; the O(1) term decides)")
for k in (1, 2, 5, 10, 20, 50, 100):
    th = arb(k) / arb(100000)
    L, R, Q, Qp = LHS_RHS(th)
    print(
        "     theta=%.5f   LHS-1/th^2 = %+12.6f    RHS-1/th^2 = %+12.6f   RHS-LHS = %+.6f"
        % (float(th), float(L - 1 / (th * th)), float(R - 1 / (th * th)), float(R - L))
    )
    sys.stdout.flush()
