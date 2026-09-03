from flint import fmpq, fmpq_poly
from math import comb
import sys


def esyms(b, t):
    e = [fmpq(1)]
    for v in b:
        e = [e[q] + (v * e[q - 1] if q else 0) for q in range(len(e))] + [e[-1] * v]
        if len(e) > t + 1:
            e = e[: t + 1]
    return e


def true_p(nn, i):
    b = [fmpq((nn - 2 * k) ** 2) for k in range(1, nn)]
    return esyms(b, i)[i] / fmpq(comb(nn - 1, i))


def p_poly(i):
    deg = 2 * i
    xs = list(range(i + 2, i + 2 + deg + 1))
    ys = [true_p(x, i) for x in xs]
    P = fmpq_poly([0])
    for a in range(len(xs)):
        L = fmpq_poly([1])
        den = fmpq(1)
        for bb in range(len(xs)):
            if bb == a:
                continue
            L = L * fmpq_poly([-xs[bb], 1])
            den *= fmpq(xs[a] - xs[bb])
        P = P + L * (ys[a] / den)
    for x in range(i + 2 + deg + 1, i + 2 + deg + 13):
        assert P(fmpq(x)) == true_p(x, i)
    return P


print("   SHIFT CERTIFICATES for the ladder in i")
print("   smallest c with all coefficients of W_i(n+c) nonnegative, against the regime bound 2i+5")
print("     i    deg    c      2i+5    certificate covers the regime?     c - i")
P = {}
for i in range(61, 141):
    for j in (i, i + 1, i + 2, i + 3):
        if j not in P:
            P[j] = p_poly(j)
    W = P[i + 2] * P[i + 2] * P[i + 2] * P[i] - P[i + 1] * P[i + 1] * P[i + 1] * P[i + 3]
    c = None
    for cc in range(0, 3 * i + 30):
        Q = W(fmpq_poly([cc, 1]))
        if all(Q.coeffs()[q] >= 0 for q in range(Q.degree() + 1)):
            c = cc
            break
    ok = c is not None and c <= 2 * i + 5
    print(
        "   %3d   %5d   %-5s  %4d       %-4s                     %s"
        % (
            i,
            W.degree(),
            str(c),
            2 * i + 5,
            "YES" if ok else "no",
            (c - i) if c is not None else "?",
        )
    )
    sys.stdout.flush()
