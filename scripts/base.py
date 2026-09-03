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


p0, p1, p2 = p_poly(0), p_poly(1), p_poly(2)
print("   p_0 = %s" % p0)
print("   p_1 = %s" % p1)
print("   p_2 = %s" % p2)
# M(n,1) = n ( p_1^2/(p_0 p_2) - 1 ) > 4/5   <=>  5 n (p_1^2 - p_0 p_2) - 4 p_0 p_2 > 0
X = fmpq_poly([0, 5]) * (p1 * p1 - p0 * p2) - fmpq_poly([4]) * (p0 * p2)
print()
print("   X(n) = 5n(p_1^2 - p_0 p_2) - 4 p_0 p_2      must be > 0 for every odd n >= 3")
print("   X = %s" % X)
print("   degree %d" % X.degree())
for c in range(0, 12):
    Q = X(fmpq_poly([c, 1]))
    cs = [Q.coeffs()[q] for q in range(Q.degree() + 1)]
    nneg = sum(1 for v in cs if v < 0)
    if nneg == 0:
        print(
            "   shift c = %d : ALL COEFFICIENTS NONNEGATIVE  =>  M(n,1) > 4/5 for every n >= %d"
            % (c, c)
        )
        print("   coefficients: %s" % [str(v) for v in cs])
        break
    else:
        print("   shift c = %d : %d negative" % (c, nneg))
    sys.stdout.flush()
print()
for nn in (3, 5, 7, 9, 11, 21, 41, 101):
    m = nn * (true_p(nn, 1) ** 2 / (true_p(nn, 0) * true_p(nn, 2)) - 1)
    print(
        "   n=%3d   M(n,1) = %s = %.10f   > 4/5 : %s"
        % (nn, m, float(m.numer().str()) / float(m.denom().str()), m > fmpq(4, 5))
    )
