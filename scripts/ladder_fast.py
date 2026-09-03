from flint import fmpq, fmpq_poly, fmpz_poly, fmpz
from math import comb, gcd
import sys, time


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
    return [
        (faulhaber(2 * k, B)(fmpq_poly([0, 2])) - faulhaber(2 * k, B) * fmpq(4**k)) * 2
        for k in range(1, K + 1)
    ]


def binomM(i):
    P = fmpq_poly([1])
    for r in range(i):
        P = P * fmpq_poly([-r, 2])
    f = fmpq(1)
    for r in range(1, i + 1):
        f *= r
    return P * (fmpq(1) / f)


def to_z(P):
    "clear denominators: return (fmpz_poly, positive rational scale) with P = scale * Z"
    cs = [P.coeffs()[q] for q in range(P.degree() + 1)]
    den = 1
    for c in cs:
        d = int(c.denom())
        den = den * d // gcd(den, d)
    Z = fmpz_poly([int((c * fmpq(den)).numer()) for c in cs])
    g = 0
    for q in range(Z.degree() + 1):
        g = gcd(g, abs(int(Z.coeffs()[q])))
    if g > 1:
        Z = fmpz_poly([int(Z.coeffs()[q]) // g for q in range(Z.degree() + 1)])
    return Z, fmpq(g, den)


IMAX = int(sys.argv[1])
START = int(sys.argv[2])
K = IMAX + 4
t0 = time.time()
B = bernoulli(2 * K + 2)
ps = powersums(K, B)
print("setup %.1fs" % (time.time() - t0))
sys.stdout.flush()
e = [fmpq_poly([1])]
PZ = {}


def getpz(i):
    while len(e) <= i:
        m = len(e)
        s = fmpq_poly([0])
        for k in range(1, m + 1):
            t = e[m - k] * ps[k - 1]
            s = s + t if (k % 2 == 1) else s - t
        e.append(s * (fmpq(1) / fmpq(m)))
    if i not in PZ:
        q, r = divmod(e[i], binomM(i))
        assert r == fmpq_poly([0])
        PZ[i] = to_z(q)
    return PZ[i]


print("   INTEGER LADDER   i    deg W    c_M    covers n <=    time/index")
tprev = time.time()
for i in range(START, IMAX + 1):
    (A0, s0) = getpz(i)
    (A1, s1) = getpz(i + 1)
    (A2, s2) = getpz(i + 2)
    (A3, s3) = getpz(i + 3)
    # W = s2^3 s0 * A2^3 A0 - s1^3 s3 * A1^3 A3 ; put over a common rational factor
    L = s2**3 * s0
    R = s1**3 * s3
    # scale both by the lcm of denominators / gcd of numerators -> keep integers
    num = L / R
    a = int(num.numer())
    b = int(num.denom())
    WZ = A2**3 * A0 * fmpz(a) - A1**3 * A3 * fmpz(b)
    pred = (i + 5) // 2 if i % 2 else (i + 2) // 2
    c = None
    for cc in [pred] + [pred + d for d in (1, -1, 2, -2, 3, -3)]:
        if cc < 0:
            continue
        Q = WZ(fmpz_poly([cc, 1]))
        if all(Q.coeffs()[q] >= 0 for q in range(Q.degree() + 1)):
            c = cc
            break
    now = time.time()
    print(
        "   %4d  %6d   %5s   n <= %5d   %6.1fs" % (i, WZ.degree(), str(c), 2 * i + 5, now - tprev)
    )
    tprev = now
    sys.stdout.flush()
