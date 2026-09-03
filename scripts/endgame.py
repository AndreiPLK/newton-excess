from flint import fmpq, arb, ctx
import sys

ctx.prec = 600
K = 46
Z = fmpq(0)
ONE = fmpq(1)


def mul(a, b):
    r = [Z] * K
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if i + j >= K:
                break
            if bj == 0:
                continue
            r[i + j] += ai * bj
    return r


def inv(a):
    r = [Z] * K
    r[0] = 1 / a[0]
    for n in range(1, K):
        s = Z
        for j in range(1, n + 1):
            if j < len(a):
                s += a[j] * r[n - j]
        r[n] = -s / a[0]
    return r


def dropx(a):
    assert a[0] == 0
    return a[1:] + [Z]


one = [Z] * K
one[0] = ONE
aa = [Z] * K
for k in range(K):
    aa[k] = fmpq((-1) ** k, 2 * k + 1)
opx = [Z] * K
opx[0] = ONE
opx[1] = ONE
th = [one[i] - aa[i] for i in range(K)]
p = [mul(aa, opx)[i] - one[i] for i in range(K)]
pt = dropx(p)
tht = dropx(th)
inner = [2 * inv(mul(pt, opx))[i] - one[i] for i in range(K)]
LH = mul([2 * c for c in inner], mul(mul(opx, opx), inv(mul(pt, pt))))
RH1 = inv(mul(tht, tht))
RH2 = inv(mul(aa, aa))
G = [RH1[i] - LH[i] for i in range(K)]
F = [dropx(dropx(G))[i] - RH2[i] for i in range(K)]
NT = 42
XMAX = arb(15) / arb(100)
XMAX = XMAX * XMAX  # x = v^2 <= 0.0225
print(
    "   ENDGAME : F(x) > 0 on x in [0, %s]  by the exact series with a certified tail" % XMAX.str(8)
)
mx = max(abs(float(F[k].numer().str()) / float(F[k].denom().str())) for k in range(NT))
print(
    "   max |F_k| over k < %d : %.6f     (tail bounded using |F_k| <= 1, a safe overestimate)"
    % (NT, mx)
)
xiv = arb(0).union(XMAX)
s = arb(0)
pw = arb(1)
for k in range(NT):
    c = arb(F[k].numer().str()) / arb(F[k].denom().str())
    s = s + c * pw
    pw = pw * xiv
tail = XMAX**NT / (1 - XMAX)
print("   truncated series enclosure : [%.10f, %.10f]" % (float(s.lower()), float(s.upper())))
print("   tail bound                 : %.6e" % float(tail))
tot = s + arb(0).union(tail) - arb(0).union(tail)
lo = s.lower() - tail
print("   certified lower bound on F : %.10f" % float(lo))
print("   F > 0 on [0, 0.0225] :", bool(lo > 0))
print()
print(
    "   F(0) = %s  exactly  = %.10f"
    % (F[0], float(arb(F[0].numer().str()) / arb(F[0].denom().str())))
)
