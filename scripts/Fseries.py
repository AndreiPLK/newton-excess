from flint import fmpq
import sys

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
x = [Z] * K
x[1] = ONE
a = [Z] * K
for k in range(K):
    a[k] = fmpq((-1) ** k, 2 * k + 1)  # arctan(v)/v  in x=v^2
opx = [Z] * K
opx[0] = ONE
opx[1] = ONE  # 1+x
th = [one[i] - a[i] for i in range(K)]  # theta = 1 - a
p = [mul(a, opx)[i] - one[i] for i in range(K)]  # p = a(1+x) - 1  ~ (2/3)x
pt = dropx(p)  # p/x   -> 2/3 at 0
tht = dropx(th)  # theta/x -> 1/3 at 0
print("   p/x at 0 = %s (want 2/3)     theta/x at 0 = %s (want 1/3)" % (pt[0], tht[0]))
sys.stdout.flush()
# LHS = 2[-1 + 2x/(p(1+x))] (1+x)^2 / p^2
#     = (1/x^2) * 2[ -1 + 2/(pt (1+x)) ] (1+x)^2 / pt^2
inner = [2 * inv(mul(pt, opx))[i] - one[i] for i in range(K)]  # -1 + 2/(pt(1+x))   -> -1+3 = 2
LH = mul([2 * c for c in inner], mul(mul(opx, opx), inv(mul(pt, pt))))  # times 1/x^2
# RHS = 1/theta^2 - 1/a^2 = (1/x^2)/tht^2  -  1/a^2
RH1 = inv(mul(tht, tht))  # times 1/x^2
RH2 = inv(mul(a, a))  # no 1/x^2
# F = RHS - LHS = (1/x^2)(RH1 - LH) - RH2
G = [RH1[i] - LH[i] for i in range(K)]
print("   G[0] = %s   G[1] = %s   (both must vanish for F to be analytic)" % (G[0], G[1]))
sys.stdout.flush()
G2 = dropx(dropx(G))  # (RH1-LH)/x^2
F = [G2[i] - RH2[i] for i in range(K)]
print()
print("   F(x) = RHS - LHS  as a series in x = v^2 :   ARE ALL COEFFICIENTS >= 0 ?")
neg = 0
for k in range(K - 4):
    if F[k] < 0:
        neg += 1
    if k < 16 or k >= K - 8:
        val = float(F[k].numer().str()) / float(F[k].denom().str()) if F[k] != 0 else 0.0
        print(
            "     F_%-2d = %-30s = %+.10f %s" % (k, F[k], val, "  <-- NEGATIVE" if F[k] < 0 else "")
        )
print()
print("   negative coefficients among the first %d : %d" % (K - 4, neg))
print("   F_0 = %s   (176/175 = %s)" % (F[0], fmpq(176, 175)))
