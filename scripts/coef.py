from flint import fmpq
import sys

K = 34  # series order in x = v^2
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
    assert a[0] != 0
    r = [Z] * K
    r[0] = 1 / a[0]
    for n in range(1, K):
        s = Z
        for j in range(1, n + 1):
            if j < len(a):
                s += a[j] * r[n - j]
        r[n] = -s / a[0]
    return r


def div(a, b):
    return mul(a, inv(b))


def comp(a, b):  # a(b(x)), b[0]=0
    r = [Z] * K
    p = [Z] * K
    p[0] = ONE
    for n in range(K):
        if a[n] != 0:
            for i in range(K):
                r[i] += a[n] * p[i]
        p = mul(p, b)
    return r


def revert(b):  # b[0]=0, b[1]!=0 ; return q with b(q(y))=y
    q = [Z] * K
    q[1] = 1 / b[1]
    for n in range(2, K):
        t = comp(b, q)
        q[n] -= t[n] / b[1]
    return q


def sh(a, k):  # multiply by x^k  (k may be negative)
    r = [Z] * K
    for i in range(K):
        j = i + k
        if 0 <= j < K:
            r[j] = a[i]
    return r


x = [Z] * K
x[1] = ONE  # x = v^2
one = [Z] * K
one[0] = ONE
# ---- I_n / v  as series in x  (all I_n are odd in v, so I_n/v is a series in x) ----
at = [Z] * K  # arctan(v)/v  = sum (-1)^k x^k/(2k+1)
for k in range(K):
    at[k] = fmpq((-1) ** k, 2 * k + 1)
J = [None] * 6
J[0] = one[:]  # I_0/v = 1
J[1] = at[:]  # I_1/v
opx = [Z] * K
opx[0] = ONE
opx[1] = ONE  # 1+x
for n in (1, 2, 3):
    inv_n = inv(opx)
    p = one[:]
    for _ in range(n):
        p = mul(p, inv_n)
    term = [c / (2 * n) for c in p]  # v/(2n(1+v^2)^n) / v
    J[n + 1] = [term[i] + fmpq(2 * n - 1, 2 * n) * J[n][i] for i in range(K)]
# s_k = (combination of I)/v
S2 = [J[1][i] - J[2][i] for i in range(K)]
S3 = [-J[1][i] + 3 * J[2][i] - 2 * J[3][i] for i in range(K)]
S4 = [J[1][i] - 7 * J[2][i] + 12 * J[3][i] - 6 * J[4][i] for i in range(K)]


# each ~ x * (...)  -> divide by x
def dropx(a):
    assert a[0] == 0, a[0]
    return a[1:] + [Z]


S2 = dropx(S2)
S3 = dropx(S3)
S4 = dropx(S4)
print("S2(0)=%s  S3(0)=%s  S4(0)=%s   (all should be 1/3)" % (S2[0], S3[0], S4[0]))
G = [S4[i] for i in range(K)]
G = div(G, [8 * c for c in mul(S2, S2)])
G2 = div([5 * c for c in mul(S3, S3)], [24 * c for c in mul(mul(S2, S2), S2)])
G = [G[i] - G2[i] for i in range(K)]
print("G(0) = %s   (should be -1/4)" % G[0])
sys.stdout.flush()
TH = [Z] * K  # theta = 1 - arctan(v)/v  = x/3 - x^2/5 + ...
for k in range(1, K):
    TH[k] = -at[k]
THt = dropx(TH)  # theta / x
print("theta/x at 0 = %s  (should be 1/3)" % THt[0])
num = [one[i] - TH[i] + mul(TH, TH)[i] for i in range(K)]
den = [12 * c for c in mul(THt, [one[i] - TH[i] for i in range(K)])]
B = div(num, den)  # (1-th+th^2)/(12*(theta/x)*(1-theta))  = -x*tau_bin
H = [G[i] + B[i] for i in range(K)]
print("H(0) = %s   (must be 0 : the poles cancel)" % H[0])
sys.stdout.flush()
TAU_x = dropx(H)  # tau as a series in x
# ---- convert to a series in theta ----
xoftheta = revert(TH)
TAU = comp(TAU_x, xoftheta)
print()
print("   tau(theta) = sum t_m theta^m   (exact rationals)")
for m in range(0, 12):
    print(
        "     t_%-2d = %-28s = %+.10f"
        % (
            m,
            TAU[m],
            float(TAU[m].numer().str()) / float(TAU[m].denom().str()) if TAU[m] != 0 else 0.0,
        )
    )
import json

open("tau_coeffs.json", "w").write(
    json.dumps([[TAU[m].numer().str(), TAU[m].denom().str()] for m in range(K)])
)
print("\n   wrote %d exact coefficients" % K)
