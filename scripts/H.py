from flint import fmpq
import sys

K = 42
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


def div(a, b):
    return mul(a, inv(b))


def comp(a, b):
    r = [Z] * K
    p = [Z] * K
    p[0] = ONE
    for n in range(K):
        if a[n] != 0:
            for i in range(K):
                r[i] += a[n] * p[i]
        p = mul(p, b)
    return r


def revert(b):
    q = [Z] * K
    q[1] = 1 / b[1]
    for n in range(2, K):
        t = comp(b, q)
        q[n] -= t[n] / b[1]
    return q


def dropx(a):
    assert a[0] == 0
    return a[1:] + [Z]


one = [Z] * K
one[0] = ONE
at = [Z] * K
for k in range(K):
    at[k] = fmpq((-1) ** k, 2 * k + 1)  # arctan(v)/v  in x = v^2
TH = [Z] * K
for k in range(1, K):
    TH[k] = -at[k]  # theta(x)
opx = [Z] * K
opx[0] = ONE
opx[1] = ONE  # 1 + x  = 1 + v^2
# bracket = (1 - theta) - 1/(1+x)
br = [(one[i] - TH[i]) - inv(opx)[i] for i in range(K)]
print("bracket[0] = %s   (must be 0)" % br[0])
BR = dropx(br)  # bracket / x
print("bracket/x at 0 = %s   (should be 2/3)" % BR[0])
sys.stdout.flush()
THt = dropx(TH)  # theta / x   -> 1/3 at 0
# H = 2/bracket - 1/(theta(1-theta))  = (1/x)[ 2/BR - 1/(THt (1-theta)) ]
Ht = [2 * c for c in inv(BR)]
Hb = inv(mul(THt, [one[i] - TH[i] for i in range(K)]))
Hx = [Ht[i] - Hb[i] for i in range(K)]
print("numerator[0] = %s   (must be 0 : the 1/theta poles cancel)" % Hx[0])
sys.stdout.flush()
H_in_x = dropx(Hx)
Hth = comp(H_in_x, revert(TH))  # H as a series in theta
print()
print("   H(theta) = sum h_k theta^k     ARE ALL COEFFICIENTS POSITIVE ?")
neg = 0
for k in range(0, K - 2):
    val = float(Hth[k].numer().str()) / float(Hth[k].denom().str()) if Hth[k] != 0 else 0.0
    if Hth[k] < 0:
        neg += 1
    if k < 14 or k >= K - 6:
        print(
            "     h_%-2d = %-34s = %+.10f %s"
            % (k, Hth[k], val, "   <-- NEGATIVE" if Hth[k] < 0 else "")
        )
print()
print("   negative coefficients among the first %d :  %d" % (K - 2, neg))
if neg == 0:
    rr = [
        float(Hth[k + 1].numer().str())
        / float(Hth[k + 1].denom().str())
        / (float(Hth[k].numer().str()) / float(Hth[k].denom().str()))
        for k in range(K - 8, K - 3)
    ]
    print("   consecutive ratios h_{k+1}/h_k near the top: %s" % ["%.6f" % x for x in rr])
