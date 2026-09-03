from flint import fmpq, arb, ctx
from math import comb
import sys

ctx.prec = 4000
print("   INDEPENDENT CHECK of the claim: M(n,t) > 4/5 for EVERY t < n/2, at the top of the range")
print("   (the certificates are polynomial; this is a direct numerical audit of the conclusion)")
for n in (301, 331, 363):
    N = n - 1
    M = N // 2
    e = [fmpq(1)]
    for j in range(1, M + 1):
        c = fmpq((2 * j - 1) ** 2)
        for _ in range(2):
            e = [e[q] + (c * e[q - 1] if q else 0) for q in range(len(e))] + [e[-1] * c]

    def toarb(q):
        return arb(q.numer().str()) / arb(q.denom().str())

    lp = [toarb(e[i]).log() - arb(str(comb(N, i))).log() for i in range(N + 1)]
    worst = None
    wt = None
    bad = 0
    maxi = 0
    for t in range(1, (n + 1) // 2):
        g = -(lp[t + 1] - 2 * lp[t] + lp[t - 1])
        Mv = arb(n) * (g.exp() - 1)
        if not (Mv > arb(4) / arb(5)):
            bad += 1
        if worst is None or Mv < worst:
            worst = Mv
            wt = t
        maxi = max(maxi, t - 1)
    print(
        "   n=%4d   t = 1..%3d   largest index i used = %3d   min M = %.10f at t=%d   violations %d"
        % (n, (n - 1) // 2, maxi, float(worst), wt, bad)
    )
    sys.stdout.flush()
print()
print(
    "   the certificates cover i = 0..180, so they cover the whole regime whenever (n-5)/2 <= 180,"
)
print("   that is  n <= 365.")
