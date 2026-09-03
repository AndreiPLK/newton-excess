from flint import fmpq, arb, ctx
import sys

ctx.prec = 1200


def toarb(q):
    return arb(q.numer().str()) / arb(q.denom().str())


def make(M):
    N = 2 * M
    b = [toarb(fmpq((2 * j - 1) ** 2)) for j in range(1, M + 1)] * 2

    def qs(z):
        return [v * z / (1 + v * z) for v in b]

    def A(z):
        return sum(qs(z), arb(0))

    def cum(z):
        s2 = arb(0)
        s3 = arb(0)
        for u in qs(z):
            w = u * (1 - u)
            s2 += w
            s3 += w * (1 - 2 * u)
        return s2, s3

    def saddle(x):
        lo, hi = arb(1) / arb(10**14), arb(10**14)
        for _ in range(220):
            mid = (lo * hi).sqrt()
            if A(mid) < x:
                lo = mid
            else:
                hi = mid
        return (lo * hi).sqrt()

    return N, saddle, cum


print("   CLOSED FORM :  N^2 D' = (1-2t)/(t^2(1-t)^2) * [ 1 - k3ratio/rho^3 ]")
print("     n     theta    N^2 D' numeric    closed form      ratio       A3 g2^3 < g3 A2^3 ?")
for M in (20, 40, 80):
    N, saddle, cum = make(M)
    n = N + 1

    def D(x):
        z = saddle(x)
        s2, _ = cum(z)
        th = x / arb(N)
        return 1 / s2 - 1 / (arb(N) * th * (1 - th))

    hh = arb(1) / arb(500)
    for frac in (0.05, 0.15, 0.25, 0.35):
        x = arb(N) * arb(int(frac * 100)) / arb(100)
        if float(x) < 2.5:
            continue
        num = ((D(x + hh) - D(x - hh)) / (2 * hh)) * arb(N) * arb(N)
        z = saddle(x)
        s2, s3 = cum(z)
        th = x / arb(N)
        b2 = arb(N) * th * (1 - th)
        b3 = b2 * (1 - 2 * th)
        rho = s2 / b2
        k3r = s3 / b3
        cf = (1 - 2 * th) / (th * th * (1 - th) * (1 - th)) * (1 - k3r / rho**3)
        A2 = s2 / arb(N)
        A3 = s3 / arb(N)
        g2 = th * (1 - th)
        g3 = g2 * (1 - 2 * th)
        print(
            "   %4d   %.3f   %+12.6f     %+12.6f    %.8f     %s"
            % (
                n,
                float(th),
                float(num),
                float(cf),
                float(num / cf),
                "YES" if A3 * g2**3 < g3 * A2**3 else "NO",
            )
        )
        sys.stdout.flush()
