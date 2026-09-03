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

    def cums(z):
        q = qs(z)
        s2 = arb(0)
        s3 = arb(0)
        for u in q:
            s2 += u * (1 - u)
            s3 += u * (1 - u) * (1 - 2 * u)
        return s2, s3

    def saddle(x):
        lo, hi = arb(1) / arb(10**14), arb(10**14)
        for _ in range(200):
            mid = (lo * hi).sqrt()
            if A(mid) < x:
                lo = mid
            else:
                hi = mid
        return (lo * hi).sqrt()

    return N, saddle, cums


print("   IDENTITY CHECK :  D'(x) = kappa3_bin/sigma_bin^6 - kappa3/sigma^6 ?")
print("     n     x      D' numeric       D' from cumulants     ratio")
for M in (20, 40):
    N, saddle, cums = make(M)
    n = N + 1

    def D(x):
        z = saddle(x)
        s2, _ = cums(z)
        sb = arb(x) * (arb(N) - arb(x)) / arb(N)
        return 1 / s2 - 1 / sb

    hh = arb(1) / arb(1000)
    for x in (5, M // 2, M - 2):
        num = (D(arb(x) + hh) - D(arb(x) - hh)) / (2 * hh)
        z = saddle(arb(x))
        s2, s3 = cums(z)
        th = arb(x) / arb(N)
        sb = arb(N) * th * (1 - th)
        k3b = arb(N) * th * (1 - th) * (1 - 2 * th)
        cum = k3b / sb**3 - s3 / s2**3
        print(
            "     %4d %5.1f   %+14.8e   %+14.8e   %.8f"
            % (n, x, float(num), float(cum), float(num / cum))
        )
        sys.stdout.flush()
print()
print("   THE ELEMENTARY FORM :  kappa3/sigma^6 vs its value at uniform q  (same sum)")
print("     n     x      spectrum k3/s^6      uniform k3/s^6      spectrum smaller ?")
for M in (10, 20, 40):
    N, saddle, cums = make(M)
    n = N + 1
    for x in (2, M // 2, M - 2):
        z = saddle(arb(x))
        s2, s3 = cums(z)
        th = arb(x) / arb(N)
        sb = arb(N) * th * (1 - th)
        k3b = arb(N) * th * (1 - th) * (1 - 2 * th)
        a = s3 / s2**3
        c = k3b / sb**3
        print(
            "     %4d %5.1f   %+14.8e   %+14.8e    %s"
            % (n, x, float(a), float(c), "YES" if a < c else "NO")
        )
        sys.stdout.flush()
