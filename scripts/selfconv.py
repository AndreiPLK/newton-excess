from flint import fmpq
from math import comb
import sys


def esyms(b):
    e = [fmpq(1)]
    for v in b:
        e = [e[q] + (v * e[q - 1] if q else 0) for q in range(len(e))] + [e[-1] * v]
    return e


def rlcfail(b):
    N = len(b)
    e = esyms(b)
    p = [e[j] / fmpq(comb(N, j)) for j in range(N + 1)]
    return [
        N - i for i in range(0, N - 2) if not (p[i + 2] ** 3 * p[i] >= p[i + 1] ** 3 * p[i + 3])
    ]


print("   THE SPECTRUM IS A SELF-CONVOLUTION :  the doubled polynomial is P(x)^2")
print("   does RLC hold for the SINGLE odd-square spectrum, and survive the squaring ?")
print("     M    single {(2k-1)^2}          doubled = self-convolution")
for M in (4, 6, 8, 12, 20, 30, 40):
    single = [fmpq((2 * j - 1) ** 2) for j in range(1, M + 1)]
    doubled = single * 2
    fs = rlcfail(single)
    fd = rlcfail(doubled)
    print(
        "   %4d   fails at m = %-22s  fails at m = %s"
        % (
            M,
            (", ".join(map(str, fs)) if fs else "NONE"),
            (", ".join(map(str, fd)) if fd else "NONE"),
        )
    )
    sys.stdout.flush()
