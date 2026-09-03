from flint import fmpq
from math import comb
import sys

print("   RLC ON THE WHOLE SCALE for the centred-square spectrum :  where exactly does it fail ?")
print("     n     i tested      failures (as m = N-i)")
for M in (5, 10, 20, 40, 60, 80):
    N = 2 * M
    n = N + 1
    e = [fmpq(1)]
    for j in range(1, M + 1):
        c = fmpq((2 * j - 1) ** 2)
        for _ in range(2):
            e = [e[q] + (c * e[q - 1] if q else 0) for q in range(len(e))] + [e[-1] * c]
    p = [e[j] / fmpq(comb(N, j)) for j in range(N + 1)]
    bad = []
    for i in range(0, N - 2):
        if not (p[i + 2] ** 3 * p[i] >= p[i + 1] ** 3 * p[i + 3]):
            bad.append(N - i)
    print(
        "   %4d   i = 0..%-4d    %s"
        % (n, N - 3, ("m = " + ", ".join(str(x) for x in bad)) if bad else "NONE")
    )
    sys.stdout.flush()
