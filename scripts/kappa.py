from flint import fmpq, arb, ctx
from math import comb, factorial
import sys

ctx.prec = 800
print("   THE KAPPA ARGUMENT :  ehat_0 = 1 carries no doubling factor, the others carry 1/2")
print()
print("   step 1 : the doubled reciprocal generating function")
print("      prod (1+z/(2k-1)^2)^2 = cosh^2(pi sqrt z/2) = (1+cosh(pi sqrt z))/2")
print("      => ehat_0 = 1 ,  ehat_m = pi^(2m)/(2 (2m)!)  for m >= 1")
print()
print("   step 2 : the predicted limit ratio, with the kappa counted correctly")
for m in range(3, 10):

    def eh(j):
        return fmpq(1) if j == 0 else fmpq(1, 2 * factorial(2 * j))  # pi^{2j} tracked separately

    A = eh(m - 2) ** 3 * eh(m) * fmpq(factorial(m - 2) ** 3 * factorial(m))
    B = eh(m - 1) ** 3 * eh(m - 3) * fmpq(factorial(m - 1) ** 3 * factorial(m - 3))
    naive = fmpq((2 * m - 3) ** 2, (2 * m - 1) * (2 * m - 5))
    print(
        "     m=%2d   corrected A/B = %-12s   naive (2m-3)^2/((2m-1)(2m-5)) = %-12s   %s"
        % (
            m,
            A / B,
            naive,
            "SAME" if A / B == naive else "HALVED" if A / B == naive / 2 else "other",
        )
    )
print()
print("   step 3 : against the exact object at m = 3")
for M in (20, 40, 80, 160):
    N = 2 * M
    e = [fmpq(1)]
    for j in range(1, M + 1):
        c = fmpq((2 * j - 1) ** 2)
        for _ in range(2):
            e = [e[q] + (c * e[q - 1] if q else 0) for q in range(len(e))] + [e[-1] * c]
    p = [e[j] / fmpq(comb(N, j)) for j in range(N + 1)]
    i = N - 3
    r = (p[i + 2] ** 3 * p[i]) / (p[i + 1] ** 3 * p[i + 3])
    print(
        "     n=%4d   exact ratio at m=3 = %.10f     prediction 9/10 = 0.9000000000"
        % (N + 1, float(arb(r.numer().str()) / arb(r.denom().str())))
    )
    sys.stdout.flush()
