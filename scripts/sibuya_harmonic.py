"""SIBUYA'S CONJECTURE (1988, eq. 3.4), the few-missing-roots corner:  j' = N - j in [1, JMAX-1], every N >= 1001 + j'.

For the spectrum {1..N} write E_k := e_k(1, 1/2, ..., 1/N).  Since e_j(1..N) = N! E_{N-j}, the normalised ratio at
index j = N - j' is

    p_j^2/(p_{j-1} p_{j+1}) = [E_{j'}^2 / (E_{j'+1} E_{j'-1})] * [j'(N-j')/((j'+1)(N-j'+1))],

and Sibuya's target 1 + 1/(3n-j) = 1 + 1/(2N+3+j').  So the statement is  F := E_{j'}^2 - T E_{j'+1} E_{j'-1} >= 0 with

    T(N, j') = (1 + 1/(2N+3+j')) (j'+1)(N-j'+1) / (j'(N-j')),      T -> (j'+1)/j'  as N -> infinity.

The E_k are polynomials in the power sums P_r = H_N^{(r)} (Newton's identities), P_1 = H_N =: H, and for r >= 2
P_r in [zeta(r) - 1/((r-1) N^{r-1}), zeta(r)) (integral bounds on the tail).  Hence, for a band N in [N_a, N_b]:
H in [H(N_a), H(N_b)], P_r in an interval, T in [T(N_b), T(N_a)], and F is a polynomial in H with interval
coefficients whose positivity on the H-interval is checked by evaluation with subdivision (the dependence of H, P_r, T
on the same N is dropped -- conservative).  Bands double from N = 1001 + j' to N_TAIL = 2^48.

The tail N >= N_TAIL: with G_0 := E_{j'}^2 - ((j'+1)/j') E_{j'+1} E_{j'-1} (exact leading cancellation: G_0 has degree
2j'-2, leading coefficient a positive multiple of P_2) and the deficit T - (j'+1)/j' <= K/N, K = 1.7 (j'+1)/j' (N >= 10 j'),
1/N < e^{gamma - H} (H_N > log N + gamma), and E_{j'+1} E_{j'-1} <= H^{2j'}/((j'+1)!(j'-1)!):

    F >= G_0(H) - K e^{gamma} e^{-H} H^{2j'} / ((j'+1)!(j'-1)!),

and dividing by H^{2j'-2}: G_0/H^{2j'-2} >= c_lead - SUM_{i<2j'-2} |g_i| H_T^{i-(2j'-2)} (each term decreasing in H >= H_T),
while e^{-H} H^2 is decreasing for H >= 2; one inequality at H = H_T closes the tail.

Run:  uv run python projects/qg-bootstrap/release/scripts/sibuya_harmonic.py [JMAX]
"""

from __future__ import annotations

import sys
import time
from math import factorial

from flint import arb, ctx, fmpq

ctx.prec = 400
JMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 100
J_MIN_INDEX = 1001  # the ladder closes j <= 1000
N_TAIL_LOG2 = 48
EULER = arb.const_euler()


def harmonic(N):
    """H_N = psi(N+1) + gamma, certified."""
    return arb(N + 1).digamma() + EULER


def power_sum_interval(r, Na, Nb):
    """P_r = H_N^(r) for N in [Na, Nb]: zeta(r) - tail, tail in [1/((r-1)(Nb+1)^(r-1)), 1/((r-1) Na^(r-1))]."""
    z = arb(r).zeta()
    lo = z - arb(1) / ((r - 1) * arb(Na) ** (r - 1))
    hi = z - arb(1) / ((r - 1) * arb(Nb + 1) ** (r - 1))
    return arb((lo + hi) / 2, ((hi - lo) / 2).upper() + 1e-90)


# ---- polynomials in H as lists of arb coefficients (low to high)
def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else arb(0)) + (b[i] if i < len(b) else arb(0)) for i in range(n)]


def pscale(a, s):
    return [c * s for c in a]


def pmul(a, b):
    out = [arb(0)] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, w in enumerate(b):
            out[i + j] += u * w
    return out


def peval(a, H):
    tot = arb(0)
    for c in reversed(a):
        tot = tot * H + c
    return tot


def elementary_polys(kmax, P):
    """E_k as polynomials in H via Newton's identities: k E_k = SUM_{r=1}^k (-1)^{r-1} P_r E_{k-r}, P_1 = H."""
    E = [[arb(1)]]
    for k in range(1, kmax + 1):
        acc = [arb(0)]
        for r in range(1, k + 1):
            term = pmul([arb(0), arb(1)], E[k - 1]) if r == 1 else pscale(E[k - r], P[r])
            acc = padd(acc, pscale(term, (-1) ** (r - 1)))
        E.append(pscale(acc, arb(1) / k))
    return E


def T_of(N, jp):
    N = arb(N)
    return (1 + 1 / (2 * N + 3 + jp)) * (jp + 1) * (N - jp + 1) / (jp * (N - jp))


def positive_on_fn(f, Hlo, Hhi, depth=0):
    """Is f(H) > 0 on [Hlo, Hhi]?  Interval evaluation with bisection; f takes an arb and returns an arb."""
    H = arb((Hlo + Hhi) / 2, ((Hhi - Hlo) / 2).upper())
    val = f(H)
    if val.lower() > 0:
        return True
    if depth > 60:
        return False
    mid = (Hlo + Hhi) / 2
    return positive_on_fn(f, Hlo, mid, depth + 1) and positive_on_fn(f, mid, Hhi, depth + 1)


def band_ok(jp, Na, Nb, depth=0):
    """Is F > 0 for every integer N in [Na, Nb]?  Interval evaluation; on failure the band is split
    geometrically (the loss comes from treating H, P_r and T as independent over the band)."""
    P = [None, None] + [power_sum_interval(r, Na, Nb) for r in range(2, jp + 2)]
    E = elementary_polys(jp + 1, P)
    # Each E_k is evaluated SEPARATELY at H and the products are formed numerically: multiplying the
    # polynomials first and evaluating afterwards loses the correlation between the factors and inflates
    # the enclosure by two orders of magnitude (measured at j' = 14: +/- 0.47 against a value of 8).
    Tmax = T_of(
        Na, jp
    ).upper()  # T decreases in N; E_{j'+1} E_{j'-1} > 0, so T_max is the conservative choice

    def F(H):
        return peval(E[jp], H) ** 2 - arb(Tmax) * peval(E[jp + 1], H) * peval(E[jp - 1], H)

    if positive_on_fn(F, arb(harmonic(Na).lower()), arb(harmonic(Nb).upper())):
        return True, 1
    if Nb - Na <= 1 or depth > 60:
        return False, 1
    Nm = int((Na * Nb) ** 0.5)
    Nm = min(max(Nm, Na + 1), Nb - 1)
    ok1, c1 = band_ok(jp, Na, Nm, depth + 1)
    if not ok1:
        return False, c1
    ok2, c2 = band_ok(jp, Nm, Nb, depth + 1)
    return ok2, c1 + c2


N_EXACT = 1000000  # below this the 1/N corrections are of the size of the margin: exact recursion, not bands


def exact_pass(jmax, N_hi=N_EXACT):
    """One incremental pass over N: E_k(N) = E_k(N-1) + E_{k-1}(N-1)/N for the spectrum {1/k}, all k <= jmax+1
    at once, in certified arb.  Checks F = E_j'^2 - T E_{j'+1} E_{j'-1} > 0 for every j' < jmax with
    N >= J_MIN_INDEX + j'.  Returns (ok, worst relative margin, where)."""
    E = [arb(1)] + [arb(0)] * (jmax + 1)
    worst = None
    Tq = [None] + [fmpq(jp + 1, jp) for jp in range(1, jmax + 1)]
    for N in range(1, N_hi + 1):
        inv = arb(1) / N
        for k in range(min(jmax + 1, N), 0, -1):
            E[k] = E[k] + inv * E[k - 1]
        if N < J_MIN_INDEX + 1:
            continue
        for jp in range(1, jmax):
            if N < J_MIN_INDEX + jp:
                continue
            T = (1 + arb(1) / (2 * N + 3 + jp)) * arb(Tq[jp]) * (arb(N - jp + 1) / (N - jp))
            F = E[jp] ** 2 - T * E[jp + 1] * E[jp - 1]
            if not (F.lower() > 0):
                return False, None, f"exact pass FAILED at N = {N}, j' = {jp}: F = {F.str(6)}"
            rel = float(F / (E[jp] ** 2))
            if worst is None or rel < worst[0]:
                worst = (rel, N, jp)
    return True, worst, ""


def tail_ok(jp, K):
    """The tail N >= 2^K: F >= G_0(H) - K_delta e^{-H} H^{2j'}/((j'+1)!(j'-1)!), and
    G_0(H)/H^{2j'-2} >= c_lead - SUM_{i<deg} |g_i| H^{i-deg} at H = H(2^K) (each term decreasing in H)."""
    N_tail = 2**K
    P = [None, None] + [power_sum_interval(r, N_tail, 10**400) for r in range(2, jp + 2)]
    E = elementary_polys(jp + 1, P)
    G0 = padd(pmul(E[jp], E[jp]), pscale(pmul(E[jp + 1], E[jp - 1]), -arb(fmpq(jp + 1, jp))))
    assert abs(G0[2 * jp]).upper() < 1e-40 and abs(G0[2 * jp - 1]).upper() < 1e-40, (
        "leading terms did not cancel"
    )
    deg = 2 * jp - 2
    c_lead = G0[deg].lower() if jp >= 2 else G0[0].lower()
    if not (c_lead > 0):
        return False, "tail: leading coefficient of G_0 not positive"
    HT = arb(harmonic(N_tail).lower())
    c_min = arb(c_lead)
    for i in range(deg):
        c_min -= abs(G0[i]).upper() / HT ** (deg - i)
    Kd = arb("1.7") * (jp + 1) / jp * EULER.exp()
    q = arb(1) / (factorial(jp + 1) * factorial(jp - 1))
    tail_term = Kd * (-HT).exp() * HT * HT * q
    # (E_{j'+1} E_{j'-1} <= H^{2j'}/((j'+1)!(j'-1)!) since E_k <= H^k/k! for positive roots)
    assert N_tail >= 10 * jp
    return (
        c_min - tail_term
    ).lower() > 0, f"c_min = {float(c_min):.3e}, e^-H term = {float(tail_term):.3e}"


def certify_jprime(jp, N_start):
    """One wide band from N_start (>= N_EXACT, where the 1/N corrections are far below the margin) up to
    2^K, with K raised until the analytic tail closes."""
    t0 = time.time()
    K = 48
    while K <= 3072:
        ok, msg = tail_ok(jp, K)
        if ok:
            break
        K *= 2
    else:
        return False, "tail does not close even at 2^3072"
    ok, c = band_ok(jp, N_start, 2**K)
    if not ok:
        return False, f"band N in [{N_start}, 2^{K}] FAILED after {c} splits"
    return True, f"{c} band(s) {N_start} -> 2^{K}; tail at 2^{K}: {msg}  [{time.time() - t0:.1f} s]"


def main():
    t0 = time.time()
    print(f"Sibuya (3.4), the corner j' = N - j in [1, {JMAX - 1}], every N >= {J_MIN_INDEX} + j':")
    ok, worst, msg = exact_pass(JMAX)
    if not ok:
        print(f"VERDICT: {msg}")
        sys.exit(1)
    print(
        f"  exact recursion N <= {N_EXACT}: all j' ok; worst relative margin {worst[0]:.4e} at N = {worst[1]}, j' = {worst[2]}"
        f"  [{time.time() - t0:.0f} s]",
        flush=True,
    )
    for jp in range(1, JMAX):
        ok, msg = certify_jprime(jp, N_EXACT)
        if jp <= 3 or jp % 20 == 0 or not ok:
            print(f"  j' = {jp:3d}: {'ok' if ok else 'FAIL'}  {msg}", flush=True)
        if not ok:
            print(f"VERDICT: FAIL at j' = {jp}: {msg}")
            sys.exit(1)
    print()
    print(
        f"VERDICT: Sibuya (3.4) holds for j' = N - j in [1, {JMAX - 1}] and every N >= {J_MIN_INDEX} + j'"
        f" (exact recursion to {N_EXACT}, bands to 2^{N_TAIL_LOG2}, analytic tail): True  [{time.time() - t0:.0f} s]"
    )


if __name__ == "__main__":
    main()
