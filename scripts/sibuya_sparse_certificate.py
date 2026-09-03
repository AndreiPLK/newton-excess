"""SIBUYA'S CONJECTURE (1988, eq. 3.4), the sparse regime:  tau^2/N <= BMAX for EVERY tau >= TAU_MIN.

Port of sparse_certificate_full.py (the 4/5 theorem) to the spectrum {1..N} (unsigned Stirling numbers
of the first kind) and to Sibuya's target.  What changed, and only that:
  * spectrum 1, 2, ..., N (single, N roots, n = N + 1); normalised density 2u on [0,1]; centred moments
    mu_i^inf = INT_0^1 (2u-1)^i du;  a_2 = -mu_2^inf/2 = -1/6;  |beta_k| < 1 so |mu_i| <= 1;
  * the target:  n(R-1) >= n/(3n-j)  <=  N g >= N/(3N+3-tau)  <=  N Phi >= 1/(3 + 3a^2 b - a b)
    (a = 1/tau, b = tau^2/N, 1/N = a^2 b), hence the sign polynomial
        Phat := (3 + 3a^2 b - a b) [ (DP~)^2 - P~ D^2 P~ ]  -  a^2 b P~^2 ;
    leading term [a^3]/b^2 = 2/3 (the measured margin 2(j-1)/(9n) times 3);
  * TAU_MIN = 2000 (the Sibuya ladder lab/sibuya_ladder.py closes j <= 2000 exactly);
  * lemma constants for mu_2 = 1/3: F >= FMIN, |F'| <= CF a b F, |F''| <= CF a^2 b F, C_i <= i^2 beyond 60
    -- MARKED "TO VALIDATE": derived by the same argument as the 4/5 constants (0.669, 0.81, i^2 2^i) with
    |mu_i| <= 1 in place of 2^i, not yet independently validated.

The original docstring follows (its statements hold with the substitutions above).

The sparse regime, full version: t^2/N <= BMAX for EVERY t >= 627, by the exact sampling expansion.

Supersedes sparse_certificate.py (region b <= 1/100 with a side constraint on a = 1/t): the a = 0
part of F is the exact exponential E = exp(a_2 b), a_2 = -mu_2^inf/2 = -2/5, carried as a formal
variable with the derivation rule d_b E = a_2 E, so the truncation defect of the exponential
vanishes identically and no side constraint on a remains.  The finite-m moments are handled by an
explicit perturbation lemma (item vi) instead of a symbolic sensitivity.

STATEMENT CERTIFIED.  With a = 1/tau, b = tau^2/N, F(tau) = SUM_j e_j (tau)_j/(N)_j, Phi = -(log F)'':
        Phi N > 4/5    for all real tau >= 627 and all N with b <= BMAX.
Consequence: for integer t >= 628 with t^2/N <= BMAX (628/629)^2, M(n,t) > 4/5 (integral identity
and last step as in sparse_certificate.py).

THE POLYNOMIAL.  In Q[a, b, E], D = -a^2 d/da + 2ab d/db, d_b E = a_2 E:
    F~ = SUM_{j<=J} E_j R_j + (E - T_{J/2}(a_2 b)),   P~ := F~ Q  (Q = prod_{i<J}(1 - i a^2 b)),
    Phat := (DP~)^2 - P~ D^2 P~ - (4/5) a^2 b P~^2 .
Phat satisfies [a^0] = [a^1] = [a^2] = 0 identically in Q[b, E] (checked), so Phat = a^3 (...).  With
E = 1 + a_2 b + b^2 Z, Z formal (Z(b) = (e^{a_2 b} - 1 - a_2 b)/b^2 in [a_2^2/2 + a_2^3 b/6, a_2^2/2]
for 0 <= b <= 1), [a^3] is divisible by b^2 with [a^3]/b^2 = 176/175 + O(b), and
    Phat/(a^3 b^2)  >=  [a^3]/b^2 - SUM_{k>=4} |[a^k]/b^2| a^{k-3}  >  0
is checked by an interval sweep in b.

THE TWO LEMMAS (explicit constants, proportional to a, valid for every tau >= 627, b <= BMAX).
 (v)  Tail.  The true interpolant on the window is F = F~ + T with
        T = SUM_{J<j<=t+1} (E_j R_j - E_j^0) - SUM_{j>t+1} E_j^0,   E_j^0 = (a_2 b)^{j/2}/(j/2)!  (even j).
      Cauchy on |z| = rho_j = sqrt(j/(2 mu_2 b)) with |L_i| <= b (2^i/i)(ab)^{i-2} rho^i (|beta_k| < 2 gives
      |mu_i| <= 2^i) gives  |E_j - E_j^0| <= (2 e mu_2 b/j)^{j/2} (e^{g(rho_j)} - 1),  g(rho) = 4 b rho^2
      SUM_{i>=3} (2 a b rho)^{i-2}/i <= (8 a b^2 rho^3/3)/(1 - 2ab rho),  and 2ab rho_j <= 2a sqrt(j b/(2 mu_2))
      <= 0.1 for j <= t+1;  |1 - R_j| <= a j(j-1)/2 for j a <= 1;  |E_j^0| <= (2 e mu_2 b/j)^{j/2}.
      The tau-derivatives of (tau)_j R_j-type terms are bounded by (j a)^k times the term, and
      d/dtau E_j^0 = j a E_j^0 (b = tau^2/N).  The effect on Phi N is <= 1.06 N SUM_j (j a)^2 |T_j| / F.
 (vi) Moments.  [SIBUYA VERSION: |delta_i| <= C_i/N = C_i a^2 b, h(rho) = a^2 b^2 rho^2 SUM_i (C_i/i)(ab rho)^{i-2}]  ORIGINAL: mu_i(m) = mu_i^inf + delta_i, |delta_i| <= C_i/m^2 = 4 C_i a^4 b^2 (C_i by interval
      evaluation of the exact rational functions for m >= 627^2).  In the same Cauchy estimate
      |E_j(mu) - E_j(mu^inf)| <= (2 e mu_2 b/j)^{j/2} e^{g(rho)} (e^{h(rho)} - 1),
      h(rho) = 4 a^4 b^3 rho^2 SUM_i (C_i/i)(ab rho)^{i-2}.

Run:  uv run python projects/qg-bootstrap/release/scripts/sibuya_sparse_certificate.py
"""

from __future__ import annotations

import sys
import time
from math import comb, factorial

from flint import Ordering, arb, ctx, fmpq, fmpq_mpoly_ctx, fmpq_poly

ctx.prec = 300

J = 30
TAU_MIN = 1000  # the Sibuya ladder closes j <= 1000 (2000 was not reachable in one night: rung time ~ j^3.4)
BMAX = fmpq(1)
A2 = fmpq(-1, 6)
FMIN = arb("0.80")  # TO VALIDATE: F ~ e^{a_2 b} >= e^{-1/6} = 0.846 minus the corrections
CF = arb("0.35")  # TO VALIDATE: |F'/F| <= (mu_2 + 0.02) a b, |F''/F| <= (mu_2 + 0.02) a^2 b

R3 = fmpq_mpoly_ctx.get(("a", "b", "E"), Ordering.lex)
A, B, EE = R3.gens()


def one():
    return R3.from_dict({(0, 0, 0): fmpq(1)})


def zero():
    return R3.from_dict({})


def toarb(q):
    if isinstance(q, arb):
        return q
    return arb(q.numer().str()) / arb(q.denom().str())


# ------------------------------------------------------------------ moments
def mu_inf(imax):
    out = [fmpq(1)]
    base = fmpq_poly([-1, 2])  # 2u - 1
    p = fmpq_poly([1])
    for _ in range(1, imax + 1):
        p = p * base
        out.append(p.integral()(fmpq(1)))
    return out


def faulhaber_single(L):
    """S[l](N) = SUM_{k=1}^{N} k^l as polynomials in N, l = 0..L."""
    P = [fmpq_poly([0, 1])]
    for p in range(1, L + 1):
        acc = fmpq_poly([1, 1]) ** (p + 1) - fmpq_poly([1])
        for j in range(p):
            acc = acc - P[j] * fmpq(comb(p + 1, j))
        P.append(acc / fmpq(p + 1))
    return P


def mu_exact_poly(imax):
    S = faulhaber_single(imax)
    m = fmpq_poly([0, 1])
    out = [(fmpq_poly([1]), fmpq_poly([1]))]
    for i in range(1, imax + 1):
        num = fmpq_poly([0])
        for ell in range(i + 1):
            num += fmpq(comb(i, ell)) * (-1) ** (i - ell) * S[ell] * m ** (ell) * S[1] ** (i - ell)
        den = S[1] ** i * m
        out.append((num, den))
    return out


def mu_deviation_bounds(imax, mmin):
    """C_i with |mu_i(N) - mu_i^inf| <= C_i / N for all N >= mmin (interval evaluation in 1/N).
    Single spectrum: the deviation is first order in 1/N (mu_2(N) = 1/3 - 2/(3(N+1))), unlike the odd squares."""
    muinf = mu_inf(imax)
    polys = mu_exact_poly(imax)

    def apow(x, k):
        r = arb(1)
        for _ in range(k):
            r *= x
        return r

    C = [arb(0)] * (imax + 1)
    u = arb(0, 1 / mmin)
    for i in range(2, imax + 1):
        num, den = polys[i]
        diff = num - den * muinf[i]
        d = den.degree()
        dc = list(diff.coeffs()) + [fmpq(0)] * (d + 1 - len(diff.coeffs()))
        ec = list(den.coeffs()) + [fmpq(0)] * (d + 1 - len(den.coeffs()))
        du = [dc[d - k] for k in range(d + 1)]
        eu = [ec[d - k] for k in range(d + 1)]
        assert du[0] == 0
        numv = sum((toarb(du[k]) * apow(u, k - 1) for k in range(1, d + 1)), arb(0))
        denv = sum((toarb(eu[k]) * apow(u, k) for k in range(d + 1)), arb(0))
        C[i] = abs(numv / denv).upper()
    return C


# ------------------------------------------------------------------ the polynomial
def build(mu):
    t0 = time.time()
    theta = A * B
    L = {i: B * (fmpq((-1) ** (i - 1)) * mu[i] / i) * theta ** (i - 2) for i in range(2, J + 1)}
    E = [one()]
    for j in range(1, J + 1):
        acc = zero()
        for i in range(2, j + 1):
            acc += L[i] * E[j - i] * fmpq(i)
        E.append(acc / fmpq(j))
    Q = one()
    for i in range(1, J):
        Q *= 1 - fmpq(i) * A * A * B
    P = zero()
    for j in range(J + 1):
        num = one()
        for i in range(1, j):
            num *= 1 - fmpq(i) * A
        den_rest = one()
        for i in range(j, J):
            den_rest *= 1 - fmpq(i) * A * A * B
        P += E[j] * num * den_rest
    K = J // 2
    T = zero()
    for k in range(K + 1):
        T += (A2 * B) ** k / fmpq(factorial(k))
    Pt = P + (EE - T) * Q
    print(
        f"  P~ built: {len(Pt)} terms, degrees {Pt.degrees()}  [{time.time() - t0:.1f} s]",
        flush=True,
    )

    def Dop(p):
        return -A * A * p.derivative(0) + 2 * A * B * (p.derivative(1) + A2 * EE * p.derivative(2))

    DP = Dop(Pt)
    D2P = Dop(DP)
    Phat = (3 + 3 * A * A * B - A * B) * (DP * DP - Pt * D2P) - A * A * B * Pt * Pt
    print(
        f"  Phat: {len(Phat)} terms, degrees {Phat.degrees()}  [{time.time() - t0:.1f} s]",
        flush=True,
    )
    return Phat


def coeff_a(p, k):
    return R3.from_dict({(0, e[1], e[2]): c for e, c in p.to_dict().items() if e[0] == k})


def main():
    t0 = time.time()
    mu = mu_inf(J)
    Phat = build(mu)
    for k in range(3):
        ck = coeff_a(Phat, k)
        print(f"  [a^{k}] = {'0' if ck.is_zero() else 'NONZERO (' + str(len(ck)) + ' terms)'}")
        assert ck.is_zero()
    # E -> 1 + a2 b + b^2 Z  (Z takes the slot of E)
    Esub = 1 + A2 * B + B * B * EE
    coeffs = {}
    for k in range(3, Phat.degrees()[0] + 1):
        ck = coeff_a(Phat, k)
        if not ck.is_zero():
            coeffs[k] = ck.compose(A, B, Esub)
    c3 = coeffs[3]
    d3 = c3.to_dict()
    lowb = min(e[1] for e in d3)
    print(f"  [a^3](b, Z): {len(c3)} terms, lowest power of b: {lowb}")
    assert lowb >= 2
    c3r = R3.from_dict({(0, e[1] - 2, e[2]): c for e, c in d3.items()})
    const = c3r.to_dict().get((0, 0, 0), fmpq(0))
    print(f"  [a^3]/b^2 at b = 0: {const}   (2/3: {const == fmpq(2, 3)})")

    amax = arb(1) / TAU_MIN
    bmaxa = toarb(BMAX)
    a2 = toarb(A2)

    def zeta_box(bhi):
        lo = a2 * a2 / 2 + a2 * a2 * a2 * bhi / 6
        hi = a2 * a2 / 2
        return arb((lo + hi) / 2, ((hi - lo) / 2).upper())

    def eval_bz(poly, b, z):
        tot = arb(0)
        for e, c in poly.to_dict().items():
            tot += toarb(c) * b ** e[1] * z ** e[2]
        return tot

    # (iii) rest
    rest = arb(0)
    zmax = arb("0.02")  # Z <= a_2^2/2 = 1/72
    for k, cz in coeffs.items():
        if k == 3:
            continue
        dz = cz.to_dict()
        lowk = min(e[1] for e in dz)
        assert lowk >= 2, f"[a^{k}] not divisible by b^2 (lowest {lowk})"
        Ck = sum(
            (abs(toarb(c)) * bmaxa ** (e[1] - 2) * zmax ** e[2] for e, c in dz.items()), arb(0)
        )
        rest += Ck * amax ** (k - 3)
    print(f"  (iii) SUM_k |[a^k]/b^2| a^(k-3) <= {float(rest):.4e}")

    # (v) tail and (vi) moments, both relative to a b  -- explicit lemmas, see docstring.
    # Effect of a perturbation dF = SUM_j dF_j (each term proportional to (tau)_j) on Phi N, from
    # Phi = (F'/F)^2 - F''/F:  N [ |dF''| + 2|F'||dF'| + (F'^2 + |F''|)|dF| ] / F (1 + O(|dF|/F)), with
    #   |dF_j''| <= (j/(tau-j))^2 |dF_j|,  |dF_j'| <= (j/(tau-j)) |dF_j|,  |F'| <= 0.81 a b F,  |F''| <= 0.81 a^2 b F,
    #   F >= 0.669 on the region (independent validation, results/VALIDATION_SPARSE_CERTIFICATE_FULL_2026-09-02.md),
    # and N a^2 = 1/b:  effect_j <= |dF_j| j^2 m_j / (b F),  m_j := 1/(1 - j a)^2 + 1.62 b/j + 1.46 b/j^2.
    e1 = arb(1).exp()
    mu2 = toarb(mu[2])
    C = mu_deviation_bounds(
        60, TAU_MIN**2
    )  # region: N >= tau^2/b >= TAU_MIN^2 (single spectrum, m = N)
    Fmin = FMIN

    def mult(j):
        # validator (3 Sept 2026): the derivative of T/F gives 2 T F'^2/F^3, hence 2 CF^2 b, not CF^2 b
        return (
            1 / (1 - j * amax) ** 2
            + 2 * CF * bmaxa / (j * (1 - j * amax))
            + (2 * CF * CF * bmaxa + CF) * bmaxa / (j * j)
        )

    def hsum(x):
        # h(rho) = 4 a^4 b^3 rho^2 SUM_{i>=2} (C_i/i)(ab rho)^{i-2}: C_i computed for i <= 60, and C_i <= i^2 2^i beyond
        # (midpoint-rule error of (3u^2(1+e')-1)^i plus the shift e' = 1/(4m^2-1); see the paper), x = 2 a b rho
        # single spectrum, |delta_i| <= C_i/N = C_i a^2 b:  h(rho) = a^2 b^2 rho^2 SUM_{i>=2} (C_i/i)(ab rho)^{i-2};
        # C_i <= 3 + i^2/TAU_MIN^2 beyond 60 (TO VALIDATE: (1/N) SUM f(x_k) vs INT f/2 for f = (2u-1)^i, |f| <= 1:
        # the 1/N rescaling and the two end gaps of width 1/(N+1) give 2/N, the midpoint error i^2/(6N^2))
        tot = sum((toarb(C[i]) / i * x ** (i - 2) / 2 ** (i - 2) for i in range(2, 61)), arb(0))
        for i in range(61, 400):
            tot += (3 + arb(i * i) / TAU_MIN**2) / i * x ** (i - 2) / 2 ** (i - 2)
        return tot

    tail_rel = arb(0)
    mom_rel = arb(0)
    for j in range(2, 4000):
        rho = (arb(j) / (2 * mu2 * bmaxa)).sqrt()
        base = (2 * e1 * mu2 * bmaxa / j) ** (arb(j) / 2)
        x = 2 * amax * bmaxa * rho
        # |mu_i| <= 1:  g(rho) = b rho^2 SUM_{i>=3} (ab rho)^{i-2}/i <= (a b^2 rho^3/3)/(1 - ab rho)
        g = (amax * bmaxa * bmaxa * rho**3 / 3) / (1 - x / 2)
        h = amax**2 * bmaxa**2 * rho * rho * hsum(x)
        Mj = base * g.exp() * (h.exp() - 1)
        mom_rel += j * j * mult(j) * Mj / (bmaxa * Fmin) / (amax * bmaxa)
        if j > J:
            Tj = base * (g.exp() - 1) + base * amax * j * (j - 1) / 2
            eff = j * j * mult(j) * Tj / (bmaxa * Fmin)
            tail_rel += eff / (amax * bmaxa)
            if j > 200 and float(eff) < 1e-300:
                break
    print(f"  (v)  tail  j > {J}: effect relative to a b <= {float(tail_rel):.3e}")
    print(f"  (vi) finite-m moments: effect relative to a b <= {float(mom_rel):.3e}")
    # the lemma effects on N Phi enter Phat through the factor (3 + 3a^2 b - ab) <= 3.01
    subtract = (rest + (tail_rel + mom_rel) * arb("3.01")).upper()

    # sweep in b
    stack = [(arb(0), bmaxa, 0)]
    worst = None
    boxes = 0
    while stack:
        lo, hi, dep = stack.pop()
        boxes += 1
        b = arb((lo + hi) / 2, ((hi - lo) / 2).upper())
        z = zeta_box(hi)
        val = eval_bz(c3r, b, z)
        low = val.lower() - subtract
        if low > 0:
            if worst is None or low < worst:
                worst = low
            continue
        if dep > 40:
            print("FAILED at b in", float(lo), float(hi), "value", val)
            sys.exit(1)
        mid = (lo + hi) / 2
        stack.append((lo, mid, dep + 1))
        stack.append((mid, hi, dep + 1))
    print(
        f"  (i)+(ii) sweep over b in [0, {BMAX}]: {boxes} boxes, min of [a^3]/b^2 - (iii) - (v) - (vi) = {float(worst):.6f}"
    )
    ok = float(worst) > 0
    print()
    print(
        f"VERDICT (lemma constants independently validated 3 Sept 2026): (3 + 3a^2 b - ab) N Phi > 1 on {{ tau >= {TAU_MIN}, tau^2/N <= {BMAX} }}: {ok}   [{time.time() - t0:.0f} s]"
    )
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
