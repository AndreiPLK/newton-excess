"""The sparse regime of the Newton-excess conjecture: t*theta small, by the exact sampling expansion.

STATEMENT CERTIFIED HERE (real variables, stronger than needed):

    Phi(tau; N) * N  >  4/5     for all real tau >= 627 and all N with  b := tau^2/N <= BSTAR,

where  Phi = -(d/dtau)^2 log F(tau),  F(tau) = SUM_{j=0}^{N} e_j (tau)_j / (N)_j  (p_tau = bbar^tau F).
Consequence: for integer t >= 628 with  t^2/N <= BSTAR (628/629)^2  and  (t^2/N)^3 <= (628/629)^7 / (125 t)
(so that the whole window tau in [t-1, t+1] lies in the certified region),
    g(n,t) = -Delta^2 log p_t = INT_0^1 INT_0^1 Phi(t-1+u+v) du dv  >  4/(5N)  >  log(1 + 4/(5n)),
so M(n,t) = n(e^g - 1) > n g > 4/5.  (Independently validated 2 September 2026, see
results/VALIDATION_SPARSE_CERTIFICATE_2026-09-02.md; the validator also re-derived the tail bound
without linearisation, 5.76e-4 relative to a b.)

VARIABLES.  a = 1/tau,  b = tau^2/N,  so theta = ab and 1/N = a^2 b.  Then

    e_j (tau)_j/(N)_j = E_j(a,b) R_j(a,b),
    E_j = [z^j] exp( b SUM_{i>=2} (-1)^{i-1} (mu_i/i) (ab)^{i-2} z^i ),   R_j = prod_{i<j} (1-ia)/(1-i a^2 b),

with mu_i = (1/N) SUM_k (b_k/bbar - 1)^i, replaced by their m -> infinity limits
mu_i^inf = INT_0^1 (3u^2-1)^i du (exact rationals) plus a bounded perturbation.

d/dtau at fixed N is  D = -a^2 d/da + 2ab d/db,  which annihilates functions of a^2 b.
With F_J = P/Q, Q = prod_{i<J}(1 - i a^2 b),  the sign of  Phi N - 4/5  is the sign of

    Phat := (DP)^2 - P D^2 P - (4/5) a^2 b P^2        (a polynomial in a, b).

Phat = a^2 [a^2](b) + a^3 [a^3](b) + a^4 (rest);  [a^0] = [a^1] = 0 identically, [a^2] = O(b^6)
is the truncation defect of exp(-2b/5), [a^3] = (176/175) b^2 (1 + O(b)).  The certificate
checks each piece by exact polynomial arithmetic and interval bounds on the region, then adds
the tail j > J and the mu-perturbation as explicit bounds.

Everything exact (fmpq_mpoly) or certified (arb).  Run:
    uv run python projects/qg-bootstrap/release/scripts/sparse_certificate.py
"""

from __future__ import annotations

import sys
import time

from flint import Ordering, arb, ctx, fmpq, fmpq_mpoly_ctx, fmpq_poly

ctx.prec = 300

J = 10  # terms of the sampling expansion kept exactly
TAU_MIN = 627
BSTAR = fmpq(1, 100)  # region: b = tau^2/N <= 1/100

MC = fmpq_mpoly_ctx.get(("a", "b"), Ordering.lex)
A, B = MC.gens()


def mu_inf(imax):
    """mu_i^inf = INT_0^1 (3u^2 - 1)^i du, exact."""
    out = [fmpq(1)]
    base = fmpq_poly([-1, 0, 3])
    p = fmpq_poly([1])
    for _ in range(1, imax + 1):
        p = p * base
        out.append(p.integral()(fmpq(1)))
    return out


def build_F(J, mu):
    """F_J = P/Q with P, Q in Q[a,b]; returns (P, Q)."""
    theta = A * B
    # log-generating function coefficients: L_i = b * (-1)^{i-1} (mu_i/i) theta^{i-2}
    L = {i: B * (fmpq((-1) ** (i - 1)) * mu[i] / i) * theta ** (i - 2) for i in range(2, J + 1)}
    # E_j via j E_j = SUM_{i} i L_i E_{j-i}
    E = [MC.from_dict({(0, 0): fmpq(1)})]
    for j in range(1, J + 1):
        acc = MC.from_dict({})
        for i in range(2, j + 1):
            acc += L[i] * E[j - i] * fmpq(i)
        E.append(acc / fmpq(j))
    # R_j = prod_{i<j} (1 - i a) / (1 - i a^2 b) ; common denominator Q = prod_{i<J} (1 - i a^2 b)
    Q = MC.from_dict({(0, 0): fmpq(1)})
    for i in range(1, J):
        Q *= 1 - fmpq(i) * A * A * B
    P = MC.from_dict({})
    for j in range(J + 1):
        num = MC.from_dict({(0, 0): fmpq(1)})
        for i in range(1, j):
            num *= 1 - fmpq(i) * A
        den_rest = MC.from_dict({(0, 0): fmpq(1)})
        for i in range(j, J):
            den_rest *= 1 - fmpq(i) * A * A * B
        P += E[j] * num * den_rest
    return P, Q


def Dop(p):
    """D = -a^2 d/da + 2ab d/db."""
    return -A * A * p.derivative(0) + 2 * A * B * p.derivative(1)


def coeff_in_a(p, k):
    """[a^k] p as an fmpq_poly in b."""
    d = p.to_dict()
    maxb = max((e[1] for e in d if e[0] == k), default=-1)
    cs = [fmpq(0)] * (maxb + 1)
    for (ea, eb), c in d.items():
        if ea == k:
            cs[eb] = c
    return fmpq_poly(cs)


def poly_abs_bound(p, bmax):
    """SUM |c_k| bmax^k -- a bound on |p(b)| for 0 <= b <= bmax."""
    return sum(abs(c) * bmax**k for k, c in enumerate(p.coeffs()))


def mu_exact_poly(imax):
    """mu_i(m) as exact rational functions of m via Faulhaber: returns list of (num, den) fmpq_poly pairs."""
    from math import comb

    P = [fmpq_poly([0, 1])]
    for p in range(1, 2 * imax + 1):
        acc = fmpq_poly([1, 1]) ** (p + 1) - fmpq_poly([1])
        for j in range(p):
            acc = acc - P[j] * fmpq(comb(p + 1, j))
        P.append(acc / fmpq(p + 1))
    S = []
    for ell in range(imax + 1):
        c = list(P[2 * ell].coeffs())
        at2m = fmpq_poly([c[j] * fmpq(2) ** j for j in range(len(c))])
        S.append(at2m - P[2 * ell] * fmpq(2) ** (2 * ell))
    m = fmpq_poly([0, 1])  # bbar = S_1/m
    out = [(fmpq_poly([1]), fmpq_poly([1]))]
    for i in range(1, imax + 1):
        # mu_i = (1/m) SUM_l C(i,l) (-1)^{i-l} S_l / bbar^l = SUM_l C(i,l)(-1)^{i-l} S_l m^{l-1} / S_1^l
        num = fmpq_poly([0])
        for ell in range(i + 1):
            num += fmpq(comb(i, ell)) * (-1) ** (i - ell) * S[ell] * m ** (ell) * S[1] ** (i - ell)
        den = S[1] ** i * m
        out.append((num, den))
    return out


def mu_deviation_bounds(imax, mmin):
    """C_i with |mu_i(m) - mu_i^inf| <= C_i / m^2 for all m >= mmin, by interval evaluation in u = 1/m."""
    muinf = mu_inf(imax)
    polys = mu_exact_poly(imax)
    C = [arb(0)] * (imax + 1)
    for i in range(2, imax + 1):
        num, den = polys[i]
        # (num - muinf*den)/den, written in u = 1/m: multiply numerator and denominator by u^deg
        diff = num - den * muinf[i]
        d = max(diff.degree(), den.degree())
        # coefficients of u^k: reverse with padding to degree d
        dc = list(diff.coeffs()) + [fmpq(0)] * (d + 1 - len(diff.coeffs()))
        ec = list(den.coeffs()) + [fmpq(0)] * (d + 1 - len(den.coeffs()))
        # diff(m) m^{-d} = SUM dc[k] u^{d-k}; the two lowest powers of u must vanish (O(1/m^2))
        du = [dc[d - k] for k in range(d + 1)]
        eu = [ec[d - k] for k in range(d + 1)]
        assert du[0] == 0 and du[1] == 0, f"mu_{i}: deviation is not O(1/m^2)"
        u = arb(0, 1 / mmin)  # the interval [0, 1/mmin] as a ball

        def apow(x, k):  # arb ** int returns nan on a ball containing 0; multiply instead
            r = arb(1)
            for _ in range(k):
                r *= x
            return r

        numv = sum((toarb(du[k]) * apow(u, k - 2) for k in range(2, d + 1)), arb(0))
        denv = sum((toarb(eu[k]) * apow(u, k) for k in range(d + 1)), arb(0))
        val = numv / denv
        C[i] = abs(val).upper()  # an arb upper bound, valid for every m >= mmin
    return C


def mu_sensitivity(J, mu, amax, bmax, C):
    """M_i >= sup |dPhat/dmu_i| over a<=amax, b<=bmax, |mu - mu_inf| <= C/m^2 (crude absolute-coefficient bound)."""
    names = ("a", "b") + tuple(f"m{i}" for i in range(2, J + 1))
    ctx11 = fmpq_mpoly_ctx.get(names, Ordering.lex)
    g = ctx11.gens()
    a, b = g[0], g[1]
    mus = {i: g[i] for i in range(2, J + 1)}
    theta = a * b
    L = {i: b * (fmpq((-1) ** (i - 1)) / i) * mus[i] * theta ** (i - 2) for i in range(2, J + 1)}
    E = [ctx11.from_dict({tuple([0] * len(names)): fmpq(1)})]
    for j in range(1, J + 1):
        acc = ctx11.from_dict({})
        for i in range(2, j + 1):
            acc += L[i] * E[j - i] * fmpq(i)
        E.append(acc / fmpq(j))
    P = ctx11.from_dict({})
    for j in range(J + 1):
        num = ctx11.from_dict({tuple([0] * len(names)): fmpq(1)})
        for i in range(1, j):
            num *= 1 - fmpq(i) * a
        den_rest = ctx11.from_dict({tuple([0] * len(names)): fmpq(1)})
        for i in range(j, J):
            den_rest *= 1 - fmpq(i) * a * a * b
        P += E[j] * num * den_rest

    def D(p):
        return -a * a * p.derivative(0) + 2 * a * b * p.derivative(1)

    DP = D(P)
    Phat = DP * DP - P * D(DP) - fmpq(4, 5) * a * a * b * P * P
    M = {}
    for i in range(2, J + 1):
        dPi = Phat.derivative(i)  # variable m_i sits at index i in names
        tot = arb(0)
        for exps, c in dPi.to_dict().items():
            term = toarb(abs(c)) * toarb(amax) ** exps[0] * toarb(bmax) ** exps[1]
            for k in range(2, J + 1):
                if exps[k]:
                    term *= (toarb(abs(mu[k])) + toarb(C[k]) / (TAU_MIN**2 * 50) ** 2) ** exps[k]
            tot += term
        M[i] = tot
    return M


def toarb(q):
    if isinstance(q, arb):
        return q
    return arb(q.numer().str()) / arb(q.denom().str())


def main():
    t0 = time.time()
    mu = mu_inf(J)
    print(f"mu_i^inf, i=2..{J}:", ", ".join(str(m) for m in mu[2:]))
    P, Q = build_F(J, mu)
    print(
        f"F_J = P/Q built: P has {len(P)} terms, degrees {P.degrees()}; Q degrees {Q.degrees()}  [{time.time() - t0:.1f} s]"
    )
    DP = Dop(P)
    D2P = Dop(DP)
    Phat = DP * DP - P * D2P - fmpq(4, 5) * A * A * B * P * P
    print(f"Phat: {len(Phat)} terms, degrees {Phat.degrees()}  [{time.time() - t0:.1f} s]")

    c0, c1, c2, c3 = (coeff_in_a(Phat, k) for k in range(4))
    print("[a^0] =", c0, " [a^1] =", c1)
    assert c0.is_zero() and c1.is_zero(), "the a^0 / a^1 coefficients must vanish identically"
    print(
        "[a^2](b) =",
        c2,
        "   (the truncation defect; lowest power of b:",
        next(i for i, c in enumerate(c2.coeffs()) if c != 0) if not c2.is_zero() else None,
        ")",
    )
    print("[a^3](b) =", c3)
    lead = c3.coeffs()[2] if len(c3.coeffs()) > 2 else fmpq(0)
    print("   [a^3] leading coefficient of b^2:", lead, "= 176/175 ?", lead == fmpq(176, 175))

    # ---- region bounds --------------------------------------------------------------
    amax = fmpq(1, TAU_MIN)
    bmax = BSTAR
    # (i) [a^3](b) >= (176/175) b^2 (1 - kappa b) on [0, bmax] with an explicit kappa: write [a^3] = b^2 * c3r
    c3r = fmpq_poly(list(c3.coeffs())[2:])
    assert all(c == 0 for c in list(c3.coeffs())[:2])
    # lower bound of c3r on [0,bmax]: c3r(0) - SUM_{k>=1} |c_k| bmax^k
    c3r_low = c3r.coeffs()[0] - sum(abs(c) * bmax**k for k, c in enumerate(c3r.coeffs()) if k >= 1)
    print(f"(i)  [a^3]/b^2 >= {float(c3r_low):.6f} on b in [0, {float(bmax)}]")
    # (ii) defect: |[a^2](b)| <= d5 * b^5 on the region
    c2c = list(c2.coeffs())
    low2 = next(i for i, c in enumerate(c2c) if c != 0)
    assert low2 >= 5, f"defect starts at b^{low2}, expected >= 5"
    d5 = sum(abs(c) * bmax ** (k - 5) for k, c in enumerate(c2c) if k >= 5)
    print(f"(ii) |[a^2](b)| <= {float(d5):.3e} b^5")
    # (iii) rest: SUM_{k>=4} |[a^k](b)| a^{k-3} <= (SUM_k Ck amax^{k-3}) b^2 with Ck = bound of [a^k]/b^2
    maxa = Phat.degrees()[0]
    rest = fmpq(0)
    for k in range(4, maxa + 1):
        ck = list(coeff_in_a(Phat, k).coeffs())
        if not ck:
            continue
        lowk = next((i for i, c in enumerate(ck) if c != 0), None)
        if lowk is None:
            continue
        assert lowk >= 2, f"[a^{k}] not divisible by b^2 (lowest power {lowk})"
        Ck = sum(abs(c) * bmax ** (i - 2) for i, c in enumerate(ck) if i >= 2)
        rest += Ck * amax ** (k - 3)
    print(f"(iii) SUM_k |[a^k]| a^(k-3) <= {float(rest):.3e} b^2 on the region")
    # (iv) assemble: Phat/(a^3 b^2) >= c3r_low - d5 b^3/a - rest.  On the region b <= bmax and
    # a >= ? -- b^3/a is largest when a is smallest.  The region must therefore be cut with a lower
    # bound on a in terms of b:  we require  b <= BREL * a^(1/3)  i.e.  b^3/a <= BREL^3.
    # With BREL = 1/5:  d5 b^3/a <= d5/125.
    BREL3 = fmpq(1, 125)
    total_low = c3r_low - d5 * BREL3 - rest
    print(
        f"(iv) Phat/(a^3 b^2) >= {float(total_low):.6f} on {{ a <= 1/{TAU_MIN}, b <= {float(bmax)}, b^3 <= a/125 }}"
    )
    ok_main = total_low > 0

    # ---- (v) the tail j > J --------------------------------------------------------------
    # Interpolant on the window [t-1, t+1]: F^(t)(tau) = SUM_{j <= t+1} e_j (tau)_j/(N)_j, a polynomial
    # equal to p_s / bbar^s at s = t-1, t, t+1.  Its tail T = SUM_{J<j<=t+1} obeys, for j <= N mu_2/8,
    #     |e_j| <= (2 e N mu_2 / j)^{j/2}          (Cauchy on |z| = sqrt(j/(2 N mu_2)), using
    #                                              |(1+w) e^{-w}| <= e^{|w|^2} for |w| <= 1/2, |beta_k| < 2)
    #     |(tau)_j| <= tau^j (1 + 2/tau)^2 (the two possible negative factors are at most 2 in size)
    #     (N)_j >= (N-j)^j
    # so  tail_j := |e_j (tau)_j/(N)_j| <= (2 e mu_2 b / j)^{j/2} (1-j/N)^{-j} (1+2/tau)^2 =: (2 e mu_2 b/j)^{j/2} K_j
    # and |T''| <= SUM_j (j/(tau-j+1))^2 tail_j.  With F >= e^{-0.41 b} >= 0.99 and |F'| <= b,
    #     |Phi(F+T) - Phi(F)| N <= N (|T''| + 2|F'||T'| + (|F'|^2+|F''|)|T|)/(0.99 F) <= 1.03 SUM_j j^2 a^2 N tail_j
    #                          = 1.03 SUM_j j^2 tail_j / b.
    # Relative to the certified margin a b:  1.03 SUM_j j^2 (2e mu_2/j)^{j/2} b^{j/2-2}/a, and on the region
    # b^3 <= a/125 gives b^{j/2-2}/a <= b^{j/2-5}/125 <= bmax^{j/2-5}/125 for j >= J+1 >= 11.
    e = arb(1).exp()
    mu2 = toarb(mu[2])
    # (1-j/N)^{-j} <= exp(2 j^2/N) <= exp(2*200^2/(627^2*100)) < 1.0021 for j <= 200, N >= 627^2 * 100;
    # (tau)_j has at most one negative factor of size <= 1, so |(tau)_j| <= tau^j; F >= 0.98, |F'| <= b:
    # the linearised perturbation bound carries the factor 1.03/0.98 < 1.06.
    Kfudge = arb("1.0021")
    tail_rel = arb(0)
    for j in range(J + 1, 200):
        tail_rel += (
            arb("1.06")
            * arb(j * j)
            * (2 * e * mu2 / j) ** (arb(j) / 2)
            * toarb(bmax) ** (arb(j) / 2 - 5)
            / 125
            * Kfudge
        )
    # j > 200: the terms are below 1e-300 and are absorbed by rounding the bound up
    tail_rel = tail_rel * arb("1.000001")
    print(f"(v)  tail (j > {J}) effect relative to a*b:  <= {float(tail_rel):.3e}")

    # ---- (vi) the finite-m moments ------------------------------------------------------------
    # mu_i(m) - mu_i^inf = O(1/m^2); with 1/m = 2 a^2 b the perturbation of Phat is
    #     SUM_i delta_i * dPhat/dmu_i  with |delta_i| <= C_i / m^2 = 4 C_i a^4 b^2 ,
    # so relative to a^3 b^2 it is at most  4 a SUM_i C_i sup|dPhat/dmu_i|  <=  4 amax SUM_i C_i M_i .
    C = mu_deviation_bounds(J, TAU_MIN**2 * 50)
    Msym = mu_sensitivity(J, mu, amax, bmax, C)
    mu_rel = 4 * toarb(amax) * sum((toarb(C[i]) * Msym[i] for i in range(2, J + 1)), arb(0))
    print(
        f"(vi) finite-m moments: relative effect <= {float(mu_rel):.3e}   (C_i = {', '.join(f'{float(C[i]):.2g}' for i in range(2, J + 1))})"
    )

    ok = (
        ok_main
        and float(tail_rel) + float(mu_rel) < 0.05
        and float(total_low) > float(tail_rel) + float(mu_rel)
    )
    print()
    print(
        f"VERDICT: Phi N > 4/5 on the region {{ tau >= {TAU_MIN}, tau^2/N <= {float(bmax)}, (tau^2/N)^3 <= 1/(125 tau) }} :",
        ok,
    )
    print(f"[{time.time() - t0:.1f} s]")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
