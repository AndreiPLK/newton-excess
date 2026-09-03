"""The dense regime, part (a): 0.05 <= theta <= 1/2, t >= 628, N >= 1260 -- by the fixed-tilt Fourier
weight, the exact Edgeworth polynomial of its second cumulant, and polygamma closed forms of the
tilted cumulants with a Stirling remainder.  Certified interval arithmetic throughout (arb/acb).

STATEMENT CERTIFIED.  For every real x (the tilted mean) and N = 2m with
    N >= 1260,  x >= 628,  0.05 <= theta := x/N <= 1/2 :
        N * g_lower(x, N)  >  4/5 ,
where g_lower <= g(n,t) := -Delta^2 log p_t at t = x whenever x is an integer.  Since
log(1 + 4/(5n)) < 4/(5n) < 4/(5N), this gives M(n,t) = n(e^g - 1) > n g > 4/5.

THE OBJECTS.  Tilt r fixed at the saddle of x (SUM_k q_k = x).  For real tau,
    P_r(tau) = (1/2pi) INT_{-pi}^{pi} e^{K(phi) - i tau phi} dphi,   K = SUM_k 2 log(1 + q_k(e^{i phi}-1)),
agrees with e_tau r^tau / E(r) at integers, and (log P_r)'' = -kappa_2[phi], the second cumulant of
phi under the complex weight.  Hence
    -Delta^2 log e_t = INT_0^1 INT_0^1 kappa_2[phi](t-1+u+v) du dv  >=  min_{|tau - t| <= 1} kappa_2[phi](tau),
    g = -Delta^2 log e_t + Delta^2 log C(N,t) = INT INT kappa_2[phi] - L_bin,
    L_bin := log[ (t+1)(N-t+1) / (t (N-t)) ]     (exact).
With u = sigma phi, sigma^2 = k_2, c_1 = (x - tau)/sigma in [-1/sigma, 1/sigma], c_j = k_j/(j! sigma^j):
    kappa_2[phi] = kappa_2[u]/sigma^2,   kappa_2[u] = M_2/M_0 - (M_1/M_0)^2,
    M_j = INT_{|u| <= pi sigma} u^j e^{-u^2/2 + rho(u)} du / sqrt(2 pi),  rho(u) = i c_1 u + SUM_{j>=3} c_j (iu)^j.
M_j is enclosed by the weight-<= W truncation of exp(rho) integrated against the Gaussian (exact
polynomial in the c's, weights wt(c_1) = 1, wt(c_j) = j-2) plus an explicit remainder:
truncation of exp inside |u| <= L, the tail of the cumulant series, and the region |u| > L where
|e^{K(phi)}| <= exp(-k_2 (1 - cos phi)).

THE CUMULANTS.  With a = m + 1/2, y = m/v, z = a + i y, w = eps z = 1 + i/v + eps/2, Y = 1/v:
    k_1 = 2[ m + y Im psi(z) - (pi y/2) tanh(pi y) ],     k_{p+1} = -(y/2) d/dy k_p ,
    kt_{p+1} := k_{p+1}/N = delta_{p0} + SUM_q c_{p,q} [ Y^{q+1} Im(i^q Phi_q) + q Y^q Im(i^{q-1} Phi_{q-1}) ]
                            - (pi/2)(-1/2)^p Y ,
    Phi_q := psi^{(q)}(z) / eps^q  (Phi_0 := psi(z) + log eps), the c_{p,q} from T^p = SUM_q c_{p,q} y^q d^q,
    T = -(y/2) d/dy:  c_{p+1,q} = -(q/2) c_{p,q} - (1/2) c_{p,q-1}.
The Phi_q are Stirling partial sums in powers of eps with the remainder
    |R_n| <= |B_{2K+2}| (n+2K+1)! / ((2K+2)! (Re z)^{n+2K+2}),
from psi^{(n)}(z) = (-1)^{n+1} INT_0^inf t^n e^{-zt}/(1-e^{-t}) dt and the Bernoulli-remainder lemma
| t/(1-e^{-t}) - 1 - t/2 - SUM_{k<=K} B_{2k} t^{2k}/(2k)! | <= |B_{2K+2}| t^{2K+2}/(2K+2)!  (t > 0),
which this script certifies by a 1-D interval sweep before use.

The sweep is an adaptive bisection over (eps, v) boxes with eps = 1/m in (0, 1/630] and v in the
range where theta_0(v) = 1 - arctan(v)/v covers [0.05, 1/2]; boxes entirely outside the region
(x_hi < 628, theta_hi < 0.05 or theta_lo > 1/2) are skipped.

Run:  uv run python projects/qg-bootstrap/release/scripts/dense_certificate_a.py [max_depth]
"""

from __future__ import annotations

import sys
import time
from math import factorial

from flint import acb, arb, ctx, fmpq

ctx.prec = 200
PI = arb.pi()

W = 6  # exact Edgeworth weight kept
JC = 8  # cumulants k_1..k_JC by polygamma; c_9..c_13 by sup bounds; j >= 14 by Cauchy
KST = 6  # Stirling terms
BERN = {
    2: fmpq(1, 6),
    4: fmpq(-1, 30),
    6: fmpq(1, 42),
    8: fmpq(-1, 30),
    10: fmpq(5, 66),
    12: fmpq(-691, 2730),
    14: fmpq(7, 6),
    16: fmpq(-3617, 510),
}
SUPK = {  # sup_{q in [0,1]} |kappa_j(q)| (certified separately in lab/, rounded up)
    3: 0.0963,
    4: 0.1251,
    5: 0.1277,
    6: 0.2501,
    7: 0.4084,
    8: 1.0626,
    9: 2.3901,
    10: 7.7501,
    11: 22.253,
    12: 86.376,
    13: 302.05,
}


def toarb(q):
    if isinstance(q, arb):
        return q
    return arb(q.numer().str()) / arb(q.denom().str())


def apow(x, e):
    """x**e by repeated multiplication: arb ** int returns nan on a ball containing 0."""
    r = arb(1)
    for _ in range(e):
        r *= x
    return r


# ------------------------------------------------------------------ Bernoulli remainder lemma
def bernoulli_lemma(K, X0=60.0, depth=40):
    """Certify  |(t/2)coth(t/2) - SUM_{k<=K} B_{2k} t^{2k}/(2k)!| <= |B_{2K+2}| t^{2K+2}/(2K+2)!  on (0, X0],
    and prove it for t > X0 by the crude bound (t/2)coth(t/2) <= t/2 + 1."""
    B = [fmpq(1)] + [BERN[2 * k] for k in range(1, K + 2)]

    def S(t):
        return sum((toarb(B[k]) * t ** (2 * k) / factorial(2 * k) for k in range(K + 1)), arb(0))

    def bound(t):
        return abs(toarb(B[K + 1])) * t ** (2 * K + 2) / factorial(2 * K + 2)

    def f(t):
        return (t / 2) / (t / 2).tanh()

    stack = [(arb(0), arb(X0), 0)]
    boxes = 0
    while stack:
        lo, hi, d = stack.pop()
        boxes += 1
        t = arb((lo + hi) / 2, (hi - lo) / 2)
        if t.lower() <= 0:
            # near 0 use the series bound: the function minus S is the tail SUM_{k>K} B_2k t^2k/(2k)!,
            # alternating with decreasing terms for t < 2 pi, hence bounded by its first term.
            if hi < 6:
                continue
        ok = False
        if t.lower() > 0:
            ok = bool(abs(f(t) - S(t)) < bound(t))
        if ok:
            continue
        if d > depth:
            return False, boxes
        mid = (lo + hi) / 2
        stack.append((lo, mid, d + 1))
        stack.append((mid, hi, d + 1))
    # t > X0:  |f - S| <= t/2 + 1 + SUM_k |B_2k| t^{2k}/(2k)!  <=  bound(t)  -- check at X0 with monotone growth
    t = arb(X0)
    crude = (
        t / 2
        + 1
        + sum((abs(toarb(B[k])) * t ** (2 * k) / factorial(2 * k) for k in range(K + 1)), arb(0))
    )
    # the ratio bound(t)/crude(t) is increasing in t (bound has the highest power); enough at t = X0
    return bool(crude < bound(t)), boxes


# ------------------------------------------------------------------ Stirling forms of psi^{(n)}(z)/eps^n
def phi_stirling(n, eps, w, K):
    """Phi_n = psi^{(n)}(z)/eps^n with z = w/eps, Re z = (1+eps/2)/eps; includes the remainder ball."""
    if n == 0:
        val = w.log() - eps / (2 * w)
        for k in range(1, K + 1):
            val -= toarb(BERN[2 * k]) * apow(eps, 2 * k) / (2 * k * w ** (2 * k))
    else:
        val = acb(factorial(n - 1)) / w**n + acb(factorial(n)) * eps / (2 * w ** (n + 1))
        for k in range(1, K + 1):
            val += (
                toarb(BERN[2 * k])
                * factorial(2 * k + n - 1)
                * apow(eps, 2 * k)
                / (factorial(2 * k) * w ** (2 * k + n))
            )
        val = val * ((-1) ** (n + 1))
    # remainder: |R_n| <= |B_{2K+2}| (n+2K+1)! / ((2K+2)! Re(z)^{n+2K+2}),  Re z = (1 + eps/2)/eps
    # scaled by eps^{-n}:  |R_n|/eps^n <= |B| (n+2K+1)!/(2K+2)! * eps^{2K+2} / (1+eps/2)^{n+2K+2}   (no division by eps)
    rem = (
        abs(toarb(BERN[2 * K + 2]))
        * factorial(n + 2 * K + 1)
        / factorial(2 * K + 2)
        * apow(eps, 2 * K + 2)
        / apow(1 + eps / 2, n + 2 * K + 2)
    )
    r = rem.upper()
    return val + acb(arb(0, r), arb(0, r))


def scaled_cumulants(eps, v, jmax, K=KST):
    """kt_j = k_j / N for j = 1..jmax as arb enclosures, for eps = 1/m and v in boxes."""
    w = acb(1 + eps / 2, 1 / v)
    Y = 1 / v
    Phi = [phi_stirling(q, eps, w, K) for q in range(jmax)]

    # I_q := Im(i^q Phi_q), except that I_0 carries the -pi/2 of the -(pi y/2) term:
    # Im Phi_0 - pi/2 = arg(w) - pi/2 + Im(Phi_0 - log w) = -arctan(v (1 + eps/2)) + Im(Phi_0 - log w),
    # which removes the cancellation between y Im psi and pi y/2 that ruins interval enclosures.
    def Im_iq(q):
        if q == 0:
            return -(v * (1 + eps / 2)).atan() + (Phi[0] - w.log()).imag
        r = acb(0, 1) ** q * Phi[q]
        return r.imag

    Iq = [Im_iq(q) for q in range(jmax)]
    # c_{p,q}
    c = {(0, 0): fmpq(1)}
    for p in range(0, jmax - 1):
        for q in range(0, p + 2):
            val = fmpq(0)
            if (p, q) in c:
                val += -fmpq(q, 2) * c[(p, q)]
            if (p, q - 1) in c:
                val += -fmpq(1, 2) * c[(p, q - 1)]
            if val != 0:
                c[(p + 1, q)] = val
    kt = []
    for p in range(jmax):
        acc = arb(1) if p == 0 else arb(0)
        for q in range(p + 1):
            if (p, q) not in c:
                continue
            term = Y ** (q + 1) * Iq[q]
            if q >= 1:
                term += q * Y**q * Iq[q - 1]
            acc += toarb(c[(p, q)]) * term
        # (the -(pi/2)(-1/2)^p Y term is already inside I_0, see Im_iq)
        # tanh correction: (pi y/2)(1 - tanh(pi y)) and its T-derivatives are below e^{-1600}
        acc += arb(0, 1e-300)
        kt.append(acc)
    return kt


# ------------------------------------------------------------------ Edgeworth engine (weights)
JMAX = 13
VARS = [1] + list(range(3, JMAX + 1))
WT = {1: 1, **{j: j - 2 for j in range(3, JMAX + 1)}}
NV = len(VARS)
ZERO = tuple([0] * NV)


def weight(mon):
    return sum(WT[VARS[i]] * mon[i] for i in range(NV))


def padd(a, b):
    r = dict(a)
    for k, val in b.items():
        r[k] = r.get(k, fmpq(0)) + val
        if r[k] == 0:
            del r[k]
    return r


def pmul(a, b, Wmax):
    r = {}
    for ka, va in a.items():
        wa = weight(ka)
        for kb, vb in b.items():
            if wa + weight(kb) > Wmax:
                continue
            k = tuple(x + y for x, y in zip(ka, kb, strict=True))
            r[k] = r.get(k, fmpq(0)) + va * vb
    return {k: val for k, val in r.items() if val != 0}


def var(j):
    m = [0] * NV
    m[VARS.index(j)] = 1
    return tuple(m)


def truncated_exp(Wmax):
    """exp(rho) with rho = SUM c_j U^j (U = i u), truncated by weight <= Wmax: {U-degree: c-poly}."""
    rho = {1: {var(1): fmpq(1)}}
    for j in range(3, JMAX + 1):
        if WT[j] <= Wmax:
            rho[j] = {var(j): fmpq(1)}
    expo = {0: {ZERO: fmpq(1)}}
    power = {0: {ZERO: fmpq(1)}}
    for n in range(1, Wmax + 1):
        newp = {}
        for da, pa in power.items():
            for db, pb in rho.items():
                p = pmul(pa, pb, Wmax)
                if p:
                    newp[da + db] = padd(newp.get(da + db, {}), p)
        power = newp
        if not power:
            break
        for d, p in power.items():
            expo[d] = padd(expo.get(d, {}), {k: val / factorial(n) for k, val in p.items()})
    return expo


def gauss_abs_moment(d):
    """INT_R |u|^d e^{-u^2/2} du / sqrt(2 pi) = 2^{d/2} Gamma((d+1)/2) / sqrt(pi)."""
    return arb(2) ** (arb(d) / 2) * (arb(d + 1) / 2).gamma() / PI.sqrt()


def gauss_abs_moment_tail(d, L):
    """INT_{|u|>L} |u|^d e^{-u^2/2} du / sqrt(2 pi) = 2^{d/2} Gamma((d+1)/2, L^2/2) / sqrt(pi)."""
    # python-flint: x.gamma_upper(s) == Gamma(s, x)  (checked: arb(2).gamma_upper(arb(3)) = Gamma(3,2))
    return arb(2) ** (arb(d) / 2) * (L * L / 2).gamma_upper(arb(d + 1) / 2) / PI.sqrt()


def dfact2(k):
    r = 1
    for i in range(k - 1, 0, -2):
        r *= i
    return r


def eval_cpoly(p, cabs):
    """SUM coef * prod c_j^e with c_j given as arb (signed or absolute)."""
    tot = arb(0)
    for k, val in p.items():
        term = toarb(val)
        for i, e in enumerate(k):
            if e:
                term *= apow(cabs[VARS[i]], e)
        tot += term
    return tot


def eval_cpoly_abs(p, cabs):
    tot = arb(0)
    for k, val in p.items():
        term = abs(toarb(val))
        for i, e in enumerate(k):
            if e:
                term *= apow(cabs[VARS[i]], e)
        tot += term
    return tot


EXPO = truncated_exp(W)
DEBUG = False
kappa2_enclosure_last_delta = None


def moments_in_c1(cvals):
    """The weight-<=W Gaussian moments M_0, A_1 (M_1 = i A_1), M_2 as polynomials in c_1:
    returns three lists indexed by the power of c_1, coefficients arb (the other c_j from cvals)."""
    out = []
    for j in range(3):
        acc = {}
        for d, p in EXPO.items():
            if (d + j) % 2:
                continue
            g = arb(dfact2(d + j))
            k = d % 4
            sign = {0: 1, 2: -1, 1: 1, 3: -1}[k]
            for mon, val in p.items():
                e1 = mon[0]
                term = toarb(val) * g * sign
                for i, e in enumerate(mon):
                    if i == 0 or e == 0:
                        continue
                    term *= apow(cvals[VARS[i]], e)
                acc[e1] = acc.get(e1, arb(0)) + term
        out.append(acc)
    return out


def poly_avg(coeffs_by_power, inv_sigma):
    """Window average of SUM_e coeff_e (delta inv_sigma)^e over delta with density (1-|delta|) on [-1,1]:
    E[delta^e] = 2/((e+1)(e+2)) for even e, 0 for odd e."""
    tot = arb(0)
    for e, c in coeffs_by_power.items():
        if e % 2:
            continue
        tot += c * apow(inv_sigma, e) * fmpq(2, (e + 1) * (e + 2))
    return tot


def kappa2_window_average(cvals, cabs, sigma_lo, half_eps, k2t, L, inv_sigma, geo_K=None):
    """Enclosure of the window average of kappa_2[u], with the c_1 dependence integrated exactly.

    kappa_2[u] = M_2/M_0 + A_1^2/M_0^2 with M_j = M_j^(W)(c_1) + e_j, |e_j| <= delta_j (from kappa2_enclosure's
    remainder with c_1 = the full ball).  Write M_0 = 1 + h, |h| <= hbar < 1, 1/M_0 = SUM_{k<=3} (-h)^k + R_4,
    |R_4| <= hbar^4/(1 - hbar).  Then kappa_2[u] is a polynomial in delta = c_1 sigma of degree <= 4W plus
    ball terms, and the average over the window is taken term by term.
    """
    # remainders delta_j with c_1 as the full ball (they are monotone in |c_1|, so this is an upper bound)
    _, M0ball = kappa2_enclosure(cvals, cabs, sigma_lo, half_eps, k2t, L, geo_K=geo_K)
    dj = kappa2_enclosure.last_delta  # set by kappa2_enclosure
    M0c, A1c, M2c = moments_in_c1(cvals)

    # polynomial arithmetic in delta (dict power -> arb), c_1 = delta * inv_sigma folded into the coefficients
    def scale(cd):
        return {e: c * apow(inv_sigma, e) for e, c in cd.items()}

    M0p, A1p, M2p = scale(M0c), scale(A1c), scale(M2c)
    h = dict(M0p)
    h[0] = h.get(0, arb(0)) - 1

    def pmul_d(a, b):
        r = {}
        for ea, ca in a.items():
            for eb, cb in b.items():
                r[ea + eb] = r.get(ea + eb, arb(0)) + ca * cb
        return r

    def padd_d(a, b, s=1):
        r = dict(a)
        for e, c in b.items():
            r[e] = r.get(e, arb(0)) + s * c
        return r

    # bound hbar = max over |delta| <= 1 of |h(delta)| + delta_0
    hbar = sum((abs(c) for c in h.values()), arb(0)) + dj[0]
    if not (hbar.upper() < 0.5):
        return None
    inv = {0: arb(1)}
    hp = {0: arb(1)}
    for k in range(1, 4):
        hp = pmul_d(hp, h)
        inv = padd_d(inv, hp, (-1) ** k)
    # 1/M_0 with M_0 = 1 + h + e_0: inv3 is built from h alone, so the omitted piece
    # inv3(h + e_0) - inv3(h) is bounded by delta_0/((1-hbar)(1-hbar-delta_0)) (validator finding, 3 Sept 2026)
    R4 = apow(hbar, 4) / (1 - hbar) + dj[0] / ((1 - hbar) * (1 - hbar - dj[0]))
    # kappa_2 polynomial part:  M2 * inv + (A1 * inv)^2
    A1inv = pmul_d(A1p, inv)
    kpoly = padd_d(pmul_d(M2p, inv), pmul_d(A1inv, A1inv))
    avg = poly_avg(kpoly, arb(1))  # inv_sigma already folded in
    # ball terms: |M2| |R4| + |e_2| |1/M0| + |A1 inv|^2 corrections ... bounded crudely with sup norms over |delta|<=1
    sup = lambda cd: sum((abs(c) for c in cd.values()), arb(0))  # noqa: E731
    invsup = sup(inv) + R4
    A1sup = sup(A1p) + dj[1]
    err = (
        sup(M2p) * R4
        + dj[2] * invsup
        + (2 * A1sup * dj[1] + dj[1] * dj[1]) * invsup * invsup
        + (sup(A1p) ** 2) * (2 * invsup * R4 + R4 * R4)
    )
    return avg + arb(0, err.upper())


def kappa2_enclosure(cvals, cabs, sigma_lo, half_eps, k2t, L, geo_K=None):
    """Enclosure of kappa_2[u] from the weight-<=W Gaussian moments plus explicit remainders.

    cvals: dict j -> arb enclosure of c_j (j = 1 interval [-1/sigma, 1/sigma], 3..JC from polygamma,
           JC+1..13 zero-centred balls from sup bounds);  cabs: absolute bounds |c_j| for j <= 13.
    geo_K: optional (Kg, rho) for the j >= 14 tail in the sigma-free form
           SUM_{j>=14} |c_j| |u|^j <= Kg sigma^2 SUM_{j>=14} (|u|/(rho sigma))^j
           (used by the dense-(b) corner, where the default N/sigma^14 form is not sigma-free; the default
           is kept for the boxes of this script and the bands of dense_certificate_b).
    """

    # M_j^(W) for j = 0,1,2: SUM_d [i^d (d+j-1)!!] coef_d  (real for even d+j; M_1 = i * A_1)
    def MW(j):
        re, im = arb(0), arb(0)
        for d, p in EXPO.items():
            if (d + j) % 2:
                continue
            g = arb(dfact2(d + j))
            val = eval_cpoly(p, cvals) * g
            k = d % 4
            if k == 0:
                re += val
            elif k == 2:
                re -= val
            elif k == 1:
                im += val
            else:
                im -= val
        return re, im

    M0, _ = MW(0)
    _, A1 = MW(1)
    M2, _ = MW(2)

    # ---- remainder delta_j = |T_W outer| + |R_W inner| + |outer true|
    # rho-bar(u) = SUM |c_j| |u|^j (j = 1, 3..13) + 2.7 N SUM_{j>=14} (|u|/(3 sigma))^j
    # (the constant 2.7/3^j for |kappa_j|/N comes from the Hadamard-zero identity
    #  kappa_j(q) = -(j-1)! SUM_k s_k^{-j}, giving sup_q |kappa_j| <= 2 (j-1)! (1 - 2^-j) zeta(j)/pi^j,
    #  36x below 2.7/3^j for j >= 14; the circle Cauchy estimate alone does not give it -- validator, 3 Sept 2026)
    def rhobar_poly_coeffs():
        co = {1: cabs[1]}
        for j in range(3, JMAX + 1):
            co[j] = cabs[j]
        return co

    co = rhobar_poly_coeffs()
    # 2.7 N (L/(3 sigma))^14 / (1 - L/(3 sigma)) with N/sigma^14 = (1/N)^6 / kt_2^7 = half_eps^6 / k2t^7
    # upper bound: the largest 1/N and the smallest kt_2 on the box (a wide k2t ball raised to the 7th power
    # would otherwise swallow 0)
    if geo_K is None:
        geo_tail_at_L = (
            arb("2.7")
            * apow(arb(half_eps.upper()), 6)
            / apow(arb(k2t.lower()), 7)
            * (L / 3) ** 14
            / (1 - L / (3 * sigma_lo))
        )
    else:
        # sigma-free form: Kg sigma^2 (L/(rho sigma))^14 / (1 - L/(rho sigma)), with the smallest sigma
        # giving the largest value when L <= sigma/2 (the power 2 - 14 < 0 of sigma)
        Kg, rho = geo_K
        ratio = L / (rho * sigma_lo)
        assert ratio.upper() < 1
        geo_tail_at_L = Kg * sigma_lo * sigma_lo * apow(ratio, 14) / (1 - ratio)

    # (a) truncation of exp inside |u| <= L: for n <= W the weight > W part of rhobar^n, for n > W all.
    # weight>W part of rhobar^n  =  rhobar^n - (rhobar restricted to c_1..c_{W+2})^n truncated at weight W.
    # Both are polynomials in |u|; we integrate them against |u|^j e^{-u^2/2} over R (an upper bound).
    def poly_pow_coeffs(coeffs, n):
        # coeffs: dict degree -> arb ; returns dict for the n-th power
        out = {0: arb(1)}
        for _ in range(n):
            new = {}
            for da, ca in out.items():
                for db, cb in coeffs.items():
                    new[da + db] = new.get(da + db, arb(0)) + ca * cb
            out = new
        return out

    delta = [arb(0), arb(0), arb(0)]
    full_co = dict(co)
    full_co[14] = (
        geo_tail_at_L / L**14
    )  # the geometric tail bounded by one monomial of degree 14 on |u|<=L
    # weight-<=W part of exp with absolute coefficients, as {degree: abs-c-poly value}
    for n in range(1, W + 1):
        pw = poly_pow_coeffs(full_co, n)
        # the kept part: monomials of weight <= W built from c_1..c_{W+2} -- take from EXPO-like structure:
        # we recompute the weight-truncated n-th power with the dict engine
        rho_small = {1: {var(1): fmpq(1)}}
        for j in range(3, W + 3):
            rho_small[j] = {var(j): fmpq(1)}
        power = {0: {ZERO: fmpq(1)}}
        for _ in range(n):
            newp = {}
            for da, pa in power.items():
                for db, pb in rho_small.items():
                    p = pmul(pa, pb, W)
                    if p:
                        newp[da + db] = padd(newp.get(da + db, {}), p)
            power = newp
        for j in range(3):
            tot = arb(0)
            for d, cf in pw.items():
                kept = eval_cpoly_abs(power.get(d, {}), cabs) if d in power else arb(0)
                dropped = cf - kept  # >= 0 : the weight > W part (absolute coefficients)
                tot += dropped * gauss_abs_moment(d + j)
            delta[j] += tot / factorial(n)
            if DEBUG:
                print(f"   (a) n={n} j={j}: {(tot / factorial(n)).str(4)}")
    # n > W: SUM_{n>W} rhobar^n/n! <= rhobar(u)^{W+1} e^{rhobar(u)}/(W+1)!.  On |u| <= L every term of degree
    # j >= 3 obeys |u|^j <= L^{j-2} u^2, so rhobar(u) <= c_1 |u| + beta u^2 with beta = SUM_{j>=3} |c_j| L^{j-2},
    # and e^{rhobar} <= e^{c_1 L} e^{beta u^2}.  The polynomial rhobar^{W+1} is kept EXACTLY and integrated against
    # the wider Gaussian e^{-(1/2 - beta) u^2}:  INT |u|^k e^{-gamma u^2} du / sqrt(2pi) = Gamma((k+1)/2) / (gamma^{(k+1)/2} sqrt(2 pi)).
    beta = sum((co[j] * apow(L, j - 2) for j in co if j >= 3), arb(0)) + geo_tail_at_L / (L * L)
    gamma = arb(1) / 2 - beta
    if not (gamma.lower() > arb("0.1")):
        for j in range(3):
            delta[j] += arb("inf")
    else:
        pwW1 = poly_pow_coeffs(full_co, W + 1)
        pref = (cabs[1] * L).exp() / factorial(W + 1)
        for j in range(3):
            tot = arb(0)
            for d, cf in pwW1.items():
                k = d + j
                tot += cf * (arb(k + 1) / 2).gamma() / (gamma ** (arb(k + 1) / 2) * (2 * PI).sqrt())
            delta[j] += pref * tot
            if DEBUG:
                print(f"   (a') n>W j={j}: {(pref * tot).str(4)}   (beta = {beta.str(4)}, L = {L})")
    # (b) the truncated polynomial T_W outside |u| > L
    for j in range(3):
        tot = arb(0)
        for d, p in EXPO.items():
            tot += eval_cpoly_abs(p, cabs) * gauss_abs_moment_tail(d + j, L)
        delta[j] += tot
        if DEBUG:
            print(f"   (b) j={j}: {tot.str(4)}")
    # (c) the true integrand outside |u| > L:  |e^{K}| <= exp(-k_2 (1 - cos phi)), phi = u/sigma.
    #     for |phi| <= 1: 1 - cos phi >= 0.4583 phi^2 ; for 1 < |phi| <= pi: 1 - cos phi >= 0.4596.
    for j in range(3):
        # |u| in (L, sigma]:  INT |u|^j e^{-0.4583 u^2} du / sqrt(2pi) ; substitute s = sqrt(0.9166) u
        s0 = arb("0.9166").sqrt()
        part1 = gauss_abs_moment_tail(j, L * s0) / s0 ** (j + 1)
        # |u| in (sigma, pi sigma]: length <= 2 pi sigma, |u|^j <= (pi sigma)^j, weight e^{-0.4596 k_2};
        # sigma^{j+1} e^{-0.4596 sigma^2} is decreasing for sigma^2 > (j+1)/0.9192, so sigma_lo is the worst
        part2 = (
            2
            * PI
            * sigma_lo
            * (PI * sigma_lo) ** j
            * (-arb("0.4596") * sigma_lo * sigma_lo).exp()
            / (2 * PI).sqrt()
        )
        delta[j] += part1 + part2
        if DEBUG:
            print(f"   (c) j={j}: part1 {part1.str(4)} part2 {part2.str(4)}")
    kappa2_enclosure.last_delta = [d.upper() for d in delta]
    M0e = M0 + arb(0, delta[0].upper())
    A1e = A1 + arb(0, delta[1].upper())
    M2e = M2 + arb(0, delta[2].upper())
    # kappa_2[u] = M2/M0 - (i A1/M0)^2 = M2/M0 + A1^2/M0^2
    return M2e / M0e + (A1e / M0e) ** 2, M0e


# ------------------------------------------------------------------ one box
def certify_box(eps, v, verbose=False):
    """Return (status, info): status in {'ok', 'skip', 'fail'} for the box eps x v (arb balls).

    Every quantity is written through eps = 1/m and the scaled cumulants kt_j = k_j/N, so that a
    box with eps -> 0 (N -> infinity) is handled without infinities: N^{1-j/2} = (eps/2)^{j/2-1}.
    """
    # eps box clamped to [0, upper]
    if eps.lower() < 0:
        eps = arb(eps.upper() / 2, eps.upper() / 2)
    # mean-value enclosure in v:  d kt_j / dv = (2/v) kt_{j+1}  (exact: d/dv = (2/v) T at fixed m),
    # so kt_j(v) = kt_j(v_c) + (v - v_c) (2/xi) kt_{j+1}(xi).  This kills the dependency blow-up of the
    # T-recursion when v is a box: the centre is evaluated on a thin interval.
    v_c = arb(v.mid())
    kt_c = scaled_cumulants(eps, v_c, JC + 2)
    raw = scaled_cumulants(
        eps, v, JC + 2
    )  # plain box evaluation, loose but only used under a factor dv
    dv = arb(0, v.rad())
    kt_b = [kt_c[j] + dv * (2 / v) * raw[j + 1] for j in range(JC + 1)]  # first mean-value pass
    kt = [kt_c[j] + dv * (2 / v) * kt_b[j + 1] for j in range(JC)]  # second pass, tight
    theta = kt[0]
    eps_lo = eps.lower()
    N_lo = 2 / eps.upper()
    # x = N theta >= N_lo * theta_lo ; skip only if even the largest x in the box is below 628
    x_hi = (2 / eps_lo) * theta.upper() if eps_lo > 0 else arb("inf")
    if x_hi < 628 or theta.upper() < arb("0.05") or theta.lower() > arb("0.5"):
        return "skip", {}
    k2t = kt[1]
    if not (k2t.lower() > 0):
        return "fail", {"why": "k2 not positive"}
    half_eps = eps / 2  # = 1/N
    sigma_lo = (N_lo * k2t.lower()).sqrt()
    # 1/sigma = sqrt(half_eps / k2t): the ball may touch 0 (eps -> 0), so enclose [0, sqrt(upper)] by hand
    q_up = (half_eps / k2t).upper()
    if (half_eps / k2t).lower() > 0:
        inv_sigma = (half_eps / k2t).sqrt()
    else:
        r_up = q_up.sqrt().upper()
        inv_sigma = arb(r_up / 2, r_up / 2)
    cvals = {1: arb(0, inv_sigma.upper())}
    cabs = {1: arb(inv_sigma.upper())}
    for j in range(3, JC + 1):
        # c_j = N^{1-j/2} kt_j / (j! kt_2^{j/2}) = (eps/2)^{j/2-1} kt_j / (j! kt_2^{j/2})
        cj = apow(inv_sigma, j - 2) * kt[j - 1] / (factorial(j) * k2t)
        cvals[j] = cj
        cabs[j] = abs(cj).upper()
    for j in range(JC + 1, JMAX + 1):
        bnd = (apow(inv_sigma, j - 2) * SUPK[j] / (factorial(j) * k2t)).upper()
        cvals[j] = arb(0, bnd)
        cabs[j] = arb(bnd)
    # inner cut L: every candidate gives a valid enclosure; keep the tightest (all bounds are rigorous).
    # The choice is made with c_1 as the full ball [-1/sigma, 1/sigma].
    best = None
    for Lc in (4, 5, 6, 7, 8, 10, 12, 15, 20, 30):
        if Lc > sigma_lo.lower():
            break
        k2u, M0c = kappa2_enclosure(cvals, cabs, sigma_lo, half_eps, k2t, arb(Lc))
        if not (M0c.lower() > 0):
            continue
        if best is None or k2u.rad() < best[0].rad():
            best = (k2u, M0c, Lc)
    if best is None:
        return "fail", {"why": "no valid split", "sigma_lo": sigma_lo}
    _, M0, L = best
    kappa2u = kappa2_window_average(cvals, cabs, sigma_lo, half_eps, k2t, arb(L), inv_sigma)
    if kappa2u is None:
        return "fail", {"why": "window average: h too large"}
    # N g >= (N/k_2) kappa_2[u] - N L_bin,  N L_bin = (1/theta) G(x) + (1/(1-theta)) G(N-x),
    # G(x) := x log(1 + 1/x) in [1 - 1/(2x), 1], increasing; x >= 628, N - x >= N/2.  Written without
    # the cancellation between N/k_2 and N/(x(N-x)) through the variance identity:
    #   N g = NV/(kt_2 Q) + (kappa_2[u] - 1)/kt_2 + (1 - G(x))/theta + (1 - G(N-x))/(1 - theta),
    #   NV := theta(1-theta) - kt_2 = Var(q) >= 0,  Q = theta(1-theta).
    # NV is enclosed by the mean value theorem in v with its own derivative (2/v)[kt_2(1-2kt_1) - kt_3],
    # which is O(theta^2), so the enclosure does not lose the factor 1/theta^2.
    NV_c = kt_c[0] * (1 - kt_c[0]) - kt_c[1]
    NV = NV_c + dv * (2 / v) * (kt_b[1] * (1 - 2 * kt_b[0]) - kt_b[2])
    Q = theta * (1 - theta)
    x_lo = arb(max(628.0, float(N_lo * theta.lower())))
    one_minus_Gx = arb(0, (1 / (2 * x_lo)).upper())  # 1 - G(x) in [0, 1/(2 x_lo)]
    one_minus_Gnx = arb(0, (1 / N_lo).upper())  # 1 - G(N-x) in [0, 1/(2 (N/2))]
    Ng = NV / (k2t * Q) + (kappa2u - 1) / k2t + one_minus_Gx / theta + one_minus_Gnx / (1 - theta)
    info = {"theta": theta, "sigma_lo": sigma_lo, "NV": NV, "kappa2u": kappa2u, "Ng": Ng, "L": L}
    if Ng.lower() > arb("0.8"):
        return "ok", info
    return "fail", info


def main():
    max_depth = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    t0 = time.time()
    ok, nb = bernoulli_lemma(KST)
    print(f"Bernoulli remainder lemma (K={KST}): {ok}  [{nb} boxes]")
    assert ok

    # v range: theta_0(v) in [0.049, 0.501]
    def theta0(v):
        return 1 - v.atan() / v

    vlo = arb("0.38")  # theta_0(0.38) = 0.0466 < 0.05
    vhi = arb("2.34")  # theta_0(2.34) = 0.5003 > 0.5
    assert theta0(vlo) < arb("0.05") and theta0(vhi) > arb("0.5")
    eps_max = arb(1) / 630
    stack = [(arb(0), eps_max, vlo, vhi, 0)]
    n_ok = n_skip = 0
    worst = None
    while stack:
        e1, e2, v1, v2, d = stack.pop()
        eps = arb((e1 + e2) / 2, ((e2 - e1) / 2).upper())
        v = arb((v1 + v2) / 2, ((v2 - v1) / 2).upper())
        try:
            status, info = certify_box(eps, v)
        except Exception as ex:  # noqa: BLE001 -- an arithmetic failure is a failed box, bisect it
            status, info = "fail", {"why": repr(ex)}
        if status == "ok":
            n_ok += 1
            if worst is None or info["Ng"].lower() < worst[0]:
                worst = (info["Ng"].lower(), float(eps.mid()), float(v.mid()))
            continue
        if status == "skip":
            n_skip += 1
            continue
        if d >= max_depth:
            print("FAILED box:", float(e1), float(e2), float(v1), float(v2), info)
            sys.exit(1)
        # bisect the wider relative side: v (its width matters most), then eps
        if (v2 - v1) > (e2 - e1) * 500 or d % 3 != 2:
            vm = (v1 + v2) / 2
            stack.append((e1, e2, v1, vm, d + 1))
            stack.append((e1, e2, vm, v2, d + 1))
        else:
            em = (e1 + e2) / 2
            stack.append((e1, em, v1, v2, d + 1))
            stack.append((em, e2, v1, v2, d + 1))
        if (n_ok + n_skip) % 500 == 0 and n_ok + n_skip > 0:
            print(
                f"  ... {n_ok} ok, {n_skip} skipped, stack {len(stack)}  [{time.time() - t0:.0f} s]",
                flush=True,
            )
    print(
        f"certified: {n_ok} boxes ok, {n_skip} skipped, worst N g lower bound {float(worst[0]):.6f} at eps={worst[1]:.2e}, v={worst[2]:.4f}"
    )
    print(
        f"VERDICT: N g > 4/5 on the dense region (a) {{ 0.05 <= theta <= 1/2, x >= 628, N >= 1260 }} : True  [{time.time() - t0:.0f} s]"
    )


if __name__ == "__main__":
    main()
