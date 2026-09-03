"""SIBUYA'S CONJECTURE (1988, eq. 3.4), the dense regime (b'): theta <= 0.1 and t^2/N >= 1, uniformly as
theta -> 0, by an exact power series in v (the tilt) whose coefficients are polynomials in zeta = eta/v.

Port of dense_certificate_b.py (the 4/5 theorem, series in V = v^2 for the doubled odd squares) to the
spectrum {1..N}.  What changed, and only that:
  * tilt r = v/N, q_k = (k v/N)/(1 + k v/N);  eps = 1/N (NOT 1/m = 2/N as in the original: every "2" that
    came from N = 2/eps is gone: 1/N = eps = zeta v^2, t = khat_1/(zeta v) = khat_1/eta, sigma^2 = khat_2/eta,
    1/x = zeta v/khat_1, 1/(N-x) = zeta v^2/(1 - v khat_1), c_j = eta^{j/2-1} khat_j/(j! khat_2^{j/2}));
  * the series variable is v itself (kt_j = v khat_j, khat_j(0) = 1/2), the Faulhaber factor is the single
    sum  sigma_i(eps) = eps^{i+1} SUM_{k<=N} k^i  (a polynomial in eps with a FIRST-order term), and the
    Riemann bound is sigma_i(eps) <= (1+eps)^{i+1}/(i+1);
  * the target: (3 + 3 eps - theta) N g > 1 (Sibuya's n(R-1) >= n/(3n-t), see sibuya_dense_a.py), so the
    regular series is  S := ((3 + 3 zeta v^2 - v khat_1) N g - 1)/v  whose constant term must vanish exactly;
    every loss on N g enters S through the factor |3 + 3 eps - theta| <= 3.01;
  * REGION: 0 < v <= VMAX = 0.24 (theta <= 0.1 => v <= 0.24 since theta(0.24, N) >= 0.1037 for every N),
    0 <= zeta <= ZMAX = 0.26 (N < (j+1)^2 => zeta = khat_1^2/b <= 0.2506 with khat_1 <= (1+eps)/2), zeta v <= ETAMAX = 5.1e-4
    (j >= 1001 => eta = khat_1/j <= 4.9955e-4).  Constants independently validated 3 Sept 2026
    (results/VALIDATION_SIBUYA_DENSE_B_2026-09-03.md, seven items PASS).

The original docstring follows (its structure holds with the substitutions above).

The dense regime, part (b): theta <= 0.05 and t^2/N >= 0.95, uniformly as theta -> 0,
by an exact power series in V = v^2 whose coefficients are polynomials in zeta = eta/V.

VARIABLES.  v = 2 m sqrt(r) (the tilt), V = v^2, eps = 1/m, eta = eps/V, zeta = eta/V, so
    eps = zeta V^2,   1/N = eps/2 = zeta V^2/2,   theta = kt_1 = V khat_1,   t = x = N theta = 2 khat_1/(zeta V),
    b = t theta = 2 khat_1^2/zeta,   sigma^2 = N kt_2 = 2 khat_2/(zeta V) .
REGION certified:  0 < V <= VMAX = 0.17,  0 <= zeta <= ZMAX = 0.24,  zeta V <= ETAMAX = 1.07e-3.
It contains {theta <= 0.05, t^2/N >= 0.95, t >= 628}: theta <= 0.05 gives V <= 0.17 (theta <= V/3 hmm,
theta = V khat_1 with khat_1 <= 1/3), b >= 0.95 gives zeta <= 2 khat_1^2/0.95 <= 0.234, t >= 628 gives
eta = 2 khat_1/t <= 1.07e-3.

STATEMENT CERTIFIED.   N g > 4/5  on the region, hence M(n,t) > 4/5 there (g >= min over the window
of kappa_2[phi] minus the exact binomial term, as in dense_certificate_a.py).

THE SERIES.  Every ingredient is a power series in V with coefficients in Q[zeta]:
  * kt_j = k_j/N = SUM_i [SUM_p a_{j,p} c_{p,i}] V^i sigma_i(eps),  kappa_j(q) = SUM_p a_{j,p} q^p (Bernoulli
    cumulant polynomials), q^p = SUM_i c_{p,i} u^{2i} with c_{p,i} = (-1)^{i-p} C(i-1, p-1), and
    sigma_i(eps) = 2^{-2i} eps^{2i+1} S_i(1/eps) with S_i(m) = SUM_{k<=m} (2k-1)^{2i} (Faulhaber, exact),
    a polynomial in eps = zeta V^2.   khat_j := kt_j / V  has constant term 1/3.
  * NV := theta(1-theta) - kt_2 = V^2 nhat,  nhat = (khat_1 - khat_2)/V - khat_1^2,  nhat(0) = 4/45;
    NV/(kt_2 Q) = nhat / (khat_2 khat_1 (1 - V khat_1)),  constant term 4/5.
  * Edgeworth: kappa_2[u] - 1 = SUM over weights w = 2, 4, ..., W of the exact polynomials in the c_j
    (lab/edgeworth_exact.py engine), with  c_j = (eta/2)^{j/2-1} khat_j / (j! khat_2^{j/2})  (j >= 3),
    c_1 = delta sqrt(eta/(2 khat_2)),  delta = x - tau in [-1, 1].  In a monomial of weight w with
    e_1 powers of c_1 and factors c_{j_l}^{e_l} (# = SUM e_l) the half powers pair up:
        delta^{e_1} zeta^{w/2} V^{w/2} 2^{-w/2} khat_2^{-(w/2 + #)} PROD (khat_{j_l}/j_l!)^{e_l},
    and the window average replaces delta^{e_1} by 2/((e_1+1)(e_1+2)) (even e_1) or 0.
    Then (kappa_2[u]_avg - 1)/kt_2 = (...)/(V khat_2).
  * binomial term:  N L_bin = (1/theta) G(x) + (1/(1-theta)) G(N-x),  G(x) = x log(1 + 1/x);
    1 - G(x) = SUM_{k>=1} (-1)^{k+1} (1/x)^k/(k+1),  1/x = zeta V/(2 khat_1);  1/(N-x) = zeta V^2/(2(1 - V khat_1)).
  N g = NV/(kt_2 Q) + (kappa_2 avg - 1)/kt_2 + (1 - G(x))/theta + (1 - G(N-x))/(1-theta).
  The zeta^k V^0 terms cancel exactly against the binomial term for k <= W/2 (checked: the constant
  term of N g is 4/5 in Q[zeta]), so  S := (N g - 4/5)/V  is a regular series, S(0, zeta) = 176/525 (zeta-free).

TAILS.  Each series carries a bound M(r) on its supremum over the disc |V| <= r, |zeta| <= ZMAX
(r = 0.3 > VMAX): sums, products, 1/(c + g) with M_g < |c|, and division by V (maximum principle)
propagate it.  The truncation of S at degree I is then bounded by M_S(r) (V/r)^{I+1}/(1 - V/r) (Cauchy).
The Edgeworth remainder (weights > W, the exp tail, the cumulant tail, the outer region) is bounded
with dense_certificate_a.kappa2_enclosure on V-bands, as a decreasing function of sigma, and shown to
be below 0.05 V^2 khat_2 on every band; the binomial log series is cut with its first omitted term.

Run:  uv run python projects/qg-bootstrap/release/scripts/sibuya_dense_b.py
"""

from __future__ import annotations

import os
import sys
import time
from math import comb, factorial

from flint import arb, ctx, fmpq, fmpq_poly

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dense_certificate_a as DA  # noqa: E402  (the Edgeworth engine and the remainder machinery)

ctx.prec = 300

I_TRUNC = (
    150  # v-degree kept (VMAX/R_DISC = 0.767 needs (0.767)^(I-2) M_S << 1, M_S ~ 2.5e12: I >= 140)
)
W = 10  # Edgeworth weight kept exactly
# the remainder machinery of dense_certificate_a must truncate at the SAME weight as the series here
DA.W = W
DA.EXPO = DA.truncated_exp(W)
VMAX = fmpq(24, 100)  # theta <= 0.1 => v <= 0.23 (khat_1 >= 0.4365 there); 0.24 leaves room
ZMAX = fmpq(26, 100)  # b >= 1 => zeta = khat_1^2/b <= (1+eps)^2/4 <= 0.2502
ETAMAX = fmpq(51, 100000)  # t >= 1001 => eta = khat_1/t <= (1+eps)/2002 <= 4.997e-4
R_DISC = fmpq(
    3, 10
)  # Cauchy radius in v (0.4 fails the inverse gap of khat_1: sup of the rest 0.53 > 1/2)
FAC = arb("3.01")  # |3 + 3 eps - theta| <= 3.01 on the region
JCUM = 13  # cumulants needed: c_j up to j = JMAX of the Edgeworth engine


def toarb(q):
    if isinstance(q, arb):
        return q
    return arb(q.numer().str()) / arb(q.denom().str())


def poly_sup(p, zmax):
    """sup_{|zeta| <= zmax} |p(zeta)| bounded by SUM |c_k| zmax^k."""
    return sum((abs(toarb(c)) * toarb(zmax) ** k for k, c in enumerate(p.coeffs())), arb(0))


# ------------------------------------------------------------------ series with a sup bound
class Ser:
    """Truncated power series in V with fmpq_poly (in zeta) coefficients, a certified sup bound M on the
    polydisc |V| <= R_DISC, |zeta| <= ZMAX (an arb upper bound), the exact degree d (coefficients of V^n are
    exact for n <= d; d drops by one under divV, and a product is exact to the smaller d), and a real-point
    perturbation P: a polynomial in V with non-negative arb coefficients such that, on the real region, the
    TRUE quantity differs from the analytic object described by (c, M) by at most P(V).  P starts with the
    i > I remainder of the base series kt_j (order V^{I+1}) and propagates through every operation
    (validator findings of 3 Sept 2026).  Only the operations used below are provided."""

    __slots__ = ("c", "M", "d", "P")
    DINF = 10**9

    def __init__(self, coeffs, M, d=DINF, P=None):
        self.c = list(coeffs) + [fmpq_poly([0])] * (I_TRUNC + 1 - len(coeffs))
        self.M = M
        self.d = d
        self.P = dict(P) if P else {}

    @staticmethod
    def const(q):
        return Ser([fmpq_poly([q])], abs(toarb(q)))

    @staticmethod
    def from_terms(terms, extra_tail=None, P=None):
        """terms: dict (V-power) -> fmpq_poly in zeta, exact for powers <= I_TRUNC; extra_tail: the polydisc
        sup of the terms of V-power > I_TRUNC of the analytic base object; P: real-point perturbation."""
        if extra_tail is None:
            extra_tail = arb(0)
        cs = [fmpq_poly([0])] * (I_TRUNC + 1)
        for i, p in terms.items():
            if i <= I_TRUNC:
                cs[i] = cs[i] + p
        M = (
            sum((poly_sup(cs[i], ZMAX) * toarb(R_DISC) ** i for i in range(I_TRUNC + 1)), arb(0))
            + extra_tail
        )
        return Ser(cs, M, I_TRUNC, P)

    # ---- perturbation polynomials (dict power -> arb, non-negative)
    @staticmethod
    def _padd(a, b):
        out = dict(a)
        for k, v in b.items():
            out[k] = out[k] + v if k in out else v
        return out

    @staticmethod
    def _pscale(a, s):
        return {k: v * s for k, v in a.items()}

    @staticmethod
    def _pmul(a, b):
        out = {}
        for i, u in a.items():
            for j, w in b.items():
                out[i + j] = out[i + j] + u * w if (i + j) in out else u * w
        return out

    @staticmethod
    def _peval(a, Vhi):
        return sum((v * DA.apow(arb(Vhi), k) for k, v in a.items()), arb(0))

    def __add__(self, o):
        if not isinstance(o, Ser):
            o = Ser.const(o)
        return Ser(
            [a + b for a, b in zip(self.c, o.c, strict=True)],
            self.M + o.M,
            min(self.d, o.d),
            Ser._padd(self.P, o.P),
        )

    __radd__ = __add__

    def __neg__(self):
        return Ser([-a for a in self.c], self.M, self.d, self.P)

    def __sub__(self, o):
        if not isinstance(o, Ser):
            o = Ser.const(o)
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, Ser):
            q = fmpq(o) if not isinstance(o, fmpq) else o
            return Ser(
                [a * q for a in self.c],
                self.M * abs(toarb(q)),
                self.d,
                Ser._pscale(self.P, abs(toarb(q))),
            )
        out = [fmpq_poly([0]) for _ in range(I_TRUNC + 1)]
        for i, a in enumerate(self.c):
            if a.is_zero():
                continue
            for j in range(I_TRUNC + 1 - i):
                b = o.c[j]
                if b.is_zero():
                    continue
                out[i + j] = out[i + j] + a * b
        # |(f+df)(g+dg) - fg| <= |f| dg + |g| df + df dg, with |f| <= M_f on the real region (in the polydisc)
        P = Ser._padd(
            Ser._padd(Ser._pscale(o.P, self.M), Ser._pscale(self.P, o.M)), Ser._pmul(self.P, o.P)
        )
        return Ser(out, self.M * o.M, min(self.d, o.d), P)

    __rmul__ = __mul__

    def inv(self):
        """1/self, requiring a nonzero rational constant term c0 (zeta-free) and M(self - c0) < |c0|."""
        c0 = self.c[0]
        assert c0.degree() == 0 and not c0.is_zero(), "inverse needs a nonzero constant"
        c0q = c0.coeffs()[0]
        rest = Ser([fmpq_poly([0])] + self.c[1:], arb(0))
        Mrest = sum(
            (poly_sup(rest.c[i], ZMAX) * toarb(R_DISC) ** i for i in range(1, I_TRUNC + 1)), arb(0)
        )
        # tail of self beyond I contributes to M(self - c0) as well: bound by self.M - |c0| ... conservative:
        Mrest = Mrest + (
            self.M
            - sum(
                (poly_sup(self.c[i], ZMAX) * toarb(R_DISC) ** i for i in range(I_TRUNC + 1)), arb(0)
            )
        )
        gap = abs(toarb(c0q)) - Mrest
        assert gap.lower() > 0, f"inverse: sup of the rest {Mrest} not below |c0| = {c0q}"
        # exact coefficients by recursion: (c0 + g) h = 1  ->  h_0 = 1/c0,  h_i = -(1/c0) SUM_{k=1}^{i} g_k h_{i-k}
        h = [fmpq_poly([0]) for _ in range(I_TRUNC + 1)]
        h[0] = fmpq_poly([1 / c0q])
        for i in range(1, I_TRUNC + 1):
            acc = fmpq_poly([0])
            for k in range(1, i + 1):
                if not self.c[k].is_zero():
                    acc = acc + self.c[k] * h[i - k]
            h[i] = -acc * (1 / c0q)
        # |1/(f+df) - 1/f| <= |df| / (|f| (|f| - |df|)) with |f| >= gap and |df| <= P(VMAX) on the real region
        pmax = Ser._peval(self.P, float(VMAX))
        assert (gap - pmax).lower() > 0, "inverse: perturbation not below the gap"
        P = Ser._pscale(self.P, 1 / (gap * (gap - pmax)))
        return Ser(h, 1 / gap, self.d, P)

    def divV(self):
        """self / V, requiring a zero constant term; M by the maximum principle: M/r; exact degree drops by one;
        the perturbation P(V)/V is again a polynomial with non-negative coefficients when P has no constant term."""
        assert self.c[0].is_zero(), "division by V needs a zero constant term"
        assert all(k >= 1 for k in self.P), "division by V of a perturbation with a constant term"
        return Ser(
            self.c[1:] + [fmpq_poly([0])],
            self.M / toarb(R_DISC),
            self.d - 1,
            {k - 1: v for k, v in self.P.items()},
        )

    def pow(self, e):
        r = Ser.const(fmpq(1))
        for _ in range(e):
            r = r * self
        return r

    def evaluate(self, V, zeta):
        """Interval enclosure at V (arb box in [0, VMAX]) and zeta (arb box): the exact polynomial part to degree
        d, the Cauchy tail M (V/r)^{d+1}/(1 - V/r) of the analytic object, and the real-point perturbation P(V)."""
        dd = min(self.d, I_TRUNC)
        tot = arb(0)
        for i in range(dd, -1, -1):
            pv = arb(0)
            for ck in reversed(list(self.c[i].coeffs())):
                pv = pv * zeta + toarb(ck)
            tot = tot * V + pv
        ratio = V / toarb(R_DISC)
        tail = self.M * DA.apow(ratio, dd + 1) / (1 - ratio)
        pert = Ser._peval(self.P, V.upper())
        return tot + arb(0, tail.upper()) + arb(0, pert.upper())


# ------------------------------------------------------------------ ingredients
def faulhaber_S(L):
    P = [fmpq_poly([0, 1])]
    for p in range(1, L + 1):
        acc = fmpq_poly([1, 1]) ** (p + 1) - fmpq_poly([1])
        for j in range(p):
            acc = acc - P[j] * fmpq(comb(p + 1, j))
        P.append(acc / fmpq(p + 1))
    return P[: L + 1]  # single spectrum: S_l(N) = SUM_{k<=N} k^l, all l


def sigma_polys(imax):
    """sigma_i(eps) = eps^{i+1} S_i(1/eps) as polynomials in eps (fmpq_poly), i <= imax (single spectrum)."""
    S = faulhaber_S(imax)
    out = []
    for i in range(imax + 1):
        cs = list(S[i].coeffs())  # S_i(N) = SUM_d cs[d] N^d, degree i+1, cs[0] = 0
        pe = [fmpq(0)] * (i + 2)
        for d, c in enumerate(cs):
            if c != 0:
                pe[i + 1 - d] += c
        out.append(fmpq_poly(pe))
    return out


def bernoulli_cumulant_polys(jmax):
    kp = [None, fmpq_poly([0, 1])]
    q1q = fmpq_poly([0, 1, -1])
    for j in range(1, jmax):
        kp.append(q1q * kp[j].derivative())
    return kp


def cumulant_series(jmax):
    """kt_j as Ser (j = 1..jmax), with the tail bound from |coef_i| <= A_j 2^i sup|sigma_i|."""
    sig = sigma_polys(I_TRUNC + 1)
    kp = bernoulli_cumulant_polys(jmax)
    out = []
    zmax, rr = toarb(ZMAX), toarb(R_DISC)
    for j in range(1, jmax + 1):
        a = list(kp[j].coeffs())  # a[p] = coefficient of q^p, p = 1..j
        terms = {}
        Aj = sum((abs(toarb(c)) for c in a), arb(0))
        for i in range(1, I_TRUNC + 1):
            coef = fmpq(0)
            for p in range(1, min(i, j) + 1):
                coef += a[p] * (-1) ** (i - p) * comb(i - 1, p - 1)
            if coef == 0:
                continue
            # coef * V^i * sigma_i(eps), eps = zeta V^2:  sigma_i(eps) = SUM_k s_k eps^k -> s_k zeta^k V^{i+2k}
            for k, sk in enumerate(sig[i].coeffs()):
                if sk == 0 or i + 2 * k > I_TRUNC:
                    continue
                terms[i + 2 * k] = terms.get(i + 2 * k, fmpq_poly([0])) + fmpq_poly(
                    [fmpq(0)] * k + [coef * sk]
                )
        # THE BASE OBJECT IS THE POLYNOMIAL G_j := SUM_{i<=I} coef_i V^i sigma_i(zeta V^2), sigma_i the Faulhaber
        # polynomials (validator finding of 3 Sept 2026: the full series with Faulhaber polynomials is NOT bounded
        # on a complex eps-disc for large i, so it is never used as an analytic object).  G_j is a polynomial in
        # (V, zeta), analytic everywhere; its polydisc sup is the sum of the kept coefficients' sups (from_terms)
        # plus tail 2 below (its own terms of V-degree > I).  On the real region the true kt_j equals G_j plus
        # T_j := SUM_{i>I} coef_i V^i sigma_i(eps) with eps = 1/m an integer point, where the Riemann-sum bound
        # S_i(m) = SUM_{k<=m} (2k-1)^{2i} <= INT_{1/2}^{m+1/2} (2x)^{2i} dx <= (2m+1)^{2i+1}/(2(2i+1)) gives
        # sigma_i(eps) <= (1+eps/2)^{2i+1}/(2i+1), and |coef_i| <= A_j 2^{i-1}:
        #     |T_j(V)| <= A_j (1+eps/2)/(2(2I+3)) (2(1+eps/2)^2)^{I+1} V^{I+1} / (1 - q),  q = 2 VMAX (1+eps/2)^2,
        # which is carried as the real-point perturbation P (a multiple of V^{I+1}) and propagated by Ser.
        eps_max = zmax * rr * rr  # >= every real eps = zeta V^2 on the region
        # single spectrum: sigma_i(eps) <= (1+eps)^{i+1}/(i+1)  (S_i(N) <= (N+1)^{i+1}/(i+1))
        qq = 2 * toarb(VMAX) * (1 + eps_max)
        assert qq.upper() < 1
        P_base = {
            I_TRUNC + 1: Aj
            * (1 + eps_max)
            / (2 * (I_TRUNC + 2))
            * (2 * (1 + eps_max)) ** (I_TRUNC + 1)
            / (1 - qq)
        }
        tail = arb(0)
        # tail 2: the terms of G_j of V-degree > I (i <= I, i + 2k > I), polydisc sup |coef_i| |s_{i,k}| ZMAX^k r^{i+2k}
        # = |coef_i| r^i |s_{i,k}| eps_max^k  (the earlier version carried a spurious extra factor r^{2k};
        # removed 3 Sept 2026, validator finding)
        for i in range(1, I_TRUNC + 1):
            coef = fmpq(0)
            for p in range(1, min(i, j) + 1):
                coef += a[p] * (-1) ** (i - p) * comb(i - 1, p - 1)
            if coef == 0:
                continue
            for k, sk in enumerate(sig[i].coeffs()):
                if sk == 0 or i + 2 * k <= I_TRUNC:
                    continue
                tail += abs(toarb(coef)) * rr**i * abs(toarb(sk)) * eps_max**k
        out.append(Ser.from_terms(terms, extra_tail=tail, P=P_base))
    return out


def window_moment(e):
    return fmpq(2, (e + 1) * (e + 2)) if e % 2 == 0 else fmpq(0)


def build_Ng(kt):
    """N g as a Ser, and the pieces."""
    V = Ser([fmpq_poly([0]), fmpq_poly([1])], toarb(R_DISC))
    khat = [k.divV() for k in kt]  # khat[j-1] = kt_j / V
    k1, k2 = khat[0], khat[1]
    theta = V * k1
    one = Ser.const(fmpq(1))
    # NV/(kt_2 Q)
    nhat = (k1 - k2).divV() - k1 * k1
    inv_k1 = k1.inv()
    inv_k2 = k2.inv()
    inv_1mt = (one - theta).inv()
    piece1 = nhat * inv_k2 * inv_k1 * inv_1mt
    # Edgeworth (weights <= W), window-averaged
    k2_derive, _ = DA_derive_kappa2()
    edge = Ser.const(fmpq(0))
    zeta = Ser([fmpq_poly([0, 1])], toarb(ZMAX))
    for mon, coef in k2_derive.items():
        if mon == DA.ZERO:
            continue
        e1 = mon[0]
        w = DA.weight(mon)
        if w > W or w % 2:
            continue
        me = window_moment(e1)
        if me == 0:
            continue
        nf = sum(mon[1:])
        term = Ser.const(
            coef * me
        )  # no 2^{-w/2}: here 1/N = eps, so c_j carries eta^{j/2-1}, not (eta/2)^{j/2-1}
        term = term * zeta.pow(w // 2) * V.pow(w // 2) * inv_k2.pow(w // 2 + nf)
        for idx, e in enumerate(mon):
            if idx == 0 or e == 0:
                continue
            j = DA.VARS[idx]
            term = term * (khat[j - 1] * fmpq(1, factorial(j))).pow(e)
        edge = edge + term
    piece2 = edge.divV() * inv_k2
    # binomial: (1 - G(x))/theta, 1/x = zeta V/(2 khat_1)
    invx = zeta * V * inv_k1  # 1/x = 1/t = eps/theta = zeta v/khat_1
    oneG = Ser.const(fmpq(0))
    for k in range(1, 6):
        oneG = oneG + invx.pow(k) * fmpq((-1) ** (k + 1), k + 1)
    piece3 = oneG.divV() * inv_k1
    invNx = zeta * V * V * inv_1mt  # 1/(N - x) = eps/(1 - theta)
    oneGn = Ser.const(fmpq(0))
    for k in range(1, 4):
        oneGn = oneGn + invNx.pow(k) * fmpq((-1) ** (k + 1), k + 1)
    piece4 = oneGn * inv_1mt
    Ng = piece1 + piece2 + piece3 + piece4
    return Ng, (piece1, piece2, piece3, piece4), khat, invx, invNx


_cache = {}


def DA_derive_kappa2():
    if "k2" not in _cache:
        sys.path.insert(
            0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lab")
        )
        import edgeworth_exact as EE  # noqa: E402

        k2, k4 = EE.derive(W)
        # remap monomial tuples from EE's VARS to DA's VARS (identical layout: [1, 3..13])
        assert EE.VARS == DA.VARS
        _cache["k2"] = (k2, k4)
    return _cache["k2"]


def main():
    t0 = time.time()
    kt = cumulant_series(JCUM)
    print(f"cumulant series kt_1..kt_{JCUM} built, I = {I_TRUNC}  [{time.time() - t0:.1f} s]")
    print(
        "  kt_1 =",
        " + ".join(f"({kt[0].c[i]})V^{i}" for i in range(1, 4)),
        "...   M =",
        kt[0].M.str(5),
    )
    print("  kt_2 =", " + ".join(f"({kt[1].c[i]})V^{i}" for i in range(1, 4)), "...")
    Ng, pieces, khat, invx, invNx = build_Ng(kt)
    print(f"N g built  [{time.time() - t0:.1f} s]")
    c0 = Ng.c[0]
    print("  constant term of N g (must be 1/3 exactly, zeta-free):", c0)
    assert c0.degree() == 0 and c0.coeffs()[0] == fmpq(1, 3), "the zeta^k v^0 terms did not cancel"
    Vs = Ser([fmpq_poly([0]), fmpq_poly([1])], toarb(R_DISC))
    zs = Ser([fmpq_poly([0, 1])], toarb(ZMAX))
    target_fac = fmpq(3) + zs * Vs * Vs * fmpq(3) - Vs * khat[0]  # 3 + 3 eps - theta
    Ngs = target_fac * Ng - fmpq(1)
    c0s = Ngs.c[0]
    print("  constant term of (3 + 3 eps - theta) N g - 1 (must be 0 exactly):", c0s)
    assert c0s.is_zero(), "the constant term of the Sibuya sign series did not cancel"
    S = Ngs.divV()
    print(
        "  S = ((3 + 3 eps - theta) N g - 1)/v:  S(0, zeta) =",
        S.c[0],
        "   v^1 coefficient:",
        S.c[1],
    )
    print(
        "  sup bound M_S(r = 0.3) =",
        S.M.str(6),
        "  exact degree d =",
        S.d,
        "  real-point perturbation P_S(VMAX) =",
        Ser._peval(S.P, float(VMAX)).str(3),
    )
    # ---- remainders not in the series ----------------------------------------------------------
    # (r1) binomial log series cut: |(1/x)^6/7| / theta and |(1/(N-x))^4/5| / (1-theta), relative to V
    # (r2) Edgeworth remainder on V-bands via DA.kappa2_enclosure, expressed as a bound on |kappa_2 - 1 - poly|/(kt_2 V)
    bands = []
    v_hi = float(VMAX)
    while v_hi > 1e-6:
        bands.append((v_hi / 2, v_hi))
        v_hi /= 2
    worst_S = None
    for lo, hi in bands:
        Vb = arb((lo + hi) / 2, (hi - lo) / 2)
        # region in zeta: zeta <= min(ZMAX, ETAMAX / V); V-band -> zeta <= min(ZMAX, ETAMAX/lo)
        zhi = min(float(ZMAX), float(ETAMAX) / lo)
        zb = arb(zhi / 2, zhi / 2)
        # Edgeworth remainder: worst case at the largest eta on the band: t_min = 2 khat_1_min /(eta_max)
        eta_max = min(float(ETAMAX), zhi * hi)
        # bounds on khat_j over the band (interval evaluation of the series)
        kh = [khat[j].evaluate(Vb, zb) for j in range(JCUM)]
        k1lo = kh[0].lower()
        k2lo = kh[1].lower()
        assert k1lo > 0 and k2lo > 0
        t_min = k1lo / eta_max
        half_eps = arb(eta_max * hi / 2, eta_max * hi / 2)  # 1/N = eps = eta v in [0, eta_max v_hi]
        k2t = arb(
            (kh[1].lower() * lo + kh[1].upper() * hi) / 2,
            (kh[1].upper() * hi - kh[1].lower() * lo) / 2,
        )
        sigma_lo = (k2lo / eta_max).sqrt()  # sigma^2 = N kt_2 = khat_2/eta
        cvals = {1: arb(0, (1 / sigma_lo).upper())}
        cabs = {1: arb((1 / sigma_lo).upper())}
        for j in range(3, DA.JMAX + 1):
            # |c_j| <= eta^{j/2-1} |khat_j| / (j! khat_2^{j/2}) <= eta_max^{j/2-1} sup|khat_j| / (j! k2lo^{j/2})
            supk = abs(kh[j - 1]).upper()
            bnd = (
                arb(eta_max) ** (arb(j) / 2 - 1) * supk / (factorial(j) * k2lo ** (arb(j) / 2))
            ).upper()
            cvals[j] = arb(0, bnd)
            cabs[j] = arb(bnd)
        best = None
        for Lc in (6, 8, 10, 12, 15, 20, 30, 50):
            if Lc > sigma_lo.lower():
                break
            try:
                _, M0c = DA.kappa2_enclosure(cvals, cabs, sigma_lo, half_eps, k2t, arb(Lc))
            except Exception:  # noqa: BLE001
                continue
            dj = DA.kappa2_enclosure.last_delta
            tot = dj[0] + dj[1] + dj[2]
            if best is None or tot < best:
                best = tot
        # remainder in kappa_2[u]: the delta_j enter kappa_2 = M2/M0 + A1^2/M0^2 with M0 ~ 1: |d kappa_2| <= 3 (d0 + d1 + d2)
        R_edge = 3 * best
        # effect on S: R_edge/(kt_2 V) <= R_edge / (k2lo lo^2)
        S_loss = R_edge / (k2lo * lo * lo) * FAC  # times |3 + 3 eps - theta| <= FAC
        # binomial cut: (1/x)^6/(7 theta V) with 1/x <= zeta V/(2 k1lo): <= (zhi hi/(2 k1lo))^6 /(7 k1lo lo^2)
        r1 = (
            (arb(zhi * hi) / k1lo) ** 6 / (7 * k1lo * lo * lo)
            + (arb(zhi * hi * hi) / (1 - 0.11)) ** 4 / (5 * 0.89 * lo)
        ) * FAC
        # the series S itself on the band: subdivide in (V, zeta) until S - losses > 0 on every sub-box
        subtract = S_loss.upper() + r1.upper()
        stack = [(lo, hi, 0.0, zhi, 0)]
        nsub = 0
        band_min = None
        failed = False
        while stack:
            v1, v2, z1, z2, dep = stack.pop()
            nsub += 1
            Sval = S.evaluate(arb((v1 + v2) / 2, (v2 - v1) / 2), arb((z1 + z2) / 2, (z2 - z1) / 2))
            low = Sval.lower() - subtract
            if low > 0:
                band_min = low if band_min is None else min(band_min, low)
                continue
            if dep > 30:
                failed = True
                print("    FAILED sub-box", v1, v2, z1, z2, Sval)
                break
            if (v2 - v1) / v2 > (z2 - z1) / max(z2, 1e-9):
                vm = (v1 + v2) / 2
                stack.append((v1, vm, z1, z2, dep + 1))
                stack.append((vm, v2, z1, z2, dep + 1))
            else:
                zm = (z1 + z2) / 2
                stack.append((v1, v2, z1, zm, dep + 1))
                stack.append((v1, v2, zm, z2, dep + 1))
        low = band_min if not failed else arb(-1)
        print(
            f"  band V in [{lo:.3e}, {hi:.3e}], zeta <= {zhi:.3f}, t >= {float(t_min):.0f}, sigma >= {float(sigma_lo):.1f}: "
            f"{nsub} sub-boxes, edge loss {float(S_loss):.2e}, cut {float(r1):.1e}  -> min(S - losses) = {float(low):+.4f}"
        )
        if worst_S is None or low < worst_S:
            worst_S = low
    # ---- the corner V -> 0 (all V below the last band) ---------------------------------------------
    # With the inner cut L = sigma/2, every inner remainder piece is a sum of terms C_w sigma^{-w} with w >= 11
    # (|c_j| <= K_j sigma^{2-j} for j <= 13 from the series, c_1 <= 1/sigma, the Gaussian moments are sigma-free,
    # beta and e^{c_1 L} are sigma-free), and the outer pieces are bounded by explicit functions decreasing faster
    # than any power of sigma.  For j >= 14 the sigma-free K_j comes from ONE Bernoulli cumulant (3 Sept 2026,
    # after the validator showed the N/sigma^14 form is not sigma-free): Cauchy on |s| = rho = 3/2 for
    # log(1 + q(e^s - 1)) with q <= V/(1+V) <= 0.1453 gives |kappa_j(q)| <= j! rho^{-j} q (e^rho - 1)/(1 - q(e^rho - 1))
    # <= 7.05 j! (2/3)^j q; summing over the roots (SUM 2 q_k = t = k_1) and dividing by j! sigma^j,
    #     |c_j| <= 7.05 (2/3)^j (t/sigma^2) sigma^{2-j} = 7.05 (2/3)^j (khat_1/khat_2) sigma^{2-j},
    # so SUM_{j>=14} |c_j| |u|^j <= Kg sigma^2 SUM_{j>=14} (|u|/(rho sigma))^j with Kg = 7.05 khat_1/khat_2, rho = 3/2.
    # of sigma.  Hence for sigma >= sigma_last:  R_inner(sigma) <= R_inner(sigma_last) (sigma_last/sigma)^11,
    # and with sigma^2 >= 2 khat_2/(zeta V) >= (2 k2lo/ZMAX)/V the loss R/(k2 V^2) on (0, V_last] is at most
    #   loss_inner_last (V/V_last)^{3.5} + R_outer(sigma(V))/(k2lo V^2),
    # both decreasing as V -> 0 (the second because e^{-c/V}/V^2 decreases for V < c/2).  So the corner is
    # covered by the value at V_last.  Numbers with L = sigma_last/2:
    lo, hi = bands[-1]
    Vb = arb(
        hi / 2, hi / 2
    )  # the whole corner (0, V_last], not only the last band (validator, 3 Sept 2026)
    zb = arb(float(ZMAX) / 2, float(ZMAX) / 2)
    kh = [khat[j].evaluate(Vb, zb) for j in range(JCUM)]
    k2lo = kh[1].lower()
    eta_max = min(float(ETAMAX), float(ZMAX) * hi)
    sigma_lo = (k2lo / eta_max).sqrt()
    half_eps = arb(eta_max * hi / 2, eta_max * hi / 2)
    k2t = arb(
        (kh[1].lower() * lo + kh[1].upper() * hi) / 2, (kh[1].upper() * hi - kh[1].lower() * lo) / 2
    )
    cvals = {1: arb(0, (1 / sigma_lo).upper())}
    cabs = {1: arb((1 / sigma_lo).upper())}
    for j in range(3, DA.JMAX + 1):
        bnd = (
            arb(eta_max) ** (arb(j) / 2 - 1)
            * abs(kh[j - 1]).upper()
            / (factorial(j) * k2lo ** (arb(j) / 2))
        ).upper()
        cvals[j] = arb(0, bnd)
        cabs[j] = arb(bnd)
    Lc = arb(sigma_lo.lower() / 2)
    # q_k <= v/(1+v) <= VMAX/(1+VMAX); Cauchy on |s| = 3/2: |kappa_j(q)| <= j! (2/3)^j q (e^{3/2}-1)/(1 - q (e^{3/2}-1))
    qmax = toarb(VMAX) / (1 + toarb(VMAX))
    e32 = (arb(3) / 2).exp() - 1
    assert (qmax * e32).upper() < 1
    Kfac = (e32 / (1 - qmax * e32)).upper()
    Kg = arb(Kfac) * kh[0].upper() / k2lo
    DA.kappa2_enclosure(cvals, cabs, sigma_lo, half_eps, k2t, Lc, geo_K=(Kg, arb(3) / 2))
    dj = DA.kappa2_enclosure.last_delta
    R_last = 3 * (dj[0] + dj[1] + dj[2])
    loss_corner = R_last / (k2lo * lo * lo) * FAC
    print(
        f"  corner V <= {hi:.3e}: with L = sigma/2 = {float(Lc):.1f}, remainder {float(R_last):.2e}, loss <= {float(loss_corner):.2e}"
        f"  (S >= {float(S.evaluate(arb(hi / 2, hi / 2), zb).lower()):.4f} there; the corner is covered by monotonicity, see text)"
    )
    ok = worst_S is not None and worst_S > 0 and float(loss_corner) < 0.1
    print()
    print(
        f"VERDICT: (3 + 3/N - theta) N g > 1 on the dense region (b') {{ v <= {VMAX}, zeta <= {ZMAX}, zeta v <= {ETAMAX} }}: {ok}   [{time.time() - t0:.0f} s]"
    )
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
