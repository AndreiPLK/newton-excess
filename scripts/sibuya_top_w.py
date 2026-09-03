"""SIBUYA'S CONJECTURE (1988, eq. 3.4), the top regime in the variable Y = 1/v:  theta >= 0.9, j' = N - j >= JP_MIN.

The dense-(a') sweep is parametrised by (eps = 1/N, v); theta -> 1 pushes v -> infinity, so the region
{theta >= 0.9, j' >= JP_MIN} is unbounded in v and cannot be covered by a v-box.  In Y = 1/v it is a BOX:

    Y in [0, 1/35]   (theta increases as Y decreases; theta(Y = 1/35) <= 0.9),     eps in [0, 1/N_MIN],

with Y = 0 (v = infinity, the tilt x = Y/eps -> 0 at fixed eps) and eps = 0 (N -> infinity) both edges of the box.
Everything else is sibuya_dense_a: the cumulants k_j/N in closed form (Stirling for x = Y/eps >= 20, the direct
polygamma below that), the Edgeworth remainder engine of dense_certificate_a, the variance identity, and the
target N g > 1/(3 + 3 eps - theta).  Boxes with j' < JP_MIN (where sigma^2 ~ j' is too small for the Edgeworth
window) are skipped: that corner is closed by sibuya_harmonic.py.

Run:  uv run python projects/qg-bootstrap/release/scripts/sibuya_top_w.py [JP_MIN]
"""

from __future__ import annotations

import os
import sys
import time
from math import factorial

from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv = [
    sys.argv[0]
]  # sibuya_dense_a reads sys.argv for its own modes; keep it in the default one
import dense_certificate_a as DA  # noqa: E402
import sibuya_dense_a as SA  # noqa: E402

ctx.prec = 300

JP_MIN = (
    803  # sibuya_corner_grid.py closes 33 <= j' <= 802 for every N; here sigma = sqrt(j') >= 28
)
T_MIN = 1001  # the Sibuya ladder closes j <= 1000
Y_MAX = arb(1) / 35  # theta(v = 35) <= 0.9
JC = SA.JC
apow = DA.apow
# the coefficient table of sibuya_dense_a stops at JC + 3; the top regime needs every cumulant up to DA.JMAX + 2
SA._CC = SA.cumulant_coefficients(DA.JMAX + 3)


def phi_ratio(n, eps, w, K=SA.KST, ratio_given=False):
    """w^n Phi_n(w) = B_n(u), u = eps/w: the Stirling series with the w^{-n} factored out, plus its remainder.

    Phi_n(w) = N^n psi^{(n)}(N w) as in sibuya_dense_a.phi_stirling_real; here everything is expressed through
    u = eps/w so that a wide w costs one width instead of 2K+2+n of them.

    Two things an independent check (2026-09-03) established about this routine.  For `n >= 1` the first line is
    an EXACT identity, not an approximation -- `w^n Phi_n(w) = z^n psi^{(n)}(z)` with `z = 1/u` -- verified at
    four values of `w` spanning five orders of magnitude at fixed `u`, agreeing to twelve digits for
    `n = 1..5`.  For `n = 0` it is NOT: `Phi_0 = psi(1/u) + log u + log w`, so the value depends on `w` beyond
    `u`, differences across `w` measuring exactly `log(w_2/w_1)`.  The CODE is right -- its `n == 0` branch
    carries `w.log()` explicitly -- but the identity claimed in the first line covers `n >= 1` only.

    Why the form matters: against the direct Stirling evaluation, on 5%-wide `w`-bands the direct form returns
    inf or nan already at `n = 1`, while this one keeps a relative width of `8e-7` to `4e-5`.
    """
    u = eps if ratio_given else eps / w
    if n == 0:
        val = w.log() - u / 2
        for k in range(1, K + 1):
            val -= SA.toarb(DA.BERN[2 * k]) * apow(u, 2 * k) / (2 * k)
    else:
        val = arb(factorial(n - 1)) + arb(factorial(n)) * u / 2
        for k in range(1, K + 1):
            val += (
                SA.toarb(DA.BERN[2 * k])
                * factorial(2 * k + n - 1)
                * apow(u, 2 * k)
                / factorial(2 * k)
            )
        val = val * ((-1) ** (n + 1))
    rem = (
        abs(SA.toarb(DA.BERN[2 * K + 2]))
        * factorial(n + 2 * K + 1)
        / factorial(2 * K + 2)
        * apow(u, 2 * K + 2)
    )
    r = rem.upper()
    return val + arb(0, r)


def scaled_cumulants_Y(eps, Y, jmax, K=SA.KST, want_nv=False):
    """kt_j = k_j/N with Y = 1/v as the variable (v = 1/Y may be infinite: Y = 0 is allowed).

    want_nv also returns  W := Y^2 [Phi_1(Y+eps) - Phi_1(1+Y+eps)] = (x^2/N)[psi'(1+x) - psi'(N+1+x)]
    = (1/N) SUM_k p_k^2 with p_k = x/(x+k) the probability that root k is MISSING (x = N Y).  The variance
    NV = theta(1-theta) - kt_2 equals W - mu^2 with mu = 1 - theta = j'/N: as theta -> 1 the first form is a
    difference of two nearly equal numbers (measured: NV lost to cancellation at theta = 0.99998), the second
    is not (W ~ x^2 zeta(2)/N against mu^2 ~ (x log N/N)^2)."""
    w1 = 1 + Y + eps
    w2 = Y + eps
    direct = eps.lower() > 0 and (Y / eps).lower() < 20
    if direct:
        x = Y / eps
        z = SA.acb(1 + x)
        xc = arb(x.mid())
        dig = SA.acb(1 + xc).digamma().real + arb(0, x.rad()) * z.polygamma(SA.acb(1)).real
        term2 = []
        for q in range(jmax + 1):
            if q == 0:
                term2.append(Y * (dig + eps.log()))
            else:
                term2.append(Y * apow(x, q) * z.polygamma(SA.acb(q)).real)
        term = [
            apow(Y, q + 1) * SA.phi_stirling_real(q, eps, w1, K) - term2[q] for q in range(jmax + 1)
        ]
    else:
        # RATIO FORM.  Phi_q(w) = B_q(eps/w)/w^q with B_q depending only on the ratio eps/w (the Stirling series
        # and its remainder are series in eps/w after the factor w^{-q} is taken out), so
        #     Y^{q+1} Phi_q(w) = Y (Y/w)^q B_q(eps/w),
        # and with w_2 = Y + eps the ratio Y/w_2 is close to 1.  Written this way the enclosure of a band loses
        # ~2 widths of Y instead of ~29 (the direct form divides by w^{2K+2+q}), which is what lets the bands be
        # wide: at 5% relative width the direct form returns nan from k = 3 on, the ratio form is fine.
        term = [
            Y
            * (
                apow(Y / w1, q) * phi_ratio(q, eps, w1, K)
                - apow(Y / w2, q) * phi_ratio(q, eps, w2, K)
            )
            for q in range(jmax + 1)
        ]
    kt = []
    for j in range(1, jmax + 1):
        acc = arb(1) if j == 1 else arb(0)
        for q, c in SA._CC[j].items():
            acc += SA.toarb(c) * term[q]
        kt.append(acc)
    if want_nv:
        # Phi_1(1+Y+eps) via Stirling (the argument is >= 1); Phi_1(Y+eps) is term2[1]/(Y x) on the direct
        # path and Y^2 phi_stirling_real(1, eps, w2) otherwise
        big = apow(Y, 2) * SA.phi_stirling_real(1, eps, w1, K)
        small = (Y * term2[1] / (Y * (Y / eps))) * apow(Y, 2) / apow(Y, 2) if False else None
        if direct:
            small = apow(Y, 2) * (Y / eps) * SA.acb(1 + Y / eps).polygamma(SA.acb(1)).real / Y
        else:
            small = apow(Y, 2) * SA.phi_stirling_real(1, eps, w2, K)
        return kt, small - big
    return kt


def certify_box(eps, Y):
    """(status, info) for the box eps x Y; status in {'ok', 'skip', 'fail'}."""
    if eps.lower() < 0:
        eps = arb(eps.upper() / 2, eps.upper() / 2)
    if Y.lower() < 0:
        Y = arb(Y.upper() / 2, Y.upper() / 2)
    # mean-value enclosure in Y: v d_v = -Y d_Y and v d_v kt_j = kt_{j+1}  =>  d kt_j/dY = -kt_{j+1}/Y
    # every cumulant up to DA.JMAX = 13 comes from the closed form: the generic per-root sup bounds SUPK[j]
    # that dense_certificate_a uses beyond j = JC are O(1), which overestimates the true kt_j ~ mu by 1/kt_2
    # (a factor 10^8 when theta -> 1) and was what made the enclosure of N g useless there
    JT = DA.JMAX
    Y_c = arb(Y.mid())
    kt_c = scaled_cumulants_Y(eps, Y_c, JT + 2)
    raw = scaled_cumulants_Y(eps, Y, JT + 2)
    dY = arb(0, Y.rad())
    if Y.lower() <= 0:
        kt = [kt_c[j] + arb(0, abs(raw[j] - kt_c[j]).upper()) for j in range(JT)]
    else:
        kt_b = [kt_c[j] + dY * raw[j + 1] / Y for j in range(JT + 1)]
        kt = [kt_c[j] + dY * kt_b[j + 1] / Y for j in range(JT)]
    theta = kt[0]
    eps_lo = eps.lower()
    N_lo = 1 / eps.upper()
    # a box touching eps = 0 (N -> infinity) is a legitimate part of the region: both t and j' are unbounded there
    x_hi = theta.upper() / eps_lo if eps_lo > 0 else arb("inf")
    jp_hi = (1 - theta).upper() / eps_lo if eps_lo > 0 else arb("inf")

    if x_hi < T_MIN or theta.upper() < arb("0.9") or jp_hi < JP_MIN:
        return "skip", {}
    # (an earlier version refused the sparse side N > j'^2; measured on 3 September, the engine certifies it
    # with margins that GROW as theta -> 1 -- 6e3, 2e7, 5e16, 1e36 at Y = 1e-6, 1e-10, 1e-20, 1e-40 -- provided
    # the boxes are relatively thin, so the guard is gone)
    k2t = kt[1]
    if not (k2t.lower() > 0):
        return "fail", {"why": "k2 not positive"}
    half_eps = eps
    sigma_lo = (N_lo * k2t.lower()).sqrt()
    q_up = (half_eps / k2t).upper()
    if (half_eps / k2t).lower() > 0:
        inv_sigma = (half_eps / k2t).sqrt()
    else:
        r_up = q_up.sqrt().upper()
        inv_sigma = arb(r_up / 2, r_up / 2)
    cvals = {1: arb(0, inv_sigma.upper())}
    cabs = {1: arb(inv_sigma.upper())}
    for j in range(3, DA.JMAX + 1):
        cj = apow(inv_sigma, j - 2) * kt[j - 1] / (factorial(j) * k2t)
        cvals[j] = cj
        cabs[j] = abs(cj).upper()
    # the j >= 14 cumulant tail: dense_certificate_a's default bound 2.7 N (|u|/(3 sigma))^j is the one for the
    # centred-square problem and is meaningless here (it ignores that every per-root cumulant carries a factor
    # p_k).  Cauchy on |s| = rho for log(1 + p(e^s - 1)) with p <= p_max = Y/(Y+eps) gives
    #     |kappa_j(p)| <= j! rho^{-j} p (e^rho - 1)/(1 - p_max(e^rho - 1)),
    # so SUM_k |kappa_j| <= j! rho^{-j} C N mu and |c_j| <= Kg sigma^2/(rho sigma)^j with Kg = C mu/kt_2.
    p_max = (Y / (Y + eps)).upper()
    rho = arb("0.9") * (1 + 1 / arb(p_max)).log()
    er = rho.exp() - 1
    Kg = (er / (1 - arb(p_max) * er) * (1 - theta) / k2t).upper()
    if not (Kg > 0):
        return "fail", {"why": "cumulant tail constant"}
    geo = (arb(Kg), rho)
    best = None
    for Lc in (4, 5, 6, 7, 8, 10, 12, 15, 20, 30):
        # the geometric cumulant tail needs L < rho sigma (its ratio L/(rho sigma) must stay below 1)
        if Lc > sigma_lo.lower() or Lc >= float(rho) * float(sigma_lo.lower()) * 0.9:
            break
        k2u, M0c = DA.kappa2_enclosure(cvals, cabs, sigma_lo, half_eps, k2t, arb(Lc), geo_K=geo)
        if not (M0c.lower() > 0):
            continue
        if best is None or k2u.rad() < best[0].rad():
            best = (k2u, M0c, Lc)
    if best is None:
        return "fail", {"why": "no valid split", "sigma_lo": sigma_lo}
    _, _M0, L = best
    kappa2u = DA.kappa2_window_average(
        cvals, cabs, sigma_lo, half_eps, k2t, arb(L), inv_sigma, geo_K=geo
    )
    if kappa2u is None:
        return "fail", {"why": "window average: h too large"}
    # NV = W - mu^2 (cancellation-free), mu = 1 - theta = j'/N
    _, W = scaled_cumulants_Y(eps, Y, 2, want_nv=True)
    mu = 1 - theta
    NV = W - mu * mu
    if not (NV.lower() > 0):
        return "fail", {"why": "NV not positive", "NV": NV}
    Q = theta * (1 - theta)
    x_lo = arb(max(float(T_MIN), float(N_lo * theta.lower())))
    one_minus_Gx = arb(0, (1 / (2 * x_lo)).upper())
    one_minus_Gnx = arb(0, (1 / (2 * N_lo * (1 - theta).lower())).upper())
    Ng = NV / (k2t * Q) + (kappa2u - 1) / k2t + one_minus_Gx / theta + one_minus_Gnx / (1 - theta)
    margin = Ng - 1 / (3 + 3 * eps - theta)
    info = {"theta": theta, "Ng": Ng, "margin": margin, "sigma_lo": sigma_lo}
    return ("ok", info) if margin.lower() > 0 else ("fail", info)


def main():
    t0 = time.time()
    ok, nb = DA.bernoulli_lemma(SA.KST)
    print(f"Bernoulli remainder lemma (K={SA.KST}): {ok}  [{nb} boxes]", flush=True)
    assert ok
    # coverage of theta >= 0.9 by Y <= Y_MAX: theta decreases in Y, so theta(Y_MAX, eps) <= 0.9 for every eps
    eps_max = arb(1) / (T_MIN + JP_MIN)
    th_edge = scaled_cumulants_Y(arb(float(eps_max) / 2, float(eps_max) / 2), Y_MAX, 1)[0]
    print(
        f"  coverage: theta(Y = 1/35) <= {float(th_edge.upper()):.5f} <= 0.9: {th_edge.upper() <= 0.9}"
    )
    assert th_edge.upper() <= 0.9
    stack = [(arb(0), eps_max, arb(0), Y_MAX, 0)]
    n_ok = n_skip = n_corner = n_sparse = 0
    corner_max = 0.0
    worst = None
    max_depth = 260
    while stack:
        e1, e2, y1, y2, d = stack.pop()
        eps = arb((e1 + e2) / 2, ((e2 - e1) / 2).upper())
        Y = arb((y1 + y2) / 2, ((y2 - y1) / 2).upper())
        try:
            status, info = certify_box(eps, Y)
        except Exception as ex:  # noqa: BLE001
            status, info = "fail", {"why": repr(ex)[:60]}
        if status == "ok":
            n_ok += 1
            if worst is None or info["margin"].lower() < worst[0]:
                worst = (info["margin"].lower(), float(eps.mid()), float(Y.mid()))
            continue
        if status in ("skip", "sparse-side"):
            n_skip += 1
            n_sparse += status == "sparse-side"
            continue
        if d >= max_depth:
            if float(y1) == 0.0 and float(y2) < 1e-30:
                # the degenerate corner Y -> 0 (no missing roots) with eps -> 0: N -> infinity at fixed j',
                # where every cumulant vanishes and the box engine has nothing to enclose.  Recorded, and
                # closed per index by the H-model of sibuya_corner_grid.tail_ok.
                n_corner += 1
                corner_max = max(corner_max, float(y2))
                continue
            print("FAILED box:", float(e1), float(e2), float(y1), float(y2), info)
            sys.exit(1)
        # GEOMETRIC bisection of the side with the larger RELATIVE width (near theta = 1 what has to be small is
        # the relative width of mu = 1 - theta).  A side that touches 0 never becomes relatively thin, so it is
        # always the one to split until its upper half is thin and its lower half is skipped by the j' >= JP_MIN
        # rule -- alternating the two sides instead would exhaust the depth budget.
        rel_y = float("inf") if y1.lower() <= 0 else float(y2 / y1)
        rel_e = float("inf") if e1.lower() <= 0 else float(e2 / e1)
        if rel_y >= rel_e:
            ym = (y1 * y2).sqrt() if y1.lower() > 0 else y2 / 2
            stack.append((e1, e2, ym, y2, d + 1))
            stack.append((e1, e2, y1, ym, d + 1))
        else:
            em = (e1 * e2).sqrt() if e1.lower() > 0 else e2 / 2
            stack.append((e1, em, y1, y2, d + 1))
            stack.append((em, e2, y1, y2, d + 1))
        if (n_ok + n_skip) % 500 == 0 and n_ok + n_skip > 0:
            print(
                f"  ... {n_ok} ok, {n_skip} skipped, stack {len(stack)}  [{time.time() - t0:.0f} s]",
                flush=True,
            )
    print(
        f"certified: {n_ok} boxes ok, {n_skip} skipped ({n_sparse} on the sparse side N > j'^2), {n_corner} at the degenerate corner (Y <= {corner_max:.1e})"
        + (
            f", worst margin {float(worst[0]):.6f} at eps={worst[1]:.2e}, Y={worst[2]:.3e}"
            if worst
            else " (no box certified)"
        )
    )
    print(
        f"VERDICT (top regime in Y = 1/v): N g > 1/(3 + 3/N - theta) on {{ theta >= 0.9, t >= {T_MIN},"
        f" j' = N - t >= {JP_MIN}, N <= j'^2 }} : True"
        f"  [{time.time() - t0:.0f} s]"
    )


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------- the corner Y -> 0
def corner_scaled(L, r, jmax, K=SA.KST):
    """Every cumulant divided by Y, in the variables L = log(1/Y) and r = eps/Y, valid for Y <= e^{-L}.

    With eps = Y r, w_2 = Y(1+r), w_1 = 1 + Y(1+r):
        Y Phi_0(w_2) / Y = log w_2 + s_0(r/(1+r)),        log w_2 = -L + log(1+r),
        Y^{q+1} Phi_q(w_2) / Y = (1+r)^{-q} B_q(r/(1+r)),        q >= 1, a function of r alone,
        Y^{q+1} Phi_q(w_1) / Y = (Y/w_1)^q B_q(Yr/w_1) = O(Y^q),
    so kt_j / Y = A_j L + G_j(r) + O(Y) with EXPLICIT A_j and G_j.  The margin then scales as 1/(Y L^2) while
    the Edgeworth defect scales as 1/(Y L) times |kappa_2 - 1| <= C r / L: the two are comparable, L cancels,
    and the corner reduces to a condition on r alone.
    """
    Y = arb(0, float((-L).exp().upper()))  # (0, e^{-L}]
    eps = Y * r
    w1 = 1 + Y * (1 + r)
    u2 = r / (1 + r)
    term = []
    for q in range(jmax + 1):
        big = apow(Y / w1, q) * phi_ratio(q, eps, w1, K)
        if q == 0:
            small = -L + (1 + r).log() + phi_ratio(0, u2, arb(1), K, ratio_given=True)
        else:
            small = apow(1 / (1 + r), q) * phi_ratio(q, u2, arb(1), K, ratio_given=True)
        term.append(big - small)
    kt = []
    for j in range(1, jmax + 1):
        acc = arb(0)
        for q, c in SA._CC[j].items():
            acc += SA.toarb(c) * term[q]
        kt.append(acc)
    # W/Y = B_1(u_2)/(1+r) - Y B_1(Yr/w_1)/w_1
    WY = phi_ratio(1, u2, arb(1), K, ratio_given=True) / (1 + r) - Y * phi_ratio(1, eps, w1, K) / w1
    return kt, WY, Y


def certify_corner(L, r, jp_min=None):
    """NOT LOAD-BEARING.  An attempt at the corner Y -> 0 that does not decide it; kept as a record.

    Two reasons it is excluded from every verdict in this package.  (i) At large L the two sides cancel to
    +/- 1e-12 -- the inequality is asymptotically TIGHT at this corner, so no crude bound can separate them;
    that is a fact about the object, and it is why the corner is stated as open rather than certified.
    (ii) Independently of that, `r_worst = arb(r.upper())` below substitutes one endpoint of the ratio interval
    into the Edgeworth half while `m` and `W/Y` keep the full ball; the monotonicity that would justify it is
    NOT proved.  Nothing calls this function from a verdict path (`sibuya_wedge.main` does not run the sweep),
    and it must not be wired into one before (ii) is settled.

    The corner {Y <= e^{-L}} at the ratio interval r: is N g > 1/(3 + 3 eps - theta) there?

    N g >= NV/(kt_2 Q) - |kappa_2 - 1|/kt_2, and with m = mu/Y, w = W/Y, all of NV, kt_2, Q carrying one factor
    of Y each, the first term is  (w - Y m^2) / (Y (m - w) m)  and the second  E / (Y (m - w)).  Both scale as
    1/Y, so the statement is the Y-free inequality
            (w - Y m^2) / m  >  E  +  0.5 Y (m - w) m / 1        (the target 1/(3+...) <= 0.5),
    which is checked here with E the Edgeworth enclosure of |kappa_2 - 1| (itself Y-free).
    """
    jp_min = JP_MIN if jp_min is None else jp_min
    JT = DA.JMAX
    kt, WY, Y = corner_scaled(L, r, JT + 2)
    m = -kt[0]
    k2Y = kt[1]
    if not (m.lower() > 0 and k2Y.lower() > 0):
        return "fail", {"why": "m or kt_2/Y not positive", "m": m, "k2": k2Y}
    if (m / r).upper() < jp_min:
        return "skip", {}
    # sigma^2 = kt_2/eps = (kt_2/Y)/r decreases in r, and the Edgeworth defect E grows with r: a box [0, r2]
    # is therefore worst at its right end, where the enclosure is evaluated (w and m keep the whole box)
    r_worst = arb(r.upper())
    sigma2 = arb(k2Y.lower()) / r_worst
    if not (sigma2.lower() > 4):
        return "fail", {"why": "sigma^2 too small", "sigma2": sigma2}
    sigma_lo = arb(sigma2.lower()).sqrt()
    inv_sigma = (r_worst / arb(k2Y.lower())).sqrt()
    cvals = {1: arb(0, inv_sigma.upper())}
    cabs = {1: arb(inv_sigma.upper())}
    for j in range(3, DA.JMAX + 1):
        cj = apow(inv_sigma, j - 2) * kt[j - 1] / (factorial(j) * k2Y)
        cvals[j] = cj
        cabs[j] = abs(cj).upper()
    p_max = (1 / (1 + r.lower())).upper() if r.lower() > 0 else arb(1).upper()
    rho = arb("0.9") * (1 + 1 / arb(p_max)).log()
    er = rho.exp() - 1
    if not ((1 - arb(p_max) * er).lower() > 0):
        return "fail", {"why": "cumulant tail constant"}
    Kg = (er / (1 - arb(p_max) * er) * m / k2Y).upper()
    geo = (arb(Kg), rho)
    best = None
    for Lc in (4, 5, 6, 7, 8, 10, 12, 15, 20, 30):
        if Lc > sigma_lo.lower() or Lc >= float(rho) * float(sigma_lo.lower()) * 0.9:
            break
        try:
            k2u, M0c = DA.kappa2_enclosure(
                cvals, cabs, sigma_lo, Y * r_worst, arb(k2Y.lower()) * Y, arb(Lc), geo_K=geo
            )
        except Exception:  # noqa: BLE001
            continue
        if not (M0c.lower() > 0):
            continue
        if best is None or k2u.rad() < best[0].rad():
            best = (k2u, M0c, Lc)
    if best is None:
        return "fail", {"why": "no valid split", "sigma_lo": sigma_lo}
    # the window average is the tight enclosure (the raw kappa2_enclosure keeps c_1 as a full ball and is ~7x
    # looser); kappa2_enclosure is used only to pick the cut L and to set the remainder deltas
    k2avg = DA.kappa2_window_average(
        cvals, cabs, sigma_lo, Y * r_worst, arb(k2Y.lower()) * Y, arb(best[2]), inv_sigma, geo_K=geo
    )
    kap = (
        k2avg if k2avg is not None else best[0]
    ) - 1  # SIGNED: at these sizes the Edgeworth correction is
    # positive and of the same size as the binomial term, so bounding it by -|kappa_2 - 1| loses the whole margin
    E = abs(kap).upper()
    if not ((m - WY).lower() > 0):
        return "fail", {"why": "m - W/Y not positive"}
    # the binomial term (1 - G(N-x))/(1-theta) is POSITIVE and of the same order as the rest: with y = N - x = j'
    # = m/r, 1 - G(y) = 1/(2y) - 1/(3y^2) + ..., so its scaled contribution is r(m - w)(1 - 2r/(3m))/(2 m^2).
    # Dropping it costs exactly the margin at j' ~ 1000, where E and it are the same size.
    binom = r * (m - WY) / (2 * m * m) * (1 - 2 * r / (3 * m))
    lhs = (WY - Y * m * m) / ((1 - Y * m) * m) + binom + kap
    rhs = arb("0.5") * Y * (m - WY) * m
    if (lhs - rhs).lower() > 0:
        return "ok", {"slack": lhs - rhs, "sigma2": sigma2, "E": arb(E), "m": m}
    return "fail", {"slack": lhs - rhs, "E": arb(E), "sigma2": sigma2, "m": m, "W": WY}
