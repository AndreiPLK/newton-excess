"""SIBUYA'S CONJECTURE (1988, eq. 3.4), the dense regime (a'):  0.1 <= theta <= 1/2,  t >= 2001,  N >= 4002.

Spectrum 1, 2, ..., N (unsigned Stirling numbers of the first kind, N roots, n = N + 1).  Statement certified:

    N g  >  N / (3N + 3 - t)  =  1 / (3 + 3 eps - theta),      eps = 1/N,  theta = t/N,

which gives  n (R_t - 1) > n g > n/(3n - t), i.e. Sibuya's  p_t^2/(p_{t-1} p_{t+1}) >= 1 + 1/(3n - t).

The machine is dense_certificate_a.py (the 4/5 theorem) with the ONLY spectrum-dependent part replaced:
the tilted cumulants.  With the tilt r = v/N,  q_k = k r/(1 + k r)  and  x := 1/r = N/v,

    k_1 = SUM_k q_k = N - x [ psi(N + 1 + x) - psi(1 + x) ],        k_{j+1} = v d_v k_j   (fixed N),

so with  Y = 1/v,  w_1 = 1 + 1/v + eps,  w_2 = 1/v + eps,  Phi_q(w) := N^q psi^{(q)}(N w)  (Stirling in eps
with an explicit remainder, real arguments) and  Delta_q := Phi_q(w_1) - Phi_q(w_2):

    kt_j := k_j / N = [j == 1] + SUM_q C_{j,q} Y^{q+1} Delta_q,     C_{1,0} = -1,
    C_{j+1,q} = -(q+1) C_{j,q} - C_{j,q-1}          (from v d_v Y = -Y,  v d_v Delta_q = -Y Delta_{q+1}).

Everything else -- the exact Edgeworth polynomial, the window average, the remainder engine
kappa2_enclosure, the variance identity  N g = NV/(kt_2 Q) + (kappa_2 - 1)/kt_2 + (1 - G(x))/theta
+ (1 - G(N-x))/(1 - theta)  with  NV = theta(1 - theta) - kt_2 = Var(q) -- is imported unchanged.

Run:  uv run python projects/qg-bootstrap/release/scripts/sibuya_dense_a.py [--selftest] [--top]
      --top certifies the upper half 1/2 <= theta <= 9/10 (v in [2.5, 45], N >= 2224) instead of [1/10, 1/2].
"""

from __future__ import annotations

import os
import sys
import time
from math import factorial

from flint import acb, arb, ctx, fmpq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dense_certificate_a as DA  # noqa: E402

ctx.prec = 300

T_MIN = 1001  # the Sibuya ladder closes j <= 1000
N_MIN = 2 * T_MIN  # t >= 2001 and theta <= 1/2
if "--top" in sys.argv:
    N_MIN = (10 * T_MIN) // 9  # t >= 2001 and theta <= 9/10
THETA_LO = fmpq(1, 10)
THETA_HI = fmpq(1, 2)
if (
    "--top" in sys.argv
):  # the upper half: 1/2 <= theta <= 9/10 (Sibuya's (3.4) is stated for every j <= n-2)
    THETA_LO = fmpq(1, 2)
    THETA_HI = fmpq(9, 10)
JPRIME_MIN = 0  # --top2: boxes whose every N has fewer than JPRIME_MIN missing roots are skipped (the (T2) corner)
if (
    "--top2" in sys.argv
):  # EXPLORATORY: the top regime 9/10 <= theta <= 99/100 with j' = N - j >= JPRIME_MIN
    THETA_LO = fmpq(9, 10)
    THETA_HI = fmpq(99, 100)
    N_MIN = (100 * T_MIN) // 99
    _i = sys.argv.index("--top2")
    JPRIME_MIN = int(sys.argv[_i + 1]) if len(sys.argv) > _i + 1 else 100
KST = DA.KST
JC = DA.JC
apow = DA.apow
toarb = DA.toarb


# ------------------------------------------------------------------ real Stirling forms
def phi_stirling_real(n, eps, w, K=KST):
    """Phi_n(w) = N^n psi^{(n)}(N w) with N = 1/eps (n = 0: psi(N w) - log N), real w > 0, eps in [0, eps_hi];
    Stirling series to K terms with the remainder |B_{2K+2}| (n+2K+1)!/(2K+2)! eps^{2K+2} / w^{n+2K+2}
    (the classical bound for real positive arguments: the error is bounded by the first omitted term)."""
    if n == 0:
        val = w.log() - eps / (2 * w)
        for k in range(1, K + 1):
            val -= toarb(DA.BERN[2 * k]) * apow(eps, 2 * k) / (2 * k * apow(w, 2 * k))
    else:
        val = arb(factorial(n - 1)) / apow(w, n) + arb(factorial(n)) * eps / (2 * apow(w, n + 1))
        for k in range(1, K + 1):
            val += (
                toarb(DA.BERN[2 * k])
                * factorial(2 * k + n - 1)
                * apow(eps, 2 * k)
                / (factorial(2 * k) * apow(w, 2 * k + n))
            )
        val = val * ((-1) ** (n + 1))
    rem = (
        abs(toarb(DA.BERN[2 * K + 2]))
        * factorial(n + 2 * K + 1)
        / factorial(2 * K + 2)
        * apow(eps, 2 * K + 2)
        / apow(w, n + 2 * K + 2)
    )
    return val + arb(0, rem.upper())


def cumulant_coefficients(jmax):
    """C_{j,q} for j = 1..jmax."""
    C = {1: {0: fmpq(-1)}}
    for j in range(1, jmax):
        nxt = {}
        for q, c in C[j].items():
            nxt[q] = nxt.get(q, fmpq(0)) - (q + 1) * c
            nxt[q + 1] = nxt.get(q + 1, fmpq(0)) - c
        C[j + 1] = {q: c for q, c in nxt.items() if c != 0}
    return C


_CC = cumulant_coefficients(JC + 3)


def scaled_cumulants(eps, v, jmax, K=KST):
    """kt_j = k_j / N for j = 1..jmax as arb enclosures, eps = 1/N and v in boxes."""
    Y = 1 / v
    w1 = 1 + Y + eps
    w2 = Y + eps
    # Phi_q(w_2) = N^q psi^(q)(N w_2) with N w_2 = 1 + x, x = N/v = Y/eps.  The Stirling series in eps/w_2 = 1/(1+x)
    # is usable only for x >= 20; on boxes with eps bounded away from 0 and small x the polygamma is evaluated
    # DIRECTLY in arb (a certified enclosure on the ball 1 + x): Phi_q = eps^(-q) psi^(q)(1+x), Phi_0 = psi(1+x) + log eps.
    direct = eps.lower() > 0 and (Y / eps).lower() < 20
    if direct:
        # well-conditioned form: Y^(q+1) Phi_q(w_2) = Y x^q psi^(q)(1+x)  (q >= 1),  Y (psi(1+x) + log eps)  (q = 0),
        # so that eps^(-q) and x = Y/eps are never multiplied as independent balls
        x = Y / eps
        z = acb(1 + x)
        # arb's digamma of a wide ball is useless (radius 0.33 -> [+/- 21]); enclose it by the mean value theorem:
        # psi(1 + x) in psi(1 + x_c) + (x - x_c) psi'(1 + xi), xi in the ball (psi' of the ball is well enclosed)
        xc = arb(x.mid())
        dig = acb(1 + xc).digamma().real + arb(0, x.rad()) * z.polygamma(acb(1)).real
        term2 = []
        for q in range(jmax + 1):
            if q == 0:
                term2.append(Y * (dig + eps.log()))
            else:
                term2.append(Y * apow(x, q) * z.polygamma(acb(q)).real)
        term = [
            apow(Y, q + 1) * phi_stirling_real(q, eps, w1, K) - term2[q] for q in range(jmax + 1)
        ]
    else:
        term = [
            apow(Y, q + 1) * (phi_stirling_real(q, eps, w1, K) - phi_stirling_real(q, eps, w2, K))
            for q in range(jmax + 1)
        ]
    kt = []
    for j in range(1, jmax + 1):
        acc = arb(1) if j == 1 else arb(0)
        for q, c in _CC[j].items():
            acc += toarb(c) * term[q]
        kt.append(acc)
    return kt


def exact_cumulants(N, v, jmax):
    """Exact kt_j from the sums over k (for the self-test): kappa_j of independent Bernoulli(q_k) via the
    recursion kappa_{j+1}(q) = q(1-q) d_q kappa_j(q) evaluated as polynomials."""
    from flint import fmpq_poly

    kp = [None, fmpq_poly([0, 1])]
    q1q = fmpq_poly([0, 1, -1])
    for j in range(1, jmax):
        kp.append(q1q * kp[j].derivative())
    coeffs = [None] + [[toarb(c) for c in kp[j].coeffs()] for j in range(1, jmax + 1)]
    r = v / N
    tot = [arb(0)] * (jmax + 1)
    for k in range(1, N + 1):
        q = k * r / (1 + k * r)
        for j in range(1, jmax + 1):
            pv = arb(0)
            for c in reversed(coeffs[j]):
                pv = pv * q + c
            tot[j] += pv
    return [tot[j] / N for j in range(1, jmax + 1)]


def selftest():
    for N, v in ((5000, arb("0.7")), (4002, arb("2.5")), (20000, arb("0.22"))):
        eps = arb(1) / N
        kt = scaled_cumulants(eps, v, 6)
        ex = exact_cumulants(N, v, 6)
        mid = max(float(abs(arb(a.mid()) - arb(b.mid()))) for a, b in zip(kt, ex, strict=True))
        rad = max(float(a.rad()) for a in kt)
        contained = all(
            (a - b).lower() <= 0 <= (a - b).upper() for a, b in zip(kt, ex, strict=True)
        )
        print(
            f"  self-test N = {N}, v = {v}: max |mid(closed form) - mid(exact sum)| = {mid:.2e}, max radius {rad:.2e}, exact inside the enclosure: {contained}"
        )
        assert contained and mid < 1e-30, "closed form disagrees with the exact sums"
    print("  self-test PASS")


# ------------------------------------------------------------------ one box
def certify_box(eps, v, verbose=False):
    """Return (status, info): status in {'ok', 'skip', 'fail'} for the box eps x v (arb balls), eps = 1/N."""
    if eps.lower() < 0:
        eps = arb(eps.upper() / 2, eps.upper() / 2)
    # mean-value enclosure in v:  d kt_j / dv = (1/v) kt_{j+1}  (exact: v d_v k_j = k_{j+1} at fixed N)
    v_c = arb(v.mid())
    kt_c = scaled_cumulants(eps, v_c, JC + 2)
    raw = scaled_cumulants(eps, v, JC + 2)
    dv = arb(0, v.rad())
    kt_b = [kt_c[j] + dv * (1 / v) * raw[j + 1] for j in range(JC + 1)]
    kt = [kt_c[j] + dv * (1 / v) * kt_b[j + 1] for j in range(JC)]
    theta = kt[0]
    eps_lo = eps.lower()
    N_lo = 1 / eps.upper()
    x_hi = theta.upper() / eps_lo if eps_lo > 0 else arb("inf")
    if x_hi < T_MIN or theta.upper() < toarb(THETA_LO) or theta.lower() > toarb(THETA_HI):
        return "skip", {}
    if JPRIME_MIN and eps_lo > 0 and (1 - theta).upper() / eps_lo < JPRIME_MIN:
        return (
            "skip",
            {},
        )  # every N in the box has j' = (1 - theta) N < JPRIME_MIN: the (T2) corner, not this sweep
    k2t = kt[1]
    if not (k2t.lower() > 0):
        return "fail", {"why": "k2 not positive"}
    half_eps = eps  # = 1/N  (the name follows dense_certificate_a, where eps = 1/m = 2/N)
    sigma_lo = (N_lo * k2t.lower()).sqrt()
    q_up = (half_eps / k2t).upper()
    if (half_eps / k2t).lower() > 0:
        inv_sigma = (half_eps / k2t).sqrt()
    else:
        r_up = q_up.sqrt().upper()
        inv_sigma = arb(r_up / 2, r_up / 2)
    cvals = {1: arb(0, inv_sigma.upper())}
    cabs = {1: arb(inv_sigma.upper())}
    for j in range(3, JC + 1):
        cj = apow(inv_sigma, j - 2) * kt[j - 1] / (factorial(j) * k2t)
        cvals[j] = cj
        cabs[j] = abs(cj).upper()
    for j in range(JC + 1, DA.JMAX + 1):
        bnd = (apow(inv_sigma, j - 2) * DA.SUPK[j] / (factorial(j) * k2t)).upper()
        cvals[j] = arb(0, bnd)
        cabs[j] = arb(bnd)
    best = None
    for Lc in (4, 5, 6, 7, 8, 10, 12, 15, 20, 30):
        if Lc > sigma_lo.lower():
            break
        k2u, M0c = DA.kappa2_enclosure(cvals, cabs, sigma_lo, half_eps, k2t, arb(Lc))
        if not (M0c.lower() > 0):
            continue
        if best is None or k2u.rad() < best[0].rad():
            best = (k2u, M0c, Lc)
    if best is None:
        return "fail", {"why": "no valid split", "sigma_lo": sigma_lo}
    _, M0, L = best
    kappa2u = DA.kappa2_window_average(cvals, cabs, sigma_lo, half_eps, k2t, arb(L), inv_sigma)
    if kappa2u is None:
        return "fail", {"why": "window average: h too large"}
    NV_c = kt_c[0] * (1 - kt_c[0]) - kt_c[1]
    NV = NV_c + dv * (1 / v) * (kt_b[1] * (1 - 2 * kt_b[0]) - kt_b[2])
    Q = theta * (1 - theta)
    x_lo = arb(max(float(T_MIN), float(N_lo * theta.lower())))
    one_minus_Gx = arb(0, (1 / (2 * x_lo)).upper())
    one_minus_Gnx = arb(
        0, (1 / (2 * N_lo * (1 - theta).lower())).upper()
    )  # 1 - G(N-x) <= 1/(2(N-x))
    Ng = NV / (k2t * Q) + (kappa2u - 1) / k2t + one_minus_Gx / theta + one_minus_Gnx / (1 - theta)
    target = 1 / (3 + 3 * eps - theta)
    margin = Ng - target
    info = {
        "theta": theta,
        "sigma_lo": sigma_lo,
        "NV": NV,
        "kappa2u": kappa2u,
        "Ng": Ng,
        "L": L,
        "margin": margin,
    }
    if margin.lower() > 0:
        return "ok", info
    return "fail", info


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    max_depth = 30
    t0 = time.time()
    ok, nb = DA.bernoulli_lemma(KST)
    print(f"Bernoulli remainder lemma (K={KST}): {ok}  [{nb} boxes]")
    assert ok
    selftest()
    eps_max = arb(1) / N_MIN
    eps_box = arb(eps_max / 2, eps_max / 2)
    # coverage of theta in [0.1, 1/2]: theta(eps, v) is increasing in v; at the v-edges, for every eps in [0, eps_max]
    vlo, vhi = (arb("2.5"), arb("45")) if "--top" in sys.argv else (arb("0.21"), arb("2.6"))
    if "--top2" in sys.argv:
        vlo, vhi = arb("35"), arb("800")  # theta(800; N >= 10000) >= 0.9917 > 0.99
    th_lo = scaled_cumulants(eps_box, vlo, 1)[0]
    if "--top2" in sys.argv:
        # coverage is needed only where the sweep accepts: j' >= JPRIME_MIN, i.e. eps <= (1 - THETA_HI)/JPRIME_MIN;
        # eps -> 0 needs the Stirling form (x >= 20), the rest uses the direct polygamma
        e_cov = min(float(eps_max), (1 - float(THETA_HI)) / JPRIME_MIN)
        e_st = min(e_cov, 1 / (20 * float(vhi)))
        th_a = scaled_cumulants(arb(e_st / 2, e_st / 2), vhi, 1)[0]
        th_hi = arb(th_a.lower())
        if e_st < e_cov:
            th_b = scaled_cumulants(arb((e_st + e_cov) / 2, (e_cov - e_st) / 2), vhi, 1)[0]
            th_hi = arb(min(th_a.lower(), th_b.lower()))
    else:
        th_hi = scaled_cumulants(eps_box, vhi, 1)[0]
    print(
        f"  coverage: theta(eps, {vlo}) <= {float(th_lo.upper()):.5f} < {THETA_LO};  theta(eps, {vhi}) >= {float(th_hi.lower()):.5f} > {THETA_HI}"
    )
    assert th_lo.upper() < toarb(THETA_LO) and th_hi.lower() > toarb(THETA_HI)
    stack = [(arb(0), eps_max, vlo, vhi, 0)]
    n_ok = n_skip = 0
    worst = None
    while stack:
        e1, e2, v1, v2, d = stack.pop()
        if JPRIME_MIN and float(e1) == 0.0 and float(e2 * v2) > 0.05 and d < max_depth:
            # a box touching eps = 0 needs the Stirling form at 1 + x, valid only for x = 1/(eps v) >= 20:
            # split off the part eps <= 1/(20 v_hi); the rest is handled by the direct polygamma
            em = arb(1) / (20 * v2)
            if float(em) < float(e2):
                stack.append((e1, em, v1, v2, d + 1))
                stack.append((em, e2, v1, v2, d + 1))
                continue
        eps = arb((e1 + e2) / 2, ((e2 - e1) / 2).upper())
        v = arb((v1 + v2) / 2, ((v2 - v1) / 2).upper())
        try:
            status, info = certify_box(eps, v)
        except Exception as ex:  # noqa: BLE001 -- an arithmetic failure is a failed box, bisect it
            status, info = "fail", {"why": repr(ex)}
        if status == "ok":
            n_ok += 1
            if worst is None or info["margin"].lower() < worst[0]:
                worst = (
                    info["margin"].lower(),
                    float(eps.mid()),
                    float(v.mid()),
                    info["Ng"].lower(),
                )
            continue
        if status == "skip":
            n_skip += 1
            continue
        if d >= max_depth:
            print("FAILED box:", float(e1), float(e2), float(v1), float(v2), info)
            sys.exit(1)
        if (v2 - v1) > (e2 - e1) * 500 or d % 3 != 2:
            vm = (v1 + v2) / 2
            # push the upper half first so the lower-v half (where the top regime's accepted boxes are) pops first
            stack.append((e1, e2, vm, v2, d + 1))
            stack.append((e1, e2, v1, vm, d + 1))
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
        f"certified: {n_ok} boxes ok, {n_skip} skipped, worst margin N g - 1/(3 + 3 eps - theta) >= {float(worst[0]) if worst else float('nan'):.6f}"
        f" at eps={worst[1]:.2e}, v={worst[2]:.4f} (N g >= {float(worst[3]):.6f})"
    )
    if JPRIME_MIN:
        print(
            f"VERDICT (top regime): N g > 1/(3 + 3/N - theta) on {{ {THETA_LO} <= theta <= {THETA_HI}, t >= {T_MIN}, j' = N - t >= {JPRIME_MIN} }} : True"
            f"   (the corner j' < {JPRIME_MIN} is the harmonic regime, OBSTRUCTION.md)  [{time.time() - t0:.0f} s]"
        )
        return
    print(
        f"VERDICT: N g > 1/(3 + 3/N - theta) on the dense region (a') {{ {THETA_LO} <= theta <= {THETA_HI}, t >= {T_MIN}, N >= {N_MIN} }} : True  [{time.time() - t0:.0f} s]"
    )


if __name__ == "__main__":
    main()
