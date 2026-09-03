> **FINAL STATE OF THIS VALIDATION (read Addendum 2 at the end): all seven items PASS.** The two FAIL
> verdicts in the body below are the FIRST pass, before the repairs of 3 September; Addendum 1 closed
> item 5 and Addendum 2 closed item 4. The recommendation as of the last addendum is
> `independently-validated`.

# Independent validation: `release/scripts/dense_certificate_b.py` (dense regime (b), theta <= 0.05)

**Time (measured with `date`):** Thu Sep 3 00:03:41 RDT 2026 (work started Wed Sep 2 23:03:00 RDT 2026)
**Validator:** independent-validator (did not write the certificate; did not import it, nor `dense_certificate_a.py`, nor `lab/edgeworth_exact.py`; read all three)
**Object validated:** `projects/qg-bootstrap/release/scripts/dense_certificate_b.py`, last changed in commit `c959de262225f6865b641ee72ce4139636dc4ac3` (HEAD at validation time `a24097af1cd5bdba0e8854d6df3c79cb43281b0c`, working tree clean for this file)
**Engine:** python-flint 0.9.0 (fmpq / fmpq_poly / fmpz / arb / acb), Python 3.12.10. No sympy, no `fractions.Fraction`.
**Independence:** separate formulation. Cumulants from the moment-to-cumulant recursion (not the `q(1-q) d/dq` recursion), Faulhaber via `fmpq.bernoulli`, own bivariate series class (list of `fmpq_poly` in zeta), own weight-truncated Edgeworth engine, own Cauchy-integral evaluation of `M(n,t)` from a closed form of the generating function (no Edgeworth at all). Shared dependency: python-flint itself (unavoidable, documented).
**Validation code (session scratchpad):** `v1_cumulants.py`, `v23_series.py`, `edge_engine.py`, `v3b_bands.py`, `v4_tails.py`, `v4b_M.py`, `v5_remainder.py`, `v6_coverage.py`, `v7_exact.py`, `v7_contour.py`. Load-bearing excerpts are quoted below.
**Script under test, own output (run as a subprocess, not imported):** constant term `4/5`, `S(0, zeta) = 176/525`, V^1 coefficient `2/35 x - 704/7875`, `M_S(0.3) = 1.25327e9`, 18 bands all with `min(S - losses) >= +0.3206`, corner loss `3.06e-15`, `VERDICT ... True`, exit 0, 36 s.

## Claim under test

For the doubled odd-square spectrum `(2k-1)^2` (k <= m, each twice), `N = 2m`, `n = N + 1`, `p_t = e_t/C(N,t)`, `g = -Delta^2 log p_t`:
`N g > 4/5` on `{ 0 < V <= 0.17, 0 <= zeta <= 0.24, zeta V <= 1.07e-3 }`, where `V = v^2` is the squared tilt at the saddle of `x = t`, `eps = 1/m`, `eta = eps/V`, `zeta = eta/V`; hence `M(n,t) = n(e^g - 1) > 4/5` for every integer pair with `t >= 628`, `N >= 1260`, `theta = t/N <= 0.05`, `t^2/N >= 1`.

## Validation plan (declared before running)

| # | item | pass criterion (deterministic) |
|---|---|---|
| 1 | cumulant series `kt_j` vs exact sums | `|exact - series| < 1e-30` at (m, v) = (2000, 0.3), (1e5, 0.05), (1e6, 0.4), j = 1..6; leading terms exact in `fmpq` |
| 2 | the identities (NV, nhat, c_j rewrite, monomial rewrite, binomial series) | exact `fmpq` equality where symbolic; numeric agreement at one point to the arb radius |
| 3 | the cancellation: constant term, `S(0)`, V^1 coefficient of `S`, weight-2 Edgeworth part | exact `fmpq` equality from an independent derivation |
| 4 | Ser sup-bound bookkeeping | each rule either proved valid for the objects it is applied to, or a counterexample exhibited with numbers |
| 5 | per-band remainder inputs, factor 3, corner lemma scaling | each scaling claim either derived or a counterexample exhibited with numbers |
| 6 | coverage `khat_1 in [0.30274, 0.33334]` and the three-way union | own bounds contain the script's; exhaustive case split; random grid with zero misses |
| 7 | end-to-end `M(n,t)` at four points | certified arb balls with `M > 4/5`; margin vs `(176/175) theta` reported; two points cross-checked exactly |

## Results

### Item 1 -- cumulant series: PASS

Own Bernoulli cumulants from `kappa_n = mu_n - SUM_{k<n} C(n-1,k-1) kappa_k mu_{n-k}` with all moments `q`; `S_i(m) = P_{2i}(2m) - 4^i P_{2i}(m)` with `P_p` from `fmpq.bernoulli`; exact sums `(1/m) SUM_k kappa_j(q_k)`, `q_k = x/(1+x)`, `x = V (2k-1)^2/(4m^2)`, in arb at 400 bits.

| (m, v) | worst `|exact - series|` over j = 1..6 | PASS |
|---|---|---|
| (2000, 3/10) | 2.99e-119 | yes |
| (1e5, 1/20) | 1.79e-119 | yes |
| (1e6, 2/5) | 4.57e-89 (series truncated at I = 120; V = 0.16) | yes |

Leading terms (`eps -> 0`, `sigma_i -> 1/(2i+1)`): `kt_1 = [1/3, -1/5, 1/7]`, `kt_2 = [1/3, -2/5, 3/7]` -- identical to the script's print. Sample values: `kt_1(2000, 0.3) = 0.028477350162293029788`, `kt_6(1e6, 0.4) = -0.0207775224906194861`.

### Item 2 -- identities: PASS

At m = 1e5, V = 0.09, zeta = 1.2346e-3 (exact sums for `kt_j`):

- `NV = theta(1-theta) - kt_2 = 6.20664378722564e-4`, `V^2 nhat` the same, difference `[+/- 1.7e-117]`; `nhat(0) = 4/45` and `NV/(kt_2 Q)` constant term `4/5` exactly in the validator's series (own derivation).
- `c_j = k_j/(j! sigma^j)` vs `(eta/2)^{j/2-1} khat_j/(j! khat_2^{j/2})`: differences `[+/- 1.5e-118]` (j = 3), `[+/- 5e-121]` (j = 4), `[+/- 1e-123]` (j = 5), `[+/- 7e-127]` (j = 6); `c_1 = delta/sigma = delta sqrt(eta/(2 khat_2))` agrees to 12 digits at `delta = 0.37`.
- Monomial `c_1^2 c_3^2 c_4` (w = 6, # = 3): `5.83580800023e-16` both ways, difference `[+/- 1.6e-130]`. Derivation checked by hand: the exponent of `khat_2` is `-(e_1/2 + SUM j_l e_l/2) = -(w/2 + #)`.
- `1 - x log(1+1/x) = SUM (-1)^{k+1} x^{-k}/(k+1)` at x = 700: both `7.13606170207849e-4`.
- The assembly `N g = NV/(kt_2 Q) + (kappa_2[u]-1)/kt_2 + (1-G(x))/theta + (1-G(N-x))/(1-theta)` was re-derived from `(log P_r)'' = -kappa_2[phi]`, `kappa_2[phi] = kappa_2[u]/(N kt_2)`, `N L_bin = G(x)/theta + G(N-x)/(1-theta)`, `1/kt_2 - 1/Q = NV/(kt_2 Q)`; the window density `(1-|delta|)` gives `E[delta^e] = 2/((e+1)(e+2))` for even e (checked).

### Item 3 -- the cancellation: PASS

Own Edgeworth engine (`edge_engine.py`): `exp(rho)` as the Taylor sum with weight truncation after every product, `M_j = SUM_d i^d (d+j-1)!! [U^d]`, `kappa_2 = M_2/M_0 + A_1^2/M_0^2` with `1/M_0 = SUM (-h)^k`. Monomial counts per weight `{0: 1, 2: 3, 4: 10, 6: 26, 8: 60, 10: 127}`.

- weight-2 part: `12 c4 - 36 c3^2 + 6 c1 c3` -- exact match (`== True`).
- weight-4 part (for the record): `-90 c6 + 384 c4^2 + 750 c3 c5 - 3924 c3^2 c4 + 4050 c3^4 - 60 c1 c5 + 504 c1 c3 c4 - 864 c1 c3^3 - 12 c1^2 c4 + 54 c1^2 c3^2`.
- `N g`: V^0 = `4/5` (zeta-free), V^1 = `176/525`, V^2 = `2/35 zeta - 704/7875`, V^3 = `-8/2625 zeta + 84752/1819125`. So `S(0) = 176/525` and `[V^1] S = 2 zeta/35 - 704/7875` -- both exact matches with the script's prints, and `3 * 176/525 = 176/175` (the sparse constant).
- The cancellation made explicit, piece by piece: `piece2 (Edgeworth)`: V^0 `-9/4 zeta`, V^1 `9/4 zeta^2 - 27/10 zeta`; `piece3 ((1-G(x))/theta)`: V^0 `+9/4 zeta`, V^1 `-9/4 zeta^2 + 27/10 zeta`; `piece1`: `4/5 + 176/525 V`; `piece4`: 0 at these orders.

```python
# validator, item 3 (excerpt): the weight-2 check and the assembly's constant term
w2_expected = {
    cvar(4): fmpq(12),
    tuple(2 * x for x in cvar(3)): fmpq(-36),
    tuple(x + y for x, y in zip(cvar(1), cvar(3))): fmpq(6),
}
print("weight-2 part == 12 c4 - 36 c3^2 + 6 c1 c3 :", byw.get(2, {}) == w2_expected)  # True
S = (Ng - fmpq(4, 5)).divV()
print(
    S.c[0] == fmpq_poly([fmpq(176, 525)]), S.c[1] == fmpq_poly([fmpq(-704, 7875), fmpq(2, 35)])
)  # True True
```

**Why the `zeta^k V^0` terms must cancel.** In the variables `(V, zeta)`, letting `V -> 0` at fixed `eta = zeta V` means `zeta -> infinity` with `t = 2 khat_1/eta` fixed and `N = 2/(eta V) -> infinity`. The `V^0` part of `N g` is a polynomial `SUM_k g_k zeta^k`; a finite limit at fixed `t` forces `g_k = 0` for `k >= 1` and `g_0 = lim_{N -> inf, t fixed} N g`. That limit is the same for every `t`: the sparse expansion has `N g - 4/5 = O(t^2/N)` (the `[a^0]` and `[a^1]` parts of `Phat` vanish identically), and at `t = 1` it is elementary, `N g = N log[(N-1)/(N - m_2/m_1^2)] -> m_2/m_1^2 - 1 = Var(lambda)/E[lambda]^2 = (1/5)/(1/9) - 1 = 4/5` for the limiting law `lambda = s^2`, `s` uniform on `[0,1]`. Term by term: a monomial `zeta^k V^n` corresponds to `a^n b^{n-k}` (`a = 1/t`, `b = t^2/N`), and `N g - 4/5 = O(b)` at every fixed `a` forces `k <= n - 1`; at `V^0` only the weight-2 Edgeworth term (`-9/4 zeta`, from `12 c_4 - 36 c_3^2` with `khat_j(0) = 1/3`) and the `k = 1` binomial term (`+zeta/(4 khat_1^2) = +9/4 zeta`) appear, and they cancel; at `V^1` the `zeta^2` and `zeta^1` terms cancel likewise (table above), leaving `176/525`.

Independent evaluation of the truncated `S` (degree 34, polynomial part only) on the script's bands gives point minima 0.32126, 0.32796, 0.33152, 0.33336, 0.33429, 0.33476, 0.33500, 0.33512 (bands 1-8, always at `V = hi`, `zeta = 0`), each above the script's interval lower bounds 0.3206, 0.3278, 0.3315, 0.3333, 0.3342, 0.3343, 0.3349, 0.3351 -- consistent (the script's numbers are box lower bounds minus losses).

### Item 4 -- the Ser sup-bound bookkeeping: FAIL as written (numbers survive a stated repair)

The algebraic rules, taken one by one, for functions analytic on the polydisc `|V| <= r = 0.3`, `|zeta| <= 0.24`:

- product: `sup|fg| <= sup|f| sup|g|` pointwise -- valid.
- `1/(c + g)` with `M(g) < |c|`: `sup|1/(c+g)| <= 1/(|c| - M(g))` -- valid, and `c + g` has no zero on the polydisc so the inverse is analytic there. The script's `Mrest` formula (`SUM_{i>=1} poly_sup r^i + (M_self - SUM_{i<=I} poly_sup r^i)`) is **not** valid for a general Ser (for a product, `M_a M_b` can be below the coefficient sum and the bracket goes negative), but it is valid for the three series actually inverted: for `k1`, `k2` (built by `from_terms` then `divV`) the bracket equals `extra_tail/r`, a genuine sup bound of the tail by the maximum principle; for `1 - theta` it equals `M(kt_1)`.
- `f/V` with `f(0) = 0`: `M/r` by the maximum principle -- valid for `f` analytic on the closed disc.
- tail `M (V/r)^{I+1}/(1 - V/r)` -- valid by Cauchy for `f` analytic on `|V| <= r` with `|f| <= M`.

**The gap is in the base objects.** The script bounds `kt_j` as the function `F_j(V, zeta) = SUM_i coef_i V^i sigma_i(zeta V^2)` with `sigma_i` the Faulhaber *polynomials*. Off the integers `1/eps` the Faulhaber polynomial is not a Riemann sum: its `eps^{2i}` coefficient is `~ -B_{2i}`, so on the polydisc (`|eps| <= 0.24 * 0.09 = 0.0216`) `sup|sigma_i|` grows super-exponentially and the `i`-sum **diverges** -- `F_j` is not an analytic function on the polydisc, and the docstring's "bound on its supremum over the disc" refers to an object that does not exist. Exact numbers (own Faulhaber, absolute coefficient sums):

| i | `sup_{|eps|<=0.0216} |sigma_i|` | `1/(2i+1)` | `|s_{i,2i}|` | `|s_{i,2i}| 0.0216^{2i}` |
|---|---|---|---|---|
| 60 | 0.011131 | 0.00826 | 2.2e103 | 3.0e-97 |
| 100 | 0.012726 | 0.00498 | 3.6e215 | 2.8e-118 |
| 166 | first `> 2` | | | |
| 200 | 1.5e7 | 0.00249 | 6.9e549 | 4.1e-117 |
| 400 | 5.9e129 | 0.00125 | 4.4e1338 | 1.6e6 |
| 460 | 1.8e186 | 0.00109 | 9.8e1594 | 4.9e62 |

Consequences for the script's `cumulant_series` tail (lines 262-286):

1. "factor 2 covers `sup|sigma_i(eps)| <= 2` for eps <= 0.03" is false for `i >= 166`; with the script's own truncation `i < 460` and the true sups, `tail1` would be `9.6e83` (j = 1) to `5.4e94` (j = 13) instead of `7.3e-14` to `4.1e-3`. The number the script calls `M` is therefore not a sup bound of the object it names.
2. The dropped-Faulhaber-term bound (line 280-286) multiplies by a spurious extra `rr**(2k)`: `zmax^k r^{i+4k}` instead of the correct `zmax^k r^{i+2k}`. Corrected values: `2.4e-28` (j = 1) ... `2.5e-9` (j = 13) versus the script's `4.7e-35` ... `5.5e-14`.
3. Nothing bounds `i > 460` at all.

**Repair, verified numerically.** Take the base object to be the *polynomial* `G_j := SUM_{i<=60} coef_i V^i sigma_i(zeta V^2)` (full Faulhaber polynomials, no `tail1`). It is analytic; every derived Ser (products, the three inverses, V-divisions of series vanishing at `V = 0`) is analytic on the polydisc; the sup bounds and the Cauchy tails are then rigorous as coded. With the corrected dropped-term bound, `M(G_j)` equals the script's `M(kt_j)` to every printed digit (j = 1: `0.123078`, matching the script's `0.12308`; j = 2: `0.152803`; j = 13: `2.05846e6`), because the tails are negligible in both accountings. The step the script is missing is the passage back to the true cumulants at the real points: for integer `m`, `sigma_i(1/m) in [0, 1]`, so `|kt_j - G_j(V, 1/(m V^2))| <= SUM_{i>60} |coef_i| V^i`, exactly `1.4e-47` (j = 1) ... `3.8e-26` (j = 13) at `V = 0.17`, i.e. at most `1.3e-24` after the worst division by `V^2` in `nhat`, and `O(V^{59})` below; the effect on `S` is below `1e-20` on every band against margins `>= 0.32`. This step must be written (a bound on the sensitivity of `N g` to the `khat_j`), not assumed. Verdict: the numerical conclusion of the band sweep stands; the certificate text as written does not prove it.

```python
# validator, item 4 (excerpt): the true sup of sigma_i on the polydisc, absolute coefficient sum
pe = sigma_coeffs(i)  # eps-coefficients of 4^{-i} eps^{2i+1} S_i(1/eps), exact fmpq
s = sum(abs(toarb(ck)) * epsmax**k for k, ck in enumerate(pe) if ck != 0)  # i = 460: 1.7959e+186
```

### Item 5 -- remainder per band and the corner lemma: bands PASS, corner FAIL

**Bands.** The inputs to `DA.kappa2_enclosure` were checked against the identities of item 2: `|c_j| <= (eta_max/2)^{j/2-1} sup|khat_j|/(j! k2lo^{j/2})`, `c_1 <= 1/sigma_lo`, `sigma_lo^2 = 2 k2lo/eta_max`, `1/N in [0, eta_max hi/2]`, `kt_2 in [k2lo lo, k2hi hi]`, `W = 10` forced on the DA engine (`DA.W`, `DA.EXPO` reset) -- all consistent with the series' weight-10 truncation. `R_edge/(k2lo lo^2)` is the correct conversion of an absolute error in `kappa_2[u]` into `S`. The binomial cuts are first omitted terms of alternating series with decreasing terms (x >= 563), valid.

**Factor 3.** At the worst band (sigma >= 22.5) the validator's own Edgeworth polynomials with the band's c-bounds give `|M_0 - 1| <= 0.0040`, `|A_1| <= 0.081`, `|M_2 - 1| <= 0.020`; the first-order sensitivities of `M_2/M_0 + A_1^2/M_0^2` are `1.042 d_0 + 0.163 d_1 + 1.004 d_2`, and the non-linearised bound at `d_j = 1e-6` (the bands have `<= 1e-8`) is `2.22e-6 < 9e-6`. So `3(d_0 + d_1 + d_2)` is justified -- but the script asserts none of `M_0 ~ 1`, `|A_1| << 1`; it should.

**Corner (V below the last band).** Scaling claims checked: `c_j = khat_j/(j! khat_2) sigma^{2-j}` for `j <= 13` (from `sigma^2 = 2 khat_2/eta`), `K_j` sigma-free -- valid; `c_1 <= 1/sigma`, `e^{c_1 L} <= e^{1/2}` at `L = sigma/2`, `beta_{j<=13} = SUM K_j 2^{2-j}` sigma-free, Gaussian moments sigma-free -- valid; the outer pieces decay as `exp(-0.1146 sigma^2)` and `exp(-0.4596 sigma^2)` -- valid; `sigma^{-11} -> (V/V_last)^{3.5}` for the loss -- valid arithmetic. **But the `j >= 14` cumulant tail used by `dense_certificate_a` is `c_j <= 2.7 N (3 sigma)^{-j}`, which carries a factor `N = sigma^2/kt_2`.** Its contribution to `beta` with `L = sigma/2` is

    beta_geo = 2.7 (6/5) 4 6^{-14} N/sigma^2 = 1.654e-10 / kt_2 ,   kt_2 = V khat_2 -> 0 ,

not sigma-free: `beta_geo = 3.8e-4` at `V = 1.3e-6` (last band), `0.0496` at `V = 1e-8`, `0.397` at `V = 1.25e-9`, `4.96` at `V = 1e-10`. The DA machinery requires `gamma = 1/2 - beta > 0.1`, which fails below `V* = 1.24e-9`, i.e. for `t = 2 khat_1/(zeta V) >= 2.2e9` (with `b >= 1`). The corner lemma's sentence "beta and `e^{c_1 L}` are sigma-free" is false for this term, and the monotone extrapolation `R_inner(sigma) <= R_inner(sigma_last)(sigma_last/sigma)^11` does not follow: the piece (a') is not even finite there. Two lesser inaccuracies in the same lemma: the `|u|^14` monomial `full_co[14] = const sigma^{-12}/kt_2 = const' zeta sigma^{-10}` has effective weight 10, not `>= 11` (harmless: `(V/V_last)^3` is still decreasing); and the `K_j` are taken from the last band `[V_last, 2 V_last]` rather than the corner `(0, V_last]` (harmless at `O(1e-6)` relative, but should be the corner sups from the series evaluated on `[0, V_last]`).

**Repair.** Replace the `j >= 14` tail bound by one of the form `K_j sigma^{2-j}` with sigma-free `K_j` (e.g. from `|kappa_j(q)| <= q C_j`, so `|khat_j| <= C_j khat_1`, or from the radius `pi` of the Bernoulli cumulant generating function), or choose `L = L(V)` such that `N (L/(3 sigma))^{14}/L^2` stays bounded (`L <= sigma kt_2^{1/12}/C` works and keeps the outer pieces decaying like `exp(-c V^{-5/6})`). Either way the corner lemma has to be rewritten and re-run. Until then the certificate covers `V >= 6.485e-7` (the bands) rigorously modulo item 4, and the corner `V < 6.485e-7` only for `V >= 1.24e-9` by the arithmetic of the existing lemma with the exponent corrected to 3.

### Item 6 -- coverage: PASS

Own derivation without the series: `khat_1 = (1/m) SUM_k s_k^2/(1 + V s_k^2)`, `s_k = (2k-1)/(2m)`, decreasing in `V`; upper bound `khat_1(0) = 1/3 - 1/(12 m^2) < 0.33334`; lower bound by the midpoint rule (`|f''| <= 2`): `khat_1(0.17) >= (1 - arctan(sqrt 0.17)/sqrt 0.17)/0.17 - 1/(12 * 630^2) = 0.30298221 > 0.30274` (exact sums: `0.3029822683` at m = 630, `0.3029824217` at m = 1e5). The script's interval `[0.30274, 0.33334]` is contained in the validator's `[0.30298, 0.33334)`. Implications: `theta <= 0.05 => V <= 0.165158 <= 0.17` (not circular: `V > 0.17` would give `theta = kt_1(V) > 0.17 * 0.30274 = 0.0515`); `b >= 1 => zeta = 2 khat_1^2/b <= 0.222231 <= 0.24`; `t >= 628 => eta = 2 khat_1/t <= 1.06159e-3 <= 1.07e-3`. Case split for integer `(N, t)`, `t >= 628`, `N >= 1260` even, `t <= N/2`: `b <= 1` -> sparse (every `t >= 627`); `b >= 1, theta >= 0.05` -> dense (a); `b >= 1, theta <= 0.05` -> dense (b). Exhaustive by construction; 200 000 random pairs, zero outside all three. `N in {1256, 1258}` handled exactly in `theorem.py` (B).

### Item 7 -- end-to-end `M(n,t)`: PASS

Two independent computations. (i) Exact: power sums by Faulhaber, `e_j` by Newton's identities in `fmpz`, `M` as an exact `fmpq` (8796-digit `e_630`, 29309-digit `e_2000`). (ii) Certified arb by a separate formulation: `E(z) = PROD (1 + z(2k-1)^2)^2` in closed form,

    log E(z) = 2m log(1+X) - 2(2m+1) theta0(X) + 2[St(w+) + St(w-)] + 2 log(1 + e^{-2 pi y}),
    X = 4 z (m+1/2)^2,  y = (m+1/2)/sqrt X,  w+- = m + 1/2 +- i y,  theta0(X) = 1 - arctan(sqrt X)/sqrt X

(Stirling with the explicit complex remainder `|R_4| <= |B_10|/(90 |w|^9 cos^10(arg w/2))`, reflection formula for the small gammas; checked against the direct product at m = 50), `e_k r^k/E(r)` as `acb.integral` over `|phi| <= phi0` with a mean-value enclosure of the large part (`P(phi) in P(phi_c) + P'(ball)(phi - phi_c)`), holomorphy checked on every ball, tail `2(pi - phi0) exp(-k_2(1 - cos phi0))` with `k_2` bounded below by the midpoint rule. Self-tests against the exact integer DP: `(2000, 100)`: `0.853596274420188 +/- 1.9e-16` contains the exact value; `(4000, 300)`: `0.881887184787830 +/- 2.0e-16` contains it.

| (N, t) | theta | b | `M(n,t)` (certified ball) | exact `M` | `M - 4/5` | `(176/175) theta` |
|---|---|---|---|---|---|---|
| (400000, 630) | 1.575e-3 | 0.992 | `0.801589591487799 +/- 3.8e-16` | 0.801589591488 | 1.5896e-3 | 1.5840e-3 |
| (1e6, 2000) | 2.000e-3 | 4.0 | `0.802016697316824 +/- 3.0e-16` | 0.802016697317 | 2.0167e-3 | 2.0114e-3 |
| (1e8, 20000) | 2.000e-4 | 4.0 | `0.800201195440152 +/- 2.5e-16` | -- | 2.0120e-4 | 2.0114e-4 |
| (2e7, 1e6) | 5.000e-2 | 5e4 | `0.852931556868305 +/- 1.1e-16` | -- | 5.2932e-2 | 5.0286e-2 |

`M > 4/5` at all four; the margin is `(176/175) theta` to 0.3 percent at small `V` and, at `(2e7, 1e6)` (`V = 0.1646`), equals `V S(V)` with the validator's `S(0.1646, ~0) = 0.3216`: `0.1646 * 0.3216 = 0.0529` -- the series and the direct computation agree there to three digits.

```python
# validator, item 7 (excerpt): mean-value enclosure of the large part of the exponent on a ball phi
Pc = 2 * m * (1 + Xc).log() - 2 * (2 * m + 1) * th0c - acb(0, 1) * k * phic  # thin centre
dP = acb(0, 1) * (2 * m * X / (1 + X) - 2 * (2 * m + 1) * xth1 - k)  # P' on the whole ball
P = Pc + dP * (phi - phic)  # rigorous: the ball is convex, P analytic
```

## Failure map

| region / step | status | size of the defect |
|---|---|---|
| series algebra, cancellation, `S` coefficients | reproduced exactly | -- |
| bands `V in [6.485e-7, 0.17]`, polynomial part of `S` | reproduced (point minima 0.3213..0.3351 vs script 0.3206..0.3352) | -- |
| Cauchy tails of every Ser (item 4) | invalid as stated: the limiting object is not analytic; `sup|sigma_i| <= 2` false for `i >= 166`; spurious `r^{2k}`; no bound for `i > 460` | numerically nil after restating with the polynomial `G_j` (`M` unchanged to 6 digits; real-point remainder `<= 1.3e-24`); the restatement and the sensitivity step are missing from the text |
| Edgeworth remainder on the bands (item 5) | inputs and the factor 3 verified | `M_0 ~ 1`, `A_1 << 1` unasserted |
| corner `1.24e-9 <= V < 6.485e-7` | lemma arithmetic holds with exponent 3 instead of 3.5 | documentation |
| corner `V < 1.24e-9` (`t >= 2.2e9`, `b >= 1`) | **not covered**: `beta_geo = 1.65e-10/kt_2` exceeds the `gamma > 0.1` threshold | genuine gap; repair specified |
| coverage of the integer region | verified independently | -- |
| docstring line 33 (`S(0, zeta) = 0.64 + 0.54 zeta`) | stale, contradicts the print `176/525` | text only |

## Claim-state recommendation

`dense_certificate_b.py`: **not** `independently-validated`. Recommended state: `experimentally-supported` (every number the script prints is reproduced by an independent formulation; the end-to-end values at four points confirm `M > 4/5` with the predicted margin) with two blockers before promotion:

1. Item 4: restate the Ser objects as the polynomials `G_j` (drop `tail1`, fix the `r^{2k}` factor) and add the bound on the effect of `SUM_{i>60} coef_i sigma_i(1/m) V^i` on `N g` at the real points (order `1e-24`).
2. Item 5: rewrite the corner lemma with a sigma-free cumulant tail for `j >= 14` (or a `V`-dependent `L`), and re-run; the present lemma leaves `t >= 2.2e9`, `theta <= 0.05`, `t^2/N >= 1` unproved.

Neither blocker changes any printed number; both change what the text proves. Re-validate after the fixes (items 4 and 5 only).

```yaml
# templates/validation.yaml schema
id: VAL-DENSE-B-2026-09-02
experiment_id: dense_certificate_b
validator: independent-validator
independence_level: separate-formulation   # shared dependency: python-flint 0.9.0
checks:
  known_answer: pass          # items 1-3, 6, 7 (exact e_t at two points; DP self-tests)
  negative_control: pass      # item 4: sup|sigma_i| computed exactly refutes the stated bound at i >= 166
  hidden_points: pass         # item 7: (1e8, 20000), (2e7, 1e6) by contour, not by the certificate's method
  precision_convergence: pass # item 1 to 1e-89; item 7 balls 1e-16
  boundary_map: fail          # item 5: corner V < 1.24e-9 not covered by the lemma as written
  signature: pass             # constant 4/5, 176/525, 2 zeta/35 - 704/7875 reproduced exactly
decision: fail                # as a proof text; numerical content reproduced
allowed_claim_state: experimentally-supported
blockers:
  - "item 4: Cauchy-tail bookkeeping refers to a non-analytic object; restate with the polynomial G_j and add the real-point remainder step"
  - "item 5: corner lemma: beta_geo = 1.65e-10/kt_2 is not sigma-free; gamma > 0.1 fails below V = 1.24e-9 (t >= 2.2e9)"
```

**Addendum time (measured with date):** Thu Sep  3 00:23:36 RDT 2026

---

## Addendum: re-check of items 4 and 5 after the repairs (commit `43d2775`)

**Object re-checked:** `dense_certificate_b.py` and the `geo_K` branch of `dense_certificate_a.py` at commit `43d2775` (rerun log `results/dense_certificate_b_rerun_2026-09-03.txt`: bands unchanged, corner remainder 6.40e-29, loss 4.57e-16, VERDICT True, exit 0). Validation code: `v45_recheck.py` (scratchpad), same independence rules as above (nothing imported).

### Item 4 (sup bounds at fixed integer m): M bounds PASS; `evaluate()` bookkeeping FAIL (residual, numerically harmless)

What now holds, verified:

- Fixed `m`: `kt_j(V; m)` is a rational function of `V` with poles at `V = -4m^2/(2k-1)^2`, `|.| > 1`, so it is analytic on `|V| <= 0.3`; its Taylor coefficients are the numbers `coef_i sigma_i(1/m)`. Valid.
- Riemann bound `sigma_i(1/m) <= (1 + eps/2)^{2i+1}/(2i+1)`: derivation checked (`(2k-1)^{2i} <= INT_{k-1/2}^{k+1/2} (2x)^{2i} dx`), and exact values tested: (i, m) = (1, 5): 0.3300 vs 0.4437; (60, 46): 6.29e-3 vs 3.06e-2; (60, 630): 8.252e-3 vs 9.097e-3; (200, 100): 1.37e-3 vs 1.84e-2 -- all `ok`.
- Tail 1 closed form `A_j (1+eps/2)/(2(2I+3)) q^{I+1}/(1-q)`, `q = 2r(1+eps/2)^2 = 0.613`: dominates the direct sum (j = 1: 1.155e-15 vs 1.127e-15; j = 13: 6.49e-5 vs 6.33e-5). `eps <= ZMAX r^2 = 0.0216` holds on the region with room (`eps = eta V <= 1.82e-4`).
- Tail 2 without the spurious `r^{2k}`: identical to the validator's corrected values (2.40e-28 for j = 1, 2.46e-9 for j = 13).
- Domination of the kept part: a kept monomial `coef_i s_ik eps^k V^i` has sup `|coef_i s_ik| (eps/r^2)^k r^n` on `|V| = r` (`n = i + 2k`), which is at most `|coef_i s_ik| ZMAX^k r^n` because `eps/r^2 <= ZMAX`. Valid, and `M(kt_1) = 0.123078`, `M(kt_13) = 2.05846e6` agree with the validator's `M(G_j)` to all digits.
- Analyticity in zeta: nothing in the Ser algebra or in `evaluate()` needs it. Each Ser is, for each fixed `m`, a function analytic in `V` on the disc with sup `<= M` (products, the three inverses with the gap check, V-divisions by the maximum principle all preserve this), and `evaluate()` applies Cauchy in `V` only. The polynomials in zeta are bookkeeping.

What still does not hold:

1. **The fixed-`m` Cauchy tail does not cover the dropped terms.** At fixed `m` the polynomial part `SUM_{n<=I} c_n(eps/V^2) V^n` keeps the monomials `eps^k V^i` with `i + 2k <= 60`; the omitted monomials with `i <= 60`, `i + 2k > 60` have V-degree `<= 60` and are **not** part of the Taylor tail beyond degree 60, so `M (V/r)^{61}/(1 - V/r)` does not bound them (it reaches them only through `tail2 * (V/r)^{61}/(1-V/r)`, short by the factor `(V/r)^{61}/(1-V/r) ~ 2e-15`). Their true size at the real points is tiny -- `D <= 2.6e-51` (j = 1), `3.2e-30` (j = 13) at `V = 0.17`, `eps <= 1.82e-4` -- but the text has no bound for them. The polynomial-`G_j` view of the first report does cover them (there they have total degree `n > 60` and the two-variable Cauchy tail at fixed zeta applies); the new fixed-`m` comment does not. Either state the `G_j` view, or add `D` explicitly to `evaluate()`.
2. **`divV` loses one exactly-known degree each time and `evaluate()` does not track it.** `divV` sets the coefficient of `V^60` to zero although the true coefficient (`V^61` of the parent) is unknown; the tail formula still starts at `(V/r)^{61}`. `khat_j` is exact to degree 59, `nhat` to 58, `S` to 57. At `V = 0.17` the missing degrees of `S` alone are bounded by `M_S [(V/r)^58 + (V/r)^59 + (V/r)^60] = 1.17e-5`, versus the tail used, `2.60e-6`; the correct tail is `M_S (V/r)^{58}/(1-V/r) = 1.43e-5` (band-1 margin 0.3206, so harmless). For `khat_13` on band 1 the missing degree-60 term is `<= 1.09e-8` versus the tail used, `1.42e-8` (covered by accident only because `V/r >= 1/2` there; not covered on the lower half of the band). Fix: carry an exact-degree field, decremented by `divV`, and use `(V/r)^{d+1}/(1-V/r)` in `evaluate()`.
3. **The `i > 60` real-point remainder is still not propagated in writing.** Folding tail 1 into `M` covers it for `kt_j` itself (`tail1 (V/r)^{61}/(1-V/r) >= delta_j` since `(1-q_r)(1-V/r) <= 1 - q_V`), but after each `divV` the true effect `delta_j/V` exceeds the folded route by `r(1-q_r)(1-V/r)/V`, i.e. for `V < 0.116`. Size: `delta_13 <= 9.2e-21` at `V = 0.17` and `~(2V)^{61}` below, so the effect on `S` is below `1e-10` everywhere. One sentence with a sensitivity constant closes it.

Total under-accounted tail at the worst point (`V = 0.17`): `<= 1.5e-5` against `min(S - losses) = 0.3206`. The numbers of the certificate are unaffected; the proof text is not yet closed on item 4.

### Item 5 (corner with `geo_K`): PASS

- Constant 7.05: `q_max = 0.17/1.17 = 0.14529915 <= 0.1453`; `q(e^{3/2} - 1) = 0.505886 < 1`; `(e^{3/2}-1)/(1 - q(e^{3/2}-1)) = 7.04633 <= 7.05`. The chain `|log(1+w)| <= -log(1-|w|) <= |w|/(1-|w|)`, Cauchy on `|s| = 3/2` for the Bernoulli cumulant generating function `K(s) = log(1 + q(e^s - 1))` (analytic there since `|q(e^s-1)| < 1`), gives `|kappa_j(q)| <= 7.05 j! (2/3)^j q`. Verified.
- `q` bound: `q_k = b_k r/(1 + b_k r)` with `b_k r <= (2m-1)^2 V/(4m^2) < V`, so `q_k < V/(1+V)`. Verified.
- `SUM_k 2 q_k = k_1 = t` at the saddle, `c_j = k_j/(j! sigma^j)`, `t/sigma^2 = N kt_1/(N kt_2) = khat_1/khat_2`. Verified; `Kg = 7.05 khat_1/khat_2`, `rho = 3/2`.
- With `L = sigma/2`: `full_co[14] = Kg (3/2)(2/3)^{14} sigma^{-12}` (weight 12, sigma-free constant 0.0362) and `beta_geo = 6 Kg 3^{-14} = 8.84e-6`, sigma-free; the `N/sigma^{14}` term that broke the lemma is gone. Every inner piece is now `C_w sigma^{-w}` with `w >= 11` and sigma-free `C_w` (`c_1 <= 1/sigma`, `e^{c_1 L} = e^{1/2}`, `beta = SUM K_j 2^{2-j} + 8.84e-6`, `gamma > 0.1`, Gaussian moments sigma-free), the outer pieces decay as `exp(-0.1146 sigma^2)` and `exp(-0.4596 sigma^2)`, so `R_inner(sigma) <= R_inner(sigma_last)(sigma_last/sigma)^{11}`, `loss(V) <= loss_last (V/V_last)^{3.5} + R_outer/(k2lo V^2)`, both decreasing on `(0, V_last]` with **no lower limit on V**. The DA `geo_K` branch computes exactly this (`Kg sigma_lo^2 ratio^14/(1 - ratio)`, `ratio = L/(rho sigma_lo) = 1/3`).
- One bookkeeping slip, harmless: the corner `K_j` must be sups over `(0, V_last]`, but `kh` is evaluated on the box `[V_last/2, V_last]`. `khat_j = 1/3 - 2^{j-1} V/5 + ...` is decreasing, so the corner sup is `1/3` and the band box under-reads it by a relative `3 * 2^{j-1} V_last/5`: 1.6e-6 (j = 3), 5.0e-5 (j = 8), 1.6e-3 (j = 13). `khat_1/khat_2` (increasing) and `k2lo` (decreasing) take their extremes at `V_last`, inside the box. The remainder `6.4e-29` would grow by well under a factor 2 with the correct sups, against a margin of 0.3352; evaluate `kh` on `[0, V_last]` to close it.

### Claim-state recommendation (updated)

Item 5 closed. Item 4 remains a text-level FAIL with a bounded, harmless residual (`<= 1.5e-5` at `V = 0.17`, far less elsewhere): (i) exact-degree tracking through `divV` in `evaluate()`, (ii) either the polynomial-`G_j` justification or an explicit `D` term for the dropped monomials, (iii) one written sentence propagating the `i > 60` remainder. Recommended state unchanged: `experimentally-supported`; promote to `independently-validated` after (i)-(iii) and a rerun (item 4 only; no other item needs re-validation).

**Addendum 2 time (measured with date):** Thu Sep  3 00:35:49 RDT 2026

---

## Addendum 2: re-check of item 4 (and the item-5 corner window) after commit `6393d87`

**Object re-checked:** `dense_certificate_b.py` at commit `6393d87` (Ser with exact degree `d` and real-point perturbation `P`; base object the polynomial `G_j`; corner window `(0, V_last]`). Rerun log `results/dense_certificate_b_rerun2_2026-09-03.txt`: `d(S) = 57`, `P_S(VMAX) = 3.70e-18`, bands unchanged (`min(S - losses) >= 0.3206`), corner remainder 6.40e-29, VERDICT True, exit 0. Validation code: `v4d_pbase.py` (scratchpad); the rules below were checked by derivation, nothing imported.

### Item 4: PASS

The certificate now has a two-layer structure and each layer was checked separately.

**Layer 1 -- the analytic objects.** Base `G_j := SUM_{i<=60} coef_i V^i sigma_i(zeta V^2)` (Faulhaber polynomials, finite in `i`), a polynomial in `(V, zeta)`. `from_terms` stores its coefficients of V-degree `<= 60` exactly (`d = 60`) and bounds its polydisc sup by the kept-coefficient sums plus tail 2 (the terms of `G_j` of V-degree `> 60`, `|coef_i| |s_ik| ZMAX^k r^{i+2k}`); tail 1 no longer enters `M`, so `M` is a sup of a genuine analytic function. Derived objects are the exact algebraic images (sums, products, `1/(c0+g)` with `M(g) < |c0|` so the inverse is analytic on the polydisc, `f/V` with `f(0, zeta) = 0` exactly), all analytic on the polydisc with the propagated `M` as sup bound. Exact-degree bookkeeping: `from_terms: 60`, `const/V/zeta: infinite`, `add/mul: min`, `inv: d` (`h_i` uses `c_k`, `k <= i`), `divV: d - 1`, scalar: `d`; `S` ends at `d = 57` as I derived in Addendum 1. `evaluate()` sums `n <= min(d, 60)` and takes the Cauchy tail `M (V/r)^{d+1}/(1 - V/r)` at fixed zeta -- valid because every omitted coefficient of the analytic object is a genuine Taylor coefficient of degree `> d` in `V` (at fixed zeta) bounded by `M/r^n`. This answers the second question: in the `G_j` formulation every dropped monomial of every derived object (including the base tail-2 terms, which are the terms of `G_j` itself of V-degree `> 60`) is part of the Taylor tail of the analytic object beyond its exact degree; no monomial of degree `<= d` is missing anywhere on the path from the `G_j` to `S`.

**Layer 2 -- the real-point perturbation.** On the real region the true `N g` is the same rational expression in `(kt_j, V, zeta)` that the Ser algebra evaluates in the `G_j` (the variance identity for piece 1, the weight-`<= 10` Edgeworth polynomial for piece 2, the cut binomial series for pieces 3-4 are exact identities; the Edgeworth remainder and the binomial cut are accounted separately as before). Hence the only discrepancy between the true value and the analytic object is `kt_j - G_j(V, eps/V^2) = T_j = SUM_{i>60} coef_i sigma_i(1/m) V^i` at the integer `m`, which is what `P` carries. Checks:

- Base: `|T_j(V)| <= A_j (1+eps/2)/(2(2I+3)) (2(1+eps/2)^2)^{61} V^{61}/(1 - q)`, `q = 2 VMAX (1+eps/2)^2 = 0.347 < 1`, monotone in `eps <= eps_max`. Recomputed: `P_base = 5.383e16 V^61` (j = 1), `1.077e17 V^61` (j = 2), `1.747e26 V^61` (j = 12), `3.025e27 V^61` (j = 13); at `V = 0.17` these are `6.14e-31`, `1.23e-30`, `1.99e-21`, `3.45e-20`, each above the direct 2000-term sum of the Riemann-bounded terms (`6.09e-31`, `1.22e-30`, `1.98e-21`, `3.42e-20`).
- Product rule `P_fg = M_f P_g + M_g P_f + P_f P_g`: from `(f+df)(g+dg) - fg = f dg + g df + df dg` with `f`, `g` the analytic objects. Using the polydisc sup `M_f` for `|f|` at the real points is correct: the real region `V in [0, 0.17]`, `zeta in [0, 0.24]` lies inside the polydisc `|V| <= 0.3`, `|zeta| <= 0.24`, so `|f_G| <= M_f` there (conservative, never wrong). First question answered: yes.
- Inverse rule `P_{1/f} = P_f/(gap (gap - P_f(VMAX)))`: from `|1/(f+df) - 1/f| = |df|/(|f| |f+df|)` with `|f_G| >= |c0| - M(g) = gap` on the polydisc and `|f_true| >= gap - P_f(V) >= gap - P_f(VMAX)` (coefficients of `P` non-negative, `V <= VMAX`). Third question answered: `P` must enter the real-point gap condition and it does (`assert (gap - pmax) > 0`); it must not enter the analytic gap condition `M(g) < |c0|`, and it does not. Both asserts are present. The `Mrest` formula remains valid for the three inverted series (`k1`, `k2`: bracket = tail2/r, a sup of the V-divided tail by the maximum principle; `1 - theta`: bracket = `M(kt_1)`).
- `divV`: `P_{f/V} = P_f/V` (powers shift down), requires `P` without constant term (asserted) and `f(0, zeta) = 0` exactly (asserted on an exact coefficient since `d >= 0` everywhere on the path).
- `evaluate()` adds `P(V_hi)`; `P` has non-negative coefficients so it is monotone in `V`.
- Sanity of `P_S(VMAX) = 3.7e-18`: the base constants are of order `A_j 5e16 V^61`, the `j!` divisions in the Edgeworth monomials and the `M`-multipliers of the inverses (`6.4` per power) leave the total in the `1e-18` range; against the band-1 margin `0.3206` it is immaterial, and all three residuals of Addendum 1 are now explicit terms of the enclosure (degree slippage: `(V/r)^{58}`; dropped monomials: inside the Cauchy tail of `G_j`; `i > 60` remainder: `P`).

Nothing in the Ser algebra or in `evaluate()` uses analyticity in zeta beyond polynomial dependence; the Cauchy estimate is in `V` at fixed real zeta.

### Item 5, corner window: PASS

`kh` is now evaluated on `Vb = [0, V_last]` (`arb(hi/2, hi/2)`), so the `K_j` (sup of `|khat_j|`, attained at `V -> 0`), `k2lo` and `Kg = 7.05 khat_1/khat_2` are sups/infs over the whole corner; `half_eps = [0, eta_max V_last/2]` and `sigma_lo^2 = 2 k2lo/(ZMAX V_last)` are valid on `(0, V_last]`; the `k2t` box (whose lower end `k2lo * lo` is not a corner-wide bound) is not used in the `geo_K` branch. The remainder stays `6.40e-29` at the printed precision (the corrected `K_13` is 0.16 percent larger, invisible against the dominant low-weight terms), loss `4.57e-16` against `S >= 0.3352`.

### Claim-state recommendation (final for this validation)

All seven items now PASS (items 1-3, 6, 7 from the first report; item 5 from Addendum 1 with the window slip closed here; item 4 here). The certificate text and code, as of commit `6393d87`, prove `N g > 4/5` on the stated region with the two shared dependencies documented (python-flint; `dense_certificate_a.kappa2_enclosure` as the Edgeworth-remainder engine, validated separately by the dense-(a) validator with the `geo_K` branch checked here). Recommended state for `dense_certificate_b`: `independently-validated`. The human still approves any public claim.
