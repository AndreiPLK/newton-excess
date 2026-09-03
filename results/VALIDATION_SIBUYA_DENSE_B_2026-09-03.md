# Independent validation: `release/scripts/sibuya_dense_b.py` (Sibuya dense regime (b'): theta <= 0.1, j^2/N >= 1)

**Time (measured with `date`):** Thu Sep 3 02:44:16 RDT 2026 (work started Thu Sep 3 02:09:08 RDT 2026)
**Validator:** independent-validator (did not write the certificate; imported nothing from it -- its numbers come from its log and from one subprocess re-run)
**Object validated:** `projects/qg-bootstrap/release/scripts/sibuya_dense_b.py`, blob `8ca3fed586f42fde907a2cd13dd07ffa78cf7ae4` = the version committed in `bb46ccee` (HEAD `942dcff7`, working tree clean for this file). Log validated: `results/sibuya_dense_b_2026-09-03.txt` (VERDICT True, EXIT 0, 539 s, constants `v <= 6/25, zeta <= 13/50, zeta v <= 51/100000`, `M_S = 3.29401e12`, `d = 147`, `P_S(VMAX) = 7.36e-34`, band-1 `min(S - losses) = +0.2398`, other bands `>= +0.3234`, corner remainder `3.24e-27`, loss `2.32e-14`). The earlier logs `_run1.txt` (v <= 0.23, zeta <= 0.25, eta <= 2.5e-4) and `_tmin2001.txt` (eta <= 2.6e-4) are superseded by the final constants and were only read for consistency (same series prints, same structure).
**Engine:** python-flint 0.9.0 (fmpq / fmpq_poly / fmpz / fmpz_poly / arb), Python 3.12.10. No sympy, no `fractions.Fraction`, no float in exact comparisons.
**Independence:** separate formulation. Own series engine (list of `fmpq_poly` in zeta, index = v-power; lineage: the validator's own `v23_series.py` of 2 Sept re-derived for the single spectrum), own Bernoulli cumulants (moment-to-cumulant recursion, not the `q(1-q) d/dq` recursion), own Faulhaber via `fmpq.bernoulli` with `B_1 = +1/2`, own Edgeworth engine (`edge_engine.py`, validator's own), own exact `e_j(1..N)` by TWO routes (fmpz_poly product tree with `mul_low`; Eulerian numbers of the second kind, `c(n, n-k) = SUM_m <<k,m>> C(n+m, 2k)`), own tilt solver (arb digamma). Shared dependencies, declared: python-flint; `dense_certificate_a.kappa2_enclosure` is called ONLY in item 5, with the validator's own inputs, to re-run the corner remainder (the engine was validated 2 Sept, its `geo_K` branch 3 Sept).
**Script under test, own subprocess re-run (finished Thu Sep 3 02:52:27 RDT 2026, `script_rerun.log`, 491 s, EXIT 0):** identical to the logged run line by line -- constant `1/3`, `0`, `S(0) = 1/3`, `[v^1] S = (-5/6)*x + (-7/90)`, `M_S = 3.29401e+12`, `d = 147`, `P_S(VMAX) = 7.36e-34`, band 1 `+0.2398`, last band `+0.3333`, corner `3.24e-27 / 2.32e-14`, VERDICT True.
**Validation code (session scratchpad `...\scratchpad\sib_b\`):** `PLAN.md` (frozen before any result), `sb_series.py`, `sb_item1.py`, `sb_item2.py`, `sb_item34.py`, `sb_item5.py`, `sb_item6.py`, `sb_item7.py`; logs `item1.log`, `item2.log`, `item34.log` (I = 40), `item34_150.log` (I = 150), `item5.log`, `item6.log`, `item7.log`, `script_rerun.log`.

## Claim under test

Spectrum `1, 2, ..., N` (unsigned Stirling numbers of the first kind), `n = N + 1`, `p_j = e_j/C(N, j)`, `g = -Delta^2 log p_j`. With the tilt `v` at the saddle (`SUM_k q_k = j`, `q_k = (kv/N)/(1 + kv/N)`), `eps = 1/N`, `eta = eps/v`, `zeta = eta/v`, `theta = j/N = v khat_1`:
`(3 + 3 eps - theta) N g > 1` on `{ 0 < v <= 0.24, 0 <= zeta <= 0.26, zeta v <= 5.1e-4 }`, hence `N g > N/(3N + 3 - j)`, hence Sibuya's (3.4) `p_j^2/(p_{j-1}p_{j+1}) >= 1 + 1/(3n - j)` for every integer pair in the region.

## Validation plan (declared before running; `PLAN.md`)

| # | item | pass criterion (deterministic) |
|---|---|---|
| 1 | series in v; `sigma_i(eps) = eps^{i+1} SUM_{k<=N} k^i` with first-order term; Riemann bound; coefficient bound; P_base | exact fmpq identities; series vs exact sums `< 1e-30` at three (N, v); P_base constant dominates the direct tail |
| 2 | the eps = 1/N conventions | numeric identities at N = 1e5, v = 1/5 (exact sums) to `< 1e-60`; theta(0.24) <= 0.104 so 0.89 and FAC = 3.01 hold |
| 3 | target, cancellation, `S(0) = 1/3`, `[v^1] S`; decoding of (3.4) from the paper's text | exact fmpq equalities from the own series; exact identity for 50 random (n, m) |
| 4 | I = 150, r = 0.3, Cauchy tail at v = 0.24 | own coefficients satisfy `sup|S_n| r^n <= M_S`; tail arithmetic reproduced; own `M(kt_1)` within 1e-3 of 0.19395; khat_1 gap > 0 |
| 5 | corner: Kfac from q_max = 0.24/1.24 on |s| = 3/2; sigma-free K_j; corner remainder with own inputs | `q e32 < 1`; own remainder within a factor 2 of 3.24e-27, loss < 1e-10; exponent 3.5 derived |
| 6 | coverage of {theta <= 0.1, j >= 1001, N <= j^2 or j^2 <= N < (j+1)^2} | exact bounds; grid + 1000 random integer pairs, zero misses |
| 7 | end-to-end exact p_j at integer pairs | known answer, two exact routes agree, negative control fails, every pair passes exactly, `|S_true - S_own|` below the band loss, `S_true` above the band's logged minimum |

## Results

### Item 1 -- the series in v, the Faulhaber factor, the Riemann bound, P_base: PASS

Own derivation: `q = u/(1+u)`, `u = k v/N = k eps v`, `q^p = SUM_i (-1)^{i-p} C(i-1, p-1) u^i`, so `kt_j = (1/N) SUM_k kappa_j(q_k) = SUM_i coef_{j,i} v^i eps^{i+1} SUM_{k<=N} k^i = SUM_i coef_{j,i} v^i sigma_i(eps)` -- one Faulhaber sum, no doubling, no `2^{-2i}`.

- `sigma_i(eps) = 1/(i+1) + eps/2 + ...` for i = 1..12 exactly (`[eps^2] = (i)/12`; e.g. `sigma_2 = 1/3 + eps/2 + eps^2/6`). The first-order term `eps/2` is what distinguishes this port from the odd-square original (which had only even powers). PASS.
- Riemann bound `sigma_i(1/N) <= (1 + 1/N)^{i+1}/(i+1)` (from `k^i <= INT_k^{k+1} x^i dx`): exact fmpq at (i, N) = (1,5), (5,1), (60,46), (60,630), (150,1001), (151,8170), (200,100), (150,10); the Faulhaber polynomial equals the direct sum at every N <= 1001. E.g. (150, 1001): `7.134e-3 <= 7.700e-3`; (151, 8170): `6.640e-3 <= 6.702e-3`. PASS.
- `|coef_{j,i}| <= A_j 2^{i-1}` with `A_j = SUM_p |a_{j,p}|` (`A_1..A_4 = 1, 2, 6, 26`, `A_13 = 5.6e10`): exact for i <= 200, j <= 13. PASS.
- Series vs exact sums (arb, 400 bits, series to i = 220 with exact `sigma_i(1/N)`): worst `|exact - series|` over j = 1..6 is `5.3e-107` at (2000, 3/10), `6.9e-118` at (1e5, 1/20), `6.2e-117` at (2e5, 6/25). Leading terms (eps -> 0): `kt_1 = [1/2, -1/3, 1/4]`, `kt_2 = [1/2, -2/3, 3/4]` -- identical to the script's prints `(1/2)V + (-1/3)V^2 + (1/2 x + 1/4)V^3` and `(1/2)V + (-2/3)V^2 + (1/2 x + 3/4)V^3` (the `x/2` is the `eps/2` term with `eps = zeta v^2`). PASS.
- P_base: `|T_j(v)| = |SUM_{i>150} coef_i v^i sigma_i(eps)| <= A_j (1+e)/(2(I+2)) (2(1+e))^{I+1} v^{I+1}/(1 - 2(1+e)v)`, with `e = eps_max = ZMAX r^2 = 0.0234` (`>=` every real eps, which is `<= 1.2e-4`) and `qq = 2 VMAX (1+e) = 0.4912 < 1`. Recomputed constants `6.209e44` (j = 1), `1.242e45` (j = 2), `2.015e54` (j = 12), `3.489e55` (j = 13); at v = 0.24 they give `1.603e-49, 3.206e-49, 5.202e-40, 9.006e-39`, each above the direct 2850-term sum of the Riemann-bounded terms (`1.593e-49, 3.186e-49, 5.170e-40, 8.950e-39`). PASS.

### Item 2 -- the eps = 1/N conventions: PASS

Numeric identities at N = 1e5, v = 1/5 (`theta = 0.08839305`, `x = t = 8839.30`, `zeta = 2.5e-4`, `eta = 5e-5`), cumulants by exact sums, all differences `<= 4.3e-120` or smaller: `1/x = zeta v/khat_1`; `1/(N-x) = zeta v^2/(1 - v khat_1)`; `t = khat_1/eta`; `sigma^2 = N kt_2 = khat_2/eta`; `c_j = k_j/(j! sigma^j) = eta^{j/2-1} khat_j/(j! khat_2^{j/2})` for j = 3..6 (no `/2` inside the power); `c_1 = delta/sigma = delta sqrt(eta/khat_2)`; the monomial `c_1^2 c_3^2 c_4` (w = 6, # = 3) `= delta^2 zeta^3 v^3 khat_2^{-6} (khat_3/6)^2 (khat_4/24)` with NO `2^{-w/2}` (`8.0939e-17` both ways); `NV = theta(1-theta) - kt_2 = v^2 nhat`; `NV/(kt_2 Q) = nhat/(khat_2 khat_1 (1-theta))`; `1 - G(2000) = SUM (-1)^{k+1} x^{-k}/(k+1)`. Hence the script's `sigma_lo = sqrt(k2lo/eta_max)`, `half_eps = 1/N = eta v in [0, eta_max v_hi]`, `t_min = k1lo/eta_max`, `invx = zeta V inv_k1`, `invNx = zeta V^2 inv_1mt` and the Edgeworth monomial without `2^{-w/2}` are the correct single-spectrum forms. The `half_eps` name is now literally `1/N`, which is what `kappa2_enclosure` uses in `N/sigma^14 = half_eps^6/k2t^7` (N Bernoulli variables, correct here).

Binomial cut `r1`: first omitted terms `(1/x)^6/(7 theta v)` and `(1/(N-x))^4/(5 (1-theta) v)` of alternating series with decreasing terms (`x >= 1001`); at the test point they are `1.7e-23` and `1.6e-20`. The constants `1 - 0.11` and `0.89` need `theta <= 0.11` on the region: `theta(0.24, N) = 0.10379926` (N = 1001), `0.10371445` (8170), `0.10370355` (1e5), `0.10370259` (1e7), and in general `theta(0.24, N) <= theta_inf(0.24) + (0.24/1.24)/N <= 0.103896` for N >= 1001 (`theta_inf = 1 - log 1.24/0.24 = 0.10370258`). So `1 - theta >= 0.896 >= 0.89` and `3 + 3 eps - theta in [2.896, 3.00037]`, `|.| <= FAC = 3.01`; the lower end is positive, so `(3 + 3 eps - theta) N g - 1 > 0` does imply `N g > 1/(3 + 3 eps - theta)`. PASS.

### Item 3 -- the target and the cancellation: PASS

Own series (I = 12, 40 and 150 -- identical low-order coefficients): `nhat(0) = 1/12`; `N g = 1/3 + v/6 + (-11/18 zeta - 19/540) v^2 + ...` (constant `1/3` zeta-free); `(3 + 3 zeta v^2 - v khat_1) N g - 1` has constant term `0` exactly; `S(0, zeta) = 1/3`; `[v^1] S = -5/6 zeta - 7/90`; `[v^2] S = -47/60 zeta + 23/540`. All four match the script's prints exactly (`1/3`, `0`, `1/3`, `(-5/6)*x + (-7/90)`). The cancellation, piece by piece: `piece2` (Edgeworth) `v^0: -2 zeta`, `v^1: 8/3 zeta^2 - 8/3 zeta`; `piece3` (`(1-G(x))/theta`) `v^0: +2 zeta`, `v^1: -8/3 zeta^2 + 8/3 zeta`; `piece1 = 1/3 + v/6 + ...`; `piece4 = zeta v^2/2 + ...`. Why `S(0) = 1/3`: with `N g = 1/3 + v/6 + O(v^2)` and `theta = v/2 + O(v^2)`, `(3 - theta) N g - 1 = v(1/2 - 1/6) = v/3`; equivalently `N g - 1/(3 - theta) ~ (2/9) theta`, the measured margin `2(j-1)/(9n)` of `results/sibuya_margin_probe_2026-09-02.txt`. Own weight-2 Edgeworth part `12 c_4 - 36 c_3^2 + 6 c_1 c_3` exact; monomial counts `{0:1, 2:3, 4:10, 6:26, 8:60, 10:127}`.

The target rederived from the paper's text (scratchpad `sibuya1988.txt`, OCR of the original, p. 699): "(3.4): `(m-1)(2n+m)/(n-m+1) [n,m]/[n,m-1]`, m = 2, ..., n-1, is strictly decreasing and the same for m = n-1 and n". With `e_j = c(n, n-j)`, `m = n - j`, `p_j = e_j/C(N, j)`, `N = n - 1`: `s_m > s_{m+1}` is `e_j^2/(e_{j-1}e_{j+1}) > m(2n+m+1)(n-m+1)/((m-1)(2n+m)(n-m))`, and `C(N,j-1)C(N,j+1)/C(N,j)^2 = j(N-j)/((j+1)(N-j+1)) = (n-m)(m-1)/((n-m+1)m)`, whose product is exactly `1 + 1/(2n+m) = 1 + 1/(3n-j)` -- verified as an fmpq identity for 50 random (n, m); "the same for m = n-1 and n" is the equality at j = 1. The chain: `R = e^g`, `e^g - 1 >= g`, `g >= 1/(3n-j) <=> N g >= N/(3N+3-j) = 1/(3 + 3/N - j/N)` (fmpq identity, 50 random (N, j)). This closes the "transcription not validated" limitation recorded by the dense-(a') and sparse validators. PASS.

### Item 4 -- I_TRUNC = 150, R_DISC = 0.3, the Cauchy tail: PASS

- Arithmetic: `M_S (0.24/0.3)^{148}/(1 - 0.8) = 3.29401e12 * 4.543e-15 / 0.2 = 0.07482`. With the original `I = 60` (d = 57) the tail would be `3.9e7` and with d = 137 it would be `0.70`: the increase to I = 150 is necessary, and the band-1 result `0.2398` is `S_poly - 0.0748 - losses` (own polynomial part at v = 0.24: `0.316781` at zeta = 0, `0.315793` at zeta = 0.004; `0.3158 - 0.0748 = 0.2410`, and the script's box evaluation over `[0.12, 0.24] x [0, 0.004]` reads `0.2398` -- consistent).
- Necessary condition on the own exact coefficients (I = 150): `sup_{|zeta|<=0.26} |S_n| 0.3^n <= M_S` for every n <= 150, with maximum ratio `1.01e-13` (n = 0). The own coefficients grow (`sup|S_n| = 0.33, 0.29, 0.25, 0.25, 0.83, 54, 5666, 4.8e5` at n = 0, 1, 2, 5, 10, 20, 30, 40) but `0.3^n` wins by 13 orders: `M_S` is a valid (very conservative) sup bound as far as the coefficients can tell. At the real points the series converges fast: `|S_n| 0.24^n ~ 1e-91` at n = 145..150, and the degree-40 and degree-150 polynomial parts agree at v = 0.24 to `1e-13`.
- Base object with the new Faulhaber factor: own `M(kt_1)` on `|v| <= 0.3, |zeta| <= 0.26` = `0.193945` (script: `0.19395`); tail 2 (terms with i <= 150, i + 2k > 150) `2.5e-65`; khat_1 inverse: sup of the rest `0.14648`, gap `1/2 - 0.14648 = 0.35352 > 0`. The eps-terms of the Faulhaber polynomials (Bernoulli growth) are harmless on this polydisc: the largest kept monomials carry `0.26^k 0.3^{i+2k}` against `C(i+1,k)|B_k|/(i+1)`.
- The rules of `Ser` (product, inverse with gap, `divV`, exact degree `d`, real-point perturbation `P`) are unchanged from the version validated in Addendum 2 of `VALIDATION_DENSE_CERTIFICATE_B_2026-09-02.md`; the only new ingredient is the base object `G_j = SUM_{i<=150} coef_i v^i sigma_i(zeta v^2)` with a first-order eps term, which is still a polynomial in (v, zeta), so nothing in that proof changes. Exact-degree bookkeeping traced: kt 150 -> khat 149 -> nhat 148 -> N g 148 -> `(3 + 3 eps - theta) N g - 1` 148 -> S 147, matching the printed `d = 147`.
- Not reproduced (harmless): the comment "0.4 fails the inverse gap of khat_1: sup of the rest 0.53 > 1/2". With the current base object the own rest at r = 0.4 is `0.2276` (gap `0.272 > 0`); the comment probably refers to an earlier bookkeeping (tail 1 in M). The certificate uses r = 0.3, where the gap is `0.3535`. Documentation only.

### Item 5 -- the corner v -> 0: PASS

- `q_k = (kv/N)/(1 + kv/N) <= v/(1+v) <= q_max = 0.24/1.24 = 6/31 = 0.193548` (k <= N, `u/(1+u)` increasing). `e^{3/2} - 1 = 3.4816891`, `q_max (e^{3/2} - 1) = 0.6738753 < 1`, so `K(s) = log(1 + q(e^s - 1))` is analytic on `|s| <= 3/2`, `|K| <= |w|/(1-|w|)`, and Cauchy gives `|kappa_j(q)| <= j! (2/3)^j q Kfac`, `Kfac = 3.4816891/(1 - 0.6738753) = 10.675944`. Summing over k (`SUM q_k = k_1 = t`, single spectrum: no factor 2, and the ratio `t/sigma^2 = khat_1/khat_2` is the same either way): `|c_j| <= Kfac (2/3)^j (khat_1/khat_2) sigma^{2-j}` -- `K_j` sigma-free (`K_14 <= 0.0366`, `K_15 <= 0.0244`). The DA `geo_K` branch computes `Kg sigma_lo^2 (1/3)^{14}/(1 - 1/3)` with `L = sigma_lo/2`; this is the term of weight 12 identified in the earlier validation, so every inner piece is `C_w sigma^{-w}` with w >= 11 and sigma-free `C_w`.
- Monotonicity in the new variable: `sigma^2 = khat_2/(zeta v) >= k2lo/(ZMAX v)`, so `sigma^{-11} ~ v^{5.5}` and the loss `R/(khat_2 v^2) ~ v^{3.5}`, decreasing on `(0, v_last]`; the outer pieces `exp(-0.4583 L^2)`, `exp(-0.4596 sigma^2)` with `L = sigma/2` decay like `e^{-c/v}`. Same exponent 3.5 as in the original (there `V = v^2` and `sigma^2 ~ 1/V`).
- Corner remainder re-run with the validator's OWN inputs (own series on `[0, v_last] x [0, 0.26]`, `v_last = 1.8311e-6`: `khat_1 in 0.5 +/- 6e-7`, `k2lo = 0.499999`, `eta_max = 4.761e-7`, `sigma_lo = 1024.82`, `L = 512.41`, `Kg = 10.676`) through the shared engine: `R = 3(d_0 + d_1 + d_2) = 3.237e-27` (script `3.24e-27`), `loss = R/(k2lo v_lo^2) * 3.01 = 2.325e-14` (script `2.32e-14`).
- Bookkeeping notes (harmless): (i) the binomial cut `r1` is not subtracted in the corner; it is `<= 7.7e-25` there and decreases like `v^4`; (ii) the verdict tests `loss_corner < 0.1` rather than `S_corner.lower() - loss_corner > 0`; the printed `S >= 0.3333` on `(0, v_last]` makes the conclusion true, but the mechanical check should be the difference.

### Item 6 -- coverage: PASS

Own derivation without the series: (a) `theta = (1/N) SUM q_k <= (1/N) SUM kv/N = v(N+1)/(2N)`, so `khat_1 = theta/v <= (1 + eps)/2` exactly. (b) `theta(v, N) = 1 - (1/N) SUM f(k/N)`, `f(s) = 1/(1+vs)` decreasing, so `theta(v, N) >= theta_inf(v) = 1 - log(1+v)/v` for EVERY N; `theta_inf(0.24) = 0.1037025849 > 0.1`; and `d theta/dv = kt_2/v > 0`, so `theta <= 0.1 => v < 0.24` (checked numerically: `theta(0.24, N) = 0.10379926, 0.10371225, 0.10370355, 0.10370259, 0.10370259` at N = 1001, 10010, 1e5, 1e7, 1e9, all `>= theta_inf`). (c) `theta <= 0.1` gives `N >= 10 j`, `eps <= 1/(10 j)`; with `b = j^2/N > (j/(j+1))^2` on the strip, `zeta = khat_1^2/b <= ((1+eps)/2)^2 ((j+1)/j)^2 = 0.250550` at j = 1001 `<= 0.26`. (d) `eta = zeta v = 1/(Nv) = khat_1/j <= (1 + 1/10010)/2002 = 4.9955e-4 <= 5.1e-4`. Closed form `k_1 = N - x[psi(N+1+x) - psi(1+x)]` vs direct sum at (1001, 0.24): `|diff| <= 9.2e-116`.

Integer pairs (tilt solved by bisection on the closed form, 300 bits): grid of 120 pairs over j in {1001, 1002, 1100, 1500, 2000, 3000, 5000, 1e4, 3e4, 1e5} and N in {10j, 10j+1, 10j+7, 11j, 20j, 50j, 100j, j^2/3+5, j^2/2, j^2-1, j^2, (j+1)^2-1} (theta <= 0.1, N <= (j+1)^2 - 1): zero misses; extremes `v = 0.23016` (theta = 0.1), `zeta = 0.249998` (N = (j+1)^2 - 1, j = 1e5), `eta = 4.98837e-4` (j = 1001). 1000 random pairs (j log-uniform in [1001, 1e5], N uniform in [10 j, (j+1)^2 - 1]): zero misses, worst `v = 0.1515`, `zeta = 0.2497`, `eta = 4.964e-4`. The strip `j^2 <= N < (j+1)^2` IS inside this script's region (zeta <= 0.2506 at its worst point). PASS for this script; see the assembly note below.

### Item 7 -- end-to-end at integer pairs: PASS

Two exact routes for `e_j(1..N)`: (A) `fmpz_poly` product tree with `mul_low` at degree j + 2; (B) `e_k(1..N) = c(N+1, N+1-k) = SUM_{m=0}^{k-1} <<k, m>> C(N+1+m, 2k)` with the Eulerian numbers of the second kind from `<<n,m>> = (m+1)<<n-1,m>> + (2n-m-1)<<n-1,m-1>>` (own rows `1; 1; 1,2; 1,8,6; 1,22,58,24; ...`) and the binomials updated exactly. Known answer: N = 4 gives `[1, 10, 35, 50, 24]` by both routes; they agree at N = 300 for all k <= 61 and at (N, j) = (1e6, 1001) for `e_{j-1}, e_j, e_{j+1}` (9140-digit integers; tree 182 s, Eulerian 1 s). Negative control: all roots equal gives `R = 1` and fails `n(R-1) >= n/(3n-j)`, as it must. Then `p_j = e_j/C(N,j)`, `R = p_j^2/(p_{j-1}p_{j+1})` exact in fmpq, Sibuya's (3.4) tested exactly in fmpq, `N g = N log R` in arb (400 bits), the tilt from `k_1(v) = j`, `S_true = ((3 + 3/N - theta) N g - 1)/v`, and the validator's own series `S_own` (polynomial part, degree 40):

| (N, j) | theta | b | v | zeta | band | `n(R-1) - n/(3n-j)` (exact > 0) | `N g` | `1/(3+3eps-theta)` | `S_true` | `|S_true - S_own|` | band edge loss | band min (log) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| (10010, 1001) | 0.1 | 100.1 | 0.230138 | 1.886e-3 | 1 | 2.5159e-2 yes | 0.369941512732 | 0.3447919537 | 0.3169456184 | 1.8e-17 | 1.9e-6 | 0.2398 |
| (20000, 2000) | 0.1 | 200 | 0.230150 | 9.44e-4 | 1 | 2.5174e-2 yes | 0.369979466507 | 0.3448097512 | 0.3171662662 | 2.8e-19 | 1.9e-6 | 0.2398 |
| (1e5, 1500) | 0.015 | 22.5 | 0.030610 | 1.067e-2 | 3 | 3.3919e-3 yes | 0.338396331344 | 0.3350050083 | 0.3307115116 | 8.2e-17 | 8.1e-7 | 0.3279 |
| (1e6, 1001) | 1.0e-3 | 1.002 | 2.0047e-3 | 0.24884 | 7 | 2.2249e-4 yes | 0.333666692202 | 0.3334442591 | 0.3327611071 | 2.1e-13 | 1.4e-4 | 0.3321 |
| (1002001, 1001) | 9.99e-4 | 1.000 | 2.0007e-3 | 0.24934 | 7 | 2.2204e-4 yes | 0.333666025878 | 0.3334440375 | 0.3327614187 | 2.1e-13 | 1.4e-4 | 0.3321 |
| (1004003, 1001) strip | 9.97e-4 | 0.998 | 1.9967e-3 | 0.24983 | 7 | 2.2160e-4 yes | 0.333665361881 | 0.3334438167 | 0.3327617292 | 2.1e-13 | 1.4e-4 | 0.3321 |
| (1e7, 3200) | 3.2e-4 | 1.024 | 6.4027e-4 | 0.24393 | 9 | 7.1116e-5 yes | 0.333439969956 | 0.3333688593 | 0.3331533207 | 2.0e-15 | 2.7e-5 | 0.3330 |
| (9e6, 3000) | 3.33e-4 | 1.000 | 6.6696e-4 | 0.24978 | 9 | 7.4079e-5 yes | 0.333444410251 | 0.3333703374 | 0.3331425632 | 2.6e-15 | 2.7e-5 | 0.3330 |
| (2.5e7, 5000) | 2.0e-4 | 1.000 | 4.0011e-4 | 0.24987 | 10 | 4.4446e-5 yes | 0.333399987696 | 0.3333555437 | 0.3332188783 | 3.4e-16 | 2.1e-6 | 0.3332 |
| (25010000, 5000) strip | 2.0e-4 | 0.9996 | 3.9995e-4 | 0.24997 | 10 | 4.4428e-5 yes | 0.333399961040 | 0.3333555348 | 0.3332188907 | 3.4e-16 | 2.1e-6 | 0.3332 |

Every pair satisfies (3.4) exactly; `S_true` exceeds the script's logged band minimum at every pair (band 1: `0.3169 > 0.2398`); the own series reproduces the exact value to `1e-13..1e-19`, i.e. the Edgeworth remainder and cuts are 5-10 orders below the script's per-band allowances (the script's losses are upper bounds, not estimates). The margin `S_true ~ 1/3` at small v and `0.317` at theta = 0.1 is the `2(j-1)/(9n)` law with its first correction (`[v^1] S = -7/90 - 5/6 zeta`).

```python
# validator, item 7 (excerpt): the second exact route, e_k(1..N) = c(N+1, N+1-k)
def e_eulerian(N, k, row):  # row = <<k, 0..k-1>> (Eulerian numbers of the second kind)
    n = N + 1
    C = binom(n, 2 * k)  # C(n, 2k)
    tot = fmpz(0)
    for m in range(k):
        tot += row[m] * C
        M1 = n + m + 1
        C = binom(M1, 2 * k) if M1 - 2 * k <= 0 else C * M1 // (M1 - 2 * k)  # C(n+m+1, 2k), exact
    return tot
```

## Failure map

| region / step | status | size of the defect |
|---|---|---|
| series algebra in v, first-order Faulhaber term, cancellation, `S` coefficients | reproduced exactly (own formulation) | -- |
| eps = 1/N conventions, `FAC = 3.01`, `0.89` | verified against exact sums and `theta(0.24) <= 0.1039` | -- |
| Cauchy tails, `M_S`, `I = 150`, `d = 147` | own coefficients consistent with `M_S` (ratio `1e-13`); tail 0.0748 accounts for the band-1 drop | comment "r = 0.4 fails (0.53)" not reproduced (own 0.228); documentation |
| bands `v in [9.2e-7, 0.24]` | own S polynomial and 10 exact pairs consistent with every logged band minimum | -- |
| corner `v < 9.2e-7` (all N, t up to infinity) | own-input remainder `3.237e-27`, loss `2.3e-14`, K_j sigma-free with `Kfac = 10.676` | `r1` not subtracted in the corner (`<= 7.7e-25`); verdict tests `loss < 0.1` instead of `S - loss > 0`; both harmless, bookkeeping |
| coverage of the integer region incl. the strip `j^2 <= N < (j+1)^2` | verified independently, 1120 pairs, zero misses | -- |
| docstring header (VMAX 0.23, ZMAX 0.25, ETAMAX 2.5e-4, t >= 2001) and the odd-square comment block in `cumulant_series` | stale relative to the code (0.24, 0.26, 5.1e-4, t >= 1001) | text only |
| assembly (`sibuya_theorem.py`, out of scope here) | its coverage lemma states `j^2/N >= 1` for C3' while the sparse piece covers `N >= (j+1)^2` (sparse validator's finding); this script DOES cover the strip (item 6), but the lemma text must say `b >= (j/(j+1))^2` and use `hi_all^2 ((j+1)/j)^2 <= ZMAX` | theorem-level documentation, not a defect of this script |

## Claim-state recommendation

Statement "`(3 + 3/N - theta) N g > 1` on `{v <= 0.24, zeta <= 0.26, zeta v <= 5.1e-4}` for the spectrum `{1..N}`, hence Sibuya's (3.4) for every integer pair with `theta <= 0.1`, `j >= 1001`, `N <= (j+1)^2 - 1`": **independently-validated** for `sibuya_dense_b.py` at blob `8ca3fed5` (commit `bb46ccee`) with log `results/sibuya_dense_b_2026-09-03.txt`. All seven items PASS by separate formulation; the shared dependencies are python-flint and, in item 5 only, the previously validated `kappa2_enclosure` engine fed with the validator's own inputs. The three bookkeeping remarks (corner `r1`, corner verdict test, stale docstring/comments) do not change any number or any step of the proof; fixing them needs no re-validation. Promotion itself goes through `allowed_claim_promotion`; the human approves any public claim. The assembled theorem is NOT covered by this report (the strip wording in `sibuya_theorem.py` and the top regime `theta > 0.9` are open at the assembly level).

```yaml
# templates/validation.yaml schema
id: VAL-SIBUYA-DENSE-B-2026-09-03
experiment_id: sibuya_dense_b (results/sibuya_dense_b_2026-09-03.txt)
validator: independent-validator
independence_level: separate-formulation   # shared: python-flint 0.9.0; dense_certificate_a.kappa2_enclosure in item 5 only, own inputs
checks:
  known_answer: pass          # Stirling row N = 4; two exact routes agree at (1e6, 1001); series vs exact sums 1e-107..1e-118
  negative_control: pass      # all-roots-equal spectrum fails the target; I = 60 would leave a Cauchy tail of 3.9e7 (the port needed I = 150)
  hidden_points: pass         # 10 exact integer pairs incl. the strip N = (j+1)^2 - 1 and the band-1 point (10010, 1001)
  precision_convergence: pass # own S at degree 40 and 150 agree to 1e-13 at v = 0.24; exact pairs reproduced to 1e-13..1e-19
  boundary_map: pass          # theta(0.24, N) >= 0.10370 > 0.1 for every N; zeta <= 0.2506 on the strip; eta <= 4.9955e-4; corner remainder 3.24e-27 with sigma-free K_j
  signature: pass             # constant 1/3, S(0) = 1/3, [v^1] S = -5/6 zeta - 7/90, kt_1, kt_2 leading terms, M(kt_1) = 0.19394, corner 3.24e-27 / 2.32e-14 all reproduced
decision: pass
allowed_claim_state: independently-validated
blockers: []
notes:
  - "bookkeeping (no number changes): corner verdict should test S_corner - loss > 0 and subtract r1; docstring constants and the odd-square comment in cumulant_series are stale"
  - "assembly-level (sibuya_theorem.py): state the strip j^2 <= N < (j+1)^2 explicitly (b >= (j/(j+1))^2, zeta <= 0.2506 <= ZMAX); top regime theta > 0.9 open"
```
