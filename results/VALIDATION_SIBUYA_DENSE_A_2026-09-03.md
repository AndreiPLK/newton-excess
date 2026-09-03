# Independent validation of `sibuya_dense_a.py` (Sibuya dense regime (a'): 0.1 <= theta <= 1/2, t >= 2001, N >= 4002)

**Validator run: started Thu 3 Sep 2026 01:41 RDT, report written Thu 3 Sep 2026 01:53 RDT** (both read from `date`).
Object under test: `projects/qg-bootstrap/release/scripts/sibuya_dense_a.py` (read in full; its functions were NOT
imported -- its numbers come only from running it as a subprocess). The shared engine `dense_certificate_a.py`
(`kappa2_enclosure`, `kappa2_window_average`, the Edgeworth polynomial, the Bernoulli lemma, SUPK, the 2.7/3^j Cauchy
constant) was validated on 2 Sept (`VALIDATION_DENSE_CERTIFICATE_A_2026-09-02.md`); the R4 fix from that report is
confirmed present in the current source (line 407) and the Sibuya script calls the fixed version. This report
concentrates on what changed. Validator code (own formulations, python-flint only, no sympy, no Fraction, no float in
exact comparisons), scratchpad
`C:\Users\user\AppData\Local\Temp\claude\C--Users-user-ScienceBro\ac66e2dc-eaec-44ab-b12d-0b61f841fa72\scratchpad\val2\`:
`s1_cumulants.py` (item 1), `s2_box.py` (own per-box certificate and own sweep, items 2 and 5; built by `build_s2.py`
from the validator's own 2-Sept engine `val/item567_box_fast.py`, lines 99-327, with new cumulants and conventions),
`s3_exact.py` (exact M(n,t)), `s4_coverage.py` + `s4b_fd.py` (items 3-4); logs `s1.log`, `s2_point.log`, `s2_sweep.log`,
`s2_compare.log`, `s2_compare2.log`, `s3.log`, `s4.log`, `s4b.log`, `script_run.log`.

**Statement certified by the script.** For real x and N with N >= 4002, x >= 2001, 0.1 <= theta = x/N <= 1/2, spectrum
b_k = k (k = 1..N, once each), p_t = e_t/C(N,t): `N g > 1/(3 + 3/N - theta) = N/(3N + 3 - t)` with g = -Delta^2 log p_t
at t = x; hence with n = N + 1: `p_t^2/(p_{t-1}p_{t+1}) - 1 = e^g - 1 >= g > 1/(3n - t)`, Sibuya's (3.4) as transcribed
in `sibuya_theorem.py`.

**Reproduction of the script's own run (subprocess, 136 s on the loaded machine, exit 0):**
`certified: 266 boxes ok, 43 skipped, worst margin N g - 1/(3 + 3 eps - theta) >= 0.000010 at eps=4.69e-05, v=0.2444 (N g >= 0.345496)`
-- identical to the logged `results/sibuya_dense_a_2026-09-03.txt` (which took 95 s).

## Validation record (templates/validation.yaml)

```yaml
id: VAL-SIBUYA-DENSE-A-2026-09-03
experiment_id: sibuya_dense_a (results/sibuya_dense_a_2026-09-03.txt)
validator: independent-validator
independence_level: separate-formulation     # cumulants: own Leibniz operator recursion; box engine: validator's own
                                             # 2-Sept code; shared dependencies documented below
checks:
  known_answer: pass          # N = 4 row of the unsigned Stirling triangle [1, 10, 35, 50, 24]; exact-sum cumulants
  negative_control: pass      # all-roots-equal spectrum (R = 1) fails the target, as it must
  hidden_points: pass         # 7 integer pairs inside the region, exact, all inside own enclosures and above target
  precision_convergence: pass # closed form vs exact sums to 1e-36..1e-60 (Stirling), 1e-174 (direct polygamma)
  boundary_map: pass          # binding box at theta ~ 0.105 (lower theta edge), N in [16008, 32020]; margin 1.0e-5
  signature: pass             # own sweep reproduces 266 ok / 43 skipped and the same worst box
decision: pass
allowed_claim_state: independently-validated     # for the script's statement only; see the out-of-scope note
blockers: []                                     # for this script; one assembly-level note below
```

## Verdict summary

| # | item | verdict | key number |
|---|---|---|---|
| 1 | tilted cumulants closed form, recursion, real Stirling remainder | PASS | own Leibniz recursion vs exact sums: 1e-36..1e-60 (K = 6 Stirling), 5e-174 (direct polygamma); containment at all 6 points; remainder ratio max 1 - O(1e-9) < 1 for n = 0..11 |
| 2 | eps = 1/N conventions | PASS | own certificate with these conventions reproduces the script's box structure and worst bound (0.345496) |
| 3 | target 1/(3 + 3 eps - theta) from (3.4) | PASS | exact fmpq identities; e^g - 1 >= g |
| 4 | coverage, skip rules, x_lo, G bounds | PASS | theta(N, 0.21) <= 0.0923274 < 0.1 and theta(N, 2.6) >= 0.5073331 > 1/2 for every N >= 4002 (own Riemann bound); d theta/dv = kt_2/v to 3.8e-83 |
| 5 | own sweep + exact M(n,t) | PASS | own sweep: 266 ok, 43 skipped, worst margin 1.018e-5 (same box); 7 exact pairs inside own enclosures, min exact margin N g - target = 0.02517 at (20010, 2001) |

**Overall: the certified statement holds** -- independent implementation, same box structure, worst own margin
1.018e-5 > 0 at the same binding box as the script.

---

## 1. Tilted cumulants for q_k = k r/(1 + k r) -- PASS

**Own derivation** (`s1_cumulants.py`, different bookkeeping from the script's C_{j,q}): with x = N/v,
`k_1 = SUM (1 - 1/(1 + k r)) = N - x SUM_{k=1}^N 1/(k + x) = N - x [psi(N+1+x) - psi(1+x)]` (the psi difference is the
finite sum by psi(z+N) - psi(z) = SUM_{i<N} 1/(z+i) at z = 1+x). The tilt derivative: d q_k/d log r = q_k(1 - q_k) and
kappa_{j+1}(q) = q(1-q) kappa_j'(q) for a Bernoulli, so k_{j+1} = r d_r k_j = v d_v k_j at fixed N. Writing
F(x) = -x G(x), G = psi(N+1+x) - psi(1+x), D = v d_v = -x d_x, D^p = SUM_q a_{p,q} x^q d_x^q with
`a_{p+1,q} = -q a_{p,q} - a_{p,q-1}`, and Leibniz `d_x^q (x G) = x G^{(q)} + q G^{(q-1)}`, then with Y = 1/v,
Phi_q(w) = N^q psi^{(q)}(N w), w_1 = 1 + Y + eps, w_2 = Y + eps, Delta_q = Phi_q(w_1) - Phi_q(w_2):
`x^q G^{(q-1)} = N Y^q Delta_{q-1}`, `x^{q+1} G^{(q)} = N Y^{q+1} Delta_q`, hence
`kt_{p+1} = [p = 0] - SUM_q a_{p,q} (q Y^q Delta_{q-1} + Y^{q+1} Delta_q)`. This must agree with the script's
`kt_j = [j = 1] + SUM_q C_{j,q} Y^{q+1} Delta_q`, C_{1,0} = -1, C_{j+1,q} = -(q+1) C_{j,q} - C_{j,q-1}; the script's
recursion is also checked by hand: v d_v Y = -Y and v d_v Delta_q = -Y Delta_{q+1} (from d_v w_i = -Y/v and
d_w Phi_q(w) = N^q psi^{(q+1)}(Nw) N = Phi_{q+1}(w)) give exactly the terms -(q+1) C_{j,q} and -C_{j,q-1}.

**Three routes compared** (600 bits): (A) exact sums of the Bernoulli cumulant polynomials over k = 1..N (arb_poly);
(B) own recursion with flint digamma/polygamma at the real points N w (no Stirling); (C) own recursion with own
K = 6 Stirling series and the remainder |B_14| (n+13)!/14! eps^14/w^{n+14}:
```
(N,v)=(4002,2.5):    theta 0.498984041750  kt_2 0.215416385676  |A-B| <= 3.2e-175  |A-C mid| <= 6.2e-36  C radius <= 1.3e-33  contained: True
(N,v)=(4002,0.21):   theta 0.0923056848904 kt_2 0.0812876366824 |A-B| <= 2.5e-174  |A-C mid| <= 6.5e-50  C radius <= 1.4e-47  contained: True
(N,v)=(5000,0.7):    theta 0.242000816284  kt_2 0.169829283805  |A-B| <= 1.6e-175  |A-C mid| <= 1.8e-44  C radius <= 3.9e-42  contained: True
(N,v)=(20000,0.22):  theta 0.0961369684309 kt_2 0.0841991038203 |A-B| <= 7.7e-175  |A-C mid| <= 2.0e-59  C radius <= 4.3e-57  contained: True
(N,v)=(20010,0.2444):theta 0.105350709320  kt_2 0.0910580133004 |A-B| <= 5.3e-174  |A-C mid| <= 7.7e-59  C radius <= 1.7e-56  contained: True
(N,v)=(100000,1.0):  theta 0.306855319434  kt_2 0.193148430552  |A-B| <= 2.1e-174  |A-C mid| <= 1.2e-60  C radius <= 2.5e-58  contained: True
```
(kt_1..kt_11 at every point; "contained" = the exact-sum balls and the direct-polygamma balls both overlap the
Stirling enclosures for all eleven cumulants.) The script's own self-test numbers (4.9e-50 / 1.7e-41 / 5.3e-65 at its
three points) are consistent with these.

**Stirling remainder for real positive arguments (n >= 1 and n = 0).** Valid: psi^{(n)}(z) = (-1)^{n+1} INT t^n
e^{-zt}/(1 - e^{-t}) dt and 1/(1 - e^{-t}) = 1/2 + (1/t)(t/2)coth(t/2); the Bernoulli lemma
|(t/2)coth(t/2) - SUM_{k<=K} B_2k t^2k/(2k)!| <= |B_{2K+2}| t^{2K+2}/(2K+2)! holds for EVERY t > 0 (own partial-fraction
proof, 2-Sept report item 3; the script also re-certifies it on (0, 60] with 2185 boxes), so integrating termwise against
t^{n-1} e^{-zt} gives |R_n| <= |B_{2K+2}| (n+2K+1)!/((2K+2)! z^{n+2K+2}) for n >= 1, and against t^{-1} e^{-zt} (from
psi(z) = log z + INT (1/t - 1/(1-e^{-t})) e^{-zt} dt) gives |R_0| <= |B_{2K+2}|/((2K+2) z^{2K+2}) = the same formula at
n = 0. Scaling by N^n with z = N w produces exactly the script's eps^{2K+2}/w^{n+2K+2}. Measured against direct polygamma
at N in {4002, 20000}, w in {0.3847, 0.5, 1, 1.25, 5.762} (the extreme w_2 and w_1 of the sweep), n = 0..11: the ratio
|Stirling_6 - polygamma| / bound has maximum `1.0000 - O(1e-9)`, strictly below 1 (the ball's upper bound is < 1). For real
z the remainder equals the first omitted term up to a factor 1 - O(1/z^2), so the bound is sharp but never violated.

**Derivative identity** v d_v kt_j = kt_{j+1} at N = 4002, v = 0.7: central differences agree with the closed form to
3.9e-61 (Stirling route, h = 1e-30, j = 1..7).

## 2. Conventions eps = 1/N -- PASS

Checked line by line against the mathematics (the original had eps = 1/m = 2/N):
- `half_eps = eps` = 1/N: `kappa2_enclosure` uses it only in N/sigma^14 = N^{-6} kt_2^{-7} = half_eps^6/k2t^7 (grep: lines
  471-477 of the shared engine); with N Bernoulli terms this is correct. `N_lo = 1/eps.upper()`, `sigma_lo = sqrt(N_lo kt_2.lower)`
  (sigma^2 = k_2 = N kt_2), `inv_sigma = sqrt(eps/kt_2)`, `c_j = sigma^{2-j} kt_j/(j! kt_2) = eps^{j/2-1} kt_j/(j! kt_2^{j/2})`
  (script: `apow(inv_sigma, j-2) * kt[j-1] / (j! k2t)`) -- all consistent.
- Mean-value passes with `(1/v)`: d kt_j/dv = (1/v) kt_{j+1} (item 1); `kt_b = kt_c + dv (1/v) raw[j+1]` with `1/v` and
  `raw` both evaluated on the box encloses (v - v_c) kt_{j+1}(xi)/xi for every xi in the box; second pass tightens.
- NV = theta(1 - theta) - kt_2 = mean(q)(1 - mean(q)) - mean(q(1-q)) = Var_k(q_k) >= 0;
  d NV/dv = (1/v)[kt_2(1 - 2 kt_1) - kt_3] (chain rule with the identity above) -- as coded.
- `x_hi = theta.upper()/eps_lo` = N_hi theta_hi, `x_lo = max(2001, N_lo theta_lo)`, `1 - G(N-x) <= 1/N_lo` -- item 4.
- The target is evaluated as a ball over the box and the margin ball `Ng - target` is required to have positive lower
  bound: this implies N g(p) > target(p) at every point p of the box.
Own certificate with exactly these conventions (`s2_box.py`) at the script's worst point, thin box eps = 4.69e-5,
v = 0.2444: theta = 0.1053504074, kt_2 = 0.09105777, sigma_lo = 44.06, L = 8, delta = (3.96e-11, 1.60e-10, 6.74e-10),
kappa_2 window average = 0.999805262 +/- 9.0e-10, NV = 0.003193928, N g in 0.37 +/- 2.2e-3 (the +/- 2.2e-3 is the
1 - G(x) in [0, 1/(2x)] slack divided by theta at x ~ 2246), target 0.3454481668, margin 0.0224. The full-box version
in the sweep (item 5) reproduces the script's 0.345496 lower bound.

## 3. Target 1/(3 + 3 eps - theta) from (3.4) -- PASS

Exact fmpq (`s4_coverage.py`): `N/(3N + 3 - t) == 1/(3 + 3/N - t/N)` and `N/(3n - t) == N/(3N + 3 - t)` with n = N + 1 at
(N,t) = (4002, 2001), (20010, 2001), (12345, 4321): True. Chain: g > 1/(3n - t) <=> n g > n/(3n - t) <=> N g > N/(3n - t)
(multiplying by the positive N/n), and R - 1 = e^g - 1 >= g (e^g - 1 - g >= 0, minimum 0 at g = 0, checked on a grid of
601 points in [-3, 3]). So the script's statement implies `p_t^2/(p_{t-1}p_{t+1}) >= 1 + 1/(3n - t)` for the pairs in
its region. **Limitation:** the form of (3.4) was taken from `sibuya_theorem.py` (p_j = e_j(1..n-1)/C(n-1, j),
1 <= j <= n-2); the repository's literature verdict records Sibuya 1988 as NOT ACCESSED, so the transcription itself is
not validated here.

## 4. Coverage -- PASS

theta(N, v) = 1 - (1/N) SUM_{k=1}^N f(k/N), f(s) = 1/(1 + v s) decreasing, so INT_0^1 f - (f(0) - f(1))/N <= (1/N) SUM f(k/N)
<= INT_0^1 f, i.e. `theta_inf(v) <= theta(N, v) <= theta_inf(v) + (v/(1+v))/N`, theta_inf = 1 - log(1+v)/v:
```
v = 0.21: theta_inf = 0.09228400186, slack(N = 4002) = 4.337e-5  -> theta(N, 0.21) <= 0.0923273687 < 0.1  for all N >= 4002
v = 2.6:  theta_inf = 0.5073331364 > 1/2                          -> theta(N, 2.6) >= 0.5073331 > 1/2   for all N
exact:    theta(4002, 0.21) = 0.09230568489, theta(20000, 0.21) = 0.09228834069, theta(4002, 2.6) = 0.5074233565,
          theta(20000, 2.6) = 0.5073511914  -- all inside [theta_inf, theta_inf + slack]
```
The script's interval evaluation over eps in [0, 1/4002] (0.09264 and 0.50714) is consistent: my own box evaluation
([0.092 +/- 6.4e-4], [0.507 +/- 6.2e-4]) contains both the N = infinity and N = 4002 values. Monotonicity:
d theta/dv = kt_2/v with kt_2 = (1/N) SUM q_k(1 - q_k) > 0; central differences at eps = 0 and eps = 1/4002, v in
{0.21, 1, 2.6} agree with kt_2/v to 1e-50 (eps = 0) and 3.8e-83 (N = 4002, direct polygamma at 600 bits; an earlier
Stirling-ball finite difference with h = 1e-40 showed 1.5e-10, an artifact of the 7.6e-47 remainder radius divided by
2h -- diagnosed in `s4b.log`, not a property of the object).
Hence every integer pair (N >= 4002, 2001 <= t <= N/2, t/N >= 0.1) is the point (eps, v) = (1/N, v(t/N)) with
eps in (0, 1/4002] and v in (0.21, 2.6), theta in [0.1, 1/2], x = t >= 2001; the skip rules (`x_hi < 2001`,
`theta_hi < 0.1`, `theta_lo > 1/2`, each computed with outward-rounded balls) are strict exclusions that cannot remove
such a point, and the bisection tiles the rectangle with closed leaf boxes that are all `ok` or `skip`. `x_lo`:
x = N theta >= N_lo theta_lo and x = t >= 2001. G(x) = x log(1 + 1/x) in [1 - 1/(2x), 1] checked at x = 2001, 2246, 1e5,
1e9 (e.g. x = 2001: 1 - G = 2.4979e-4 <= 2.4988e-4); 1 - G(N - x) <= 1/N at (4002, 2001), (20000, 2000), (20010, 2001)
(2.4979e-4 <= 2.4988e-4; 2.78e-5 <= 5.0e-5; 2.78e-5 <= 5.0e-5).

## 5. Own sweep and exact M(n,t) -- PASS

**Own sweep** (`s2_box.py sweep`: own cumulants, own Edgeworth/remainder/window engine with the corrected 1/M_0 term,
same rectangle, same skip rules, same bisection order; 83 s):
```
OWN SWEEP: 266 boxes ok, 43 skipped; worst margin N g - 1/(3 + 3 eps - theta) >= 1.018e-05 at eps=4.685e-05, v=0.2444
   (N g >= 0.345496, theta = 0.105 +/- 6.2e-4; box eps [3.123e-05, 6.247e-05] (N in [16008, 32020]), v [0.24384, 0.24501],
    L = 8, delta = (1.45e-10, 5.76e-10, 2.39e-9), hb = 6.7e-4, Rb = 1.45e-10, err = 2.58e-9)
```
Same counts and the same binding box as the script (266/43, 0.345496, eps 4.69e-5, v 0.2444).

**Exact M(n,t)** (`s3_exact.py`: e_t by binary-splitting `fmpz_poly.mul_low` products of (1 + k x), p_t = e_t/C(N,t),
R = p_t^2/(p_{t-1}p_{t+1}) exact in fmpq, Sibuya's inequality R - 1 >= 1/(3n - t) tested exactly in fmpq, N g = N log R in
arb against the fmpq target). Baseline N = 4: e = [1, 10, 35, 50, 24] (row 5 of the unsigned Stirling triangle) -- PASS.
Negative control: all roots equal ({1,...,1}, N = 4002, t = 2001) gives R = 1 and fails both tests, as it must -- PASS.
```
(N,t)      theta    M(n,t)          N g             target N/(3N+3-t)  N g - target  (3.4) exact   own thin-box enclosure of N g          contains
(4002,2001)  0.5000  0.643586962277  0.643374467867  0.3998800959       0.243494      True          [0.641375799049, 0.643374801934]       yes
(20010,2001) 0.1000  0.370001395846  0.369979485514  0.3448097601       0.0251697     True          [0.364896427316, 0.370005011580]       yes
(20000,2001) 0.1001  0.370021576750  0.369999654081  0.3448156960       0.0251840     True          [0.364919043274, 0.370025191715]       yes
(10000,3000) 0.3000  0.471337037655  0.471278803304  0.3703292227       0.100950      True          [0.469922927023, 0.471319753127]       yes
(10000,2001) 0.2001  0.414957488893  0.414907389767  0.3571173488       0.0577900     True          [0.412207147746, 0.414954685642]       yes
(10000,5000) 0.5000  0.643831329375  0.643746231774  0.3999520058       0.243794      True          [0.642946285051, 0.643746285148]       yes
(20000,4000) 0.2000  0.414930850829  0.414905801623  0.3571237255       0.0577821     True          [0.413554344692, 0.414929345103]       yes
own only:  (100000,20000) theta 0.2: N g in [0.414674233896, 0.414949233899], target 0.357139030653
outside the region (information only): (4002,1000) theta 0.25, t < 2001: N g - target = 0.0776, (3.4) True;
   (4002,3000) theta 0.75: 0.7395, True;  (4002,3900) theta 0.9745: 7.231, True
```
The exact values sit inside every own enclosure (at (4002, 2001) and (10000, 5000) within 3e-7 and 5e-8 of the upper
end); the lower ends are 0.002-0.005 below the truth, which is the 1/(2x) slack of 1 - G(x). The smallest exact margin
inside the region among these points is 0.02517 at theta = 0.1 (t = 2001, N = 20010): the certificate's binding margin
(1.0e-5) is a box-width effect at theta ~ 0.105, N ~ 16000-32000, not a near-failure of the inequality.

## Failure map / limitations

- The binding box (eps in [3.12e-5, 6.25e-5], v in [0.24384, 0.24501]) has margin 1.0e-5 against a thin-box margin of
  0.022 at its centre: 99.9% of the room is consumed by the box widths. Any change of a constant in the shared engine
  must be re-run on this box first.
- The (3.4) transcription (index convention, normalisation by C(n-1, j)) is taken from the repository; Sibuya 1988 is
  recorded as NOT ACCESSED.
- Exact end-to-end checks cover N <= 20010 (t <= 5000); larger N are covered by the certificate only.
- Shared, unavoidable dependencies: python-flint; the Bernoulli numbers B_2..B_16; the spectrum-independent constants
  SUPK[9..13] and 2.7/3^j (re-certified 2 Sept from the Hadamard-zero identity, valid for any Bernoulli q); the
  validator's own box engine is the 2-Sept validator code, not the script's.
- **Out of scope for this script, relevant to the theorem assembly** (`sibuya_theorem.py`): C' is labelled
  "j >= 2001 (so N >= 4002, theta <= 1/2)", but j >= 2001 does not imply j <= N/2. Pairs with j >= 2001 and j > N/2
  (theta > 1/2, e.g. (N, j) = (4002, 3000)) are covered by neither the ladder (j <= 2000) nor the three certificates as
  stated. The exact checks above show large margins there (0.74 at theta = 0.75, 7.2 at theta = 0.97), which is evidence,
  not a proof. This should be resolved (a certificate for theta > 1/2, or an argument why that range is not needed)
  before the assembled theorem is promoted.

## Claim-state recommendation

Statement "N g > 1/(3 + 3/N - theta) on {0.1 <= theta <= 1/2, t >= 2001, N >= 4002} for the spectrum {1..N}":
**independently-validated** (separate formulation of the cumulants, separate implementation of the certificate, same box
structure, own worst margin 1.018e-5 > 0, seven exact integer pairs inside the own enclosures and above the target,
known-answer baseline and negative control passed). Promotion itself goes through `allowed_claim_promotion`; this report
supplies the validation artifact. The assembled Sibuya theorem should NOT be promoted on this report alone: the
theta > 1/2, j >= 2001 range noted above is outside every listed piece.
