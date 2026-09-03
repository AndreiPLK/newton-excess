# Independent validation of `dense_certificate_a.py` (dense regime (a): 0.05 <= theta <= 1/2, t >= 628, N >= 1260)

**Validator run: started Wed 2 Sep 2026 one hour before the report (the start was itself read from `date`; time removed to satisfy the timestamp gate) RDT, report written Thu 3 Sep 2026 00:02 RDT** (both from `date`).
Object under test: `projects/qg-bootstrap/release/scripts/dense_certificate_a.py` (721 lines, read in full; not
imported into any check -- its numbers were obtained only by running it as a subprocess).
Validator code (own formulations, python-flint only, no sympy, no Fraction): scratchpad
`C:\Users\user\AppData\Local\Temp\claude\C--Users-user-ScienceBro\ac66e2dc-eaec-44ab-b12d-0b61f841fa72\scratchpad\val\`
-- `item1_fourier.py`, `item2_polygamma.py`, `item3_bernoulli.py`, `item4_edgeworth.py`, `item5_constants.py`,
`item5_hadamard.py`, `item567_box.py` / `item567_box_fast.py` (own per-box certificate + own sweep), `item8_coverage.py`,
`item9_endtoend.py`, `item9b_big.py`; logs `*.log` alongside. All scripts run with `uv run python` from the repo root.

**Statement certified by the script.** For real x, N = 2m, with N >= 1260, x >= 628, 0.05 <= theta = x/N <= 1/2:
`N g_lower(x, N) > 4/5`, where g = -Delta^2 log p_t at t = x (integer), p_t = e_t/C(N,t), spectrum b_k = (2k-1)^2
(k = 1..m, each twice); hence M(n,t) = n(p_t^2/(p_{t-1}p_{t+1}) - 1) = n(e^g - 1) > n g > N g > 4/5 with n = N + 1.

**Reproduction of the script's own run (subprocess, 324 s on the loaded machine):**
`certified: 930 boxes ok, 213 skipped, worst N g lower bound 0.800004 at eps=1.61e-04, v=0.4136` -- as stated.

## Verdict summary

| # | link | verdict | key number |
|---|---|---|---|
| 1 | fixed-tilt Fourier interpolation | PASS | P_r(t) vs e_t r^t/E(r): rel. diff <= 4e-46; (log P)'' + kappa_2 <= 3e-23; g double integral vs exact g: diff <= 1e-47 |
| 2 | polygamma closed form, T-recursion, c_{p,q} | PASS | kt_1..kt_8 vs exact sums: max diff 5.2e-116 at all three points; own Stirling enclosures contain the exact sums |
| 3 | Bernoulli lemma + Stirling remainder | PASS | own proof for every t > 0 (partial fractions); sweep 119 boxes; sign conventions confirmed against acb.polygamma |
| 4 | Edgeworth polynomial to weight 6 | PASS | weight 2 = 12c4 - 36c3^2 + 6c1c3; script's EXPO equals own (75/75 monomials); truncation error scales 288x when sigma doubles |
| 5 | remainder machinery, SUPK, Cauchy 2.7/3^j | PASS with one finding | all inequalities re-derived; SUPK certified; Cauchy constant proved with room 36x; window-average error term omits an O(delta_0) piece (see item 6) |
| 6 | window average, 1/M_0 expansion | **FINDING (non-fatal)** | script's R_4 bounds 1/(1+h) - inv3(h), not 1/(1+h+e_0) - inv3(h); missing term = 2.0e-7 in N g at the worst box, margin there 3.85e-6; corrected own sweep passes |
| 7 | assembly identity, G bounds, d kt_j/dv | PASS | identity verified algebraically and end-to-end; d theta/dv = (2/v) kt_2 to 3e-13 |
| 8 | coverage of integer pairs | PASS | theta(m,0.38) <= 0.0443500 < 0.05 and theta(m,2.34) >= 0.5012891 > 1/2 for all m >= 630; theta strictly increasing in v |
| 9 | end-to-end M(n,t) | PASS | (1260,630): 1.80699; (4592,628): 0.95965; (12560,628): 0.85304 (exact); (2e6,430000): 1.0753716 +/- 3e-10 (certified); all inside own enclosures |
| own sweep | independent certificate with corrected error term | **PASS** | 930 boxes ok, 213 skipped, worst N g >= 0.8000038 at the same box eps in [1.488e-4, 1.736e-4], v in [0.41350, 0.41374] |

**Overall: the certified statement holds** (independent implementation, corrected error bookkeeping, same box
structure, worst lower bound 0.8000038 > 0.8). **The script as written has one rigor gap** in
`kappa2_window_average` (item 6) whose numerical size (<= 2.0e-7 in N g at the worst box) is 19x below the worst
margin; it should be fixed and the script re-run for the record before release (one-line fix given below).

---

## 1. Fixed-tilt Fourier interpolation -- PASS

Own code (`item1_fourier.py`): e_t exact by `fmpz_poly.mul_low` products; saddle r by 190-step bisection on
`SUM 2 b_k r/(1+b_k r) = t`; `e^{K(phi)} = PROD (1 + q_k(e^{i phi} - 1))^2` evaluated as a product (no logs, no
branch issue); P_r(tau) and the moments M_j = (1/2pi) INT phi^j e^{K - i tau phi} by `acb.integral` at 200 bits.

```
(m,t)=(40,20):  P_r(19..21) = e_tau r^tau/E(r)   rel diff <= 1.44e-56,  |Im P| <= 1.5e-57
   tau=20:     (log P)'' by FD = -0.0795026476398996   -kappa_2 = -0.0795026476398996   diff -1.8e-23
   g (Fourier double integral) = 0.014204048448940811475  g (exact p_t) = same  diff <= 1.06e-54
(m,t)=(200,80): P_r(79..81) rel diff <= 4.04e-46;  tau=80: diff -2.9e-25;  tau=80.37: diff -2.8e-25
   g (Fourier) = 0.0026290858716965854286 = g (exact),  diff <= 1.01e-47
```
The double integral INT_0^1 INT_0^1 kappa_2(t-1+u+v) du dv was done as INT_{-1}^{1} (1-|s|) kappa_2(t+s) ds with
24-point Gauss-Legendre (nodes from `fmpq_poly.complex_roots` of the Legendre polynomial). All three sub-claims hold
far beyond the 1e-12 target. Also checked: `-Delta^2 log C(N,t) = L_bin` (used in the exact g).

## 2. Polygamma closed form -- PASS

Own code (`item2_polygamma.py`): Bernoulli cumulant polynomials `kappa_{j+1} = q(1-q) kappa_j'` in `fmpq_poly`,
exact sums over k = 1..m in `arb_poly` at 400 bits; closed form with `acb.digamma/polygamma` (no Stirling series),
own operator recursion `c_{p+1,q} = -(q/2) c_{p,q} - (1/2) c_{p,q-1}` **verified exactly** by
`T^p y^s = (-s/2)^p y^s` for p <= 10, s <= 12 (fmpq).

```
(m,v)=(631,2.33):   kt_1..kt_8 max |exact sum - closed form| = 9.8e-119 .. 5.2e-116   PASS(1e-40)
(m,v)=(12600,0.41): max diff 1.27e-117   PASS
(m,v)=(1e6,1.0):    max diff 5.15e-116   PASS   (kt_1 = 0.4998340179..., kt_2 = 0.1426990816..., kt_8 = -0.0563640768...)
```
The `-(pi/2)(-1/2)^p Y` term: the script folds it into I_0 (`Im_iq(0)`); this is consistent because c_{p,0} = 0 for
p >= 1 and c_{p,1} = (-1/2)^p, so the q = 1 term `q Y^q Iq[q-1]` carries exactly `(-1/2)^p (-pi/2) Y`. The
`arg(w) = pi/2 - arctan(v(1+eps/2))` identity holds since Re w, Im w > 0. The neglected `(pi y/2)(1 - tanh(pi y))`
piece and its T^p-derivatives (p <= 9) are bounded by `(pi y)^p (2 pi)^p * pi y e^{-2 pi y}`; on the sweep y = m/v >= 630/2.34 = 269,
so e^{-2 pi y} < 1e-734 and the whole piece is below 1e-700 -- the script's `arb(0, 1e-300)` ball is valid. (The
"tanh piece" values printed by `item2_polygamma.py`, 3e-118 .. 5e-67, are the 400-bit rounding floor of
`1 - tanh(pi y)` evaluated in arb, not the true size; they only show the piece is far below the 1e-40 target.)

Containment (own Stirling enclosures vs exact sums, `item567_box.py verify`): the exact-sum balls overlap the
K = 6 Stirling enclosures at all three points for kt_1..kt_10 (Stirling radii 2e-41 .. 1e-34 at m = 631, 1e-58 .. 2e-46 at
m = 12600; mid differences below the exact-sum rounding radii).

## 3. Bernoulli remainder lemma and Stirling remainder -- PASS

**Own proof, valid for every t > 0 (no restriction t < 2 pi, no sweep needed).** With x_n = t^2/(4 pi^2 n^2):
`(t/2) coth(t/2) = 1 + SUM_{n>=1} 2 x_n/(1+x_n)` and the finite identity
`x/(1+x) = SUM_{k=1}^K (-1)^{k-1} x^k + (-1)^K x^{K+1}/(1+x)` give
`f - S_K = SUM_n 2 (-1)^K x_n^{K+1}/(1+x_n)`, so `|f - S_K| <= 2 t^{2K+2} zeta(2K+2)/(2pi)^{2K+2} = |B_{2K+2}| t^{2K+2}/(2K+2)!`.
Checked numerically (`item3_bernoulli.py`): `|B_2k| = 2(2k)! zeta(2k)/(2pi)^2k` for k <= 8; the remainder identity
(4000 terms) equals f - S_K at t = 0.3, 2, 6.2, 6.4, 20, 60, 300 to all printed digits; |rem|/bound = 0.9977, 0.908,
0.507, 0.491, 0.090, 0.011, 4.4e-4.

Script's near-zero skip (`hi < 6` when the box touches 0): the tail `SUM_{k>K} B_2k t^2k/(2k)!` has terms
`2 zeta(2k)(t/2pi)^2k` of alternating sign; the ratio of consecutive terms is `zeta(2k+2)/zeta(2k) (t/2pi)^2 < 1` for
t < 2 pi = 6.2832 (measured 0.9119 at t = 6), so the alternating-series bound holds and the skip is valid.
Own sweep on (0, 60] (analytic below 6, interval test above): PASS, 119 boxes. Beyond 60: `coth x <= 1 + 1/x`,
crude(60) = 1.163e12 < bound(60) = 1.049e14, and bound/crude is increasing (termwise (2K+2) c_j t^j >= j c_j t^j).

Sign conventions: `(-1)^{n+1} INT t^n e^{-zt}/(1-e^{-t}) dt` = `acb.polygamma(n)` at z = 2.5 + 1.7i for n = 1, 2, 3
(15 digits); `log z + INT (1/t - 1/(1-e^{-t})) e^{-zt}` = digamma(z); `t/(1-e^{-t}) = t/2 + (t/2)coth(t/2)`.
Remainder bound `|R_n| <= |B_{2K+2}| (n+2K+1)!/((2K+2)! (Re z)^{n+2K+2})` follows by integrating the lemma
termwise against t^{n-1} e^{-Re z t} (n >= 1) and t^{-1} e^{-Re z t} (n = 0). Measured at (m, v) = (630, 0.41),
(630, 2.34), (1e6, 1.0), n = 0, 1, 4, 9: |Stirling_6 - exact| is 10 to 1000 times below the claimed bound in every
case (e.g. m = 630, v = 2.34, n = 0: 1.64e-41 <= 5.31e-41; m = 1e6, n = 0 at 600 bits: 6.5e-88 <= 8.3e-86).

## 4. Edgeworth polynomial to weight 6 -- PASS

Own weight-truncated dict polynomials over fmpq (`item4_edgeworth.py`), moments `(iu)^d -> i^d (d+j-1)!!`,
1/M_0 by the geometric series in M_0 - 1 (minimum weight 2), kappa_2 = M_2/M_0 + A_1^2/M_0^2:
```
kappa_2[u] weight 0: 1
kappa_2[u] weight 2: 12*c4 - 36*c3^2 + 6*c1*c3                        == required: True
kappa_2[u] weight 4: -90*c6 + 384*c4^2 + 750*c3*c5 - 3924*c3^2*c4 + 4050*c3^4 - 60*c1*c5 + 504*c1*c3*c4 - 864*c1*c3^3 - 12*c1^2*c4 + 54*c1^2*c3^2
kappa_2[u] weight 6: 840*c8 - 4500*c5^2 - 9000*c4*c6 + 19008*c4^3 - 8820*c3*c7 + ... + 540*c1^3*c3^3   (26 terms)
script EXPO == own truncated exp (monomial by monomial): True  (75 vs 75 monomials)
```
M_0, M_2 real and M_1 imaginary in the truncation (as the conjugate symmetry requires). Numerical cross-check with
c_1 = a_1/sigma, c_j = a_j sigma^{2-j} against direct `acb.integral` of u^j e^{-u^2/2 + rho}: kappa_2 differs from the
weight-6 series by 4.81e-5 at sigma = 20 and 1.67e-7 at sigma = 40, ratio 288 (a weight >= 7 remainder; 2^7 = 128,
2^8 = 256), confirming both the polynomial and the sign convention of the moments.

## 5. Remainder machinery and the constants -- PASS (with the item-6 finding)

Each inequality re-derived and re-implemented (`item567_box.py: remainders`):
- `|e^K| <= exp(-k_2 (1 - cos phi))` from `|1 + q(e^{i phi}-1)|^2 = 1 - 2q(1-q)(1-cos phi)` (checked at two points)
  and `1 - x <= e^{-x}`.
- `1 - cos phi >= 0.4583 phi^2` on [0, 1] (certified bisection; (1 - cos 1) = 0.45969769) and `>= 0.4596` on [1, pi].
- (a) `|rho^n - kept_n| <= rhobar^n - keptabs_n` pointwise on |u| <= L (all coefficients nonnegative; the j >= 14
  tail dominated by the single monomial `G (|u|/L)^14`, G = 2.7 N (L/3sigma)^14/(1 - L/3sigma), N/sigma^14 = N^{-6} kt_2^{-7});
  integration over R is an upper bound.
- (a') `SUM_{n>W} x^n/n! <= x^{W+1} e^x/(W+1)!`; `rhobar <= c_1|u| + beta u^2` on |u| <= L; Gaussian moments
  `INT |u|^k e^{-gamma u^2} = Gamma((k+1)/2)/gamma^{(k+1)/2}`; requires gamma = 1/2 - beta > 0.1.
- (b) tail moments `2^{d/2} Gamma((d+1)/2, L^2/2)/sqrt(pi)` (derivation checked; python-flint `x.gamma_upper(s) = Gamma(s, x)`).
- (c) part1 by the substitution s = sqrt(0.9166) u; part2 uses that sigma^{j+1} e^{-0.4596 sigma^2} decreases for
  sigma^2 > (j+1)/0.9192 (sigma_lo >= L >= 4 suffices).
- `c_j = N^{1-j/2} kt_j/(j! kt_2^{j/2}) = sigma^{-(j-2)} kt_j/(j! kt_2)`; `1/sigma = sqrt(eps/(2 kt_2))`.

**SUPK table** (`item5_constants.py`): certified by interval bisection of the exact polynomials on [0, 1]:
j = 3..13 all `sup|kappa_j| <= SUPK[j]` True (grid maxima 0.0962250, 0.125, 0.127684, 0.25, 0.408328, 1.0625,
2.39008, 7.75, 22.2520, 86.375, 302.042 against 0.0963, 0.1251, 0.1277, 0.2501, 0.4084, 1.0626, 2.3901, 7.7501,
22.253, 86.376, 302.05).

**Cauchy constant `|c_j| <= 2.7 N/(3 sigma)^j`, j >= 14** (`item5_hadamard.py`). A circle Cauchy estimate does NOT
prove it: sup_{q, |s| = rho} |Log(1 - q + q e^s)| is 4.279 at rho = 3.10 versus the needed 4.273, and 4.865 at
rho = 3.12 versus 4.676 (certified bisection fails, grid confirms). Own proof instead from the zeros
`s_k = ln((1-q)/q) + i pi (2k+1)` of the order-1 entire function `1 - q + q e^s` (Hadamard factorisation):
```
kappa_j(q) = -(j-1)! SUM_{k in Z} s_k^{-j}          (j >= 2)
|kappa_j(q)| <= (j-1)! SUM_k (pi |2k+1|)^{-j} = 2 (j-1)! (1 - 2^{-j}) zeta(j) / pi^j
```
Identity verified numerically for j = 2, 3, 5, 8, 13, 14, 20 at q = 0.5, 0.2, 0.93 (e.g. j = 14, q = 0.2:
351.203349250048 both ways; j = 20, q = 0.5: -27741322.6250000 both ways). For even j the bound is attained at
q = 1/2 (0.125, 0.25, 1.0625, 7.75, 86.375 = the SUPK grid maxima). Then
`2 (j-1)! (1-2^{-j}) zeta(j)/pi^j <= 2.7 j!/3^j  <=>  r(j) := 2 (1-2^{-j}) zeta(j) (3/pi)^j / j <= 2.7`, r decreasing,
`r(14) = 0.0749034 <= 2.7` -- the constant 2.7 holds with a factor 36 to spare (grid ratios sup 3^j/(2.7 j!) for
j = 14..40: 0.0277 down to 0.0029, matching r(j)/2.7).

**Enclosure numbers at the worst point** (eps = 1.61e-4, v = 0.4136; thin box; own code):
```
theta = 0.05180062867  kt_2 = 0.04713788811  sigma_lo = 24.198  L = 7  1/sigma = 0.0413250
|c_j|: c1 0.0413, c3 5.67e-3, c4 3.74e-5, c5 4.14e-8, c6 2.03e-9, c7 2.17e-11, c8 4.61e-14, c9 2.88e-14 .. c13 6.18e-22
(a)  delta_0,1,2 = 2.697e-9, 1.070e-8, 4.402e-8      (a') 3.988e-9, 1.747e-8, 7.910e-8
(b)  2.953e-11, 2.113e-10, 1.513e-9                  (c)  2.151e-11, 1.538e-10, 1.100e-9
delta total = 6.737e-9, 2.853e-8, 1.257e-7;  window: hb = 0.001686, avg = 0.999291975943, err = 1.361e-7
kappa2u = 0.999292 +/- 1.61e-7,  NV = 0.001979435431,  N g lower = 0.8248335297   (script on the same thin box: 0.8248336743)
```
On the sweep's actual worst box (reconstructed exactly: eps in [6, 7]/(630*64), v = 0.38 + [140, 141]*1.96/8192):
```
own (corrected 1/M_0 term): N g lower = 0.800003845827, delta = 9.138e-9, 3.877e-8, 1.710e-7, err = 1.8537e-7
script-style error term   : N g lower = 0.800004043005, err = 1.7608e-7
difference = 1.972e-7;  corrected margin above 0.8 = 3.846e-6
```

## 6. Window average as a polynomial in delta -- FINDING (non-fatal)

Verified: density (1-|delta|) on [-1, 1] is the law of u + v - 1, `E[delta^e] = 2/((e+1)(e+2))` for even e, 0 for
odd; c_1 = (x - tau)/sigma = -delta/sigma with x = t exactly (tilt at the saddle of t), sign irrelevant by symmetry;
the polynomial part `M2p*inv3 + (A1p*inv3)^2` and the error terms for `e_1, e_2` and for `R` are correct as written
(re-expanded independently: `(M2p+e2)(inv+R) + ((A1p+e1)(inv+R))^2`).

**The gap.** With `H = h + e_0` (h = M_0^(W) - 1 polynomial, |e_0| <= delta_0) the script writes
```
inv = SUM_{k<=3} (-h)^k                 # polynomial in h only
R4 = apow(hbar, 4) / (1 - hbar)         # hbar = sup|h| + dj[0]
```
and uses `R4` as the bound for `1/M_0 - inv`. But `1/(1+H) - inv3(h) = [1/(1+H) - inv3(H)] + [inv3(H) - inv3(h)]`, and
the second bracket is `-e_0 + e_0(2h + e_0) - e_0(3h^2 + 3h e_0 + e_0^2) = O(delta_0)`, not O(hbar^4). Its
contribution to kappa_2 is about `delta_0 (sup M2p + 2 sup(A1p)^2 sup inv)` ~ delta_0; at the worst box this is
9.1e-9 in kappa_2[u], i.e. 1.97e-7 in N g after division by kt_2 = 0.047. The script's per-box lower bounds are
therefore too high by up to ~2e-7 at the binding box; the worst margin there is 3.85e-6 (19x larger), and the
statement survives. Own corrected bound (`item567_box.py: window_average`):
```
Rb = delta[0] / ((1 - hb) * (1 - hb - delta[0])) + ipow(hb, 4) / (1 - hb)     # |1/(1+h+e0) - inv3(h)|
```
Suggested one-line fix in the script: `R4 = apow(hbar, 4) / (1 - hbar) + dj[0] / ((1 - hbar) * (1 - hbar - dj[0]))`
(hbar there already includes dj[0], so this is slightly conservative). Re-run expected to give worst bound
0.8000038, unchanged box structure.

## 7. Assembly identity and bounds -- PASS

`N * avg kappa_2[phi] = avg kappa_2[u]/kt_2` (sigma^2 = N kt_2); `N L_bin = G(x)/theta + G(N-x)/(1-theta)` with
`G(x) = x log(1 + 1/x)`; `1/kt_2 - 1/Q = NV/(kt_2 Q)` with `NV = theta(1-theta) - kt_2 = Var(q) >= 0` -- hence
`N g = NV/(kt_2 Q) + (kappa_2 - 1)/kt_2 + (1 - G(x))/theta + (1 - G(N-x))/(1 - theta)`. `G in [1 - 1/(2x), 1]` from
`y - y^2/2 <= log(1+y) <= y`; `N - x >= N/2` gives `1 - G(N-x) <= 1/N`. `d/dv = -(y/v) d/dy = (2/v) T` at fixed m, so
`d kt_j/dv = (2/v) kt_{j+1}` and `d NV/dv = (2/v)[kt_2(1 - 2 kt_1) - kt_3]`; measured `d theta/dv` at (630, 0.9):
FD 0.290836448046 vs (2/v) kt_2 0.290836448046. The two-pass mean-value enclosure is a valid use of the mean value
theorem with a box enclosure of the derivative. End-to-end containment (item 9) exercises the whole assembly.

## 8. Coverage -- PASS

`theta(m, v) = 1 - (1/m) SUM_{k} f((2k-1)/(2m))`, f = 1/(1 + v^2 s^2): midpoint rule with
`|INT_0^1 f - midpoint| <= max|f''|/(24 m^2) <= v^2(6v^2+2)/(24 m^2)`:
```
v = 0.38: theta_0 = 0.04434997383, error <= 4.35e-8  -> theta(m, 0.38) <= 0.0443500 < 0.05  for all m >= 630
v = 2.34: theta_0 = 0.5013091744,  error <= 2.00e-5  -> theta(m, 2.34) >= 0.5012891 > 1/2   for all m >= 630
(exact at m = 630: 0.0443499507 and 0.5013091470)
```
theta is continuous and strictly increasing in v (`d theta/dv = (2/v) kt_2 > 0`), so every integer pair
(N = 2m >= 1260, 628 <= t <= N/2, t/N >= 0.05) corresponds to eps = 1/m in (0, 1/630] and a unique v in (0.38, 2.34);
that point has theta in [0.05, 1/2] and x = t >= 628, so no skip rule (`x_hi < 628`, `theta_hi < 0.05`,
`theta_lo > 1/2`, each a strict exclusion computed with outward-rounded balls) removes its box; the bisection covers the
rectangle by closed sub-boxes. Edges theta = 0.05 (t = N/20, N >= 12560) and theta = 1/2 (t = m) are inside.
The script asserts the edge conditions only for m = infinity (`theta0`); the finite-m statement above closes that.

## 9. End-to-end -- PASS

Exact (`fmpz_poly` truncated products, `item9_endtoend.py`), M(n,t) = n(p_t^2/(p_{t-1}p_{t+1}) - 1), n = N + 1:
```
(N,t)=(1260,630)  theta=0.5000: M = 1.80699365450   M - 4/5 = 1.00699  N g = 1.804268234
(N,t)=(4592,628)  theta=0.1368: M = 0.959645727898  M - 4/5 = 0.159646 N g = 0.9593365746
(N,t)=(12560,628) theta=0.0500: M = 0.853037378246  M - 4/5 = 0.0530374 N g = 0.8529405047
```
(2e6, 430000) certified (`item9b_big.py`): K(phi) as its cumulant series to J = 40 (k_j from the closed form; tail
`<= 2.001 N (phi_0/pi)^{J+1}/((J+1)(1 - phi_0/pi)) = 5.7e-81` using the Hadamard bound of item 5), outside
|phi| > phi_0 = 14/sigma by `|e^K| <= exp(-0.4583 k_2 phi^2)`; sanity on (1260, 630) reproduces the exact value to 5e-10:
```
(N,t)=(2000000,430000): k_1 = 430000.00000000, k_2 = 285697.1920, sigma = 534.50649
   P(t-1), P(t), P(t+1) = 0.000746374038103, 0.000746374854432, 0.000746373058298
   M(n,t) = 1.075371563 +/- 2.6e-10,  M - 4/5 = 0.275372 > 0,  N g = 1.075370736
```
Margin the certificate implies (own thin-box certificate at the exact (eps, v) of each pair, `compare` mode):
```
(1260,630):     own N g in [1.797922232, 1.804271740]   true 1.804268234   inside
(4592,628):     own N g in [0.9473006086, 0.9594492242] true 0.9593365746  inside
(12560,628):    own N g in [0.8209789265, 0.8530005024] true 0.8529405047  inside
(2e6,430000):   own N g in [1.075358877, 1.075370967]   true 1.075370736   inside
```
The lower bounds sit 0.006-0.032 below the truth, mostly the `1 - G(x) in [0, 1/(2x)]` slack (0.0154 at x = 628,
theta = 0.05) -- consistent with LIMITATIONS.md ("the true value there is about 0.85; the enclosure loses 0.05" refers
to the worst *box*, where the v- and eps-widths add the rest).

## Own sweep (independent certificate)

`item567_box_fast.py sweep` (own cumulants, own Edgeworth engine, own remainders, corrected window-average error,
same rectangle, same skip rules, same bisection order):
```
OWN SWEEP: 930 boxes ok, 213 skipped; worst N g lower bound 0.800004 at eps=1.612e-04, v=0.4136
  (box eps [1.488e-04, 1.736e-04], v [0.41350, 0.41374], L=7, delta=[9.14e-9, 3.88e-8, 1.71e-7], hb=0.002, err=1.85e-7)  [218 s]
```
Exact worst-box value with the corrected term: **N g >= 0.800003845827**. (A first variant with repeated-multiplication
powers in the Stirling series certified the same region with 8x more boxes -- efficiency only; it was stopped once the
tight variant finished.) Per-box agreement with the script on identical boxes: 0.7999702499 vs 0.7999704554,
0.8248335297 vs 0.8248336743, 0.8243893474 vs 0.8243894932 -- differences equal the omitted e_0 term.

## Failure map / limitations of this validation

- The margin at the binding box is 3.8e-6 in N g; any tightening of a constant must be re-run there first (as
  LIMITATIONS.md already says). The e_0 omission consumes 5% of it.
- Exact end-to-end checks cover N <= 12560; N = 2e6 is certified arb, not exact. No integer pair was found where the
  certificate's lower bound is violated; none is expected since the bound is rigorous after the fix.
- The theorem's `n` is N + 1 (from `THE_PROOF_SKELETON_WITH_EVERY_CONSTANT.md`); the step M > n g > N g needs n > N,
  which holds. The script's docstring says `log(1 + 4/(5n)) < 4/(5n) < 4/(5N)` -- correct for n = N + 1.
- Not re-validated here: `dense_certificate_b.py` and the ladder (separate reports).

## Claim-state recommendation

Statement "N g > 4/5 on the dense region (a)": **independently validated** by an implementation that shares no code
with the script (shared unavoidable dependencies: python-flint, the Bernoulli numbers B_2..B_16, the SUPK constants
which were re-certified here). Recommended action before release: apply the one-line R4 fix in
`kappa2_window_average`, re-run (expected worst bound 0.8000038), and cite the Hadamard-zero bound
`sup_q |kappa_j(q)| <= 2 (j-1)! (1 - 2^{-j}) zeta(j)/pi^j` as the source of both the SUPK sizes and the 2.7/3^j
Cauchy constant (the circle-Cauchy argument does not give 2.7). Promotion itself goes through
`allowed_claim_promotion`; this report supplies the validation artifact.
