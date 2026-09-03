# Independent validation: `release/scripts/sibuya_sparse_certificate.py`

**Time (measured with `date`):** Thu Sep 3 02:02:52 RDT 2026
**Validator:** independent-validator (did not write the certificate; did not import it -- run as a subprocess only, plus one out-of-process call of its `mu_deviation_bounds` for a method-agreement table)
**Object validated:** `projects/qg-bootstrap/release/scripts/sibuya_sparse_certificate.py`. **The file changed on disk during the validation**: HEAD `8be2864` tracks the version with `TAU_MIN = 2000` (sha256 `20dcb92f0b05d29cf070640fbe47bdd8b20411a6779a6ca979878e9758d7d540`); the working copy now has `TAU_MIN = 1000` (sha256 `ef9ae383cdfce9565ebe098c8545d78e5441ce92236d289c8b4219191306de50`, one-line diff, uncommitted). **Every a-dependent number below was computed for both; the verdict is for the working copy (`TAU_MIN = 1000`), the 2000 numbers are given in brackets.**
**Engine:** python-flint (fmpq / fmpq_poly / fmpz / arb), 300 bits (2000 bits for the end-to-end cells). No sympy, no `fractions.Fraction`, no float in any comparison.
**Independence:** separate formulation. The script uses `fmpq_mpoly` in (a, b, E) and the recurrence `j E_j = SUM i L_i E_{j-i}`; the validator uses its own trivariate class (dict `(power of a, power of E) -> fmpq_poly in b`), the Taylor-sum exponential `SUM S^k/k!`, its own `D` with `d_b E = a_2 E`, its own substitution `E -> 1 + a_2 b + b^2 Z`, its own sweep; own Faulhaber (Bernoulli recursion) and own rational functions in `u = 1/N` for the moment deviations; own Newton identities + own Cauchy tail for the exact `p_j`. Shared dependency: python-flint only.
**Validation code (scratchpad `sib/`):** `w1_phat.py` (items 1, 2, 3, 6), `ci_lib.py` + `w2_ci.py` (item 4), `w3_lemmas.py` (items 5, 8), `w4_end2end.py` (items 5, 7); logs `w1_1000.log`, `w2_1000.log`, `w3_1000.log`, `w4_1000.log` (and the `TAU_MIN = 2000` runs `w1.log`..`w4.log`).
**Script under test, own output (subprocess):** `TAU_MIN = 1000`: `[a^0]=[a^1]=[a^2]=0`, `[a^3]/b^2(0) = 2/3`, (iii) `2.8142e-03`, (v) `6.951e-14`, (vi) `4.221e-02`, sweep min `0.291732` (1 box), VERDICT True, exit 0, 3 s. [`TAU_MIN = 2000`: (iii) `1.2037e-03`, (v) `6.696e-14`, (vi) `2.097e-02`, min `0.357287`, True.]

## Claim under test

Real-variable certificate: with `a = 1/tau`, `b = tau^2/N`, spectrum `beta_k = (2k - N - 1)/(N + 1)`, `F(tau) = SUM_i e_i(beta) (tau)_i/(N)_i`, `Phi = -(log F)''`:
**`(3 + 3a^2 b - ab) N Phi > 1` for all real `tau >= TAU_MIN` and all `N` with `b <= 1`.**
Integer consequence (Sibuya 1988, eq. 3.4, `p_j = e_j(1..N)/C(N,j)`, `n = N + 1`): `p_j^2/(p_{j-1}p_{j+1}) > 1 + 1/(3n - j)` for every integer `j >= TAU_MIN + 1` and every `N >= (j + 1)^2`.

## Validation plan (declared before running)

| # | item | pass criterion (deterministic) |
|---|---|---|
| 1 | moments: `mu_i^inf = INT_0^1 (2u-1)^i du`, `a_2 = -1/6`, `\|beta_k\| < 1` hence `\|mu_i\| <= 1` | exact identities |
| 2 | target polynomial in own representation: `[a^0..a^2] = 0`, `[a^3]` divisible by `b^2`, `[a^3]/b^2(0) = 2/3`; the chain (3.4) -> `(3+3a^2b-ab) N Phi >= 1` rederived | exact |
| 3 | `FMIN = 0.80`, `CF = 0.35`: certified sweep of `P~`, `DP~/(ab P~)`, `D^2P~/(a^2 b P~)`; the multiplier `m_j` rederived from `-(log(1+u))''` | certified balls |
| 4 | own `C_i` (i <= 60); exact deviations at 5 values of N `<= C_i`; own bound for i > 60 checked against exact deviations to i = 400; `h` formula | exact / interval |
| 5 | tail lemma: `g`, base rederived; ingredients against exact `E_j` at a region cell | certified |
| 6 | factor 3.01 `>= (3+3a^2b-ab) max P~^2` | certified |
| 7 | end-to-end `n(R-1) - n/(3n-j) > 0` at cells with `j >= 1001`, `(j+1)^2 <= N`; known answers against direct Stirling DP and the probe file | certified balls |
| 8 | own sweep minus corrected (iii), (v), (vi): margin > 0; monotonicity in b | certified |

## Results

### Item 1 -- moments: PASS

`mu_i^inf = 0` (odd), `1/(i+1)` (even), checked against the integral of the polynomial `(2u-1)^i` for i <= 30; `mu_2 = 1/3`, `a_2 = -mu_2/2 = -1/6`. `|beta_k| = |2k - N - 1|/(N + 1) <= (N - 1)/(N + 1) < 1`, so `|mu_i(N)| <= 1` (checked exactly for i <= 170 at N = 1002^2 and i <= 60 at N = 10^6) and `|mu_i^inf| <= 1`. Odd finite-N moments vanish exactly (symmetric spectrum). `mu_2(N) = 1/3 - 2/(3(N+1))` exactly. The normalisation is forced: `1 + beta_k = 2k/(N+1)`, and the geometric factor `((N+1)/2)^j` cancels in `p_j^2/(p_{j-1}p_{j+1})`, so `R = F(j)^2/(F(j-1)F(j+1))` exactly with `F(s) = e_s(1+beta)/C(N,s) = SUM_i e_i(beta) C(s,i)/C(N,i)` (random-subset identity `C(N-i, s-i)/C(N,s) = C(s,i)/C(N,i)`).

### Item 2 -- the target polynomial: PASS (exact)

Own construction: `P~` degrees (a, b, E) = (58, 44, 1), `Phat` degrees (120, 89, 2) -- same as the script. `[a^0] P~ = E` exactly. `[a^0] Phat = [a^1] Phat = [a^2] Phat = 0` identically in Q[b, E].

**The E^2 part does not cancel here** (it did for the 4/5 target): validator's derivation `D(EQ) = 2a_2 ab EQ`, `D^2(EQ) = 2a_2 a^2 b(1 + 2a_2 b) EQ` gives the E^2 coefficient `Q^2 a^2 b [(3+3a^2b-ab)(-2a_2) - 1] = Q^2 a^3 b^2 (a - 1/3)`; confirmed exactly in the polynomial (`E^2 part == Q^2 E^2 a^3 b^2 (a - 1/3): True`). Hence `[a^3]` has Z-degree 2 after the substitution.

After `E -> 1 + a_2 b + b^2 Z`: `[a^3]` lowest b-power 2; `[a^3]/b^2 |_{b=0} = 2/3` exactly; every `[a^k]`, k >= 4, divisible by `b^2`. Own output:

```
[a^3]/b^2 = 2/3 - (2/9) b + (7/216) b^2 - (1/324) b^3 + (5/31104) b^4 - (1/155520) b^5 + (7/33592320) b^6 - (1/176359680) b^7 + ... (to b^16)
          + Z ( b^2/3 - b^3/18 + b^4/72 - b^5/1296 + b^6/31104 - b^7/933120 + ... (to b^17) )  - (1/3) Z^2 b^4
```
(full coefficient list in `w1_1000.log`; the line above is the leading part, the exact polynomial is printed there.)

**Structural fact.** With the exact `Z(b)`, `[a^3]/b^2 = (2/3) e^{2 a_2 b} = (2/3) e^{-b/3}` to 12 digits at b = 0, 0.25, 0.5, 0.75, 1 (`0.666666666667, 0.613362943086, 0.564321149927, 0.519200522048, 0.477687540383`; the residual is the O(b^14) truncation). Since `P~ -> E` as `a -> 0`, `((3+3a^2b-ab) N Phi - 1)/(ab) -> [a^3]/(b^2 E^2) = 2/3` for every b: in Sibuya's units `n(R-1) - n/(3n-j) -> (2/9) j/n`, the measured margin `2(j-1)/(9n)`. Item 7 confirms this at every cell.

**The translation from (3.4), rederived.** `R = F(j)^2/(F(j-1)F(j+1)) = e^g`, `g = -Delta^2 log F(j-1) = INT_0^1 INT_0^1 Phi(j-1+u+v) du dv` (exact, any C^2 interpolant with `F > 0` on `[j-1, j+1]`). `R - 1 = e^g - 1 > g` for `g > 0`. If `Phi(tau) > 1/(3n - tau)` on the window then, `x -> 1/(3n-x)` being convex and `E[j-1+u+v] = j`, Jensen gives `g > 1/(3n - j)`. And `N/(3n - tau) = N/(3N + 3 - tau) = 1/(3 + 3/N - tau/N) = 1/(3 + 3a^2 b - ab)` exactly since `1/N = a^2 b`, `tau/N = ab`. So `(3+3a^2b-ab) N Phi > 1` on `[j-1, j+1]` implies `n(R-1) > n/(3n-j)`, i.e. (3.4) with strict inequality. The j = 1 equality of (3.4) (`R_1 = 3(N+1)/(3N+2) = 1 + 1/(3n-1)`) is reproduced exactly at N = 399 and 4000 by direct Stirling numbers, which fixes the transcription of (3.4) used here (the form itself was read from the full text of Sibuya 1988 in `LITERATURE_VERDICT_STIRLING_NEWTON_GAP_2026-09-02.md`; not re-read by this validator).

Z bracket (Leibniz, ratio `|a_2| b/(k+1) <= 1/18`): `Z in [1/72 - b/1296, 1/72]`, checked with the exact Z at b = 0.001, 0.1, 0.5, 0.9, 1 (e.g. b = 1: `0.0131484 in [0.0131173, 0.0138889]`). `|Z| <= 1/72 = 0.01389`; the script's 0.02 is conservative.

(iii) recomputed: `SUM_{k>=4} |[a^k]/b^2| a^{k-3} <= 2.81422e-03` at a = 1/1000, |Z| <= 0.02 (script `2.8142e-03`; with |Z| <= 1/72: `2.7947e-03`) [a = 1/2000: `1.20369e-03`, script `1.2037e-03`].

Own sweep, 256 boxes in b with the Z box: `min_b [a^3]/b^2 = 0.4775246` (true min at b = 1: `0.4776875`); minus (iii): `0.4747104` [2000: `0.4763209`].

### Item 3 -- `FMIN`, `CF`, the multiplier: PASS on the constants; the multiplier formula has one factor-2 slip (non-blocking)

Certified sweep (256 boxes in b, a as the ball `[0, 1/TAU_MIN]`, `E = exp(a_2 b)` as a ball), `F~ = P~/Q`:

| constant | validator, `TAU_MIN = 1000` | validator, `TAU_MIN = 2000` | script |
|---|---|---|---|
| `F~ >=` | **0.845928** | 0.846343 | `FMIN = 0.80` (justified; `e^{-1/6} = 0.846482`) |
| `\|F~'/F~\| <=` | **0.333806** `ab` | 0.333624 `ab` | `CF = 0.35` (justified) |
| `\|F~''/F~\| <=` | **0.333553** `a^2 b` | 0.333552 `a^2 b` | `CF = 0.35` (justified) |
| `P~ <=` | 1.0000019 | 1.00000065 | -- |

(`DP~` is divisible by `ab` and `D^2P~` by `a^2 b` exactly, so the ratios are polynomial and the sweep is rigorous.)

**The multiplier.** From `Phi(F~ + T) - Phi(F~) = -(log(1+u))''`, `u = T/F~`: `|Delta Phi| <= [ |T''|/F + 2|T'||F'|/F^2 + |T||F''|/F^2 + 2|T|F'^2/F^3 ]/(1-|u|) + |u'|^2/(1-|u|)^2`. With `|T_j^{(k)}| <= (ja)^k W_j` (each `|tau - l| <= tau` for `l <= j - 1 <= 2 tau`), `|F'| <= c_1 ab F`, `|F''| <= c_2 a^2 b F`: `effect_j <= (j^2/b) W_j m_j/F` with
`m_j = 1 + 2 c_1 b/j + (c_2 + 2 c_1^2 b) b/j^2`.
The script's `mult(j) = 1/(1-ja)^2 + 2 CF b/(j(1-ja)) + (CF^2 b + CF) b/j^2` has the `(F'/F)^2 |T|` term with coefficient 1 where the derivative of `T/F` gives 2 (`(T/F)'' = T''/F - 2T'F'/F^2 - TF''/F^2 + 2T F'^2/F^3`). At b = 1: `m_2` script 1.4728 vs derived 1.4988 (1.8% low), `m_3` 1.2926 vs 1.2994; at j = 31 the script's `1/(1-ja)^2` makes it 1.0888 vs derived 1.0232 (conservative). The script also omits `K_j = 1/prod_{l<j}(1 - l/N) <= exp(j(j-1)/(2(N-j)))` on the moment terms (`K_30 - 1 <= 4.4e-4` at N >= 10^6), the factor `1/(1-|u|)` (`|u| <= 2.9e-6`) and the `|u'|^2` term (`8.6e-8` relative to ab). Net effect on (vi): +0.1%; see item 8.

Fragility, not an error: `mult(j)` is singular at `j = TAU_MIN` (`1 - j amax = 0`) and negative beyond, and the loop runs `j < 4000`; it is never reached only because the loop breaks at `eff < 1e-300` near j = 285. No remainder for `j > 285` is bounded by the script (validator's bound: `< 3.2e-458`, see item 8).

### Item 4 -- moment lemma constants `C_i`: PASS

Own `C_i = sup_{N >= TAU_MIN^2} N |mu_i(N) - mu_i^inf|` by the exact rational function `(SUM_k c_k u^{i-k})/(1+u)^i`, `u = 1/N` on `[0, 1/TAU_MIN^2]` (own Faulhaber; power sums checked against direct summation at N = 37, 1000 for i <= 60 -- this known-answer test caught a bug in the validator's own Faulhaber before it ran, as intended):

| i | `C_i` (N >= 10^6) | `C_i` (N >= 4e6) | script's function, N >= 4e6 (subprocess) | `1 - 1/(i+1)` |
|---|---|---|---|---|
| 2 | 0.666668668389 | 0.666667166674 | 0.666667166674 | 0.666667 |
| 4 | -- | 0.800001269206 | 0.800001269206 | 0.800000 |
| 10 | -- | 0.909094623290 | 0.909094623290 | 0.909091 |
| 20 | 0.952412384562 | 0.952388811857 | 0.952388811857 | 0.952381 |
| 30 | 0.967790003866 | 0.967753953300 | 0.967753953300 | 0.967742 |
| 60 | 0.983704602346 | 0.983631069772 | 0.983631069772 | 0.983607 |

Two different formulations, twelve-digit agreement. Odd `C_i = 0` exactly. Exact deviations at N = 10^6, 1002^2, 4e6, 10^7, 10^8: `N|mu_i(N) - mu_i^inf| <= C_i` for all i <= 60 (max ratio 0.999997). Domain: the region gives `N >= tau^2/b >= TAU_MIN^2` directly (single spectrum, `m = N`; no factor-2 issue as in the 4/5 case).

**Beyond i = 60.** Validator's bound, by the trapezoid rule with the Peano kernel: `N dev_i = INT f - (f(0)+f(1))/2 - (N+1) E_T` with `f = (2u-1)^i`, `|E_T| <= (h^2/8) INT_0^1 |f''| = i h^2/2` (`INT|f''| = 4i(i-1) INT|2u-1|^{i-2} = 4i`), `h = 1/(N+1)`; hence
**`C_i <= 2 + i/(2(N+1))`** for every i >= 2 and every N. Checked against the exact `N|dev_i|` for all i <= 400 at N = 10^6 (max ratio 0.4986 -- the true value is `~ 1 - 1/(i+1)` for even i). This implies the script's `3 + i^2/TAU_MIN^2` for every i, so the script's constant is valid (its comment's "midpoint error i^2/(6N^2)" is not the argument that proves it, but the claimed bound is weaker than the proved one). `h(rho) = a^2 b^2 rho^2 SUM_i (C_i/i)(ab rho)^{i-2}` is the correct Cauchy bound for `|SUM_i b (delta_i/i)(ab)^{i-2} z^i|` on `|z| = rho` with `|delta_i| <= C_i a^2 b`; the script's `x**(i-2)/2**(i-2)` is `(ab rho)^{i-2}`. The script stops the i-sum at 399 without a remainder; validator's remainder `<= (2/61 + 1/(2(N+1))) x^59/(1-x) < 1e-95` at `x = ab rho <= 0.0245`.

### Item 5 -- tail lemma: PASS

Derivation: `|E_j - E_j^0| <= rho^{-j} e^{|L_2| rho^2} (e^{g} - 1)` with `|L_2| rho^2 = j/4` on `rho_j = sqrt(j/(2 mu_2 b))`, so `<= (2 mu_2 b/j)^{j/2} e^{j/4}(e^g - 1) <= (2e mu_2 b/j)^{j/2}(e^g - 1)`; `g(rho) = SUM_{i>=3} |L_i| rho^i <= b rho^2 SUM_{i>=3}(ab rho)^{i-2}/i <= (a b^2 rho^3/3)/(1 - ab rho)` using `|mu_i| <= 1`. Both as in the script. Convergence: `ab rho_j <= sqrt((1+2a) ab/(2 mu_2)) = 0.03877` for `j <= tau + 2`, `tau >= 1000`, `b <= 1` [0.02740 at 2000] -- inside the disc.

Ingredients with the exact `e_j` (own Newton from exact power sums, j <= 120) at the region cell `N = 1002^2, tau = 1001` (b = 0.999, a = 1/1001):

| j | `\|E_j(N) - E_j^0\|` / Cauchy bound (true moments + `C_2`) | `\|E_j(inf) - E_j^0\|` / Cauchy bound | `\|E_j(N) - E_j(inf)\|` / moment bound |
|---|---|---|---|
| 2 | 2.11e-4 | 0 | 0.3672 (= 1/e: the `e^{j/4} -> e^{j/2}` slack) |
| 4 | 1.04e-4 | 4.9e-5 | 0.1347 |
| 10 | 1.28e-5 | 1.0e-5 | 0.0108 |
| 30 | 1.19e-8 | 1.1e-8 | 5.6e-6 |
| 60 | 3.3e-13 | 3.2e-13 | 1.0e-10 |
| 120 | 2.3e-22 | 2.2e-22 | -- |

All hold for j = 2..120 (odd j are exactly 0 on both sides). `|1 - R_j| < a j(j-1)/2` for 3 <= j <= 400 (j = 2 is an exact equality; ratio 0.9997 at j = 400). `|mu_i(N)| <= 1` (i <= 170), `mu_2(N) <= 1/3` (so the radius `rho_j` with `mu_2^inf` is valid for the true `L_2(N)`). The excluded case `j = t+1`, `tau in [t-1, t)`: `|1 - R_j| <= 2 <= a j(j-1)/2 (>= 500)`. Terms `j > t+1` (`-E_j^0` only): `d/dtau E_j^0 = j a E_j^0` exactly, `|E_j^0| <= (e mu_2 b/j)^{j/2}`. The script's `Tj = base(e^g - 1) + base a j(j-1)/2` bounds `E_j^0(R_j - 1)` and its derivatives by `(j/(tau-j))^k` times the value; for the difference `(tau)_j - tau^j` this needs the sign argument (the two parts of its derivative have opposite signs) and holds where `a j(j-1)/2 <= 1`, i.e. j <= 45 [63]; undocumented. The validator's derivative-safe form uses `(e^{x_j} - 1) + (K_j - 1)`, `x_j = a j(j-1)/2`, or `K_j + 1`, whichever is smaller; it gives (v) = `8.0e-14` vs the script's `6.95e-14` (the script's (v) is 13% low, on a quantity 12 orders below the margin).

### Item 6 -- the factor 3.01: PASS

`(3 + 3a^2 b - ab) P~^2 <= 3.0000135` on the region [3.0000049 at 2000]; the conversion `Phat/(a^3 b^2) = (Q^2 F~^2/(ab)) [(3+3a^2b-ab) N Phi - 1]` means a lemma effect `delta ab` on `N Phi` moves `Phat/(a^3 b^2)` by at most `(3+3a^2b-ab) P~^2 delta <= 3.01 delta`. Justified.

### Item 7 -- end-to-end, Sibuya's inequality at integer cells: PASS

Exact `e_i(beta)` (Newton, i <= 170), `F(s) = SUM_{i<=170} e_i (s)_i/(N)_i` exact in fmpq plus a certified Cauchy tail ball for `170 < i <= s` (own bound `|e_i| <= (e N mu_2/i)^{i/2} exp(N rho^3/(3(1-rho)))`, `rho = sqrt(3i/N)`), 2000 bits:

| N | j | `(j+1)^2/N` | `n(R-1) - n/(3n-j)` | `D n/j` | `D 9n/(2(j-1))` | tail radius |
|---|---|---|---|---|---|---|
| 1002^2 | 1001 | 1.00000 | +2.21600613911976e-4 | 0.222266 | 1.001197 | 1.5e-193 |
| 1100000 | 1001 | 0.91273 | +2.02240632372743e-4 | 0.222243 | 1.001092 | 4.9e-197 |
| 2100000 | 1001 | 0.47810 | +1.05880560064626e-4 | 0.222127 | 1.000572 | 1.6e-221 |
| 1502^2 | 1501 | 1.00000 | +1.47871786683951e-4 | 0.222251 | 1.000798 | 3.8e-194 |
| 2002^2 | 2001 | 1.00000 | +1.10955632196911e-4 | 0.222244 | 1.000599 | 1.9e-194 |
| 4100000 | 2001 | 0.97756 | +1.08464546900760e-4 | 0.222241 | 1.000586 | 2.7e-195 |
| 8100000 | 2001 | 0.49482 | +5.48859403983707e-5 | 0.222177 | 1.000296 | 8.0e-221 |
| 10^7 | 2500 | 0.62550 | +5.55499860929170e-5 | 0.222200 | 1.000300 | 3.5e-212 |
| 3000^2 | 2999 | 1.00000 | +7.40542844024388e-5 | 0.222237 | 1.000400 | 1.0e-194 |
| 10^7 | 3000 | 0.90060 | +6.66684299946578e-5 | 0.222228 | 1.000360 | 1.2e-198 |
| 10^8 | 3000 | 0.09006 | +6.66468423403406e-6 | 0.222156 | 1.000036 | 1.6e-284 |

`D > 0` as a certified ball at all 11 cells; `D n/j -> 2/9` and `D 9n/(2(j-1)) = 1 + O(1/j)`: the margin `2(j-1)/(9n)` is reproduced, b-independent, as item 2 predicts. Known answers: at N = 4000 the ratio `R` from a direct integer DP for `e_j(1..4000)` equals the Newton-on-beta pipeline **exactly** (fmpq equality) at j = 2, 3, 5, 10, 20, 32, and `D` matches `results/sibuya_margin_probe_2026-09-02.txt` (n = 4001: j = 2 `5.55602e-5`, j = 32 `0.00173799`); `(R-1)(3n-j) - 1 = 0.001669` at (n, j) = (400, 2) as in `sibuya_exact_check_2026-09-02.txt`; `p_1^2/(p_0 p_2) = 1 + 1/(3n-1)` exactly at N = 399, 4000.

### Item 8 -- the lemma totals and the corrected margin: PASS

Own recomputation with the derived multipliers, `K_j`, the i-remainder, both j-remainders, `1/(1-|u|)` and `|u'|^2`, at a = 1/1000 (relative to ab):

| b | (v) tail | (vi) moments | `\|u\| <=` | `\|u'\|^2` term | total | script (v) / (vi) |
|---|---|---|---|---|---|---|
| 0.25 | 4.9e-22 | 1.258e-2 | 1.1e-7 | 1.2e-9 | 1.2577e-2 | -- |
| 0.5 | 6.1e-18 | 2.043e-2 | 5.3e-7 | 8.5e-9 | 2.0431e-2 | -- |
| 0.75 | 1.6e-15 | 3.015e-2 | 1.4e-6 | 3.1e-8 | 3.0150e-2 | -- |
| **1** | **8.0e-14** | **4.226e-2** | 2.9e-6 | 8.6e-8 | **4.2260e-2** | 6.951e-14 / 4.221e-2 |

(script's constants FMIN = 0.80, CF = 0.35; with the validator's constants the total is `3.9573e-2`). Remainders: `j in (400, t+1]`: `<= 3.2e-458` (geometric, ratio 0.069, `K_j <= 1.653`, `g_j <= 0.0202 j`); `j >= t+2`: `<= 2.5e-1516`. The totals increase with b and with a (each term is `b^{j/2-1}` times increasing factors, and proportional to a), so b = 1, a = 1/TAU_MIN is the sup. Split (A) (script's: moments for all j) and split (B) (validator's: moments for j <= 30 only, tail with the true moments and `C_2`) give the same totals to 4 digits.

**Assembly:** `subtract = (iii) + 3.0000135 x total = 0.129594`; **margin = `0.4775 - 0.1296 = 0.3479`** (validator's constants: `0.3560`) [2000: `0.4131` / `0.4171`]. The script's `0.291732` (one box, the polynomial evaluated on the whole `b in [0,1]` as a ball) is a valid, coarser lower bound. The dominant correction is (vi), which is genuinely first order in `1/N` here (`mu_2(N) - 1/3 = -2/(3(N+1))`); the true effect measured at the cell is 60x below the Cauchy bound (the `1/e` ratio at j = 2, odd j contribute exactly 0).

## What the certificate implies for integers, and a coverage gap in `sibuya_theorem.py`

For integer `j`, the window is `tau in [j-1, j+1]`, so `a <= 1/(j-1)` and `b <= (j+1)^2/N`. The certificate therefore gives (3.4) for **`j >= TAU_MIN + 1` and `N >= (j+1)^2`**, not for `j^2 <= N`. `sibuya_theorem.py` states C1' as "`j^2/N <= 1`" and C3' as "`j^2/N >= 1`": the strip `j^2 <= N < (j+1)^2` is claimed by neither as written. Arithmetic on the dense-(b) constants (`b >= (j/(j+1))^2 = 0.998` gives `zeta <= 0.2502/0.998 = 0.2507 <= ZMAX = 0.26`) suggests C3' can absorb it, but its coverage lemma must say so; the alternative is to rerun this script with `BMAX = ((TAU_MIN+2)/(TAU_MIN+1))^2` (the margin is flat in b at the 1e-3 level, `C_i` would be recomputed for `N >= TAU_MIN^2/BMAX`). Second, **the change to `TAU_MIN = 1000` moves the boundary of the sparse piece to `j >= 1001`, while the dense pieces and the coverage lemma were set up for `j >= 2001`** (`ETAMAX = 2.6e-4` covers `eta <= khat_1/j` only for j >= 2001; for j = 1001 one needs `5.0e-4`; C2' states `N >= 4002`). Both are theorem-level blockers, not defects of this script.

## Failure map

| boundary | status |
|---|---|
| `tau -> TAU_MIN` (a max) | certified; (iii) = 2.8e-3 and (vi) = 3 x 4.2e-2 are both proportional to a; margin 0.348 at 1000, 0.413 at 2000 |
| `b -> 1` | certified; `[a^3]/b^2 = 0.4777` (min), corrections maximal there (monotone in b) |
| `b -> 0` | certified; `[a^3]/b^2 -> 2/3`, all corrections vanish as powers of b |
| `b > 1` | not certified; the integer consequence needs `b` up to `(j+1)^2/N`, hence `N >= (j+1)^2` (see above) |
| `j = t+1`, `tau < t` | `\|1 - R_j\| <= 2 <= a j(j-1)/2`; holds, undocumented |
| `j > 30` moment terms | covered: the script's (vi) sums all j with `h` over all i (with the trapezoid `C_i` beyond 60) -- the "Gap B" of the 4/5 report does not recur |
| `j = 2` multiplier | factor 2 missing on the `(F'/F)^2` term (1.8% at j = 2, 0.1% on (vi)) |
| `j >= 285` (loop break) | script bounds no remainder; validator: `< 3.2e-458` |
| `j = TAU_MIN` in `mult(j)` | singular (`1/(1-ja)^2`); unreachable at the current parameters because the loop breaks near j = 285 |
| `mu_2(N) < mu_2^inf` | needed for the Cauchy radius; holds (`1/3 - 2/(3(N+1))`) |

## Validation record (schema: templates/validation.yaml)

```yaml
id: VAL-SIBUYA-SPARSE-2026-09-03
experiment_id: sibuya_sparse_certificate.py@sha256:ef9ae383 (working copy, TAU_MIN = 1000, uncommitted; HEAD 8be2864 tracks sha256:20dcb92f with TAU_MIN = 2000)
validator: independent-validator
independence_level: separate-formulation   # shared dependency: python-flint only
checks:
  known_answer: pass          # [a^0..2] = 0 and 2/3 exact; E^2 part = Q^2 E^2 a^3 b^2 (a - 1/3) as derived; C_i to 12 digits by two methods; direct Stirling DP == Newton pipeline exactly at N = 4000; j = 1 equality of (3.4) exact; probe and exact-check values reproduced
  negative_control: pass      # own Faulhaber bug caught by the direct-summation test before use; Z bracket contains the exact Z; odd moments and odd C_i vanish exactly; the deviation is O(1/N) (degree check asserted)
  hidden_points: pass         # 11 integer cells (j = 1001..3000, (j+1)^2/N from 0.09 to 1) and a region cell for the ingredients, none chosen by the implementer
  precision_convergence: pass # 300-2000 bits; exact fmpq wherever possible; certified balls with explicit radii (tail 1.5e-193 and below)
  boundary_map: pass          # table above
  signature: pass             # (n(R-1) - n/(3n-j)) 9n/(2(j-1)) = 1 + O(1/j) at every cell, b-independent, as [a^3]/(b^2 E^2) = 2/3 predicts
decision: pass
allowed_claim_state: independently-validated
blockers:
  - "theorem-level (sibuya_theorem.py), not this script: the integer consequence is j >= TAU_MIN + 1 AND N >= (j+1)^2, not j^2 <= N; the strip j^2 <= N < (j+1)^2 must be assigned to C3' (zeta <= 0.2507 <= ZMAX = 0.26 makes it possible) or BMAX raised to ((TAU_MIN+2)/(TAU_MIN+1))^2"
  - "theorem-level: TAU_MIN = 1000 moves the sparse boundary to j >= 1001 while C2', C3' and the coverage lemma assume j >= 2001 (ETAMAX = 2.6e-4 < 5.0e-4 needed at j = 1001; C2' states N >= 4002); either the ladder reaches 2000 or the dense pieces are re-validated for j >= 1001"
  - "rigor (non-blocking numerically): mult(j) should carry 2 CF^2 b, not CF^2 b, on the (F'/F)^2 term; K_j on the moment terms; the j-remainder beyond the break; the i-remainder beyond 399; and the docstring should replace 'C_i <= i^2 beyond 60' / 'midpoint error' by the proved C_i <= 2 + i/(2(N+1))"
  - "housekeeping: the TAU_MIN = 1000 change is uncommitted; the docstring's ported paragraphs still quote the 4/5 constants (627, 0.669, 0.81, 4/5, 176/175) as if they applied"
```

## Verdict

| item | verdict |
|---|---|
| 1. moments: `mu_i^inf`, `a_2 = -1/6`, `\|mu_i\| <= 1` | **PASS** (exact) |
| 2. target polynomial: `[a^0..a^2] = 0`, `[a^3]/b^2(0) = 2/3`, chain from (3.4) | **PASS** (exact; own min `0.4775`, `[a^3]/b^2 = (2/3)e^{-b/3}`) |
| 3. `FMIN = 0.80`, `CF = 0.35`, multiplier | **PASS** (`F~ >= 0.8459`, `c_1 <= 0.3339`, `c_2 <= 0.3336`); multiplier: factor 2 on one term, K_j, remainders -- corrected in item 8 |
| 4. `C_i`, i > 60 bound, `h` | **PASS** (12-digit agreement; proved `C_i <= 2 + i/(2(N+1))`, which implies the script's `3 + i^2/TAU_MIN^2`) |
| 5. tail lemma `g`, base, ingredients | **PASS** (all ingredients at the region cell, worst ratio 2.1e-4; `ab rho <= 0.0388`) |
| 6. factor 3.01 | **PASS** (`3.0000135`) |
| 7. end-to-end `n(R-1) > n/(3n-j)` | **PASS** (11 cells, certified balls, margin `2(j-1)/(9n)` reproduced) |
| 8. corrected margin | **PASS**: `0.4775 - 0.1296 = 0.348` at `TAU_MIN = 1000` (script's `0.2917` is a valid coarser bound) |

**Claim-state recommendation:** promote the script's own claim -- `(3 + 3a^2 b - ab) N Phi > 1` on `{tau >= 1000, tau^2/N <= 1}` and its integer consequence `p_j^2/(p_{j-1}p_{j+1}) > 1 + 1/(3n-j)` for `j >= 1001`, `N >= (j+1)^2` -- to **independently-validated**. The certificate survives the corrected lemma constants with a margin of 0.348 against corrections of 0.130. Do **not** promote the assembled theorem (`sibuya_theorem.py`) until the two coverage blockers above are closed: the `(j+1)^2` strip and the `j >= 1001` vs `j >= 2001` mismatch introduced by the `TAU_MIN` change. Human approval required for any public claim.
