# Independent validation: `release/scripts/sparse_certificate.py`

**Time (measured with `date`):** Wed Sep 2 21:06:19 RDT 2026
**Validator:** independent-validator (did not write the certificate; did not import it)
**Object validated:** `projects/qg-bootstrap/release/scripts/sparse_certificate.py` at commit `226af459657829149a73bc640345ea11befacf73`
**Engine:** python-flint 0.9.0 (fmpq / fmpq_poly / fmpz / arb), Python 3.12.10. No sympy, no `fractions.Fraction`.
**Independence:** separate formulation. The certificate uses `fmpq_mpoly` and the recurrence `j E_j = SUM i L_i E_{j-i}`; the validation uses its own bivariate class (list of `fmpq_poly` in b indexed by the power of a), the Taylor-sum exponential, Bernoulli-number Faulhaber sums, and Newton's identities. Shared dependency: python-flint itself (unavoidable, documented).
**Validation code:** six scripts, run in the session scratchpad (`v1_identity.py`, `v23_phat.py`, `v4_tail.py`, `v5_mu.py` + `v5b_Ci.py` + `v5c_compare.py`, `v6_cell.py`, `v6_tailbound.py`, `v6_mu_spot.py`); load-bearing excerpts are quoted below.
**Script under test, own output (run as a subprocess, not imported):** `total_low = 0.997424`, `tail_rel <= 6.537e-04`, `mu_rel <= 9.974e-11`, VERDICT True, exit 0, 1.9 s.

## Claim under test

Integer statement: for the doubled centred-square spectrum `b_k = (2k-1)^2` (each twice), `N = 2m`, `p_t = e_t/C(N,t)`, `n = N+1`,
`M(n,t) = n(p_t^2/(p_{t-1}p_{t+1}) - 1) > 4/5` for every integer `t >= 628` and every `N` with `t^2/N <= 1/100` and `(t^2/N)^3 <= 1/(125 t)`.
Real-variable certificate: `Phi(tau;N) N > 4/5` for `tau >= 627` on `{a <= 1/627, b <= 1/100, b^3 <= a/125}`, `a = 1/tau`, `b = tau^2/N`.

## Validation plan (declared before running)

| # | item | pass criterion (deterministic) |
|---|---|---|
| 1 | sampling identity, product form of `(t)_j/(N)_j`, `E_j` = coefficient of `exp(...)` | exact `fmpq` equality at 4 cells, `m <= 400`, `t <= 40` |
| 2 | operator `D`, annihilation of `a^2 b`, sign/value of `Phi N - 4/5` vs `Phat/(a^2 b P^2)` | symbolic zero; relative difference of finite-difference `Phi` (h = 1e-6, 400 bits) vs `Phat` below 1e-9 at 5 points |
| 3 | `[a^k] Phat`, k = 0..3 | `[a^0] = [a^1] = 0`, `[a^2]` starts at `b^6`, `[a^3]` coefficient of `b^2` is exactly `176/175` |
| 4 | Cauchy tail bound on `e_j(beta)`, and `|beta_k| < 2` | ratio `|e_j|/bound < 1` for j = 2..30 at m = 200, 400 |
| 5 | `|mu_i(m) - mu_i^inf| <= C_i/m^2` | exact deviations vs independently computed sup constants; agreement with the printed 2-s.f. `C_i` |
| 6 | gap between integer claim and real certificate | each sub-question answered by an exact identity or a certified bound; the integer claim checked end-to-end at real cells |

## Results

### Item 1 -- exact sampling identity: PASS

`p_t = bbar^t SUM_{j<=t} e_j(beta) (t)_j/(N)_j` holds with exact equality (`fmpq`) at (m, t) = (50, 10), (100, 25), (200, 40), (400, 40), with `e_t(b)` computed directly by integer DP over all N roots and `e_j(beta) = e_j(3 b_k - (4m^2-1)) / (4m^2-1)^j`. Also exact: `(t)_j/(N)_j = theta^j prod_{i<j} (1-i/t)/(1-i/N)` for all `j <= t`, and `e_j theta^j = [z^j] exp(N SUM_{i>=2} (-1)^{i-1} mu_i theta^i z^i / i)` for all `j <= t`, where the exponential was computed as the truncated Taylor sum `SUM_k S^k/k!` (not the script's recurrence). Values: `p_10 = 1.159387e35` (m=50), `p_25 = 3.727277e102` (m=100), `p_40 = 2.352129e188` (m=200), `p_40 = 6.547119e212` (m=400).

```python
# validator, item 1 (excerpt)
S = fmpq_poly([0, 0] + [N * fmpq((-1) ** (i - 1)) * mu[i] * theta**i / i for i in range(2, t + 1)])
expS = fmpq_poly([1])
term = fmpq_poly([1])
for k in range(1, t // 2 + 1):
    term = fmpq_poly(list(((term * S) / k).coeffs())[: t + 1])
    expS += term
ok3 = all(expS[j] == e_beta[j] * theta**j for j in range(t + 1))
```

### Item 2 -- operator and sign equivalence: PASS

- Chain rule at fixed N: `da/dtau = -a^2`, `db/dtau = 2 tau/N = 2ab`; numeric residuals at (700, 1e8): `-4.2e-28`, `1e-114`. Hence `d/dtau = -a^2 d/da + 2ab d/db`.
- `D (a^2 b)^k = 0` exactly for k = 1..5 in the validator's own polynomial class (`D(a^2 b) = -2a^3 b + 2a^3 b`).
- `N Phi_J - 4/5 = Phat/(a^2 b P^2)` exactly (since `Q` is a function of `a^2 b`, `DQ = 0`, and `F = P/Q`). Numeric comparison of `N * Phi_fd - 4/5` (central second difference of `log(P/Q)` in tau, h = 1e-6, 400 bits) with `Phat/(a^2 b P^2)`:

| tau | N | b | `N Phi_fd - 4/5` | `Phat/(a^2 b P^2)` | rel. diff | `(N Phi - 4/5)/(ab)` |
|---|---|---|---|---|---|---|
| 627 | 3.931e7 | 1.00e-2 | 1.60432632215e-5 | 1.60432632215e-5 | 6.9e-24 | 1.0059126 |
| 1000 | 1e8 | 1.00e-2 | 1.00583863580e-5 | 1.00583863580e-5 | 1.7e-24 | 1.0058386 |
| 2000 | 1e9 | 4.00e-3 | 2.01154688073e-6 | 2.01154688073e-6 | 8.4e-26 | 1.0057734 |
| 627 | 1e10 | 3.93e-5 | 6.30697182394e-8 | 6.30697182394e-8 | 2.7e-26 | 1.0058966 |
| 5000 | 2.5e9 | 1.00e-2 | 2.01147838354e-6 | 2.01147838354e-6 | 2.5e-26 | 1.0057392 |

Same sign and same value to 24 digits at every point.

### Item 3 -- Taylor coefficients of Phat: PASS

Recomputed with the validator's own representation (P: a-degree 18, b-degree 9; Phat: a-degree 38, b-degree 19 -- same as the script reports):

- `[a^0] = 0`, `[a^1] = 0` identically.
- `[a^2](b) = 176/234375 b^6 - 96/390625 b^7 + 224/5859375 b^8 - 64/17578125 b^9 + 32/146484375 b^10 - 64/10986328125 b^11` -- lowest power **b^6** (the script asserts only `>= 5` and bounds with `d5 b^5`; both consistent).
- `[a^3](b) = 176/175 b^2 - 704/875 b^3 + 1408/4375 b^4 - 5632/65625 b^5 - 176/28125 b^6 + 136256/24609375 b^7 - 137824/123046875 b^8 + 242432/1845703125 b^9 - 88544/9228515625 b^10 + 13952/46142578125 b^11`; coefficient of `b^2` is exactly **176/175**.

Every coefficient agrees with the script's printed polynomials.

Region bounds (i)-(iv) recomputed by the validator (own code, same interval-free absolute-coefficient method):

| region | `[a^3]/b^2 >=` | `|[a^2]| <=` | rest `<=` | `Phat/(a^3 b^2) >=` |
|---|---|---|---|---|
| script: `a<=1/627, b<=1/100, b^3<=a/125` | 0.997636 | 7.5339e-6 b^5 | 2.1229e-4 b^2 | **0.997424** (matches script) |
| integer-claim: `b <= (629/628)^2/100, b^3 <= (629/628)^7 a/125` | 0.997610 | 7.5580e-6 b^5 | 2.1238e-4 b^2 | **0.997398** |

### Item 4 -- Cauchy tail bound: PASS

`|beta_k| < 2`: max `|beta_k|` = 1.985037 (m = 200), 1.992509 (m = 400); analytically `max beta = 3(2m-1)/(2m+1) - 1 < 2`, `min beta = 3/(4m^2-1) - 1 > -1`. `N mu_2/8` = 40.0 (m = 200), 80.0 (m = 400), so `j <= 30` is inside the bound's domain. Exact `e_j(beta)` vs `(2 e N mu_2/j)^{j/2}` (with the exact finite-m `mu_2`):

| m | j=2 | j=3 | j=5 | j=10 | j=20 | j=30 | worst ratio |
|---|---|---|---|---|---|---|---|
| 200 | 1.84e-1 | 4.36e-3 | 4.26e-3 | 4.78e-3 | 4.95e-5 | 1.89e-7 | 0.184 |
| 400 | 1.84e-1 | 3.09e-3 | 3.03e-3 | 5.12e-3 | 7.99e-5 | 6.90e-7 | 0.184 |

The bound holds with a factor >= 5.4 to spare at every j. Derivation checked: `|(1+w)e^{-w}| <= e^{|w|^2}` for `|w| <= 1/2` (`SUM_{k>=2} |w|^k/k <= |w|^2 (1/2 + 1/6 + 1/16 + ...) < |w|^2`), radius `r^2 = j/(2N mu_2)`, `r <= 1/4` iff `j <= N mu_2/8`.

### Item 5 -- finite-m moments: PASS (with a precision note)

Independent constants: `mu_i(m)` as exact rational functions of m (Faulhaber via `fmpq.bernoulli`, checked against direct summation at m = 100), written in `u = 1/m`; the two lowest powers of u vanish for every i (deviation is `O(1/m^2)`); `C_i := sup_{m >= m_min} m^2 |mu_i(m) - mu_i^inf|` by arb evaluation on the ball `[0, 1/m_min]`, `m_min = 627^2 * 50`:

| i | validator `C_i` | exact `m^2 (mu_i - mu_i^inf)` at m = 1000 | at m = 10000 | script print |
|---|---|---|---|---|
| 2 | 0.60000000 | -0.60000 | -0.600000 | 0.6 |
| 3 | 2.0571429 | -2.05714 | -2.057143 | 2.1 |
| 4 | 6.1714286 | -6.17142 | -6.171428 | 6.2 |
| 5 | 16.207792 | -16.20774 | -16.207792 | 16 |
| 6 | 40.423576 | -40.42335 | -40.423574 | 40 |
| 7 | 96.671329 | -96.67053 | -96.671321 | 97 |
| 8 | 225.25051 | -225.24791 | -225.250488 | 2.3e+02 |
| 9 | 514.18372 | -514.17581 | -514.183644 | 5.1e+02 |
| 10 | 1156.0155 | -1155.99268 | -1156.015273 | 1.2e+03 |

All nine printed constants are the 2-s.f. roundings of the validator's sup constants. `|m^2 dev_i| <= C_i` holds at m = 1000, 10000, m_min, 10 m_min for i = 3..10. For i = 2 the exact identity `mu_2(m) = 4/5 - 12/(5(4m^2-1))` (verified as a polynomial identity) gives `m^2 |dev_2| = 0.6 (1 + 1/(4m^2-1))`, which exceeds 0.6 by 1.5e-7 at m = 1000 and 1.5e-9 at m = 10000 -- both **outside** the certificate's domain `m >= m_min = 19,656,450`, where the excess is 3.9e-16 and is covered by the arb bound. So the literal test "at m = 1000, 10000 against the printed C_i" is not meaningful at the 1e-9 level; within the certificate's stated domain the bound holds. Not a defect.

Spot check of step (vi) (measurement, not a bound): `Phat/(a^3 b^2)` at the corner `(1/627, 1/100)` with the exact `mu_i(m)` instead of `mu_i^inf` changes by `-9.68e-11` (m = m_min) and `-9.35e-11` (m = 2e7), relative; the script's bound is `9.97e-11`. The bound is valid and tight to 3 percent at the corner (as expected, the perturbation is linear at leading order and the corner attains `1/m^2 = 4 a^4 b^2`).

### Item 6 -- attack on the integer-to-real step

**6a. Is the interpolant `F^(t)` positive on `[t-1, t+1]`?** PASS. `F^(t) = F_J + T`, and on the integer-claim region `F_J >= 1 - |P - 1| >= 0.995972` (absolute-coefficient bound `|P-1| <= 0.00403`, `Q <= 1`) while `|T| <= 6.49e-14` (non-linearised tail bound below). Measured at the five real cells: `min F^(t)` on 17 points of `[t-1, t+1]` is 0.996002 (t = 628, N = 100*628^2), 0.999984 (t = 628, N = 1e10). `F^(t)` is a polynomial, hence `C^2`.

**6b. Is the integral identity applied correctly?** PASS. `INT_0^1 INT_0^1 f''(x+u+v) du dv = f(x+2) - 2f(x+1) + f(x)`; with `x = t-1`, `f = log(bbar^tau F^(t)(tau))`: `Delta^2 f(t-1) = log p_{t+1} - 2 log p_t + log p_{t-1} = -g`, and `f'' = (log F^(t))'' = -Phi^{(t)}` (the `tau log bbar` term is linear). Hence `g = INT INT Phi^{(t)}(t-1+u+v) du dv` exactly, where `Phi^{(t)}` is built from the interpolant (sum to `j = t+1`), not from the full `F` (sum to N) -- this is the correct choice, because the Cauchy bound only covers `j <= N mu_2/8`, and `t+1 <= sqrt(N)/10 + 1 << N/10` on the region.

**6c. Does `M > 4/5` follow from `Phi N > 4/5`?** PASS. `g = INT INT Phi^{(t)} > 4/(5N)`; `4/(5N) > 4/(5n)` since `n = N+1`; `log(1+x) <= x` gives `4/(5n) >= log(1 + 4/(5n))`; so `e^g > 1 + 4/(5n)` and `M = n(e^g - 1) > 4/5`. This requires `Phi^{(t)} N > 4/5` at every real `tau in [t-1, t+1]`, i.e. the real certificate must hold on the image of the window -- see 6d.

**6d. Gap in the stated region (statement-level, closed here by recomputation).** The task's claim statement ("`t^2/N <= 1/100` and `(t^2/N)^3 <= 1/(125 t)`") does not literally follow from the certificate on `{b <= 1/100, b^3 <= a/125}`: for `tau = t+1` one has `b = (t+1)^2/N <= (629/628)^2/100 > 1/100` and `b^3/a = tau^7/N^3 <= (629/628)^7/125 > 1/125`. The docstring shrinks the b-constraint (`t^2/N <= BSTAR (628/629)^2`) but says nothing about the `b^3 <= a/125` constraint in its "Consequence" line. The validator therefore recomputed the whole real certificate on the enlarged region `{a <= 1/627, b <= (629/628)^2/100, b^3 <= (629/628)^7 a/125}`, which is exactly the image of the window for the task's integer statement: `Phat/(a^3 b^2) >= 0.997398 > 0`, tail effect `<= 5.76e-4` (below), mu effect `<= 9.97e-11` (script; scale unchanged). So the integer claim **as stated in the task** holds on its stated region; the docstring's "Consequence" sentence should add `(t^2/N)^3 <= (628/629)^7/(125 t)` or equivalently state the enlarged real region. In (t, N) terms the certified region is `N >= max(100 t^2, 5 t^{7/3})` up to these factors; the second constraint binds only for `t > 8000`.

**6e. Hidden assumptions in the tail step (v).** The script's comment derives `|T''| <= SUM (j/(tau-j+1))^2 tail_j <= 1.03 SUM j^2 a^2 N tail_j`; the per-term factor `(tau/(tau-j+1))^2` at `tau = 627` is 1.033 at j = 11 and 1.036 at j = 12, so the per-term "1.03" is false for `j >= 12` (immaterial: the j = 12 term is 40x smaller than j = 11), and the step `Phi(F+T) - Phi(F)` is linearised with an unproved fudge 1.03/0.98. The validator replaced this with an exact, non-linearised bound: `Phi(F+T) - Phi(F) = -(log(1+u))''`, `u = T/F`, `|(log(1+u))''| <= |u''|/(1-u_0) + u_1^2/(1-u_0)^2`, with `|(tau)_j^{(k)}| <= j^{(k)} tau^{j-k}` (each factor `|tau - i| <= tau`), `|F_J| >= 0.995972`, `|F_J'| <= 0.8039 ab`, `|F_J''| <= 0.8097 a^2 b` (absolute-coefficient sums of `DP`, `D^2P` over `Q_min`), `mu_2 <= 0.8 + 1e-12`, `(1-j/N)^{-j} <= exp(j^2/(N-j)) <= 1.01008` uniformly for `j <= t+1 <= sqrt(N)/10 + 1`, sum to j = 400 plus a geometric remainder. Result on the integer-claim region: **tail effect on `N Phi` relative to `ab` `<= 5.759e-4`** (script: 6.537e-4). Since `N Phi_J - 4/5 = Phat/(a^2 b P^2) >= 0.997398 ab / P^2` and `P^2 <= 1.00807`, the margin after the tail is `0.98941 - 0.00058 = 0.98883 > 0`. PASS. Measured true tail effect (J = 40 vs J = 10 series, exact `e_j(m)`, at `tau = t - 1/2`): `-4.5e-9` (t = 628, b = 0.01), `-3.7e-8` (t = 5000, b = 0.01), `-1.1e-18` (t = 628, b = 3.9e-5) -- five orders below the bound; it scales as `b^4/a` (the `exp(-2b/5)` truncation defect of `E_12 ∋ L_2^6/6!`), which is why some coupling constraint between a and b is intrinsic at fixed J; the specific `b^3 <= a/125` is an artefact of J = 10 and of the Cauchy bound, not of the theorem.

**6f. Hidden assumption on J.** `Q = prod_{i<10}(1 - i a^2 b) = prod (1 - i/N) > 0`, so `sign(N Phi_J - 4/5) = sign(Phat)`; `E_j` for `j <= 10` involves only `mu_2..mu_10`; the tail bound `bmax^{j/2-5}` needs `j >= 10` and is used for `j >= 11`. No hidden dependence found.

**6g. End-to-end check of the INTEGER claim at real cells** (exact `e_j(beta)` for `j <= 40` from Faulhaber power sums and Newton's identities at `m ~ 2e7`, Cauchy ball for `40 < j <= t+1`, 500 bits):

| t | N | b | g | `M - 4/5` | `(M - 4/5)/(ab)` | tail radius |
|---|---|---|---|---|---|---|
| 628 | 100*628^2 | 1.00e-2 | 2.02852047170150e-8 | +1.604611097e-5 | 1.0076958 | 1.1e-61 |
| 628 | 4e7 | 9.86e-3 | 2.00003948204831e-8 | +1.582082003e-5 | 1.0076955 | 8.1e-62 |
| 1000 | 1e8 | 1.00e-2 | 8.00010058386285e-9 | +1.006958647e-5 | 1.0069586 | 1.1e-61 |
| 5000 | 2.5e9 | 1.00e-2 | 3.20000804591323e-10 | +2.011926310e-6 | 1.0059632 | 1.1e-61 |
| 628 | 1e10 | 3.94e-5 | 8.00000063170290e-11 | +6.328228969e-8 | 1.0076798 | 5.5e-111 |

`M > 4/5` at all five cells, with `(M - 4/5)/(ab) = 1.006..1.008`, consistent with the certified `(N Phi - 4/5)/(ab) >= 0.989` plus the `O(1/N)` terms of `n(e^g - 1) - 4/5`. PASS.

```python
# validator, item 6g (excerpt): power sums by Bernoulli numbers, Newton identities, tail ball
def S_pow(p, M):  # SUM_{n=1}^{M} n^p
    if p == 0:
        return fmpq(M)
    M1 = fmpq(M + 1)
    return sum(
        fmpq(comb(p + 1, j)) * fmpq.bernoulli(j) * M1 ** (p + 1 - j) for j in range(p + 1)
    ) / (p + 1)


def esym_from_moments(mu, N, jmax):  # j e_j = SUM_{i=1}^{j} (-1)^{i-1} e_{j-i} N mu_i
    e = [fmpq(1)]
    for j in range(1, jmax + 1):
        e.append(sum(fmpq((-1) ** (i - 1)) * e[j - i] * N * mu[i] for i in range(1, j + 1)) / j)
    return e


g = 2 * Fs[t].log() - Fs[t - 1].log() - Fs[t + 1].log()
M = n * (g.exp() - 1)
```

### Minor observations (do not affect validity)

1. `float(total_low) > float(tail_rel) + float(mu_rel)` compares arb midpoints, not upper bounds (radii are ~1e-10 so harmless; formally against the project's no-float-comparison law).
2. `arb(0, 1/mmin)` uses a float radius; rounding is 1e-24 relative and is swallowed by the arb evaluation width (2e-13). Harmless.
3. The Cauchy bound in (v) uses `mu_2^inf = 4/5` where the true `mu_2(m) = 4/5 - 12/(5(4m^2-1)) < 4/5`, so it is on the safe side.
4. The `j > 200` remainder comment ("below 1e-300") is true but the `K_j <= 1.0021` statement is proved only for `j <= 200`; the validator's uniform `K_j <= 1.01008` for all `j <= t+1` closes this.

## Failure map

| boundary | status |
|---|---|
| `tau -> 627` (a max) | certified; margin 0.997 in `Phat/(a^3 b^2)` is essentially flat in a |
| `b -> 1/100` (dense end of the sparse regime) | certified; tail bound 5.8e-4 relative is the largest correction; true tail effect 4.5e-9 |
| `b -> 0` (N -> infinity at fixed t) | certified; all corrections vanish as powers of b |
| `b^3 -> a/125` (large t at fixed b) | certified; beyond it the Cauchy tail bound (not the theorem) fails -- true tail effect grows like `b^4/a` for J = 10 |
| `t+1 <= N mu_2/8` (Cauchy domain) | satisfied by a factor >= 6e3 on the region |
| window `tau in [t-1, t+1]` | requires the enlarged region of 6d; verified there |

## Validation record (schema: templates/validation.yaml)

```yaml
id: VAL-SPARSE-CERT-2026-09-02
experiment_id: sparse_certificate.py@226af459
validator: independent-validator
independence_level: separate-formulation   # shared dependency: python-flint 0.9.0 only
checks:
  known_answer: pass          # exact identities (item 1), exact Taylor coefficients (item 3), closed form mu_2(m)
  negative_control: pass      # [a^0]=[a^1]=0 must vanish and do; deviation must be O(1/m^2) and is
  hidden_points: pass         # 5 (tau,N) points for item 2, 5 real cells for item 6g, not chosen by the implementer
  precision_convergence: pass # 300-500 bits, finite-difference agreement 1e-24, arb balls with explicit radii
  boundary_map: pass          # table above
  signature: pass             # (M-4/5)/(ab) = 1.006..1.008 vs certified >= 0.989 and predicted 176/175 leading term
decision: pass
allowed_claim_state: independently-validated
blockers:
  - "statement-level: the docstring 'Consequence' line must add the constraint (t^2/N)^3 <= (628/629)^7/(125 t) (or state the enlarged real region); the validator verified the certificate on that enlarged region, so this is a wording fix, not a mathematical gap"
  - "rigor-level (non-blocking): tail step (v) as written is linearised with an unproved 1.03/0.98 fudge; an exact non-linearised bound (5.76e-4) is given in this report and could replace it"
```

## Verdict

| item | verdict |
|---|---|
| 1. exact sampling identity, product form, `E_j` exponential | **PASS** (exact equality, 4 cells) |
| 2. operator `D`, annihilation, sign equivalence | **PASS** (symbolic; numeric to 1e-24) |
| 3. `[a^0]=[a^1]=0`, `[a^2] = O(b^6)`, `[a^3] ∋ 176/175 b^2` | **PASS** (exact) |
| 4. Cauchy tail bound, `|beta_k| < 2` | **PASS** (worst ratio 0.184) |
| 5. `|mu_i(m) - mu_i^inf| <= C_i/m^2` | **PASS** on the certificate's domain `m >= 627^2*50`; the 2-s.f. prints are roundings of sup constants reproduced to 8 digits |
| 6. integer claim vs real certificate | **PASS** after enlarging the real region to the image of the window (6d); interpolant positive; integral identity exact; `M > 4/5` follows; end-to-end `M > 4/5` at 5 real cells |

**Claim-state recommendation:** promote the real-variable statement `Phi N > 4/5` on `{tau >= 627, tau^2/N <= (629/628)^2/100, (tau^2/N)^3 <= (629/628)^7/(125 tau)}` and its integer consequence `M(n,t) > 4/5` for `t >= 628`, `t^2/N <= 1/100`, `(t^2/N)^3 <= 1/(125 t)` to **independently-validated**, conditional on the one-line wording fix in the docstring (blocker 1). Human approval required for any public claim.
