# Independent validation: `release/scripts/sparse_certificate_full.py`

**Time (measured with `date`):** Wed Sep 2 21:47:06 RDT 2026
**Validator:** independent-validator (did not write the certificate; did not import it)
**Object validated:** `projects/qg-bootstrap/release/scripts/sparse_certificate_full.py`, **untracked in git** at HEAD `ee4a760`; identified by sha256 `3726339a442c0cd7478543d6c85499c96cfa747411aa1aa85e1d97e971476ce7`
**Engine:** python-flint 0.9.0 (fmpq / fmpq_poly / fmpz / arb), Python 3.12.10. No sympy, no `fractions.Fraction`.
**Independence:** separate formulation. The script uses `fmpq_mpoly` in (a, b, E) and the recurrence `j E_j = SUM i L_i E_{j-i}`; the validator uses its own trivariate class (dict `(power of a, power of E) -> fmpq_poly in b`), the Taylor-sum exponential, its own `D` with `d_b E = a_2 E`, its own substitution `E -> 1 + a_2 b + b^2 Z`, its own sweep, and for the lemmas its own (non-linearised) derivation with constants computed by interval evaluation. Shared dependency: python-flint itself (unavoidable, documented). Predecessor report: `VALIDATION_SPARSE_CERTIFICATE_2026-09-02.md`; its Faulhaber / Newton / Cauchy-ball code was reused (validator's own code).
**Validation code:** `scratchpad/full/w1_phat.py` (items 1, 2, and the derived multipliers), `w4_Ci.py` (item 4 constants), `w3_cells.py` + `w3b_assembly.py` (item 3), `w5_cells.py` (item 5), `w1c_identity.py`. Load-bearing excerpts quoted below.
**Script under test, own output (subprocess, not imported):** `[a^0]=[a^1]=[a^2]=0`, `[a^3]/b^2(0) = 176/175`, (iii) `1.2469e-02`, (v) `9.256e-08`, (vi) `2.431e-06`, sweep min `0.141702` (3 boxes), VERDICT True, exit 0, 2 s.

## Claim under test

Real-variable certificate: with `a = 1/tau`, `b = tau^2/N`, `F = SUM_j e_j (tau)_j/(N)_j`, `Phi = -(log F)''`: **`Phi N > 4/5` for all real `tau >= 627` and all `N` with `b <= 1`**, no side constraint on `a`. Integer consequence: `M(n,t) > 4/5` for integer `t >= 628`, `t^2/N <= (628/629)^2`.

## Validation plan (declared before running)

| # | item | pass criterion (deterministic) |
|---|---|---|
| 1 | Phat in Q[a,b,E] with own representation | `[a^0]=[a^1]=[a^2]=0` exactly; E-degree of Phat = 1 (E^2 cancels); `[a^3]` divisible by `b^2`; `[a^3]/b^2(0) = 176/175` exactly |
| 2 | Z bracket and sweep | Leibniz bracket holds (exact Z vs bracket, arb, 5 values of b); own sweep min > 0 after subtracting (iii), (v), (vi) |
| 3 | tail lemma (v) | each ingredient holds with exact `E_j` at 3 cells x 5 tau; the assembly constants 1.06 / 0.66 derived from an exact identity, or replaced |
| 4 | moment lemma (vi) | own `C_i`, i <= 30, agree with the script's method; `C_i >= m^2 |dev_i|` at 3 values of m; `h` formula checked; domain of m checked |
| 5 | end-to-end | `M > 4/5` as a certified arb ball at 6 cells, `b in {0.5, 1}`, `t in {628, 2000, 20000}` |
| 6 | state-level | window `[t-1, t+1]` maps into the region by an exact inequality |

## Results

### Item 1 -- Phat in Q[a, b, E]: PASS

Own representation, same sizes as the script: `P~`: 795 terms, degrees (58, 44, 1); `Phat`: 5119 terms, degrees (118, 89, 1).

- `[a^0] P~ = E` exactly (the a = 0 part of the truncated F~ is the exact exponential, by construction).
- `[a^0] Phat = [a^1] Phat = [a^2] Phat = 0` identically in Q[b, E].
- **E^2 cancels identically.** Derivation (validator): with `DQ = 0`, `D(EQ) = 2 a_2 ab EQ`, `D^2(EQ) = 2 a_2 a^2 b (1 + 2 a_2 b) EQ`, so the E^2 coefficient of Phat is `a^2 b Q^2 [4 a_2^2 b - 2 a_2 (1 + 2 a_2 b) - 4/5] = a^2 b Q^2 (-2 a_2 - 4/5) = 0` for `a_2 = -2/5`. Confirmed exactly in the polynomial (`E^2 part of Phat is zero: True`).
- After `E -> 1 + a_2 b + b^2 Z`: `[a^3]` has 33 terms, lowest b-power 2, Z-degree 1; `[a^3]/b^2 |_{b=0} = 176/175` exactly; every `[a^k]`, k >= 4, is divisible by `b^2` (lowest power 2).

`[a^3]/b^2` as a polynomial in (b, Z), validator's own output:

```
176/175 - 704/875 b + 1056/4375 b^2 - 2816/65625 b^3 + 352/65625 b^4 - 1408/2734375 b^5 + 704/17578125 b^6
 - 11264/4306640625 b^7 + 352/2392578125 b^8 - 1408/193798828125 b^9 + 7744/24224853515625 b^10
 - 512/40374755859375 b^11 + 832/1816864013671875 b^12 - 256/16870880126953125 b^13
 - 10048/2165096282958984375 b^14 + 185728/97429332733154296875 b^15 + 24064/487146663665771484375 b^16
 + Z (176/175 b^2 - 352/875 b^3 + 352/4375 b^4 - 704/65625 b^5 + 352/328125 b^6 - 704/8203125 b^7
      + 704/123046875 b^8 - 1408/4306640625 b^9 + 352/21533203125 b^10 - 704/968994140625 b^11
      + 704/24224853515625 b^12 - 128/121124267578125 b^13 + 64/1816864013671875 b^14
      - 128/118096160888671875 b^15 - 6592/1299057769775390625 b^16 - 12032/97429332733154296875 b^17)
```

```python
# validator, item 1 (excerpt): own operator with the E rule, own substitution
def dbE(s):
    return TP(
        {(ka, ke): p.derivative() + p * (A2 * ke) for (ka, ke), p in s.d.items()}
    )  # d_b + a2 E d_E


def D(p):
    return (
        TP.mono(2, 0, 0, -1) * p.da() + TP.mono(1, 1, 0, 2) * p.dbE()
    )  # -a^2 d_a + 2ab(d_b + a2 E d_E)


assert all(D(TP.mono(2 * k, k)).is_zero() for k in range(1, 6))
assert (D(Evar) - TP.mono(1, 1, 1, 2 * A2)).is_zero()
```

**Structural observation (not required, but load-bearing for item 5).** With the exact `Z(b)`, `[a^3]/b^2 = (176/175) e^{2 a_2 b} + r(b)` with `r = -1.7e-23, -2.6e-19, -3.5e-15` at `b = 0.25, 0.5, 1` (`O(b^14)`, the J = 30 / T_15 truncation). Since `P~ -> E` as `a -> 0`, `(N Phi - 4/5)/(ab) -> [a^3]/(b^2 E^2) = 176/175 + O(b^14)` for **every** b: the leading correction is b-independent. This is what item 5 measures (1.0058 at t = 20000 at both b = 0.5 and b = 1).

### Item 2 -- Z bracket and sweep: PASS

`Z(b) = SUM_{k>=2} a_2^k b^{k-2}/k!` is alternating (`a_2 < 0`) with `|term_{k+1}/term_k| = |a_2| b/(k+1) <= 0.4/3 < 1` for `b <= 1`, so by Leibniz `Z in [a_2^2/2 + a_2^3 b/6, a_2^2/2]`. Checked with the exact `Z` in arb at b = 0.001, 0.1, 0.5, 0.9, 1: in bracket at all five (e.g. b = 1: `Z = 0.0703200460`, bracket `[0.0693333, 0.08]`). Hence `|Z| <= 0.08` as the script uses.

(iii) recomputed: `SUM_{k>=4} |[a^k]/b^2| a^{k-3} <= 1.246862e-02` at `a = 1/627, b <= 1, |Z| <= 0.08` (script `1.2469e-02`).

Own sweep, 256 boxes in b with the Z box `[Zlo(b_hi), 0.08]`:

| quantity | validator | script |
|---|---|---|
| `min_b [a^3]/b^2` | `0.44999633` (true min is at b = 1: `0.4518965582`) | -- |
| minus (iii) | `0.43752771` | -- |
| minus (iii), (v), (vi) x `P~^2` | **`0.43752442`** | `0.141702` (3 boxes) |

The script's `0.141702` is the lower end of a crude 2-box enclosure (the polynomial evaluated on `b in [0, 1/2]` and `[1/2, 1]` as balls), not the minimum; it is a valid lower bound. The true margin is three times larger. Values of `[a^3]/b^2` at b = 0, 0.25, 0.5, 0.75, 1: `1.0057143, 0.8234092, 0.6741504, 0.5519477, 0.4518966`, and `x e^{0.8 b}` = `1.0057143` at each.

### Item 3 -- tail lemma (v): PASS, with the assembly constants now derived (they were asserted)

**Ingredients, numerically with exact `E_j`** (exact `e_j(beta)` by integer DP over the full spectrum at (m, t) = (300, 20), (2000, 60); by Newton from Faulhaber power sums at the region cell (m, t) = (628^2/2, 628), b = 1, j <= 120); at tau = t-1, t-1/2, t, t+1/2, t+1:

| cell | `\|mu_i\| <= 2^i` | Cauchy `\|E_j - E_j^0\| <= (2e mu_2 b/j)^{j/2}(e^g - 1)` | worst ratio | `\|1-R_j\| <= a j(j-1)/2` | worst ratio ((j-1)a <= 1) | max `2ab rho_j` |
|---|---|---|---|---|---|---|
| (300, 20), b = 2/3 | True (i <= 21) | True, j = 2..21 | 0.0096 | True (all j incl. j = t+1, tau < t) | 0.970 | 0.296 |
| (2000, 60), b = 0.9 | True (i <= 61) | True, j = 2..61 | 0.0115 | True | 0.986 | 0.195 |
| (197192, 628), b = 1 | True (i <= 120) | True, j = 2..120 | 0.0126 | True | 0.998 | 0.0276 |

- `d/dtau E_j^0 = j a E_j^0`: exact (`E_j^0 = (a_2/N)^{j/2} tau^j/(j/2)!`); numeric relative residual `1.8e-26`.
- `2ab rho_j <= 0.1` on the region: `j <= t+1 <= tau + 2` gives `2ab rho_j = 2 sqrt((1+2a) a b/(2 mu_2)) <= 0.0632452` for `tau >= 627, b <= 1` (analytic; also the series `SUM_i L_i z^i` converges on `|z| = rho_j` since `2 theta rho_j < 1`).
- `|1 - R_j| <= a j(j-1)/2` is proved by the Weierstrass product inequality only for `(j-1) a <= 1`; the one excluded case is `j = t+1` with `tau in [t-1, t)`, where `|1 - R_j| <= 1 + K_j <= 2.7 <= a j(j-1)/2 (>= 314)`. Holds; the docstring should say so.
- `|E_j^0| <= (2 e mu_2 b/j)^{j/2}`: from `(j/2)! >= (j/(2e))^{j/2}`, in fact `|E_j^0| <= (e mu_2 b/j)^{j/2}`.
- The `e` in the base: the validator's Cauchy estimate gives `e^{j/4}` (from `exp(|L_2| rho^2)`, `|L_2| rho^2 = j/4`); the script's `e^{j/2}` is a safe overestimate.
- Terms `j > t+1` (`-E_j^0`, no `R_j`): the script bounds them with the same formula; valid because `a j(j-1)/2 >= (j-1)/2 > 1` there. Not stated in the docstring.

**The assembly.** The docstring asserts "effect on Phi N <= 1.06 N SUM (ja)^2 |T_j| / F" with `F >= 0.66`. Neither constant is derived. Validator's derivation from the exact identity `Phi(F~ + T) - Phi(F~) = -(log(1+u))''`, `u = T/F~`:

`|Delta Phi| <= [ |T''|/F + 2|T'||F'|/F^2 + |T||F''|/F^2 + 2|T||F'|^2/F^3 ]/(1-|u|) + |u'|^2/(1-|u|)^2`,

with, for each `j <= t+1`, `|(tau)_j^{(k)}| <= j^{(k)} tau^{j-k}` (all `|tau - i| <= tau` since `j - 1 <= 2 tau`), so `|d^k[(E_j - E_j^0) R_j]| <= (ja)^k |E_j - E_j^0| K_j` and `|d^k[E_j^0 (R_j - 1)]| <= (ja)^k |E_j^0| K_j a j(j-1)/2 (1 + ab)`, where `K_j = 1/prod_{i<j}(1 - i/N) <= exp(j(j-1)/(2(N_min(j) - j)))`, `N_min(j) = max(627^2, (j-2)^2)` (`K_31 <= 1.0012`; `K_j -> e^{0.5}` only where the terms are `< 1e-600`). The `F` constants by interval sweep (512 boxes in b, `a in [0, 1/627]`, `E = exp(a_2 b)` as a ball):

| constant | validator (certified) | script |
|---|---|---|
| `F~ >=` | **0.66920516** | 0.66 (justified; `e^{-0.41} = 0.6637` is below the true min) |
| `\|F~'\| <=` | 0.802424 `a b F~` | -- |
| `\|F~''\| <=` | 0.800635 `a^2 b F~` | -- |
| `P~ <=`, `P~^2 <=` | 1.0000029, 1.0000057 | multiplies by 176/175 = 1.005714 (conservative, unexplained) |

Hence the multiplier per term is `m_j = [1 + 2 c_1 b/j + c_2 b/j^2 + 2 c_1^2 b^2/j^2] K_j (1+ab)/(1-|u|)` (+ a quadratic `u'^2` term). At the worst tail term `j = 31, b = 1`: `m_31 = 1.05394 x 1.0012 x 1.0016 = 1.0569`. **So 1.06 is (barely) justified for the tail (v)**, by a derivation the script does not contain. Recomputed (v) with the derived multipliers, `K_j`, the `(1+ab)` factor and a geometric remainder for `j >= 600` (the script simply breaks at `eff < 1e-300` without bounding the remainder): **`tail_rel <= 9.092e-08`** at b = 1 (script `9.256e-08`); at b = 0.25, 0.5, 0.75: `4.4e-16, 6.1e-12, 1.7e-9` -- increasing in b, so evaluating at b = 1 is the sup. `|T/F| <= 1.69e-10`, so the interpolant `F = F~ + T > 0` on the window.

**Measured at the region cell** (exact `e_j`, tau = 627.5 and 628.5, `F` vs `F~_m` = same finite-m moments for `j <= 30` plus `E - T_15`, so that only the tail differs):

| tau | `T = F - F~_m` | measured `N(Phi(F) - Phi(F~_m))` | relative to ab | script assembly with the TRUE `T_j` | measured/assembled |
|---|---|---|---|---|---|
| 627.5 | -7.45e-21 | 1.0966e-17 | 6.89e-15 | 2.383e-17 | 0.460 |
| 628.5 | -7.82e-21 | 1.1502e-17 | 7.22e-15 | 2.500e-17 | 0.460 |

The exact identity reproduces the measured effect to `1e-8`; the assembled bound holds with factor 0.46; the Cauchy-bounded script number (9.256e-8 relative) is 7 orders above the measured effect. (A first measurement that compared `F` against the certificate's `F~` gave `-9.7e-9` relative -- that is the moment perturbation (vi), see item 4, not the tail.)

```python
# validator, item 3 (excerpt): the derived per-term multiplier and the K_j factor
Nmin = max(627**2, (j - 2) ** 2)
Kj = (arb(j) * (j - 1) / (2 * (Nmin - j))).exp()  # 1/prod(1-i/N) <= exp(j(j-1)/(2(N-j)))
mj = 1 + 2 * c1 * b / j + c2 * b / j**2 + 2 * c1 * c1 * b * b / j**2
tau_j = base * ((g.exp() - 1) + amax * j * (j - 1) / 2 * (1 + amax * b)) * Kj
S_eff += j * j * tau_j * mj / (b * F)  # N (ja)^2 = j^2/b
```

### Item 4 -- moment lemma (vi): PASS on the numbers, two rigor gaps named

**Constants.** Own `C_i = sup_{m >= m_min} m^2 |mu_i(m) - mu_i^inf|`, i = 2..30, by Bernoulli-Faulhaber rational functions in `u = 1/m` on the ball `[0, 1/m_min]` (the deviation is `O(1/m^2)`: the `u^0, u^1` coefficients vanish for every i, asserted). Exact deviations at m = 196565, 393129, 1e7 are all `<= C_i` (i = 2..30). Excerpt of the table (full table in `w4.log`):

| i | `C_i` (m >= 627^2/2) | `C_i` (m >= 627^2, script's domain) | `m^2 dev` at m = 1e7 |
|---|---|---|---|
| 2 | 0.600000000010 | 0.600000000002 | -0.6000000 |
| 3 | 2.05714285721 | 2.05714285716 | -2.0571429 |
| 10 | 1156.01550193 | 1156.01550163 | -1156.0155 |
| 20 | 2492481.52409 | 2492481.52150 | -2492481.5 |
| 30 | 3893772838.32 | 3893772829.27 | -3893772804 |

`|mu_i(m)| <= 2^i` at m = 100, 196565 and `|mu_i^inf| <= 2^i`, i <= 30: True.

**Gap A (domain of m).** The script computes `C_i` for `m >= TAU_MIN^2 = 627^2`, but the region `tau >= 627, b <= 1` only gives `N = 2m = tau^2/b >= 627^2`, i.e. `m >= 627^2/2 = 196564.5`. On `[627^2/2, 627^2)` the script's `C_2` is exceeded (`m^2 |dev_2| = 0.6(1 + 1/(4m^2-1))`, excess `2.4e-12` relative at m = 196565); for i >= 3 the script's constants happen to hold there (`-1.8e-11` to `-6.4e-9` relative slack). Numerically nothing; formally the stated domain is wrong by a factor 2 in m. Fix: `mu_deviation_bounds(J, TAU_MIN**2 // 2)`.

**Gap B (the i-sum in h for j > J).** For `j > J` the script bounds `|E_j(mu) - E_j(mu^inf)|` with `h(rho) = 4 a^4 b^3 rho^2 SUM_{i=2}^{30} (C_i/i)(ab rho)^{i-2}`, but `E_j` for `j > 30` depends on `mu_i` for all `i <= j`, and no `C_i` for `i > 30` exists in the script. Validator's replacement, needing only `C_2`: for `j > J` split `T_j = [E_j(m) - E_j^0(m)] R_j + [E_j^0(m) - E_j^0(inf)] R_j + E_j^0(inf)(R_j - 1)`, where the first bracket is the Cauchy estimate with the TRUE moments (only `|mu_i(m)| <= 2^i` and `mu_2(m) <= 4/5` are used, so it is the same `g`), and `|E_j^0(m) - E_j^0(inf)| <= (j/2)(C_2/m^2)(b/2)(0.4 b)^{j/2-1}/(j/2)!`. Recomputed this way: the `j > J` moment contribution is `< 1e-12` relative.

**Gap C (the multiplier for j <= J).** The script applies the same `1.06 N (ja)^2 |M_j|/F` assembly to the moment perturbation of the terms `j = 2..30`. For `j = 2` the derived multiplier is `m_2 = 1 + 0.8 b + 0.2 b + 0.32 b^2 = 2.32` at b = 1 (the `2|T'||F'|/F^2` term is comparable to `|T''|/F` when `j ~ b`), so **1.06 is not a valid constant for (vi)**. With the derived multipliers (and gaps A, B closed): **`mom_rel <= 3.201e-06`** at b = 1 (script `2.431e-06`); b = 0.25, 0.5, 0.75: `1.05e-7, 4.7e-7, 1.36e-6` (increasing in b). Measured at the region cell (`F~_m` vs `F~`): `-9.67e-9` relative to ab at tau = 627.5, `-9.65e-9` at 628.5 -- 330x below the corrected bound.

**h formula.** `Delta S = SUM_i b (-1)^{i-1} (delta_i/i)(ab)^{i-2} z^i`, `|Delta S| <= SUM_i b (C_i/i) 4a^4 b^2 (ab)^{i-2} rho^i = h(rho)` on `|z| = rho`; `|exp(S_inf)| <= e^{j/4} e^{g}`; so `|E_j(mu) - E_j(mu^inf)| <= rho^{-j} e^{j/4} e^g (e^h - 1) <= (2e mu_2 b/j)^{j/2} e^g (e^h - 1)`. Correct (with the same `e^{j/4} -> e^{j/2}` slack). The script's `x**(i-2)/2**(i-2)` is `(ab rho)^{i-2}`: correct.

### Item 5 -- end-to-end, integer claim at real cells: PASS

Exact `mu_i(m)` (Faulhaber), exact `e_j` (Newton) for `j <= 170`, `F(s) = SUM_{j<=170} e_j (s)_j/(N)_j` exact in fmpq plus a Cauchy tail ball for `170 < j <= s` (needs `s <= N mu_2/8`: 20001 <= 4e7), `g = 2 log F(t) - log F(t-1) - log F(t+1)`, `M = n(e^g - 1)`, 2000 bits:

| t | N | b | `M - 4/5` | `(M - 4/5)/(ab)` | tail radius |
|---|---|---|---|---|---|
| 628 | 628^2 = 394384 | 1.0 | +1.60714689653e-3 | 1.009288251 | 5.5e-137 |
| 628 | 2 x 628^2 | 0.5 | +8.02932597616e-4 | 1.008483343 | 9.4e-163 |
| 2000 | 4e6 | 1.0 | +5.03417491904e-4 | 1.006834984 | 5.1e-137 |
| 2000 | 8e6 | 0.5 | +2.51645785904e-4 | 1.006583144 | 9.1e-163 |
| 20000 | 4e8 | 1.0 | +5.02913146346e-5 | 1.005826293 | 5.1e-137 |
| 20000 | 8e8 | 0.5 | +2.51450286430e-5 | 1.005801146 | 9.0e-163 |

`M > 4/5` at all six cells (the displayed `+/-` in the log is display rounding; the certified radii are the tail radii listed). `(M - 4/5)/(ab) = 176/175 + ~2.2 a` (1.0093 at a = 1/628, 1.0068 at 1/2000, 1.0058 at 1/20000), at both b = 0.5 and b = 1 -- exactly the b-independent leading term found in item 1. The task's expectation "close to 1.006 for t large" is confirmed; note the certificate's guaranteed value at b = 1 is only `0.4375/P~^2 >= 0.4375` because it bounds `[a^3]/b^2 = 0.452` and `P~^2` separately, whereas the true ratio is `[a^3]/(b^2 E^2) = 176/175`.

```python
# validator, item 5 (excerpt): exact main sum, Cauchy tail ball, M
main = sum((e[j] * ff_exact(fmpq(s), j) / ff_exact(fmpq(N), j) for j in range(jmax + 1)), fmpq(0))
for j in range(jmax + 1, jtop + 1):
    tail += (2 * e1 * N * qa(mu2) / j) ** (arb(j) / 2) * (sA / (N - j)) ** j
g = 2 * Fs[t].log() - Fs[t - 1].log() - Fs[t + 1].log()
M = n * (g.exp() - 1)
```

### Item 6 -- state-level, window to region: PASS

For integer `t >= 628`, `t^2/N <= (628/629)^2`, and `tau in [t-1, t+1]`: `b = tau^2/N <= ((t+1)/t)^2 (t^2/N) <= (629/628)^2 (628/629)^2 = 1`, and `a = 1/tau <= 1/(t-1) <= 1/627`. So the whole window lies in the certified region; no second constraint is needed (the predecessor's `b^3 <= a/125` is gone because the exponential is now exact). The integral identity `g = INT_0^1 INT_0^1 Phi^{(t)}(t-1+u+v) du dv` and the last step `M = n(e^g - 1) > 4/5` are unchanged from the predecessor report (6b, 6c) and rely only on `Phi^{(t)} N > 4/5` on the window and on `F^{(t)} > 0` there (`F~ >= 0.6692`, `|T| <= 1.7e-10`). The Cauchy domain `t+1 <= N mu_2/8` holds since `N >= t^2 >= 628^2`. The docstring's "Consequence" line is correct as written.

### Minor observations (do not affect validity)

1. `ok = float(worst) > 0` and `if j > 200 and float(eff) < 1e-300` compare floats of arb midpoints (project law: exact comparisons); harmless here (margin 0.14 / 0.44).
2. The remainder of the j-sums beyond the break point is not bounded; the validator's geometric remainder is `< 1e-600`.
3. The sweep's 3 boxes give a lower bound 3x below the true minimum; 256 boxes give `0.4375`.
4. The factor `176/175` applied to `(v)+(vi)` stands in for `P~^2 <= 1.0000057`; it is conservative but the docstring does not say what it is.

## Failure map

| boundary | status |
|---|---|
| `tau -> 627` (a max) | certified; margin flat in a; (iii) = 0.0125 is the only a-dependent correction |
| `b -> 1` (dense end) | certified; `[a^3]/b^2 = 0.452`, corrections `1e-7`; true `(N Phi - 4/5)/(ab) = 1.0093` at (628, 628^2) |
| `b -> 0` | certified; `[a^3]/b^2 -> 176/175`, all lemma bounds vanish as powers of b |
| `b > 1` | not certified; Leibniz bracket for Z needs `\|a_2\| b/3 < 1` (fine to b = 7.5) but (iii) and the tail base `(2e mu_2 b/31)^{15.5}` grow; not tested |
| `j = t+1`, `tau < t` (`(j-1)a > 1`) | `\|1 - R_j\| <= 2.7 <= a j(j-1)/2`; holds, undocumented |
| `m in [627^2/2, 627^2)` | script's `C_i` domain excludes it; `C_2` exceeded by 2.4e-12 relative there (gap A) |
| `j > 30` moment terms | script's `h` omits `delta_i`, i > 30 (gap B); replaced by a `C_2`-only split, `< 1e-12` |
| `j = 2` moment term | script multiplier 1.06 vs derived 2.32 (gap C); corrected (vi) = 3.2e-6 |

## Validation record (schema: templates/validation.yaml)

```yaml
id: VAL-SPARSE-CERT-FULL-2026-09-02
experiment_id: sparse_certificate_full.py@sha256:3726339a (untracked at HEAD ee4a760)
validator: independent-validator
independence_level: separate-formulation   # shared dependency: python-flint 0.9.0 only
checks:
  known_answer: pass          # [a^0..2] = 0, E^2 cancellation (-2 a_2 - 4/5 = 0), 176/175 exact, C_i to 10 digits
  negative_control: pass      # Z bracket must contain exact Z and does; deviation must be O(1/m^2) and is; E^2 must vanish and does
  hidden_points: pass         # 6 integer cells (b in {0.5,1}, t in {628,2000,20000}) and 3 ingredient cells x 5 tau, not chosen by the implementer
  precision_convergence: pass # 300-2000 bits; exact fmpq where possible; certified balls with explicit radii (tail 5e-137)
  boundary_map: pass          # table above
  signature: pass             # (M-4/5)/(ab) -> 176/175 at every b, as the polynomial predicts
decision: pass
allowed_claim_state: independently-validated
blockers:
  - "rigor (non-blocking numerically): the (vi) assembly constant 1.06 is invalid for j = 2 (derived multiplier 2.32 at b = 1); corrected bound 3.2e-6 vs margin 0.4375"
  - "rigor (non-blocking numerically): C_i domain must be m >= 627^2/2, not 627^2; and for j > 30 the h-sum must not stop at i = 30 -- use the C_2-only split given in item 4"
  - "wording: the 1.06 / 0.66 constants of (v) need the derivation given in item 3 (F~ >= 0.6692, |F~'| <= 0.8025 ab F~, |F~''| <= 0.8007 a^2 b F~, K_j, (1+ab)); the j = t+1, tau < t case of |1 - R_j| and the j > t+1 terms need one sentence each"
  - "housekeeping: the script is untracked in git; commit before any claim promotion"
```

## Verdict

| item | verdict |
|---|---|
| 1. Phat in Q[a,b,E]: `[a^0]=[a^1]=[a^2]=0`, E^2 cancels, `[a^3]/b^2(0) = 176/175` | **PASS** (exact) |
| 2. Z bracket; sweep of `[a^3]/b^2 - (iii)` | **PASS**; own min `0.4375` (script's `0.1417` is a valid coarse enclosure) |
| 3. tail lemma (v): ingredients, `2ab rho_j <= 0.063`, `F >= 0.6692`, assembly | **PASS**; 1.06 justified only after the validator's derivation (`m_31 = 1.0569`); own (v) = 9.09e-8 |
| 4. moment lemma (vi): `C_i`, `h`, assembly | **PASS on the numbers**; three rigor gaps (domain of m, i > 30 for j > 30, multiplier at j = 2); corrected (vi) = 3.2e-6 |
| 5. end-to-end `M > 4/5`, b in [0.5, 1], t up to 20000 | **PASS** (6 cells, certified balls, `(M-4/5)/(ab) -> 176/175`) |
| 6. window in region, integer consequence | **PASS** (exact inequality) |

**Claim-state recommendation:** promote `Phi N > 4/5` on `{tau >= 627, tau^2/N <= 1}` and its integer consequence `M(n,t) > 4/5` for `t >= 628`, `t^2/N <= (628/629)^2` to **independently-validated**. The certificate's conclusion survives every attack: the corrected lemma constants change the subtracted corrections from `2.5e-6` to `3.3e-6` against a margin of `0.4375`. The three rigor blockers are text/constant fixes to the lemmas, not mathematical gaps in the conclusion; they should be applied before release-ready. Human approval required for any public claim.
