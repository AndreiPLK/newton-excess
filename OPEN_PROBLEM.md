# The obstruction

**None left in the mathematics. The theorem is assembled** (`release/scripts/theorem.py`, 21:54,
2 September 2026, clock measured): `M(n,t) > 4/5` for every odd `n >= 5` and every `t < n/2`.

    A.  t <= 627, every n              the certificate ladder (627 rungs; being re-run into
                                       results/m_ladder_log_2026-09-02.txt as the durable artifact)
    B.  n = 1257, 1259, t >= 628       exact:  min M = 1.807, 1.804
    C1. t^2/N <= 1, every t >= 627     sparse_certificate_full.py   exit 0, validated (six PASS)
    C2. 0.05 <= theta <= 1/2           dense_certificate_a.py       exit 0, 930 boxes, validation running
    C3. theta <= 0.05, t^2/N >= 1      dense_certificate_b.py       exit 0, series, validation running
    coverage of C by C1 u C2 u C3      from khat_1 = theta/V in [0.30274, 0.33334] on V <= 0.17:
                                       V <= 0.1652, zeta <= 0.2223, eta <= 1.062e-3 -- inside the
                                       certified boxes (0.17, 0.24, 1.07e-3)

**What still has to happen before the word "proved" is used in public:** the two validation
reports for C2 and C3 (independent agents, no imports of the certificates), the ladder log
completing (it is a long run; the previous session's 627 rungs were seen only on a terminal),
and a written proof section in the paper that a referee can follow without running anything.

**The five constants of the day, each falling out of a computation aimed at something else:**
`4/5` as the exact constant term of `N g` in the sparse series; `176/175` as the exact `b^2`
coefficient of `[a^3]` and again as `S(0) = 176/525 = H'(0)/3` in the dense series; `E^2`
cancelling identically by `-2 a_2 - 4/5 = 0`; and `[a^3]/b^2 = (176/175) e^{2 a_2 b} + O(b^14)`
(found by the validator), so that `(N Phi - 4/5)/(a b) -> 176/175` for every `b`.

---

## Head of 21:30, kept

**One finite certificate left: the dense regime `t^2/N >= 0.95`, `theta <= 0.05`, uniformly as
`theta -> 0`, by a formal series in `v^2` with exact coefficients.** (21:30, 2 September 2026,
clock measured.)

What is certified and independently validated since 20:36:

- **Sparse regime, closed for every `t`:** `M(n,t) > 4/5` for all `t >= 628` with `t^2/N <= 1`
  (`release/scripts/sparse_certificate_full.py`, 2 s). The `a = 0` part of the sampling expansion
  is the exact exponential `E = e^{-2b/5}` carried as a formal variable with `d_b E = a_2 E`; then
  `[a^0] = [a^1] = [a^2] = 0` identically, `E^2` cancels identically (the identity
  `-2a_2 - 4/5 = 0`), `[a^3]/b^2 |_{b=0} = 176/175` exactly, and the sweep over `b in [0, 1]` has
  minimum `0.1417` after the tail (`9e-8`) and finite-`m` (`2.4e-6`) lemmas. The earlier version
  with the side constraint (`b <= 1/100`, `b^3 <= a/125`) was validated by an independent agent,
  six items PASS (`results/VALIDATION_SPARSE_CERTIFICATE_2026-09-02.md`); the full version is
  under validation now.
- **Dense regime (a), `0.05 <= theta <= 1/2`:** `release/scripts/dense_certificate_a.py` -- the
  fixed-tilt Fourier weight, the exact Edgeworth polynomial of `kappa_2[u]` to weight 6 with an
  explicit remainder (inner Taylor tail with `rhobar(u) <= (|u|/L) rhobar(L)`, cumulant tail by
  Cauchy, outer region by `|e^K| <= e^{-k_2(1-cos phi)}`), cumulants by polygamma with the
  Stirling remainder from the Bernoulli lemma (itself certified by a 1-D sweep), everything
  written through `eps = 1/m` so that `N -> infinity` is a box, two-pass mean-value enclosure in
  `v` (`d kt_j/dv = (2/v) kt_{j+1}` exactly), and the variance identity to remove the `1/theta`
  cancellation. Test boxes pass with `N g` in `[0.84, 1.9]` against the target `0.8`; the adaptive
  sweep over `(eps, v)` is running (about 2000 boxes certified per 8 minutes, no failure so far).

**Why the remaining piece needs a formal series and not boxes.** In the dense regime at fixed
`b = t theta >= 1` and `t -> infinity`, `theta -> 0` and the margin `theta/2` vanishes while the
pieces of `N g` stay `O(1)`; any box method needs `theta^{-1/k}` boxes per e-fold of `v` for a
`k`-th order enclosure, which diverges as `theta -> 0`. The cure is the same as for `H''` at the
edge: an exact power series in `V = v^2` whose coefficients are polynomials in `eta = eps/V`
(`t = 2 khat_1/eta`), in which the `1/V` terms of `(kappa_2[u]-1)/kt_2` and of `(1-G(x))/theta`
cancel symbolically, and whose constant term is `4/5` exactly. This is a finite computation on
the fast engine (the cumulants are Faulhaber sums in `V` with polynomial coefficients in
`eps = V eta`, the Edgeworth weights pair every half-power of `eta`).

---

## Head of 20:36, kept

**Two finite certificates for the DIRECT statement `g(n,t) := -Delta^2 log p_t > log(1 + 4/(5n))`
on `t >= 628, n >= 1261`:** (20:36, 2 September 2026, clock measured)

    (II) sparse, t*theta <= 2 (N >= t^2/2):   the exact finite sampling expansion of p_t, and the
         inequality  Xi(a, b) > 0  for the rational function  Xi = (N Phi - 4/5)/(a b),
         a = 1/t, b = t^2/N, Phi = -(log F)'' , on the box  a <= 1/627, b <= 2.01.
    (I)  dense, t*theta >= 1:                  Edgeworth for kappa_2[phi] through weight 9 with a
         crude weight-10 remainder (cushion 10^4 at the sparse end), cumulants by polygamma with
         a Stirling remainder, certificate in (v, eps).

**Why this is smaller than the head of 19:52 (the `Delta^4` / RLC chain), not merely different.**
Measured this hour, exact arithmetic:

- The RLC route needs `-N^3 Delta^4 log p = H''(theta)` with `H'' >= 2`; in the sparse regime the
  pieces of the fourth difference are `O(1/N^2)` and cancel to `O(1/N^3)` -- a cancellation of
  order `N/t^2`... measured: partial sums of `-N^3 Delta^4 log F` at `t theta = 4` swing through
  `10^6` before settling (lab/sampling_expansion.py). It would need Edgeworth through weight 13
  and the binomial half expanded too.
- The direct route's leading term is the proved `H(theta) >= 4/5 + 1.0077 theta` (variance
  identity), the binomial half is exact (`Delta^2 log C(N,t)` in closed form), and the margin is
  `(176/175) t / N^2 + (1.5 + 0.9 t theta)/N^2`, positive and growing in `t`, with no
  cancellation beyond a factor `t` among `O(1)` pieces. Measured on 12 cells `t = 628..2000`,
  `t theta = 0.25..2` (lab/sparse_direct.py): `-N^2 Delta^2 log F = H(theta)` to six digits,
  `margin * N^2 / t = 1.0065..1.0109`.
- Edgeworth for `kappa_2[u]` through weight 10 is derived EXACTLY on the fast engine in 1.0 s
  (lab/edgeworth_exact.py, sparse fmpq polynomials with weight truncation; the sympy version ate
  2 GB and was killed) and matches the exact second cumulant of the fixed-tilt Fourier weight to
  `5e-14` at `N = 400`. Its weight-2 term `12 c_4 - 36 c_3^2 + 6 c_1 c_3` is `(k_2/2)(log k_2)''`,
  the saddle-point `-(1/2) log V` term, recovered.
- The sampling expansion is an exact finite identity: `g` from it agrees with `g` from the exact
  `e_t` to `1e-118` (three cells). Its tail is bounded by `(2e s_2 t theta / j)^{j/2}` from
  `|(1+w)e^{-w}| <= e^{|w|^2}` on `|w| <= 1/2`.
- The suprema of the Bernoulli cumulant polynomials, needed for every crude remainder, are
  `sup |kappa_j(q)| <= 2.5 (j-1)!/pi^j` for `j <= 13` (exact roots, flint), and `j! 2.7/3^j` for
  all `j` by Cauchy on `|s| = 3`.

**Not done:** the certificate builders themselves -- (II) first (rational function in `(a,b)`,
exact division by `ab`, interval sweep), then (I).

Everything below this line is the state before the direct route was taken up.

---

## Head of 19:52, kept

**An explicit remainder for the fourth `t`-derivative of `log P(t)`, where `P` is the
Poisson-binomial mass function of the tilted spectrum, uniform in `theta in (0, 1/2]` for
`N >= 1260`, `t >= 628`, with constant below `2N`.** (19:52, 2 September 2026, clock measured.)

Everything else in the chain `A1 + A2 + B + C => conjecture` is either proved or is now an
explicit function of two variables with a known expansion. What moved since 17:45, and why it
is smaller and not merely different:

1. **Every tilted cumulant is a polygamma function of complex argument.** With `a = m + 1/2`,
   `y = 1/(2 sqrt r)`,

       x = k_1 = 2 [ m + y Im psi(a + i y) - (pi y / 2) tanh(pi y) ] ,
       k_{j+1} = -(y/2) d/dy k_j .

   Checked against the exact sums to 80 digits at `m = 1, 3, 10, 40, 200` (scratch script
   `digamma_cumulants.py`; the literature agent found the same identity independently, to
   `1e-31` in nine cells). So **item 1 (the cumulant inequality) and item 2 (the `(1/2) log V`
   term) are explicit functions of `(m, v)`**, and the Stirling series of `psi` carries an
   explicit remainder (DLMF 5.11.11) uniform in `|arg z| <= pi/2`, which holds here since
   `Re z = a > 0`. No midpoint rule, no quadrature tail. The expansion in `1/m` at fixed `v`
   reproduces the midpoint-rule boundary terms (checked at leading order: `theta = theta_0(v)
   + O(1/m^2)`, the `O(1)` terms cancelling exactly).

2. **The `t`-derivatives are exact, so no remainder is ever differentiated.** Fix `r` for the
   whole stencil. `P_r(t) = (1/2 pi) INT E(r e^{i phi}) e^{-i t phi} d phi / E(r)` is smooth in
   real `t`, agrees with the mass function at integers, and its fourth logarithmic derivative
   is `kappa_4[phi]`, the fourth cumulant of `phi` under the complex weight
   `exp(K_r(phi) - i t phi)` on `[-pi, pi]` -- an algebraic combination of the five moments
   `M_j = INT phi^j exp(K_r - i t phi) d phi`. Then `Delta^4 log p = INT^4 kappa_4[phi]` with
   no remainder, and each `M_j` is a Laplace integral whose expansion has an explicit
   remainder. Leading term: `-(1/k_2)''(x) = -(3k_3^2 - k_2 k_4)/k_2^5`, i.e. item 1 -- the
   same object seen from the dual side.

3. **The sub-obstruction that was never named: uniformity as `theta -> 0`.** All of yesterday's
   cells had `theta >= 0.125`. The open region is `t >= 628, N >= 1260`, and at fixed `t` it
   contains every `N`, so `theta -> 0` inside it. There the spectrum and the binomial are both
   nearly Poisson with mean `t`, each `kappa_4 ~ -2/t^3`, and their difference is `-H''/N^3`.
   A bound on each half separately can never reach the difference (relative accuracy
   `theta^3` would be required); the two halves must be treated by the identical procedure,
   ERR-0033 for the fourth time. Two regimes, both finite:

       (beta)  theta in [0.2, 1/2], N >= 1260:  Edgeworth expansion of `M_0..M_4` for each
               half to relative order `1/N`, explicit remainder `C/N^2`; the split is
               harmless here because the required relative accuracy is `theta^3 >= 0.008`.
       (gamma) theta in (0, 0.2], t >= 628:     exact perturbation of the spectrum around the
               binomial weight. With `dK = K_spec - K_bin = SUM_j (-1)^{j-1}(S_j - S_j^bin)
               (e^{i phi} - 1)^j / j`, `S_j = 2 SUM q_k^j`, one has `|dK| <= (2 theta/5)
               |log|e^{K_bin}|| (1 + O(theta))`, so `log E_w[e^{dK}]` expands in cumulants of
               `dK` with every moment `E_w[(e^{i phi}-1)^J]` an exact finite combination of
               binomial coefficients, and `Delta^4` acts on those exactly. Geometric
               convergence with ratio about `theta`, uniform as `theta -> 0`.

4. **Literature verdict, negative, recorded** in
   `results/LITERATURE_VERDICT_RATIO_LOG_CONCAVITY_2026-09-02.md`: no theorem gives a third
   difference of `log` of normalised coefficients for any class of real-rooted polynomials.
   Structural reason found there: coefficient reversal preserves real-rootedness and negates
   RLC (index `i -> N-1-i`), so no reversal-invariant theorem (Newton, Rosset, Sylvester,
   Jensen hyperbolicity) can ever give it. The observed unique failure at `i = N-2` is the
   image of `i = 1`.

**Not done, in one sentence each:** the Stirling-series form of item 1 with its `1/m^2`
coefficient and remainder; the Edgeworth remainder constants for `kappa_4[phi]` in regime
(beta); the perturbation-series bounds in regime (gamma); the assembly.

---

## Previous head (the ERR-0044 correction), kept

**`(3 k_3^2 - k_2 k_4) / k_2^5` is larger for the spectrum than for the binomial**, for every
`x` on the theorem's range, where `k_j` are the cumulants of the tilted vector
`q_k = b_k r(x)/(1+b_k r(x))` at the saddle `sum q_k = x`.

Four cumulants of an explicit probability vector against the constant vector. Derived from
`dk_j/dx = k_{j+1}/k_2`, which follows from `dq_k/dr = q_k(1-q_k)/r`. Verified in 16 of 16 cells;
`N^3` times the difference converges to `H''(theta)`, proved positive on all of `[0, 1/2]`.

**Corrected 17:45 (ERR-0044).** This file previously said "everything else is proved". It is
not: two links of the chain are measured, not proved --

    A2  ->  G'''' < 0            ratio 0.99995     (the (1/2) log V term in the saddle F)
    Delta^4 G  ->  Delta^4 log p ratio 1.000005    (the saddle expansion's own error)

So the open work is **A2 plus those two remainder bounds**. What IS proved, unconditionally:
`H'' >= 352/175 > 2` on all of `[0, 1/2]`, `r_1 > 1` for every `m`, `M(n,1) > 4/5` for every
odd `n >= 5`, and the certificate ladders (`t <= 627`, `n <= 1259`).



## Honest bookkeeping, 16:44 -- the obstruction did NOT move this evening

Since 16:44 the head has read the same sentence: **prove RLC at every index**. What the evening
produced, and what it is worth:

**Real, and not movement.** The rung factors as a trivially positive linear part times ONE
irreducible of degree `8i-3` -- half the rung. That halves the per-rung work and describes the
object better. It does not make the infinity smaller: there is still one new irreducible
polynomial per index.

**Two routes killed, which narrows the search and moves nothing.**
- a closed form for the leading coefficient of the hard factor -- the 2-part is exactly
  `2^(8i-6)`, the odd part is not a fixed product (new primes at `i = 4, 6, 7`, exponents not
  monotone). Dead.
- monotonicity of `r_i` in the index, which would have reduced RLC everywhere to RLC at one
  index. Not monotone at any of five sizes. Dead.

**One overstatement of mine corrected** (ERR-0041): RLC does not hold on the full index range;
it fails at `i = N-2` at every size, exactly one index, always the last. Sharper than what I
wrote, and irrelevant to the theorem, which needs `i <= (N-1)/2 - 1`.

**Named as required by the law:** this evening's work after 16:44 is structure and dead ends.
The blocking sentence is unchanged. The ladder climbed from 89 to 132 rungs, which is coverage,
not closure.

*Previous head, kept:* a one-sided bound `N^2 r >= -706 110` on `t >= 591, N >= 1186`; measured
`-0.09515` at the corner, slack 417 816. That route is still open and still valid.

**A one-sided bound `N^2 r >= -706 110` on the region `t >= 591, N >= 1186`.** That region is
what the first 590 certificates left open when the bound was computed; the ladder has since
passed 601, which only shrinks the region further and weakens the requirement, and on it `H(theta) >= 4/5 + (176/175) t/N` (proved slope,
`min H' = 1.007723`) supplies the cushion. Measured: `N^2 r >= -1.69` for every `n >= 203`, a
slack of 417 816. The same requirement read at the smallest size in the problem gave a threshold
of 5.41 against 2.86 -- which is why a near-sharp analysis seemed necessary, and is not.

**A uniform bound `|N^2 r| <= 1070` for `N >= 1186`.** Any bound of that size closes the
conjecture -- the constant is allowed to be 635 times worse than the measured value (1.69).
Reading the same requirement at the smallest size gave a threshold of 5.41 against 2.86, a
slack of 1.9, which is why a near-sharp analysis was needed; the 590 unconditional
certificates make `n > 1185` free, and the threshold grows linearly in `N`.

**OPEN, 2 September 15:06.** One bounded number:

    |N^2 r|  <=  2.86     for every n >= 7 ,
    r := n g - H(theta) - [H(theta) + (1/2)R''(theta)]/N .

The conjecture needs `|N^2 r| <= 0.9021555 N`, an allowance that grows linearly while the measured
`|N^2 r|` **decreases** in `n` (`2.856, 2.463, 2.263, 2.142, ... , 1.686` at `n = 7 .. 201`). The
margin runs from `1.90` at `n = 7` to over `100` by `n = 201` and keeps growing. `n = 5` is the sole
exception and is proved directly; `n = 3` is where the conjecture is false.

Checked exhaustively: every index of every odd `n` from 7 to 501 — **31 372 cells, no failure**.

Everything else is proved: `H >= 4/5` (twice today), `(1/2)R'' >= 0.1021555` (majorant, Cauchy tail
`9.8e-32`), `M = n(e^g-1) > n g`, and the range `5 <= n <= 1139` by exact certificates. The worst
cell is `t = 1`, which is itself proved unconditionally.

See `results/THE_UNIFORM_EXPANSION.md`.


**OPEN, 2 September 14:59.** One bounded remainder:

    |N^2 r|  <=  0.9021555 . N ,    r := n g - H(theta) - [H(theta) + (1/2)R''(theta)]/N .

**MEASURED: `|N^2 r| <= 1.4264`** across the regime (it even changes sign, `+0.095` at
`theta = 1/2`). The requirement is satisfied for every `N >= 1.58`; at the edge of the certified
range (`n = 1119`) the allowance is `1009` against `1.43` — **a margin of about 705**.

Everything else is proved: `H >= 4/5` (twice today), `(1/2)R'' >= 0.1021555` (majorant with a
Cauchy tail of `9.8e-32`), and `M = n(e^g-1) > n g`. See `results/THE_UNIFORM_EXPANSION.md`.


**OPEN, 2 September 14:55.** A rigorous bound on the remainder of the **uniform** expansion

    n g(n,t)  =  H(theta)  +  [ H(theta) + (1/2)R''(theta) ] / N  +  remainder ,
    g = -Delta^2 log p ,  theta = t/N ,  N = n-1 .

**Both leading coefficients are PROVED positive**: `H >= 4/5` (proved twice today) and
`(1/2)R'' >= 0.1021555` (majorant, Cauchy tail `9.8e-32`), so `H + (1/2)R'' >= 0.9021555` and

    n g  >=  4/5 + 0.9021555/N - |remainder| .

The conjecture needs `|remainder| < 0.9021555/N`, i.e. a `1/N^2` coefficient below `0.9 N`.
**The measured coefficient is about `0.90`** — satisfied for every `N >= 1`, with a factor of about
`1120` to spare at the edge of the certified range.

The formula is matched to six digits at five values of `theta` and reproduces the derived fixed-`t`
expansion `n g = 4/5 + (176t+160)/(175N)` in the limit `theta -> 0`. See
`results/THE_UNIFORM_EXPANSION.md`.


**OPEN, 2 September 14:48.** A rigorous bound on the remainder of

    n g(n,t)  =  4/5  +  (176 t + 160)/(175 N)  +  remainder ,     g = -Delta^2 log p .

The conjecture needs `|remainder| <= (176t+160)/(175N)`, i.e. a `1/N^2` coefficient below `1.92 N`
at `t = 1`. **The measured value is about `1/3`** — a margin of roughly `5.8 N`, which at `N = 1120`
(where the certificates stop) is a factor of about `6500`.

**Why this is the right target.** The leading term is `4/5` exactly, so any bound correct only to
leading order yields `M >= 4/5` without strictness; the conjecture lives in the next term, and that
term is now explicit: `c_1(t) = (176t+160)/175`, MEASURED to between `4e-07` and `5e-04` at eight
values of `t`, with slope `176/175` — the proved minimum gap of the limiting inequality, one gap per
index. See `results/THE_SECOND_ORDER_IN_CLOSED_FORM.md`.


**OPEN, 2 September 14:41.** One explicit inequality, with no construction in it:

    n g(n,t)  >=  (1 + 1/N) . [ 4/5 - 12/(5(4m^2-1)) ] ,
    g := log[ p_t^2/(p_{t-1}p_{t+1}) ] ,   N = n-1 = 2m .

Left: an ordinary second difference of `log p`. Right: an explicit rational number.

**It gives the conjecture in two elementary steps**: `M = n(e^g-1) > n g`, and the right-hand side
exceeds `4/5` exactly when `4m^2 - 1 > 6m`, i.e. `m >= 2` (pure algebra). It arises by letting the
free parameter `r` of the Poisson-binomial identity go to its endpoint, where the constant becomes
`Var(b)/bbar^2` — the spectrum's own relative spread.

**MEASURED at every index of eight sizes**, minimum always at `t = 1`, margin `~1.12/N` against a
threshold excess of `~0.8/N` — a factor of **2.4** at the tight point. See
`results/THE_CONJECTURE_AS_ONE_EXPLICIT_INEQUALITY.md`.


**OPEN, 2 September 14:15.** One deviation bound:

    | N^2 (E + g^2/2)  -  limit |  <=  0.4222     for  N > 1056 ,

where `E := g - (sigma_bin^2-sigma^2)/(sigma^2 sigma_bin^2)`, `g := -Delta^2 log p`, and the limit
is `(1/2)R'' + H^2/2 >= 0.1021555 + 0.32 = 0.4222`, **both constants PROVED**.

**The measured deviation is 0.0071. Margin: about 59.**

Everything else is done: `M = n(e^g-1) = P + n(E+g^2/2) + n g^3/6 + ...` with every later term
positive; `P > 4/5` from the PROVED limit inequality plus (C) (a pure condition `t >= 3`, uniform
in `m`); and the exact certificates cover `5 <= n <= 1057`. `E + g^2/2 >= 0` was checked at every
index of every size from `n = 5` to `501` with **no exception** (global minimum `+0.17448` at
`n = 5`).

See `results/THE_THEOREM_IS_A_VARIANCE_INEQUALITY.md`.


**OPEN, 2 September 14:10.** One statement, in its weakest sufficient form:

    (B)   E  >=  - g^2/2 ,      E := g - (sigma_bin^2 - sigma^2)/(sigma^2 sigma_bin^2) ,
                                g := -Delta^2 log p .

In scaled terms: `N^2 E >= -0.322`. **The measured value is `+0.105`** — so the bound may be wrong
by `0.427` in absolute terms and the conjecture still closes. (`M = n(e^g-1) = P + nE + n g^2/2 +
...`; the `n g^2/2` term is positive and had been discarded by using `M >= n g`.)

Everything else in the chain is done: the limit inequality is **PROVED**, `(C)` is settled as a
pure condition `t >= 3` uniform in `m`, `(B)`'s leading order is **PROVED**
(`R'' >= 0.2043110084`), and the exact certificates cover `5 <= n <= 1057`.

See `results/THE_THEOREM_IS_A_VARIANCE_INEQUALITY.md`.


**OPEN, 2 September 14:01.** One statement:

    (B)   g  >=  (sigma_bin^2 - sigma^2)/(sigma^2 sigma_bin^2) ,      g = -Delta^2 log p .

Everything else in the chain is now done. `(A)`, the finite variance inequality
`V(5n+4NQ) > 4NQ^2`, follows from the **PROVED** limit inequality plus `(C)`, and `(C)` is settled:
the exact midpoint bound `|INT f - (1/m)SUM f(y_k)| <= (1/24m^2) INT |f''|` gives a threshold
`m > 2 coef/(5V)`, which in the variable `t` reads `m/threshold = 0.4 t` **independently of `m`** —
so it holds for every `t >= 3`, and `t <= 2` is covered by the exact certificates.

`(A)` and `(B)` give `M = n(e^g-1) >= n g >= P > 4/5`.

`(B)`'s leading order is PROVED (`R''(theta) >= 0.2043110084 > 0` on `[0,1/2]`); what is missing is
`(B)` exactly. See `results/THE_THEOREM_IS_A_VARIANCE_INEQUALITY.md`.


**OPEN, 2 September 13:58.** Two "next order" statements, both MEASURED:

    (B')  E := g - (sigma_bin^2 - sigma^2)/(sigma^2 sigma_bin^2)  >=  0   exactly
          (its leading term (1/2)R''/N^2 is now PROVED >= 0.1021555/N^2)
    (C)   the 1/m^2 midpoint-rule terms are dominated by the +5V/(2m) leading correction
          (that leading term is verified positive at 7 values of v and 3 of m)

With those the conjecture follows: the finite variance inequality `V(5n+4NQ) > 4NQ^2` comes from
the **PROVED** limit inequality plus (C), and it with (B') gives `M = n(e^g-1) >= n g > 4/5`.

**PROVED today:** `H(theta) > 4/5` by a short new route (series + 359 637 boxes + analytic tail,
`release/scripts/variance_limit.py`, exit 0); `4/5 = lim Var(b)/bbar^2` for the centred-square
spectrum; `R''(theta) >= 0.2043110084 > 0` on `[0,1/2]`; `c_1 >= 0.993584 > 0` with
`c_1(0) = 216/175`; `L_i = (176/175) 2^(8i+10)/3^(4i+6)`; the dyadic rungs `t = 1..16`.

See `results/THE_THEOREM_IS_A_VARIANCE_INEQUALITY.md`.


**OPEN, 2 September 13:55.** Two elementary statements, both MEASURED, neither proved:

    (B)  g  >=  (sigma_bin^2 - sigma^2)/(sigma^2 sigma_bin^2) ,   g = -Delta^2 log p
    (C)  the 1/m^2 midpoint-rule terms are dominated by the +5V/(2m) leading correction

With those, the conjecture follows: `(A)` — the finite variance inequality `V(5n+4NQ) > 4NQ^2` —
comes from the **PROVED** limit inequality plus `(C)`, and `(A)+(B)` give
`M = n(e^g-1) >= n g >= P > 4/5`.

**PROVED today (2 September):** the limit inequality itself, i.e. `H(theta) > 4/5`, by a second and
much shorter route — an exact series on `(0,0.2]` with leading term `(352/945)v^6`, 359 637
certified interval boxes on `[0.2,3.5]`, and an analytic argument for `v >= 3.5`
(`release/scripts/variance_limit.py`, exit 0). And the constant is explained: `4/5` is
`lim Var(b)/bbar^2` for the centred-square spectrum.

See `results/THE_THEOREM_IS_A_VARIANCE_INEQUALITY.md`.


**OPEN, 2 September 13:56.** An explicit constant `D <= 478` with

    | g(n,t) - H(theta)/N |  <=  D / N^2 ,     g = -Delta^2 log p ,  N = n-1 ,  theta = t/N .

**Measured value of the constant: about 0.11.** (`lim N(Ng - H) = (1/2)R''`, measured at
`+0.113, +0.112, +0.109, +0.106, +0.103` for `theta = 1/16 .. 8/16`, agreeing with `(1/2)R''` to six
digits.) **Margin: about 4300.**

**Why it closes the conjecture.** `M = n(e^g - 1) >= n g`, so `M >= H(1 + 1/n) - D/n`, and
`H - 4/5 >= (176/175) theta` gives `M - 4/5 >= (176/175) t/N + (H - D)/n`. The exact certificates
cover `t <= 474`, so only `t >= 475` needs it, and there `(176/175) . 475 = 477.7`, hence any
`D <= 478` suffices.

This is the simplest and most generous form the obstruction has taken: **one two-sided bound on how
far the second difference of `log p` sits from its limit**, on a quantity that is elementary and
finite at every `n`, with three orders of magnitude of slack.


**OPEN, 2 September 13:38.** An explicit constant `C <= 477` with

    | M(n,t) - H(t/(n-1)) |  <=  C / n        for every n and every t in the regime.

That single crude bound closes the conjecture: `H - 4/5 >= (176/175) theta` is proved, the exact
certificates cover `t <= 474`, so only `t >= 475` needs it and there `(176/175).475 = 477.7`.
**The measured deviation is `n|M-H| <= 3.70` over the whole regime at seven sizes — a margin of
129.**

Nothing about the `c_1` expansion is needed for the closing argument; today's work on `c_1` is what
fixes the sign of the deviation (`c_1 >= 0.9936 > 0`, proved) and predicts its size
(`c_1(0) = 216/175` exactly). See `results/THE_PROOF_SKELETON_WITH_EVERY_CONSTANT.md`.

This is the smallest the obstruction has been: one two-sided envelope, on a quantity measured
everywhere, with two orders of magnitude of slack.


**OPEN, 2 September 13:17.** A bound `|M(n,t) - H(theta) - c_1(theta)/n| <= K_max/n^2` with an
explicit `K_max <= 1287`.

That closes the conjecture. `H >= 4/5` is PROVED; `c_1 = H^2/2 + H + (1/2)(log[s/(theta(1-theta))])''`
is a closed form matched to 6 digits at 8 values of theta, **and every one of its terms is
positive**; exact certificates already cover `n <= 953`, which is what makes `1287` the threshold.
The measured second-order coefficient runs `0.51 .. 6.42`, so **the bound may be crude by a factor
of about 200**. See `results/THE_FIRST_FINITE_SIZE_CORRECTION_IN_CLOSED_FORM.md`.

This is SMALLER than every previous form of the obstruction: it is one remainder bound, on a
quantity whose leading term is now known in closed form and positive, with two orders of magnitude
of slack.


**OPEN, 2 September 12:38.** An explicit constant `C <= 466` with

    M(n,t) >= H(t/(n-1)) - C/n        for every n and every t in the regime.

That single bound closes the theorem: the bulk is then handled by `H - 4/5 >= (176/175) theta`
(proved) and the edge `t <= 464` by the existing exact certificates. **The measured value of the
constant is about 3.7**, so the bound may be crude by a factor of 125 and still suffice. See
`results/THE_CROSSOVER_ARGUMENT_THE_REMAINING_BOUND_MAY_BE_CRUDE.md`.

This is SMALLER than the previous obstruction ("every coefficient of the shifted certificate is
non-negative", an infinite family): it is one inequality with an explicit, generous constant, and
the quantity it bounds was measured in `[1.238, 3.70]` over the whole regime at seven sizes.

Independent second route, measured only: the dyadic refinement `(n,t) -> (2n-1,2t)` preserves theta
exactly and `M` decreases along it in 2232 of 2233 pairs (the exception being n=5, t=1). If that
holds, `M > H >= 4/5` with no remainder bound at all.


**OPEN, 2 September 12:32.** Every coefficient of `Z_i(n + c(i))` is non-negative, for every `i`.

What is now PROVED, and no longer part of the obstruction: the top coefficient,
`L_i = (176/175) . 2^(8i+10)/3^(4i+6) > 0` for every `i` at once, derived from four constants
`f_1=8/3, f_2=-16/5, f_3=128/21, g_1=-2/3` of the cumulant expansion plus the binomial
normalisation. See `results/PROVED_THE_LEADING_COEFFICIENT_AND_THE_ORIGIN_OF_176_OVER_175.md`.
The same derivation gives `a_i=(4/3)^i` and `u_i=-i(i-1)/5`, and explains why two orders cancel:
`a_i^2 = a_{i-1}a_{i+1}` is Newton's equality case.

WITHDRAWN the same hour: "log-concavity of the coefficient sequence plus positive endpoints implies
all coefficients positive". False -- `1,-1,1` is a counterexample. The log-concavity measurement
(i = 0..24) stands as a measurement; the implication does not.


**CONJECTURED, 2 September 12:25.** The coefficient sequence of the shifted certificate
`Z_i(n + c(i))` is log-concave, for every `i`.

With that, the certificate follows: log-concavity puts every coefficient above the chord between
the two end values, the top end is `L_i = (176/175) . 2^(8i+10)/3^(4i+6) > 0` (closed form, exact at
`i = 0..20`, positive for every `i` at once) and the bottom end is `Z_i(c(i)) > 0` (verified
`i = 0..10`). Log-concavity itself is verified at `i = 0..24` with no exception.

It does NOT come from real-rootedness — `Z_i` has only about half its roots real or negative — so
Newton's inequality is unavailable and this is a property of the family, not a standard theorem.

See `results/THE_CERTIFICATE_REDUCES_TO_LOG_CONCAVITY_OF_ITS_OWN_COEFFICIENTS.md`.

This REPLACES "every coefficient is at least the leading one", which was the same statement in a
weaker form; the log-concave version has both endpoints already handled.


**CONJECTURED, 2 September 12:19.** For every `i`, every coefficient of `Z_i(n + c(i))` is at least
the leading one, `L_i = (176/175) . 2^(8i+10) / 3^(4i+6)`.

`L_i` is now a closed form, exact at `i = 0..20`, and positive for every `i` at once. "Leading is
the minimum after the shift" holds at every index checked (0..20). The certificate also always
carries the positive factor `(n+1)(2n+1)^4`, confirmed by exact division at `i = 0..12`.

The same constant `176/175` is the proved minimum gap of the limiting inequality, so the finite
certificates and the limit shape now agree on one number computed two independent ways. See
`results/WHY_176_OVER_175_THE_LEADING_COEFFICIENT_IN_CLOSED_FORM.md`.

This REPLACES the induction-step obstruction, which died by counterexample the same hour.


**CONJECTURED LEMMA (the induction step), 2 September 12:12.** Let `b_1 < ... < b_m` have a doubled
product that is ratio log-concave in the bulk. If `b_{m+1}/b_m <= 3/2`, the doubled product over
`b_1 .. b_{m+1}` is ratio log-concave in the bulk. Our step is `((2m+1)/(2m-1))^2 <= 3/2` for
`m >= 5`, and the base is proved, so this lemma closes the theorem for every `n`. The hypothesis
`P` is not yet correctly identified: "RLC in the bulk" alone is probably too weak.

Measured across six unrelated base families; the minimum threshold rises (1.594 -> 1.708 over
`m = 6..20`) while our required step falls (1.397 -> 1.105). See
`results/THE_INDUCTION_STEP_HAS_A_UNIVERSAL_THRESHOLD.md`.

This REPLACES the previous obstruction (a uniform remainder bound on the saddle expansion). The new
one is local: one appended root, one inequality.


*One sentence. Nothing else belongs in this file. Changing it is the only thing that counts as
progress on the theorem.*

## Stated

    (log p)''' < 0 pointwise on the regime, for an interpolant one can control.
    The triple-integral identity makes the pointwise form sufficient (a real gain -- interval
    arithmetic certifies pointwise, not averaged, statements), but the third derivative of the
    saddle interpolant still carries the Edgeworth term. Barrier relocated, not removed.

## Why it is the obstruction

Links 1-5 of the chain are proved (see `results/STATE_OF_THE_THEOREM.md`). This is the only
unproved link, and the requirement on it is now sharp: not accuracy, but accuracy *after* the
cancellation.

## Routes closed by test, so they are not tried again

**Recurrence-based ratio log-concavity criteria (closed 2 September 11:57).** The literature proves
ratio log-concavity for sequences satisfying linear recurrences of BOUNDED order, via bounding
functions. Tested directly for a P-recursion in `i`: orders 2, 3, 4 against coefficient degrees
2, 4, 6 -- **no recurrence exists in any of the nine cases**. The structural reason: for fixed `M`
the generating function gives `E' Q = P E` with `Q = PROD(1+b_k z)` of degree `M`, so the sequence
is P-recursive only at order `~M`, which grows. Criteria needing bounded order cannot apply.

**Jensen / Hermite / Petrov (closed 2 September 11:29).** The machinery proves hyperbolicity of the
Jensen polynomials, i.e. the Turan hierarchy. Tested on 600 random real-rooted spectra: RLC fails in
212 of them while Jensen hyperbolicity at d = 4, 5 fails in **none**. Jensen hyperbolicity is
strictly weaker than what is needed, so the route cannot reach the theorem -- and this explains why
no log-concavity literature has applied at any point: `Newton < higher Turan < Jensen at every
degree < RLC`. The statement is not in the Turan hierarchy.

## Where it is hard

It is no longer hard in kind, only unfinished: a remainder bound for a saddle-point expansion whose
integrand is a product of Gamma functions, where Stirling remainders are classical and explicit.
The two constants ahead of it are already computed and tiny.

## Detours since the last move, logged by name

**2026-09-02 09:11 -- the root-location route, which collapses back to the recorded obstruction.**
I proposed deriving the shift law from "all roots have Re < c", which would make it a classical
root-location question. For REAL `n` the complex roots are irrelevant -- a conjugate pair
contributes `|n - alpha|^2 > 0` -- so only the largest REAL root matters, which is exactly the
obstruction as recorded at 01:13. New data all the same: the maximum real part over ALL roots is
`~ i + 4` uniformly for `i = 0..18` (2.88, 5.32, 5.81, 7.20, ..., 22.00), which explains the size
of the shift `c(i) = i+2 / i+5` but does not advance the statement. Logged as a detour.


**2026-09-02 08:22 -- the leading term made exact.** `D' = kappa_3^bin/sigma_bin^6 - kappa_3/sigma^6`
is an identity at the finite-n saddle (verified 6 of 6 to eight digits), so the main term no longer
needs the limit shape or any asymptotics. The general Schur-type version of the resulting
inequality is false (1956 violations of 4000). This strengthens the foundation but does not by
itself move the blocking statement, which is still the bound on X. Logged as a detour.


**2026-09-02 04:32 -- both edge cases derived from one argument.** The single spectrum has no
exception because `ehat_0 = 1` fits its formula; the doubled one has exactly one because it does
not. Verified at four sizes each. Explanatory; the blocking statement is unchanged. Detour.

**2026-09-02 03:59 -- the self-convolution link.** The single odd-square spectrum satisfies ratio
log-concavity with NO failures at any index; the doubled spectrum, which is its self-convolution,
fails only at m = 3. This connects the theorem to the project's published counterexample showing
self-convolution does not preserve the property in general. Structural, and does not move the
blocking statement. Logged as a detour.


**2026-09-02 03:26 -- the far edge, derived but not moving the obstruction.** Three results were
recorded without the blocking statement changing, and the gate caught it:

- the reciprocal generating function `cosh^2(pi sqrt z/2)` and the exact edge ratio
  `(2m-3)^2/((2m-1)(2m-5))`;
- ratio log-concavity fails at **exactly one index, m = 3**, at every n tested;
- that exception derived in closed form as `9/10`, because `ehat_0 = 1` carries no doubling factor.

These explain the structure and pin the boundary of the phenomenon exactly. They do **not** weaken
`Delta g_i >= log(1 + 4/((2m-1)(2m-5)))`, which is still what must be proved. Logged as a detour,
not as progress.

## Last changed

2026-09-02 10:48   (the difference-versus-derivative barrier is removed by an exact identity:
                    Delta^3 log p is a triple integral of h', so pointwise h' > 0 suffices and
                    no remainder bound appears anywhere)

2026-09-02 05:43   (the object fully factored: a bracket with computed constants times
                    (1 + X/N^2) with X measured at 2.85..5.5, needing only |X| < N^2)

2026-09-02 05:08   (reduction to Delta^3 makes derivatives legitimate again; every constant of the
                    expansion computed -- max|R'''| = 0.051355, max|tau'''| = 4.124872 -- giving
                    positivity at every N >= 3 from three explicit terms)

2026-09-02 02:18   (the far edge computed in closed form via cosh(pi sqrt z/2); the resulting
                    elementary fraction turns out to predict the slack throughout the regime,
                    always from below)

## History of statements

- 2026-09-02 05:43  a bound |X| < N^2 on the residual coefficient (current). The object factors
  as a bracket times (1 + X/N^2); the bracket's positivity at every N >= 3 is established from
  computed constants, and X is measured at +2.85 to +5.5, converging in N, where 9 suffices at
  N = 3. See `results/THE_CONSTANTS_ARE_COMPUTED_AND_TINY.md`.

- 2026-09-02 05:08  bound the third derivative of the 1/N^2 remainder. Since the theorem
  needs only Delta^3, the mean value theorem applies at a fixed low order and the derivative route
  is legitimate. The expansion h' = H'/N^2 + R'''/(2N^3) - tau'''/N^4 - ... has every constant
  computed: H' >= 176/175 proved, max|R'''| = 0.051355, max|tau'''| = 4.124872, giving
  H' - |R'''|/(2N) - |tau'''|/N^2 > 0 at every N >= 3. See
  `results/THE_CONSTANTS_ARE_COMPUTED_AND_TINY.md`.

- 2026-09-02 02:18  Delta g_i >= log(1 + 4/((2m-1)(2m-5))). The reciprocal spectrum has
  the Weierstrass product cosh(pi sqrt z/2), the far-edge ratio collapses to
  (2m-3)^2/((2m-1)(2m-5)) exactly, and that fraction predicts the slack throughout the regime with
  the exact object always above it. See `results/THE_FAR_EDGE_IN_CLOSED_FORM.md`.
  (A retraction is recorded there: the appearance of 4/5 at m = 3 is outside the formula's domain
  and is not a derivation of the theorem's constant.)

- 2026-09-02 01:41  T_i >= 1 at finite n. The limit form is now PROVED outright:
  2Q'/Q^2 < 1/theta^2 - 1/(1-theta)^2 on [0,1/2], by a series with a 6.33e-70 tail near zero and
  2 163 485 certified interval boxes on the rest. The exact finite-n ratio T_i is never below
  2.012 while 1 suffices. See `results/THE_ANALYTIC_CORE_IS_ONE_ELEMENTARY_INEQUALITY.md` and
  `results/THE_EXACT_MARGIN_IS_A_FACTOR_OF_TWO.md`.

- 2026-09-02 01:30  the slack log B_i - delta_i is increasing in i. The theorem splits
  exactly into a spectral deficit and a binomial surplus; the slack runs from 1/n^2 at i = 0 to
  4/n^2 at the top of the regime, increasing at every index tested at six values of n. Its value
  at i = 0 is the proved base case. See
  `results/THE_SHARPEST_FORM_A_PROVED_BASE_AND_ONE_MONOTONICITY.md`.

- 2026-09-02 01:13  a root bound: largest real root of W_i <= i + 4.1. Running the
  ladder over i instead of k makes the degrees 8i+10 instead of 10, 29, 76, ..., 5111, and every
  member closes its own regime n >= 2i+5 with a margin growing like i. This proved t <= 32 at
  every n, up from t <= 8. See `results/THE_RLC_LADDER_IN_I_CLOSES_ITS_OWN_REGIME.md`.

- 2026-09-02 01:00  Delta^3 log p_i <= 0 for every i < N/2. The minimum of M(n,t) is at
  t = 1 for all 69 values of n tested, and t = 1 is proved for every n, so the theorem is exactly
  monotonicity in t -- one difference order at every i, instead of every order at one i. See
  `results/THE_THEOREM_NEEDS_ONE_ORDER_NOT_ALL.md`.

- 2026-09-02 00:46  non-negativity of every slice of the double series in (theta', 1/N).
  The identity N^(k+1) Delta^k g_0 = SUM_{j>=k} h_j D_{j,k} N^-(j-k), with D_{j,k} = Delta^k[(i+1)^j]_0
  >= 0, turns one property of H into positivity, monotonicity in n, and the limit at once.
  Confirmed against exact data to six figures at eleven orders, with c_1/c_0 = C(k+2,2) exactly as
  predicted. See `results/ONE_PROPERTY_GIVES_EVERYTHING.md`.

- 2026-09-02 00:35  monotonicity of N^(k+1) Delta^k g_0 in N. This removes the tail of
  the 1/n expansion from the proof entirely: limit positivity is settled, and one-sided
  convergence carries it to every finite n. Its k = 0 case is the already-proved A_1 > 0. See
  `results/THE_CLOSING_STRUCTURE_MONOTONE_IN_N.md`.

- 2026-09-02 00:30  non-negativity of the finite-n coefficients in combination.

  Why it weakened. `Delta^k[i^j] >= 0` for all i,j,k >= 0 (Stirling-number formula; 1521 cases
  checked, no negatives). Hence any G with non-negative Taylor coefficients has
  `Delta^k[G(i/N)] >= 0` for every k -- the differences never see the singularity that made the
  derivative comparison predict a 3.83^k blow-up. H's coefficients are h_0 = 4/5, h_1 = h_2 =
  176/175, ..., converging to 1, all 40 exact ones positive. See
  `results/THE_KEY_LEMMA_DIFFERENCES_NOT_DERIVATIVES.md`.

- 2026-09-02 00:25  Delta^m[S - S_bin] < 0, an inequality between two explicit closed forms
  (current). The approximation step behind it is now measured at the edge itself rather than
  extrapolated: N = 80 to m = 41 and N = 160 to m = 81, sign correct in 40 of 40 rows, relative
  error 2.461e-02 and 2.221e-02 at the two edges -- decreasing with N. See
  `results/THE_LEMMA_MEASURED_TO_THE_EDGE.md`.

- 2026-09-02 00:05  prove the measured  0.05 m^2/N^2  cancellation law

  Why it weakened. Approximating BOTH halves by the identical saddle procedure -- `S(i)` from the
  Gamma closed form, `S_bin(i)` from `(1+z)^N` -- makes the thirteen-order cancellation happen
  inside the approximation instead of after it. Sign correct in 39 of 39 cells at N = 80, 160,
  320 for m = 2..14; relative error 1.0e-4 where the previous assembly gave 1.1e+5. See
  `results/THE_CANCELLATION_MADE_ANALYTIC.md`.

- 2026-09-02 00:02  an expansion of the DIFFERENCE
- 2026-09-02 00:00  same sign as the closed-form saddle -- KILLED in three minutes.

  The sign form was true for `log e` alone: the correction to the Gamma-closed-form saddle is
  exactly the Stirling term `-1/(12 i)`, verified in 39 of 39 cells at N = 80, 160, 320 for
  m = 2..14, with the ratio converging to 1 (0.2 % at N = 320). And every term of a Stirling
  series carries the same sign under `Delta^m`, since `Delta^m[i^-k] = (-1)^m k! (positive)`
  independently of `k`. All of that stands.

  It does not transfer. The theorem lives on `log p = log e - log C`, where the two halves
  cancel to thirteen orders. Assembling the explicit pieces at N = 320 gives, against an exact
  `Delta^14 log p = -1.03e-23`, the value `-1.18e-18` -- a factor of 10^5 -- and at m = 9, 11, 13
  the assembly returns the WRONG SIGN.

  **What this rules out, permanently:** any route that expands `log e` and `log C` separately and
  subtracts. Three of the day's routes were of that kind. The expansion has to be of the
  difference.

- 2026-09-01 23:47  a uniform-in-m bound on the m-th derivative of the remainder

## Routes closed by test, 2 September

- **`log b` concave => doubled RLC (criterion C2).** DEAD. 8 failures in 1077 log-concave spectra,
  all geometric or near-geometric. Survived 776 random spectra first — the confirmations measured
  the generator.
- **Any soft class criterion.** DEAD by structure: geometric spectra fail, and their failure
  threshold in `q` falls with `m` (`q=3` fails by `m=16`, `q=5/2` by `m=26`, `q=2`'s margin decays
  to `+0.0004`). A sufficient condition must therefore be quantitative, i.e. the limiting
  inequality already proved.
- **The ladder as a route to all `n`.** Stopped by decision, not by failure: cost grows like `i^3`,
  so a 1000x speedup buys 10x in range. It proved `5 <= n <= 953`.
- **The step lemma** (*doubled base RLC + `b_{m+1}/b_m <= 3/2` implies the extension is RLC*). DEAD
  within two minutes of being stated: 6 counterexamples among 3311 bases satisfying both hypotheses,
  e.g. `b = (4,5,6,7)` with `rho = 13/10`. The failing margins are tiny (`-9e-5`, `-5e-6`), so the
  hypothesis is close to right and is missing one condition — but as stated it is false, and the
  universal-threshold table that suggested it measured six well-behaved families only.
- **Majorization (2 September 12:45).** Exact identity, worth keeping: since `log e_t = log P(S=t)
  + log E(r) - t log r` for any `r>0` and the last two terms are affine in `t`, the third difference
  kills them, so `Delta^3 log e_t = Delta^3 log P_PB(t)` for the Poisson-binomial `S` with
  `q_k = b_k r/(1+b_k r)`, and `C(N,t)` is the Poisson-binomial with constant `q = 1/2`. Hence
  **RLC is a comparison of `Delta^3 log` between a spread `q`-vector and the constant one.** The
  natural claim -- *`q` majorizes `q'` (same sum) implies `Delta^3 log P_q <= Delta^3 log P_q'`* --
  is DEAD: 5 counterexamples in the first 5 majorization pairs tested, e.g.
  `q = (2/5,1/5,1/2,1/5,1/2,3/5,1/5)` against its two-entry average. Spread does not beat uniform in
  general; our spectrum is special, consistent with everything else killed today.
- **Literature, central factorial numbers (2 September 12:45).** The r-central factorial numbers with
  even indices are known to be "strictly log-concave and Poisson-binomially distributed", which is
  our object's first difference, not its third. Verdict NEGATIVE, recorded so it is not searched
  again.
- **The bootstrap one floor up (2 September 12:56).** MEASURED: `W_t = (p_{t+2}^3 p_t)/(p_{t+1}^3
  p_{t+3})` is strictly increasing in `t` at every size and index tested (`n = 21..501`, ~860
  indices). Algebraically `W_t < W_{t+1}` is exactly `p_{t+1}^4 p_{t+3}^4 > p_{t+2}^6 p_t p_{t+4}`,
  i.e. **`Delta^4 log p < 0`** -- the next floor of the difference hierarchy. So
  `Delta^4 <= 0` plus `Delta^3 <= 0` at the single index `t = 0` (PROVED by the ladder) would give
  RLC everywhere. Logged as a DETOUR by this project's own rule: it replaces one infinite family by
  a strictly stronger one, and the escape upward (complete monotonicity / Hausdorff) is already
  recorded as closed -- the measure does not exist. Kept because it names the bootstrap exactly.
- **Newton forward-difference certificate (2 September 13:01).** `Z_i(m) = SUM_j Delta^j Z_i(m0)
  C(m-m0, j)`, so all forward differences non-negative at `m0` gives `Z_i(m) >= 0` for every integer
  `m >= m0` -- a strictly WEAKER requirement than the Taylor shift. Measured `i = 0..20`: for **even
  `i`** it holds at `m0 = 0` (through `i = 14`; then `2, 4, 5` at `i = 16, 18, 20`), i.e. no shift is
  needed at all, a gain of `+3..+8` over the Taylor shift. For **odd `i`** it gives exactly the same
  `m0 = (i+5)/2` -- **no gain, and odd `i` is therefore the binding case.** The same parity split
  appears in the Taylor shift law (`(i+2)/2` even, `(i+5)/2` odd), so the asymmetry is structural,
  not an artefact of one certificate form. `L_i` has no parity dependence, so whatever causes it
  lives in the lower coefficients.

## Displacement check on the afternoon's work (14:17)

**Was the obstruction moved, or re-described?** Both, and the honest split is:

**Moved (real):**
- the limiting inequality `H > 4/5` now has a **second, shorter proof** (359 637 boxes against
  2 163 485, no 40-term expansion of `H`, edge in closed form) — and the same computation
  **explains `4/5`** as `lim Var(b)/bbar^2`, which was unknown this morning;
- `R''(theta) >= 0.2043110084 > 0` — a new proved statement;
- `(C)` reduced to a pure condition `t >= 3`, uniform in `m`, by the exact midpoint bound;
- the remaining quantity is now **elementary** — variances of explicit probabilities — where this
  morning it was the third difference of a logarithm of elementary symmetric functions.

**Re-described (not moved):** stating the residue as *"a deviation bound with 0.4222 available"* is
**the same statement** as `E + g^2/2 >= 0`, because the limit is proved `>= 0.4222`. The `59x`
figure is an honest measure of how crude the bound may be, not a reduction of the statement.

**So the last step is still what it has been all day: a uniform asymptotic error bound.** What
changed is the object it is about and how much slack it has, not its kind.


---

*Timestamps in this file were repaired on 2026-09-02 19:33 +0200 by `tools/timestamp_repair.py`: the typed times drifted up to five hours ahead of the clock, and each was replaced by the time of the git commit that introduced it (ERR-0045).*

## Goal 2 opened (2 September 2026, 23:58): Sibuya's conjecture (1988, eq. 3.4) by the same machine

**The blocking statement.** For the spectrum {1..N}: `n(R_j - 1) >= n/(3n - j)` for all `j >= J0` and all `n`, where the true margin is only `2(j-1)/(9n)` (measured, results/sibuya_margin_probe_2026-09-02.txt): the sparse certificate must certify at least ONE THIRD of the true slope (`1/9` out of `1/3`), where the 4/5 certificate certified 14% (`0.1417` out of `1.006`). The dense regime has 3-60% of room and the exact V-series resolves it; the ladder (degree `2j-1`, shift `j+2`) closes `j <= J0` exactly and is running to `J0 = 2000`.

**What must shrink.** The sparse certificate's error budget (the `j > J` tail and the finite-m moment perturbation, both crude Cauchy bounds proportional to `a = 1/j`): with `J0 = 2000` the budget is 3x smaller by `a` alone; the rest needs sharper tail bounds, not a new language.

Goal 1 (`4/5`) has no mathematical obstruction left; its artifacts (merged ladder log, two validation reports) are being produced tonight.

## Goal 1 SEALED (3 September 2026, 01:14)

Every artifact of the theorem M(n,t) > 4/5 is in place and checked by a script, not by a sentence:
- the ladder t <= 627: results/m_ladder_log_2026-09-02_full.txt, 627 rungs, merged and checked by lab/ladder_log_merge.py (theorem.py re-runs the check as piece A);
- sparse_certificate_full.py: independently validated (six PASS);
- dense_certificate_a.py: independently validated (nine links PASS after the R4 fix; rerun 930 boxes, worst 0.800004);
- dense_certificate_b.py: independently validated (seven items PASS after the Ser repairs: base object G_j, exact degree d, real-point perturbation P, sigma-free corner tail);
- theorem.py --full: three VERDICT True, exit 0 (409 s); cross-check lab/crosscheck_certificates.py: the series equals the exact N g to 1e-19..1e-24, the boxes contain the exact value.
What "proved" still waits for: a human referee (the founder's decision on publication). No mathematical obstruction remains for goal 1.

## Goal 2 moved (3 September 2026, 01:44): the three Sibuya certificates all pass

The sparse certificate was the named obstruction (it had to certify a third of the true slope). It certifies
54% of it: `[a^3]/b^2 = 2/3` exactly, residual 0.357 after the lemmas (results/sibuya_sparse_certificate_2026-09-03.txt),
because `|mu_i| <= 1` (not `2^i`) and `tau >= 2000` shrink the error budget by more than the target grew.
Dense (a') passes with a thin margin `1.0e-5` at `theta = 0.1` (266 boxes, 95 s); dense (b') passes with
`S(0, zeta) = 1/3` exactly and `min(S - losses) = 0.3158` on 18 bands (I = 150 to control the Cauchy tail at v = 0.23).

**What blocks "proved" now (goal 2):** (i) the ladder to j = 2000 (940 done, six parallel ranges running);
(ii) independent validation of the three ports -- in particular the sparse lemma constants FMIN = 0.80,
CF = 0.35, C_i <= 3 + i^2/tau^2 beyond 60, which were derived by analogy and are marked TO VALIDATE;
(iii) the (b') rerun with the region widened by (1 + eps). No mathematical obstruction is named; the
remaining work is verification.

## Goal 2, the named obstruction (02:01): the top regime theta > 0.9 of Sibuya's (3.4)

Sibuya's (3.4) is stated for EVERY j <= n - 2; the spectrum {1..N} has no j <-> N-j symmetry (the validator of
dense (a') caught this). Certified tonight: the ladder j <= 1000 (running: 941..1000), and for j >= 1001 the
sparse regime (tau >= 1000, tau^2/N <= 1), the series regime (theta <= 0.1, j^2/N >= 1) and the boxes
0.1 <= theta <= 0.9 (two sweeps: [0.1, 1/2] and [1/2, 0.9]; the upper one has margin 0.198).

**The blocking statement:** n(R_j - 1) > n/(3n - j) for j >= 1001 and 0.9 < theta = j/N < 1 (N -> infinity),
i.e. for j' = N - j < 0.1 N missing roots. By e_j(1..N) = N! e_{j'}(1, 1/2, ..., 1/N) this is the Newton gap
of the RECIPROCAL spectrum {1/k} at the small index j', target 1 + 1/(2n + 1 + j'); the true excess there is
huge (exact: 0.74 at theta = 0.75, 7.2 at theta = 0.97; at j' = 1 it is (N H^(2) - H^2)/(N(H^2 - H^(2))) ~ zeta(2)/log^2 N,
against a target ~ 1/(2N)). Two sub-regimes: (T1) j' >= ~200: boxes in (eps, x), x = N/v the tilt parameter,
with q_k = k/(k + x), sigma^2 ~ x log(N/x) ~ j', the cumulants through psi^(q)(x + 1) directly (no Stirling
needed in that argument) -- the existing engine with a re-parametrisation; (T2) j' < ~200: the harmonic regime,
e_{j'}(1/k) as polynomials in H_N^(r), a new explicit-error lemma. Neither is a detour: the gap cannot be
closed by any of the four existing pieces.

## Goal 2 status (03:52): everything below theta = 0.9 is closed and machine-checked

sibuya_theorem.py --full: four certificates exit 0 with VERDICT True (results/sibuya_theorem_full_2026-09-03.txt),
the ladder log merged and checked (1005 rungs), all four ports independently validated
(results/VALIDATION_SIBUYA_*_2026-09-03.md). The obstruction is exactly the previous section: theta > 0.9, the
reciprocal spectrum {1/k}, sub-regimes (T1) and (T2). Nothing else stands between the work and Sibuya's full statement.

## Goal 2, the obstruction NARROWED (05:25): only the small-j' corner is left

An exploratory sweep of the box engine in (eps, v) over 0.9 <= theta <= 0.99 (results/sibuya_dense_a_top2_exploratory_2026-09-03.txt):
303 boxes ok, worst margin 1.03 (the target is ~0.5), 37 s -- wherever x = N/v >= 20 and j' = N - j >= 100 the existing
machine certifies Sibuya's target with room to spare, uniformly in N (eps -> 0 is a box edge). Not a certificate yet:
a sliver along eps v = 1/20 was skipped at depth 14 and the corner was excluded. **The obstruction is now exactly:
j >= 1001 with x = N/v < 20 or j' = N - j < 100 (the few-missing-roots corner), uniformly in N.** In that corner the
tilted distribution of the missing count is Poisson-like with mean j' ~ x H_N, the Stirling series at 1 + x is not
usable, and the Edgeworth needs sigma^2 ~ j' large: a genuinely different tool (the harmonic-number regime,
e_{j'}(1/k) as polynomials in H_N^(r)) -- sub-regime (T2) of the previous section; (T1) is closed in principle.

## (T2) route named (05:39): make Hwang's expansion effective ourselves

Literature (results/LITERATURE_VERDICT_STIRLING_SMALL_K_2026-09-03.md): no explicit-constant uniform asymptotic for
c(n, m) with m small exists; Hwang 1995 is uniform for m <= eta log n with O-terms, and its ONE non-effective step is the
remainder of Gamma(n+w)/Gamma(n+1) = n^(w-1)(1 + SUM pi_k(w)/n^k + R) uniformly on a disc |w| <= eta. That remainder is a
Stirling-series remainder on a disc -- the same object as phi_stirling (dense certificates), certifiable with the Bernoulli
lemma. So (T2) = Hwang's Theorem 1 with our explicit R, applied to the ratio at three consecutive m; the target gap is
~1/(2n) against an excess ~H^(2)/H^2, so first-order accuracy in 1/log n suffices once R is explicit. Adell 2022's
explicit remainder |E| <= 2(em/H)^3 is too weak (decides the gap only for m^3 <~ H/160). Recorded; not started.

## Goal 2, the top regime measured precisely (09:51)

Closed today: the corner j' = N - j <= 32 for EVERY N >= 1001 + j' (release/scripts/sibuya_harmonic.py, exact
incremental recursion E_k(N) = E_k(N-1) + E_{k-1}(N-1)/N up to N = 10^6, then bands in H with the exact rational
part of T, then the analytic tail at 2^48; log results/sibuya_harmonic_2026-09-03.txt, 271 s), and
0.9 <= theta <= 0.99 with j' >= 100 (the box sweep, margin 1.03).

**Two walls, both measured, not guessed.**
1. *The harmonic representation stops at j' ~ 32.* E_k(H) = SUM_m c_m H^{k-m}/(k-m)! is an alternating sum whose
   largest term exceeds the value by e^H k!/H^k; at j' = 40, N = 2*10^4 that is 1e10, so the P_r interval widths
   (~1/N^2) swamp the margin (measured: E_j' enclosure +/- 1.7e-5 against a value 2.8e-7). Raising the precision
   does not help (the width is not rounding); the cancellation disappears only when H >= j', i.e. N >= e^{j'}.
2. *The box engine stops at theta ~ 0.999.* N g contains (kappa_2 - 1)/kt_2 with kt_2 = j'/N -> 0 as theta -> 1,
   so the Edgeworth enclosure (~1e-4 at best) must be small compared with j'/N (measured: at theta = 0.99998,
   j' = 100 the enclosure of N g is +/- 9.1e3).

**The obstruction, exactly:** Sibuya (3.4) for j' >= 33 with theta > 0.99 (equivalently N > 100 j'), where neither
the alternating harmonic sum nor the Edgeworth box survives. What is needed is a stable evaluation of
e_{j'}(1, 1/2, ..., 1/N) for moderate j' and large N -- the Gamma-ratio series
Gamma(N+1+x)/(Gamma(N+1)Gamma(1+x)) computed as a power series in x at high precision, at a grid of N fine enough
for the monotone bracket E_k(N_a) <= E_k(N) <= E_k(N_b) -- or a genuinely new lower bound for the Newton gap of
the reciprocal spectrum {1/k}.

## Goal 2, the obstruction as it stands (12:06)

CLOSED today: the whole top of the row with at most 399 missing indices, at every N
(release/scripts/sibuya_corner_grid.py: exact recursion to 10^6, a third-order grid in N on exact Gamma-ratio
series to e^20, and the H-model certified for all 367 indices by one Taylor shift each; 504 s, log
results/sibuya_corner_grid_2026-09-03.txt).

**The obstruction, exactly: theta > 0.9 together with j' = N - j >= 400.** The harmonic instrument stops there
because the H-model's shift certificate fails around j' = 1000 at any starting height; the box instrument stops
because kt_2 = j'/N collapses on the sparse side N > j'^2, and the region has a cusp at (eps, Y) -> (0,0) that a
box scheme in those variables never resolves. The named fix: rerun the box machine in (mu, z) with mu = j'/N and
eps = mu^2 z, where the sparse side is the rectangle 0 <= z <= 1 rather than a cusp.

## Goal 2, the obstruction as it stands (12:06)

CLOSED today: the whole top of the row with at most 399 missing indices, at every N (sibuya_corner_grid.py:
exact recursion to 10^6, a third-order grid in N on exact Gamma-ratio series to e^20, then the H-model
certified for all 367 indices by one Taylor shift each; 504 s, log results/sibuya_corner_grid_2026-09-03.txt).

**The obstruction, exactly: theta > 0.9 together with at least 400 missing indices.** The harmonic instrument
stops there because the H-model shift certificate fails around 1000 missing indices at any starting height; the
box instrument stops because kt_2 collapses on the sparse side and the region has a cusp at the origin of its
variables. The named fix: rerun the box machine in the variables (mu, z) with mu the missing fraction and
1/N = mu^2 z, where the sparse side is a rectangle rather than a cusp.

---

## 3 September, after the debate: the obstruction shrank from 400 to 803 missing indices

**The obstruction, exactly: `theta > 0.9` together with at least 803 missing indices, at `j >= 1001`.**

Two things moved it. The window certified by the harmonic instrument doubled -- `sibuya_corner_grid.py` now
closes `33 <= j' <= 802` at every `N` (exact recursion to `10^6` from an artifact the verdict now reads rather
than asserts, a grid of 4705 steps to `N = 4.85e8` with worst relative margin `1.92e-10`, and the `H`-model
above that, all 770 indices by one Taylor shift each in 177 s). And the wedge `{theta >= 0.9, j' >= 803}` is
certified down to `Y = 1e-30` in `Y = 1/v`, which is a statement about most of the region but not all of it.

**What did NOT move, and is the whole difficulty.** At fixed `j'` with `N -> infinity` the two sides of the
inequality cancel at leading order -- the slack encloses zero to `+/- 1e-12` at `j' = 20000`. That corner is
where the certified statement stops being a statement about a margin and becomes a statement about which of two
equal leading terms wins. The Laguerre-Polya observation of the same day says the LIMIT there is classical
Newton and holds for every index; the finite-`N` factor `1 + 1/(2N+3+j')` is what remains, and no bound in this
package reaches it uniformly in `j'`.
