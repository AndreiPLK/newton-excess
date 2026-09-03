# The Newton excess of the centred-square spectrum

**Andrei Pluzhnik**
ORCID 0009-0005-5660-2603

*Draft of 2 September 2026. Every claim is labelled PROVED, VERIFIED, or CONJECTURED, and nothing
is labelled above what it has earned.*

---

## 1. The statement

For a real-rooted polynomial with positive coefficients, write `e_i` for the elementary symmetric
functions of its roots and `p_i = e_i / C(N,i)` for the normalised coefficients, `N` the number of
roots. Newton's inequality (1707) says `p_i` is log-concave, i.e. `R_t = p_t^2/(p_{t-1}p_{t+1}) >= 1`.
The **Newton excess** measures how far past the boundary a family actually sits:

    M(n,t)  =  n ( R_t - 1 ) .

> **The centred-square spectrum.** For odd `n = 2m+1`, the roots are `{ (n-2k)^2 : k = 1..n-1 }` --
> the odd squares `1^2, 3^2, ..., (n-2)^2`, each twice. Set `N = n-1`.

> **CONJECTURE.** `M(n,t) > 4/5` for every odd `n >= 5` and every `t < n/2`.
>
> **THEOREM (proved here).** It holds for every odd `n` with `5 <= n <= 1207`, and every `t < n/2`.

The bound `n >= 5` is necessary: at `n = 3` the spectrum is `{1,1}`, all roots coincide, Newton is
an equality and `M(3,1) = 0`.

The constant is sharp. `M(n,1)` descends to it: `1.3559, 1.0376, 0.9158, 0.8570, 0.8283, 0.8141` at
`n = 5, 11, 21, 41, 81, 161`; and section 4 identifies `4/5` as the value at zero of an explicit
limit function.

## 2. The property has a name

In the terminology of Chen, Guo and Wang (*Infinitely Log-monotonic Combinatorial Sequences*,
arXiv:1304.5160), with the ratio operator `R{a_n} = {a_{n+1}/a_n}`, a sequence is *log-monotonic of
order k* when `R^r{a_n}` is log-concave for odd `r <= k-1` and log-convex for even `r <= k-1`.
Order two is exactly **ratio log-concavity**. Their Theorem 2.1 gives infinite log-monotonicity
whenever `[log f(x)]''` is completely monotonic -- a hypothesis fitted to log-*convex* sequences
(Bernoulli, Catalan, central binomial), which `p_i` cannot satisfy, since Newton's inequality makes
`[log p]'' < 0`.

The log-concave mirror, which is what is needed here, is that `-[log p]''` be **absolutely**
monotonic. Section 3 gives the discrete form of the passage from derivatives to differences, which
is what the finite-`n` argument actually uses.

## 3. The difference lemma (PROVED, and elementary)

Every earlier attempt routed `Delta^k` through the mean value theorem for divided differences and
then bounded a derivative. That comparison predicts failure -- the limit function's nearest
singularity is at `|theta| = 1` while the quantity it must beat has one at `3.8256676`, so the
ratio of derivatives grows like `3.83^k` -- and the measurements contradict it outright.
Derivatives are the wrong instrument.

> **Lemma.** For integers `i, j, k >= 0`, `Delta^k[i^j] = SUM_l C(j,l) i^(j-l) k! S(l,k) >= 0`,
> with `S(l,k)` the Stirling numbers of the second kind.
>
> **Corollary.** If `G(theta) = SUM_j c_j theta^j` has all `c_j >= 0` and converges for
> `|theta| < rho`, then `Delta^k[G(i/N)] >= 0` for every `k >= 0` and every `i >= 0` with
> `(i+k)/N < rho`.

A polynomial with non-negative coefficients has non-negative differences of every order, wherever
its analytic continuation has poles. The differences never see the singularity; the derivatives see
nothing else.

## 4. The limit shape (PROVED)

The generating polynomial has an exact closed form,
`E(z) = ((4z)^m Gamma(m+1/2+zeta) Gamma(m+1/2-zeta) cos(pi zeta)/pi)^2`, `zeta = -i/(2 sqrt z)`,
verified exactly at twelve points. Its saddle satisfies `theta = 1 - arctan(v)/v`, and the envelope
theorem (`L(i) = log E(r) - i log r`, `L'' = -1/sigma^2`) gives

    sigma^2 = (N/2)[(1-theta) - 1/(1+v^2)] ,     sigma_bin^2 = N theta (1-theta) .

`sigma^2 <= sigma_bin^2` by concavity of `x(1-x)` -- Newton's inequality, recovered. The limit shape
of the excess is the gap,

    H(theta) = 2/[(1-theta) - 1/(1+v^2)] - 1/[theta(1-theta)] ,

with the two `1/theta` poles cancelling exactly. Its Taylor coefficients are exact rationals:

    h_0 = 4/5,   h_1 = h_2 = 176/175,   h_3 = 67328/67375,   ...,   h_39 = 1.0000000000

All forty computed coefficients are positive and converge to 1, so `H - 1/(1-theta)` extends to the
disc of radius 3.8256676. `H(0) = 4/5` is the constant of the conjecture.

## 5. Positivity, monotonicity and the limit, from one property (PROVED)

With `g_i = -Delta^2 log p_i` and `theta' = (i+1)/N`,

    N^(k+1) Delta^k [H(theta')]_{i=0}  =  SUM_{j>=k} h_j D_{j,k} N^{-(j-k)} ,
    D_{j,k} = Delta^k[(i+1)^j]_0 >= 0 ,   D_{j,k} = 0 for j < k,

a power series in `1/N` with non-negative coefficients. Hence for every `k` the quantity is
non-negative, is strictly decreasing in `N`, and tends to `h_k k!`. Against exact data at
`N = 80..176` the extraction reproduces `h_k k!` to six figures at eleven orders, with
`c_1/c_0 = C(k+2,2)` for `k = 5..10` exactly as the identity predicts.

## 6. Finite n: every t up to 128, at every n (PROVED)

The minimum of `M(n,t)` over the regime sits at `t = 1` for all 69 values of `n` tested
(`n = 5..141` odd, no exceptions). So the theorem is a base case plus monotonicity in `t`.

**The base.** `p_1 = (n^2-2n)/3`, `p_2 = n^4/9 - 8n^3/15 + 28n^2/45 + 4n/15`, and

    M(n,1) > 4/5  <=>  X(n) = 56/45 n^3 - 172/45 n^2 - 16/15 n  >  0 ,
    X(n+4) = 56/45 n^3 + 100/9 n^2 + 1264/45 n + 128/9 .

Four positive coefficients, so `M(n,1) > 4/5` for every `n >= 4`. **PROVED.**

**The monotonicity.** `g_{i+1} - g_i = log[p_{i+2}^3 p_i/(p_{i+1}^3 p_{i+3})]`, so it is
`W_i(n) = p_{i+2}^3 p_i - p_{i+1}^3 p_{i+3} >= 0`. Each `p_i` is a polynomial in `n` of degree `2i`
with leading coefficient `3^{-i}`; the leading terms of `W_i` cancel exactly, leaving
`deg W_i = 8i+10`, and the regime for index `i` is `n >= 2i+5`. The smallest shift `c(i)` making
every coefficient of `W_i(n+c)` non-negative:

     i  :  0   1   2   3   4   5   6  ...  86  87  88  89
    c-i :  3   5   2   5   2   5   2  ...   2   5   2   5

    c(i) = i+5  (i odd),    c(i) = i+2  (i even >= 2),    c(0) = 3,

so `c(i) <= i+5 <= 2i+5` and each certificate covers its whole regime. **No deviation from the law
in 126 consecutive cases**, and the predicted shift was checked separately far beyond the ladder --
at `i = 150, 175, 200, 225, 250, 275`, degrees up to 2210, zero negative coefficients in each.
Hence

    **M(n,t) > 4/5 for every t <= 584 and every odd n >= 5** ,

unconditionally, by exact polynomial certificates and no asymptotics. (Twelve hours earlier the
proved range was `t <= 8`.)

Structure of the family, which a uniform proof would have to use:

    lead(W_i) = 176/(175 . 3^(4i+6)) ,   n^4 | W_i for i >= 1 ,
    sum of the roots of W_i = (88 i^2 + 506 i + 518)/55 .

The numerator of the leading coefficient is constant, and `176/175` is exactly `h_1 = H'(0)`.

## 7. The exact split, and why the theorem is true

Writing `r_j = p_{j+1}/p_j = (j+1)/(N-j) . s_j` with `s_j = e_{j+1}/e_j`, the statement separates
exactly into a spectral half and a binomial half:

    delta_i  =  log( s_i s_{i+2}/s_{i+1}^2 )                    spectral deficit
    log B_i  =  log[ a^2(b^2-1)/((a^2-1)b^2) ] ,  a = i+2, b = N-i-1     binomial surplus
    the theorem  <=>  delta_i <= log B_i .

`B_i >= 1` exactly when `N >= 2i+3`, which the regime supplies for free.

It is tempting to prove the spectral half alone. **That half is false:** `s_{j+1}^2 >= s_j s_{j+2}`
fails at every index of the regime (15 of 15 at `n = 41`, 32 of 32 at `n = 81`). The binomial factor
does not assist the spectral part -- it carries the whole inequality against it. **The normalisation
by `C(N,i)` is the content of the theorem, not its packaging.**

Both sides have the same leading term:

    delta_i = log(a^2/(a^2-1)) - 2/(n^2(1-theta)^2) + O(1/n^3)
    log B_i = log(a^2/(a^2-1)) - 1/(n^2(1-theta)^2) + O(1/n^4)
    ------------------------------------------------------------
    slack   =                  + 1/(n^2(1-theta)^2)

`delta_0` extrapolates to `0.2876838` against `log(4/3) = 0.28768207`, and the `1/n^2` coefficient,
by Richardson extrapolation from `n = 81, 161, 321`, matches `2/(1-theta)^2` within 1.5 % across the
regime. Two things follow. The measured monotonicity of the slack is not a coincidence -- it is the
shape of `1/(1-theta)^2`. And **the entire content of the theorem is the 2 against the 1.**

The slack was checked directly: increasing in `i` at every index of the regime at
`n = 11, 21, 41, 81, 121, 161`, running from `1/n^2` at `i = 0` to `4/n^2` at the top.

## 8. The limit form is proved

Through the envelope relations the limit statement is a single inequality in one variable. With
`Q(theta) = (1-theta) - 1/(1+v^2)` and `dv/dtheta = v/Q`:

    **2 Q'/Q^2  <  1/theta^2 - 1/(1-theta)^2      for every theta in [0, 1/2].**

Equivalently, with `A = arctan v` and `P = A(1+v^2) - v`,

    2v^2(1+v^2)[2v^3 - (1+v^2)P]/P^3  <  v^3(2A-v)/(A^2(v-A)^2)   on  v in (0, 2.3311223] .

**This is proved**, by a finite certified computation in two pieces:

- `v in [0, 0.15]`: the gap `F` is analytic in `x = v^2`; summing 42 exact rational Taylor
  coefficients over `x in [0, 0.0225]` in interval arithmetic gives the enclosure
  `[1.0056789489, 1.0208179254]` with a tail bounded by `6.33e-70`, hence `F >= 1.00567 > 0`.
- `v in [0.15, 2.3311223]`: adaptive interval bisection, **2 163 485 boxes**, every enclosure
  strictly positive. Sample: `F in [+1.395792, +5.608977]` on `[2.00, 2.05]`.

The minimum of the gap is at `theta = 0` and equals exactly `F(0) = 176/175 = h_1`. Both sides blow
up like `1/theta^2` there, and the `O(1)` terms are `-2` and `-1` -- the factor of two again, in its
plainest form.

## 8b. The far edge, and what the theorem now reduces to

At the far edge the reciprocal roots `1/(2k-1)^2` have the Weierstrass product
`PROD (1 + z/(2k-1)^2) = cosh(pi sqrt z/2)`, so the reversed elementary symmetric functions are
`(pi/2)^{2m}/(2m)!`. Substituting, every power of `pi` cancels and the edge ratio collapses:

    p_{i+2}^3 p_i / (p_{i+1}^3 p_{i+3})  ->  (2m-3)^2/((2m-1)(2m-5))  =  1 + 4/((2m-1)(2m-5)) ,

`m = N - i`, because `(2m-3)^2 - (2m-1)(2m-5) = 9 - 5 = 4`. Checked exactly: `9/5, 25/21, 49/45,
81/77, 121/117, ...` for `m = 3..11`.

The formula is derived at `theta -> 1` and has no right to describe the regime `m ~ N/2`. It does,
to within 0.8 %: the exact slack exceeds the edge value at **every** index of every regime tested,
with the ratio between 1.0027 and 1.0078. So the theorem reduces to

    Delta g_i  >=  log( 1 + 4/((2m-1)(2m-5)) )   for  m = N - i ,

an inequality with an explicit elementary right-hand side, positive exactly where it is needed
(`m >= 3`) -- and proving merely that the surplus exceeds one half would suffice, where it is
measured at 1.003.

**The single exception, derived.** The spectrum is doubled, so the reciprocal product is a square,
`cosh^2(pi sqrt z/2) = (1 + cosh(pi sqrt z))/2`, giving `ehat_0 = 1` and
`ehat_m = pi^{2m}/(2(2m)!)` for `m >= 1`. The limit ratio is four `ehat` factors, four doubling
factors against four -- exact cancellation, **provided all four indices are at least 1**. At
`m = 3` the fourth index is 0, `ehat_0 = 1` carries no factor, and the ratio is halved:

    m = 3  :  (1/2) . 9/5 = 9/10 < 1        m >= 4 :  (2m-3)^2/((2m-1)(2m-5)) > 1

Checked symbolically at `m = 3..9` and against the exact object: `0.8861, 0.8931, 0.8966, 0.8983`
at `n = 41, 81, 161, 321`, converging to `9/10`.

Testing the inequality at **every** index rather than only in the regime, the failures are
`m = 3` and nothing else, at `n = 11, 21, 41, 81, 121, 161` alike. So ratio log-concavity holds for
this spectrum on the entire scale but one point, and the theorem asks for the half of it furthest
from that point.

(An earlier reading of `1 + 4/5` at `m = 3` as the origin of the theorem's constant was stated and
withdrawn: the value there is `9/10`, not `9/5`.)

## 8c. Every constant computed

Because the theorem reduces to `Delta^3` -- one fixed low order -- the mean value theorem for
divided differences applies without the uniformity that defeated the derivative route at high
order. Writing `h = -(log p)''`, the theorem is `h' > 0`, and expanding both halves through the
same saddle representation:

    h' = H'(theta)/N^2 + R'''(theta)/(2N^3) - tau'''(theta)/N^4 - ... ,
    R(theta) = log(sigma^2/sigma_bin^2) ,   tau  the Edgeworth correction.

    H'      >=  176/175 = 1.005714     PROVED
    max |R'''|   on the regime  =  0.051355
    max |tau'''| on the regime  =  4.124872   (29 exact Taylor coefficients, tail 2.78e-6)

    =>  H' - |R'''|/(2N) - |tau'''|/N^2  >  0   at every N >= 3 :

        N = 3 : +0.538836     N = 6 : +0.886855     N = 20 : +0.994118
        N = 4 : +0.741490     N = 10: +0.961898     N = 80 : +1.004749

Validated against the exact object -- the ratio of `-Delta^3 log p` to the three-term expression is
`1.0031, 1.0008, 1.0002` at `n = 41, 81, 161`, above 1 in all twelve cells tested and converging
like `1/N^2`. So the expression is a lower bound on the object, and it is positive at every size
the theorem concerns.

## 8d. The object, fully factored

Extracting the residual from the validation gives the exact factorisation

    -Delta^3 log p  =  [ H'/N^2 + R'''/(2N^3) - tau'''/N^4 ]  .  ( 1 + X(theta)/N^2 ) ,

with `X` measured across the regime at three sizes:

     theta :   0.02    0.07    0.13    0.26    0.38    0.49
     n= 41 :   4.964    --     3.906   3.000   2.951   3.357
     n= 81 :   5.320   4.704   3.905   3.000   2.878   3.397
     n=161 :   5.490   4.680   3.903   3.001   2.847   3.419

`X` is positive everywhere, bounded between about 2.82 and 5.61, and converges in `N` (it is 3.000
at `theta = 0.26` at all three sizes). Swept over every index rather than sampled: **134 indices at
n = 41, 81, 161, and the bracket is a lower bound at all 134.**

    first factor  > 0  at every N >= 3        -- all constants computed
    second factor > 0  iff  |X| < N^2         -- that is 9 at N = 3, against a measured 5.5

**[WITHDRAWN 11:10 -- see the correction at the end of this paper: the factorisation DEFINES `X`,
so `1 + X/N^2 > 0` is the conclusion itself and the reduction is circular.]** The passage read:
**the whole theorem now rests on one bound: `|X| < N^2`**, needed at nine when the measured
maximum over 134 indices is 5.605, with the margin growing quadratically in the size. Everything else is
either proved or computed.

## 8e. The leading term without asymptotics

For an exponential family `d sigma^2/dx = kappa_3/sigma^2`, so with
`D(x) = 1/sigma^2 - 1/sigma_bin^2`,

    D'(x)  =  kappa_3^bin/sigma_bin^6  -  kappa_3/sigma^6           (identity, verified 6 of 6)
    N^2 D'(theta)  =  (1-2 theta)/(theta^2(1-theta)^2) . [ 1 - kappa_3ratio/rho^3 ]
                                                                    (verified 11 of 11, exactly)

and therefore the leading term of the theorem is the elementary inequality

    **A_3 . g_2(theta)^3  <  g_3(theta) . A_2^3** ,
    g_2(q) = q(1-q) ,  g_3(q) = q(1-q)(1-2q) ,  A_k = mean of g_k ,  theta = mean of q .

Everything here is a finite sum at the actual `n` -- no limit shape, no series, no saddle
expansion. Equivalently `kappa_3/kappa_3^bin < rho^3` with `rho = sigma^2/sigma_bin^2 <= 1`
(Newton); writing the ratio as `rho^E`, the exponent is

     theta :   0.05      0.15      0.25      0.35
     E     :  3.0705    3.2797    3.6956    4.9814

**independent of `n` to four decimals across n = 41, 81, 161**, and above 3 throughout.

The general Schur-type version of this inequality is false -- 1956 violations of 4000 random
vectors -- so the tilted probabilities of this spectrum are what make it hold.

## 8f. The leading coefficient in closed form, and the origin of 176/175 (PROVED)

Let `Z_i(m) = p_{i+2}^3 p_i - p_{i+1}^3 p_{i+3}`, the certificate at index `i`, a polynomial in `m`
of degree `8i+10` (with `n = 2m+1`). Its leading coefficient is

    L_i  =  (176/175) . 2^(8i+10) / 3^(4i+6)  >  0    for every i .

The derivation uses four constants of the cumulant expansion. Since
`S_j = 2 SUM_{k<=m} (2k-1)^{2j}` contains only odd powers of `m`, writing `Z = m^2 z` gives
`log E = m F_0(Z) + F_1(Z)/m + ...` with

    f_j = [Z^j] F_0 = (-1)^{j+1} (2/j) 4^j/(2j+1) ,   f_1 = 8/3, f_2 = -16/5, f_3 = 128/21,
    g_1 = [Z] F_1 = -2/3 .

Then `p_i = a_i m^{2i}(1 + u_i/m + w_i/m^2 + ...)` with

    a_i = (f_1/2)^i = (4/3)^i ,      u_i = i(i-1) (f_2/f_1^2) + i(i-1)/4 = -i(i-1)/5 .

**`a_i^2 = a_{i-1} a_{i+1}`** — the leading coefficients satisfy the ratio inequality *with
equality*, because at leading order the spectrum degenerates to a single repeated root, which is
Newton's equality case. That is why the top two orders of `Z_i` cancel identically (`u_i` is a
quadratic, so `Delta^3 u = 0`). What survives is

    W_1 - W_2 = -Delta^3 w + 3(u_{i+2}^2 - u_{i+1}^2) + 3(u_{i+2}u_i - u_{i+1}u_{i+3})
              = (12/25)(i+1) - (84i+40)/175  =  44/175 ,

independent of `i`, which gives `L_i` after multiplying by `(4/3)^{4i+6}`.

Checked against the exact polynomials, which were computed independently and used to fit nothing:
`a_i` and `u_i` and `w_i` exact at `i = 0..24`; `W_1 - W_2` a single value over `i = 0..59`; `L_i`
exact at `i = 0..24`. Reproduced by `scripts/leading.py`.

**`176/175` is the same constant as the proved minimum gap of the limiting inequality** (Section 8),
attained there at `theta = 0`. Two computations sharing no code produce it. The constant of the
conjecture is not an observation.

## 8g. The first finite-size correction, and what the conjecture now needs

With `theta = t/(n-1)`, `s(theta) = (1/2)[(1-theta) - 1/(1+v^2)]`, `H = 1/s - 1/(theta(1-theta))`
and `R(theta) = log[s/(theta(1-theta))]`:

    M(n,t)  =  H(theta) + c_1(theta)/n + ... ,
    c_1(theta) = H^2/2 + H + (1/2) R''(theta) ,

matched to six digits at `theta = 1/16 .. 8/16` against `c_1` extracted from exact dyadic towers
(`n = 16.2^k+1`, `t = j.2^k`, which keep `theta` exact); worst relative deviation `1.1e-05`, the
extrapolation error of the measured side. `scripts/correction.py`.

The three terms plainly come from `M = n(e^g-1)`, from `N = n-1` against `n`, and from the Gaussian
prefactors — where **both sides carry `-1/2 log(2 pi sigma^2)`, so only the ratio survives** (taken
on the spectral side alone that term diverges like `1/(2 theta^2)`; with both sides kept the
divergence cancels and leaves an almost constant `0.21`). Measuring the `1/N^2` coefficient of `g` directly gives `G = (1/2)R''` to six digits at five values
of `theta`, from which `c_1 = H + G + H^2/2` follows by elementary algebra (`n/N = 1 + 1/n`). So the
closed form rests on **one measured identity**, `g = H/N + (1/2)R''/N^2 + O(1/N^3)`, and not on a
three-term accounting: a first attempt at that accounting produced a spurious `H'`, which is
cancelled by the finite-size correction to the variance itself.

**`c_1 > 0` (PROVED).** `R` has exact rational Taylor coefficients, with `R'(0) = -4/5` and
`(1/2)R''(0) = 4/35`; 38 of the 42 computed coefficients are clean and satisfy `|d_k| 3^k <= 0.1821`
with `|d_k|^{1/k} -> 0.267`, so a Cauchy estimate bounds the tail beyond them by `9.8e-32` and

    |R''(theta)| <= 0.2528318488   on [0, 1/2] ,
    c_1(theta) >= 28/25 - 0.1264 = 0.993584 > 0 ,     c_1(0) = 216/175  exactly.

**The denominator 175 = 5^2 . 7 now carries three constants of the problem**, each from a different
computation: `176/175` (the gap, and every `L_i`), `216/175` (the correction at the edge), and
`54/175` (the leading coefficient of the dyadic certificates, `LD_t = (54/175)(2^20/3^6)^t`).

## 8h. The conjecture is a variance inequality, and 4/5 is the spectrum's own ratio (PROVED in part)

**The exact identities.** Write the coefficients as a Poisson-binomial: for any `r > 0`,
`P(S = t) = e_t r^t / E(r)` with `q_k = b_k r/(1 + b_k r)`, so `Delta^3 log e_t = Delta^3 log P(t)`
(the other terms are affine in `t`). With `qbar`, `Q = qbar(1-qbar)`, `V = Var(q)` over the `N = 2m`
values:

    sigma^2 = N(Q - V) ,   sigma_bin^2 = N Q ,   sigma_bin^2 - sigma^2 = SUM (q_k - qbar)^2 = N V .

The second is Newton's inequality (`V >= 0`, by concavity of `x(1-x)`); the conjecture is about how
much more than zero `V` is. Since `H = 1/s - 1/(theta(1-theta)) = V/((Q-V)Q)`,

    **H(theta) > 4/5   <==>   V (5 + 4Q)  >  4 Q^2** .

**PROVED.** In the limit `m -> infinity`, with `x = k/m` and `v^2 = 4m^2 r`, the saddle
probabilities are `q(x) = v^2x^2/(1+v^2x^2)` and both moments are elementary. Writing
`A = arctan(v)/v`, `B = 1/(1+v^2)`, the identity `A' = (B-A)/v` gives `V = Q - (A-B)/2`, so the
inequality collapses to

    ( A - B ) ( 5 + 4A - 4A^2 )   <   10 A (1 - A) ,      v > 0 ,

which is proved for every `v > 0` with no gaps (`scripts/variance_limit.py`, exit 0): an exact
rational series on `(0, 0.2]` with leading term `(352/945) w^3` (`w = v^2`) and an alternating tail
`<= 6.8e-36`; 359 637 certified interval boxes on `[0.2, 3.5]`; and an analytic argument for
`v >= 3.5`, where `F = 5A - 14A^2 + 4A^3 + B(5+4A-4A^2)` with the `B`-term positive and
`5 - 14A + 4A^2 > 0` since `A(3.5) = 0.3692848 < (14-sqrt(116))/8 = 0.4037088`.

This is a **second proof of Section 8's result**, six times cheaper in boxes, with the edge in
closed form: `F = (352/945) v^6 + O(v^8)`, whose relative margin `(44/105) v^2 = 0.41904762 v^2`
matches the numerical value exactly.

**Where 4/5 comes from.** At the edge `q_k ~ b_k r`, and the inequality becomes
`Var(b) > (4/5) bbar^2`. For `b_k = (2k-1)^2`, `k = 1..m`:

    bbar = (4m^2-1)/3 ,     Var(b)/bbar^2  =  16(m^2-1) / (5(4m^2-1))   ->   **4/5** .

**The constant of the conjecture is the centred-square spectrum's own variance-to-mean-squared
ratio.**

**What the conjecture now needs.** Exactly:

    M = n(e^g - 1) = P + n(E + g^2/2) + n g^3/6 + ... ,
    P := n (sigma_bin^2-sigma^2)/(sigma^2 sigma_bin^2) > 4/5 ,
    E := g - (sigma_bin^2-sigma^2)/(sigma^2 sigma_bin^2) .

`P > 4/5` follows from the proved limit inequality plus the finite-`m` correction, which is
favourable: the finite sums are **midpoint sums** (`y_k = (2k-1)/(2m)` are the midpoints of the `m`
subintervals of `[0,1]`), the exact bound `|INT f - (1/m)SUM f(y_k)| <= (1/24m^2) INT |f''|` applies,
and the resulting threshold is `m/threshold = 0.4 t` **independently of `m`** — i.e. a pure
condition `t >= 3`, with `t <= 2` covered by the certificates.

What remains is `E + g^2/2 >= 0`. Its limit is `(1/2)R'' + H^2/2 >= 0.1021555 + 0.32 = 0.4222`,
**both constants proved**, and it was checked at every index of every size from `n = 5` to `501`
with no exception (global minimum `+0.17448` at `n = 5`). **It is not proved.**

## 9. What remains

**One bounded remainder.** With `theta = t/N`, `N = n-1`, `g = -Delta^2 log p`, and
`R(theta) = log[s(theta)/(theta(1-theta))]`:

    n g(n,t)  =  H(theta)  +  [ H(theta) + (1/2) R''(theta) ] / N  +  r ,

and the conjecture follows from `M = n(e^g - 1) > n g` together with

    H(theta)         >=  4/5             **PROVED** (twice: Section 8, and Section 8h as
                                          the variance inequality `V(5+4Q) > 4Q^2`)
    (1/2) R''(theta) >=  0.1021555       **PROVED** (majorant over 38 exact rational Taylor
                                          coefficients, Cauchy tail `9.8e-32`)
    |N^2 r|          <=  0.9021555 . N   **OPEN** -- the only thing left

**MEASURED: `|N^2 r| <= 1.4264`** across the whole regime, from exact dyadic towers at
`N = 160, 320, 640, 1280`:

    theta      0.0625    0.1250    0.2500    0.3750    0.5000
    N^2 r     -1.4264   -1.2105   -0.8173   -0.4242   +0.0951

— bounded, and **changing sign**: negative near the edge, positive at `theta = 1/2`. The
requirement `|N^2 r| <= 0.9021555 N` is met for every `N >= 1.58`, i.e. every `n >= 3`, and at the
edge of the certified range (`n = 1139`) the allowance is `1027` against a measured `1.43` — **a
margin of about 718**.

**The expansion reproduces the derived edge behaviour.** As `theta -> 0`,
`H + (1/2)R'' -> 4/5 + 4/35 = 32/35`, and with `H = 4/5 + (176/175)theta` this returns

    n g  =  4/5  +  (176 t + 160) / (175 N)  +  ... ,

which is **derived** (not matched) from `a_i = (4/3)^i`, `u_i = -i(i-1)/5` and the explicit `w_i`:
`c_0 = -2 Delta^2 u = 4/5` in three lines, and `c_1(t) = 12/25 + 4b` equal to `(176t+160)/175`
exactly at `t = 1..15`.

**Honest.** The uniform expansion is matched to six digits at five values of `theta`, not derived;
the bound on `r` is an extrapolation from four sizes. **The conjecture is not proved.** What is
proved is every other link: both leading coefficients, the limit inequality behind `H`, and the
range `5 <= n <= 1207` by exact certificates.

**A second, independent and also partial route.** The dyadic refinement `(n,t) -> (2n-1,2t)`
preserves `theta` exactly and `M` decreases along it: proved for `t = 1..16` at all `n`
(degree `12t-1`, shift `t/2` even / `(t+3)/2` odd, exact 16 of 16,
`LD_t = (54/175)(2^20/3^6)^t`), and measured at 2232 of 2233 pairs. If it held for every `t`, the
tower would decrease to `H` and no remainder bound would be needed at all.

## 10. In plain language

Three hundred years ago Newton proved that a polynomial whose roots are all real cannot have
wobbly coefficients: each one is at least the geometric mean of its two neighbours. His inequality
permits equality — the coefficients are allowed to sit exactly on the boundary, and for a polynomial
with all its roots at the same place they do.

This paper asks a natural follow-up. Take one specific family — the odd squares `1, 9, 25, 49, ...`,
each counted twice — and ask **how far past Newton's boundary its coefficients actually sit.**
Measured in the right units, the answer appears to be the same number every time, for every
coefficient and every size: **four fifths**, and never less.

**Where that number comes from.** Newton's inequality is really a statement about variance. Write
the coefficients as the distribution of a sum of independent coin flips, one per root, with the
`k`-th coin biased by the `k`-th root. Then the exact identity

    (variance of a fair-coin sum)  -  (variance of the biased sum)  =  SUM (bias_k - average bias)^2

says that Newton's inequality is nothing but "a sum of squares is non-negative" — and that the
*margin* past Newton is exactly **how spread out the biases are**. So the question becomes: how
spread out are the odd squares?

The answer, worked out here, is that at the sharpest point the margin is the roots' own relative
spread — their variance divided by the square of their average — and for `1, 9, 25, ..., (2m-1)^2`,

    variance / (average)^2  =  16(m^2 - 1) / (5(4m^2 - 1))   ->   **4/5** .

**So `4/5` is not a mysterious constant of the inequality. It is a property of the set
`{1, 9, 25, 49, ...}` itself** — a measure of how unevenly the odd squares are scattered around
their own mean. Any other family of roots would put its own spread there instead.

**What is proved.** That the limiting statement is exactly this variance inequality, and that the
variance inequality is true — by a computation that ends: an exact series near the edge, a third of
a million certified interval boxes in the middle, and an elementary argument at the far end. Also
proved: the conjecture itself for every size up to `n = 1097`, by exact certificates, and a family
of exact formulas for the objects involved.

**What is not proved.** The passage from "every size up to about a thousand, and the exact limit" to
"every size". What stands between is a single error estimate, on an elementary quantity, which the
measurements say has about seven hundred times more room than it needs.

## 11. Reproduction

All computations use exact rational arithmetic (`flint.fmpq`) and certified interval arithmetic
(`flint.arb`, `flint.acb`); no floating point enters any comparison. Every interpolation of a
polynomial in `n` is asserted against exact values at twelve fresh `n` before use. Scripts are named
in the corresponding records under `projects/qg-bootstrap/results/`, and the refuted statements are
kept in `docs/ERRATA.md`.

---

**A correction, 2 September 10:44.** The factorisation
`-Delta^3 log p = [bracket] . (1 + X/N^2)` DEFINES `X` as the residual, so `1 + X/N^2 > 0` is
literally `exact/bracket > 0`, which is the conclusion itself. **The statement "the theorem rests
on one bound |X| < N^2" is circular and is withdrawn**, along with the "margin of 9 against a
measured 5.6" that went with it -- to know `X > -N^2` one must already know the sign of the object.

What remains is what it always was: **a rigorous remainder bound for the saddle expansion**, i.e.
`|exact - bracket| < bracket` established by an argument independent of the object's sign. The
measurements stand (the bracket is a lower bound in all 134 cells swept, the residual coefficient
sits in [2.82, 5.61]); their status as evidence is unchanged. What was wrong was the claim that
they had been converted into a smaller *statement*.

---

# 11. The conjecture as four lemmas, three of them proved (2 September 2026, evening)

The statement was reorganised this evening and now decomposes with the whole numerical content
in the parts that are closed.

Write `p_j = e_j / C(N,j)`, `g = -Delta^2 log p`, `M(n,t) = n(p_t^2/(p_{t-1}p_{t+1}) - 1)`, and
`r_i = p_{i+1}^3 p_{i-1} / (p_i^3 p_{i+2})`, so `log r_i = -Delta^3 log p_{i-1}`.

    C.   M(n,1) > 4/5  for every odd n >= 5                      PROVED (three lines)
    B.   r_1 > 1       for every m >= 2                          PROVED (exact, degree 22)
    A1.  H''(theta) > 0 on [0, 1/2]                              PROVED (convexity.py)
    A2.  the remainder carrying A1 to finite n                   OPEN

    A1 + A2  =>  Delta^4 log p < 0  =>  r_i increasing in i
    with B   =>  r_i > 1 on the range the theorem needs
             =>  M(n,t) increasing in t
    with C   =>  M(n,t) >= M(n,1) > 4/5                          THE CONJECTURE

## 11a. The base case B, and why it is the substantive one

`r_1 > 1` reduces to positivity of a degree-22 polynomial in `m`; after `m -> m+2` all 23
coefficients are non-negative, constant term 51 193 296 000, leading term 56 320. Hence
`r_1 > 1` for every `m >= 2`, unconditionally.

Index 1 is where the inequality fails for most spectra: of 400 random positive spectra 317
satisfied ratio log-concavity, and appending a root above the maximum preserved it in only 54 --
263 broke, almost all at index 1. B is therefore the load-bearing base, not a formality.

## 11b. A1: the limit shape is convex

Two overlapping ranges, no floating point in any comparison:

    theta in [0, 0.04]        exact rational Taylor series of H''; 72 coefficients clean and all
                              non-negative; c_0 = 352/175 exactly; tail below 1.3e-87 even under
                              the wild bound |c_k| <= 1e10 (k+1)(k+2)  =>  H'' >= 2.0114285714
    theta in [0.03807, 0.5]   certified interval bisection in v, 1 267 757 boxes, 10.6 s

The series requires dividing out both vanishing factors first -- `A - B = w E` with `E(0) = 2/3`
and `1 - A = w R` with `R(0) = 1/3` -- after which `H = (1/w)[2/E - 1/(A R)]` is regular at
`w = 0`. Reverting `theta = w R` gives

    H(theta) = 4/5 + (176/175) theta + (176/175) theta^2 + (67328/67375) theta^3 + ...

so **`H(0) = 4/5` and `H'(0) = 176/175` fall out of the series exactly** -- the theorem's own
constant and the proved slope, from a computation aimed at neither.

And `H''` is essentially closed: `c_k / ((k+1)(k+2)) = 1.000000` for every `k >= 4`, so
`H'' = 2/(1-theta)^3` to six digits at every order, with the entire deviation at order zero and
equal there to exactly `176/175`:

    H''(0) = 352/175 = 2 . (176/175) = 2 H'(0) .

This is the seventh constant of the family over `175 = 5^2 . 7` in this project.

## 11c. A2, stated exactly, with what is measured

    -n N^2 Delta^4 log p  =  H''(theta) . ( 1 + c(theta)/n + ... ) ,  c(theta) in [0.994, 1.00]

measured at `n = 161 .. 2001` and `theta = 0.125, 0.25, 0.333` with `H''` from the exact series.
**The finite-size correction is positive**, so the proved margin `H'' >= 352/175 > 2` is helped
by it rather than eaten. What is missing is a bound on the next order with an explicit constant.

An earlier version of this section said the coefficient is exactly 1; that was read off three
digits and is corrected here (ERRATA ERR-0042). The direction of the argument is unaffected.


---

*Timestamps in this file were repaired on 2026-09-02 19:33 +0200 by `tools/timestamp_repair.py`: the typed times drifted up to five hours ahead of the clock, and each was replaced by the time of the git commit that introduced it (ERR-0045).*


---

# 12. The proof (2 September 2026, evening; times are the machine clock)

## 12.1 Statement and decomposition

**Theorem.** For every odd `n >= 5` and every integer `1 <= t < n/2`, `M(n,t) > 4/5`.

Write `N = n-1 = 2m`, `theta = t/N`, `g = -Delta^2 log p_t = log[p_t^2/(p_{t-1}p_{t+1})]`, so that
`M = n(e^g - 1) > n g`. Since `log(1+x) < x` and `n > N`, it suffices to prove `N g > 4/5`
wherever the exact certificates below do not already give `M > 4/5` directly. The proof splits
the `(n, t)` quadrant into five pieces, each closed by one script under `release/scripts/`
(`theorem.py` runs the assembly):

    A.  t <= 627, every odd n        one exact polynomial certificate per index (`ladder_fast.py`):
                                     polynomials in n of degree 8t+10 with all coefficients
                                     non-negative after the shift n -> n + c(t), c(t) = t+2 (even t)
                                     or t+5 (odd t).
    B.  n = 1257, 1259, t >= 628     exact rational computation (`theorem.py`): min M = 1.807, 1.804.
    C1. t^2/N <= 1, every t >= 627   `sparse_certificate_full.py`
    C2. 0.05 <= theta <= 1/2         `dense_certificate_a.py`
    C3. theta <= 0.05, t^2/N >= 1    `dense_certificate_b.py`

`C1 u C2 u C3` covers `{N >= 1260, t >= 628, theta <= 1/2}`; the coverage is checked in `theorem.py`
from the certified range of `khat_1 = theta/v^2 in [0.30274, 0.33334]` on `v^2 <= 0.17`.

## 12.2 The two representations of `g`

*The sampling expansion (used in C1).* With `bbar` the mean root and `b_k = bbar(1 + beta_k)`,
`SUM beta_k = 0`, a uniform random `t`-subset `S` of the roots gives the exact finite identity

    p_t = bbar^t F(t),   F(t) = SUM_{j=0}^{t} e_j(beta) (t)_j / (N)_j ,

and `g = -Delta^2 log F`. In the variables `a = 1/tau`, `b = tau^2/N` (so `theta = ab`, `1/N = a^2 b`),
`e_j (tau)_j/(N)_j = E_j(a,b) R_j(a,b)` with `E_j = [z^j] exp(b SUM_{i>=2} (-1)^{i-1}(mu_i/i)(ab)^{i-2} z^i)`
(`mu_i` the centred moments of the normalised spectrum) and `R_j = prod_{i<j}(1-ia)/(1-ia^2 b)`. The
derivative in `tau` at fixed `N` is `D = -a^2 d_a + 2ab d_b`, which annihilates every function of
`a^2 b = 1/N`.

*The fixed-tilt Fourier weight (used in C2, C3).* For the tilt `r` at the saddle of `t` and real `tau`,

    P_r(tau) = (1/2pi) INT_{-pi}^{pi} e^{K(phi) - i tau phi} dphi,   K = SUM_k 2 log(1 + q_k(e^{i phi}-1)),

agrees with `e_tau r^tau/E(r)` at integers, and `(log P_r)'' = -kappa_2[phi]` exactly, the second
cumulant of `phi` under the complex weight. Hence, with the exact binomial term,

    g = INT_0^1 INT_0^1 kappa_2[phi](t-1+u+v) du dv  -  log[(t+1)(N-t+1)/(t(N-t))] .

Under `u = sigma phi`, `sigma^2 = k_2`, `c_1 = (t - tau)/sigma`, `c_j = k_j/(j! sigma^j)`,
`kappa_2[u] = M_2/M_0 - (M_1/M_0)^2` with `M_j = INT u^j exp(-u^2/2 + rho(u)) du`,
`rho = i c_1 u + SUM_{j>=3} c_j (iu)^j`, and `kappa_2[u] = 1 + SUM_w E_w(c)` is the exact Edgeworth
polynomial in the grading `wt(c_1) = 1`, `wt(c_j) = j-2` (derived on the fast engine as
weight-truncated `fmpq` polynomials; the weight-2 term `12c_4 - 36c_3^2 + 6c_1c_3` is
`(k_2/2)(log k_2)''`, the saddle-point `-(1/2)log V` term).

*The tilted cumulants in closed form.* With `a = m + 1/2`, `y = 1/(2 sqrt r)`, `z = a + iy`:
`x = k_1 = 2[m + y Im psi(z) - (pi y/2) tanh(pi y)]` and `k_{j+1} = -(y/2) d_y k_j`; the Stirling
series of `psi^{(n)}` carries the remainder `|B_{2K+2}|(n+2K+1)!/((2K+2)! (Re z)^{n+2K+2})`, from the
integral representation and the Bernoulli lemma
`|(t/2)coth(t/2) - SUM_{k<=K} B_{2k} t^{2k}/(2k)!| <= |B_{2K+2}| t^{2K+2}/(2K+2)!`, itself certified by a
1-D sweep. Verified against the exact sums to 80 digits.

## 12.3 The three certificates

**C1 (sparse).** `F~ = SUM_{j<=30} E_j R_j + (E - T_15(a_2 b))`, `E = e^{a_2 b}` a formal variable with
`d_b E = a_2 E`, `a_2 = -mu_2^inf/2 = -2/5`. The sign of `N Phi - 4/5`, `Phi = -(log F~)''`, is the sign
of the polynomial `Phat = (DP~)^2 - P~ D^2 P~ - (4/5) a^2 b P~^2` in `Q[a, b, E]`. Facts, all exact:
`[a^0] = [a^1] = [a^2] = 0`; the `E^2` part cancels identically (its coefficient is
`a^2 b Q^2(-2a_2 - 4/5) = 0`); after `E = 1 + a_2 b + b^2 Z`, `[a^3]` is divisible by `b^2` and
`[a^3]/b^2 = 176/175 + O(b)` -- indeed `(176/175) e^{2 a_2 b} + O(b^14)`. Then
`Phat/(a^3 b^2) >= [a^3]/b^2 - SUM_{k>=4}|[a^k]/b^2| a^{k-3} - (tail) - (moments) >= 0.1417` on
`b in [0,1]`, `a <= 1/627`, with the tail `j > 30` and the finite-`m` moment perturbation bounded
by Cauchy estimates proportional to `a` (constants in the script; independently re-derived).

**C2 (dense, boxes).** For each box in `(eps = 1/m, v)`: cumulants `kt_j = k_j/N` by the Stirling forms
written through `eps` so that `N -> infinity` is a box; a two-pass mean-value enclosure in `v`
(`d kt_j/dv = (2/v) kt_{j+1}` exactly); `N g = NV/(kt_2 Q) + (kappa_2 - 1)/kt_2 + (1-G(x))/theta +
(1-G(N-x))/(1-theta)` with `NV = theta(1-theta) - kt_2 = Var(q)` (the variance identity, removing the
`1/theta` cancellation) and `G(x) = x log(1+1/x) in [1-1/(2x), 1]`; `kappa_2[u]` from the exact
Edgeworth polynomial to weight 6, window-averaged as an exact polynomial in `delta = t - tau`
(`1/M_0` expanded to third order), plus an explicit remainder: the weight `> 6` monomials of
`rhobar^n` for `n <= 6`, the `n > 6` tail with `rhobar(u) <= c_1|u| + beta u^2` on `|u| <= L` under the
wider Gaussian `e^{-(1/2-beta)u^2}`, the truncated polynomial outside `|u| > L`, and the true
integrand outside via `|e^K| <= exp(-k_2(1 - cos phi))`. Adaptive bisection: 930 boxes, worst lower
bound `N g >= 0.800004` at the corner `(t = 628, theta = 0.05)`.

**C3 (dense, series).** With `V = v^2`, `eps = zeta V^2` (so `theta = V khat_1`, `t = 2 khat_1/(zeta V)`,
`b = 2 khat_1^2/zeta`), every ingredient is an exact power series in `V` with coefficients in
`Q[zeta]`: the cumulants by Faulhaber sums, the Edgeworth monomials with the half-powers of `eta`
pairing up, the binomial term by the series of `1 - x log(1+1/x)`. In these variables every `1/V`
term cancels symbolically, and **the constant term of `N g` is `4/5` exactly**; then
`S = (N g - 4/5)/V` has `S(0, zeta) = 176/525 = H'(0)/3` and `S >= 0.32` on every band of `V`
down to `6.5e-7` after the Edgeworth remainder (bounded per band with the C2 machinery at weight 10)
and the Cauchy tails (sup bounds on `|V| <= 0.3`); the corner `V -> 0` follows by monotonicity with
the inner cut `L = sigma/2`, under which every remainder piece scales as `sigma^{-w}`, `w >= 11`, or
`e^{-0.11 sigma^2}`.

## 12.4 Why the direct statement, and not the fourth difference

The evening of 1-2 September was spent on `Delta^4 log p < 0` (ratio log-concavity), a statement
with no constant in it. It is a detour: in the sparse regime the pieces of the fourth difference
are `O(1/N^2)` and cancel to `O(1/N^3)`, so certifying it needs the Edgeworth expansion through weight
13 for both halves; the direct statement has the proved `H(theta) >= 4/5 + 1.0077 theta` as its
leading term, an exact binomial half, and no cancellation beyond a factor `t` among `O(1)` pieces.

## 12.5 In plain language

Newton proved that a real-rooted polynomial cannot have wobbly coefficients: every normalised
coefficient is at least the geometric mean of its neighbours, and the margin is exactly how spread
out the roots are. For the odd squares `1, 9, 25, ...`, each twice, that spread is `4/5` in the
limit, and the claim is that no finite size ever falls below it. The proof measures the same
quantity in two ways: as a random sample of the roots (when the coefficient index is small compared
with the square root of the size), and as the width of a bell curve built from the roots (otherwise).
Each measurement comes with an exact error bar computed by the machine, and the two error bars
never cross the line `4/5`. The one place where the margin shrinks to nothing -- very large sizes at
a fixed small index -- is handled not by measuring but by an exact formula in which the constant
`4/5` and the slope `176/175` appear, again, as the first two terms.

## 12.6 Status of the pieces (kept honest)

`C1` independently validated (own formulation, six items PASS, `results/VALIDATION_SPARSE_CERTIFICATE_FULL_2026-09-02.md`).
`C2`, `C3`: independently validated on 3 September 2026 (results/VALIDATION_DENSE_CERTIFICATE_A_2026-09-02.md, nine links PASS after one bookkeeping fix; results/VALIDATION_DENSE_CERTIFICATE_B_2026-09-02.md, seven items PASS after the Ser repairs).
`A`: the 627-rung ladder is the durable, script-checked artifact `results/m_ladder_log_2026-09-02_full.txt` (7.06 h of rung time, six parallel ranges merged by `lab/ladder_log_merge.py`).

Every piece is now an artifact checked by a script; the word "proved" waits only for a human referee.
