# Quantitative Newton inequalities: a sharp constant for the centred-square spectrum, and Sibuya's 1988 conjecture

**Andrei Pluzhnik** (ORCID 0009-0005-5660-2603)

*Version 1.0.0, 3 September 2026. Archived at Zenodo, DOI [10.5281/zenodo.22282840](https://doi.org/10.5281/zenodo.22282840) (concept DOI 10.5281/zenodo.22282839). Every numbered claim below is produced by a script in `release/scripts/`; every
script writes a log under `results/`. Sections 6 and 7 state exactly what is proved and what is not.*

---

## Abstract

Newton's inequality (1707) says that the normalised coefficients `p_j = e_j/C(N,j)` of a real-rooted polynomial
are log-concave: `R_j := p_j^2/(p_{j-1}p_{j+1}) >= 1`. We prove two quantitative strengthenings.

**(1)** For the centred-square spectrum `{(n-2k)^2 : k = 1..n-1}` (the odd squares, each twice), which is the
spectrum appearing in the partial-wave positivity of deformed Veneziano amplitudes,

        M(n,t) := n (R_t - 1) > 4/5      for every odd n >= 5 and every t < n/2,

with the constant `4/5 = lim Var(b)/E(b)^2` sharp: `M(n,1)` decreases to it. The proof is a finite decomposition
into five machine-checked pieces and is complete.

**(2)** For the unsigned Stirling numbers of the first kind (`b_k = k`, `k = 1..N`), Sibuya conjectured in 1988
(Ann. Inst. Statist. Math. 40, eq. (3.4), unproved and, as far as our literature search reaches, untouched since)
that `R_j >= 1 + 1/(3n-j)` with equality at `j = 1`. We prove it for every `j <= 1000` and every `n`; for every
`j >= 1001` with `j/(n-1) <= 0.9`; and, at the top of the row, for every `n` and every `j` with at most 802
indices missing. The residual region is stated exactly in section 7.

Both proofs use the same instrument: the excess is written as an exact quantity with two representations -- a
sampling expansion of the elementary symmetric functions when the index is small, and the second cumulant of an
exactly tilted Fourier weight when it is not -- and every error term is an explicit interval computed in
certified arithmetic (`arb`), never an `O(.)`.

---

## 0. The result in plain words

Take any list of positive numbers and build from it the quantities that count how many ways you can pick one of
them, two of them, three of them, and so on, multiplying the picks together and adding up. Newton showed in 1707
that this sequence never has a dent in it: the middle of any three consecutive terms is at least as large as
what the two neighbours would predict. That is one of the oldest inequalities in algebra, and it is stated as
"at least" -- Newton's proof leaves open how much room is actually there.

**This paper measures the room.** For one particular list of numbers, the answer turns out to be a specific
number and not merely "some", and that number is `4/5`. It is the smallest possible: the room shrinks toward
`4/5` and never goes below it. Where the list comes from matters -- it is the list that decides whether certain
models of gravity and strings are physically consistent -- and a floor under it is a statement about how far
those models can be deformed before they break.

The second result concerns a guess made by Masaaki Sibuya in 1988 about a different, very classical list: the
whole numbers `1, 2, 3, ..., N`. He wrote down exactly how much room he believed there was, and, as far as we
can find, nobody proved or disproved it in the thirty-eight years since. We prove it over most of its range and
say precisely where our proof stops.

**How the room is measured.** Two pictures of the same quantity, used in different places. When few of the
numbers are involved, we expand it directly and control every leftover term by hand. When many are, we turn the
list into a probability distribution, tilt it until its centre sits exactly where we need it, and read the room
off the spread of the tilted distribution. Neither picture is an approximation: every leftover is a numerical
interval with proved endpoints, computed in arithmetic that carries its own error bars, so a claim is only made
when the interval sits strictly on the right side of zero.

**What a reader should be sceptical about.** The proofs are finite but large, and a machine checked the
arithmetic. Every step writes a log, every log is in this repository, and every script can be re-run. What we
ask a reader to verify is not the arithmetic but the setup: that the quantities certified are the ones the
theorems are about. Section 7 says what is not proved, and `LIMITATIONS.md` says what would change our minds.

---

## 1. The two statements

Let `b_1, ..., b_N > 0`, `e_j` their elementary symmetric functions, `p_j = e_j/C(N,j)`. Newton: `R_j >= 1`.
The **Newton excess** `M = n(R_j - 1)`, `n = N+1`, measures the distance past that boundary. Two elementary facts
fix the scale:

- at the first index the excess is exactly the relative variance of the roots,
  `M(n,1) = rho n/(N - 1 - rho)`, `rho = Var(b)/E(b)^2` (Newton's identities);
- for `b_k = k^p` on `k = 1..m`, `rho -> Var(u^p)/E(u^p)^2 = p^2/(2p+1)` -- so `1/3` for `p = 1`, `4/5` for
  `p = 2`, `9/7` for `p = 3`.

**Theorem A.** For the centred-square spectrum, `M(n,t) > 4/5` for every odd `n >= 5` and every `t < n/2`.

**Theorem B (Sibuya's conjecture, partial).** For `b_k = k`, `k = 1..n-1`,
`R_j >= 1 + 1/(3n-j)` in the ranges listed in the abstract; equality holds at `j = 1` for every `n`.

The universal form of these statements is **false**: for the density `1/(1+u)` the minimum of `M(n,t)` over `t`
falls below `rho` (measured, `results/FAMILY_POWER_SPECTRA_2026-09-02.md`). Both theorems are statements about a
family, not about all real-rooted polynomials.

## 2. Where the centred-square spectrum comes from

It is the spectrum of the graviton partial-wave problem for the three-parameter Veneziano deformations of
Cheung, Hillman and Remmen (arXiv:2406.02665). The boundary of the positivity island of that family -- the
region where the deformation is still unitary -- is governed by the same odd squares; the exact edge laws are in
the companion draft `article/qg-island-draft.md`. Theorem A is therefore a statement about how sharply the
coefficients of that boundary problem are log-concave. Applying it back to the island is not done here.

## 3. The instrument

Write `g = -Delta^2 log p_t`, so `R_t - 1 = e^g - 1 > g`. Two exact representations of `g`:

**(a) Sampling.** With `b_k = bbar(1+beta_k)`, a uniform random `t`-subset of the roots gives the finite identity
`p_t = bbar^t SUM_j e_j(beta) (t)_j/(N)_j`. In the variables `a = 1/tau`, `b = tau^2/N` the derivative in `tau`
at fixed `N` is `D = -a^2 d_a + 2ab d_b`, which annihilates every function of `1/N`; the sign of `N g - c` is the
sign of a polynomial in `Q[a, b, E]` with `E = e^{a_2 b}` carried as a formal variable (`d_b E = a_2 E`). Its
first three coefficients in `a` vanish identically and the fourth is the constant of the theorem.

**(b) Tilting.** For the tilt `r` at the saddle of `t`, `P_r(tau) = (1/2pi) INT e^{K(phi) - i tau phi} dphi`
agrees with `e_tau r^tau/E(r)` at integers and satisfies `(log P_r)'' = -kappa_2[phi]` **exactly**. Hence
`g = INT_0^1 INT_0^1 kappa_2[phi](t-1+u+v) du dv - log[(t+1)(N-t+1)/(t(N-t))]`, the second term exact.
`kappa_2` is expanded by an exact Edgeworth polynomial in the scaled cumulants (derived on weight-truncated
`fmpq` dictionaries), and the tilted cumulants themselves are closed forms: polygamma at complex argument for
the centred squares, `k_1 = N - x[psi(N+1+x) - psi(1+x)]` with `k_{j+1} = v d_v k_j` for `{1..N}`.

Three devices make the intervals survive: the **variance identity** `NV = theta(1-theta) - k_2/N = Var(q)`,
which removes the `1/theta` cancellation; a **two-pass mean-value enclosure** in the tilt, using the exact
derivative rule for the cumulants; and, near the degenerate corner, an **exact power series** in the tilt whose
coefficients live in `Q[zeta]` and in which every `1/V` term cancels symbolically.

## 4. Theorem A: the five pieces

| piece | region | script | result |
|---|---|---|---|
| A | `t <= 627`, every odd `n` | ladder | 627 exact polynomial certificates, all coefficients non-negative after the shift `n -> n + c(t)`; log `results/m_ladder_log_2026-09-02_full.txt` (7.06 h), merged and checked by `lab/ladder_log_merge.py` |
| B | `n = 1257, 1259`, `t >= 628` | `theorem.py` | exact rationals, `min M = 1.807, 1.804` |
| C1 | `t^2/N <= 1`, `t >= 627` | `sparse_certificate_full.py` | `[a^3]/b^2 = 176/175` exactly, residual `0.1417` |
| C2 | `0.05 <= theta <= 1/2` | `dense_certificate_a.py` | 930 boxes, worst `N g >= 0.800004` |
| C3 | `theta <= 0.05`, `t^2/N >= 1` | `dense_certificate_b.py` | series in `V`: constant term `4/5` exactly, `S(0) = 176/525`, 18 bands, `min >= 0.3206` |

`theorem.py` re-runs the ladder check and the coverage lemma (1 s) and, with `--full`, the three certificates
(409 s, three `VERDICT True`, exit 0). Each certificate has an independent validation report under
`results/VALIDATION_*.md`, written by an agent that did not import the code it validated; three defects found
there (a missing term in a `1/M_0` bound, a non-analytic base object, a non-sigma-free tail) were fixed and the
scripts re-run.

## 5. Theorem B: the same machine on `{1..N}`

| piece | region | script | result |
|---|---|---|---|
| A' | `j <= 1000`, every `n` | `lab/sibuya_ladder.py` | 1000 exact certificates, degree `2j-1`, shift `j+2`; `j = 1` is an identity; log `results/sibuya_ladder_log_full.txt` |
| C1' | `j^2/N <= 1`, `j >= 1001` | `sibuya_sparse_certificate.py` | `[a^3]/b^2 = 2/3` exactly, residual `0.29` |
| C2' | `0.1 <= theta <= 0.9` | `sibuya_dense_a.py` (`--top` for the upper half) | 756 + 84 boxes, margins `5e-6` and `0.163` (logged) |
| C3' | `theta <= 0.1` | `sibuya_dense_b.py` | series in `v`: `S(0, zeta) = 1/3` exactly, 18 bands, `min >= 0.2398` |
| D' | `N - j <= 32`, every `N` | `sibuya_harmonic.py` | exact recursion to `N = 10^6`, bands in `H`, analytic tail at `2^48` |
| D'' | `33 <= N - j <= 802`, every `N` | `sibuya_corner_grid.py` | exact recursion to `10^6`; a third-order grid in `N` on exact Gamma-ratio series to `e^20`; beyond that the `H`-model with an explicit `O(e^{-H})` defect, certified for all 770 indices by ONE Taylor shift each (177 s), the grid taking 4705 steps to `N = 4.85e8` with worst relative margin `1.92e-10` |
| E' (exploratory, not a certificate) | `0.9 <= theta <= 0.99`, `N - j >= 100` | `sibuya_dense_a.py --top2` | 303 boxes, margin `1.03` (exploratory sweep: a sliver along the validity curve is skipped) |

**Why the `H`-model half needs no certificates of its own.** The limit of `E_k(N) = e_k(1, 1/2, ..., 1/N)` at
fixed harmonic number is `Ehat_k(H) = [x^k] Phi(x)` with

    Phi(x) = e^{Hx} PROD_{n>=1} (1 + x/n) e^{-x/n} = e^{(H - gamma) x} / Gamma(1 + x) ,

a Laguerre-Polya function: genus one, no Gaussian factor, zeros exactly at `x = -1, -2, -3, ...`, all real,
negative and simple. Such functions are locally uniform limits of real-rooted polynomials, so Newton's
inequality passes to their coefficients by Hurwitz, in the form `a_k^2 >= ((k+1)/k) a_{k-1} a_{k+1}` that the
finite Newton factor `T_inf = (j'+1)/j'` tends to. That is exactly the statement the shift certificates verify
one index at a time -- they were re-deriving a classical inequality. The mathematics that remains is the
finite-`N` part: Sibuya's inequality is this limit statement TIMES the extra factor `1 + 1/(2N+3+j')`, and it is
that factor, and only that factor, which the machine work has to reach.

**A repair the debate forced.** The `O(e^{-H})` defect above was, until 3 September, bounded through a ratio
`Ehat_{k-2}/Ehat_k <= k^2/H^2` -- the Poisson estimate, applied where `k/H` reaches 40. Measured, that ratio
reaches 1.70 at `H = 20` and 2.79 at `H = 12`: the step is false, and the bound it supported survived only on
unstated slack. The same derivation substituted `1/N <= e^{gamma-H}`, whose true direction is the reverse. The
defect is now carried as a POLYNOMIAL in `H`, from one convolution shared by every index, with no ratio and no
asymptotic form anywhere; it is about 2.8 times tighter than what it replaces. Details and the measurements are
in `results/DEBATE_REPAIR_TAIL_2026-09-03.md`.

`sibuya_theorem.py` assembles A', C1', C2', C3' and checks their coverage; its `--full` run is logged
(`results/sibuya_theorem_full_2026-09-03.txt`, four certificates, exit 0). The sparse and both dense
certificates have independent validation reports (`results/VALIDATION_SIBUYA_*.md`, 8/8, 5/5 and 7/7 items).

## 6. What is proved

- **Theorem A: complete.** Every piece is an artifact checked by a script, and every certificate is
  independently validated (`results/VALIDATION_SPARSE_CERTIFICATE_FULL_2026-09-02.md`,
  `VALIDATION_DENSE_CERTIFICATE_A_2026-09-02.md`, `VALIDATION_DENSE_CERTIFICATE_B_2026-09-02.md` -- the last one
  after two addenda, both recorded in the file). The remaining step is human refereeing.
- **Theorem B: the ranges of section 5.** Together they give Sibuya's inequality for every `n` and every
  `j <= 1000`; for every `j >= 1001` with `theta <= 0.9`; and, at the top of the row, for every `n` and every `j`
  with at most 802 indices missing. The sweep E' is exploratory and is not counted in that list.

## 7. What is not proved

Sibuya's inequality in the top decile of the row when many indices are missing:

    j >= 1001,   theta = j/(n-1) > 0.9,   and   n - 1 - j >= 803 .

Everything else is covered: the ladder takes `j <= 1000` at every `n`, the three certificates take `theta <= 0.9`
at every `j >= 1001`, and the harmonic instruments take the top of the row whenever at most 802 indices are
missing, at every `n`.  (The exploratory sweep E' of section 5 also covers part of the residual, but it skipped
473 of its 776 boxes and is not counted as a certificate anywhere in this article.)

Two instruments meet at that boundary and neither crosses it, for reasons that were measured, not guessed:

1. The **harmonic side** (exact `e_k(1, 1/2, ..., 1/N)` and its `H`-model) is an alternating sum whose
   cancellation grows with the index; the `H`-model's Taylor-shift certificate stops working around 1000 missing
   indices at any starting height, because `Fhat` is not yet monotone for `H` far below the index.
2. The **box side** (the tilted Fourier weight) needs `kt_2 = j'/N` to stay above the Edgeworth enclosure; on the
   sparse side `N > j'^2` it collapses, and the region has a cusp at `(1/N, 1/v) -> (0,0)` that a box scheme in
   those variables walks down without ever becoming relatively thin.  The instrument written for it,
   `release/scripts/sibuya_top_w.py`, certifies boxes away from that cusp and **exits 1 at it**; it is shipped as
   the starting point for closing the region, not as a certificate.

What would close it: the same machine in the variables `(mu, z)` with `mu = j'/N` and `1/N = mu^2 z`, in which the
sparse side is the rectangle `0 <= z <= 1` instead of a cusp; or a lower bound for the Newton gap of the
reciprocal spectrum `{1/k}` that is uniform in the index.

## 8. Reproducing

```
uv sync
uv run python projects/qg-bootstrap/release/scripts/theorem.py --full          # Theorem A, ~7 min
uv run python projects/qg-bootstrap/release/scripts/sibuya_theorem.py --full   # Theorem B, ~14 min
uv run python projects/qg-bootstrap/release/scripts/sibuya_harmonic.py 33      # the top corner, ~5 min
```

Everything that decides a claim is exact or certified interval arithmetic (`python-flint`: `fmpq`, `fmpq_poly`,
`fmpz_poly`, `arb`, `acb`). Floating point appears only in printing, in the adaptive choice of box and step
sizes, and in loop bookkeeping -- never in a comparison that a conclusion rests on. Logs of every run quoted above are in `results/`.

## 9. Sources

Newton, *Arithmetica Universalis* (1707). Sibuya, *Log-concavity of Stirling numbers and unimodality of Stirling
distributions*, Ann. Inst. Statist. Math. 40 (1988) 693-714, eq. (3.4) (the conjecture; full text read).
Lieb, J. Combin. Theory 5 (1968) 203-206 (strict Newton for both kinds; abstract only). Hwang, J. Combin.
Theory Ser. A 71 (1995) (uniform asymptotics for `s(n,m)`, `O`-terms only; full text read). Cheung, Hillman,
Remmen, arXiv:2406.02665 (the Veneziano deformations). Literature verdicts with exact locations:
`results/LITERATURE_VERDICT_STIRLING_NEWTON_GAP_2026-09-02.md`,
`results/LITERATURE_VERDICT_STIRLING_SMALL_K_2026-09-03.md`.

## 10. AI disclosure

The proofs were constructed and the code written in collaboration with an AI assistant (Claude), under the
author's direction; every certificate is machine-checked and every claim is traceable to a logged run. See
`release/AI_DISCLOSURE.md`.
