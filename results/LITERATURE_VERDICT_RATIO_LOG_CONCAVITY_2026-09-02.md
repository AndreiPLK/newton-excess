# Literature verdict: ratio log-concavity of the normalised coefficients of `prod_{k=1}^m (1 + (2k-1)^2 z)^2`

*Acquisition date: 2026-09-02, 19:32-19:50 (the end time as first typed ran a quarter hour ahead of the clock and is corrected here) (search performed by the literature agent; no computation heavier
than a 30-digit mpmath check of one identity). Every entry states what was actually read:
`FULL TEXT` (LaTeX source or PDF text extracted and quoted), `ABSTRACT ONLY -- theorem not verified`,
or `NOT ACCESSED` (metadata only). Nothing below is cited from memory without one of these labels.*

## 0. The statement searched for

Object: `N = 2m`, `b_k = (2k-1)^2` each doubled, `e_j` the elementary symmetric functions,
`p_j = e_j / C(N,j)`. Needed ("RLC", ratio log-concavity, first half):

    p_{i+1}^3 p_{i-1} > p_i^3 p_{i+2}      for 1 <= i <= N/2
    <=>  Delta^3 log p_i < 0
    <=>  the Newton ratio  R_i = p_i^2/(p_{i-1} p_{i+1})  is increasing in i on the first half.

Known in this repository before the search: RLC is FALSE for general real-rooted polynomials (97 of 400
random spectra), and it fails for this spectrum at `i = N-2` (ERR-0041), i.e. it is a first-half statement.

## 1. Verdict table

| # | lead | source (read status) | verdict |
|---|---|---|---|
| 0 | any inequality valid for ALL real-rooted polynomials | elementary argument, §2 | DOES NOT APPLY, structurally: such inequalities are invariant under coefficient reversal, RLC is reversal-odd |
| 1a | Newton / Maclaurin | Niculescu 2000 Thm 1.1 (FULL TEXT) | DOES NOT APPLY: second difference only |
| 1b | Rosset 1989 cubic ("higher-order Newton") | Rosset primary NOT ACCESSED (paywalled); statement as eq. (1.5) of Niculescu 2000 (FULL TEXT) | DOES NOT APPLY: reversal-symmetric discriminant condition, cannot imply a reversal-odd statement |
| 1c | Niculescu 2000 quartic / Sylvester discriminating families | Niculescu 2000 §1, §3, Prop. 5.1-5.3 (FULL TEXT) | DOES NOT APPLY: all are discriminants of sections, hence reversal-symmetric |
| 1d | Newton-type inequalities under extra hypotheses (roots in arithmetic / geometric progression, "spread" sets) | searches, §3.4 | NOT FOUND: no theorem giving Newton-ratio monotonicity under any root-configuration hypothesis |
| 1e | Guo 2020, arXiv:2012.03530 (two-sided bound on `D_{k+1}/D_k`) | FULL TEXT, Thm 1 quoted | DOES NOT APPLY: lower bound is below the RLC threshold by the factor `(1-sqrt(1-c_k))^2 < 1` at every `k`; the "simple sufficient condition" is a hypothesis, not a theorem about real-rooted polynomials |
| 2a | Arizmendi-Perales 2018, finite free cumulants | arXiv:1611.06598 (FULL TEXT source) | DOES NOT APPLY: no inequality among normalised coefficients; only conditional positive definiteness of cumulants (infinite divisibility) |
| 2b | Garza-Vargas-Srivastava 2026, finite free information inequalities | arXiv:2602.15822 (FULL TEXT source) | DOES NOT APPLY: theorems are about `Phi_n` (Coulomb potential of roots) and `chi_n` (log discriminant); the strings "log-concav" and "ratio" do not occur in the source. A tool summary claiming a "Theorem 8 on ratio log-concavity" was FALSE and is recorded as such |
| 2c | Martinez-Finkelshtein-Morales-Perales, finite free multiplicative convolution | arXiv:2309.10970, already in `results/FFP_LITERATURE_PASS.md` | unchanged: Schur-Szego composition preserves real-rootedness, gives nothing on third differences |
| 3a | Blackwell-Hodges 1959 lattice large-deviation local theorem | Ann. Math. Statist. 30, 1113-1120, DOI 10.1214/aoms/1177706094 (abstract + main formula visible on Project Euclid; FULL TEXT not extracted) | DOES NOT APPLY: i.i.d. only; error `O(n^-2)` relative with no explicit constant |
| 3b | Chaganty-Sethuraman 1985, arbitrary sequences, lattice case | Ann. Probab. 13, 97-114, DOI 10.1214/aop/1176993069 (ABSTRACT ONLY -- PDF is a scanned image with no text layer) | UNCLEAR / very probably DOES NOT APPLY: gives `P(T_n = k) ~` saddlepoint form with `o(1)` or `O(1/n)`-type relative error and no explicit constants (abstract wording only) |
| 3c | Deheuvels-Puri-Ralescu 1989, Edgeworth for non-identical Bernoulli sums | J. Multivariate Anal. 28, 282-303, DOI 10.1016/0047-259X(89)90111-5 (ABSTRACT ONLY) | DOES NOT APPLY as stated: expansion of the DISTRIBUTION FUNCTION to `O(sigma_n^-2)`, central regime, constants not explicit |
| 3d | Roos 2000 Krawtchouk expansion of Poisson-binomial around the binomial | Theory Probab. Appl. 45:2, 258-272 (2001), DOI 10.4213/tvp466 (ABSTRACT ONLY) | UNCLEAR: an exact finite expansion around the binomial with explicit bounds in total variation and point metric; bounds are on `P`, not on `Delta^3 log P`; the size parameter is `sum (p_k - pbar)^2`, which is `O(N)` here, not small. Worth one look only if the tilted comparison against the binomial is re-attacked |
| 3e | Anade-Gorce-Mary-Perlaza 2020 explicit saddlepoint error bound | arXiv:2007.05319 (ABSTRACT ONLY) | DOES NOT APPLY: CDF, i.i.d. |
| 3f | Daniels 1954; Jensen 1995 ch. 2-3, ch. 6; Butler 2007; Petrov 1975 ch. VII; Kolassa | NOT ACCESSED (books / paywalled) | UNCLEAR by construction: known to give relative error `1 + O(1/n)` with unspecified constants for lattice sums (Jensen ch. 6 is "Uniform saddlepoint approximations" per the OUP table of contents); no statement about differences of `log P` is expected from the genre |
| 3g | Hermite-Jensen limits via Petrov's method | arXiv:2511.02628 (FULL TEXT source), Thm `HermiteLimit` quoted | DOES NOT APPLY (already refuted in `results/THE_LITERATURE_ROUTE_JENSEN_HERMITE_PETROV.md`): central window `|m - mu| <= C sigma`, and Jensen hyperbolicity is strictly weaker than RLC |
| 3h | NEW, from DLMF: `E(r)` is a Gamma-function product; the tilted cumulants are polygamma functions of complex argument; DLMF 5.11.11 gives an explicit remainder bound | DLMF 5.4.4, 5.11.1, 5.11.10-11 (FULL TEXT of the DLMF section); identity checked numerically here to 1e-31 | APPLIES as a TOOL for items 2 and 3 of the open chain (remainder bounds), not as a theorem; see §4.3 |
| 4a | OEIS A008956 | FULL TEXT of entry | confirms the identification `e_j((2i+1)^2 : i = 0..n-1) = |4^k t(2n+1, 2n+1-2k)|`; no comment on log-concavity or ratios |
| 4b | Butzer-Schmidt-Stark-Vogt 1989 | Numer. Funct. Anal. Optim. 10, 419-488, DOI 10.1080/01630568908816313 (ABSTRACT ONLY) | UNCLEAR: "systematic treatment" of properties; no ratio result visible in the abstract; the paper is 70 pages and NOT ACCESSED |
| 4c | Chow 2024, "Central factorial numbers are Polya frequency sequences" | J. Math. Anal. Appl. 535(1) 128077, DOI 10.1016/j.jmaa.2023.128077 (ABSTRACT ONLY, second-hand via search summary) | DOES NOT APPLY: real-rootedness of the generating functions (odd and even indices), i.e. log-concavity level; nothing on third differences |
| 4d | Andrews-Gawronski-Littlejohn, Legendre-Stirling numbers | PDF (FULL TEXT extracted), Thms 5.5, 5.6, 5.10 | DOES NOT APPLY but IDENTIFIES the object: the single spectrum is an affine image of the Legendre-Stirling first-kind spectrum, `(2k-1)^2 = 4 k(k-1) + 1`; their first-kind results are the generating function, a recurrence and unimodality |
| 4e | Stirling numbers of the first kind: Hammersley 1951, Erdos 1953, Lieb 1968, Sibuya 1988 | NOT ACCESSED (metadata only) | NOT FOUND: no result on monotonicity of the Newton ratio of `s(n,k)` in `k`; the literature proves log-concavity (Newton), unique maximum and mode location |
| 4f | Heim-Neuhauser, CLT/LLT interpolating binomial and `s(n,k)` | arXiv:2208.09928 (ABSTRACT ONLY) | DOES NOT APPLY: central / local limit theorems |
| 4g | Shankar 2025, rows of super-recurrence triangles | arXiv:2508.12467 (ABSTRACT ONLY) | DOES NOT APPLY: ordinary log-concavity, gamma-nonnegativity, real-rootedness of rows; central factorial numbers not mentioned |
| 5a | Branden-Huh, Lorentzian polynomials | arXiv:1902.03719 (FULL TEXT source grepped) | DOES NOT APPLY: "Newton" occurs only in the bibliography; no quantitative gap in Newton's inequality |
| 5b | Strongly Rayleigh / ULC as used for Poisson-binomial | Tang-Tang survey arXiv:1908.10024 (FULL TEXT source) | DOES NOT APPLY: Newton's inequality stated as the necessary condition, no gap; the survey cites Rosset only as "higher order Newton's inequalities" |
| 5c | a general gap `p_j^2/(p_{j-1}p_{j+1}) >= 1 + c/N` | searches, §3.5 | NOT FOUND, and cannot exist without a spread hypothesis (equality iff all roots equal, Niculescu Thm 1.1) |
| 6a | Chen-Guo-Wang 2014, infinitely log-monotonic sequences | arXiv:1304.5160 (FULL TEXT source), Thm `rlcm` quoted | DOES NOT APPLY: requires `[log f]''` completely monotonic on `[1, inf)` for an `n`-indexed sequence `a_n = f(n)`, and produces sequences that are log-CONVEX with log-concave ratio; ours is a finite log-concave row |
| 6b | McNamara-Sagan 2010, infinite log-concavity | arXiv:0808.1065 (FULL TEXT source), §"real-rooted polynomials" read | DOES NOT APPLY: the `L`-operator tower (`b_k = a_k^2 - a_{k-1}a_{k+1}`), conjectures on `L` preserving real-rootedness (Stanley, Fisk); no statement on the ratio sequence |
| 7 | fermionic canonical partition functions (`Z_N = e_N(e^{-beta eps_k})`), convexity of the chemical potential in `N` | one search, §3.7 | NOT FOUND: no theorem beyond log-concavity located in one search; recorded so it is not searched twice without a better query |

**Net verdict.** No published theorem gives the third-difference statement for this spectrum, for the
Legendre-Stirling / central-factorial family, for Stirling numbers of the first kind, or for any class of
real-rooted polynomials. The literature stops at log-concavity (Newton), Turan-type and Jensen-hyperbolicity
statements, all of which this repository has already shown to be strictly weaker than RLC. The negative
answer is structural (§2), not an accident of coverage.

## 2. Why no general real-rooted theorem can work (elementary, not from the literature)

*Own observation, three lines, easily checked; not attributed to any source.*

Reversal `f(x) -> x^N f(1/x)` preserves real-rootedness (roots `b_k -> 1/b_k`, all nonzero here) and maps
`p_j -> p'_j = p_{N-j}`. RLC for `p'` at index `i'` reads `p_{N-i'-1}^3 p_{N-i'+1} > p_{N-i'}^3 p_{N-i'-2}`;
with `i = N-1-i'` this is exactly the strict NEGATION of RLC for `p` at index `i`. Hence:

* any inequality valid for every real-rooted polynomial, evaluated on consecutive normalised coefficients,
  is invariant under reversal, while RLC changes sign under reversal;
* at the middle indices `i in {N/2 - 1, N/2}` both `i` and `N-1-i` are in the theorem's range, so a general
  theorem giving RLC on `i <= N/2` would contradict itself.

Every discriminant-type inequality (Newton = quadratic sections, Rosset = cubic sections, Niculescu = quartic
sections, Sylvester's discriminating families) is of this reversal-invariant kind. This is the same fact as
the repository's "97 of 400 random spectra fail RLC", seen from the symmetry side, and it explains why the
observed failure sits at `i = N-2` (the reversal image of `i = 1`, where the margin is largest).

Useful reformulation that fell out (with `a = p_{k-1}, b = p_k, c = p_{k+1}, d = p_{k+2}`,
`D_k = b^2 - ac`, `D_{k+1} = c^2 - bd` the Newton defects):

    a c^3 - b^3 d  =  b^2 D_{k+1} - c^2 D_k ,    so    RLC at k  <=>  D_{k+1}/D_k  >  (p_{k+1}/p_k)^2 .

The threshold `(p_{k+1}/p_k)^2` lies strictly between the two thresholds `p_{k+2}/p_k` and `p_{k+1}/p_{k-1}`
of the "simple sufficient condition" in arXiv:2012.03530 (by Newton at `k+1` and at `k` respectively), which
is the precise form of the "misses by `alpha^2`" finding in `results/THE_PROBLEM_HAS_A_LITERATURE_AND_A_NAME.md`.

## 3. Lead-by-lead record

### 3.1 Lead 1 -- higher-order Newton inequalities

**Niculescu, C. P., "A new look at Newton's inequalities", J. Inequal. Pure Appl. Math. 1(2) (2000), Art. 17.**
FULL TEXT: PDF obtained from the EMIS Poland mirror
(`http://emis.icm.edu.pl/journals/JIPAM/images/014_99_JIPAM/014_99.pdf`, 186 931 bytes,
sha256 `6c6cdd0342459b7c471a81274b78bc264824fb916b44295a6ecc15cd60965c5a`; emis.de itself now redirects to zbMATH).
Local copy: scratchpad `niculescu2000.pdf` / `niculescu2000.txt`.

* Notation: `E_k = e_k / C(n,k)` -- exactly our `p_k`.
* Theorem 1.1 (Newton, Maclaurin), as written: "Let F be an n-tuple of non-negative numbers. Then:
  (1.1) `E_k^2(F) > E_{k-1}(F) E_{k+1}(F)`, `1 <= k <= n-1` unless all entries of F coincide;
  (1.2) `E_1 > E_2^{1/2} > ... > E_n^{1/n}` unless all entries of F coincide." The text adds that (1.1)
  holds for real, not necessarily positive, entries.
* Rosset's cubic Newton inequalities, eq. (1.5), as written:
  `6 E_k E_{k+1} E_{k+2} E_{k+3} + 3 E_{k+1}^2 E_{k+2}^2 >= 4 E_k E_{k+2}^3 + E_k^2 E_{k+3}^2 + 4 E_{k+1}^3 E_{k+3}`,
  `k = 0, ..., n-3`, equivalent to the non-negativity of the discriminant of
  `E_k x^3 - 3E_{k+1} x^2 y + 3E_{k+2} x y^2 - E_{k+3} y^3`; rewritten by Niculescu as
  `4 (E_{k+1}E_{k+3} - E_{k+2}^2)(E_k E_{k+2} - E_{k+1}^2) >= (E_{k+1}E_{k+2} - E_k E_{k+3})^2`.
  Hypothesis: all roots real. Niculescu attributes the inductive proof and the observation that (1.5) is
  strictly stronger than (1.1) to Rosset [20].
* The general order-`n` Newton inequalities (his `(N_n)`): non-negativity of the discriminant
  `D_n(1, -C(n,1)E_{k+1}/E_k, C(n,2)E_{k+2}/E_k, ..., (-1)^n E_{k+n}/E_k) >= 0`, `k = 0..m-n`.
* Theorem 2.1: for `alpha + beta = 1`, `j alpha + k beta in {0..n}`: `E_{j alpha + k beta} >= E_j^alpha E_k^beta`
  (log-concavity of `k -> E_k`, proved by Rosset's induction). Hypothesis: non-negative entries.
* Section 3: quartic Newton inequalities (Lemmas 3.1, 3.2 -- Euler's resolvent cubic). Section 5:
  Propositions 5.1-5.3 transfer the same inequalities to `e_k`, to `(-1)^k k! P^{(n-k)}`, and to
  falling-factorial coefficients (via Brenti, Thm 2.4.2).
* **Our parameters satisfy every hypothesis (real, positive roots).** What fails is the conclusion's shape:
  all statements are discriminants of sections, reversal-invariant (§2); none concerns the third difference
  of `log E_k` or the monotonicity of `E_k^2/(E_{k-1}E_{k+1})`. The paper contains no statement of that kind
  under any extra hypothesis on the roots.
* Later strengthening: none found that changes this.

**Rosset, S., "Normalized symmetric functions, Newton's inequalities and a new set of stronger inequalities",
Amer. Math. Monthly 96 (1989), no. 9, 815-819, DOI 10.1080/00029890.1989.11972286 (also 10.2307/2324844).**
NOT ACCESSED (paywalled; Crossref metadata verified). The statement is taken from Niculescu (1.5) above and
is cited as "higher order Newton's inequalities" in Tang-Tang arXiv:1908.10024 (line 677 of the source).
Verdict: DOES NOT APPLY, by §2.

**Sylvester-type cubic discriminant inequality.** This is (1.5) itself (Niculescu: "the equation
`E_k x^3 - 3E_{k+1}x^2 y + 3E_{k+2} x y^2 - E_{k+3} y^3 = 0` has all its roots `x/y` real if and only if
(1.5) holds"). Same verdict.

**Guo, "An inequality for coefficients of the real-rooted polynomials", arXiv:2012.03530.** FULL TEXT (LaTeX
source). Theorem `mainthm`, as written: "For a real-rooted polynomial `f(x) = sum_{k=0}^n C(n,k) a_k x^k`, if
`a_k a_{k+1}(a_k^2 - a_{k-1}a_{k+1}) != 0`, then
`(a_{k+1}^2/a_k^2)(1 - sqrt(1-c_k))^2 <= (a_{k+1}^2 - a_k a_{k+2})/(a_k^2 - a_{k-1}a_{k+1}) <= (a_{k+1}^2/a_k^2)(1 + sqrt(1-c_k))^2`
holds for `1 <= k <= n-2`, where `c_k = a_k a_{k+2}/a_{k+1}^2`." Hypotheses: real-rooted; the stated
non-vanishing. Ours satisfy both. By §2's reformulation RLC needs `D_{k+1}/D_k > (a_{k+1}/a_k)^2`; the
theorem's lower bound is that threshold times `(1 - sqrt(1-c_k))^2 < 1`, so the theorem never reaches RLC.
The abstract's "simple sufficient condition" `a_{k+2}/a_k <= D_{k+1}/D_k <= a_{k+1}/a_{k-1}` is a hypothesis
on the sequence (proved there for `C(n,k)` and for sequences with `a_k^2 >= 2 a_{k-1}a_{k+1}`), not a
consequence of real-rootedness. Verdict: DOES NOT APPLY.

**Extra-hypothesis versions (arithmetic / geometric progressions, spread sets).** Searched: "Newton's
inequalities" with "arithmetic progression", "geometric progression", "roots in an interval", "spread", "improved
constant"; Gaussian-binomial rows (roots in geometric progression) return only the `q`-log-concavity literature
(Krattenthaler 1989, Sagan 1992, Butler 1990 -- NOT ACCESSED, titles only), which is log-concavity in the
`q`-polynomial sense. Verdict: NOT FOUND.

### 3.2 Lead 2 -- finite free probability

**Arizmendi, Perales, "Cumulants for finite free convolution", J. Combin. Theory A 155 (2018) 244-266,
arXiv:1611.06598.** FULL TEXT (LaTeX source, file `finitecumulants28022017.tex`). The cumulants are defined
from the coefficients `a_k` of `p(x) = sum x^{d-k}(-1)^k a_k` by the moment-coefficient lemma
`a_n = (1/n!) sum_{pi in P(n)} d^{|pi|} mu(0,pi) m_pi` and the moment-cumulant theorem (their eq.
`cum-mom formula`). The only inequalities in the paper: Proposition `divinf` -- `p` is infinitely divisible
w.r.t. `boxplus_d` iff the sequence `{kappa_n}` (equivalently `{tilde kappa_n}`) is conditionally positive
definite -- and the corollary that the only infinitely divisible polynomials are the rescaled Hermite
polynomials. No Newton-type, log-concavity or ratio inequality. Verdict: DOES NOT APPLY.

**Garza-Vargas, Srivastava, "Finite free information inequalities", arXiv:2602.15822 (2026).** FULL TEXT
(LaTeX source `paper.tex`). Abstract, as written: "We develop finite free information theory for real-rooted
polynomials, establishing finite free analogues of entropy and Fisher information monotonicity, as well as
the Stam and entropy power inequalities. These results resolve conjectures by Shlyakhtenko and Gribinski ...".
Definition (their `def:finite_free_entropies`): `Phi_n(p) = (1/n) sum_i ((2/(n-1)) sum_{j != i} 1/(a_i - a_j))^2`
(a scaling of `sum_{i<j} (a_i - a_j)^{-2}`), and `chi_n[p] = (1/(2 C(n,2))) log Disc(p)`. Main results:
`Phi_{n-1}(tilde p') <= Phi_n(p)`; `1/Phi_n(p) + 1/Phi_n(q) <= 1/Phi_n(p boxplus_n q)`; entropy monotonicity
under differentiation; the finite free EPI. **The strings "log-concav", "ratio", "Newton", "ultra" do not occur
anywhere in the source.** A WebFetch summary of the PDF asserted a "Theorem 8: ratio log-concavity of
`e_k/C(n,k)`"; that assertion is FALSE (no such theorem exists in the source) and is recorded here so it is
not propagated. Verdict: DOES NOT APPLY. (Read beyond the question: the quantities are functions of the ROOTS,
`sum (a_i - a_j)^{-2}` and `prod (a_i - a_j)^2`, both explicit for our spectrum; they carry no information
about coefficient ratios.)

**Martinez-Finkelshtein, Morales, Perales, arXiv:2309.10970** -- already recorded in
`results/FFP_LITERATURE_PASS.md`; nothing new.

**Even hypergeometric polynomials and finite free commutators, arXiv:2502.00254.** Source grepped for
"central factorial", "(2k-1)": no occurrence. Our even polynomial `prod (x^2 - (2k-1)^2)` is not treated.
Verdict: DOES NOT APPLY.

### 3.3 Lead 3 -- local expansions with explicit constants for differences of `log P(S = j)`

Already recorded negative (do not repeat): Korolev-Zhukov 2000, Giuliano-Weber 2017,
Siripraparat-Neammanee 2021, J. Inequal. Appl. 2024 -- density bounds only.

**Blackwell, Hodges, "The probability in the extreme tail of a convolution", Ann. Math. Statist. 30 (1959)
1113-1120, DOI 10.1214/aoms/1177706094.** Abstract and main formula read on Project Euclid (the page shows
them in full). Hypotheses: i.i.d. integer-valued with g.c.d. of differences 1, m.g.f. finite near 0,
`E X_1 < a < sup X_1`. Result: `Pr{X_1 + ... + X_n = na} = pi_n^{**}[1 + O(n^{-2})]` with
`pi_n^{**} = ([m(a)]^n / (sigma sqrt(2 pi n))) [1 + (1/(8n))(mu_4/mu_2^2 - 3 - (5/3) mu_3/mu_2)]`
(moments of the tilted law). **What fails:** i.i.d. -- our summands are non-identical; and the `O(n^{-2})`
carries no explicit constant. It is the exact shape of expansion the project's "saddle expansion" uses,
so it is the right ancestor to cite, not a tool. Verdict: DOES NOT APPLY.

**Chaganty, Sethuraman, "Large deviation local limit theorems for arbitrary sequences of random variables",
Ann. Probab. 13 (1985) 97-114, DOI 10.1214/aop/1176993069.** ABSTRACT ONLY -- theorem not verified (the
Project Euclid PDF is a CCITT-fax scan with no text layer; local copy in the session tool-results folder).
Abstract, as written: "The results of W. Richter (Theory Probab. Appl. (1957) 2 206-219) on sums of
independent, identically distributed random variables are generalized to arbitrary sequences of random
variables `T_n`. Under simple conditions on the moment generating function of `T_n`, which imply that
`T_n/n` converges to zero, it is shown ... that `k_n(m_n)`, the probability density function of `T_n/n` at
`m_n`, is asymptotic to an expression involving the large deviation rate of `T_n/n`. Analogous results for
lattice valued random variables are also given." This covers non-identical Bernoulli sums (the m.g.f.
conditions are on `log E(r e^s)`), but the conclusion is "asymptotic to", i.e. relative error `o(1)`, with no
constant. Verdict: UNCLEAR on the exact error order (full text not read), DOES NOT APPLY as an explicit bound.
Follow-up by the same authors: "Strong large deviation and local limit theorems", Ann. Probab. 21 (1993),
DOI 10.1214/aop/1176989136 -- NOT ACCESSED.

**Deheuvels, Puri, Ralescu, "Asymptotic expansions for sums of nonidentically distributed Bernoulli random
variables", J. Multivariate Anal. 28 (1989) 282-303, DOI 10.1016/0047-259X(89)90111-5.** ABSTRACT ONLY
(ScienceDirect blocked; wording from the search index): expansion for the DISTRIBUTION of the sum of
independent zero-one variables with `sigma_n^2 -> infinity`, to order `O(sigma_n^{-2})`. Distribution
function, central regime, no explicit constants visible. Verdict: DOES NOT APPLY.

**Roos, B., "Binomial approximation to the Poisson binomial distribution: the Krawtchouk expansion",
Teor. Veroyatnost. i Primenen. 45:2 (2000) 328-344; Theory Probab. Appl. 45:2 (2001) 258-272,
DOI 10.4213/tvp466.** ABSTRACT ONLY, as written on mathnet.ru: "The Poisson binomial distribution is
approximated by a binomial distribution and also by finite signed measures resulting from the corresponding
Krawtchouk expansion. Bounds and asymptotic relations for the total variation distance and the point metric
are given." Relevance: the project's obstruction IS a comparison of the tilted Poisson-binomial with the
binomial at the same mean; Roos's expansion is exactly that comparison, exact and finite. What fails: the
bounds are on `P` (point metric), not on `Delta^3 log P`, and the expansion parameter is
`sum_k (q_k - qbar)^2`, which for our tilted vector is of order `N`, not small. Verdict: UNCLEAR, low priority.

**Anade, Gorce, Mary, Perlaza, "An upper bound on the error induced by saddlepoint approximations --
applications to information theory", arXiv:2007.05319.** ABSTRACT ONLY: bound on the CDF for i.i.d. sums
with finite third absolute moment. Verdict: DOES NOT APPLY (CDF, i.i.d.).

**Daniels 1954 (Ann. Math. Statist. 25, 631-650); Jensen, "Saddlepoint Approximations", OUP 1995 (ch. 2-3,
ch. 6 "Uniform saddlepoint approximations" per the publisher's table of contents); Butler, "Saddlepoint
Approximations with Applications", CUP 2007; Petrov, "Sums of Independent Random Variables", 1975, ch. VII;
Kolassa, "Series Approximation Methods in Statistics".** NOT ACCESSED. Recorded expectation, not verified:
these give lattice saddlepoint/Edgeworth expansions with relative error `O(1/n)` and unspecified constants;
none is known to state anything about differences of `log P`. Petrov ch. VII is the method named by
arXiv:2511.02628 (source line 718: "we employ Petrov's method (see Ch. VII of [Petrov75])").

**arXiv:2511.02628, "Hermite-Jensen limits and d-log-concavity of q-multinomials".** FULL TEXT (LaTeX
source). Theorem `HermiteLimit`, as written: "Fix `d >= 1` and `lambda in (0,1)`, and suppose `a, b -> +inf`
with `a/(a+b) -> lambda`. Then for every `C > 0`, uniformly for integers `m` with `|m - mu_{a,b}| <= C sigma_{a,b}`,
coefficientwise we have `J^{d,m}_{a,b}(X) = H_d(X) + O_{d,lambda,C}((a+b)^{-1/2})`." What fails for us: the
central window, and (already shown in this repository) Jensen hyperbolicity does not imply RLC. Verdict:
DOES NOT APPLY.

### 3.4 Lead 3, continued -- what the search found beyond the question (APPLIES as a tool)

The generating function of the single spectrum is a Gamma-function product. With `c = 1/(2 sqrt r)`,

    E_1(r) := prod_{k=1}^m (1 + (2k-1)^2 r)  =  (4r)^m |(1/2 + i c)_m|^2
                                            =  (4r)^m |Gamma(m + 1/2 + i c)|^2 / |Gamma(1/2 + i c)|^2 ,

because `1 + (2k-1)^2 r = 4r ((k - 1/2)^2 + c^2) = 4r |k - 1/2 + ic|^2`. Checked numerically here
(mpmath, 30 digits) at `m = 3, 8, 20` and `r = 0.013, 0.7, 5`: relative discrepancy `<= 3e-31` in all nine
cells (`scratchpad/gammacheck.py`). The doubled spectrum has `E(r) = E_1(r)^2`. Consistency with the known
limit: DLMF 5.4.4, as written, "`Gamma(1/2 + iy) Gamma(1/2 - iy) = |Gamma(1/2 + iy)|^2 = pi / cosh(pi y)`",
which is the `m -> infinity` form `prod (1 + (2k-1)^2 r) -> cosh(pi sqrt r / 2)` recorded in HANDOFF §7.1.
The same representation appears, in Beta-function form, in the OEIS A008956 comment by J. W. Meijer
(`(Pi/2) Beta(n-1/2-z/2, n-1/2+z/2)/Beta(n-1/2, n-1/2)` generating the rows), FULL TEXT of the entry.

Consequences for the open items 2 and 3 of the chain (remainder bounds), stated as leads, not results:

* the tilted cumulant generating function is `K(s) = log E(r e^s) - log E(r)`, so every cumulant `k_j` of
  the tilted vector is a finite combination of polygamma functions `psi^{(n)}` at the complex points
  `m + 1/2 + ic` and `1/2 + ic` with `c = 1/(2 sqrt(r))` -- a closed form for the midpoint sums of HANDOFF §4,
  with the Euler-Maclaurin boundary terms being the Stirling-series terms;
* DLMF 5.11.1 (Stirling's series for `Ln Gamma(z)`, sector `|ph z| <= pi - delta`), and the explicit bound
  DLMF 5.11.11, as written: "`|R_K(z)| <= ((1 + zeta(K)) Gamma(K) / (2 (2 pi)^{K+1} |z|^K)) (1 + min(sec(ph z), 2 K^{1/2}))`"
  for `|ph z| <= pi/2` in the expansion 5.11.10 `Gamma(z) = e^{-z} z^z (2 pi / z)^{1/2} (sum_{k=0}^{K-1} g_k/z^k + R_K(z))`
  (for `K = 1` the factor `1 + zeta(K)` is replaced by 4). DLMF 5.11(ii) also states that for the `psi`
  expansion 5.11.2 the remainder is bounded by `sec^{2n+1}(ph z / 2)` times the first neglected term (in the
  sector stated there); for real positive `z` the remainder is bounded by the first neglected term with the
  same sign. **DLMF 5.15 gives the polygamma expansions 5.15.8-5.15.9 only as Poincare expansions, with no
  error bound**, so bounds for `k_3, k_4` (which need `psi'', psi'''`) would have to be derived (Cauchy
  estimates from 5.11.11 on a disc, or the integral representations of DLMF §5.9 -- equation numbers not
  verified here).
* our argument has `ph z = arctan(c/(m + 1/2))`, well inside `|ph z| <= pi/2`, so 5.11.11 applies at the
  upper point for every `m`; at the lower point `1/2 + ic` the phase approaches `pi/2` as `c -> infinity`
  (i.e. `r -> 0`, the `t = 1` corner), which is exactly where the exact certificates already cover the
  statement.

This is the "classical object with classical identities" that HANDOFF §7.1 asked for. It is a tool for the
remainder bounds, not for item 1 (the cumulant inequality itself).

### 3.5 Lead 4 -- central factorial numbers and Stirling numbers of the first kind

**OEIS A008956** (FULL TEXT of the entry, retrieved 2026-09-02). Name, as written: "Triangle of central
factorial numbers `|4^k t(2n+1,2n+1-2k)|` read by rows (n>=0, k=0..n)." Comment (Sloane, 2011): "The n-th row
gives the coefficients in the expansion of `Product_{i=0..n-1} (x+(2i+1)^2)`, highest powers first". Recurrence
(Meijer): `t2(n,k) = (2n-1)^2 t2(n-1,k-1) + t2(n-1,k)`, `t2(n,0) = 1`, `t2(n,n) = ((2n-1)!!)^2`. References given:
Butzer et al. 1989; Riordan, Combinatorial Identities, Wiley 1968, p. 217. No comment on log-concavity,
real roots, ratios or asymptotics. So `e_j` of the single spectrum is `A008956(m, j)` exactly.

**Butzer, Schmidt, Stark, Vogt, "Central factorial numbers; their main properties and some applications",
Numer. Funct. Anal. Optim. 10 (1989) 419-488, DOI 10.1080/01630568908816313.** ABSTRACT ONLY (OpenAlex),
as written: "The purpose of this paper is to present a systematic treatment of central factorial numbers (cfn),
including their main properties, as well as to employ them in a variety of applications ...". 70 pages, NOT
ACCESSED. Verdict: UNCLEAR; if anyone obtains it, look for a section on asymptotics of `t(n,k)` in `k`.

**Chow, C.-O., "Central factorial numbers are Polya frequency sequences", J. Math. Anal. Appl. 535 (2024)
128077, DOI 10.1016/j.jmaa.2023.128077.** ABSTRACT ONLY, and second-hand (search-engine summary; ScienceDirect
and OpenAlex both withheld the abstract): the generating functions of the central factorial numbers of odd and
even indices are shown to be real-rooted, hence PF. For the first kind with fixed `n` this is immediate from the
product form; the content presumably concerns the second kind and/or the column direction. Nothing on third
differences is claimed in what was visible. Verdict: DOES NOT APPLY.

**Andrews, Gawronski, Littlejohn, "The Legendre-Stirling numbers" (Discrete Math. 311 (2011), per the search
index; the PDF read is the authors' preprint, sha256 `e601e6967777e53ef0bf0419051e603fe72f3a3aa6ed00ad2d42406da6eb9125`).**
FULL TEXT extracted. Theorem 5.5, as written: "The Legendre-Stirling numbers of the first kind `{Ps_n^{(j)}}`
satisfy the horizontal generating function (5.22) `<x>_n = sum_{j=0}^n Ps_n^{(j)} x^j`", with the worked example
"`<x>_4 = x(x-2)(x-6)(x-12) = -144x + 108x^2 - 20x^3 + x^4`" -- i.e. the roots are `k(k+1)`, `k = 0..n-1`.
Theorem 5.6: recurrence `Ps_n^{(j)} = Ps_{n-1}^{(j-1)} - n(n-1) Ps_{n-1}^{(j)}`. Theorem 5.10: "The unsigned
Legendre-Stirling numbers of the first kind ... are unimodal with either a peak or a plateau of 2 points."
No statement on ratios or higher-order log-concavity. **Identification (own algebra, one line):**
`(2k-1)^2 = 4k(k-1) + 1`, so `prod_{k=1}^m (1 + (2k-1)^2 z) = prod_{j=0}^{m-1} ((1+z) + 4 z j(j+1))`: the single
spectrum is the Legendre-Stirling first-kind spectrum `{j(j+1)}` under the affine map `c -> 4c + 1`. The
Legendre-Stirling / Jacobi-Stirling literature (Gelineau-Zeng arXiv:0905.2899; Mongelli, total positivity;
Andrews-Egge-Gawronski-Littlejohn arXiv:1112.6111 -- all NOT ACCESSED beyond titles/abstract snippets) is
therefore the correct neighbouring theory; what is known there is real-rootedness, total positivity of the
triangles, unimodality and (second kind) asymptotic normality (arXiv:1408.0477, title only). Verdict:
DOES NOT APPLY to RLC; APPLIES as the standard name of the object.

**Stirling numbers of the first kind.** Hammersley, "The sums of products of the natural numbers", Proc.
London Math. Soc. (3) 1 (1951) 435-452; Erdos, "On a conjecture of Hammersley", J. London Math. Soc. 28 (1953)
232-236 (DOI 10.1112/jlms/s1-28.2.232); Lieb, "Concavity properties and a generating function for Stirling
numbers", J. Combin. Theory 5 (1968); Sibuya, "Log-concavity of Stirling numbers and unimodality of Stirling
distributions", Ann. Inst. Statist. Math. 40 (1988) 693-714 (DOI 10.1007/BF00049427). All NOT ACCESSED
(metadata and search snippets only). What the snippets state: strong log-concavity in `k` (Newton, since the
row polynomial is `(x+1)...(x+n-1)`), unique maximum for `n >= 3`, log-concavity in three directions of the
triangle. **No source found stating that `s(n,k)^2/(s(n,k-1)s(n,k+1))` is monotone in `k`, or any third
difference of `log s(n,k)`.** The repository's own check (HANDOFF §3: `1, ..., m` satisfies RLC) is, as far as
this search can tell, unpublished. Verdict: NOT FOUND -- a genuine gap, not a wall.

**arXiv:2208.09928 (Heim-Neuhauser)** ABSTRACT ONLY: central and local limit theorems for a family interpolating
`C(n-1,k-1)` and `s(n,k)/n!`. DOES NOT APPLY. **arXiv:2508.12467 (Shankar)** ABSTRACT ONLY: log-concavity of
rows of super-recurrence triangles; central factorial numbers not mentioned. DOES NOT APPLY.

### 3.6 Lead 5 -- ULC, strongly Rayleigh, Lorentzian: a quantitative Newton gap

**Branden, Huh, "Lorentzian polynomials", arXiv:1902.03719.** FULL TEXT source grepped: the word "Newton"
occurs only in two bibliography entries ("On multivariate Newton-like inequalities", "M-matrices satisfy
Newton's inequalities"). No quantitative strengthening of Newton's inequality. DOES NOT APPLY.

**Tang, Tang, "The Poisson binomial distribution -- old & new", arXiv:1908.10024.** FULL TEXT source. Newton's
inequality stated (their eq. `Newton`) as the necessary condition for strong Rayleigh,
`a_i^2 >= a_{i-1}a_{i+1}(1 + 1/i)(1 + 1/(n-i))`, `1 <= i <= n-2`, "also said to be ultra-logconcave [Pem00]";
Darroch's rule quoted with the exact mode formula (Darroch 1964, Thm 4); Roos 2000 cited for the Krawtchouk
expansion; Rosset cited for "higher order Newton's inequalities". No gap, no third difference. DOES NOT APPLY.

**A gap `p_j^2/(p_{j-1}p_{j+1}) >= 1 + c/N` for a class.** Searches for "sharpened / refined / strengthened
Newton's inequalities", "Newton's inequalities variance of roots", "roots in an interval": nothing beyond
Kurtz's sufficient condition `4 a_{k-1}a_{k+1} < a_k^2` (a criterion for real-rootedness, the opposite direction)
and Niculescu's equality case. Since equality holds when all roots coincide (Niculescu Thm 1.1), any gap must
carry a spread hypothesis; none is in print that this search could find. Verdict: NOT FOUND. (Note the
project's own `M(n,1) > 4/5`, i.e. `R_1 >= 1 + 4/(5N)` for this spectrum, is such a gap and appears to be new.)

### 3.7 Lead 6 -- infinite log-concavity / log-monotonicity

**Chen, Guo, Wang, "Infinitely log-monotonic combinatorial sequences", Adv. Appl. Math. 52 (2014),
arXiv:1304.5160.** FULL TEXT source. Definition, as written: "A sequence `{a_n}_{n>=0}` is said to be ratio
log-concave if `{a_{n+1}/a_n}_{n>=0}` is log-concave." Theorem `rlcm`, as written: "Assume that `f(x)` is a
function such that `[log f(x)]''` is completely monotonic for `x >= 1`. Let `a_n = f(n)` for `n >= 1`. Then the
sequence `{a_n}_{n>=1}` is infinitely log-monotonic." Log-monotonic of order 2 is (their §1) "log-convex and
the ratio sequence log-concave". Theorem `rlcav`: ratio log-concavity plus an initial condition gives strict
log-concavity of `a_n^{1/n}`. **What fails:** (i) our sequence is a finite row indexed by `k` at fixed `N`;
(ii) it is log-concave (Newton), whereas the theorem manufactures log-convex sequences; (iii) the hypothesis is
on a generating function `f` with `[log f]''` completely monotonic, which no finite row supplies. Verdict:
DOES NOT APPLY. The `n`-indexed applications (Motzkin, derangement, Domb, ...) are irrelevant here.

**McNamara, Sagan, "Infinite log-concavity: developments and conjectures", Adv. Appl. Math. 44 (2010),
arXiv:0808.1065.** FULL TEXT source, section on real-rooted polynomials read. Tower: `L(a_k) = a_k^2 - a_{k-1}a_{k+1}`
iterated. Conjecture (Stanley; Fisk): "If `p[a_k]` has only real roots then the same is true of `p[L(a_k)]`";
Aissen-Schoenberg-Whitney quoted (Toeplitz total non-negativity). No statement about the ratio sequence
`a_{k+1}/a_k`. Verdict: DOES NOT APPLY (confirms HANDOFF §5: a different tower, agreeing at level 2 only).

### 3.8 Physics neighbour (one search only)

`e_N(e^{-beta eps_k})` is the canonical partition function of `N` free fermions; RLC in `N` is convexity of the
chemical potential `mu_N = -(1/beta) log(Z_{N+1}/Z_N)` in `N`. One search ("canonical ensemble fermions
partition function log-concave in particle number chemical potential convex elementary symmetric") returned
lecture notes and the recursion literature (Borrmann-Franke type), nothing on convexity of `mu_N`. NOT FOUND;
not pursued further. A better query would target "canonical ensemble" + "convexity of the free energy in N"
+ "finite systems".

## 4. What this search did NOT do (so the next one starts here)

* Did not open Rosset 1989, Butzer et al. 1989, Jensen 1995, Butler 2007, Petrov 1975, Kolassa, Riordan p. 217,
  Hammersley 1951, Erdos 1953, Lieb 1968, Sibuya 1988, Chaganty-Sethuraman 1985 full text, Deheuvels et al. 1989
  full text, Roos 2000 full text, Chow 2024 full text. All are either paywalled or books. None is expected,
  from its abstract, to contain a third-difference statement; Butzer et al. and Roos are the two worth a library
  request.
* Did not search MathSciNet/zbMATH reviews (no access from here).
* Did not run any computation beyond the nine-cell identity check.

## 5. Local copies (session scratchpad, not committed)

`C:\Users\user\AppData\Local\Temp\claude\C--Users-user-ScienceBro\ac66e2dc-eaec-44ab-b12d-0b61f841fa72\scratchpad\`:
`niculescu2000.pdf` (+ `.txt`), `agl_legendre_stirling.pdf` (+ `agl.txt`), `gammacheck.py`, and under `src/`
the arXiv LaTeX sources of 1611.06598, 2602.15822, 0808.1065, 1908.10024, 2511.02628, 2012.03530, 1604.05148,
1702.04665, 1010.2043, 1304.5160, 1902.03719, 2502.00254. Move to `research/literature/pdfs/` if they are to be
kept (the task forbade touching any file other than this one).
