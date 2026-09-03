# Literature verdict: quantitative Newton gap for Stirling numbers of the first kind

Written: Wed Sep  2 23:30:18 RDT 2026 (output of `date` at time of writing)
Branch: claude/handoff-markdown-review-anzwo5
Author: literature search agent for the qg-bootstrap Newton-excess programme

**Question (from the structural agent, 2026-09-02 evening).** For the unsigned Stirling numbers of the
first kind c(n,m) (row = elementary symmetric functions of 1, ..., n-1; real-rooted), Newton's inequality
gives log-concavity of the normalised coefficients p_j = e_j / C(N, j). Is there ANY published quantitative
strengthening of the form p_j^2/(p_{j-1} p_{j+1}) >= 1 + c/N (equivalently N (ratio - 1) >= const), for
Stirling numbers of either kind, binomial-type triangles, or Eulerian numbers? Does Var/mean^2 of the
limiting root density appear anywhere as a Newton-gap constant?

**Our finding to be checked against the record.** For the spectrum {1, ..., m}, min over j of
(m+1)(p_j^2/(p_{j-1}p_{j+1}) - 1) is attained at j = 1 and decreases to exactly 1/3 as m grows; for {k^2}
the limit is 4/5, for {k^3} it is 9/7; in general q^2/(2q+1) = Var(u^q)/E(u^q)^2 on [0, 1].

**Notation dictionary used throughout (needed because no source uses ours).**
Sibuya writes [n m] for c(n,m) and {n m} for S(n,m), with x(x+1)...(x+n-1) = sum_m c(n,m) x^m.
Hence c(n,m) = e_{n-m}(1, ..., n-1): the row has N = n-1 roots, and OUR spectrum {1, ..., m} is
Sibuya's n = m+1. Our p_j = e_j(1..N)/C(N,j) = c(n, n-j)/C(n-1, j), j = n-m. Conversion:
p_j^2/(p_{j-1}p_{j+1}) = [c(n,m)^2/(c(n,m-1)c(n,m+1))] * (m-1)(n-m)/(m(n-m+1)) with j = n-m.
Newton for the degree-(n-1) polynomial (x+1)...(x+n-1) reads, in Sibuya's indices,
c(n,m)^2/(c(n,m-1)c(n,m+1)) >= m(n-m+1)/((m-1)(n-m)).

---

## VERDICT IN FIVE LINES

1. **The j = 1 value is a classical identity, not a discovery**: p_1^2 - p_2 = sigma^2/(N-1) (Newton's
   identities), so p_1^2/p_2 - 1 = sigma^2 / ((N-1) mu^2 - sigma^2); for {1..N} this is exactly 1/(3N+2),
   which gives (N+1)(ratio - 1) -> 1/3, and for {k^q} it gives Var(u^q)/E(u^q)^2 in the limit.
2. **The full statement (minimum over j at j = 1, constant 1/3) was conjectured in print by Sibuya (1988),
   eq. (3.4), as "numerically suggested", in the slightly stronger j-dependent form**
   p_j^2/(p_{j-1}p_{j+1}) >= 1 + 1/(3n - j), with EQUALITY at j = 1 for every n (translation and equality
   verified exactly here for n = 4..12). Sibuya did not prove it.
3. **No later paper proves, restates, or even mentions Sibuya's (3.4).** All 33 works citing Sibuya in
   OpenAlex were listed; the ones that use his inequalities (Hong-Zhang 2020, Pitman 1997, Qi 2014/2016)
   use the proved ones ((3.1), (3.6), second-kind results). The conjecture is open as far as this search
   can determine.
4. **For Stirling numbers of the SECOND kind a strengthening IS proved** (Sibuya 1988, Theorem 3.1 second
   half, eq. (3.2)): the Newton ratio exceeds Newton's bound by the factor (m+1)/m, i.e. in normalised
   form p_j^2/(p_{j-1}p_{j+1}) > 1 + 1/(j+1). That is a 1/j gap, not a 1/N gap: a different regime.
5. **Nothing of the form 1 + c/N with c = Var/mean^2 of the root density exists in the literature for
   j >= 2**, in Stirling, Eulerian, binomial-type, sharpened-Newton, finite-free-cumulant, or
   exponential-profile sources. Finite free cumulants reproduce only the j = 1 identity.

Consequence for the programme: a PROOF of the uniform-in-j gap (either Sibuya's 1 + 1/(3n-j) or our weaker
"min at j = 1") would be new; the CONJECTURE has priority Sibuya 1988 and must be cited as such.

---

## CANDIDATE RESULTS, ONE BY ONE

### 1. Sibuya, M. (1988). Log-concavity of Stirling numbers and unimodality of Stirling distributions.
Ann. Inst. Statist. Math. 40(4), 693-714. DOI 10.1007/BF00049427. OpenAlex W2094158818.
**FULL TEXT READ** (scanned PDF from https://www.ism.ac.jp/editsec/aism/pdf/040_4_0693.pdf, text layer plus
page images pp. 693-705 rendered and read; acquired 2026-09-02, local copy in the session scratchpad
`sibuya1988.txt`, `sibuya_p-03..13.png`).

* **Theorem 3.1 (credited to Lieb (1968)).** "The following sequences are strictly decreasing for any
  n = 3, 4, ...;
  (A1)  ((m-1)/(n-m+1)) [n m]/[n m-1],  m = 2, ..., n;
  (A1') (m/(n-m)) [n n-m]/[n n-m+1],  m = 1, ..., n-1;
  (A2)  ((m-1)m/(n-m+1)) {n m}/{n m-1},  m = 2, ..., n;
  (A2') (m/((n-m)(n-m+1))) {n n-m}/{n n-m+1},  m = 1, ..., n-1."
  Proof text: "The first half (A1) and (A1') is equivalent to
  (3.1)  (m-1)(n-m) [n m]^2 > m(n-m+1) [n m-1][n m+1],  m = 2, ..., n-1."
  This is exactly strict Newton for (x+1)...(x+n-1); in our normalisation it is p_j^2 > p_{j-1}p_{j+1}
  with no gap. "The second half (A2) and (A2') is equivalent to
  (3.2)  (m-1)(n-m) {n m}^2 > (m+1)(n-m+1) {n m-1}{n m+1},  m = 2, ..., n-1."
  Newton for the Bell polynomial/x would give the right side with factor m(n-m+1); (3.2) has (m+1)(n-m+1),
  i.e. **Newton times (m+1)/m, proved by induction on n via the recurrence** (pp. 697-698; both proofs
  are displayed as explicit sums of nonnegative terms). Checked here exactly for n = 4..10: the normalised
  ratio minus (m+1)/m is positive for every m (smallest margin 0.089 at n = 10).
  Remarks pp. 698-699: "Lieb (1968) showed (3.1) from Newton's inequality and the fact that the generating
  function of ([n m])_{m=1}^n has only real roots. He showed also the sequence
  (3.3)  ((m-1)/(n-m+1)) {n m}/{n m-1},  m = 2, ..., n,
  strictly decreasing. The inequality (A2) slightly improves this result."
* **Eq. (3.4), p. 699, NOT PROVED -- the hit.** Verbatim: "Numerically, it is suggested that the sequence
  (3.4)  ((m-1)(2n+m)/(n-m+1)) [n m]/[n m-1],  m = 2, 3, ..., n-1,
  is strictly decreasing and the same for m = n-1 and n, and that
  (3.5)  ((m-1)^2/(n-m+1)) {n m}/{n m-1},  m = 2, 3, ..., n,
  is strictly decreasing."
  Translation of (3.4) (this note, verified exactly for n = 4..12 with the row recurrence in exact
  rationals, scratchpad `check_sibuya.py`): s_m > s_{m+1} is equivalent to
  c(n,m)^2/(c(n,m-1)c(n,m+1)) > [m(n-m+1)/((m-1)(n-m))] * (2n+m+1)/(2n+m), i.e. with j = n-m
      p_j^2/(p_{j-1} p_{j+1}) >= 1 + 1/(3n - j),   j = 1, ..., n-2,
  and the last pair (m = n-1, n; i.e. j = 1) is an EQUALITY for every n tested (s_{n-1} = s_n), which is
  why Sibuya separates "the same for m = n-1 and n". Hence Sibuya's (3.4) implies
      n (p_j^2/(p_{j-1}p_{j+1}) - 1) >= n/(3n-1) > 1/3 for all j, with the minimum at j = 1,
  which is precisely tonight's finding (our "m+1" is Sibuya's n). Sibuya's form is slightly stronger for
  j >= 2 (bound 1/(3n-j) instead of 1/(3n-1)).
  Hypotheses: n >= 3 integer; nothing else. Our parameters: spectrum {1..m} = Sibuya's n = m+1 --
  satisfied. Spectra {k^2}, {k^3}: NOT covered by any statement in the paper.
* **Corollary 3.3, eq. (3.6), PROVED** (two-sided ratio bound worth keeping for the ladder): for 2 <= m <= n,
  ((n-m+1)/((m-1)(n-1))) H_{n-1} >= [n m]/[n m-1] >= 2(n-m+1)/((m-1)n), H_k the harmonic number; left
  equality at m = 2, right equality at m = n, strict for 2 < m < n. (3.7) is the analogue for {n m} with
  (2^n - 2)/(n-1) on the left; remark p. 701: (3.7), (3.7') improve (4.17)-(4.20) of Neuman (1985).
* Theorems 3.2-3.5 (pp. 701-705): strict TP2 of both triangles (B1), (B2); (C1) (1/n)[n+1 m]/[n m] and (C2)
  {n+1 m}/{n m} strictly decreasing in n; (D1) (n-1)[n-1 m][n m] > n[n-1 m-1][n m+1]; (E1)
  (n-1)^2[n-1 m][n m-1] > n^2[n-1 m-2][n m+1]; (F1), (G1) [n m]^2 >= [n-1 m-1][n+1 m+1]; (D2), (E2) second-kind
  analogues. These are mixed-index inequalities, not row-Newton gaps.
* Later strengthening: **none found.** Forward citations (OpenAlex, 33 works, all listed and scanned by
  title; acquired 2026-09-02): Pitman 1997; Qi 2014, 2016 (second kind); Ferroni 2022; Yamato 1997, 2001;
  Sibuya 1991, 1993, 1997, 2005, 2014; Kabluchko-Marynych-Sulzbach 2016 (arXiv 1609.03798, full text: Sibuya
  appears only in the bibliography); Kabluchko-Steigenberger 2022, 2023; Genitrini et al. 2013, 2014, 2016;
  Tsukuda 2018, 2020; Finner-Roters 1993; Hong-Zhang 2020/2021 (arXiv 2008.10069, full text: they quote
  and use (3.6), the harmonic-number ratio bound, for the Nekrasov-Okounkov unimodality problem); Cramer
  2000; Nishimura-Sibuya 1997; Port 1994; Luca-Stanica 2011; Heim-Neuhauser 2021; Barranco-Chamorro 2015;
  Miegielsen 2014; Mansour 2012; Basar 2025; Ihm 2000. None restates or proves (3.4).

### 2. Lieb, E. H. (1968). Concavity properties and a generating function for Stirling numbers.
J. Combin. Theory 5, 203-206. DOI 10.1016/S0021-9800(68)80057-2.
**ABSTRACT ONLY -- theorem not verified.** ScienceDirect returned a bot wall (HTML) for both PDF endpoints;
Springer Selecta reprint redirected to an auth page; Semantic Scholar has the abstract elided by the
publisher. Abstract as returned by search: "The Stirling numbers of the first kind, S(N,k), and of the
second kind, sigma(N,k), are shown to be strongly logarithmically concave as functions of k for fixed N.
This result is stronger than the unimodality conjecture which was heretofore proved only for the Stirling
numbers of the second kind by Harper." Content of the theorems is known here SECOND-HAND through Sibuya's
Theorem 3.1 attribution: Lieb's "strongly log-concave" for the first kind is (3.1), i.e. strict Newton
with the binomial factor and no further gap; for the second kind it is (3.3). Cited also by Jalowy-
Kabluchko-Marynych 2025 as "[41], formula (3)" for the log-concavity of Bernoulli-sum laws.

### 3. Hammersley (1951), Erdos (1953), Kurtz (1972), Neuman (1985), Bona (textbook), Sagan (1988).
**NOT CONSULTED (primary texts not retrieved); second-hand only.** Per Sibuya, Corollary 3.2 remark:
"Hammersley (1951) and Erdos (1953) showed [n m]/[n m-1] != 1 for any n and m" (no two consecutive
first-kind numbers are equal), and Kurtz (1972) "gives log-concavity conditions for triangular arrays
defined by a recurrence formula" (a sufficient condition, no constant). Neuman (1985) proved (C2) via
log-concavity of symmetric means and had ratio bounds (4.17)-(4.20) for the second kind, improved by
Sibuya (3.7). Sagan 1988 and Bona are cited in Abdesselam 2020 for injective/real-root proofs of plain
log-concavity. None is reported anywhere as containing a c/N gap.

### 4. Niculescu, C. P. (2000). A new look at Newton's inequalities. JIPAM 1(2), Art. 17.
**FULL TEXT READ** (scratchpad `niculescu2000.txt`). Theorem 1.1 (Newton, Maclaurin): E_k^2 > E_{k-1}E_{k+1}
"unless all entries of F coincide" (strictness only). Eq. (1.5) (credited to Rosset [20], "cubic Newton
inequalities"): 6E_kE_{k+1}E_{k+2}E_{k+3} + 3E_{k+1}^2E_{k+2}^2 >= 4E_kE_{k+2}^3 + E_k^2E_{k+3}^2 +
4E_{k+1}^3E_{k+3}, rewritten as 4(E_{k+1}E_{k+3} - E_{k+2}^2)(E_kE_{k+2} - E_{k+1}^2) >= (E_{k+1}E_{k+2} -
E_kE_{k+3})^2 -- "strictly stronger than (1.1)", but it is the higher-order Turan condition, not a gap
constant. Theorem 2.1: E_{j alpha + k beta} >= E_j^alpha E_k^beta (interpolation). Section 5: P'^2 >=
(n/(n-1)) P P'' for real-rooted P. **No quantitative gap in terms of the roots anywhere in the paper.**

### 5. Rosset, S. (1989). Normalized symmetric functions, Newton's inequalities and a new set of stronger
inequalities. Amer. Math. Monthly 96(9), 815-819. DOI 10.1080/00029890.1989.11972286.
**ABSTRACT ONLY (403 at publisher)**; content known through Niculescu's quotation of (1.5) above. Cubic
(discriminant) inequalities; no c/N constant.

### 6. Sharma, R., Sharma, A., Saini, R., Kapoor, G. (2017). Means, moments and Newton's inequalities.
arXiv:1702.04665 [math.ST]. **FULL TEXT READ** (scratchpad `sharma2017.txt`). Abstract: "Newton's
inequalities and the related Maclaurin's inequalities provide several refinements of the fundamental
Arithmetic mean - Geometric mean - Harmonic mean inequality in terms of the means and variance of positive
real numbers." Proof of Theorem 2.2 uses 2C_2 = C_1^2 - sum x_i^2 (Newton's identity (1.6)), i.e. exactly
S_1^2 - S_2 = s^2/(n-1) with s^2 the population variance -- the j = 1 gap. Companion: Sharma-Bhandari
(2015), Skewness, kurtosis and Newton's inequality, Rocky Mountain J. Math. 45(5), 1639-1643, DOI
10.1216/RMJ-2015-45-5-1639, **ABSTRACT ONLY**: "an inequality related to Newton's inequality provides one
more relation between skewness and kurtosis" (the k = 2 Newton inequality in terms of the third and
fourth moments). **These are the only sources that write a Newton gap in terms of the root variance, and
only for j = 1 (and the moment form of j = 2).** Nothing uniform in j.

### 7. Arizmendi, O., Perales, D. (2018). Cumulants for finite free convolution. J. Combin. Theory Ser. A
155, 244-266. arXiv:1611.06598. **FULL TEXT READ** (scratchpad `arizmendi_perales.txt`). Proposition 3.4
(coefficient-cumulant formula), eqs. (3.2)-(3.3), for p(x) = sum_i x^{d-i}(-1)^i a_i with a_i normalised by
(d)_i/i! (their a-tilde_i = e_i/C(d,i), our p_i up to sign); Theorem 4.2, eqs. (4.4)-(4.5)
(moment-cumulant). **No statement about Newton's inequalities or log-concavity of the a-tilde_i appears in
the paper.** The n = 2 specialisation of (3.3) gives kappa_2 proportional to a-tilde_1^2 - a-tilde_2, i.e.
the j = 1 Newton gap; this reduction is mine and was not checked line by line against the paper. Marcus,
"Polynomial convolutions and (finite) free probability" (ff_main.pdf, full text grep for Newton /
log-concave / ultra / variance): only Newton's identities used; no coefficient inequality with a gap.
Arizmendi-Perales-et-al. later works (2108.08489, 2412.20488, 2606.10870): titles/abstracts only; nothing
on Newton gaps found in search.

### 8. Jalowy, J., Kabluchko, Z., Marynych, A. (2025). Zeros and exponential profiles of polynomials I.
arXiv:2504.11593. **FULL TEXT GREP** (scratchpad `jkm1.txt`). Newton's inequality is used only
qualitatively ("(P[S_n = k]) is a log-concave sequence ... This follows from Newton's inequality as
explained in [41], formula (3)", [41] = Lieb 1968); the exponential profile g is proved "infinitely
differentiable and strictly concave" in the limit (Theorem 2.2 area). Asymptotic strict concavity of g is
the limiting shadow of a Newton gap but carries no finite-N constant. Part II (arXiv 2509.11248) not read.

### 9. Pitman, J. (1997). Probabilistic bounds on the coefficients of polynomials with only real zeros.
J. Combin. Theory Ser. A 77, 279-303. DOI 10.1006/jcta.1997.2747. **FULL TEXT READ** (UC Berkeley Tech.
Report 453, https://statistics.berkeley.edu/sites/default/files/tech-reports/453.pdf, pp. 7-8 as images).
Section "Bounds for Consecutive Ratios", eq. (20): if A(z) = sum a_k z^k has non-negative coefficients and
only real zeros then theta(k + 1/(k+2)) <= a_k/a_{k+1} <= theta(k + 1 - 1/(n-k+1)), where theta(x) is the
unique positive solution of theta A'(theta)/A(theta) = x, derived from Darroch's rule for the mode, eq.
(13). This is a QUANTITATIVE ratio control (implicit through the tilt) valid for both Stirling kinds and
Eulerian numbers, but it bounds each ratio a_k/a_{k+1} separately, not the second difference of
log a_k; no Newton-gap constant. Cites Sibuya [110] and Lieb [82].

### 10. Abdesselam, A. (2020). A local injective proof of log-concavity for increasing spanning forests.
arXiv:2012.14330. **FULL TEXT READ (introduction and Theorems 1.1-1.2).** Only log-concavity and strong
x-log-concavity; ultra-log-concavity is cited as "[25, Theorem 2]" (Stanley 1989). No constant.

### 11. Higher-order Turan inequalities.
arXiv API, all fields "higher order Turan" AND Stirling: **0 results** (2026-09-02). Chen-Jia-Wang 2019
(Trans. AMS, DOI 10.1090/tran/7707) is the partition function; L. X. W. Wang, "Higher order Turan
inequalities for combinatorial sequences", Adv. Appl. Math. 110 (2019) (abstract-only via search): Motzkin,
Fine, Domb-type sequences, not Stirling. Negative.

### 12. Eulerian numbers, binomial-type triangles, second kind beyond Sibuya.
Searches ("Eulerian numbers" + Newton / ultra-log-concave + strict/quantitative; Bell polynomials + Turan +
ratio bound; OpenAlex "Turan inequalities" + "Stirling numbers" -> 1 irrelevant hit; OpenAlex "Newton's
inequalities" + "Stirling numbers" -> 25 works, all listed, none with a gap; arXiv all-fields "Stirling
numbers" AND Newton AND inequalit* -> 0) returned only plain log-concavity / ultra-log-concavity
(= Newton) statements. Johnson-Goldschmidt 2006 (arXiv math/0502548) reprove log-concavity of S(n,k) and
Eulerian numbers, no constant. Adell 2022 (JCTA 192, 105669; abstract only): explicit, asymptotically sharp
UPPER bounds on the values c(n,m), not on ratios. **Negative for all three families.**

---

## WHAT THIS MEANS FOR THE CLAIM

* Cite as prior art: **Sibuya 1988, eq. (3.4)** (conjecture, numerically suggested, first kind) and the
  classical identity for j = 1. Do not describe the 1/3 constant as new; describe a PROOF of the
  uniform-in-j statement as new, if and when it exists.
* The literature form of the conjecture is stronger than ours by the j-dependence 1/(3n-j); a proof of our
  weaker "min at j = 1" form would still be the first proof of any first-kind Newton gap.
* Sibuya's proved (3.2) shows the second-kind gap is O(1/j), not O(1/N): the two Stirling kinds are not
  analogous here, and a second-kind argument will not transfer.
* Two proved tools worth carrying: Sibuya (3.6) (harmonic-number two-sided ratio bound for c(n,m)) and
  Pitman (20) (tilt bounds on consecutive ratios of any real-rooted non-negative polynomial).

## ACQUISITION RECORD (all 2026-09-02)
ISM scan of Sibuya 1988 (full); arXiv 1702.04665, 1611.06598, 2504.11593, 2012.14330, 1609.03798,
2008.10069 (full text via pdftotext); Pitman TR 453 (full); Niculescu JIPAM 2000 (full, file already in
scratchpad); Marcus ff_main.pdf (full grep). Abstract-only: Lieb 1968, Rosset 1989, Sharma-Bhandari 2015,
Adell 2022, Wang 2019. Not retrieved: Hammersley 1951, Erdos 1953, Kurtz 1972, Neuman 1985, Bona, Sagan
1988. Exact checks: scratchpad `check_sibuya.py` (row recurrences in exact rationals, n <= 12; light).
