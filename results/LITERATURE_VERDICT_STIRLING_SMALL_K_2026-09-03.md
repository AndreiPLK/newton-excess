# Literature verdict: uniform, explicit-error asymptotics for c(n,k), k small, and Newton-type gaps

Written: Thu Sep  3 05:35:12 RDT 2026 (output of `date` at time of writing)
Branch: claude/handoff-markdown-review-anzwo5 (HEAD 4be6cc8)
Author: literature search agent for the qg-bootstrap Newton-excess programme
Companion to: `LITERATURE_VERDICT_STIRLING_NEWTON_GAP_2026-09-02.md` (Sibuya (3.4) is the conjecture;
nothing proves it) and `LITERATURE_VERDICT_RATIO_LOG_CONCAVITY_2026-09-02.md`.

**Question (structural agent, 2026-09-03).** For the unsigned Stirling numbers of the first kind
c(n,k) = [n k] = (n-1)! e_{k-1}(1, 1/2, ..., 1/(n-1)) in the regime k small (k <= ~200, or
k <= C log n), is there a uniform-in-n asymptotic with EXPLICIT (numerical-constant) error terms, and
does anyone prove a quantitative Newton gap c(n,m)^2/(c(n,m-1)c(n,m+1)) >= Newton * (1 + c/n)
from the harmonic-number (log-gamma) representation of the row? Five sub-questions, answered below.

**Notation dictionary (needed: no source uses ours).** All sources below write s(n,m) or [n m] for the
number of permutations of n with m cycles, x(x+1)...(x+n-1) = sum_m s(n,m) x^m, so s(n,m) =
e_{n-m}(1,...,n-1) = (n-1)! e_{m-1}(1, 1/2, ..., 1/(n-1)). "k small" in the question is m small in every
source (few cycles): this is the Jordan / Moser-Wyman / Wilf / Hwang regime. In the p_j normalisation of
the 2026-09-02 verdict (N = n-1 roots {1..N}, p_j = e_j(1..N)/C(N,j)), m small is j = N-(m-1) near N; the
Newton ratio is invariant under x -> 1/x, so p_{N-k}^2/(p_{N-k-1}p_{N-k+1}) = q_k^2/(q_{k-1}q_{k+1}) with
q_k := e_k(1/1, ..., 1/N)/C(N,k), k = m-1. Sibuya's (3.4) in this regime reads
q_k^2/(q_{k-1}q_{k+1}) >= 1 + 1/(3n - (N-k)) = 1 + 1/(2n+1+k). At k = 1 it is an exact identity:
q_1^2/q_2 = (1 - 1/N)/(1 - H_N^{(2)}/H_N^2) (from e_1 = H_N, e_2 = (H_N^2 - H_N^{(2)})/2), so the TRUE
excess in the few-cycle regime is of order zeta(2)/(log n)^2, enormously larger than the 1/n that (3.4)
asks for. The difficulty in this regime is therefore not the size of the gap but controlling the
expansion remainder UNIFORMLY in m with explicit constants. That is exactly what the literature does not
supply -- see the verdict.

---

## VERDICT IN SEVEN LINES

1. **Hwang 1995 (full text read, preprint version): uniform in 1 <= m <= eta log n for any eta > 0,
   but every error term is an O-term with unspecified constants.** Theorem 1: s(n,m)/n! = (1/n)
   sum_{0<=k<=nu} Pi_{m,k}(log n) n^{-k} + O((log n)^m /(m! n^{nu+2})). Theorem 2: a Selberg-type uniform
   expansion with remainder Z_mu(m,n) = O(K_mu (m-1)^{mu/2}/(log n)^mu + (log n)^m/(m! n)). No numerical
   constant anywhere. The published JCTA version (DOI 10.1016/0097-3165(95)90010-1) could not be
   retrieved (ScienceDirect 403); the September 1994 LIX preprint (7 pp.) was read in full.
2. **The only EXPLICIT-constant small-m expansion in print is Moser-Wyman 1958 as quoted by Adell 2022,
   eq. (2):** s(n+1,m+1) = ((-1)^{n-m} n! (H_n^{(1)})^m / m!) (1 - m(m-1)H_n^{(2)}/(2 (H_n^{(1)})^2) + E_{n,m}),
   |E_{n,m}| <= 2 (e m / H_n^{(1)})^3. SECOND-HAND (Moser-Wyman is paywalled, 403 at both OUP and Wiley).
   For a Newton gap this remainder is useless in practice: the gap it must resolve is ~ H^{(2)}/H^2, the
   remainder is ~ 2e^3 m^3/H^3, so it decides only when m^3 < ~0.04 H_n, i.e. for m = 2 only for n beyond
   e^200 (this note's arithmetic, section 2 below).
3. **Wilf 1993 (full text read, author's PostScript): complete asymptotic series for FIXED k; error
   O((log n)^{k-2}/n); Theorem 2 gives the full series to any order in 1/n; no explicit constants;
   validity "in no event can ... extend beyond about k ~ C log n / log log n" (his words).**
4. **Arratia-DeSalvo 2016/17 (full text): the only completely effective (explicit, all finite n)
   two-sided bounds for the first kind -- but at the OTHER end, s(n, n-k) with k = O(sqrt n)
   (many cycles).** Not our regime. Adell 2022 (full text): explicit UPPER bounds only, valid for all
   m, e.g. |s(n+1,m+1)| <= n! (log n)^m/m! (1 + m/log n) (Theorem 2.1, eq. (8)); no lower bound, no ratio.
5. **Temme 1993 (full text): one-term uniform saddle-point approximation (3.5) for the first kind;
   accuracy shown by tables (max relative error 0.0082 at n = 10), no error bound of any kind.**
   Louchard 2010 (full text): central region (Theorem 2.1, Edgeworth-type in powers of 1/sqrt(A_n)) and
   j = n - n^alpha, alpha > 1/2 (Theorem 3.1); all O-terms; not the small-m regime.
6. **No paper proves log-concavity or a Newton-type gap from the harmonic-number (Bell-polynomial)
   representation e_k(1/i) = Y_k(H, -H^{(2)}, 2! H^{(3)}, ...)/k!.** The representation is standard
   (Adamchik 1997; Wilf's (1)-(4) is exactly it) but is used only for identities and asymptotics.
   Sibuya (3.6) remains the only proved two-sided ratio bound with harmonic numbers; nothing sharper found.
   In the n-direction the state of the art is a CONJECTURE (Liang-Sagan 2024, Conjecture 8.9).
7. **Consequence.** An explicit-constant, uniform-in-m (m <= C log n) expansion of c(n,m) with a
   remainder small compared to H^{(2)}/H^2 does not exist in the literature; deriving one (from Hwang's
   Theorem 1 machinery, whose only non-explicit step is the singularity-analysis remainder in his eq. (6)
   for Gamma(n+w)/(Gamma(w) n!), or from the exact Bell-polynomial form) would be new, and the gap
   inequality proved from it would be the first proof of any first-kind Newton gap (cf. 2026-09-02 verdict).

---

## 1. Hwang, H.-K. (1995). Asymptotic expansions for the Stirling numbers of the first kind.
J. Combin. Theory Ser. A 71(2), 343-351. DOI 10.1016/0097-3165(95)90010-1.
**FULL TEXT READ -- of the LIX Ecole Polytechnique preprint "NOTE, September 1, 1994" (7 pp.), obtained
from CiteSeerX (https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=05a0f4171d12cd494a80267b2924bb3fe6d81d77),
which is the document Semantic Scholar lists as the paper's PDF.** The journal version was not
retrievable (ScienceDirect returned 403 for the PDF endpoint); theorem numbering and statements below
are from the preprint and may differ in the published version by editorial changes only (the preprint
already cites Temme 1993 and Wilf 1993 with journal data). Local copy: scratchpad `hwang1995_preprint.pdf`,
sha256 a604153e..., acquired 2026-09-03.

Setting: s(n,m) := [w^m](w(w+1)...(w+n-1)), 1 <= m <= n. g_k := [w^k] 1/Gamma(1+w) (g_0 = 1, g_1 = gamma).

* **Theorem 1 (verbatim).** "For any eta > 0 and nu in N, the Stirling numbers of the first kind s(n,m)
  satisfy asymptotically
      s(n,m)/n! = (1/n) sum_{0<=k<=nu} Pi_{m,k}(log n)/n^k + O( (log n)^m / (m! n^{nu+2}) )   (n -> infinity),  (2)
  uniformly for 1 <= m <= eta log n, where Pi_{m,k}(X) are polynomials in X of degree m-1."
  With (3): Pi_{m,k}(X) := sum_{0<=l<=m-1} g_{k,m-1-l} X^l / l!, g_{k,l} := [w^l] pi_k(w)/Gamma(1+w),
  pi_k(w) = ((-1)^k/k!) B_k^{(w)} (w-1)(w-2)...(w-k), B_k^{(w)} := k! [t^k] (t/(e^t-1))^w (generalized
  Bernoulli numbers); pi_1(w) = w(w-1)/2, pi_2(w) = w(w-1)(w-2)(3w-1)/24, pi_3(w) = w^2(w-1)^2(w-2)(w-3)/48.
* **Corollary 1 (verbatim).** "Formula (1) holds for m in the range 1 <= m <= eta log n for any eta > 0,
  provided that the error term is replaced by O((log n)^m m!^{-1} n^{-2})." Here (1) is Wilf's
  s(n,m)/n! = (1/n)( g_0 (log n)^{m-1}/(m-1)! + g_1 (log n)^{m-2}/(m-2)! + ... + g_{m-1} ) + O((log n)^{m-2}/n^2).
* **Theorem 2 (verbatim).** "The Stirling numbers of the first kind satisfy, uniformly with respect to m,
  2 <= m <= eta log n, eta > 0, (r := (m-1)/log n)
      s(n,m)/n! = ((log n)^{m-1}/(n (m-1)!)) ( 1/Gamma(1+r) + sum_{2<=k<mu} g^{(k)}(r) tau_k(m)/(log n)^k + Z_mu(m,n) )  (5)
  where tau_k(m) is a polynomial in m of degree floor(k/2), cf. (10), tau_k(m) := [w^k](e^{-w}(1+w))^{m-1}
  (k = 0, 1, ..., m-1), g(w) := Gamma^{-1}(1+w), mu is any fixed positive integer such that 2 <= mu <= m, and
      Z_mu(m,n) = O( K_mu (m-1)^{mu/2}/(log n)^mu + (log n)^m/(m! n) ),
  with |g^{(mu)}(z)| <= K_mu for |z| <= eta."
  tau_2(m) = -(m-1)/2, tau_3(m) = (m-1)/3, tau_4(m) = (m-1)(m-3)/8, tau_5(m) = -(m-1)(5m-11)/30,
  tau_6(m) = -(m-1)(3m^2-32m+53)/144; tau_k(m) = L_k^{(-k-1)}(m-1) (Laguerre), with recurrence
  (k+1) tau_{k+1}(m) = -k tau_k(m) + (1-m) tau_{k-1}(m).
* Proof structure (Section 2): from singularity analysis, eq. (6),
  w(w+1)...(w+n-1)/n! = [z^n](1-z)^{-w} = (n^{w-1}/Gamma(w)) (1 + sum_{1<=k<=nu} pi_k(w)/n^k + O(n^{-nu-1})),
  "the O being uniform with respect to w, |w| <= eta"; then Cauchy's formula on |w| = m/log n and an
  explicit Bessel-type estimate (7): int_{-pi}^{pi} e^{m cos t} dt < 2 sqrt(pi) m^{-1/2} e^m + pi.
  Section 3 (Theorem 2): Selberg's saddle-point variant; the remainder Y_mu(X) is bounded explicitly
  through K_mu, Gamma((mu+1)/2), and Stirling's inequality e^{m-1}/(m-1)^{m-1/2} < e^{1/12}sqrt(2pi)/(m-1)!.
* **Are the error terms explicit?** NO. Both theorems carry O(.) with unspecified absolute constants. The
  ONLY non-effective ingredient in the proof of Theorem 1 is the constant in the singularity-analysis
  remainder of eq. (6) (Flajolet-Odlyzko 1990, p. 220, proved there with z = 1 + t/n; Hwang uses z = e^{-t});
  everything downstream is written with explicit inequalities. In Theorem 2 the implied constants in
  Z_mu come from (8)-(9) and are likewise traceable but not stated.
* Hypotheses vs ours: 1 <= m <= eta log n -- satisfied for k <= C log n (any C, but the O-constants
  depend on eta); k <= 200 fixed with n -> infinity -- satisfied. What fails: no numerical constant, so
  nothing can be CERTIFIED at a finite n from this paper as written.
* Later strengthening: Hwang's own "Large deviations for combinatorial distributions II: local limit
  theorems" (in preparation in 1994; Ann. Appl. Probab. 1998) and the quasi-power / mod-Poisson framework
  (Kabluchko-Marynych-Pitters 2022, below) restate the same expansion with O(1/n) speed; none makes the
  constants explicit. Louchard 2010 p. 167: "An asymptotic expansion for j = O(1) is given in Wilf [14],
  which has been extended to the range j = O(ln n) by Hwang [6]" -- still the reference of record in 2010
  and in Adell 2022.

Remark from the preprint worth keeping (p. 2): "Theorem 1 says merely that sum_k Pi_{m,k}(log n) n^{-k-1}
is an asymptotic expansion for s(n,m)/n!, thus whenever we use (2), each 'coefficient' Pi_{m,k}(log n)
should be calculated as a whole ... since the Pi_{m,k}'s are themselves asymptotic expansions only when
m = o(log n)" -- i.e. for m ~ c log n the expansion in powers of 1/n is NOT an expansion in powers of
1/log n; the two scales must be separated (this is why Theorem 2 exists).

## 2. Moser, L., Wyman, M. (1958). Asymptotic development of the Stirling numbers of the first kind.
J. London Math. Soc. 33, 133-146. DOI 10.1112/jlms/s1-33.2.133.
**NOT RETRIEVED (403 at academic.oup.com and londmathsoc.onlinelibrary.wiley.com). SECOND-HAND from three
independent full texts:**
* Wilf 1993 (Section 4, eq. (9)), "their formula for the range of small k begins as
  (1/(n-1)!) [n k] = sigma_n(1)^{k-1}/(k-1)! { 1 - C(k-1,2) sigma_n(2)/sigma_n(1)^2 + 2 C(k-1,3) sigma_n(3)/sigma_n(1)^3
  + 3 C(k-1,4) sigma_n(2)^2/sigma_n(1)^4 - 2 C(k-1,4)?? sigma_n(4)/sigma_n(1)^4 + ... }" [the last binomial is
  unreadable in the PS-to-text conversion; sigma_n(s) := sum_{j<=n} j^{-s}]. Wilf: "Moser and Wyman gave a
  comprehensive study ... covered the entire range 1 <= k <= n by breaking it into three subranges ...
  The formulas that they obtain for the low end of the range, i.e., for fixed k, are however not phrased in
  terms of ... powers of n and of log n. Instead the formulas are expressed as series of polynomials in the
  partial sums of the zeta function, containing p(k) monomials."
* Arratia-DeSalvo (Section 6.2): "Jordan [9] is credited by Moser and Wyman [14] for |s(n,m)| ~
  (n-1)!(log n + gamma)^{m-1}/(m-1)!. Moser and Wyman [14] extended this first-order asymptotic formula to
  values of m = o(log n)." And (p. 2-3) for the middle range, quoting MW pp. 142-143: "One of the defects of
  formula (3) is the fact that we have given no estimate of the error ... we have been unable to give even
  such a crude estimate of the error involved in the use of (3)" -- (3) being |s(n,m)| = Gamma(n+R)/
  (R^m Gamma(R) sqrt(2 pi H)) (1 + O(m^{-1})), R the tilt with sum_{k<n} R/(R+k) = m, H = m - sum R^2/(R+k)^2.
* **Adell 2022, eq. (2), read from the page image (p. 2):** "Moser and Wyman [14] gave the estimate
      s(n+1, m+1) = ((-1)^{n-m} n! (H_n^{(1)})^m / m!) ( 1 - m(m-1) H_n^{(2)} / (2 (H_n^{(1)})^2) + E_{n,m} ),   (2)
  where |E_{n,m}| <= 2 (e m / H_n^{(1)})^3. Of course, formula (2) works in the range m = o(H_n^{(1)}) =
  o(log n). The same authors showed that more exact terms can be added in formula (2)." (H_n^{(a)} :=
  sum_{j<=n} j^{-a}.)
  **This is the only explicit-constant small-m statement located anywhere.** It must be confirmed against
  the primary before use (the constant 2e^3 and the power 3 are Adell's transcription).
* **What it buys for the Newton gap (this note's arithmetic, not from any source).** Write H = H_n^{(1)},
  H2 = H_n^{(2)}, A_m := 1 - m(m-1)H2/(2H^2) + E_{n,m}. Then
      s(n+1,m+1)^2 / (s(n+1,m) s(n+1,m+2)) = ((m+1)/m) * A_m^2/(A_{m-1}A_{m+1}),
  and since the second difference of m(m-1) is 2, A_m^2/(A_{m-1}A_{m+1}) = 1 + H2/H^2 + O(H2^2 m^2/H^4) + (error
  terms), i.e. the excess over Newton's bound ((m+1)/m)(1 + 1/(n-m)) is ~ H2/H^2 ~ zeta(2)/(log n)^2, while
  Sibuya (3.4) asks only for the factor 1 + 1/(2n+m+3). The obstruction is the remainder: with
  |E| <= 2(e(m+1)/H)^3 in each of the three factors one needs roughly 8 e^3 (m+1)^3/H^3 < H2/H^2, i.e.
  (m+1)^3 < ~ H/160 -- false for every m >= 1 at every n below e^{300}. So Moser-Wyman's explicit remainder,
  even if confirmed, cannot certify the gap for any m >= 2 at any n reachable by the ladder. An explicit
  remainder of relative size o(H2/H^2) uniformly in m is what is needed and is not in print.

## 3. Wilf, H. S. (1993). The asymptotic behavior of the Stirling numbers of the first kind.
J. Combin. Theory Ser. A 64, 344-349. DOI 10.1016/0097-3165(93)90103-F.
**FULL TEXT READ** (author's PostScript from http://www.math.upenn.edu/~wilf/website/stirling.ps, linked
from his reprints page; converted with ps2pdf/pdftotext; 6 pp.; scratchpad `wilf1993.ps`, sha256 7fd9c0de...).
* **Theorem 1 (verbatim).** "If [n k] is the (signless) Stirling number of the first kind, then for each
  fixed integer k > 1 we have
      (1/(n-1)!) [n k] = alpha_1 (log n)^{k-1}/(k-1)! + alpha_2 (log n)^{k-2}/(k-2)! + ... + alpha_k + O((log n)^{k-2}/n)
  where the alpha_j are the coefficients in the expansion 1/Gamma(z) = sum_{j>=1} alpha_j z^j" (alpha_1 = 1,
  alpha_2 = gamma, alpha_3 = (6gamma^2 - pi^2)/12, ...).
* Proof (Section 2, eqs. (1)-(4)): (1/n!) sum_k [n+1 k+1] x^k = prod_{j<=n}(1 + x/j) = exp(sum_s (-1)^{s-1}
  x^s sigma_n(s)/s), sigma_n(s) expanded by Euler-Maclaurin (eq. (2), with Bernoulli-number coefficients
  d(s,j)), and the identity exp(sum_{s>=2} (-1)^{s-1} x^s zeta(s)/s) = e^{-gamma x}/Gamma(1+x) (eq. (4)).
  **This IS the harmonic-number / log-gamma representation asked about in sub-question (3)**: the row
  generating function is exactly Gamma(n+1+x)/(Gamma(n+1)Gamma(1+x)) and its coefficients are polynomials
  in the sigma_n(s) = H_n^{(s)}. Wilf uses it only to derive the asymptotic series.
* **Theorem 2 (verbatim, structure).** "The complete asymptotic series for the Stirling numbers of the
  first kind, for fixed k, is obtained by equating the coefficients of x^k on both sides of
      (1/n!) sum_k [n+1 k+1] x^k ~ (e^{x log n}/(x Gamma(x))) exp( sum_{j>=1} phi_j(x)/n^j ),   (6)
  where phi_j(x) = (-1)^j x^{j+1}/(j(j+1)) + sum_{s=1}^{j} C(j-1,s-1) (-B_{j+1-s}/(s(j+1-s))) x^s ... (j = 1,2,...)"
  [the exact form of phi_j is partly garbled by the conversion; the displayed order-1/n^2 result (7) is
  (1/n!)[n+1 k+1] = sum_{r=0}^{k} (log n)^r/r! { alpha_{k+1-r} + (1/(2n))(-alpha_{k-r} + alpha_{k-r-1}) +
  (1/(24 n^2))(-2alpha_{k-r} - 3alpha_{k-r-1} + 2alpha_{k-r-2} + 3alpha_{k-r-3}) } + O(n^{-3+epsilon})].
* Explicit constants: NONE; fixed k only. Wilf (p. 5): "Our formulas above are considerably more explicit,
  though they have been proved here to hold only for constant k, and in no event can their validity extend
  beyond about k ~ C log n / log log n." Numerical table p. 5: [100 5]/99! = 21.1204415, eq. (7) to order
  1/n^2 gives 21.1204409 (error 3e-6 %); [50 9]/49! = 0.5591960, Theorem 1 alone is off by 9.2 %.
* Side remark (p. 5): "The first several of these polynomials v_k(y) [e^{xy}/(x Gamma(x)) = sum v_k(y) x^k]
  appear to have only real zeros, though I am not aware of any theorems that would assure this."
  Unresolved in the paper; not pursued here.
* Hypotheses vs ours: fixed k -- fails for k <= C log n (Hwang extends it); no constants -- fails for
  certification.

## 4. Arratia, R., DeSalvo, S. (2016/2017). Completely effective error bounds for Stirling numbers of the
first and second kinds via Poisson approximation. arXiv:1404.3007v4; Ann. Comb. 21(1) (2017) 1-24,
DOI 10.1007/s00026-017-0339-z. **FULL TEXT READ** (arXiv v4, 19 pp.; scratchpad `arratia_desalvo.pdf`).
* **Theorem 5 (verbatim, first-kind half).** "Suppose n >= 3 and n >= k >= 2. Let N := C(n,2), and define
  mu = mu_{n,k} := C(k,2) C(n,3) / C(N,2), [d_1, d_2 explicit rational functions of n, k, N -- eq. block p. 7],
  D_n^{(1,k)} := min(d_1, mu d_1, 1) ... Then ... C(N,k) e^{-mu}(1 - e^{mu} D_n^{(1,k)}) <= |s(n, n-k)| <=
  C(N,k) e^{-mu}(1 + e^{mu} D_n^{(1,k)})." "Note that the error term is D = O(k^3/n^2 ...)" [rendering
  garbled], "goes to 0 for k = O(sqrt n)".
* Proposition 1 (eqs. (8), (10)): for m = n - t n^a, 0 <= a <= 1/2, |s(n,m)| ~ C(n,2)^{n-m}/(n-m)! e^{-2 t^2
  n^{2a-1}/3 ...} matching Moser-Wyman for a < 1/2 and Louchard for a > 1/2 (Theorem 3, from Louchard).
* Regime: s(n, n-k), k <= O(sqrt n), i.e. MANY cycles (m = n-k). **Our regime (m small) is not covered:
  the error D is O(1) or larger there.** Recorded because it is the only "completely effective" first-kind
  bound in existence, and because its Section 6.2 is the cleanest historical map of the small-m literature
  (Jordan -> Moser-Wyman [m = o(log n)] -> Wilf [m = O(1), full series] -> Hwang [m = O(log n)]).
* Also from the same paper, second kind, eq. (1): Moser-Wyman's explicit hard error for S(n,m) in the
  range m = n - o(sqrt n) -- again the many-blocks end.

## 5. Adell, J. A. (2022). Explicit upper bounds for the Stirling numbers of the first kind.
J. Combin. Theory Ser. A 192, 105669. DOI 10.1016/j.jcta.2022.105669 (CC-BY).
**FULL TEXT READ** (Zaragoza repository https://zaguan.unizar.es/record/119973/files/texto_completo.pdf,
17 pp., pages 2-4 also read as images; scratchpad `adell2022.pdf`, sha256 8f961785...).
* Definitions (6)-(7): mu_n(t) = sum_{j<=n} t/(t+j), sigma_n^2(t) = sum_{j<=n} (t/(t+j))(1 - t/(t+j));
  tau = tau(n,m) > 0 the unique solution of mu_n(tau) = m (m = 1..n-1). "Such functions were introduced by
  Moser and Wyman [14] (see also Temme [18], Chelluri et al. [6], and Louchard [13]). Their respective
  probabilistic meaning, as the mean and the variance of a certain random variable having the
  Poisson-binomial distribution, will be established in Section 3."
* **Theorem 2.1 (verbatim).** "Let m = 1, ..., n-1 and let tau be as in (7). Then,
      |s(n+1, m+1)| <= (n! (log n)^m / m!) (1 + m/log n)                                              (8)
  and
      |s(n+1, m+1)| <= (Gamma(tau+n+1)/(Gamma(tau) tau^{m+1})) min( 1, 1/(sigma_n(tau) sqrt(2 pi (1 - log sigma_n(tau)/sigma_n^2(tau)))) + 1/sigma_n^2(tau) ).   (9)
  In addition, whenever n - m <= sqrt(m+1)/4, we have
      |s(n+1, m+1)| <= C(n+1, m+1) ((m+1)/2)^{n-m} (1 + 32 e^{1/6} (n-m)^2/(m+1))."                     (10)
* **Theorem 2.2 (verbatim).** "Let m = 1, ..., n-1. If n - m <= sqrt(m+1)/4, then
  | (-1)^{n-m} s(n+1,m+1) - C(n+1,m+1) ((m+1)/2)^{n-m} (1 + 5(n-m)(n-m-1)/(6(m+1))) |
  <= 2^7 e^{1/6} C(n+1,m+1) ((m+1)/2)^{n-m} ((n-m)/sqrt(m+1))^3."
* Text after 2.2: "(8) is asymptotically sharp, provided that m = o(log n). We point out that Hwang [11]
  extended formula (2) into an asymptotic series valid for m = O(log n)." And Chelluri-Richmond-Temme
  (2000), eq. (11): |s(n+1,m+1)| = (Gamma(tau+n+1)/(Gamma(tau) tau^{m+1} sigma_n(tau) sqrt(2 pi))) (1 + O(1/m))
  for n/log n <= m <= n - n^{1/3}; and Moser-Wyman [14, Lemma 4.1]: sigma_n(tau) -> infinity there.
  Practical guide (p. 5): use (8) for m <= (log n)^r, (9) for (log n)^r <= m <= n - n^s, (10) beyond,
  1/2 <= r < 1, 1/3 <= s < 1/2.
* Hypotheses vs ours: (8) valid for all n >= 2, m <= n-1 -- satisfied; **it is an UPPER bound only, and a
  ratio inequality needs a matching LOWER bound of relative precision o(H2/H^2); none is given.** (10) and
  Theorem 2.2 are the many-cycles end (n - m <= sqrt(m+1)/4), not ours; their explicit constants (32e^{1/6},
  2^7 e^{1/6}) show the style of bound that WOULD be needed at our end.
* Also read: Proposition 7.1 (zeta(m+2) tail bound from (8)); Theorem 7.2 (Comtet numbers). Method:
  Poisson-binomial representation s(n+1,m+1)/(...) = P(W(p_n) = m) with p_j = tau/(tau+j), and
  Shorgin/Roos-type explicit Poisson-binomial estimates (Lemmas 4.1, 4.2, 5.1, 5.2). **This is the same
  Poisson-binomial local-limit machinery the 2026-09-02 STUCK-protocol search found (Korolev-Zhukov,
  Roos, Shorgin): it bounds the density, not the second difference of its logarithm.**

## 6. Temme, N. M. (1993). Asymptotic estimates of Stirling numbers. Stud. Appl. Math. 89, 233-243.
DOI 10.1002/sapm1993893233. **FULL TEXT READ** (CWI repository https://ir.cwi.nl/pub/2304/2304D.pdf,
11 pp. scan; scratchpad `temme1993.pdf`).
* Section 3 (first kind): integral representation via phi(x) = ln Gamma(x+n+1) - ln Gamma(x+1) - m ln x,
  transformation (3.2) phi(x) = n ln(1+t) - m ln t + B with t_0 = m/(n-m), one-term approximation
  (3.5): (-1)^{n-m} s(n+1, m+1) ~ C(n,m) e^B g(t_0) ... with g(t_0) = sqrt( m(n-m)/(n x_0^2 phi''(x_0)) )
  [scan partly garbled]; and an exact finite "curious" expansion (-1)^{n-m} s(n+1,m+1) = e^B sum_k c_k C(n, m-k).
* Accuracy: tables only. "The maximal relative error now occurs at m = 3, and is 0.0082. For n = 20,
  n = 30, the maximal errors are: 0.0063 and 0.0053" (p. 239). Section 4: higher-order terms worked out for
  the SECOND kind only ("the treatment for the numbers of the first kind is more complicated").
* Verdict: uniform in m, one term, NO error bound. Not usable for certification.

## 7. Louchard, G. (2010). Asymptotics of the Stirling numbers of the first kind revisited: a saddle point
approach. DMTCS 12(2), 167-184. DOI 10.46298/dmtcs.501. **FULL TEXT READ** (episciences PDF, 18 pp.).
* Theorem 2.1 (central region, x := (j - M)/sigma, B_n := sqrt(A_n), A_n := ln n - pi^2/6 + gamma):
  Z_n(j) := [n j]/n! ~ (e^{-x^2/2}/(sqrt(2 pi) B_n)) { 1 + (x^3/6 - x/2)/B_n + (3x^2/8 - x^4/6 - 1/12 +
  x^6/72)/B_n^2 + (...)/B_n^3 + ... } -- Edgeworth-type; no error bound.
* Theorem 3.1 (j = n - n^alpha, 1 > alpha > 1/2; x = n^alpha, y = n^{1-alpha}): an explicit double series in
  1/x, 1/y with O(1/x^3) remainder. Many-cycles end, O-terms.
* p. 173 on Moser-Wyman: "For alpha < 1/2, Moser and Wyman (6.9) give an explicit asymptotic expression. For
  alpha > 1/2, they first compute in (4.52) the numerical solution z_n of S'(z_n) = 0 and give in (4.51) an
  asymptotic expression. This is rather precise: for n = 50, this gives a precision of order 10^{-4}."
* Verdict: not the small-m regime; no explicit constants.

## 8. Erdos, P. (1953). On a conjecture of Hammersley. J. London Math. Soc. 28, 232-236.
DOI 10.1112/jlms/s1-28.2.232. **FULL TEXT READ** (Renyi Institute archive https://users.renyi.hu/~p_erdos/1953-04.pdf,
scratchpad `erdos1953.pdf`). Notation: Sigma_{n,s} = e_s(1, ..., n); f(n) = largest maximising s.
* Theorem 1: "For sufficiently large n all the integers Sigma_{n,s}, 1 <= s <= n, are different." Proof by
  a prime p_k in (n/(k+1), n/k) dividing Sigma_{n,n-r} for r < k but not Sigma_{n,n-k}; second, elementary
  proof valid for n > 10^8 with k < e(log n + 1).
* Theorem 2 (Hammersley's conjecture): "The value of s which maximises Sigma_{n,s} is unique," proved for
  all n via primes in (n/(t+2), n/(t+1)), t = n - f(n), log n - 2 < t < log n, plus Hammersley's tables.
* From Hammersley 1951 (Proc. London Math. Soc. (3) 1, 435-452), quoted by Erdos: Hammersley proves
  uniqueness for n <= 188 and the mode formula (3): f(n) = n - [log(n+1) + gamma - 1 - zeta(2)/(log(n+1)
  + gamma - 1/2) + h/(log(n+1) + gamma - 1/2)^2] ... with -1.1 < h < 1.5 (scan partly garbled at the
  zeta term); hence (4): [log n - 2] < n - f(n) < [log n] for n > 188. Hammersley 1951 itself NOT RETRIEVED.
* Closing remark (p. 236, Erdos-Stone theorem): for any positive reals u_1 < u_2 < ... with sum 1/u_i =
  infinity and sum 1/u_i^2 < infinity, the maximising s satisfies f(n) = n - [sum_{i<=n} 1/u_i - sum 1/(u_i^2 ...)]
  (garbled) -- a mode formula for general spectra, relevant to the {k^2}, {k^3} spectra of the programme.
* Verdict: existence/uniqueness of the mode via divisibility; NOTHING about the size of the Newton ratio.
  (Heim-Neuhauser 2022, arXiv:2208.09928, Section 6.2, give a new proof of one mode "for infinitely many n"
  via Darroch's theorem -- also no quantitative gap.)

## 9. Log-concavity in the n-direction; ratio monotonicity; sharper ratio bounds (sub-questions 4-5)
* **Liang, J., Sagan, B. E. (2024). Log-concavity and log-convexity via distributive lattices.
  arXiv:2408.02782 (published Order 2026, DOI 10.1007/s11083-026-09735-2, not retrieved). FULL TEXT READ
  (arXiv v1).** Section 8.2, verbatim: "We have checked the following conjecture for 1 <= k <= n <= 100.
  Conjecture 8.9. Given k, there is an integer N_k such that (c(n,k))_{n>=0} is log-concave for n < N_k and
  log-convex for n >= N_k." **Open as of the search (no proof found; WebSearch 2026-09-03).** So in the
  n-direction not even qualitative log-concavity/log-convexity of c(n,k) at fixed k is settled.
* **Sibuya 1988, Corollary 3.3, eq. (3.6)** (read in full on 2026-09-02): for 2 <= m <= n,
  ((n-m+1)/((m-1)(n-1))) H_{n-1} >= [n m]/[n m-1] >= 2(n-m+1)/((m-1)n). Left equality at m = 2, right at
  m = n. **No sharper published two-sided bound on [n m]/[n m-1] was found** (OpenAlex searches "ratio
  Stirling numbers first kind inequality harmonic", "Stirling numbers first kind explicit bounds small m";
  WebSearch for Qi/Guo/Chen-type ratio bounds; the 33 forward citations of Sibuya listed on 2026-09-02).
  Hong-Zhang 2020 (arXiv:2008.10069, Section 3) use exactly (3.6), which gives for m >= 2H_n + 1
  [n+1, m+t+1]/[n+1, m+1] <= (H_n/m)^t <= 2^{-t}; their own Lemma 2.2 (sum_{m<=n} c_{m,k} = (pi^2/6)^k
  C(n,k)(1 + O(k^2 log n / n)), n >= k^2, for a partition-related c_{m,k}, NOT Stirling) is an O-term.
* **Pitman 1997 eq. (20)** (read 2026-09-02) remains the only general implicit ratio control:
  theta(k + 1/(k+2)) <= a_k/a_{k+1} <= theta(k + 1 - 1/(n-k+1)) for real-rooted A(z), theta(x) the tilt with
  theta A'(theta)/A(theta) = x. Per-ratio, no second-difference constant.
* **Monotonicity of c(n,k)/c(n,k+1) in k** is Sibuya's (A1) = strict Newton (Lieb 1968); in n it is
  Sibuya (C1): (1/n)[n+1 m]/[n m] strictly decreasing in n. Nothing beyond these.
* Second kind, for contrast: Canfield, E. R., Pomerance, C., "On the problem of uniqueness for the maximum
  Stirling number(s) of the second kind" (preprint header: "Integers: Electronic Journal of Combinatorial
  Number Theory 1 (2001), #Axx", received 2001-09-10, accepted 2001-12-05; final volume/article number NOT
  resolved this session -- OpenAlex search returned only later papers; full text read from
  https://math.dartmouth.edu/~carlp/MaxStirCorr.pdf),
  eq. (2): "S(n,k)^2 >= (1 + 3/k) S(n,k-1) S(n,k+1), 1 <= k <= n" by induction on the recurrence -- a 1/k
  gap of the same type as Sibuya's (3.2). No first-kind analogue exists (the first-kind recurrence
  [n+1 m] = n[n m] + [n m-1] does not propagate a constant gap the same way; this is the 2026-09-02
  finding that the two kinds are not analogous here).

## 10. Harmonic-number / Bell-polynomial representation (sub-question 3)
* The representation c(n, k+1)/(n-1)! = e_k(1, 1/2, ..., 1/(n-1)) = (1/k!) Y_k(H, -1! H^{(2)}, 2! H^{(3)},
  ..., (-1)^{k-1}(k-1)! H^{(k)}) (H^{(r)} := H_{n-1}^{(r)}; Y_k complete exponential Bell polynomial) is the
  content of Wilf's eqs. (1)-(4) and of Adamchik, V. (1997), On Stirling numbers and Euler sums, J. Comput.
  Appl. Math. 79, 119-130 (metadata only; not retrieved); also arXiv:1001.2835 (Connon, "Various
  applications of the (exponential) complete Bell polynomials", ABSTRACT ONLY via search), and Adell-Lekuona
  2019 (Adv. Differ. Equ. 398, "Explicit expressions and integral representations for the Stirling numbers.
  A probabilistic approach", metadata only).
* **No source uses this representation to prove log-concavity, a Turan inequality, or a Newton gap.**
  Searches: OpenAlex "elementary symmetric functions reciprocals integers inequality log-concave" (0 relevant
  of 12); arXiv API returned empty result sets for every query this session (the export API answered with
  0 entries even for known titles -- treated as a tool failure, not as evidence); WebSearch on
  "Stirling numbers of the first kind" + "complete Bell polynomials" + "harmonic" returned only identity
  papers (Connon 2010, Kim-Kim degenerate Stirling 2018, Integers 15 (2015) #A8 "Representations of
  Stirling numbers of the first ..." -- none with inequalities).
* Mod-Poisson framing (Kabluchko, Marynych, Pitters 2022/23, arXiv:2209.06808, FULL TEXT GREP, eq. (8)):
  E e^{z omega_n} / e^{(log n)(e^z - 1)} -> 1/Gamma(e^z) "locally uniformly in z in C. Moreover, the speed of
  convergence in (8) is O(1/n), again locally uniformly" -- this is Hwang's (6) in probabilistic clothing
  (Feray-Meliot-Nikeghbali, Example 2.1.3), again with an O(1/n) and no constant. Their Theorems 2.1-2.3 are
  for tilted (theta-weighted) Stirling distributions; local limit theorems via [27, Theorem 2.7]; nothing
  explicit.

---

## WHAT THIS MEANS FOR THE PROGRAMME

1. **Uniform-in-m explicit-error expansion for m <= C log n: does not exist.** Hwang 1995 is the reference of
   record and is O-only. The single explicit remainder (Moser-Wyman via Adell (2)) is (a) second-hand and
   (b) three orders too weak in m for the ratio (needs m^3 < ~H/160).
2. **The gap in the few-cycle regime is huge (zeta(2)/(log n)^2 against a needed 1/n)**; the whole
   difficulty is an explicit, m-uniform remainder. Two routes are visible and neither is in print:
   (i) make Hwang's eq. (6) effective -- i.e. an explicit bound for the singularity-analysis / Gamma-ratio
   remainder in Gamma(n+w)/(Gamma(w) Gamma(n+1)) = n^{w-1}(1 + sum pi_k(w) n^{-k} + R) uniformly on |w| <= eta
   (lead, NOT verified: DLMF 5.11(iii) gives the expansion of Gamma(z+a)/Gamma(z+b) without a remainder
   constant; explicit remainder bounds for that ratio exist in the literature on Tricomi-Erdelyi expansions
   (Frenzen 1987 is the usual citation) -- to be checked before use); (ii) work directly with the exact
   Bell-polynomial form and bound the tail of Y_k by the explicit inequalities H^{(r)} <= zeta(r), which is
   the programme's own "three certificates" style.
3. **Prior art to cite in the paper:** Jordan (first-order), Moser-Wyman 1958 (small-m series in partial
   zeta sums, explicit remainder per Adell), Wilf 1993 (fixed k, full series), Hwang 1995 (m <= eta log n,
   O-terms), Temme 1993 and Louchard 2010 (other regimes), Arratia-DeSalvo 2017 and Adell 2022 (the only
   explicit bounds, other end / upper only), Sibuya 1988 (3.4) (the conjecture) and (3.6) (the ratio bound),
   Erdos 1953 (unique mode), Liang-Sagan Conjecture 8.9 (n-direction open).
4. **Negative results recorded so nobody searches again:** no Turan/Newton gap from harmonic-number
   polynomials; no sharper ratio bound than Sibuya (3.6); no explicit lower bound for c(n,m) at small m;
   n-direction log-concavity/convexity is conjectural; arXiv export API unusable this session.

## ACQUISITION RECORD (all 2026-09-03, scratchpad
`C:\Users\user\AppData\Local\Temp\claude\C--Users-user-ScienceBro\ac66e2dc-eaec-44ab-b12d-0b61f841fa72\scratchpad`)
| file | source | status | sha256 (prefix) |
|---|---|---|---|
| hwang1995_preprint.pdf/.txt/.png | CiteSeerX (LIX preprint, Sept 1994) | FULL TEXT | a604153e |
| wilf1993.ps/.pdf/.txt | math.upenn.edu/~wilf/website/stirling.ps | FULL TEXT | 7fd9c0de |
| adell2022.pdf/.txt, adell_p-02..04.png | zaguan.unizar.es/record/119973 (CC-BY) | FULL TEXT | 8f961785 |
| arratia_desalvo.pdf/.txt | arXiv 1404.3007v4 | FULL TEXT | ab6b2921 |
| louchard2010.pdf/.txt | dmtcs.episciences.org/501/pdf | FULL TEXT | a183d963 |
| temme1993.pdf/.txt | ir.cwi.nl/pub/2304/2304D.pdf (scan) | FULL TEXT (OCR garbled in formulas) | 95351182 |
| erdos1953.pdf/.txt | users.renyi.hu/~p_erdos/1953-04.pdf | FULL TEXT (OCR garbled in formulas) | cf4dc720 |
| lc_lattices.pdf/.txt | arXiv 2408.02782v1 (Liang-Sagan) | FULL TEXT | be5834f9 |
| kab2022_modphi.pdf/.txt | arXiv 2209.06808v2 | FULL TEXT GREP | 9119a85b |
| clt_stirling.pdf/.txt | arXiv 2208.09928v1 (Heim-Neuhauser) | FULL TEXT GREP | d02d1488 |
| pomerance_maxstir.pdf/.txt, cp_p-01.png | math.dartmouth.edu/~carlp/MaxStirCorr.pdf | FULL TEXT (intro) | fbd940af |
| hz2020.pdf/.txt | arXiv 2008.10069 (already present) | FULL TEXT GREP | 9643f516 |
Not retrieved (paywalled, 403): Moser-Wyman 1958 (OUP, Wiley), Hammersley 1951, Hwang 1995 journal
version (ScienceDirect), Chelluri-Richmond-Temme 2000, Adamchik 1997, Liang-Sagan published version.
Abstract/metadata only: Connon 2010, Adell-Lekuona 2019, Frenzen 1987 (lead only).
