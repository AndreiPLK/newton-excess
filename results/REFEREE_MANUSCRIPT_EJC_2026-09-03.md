# Referee report: `projects/qg-bootstrap/release/manuscript/main.tex` (Electronic Journal of Combinatorics)

**Referee role:** external referee, reading as an EJC editor would, checking the manuscript text
against the package it describes.
**Report written:** 2026-09-03, clock checked with `date` at 19:46:38 RDT (session start) and
19:56:01 RDT (report time).
**Scope reviewed:** `release/manuscript/main.tex` (9 pages, compiles to `main.pdf`), cross-checked
against `release/ARTICLE.md`, `release/README.md`, `release/LIMITATIONS.md`, `release/AI_DISCLOSURE.md`,
`release/OUTREACH.md` §4.2, every script named in the manuscript under `release/scripts/`, and the
logs/reports in `results/` (`theorem_full_2026-09-03.txt`, `sibuya_theorem_full_recheck_2026-09-03.txt`,
`sibuya_corner_grid802_repaired_2026-09-03.txt`, `DEBATE_*`, `VALIDATION_*`,
`REFEREE_REVIEW_RELEASE_2026-09-03.md`). Citations verified live against Crossref and the arXiv API
(internet access confirmed available in this session; results below).
**Runs executed by me:** `release/scripts/variance_limit.py`, `release/scripts/certify2.py` (both to
reproduce numbers quoted in the manuscript's Proposition proof; both exit 0, both reproduced exactly).
No certificate script (`theorem.py --full`, `sibuya_theorem.py --full`, `sibuya_corner_grid.py`) was
re-run; their logged verdicts were read.
**Files modified:** this report only.

---

## VERDICT: NEEDS THESE FIXES FIRST — not a clean "ready to submit"

The mathematics and the numerical bookkeeping are in genuinely strong shape: every headline number I
checked (about 40 distinct quantities, several independently recomputed from scratch) traced correctly
to a script or a logged run, including two I reproduced myself from a cold start
(`variance_limit.py`: 359637 boxes, tail `6.79e-36`; `certify2.py`: 2163485 boxes). The package-level
blockers found by the 2026-09-03 12:03 release review (`REFEREE_REVIEW_RELEASE_2026-09-03.md`) appear
to have been substantially closed since — the dense-certificate-B validation is now full PASS
(`independently-validated`), the residual region is now stated as the single consistent `n-1-j >= 803`
window, and `LIMITATIONS.md` is current.

But I found one clear, previously unflagged **numerical error in the manuscript's own prose** (the
empirical shift formula for Piece A, item 1 below — off by a factor of exactly 2 against the logged
data), one **inconsistency in the "how it was checked" bookkeeping** relative to `README.md` (three
defects claimed vs. four actually found and fixed), a **structural defect** an EJC referee will
notice immediately (five of eight bibliography entries are never `\cite{}`d in the running text), a
genuine **notation collision** (`H` denotes two different load-bearing objects in the two theorems),
and several places where the manuscript satisfies the letter of EJC's "enough detail for a human to
check it" policy only by pointing at code, not by giving the derivation. None of this touches the
correctness of the theorems as stated; all of it is fixable without new computation.

---

## Findings, most serious first

### 1. FAIL — Piece A's stated empirical shift formula is wrong by a factor of 2 against the logged data

`main.tex:127-131` (Piece A, the ladder):

> "Empirically $c(t) = t+5$ for odd $t$ and $c(t) = t+2$ for even $t\ge2$, with $c(0)=3$"

I checked this against the actual computation. `release/scripts/ladder_fast.py:99` computes the
predicted shift as
```
pred = (i + 5) // 2 if i % 2 else (i + 2) // 2
```
and the logged shifts in `results/m_ladder_log_2026-09-02_full.txt` and
`results/theorem_full_2026-09-03.txt` confirm this, not the manuscript's formula:

| $t$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| logged $c(t)$ | 3 | 2 | 4 | 3 | 5 | 4 | 6 | 5 | 7 | 6 | 8 | 7 |
| manuscript's $t+5$ / $t+2$ | 6 | 4 | 8 | 6 | 10 | 8 | 12 | 10 | 14 | 12 | 16 | 14 |
| actual formula $(t+5)/2$ / $(t+2)/2$ | 3 | 2 | 4 | 3 | 5 | 4 | 6 | 5 | 7 | 6 | 8 | 7 |

The manuscript is missing the "`/2`" (equivalently, floor division) in both branches. This is a
one-character-class fix (`c(t) = \lceil (t+5)/2 \rceil$ for odd $t$, `c(t) = (t+2)/2$ for even $t$),
and the text is careful to say the certificate itself does not depend on the pattern being right ("the
certificate does not rely on this pattern, only on the non-negativity of the shifted coefficients") —
so no certificate is affected — but the stated closed form is simply arithmetically false as printed,
and a referee who checks it against the reproducing instructions will find the discrepancy in under a
minute. `c(0) = 3` I could not independently confirm either way (the ladder log starts at $t=1$;
$t=0$ is the trivial index $p_0 \equiv 1$ and may not need a certificate at all — if so, drop the
clause rather than leave an unverifiable value).

### 2. FAIL — "three defects" undercounts the package's own count of four

`main.tex:216` ("How the computations were checked"):

> "three defects were found this way (a missing term in a $1/M_0$ bound, a non-analytic base object, a
> non-$\sigma$-free tail), fixed, and the scripts re-run."

`release/README.md:82-83` (same package, same event) says:

> "**Four** defects found there (a missing term in a `1/M_0` bound, a non-analytic base object, a
> non-`σ`-free tail, **a multiplier off by a factor 2**) were fixed and the scripts re-run."

The fourth defect is real and traceable: `results/VALIDATION_SIBUYA_SPARSE_2026-09-03.md:60,167,201`
— "the multiplier formula has one factor-2 slip" in `sibuya_sparse_certificate.py`'s $j=2$ multiplier
derivation (`-(\log(1+u))''$ term), corrected in that report's item 8 and consistent with
`LIMITATIONS.md`'s "eight items PASS, one multiplier corrected." The manuscript's list of three is
simply short one item relative to its own package's README. Either add the fourth defect to the
manuscript's list or explain why it is excluded (e.g. if it is judged not "independent-validator"
class) — as written the two public-facing documents of the same release disagree on a factual count.

### 3. MAJOR — five of eight bibliography entries are never cited in the running text

`\cite{...}` appears exactly three times in `main.tex` (lines 67, 72, 214), citing `Sibuya1988`,
`CHR2024`, and `Arb2017`. The bibliography (`main.tex:233-251`) lists eight entries. The other five —
`Newton1707`, `Lieb1968`, `Hwang1995`, `LiangSagan2024`, `BrandenHuh2020` — are never referenced by a
`\cite{}` command anywhere in the body. (Newton is *named* in prose at line 39, "Newton's inequality
(1707)", but not tied to the bibliography entry by a citation command; the LaTeX will compile without
complaint, since `thebibliography` does not require every entry to be cited, but a referee reading the
PDF will find four names — Lieb, Hwang, Liang–Sagan, Brändén–Huh — that never appear in a single
sentence of the paper.) I traced all five against Crossref/arXiv (details in Finding 6) and they are
genuine, correctly-transcribed papers on adjacent topics, but their relevance to *this* paper is never
stated. This is very likely leftover from `ARTICLE.md §9` ("Sources"), which does discuss Lieb and
Hwang in one sentence each, and from an earlier draft's related-work paragraph that did not survive
into the LaTeX. **Fix:** either add one sentence per reference tying it to the text (Lieb 1968 and
Hwang 1995 belong naturally in the paragraph after eq. (3.4) discussing what is and is not known about
Stirling-number log-concavity; Liang–Sagan and Brändén–Huh belong in a related-work remark about
quantitative/Lorentzian log-concavity), or drop the entries that are not going to be discussed.

### 4. MAJOR — notation collision: `H` denotes two different load-bearing objects

- `main.tex:157-164` (Section 4, Theorem A): `H(\theta)` is *the limit function* of the Newton excess,
  with `H(0) = 4/5`, Taylor coefficients `h_0, h_1, h_2, h_3`, and Proposition 1 is the statement
  `H(\theta) > 4/5`.
- `main.tex:195-198` (Section 5, Theorem B): `H = H_N` is *the harmonic number* `H_N = \sum_{k=1}^N 1/k`,
  the variable the entire "$H$-model" (the load-bearing machine for pieces D′, D″) is built around.

Both are central, both are used repeatedly over several pages, and both are given the bare capital
letter `H` with no subscript convention distinguishing them (the harmonic number is sometimes written
`H_N`, sometimes bare `H`; the limit function is always bare `H(\theta)` or bare `H`). A referee
reading straight through, or one who has the two sections open side by side while checking the
`H`-model section against Section 4's `H(\theta)`, will momentarily read the wrong object. Cheap fix:
rename one of them (e.g. keep `H(\theta)` for the limit function and use `\mathcal{H}_N` or `h_N` for
the harmonic number, which is in any case the more standard symbol clash to avoid since `H` for
"harmonic number" is common in the literature and worth keeping — so better to rename the limit
function, e.g. to `\Phi(\theta)` or `\Lambda(\theta)`).

### 5. MINOR (carried over, confirm still applies) — eq. (3.4) is presented as Sibuya's own displayed inequality; it is this project's translation of his monotonicity statement

`main.tex:67`: "The inequality in Theorem 2 is equation (3.4) of Sibuya, stated there as a conjecture."
This exact point was raised in the prior package-level review (`REFEREE_REVIEW_RELEASE_2026-09-03.md`,
finding m1) against `ARTICLE.md`/`README.md`, and it applies identically to the manuscript, which was
apparently not touched on this point. Per `results/LITERATURE_VERDICT_STIRLING_NEWTON_GAP_2026-09-02.md`
(verbatim quotation of Sibuya 1988, p. 699), equation (3.4) in the original paper is the assertion that
a specific sequence `s_m` is *strictly decreasing*, not the displayed ratio inequality
`p_j^2/(p_{j-1}p_{j+1}) \ge 1+1/(3n-j)`. The translation is algebraically verified sound (both by the
prior review and independently by `VALIDATION_SIBUYA_DENSE_B_2026-09-03.md:50`, an OCR re-derivation),
so this is not a mathematical problem — but stating flatly "the inequality ... is equation (3.4)" invites
a referee who pulls the 1988 paper to conclude the displayed formula was misquoted. One clause fixes it:
"...is the equivalent ratio form of equation (3.4) of Sibuya (a monotonicity statement in his notation;
the equivalence is elementary and stated in the supplementary material)."

### 6. Citation verification (live, against Crossref / arXiv API)

All entries **PASS** on title/author/venue/volume/pages/year, with two notes:

| key | manuscript claim | verified | verdict |
|---|---|---|---|
| `Sibuya1988` | Sibuya, *Log-concavity of Stirling numbers...*, Ann. Inst. Statist. Math. 40 (1988), 693–714 | Crossref: exact title/author/journal/volume/pages match, DOI `10.1007/bf00049427` | PASS |
| `CHR2024` | Cheung, Hillman, Remmen, *Bootstrap Principle...*, Phys. Rev. Lett. 133 (2024), 251601; arXiv:2406.02665 | Crossref + arXiv API: exact title/authors/journal/volume/article-number/year, DOI `10.1103/PhysRevLett.133.251601` | PASS |
| `Lieb1968` | Lieb, *Concavity properties...*, J. Combin. Theory 5 (1968), 203–206 | Crossref: exact match, DOI `10.1016/s0021-9800(68)80057-2` | PASS (uncited — Finding 3) |
| `Hwang1995` | Hwang, *Asymptotic expansions...*, J. Combin. Theory Ser. A 71 (1995), 343–351 | Crossref: exact match, DOI `10.1016/0097-3165(95)90010-1` | PASS (uncited — Finding 3) |
| `LiangSagan2024` | Liang, Sagan, arXiv:2408.02782 (2024); Conjecture 8.9 | arXiv API: exact title/authors/date. I downloaded the PDF and located Conjecture 8.9 (p. 33): "Given $k$, there is an integer $N_k$ such that $(c(n,k))_{n\ge0}$ is log-concave for $n<N_k$ and log-convex for $n\ge N_k$." This exists and is correctly numbered, but it concerns a log-concave-to-log-convex *transition point*, not a quantitative margin — a different phenomenon from Sibuya's inequality. It is never discussed in the running text (Finding 3), so its relevance is asserted only by co-location in the bibliography. | PASS on transcription, uncited (Finding 3), relevance unstated |
| `BrandenHuh2020` | Brändén, Huh, *Lorentzian polynomials*, Ann. of Math. 192 (2020), 821–891 | Crossref: exact title/authors/journal/volume/year, DOI `10.4007/annals.2020.192.3.4`; Crossref did not return a page range in its metadata (issue 3 confirmed, consistent with the stated pages) | PASS, page range not independently re-confirmed via API (uncited — Finding 3) |
| `Arb2017` | Johansson, *Arb: efficient arbitrary-precision midpoint-radius interval arithmetic*, IEEE Trans. Comput. 66 (2017), 1281–1292 | Crossref: exact match, DOI `10.1109/tc.2017.2690633` | PASS |
| `Newton1707` | Newton, *Arithmetica Universalis*, 1707 | Pre-DOI historical text; standard citation form, nothing to verify against an API | PASS (as far as checkable), uncited (Finding 3) |

No fabricated or unverifiable citation was found — a real positive, and directly relevant to EJC's AI
policy warning about "fictitious citations" (Finding 7).

### 7. EJC's AI policy: "enough detail ... to allow another human to check it"

The exact clause, quoted from `release/OUTREACH.md §4.2` (EJC About page, verbatim in full there):

> "If you use AI to assist in finding a proof or an argument, you must check it yourself and provide
> enough detail in your paper to allow another human to check it. You must also ensure that prior work
> is properly referenced, since AI often takes ideas from previous publications without attribution and
> sometimes suggests fictitious citations."

**Where the manuscript satisfies this well:** the AI-disclosure paragraph is present and honest
(`main.tex:229-231`); Sibuya's priority is credited correctly and explicitly (Finding 6 found no
fictitious citation — a real point in the manuscript's favour against the specific risk the policy
names); the sharpness argument (§ "Sharpness", the $X(n)$ cubic and its shift) and Proposition 1's proof
sketch are given in enough closed form that I reproduced them independently without running any script
(hand/`python -c` verification of $M(n,1)$ at $n=5,11,21,41,81,161$, and of $X(n+4)$'s four positive
coefficients, both confirmed exactly).

**Where a referee would honestly say "too thin, points to code":**

- **Piece C1 (sparse certificate), `main.tex:137`.** "the remaining terms are bounded on the region by
  explicit lemmas (a tail lemma for the $e_j(\beta)$, a moment lemma for the deviations $\beta_k$, and a
  bound $C_i \le i^2 2^i$ on the coefficient growth beyond $i=60$, all stated with their constants **in
  the script**...)" — the lemmas are named, not given. A reader cannot check the $0.1417$ residual
  without opening `sparse_certificate_full.py`.
- **Pieces C2, C3, C2′, C3′ (dense certificates), `main.tex:139-143,192`.** The Edgeworth remainder
  bound, the "explicit constants," and the box-adaptive-subdivision rule are described only by name
  ("the Edgeworth remainder and the cumulant tail are enclosed by explicit constants"); none of the
  actual bounding inequalities appear in the paper.
- **The $H$-model tail defect, `main.tex:202`.** "$|r_k| \le t^2 R_k(H)$, where $R_k$ is an explicit
  polynomial in $H$ obtained from one convolution" — $R_k$ itself is never written down, even though
  this is precisely the object the 3 September debate found and fixed a factor-2.8 error in
  (`results/DEBATE_REPAIR_TAIL_2026-09-03.md`). Given that this exact term was the site of a real bug
  five hours before the manuscript's own writing, a referee would specifically want to see it in the
  paper, not just hear that it exists.
- **The literature search for the open region, `main.tex:210`.** "A literature search ... returned
  nothing usable; the verdict, with the nearest results and why they do not apply, is recorded in the
  repository." This is a pointer, not content — the nearest results and the reason they fail are
  exactly the kind of detail the AI policy asks to be *in the paper*.

**Recommended minimal additions** (cheapest fixes, roughly in order of value): (a) state the
$C_i \le i^2 2^i$ bound's one-line derivation or its source lemma name in the text, not only "in the
script"; (b) write out $R_k(H)$'s degree and leading behaviour (even without the full polynomial,
stating its degree and that it comes from "one convolution of the geometric tail with itself" — one
sentence) since this is the exact object a June-Sept 2026 defect lived in; (c) name the one or two
nearest literature results for the open region (the repository's own verdict file already has this —
Korolev–Zhukov 2000 / Giuliano–Weber 2017 / Siripraparat–Neammanee 2021 class results are named
elsewhere in this project's history as the negative verdict) instead of only pointing at "the
repository." None of these require new computation — they are transcription from files that already
exist.

### 8. Theorem 2 statement and index conventions — PASS, checked line by line

- $e_j(1,\dots,n-1) = c(n,n-j)$ (`main.tex:60`): the standard rising-factorial identity
  $\prod_{k=1}^{n-1}(x+k) = \sum_j c(n,n-j)\,x^j$. Correct.
- $1+1/(3n-j)$ (`main.tex:62`, Theorem 2) vs. the script comment's $1+1/(2N+3+j')$
  (`sibuya_theorem.py`'s docstring: "Sibuya's target `1 + 1/(3n-j) = 1 + 1/(2N+3+j')`") and the
  manuscript's own eq. (sibuyaE) (`main.tex:187`): with $N=n-1$ and $j'=N-j$,
  $3n-j = 3N+3-j = 2N+3+(N-j) = 2N+3+j'$ — an exact algebraic identity, confirmed by direct
  substitution. **Consistent.**
- The reciprocal-duality step $e_j(1,\dots,N) = N!\,e_{N-j}(1,\tfrac12,\dots,\tfrac1N)$
  (`main.tex:185`) matches `sibuya_harmonic.py`'s docstring verbatim ("Since
  `e_j(1..N) = N! E_{N-j}`"). **Consistent.**
- $T$ (`main.tex:187`) matches `sibuya_harmonic.py`'s $T(N,j')$ formula term for term, including the
  $N\to\infty$ limit $T\to(j'+1)/j'$ stated in both. **Consistent.**
- One notation gap, not an error: `main.tex:180` (start of §5, "Theorem 2: the same machine on
  $\{1,\dots,N\}$") uses $N$ immediately without redefining it in that section; the reader must recall
  $N=n-1$ from the general setup in §1 (`main.tex:39`), four pages and one full theorem earlier, and
  Theorem 2's own statement (`main.tex:59-65`) never uses the letter $N$ at all. A one-clause reminder
  ("with $N=n-1$" at the start of §5) would remove the need to backtrack.

### 9. Numeric cross-check — PASS on every value checked (about 40 quantities; failures noted above)

Every number I could locate a source for matched the manuscript exactly, including several I
recomputed independently from a cold start rather than trusting a log file:

| manuscript claim | location | check |
|---|---|---|
| $M(n,1)$ at $n=5,11,21,41,81,161$: $1.3559,1.0376,0.9158,0.8570,0.8283,0.8141$ | Sharpness, `main.tex:177` | recomputed independently (exact `Fraction` arithmetic on the closed-form $\rho$): matches to the last printed digit at every point |
| $X(n+4)$ has four positive coefficients | Sharpness, `main.tex:177` | recomputed: coefficients $56/45,\,100/9,\,1264/45,\,640/45$, all positive |
| $(0,0.2]$: $(352/945)v^6$, tail $\le 6.8\times10^{-36}$; $[0.2,3.5]$: $359{,}637$ boxes; $A(3.5)=0.36928$ | Prop. 1 proof, `main.tex:171` | reran `variance_limit.py` cold: `359637 boxes`, tail `6.79e-36`, `A(3.5)=0.36928476` — exact match |
| second proof, $2{,}163{,}485$ boxes | `main.tex:171` | reran `certify2.py` cold: `boxes 2163485` on `v in [0.15, 2.331122]` — exact match |
| $[a^3]/b^2=176/175$, residual $0.1417$ (piece C1) | `main.tex:137` | `results/sparse_certificate_full_2026-09-03.txt`: `176/175`, `0.141701` |
| 930 boxes, $0.800004$ (piece C2) | `main.tex:140` | `results/dense_certificate_a_rerun_2026-09-03.txt`: `930 boxes ok, 213 skipped`, `0.800004` |
| $S(0)=176/525$, 18 bands, $\min\ge0.3206$ (piece C3) | `main.tex:143` | `results/dense_certificate_b_rerun3_2026-09-03.txt`: `176/525`, 18 band lines, first band `+0.3206` |
| 627 rungs, 7.06 h | `main.tex:131` | `results/m_ladder_log_2026-09-02_full.txt`: `627 5026 316 n<=1259`, `total CPU time 7.06 h` |
| $[a^3]/b^2=2/3$, residual $0.29$ (piece C1′) | `main.tex:192` | `results/sibuya_sparse_certificate_2026-09-03.txt`: `2/3`, `0.290647` |
| $756+84$ boxes, margins $5\times10^{-6}$ and $0.163$ (piece C2′) | `main.tex:192` | `sibuya_dense_a_2026-09-03.txt`: `756 boxes ok`, margin `0.000005`; `sibuya_dense_a_top_2026-09-03.txt`: `84 boxes ok`, margin `0.162638` — note this corrects a mistranscription the 2026-09-03 12:03 package review found in `ARTICLE.md` (which had `1.0e-5`, the validator's number, not the script's `5e-6`); the manuscript now has the right number |
| $S(0,\zeta)=1/3$, 18 bands, $\min\ge0.2398$ (piece C3′) | `main.tex:192` | `sibuya_dense_b_2026-09-03.txt`: `1/3`, 18 bands, first band `+0.2398` |
| 4705 grid steps, worst relative margin $1.92\times10^{-10}$, $N\approx4.85\times10^8$; 770 indices | `main.tex:195` | `sibuya_corner_grid802_repaired_2026-09-03.txt`: `4705 grid steps ... worst relative margin 1.924e-10`; `770/770 indices certified` |
| 280,756-grid coverage check, C1 reaches $b\le(628/629)^2$, C3 reaches $b\approx0.926$ | `main.tex:146` | `theorem.py:22` docstring and `results/DEBATE_A_VERIFIER_2026-09-03.md:198` (280,756-point grid, 0 gaps) — match |
| the 2.8× Poisson-ratio defect, $1/N\le e^{\gamma-H}$ substituted in the unsafe direction | `main.tex:216` | `results/DEBATE_REPAIR_TAIL_2026-09-03.md`: ratio `1.697` at $H=20$, `2.793` at $H=12$ ("false by up to a factor of 2.8"); "new bound is about 2.8x TIGHTER" — match |
| the loss-vs-threshold fix in Theorem A's debate | `main.tex:216` | `results/DEBATE_A_SCEPTIC_2026-09-03.md` (the `loss_corner < 0.1` "magic threshold" finding) and `release/scripts/dense_certificate_b.py:600-605` (now `corner_ok = (S_corner - loss_corner).lower() > 0`) — confirmed fixed in the shipped script, not just recommended |
| "seven blockers, all closed" | `main.tex:216` | `results/REFEREE_REVIEW_RELEASE_2026-09-03.md` names exactly seven blockers B1–B7; spot-checked B1 (dense-certificate-B validation) closed via the report's own Addendum 2 ("Recommended state ... `independently-validated`") |
| two adversarial reviews per theorem (sceptic + verifier) | `main.tex:216` | `DEBATE_A_SCEPTIC`, `DEBATE_A_VERIFIER`, `DEBATE_SCEPTIC`, `DEBATE_VERIFIER` (4 files = 2 per theorem) all present and dated 2026-09-03 |
| DOI `10.5281/zenodo.22282840` | `main.tex:220` | matches `release/README.md`, `release/ARTICLE.md` throughout |

**No number I checked was wrong except the shift formula in Finding 1.** I did not independently
re-derive the Edgeworth-remainder machinery inside `dense_certificate_a/b.py` or
`sibuya_dense_a/b.py` (I read the logs and the `VALIDATION_*` reports, as the prior package review
also did) — that is the one class of claim in the table above I am trusting rather than re-deriving.

### 10. What is proved vs. what is claimed — no gate violation found in the manuscript text

Cross-checked every use of "proved"/"proves" (7 instances) against `ARTICLE.md §6-7` and the current
`LIMITATIONS.md`: all are correctly scoped (Theorem A: unconditionally "complete," now earned — see
below; Theorem B: hedged to the three stated ranges, with §"What is not proved" giving the residual
region exactly, matching `LIMITATIONS.md`'s current Theorem B section). None of the forbidden words
"discovered / novel / first / confirmed / refuted" appears anywhere in the manuscript. The claim
"Theorem A: complete... every certificate is independently validated" (implicit in "five machine-checked
pieces and is complete", `main.tex:33`) was a blocking claim-gate violation in the 12:03 package review
(B1: the dense-certificate-B validator had recommended `experimentally-supported`, not
`independently-validated`) — I checked this specifically and found it **resolved**:
`results/VALIDATION_DENSE_CERTIFICATE_B_2026-09-02.md`'s final "Addendum 2" (dated
"Thu Sep 3 ... after commit `6393d87`") now reads "All seven items now PASS... Recommended state for
`dense_certificate_b`: `independently-validated`." The manuscript's "complete" claim is now earned by
the gate as I read the current validation report. **This is a genuine, verified improvement over the
package state the 12:03 review saw**, and the manuscript is not overclaiming here.

---

## Summary table

| # | severity | location | defect |
|---|---|---|---|
| 1 | FAIL | `main.tex:127-131` | empirical shift formula for Piece A is off by a factor of 2 against the logged data (`t+5`/`t+2` should be `(t+5)/2`/`(t+2)/2`) |
| 2 | FAIL | `main.tex:216` vs `README.md:82-83` | "three defects" undercounts the package's own "four defects" (missing the sparse-Sibuya multiplier factor-2 slip) |
| 3 | MAJOR | `main.tex:233-251` | five of eight bibliography entries never `\cite{}`d in the running text |
| 4 | MAJOR | `main.tex:157-164` vs `195-198` | `H` denotes two different objects (limit function of Thm A vs. harmonic number of Thm B) |
| 5 | minor (carried over) | `main.tex:67` | states eq. (3.4) *is* the displayed inequality; it is this project's equivalent translation of Sibuya's monotonicity statement |
| 6 | PASS, one caveat | bibliography | all 8 entries verified against Crossref/arXiv; no fictitious citation; 5 uncited (see #3) |
| 7 | mixed | AI-policy compliance | disclosure and priority-crediting are strong; four specific spots ("explicit lemmas ... in the script", the Edgeworth remainder, $R_k(H)$, the literature verdict) are named but not given, and a referee would ask for at least one sentence more on each |
| 8 | PASS | Theorem 2 / index conventions | $3n-j=2N+3+j'$, the reciprocal-duality identity, and $T$ all check out exactly against the scripts; one un-reintroduced symbol ($N$ at the start of §5) |
| 9 | PASS (with #1 as the exception) | ~40 numeric claims | every other number traced to a log or was independently reproduced from scratch |
| 10 | PASS | claim-gate compliance | no forbidden word misused; "Theorem A: complete" is now earned (dense-certificate-B validation closed to full PASS since the 12:03 package review) |

---

## What is genuinely strong and should not be lost in revision

The theorem statements, the instrument (two exact representations of the excess), the piece-by-piece
decomposition, and the coverage argument are mathematically sound and match the underlying computation
with unusual fidelity for a 9-page paper describing this much machine work — I was able to reproduce
two of the harder numbers (`359,637` and `2,163,485` boxes) from a cold `uv run` rather than trusting
the manuscript, and both matched exactly. The claim-gate discipline (no "proved" without an earned
gate, the residual region stated honestly, the debate-found defects disclosed rather than hidden) is
real and is exactly what an EJC referee wants to see from an AI-assisted submission. The two failing
findings (1, 2) are both small, mechanical, and fixable by editing the manuscript text against numbers
that already exist in the repository — no new computation is required for any of the ten findings
above.

**Publication remains a deliberate human action.** This report is a referee's list of what to fix
before submission; nothing here authorises submitting, and the author decides.
