# Referee review of `projects/qg-bootstrap/release/`

**Reviewer role:** release-reviewer, refereeing as an experienced combinatorialist sceptical of
computer-assisted proofs, with one hour.
**Written:** 2026-09-03 12:03:05 RDT (clock measured with `date`).
**Scope reviewed:** `release/README.md`, `ARTICLE.md`, `LIMITATIONS.md`, `AI_DISCLOSURE.md`,
`.zenodo.json`, `CITATION.cff`; the scripts `theorem.py`, `sibuya_theorem.py`, `sibuya_harmonic.py`,
`sibuya_corner_grid.py`, `sibuya_top_w.py`, `sparse_certificate_full.py`, `dense_certificate_a.py`,
`dense_certificate_b.py`, `sibuya_sparse_certificate.py`, `sibuya_dense_a.py`, `sibuya_dense_b.py`,
`fig3d_theorem.py`, `fig3d_sibuya.py`; the logs and `VALIDATION_*.md` in `results/`.
**Runs executed by the reviewer:** `theorem.py` (no `--full`) and `fig3d_sibuya.py` only, as instructed.
**Files modified by the reviewer:** this report only.

---

## Verdict, up front

**DO NOT PUBLISH AS IT STANDS.** The mathematics is in far better shape than the package around it,
but there are two blocking defects that a referee will find in the first twenty minutes — one is a
claim-gate violation (a certificate is called "independently validated" when its own validator
explicitly refused to promote it), and one is a wrong number printed inside a published figure. There
is also an unresolvable-as-written inconsistency about what Theorem B actually covers: README, ARTICLE
§6 and ARTICLE §7 state three different residual regions, none of which equals the region the scripts
actually leave open.

Seven blockers, six majors, nine minors below. The headline numbers themselves are, with two
exceptions, correctly transcribed from logs that exist — the arithmetic is not the problem, the
bookkeeping around it is.

---

# BLOCKERS (must change before any publication)

## B1. `dense_certificate_b.py` is called "independently validated"; its validator refused to validate it

This is a claim-gate violation, and it is the one that would end the referee's goodwill.

`results/VALIDATION_DENSE_CERTIFICATE_B_2026-09-02.md`, final line of the re-check (line 235), reads
verbatim:

> Item 4 remains a text-level FAIL with a bounded, harmless residual (`<= 1.5e-5` at `V = 0.17`, far
> less elsewhere): (i) exact-degree tracking through `divV` in `evaluate()`, (ii) either the
> polynomial-`G_j` justification or an explicit `D` term for the dropped monomials, (iii) one written
> sentence propagating the `i > 60` remainder. **Recommended state unchanged: `experimentally-supported`;
> promote to `independently-validated` after (i)-(iii) and a rerun** (item 4 only; no other item needs
> re-validation).

The same report, line 79, records `### Item 4 -- the Ser sup-bound bookkeeping: FAIL as written`.

Against that, the package asserts the opposite in five places:

- `release/ARTICLE.md:119-120` — "**Theorem A: complete.** Every piece is an artifact checked by a
  script, and **every certificate is independently validated**."
- `release/README.md:78-79` — "Each certificate was re-derived and re-run by an independent agent that
  did not import the code it validated".
- `release/.zenodo.json` (description) — "each was validated by an independent re-implementation."
- `release/CITATION.cff` (abstract) — "each was validated by an independent re-implementation."
- `release/scripts/theorem.py:141` — printed output: "all three independently validated."

Under `.claude/rules/claim-gates.md` the word **proved** requires the matching gate; Theorem A is
declared *complete* and *proved* on the strength of a certificate the validator placed at
`experimentally-supported`. **Either close items (i)-(iii) and rerun, or downgrade the wording
everywhere to "two of three certificates independently validated; the third validated on its numbers
with a named text-level bookkeeping gap (residual <= 1.5e-5)".** The residual really does look
harmless, which is exactly why the honest sentence is cheap to write and the overclaim is inexcusable.

## B2. The Sibuya figure prints the wrong constant

`release/scripts/fig3d_sibuya.py:118` emits, and `release/data/sibuya_3d.svg` therefore contains:

    the lowest edge: t = 1, sinking towards 4/5

The Sibuya floor is **1/3**, not 4/5. The label is copy-pasted from `fig3d_theorem.py` and never
updated; every other label in the same figure says 1/3 ("the floor M = 1/3 — never touched", "the
limit 1/3 = Var/mean² of the roots"). This figure is the second image in `README.md:38` and would be
the first thing on the Zenodo landing page. Exact values I recomputed independently
(`e_j(1..N)` by the row recurrence in `fmpq`): at `n = 513`, `M(n,1) = 0.33355006...`, descending to
1/3 — confirming 1/3 and refuting the label.

## B3. The Sibuya figure draws the floor plane at 0.32, not at 1/3, and thereby inflates the margin ~60x

`release/scripts/fig3d_sibuya.py:63`:

    floor = [proj(0, 0, 0.32), proj(1, 0, 0.32), proj(1, 1, 0.32), proj(0, 1, 0.32)]

The plane is drawn at height `0.32` while its label (line 69) says `M = 1/3` (= 0.33333). The whole
point of the picture is that the yellow path approaches the floor without touching it. The true
clearance at `n = 513` is `0.33355 - 0.33333 = 2.2e-4`; the drawn clearance is
`0.33355 - 0.32 = 1.4e-2`, about **60 times larger than reality**. A sceptical reader who measures the
figure will conclude the margin is comfortable when it is razor-thin.

For contrast, `fig3d_theorem.py:61` draws its floor at exactly `0.8 = 4/5`, correctly. Only the Sibuya
figure is wrong. Fix the constant to `1/3` (and, if the gap then becomes invisible, say so in the
caption rather than moving the plane).

## B4. ARTICLE §7 states a residual region that is neither what the scripts leave open nor what README says

This is the single most important mathematical bookkeeping defect. I reconstructed the covered region
independently from the scripts' own region constants, not from the prose.

**Region constants, read from source:**

| piece | script | region as coded |
|---|---|---|
| A' | `lab/sibuya_ladder.py` | `j <= 1000`, every `n` (`sibuya_ladder_log_full.txt`: "all 1000 indices present") |
| C1' | `sibuya_sparse_certificate.py:68-69` | `TAU_MIN = 1000`, `BMAX = fmpq(1)` → `tau >= 1000`, `tau^2/N <= 1` |
| C2' | `sibuya_dense_a.py:46-47` | `THETA_LO = 1/10`, `THETA_HI = 1/2`, `t >= 1001`, `N >= 2002` |
| C2't | `sibuya_dense_a.py:51-52` (`--top`) | `THETA_LO = 1/2`, `THETA_HI = 9/10`, `t >= 1001`, `N >= 1112` |
| C3' | `sibuya_dense_b.py` | `VMAX/ZMAX/ETAMAX`; `theta <= 0.1`, `j^2/N >= 1` |
| D' | `sibuya_harmonic.py` | `j' = N - j in [1, 32]`, every `N >= 1001 + j'` |
| D'' | `sibuya_corner_grid.py` | `33 <= j' <= 399`, every `N >= 1001 + j'` |
| E' | `sibuya_dense_a.py --top2 100` | `9/10 <= theta <= 99/100`, `t >= 1001`, `j' >= 100` — **exploratory** |

**Union, worked through.** `j <= 1000` → A'. For `j >= 1001`: `theta in [0.1, 0.9]` → C2'/C2't;
`theta <= 0.1` splits on `j^2/N` into C3' / C1'. So certificates give
`{j <= 1000} ∪ {j >= 1001, theta <= 0.9}`. The corner scripts add `{j' <= 399}` (and for
`N < 1001 + j'` we have `j <= 1000`, already in A' — that seam is clean). Hence

> **True residual, counting only certificates: `j >= 1001` **and** `j/(n-1) > 0.9` **and**
> `n - 1 - j >= 400`.**

ARTICLE.md:127 instead says:

> Sibuya's inequality for `n - 1 - j >= 33` together with `j/(n-1) > 0.99`.

Both halves are wrong. `>= 33` is stale — `sibuya_corner_grid.py` closes `j' <= 399`
(`results/sibuya_corner_grid_2026-09-03.txt`: "VERDICT: ... holds for 33 <= j' <= 399 and every N ...:
True, EXIT 0"), so the correct bound is `>= 400`. And `> 0.99` silently promotes the **exploratory**
sweep E' to a certificate: ARTICLE §6 line 123 itself admits "(that last sweep is exploratory)", and
the log `results/sibuya_dense_a_top2_exploratory_2026-09-03.txt` self-labels
"EXPLORATORY (not a certificate)". §7 must not consume it. `OBSTRUCTION.md:835` agrees with me and not
with §7: "the named obstruction (02h01 (validator clock, quoted from the report it cites)): **the top regime theta > 0.9** of Sibuya's (3.4)".

**Rewrite §7 to the region above.** As written, §7 claims proof of the band
`0.9 <= theta <= 0.99, 33 <= j' < 400` that no certificate covers.

## B5. README, ARTICLE §6 and ARTICLE §7 state three mutually inconsistent versions of Theorem B

- `README.md:34-35`: "for every `j ≤ 1000` and every `n`; for every `j ≥ 1001` with `j/(n−1) ≤ 0.9`;
  and, at the top of the row, **for every `n` with at most 399 indices missing**." (No mention of the
  0.9–0.99 band at all.)
- `ARTICLE.md:26-27` (abstract): "... for every `n` with **at most 32 indices from the top**; and for
  `0.9 <= j/(n-1) <= 0.99` with `n - 1 - j >= 100`."
- `ARTICLE.md:127` (§7): residual `j' >= 33` and `theta > 0.99`.

The README version is the closest to correct but omits E' entirely; the abstract is stale on the
corner (32 vs 399) and launders E' into the theorem statement without the exploratory caveat that
appears only 96 lines later at §6. A referee reading the abstract and then §7 gets two different
theorems. Pick one statement — the B4 region — and use it verbatim in all three places.

## B6. `sibuya_top_w.py` is shipped in the release and it **fails**

`results/sibuya_top_w_2026-09-03.txt` (timestamped 11:52 today, i.e. current):

    Bernoulli remainder lemma (K=6): True  [2185 boxes]
      coverage: theta(Y = 1/35) <= 0.89833 <= 0.9: True
    FAILED box: 0.0 0.0007137758743754461 1.542172956266865e-80 3.08434591253373e-80 {'why': 'k2 not positive'}
    EXIT 1

This script's own docstring (line 1) says it is exactly the piece that would close the residual region
of B4: "the top regime in the variable Y = 1/v: theta >= 0.9, j' = N - j >= JP_MIN", with
`JP_MIN = 400  # sibuya_corner_grid.py closes 33 <= j' <= 399`. It is referenced by **nothing** in
README.md or ARTICLE.md. So the release ships a non-exit-0 script, silently, whose failure is precisely
the open problem — while §7 describes the open problem as something else (B4).

Keeping the failing run visible is correct under the experiment contract, and I credit it. But it must
be *named*: add `sibuya_top_w.py` to the scripts table with its status ("attempts the residual region;
currently fails at the `k2 not positive` box at `Y -> 0`"), or remove it from the package. As it
stands, a referee who runs everything in `scripts/` gets a non-zero exit and no explanation.

## B7. `LIMITATIONS.md` is stale by a full day and contradicts the rest of the package

`LIMITATIONS.md` was last touched 2 Sept 23:03 and still describes a pre-assembly state. It directly
contradicts README and ARTICLE:

- Line 7-9: "**The conjecture itself, for all `n`.** ... As of 2 September 2026 the certificates cover
  indices `i = 0..584`, giving `5 <= n <= 1207`" — the ladder is at 627 rungs / `n <= 1259`, and
  Theorem A is claimed complete for **all** odd `n >= 5`.
- Lines 135-137: "**The conjecture is assembled but not yet called proved.**"
- Lines 138-142: "**Two validation reports are pending**: ... `dense_certificate_a.py` ... and
  `dense_certificate_b.py`. Independent agents were started at 23:02; their first run was killed ...
  before writing." Both reports now exist (`VALIDATION_DENSE_CERTIFICATE_A/B_2026-09-02.md`).
- Lines 143-146: "**The ladder artifact.** ... the durable log ... is being regenerated (rung 186
  reached at 21:57 ...)". It completed: `m_ladder_log_2026-09-02_full.txt` ends
  "# total CPU time 7.06 h; all 627 rungs present, every shift equals the prediction, no None".
- Lines 159-162: "**What is claimed at what level.** ... the whole theorem: assembled, awaiting the two
  validations and the ladder log."

`LIMITATIONS.md` also contains **no Sibuya section at all** — Theorem B, the half of the package with a
genuinely open region, has no limitations entry. The document a referee is told (README:70) separates
"the present state from a referee's 'proved'" is the document most out of date. Rewrite it, and give
Theorem B its own section stating the B4 residual and the E'-is-exploratory caveat.

---

# MAJOR DEFECTS

## M1. ARTICLE §4 cites a `theorem.py --full` run whose log does not exist

`ARTICLE.md:96-97`:

> `theorem.py` re-runs the ladder check and the coverage lemma (1 s) and, with `--full`, the three
> certificates (**409 s, three `VERDICT True`, exit 0**).

There is no such log anywhere in the repository. I searched `results/` (39 `.txt` files) and the whole
tree: the only `*theorem*full*` artifact is `sibuya_theorem_full_2026-09-03.txt`, which is the *other*
theorem. The "409 s, three VERDICT True, exit 0" is therefore an unlogged assertion in a package whose
README (line 62) promises "Each run writes a log under `../results/`". Theorem B's equivalent claim is
properly logged; Theorem A's is not. **Produce `results/theorem_full_2026-09-03.txt` or delete the
parenthesis.**

## M2. `sparse_certificate_full.py` has no run log; its headline `0.1417` survives only second-hand

`ARTICLE.md:92` states for piece C1: "`[a^3]/b^2 = 176/175` exactly, residual `0.1417`". There is no
`results/sparse_certificate_full_*.txt`. The number is recoverable only from two indirect sources:

- `VALIDATION_SPARSE_CERTIFICATE_FULL_2026-09-02.md:9` — "sweep min `0.141702` (3 boxes), VERDICT True,
  exit 0, 2 s";
- `OBSTRUCTION.md:805` — "the 4/5 certificate certified 14% (`0.1417` out of `1.006`)".

That the validator reproduces it is reassuring, and note the same report (line 234) records "own min
`0.4375` (script's `0.1417` is a valid coarse enclosure)" — so the script is conservative, which is the
right direction. But the primary artifact is missing while its sibling certificates
(`dense_certificate_a_rerun_2026-09-03.txt`, `dense_certificate_b_rerun2_2026-09-03.txt`) both have one.
Generate the log.

## M3. Every reproduction command in README.md is wrong from the repository root

`README.md:52-57` gives:

    uv run python release/scripts/theorem.py --full
    uv run python release/scripts/sibuya_theorem.py --full
    uv run python release/scripts/sibuya_harmonic.py 33
    uv run python release/scripts/sibuya_corner_grid.py 399
    uv run python release/scripts/fig3d_theorem.py
    uv run python release/scripts/fig3d_sibuya.py

Run from the repository root (which is where `uv sync` on line 51 puts you), all six fail. Verified:

    $ uv run python release/scripts/theorem.py
    can't open file 'C:\Users\user\ScienceBro\release\scripts\theorem.py': [Errno 2] No such file or directory

The correct prefix is `projects/qg-bootstrap/release/scripts/`, which `ARTICLE.md:142-144` and the
scripts' own docstrings use correctly. Only README is wrong — and README is the file a stranger reads
first. With the corrected path both fast scripts work:

- `theorem.py` → exit 0, 2 s, `B. n = 1257: min ... = 1.807003`, `n = 1259: 1.803808`, coverage
  `khat_1 ... in [0.30274, 0.33334]`, `Coverage of C by C1 u C2 u C3: True`, ladder check OK.
- `fig3d_sibuya.py` → exit 0, "8 rows, n = 5..513; min at t=1: 0.3336", and it regenerates
  `sibuya_3d.svg` byte-identically (no git diff) — deterministic, which is good.

## M4. README.md links to `OUTREACH.md`, which does not exist

`README.md:71`, in the "What is in here" table:

    | [OUTREACH.md](OUTREACH.md) | who to send it to, and the arXiv endorsement mechanics |

`release/` contains no `OUTREACH.md`. A dead link in the package's own contents table. (Also: `data/`
holds `newton_excess_floor.svg` and `variance_is_the_constant.svg`, which no document references — either
cite them or drop them.)

## M5. "No floating-point number enters a comparison" is false in the coverage lemma of both theorems

`README.md:61` (bold) and `ARTICLE.md:31-32,147-148` assert this absolutely. The coverage lemma —
which is *part of the theorem*, not a printout — compares floats. `theorem.py:87-97`:

    imp1 = 0.05 / float(lo_all)          # theta <= 0.05 => V <= this
    imp2 = 2 * float(hi_all) ** 2 / 1.0  # b >= 1 => zeta <= this
    imp3 = 2 * float(hi_all) / 628       # t >= 628 => eta <= this
    ...
    assert imp1 <= float(DB.VMAX) and imp2 <= float(DB.ZMAX) and imp3 <= float(DB.ETAMAX)

`sibuya_theorem.py:82-91` does the same (`theta_edge`, `imp2`, `imp3`, all `float`). The inputs
`lo_all`/`hi_all` come from certified `arb` enclosures, so the *bounds* are rigorous; but the division,
squaring and final comparison are double-precision. And the margins are not generous: `theorem.py`
prints `eta <= 0.001062 <= 0.00107` — a 0.8% gap decided in floating point. This is a small fix (do it
in `fmpq`/`arb` and compare exactly) and a large credibility item, because the absolute claim is the
package's main selling point against a sceptic. Fix the code, not the sentence.

## M6. Coverage figures omit the skipped boxes, and for piece E' the omission is 61%

The box counts are quoted as if they were coverage:

- `ARTICLE.md:93` — "930 boxes, worst `N g >= 0.800004`". Log
  (`dense_certificate_a_rerun_2026-09-03.txt`): "certified: 930 boxes ok, **213 skipped**".
- `ARTICLE.md:108` — "756 + 84 boxes". Logs: "756 boxes ok, **161 skipped**"; "84 boxes ok, **6 skipped**".
- `ARTICLE.md:111` — "303 boxes, margin `1.03`". Log
  (`sibuya_dense_a_top2_exploratory_2026-09-03.txt`): "certified: 303 boxes ok, **473 skipped**",
  i.e. 61% of the sweep skipped, and the script's own line names what: "EXCEPT the corner `j' < 100` or
  `x < 20` and a sliver along `eps v = 1/20` (skipped)".

Skipped boxes are presumably handled by the subdivision/skip rules (the dense-A validator checked those,
`VALIDATION_DENSE_CERTIFICATE_A:53` "coverage, skip rules, x_lo, G bounds | PASS"), but the reader is
not told they exist. `ARTICLE.md:111` reduces the E' exclusions to "a sliver along the validity curve",
which understates "473 of 776 boxes, including the whole `x < 20` corner". State the skipped counts and
say where the skipped region is closed.

---

# MINOR DEFECTS

## m1. Sibuya's (3.4) is presented as an inequality; Sibuya wrote a monotonicity statement

Checked against `results/LITERATURE_VERDICT_STIRLING_NEWTON_GAP_2026-09-02.md:80-93`, which quotes the
paper verbatim (p. 699):

> "Numerically, it is suggested that the sequence (3.4) `((m-1)(2n+m)/(n-m+1)) [n m]/[n m-1]`,
> m = 2, 3, ..., n-1, is strictly decreasing and the same for m = n-1 and n, and that (3.5) ... is
> strictly decreasing."

The inequality `p_j^2/(p_{j-1}p_{j+1}) >= 1 + 1/(3n-j)` is **this project's translation** of that
monotonicity statement — the verdict file says so explicitly ("Translation of (3.4) (this note,
verified exactly for n = 4..12 ...)"). I re-derived the equivalence and it is correct: with `m = n-j`,
`s_m > s_{m+1}` combined with `C(N,j-1)C(N,j+1)/C(N,j)^2 = j(N-j)/((j+1)(N-j+1))` gives exactly the
factor `1 + 1/(2n+m) = 1 + 1/(3n-j)`; `VALIDATION_SIBUYA_DENSE_B_2026-09-03.md:50` independently
re-derives it from an OCR of the original and verifies it as an `fmpq` identity at 50 random `(n,m)`.
So the transcription is **sound** — but `ARTICLE.md:24-25` and `README.md:27-32` present the inequality
as Sibuya's own wording. One clause fixes it: "in the equivalent form (translation verified here; Sibuya
states the monotonicity of the sequence (3.4))". A combinatorialist who pulls the 1988 paper will not
find the displayed formula in it and will wonder what else was paraphrased.

**PASS on substance:** the range `j = 1, ..., n-2`, the equality at `j = 1`, and the credit are all
correct. `>=` with equality at `j = 1` matches Sibuya's "and the same for m = n-1 and n".

## m2. The chain `n(R-1) > n g > n/(3n-j)` — PASS

Verified as written. `g = -Delta^2 log p_j` gives `log R_j = 2 log p_j - log p_{j-1} - log p_{j+1} = g`,
so `R_j = e^g` exactly and `R_j - 1 = e^g - 1 > g` for `g > 0`. The certificates prove `N g > 1/(3 + 3/N - theta)`,
and `N/(3n-j) = N/(3N+3-j) = 1/(3 + 3/N - j/N)` is an exact identity. Chain holds.
`ARTICLE.md:66` states it correctly; `VALIDATION_SIBUYA_SPARSE:52` and `VALIDATION_SIBUYA_DENSE_A:53`
both re-derive it and mark PASS. The one thing worth adding to §3: the step needs `g > 0`, which is
Newton — say so, it costs four words.

## m3. Priority and novelty — PASS, with one word to soften

The novelty claims are *under*-stated relative to what the verdicts support, which is the right
direction. `LITERATURE_VERDICT_STIRLING_NEWTON_GAP:33-36,49` supports "conjectured in print by Sibuya
(1988), eq. (3.4) ... Sibuya did not prove it" and "the CONJECTURE has priority Sibuya 1988 and must be
cited as such". Sibuya is credited in every place the conjecture appears: `README.md:27,28`,
`ARTICLE.md:23-24,49,152-153`, `.zenodo.json` title and description, `CITATION.cff` title and abstract,
and the keyword list. Good.

`ARTICLE.md:24-25` — "unproved and, as far as our literature search reaches, untouched since" — matches
verdict item 3 ("No later paper proves, restates, or even mentions Sibuya's (3.4)"; 33 OpenAlex citers
enumerated). Fine as hedged.

The one word: `README.md:6` says the repository "**proves** that for two classical families the
inequality holds with a definite margin". For Theorem B that is false as an unqualified sentence — it
is proved on a large but incomplete region. Line 27 does say "large parts", but line 6 is the sentence
that gets quoted. Change to "proves ... for one classical family, and for large parts of a second".

## m4. Theorem A's other headline numbers — PASS

Each traced to a log line that says what the text says:

| claim | location | log line |
|---|---|---|
| 627-rung ladder, 7.06 h | `ARTICLE.md:90` | `m_ladder_log_2026-09-02_full.txt`: `627  5026  316  n <=  1259` and `# total CPU time 7.06 h; all 627 rungs present, every shift equals the prediction, no None` |
| `min M = 1.807, 1.804` | `ARTICLE.md:91` | reproduced live by me: `1.807003`, `1.803808` |
| 930 boxes, `0.800004` | `ARTICLE.md:93` | `dense_certificate_a_rerun_2026-09-03.txt`: `certified: 930 boxes ok, 213 skipped, worst N g lower bound 0.800004 at eps=1.61e-04, v=0.4136` |
| `4/5` exact, `S(0)=176/525`, 18 bands, `min >= 0.3206` | `ARTICLE.md:94` | `dense_certificate_b_rerun2_2026-09-03.txt`: `constant term of N g (must be 4/5 exactly, zeta-free): 4/5`; `S(0, zeta) = 176/525`; first band `min(S - losses) = +0.3206`; 18 band lines counted |
| `[a^3]/b^2 = 176/175` | `ARTICLE.md:92` | only via `VALIDATION_SPARSE_CERTIFICATE_FULL:9,233` — see M2 |

## m5. Theorem B's headline numbers — PASS except one

| claim | location | log line |
|---|---|---|
| 1000 rungs, degree `2j-1`, shift `j+2` | `ARTICLE.md:106` | `sibuya_ladder_log_full.txt`: `1000  1999  1002  ok` and `# total time 3.53 h; all 1000 indices present, j = 1 the identity, every shift j + 2, every check ok` |
| four certificates, exit 0 | `ARTICLE.md:114` | `sibuya_theorem_full_2026-09-03.txt`: four `C'.` lines each `exit 0 ... : True`, closing `EXIT 0` |
| `[a^3]/b^2 = 2/3`, residual `0.29` | `ARTICLE.md:107` | `sibuya_sparse_certificate_2026-09-03.txt`: `[a^3]/b^2 at b = 0: 2/3   (2/3: True)`; `min of [a^3]/b^2 - (iii) - (v) - (vi) = 0.290647` |
| `84 boxes`, margin `0.163` | `ARTICLE.md:108` | `sibuya_dense_a_top_2026-09-03.txt`: `84 boxes ok, 6 skipped, worst margin ... >= 0.162638` |
| `S(0,zeta) = 1/3`, 18 bands, `min >= 0.2398` | `ARTICLE.md:109` | `sibuya_dense_b_2026-09-03.txt`: `S(0, zeta) = 1/3`; first band `min(S - losses) = +0.2398`; 18 bands counted |
| corner window verdict | `README.md:55` | `sibuya_corner_grid_2026-09-03.txt`: `VERDICT: Sibuya (3.4) holds for 33 <= j' <= 399 and every N (exact recursion N <= 1000000 by sibuya_harmonic.py, grid to 485165195, H-model beyond): True`, `EXIT 0` |
| top corner `j' in [1,32]` | `ARTICLE.md:110` | `sibuya_harmonic_2026-09-03.txt`: `VERDICT: ... j' = N - j in [1, 32] and every N >= 1001 + j' ...: True`, `EXIT 0` |
| `303 boxes`, margin `1.03` | `ARTICLE.md:111` | `sibuya_dense_a_top2_exploratory_2026-09-03.txt`: `303 boxes ok, 473 skipped, worst margin ... >= 1.033092` — number PASS, framing per M6 |

**The one miss:** `ARTICLE.md:108` gives C2' the margin "`1.0e-5`". The script's own log says
`worst margin N g - 1/(3 + 3 eps - theta) >= 0.000005` (5e-6). `1.018e-5` is the **validator's**
independent sweep (`VALIDATION_SIBUYA_DENSE_A:54`, "own sweep: 266 ok, 43 skipped, worst margin
1.018e-5"). The article attributes the script's box count and the validator's margin to the same run.
Use `5e-6` for `sibuya_dense_a.py`, or attribute both numbers.

**"1005 Sibuya rungs" appears nowhere** in README.md or ARTICLE.md — both say 1000, which the log
supports. No defect; noting it since the number was raised.

## m6. `theorem.py` rewrites an artifact in `results/` on every run

`theorem.py:104-113` shells out to `lab/ladder_log_merge.py 627`, whose printed result is
`OK: wrote ...\results\m_ladder_log_2026-09-02_full.txt`. So merely *verifying* the theorem overwrites
a file in `results/`. It is a deterministic re-merge of immutable per-range logs, so nothing is lost —
but it collides with the experiment contract's "preserve raw outputs as immutable artifacts under
results/raw/ (never rewrite)", and it means the merged log's mtime tells you nothing about when the
7.06 h of computation happened. Have the checker verify without writing, or write to a scratch path.

## m7. A stale, verdictless 9715-second log sits in `results/`

`results/sibuya_dense_a_top2_2026-09-03.txt` (86 KB, 11:51 today) ends mid-sweep:

    ... 23921 ok, 54579 skipped, stack 16  [9715 s]

No `VERDICT`, no `EXIT`. It is a killed or still-running `--top2` attempt, not cited by any document,
sitting next to the log that *is* cited. Keeping failed runs is correct; leaving an ambiguous one
unlabelled in a publication package is not. Add a one-line header saying it was killed, or move it out
of the release-adjacent set.

## m8. Metadata — mostly PASS, three fixes

Both files parse (`json.load` OK, 9 keys; `yaml.safe_load` OK, 8 keys) and CFF 1.2.0's required fields
(`cff-version`, `message`, `title`, `authors`) are present.

- **ORCID: PASS.** `0009-0005-5660-2603` in both files, matching `README.md:3`, `ARTICLE.md:3` and the
  repository's established identity (`article/DATA_LOG.md:114`, `article/PAPER_*.md:4`).
- **Licence: PASS.** MIT in `.zenodo.json`, `CITATION.cff` and `README.md:96`, consistent with the
  repository `LICENSE` (MIT). `README.md:96`'s relative path `../../../LICENSE` resolves correctly from
  `release/`.
- **Fix 1 — type mismatch.** `.zenodo.json` declares `"upload_type": "publication"`,
  `"publication_type": "preprint"`; `CITATION.cff` declares `type: dataset`. Same artifact, two
  categories. Pick one (preprint, given the content).
- **Fix 2 — wrong relation.** `.zenodo.json` `related_identifiers` asserts this record
  `isPartOf` DOI `10.5281/zenodo.21915627`. That DOI is the **spacetime-verifier / AInstein audit**
  software record — a different, completed project (per `CLAUDE.md`). Newton's inequalities are not part
  of it. Either drop the relation or change it to `isSupplementTo`/`references` with a correct target.
- **Fix 3 — name spelling.** `LICENSE` says "Copyright (c) 2026 **Andrey** Pluzhnik"; every release file
  says "**Andrei** Pluzhnik". Harmonise before a DOI freezes the transliteration.
- Missing but expected on a Zenodo/CFF record: `version`, `date-released`, `repository-code`, and
  (post-publication) `doi` in `CITATION.cff`.

## m9. AI disclosure — adequate in substance, one gap for a journal

`AI_DISCLOSURE.md` is better than most: it names the model ("Claude Opus 5, via Claude Code"), splits
roles concretely, states "The model was not permitted to decide scientific truth", defines the
PROVED/VERIFIED/CONJECTURED vocabulary, and commits to keeping refuted claims visible. That satisfies
the substance of ICMJE/Nature-style AI policies.

The gap: those policies require the **human author to take full responsibility for the entire content**,
including AI-generated parts, and to state that the AI is **not an author**. Neither sentence appears.
Add both verbatim. Also: line 169 points at "the parent repository's night report" without a path —
give the file. And the disclosure claims "No floating-point number enters a comparison anywhere in the
chain", which M5 shows is false.

---

# Summary table

| # | severity | file:line | defect |
|---|---|---|---|
| B1 | blocker | `ARTICLE.md:119`, `README.md:78`, `.zenodo.json`, `CITATION.cff`, `theorem.py:141` | "independently validated" / "complete" contradicts `VALIDATION_DENSE_CERTIFICATE_B:235` (recommends `experimentally-supported`) |
| B2 | blocker | `fig3d_sibuya.py:118` → `data/sibuya_3d.svg` | label says "sinking towards **4/5**"; the Sibuya floor is 1/3 |
| B3 | blocker | `fig3d_sibuya.py:63` | floor plane drawn at 0.32, labelled 1/3; inflates the margin ~60x |
| B4 | blocker | `ARTICLE.md:127` | residual region wrong; true residual is `j >= 1001` ∧ `theta > 0.9` ∧ `j' >= 400` |
| B5 | blocker | `README.md:34`, `ARTICLE.md:26`, `ARTICLE.md:127` | three inconsistent statements of Theorem B |
| B6 | blocker | `scripts/sibuya_top_w.py`, `results/sibuya_top_w_2026-09-03.txt` | shipped script exits 1, uncited, and is the piece that would close B4 |
| B7 | blocker | `LIMITATIONS.md:7,135,138,143,159` | stale by a day; contradicts README/ARTICLE; no Sibuya section |
| M1 | major | `ARTICLE.md:97` | "409 s, three VERDICT True, exit 0" — no such log exists |
| M2 | major | `ARTICLE.md:92` | `sparse_certificate_full.py` has no run log; `0.1417` only second-hand |
| M3 | major | `README.md:52-57` | all six reproduction commands fail from the repository root |
| M4 | major | `README.md:71` | dead link to `OUTREACH.md` |
| M5 | major | `README.md:61`, `theorem.py:87-97`, `sibuya_theorem.py:82-91` | "no float enters a comparison" is false in both coverage lemmas |
| M6 | major | `ARTICLE.md:93,108,111` | skipped-box counts omitted; E' is 303 ok / 473 skipped |
| m1 | minor | `ARTICLE.md:24`, `README.md:27` | (3.4) presented as Sibuya's inequality; it is this project's (correct) translation of his monotonicity claim |
| m2 | minor | `ARTICLE.md:66` | chain PASS; add that it needs `g > 0` |
| m3 | minor | `README.md:6` | "proves ... for two classical families" overstates Theorem B |
| m4 | pass | — | Theorem A headline numbers all trace to logs |
| m5 | minor | `ARTICLE.md:108` | C2' margin `1.0e-5` is the validator's; the script logs `5e-6` |
| m6 | minor | `theorem.py:104` | verification rewrites `results/m_ladder_log_2026-09-02_full.txt` |
| m7 | minor | `results/sibuya_dense_a_top2_2026-09-03.txt` | 9715 s log with no verdict, unlabelled |
| m8 | minor | `.zenodo.json`, `CITATION.cff`, `LICENSE` | type mismatch; wrong `isPartOf` DOI; Andrei/Andrey |
| m9 | minor | `AI_DISCLOSURE.md` | missing "author takes full responsibility" + "AI is not an author" |

---

# What I could not check in the hour

- The five certificates were not re-run (`--full` runs total ~21 min plus `sibuya_corner_grid` ~9 min);
  I read their logs and their independent validation reports instead. The mathematics inside
  `dense_certificate_a/b.py` and `sibuya_dense_a/b.py` (Edgeworth remainder machinery, Bernoulli lemma,
  `Ser` sup-bound bookkeeping) I did not re-derive — the `VALIDATION_*` reports do that, and I took
  their PASS/FAIL findings at face value, which is how B1 was found.
- I did not obtain Sibuya 1988 directly; m1 rests on the verbatim quotation in
  `LITERATURE_VERDICT_STIRLING_NEWTON_GAP_2026-09-02.md:80-84` and the independent OCR re-derivation in
  `VALIDATION_SIBUYA_DENSE_B_2026-09-03.md:50`. Two independent readings of the same page agree, which
  is reasonable but not the same as my own copy.
- `PAPER.md` (41 KB working notebook) was out of scope and unread.

---

# FINAL VERDICT

**NOT READY to be published as a Zenodo record or sent to researchers.** Seven blockers stand between
this package and a defensible submission, and two of them (B1, B2) are the kind a hostile referee finds
immediately and never forgets: a claim-gate violation asserting independent validation that the
validator explicitly withheld, and a headline figure printing the wrong constant.

**Minimum set of changes before publication:**

1. **B1** — either close items (i)-(iii) of `VALIDATION_DENSE_CERTIFICATE_B` and rerun, or strike
   "every certificate is independently validated" / "each was validated by an independent
   re-implementation" from `ARTICLE.md:119`, `README.md:78`, `.zenodo.json`, `CITATION.cff` and
   `theorem.py:141`, and downgrade "Theorem A: complete" to the honest form. **The word "proved" for
   Theorem A is not currently earned by the gate.**
2. **B2, B3** — fix `fig3d_sibuya.py:118` (4/5 → 1/3) and `:63` (0.32 → 1/3), regenerate `sibuya_3d.svg`.
3. **B4, B5** — adopt one statement of Theorem B, the true residual
   `j >= 1001` ∧ `j/(n-1) > 0.9` ∧ `n-1-j >= 400`, and use it identically in `README.md:34`,
   `ARTICLE.md:26`, `ARTICLE.md:127`; mark E' exploratory wherever it is mentioned, including the
   abstract.
4. **B6** — document or remove `sibuya_top_w.py` and its `EXIT 1` log.
5. **B7** — rewrite `LIMITATIONS.md` to today's state and give Theorem B a section.
6. **M1, M2** — generate the two missing logs (`theorem.py --full`, `sparse_certificate_full.py`).
7. **M3, M4** — fix the six README commands and the `OUTREACH.md` link.
8. **M5** — move the two coverage lemmas off `float` and onto `fmpq`/`arb`, then the absolute claim is
   true and can stay.

**What is genuinely strong, and should not be lost in the rewrite:** the region decomposition tiles
correctly with no circularity (I reconstructed both unions from the scripts' own constants and they
close — Theorem A's `A ∪ B ∪ C1 ∪ C2 ∪ C3` is airtight, and Theorem B's gap is one honest region, not a
hidden one); the chain `n(R-1) > n g > n/(3n-j)` is correct as written; Sibuya's (3.4) is transcribed
faithfully and credited everywhere; the priority claims are, if anything, under-stated; every headline
number but one traces to a log line that says what the text says; the two fast scripts run and are
deterministic; and the practice of keeping refuted routes and failing runs visible is exactly right.
This is a package with a bookkeeping problem, not a mathematics problem. Fix the eight items and it is
publishable.

**Publication remains a deliberate human action.** Nothing here authorises submission; this report is a
referee's list, and the author decides.
