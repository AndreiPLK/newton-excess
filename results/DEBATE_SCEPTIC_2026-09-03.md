# Sceptic's brief: the Sibuya (3.4) coverage claim, 3 September 2026

**Header time.** `date` could not be run: the Bash tool is disabled in this session (tool error
`No such tool available: Bash`). The date is taken from the session clock, **2026-09-03**, and the
time of writing is therefore *not* a measurement. Treat the timestamp as unverified.

**Role.** Sceptic in an internal debate. Everything below is an attack on the claim, not a review of
the work's value. I am not an accredited mathematician; where I could not settle a point by reading
the code I say so instead of guessing.

**Claim under attack.** Sibuya's 1988 conjecture `p_j^2/(p_{j-1}p_{j+1}) >= 1 + 1/(3n-j)` is proved
in the ranges of `release/ARTICLE.md` §5 **plus** the top region, by
`release/scripts/sibuya_harmonic.py`, `sibuya_corner_grid.py`, `sibuya_wedge.py`.

**Verdict in one sentence.** The `+ the top region` half of the claim does not hold: three of the
four legs it rests on have no completed exit-0 run, the fourth (the wedge) does not close its region
but pushes its boundary to `n ~ 1.2e31` and leaves an infinite family behind, and the referral that
is supposed to absorb the remainder points at an instrument the project's own ARTICLE says fails
exactly there. The `ARTICLE.md` §5 + §7 statement as *committed* (open region: `j >= 1001`,
`theta > 0.9`, `n-1-j >= 400`) is not damaged by anything I found.

---

## Ordered by how much of the claim each finding destroys

### 1. FATAL to "the top region is proved": the wedge closes no region, it moves a boundary

`release/scripts/sibuya_wedge.py:141-177`. `main()` sweeps geometric bands in `Y = 1/v` from
`Y_MAX = 1/35` down to `10**-Y_MIN_EXP` and then prints

```
VERDICT (wedge): N g > 1/(3 + 3/N - theta) on { theta >= 0.9, j' >= 803, Y >= 1e-{Y_MIN_EXP} } : True
```

The verdict is honestly conditioned on `Y >= 1e-Y_MIN_EXP`, and that condition is not a technicality:
it is an unbounded exclusion in `n`.

Translation, from the script's own variables (`eps = 1/N`, `Y = x/N`, `x` the tilt, `N = n-1`):

```
j' = SUM_{k=1..N} x/(x+k) = x [ psi(N+1+x) - psi(1+x) ] = N Y m ,   m := psi(N+1+x) - psi(1+x)
```

and `m -> log(1/Y)` for `x >> 1` (this is confirmed by the script's own corner run, which reports
`m = 1237.4887` at `L = log(1/Y) = 1238`). The region's constraint `j' >= 803` therefore forces

```
N  >=  803 / (Y * m)  ~  803 / (Y * log(1/Y)) .
```

At the cut `Y = 1e-30` (`results/sibuya_wedge_2026-09-03.txt:1`, the only band sweep that ran to
completion) this gives

```
N  >  803 / (1e-30 * 69.08)  =  1.16e31 .
```

**So the completed wedge run leaves uncovered, exactly: every `(n, j)` with `j >= 1001`,
`theta = j/(n-1) >= 0.9`, `j' = n-1-j >= 803` and `n - 1 > ~1.16e31`.** That is an infinite family,
not a finite residue. Running with `Y_MIN_EXP = 300` (`results/sibuya_wedge300_2026-09-03.txt`,
itself truncated after 400 of ~13 400 bands) would move the same boundary to `n ~ 1.2e300` and leave
an infinite family behind it.

This is the `per-instance grinding is not progress` failure named in the project's own memory: the
number of bands measures the run, not the theorem. A `Y >= 1e-30` cut is a coverage statement of the
form "for all `n` below `1.16e31`", and the conjecture is a statement for all `n`.

### 2. FATAL to the same half: `400 <= j' <= 802` has no completed certificate at any `N`

`sibuya_corner_grid.py` covers a window `JPMIN = 33 <= j' <= JPMAX` (line 36-37). The runs on disk:

| log | window | outcome |
|---|---|---|
| `sibuya_corner_grid_2026-09-03.txt` | `33..399` | VERDICT True, EXIT 0 |
| `sibuya_corner_grid800_2026-09-03.txt` | `33..800` | `grid FAILED at N = 1000000 (step could not be shrunk enough)`, **EXIT 1** |
| `sibuya_corner_grid802_2026-09-03.txt` | `33..802` | **truncated**: last line `grid: N = 1.9481e+07 ... 2500 steps` at `2501 s`, target `N_hi = 4.85e8`. No VERDICT line, no EXIT line |
| `sibuya_corner_grid_test45.txt` | `33..45` | `grid FAILED at N = 1000000`, EXIT 1 |

So for `400 <= j' <= 802` the only thing certified is the H-model tail `N >= e^20 = 4.85e8`
(`sibuya_corner_grid802_2026-09-03.txt:2`, 770/770 indices). The band
`10^6 < N <= 4.85e8` was covered only to `N = 1.95e7`, i.e. **`1.95e7 < N <= 4.85e8` with
`400 <= j' <= 802` is uncertified**, and the run that would close it is the same code path that
returned EXIT 1 at `jpmax = 800`. The `800` failure and the `802` non-failure differ by two indices;
without a completed `802` log I have no reason to believe the second run would have finished, and one
concrete reason to expect trouble (the grid step is shrunk by a shared `for _ in range(12)` loop over
*all* `j'` at once, `sibuya_corner_grid.py:139-156`, so one bad index stalls every index).

### 3. FATAL as written: the `N <= 10^6` leg of the `33 <= j' <= 399` window was never executed

`sibuya_corner_grid.py:286-288` prints

```
VERDICT: Sibuya (3.4) holds for 33 <= j' <= 399 and every N (exact recursion N <= 1000000 by
sibuya_harmonic.py, grid to 485165195, H-model beyond): True
```

**`corner_grid.py` never runs that exact recursion and never imports the module that contains it.**
`main()` (lines 259-288) calls `u_coeffs`, `tail_ok` per index, and `grid_ok(N_EXACT, N_hi)` — nothing
else. `sys.path.insert(0, HERE)` at line 34 is the only trace of `sibuya_harmonic`; there is no
`import sibuya_harmonic` anywhere in the file. The claim is a sentence in an f-string.

And the run it defers to was never made at that window. The only completed
`sibuya_harmonic.py` log (`results/sibuya_harmonic_2026-09-03.txt:1-2`) reads

```
Sibuya (3.4), the corner j' = N - j in [1, 32], every N >= 1001 + j':
  exact recursion N <= 1000000: all j' ok; ... [105 s]
```

i.e. `exact_pass(33)`, which loops `for jp in range(1, jmax)` = `1..32` (`sibuya_harmonic.py:156`).
The one attempt at a wider pass, `results/sibuya_exact_pass_803_2026-09-03.txt`, is a **zero-byte
file** — `main()`'s first `print` has no `flush=True` (line 217), so a block-buffered redirect
produces an empty file exactly when the process is killed before `exact_pass` returns.

**Consequence.** For `33 <= j' <= 802` and `1001 + j' <= N <= 10^6` there is no completed run at all.
That range is inside `theta > 0.9` (e.g. `j' = 399`, `N = 10^4` gives `theta = 0.96`, `j = 9601`), so
it is not picked up by C1'/C2'/C3', and `j >= 1001` puts it outside the ladder. This punctures the
`33 <= j' <= 399` window **that the committed `release/README.md:34-36` and `ARTICLE.md` §5 row D''
already claim**, not only today's addition.

This one is probably cheap to repair — it is a missing 20-40 minute run, not a wrong argument — but
until the log exists, row D'' is not an artifact-backed claim, and `REFEREE_REVIEW_RELEASE_2026-09-03.md:373`
accepted the f-string as evidence.

### 4. SEVERE: the wedge's own artifact ends `EXIT 1`, and the failing check was removed from `main()`

`results/sibuya_wedge_2026-09-03.txt:9-10`:

```
FAILED (corner): corner FAILED at L = 1238, r in [1.201e+00, 1.201e+00]:
  {'slack': [+/- 1.78e-7], 'E': [0.000484545678... ], 'sigma2': [1030.13...], 'm': [1237.4887...], 'W': [0.600...]}
EXIT 1
```

The current `sibuya_wedge.py:141-177` `main()` **does not call `certify_corner_sweep`** (defined at
line 89 and now unreachable from `main`). In its place, lines 170-173:

```python
# NOTE on the corner Y -> 0 at FIXED j': there the leading orders of N g and of the target cancel exactly
# (measured: the slack encloses 0 to +/- 1e-12 at j' = 20000, L = 10^9), so the crude corner lemma cannot
# decide it -- that limit is the H-model's territory, certified per index by sibuya_corner_grid.tail_ok.
```

Three objections, in increasing order of seriousness.

(a) **The artifact and the script disagree.** The only wedge log on disk is from the version that ran
the corner and exited 1. A future reader who runs the current script gets EXIT 0 and a verdict; the
log says EXIT 1. Whichever is right, they cannot both be the evidence.

(b) **The referral is to an instrument documented to fail there.** `ARTICLE.md:141-143`, in the
project's own words: *"the `H`-model's Taylor-shift certificate stops working around 1000 missing
indices at any starting height, because `Fhat` is not yet monotone for `H` far below the index."*
The failing corner box is at `m = 1237.5`, `r = 1.201`, i.e. `j' = m/r = 1030` — inside the band where
the H-model is said to stop. And `tail_ok` has in fact never been run above `j' = 802`.

(c) **The corner is not a bounded set of indices.** `certify_corner` is invoked over `r in [0, L/JP_MIN]`
for `L -> infinity`, so `j' = m/r` is unbounded in it. "Certified per index by
`sibuya_corner_grid.tail_ok`" is a per-instance promise over an infinite index set even if every
individual index passed.

**What the failing number actually says, and it is worth more than the failure.** In that box
`W/Y = 0.600`, `m = 1237.5`, so the positive leading term is `W/(Y m) = 4.85e-4`; the reported
Edgeworth enclosure is `E = 4.845e-4`. The two agree to three digits. That is not a precision
problem — it is the statement that at `Y -> 0` with `j' ~ 10^3` the margin in `N g - 1/(3+3eps-theta)`
**cancels at leading order**, so the inequality is asymptotically tight in that corner and any method
that keeps only the leading order will report `[+/- 1e-7]` for ever. This is a substantive finding
about the object, and it is the reason the corner is hard. It should be recorded as such rather than
routed around.

### 5. MODERATE: the wedge's `eps` cover is built at the band midpoint, so a sliver of each band is never visited

`sibuya_wedge.py:49-53`:

```python
def eps_max_for(Y):
    mu = mu_of(arb(0, 1e-300), arb(Y.mid()))
    return arb(float(mu.lower()) / JP_MIN)
```

Two substitutions, neither justified in the file:

* **`Y.mid()` instead of the band.** `theta` decreases in `Y` (the script's own premise,
  `sibuya_top_w.py:236`), so `mu = 1 - theta` increases in `Y`, so the true
  `eps_max(Y) = mu(Y)/JP_MIN` at the band's top `y_hi` exceeds the value computed at the midpoint by
  the same relative amount as `mu` varies over the band. With `BAND = 1.05` (line 41) that is about
  2.5 %. Points with `Y` near `y_hi` and `eps` in `(mu(mid)/803, mu(y_hi)/803]` satisfy `j' >= 803`,
  lie in the claimed region, and are outside the box the script actually covers. In index terms the
  uncovered sliver sits at `803 <= j' <~ 843` near the top of every one of the 1343 bands.
* **`eps = arb(0, 1e-300)` instead of the `eps` range.** `mu` also depends on `eps` at fixed `Y`;
  evaluating at `eps ~ 0` and using `mu.lower()` is asserted to be conservative but the sign of
  `d mu / d eps` is never established in the file. If `mu` increases in `eps`, the sliver above is
  wider than 2.5 %.

Also `float(mu.lower())` drops from a certified `arb` to a double before dividing: a float decides
the edge of a covering.

### 6. MODERATE: the wedge's summary cannot distinguish "certified" from "skipped"

`certify_band` (`sibuya_wedge.py:55-86`) returns `f"{n_ok} ok, {n_skip} skipped"` but `main()`
discards that string (line 153 binds it to `msg` and never prints it); the run-level output is only
the band count and the worst margin. A band in which every box was *skipped* is indistinguishable in
the log from a band in which every box was *certified*. Compounding this, `certify_band` calls
`TW.certify_box`, whose skip rule uses `TW.JP_MIN = 400` (`sibuya_top_w.py:36, 163`), not the wedge's
`JP_MIN = 803` — the direction is conservative (fewer skips) but it means the wedge's stated region
and the region its skip rule enforces are two different sets, and neither is reported per band.

The constant `worst margin 1.628e+00` across all 1343 bands (`sibuya_wedge_2026-09-03.txt:2-8`) is
consistent with the worst box being found in the first band and never revisited, which is plausible;
it is also consistent with almost everything being skipped. The log cannot tell them apart.

### 7. The interfaces I was asked to check that **do** tile (one sentence each)

* **`t <= 1000` vs `t >= 1001`.** Sound. `certify_box` skips only when `x_hi = theta.upper()/eps_lo < T_MIN = 1001`
  (`sibuya_top_w.py:160,163`), i.e. only when every integer `t` in the box is `<= 1000`; boxes touching
  `eps = 0` get `x_hi = inf` and are never skipped by this rule.
* **`theta <= 0.9` vs `theta >= 0.9`.** Sound, and actually checked: `sibuya_top_w.py:236-242` asserts
  `theta(Y = 1/35, eps in [0, 1/1401]).upper() <= 0.89833 <= 0.9`, so `{theta >= 0.9}` is inside
  `{Y <= 1/35}` where the sweep starts. Caveat in finding 9 below.
* **Grid end vs H-model start in `corner_grid`.** Sound with overlap: the grid runs to
  `N_hi = int(e^20) = 485165195` and `tail_ok` covers `H_N >= 20`, i.e. `N >= e^{20-gamma} = 2.72e8`,
  so `[2.72e8, 4.85e8]` is double-covered and there is no seam.
* **`exact_pass` end vs band start in `harmonic`.** Sound: `exact_pass` covers `N <= 10^6` and
  `certify_jprime(jp, N_EXACT)` starts its bands at `10^6`, with `D = 0` / `Na = 10^6` included.
* **Corner-grid window vs wedge floor.** *Not* an interface but a gap: `corner_grid` completed at
  `j' <= 399`, the wedge starts at `j' >= 803`. Findings 2 and 3.

---

## The instruments: the assumption that would silently break each

### 8. `corner_grid.tail_ok`: the constant `0.6 (2j')^2 e^{gamma-H}/H^2` is asserted, never checked, and is asserted in the one regime where its derivation is invalid

`sibuya_corner_grid.py:206-217, 245`. The whole H-model leg — which is the *only* completed leg for
`400 <= j' <= 802`, and the leg the wedge's NOTE defers its corner to — rests on

```
E_k(N) = Ehat_k(H)(1 + eta_k),   |eta_k| <= eps(H) := 0.6 (2j')^2 e^{gamma-H} / H^2
```

The derivation sketched in the docstring is: `Gamma(N+1+x)/Gamma(N+1) = e^{psi(N+1)x} exp(SUM_{r>=2} ...)`,
`|v_2| <= 1/(2N)`, so `E_k - Ehat_k ~ v_2 Ehat_{k-2}` — and to turn that into a *relative* bound one
must divide by `Ehat_k`, i.e. one must bound `Ehat_{k-2}/Ehat_k`. The factor `(2j')^2/H^2` is exactly
the estimate `Ehat_{k-2}/Ehat_k ~ k(k-1)/H^2`, which is the *Poisson* estimate `Ehat_k ~ H^k/k!`.

That estimate is used at `h1 = 20` with `j'` up to `802`, i.e. at `k/H = 40`. In that regime
`Ehat_k(H) = SUM_m u_m H^{k-m}/(k-m)!` is dominated by large `m`, and the file itself says so, at
line 191-193:

> a degree-800 polynomial with alternating coefficients evaluated directly at `H ~ 10^3` loses every digit

An alternating, heavily cancelling sum is precisely a sum for which `Ehat_k ~ H^k/k!` is not a valid
proxy. Nothing in the code checks `Ehat_{k-2}/Ehat_k <= 2 * 0.6 * (2j')^2 / (k(k-1)) * k(k-1)/H^2`,
or checks `Ehat_k(H) > 0` for `H >= 20, k = 800`, which the relative formulation silently presupposes.
The nominal headroom over the naive `0.5 k^2/H^2` is a factor `~4.8`, which is not a lot of room for
an unverified ratio.

**This is the single assumption whose failure would destroy the most: rows D'' of ARTICLE §5 and the
whole `400 <= j' <= 802` tail.** I could not settle it by reading; it needs the numerical test in §14.

### 9. `corner_grid.step_margin`: `h_3 <= D^3/6` is used as an equality, and the deficit is the same size as the deficit that *is* carried

`sibuya_corner_grid.py:88-119`. The bracket is

```python
def br(k):
    pol = [E[k], E[k - 1], E[k - 2] / 2, E[k - 3] / 6]
    rest = arb(0, float(S2 / 2 * Eb[k - 2])) + D4 * Eb[max(k - 4, 0)]
```

with the docstring stating `h_2 = (D^2 - S_2)/2`, `h_3 <= D^3/6`, tail `<= 1.5 (D^4/24) E_{k-4}`.
The `h_2` deficit `S_2/2 * E_{k-2}` is carried in `rest`. The **`h_3` deficit is carried nowhere.**
Exactly,

```
h_3 = (D^3 - 3 D S_2 + 2 S_3)/6   =>   D^3/6 - h_3 = (3 D S_2 - 2 S_3)/6 ~ D S_2 / 2 ,
```

so the uncarried offset in each bracket is about `(D S_2/2) E_{k-3}`, against the carried
`(S_2/2) E_{k-2}`. Their ratio is `D * E_{k-3}/E_{k-2} ~ D k/H`. At the grid's own operating point
(`N ~ 10^6`, `H ~ 14.4`, `D ~ 0.02` from `Dm_guess`, `k ~ 800`) that ratio is
`0.02 * 800/14.4 ~ 1.1`. **The unaccounted term is the same size as the accounted one.**

I am not claiming the certificate is thereby false: all three brackets `E_{j'-1}, E_{j'}, E_{j'+1}`
are biased the same way and `F = A^2 - T B C` cancels first order, which is the whole design of the
step. I am claiming that the cancellation of the `h_3` bias is *assumed*, not computed, while the
structurally identical `h_2` bias is computed — and that an "upper bound used where the sign matters"
is exactly the failure mode this bracket was built to avoid. The fix is one line: put
`arb(0, (Dm * S2 / 2 * Eb[k-3]).upper())` into `rest`, re-run, and see if the margins survive
(worst relative margin logged is `2.08e-9` at `j' <= 399` and `4.18e-10` at `j' <= 802` — these are
not large).

Related and unproved in the same function: the `m >= 4` tail bound uses `Eb[max(k-4,0)]` as a majorant
of every `E_{k-m}`, `m >= 4`, which needs `E_j` decreasing in `j` over `j <= k-4` — true when `k > H`
but never checked, and false near `j = 0` where `E_0 = 1` dwarfs `E_{k-4}` (the compensation is that
`h_k` is astronomically small; the argument is fine but is nowhere written).

### 10. Both `tail_ok`s use `1/N <= e^{gamma-H}`, which is **false** as stated

`sibuya_harmonic.py:20` (`1/N < e^{gamma - H} (H_N > log N + gamma)`) and
`sibuya_corner_grid.py:211-212` (`1/N <= e^{gamma-H}`), both load-bearing (they convert the `1/N`
defect into a function of `H` so it can be beaten by a polynomial inequality).

From `H_N = log N + gamma + 1/(2N) - 1/(12N^2) + ...` one gets `H_N > log N + gamma`, hence
`gamma - H_N < -log N`, hence **`e^{gamma - H} < 1/N`** — the opposite of what is used. The parenthetical
`(H_N > log N + gamma)` is the correct fact and the inequality drawn from it is the wrong way round.

Magnitude: the true statement is `1/N = e^{gamma-H} e^{1/(2N) - ...} <= e^{gamma-H}(1 + 1/(2N)*1.001)`,
a relative overshoot of `< 1e-9` at the smallest `N` in either tail (`2^48` and `e^20`). Against the
headroom in the logs (`c_min = 1.5e-37` vs `e^-H` term `1.6e-62` in `sibuya_harmonic_2026-09-03.txt:6`)
this changes nothing numerically. **But the lemma as written in the released scripts is false, and a
referee who checks it will stop there.** Both docstrings need the `e^{1/(2N)}` factor, or the
substitution `1/N <= 1.001 e^{gamma-H}`.

### 11. The wedge never executes the two preconditions of the engine it imports

`sibuya_wedge.py:32` imports `sibuya_top_w` with `sys.argv` stubbed, which (correctly) does not run
`TW.main()`. But `TW.main()` is where the two global preconditions live:

* `DA.bernoulli_lemma(SA.KST)` (`sibuya_top_w.py:233-235`) — the certificate for the Stirling
  remainder that every call to `phi_ratio`/`phi_stirling_real` depends on;
* the coverage assertion `theta(Y = 1/35).upper() <= 0.9` (lines 236-242) — the statement that
  `{theta >= 0.9}` is contained in the swept `Y` range at all.

Neither runs in a wedge run. Both hold (they are in `results/sibuya_top_w_2026-09-03.txt:1-2`, `True`
with 2185 boxes and `0.89833`), so this is a bookkeeping hole rather than a mathematical one — but the
wedge's exit-0 artifact does not contain them, and the `top_w` artifact that does contain them **exits 1**
three lines later.

### 12. `certify_corner`'s "everything scales like `Y`": the log *is* handled; the cancellation is the real problem

The prompt suggested the corner lemma wrongly assumes pure `Y`-scaling when there is a log. Reading
`sibuya_top_w.py:314-345`, that particular charge does not stick: `corner_scaled` explicitly carries
`log w_2 = -L + log(1+r)` (line 333) and the docstring states `kt_j/Y = A_j L + G_j(r) + O(Y)`. The
`L` is there.

What is wrong is different and worse: the lemma is built so that `L` *cancels* between the margin
(`~1/(Y L^2)`) and the defect (`~1/(Y L)` times `|kappa_2 - 1| <= C r/L`), i.e. the whole conclusion
sits on the ratio of two quantities that are equal to leading order — and the run confirms they are
equal to three digits (finding 4). A lemma whose design principle is "the two sides have the same
leading order and I keep only the leading order" cannot conclude anything; the sign lives in the next
order, which the lemma discards.

One genuine inconsistency inside it: `Y` is kept as the ball `(0, e^{-L}]` (line 328) while
`log Y` is replaced by the single point `-L` (line 333), i.e. the same variable is treated as an
interval in one term and as its own supremum in another, inside one expression.

### 13. `phi_ratio`'s Stirling series at small argument: **I attacked this and it holds**

`sibuya_top_w.py:46-74` calls the Stirling series through `u = eps/w`, i.e. at effective argument
`z = 1/u`. In `corner_scaled` (line 335) `z = (1+r)/r` can be as small as `~1.8` (the failing box has
`r = 1.201`), and `DA.bernoulli_lemma` is stated with `X0 = 60`, which looks fatal.

It is not. `bernoulli_lemma` (`dense_certificate_a.py:99-142`) certifies
`|(t/2)coth(t/2) - SUM_{k<=K} B_{2k}t^{2k}/(2k)!| <= |B_{2K+2}| t^{2K+2}/(2K+2)!` in the *Laplace
variable* `t`, not in the argument `z`. Through
`psi^{(n)}(z) = (-1)^{n+1} INT_0^inf t^{n-1} e^{-zt} [1 + t/2 + SUM B_{2k}t^{2k}/(2k)!] dt`
and `e^{-zt} t^{n-1} >= 0`, the pointwise bound in `t` transfers to the remainder for every `z > 0`,
and `phi_ratio`'s `rem = |B_{2K+2}| (n+2K+1)!/(2K+2)! u^{2K+2}` matches the `k = K+1` term of that
integral exactly. So the Stirling remainder is uniform in `z`, and `X0 = 60` is irrelevant to it.
One line of attack that found nothing.

---

## Interval arithmetic

### 14. Floats deciding coverings and bounds (small, real, and everywhere)

* `sibuya_corner_grid.py:141` — `Dm_eff = arb(float(harmonic(N2) - H))`. This **discards the radius**
  of a certified `arb` difference and rounds the midpoint to a double. `Dm_eff` is then used as the
  right end of the bracket `D in [0, Dm]`. If the true `H_{N2} - H_N` exceeds the rounded double by
  one ulp, the bracket does not contain `D` and the step is not an enclosure. Same at line 112
  (`Db = arb(float(Dm.upper())/2, float(Dm.upper())/2)`) and 99 (`D4`).
* `sibuya_corner_grid.py:142` — `S2 = float(...)`, used as an **upper** bound for `SUM 1/i^2` at
  line 103; the midpoint of a ball is not an upper bound.
* `sibuya_wedge.py:57, 64, 83` and `sibuya_top_w.py:250-251, 279-288` — band and box endpoints are
  formed as `arb(float(a)+float(b))/2, (float(b)-float(a))/2)`. Both midpoint and radius are rounded
  in `float`, so `mid +/- rad` may fail to cover `[a, b]` by ~1 ulp at each end, at every one of the
  ~1343 bands and every bisection.

Each of these is `O(1e-16)` relative against margins of `1e-9` to `1e-10`, so I do not claim any of
them changes an outcome. I do claim that a computer-assisted proof cannot have a float deciding the
edge of a covering, and that a referee will say the same. The repair is mechanical: use
`.upper()`/`.lower()` on `arb`s and never round through `float` on a bound.

(`arb(mid, rad)` with an `arb` midpoint is, as far as I can determine from python-flint's semantics,
`arb_set` followed by `arb_add_error`, i.e. it *preserves* the midpoint's own radius; the idiom at
`sibuya_harmonic.py:55` is therefore safe. I did not verify this against the installed python-flint
and it deserves a two-line unit test, because if it is wrong every band in `harmonic.band_ok` loses
its power-sum enclosure.)

### 15. Monotonicity assumed rather than proved

* **Proved on inspection, no objection:** `E_k` increasing in `N` (`corner_grid.py:138` — true, adding
  a positive root increases every elementary symmetric function; and the comment at line 154 that `Eb`
  stays valid after the step shrinks is correct). `T(N, j')` decreasing in `N`
  (`sibuya_harmonic.py:120-122` — true, both factors decrease). Both are used in the right direction.
* **Assumed:** `sibuya_top_w.py:365-373` — `r_worst = arb(r.upper())` is substituted for the whole `r`
  interval in the Edgeworth half, on the stated grounds that *"sigma^2 = kt_2/eps decreases in r and
  the Edgeworth defect E grows with r"*. Neither is established. Meanwhile `m` and `WY` keep the full
  `r` ball, so the same variable is a point in one half of the final inequality and an interval in the
  other. **This is the place where a wrong monotonicity produces `"ok"` on a false statement**, not a
  spurious `"fail"`.
* **Assumed:** `sibuya_wedge.py:51` — `theta` monotone in `Y` and `mu` monotone in `eps` (finding 5).
* **Assumed:** `sibuya_top_w.py:415-418` — `1 - G(y) = 1/(2y) - 1/(3y^2) + ...` treated as positive and
  truncated after two terms; the alternating-decreasing property is not shown.

### 16. Division by a ball containing zero: **checked, and it is harmless where I expected it**

`certify_corner` passes `arb(k2Y.lower()) * Y` as `k2t` and `Y * r_worst` as `half_eps` into
`DA.kappa2_enclosure` / `kappa2_window_average`, with `Y = arb(0, e^{-L})` — balls containing zero,
and `k2t` appears in `dense_certificate_a.py:478` as `apow(arb(k2t.lower()), 7)` in a denominator.
That would be fatal, **except** that both call sites pass `geo_K=(...)` (`sibuya_top_w.py:202, 392`),
which takes the `else` branch at line 482-488 where neither `half_eps` nor `k2t` appears at all. So
the zero-containing balls are dead arguments. Worth a comment in the code, because the next person to
drop `geo_K` gets a silent `nan`.

Elsewhere `M2e/M0e` (line 587) divides by `M0e`, which can contain zero; the guard `M0c.lower() > 0`
is applied by the callers *after* the division, so the division itself produces `nan`/`inf` and the
guard then rejects it. Safe, but by luck rather than by design.

### 17. An upper bound where a lower bound is needed — one instance, defensible only under a reading that is not written down

`sibuya_top_w.py:222-223`:

```python
x_lo = arb(max(float(T_MIN), float(N_lo * theta.lower())))
one_minus_Gx = arb(0, (1 / (2 * x_lo)).upper())
```

When the box's true smallest `x` is below `1001`, the code **substitutes `1001`**, which *shrinks* the
uncertainty ball on a term that is added to `N g`. The skip rule only guarantees `x_hi >= 1001`, not
`x_lo >= 1001`, so straddling boxes exist. The enclosure is then valid only at the points of the box
with `x >= 1001` — which is exactly the claimed region, so the step is defensible. But the script's
own output says `N g > ... on { theta >= 0.9, t >= 1001, ... }` without ever stating that the box
conclusion is conditional on `t >= 1001` rather than uniform on the box. That distinction has to be
written down, or a reader will (correctly) read the code as unsound.

---

## The question I was told matters most

### 18. In `sibuya_top_w.certify_corner`, is the near-cancellation regime resolved as "fail" (safe) or could rounding make it "ok" (unsafe)?

**The decision itself is safe.** `sibuya_top_w.py:421-423`:

```python
if (lhs - rhs).lower() > 0:
    return "ok", {...}
return "fail", {...}
```

`arb` arithmetic is rigorous, so `lhs - rhs` always contains the true value; a ball straddling zero has
`.lower() <= 0` and falls through to `"fail"`. The logged failure
(`slack: [+/- 1.78e-7]`, a ball centred on zero) is exactly that path. No rounding accident can turn a
straddling ball into `"ok"` **through this comparison**, and the caller
(`sibuya_wedge.py:117-122`) escalates `"fail"` to a hard `sys.exit(1)` after 40 bisections. The
`"skip"` branch (line 365, `(m/r).upper() < jp_min`) is also sound — it fires only when every point of
the box has `j' < jp_min`.

**But "safe here" is not "safe".** The unsoundness risk in this function is not in the comparison, it
is upstream of it, in three deliberate departures from interval arithmetic that would each let a
genuinely-negative slack be enclosed by a positive interval:

1. `r_worst = arb(r.upper())` for the Edgeworth half only (finding 15) — the load-bearing one.
2. `sigma2 = arb(k2Y.lower())/r_worst` and `inv_sigma = (r_worst/arb(k2Y.lower())).sqrt()`
   (lines 369, 373): interval endpoints promoted to exact points. Both happen to be conservative
   *if* larger `inv_sigma` means larger Edgeworth coefficients, which holds for `cvals[j]`
   (`|c_j| ~ inv_sigma^{j-2}`) and for the `c_1` ball — I checked this one and it is fine.
3. `Y = arb(0, e^{-L})` as an interval but `log Y = -L` as a point, in one expression (finding 12).

And there is a fourth, in the caller: `r_hi = float((Lb/JP_MIN).upper())` (`sibuya_wedge.py:101`) uses
`L` as a proxy for `m` when deciding how far `r` must be covered, but the region is `r <= m/JP_MIN`.
In the logged box `m = 1237.4887 < L = 1238`, so the cover happened to be generous there; the sign of
`m - L` is not established, and where `m > L` the `r`-cover is short and the corner sweep would
silently miss the boxes at its own edge.

So: **the answer to the question as asked is "fail, safely".** The corner cannot be certified by
accident. It also cannot be certified at all by this lemma, for the structural reason in finding 4,
and the surrounding code has three point-for-interval substitutions that are *not* protected by `arb`.

---

## Known-answer tests I would demand before any of this is called proved

Cheap, deterministic, each one decisive against a specific finding above. I could not run them: Bash
is disabled in this session.

1. **The `eps(H)` constant (finding 8, the highest-value test).** At `H1 = 20` compute, in `arb`,
   the true `|E_k(N)/Ehat_k(H_N) - 1|` for `k = 100, 400, 800` at the `N` with `H_N` nearest 20, using
   `series_E` for `E_k` and the `u`-polynomial for `Ehat_k`. Compare against
   `0.6 (2k)^2 e^{gamma-H}/H^2`. If the measured value exceeds the bound at any `k`, the entire
   `33 <= j' <= 802` H-model tail falls, including committed row D''. Also print `Ehat_k(20)` and check
   it is positive and not the residue of a cancellation of more than the working precision.
2. **The `h_3` deficit (finding 9).** Add `arb(0, (Dm_eff * S2 / 2 * Eb[k-3]).upper())` to `rest` in
   `step_margin` and re-run `sibuya_corner_grid.py 45` over one decade of `N`. If the worst relative
   margin (`2.08e-9` logged) survives, finding 9 is closed; if it goes negative, the grid is not a
   certificate as it stands.
3. **The missing exact pass (finding 3).** Run `sibuya_harmonic.exact_pass(803)` to `N = 10^6` with
   `flush=True` on every print, and keep the log. Estimated cost from the `jmax=33` timing (105 s):
   roughly 40-50 min. This is the cheapest of all the repairs and it closes the largest committed hole.
4. **Negative control for the wedge.** Feed `certify_box` a target that is known to be FALSE (e.g.
   replace `1/(3 + 3*eps - theta)` by `1/(2 - theta)`) and confirm the sweep exits 1 quickly. A sweep
   that certifies 1343 bands and has never been shown to reject anything is not yet an instrument.
   The experiment contract in `.claude/rules/experiment-contract.md` requires this and I see no such
   run on disk for any of the three top-region scripts.
5. **`arb(mid, rad)` semantics (finding 14).** Two lines: `x = arb(1, 1e-3); y = arb(x, 1e-9);
   assert y.rad() >= 1e-3`. If it fails, `power_sum_interval` is not an enclosure and `sibuya_harmonic`
   falls.
6. **Convergence check on the grid.** Re-run one decade of `sibuya_corner_grid.grid_ok` at
   `prec_for(kmax)*2` and confirm the worst relative margin moves by less than itself. `PREC = 4000`
   and the `8*kmax` rule are justified by a single measurement in a comment (line 47); the
   experiment contract requires a precision-convergence check.

## Questions for an external expert (I could not settle these by reading)

1. Is there a *uniform in `k`* two-sided bound on `e_k(1, 1/2, ..., 1/N) / Ehat_k(H_N)` in the regime
   `k >> H_N` — i.e. does the literature give the constant that finding 8 needs, or is it genuinely new?
   (The project's recorded literature verdict on Poisson-binomial local limit theorems,
   `results/LITERATURE_VERDICT_STIRLING_SMALL_K_2026-09-03.md`, should be re-read for the reverse
   regime `k >> log N`.)
2. Is there a known uniform lower bound for the Newton gap of the reciprocal spectrum `{1/k}` — the
   route `ARTICLE.md:150-152` names as what would close the region — and does it survive the
   `theta -> 1` corner where the margin is measured to vanish at leading order (finding 4)?
3. Sibuya's own paper (Ann. Inst. Statist. Math. 40, 1988, eq. 3.4): does he state the conjecture for
   all `1 <= j <= n-2`, or with a restriction that would make the `theta > 0.9` corner vacuous? Nobody
   should certify a region the conjecture does not claim. I did not have the paper and could not check.

## What I could not damage

* The ladder leg `j <= 1000` and the three `theta <= 0.9` certificates (ARTICLE §5 rows A', C1', C2',
  C3') — I did not attack them; they have completed logs and independent validation reports on disk,
  and the coverage lemma in `sibuya_theorem.py:64-92` is an executed assertion rather than a comment.
* The `harmonic` leg for `j' <= 32`: exact pass, bands, and analytic tail all have a completed exit-0
  log, the monotonicity substitutions in `band_ok` are in the correct direction, and the `1/N`
  inequality error of finding 10 is swamped by 25 orders of headroom.
* The four interfaces listed in finding 7 tile.
* The Stirling-remainder objection (finding 13) and the divide-by-zero-ball objection (finding 16)
  both failed on inspection. Recording them so nobody spends the hour again.
* `ARTICLE.md` §7 and `release/README.md:34-36` as committed. They state the open region as
  `j >= 1001, theta > 0.9, n-1-j >= 400` and say `sibuya_top_w.py` "exits 1 at its cusp; it is not a
  certificate". That is accurate, honest, and — findings 2, 3 and 5 aside — matches the evidence.
  **The claim I was asked to attack is stronger than what the release says, and it is the excess that
  fails.**

## One-line summary for the founder

The proof of Sibuya's conjecture stands where the article says it stands (`j <= 1000` everywhere,
`theta <= 0.9` everywhere above that, and the very top of the row when fewer than 400 indices are
missing) — but today's three scripts do not extend it: two of them have no completed run, and the
third does not close its region, it only moves the boundary out to `n` about `10^31` and leaves
infinitely many cases on the far side. There is also one hole in what was already published: a
40-minute run that row D'' cites was never actually made.
