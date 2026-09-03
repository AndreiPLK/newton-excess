# Literature verdict: is there a quantitative ultra-log-concavity gap bound?

*Searched 2026-09-03 15:48 (stamp from `date`). Two attempts to delegate this search died on server errors (500, 529); the
search below was done directly, four web queries and three abstract fetches. Verdict recorded either way, as
the stuck protocol requires.*

## Verdict, one sentence

**No.** Nothing in the reachable literature bounds the ULC excess
`P(k)^2/(P(k-1)P(k+1)) - (k+1)(N-k+1)/(k(N-k))` from below in terms of the spread of the Bernoulli
parameters; the corner of Sibuya's conjecture cannot be closed by citation.

## What was found, and why each does not apply

| item | location | read | why not |
|---|---|---|---|
| Newton's inequality, strict form for distinct roots: `a_k^2 > (1+1/k)(1+1/(n-k)) a_{k-1}a_{k+1}` | standard (e.g. Wikipedia "Newton's inequalities", Stanford lecture notes on real-rooted polynomials) | full | this IS the ULC bound; strictness with no constant — the gap is exactly what is unquantified |
| J.J.F. Guo, *An inequality for coefficients of real-rooted polynomials*, arXiv:2012.03530 | abstract, main theorem | abstract_only: true | bounds the RATIO of two consecutive log-concavity gaps by `(1 ± sqrt(1-c_k))^2`, not the gap itself; useful for 2-log-concavity, not for a lower bound on one gap |
| Shenfeld & van Handel, *Extremals of the Alexandrov–Fenchel inequality for convex polytopes*, Acta Math. 231 (2023); Chan–Pak–Panova, *Equality cases of AF are not in the polynomial hierarchy*, Forum Math. Pi | abstracts | abstract_only: true | characterise EQUALITY cases only; no quantitative deficit; the complexity result says an effective description is not to be expected |
| Johnson–Kontoyiannis–Madiman, *Log-concavity, ultra-log-concavity, and a maximum entropy property* (2009/2013); Liggett 1997 | abstracts | abstract_only: true | ULC is preserved under convolution; entropy is maximal for Poisson; qualitative, no gap constant |
| *Quantitative limit theorems via relative log-concavity*, arXiv:2210.11632 | abstract | abstract_only: true | total-variation distance between measures, one log-concave relative to the other; no second-difference bound |
| Chávez & Sheng, *Stability of Khintchine-type inequalities via log-monotonicity*, arXiv:2606.19313 (2026) | abstract | abstract_only: true | log-monotonicity of moment ratios of symmetric random variables; different inequality, no bridge to elementary symmetric functions stated |
| Sibuya 1988, AISM 40(4) 693–714 | full text (earlier session) | full | the conjecture itself; **no work found citing eq. (3.4) as a conjecture, in 38 years** — the only citations are to the qualitative log-concavity results |

## Dead ends, named so nobody searches them again

- Explicit-constant local limit theorems (already recorded negative on 2 September): bound the density, not the
  second difference of its logarithm.
- Stability of Alexandrov–Fenchel: only equality cases are characterised; the quantitative deficit results that
  exist (arXiv:2503.15884 gives UPPER bounds on the deficit) go the wrong way.
- Schur-convexity of ratios of elementary symmetric functions (Springer J. Inequal. Appl. 2012; MDPI Symmetry
  2021; RGMIA notes): monotonicity under majorisation, never an explicit lower bound on the Newton excess.
- "Quantitative ULC" / "strict ULC with constant": no such phrase occurs in any indexed title or abstract.

## What this means for the corner

The corner stays open and the reason is now recorded: the general theory has the inequality but not its
size, and the size is the whole content of Sibuya's improvement. The nearest positive fact is our own
(`FIRST_ORDER_CANCELLATION_2026-09-03.md`, `LAGUERRE_POLYA_EXPLAINS_THE_H_MODEL_2026-09-03.md`): in the
limit the excess is `psi'(1+k/H)/H^2` and the finite-`N` correction enters through a bracket where the
leading term cancels. A uniform bound on that excess would be a new theorem, not a citation.

## A second observation worth stating

The absence of citations is itself evidence about the conjecture's history: it was not attacked and abandoned,
it was not noticed. It sits as one displayed line in a paper about unimodality of Stirling distributions.
