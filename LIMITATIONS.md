# Limitations

*What separates the present state from a referee's "proved". Written 3 September 2026; every statement here is
checked against the logs in `../results/`.*

## Theorem A (`M(n,t) > 4/5`)

Nothing in the mathematics is outstanding. The proof is a finite decomposition into five pieces, each an exit-0
script with a logged run, each independently validated by an agent that did not import the code it validated.
What remains is what remains for any computer-assisted proof:

- **A referee has not read it.** The certificates are exact and self-checking, but no human has audited the
  derivation of the lemma constants in `sparse_certificate_full.py` (the tail and moment lemmas) beyond the
  independent validator's own re-derivation.
- **The ladder is an artifact, not a proof by hand.** 627 polynomial certificates, one per index, checked by
  `lab/ladder_log_merge.py`. A reader who does not run the script has to trust the log.
- **The dense certificates rest on a shared engine.** `dense_certificate_a.py` supplies the Edgeworth remainder
  machinery used by `dense_certificate_b.py` and by the Sibuya ports; a defect there would touch several pieces
  at once. It has been validated once, independently.

## Theorem B (Sibuya's conjecture)

- **One region is open:** `j >= 1001` with `j/(n-1) > 0.9` and at least 803 indices missing. Section 7 of
  `ARTICLE.md` states it exactly and names the two instruments that stop at its boundary.
- **The `--top2` sweep is exploratory.** It skipped 473 of 776 boxes and self-labels as not a certificate; it is
  quoted in section 5 as evidence, never as part of the proved range.
- **`sibuya_top_w.py` exits 1.** It is the instrument written for the open region; it certifies boxes away from
  the cusp of that region and stops at the cusp. It is shipped as a starting point, not as a certificate.
- **The sparse port's lemma constants** (`FMIN`, `CF`, the `C_i` bound beyond `i = 60`) were derived by analogy
  with the 4/5 certificate and then validated independently (`results/VALIDATION_SIBUYA_SPARSE_2026-09-03.md`,
  eight items PASS, one multiplier corrected).

## Both

- **Floating point is used for bookkeeping.** Box and step sizes, loop bounds and printing use doubles. Every
  comparison a conclusion rests on is exact (`fmpq`) or certified interval (`arb`/`acb`). The coverage lemmas in
  `theorem.py` and `sibuya_theorem.py` were converted to interval comparisons on 3 September after a referee
  review flagged them.
- **The physics application is not done.** The centred-square spectrum comes from the partial-wave positivity of
  deformed Veneziano amplitudes; nothing here is carried back to that problem.
- **AI assistance.** See `AI_DISCLOSURE.md`. The proofs and code were produced in collaboration with an AI
  assistant under the author's direction.
