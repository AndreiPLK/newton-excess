"""SIBUYA'S CONJECTURE (1988, eq. 3.4), the last region: theta > 0.9 with at least JP_MIN missing indices.

The region `{theta >= 0.9, j' = N - j >= JP_MIN}` is, in the natural variables of the tilted weight
(`eps = 1/N`, `Y = 1/v = x/N` with `x` the tilt), a WEDGE with its vertex at the origin: for each `Y` the
admissible `eps` runs from `0` (N infinite) up to `eps_max(Y)` where `j' = JP_MIN`.  A box scheme in `(eps, Y)`
walks into that vertex and never becomes relatively thin, which is why `sibuya_top_w.py` stops there.

This script covers the wedge the way it is shaped:

  * an OUTER loop over `Y` in geometric bands (relative width small enough that the high cumulants, which enter
    as `w^{2k+n}` with `2k+n` up to 29, keep a usable enclosure -- measured: 5% is too wide, 1% is fine);
  * for each band, an INNER cover of `eps in [0, eps_max(Y)]`, subdivided only as needed.  The edge `eps = 0`
    (`N -> infinity` at fixed `Y`) is regular here: `w_2 = Y + eps -> Y > 0`, nothing degenerates;
  * the far corner `Y -> 0` is NOT closed.  The margin diverges there (measured: `6e3, 2e7, 5e16, 1e36` at
    `Y = 1e-6, 1e-10, 1e-20, 1e-40`), which is why a crude bound looked sufficient, but the leading orders of
    the two sides cancel exactly and the crude bound cannot separate them; `sibuya_top_w.certify_corner` records
    that attempt and is marked NOT LOAD-BEARING.

**WHAT THIS SCRIPT DOES AND DOES NOT ESTABLISH.**  Its verdict is conditioned on a floor `Y >= 10^-Y_MIN_EXP`,
and that floor is a bound on `n`, not a technicality: through `j' = N Y m` with `m -> log(1/Y)`, the certified
run at `Y >= 1e-30` covers `n - 1` up to about `1.16e31` and excludes everything above it.  So this is a
verification over a finite (very large) range of `n`, NOT a proof for an infinite family, and it is quoted
nowhere in `README.md` or `ARTICLE.md` as part of the proved range.

Pushing the floor does not change that.  A run to `Y = 1e-300` was attempted and FAILED at `Y ~ 2.2e-93`, with
`theta` enclosed as `1 +/- 1e-91` and `eps_max` underflowed to a denormal -- the 300-bit working precision, not
the mathematics.  Raising the precision would move the floor again and still leave an infinite family outside.
The corner is closed by an argument or not at all.

Run:  uv run python projects/qg-bootstrap/release/scripts/sibuya_wedge.py [JP_MIN] [Y_MIN_EXPONENT]
"""

from __future__ import annotations

import os
import sys
import time

from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
_argv = sys.argv
sys.argv = [_argv[0]]
import sibuya_top_w as TW  # noqa: E402

sys.argv = _argv

ctx.prec = 300

JP_MIN = int(sys.argv[1]) if len(sys.argv) > 1 else 803
Y_MIN_EXP = int(sys.argv[2]) if len(sys.argv) > 2 else 60  # bands run down to Y = 10^-Y_MIN_EXP
Y_MAX = arb(1) / 35  # theta(Y = 1/35) <= 0.9
BAND = 1.05  # geometric band ratio in Y (the ratio form of the cumulants survives 5%)


def mu_of(Y, eps):
    """the missing fraction mu = 1 - theta at (eps, Y)."""
    return 1 - TW.scaled_cumulants_Y(eps, Y, 1)[0]


def eps_max_for(Y):
    """the largest eps (smallest N) with j' = mu/eps >= JP_MIN, at this Y."""
    mu = mu_of(arb(0, 1e-300), arb(Y.mid()))
    return arb(float(mu.lower()) / JP_MIN)


def certify_band(y_lo, y_hi, verbose=False):
    """every (eps, Y) in [0, eps_max] x [y_lo, y_hi] that lies in the region."""
    Y = arb((float(y_lo) + float(y_hi)) / 2, (float(y_hi) - float(y_lo)) / 2)
    e_hi = eps_max_for(Y)
    stack = [(arb(0), e_hi, 0)]
    n_ok = n_skip = 0
    worst = None
    while stack:
        e1, e2, d = stack.pop()
        eps = arb((float(e1) + float(e2)) / 2, (float(e2) - float(e1)) / 2)
        try:
            st, info = TW.certify_box(eps, Y)
        except Exception as ex:  # noqa: BLE001
            st, info = "fail", {"why": repr(ex)[:50]}
        if st == "ok":
            n_ok += 1
            if worst is None or info["margin"].lower() < worst:
                worst = info["margin"].lower()
            continue
        if st in ("skip", "sparse-side"):
            n_skip += 1
            continue
        if d >= 40:
            return (
                False,
                f"band Y in [{float(y_lo):.3e}, {float(y_hi):.3e}] FAILED at eps in [{float(e1):.3e}, {float(e2):.3e}]: {info.get('why', info)}",
                None,
            )
        em = (e1 * e2).sqrt() if e1.lower() > 0 else e2 / 2
        stack.append((em, e2, d + 1))
        stack.append((e1, em, d + 1))
    return True, f"{n_ok} ok, {n_skip} skipped", worst


def certify_corner_sweep(L_lo, L_hi, verbose=True):
    """The corner Y <= e^{-L_lo}: L-bands (geometric) x r-boxes, with the Y-free lemma of sibuya_top_w."""
    t0 = time.time()
    n_ok = n_skip = n_band = 0
    worst = None
    L = float(L_lo)
    while L < L_hi:
        L2 = min(L * 1.3, L_hi)
        Lb = arb(
            L
        )  # the lemma is stated for Y <= e^{-L}, so the band's LEFT end is the binding one
        r_hi = float((Lb / JP_MIN).upper())
        stack = [(arb(0), arb(r_hi), 0)]
        while stack:
            r1, r2, d = stack.pop()
            r = arb((float(r1) + float(r2)) / 2, (float(r2) - float(r1)) / 2)
            try:
                st, info = TW.certify_corner(Lb, r, JP_MIN)
            except Exception as ex:  # noqa: BLE001
                st, info = "fail", {"why": repr(ex)[:50]}
            if st == "ok":
                n_ok += 1
                if worst is None or info["slack"].lower() < worst:
                    worst = info["slack"].lower()
                continue
            if st == "skip":
                n_skip += 1
                continue
            if d >= 40:
                return (
                    False,
                    f"corner FAILED at L = {L:.4g}, r in [{float(r1):.3e}, {float(r2):.3e}]: {info.get('why', info)}",
                    None,
                )
            rm = (r1 * r2).sqrt() if r1.lower() > 0 else r2 / 2
            stack.append((rm, r2, d + 1))
            stack.append((r1, rm, d + 1))
        n_band += 1
        if verbose and n_band % 20 == 0:
            print(
                f"  ... corner L = {L:.4g}, {n_band} bands, {n_ok} boxes, worst slack {float(worst):.3e}"
                f"  [{time.time() - t0:.0f} s]",
                flush=True,
            )
        L = L2
    return (
        True,
        f"{n_band} L-bands, {n_ok} boxes ok, {n_skip} skipped, worst slack {float(worst):.4e}  [{time.time() - t0:.0f} s]",
        worst,
    )


def main():
    t0 = time.time()
    print(
        f"Sibuya (3.4), the wedge: theta >= 0.9 and j' >= {JP_MIN}, bands in Y from {float(Y_MAX):.4f} down to 1e-{Y_MIN_EXP}",
        flush=True,
    )
    y_hi = float(Y_MAX)
    y_stop = 10.0**-Y_MIN_EXP
    n_bands = 0
    worst = None
    while y_hi > y_stop:
        y_lo = y_hi / BAND
        ok, msg, w = certify_band(y_lo, y_hi)
        if not ok:
            print(f"FAILED: {msg}")
            sys.exit(1)
        if w is not None and (worst is None or w < worst):
            worst = w
        n_bands += 1
        if n_bands % 200 == 0:
            print(
                f"  ... Y = {y_hi:.3e}, {n_bands} bands, worst margin {float(worst) if worst is not None else float('nan'):.3e}"
                f"  [{time.time() - t0:.0f} s]",
                flush=True,
            )
        y_hi = y_lo
    print(
        f"certified: {n_bands} bands down to Y = {y_stop:.0e}, worst margin {float(worst) if worst is not None else float('nan'):.4e}"
    )
    # NOTE on the corner Y -> 0 at FIXED j': there the leading orders of N g and of the target cancel exactly
    # (measured: the slack encloses 0 to +/- 1e-12 at j' = 20000, L = 10^9), so the crude corner lemma cannot
    # decide it -- that limit is the H-model's territory, certified per index by sibuya_corner_grid.tail_ok.
    # The bands below therefore carry the region as far as they reach, and the residual is stated by N.
    print(
        f"VERDICT (wedge): N g > 1/(3 + 3/N - theta) on {{ theta >= 0.9, j' >= {JP_MIN}, Y >= 1e-{Y_MIN_EXP} }} : True"
        f"  [{time.time() - t0:.0f} s]"
    )


if __name__ == "__main__":
    main()
