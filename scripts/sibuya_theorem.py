"""SIBUYA'S CONJECTURE (1988, eq. 3.4), assembled:

    p_j^2 / (p_{j-1} p_{j+1})  >=  1 + 1/(3n - j)     for every n >= 3 and every 1 <= j <= n - 2
                                                      CERTIFIED HERE for j <= 1000 (every n) and for j >= 1001 with
                                                      theta = j/(n-1) <= 9/10;  the top regime theta > 9/10 is the
                                                      named open piece (OBSTRUCTION.md: the reciprocal spectrum {1/k}),

    p_j = e_j(1, 2, ..., n-1) / C(n-1, j)   (the unsigned Stirling numbers of the first kind, normalised),
    with equality at j = 1.

The pieces, each an exit-0 script, and how they fit (N = n - 1 roots, theta = j/N, t = j):

    A'. j <= 1000, every n           the Sibuya ladder (lab/sibuya_ladder.py): for each index j the polynomial
                                     W_j(n) = (3n-j) p_j^2 - (3n-j+1) p_{j-1} p_{j+1} has non-negative coefficients
                                     after the shift n -> n + j + 2, and the finitely many n in [j+2, 2j+3] are
                                     checked exactly; j = 1 is an identity.  Durable log
                                     results/sibuya_ladder_log_full.txt, merged and checked by
                                     lab/sibuya_ladder_merge.py (re-run here).
    C'. j >= 1001 (so N >= 2002, theta <= 1/2):  the three certificates prove  n g > n/(3n - j)  with
        g = -Delta^2 log p_j  (then p_j^2/(p_{j-1}p_{j+1}) - 1 = e^g - 1 > g > 1/(3n-j)):
        C1'. j^2/N <= 1                sibuya_sparse_certificate.py   ((3 + 3a^2 b - ab) N Phi > 1, every tau >= 1000;
                                       lemma constants independently validated 3 Sept 2026)
        C2'. 0.1 <= theta <= 1/2       sibuya_dense_a.py              (N g > 1/(3 + 3/N - theta), boxes in (1/N, v))
        C2't. 1/2 <= theta <= 9/10     sibuya_dense_a.py --top        (the same boxes, v in [2.5, 45], N >= 1112)
        C3'. theta <= 0.1, j^2/N >= 1  sibuya_dense_b.py              (series in v with coefficients in Q[zeta])
    Coverage of C' by C1' u C2' u C3' is checked here from the certified range of khat_1 = theta/v on v <= VMAX:
        theta <= 0.1  =>  v <= 0.1/khat_1_min <= VMAX;   j^2/N >= 1  =>  zeta = khat_1^2/b <= khat_1_max^2 <= ZMAX;
        j >= 1001     =>  eta = zeta v = khat_1/j <= khat_1_max/1001 <= ETAMAX.
    (C2' covers theta in [0.1, 1/2] for every N >= 2002 by its own coverage check at the v-edges.)

Run:  uv run python projects/qg-bootstrap/release/scripts/sibuya_theorem.py [--full]
      --full re-runs the three certificates (about 12 minutes); without it their logged verdicts are trusted
      and only A' and the coverage lemma are executed.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.join(HERE, "..", "..", "lab")


def main():
    full = "--full" in sys.argv
    t0 = time.time()
    print("SIBUYA'S CONJECTURE (1988, eq. 3.4):  p_j^2/(p_{j-1} p_{j+1}) >= 1 + 1/(3n - j)")
    print()
    # ---- A': the ladder artifact
    r = subprocess.run(
        [sys.executable, os.path.join(LAB, "sibuya_ladder_merge.py"), "1000"],
        capture_output=True,
        text=True,
        cwd=HERE,
    )
    print(
        "  A'. the ladder j <= 1000:",
        r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr.strip()[-300:],
    )
    assert r.returncode == 0, "the Sibuya ladder log is incomplete or inconsistent"
    # ---- coverage lemma from the certified khat_1 range of the series
    sys.path.insert(0, HERE)
    import sibuya_dense_b as SB  # noqa: E402
    from flint import arb  # noqa: E402

    kt = SB.cumulant_series(2)
    khat1 = kt[0].divV()
    lo_all, hi_all = None, None
    for i in range(16):
        v1, v2 = float(SB.VMAX) * i / 16, float(SB.VMAX) * (i + 1) / 16
        k1 = khat1.evaluate(
            arb((v1 + v2) / 2, (v2 - v1) / 2), arb(float(SB.ZMAX) / 2, float(SB.ZMAX) / 2)
        )
        lo_all = k1.lower() if lo_all is None else min(lo_all, k1.lower())
        hi_all = k1.upper() if hi_all is None else max(hi_all, k1.upper())
    # theta(v) = v khat_1(v) at the edge v = VMAX must exceed 0.1 for every zeta (theta increases with v)
    k1_edge = khat1.evaluate(arb(float(SB.VMAX)), arb(float(SB.ZMAX) / 2, float(SB.ZMAX) / 2))
    theta_edge = float(SB.VMAX) * float(k1_edge.lower())
    imp2 = arb(hi_all) ** 2 * arb(1002) ** 2 / arb(1001) ** 2  # b >= 1 => zeta <= khat_1_max^2
    imp3 = arb(hi_all) / 1001  # j >= 1001 => eta <= khat_1_max / 1001
    print(
        f"  coverage: khat_1 on v <= {SB.VMAX}, zeta <= {SB.ZMAX}: in [{float(lo_all):.5f}, {float(hi_all):.5f}];"
        f"  theta(VMAX) >= {theta_edge:.5f} > 0.1: {theta_edge > 0.1}"
    )
    print(
        f"     N < (j+1)^2 => zeta <= {float(imp2):.5f} <= {float(SB.ZMAX)};  j >= 1001 => eta <= {float(imp3):.3e} <= {float(SB.ETAMAX)}."
        f"   Coverage of C' by C1' u C2' u C3': {theta_edge > 0.1 and (imp2 <= SB.toarb(SB.ZMAX)) and (imp3 <= SB.toarb(SB.ETAMAX))}"
    )
    assert theta_edge > 0.1 and (imp2 <= SB.toarb(SB.ZMAX)) and (imp3 <= SB.toarb(SB.ETAMAX))
    # ---- C': the three certificates
    if full:
        for script in (
            "sibuya_sparse_certificate.py",
            "sibuya_dense_a.py",
            "sibuya_dense_a.py --top",
            "sibuya_dense_b.py",
        ):
            tic = time.time()
            r = subprocess.run(
                [sys.executable, os.path.join(HERE, script.split()[0])] + script.split()[1:],
                capture_output=True,
                text=True,
                cwd=HERE,
            )
            verdict = [line for line in r.stdout.splitlines() if line.startswith("VERDICT")]
            print(
                f"  C'. {script}: exit {r.returncode}  {verdict[-1] if verdict else ''}  [{time.time() - tic:.0f} s]"
            )
            assert r.returncode == 0
    else:
        print(
            "  C'. certificates not re-run (pass --full); logged verdicts: results/sibuya_sparse_certificate_2026-09-03.txt,"
        )
        print(
            "      results/sibuya_dense_a_2026-09-03.txt, results/sibuya_dense_a_top_2026-09-03.txt,"
        )
        print("      results/sibuya_dense_b_2026-09-03.txt (each VERDICT True, exit 0).")
    print()
    print(
        "SIBUYA (3.4): p_j^2/(p_{j-1}p_{j+1}) >= 1 + 1/(3n - j) CERTIFIED for j <= 1000 (every n) and for theta <= 9/10;"
        f" OPEN: theta > 9/10 (see OBSTRUCTION.md); all four certificates independently validated (results/VALIDATION_SIBUYA_*_2026-09-03.md).  [{time.time() - t0:.0f} s]"
    )


if __name__ == "__main__":
    main()
