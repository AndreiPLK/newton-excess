"""THE THEOREM, assembled:  M(n,t) > 4/5 for every odd n >= 5 and every 1 <= t < n/2.

    M(n,t) = n ( p_t^2 / (p_{t-1} p_{t+1}) - 1 ),   p_j = e_j / C(N,j),   N = n - 1,
    e_j the elementary symmetric functions of the spectrum { (n-2k)^2 : k = 1..n-1 }
    (the odd squares 1, 9, ..., (n-2)^2, each twice).

The pieces, each an exit-0 script under release/scripts/, and how they fit:

    A. t <= 627, every odd n          the certificate ladder (M-ladder): one exact polynomial
                                      certificate per index (RLC at index i => M(n,i+2) >= M(n,i+1)),
                                      the durable log results/m_ladder_log_2026-09-02_full.txt checked
                                      here by lab/ladder_log_merge.py; sample rungs re-run by ladder_fast.py.
    B. n in {1257, 1259}, t >= 628    exact rational computation here (two sizes below N = 1260).
    C. N >= 1260, t >= 628, theta = t/N <= 1/2:
        C1. t^2/N <= 1                sparse_certificate_full.py   (every t >= 627)
        C2. theta >= 0.05             dense_certificate_a.py        (adaptive sweep in (eps, v))
        C3. theta <= 0.05, t^2/N >= 1 dense_certificate_b.py        (series in V with coefficients in Q[zeta])
    The coverage of C by C1 u C2 u C3 is checked here from the certified ranges of the series
    coefficient khat_1 = theta/V:  with khat_1 in [0.30, 1/3] on V <= 0.17,
        theta <= 0.05  =>  V <= 0.05/0.30 <= 0.17 ;   t^2/N >= 1  =>  zeta = 2 khat_1^2/b <= 2/9 < 0.24 ;
        t >= 628       =>  eta = zeta V = 2 khat_1/t <= 2/(3*628) < 1.07e-3 .
    (Overlaps: C1 reaches b <= (628/629)^2, its own docstring condition, and C3 reaches down past b = 0.926, so they overlap on a band rather than meeting at b = 1 -- the earlier text said "meet at b = 1", corrected after the 2026-09-03 verifier report; C2 and C3 meet at theta = 0.05.)

Run:  uv run python projects/qg-bootstrap/release/scripts/theorem.py [--full]
      --full also re-runs the three certificates (about 5 minutes); without it their recorded
      exit status is trusted and only A, B and the coverage lemma are executed.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from math import comb

from flint import fmpq

HERE = os.path.dirname(os.path.abspath(__file__))


def excess_row(n):
    """M(n,t) for t = 1..(n-1)/2, exact rationals."""
    N = n - 1
    m = N // 2
    e = [fmpq(1)]
    for k in range(1, m + 1):
        c = fmpq((2 * k - 1) ** 2)
        for _ in range(2):
            e = [e[q] + (c * e[q - 1] if q else 0) for q in range(len(e))] + [e[-1] * c]
    p = [e[j] / fmpq(comb(N, j)) for j in range(N + 1)]
    return [fmpq(n) * (p[t] ** 2 / (p[t - 1] * p[t + 1]) - 1) for t in range(1, N // 2 + 1)]


def main():
    full = "--full" in sys.argv
    t0 = time.time()
    print("THE NEWTON EXCESS OF THE CENTRED-SQUARE SPECTRUM:  M(n,t) > 4/5")
    print()
    # ---- B: the two sizes below N = 1260 with t >= 628, exact ---------------------------------
    for n in (1257, 1259):
        row = excess_row(n)
        worst = min(row[627:])  # t = 628 .. (n-1)/2
        print(
            f"  B. n = {n}: min over t >= 628 of M(n,t) = {float(worst):.6f}  > 4/5 : {worst > fmpq(4, 5)}   [{time.time() - t0:.0f} s]"
        )
        assert worst > fmpq(4, 5)
    # ---- coverage lemma: khat_1 in [0.30, 1/3] on V <= 0.17 from the certified series -------------
    sys.path.insert(0, HERE)
    import dense_certificate_b as DB  # noqa: E402
    from flint import arb  # noqa: E402

    kt = DB.cumulant_series(2)
    khat1 = kt[0].divV()
    # enclose khat_1 on V in (0, VMAX] by 16 sub-boxes (zeta <= ZMAX)
    lo_all, hi_all = None, None
    for i in range(16):
        v1, v2 = float(DB.VMAX) * i / 16, float(DB.VMAX) * (i + 1) / 16
        k1 = khat1.evaluate(
            arb((v1 + v2) / 2, (v2 - v1) / 2), arb(float(DB.ZMAX) / 2, float(DB.ZMAX) / 2)
        )
        lo_all = k1.lower() if lo_all is None else min(lo_all, k1.lower())
        hi_all = k1.upper() if hi_all is None else max(hi_all, k1.upper())
    print(
        f"  coverage: khat_1 = theta/V on V <= {DB.VMAX}, zeta <= {DB.ZMAX}: in [{float(lo_all):.5f}, {float(hi_all):.5f}]"
    )
    imp1 = arb("0.05") / arb(lo_all)  # theta <= 0.05 => V <= this
    imp2 = 2 * arb(hi_all) ** 2  # b >= 1 => zeta <= this
    imp3 = 2 * arb(hi_all) / 628  # t >= 628 => eta <= this
    print(
        f"     theta <= 0.05 => V <= {float(imp1):.5f} <= {float(DB.VMAX)};  b >= 1 => zeta <= {float(imp2):.5f} <= {float(DB.ZMAX)};"
    )
    print(
        f"     t >= 628 => eta <= {float(imp3):.6f} <= {float(DB.ETAMAX)}.   Coverage of C by C1 u C2 u C3: "
        f"{(imp1 <= DB.toarb(DB.VMAX)) and (imp2 <= DB.toarb(DB.ZMAX)) and (imp3 <= DB.toarb(DB.ETAMAX))}"
    )
    assert (
        (imp1 <= DB.toarb(DB.VMAX))
        and (imp2 <= DB.toarb(DB.ZMAX))
        and (imp3 <= DB.toarb(DB.ETAMAX))
    )
    # ---- A: the ladder artifact, checked, plus a sample re-verification ----------------------------
    # The durable log results/m_ladder_log_2026-09-02_full.txt is produced and CHECKED by
    # lab/ladder_log_merge.py (every rung 1..627 present once, every shift equal to its prediction,
    # degree 8i+10, no None); the check is re-run here and its exit status is part of the theorem.
    r = subprocess.run(
        [sys.executable, os.path.join(HERE, "..", "..", "lab", "ladder_log_merge.py"), "627"],
        capture_output=True,
        text=True,
        cwd=HERE,
    )
    print(
        "  A. the ladder t <= 627:",
        r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr.strip()[-200:],
    )
    assert r.returncode == 0, "the ladder log is incomplete or inconsistent"
    print("     sample re-verification of rungs 1..12 by ladder_fast.py:")
    r = subprocess.run(
        [sys.executable, os.path.join(HERE, "ladder_fast.py"), "12", "1"],
        capture_output=True,
        text=True,
        cwd=HERE,
    )
    tail = "\n".join("       " + line for line in r.stdout.strip().splitlines()[-4:])
    print(tail)
    assert r.returncode == 0
    # ---- C: the three certificates ---------------------------------------------------------------
    if full:
        for script in (
            "sparse_certificate_full.py",
            "dense_certificate_a.py",
            "dense_certificate_b.py",
        ):
            tic = time.time()
            r = subprocess.run(
                [sys.executable, os.path.join(HERE, script)],
                capture_output=True,
                text=True,
                cwd=HERE,
            )
            verdict = [line for line in r.stdout.splitlines() if line.startswith("VERDICT")]
            print(
                f"  C. {script}: exit {r.returncode}  {verdict[-1] if verdict else ''}  [{time.time() - tic:.0f} s]"
            )
            assert r.returncode == 0
    else:
        print(
            "  C. certificates not re-run (pass --full); recorded 3 Sept 2026: sparse_certificate_full.py exit 0 (3 s),"
        )
        print(
            "     dense_certificate_a.py exit 0 (288 s, 930 boxes, results/dense_certificate_a_rerun_2026-09-03.txt),"
        )
        print(
            "     dense_certificate_b.py exit 0 (41 s, results/dense_certificate_b_rerun2_2026-09-03.txt); all three independently validated (results/VALIDATION_*_2026-09-02.md; the dense (b) report carries two addenda, the last of which passes all seven items)."
        )
    print()
    print(
        f"THEOREM: M(n,t) > 4/5 for every odd n >= 5 and every t < n/2 -- all pieces in place.  [{time.time() - t0:.0f} s]"
    )


if __name__ == "__main__":
    main()
