"""The 3D picture of Sibuya's conjecture: the Newton excess of the Stirling numbers of the first kind and its floor 1/3.

Exact values (fmpq) on a grid of odd n and t < n/2, rendered as a hand-written isometric SVG (the repository's
figures carry no plotting dependency).  The reading, in the founder's standard: the camera looks at the wall,
the zones are named in words, the yellow path is where the surface comes closest to the floor, and the diamond
sits at the limit that the theorem proves is never reached.

Run:  uv run python projects/qg-bootstrap/release/scripts/fig3d_theorem.py [out.svg]
"""

from __future__ import annotations

import sys
from math import comb, log

from flint import fmpq

OUT = sys.argv[1] if len(sys.argv) > 1 else "projects/qg-bootstrap/release/data/sibuya_3d.svg"
W, HT = 1180, 760


def excess_row(n):
    """M(n,t) = n (p_t^2/(p_{t-1}p_{t+1}) - 1) for the spectrum 1, 2, ..., n-1, exact rationals.

    p_j = e_j(1..n-1)/C(n-1, j) are the normalised unsigned Stirling numbers of the first kind, and Sibuya's
    conjecture (1988, eq. 3.4) says the excess never falls below the floor n/(3n-t) -> 1/3.
    """
    N = n - 1
    e = [fmpq(1)]
    for k in range(1, N + 1):
        c = fmpq(k)
        e = [e[q] + (c * e[q - 1] if q else 0) for q in range(len(e))] + [e[-1] * c]
    p = [e[j] / fmpq(comb(N, j)) for j in range(N + 1)]
    return [float(fmpq(n) * (p[t] ** 2 / (p[t - 1] * p[t + 1]) - 1)) for t in range(1, N // 2 + 1)]


# ---------------------------------------------------------------- projection
def proj(u, v, h):
    """isometric: u across (theta), v into the page (log n), h up (the excess)."""
    x = 210 + 640 * u + 250 * v
    y = 600 - 300 * v - 300 * (h - 1 / 3) / 0.40
    return x, y


def path(points, close=False):
    d = "M " + " L ".join(f"{x:.1f} {y:.1f}" for x, y in points)
    return d + (" Z" if close else "")


def main():
    ns = [5, 9, 17, 33, 65, 129, 257, 513]
    rows = []
    for n in ns:  # b_k = k: the unsigned Stirling numbers of the first kind
        r = excess_row(n)
        v = log(n / ns[0]) / log(ns[-1] / ns[0])
        rows.append((n, v, r))
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {HT}" width="{W}" height="{HT}" '
        f'font-family="Georgia, serif">',
        f'<rect width="{W}" height="{HT}" fill="#fbfaf7"/>',
    ]
    # ---- the floor: the plane M = 1/3, drawn as the wall the surface never crosses
    floor = [proj(0, 0, 1 / 3), proj(1, 0, 1 / 3), proj(1, 1, 1 / 3), proj(0, 1, 1 / 3)]
    out.append(
        f'<path d="{path(floor, True)}" fill="#f0d9c8" fill-opacity="0.55" stroke="#c9855a" stroke-width="1.5"/>'
    )
    out.append(
        f'<text x="{proj(0.52, 0.5, 1 / 3)[0]:.0f}" y="{proj(0.52, 0.5, 1 / 3)[1] + 26:.0f}" font-size="20" '
        f'fill="#b4653a" text-anchor="middle">the floor  M = 1/3  —  never touched</text>'
    )
    # ---- the surface, drawn back to front as filled ribbons between consecutive n
    for i in range(len(rows) - 1, 0, -1):
        n1, v1, r1 = rows[i - 1]
        n2, v2, r2 = rows[i]
        m = 60
        top, bot = [], []
        for s in range(m + 1):
            u = s / m
            j1 = min(int(u * (len(r1) - 1)), len(r1) - 1)
            j2 = min(int(u * (len(r2) - 1)), len(r2) - 1)
            top.append(proj(u, v2, min(r2[j2], 0.72)))
            bot.append(proj(u, v1, min(r1[j1], 0.72)))
        band = top + bot[::-1]
        shade = 0.35 + 0.5 * (i / len(rows))
        col = f"rgb({int(120 + 90 * shade)},{int(150 + 70 * shade)},{int(185 + 55 * shade)})"
        out.append(
            f'<path d="{path(band, True)}" fill="{col}" fill-opacity="0.9" stroke="#54708c" stroke-width="0.7"/>'
        )
    # ---- the ridge lines n = const, and the labels
    for n, v, r in rows:
        pts = [
            proj(s / 60, v, min(r[min(int((s / 60) * (len(r) - 1)), len(r) - 1)], 0.72))
            for s in range(61)
        ]
        out.append(f'<path d="{path(pts)}" fill="none" stroke="#33506b" stroke-width="1.2"/>')
        x, y = pts[-1]
        out.append(
            f'<text x="{x + 8:.0f}" y="{y + 4:.0f}" font-size="14" fill="#33506b">n = {n}</text>'
        )
    # ---- the yellow path: t = 1, where the surface dives closest to the floor
    yel = [proj(0, v, min(r[0], 0.72)) for _, v, r in rows]
    out.append(
        f'<path d="{path(yel)}" fill="none" stroke="#e0a300" stroke-width="4.5" stroke-linecap="round"/>'
    )
    for x, y in yel:
        out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="#e0a300"/>')
    dx, dy = proj(0, 1.18, 1 / 3)
    out.append(
        f'<path d="M {dx:.1f} {dy - 11:.1f} L {dx + 11:.1f} {dy:.1f} L {dx:.1f} {dy + 11:.1f} L {dx - 11:.1f} {dy:.1f} Z" '
        f'fill="#e0a300" stroke="#8a6400" stroke-width="1.5"/>'
    )
    out.append(
        f'<text x="{dx - 18:.0f}" y="{dy + 6:.0f}" font-size="17" fill="#8a6400" text-anchor="end">'
        f"the limit 1/3 = Var/mean² of the roots</text>"
    )
    out.append(
        f'<text x="{yel[3][0] + 14:.0f}" y="{yel[3][1] - 14:.0f}" font-size="17" fill="#8a6400">'
        f"the lowest edge: t = 1, sinking towards 1/3</text>"
    )
    # ---- zones, in words
    out.append(
        f'<text x="{proj(0.75, 0.15, 0.68)[0]:.0f}" y="{proj(0.75, 0.15, 0.68)[1]:.0f}" font-size="17" '
        f'fill="#33506b">the middle of the row rises to 0.63 and stays there</text>'
    )
    # ---- axes in words
    out.append(
        '<text x="200" y="712" font-size="18" fill="#5c584d">left to right: the coefficient index, from the '
        "first one to the middle of the row</text>"
    )
    out.append(
        '<text x="820" y="640" font-size="18" fill="#5c584d" text-anchor="middle">into the page: the size n, '
        "doubling each step from 5 to 513</text>"
    )
    out.append(
        '<text x="60" y="150" font-size="18" fill="#5c584d">up: the Newton excess M(n,t)</text>'
    )
    out.append(
        '<text x="60" y="58" font-size="26" fill="#2c2a25">Newton\'s inequality has a floor, and for the Stirling numbers it is 1/3</text>'
    )
    out.append(
        '<text x="60" y="88" font-size="17" fill="#6b6659">M(n,t) = n (p_t² / p_{t-1}p_{t+1} − 1) for the '
        "spectrum 1, 2, ..., n-1 (Stirling numbers of the first kind); every value shown is exact arithmetic</text>"
    )
    out.append("</svg>")
    open(OUT, "w", encoding="utf-8").write("\n".join(out))
    print(
        f"wrote {OUT}  ({len(rows)} rows, n = {ns[0]}..{ns[-1]}; min at t=1: {min(r[0] for _, _, r in rows):.4f})"
    )


if __name__ == "__main__":
    main()
