"""Figure: the conjecture's constant 4/5 is the odd squares' own relative spread.

Draws two panels:
  left  - the odd squares 1, 9, 25, ... scaled to their mean, showing how they scatter
  right - Var(b)/bbar^2 = 16(m^2-1)/(5(4m^2-1))  converging to 4/5

Writes data/variance_is_the_constant.svg .   Exact rational arithmetic (flint fmpq).
"""

from flint import fmpq
import pathlib


def ratio(m):
    return fmpq(16 * (m * m - 1), 5 * (4 * m * m - 1))


W, H = 900, 380
L, R, T, B = 60, 40, 40, 50
pw = (W - L - R - 60) / 2

out = []
out.append(
    f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="Georgia,serif">'
)
out.append(f'<rect width="{W}" height="{H}" fill="#fcfcfa"/>')

# ---- left panel: the spectrum, scaled by its mean -------------------------------
m = 14
bs = [(2 * k - 1) ** 2 for k in range(1, m + 1)]
mean = sum(bs) / m
x0, y0, ph = L, T + 30, H - T - B - 30
maxv = max(bs) / mean
out.append(
    f'<text x="{x0}" y="{T + 12}" font-size="15" fill="#222">the odd squares, divided by their mean</text>'
)
out.append(
    f'<line x1="{x0}" y1="{y0 + ph}" x2="{x0 + pw}" y2="{y0 + ph}" stroke="#888" stroke-width="1"/>'
)
ymean = y0 + ph - ph * (1 / maxv) * 0.92
out.append(
    f'<line x1="{x0}" y1="{ymean:.1f}" x2="{x0 + pw}" y2="{ymean:.1f}" stroke="#c0392b" stroke-width="1.2" stroke-dasharray="5 4"/>'
)
out.append(
    f'<text x="{x0 + pw - 8}" y="{ymean - 6:.1f}" font-size="12" fill="#c0392b" text-anchor="end">mean</text>'
)
for i, b in enumerate(bs):
    xx = x0 + (i + 0.5) * pw / m
    v = b / mean
    yy = y0 + ph - ph * (v / maxv) * 0.92
    out.append(
        f'<line x1="{xx:.1f}" y1="{y0 + ph}" x2="{xx:.1f}" y2="{yy:.1f}" stroke="#2c3e50" stroke-width="2.4"/>'
    )
    out.append(f'<circle cx="{xx:.1f}" cy="{yy:.1f}" r="3.1" fill="#2c3e50"/>')
out.append(
    f'<text x="{x0}" y="{y0 + ph + 22}" font-size="12" fill="#555">1, 9, 25, ... spread far above and far below their own average</text>'
)

# ---- right panel: the ratio ------------------------------------------------------
x1 = L + pw + 60
ms = [2, 3, 4, 6, 9, 14, 22, 35, 55, 90, 150, 250, 400]
vals = [float(ratio(mm).numer()) / float(ratio(mm).denom()) for mm in ms]
lo, hi = 0.55, 0.815
out.append(
    f'<text x="{x1}" y="{T + 12}" font-size="15" fill="#222">Var(b) / mean(b)^2  =  16(m^2-1) / (5(4m^2-1))</text>'
)
out.append(
    f'<line x1="{x1}" y1="{y0 + ph}" x2="{x1 + pw}" y2="{y0 + ph}" stroke="#888" stroke-width="1"/>'
)
y45 = y0 + ph - ph * (0.8 - lo) / (hi - lo)
out.append(
    f'<line x1="{x1}" y1="{y45:.1f}" x2="{x1 + pw}" y2="{y45:.1f}" stroke="#c0392b" stroke-width="1.6"/>'
)
out.append(
    f'<text x="{x1 + pw - 6}" y="{y45 - 8:.1f}" font-size="14" fill="#c0392b" text-anchor="end">4/5</text>'
)
pts = []
for i, (mm, v) in enumerate(zip(ms, vals)):
    xx = x1 + i * pw / (len(ms) - 1)
    yy = y0 + ph - ph * (v - lo) / (hi - lo)
    pts.append(f"{xx:.1f},{yy:.1f}")
out.append(
    '<polyline points="' + " ".join(pts) + '" fill="none" stroke="#2c3e50" stroke-width="2"/>'
)
for p in pts:
    a, b2 = p.split(",")
    out.append(f'<circle cx="{a}" cy="{b2}" r="3" fill="#2c3e50"/>')
out.append(
    f'<text x="{x1}" y="{y0 + ph + 22}" font-size="12" fill="#555">m = 2 .. 400  (number of distinct roots)</text>'
)

out.append(
    f'<text x="{L}" y="{H - 14}" font-size="13" fill="#333">'
    "Newton's margin is the roots' own spread: the conjecture's constant 4/5 is a property of {1, 9, 25, ...}, not of the inequality."
    "</text>"
)
out.append("</svg>")

pathlib.Path("data").mkdir(exist_ok=True)
pathlib.Path("data/variance_is_the_constant.svg").write_text("\n".join(out), encoding="utf-8")
print("wrote data/variance_is_the_constant.svg")
print("ratio at m = 2, 14, 400 :", ratio(2), ratio(14), ratio(400))
print(
    "as floats               : %.6f %.6f %.6f"
    % tuple(float(ratio(x).numer()) / float(ratio(x).denom()) for x in (2, 14, 400))
)
