"""The theorem as terrain, in the lab's approved 3D style: the excess surface M(n,t) from EXACT values, the
floor 4/5 as a plane it never touches, the lowest ridge (t = 1) as the yellow lane, and the diamond at the
limit the theorem proves is never reached.  Rendered with plotly + kaleido at 2x.

Run:  uv run python projects/qg-bootstrap/release/scripts/fig3d_theorem_plotly.py [out.png] [n_max]
"""
from __future__ import annotations

import os
import sys

import numpy as np
import plotly.graph_objects as go

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
_argv = sys.argv
sys.argv = [_argv[0]]
from fig3d_theorem import excess_row  # noqa: E402

sys.argv = _argv
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "data", "theorem_3d.png")
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 121

ns = list(range(5, NMAX + 1, 2))
tmax = (NMAX - 1) // 2
Z = np.full((tmax, len(ns)), np.nan)
for j, n in enumerate(ns):
    row = excess_row(n)
    for i, v in enumerate(row):
        Z[i, j] = float(v)
T = np.arange(1, tmax + 1)
N = np.array(ns)
zmin = np.nanmin(Z)
zmax = np.nanmax(Z)
print(f"grid {Z.shape}, M in [{zmin:.4f}, {zmax:.4f}], floor 4/5 = 0.8")

# compress the dynamics so the floor region and the high ridge both read
def comp(z):
    return np.log10(z)

floor = 0.8
surf = go.Surface(
    x=N, y=T, z=comp(Z),
    colorscale=[[0.0, "#ffb347"], [0.18, "#e8743b"], [0.4, "#a8437a"], [0.7, "#1682a8"], [1.0, "#0a3450"]],
    cmin=comp(floor), cmax=comp(zmax),
    showscale=False, opacity=1.0,
    lighting=dict(ambient=0.72, diffuse=0.45, specular=0.06, roughness=0.95),
    hoverinfo="skip",
)
# the floor: a plane at 4/5 under the whole grid
FX, FY = np.meshgrid(N, T)
plane = go.Surface(
    x=N, y=T, z=np.full_like(Z, comp(floor)),
    colorscale=[[0, "#ff5f95"], [1, "#ff5f95"]], showscale=False, opacity=0.55, hoverinfo="skip",
    lighting=dict(ambient=0.9, diffuse=0.2, specular=0.0),
)
# the yellow lane: t = 1, the lowest ridge, creeping down to the floor
lane = go.Scatter3d(
    x=N, y=np.ones_like(N), z=comp(Z[0, :]) + 0.004, mode="lines",
    line=dict(color="#f9f871", width=9), hoverinfo="skip",
)
diamond = go.Scatter3d(
    x=[N[-1]], y=[1], z=[comp(Z[0, -1]) + 0.006], mode="markers",
    marker=dict(size=9, color="white", symbol="diamond", line=dict(color="#f9f871", width=2)), hoverinfo="skip",
)
labels = go.Scatter3d(
    x=[N[len(N) // 2], N[-1] * 0.72, N[-1] * 0.82],
    y=[tmax * 0.55, tmax * 0.78, 4],
    z=[comp(zmax) * 0.98, comp(floor) - 0.03, comp(Z[0, -1]) + 0.16],
    mode="text",
    text=["ALIVE  —  ABOVE THE FLOOR", "THE FLOOR  4/5", "4/5 at infinity — never reached"],
    textfont=dict(color=["#9ff5ff", "#ff5f95", "#ffe94a"], size=[20, 20, 16], family="Segoe UI, Arial"),
    hoverinfo="skip",
)
fig = go.Figure([surf, plane, lane, diamond, labels])
fig.update_layout(
    paper_bgcolor="#05030c", plot_bgcolor="#05030c", showlegend=False,
    width=1400, height=900, margin=dict(l=0, r=0, t=90, b=0),
    title=dict(text="<b>FLOOR</b>", x=0.04, y=0.95, font=dict(size=36, color="#ffffff", family="Segoe UI, Arial")),
    annotations=[dict(x=0.04, y=0.88, xref="paper", yref="paper", showarrow=False, align="left",
                      text="height = Newton excess M(n,t), log scale  ·  pink plane = the floor 4/5  ·  yellow lane = the lowest ridge t = 1",
                      font=dict(size=14, color="#c9d3e6", family="Segoe UI, Arial"))],
    scene=dict(
        xaxis=dict(title="n  (size of the list)", color="#9aa6be", gridcolor="#1c2333", backgroundcolor="#05030c", showbackground=True),
        yaxis=dict(title="t  (which step of the staircase)", color="#9aa6be", gridcolor="#1c2333", backgroundcolor="#05030c", showbackground=True),
        zaxis=dict(title="room above Newton", color="#9aa6be", gridcolor="#1c2333", backgroundcolor="#05030c", showbackground=True,
                   tickvals=[comp(0.8), 0, 1, 2], ticktext=["4/5", "1", "10", "100"]),
        camera=dict(eye=dict(x=1.85, y=-1.55, z=0.75), center=dict(x=0, y=0, z=-0.12)),
        aspectratio=dict(x=1.5, y=1.0, z=0.6),
    ),
)
fig.write_image(OUT, scale=2)
print("wrote", OUT)
