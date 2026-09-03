from flint import fmpq, arb, ctx
from math import comb
import sys

ctx.prec = 2000
K = 42
Z = fmpq(0)
ONE = fmpq(1)


def mul(a, b):
    r = [Z] * K
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if i + j >= K:
                break
            if bj == 0:
                continue
            r[i + j] += ai * bj
    return r


def inv(a):
    r = [Z] * K
    r[0] = 1 / a[0]
    for n in range(1, K):
        s = Z
        for j in range(1, n + 1):
            if j < len(a):
                s += a[j] * r[n - j]
        r[n] = -s / a[0]
    return r


def comp(a, b):
    r = [Z] * K
    p = [Z] * K
    p[0] = ONE
    for n in range(K):
        if a[n] != 0:
            for i in range(K):
                r[i] += a[n] * p[i]
        p = mul(p, b)
    return r


def revert(b):
    q = [Z] * K
    q[1] = 1 / b[1]
    for n in range(2, K):
        t = comp(b, q)
        q[n] -= t[n] / b[1]
    return q


def dropx(a):
    return a[1:] + [Z]


one = [Z] * K
one[0] = ONE
at = [Z] * K
for k in range(K):
    at[k] = fmpq((-1) ** k, 2 * k + 1)
TH = [Z] * K
for k in range(1, K):
    TH[k] = -at[k]
opx = [Z] * K
opx[0] = ONE
opx[1] = ONE
br = [(one[i] - TH[i]) - inv(opx)[i] for i in range(K)]
Hx = [
    2 * inv(dropx(br))[i] - inv(mul(dropx(TH), [one[i] - TH[i] for i in range(K)]))[i]
    for i in range(K)
]
Hth = comp(dropx(Hx), revert(TH))
Hc = [float(Hth[k].numer().str()) / float(Hth[k].denom().str()) for k in range(K - 2)]


def Hval(t):
    s = 0.0
    p = 1.0
    for c in Hc:
        s += c * p
        p *= t
    return s


def Mcurve(M):
    N = 2 * M
    n = N + 1
    e = [fmpq(1)]
    for j in range(1, M + 1):
        c = fmpq((2 * j - 1) ** 2)
        for _ in range(2):
            e = [e[q] + (c * e[q - 1] if q else 0) for q in range(len(e))] + [e[-1] * c]

    def toarb(q):
        return arb(q.numer().str()) / arb(q.denom().str())

    lp = [toarb(e[i]).log() - arb(str(comb(N, i))).log() for i in range(N + 1)]
    pts = []
    for t in range(1, N // 2):
        g = -(lp[t + 1] - 2 * lp[t] + lp[t - 1])
        R = (g).exp()
        pts.append((t / n, float(n * (R - 1))))
    return pts


curves = {M: Mcurve(M) for M in (10, 20, 40, 80)}
Hpts = [(x / 500.0 * 0.5, Hval(x / 500.0 * 0.5)) for x in range(0, 501)]
W, Hh = 980, 560
L, Rm, T, B = 90, 40, 50, 70
x0, x1 = 0.0, 0.5
y0, y1 = 0.70, 2.10


def X(t):
    return L + (t - x0) / (x1 - x0) * (W - L - Rm)


def Y(v):
    return T + (y1 - v) / (y1 - y0) * (Hh - T - B)


def path(pts):
    return "M " + " L ".join("%.2f %.2f" % (X(a), Y(b)) for a, b in pts if y0 <= b <= y1)


cols = {10: "#8ecae6", 20: "#4d96c9", 40: "#2a6f97", 80: "#013a63"}
out = []
out.append(
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" font-family="Georgia, serif">'
    % (W, Hh, W, Hh)
)
out.append('<rect width="%d" height="%d" fill="#fbfaf7"/>' % (W, Hh))
for gv in [0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]:
    out.append(
        '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#e6e2d8" stroke-width="1"/>'
        % (L, Y(gv), W - Rm, Y(gv))
    )
    out.append(
        '<text x="%.1f" y="%.1f" font-size="13" fill="#8a8577" text-anchor="end">%.1f</text>'
        % (L - 10, Y(gv) + 4, gv)
    )
for tv in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
    out.append(
        '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#e6e2d8" stroke-width="1"/>'
        % (X(tv), T, X(tv), Hh - B)
    )
    out.append(
        '<text x="%.1f" y="%.1f" font-size="13" fill="#8a8577" text-anchor="middle">%.1f</text>'
        % (X(tv), Hh - B + 22, tv)
    )
# the 4/5 floor
out.append(
    '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#c1121f" stroke-width="2.5" stroke-dasharray="7 5"/>'
    % (L, Y(0.8), W - Rm, Y(0.8))
)
out.append(
    '<text x="%.1f" y="%.1f" font-size="15" fill="#c1121f" font-style="italic">4/5 &#8212; the floor, never crossed</text>'
    % (L + 12, Y(0.8) - 10)
)
out.append('<path d="%s" fill="none" stroke="#111" stroke-width="3"/>' % path(Hpts))
for M in (10, 20, 40, 80):
    out.append(
        '<path d="%s" fill="none" stroke="%s" stroke-width="1.8" opacity="0.95"/>'
        % (path(curves[M]), cols[M])
    )
    p = [q for q in curves[M] if y0 <= q[1] <= y1]
    if p:
        out.append(
            '<circle cx="%.1f" cy="%.1f" r="3.2" fill="%s"/>' % (X(p[0][0]), Y(p[0][1]), cols[M])
        )
out.append(
    '<text x="%.1f" y="26" font-size="19" fill="#111">The Newton excess of the centred-square spectrum, and the floor it never crosses</text>'
    % L
)
out.append(
    '<text x="%.1f" y="44" font-size="13" fill="#6b665a">M(n,t) = n(p_t&#178;/p_{t-1}p_{t+1} - 1)   against   t/n.  Black: the limit shape H. Blue: n = 21, 41, 81, 161.</text>'
    % L
)
out.append(
    '<text x="%.1f" y="%.1f" font-size="14" fill="#6b665a" text-anchor="middle">t / n</text>'
    % ((L + W - Rm) / 2, Hh - 18)
)
lx = W - Rm - 250
out.append(
    '<text x="%.1f" y="%.1f" font-size="13" fill="#111">limit shape H(&#952;)</text>'
    % (lx + 22, T + 34)
)
out.append(
    '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#111" stroke-width="3"/>'
    % (lx, T + 30, lx + 16, T + 30)
)
for q, M in enumerate((10, 20, 40, 80)):
    out.append(
        '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.8"/>'
        % (lx, T + 52 + q * 19, lx + 16, T + 52 + q * 19, cols[M])
    )
    out.append(
        '<text x="%.1f" y="%.1f" font-size="13" fill="#6b665a">n = %d</text>'
        % (lx + 22, T + 56 + q * 19, 2 * M + 1)
    )
out.append("</svg>")
open("newton_excess_floor.svg", "w").write("\n".join(out))
print("wrote newton_excess_floor.svg")
for M in (10, 20, 40, 80):
    p = curves[M]
    print(
        "   n=%3d   M at t=1 : %.6f    min over the regime : %.6f"
        % (2 * M + 1, p[0][1], min(v for _, v in p))
    )
print("   H(0) = %.10f" % Hval(0.0))
