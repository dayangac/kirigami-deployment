#!/usr/bin/env python3
"""plot_teaser.py -- render results/final/figures/teaser_cat.json (from scripts/teaser_cat.jl)
as the tech report's head figure: planar graph with sigma, the Eq. (6) projection jammed at
0+, and the range-maximising embedding at Theta_max/2.

    arch -arm64 /usr/local/bin/python3 Kirigami/scripts/plot_teaser.py [json] [out.png]
"""
import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch

src = sys.argv[1] if len(sys.argv) > 1 else "results/final/figures/teaser_cat.json"
out = sys.argv[2] if len(sys.argv) > 2 else "results/final/figures/fig_teaser.png"
FS = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0   # font scale (1.4 for the two-column abstract)
d = json.load(open(src))
T, C, A = 11 * FS, 9 * FS, 8.5 * FS
WIDE = FS > 1.0   # compact wording for the two-column abstract

PINK, BLUE = "#f4b6c9", "#b9d4ef"      # the two orientations, as in the 2026 paper
FACE, EDGE = "#c9a9b3", "#2b2b2b"
OURS = "#a9c8b8"
RED = "#c8102e"

def draw(ax, polys, colors, lw=0.6, alpha=1.0, ec=EDGE):
    for p, c in zip(polys, colors):
        ax.add_patch(Polygon(np.asarray(p), closed=True, facecolor=c, edgecolor=ec,
                             linewidth=lw, alpha=alpha, joinstyle="round"))

def frame(ax, polys, pad=0.12):
    P = np.vstack([np.asarray(p) for p in polys])
    lo, hi = P.min(0) - pad, P.max(0) + pad
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1])
    ax.set_aspect("equal"); ax.axis("off")

fig, axs = plt.subplots(1, 4, figsize=(14.0, 4.9) if WIDE else (15.2, 4.5))
sig = d["sigma"]

# (a) the planar graph with its orientation
ax = axs[0]
draw(ax, d["faces_M"], [PINK if s > 0 else BLUE for s in sig])
frame(ax, d["faces_M"])
ax.set_title(r"graph $M$, orientations $\sigma$" if WIDE else r"planar graph $M$, orientations $\sigma$", fontsize=T, style="italic")
ax.text(0.5, -0.02, f"{d['n_faces']} faces, {d['n_split']} split cuts,\nshape space of dim. {d['dim_null']}",
        transform=ax.transAxes, ha="center", va="top", fontsize=C)

# (b) the Eq. (6) projection, opened: faces penetrate
ax = axs[1]
bad = set(d["X0_overlapping"])
cols = [RED if i + 1 in bad else FACE for i in range(len(d["X0_deployed"]))]
order = sorted(range(len(cols)), key=lambda i: cols[i] == RED)   # red on top
draw(ax, [d["X0_deployed"][i] for i in order], [cols[i] for i in order], lw=0.5, alpha=0.95)
frame(ax, d["X0_deployed"])
ax.set_title((r"projection $X_0$ (Eq. 6), $\theta=%.2f$" if WIDE else r"published projection $X_0$ (Eq. 6), $\theta=%.2f$") % d["theta_jam"],
             fontsize=T, style="italic")
ax.text(0.5, -0.02, r"$\Theta_{\max}=0$: %d faces penetrate;" % len(bad)
        + "\nsplit cuts close inward at $0^+$",
        transform=ax.transAxes, ha="center", va="top", fontsize=C, color=RED)

# (c) our embedding, flat
ax = axs[2]
draw(ax, d["ours_flat"], [PINK if s > 0 else BLUE for s in sig])
frame(ax, d["ours_flat"])
ax.set_title(r"our embedding $X$, $\theta=0$" if WIDE else r"range-maximising embedding $X$, $\theta=0$", fontsize=T, style="italic")
ax.text(0.5, -0.02, "same null space; every corner convex,\nevery split cut opens outward",
        transform=ax.transAxes, ha="center", va="top", fontsize=C)

# (d) ours, opened
ax = axs[3]
draw(ax, d["ours_deployed"], [OURS] * len(d["ours_deployed"]), lw=0.5)
frame(ax, d["ours_deployed"])
ax.set_title(r"deployed, $\theta=%.2f$" % d["theta_half"], fontsize=T, style="italic")
ax.text(0.5, -0.02, r"exact $\Theta_{\max}=%.2f$ rad," % d["ours_theta_max"]
        + "\n" + r"certified $\varepsilon_{\max}=%.2f$ rad" % d["ours_eps_max"],
        transform=ax.transAxes, ha="center", va="top", fontsize=C)

# arrows between panels
labels = ["cut, project\n(Eq. 6)", "maximise the\n$0^+$ margin", "rotate every face\nby $\\sigma\\theta/2$"]
for i, label in enumerate(labels):
    x0 = 0.245 + 0.2475 * i
    a = FancyArrowPatch((x0 - 0.012, 0.5), (x0 + 0.012, 0.5),
                        transform=fig.transFigure, arrowstyle="-|>", mutation_scale=16,
                        color="#7a9cc6", linewidth=1.6)
    fig.patches.append(a)
    fig.text(x0, 0.545, label, ha="center", va="bottom", fontsize=A,
             style="italic", color="#3b5a80")

plt.subplots_adjust(left=0.005, right=0.995, top=0.9, bottom=0.11, wspace=0.16)
fig.savefig(out, dpi=220)
fig.savefig(out.rsplit(".", 1)[0] + ".pdf")
print("wrote", out)
