#!/usr/bin/env python3
"""plot_headline.py -- the two reader-facing result figures of the papers, drawn from the
canonical CSVs with large type for column width:

  fig_yield_all.{png,pdf}   yield against |F| from 24 to 793 faces: the paper-regime
                            population (results/regime/regime.csv, 200 designs) and the
                            K9c population (results/experiments/k9c/k9c.csv, 400 designs),
                            Wilson 95 % intervals; least-norm projection, published pipeline
                            where it was run, and the range-maximising embedding
  fig_eps_all.{png,pdf}     certified eps_max over the 400 K9c designs against K9, with the
                            medians and the teaser design marked

    arch -arm64 /usr/local/bin/python3 Kirigami/scripts/plot_headline.py [outdir]
"""
import csv, json, math, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = sys.argv[1] if len(sys.argv) > 1 else "results/final/figures"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
                     "legend.fontsize": 9.5, "xtick.labelsize": 10, "ytick.labelsize": 10})
C_PROJ, C_NAT, C_OURS = "#8c8c8c", "#c8102e", "#2b6f8f"
FAM = {"delaunay": "#2b6f8f", "voronoi": "#3f9f6f", "quad_random": "#c98a2b"}

def wilson(k, n, z=1.96):
    if n == 0: return (float("nan"),) * 3
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, c - h, c + h

reg = [r for r in csv.DictReader(open("results/regime/regime.csv")) if r["set"] == "pop"]
k9c = list(csv.DictReader(open("results/experiments/k9c/k9c.csv")))
rows = [dict(F=int(r["F"]), kind=r["kind"], proj=float(r["base_theta"]) > 0,
             ours=float(r["k9c_theta"]) > 0,
             nat=(float(r["natcol_theta"]) > 0) if r["natcol_status"] == "completed" else None)   # their own colouring
        for r in reg]
rows += [dict(F=int(r["F"]), kind=r["kind"], proj=False, ours=float(r["k9c_theta"]) > 0, nat=None)
         for r in k9c]   # K1a: the projection has Theta_max = 0 on all 400 (EXPERIMENTS.md)

edges = [20, 50, 100, 230, 420, 800]
fig, ax = plt.subplots(figsize=(6.4, 3.5))
def series(key, color, label, subset=None, lw=2.2, ls="-", marker="o"):
    xs, ps, lo, hi = [], [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = [r for r in rows if a <= r["F"] < b and (subset is None or r["kind"] == subset)
               and r[key] is not None]
        if len(sel) < 8: continue
        p, l, h = wilson(sum(r[key] for r in sel), len(sel))
        xs.append(math.sqrt(a * b)); ps.append(p); lo.append(l); hi.append(h)
    if not xs: return
    ax.plot(xs, ps, ls, color=color, lw=lw, marker=marker, ms=5, label=label)
    if lw > 1.5:
        ax.fill_between(xs, lo, hi, color=color, alpha=0.15, lw=0)
series("proj", C_PROJ, "least-norm projection (Eq. 6)")
series("nat", C_NAT, "published pipeline (its own colouring), completed cells")
series("ours", C_OURS, "range-maximising embedding (ours)")
for k, c in FAM.items():
    series("ours", c, None, subset=k, lw=1.0, ls="--", marker="")
ax.axvline(100, color="k", lw=0.8, ls=":")
ax.text(97, 0.5, "sizes of prior work's examples", rotation=90, ha="right", va="center", fontsize=9)
ax.set_xscale("log"); ax.set_xlim(22, 800); ax.set_ylim(-0.02, 1.05)
ax.set_xticks([25, 50, 100, 200, 400, 800]); ax.set_xticklabels(["25", "50", "100", "200", "400", "800"])
ax.set_xlabel("faces $|F|$"); ax.set_ylabel("fraction deployable (exact $\\Theta_{\\max}>0$)")
ax.legend(loc="center right", frameon=False)
ax.text(0.02, 0.03, "dashed: ours per family (Delaunay, Voronoi, quad-dominant)", transform=ax.transAxes,
        fontsize=8.5, color="#555")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(out, "fig_yield_all.png"), dpi=220); fig.savefig(os.path.join(out, "fig_yield_all.pdf"))

# certified eps_max
eps_c = np.array([float(r["k9c_eps"]) for r in k9c]); eps_9 = np.array([float(r["k9_eps"]) for r in k9c])
fig, ax = plt.subplots(figsize=(6.4, 2.7))
bins = np.linspace(0, math.pi, 33)
ax.hist(eps_c[eps_c > 0], bins=bins, color=C_OURS, alpha=0.85, label=f"ours, {int((eps_c>0).sum())} of 400 certified")
ax.hist(eps_9[eps_9 > 0], bins=bins, color="#e0a030", alpha=0.85, label=f"proximity objective (K9), {int((eps_9>0).sum())} of 400")
m = np.median(eps_c[eps_c > 0])
ax.axvline(m, color=C_OURS, ls="--", lw=1.2); ax.text(m + 0.04, ax.get_ylim()[1] * 0.92, f"median {m:.2f} rad", color=C_OURS, fontsize=9.5)
try:
    t = json.load(open(os.path.join(out, "teaser_cat.json")))
    ax.axvline(t["ours_eps_max"], color="k", lw=1.0)
    ax.text(t["ours_eps_max"] - 0.04, ax.get_ylim()[1] * 0.6, "Fig. 1 design", rotation=90, ha="right", va="center", fontsize=9)
except FileNotFoundError:
    pass
ax.set_xlabel("certified opening $\\varepsilon_{\\max}$ (rad)"); ax.set_ylabel("designs")
ax.set_xlim(0, math.pi); ax.set_xticks([0, 0.5, 1, 1.5, 2, 2.5, math.pi]); ax.set_xticklabels(["0", "0.5", "1", "1.5", "2", "2.5", "$\\pi$"])
n_not = int((eps_c == 0).sum())
ax.text(0.40, 0.50, f"{n_not} designs not certified\n(exact $\\Theta_{{\\max}}=0$ on 88)", transform=ax.transAxes, ha="left", fontsize=9, color="#555")
ax.legend(frameon=False, loc="upper right")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(out, "fig_eps_all.png"), dpi=220); fig.savefig(os.path.join(out, "fig_eps_all.pdf"))
print("wrote fig_yield_all, fig_eps_all; medians ours %.3f K9 %.3f" % (m, np.median(eps_9[eps_9 > 0])))
