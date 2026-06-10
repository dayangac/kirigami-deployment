"""Compact single-column yield figure for the SIGGRAPH poster paper.

Reuses code/scripts/plot_final.py's data loading, |F| binning and Wilson interval
verbatim (imported, not re-implemented). Difference from fig_yield.png: only the
baseline and the K9c range-maximising arm are drawn, only the exact Theta_max > 0
metric, and the two orientation rules are POOLED per (family, bin) so that one
panel fits a single column.  Run from the repo root with
  arch -arm64 /usr/local/bin/python3 docs/poster/plot_poster_yield.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "code", "scripts"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plot_final as pf

rows = pf.load_k9c_merged()
edges = [100, 230, 420, 800]
labels = ["|F|<230", "230-419", "|F|>=420"]

def bin_of(F):
    for i in range(len(edges) - 1):
        if edges[i] <= F < edges[i + 1]:
            return i
    return len(edges) - 2

fig, ax = plt.subplots(figsize=(3.32, 1.95))
colors = {"delaunay": "#2b5fb5", "voronoi": "#63a0d8", "quad_random": "#8fbf6f"}
markers = {"delaunay": "o", "voronoi": "s", "quad_random": "^"}
offs = {"delaunay": -0.09, "voronoi": 0.0, "quad_random": 0.09}
total_n = 0
for fam in pf.FAMILIES:
    sub = [r for r in rows if r["kind"] == fam]
    total_n += len(sub)
    xs, ys, lo, hi = [], [], [], []
    for b in range(len(edges) - 1):
        rs = [r for r in sub if bin_of(int(r["F"])) == b]
        if not rs:
            continue
        k = sum(1 for r in rs if float(r["k9c_theta"]) > 1e-9)
        p, l, h = pf.wilson(k, len(rs))
        xs.append(b + offs[fam]); ys.append(p)
        lo.append(max(0.0, p - l)); hi.append(max(0.0, h - p))
        print(fam, labels[b], "%d/%d" % (k, len(rs)))
    ax.errorbar(xs, ys, yerr=[lo, hi], fmt=markers[fam], ms=3.6, color=colors[fam],
                linestyle="-", linewidth=1.1, capsize=2, elinewidth=0.9,
                label=pf.FAMILY_LABEL[fam])
ax.plot([0, 1, 2], [0, 0, 0], "x--", color="0.45", ms=4, linewidth=1.0,
        label="Eq.(6)+repairs (all)")
ax.set_xticks(range(3)); ax.set_xticklabels(labels, fontsize=7.5)
ax.set_ylim(-0.04, 1.08); ax.tick_params(labelsize=7.5)
ax.set_xlabel("face count $|F|$", fontsize=8)
ax.set_ylabel("fraction deployable", fontsize=8)
ax.legend(fontsize=6.6, frameon=False, ncol=2, loc="lower left",
          handlelength=1.6, columnspacing=1.0, borderaxespad=0.2)
fig.tight_layout(pad=0.25)
out = "docs/poster/figs/fig_yield_poster.png"
fig.savefig(out, dpi=400)
print("wrote", out, "N =", total_n)
