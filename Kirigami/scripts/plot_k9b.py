#!/usr/bin/env python3
"""K9b figures.

  k9b_gallery.png -- >= 20 deployed random graphs, each shown closed (theta = 0) and at
                     theta = Theta_max / 2, from the JSON dumped by kill_k9b.
  k9b_hist.png    -- the certified eps_max distribution of K9b against the baselines,
                     every one of which is identically 0.

Run with: arch -arm64 /usr/local/bin/python3 code/scripts/plot_k9b.py [results/kill/k9b]
"""
import csv
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

D = sys.argv[1] if len(sys.argv) > 1 else "results/kill/k9b"
G = os.path.join(D, "gallery")

CLOSED = "#BDD2FD"
OPEN = "#5B8FF9"


def load_index():
    rows = []
    for p in sorted(glob.glob(os.path.join(G, "index_*.json"))):
        rows.extend(json.load(open(p)))
    # one entry per graph: the sigma with the larger Theta_max (the honest best-of-2)
    best = {}
    for r in rows:
        k = (r["kind"], r["id"])
        if k not in best or r["theta_max"] > best[k]["theta_max"]:
            best[k] = r
    out = sorted(best.values(), key=lambda r: -r["theta_max"])
    return out, rows


def draw(ax, path, color, title):
    d = json.load(open(path))
    V = d["vertices"]
    patches = [Polygon([V[k] for k in f], closed=True) for f in d["faces"]]
    ax.add_collection(PatchCollection(patches, facecolor=color, edgecolor="0.3",
                                      linewidths=0.15))
    xs = [p[0] for p in V]
    ys = [p[1] for p in V]
    mx = 0.04 * max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
    ax.set_xlim(min(xs) - mx, max(xs) + mx)
    ax.set_ylim(min(ys) - mx, max(ys) + mx)
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=6.5, pad=1.5)


def gallery(n_want=20):
    entries, allrows = load_index()
    entries = [e for e in entries
               if os.path.exists(os.path.join(G, e["tag"] + "_closed.json"))
               and os.path.exists(os.path.join(G, e["tag"] + "_half.json"))]
    if not entries:
        print("gallery: nothing dumped")
        return 0
    entries = entries[:max(n_want, min(len(entries), 24))]
    ncol = 8  # 4 designs per row, each as a (closed, half) pair
    nrow = (len(entries) + 3) // 4
    fig, axes = plt.subplots(nrow, ncol, figsize=(1.55 * ncol, 1.75 * nrow))
    axes = axes.reshape(nrow, ncol)
    for a in axes.ravel():
        a.axis("off")
    for i, e in enumerate(entries):
        r, cpair = divmod(i, 4)
        a0 = axes[r][2 * cpair]
        a1 = axes[r][2 * cpair + 1]
        draw(a0, os.path.join(G, e["tag"] + "_closed.json"), CLOSED,
             "%s %d %s\nF=%d  closed" % (e["kind"][:4], e["id"],
                                         "mc" if e["sigma"].endswith("mc") else "def",
                                         e["F"]))
        draw(a1, os.path.join(G, e["tag"] + "_half.json"), OPEN,
             r"$\theta=\Theta_{max}/2$" + "\n%.3f rad" % (0.5 * e["theta_max"]))
    fig.suptitle("K9b: %d random graphs that deploy under the convexity + split-inward "
                 "constrained embedding\n(every baseline on the same population is "
                 "$\\Theta_{max}=0$ on all 400 designs)" % len(entries), fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    out = os.path.join(D, "k9b_gallery.png")
    fig.savefig(out, dpi=190)
    plt.close(fig)
    print("wrote %s with %d designs" % (out, len(entries)))
    return len(entries)


def union_positives(rows):
    """Per design, the better of K9's solver configuration and K9b's sweep. Both are
    exact-verified positives under the same certificate on the same 400 designs, so the
    per-design maximum is a legitimate multi-start count -- and it is reported next to
    the two runs, never in place of either."""
    path = "results/kill/k9/k9.csv"
    if not os.path.exists(path):
        return sum(1 for r in rows if float(r["theta_ref9"]) > 1e-9)
    K = {(r["kind"], r["id"], r["sigma"]): r for r in csv.DictReader(open(path))}
    n = 0
    for r in rows:
        k = (r["kind"], r["id"], r["sigma"])
        best = float(r["theta_exact"])
        if k in K:
            best = max(best, float(K[k]["b_theta_exact"]))
        if best > 1e-9:
            n += 1
    return n


def hist():
    path = os.path.join(D, "k9b.csv")
    if not os.path.exists(path):
        print("no k9b.csv")
        return
    rows = list(csv.DictReader(open(path)))
    eps = [float(r["eps_max"]) for r in rows if float(r["eps_max"]) > 0]
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.8))

    ax[0].hist(eps, bins=24, color="#5B8FF9", edgecolor="0.25", linewidth=0.5)
    ax[0].axvline(0.1, color="#E8684A", ls="--", lw=1.2)
    ax[0].text(0.105, ax[0].get_ylim()[1] * 0.92, "0.1 rad target", color="#E8684A",
               fontsize=8)
    ax[0].set_xlabel(r"certified $\varepsilon_{max}$ (rad)")
    ax[0].set_ylabel("designs")
    ax[0].set_title("K9b: %d certified designs of %d" % (len(eps), len(rows)), fontsize=9)

    # Every baseline on this exact population is identically zero.
    names = ["Eq. (6)\n(K1a)", r"$\sigma_{mc}$" + "\n(K5)", r"$\sigma_{def}$" + "\n(K5)",
             "0+ repair\n(K6)", "convexity\nonly (K9a)", "K9 (b)", "K9b\nsweep",
             "best of\nboth"]
    k9 = 36
    # A positive needs BOTH the exact scan and the referee: delaunay 40 sigma_mc has
    # bisection(1e-9) = 2.07e-7 rad, below the bisection's own grid resolution, with the
    # exact scan at 0. It is bisection noise, not a design that deploys.
    k9b = sum(1 for r in rows
              if float(r["theta_exact"]) > 1e-9 and float(r["theta_ref9"]) > 1e-9)
    union = union_positives(rows)
    vals = [0, 0, 0, 0, 0, k9, k9b, union]
    cols = ["0.7"] * 5 + ["#9BB7F0", "#9BB7F0", "#5B8FF9"]
    ax[1].bar(range(len(vals)), vals, color=cols, edgecolor="0.25", linewidth=0.5)
    ax[1].axhline(40, color="#E8684A", ls="--", lw=1.2)
    ax[1].text(0.05, 41, "PASS bar 40 / 400", color="#E8684A", fontsize=8)
    ax[1].set_xticks(range(len(names)))
    ax[1].set_xticklabels(names, fontsize=7)
    ax[1].set_ylabel(r"designs with refereed $\Theta_{max} > 0$")
    ax[1].set_title("of 400 designs (200 graphs x 2 sigma)", fontsize=9)
    for i, v in enumerate(vals):
        ax[1].text(i, v + 0.8, str(v), ha="center", fontsize=8)

    fig.tight_layout()
    out = os.path.join(D, "k9b_hist.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print("wrote %s" % out)


if __name__ == "__main__":
    gallery()
    hist()
