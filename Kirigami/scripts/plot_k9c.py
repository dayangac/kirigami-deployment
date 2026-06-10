#!/usr/bin/env python3
"""K9c figures.

  k9c_gallery.png -- the 20 best designs by certified eps_max, each shown closed
                     (theta = 0) and at theta = Theta_max / 2, from the JSON dumped by
                     kill_k9c. Sorted by eps_max, descending.
  k9c_hist.png    -- the certified eps_max distribution of the range-maximising embedding
                     against every baseline on the same 400 designs, plus the per-solver
                     counts (K9's proximity arm, K9b's proximity arm, K9c, best of three).

Run with: arch -arm64 /usr/local/bin/python3 code/scripts/plot_k9c.py [results/kill/k9c]
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

D = sys.argv[1] if len(sys.argv) > 1 else "results/kill/k9c"
G = os.path.join(D, "gallery")

CLOSED = "#BDD2FD"
OPEN = "#5B8FF9"


def load_index():
    rows = []
    for p in sorted(glob.glob(os.path.join(G, "index_*.json"))):
        rows.extend(json.load(open(p)))
    if not rows:
        # The shard index files are written only when a shard finishes. Rebuild the same
        # records from the dumped geometry plus k9c.csv, so a partially aggregated run
        # still produces the gallery.
        csvp = os.path.join(D, "k9c.csv")
        if os.path.exists(csvp):
            by_tag = {}
            for r in csv.DictReader(open(csvp)):
                by_tag["%s_%s_%s" % (r["kind"], r["id"], r["sigma"])] = r
            for q in sorted(glob.glob(os.path.join(G, "*_closed.json"))):
                tag = os.path.basename(q)[:-len("_closed.json")]
                r = by_tag.get(tag)
                if r is None:
                    continue
                rows.append({"tag": tag, "id": int(r["id"]), "kind": r["kind"],
                             "sigma": r["sigma"], "F": int(r["F"]),
                             "theta_max": float(r["theta_exact"]),
                             "eps_max": float(r["eps_max"]),
                             "src": r["best_src"]})
    # one entry per graph: the sigma with the larger certified eps_max (best-of-2)
    best = {}
    for r in rows:
        k = (r["kind"], r["id"])
        if k not in best or r["eps_max"] > best[k]["eps_max"]:
            best[k] = r
    return sorted(best.values(), key=lambda r: -r["eps_max"])


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
    entries = [e for e in load_index()
               if os.path.exists(os.path.join(G, e["tag"] + "_closed.json"))
               and os.path.exists(os.path.join(G, e["tag"] + "_half.json"))]
    if not entries:
        print("gallery: nothing dumped")
        return 0
    entries = entries[:n_want]
    ncol = 8  # 4 designs per row, each as a (closed, half) pair
    nrow = (len(entries) + 3) // 4
    fig, axes = plt.subplots(nrow, ncol, figsize=(1.55 * ncol, 1.75 * nrow),
                             squeeze=False)
    for a in axes.ravel():
        a.axis("off")
    for i, e in enumerate(entries):
        r, cpair = divmod(i, 4)
        draw(axes[r][2 * cpair], os.path.join(G, e["tag"] + "_closed.json"), CLOSED,
             "%s %d %s\nF=%d  closed" % (e["kind"][:4], e["id"],
                                         "mc" if e["sigma"].endswith("mc") else "def",
                                         e["F"]))
        draw(axes[r][2 * cpair + 1], os.path.join(G, e["tag"] + "_half.json"), OPEN,
             r"$\theta=\Theta_{max}/2$" + "\n%.3f rad, $\\varepsilon$=%.2f"
             % (0.5 * e["theta_max"], e["eps_max"]))
    fig.suptitle("K9c: the %d best random graphs by certified $\\varepsilon_{max}$ under "
                 "the RANGE-MAXIMISING constrained embedding\n(closed and half deployed; "
                 "every baseline on the same 400 designs is $\\Theta_{max}=0$)"
                 % len(entries), fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    out = os.path.join(D, "k9c_gallery.png")
    fig.savefig(out, dpi=190)
    plt.close(fig)
    print("wrote %s with %d designs" % (out, len(entries)))
    return len(entries)


def hist():
    path = os.path.join(D, "k9c.csv")
    if not os.path.exists(path):
        print("no k9c.csv")
        return
    rows = list(csv.DictReader(open(path)))
    eps = [float(r["eps_max"]) for r in rows if float(r["eps_max"]) > 0]
    fig, ax = plt.subplots(1, 2, figsize=(10.4, 3.8))

    ax[0].hist(eps, bins=24, color="#5B8FF9", edgecolor="0.25", linewidth=0.5)
    ax[0].axvline(0.1, color="#E8684A", ls="--", lw=1.2)
    ax[0].text(0.105, ax[0].get_ylim()[1] * 0.92, "0.1 rad target", color="#E8684A",
               fontsize=8)
    ax[0].set_xlabel(r"certified $\varepsilon_{max}$ (rad)")
    ax[0].set_ylabel("designs")
    ax[0].set_title("K9c best-of-3: %d certified designs of %d" % (len(eps), len(rows)),
                    fontsize=9)

    # Every baseline on this exact population is identically zero (K1a, K5, K6, K9a, B4,
    # K8a, and the authors' native `prevent`). The three solver arms are counted from the
    # same CSV; only the best-of-3 column is refereed by bisection.
    k9 = sum(1 for r in rows if float(r["k9_theta"]) > 1e-9)
    k9b = sum(1 for r in rows if float(r["k9b_theta"]) > 1e-9)
    k9c = sum(1 for r in rows if float(r["k9c_theta"]) > 1e-9)
    best = sum(1 for r in rows
               if float(r["theta_exact"]) > 1e-9 and float(r["theta_ref9"]) > 1e-9)
    names = ["Eq. (6)\n(K1a)", r"$\sigma_{mc}$" + "\n(K5)", r"$\sigma_{def}$" + "\n(K5)",
             "0+ repair\n(K6)", "convexity\nonly (K9a)",
             "prox. arm\n(K9)", "prox. arm\n(K9b-like)", "K9c\nrange-max",
             "best of\nthree"]
    vals = [0, 0, 0, 0, 0, k9, k9b, k9c, best]
    cols = ["0.7"] * 5 + ["#9BB7F0", "#9BB7F0", "#5B8FF9", "#3363C7"]
    ax[1].bar(range(len(vals)), vals, color=cols, edgecolor="0.25", linewidth=0.5)
    ax[1].axhline(40, color="#E8684A", ls="--", lw=1.2)
    ax[1].text(0.05, 41, "PASS bar 40 / 400", color="#E8684A", fontsize=8)
    ax[1].set_xticks(range(len(names)))
    ax[1].set_xticklabels(names, fontsize=6.2)
    ax[1].set_ylabel(r"designs with exact $\Theta_{max} > 0$")
    ax[1].set_title("of the %d designs measured (population: 200 graphs x 2 sigma = 400)"
                    % len(rows), fontsize=9)
    for i, v in enumerate(vals):
        ax[1].text(i, v + 0.8, str(v), ha="center", fontsize=8)

    fig.tight_layout()
    out = os.path.join(D, "k9c_hist.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print("wrote %s" % out)


if __name__ == "__main__":
    gallery()
    hist()
