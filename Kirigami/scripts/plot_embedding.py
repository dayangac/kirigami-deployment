#!/usr/bin/env python3
"""Plot the uncut mesh M and the kirigami structure M' at several opening angles.

Reads the JSON produced by the C++ tools (kiri_analyze / kiri_deploy / kiri_reference).
Run under Rosetta-free Python:  arch -arm64 /usr/local/bin/python3 plot_embedding.py ...

usage:
  plot_embedding.py case <case_dir> [-o out.png]     # M.json + deploy_*.json panels
  plot_embedding.py hist <sweep.csv> [-o out.png]    # dim_null vs (#interior - H)
"""
import json
import sys
import os
import csv
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

CW = "#5B8FF9"   # sigma = +1  (clockwise)
CCW = "#F6BD16"  # sigma = -1  (counter-clockwise)
HOLE = "#E8684A"


def _axis(ax, title):
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=9)


def draw_mesh(ax, path, title):
    d = json.load(open(path))
    V = d["vertices"]
    sig = d.get("orientation")
    patches, colors = [], []
    for i, f in enumerate(d["faces"]):
        patches.append(Polygon([V[k] for k in f], closed=True))
        colors.append(CW if (sig and sig[i] == 1) else CCW)
    ax.add_collection(PatchCollection(patches, facecolor=colors, edgecolor="0.25", linewidths=0.4))
    ax.autoscale_view()
    ax.relim()
    xs = [p[0] for p in V]
    ys = [p[1] for p in V]
    m = 0.05 * max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
    ax.set_xlim(min(xs) - m, max(xs) + m)
    ax.set_ylim(min(ys) - m, max(ys) + m)
    _axis(ax, title)


def draw_deployment(ax, path, title):
    d = json.load(open(path))
    V = d["vertices"]
    patches = [Polygon([V[k] for k in f], closed=True) for f in d["faces"]]
    ax.add_collection(PatchCollection(patches, facecolor="#BDD2FD", edgecolor="0.25",
                                      linewidths=0.4))
    holes = [Polygon([V[k] for k in h], closed=True) for h in d.get("holes", []) if len(h) >= 3]
    if holes:
        ax.add_collection(PatchCollection(holes, facecolor=HOLE, alpha=0.45, edgecolor="none"))
    xs = [p[0] for p in V]
    ys = [p[1] for p in V]
    m = 0.05 * max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
    ax.set_xlim(min(xs) - m, max(xs) + m)
    ax.set_ylim(min(ys) - m, max(ys) + m)
    _axis(ax, "%s  (%d holes)" % (title, len(d.get("holes", []))))


def plot_case(case_dir, out):
    panels = [("M.json", "M (input, sigma coloured)"),
              ("X0.json", "X0 (Eq. 6 projection)"),
              ("deploy_30deg.json", "M' at theta = 30 deg"),
              ("deploy_60deg.json", "M' at theta = 60 deg"),
              ("deploy_thetamax.json", "M' at theta_max")]
    have = [(f, t) for f, t in panels if os.path.exists(os.path.join(case_dir, f))]
    fig, axes = plt.subplots(1, len(have), figsize=(3.1 * len(have), 3.3))
    if len(have) == 1:
        axes = [axes]
    for ax, (f, t) in zip(axes, have):
        p = os.path.join(case_dir, f)
        if f.startswith("deploy"):
            draw_deployment(ax, p, t)
        else:
            draw_mesh(ax, p, t)
    fig.suptitle(os.path.basename(os.path.normpath(case_dir)), fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("wrote", out)


def plot_hist(csv_path, out):
    dim, claim = [], []
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            dim.append(int(row["dim_null"]))
            claim.append(int(row["n_interior_minus_H"]))
    diff = [a - b for a, b in zip(dim, claim)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].scatter(claim, dim, s=14, alpha=0.7, color=CW)
    lim = [0, max(max(claim, default=1), max(dim, default=1)) * 1.05]
    axes[0].plot(lim, lim, "k--", lw=0.8)
    axes[0].set_xlabel("#interior vertices - H")
    axes[0].set_ylabel("dim null(Eq. 4)")
    axes[0].set_title("Rank claim, %d sweep graphs" % len(dim))
    axes[1].hist(diff, bins=range(min(diff, default=0) - 1, max(diff, default=0) + 3),
                 color=CCW, edgecolor="0.3")
    axes[1].set_xlabel("dim_null - (#interior - H)")
    axes[1].set_ylabel("count")
    axes[1].set_title("Deviation (0 = claim holds)")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("wrote", out, "| violations:", sum(1 for d in diff if d != 0), "of", len(diff))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["case", "hist"])
    ap.add_argument("path")
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args()
    if a.mode == "case":
        plot_case(a.path, a.out or os.path.join(a.path, "figure.png"))
    else:
        plot_hist(a.path, a.out or "rank_claim.png")
