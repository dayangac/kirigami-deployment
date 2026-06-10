#!/usr/bin/env python3
"""fig_regime.png -- deployment yield against |F| from 20 to 800, three arms.

Joins WP7a's paper-regime population (results/regime/regime.csv, |F| in [20,100]) with the
population every earlier kill experiment used (|F| in [101,793]), so the two halves of the
face-count axis are read off one figure:

  Eq. (6) baseline   regime.csv base_theta                | results/kill/k1a/k1a.csv
                                                            theta_max_X0 (the Eq. (6)
                                                            projection, one row per graph)
  native pipeline    regime.csv natour_theta / natcol_theta | results/kill/native200/
                                                            native200.csv our_theta_exact,
                                                            completed cells only
  k9c range-max      regime.csv k9c_theta                  | results/kill/k9c/k9c.csv
                                                            k9c_theta

A design counts as deployable when the exact T4.2'' Theta_max is > 1e-9; every column above
is that same exact scan. Bars are Wilson 95 % intervals. Run under
`arch -arm64 /usr/local/bin/python3` (D3).
"""
import csv, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "results", "regime", "fig_regime.png")
EDGES = [20, 35, 50, 71, 100, 150, 220, 330, 480, 800]


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def wilson(k, n):
    if n == 0:
        return (math.nan, math.nan, math.nan)
    z, p = 1.959963984540054, k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    r = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return p, max(0.0, min(p, (c - r) / d)), min(1.0, max(p, (c + r) / d))


def collect():
    """-> {arm: [(F, deployable0or1), ...]}"""
    arms = {"Eq. (6) baseline": [], "authors' native pipeline": [], "k9c range-max": []}
    reg = rows(os.path.join(ROOT, "results", "regime", "regime.csv"))
    for r in reg:
        if r["set"] != "pop":
            continue
        F = int(r["F"])
        if r["base_status"] != "not_run":
            arms["Eq. (6) baseline"].append((F, float(r["base_theta"]) > 1e-9))
        if r["k9c_status"] != "not_run":
            arms["k9c range-max"].append((F, float(r["k9c_theta"]) > 1e-9))
        for pre in ("natour", "natcol"):
            if r[pre + "_status"] == "completed":
                arms["authors' native pipeline"].append(
                    (F, float(r[pre + "_theta"]) > 1e-9))
    k1a = os.path.join(ROOT, "results", "kill", "k1a", "k1a.csv")
    if os.path.exists(k1a):
        for r in rows(k1a):
            arms["Eq. (6) baseline"].append(
                (int(r["F"]), float(r["theta_max_X0"]) > 1e-9))
    k9c = os.path.join(ROOT, "results", "kill", "k9c", "k9c.csv")
    if os.path.exists(k9c):
        for r in rows(k9c):
            arms["k9c range-max"].append((int(r["F"]), float(r["k9c_theta"]) > 1e-9))
    nat = os.path.join(ROOT, "results", "kill", "native200", "native200.csv")
    if os.path.exists(nat):
        for r in rows(nat):
            if r["status"] == "completed" and r["embedding_ok"] == "1":
                arms["authors' native pipeline"].append(
                    (int(r["F"]), float(r["our_theta_exact"]) > 1e-9))
    return arms


def main():
    arms = collect()
    colors = {"Eq. (6) baseline": "#b03030",
              "authors' native pipeline": "#7a5cc0",
              "k9c range-max": "#1f6f4a"}
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    for name, data in arms.items():
        xs, ys, lo, hi, ns = [], [], [], [], []
        for a, b in zip(EDGES[:-1], EDGES[1:]):
            sel = [d for F, d in data if a <= F < b]
            if len(sel) < 3:
                continue
            p, l, h = wilson(sum(1 for s in sel if s), len(sel))
            xs.append(math.sqrt(a * b))
            ys.append(100 * p)
            lo.append(100 * (p - l))
            hi.append(100 * (h - p))
            ns.append(len(sel))
        if not xs:
            continue
        ax.errorbar(xs, ys, yerr=[lo, hi], marker="o", ms=5, lw=1.8, capsize=3,
                    color=colors[name], label="%s (n=%d)" % (name, sum(ns)))
    ax.axvline(100, color="0.6", ls="--", lw=1)
    ax.text(103, 52, "WP7a population  |  earlier kill population", fontsize=8,
            color="0.35", rotation=90, va="center")
    ax.set_xscale("log")
    ax.set_xlabel(r"faces $|F|$  (log scale)")
    ax.set_ylabel(r"designs with exact $\Theta_{\max} > 0$   [%]")
    ax.set_ylim(-4, 104)
    ax.set_xlim(18, 900)
    ax.set_xticks([20, 50, 100, 200, 400, 800])
    ax.set_xticklabels(["20", "50", "100", "200", "400", "800"])
    ax.grid(alpha=0.25)
    ax.legend(loc="center left", fontsize=9, framealpha=0.95)
    ax.set_title("Deployment yield vs. face count, three arms\n"
                 "bars: Wilson 95 % intervals; left of the dashed line is the 2026 "
                 "paper's own scale", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT, dpi=170)
    print("wrote", OUT)
    for name, data in arms.items():
        print(name, "n =", len(data))


if __name__ == "__main__":
    main()
