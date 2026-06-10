#!/usr/bin/env python3
"""K6 figure: inward split edges at t = 0, and the certified Theta_max distributions.

Run as:  arch -arm64 /usr/local/bin/python3 code/scripts/plot_k6.py results/kill/k6/k6.csv \
             -o results/kill/k6/k6_zero_plus.png
"""
import argparse
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FAMS = ["sigma_mc", "sigma_def"]
COL = {"sigma_mc": "#3b6ea5", "sigma_def": "#c2622d"}


def load(path):
    rows = {f: [] for f in FAMS}
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["sigma"] in rows:
                rows[r["sigma"]].append(r)
    return rows


def col(rs, key, cast=float):
    return np.array([cast(r[key]) for r in rs]) if rs else np.array([])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("-o", "--out", default="k6.png")
    a = ap.parse_args()
    rows = load(a.csv)
    if not any(rows.values()):
        sys.exit("no rows")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

    # (a) inward split edges at t = 0, as a fraction of |E_split|
    for f in FAMS:
        v = col(rows[f], "inward_frac0")
        if v.size:
            ax[0].hist(v, bins=np.linspace(0, 1, 26), alpha=0.6, label=f, color=COL[f])
    ax[0].set_xlabel(r"inward split edges $q_e\leq 0$ / $|E_{split}|$ at $t=0$")
    ax[0].set_ylabel("graphs")
    ax[0].set_title("(a) why K5 saw $\\Theta_{max}=0$")
    ax[0].legend()

    # (b) absolute counts vs |E_split|
    for f in FAMS:
        x, y = col(rows[f], "n_split"), col(rows[f], "inward0")
        if x.size:
            ax[1].scatter(x, y, s=12, alpha=0.6, label=f, color=COL[f])
    lim = max([col(rows[f], "n_split").max() if rows[f] else 0 for f in FAMS] + [1])
    ax[1].plot([0, lim], [0, lim], "k--", lw=0.8, label="all inward")
    ax[1].set_xlabel(r"$|E_{split}|$")
    ax[1].set_ylabel("inward split edges at $t=0$")
    ax[1].set_title("(b) count vs split-cut count")
    ax[1].legend()

    # (c) certified Theta_max after the repair, plus the refereeing bisection
    labels, data, colors = [], [], []
    for f in FAMS:
        rs = [r for r in rows[f] if int(r["feasible"])]
        for key, tag in (("eps_max", "certified"), ("theta_bisect", "bisection")):
            v = np.array([float(r[key]) for r in rs]) if rs else np.array([])
            labels.append(f"{f}\n{tag}\n(n={v.size})")
            data.append(v if v.size else np.array([0.0]))
            colors.append(COL[f])
    pos = np.arange(len(data))
    for p, d, c in zip(pos, data, colors):
        jitter = (np.random.RandomState(0).rand(d.size) - 0.5) * 0.3
        ax[2].scatter(p + jitter, d, s=14, alpha=0.7, color=c)
    ax[2].set_xticks(pos)
    ax[2].set_xticklabels(labels, fontsize=7)
    ax[2].set_ylabel(r"$\Theta_{max}$ (rad)")
    ax[2].set_title("(c) after the 0+ repair, at the $0^+$-feasible points")
    ax[2].set_ylim(-0.05, max(0.35, float(np.max([d.max() for d in data])) * 1.1))

    fig.suptitle("K6 -- 0+ repair in the Tutte auxetic null space, 200 random graphs")
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
