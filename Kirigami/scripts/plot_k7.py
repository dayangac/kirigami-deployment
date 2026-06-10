#!/usr/bin/env python3
"""K7 figures: Poisson-ratio curves (predicted vs measured) and the achievable-K
dimension histogram. Run with `arch -arm64 /usr/local/bin/python3`."""
import csv
import os
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = sys.argv[1] if len(sys.argv) > 1 else "results/kill/k7"


def nu_figure():
    names = ["hexagons_2x2", "snub_square_2x2", "t3_4_3_12_2x2", "voronoi_torus_0_n20"]
    files = [(n, os.path.join(D, f"k7_nu_{n}.csv")) for n in names]
    files = [(n, p) for n, p in files if os.path.exists(p)]
    if not files:
        return
    fig, axes = plt.subplots(1, len(files), figsize=(4.2 * len(files), 3.6), sharey=True)
    if len(files) == 1:
        axes = [axes]
    for ax, (name, path) in zip(axes, files):
        rows = list(csv.DictReader(open(path)))
        series = defaultdict(lambda: ([], [], []))
        for r in rows:
            if r["design"] != "random_in_K":
                continue
            th, p, m = float(r["theta"]), float(r["nu_pred"]), float(r["nu_meas"])
            if abs(p) > 20:
                continue
            s = series[int(r["dir_deg"])]
            s[0].append(th)
            s[1].append(p)
            s[2].append(m)
        for k in sorted(series):
            th, p, m = series[k]
            ln, = ax.plot(th, p, "-", lw=1.6, label=f"predicted, d = {k}°")
            ax.plot(th[::4], m[::4], "o", ms=3.2, mfc="none", color=ln.get_color())
        # the conformal design of the same pattern
        cth, cnu = [], []
        for r in rows:
            if r["design"] == "conformal" and int(r["dir_deg"]) == 0:
                cth.append(float(r["theta"]))
                cnu.append(float(r["nu_pred"]))
        if cth:
            ax.plot(cth, cnu, "k--", lw=1.2, label="conformal design")
        ax.set_title(name, fontsize=9)
        ax.set_xlabel(r"$\theta$")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel(r"Poisson ratio $\nu(\theta)$")
    axes[0].legend(fontsize=6.5)
    fig.suptitle("K7 C3: Poisson ratio from the closed-form K (lines) vs forward "
                 "kinematics (circles)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(D, "k7_nu_curves.png"), dpi=160)


def dim_figure():
    path = os.path.join(D, "k7_main.csv")
    if not os.path.exists(path):
        return
    rows = [r for r in csv.DictReader(open(path)) if r.get("dimK")]
    c = Counter(int(r["dimK"]) for r in rows)
    ub = Counter(int(r["dimK_min4_2k"]) for r in rows)
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
    ks = [0, 1, 2, 3, 4]
    ax[0].bar([k - 0.18 for k in ks], [c.get(k, 0) for k in ks], width=0.36,
              label=r"measured $\dim\mathcal{K}$")
    ax[0].bar([k + 0.18 for k in ks], [ub.get(k, 0) for k in ks], width=0.36,
              label=r"bound $\min(4, 2\,\dim\,\mathrm{null})$")
    ax[0].set_xticks(ks)
    ax[0].set_xlabel(r"$\dim\mathcal{K}$")
    ax[0].set_ylabel("patterns")
    ax[0].legend(fontsize=8)
    ax[0].set_title(f"achievable-set dimension ({len(rows)} periodic patterns)", fontsize=9)
    ax[1].plot([2 * int(r["rankD"]) for r in rows], [int(r["dimK"]) for r in rows], "o",
               ms=5, alpha=0.6)
    ax[1].plot([0, 4], [0, 4], "k--", lw=1)
    ax[1].set_xlabel(r"$2\,\mathrm{rank}(D)$")
    ax[1].set_ylabel(r"$\dim\mathcal{K}$")
    ax[1].set_title(r"$\dim\mathcal{K} = 2\,\mathrm{rank}(D)$ on every pattern", fontsize=9)
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(D, "k7_dimK_hist.png"), dpi=160)


def c4_figure():
    """Where the hole-area harmonic and forward kinematics part company: only on the
    theta-interval where the deployed cell has inverted (det P_theta < 0), because the
    measured area uses |det P| while the harmonic is the signed det P_0 (det J - 1)."""
    path = os.path.join(D, "k7_c4_sweep.csv")
    if not os.path.exists(path):
        return
    rows = list(csv.DictReader(open(path)))
    names = []
    for r in rows:
        if r["name"] not in names:
            names.append(r["name"])
    fig, axes = plt.subplots(1, len(names), figsize=(3.6 * len(names), 3.2), sharex=True)
    if len(names) == 1:
        axes = [axes]
    for ax, name in zip(axes, names):
        rr = [r for r in rows if r["name"] == name]
        th = [float(r["theta"]) for r in rr]
        dev = [max(float(r["rel_dev"]), 1e-17) for r in rr]
        det = [float(r["detP"]) for r in rr]
        ax.semilogy(th, dev, "-", lw=1.4, color="C3")
        ax.set_title(name, fontsize=9)
        ax.set_xlabel(r"$\theta$")
        ax.grid(alpha=0.3)
        # shade where det P_theta < 0
        for i in range(len(th) - 1):
            if det[i] < 0:
                ax.axvspan(th[i], th[i + 1], color="C0", alpha=0.18, lw=0)
    axes[0].set_ylabel("relative $|A_{\\mathrm{fk}} - A_{\\mathrm{harmonic}}|$")
    fig.suptitle(r"K7 C4: the harmonic is exact except where $\det P_\theta < 0$ "
                 r"(shaded), an artefact of measuring $|\det P_\theta|$", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(D, "k7_c4_area_fit.png"), dpi=160)


def c3_figure():
    """C3 outcome per pattern: the target is always hit to machine precision; what
    fails is certifying a positive deployment range at the hitting design."""
    # prefer the post-F32 certificate run when it is present
    path = os.path.join(D, "k7_c3_all_v2.csv")
    if not os.path.exists(path):
        path = os.path.join(D, "k7_c3_all.csv")
    if not os.path.exists(path):
        return
    rows = list(csv.DictReader(open(path)))
    fams = []
    for r in rows:
        if r["family"] not in fams:
            fams.append(r["family"])
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.4))
    for tgt, a in zip(["random_in_K", "diag_1_-0.5"], ax):
        tot = [sum(1 for r in rows if r["family"] == f and r["target"] == tgt) for f in fams]
        cer = [sum(1 for r in rows if r["family"] == f and r["target"] == tgt
                   and r["cert"] == "1") for f in fams]
        fea = [sum(1 for r in rows if r["family"] == f and r["target"] == tgt
                   and r["zp_feasible"] == "1") for f in fams]
        x = range(len(fams))
        a.bar([i - 0.26 for i in x], tot, width=0.26, label="patterns", color="0.75")
        a.bar(x, fea, width=0.26, label=r"$0^+$ feasible", color="C0")
        a.bar([i + 0.26 for i in x], cer, width=0.26, label=r"certified $\Theta_{\max}>0$",
              color="C2")
        a.set_xticks(list(x))
        a.set_xticklabels(fams, rotation=35, ha="right", fontsize=7)
        a.set_title(f"target = {tgt}", fontsize=9)
        a.set_ylabel("patterns")
        a.grid(axis="y", alpha=0.3)
    ax[0].legend(fontsize=7)
    fig.suptitle("K7 C3: every target is hit exactly; certification is what fails",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(D, "k7_c3_certified.png"), dpi=160)


nu_figure()
dim_figure()
c4_figure()
c3_figure()
print("wrote", D)
