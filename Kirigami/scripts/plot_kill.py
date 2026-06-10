#!/usr/bin/env python3
"""Figures for results/kill/. Reads only CSV dumped by the C++ drivers (directive D3).
Run: arch -arm64 /usr/local/bin/python3 code/scripts/plot_kill.py"""
import csv, math, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
K = os.path.join(ROOT, "results", "kill")


def rows(p):
    with open(p) as f:
        return list(csv.DictReader(f))


def fnum(r, k, d=float("nan")):
    try:
        return float(r[k])
    except (KeyError, ValueError, TypeError):
        return d


def k2a():
    p = os.path.join(K, "k2a", "k2a.csv")
    if not os.path.exists(p):
        return
    R = rows(p)
    roots = np.array([fnum(r, "theta_roots") for r in R])
    t422 = np.array([fnum(r, "theta_t422") for r in R])
    bis = np.array([fnum(r, "bis_1e12") for r in R])
    graze = t422 > roots + 1e-7
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    ax[0].plot([0, math.pi], [0, math.pi], "k--", lw=0.8, label="exact agreement")
    ax[0].scatter(bis[~graze], t422[~graze], s=22, c="#2b6cb0", label="T4.2$''$ scan")
    ax[0].scatter(bis[graze], t422[graze], s=60, facecolors="none", edgecolors="#c53030",
                  lw=1.6, label="graze cases (%d)" % graze.sum())
    ax[0].scatter(bis[graze], roots[graze], s=40, marker="x", c="#c53030",
                  label="min over roots (wrong)")
    for i in np.where(graze)[0]:
        ax[0].annotate("", xy=(bis[i], t422[i]), xytext=(bis[i], roots[i]),
                       arrowprops=dict(arrowstyle="->", color="#c53030", lw=1.0))
    ax[0].set_xlabel(r"$\Theta_{max}$ by collision bisection (shrink $10^{-12}$)")
    ax[0].set_ylabel(r"closed-form $\Theta_{max}$")
    ax[0].set_title("K2a: closed form vs bisection, %d configurations" % len(R))
    ax[0].legend(fontsize=8, loc="upper left")
    ax[0].grid(alpha=0.3)

    e422 = np.abs(t422 - bis)
    eroot = np.abs(roots - bis)
    lo = 1e-16
    ax[1].hist(np.log10(np.maximum(e422, lo)), bins=30, alpha=0.75, color="#2b6cb0",
               label=r"T4.2$''$ scan (187/187 pass)")
    ax[1].hist(np.log10(np.maximum(eroot, lo)), bins=30, alpha=0.6, color="#c53030",
               label="min over roots (182/187)")
    ax[1].axvline(math.log10(1e-5), color="k", ls="--", lw=1.0, label=r"rule: $10^{-5}$ rad")
    ax[1].set_xlabel(r"$\log_{10}\,|\Theta_{max}^{\rm closed}-\Theta_{max}^{\rm bisect}|$")
    ax[1].set_ylabel("configurations")
    ax[1].set_title("K2a: error distribution")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(K, "k2a", "k2a_agreement.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


def k1c():
    p = os.path.join(K, "k1c", "k1c.csv")
    if not os.path.exists(p):
        return
    R = rows(p)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    labels, naive, samerange, corr, binding = [], [], [], [], []
    for src, lab in (("native_eq9", "authors' native\nprevent"),
                     ("ours_eq9", "our\noptimize_collision_sweep")):
        S = [r for r in R if r["source"] == src and r["valid_Y9"] == "1"]
        if not S:
            continue
        n = len(S)
        labels.append("%s\n(%d certified)" % (lab, n))
        naive.append(100.0 * sum(1 for r in S if int(r["n_naive"]) > 0) / n)
        samerange.append(100.0 * sum(1 for r in S if int(r["n_corr_tmax"]) > 0) / n)
        corr.append(100.0 * sum(1 for r in S if int(r["n_below_beta"]) > 0) / n)
        binding.append(100.0 * sum(1 for r in S if int(r["n_binding"]) > 0) / n)
    x = np.arange(len(labels))
    w = 0.2
    ax.bar(x - 1.5 * w, naive, w, color="#c53030",
           label=r"rule as written: $\theta^*<\Theta_{max}$, no clauses")
    ax.bar(x - 0.5 * w, samerange, w, color="#dd6b20",
           label=r"clauses + same range $\Theta_{max}$ (empty by definition)")
    ax.bar(x + 0.5 * w, corr, w, color="#2b6cb0",
           label=r"corrected: re-closure below $\min\beta$")
    ax.bar(x + 1.5 * w, binding, w, color="#2f855a",
           label=r"corrected: re-closure is the binding contact")
    ax.axhline(3.0, color="k", ls="--", lw=1.0, label="PASS threshold 3%")
    for off, vals in ((-1.5 * w, naive), (-0.5 * w, samerange), (0.5 * w, corr),
                      (1.5 * w, binding)):
        for xi, v in zip(x + off, vals):
            ax.text(xi, v + 1, "%.1f%%" % v, ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("designs with at least one Eq. (9) false negative (%)")
    ax.set_title("K1c: the rule as written counts a set that is empty once the\n"
                 "interval and crossing clauses are applied")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    out = os.path.join(K, "k1c", "k1c_rates.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


def k3a():
    p = os.path.join(K, "k3a", "k3a.csv")
    if not os.path.exists(p):
        return
    R = [r for r in rows(p) if r["kind"] == "delaunay"]
    if not R:
        return
    mf = np.array([fnum(r, "m_full") for r in R])
    mc = np.array([fnum(r, "m_core") for r in R])
    F = np.array([fnum(r, "F") for r in R])
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    hi = max(mf.max(), mc.max())
    ax[0].plot([0, hi], [0, hi], "k--", lw=0.8, label="no change")
    ax[0].scatter(mf, mc, s=18, c=F, cmap="viridis")
    ax[0].axhline(5, color="#c53030", lw=1.2, ls=":",
                  label=r"adversary's prediction $m_{core}\leq 5$")
    ax[0].set_xlabel(r"$m_{full}$ (whole hinge graph)")
    ax[0].set_ylabel(r"$m_{core}$ (2-core only)")
    ax[0].set_title("K3a: deleting dangling faces removes ~5%% of mobility (n=%d)" % len(R))
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)
    ax[1].scatter(F, mc, s=18, c="#2b6cb0", label=r"$m_{core}$")
    ax[1].axhline(5, color="#c53030", lw=1.2, ls=":", label=r"$m_{core}\leq 5$: 0 of %d" % len(R))
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("faces")
    ax[1].set_ylabel(r"$m_{core}$")
    ax[1].set_title("K3a: intrinsic mobility grows with the patch")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    out = os.path.join(K, "k3a", "k3a_mobility.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


def k5():
    p = os.path.join(K, "k5", "k5.csv")
    if not os.path.exists(p):
        return
    R = rows(p)
    if not R:
        return
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
    Dmc = np.array([fnum(r, "mc_D") for r in R])
    Ddef = np.array([fnum(r, "def_D") for r in R])
    m = (Dmc > 0) & (Ddef > 0)
    ax[0].loglog(Dmc[m], Ddef[m], "o", ms=4, color="#2b6cb0")
    lo, hi = min(Dmc[m].min(), Ddef[m].min()), max(Dmc[m].max(), Ddef[m].max())
    ax[0].plot([lo, hi], [lo, hi], "k--", lw=0.8, label="no change")
    ax[0].set_xlabel(r"$D(\sigma_{mc})$  (max-cut, Eq. 1)")
    ax[0].set_ylabel(r"$D(\sigma_{def})$  (defect search)")
    ax[0].set_title("K5: deployability defect, %d graphs" % len(R))
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3, which="both")

    pmc = np.array([fnum(r, "mc_proj_rel") for r in R])
    pdf = np.array([fnum(r, "def_proj_rel") for r in R])
    ax[1].hist(pmc, bins=30, alpha=0.7, color="#c53030", label=r"$\sigma_{mc}$")
    ax[1].hist(pdf, bins=30, alpha=0.7, color="#2b6cb0", label=r"$\sigma_{def}$")
    ax[1].set_xlabel(r"$\|X_0-X_{ini}\|_\infty$ / median edge length")
    ax[1].set_ylabel("graphs")
    ax[1].set_title("K5: how far Eq. (6) has to move the mesh")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    n = len(R)
    cats = ["overlap-free\nat $\\theta=0$", "certificate\nPOS/NOOVL/NOROOT",
            r"$\Theta_{max}>0$"]
    mcv = [sum(1 for r in R if r["mc_overlap_free0"] == "1"),
           sum(1 for r in R if r["mc_cert"] == "1"),
           sum(1 for r in R if fnum(r, "mc_theta") > 1e-9)]
    dfv = [sum(1 for r in R if r["def_overlap_free0"] == "1"),
           sum(1 for r in R if r["def_cert"] == "1"),
           sum(1 for r in R if fnum(r, "def_theta") > 1e-9)]
    x = np.arange(3)
    w = 0.36
    ax[2].bar(x - w / 2, [100.0 * v / n for v in mcv], w, color="#c53030",
              label=r"$\sigma_{mc}$")
    ax[2].bar(x + w / 2, [100.0 * v / n for v in dfv], w, color="#2b6cb0",
              label=r"$\sigma_{def}$")
    for xi, v in zip(x - w / 2, mcv):
        ax[2].text(xi, 100.0 * v / n + 1, str(v), ha="center", fontsize=8)
    for xi, v in zip(x + w / 2, dfv):
        ax[2].text(xi, 100.0 * v / n + 1, str(v), ha="center", fontsize=8)
    ax[2].set_xticks(x)
    ax[2].set_xticklabels(cats, fontsize=8)
    ax[2].set_ylabel("% of 200 graphs")
    ax[2].set_ylim(0, 70)
    ax[2].set_title("K5: it fixes the flat sheet, not the deployment")
    ax[2].legend(fontsize=8)
    ax[2].grid(alpha=0.3, axis="y")
    fig.tight_layout()
    out = os.path.join(K, "k5", "k5_orientation.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


def k2c_drift():
    p = os.path.join(K, "k2c", "k2c_drift.csv")
    if not os.path.exists(p):
        return
    # the growing-patch table lives in summary.txt; parse it
    sp = os.path.join(K, "k2c", "summary.txt")
    diam, drift = [], []
    for line in open(sp):
        t = line.split()
        if len(t) >= 3 and t[0] == "squares":
            diam.append(float(t[1]))
            drift.append(float(t[2]))
    if len(diam) < 3:
        return
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(diam, drift, "o-", ms=5, color="#c53030", label="growing square patch")
    b, a = np.polyfit(diam, drift, 1)
    xs = np.linspace(min(diam), max(diam), 10)
    ax.plot(xs, a + b * xs, "k--", lw=1.0, label="linear fit, slope %.4f" % b)
    ax.axhline(1.0, color="#2f855a", lw=1.2, ls=":",
               label=r"H-LOC would need a constant $\kappa$")
    ax.set_xlabel("patch diameter")
    ax.set_ylabel(r"$\max_f \|\gamma_f(\theta)-\bar{x}_f\| / r_f$")
    ax.set_title("K2c: centroid drift grows linearly - H-LOC is refuted")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(K, "k2c", "k2c_hloc.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


def k2c():
    p = os.path.join(K, "k2c", "k2c.csv")
    if not os.path.exists(p):
        return
    R = rows(p)
    key_pairs = next((k for k in R[0] if "pairs" in k and "per" in k), None)
    F = np.array([fnum(r, "F") for r in R])
    rho = np.array([fnum(r, "max_rho_over_r") for r in R]) if "max_rho_over_r" in R[0] else None
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    if rho is not None and np.isfinite(rho).any():
        ax[0].scatter(F, rho, s=14, c="#2b6cb0")
        ax[0].axhline(1.0, color="#c53030", lw=1.0, ls="--", label=r"$\rho_f/r_f = 1$ exactly")
        ax[0].set_xscale("log")
        ax[0].set_ylim(0.999, 1.001)
        ax[0].legend(fontsize=8)
    ax[0].set_xlabel("faces")
    ax[0].set_ylabel(r"$\max_f \rho_f / r_f$")
    ax[0].set_title("K2c: the swept disc is exactly the flat circumdisc")
    ax[0].grid(alpha=0.3, which="both")
    if key_pairs:
        pp = np.array([fnum(r, key_pairs) for r in R])
        m = np.isfinite(pp) & (pp > 0) & (F > 0)
        ax[1].scatter(F[m], pp[m], s=14, c="#2f855a")
        if m.sum() > 2:
            b, a = np.polyfit(np.log(F[m]), np.log(pp[m]), 1)
            xs = np.linspace(np.log(F[m].min()), np.log(F[m].max()), 20)
            ax[1].plot(np.exp(xs), np.exp(a + b * xs), "k--", lw=1.2,
                       label="slope %.4f (rule $\\leq 0.3$)" % b)
            ax[1].legend(fontsize=8)
        ax[1].set_xscale("log")
        ax[1].set_yscale("log")
        ax[1].set_xlabel("faces")
        ax[1].set_ylabel("surviving broad-phase pairs per face")
        ax[1].set_title("K2c: the active set is O(n)")
        ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    out = os.path.join(K, "k2c", "k2c_locality.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


def k2b():
    p = os.path.join(K, "k2b", "k2b.csv")
    if not os.path.exists(p):
        return
    R = rows(p)
    if not R:
        return
    base = np.array([fnum(r, "theta_native_best") for r in R])
    ours = np.array([fnum(r, "theta_ours") for r in R])
    eq9 = np.array([fnum(r, "theta_ours_eq9") for r in R])
    names = [r["name"] for r in R]
    order = np.argsort(base)
    x = np.arange(len(R))
    fig, ax = plt.subplots(figsize=(max(7.5, 0.34 * len(R)), 4.6))
    ax.plot(x, base[order], "o-", ms=4, color="#c53030", label="authors' native prevent (best of ladder)")
    ax.plot(x, ours[order], "s-", ms=4, color="#2b6cb0", label="ours (softmin of closed-form roots)")
    ax.plot(x, eq9[order], "^", ms=4, color="#a0aec0", label="our Eq. (9) reimplementation")
    ax.set_xticks(x)
    ax.set_xticklabels([names[i] for i in order], rotation=90, fontsize=6)
    ax.set_ylabel(r"$\Theta_{max}$ (our bisection, shrink $10^{-12}$)")
    ax.set_title("K2b: range against the authors' native baseline (n=%d)" % len(R))
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(K, "k2b", "k2b_margin.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    for fn in (k2a, k1c, k3a, k2c, k2c_drift, k2b, k5):
        try:
            fn()
        except Exception as e:
            print("skip %s: %s" % (fn.__name__, e), file=sys.stderr)
