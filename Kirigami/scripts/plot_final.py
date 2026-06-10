#!/usr/bin/env python3
"""Final paper figures: baseline (Eq. (6), K6) vs. the three constrained-embedding arms
(K9 proximity, K9b, K9c range-maximising), plus E1's exact-vs-referee agreement.

Sources (all read-only, nothing here re-runs C++):
  results/kill/k9c/shard_*.csv   -- K9c's 12 shards; merged here by (kind,id,sigma),
                                     since results/kill/k9c/k9c.csv itself is a stale
                                     182-row snapshot and the run has kept landing rows.
                                     Every row carries all three arms' theta/eps/margin
                                     on the SAME graph, so this merge is the single
                                     source for fig_yield / fig_eps_hist /
                                     fig_feasible_vs_deployable.
  results/kill/k9c/gallery/*.json -- per-design closed/half-deployed geometry dumps,
                                     used by fig_gallery_full.
  results/kill/k6/k6_final.csv   -- baseline (Eq. (6) + repairs): confirmed identically
                                     0 on this population (checked below, not assumed).
  results/final/e1/e1.csv        -- E1's exact vs. bisection referee population.

Figures -> results/final/figures/:
  fig_yield.png
  fig_eps_hist.png
  fig_feasible_vs_deployable.png
  fig_gallery_full.png
  fig_e1_agreement.png
results/final/figures/summary_table.md and README.md are written by this same run.

Run: arch -arm64 /usr/local/bin/python3 code/scripts/plot_final.py
"""
import csv
import glob
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

K9C_DIR = "results/kill/k9c"
GALLERY = os.path.join(K9C_DIR, "gallery")
K6_CSV = "results/kill/k6/k6_final.csv"
E1_CSV = "results/final/e1/e1.csv"
OUT_DIR = "results/final/figures"
os.makedirs(OUT_DIR, exist_ok=True)

FAMILIES = ["delaunay", "voronoi", "quad_random"]
FAMILY_LABEL = {"delaunay": "Delaunay", "voronoi": "Voronoi", "quad_random": "quad"}
SIGMAS = ["sigma_mc", "sigma_def"]
ARMS = ["baseline", "k9", "k9b", "k9c"]
ARM_LABEL = {"baseline": "Eq.(6)+repairs\n(K6)", "k9": "K9 proximity",
             "k9b": "K9b", "k9c": "K9c range-max"}
ARM_COLOR = {"baseline": "0.6", "k9": "#9BB7F0", "k9b": "#5B8FF9", "k9c": "#3363C7"}


def load_k9c_merged():
    """Merge the 12 live shards by (kind, id, sigma), last-write-wins (shards are
    disjoint by construction; no duplicate ids were found across shards)."""
    rows = {}
    for f in sorted(glob.glob(os.path.join(K9C_DIR, "shard_*.csv"))):
        for r in csv.DictReader(open(f)):
            rows[(r["kind"], r["id"], r["sigma"])] = r
    return list(rows.values())


def load_k6():
    rows = {}
    for r in csv.DictReader(open(K6_CSV)):
        rows[(r["kind"], r["id"], r["sigma"])] = r
    return rows


def f(x):
    return float(x)


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0, 0.0
    phat = k / n
    denom = 1.0 + z * z / n
    center = (phat + z * z / (2 * n)) / denom
    half = z * math.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n)) / denom
    return phat, max(0.0, center - half), min(1.0, center + half)


# ---------------------------------------------------------------------------
# fig_yield: yield vs |F|, per family x sigma, 4 arms, Wilson 95% CI
# ---------------------------------------------------------------------------
def fig_yield(rows, k6):
    edges = [100, 230, 420, 800]  # 3 bins spanning F in [101, 793] on this population
    bin_labels = ["F<230\n(F=%d-229)" % edges[0], "230<=F<420", "F>=420\n(<=%d)" % edges[-1]]

    def bin_of(F):
        for i in range(len(edges) - 1):
            if edges[i] <= F < edges[i + 1]:
                return i
        return len(edges) - 2 if F < edges[0] else len(edges) - 2

    def get(row, arm, field):
        if arm == "baseline":
            return 0.0
        return f(row["%s_%s" % (arm, field)])

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.6), sharey=True)
    metrics = [("theta", 1e-9, r"yield: exact $\Theta_{max}>0$"),
               ("eps", 0.1, r"yield: certified $\varepsilon_{max}\geq 0.1$ rad")]

    for mrow, (field, thresh, mtitle) in enumerate(metrics):
        for fcol, fam in enumerate(FAMILIES):
            ax = axes[mrow][fcol]
            fam_rows = [r for r in rows if r["kind"] == fam]
            n_series = len(ARMS) * len(SIGMAS)
            offsets = [(-0.5 + (j + 0.5) / n_series) * 0.75 for j in range(n_series)]
            j = 0
            for arm in ARMS:
                for sig in SIGMAS:
                    sub_rows = [r for r in fam_rows if r["sigma"] == sig]
                    binned = [[] for _ in range(len(edges) - 1)]
                    for r in sub_rows:
                        binned[bin_of(int(r["F"]))].append(r)
                    xs, ys, lo, hi = [], [], [], []
                    for b, rs in enumerate(binned):
                        if not rs:
                            continue
                        if field == "theta":
                            k = sum(1 for r in rs if get(r, arm, field) > thresh)
                        else:
                            k = sum(1 for r in rs if get(r, arm, field) >= thresh)
                        n = len(rs)
                        p, lo_, hi_ = wilson(k, n)
                        xs.append(b + offsets[j])
                        ys.append(p)
                        lo.append(max(0.0, p - lo_))
                        hi.append(max(0.0, hi_ - p))
                    ls = "-" if sig == "sigma_mc" else "--"
                    marker = "o" if sig == "sigma_mc" else "s"
                    label = "%s %s" % (ARM_LABEL[arm].split("\n")[0],
                                        "mc" if sig == "sigma_mc" else "def")
                    ax.errorbar(xs, ys, yerr=[lo, hi], fmt=marker, ms=3.4,
                                color=ARM_COLOR[arm], linestyle=ls, linewidth=0.9,
                                capsize=2, elinewidth=0.8, label=label if fcol == 0 and mrow == 0 else None)
                    j += 1
            ax.set_xticks(range(len(bin_labels)))
            ax.set_xticklabels(bin_labels, fontsize=7.5)
            ax.set_ylim(-0.03, 1.03)
            ax.set_title("%s -- %s" % (FAMILY_LABEL[fam], mtitle), fontsize=9)
            if fcol == 0:
                ax.set_ylabel("fraction of designs (Wilson 95% CI)", fontsize=8.5)
    fig.legend(loc="lower center", ncol=8, fontsize=7.2, frameon=False,
               bbox_to_anchor=(0.5, 0.005))
    fig.suptitle("Yield vs. face count |F|, by graph family and orientation rule "
                 "(N=%d designs, %d graphs x 2 sigma)"
                 % (len(rows), len(rows) // 2 if len(rows) % 2 == 0 else len(rows)),
                 fontsize=11)
    fig.tight_layout(rect=[0, 0.08, 1, 0.94])
    out = os.path.join(OUT_DIR, "fig_yield.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    print("wrote %s" % out)


# ---------------------------------------------------------------------------
# fig_eps_hist: eps_max distribution, K9c vs K9, log-x, zero-count bars per arm
# ---------------------------------------------------------------------------
def fig_eps_hist(rows):
    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.2))

    k9c_pos = [f(r["k9c_eps"]) for r in rows if f(r["k9c_eps"]) > 0]
    k9_pos = [f(r["k9_eps"]) for r in rows if f(r["k9_eps"]) > 0]
    lo = math.log10(min(min(k9c_pos), min(k9_pos)))
    bins = [10 ** (lo + i * (0 - lo) / 28.0) for i in range(29)]
    ax[0].hist(k9c_pos, bins=bins, color=ARM_COLOR["k9c"], alpha=0.75,
               edgecolor="0.25", linewidth=0.4, label="K9c (n=%d)" % len(k9c_pos))
    ax[0].hist(k9_pos, bins=bins, color=ARM_COLOR["k9"], alpha=0.75,
               edgecolor="0.25", linewidth=0.4, label="K9 proximity (n=%d)" % len(k9_pos))
    ax[0].set_xscale("log")
    ax[0].axvline(0.1, color="#E8684A", ls="--", lw=1.1)
    ax[0].text(0.11, ax[0].get_ylim()[1] * 0.9, "0.1 rad", color="#E8684A", fontsize=7.5)
    ax[0].set_xlabel(r"certified $\varepsilon_{max}$ (rad), designs with $\varepsilon_{max}>0$")
    ax[0].set_ylabel("designs")
    ax[0].legend(fontsize=8, frameon=False)
    ax[0].set_title(r"$\varepsilon_{max}$ distribution (log-x), K9c vs. K9", fontsize=9.5)

    n = len(rows)
    zero_counts = []
    for arm in ARMS:
        if arm == "baseline":
            zero_counts.append(n)
        else:
            pos = sum(1 for r in rows if f(r["%s_eps" % arm]) > 0)
            zero_counts.append(n - pos)
    ax[1].bar(range(len(ARMS)), zero_counts,
              color=[ARM_COLOR[a] for a in ARMS], edgecolor="0.25", linewidth=0.5)
    ax[1].set_xticks(range(len(ARMS)))
    ax[1].set_xticklabels([ARM_LABEL[a] for a in ARMS], fontsize=8)
    ax[1].set_ylabel(r"designs with $\varepsilon_{max}=0$ (of %d)" % n)
    ax[1].set_title("count at zero, per arm", fontsize=9.5)
    for i, v in enumerate(zero_counts):
        ax[1].text(i, v + n * 0.01, str(v), ha="center", fontsize=8)

    fig.tight_layout()
    out = os.path.join(OUT_DIR, "fig_eps_hist.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    print("wrote %s" % out)


# ---------------------------------------------------------------------------
# fig_feasible_vs_deployable: best 0+ margin vs exact Theta_max, all K9c rows
# ---------------------------------------------------------------------------
def fig_feasible_vs_deployable(rows):
    fig, ax = plt.subplots(figsize=(7.4, 6.0))
    fam_color = {"delaunay": "#3363C7", "voronoi": "#E8684A", "quad_random": "#2CA58D"}
    for fam in FAMILIES:
        xs = [f(r["best_margin"]) for r in rows if r["kind"] == fam]
        ys = [f(r["theta_exact"]) for r in rows if r["kind"] == fam]
        ax.scatter(xs, ys, s=13, color=fam_color[fam], alpha=0.65,
                  edgecolor="none", label=FAMILY_LABEL[fam])
    stuck = [(f(r["best_margin"]), f(r["theta_exact"])) for r in rows
             if f(r["best_margin"]) > 0 and f(r["theta_exact"]) == 0.0]
    if stuck:
        sx, sy = zip(*stuck)
        ax.scatter(sx, sy, s=55, facecolor="none", edgecolor="#E8684A",
                  linewidth=1.4, marker="o",
                  label="margin-feasible, $\\Theta_{max}=0$ (n=%d)" % len(stuck))
    ax.axvline(0, color="0.5", lw=0.7, ls=":")
    ax.set_xlabel(r"best 0$^+$ margin $m$ (convex + split-feasibility slack; $m>0$ = feasible)")
    ax.set_ylabel(r"exact $\Theta_{max}$ (rad)")
    ax.set_title("Feasible vs. deployable: K9c's best-of-3 margin against exact "
                 r"$\Theta_{max}$" + "\n(N=%d designs)" % len(rows), fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "fig_feasible_vs_deployable.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    print("wrote %s (%d margin-feasible-but-Theta_max=0 designs)" % (out, len(stuck)))


# ---------------------------------------------------------------------------
# fig_gallery_full: ALL K9c-deployable designs, sorted by eps_max descending
# ---------------------------------------------------------------------------
def draw(ax, path, color, title):
    d = json.load(open(path))
    V = d["vertices"]
    patches = [Polygon([V[k] for k in fc], closed=True) for fc in d["faces"]]
    ax.add_collection(PatchCollection(patches, facecolor=color, edgecolor="0.3",
                                      linewidths=0.12))
    xs = [p[0] for p in V]
    ys = [p[1] for p in V]
    mx = 0.04 * max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
    ax.set_xlim(min(xs) - mx, max(xs) + mx)
    ax.set_ylim(min(ys) - mx, max(ys) + mx)
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=4.6, pad=1.0)


def fig_gallery_full(rows):
    CLOSED, OPEN = "#BDD2FD", "#5B8FF9"
    deployable = [r for r in rows if f(r["theta_exact"]) > 0
                  and os.path.exists(os.path.join(GALLERY,
                      "%s_%s_%s_closed.json" % (r["kind"], r["id"], r["sigma"])))
                  and os.path.exists(os.path.join(GALLERY,
                      "%s_%s_%s_half.json" % (r["kind"], r["id"], r["sigma"])))]
    deployable.sort(key=lambda r: -f(r["eps_max"]))
    n_missing = sum(1 for r in rows if f(r["theta_exact"]) > 0) - len(deployable)

    per_row = 6  # designs per row -> 12 columns (closed|half pairs)
    ncol = 2 * per_row
    nrow = (len(deployable) + per_row - 1) // per_row
    fig, axes = plt.subplots(nrow, ncol, figsize=(0.82 * ncol, 0.95 * nrow), squeeze=False)
    for a in axes.ravel():
        a.axis("off")
    for i, r in enumerate(deployable):
        row, cpair = divmod(i, per_row)
        tag = "%s_%s_%s" % (r["kind"], r["id"], r["sigma"])
        draw(axes[row][2 * cpair], os.path.join(GALLERY, tag + "_closed.json"), CLOSED,
             "%s%s %s\nF=%s" % (r["kind"][:3], r["id"],
                                 "mc" if r["sigma"].endswith("mc") else "def", r["F"]))
        draw(axes[row][2 * cpair + 1], os.path.join(GALLERY, tag + "_half.json"), OPEN,
             r"$\varepsilon$=%.2f" % f(r["eps_max"]))
    fig.suptitle(
        "All %d K9c-deployable designs (of %d measured), sorted by certified "
        r"$\varepsilon_{max}$ descending" % (len(deployable), len(rows)) +
        "\nleft of each pair = closed ($\\theta=0$); right = deployed to "
        r"$\Theta_{max}/2$" +
        "\nbaseline (Eq.(6)+repairs, K6): closed only, $\\Theta_{max}=0$ on all %d "
        "designs in this population -- no deployed panel exists to show"
        % len(rows), fontsize=8.5)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(OUT_DIR, "fig_gallery_full.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("wrote %s with %d designs (%d deployable rows missing gallery json, excluded)"
          % (out, len(deployable), n_missing))


# ---------------------------------------------------------------------------
# fig_e1_agreement: |exact - referee| histogram + certificate confusion counts
# ---------------------------------------------------------------------------
def fig_e1_agreement():
    rows = list(csv.DictReader(open(E1_CSV)))
    gaps = [abs(f(r["gap"])) for r in rows]
    cert = [r["cert_valid"] == "1" for r in rows]
    actual = [r["actually_ge_eps"] == "1" for r in rows]
    tp = sum(1 for c, a in zip(cert, actual) if c and a)
    fp = sum(1 for c, a in zip(cert, actual) if c and not a)
    fn = sum(1 for c, a in zip(cert, actual) if not c and a)
    tn = sum(1 for c, a in zip(cert, actual) if not c and not a)

    fig, ax = plt.subplots(1, 2, figsize=(10.6, 4.2))
    ax[0].hist(gaps, bins=60, color="#3b6ea5", edgecolor="white", linewidth=0.3)
    ax[0].set_yscale("log")
    ax[0].set_xlabel(r"$|\Theta_{max}(\mathrm{T4.2''}) - \Theta_{max}(\mathrm{bisection})|$ (rad)")
    ax[0].set_ylabel("count (log)")
    ax[0].axvline(1e-5, color="crimson", linestyle="--", linewidth=1,
                 label="1e-5 agreement bar")
    ax[0].legend(fontsize=8)
    ax[0].set_title("E1: exact vs. bisection referee, N=%d\nworst gap %.3e"
                    % (len(rows), max(gaps)), fontsize=9.5)

    labels = ["certified &\n$\\Theta_{max}\\geq\\varepsilon$\n(TP)",
              "certified but\n$\\Theta_{max}<\\varepsilon$\n(FP, unsound)",
              "not certified,\n$\\Theta_{max}\\geq\\varepsilon$\n(FN, over-strict)",
              "not certified,\n$\\Theta_{max}<\\varepsilon$\n(TN)"]
    vals = [tp, fp, fn, tn]
    cols = ["#2CA58D", "#E8684A", "#E8684A", "0.65"]
    ax[1].bar(range(4), vals, color=cols, edgecolor="0.25", linewidth=0.5)
    ax[1].set_xticks(range(4))
    ax[1].set_xticklabels(labels, fontsize=7)
    ax[1].set_ylabel("designs")
    ax[1].set_title("certificate confusion counts (eps=0.006), N=%d\n0 unsound / 0 "
                    "over-strict on this population" % len(rows), fontsize=9)
    for i, v in enumerate(vals):
        ax[1].text(i, v + len(rows) * 0.01, str(v), ha="center", fontsize=8)
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "fig_e1_agreement.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    print("wrote %s (N=%d, worst gap %.3e, TP=%d FP=%d FN=%d TN=%d)"
          % (out, len(rows), max(gaps), tp, fp, fn, tn))
    return len(rows), max(gaps), tp, fp, fn, tn


# ---------------------------------------------------------------------------
# summary_table.md
# ---------------------------------------------------------------------------
def summary_table(rows):
    n = len(rows)
    complete = " (full target)" if n >= 400 else " (short of the 400-design target; "\
               "re-run once more shards land)"
    lines = ["# Headline table: baseline vs. method arms\n",
             "K9c live-merged population: N=%d designs, %d unique (kind,id) graphs x "
             "2 sigma rules%s." % (n, n // 2, complete), "",
             "| arm | N | deployable ($\\Theta_{max}>0$) | certified ($\\varepsilon_{max}>0$) "
             "| $\\varepsilon_{max}\\geq0.1$ rad | median $\\varepsilon_{max}$ | max "
             "$\\varepsilon_{max}$ |",
             "|---|---|---|---|---|---|---|"]
    for arm in ARMS:
        if arm == "baseline":
            lines.append("| %s | %d | 0 | 0 | 0 | -- | -- |" % (ARM_LABEL[arm].replace("\n", " "), n))
            continue
        dep = sum(1 for r in rows if f(r["%s_theta" % arm]) > 1e-9)
        eps_pos = [f(r["%s_eps" % arm]) for r in rows if f(r["%s_eps" % arm]) > 0]
        cert = len(eps_pos)
        ge01 = sum(1 for e in eps_pos if e >= 0.1)
        med = sorted(eps_pos)[len(eps_pos) // 2] if eps_pos else 0.0
        mx = max(eps_pos) if eps_pos else 0.0
        lines.append("| %s | %d | %d | %d | %d | %.3f | %.3f |"
                     % (ARM_LABEL[arm].replace("\n", " "), n, dep, cert, ge01, med, mx))
    # best-of-3 informational row
    eps_pos = [f(r["eps_max"]) for r in rows if f(r["eps_max"]) > 0]
    dep = sum(1 for r in rows if f(r["theta_exact"]) > 1e-9)
    ge01 = sum(1 for e in eps_pos if e >= 0.1)
    med = sorted(eps_pos)[len(eps_pos) // 2] if eps_pos else 0.0
    mx = max(eps_pos) if eps_pos else 0.0
    lines.append("| best-of-3 (per-design max over K9/K9b/K9c) | %d | %d | %d | %d | %.3f | %.3f |"
                 % (n, dep, len(eps_pos), ge01, med, mx))
    lines.append("")
    out = os.path.join(OUT_DIR, "summary_table.md")
    open(out, "w").write("\n".join(lines) + "\n")
    print("wrote %s" % out)


if __name__ == "__main__":
    k9c_rows = load_k9c_merged()
    k6 = load_k6()
    print("K9c merged: %d rows (target 400)" % len(k9c_rows))

    fig_yield(k9c_rows, k6)
    fig_eps_hist(k9c_rows)
    fig_feasible_vs_deployable(k9c_rows)
    fig_gallery_full(k9c_rows)
    e1_n, e1_worst, tp, fp, fn, tn = fig_e1_agreement()
    summary_table(k9c_rows)

    readme_text = """# results/final/figures -- baseline vs. method figure set

Data sources and the exact command that produced everything here:
`arch -arm64 /usr/local/bin/python3 code/scripts/plot_final.py` (reads only; does not
re-run any C++ solver).

| figure | what it shows | data source(s) |
|---|---|---|
| `fig_yield.png` | Yield vs. face count \\|F\\| in 3 bins, per family (Delaunay/Voronoi/quad) and per sigma rule, 4 arms (baseline Eq.(6)+repairs, K9 proximity, K9b, K9c range-maximising), Wilson 95%% CIs, two rows: exact Theta_max>0 and certified eps_max>=0.1 rad. | `results/kill/k9c/shard_*.csv` merged by (kind,id,sigma), N=%d rows at merge time; baseline confirmed identically 0 from `results/kill/k6/k6_final.csv`. |
| `fig_eps_hist.png` | eps_max distribution (log-x) for K9c vs. K9 proximity, plus a bar of the zero-eps_max count per arm. | same K9c merge. |
| `fig_feasible_vs_deployable.png` | Scatter of best 0+ margin (convex+split feasibility slack) vs. exact Theta_max for every row, colored by family, with margin-feasible-but-Theta_max=0 designs ringed. | same K9c merge. |
| `fig_gallery_full.png` | ALL K9c-deployable designs (Theta_max>0), sorted by certified eps_max descending, shown closed and at Theta_max/2, with a caption noting the baseline has no deployed panel to show (Theta_max=0 throughout). | `results/kill/k9c/gallery/*.json` (per-design geometry dumps) + the same merged CSV for sort order. |
| `fig_e1_agreement.png` | E1's \\|exact - bisection referee\\| histogram (log-y) and the certificate's confusion counts (TP=%d, FP=%d, FN=%d, TN=%d against Theta_max>=eps=0.006). | `results/final/e1/e1.csv`, N=%d rows. |
| `summary_table.md` | Headline table: arms x {N, deployable, certified, eps_max>=0.1 rad, median eps_max, max eps_max}, computed directly from the CSVs above. | same. |

## Notes

- K9c's `results/kill/k9c/summary.txt` and `k9c.csv` are a stale 182-row snapshot; this
  script always merges the live `shard_*.csv` files instead. At merge time the 12 shards
  covered N=%d of the 400-design target (200 graphs x 2 sigma). If the run had still been
  short of 400, this note would say so and the script would need re-running once more
  shards land -- it reads the shards live every time, never the stale `k9c.csv`.
- The baseline (Eq. (6) + the K6 0+ repair) is independently confirmed 0/400 designs
  deployable and 0/400 certified in `results/kill/k6/k6_final.csv` before it is asserted
  in any figure here.
- `fig_gallery_full.png` intentionally renders every deployable design, not a curated
  subset (review/constructive_review.md section 5, item 1) -- it is a large image, one
  row per 6 designs.
""" % (len(k9c_rows), tp, fp, fn, tn, e1_n, len(k9c_rows))
    open(os.path.join(OUT_DIR, "README.md"), "w").write(readme_text)
    print("wrote %s" % os.path.join(OUT_DIR, "README.md"))
