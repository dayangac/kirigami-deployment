#!/usr/bin/env python3
"""A3 jitter-transition figure: deployable / certified fraction against jitter
amplitude, one curve per authored tiling per sigma rule, drawn from ladder.csv
written by `kill_jitter --analyze`. Run with `arch -arm64 /usr/local/bin/python3`."""
import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = sys.argv[1] if len(sys.argv) > 1 else "results/kill/jitter"

rows = list(csv.DictReader(open(os.path.join(D, "ladder.csv"))))
series = defaultdict(lambda: ([], [], [], []))  # (tiling, rule) -> a, f_cert, f_theta, med
nsplit = {}
for r in rows:
    if int(r["amp_idx"]) < 0:  # the a = 0 baseline is drawn as a marker, not a point
        continue
    k = (r["tiling"], r["rule"])
    s = series[k]
    s[0].append(float(r["amp"]))
    s[1].append(float(r["frac_cert"]))
    s[2].append(float(r["frac_theta_pos"]))
    s[3].append(float(r["median_theta"]) if r["median_theta"] not in ("", "nan") else float("nan"))
    if r["rule"] == "mc":  # |E_split| is sigma-dependent; label with the max-cut value
        nsplit[r["tiling"]] = int(r["n_split"])

tilings = sorted({t for t, _ in series})
# split-free tilings first, in a muted grey; split-bearing ones get the colour cycle
free = [t for t in tilings if nsplit[t] == 0]
bear = [t for t in tilings if nsplit[t] > 0]
cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
colour = {t: "0.65" for t in free}
for i, t in enumerate(bear):
    colour[t] = cycle[i % len(cycle)]

fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3))
panels = [
    (axes[0], 2, "fraction with exact $\\Theta_{max} > 0$"),
    (axes[1], 1, "fraction CERTIFIED (POS $\\wedge$ NOOVERLAP $\\wedge$ NOROOT, $\\epsilon$ = 0.006)"),
]
for ax, idx, title in panels:
    for t in free + bear:
        for rule, ls in (("mc", "-"), ("def", "--")):
            k = (t, rule)
            if k not in series:
                continue
            a, *cols = series[k]
            lbl = t if rule == "mc" else None
            ax.plot(a, cols[idx - 1], ls, color=colour[t], lw=1.8 if rule == "mc" else 1.2,
                    marker="o" if rule == "mc" else None, ms=3, label=lbl, alpha=0.95)
    ax.set_xscale("log")
    ax.set_xlabel("jitter amplitude $a$  (median edge lengths, interior vertices)")
    ax.set_ylim(-0.04, 1.04)
    ax.axhline(0.5, color="0.8", lw=0.8, zorder=0)
    ax.set_title(title, fontsize=9)
    ax.grid(alpha=0.25, lw=0.5)

ax = axes[2]
for t in free + bear:
    a, fc, ft, med = series[(t, "mc")]
    ax.plot(a, med, "-o", color=colour[t], lw=1.8, ms=3, label=t)
ax.set_xscale("log")
ax.set_xlabel("jitter amplitude $a$")
ax.set_title("median exact $\\Theta_{max}$ (rad), $\\sigma_{mc}$", fontsize=9)
ax.grid(alpha=0.25, lw=0.5)

axes[0].set_ylabel("fraction of 20 seeds")
h, l = axes[2].get_legend_handles_labels()
lab = [f"{t}  ($|E_{{split}}|$ = {nsplit[t]})" for t in l]
fig.legend(h, lab, loc="lower center", ncol=4, fontsize=8, frameon=False,
           bbox_to_anchor=(0.5, -0.02))
fig.suptitle("A3 jitter transition: authored tilings, combinatorics fixed, interior geometry "
             "perturbed\n(solid = $\\sigma$ max-cut, dashed = $\\sigma$ defect-minimising; "
             "grey = split-free patterns, whose Eq. (6) projection undoes the jitter exactly)",
             fontsize=10)
fig.tight_layout(rect=(0, 0.06, 1, 0.90))
out = os.path.join(D, "transition.png")
fig.savefig(out, dpi=160, bbox_inches="tight")
print("wrote", out)
