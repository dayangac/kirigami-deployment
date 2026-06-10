#!/usr/bin/env python3
# E1 -- histogram of |Theta_max(T4.2") - bisect| over the extended population.
# Run: arch -arm64 /usr/local/bin/python3 code/scripts/plot_e1.py
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = "results/final/e1/e1.csv"
OUT = "results/final/e1/e1_gap.png"

gaps = []
with open(CSV) as f:
    for row in csv.DictReader(f):
        try:
            gaps.append(float(row["gap"]))
        except ValueError:
            continue

fig, ax = plt.subplots(figsize=(7, 4.5))
import math
nz = [g for g in gaps if g > 0]
bins = [0] + list(sorted(set([10 ** (e / 4.0) for e in range(-4 * 16, 4 * 2)])))
ax.hist(gaps, bins=60, color="#3b6ea5", edgecolor="white", linewidth=0.3)
ax.set_yscale("log")
ax.set_xlabel(r"$|\Theta_{max}(\mathrm{T4.2''}) - \Theta_{max}(\mathrm{bisection})|$ (rad)")
ax.set_ylabel("count (log)")
ax.set_title(f"E1: exact vs. bisection referee, N={len(gaps)}")
ax.axvline(1e-5, color="crimson", linestyle="--", linewidth=1, label="1e-5 agreement bar")
ax.legend()
fig.tight_layout()
fig.savefig(OUT, dpi=150)
print(f"wrote {OUT}, N={len(gaps)}, worst={max(gaps) if gaps else float('nan')}")
