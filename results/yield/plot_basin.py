# WP2b figure. Reads results/yield/basin.csv (written by Kirigami/apps/kill_basin.jl and
# merged by Kirigami/apps/kill_basin_agg.jl) and writes results/yield/fig_basin.png.
# Run as: arch -arm64 /usr/local/bin/python3 results/yield/plot_basin.py
# Plotting only -- every number quoted in BASIN.md comes from basin_stats.txt.
import csv
import collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = "results/yield/basin.csv"
OUT = "results/yield/fig_basin.png"
EPS = 1e-9

rows = list(csv.DictReader(open(CSV)))
D = collections.OrderedDict()
for r in rows:
    key = (int(r["id"]), r["sigma"])
    D.setdefault(key, {"cls": r["orig_class"], "kind": r["kind"], "seeds": {}})
    D[key]["seeds"][int(r["k"])] = float(r["theta_exact"])

fails = [(k, v) for k, v in D.items() if v["cls"] == "fail"]
succs = [(k, v) for k, v in D.items() if v["cls"] == "success"]

# failures sorted by best Theta_max over the eight new seeds, descending
fails.sort(key=lambda kv: -max([t for j, t in kv[1]["seeds"].items() if j > 0] or [0.0]))
# successes sorted by median Theta_max over the nine seeds
def median(a):
    a = sorted(a)
    n = len(a)
    return 0.0 if n == 0 else (a[n // 2] if n % 2 else 0.5 * (a[n // 2 - 1] + a[n // 2]))
succs.sort(key=lambda kv: -median(list(kv[1]["seeds"].values())))

fig, ax = plt.subplots(2, 1, figsize=(11, 8.2))

# ---- panel 1: the 93 originally failing designs ------------------------------------
a = ax[0]
nflip = 0
for i, (key, v) in enumerate(fails):
    new = [t for j, t in v["seeds"].items() if j > 0]
    a.scatter([i] * len(new), new, s=9, color="#4477aa", alpha=0.75, zorder=3,
              label="one of the 8 new seeds" if i == 0 else None)
    a.scatter([i], [v["seeds"].get(0, 0.0)], s=14, marker="_", color="#cc3311", zorder=4,
              label="archived K9c seed (k = 0)" if i == 0 else None)
    if max(new or [0.0]) > EPS:
        nflip += 1
a.axhline(0, color="0.6", lw=0.7)
a.set_xlabel("originally failing design, sorted by the best of its 8 new seeds "
             "(%d designs)" % len(fails))
a.set_ylabel(r"exact $\Theta_{\max}$  [rad]")
a.set_title("A. The 93 K9c failures re-seeded: %d of %d flip (at least one new seed "
            "deploys)" % (nflip, len(fails)))
a.legend(loc="upper right", frameon=False, fontsize=9)
a.grid(axis="y", alpha=0.25)

# ---- panel 2: the sampled successes ------------------------------------------------
b = ax[1]
for i, (key, v) in enumerate(succs):
    th = sorted(v["seeds"].values())
    n = len(th)
    q25 = th[int(0.25 * (n - 1))]
    q75 = th[min(n - 1, int(0.75 * (n - 1) + 0.999))]
    b.vlines(i, th[0], th[-1], color="0.65", lw=1.0, zorder=2)
    b.vlines(i, q25, q75, color="#4477aa", lw=4.0, zorder=3)
    b.scatter([i], [v["seeds"].get(0, 0.0)], s=18, marker="x", color="#cc3311", zorder=4,
              label="archived K9c seed (k = 0)" if i == 0 else None)
b.set_xlabel("sampled successful design, sorted by median $\\Theta_{\\max}$ over its 9 "
             "seeds (%d designs)" % len(succs))
b.set_ylabel(r"exact $\Theta_{\max}$  [rad]")
b.set_title("B. Seed-to-seed spread on the successes: min-max bar, IQR box, "
            "the archived seed marked")
b.legend(loc="upper right", frameon=False, fontsize=9)
b.grid(axis="y", alpha=0.25)

fig.tight_layout()
fig.savefig(OUT, dpi=160)
print("wrote", OUT, "-- %d failures, %d successes" % (len(fails), len(succs)))
