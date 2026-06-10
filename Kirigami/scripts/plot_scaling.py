#!/usr/bin/env python3
"""fig_scaling.png + the fitted exponents of SCALING.md, from results/scaling/scaling.csv.

Left panel: log-log wall time against |F| for each routine, with an ordinary least-squares
fit of log10(t) on log10(|F|). The reported exponent is that slope; its 95 % interval is
slope +- t(0.975, n-2) * SE(slope), with t read from a small table (numpy only, no scipy;
1.96 is used beyond df = 30). Rows whose wall time is below the 1 ms clock resolution of
the driver are excluded from the fit and drawn hollow, because log10(0) is not a datum.

Right panel: peak resident set (getrusage ru_maxrss at process exit) against |F|.

Run under `arch -arm64 /usr/local/bin/python3` (D3).
"""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SC = os.path.join(ROOT, "results", "scaling")
TCRIT = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
         8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145,
         15: 2.131, 16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093, 20: 2.086, 21: 2.080,
         22: 2.074, 23: 2.069, 24: 2.064, 25: 2.060, 26: 2.056, 27: 2.052, 28: 2.048,
         29: 2.045, 30: 2.042}
ROUTINES = [("solve_dense", "assemble + solve_system (dense SVD)", "#b03030"),
            ("rank_sparse", "rank_only_sparse (sparse QR)", "#c07a20"),
            ("characterize", "characterize (exact scan + certificate)", "#1f6f4a"),
            ("design_range_max", "design_range_max (incl. its own solve)", "#3050b0")]
MIN_T = 0.002   # the driver's clock is milliseconds; below this the datum is quantised


def load():
    with open(os.path.join(SC, "scaling.csv")) as f:
        return [r for r in csv.DictReader(f)]


def fit(F, t):
    """OLS slope of log10(t) on log10(F); returns (slope, lo, hi, n, r2)."""
    x, y = np.log10(np.asarray(F, float)), np.log10(np.asarray(t, float))
    n = len(x)
    if n < 3:
        return (math.nan,) * 3 + (n, math.nan)
    A = np.vstack([x, np.ones_like(x)]).T
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    dof = n - 2
    s2 = float(resid @ resid) / dof
    cov = s2 * np.linalg.inv(A.T @ A)
    se = math.sqrt(cov[0, 0])
    tc = TCRIT.get(dof, 1.96)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - float(resid @ resid) / ss_tot if ss_tot > 0 else math.nan
    return beta[0], beta[0] - tc * se, beta[0] + tc * se, n, r2


def main():
    rows = load()
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.0))
    lines = []
    for key, label, color in ROUTINES:
        sel = [r for r in rows if r["routine"] == key and r["status"] == "ok"]
        if not sel:
            continue
        F = np.array([float(r["F"]) for r in sel])
        t = np.array([float(r["secs"]) for r in sel])
        rss = np.array([float(r["rss_mb"]) for r in sel])
        good = t >= MIN_T
        axes[0].plot(F[~good], np.maximum(t[~good], 1e-4), "o", ms=4, mfc="none",
                     color=color, alpha=0.6)
        axes[0].plot(F[good], t[good], "o", ms=5, color=color)
        s, lo, hi, n, r2 = fit(F[good], t[good])
        if not math.isnan(s):
            xx = np.array([F[good].min(), F[good].max()])
            c = np.log10(t[good]).mean() - s * np.log10(F[good]).mean()
            axes[0].plot(xx, 10 ** (c + s * np.log10(xx)), "-", lw=1.6, color=color,
                         label=r"%s: $\alpha$ = %.2f [%.2f, %.2f], n=%d, $R^2$=%.3f"
                               % (label, s, lo, hi, n, r2))
        lines.append((key, label, s, lo, hi, n, r2, int(F.max()), float(t.max()),
                      float(rss.max())))
        axes[1].plot(F, rss, "o-", ms=4, lw=1.2, color=color, label=label)
    for a in axes:
        a.set_xscale("log")
        a.set_yscale("log")
        a.set_xlabel(r"faces $|F|$")
        a.grid(alpha=0.25, which="both")
    axes[0].set_ylabel("wall time [s]")
    axes[0].set_title(r"wall time vs. $|F|$, fitted $t \propto |F|^{\alpha}$", fontsize=10)
    axes[0].legend(fontsize=7.5, loc="upper left")
    axes[1].set_ylabel("peak RSS [MB]  (ru_maxrss at exit)")
    axes[1].set_title("peak resident set vs. $|F|$", fontsize=10)
    axes[1].legend(fontsize=7.5, loc="upper left")
    fig.tight_layout()
    out = os.path.join(SC, "fig_scaling.png")
    fig.savefig(out, dpi=170)
    print("wrote", out)

    with open(os.path.join(SC, "fits.csv"), "w") as f:
        f.write("routine,label,exponent,ci_lo,ci_hi,n_fit,r2,max_F_completed,"
                "max_secs,max_rss_mb\n")
        for r in lines:
            f.write("%s,%s,%.4f,%.4f,%.4f,%d,%.4f,%d,%.4f,%.2f\n" % r)
    print("wrote", os.path.join(SC, "fits.csv"))
    for r in lines:
        print(r)


if __name__ == "__main__":
    main()
