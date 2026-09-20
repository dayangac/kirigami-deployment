#!/usr/bin/env python3
"""WP2 -- why the K9c yield is not scale-free.

Reads only CSVs written by Kirigami/apps/exp_yield_features.jl:
  results/yield/features.csv  400 K9c designs, input-side features + the K9c outcomes
                              joined from results/experiments/k9c/k9c.csv
  results/yield/fresh.csv     held-out designs (ids 1000..), features and outcomes both
                              produced by the same binary via method::design_range_max

Writes:
  results/yield/analysis_out.txt      every number YIELD.md cites
  results/yield/fig_yield_predictor.png
  results/yield/fig_yield_matched.png
  results/yield/fig_yield_auc.png

No scikit-learn (not installed). AUC, Spearman, Wilson and the logistic fit are all
implemented here on numpy.

Run: arch -arm64 /usr/local/bin/python3 results/yield/analyse.py
"""
import csv
import itertools
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = open(os.path.join(HERE, "analysis_out.txt"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.write(s + "\n")
    OUT.flush()


# --- statistics --------------------------------------------------------------
def wilson(k, n, z=1.96):
    """Wilson score interval, identical to Kirigami/scripts/plot_common.jl's."""
    if n == 0:
        return 0.0, 0.0, 0.0
    phat = k / n
    denom = 1.0 + z * z / n
    center = (phat + z * z / (2 * n)) / denom
    half = z * math.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n)) / denom
    return phat, max(0.0, center - half), min(1.0, center + half)


def ranks(x):
    """Average ranks, ties shared."""
    x = np.asarray(x, float)
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x), float)
    sx = x[order]
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and sx[j + 1] == sx[i]:
            j += 1
        r[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return r


def auc(score, label):
    """Mann-Whitney U / (n_pos n_neg), ties counted at 0.5."""
    score = np.asarray(score, float)
    label = np.asarray(label, int)
    npos = int(label.sum())
    nneg = len(label) - npos
    if npos == 0 or nneg == 0:
        return float("nan")
    r = ranks(score)
    return (r[label == 1].sum() - npos * (npos + 1) / 2.0) / (npos * nneg)


def spearman(x, y):
    rx, ry = ranks(x), ranks(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    d = math.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float("nan") if d == 0 else float((rx * ry).sum() / d)


def logistic_fit(X, y, ridge=1e-3, iters=200):
    """IRLS / Newton with a ridge, on standardised columns. Returns (w, b, mu, sd)."""
    X = np.asarray(X, float)
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    Z = np.hstack([(X - mu) / sd, np.ones((len(X), 1))])
    w = np.zeros(Z.shape[1])
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-np.clip(Z @ w, -30, 30)))
        W = np.clip(p * (1 - p), 1e-9, None)
        g = Z.T @ (y - p) - ridge * np.r_[w[:-1], 0.0]
        H = (Z.T * W) @ Z + ridge * np.diag(np.r_[np.ones(Z.shape[1] - 1), 0.0])
        try:
            step = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            break
        w = w + step
        if np.max(np.abs(step)) < 1e-10:
            break
    return w, mu, sd


def logistic_predict(w, mu, sd, X):
    Z = np.hstack([(np.asarray(X, float) - mu) / sd, np.ones((len(X), 1))])
    return 1.0 / (1.0 + np.exp(-np.clip(Z @ w, -30, 30)))


# --- data --------------------------------------------------------------------
def load(path):
    with open(path) as f:
        return list(csv.DictReader(f))


feat_path = os.path.join(HERE, "features.csv")
fresh_path = os.path.join(HERE, "fresh.csv")
rows = load(feat_path)
fresh = load(fresh_path) if os.path.exists(fresh_path) else []

# Derived columns. Size-like counts span a decade, and a logistic link is linear in the
# feature, so their logs are added as candidates. Both are functions of the input only.
def add_derived(rs):
    for r in rs:
        for c in ("F", "N", "E_int", "n_incid", "n_constraints_0plus", "n_hinge",
                  "n_split"):
            r["log_" + c] = "%.10g" % math.log(max(1.0, float(r[c])))


add_derived(rows)
add_derived(fresh)

TEXT = {"id", "kind", "sigma", "best_src", "binding"}
OUTCOME = {"theta_exact", "eps_max", "best_feas", "best_margin", "dist", "secs"}
SOLVER = ["best_margin", "dist", "secs", "best_feas"]

all_cols = [c for c in rows[0].keys() if c not in TEXT]
FEATURES = [c for c in all_cols if c not in OUTCOME]
INPUT_SIDE = [c for c in FEATURES if not c.startswith("proj_")]
PROJ_SIDE = [c for c in FEATURES if c.startswith("proj_")]


def mat(rs, cols):
    return np.array([[float(r[c]) for c in cols] for r in rs], float)


def col(rs, c):
    return np.array([float(r[c]) for r in rs], float)


def ycol(rs):
    return (col(rs, "theta_exact") > 1e-9).astype(int)


say("WP2 yield analysis")
say("features.csv rows: %d   fresh.csv rows: %d" % (len(rows), len(fresh)))
y = ycol(rows)
p, lo, hi = wilson(int(y.sum()), len(y))
say("K9c population yield (theta_exact > 0): %d/%d = %.4f  Wilson95 [%.4f, %.4f]"
    % (y.sum(), len(y), p, lo, hi))
if fresh:
    yf = ycol(fresh)
    pf, lof, hif = wilson(int(yf.sum()), len(yf))
    say("FRESH yield (theta_exact > 0): %d/%d = %.4f  Wilson95 [%.4f, %.4f]"
        % (yf.sum(), len(yf), pf, lof, hif))
    say("FRESH eps_max >= 0.1: %d/%d" % (int((col(fresh, "eps_max") >= 0.1).sum()), len(yf)))

# ---------------------------------------------------------------------------
# 1. univariate table
say("")
say("=" * 78)
say("1. UNIVARIATE (400 K9c designs): AUC for theta_exact > 0, Spearman with eps_max")
say("   AUC > 0.5 means LARGER feature value -> MORE likely to deploy.")
say("=" * 78)
eps = col(rows, "eps_max")
uni = []
for c in FEATURES:
    v = col(rows, c)
    if np.all(v == v[0]):
        uni.append((c, float("nan"), float("nan"), 0.0))
        continue
    a = auc(v, y)
    uni.append((c, a, spearman(v, eps), abs(a - 0.5)))
uni.sort(key=lambda t: -t[3])
say("%-26s %8s %8s   %s" % ("feature", "AUC", "rho(eps)", "group"))
for c, a, s, _ in uni:
    grp = "projection" if c.startswith("proj_") else "input"
    say("%-26s %8.4f %8.4f   %s" % (c, a, s, grp))

say("")
say("solver-side (EXCLUDED from the pre-registered bar; restatements of the outcome):")
for c in SOLVER:
    v = col(rows, c)
    say("%-26s %8.4f %8.4f   solver" % (c, auc(v, y), spearman(v, eps)))

# ---------------------------------------------------------------------------
# 2. the size / family confound
EDGES = [100, 230, 420, 800]
BINLAB = ["F<230", "230<=F<420", "F>=420"]
FAMS = ["voronoi", "delaunay", "quad_random"]


def bin_of(F, edges=EDGES):
    for i in range(len(edges) - 1):
        if edges[i] <= F < edges[i + 1]:
            return i
    return 0 if F < edges[0] else len(edges) - 2


say("")
say("=" * 78)
say("2. THE SIZE / FAMILY CONFOUND")
say("=" * 78)
say("yield (theta_exact > 0) per family at matched |F| bins  [k/n, Wilson95]")
say("%-13s %s" % ("family", "  ".join("%-24s" % b for b in BINLAB)))
fam_bin = {}
for fam in FAMS:
    cells = []
    for b in range(3):
        sub = [r for r in rows if r["kind"] == fam and bin_of(float(r["F"])) == b]
        k = sum(1 for r in sub if float(r["theta_exact"]) > 1e-9)
        pp, l, h = wilson(k, len(sub))
        fam_bin[(fam, b)] = (k, len(sub), pp, l, h)
        cells.append("%3d/%-3d %.2f [%.2f,%.2f]" % (k, len(sub), pp, l, h))
    say("%-13s %s" % (fam, "  ".join("%-24s" % c for c in cells)))

# split-density bins, same treatment
sd_all = col(rows, "split_density")
sd_edges = list(np.quantile(sd_all, [0.0, 1 / 3, 2 / 3, 1.0]))
sd_edges[-1] += 1e-9
say("")
say("split_density terciles of the 400: edges %s" % ["%.4f" % e for e in sd_edges])
say("yield per family at matched split_density bins")
say("%-13s %s" % ("family", "  ".join("%-24s" % ("sd bin %d" % b) for b in range(3))))
for fam in FAMS:
    cells = []
    for b in range(3):
        sub = [r for r in rows if r["kind"] == fam
               and sd_edges[b] <= float(r["split_density"]) < sd_edges[b + 1]]
        k = sum(1 for r in sub if float(r["theta_exact"]) > 1e-9)
        pp, l, h = wilson(k, len(sub))
        cells.append("%3d/%-3d %.2f [%.2f,%.2f]" % (k, len(sub), pp, l, h))
    say("%-13s %s" % (fam, "  ".join("%-24s" % c for c in cells)))

say("")
say("family mean split_density (400 designs):")
for fam in FAMS:
    v = np.array([float(r["split_density"]) for r in rows if r["kind"] == fam])
    F = np.array([float(r["F"]) for r in rows if r["kind"] == fam])
    say("  %-13s split_density mean %.4f  median |F| %.0f  n=%d"
        % (fam, v.mean(), np.median(F), len(v)))
say("family medians of the leading features (400 designs):")
say("  %-13s %10s %10s %10s %10s %10s" % ("family", "min_mu_ini", "ang_min", "asp_med",
                                          "dof_ratio", "n_incid"))
for fam in FAMS:
    sub = [r for r in rows if r["kind"] == fam]
    say("  %-13s %10.3f %10.2f %10.3f %10.4f %10.0f"
        % (fam, np.median([float(r["min_mu_ini"]) for r in sub]),
           np.median([float(r["ang_min"]) for r in sub]),
           np.median([float(r["asp_med"]) for r in sub]),
           np.median([float(r["dof_ratio"]) for r in sub]),
           np.median([float(r["n_incid"]) for r in sub])))
say("yield in |F| bins, all families pooled:")
for b in range(3):
    sub = [r for r in rows if bin_of(float(r["F"])) == b]
    k = sum(1 for r in sub if float(r["theta_exact"]) > 1e-9)
    pp, l, h = wilson(k, len(sub))
    say("  %-12s %3d/%-3d %.3f [%.3f,%.3f]" % (BINLAB[b], k, len(sub), pp, l, h))

# ---------------------------------------------------------------------------
# 3. held-out logistic models
say("")
say("=" * 78)
say("3. LOGISTIC MODELS -- fit on one sigma half, evaluated on the other and on fresh")
say("   Model choice is made on the FIT set only; held-out AUC is reported, never used")
say("   to choose.")
say("=" * 78)

mc = [r for r in rows if r["sigma"] == "sigma_mc"]
de = [r for r in rows if r["sigma"] == "sigma_def"]
CAND = FEATURES  # every input-side and projection-side feature; no solver output


def cv_auc(fit_rows, cols, yv, folds=5, seed=0):
    """Mean out-of-fold AUC INSIDE the fit set. Selection uses this and nothing else."""
    n = len(fit_rows)
    rng = np.random.RandomState(seed)
    idx = rng.permutation(n)
    X = mat(fit_rows, cols)
    scores = np.zeros(n)
    for k in range(folds):
        te = idx[k::folds]
        tr = np.setdiff1d(idx, te)
        if len(np.unique(yv[tr])) < 2:
            return float("nan")
        w, mu, sd = logistic_fit(X[tr], yv[tr])
        scores[te] = logistic_predict(w, mu, sd, X[te])
    return auc(scores, yv)


def eval_split(fit_rows, ho_rows, fit_name, ho_name):
    yf_, yh = ycol(fit_rows), ycol(ho_rows)
    res = {}
    # -- best single feature, chosen by 5-fold CV AUC INSIDE the fit set
    singles = []
    for c in CAND:
        v = col(fit_rows, c)
        if np.all(v == v[0]):
            continue
        singles.append((cv_auc(fit_rows, [c], yf_), c))
    singles = [s for s in singles if not math.isnan(s[0])]
    singles.sort(reverse=True)
    best1 = singles[0][1]
    w, mu, sd = logistic_fit(mat(fit_rows, [best1]), yf_)
    p_fit = logistic_predict(w, mu, sd, mat(fit_rows, [best1]))
    p_ho = logistic_predict(w, mu, sd, mat(ho_rows, [best1]))
    res["single"] = (best1, auc(p_fit, yf_), auc(p_ho, yh), (w, mu, sd))
    say("[%s -> %s] best single feature by fit-set CV AUC: %s (CV %.4f)"
        % (fit_name, ho_name, best1, singles[0][0]))
    say("    in-sample fit AUC %.4f    held-out AUC %.4f"
        % (res["single"][1], res["single"][2]))
    say("    fit-set top 8 singles by CV AUC: %s"
        % ", ".join("%s %.3f" % (c, s) for s, c in singles[:8]))
    # -- best pair, chosen by the same CV, over the top 14 singles
    top = [c for _, c in singles[:14]]
    pairs = []
    for a_, b_ in itertools.combinations(top, 2):
        pairs.append((cv_auc(fit_rows, [a_, b_], yf_), a_, b_))
    pairs.sort(key=lambda t: -t[0])
    fa, fb = pairs[0][1], pairs[0][2]
    w2, mu2, sd2 = logistic_fit(mat(fit_rows, [fa, fb]), yf_)
    p_ho2 = logistic_predict(w2, mu2, sd2, mat(ho_rows, [fa, fb]))
    res["pair"] = ((fa, fb), auc(logistic_predict(w2, mu2, sd2, mat(fit_rows, [fa, fb])),
                                 yf_), auc(p_ho2, yh), (w2, mu2, sd2))
    say("[%s -> %s] best pair by fit-set CV AUC: (%s, %s) (CV %.4f)"
        % (fit_name, ho_name, fa, fb, pairs[0][0]))
    say("    in-sample fit AUC %.4f    held-out AUC %.4f"
        % (res["pair"][1], res["pair"][2]))
    say("    runner-up pairs (CV): %s"
        % "; ".join("(%s,%s) %.3f" % (a_, b_, s) for s, a_, b_ in pairs[1:4]))
    return res


res_mc = eval_split(mc, de, "fit sigma_mc", "held-out sigma_def")
say("")
res_de = eval_split(de, mc, "fit sigma_def", "held-out sigma_mc")

say("")
say("--- additional fit: ALL 400 K9c designs, selected by 5-fold CV on those 400 only,")
say("    with the FRESH graphs as the held-out set (the strongest legitimate test: the")
say("    fresh graphs are disjoint from the K9 population, fresh_check.txt) ---")
y400 = ycol(rows)
sing400 = sorted(((cv_auc(rows, [c], y400), c) for c in CAND
                  if not np.all(col(rows, c) == col(rows, c)[0])), reverse=True)
sing400 = [s for s in sing400 if not math.isnan(s[0])]
b1 = sing400[0][1]
w1a, mu1a, sd1a = logistic_fit(mat(rows, [b1]), y400)
top400 = [c for _, c in sing400[:14]]
pr400 = sorted(((cv_auc(rows, [a_, b_], y400), a_, b_)
                for a_, b_ in itertools.combinations(top400, 2)), reverse=True)
pa, pb = pr400[0][1], pr400[0][2]
w2a, mu2a, sd2a = logistic_fit(mat(rows, [pa, pb]), y400)
res_all = {"single": (b1, auc(logistic_predict(w1a, mu1a, sd1a, mat(rows, [b1])), y400),
                      float("nan"), (w1a, mu1a, sd1a)),
           "pair": ((pa, pb),
                    auc(logistic_predict(w2a, mu2a, sd2a, mat(rows, [pa, pb])), y400),
                    float("nan"), (w2a, mu2a, sd2a))}
say("    best single by CV on the 400: %s (CV %.4f, in-sample %.4f)"
    % (b1, sing400[0][0], res_all["single"][1]))
say("    top 8 singles by CV: %s" % ", ".join("%s %.3f" % (c, s) for s, c in sing400[:8]))
say("    best pair by CV on the 400: (%s, %s) (CV %.4f, in-sample %.4f)"
    % (pa, pb, pr400[0][0], res_all["pair"][1]))
say("    runner-up pairs (CV): %s"
    % "; ".join("(%s,%s) %.3f" % (a_, b_, s) for s, a_, b_ in pr400[1:4]))

if fresh:
    say("")
    say("--- the models applied to the FRESH graphs (%d rows) ---" % len(fresh))
    yfr = ycol(fresh)
    summary_fresh = []
    for nm, res in (("fit sigma_mc", res_mc), ("fit sigma_def", res_de),
                    ("fit all 400", res_all)):
        c1, _, _, (w, mu, sd) = res["single"]
        p1 = logistic_predict(w, mu, sd, mat(fresh, [c1]))
        a1 = auc(p1, yfr)
        (ca, cb), _, _, (w2, mu2, sd2) = res["pair"]
        p2 = logistic_predict(w2, mu2, sd2, mat(fresh, [ca, cb]))
        a2 = auc(p2, yfr)
        say("[%s] single %-22s fresh AUC %.4f   predicted yield %.4f (obs %.4f)"
            % (nm, c1, a1, p1.mean(), yfr.mean()))
        say("[%s] pair   (%s, %s)  fresh AUC %.4f   predicted yield %.4f"
            % (nm, ca, cb, a2, p2.mean()))
        summary_fresh.append((nm, c1, a1, p1.mean(), (ca, cb), a2, p2.mean()))
    pf, lof, hif = wilson(int(yfr.sum()), len(yfr))
    say("observed fresh yield %.4f  Wilson95 [%.4f, %.4f]" % (pf, lof, hif))
    say("")
    say("PASS BAR: held-out AUC >= 0.85 on fresh AND predicted yield inside the Wilson"
        " interval.")
    best_auc = max(max(t[2], t[5]) for t in summary_fresh)
    ok = False
    for nm, c1, a1, y1, cp, a2, y2 in summary_fresh:
        for lab, a_, yp in (("single " + c1, a1, y1),
                            ("pair " + str(cp), a2, y2)):
            inside = lof <= yp <= hif
            say("  %-14s %-40s AUC %.4f  pred %.4f  in-interval %s  -> %s"
                % (nm, lab, a_, yp, inside,
                   "PASS" if (a_ >= 0.85 and inside) else "FAIL"))
            ok = ok or (a_ >= 0.85 and inside)
    say("VERDICT: %s   (best fresh AUC obtained: %.4f)"
        % ("PASS" if ok else "FAIL", best_auc))

# ---------------------------------------------------------------------------
# 4. binding modes of the failures vs the top feature
TOP = res_all["single"][0]   # the univariate / CV winner on the whole 400
say("")
say("=" * 78)
say("4. BINDING MODE OF THE FAILURES vs %s" % TOP)
say("=" * 78)
tv = col(rows, TOP)
q = np.quantile(tv, [0.25, 0.5, 0.75])
say("%s quartiles over the 400: %s" % (TOP, ["%.4f" % x for x in q]))
modes = sorted({r["binding"] for r in rows})
for md in modes:
    sub = [r for r in rows if r["binding"] == md]
    v = np.array([float(r[TOP]) for r in sub])
    sdv = np.array([float(r["split_density"]) for r in sub])
    say("  %-14s n=%3d   median %s %.4f   median split_density %.4f   median |F| %.0f"
        % (md, len(sub), TOP, np.median(v), np.median(sdv),
           np.median([float(r["F"]) for r in sub])))
if fresh:
    say("fresh binding modes: "
        + ", ".join("%s=%d" % (m, sum(1 for r in fresh if r["binding"] == m))
                    for m in sorted({r["binding"] for r in fresh})))

# ---------------------------------------------------------------------------
# 5. mechanism probes: is it the constraint count, the null-space dimension, or the
#    extreme value of a min over a growing number of constraints?
say("")
say("=" * 78)
say("5. MECHANISM PROBES")
say("=" * 78)
lF = np.log(col(rows, "F"))


def slope(x, yv):
    A = np.vstack([x, np.ones(len(x))]).T
    b = np.linalg.lstsq(A, yv, rcond=None)[0]
    pred = A @ b
    ss = 1 - ((yv - pred) ** 2).sum() / ((yv - yv.mean()) ** 2).sum()
    return b[0], ss


for nm in ("n_constraints_0plus", "n_incid", "n_split", "n_corner", "dim_null", "N"):
    sl, r2 = slope(lF, np.log(np.maximum(1.0, col(rows, nm))))
    say("  log %-20s ~ %.3f log|F|   (R2 %.3f)" % (nm, sl, r2))
say("  2*dim_null vs n_constraints_0plus: dof_ratio median %.4f, Spearman with |F| %.4f"
    % (np.median(col(rows, "dof_ratio")), spearman(col(rows, "F"),
                                                   col(rows, "dof_ratio"))))
say("  => dof_ratio is %s in |F|; a 'constraints outgrow the null space' story needs a"
    % ("essentially flat" if abs(spearman(col(rows, "F"), col(rows, "dof_ratio"))) < 0.3
       else "trending"))
say("     clearly negative trend here.")
say("")
say("  Spearman(min_mu_ini, log|F|) = %.4f ; Spearman(min_q_ini, log|F|) = %.4f"
    % (spearman(col(rows, "min_mu_ini"), lF), spearman(col(rows, "min_q_ini"), lF)))
say("  median min_mu_ini by |F| bin:")
for b in range(3):
    sub = [r for r in rows if bin_of(float(r["F"])) == b]
    say("    %-12s n=%3d  median min_mu_ini %.5f  median min_q_ini %.5f  median n_incid %.0f"
        % (BINLAB[b], len(sub), np.median([float(r["min_mu_ini"]) for r in sub]),
           np.median([float(r["min_q_ini"]) for r in sub]),
           np.median([float(r["n_incid"]) for r in sub])))
say("")
say("  residual power: AUC of each candidate WITHIN matched |F| bins (pooled over bins)")
for c in ("min_mu_ini", "min_q_ini", "proj_min_q_x0", "ang_min", "split_density",
          "dof_ratio", "asp_med"):
    aucs, ns = [], []
    for b in range(3):
        sub = [r for r in rows if bin_of(float(r["F"])) == b]
        a_ = auc(np.array([float(r[c]) for r in sub]), ycol(sub))
        if not math.isnan(a_):
            aucs.append(a_)
            ns.append(len(sub))
    say("    %-18s within-bin AUC %s   n-weighted %.4f"
        % (c, ["%.3f" % a_ for a_ in aucs],
           float(np.average(aucs, weights=ns))))
say("  and AUC of |F| WITHIN terciles of min_mu_ini:")
mq = np.quantile(col(rows, "min_mu_ini"), [0, 1 / 3, 2 / 3, 1.0])
mq[-1] += 1e-12
aucs, ns = [], []
for b in range(3):
    sub = [r for r in rows if mq[b] <= float(r["min_mu_ini"]) < mq[b + 1]]
    a_ = auc(np.array([float(r["F"]) for r in sub]), ycol(sub))
    aucs.append(a_)
    ns.append(len(sub))
say("    F within min_mu_ini terciles: %s   n-weighted %.4f"
    % (["%.3f" % a_ for a_ in aucs], float(np.average(aucs, weights=ns))))

# ---------------------------------------------------------------------------
# 6. does the picture reproduce on the fresh graphs?
if fresh:
    say("")
    say("=" * 78)
    say("6. REPRODUCTION ON THE FRESH GRAPHS")
    say("=" * 78)
    yfr = ycol(fresh)
    say("  fresh yield in |F| bins:")
    for b in range(3):
        sub = [r for r in fresh if bin_of(float(r["F"])) == b]
        k = sum(1 for r in sub if float(r["theta_exact"]) > 1e-9)
        pp, l, h = wilson(k, len(sub))
        say("    %-12s %3d/%-3d %.3f [%.3f,%.3f]" % (BINLAB[b], k, len(sub), pp, l, h))
    say("  fresh yield per family:")
    for fam in FAMS:
        sub = [r for r in fresh if r["kind"] == fam]
        k = sum(1 for r in sub if float(r["theta_exact"]) > 1e-9)
        pp, l, h = wilson(k, len(sub))
        say("    %-13s %3d/%-3d %.3f [%.3f,%.3f]" % (fam, k, len(sub), pp, l, h))
    say("  univariate AUC on the 400 vs on the 200 fresh, top features:")
    say("    %-24s %8s %8s" % ("feature", "AUC 400", "AUC fresh"))
    for c in ("min_mu_ini", "log_F", "min_q_ini", "n_incid", "proj_min_q_x0",
              "n_constraints_0plus", "ang_min", "hole_size_mean", "mean_split_comp",
              "dof_ratio", "split_density"):
        say("    %-24s %8.4f %8.4f"
            % (c, auc(col(rows, c), y), auc(col(fresh, c), yfr)))
    say("  Spearman(min_mu_ini, log|F|) on fresh = %.4f"
        % spearman(col(fresh, "min_mu_ini"), np.log(col(fresh, "F"))))

# ---------------------------------------------------------------------------
# figures
COL_400 = "#2b4f81"
COL_FRESH = "#c1440e"


def fig_predictor():
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.0))
    tv400 = col(rows, TOP)
    ax = axes[0]
    # binned yield with Wilson bars
    for rs, cl, lab in ((rows, COL_400, "K9c 400"),
                        (fresh, COL_FRESH, "fresh %d" % len(fresh)) if fresh else (None,) * 3):
        if not rs:
            continue
        v = np.array([float(r[TOP]) for r in rs])
        yy = ycol(rs)
        edges = np.quantile(col(rows, TOP), np.linspace(0, 1, 7))
        edges[-1] += 1e-9
        xs, ps, lo_, hi_ = [], [], [], []
        for i in range(len(edges) - 1):
            sel = (v >= edges[i]) & (v < edges[i + 1])
            if sel.sum() < 3:
                continue
            pp, l, h = wilson(int(yy[sel].sum()), int(sel.sum()))
            xs.append(v[sel].mean())
            ps.append(pp)
            lo_.append(pp - l)
            hi_.append(h - pp)
        ax.errorbar(xs, ps, yerr=[lo_, hi_], fmt="o-", color=cl, ms=4, capsize=2,
                    lw=1.2, label=lab)
    ax.set_xlabel(TOP)
    ax.set_ylabel(r"yield: exact $\Theta_{max} > 0$ (Wilson 95%)")
    ax.set_ylim(-0.03, 1.03)
    ax.legend(fontsize=8, frameon=False)
    ax.set_title("yield vs %s" % TOP, fontsize=10)

    for ax, field, lab in ((axes[1], "theta_exact", r"$\Theta_{max}$ (rad)"),
                           (axes[2], "eps_max", r"$\varepsilon_{max}$ (rad)")):
        ax.scatter(tv400, col(rows, field), s=9, alpha=0.55, color=COL_400,
                   label="K9c 400", edgecolors="none")
        if fresh:
            ax.scatter([float(r[TOP]) for r in fresh], col(fresh, field), s=13,
                       alpha=0.8, color=COL_FRESH, label="fresh", edgecolors="none")
        ax.set_xlabel(TOP)
        ax.set_ylabel(lab)
        ax.legend(fontsize=8, frameon=False)
        ax.set_title("%s vs %s" % (lab, TOP), fontsize=10)
    fig.suptitle("WP2: the top input-side predictor of K9c yield", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(HERE, "fig_yield_predictor.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    say("wrote %s" % out)


def fig_matched():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), sharey=True)
    colors = {"voronoi": "#2b4f81", "delaunay": "#1a7a4c", "quad_random": "#c1440e"}
    ax = axes[0]
    for fam in FAMS:
        xs, ps, lo_, hi_ = [], [], [], []
        for b in range(3):
            k, n, pp, l, h = fam_bin[(fam, b)]
            if n == 0:
                continue
            xs.append(b)
            ps.append(pp)
            lo_.append(pp - l)
            hi_.append(h - pp)
        ax.errorbar(xs, ps, yerr=[lo_, hi_], fmt="o-", color=colors[fam], ms=5,
                    capsize=3, lw=1.3, label=fam)
    ax.set_xticks(range(3))
    ax.set_xticklabels(BINLAB)
    ax.set_ylabel(r"yield: exact $\Theta_{max}>0$ (Wilson 95%)")
    ax.set_title("matched on |F|", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.set_ylim(-0.03, 1.03)

    ax = axes[1]
    for fam in FAMS:
        xs, ps, lo_, hi_ = [], [], [], []
        for b in range(3):
            sub = [r for r in rows if r["kind"] == fam
                   and sd_edges[b] <= float(r["split_density"]) < sd_edges[b + 1]]
            if len(sub) < 3:
                continue
            k = sum(1 for r in sub if float(r["theta_exact"]) > 1e-9)
            pp, l, h = wilson(k, len(sub))
            xs.append(b)
            ps.append(pp)
            lo_.append(pp - l)
            hi_.append(h - pp)
        ax.errorbar(xs, ps, yerr=[lo_, hi_], fmt="s--", color=colors[fam], ms=5,
                    capsize=3, lw=1.3, label=fam)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["sd<%.3f" % sd_edges[1], "%.3f-%.3f" % (sd_edges[1], sd_edges[2]),
                        "sd>=%.3f" % sd_edges[2]], fontsize=8)
    ax.set_title("matched on split_density", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle("Family yield at matched bins (400 K9c designs)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(HERE, "fig_yield_matched.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    say("wrote %s" % out)


def fig_auc():
    top = [t for t in uni if not math.isnan(t[1])][:20][::-1]
    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    names = [t[0] for t in top]
    vals = [t[1] for t in top]
    cols = ["#8a8f98" if n.startswith("proj_") else "#2b4f81" for n in names]
    ax.barh(range(len(top)), [v - 0.5 for v in vals], left=0.5, color=cols)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(names, fontsize=8)
    ax.axvline(0.5, color="k", lw=0.8)
    ax.set_xlabel(r"AUC for exact $\Theta_{max}>0$ (400 K9c designs)")
    ax.set_title("Top 20 features by |AUC - 0.5|\n"
                 "blue = input-side, grey = Eq. (6) projection-side", fontsize=10)
    fig.tight_layout()
    out = os.path.join(HERE, "fig_yield_auc.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    say("wrote %s" % out)


say("")
fig_predictor()
fig_matched()
fig_auc()
OUT.close()
