#!/usr/bin/env python3
"""Summarise results/kill/k8a/k8a.csv (written by code/apps/kill_k8a.cpp).

Run as:  arch -arm64 /usr/local/bin/python3 code/scripts/summarise_k8a.py results/kill/k8a
No numpy / matplotlib: this is a text reduction of the merged CSV, printed to stdout so
the caller can tee it into summary.txt.
"""
import collections
import csv
import sys

d = sys.argv[1] if len(sys.argv) > 1 else "results/kill/k8a"
rows = list(csv.DictReader(open(d + "/k8a.csv")))
f = lambda r, k: float(r[k])
med = lambda v: sorted(v)[len(v) // 2]

print("K8a -- expansive-cone LP")
print("configurations in k8a.csv:", len(rows))

ref = [r for r in rows if r["kind"] == "reference"]
dep = [r for r in rows if r["kind"].startswith("deployable_")]

print("\n-- HARD SOUNDNESS CONTROL: the four split-bearing authored tilings")
for r in ref:
    print(
        f"   {r['id']:22s} F={r['F']:>3} split={r['n_split']:>3} dim_ker_A={r['dim_ker_A']:>3}"
        f" dim_flex={r['dim_flex']:>3}  sigma_in_cone={r['sigma_in_cone']}"
        f"  min q={r['sigma_min_q']}  min mu={r['sigma_min_mu']}"
        f"  LP margin={r['margin_l2']}  euler_eps={r['euler_eps']}"
    )
print(
    "   sigma in P(X): %d/%d ; LP margin > 0: %d/%d"
    % (
        sum(int(r["sigma_in_cone"]) for r in ref),
        len(ref),
        sum(1 for r in ref if f(r, "margin_l2") > 1e-8),
        len(ref),
    )
)

print("\n-- DEPLOYABLE POPULATION (kill_common::deployable_population)")
print(
    "   configs %d ; sigma in P(X) %d ; LP margin > 0 %d ; Euler step tested %d, collision-free %d"
    % (
        len(dep),
        sum(int(r["sigma_in_cone"]) for r in dep),
        sum(1 for r in dep if f(r, "margin_l2") > 1e-8),
        sum(int(r["euler_tested"]) for r in dep),
        sum(1 for r in dep if f(r, "euler_eps") > 0),
    )
)
print(
    "   euler_eps histogram (median-edge units):",
    dict(collections.Counter(r["euler_eps"] for r in dep if int(r["euler_tested"]))),
)
for r in dep:
    if not int(r["sigma_in_cone"]):
        print(
            f"   sigma OUTSIDE P(X): {r['id']}  min q={r['sigma_min_q']}"
            f"  min mu={r['sigma_min_mu']}  LP margin={r['margin_l2']}"
        )

for w in ("X_ini", "X0"):
    S = [r for r in rows if r["where"] == w]
    if not S:
        continue
    print(f"\n-- K1a POPULATION at {w}  (n={len(S)})")
    print(
        f"   LP margin > 0                                     : "
        f"{sum(1 for r in S if f(r,'margin_l2')>1e-8)}/{len(S)}   (PASS bar: >= 20/100)"
    )
    for t in (1e-6, 1e-5):
        print(
            f"   dual_max < {t:g} (every branch chart tried certified) : "
            f"{sum(1 for r in S if f(r,'dual_max')<t)}/{len(S)}"
        )
    print(
        f"   dual_max: median {med([f(r,'dual_max') for r in S]):.2e}"
        f"  max {max(f(r,'dual_max') for r in S):.2e}"
    )
    print(f"   sigma is a flex of the framework : {sum(int(r['sigma_is_flex']) for r in S)}/{len(S)}")
    print(f"   sigma in P(X)                    : {sum(int(r['sigma_in_cone']) for r in S)}/{len(S)}")
    print(
        f"   dim ker A median {med([int(r['dim_ker_A']) for r in S])},"
        f" dim flex median {med([int(r['dim_flex']) for r in S])},"
        f" c(Gamma) values {sorted(set(int(r['components']) for r in S))}"
    )
    print(
        f"   rows median {med([int(r['n_rows']) for r in S])}"
        f" (split {med([int(r['n_split']) for r in S])},"
        f" corner incidences {med([int(r['n_corner']) for r in S])},"
        f" all convex: {sum(int(r['n_convex']) for r in S)==sum(int(r['n_corner']) for r in S)})"
    )
    for kind in ("voronoi", "delaunay", "quad_random"):
        T = [r for r in S if r["kind"] == kind]
        if not T:
            continue
        print(
            f"     {kind:12s} n={len(T):3d} dim_ker_A med {med([int(r['dim_ker_A']) for r in T]):4d}"
            f"  margin>0 {sum(1 for r in T if f(r,'margin_l2')>1e-8):3d}"
            f"  dual_max<1e-6 {sum(1 for r in T if f(r,'dual_max')<1e-6):3d}"
            f"  dual_max med {med([f(r,'dual_max') for r in T]):.1e}"
        )

print("\n-- SELF-CHECKS")
print("   max |R N| / max|R| over every configuration : %.2e" % max(f(r, "flex_resid") for r in rows))
print("   passes solved per configuration            :",
      dict(collections.Counter(int(r["passes"]) for r in rows)))
