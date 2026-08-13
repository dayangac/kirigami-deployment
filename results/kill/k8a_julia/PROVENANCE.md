# K8a -- provenance of the three CSVs (archived C++ / fresh C++ / Julia)

Three-way comparison on the 82 rows of `k8a_julia/k8a.csv` (arm64 Julia 1.12.7,
`kill_k8a.jl --limit 5`: the 4 split-bearing reference tilings, the 73 configs of
`deployable_population(2, 0.2)`, and K1a graphs 0..4 at X_ini), against

* **archived**: `results/kill/k8a/k8a.csv` (277 rows, merged shards; committed 2026-09-04 in
  `dd74944` of kirigami-experiments, 31 columns);
* **fresh C++**: `code/build/kill_k8a --n 5` run 2026-09-19 with the current binary
  (`kill_k8a.cpp` / `expansive_cone.cpp` at `d38b0a4`, 39 columns).

Rows were matched on (`id`, `where`); every Julia row has a partner in both C++ files.

## Combinatorial columns

| column | archived != fresh C++ | archived != Julia | fresh C++ != Julia |
|---|---|---|---|
| N, F, n_split, n_corner, n_rows, dim_ker_A, dim_flex, components, sigma_is_flex | 0 / 82 | 0 / 82 | 0 / 82 |
| n_convex, sigma_in_cone, sign of sigma_min_q, sigma_bad_q, sigma_bad_mu, pass_feasible, sign of margin_l2 | 0 / 82 | 0 / 82 | 0 / 82 |
| euler_eps | 1 / 82 | 1 / 82 | 0 / 82 |

**Every combinatorial / geometric column agrees three ways** (the one `euler_eps` row,
`snub_square_R40_s1`, is archived 0.2 vs 0.5 in both the fresh C++ and Julia: an LP-witness
dependent quantity). On the 73 deployable rows Julia also matches the fresh C++ at rtol
1e-9 on everything except the LP solver outputs (`margin_l2`, `margin_inf`, `dual_bound`,
`dual_max`, `n_active`, `n_dual_support`, `farkas_*`, `lp_gap`); `sigma_chart` matches
there too.

History of this table: a first Julia run rebuilt `deployable_population(2, 0.2)` in Julia
and differed from both C++ files on 11 shape-space SAMPLE rows (`n_convex` 9, `sigma_in_cone`
5, ...). A sample is `X = X0 + Phi * T` with `T` from the bit-exact `mt19937(20260903)`
stream but `Phi` the null-space basis of `solve_system`, which Julia's LAPACK SVD returns in
a different orthogonal frame than Eigen's `JacobiSVD`; the base rows never involve `Phi`
and always matched. Since 2026-09-20 the C++ Eigen-basis samples are frozen as
`data/corpus/deployable_population_{2_0.2,20_0.2,20_0.35}.json` and `kill_k8a.jl` /
`kill_e1.jl` read them (`--regenerate` rebuilds in Julia and brings the sample-row
differences back). With the frozen files the E1 authored rows also match the fresh C++
`kill_e1` at rtol 1e-9 (shards 0 and 1 of 200: 15 + 22 rows, all columns).
(On the x86_64/Rosetta Julia used even earlier, two truncated-square BASES differed as
well, through the orientation relaxation's trig; on arm64 all bases match.)

## Solver columns: what the archived CSV does NOT reproduce

Archived vs fresh C++ differ on `margin_l2` (44/82), `margin_inf` (43), `dual_bound`
(30), `dual_max` (30), `n_active` (41), `n_dual_support` (55), and the fresh file carries 8
columns the archived one lacks (`sigma_chart`, `sigma_span`, `farkas_resid`, `farkas_lmin`,
`farkas_sumerr`, `cert_1e9`, `cert_1e6`, `lp_gap`). Typical: `delaunay_1/X_ini`
`dual_max` archived `9.999e-08` (the old stopping threshold) vs fresh `4.448e-03`;
`triangles_R40` `margin_l2` `5.698e-03` vs `5.095e-02`.

Explanation from the C++ history: the archived CSV was produced by the binary at
`dd74944` (2026-09-04, "Round 2 kills ... K8a FAIL (dual-certified)"). The next commit
`d38b0a4` (2026-09-04, "Hostile reviews ..., STATE/KILL_REPORT wording corrections")
changed `code/src/method/expansive_cone.cpp` (+151 lines) and `kill_k8a.cpp` (+38) and
added `kill_k8a_recheck.cpp`: the Frank-Wolfe dual loop's stopping rule was replaced
(`if (gn < opt.dual_tol) break;` -> a duality-gap test `gap_fw <= gap_tol`, with
`dual_tol` now applied together with `gap_fw <= dual_tol^2`), the sigma-chart lower
bound and the independent Farkas re-verification were added, and the driver gained the 8
recheck columns. This is the "the solver was wrong; the certificates were not" recheck
documented in `results/kill/k8a/recheck.md` and KILL_REPORT.md section "K8a -- Recheck".
`results/kill/k8a/k8a.csv` was not regenerated after that commit, so its LP values are
those of the pre-recheck solver; the archived `recheck_k9*.csv` files are post-`d38b0a4`.

Consequence for the paper: any number quoted from `results/kill/k8a/k8a.csv` that is a
COUNT of combinatorial facts (sigma in P(X), n_convex, dim_flex, components, split
counts) is reproduced by both the current C++ and Julia on the bases, reference tilings
and K1a graphs. Any quoted LP value (margins, dual bounds, "dual-certified" tallies,
n_active / n_dual_support) comes from the superseded solver and must be re-derived
from a fresh `kill_k8a` run (C++) -- the Julia LP (`expansive_cone.jl`) agrees with the
fresh C++ on feasibility flags and on `sigma_chart` for the deployable rows, but not to
better than a few percent on the margin/dual values, and on the K1a graphs at X_ini
`sigma_chart` (documented as solver-free) still differs (e.g. voronoi_3 -0.83 vs -0.57),
which is an open item for the method port.

Files: fresh C++ run at the session scratchpad `cpp_k8a/k8a.csv` (re-runnable with
`cd code/build && ./kill_k8a --n 5 --out <dir>`); comparison script `threeway.py` there.
