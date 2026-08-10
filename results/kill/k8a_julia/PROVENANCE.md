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
| n_convex | 0 / 82 | 9 / 82 | 9 / 82 |
| sigma_in_cone | 0 / 82 | 5 / 82 | 5 / 82 |
| sign of sigma_min_q | 0 / 82 | 5 / 82 | 5 / 82 |
| sigma_bad_q, sigma_bad_mu | 0 / 82 | 6 / 82 | 6 / 82 |
| pass_feasible, sign of margin_l2 | 0 / 82 | 4 / 82 | 4 / 82 |
| euler_eps | 1 / 82 | 11 / 82 | 10 / 82 |

**The archived combinatorial columns are reproduced exactly by the current C++.** The
Julia differences are the same set in both comparisons and all 11 rows are shape-space
SAMPLES of the deployable population (`snub_square_R20_s1`, `snub_square_R40_s0`,
`t3_4_3_12_R{25,30,35,40}_s0`, `t3_4_3_12_R40_s1`, `trunc_square_R30_s{0,1}`,
`trunc_square_R40_s{0,1}`). Every base configuration (`*_R*` without `_sN`), every
reference tiling and every K1a graph agrees on all of these columns.

Cause of the sample rows: a sample is `X = X0 + Phi * T` with `T` drawn from the shared
`mt19937(20260903)` stream (reproduced bit-exactly), but `Phi` is the null-space basis
returned by `solve_system`, which is only defined up to an orthogonal change of basis;
Julia's LAPACK SVD returns a different basis than Eigen's `JacobiSVD`, so the same `T`
is a different point of the same shape space. The base rows do not involve `Phi` and
match. (On the x86_64/Rosetta Julia used earlier, two truncated-square BASES also
differed because the orientation relaxation's trig landed on a different sigma; on arm64
all bases match the C++ and the frozen `deployable_population.json`.) The `(2, 0.2)`
population is not frozen in `data/corpus` -- only the `(8, 0.2)` one is -- so these
sample rows can only be made bit-comparable by freezing that population from the C++.

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
from a fresh `kill_k8a` run (C++) -- the Julia LP (`expansive_cone.jl`) currently agrees
with the fresh C++ on feasibility flags but not to better than a few percent on the
margin/dual values, and `sigma_chart` (documented as solver-free) differs, which is an
open item for the method port.

Files: fresh C++ run at the session scratchpad `cpp_k8a/k8a.csv` (re-runnable with
`cd code/build && ./kill_k8a --n 5 --out <dir>`); comparison script `threeway.py` there.
