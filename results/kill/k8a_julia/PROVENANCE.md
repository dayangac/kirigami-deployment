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

## LP layer: which K8a columns are reproducible across implementations (2026-09-20)

Probe by port-method-3 on bit-identical inputs (`lp_probe/`: `probe_k1a.jl`,
`probe_inputs.jl`, `probe.cpp` linked against `code/build/libkiri_core.a`, Xcode clang++
arm64 -O2; outputs `julia_*.txt` / `cpp_*.txt`). For each (cut structure, X) both sides
compute the flex basis `N`, the cone rows `A`, `M = A N`, the row norms of `M` BEFORE
normalisation, the sigma chart, the LP on that chart at 600 and at 20000 Frank-Wolfe
iterations, and the same LP with the near-zero rows (norm < 1e-12) dropped before
normalising ("clean"). `A` and the flex projector `N N^T` agree between C++ and Julia to
1e-15 on every input tested; the differences below are entirely downstream of them.

**Mechanism 1 -- noise rows (not reproducible).** `expansive_cone.cpp::normalise_rows`
skips only rows of norm <= 0. A row whose exact value on the flex space is zero (two
copies of a vertex that never separate at first order) arrives as ~1e-16 rounding noise
and is scaled to a unit row whose direction is the rounding of the flex basis (Eigen
Householder vs LAPACK). Counts of such rows (norm < 1e-12, out of n_rows):

| input | noise rows | in the sigma chart |
|---|---|---|
| voronoi_0 / X_ini | 5095 / 5169 | 2739 |
| delaunay_1 / X_ini | 80 / 1641 | 40 |
| quad_random_2 / X_ini | 6594-6635 / 8905 | ~3370 |
| voronoi_3 / X_ini | 6138-6165 / 6201 | ~3310 |
| delaunay_4 / X_ini | 212-214 / 5966 | ~106 |
| trunc_square_R20 / R35 / R40, truncated_square_488 | 20 / 56 / 9 / 56 | 1-2 |
| snub_square_R20 | 36 / 196 | 2-4 |
| t3_4_3_12_R25, _s0 | 2 / 154 | 0-1 |
| every other deployable row and reference tiling probed | 0 | 0 |

(At a non-deployable X_ini the flex space is essentially the rigid motions, so nearly
every row vanishes on it; the C++/Julia counts differ by a few rows sitting at the 1e-12
threshold.) Consequences, confirmed by dropping those rows:

* `sigma_chart` on the K1a X_ini rows: C++ -0.8186 / -0.2357 / -0.7333 / -0.8349 / -0.0498
  vs Julia -0.3677 / -0.1314 / -0.6385 / -0.3715 / -0.0629 (voronoi_0, delaunay_1,
  quad_random_2, voronoi_3, delaunay_4). With the noise rows dropped both give
  -0.166538 / -0.131376 / -0.300918 / -0.0531013 / -0.0497877. So the documented
  "solver-free" quantity is solver-free but not noise-free; on those rows it is a random
  number, and its sign carries no information. On every deployable row (0 noise rows in
  the chart, or the noise rows not attaining the min) `sigma_chart` agrees to 1e-9.
* `margin_l2` / `margin_inf` / `dual_bound` / `n_active` / `n_dual_support` on the
  truncated-square family and the other rows above: closed brackets on DIFFERENT
  matrices, so e.g. trunc_square_R20 0.9006 (C++) vs 0.9574 (Julia), trunc_square_R35
  0.8085 vs 0.7735, trunc_square_R40 0.7016 vs 0.7423, truncated_square_488 0.8523 vs
  0.8953, snub_square_R20 0.1646 vs 0.1694, t3_4_3_12_R25 0.2701 vs 0.2833 -- unchanged
  at 20000 iterations. With the noise rows dropped: 1.0 / 1.0 / 0.74235 / 1.0 / 0.169401 /
  0.283291 on both sides. The feasibility flags (`pass_feasible`, sign of `margin_l2`)
  were never affected on the rows tested, but a noise row CAN in principle flip them.

**Mechanism 2 -- Frank-Wolfe stopping path (reproducible at convergence only).** The
drivers run `--dual-iters 600`; on the larger snub_square / kagome / snub_square_33434
rows the bracket is still open at 600 (lp_gap 1e-3 .. 5e-2; kagome_R40 margin 0.078 vs
dual 0.128) and the two implementations stop at different points of the path:
margins differ by 0.1-2 % (snub_square_R35_s0 0.05376 vs 0.05482, snub_square_33434
0.05954 vs 0.05913). At 20000 iterations both close the bracket and agree to ~1e-7
(kagome_R40 0.115938, snub_square_33434 0.0657838, snub_square_R30 0.100594, ...).
`dual_bound`, `dual_max`, `farkas_resid` and `lp_gap` on the infeasible K1a rows
(1e-10 .. 1e-4) and the tiny `farkas_sumerr` values are likewise stopping-path numbers.

**Verdict per column** (C++ vs Julia, bit-identical inputs):

| column | reproducible? |
|---|---|
| N, F, n_split, n_corner, n_convex, n_rows, dim_ker_A, dim_flex, components, flex_resid, sigma_resid, sigma_is_flex, sigma_min_q, sigma_min_mu, sigma_bad_*, sigma_in_cone, sigma_span, euler_eps | yes (all 82 rows) |
| pass_feasible, sign of margin_l2, cert_1e9 | yes on all 82 rows, but not guaranteed on inputs with noise rows |
| sigma_chart | yes on the 77 deployable rows (1e-9); NOT on the 5 K1a X_ini rows (noise rows) |
| margin_l2, margin_inf, dual_bound, dual_max, n_active, n_dual_support | yes to ~1e-6 on rows without noise rows once the bracket is closed (20000 iterations); at the driver's 600 iterations only to 0.1-2 % on the large snub/kagome rows; NOT on the truncated-square family, snub_square_R20, t3_4_3_12_R25 (noise rows, 3-6 %) |
| farkas_resid, farkas_lmin, farkas_sumerr, lp_gap, cert_1e6 | stopping-path numbers; only their order of magnitude is comparable |

**Intent and fix (not applied).** The C++ header says "rows of norm <= 0 are left
alone and reported through n_degenerate"; an exactly-zero row left in the matrix would
force the LP margin to 0, so the degenerate case was never handled consistently.
A row that is identically zero on the flex space is no constraint and should be dropped
before normalising (e.g. norm < 1e-9 max_i |a_i N|); after that every column above is
reproducible to solver tolerance. The Julia port stays faithful so that the columns
quoted from the C++ remain comparable; the docstring of `normalise_rows!` carries this
verdict.
