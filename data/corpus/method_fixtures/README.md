# method_fixtures/ -- frozen inputs and C++ expected numbers of the test_method.cpp cases

Each `<test>.json` replays the inputs of a group of `code/tests/test_method.cpp` cases
exactly (same generator calls, same `std::mt19937` seeds, same sigma method, same Eq. (6)
projection) and stores, next to the inputs, every number the C++ computed. Producer: the
`freeze_<test>.cpp` next to it, compiled with the Xcode clang (arm64, `-O2 -std=c++20`,
default `-ffp-contract=on`) against `code/build/libkiri_core.a` on 2026-09-19. Julia tests:
`Kirigami/test/<test>.jl`.

Index base inside the files is the C++ one (0-based faces, M'-vertices, face pairs,
`bad_*` fields; `-1` = none); the Julia loaders convert at the boundary.

| file | C++ cases replayed | contents |
|---|---|---|
| `test_method_1.json` | the `make_case` tilings of the deploy_basis / contact / certificate cases (10 tilings: `squares_checker`, `triangles_checker`, `kagome_checker`, `snub_{2_6,3_0,3_2}`, `trunc_{3_0,4_0}`, `hexagons_3_0`, `t3_4_3_12_4_2`), the `delaunay_of_random_points` meshes of the two mobility cases (`rng(31337)`, `rng(4711)`), the range-objective case (`trunc_3_0`, `rng(5)` normal perturbation) | per tiling: mesh + sigma + X, deploy basis C/S, swept radii, pruned pair list, `exact_theta_max` report, `exact_theta_max_overlap` report with C(X), bisection and shipped theta_max, `validity_certificate` at eps = 0.006 (both overloads) and eps = 3.0 with the first_root probes, `contact_angles(0, 0.006)`; per mobility mesh: hinge graph, X0, `MobilityReport` at X_ini and X0, rigidity rank; range objective: Phi, t, f0, gradient, the 12 central-difference probes, a 2x5 `maximize_range` run |

Provenance strings (`_provenance`, per-block `provenance`) name the exact C++ calls.

Other files in this directory belong to sibling ports and document themselves:
`test_method_2.json` / `freeze_method_2.cpp` (port-method-2), `method_3_reference.json`
(port-method-3), `design_intermediates_148.json`, `k9c_calls_*.json`, `freeze_design.cpp`,
`freeze_k9c.cpp` with `README_design.md` (port-method-4).

## `test_method_2.json` (test_method_2.jl: zero_plus / convex_embed; test_range_embed.jl)

Produced by `freeze_method_2.cpp` (same toolchain and linkage as above, run 2026-09-19; it
makes the SAME calls as `tests/test_method.cpp::make_case` / `tests/test_range_embed.cpp::
make_rcase` and the test cases that follow). One entry per reference tiling under `cases`
(`hexagons_3.0`, `snub_square_3.2`, `truncated_square_{4.0,3.0}`, `3_4_3_12_4.2`,
`kagome_2.5_checker`):

| key | contents |
|---|---|
| `mesh`, `X0`, `X_case`, `Phi`, `dim_null`, `med` | mesh with sigma (`mt19937(2026)` relaxation or checkerboard), the Eq. (6) projection, `make_case`'s X, the C++ (Eigen JacobiSVD) null basis, the median edge |
| `q_at_case`, `mu_at_case`, `crosses_at_case`, `margin_*_at_X0`, `modes_at_X0` | `zero_plus_q`, `zero_plus_corner_margin`, `corner_crosses`, `zero_plus_margin`, `range_embed_modes` |
| `probe_t`, `*_at_probe` | at `t_i = 0.02 med sin(i+1)`: `zero_plus_eval`, `zero_plus_objective` (lambda 1e-3, w_corner 1, w_prox 0.5), `convex_embed_objective` (a: convexity only; b: + split 1e-3), `range_embed_objective` (delta 1e-4, kappa 2e-2, bw 0 and 1e-3) with values and gradients |
| `convex_embed_solve`, `range_embed_solve`, `zero_plus_repair_solve` | the full solves with the options of the corresponding C++ test cases |

The Julia null basis (LAPACK SVD) spans the same space as Eigen's (projector difference
~1e-16) but its columns differ, so the C++ tests that draw random `t` in null-space
coordinates are replayed on the frozen `Phi`. Exact quantities, objectives and gradients
agree with the C++ to ~1e-14; the L-BFGS solves agree to ~1e-8 except where the path is
chaotic (snub_square range_embed: 1e-4 after 3 x 60 iterations, see test_range_embed.jl).

## `method_3_reference.json` (test_method_3.jl: periodic_jacobian / budget / expansive_cone)

Produced by a scratchpad freezer (`freeze_method_3.cpp`, linked against
`code/build/libkiri_core.a`, Apple clang arm64 `-O2 -std=c++20`, `-ffp-contract=on`) that
makes the SAME calls as `tests/test_method.cpp`; nothing is invented. Indices are 0-based
as the C++ wrote them; meshes are `mesh_to_json_string` output.

| key | contents |
|---|---|
| `lattice` | `detect_lattice` on the seven `make_tiling_pattern` big tilings: t1, t2 (the libc++ `std::sort` tie order decides the representative) |
| `c4` | the two K7 rows of the C4 test: cell (with sigma), T, P0, Q, K, theta_c, face_area |
| `periodic_budget` | the five 2x2 families of the periodic budget test: sigma, counts, W, R, det P0 tr K |
| `patch_budget` | `make_budget_case` on five generator kinds: mesh+sigma, X, W, R, bW, bR, border functional, void area at 0.7 |
| `cone` | `make_case` tilings: mesh+sigma, X, flex/cone dimensions, LP margins, sigma chart margin |
| `cone_lp_1200x60` | the random 1200x60 LP: witness, margin, dual bound, gap, first row sample |

Design intermediates (`k9c_calls_*.json`, `design_intermediates_148.json`): see `README_design.md`.
