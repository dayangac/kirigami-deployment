# reference_patterns/ -- stored pattern files the tests and reference apps read

Archived outputs of `kiri_reference` (Phase 2, Gate 2, 2026-09-03) and of one derivation
test, kept verbatim (same basenames). None of these is a `make_graph` population.

| File | Origin | Consumed by | What it pins |
|---|---|---|---|
| `reference_cases.json` | `results/core_validation/reference_cases.json`, written by `kiri_reference` (`write_json(jrows, outdir * "/reference_cases.json")`) | `Kirigami/test/test_reference_cases.jl` (both testsets: "reference_cases.json is present and complete" and the per-case regression) | 8 rows `{name, N, F, H, H_geo, interior, hinge, split, rank_L, dim_null, res_before, res_after, min_beta, theta_max, theta_max_opt, claim, claim_applicable, claim_ok}`. The test asserts `F`, `split`, `dim_null` exactly and cross-checks `theta_max` only at 3e-6 (it is a pre-F34 bisection value; see the test header). The meshes themselves are NOT in this file -- the test rebuilds them from the seeded generator calls; the frozen copies are `../reference_cases_8.json`. |
| `cases/<name>/M.json` | `results/core_validation/cases/<name>/M.json` (`mesh_json(g)`: uncut graph with `orientation`) | `Kirigami/apps/kiri_sweep.jl` (`--cases`, default `../results/core_validation/cases`); not read by `Pkg.test()` | Input mesh of each of the 8 cases with its sigma. Should equal `../reference_cases_8.json[i].mesh` (see verification note below). |
| `cases/<name>/X0.json` | same dir, `mesh_json(m0)` with `X = X0` | `kiri_sweep` | Eq. (6) projected flat state X0 of each case. |
| `cases/<name>/deploy_30deg.json`, `deploy_60deg.json`, `deploy_thetamax.json` | same dir, `deployment_json(c, X0, theta)` | `Kirigami/scripts/plot_embedding.jl` (figures only) | Deployed vertex positions `Y`, `prime_faces`, geometric hole cycles at theta = 30deg, 60deg and 0.999 theta_max. Useful as end-to-end regression targets for `deploy()`. |
| `T5_1_positive_orientation_zero_range_hexagons.json` | `derivations/check_failures/`, written by `dump_failure()` in the derivation tests | nothing reads it; it is the saved failing case of the T5.1 check (hexagons_auto: positive orientation but zero range) | The hexagons_auto mesh + sigma that violated derivation claim T5.1 -- a known counterexample worth keeping as a test fixture. |

Case names (`cases/`): `rotating_squares`, `triangles_alternating`, `kagome_3636`,
`hexagons_auto`, `truncated_square_488`, `snub_square_33434`, `tiling_3_4_3_12`,
`periodic_squares_4x4`.

Other test fixtures are NOT files: `test_design.jl` rebuilds its named K9/K9c designs
with `make_graph(id, 100, 800, 1400)` (rows of `../k1a_200.json` / `../native200.json`),
`derivation_tests.jl` builds its L1 corpus from the 7 non-periodic reference cases plus
`make_graph(id, 18, 46, 220)` (`../derivation_l1_small.json`), and the remaining tests call
the tiling generators directly (`tiling_squares(rect(...))`, `generate(kind, {R}, rng)`).

Verification performed when freezing (2026-09-19): `reference_cases.json` rows agree with
`../reference_cases_8.json` on `name`, `N`, `F` for all 8 cases; `cases/*/M.json` vertices
and faces are compared against `../reference_cases_8.json` in `../README.md`.

Unit-test reference fixtures (`test_fixtures_*.json`: mesh_cut, holes, system, kinematics,
collision, rank_checks) replay the procedurally built test cases with their seeds and store
the reference numbers next to the inputs. See `../README.md` for the per-file contents.
