# data/corpus -- frozen inputs of the C++ experiments

Every input population that the C++ reference repo (`kirigami-experiments`) builds
procedurally, written to JSON once so the Julia port can rerun each experiment on
byte-identical inputs and so a bit-exact `std::mt19937` port can be tested.

Producer: `code/apps/freeze_corpus.cpp` in `kirigami-experiments` (build dir
`code/build-freeze`, target `freeze_corpus`), run 2026-09-19 from the repo root with
default paths. The app calls the SAME functions with the SAME arguments as the owning
experiment app (table below); it invents nothing. Archived orientations are read from the
same directories the apps read them from. Files in `reference_patterns/` are plain copies
(see `reference_patterns/INDEX.md`).

Toolchain that produced the numbers: Apple clang 16.0.0 (clang-1600.0.26.3), arm64,
`-O2 -std=c++20`, libc++ `_LIBCPP_VERSION = 180100`, Eigen 5.0.1, nlohmann-json 3.12.
Doubles are serialised by nlohmann-json (shortest round-trip representation), so parsing
them back gives the exact binary value.
Floating-point contraction: clang's arm64 default `-ffp-contract=on` was in effect, so
`a*b+c` expressions in libc++'s distributions and in the generators compile to fused
multiply-adds. The Julia port's bit-exact RNG test (port-mt19937, 427/427 values) matches
only under FMA; rebuilding `freeze_corpus` with `-ffp-contract=off` would change the
`uniform_real_*` / `normal_0_1` vectors and, through them, every random mesh here.

## Population files

All population files are a JSON **array**, one object per graph, ids in ascending order.
Every random graph is a deterministic function of its `id` through
`kill::make_graph(id, min_faces, max_faces, n_cap)` (`code/apps/kill_common.hpp`):

* `kind = {"voronoi","delaunay","quad_random"}[id % 3]`
* `rng = std::mt19937(1000003*id + 20260903)` (stored as `seed`)
* `target = min_faces * (max_faces/min_faces)^U(0,1)`, sites = `max(20, target/2)` (delaunay),
  `max(20, target*1.15)` (voronoi), `max(20, target*1.6)` (quad_random); up to 8 attempts,
  `generate(kind, {sites, 40.0}, rng)` on the SAME rng, accept iff
  `min_faces <= F <= max_faces` and `N <= n_cap`, then
  `sigma = assign_orientation_relaxation(m, rng, 4, 300, 90).sigma` (Eq. (1) max-cut relaxation).

Common fields of a graph object:

| field | meaning |
|---|---|
| `id`, `kind`, `ok` | graph id, family, whether `make_graph` succeeded (rows with `ok=false` have an empty mesh and are kept so the id -> row map is total; only `k2c_500.json` has any, 38 of them) |
| `make_graph_args` | `[min_faces, max_faces, n_cap]` used |
| `seed` | the `std::mt19937` seed `1000003*id + 20260903` |
| `N`, `F`, `E`, `n_interior` | vertex / face / edge / interior-vertex counts after `build_topology()` |
| `mesh` | `mesh_to_json_string` output parsed: `{vertices: [[x,y],...], faces: [[v,...],...], orientation: [...]}`; `orientation` is `sigma_mc` |
| `sigma_mc` | the Eq. (1) orientation `make_graph` leaves in `mesh.sigma` (+1 clockwise, -1 counter-clockwise) |
| `sigma_checker` | `kill::checkerboard(m)`: 2-colouring of the dual graph by DFS from face 0 (what K2c uses as its alternative sigma); computed for every graph |
| `sigma_def` / `sigma_def_source` | K5's defect-minimising orientation, or `null`. Source is the archived file it was read from, or the call that recomputed it |
| `gidx` | (native200, k2b only) the population index the sharded apps use |

| file | graphs | call | sigma variants | owning apps / tests | notes |
|---|---|---|---|---|---|
| `k1a_200.json` | 200 (ids 0..199, all ok) | `make_graph(id, 100, 800, 1400)` | `sigma_mc`, `sigma_checker`, `sigma_def` (all 200, from `results/kill/k5/sigma/<kind>_<id>.json`) | `kill_k1a` (ids 0..199), `kill_k1b`, `kill_k8a` (ids 0..99), `kill_yield --mode check`, `test_design.cpp` (`k9_design(id)`, `k9c_row(id)`) | F in [101, 793], N in [56, 1338]. Verified bit-identical (vertices, faces, sigma_def) to the 200 archived K5 files. |
| `native200.json` | 200 (ids 0..199, `gidx` 0..199) | `make_graph(id, 100, 800, 1400)`; walk ids 0..399, keep iff `ok` AND an archived K5 sigma_def matches, stop at 200 | same as k1a | `kill_k5`, `kill_k6`, `kill_k9`, `kill_k9b`, `kill_k9c`, `kill_t1`, `kill_b4`, `kill_k8a_recheck`, `kill_native200`, `native200_cell`, `kill_yield --mode k9c`, `kill_basin` (subset via k9c.csv) | Turns out to be the SAME 200 graphs as `k1a_200.json` (every id 0..199 builds and has a sigma_def), so `gidx == id`. Kept as its own file because the apps define the population by the filtered walk. |
| `k3a_500.json` | 500 (ids 0..499, all ok) | `make_graph(id, 100, 5000, 1 << 30)` | `sigma_mc`, `sigma_checker` | `kill_k3a`, `kill_k3a_recheck` (15 of these ids) | No vertex cap. F in [100, 4918], N up to 9838 (212 graphs have N > 1400). (id, N, F) match the archived `results/kill/k3a/k3a.csv` for all 500 ids. 165 MB. |
| `k2c_500.json` | 500 rows, **462 ok** (ids 0..499; 38 ids have `ok=false` because all 8 attempts exceeded the 1400-vertex cap -- `kill_k2c` skips them with `if (!g.ok) continue`) | `make_graph(id, 100, 5000, 1400)` | `sigma_mc`, `sigma_checker` (K2c's second sigma) | `kill_k2c` | Same seeds as K3a; the graph differs from `k3a_500.json` on exactly the 212 ids where the K3a graph has N > 1400 (cap forced a retry), and is identical on the other 288. F in [100, 2775]. The 462 ok ids and their (N, F) match the archived `results/kill/k2c/k2c.csv` (462 rows) exactly. |
| `regime_100.json` | 100 (ids 2000..2099) | `make_graph(id, 20, 100, 1400)` | `sigma_mc`, `sigma_checker`, `sigma_def` = `method::orientation_defect(m0, sigma_mc, 20*F, 7000 + id).sigma` (recomputed exactly as `kill_regime.cpp` does; extra fields `sigma_def_ok`, `sigma_def_defect`) | `kill_regime` | All 100 defect searches succeeded (`sigma_def_ok = true`); sigma_def differs from sigma_mc on all 100. F in [24, 100]. |
| `yield_fresh_100.json` | 100 (ids 1000..1099, all ok) | `make_graph(id, 100, 800, 1400)` | `sigma_mc`, `sigma_checker`, `sigma_def` from `results/yield/fresh_sigma/<kind>_<id>.json` (written by `kill_yield --mode fresh` with `orientation_defect(m0, sigma_mc, 20*F, 7000 + id)`); plus `sigma_def_recomputed`, `sigma_def_recompute_ok`, `sigma_def_recompute_matches` from rerunning that call here | `kill_yield --mode fresh` (held-out set) | F in [101, 769]. Recompute reproduced the archived sigma_def on 100/100 graphs (so `orientation_defect` is deterministic given the RNG; 0.2-2 s per graph). |
| `e1_900.json` | 900 (ids 0..899) | `make_graph(id, 100, 800, 1400)` | `sigma_mc`, `sigma_checker`; `sigma_def` only for ids 0..199 (K5 archive) | `kill_e1` part (B) | E1's own sigma_def uses `defect_search_at(m, m.X, sigma_mc, 20*F, 5000000 + id)`, a function local to `kill_e1.cpp` -- NOT frozen (see "not frozen"). |
| `k2b_random_8.json` | 8 (`gidx` 0..7 = ids 0..7) | `make_graph(id, 100, 160, 1400)`, first 8 of ids 0..399 that build | `sigma_mc`, `sigma_checker` | `kill_k2b` | |
| `derivation_l1_small.json` | 80 (ids 0..79, all ok) | `make_graph(id, 18, 46, 220)` | `sigma_mc`, `sigma_checker` | `tests/derivation_tests.cpp` `l1_corpus()` (uses these after the 7 non-periodic reference cases until 60 cases) | F in [21, 46]. |
| `derivation_l1_140.json` | 140 (ids 0..139, all ok) | `make_graph(id, 18, 46, 220)` | `sigma_mc`, `sigma_checker` | `derivations/scratch/check_l1` (`check_l1 140 20`; 84 of these survive its own filter) | Same call; rows 0..79 are byte-identical to `derivation_l1_small.json`. 1.0 MB. |
| `scaling_42.json` | 42 cells | `kill_scaling::build_graph(kind, sites, seed, true)`: `rng = mt19937(1000003*seed + 20260908 + 7919*sites)`, `generate(kind, {sites, 40.0}, rng)`, then `assign_orientation_relaxation(m, rng, 4, 300, 90)`; grid of `code/scripts/run_scaling.sh`: targets {50,100,200,500,1000,2000,5000}, delaunay sites = target/2, voronoi sites = target, seeds 0,1,2 | `sigma_mc` | `kill_scaling` | Fields `kind, target_faces, sites, seed, rng_seed, N, F, E, mesh, sigma_mc`. N, F agree with all 39 cells present in the archived `results/scaling/scaling.csv`. |
| `reference_cases_8.json` | 8 | `kill::reference_cases()`: the Phase-2 tilings with ONE shared `mt19937(20260903)` for the four relaxation calls (`assign_orientation_relaxation(g, rng, 8, 500, 180)`, in file order), checkerboard for the rest | `sigma` (per case) | `kiri_reference`, `test_reference_cases.cpp`, `test_design.cpp`, `derivation_tests.cpp` | Fields `name, periodic, N, F, mesh, sigma`. `periodic_squares_4x4` carries `mesh.periodic`. Verified equal to `reference_patterns/cases/<name>/M.json` for all 8. |
| `deployable_population.json` | 187 configs (35 bases + 152 shape-space samples) | `kill::deployable_population(8, 0.2)`: 7 tiling families x radii {2.0,2.5,3.0,3.5,4.0}, sigma checkerboard (squares, triangles, kagome) or `assign_orientation_relaxation(m, mt19937(777 + 13*ri + 101), 8, 500, 180)`, X = Eq. (6) projection (dense SVD), samples from `mt19937(20260903)` normal draws filtered for no inverted face and no collision | `mesh.orientation` | `kill_k1a`, `kill_k2b`, `kill_e1` part (A) | Fields `name, family, dim_null, sample (-1 = base), mesh, X`. NOTE: unlike the others this contains solver OUTPUT (`X`), so a port should match it to ~1e-9, not bit-exactly. |
| `deployable_population_2_0.2.json` | 73 configs (35 bases + 38 samples) | `kill::deployable_population(2, 0.2)` | as above | `kill_k8a` (its 73 "deployable_*" rows) | Bases identical to `deployable_population.json` (35/35). Samples are NOT a prefix of the (8, 0.2) samples beyond the first base: `rng0` is one shared stream and each base consumes a number of normal draws that depends on `samples_per_base`, so every (samples, radius) pair is its own population and must be frozen from C++. 0.7 MB. |
| `deployable_population_20_0.2.json` | 415 configs (35 bases + 380 samples) | `kill::deployable_population(20, 0.2)` | as above | `kill_e1` part (A), tag `a1` | Bases identical to the (8, 0.2) file (35/35). 3.9 MB. |
| `deployable_population_20_0.35.json` | 415 configs (35 bases + 380 samples) | `kill::deployable_population(20, 0.35)` | as above | `kill_e1` part (A), tag `a2` | Bases identical to the (8, 0.2) file (35/35). 3.9 MB. |

## `mt19937_vectors.json`

For seeds `[0, 1, 5489, 20260903, 27260924]` (the last is `1000003*7 + 20260903`, i.e.
`make_graph` id 7), each drawn from a FRESH engine and a fresh distribution object:

* `raw_u32`: first 20 `operator()` outputs (seed 5489 starts `3499211612, 581869302, ...`).
* `uniform_real_0_1`, `uniform_real_0_40`, `uniform_real_0_2pi`: first 20 draws of
  `std::uniform_real_distribution<double>(a, b)`, as JSON numbers and as `%.17g` strings
  (`*_str`). libc++ builds these from `generate_canonical<double, 53>` = two 32-bit draws
  per double, so a port that uses one draw per double will NOT match.
* `uniform_int_0_99`, `uniform_int_0_1`: first 20 draws of `std::uniform_int_distribution<int>`
  (libc++ rejection sampling; `(0,1)` is what `orientation_defect` uses for random starts).
* `normal_0_1`: first 20 draws of `std::normal_distribution<double>(0,1)` (libc++ Marsaglia
  polar method with a cached second value; used by `deployable_population`).
* `shuffle_0_19`: `std::shuffle` of `[0..19]` from a fresh engine (libc++ uses
  `uniform_int_distribution` on the decreasing range; used by `quad_dominant_random` and
  `orientation_defect`).
* `__cplusplus = 202002`, `_LIBCPP_VERSION = 180100`, `compiler`.

These are the ONLY RNG primitives the corpus generators use (`random_points`: U(0,box);
`assign_orientation_relaxation`: U(0, 2pi); `quad_dominant_random`: shuffle;
`orientation_defect`: uniform_int(0,1) + shuffle; `deployable_population`: normal).

## `reference_patterns/`

Copies of the stored pattern files the C++ tests and reference apps read; see
`reference_patterns/INDEX.md` for the per-file consumer table.

## Verification done at freeze time (python3, 2026-09-19)

* Every population file reloads; row counts as in the table; for every graph
  `len(mesh.vertices) == N`, `len(mesh.faces) == F`, `mesh.orientation == sigma_mc`,
  `len(sigma_*) == F`, all sigma values in {-1, +1}; F within `[min_faces, max_faces]` and
  N <= n_cap for every graph.
* `k1a_200.json` vertices, faces and `sigma_def` are bit-identical to the 200 archived
  `results/kill/k5/sigma/*.json` files (which store the graph alongside sigma_def) -- i.e.
  today's `make_graph` reproduces the graphs the 2026-09 kill runs used.
* `reference_cases_8.json` meshes equal `reference_patterns/cases/*/M.json` (8/8) and its
  `name/N/F` match `reference_patterns/reference_cases.json`.
* `scaling_42.json` N, F equal the archived `results/scaling/scaling.csv` for all 39 cells
  that CSV contains (3 cells were capped in the original run and have no row).
* `yield_fresh_100.json`: rerunning `method::orientation_defect(m0, sigma_mc, 20*F, 7000+id)`
  reproduced the archived `results/yield/fresh_sigma` sigma_def on 100/100 graphs.
* `k3a_500.json`: (id, N, F) equal the archived `results/kill/k3a/k3a.csv` for all 500 ids;
  `k2c_500.json`: the 462 ok ids and their (N, F) equal `results/kill/k2c/k2c.csv` (462 rows).
* `k2c_500.json` vs `k3a_500.json`: vertices differ on 212 ids (exactly those where the
  uncapped K3a graph has N > 1400) and are identical on the other 288.

## Other files in this directory

* `mt19937_vectors_local.json` was NOT written by `freeze_corpus` (it predates the freeze run;
  another teammate's local RNG check). `mt19937_vectors.json` is the one described above.

## Not frozen, and why

* `kill_e1.cpp`'s `sigma_def` (`defect_search_at`, seed `5000000 + id`, and its authored-
  tiling variant with seed `900001 + idx` and cap 3F) is a function in the app's anonymous
  namespace, not in the library, so it cannot be called from `freeze_corpus` without
  copying code. Recreate it from `kill_e1.cpp` lines 82-130 if E1 is ported.
* `kill_k5.cpp`'s own `defect_search` copy (the run that produced `results/kill/k5/sigma`)
  is likewise app-local; what is frozen is its ARCHIVED result (`sigma_def` in
  `k1a_200.json` / `native200.json`), which is what every downstream app consumes.
* K6's shape-space cache (`results/kill/k6/cache/*/shape_<id>.bin`, dense-SVD X0 and Phi)
  is solver output, not input; `kill_yield --mode k9c`, `kill_t1` and `kill_basin` read it.
  Regenerate with `assemble_system` + `solve_system` (BoundaryMode::Fixed) if needed.
* `kill_basin`'s design list comes from `results/kill/k9c/k9c.csv` rows (theta_exact <=
  1e-9 plus every 10th success); it is a subset of `native200.json` selected by an archived
  result CSV, so it is not a separate input.
* `kill_jitter`, `kill_k2a`, `kill_k7`, `kill_f23`, `kill_b3`, `kill_t1`'s extra inputs and
  the `kiri_*` apps use authored tilings via `generate(kind, {R}, rng)` or the reference
  cases, all reproducible from `reference_cases_8.json` / the tiling generators; they were
  not enumerated separately.

## Unit-test fixtures for `test_mesh_cut.jl` / `test_holes.jl` (`reference_patterns/test_fixtures_*.json`)

`code/tests/test_mesh_cut.cpp` and `code/tests/test_holes.cpp` load no files: every case
builds its meshes procedurally (`generate`, `delaunay_of_random_points`,
`voronoi_of_random_points`, `quad_dominant_random`, `tiling_squares`) with a case-local
`std::mt19937` seed, and draws sigma with `test::random_sigma` (`bernoulli(0.5)`),
`test::checkerboard_sigma` or `assign_orientation_relaxation`. The geometric hole checks
additionally run `assemble_system` / `solve_system` / `deploy` / `has_collision`. Until those
units are ported, the Julia tests read frozen replays:

| File | Producer | Contents |
|---|---|---|
| `reference_patterns/test_fixtures_mesh_cut.json` | `reference_patterns/freeze_fixtures.cpp` (built against `code/build/libkiri_core.a`, same toolchain as above, run 2026-09-19) | `remark_A1` (7 families: base mesh + 30 sigmas each; 40 Delaunay graphs with sigma as `orientation`), `A2A3_kagome` (kagome + 20 sigmas). |
| `reference_patterns/test_fixtures_holes.json` | same | `case1` (65 graphs: mesh with sigma, `forest`, `connected`, and where the C++ reached the geometric stage the deployed M'-vertex positions `Yd` and the C++ `holes_geometric` result `geo`), `rotating_squares`, `case3` (4 families x 8 sigmas), `case4` (107 `try_check` instances: mesh, embedding `X`, `Yd`/`geo` where deployed), plus the C++ run's own tallies under `counts`. |

Every block carries a `provenance` object (test case name, seed, generator call + args in
call order, sigma method, how `X`/`Yd` were produced). Meshes are in the standard mesh JSON
format (0-based); `geo` edge ids are 0-based C++ edge indices (edge numbering is the
`build_topology` traversal order, identical in the port). The freezer replays the C++ test
bodies verbatim, so `counts` equals the numbers the C++ tests print (65 graphs / 7 geometric
agreements / 37 split cycles / 53 disconnected / 0 mismatches; 60 checked / 1 of 30 / 47
skipped). The Julia tests recompute `holes_geometric(c, Yd)` and compare it against `geo`
and against the combinatorial preimages; they do not just echo the stored value.
`TODO(generators)` comments mark the load sites to swap for direct generator calls.

## Unit-test fixtures for `test_system.jl` / `test_kinematics.jl` / `test_collision.jl` / `test_rank_checks.jl`

Same situation as above: `code/tests/test_system.cpp`, `test_kinematics.cpp`, `test_collision.cpp`
and `test_rank_checks.cpp` build every input procedurally (tiling_* / periodic_* / torus_* /
generate / delaunay_of_random_points, sigma from `checkerboard_sigma`, `random_sigma` or
`assign_orientation_relaxation` on a case-local `std::mt19937`) and compute their expected
numbers inline. `reference_patterns/freeze_fixtures_2a.cpp` (built against
`code/build/libkiri_core.a`, Xcode clang++ `-O2 -std=c++20 -arch arm64`, run 2026-09-19)
replays the four test bodies verbatim and writes:

| File | Contents |
|---|---|
| `reference_patterns/test_fixtures_system.json` | one block per TEST_CASE: mesh (with sigma), the C++ `solve_system` report (`rank_full`, `rank_L`, `dim_null`, row counts, `projection_ok`, `sv_tol`, `X0`), the exhaustive-hexagon tallies (39/39/0), the (3,4,3,12) residual data (12 of 16, `per_hole`, `max_norm`), the 100 affine maps `[a00,a01,a10,a11,bx,by]` drawn from mt19937(19), the null-space offsets `[i][t] = [ux,uy]` drawn from mt19937(23) and the C++ `Phi`, the brute-force / relaxation orientation reports of the small hexagon patch, and the periodic_kagome(3,3) mesh. |
| `reference_patterns/test_fixtures_kinematics.json` | per family: mesh, sigma, the Eq. (6) embedding `X`, `max_mismatch` per theta; the 8 shuffled face orders (0-based) per kind of the BFS-order case; the FD case's per-trial sigma + 6 theta draws and its tally (20532 samples, worst 1.0875e-9); the Remark A.4 case (76 pairs). |
| `reference_patterns/test_fixtures_collision.json` | per kind: mesh, sigma, `X0`, `theta_max_geometric`, `min_beta`, `beta`, `collided`; the rotating-squares report; the deployment velocity `z`; the collision-sweep results (`theta_max_before/after`, `gamma_used`, full gamma ladder). |
| `reference_patterns/test_fixtures_rank_checks.json` | the 12 bounded and 3 torus cases (name, mode, mesh with sigma) and every field of the C++ `RowSumReport` (incl. `r`), `FactorizationReport` and `HingeGraphReport`, plus `solve_system(...).rank_L`. |

Every mesh object carries `faces_roundtrip_stable`: whether `mesh_from_json(mesh_to_json(m))`
returns the same face lists. It is `false` for the three torus meshes (their wrap-around faces
are geometrically degenerate, so the CCW normalisation flips them); the Julia loader
(`test/helpers.jl::fixture_mesh_raw`) therefore takes faces AS STORED and only builds topology.
The collision-sweep numbers are OPTIMIZER output (L-BFGS over a null-space basis whose columns
differ between Eigen and LAPACK) and are recorded for comparison, not asserted.
`TODO(generators)` comments mark the load sites to swap for direct generator calls.

## Export-writer fixtures for `test_export.jl` (`export_fixtures/`)

`code/tests/test_export.cpp` builds its meshes in place (a 3x3 `squares_patch`) except
one case that needs `tiling_3_4_3_12(disk(...))` + `assign_orientation_relaxation`. The
Julia port compares every writer against the C++ **bytes**, so the reference outputs were
frozen by `export_fixtures/freeze_export.cpp` (linked against `code/build/src/export/
libkiri_export.a` + `libkiri_core.a`, same toolchain as above, run 2026-09-19):

| File | What |
|---|---|
| `squares_3x3.json`, `squares_3x3_split.json` | the test's `squares_patch(3,3)` and its `make_split` variant (mesh JSON with `orientation`) |
| `t34312_disk2.5.json` | `tiling_3_4_3_12(disk(Vec2(0,0),2.5))` with `assign_orientation_relaxation(m, std::mt19937(12345))` sigma |
| `squares_felt_theta0.svg`, `squares_split_felt_theta0.svg`, `squares_felt_deployed.svg`, `squares_pla_theta0.svg`, `t34312_frac0.svg`, `t34312_frac1.svg` | `svg_string(build_layout(...), SvgOptions{})` at theta 0 / 0.8 theta_max, scale 10, felt_laser or pla_print |
| `squares_felt.stl`, `squares_pla.stl` | `write_stl_binary(build_solid(L))`, default header |
| `squares_felt.3mf`, `squares_pla.3mf` | `write_3mf(build_solid(L))`, default title |
| `index.json` | per-case layout summaries (bbox, hinge necks/setbacks/points, per-face outline sizes and areas, warnings), solid summaries (`check_manifold` report, `split_components` count), the C++ `theta_max` values, and a `provenance` block |

Comparison used by the Julia tests (`test_export.jl`, policy comment at its top): SVG text
byte-for-byte, always. STL and 3MF byte-for-byte first (the zip writer is STORE-only with a
fixed 1980-01-01 stamp, so the archive is deterministic); when the bytes differ the files
must still carry the identical vertex set / indexed vertex multiset, the same triangle and
object counts and the same manifold report, and on a native arm64 Julia the byte identity
is asserted outright.

Why the fallback exists: the reference ran natively on arm64 and its trig goes through
Apple's `__sincos_stret`; the Julia on this machine is x86_64 under Rosetta, where neither
Base nor the x86 libm reproduces those values to the last ulp. Vertex coordinates still
agree to 9+ digits, but two things move: the ~1e-18 residual of a float32 wall normal
(pad/hole arcs), and which of several exactly-degenerate needles ear clipping removes
first in a stalled face (same vertices, a few differently split triangles). Measured on
2026-09-19 (x86_64 Julia 1.12.1): 23/23 migrated SVGs byte-identical; of the 23 3D files
14 byte-identical, 9 structurally identical only (hero open_0.9tm 3MF/STL, the pin-pad
STL, and float32-normal-only diffs of 6-16 bytes in five STLs).

The reference build also contracts `x*y - z*w` into `fma(x, y, -(z*w))` (Apple clang,
arm64). The ear-clipping tie-breaks and the float32 STL normals depend on those sub-ulp
residuals, so `solid.jl` reproduces the contraction explicitly (`fms`, `cross_cpp`);
without it four triangles of the hero solid differ even at theta = 0.

The migrated exports under `export/{hero,hero2,samples}` (and the `sequence/` frames) are
also regenerated by `test_export.jl` from their input graphs (`hero*_graph.json`,
`reference_patterns/cases/<name>/M.json`) at the theta recorded in each report JSON and
compared byte for byte with the C++ files.

## Inputs of `derivation_tests.jl` (`reference_patterns/derivation_inputs{,_l2}.json`)

`code/tests/derivation_tests.cpp` builds every population procedurally.
`reference_patterns/freeze_derivation_inputs.cpp` (same toolchain as above, run 2026-09-19)
replays them: `derivation_inputs.json` holds `corpus` (the 16 cases of `corpus()`: mesh with
sigma from one `mt19937(20260903)` shared over the relaxation calls, plus the C++
`solve_system` X0 and null-space basis Phi, and k/H/n_hinge/n_split), `hloc` (6 square patches,
relaxation seed 5), `tori` (3 torus patches, relaxation seed 9) and `l1_corpus` (7 non-periodic
reference cases + `make_graph(id, 18, 46, 220)` until 60 usable, with X0/Phi);
`derivation_inputs_l2.json` holds the L2 periodic patterns `k7` (21 `make_tiling_pattern` + 12
Voronoi tori), `fresh110` and `ce2000` (cell mesh with sigma from `quotient_sigma`, T row-major).
The Julia test rebuilds what it can with the Julia generators and checks it against these
(on the arm64 Julia every mesh, sigma and Voronoi torus is bit-identical; `make_tiling_pattern`
is not for squares/kagome), then runs on the frozen copies and samples the shape space with
the C++ Phi, so every printed C++ tally is reproduced exactly (see the test file header).
