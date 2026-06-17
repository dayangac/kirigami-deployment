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
| `scaling_42.json` | 42 cells | `kill_scaling::build_graph(kind, sites, seed, true)`: `rng = mt19937(1000003*seed + 20260908 + 7919*sites)`, `generate(kind, {sites, 40.0}, rng)`, then `assign_orientation_relaxation(m, rng, 4, 300, 90)`; grid of `code/scripts/run_scaling.sh`: targets {50,100,200,500,1000,2000,5000}, delaunay sites = target/2, voronoi sites = target, seeds 0,1,2 | `sigma_mc` | `kill_scaling` | Fields `kind, target_faces, sites, seed, rng_seed, N, F, E, mesh, sigma_mc`. N, F agree with all 39 cells present in the archived `results/scaling/scaling.csv`. |
| `reference_cases_8.json` | 8 | `kill::reference_cases()`: the Phase-2 tilings with ONE shared `mt19937(20260903)` for the four relaxation calls (`assign_orientation_relaxation(g, rng, 8, 500, 180)`, in file order), checkerboard for the rest | `sigma` (per case) | `kiri_reference`, `test_reference_cases.cpp`, `test_design.cpp`, `derivation_tests.cpp` | Fields `name, periodic, N, F, mesh, sigma`. `periodic_squares_4x4` carries `mesh.periodic`. Verified equal to `reference_patterns/cases/<name>/M.json` for all 8. |
| `deployable_population.json` | 187 configs (35 bases + 152 shape-space samples) | `kill::deployable_population(8, 0.2)`: 7 tiling families x radii {2.0,2.5,3.0,3.5,4.0}, sigma checkerboard (squares, triangles, kagome) or `assign_orientation_relaxation(m, mt19937(777 + 13*ri + 101), 8, 500, 180)`, X = Eq. (6) projection (dense SVD), samples from `mt19937(20260903)` normal draws filtered for no inverted face and no collision | `mesh.orientation` | `kill_k1a`, `kill_k2b`, `kill_e1` part (A) | Fields `name, family, dim_null, sample (-1 = base), mesh, X`. NOTE: unlike the others this contains solver OUTPUT (`X`), so a port should match it to ~1e-9, not bit-exactly. |

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
