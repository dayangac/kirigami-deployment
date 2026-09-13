# Kirigami.jl -- package tour

One module, `Kirigami`, that `include`s one file per C++ translation unit of the reference
implementation (`src/Kirigami.jl` lists the order). Nothing is exported: call
`Kirigami.f(...)`. C++ function and field names are kept so that every equation, theorem
and experiment reference in `derivations/`, `results/` and the papers still resolves; the
deviations are listed at the end. Conventions are in `../PORTING.md`; tests run with
`julia --project=Kirigami -e 'using Pkg; Pkg.test()'` (187 test sets / 176,692 assertions (23 marked broken), 0 failures).

Types: `Vec2 = SVector{2,Float64}`; dense matrices are `Matrix{Float64}`; shape-space
coefficients `t` are `Vector{Float64}` of length `2 * dim_null`; the RNG is the bit-exact
`MT19937`. Indices are 1-based everywhere in memory; JSON files keep the C++ 0-based format
and convert at the I/O boundary.

## `src/core/` -- the 2026 pipeline

**`mesh.jl`** -- the uncut embedded planar graph `Mesh` (`X`, `faces`, `sigma`, optional
`periodic`), its half-edge topology (`build_topology!`), `face_signed_area`,
`corner_angle`, and the sigma-induced half-edge direction `half_edge_direction` of 2026
Sec. 3 (sigma = +1 clockwise, -1 counter-clockwise). JSON I/O: `load_mesh_json`,
`save_mesh_json`, `mesh_from_json`, `mesh_to_json` (the README graph contract: vertices,
faces, orientation, periodic).

**`cut.jl`** -- `make_cut(m) -> CutStructure`: edge classification into
`Hinge`/`Split`/`Border` and the kirigami structure M' = (X', F') (2026 Sec. 3); a hinge
cut keeps the hinge at the SOURCE vertex and duplicates the target. `check_remark_A1`
(hinge in-degree = out-degree at interior vertices), `count_components`, `prime_vertex`.

**`holes.jl`** -- hole preimages: `holes_seed_growing` (2026 Sec. 4.1, Algorithm 1),
`holes_partition` (the independent split-forest partition formulation), `same_hole_sets`,
`holes_partition_edges`, `split_subgraph_is_forest` (Remark A.4), and the geometric
cross-check `holes_geometric` / `holes_geometric_cycles` on the deployed M' (Definition 4.1).

**`tutte_auxetic.jl`** -- `hole_residuals` (the deployability residual, Eq. (2)),
`assemble_system(c, hs, X_ini, mode::BoundaryMode)` (the Tutte auxetic linear system,
Eqs. (3)-(5), with `Fixed`, `Periodic` or `None` boundary rows) and `solve_system` (the
Eq. (6) projection `X0` and the null-space basis `Phi`, dense SVD with the C++ rank rule);
`rank_only_sparse` for large graphs; `matrix_to_points` / `points_to_matrix`.

**`kinematics.jl`** -- `deploy(c, X, theta) -> Deployment`: forward kinematics of M' (2026
Sec. 4.5), every face rotated by `-sigma(f) theta/2` and placed by BFS over the hinge
adjacency, with `max_mismatch` over non-tree hinges; `rigid_align_residual`.

**`collision.jl`** -- exact polygon predicates (`polygons_overlap`, `has_collision`), the
grid-plus-bisection `theta_max` referee and `hinge_beta` (2026 Sec. 5.2), and the
collision-aware optimisation of Eqs. (7)-(9): `optimize_collision`,
`optimize_collision_sweep` (a gamma ladder, never worse than the input).

**`orientation.jl`** -- `assign_orientation_relaxation(m, rng, restarts, iters,
n_diameters)`: the Eq. (1) max-cut relaxation (2026 Sec. 4.2), bit-exact for a given
`MT19937`; `brute_force_orientation`; `dual_graph`; `describe_orientation`.

**`generators.jl`** -- authored tilings clipped to a `disk`/`rect` (`tiling_squares`,
`tiling_triangles`, `tiling_hexagons`, `tiling_kagome`, `tiling_3_4_3_12`,
`tiling_snub_square`, `tiling_truncated_square`, the `periodic_*` and `torus_*` patches)
and the random families (`delaunay_of_random_points`, `voronoi_of_random_points`,
`quad_dominant_random`), all behind `generate(kind, params, rng)`. The random meshes are
bit-exact with the C++ for a given `MT19937`; the `libm_*` shims call the system libm.

**`mt19937.jl`** -- `MT19937`, `next_u32`, `uniform_real`, `uniform_int`, `shuffle!`,
`NormalDist`/`normal`: `std::mt19937` and the libc++ distributions the corpus was drawn
with, bit-exact (`data/corpus/mt19937_vectors.json`).

**`import_soup.jl`** -- `weld_segments`, `read_svg_segments`, `import_svg_soup`: the welding
importer that turns a polygon soup (the 2025 paper's fabrication SVGs) into a `Mesh`.

**`optimize.jl`** -- `lbfgs_minimize(fg, x0, LbfgsOptions())`: the self-contained L-BFGS
with backtracking line search every optimiser in `method/` uses (no external solver).

**`rank_checks.jl`** -- the three measurements on the hole-constraint matrix L behind
theorem T7 / the 2026 Sec. 4.4 erratum: `check_row_sum` (1ᵀL and the notch-restricted
degree identity), `check_factorization` (L = R·D and the out-harmonic left null space),
`check_hinge_graph` (H = |E_hinge| - |F| + c(Γ)).

**`kill_common.jl`** -- the graph populations: `make_graph(id, min_faces, max_faces,
n_cap)` (the kill experiments' population, id-deterministic through `MT19937`),
`checkerboard`, `reference_cases()` (the eight Phase-2 cases).

## `src/method/` -- this project's theory and algorithms

**`deploy_basis.jl`** -- `deploy_basis(c, X) -> DeployBasis`: the trig-linear form
`Y(θ) = cos(θ/2) C + sin(θ/2) S` with `C`, `S` linear in `X` (T1, `derivations/core.md`),
and the harmonic calculus on top of it: `Harmonic` (`p + q cos θ + r sin θ`),
`orient_harmonic` / `dot_harmonic` (T3: the orientation determinant and the two projection
inequalities of a vertex-edge pair), `harmonic_roots`, `classify_harmonic` and
`harmonic_roots_deflated` (the tau = 0 deflation of the systematic degeneracies).

**`contact.jl`** -- the exact first-contact calculus. `swept_discs` and `candidate_pairs`
(the sound moving broad phase, T4.5a), `contact_angles` (the complete contact set C(X),
Lemma T4.2 / Corollary T4.2'), `exact_theta_max_overlap` (Θ_max by T4.2'': first interval
whose midpoint has an interior overlap; `zero_range` for penetration at 0⁺), and
`validity_certificate` (POS ∧ NOOVERLAP(ε/2) ∧ NOROOT(ε), Proposition T5.2b' with the
interval inequalities of T4.1b; `valid(cert)`), `penetrates_immediately`.

**`design.jl`** -- the deliverable API. `characterize(mesh, sigma, X)` (exact range,
contact list, certificate, largest certified `eps_max`, binding type);
`design_constrained` (K9 variant (b): convexity + split-outward feasibility by proximity
to `X_ini`), `design_baseline` (the Eq. (6) projection alone, t = 0),
`design_range_max` (K9c: arms k9/k9b, stage A `range_embed` from each start, stage B
`maximize_margin_range`, best-of-three by exact Θ_max with `provenance`);
`orientation_maxcut` and `orientation_defect` (K5's defect-minimising sigma). Defaults are
the K9 / K9c runs' parameters.

**`zero_plus.jl`** -- the 0⁺ separation calculus of a split cut: `zero_plus_q` (the split
signs `q_e`), `zero_plus_corner_margin` (the corner incidences `mu_j`), `zero_plus_form`
(their quadratic form in `t`), and `zero_plus_repair` (K6's split-only repair).

**`convex_embed.jl`** -- `corner_crosses`, `convex_embed(c, X0, Phi, X_ini, med, opt)`:
the convexity-constrained (and split-constrained) point of the shape space found by
penalty continuation plus log-barrier proximity stages (experiments K9, K9b).

**`range_embed.jl`** -- the paper's headline algorithm (K9c). `zero_plus_margin` (the exact
margin `m(X) = min(min q_e, min mu_j) / med²`), `range_embed` (stage A: log-sum-exp softmin
of the margin inside the barriers, mode refresh per stage) and `maximize_margin_range`
(stage B: push the exact Θ_max under a shrinking trust region with exact rejection).

**`range_opt.jl`** -- `maximize_range` (K2b): softmin of the closed-form first-contact
angles over the null space with the T6 implicit-differentiation gradients and an
active-set refresh; `RangeObjective` exposes objective and gradient for finite-differencing.

**`mobility.jl`** -- `build_hinge_graph` (Γ = (F, E_hinge) with its cycle space),
`build_rigidity`, `two_core`, `mobility_at`: the angular-velocity matrix and the mobility
identity `dim ker A = |F \ core2| + dim ker A|core2` of experiment K3a.

**`periodic_jacobian.jl`** -- the torus quotient of a periodic pattern (`detect_lattice`,
`fundamental_domain`, `build_quotient`, `quotient_system`) and `periodic_jacobian`: the
deployment Jacobian of the fundamental parallelogram, 2026 Sec. 5.1 Eqs. (11)-(13), whose
exact conformality is contribution C5 (K7); `fk_period_matrix` cross-checks it by forward
kinematics.

**`expansive_cone.jl`** -- `cone_system`, `cone_lp`, `expansive_cone`: the polyhedral cone
of first-order expansive flexes of a cut structure and the LP (Frank-Wolfe on the min-norm
point, `farkas_residual` dual certificate) that decides it, experiment K8a; `flex_basis`,
`sigma_flex`, `euler_step`.

**`budget.jl`** -- `face_potential`, `budget_terms`, `periodic_cell_edges`: the
expansion-budget identity of `ideas/round2_theorist_b.md` (B1/B3), the first-order
hole-opening rate as a border functional split into hinge and split parts.

## `src/export/` -- fabrication

**`material.jl`** -- `MaterialProfile` (`profile_felt_laser`, `profile_paper_laser`,
`profile_pla_print`, `profile_by_name`; `HingeType` `LivingHingeNeck` / `PinPad`), the only
place fabrication constants live. **`layout.jl`** -- `build_layout(m, c, X, theta, scale,
profile) -> Layout`: the deployed, millimetre-scaled, inset face outlines with a neck at
every hinge; consumed by every writer. **`svg.jl`** -- `write_svg` / `svg_string`
(cut / score / engrave groups). **`solid.jl`** -- `build_solid(L) -> (TriMesh, warnings)`
(face extrusion, one closed manifold per living-hinge sheet), `check_manifold`, `summary`.
**`stl.jl`** -- `write_stl_binary`, `read_stl_binary`. **`threemf.jl`** -- `write_3mf`,
`split_components`. **`zip.jl`**, **`xml.jl`** -- the STORE-only zip writer and the
well-formedness checker the 3MF/SVG tests use.

## Reproducibility of the numbers

`characterize`, the certificate, the generators, the RNG and the linear algebra reproduce
the C++ to the tolerances the tests state (bit-exact for RNG streams, generated meshes and
the closed-form contact calculus at a given point). The optimiser outputs
(`convex_embed`, `zero_plus_repair`, `range_embed`, `maximize_margin_range`) are
reproduced to path level: same objective, gradient and start points to 1e-14, converged
solves agree to solver tolerance, and the fixed-budget `range_embed` diverges by rounding
amplification over its 6 x 120 iterations, so the K9 / K9c CSV rows are reproduced at CSV
precision by the converged arms and NOT to 1e-9 through stage A (those C++ locks are kept
as `@test_broken` in `test/test_design.jl`; call-by-call fixtures in
`../data/corpus/method_fixtures/`).

## Deviations from the C++ API

| C++ | Julia | where |
|---|---|---|
| mutating members `m.build_topology()`, `m.normalize_face_ccw()`, `std::shuffle`, `hs.finalize()`, `g.build_forest()` | Julia `!` suffix: `build_topology!(m)`, `normalize_face_ccw!(m)`, `shuffle!(v, rng)`, `finalize!(hs, c)`, `build_forest!(g)` | `mesh.jl`, `mt19937.jl`, `holes.jl`, `mobility.jl` |
| out-pointers (`std::vector<int>* bad`, `std::string* diff`, `double* min_q`, ...) | extra return values, as tuples: `check_remark_A1(c) -> (ok, bad)`, `same_hole_sets -> (same, msg)`, `holes_partition_edges -> (ok, multiply_covered, uncovered)`, `split_subgraph_is_forest -> (ok, cycle_vertices)`, `zero_plus_margin -> (margin, min_q, min_mu)`, `build_solid -> (TriMesh, warnings)`, `subgraph -> (HingeGraph, sub_pins)`, `matrix_rank -> (rank, used_sparse)`, `face_potential -> (u, worst_closure)`, `periodic_cell_edges -> (keep, n_preimage)`, `fk_period_matrix -> (P, spread)`, `farkas_residual -> (resid, lambda_min, sum_err)`, `harmonic_fit -> (h, max_res, scale)`, `read_svg_segments -> (segments, n_path_ignored)`, `xml_well_formed -> (ok, err)`, `*_objective(...) -> (F, grad)` | throughout |
| `enum class EdgeType { Border, Hinge, Split }`, `BoundaryMode::{None, Fixed, Periodic}`, `HingeType`, `HarmonicClass`, `ConeRowKind` | `@enum` with the same member names, unscoped: `Kirigami.Hinge`, `Kirigami.Fixed`, ... | `cut.jl`, `tutte_auxetic.jl`, `material.jl`, `deploy_basis.jl`, `expansive_cone.jl` |
| 0-based indices; "-1 = none" sentinels (`parent`, `pedge`, `ContactWitness.pv`, `bad_pv`) | 1-based; "0 = none" for index-valued sentinels. Exceptions kept as the C++ wrote them: `OverlapRangeReport.i_star` and `Characterization.i_star` (C++ index into `candidates`, 0 => zero range, -1 => none), `DefectOrientationResult.restart_used`, `best_start` (0 = t = 0) | `contact.jl`, `mobility.jl`, `design.jl` |
| member functions `cert.valid()`, `B.c(i)`, `B.s(i)`, `B.n_prime()`, `r.deployable(tol)`, `L.width()`, `mr.summary()` | free functions `valid(cert)`, `basis_c(B, i)`, `basis_s(B, i)`, `n_prime(B)`, `deployable(r, tol)`, `width(L)`, `summary(mr)` | `contact.jl`, `deploy_basis.jl`, `tutte_auxetic.jl`, `layout.jl`, `solid.jl` |
| options structs with default member initialisers; `Eigen::VectorXd t_init` empty = t = 0 | `Base.@kwdef mutable struct` (`ConvexEmbedOptions(delta_rel = 1e-2)`); `t_init::Vector{Float64} = Float64[]`, empty = t = 0 | `method/*.jl` |
| `std::optional<T>`, nullable pointers | `Union{T,Nothing}` | throughout |
| `std::invalid_argument` / `std::runtime_error` | `ArgumentError` / `error(...)` | `design.jl` and elsewhere |
| `Eigen::VectorXd t`, `MatrixXd X0` (N x 2) | `Vector{Float64}`; `SolveReport.X0` stays `Matrix{Float64}` N x 2, converted with `matrix_to_points` | `tutte_auxetic.jl` |
| `CutStructure` holds `const Mesh*` | `CutStructure.mesh::Mesh` (a reference to the same object; the mesh must not be mutated after `make_cut`) | `cut.jl` |
| `CollisionSweepResult : CollisionOptResult` (inheritance) | one flat struct with all fields | `collision.jl` |
| app-side `kill::` helpers (`Shape` cache, `deployable_population`, `sci`/`fx`, `Timer`) | `apps/kill_common.jl` (an included file, not part of the module), reading `data/corpus/` by default | `apps/` |
