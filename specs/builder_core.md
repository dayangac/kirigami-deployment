# Builder-Core: C++ reimplementation of the Segall 2026 pipeline

Output: code/ (CMake project), all tests passing, plus code/README.md documenting the API and the CLI tools.
Read first: notes/paper_2026.md (full), notes/paper_2025.md (Sec 3–4), and specs/common_preamble.md. If notes are missing, read papers/paper_2026_tutte_flow.txt.

## Scope (reference pipeline; do NOT invent new methods)
Implement, in this order, committing progress to files continuously:

### 0. Project skeleton
- code/CMakeLists.txt (C++20, `-O2`, `CMAKE_OSX_ARCHITECTURES=arm64`, warnings on). Targets: `kiri_core` (static lib), `kiri_tests` (doctest), CLI apps under code/apps/.
- Headers under code/src/core/, one concern per file. Namespace `kiri`.
- Dependencies: Eigen (dense+sparse), nlohmann/json, doctest — all header-only in /opt/homebrew/include.
- Build command must be documented and verified: `cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j && ./code/build/kiri_tests`.

### 1. Graph I/O and data model (`mesh.hpp`)
- JSON format (this is the contract for all later agents — implement exactly):
  ```
  { "vertices": [[x,y], ...],
    "faces":    [[i0,i1,...], ...],        // each face CCW in the given geometry, simple polygon
    "orientation": [+1|-1, ...],            // optional, one per face (σ). +1 = clockwise per paper, -1 = CCW
    "periodic": { "th":[tx,ty], "tv":[tx,ty], "pairs_h":[[p,p'],...], "pairs_v":[[q,q'],...] }  // optional
  }
  ```
- Build a half-edge structure: for each face, half-edges in the face's cyclic order. Interior edge = shared by exactly 2 faces; boundary edge = 1 face. Reject non-manifold input (edge shared by >2 faces) with a clear error. Detect interior vs boundary vertices.
- The paper's orientation σ induces a direction on each half-edge: faces with σ=+1 traverse their boundary clockwise, σ=−1 counter-clockwise (2026 Sec 3). Implement `half_edge_direction(face, edge)` accordingly and CHECK against 2026 Fig. 4: for adjacent faces with opposite σ, the two half-edges of the shared edge have the SAME direction (v_i→v_j); with equal σ they are opposite.

### 2. Edge classification and cutting (`cut.hpp`)
- Classify interior edges into E_hinge (σ differs; directed v_i→v_j with source v_i = hinge vertex, target v_j duplicated) and E_split (σ equal; both endpoints duplicated). Boundary edges untouched.
- Build the cut structure M' = (X', F'): for each face, new vertex indices. Rule: an original vertex v is shared between two faces adjacent around v iff their shared edge is a hinge cut whose SOURCE is v. Implement by walking around each vertex and grouping its incident face-corners into equivalence classes joined by hinge cuts with source v. Each class becomes one vertex of M'. Store the map corner → M'-vertex and M'-vertex → original vertex.
- Verify Remark A.1 numerically on all test graphs: at every interior vertex, #hinge edges directed in == #hinge edges directed out. Fail loudly otherwise (this is a correctness check on the direction convention).

### 3. Hole preimage detection (Algorithm 1, 2026 Sec 4.1) (`holes.hpp`)
- Implement the seed-growing algorithm exactly as in the paper (rules: split–split always merge; hinge–hinge merge iff both point INTO the shared vertex; hinge–split merge iff the hinge points INTO the shared vertex; border edges excluded). "Neighboring edges" of an edge = edges sharing an endpoint with it.
- IMPORTANT clarification from the orchestrator (the paper's pseudocode line 10 is incomplete and the paper never states whether preimages partition the edges): implement the three PROSE rules, where for a hinge edge only its TARGET endpoint is the shared vertex that matters. Equivalent clean formulation to implement as a SECOND independent routine and cross-check against the seed-growing one: let S be the subgraph of M formed by split-cut edges (expected to be a forest; assert and report if a cycle appears). For every connected component K of S, and for every isolated interior vertex v not touched by any split edge (treat as a trivial component {v}), the hole preimage is C_K = E(K) ∪ { hinge edges whose target vertex ∈ V(K) }. Components consisting only of boundary vertices are skipped. The preimages must partition E_hinge ∪ E_split — assert it.
- Enumerate ALL hole preimages by seeding from every not-yet-assigned hinge edge. Then handle split-cut edges not reached from any hinge seed (they should be reached; assert and report if not).
- Cross-check: independently compute holes GEOMETRICALLY on M' deployed at a small θ (after step 5): each hole is a boundary cycle of M' made only of duplicated interior edges; map back to original edges and compare with the combinatorial preimages. Both must agree on every test graph. Report the number of holes H.

### 4. Deployability test and the linear system (Eq. 2–6) (`tutte_auxetic.hpp`)
- `residual(C)` = Σ_{hinge e∈C} (x_target − x_source) for each hole preimage; deployable iff all ≈ 0 (Prop 4.1).
- Assemble L (H×N, sparse): row per hole preimage; +1 at target, −1 at source for each hinge edge in it (entries accumulate).
- Boundary: (a) fixed — rows selecting boundary vertices, rhs = given positions; (b) periodic — rows x_p − x_p' = t_h etc. per Eq. 3b–c, plus pin ONE vertex to kill translation.
- Solve: dense for N ≤ ~5000 unknowns (Eigen ColPivHouseholderQR / BDCSVD), sparse SparseQR otherwise. Provide: rank, null-space basis {φ_i} (columns, N×k) via SVD of the constraint matrix, particular solution X0 = argmin ‖X − X_ini‖ s.t. constraints (Eq. 6; implement as projection onto the affine solution set using the SVD).
- Report `dim_null = N_interior − rank` and compare with the paper's prose claim "one hole per interior vertex ⇒ full rank when E_split = ∅" and "#holes ≤ #interior vertices". Log both numbers for every test graph into a CSV (see step 8). Do NOT assume the claim; measure.

### 5. Forward kinematics (`kinematics.hpp`)
- Given M', σ, and θ, compute Y_θ: pick a seed face (fixed), BFS over the hinge graph (faces adjacent via shared M'-vertices), rotate each newly reached face rigidly about the shared hinge vertex by ±θ relative to its parent, sign determined so that the two duplicated copies of the hinge edge open by exactly θ (check the sign convention against 2026 Fig. 5(c)/Fig. 8: all hinges open by the same θ; in a uniformly deployable structure the result must be independent of BFS order — verify this: recompute with a different seed / order and compare, tolerance 1e-9 after aligning by a rigid motion).
- Consistency check: for a deployable embedding, every M'-vertex reached via two different BFS paths lands at the same position (this is the geometric content of Eq. 2). For a NON-deployable embedding, report the max mismatch; it must be > 0 (use 2026 Fig. 21's (3,4,3,12) tiling as the test).
- Also implement `derivative_wrt_theta` (∂Y_θ/∂θ) by differentiating the rotation chain (needed for Eq. 7–9) and validate against finite differences (relative error < 1e-6 on ≥ 1000 random samples).

### 6. Collision handling (Eq. 7–9) and deployment range
- Polygon–polygon overlap test (exact for convex, SAT or clipping; for non-convex use triangulation or a robust segment-intersection + point-in-polygon test). `theta_max(Y0)`: bisection on θ ∈ (0, π] for the first overlap between any two faces (use a uniform grid + bbox culling). Also compute the analytic per-hinge bound β_i = 2π − α_lik of 2026 Sec 5.2 and report min_i β_i alongside; note when they differ (they differ when collisions happen at split-cut edges earlier, cf. 2026 Fig. 14–15).
- Collision energy Eq. (9): z_i = row of ∂g/∂α at α=0 (== ∂Y_θ/∂θ at θ=0), barrier B(s) = soft penalty for s<0 (e.g. log-barrier or smooth hinge), summed over split-cut edges, plus γ‖Y−Y0‖², minimized over Y ∈ X (parametrize Y = X0 + Φ t and optimize t with L-BFGS or Newton; write your own small L-BFGS or gradient descent with line search — no external libs). Validate: on a case where θ_max is small before optimization, θ_max after optimization must increase; print both.

### 7. Orientation assignment (Eq. 1, 2026 Sec 4.2)
- Dual graph; relaxation: minimize Σ_{(i,j)∈E_d} x_i·x_j s.t. ‖x_i‖=1 in R² (projected gradient descent from random init, several restarts), then split by the best of a set of candidate diameters (maximize cut). Also provide `brute_force_orientation` for ≤ 20 faces (enumerate all 2^F, filter: connected after cutting) for tests. Also allow reading σ from JSON.
- Report #split edges, connectivity of M', number of holes.

### 8. Test-graph generators (`generators.hpp`) and CLI apps
- Generators: regular tilings clipped to a disk or rectangle: triangles, squares, hexagons, kagome (3.6.3.6), (3,4,3,12) 2-uniform tiling (2026 Fig. 21), snub square (3.3.4.3.4), truncated square (4.8.8); random planar graphs: Delaunay triangulation of n random points, Voronoi diagram of n random points clipped to a box (polygonal faces), and random quad-dominant meshes via edge-collapse of Delaunay. Optional periodic versions of the regular tilings (n×m cells with pairs_h/pairs_v filled in).
- Apps: `kiri_analyze <graph.json> [--orient auto|json|brute] --out <dir>`: writes cut structure, hole preimages, rank/null dim, residuals, X0, θ_max, and a CSV summary; `kiri_deploy <graph.json> --theta t --out y.json`; `kiri_gen <kind> <params> --out g.json`; `kiri_sweep --n 100 --out sweep.csv` running the whole pipeline on random graphs (100–5000 faces) and logging: N, F, #interior, #hinge, #split, H, rank, dim_null, deployable?, θ_max_geometric, min β_i, timings.
- Plotting: code/scripts/plot_embedding.py (reads the JSON dumps; draws M, M' at several θ with holes shaded) — run with `arch -arm64 /usr/local/bin/python3`.

### 9. Tests (doctest, code/tests/) — all must pass; write failing-first where possible
- Direction convention (Fig. 4 check), Remark A.1 in/out balance on ≥ 200 random (graph, σ) pairs.
- Rotating squares (checkerboard σ on square grid): deployable for the regular grid; unique solution (dim_null = 0); FK produces the classic rotating-squares pattern (hole = square of side ∝ sin(θ/2)…—just check all hinges open by θ and no mismatch).
- Equilateral triangle tiling with alternating σ (kagome-type cut of 2025 Fig 2 row 1): deployable, dim_null 0.
- Hexagonal tiling: with the 4 σ assignments of 2026 Fig 6 (reconstruct them from notes; if the exact σ is unrecoverable, use: (a) all-alternating impossible since hex tiling dual is triangular ⇒ every σ has split edges; enumerate σ on a small patch by brute force, keep connected ones) — check that dim_null == #interior − H on every case and report any counterexample to the paper's prose claim.
- (3,4,3,12) tiling: NOT deployable as Euclidean tiling (residual of the preimage C1 ≠ 0, 2026 Fig 21b); after projection Eq. 6 it IS deployable; FK mismatch < 1e-9.
- Affine invariance (Remark 4.1): apply 100 random affine maps to a deployable embedding; residuals stay 0 (to 1e-12 relative).
- Remark A.4 numerically: on deployable embeddings with split edges, duplicated copies of every split-cut edge remain parallel for θ ∈ {0.1,…,θ_max} (cross product < 1e-9).
- Hole preimage: combinatorial vs geometric agreement on ≥ 50 random graphs.
- FK derivative vs finite differences.
- Null-space: for each φ_i, X0 + φ_i t^T remains deployable for 100 random t.
- Collision: θ_max via geometry ≤ min β_i always; equality on rotating squares.

### 10. Phase-2 gate artifacts (write to results/core_validation/)
- `reference_cases.md`: for each reference case above, the numbers (rank, dim_null, residual norms, θ_max) and a PNG of M and M' at θ = 30°, 60°, θ_max.
- `sweep.csv` from `kiri_sweep --n 100` on random Voronoi/Delaunay graphs 100–5000 faces + a histogram PNG of dim_null vs (#interior − H).
- `rank_claim.md`: table stating for how many of the sweep graphs dim_null == #interior − H held, and listing every violation with its JSON saved under results/core_validation/violations/.

## Acceptance criteria
1. `cmake --build` succeeds with zero warnings-as-errors disabled but no errors; `kiri_tests` reports all passed.
2. All step-10 artifacts exist. Any failing reference case is documented, not hidden.
3. code/README.md explains the JSON contract, the direction/σ convention with a small ASCII diagram, and every CLI.
Report back: build/test output tail, the numbers in reference_cases.md, and every place where your measurements disagree with the paper's text.
