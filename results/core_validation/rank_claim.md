# Rank claim on the sweep (Builder-Core)

Source: `sweep.csv`, produced by `kiri_sweep --n 50 --seed 20260903`.
Random planar graphs (Voronoi cells of random points, Delaunay triangulations,
quad-dominant edge-collapsed Delaunay), face orientations from the paper's own Eq. (1)
relaxation, fixed-boundary system of Eq. (4).

Graphs analysed: **50**, faces from **115** to **2734**.

## The paper's Sec. 4.4 claims, measured

| claim (2026 Sec. 4.4, prose, unproved in the paper) | holds on |
|---|---|
| `dim null(Eq. 4) == #interior vertices - H` | **50 / 50** |
| `rank(L) == H` (the rows of L are independent) | **50 / 50** |
| `H <= #interior vertices` | **50 / 50** |
| stronger relation found here: `dim null == |E_split|` | **50 / 50** |

## Structural checks

| check | holds on |
|---|---|
| Remark A.1 (hinge in-degree == out-degree at interior vertices) | 50 / 50 |
| Algorithm 1 (seed-growing) == split-forest partition formulation | 50 / 50 |
| hole preimages partition E_hinge union E_split | 50 / 50 |
| split-cut subgraph is a forest (Remark A.4) | 50 / 50 |

## Violations

None. `dim_null == #interior - H` held on every graph in the sweep.

## Why the geometric hole cross-check is empty on this sweep

The geometric check traces the bounded complement components of the deployed M', which
requires the Tutte auxetic embedding to be a valid (non-self-intersecting) straight-line
embedding. On these random graphs it never is:

| quantity | value |
|---|---|
| graphs where the dense null space / X0 was computed | 27 / 50 |
| of those, X0 free of face-face overlap at theta = 0 | **0 / 27** |
| mean fraction of faces whose orientation flips in X0 | 0.058 |
| largest vertex displacement of the Eq. (6) projection (box side 40) | 11.095 |
| geometric cross-checks actually run / agreed | 0 / 0 |

This is the paper's own open problem (2026 Sec. 6): the shape space X contains
embeddings with self-intersections and the paper offers no certificate for the valid
subset. The combinatorial-vs-geometric agreement is instead verified on 60 valid
instances by `Pkg.test()` (see `test/test_holes.jl`).

Worst dense solve time in the sweep: 2382.611 ms.

## Orchestrator checks (added)

Three measurements on the hole-constraint matrix `L` that nothing else in the
pipeline needed: the row sum `1^T L`, the factorization `L = R D`, and the
hinge-graph Euler count. `L` is taken verbatim from `assemble_system` (its first
`H` rows), so what is measured is the matrix the solver actually uses. Definitions:

- `D` is the `|E_hinge| x N` signed hinge incidence (+1 target, -1 source);
  `R` is the `H x |E_hinge|` 0/1 matrix saying which hole preimage owns each hinge edge.
- `Z` is the left null space of `L`; `g(v) = y_{K(v)}` with `K(v)` the preimage whose
  split-forest component contains `v`, and `g(v) = 0` when `K(v)` is a notch (no row of L).
- `Gamma = (nodes = faces, edges = hinge cuts)`; `c(Gamma)` is its component count.
- `H` counts ONLY all-interior preimages, as Builder-Core defines it; the
  boundary-touching preimages are counted separately as notches.

### Reference cases (read back from `cases/*/M.json`, boundary rows omitted)

| case | F | \|E_hinge\| | H | notches | c(G) | pred = \|E_hinge\|-\|F\|+c(G) | rank(L) | dim Z | 1^T L = 0 | naive row sum | restricted row sum | L = R D | rank = H - dim Z | Z out-harmonic | H = pred | H_all = pred |
|---|--:|--:|--:|--:|--:|--:|--:|--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| `hexagons_auto` | 10 | 13 | 4 | 5 | 1 | 4 | 4 | 0 | **no** | **no** | yes | yes | yes | yes | yes | **no** |
| `kagome_3636` | 97 | 174 | 78 | 18 | 1 | 78 | 78 | 0 | **no** | **no** | yes | yes | yes | yes | yes | **no** |
| `periodic_squares_4x4` | 16 | 24 | 9 | 6 | 1 | 9 | 9 | 0 | **no** | **no** | yes | yes | yes | yes | yes | **no** |
| `rotating_squares` | 25 | 40 | 16 | 8 | 1 | 16 | 16 | 0 | **no** | **no** | yes | yes | yes | yes | yes | **no** |
| `snub_square_33434` | 53 | 64 | 12 | 14 | 1 | 12 | 12 | 0 | **no** | **no** | yes | yes | yes | yes | yes | **no** |
| `tiling_3_4_3_12` | 25 | 40 | 16 | 8 | 1 | 16 | 16 | 0 | **no** | **no** | yes | yes | yes | yes | yes | **no** |
| `triangles_alternating` | 89 | 121 | 33 | 19 | 1 | 33 | 33 | 0 | **no** | **no** | yes | yes | yes | yes | yes | **no** |
| `truncated_square_488` | 21 | 32 | 12 | 8 | 1 | 12 | 12 | 0 | **no** | **no** | yes | yes | yes | yes | yes | **no** |

### Boundary-free (torus) patches -- where `1^T L = 0` and `dim Z > 0`

The generators' "periodic" patches are finite patches with a real boundary (Eqs. 3b-3c
are constraints, not a topological identification), so they do NOT test the
boundary-free case. `torus_squares` / `torus_triangles` do: every vertex is interior.
Their wrap-around faces are geometrically degenerate, which is harmless here because
all three checks are purely combinatorial.

| patch | F | \|E_hinge\| | H | notches | c(G) | pred | rank(L) | dim Z | 1^T L = 0 | rank = H - 1 | L = R D | Z out-harmonic | max harmonic residual |
|---|--:|--:|--:|--:|--:|--:|--:|--:|:--:|:--:|:--:|:--:|--:|
| `torus_squares_4x4` | 16 | 32 | 16 | 0 | 1 | 17 | 15 | 1 | yes | yes | yes | yes | 3.331e-16 |
| `torus_squares_6x4` | 24 | 48 | 24 | 0 | 1 | 25 | 23 | 1 | yes | yes | yes | yes | 7.216e-16 |
| `torus_triangles_4x4` | 32 | 48 | 16 | 0 | 1 | 17 | 15 | 1 | yes | yes | yes | yes | 5.551e-16 |

Note `pred = |E_hinge| - |F| + c(Gamma)` overshoots `H` by one on every torus patch:
`Gamma` is embedded on a surface of Euler characteristic 0, not 2.

### Sweep aggregate (50 random graphs)

| check | holds on |
|---|---|
| 1a. `1^T L` = naive degree difference `indeg(v) - outdeg(v)` | 0 / 50 |
| 1b. `1^T L` = notch-restricted degree difference (corrected form) | **50 / 50** |
| 1c. `1^T L` supported on boundary vertices only | 0 / 50 |
| 2a. `L == R D` exactly (max abs difference 0) | **50 / 50** |
| 2b. `rank(L) == H - dim Z` | **27 / 27** (dense only) |
| 2c. every `y` in `Z` gives an out-harmonic `g`, tol 1e-9 | 27 / 27 **vacuously** (dense only; `dim Z = 0` on all 27, so 0 residuals were evaluated) |
| 3a. `H == |E_hinge| - |F| + c(Gamma)` | 50 / 50 |
| 3b. `H_all == |E_hinge| - |F| + c(Gamma)` (holes + notches) | **0 / 50** |
| Eq. (1) auto-orientation yields `c(Gamma) == 1` | 50 / 50 |

Out-harmonic residuals actually evaluated in the sweep: 0 (worst 0.000e+00).

**Caveat: check 2c is vacuous on the sweep and on the 8 reference cases.** `dim Z == 0` on
every one of the 27 sweep graphs where the left null space was computed, and on all 8
reference cases (the `dim Z` column above), so there is no basis vector to test and the
`27 / 27` in row 2c counts graphs on which nothing was checked, not residuals that passed.
The row is consistent with row 2b: `dim Z == 0` is exactly `rank(L) == H`. The check has
real content only where `dim Z > 0`, which in everything measured means the boundary-free
patches: the three torus rows above and `test/test_rank_checks.jl`, where `dim Z == 1`
and the out-harmonic residual is `3.3e-16 / 7.2e-16 / 5.6e-16`.

### Identities that failed, and the corrected form

**Check 1 (row sum), as first stated, is false on any patch with a boundary.**
`L` carries a row only for an all-interior preimage, so every hinge edge owned by a
boundary-touching preimage (a notch) contributes nothing to `1^T L`. Its target loses
the `+1` and its source loses the `-1`, and the source can be an interior vertex --
so `1^T L` is NOT supported on the boundary either. The corrected identity, which
does hold everywhere measured, restricts both degrees to hinge edges owned by hole
rows. With `I(v) = 1` iff `K(v)` is a row of `L`:

```
    (1^T L)_v  ==  I(v) * indeg_hinge(v)  -  sum_{hinge v->w} I(w)
```

(All in-edges of `v` belong to `K(v)`, so the in-degree is all-or-nothing; the
out-edges are distributed over the preimages of their targets.) This is exactly the
`y = 1` instance of the out-harmonic characterization of check 2. When there are no
notches -- the torus patches in `test/test_rank_checks.jl` -- it collapses to the
naive form, `1^T L == 0` exactly, and `rank(L) == H - 1` (measured, not assumed).

**Check 3 (hinge-graph Euler count) holds exactly as stated**, on all 50 / 50 sweep graphs and all 8 reference cases, with `H` counting only the all-interior preimages:

```
    H  ==  |E_hinge| - |F| + c(Gamma)          (H = holes only, notches excluded)
```

So the notches are NOT bounded faces of `Gamma`: counting them in breaks the identity
(`H_all == pred` on 0 / 50 sweep graphs). This
confirms unverified claim U1 numerically and, with `c(Gamma) == 1` on 50 / 50 sweep graphs under the Eq. (1) auto-orientation,
reduces to Scout-c's `H = |E_hinge| - |F| + 1`.

The one place the plane count needs a correction is a boundary-free patch, where
`Gamma` lives on a surface of Euler characteristic 0 rather than 2 and the count
overshoots by exactly one: `H == |E_hinge| - |F| + c(Gamma) - 1`, verified in
`test/test_rank_checks.jl` on the 4x4 and 6x4 square tori and the 4x4 triangle torus.

Checks 2a and 2b held on every graph and every reference case measured. Check 2c held
wherever it had content -- the three torus patches, `dim Z == 1`, residuals `<= 7.3e-16` --
and was vacuous everywhere else, since `dim Z == 0` on all 27 dense sweep graphs and all 8
reference cases. No counterexample to the factorization or to the out-harmonic
characterization was found, but the out-harmonic characterization has been exercised on
three patches only, and that is the honest strength of the evidence for it.

Counterexamples saved (graph JSON, loadable by `kiri_analyze`):

| case | identity that failed | file |
|---|---|---|
| `hexagons_auto` | rowsum_naive_degree_identity | `violations/hexagons_auto__rowsum_naive_degree_identity.json` |
| `kagome_3636` | rowsum_naive_degree_identity | `violations/kagome_3636__rowsum_naive_degree_identity.json` |
| `periodic_squares_4x4` | rowsum_naive_degree_identity | `violations/periodic_squares_4x4__rowsum_naive_degree_identity.json` |
| `rotating_squares` | rowsum_naive_degree_identity | `violations/rotating_squares__rowsum_naive_degree_identity.json` |
| `snub_square_33434` | rowsum_naive_degree_identity | `violations/snub_square_33434__rowsum_naive_degree_identity.json` |
| `tiling_3_4_3_12` | rowsum_naive_degree_identity | `violations/tiling_3_4_3_12__rowsum_naive_degree_identity.json` |
| `triangles_alternating` | rowsum_naive_degree_identity | `violations/triangles_alternating__rowsum_naive_degree_identity.json` |
| `truncated_square_488` | rowsum_naive_degree_identity | `violations/truncated_square_488__rowsum_naive_degree_identity.json` |

