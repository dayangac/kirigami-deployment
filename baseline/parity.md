# Parity: our `code/src/core` vs. the authors' published implementation

Everything below is measured with `baseline/native/build/tuttekiri_cli` running the
authors' own C++ on the eight reference graphs in `results/core_validation/cases/*/M.json`,
with our `code/build/kiri_analyze` / `kiri_deploy` on the same inputs. Our reference
numbers are the ones already recorded in `results/core_validation/reference_cases.md`.
Build details and the sigma <-> colour correspondence are in `baseline/README.md`.

## 1. Cuts, holes, rank, projection

`bnd` = boundary vertices. `E_h`/`E_s` = hinge / split edge counts, `set=` is set equality of
`E_hinge` as unordered vertex pairs (our `cut.json` against their classification).
`rows` is the number of rows their `make_deployable` builds (their printed `num_holes`);
`rank` is `rank(L)` for us and `nv - dim ker A` for them.
`res in theirs` is the residual of **our** input `M` in **their** assembled system;
`res after` is the residual of their own projected `X0` in their own system.
`FK` is the max deviation, after rigid alignment, between the two deployed `M'` at 30 deg.

| case | nv | nf | int | bnd | E_h ours | E_h theirs | E_s ours | E_s theirs | set= | H ours | rows theirs | rank ours | rank theirs | res in theirs | res after | FK |
|---|--:|--:|--:|--:|--:|--:|--:|--:|:--:|--:|--:|--:|--:|--:|--:|--:|
| hexagons_auto | 33 | 10 | 9 | 24 | 13 | 13 | 5 | 5 | yes | 4 | 28 | 4 | 6 | 1.73e+00 | 1.66e-15 | 1.3e-08 |
| kagome_3636 | 123 | 97 | 78 | 45 | 174 | 174 | 0 | 0 | yes | 78 | 123 | 78 | 79 | 2.84e-15 | 2.30e-14 | 1.6e-08 |
| periodic_squares_4x4 | 25 | 16 | 9 | 16 | 24 | 24 | 0 | 0 | yes | 9 | 25 | 9 | 10 | 0.00e+00 | 3.33e-15 | 1.6e-08 |
| rotating_squares | 36 | 25 | 16 | 20 | 40 | 40 | 0 | 0 | yes | 16 | 36 | 16 | 17 | 0.00e+00 | 1.05e-14 | 2.0e-08 |
| snub_square_33434 | 50 | 53 | 25 | 25 | 64 | 64 | 13 | 13 | yes | 12 | 37 | 12 | 14 | 3.24e+00 | 5.99e-15 | 2.7e-01 |
| tiling_3_4_3_12 | 56 | 25 | 20 | 36 | 40 | 40 | 4 | 4 | yes | 16 | 52 | 16 | 17 | 3.59e+00 | 6.44e-15 | 3.6e-01 |
| triangles_alternating | 58 | 89 | 33 | 25 | 121 | 121 | 0 | 0 | yes | 33 | 58 | 33 | 34 | 2.04e-15 | 1.48e-14 | 3.9e-08 |
| truncated_square_488 | 56 | 21 | 24 | 32 | 32 | 32 | 12 | 12 | yes | 12 | 44 | 12 | 13 | 8.31e-16 | 3.72e-15 | 1.7e-08 |

### What matches

* **`E_hinge` and `E_split` agree exactly, as sets, on 8/8 cases.** Symmetric difference 0.
  The cut classification is not in doubt.
* **The Eq. (6) projection is exact in both.** Their `X0` satisfies their system to 6e-15
  and satisfies **our** system to 4e-15 on 8/8 (measured by feeding their `X0` to
  `kiri_analyze`). Their KKT formulation `[[I, A^T],[A, 0]]` solved with `SparseQR` is the
  same least-norm projection our dense-SVD route computes.
* **Forward kinematics agrees to 4e-8** (their `deploy` takes a `float` angle, which is the
  size of the residual) on **8/8** cases when both are applied to a *deployable* embedding,
  i.e. to `X0`. The two cases with a 0.27-0.36 FK gap in the table above are exactly the two
  our reference file already marks "deployable as given: no" (`snub_square_33434`,
  `tiling_3_4_3_12`); for a non-deployable `M` the BFS visiting order changes the answer, so
  the disagreement is the expected symptom, not an implementation difference. Re-run on
  `X0` the same two cases give 3.0e-8 and 1.8e-8.

### What does not match, and why

**Their constraint system is strictly stronger than ours on 2 of 8 cases.** `rank theirs`
is `rank ours + 1` on 6/8 (the +1 is their pin of vertex 0), but `rank ours + 2` on
`hexagons_auto` and `snub_square_33434`. The cross-residuals localise it exactly:

| direction | hexagons_auto | snub_square_33434 | other six |
|---|--:|--:|--:|
| their residual of **our** `X0` | 1.73e+00 | 1.00e+00 | <= 4e-15 |
| our residual of **their** `X0` | 4.00e-16 | 1.13e-15 | <= 4e-15 |

Our embeddings do **not** satisfy their system; theirs always satisfy ours. The cause is
our documented Deviation 2 (`code/README.md`): we drop every hole preimage that touches the
mesh boundary, on the argument that those are notches open to the exterior rather than
holes. The authors do not drop them. Their `make_deployable` builds one row per *vertex
class* of the split-edge subgraph (hence `rows theirs` = 36, 44, 52, ... , far more than
`H`), and a class that contains both interior and boundary vertices still collects the
hinge edges that cut into its **interior** vertices. On these two patterns that produces one
extra, independent constraint that our `L` omits. Both cases are exactly the ones where a
split-cut component reaches the boundary.

This is a real disagreement about the mathematics, not a numerical one, and it means our
shape space is one dimension too large on such patterns. It should be resolved before any
rank claim built on `L` is asserted.

**Everything else about the rank difference is bookkeeping.** Their system pins a single
vertex instead of fixing the whole boundary, so `dim ker A` per coordinate is
`nv - rank(L) - 1` and is dominated by free boundary motion; it is not comparable to our
`dim_null` under fixed boundary conditions and no conclusion should be drawn from the raw
numbers. Their `num_holes` is a count of vertex classes, not of holes: rows belonging to
classes made only of boundary vertices come out identically zero.

**Their hole tracer could not be run at all.** `UnitPattern::get_holes()` dereferences
`cur_edge->prev()->twin()` without a null check, so it segfaults on any patch with a
boundary; all eight cases have one. `n_traced_holes` is `-1` throughout. Their hole *count*
is therefore not directly comparable; only the constraint rows are.

## 2. Face-orientation assignment

Their `coloring::initialized_two_face_coloring` (circular relaxation of Eq. (1), then the
best of 100 evenly spaced diameters against a greedy 2-colouring) versus the `sigma` we
stored. `agree` is up to a global flip.

| case | agree | their hinge | our hinge | their split | our split |
|---|--:|--:|--:|--:|--:|
| hexagons_auto | 9/10 | 13 | 13 | 5 | 5 |
| kagome_3636 | 97/97 | 174 | 174 | 0 | 0 |
| periodic_squares_4x4 | 16/16 | 24 | 24 | 0 | 0 |
| rotating_squares | 25/25 | 40 | 40 | 0 | 0 |
| snub_square_33434 | 43/53 | 64 | 64 | 13 | 13 |
| tiling_3_4_3_12 | 25/25 | 40 | 40 | 4 | 4 |
| triangles_alternating | 89/89 | 121 | 121 | 0 | 0 |
| truncated_square_488 | 21/21 | 32 | 32 | 12 | 12 |

The two implementations reach **the same objective value** on 8/8 (identical hinge and
split counts), including the two where the assignment itself differs. Our Eq. (1) solver is
not worse than theirs on any case measured. Note their routine seeds from `Eigen::Random`
without seeding `srand`, so it is deterministic per process but arbitrary.

## 3. theta_max

`kin` is their `kirigami::max_opening_angle` (the kinematic bound only); `coll` adds their
1%-step forward collision scan, i.e. `theta_max` exactly as `bind.cpp` computes it for the
web UI. Both are evaluated on the projected `X0` so that the pattern is actually deployable.
`ours` is the `theta_max` column of `reference_cases.md`, `min beta` its neighbour.

| case | their kin | our min beta | their coll | our theta_max |
|---|--:|--:|--:|--:|
| hexagons_auto | 2.094395 | 2.094395 | 2.094395 | 2.094397 |
| kagome_3636 | 3.141593 | 3.141593 | 3.141593 | 3.141593 |
| periodic_squares_4x4 | 3.141593 | 3.141593 | 3.141593 | 3.141593 |
| rotating_squares | 3.141593 | 3.141593 | 3.141593 | 3.141593 |
| snub_square_33434 | 3.383898 | 3.383898 | **0.067307** | 1.646137 |
| tiling_3_4_3_12 | 2.617994 | 2.380636 | 2.356194 | 2.380638 |
| triangles_alternating | 4.188790 | 4.188790 | 4.188790 | 3.141593 |
| truncated_square_488 | 2.356194 | 2.356194 | 2.356194 | 2.356196 |

* **Their kinematic bound is our `min beta`, to 2e-6, on 7/8.** The two quantities are the
  same object. The only exception is `tiling_3_4_3_12` (2.617994 vs 2.380636), where their
  bound is the looser one; our `theta_max` there (2.380638) is confirmed by their own
  collision scan (2.356194, one grid step below).
* **We cap `theta_max` at pi; they do not.** `triangles_alternating` and
  `snub_square_33434` are the cases where the bound exceeds pi.
* **Their collision-aware `theta_max` is wrong on `snub_square_33434` (0.0673 vs 1.6461).**
  Their scan applies `merge_close_verts()` with tolerance `0.1 * avg_edge_len` to the
  deployed mesh before the overlap test, which fuses hinge duplicates that have opened by
  less than that tolerance and creates a spurious overlap. A 60-point scan of their own
  predicate over `(0, 3.3839]` (`tuttekiri_cli collide`) is clean up to 1.6919 and first
  reports contact in the bracket **1.6919-1.7483**, which contains our 1.646137. The 0.0673
  hit is an isolated artefact at step 3 of a 200-point scan. Our number is the correct one.

## 4. Their collision-prevention optimization (Eq. 9)

`opt::prevent_intersections` run on their own `X0`, with the web UI's defaults
(`src/state.js`: barrier 0.1, barrier strength 10, close-to-original 0.1) and over the
84-point grid `barrier x strength x close = {0,0.05,0.1,0.2,0.3,0.5,1} x {1,10,50,100} x
{0,0.1,1}`. Because their own collision detector is unreliable (above), the achieved
`theta_max` is re-measured with **our** independent detector by feeding their optimised
pattern back through `kiri_analyze --collision`.

| case | start (their X0) | their defaults | best of 84 | our `theta_max` of their result | our `min beta` of their result | ours (`reference_cases.md`) |
|---|--:|--:|--:|--:|--:|--:|
| snub_square_33434 | 3.365357 | 3.143889 | 3.377903 (b 0.1, s 10, w 1.0) | **3.141593** (pi cap) | 3.143889 | 1.646137 |
| truncated_square_488 | 2.356194 | 2.353852 | 2.356721 (b 0.3, s 50, w 0.0) | 2.356721 | 2.356721 | 2.356196 |
| hexagons_auto | 1.640229 | 1.720933 | not swept | -- | -- | 2.094397 |
| tiling_3_4_3_12 | 2.274526 | 2.281580 | not swept | -- | -- | 2.415428 |

All 168 sweep runs produced finite vertices; total optimiser time 16.9 s (snub) and 9.1 s
(truncated). The residual of every optimised pattern in **our** constraint system is
<= 1.4e-15, so their Eq. (9) stays inside the deployable subspace, as it must.

### The headline result

**On `snub_square_33434` the authors' Eq. (9) roughly doubles `theta_max` where ours could
not improve it at all.** Our `code/README.md` records "on the snub-square pattern no value
of gamma improves it", with `theta_max` stuck at 1.646137. Their optimisation, at its
published default parameters, produces a pattern that **our own** collision detector
certifies to `theta_max = 3.141593` (our pi cap; the kinematic bound of that pattern is
3.143889). This is a genuine gap in our implementation of Eq. (9), not a measurement
artefact: the pattern was produced by their code, checked by our code, and the two
independent detectors agree on it.

On `truncated_square_488` the two are level: their defaults give 2.353852 against our
2.356196, and their best sweep point gives 2.356721 — a 0.02% improvement. Our
README's claim that Eq. (9) "raises `theta_max` from 0.059 to 2.027" on this pattern refers
to a different instance than the reference case, whose baseline is already 2.356196; the
`2.027` figure could not be reproduced here and should be re-checked.

The likely reason for the snub-square gap is a difference in what Eq. (9) optimises. Theirs
is parameterised by the null-space coefficients and its barrier acts on `trace_hole`, the
sign of the gradient of the relative motion of the two duplicates of each **split** edge —
a first-order, collision-free quantity. Ours (Deviation 6/7 in `code/README.md`) invents a
barrier `B(s) = kappa log(1 + exp(-s/kappa))` on a different argument and sweeps a `gamma`
ladder. The authors' formulation is the one that works here.

## 5. Fully closed

`opt::optimize_for_fully_closed` on their `X0`. `theta_max` before/after is their kinematic
bound.

| case | dim ker (per coord) | before | after | seconds |
|---|--:|--:|--:|--:|
| hexagons_auto | 27 | 1.640229 | 2.178433 | 1.2 |
| snub_square_33434 | 36 | 3.365357 | 3.663208 | 41.1 |
| truncated_square_488 | 43 | 2.356194 | 2.356194 | 1.5 |
| tiling_3_4_3_12 | 39 | 2.274526 | 2.276595 | 8.2 |
| rotating_squares | 19 | 3.141593 | 3.141593 | 0.1 |

All runs finite. We have no counterpart implementation, so this is a baseline record only.

## 6. Timings

`make_deployable` (Eq. 6, their `SparseQR` KKT solve) on these sizes: 0.25-1.4 ms.
`prevent_intersections`: 26-160 ms per run at defaults. `optimize_for_fully_closed`:
1.2-41 s. Their collision scan dominates every `theta_max` query at 0.1-3 s.

## 7. What could not be checked

* Their hole tracer (`get_holes`) — segfaults on bounded patches, see above.
* Periodic boundary conditions in their formulation. `periodic_squares_4x4/M.json` carries
  no `periodic` block, so the baseline treated it as a bounded patch with a single pinned
  vertex; their Eq. (3b)-(3c) path (`make_periodic` + the periodic rows of `A`) is built and
  compiles but was not exercised on a case where we have a number to compare against.
* The inverse-design pipeline (`param::lift`, `opt::optimize_rigidity`), which needs a 3D
  target mesh.
* Their `conformalize` (Eq. 13), which we have no counterpart for.
