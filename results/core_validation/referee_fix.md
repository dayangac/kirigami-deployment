# The bisection referee's collision predicate: the F34 hinge-vertex artefact, fixed

Scope: `Kirigami/src/core/collision.jl` -- `polygons_overlap` / `has_collision`, the
primitive under every `theta_bisect` column in the experiments and under `theta_max()`.
All numbers below come from a full re-run of K2a, K5, K6 and K9 against a matched
reverted-predicate build of the same source tree; nothing is quoted from the recorded runs.

## 1. The defect (F34, derivation check R6-c; K9 post-hoc audit)

The old predicate shrank each polygon toward its OWN centroid by a fixed relative amount
(`shrink`, 1e-12 at every experiment call site) and then ran a STRICT segment-crossing test on
the shrunk copies.

Two faces incident to the same hinge edge share the pin vertex: both read the same
`Y[pv]`, so the two copies of that point are bit-for-bit equal, and every orientation
determinant that uses the pin twice is exactly zero -- i.e. the strict crossing test is
exactly right on the raw polygons. `shrink_poly` destroys that: it moves the two copies of
the pin by `-s * (p - centroid_A)` and `-s * (p - centroid_B)`, two DIFFERENT
displacements, so the hinge-incident edge pairs acquire determinants of size 1e-11 ...
1e-17 whose SIGNS then decide the crossing. The sign depends on `s` and on `theta`, so the
misfire is non-monotone in `theta` (K9 delaunay id 13: collision reported true at
`theta` = 1e-8, 1e-7, 1e-5 and false at 1e-9, 1e-6, 1e-4, ...) and it survives arbitrarily
small `s`. Shrinking toward the centroid is also simply wrong when the hinge sector is
reflex, which is the `voronoi_93` case (`alpha_g = 4.105228 > pi`).

Known instances: B4 `voronoi_93` faces (94, 184), where the referee returned 0 against the
exact scan's 0.2484 (a 1200 x 1200 grid finds 0 of 1 442 401 probes inside both faces);
K9 delaunay ids 13, 40, 103, 190, the 4 designs of the 36 exact positives that the referee
rejected at shrink 1e-12 but accepted at 1e-9.

## 2. The fix

`polygons_overlap` no longer shrinks anything. It now decides the question it is supposed
to decide -- do the OPEN interiors meet -- exactly, keeping T4.1b's closed-segment
convention:

1. bounding-box rejection (unchanged);
2. each edge of one polygon is split at every intersection with the other boundary
   (proper crossings, endpoint incidences and collinear overlaps all contribute
   parameters), and the MIDPOINT of every resulting piece is classified against the other
   polygon with a three-state point-in-polygon test (inside / on the boundary / outside).
   A piece strictly inside means the two interiors meet. Done in both directions;
3. if no boundary arc enters the other interior, the only remaining overlap is containment
   with coincident boundaries, so a point constructed strictly inside each polygon (an ear
   construction, no offsets) is classified against the other.

The orientation predicate is accumulated in Float64 in a fixed contraction form
(`core/collision.jl`, `orient_raw`) and snapped to zero by a RELATIVE tolerance: `c` counts as ON the line `ab` when `|det| <= eps * |u| * S` with
`S = max(|c - a|, |a|, |b|, |c|)`, i.e. when its distance from the line is at most `eps`
times the scale of the configuration. It is exact -- `eps` plays no role at all -- whenever
two of the three points coincide, which is exactly the shared-pin case, so the hinge
incidence (`h == 0`, the T3.H.1 exclusion the certificate uses) is handled by exact
cancellation and not by a threshold.

The `|a|, |b|, |c|` terms in `S` are load-bearing and were added after a second, opposite
defect showed up in the re-referee (K5 `voronoi_75`, faces 106 and 390): a tolerance
relative to the edge lengths ALONE rejects a constructed midpoint from its own edge, because
that midpoint carries an absolute ulp of the world coordinate (~4e-15 on a sheet 30 units
across) while the edge is 6e-3 long. Two faces adjacent across a split edge, whose two
copies coincide exactly at `theta = 0`, were then reported as overlapping. Both defects have
the same root cause -- a tolerance whose scale does not match the geometry -- and both are
now pinned by tests.

The third argument keeps its position but changes meaning: it is now that relative tolerance
`tol_rel`, clamped to `[1e-14, 1e-9]`, so all historical call sites (1e-12, 1e-9 and the
1e-6 defaults) agree with each other. The floor exists because below 1e-14 the predicate
would be deciding the rounding noise of its own midpoints. Defaults are 1e-12.

Cost: the predicate is ~3x slower per pair (K2a wall 13.6 s -> 42.4 s for the whole
187-design run, of which the referee bisection is the bulk).

## 3. Tests

New, in `Kirigami/test/test_collision.jl` (written first, failing against the old predicate):

* `polygons_overlap: faces sharing a hinge vertex with a reflex sector do not overlap` --
  the two `voronoi_93` polygons verbatim, at tolerances 0, 1e-15, 1e-12, 1e-9, 1e-6, both
  argument orders. Old predicate: 10 of 10 wrong. New: all correct.
* `polygons_overlap: shared-vertex touching vs. genuine overlap through the pin` -- the
  minimal synthetic version (a reflex sector meeting a convex one at one shared vertex),
  plus the positive control obtained by translating one polygon across the pin.
* `polygons_overlap: collinear and boundary-only degeneracies` -- an exactly shared edge and
  a partially shared edge (no overlap), and a triangle whose vertices ALL lie on the other
  polygon's boundary while the interiors do intersect (must fire).
* `polygons_overlap: faces meeting along a coincident split edge far from the origin` --
  the K5 `voronoi_75` pair verbatim, at five tolerances and in both orders, plus a
  translation invariance check. This is the regression test for the second defect above.

`derivation_tests.jl` R6-c, the pinned reproduction from derivation-check round 6, is INVERTED:
`R6-c2` now asserts that `polygons_overlap` agrees with the dense probe at every tolerance
(0 misfires of 6 values, both orders), and its comment records the fix.

Suites (both re-run after the fix, from the repo root; the `Pkg.test()` totals include test
cases other agents added to the same working tree while this ran):

| suite | cases | assertions | result |
|---|--:|--:|---|
| `Pkg.test()` (whole suite) | 187 test sets | 176,692 (23 marked broken) | 0 failures |
| `derivation_tests.jl` | 38 test sets | 150,191 | 0 failures |

## 4. Re-referee: agreement before / after

Every "before" number here was produced by a REVERTED-PREDICATE BUILD: the current source
tree with only `Kirigami/src/core/collision.jl` reverted to `HEAD`, run from a copy in
the scratchpad, so the comparison isolates this fix from the round-6 corrections other
agents made to `contact.jl` / `zero_plus.jl` in the same working tree. The baseline K5 run
reproduces `results/experiments/k5/summary.txt` line for line (only the wall time differs), which
checks that control.

Outputs live in `results/core_validation/referee_fix/{k2a,pre_k2a,k5,pre_k5,k6,pre_k6,k9,pre_k9}`.

### K2a -- 187 tightly embedded designs, exact T4.2" scan vs the referee

| quantity | before | after |
|---|--:|--:|
| `|Theta_max(T4.2") - bisect(1e-9)| <= 1e-5` | 187 / 187 | 187 / 187 |
| worst gap at 1e-9 | 1.956e-07 | 1.817e-07 |
| `|Theta_max(T4.2") - bisect(1e-12)| <= 1e-5` | 187 / 187 | 187 / 187 |
| worst gap at 1e-12 | 1.956e-10 | 1.817e-10 |
| min-over-roots vs bisection, tolerance 1e-6 | 149 / 187 | 182 / 187 |
| min-over-roots vs bisection, tolerance 1e-12 | 182 / 187 | 182 / 187 |
| certificate columns (holds / POS / NOOVERLAP / NOROOT) | 173 / 187 / 173 / 187 | unchanged |

The headline 187/187 is unchanged and the referee has become tolerance-independent: 1e-6, 1e-9 and 1e-12 now give the same 182/187 (the 5
remaining are the known grazing designs, where the closed-form min-over-roots rule -- not
the referee -- is the quantity being tested). Both verdicts unchanged.

### B4 -- the F34 case itself (`derivations/scratch/check_b4_93.jl`, ids 93 and 96)

| design | exact scan | referee before | referee after |
|---|--:|--:|--:|
| `voronoi_93` (sigma_def, free) | 0.248400 | 0.000000 | **0.248400** |
| `voronoi_96` (sigma_def, free) | 0.241884 | 0.241884 | 0.241884 |

For id 93 the old predicate fired on face pair (94, 184) at every probe from 1e-7 to 0.01
and on a second pair (476, 583) at 0.05; after the fix it fires nowhere below the exact
first contact, and the first pair it reports at 0.2484 is (151, 296), the pair the scan
names. F34 is closed: the referee now reproduces the exact scan on the case that exposed it.

### K5 -- 200 graphs, the orientation-defect experiment

| quantity | before | after |
|---|--:|--:|
| `sigma_def` gives an overlap-free `X0` at `theta = 0` | 114 / 200 | **125 / 200** |
| ... as a rate on the graphs where `sigma_mc` fails | 57.00 % | 62.50 % |
| `sigma_mc` gives an overlap-free `X0` | 0 / 200 | 0 / 200 |
| smallest ladder angle with an interior overlap, worst over graphs | 1e-12 | 1e-10 |
| every other line of the summary | -- | unchanged |

Eleven designs whose flat state the old predicate called self-overlapping are in fact
embedded, and none goes the other way. The ladder line changes too: the old predicate found
an interior overlap at the first rung, `1e-12`, on every graph that overlaps at all; the new
one finds it at `1e-12` on 36 of the 125 and at `1e-10` on the other 89. The conclusion is
unaffected -- the flat sheet still self-intersects as soon as the cuts open -- and the
PASS/FAIL verdict (bar: >= 20 %) does not move. `theta_max` itself (`def_theta`, `mc_theta`)
is 0 on all 400 designs before and after.

### K6 -- 400 designs, the 0+ repair (primary variant, `--no-native`, K6's own cache)

| quantity | before | after |
|---|--:|--:|
| exact `Theta_max > 0` | 0 / 400 | **2 / 400** |
| refereed `Theta_bisect > 0` | 0 / 400 | **2 / 400** |
| `|theta_exact - theta_bisect| <= 1e-5` | 400 / 400 | 400 / 400 |
| worst gap | 0 | 2.0e-08 |
| soundness `eps_max <= theta_bisect` | 400 / 400 | 400 / 400 |

The two designs are `voronoi 12 / sigma_def` (`Theta_max = 0.002937`) and
`voronoi 93 / sigma_def` (`Theta_max = 0.238947`, the F34 design itself). They move because
`polygons_overlap` is called INSIDE the exact scan too (`exact_theta_max_overlap`'s interval
probe, `contact.jl`), not only in the referee: the same hinge artefact was rejecting the
first interval. Both are still uncertified (`cert_noroot = 0`), and 2 / 400 = 0.5 % is far
under K6's 20 % bar, so the verdict does not move -- but the sentence "`Theta_max = 0` on all
400" is no longer true and is corrected in `EXPERIMENTS.md`.

Control: the reverted-predicate build reproduces `results/experiments/k6/k6_final.csv` on every column
except `sec_eps_max` on 15 rows, which is another agent's round-6 change to
`method/zero_plus.jl` in the same working tree, not this fix.

### K9 -- 400 designs, the convexity + split-inward embedding

| quantity | before | after |
|---|--:|--:|
| exact `Theta_max > 0`, arm (b) | 36 / 400 | 36 / 400 |
| **refereed `Theta_bisect > 0`, arm (b)** | **32 / 400** | **36 / 400** |
| refereed, `sigma_mc` / `sigma_def` | 20 / 12 | 24 / 12 |
| `|theta_exact - theta_bisect| <= 1e-5` | 396 / 400 | **400 / 400** |
| worst gap | 3.020e-01 | **0** |
| soundness `eps_max <= theta_bisect` | 396 / 400 | 400 / 400 |
| arms `X0` and (a) | 0 exact, 0 refereed, 400/400 agree | unchanged |

The only four cells that changed in the whole 400-row table are the `b_theta_bisect` entries
of delaunay 13, 40, 103 and 190 (all `sigma_mc`), i.e. exactly the four F36 disagreements, and
each now returns the exact value: 0.219029, 0.177931, 0.301957, 0.103354. No other column of
any row differs.

## 5. Verdicts that move

**None.** Every experiment verdict is unchanged: K2a PASS with the T4.2" scan (187 / 187),
K5 PASS on its own bar, K6 FAIL (2 / 400 = 0.5 % against a 20 % bar), K9 FAIL (36 / 400 =
9.0 % against a 10 % bar), B4 unchanged (`voronoi_93` was already reported at the scan's
0.2484; only the referee column moves).

What does change is the wording that has to be used about them:

1. K9's headline is now **36 / 400 exact and 36 / 400 refereed** -- the "refereed 32 at shrink
   1e-12, 36 at 1e-9" split disappears, and with it the shrink-dependence caveat. The
   `9.0 %` figure and the `FAIL by 4` conclusion stand.
2. K6 is no longer "`Theta_max = 0` on all 400": it is 2 / 400 with `Theta_max > 0`, both
   uncertified, both `sigma_def` voronoi.
3. K5's `sigma_def` embedding rate is 125 / 200 = 62.5 %, not 114 / 200 = 57.0 %, and its
   overlap ladder bites at `1e-10` rather than at `1e-12` on 89 of those designs.
4. The referee is now tolerance-independent over the whole historical range (1e-12 to 1e-6),
   so it may be described as an independent cross-check without the shrink caveat -- though
   the exact T4.2" scan remains the ground truth and the referee remains a bisection on a
   4000-point grid, i.e. resolution-limited, not exact.
5. F34 is closed. The old predicate's bias was one-sided in the runs measured here -- every
   one of its 4 + 2 + 11 disagreements is a collision it reported that is not there -- so
   every FAIL verdict it produced was conservative, which is why nothing flips. That is a
   measurement, not a theorem: shrinking a polygon can in principle also HIDE an overlap
   thinner than the shrink, and the interim version of this fix (before the scale term was
   added to the tolerance) produced exactly one such spurious flip, which is now a test.
