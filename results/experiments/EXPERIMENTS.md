# Experiments

Each section below is one experiment: the question, how it was run, what was measured,
and what it settles. All numbers were produced by the drivers `Kirigami/apps/exp_*.jl` on
Apple silicon with the numerical conventions of `docs/NUMERICS.md`. Populations are
deterministic functions of the graph id and are frozen in `data/corpus/`
(`data/corpus/README.md`), so every experiment sees the same graphs. Where an experiment
had a pre-registered pass bar, the bar and the outcome against it are both stated.

Labels: `K1a … K9c`, `A3`, `B3`, `B4`, `F23`, `T-1`, `E1`, `Native200` are the experiment
ids used by the CSVs, figures and the papers. `F<n>` labels refer to defects found and fixed
during the project (F24, F32, F34, F42 are the ones that matter here).

Contents: conventions · summary · K3a · K1b · K1a · K2c · K2a · K1c · K5 · K2b · F23 ·
rank checks · K6 · A3 · K7 · B4 · K8a · B3 · T-1 · K9 · K9b · K9c · Native200 · reproducing.

---

## Conventions and corrections that affect the numbers

These were established during the derivation checks (`notes/derivations/core.md`, `check.md`) and are implemented in `Kirigami/src/method/` with unit tests in
`Kirigami/test/test_method.jl`.

1. **`Θ_max` is not the minimum over harmonic roots** (T4.2″). A root can be a *graze*, a
   tangency without a crossing: `hexagons_auto` grazes at `π/3 = 1.047198` and does not
   overlap until `2π/3 = 2.094395`. `Θ_max` is computed by the T4.2″ interval scan
   (`exact_theta_max_overlap`): take the complete contact set `C(X)`, walk the gaps, return
   the left endpoint of the first gap whose midpoint has an interior overlap. This changed
   K2a from fail to pass.

2. **The `τ = 0` deflation** (T5.2b.2). With `A = p − q`, `B = 2r`, `C = p + q`,
   `g(τ) = C + Bτ + Aτ²`, `τ = tan(θ/2)`: `C ≠ 0` is class 1; `C = 0, B ≠ 0` is class 2, a
   simple root at `θ = 0`; `C = 0, B = 0, A ≠ 0` is class 3, `h = p(1 − cos θ)`, a double root
   at `θ = 0` and constant sign on `(0, π)` — never a contact; `A = B = C = 0` is a permanent
   incidence. `C = 0` holds identically on the whole shape space for every permanent
   incidence and every split-edge duplicate pair. Without the deflation the flat-state root
   leaks in at `~1e−8` and every pattern with split cuts reports `Θ_max = 0`. On K2a's 187
   configurations the deflation removes 13 630 spurious roots in `(0, ε)`; the deflated root
   list is exact on 60 000 synthetic harmonics of all three classes.

3. **The swept radius is exact, without a `√2`.** In the face's own frame the per-face
   translation cancels and `S_u − ḡ_f = −σ_f J (C_u − ḡ_f)`, so the trajectory is a circle
   of radius `‖x_u − x̄_f‖` exactly. The broad phase uses that radius with the
   moving-centroid distance test (T4.5a). Surviving pairs per face drop from a median of
   22.23 to 12.69, `Θ_max` unchanged.

4. **The flat-centroid pruning is unsound in principle.** Deployment contracts centroid
   distances while face radii stay fixed, so a flat-centroid test can discard a pair whose
   faces approach (19 842 such pairs measured). The moving test is the default.

5. **`Θ_max = min(min_e β_e, π)` on split-free patterns**, not `min_e β_e`: the triangle
   tiling has `min β = 4π/3 > π`. 16/16 split-free configurations satisfy the capped identity
   to `8.882e−16`; the cap bites on 5.

6. **The `τ` chart is singular at `θ = π`.** `harmonic_roots` uses the amplitude/phase
   form, and the class-2 branch handles `A = 0` explicitly as a root at `θ = π`. All searches
   are capped at `π`.

7. **The validity certificate** (T5.2b′): `VALID(ε) = POS ∧ NOOVERLAP(θ₁ = ε/2) ∧ NOROOT(ε)`.
   `POS`: all face signed areas `> 0` at `θ = 0`; `NOOVERLAP`: one exact polygon–polygon test
   at the single interior angle `ε/2`; `NOROOT`: no admissible deflated root in `(0, ε)`. A
   `θ = 0` overlap test is not used: at `θ = 0` the two copies of every split edge coincide
   and adjacent faces touch along whole edges, so the test is degenerate. On K2a's population
   the certificate implies `Θ_max ≥ ε` on 173/173 configurations where it holds, 0
   violations, and is exact there (173/187 have `Θ_max ≥ ε`, 173 are certified). Before F32
   (`NOROOT` lacked the two interval tests) it certified only 65/187; see §A3 and
   `jitter/cert_diagnosis.md`.

8. **H-LOC is refuted.** No `O(n)` active-set claim is made; see K2c.

9. **The bisection referee needs a fine grid.** At 180 samples the bisection steps over a
   genuine overlap window of width `2.9e−3` rad at `θ = 0.9467` on `trunc_square_R20_s2` and
   reports `2.025` where the closed form says `0.9467` (confirmed by a `1e−4` scan,
   `apps/dbg_k2a_out.jl`). All referee bisections use 4000 grid points.

### Which predicate each number uses

| number | predicate |
|---|---|
| K2a's certificate columns, the T5.2b′ soundness check | `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT`, `ε = 0.006` |
| every `Θ_max` reported anywhere | T4.2″ interval scan, refereed by bisection on the exact overlap predicate, 4000-point grid. After F34 (`results/core_validation/referee_fix.md`) the predicate no longer shrinks the faces and its tolerance is a relative degeneracy tolerance, so the referee agrees at `1e−12`, `1e−9` and `1e−6` alike |
| K1a's `p_valid` and "X0 injective" | `n_inv = 0` and no overlap at `θ = 0` — the earlier, since-withdrawn predicate. K1a is quoted as measured |
| K2b's "valid ladder point", K1c's "certified design" | no overlap at `θ = 0` — the withdrawn predicate (both ran before the certificate was defined) |
| K5's "overlap-free X0" | no overlap at `θ = 0`, for comparability with K1a; K5 also reports the certificate, and the difference is the finding |

## Summary

| experiment | outcome | key number |
|---|---|---|
| K3a | identity holds; second clause fails | identity 505/508 as measured, every exception a sparse-rank estimate above the dense limit, 7/7 rechecked ids hold exactly with a dense rank. `m_core ≤ 5` on 0 of 167 Delaunay patches (median `m_core = 305`). `σ ∈ ker A` on 296/296 at `X0`, 0/212 at `X_ini` |
| K1b | pass | worst relative fit residual `1.450e−12` over 207 096 triples (bar `< 1e−10`) |
| K1a | the shape space of a random graph is empty near the projection | `n_inv(X0) > 0` on 142/200; `p_valid = 0` on 134/200; an injective sample on 1/200 graphs, `Θ_max > 0` on 0/200 |
| K2c | pass on the corrected rule; H-LOC refuted | pruned == unpruned exact `Θ_max` on 285/285, worst difference 0. The original gate statistic is identically 1 by algebra |
| K2b | fail | median relative gain over the authors' native baseline 0.0000 on 30 live designs (bar `≥ 25 %`); ours better on 7, native better on 15, tie on 16 |
| K1c | pass | corrected re-closure rate 9.47 % on the authors' native Eq. (9) (bar `≥ 3 %`). The rule as first written gives 65.09 %, and with the crossing clause at its own reference range exactly 0, because that rule is self-contradictory |
| K2a | pass with T4.2″ | `\|Θ_max − bisection\| ≤ 1e−5` on 187/187, worst `1.817e−10`; min over roots gets 182/187, worst error `1.047` rad. Certificate implies `Θ_max ≥ ε` on 173/173 |
| K5 | fail | `σ_def` gives a certified `X0` (`ε = 0.3`) on 0/200, as does `σ_mc`. It does fix the flat sheet: defect down 122×, `POS` from 58/200 to 186/200, projection distance down 4.7×, 62.5 % pass the withdrawn `θ = 0` test |
| K6 | fail | certified `Θ_max > 0` on 2/400 under the primary repair, 15/200 `σ_mc` (7.5 %) under the split-only secondary; 17/400 taking the best per design, 15 with `ε_max ≥ 0.1` rad, every one by leaving the projection's neighbourhood |
| A3 | mixed | exact `Θ_max > 0` transition sharp at `a* = 0.16–0.52` median edges on the 4 split-bearing tilings (4 of 8 are split-free: shape space is a point); certificate certifies 1856/1859 rows with `Θ_max ≥ 1` rad after F32 (1622 false negatives before it); best geometric predictor AUC 0.778 (bar 0.9) |
| K7 | C1/C2/C4 pass, C3 fail | `J(θ) = cos(θ/2)I + sin(θ/2)K` exact to 7.3e−14, `K` affine to 4.8e−14, `dim 𝒦 = 2·rank D` on 33/33; conformal at every θ (≤1e−15, 25/25); `θ_c = 2·atan2(tr K, 1 − det K)` to 5.6e−14 on 39/39; target-driven certified design 13/25 = 52 % (bar 80 %) — the obstruction is the 0⁺ collision |
| K9 | fail by 4 designs; first non-zero result near the projection | convexity + split-inward constrained embedding: exact `Θ_max > 0` on 36/400 = 9.0 % (refereed 36/36 after F34), certified `ε_max` median 0.25 rad, max 2.0, `≥ 0.1` on 33/36; bar ≥ 40. Every baseline on the same 400: Eq. (6) 0, `σ_mc` 0, `σ_def` 0, K6 2 (primary) / 15 (split-only, far from the projection), B4 2, K8a 0, native `prevent` 0. Convexity alone: 118 feasible, 0 deployable |
| K9b | fail | the solver push gives 28/400 (7.0 %), below K9's 36; best over both 37/400, 34 with `ε_max ≥ 0.1` rad. Only the `δ` sweep paid. The 9 K9 positives lost are exactly margin-feasible and still `Θ_max = 0`: the cap is basin selection, not budget |
| K9c | pass | replacing the proximity objective by the 0⁺ margin inside the same barriers: exact `Θ_max > 0` on 312/400 (78.0 %); 310 certified, 192 with `ε_max ≥ 0.1` rad (201 best-of-3), `ε_max` median 0.162 (0.172 best-of-3), max π. Optimiser-path dependent at ±2 % (an earlier pass gave 307/306/214). Proximity arms on the same 400: 36 and 27. All 9 of K9b's lost designs recovered. Soundness 400/400 |
| Native200 | fail (partial coverage) | the authors' full pipeline on the K9 population: exact `Θ_max > 0` on 0/195 completed cells; their own test says `Θ > 0` on 27/195, all false negatives of the F24 mechanism. 574/600 dispatched: 195 completed, 364 timed out, 15 crashed. Median `\|F\|` 202 (completed) vs 332.5 (population) vs ≤ 97 (paper's examples) |
| F23 | confirmed | forward kinematics from all 53 seed faces and 20 BFS orders agree to `6.24e−15` |

---

## K3a — the 2-core mobility identity

**Question.** Does `dim ker A = |F \ core2(Γ)| + dim ker A|core2` hold as an integer
identity, and is the mobility of a free-boundary Delaunay patch a dangling-face artefact
(`m_core = dim ker A|core2 − c(core2) ≤ 5` on ≥ 90 % of patches)?

**Measured** (`k3a/{k3a.csv,summary.txt}`, 508 graphs, `F ∈ [10, 4918]`, wall 44.5 s; ranks
by column-pivoted QR with threshold `1e-10 ‖A‖_2`):

| quantity | value |
|---|---|
| identity holds | 505 / 508 |
| identity at `θ = 0.3·Θ_max` | 8 / 8 |
| `σ ∈ ker A` (Eq. (2)) | 296 / 508 |
| Delaunay patches with `m_core ≤ 5` | 0 / 167 (0.00 %) |
| Delaunay median `m_full` / median `m_core` | 322 / 305 |
| Delaunay `m_full` range / `m_core` range | [39, 1875] / [36, 1675] |

The three violators (`delaunay_88`, `F = 1644`; `delaunay_247`, `F = 3103`;
`quad_random_266`, `F = 1079`) all sit above `mobility.jl`'s `dense_limit = 700`, where
`matrix_rank` switches to a sparse QR whose rank estimate is unreliable near the threshold;
which large graphs trip it depends on the backend (an earlier pass reported 15 exceptions,
every one off by exactly ±1). `exp_k3a_core_mobility_recheck.jl` recomputes both ranks
densely on the ids of that earlier list, with K3a's exact configuration (`X0` where
`N ≤ 1400`, else `X_ini`):

| id | kind | F | X | sparse-QR defect | dense defect | smallest kept `σ_i / σ_0` |
|---|---|--:|---|--:|--:|--:|
| 257 | quad_random | 1396 | X0 | 0 | 0 | 5.913e−05 |
| 266 | quad_random | 1079 | X0 | −2 | 0 | 1.886e−04 |
| 320 | quad_random | 1358 | X0 | 0 | 0 | 6.678e−04 |
| 425 | quad_random | 1428 | X0 | 0 | 0 | 5.869e−05 |
| 445 | delaunay | 1138 | X0 | 0 | 0 | 1.698e−03 |
| 467 | quad_random | 911 | X0 | 0 | 0 | 6.131e−05 |
| 494 | quad_random | 879 | X0 | 0 | 0 | 2.632e−04 |

7 of 7 hold exactly with a dense rank. The singular gap at the cut is `5.9e−5` to
`1.7e−3`, three to eight orders above the threshold, so the sparse estimate is simply wrong
there. The violators with `F` from 1644 to 4902 are too large for a dense QR and are not
claimed either way.

**Reading.**
* The identity is true and `A` is correctly assembled; the exceptions are the rank
  estimator. No rank above 700 faces should be quoted without a dense recheck; every
  rank-derived number in the table above `F = 700` carries that caveat.
* The dangling-face hypothesis is refuted outright: deleting the dangling faces moves the
  median mobility from 322 to 305 (about 5 %). The mobility of a free-boundary Delaunay
  patch is intrinsic to its 2-core.
* `σ ∈ ker A` is 296/508 only because two configurations were mixed: 296 of 296 where `A`
  was evaluated at the Eq. (6) projection `X0` (worst residual `3.595e−12`), 0 of 212 where
  the solve was skipped (`N > 1400`) and `A` was built at `X_ini` (residuals up to
  `6.373e+02`). `σ ∈ ker A` *is* Eq. (2), so it holds exactly when the configuration is
  uniformly deployable. The correct statement is 296/296 at `X0`.

Artifacts: `k3a/{k3a.csv, summary.txt, k3a_recheck.csv, recheck_summary.txt}`. Drivers
`exp_k3a_core_mobility.jl`, `exp_k3a_core_mobility_recheck.jl`.

---

## K1b — the harmonic identity

**Question.** Is every orientation determinant along the uniform deployment a single first
harmonic `p + q cos θ + r sin θ`, and are faces rigid? This gates the whole closed-form
contact calculus.

**Measured** (208 graphs, `F ∈ [10, 793]`, 20 random `X` in the shape space per graph, 50
random (vertex, edge) triples per design sampled at 200 angles on `(0, min(π, θ_max))`,
207 096 triples, wall 42.2 s):

| quantity | value | bar |
|---|---|---|
| worst relative LS residual of `p + q cos + r sin` | 1.450e−12 (`voronoi_129`) | `< 1e−10` |
| triples with residual / max\|det\| ≥ 1e−10 | 0 | |
| worst residual / (\|AB\| \|AP\|) | 1.425e−12 | |
| worst \|closed form − `deploy()`\| / (\|AB\| \|AP\|) | 1.442e−12 | |
| worst \|fitted − closed-form\| coefficients | 3.141e−13 | |
| worst relative \|`deploy` − (cos C + sin S)\| | 1.854e−15 | |
| worst relative face signed-area drift | 1.002e−13 (`quad_random_5`) | `1e−12` |
| triples skipped as numerically coincident | 14 (both \|AB\| and \|AP\| below `1e−6` of the diameter) | |

**Pass.** The trig-linear form `Y(θ) = cos(θ/2) C + sin(θ/2) S` and the harmonic structure
of every orientation determinant hold to `1.5e−12`; faces are rigid to `1.0e−13`.

Artifacts: `k1b/{k1b.csv,summary.txt}`. Driver `exp_k1b_harmonic_identity.jl`.

---

## K1a — is the Eq. (6) projection valid, and is the shape space usable?

**Question.** On random planar graphs, is the authors' Eq. (6) projection `X0` a valid
embedding, and does isotropic sampling of the null space around it find valid, deployable
designs?

**Measured** (200 graphs, `F ∈ [101, 793]`, `10^4` Gaussian samples each plus a
trust-region ladder, wall 60.6 s):

| quantity | value |
|---|---|
| `n_inv(X_ini) > 0` (asserted 0) | 0 |
| `n_inv(X0) > 0` | 142 / 200 |
| `n_inv(X0) = 0` | 58 / 200 |
| mean inverted fraction of `X0` | 0.0410 (max 0.1655) |
| `p_valid > 0` (fraction of samples with `n_inv = 0`) | 66 / 200 |
| `p_valid = 0` | 134 / 200 |
| any trust-region radius valid | 69 / 200 |
| `X_ini` injective (control) | 200 / 200 |
| `X0` injective (no face overlap) | 0 / 200 |
| some sample injective | 1 / 200 |
| some sample with `Θ_max > 0` | 0 / 200 |

**Reading.** Measured by inverted faces the projection is invalid on 142/200. Under the
real criterion, no face–face overlap in the flat state, it is invalid on 200/200, and of
the two million samples drawn a single one on a single graph is a valid embedding, with
`Θ_max = 0`. 58 graphs have every face positively oriented and still self-intersect, which
is T5.1's point that positive orientation certifies nothing. On this population the usable
region `U(ε)` is not reachable by isotropic sampling of the null space around `X0`. It is
not empty in general: the authored tilings supply 187 valid deployable configurations
(`deployable_population()`), the population K2a and K1c run on.

Artifacts: `k1a/{k1a.csv,summary.txt}`. Driver `exp_k1a_projection_validity.jl`.

---

## K2c — the broad phase, and the locality hypothesis

**Question.** Is the swept-disc broad phase sound (pruned and unpruned exact `Θ_max`
agree), and is there an `O(n)` active set (H-LOC)?

The originally planned gate — a log-log slope of `max_f ρ_f / r_f` against `log n` — is
vacuous: in the face's own frame `χ_u = −σ_f J x_u`, hence `‖χ_u‖ = ‖x_u‖` and
`ρ_f = r_f` for every face at every size (measured `1.000000000001`). The rule used
instead: pruned and unpruned exact `Θ_max` agree to `1e−12` on every graph, and the growth
of the candidate set with `n` is reported without an `O(n)` claim.

**Measured** (462 graphs, 285 with `F ≤ 800` for the `Θ_max` half; wall 105 s):

| quantity | value |
|---|---|
| pruned == unpruned exact `Θ_max` to `1e−12` | 285 / 285, worst difference 0.000e+00 |
| the same for the flat-centroid variant | 285 / 285, worst difference 0.000e+00 |
| worst \| ‖Y_pv − m_f(θ)‖ − r_f \| (exactness of the moving disc) | 6.729e−13 |
| surviving pairs per face | median 12.69, max 18.40 |
| log-log slope of pairs per face vs `log F` | 0.1695 (R² 0.32, n = 285) |

The exact radius is worth a factor of 1.75 over the loose `√2` radius (median 22.23 → 12.69
pairs per face), `Θ_max` unchanged on every graph.

**H-LOC is refuted.** It would need `max_f ‖γ_f(θ) − x̄_f‖ ≤ κ · r_f` with `κ` independent
of the patch. On a growing square tiling:

| F | patch diameter | `max_f ‖γ_f − x̄_f‖ / r_f` |
|--:|--:|--:|
| 4 | 2.828 | 1.414 |
| 16 | 5.657 | 2.316 |
| 36 | 8.485 | 3.924 |
| 64 | 11.314 | 5.560 |
| 100 | 14.142 | 7.205 |
| 144 | 16.971 | 8.854 |
| 196 | 19.799 | 10.505 |
| 256 | 22.627 | 12.157 |
| 289 | 24.042 | 12.984 |
| 361 | 26.870 | 14.637 |

`drift/r = a + 0.5656 · diameter`, a straight line. On the fixed-box random population the
drift ranges from 25.8 to 105.0 circumradii. The swept region of a face is not contained in
any fixed multiple of its flat circumdisc, so the packing argument does not close. What
survives is an exact broad phase, sound with no hypothesis and empirically tight.

Artifacts: `k2c/{k2c.csv,k2c_drift.csv,summary.txt,k2c_locality.png}`. Driver
`exp_k2c_active_set_locality.jl`.

---

## K2a — exact `Θ_max` versus bisection

**Question.** Does the closed-form `Θ_max` agree with an independent bisection on the
exact overlap predicate to `1e−5` rad on every design?

**Population.** Not random graphs: K1a shows no random graph in this population has an
embedded Eq. (6) projection, and on a non-deployable embedding `deploy()` disagrees with
itself at the shared hinge pins. The population is `deployable_population()`: seven authored
tilings at five clip radii, each with `X0` and up to eight still-embedded shape-space
samples, 187 configurations, every one uniformly deployable and embedded.

**Measured** (187 configurations, wall 40.6 s):

| method | agreement with bisection (4000-point grid) | worst error |
|---|--:|--:|
| T4.2″ interval scan | 187 / 187 | 1.817e−10 |
| min over roots | 182 / 187 | 1.047 rad (`hexagons_R40`) |
| min over roots, `1e−6` shrink | 182 / 187 | 1.047 rad |

Since the F34 predicate fix the referee is tolerance-independent: the min-over-roots row
reads 182/187 at `1e−6`, `1e−9` and `1e−12` alike (before the fix the `1e−6` shrink gave
149/187).

| quantity | value |
|---|---|
| configurations with `Θ_max = 0` (`i* = 0`, immediate penetration) | 14 |
| configurations where min-over-roots < `Θ_max` (a graze) | 5 |
| mean \|C(X)\| / mean intervals probed | 40.1 / 1.9 |
| flat state embedded (overlap predicate at `1e−12`) | 187 / 187 |
| binding edge is a split duplicate / hinge / no contact | 130 / 42 / 10 |
| binding faces not adjacent in `M` | 1 |

The graze correction is worth 5 configurations and, on `hexagons`, `1.047` rad — a third of
the range. The two duplicates of a split edge are congruent translates, so when they become
collinear the contact is vertex-to-vertex: the faces touch at a point and separate. This
holds identically on the whole shape space of a symmetric tiling. On the reference patterns
(`apps/dbg_t422.jl`):

| pattern | \|C\| | `θ₁` (first contact) | `Θ_max` (T4.2″) | bisection | `min β` |
|---|--:|--:|--:|--:|--:|
| rotating_squares | 1 | — | 3.141593 | 3.141593 | 3.141593 |
| triangles_alternating | 0 | — | 3.141593 | 3.141593 | 4.188790 |
| kagome_3636 | 1 | 3.141593 | 3.141593 | 3.141593 | 3.141593 |
| hexagons_auto | 3 | 1.047198 | 2.094395 | 2.094395 | 2.094395 |
| truncated_square_488 | 12 | 2.356194 | 2.356194 | 2.356194 | 2.356194 |
| snub_square_33434 | 18 | 1.646136 | 1.646136 | 1.646136 | 3.383898 |
| tiling_3_4_3_12 | 10 | 2.380636 | 2.380636 | 2.380636 | 2.380636 |

The corrections of the previous section, on this population:

| quantity | value |
|---|---|
| split-free configurations with `Θ_max = min(min β, π)` to `1e−12` | 16 / 16, worst `8.882e−16` |
| of those, patterns where the `π` cap bites | 5 |
| candidate harmonics by class 1 / 2 / 3 + identically zero | 728 783 / 98 376 / 631 + 0 |
| spurious roots in `(0, ε)` removed by the `τ = 0` deflation | 13 630 |
| certificate `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT` holds (`ε = 0.006`) | 173 / 187 (65 / 187 before F32) |
| of which `POS` / `NOOVERLAP` / `NOROOT` alone | 187 / 173 / 187 |
| T5.2b′: certificate holds ⟹ bisection `Θ_max ≥ ε` | 173 / 173, 0 violations |
| configurations that actually have `Θ_max ≥ ε` | 173 / 187 |

Class-2 candidates are 12 % of the population (98 376 of 827 790), and every one has its
`θ = 0` root deflated; without that, 13 630 would be reported as contacts inside `(0, ε)`.
The certificate is sound with no exceptions and, after F32, exact on this population;
`NOOVERLAP` is the binding clause.

One binding pair is not adjacent in `M`: the paper's stated empirical assumption (§4.5)
that the first contact is between faces sharing a split edge is false on at least one of
187 configurations, and true on 130 — a good heuristic, not a theorem.

Artifacts: `k2a/{k2a.csv,summary.txt}`. Drivers `exp_k2a_thetamax_bisection.jl`,
`apps/dbg_t422.jl`, `apps/dbg_k2a_out.jl`.

---

## K1c — Eq. (9)'s false-negative rate

**Question.** Eq. (9) of the 2026 paper certifies `r = h′(0) > 0` at each split edge and
is blind to `p`. How often does a design it accepts have a split-edge re-closure inside its
deployment range?

A root of the separation harmonic is a re-closure only if the duplicates meet as segments
and the faces cross there (T5.3):

```
   r = <d, du> > 0                        Eq. (9) certifies this, and nothing else
   p = -sigma_f det(d, du) < 0            so tau* = -r/p > 0 exists
   theta* = 2 arctan(-r/p) in range
   INTERVAL clause  0 <= dot <= |.|^2 at theta*      (segments, not lines)
   CROSSING clause  the two faces' interiors overlap just after theta* and not just before
```

The T5.2 closed forms were checked against the generic harmonic on every split edge:
`p + q = 0` exactly, `|closed − generic| / scale` worst `1.899e−08` (carried by nearly
collapsed split edges; median `2.0e−15`). Recovering the derivation's `du` from this basis
needs `du = −½ J (S(a″) − S(a′))`.

**Measured** (171 deployable embedded designs with split cuts, wall 61.1 s):

| `Y_9` source | certified designs | `θ* < Θ_max`, no clauses | clauses applied, same range `Θ_max` | corrected: `θ* < min β` | corrected: `θ*` is the binding contact |
|---|--:|--:|--:|--:|--:|
| authors' native `prevent` (published defaults) | 169 | 110 (65.09 %) | 0 (0.00 %) | 16 (9.47 %) | 12 (7.10 %) |
| our `optimize_collision_sweep` (γ ladder) | 163 | 51 (31.29 %) | 0 (0.00 %) | 46 (28.22 %) | 12 (7.36 %) |

Worst `|θ* − per-pair bisection|` over qualifying re-closures: `9.390e−11`.

The native column is a fixed function of the authors' binary; the `optimize_collision_sweep`
column is optimiser-path dependent (an earlier pass ended at 115 certified designs, 40 / 9 /
32 in the last three columns); the conclusion does not depend on which endpoint is used.

**Reading.** Column 4 is zero on both sources and has to be: the crossing clause says the
faces' interiors overlap just after `θ*`, and `Θ_max` is by definition the first angle at
which any two interiors overlap, so "`θ* < Θ_max`" and "`θ*` is a genuine crossing" are
contradictory. The 65.09 % in column 3 is made entirely of roots that are not collisions
(grazes and collinearities). Two well-posed references remain:

* `θ* < min_e β_e` — did Eq. (9) leave a re-closure inside the range the pattern would
  otherwise reach? `min β` is a kinematic bound computed without reference to collisions.
  9.47 % on the authors' code, 28.22 % on ours.
* `θ* == Θ_max` — is a split-edge re-closure what actually stops the deployment? 7.10 % and
  7.36 %.

Eq. (9) is blind to `p` by construction, and on about one design in ten of the authors' own
output that blindness costs range. The 65.09 % figure is not a result.

Our `optimize_collision_sweep` collapses 96 split edges to zero length across this
population; the native `prevent` collapses 0. That is why 163 of 171 of our outputs are valid
flat embeddings against 169 of 171 of theirs, and a reason not to quote our Eq. (9)
reimplementation as a baseline.

Artifacts: `k1c/{k1c.csv,summary.txt}`. Driver `exp_k1c_eq9_false_negatives.jl`.

---

## K5 — is the orientation why the shape space is empty?

**Question.** K1a found that with `σ` from Eq. (1), a max-cut relaxation on the dual graph,
not one of two million shape-space samples over 200 random graphs is a valid embedding.
Eq. (1) maximises the number of hinge cuts and knows nothing about whether the hole-closure
system has a solution near `X_ini`. Does choosing `σ` to minimise the deployability defect
instead fix that?

**Definitions.** `D(σ) = Σ_holes ‖Σ_{hinge e ∈ C} (x_target − x_source)‖²` at `X_ini`,
exactly Eq. (2)'s residual, so `D = 0` iff `X_ini` is already uniformly deployable.
`σ_def`: greedy local search over face flips minimising `D` (single-face and adjacent-pair
flips, strictly improving, any flip rejected that makes `c(Γ) > 1` or detaches a face),
started from `σ_mc` and 3 random `σ`, best kept, at most `20·|F|` attempts per start, seed
`7000 + id`.

**Bar:** `σ_def` gives a certified `X0` (`POS ∧ NOOVERLAP(ε/2) ∧ NOROOT`, `ε = 0.3`) on
≥ 20 % of graphs.

**Measured** (the 200 K1a graphs, `F ∈ [100, 800]`, wall 152 s; `c(Γ) = 1` held on every
accepted `σ`):

| quantity | `σ_mc` (Eq. 1) | `σ_def` |
|---|--:|--:|
| certified `X0`, `ε = 0.3` | 0 / 200 | 0 / 200 |
| of which `POS` alone | 58 | 186 |
| `NOOVERLAP(0.15)` alone | 0 | 0 |
| `NOROOT(0.3)` alone | 0 | 0 |
| certified `Θ_max > 0` | 0 / 200 | 0 / 200 |
| median `D(σ)` | 2.245e+03 | 1.840e+01 |
| median inverted faces of `X0` | 9.5 | 0 |
| median `‖X0 − X_ini‖_∞` / median edge length | 1.587 | 0.334 |
| median `\|E_split\|` | 126.5 | 321.5 |
| median `dim_null` | 126.5 | 321.5 |
| graphs with `D = 0` exactly | 0 | 16 |
| overlap-free `X0` at `θ = 0` (the withdrawn predicate) | 0 / 200 | 125 / 200 |

Correlation of `D` with `|E_split|`: −0.6446 under `σ_mc`, −0.0346 under `σ_def`. Median
flips accepted 138.5.

**Fail: 0 % against 20 %.** Under the `θ = 0` overlap test `σ_def` would score 62.5 %; the
gap between the two predicates is the result.

`σ_def` fixes the flat sheet: defect down 122× in the median, the projection moves the mesh
4.7× less far, inverted faces from 9.5 to zero, `POS` from 58/200 to 186/200, and on 16
graphs `X_ini` is already deployable. It does not touch the deployment: `NOOVERLAP(0.15)` is
false and `Θ_max = 0` on all 400 designs. A log-spaced ladder (`1e−12 … 1.0`) finds an
interior overlap at the first or second rung on every graph that overlaps at all, and a
direct scan (`apps/dbg_k5.jl`) confirms overlap at `1e−2`, `0.1` and `0.5`. The flat sheet
is a valid cut pattern that self-intersects the instant the cuts open. At `θ = 0` the copies
of every split edge coincide, so a `θ = 0` test cannot see this — the reason the certificate
does not use one.

Mechanically, `σ_def` buys its low defect by more than doubling the split cuts (median 126.5
→ 321.5): a split edge imposes no closure constraint, so cutting more edges apart makes
Eq. (2) easy. The hinge graph is then sparse and the structure collapses into itself
immediately. A defect-only objective fails for the same reason max-cut does: neither knows
about collisions. `σ` is a real and badly chosen design variable, but no purely combinatorial
objective on it produces a deployable random graph; the objective has to carry a deployment
term.

The planned follow-on comparison against the range optimiser was not run: both sides would
start from `Θ_max = 0`. `σ_def` orientations are saved as JSON graphs in `k5/sigma/<kind>_<id>.json`
(200 files) and are reused by K6, B4, K9, K9b, K9c and Native200.

Artifacts: `k5/{k5.csv,summary.txt,sigma/,k5_orientation.png}`. Driver
`exp_k5_orientation.jl`, diagnostic `apps/dbg_k5.jl`.

---

## K2b — the range-optimisation margin, against the authors' native code

**Question.** Does `maximize_range` (softmin of the closed-form first-contact roots over
the null space, T6 gradients) beat the authors' own `prevent` (their Eq. 9, run natively)
on deployment range? Bar: median relative gain ≥ 25 % over the best baseline, with
certified validity on ≥ 20 graphs; secondary, snub square closes ≥ 30 % of its gap from
1.646 toward `min β = 3.384`.

**Setup.**
* Baseline: `tuttekiri_cli prevent --sweep` from their Eq. (6) projection, every ladder
  point dumped (`--dumpdir`, added in `baseline/native/src/main.cpp`; nothing in
  `baseline/tuttekiri/` touched), plus their `X0` as the `T = 0` fallback. Full 84-point grid
  `{0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0} × {1, 10, 50, 100} × {0, 0.1, 1.0}` on the four
  reference cases; the 36-point sub-ladder `{0, 0.1, 0.3, 1.0} × {1, 10, 100} × {0, 0.1, 1.0}`
  elsewhere. Their own `theta_max` is never used (F24: `merge_close_verts` fuses hinge
  duplicates before their overlap test).
* Referee: our bisection on the true polygon overlap, 4000-point grid, applied identically
  to every ladder point, our Eq. (9) output and ours. A ladder point counts only if it is a
  valid flat embedding.
* Ours: `maximize_range` with a signed-area log barrier, active set refreshed every 10
  iterations, `T = 0` fallback.
* Population: 38 designs — the 4 reference cases with split cuts, 8 random graphs, 26
  members of `deployable_population()` with `dim_null ≥ 2`. Wall 1533 s.

Two deviations from the plan. The native `prevent` is steep in the null-space dimension
(one ladder point: `0.95 s` at 60 faces, `2.3 s` at 100, `8.8 s` at 150, `27 s` at 200,
`63 s` at 300), so random graphs are capped at 160 faces. And on the random graphs the
experiment is vacuous: 0 of 37 native ladder points is a valid flat embedding (one graph had
4 valid points, all penetrating immediately), so baseline, our Eq. (9) and our optimiser
all return `Θ_max = 0`. Eight such graphs are kept to document it; the other slots went to
`deployable_population()`.

**Measured**:

| quantity | value | bar |
|---|--:|---|
| designs | 38 (30 live, 8 vacuous) | |
| median relative gain over the native baseline | 0.0000 | `≥ 0.25` ✗ |
| median gain restricted to live designs | 0.0000 (n = 30) | ✗ |
| designs with certified validity and range > 0 | 30 / 38 | `≥ 20` ✓ |
| ours better / native better / tie | 7 / 15 / 16 | |
| gain quartiles over the 30 live designs | −0.305, −0.018, 0.000, 0.000, +0.106 | |
| snub square | native `3.1416`, ours `3.1416`, `min β = 3.3839` | gap closed 0.0 %, bar `≥ 30 %` ✗ |

**Fail on both clauses.** The authors' code is strong: on `snub_square_33434` their ladder
reaches `Θ_max = π`, certified by our bisection, so there is no headroom. Our optimiser
mostly does not move — on 16 of 30 live designs the `T = 0` fallback wins. Where it loses,
it loses badly: on `hexagons_R20_s5/_s6` the native reaches `3.0148` and `2.9530`, ours
stays at `2.0944 = 2π/3`, the graze-limited hexagon value, while their penalty walks the
design out of the hexagon geometry. Best gain `+10.6 %` on `tiling_3_4_3_12`. Our Eq. (9)
reimplementation is the weakest of the three (`Θ_max = 0` on every `snub_square_R20*`
design and two `trunc_square` designs, consistent with K1c's 96 collapsed split edges) and
is not cited as a baseline.

The margin claim over the authors' method is therefore not supported by range optimisation
in the null space alone; the diagnostic results (K1a, K2a, K1c and the vacuity finding)
stand.

Artifacts: `k2b/{k2b.csv,summary.txt,k2b_margin.png}`. Driver `exp_k2b_range_margin.jl`.

---

## F23 — is the authors' extra boundary-component row an over-constraint?

**Question.** The authors' code adds a constraint row for split-forest components that
touch the mesh boundary; we drop them as notches. Their `X0` satisfies our system 8/8; our
`X0` has residual 1.73 (`hexagons_auto`) / 1.00 (`snub_square`) in theirs, and our null
space is one dimension larger on those two. The faces around a boundary-touching component
form a path in `Γ`, so no closure constraint should be needed. Test: deploy our `X0` by
forward kinematics from every seed face and 20 randomised BFS orders, at four angles, and
check that all `M′` vertex positions agree to `< 1e−9` and all hinges open by `θ`.

**Measured** (`exp_f23_boundary_row.jl`):

| pattern | F | n_split | H | rank(L) | dim_null | Eq. (2) residual of our `X0` |
|---|--:|--:|--:|--:|--:|--:|
| `hexagons_auto` | 10 | 5 | 4 | 4 | 5 | 2.483e−16 |
| `snub_square_33434` | 53 | 13 | 12 | 12 | 13 | 1.351e−15 |

| pattern | θ | seeds | orders | worst seed disagreement | worst order disagreement | worst hinge-angle error | worst split-duplicate parallel error |
|---|--:|--:|--:|--:|--:|--:|--:|
| hexagons | 0.30 | 10 | 20 | 2.245e−15 | 1.778e−15 | 3.886e−16 | 3.377e−16 |
| hexagons | 0.70 | 10 | 20 | 2.667e−15 | 8.882e−16 | 4.441e−16 | 3.331e−16 |
| hexagons | 1.20 | 10 | 20 | 1.807e−15 | 8.882e−16 | 6.661e−16 | 2.483e−16 |
| hexagons | 2.00 | 10 | 20 | 1.691e−15 | 8.889e−16 | 6.661e−16 | 4.578e−16 |
| snub_square | 0.30 | 53 | 20 | 3.560e−15 | 1.986e−15 | 8.327e−16 | 5.097e−16 |
| snub_square | 0.70 | 53 | 20 | 4.572e−15 | 8.882e−16 | 1.554e−15 | 6.184e−16 |
| snub_square | 1.20 | 53 | 20 | 5.348e−15 | 2.220e−15 | 2.442e−15 | 4.647e−16 |
| snub_square | 2.00 | 53 | 20 | 6.242e−15 | 1.332e−15 | 3.109e−15 | 1.768e−15 |

Worst over both patterns: seed/order disagreement `6.242e−15` (after rigid alignment),
hinge-angle error `3.109e−15`, split-duplicate parallel error `1.768e−15`.

**Yes.** Our `X0` closes every hole and its forward kinematics is path-independent from
every seed face and BFS order, so there is nothing for an extra closure equation to enforce.
Their row is harmless in that their own `X0` still satisfies our system; it only removes
admissible designs, one dimension of shape space on each of these patterns.

Artifacts: `f23/{f23.csv,summary.txt}`. Driver `exp_f23_boundary_row.jl`.

---

## Rank checks on the hole-constraint matrix

From `results/core_validation/rank_claim.md`; final numbers.

1. **`H == |E_hinge| − |F| + c(Γ)`** holds on 50/50 sweep graphs and 8/8 reference cases,
   with `H` counting only all-interior preimages. Counting the notches in breaks it
   (`H_all == pred` on 0/50), so notches are not bounded faces of `Γ`. On a boundary-free
   patch `Γ` lives on a surface of Euler characteristic 0 and the count overshoots by one:
   `H == |E_hinge| − |F| + c(Γ) − 1`, verified on the 4×4 and 6×4 square tori and the 4×4
   triangle torus. Eq. (1) auto-orientation gives `c(Γ) == 1` on 50/50.
2. **`L = R·D`** exactly (max absolute difference 0) on 50/50 sweep graphs and all 8
   reference cases, `D` the `|E_hinge| × N` signed hinge incidence, `R` the `H × |E_hinge|`
   0/1 preimage-ownership matrix. `rank(L) == H − dim Z` on 27/27 (dense only).
3. **Periodic `rank(L) = H − 1`** on the three boundary-free tori: `1ᵀL = 0` exactly,
   `dim Z = 1`, out-harmonic residual `3.3e−16 / 7.2e−16 / 5.6e−16`. So the 2026 paper's
   §4.4 statement that the number of independent equations equals the number of holes is
   false as stated for periodic patterns. With a boundary the naive row-sum identity fails on
   50/50; the corrected form `(1ᵀL)_v = I(v)·indeg(v) − Σ_{v→w} I(w)` holds on 50/50 and 8/8.

Check 2c (every `y ∈ Z` gives an out-harmonic `g`) is vacuous on all 27 dense sweep graphs
and all 8 reference cases (`dim Z == 0`); it has content only on the three tori, where
`dim Z == 1` and the residuals are `≤ 7.3e−16`. `rank_claim.md` records it as
`27 / 27 vacuously`.

---

## K6 — 0⁺ repair inside the shape space

**Question.** K1a and K5 end at `Θ_max = 0` and the exact scan always names the same
culprit: at `θ = 0⁺` the two copies of a split cut translate into each other. The null space
is large (median `dim_null` 126 under `σ_mc`, 321 under `σ_def`). Is there a point in it at
which every cut opens outward?

**The 0⁺ calculus** (`src/method/zero_plus.jl`). By T1.B the two copies of a split edge `e`
differ by a pure translation `sin(θ/2)·dS_e`, so `q_e := det(dS_e, d_e)` decides the side:
`q_e > 0` the cut opens, `q_e ≤ 0` the copies move into each other and `Θ_max = 0`. The
same first-order argument for a corner — a copy of vertex `v` entering the corner cone of a
face at `v` — gives a margin `μ`. Both are exact quadratics in the shape-space coefficients
`t`.

**Bar:** certified `Θ_max > 0` (`ε = 0.3`) on ≥ 20 % of designs for at least one `σ` family.

**Setup.**
* Population: the 200 K1a/K5 graphs (`voronoi` 67, `delaunay` 67, `quad_random` 66; `F`
  from 101 to 793, median 332.5), under `σ_mc` and K5's saved `σ_def` — 400 designs.
* Repair: `F(t) = Σ_split softplus(−q_e/s) + Σ_faces softplus(−a_f/s) [+ w_corner Σ softplus(−μ/s)]
  + w_prox‖X(t) − X0‖²/(N s) + λ‖t‖²`, `s = (median edge)²`, L-BFGS from `t = 0` and 3
  Gaussian starts of scale `0.1·med`, `max_iter = 1200`, seed `6000 + 7·id`. Feasibility is
  always decided by the exact constraint values.
* Primary variant: split + vertex-edge with proximity, `w_corner = 0.05`, `w_prox = 1e3`
  (weights from a 3-graph pilot; `w_corner = 1` lets the corner terms swamp the area terms).
  Secondary: split-only, `w_corner = w_prox = 0`. Also reported: a `λ` ladder
  `{1e−6, 1e−3, 1e−1, 1}` and a stage 2 handing a feasible point to `maximize_range`.
* 12 shards, 1778 s of solver time. The authors' `prevent` from the same `X0` was run on
  designs that are 0⁺-feasible with `F ≤ 550`, at most 12 per family, 420 s timeout, in a
  separate 40-graph capped run (`k6/native/`).

**Measured** (`k6/k6.csv`, 400 rows):

| quantity | `σ_mc` | `σ_def` |
|---|--:|--:|
| designs | 200 | 200 |
| median `\|E_split\|` = median `dim_null` | 126.5 | 321.5 |
| inward split edges at `t = 0` (q1 / median / q3 / max) | 31.8 / 65.5 / 119.8 / 265 | 58 / 103.5 / 192.5 / 340 |
| median inward fraction of split edges | 0.500 | 0.342 |
| designs with zero inward split edges at `t = 0` | 0 / 200 | 0 / 200 |
| inward corner incidences at `t = 0` (q1 / median / q3 / max) | 120.5 / 215.5 / 406.8 / 670 | 149 / 252.5 / 468.3 / 798 |
| designs with zero inward corners | 0 / 200 | 0 / 200 |
| median `min q` at `t = 0` | −77.81 | −130.43 |
| 0⁺-feasible, primary | 0 / 200 | 2 / 200 (1.0 %) |
| 0⁺-feasible, secondary (split-only) | 59 / 200 | 67 / 200 |
| 0⁺-feasible under the `λ` ladder | 0 / 200 | 2 / 200 |
| certified `Θ_max > 0` (`eps_max > 0`), primary | 0 / 200 | 2 / 200 (1.0 %) |
| certified `Θ_max > 0`, secondary | 15 / 200 (7.5 %) | 0 / 200 |
| certified `Θ_max > 0`, ladder / stage 2 | 0 / 200 each | 2 / 200 each (the same two) |
| certificate at `ε = 0.3` | 0 / 200 | 0 / 200 |
| of which `POS` alone | 0 | 2 |
| bisection `Θ_max > 0` on the secondary points | 15 / 200 | 0 / 200 |
| secondary certified `ε_max`, median / max | 2.466 / π | — |
| binding contact at the primary failures | `split-inward` 200 | `split-inward` 188, `vertex-edge` 9, `inverted` 1 |

Which designs the L-BFGS repair lands on is optimiser-path dependent (an earlier pass
reported the primary at 0/400 and the secondary at 0/400 exact, with the bisection positive
on 13 of these `σ_mc` designs); the tallies are those of the canonical `k6.csv`.

*Primary.* The two feasible designs are `voronoi_12` (`F = 521`, 916 split edges) and
`voronoi_93` (`F = 635`, 1124 split edges): `POS`, every `q_e > 0`, every `μ > 0`, and since
F34 a positive exact range (`Θ_max = 0.002354` and `0.239603`, agreed to `2.6e−8` by the
bisection). Neither is certified at `ε = 0.3`. They sit 10.6 and 12.8 median edge lengths
from `X0`: where the repair succeeds, it does so by leaving the projection's neighbourhood.

*Secondary.* The split-only objective reaches a 0⁺-feasible point on 59 `σ_mc` designs, and
on 15 of them — all Delaunay, `F` from 101 to 268 — the certificate is positive with a large
range (`ε_max` median 2.466 rad, max π); the bisection returns the same `Θ_max` on 14 of 15
(on `delaunay_127` the certificate gives 0.058 rad against a bisection of 3.023). With
`w_prox = 0` these points are not near `X0`.

*Capped native comparison* (`k6/native/`, 40 graphs, 80 designs): 8 designs qualified;
ours better / theirs better / tie = 1 / 0 / 3 under `σ_mc`, 0 / 0 / 4 under `σ_def`; median
`Θ_max` 0.0000 for both sides; `native_theta > 0` on 0/80.

*History of the secondary 15.* An early pass reported the exact `Θ_max` positive on 15/200
`σ_mc` designs while the referee agreed on only 13 and disagreed by up to 2.466 rad on three.
Two defects were entangled: `validity_certificate` set `cert.first_root` to the first
admissible root the loop met rather than the smallest (fixed by the guard at
`src/method/contact.jl:526`; with that fix alone the count fell to 0/400), and the pre-F34
overlap predicate rejected the first interval at hinge vertices. With both fixed the exact
scan and the bisection agree on the same 15 designs. The `worst 2.466e+00` line in
`k6/native/summary.txt` is from the 40-graph run with the early referee and affects that
column only.

**Fail: 7.5 % at best against 20 %.** The emptiness of the random-graph design space is not
an artefact of nobody having tried to repair the 0⁺ jam. The jam is first-order and exactly
characterised; the repair problem is solvable on a third of the designs (every `q_e > 0` and
every face area positive on 126 of 400 under the split-only objective), and solving it buys
range on only 17 of 400, all by leaving the projection's neighbourhood, because on the rest a
fresh contact waits behind the one repaired. Fixing `q` moves the binding event from
`split-inward` to `vertex-edge` and then to `inverted`; the obstruction is the whole
first-order cone being empty, not one bad sign. This result ended the random-graph
constructive programme in its original form and moved it to periodic patterns (K7) and then
to constrained embeddings (K9).

Artifacts: `k6/{k6.csv, summary.txt, cache/, native/}`. Driver `exp_k6_zero_plus_repair.jl`,
method `src/method/zero_plus.jl`, diagnostic `apps/dbg_k6.jl`.

---

## A3 — jitter transition

**Question.** The emptiness result reads "0/200 random graphs, 8/8 authored tilings", and
random Voronoi/Delaunay patches differ from authored tilings in every respect at once. Does
deployability die as a sharp function of one geometric ratio when an authored tiling's
interior vertices are jittered continuously, with the combinatorics held fixed? And does a
scalar computable from `(G, σ, X_ini)` predict the outcome (bar: ROC AUC ≥ 0.9)?

**Setup.**
* Driver `exp_a3_jitter.jl`; no `core/` or `method/` source modified.
* Population: the eight `reference_cases()` tilings with `X[v] ← X[v] + a·(median edge)·N(0, I)`
  on interior vertices only, `a` on 20 log-spaced steps in `[0.005, 1.0]`, 20 seeds each,
  plus one `a = 0` baseline per tiling: 3208 designs, each under two σ rules = 6416 rows.
  Boundary vertices are not jittered because the pipeline solves Eq. (4)/(6) with
  `BoundaryMode::Fixed`.
* σ rules as in K5: `σ_mc` depends on the dual graph only and is constant along the ladder;
  `σ_def` is K5's greedy search re-run on each jittered geometry (cap `3|F|` attempts).
* Measurement: `make_cut → holes_partition → assemble_system(Fixed) → solve_system → X0`,
  then the certificate at `ε = 0.006` and `ε = 0.3`, and the exact `Θ_max` by T4.2″.
* 12 shards, wall 2.2 s.
* Predictors from `(G, σ, X_ini)`: `rho_geom` (median over split-forest components of
  geometric diameter / median face inradius), `rho_hop`, `rho_max`, `psplit = |E_split|/|F|`,
  `badq_ini` (fraction of split edges with `q_e ≤ 0`), `aniso`, `D(σ)`, `amp` as positive
  control; `X0`-derived controls `badq_X0`, `margin_X0` reported separately.

**Measured** (`jitter/{jitter.csv, ladder.csv, summary.txt, transition.png}`).

The population splits in two by `|E_split|`. Four of the eight reference cases are
split-free under their own σ (`rotating_squares`, `triangles_alternating`, `kagome_3636`,
`periodic_squares_4x4`): `dim_null = 0`, the fixed-boundary Eq. (4) system has a unique
solution, and the Eq. (6) projection sends any jittered interior geometry back to the
authored tiling exactly (`‖X0 − X_base‖_∞ = 0.0000` at every amplitude up to `a = 1.0`,
while `‖X0 − X_ini‖_∞` reaches 2.8 median edges). Those four curves are flat at 1.0 because
they never left the start point (unit test "split-free tiling: Eq. (6) projection undoes an
arbitrary interior jitter"). The ladder is informative only on the four split-bearing
tilings (`hexagons_auto` 5, `truncated_square_488` 12, `snub_square_33434` 13,
`tiling_3_4_3_12` 4 split cuts under `σ_mc`).

On those four the exact `Θ_max` has a monotone, reasonably sharp transition:

| tiling (σ_mc) | `a*` (Θ_max > 0) | `a` at frac ≥ 0.9 | `a` at frac ≤ 0.1 | width |
|---|--:|--:|--:|--:|
| `snub_square_33434` | 0.163 | 0.107 | 0.248 | 0.36 dec |
| `hexagons_auto` | 0.277 | 0.142 | 0.573 | 0.61 dec |
| `truncated_square_488` | 0.320 | 0.248 | 0.573 | 0.36 dec |
| `tiling_3_4_3_12` | 0.522 | 0.248 | 0.757 | 0.48 dec |

`a*` under `σ_def` is uniformly smaller (0.071 / 0.188 / 0.096 / already dead at 0.005),
and `σ_def` doubles to quintuples `|E_split|` (`snub` 13 → 25, `3_4_3_12` 4 → 20). Over all
split-bearing rows `Θ_max > 0` holds on 1225/1604 (76.4 %) under `σ_mc` and 634/1604
(39.5 %) under `σ_def`.

**The certificate, and F32.** As first run, at `a = 0.005` every split-bearing tiling had
`POS 20/20`, `NOOVERLAP 20/20` and exact `Θ_max` between 1.0 and 2.4 rad on 20/20 seeds, yet
the certificate passed on 0/20 (`hexagons_auto`, `truncated_square_488`), 5/20
(`snub_square`), 10/20 (`tiling_3_4_3_12`); over the ladder 1622 of the 1859 rows with
`Θ_max > 0` were not certified (87.3 % false negatives). Diagnosis
(`jitter/cert_diagnosis.md`): an infinitesimal perturbation of an exactly symmetric tiling
splits a degenerate contact and pushes a collinearity root — the duplicate edge on the *line*
of its partner, not on the segment — into `(0, 0.006)`, and `NOROOT` counted it because it
lacked the two interval tests. With the interval tests restored (F32, `contact.jl`) the
canonical `jitter.csv` reads: certificate holds on 20/20 seeds of all four split-bearing
tilings at `a = 0.005`, the `a = 0` baselines certify 8/8, and 1856 of the 1859 rows with
`Θ_max > 0` are certified; the three exceptions have `Θ_max` between 0.0018 and 0.0060 rad,
below `ε = 0.006`, and are correct rejections. `a*(certified)` coincides with
`a*(Θ_max > 0)` on every tiling.

**Predictors** (directed AUC = max(AUC, 1−AUC); the split-bearing rows are the right
population, the all-rows column is contaminated by the split-free dichotomy):

| predictor | AUC(Θ_max > 0), all rows | AUC(Θ_max > 0), split-bearing | AUC(certified), split-bearing |
|---|--:|--:|--:|
| `rho_geom` | 0.918 | 0.778 | 0.777 |
| `rho_hop` | 0.892 | 0.706 | 0.706 |
| `rho_max` | 0.952 | 0.870 | 0.870 |
| `psplit` | 0.874 | 0.657 | 0.656 |
| `badq_ini` | 0.987 | 0.976 | 0.976 |
| `aniso` | 0.745 | 0.816 | 0.816 |
| `D(σ)` | 0.585 | 0.504 | 0.504 |
| `amp` (positive control) | 0.743 | 0.831 | 0.831 |
| `margin_X0` (X0 control) | 1.000 | 1.000 | 1.000 |
| `badq_X0` (X0 control) | 1.000 | 1.000 | 0.999 |

The AUC routine was cross-checked against a brute-force Mann–Whitney count on `badq_ini`
(0.012789 vs 0.0128).

**Reading.** The transition exists, is monotone and sharp (0.36–0.61 decades), with
`a* ∈ [0.16, 0.52]` under `σ_mc`: deployability survives perturbations of about one sixth
of an edge length and dies around one third to one half. Emptiness is a function of measured
geometry, reached continuously from patterns whose combinatorics were never randomised. The
predictor bar is met only by `badq_ini` (0.976), which restates the K5/K6 failure mode rather
than giving an independent geometric characteristic; the proposed ratio `rho_geom` reaches
0.778 and `rho_max` 0.870. Two consequences: split-free authored patterns are outside the
problem (their shape space is a point), so "authored tilings pass 8/8" is "4 trivially, 4
substantively"; and the certificate needed the interval tests before it could label authored
geometry.

Artifacts: `jitter/{jitter.csv, ladder.csv, summary.txt, transition.png, cert_diagnosis.md}`.
Driver `exp_a3_jitter.jl`, plot `scripts/plot_jitter.jl`.

---

## K7 — periodic deployment Jacobians

**Question.** On a torus, is the deployment of a periodic pattern an affine family? With
`P_θ = [p_x(θ) p_y(θ)]` the deployed period matrix and `J(θ) = P_θ P_0⁻¹`:

* C1: `J(θ) = cos(θ/2) I + sin(θ/2) K(X)` exactly, `K` constant in `θ` and affine in
  `X = X_0 + Φ t`, so the achievable set `𝒦 = {K(X)}` is an affine subspace and inverse
  design of the whole `θ`-history is one linear solve. Bar: all tolerances `1e-9` on ≥ 30
  patterns.
* C2: the conformal designs are the solutions of two linear equations in `t`
  (`K₁₁ = K₂₂`, `K₁₂ = −K₂₁`), and such a design is conformal at every `θ` — the sentence
  the 2026 paper's §5.1 asserts empirically. Bar: `< 1e-9` at 200 angles.
* C3: a target `K*` can be hit by least squares with leftover freedom spent on certified
  `Θ_max`. Bar: hit `< 1e-8` with a certified `Θ_max > 0` design on ≥ 80 % of patterns
  with `dim 𝒦 ≥ 1`.
* C4: the per-cell hole area is a first harmonic, so the second fully-closed angle has a
  closed form. Bar: `|θ_c − measured| < 1e-6` on every pattern that has one.

**Setup.**
* Driver `exp_k7_periodic_jacobian.jl` (stages `main`, `c3`, `c4bounded`, `nu`, `c4sweep`,
  `perdbg`); periodic layer `src/method/periodic_jacobian.jl` (lattice detection, torus
  quotient, super patch, `periodic_jacobian`, achievable set).
* Why a quotient: the finite-patch route (`Kirigami/README.md` judgement call 11) imposes
  Eqs. (3b)–(3c) on a patch whose boundary edges are not identified, so every hole preimage
  straddling the seam contributes no row of `L` and the wrap-around conditions go missing.
  K7 builds the genuine quotient: one face per translation class, edges keyed by (class pair,
  lattice offset), hole rows carrying the integer offset sum. The finite-patch `dim_null` is
  reported as `dim_null_patch` and is systematically larger (`voronoi_torus_3_n45`: 45 on
  the quotient, 63 on the patch).
* Population: 33 patterns — seven periodic families (`squares`, `triangles`, `hexagons`,
  `kagome`, `snub_square`, `trunc_square_488`, `t3_4_3_12`) at 2×2, 3×2, 3×3 cells, plus 12
  periodic Voronoi tori of 20–200 faces per cell; σ from Eq. (1) solved as max-cut on the
  quotient dual.
* `K` comes from the closed-form deploy basis; every check uses an independent `deploy()`
  forward-kinematics call. Certification is the standard certificate with the `first_root`
  fix. 12 shards; the `main` stage is set by `voronoi_torus_11_n200` (1824 s).

**Measured** (`k7/{k7_main.csv (33), k7_c3_all.csv (50), k7_c4_bounded.csv (7), k7_periodicity.csv (33), k7_nu_*.csv, k7_c4_sweep.csv}`).

**C1 — pass.**

| quantity | bar | worst over 33 patterns |
|---|--:|--:|
| `‖P_θ^fit − P_θ^FK‖ / ‖P_θ^FK‖`, 50 angles | `1e-9` | 7.29e−14 |
| `‖K_fit − K_basis‖_∞` | `1e-9` | 5.99e−14 |
| second difference of `t ↦ K` (affinity), 31 triples/pattern | `1e-9` | 4.79e−14 |
| `‖P_0 − T‖_∞` (measured period vs the lattice) | — | 6.39e−14 |

`dim 𝒦 = 2·rank(D)` on 33/33, `D` the `k × 2` matrix of period-potential increments. The
guess `min(4, 2·dim_null)` fails on 1/33: `squares_3x3` has `dim_null = 4` (6 split cuts)
yet `dim 𝒦 = 0` — its `K` is frozen at a rotation (`tr K = 0`, `det K = 1`) over the whole
shape space, so the cell never opens. Shape-space dimension does not bound designability;
the rank of the period-potential map does. Where `dim 𝒦 = 4` (25 of 33) the achievable set
is all of the 2×2 matrices. Figure `k7_dimK_hist.png`.

**C2 — pass, and it settles §5.1.** On all 25 patterns with `dim 𝒦 ≥ 1` the two
conformality equations are consistent (residual ≤ 9.19e−16), and at the solution the
conformal distortion of `J(θ)` over 200 angles in `[0, π]` is ≤ 2.24e−14. Since `J(θ)` is
an affine path through `I` and the similarities are the span of `{I, Jrot}`, `J(θ)` is a
similarity for every `θ` exactly when `K` is one; the paper's empirical statement is a
one-line consequence of C1.

**C3 — fail, 52 % against 80 %.**

| target | patterns | max hit error | `0⁺` feasible | certified `Θ_max > 0` |
|---|--:|--:|--:|--:|
| random `K* ∈ 𝒦` | 25 | 1.44e−15 | 21/25 | 13/25 (52 %) |
| anisotropic `K* = diag(1, −0.5)` | 25 | 5.97e−15 (24) / 1.0 (1) | 0/25 | 0/25 |

The certified `ε` agrees with the exact `Θ_max` to `2.2e−9` rad on all 13. The random
target is drawn from a null-space basis, so which target is drawn and which design realises
it depends on the linear-algebra backend (earlier passes gave 11/25 and 12/25); the
anisotropic row is fixed. The `c3_cert` columns of `k7_main.csv` are from the `main` stage's
own draw; read `k7_c3_all.csv`.

The linear solve does everything C1 promises: on 24 of 25 patterns even the anisotropic
target is hit to `8.9e−16`. The one unreachable case is `squares_3x2`, the only pattern with
`rank(D) = 1` (hit error exactly `1.0`). Every other failure is the same 0⁺ obstruction as
K5/K6: the design realising `K*` places at least one split cut so that its duplicates
translate into each other, `min_e q_e < 0`, `Θ_max = 0`.

| family | `nsplit` per cell | `dim 𝒦` | `rank D` | certified (random target) | `min q` at the anisotropic design |
|---|--:|--:|--:|--:|--:|
| `squares` 2×2, `triangles`, `kagome` (8 patterns) | 0 | 0 | 0 | — (`𝒦` is a point) | — |
| `squares_3x3` | 6 | 0 | 0 | — (`K` is a rotation) | — |
| `squares_3x2` | 2 | 2 | 1 | 0/1 | −33.5, target unreachable |
| `hexagons` | 4–9 | 4 | 2 | 2/3 | −5.91 … −0.30 |
| `snub_square` | 8–18 | 4 | 2 | 2/3 | −41.5 … −14.9 |
| `trunc_square_488` | 8–18 | 4 | 2 | 2/3 | −2.66 … −0.35 |
| `t3_4_3_12` | 8–24 | 4 | 2 | 3/3 | −49.0 … −38.7 |
| `voronoi_torus` | 20–204 | 4 | 2 | 4/12 | −15.0 … −0.29 |

`kagome`, `triangles` and `squares_2x2` are split-free under the max-cut σ on the quotient
dual, so their shape space is a point — the same dichotomy A3 found on patches.
`squares_3x3` is the opposite degeneracy. On the remaining 25, `dim 𝒦 = 4` means
reachability is not the question; validity is. On `snub_square` and `t3_4_3_12` the
anisotropic design overshoots by `min q ≈ −39` to `−49`; on `hexagons`, 4.8.8 and Voronoi
the overshoot is `O(1)` and sometimes only `−0.29` (`voronoi_torus_0_n20`), where a
margin-aware search might recover a design.

Where a design is certified, the closed-form `ν(θ)` is right: over the 13 certified designs
the largest discrepancy between `ν` predicted from `K` and forward kinematics at 20 angles
and 2 directions is `1.15e−09`; both the certified `ε` and the exact `Θ_max` reach 2.18 rad
(`t3_4_3_12_2x2`). Figures `k7_nu_curves.png`, `k7_c3_certified.png`.

**C4 — pass.** `A(θ) = det(P_0)(det J(θ) − 1)`; substituting C1 and the half-angle
identities gives the first harmonic `A(θ) = q(cos θ − 1) + r sin θ` with

```
q = -det(P_0) (det K - 1) / 2 ,    r = det(P_0) tr K / 2 ,
theta_c = 2 atan2(r, q) = 2 atan2(tr K, 1 - det K) .
```

`det P_0` cancels, so the right-hand form holds identically; `p + q = 0` automatically, so
`θ = 0` is always a root and the second root is the content.

| population | rows | worst `|θ_c − θ_measured|` | bar |
|---|--:|--:|--:|
| periodic (bisection on FK hole area) | 32 (of 33) | 5.62e−14 | 1e−6 |
| bounded patches, `k7_c4_bounded.csv` | 7 | 4.44e−16 | 1e−6 |

The excluded row is `squares_3x3` (`q = r = 0`, `A(θ) ≡ 0`, no second closed state). On the
bounded patches the closed state is overlap-free (`θ_c ≤ Θ_max`) on 4/7 — `rotating_squares`
(π), `kagome_3636` (π), `hexagons_auto` (2π/3), `truncated_square_488` (3π/4) — and
unreachable on `triangles_alternating` (`θ_c = 4π/3 > π`), `snub_square_33434` (`3.768` vs
`Θ_max = 1.646`) and `tiling_3_4_3_12` (`2.818` vs `2.381`).

Two anomalies. The area-harmonic fit error over `θ ∈ (0, 2π)` reaches `2.3e−2` on the
Voronoi tori: `k7_c4_sweep.csv` / `k7_c4_area_fit.png` localise it to the interval where
`det P_θ < 0` (the cell has inverted), an artefact of measuring `|det P_θ|` where the
harmonic is the signed quantity; off that interval the agreement is `1e−15`. And
`periodic_jacobian`'s consistency statistic is `≤ 1.1e−13` on 32/33 patterns and 27.5 on
`snub_square_3x3`: its super patch's cut structure has 9 connected components, each threading
all 9 cell copies, and `deploy()` gauges each independently, so the period read off one
component differs from another by a gauge offset growing with `θ`. `snub_square_3x3`'s `K` is
gauge-dependent and should not be quoted; removing it changes none of the worst-case figures
above.

**Verdict.**

| claim | bar | measured | verdict |
|---|---|---|---|
| C1 affine Jacobian | `1e-9`, ≥ 30 patterns | 7.29e−14 / 4.79e−14 on 33 | pass |
| C2 exact conformality | `< 1e-9` at 200 angles | 2.24e−14, 25/25 | pass |
| C3 designable Poisson family | ≥ 80 % certified | 13/25 = 52 % (basis-dependent draw) | fail |
| C4 closed-form `θ_c` | `< 1e-6` | 5.62e−14 on 32 + 7 | pass |

The algebra is right and better than claimed: `J(θ)` is exactly affine, the achievable set
is on 25 of 33 patterns the whole 2×2 matrix space, conformal design is two linear equations,
and the second closed angle is `2 atan2(tr K, 1 − det K)`. The geometry fails in the same way
as K5, K6 and F30: the solution set of the linear design problem is dominated by embeddings
whose split cuts open inward at `θ = 0⁺`. That K7 reproduces the 0⁺ obstruction on a
population sharing none of K5/K6's construction (periodic, tori, tilings as well as random
Voronoi) is the strongest evidence that it is a property of the problem and not of one
generator.

Unit test: "C4: theta_c = 2 atan2(tr K, 1 − det K) is where forward kinematics recloses"
(`test_method.jl`), on `hexagons_2x2` and `t3_4_3_12_2x2` (recorded values 2.09439510239 and
2.72771226181).

Artifacts as listed above; figures `k7_nu_curves.png`, `k7_dimK_hist.png`,
`k7_c3_certified.png`, `k7_c4_area_fit.png`. Tiling seeds `20260904 + 7919·i`, Voronoi seeds
`9100001 + 104729·i`; a single pattern `i` re-runs with `--nshard 33 --shard i`.

---

## B4 — does the fixed boundary empty the space?

**Question.** K6 measured 2/400 certified under its primary repair on fixed-boundary
patches; K7 ran the same Voronoi generator on tori and got 11/12 zero-plus feasible and 4/12
certified. The two populations differ in one thing: the patch carries Eq. (4)'s fixed
boundary rows `B`, the torus does not. Is `B`, rather than graph randomness, what empties
the space? Bar: exact `Θ_max > 0` on ≥ 10 % of the 400 `free` designs.

**Three boundary variants**, all on K6's population verbatim (same graphs, `σ_mc`, saved
`σ_def`, `X_ini`, repair objective, weights `w_corner = 0.05`, `w_prox = 1e3`, `λ = 1e−6`,
`max_iter = 1200`, 3 Gaussian starts, seeds `6000 + 7·id + which`):

| variant | system solved | note |
|---|---|---|
| `fixed` | `[L; B] X = [0; T]` | the K6 system, run as a control from K6's own cache |
| `free` | `[L; e_0ᵀ] X = [0; x_0]` | `B` dropped. The only gauge `L` cannot see is global translation (every `L` row sums to zero), killed by pinning one vertex. Rotation and scaling are left in: `q_e` is rotation-invariant and positively homogeneous, so they change no label |
| `split_free` | `B` restricted to boundary vertices touching no split edge | free exactly the border where the cuts live |

Downstream is K6's pipeline unchanged: Eq. (6) projection, `zero_plus_repair` (split-only
first, then primary warm-started from it), the exact `Θ_max` of T4.2″ as primary label, the
certificate as secondary, the bisection as referee. 12 shards.

**Measured** (`b4/b4.csv`, 1200 rows = 200 graphs × 2 σ × 3 variants):

| variant \| σ | designs | median `dim_null` | median `min q` at `t=0` | 0⁺-feasible (split-only) | 0⁺-feasible (primary) | exact `Θ_max > 0` | certified `eps_max > 0` |
|---|--:|--:|--:|--:|--:|--:|--:|
| `fixed` \| `σ_mc` | 200 | 126.5 | −77.81 | 59 | 0 | 0 | 0 |
| `split_free` \| `σ_mc` | 200 | 131 | −66.06 | 85 | 0 | 0 | 0 |
| `free` \| `σ_mc` | 200 | 157 | −51.49 | 176 | 0 | 0 | 0 |
| `fixed` \| `σ_def` | 200 | 321.5 | −130.43 | 67 | 2 | 2 | 2 |
| `split_free` \| `σ_def` | 200 | 330 | −130.43 | 180 | 1 | 1 | 1 |
| `free` \| `σ_def` | 200 | 339.5 | −130.43 | 200 | 2 | 2 | 2 |

Which one or two designs the primary repair lands on is a knife-edge outcome of an
unconverged 1200-iteration L-BFGS run (an earlier pass had 3 feasible and 2 positive under
`free`, with a different pair); the `t = 0` columns, the split-only counts and the dimension
identity below are not path-dependent.

**Fail: 2/400 = 0.50 % under `free`, the same rate as the `fixed` control, against 10 %.**

Three things the run establishes:

1. The control reproduces K6 exactly: under `fixed` the primary repair is feasible on
   `voronoi_12` and `voronoi_93` under `σ_def`, the same two as K6, with `Θ_max` 0.002354 and
   0.238891 (K6: 0.002354 and 0.239603; the drivers share `X0` and differ only by the
   repair's floating-point path).
2. The dimension count is exactly right: `dim_null(free) − dim_null(fixed) = |V_∂| − 1` on
   400/400 rows (mean gain 32.8; the `−1` is the pin). `split_free` buys a mean of 12.4.
3. The boundary does constrain the 0⁺ cone, and that buys almost no range. The split-only
   repair goes from 59/200 to 176/200 under `σ_mc` and from 67/200 to 200/200 under `σ_def`
   when `B` is dropped, monotonically through `split_free`; the median `min q` at `t = 0`
   improves under `σ_mc` (−77.8 → −51.5). None of it survives the vertex-edge margin: the
   primary objective is feasible on 2 of 200 at best, and the binding contact migrates
   `split-inward` 185, `vertex-edge` 12, `inverted` 1 under `free`/`σ_def`. The first-order
   cone is emptied by the number of simultaneous sign conditions, and `|V_∂| ≈ 33` extra
   dimensions against a median of 126 to 321 split-sign constraints is not enough.

The positives, all refereed:

| design | variant | `F` | `n_split` | `dim_null` | `min q` | exact `Θ_max` | certified `eps_max` | bisection |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| `voronoi_96`, `σ_def` | `free` | 497 | 866 | 947 | 0.1640 | 0.24882 | 0.24882 | 0.24882 |
| `voronoi_12`, `σ_def` | `free` | 521 | 916 | 995 | 0.1058 | 0.07594 | 0.07594 | 0.07594 |
| `voronoi_12`, `σ_def` | `split_free` | 521 | 916 | 960 | 0.0302 | 0.02128 | 0.02128 | 0.02128 |
| `voronoi_93`, `σ_def` | `fixed` | 635 | 1124 | 1124 | 0.1957 | 0.23889 | 0.23889 | 0.23889 |
| `voronoi_12`, `σ_def` | `fixed` | 521 | 916 | 916 | 0.0035 | 0.00235 | 0.00235 | 0.00235 |

`|Θ_exact − bisection| ≤ 1e−5` on 5/5, worst gap `2.0e−8`. The `free` positives sit 10.3 and
10.7 median edge lengths from `X0`: where the free boundary helps, it helps by leaving the
projection's neighbourhood, with nothing holding the border in place.

The sentence "the published boundary condition empties the space" is not licensed: the
space is essentially empty either way, and the torus advantage is not explained by the
absence of `B`.

Artifacts: `b4/{b4.csv, shard_*.csv, summary_shard_*.txt, cache/}`. Driver
`exp_b4_boundary_vs_graph.jl` (builds its variant systems from
`assemble_system(..., BoundaryMode::None)` and appends the pin rows itself).

---

## K8a — the expansive-cone LP

**Question.** F30 says the *uniform* deployment collides at `0⁺` on random graphs. Is
there any first-order flex — uniform or not — that separates every cut and every corner?
This is the polyhedral expansion cone of Rote–Santos–Streinu 2003 (Lemma 3.2) and
Connelly–Demaine–Rote 2003, instantiated on a body-and-pin framework with a cut structure,
where the rows are per split edge and per corner incidence. Bar: a strictly feasible flex
on ≥ 20 of 100 K1a graphs; hard soundness control on four split-bearing tilings whose
uniform ray must be inside the cone; Euler-step sanity on ≥ 10.

**Setup.** A flex gives every face an angular velocity and a translation,
`V_f(y) = ω_f J y + w_f`, and hinge edges impose `V_f(p_i) = V_g(p_i)` at the pin
(`mobility.jl::build_rigidity`). With `X` held fixed every 0⁺ separation condition of
`zero_plus.jl` becomes a strict linear inequality in the flex:

```
split edge e :  row_e = det( V_{f1}(X_v) - V_{f0}(X_v), d_e )   = q_e / 2
corner (p into the corner of f at v) :  -g1 = det(dV, e1),  -g2 = -det(dV, e2),
                                        mu = max(-g1,-g2) convex, min reflex,  = mu_0+/2
```

Both identities are asserted to `1e-9` relative against `zero_plus_q` and
`zero_plus_corner_margin` in a unit test. The flex basis is `ker A` plus the integrated
translations plus one free translation pair per component of `Γ`; a unit test checks
`dim ker R = dim ker A + 2 c(Γ)` against an independent dense QR rank. Over all 277
configurations, `max |R N| / max|R| = 1.2e-13`.

With rows normalised the quantity solved is the matrix-game value
`val = max_{‖z‖₂ ≤ 1} min_i a_i · z = min_{λ ∈ simplex} ‖Aᵀ λ‖₂`. Both sides are `≥ 0` and
`val > 0` iff the strict system is feasible. Every iterate of either loop is feasible for
its side, so the reported `margin_l2` is a rigorous lower bound and `dual_bound` a rigorous
upper bound on the LP value; a small dual bound is an approximate Farkas certificate. Three
hand-solved LPs pin the solver in a unit test. The `l2` game is used rather than the
`l_∞`/multiplicative-weights form originally proposed; the cone is homogeneous so the two
are positive iff each other, and `margin_inf` is reported anyway.

Every one of the 244 931 corner incidences on the K1a population at `X_ini` is at a convex
corner, so `P(X)` is a union of `2^{n_convex}` polyhedra. Four passes are run: the
conjunctive relaxation (a subset of `P(X)`), the branch the uniform ray selects, and two
repairs from the best witness. Feasibility of any pass proves `P(X) ≠ ∅`; infeasibility of
the passes tried proves nothing — the negative reading rests on the dual bound.

*Reproducibility.* The combinatorial columns of `k8a.csv` reproduce bit for bit. The LP
values are the endpoints of an iterative solver (away-step Frank–Wolfe, 100 000 dual steps)
and move at the last digit or two between environments; the σ-chart bound on a
non-deployable `X_ini` and the margin magnitudes on the truncated-square family depend on a
degenerate-row normalisation and are not reproducible across environments.

**Hard soundness control: 4 of 4.**

| tiling | F | split | dim ker A | dim flex | σ ∈ P(X) | min `q_e` | min `mu` | LP margin |
|---|---|---|---|---|---|---|---|---|
| `hexagons_auto` | 10 | 5 | 3 | 5 | yes | 1.000 | 1.000 | 6.21e−1 |
| `truncated_square_488` | 21 | 12 | 2 | 4 | yes | 1.414 | 1.414 | 9.25e−1 |
| `snub_square_33434` | 53 | 13 | 29 | 31 | yes | 0.688 | 0.509 | 6.58e−2 |
| `tiling_3_4_3_12` | 25 | 4 | 5 | 7 | yes | 1.742 | 1.216 | 3.12e−1 |

(`tiling_3_4_3_12` at its Eq. (6) projection, the others at `X_ini`: "X_ini if Eq. (2)
holds there, else the projection".) Extended control on `deployable_population()` (73
configurations): σ ∈ `P(X)` on 71 of 73, LP margin > 0 on 72 of 73.

**The K1a population.**

| | at `X_ini` | at `X0` |
|---|---|---|
| LP margin > 0 | 0 / 100 | 0 / 100 |
| `dual_max < 1e-6` (every branch chart tried certified) | 95 / 100 | 99 / 100 |
| `dual_max < 1e-5` | 99 / 100 | 100 / 100 |
| `dual_max` median / max | 6.7e−10 / 1.3e−5 | 6.7e−10 / 1.7e−6 |
| σ is a flex of the framework | 0 / 100 | 100 / 100 |
| σ ∈ `P(X)` | 0 / 100 | 0 / 100 |
| dim ker A, median | 17 | 65 |
| rows, median | 4 911 (142 split, 2 338 corner incidences) | same |

Per family at `X_ini`: Voronoi `dim ker A` median 3 (34/34 dual-certified), quad-random 17
(33/33), Delaunay 72 (28/33, `dual_max` median 4.9e−9). The Delaunay residue is Frank–Wolfe
convergence (at 20 000 dual steps the same graphs sat at `1e−4`, at 100 000 at `1e−6` or
below); no Delaunay graph has a positive primal margin at any budget.

**Reading.** On this population the emptiness is not "one ray of a 305-dimensional cone
points the wrong way". On 95 of 100 graphs there is a non-negative combination of the
split-edge and corner rows summing to within `1e-6` of zero — an approximate Farkas
certificate that no flex at all separates every cut and corner at first order, over every
branch chart reached. The failure is a property of the cut structure, not of the two
orientation heuristics. `σ ∈ ker A` on 100/100 at `X0` and 0/100 at `X_ini` reproduces K3a's
split from independent code; at `X0`, σ is a flex and is nevertheless outside `P(X0)` on
100/100 (median `min q_e = −79.2`).

One non-uniform rescue exists: `snub_square_R20_s1` has σ outside `P(X)`
(`min q_e = −6.05e-2`, `Θ_max = 0`) and the LP certifies a strictly positive non-uniform
margin of `5.93e-2`, with a collision-free Euler step at the largest tested amplitude.
`trunc_square_R40_s0` is the control: σ outside `P(X)` and the LP finds nothing.

**Euler-step sanity.** On every configuration with a positive margin the witness flex is
integrated one explicit Euler step at `h = eps · (median edge) / max_p |V_p|`, `eps`
descending through `{0.5, 0.2, 0.1, 0.03, 0.01, 0.003, 0.001}`, and tested for face–face
overlap: 72 of 72 collision-free, 65 already at `eps = 0.5`, 4 at 0.2, 2 at 0.1, 1 at 0.03.
The T4.2″ scan is parameterised by the uniform angle and has no meaning along a non-uniform
flex, so the predicate is `has_collision` / `polygons_overlap`; a unit test checks that along
the uniform ray the Euler step agrees with `deploy(c, X, h)` to `1e-3` median edges.

**Caveats.** Infeasibility is certified to `1e-6`, not proved; 5 of 100 graphs at `X_ini`
have `dual_max ≥ 1e-6` (all Delaunay, below `1.3e-5`), and there the negative rests on "no
witness found". The branch enumeration (four of `2^{2338}` charts) is not exhaustive and
cannot be. `dim ker A` medians here (Voronoi 3, quad-random 17, Delaunay 72 at `X_ini`) are
the full graph at a different size distribution, not K3a's 2-core median of 305. The margin
is scale-dependent through the row normalisation; only its sign is invariant.

### Recheck — the solver was wrong; the certificates were not

Full write-up `k8a/recheck.md`; artifacts `k8a/recheck_k9.csv`, `recheck_k9_summary.txt`,
`k8a/recheck/` (20 000 dual steps) and `k8a/recheck100k/`. Driver
`exp_k8a_expansive_cone_recheck.jl`.

`cone_lp` computed the two ends of its bracket with unrelated algorithms. The dual end
(Frank–Wolfe on `min_λ ‖Aᵀ λ‖`) was sound. The primal end came from a separate log-sum-exp
smoothing loop (about 1 000 projected-gradient steps of size `≈ 1e−5`, from `z = 0`, in up
to 300 dimensions against 5 000 rows) and was the only source of the witness. The
Frank–Wolfe iterate `w = Aᵀ λ` was never evaluated as a primal direction, although it is the
optimum: at the minimum-norm point of `conv{a_i}` one has `a_i · w ≥ ‖w‖²` for every `i`,
so `z = w/‖w‖` attains `min_i a_i · z ≥ ‖w‖` and the bracket closes. Fix: away-step
Frank–Wolfe (Lacoste-Julien & Jaggi 2015) on the min-norm point with the primal witness read
off the same iterate; new `ConeLPResult::gap`, `ExpansiveConeReport::sigma_chart_margin` (a
solver-free lower bound), and `farkas_residual()`.

On K9's 30 variant-(b) embeddings, regenerated bit-identically:

| | old solver | corrected |
|---|---|---|
| `σ ∈ P(X)` by direct `zero_plus` measurement | 23 | 23 |
| of which LP margin `> 0` | 0 / 23 | 23 / 23 |
| of which margin ≥ the solver-free `sigma_chart_margin` | — | 23 / 23 |
| of which a non-uniform flex strictly beats `σ` | — | 23 / 23 |
| `σ` outside `P(X)`, LP margin `> 0` anyway | 0 / 7 | 6 / 7 |
| feasible overall | 0 / 30 | 29 / 30 |

`margin / sigma_chart_margin`: median 8.73, max 56.9. The `feasible 0` lines of the original
K9 pass were solver artefacts.

The certificates were real: every dual multiplier vector is a genuine simplex point (worst
`min_i λ_i` exactly `0`, worst `|Σλ − 1|` `2.5e−13` over 277 configurations). Independently
verified, 100 000 dual steps:

| | at `X_ini` | at `X0` |
|---|---|---|
| reported `dual_max < 1e−6` (original run) | 84 / 100 | 68 / 100 |
| verified `‖Aᵀ λ‖ < 1e−6` | 96 / 100 | 99 / 100 |
| verified `< 1e−9` | 65 / 100 | 57 / 100 |
| verified `< 1e−5` | 98 / 100 | 100 / 100 |
| worst verified residual | `3.9e−5` | `1.8e−6` |

Per family at `X_ini`, verified `< 1e−9`: Voronoi 34/34, quad-random 31/33, Delaunay 0/33
(29/33 at `1e−6`; `dim_flex` up to 156).

LP margin `> 0` on the K1a population with the corrected solver: 0/100 at `X_ini`, 0/100 at
`X0` — identical to the original, so the 0/100 headline stands. The controls improve, the
second signature of the old primal's under-convergence: `hexagons_auto` `5.40e−1 → 6.21e−1`,
`truncated_square_488` `7.35e−1 → 9.25e−1` (this family's magnitude depends on the
degenerate-row normalisation; the recheck pass printed `8.52e−1`), `snub_square_33434`
`3.53e−2 → 6.58e−2`, `tiling_3_4_3_12` unchanged at `3.12e−1`. The main table above is the
canonical `k8a.csv` with the corrected solver at 100 000 dual steps. A verified residual `r`
with `margin = 0` pins the chart's LP value to `[0, r]` — a numerical bound of those charts,
never of `P(X)`.

Artifacts: `k8a/{k8a.csv, summary.txt, summary_*.txt, recheck.md, recheck/, recheck100k/}`.
Drivers `exp_k8a_expansive_cone.jl`, `exp_k8a_expansive_cone_recheck.jl`, summary
`scripts/summarise_k8a.jl`; method `src/method/expansive_cone.jl`.

---

## B3 — `tr K` as the expansion budget

Driver `exp_b3_expansion_budget.jl`, library `src/method/budget.jl`; outputs
`b3/{b3_b1_global.csv, b3_retro.csv, b3_constrained.csv, summary.txt}`.

**Question.** Is the total first-order hole-opening rate a functional of the border alone
(B1), does `tr K` bound it on a torus such that `tr K ≥ τ*` is a linear constraint on the
achievable set (B3), and does imposing that constraint raise K7's certified count?

**B1, global half.** On 10 fixed-boundary patches the total void area (the shoelace of the
deployed border half-edges minus the θ-independent face-area sum, so holes and notches with
no walk tracing) equals the border functional's first harmonic to `3.41e-15` at 6 angles,
and the edge-wise budget satisfies `W + R = 2 B(X)` to `2.92e-16`. Bar `1e-10`: pass.
Notches are not a correction: `Σ_C a_C` is 28–46 % below the hole-only sum on every pattern.

**B3's identity needs a repair.** Both copies of an edge contribute
`⟨x_a − x_b, u_{f1} − u_{f0}⟩`. For a split edge `σ_{f0} = σ_{f1}` so the term is
lattice-translation invariant, but for a hinge edge `σ_{f0} = −σ_{f1}` and
`u_{f1+t} − u_{f0+t} = u_{f1} − u_{f0} + σ_{f1} t`, so the term depends on which lift the
preimage walk uses. Summing one arbitrary representative per quotient edge misses a defect
of `0.38–1.88 × det P₀` on the 33 K7 patterns. Summing one representative preimage per
translation class at its own lift (`method::periodic_cell_edges`) gives
`W + R = det P₀ · tr K` to `3.9e-14` on 46 of 50 C3 rows — every row of 23 of the 25
patterns with `dim 𝒦 ≥ 1` (exceptions: `squares_3x2`, whose 3×3 super patch carries no full
representative set, and `snub_square_3x3`, whose `K` is gauge-dependent). Consequence:
`tr K − τ* = R / det P₀` exactly, so `τ*` is a function of the design rather than of the
pattern, and `tr K ≥ τ*` is not a linear constraint on the achievable set.

**Retrodiction.** The baseline reproduces `k7_c3_all.csv` bit for bit (0 mismatches over 50
rows; 13 certified, 21 with `min q > 0`).

| quantity | value | bar |
|---|---|---|
| AUC(`tr K − τ*`) for `min q > 0` | 0.989 | 0.80 |
| AUC(`tr K − τ*`) for certified | 0.906 | 0.80 |
| rows with `min q > 0` at `tr K < τ*` | 0 | must be 0 |
| same, restricted to the 46 identity-exact rows | 0.987 / 0.895 | — |
| AUC(`R/(det P₀ n_split)`) for `min q > 0` / certified | 0.982 / 0.963 | — |
| spread(W) / spread(det P₀ tr K) over 𝒦, 20 points × 25 patterns | median 0.776, max 1.401 | hard fail |

The last row is the clause that decides B3: `W` moves over the achievable set by ~0.78 of
what `det P₀ tr K` moves, so `τ*` is a penalty, not a threshold. (Random-target rows
inherit K7 C3's basis-dependent draw; the ratio and the identity checks do not.)

**Constrained re-solve.** Minimising `‖K − K*‖` over 𝒦 subject to `tr K ≥ τ* + δ|τ*|`; at
`dim 𝒦 = 4` the minimiser is `K* + ((τ*(1+δ) − tr K*)/2) I`, an isotropic expansion added to
the target; `τ*` refreshed by a 4-step fixed point.

| δ | certified | flips | `min q > 0` | max ‖K − K*‖ |
|---|---|---|---|---|
| 0.0 | 13 → 15 | +2, −0 | 21 → 37 | 10.4 |
| 0.1 | 13 → 14 | +2, −1 | 21 → 37 | 19.8 |
| 0.5 | 13 → 2 | +1, −12 | 21 → 25 | 588.9 |

`+2` and `+1` of 50 are inside the binomial error bar (±4.7 at `p = 0.26`). The constraint
nearly doubles the designs whose split cuts all open at `0⁺` (21 → 37), but that does not
convert into certified range — the K5/K6 obstruction again — and the target is abandoned by
10–589 in Frobenius norm.

**Verdict:** theorem pass after the repair to its derivation; algorithm fail; B3's hard-fail
clause triggered.

---

## T-1 — the vertex balance defect

**Question.** At an interior vertex `v` touched by no split edge (a *pure* vertex) the
barycentric row of Eq. (2) reads `δ_v := Σ_{u ∈ In(v)} (x_u − x_v) = 0`. Read
kinematically, the `k_v` copies of `v` separate at `0⁺` with velocities proportional to
`J(x_u − x_v)`, which sum to zero and for `k_v ≥ 3` positively span the plane. Does that
keep every copy of `v` out of every neighbouring face's material, so that the `0⁺` collapse
is localised on split-touched vertices? The competing reading: the `0⁺` corner predicate is
`mu = max(−g1, −g2)` at a convex corner, a disjunction, so positively spanning velocities do
not by themselves keep every copy out of every corner cone. Test: on the K6 population at
`X0`, count the pure interior vertices and those with `mu_v ≤ 0`. Pass iff `N_bad = 0`.

**Measured** (`t1/t1_summary.txt`, 200 graphs × 2 σ = 400 designs at K6's cached `X0`,
wall 1.5 s):

| quantity | value |
|---|---|
| interior vertices with ≥ 2 copies | 165 788 |
| pure (no incident split edge) | 14 886, on 250 / 400 designs, median pure fraction 0.0121 |
| touched (≥ 1 incident split edge) | 150 902 |
| `N_bad` (pure, `mu_v ≤ 0`) | 1 970, on 142 / 400 designs |
| touched with `mu_v ≤ 0` | 66 521 |
| restricted to `k_v ≥ 3` | pure 4 720, bad 1 014 |
| restricted to `k_v = 2` | pure 10 166, bad 956 |
| binding incidence at a bad pure vertex | convex corner 1, reflex corner 1 969 |
| pure vertices with every incidence at a convex corner | 12 916, bad 0 |
| faces non-convex at `X0` / at the input embedding | 22 107 / 156 220 ; 0 |
| collapse rate, pure vs touched | 0.1323 vs 0.4408 |
| pure vertices with `\|δ_v\| > tol` | 0 (max `3.425e−14`, the F11 self-check) |
| corr(`\|δ_v\|`, `mu_v`) on touched vertices | Spearman 0.0617, Pearson −0.1016 |
| culprit split edge incident to the argmax-`\|δ\|` vertex / to a top-10 % vertex | 61 / 400 ; 254 / 400 |

**Fail as stated, and the failure is the finding.** Pure vertices do collapse (1 970 of
14 886), so the balance row does not protect them — but 1 969 of the 1 970 bad pure
vertices sit at a *reflex* corner, and among the 12 916 pure vertices whose every incidence
is at a convex corner there are none. The input embeddings have zero non-convex faces;
`X0` has 22 107 non-convex corners out of 156 220. The reflex corners are manufactured by
the Eq. (6) projection, and they are where the corner disjunction fails. The balance defect
itself is a weak predictor (correlation ≈ 0). This is the observation K9 acts on.

Artifacts: `t1/{t1.csv, t1_graph.csv, t1_summary.txt}`. Driver `exp_t1_balance_defect.jl`
(reads K6's shape-space cache and K5's saved `σ_def`; nothing recomputed).

---

## K9 — convexity-constrained embedding

**Question.** T-1 measured that 1 969 of 1 970 balanced pure vertices with a non-positive
`0⁺` margin sit at a reflex corner, and that the input embeddings have zero non-convex
faces while `X0` has thousands of non-convex corners. At a convex corner the margin is
`μ = max(−g₁, −g₂)`, a disjunction; at a reflex corner it is `min(−g₁, −g₂)`, a conjunction,
and that is the one that fails. The reflex corners are created by the unconstrained Eq. (6)
projection, not by the graph. Does putting them back inside the shape space buy range that
K6's `q`-only repair could not? Bar: refereed exact `Θ_max > 0` on ≥ 10 % of 400 designs.

**The problem** (`src/method/convex_embed.jl`). Inside `X(t) = X0 + Φ t`,

```
minimise ‖X(t) − X_ini‖²_F   subject to   cross_i(t) ≥ δ   for every corner i of every face
                                          [ and q_e(t) ≥ δ′ for every split cut ]
```

with `cross = det(X[v] − X[v_prev], X[v_next] − X[v])`. Every corner strictly convex implies
every face convex and positively oriented, so this subsumes K6's face-area constraint. Both
families are quadratic in `t` and the feasible set is not convex; phase A minimises a
softplus penalty under a continuation over the target (`0.01, 0.1, 1` × the solve target),
phase B refines the proximity behind a log barrier keeping only strictly feasible iterates.
FEASIBLE always means the exact constraint values at the returned point.

**Setup.** K6's 200 graphs under both orientation rules, `σ_def` read from `k5/sigma/`.
`δ = δ′ = 1e−3·med²`, 3 Gaussian starts of scale `0.1·med` plus `t = 0`, L-BFGS 600
iterations per continuation stage, 6 barrier stages, certificate `ε = 0.3`. Three arms:
`X0` (control), (a) convexity only, (b) convexity + split-inward, (a)'s minimiser warm-starting
(b). 12 shards, 1.60 h of solver time.

**Measured** (`k9/{k9.csv, summary.txt}`, 400 rows):

| quantity | `X0` (control) | (a) convexity | (b) convexity + split |
|---|--:|--:|--:|
| convex-feasible designs | 79 / 400 (19.8 %) | 118 / 400 (29.5 %) | 30 / 400 (7.5 %) |
| exact `Θ_max > 0` | 0 | 0 | 36 / 400 (9.0 %) |
| refereed (bisection, fixed predicate) | 0 | 0 | 36 |
| refereed (pre-F34 predicate at shrink `1e−12`, earlier pass) | 0 | 0 | 32 |
| certified (`ε_max > 0`) | 0 | 0 | 36 |
| all `0⁺` corner margins `μ > 0` at the returned `X` | 0 / 79 | 0 / 118 | 23 / 30 |
| median `min μ / med²` over feasible | −8.79 | −10.06 | +4.96e−2 |
| `‖X − X_ini‖` per vertex, median / q90 (median edges) | 0 / 0 | 0.275 / 0.399 | 0.376 / 0.567 |
| expansive-cone LP feasible (corrected solver) | — | 0 / 118 | 29 / 30 |
| binding contact at `Θ_max = 0` | split-inward 400 | split-inward 400 | split-inward 259, inverted 73, vertex-edge 32 |

(The (a) arm's convex-feasible count and the (b) arm's endpoints are optimiser-path
dependent — an earlier pass gave 115 under (a) — while the 36 positives, their `Θ_max` and
the `t = 0` columns reproduce.)

Per orientation, arm (b):

| quantity | `σ_mc` | `σ_def` |
|---|--:|--:|
| convex-feasible | 13 / 200 | 17 / 200 |
| exact `Θ_max > 0` | 24 / 200 (12.0 %) | 12 / 200 (6.0 %) |
| refereed (fixed predicate; pre-F34 in brackets) | 24 (20) | 12 (12) |
| `ε_max` median / max | 0.278 / 1.997 | 0.152 / 0.381 |
| `‖X − X_ini‖` per vertex, median | 0.523 | 0.297 |
| expansive-cone LP feasible | 13 / 13 | 16 / 17 |

Over the 36 positives: `Θ_max` median 0.249, min 0.027, max 1.997 rad; 33 of 36 have
`ε_max ≥ 0.1` rad; by family delaunay 25, voronoi 10, quad_random 1. Taking the better `σ`
per graph, 33 of 200 graphs carry a deployable design. Only 23 of the 36 positives are
`b_feasible` at the full margin `δ′`: the margin is a sufficient handle on the 0⁺ cone, not a
necessary one.

Convexity alone is necessary, not sufficient: arm (a) raises convex feasibility from 79 to
118 designs and buys zero range — the binding contact stays `split-inward` on all 400. Both
constraints together move the count off zero, the lesson K6 reached from the other side.
Non-convex corners at `X0`: 35 336 / 656 736 = 5.38 %.

**The referee audit (F34).** Four of the 36 exact positives (delaunay 13, 40, 103, 190, all
`σ_mc`) refereed as `Θ_bisect = 0` at the historical shrink `1e−12`. `apps/dbg_k9.jl` traces
each to the bisection's early-exit `has_collision(1e−7)` firing on a pair of hinge-adjacent
faces, which touch at the pin by construction, non-monotonically in `θ` (true at
`1e−8, 1e−7, 1e−5`, false at `1e−9, 1e−6, 1e−4, …, 0.1`): the shrink displaced the two copies
of the shared pin by different amounts. `polygons_overlap` was rewritten
(`results/core_validation/referee_fix.md`) and K9 re-run in full: the referee agrees with
the exact scan on 400/400 with worst gap 0, the refereed count is 36/36, and the only four
cells that changed are these four (0.219029, 0.177931, 0.301957, 0.103354, matching the
`1e−9` re-referee). Soundness (`ε_max ≤ bisection Θ_max`) held on 396/400 before the fix and
400/400 after.

**The expansive-cone LP control.** As first run, on the 30 (b)-feasible designs the LP found
no witness (0/30) although on 23 of them `σ` is in `P(X)` by direct measurement, so a
feasible point provably exists. That was `cone_lp` non-convergence (§K8a recheck); with the
corrected solver the same 30 embeddings give 29 feasible, all 23 with `σ ∈ P(X)` among them,
at a median 8.7× the uniform ray's margin, plus 6 of the 7 where `σ` is outside `P(X)`. The
0/30 figure is not a result.

**Fail: 9.0 % (36 of 400) against 10 %**, missed by 4 designs while every comparable baseline
near the projection reports zero or nearly so: Eq. (6) alone 0/200 (K1a), `σ_def` and `σ_mc`
0/200 (K5), K6's primary repair 2/400 (split-only 15/400, by leaving the neighbourhood), B4
2/400, K8a 0/100, native `prevent` 0/80. Against that, 36 designs with a certified `ε_max` up
to 2.0 rad is the first constructive result on random graphs that stays near the
projection. The random-graph design space is empty under unconstrained Eq. (6), and sparsely
non-empty when convexity and the split signs are imposed together inside the null space.

Artifacts: `k9/{k9.csv, summary.txt, shard_*.csv, shard_*.txt}`. Driver
`exp_k9_convex_embedding.jl`, method `src/method/convex_embed.jl`, diagnostic
`apps/dbg_k9.jl`; unit tests in `test_method.jl` (`corner_crosses`, a hand-known constrained
minimiser on the 2×2 grid, the phase-A gradient against a central difference, the
FEASIBLE-is-exact invariant, barrier-stage monotonicity).

---

## K9b — pushing the constrained embedding

**Question.** K9 missed its bar by 4 designs. Holding the population, the shape space and
the certificate fixed, is 9 % a property of the design space or of the search?

**Four levers** (`exp_k9b_convex_embedding_search.jl`):
* (i) per-graph best of `σ` — a reporting convention, given as "best-of-2 over graphs"
  beside the per-design count, never mixed into it.
* (ii) 8 Gaussian starts + `t = 0` (K9: 3 + 0), plus the two warm starts K9 used; 10
  barrier stages (K9: 6).
* (iii) a sweep of `δ = δ′` over `{1e−4, 1e−3, 3e−3, 1e−2}·med²`, best-first, cut short
  once a design is certified past `ε = 0.1` or a 200 s per-design budget is spent.
* (iv) stage 2: `range_opt.jl`'s softmin-of-first-contact objective re-centred at the
  feasible point, accepted only if the exact constraints still hold and the exact `Θ_max`
  improved. The barriers are enforced by exact rejection, not inside the stage-2 objective.

Every design refereed at both shrinks, `1e−12` and `1e−9`. 12 shards, 7.69 h of solver time.

**Measured** (`k9b/{k9b.csv, summary.txt}`, 400 rows):

| quantity | K9 (b) | K9b sweep | best over both |
|---|--:|--:|--:|
| exact `Θ_max > 0` | 36 / 400 (9.0 %) | 28 / 400 (7.0 %) | 37 / 400 (9.25 %) |
| refereed at shrink `1e−12` | 36 (32 pre-F34) | 28 | 37 |
| refereed at shrink `1e−9` | 36 | 28 | 37 |
| certified (`ε_max > 0`) | 36 | 28 | 37 |
| of which `ε_max ≥ 0.1` rad | 33 | 16 | 34 |
| `ε_max` median / max (rad) | 0.253 / 1.997 | 0.147 / 1.828 | 0.253 / 1.997 |
| convexity + split feasible | 30 | 28 | — |
| graphs with ≥ 1 deployable `σ` (of 200) | 33 | 26 | 34 |

Per orientation: `σ_mc` 22/200 (11.0 %), `σ_def` 6/200 (3.0 %). Worst
`|exact − bisection(1e−9)|` `4.1e−06` over all 400; soundness 400/400. Sixteen rows have a
bisection value below `5e−6` with the exact scan at 0; that is below the bisection's grid
resolution and is counted as zero.

**Which lever mattered.** Only the `δ` sweep: among the 28 positives the winning margin is
`δ = 1e−2·med²` on 15, K9's `1e−3` on 12, `1e−4` on 1 — the constraint should be posed with
a generous margin. Restarts and barrier stages did not pay; stage 2 supplied the accepted
point on 2 of 28 and, being gated on `Θ_max > 0`, never ran on the feasible-but-jammed
designs where it would have mattered (a design error in this driver).

**The 9 lost designs.** K9b keeps 27 of K9's 36, adds 1 (`delaunay 22 σ_def`,
`Θ_max = 0.151`), loses 9 — all `delaunay`, and on all nine K9b's returned point is exactly
feasible at the full margin and still has `Θ_max = 0`, binding on `vertex-edge`. Two solvers
found different points inside the same feasible set, one deployable and one not. Margin
feasibility does not determine deployability, and the proximity objective that selects the
point inside the feasible set is not aligned with range. That, not the search budget, caps
the count near 9 %.

**Fail:** 28 (K9b alone), 37 (best over both) against ≥ 40. The constructive claim survives
unchanged in kind: 37 random graphs deploy from a constrained embedding near the projection,
against 2 without it.

Figures: `k9b/k9b_gallery.png` (24 random graphs closed and at `θ = Θ_max/2`, `F` from 101
to 442), `k9b/k9b_hist.png`. Artifacts `k9b/{k9b.csv, summary.txt, shard_*, gallery/}`;
plots `scripts/plot_k9b.jl`.

---

## K9c — range-maximising embedding

**Premise.** K9b's nine lost designs show that the feasible set is not what caps the yield;
the objective that selects the point inside it — proximity, `‖X − X_ini‖` — is. K9c keeps
the population (400 designs), the shape space, the barriers (`cross_i ≥ δ`, `q_e ≥ δ′`) and
the certificate fixed, and replaces only the objective.

**What is maximised** (`src/method/range_embed.jl`). The 0⁺ deployment margin, in units of
`med²`, `m(X) = min( min_e q_e(X), min_j mu_j(X) ) / med²` — the joint minimum of the split
sign `q_e` and the corner margin `mu_j`. Both are quadratic in `t`; `mu` carries the max/min
disjunction whose branch is the sign of the same `cross` that `convex_embed.jl` constrains,
so inside the convexity barrier the reflex branch cannot occur. The margin list is flattened
into smooth entries and the objective is the log-sum-exp softmin
`M_kappa(t) = −kappa · log Σ_i exp(−c_i(t)/kappa) ≤ min_i c_i(t)`, a lower bound on `m` for
every `t` and tight where the branch modes were read (both claims unit-tested in
`test_range_embed.jl`, plus finite-difference checks of the analytic gradient).

Stage A minimises `−M_kappa + w·(convexity + split log barriers)` with `w` decreased
geometrically over 6 continuation stages, branch modes refreshed between them; the point kept
is the one with the largest exact `m`. Stage B, from a point with `m > 0`, pushes the exact
`Θ_max` with `range_opt.jl`'s objective under a shrinking trust region (caps `0.25, 0.10,
0.04` median edges), accepting a step only if the exact convexity, the exact split signs, a
margin floor (`m ≥ ½ m₀`) and an actual `Θ_max` improvement all hold.

**Setup.** Three arms on every design, the best recorded with its provenance:
* arm `k9` — `convex_embed` with K9's settings, proximity objective.
* arm `k9b` — `convex_embed` with `δ = δ′ = 1e−2·med²`, 8 starts, 10 barrier stages, warm
  started from arm `k9` (K9b-like, not K9b: 27 of K9b's 28 positives won at `1e−3` or
  `1e−2`).
* arm `k9c` — stage A from arm `k9`'s point, arm `k9b`'s point and `t = 0` (6 stages × 120
  L-BFGS iterations), then stage B.

Only the winning point of each design is refereed, at both shrinks. 12 shards over the graph
index; all 400 designs completed (the run was first reported at an interim 182 under heavy
machine contention and later finished and re-merged with `--aggregate --nshards 12`).

**Measured** (`k9c/{k9c.csv, summary.txt}`, 400 of 400):

| quantity | arm `k9` (proximity) | arm `k9b` (proximity, `δ=1e−2`) | arm `k9c` (range-max) | best of three |
|---|--:|--:|--:|--:|
| convexity + split feasible | 30 (7.5 %) | 24 (6.0 %) | 258 (64.5 %) | 212 |
| of which 0⁺ margin `m > 0` | 23 | 12 | 239 | 210 |
| exact `Θ_max > 0` | 36 (9.0 %) | 27 (6.8 %) | 312 (78.0 %) | 312 |
| refereed at shrink `1e−12` | — | — | — | 311 |
| refereed at shrink `1e−9` | — | — | — | 312 |
| certified (`ε_max > 0`) | 36 | 27 | 310 | 310 |
| of which `ε_max ≥ 0.1` rad | 33 | 14 | 192 | 201 |
| `ε_max` median / q90 / max (rad) | 0.253 / 0.540 / 2.00 | 0.115 / 0.470 / 1.83 | 0.162 / 0.958 / 3.14 | 0.172 / 0.958 / π |
| 0⁺ margin `m` over feasible, median (med²) | 0.050 | 0.004 | 0.529 | 0.503 |

The K9c arm's endpoints are optimiser-path dependent: which designs cross `Θ_max = 0`, how
far each opens and where the winning point comes from move between floating-point
environments (an earlier pass gave 307 positives, 306 certified, 214 at `ε_max ≥ 0.1` rad
and a median `ε_max` of 0.278 rad); the proximity arms, the `t = 0` columns and the referee
do not. The rate is stable to about ±2 % of the population.

Per orientation, best-of-three: `σ_mc` 156/200 (78.0 %) with 80 at `ε_max ≥ 0.1` rad,
`σ_def` 156/200 with 121. Per family, Voronoi 126/134, Delaunay 109/134, quad-random 77/132.
Taking the better `σ` per graph, 187 of 200 graphs carry a deployable design and 147 one with
`ε_max ≥ 0.1` rad. `‖X − X_ini‖` per vertex is median 0.927 median edges against K9's 0.38
and K9b's 0.21: the range-maximising point is farther from the input embedding, which is the
price of the range.

Yield is not scale-free: `results/final/figures/fig_yield.png` bins by `|F|` with Wilson
95 % intervals; Delaunay falls from about 0.94 to 0.65 across the bins and quad-random from
0.88 to 0.44.

Provenance of the winning point: arm `k9` on 102 (88 of them designs where all arms return
`Θ_max = 0` and the tie-break keeps the first), arm `k9b` on 2, a K9c point on 296
(`k9c/x0` 154, `k9c/x0+B` 127, `k9c/k9+B` 9, `k9c/k9b+B` 4, `k9c/k9` 1, `k9c/k9b` 1). Within
K9c the `t = 0` start supplies the winner on 281 and the warm starts on 15: the range
objective does not need K9's feasible point to start from. Stage B supplied the winning
point on 140 of the 312 positives; stage A alone gets 239 of 400 designs into `m > 0`
(against 23 for the proximity arm) and stage B converts that margin into range.

**Referee.** Soundness (`ε_max ≤ bisection at 1e−9`) 400/400; the two shrinks agree on
every positive but `delaunay 136 σ_def` (`Θ_max = 1.9e−7` rad, below the grid resolution;
not counted). On `delaunay 103 σ_def` the exact scan reports `Θ_max = 0.2986` and the
bisection `0.3294` — the bisection's grid of `π/4000 ≈ 7.9e−4` rad can step over a
collision interval narrower than that; the exact scan tests the overlap status at the
closed-form contact angles themselves and is the smaller, conservative number.

Of the 88 designs with `Θ_max = 0`: `split-inward` 58, `inverted` 24, `vertex-edge` 6. The
`inverted` count is new: when no start reaches feasibility the driver records the
least-infeasible point (`best_feas = 0`) rather than discarding it.

K9b's nine lost designs (`delaunay 10, 37, 43, 64, 112, 163, 190 (σ_def)`; `delaunay 40,
103 (σ_mc)`) are recovered by the K9c arm alone, 9/9 — the direct test of K9b's diagnosis.

**Pass.**

| target | achieved |
|---|---|
| ≥ 40 / 400 refereed exact `Θ_max > 0` | 312 / 400 |
| as many as possible with `ε_max ≥ 0.1` rad | 192 (201 best-of-3) |
| recover K9b's 9 lost designs | 9 / 9 |

Most random graphs deploy from a constrained embedding once the embedding is chosen to
maximise the 0⁺ margin rather than to stay near the input: 78.0 % of 400 designs, against 0
for Eq. (6) alone, `σ_mc`, `σ_def`, convexity alone, K8a and the native `prevent`, and 2/400
for K6's primary repair and B4 (17/400 counting K6's split-only relaxation), all on this
population. The two proximity arms measured here reproduce K9's 9.0 % and K9b's 7.0 %, which
is the control that says the jump is the objective and not a change of population or referee.

Figures: `results/final/figures/` (`fig_gallery_full.png`, `fig_yield.png`,
`fig_feasible_vs_deployable.png`, `summary_table.md`); `k9c/k9c_gallery.png`,
`k9c/k9c_hist.png`. Artifacts `k9c/{k9c.csv, summary.txt, shard_*, gallery/}`. Driver
`exp_k9c_range_embedding.jl`, method `src/method/range_embed.jl`, tests
`test_range_embed.jl`, plots `scripts/plot_k9c.jl`.

---

## Native200 — the authors' full pipeline on the K9 population

**Question.** Every prior comparison against the authors' code was capped (K2b: 30 live
designs; K6: 8 of 400). Run their full pipeline — their `initialized_two_face_coloring`,
their Eq. (6), their Eq. (9), their FK collision test — unconditionally on the same 200-graph
population K5/K6/K9/K9b/K9c use. All native numbers come from their unchanged binary
(`baseline/native/`); only our side of each comparison is recomputed.

**Setup.** `exp_native200.jl`, three variants per graph: *native* (their coloring via the
CLI's `color` subcommand, fed to `prevent`), *sigma_mc*, *sigma_def*. Each (graph, variant)
cell is `tuttekiri_cli prevent … --collisions --dump …` under a 600 s alarm, so a cell is
`completed`, `timed_out` or `crashed` (`Segmentation fault: 11` in every case observed). The
returned embedding is re-evaluated with the exact T4.2″ scan, the bisection referee and the
`ε = 0.3` certificate; the authors' `native_theta_collisions` is recorded as a diagnostic
column only (F24). 12 shards.

**Coverage.** The 600 s pass was stopped at 573/600 cells (172 completed, 362 timed out, 39
crashed). The 39 crashed cells were rerun against the F42-patched binary at 3 600 s
(`crashfix3600/`, 26 rows; F42 traced the crash to `merge_close_verts` being called on a
deployed configuration that had collapsed to zero faces — `native200/crashfix.patch`; on 3
completed cells the patched and original binaries agree exactly), and 4 timed-out cells at
3 600 s (`rerun3600/`). `native200_merge.jl` merges the passes into `native200_final.csv`.
On the 600-cell grid: 574 dispatched (26 never dispatched), 195 completed, 364 timed out, 15
crashed. Timeouts are concentrated where `|F|` is largest — median `|F|` 202 for completed
cells, 482 timed out, 431 crashed, against 332.5 for the population and ≤ 97 for the 2026
paper's own worked examples. `quad_random` is worst hit: 150/198 timed out, 29/198 completed.

| | `native` | `sigma_mc` | `sigma_def` | total |
|---|--:|--:|--:|--:|
| completed | 58 | 64 | 73 | 195 |
| timed_out | 130 | 125 | 109 | 364 |
| crashed | 1 | 4 | 10 | 15 |
| never dispatched | 11 | 7 | 8 | 26 |
| dispatched | 189 | 193 | 192 | 574 |

| family | completed | timed_out | crashed | never dispatched |
|---|--:|--:|--:|--:|
| delaunay | 120 | 67 | 5 | 9 |
| voronoi | 46 | 147 | 1 | 7 |
| quad_random | 29 | 150 | 9 | 10 |

**Measured, on the 195 completed cells.** Their own FK collision test reports `Θ > 0` on
27/195; our exact scan on the same dumped embedding gives `Θ_max > 0` on 0/195, the
bisection 0/195, the certificate 0/195. The 27 are 25 `sigma_def` and 2 `native`; by family
`quad_random` 16, `voronoi` 8, `delaunay` 3. All 27 fail `POS` at `ε = 0.3`: the exact scan
finds an immediate `0⁺` collision that `merge_close_verts` fuses away before their FK test
sees it — the F24/K1c mechanism, now measured unconditionally. The full 27-row table is in
`NATIVE200_FINAL.md` §3; the five from the 600 s pass:

| id | kind | variant | `F` | `dim_null` | native `Θ` | our exact | our bisect |
|---|---|---|--:|--:|--:|--:|--:|
| 160 | delaunay | native | 292 | 59 | 0.0880 | 0.0000 | 0.0000 |
| 36 | voronoi | sigma_def | 131 | 235 | 0.6726 | 0.0000 | 0.0000 |
| 102 | voronoi | sigma_def | 336 | 600 | 0.3508 | 0.0000 | 0.0000 |
| 168 | voronoi | sigma_def | 138 | 238 | 1.0248 | 0.0000 | 0.0000 |
| 156 | voronoi | sigma_def | 308 | 566 | 0.5704 | 0.0000 | 0.0000 |

Timing: completed cells median 13.5 s, max 570.7 s; timed-out cells at the alarm (median
600.2 s; the tail to 3 711 s is the four 3 600 s reruns); crashed cells median 264 s.

**Reading.** The authors' full pipeline, run unconditionally on the K9 population, does not
escape the `0⁺` jam: 0/195 by exact scan, 27/195 by their own predicate, all 27 false
negatives. Caveat: this is 0 % of what finished, on a population whose median `|F|` is 3–4×
the paper's worked examples; the pipeline's practical operating range on this hardware is
smaller graphs than this population supplies, which is itself a finding about scaling,
independent of the `0⁺` obstruction.

Artifacts: `native200/{native200_final.csv, NATIVE200_FINAL.md, native200.csv, summary.txt,
shards/, crashfix3600/, rerun3600/, logs/, run_all.log}`. Driver `exp_native200.jl`, merge
`native200_merge.jl`.

---

## Reproducing

`Kirigami/scripts/rerun_all.sh` runs every driver in dependency order with the committed
flags and wall-time estimates; individual commands:

```
julia --project=Kirigami -e 'using Pkg; Pkg.test()'   # 187 test sets, 176,692 assertions, 23 broken, 0 failures
A=Kirigami/apps; R=results/experiments
julia --project=Kirigami $A/exp_k3a_core_mobility.jl          --out $R/k3a                  # ~25 min
julia --project=Kirigami $A/exp_k3a_core_mobility_recheck.jl  --out $R/k3a                  # ~4 min
julia --project=Kirigami $A/exp_k1b_harmonic_identity.jl      --out $R/k1b                  # ~34 s
julia --project=Kirigami $A/exp_k1a_projection_validity.jl    --out $R/k1a                  # ~40 s
julia --project=Kirigami $A/exp_k5_orientation.jl             --out $R/k5 --n 200 --maxf 800   # ~5 min
julia --project=Kirigami $A/exp_k2c_active_set_locality.jl    --out $R/k2c                  # ~3 min
julia --project=Kirigami $A/exp_k2a_thetamax_bisection.jl     --out $R/k2a                  # ~7 s
julia --project=Kirigami $A/exp_k1c_eq9_false_negatives.jl    --out $R/k1c                  # ~31 s   (needs baseline/native)
julia --project=Kirigami $A/exp_k2b_range_margin.jl           --out $R/k2b --maxf 160       # ~26 min (needs baseline/native)
julia --project=Kirigami $A/exp_f23_boundary_row.jl           --out $R/f23                  # ~1 s
julia --project=Kirigami/scripts Kirigami/scripts/plot_experiments.jl
```

Sharded experiments (A3, K6, K7, B4, K8a, K9, K9b, K9c, Native200) take `--shard I --nshards N`
(`--nshard` for K7/K8a) and are merged by concatenating the shard CSVs; every row is a
deterministic function of its (graph, σ, variant) so shards may be run in any order. The
exact invocations, including the K7 stage sequence and the Native200 alarm wrapper, are in
`rerun_all.sh`. Figures: `k2a/k2a_agreement.png`, `k1c/k1c_rates.png`,
`k3a/k3a_mobility.png`, `k2c/k2c_locality.png`, `k2c/k2c_hloc.png`, `k2b/k2b_margin.png`,
`k5/k5_orientation.png`, `jitter/transition.png`, `k7/k7_*.png`, `k9b/k9b_*.png`,
`k9c/k9c_*.png`, and `results/final/figures/`.
