# WP2 — why the K9c yield is not scale-free

**Verdict: FAIL against the pre-registered bar.** Best held-out AUC on the fresh graphs is
**0.751** (two-feature model `n_incid` + `n_constraints_0plus`), against a bar of 0.85. The
yield half of the bar is met: that model predicts a fresh yield of 0.785 against an observed
157/200 = 0.785, Wilson 95 % [0.723, 0.836]. The AUC half is not.

The outcome columns of `features.csv` are K9c's, and K9c's endpoints are optimiser-path
dependent (`results/kill/KILL_REPORT.md` §K9c): the labels of about 5 % of the 400 designs
and the exact AUC values and selected feature pairs below move between floating-point
environments (an earlier pass of the same analysis, on an earlier K9c pass with 307
positives, selected `log_F` + `hole_size_mean` at a held-out AUC of 0.739). The verdict, the
size and family trends and the mechanism ranking do not.

Every number below is reproduced by
`arch -arm64 /usr/local/bin/python3 results/yield/analyse.py`, which writes
`results/yield/analysis_out.txt`; section references are to that file.

## What was measured

`Kirigami/apps/kill_yield.jl` (added to `Kirigami/Project.toml`) computes, per design, 45 raw
features (52 candidate columns after `analyse.py` adds the log transforms) that depend only on the input — the graph, the orientation σ, the input embedding
`X_ini`, and the cut structure derived from them — plus a separately labelled group of
**projection-side** quantities produced by the Eq. (6) solve (`proj_*`). No feature is a
K9c solver output; the solver outputs sit in the outcome columns and are excluded from the
bar.

- `results/yield/features.csv` — **400 rows**, the K9c population (ids 0..199, both σ),
  same graphs and same σ as `Kirigami/apps/kill_k9c.jl`, the Eq. (6) projection read from the
  `results/kill/k6/cache` shape cache, outcomes joined from `results/kill/k9c/k9c.csv` on
  (id, kind, sigma). `N`, `F`, `n_split`, `dim_null`, `n_corner`, `med_edge` and
  `nonconvex0` were checked to agree with `k9c.csv` on all 400 rows. K9c was not rerun.
- `results/yield/fresh.csv` — **200 rows**, held out: ids 1000..1099, σ_mc from
  `make_graph` and σ_def from `orientation_defect` (cap 20|F|, seed 7000+id, K5's
  rule), solved by `design_range_max` (`Kirigami/src/method/design.jl`) at the run defaults with seed
  9300 + 7·id + which. `kill_yield.jl --mode check` records 0 fingerprint collisions
  with the 200 K9 graphs.
- `results/yield/analyse.py` — AUC, Spearman, Wilson and the logistic fit implemented on
  numpy; scikit-learn is not installed and is not used.

**Fresh-graph outcome:** exact `Θ_max > 0` on **157/200 = 0.785**, Wilson 95 %
[0.723, 0.836]; `eps_max ≥ 0.1` on 107/200. The K9c population's 312/400 = 0.780,
[0.737, 0.818], sits inside that interval, so the held-out population reproduces the
headline yield.

Verified along the way: `frac_convex_ini` is exactly 1.0 on all 400 designs — every face
corner is convex at `X_ini`. That is the K1a claim: the reflex corners are created by the
Eq. (6) projection, not by the input graph. `min_mu_ini` and `min_q_ini` are negative on
all 400, so no input embedding is 0⁺ feasible before the solver runs.

## 1. Univariate table (400 designs)

AUC is for `theta_exact > 0`; ρ is Spearman with `eps_max`. AUC > 0.5 means a larger
feature value makes deployment more likely. Log transforms of the size counts are added by
`analyse.py` and are rank-identical to the raw counts, so they share a row.

| feature | AUC | ρ(eps_max) | group |
|---|---|---|---|
| `min_mu_ini` — worst corner-incidence 0⁺ margin at `X_ini`, in med² | **0.747** | +0.41 | input |
| `ang_min` — smallest corner angle at `X_ini`, degrees | 0.745 | +0.26 | input |
| `F` (= `log_F`) | 0.269 | −0.52 | input |
| `proj_min_q_x0` — worst split margin at `X0` | 0.730 | +0.25 | projection |
| `min_q_ini` — worst split-edge 0⁺ margin at `X_ini` | 0.726 | +0.36 | input |
| `n_incid` — corner incidences | 0.280 | −0.45 | input |
| `n_constraints_0plus` = `n_split` + `n_incid` | 0.297 | −0.43 | input |
| `n_hinge` | 0.304 | −0.56 | input |
| `hole_size_mean` | 0.644 | +0.26 | input |
| `mean_split_comp` | 0.359 | −0.05 | input |
| `E_int`, `n_corner`, `E` | 0.366 / 0.370 / 0.374 | ≈ −0.44 | input |
| `frac_q_pos_ini` | 0.633 | +0.10 | input |
| `proj_frac_q_pos_x0` | 0.630 | +0.17 | projection |
| `asp_min` | 0.628 | +0.12 | input |
| `dof_ratio` = 2·dim_null / n_constraints_0plus | 0.610 | +0.02 | input |
| `asp_med` | 0.595 | −0.09 | input |
| `split_density` = n_split / E_interior | 0.551 | −0.00 | input |
| `asp_max` | 0.558 | −0.11 | input |
| `med_edge` | 0.534 | +0.31 | input |
| `hinge_diam` | 0.534 | −0.21 | input |
| `proj_nonconvex0` | 0.491 | −0.27 | projection |
| `frac_mu_pos_ini` | 0.452 | +0.14 | input |
| `proj_dist_x0` = ‖X0 − X_ini‖ per vertex, in med | 0.451 | −0.17 | projection |
| `frac_convex_ini` | constant 1.0 | — | input |

The complete 52-row table is `analysis_out.txt` §1; `results/yield/fig_yield_auc.png` plots
the top 20.

**Solver-side, reported and excluded from the bar:** `best_margin` AUC 0.999 (the outcome
restated — K9c returns `Θ_max > 0` almost exactly when its point has a positive 0⁺ margin),
`best_feas` 0.833, `secs` 0.502, `dist` 0.148.

## 2. The size / family confound

Yield per family at the `plot_final.jl` bins `[100, 230, 420, 800]`, Wilson 95 %:

| family | F<230 | 230≤F<420 | F≥420 |
|---|---|---|---|
| voronoi | 33/34 = 0.97 [0.85, 0.99] | 44/46 = 0.96 [0.85, 0.99] | 49/54 = 0.91 [0.80, 0.96] |
| delaunay | 60/64 = 0.94 [0.85, 0.98] | 27/36 = 0.75 [0.59, 0.86] | 22/34 = 0.65 [0.48, 0.79] |
| quad_random | 14/16 = 0.88 [0.64, 0.97] | 32/46 = 0.70 [0.55, 0.81] | 31/70 = 0.44 [0.33, 0.56] |
| pooled | 107/114 = 0.939 [0.879, 0.970] | 103/128 = 0.805 [0.728, 0.864] | 102/158 = 0.646 [0.568, 0.716] |

**The family gap survives matching on |F|.** In the largest bin the families are
0.91 / 0.65 / 0.44, and the voronoi and quad-random intervals do not overlap. Size alone
does not explain the family ordering, and family alone does not explain the size trend:
yield falls with |F| inside every family. Both patterns reappear on the fresh graphs
(`analysis_out.txt` §6): 0.90 / 0.85 / 0.66 by bin, and 0.89 / 0.79 / 0.67 by family in the
same order.

Matched on `split_density` terciles (edges 0.075 / 0.301 / 0.334 / 0.631) the families
barely share a bin — voronoi occupies the top two terciles (mean density 0.471), delaunay
the bottom two (0.215), quad-random the outer two (0.307). That matching is therefore
uninformative, and what it does show points the wrong way: the family with the **highest**
split density has the **highest** yield, and `split_density` on its own has AUC 0.551.
Split-cut density is not the cause.

Family medians of the leading features:

| family | `min_mu_ini` | `ang_min` (°) | `asp_med` | `dof_ratio` | `n_incid` |
|---|---|---|---|---|---|
| voronoi | −7.30 | 20.41 | 7.85 | 0.317 | 2426 |
| delaunay | −8.43 | 0.46 | 2.01 | 0.076 | 1817 |
| quad_random | −16.79 | 0.28 | 2.59 | 0.137 | 3644 |

The family ordering of the yield is the family ordering of `min_mu_ini`. It is **not** the
ordering of aspect ratio (voronoi has the worst median aspect ratio and the best yield) and
not the ordering of split density.

Figure: `results/yield/fig_yield_matched.png`.

## 3. Held-out models

Selection uses 5-fold cross-validated AUC **inside the fit set only**; the held-out AUC is
reported and never used to choose. Solver outputs are not candidates.

| fit set | model | fit CV AUC | AUC on the other σ half | AUC on fresh (200) | predicted fresh yield |
|---|---|---|---|---|---|
| σ_mc (200) | single `split_depth` | 0.738 | 0.664 | 0.584 | 0.390 ✗ |
| σ_mc (200) | pair `proj_dist_x0` + `log_n_hinge` | 0.785 | 0.663 | 0.564 | 0.892 ✗ |
| σ_def (200) | single `mean_split_comp` | 0.812 | 0.742 | 0.620 | 0.882 ✗ |
| σ_def (200) | pair `log_n_incid` + `hinge_frac_largest` | 0.838 | 0.717 | 0.739 | 0.785 ✓ |
| all 400 | single `min_mu_ini` | 0.742 | — | 0.633 | 0.789 ✓ |
| all 400 | pair `n_incid` + `n_constraints_0plus` | 0.789 | — | **0.751** | 0.785 ✓ |

✓/✗ is whether the predicted yield lies inside the observed fresh Wilson interval
[0.723, 0.836]. The "all 400" rows are an extra fit the spec does not require, added
because it is the strongest legitimate test available: the fresh graphs are disjoint from
the K9 population, so a model selected by CV on all 400 is still evaluated out of sample.
It did not help.

Two things the numbers say clearly. **Selection transfers badly across σ**: every
two-feature model loses 0.07–0.12 AUC going from its own CV estimate to the other σ half,
and the σ_mc pair collapses on the fresh graphs (0.785 → 0.663 → 0.564). **The margin
features transfer worse than the size features**: `min_mu_ini` drops 0.747 → 0.633 from the
K9c population to the fresh graphs, while `log_F` moves 0.731 → 0.705 and `ang_min` only
0.745 → 0.686 (`analysis_out.txt` §6).

Figure: `results/yield/fig_yield_predictor.png` (top feature vs yield, `theta_exact` and
`eps_max`, fresh graphs overplotted in orange).

## 4. Binding modes of the 88 failures

| binding | n | median `min_mu_ini` | median `split_density` | median \|F\| |
|---|---|---|---|---|
| n/a (deploys) | 312 | −8.25 | 0.326 | 308 |
| split-inward | 58 | −17.50 | **0.175** | 622 |
| inverted | 24 | −16.65 | 0.444 | 574 |
| vertex-edge | 6 | −15.23 | 0.446 | 480 |

**Answer to the question the spec asks: no.** Split-inward binding sits at *low* split
density — median 0.175 against 0.326 for the designs that deploy — and at large |F| (622
against 308). Every failure mode is concentrated at large |F| and at a much deeper worst
0⁺ violation (`min_mu_ini` −15 to −18 against −8.3); what the density selects is *which*
mode binds, in the opposite direction to the naive hypothesis. The fresh graphs give the
same mode mixture (`analysis_out.txt` §4): split-inward 29, inverted 11, vertex-edge 3 of
the 43 failures.

## 5. Mechanism probes

Power-law fits over the 400 designs (`analysis_out.txt` §5):

| quantity | exponent in log \|F\| | R² |
|---|---|---|
| `n_constraints_0plus` | 1.021 | 0.923 |
| `n_incid` | 1.009 | 0.928 |
| `n_corner` | 1.070 | 0.846 |
| `dim_null` | 1.230 | 0.482 |

`dof_ratio` = 2·dim_null / n_constraints_0plus has median 0.139 and Spearman **+0.111**
with |F| — flat, if anything slightly increasing — and AUC 0.610.

What does degrade with size is the depth of the worst violation:

| bin | n | median `min_mu_ini` | median `min_q_ini` | median `n_incid` |
|---|---|---|---|---|
| F<230 | 114 | −5.44 | −10.53 | 1094 |
| 230≤F<420 | 128 | −9.07 | −14.94 | 2197 |
| F≥420 | 158 | −14.06 | −23.18 | 4527 |

Spearman(`min_mu_ini`, log|F|) = −0.638 on the 400 and −0.649 on the fresh 200. Within
matched |F| bins `min_mu_ini` keeps an n-weighted AUC of 0.660 and `ang_min` 0.716, while
`split_density` keeps 0.595 and `dof_ratio` 0.625. Conversely |F| keeps AUC 0.362 (0.638
with the sign flipped) inside terciles of `min_mu_ini`. Size and margin depth both act, and
neither is a proxy for the other.

## Mechanism

**The yield falls with size not because the 0⁺ constraints outgrow the degrees of freedom —
they do not, `n_constraints_0plus ∝ |F|^1.02` against `dim_null ∝ |F|^1.23`, with
`dof_ratio` flat in |F| (ρ = +0.11) and near-useless as a predictor (AUC 0.610) — but
because 0⁺ feasibility is a *minimum* over a constraint set whose size grows linearly with
|F|, so the single deepest violation the solver must repair deepens with the count: the
median worst corner-incidence margin at `X_ini` falls from −5.4 med² to −14.1 med² as the
median incidence count goes 1094 → 4527, and that depth (`min_mu_ini`, AUC 0.747, with the
smallest corner angle `ang_min` at 0.745 right behind it) is the strongest input-side
predictor there is, while split-cut density (0.551) and face aspect ratio (0.558 / 0.595)
have almost none.** The same quantity orders the families: voronoi −7.3,
delaunay −8.4, quad-random −16.8, exactly the yield ordering, and not the ordering of
aspect ratio or density. This is stated as the mechanism the numbers *support*, not as a
sufficient explanation — see the verdict.

## Verdict: FAIL

Pre-registered bar: held-out AUC ≥ 0.85 on the fresh graphs for one feature or a
two-feature model, **and** the model's predicted fresh yield inside the Wilson 95 %
interval of the observed fresh yield.

- Best held-out AUC obtained on the fresh graphs: **0.751**, from the two-feature model
  `n_incid` + `n_constraints_0plus` fitted on all 400 K9c designs. Best single feature on
  fresh among the selected models: `min_mu_ini` at 0.633 (fitted on all 400).
- That model's predicted fresh yield is 0.785, inside [0.723, 0.836]. Three of the six
  models pass the yield half of the bar; none passes the AUC half.
- The in-population ceiling is itself below the bar: the best cross-validated AUC on the
  400 is 0.789 for any pair of these 52 candidate features. No feature set measured here separates
  deployable from non-deployable designs at 0.85.

So the reviewer's question is answered in the negative direction as well as the positive
one: the trend is real, reproduces on 200 held-out designs, and is driven by size and by
the depth of the worst 0⁺ violation rather than by split density, aspect ratio or the
constraint/DOF budget — but the input alone does not determine the outcome at the accuracy
the bar demanded. Roughly a quarter of the ranking information is missing from every
input-side and projection-side quantity measured. The remaining candidate named in
REPORT.md "Next steps" 7 is the solver's basin — whether the same design flips outcome
under a different seed or start set; `results/yield/BASIN.md` (WP2b) tests it and finds
the failure set predominantly geometric.

## Features that did **not** predict (AUC within 0.05 of 0.5, on the 400)

`hinge_diam` 0.534, `med_edge` 0.534, `n_split_comp` 0.537, `n_border` 0.549,
`frac_mu_pos_ini` 0.452, `proj_dist_x0` 0.451, `H` 0.480, `proj_nonconvex0` 0.491,
`proj_frac_convex_x0` 0.492, `frac_convex_ini` (constant 1.0, no information). Notably: the Eq. (6) projection distance
and the number of reflex corners it creates — the two quantities the earlier reports point
at as the obstruction — carry **no** signal about which designs the K9c solver can then
repair.

## Files

| file | what |
|---|---|
| `Kirigami/apps/kill_yield.jl` | feature extractor, three modes (`k9c`, `fresh`, `check`) |
| `results/yield/features.csv` | 400 rows, features + joined K9c outcomes |
| `results/yield/fresh.csv` | 200 held-out rows, features + new outcomes |
| `results/yield/fresh_sigma/` | the 100 held-out `σ_def` orientations, one JSON per graph |
| `results/yield/analyse.py`, `analysis_out.txt` | the analysis and its full output |
| `results/yield/fig_yield_predictor.png` | top feature vs outcome, both populations |
| `results/yield/fig_yield_matched.png` | family yield at matched \|F\| and density bins |
| `results/yield/fig_yield_auc.png` | AUC table as a bar chart |
