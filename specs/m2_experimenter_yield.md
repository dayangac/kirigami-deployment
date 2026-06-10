# Mission 2 / WP2 — Experimenter-Y: explain the K9c yield trend

Paste of specs/common_preamble.md applies (read it first).

## Why this exists
REPORT.md §Results says the K9c method's yield "is not scale-free": exact `Θ_max > 0` on 307/400 designs overall, but Delaunay falls from about 1.0 to 0.8 and quad-random from 0.75 to 0.43 across `|F|` bins (results/final/figures/fig_yield.png, code/scripts/plot_final.py `fig_yield`, bins [100, 230, 420, 800]); per family Voronoi 122/134, Delaunay 108/134, quad-random 77/132. REPORT.md "Next steps" item 7: "Whether that is split-cut density, aspect ratio, or the solver's basin is unknown, and it is the first question a reviewer will ask." Your job is to answer it with a measured, out-of-sample-validated predictor and a mechanism sentence.

## Read first
- results/kill/k9c/k9c.csv (400 rows; header: id,kind,sigma,N,F,med_edge,n_split,dim_null,n_corner,nonconvex0,k9_*,k9b_*,k9c_feas,k9c_margin,k9c_theta,k9c_eps,k9c_stage_b,best_src,best_feas,best_margin,best_min_cross,best_min_q,best_min_mu,eps_max,theta_exact,theta_ref12,theta_ref9,binding,dist,secs) and results/kill/k9c/summary.txt.
- results/kill/KILL_REPORT.md §K9c (lines ~2188–2430) — what the arms are, the `binding` classification (split-inward 57 / inverted 27 / vertex-edge 9 of the 93 failures), the "yield is not scale-free" paragraph.
- code/apps/kill_common.hpp (`make_graph(id, min_faces, max_faces, n_cap)`: population is a deterministic function of id; kind by id % 3; the K9 population is ids 0..199 with min 100, max 800, n_cap per kill_k9c.cpp — read that file's call to confirm) and code/apps/kill_k9c.cpp (how σ_mc and σ_def are produced per graph; the results/kill/k5 sigma cache it reads).
- code/src/method/design.hpp — `design_range_max` (the K9c pipeline as a library call; defaults reproduce the run), `characterize`, `orientation_maxcut`, `orientation_defect`.
- code/src/core/{cut.hpp, holes.hpp, tutte_auxetic.hpp} for the input-side quantities; code/src/method/zero_plus.hpp (`zero_plus_q`, `zero_plus_form`) and convex_embed.hpp for the 0⁺ constraint counts.
- code/scripts/plot_final.py for the existing figure style and the Wilson-interval code.

## Task
1. **Input-side features, per design** (computed at `X_ini` and from the cut structure only — nothing that the solver produces). Write code/apps/kill_yield.cpp that, for each of the 400 K9c designs (same ids, same σ_mc / σ_def as kill_k9c.cpp), regenerates the graph and σ and emits results/yield/features.csv with at least:
   - `F`, `N`, `E`, `n_split`, `n_hinge`, `split_density = n_split / E_interior`, `dim_null`;
   - `n_constraints_0plus` = number of split-sign constraints + number of corner incidences (the count `range_embed` actually enforces), and `dof_ratio = 2·dim_null / n_constraints_0plus`;
   - face geometry at `X_ini`: min / median corner angle, min / median face aspect ratio (circumradius over inradius, or longest over shortest edge), fraction of faces with a corner below 20°;
   - cut structure: number of split-forest components, size of the largest split component, split-forest depth, hinge-graph diameter, fraction of interior vertices touched by no split edge ("pure" vertices);
   - hole structure: `H`, mean and max preimage size, number of notches (boundary-touching components);
   - the 0⁺ signs **at X_ini** (not at X0): fraction of split edges with `q_e(X_ini) > 0`, and `min q_e(X_ini) / med²`; likewise fraction of convex corners at `X_ini` (should be 1.0 by K1a — verify and report);
   - the Eq. (6) projection distance `‖X0 − X_ini‖` per vertex in median edges, and `nonconvex0` (already in k9c.csv) — these are *projection*-side, keep them in a separate column group labelled as such.
   Join with k9c.csv on (id, kind, sigma) to attach the outcomes `theta_exact > 0`, `eps_max ≥ 0.1`, `binding`.
2. **Fresh graphs for held-out validation.** Generate 100 new graphs with `make_graph(id, 100, 800, n_cap)` for ids 1000..1099 (outside the K9 population; confirm they are not equal to any id 0..199 graph), both σ rules, run `design_range_max` with the run defaults (seed 9300 + 7·id + which, as kill_k9c.cpp does), characterize, and record the same features and outcomes to results/yield/fresh.csv (200 rows). Runtime: the K9c run took ~50 s per design on a loaded machine; budget accordingly, shard 4-way, and write incrementally. If the full 200 does not fit in ~3 h, stop at ≥ 120 rows and say so.
3. **Analysis** (Python allowed here: matplotlib for figures; use numpy for the logistic fit and AUC — implement Wilson intervals and AUC yourself or with numpy, do not assume scikit-learn is installed; it is not). Write results/yield/analyse.py:
   - univariate: for every feature, AUC for `theta_exact > 0` on the 400, and Spearman correlation with `eps_max`; table sorted by AUC.
   - the size/family confound: yield per family **at matched `|F|` bins** (the plot_final.py bins) and the same for `split_density` bins, so the reader sees whether family differences survive matching.
   - a logistic model fit on the σ_mc half (200), evaluated on the σ_def half (200) and on fresh.csv; then the reverse fit. Report held-out AUC for the best single feature and for the best two-feature model; choose by fit-set AUC, report held-out AUC — never choose on the held-out set.
   - the binding-mode breakdown of the 93 failures against the top feature (does split-inward binding sit at high split density?).
   - `theta_exact` and `eps_max` vs the top feature, scatter with the fresh graphs overplotted in a second colour.
   Figures: results/yield/fig_yield_predictor.png (top feature vs outcome, both populations), fig_yield_matched.png (family yield at matched bins), fig_yield_auc.png (AUC table as a bar chart).
4. **Report** results/yield/YIELD.md: the feature table with AUCs; the held-out numbers; the matched-bin table; one paragraph naming the mechanism the winning predictor supports (e.g. "yield falls because the number of split-sign constraints grows faster than the null-space dimension" — only if the numbers say so); the PASS/FAIL verdict against the bar below; and, honestly, the features that did **not** predict. Every number cites features.csv / fresh.csv / analyse.py output.

## PASS bar (pre-registered)
Held-out AUC ≥ 0.85 for `theta_exact > 0` on the fresh graphs, for one feature or a two-feature model, **and** the model's fresh-graph yield prediction within the Wilson 95 % interval of the observed fresh yield. Anything else is FAIL, reported with the best held-out AUC obtained. A predictor that uses solver outputs (`k9c_margin`, `best_*`, `dist`, `secs`, `binding`) is a restatement and is excluded from the bar; you may still report such correlations in a clearly labelled "solver-side" table.

## Rules
- Do not modify code/src/method or code/src/core. New code goes in code/apps/kill_yield.cpp (add to CMakeLists.txt's app list) and results/yield/analyse.py.
- Do not rerun K9c on the 400; read its CSV. The fresh-graph run is the only new solver work.
- Write features.csv, fresh.csv, YIELD.md incrementally. Keep concurrency ≤ 4 processes; another long job (Native200 rerun) shares the machine.
- Do not commit; do not edit REPORT.md/IDEA.md/STATE.md.

## Reply
≤15 lines: rows in features.csv and fresh.csv; observed fresh yield with Wilson interval; top three features by held-out AUC with their numbers; PASS/FAIL; the mechanism sentence; anything not done.
