# Is the K9c failure set a solver basin, or is it geometry?

Every number below is copied from `results/yield/basin_stats.txt`,
which `Kirigami/apps/exp_basin_agg.jl --nshards 6` derives from `results/yield/basin.csv`
(1,107 rows = 123 designs × 9 seeds). Nothing here is estimated by hand.

## The question

`results/yield/YIELD.md` could not predict K9c's 93 failures from input-side features
(best held-out AUC 0.739 against a 0.85 bar) and left untested whether those failures are
basins of the non-convex two-stage search rather than properties of the geometry. K9b
already answered the analogous question for the **proximity** objective — a 4× budget made
things worse — but the **range** objective had never been re-seeded on the failures. A
referee is entitled to ask "did you just not search hard enough?", and this file answers
exactly that.

## What was run

- Driver `Kirigami/apps/exp_basin.jl`; aggregator `Kirigami/apps/exp_basin_agg.jl`; figure
  `results/yield/plot_basin.py`. Nothing under `Kirigami/src` and nothing under `results/experiments`
  was modified.
- **Population, 123 designs.** The 93 rows of the K9c pass that was current when this
  experiment was specified (`theta_exact ≤ 1e-9`), plus 30 successes drawn deterministically
  as every 10th success in row order (the 10th, 20th, …, 300th of its 307). The canonical
  `results/experiments/k9c/k9c.csv` is a later pass of the same driver and, K9c's endpoints being
  optimiser-path dependent, has 88 failures: 75 of them are in this file's 93, and 18 of
  the 93 deploy in the canonical pass. The `k = 0` column below is the archived answer of the
  pass the population was drawn from, so that the flip statistics are internally
  consistent.
- **Seeds.** `design_range_max` (`Kirigami/src/method/design.jl`) at the K9c run's default options, with
  `seed = 9300 + 7·id + which + 1000·k` for `k = 1…8`. The archived K9c answer is copied
  out of `k9c.csv` and recorded as `k = 0`, so nine seeds sit in one file and `k = 0` is
  never a recomputation.
- **Designs are regenerated exactly as K9c and `YIELD.md` build them**: `make_graph(id, 100, 800,
  1400)`, `sigma_mc = mesh.sigma`, `sigma_def` read from `results/experiments/k5/sigma`,
  `X_ini = mesh.X`.
- **Instrument.** `characterize` — the exact T4.2″ `Θ_max`, the binding
  classification, and the largest certified `ε_max`. The independent bisection referee (at
  shrink 1e-12 through `characterize`, and at 1e-9 through K9c's own routine) is run only
  on points with `Θ_max > 0`.
- **Verification that the solver is the K9c solver**: the four locked regression cases in
  `Kirigami/test/test_design.jl` ("design_range_max reproduces the K9c CSV row …") run
  under `Pkg.test(test_args=["test_design"])` — 254 assertions pass, 0 fail, 23 are
  `@test_broken`: the structural outputs (vertex and face counts, `dim_null`, options) are
  locked exactly, while the `Θ_max`/`ε_max` endpoints of the fixed-budget range optimiser
  are path-dependent locks kept as `@test_broken` (`docs/NUMERICS.md`), which is the same
  path dependence this experiment measures.
- Output: `results/yield/basin.csv` (one row per design × seed, 6-way sharded and merged),
  `results/yield/basin_stats.txt` (every derived number), `results/yield/fig_basin.png`.

## 1. Flip rate of the failures

- **7 of 93 = 7.5 %** of the failures deploy under at least one of the 8 new seeds, Wilson
  95 % [3.7, 14.7] %. Of the 7, **4** reach `ε_max ≥ 0.1` rad at their best seed and all
  **7** are confirmed by the `1e-9` bisection referee.
- By family: Delaunay 3 / 26 (11.5 %), quad-random 2 / 55 (3.6 %), Voronoi 2 / 12 (16.7 %).
  By the original binding contact: inverted 2 / 27, split-inward 4 / 57, vertex-edge 1 / 9.
  By orientation: `σ_def` 5 / 48 (10.4 %), `σ_mc` 2 / 45 (4.4 %).
- The flipped designs (`n` = deploying seeds of 8, `Θ_max` at the best seed):
  `voronoi_39 σ_def` (8, 0.0066), `quad_random_47 σ_def` (8, 0.165),
  `delaunay_40 σ_def` (8, 0.480), `delaunay_151 σ_def` (8, 0.206),
  `delaunay_145 σ_def` (8, 0.148), `voronoi_135 σ_mc` (8, 0.031),
  `quad_random_170 σ_mc` (8, 0.062). Best `Θ_max` over the flipped: min 0.0066, median
  0.148, max 0.480 rad.
- **The distribution of deploying seeds is all-or-nothing: 86 designs at `n = 0`, 7 at
  `n = 8`, none in between.** On every one of the 123 designs the 8 new seeds return the
  *same* `Θ_max` to `1e-9`: the Gaussian restarts that the seed controls never beat the
  deterministic `t = 0` start (provenance `k9c/x0`, `k9c/x0+B`, `k9c/k9b+B` or the deterministic arm-`k9`
  point on every `k ≥ 1` row), so the seed is not a lever at all. What the 7 flips measure is
  therefore the difference between the archived `k = 0` endpoint and the current solver's
  endpoint on the same design — the optimiser-path dependence K9c already documents
  (`results/experiments/EXPERIMENTS.md` §K9c) — not a seed-selected basin. Indeed 4 of the 7 also
  deploy in the canonical `k9c.csv`, at comparable values (0.0066, 0.100, 0.222 and 0.032
  rad); the other 3 are at `Θ_max = 0` there.

## 2. Seed-to-seed spread on the successes

- All 30 sampled successes deploy under all 9 seeds; **none has a non-deploying seed.**
- Relative IQR of `Θ_max` over the 9 seeds: median 0, q90 0, max 0 — again because the 8 new
  seeds coincide. The only spread is between the archived `k = 0` value and the current
  endpoint: best / original `Θ_max` median 1.00, q90 1.25, max 2.35; the archived seed is
  the best of the 9 on 22 of 30 (73 %), the current endpoint wins on 8.
- `ε_max ≥ 0.1` rad on 21 of 30 at `k = 0`, and 21 of 30 at the best of 9.

## 3. Best-of-9 yield on the 400

- Single-seed, the pass the population was drawn from: 307 / 400 = 76.8 %.
- Best-of-9 (exact on the deploy count, since every failure was rerun and a success stays a
  success with `k = 0` among the nine): **314 / 400 = 78.5 %**, Wilson 95 % [74.2, 82.2] %.
  The canonical single-pass `k9c.csv` sits at 312 / 400 = 78.0 %, inside that interval.
- `ε_max ≥ 0.1` rad is **extrapolated** from the 30 sampled successes only (rate 0.70 at
  `k = 0` and at the best of 9): implied 214.9 over the 307 successes plus 4 from the flipped
  failures, an implied 400-design rate of about 55 %. An estimate, not a measurement.

## 4. Verdict — search floor or geometry?

Pre-registered rule, fixed before the run: **flip rate ≥ 20 % → search floor**, and the
paper must then report best-of-9 alongside single-seed; **< 20 % → predominantly
geometric**.

**7.5 % < 20 % → predominantly geometric: the failure set is not a solver basin.** The
finding is stronger than the rule anticipated, because the eight extra seeds never produce a
point the deterministic start does not: 86 of the 93 failures stay at `Θ_max = 0` under
every seed, and the 7 that move are moved by floating-point path, not by search. Re-seeding
the range optimiser is not a way to raise the K9c yield; the paper reports the single-pass
count and notes its ±2 % path dependence.

## 5. What distinguishes the flipped designs (only if the rule fires)

Not run — the rule did not fire.
