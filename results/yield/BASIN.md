# WP2b — is the K9c failure set a solver basin, or is it geometry?

**STATUS: RUN IN PROGRESS.** This file is written incrementally; every number below is
filled in from `results/yield/basin_stats.txt` once the shards finish, and the sections
that are not yet filled say so explicitly. Nothing here is estimated by hand.

## The question

`results/yield/YIELD.md` (WP2) could not predict K9c's 93 failures from input-side features
(best held-out AUC 0.739 against a 0.85 bar) and left untested whether those failures are
basins of the non-convex two-stage search rather than properties of the geometry. K9b
already answered the analogous question for the **proximity** objective — a 4× budget made
things worse — but the **range** objective had never been re-seeded on the failures. A
referee is entitled to ask "did you just not search hard enough?", and this file answers
exactly that.

## What was run

- Driver `Kirigami/apps/kill_basin.jl`; aggregator `Kirigami/apps/kill_basin_agg.jl`; figure
  `results/yield/plot_basin.py`. Nothing under `Kirigami/src` and nothing under `results/kill`
  was modified.
- **Population, 123 designs.** All 93 rows of `results/kill/k9c/k9c.csv` with
  `theta_exact ≤ 1e-9`, plus 30 successes drawn deterministically as every 10th success in
  row order (the 10th, 20th, …, 300th of the 307).
- **Seeds.** `method::design_range_max` at the K9c run's default options, with
  `seed = 9300 + 7·id + which + 1000·k` for `k = 1…8`. The archived K9c answer is copied
  out of `k9c.csv` and recorded as `k = 0`, so nine seeds sit in one file and `k = 0` is
  never a recomputation.
- **Designs are regenerated exactly as K9c and WP2 build them**: `make_graph(id, 100, 800,
  1400)`, `sigma_mc = mesh.sigma`, `sigma_def` read from `results/kill/k5/sigma`,
  `X_ini = mesh.X`.
- **Instrument.** `method::characterize` — the exact T4.2″ `Θ_max`, the binding
  classification, and the largest certified `ε_max`. The independent bisection referee (at
  shrink 1e-12 through `characterize`, and at 1e-9 through K9c's own routine) is run only
  on points with `Θ_max > 0`.
- **Verification that the solver is the K9c solver**: the four locked regression cases in
  `Kirigami/test/test_design.jl` ("design_range_max reproduces the K9c CSV row …") pass in
  this build tree — 6 test cases, 350 assertions, 0 failures.
- Output: `results/yield/basin.csv` (one row per design × seed, 3-way sharded and merged),
  `results/yield/basin_stats.txt` (every derived number), `results/yield/fig_basin.png`.

## 1. Flip rate of the failures

_Pending: filled from `basin_stats.txt` §1 when the run completes._

## 2. Seed-to-seed spread on the successes

_Pending: filled from `basin_stats.txt` §2 when the run completes._

## 3. Best-of-9 yield on the 400

_Pending: filled from `basin_stats.txt` §3 when the run completes._

## 4. Verdict — search floor or geometry?

Pre-registered rule, fixed before the run: **flip rate ≥ 20 % → search floor**, and the
paper must then report best-of-9 alongside single-seed; **< 20 % → predominantly
geometric**.

_Pending._

## 5. What distinguishes the flipped designs (only if the rule fires)

_Pending; runs only when the flip rate is ≥ 20 %._
