# Mission 2 / WP2b — Experimenter-B: is the K9c failure set a solver basin or geometry?

Paste of specs/common_preamble.md applies (read it first). Theory/computation only (D15).
Environment: build with `/Applications/Xcode.app/Contents/Developer/usr/bin/make -C code/build -j2 <target>` after `cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64` (bare `make`/`cmake --build` are broken by an xcode-select shim). Or build in your own tree `code/build-basin` as WP2 did (`cmake -S code -B code/build-basin ...`) to avoid contention with other agents. The machine runs ~14 baseline processes; keep your concurrency ≤ 3.

## Why
results/yield/YIELD.md (WP2) could not predict the 93 K9c failures from input-side features (best held-out AUC 0.739) and explicitly left untested whether those failures are basins of the non-convex two-stage search. A referee will ask "did you just not search hard enough?" — K9b already showed that a 4× budget under the proximity objective made things worse, but nobody has re-seeded the range objective on the failures.

## Read first
- results/yield/YIELD.md §Mechanism and §Verdict; results/kill/k9c/k9c.csv (columns; `theta_exact`, `binding`, `best_src`); results/kill/KILL_REPORT.md §K9c "Provenance of the winning point".
- code/src/method/design.hpp `design_range_max` / `RangeMaxOptions` (seed = 9300 + 7·id + which in the run; arms; stage settings 6 × 120) and code/apps/kill_k9c.cpp (how graph, σ_mc, σ_def and X_ini are produced per design — reuse exactly; σ_def comes from results/kill/k5/sigma).
- code/apps/kill_yield.cpp (WP2's driver regenerates the same designs; copy its design-regeneration code path).

## Task
1. Driver code/apps/kill_basin.cpp. Designs: all 93 K9c rows with `theta_exact == 0` (or ≤ 1e-9), plus 30 successes drawn deterministically (every 10th success by row order). For each design run `design_range_max` with the run defaults but 8 new seeds: `seed = 9300 + 7·id + which + 1000·k`, k = 1..8. Also record the original run's value from k9c.csv as k = 0. Per (design, k): exact `Θ_max`, `ε_max`, certified, provenance, margin, `binding`, seconds; referee (`characterize` with `referee = true`) on every point with `Θ_max > 0`. Write results/yield/basin.csv incrementally, 3-way sharded.
2. results/yield/BASIN.md: flip rate of failures (designs with ≥ 1 of 8 seeds deploying), by family and by original `binding`; distribution of the number of deploying seeds; for successes, seed-to-seed spread of `Θ_max` (median, IQR, min/max per design) and how often the original seed was the best; the total "best of 9 seeds" yield on the 400 (extrapolate only from the sampled successes, and say so); a paragraph answering: is the 77 % a search floor or geometric? Rule: flip rate ≥ 20 % → search floor, the paper must report best-of-9 alongside single-seed; < 20 % → predominantly geometric. Figure results/yield/fig_basin.png: per-design `Θ_max` over seeds for the 93 failures (sorted), successes in a second panel.
3. If flip rate ≥ 20 %: additionally report what distinguishes the flipped designs from the never-deploying ones using WP2's features.csv (join on id, sigma) — univariate AUC only, no new fitting.

## Rules
- Do not modify code/src. Do not touch results/kill/*. Do not commit; do not edit REPORT.md/IDEA.md/STATE.md.
- Concurrency ≤ 3; budget ≈ 123 designs × 8 seeds × ~50 s / 3 ≈ 4.5 h; write incrementally and report partial state if you must stop early.

## Reply
≤15 lines: designs run × seeds; flip rate overall / by family / by binding; success spread; best-of-9 extrapolated yield; verdict sentence; anything not done.
