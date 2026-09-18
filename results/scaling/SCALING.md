# WP7b — the scaling benchmark

**Status: COMPLETE.** Every number is a function of `results/scaling/scaling.csv` (165
rows, one per graph × routine) and of `results/scaling/fits.csv` (the least-squares
exponents written by `Kirigami/scripts/plot_scaling.jl`, which also draws
`results/scaling/fig_scaling.png`); `results/scaling/capped.csv` lists the cells the wall
cap killed.

## What was run

`Kirigami/apps/kill_scaling.jl` (driven by `Kirigami/scripts/run_scaling.sh`): random
Delaunay and Voronoi graphs at site counts landing near target |F| ∈ {50, 100, 200, 500,
1000, 2000, 5000}, three seeds each, the kill population's orientation rule, one process per
(graph, routine) so that the peak resident set is the routine's own. The measured |F| runs
from 40 to 5,000. Four routines per graph, each under a 1,800 s wall cap:

| routine | what it times |
|---|---|
| `solve_dense` | assemble Eq. (4) and solve it by a dense SVD: rank, `dim_null`, the Eq. (6) projection |
| `rank_sparse` | the rank alone, by sparse QR (no null-space basis, so nothing downstream can run from it) |
| `characterize` | the exact T4.2″ `Θ_max` scan plus the `ε = 0.3` certificate at `X0` |
| `design_range_max` | the K9c range-maximising embedding at run defaults, **including** its own dense solve |

## Measured

Least-squares exponents of `secs` against |F| (log–log, 95 % CI), `results/scaling/fits.csv`:

| routine | exponent | 95 % CI | R² | n | largest |F| completed | max wall over its rows (s) | peak RSS (MB) |
|---|--:|--:|--:|--:|--:|--:|--:|
| `solve_dense` | **1.86** | [1.57, 2.15] | 0.80 | 42 | 5,000 | 158 | 4,256 |
| `rank_sparse` | **0.21** | [0.12, 0.29] | 0.37 | 42 | 5,000 | 0.45 | 1,824 |
| `characterize` | **0.85** | [0.74, 0.95] | 0.86 | 42 | 5,000 | 2.1 | 1,144 |
| `design_range_max` | **1.49** | [1.25, 1.73] | 0.80 | 39 | 4,974 | 1,048 | 1,848 |

(The max-wall column is the slowest row of each routine, not the wall at the largest |F|:
for `design_range_max` it is the seed-2 Voronoi cell at |F| = 2,000; the largest completed
Delaunay cell, |F| = 4,974, took 158 s.)

Per cell (median over the three seeds): the exact characterization costs 0.02 s at |F| = 41
and 2.05 s at |F| = 5,000; the dense solve 0.005 s and 128 s; the range optimiser 0.4 s at
|F| = 41, 22 s at |F| ≈ 2,000 (Delaunay), 792 s at |F| = 2,000 (Voronoi), and 158 s at
|F| ≈ 5,000 on the Delaunay side, where it completed all three seeds. The three Voronoi
cells at 5,000 sites are the only capped cells (`capped.csv`: the range optimiser exceeded
the 1,800 s cap on all three seeds); every other (graph, routine) cell completed.

Deployability along the way (`theta_max > 0` after `design_range_max`, 3 seeds per cell):
3/3 at |F| ≤ 185 on both families; 2/3 or 3/3 up to |F| = 2,000 (Delaunay 2/3 at 980 and
1,971, Voronoi 2/3 at 200 and 500, 3/3 at 1,000, 2/3 at 2,000); 1/3 at |F| ≈ 5,000
(Delaunay). The `characterize` rows report `Θ_max = 0` on every graph, since they are
evaluated at the raw Eq. (6) projection `X0` — the K1a/K6 result at every scale measured.
The optimiser's endpoints are path dependent (`results/kill/KILL_REPORT.md` §K9c), so the
per-cell deployable counts are indicative; the wall-time exponents are not sensitive to
them.

## Reading

The exact range and its certificate are cheap and scale sub-linearly in the face count
(exponent 0.85, 2 s at 5,000 faces): the closed-form calculus is not the bottleneck at any
size in the paper's regime or an order of magnitude above it. The dense null-space solve is
the memory-bound step (4.3 GB at 5,000 faces, exponent 1.9 in time), and the range optimiser
inherits that solve plus its own iteration (exponent 1.5, 13 min median and 17 min worst at
2,000 Voronoi faces).
The sparse rank is essentially free but yields no basis, so it cannot replace the dense solve
for design. Figure: `results/scaling/fig_scaling.png`.
