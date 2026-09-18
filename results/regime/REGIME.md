# WP7a — the paper-regime population

**Status: COMPLETE.** Everything here is a function of `results/regime/regime.csv` (213
design rows), via `results/regime/summary.txt`, which
`julia --project=Kirigami Kirigami/apps/kill_regime.jl --mode aggregate --out results/regime`
derives from that CSV and from nothing held in memory. The native cells are the authors'
unchanged binary (`baseline/native/`), run with a 600 s cap per cell.

## What was run

Three sets, each under the Eq. (6) baseline (`design_baseline`) and the K9c range-maximising
embedding (`design_range_max`), plus the authors' full native pipeline where it could be
dispatched (their own colouring, and our σ handed to their pipeline):

- **pop** — 100 random graphs, `make_graph(id, 20, 100, 1400)`, ids 2000–2099, so that |F|
  sits in the paper's own regime: min 24, median 52, max 100 (Delaunay 33 graphs, median 32;
  quad-random 34, median 84; Voronoi 33, median 48). Both orientation rules: 200 designs.
- **ref** — the 8 authored reference tilings with their given σ.
- **pat2025** — the 5 fabrication SVGs of the 2025 paper that load as polygon soups
  (`patterns_load.csv`: 16 further files are line-group specs with no geometry).

## Measured — `results/regime/summary.txt`

Random population (200 designs), deployable = exact `Θ_max > 0`, refereed by the `1e-9`
bisection:

| arm | N | deployable | refereed | certified (`ε_max > 0`) | `ε_max ≥ 0.1` rad | median `ε_max` | max |
|---|--:|--:|--:|--:|--:|--:|--:|
| Eq. (6) baseline | 200 | 1 | 1 | 1 | 1 | 0.560 | 0.560 |
| native pipeline, their colouring | 52 dispatched | 12 | 12 | 12 | 12 | 0.290 | 1.109 |
| native pipeline, our σ | 45 dispatched | 5 | 5 | 4 | 4 | 0.213 | 0.646 |
| **K9c range-max** | 200 | **192** | **192** | **192** | **178** | **0.858** | π |

Yield against |F| with Wilson 95 % intervals: in [20, 50) the baseline gives 1/98 (1.0 %),
the native pipeline 12/30 (40.0 %, [24.6, 57.7] %) with their colouring and 5/27 (18.5 %)
with ours, K9c 96/98 (98.0 %, [92.9, 99.4] %); in [50, 100] the baseline 0/102, the native
pipeline 0/22 and 0/18, K9c 96/102 (94.1 %, [87.8, 97.3] %).

Native coverage: of the 97 native cells dispatched on the random population, 86 completed and
11 timed out at 600 s (`native cell status` table in `summary.txt`); the pattern and
reference sets were not dispatched to the native pipeline. Every native positive is a graph
of 27–31 faces (28–31 with their colouring, 27–29 with our σ).

Reference tilings: all 8 deploy and certify `ε_max ≥ 0.1` rad under both arms; 5 of them
reach `Θ_max = π` under both (the three arms are identical points there), the other three
open to 2.09, 2.36 and 1.65 rad under the baseline, and K9c lifts snub_square_33434 from
1.65 rad to π and tiling_3_4_3_12 from 2.38 to 2.62 rad. The 5 imported 2025
patterns have an exact `Θ_max` of 2.89–3.14 rad under both arms (identical points: the
range optimiser accepts the projection as is), but a certified `ε_max` of at most 0.067 rad,
because the imported fabrication geometry carries near-grazing contacts that the certificate
counts as roots.

## Reading

In the paper's own regime (|F| ≤ 100) the picture is the K9 population's, shifted up: the
published projection deploys 1 of 200 random designs, the authors' full pipeline 12 of the 52
designs it was dispatched on with their colouring (48 completed, 4 timed out; every positive
has 28–31 faces), and the range
objective 192 of 200 with a median certified opening of 0.86 rad. The K9c endpoints are
optimiser-path dependent (`results/kill/KILL_REPORT.md` §K9c); the baseline and native
columns are not. Figure: `results/regime/fig_regime.png`.
