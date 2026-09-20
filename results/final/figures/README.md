# results/final/figures -- baseline vs. method figure set

Data sources and the exact command that produced everything here:
`julia --project=Kirigami/scripts Kirigami/scripts/plot_final.jl` (reads only; does not
re-run any solver).

| figure | what it shows | data source(s) |
|---|---|---|
| `fig_yield.png` | Yield vs. face count \|F\| in 3 bins, per family (Delaunay/Voronoi/quad) and per sigma rule, 4 arms (baseline Eq.(6)+repairs, K9 proximity, K9b, K9c range-maximising), Wilson 95% CIs, two rows: exact Theta_max>0 and certified eps_max>=0.1 rad. | `results/experiments/k9c/shard_*.csv` merged by (kind,id,sigma), N=400 rows at merge time; the baseline arm is the best of K6's four repair variants per design from `results/experiments/k6/k6.csv` (17 of 400 deploy). |
| `fig_eps_hist.png` | eps_max distribution (log-x) for K9c vs. K9 proximity, plus a bar of the zero-eps_max count per arm. | same K9c merge. |
| `fig_feasible_vs_deployable.png` | Scatter of best 0+ margin (convex+split feasibility slack) vs. exact Theta_max for every row, colored by family, with margin-feasible-but-Theta_max=0 designs ringed. | same K9c merge. |
| `fig_gallery_full.png` | ALL K9c-deployable designs (Theta_max>0), sorted by certified eps_max descending, shown closed and at Theta_max/2; the baseline's 17 deployable designs are not rendered. | `results/experiments/k9c/gallery/*.json` (per-design geometry dumps) + the same merged CSV for sort order. |
| `fig_e1_agreement.png` | E1's \|exact - bisection referee\| histogram (log-y) and the certificate's confusion counts (TP=2150, FP=0, FN=0, TN=1393 against Theta_max>=eps=0.006). | `results/final/e1/e1.csv`, N=3543 rows. |
| `summary_table.md` | Headline table: arms x {N, deployable, certified, eps_max>=0.1 rad, median eps_max, max eps_max}, computed directly from the CSVs above. | same. |

## Notes

- This script merges the live `shard_*.csv` files rather than reading `k9c.csv`; the two
  agree row for row (N=400 of the 400-design target, 200 graphs x 2 sigma).
- The baseline arm (Eq. (6) + the K6 repairs) is read from `results/experiments/k6/k6.csv` as the
  best of the four repair variants per design: 17 of 400 deploy, all by leaving the
  projection's neighbourhood (EXPERIMENTS.md section K6).
- The K9c arm's endpoints are optimiser-path dependent at the +-2 % level
  (EXPERIMENTS.md section K9c); every count here is that of the canonical CSVs.
- `fig_gallery_full.png` intentionally renders every deployable design, not a curated
  subset (review/constructive_review.md section 5, item 1) -- it is a large image, one
  row per 6 designs.

## Teaser (tech report head figure)

`fig_teaser.{png,pdf}`, `fig_teaser_wide.{png,pdf}` (larger type, for the two-column abstract) — a 134-face cat-head planar graph (eyes and nose are single faces) (`Kirigami/scripts/teaser_cat.jl --seed 7`,
geometry in `teaser_cat.json`, graphs in `cat_input.json` / `cat_graph.json`, fabrication
export in `export/samples/cat_*`): the max-cut orientation, the Eq. (6) projection opened to
`θ = 0.35` with the 22 penetrating faces in red (`Θ_max = 0`, split-inward), the
range-maximising embedding flat, and deployed at `θ = 0.70` (`Θ_max = 2.57` rad, certified).
Rendered by `Kirigami/scripts/plot_teaser.py`.

## Headline figures for the papers

`Kirigami/scripts/plot_headline.py` draws, with large type for column width:

- `fig_yield_all.{png,pdf}` — fraction with exact `Θ_max > 0` against `|F|` from 24 to 793
  faces: the paper-regime population (`results/regime/regime.csv`, 200 designs) and the K9c
  population (`results/experiments/k9c/k9c.csv`, 400 designs), Wilson 95 % intervals; the
  least-norm projection, the published pipeline with its own colouring where it completed
  (52 cells, all at the small sizes), and the range-maximising embedding, pooled and per family.
- `fig_eps_all.{png,pdf}` — certified `ε_max` over the 400 K9c designs against K9, with the
  medians (0.160 / 0.249 rad) and the teaser design marked.
