# results/final/figures -- baseline vs. method figure set

Data sources and the exact command that produced everything here:
`julia --project=Kirigami/scripts Kirigami/scripts/plot_final.jl` (reads only; does not
re-run any solver).

| figure | what it shows | data source(s) |
|---|---|---|
| `fig_yield.png` | Yield vs. face count \|F\| in 3 bins, per family (Delaunay/Voronoi/quad) and per sigma rule, 4 arms (baseline Eq.(6)+repairs, K9 proximity, K9b, K9c range-maximising), Wilson 95% CIs, two rows: exact Theta_max>0 and certified eps_max>=0.1 rad. | `results/kill/k9c/shard_*.csv` merged by (kind,id,sigma), N=400 rows at merge time; baseline confirmed identically 0 from `results/kill/k6/k6_final.csv`. |
| `fig_eps_hist.png` | eps_max distribution (log-x) for K9c vs. K9 proximity, plus a bar of the zero-eps_max count per arm. | same K9c merge. |
| `fig_feasible_vs_deployable.png` | Scatter of best 0+ margin (convex+split feasibility slack) vs. exact Theta_max for every row, colored by family, with margin-feasible-but-Theta_max=0 designs ringed. | same K9c merge. |
| `fig_gallery_full.png` | ALL K9c-deployable designs (Theta_max>0), sorted by certified eps_max descending, shown closed and at Theta_max/2, with a caption noting the baseline has no deployed panel to show (Theta_max=0 throughout). | `results/kill/k9c/gallery/*.json` (per-design geometry dumps) + the same merged CSV for sort order. |
| `fig_e1_agreement.png` | E1's \|exact - bisection referee\| histogram (log-y) and the certificate's confusion counts (TP=2150, FP=0, FN=0, TN=963 against Theta_max>=eps=0.006). | `results/final/e1/e1.csv`, N=3113 rows. |
| `summary_table.md` | Headline table: arms x {N, deployable, certified, eps_max>=0.1 rad, median eps_max, max eps_max}, computed directly from the CSVs above. | same. |

## Notes

- K9c's `results/kill/k9c/summary.txt` and `k9c.csv` are a stale 182-row snapshot; this
  script always merges the live `shard_*.csv` files instead. At merge time the 12 shards
  covered N=400 of the 400-design target (200 graphs x 2 sigma). If the run had still been
  short of 400, this note would say so and the script would need re-running once more
  shards land -- it reads the shards live every time, never the stale `k9c.csv`.
- The baseline (Eq. (6) + the K6 0+ repair) is independently confirmed 0/400 designs
  deployable and 0/400 certified in `results/kill/k6/k6_final.csv` before it is asserted
  in any figure here.
- `fig_gallery_full.png` intentionally renders every deployable design, not a curated
  subset (review/constructive_review.md section 5, item 1) -- it is a large image, one
  row per 6 designs.
