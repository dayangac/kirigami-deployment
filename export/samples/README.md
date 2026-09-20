# export/samples

## `cat_*` — the tech report's teaser design

A 134-face cat-head graph (two hexagonal eyes and a triangular nose are single faces),
built and designed by `Kirigami/scripts/teaser_cat.jl --seed 7`: max-cut orientation, 35
split cuts, shape space of dimension 35. The Eq. (6) projection has `Θ_max = 0`
(split-inward); the range-maximising embedding has exact `Θ_max = 2.5688` rad, certified.

| file | contents |
|---|---|
| `cat_input.json` | the graph at the input embedding, with `orientation` |
| `cat_graph.json` | the same graph at the range-maximising embedding (what `kiri_export` read) |
| `cat_{closed,open_half,open_0.9tm}.{svg,stl,3mf,json,png}` | `kiri_export --profile felt_laser --scale 40` at `θ = 0`, `Θ_max/2`, `0.9 Θ_max`; closed sheet 84 × 102 mm |

All three solids report `closed = true`, 0 boundary and 0 non-manifold edges; the
ear-clipping warnings in the JSON reports are the fan fallback on a few sliver faces and
do not open the solid. At the default 10 mm/unit the hinge necks are too large for the
sliver faces and the solid is not closed, which is why the scale is 40.
