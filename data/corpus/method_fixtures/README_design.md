# method_fixtures -- frozen intermediates of `method/design.cpp` (K9 / K9c rows)

Produced 2026-09-19 by `freeze_design.cpp` and `freeze_k9c.cpp` (this directory), built
against `kirigami-experiments/code/build/libkiri_core.a` with the corpus toolchain
(Apple clang 16, `-O2 -std=c++20 -arch arm64`, default `-ffp-contract=on`, Eigen 5.0.1,
nlohmann-json 3.12):

    clang++ -std=c++20 -O2 -arch arm64 -I code/src -I code/apps -I /opt/homebrew/include \
        -I /opt/homebrew/include/eigen3 freeze_<x>.cpp code/build/libkiri_core.a -o freeze_<x>
    ./freeze_design 148 design_intermediates_148.json
    ./freeze_k9c <id> k9c_calls_<id>.json        # id in 130 148 30 42

All meshes come from `kill::make_graph(id, 100, 800, 1400)` (= `data/corpus/k1a_200.json`
row `id`, bit-identical), `X_ini = mesh.X`, orientation `sigma_mc`.

| file | what |
|---|---|
| `design_intermediates_148.json` | K9 row id 148, `DesignOptions` defaults, seed `9000 + 7*148`: `med_edge`, `dim_null`, `X0`, `Phi` (row-major N x k), the three solver calls of `design_constrained` in order (`ra` = convexity-only warm start, `sp` = split-only repair, `rb` = variant (b)) with their `t`, `X` and reports, the final `design` (X, t, margins, characterization) and `characterize` at the design X, at X0 and at X_ini. |
| `k9c_calls_<id>.json` | K9c rows (seed `9300 + 7*id`, `RangeMaxOptions` defaults): `X0`, `Phi`, arm `r9` and `r9b` (`convex_embed` results), `stage_a[]` = every `range_embed` call with its `t_init`, `drel`, result, exact `theta_max` and `margin_exact`; `stage_b` = `maximize_margin_range` replayed from the stage-A winner's X; `result` and `arms` of `design_range_max`. |

`Phi` is Eigen's null-space basis; the Julia `solve_system` spans the same subspace
(projector difference 1e-15) in a different basis, so feed THIS `Phi` (and the recorded
`t_init`) to a Julia optimiser call when comparing it call by call.

Vertex / face indices inside these files are 0-based where the C++ wrote them (none of
the fields above are index lists; `X`, `X0`, `t` are plain coordinate arrays).

## What the comparison showed (2026-09-19, Kirigami.jl at that date)

* `characterize` at the recorded C++ points is bit-identical (theta_max, eps_max,
  contact count, pair count, binding) on every fixture.
* `convex_embed` / `zero_plus_repair` on identical inputs agree to converged-solver
  tolerance (max |dX| 1.6e-10 .. 1.4e-6, same feasibility and best_start, iteration
  counts within a few %).
* `range_embed` (stage A: 6 x 120 fixed L-BFGS iterations, never converged) is
  reproduced to PATH level, not bit-exactly (port-method-2's stated intent): objective,
  gradient, modes and margin_start agree to 1e-14 at every probe, the iterates agree to
  1e-15 for the first ~10 iterations and then diverge by rounding amplification
  (1e-12 / 1e-8 / 1e-5 after 20 / 40 / 60 iterations), reaching max |dX| 1.8e-2 .. 0.26
  after 714 iterations, hence a different exact theta at the output (hero2 t = 0 start
  0.2887 vs 0.2790). Bit-exactness would need fma placement matched in the objective,
  the gradient and core/optimize.jl's L-BFGS. This is why the K9c CSV rows are not
  reproduced to 1e-9 by the Julia port (`@test_broken` in `Kirigami/test/test_design.jl`).
* `maximize_margin_range` from the same X: identical `improved`, caps and theta_after.
