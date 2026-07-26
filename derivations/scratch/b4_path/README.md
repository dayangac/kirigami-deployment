# b4_path -- why `check_b4_93` does not replay b4.csv: the optimiser path, not the port

Investigated 2026-09-20 (port-method-2). `probe_b4.cpp` (build against
`code/build/libkiri_core.a` with the corpus toolchain: Xcode clang, `-O2 -std=c++20 -arch arm64`,
`-I code/src -I code/apps`) and `probe_b4.jl` make kill_b4's exact calls on one graph id with the
archived `free_sigma_def` shape (`results/kill/b4/cache/free_sigma_def/shape_<id>.bin`, so X0 AND
Phi are bit-identical on both sides, not just the null space), and print every field of the
secondary (split-only) and primary (split + vertex-edge + proximity) `zero_plus_repair` results.
`trace_b4.jl` and `probe_b4 <id> <max_iter>` run the secondary repair at increasing iteration caps.

```
./probe_b4 93            # from ~/Documents/kirigami-experiments
./probe_b4 93 <iters>
julia --project=Kirigami derivations/scratch/b4_path/probe_b4.jl 93
julia --project=Kirigami derivations/scratch/b4_path/trace_b4.jl 93
```

## id 93, seed 6652, n_random 3, lambda_rel 1e-6 (the B4 options)

Start state, bit-identical: `f_start` 874.32613327176898 (C++) / 874.32613327176887 (Julia);
`min_q_start` -12.69905392006687 both; 275 inward split edges both; 5206 corner incidences both;
med 1.0291298269997151 both.

Secondary repair at iteration cap `it` (best of 4 starts, `f_end` = F at the returned point):

| it | best_start C++ / Jl | f_end C++ | f_end Julia | rel. diff |
|---|---|---|---|---|
| 1 | 0 / 0 | 784.08291749720274 | 784.08291749720365 | 1e-15 |
| 2 | 0 / 0 | 699.74914777431968 | 699.74914777432014 | 7e-16 |
| 3 | 0 / 0 | 484.56630303834748 | 484.56630303834697 | 1e-15 |
| 5 | 0 / 0 | 349.75716509459545 | 349.75716509459568 | 7e-16 |
| 10 | 0 / 0 | 182.66394236457253 | 182.66394236457251 | 1e-16 |
| 20 | 0 / 0 | 116.01352503941179 | 116.01352503941547 | 3e-14 |
| 40 | 1 / 1 | 82.505651793019652 | 82.505651701932621 | 1e-9 |
| 80 | 0 / 0 | 44.145523863706643 | 44.162927558953506 | 4e-4 |
| 160 | 0 / 0 | 24.270304245801121 | 23.832249255100368 | 2e-2 |
| 320 | 0 / 0 | 4.8057506691038077 | 4.9250888296344968 | 2e-2 |
| 640 | 0 / 0 | 0.33872123723383518 | 0.37483278905061623 | 1e-1 |
| 1199 (B4 cap, not converged) | 0 / 0 | 0.25074497545363716, \|t\| 495.42 | 0.2628781911250177, \|t\| 504.42 | -- |

Primary repair, warm-started from that secondary point:

| | C++ (= b4.csv row) | Julia |
|---|---|---|
| feasible | 1 | 0 |
| best_start / iters | 0 / 150 | 3 / 195 |
| min_q | 0.20245660165651105 | 0.23282833991157695 |
| min_area | 0.124206083212961 | 0.14775075468935484 |
| min_margin / n_bad_margin | 0.16840917246835591 / 0 | -0.30642026654058152 / 1 (of 5206) |

id 96 shows the same profile: `f_end` identical to 1e-16 at it = 10 (117.40787386598787 /
...788, both best_start 2), 3e-12 at it = 40, diverged by it = 160 (7.22 vs 5.61).

## Verdict

Identical start, identical first ~20 L-BFGS iterates, then rounding-path divergence
(1e-14 -> 1e-9 -> 4e-4 at 20 / 40 / 80 iterations) in a 2434-dimensional run that the B4
experiment caps at 1199 iterations without convergence. The secondary end point differs by ~2 %
in norm, the primary is warm-started from it and lands in a different local minimum, and
feasibility flips on exactly one corner margin out of 5206. This is not a port bug: objective,
gradient and every exact quantity agree to 1e-14 (see `data/corpus/method_fixtures/test_method_2.json`),
and bit-exactness over the path would require matching the fma contraction of the objective and the
reduction order of every dot/norm inside the L-BFGS (Eigen vectorised reductions vs Julia BLAS).
Consequence for the paper: B4's single refereed positive (`voronoi_96`, 0.241884 rad) and the
`voronoi_93` row are knife-edge outcomes of a non-converged optimiser and must not be quoted as
reproducible numbers; what reproduces is the start state and the verdict-level behaviour of the
test cases.
