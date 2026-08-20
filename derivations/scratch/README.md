# derivations/scratch -- the Checker's and Deriver's scratch programs

Each `check_X.cpp` here is the C++ program the derivation files quote; `check_X.jl` next to it is
its Julia port (same corpus, same seeds, same CLI arguments, same printed lines). The C++ sources
stay as documentary evidence (`derivations/check.md`, `core.md`, `lemmas.md` cite them by name)
and are not built. `ideas/rigidity_rig_check.{cpp,jl}` is the same for the rigidity ideator's check.

Run from the repo root:

```
julia --project=Kirigami derivations/scratch/check_t1_t2.jl
julia --project=Kirigami derivations/scratch/check_t4_t5.jl
julia --project=Kirigami derivations/scratch/check_r2.jl [eps=0.02]
julia --project=Kirigami derivations/scratch/check_r3.jl [eps=0.02]
julia --project=Kirigami derivations/scratch/check_b1.jl
julia --project=Kirigami derivations/scratch/check_t5_adm.jl
julia --project=Kirigami derivations/scratch/check_l1.jl [NRAND=52] [NSAMP=20]
julia --project=Kirigami derivations/scratch/check_l2.jl [NEXTRA=100]
julia --project=Kirigami derivations/scratch/check_b4_93.jl [id ...]          # default 93 96
julia --project=Kirigami ideas/rigidity_rig_check.jl
```

`corpus_common.jl` is the slice of `code/apps/kill_common.hpp` that `check_l1` / `check_b4_93` need
(`make_graph`, `reference_cases`), reading the frozen populations in `data/corpus/`
(`derivation_l1_140.json` for `make_graph(id, 18, 46, 220)`, `k1a_200.json` for `(100, 800, 1400)`)
instead of regenerating them. `TODO(julia-port)`: replace by `Kirigami/apps/kill_common.jl`.

## What can and cannot be reproduced line for line

Three kinds of printed numbers, with different reproducibility:

1. **Corpus and structural counts** (graphs used, holes, split edges, candidate pairs, class-3
   counts, rank distributions, counterexample identity): the generators and `std::mt19937` are
   ported bit-exactly, so these **must** agree. They do (table below).
2. **Roundoff-level residuals** (`1e-14`-scale maxima): Eigen's SVD/QR is replaced by LAPACK, so
   `X0`, `Phi` and everything downstream differ at the last bits. Same magnitude, not the same
   digits; the acceptance rule is the tolerance the derivation states, not digit equality.
3. **Shape-space sample counts** (`check_r2`, `check_r3`, `check_t5_adm`, `check_l1`): samples are
   `X0 + Phi t` with Gaussian `t`. `Phi` is a null-space basis, unique only up to an orthogonal
   transform inside the null space, and LAPACK's differs from Eigen's. The sample *distribution* is
   the same, the sample *points* are not, so "N samples positively oriented" and every count over
   samples moves by a few percent. Verdict-level results (0 mismatches, 0 violations, 0 failures)
   are unchanged.

`check_b4_93` is the one program that does NOT reproduce; see its row.

## Expected numbers (C++ era) vs the Julia rerun (2026-09-19)

Quoted from the derivation files (`core.md` "Programs run" table l.1984-1991 and 2261-2263,
`check.md` R2/R3/R6, `lemmas.md` Check L1/L2 outputs, `round2_theorist_b.md` §0.g,
`persona_rigidity.md`). "=" means identical; "~" means same magnitude / same verdict.

| program | quantity | C++ (quoted) | Julia rerun | |
|---|---|---|---|---|
| `check_t1_t2` | graphs solved / deployable | 30 / 30 | 30 / 30 | = |
| | C1 closure | 1.60e−14 | 3.70e−14 | ~ |
| | C2 closed form vs `deploy()` | 1.42e−14 | 1.49e−14 | ~ |
| | C3 hinge angle | 9.73e−14 | 2.47e−13 | ~ |
| | C4 split duplicates | 2.08e−14 | 4.41e−14 | ~ |
| | C5 area constant | 4.46e−14 | 1.16e−13 | ~ |
| | C6 pin equation | 1.83e−14 | 3.73e−14 | ~ |
| | C7 conic | 1.68e−10 | 2.60e−10 | ~ |
| `check_t4_t5` | graphs used (split-free) | 16 (8) | 16 (8) | = |
| | D1 exact vs bisection (rule ≤ 1e−5) | 2.09e−6 | 1.18e−11 | ~ (both pass; the Julia bisection lands closer) |
| | D2 split-free vs min β | 8.88e−16 | 8.88e−16 | = |
| | D3 pairs / worst rel. dev. | 72 / 4.37e−15 | 72 / 9.22e−15 | = / ~ |
| | D4 violations / copies, worst ratio | 719 / 1061, 1.4141 | 719 / 1061, 1.4141 | = |
| `check_r2` | graphs used | 16 | 16 | = |
| | R2-A copies, worst | 1 628, 3.11e−15 | 1 628, 2.89e−15 | = / ~ |
| | R2-B drift/r_f at diameters 5.66…19.80 | 4.62, 7.71, 10.79, 13.87, 16.95, 20.03 | 4.62, 7.71, 10.79, 13.87, 16.95, 20.03 (2-dec.) | = |
| | R2-B ρ_spec/r | 1.000000 ×6 | 1.000000 ×6 | = |
| | R2-C split-free patterns, worst | 8, 8.88e−16 | 8, 8.88e−16 | = |
| | R2-D samples / violations | 197 / 0 | 193 / 0 | ~ (sample kind 3) / = |
| | R2-E tests (withdrawn line) | 2 363 380 | 2 227 276 | ~ (kind 3) |
| `check_r3` (eps 0.02) | graphs / samples | 16 / 197 | 16 / 193 | = / ~ |
| | harmonics decided / ambiguous | 2 145 387 / 47 | 2 009 262 / 68 | ~ (kind 3) |
| | class sizes 1 / 2 / 3 | 2 070 850 / 73 446 / **1 138** | 1 937 338 / 70 854 / **1 138** | ~ / ~ / = (class 3 is structural) |
| | round-2 list mismatches (all class 3) | 153 | 176 | ~ (kind 3; still all class 3) |
| | round-3 list mismatches | **0** | **0** | = |
| `check_b1` | graphs / holes / split edges | 10 / 471 / 58 | 10 / 471 / 58 | = |
| | B.1 closed form vs fit (rel) | 8.7e−14 | 9.7e−14 | ~ |
| | u closure, fit resid, B.2s | (not quoted) | 2.7e−15, 1.5e−14, 6.4e−15 | |
| | B10 argmax − θ_c/2 | (not quoted) | 6.40e−4 (= the 2000-point grid step) | |
| `check_t5_adm` | corpus | 7 cases | 7 cases | = |
| | A: pairs tested | 1 363 | 1 365 | ~ (kind 3) |
| | A: ratio range | [0.5404, 1.0] | [0.5439, 1.0] | ~ |
| | A: max \|s/‖e‖² − L′/L\| | 1.55e−15 | 1.67e−15 | ~ |
| | A: admissibility failures | **0** | **0** | = |
| | B: counterexample | hexagons, θ₁ = π/3, Θ_max = 2π/3, ε = 1.570796 | same | = |
| `check_l1 140 20` | graphs used | 84 (7 ref + 77 random of 140 ids) | 84 | = (corpus `derivation_l1_140.json`, ids 0..139) |
| | candidate harmonics / coincident-copy pairs | 30 119 368 / 801 816 | 31 343 774 / 824 396 | ~ (kind 3: more samples survive the orientation screen) |
| | shared but C ≠ 0 (must be 0) | 0 | 0 | = |
| | class disagreements | 0 / 801 816 | 0 / 824 396 | = |
| | class-3 by shared-edge type | hinge 160 / split 24 / none 100 | 160 / 24 / 100 | = |
| | class-3 hinge ⟂ | 160 of 160 | 160 of 160 | = |
| | class-3 persistent / transient / undecidable | 58 / 152 / 1 259 | 59 / 152 / 1 259 | ~ / = / = |
| | closed-form errors p, q, r | 1.99e−12, 1.99e−12, 3.41e−13 | 2.93e−12, 2.93e−12, 3.41e−13 | ~ |
| | hinge dS pairs, err; split dS pairs, err | 479 784, 1.00e−13; 159 816, 8.64e−14 | 494 840, 2.32e−13; 166 152, 1.22e−13 | ~ |
| | persistent non-shared C = 0 by #frozen vertices 0:1:2:3 | 18 762 : 20 182 : 0 : 135 410 | 19 524 : 20 822 : 0 : 135 026 | ~ (kind 3). The earlier Julia line 0 : 0 : 9 986 : 170 738 was a port bug, fixed 2026-09-20: the frozen-vertex scan was written as one `for v, j` loop, whose `break` on the first nonzero Phi entry left BOTH loops, so every later vertex stayed marked frozen; the C++ breaks only the inner loop. The remaining spread is the shape-space sample-point effect (this row counts per-sample pairs, like the coincident-copy row) |
| `check_l1 52 20` (default) | | (not quoted) | 34 graphs, 10 482 152 pairs, 0 / 0 failures | |
| `check_l2` | patterns evaluated | 133 | 133 | = |
| | rank(A) = 2 rank(D) | 133 / 133 | 133 / 133 | = |
| | rank(D) distribution | 0→8, 1→1, 2→124 | 0→8, 1→1, 2→124 | = |
| | max \|M_j − formula(D)\| | 5.271e−14 | 5.478e−14 | ~ |
| | max vanishing row | 1.592e−14 | 2.545e−14 | ~ |
| | squares_3x3 line | max\|D\| 3.08e−18, tr 0, det 1 | 5.55e−17, 0, 1, K0(1,0) = +1 | ~ |
| | covector test | 94 ok, 0 mismatch | 94 ok, 0 mismatch | = |
| | H2 excluded (R2 line, quoted only in `check_lemmas.md`) | snub_square_3x3 (consistency 27.50) | snub_square_3x3 (27.50) + voronoi_torus_1056_n12 (27.15) + voronoi_torus_1096_n36 (23.00) | = / two extra tori; the C++ `check_l2` printout that included this line was never quoted, only the Checker's own test on a different corpus |
| `rigidity_rig_check` | graphs | 21 | 22 rows (6 tilings + 8 + 8 random), 0 BAD, 0 skipped | ~ (the persona says 21; the program runs 22) |
| | V1 max error | 6e−15 at scale ~6 | 9.1e−15 at scale 8.4 (worst row) | ~ |
| | V2 m_A = m_R | 21/21 OK | 22/22 OK (ini and sol) | = |
| | V4 mobility drop at θ = 0 | squares 5×5 10→1, kagome 18→7, periodic 4→1 | 10→1, 18→7, 4→1 | = |
| `check_b4_93` | id 93 replay | feasible, scan 0.248400, referee 0, \|C\| = 3234, pruned 3 625 of 201 295 | **repair infeasible** (min_q 0.2327) | ✗ see below |
| | id 96 replay | Θ = 0.241884 = referee | feasible, Θ = 0.200776 = referee (0.163520 with a recomputed shape) | ✗ |

### `check_b4_93`: not reproduced

The C++ read `(X0, Phi)` for the `free` system from the B4 cache
(`results/kill/b4/cache/free_sigma_def/shape_<id>.bin`), which was not migrated. The port
recomputes the shape by default, and can read the C++ cache when `KIRI_B4_CACHE=<dir>` is set
(the source repo still has it). **Even with the C++ `(X0, Phi)` loaded, `zero_plus_repair` does not
return the archived design**: id 93 ends infeasible (b4.csv: feasible, first-pass `min_q = 0.202457`,
150 iterations) and id 96 lands at `Θ_max = 0.200776` (b4.csv: `0.241884`). The corpus, sigma_def,
seeds (`6000 + 7·id + 1`) and options (`n_random 3, max_iter 1200, lambda_rel 1e-6`, then
`w_corner 0.05, w_prox 1e3`, warm start) are the C++ ones. So the divergence is inside
`Kirigami/src/method/zero_plus.jl` (`zero_plus_repair`, its L-BFGS, or `zero_plus_form`), not in
this driver -- to be investigated by the method port, not silently accepted. Until then the R6.5
argument of `check.md` (the referee's misfire at the shared hinge vertex of faces (94, 184)) cannot be
re-run on the Julia side; note also that the referee artefact itself was fixed in `collision.cpp`
afterwards (`results/core_validation/referee_fix.md`), so a faithful replay would now print referee
`0.2484`, not `0`.

Reproduce: `KIRI_B4_CACHE=~/Documents/kirigami-experiments/results/kill/b4/cache/free_sigma_def julia --project=Kirigami derivations/scratch/check_b4_93.jl 93 96`.

**Resolved 2026-09-20 (port-method-2): not a port bug.** With the C++ `(X0, Phi)` loaded the start
state is bit-identical and the first ~20 L-BFGS iterates agree to 1e-14; the 1199-iteration
secondary repair never converges and the rounding path diverges after that (1e-9 at 40 iterations,
4e-4 at 80), so the warm-started primary lands in a different local minimum and id 93's feasibility
flips on one corner margin out of 5206. Probes, the full iterate table and the verdict:
`derivations/scratch/b4_path/README.md`. The b4.csv rows for ids 93 and 96 (the 0.2484 / 0.241884 rad
values) are knife-edge outcomes of a non-converged optimiser, not reproducible numbers.
