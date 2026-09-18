# K8a recheck — were the 84/100 "near-Farkas certificates of structural emptiness" real?

**Short answer: the solver was wrong, but not in the direction that kills the certificates.**
The *primal* half of `cone_lp` was broken and reported `margin = 0` on systems whose
feasibility is provable by hand; the *dual* half was sound, and every certificate it
produced survives an independent recomputation. So the emptiness reading of K8a's
`0/100` survives — but K9's `0/30` was pure solver failure, and with the corrected solver
the same 30 embeddings give **29 strictly positive non-uniform margins**.

---

## 1. The bug

`cone_lp` solves the matrix game `val = max_{||z||<=1} min_i a_i.z = min_{lambda in
simplex} ||A^T lambda||_2`. The old implementation computed the two ends of the bracket
from **two unrelated algorithms**:

* the **dual** end by Frank–Wolfe on `min_lambda ||A^T lambda||` — correct, every iterate
  is in the simplex so every value is a rigorous upper bound;
* the **primal** end by a separate log-sum-exp smoothing loop with projected gradient
  ascent on the unit ball, annealed `mu = 1 -> 1e-5` in 17 stages of 60 iterations
  (`expansive_cone.jl`, lines 208–243 of the pre-fix file), started from `z = 0` with step
  `mu / smax^2`.

That second loop is hopelessly under-converged on the real systems (about 1 000 steps of
size `1e-5` in a space of dimension up to 300, against 5 000 rows) and it was the **only**
source of the returned witness. Nothing in the code ever evaluated the Frank–Wolfe iterate
`w = A^T lambda` as a primal direction — and that iterate *is* the answer:

> At the minimum-norm point `w*` of `conv{a_i}` the supporting-hyperplane inequality
> `a_i . w* >= ||w*||^2` holds for every `i` (otherwise the segment from `w*` towards
> `a_i` would contain a point of smaller norm). Hence `z* = w*/||w*||` is a unit vector
> with `min_i a_i . z* >= ||w*||`, so `margin >= dual bound`; combined with
> `margin <= val <= dual bound` this forces `margin = val = ||w*||`. **The primal witness
> is not a separate object — it is the normalised dual iterate.**

The consequence is exactly the K9 signature: the dual converged to `2e-2 … 5e-2` (a true
upper bound, carrying no infeasibility claim, since the threshold for that is `1e-7`) while
the primal sat at `0`, leaving an *open* bracket `[0, 3e-2]` that the driver's
`feasible = margin > 1e-8` test read as "no witness found".

## 2. The fix

`Kirigami/src/method/expansive_cone.jl`, `cone_lp`:

* the dual loop is now **away-step Frank–Wolfe** (Lacoste-Julien & Jaggi, NeurIPS 2015),
  which is linearly convergent on a polytope where vanilla FW is `O(1/t)` and zig-zags at
  the boundary; the active set is maintained through drop steps, and `lambda` is
  renormalised onto the simplex after each drop;
* **the primal witness is read off the same iterate** (`try_primal(w)` plus the per-step
  `cg(is)/||w||`, which is exactly the margin of `w/||w||`);
* `w` is recomputed as `A^T lambda` every 256 steps against drift, and once more at exit;
* the smoothing loop is kept only as an extra source of witnesses — the reported margin is
  the best of the two, so nothing that used to be found is lost;
* new `ConeLPResult.gap = dual_bound - margin_l2`, the width of the rigorous bracket.

Default `dual_iters` raised 600 -> 20 000.

Two new report fields make the solver testable without trusting it:

* `ExpansiveConeReport.sigma_chart_margin` — the uniform ray `sigma`, projected into the
  flex basis and evaluated on the branch chart it itself selects, divided by its norm.
  **No solver touches this number.** It is a rigorous lower bound on that chart's LP
  value, so whenever `sigma_in_cone` holds it must be positive and the solver must return
  at least it.
* `farkas_residual(A, lambda)` — recomputes `||A^T lambda||_2` row by row from the
  multipliers alone, and returns `min_i lambda_i` and `|sum lambda_i - 1|`. Applied to
  every pass inside `expansive_cone`; the report carries the worst residual, the worst
  multiplier, the worst simplex error, and the number of passes certified at `1e-9` and
  at `1e-6`.

## 3. Step 1 — the known-feasible test (`results/kill/k8a/recheck_k9.csv`)

`Kirigami/apps/kill_k8a_recheck.jl` regenerates K9's variant-(b) embeddings bit-identically
(same seeds, same three starts, same two fallbacks, `delta = delta' = 1e-3 med^2`,
600 L-BFGS iterations, 6 barrier stages) and runs the corrected LP at each. All **30**
designs K9 ran the cone on are reproduced.

| | old solver | corrected solver |
|---|---|---|
| designs with `sigma` in `P(X)` by direct `zero_plus` measurement | 23 | 23 |
| ... of which LP margin `> 0` | **0 / 23** | **23 / 23** |
| ... margin `>=` the solver-free `sigma_chart_margin` | — | **23 / 23** |
| ... a non-uniform flex strictly beats the uniform ray | — | **23 / 23** |
| designs with `sigma` **outside** `P(X)` | 7 | 7 |
| ... of which LP margin `> 0` (a genuine non-uniform rescue) | **0 / 7** | **6 / 7** |
| feasible overall | **0 / 30** | **29 / 30** |

`margin / sigma_chart_margin`: median **8.73**, range 2.23 – 56.9. So the answer is not
merely "the solver now recovers the uniform ray" — at every one of these embeddings the
best non-uniform flex separates the cuts and corners by roughly an order of magnitude more
than `sigma` does.

Bracket width `dual - margin` over the 30: median `6.3e-4`, max `4.5e-3` at 20 000 steps
(the LP value itself is `9e-4` … `4.4e-2` on the 23 in-cone designs, so the bracket is 1–50 % wide; it is a genuine
positive either way, and the ordering `margin >= sigma_chart` is what the test turns on).

The one design that stays infeasible is `id 46 / sigma_def`, where `sigma` is outside
`P(X)` (`min mu = -1.26`) and the dual reaches `3.6e-5`.

## 4. Step 2 — independent verification of the certificates

Every dual certificate is now recomputed by `farkas_residual`, which re-forms
`sum_i lambda_i a_i` row by row from the multipliers alone, with no reuse of any solver
state, and checks `lambda >= 0` and `sum lambda_i = 1`. A graph counts as certified only
if **every branch chart the driver solved** passes.

Over all 277 configurations of the re-run: worst `min_i lambda_i` is exactly `0`
(multipliers are non-negative; drop steps zero them out), worst `|sum lambda_i - 1|` is
`2.5e-13`. **No certificate was ever malformed.** The multipliers are genuine simplex
points and the residuals reproduce.

K1a population, 100 000 dual steps (`results/kill/k8a/recheck100k/`):

| | at `X_ini` | at `X0` |
|---|---|---|
| reported `dual_max < 1e-6` (original run) | 84 / 100 | 68 / 100 |
| **verified** `\|A^T lambda\| < 1e-6`, every pass | **96 / 100** | **99 / 100** |
| **verified** `< 1e-9`, every pass | **65 / 100** | **57 / 100** |
| verified `< 1e-5`, every pass | 98 / 100 | 100 / 100 |
| worst verified residual | `3.9e-5` | `1.8e-6` |

Per family at `X_ini`, verified `< 1e-9` (and `< 1e-6` in brackets): Voronoi **34 / 34**
(34), quad-random **31 / 33** (33), Delaunay **0 / 33** (29). The Delaunay shortfall is the same convergence residue
the original report identified, not a different phenomenon: those are the graphs with
`dim_flex` up to 306, and their residuals sit at `2.5e-9 … 3.9e-5`. Away-step Frank–Wolfe
improves them by two to three orders over the original vanilla Frank–Wolfe at the same
budget (Delaunay `dual_max` median `9.3e-7 -> 5.7e-9`), but does not reach `1e-9`.

**What a verified residual `r` proves, exactly.** `lambda` on the simplex gives
`val <= ||A^T lambda|| = r`, and the corrected primal returns `margin = 0`, so the LP value
of that chart lies in `[0, r]`. That is a bound, not exact emptiness — it is *not* a
rational Farkas certificate and cannot be promoted to one without exact arithmetic. It is
however now a **closed** bracket: `gap = dual - margin` equals `dual_max` on every K1a row
because the margin is identically 0, so the LP value is pinned to within `4e-10` (median).

## 5. Step 3 — the K1a population, re-run

`Kirigami/apps/kill_k8a.jl`, 100 graphs, both embeddings, 8 shards, corrected solver:

| | at `X_ini` | at `X0` |
|---|---|---|
| LP margin `> 0` (PASS bar was `>= 20 / 100`) | **0 / 100** | **0 / 100** |
| `sigma` in `P(X)` | 0 / 100 | 0 / 100 |
| verified emptiness certificates, every pass, `< 1e-6` | 96 / 100 | 99 / 100 |

**K8a's verdict is unchanged: FAIL on the PASS rule, 0 of 100 at both embeddings.** The
corrected solver does not find a single witness the old one missed on this population, and
the certificates it does produce are stronger and now independently verified.

The controls also re-run clean and strictly improve, which is the second signature of the
old primal's under-convergence — the margins on the four authored tilings all rise, and
the bracket now closes to `1e-9` or better:

| tiling | old margin | new margin | bracket width |
|---|---|---|---|
| `hexagons_auto` | 5.40e−1 | **6.21e−1** | 6.1e−13 |
| `truncated_square_488` | 7.35e−1 | **8.52e−1** | 9.2e−11 |
| `snub_square_33434` | 3.53e−2 | **6.58e−2** | 1.0e−9 |
| `tiling_3_4_3_12` | 3.12e−1 | **3.12e−1** | 1.0e−10 |

Hard soundness control: `sigma` in `P(X)` on 4 / 4, LP margin `> 0` on 4 / 4 — unchanged.
Deployable population (73 configs): `sigma` in `P(X)` 71 / 73, margin `> 0` **72 / 73**,
Euler step collision-free **72 / 72** — all unchanged.

## 6. Step 4 — non-uniform flexes at K9's constrained embeddings

Covered by section 3: at all 23 constrained embeddings where the uniform ray is admissible,
a non-uniform flex with a strictly larger margin exists, by a median factor of **8.7**
(range 2.2 – 56.9). At **6 of the 7** where the uniform ray is *not* admissible
(`min mu < 0`, so `Theta_max = 0` along `sigma`), a non-uniform flex nevertheless separates
every cut and every corner at first order.

This is the first population on which the expansive cone is non-empty at a *designed*
embedding of a random graph, and it says the constrained embedding buys more first-order
room than the uniform deployment uses. It is a **first-order** statement only: nothing here
integrates the flex, and K8a's own Euler-step sanity was not re-run on these 30 (the
`euler_step` predicate is available but the recheck driver does not call it).

## 7. What this does and does not change

**Changed.**

1. `cone_lp`'s primal was wrong. Any "no witness found" from the old solver is worthless on
   its own; only its dual bounds ever carried information. F36's parenthetical *"expansive-
   cone LP 0/30 feasible ... ⇒ cone_lp NON-CONVERGENCE"* diagnosed the symptom correctly and
   is now confirmed with a cause and a fix.
2. K9's `expansive-cone LP: run 30, feasible 0` line in `results/kill/k9/summary.txt`
   (and the `run 115 / 6 / 109 / 13 / 17, feasible 0` lines for variant (a)) are **solver
   artefacts and must not be quoted.** The corrected number for variant (b) is **29 / 30**.
   Variant (a)'s 115 + 6 + 109 + 13 + 17 configurations were not re-run here.
3. The K8a certificate counts rise and become independently verified: 84 -> **96** at
   `X_ini` and 68 -> **99** at `X0`, at `1e-6`, with 65 and 57 at `1e-9`.

**Unchanged.**

4. K8a's headline **FAIL (0 of 100, bar `>= 20`)** stands at both `X_ini` and `X0`.
5. Caveat 2 of the original report is untouched and remains the binding limitation:
   **the branch enumeration is not exhaustive.** Four charts of `2^{2338}` are visited, so a
   verified residual certifies *those charts*, never `P(X)`. "Structural emptiness" remains
   evidence, not a theorem. What the recheck adds is that the evidence is now real evidence
   rather than a solver failure — and, pointedly, that on the one population where `P(X)` is
   known to be non-empty the same four charts find it immediately.
6. F30 / F25 (the `0+` split-duplicate collision at the Eq. (6) projection) are untouched;
   nothing here concerns the uniform ray at `X0`.
7. Caveats 1, 3 and 4 of the original K8a section stand as written.

**F-facts.** F36's cone clause needs one correction and one addition; no other F-fact
changes. Proposed amendment: replace *"Expansive-cone LP 0/30 feasible on (b)-feasible
designs although sigma in P(X) on 23/30 ... K8a's dual 'certificates' on 84/100 must be
re-examined"* with *"Expansive-cone LP: `cone_lp`'s primal was broken (the Frank–Wolfe
iterate was never evaluated as a witness); with the away-step min-norm-point solver the LP
is feasible on **29/30** (b)-feasible designs, with a median margin **8.7x** the uniform
ray's and 6/7 rescues where sigma is outside P(X). K8a's dual certificates were sound and
are now independently verified at **96/100** (`X_ini`) and **99/100** (`X0`) at `1e-6`;
K8a's 0/100 FAIL stands."*

## Reproducing

```bash
julia --project=Kirigami -e 'using Pkg; Pkg.test()'   # 4 K8a-recheck cases included

# step 1 / step 4: K9's variant (b) embeddings, 12 shards
for i in $(seq 0 11); do julia --project=Kirigami Kirigami/apps/kill_k8a_recheck.jl --out results/kill/k8a \
    --shard $i --nshards 12 & done; wait

# steps 2-3: the K1a population with the corrected solver
for i in $(seq 0 7); do julia --project=Kirigami Kirigami/apps/kill_k8a.jl --n 100 --x0 --dual-iters 100000 \
    --out results/kill/k8a/recheck100k --cache results/kill/cache \
    --shard $i --nshard 8 & done; wait
```

Artifacts: `results/kill/k8a/recheck_k9.csv`, `recheck_k9_summary.txt`,
`results/kill/k8a/recheck/` (20 000 steps, includes the reference and deployable controls),
`results/kill/k8a/recheck100k/` (100 000 steps, K1a only). The original `k8a.csv` and
`summary.txt` are untouched.
