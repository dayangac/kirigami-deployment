# Hostile TOG review — K9/K9b constructive result (convexity + split-inward constrained embedding)

Reviewer stance: computational-design / physical-based-modeling reviewer (Pauly/Mahadevan-adjacent), skeptical of incremental barrier-method claims dressed as theory.

## 1. Is 36–37/400 (9%) a result or a curiosity?

Honest reading: **a curiosity that the authors correctly did not oversell, but is too thin to anchor a paper on its own.**

- ε_max distribution (K9b, n=28 certified; K9's 36 track the same shape): median 0.15–0.25 rad, only **13–16/28–36 designs ≥ 0.1 rad**, and only a handful (**snub-adjacent Delaunay outliers**, max 1.83–2.0 rad) open past 0.3 rad. Roughly a third of positives sit at ε_max < 0.03 rad — a barely-perceptible wobble, not a deployment.
- The gallery (`k9b_gallery.png`) makes this visible directly: rows 1–3 (0.91→0.12 rad) show real fanning-open; rows 5–7 (0.086→0.002 rad) are **visually indistinguishable from the closed state** — the "θ=Θmax/2" panel looks identical to the closed panel. That is the bottom third of the 24-design gallery, exactly as flagged.
- 9% against a 10% self-imposed bar, with the honest admission (K9b) that pushing the solver harder (8 restarts, 10 barrier stages, δ-sweep) made the raw count *worse* (28 vs 36) and only a best-of-two-runs union recovers 37. That is not a robust yield; it is a number sensitive to which heuristic run you report.

Verdict: report it as "existence, not prevalence" — a certified, non-trivial (0 in every baseline) foothold, but do not call 9% a design algorithm's yield rate in an abstract without the ε_max-weighted qualifier.

## 2. Is the method new, or a barrier method with two obvious constraints?

**Mechanically: yes, it is a standard L-BFGS/softplus/log-barrier feasibility-then-proximity solve over two quadratic-in-t constraint families** (`code/src/method/convex_embed.hpp`): corner convexity (`cross_i(t) ≥ δ`) and split-inward sign (`q_e(t) ≥ δ'`). Nothing about the optimizer is novel — continuation on target, Gaussian restarts, geometric barrier decay are textbook.

**The intellectual content, such as it is, is upstream of the solver, in three places, and only one is load-bearing:**

- **T-1's mechanism is the actual insight**: 1,966/1,968 balanced pure vertices with a non-positive 0⁺ margin sit at a *reflex* corner created by Eq.(6)'s unconstrained projection (0 non-convex faces at X_ini → 22,103/156,220 at X0). The disjunction-vs-conjunction observation (convex corner: μ = max(−g1,−g2), satisfied if either half-plane clears; reflex: μ = min(−g1,−g2), needs both) is the reason convexity is the right constraint to add, not an arbitrary regularizer. This is genuine, if modest, geometric content.
- **Constraint (b) subsumes (a)'s face-orientation/area constraint** (every corner strictly convex ⇒ face convex and CCW) — a clean but small remark, not a theorem.
- **The null-space parametrization** (`X(t) = X0 + Φt`) is the paper's own Eq.(5)/(6) apparatus, not new; using it to keep the search small is standard projected-optimization practice.

What is *missing* as insight: there is no argument for *why* convexity + split-inward together (and not some other pair) should be sufficient, no characterization of when the constrained set is non-empty, and K9's own data undercuts a clean story — convexity alone buys 0% (necessary, not sufficient, as stated), and K9b shows two solver configurations can both be exactly feasible on the same design and still differ in Θ_max (9 Delaunay cases), meaning **feasibility of the two constraints does not determine deployability**. That is the single most important (and self-damaging) finding in the whole set: it shows the proposed constraint pair is *not* the mechanism, only a necessary filter, and the real selector is unidentified. A reviewer will read K9b's own conclusion — "margin feasibility does not determine deployability" — as evidence the "method" doesn't explain its own results.

Bottom line: this is "add two constraints to a barrier solver, motivated by a genuine but narrow geometric observation (T-1), that does not fully explain its own output." Publishable as a negative/diagnostic finding bundled with the impossibility results; not on its own as an algorithmic contribution.

## 3. Is the comparison fair?

**Partially, and this is a real gap.** The paper text (line ~1959, 2057 of KILL_REPORT) states baselines "0/400" for Eq.(6) alone (K1a), both σ (K5), all four K6 repairs, B4, K8a, **and "the authors' native `prevent` 0/80 on the capped comparison."** So there *is* a same-population comparison against Eq.(6)+Eq.(9)-native, but only on **80 of the 400** designs (STATE.md: "native-prevent comparison capped at 8 graphs, 240s timeouts") — a 20% subsample, not the full 200-graph set, and capped for a *runtime* reason (timeouts), not a principled sampling choice. On the *separate* K2b population (30 designs, 8 graphs), native prevent actually reaches Θ_max = π on snub square and beats the paper's own null-space range optimizer 15/7/16 (worse/better/tie) — i.e. native prevent is not a weak baseline in general, it just wasn't run to completion on the full K9 population.

A rigorous reviewer will ask: **run native prevent to completion on all 400 (accept the wall-clock cost, or subsample honestly to a stated fraction with a stated random seed) before claiming 0% for it.** Right now "0/80, capped" reads as "we didn't finish the fair baseline," which invites the accusation that the comparison was abandoned once it became expensive, not because it was known to fail.

## 4. Would "first deployable random-graph kirigami under the Tutte auxetic framework" be accepted?

Only with heavy qualification. Mandatory qualifiers:

- **"under an added convexity + split-inward constrained embedding, not the paper's own Eq.(6)/Eq.(9) pipeline"** — this is not a claim about the published method, it's a claim about a new constrained variant of it.
- **"on 9–9.25% of a specific 200-graph, Delaunay/Voronoi/quad population, not all random graphs"** — Voronoi and quad_random are nearly absent (1 quad, 10 voronoi vs 25 delaunay); the result is Delaunay-biased.
- **"with median opening angle 0.15–0.25 rad (≈9–14°), not full deployment"** — "deployable" without this reads as full-range opening.
- **"certified sound but the solver is a non-convex heuristic with no completeness guarantee — feasibility ≠ deployability was directly observed (K9b)."**
- Drop "first deployable" unless every prior negative (K1a/K5/K6/B4/K8a/native-prevent-on-this-population) is cited in the same sentence, or a reviewer will read it as an overclaim relative to an already-published (if narrower) capability — e.g. native prevent reaching π on snub square is a stronger single-instance deployment than anything in the 9% bucket.

## 5. Figure/table set that would convince

1. **Full gallery, all 200 graphs × 2 σ, sorted by ε_max descending**, with an explicit baseline column per row (Eq.(6) alone / native-prevent-if-run / K9b) so a reader sees the 91% zero column directly beside the positives — not a curated 24-panel subset.
2. **Yield vs. face count (|F|) curve**, positives/total in bins of |F|, with 95% Wilson CI bands, split by graph family (Delaunay/Voronoi/quad) and by σ — to show (or disprove) whether yield is a real geometric trend or noise at n≈30 per family.
3. **One hero, fully fabricated**: pick the ε_max = 1.83–2.0 rad Delaunay outlier, laser-cut it, physically actuate it, and photograph both closed and open states next to the simulated panels — a single physical artifact does more to convince a TOG reviewer than any additional table.
4. A **feasible-but-non-deployable vs. feasible-and-deployable scatter** (the K9b 9-design finding) plotted as ‖X−X_ini‖ vs. Θ_max, colored by feasible/infeasible — this is the paper's most interesting negative result and deserves its own figure, not a paragraph.

## 6. Three cheapest experiments to harden it

1. **Finish native-prevent on all 400 designs** (relax the 240s cap or run overnight/parallelized) — closes the fairness gap in §3, cheap (compute-only, no new code).
2. **δ margin sweep reported per-family with a monotonicity check** — K9b already found δ=1e-2 wins on over half of positives; a systematic sweep (not best-of-4 cherry-pick) with a held-out validation split would turn "which δ wins" from post-hoc selection into a principled default.
3. **Characterize the 9 K9b "feasible-but-Θmax=0" designs** — diff their returned X against the 27 retained positives (same graphs where available, or matched face-count neighbors) to see if a simple secondary invariant (e.g. distance from a specific corner to its margin, or a spectral gap) separates them. This is the cheapest path to turning the negative finding into a real predictor, which is more publishable than the 9% number itself.

## 7. Rejection sentence and best rebuttal

**Rejection sentence:** "The paper proposes two barrier constraints inside an existing null-space projection and demonstrates a 9% yield with median opening under 15°, on an incomplete baseline comparison (native collision-avoidance not run to completion on the same population), falling short of its own pre-registered 10% bar and offering no characterization of why the constraint pair succeeds where it does — this reads as a negative result with a small constructive residue, not a deployable-design algorithm."

**Best rebuttal:** "We agree 9% is a floor, not a target, and report it as such against seven baselines that are all exactly zero on the identical 400-design population — no prior method (published or ours) produces a single certified positive on random graphs. The contribution is not the yield rate but the *mechanism*: T-1 identifies reflex corners, created by the paper's own unconstrained projection, as the obstruction, and K9b's negative finding (feasibility ≠ deployability on 9 designs) is itself new evidence about the shape of the design space, reported honestly rather than hidden. We commit to completing the native-prevent baseline on all 400 designs and will revise the headline claim to 'a constrained-embedding lower bound on the non-empty fraction of the random-graph design space' rather than any deployability claim."
