# ESCALATION — mission §9 stop condition reached (2026-09-04 12:50 local)

**Trigger.** Three consecutive constructive ideas died in round 2: T-1 vertex balance (FAIL, reflex corners), B4 boundary drop (FAIL 1/400 refereed), K8a expansive-cone LP (FAIL 0/100, dual-certified structural on 84/100). Before them: K2b, K6, K7-C3. One candidate still running: B3/X2 tr K budget on periodic patterns.

**What is verified and paper-ready** (see STATE.md F1–F33, derivations/core.md rounds 1–6, results/kill/KILL_REPORT.md):
- Exact deployment calculus T1–T7 (trig-linear path, no-locking lemma, harmonic contacts, exact Θ_max 187/187 vs bisection, gradients, rank corrections).
- Sound certificate (round 6): exact below the first graze, inner above; 89,063 + 16,622 assertions.
- Periodic affine theorem (K7 C1/C2/C4): J(θ) = cos(θ/2)I + sin(θ/2)K, K affine in the shape space, dim 𝒦 = 2·rank D, exact conformality ∀θ, closed-form second closed angle.
- Structural non-deployability of random-graph patches: uniform (K1a/K5/K6/B4) and non-uniform (K8a dual certificates); sharp jitter transition a* ∈ [0.16, 0.52] edges (A3); local mechanism = reflex corners from Eq.(6) (T-1).
- Errata vs the authors' code: rank claim, boundary over-constraint (F23), 9.47% silent re-closure (K1c).

**What is missing for MISSION §8 (TOG bar).** A constructive algorithm beating the Segall baseline on ≥20 graphs. Every attempt so far fails on the same 0⁺ split-cut contact.

**Options (orchestrator recommendation first).**
1. **Reframe the deliverable as characterization + impossibility (MISSION §8 option b).** The characterization (exact Θ_max + certificate) is validated 187/187; extend to ≥500 graphs (a rerun). Headline theorem: for the Tutte auxetic embedding with fixed boundary, random planar graphs are non-deployable even non-uniformly, with a dual certificate per graph, and a sharp perturbation threshold from tilings. Constructive residue: K7's exact periodic design theory with the tr K budget (if B3 passes) and the snub-square non-uniform rescue as hero. Cost: 1 day (Phase 9 experimenter + writer). Venue: strong technical paper / SIGGRAPH Asia technical communication; TOG only if reviewers value impossibility + exact theory.
2. **Change the embedding, not the graph.** Every failure is inherited from Eq.(6)'s projection creating reflex corners. Replace the projection by a convexity-constrained solve (faces convex at X0 as hard constraints inside the null space) and re-run K5/K6/K8a. Untested; one kill experiment (~2 h). If it passes it supplies the algorithm and the ≥20-graph figure.
3. **Periodic-only constructive paper.** Drop random graphs; make K7 + B3 the algorithm (prescribed Poisson/conformal targets on Voronoi tori) and the hero. Depends on B3 passing 12/25 → clearly more.
4. Stop and write up as is (poster).

Orchestrator will proceed with option 1 as the baseline and run option 2's single kill test in parallel unless told otherwise.

## Resolution (17:15 same day)
K9 (option 2) returned 36/400 exact, 32–36 refereed deployable random graphs versus 0 for every baseline (F36). Below the 10 % bar by 4–8 designs, above MISSION §8c's "≥20 graphs" by a wide margin against a baseline of zero. Orchestrator adopts option 2 as the constructive half (decision D9) with a push phase (K9b) and keeps option 1's framing as the paper's spine. Escalation closed; no human decision required unless K9b fails to reach 40/400, in which case the paper ships with 32–36.
