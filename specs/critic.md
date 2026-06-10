# Critic: rank all ideas, attack the top 5

Output: ideas/ranking.md

Read first: specs/common_preamble.md; STATE.md (Verified facts F1–F14 and Unverified claims U1–U9); MISSION constraints below; all four ideas/persona_*.md files; notes/paper_2026.md Sec 6 (limitations) and its errata/open-questions lists; the three notes/field_*.md tables; papers/related/isogami.txt Secs 3–4, 6.

## Constraints every idea must meet (MISSION §3, §8) — score each idea 0/1 on each
C1 computational at the core (math/algorithm about design space, hole constraints, null space, FK, configuration space).
C2 targets a stated limitation / assumption / heuristic of the Segall papers, with exact section cited.
C3 not in the wider field (field tables + IsoGami); cite rows.
C4 demonstrable in C++ within the run on random + authored planar graphs, 100–5,000 faces, on top of code/src/core.
C5 falsifiable: a concrete kill experiment in code terms, runnable before commitment, < 30 min.
C6 novel enough: a reader of both papers could not have written it as a future-work paragraph WITH THE RESULT KNOWN; delivers (a) theorem + proof + counterexample search, or (b) characterization validated on ≥ 500 random graphs vs brute force, or (c) algorithm beating the Segall baseline by a large explained margin on ≥ 20 graphs.
C7 not a fabrication method / UI / application / user study.

## Tasks
1. Build a master table of ALL ideas across the four personas (dedupe: merge ideas that are the same claim under different names, list which personas proposed each — convergence is evidence, not a tie-break). Columns: id, title, type, personas, C1..C7, kill experiment (one line), effort, one-line verdict.
2. Identify BUNDLES: sets of ideas that together form one coherent TOG-style contribution (one theorem family + one algorithm + one demo). The likely candidate is the "closed-form theory of uniform deployment" bundle (trig-linear path; closed-form collision times/θ_max; exact conformality; hole count via hinge graph Γ; mobility operator). Evaluate whether the bundle is a paper or a collection of paragraphs; say what the load-bearing theorem is and what the headline algorithm/figure would be.
3. Rank the top 5 ideas or bundles. For each of the top 5 write a hostile TOG reviewer's rejection paragraph (strongest possible: prior art, triviality, "folklore for rotating squares", boundary artefacts, does-not-generalize, demo too easy), then a fair rebuttal, then a verdict: PURSUE / PURSUE-IF-KILL-PASSES / DROP.
4. For each of the top 5, restate the kill experiment as an exact spec that an Experimenter can run against code/src/core with no further design decisions: generator(s), number of graphs, σ choice, routine to call, quantity to compute, tolerance, PASS/FAIL rule, expected wall time. Where two ideas share a kill experiment, say so.
5. Novelty holes: list every novelty claim that rests only on field tables and say which targeted literature search the Phase-6 Scout must run (exact query strings).
6. Final recommendation: which single idea or bundle the orchestrator should commit to first, and the fallback order.

Honesty: if an idea is likely folklore or trivial, say so even if it is elegant. Do not reward volume. Length: as needed (~400–700 lines).
