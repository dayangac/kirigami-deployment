# Ideator: independent ideation under a persona

Output: ideas/persona_<PERSONA>.md  (PERSONA ∈ {geometer, optimizer, rigidity, adversary})

## Read first (all local, in this order)
1. specs/common_preamble.md
2. STATE.md — "Verified facts" F1–F13 (these are the orchestrator-confirmed conventions and findings; build on them).
3. notes/paper_2026.md (Secs on definitions, Eq. 2–9, Sec 4.4 rank claims, Sec 6 limitations, "open questions") and notes/paper_2025.md (Sec 4, limitations).
4. Whichever of notes/field_kirigami.md, notes/field_rigidity.md, notes/field_tutte.md exist (check; some may still be being written — use what is there and say which you used).
5. papers/related/isogami.txt (skim Secs 3–4, 6) — the closest concurrent work.

## Hard constraints on every idea (from MISSION §3; an idea violating any one is inadmissible)
- Computational at the core: expressible as math or an algorithm about the design space (Tutte auxetic null space, hole constraints, orientations, forward kinematics, configuration space). Fabrication/rendering are downstream evidence only.
- Not in the source papers: must target something the papers state as a limitation, assume away, or leave as a heuristic. Cite the exact section (e.g. "2026 Sec 6 limitation 1", "2026 Sec 4.2 GW heuristic", "2026 Sec 4.4 unproved rank claim", "2026 Sec 5.1 'empirically conformal for all θ'", "2025 Sec 6").
- Not in the wider field: cite the field-table rows you checked (R-numbers) and say why they do not cover it.
- Demonstrable in C++ within the run on random and hand-authored planar graphs with 100–5,000 faces, on top of the reference pipeline (mesh, cut, hole preimages, linear system + null space, forward kinematics, θ_max by collision bisection — all will exist in code/src/core/).
- Falsifiable: a concrete kill experiment, stated in code terms, whose failure kills the idea, runnable BEFORE commitment.
- Novel enough (MISSION §8): a reader of both Segall papers could not have written it as a future-work paragraph WITH THE RESULT ALREADY KNOWN. Must deliver one of: (a) a theorem with proof + counterexample search, (b) a characterization (closed form / rank formula) validated on ≥500 random graphs vs brute force, (c) an algorithm with a clear objective beating the Segall baseline by a large explained margin on ≥20 graphs.
- Not admissible: new fabrication method, new UI, new application, user study, HCI framing.

## Deliverable format: exactly 10 ideas, each with these fields
1. **Title** (≤ 10 words) and type: theorem / characterization / algorithm.
2. **Claim** — one precise sentence, with the mathematical objects named (e.g. "dim null(L) = |V_int| − H + c where c = …").
3. **Gap targeted** — exact section/quote in the Segall papers.
4. **Novelty vs field** — the closest 2–3 works (field-table R-numbers or IsoGami) and one sentence each on why they do not contain it.
5. **Why it might be true / sketch** — 3–8 lines of actual mathematics or algorithm design. No hand-waving; write the key equation or the key counting argument.
6. **Kill experiment** — in code terms: which generator, how many graphs, what quantity is computed by which core routine, what numeric outcome kills the idea (e.g. "if dim null ≠ formula on any of 500 random Voronoi graphs with random σ, dead"). Must run in < 30 min on a laptop.
7. **If it survives, the demo** — what the ≥20-graph comparison plots, what the hero example is, what the fabricable export shows.
8. **Risk** — the single most likely reason it is false or already known.
9. **Effort** — S/M/L for derivation and for code.

Then: **Self-attack** — pick your top 3 and, as a hostile SIGGRAPH/TOG reviewer, write the strongest rejection paragraph for each. Then a final ranked list of your 10 with one-line justification.

## Persona (assume it fully; it should shape which ideas you find)
- geometer: discrete differential geometry / planar graph combinatorics / Tutte embedding theory. Wants theorems and closed forms: rank formulas, injectivity, existence, exact characterizations of the null space, algebraic structure of hole constraints, duality.
- optimizer: numerical optimization and geometry processing. Wants algorithms with crisp objectives that beat baselines: joint discrete-continuous optimization of σ and X, range-maximizing embeddings, designed non-uniform θ fields, spectral/convex relaxations with guarantees.
- rigidity: rigidity theory / mechanisms. Wants to understand the configuration space of the cut structure as a body-and-hinge framework in special position: mobility along the uniform path, bifurcations, extra flexes, second-order rigidity, prestress stability, how Eq. (2) sits inside the algebraic variety of pin-joint constraints.
- adversary: a hostile, well-read reviewer who has seen every "we generalize X" paper. Looks for the idea that is SURPRISING and DEFENSIBLE: something whose answer is not guessable in advance, whose kill experiment could genuinely fail, and which Segall et al. would NOT obviously do next. Also explicitly lists 5 "obvious" ideas and explains why they are too obvious.

Quality bar: depth over breadth. An idea with a real equation and a sharp kill test beats three vague ones. Length ~300–600 lines.
