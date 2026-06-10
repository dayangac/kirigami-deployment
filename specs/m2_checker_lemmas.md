# Mission 2 / WP1 — Checker-L: independent verification of derivations/lemmas.md

Paste of specs/common_preamble.md applies (read it first). No web access needed.

## Your position
You are the adversarial Checker for two new lemmas. **On your first pass you must not read derivations/lemmas.md.** You read only derivations/lemmas_statements.md (the bare statements, hypotheses and the definitions they use), the notation in derivations/core.md §0 and T1, and the code. You try to prove each statement yourself, or break it. Only after you have written your own verdict per statement do you open lemmas.md, compare the Deriver's proof step by step, and list every disagreement. This ordering is the point; record in your report that you followed it.

## Read first
1. derivations/lemmas_statements.md — the statements (L1.1, L1.2, L1.3, L1 corollary; L2.1, L2.2, L2 corollary) with their hypotheses.
2. derivations/core.md §0 (notation, sign conventions), T1 (closed form `Y_θ = cos(θ/2)C + sin(θ/2)S`, face potential `u`), T3 (harmonic predicates; T3.H.1–H.5), T5.1 and the deflation block around lines 1405–1520 (the τ-chart, classes 1/2/3, Lemma T5.1e), T5.3 (split-edge pairs have `p + q ≡ 0`).
3. derivations/check.md entries D5, R2.5, R3.1a, Round 6 — your predecessors' history on this material.
4. code/src/method/deploy_basis.hpp, contact.hpp (`harmonic_roots`, the class split, the identity test), zero_plus.hpp (header derivation of `dS_e`, `dC_e = 0`), periodic_jacobian.hpp (header derivation of `J(θ)`, `K = Q P₀⁻¹`, `AchievableSet` with `A`, `dimK`, `D`, `rankD`).
5. code/tests/derivation_tests.cpp — the existing standalone test file (33 cases / ~89k assertions); its build line is documented in derivations/check.md. Your new tests go here.
6. STATE.md F1, F11, F14, F22, F33.

## Task
**Pass 1 — blind.** For each statement: attempt your own proof sketch (numbered steps) or construct a counterexample. For L1.1 the natural attack is a pair type the list omits or includes wrongly; try all predicate types (vertex-vs-edge orientation, edge–edge crossing, squared distance) × all incidence relations (same M vertex two copies; vertex on an incident edge; split-edge duplicate pair; hinge-adjacent faces sharing the pin; faces sharing only a vertex; disjoint faces). For L1.2, attack the sign claim at θ = π and the "never a contact" reading (is a constant-sign harmonic really never a crossing? what about a permanent touching?). For L1.3, look for a pair with `p = q = r ≡ 0` that is not a listed permanent incidence. For L2.1, attack the factorisation: does the x-copy and the y-copy of a null vector really give two *independent* matrices whenever the increment row is non-zero? What if `P₀` has a special form (square lattice) — can the two collapse? For L2.2, attack (a) `rank(D) ≤ 2` trivially true by shape, so check what the statement actually asserts; (b) the iff; (c) the sufficient condition, if one is stated.
Write derivations/check_lemmas.md "Pass 1" with a verdict per statement: PROVED-INDEPENDENTLY / PLAUSIBLE-UNPROVED / COUNTEREXAMPLE (with the example) / STATEMENT-ILL-POSED (with why).

**Pass 2 — compare.** Open derivations/lemmas.md. For each proof step, AGREE / DISAGREE with the reason; list every disagreement as D-L1.., D-L2.. with the exact line. Distinguish: wrong (a false step), gap (a true step without justification), hypothesis (a needed hypothesis not stated), wording. Do not concede a point you have not verified.

**Tests.** Add to code/tests/derivation_tests.cpp (keep the file standalone as it is; follow its existing style and the build line in check.md):
- L1: on ≥ 1,000 random `(graph, σ, t)` triples across kill_common.hpp's population (≥ 50 graphs) and the reference tilings, for every candidate pair: predicted class from the combinatorial type per L1.1/L1.3 vs the numeric class from the coefficients; assert zero unexplained instances, where "explained" means the pair is on the accidental set and moving `t` by a random 1e-3 perturbation changes its numeric class. Also the L1.2 sign check on 10,000 random class-3 harmonics including θ = π.
- L2: `rank(A) == 2·rank(D)` on the 33 K7 patterns and ≥ 100 random Voronoi tori with random σ (find the generator in code/apps/kill_k7.cpp), dense SVD, tolerance 1e-10 relative; the reconstruction of each `M_j` from row `j` of `D` by the Deriver's formula vs the assembled `M_j`; `rank(D)` distribution reported. Counterexample search: 2,000 random (torus, σ) draws; assert none breaks the equality.
Run the whole file; report cases / assertions / failures.

**Statements file hygiene.** If lemmas_statements.md and lemmas.md differ in any statement, that is itself a disagreement; report it.

## Output
derivations/check_lemmas.md with: Pass 1 (blind verdicts, your own sketches), Pass 2 (step-by-step agreement table, disagreement list), Tests (what was added, the build line, the run output), and a one-paragraph overall verdict per lemma: VERIFIED / VERIFIED-UNDER-HYPOTHESIS (name it) / DISPUTED (list) / REFUTED (example).

## Rules
- Blind first pass; record that you did it.
- Do not edit derivations/lemmas.md or core.md; the Deriver answers your disagreements in a second round.
- Do not commit.

## Reply
≤15 lines: verdict per statement after each pass; the disagreement ids with one line each; test counts; anything you could not check.
