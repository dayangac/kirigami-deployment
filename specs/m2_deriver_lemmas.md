# Mission 2 / WP1 — Deriver-L: two lemmas that are currently measured, not proved

Paste of specs/common_preamble.md applies (read it first). No web access needed.

## Why this exists
REPORT.md "Next steps toward a paper" items 2 and 3, and review/theory_review.md, say a TOG referee will ask for proofs of two results that the project has only measured:
- **L1, the T3 degeneracy catalogue** — which contact harmonics have `g(0) = p + q = 0` *identically on the whole shape space*, and that the class `g(0) = g′(0) = 0` is never a contact.
- **L2, the periodic dimension formula** — `dim 𝒦 = 2·rank(D)` for the achievable set of period-Jacobian matrices `K`.
You write `derivations/lemmas.md` with full statements and proofs, plus scratch C++ checks. A separate Checker, in a fresh context, will see only your statements and the code and will try to break them.

## Read first (in this order)
1. derivations/core.md §0 (notation and sign conventions, lines ~116–190), T1 (lines ~190–395: the closed form `Y_θ = cos(θ/2)C(X) + sin(θ/2)S(X)`, T1.3 corollaries, the face potential `u`), T3 (lines ~502–612: harmonic predicates, T3.H.1–T3.H.5), T5.1 (lines ~958–1000 and the deflation block ~1405–1520: the τ-chart `g(τ) = C + Bτ + Aτ²` with `A = p − q, B = 2r, C = p + q`, the three classes, Lemma T5.1e, the side condition `A ≠ 0`), T5.3 (lines ~1617–1683: split-edge pair has `p + q ≡ 0`), and T8 (what is proved / refuted / conjectured).
2. derivations/check.md: the entries D5, R2.5, R3.1a, and Round 6 — the Checker's history on exactly this material. Do not repeat a claim it refuted.
3. code/src/method/deploy_basis.hpp (the C, S basis), code/src/method/contact.hpp (`harmonic_roots`, the three-class deflation, the identity test) and code/src/method/zero_plus.hpp (the header comment derives `dS_e`, `dC_e = 0` for split pairs — a proof ingredient).
4. results/kill/KILL_REPORT.md §K7 (lines ~1180–1465) and code/src/method/periodic_jacobian.hpp (header comment derives `J(θ) = cos(θ/2)I + sin(θ/2)K`, `K = Q P₀⁻¹`, `Q = 2 Jrot [w_h w_v]`, the `AchievableSet` struct with `A` (4 × 2k), `dimK = rank(A)`, `D` (k × 2 period-potential increments), `rankD`).
5. STATE.md F1, F11, F14, F22, F33 (the cut structure, preimage partition, `L = R·D`, the periodic rank, the K7 result).

## L1 — the degeneracy catalogue

**What is measured** (core.md T3.H.5, round-3 table): of 2,145,387 candidate harmonics over the K2a corpus, 2,070,850 are class 1 (`C ≠ 0`), 73,446 class 2 (`C = 0, B ≠ 0`), 1,138 class 3 (`C = B = 0, A ≠ 0`); 14,928 further candidates are identically zero and are struck by the identity test before classification. The text asserts that `C = 0` holds identically on the shape space for (i) every permanent incidence at `θ = 0` (T3.H.2: two M′ copies of the same M vertex, or a vertex sitting on an edge it is incident to) and (ii) every split-edge duplicate pair (T5.3). These are assertions with a one-line reason each, not a lemma.

**Required.**
- **Lemma L1.1 (structural `C ≡ 0`).** State precisely, in terms of the cut structure (F1: hinge edge with pin at the source vertex; split edge with both endpoints duplicated) and the three predicate types of T3 (orientation determinant of a vertex against an edge, edge–edge crossing determinant, squared distance), the complete list of candidate pairs whose harmonic satisfies `h(0) = p + q = 0` **for every X in the affine shape space `X0 + Φt`** (in fact, for every X satisfying Eq. (4); say which). Prove it from T1's closed form: `h(0) = 0` iff the predicate vanishes at the flat state, and the flat state has the two copies of every M vertex coincident and the two copies of every split edge coincident. Then prove the **converse direction you can actually prove**: for a candidate pair NOT on the list, `h(0)` is a non-zero polynomial in `t`, so `{t : h(0) = 0}` is a proper algebraic subset of the shape space (an "accidental" degeneracy). If this needs a genericity hypothesis on `X_ini` or on the graph (e.g. no three flat points collinear that are not forced to be), state it as H-L1 and say what happens without it.
- **Lemma L1.2 (class 3 is never a contact).** `h(θ) = p(1 − cos θ)` with `p ≠ 0` has constant sign on `(0, π]`. Restate and prove **without** the τ-chart, so the singularity at `θ = π` (check.md D10) never enters. Then state what class 3 *is* geometrically: which pairs have `h(0) = h′(0) = 0` structurally (the T5.3 mirror, `r = 0`), and prove that or give the counting argument that it is a codimension-1 accident among class-2 pairs, whichever is true — the measurement says 1,138 of 2.1M, which you must explain.
- **Lemma L1.3 (identically-zero harmonics).** Characterise the pairs with `p = q = r = 0` on the whole shape space (the 14,928): prove they are exactly the permanent incidences of T3.H.1/T3.H.2 — a hinge point shared by the two faces evaluated against an edge it lies on, etc. Give the list and the proof; this is what makes the identity test a *combinatorial* filter rather than a numeric one, which the paper will want to say.
- **Consequence stated as a corollary:** the certificate's `NOROOT` atom list (T5.2b.2, (T5.1c)/(T5.1d)/(T5.1e′)) can be built by combinatorial classification of pairs plus the numeric test only for accidental class-2 cases, and the accidental set has measure zero in `t`.

**Scratch checks you run yourself** (in derivations/scratch/check_l1.cpp, compiled against code/src as the existing scratch files are — see the header of derivations/scratch/check_r3.cpp for the build line): on ≥ 50 random graphs of code/apps/kill_common.hpp's population and the 8 reference tilings, for ≥ 20 random `t` each, every candidate pair's numerically detected class versus the class your lemma predicts from its combinatorial type. Report: number of pairs, agreement count, and every disagreement with its pair type. The lemma is not done until disagreements are zero or each one is explained as an accidental degeneracy at a specific `t` (and vanishes at a perturbed `t`).

## L2 — the periodic dimension formula

**What is measured** (K7 C1): `K(t)` is affine in `t`, `dim 𝒦 = rank(A) = 2·rank(D)` on 33/33 patterns, where `A` is the 4 × 2k matrix whose columns are `vec(M_j)`, `K(t) = K0 + Σ t_j M_j`, and `D` is the k × 2 matrix of period-potential increments (one row per null-space basis vector, the increment `w_t` of the face potential along the two periods). The project's earlier guess `min(4, 2·dim_null)` fails on `squares_3x3` (`dim_null = 4`, `dim 𝒦 = 0`).

**Required.**
- Set up from periodic_jacobian.hpp's derivation: `u_{f+t} = u_f + w_t + σ_f t/2`, `Q = 2 Jrot [w_h w_v]`, `K = Q P₀⁻¹`. The null space `Φ` acts on the x and y coordinates separately (t ∈ ℝ^{2k}). Show that the map `t ↦ K(t) − K0` is linear and factor it as (period-potential increment of the null vector) → (2 × 2 matrix). Identify exactly how a null vector `φ_j` (a scalar function on quotient vertices, applied to the x or the y coordinate) produces its increment pair `(w_h, w_v) ∈ ℝ²` = row `j` of `D`, and how the x- and y- copies of `φ_j` produce the two matrices `M_{2j}, M_{2j+1}`.
- **Lemma L2.1.** `dim 𝒦 = 2·rank(D)`. Prove it: the two matrices produced by the x- and y-copies of `φ_j` are `Jrot [w] P₀⁻¹` placed in the two rows (or columns — determine which), so the span of `{M}` is the image of `row(D) ⊗ ℝ²` under an injective linear map; conclude the rank doubles. Handle the case where `P₀` is singular (it never is for a lattice) by hypothesis.
- **Lemma L2.2 (what rank(D) is).** Give `rank(D)` a combinatorial meaning if you can: `D` is zero exactly on null vectors whose potential increments vanish, i.e. the shape directions that leave both periods' potential jumps unchanged. Prove at least: (a) `rank(D) ≤ 2`; (b) `rank(D) = 0` iff every null vector has zero period-potential increment, and give the `squares_3x3` case as the worked example (its 6 split cuts give `dim_null = 4` yet `K` is frozen at a rotation — explain why the increments vanish there); (c) a sufficient condition for `rank(D) = 2` in terms of the split-cut structure across the two periods, if one exists — otherwise state it as a conjecture with the 25/33 evidence.
- **Corollary.** `J(θ)` is a similarity for all θ iff `K` is one (already in the header); restate as a one-line corollary with the two linear equations in `t`, and note that it answers 2026 §5.1's "empirically".

**Scratch checks** (derivations/scratch/check_l2.cpp): rebuild `A` and `D` via the existing `periodic_jacobian` API on the 33 K7 patterns plus ≥ 100 random Voronoi tori (use `make_tiling_pattern` and the Voronoi torus generator used by kill_k7.cpp; find it there) with random `σ` from `quotient_sigma` at several seeds; check `rank(A) == 2·rank(D)` with dense SVD at tolerance 1e-10 relative to the largest singular value; report the distribution of `rank(D)` over the population and every case where equality fails. Also check the factorisation itself: reconstruct each `M_j` from row `j` of `D` by your formula and compare to the numerically assembled `M_j` (max abs error).

## Output format — derivations/lemmas.md
- Header: the standing hypotheses (import H1–H6 from core.md by name; add any new one as H-L1, H-L2 with a sentence on what fails without it).
- For each lemma: **Statement** (self-contained, in the notation of core.md §0), **Proof** (numbered steps; every step either [A] algebra a reader can verify by hand or [F] a cited fact with the source line), **What it replaces** (the sentence in core.md / KILL_REPORT.md / REPORT.md that currently states this as measured), **Check** (the scratch program, its output numbers, the build line).
- A final section "For the Checker": the three weakest steps, and the exact statements to attack.
- Tag every claim [D] / [F] / [A] / [N] as the technical report does. Never write "obviously".

## Rules
- Write derivations/lemmas.md incrementally (skeleton with all headings first).
- Do not modify derivations/core.md, code/src, or code/tests. Scratch programs go in derivations/scratch/.
- If a lemma turns out false or needs a hypothesis that the measured corpus happens to satisfy, say so in the statement; that is a result, not a failure.

## Reply
≤15 lines: the statement of each lemma in one line; proved / proved-under-hypothesis / could-not-prove for each; the scratch check numbers (pairs, agreements, disagreements; patterns, rank equalities, failures); the three weakest steps.
