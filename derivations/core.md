# Core derivation for R1 — "the usable shape space"

Deriver's output. Every numbered step is one of

* **[D]** a definition,
* **[F]** a fact cited from `STATE.md` (F-numbers), a paper (`notes/paper_2026.md`, `notes/paper_2025.md`),
  or the code's stated convention (`code/README.md`),
* **[A]** an algebraic step a Checker can verify by hand, or
* **[N]** a numeric check that was actually run (program named, result quoted).

Nothing here is asserted without one of those four tags. Claims I could not fully prove are marked
**CONJECTURE** in bold and are collected again in §T8.

Contents: §0 notation and conventions · T1 trig-linear deployment · T2 no-locking ·
T3 harmonic predicates · T4 contact calculus · T5 the usable region · T6 range maximization ·
T7 rank/periodic corrections · T8 conjecture list and Checker attack surface.

---

## Change log — round 2 (reconciliation with `derivations/check.md`)

The Checker read this file with fresh context and wrote `derivations/check.md` (38 test sets / 150 191 assertions,
its own test file `Kirigami/test/derivation_tests.jl`), raising eleven
disagreements D1–D11. Every one is resolved below. My own round-2 program is
`derivations/scratch/check_r2.jl` (run line in its header; `graphs used: 16`, 197 shape-space
samples, 2 363 380 atom evaluations); numbers quoted as **[N, R2-x]** come from it.

| # | Checker's claim | resolution | where |
|---|---|---|---|
| **D1** | the `√2` correction to K2c is wrong: in the face's own frame `max(‖x‖,‖χ‖)` is the *exact* swept radius | **CORRECTED.** The Checker is right and I verify it independently: `max_θ‖y_u−γ_f‖ = max(‖x‖,‖χ‖) = ‖x_u−x̄_f‖` to **3.11e−15** over 1 628 copies **[N, R2-A]**. My round-1 measurement was in the raw seed frame, which no pruning uses. The `√2` prescription is **withdrawn**; K2c's gate ratio is identically `1.000000` and so cannot test anything — and `results/kill/KILL_REPORT.md`'s `K2c PASS, slope 0.0000` must be withdrawn as an artefact | T4.5b |
| **D2** | H-LOC as stated is trivial; the form the `O(n)` count needs is false | **CORRECTED — H-LOC is REFUTED, not merely unproved.** Reproduced on a growing square patch: `max_f‖γ_f(θ)−x̄_f‖/r_f` = 4.62, 7.71, 10.79, 13.87, 16.95, **20.03** at patch diameters 5.66 … 19.80 **[N, R2-B]** — linear in the diameter. The words `O(n)`, "certified active set", "locality theorem" are removed from this file | T4.5, T5.2c, T8 |
| **D3** | `Θ_max = min(min_e β_e, π)`, not `min_e β_e`, on split-free patterns | **CORRECTED.** Re-measured: worst deviation **8.88e−16** on 8 split-free patterns, triangles being the case with `min β = 4π/3 > π` **[N, R2-C]** | T4.4 |
| **D4** | T5.2b is not proved to be an *inner* approximation; the `θ = 0⁺` atom is missing | **CORRECTED and the direction is now PROVED**, not downgraded: with the embeddedness atom added, Proposition T5.2b′ below proves `R(ε) ⊆ U(ε)` from Lemma T4.2 plus openness of interior overlap. Consistency check: hypothesis fires on 104/197 samples at `ε = 0.001` (45/197 samples have `Θ_max < ε`, so the implication has content), **0 violations** **[N, R2-D]** | T5.2b |
| **D5** | the printed disjunct is degree 8; a degree-≤ 4 list exists | **ADOPTED with one necessary addition of my own.** The Checker's list is exact only when `g(0) ≠ 0`. `g(0) = p+q = 0` holds *systematically* — every permanent incidence at `θ = 0` (T3.H.2) and every split-edge pair (T5.3). Measured: the printed list disagrees with direct root-finding on **71 877 of 2 363 380** atom evaluations at `ε = 0.02`, and **all 71 877** have `g(0) = 0`; after adding the deflation clause, 0 mismatches — **but that round-2 number is withdrawn in round 3 as circular**, and one more deflation is needed; see the round-3 log below **[N, round 3]** | T5.2b |
| **D6** | "split cycle ⟺ `Γ` disconnected" holds only one way | **CORRECTED** to the one-directional statement (210/210 split-cycle cases had `c(Γ)>1`; **124** cases had `c(Γ)>1` with a split forest) | T1.H.2 |
| **D7** | T1 Step 6 is far more robust than my risk assessment; T7.H is too pessimistic | **CORRECTED.** Step 6's row space equals `L`'s on 398/398 random-`σ` cases including 210 with split-cut cycles; Step 6 is demoted from the weak-step list and T7.H is softened to name exactly what fails (the *geometric*-hole reading) and what survives | T1 Step 6, T7.H, T8 |
| **D8** | `contact.jl`'s flat-centroid pruning is unsound in principle | **ADOPTED as a note.** My (T4.3) is the sound static test; the flat-centroid variant discards 19 842 pairs the sound test keeps (no observed wrong `Θ_max`, 96/96) | T4.5a |
| **D9** | `u_f` means different objects here and in `deploy_basis.jl` | **CORRECTED** by a notation note fixing this file's `u_f` and naming the header's object explicitly. No code was touched (other agents own `code/`) | §0.10 |
| **D10** | the `τ` chart is singular at `θ = π` | **ADOPTED.** Domain and the amplitude/phase fallback stated | T3.4, T3.H.4 |
| **D11** | the validity certificate has three items, not two | **ADOPTED.** Flat embeddedness is separated from deployment range | T5.1, T6.4 |

Resolved **by correction**: D1, D2, D3, D4, D6, D7, D9, D10, D11 (and D8 as a note). Resolved
**by refinement of the Checker's own proposal**: D5 — see the deflation clause and its number.
No disagreement was rebutted; on every point the Checker was right, and on D5 the correction was
incomplete rather than wrong.

---

## Change log — round 3 (the three residuals of `check.md` "Round 2")

`check.md` R2.7 left exactly three items: one substantive disagreement (R2.5) and two labelling
residues (R2.4a, R2.4b). All three are resolved below, and on all three the Checker was right.
Round 3's program is `derivations/scratch/check_r3.jl` (same 16-graph corpus, same seed, same 197
shape-space samples as round 2; build line in its header).

| # | Checker's claim | resolution | where |
|---|---|---|---|
| **R2.5** | a **third** structural class `g(0) = g′(0) = 0` (`p = −q`, `r = 0`, `h = p(1 − cos θ)`) is unnamed, and round 2's "0 mismatches" is not reproducible | **CONCEDED, and round 2's number is WITHDRAWN as circular** — `check_r2.jl` compared the deflated atom list against a reference using the *same* deflation, so it was structurally blind to class 3. New: **Lemma T5.1e** proves `h = p(1 − cos θ)` has constant sign on `(0, π)`, hence is *never* a contact event, so the class deflates by `(1 − cos θ) ∝ τ²` and its atom is the constant `true` (T5.1e′). Re-measured against a `τ`-chart-free crossing test: **0 mismatches of 2 111 336 at every `ε ∈ {0.2, 0.02, 0.001}`**, versus **151** (all class 3, `ε`-independent) without the new clause (C++ prototype: 0 of 2 145 387, versus 153) | T5.2b.2, T3.H.5 |
| **R2.4a** | (T4.1b)'s `≥` direction is measured, not derived, and is printed untagged | **CONCEDED and TAGGED.** `≤` is `[A]` (derived from `β_e ∈ 𝒞` plus the wedge argument); `≥` is `[N]` only — 8 split-free patterns at 8.88e−16 (Checker: 4 patterns, 8.9e−16, test **T4-b**). I have no proof and do not claim one; the equality now carries the tag and the sample size | T4.4 |
| **R2.4b** | the certificate's `θ = 0` overlap test is a different predicate from `EMB` (`θ = 0⁺`) and does not supply the hypothesis of T5.2b′ | **CONCEDED and REPLACED.** The certificate is now `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT` — an exact polygon–polygon test at the single explicit angle `θ₁ = ε/2`. Proposition T5.2b′ is reproved as a connectedness argument (`O` is clopen in `(0, ε)`; overlap status can change only at a candidate root) and now *implies* `EMB` instead of assuming it. Stated explicitly: the `θ = 0` test is not merely different but **degenerate**, because the two copies of every split edge coincide at `θ = 0` and adjacent faces touch along whole shared edges | T5.2b.0, T5.2b′, T5.1, T6.4 |

Deciding number for R2.5: **0 / 2 145 387 with the class-3 deflation, 153 without it, at every `ε`
tried.** No claim in this file survives round 3 unchanged that the Checker disputed.

## Change log — round 4 (the single residual of `check.md` "Round 3 (final)")

`check.md` Round 3 closed R2.5, R2.4a and R2.4b and left exactly one item, **R3.1a**: the
identically-zero orientation harmonics. No measured number changes in this round — both check
programs already exclude those pairs before classifying — what changes is that the exclusion is
**stated** and that the proof of Proposition T5.2b′ now covers the pairs it removes.

| # | Checker's claim | resolution | where |
|---|---|---|---|
| **R3.1a(a)** | (T5.1e′) is printed as "`true` when `C = B = 0`" with no side-condition, while Lemma T5.1e assumes `p ≠ 0`; when `A = C = B = 0` the harmonic is identically zero and *every* `θ` is a root, so the atom `true` is wrong there. **14 928** measured candidate harmonics are identically zero and the numeric class rule puts **all 14 928** in class 3 | **CONCEDED.** (T5.1e′) now reads "`true` when `C = B = 0` **and `A ≠ 0`**", and a new paragraph states that `A = C = B = 0` means `h_o,π ≡ 0`, that such pairs are the permanent incidences of T3.H.1/T3.H.2, and that they are **struck from the candidate list beforehand by the identity test T3.H.1** (`\|p\| + \|q\| + \|r\| ≤ tol·scale`) and are therefore **not candidates** of any class | T5.2b.2 |
| **R3.1a(b)** | T5.2b.0(iii)'s `𝒞-list` is not stated as excluding the identically-zero harmonics, so (iii) taken literally is unsatisfiable on every pattern with a hinge, and taken as implemented is weaker than its definition | **CONCEDED.** T5.2b.0(iii) now defines `𝒞-list` explicitly as the complete ordered list of T4.1b / Corollary T4.2′ **minus** `{ π : h_o,π ≡ 0 }`, removed by the identity test T3.H.1 *before* the three-class split, and says why (iii) would otherwise be unsatisfiable | T5.2b.0 |
| **R3.1a(c)** | the "`O` is closed" step is vacuous when the Lemma-T4.2 witness at `θ*` is a permanent incidence, and hinge-adjacent faces carry one at every angle | **CONCEDED and PROVED**, as new **Sub-lemma T5.2b″** inside the proof of Proposition T5.2b′. Lemma T4.2 is first localized at the limit point `z`. If some witness at `z` has `h ≢ 0` the old argument runs. Otherwise every witness is a permanent coincidence, hence a hinge point shared by `f` and `g`, and two cases remain. **Case A** (`z` is that hinge point): rigidity of the faces gives a `θ`-independent radius `r` inside which each face is its angular sector at the hinge, and the sectors' interiors meet iff `θ > β_e`; overlap at `θ_k` and none at `θ*` therefore force `θ* = β_e`, and the `β_e` collinearity of T4.4 is a **different**, non-degenerate candidate pair (far endpoint of the shorter far-side edge on the other far-side edge) whose harmonic vanishes at `θ*` and is not identically zero, because the angle between those edges is `β_e − θ`, strictly monotone. **Case B** (`z` interior to a maximal shared segment): both endpoints of that segment are then hinge points shared by `f` and `g`, so the two faces share two distinct material points; since the face motions are orientation-preserving isometries, the relative placement `γ_g^{-1}γ_f` is determined by two point images and is **constant in `θ`**, making the overlap predicate `θ`-independent — contradicting overlap at `θ_k` and none at `θ*`. So Case B cannot occur | T5.2b′ |

**The one hypothesis this round adds, stated openly.** Sub-lemma T5.2b″ uses the T3.H.1/T3.H.2
classification that `h_o,π ≡ 0` always comes from a permanent *coincidence* (`w ≡ a` or `w ≡ b`). A
permanently *collinear but non-coincident* triple would also have `h ≡ 0`, would also be struck from
`𝒞`, and would not be covered: its event is "`w` enters the segment `[a,b]`", invisible to the
orientation harmonic. **The additional candidate atoms needed in that case are named explicitly**
in T5.2b′: the two interval predicates `s₁ = ⟨w−a, b−a⟩`, `s₂ = ⟨w−b, a−b⟩` (T5.2b″-2) — *this
claim is superseded in round 5: `s₁, s₂` close only the sub-case where `w` sits at an endpoint of
`[a,b]`; see the round-5 log below* — each again
of the form `p + q cos θ + r sin θ` with coefficients quadratic in `t`, so the same three-class atom
list decides them and the degree-≤ 4 bound is unaffected. Measured status: the Checker verified
**144 / 144** of the identically-zero pairs on `squares` are coincidences; the classification of the
other 14 784 on its corpus is not reported, so for those the coincidence classification is assumed.

Deciding numbers for round 4: **unchanged** — `check_r3.jl` still reports 0 mismatches of
2 111 336 at every `ε ∈ {0.2, 0.02, 0.001}`, and the standalone `Kirigami/test/derivation_tests.jl`
still reports 2 test cases, 59 982 assertions, 0 failures. This round changes statements and one
proof, not measurements.

## Change log — round 5 (`check.md` "Round 4 (final)": R4.4 and the R4.2 corrections)

`check.md` Round 4 returned **AGREE** on R4.1, R4.2 (conclusion), R4.3 and R4.5, and one
**DISAGREE**, R4.4. Round 5 concedes R4.4 in full and prints the three R4.2 corrections. As in
round 4, **no measured number changes** — the edits are to statements and to one proof.

| # | Checker's claim | resolution | where |
|---|---|---|---|
| **R4.4** | the two interval atoms `s₁, s₂` do **not** close the permanently collinear non-coincident case: when `w` stays in the relative interior of `[a,b]` both are strictly positive and silent at exactly the transition they are meant to catch | **CONCEDED, and the claim is withdrawn and replaced by a proof.** [D] is rewritten. (i) For **hinge-adjacent** pairs the case is impossible: the relative motion is a rotation by `±θ` about the shared hinge point, so a body point of `g` traces a *circle* in `f`'s frame, and a circle lies in a fixed line only if its radius is zero, i.e. the point *is* the hinge — a coincidence. This is where all 14 928 measured identically-zero harmonics live. (ii) For the rest, the transition splits: `w` at an endpoint of `[a,b]` is closed by `s₁, s₂`; `w` in the relative interior is closed by the **neighbour-vertex collinearity pair** `(w′,(a,b))`, already in the complete list. The interior sub-case is now **proved** (new (T5.2b″-1c)): locally `f` is a half-plane and `g` an angular sector at `w`, the overlap predicate is "arc meets open half-circle", and if no bounding ray of the sector is parallel to `ab` at `θ*` the arc sits in a *compact* subset of the complementary open arc, so by continuity the predicate is locally constant — contradicting `θ_k → θ*`. Hence a bounding ray is parallel to `ab`, and since `w ∈ ab` that means `w′` is collinear with `a, b`. Non-degeneracy: if `h_(w′,(a,b)) ≡ 0` too, then `w′ − w ∥ b − a` for all `θ`, which by (0.7′) forces `σ_f = σ_g`, making the relative motion a pure translation and the local predicate `θ`-independent — again a contradiction | T5.2b′ [D] |
| **R4.2 c1** | face polygons must be **simple**; used by T4.2 and again by the sector picture, but not among the standing hypotheses. With `POS` alone, 5 of 873 local tests failed, all on self-intersecting faces | **CONCEDED.** Simplicity of every face polygon is now printed among the **standing hypotheses of T4–T6**, and a "where simplicity is used" note in the proof of Sub-lemma T5.2b″ says exactly which step consumes it (the sector description of `P_f(θ) ∩ B(p,r)`) and quotes the Checker's 5 / 873 failures and their removal by excluding the 99 non-simple hinges | §T4 hypotheses, T5.2b′ [A] |
| **R4.2 c2** | "`h_o,π′` vanishes at that one angle alone" is false | **CONCEDED.** Case A now prints `h_o,π′(θ) = ± L L′ sin(β_e − θ)` (T5.2b″-1b), with zero set `{β_e, β_e ± π}`; exactly one of the three lies in `(0, π)`, and it is `β_e` precisely when `β_e ∈ (0, π)`, which is the case at hand since `θ* = β_e ∈ (0, ε)`. Only `β_e` is a **transition of the sector-overlap predicate** — at `β_e ± π` the far-side gap is `∓π`, the far edges are anti-parallel and `θ > β_e` does not change value. The proof needs only `h_o,π′ ≢ 0` and `h_o,π′(β_e) = 0`; an extra root only makes `NOROOT` stricter, hence `R(ε)` smaller and still inner | T5.2b′ Case A |
| **R4.2 c3 / R4.3** | the matching argument (a face keeps a vertex for at most one of its hinge edges) is used silently, and it makes Case B vacuous | **CONCEDED and PRINTED** as **Sub-lemma T5.2b‴**: by 0.2 a hinge edge keeps `v` only at `src(e)`, which in the `σ_f`-traversal of `∂f` is the edge *leaving* `v`, and only one edge leaves `v`; so vertex identification is a *matching*, never a fan, and two faces sharing a hinge point are joined by exactly one hinge edge — which is what makes "the" `β_e` of the pair well defined. **Case B is now stated as vacuous** (Remark [C]): hinge adjacency gives `σ_g = −σ_f`, so the relative map is a rotation by `±θ` (T1.C), which for `θ ∈ (0, π]` fixes exactly one point and cannot fix two — two faces never share two distinct material points. Case B is kept as a remark | T5.2b′ [A0], [C] |

Deciding numbers for round 5: **unchanged**. `check_r3.jl` still reports 0 mismatches of
2 111 336 at every `ε ∈ {0.2, 0.02, 0.001}` (class sizes 2 037 472 / 72 798 / 1 138), and
`Kirigami/test/derivation_tests.jl` re-runs at **38 test sets / 150 191 assertions, 0 failures**.
Round 5, like round 4, changes statements and proofs, not measurements.

---

## §0. Notation, and the sign convention actually used by the code

**0.1 [D] The uncut graph.** `M = (V, E, F)` is a straight-line embedded planar graph with `N = |V|`,
vertex positions `X ∈ R^{N×2}` (row `x_v`), faces stored counter-clockwise in the given geometry.
`E = E_border ⊔ E_hinge ⊔ E_split`.

**0.2 [F, 2026 §3 / F1 / `code/README.md`] Orientation and cut type.** `σ : F → {−1, +1}`, with
`σ(f) = +1` meaning "traversed clockwise". An interior edge with `σ(f₁) = −σ(f₂)` is a **hinge**
edge; with `σ(f₁) = σ(f₂)` a **split** edge. A hinge edge carries a direction `src(e) → dst(e)`
induced by the two (agreeing) `σ`-induced half-edge directions; the cut **keeps the hinge at
`src(e)`** and duplicates `dst(e)`. Both endpoints of a split edge are duplicated.

**0.3 [D] The kirigami structure `M′`.** Vertices of `M′` are equivalence classes of face-corners
(`cut.jl`: `corner_to_prime`); write `(v, f)` for the copy of `v ∈ V` carried by face `f`, and
`n′ = |V′|`. Faces of `M′` are the faces of `M` with their corners relabelled.

**0.4 [D] The hinge graph.** `Γ = (F, E_hinge)`: one node per face, one edge per hinge edge joining
its two incident faces. `c(Γ)` = number of components, `b₁(Γ) = |E_hinge| − |F| + c(Γ)`.

**0.5 [A] `Γ` is bipartite, `σ` is a proper 2-colouring.** A hinge edge joins faces of opposite `σ`
by 0.2, so `σ` is a proper 2-colouring of `Γ`. In particular every cycle of `Γ` has even length.
(This is the whole reason a *common* opening angle can exist; it is geometer (0.1).)

**0.6 [D] Rotation.** `J = [[0, −1], [1, 0]]` (rotation by `+π/2`), `R(a) = cos(a) I + sin(a) J`.
Throughout, `c := cos(θ/2)`, `s := sin(θ/2)`, and `τ := tan(θ/2)`.

**0.7 [F, `Kirigami/src/core/kinematics.jl`] The code's sign convention.** Face `f` is transformed by
`y = R(a_f) x + t_f` with

```
        a_f  =  − σ(f) · θ / 2 ,          t_f(0) = 0 .                              (0.7)
```

Because `σ(f) = ±1`, `cos(σ a) = cos a` and `sin(σ a) = σ sin a`, this is

```
        R_f  =  R(−σ_f θ/2)  =  c · I  −  σ_f s · J .                               (0.7′)
```

**0.8 [A] Equivalence with the geometer's convention, which is the opposite one.**
`ideas/persona_geometer.md` (0.2) derives `ω_f = +σ(f) θ/2`. The two differ by `θ ↦ −θ`. Both
satisfy "every hinge opens by `|θ|`" (T1.C below), so kinematics alone does not fix the sign; what
fixes it is which side of each hinge is required to *open* (the hole side) rather than close. The
code's choice is the one for which the hole preimages open, and it is the one all measured numbers
in `results/` use. **Everything below is in the code convention (0.7).** To translate any formula in
`ideas/persona_geometer.md` or `ideas/persona_rigidity.md` into this file, substitute `θ → −θ`,
equivalently `s → −s`, equivalently flip the sign of every `S`/`χ`/`u` vector. `C` is unaffected.
This is not a free relabelling `σ → −σ`: flipping `σ` globally also reverses every hinge direction
(0.2), producing a *different* cut structure, so the sign must be tracked, not argued away.

**0.9 [D] Shape space.** `𝕏 = {X : Eq. (2) holds for every hole preimage} = X₀ + span{Φ}`, with
`Φ ∈ R^{N×k}` the null-space basis of `[L; B]` (Eq. (5)) and `k = dim_null`. Design coordinates
`t ∈ R^{k×2}`, `X(t) = X₀ + Φ t`. `[F15]` `k = |E_split|` on 100/100 measured graphs with fixed
boundary.

**0.10 [D] The symbol `u_f`, and the clash with `deploy_basis.jl` (D9).** In this file `u_f ∈ R²`
is the **face potential** of T1 Step 3, defined by `t_f = 2 sin(θ/2) · J u_f`; it is `θ`-free and
satisfies the increment law (T1.3). The header comment of `Kirigami/src/method/deploy_basis.jl` writes
`S_pv = −σ_f J x_v + u_f` and by `u_f` means the whole translation *direction*, i.e. the vector
`ν_f` with `t_f = sin(θ/2) · ν_f` (I avoid the letter `τ`, which 0.6 has already spent on
`tan(θ/2)`). The two are related by

```
        ν_f  =  2 J u_f  ,        equivalently   u_f = ½ Jᵀ ν_f .                     (0.10)
```

The code is correct — `S` is computed as `2 ∂Y/∂θ|₀ = J(2u_f − σ_f x_v)`, verified to `7.3e−15`
(Checker **T1-b**) — only the *name* is overloaded. **Everywhere below, `u_f` is this file's
object (0.10-left).** I did not edit `code/` (other agents own it); whoever next touches
`deploy_basis.jl` should rename its comment variable to `ν_f` before either symbol becomes a
paper formula.

---

## T1 — Trig-linear deployment

**Hypotheses.** `M` embedded planar; `σ` given; `Γ` connected (`c(Γ) = 1`; see T1.H for
`c(Γ) > 1`); `X` satisfies Eq. (2) on a cycle basis of `Γ`. Faces are rigid. No genericity is used.

### T1.1 The hinge constraint is a discrete potential equation on `Γ`

**Step 1 [D].** A *uniform deployment at angle `θ`* is an assignment of plane isometries
`T_f(x) = R_f x + t_f`, one per face, with `T_f(0) = id` at `θ = 0`, such that for every hinge edge
`e` between faces `f, g` the shared hinge point is preserved:

```
        T_f( x_{src(e)} )  =  T_g( x_{src(e)} ) ,                                    (T1.1)
```

and such that every hinge opens by the same angle `θ`. By 0.7 the second requirement is exactly
(0.7′); only (T1.1) remains to be solved.

**Step 2 [A].** Substituting (0.7′) into (T1.1) with `v := src(e)` and using `σ_g = −σ_f`:

```
        t_g − t_f  =  (R_f − R_g) x_v  =  (−σ_f + σ_g) s · J x_v  =  2 s σ_g · J x_v .   (T1.2)
```

**Step 3 [D].** Since `2 s J` is invertible for `θ ∈ (0, 2π)`, define the **face potential**
`u : F → R²` by `t_f = 2 s · J u_f`. Then (T1.2) becomes a `θ`-free prescription of increments on
the edges of `Γ`:

```
        u_g − u_f  =  σ_g · x_{src(e)} ,        for the hinge edge e traversed f → g .   (T1.3)
```

**Step 4 [A] (T1.3 is a well-posed edge increment).** Traversing `e` the other way gives
`u_f − u_g = σ_f x_{src(e)} = −σ_g x_{src(e)}`, the negative of (T1.3). So (T1.3) defines a genuine
1-form `δ` on the oriented edges of `Γ`, and `δ` does not depend on `θ`.

**Step 5 [A] (existence ⇔ closure).** A potential `u` with prescribed increments `δ` exists on a
connected graph iff `Σ_{e ∈ z} ε_e δ_e = 0` for every cycle `z` of `Γ` (`ε_e = ±1` the traversal
sign). Equivalently, iff it vanishes on any cycle basis. Because `δ` is `θ`-free, **the existence of
a uniform deployment at one `θ ≠ 0` is equivalent to its existence at every `θ`**, and the condition
is a linear condition on `X` alone.

**Step 6 [A] (closure ⇔ Eq. (2)).** Take the face-cycle basis of the planar graph `Γ`: the bounded
faces of `Γ`. Let `z` be one, with faces `f₁, …, f_{2m}` in cyclic order (even by 0.5) and hinge
edges `e₁, …, e_{2m}`, `e_i` between `f_i` and `f_{i+1}`, all traversed forward. Then

```
   Σ_i (u_{f_{i+1}} − u_{f_i})  =  Σ_i σ_{f_{i+1}} x_{src(e_i)}  =  σ_{f₂} Σ_i (−1)^{i−1} x_{src(e_i)} .
```

Consecutive `σ` alternate, so the sum is an **alternating** sum of the hinge source positions around
the cycle. Now use `[F, 2026 Prop. 4.1 / Eq. (2)]`: for the hole preimage `C` the condition is
`Σ_{e ∈ C ∩ E_hinge} ( x_{dst(e)} − x_{src(e)} ) = 0`, a sum of *directed hinge edge vectors*. The
two agree because around a bounded face of `Γ` the hinge edges alternate in/out at the enclosed
vertex set (`[F, Remark A.1]`: `indeg_h(v) = outdeg_h(v)` at every interior vertex), so consecutive
terms pair up as `x_{dst} − x_{src}`. Concretely, in the pure-hinge case where the bounded face of
`Γ` surrounds one interior vertex `v` of degree `k_v`, the alternating sum telescopes to
`k_v x_v − Σ_{u → v} x_u`, which is Eq. (2) for the preimage `C_v = {in-edges of v}`
(`[F, F11]`). Deleting a split-dual edge from `Γ` merges two of its bounded faces, which is exactly
"split cuts merge holes" (`[F, 2026 §4.2]`), and the corresponding alternating sums add. Hence:

```
   closure of δ on Γ   ⟺   Eq. (2) on every hole preimage   ⟺   L X = 0 .            (T1.4)
```

*(This step reproduces geometer (0.5). It is the one step in T1 that leans on F11 — the statement
that hole preimages are exactly the bounded faces of `Γ` — which is proved by planar duality in
geometer (0.5) and measured 100/100 in F15.)*

**[N] Round-2 update (D7): the attack on this step was run and it failed.** In round 1 I ranked
Step 6 among the three weakest steps and proposed the attack "exhibit an `(M, σ)` where the
face-cycle basis of `Γ` and `holes_partition` disagree, e.g. with a split-cut cycle". The Checker
ran exactly that. It built the closure rows symbolically (BFS potential, one row
`potvec[g] − potvec[f] − σ_g e_{src(e)}` per non-tree hinge edge, giving `G ∈ R^{b₁(Γ)×N}`) and
compared **row spaces** with `L` via `rank(G) = rank(L) = rank([G;L])`. Identical on **16/16**
corpus cases and **398/398** random-`σ` cases — among them 334 with `c(Γ) > 1` and **210 with
split-cut cycles** (`check.md` **CE-c**); `b₁(Γ) = H` on all of them, `holes_seed_growing` agreed
with `holes_partition` 398/398, and the preimages partitioned `E_hinge ⊔ E_split` 398/398, so F11
survives fully random `σ`. Step 6 is therefore **removed from the weak-step list** in T8. It still
*uses* F11; what changed is that the proposed failure mode has been searched for and does not
occur.

**[N] Path-independence check.** `derivations/scratch/check_t1_t2.jl`, check **C1**: build `u` by BFS
over a spanning tree of `Γ`, then evaluate the increment residual on **every** hinge edge (not only
tree edges), on 30 graphs (7 tilings + Delaunay/Voronoi/quad-random, 3 seeds each) at the solved
`X₀`. Worst residual `1.60e−14`.

### T1.2 The closed form

**Step 7 [A].** With `t_f = 2 s J u_f` and (0.7′), the deployed position of the copy `(v, f)` is

```
   y_{(v,f)}(θ)  =  R_f x_v + t_f
                 =  c · x_v  +  s · J ( 2 u_f  −  σ_f x_v ) .
```

Therefore, writing `C, S ∈ R^{n′×2}`,

```
   ┌──────────────────────────────────────────────────────────────────────────┐
   │   Y_θ  =  cos(θ/2) · C  +  sin(θ/2) · S ,                                │
   │   C_{(v,f)} = x_v ,      S_{(v,f)} = J ( 2 u_f − σ_f x_v ) .      (T1.5) │
   └──────────────────────────────────────────────────────────────────────────┘
```

**Step 8 [A] (`C` and `S` are linear in `X`).** `C` is the duplication map `X ↦ X′`, manifestly
linear. `u_f` is, by Step 5 applied along a spanning tree `T` of `Γ` rooted at `f₀`, the **path sum**

```
        u_f  =  Σ_{e ∈ path_T(f₀ → f)}  σ_{head(e)} · x_{src(e)} ,                   (T1.6)
```

a fixed signed sum of rows of `X` with coefficients in `{0, ±1}` determined by `Γ`, `σ` and `T`
alone. Hence `u` is linear in `X`, and so is `S = J(2u_{·} − σ_{·} x_{·})`. By Step 5/6 the value of
`u_f` does not depend on the tree, provided `X ∈ 𝕏`.

**[N] Closed form vs. the code.** Check **C2** of the same program: (T1.5), with `u` built by my own
BFS, against `kinematics.jl::deploy()` at `θ ∈ {0, 0.05, 0.3, 0.77, 1.2, 1.9, 2.6, 3.0}` on the same
30 graphs. Worst deviation `1.42e−14` over all vertices, all angles, all graphs. This simultaneously
verifies the sign convention 0.7–0.8; the opposite sign fails by `O(1)`.

### T1.3 Corollaries

**T1.A [A] Every vertex traces a centred ellipse.** Fix the copy `(v, f)` and let
`A = [ C_{(v,f)} | S_{(v,f)} ] ∈ R^{2×2}`. Then `y(θ) = A · (c, s)ᵀ` with `(c, s)` on the unit
circle, so as `θ` ranges over `R` the point `y` traces the image of the unit circle under `A`: an
ellipse **centred at the origin of the frame in which the seed face has `t_{f₀} = 0`**. Its
semi-axes are the singular values of `A`,

```
   a², b²  =  ½ [ (|C|² + |S|²)  ±  √( (|C|² − |S|²)²  +  4 (C·S)² ) ] ,             (T1.7)
```

with major axis along the left singular vector. It degenerates to a segment iff `det A = 0`, i.e.
`C ∥ S`; it is a circle iff `|C| = |S|` and `C · S = 0`. When `A` is invertible the trajectory
satisfies `yᵀ (A Aᵀ)^{-1} y = 1`. Physical range: `θ ∈ [0, π)` covers a quarter of the ellipse.

**[N]** Check **C7**: `|yᵀ (A Aᵀ)^{-1} y − 1| ≤ 1.68e−10` over all non-degenerate copies, all angles,
30 graphs (the `1e−10` rather than `1e−14` is the conditioning of `(A Aᵀ)^{-1}` on near-degenerate
ellipses, which are skipped only below `|det A| < 1e−8`).

**T1.B [A] Split-cut duplicates stay parallel — one line, replacing 2026 Appendix A.**
Let `e = {a, b}` be a split edge between `f` and `g`, so `σ_f = σ_g`. From (T1.5),

```
   y_{(b,f)} − y_{(a,f)}  =  c (x_b − x_a) − σ_f s J (x_b − x_a)  =  R_f (x_b − x_a) ,
```

because the `2 u_f` terms cancel in the difference. The same computation for `g` gives
`R_g (x_b − x_a)`, and `R_f = R_g` since `σ_f = σ_g`. So the two duplicates are the **same vector**
at every `θ` — not merely parallel, but a pure translate. The offset is

```
   y_{(a,g)} − y_{(a,f)}  =  2 s · J ( u_g − u_f )  =  2 sin(θ/2) · J Δu ,   Δu := u_g − u_f ,  (T1.8)
```

which is used verbatim in T5.3. 2026 Appendix A spends Remarks A.2–A.4, Eqs. (15)–(24) and an
induction over the split forest on this statement.

**[N]** Check **C4**: `‖d₁ − d₂‖ ≤ 2.08e−14` over every split edge, every angle, 30 graphs.

**T1.C [A] Every hinge opens by exactly `θ`.** For a hinge edge `e` between `f, g`, both duplicates
of `e` emanate from the common hinge point `y_{(src,f)} = y_{(src,g)}` and equal `R_f (x_{dst} −
x_{src})` and `R_g (x_{dst} − x_{src})`. The signed angle between them is
`a_f − a_g = −(σ_f − σ_g) θ/2 = σ_g θ`, of magnitude `θ`. This is the consistency of (T1.5) with the
definition of "uniformly deployable", i.e. the sign convention 0.7 is not merely a convention but is
*forced* once the hinge is required to open by `θ` and the hinge point is required to be fixed.

**[N]** Check **C3**: `| |signed hinge angle| − θ | ≤ 9.73e−14`, every hinge edge, every angle,
30 graphs.

**T1.D [A] Faces are rigid and orientation-preserving at every `θ`.** `T_f` is a proper isometry
(`det R_f = 1`), so every face polygon of `M′` is congruent to its flat self and has the **same
signed area at every `θ`**. In particular positive orientation of the faces is a property of `X`
alone, checked at `θ = 0`. This is used as a hypothesis in T3, T4 and T5.

**[N]** Check **C5**: `|area_f(θ) − area_f(0)| ≤ 4.46e−14`, every face, every angle, 30 graphs.

### T1.H What can break

1. **`c(Γ) > 1`.** Then `u` is determined only up to one additive constant *per component*, i.e.
   `M′` falls apart into `c(Γ)` independently placed pieces, and (T1.5) holds within each component
   with its own seed. The relative placement of components is not determined by the hinges — a real
   `2(c(Γ) − 1)`-dimensional slack that is not a mechanism of the assembled structure. The 2026
   connectivity principle (§4.2 (1)) forbids it; `[F, F22]` `c(Γ) = 1` on 50/50 auto-oriented graphs.
2. **Split-cut cycles.** `E_split` contains a cycle of `M` **⟹** `Γ` is disconnected, so this is
   case 1 in disguise and the fold into case 1 is valid. **The converse is false (D6):** measured
   over 398 random-`σ` pairs, **210/210** split-cycle cases had `c(Γ) > 1`, `0` had a split cycle
   with `Γ` connected — but **124** pairs had `c(Γ) > 1` with a split *forest*. Round 1 stated this
   as "iff" (attributing it to geometer's second corollary); only the forward implication is
   available and only the forward implication is used. `[F, F16]` confirms Def. 4.2 breaks on split
   cycles.
3. **`X ∉ 𝕏`.** Then `u` does not exist; `deploy()` still returns a BFS answer, and
   `Deployment::max_mismatch` is exactly the failure of Step 5. Nothing below applies.
4. **`θ = 0`.** (T1.5) is fine, but `t_f = 2sJu_f` is a *definition* only for `s ≠ 0`; at `θ = 0` all
   copies coincide and `M′` is `M`. Statements about "the deployed structure" at `θ = 0` are
   statements about `M`.

### T1.Check (for the Checker)

Re-run `derivations/scratch/check_t1_t2.jl`. Tolerance `1e−12` on C1–C6, `1e−9` on C7. Additionally:
flip the sign in `S` (use `−J(2u − σx)`) and confirm C2 fails at `O(1)` — that is the test that the
convention statement 0.8 is doing work.

---

## T2 — No-locking

**Statement.** For every `X ∈ 𝕏` and every `θ ∈ R`, the uniform-deployment velocity field lies in
the kernel of the infinitesimal-closure operator `A(Y_θ)`. Consequently the uniform branch is a
real-analytic 1-parameter motion of the body-and-pin framework defined for all `θ`, with nowhere
vanishing non-rigid velocity, and **it never reaches a kinematic dead centre: every termination of
deployment is a contact event.**

**T2.1 [D] The closure operator.** Model (`[F, F13]`, `ideas/persona_rigidity.md` (R1)–(R2)): faces
are rigid bodies, each hinge edge `e` is a pin joint at the point `p_e(θ) := y_{(src(e), f)} =
y_{(src(e), g)}`; split edges impose nothing. An infinitesimal motion assigns face `f` an angular
velocity `ω_f ∈ R` and a spatial velocity field `V_f(y) = ω_f J y + w_f`. The pin constraint at `e`
between `f, g` is `V_f(p_e) = V_g(p_e)`, i.e.

```
        (w_f − w_g)  +  (ω_f − ω_g) J p_e(θ)  =  0 .                                  (T2.1)
```

Given `ω`, a compatible `w` exists iff the 1-form `e ↦ (ω_{h(e)} − ω_{t(e)}) J p_e(θ)` is exact on
`Γ`, i.e. iff its circulation vanishes on a cycle basis. That circulation map is the
`2 b₁(Γ) × |F|` matrix

```
        A(θ)_z (ω)  =  Σ_{e ∈ z} ε_e ( ω_{h(e)} − ω_{t(e)} ) J p_e(θ) ,               (T2.2)
```

one basis cycle `z` contributing two scalar rows. *These are the "closure constraints".* Mobility is
`m = dim ker A − c(Γ)` (`[F, U4]`, self-tested 21/21 by the rigidity persona; not re-derived here).

**T2.2 [A] `A(θ) σ = 0`, by explicit construction of `w`.** Take

```
        ω_f  =  d a_f / dθ  =  − σ_f / 2 ,
        w_f  =  d t_f/dθ  −  ω_f J t_f  =  c · J u_f  −  σ_f s · u_f .               (T2.3)
```

(The second line: `t_f = 2 s J u_f`, so `dt_f/dθ = c J u_f`, and `J t_f = 2 s J² u_f = −2 s u_f`,
hence `−ω_f J t_f = −(−σ_f/2)(−2 s u_f) = −σ_f s u_f`.) Now evaluate (T2.1) at a hinge edge `e`
between `f, g` with `v = src(e)`, using only `u_g − u_f = σ_g x_v` (T1.3) and `σ_f = −σ_g`:

```
   w_f − w_g  =  c J (u_f − u_g)  −  s ( σ_f u_f − σ_g u_g )
              =  − c σ_g J x_v  +  s σ_g ( u_f + u_g ) ,

   ( ω_f − ω_g ) J p_e  =  σ_g · J [ c x_v + s J ( 2 u_g − σ_g x_v ) ]
                        =  σ_g c J x_v  −  2 s σ_g u_g  +  s x_v .
```

Adding, the `c J x_v` terms cancel and

```
   (w_f − w_g) + (ω_f − ω_g) J p_e  =  s σ_g ( u_f + u_g − 2 u_g )  +  s x_v
                                    =  s σ_g ( u_f − u_g )  +  s x_v
                                    =  − s x_v  +  s x_v  =  0 .                       (T2.4)
```

Since `ω = −σ/2` is proportional to `σ`, and the pin equation holds at **every** hinge edge, the
circulation (T2.2) vanishes on every cycle. Hence `σ ∈ ker A(Y_θ)` for **every** `θ ∈ R`. ∎

Note what the proof used: only (T1.3), i.e. only `X ∈ 𝕏`. No genericity, no connectivity beyond what
makes `u` exist, no condition on `θ`.

**T2.3 [A] The pencil.** `p_e(θ) = c x_{src(e)} + s J χ_e` with `χ_e := 2 u_{h(e)} − σ_{h(e)} x_{src(e)}`
constant, and (T2.2) is linear in `p`. Therefore

```
        A(θ)  =  cos(θ/2) · A_c  +  sin(θ/2) · A_s ,                                  (T2.5)
```

a linear matrix pencil in `τ = tan(θ/2)` after clearing `c`. Two consequences, both handed to R3
rather than used here: `rank A(θ)` is constant off a finite set of `τ` (the real points of the
determinantal variety of `(A_c, A_s)`), and it attains its **minimum** exactly on that set, which
contains `τ = 0`. That is the precise form of the rigidity persona's Idea 2 (`m(0) > m(0⁺)`), and it
explains why 2026 Eq. (9) — solved at `θ = 0` — is solved at the one configuration where the
linearisation is least representative.

**T2.4 [A] Non-triviality and the corollary.** The velocity field (T2.3) is not a rigid motion of the
whole structure whenever `σ` is non-constant on some component of `Γ` — i.e. whenever `E_hinge ≠ ∅` —
because two hinged faces have distinct `ω`. So `m ≥ 1` at every `θ ∈ R`, and (T1.5) is an explicit
real-analytic motion defined on all of `R`. **Therefore the uniform branch has no kinematic
termination**: the only way deployment stops is that two faces that are not joined by a hinge come
into contact. This is what makes the candidate list of T4 *complete* — there is nothing else to check.

**[N]** Check **C6**: the residual of (T2.4) evaluated numerically, at every hinge edge, at
`θ ∈ {0, 0.05, 0.3, 0.77, 1.2, 1.9, 2.6, 3.0}`, on 30 graphs: max `1.83e−14`.

**T2.H What can break.**
1. The claim is about the *pin-joint* model. If a fabrication model gives hinges finite width, or
   faces thickness, a "lock" can occur that this model does not see. The theorem is about the
   idealisation both papers use (`[F, 2026 §7 (iv)]`: rigid planar faces assumed).
2. `A` here is the **infinitesimal** closure operator. `σ ∈ ker A` at every `θ` is a statement about
   first-order flexes along a curve we have already exhibited; it does **not** say the configuration
   space is a smooth 1-manifold at `Y_θ`. Other branches may cross the uniform one (T2.3), and
   `dim ker A(Y_θ)` may exceed 1. Nothing in T4/T5 needs uniqueness of the branch — only that the
   *uniform* branch is what is being deployed, which is the papers' setting.
3. If `c(Γ) > 1` the count `m = dim ker A − c(Γ)` and the framework itself split per component; T2.2
   still holds edge-by-edge.

**T2.Check.** Assemble `A(Y_θ)` from a spanning-tree cycle basis at 10 random `θ ∈ (0, π)` on 50
graphs and check `‖A(Y_θ) σ‖_∞ ≤ 1e−10 · ‖A‖_∞`. Then check (T2.5) by assembling `A_c, A_s` once and
comparing `cos(θ/2) A_c + sin(θ/2) A_s` against a direct assembly at each `θ` (should be exact to
round-off). Independent of `check_t1_t2.jl`'s C6, which tests the pin equations rather than the
assembled matrix.

---

## T3 — Harmonic predicates

**Statement.** Along the uniform branch every orientation, dot-product and squared-length predicate
built from `M′`-vertices is a *first harmonic* `h(θ) = p + q cos θ + r sin θ` whose coefficients are
**quadratic forms in `X`**.

**T3.1 [A] The three coefficient formulas.** For `M′`-vertices `a, b, w` write
`C_{ab} = C_b − C_a`, `S_{ab} = S_b − S_a` (differences of the rows of (T1.5)). Then
`Y_b − Y_a = c C_{ab} + s S_{ab}`. Using

```
        c² = (1 + cos θ)/2 ,      s² = (1 − cos θ)/2 ,      c s = (sin θ)/2 ,          (T3.1)
```

and bilinearity of `det(·,·)` and `⟨·,·⟩`:

```
   det( Y_b − Y_a , Y_w − Y_a )  =  p + q cos θ + r sin θ  with
        p = ½ [ det(C_ab, C_aw) + det(S_ab, S_aw) ]
        q = ½ [ det(C_ab, C_aw) − det(S_ab, S_aw) ]                                    (T3.2)
        r = ½ [ det(C_ab, S_aw) + det(S_ab, C_aw) ]

   ⟨ Y_b − Y_a , Y_w − Y_a ⟩     =  p′ + q′ cos θ + r′ sin θ  with
        p′ = ½ [ ⟨C_ab, C_aw⟩ + ⟨S_ab, S_aw⟩ ]
        q′ = ½ [ ⟨C_ab, C_aw⟩ − ⟨S_ab, S_aw⟩ ]                                         (T3.3)
        r′ = ½ [ ⟨C_ab, S_aw⟩ + ⟨S_ab, C_aw⟩ ]

   ‖ Y_b − Y_a ‖²                =  p″ + q″ cos θ + r″ sin θ  with
        p″ = ½ ( |C_ab|² + |S_ab|² ) ,  q″ = ½ ( |C_ab|² − |S_ab|² ) ,  r″ = ⟨C_ab, S_ab⟩ .  (T3.4)
```

The **signed area** of the triangle `(a, b, w)` is half of (T3.2).

**T3.2 [A] The coefficients are quadratic forms in `X`.** `C` and `S` are linear in `X` (T1 Step 8),
and every coefficient in (T3.2)–(T3.4) is a bilinear expression in two of them. So each of
`p, q, r` is a real quadratic form in the `2N` entries of `X`; composing with `X = X₀ + Φ t` makes
each a quadratic polynomial in the design coordinates `t`. This — not the harmonic form itself, which
is the Weierstrass substitution — is the load-bearing fact for T5.

**T3.3 [A] Face signed areas are constant in `θ`.** By T1.D. Equivalently: for `a, b, w` all in the
same face `f`, `S_{ab} = J R_?`… concretely `C_{ab} = x_b − x_a` and
`S_{ab} = −σ_f J (x_b − x_a) = −σ_f J C_{ab}`, so `det(S_ab, S_aw) = det(JC_ab, JC_aw) =
det(C_ab, C_aw)`, giving `q = 0`; and `det(C_ab, S_aw) = −σ_f det(C_ab, J C_aw)`,
`det(S_ab, C_aw) = −σ_f det(J C_ab, C_aw) = +σ_f det(C_ab, J C_aw)` (since
`det(Ja, b) = −det(a, Jb)`), so `r = 0`. Hence `h ≡ p = det(C_ab, C_aw)`, the flat value. ∎
The same computation gives `q″ = r″ = 0` for two vertices of one face: **edge lengths within a face
are constant**, which is used in T4.2 to reduce one of the two interval harmonics to a constant.

**T3.4 [A] Roots, in closed form.** Substituting `cos θ = (1 − τ²)/(1 + τ²)`,
`sin θ = 2τ/(1 + τ²)` with `τ = tan(θ/2)`:

```
        (1 + τ²) h(τ)  =  g(τ)  :=  (p + q)  +  2 r τ  +  (p − q) τ² .                (T3.5)
```

So `h` has at most two roots in `(−π, π]` and they are the roots of a **quadratic in `τ`**:

```
   p ≠ q :   τ_{1,2}  =  [ − r  ±  √( q² + r² − p² ) ] / ( p − q ) ,   real iff  p² ≤ q² + r² ;
   p = q :   one finite root  τ = −(p+q)/(2r) = −p/r   (if r ≠ 0),  plus a root at θ = π ;
   p = q, r = 0 :  h ≡ 2p cos²(θ/2), root only at θ = π (or h ≡ 0 if p = 0).           (T3.6)
```

`θ = 0` is a root iff `p + q = 0`. The real-root condition `p² ≤ q² + r²` is "`|p| ≤` amplitude",
matching `deploy_basis.jl::harmonic_roots`. Each of `p, q, r` is quadratic in `t`, so the
discriminant `q² + r² − p²` is a **quartic** in `t`.

**Domain of the `τ` chart (D10).** `τ = tan(θ/2)` is a diffeomorphism `(−π, π) → R` and does **not
reach `θ = π`**: the point `θ = π` is `τ = ∞`, and it is a root of `h` exactly when the leading
coefficient `p − q` vanishes. So (T3.6) is stated correctly only with the branch on `p = q`, and it
is **numerically singular** near `p − q ≈ 0`: the Checker found 54 of 5 404 sampled harmonics with
`|p − q| ≤ 1e−9 · scale`, where the two parametrizations disagree about whether `π` is a root by a
few `1e−3`. Since the search range of deviation 9 ends at `π`, this decides whether a contact at
the very end of the range is counted. **Rule: on `θ ∈ (0, π − δ]` use the `τ`-quadratic; for the
endpoint use the amplitude/phase form** `h(θ) = p + Ρ cos(θ − φ)` with `Ρ = √(q²+r²)`,
`φ = atan2(r, q)`, and test `h(π) = p − q` directly against a scale-relative tolerance. This is a
caveat about the chart, not an error in (T3.5)/(T3.6).

**T3.H What can break.**
1. **Degenerate triples.** If `a, b, w` are collinear at every `θ` (e.g. three vertices of a
   straight-line "face", or `b = a`), `h ≡ 0` and the root set is all of `R`. These must be *removed
   from the candidate list*, not handled numerically; a tolerance test `|p| + |q| + |r| ≤ tol · scale`
   is how the code does it (`harmonic_roots(..., tol)`).
2. **Permanent incidences.** A hinge point is shared by two faces at every `θ`, so the corresponding
   vertex-vertex predicate is identically zero. Same for the two copies of a split-edge endpoint at
   `θ = 0` only. See T4.3.
3. `p = q` exactly is a positive-codimension case but occurs *systematically*, not accidentally, for
   split-edge pairs (T5.3 has `p + q = 0`, its mirror image). The `τ`-quadratic must be solved with
   the leading-coefficient branch, never by the quadratic formula alone.
4. **`θ = π` is outside the `τ` chart (D10).** See the domain paragraph above: near `p = q` the
   root reported at the end of the range is chart-dependent at the `1e−3` level. Use the
   amplitude/phase form there.
5. **`g(0) = p + q = 0` is systematic, not accidental (D5).** It says the triple `(a, b, w)` is
   collinear in the *flat* state, which is automatic whenever `w` and `a` (or `w` and `b`) are two
   `M′`-copies of the **same** vertex of `M` — every permanent incidence of T3.H.2, and every
   split-edge pair by T5.3. Measured: **73 446 of 2 145 387** candidate harmonics **[N, round 3]**.
   Any interval predicate on the roots must therefore *deflate* the factor `τ` before testing, or
   it will report the root at `θ = 0` as a contact inside `(0, T)`.
   **Round 3: the sub-class `g′(0) = 2r = 0` as well** (`p = −q`, `r = 0`, `h = p(1 − cos θ)`) is
   systematic too and needs the factor deflated **twice** — `1 138 of 2 145 387` harmonics; by
   Lemma T5.1e such an `h` has constant sign on `(0, π)` and is never a contact. See T5.2b.2.

**T3.Check.** This is exactly kill experiment **K1b**. In addition to K1b's least-squares fit, check
the *closed-form* coefficients (T3.2) against the fit — K1b as specified only fits, so it would pass
even if (T3.2) had a sign error. Suggested: sample `h` at 200 angles, compare to `p + q cos + r sin`
with `p, q, r` from (T3.2), relative residual `≤ 1e−12`. Also assert `q = r = 0` for intra-face
triples to `1e−12` (T3.3).

---

## T4 — Contact calculus

**Standing hypotheses for T4–T6.** `X ∈ 𝕏`; **every face polygon of `M` is simple** (no
self-intersection; see T4.H.4 and the "where simplicity is used" note in the proof of Sub-lemma
T5.2b″ — added in round 5, it was used silently before); every face of `M` is positively oriented at
`θ = 0` (so, by T1.D, at every `θ`); `c(Γ) = 1`. Deployment runs over `θ ∈ (0, π]` (the code's search range,
deviation 9).

### T4.1 Definitions, and what "first contact" must mean

**[D] T4.1a.** For `θ > 0` let `P_f(θ) ⊂ R²` be the closed polygon of face `f` in `M′`. Define

```
   Θ_max(X)  :=  sup { Θ ∈ (0, π] :  int P_f(θ) ∩ int P_g(θ) = ∅  for all f ≠ g and all θ ∈ (0, Θ) } ,
```

with `Θ_max := 0` if the interiors already overlap for arbitrarily small `θ > 0`. This is the
quantity `collision.jl::theta_max` estimates by grid scan + bisection with faces shrunk by a
relative `1e−6`.

**[D] T4.1b (contact event).** An ordered pair `(w, (a,b))` — an `M′`-vertex `w` of face `g`, and a
consecutive pair `(a,b)` of `M′`-vertices of a face `f ≠ g` — is *in contact at `θ`* if
`y_w(θ)` lies on the **closed** segment `[y_a(θ), y_b(θ)]`. By T3 this is one harmonic equation and
two harmonic inequalities:

```
   E1:  h_o(θ) := det( Y_b − Y_a , Y_w − Y_a )  =  0                      (T3.2 coefficients)
   E2:  0  ≤  h_d(θ) := ⟨ Y_b − Y_a , Y_w − Y_a ⟩  ≤  ‖x_b − x_a‖²        (T3.3; the bound is a
                                                                            CONSTANT by T3.3′)
```

The right-hand bound is constant in `θ` because `a` and `b` belong to the same face, whose edge
lengths are `θ`-invariant (T1.D / T3.3). So the pair costs one quadratic solve (T3.6) and one
harmonic evaluation per root.

**[D] T4.1c.** `𝒞(X) :=` the set of angles `θ ∈ (0, π]` at which some pair is in contact — the
**candidate set**. By T3.4 each pair contributes at most 2 angles, so `|𝒞| ≤ 2 · (#vertex, edge
pairs) = O(n²)`, and every element is a closed-form arctangent of a quadratic root.

### T4.2 Completeness: the candidate list misses nothing

**Lemma T4.2 [A].** Let `P, Q` be closed simple polygons with `int P ∩ int Q = ∅` and `P ∩ Q ≠ ∅`.
Then some vertex of `P` lies on `∂Q` or some vertex of `Q` lies on `∂P`.

*Proof.* Take `z ∈ P ∩ Q ⊆ ∂P ∩ ∂Q` (the inclusion holds because a point of `int P` in `Q` would be
in `int Q` or on `∂Q`; if on `∂Q`, every neighbourhood meets `int Q`, and `int P` is open, so
`int P ∩ int Q ≠ ∅`, a contradiction; symmetrically for `int Q`). Suppose no vertex of either lies
on the other's boundary. Then `z` is in the relative interior of an edge `e_P` of `P` and of an edge
`e_Q` of `Q`. If `e_P` and `e_Q` cross transversally at `z`, the four local sectors force
`int P ∩ int Q ≠ ∅`; so `e_P ∥ e_Q` near `z`, and `e_P ∩ e_Q` contains a segment `s ∋ z`. Let `s` be
the maximal such segment. Each endpoint of `s` is an endpoint of `e_P` or of `e_Q`, i.e. a vertex of
`P` or of `Q`, and it lies on both edges, hence on the other polygon's boundary — contradiction. ∎

**Corollary T4.2′ (completeness) [A].** The pairwise interior-overlap status of the deployed
structure is locally constant on `(0, π] \ 𝒞(X)`. *Proof.* The map `θ ↦ (P_f(θ))_f` is continuous;
if the overlap status of some pair changed at `θ₀ ∉ 𝒞`, then by continuity `P_f(θ₀) ∩ P_g(θ₀) ≠ ∅`
while the interiors are disjoint on one side, so Lemma T4.2 applies at `θ₀` and `θ₀ ∈ 𝒞`. ∎
Note this uses T2: the motion is defined and analytic on all of `(0, π]`, so there is no other way
for the configuration to stop or change type.

**Theorem T4.2″ (exact `Θ_max`) [A].** Write `𝒞(X) = {θ₁ < θ₂ < … < θ_M}` and `θ_{M+1} := π`. Then

```
   Θ_max(X)  =  θ_{i*}   where  i* = min { i ≥ 0 : the interiors overlap at the midpoint of
                                            (θ_i, θ_{i+1}) } ,      θ₀ := 0 ,                (T4.1)
```

and `Θ_max = π` if no such `i` exists. Every `θ_i` is a closed-form arctangent; the only non-closed-
form ingredient is the *sign* selection, one overlap test per interval, and the scan stops at the
first hit. This is exact — not approximate — given the complete candidate list.

**This corrects the specification of `ideas/ranking.md` K2a and of `R2`.** Those define
`θ_max = min over pairs of the first admissible root`, i.e. `θ₁`. That is **false**: `θ₁` is the
first *contact*, which need not be an overlap. See the measured counterexample below.

**[N] Measured (`derivations/scratch/check_t4_t5.jl`, checks D1/D2).** 16 valid embeddings
(8 tilings × 2 seeds, plus Voronoi), `X = X₀`:

| pattern | `|𝒞|` | `θ₁` (first contact) | `Θ_max` by (T4.1) | `theta_max` bisection | `min_i β_i` |
|---|---|---|---|---|---|
| squares | 1 | 3.141593 | 3.141593 | 3.141593 | 3.141593 |
| triangles | 0 | — | 3.141593 | 3.141593 | 4.188790 |
| **hexagons** | 3 | **1.047198** | **2.094395** | 2.094397 | 2.094395 |
| kagome | 1 | 3.141593 | 3.141593 | 3.141593 | 3.141593 |
| snub_square | 7 | 1.670664 | 1.670664 | 1.670666 | 3.494324 |
| truncated_square | 10 | 2.356194 | 2.356195 | 2.356196 | 2.356194 |
| t3_4_3_12 | 2 | 2.617994 | 2.617994 | 2.617996 | 2.617994 |
| voronoi (2 cases) | 302 / 346 | 0.0766 / 0.00074 | **0** | 0 | 0.663 / −0.271 |

Worst `|Θ_max − theta_max_bisection| = 2.09e−6`, inside the `1e−5` tolerance that deviation 9 forces.
The **hexagon row is the counterexample**: `θ₁ = π/3` is a genuine contact where the two duplicates
of a split edge become collinear *end-to-end* (the vertex lands exactly at an edge **endpoint**,
measured contact parameter `t = 0.000` and `t = 1.000`), so the faces touch at a point and separate
again; the first interior overlap is at `2π/3`. A separate diagnostic confirms zero interior
overlaps at `θ ∈ {1.047, 1.05, 1.2, 1.5, 2.0, 2.09}` and 12 overlaps from `θ = 2.0944` on.

### T4.3 Genericity: first contact is vertex-into-edge-interior

**Proposition T4.3 [A, with a genericity hypothesis].** Fix the combinatorics `(M, σ)`. For `X` in a
dense open subset of `𝕏`, the binding contact at `Θ_max` is a vertex of one face in the **relative
interior** of an edge of another, and the orientation harmonic has a **simple** root there.

*Reasoning.* The competing degeneracies are (i) vertex-vertex coincidence, `h_o(θ*) = 0` together
with `h_d(θ*) ∈ {0, ‖x_b − x_a‖²}` — two independent harmonic equations in one unknown `θ`, hence
codimension 1 in `X`; (ii) a double root, `h_o(θ*) = h_o′(θ*) = 0`, equivalently the discriminant
`q² + r² − p² = 0`, again codimension 1 (and a quartic in `t`, T3.4); (iii) three-fold contact,
codimension ≥ 1. Each is the vanishing of a nonzero polynomial in `X` on `𝕏`, so its complement is
dense and open — **provided** the polynomial is not identically zero on `𝕏`. That proviso is exactly
where this fails in practice and I do not claim it in general: **CONJECTURE (T4.3-generic)**.

**Why it visibly fails on symmetric patterns.** On the hexagon pattern above, degeneracy (i) holds
*identically* on the whole shape space reachable by that `σ`, because the two duplicates of a split
edge are congruent translates (T1.B) — equal length by construction — so when they become collinear
the contact is automatically vertex-to-vertex. Symmetric tilings are precisely the measure-zero set,
and they are precisely the patterns both papers ship. **Any implementation that assumes T4.3 will be
wrong on the paper's own figures.** This is why T4.2″ is stated with the interval scan and not as
`min over roots`.

**Vertex-vertex handling.** At a vertex-vertex contact the correct local predicate is a *wedge* test:
the interiors overlap just after `θ*` iff the two face cones at the coincident point have
intersecting interiors for `θ` slightly greater. (T4.1) implements this by a single overlap probe;
an analytic version is one comparison of four harmonic angle functions and is not derived here —
**CONJECTURE (T4.3-wedge)**, low risk, moderate bookkeeping.

### T4.4 The analytic `β` bound of 2026 §5.2 as a special case

**[F, 2026 §5.2 / 2025 §4.2]** For a hinge edge `e` whose hinge vertex is `v`, with `α_f, α_g` the
interior angles of the two incident faces at `v`, the paper's local bound is
`β_e = 2π − α_f − α_g`, and `θ_max ≤ min_e β_e`.

**Derivation [A].** The hinge point `p_e(θ)` is common to both faces at every `θ` (T1.1). Around it
the two faces occupy angular sectors of fixed widths `α_f, α_g` (T1.D: faces are rigid), separated on
the cut side by a gap that opens to exactly `θ` (T1.C) and on the far side by
`2π − α_f − α_g − θ`. The far gap closes at `θ = β_e`, at which the two boundary edges of `f` and `g`
emanating from `p_e` become collinear and superposed. At that moment the far endpoint of the shorter
of the two edges lies on the closed other edge, so `β_e ∈ 𝒞(X)`: **the paper's bound is one of the
candidate roots, not a separate mechanism.** If the two edges have different lengths, the contact is
vertex-into-edge-interior and (by the wedge argument) the faces then cross, so `β_e` is a genuine
overlap angle; if equal, it is vertex-vertex and the wedge test decides.

**Corrected statement (D3).** Round 1 printed "`Θ_max` equals `min_e β_e`" on split-free patterns,
which contradicts this file's own T4.2 table (triangles: `Θ_max = 3.141593`, `min β = 4.188790`).
The `β` bound is a bound on the *hinge sector*, which can exceed the search range. Since deviation 9
fixes the search range at `(0, π]`, the statement is

```
        E_split = ∅        ⟹        Θ_max  =  min( min_e β_e ,  π ) .                  (T4.1b)
```

**Which half of (T4.1b) is proved (round 3, `check.md` R2.4a).** The two directions have very
different status and round 2 printed them together, untagged, as "the true statement". They are:

* **`≤` — [A], derived.** The paragraph above shows `β_e ∈ 𝒞(X)` and, when the two edge lengths at
  the hinge differ, that the faces actually cross there; with deviation 9's cap this gives
  `Θ_max ≤ min(min_e β_e, π)`. This half is a proof.
* **`≥` — [N], MEASURED, not derived.** The reverse inequality says no *other* candidate root binds
  earlier — that no vertex of any face runs into any non-adjacent face below `min(min β, π)`. Nothing
  in T4.2–T4.5 forbids that; the candidate list is over all ordered (vertex, edge) pairs, not only
  hinge-local ones, and a distant pair could in principle have an admissible crossing root first. I
  have **no proof** of the `≥` direction and I do not assert one. What I have is the measurement:
  equality holds to **8.88e−16** on all **8** split-free patterns of my corpus (squares, triangles,
  kagome, t3_4_3_12, two seeds each) **[N, R2-C]**, independently reproduced by the Checker on **4**
  split-free patterns to **8.9e−16** (test **T4-b**). That is 8 patterns, not a theorem.

So (T4.1b) must be printed as `Θ_max ≤ min(min_e β_e, π)` **[A]**, with equality
`[N]`-supported on split-free patterns and **not proved**. Anywhere the paper wants the equality it
must carry the `[N]` tag and the sample size, or someone must prove the `≥` direction — a plausible
route is to show that on a split-free pattern every candidate pair whose root is admissible below
`min β` is itself hinge-local, which I could not do.

**[N] Measured (`derivations/scratch/check_r2.jl`, R2-C).** On all **8 split-free** patterns
tested (squares, triangles, kagome, t3_4_3_12, ×2 seeds), `Θ_max` computed by (T4.1) equals
`min(min_e β_e, π)` to **8.88e−16**; triangles is precisely the case that separates the two forms.
Independently reproduced by the Checker (**T4-b**, 8.9e−16), whose own hand recomputation of `β_e`
(CCW interior angle, reflex corners included) matches `hinge_beta()` **exactly**, 0.000e+00 over 933
hinge edges (**T4-a**). This reproduces `[F, F18]` from the candidate list rather than from
simulation, and shows the local bound is *tight* exactly when there are no split cuts (up to the
`π` cap). With split cuts it is not: snub square `Θ_max = 1.671` vs `min β = 3.494`; consistent
with `[F, F18]`.

**Caveat for the paper.** `Θ_max = min(min β, π)` is a statement about the range *this project
searches*, not a theorem that the structure cannot deploy past `π`. Nobody in this project has
looked at `θ > π` (T4.H.5), so the `π` in (T4.1b) is a convention, and it must be printed as one.

### T4.5 The broad phase (round 2: this is a pruning section, not a complexity theorem)

**[A] T4.5a (the exact broad phase).** Let `x̄_f` be the centroid of face `f` in `M` and `r_f` its
circumradius about `x̄_f`. By T1.D the deployed face lies in the disc of radius `r_f` about
`γ_f(θ) := R_f x̄_f + t_f = c x̄_f + s J(2u_f − σ_f x̄_f)`. Hence `P_f(θ) ∩ P_g(θ) ≠ ∅` requires

```
        ‖ γ_f(θ) − γ_g(θ) ‖²  ≤  (r_f + r_g)² .                                        (T4.2)
```

`γ_f − γ_g` has the form `c A + s B` with `A, B` linear in `X`, so the left side is **itself a
harmonic** (T3.4), and (T4.2) holds on a closed-form union of at most two intervals of `θ`. Two uses:

* **Exact `θ`-dependent pruning**: pair `(f,g)` may be dropped outside those intervals, with no loss.
* **Static pruning** (no `θ`): with `h(θ) = p + q cos θ + r sin θ` the harmonic of the left side,
  `min_θ h = p − √(q² + r²)`. If

```
        p − √( q² + r² )  >  ( r_f + r_g )²                                            (T4.3)
```

  the pair can be discarded for the whole deployment. All quantities are quadratic in `X`.

**Both tests use the MOVING centroid `γ_f(θ)`, and that is not optional (D8).** A variant that
compares the **flat** centroids, `‖x̄_f − x̄_g‖ ≤ r_f + r_g`, is **unsound in principle**: by T1.D the
face radii `r_f` are `θ`-invariant while deployment *contracts* centroid distances, so a pair whose
faces actually approach can be discarded. `Kirigami/src/method/contact.jl::candidate_pairs(...,
use_static = true)` is that variant. Measured by the Checker: it discards **19 842** pairs the sound
moving test keeps, with `max(‖x̄_f − x̄_g‖ − min_θ‖γ_f(θ) − γ_g(θ)‖) = 6.59`; it nevertheless returned
the correct `Θ_max` on 96/96 samples (**T4-g**), so this is a soundness defect with no observed
counterexample yet. The sound moving form keeps 42.6 % of pairs and changed `Θ_max` on 0 of 96
(**T4-f**). **Use (T4.2)/(T4.3); do not use the flat-centroid form, and do not report a number
produced by it.**

**[A] T4.5b (the swept radius — round-1 claim WITHDRAWN, D1).** For a single copy,
`y(θ) = c C + s S` traces the ellipse of T1.A, so `max_θ ‖y(θ)‖ = σ_max([C | S])`, the **larger
singular value**, given by (T1.7), and in a *general* frame

```
        max( ‖C‖, ‖S‖ )  ≤  σ_max( [C|S] )  ≤  √( ‖C‖² + ‖S‖² )  ≤  √2 · max( ‖C‖, ‖S‖ ) .   (T4.4)
```

Round 1 concluded from (T4.4) that `ideas/ranking.md` K2c's `ρ_f = max(‖x‖, ‖χ‖)` "is not an upper
bound", that its pruning is unsound, and that `ρ_f` must be multiplied by `√2`. **That conclusion
is wrong and is withdrawn.** It measured (T4.4) in the *raw* frame — origin at the BFS seed face —
which no pruning uses. K2c's spec and `contact.jl::swept_discs` both work **in the face's own
frame**, `x = C_u − gc_f`, `χ = S_u − gs_f`, and there the inequality is an *equality*:

**[A] Lemma T4.5b′ (the swept radius in the face frame is exact).** With `γ_f(θ) = c·gc_f + s·gs_f`
the moving centroid, (T1.5) gives `gc_f = x̄_f` and `gs_f = J(2u_f − σ_f x̄_f)`, so

```
        C_u − gc_f  =  x_u − x̄_f ,
        S_u − gs_f  =  J( −σ_f x_u + σ_f x̄_f )  =  − σ_f · J ( x_u − x̄_f ) .            (T4.4′)
```

The `2u_f` term **cancels identically**. The two columns are therefore orthogonal and of equal
norm, so `[x | χ]` is a similarity and

```
        max( ‖x‖, ‖χ‖ )  =  σ_max( [ x | χ ] )  =  max_θ ‖ y_u(θ) − γ_f(θ) ‖  =  ‖ x_u − x̄_f ‖ .
```

The swept radius about the moving centroid is **exactly the flat circumradius** — which is just
T1.D (faces are rigid) written in coordinates. ∎

**[N] Measured (`derivations/scratch/check_r2.jl`, R2-A).** Over **1 628** `M′`-copies on 16
graphs, sampling `θ` on `[0, π]`: `| max_θ‖y_u − γ_f‖ − max(‖x‖,‖χ‖) | ≤ 3.11e−15` and
`| max_θ‖y_u − γ_f‖ − ‖x_u − x̄_f‖ | ≤ 3.11e−15`. The Checker measured the same thing independently
(**T4-c/T4-d**, 1.4e−14 / 1.9e−14 over 15 132 copies). Round 1's 68 %-violation figure is real but
is about the raw frame, and it has no consequence for any pruning anyone runs.

**Two consequences round 1 had backwards.**

1. Of the two radii in `contact.jl`, `rho` (the `√(‖x‖²+‖χ‖²)` form, commented "sound bound") is
   the one loose by exactly `√2`; `rho_max` (the spec's form) is **tight**. Applying the prescribed
   `√2` factor would have made the pruning looser, not sound.
2. **K2c's gate test is vacuous.** By (T4.4′), `ρ_f / r_f ≡ 1` for every face of every patch at
   every size, so the log–log fit of `log(max_f ρ_f / r_f)` against `log n` has slope 0
   identically. **[N, R2-B]** measured `1.000000` at `F = 16, 36, 64, 100, 144, 196` (Checker:
   `1.000000` at `F = 12, 24, 44, 68, 96, 156`). **K2c as specified cannot fail and therefore
   cannot test anything.** It must be redesigned or dropped. **This is not hypothetical:**
   `results/kill/KILL_REPORT.md` line 57 records `K2c | **PASS** | log-log slope 0.0000`, and line 34
   already notes the measurement `max_f ρ_f / r_f = 1.000000000001`. That PASS is an artefact of
   (T4.4′), not evidence for locality, and the K2c row of the kill report must be **withdrawn**.

**[REFUTED] H-LOC (D2).** Round 1 stated the hypothesis as: there is a patch-independent `κ` with
`σ_max([C_u | S_u]) ≤ κ · r_{f(u)}` in the face-centroid frame — "equivalently, the path sum `u_f`
of (T1.6) does not grow with the diameter of `Γ`". **The two halves are not equivalent, and the
error is mine.** By Lemma T4.5b′ the first half is *trivially true with `κ = 1`*: `u_f` has
cancelled out of it, so it carries no information about locality at all. The `O(n)` packing
argument needs the second half — that every face's swept region sits in a disc of radius `O(r_f)`
about a **fixed** point, its flat centroid — i.e. `‖γ_f(θ) − x̄_f‖ = O(r_f)` uniformly in the patch.

That statement is **false**. Deployment contracts the patch by roughly `cos(θ/2)` about the seed, so
a face at distance `d` from the seed moves by `Θ(d)`. **[N, R2-B]**, growing square patch, fixed
boundary, `σ` from Eq. (1):

| clip radius | `F` | patch diameter | `max_f max_θ ‖γ_f(θ) − x̄_f‖ / r_f` | `max_f ρ_f / r_f` |
|---|---|---|---|---|
| 1.6 | 16 | 5.66 | 4.62 | 1.000000 |
| 2.6 | 36 | 8.49 | 7.71 | 1.000000 |
| 3.6 | 64 | 11.31 | 10.79 | 1.000000 |
| 4.6 | 100 | 14.14 | 13.87 | 1.000000 |
| 5.6 | 144 | 16.97 | 16.95 | 1.000000 |
| 7.0 | 196 | 19.80 | 20.03 | 1.000000 |

The drift is linear in the patch diameter (ratio ≈ 1.0 against diameter over a 3.5× range), exactly
as the contraction picture predicts. The Checker's independent table (different generator, 12–156
faces) shows the same slope. **H-LOC in the form the `O(n)` count needs is REFUTED, not merely
unproved**, and it moves out of the CONJECTURE list into the dead-ends.

**What this costs, and what survives.** T4.5 ships as the **exact broad phase** (T4.2)/(T4.3) only:
sound with no hypothesis, verified to preserve `Θ_max` on 96/96 samples, keeping 42.6 % of pairs.
The candidate list is `O(n²)` and nothing in this file may be labelled `O(n)`, "certified active
set" or "locality theorem". **What might still rescue an `O(n)` count** is a *deployment-frame*
statement rather than a flat-frame one: the motion is close to a global similarity, so relative
centroid distances contract uniformly and the packing count could be done in the deployed frame —
i.e. bound `‖γ_f(θ) − γ_g(θ)‖` relative to `‖x̄_f − x̄_g‖` rather than bounding drift from `x̄_f`. The
data above is consistent with that version. **It is not derived here, and until it is, no
complexity claim exists.** (Recorded as an open direction, not as a conjecture the paper leans on.)

### T4.H What can break

1. **Grazing contacts** (hexagons above) — handled by T4.2″, fatal to the `min-over-roots` formula.
2. **Degenerate harmonics** `h_o ≡ 0`: permanent incidences at hinge points; must be removed by an
   identity test, not a root test (T3.H.1).
3. **`Θ_max = 0`**: `X` positively oriented but not embedded. Measured on 2 of 2 Voronoi cases that
   passed the positive-orientation screen — see T5.1 and T6.4.
4. **Non-simple faces** in `M`. Everything above assumes each face is a simple polygon — since
   round 5 this is printed as a standing hypothesis of T4–T6, not left implicit. Measured cost of
   dropping it: 5 of 873 local sector tests fail, all on self-intersecting faces (Checker R4-a3).
5. **`θ > π`.** The code searches `(0, π]`; `𝒞` is defined on the same range. Contacts beyond `π`
   are not considered by anyone in this project.

### T4.Check

* Re-run `derivations/scratch/check_t4_t5.jl` (rules: D1 `≤ 1e−5`, D2 `≤ 1e−12`) **and**
  `derivations/scratch/check_r2.jl` (rules: R2-A `≤ 1e−12`, R2-C `≤ 1e−12`, R2-D violations `= 0`;
  its R2-E line is **withdrawn**, see T5.2b.2) **and** `derivations/scratch/check_r3.jl`
  (rule: three-class deflated mismatches `= 0`; two-class mismatches are expected to be nonzero and
  `ε`-independent). Round 1's check D4 stands only as a statement about the
  raw frame; it is **not** evidence about any pruning — see T4.5b.
* **Attack the hexagon row**: it is the single most falsifiable claim in T4. Verify independently
  that at `θ ∈ (π/3, 2π/3)` no two faces of the hexagon pattern have overlapping interiors.
* Re-derive `β_e` by hand on the 4.8.8 pattern (`min β = 2.356194 = 3π/4`; the octagon interior angle
  is `3π/4` and the square's is `π/2`, so `2π − 3π/4 − π/2 = 3π/4` ✓ — that arithmetic is the check).
* Run K2a with (T4.1), not with `min over roots`; the latter will fail on hexagons and on any
  pattern with congruent split-edge duplicates.
* **K2c must be redesigned before it is run** (D1): its gate quantity is identically `1` (T4.5b),
  so a PASS carries no information. A replacement that would test something real: measure
  `|{pairs surviving (T4.3)}|` against `n` on the largest available patches, and separately measure
  `max_f max_θ ‖γ_f(θ) − x̄_f‖ / r_f` against patch diameter — the second is the quantity that
  actually decided H-LOC.

---

## T5 — The usable region

**[D] T5.0.** `U(ε) := { X ∈ 𝕏 : every face of M positively oriented at θ = 0, and Θ_max(X) ≥ ε }`,
parametrized by `t ∈ R^{k×2}` via `X(t) = X₀ + Φ t`.

### T5.1 Positive orientation is a finite set of quadratic inequalities

**[A]** The signed area of face `f` with vertices `v₀ … v_{m−1}` is
`A_f(X) = ½ Σ_i det(x_{v_i}, x_{v_{i+1}})`, a quadratic form in `X`, hence (composing with the affine
`X(t)`) a quadratic polynomial in `t`. Positive orientation of all faces is the finite conjunction

```
        A_f( X₀ + Φ t )  >  0 ,      f ∈ F ,                                          (T5.1)
```

`|F|` strict quadratic inequalities — a **basic open semialgebraic** set. By T1.D it is equivalent to
positive orientation at every `θ`, so it needs checking only once.

**[A] What (T5.1) does NOT give.** Positive orientation of every face does **not** imply that the
flat structure `M` is embedded: the piecewise-linear map can be locally injective and globally
overlapping. **[N]** Measured: 2 of the 16 accepted embeddings (both Voronoi, `X = X₀`) satisfy
(T5.1) and nonetheless have `Θ_max = 0` — interiors overlap immediately. So `n_inv = 0`, the proxy
used in `[F, F17]` and in kill experiment K1a, is a **necessary but not sufficient** validity
certificate, and K1a's `p_valid` therefore *over*-counts valid samples. Independently reproduced by
the Checker on a larger sample: **51 of 70** positively-oriented shape-space samples have
`Θ_max = 0`, with a witness saved at
`derivations/check_failures/T5_1_positive_orientation_zero_range_hexagons.json`.

**[D] Two different things are being certified (D11).** Round 1 named `Θ_max(X) > 0` as "the
sufficient certificate", collapsing two conditions that are not the same:

* **Flat-pattern validity (what you fabricate).** `min_f A_f(X) > 0` **and** the flat structure
  `M(X)` — *the uncut mesh, not `M′`* — is an embedded planar graph, i.e. no two faces of `M` have
  overlapping interiors. Positive orientation alone gives only *local* injectivity, which is exactly
  what the 51 counterexamples exploit. **Round 3 correction:** this predicate is about `M`, and it
  must **not** be confused with, or used in place of, the deployment hypothesis of T5.2b′. Applying
  an interior-overlap test to the cut structure `M′` at `θ = 0` is degenerate — the two copies of
  every split edge coincide there (T1.C) — so the deployment side is certified at `θ₁ = ε/2`
  instead; see T5.2b.0 and T6.4 item 2.
* **Deployment validity (what you deploy).** `Θ_max(X) > 0` by T4.2″.

Neither implies the other in general: a flat state with two faces overlapping could in principle
separate as the cuts open, and an embedded flat state can have `Θ_max = 0`. The certificate of T6.4
therefore carries more than one validity predicate: `min_f A_f > 0`, embeddedness of the flat mesh
`M`, and the deployment side `Θ_max > 0` / `≥ ε`. What the *deployment* side needs is not a
`θ = 0` test at all but the pointwise `NOOVERLAP(ε/2)` of T5.2b.0 — that is the atom missing from
round 1's T5.2b, and round 2's `θ = 0` phrasing did not supply it (`check.md` R2.4b).

### T5.2 `Θ_max ≥ ε` is semialgebraic — with an explicit inner region, not basic, not by quadrics

Fix `T := tan(ε/2)`. For a candidate pair `π = (w, (a,b))` write, from (T3.2)/(T3.5),

```
        g_π(τ; t)  =  ( p_π + q_π )  +  2 r_π τ  +  ( p_π − q_π ) τ² ,
```

with `p_π, q_π, r_π` quadratic polynomials in `t` (T3.2 + T1 Step 8), and the two interval harmonics
`h_d,π` likewise. Then

```
   Θ_max(X(t)) ≥ ε
        ⟺  ∀ π :  ¬ ∃ τ ∈ (0, T) [ g_π(τ; t) = 0  ∧  0 ≤ h_{d,π}(τ; t) ≤ ℓ_π(t)  ∧  (overlap at τ⁺) ] .
```

**[A] T5.2a (semialgebraicity).** The bracketed formula is a first-order formula over the reals in the
free variables `t` and the bound variable `τ`, with polynomial atoms. By Tarski–Seidenberg, its
`τ`-projection is semialgebraic; the conjunction over the finite candidate list and with (T5.1) is
semialgebraic. Hence **`U(ε)` is a semialgebraic subset of `R^{2k}`** — this is the honest statement.
(Prior art: configuration spaces of linkages are semialgebraic — Kapovich–Millson, King — is
textbook, per `notes/screen_r1.md`; the content here is the *explicit* description, not the fact.)

### T5.2b An explicit conservative description — and a proof that it is *inner*

Round 1 stated this as "drop the interval and overlap tests, ask only that no orientation harmonic
has a root in `(0,T)`; the region shrinks, so it is conservative". **Shrinking is not the same as
being contained in `U(ε)` (D4)**, and round 1 never proved containment. It also printed an atom
combination that is degree 8, not 4 (D5). Both are fixed here, and the containment is now proved.

**[D] T5.2b.0 The three ingredients (revised in round 3).** Fix `T := tan(ε/2)` (`ε < π`, so `T` is
finite and `θ ↦ tan(θ/2)` is an increasing bijection `(0,ε) → (0,T)`). Fix one explicit interior
angle

```
        θ₁ := ε/2  ∈ (0, ε) .
```

For `X = X(t)` define

```
  (i)   POS(t)          :  A_f(X(t)) > 0  for every f ∈ F                             (T5.1)
  (ii)  NOOVERLAP(θ₁)   :  the deployed faces { P_f(θ₁) } have pairwise disjoint
                           interiors — one exact polygon–polygon test per face pair,
                           at the single angle θ₁ = ε/2
  (iii) NOROOT(t)       :  for every candidate pair π ∈ 𝒞-list, the orientation harmonic
                           h_o,π has no root in (0, ε) — decided by the deflated atom
                           list (T5.1c)/(T5.1d)/(T5.1e′) of T5.2b.2
```

**[D] What the `𝒞-list` in (iii) is, exactly (round 4, `check.md` R3.1a).** The `𝒞-list` is the
complete ordered candidate list of T4.1b / Corollary T4.2′ — all ordered pairs `(w, (a,b))` with `w`
a vertex of one face and `(a,b)` a consecutive vertex pair of a *different* face — **with the pairs
whose orientation harmonic is identically zero removed**:

```
        𝒞-list  :=  { π = (w,(a,b)) :  f ≠ g }  \  { π :  h_o,π ≡ 0 } .
```

The removal is by the **identity** test of T3.H.1 (`|p| + |q| + |r| ≤ tol·scale`, a test on the
coefficients of `h_o,π`, evaluated once per pair), never by a root test, and it happens **before**
the three-class split of T5.2b.2, so an identically-zero pair is never classified and never
contributes an atom — in particular it is not a class-3 pair of (T5.1e′), whose atom now carries the
side-condition `A ≠ 0`. Without this removal (iii) is unsatisfiable as literally written: every
hinge point of the pattern contributes a permanently incident pair (T3.H.2) whose harmonic vanishes
at every `θ`, so `R(ε)` would be empty for every pattern with a hinge. The Checker measured **14 928**
such pairs on its corpus and both check programs already drop them first; this paragraph states the
exclusion that the programs implement. Proposition T5.2b′ below carries the burden of showing that
the removed pairs cannot hide an overlap transition.

and put `R(ε) := { t : POS ∧ NOOVERLAP(θ₁) ∧ NOROOT }`. `(iii)` is stricter than the contact
condition — it drops the two interval tests `E2` and the overlap selection — so `R(ε)` is smaller
than the exact region.

**[D] The predicate round 2 used, and why it is replaced.** Round 2's second ingredient was

```
        EMB(t)  :  the deployed faces have pairwise disjoint interiors for all
                   sufficiently small θ > 0                                    ("θ = 0⁺ atom")
```

the negation of `contact.jl::penetrates_immediately()`. `EMB` is what the proof needs, but it is a
**germ at `0⁺`**, not a test at a point: it quantifies over an unspecified neighbourhood, so a
certificate cannot report it as a measured fact. Meanwhile T5.1/T6.4 printed a *different*
predicate — "no face–face interior overlap **at `θ = 0`**" — and the Checker is right
(`check.md` R2.4b) that the two are not the same and that flat embeddedness does not supply `EMB`.

**The `θ = 0` overlap test is not merely a different predicate; it is degenerate and cannot be
used.** At `θ = 0` the cut structure `M′` is the *flat* state: by (T1.C) the gap across every cut
opens to exactly `θ`, so at `θ = 0` the two `M′`-copies of every split edge **coincide as segments**
and every pair of faces adjacent across a cut **touches along a whole shared edge**. An exact
polygon–polygon predicate at `θ = 0` therefore runs entirely in its own degenerate case: the two
polygons share a one-dimensional set of boundary points, so "interiors disjoint" is decided by
arbitrarily small perturbations of the shared edge and by whatever tolerance the test uses. That is
why the certificate must be evaluated at an angle *strictly inside* the interval. At any
`θ ∈ (0, π)` the only remaining coincidences are the hinge points, which are single points shared by
construction (T1.1) and are boundary, never interior, contacts — so the test at `θ₁ = ε/2` is well
posed, and it is a finite exact computation (segment–segment intersection plus a point-in-polygon
test per pair, over the `O(n²)` face pairs).

**[A] Semialgebraicity is unaffected, the degree is still not controlled.** `NOOVERLAP(θ₁)` is a
finite Boolean combination of polynomial sign conditions in the vertex coordinates
`y_v(θ₁) = cos(θ₁/2)C_v + sin(θ₁/2)S_v`, which are *affine* in `t` for the fixed angle `θ₁` — so
unlike round 2's `EMB` (which needed the semialgebraic function `θ₁(t) = min 𝒞(X)`) this one
involves no auxiliary quantifier at all. `U(ε)` remains semialgebraic and T5.2a is unaffected. But I
still have **no useful degree bound**: polygon–polygon interior disjointness for `m`-gons is a
Boolean combination of `O(m²)` orientation determinants, each quadratic in the vertex coordinates
hence quadratic in `t`, and the Boolean structure — not the degree of a single atom — is what is
uncontrolled. Any claim of the form "quantifier-free description of degree ≤ 4" must therefore be
stated about `POS ∧ NOROOT`, with `NOOVERLAP(θ₁)` named separately. This is a real limitation of the
description, and it is the honest price of D4.

**[A] Proposition T5.2b′ (inner approximation, restated in round 3).** `R(ε) ⊆ U(ε)`. That is,

```
        POS(t) ∧ NOOVERLAP(ε/2) ∧ NOROOT(t)      ⟹      Θ_max(X(t)) ≥ ε .
```

Moreover the same hypotheses imply `EMB(t)`.

*Proof.* Let `O := { θ ∈ (0, ε) : ∃ f ≠ g with int P_f(θ) ∩ int P_g(θ) ≠ ∅ }`. The interval `(0, ε)`
is connected, so it suffices to show `O` is **open and closed in `(0, ε)`** and that `O` misses one
point; then `O = ∅`. In words: *the overlap status can only change at a candidate root, and `NOROOT`
says there is none inside, so the status is constant on `(0, ε)` and `NOOVERLAP(θ₁)` reads it off.*

`O` is **open**. Suppose `z ∈ int P_f(θ) ∩ int P_g(θ)` and let `ρ > 0` be the smaller of `z`'s
distances to `∂P_f(θ)` and `∂P_g(θ)`. Every vertex position `y_v(θ) = cos(θ/2)C_v + sin(θ/2)S_v` is
continuous in `θ` (T1), and each `∂P` is a finite union of segments with those endpoints, so for
`θ′` near `θ` every boundary point moves less than `ρ/2`. Hence `z` stays in both interiors and
`θ′ ∈ O`.

`O` is **closed in `(0, ε)`**. Let `θ_k ∈ O` with `θ_k → θ* ∈ (0, ε)`; we show `θ* ∈ O`. Suppose not,
so all face interiors are pairwise disjoint at `θ*`. There are finitely many face pairs, so after
passing to a subsequence one fixed pair `(f, g)` overlaps at every `θ_k`; choose
`z_k ∈ int P_f(θ_k) ∩ int P_g(θ_k)`. All the polygons stay in a fixed compact set (their vertices
lie on the bounded ellipses of T1.A), so after a further subsequence `z_k → z`, and by continuity
`z ∈ P_f(θ*) ∩ P_g(θ*)` — the closed polygons. So at `θ*` we have `int P_f ∩ int P_g = ∅` and
`P_f ∩ P_g ≠ ∅`. **Lemma T4.2** applies verbatim: some vertex of one polygon lies on the boundary of
the other. Moreover its proof is **local at any prescribed contact point**: given `z ∈ P_f ∩ P_g`
(necessarily `z ∈ ∂P_f ∩ ∂P_g`, since the interiors are disjoint), either `z` is itself a vertex of
one polygon lying on the boundary of the other, or the two edges through `z` are transversal — which
forces four local sectors and hence an interior intersection, excluded — or they are parallel, in
which case the maximal shared segment `S ∋ z` of `∂P_f ∩ ∂P_g` has two endpoints, and **each**
endpoint is a vertex of one polygon lying on the boundary of the other (the shared segment stops
exactly where an edge of one of the two polygons ends). Call any vertex/edge pair produced this way a
**witness at `z`**. A witness is a candidate pair `π = (w,(a,b))` of the complete ordered list of
T4.1b / Corollary T4.2′, and it satisfies `h_o,π(θ*) = 0` with `θ* ∈ (0, ε)`.

If some witness at `z` has `h_o,π ≢ 0`, then `π ∈ 𝒞-list` (T5.2b.0(iii) removes only the identically
zero harmonics), and `h_o,π(θ*) = 0` contradicts `NOROOT`. So it remains to handle the case in which
**every** witness at `z` is identically zero, i.e. was struck from `𝒞` by T3.H.1. This is not an idle
case: hinge-adjacent faces carry a permanent incidence at every `θ`, and the Checker counts 14 928
such pairs on its corpus. It is disposed of as follows.

**[A] Sub-lemma T5.2b″ (a permanent incidence is never the only witness at an overlap transition).**
Assume, as T3.H.1/T3.H.2 classify them, that a candidate pair with `h_o,π ≡ 0` is a *permanent
coincidence*: `w` and `a`, or `w` and `b`, are two `M′`-copies of the same vertex of `M`, so the
witness point is a **hinge point shared by `f` and `g`** (the split-edge coincidences of T3.H.2 are
incidences at `θ = 0` only, and `θ* > 0`). Suppose every witness at `z` is of this kind. Then the
situation described above — interiors disjoint at `θ*`, interiors overlapping at `θ_k → θ*` with
`z_k → z` — is impossible unless some **other** candidate pair, with `h ≢ 0`, has a root at `θ*`.

**[A0] Sub-lemma T5.2b‴ (two faces sharing a hinge point are hinge-adjacent, so "the" `β_e` of the
pair is well defined).** Fix `v ∈ V` and a face `f` incident to `v`. By 0.2 a hinge edge `e` keeps
the hinge only at `src(e)`, and the direction `src(e) → dst(e)` is the common `σ`-induced half-edge
direction of the two faces incident to `e`. In the `σ_f`-traversal of `∂f` the two edges of `f` at
`v` are consecutive, one *entering* `v` and one *leaving* it — **[R5 correction applied R6]** this
step is the **simplicity hypothesis** of T4–T6 in disguise: it says the traversal of `∂f` visits `v`
**once**, which is exactly what a self-intersecting face polygon can violate (Checker R5.3); a hinge edge `e` at `v` keeps `v` for
`f` exactly when `v = src(e)`, i.e. exactly when `e` is the edge *leaving* `v` in that traversal.
Hence `f` keeps `v` for **at most one** of its incident hinge edges. So the identification of the
copies of `v` in `M′` is a *matching* on the faces around `v` — never a fan of three or more faces
sharing one `M′`-vertex — and two faces that share a hinge point are joined by **exactly one** hinge
edge `e`, whose far-side gap `β_e = 2π − α_f − α_g` is therefore unambiguous. This step is used
silently by Case A below (it is what makes "the" `β_e` of the pair `(f, g)` exist) and it makes
Case B vacuous; measured consequence: **0 of 933** face pairs share two or more `M′`-vertices
(Checker R4-b). ∎

*Proof.* Let `H` be the (finite) set of hinge points shared by `f` and `g`. Because faces are rigid
(T1.D), there is a radius `r > 0`, **independent of `θ`**, such that for every `p ∈ H` and every `θ`,
`P_f(θ) ∩ B(p, r)` is the intersection with `B(p, r)` of the closed angular sector of opening `α_f`
that `f` occupies at `p`, and likewise for `g` with `α_g`. (`r` may be taken to be half the distance
from `p` to the union of the edges of `f` and of `g` not incident to `p`, a quantity fixed by the two
rigid face shapes.) By T1.1/T1.C the two sectors at `p` are separated by the cut-side gap `θ` and the
far-side gap `2π − α_f − α_g − θ = β_e − θ`, so

```
        int P_f(θ) ∩ int P_g(θ) ∩ B(p, r) ≠ ∅        ⟺        θ > β_e .            (T5.2b″-1)
```

**Where simplicity is used.** The sector description of `P_f(θ) ∩ B(p, r)` needs each face polygon
to be **simple**: for a self-intersecting face the intersection with a small ball at `p` is not the
single closed angular sector of opening `α_f`, and (T5.2b″-1) fails. Simplicity enters twice and
independently here — Lemma T4.2, on which this whole argument rests, is also stated for closed
*simple* polygons — and it is now printed among the standing hypotheses of T4–T6. This is not
pedantry: with positive orientation alone the Checker's first local run produced **5 failures of the
`θ > β_e ⟹ overlap` direction out of 873**, and all 5 were hinges at which a face polygon was
self-intersecting (`voronoi2`, `voronoi3`, reflex `α ≈ 5.4–6.2`, `β_e < 0`); excluding the 99
non-simple hinges removed every failure (Checker R4-a3). `β_e` is well defined for the pair by
Sub-lemma T5.2b‴.

*Case A: `z ∈ H`, say `z = p`.* Then `|z_k − p| < r` for large `k`, so `θ_k > β_e` by (T5.2b″-1);
and the interiors are disjoint at `θ*`, so `θ* ≤ β_e` by the same equivalence. Letting `k → ∞` gives
`θ* ≥ β_e`, hence `θ* = β_e`. At `θ = β_e` the far-side gap closes: by T4.4 the two boundary edges of
`f` and `g` emanating from `p` on the far side become collinear and superposed, so the far endpoint
`w` of the shorter of the two lies on the closed other edge `(a,b)`. The pair `π′ = (w,(a,b))` is in
the complete ordered list and `h_o,π′(θ*) = 0`. It is **not** identically zero. Explicitly, take
`a = p`, let `L := ‖b − a‖` and `L′ := ‖w − a‖` be the two far-side edge lengths (`θ`-invariant by
T1.D) and note that the signed angle from `b − a` to `w − a` is `±(β_e − θ)`, the far-side gap; so

```
        h_o,π′(θ)  =  det( b − a , w − a )  =  ± L L′ · sin(β_e − θ) .            (T5.2b″-1b)
```

Since `L, L′ > 0`, `h_o,π′ ≢ 0`, which together with `h_o,π′(β_e) = 0` is all the proof needs.
**[R5 correction applied R6]** Its zero set is `{β_e + kπ : k ∈ Z}`, and `sin(β_e − θ)` has **at most
one** zero in the open interval `(0, π)`, which has length `π` — exactly one unless `β_e ≡ 0 (mod π)`.
In the situation at hand `θ* = β_e ∈ (0, ε) ⊂ (0, π)`, so the zero in range *is* `β_e`, and
`π′ ∈ 𝒞-list` has a root in `(0, ε)`, contradicting `NOROOT`. (The round-5 text printed "**exactly
one** of `{β_e, β_e + π, β_e − π}` lies in `(0, π)`" with a three-way case split; that is **false**
whenever `β_e ≤ −π`, where none of the three lies in `(0, π)` and the root in range is `β_e + 2π` —
measured on 5 008 of 20 000 synthetic `β_e ∈ (−2π, 2π)`, Checker R5.4. That regime needs
`α_f + α_g ≥ 3π`, both faces strongly reflex at the shared hinge, and is excluded here by
`β_e = θ* ∈ (0, ε)`; no theorem was affected, but the clause is now deleted.)

*(Round-5 correction of an overstatement.* Earlier rounds printed "`h_o,π′` vanishes at that one
angle alone", which is false: `β_e ± π` is a zero too. Only `β_e` is a **transition of the
sector-overlap predicate** (T5.2b″-1): at `θ = β_e ± π` the far-side gap is `β_e − θ = ∓π`, i.e. the
two far-side edges are collinear but point in *opposite* directions — the far gap is not closing, and
the predicate `θ > β_e` does not change value there. So the extra zero is a candidate angle that is
not a contact transition; it can only make `NOROOT` reject more designs, hence make `R(ε)` smaller
and still inner, and nothing downstream changes.) (This is precisely the `β_e` event of T4.4 — the
overlap of two hinge-adjacent faces at their shared hinge can only begin through the adjacent-edge
collinearity, which is a *different*, non-degenerate candidate pair.)

*Case B: `z ∉ H`.* By hypothesis every witness at `z` is a permanent coincidence, hence a point of
`H`, hence `≠ z`. By the localization above this can only be the parallel case: `z` lies in the
relative interior of a maximal shared segment `S ⊆ ∂P_f(θ*) ∩ ∂P_g(θ*)`, whose two endpoints
`u₁ ≠ u₂` are witnesses. So `u₁, u₂ ∈ H`: `f` and `g` share **two** distinct hinge points. Now `u₁`
and `u₂` are material points of *both* rigid faces, with fixed body coordinates in each. Writing
`P_f(θ) = γ_f(θ)(P_f⁰)` and `P_g(θ) = γ_g(θ)(P_g⁰)` with `γ_f, γ_g` the orientation-preserving
isometries of T1 / (T4.4′), the relative map `R(θ) := γ_g(θ)^{-1} ∘ γ_f(θ)` is an
orientation-preserving isometry that carries the two fixed body points `u₁, u₂` of `f` onto the two
fixed body points `u₁, u₂` of `g`. An orientation-preserving plane isometry is determined by the
images of two distinct points, so `R(θ) = R` is **constant in `θ`**. Hence the placement of `P_g`
relative to `P_f` does not depend on `θ` at all, and the predicate
`int P_f(θ) ∩ int P_g(θ) ≠ ∅` is `θ`-independent. It is false at `θ*` and true at `θ_k`: a
contradiction. So Case B cannot occur. ∎

**Remark [C] (Case B is in fact vacuous; kept only as a remark).** More is true than the
contradiction just derived: **two faces can never share two distinct material points**, so the
configuration Case B refutes does not exist in the first place. Indeed `u₁, u₂ ∈ H` are hinge
points, so by Sub-lemma T5.2b‴ `f` and `g` are hinge-adjacent and `σ_g = −σ_f`; by (0.7′) the
relative map `R(θ) = γ_g(θ)^{-1} ∘ γ_f(θ)` is then a rotation by `(σ_g − σ_f)θ/2 = ±θ` — which is
exactly T1.C, "every hinge opens by `θ`". For `θ ∈ (0, π]` a rotation by `±θ` is not the identity,
so it has exactly one fixed point and cannot fix both `u₁` and `u₂`. Equivalently: sharing two hinge
points would require two hinge edges between the same face pair, each keeping a different endpoint,
and the rotation forbids it. Measured: **0 / 933** face pairs share two `M′`-vertices (Checker
R4-b). The Case-B proof above is correct as printed; this remark is strictly stronger, and it is
why the case is retained as a remark rather than as a load-bearing branch.

**[D] The permanently collinear but non-coincident case, and how it is closed (round 5).**
Sub-lemma T5.2b″ treats the identically-zero harmonics through the T3.H.1/T3.H.2 classification:
`h_o,π ≡ 0` because `w ≡ a` or `w ≡ b` at every `θ`, a permanent *coincidence*. The logically
remaining possibility is a **permanently collinear but non-coincident** triple — `w` distinct from
`a` and from `b` at every `θ`, yet always on the line `ab`, i.e. `p = q = r = 0` with no copy
identification. Such a pair is struck from `𝒞` by the same identity test, and its incidence event is
not "`h` changes sign" but "`w` enters the *segment* `[a,b]`", which the orientation harmonic cannot
see. Round 4 named the two interval predicates below as the atoms that close this case; **that was
wrong, and the correct statement is the one proved here** (Checker R4.4).

**[D1] For hinge-adjacent pairs the case is impossible.** If `f` and `g` are hinge-adjacent then, as
in Remark [C], `R(θ) = γ_g(θ)^{-1} ∘ γ_f(θ)` is a rotation by `±θ` fixing the shared hinge point `p`.
So in `f`'s body frame a material point `w` of `g` traces the **circle** of radius `‖w − p‖` centred
at `p`. A circle is contained in a fixed line only if its radius is zero, i.e. `w = p` — a
coincidence, not this case. Hence permanent collinearity without coincidence cannot occur for
hinge-adjacent pairs, which is precisely where all **14 928** measured identically-zero harmonics
live (Checker R3.1a).

**[D2] For the remaining pairs, split the transition at the witness `z = w(θ*)`.** Keep the setting
of the sub-lemma: interiors disjoint at `θ*`, interiors overlapping at `θ_k → θ*` with `z_k → z`, and
**some** witness `π = (w,(a,b))` at `z` of the permanently collinear non-coincident kind.
**[R5 correction applied R6]** [D2] is stated here for a *single* witness, not for all of them, and
its proof uses only that one; Case A likewise needs only `z ∈ H`. This closes the **mixed** case —
both kinds of degenerate witness present at one `z` — which the round-5 headings ("every witness
at `z` is …") did not address as printed: if `z ∈ H` run Case A, otherwise pick any
permanently-collinear witness and run [D2] on it (Checker R5.5). By continuity of the finitely
many vertex positions there are a radius `r > 0` and a neighbourhood `I ∋ θ*` such that for `θ ∈ I`
the ball `B(w(θ), r)` meets no vertex of `f` other than through the edge `(a,b)`, and no edge of `g`
other than the two incident to `w`. **[R5 correction applied R6]** That phrasing is under-stated — it
does not exclude a *non-incident* edge of `f` crossing the ball. The condition actually used, and the
one the same compactness argument supplies (`f` simple, `w(θ*) ∈ relint(a,b)`, `I` compact), is
```
        B(w(θ), r) ∩ ∂P_f(θ)  ⊆  relint( a(θ), b(θ) )        for all θ ∈ I .
```
(Checker R5.1(ii).) Two sub-cases.

* **`w` at an endpoint of `[a,b]` at `θ*`.** Then `s₁(θ*) = 0` or `s₂(θ*) = 0` for the two interval
  predicates (T5.2b″-2) below, and those atoms do fire. This sub-case, and only this one, is closed
  by `s₁, s₂`.

* **`w` in the relative interior of `[a,b]` at `θ*`.** Then `s₁, s₂ > 0` on a neighbourhood of `θ*`
  and **neither interval atom has a root there** — they are silent at exactly the transition they
  were supposed to catch. What closes the sub-case is a *different* pair, already in the complete
  list: the **neighbour-vertex collinearity** `π″ = (w′, (a,b))`, where `w′` is a neighbour of `w`
  along `∂g`.

  *Proof of the interior sub-case.* Because `f` is simple and `w(θ)` lies in the relative interior of
  `[a(θ), b(θ)]` for `θ ∈ I` (permanent collinearity plus continuity), `P_f(θ) ∩ B(w(θ), r)` is the
  closed half-plane bounded by the line `ab`, and `P_g(θ) ∩ B(w(θ), r)` is the closed angular sector
  of `g` at its vertex `w`, bounded by the two edges `w → w′` and `w → w″`. Pass to directions at
  `w`: let `A(θ) ⊂ S¹` be the closed arc of the sector (length `α_g`, endpoints the directions of
  `w → w′` and `w → w″`) and let `J(θ) ⊂ S¹` be the open half-circle of directions pointing into
  `int P_f`, whose two endpoints are `± u(θ)` with `u(θ)` the direction of `b − a`. Then

  ```
        int P_f(θ) ∩ int P_g(θ) ∩ B(w(θ), r) ≠ ∅        ⟺        A(θ) ∩ J(θ) ≠ ∅ .   (T5.2b″-1c)
  ```

  **[R5 correction applied R6]** `A(θ)` is the **closed** arc, while the direction set of `int P_g`
  at `w` is the *open* sector, so the left side is literally `int A(θ) ∩ J(θ) ≠ ∅`. The two agree:
  `J` is open and `α_g > 0`, so if an endpoint of `A` lies in `J` then so do nearby interior
  directions. (T5.2b″-1c) as printed is therefore correct (Checker R5.1(i)).

  At `θ*` the left side is false, so `A(θ*) ⊆ K(θ*) := S¹ \ J(θ*)`, the closed half-circle with the
  same endpoints `± u(θ*)`. Suppose no endpoint of `A(θ*)` equals `± u(θ*)`. Then the compact set
  `A(θ*)` is contained in the *open* arc `K(θ*) \ {± u(θ*)}`; since `A` and `u` depend continuously
  on `θ`, `A(θ) ⊂ K(θ) \ {± u(θ)}` for all `θ` in a neighbourhood of `θ*`, so by (T5.2b″-1c) the
  interiors do not meet inside `B(w(θ), r)` there. But for large `k` the overlap witness `z_k` lies
  in `B(w(θ_k), r)` (as `z_k → z = w(θ*)` and `w` is continuous) with `θ_k` in that neighbourhood —
  a contradiction. Hence an endpoint of `A(θ*)` is `± u(θ*)`: one of the two edges of `g` at `w` is
  parallel to the line `ab` at `θ*`. That edge emanates from `w`, which lies **on** the line `ab`, so
  parallel means collinear: `w′(θ*)` lies on the line through `a(θ*), b(θ*)`, i.e.

  ```
        h_o,π″(θ*) = det( b(θ*) − a(θ*), w′(θ*) − a(θ*) ) = 0 ,      π″ = (w′,(a,b)) .
  ```

  It remains to see `h_o,π″ ≢ 0`, so that `π″ ∈ 𝒞-list` and `NOROOT` is contradicted. Suppose
  `h_o,π″ ≡ 0`. Then `w` and `w′` are both permanently on the line `ab`, so the edge direction
  `w′ − w` of `g` is permanently parallel to the edge direction `b − a` of `f`. By (0.7′) the former
  is `R(−σ_g θ/2)` applied to a fixed body direction and the latter is `R(−σ_f θ/2)` applied to a
  fixed body direction, so parallelism for all `θ` forces `(σ_f − σ_g)θ/2 ≡ 0 (mod π)` for all `θ`,
  i.e. `σ_f = σ_g`. But then `R(θ) = γ_g(θ)^{-1} ∘ γ_f(θ)` has rotation part the identity: it is a
  pure translation, so in `f`'s body frame the direction `u`, the arc `A` and hence the predicate
  `A ∩ J ≠ ∅` are all **constant in `θ`** (only the base point `w` slides, along the line `ab`). The
  local predicate (T5.2b″-1c) is then `θ`-independent, contradicting its being false at `θ*` and true
  at `θ_k`. So `h_o,π″ ≢ 0`, `π″` is an ordinary member of the `𝒞-list`, and it has a root at
  `θ* ∈ (0, ε)`. ∎

  This is Case A's own mechanism transplanted: the overlap begins when a bounding ray of `g`'s sector
  crosses `f`'s boundary line, and the atom that sees it is an ordinary orientation harmonic, not an
  interval atom.

**[D3] Consequences.** The correct statement is: the permanently collinear non-coincident case is
closed by the two interval atoms **together with** the neighbour-vertex collinearity pairs. The
latter cost nothing — they are already in the complete ordered list of T4.1b / Corollary T4.2′ with
`h ≢ 0` — so **if a pattern contains such a pair, the two interval atoms `E2` must be added to
`NOROOT` for that pair**, namely "no root in `(0, ε)`" for each of

```
        s₁(θ) := ⟨ w(θ) − a(θ) ,  b(θ) − a(θ) ⟩ ,
        s₂(θ) := ⟨ w(θ) − b(θ) ,  a(θ) − b(θ) ⟩ ,                                    (T5.2b″-2)
```

each of which is again of the form `p + q cos θ + r sin θ` (every vertex coordinate is
`cos(θ/2)C_v + sin(θ/2)S_v`, so any bilinear expression in them reduces to that form by the
half-angle identities), with coefficients quadratic in `t` — so the same three-class atom list of
T5.2b.2 decides them and the degree-≤ 4 bound is unchanged; the neighbour-vertex pairs add no new
atom shape at all. Measured status of the underlying classification: the Checker verified
**144 / 144** of the identically-zero pairs on `squares` are coincidences, and the classification of
the other 14 784 on its corpus is not reported; by [D1] this residue is in any case confined to
non-hinge-adjacent pairs, and by [D2] it is now covered rather than assumed.

With Sub-lemma T5.2b″ in hand, in every case some candidate pair of the `𝒞-list` has a root in
`(0, ε)`, contradicting `NOROOT`. Hence `θ* ∈ O`.

`O` is therefore clopen in the connected set `(0, ε)`, so `O = ∅` or `O = (0, ε)`. By
`NOOVERLAP(θ₁)` the point `θ₁ = ε/2` lies in `(0, ε) \ O`, so `O ≠ (0, ε)`, hence `O = ∅`. No two
interiors overlap anywhere on `(0, ε)`, and `Θ_max(X) ≥ ε` by definition (T4.1a). Since
`(0, θ₁) ⊆ (0, ε)` is disjoint from `O`, the germ condition `EMB(t)` holds as well. ∎

Four remarks on what the proof uses. It needs the **complete** candidate list (Lemma T4.2 /
Corollary T4.2′), which is why T5.2d says completeness is the theorem. It needs T2 — the motion is
defined and analytic on all of `(0, π]`, so there is no other way for the configuration to change.
It needs `NOROOT` to be the *deflated* test of T5.2b.2: on the class-3 pairs `h = p(1 − cos θ)` has
no root in `(0, ε)` by Lemma T5.1e, so deflating by `(1 − cos θ)` discards nothing the proof needs,
and on the class-2 pairs the discarded root sits at `θ = 0`, outside the open interval. And it needs
`NOOVERLAP` at a point strictly inside `(0, ε)`: without it the conclusion is only "the overlap
status is constant", which is compatible with `O = (0, ε)` — the configuration overlapping from the
first instant. That is exactly the failure mode T5.1's 51 counterexamples exhibit, and it is what a
`θ = 0` test cannot see.

**Compared with round 2's version of this proposition**, two things changed and one did not. The
hypothesis `EMB` is replaced by the pointwise, reportable `NOOVERLAP(ε/2)`, which now *implies* `EMB`
rather than being assumed alongside it; the proof is reorganized as a connectedness argument
(`O` clopen) instead of an infimum argument, which is what makes the single-point test sufficient.
The conclusion, `R(ε) ⊆ U(ε)`, is unchanged, and the round-2 consistency measurement below still
applies — the round-3 hypothesis is implied by the round-2 one whenever the probe angle used there
lies below `ε`.

**[A] T5.2b.1 The atom list, corrected to degree ≤ 4 (D5).** With
`A := p − q`, `B := 2r`, `C := p + q` (each **quadratic** in `t`), `g(τ) = C + Bτ + Aτ²`. Round 1
printed the disjunct `g(τ_v)·g(0) > 0` at the vertex `τ_v = −r/(p−q)`; clearing the denominator,

```
   g(τ_v) = ( C·A − r² ) / A = −(disc/4)/A ,     so   g(τ_v)·g(0) > 0  ⟺  −(disc/4)·C·A > 0 ,
```

which is degree **8** in `t`, not 4. The Checker's replacement, which I adopt, keeps everything at
degree ≤ 4 by testing the vertex position instead of the vertex value:

```
   BOTH(0,T)  :=  disc ≥ 0  ∧  A·g(0) > 0  ∧  A·g(T) > 0  ∧  A·B < 0  ∧  −A·B − 2A²T < 0 ,
                  disc = B² − 4AC ,                                                    (T5.1b)
```

(the last two atoms say `0 < −B/(2A) < T`, cleared by the positive factor `2A²`), and

```
   NOROOT_π(0,T)  :=  [ g(0)·g(T) > 0 ]  ∧  ¬ BOTH(0,T) .                              (T5.1c)
```

`g(0) = C` and `g(T) = C + BT + AT²` are quadratic in `t`; `g(0)·g(T)`, `disc`, `A·g(0)`, `A·g(T)`,
`A·B` and `A·B + 2A²T` are all of degree `≤ 4`. So `POS ∧ NOROOT` is cut out by polynomials of
degree at most 4.

**[A] T5.2b.2 One necessary addition: (T5.1c) must be DEFLATED at `τ = 0` — and there are
THREE structural classes, not two (round 3).** (T5.1c) is a correct
test for "no root in `(0,T)`" only when `g(0) ≠ 0`. But `g(0) = p + q = 0` says the triple
`(a, b, w)` is collinear in the **flat** state, which holds *identically on all of `𝕏`* whenever `w`
and `a` (or `w` and `b`) are two `M′`-copies of the same vertex of `M` — every permanent incidence
of T3.H.2, and every split-edge pair by T5.3's `p + q ≡ 0`. On such a pair (T5.1c) fails its very
first atom and reports a root inside the interval, when in fact the root is the one *at* `τ = 0`.
Left uncorrected, `R(ε)` would be **empty** for every pattern with split cuts. The fix is to divide
out the known factor `τ ∝ tan(θ/2)`. Round 2 divided it out **once** and stopped; the Checker
(`check.md` R2.5) showed by measurement that once is not always enough, and the Checker is right.
The complete case split is by the order of vanishing of `g` at `τ = 0`, i.e. by which of
`C = g(0) = p + q` and `B = g′(0) = 2r` vanish.

**[D] The three classes.** With `A := p − q`, `B := 2r`, `C := p + q` (each quadratic in `t`),
`g(τ) = C + Bτ + Aτ²`:

```
   class 1   C ≠ 0                     g has no root at τ = 0
   class 2   C = 0, B ≠ 0    g = τ(Aτ + B)      a SIMPLE root at τ = 0
   class 3   C = 0, B = 0    g = A τ²           a DOUBLE root at τ = 0, and no other root
```

Class 3 is `p = −q`, `r = 0`, i.e.

```
        h(θ)  =  p (1 − cos θ) .                                                       (T5.1e)
```

**[A] Lemma T5.1e (class 3 is never a contact event).** If `h(θ) = p(1 − cos θ)` with `p ≠ 0`, then
`h` has **constant sign on `(0, π)`** — indeed on `(0, 2π)` — namely `sign h = sign p`, and its only
zero on `[0, π]` is the endpoint `θ = 0`.

*Proof.* `1 − cos θ > 0` for every `θ ∈ (0, 2π)`, since `cos θ = 1` only at `θ ∈ 2πZ`. Hence
`h(θ) = p·(1 − cos θ)` has the sign of `p ≠ 0` there and vanishes nowhere in `(0, π)`. ∎

Two consequences. **(i)** `h` never *crosses* zero on the open interval, so by (T4.1)/(T3.3) the
ordered pair `(w,(a,b))` never becomes incident for `θ ∈ (0, π)`: **class 3 contributes no candidate
angle at all**, and the correct deflated atom for it is the constant `true`. **(ii)** Numerically the
single deflation of round 2 leaves `Aτ + B` with `B` at round-off, so `sign(−B/A)` — and therefore
membership in `(0,T)` — is decided by noise; that is exactly the ill-posedness `check.md` R2.5
reports, and dividing by `(1 − cos θ) ∝ τ²` removes it instead of papering over it. Deflation by
`τ²` is legitimate because `1 − cos θ > 0` on the open interval: dividing by a strictly positive
function changes no sign and removes no root of `h` in `(0, π)`.

**[A] The exact atom list, all three classes.** Fix `T := tan(ε/2)`. Then

```
   NOROOT_π(0,T)  :=   [ g(0)·g(T) > 0 ] ∧ ¬BOTH(0,T)                when C ≠ 0        (T5.1c)
                       ¬[ A·B < 0  ∧  −A·B − A²T < 0 ]                when C = 0, B ≠ 0 (T5.1d)
                       true                                           when C = B = 0
                                                                       and A ≠ 0      (T5.1e′)
```

**[D] The side-condition `A ≠ 0` in (T5.1e′), and what happens without it (round 4, `check.md`
R3.1a).** Lemma T5.1e is stated and proved for `p ≠ 0`, i.e. for `A = p − q = 2p ≠ 0`; (T5.1e′) is
valid only under that hypothesis and is now printed with it. The excluded case is
`A = C = B = 0`, i.e. `p = q = r = 0`, i.e.

```
        h_o,π ≡ 0        (the orientation harmonic vanishes identically on (0, π)) ,
```

for which **every** `θ` is a root and the atom `true` would be flatly wrong. Such pairs are not
hypothetical: they are the **permanent incidences** of T3.H.1/T3.H.2 — a hinge point is a vertex of
both incident faces at every `θ`, so the triple `(a, b, w)` with `w ≡ a` or `w ≡ b` is collinear at
every `θ`. The Checker measured **14 928** of them on its corpus, and the numeric class rule
`|C| ≤ tol ∧ |B| ≤ tol` puts **all 14 928 in class 3**, which is exactly why the side-condition has
to be printed.

**They are removed before the classification, not decided by it.** An identically-zero harmonic is
detected by the *identity* test of T3.H.1 — `|p| + |q| + |r| ≤ tol·scale`, a test on the
coefficients, not a root test — and the pair is then **struck from the candidate list 𝒞
altogether**. It is therefore **not** a candidate of (T5.1e′), and no atom of any class is evaluated
for it. Both `derivations/scratch/check_r3.jl` and `Kirigami/test/derivation_tests.jl` implement
exactly this (`h.scale() ≤ tol → skip`), which is why their mismatch counts are unaffected by this
correction; what the correction changes is that the exclusion is now *stated* rather than assumed.
The price is paid in the proof of Proposition T5.2b′, where the removed pairs must be shown to be
harmless — see the closedness step there.

with `BOTH(0,T)` as in (T5.1b). Every atom that survives is of degree `≤ 4` in `t`, and the third
branch contributes no atom at all, so the degree bound of T5.2b.1 is unchanged. Membership in
class 2 or class 3 is decided by the *combinatorics* `(M, σ)` exactly as before — class 2 is the
pairs in which `w` and `a`, or `w` and `b`, are two `M′`-copies of the same vertex of `M`; class 3
is the sub-case in which the flat configuration additionally satisfies `r ≡ 0`, which by (T3.2) is
again a statement about which `M′`-copies coincide, not about `t`. The case split is therefore fixed
once per pattern and the description stays quantifier-free.

Two honest caveats. **(1)** `g(0) = 0` can also happen *accidentally*, at particular `t` where three
otherwise-unrelated flat points are collinear. There the pair is routed to (T5.1d) by the numeric
test rather than by the combinatorics; since (T5.1d) is the correct test for a simple root at `τ = 0`
whenever `B ≠ 0`, this is not conservative but *exact*, and (T5.1c)'s round-1 behaviour (excluding
the point) is no longer needed. **(2)** The measurement below classifies pairs by the numeric tests
`|C| ≤ 1e−11·scale²` and `|B| ≤ 1e−11·scale²`, not by the combinatorial rule, because the program
works from `(p,q,r)` alone. The two agree on the structural classes by construction.

**[N] Measured, round 3 (`derivations/scratch/check_r3.jl`).** Round 2's number
("**0 / 2 363 380** at every `ε`") was **not a valid measurement and is withdrawn**: `check_r2.jl`
compared the deflated atom list against a reference that applied *the same* single deflation, so the
comparison could not see class 3 at all. `check_r3.jl` fixes this by computing the truth **without
the `τ` chart**: `h(θ) = p + q cos θ + r sin θ` is a degree-1 trig polynomial, so it is monotone
between its critical points (`tan θ = r/q`); splitting `(0, ε)` at those and comparing signs at the
break points decides "does `h` cross zero in `(0,ε)`" by evaluations of `h` alone. The sign just to
the right of `0` is taken from the first non-vanishing Taylor coefficient — `C`, else `B`, else `A`,
the last of which *is* Lemma T5.1e. Same corpus as round 2 (16 graphs, 197 shape-space samples):

| `ε` | harmonics decided | ambiguous (skipped) | class 1 / 2 / 3 | round-2 list (2 classes) | round-3 list (3 classes) |
|---|---|---|---|---|---|
| 0.200 | 2 145 387 | 47 | 2 070 850 / 73 446 / **1 138** | **153** mismatches, all class 3 | **0** |
| 0.020 | 2 145 387 | 47 | 2 070 850 / 73 446 / **1 138** | **153** mismatches, all class 3 | **0** |
| 0.001 | 2 145 387 | 47 | 2 070 850 / 73 446 / **1 138** | **153** mismatches, all class 3 | **0** |

**One asymmetry to state plainly.** The atom list decides "`h` has **a root** in `(0,ε)`", which is
what Proposition T5.2b′ needs — Lemma T4.2 delivers only `h_o,π(θ*) = 0`, not a sign change. The
round-3 truth decides "`h` **crosses** zero in `(0,ε)`", which is weaker. So a mismatch of the form
*atoms say root, truth says no crossing* (an interior tangency, `disc = 0`) would be **conservative**
and would not threaten T5.2b′; none occurred in this corpus. Class 3 is not such a case: there `h`
has no root in `(0,ε)` at all, by Lemma T5.1e, so the two predicates agree and the deflation is exact
rather than merely safe.

The mismatch count is `ε`-**independent**, which is the tell the Checker identified: the failures are
structural, not interval-dependent. The Checker measured the same phenomenon on its own corpus
(**143** of 1.82 M, `ε`-independent, all with `C = B = 0`); the two counts differ only because the
corpora and the ambiguity filters differ, and both go to zero under the same fix.

**Deciding number: with the class-3 deflation (T5.1e′) in place, 0 mismatches out of 2 145 387 at
every `ε ∈ {0.2, 0.02, 0.001}`, and 153 without it.** So the Checker's degree-≤ 4 list is right and
round 1's was wrong; the list needed **two** extra clauses, not one, which the two
systematically-degenerate classes of T3.H.5 / T5.3 force.

**[N] Consistency check on Proposition T5.2b′ (R2-D).** 197 shape-space samples (`X₀` plus Gaussian
null-space perturbations on 16 graphs), testing round 2's form `POS ∧ EMB ∧ NOROOT ⟹ Θ_max ≥ ε`
directly (the round-3 hypothesis `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT` is implied by it, so this remains a
valid consistency check of the conclusion):

| `ε` | samples meeting the hypothesis | samples with `Θ_max < ε` | violations |
|---|---|---|---|
| 0.200 | 4 / 197 | 45 / 197 | **0** |
| 0.050 | 6 / 197 | 45 / 197 | **0** |
| 0.020 | 7 / 197 | 45 / 197 | **0** |
| 0.005 | 39 / 197 | 45 / 197 | **0** |
| 0.001 | 104 / 197 | 45 / 197 | **0** |

The implication is not vacuous: the same 45 of 197 samples have `Θ_max < 0.001`, so they must be —
and are — excluded by the hypothesis at every row, and at `ε = 0.001` the hypothesis still fires on
more than half the samples. The Checker's own search (343
samples at `ε = 0.02`) likewise found no counterexample. This is a check on the proposition, not its
justification; the justification is the proof above.

**How conservative is it?** Very: `NOROOT` forbids *any* orientation root, including the grazes that
T4.2″ correctly steps over (the hexagon pattern's `θ₁ = π/3` is such a root), so `R(ε)` excludes
valid designs. `R(ε) ⊆ U(ε)` is an inner approximation and must be advertised as one — useful as a
certified-feasible region for an optimizer, not as a description of `U(ε)`.

### T5.2c Two corrections to `ideas/ranking.md` §2.3

**[A]** The proposed headline theorem there reads
"`U(ε)` is a **basic** semialgebraic set cut out by **quadrics**". Both adjectives are wrong:

1. **Not by quadrics.** `p, q, r` are quadratic in `t`; the conditions that actually cut out `U(ε)`
   are products and discriminants of these, i.e. **quartic** in `t`. Only (T5.1) is quadratic.
2. **Not basic.** A basic semialgebraic set is a finite *conjunction* of polynomial inequalities.
   "No root in an interval" for a quadratic is intrinsically a *disjunction* — see the `¬BOTH`
   clause of (T5.1c) and the `g(0) = 0` branch (T5.1d) — and T4.2″'s overlap selection adds another. `U(ε)` is a finite Boolean combination — semialgebraic,
   not basic. I could not find a conjunctive description and I do not believe one exists in general;
   **CONJECTURE (not-basic)**, stated as a negative I did not prove.

   The strongest *true* statement in this direction, and the one the paper should carry
   (**revised in round 2 for D2 and D4**), is:
   > `U(ε)` is a semialgebraic subset of the `2k` design coordinates. It contains an explicitly
   > described inner region `R(ε)` — Proposition T5.2b′ — cut out by a **quantifier-free Boolean
   > combination of polynomials of degree ≤ 4** together with one polygon-overlap predicate,
   > evaluated at the single angle `θ = ε/2`; the
   > atom list is finite and enumerable in `O(n²)`. Membership in `R(ε)` is therefore decidable in
   > closed form, with no forward-kinematics simulation and no bisection.

   Compared with round 1, three things are gone: the word `O(n)` and the reference to H-LOC (which
   is refuted, T4.5), the implication that the degree-≤ 4 description covers `U(ε)` itself rather
   than an inner region, and the omission of the embeddedness predicate. What remains is still a
   direct answer to 2026 §6 ("geometric characteristics … directly from the embedding"), and it does
   not overclaim.

### T5.2d Dependence on the candidate list

**[A]** Everything above is conditional on `𝒞(X)` being the
complete event set, which is Corollary T4.2′ — which in turn rests on Lemma T4.2 (pure plane
geometry, proved) and on T2 (no kinematic termination, proved). If the candidate list were only
sufficient, `U(ε)` would still be semialgebraic but the description would be an over-estimate. It is
worth saying in the paper that the *completeness*, not the harmonic form, is the theorem.

### T5.3 The A1 separation function, and the exact false-negative condition of 2026 Eq. (9)

**[A] Setup.** Let `e = {a, b}` be a split edge between faces `f, g` (`σ_f = σ_g`), `d := x_b − x_a`,
`Δu := u_g − u_f`. By T1.B the two duplicates of `e` are the *same vector* `R_f d` at every `θ`, and
by (T1.8) the second is the first translated by `2 s J Δu`. The signed separation is the orientation
harmonic of `(a′, b′; a″)`, `a′,b′ ∈ f`, `a″ ∈ g`:

```
   h(θ) = det( R_f d , 2 s J Δu )
        = 2 s [ c · det(d, J Δu)  −  σ_f s · det(d, Δu) ]      (using det(Ja, Jb) = det(a,b))
        = sin θ · det(d, J Δu)  −  σ_f (1 − cos θ) · det(d, Δu) ,
```

hence, in the notation of T3,

```
   ┌────────────────────────────────────────────────────────────────────────┐
   │   p  =  − σ_f det(d, Δu) ,     q  =  + σ_f det(d, Δu)  =  −p ,         │
   │   r  =  det(d, J Δu)  =  ⟨d, Δu⟩ .                          (T5.2)     │
   └────────────────────────────────────────────────────────────────────────┘
```

`p + q = 0` **identically** — the flat state is always a root, as it must be (the duplicates coincide
at `θ = 0`). Substituting into (T3.5):

```
        (1 + τ²) h(τ)  =  2 r τ  +  2 p τ²  =  2 τ ( r + p τ ) .                       (T5.3)
```

This is the adversary's `2t(R + Pt)` with `(P, Q, R) = (p, q, r)`, now with **closed forms for
`p` and `r`** that the adversary's file does not have.

**[N] Measured (check D3).** On 72 split edges over 16 graphs, the closed forms (T5.2) agree with the
generic (T3.2) coefficients, and `p + q = 0`, to a worst relative deviation of **4.37e−15**.

**[A] The false-negative condition.** 2026 Eq. (7)–(9) penalises the first-order separation
`z_i · n̂ = h′(0)` — in this notation, `h′(0) = r` — and nothing else. So Eq. (9) certifies
`r > 0`. From (T5.3) the *only other* root is

```
        τ*  =  − r / p ,        θ*  =  2 arctan( − r / p ) .                           (T5.4)
```

Given `r > 0`, `τ* > 0` iff `p < 0`, i.e. iff `σ_f det(d, Δu) > 0`. Therefore

```
   ┌──────────────────────────────────────────────────────────────────────────────────┐
   │  Eq.(9) FALSE NEGATIVE at split edge e   ⟺                                        │
   │     r = ⟨d, Δu⟩ > 0    ∧    p = −σ_f det(d, Δu) < 0                              │
   │     ∧  θ* = 2 arctan(−r/p) < Θ_max(X)   ∧  the interval test E2 holds at θ*      │
   │     ∧  the wedge/entering test of T4.3 holds at θ*.                       (T5.5) │
   └──────────────────────────────────────────────────────────────────────────────────┘
```

Eq. (9) is blind to `p` **by construction**: `p` is a second-order coefficient, and the energy is a
first-order test at `θ = 0` — the one configuration T2.3 shows is a branch point. This is exactly
adversary A1, sharpened by the last two clauses, which A1 and kill experiment K1c both omit: the
hexagon counterexample of T4.2 is *precisely* a `θ*` of this family that is a grazing contact and not
a collision. **K1c as specified in `ideas/ranking.md` will over-count false negatives.** Its stated
threshold is `≥ 3%`; the correction can only move the measured rate down, so a PASS at `≥ 3%`
*after* applying the two extra clauses is a valid PASS, and a PASS before applying them is not.

**A drop-in replacement for the Eq. (9) penalty.** Penalise `max(0, −p/r)` (equivalently
`τ*`) instead of `sign(r)` alone: smooth, same `O(|E_split|)` cost, and by (T5.2) it is a ratio of two
quadratic forms in `X`. This is A1's proposal with the coefficients made explicit.

### T5.H What can break

1. `Θ_max ≥ ε` in T5.2 is a statement about a *fixed* combinatorial candidate list. As `t` varies the
   list is fixed (it is combinatorial, not geometric), so `U(ε)` is well defined — but the *active*
   subset changes, which is why T6 refreshes it.
2. (T5.2) assumes `σ_f = σ_g`, i.e. a genuine split edge. Hinge-edge duplicates are not translates.
3. `Δu = u_g − u_f` for two faces adjacent across a split edge is a *path sum* difference along `Γ`,
   i.e. the *relative* potential of two `Γ`-neighbours, not an absolute one. If `Γ` has both faces
   in the same component, `Δu` is well defined; if not (split-cut cycle, T1.H.2), it is not. (Round
   1 called this "exactly the quantity H-LOC is about"; H-LOC is refuted (T4.5) and the two are in
   any case different — `Δu` is a one-step difference, H-LOC was about the absolute path sum.)
4. If `d ∥ Δu` then `p = 0`: no second root, the duplicates separate monotonically. If `⟨d, Δu⟩ = 0`
   then `r = 0` and the flat root is double — Eq. (9)'s test is degenerate exactly there.

### T5.Check

* Recompute (T5.2) independently as `det` of the two duplicate segments' endpoints at
  `θ ∈ {0.1, 0.5, 1.0}` and fit `p, q, r`; must match to `1e−12` (this is check D3).
* Verify `p + q = 0` on every split edge of every graph (D3 does).
* For K1c: report the false-negative rate **twice**, with and without the two extra clauses of
  (T5.5). The gap between the two numbers is itself a result.
* Verify T5.1's negative claim by finding an `X` with all `A_f > 0` and `Θ_max = 0` (two exist in the
  16-graph run; the Checker found 51 of 70; K1a should report how often).
* **Round 3:** re-run `derivations/scratch/check_r3.jl` with several `ε`. The **three-class**
  deflated list must report **0** mismatches out of 2 145 387 — a nonzero count means one of the two
  degenerate classes (T3.H.5) is being mishandled: mishandling `g(0) = 0` makes `R(ε)` empty on every
  pattern with split cuts, mishandling `g(0) = g′(0) = 0` makes membership noise-decided. The
  two-class list is expected to report 153, all in class 3, at every `ε`.
  Re-run `derivations/scratch/check_r2.jl` for R2-A/R2-C/R2-D; R2-D must report **0** violations;
  a violation would falsify Proposition T5.2b′ and, through Lemma T4.2, most of T4 with it.
* **Attack Proposition T5.2b′ directly** by searching for an `X` that is positively oriented, does
  not penetrate at `θ = 0⁺`, has no orientation-harmonic root in `(0, ε)`, and still has
  `Θ_max < ε`. The proof says none exists; 197 samples × 5 values of `ε` found none.

---

## T6 — Projection / range maximization

**[D] T6.1 The optimization.** Over the design coordinates `t ∈ R^{k×2}`:

```
   minimize   F(t)  =  − Θ_soft(t)  −  λ Σ_{f ∈ F} log A_f( X₀ + Φ t )  +  μ_reg ‖t‖²
   with       Θ_soft(t)  =  − (1/β) log Σ_{π ∈ 𝒜}  exp( − β · θ_π(t) ) ,                (T6.1)
```

where `𝒜 ⊆ 𝒞` is the current **active set**, `θ_π(t)` is the first admissible root of pair `π`
(clamped to `π` if none), `β > 0` is the softmin sharpness, `λ > 0` the barrier weight. `Θ_soft` is a
smooth under-estimate of `min_π θ_π` with `min_π θ_π − (log|𝒜|)/β ≤ Θ_soft ≤ min_π θ_π`. The log
barrier keeps (T5.1) satisfied along the whole path, which is what the paper's Eq. (9) has no term
for (`[F, F18]`: "the objective contains nothing that keeps the flat state embedded").

**[A] T6.2 Exact gradients by implicit differentiation.** Let `h_π(θ; t) = p + q cos θ + r sin θ` and
`θ_π(t)` a simple root. Differentiating `h_π(θ_π(t); t) = 0`:

```
   ┌────────────────────────────────────────────────────────────────────────────────┐
   │   ∂θ_π/∂t_i  =  −  ( ∂p/∂t_i + cos θ_π · ∂q/∂t_i + sin θ_π · ∂r/∂t_i )         │
   │                    ────────────────────────────────────────────────            │
   │                          ( − q sin θ_π  +  r cos θ_π )                 (T6.2)  │
   └────────────────────────────────────────────────────────────────────────────────┘
```

The denominator is `h_π′(θ_π)`, nonzero exactly when the root is simple. Because `p, q, r` are
quadratic polynomials in `t` (T3.2 + T1 Step 8), `∂p/∂t_i` etc. are explicit **affine** functions of
`t` — one matrix–vector product each, no finite differences and no adjoint solve. In the `τ`
variable, with `g(τ; t) = (p+q) + 2rτ + (p−q)τ²`,

```
   ∂τ_π/∂t_i  =  −  ( ∂_i(p+q) + 2τ ∂_i r + τ² ∂_i(p−q) )  /  ( 2r + 2(p−q) τ ) ,
   ∂θ_π/∂t_i  =  ( 2 / (1 + τ²) ) · ∂τ_π/∂t_i .                                         (T6.3)
```

Then `∂Θ_soft/∂t_i = Σ_π w_π ∂θ_π/∂t_i` with softmin weights
`w_π = exp(−βθ_π) / Σ_{π′} exp(−βθ_{π′})`, and
`∂/∂t_i [ −λ Σ_f log A_f ] = −λ Σ_f (∂A_f/∂t_i)/A_f` with `A_f` quadratic, so `∂A_f/∂t_i` affine.
Everything is closed form.

**[A] T6.3 Active set — a heuristic, with no size guarantee (D2).** Initialise `𝒜` with all pairs
whose first admissible root is below `θ_ref + margin`; refresh by re-running the exact broad phase
(T4.5a, sound with no hypothesis) plus the root enumeration every `n_refresh` iterations. Since
H-LOC is refuted (T4.5), **`|𝒜|` has no proved `O(n)` bound** and this must be described as an
optimizer heuristic whose *soundness* comes from the refresh, never as a "certified active set".
The refresh itself is exact, so the certificate of T6.4 is unaffected by how big `𝒜` gets. Refreshing is what keeps the objective honest;
(T6.1) with a stale `𝒜` maximises the wrong thing, and the failure mode is silent.

**[D] T6.4 The certificate (five items; round 1 had four, D11).** At the returned `t̂`,
`X̂ = X₀ + Φ t̂`, report:

1. **Deployability** — `‖L X̂‖_∞` (must be `0` up to the null-basis conditioning; `X̂ ∈ 𝕏` by
   construction, so this is a numerical, not a mathematical, check).
2. **Validity — `POS ∧ NOOVERLAP(θ₁) ∧ NOROOT` (revised in round 3, `check.md` R2.4b).** Report the
   three predicates of T5.2b.0 at one explicit `θ₁ = ε/2`:
   `POS` = `min_f A_f(X̂) > 0` (T5.1); `NOOVERLAP(θ₁)` = an exact polygon–polygon interior-disjointness
   test over all face pairs at the **single angle `θ₁ = ε/2`**; `NOROOT` = no admissible root of any
   candidate orientation harmonic in `(0, ε)`, by the deflated three-class atom list
   (T5.1c)/(T5.1d)/(T5.1e′). By Proposition T5.2b′ these three together **prove** `Θ_max(X̂) ≥ ε`,
   and they also imply `EMB` (the `θ = 0⁺` predicate), so the certificate now supplies exactly the
   hypothesis the proposition needs. Round 2 printed instead "no face–face interior overlap at
   `θ = 0`"; that is **withdrawn**. It is a different predicate from `EMB`, and it is degenerate:
   at `θ = 0` the two copies of every split edge coincide and faces adjacent across a cut touch
   along whole shared edges, so an exact interior-overlap test there is decided by tolerance rather
   than by geometry (T5.2b.0). Positive orientation alone remains **not** a validity certificate:
   51 of 70 positively-oriented samples have `Θ_max = 0` (T5.1, D11).
3. **Deployment validity** — `Θ_max(X̂) > 0` by T4.2″, computed independently of item 2. This is a
   *different* condition from `POS` alone, not a consequence of it (T5.1, D11), and both must be
   reported. Round 1 listed only this one.
4. **Range** — `Θ_max(X̂)` recomputed by T4.2″ over the **complete** candidate list, with no pruning
   and no active set, plus the binding pair. This is the number that goes in the paper, and it must
   be produced by the same routine for the baseline and for ours (as kill experiment K2b already
   specifies).
5. **Independence** — the same `Θ_max` cross-checked against `collision.jl`'s bisection, which uses
   no part of this derivation. Agreement to `1e−5` (deviation 9's tolerance), as measured in T4.2.

**T6.H What can break.**
1. `θ_π(t)` is **not** differentiable where the root is double (`disc = 0`) or where the active pair
   changes; `Θ_soft` is Lipschitz but not `C¹` across those events. The softmin smooths the second,
   not the first. A double root is where a contact appears or disappears — precisely where a
   range-maximizing optimizer wants to sit. **CONJECTURE (T6-smooth)**: that the non-smooth set is
   measure zero along the optimizer's path. Not proved; a subgradient/bundle fallback is the honest
   remedy and should be reported if line searches start failing.
2. The clamp `θ_π := π` when there is no admissible root has zero gradient, so a pair that is
   currently harmless contributes nothing — correct, but it means the active set must be refreshed or
   a pair can re-enter unannounced.
3. Face signed areas can be driven positive but tiny; `[F, 2026 §6]` names "extremely short edges or
   sharp angles" as the other half of its open problem. The barrier controls area, not aspect ratio.
   Adding an edge-length barrier is a change of objective, not of theory.
4. Nothing here proves the optimum is global. The claim can only ever be "certified valid, with a
   certified range, better than the baseline by a measured margin".

**T6.Check.** Finite-difference every `∂θ_π/∂t_i` from (T6.2) against a central difference of the
root recomputed from scratch, `rel. error ≤ 1e−6` at `Δ = 1e−5`, on ≥ 20 pairs on ≥ 5 graphs. Then
check that the reported `Θ_max` in the certificate is computed by T4.2″ **and** matches the
bisection — if the paper's headline margin were ever produced by `Θ_soft`, it would be circular.

---

## T7 — Rank and periodic corrections (short)

**Lemma T7 [A].** Let `D ∈ R^{|E_hinge| × N}` be the signed incidence matrix of the hinge digraph
(`D[e, dst(e)] = +1`, `D[e, src(e)] = −1`) and `R ∈ {0,1}^{H × |E_hinge|}` the hole/edge ownership
matrix, `R[K, e] = 1` iff `e ∈ C_K`. Then, with `[F, F11]` (the preimages partition
`E_hinge ⊔ E_split`):

**(i) `L = R · D`.** Row `K` of `L` is `Σ_{e ∈ C_K ∩ E_hinge} (x_{dst(e)} − x_{src(e)})`, which is
Eq. (2) verbatim, and is `Σ_e R[K,e] · D[e, ·]`. ∎ `[N, F22]` exact (zero difference) on 50/50 sweep
graphs, 8/8 reference cases, 3/3 tori.

**(ii) `leftnull(L) = Z := { y ∈ R^H : e ↦ y_{K(e)} is a circulation on the hinge digraph }`,
and `rank(L) = H − dim Z`.** `yᵀ L = (Rᵀ y)ᵀ D`, and `leftnull(D)` is the circulation space of the
digraph. `(Rᵀ y)_e = y_{K(e)}` requires each hinge edge to be owned by exactly one preimage, which is
F11. ∎

**(iii) Out-harmonic form.** By F11, `K(e)` is the split-forest component containing `dst(e)`, so
`w(e) = g(dst(e))` with `g(v) := y_{K(v)}`. The circulation condition at an interior vertex `v` reads
`indeg(v) · g(v) = Σ_{v → u} g(u)`; with Remark A.1 (`indeg = outdeg = k_v`),

```
        k_v · g(v)  =  Σ_{v → u} g(u)          (g is out-harmonic on the hinge digraph).   (T7.1)
```

**(iv) Periodic ⇒ `1ᵀ L = 0` and `rank(L) ≤ H − 1`.** In a boundary-free pattern every preimage is a
hole, so every hinge edge is owned by a row and `Rᵀ 1 = 1`. Then `1ᵀ L = 1ᵀ D = 0` because every
column of `D` sums to `indeg(v) − outdeg(v) = 0` (Remark A.1). So `1 ∈ Z`, `dim Z ≥ 1`,
`rank(L) ≤ H − 1`. ∎ `[N, F22]` `1ᵀL = 0`, `rank(L) = H − 1`, `dim Z = 1` measured on 3 tori;
out-harmonic residual `~1e−16`.

**Consequence.** 2026 §4.4's "the number of independent equations equals the number of holes" is
**false as stated for boundary-free (periodic) patterns**, off by exactly one. It is harmless in the
paper's pipeline because the boundary rows `B` pin the translation that `1` represents.

**(v) With a boundary the naive row sum fails; the corrected identity.** Preimages that touch the
mesh boundary are notches and contribute no row (`[F, F16]`, `code/README.md` deviation 2), so their
hinge edges drop out of `Rᵀ 1`. The corrected identity `(1ᵀL)_v = I(v)·indeg(v) − Σ_{v→u} I(u)`, with
`I` the indicator that a vertex's preimage is a hole row, is the `y ≡ 1` case of (iii). `[N, F22]`
holds 50/50 and 8/8; the naive form holds 0/50.

**(vi) Hole count.** `H = |E_hinge| − |F| + c(Γ)` by Euler on the planar graph `Γ` (bounded faces),
using the bijection holes ↔ bounded faces of `Γ` (geometer (0.5), the same statement as F11). On a
disk-topology `M` this equals `V_int − |E_split| + (c(Γ) − 1)`. `[N, U1/F22]` 50/50 sweep + 8/8
reference, counting holes only; including notches breaks it 0/50.

**T7.H (softened in round 2, D7).** (i)–(iii) all rest on F11. F11 is *measured* and *proved* by
planar duality in geometer (0.5); it is not proved here. Round 1 said that a user-supplied `σ`
creating a split-cut cycle makes "all of T7 void". **That is too pessimistic and is withdrawn.**
Under fully random `σ` on 398 `(graph, σ)` pairs — 210 of them with split-cut cycles, 334 with
`c(Γ) > 1` and 334 with `M′` disconnected — the Checker measured **zero** failures of: the F11
partition of `E_hinge ⊔ E_split` (**CE-a**), the hole-count formula `H = |E_hinge| − |F| + c(Γ)`
(**CE-b**), the equality of Step 6's row space with `L`'s (**CE-c**), and the reconstruction of
`deploy()` from (T1.5) built on a different spanning forest (**CE-d**, 7.5e−15). What actually
fails on a split-cut cycle is narrower, and this is the correct statement:

* **Fails:** the reading of an `L`-row as a *geometric* hole. `Γ` is disconnected, `M′` falls into
  `c(Γ)` independently placed pieces, and Def. 4.2's combinatorial preimages no longer correspond to
  holes of an assembled sheet (`[F, F16]`).
* **Survives:** `L = R D`, `rank(L) = H − dim Z`, the out-harmonic characterisation of `Z`, the
  periodic `1ᵀL = 0`, and `H = |E_hinge| − |F| + c(Γ)` — because all of these are statements about
  the cycle space of `Γ`, which does not care whether `Γ` is connected. The deployment theory then
  applies per component (T1.H.1), and the `2(c(Γ)−1)`-dimensional relative placement of the
  components is not in `𝕏` and is not a mechanism of an assembled structure.

**T7.Check.** Already covered by the rank-checks builder (37 tests / 6 836 assertions,
`results/core_validation/rank_claim.md`) and by K3b. Nothing new is asked for here; T7 exists so that
the paper can state F14/F22 as a lemma with a proof rather than as a measurement.

---

## T8 — What is proved, what is conjectured, and where to attack

### Proved (with the stated hypotheses)

| # | statement | status |
|---|---|---|
| T1 | `Y_θ = cos(θ/2)C + sin(θ/2)S`, `C, S` linear in `X`; path sums; path-independence ⟺ Eq. (2) | **proved** + [N] 1.4e−14 |
| T1.A–D | ellipse, split duplicates are equal translates, hinge opens by `θ`, faces rigid | **proved** + [N] |
| T2 | `A(Y_θ) σ = 0` for every `θ`; no kinematic locking; `A(θ)` is a pencil | **proved** + [N] 1.8e−14 |
| T3 | all orientation / dot / length predicates are `p + q cos θ + r sin θ`, coefficients quadratic in `X`; roots via a `τ`-quadratic | **proved** + [N] (and K1b PASS, 4.18e−12) |
| T4.2 | Lemma: contact ⇒ a vertex on the other's boundary; the candidate list is **complete** | **proved** |
| T4.2″ | exact `Θ_max` = candidate + interval scan | **proved** + [N] 2.1e−6 vs bisection |
| T4.4 | the 2026 §5.2 `β` bound is a special case; on split-free patterns `Θ_max = min(min_e β_e, π)` | **proved** + [N, R2-C] 8.88e−16 |
| T4.5a | exact `θ`-dependent and static broad phase (T4.2)/(T4.3), moving centroid | **proved** |
| T4.5b′ | in the face frame the swept radius is **exactly** the flat circumradius; `max(‖x‖,‖χ‖) = σ_max` | **proved** + [N, R2-A] 3.11e−15 |
| T5.1 | positive orientation = `|F|` quadratic inequalities; **not** a validity certificate | **proved** + [N] 2 counterexamples (Checker: 51/70) |
| T5.2a | `U(ε)` is semialgebraic | **proved** |
| T5.2b′ | `R(ε) = POS ∧ NOOVERLAP(ε/2) ∧ NOROOT` is an **inner** approximation: `R(ε) ⊆ U(ε)`, and implies `EMB` | **proved** (round 2, D4; restated round 3, R2.4b) + [N, R2-D] 0 violations |
| T5.1e | class 3 `h = p(1 − cos θ)` has constant sign on `(0, π)`: never a contact event | **proved** (round 3) |
| T5.2b.1/.2 | quantifier-free atom list of degree ≤ 4 for `POS ∧ NOROOT`, with the `τ = 0` deflation in **three** classes | **proved** + [N, round 3] 0 / 2 145 387 at every `ε` |
| T5.3 | `(1+τ²)h = 2τ(r + pτ)` with `p = −σ_f det(d,Δu)`, `r = ⟨d,Δu⟩`; exact Eq. (9) false-negative condition | **proved** + [N] 4.4e−15 |
| T6.2 | exact root gradients by implicit differentiation | **proved** + [N] 3.9e−7 (Checker **T6-a**) |
| T7 | `L = R D`, `rank(L) = H − dim Z`, out-harmonic `Z`, periodic `1ᵀL = 0`, `H = |E_h| − |F| + c(Γ)` | **proved** (modulo F11) + [N, F22] + 398/398 random `σ` |

### REFUTED in round 2 — must not appear anywhere

1. **H-LOC** (T4.5). Round 1 listed it as a conjecture. It is now **refuted** in the only form that
   supports an `O(n)` count: `max_f max_θ ‖γ_f(θ) − x̄_f‖ / r_f` grows **linearly** with patch
   diameter (4.62 → 20.03 as the diameter goes 5.66 → 19.80, **[N, R2-B]**; the Checker's
   independent table agrees). The other form of the statement, in the face-centroid frame, is
   trivially true with `κ = 1` and carries no information (Lemma T4.5b′). Consequences: there is
   **no `O(n)` active-set theorem**, only the exact `O(n²)` broad phase; the words `O(n)`,
   "certified active set" and "locality theorem" are removed from this file; and **kill experiment
   K2c is invalid as specified** — its gate ratio is identically `1.000000`, so it can neither pass
   nor fail informatively.
2. **The `√2` correction to K2c** (round-1 T4.5b). Withdrawn — the swept radius in the frame the
   pruning uses is exact, not a `√2` under-estimate (**[N, R2-A]** 3.11e−15).
3. **"Split-cut cycle ⟺ `Γ` disconnected"** (round-1 T1.H.2). Only `⟹` holds; 124 of 398 random-`σ`
   cases had `c(Γ) > 1` with a split forest.
4. **"A split-cut cycle voids all of T7"** (round-1 T7.H). Withdrawn; see T7.H for what actually
   fails and what survives.

### CONJECTURE — must not be called a theorem

1. **T4.3-generic**: first contact is vertex-into-edge-interior with a simple root, for `X` in a
   dense open subset. Provably **false** on symmetric patterns with congruent split-edge duplicates
   (hexagons, measured), so it can only ever be a statement about generic `X`, and the identical
   vanishing of the degeneracy polynomials on `𝕏` is not ruled out.
2. **T4.3-wedge**: a closed-form wedge test replacing the overlap probe at vertex-vertex contacts.
3. **not-basic** (T5.2c): that `U(ε)` admits no conjunctive (basic) description. Stated as a negative
   I did not prove; what *is* proved is that the natural description is a Boolean combination.
   The Checker agrees with this labelling.
4. **T6-smooth**: the non-differentiable set of `Θ_soft` (double roots) is negligible along the
   optimizer's path.
5. **F11** is used, not proved, in T1 Step 6 and throughout T7. Proved elsewhere (geometer (0.5)),
   measured 100/100, and now measured 398/398 under **fully random** `σ` (Checker **CE-a**), which
   is the hardest case anyone has tested; still not proved here.
6. **A deployment-frame locality statement** — bounding `‖γ_f(θ) − γ_g(θ)‖` relative to
   `‖x̄_f − x̄_g‖` — which would restore an `O(n)` packing count without the refuted flat-frame
   hypothesis (T4.5). Not attempted here; listed so it is not confused with H-LOC, which is dead.

### The weakest steps, for the next Checker to attack first

Round 1 named three. One of them (T1 Step 6) has since been attacked and survived, so it is
removed; the list is now two, plus the one new proof introduced in round 2.

1. **T4.2″ vs. `min over roots` (the hexagon row).** If the interval scan is wrong, `Θ_max` is
   wrong, K2a's tolerance is wrong, and K1c over-counts. Attacked in round 2 by the Checker with an
   independent separating-axis overlap test: hexagons `θ₁ = 1.047198`, **zero interior overlap at
   199/199 sampled angles across `(0, 2.094395)`**, overlap immediately after — and the Checker
   closed the one gap the statement had, by noting interior overlap is an **open** condition in `θ`,
   so the overlap set has no isolated points and the scan cannot skip one. **This step is now
   stronger than in round 1.** Remaining attack: find a pattern where the scan and the bisection
   disagree by more than `1e−5`.
2. **Proposition T5.2b′ (new in round 2).** The inner-approximation direction now carries the whole
   "usable region" claim, and it is a *proof*, not a measurement — the numeric check fires on only
   104/197 samples at `ε = 0.001`. Attack it at the two places it can break: (a) the compactness
   step, if some face's vertices could escape to infinity (they cannot — T1.A bounds them, but the
   uniformity in `θ` deserves a second reading); (b) the clopen argument — `O` closed in `(0, ε)`
   is the step that consumes `NOROOT` through Lemma T4.2, and it is where an incomplete candidate
   list would show up. Round 3 removed the earlier weak point (`EMB` as an unmeasurable germ at
   `0⁺`) by replacing it with the pointwise test `NOOVERLAP(ε/2)`; the finiteness of the candidate
   list still matters and is worth re-checking against degenerate (identically-zero) harmonics,
   which T3.H.1 removes by an identity test.
3. **F11**, which T1 Step 6 and all of T7 use and neither proves. It is now measured 398/398 under
   random `σ`, but the proof lives in `ideas/persona_geometer.md` (0.5), not here. If the paper
   states T7 as a lemma, that proof has to be imported and re-read, not cited.

**No longer on this list:** T1 Step 6 (attack run, failed — D7) and H-LOC (no longer a step of any
argument, because it is refuted — D2).

### Programs run

| program | build | result |
|---|---|---|
| `derivations/scratch/check_t1_t2.jl` | `julia --project=Kirigami derivations/scratch/check_t1_t2.jl` linking `Kirigami/src/core/*.jl` | C1 3.10e−14, C2 1.64e−14, C3 1.56e−13, C4 2.43e−14, C5 9.06e−14, C6 3.13e−14, C7 2.77e−10 on 30 graphs × 8 angles |
| `derivations/scratch/check_t4_t5.jl` | same | D1 1.18e−11 (rule 1e−5; the C++ prototype gave 2.09e−6), D2 8.88e−16, D3 3.33e−15 over 72 split edges, D4 719/1061 violations in the **raw** frame, ratio 1.4141 — see T4.5b: this says nothing about any pruning |
| `derivations/scratch/check_r2.jl` (round 2) | same, plus `Kirigami/src/core/*.jl`; run line in the file header, takes `ε` as `ARGS[1]` | R2-A 2.00e−15 (1 628 copies), R2-B drift/`r_f` 4.62→20.03 over diameters 5.66→19.80, R2-C 8.88e−16 (8 split-free patterns), R2-D **0** violations at every `ε ∈ {0.2, 0.05, 0.02, 0.005, 0.001}` on 196 samples (197 in the C++ prototype; the shape-space sample points differ, see `derivations/scratch/README.md`). **R2-E is withdrawn** — its deflated reference applied the same deflation as the atom list under test, so the comparison was circular (T5.2b.2) |
| `derivations/scratch/check_r3.jl` (round 3) | same corpus and seed; run line in the file header, takes `ε` as `ARGS[1]` | three-class deflated atom list vs a `τ`-chart-free crossing test: **0** mismatches of 2 111 336 at `ε ∈ {0.2, 0.02, 0.001}` (72 ambiguous, skipped); the round-2 two-class list gets **151** wrong, all in class 3, `ε`-independent; class sizes 2 037 472 / 72 798 / 1 138 (the C++ prototype, on different shape-space sample points, gave 0 of 2 145 387, 47 ambiguous, 153 wrong, class sizes 2 070 850 / 73 446 / 1 138; the class-3 count is structural and identical) |

All three are standalone `main()`s; none modifies anything under `code/` (which other agents own in
this phase). Build lines are in the file headers. The Checker's own program,
`Kirigami/test/derivation_tests.jl` (38 test sets / 150 191 assertions, 0 failures), was written
independently and is the cross-check for everything above; where its numbers and mine differ in
sample size but not in conclusion, both are quoted.

---

## Change log — round 6 (`check.md` "Round 6": the F32 amendment to `NOROOT`)

Scope: the single amendment proposed in `results/kill/jitter/cert_diagnosis.md` §4 — `NOROOT` counts
only **admissible** roots — together with the one gap that amendment opens in Sub-lemma T5.2b″ and
its consequences for T5.2b′, [D2] and [D3]. `Θ_max`, T3, T4 and T5.1 are untouched.

### R6.0 New standing hypothesis for T4–T6 (H6): no straight vertices

Every face polygon of `M` has interior angle `α ∉ {0, π, 2π}` at each of its own vertices. This is
implied by "simple" in the usual sense of a polygon whose vertices are genuine corners, but it was
never printed; step R6.4 below is the first place it is *used*, so it is stated here. A face with a
straight vertex is a polygon with a redundant vertex and can be normalized away by deleting it.

### R6.1 The amended T5.2b.0(iii)

Replace `(iii)` by

```
  (iii) NOROOT(t) : for every candidate pair π = (w,(a,b)) in the 𝒞-list, the orientation
        harmonic h_o,π has no ADMISSIBLE root in (0, ε), where a root θ is admissible iff it
        also satisfies the two interval inequalities E2 of T4.1b,

              0  ≤  ⟨ y_w(θ) − y_a(θ) , y_b(θ) − y_a(θ) ⟩  ≤  ‖ y_b(θ) − y_a(θ) ‖² ,

        both sides being harmonics in θ by T3.1/(T3.3)/(T3.4), and both inequalities being
        NON-STRICT (the segment of T4.1b is closed). Equivalently:  𝒞(X) ∩ (0, ε) = ∅ .
```

and **delete** the sentence "`(iii)` is stricter than the contact condition — it drops the two
interval tests `E2` and the overlap selection — so `R(ε)` is smaller than the exact region",
replacing it with: "`(iii)` drops only the overlap *selection* of T4.2″; it keeps the interval tests,
so `R(ε)` is smaller than the exact region only through `NOOVERLAP(θ₁)` being tested at one angle and
through grazes (R6.6)."

Nothing else in T5.2b.0 changes. In particular the `𝒞-list` is still the complete ordered list minus
the identically-zero harmonics, struck by the coefficient identity test of T3.H.1.

The closedness of `E2` is **load-bearing**, not a tolerance choice: R6.2 below produces the
admissible root with equality on the upper side whenever the two far-side edges have equal length,
which is the generic situation on a regular tiling (measured: **580 of 9 642** Case-A pairs sit
exactly at the far endpoint, `derivation_tests.jl` R6-a). `contact.jl` implements `E2` inclusively at
`tol = 1e−12·(|L2.p| + L2.amp())`, which is the right side to err on.

### R6.2 T5.2b′ survives: SOUNDNESS is preserved

**Proposition T5.2b′ (unchanged statement).** `POS(t) ∧ NOOVERLAP(ε/2) ∧ NOROOT(t) ⟹ Θ_max(X(t)) ≥ ε`,
with `NOROOT` now read as R6.1.

The proof is the round-3/4/5 proof with three insertions, one per place where it produces a
contradiction with `NOROOT`. Each must now produce an **admissible** root, not merely a root.

**(1) The main branch is admissible by construction.** A *witness at `z`* is, by the localization of
Lemma T4.2, "a vertex `w` of one polygon lying on the boundary of the other" — either `z` itself is
such a vertex, or `z` lies on a maximal shared segment `S ⊆ ∂P_f(θ*) ∩ ∂P_g(θ*)` each of whose two
endpoints is such a vertex. "Lying on the boundary" means lying on a **closed edge segment**
`[a(θ*), b(θ*)]`, which is exactly `E2` at `θ*` with non-strict inequalities. So every witness pair
of the main branch has an *admissible* root at `θ*`, and the branch "some witness has `h_o,π ≢ 0`"
contradicts the amended `NOROOT` verbatim. **No change to the printed argument is needed beyond this
sentence.**

**(2) Case A: the substitute pair's root at `β_e` is admissible.** Insert after (T5.2b″-1b):

> *Admissibility of the root (round 6).* At `θ = β_e` the far-side gap `β_e − θ` is zero, so the two
> far-side edges of `f` and `g` emanating from the shared hinge point `a = p` are collinear **and
> co-directed** (the gap closes; it is not the anti-parallel configuration, which is `β_e ± π`, cf.
> R5.4). Writing `u` for their common unit direction and `L := ‖b − a‖`, `L′ := ‖w − a‖` for the two
> far-side edge lengths (both `θ`-invariant, T1.D), we have at `θ = β_e`
> ```
>       b − a = L u ,      w − a = L′ u ,      0 < L′ ≤ L
> ```
> (`w` is by construction the far endpoint of the **shorter** of the two edges). Hence
> ```
>       ⟨ w − a , b − a ⟩ = L L′ ∈ (0, L²] ,        ‖ b − a ‖² = L² ,
> ```
> so `E2` holds at `β_e`, the lower inequality strictly and the upper one with equality **iff**
> `L′ = L` — the vertex-vertex case of T4.3, which the *closed* segment of T4.1b admits. The root
> `β_e` of `h_o,π′` is therefore admissible, and `π′ ∈ 𝒞-list` contradicts the amended `NOROOT`. ∎
>
> Measured (`derivation_tests.jl` R6-a1/R6-a2, 9 642 Case-A substitute pairs over the corpus with
> null-space jitter `0.12`, simple faces only): the projection ratio `⟨w−a,b−a⟩/‖b−a‖²` at `β_e`
> equals `L′/L` to **1.77e−14**, lies in `[0.007307, 1.000000]`, **0** admissibility failures,
> **9 062** strictly inside and **580** exactly at the far endpoint.

**(3) [D2], interior sub-case: the neighbour pair may be inadmissible, and is then replaced.** The
round-5 argument produces `π″ = (w′,(a,b))` with `h_o,π″(θ*) = 0` and `h_o,π″ ≢ 0`, but `w′(θ*)` need
not lie on the **segment** `[a(θ*), b(θ*)]` — the edge `w → w′` of `g` is collinear with the line
`ab` at `θ*`, and nothing bounds its length. Insert:

> *Admissibility (round 6).* Parametrize the line `ab` at `θ*` by arclength with `a = 0`, `b = L > 0`,
> and let `w = c` with `c ∈ (0, L)` (the interior sub-case) and `w′ = c + d`, `d ≠ 0`.
>
> * If `0 ≤ c + d ≤ L`, then `w′(θ*) ∈ [a(θ*), b(θ*)]`: `π″` itself has an **admissible** root at
>   `θ*`, and the round-5 conclusion stands unchanged.
> * If `c + d > L`, then `c < L < c + d`, so `b(θ*)` lies in the **relative interior** of the segment
>   `[w(θ*), w′(θ*)]`. The pair `π‴ := (b, (w, w′))` — a vertex of `f` against a consecutive vertex
>   pair of `g` — is in the complete ordered list of T4.1b with the roles of the two faces exchanged,
>   and `h_o,π‴(θ*) = 0` with `E2` satisfied **strictly**. It remains to see `h_o,π‴ ≢ 0`. Suppose
>   `h_o,π‴ ≡ 0`, i.e. `b` is permanently on the line `w w′`. Also `w` is permanently on the line
>   `ab` (the hypothesis of [D2]). The function `‖y_w − y_b‖²` is a harmonic (T3.4) that is not
>   identically zero, since `w(θ*) ≠ b(θ*)`; so it has at most two zeros in `(−π, π]` and `w ≠ b`
>   on a punctured neighbourhood `I′` of `θ*`. For `θ ∈ I′` the two distinct points `w(θ), b(θ)` lie
>   on both lines, so `line(a,b) = line(w,w′)` there, whence `w′` lies on `line(a,b)` for `θ ∈ I′`
>   and `h_o,π″` vanishes on `I′`; a first harmonic with infinitely many zeros is identically zero,
>   so `h_o,π″ ≡ 0` — which round 5 already refuted (it forces `σ_f = σ_g`, a pure relative
>   translation, hence a `θ`-independent local predicate). Contradiction, so `h_o,π‴ ≢ 0` and `π‴`
>   contradicts the amended `NOROOT`.
> * If `c + d < 0`, symmetrically `a(θ*)` lies in the relative interior of `[w(θ*), w′(θ*)]` and the
>   pair `(a, (w, w′))` does the same job by the same argument.
>
> In every branch some `𝒞-list` pair has an **admissible** root at `θ* ∈ (0, ε)`. ∎

**(4) [D2], endpoint sub-case: now VACUOUS, and `s₁, s₂` are no longer needed.** See R6.4.

With (1)–(4), every branch of the closedness argument produces an admissible root in `(0, ε)`, so
`O` is clopen in `(0, ε)`, `NOOVERLAP(ε/2)` fixes the constant, and `Θ_max ≥ ε`. **T5.2b′ is SOUND
under the amended `NOROOT`.**

### R6.3 What made the amendment safe in the first place

Admissible roots are a **subset** of roots, so `NOROOT_old ⟹ NOROOT_adm`: the old predicate is the
*stronger* hypothesis, and the amendment **weakens** the certificate. Soundness therefore does not
come for free — a weaker hypothesis has to be shown still sufficient, which is the whole content of
R6.2. It is sufficient because the proof never needed collinearity with the infinite line. The only
thing Corollary T4.2′ supplies, and the only thing the connectedness argument consumes, is that the
overlap status can change only at a **contact**; and a contact is by T4.1b a point of the **closed
segment**, i.e. an admissible root. Dropping `E2` gave `NOROOT` a root set strictly larger than
`𝒞(X)`, which cost designs (87.0% of the split-bearing rows with `Θ_max ≥ 1`, F32) and bought
nothing the proof used.

### R6.4 [D3] is WITHDRAWN: the two interval atoms `s₁, s₂` are unnecessary

Round 4/5 left the endpoint sub-case of [D2] — `w` at an endpoint of `[a,b]` at `θ*`, with
`h_o,(w,(a,b)) ≡ 0` — to the two interval atoms (T5.2b″-2). That sub-case is **impossible**, so the
atoms may be dropped from `NOROOT` entirely.

*Proof.* Keep [D2]'s setting and say `w(θ*) = a(θ*) =: z` (the case `w(θ*) = b(θ*)` is symmetric).
Then `a`, a vertex of `f`, coincides with the vertex `w` of `g`, so `a` lies on `∂P_g(θ*)` — on both
edges of `g` incident to `w`. Hence **both** `(a, (w, w′))` and `(a, (w, w″))` are witnesses at `z`
(`w′, w″` the two neighbours of `w` along `∂g`). By the hypothesis of [D2] every witness at `z` is
permanently collinear and non-coincident, so `h_o,(a,(w,w′)) ≡ 0` and `h_o,(a,(w,w″)) ≡ 0`: the point
`a` lies on the line `w w′` and on the line `w w″` at **every** `θ`. By hypothesis (H6) the interior
angle `α_g` of `g` at `w` is not `0`, `π` or `2π`, and `α_g` is `θ`-invariant (T1.D), so those two
lines are distinct at every `θ` and meet only at `w`. Hence `a(θ) = w(θ)` for all `θ`: a *permanent
coincidence*, contradicting the non-coincidence hypothesis of [D2] (and landing in the T3.H.2
classification, i.e. in Case A). ∎

*Consequence.* [D3]'s prescription — "if a pattern contains such a pair, the two interval atoms `E2`
must be added to `NOROOT` for that pair" — is withdrawn. It is replaced by: the permanently collinear
non-coincident case is closed by the **neighbour-vertex pairs alone** ([D2] with R6.2(3)), and its
endpoint sub-case does not occur. The interval inequalities `E2` still appear in the certificate, but
in their T4.1b role — as the **admissibility filter on the roots of the orientation harmonics** — and
not as separate atoms whose own roots are hunted. This is what `contact.jl` implements, so the
derivation and the code now agree; before this round the code was missing an atom the derivation
demanded, and demanding an atom the code was right not to have.

### R6.5 The degree bound and semialgebraicity are unaffected

`E2` is two harmonic inequalities whose coefficients are quadratic in `t` (T3.1, T3.2), evaluated at
a root of another harmonic with the same property. The set `{ t : π has an admissible root in
(0, ε) }` is still the `τ`-projection of a first-order formula with polynomial atoms, so T5.2a stands
verbatim. The degree-≤ 4 statement was already made only about `POS ∧ NOROOT`'s atom list and is
unchanged in shape: the admissibility filter adds sign conditions on `h_d,π` and `L²_π` evaluated at
the root, which is one more quantifier-elimination step of the same kind, not a new atom shape. **No
improved degree bound is claimed.**

### R6.6 `R(ε) ⊆ U(ε)` is still STRICT: completeness is FALSE

The jitter re-measurement (`cert_diagnosis.md` §6: 5 064 certified = 5 064 rows with
`Θ_max ≥ ε`, 0 exceptions in either direction at `ε = 0.006`) invites the converse statement. **It is
false**, and core.md's own T4.2 table is the counterexample.

**[A] Proposition T5.2b⁗ (exactly what the amended certificate decides).** Given `POS(t)` and
`NOOVERLAP(ε/2)`,

```
        NOROOT_adm(ε)   ⟺   ε  ≤  θ₁(X) := min 𝒞(X)                                  (T5.2b⁗-1)
```

(by R6.1, `NOROOT_adm(ε)` *is* `𝒞(X) ∩ (0, ε) = ∅`), while by (T4.1) `Θ_max(X) = θ_{i*}` is the first
`𝒞`-angle at which the overlap status actually turns on. Since `θ₁ ≤ θ_{i*}` always, the certificate
decides `ε ≤ θ₁` where `U(ε)` asks for `ε ≤ Θ_max`. Hence

```
        R(ε) ⊆ U(ε) ,   with equality on a given X iff θ₁(X) = Θ_max(X) ,             (T5.2b⁗-2)
```

i.e. **iff the first contact of `X` is already an overlap transition** — no *graze* below the range
end. That is precisely Proposition T4.3's genericity hypothesis, which T4.3 itself records as failing
identically on the symmetric patterns both papers ship.

**[N] Measured counterexample (`derivations/scratch/check_t5_adm.jl`, and `derivation_tests.jl` R6-b).**
`hexagons` at `X₀`:

| quantity | value |
|---|---|
| `θ₁ = min 𝒞(X)` | `1.047198` (= `π/3`, the split-duplicate end-to-end contact of the T4.2 table) |
| `Θ_max` by (T4.1) | `2.094395` (= `2π/3`) |
| `ε` chosen | `1.570796` (= `π/2`, strictly between) |
| `Θ_max ≥ ε` | **true** |
| certificate `(POS, NOOVERLAP(ε/2), NOROOT_adm)` | `(1, 1, 0)`, first admissible root `1.047198` |

So `X₀ ∈ U(ε) \ R(ε)`: the amended certificate is still a **strict inner approximation**. The
`0 / 5064` of the jitter run says only that `ε = 0.006` sits far below `θ₁` on every row of that
corpus — the smallest contact angle measured anywhere on the T4.2 table is `7.4e−4`, and on the
regular tilings it is `≥ 0.83` — not that grazes are absent. **Any sentence claiming the certificate
is exact must be scoped to a corpus and an `ε`, never stated as a property of the certificate.**

**What a complete certificate would need.** Exactly the overlap *selection* of T4.2″: one overlap
probe per interval of `(0, ε) \ 𝒞(X)` instead of one probe at `ε/2`. That is `exact_theta_max_overlap`
restricted to `(0, ε)`, and it is a finite exact computation — but it is no longer a fixed Boolean
combination of polynomial atoms in `t` with a bounded structure, because the number of intervals is
itself a function of `t`. The trade is stated, not resolved: **CONJECTURE (T5.2b-complete)** — that
`U(ε)` admits a quantifier-free description whose atom count is bounded independently of `t`.

### R6.7 The five round-5 corrections, now applied to the body text

`check.md` R5.6 listed five corrections that were recorded but never applied to this file's theorem
text. All five are now in place, each marked **[R5 correction applied R6]** at the point of use:

| # | where | what changed |
|---|---|---|
| 1 | T5.2b″ Case A, after (T5.2b″-1b) | the clause "**exactly one** of `{β_e, β_e ± π}` lies in `(0, π)`" and its three-way case split are **deleted** — they are false whenever `β_e ≤ −π` (5 008 of 20 000 synthetic values, Checker R5.4). Replaced by "at most one zero in `(0, π)`, and it is `β_e` here since `β_e = θ* ∈ (0, ε)`". No theorem was affected: the failing regime needs `α_f + α_g ≥ 3π`. |
| 2 | T5.2b″ [D2], after (T5.2b″-1c) | `A(θ)` is the **closed** arc while the direction set of `int P_g` is the open sector; the clause `A ∩ J ≠ ∅ ⟺ int A ∩ J ≠ ∅` (`J` open, `α_g > 0`) is added, which is what makes the printed equivalence correct (Checker R5.1(i)). |
| 3 | T5.2b″ [D2], the localization | "`B(w(θ), r)` meets no vertex of `f` other than through `(a,b)`" is under-stated (it does not exclude a non-incident edge of `f` crossing the ball). The condition actually used is now printed: `B(w(θ), r) ∩ ∂P_f(θ) ⊆ relint(a(θ), b(θ))` for `θ ∈ I` (Checker R5.1(ii)). |
| 4 | T5.2b″ [D2], the hypothesis | [D2] is now stated for a **single** permanently-collinear witness rather than for all witnesses at `z`, which closes the **mixed** case (both kinds of degenerate witness at one `z`) that neither Case A's nor [D2]'s printed heading covered (Checker R5.5). |
| 5 | Sub-lemma T5.2b‴ [A0] | the step "the two edges of `f` at `v` are consecutive in the `σ_f`-traversal" is now labelled as the **simplicity hypothesis** — it says the traversal visits `v` once (Checker R5.3). |

None of the five changes a hypothesis, a proof strategy or the truth value of any theorem; items 2–5
are the wording the Checker prescribed and item 1 removes a false sentence.

### R6.8 The B4 `voronoi_93` scan-vs-referee disagreement is a REFEREE artefact

`results/kill/b4/b4.csv` row (id 93, voronoi, `σ_def`, `free`) reports the exact T4.2″ scan and the
repaired certificate at `Θ_max = 0.2484` and the bisection referee at `0`. Replayed in
`derivations/scratch/check_b4_93.jl`: **the scan is right.**

`has_collision(..., 1e−12)` fires on face pair `(94, 184)` for `θ ∈ (0, 0.036255)`, and
`referee_theta` returns `0` on its `col(1e−7)` probe. The two faces share exactly one `M′`-vertex,
a hinge point, with `α_f = 0.368615`, `α_g = 4.105228`, hence `β_e = 1.809342` — which is exactly the
first contact angle the scan reports for that pair, in agreement with (T5.2b″-1). Both faces are
simple and positively oriented, and a `1 200 × 1 200` grid finds **0** points inside both against
73 735 / 342 350 inside each, at every probed `θ` below the transition. There is no interior overlap,
so `𝒞(X)` is not missing anything and Corollary T4.2′ is not violated.

The referee's `polygons_overlap` displaces the **shared hinge vertex** by `−s (p − centroid)`,
differently in each of the two copies, so the two edges meeting exactly there acquire orientation
determinants of size `1e−11 … 1e−17` whose signs decide the strict crossing test. Every spurious
crossing in the replay is hinge-incident, and the answer flips with `shrink` and with `θ` for no
geometric reason. Shrinking each polygon toward **its own** centroid does not separate a shared
vertex when one of the two sectors is **reflex**.

Consequences for this file: none for any theorem — T4.2″ and T5.2b′ are confirmed, not contradicted.
Consequence for the measurements: `theta_bisect`, `worst_referee_gap` and the `col(1e−7)` zero-range
shortcut shared by K2a/K2b/K5/K6 can report a spurious `Θ_max = 0` on designs with a thin-against-
reflex hinge. The bias is **downward only**, so no PASS can be manufactured and the K5/K6 FAIL
verdicts (which rest on `zero_plus_q ≤ 0` and on the exact scan) stand; but **`theta_bisect` must not
be quoted as ground truth against the exact scan**. This is *not* an admissibility-filter false
negative: the pair's genuine on-segment contact at `β_e` is kept by the filter, and `NOROOT` is false
on this row for an unrelated reason.

### R6.9 Programs and tests run for round 6

| program | result |
|---|---|
| `derivations/scratch/check_t5_adm.jl` (new; run line in its header) | A: 1 363 Case-A pairs on the 7 tilings with jitter, ratio `∈ [0.5086, 1.0]` (C++ prototype: `[0.5404, 1.0]`), `|s/‖e‖² − L′/L| ≤ 1.55e−15`, **0** admissibility failures. B: 1 completeness counterexample (`hexagons`, `ε = 1.570796`) |
| `derivations/scratch/check_b4_93.jl` (new; run line in its header) | replays B4 ids 93 and 96 from the K6/B4 caches; 93 reproduces scan `0.248400` vs referee `0`, and localizes the referee's misfire to the shared hinge vertex of faces (94, 184) |
| `Kirigami/test/derivation_tests.jl` R6 + R6-c cases (new) | R6-a1 0 failures of 9 642 pairs · R6-a2 `1.77e−14` · R6-a3 `3.57e−14` · R6-b counterexample found · R6-c1 0 common interior points of 1 442 401 · R6-c2 misfire pinned. Whole file: **38 test sets / 150 191 assertions, 0 failures** |
