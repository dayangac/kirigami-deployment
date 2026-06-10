# check.md — adversarial verification of `derivations/core.md` (Phase 7 Checker)

Independent review with fresh context. I saw only `derivations/core.md`, the code, `STATE.md`
F1–F24 and `notes/paper_2026.md` — not the Deriver's reasoning. Every algebraic step was
re-derived by hand from the code's own sign convention (`code/src/core/kinematics.cpp`), and
every identity was re-tested by a program I wrote from scratch,
`code/tests/derivation_tests.cpp`, which recomputes the face potential `u`, the harmonic
coefficients `(p,q,r)`, the `τ`-quadratic, the T5.3 closed forms, the T6.2 gradient and the
`Γ`-cycle closure rows **without using the Deriver's scratch programs** (`check_t1_t2.cpp`,
`check_t4_t5.cpp` were never compiled or run by me).

Build and run — one command, from the repo root:

```
clang++ -std=c++20 -O2 -arch arm64 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
  -Icode/src code/tests/derivation_tests.cpp code/build/libkiri_core.a \
  -o code/build/derivation_tests && ./code/build/derivation_tests
```

`code/CMakeLists.txt` and every existing file under `code/` are untouched. Result:
**25 test cases, 29 033 assertions, 0 failures**; 41 recorded identities, all PASS.
Corpus: 16 `(graph, σ)` cases (7 tilings + 9 random Delaunay/Voronoi/quad, 9–102 hinge edges,
0–30 split cuts, `dim_null` 0–30), plus 398 random-`σ` pairs and 3 torus patches.

---

## 0. Verdict summary

| theorem | verdict |
|---|---|
| T1 (trig-linear deployment), T1.A–T1.D | **VERIFIED** |
| T1 Step 6 (closure ⟺ Eq. 2) | **VERIFIED**, and far more robustly than claimed — see D7 |
| T2 (no locking, `A(θ)σ = 0`, pencil) | **VERIFIED** |
| T3 (harmonic predicates, `τ`-quadratic) | **VERIFIED** (numerical caveat D10) |
| T4.2 / T4.2′ (completeness of the candidate list) | **VERIFIED** (proof re-checked by hand) |
| T4.2″ (exact `Θ_max` by interval scan) | **VERIFIED**; graze reproduced independently |
| T4.3-generic, T4.3-wedge | **CONJECTURE** (as the Deriver states) |
| T4.4 (the `β` bound) | **VERIFIED with a correction** — D3 |
| T4.5a (broad phase) | **VERIFIED** for the moving-centroid form; the flat-centroid variant is unsound in principle — D8 |
| T4.5b (swept radius, the `√2` correction to K2c) | **DISPUTED** — D1 |
| H-LOC | **DISPUTED** — the stated equivalence is false, and the form the `O(n)` argument needs is **refuted** — D2 |
| T5.1 (positive orientation quadratic; not a certificate) | **VERIFIED**; certificate wording incomplete — D11 |
| T5.2a (`U(ε)` semialgebraic) | **VERIFIED** |
| T5.2b (explicit degree-≤ 4 conservative description) | **DISPUTED** — D4 (missing atom), D5 (atom degrees) |
| T5.2c (not basic, not quadrics) | **AGREE**; "not basic" stays a **CONJECTURE** |
| T5.3 (split-edge closed forms, Eq. (9) false negatives) | **VERIFIED** |
| T6.2 (implicit-differentiation gradients) | **VERIFIED** (finite differences, 2 160 samples) |
| T6-smooth | **CONJECTURE** |
| T7 (i)–(vi) | **VERIFIED**; the F11 caveat is weaker than stated — D7 |

**Eleven disagreements**, D1–D11 below. Two are substantive enough to change what the paper can
claim: **D1/D2 (the locality story) and D4 (the conservative region is not proved to be inner).**

---

## 1. Step-by-step verification

Notation as in `core.md`. `J` is rotation by `+π/2`; I use `det(Ja, b) = −⟨a,b⟩ = −det(a, Jb)`
and `det(Ja, Jb) = det(a, b)` throughout, both verified by direct expansion.

### §0 conventions

| step | verdict | reason |
|---|---|---|
| 0.2, 0.4–0.6 | AGREE | definitions; match `cut.hpp` / `mesh.hpp` verbatim |
| 0.7, 0.7′ | AGREE | `kinematics.cpp` sets `a = -sigma[f]*theta*0.5`; `R(a) = cI + sJ` expands to `cI − σ_f s J`. Confirmed numerically: the closed form reproduces `deploy()` to 1.4e−14, the sign-flipped form fails by 45 (test **T1-a**) |
| 0.8 | AGREE | the geometer/rigidity persona convention is `θ → −θ`; the note that `σ → −σ` is *not* the same relabelling is correct, because `hinge_dir` is recomputed from `σ` |

### T1

| step | verdict | reason |
|---|---|---|
| 1 | AGREE | definition |
| 2 (T1.2) | AGREE | re-derived: `t_g − t_f = (R_f − R_g)x_v = s(σ_g − σ_f) J x_v = 2 s σ_g J x_v` using `σ_f = −σ_g` |
| 3 (T1.3) | AGREE | `t_f = 2sJu_f` is a legitimate substitution for `s ≠ 0`; gives `u_g − u_f = σ_g x_{src}`, `θ`-free |
| 4 | AGREE | reversing gives `σ_f x_{src} = −σ_g x_{src}`; a genuine 1-form |
| 5 | AGREE | standard: a prescribed edge increment on a connected graph integrates iff its circulation vanishes on a cycle basis. The "one `θ` ⟹ every `θ`" corollary follows because `δ` is `θ`-free |
| 6 | **AGREE, and stronger than claimed** | The alternating-sum telescoping is correct (`Γ` bipartite by 0.5, `indeg_h = outdeg_h` by Remark A.1). The Deriver flags this as leaning on F11 and lists it as a top-3 weak step. I tested it directly by building the `Γ`-cycle closure rows symbolically and comparing **row spaces** with `L`: identical on 16/16 corpus cases and **398/398 random-`σ` cases**, including 334 with `c(Γ) > 1` and 210 with split-cut cycles (**CE-c**). The proposed attack fails. See D7 |
| 7 (T1.5) | AGREE | `y = R_f x_v + t_f = c x_v + sJ(2u_f − σ_f x_v)`. Verified 85 806 samples, 1.4e−14 (**T1-a**), and equals `deploy_basis()` to 7.3e−15 (**T1-b**) |
| 8 (T1.6) | AGREE | path sum with coefficients in `{0,±1}`; linearity in `X` used everywhere downstream |
| T1.A | AGREE | `y = A(c,s)ᵀ` on the unit circle; (T1.7) is the standard singular-value formula. Ellipse identity holds to 2.0e−9 over 46 890 samples (**T1-A**) |
| T1.B / (T1.8) | AGREE | the `2u_f` terms cancel in the intra-face difference; offset `2 s J Δu`. 2.2e−14 / 1.9e−14 (**T1-B1/B2**) |
| T1.C | AGREE | `R_g R_f^{-1} = R(σ_f θ)`. Magnitude `θ` verified to 1.5e−13 over 27 990 hinge evaluations (**T1-C**) |
| T1.D | AGREE | proper isometry; areas constant to 5.0e−14 (**T1-D**) |
| T1.H.1 | AGREE | verified: with `c(Γ) > 1`, (T1.5) still reproduces `deploy()` per component — 398/398, 7.5e−15, using a **deliberately different spanning forest** (**CE-d**) |
| T1.H.2 | **DISAGREE (D6)** | "`E_split` contains a cycle iff `Γ` is disconnected" — only ⟹ holds. 210/210 split-cycle cases had `Γ` disconnected, but **124** cases had `Γ` disconnected with a split *forest* |
| T1.H.3 | AGREE | verified as a control: off the shape space the closure defect is ≥ 0.41 and `deploy()`'s `max_mismatch` equals `2|sin(θ/2)|` × that defect to 1.9e−15 (**T1-c/T1-d**) |

### T2

| step | verdict | reason |
|---|---|---|
| T2.1 (T2.2) | AGREE | the pin constraint and the circulation formulation are the standard body-and-pin model; matches `mobility.hpp::build_A` |
| T2.2 (T2.3, T2.4) | AGREE | re-derived independently. `V(y) = ω_f J y + w_f` with `w_f = t_f' − ω_f J t_f`; `t_f' = cJu_f`, `J t_f = −2s u_f`, so `w_f = cJu_f − σ_f s u_f`. The residual collapses to `s σ_g(u_f − u_g) + s x_v = 0` using only (T1.3). No genericity, no `θ` restriction — correct. **Task 1(d): verified at every `θ` with the code's sign convention**: pin residual 1.3e−14 over 22 392 hinge evaluations at random `θ ∈ [−4,4]` (**T2-a**), and `‖A(Y_θ)σ‖∞/‖A‖∞ ≤ 1.0e−14` over 480 assembled matrices including `θ = 0` (**T2-b**) |
| T2.3 (T2.5) | AGREE | `p_e(θ)` is affine in `(c,s)` and `A` is linear in `p`. Pencil identity 4.0e−15 (**T2-c**) |
| T2.4 | AGREE | two hinged faces have distinct `ω`, so the flex is non-rigid whenever `E_hinge ≠ ∅`; the "every termination is a contact" corollary is what makes T4's candidate list complete |
| T2.H.1–3 | AGREE | correctly scoped to the pin-joint idealisation and to infinitesimal flexes |

### T3

| step | verdict | reason |
|---|---|---|
| T3.1 (T3.2)–(T3.4) | AGREE | re-derived by expanding `det(cC_ab + sS_ab, cC_aw + sS_aw)` with `c² = (1+cos)/2`, `s² = (1−cos)/2`, `cs = sin/2`. All three coefficient triples match. Verified against direct evaluation to 3.6e−15 relative over 12 345 samples (**T3-a/b/c**) and against `orient_harmonic()` to 4.5e−14 (**T3-d**) |
| T3.2 | AGREE | `C, S` linear in `X` (T1 step 8) ⟹ `p,q,r` quadratic forms. Independently confirmed downstream: the third central difference of `A_f(X₀+Φt)` in `t` vanishes to 3.4e−14 (**T5-a**) |
| T3.3 | AGREE | `S_ab = −σ_f J C_ab` inside a face; `q = 0` and `r = 0` follow from `det(Ja,Jb) = det(a,b)` and `det(Ja,b) = −det(a,Jb)`. Measured 1.7e−14 / 1.8e−14, and `p` equals the flat value (**T3-e/f/g**); `q″ = r″ = 0` likewise (**T3-h**) |
| T3.4 (T3.5), (T3.6) | AGREE | `(1+τ²)h = (p+q) + 2rτ + (p−q)τ²`; roots `[−r ± √(q²+r²−p²)]/(p−q)`, real iff `p² ≤ q²+r²`; the `p=q` branch with the extra root at `θ = π`; `p=q, r=0` gives `2p cos²(θ/2)`. All four cases re-derived and verified: every root satisfies `h = 0` to 8.4e−14, and the real-root criterion is exact (971 of 3 682 harmonics correctly reported rootless) (**T3-i**). Discriminant quartic in `t`: correct |
| T3.H.1–3 | AGREE | degenerate triples must be removed by an identity test, not a root test — confirmed by the code's `tol` argument |
| — | **caveat (D10)** | the `τ` chart is singular at `θ = π`. In 54 of 5 404 sampled harmonics `|p−q| ≤ 1e−9·scale`, and the two parametrizations then disagree about whether `π` is a root by a few `1e−3` |

### T4

| step | verdict | reason |
|---|---|---|
| T4.1a–c | AGREE | definitions; the interval bound `‖x_b − x_a‖²` is constant by T3.3, correctly justified |
| Lemma T4.2 | AGREE | proof re-checked line by line. `P ∩ Q ⊆ ∂P ∩ ∂Q` argument is right; transversal crossing forces an interior sector to be shared (four local sectors, each polygon occupies two, they must meet); the maximal-overlap-segment argument correctly produces a vertex on the other boundary because `e_P ∩ e_Q` is a segment whose endpoints are endpoints of `e_P` or `e_Q`. Correct |
| Cor. T4.2′ | AGREE | continuity + the lemma. Uses T2 to rule out any other way for the configuration to change — correct and worth the emphasis in T5.2d |
| Thm T4.2″ (T4.1) | **AGREE** | I checked the one gap the statement could have: an *isolated* overlap angle. Interior overlap is an **open** condition in `θ` (an open ball inside both interiors persists under a continuous motion), so the overlap set has no isolated points and `Θ_max = θ_{i*}` exactly. **Task 1(a): the graze treatment is right.** Independently reproduced with my own separating-axis overlap test (not `collision.hpp`): hexagons `θ₁ = 1.047198`, **zero interior overlap at 199 sample angles across `(0, 2.094395)`**, overlap immediately after. `min`-over-roots would report `π/3`, the true range is `2π/3` |
| T4.3-generic | AGREE with the Deriver | the codimension counting is right, and the proviso ("provided the polynomial is not identically zero on 𝕏") is exactly where it fails. Correctly labelled CONJECTURE and correctly said to be false on the papers' own symmetric figures |
| T4.3-wedge | AGREE | correctly labelled CONJECTURE |
| T4.4 | **AGREE with a correction (D3)** | the derivation of `β_e = 2π − α_f − α_g` from the fixed hinge point and rigid sectors is correct; `β_e ∈ 𝒞(X)` is correct. My independent recomputation of `β_e` (CCW interior angle, reflex corners included) matches `hinge_beta()` **exactly** (0.000e+00, 933 hinge edges, **T4-a**). But the `[N]` line "`Θ_max` equals `min_e β_e` to 8.88e−16 on all 8 split-free patterns" contradicts the Deriver's own table (triangles: `Θ_max = 3.141593`, `min β = 4.188790`). The true statement is `Θ_max = min(min_e β_e, π)` on the search range, verified to 8.9e−16 (**T4-b**) |
| T4.5a (T4.2)/(T4.3) | AGREE | `γ_f(θ) − γ_g(θ) = cA + sB`, so the squared distance is a harmonic — correct, and the pruning it induces is sound. Measured: keeps 42.6 % of pairs and changes `Θ_max` on 0 of 96 samples (**T4-f**) |
| T4.5b (T4.4) | **DISAGREE (D1)** | the inequality chain is true, but the application to K2c is wrong — see D1 |
| H-LOC and the `O(n)` argument | **DISAGREE (D2)** | the two statements called equivalent are not; one is trivial, the other is false |
| T4.H.1–5 | AGREE | |

### T5

| step | verdict | reason |
|---|---|---|
| T5.1 (T5.1) | AGREE | signed area is a quadratic form; `|F|` strict inequalities; basic open semialgebraic. Third-difference test confirms exact quadratic dependence on `t` to 3.4e−14 (**T5-a**) |
| T5.1 negative claim | AGREE, reproduced | **Task 1(f): T5.1 does NOT claim positive orientation suffices** — it says explicitly that it is necessary and not sufficient, and that `n_inv = 0` makes K1a over-count. Reproduced: **51 of 70** positively-oriented shape-space samples have `Θ_max = 0`; one saved to `derivations/check_failures/T5_1_positive_orientation_zero_range_hexagons.json`. See D11 for the certificate wording |
| T5.2a | AGREE | Tarski–Seidenberg applied to a first-order formula with polynomial atoms; the honest framing (the content is the explicit description, not semialgebraicity) is right |
| T5.2b | **DISAGREE (D4, D5)** | not proved to be an inner approximation, and one stated atom is not degree 4 |
| T5.2c | AGREE | "not by quadrics" is right (products and discriminants of the quadratics are quartic); "not basic" is correctly labelled an unproved negative. The replacement headline statement is defensible |
| T5.2d | AGREE | completeness, not the harmonic form, is the theorem — correct emphasis |
| T5.3 (T5.2)–(T5.5) | AGREE | re-derived: `h = det(R_f d, 2sJΔu) = 2s[c·det(d,JΔu) − σ_f s·det(d,Δu)] = sin θ·det(d,JΔu) − σ_f(1−cos θ)·det(d,Δu)`, so `p = −σ_f det(d,Δu)`, `q = −p`, `r = det(d,JΔu) = ⟨d,Δu⟩`. `(1+τ²)h = 2τ(r + pτ)`, second root `τ* = −r/p`. `h′(0) = r`, so Eq. (9) certifies `r > 0`, and `τ* > 0` iff `p < 0` iff `σ_f det(d,Δu) > 0`. All verified: `p` 3.0e−13, `q = −p` 3.0e−13, `r` 1.9e−14, `h(2 arctan(−r/p)) = 0` 5.4e−14, over **4 300** split-edge samples (**T5-b..e**). The algebraic precondition `r>0 ∧ p<0` fires on **871 of 4 300** samples (20 %) |
| T5.3 extra clauses | AGREE | the two extra clauses (interval test at `θ*`, wedge/entering test) are necessary, and the direction of the correction is right: K1c over-counts without them, so a PASS at ≥ 3 % is only valid after applying them |
| T5.H.1–4 | AGREE | |

### T6, T7

| step | verdict | reason |
|---|---|---|
| T6.1 | AGREE | softmin bound `min − log|𝒜|/β ≤ Θ_soft ≤ min` is standard |
| T6.2 (T6.2)/(T6.3) | AGREE | re-derived: `dh/dt_i = ∂_ip + cos θ ∂_iq + sin θ ∂_ir + h′(θ)·∂θ/∂t_i = 0`, `h′(θ) = −q sin θ + r cos θ`. In `τ`: `∂_τ g = 2r + 2(p−q)τ` and `dθ/dτ = 2/(1+τ²)`. Both correct. Verified against central differences of the root recomputed from scratch: worst relative error **3.9e−7** at `Δ = 1e−5` over **2 160** (pair, coordinate) samples on 12 graphs (**T6-a**) |
| T6.3, T6.4 | AGREE | the certificate list is right, and item 2 correctly requires both positive orientation and `Θ_max > 0` (see D11 for the third item I would add) |
| T6.H.1–4 | AGREE | T6-smooth correctly labelled CONJECTURE |
| T7 (i) | AGREE | `L = R D` recomputed from my own `R` and `D`: **exactly 0** difference on 16/16 (**T7-a**) |
| T7 (ii) | AGREE | `yᵀL = (Rᵀy)ᵀD`; `leftnull(D)` is the circulation space. `rank(L) = H − dim Z` on 16/16 |
| T7 (iii) (T7.1) | AGREE | out-harmonic residual 0 on every left-null vector (vacuous on bounded patches, `dim Z = 0`); non-vacuous on the tori: residual ≤ 1.1e−15, and the left-null vector is constant to 8.0e−16 (**T7-b** + torus block) |
| T7 (iv) | AGREE | `1ᵀL = 0` **exactly** and `rank(L) = H − 1` on 3/3 torus patches. The consequence — 2026 §4.4 is off by one for boundary-free patterns — is correct |
| T7 (v) | AGREE | corrected row-sum identity holds 16/16, naive identity holds **0/16**, and all 16 corpus cases have notch-owned hinge edges |
| T7 (vi) | AGREE | `H = |E_hinge| − |F| + c(Γ)` on 16/16 corpus **and 398/398 random-`σ`** cases |
| T7.H | **partly DISAGREE (D7)** | "if a user σ creates a split-cut cycle … all of T7 is then void" is too pessimistic: F11's partition property, the `H` formula and Step 6's row space all survived 398/398 random-`σ` cases, 210 of them with split-cut cycles |

---

## 2. Disagreements, with my own derivation

### D1 — T4.5b: the `√2` correction to K2c is wrong for the frame K2c actually uses

`core.md` T4.5b says `ideas/ranking.md` K2c's `ρ = max(‖C‖, ‖S‖)` "is not an upper bound",
that K2c's pruning "as written is *unsound*", and prescribes multiplying `ρ_f` by `√2`.

K2c's spec says: "per face compute `ρ_f = max over its vertices of max(‖x_u‖, ‖χ_u‖)`
**in the face's own frame**", and `contact.hpp::swept_discs` implements exactly that, with
`x = C_u − gc_f`, `χ = S_u − gs_f`. In that frame, using (T1.5),

```
   gc_f = x̄_f ,   gs_f = J(2u_f − σ_f x̄_f) ,
   C_u − gc_f = x_u − x̄_f ,
   S_u − gs_f = J(−σ_f x_u + σ_f x̄_f) = − σ_f J (x_u − x̄_f) .
```

The `2u_f` term **cancels identically**. The two columns are therefore orthogonal and of equal
norm, so

```
   max(‖x‖, ‖χ‖)  =  σ_max([x | χ])  =  ‖x_u − x̄_f‖   exactly.
```

So in the frame K2c and the code use, `ρ_max` is the **exact** swept radius, not a `√2`
under-estimate. Measured: `ρ_max = circumradius` to 1.4e−14 (**T4-c**), and
`max_θ ‖y_u(θ) − γ_f(θ)‖ = ‖x_u − x̄_f‖` to 1.9e−14 (**T4-d**), 15 132 copies. The Deriver's
D4 measurement (68 % violations, ratio `√2`) is about the **raw** frame with the origin at the
seed face — which no pruning uses. I reproduce that number too (61 % here, worst ratio 1.4142,
and `σ_max` is a sound bound to 7.1e−15, **T4-e**), but it has no consequence for K2c.

Two consequences the Deriver has backwards:

1. `contact.hpp`'s `rho` (the `√2` form, commented "sound bound") is the one that is loose by
   exactly `√2`; `rho_max` (the "spec's form") is tight. Applying the prescribed `√2` factor
   would make the pruning looser, not sound.
2. **K2c's gate test is vacuous.** `max_f ρ_f / r_f ≡ 1` by the algebra above, so the log–log
   fit against `log n` has slope 0 for every patch, at every size. Measured `1.000000` at
   `F = 12, 24, 44, 68, 96, 156`. K2c as specified **cannot fail and cannot test H-LOC**.

### D2 — H-LOC: the "equivalently" is false; the version the `O(n)` argument needs is refuted

T4.5 states H-LOC as

> `σ_max([C_u | S_u]) ≤ κ · r_{f(u)}` … where the frame is centred at the face centroid —
> equivalently, the path sum `u_f` of (T1.6) does not grow with the diameter of `Γ`.

These are not equivalent. By D1 the first is **trivially true with `κ = 1`** — it is T1.D
(faces are rigid) restated, and `u_f` has cancelled out of it. It carries no information.

The `O(n)` argument in the next paragraph needs something else, and says so: "if every face's
swept region is contained in a disc of radius `κ r_f` about a **fixed** point (its flat
centroid)". That requires `‖γ_f(θ) − x̄_f‖ = O(r_f)` uniformly, i.e. it *is* the statement
about `u_f`. Measured on a growing square patch:

| clip radius | `F` | patch diameter | `max_f ‖γ_f(θ) − c_f‖ / r_f` | `max_f ρ_f / r_f` |
|---|---|---|---|---|
| 1.6 | 12 | 3.16 | 3.60 | 1.000000 |
| 2.6 | 24 | 5.10 | 5.75 | 1.000000 |
| 3.6 | 44 | 7.07 | 8.05 | 1.000000 |
| 4.6 | 68 | 9.06 | 10.23 | 1.000000 |
| 5.6 | 96 | 11.05 | 12.40 | 1.000000 |
| 7.0 | 156 | 13.93 | 15.88 | 1.000000 |

The drift grows linearly with patch diameter, exactly as it must: deployment contracts the
patch by roughly `cos(θ/2)` about the seed, so a face at distance `d` from the seed moves by
`Θ(d)`. **H-LOC in the form the `O(n)` count needs is FALSE, not merely unproved.**

What survives: the `O(n)` conclusion may still hold by a different argument — the motion is
close to a global similarity, so *relative* centroid distances contract uniformly and the
packing count can be done in the deployed frame. That argument is not in `core.md`. Until it
is written, T4.5 ships as the **exact broad phase** `(T4.2)/(T4.3)` (which is sound with no
hypothesis, and which I verified preserves `Θ_max` on 96/96 samples) and the words "`O(n)`",
"certified active set" and "locality theorem" must not appear.

*What would refute my refutation:* an H-LOC statement in a frame that follows the deployment
(e.g. bounding `‖γ_f(θ) − γ_g(θ)‖` relative to `‖c_f − c_g‖`) rather than the flat frame.
That is the version I recommend the Deriver attempt; the data above is consistent with it.

### D3 — T4.4's `[N]` line contradicts the Deriver's own table

"On all 8 split-free patterns tested, `Θ_max` computed by (T4.1) equals `min_e β_e` to
8.88e−16" is inconsistent with the T4.2 table, which shows triangles at `Θ_max = 3.141593`
and `min β = 4.188790`. The correct statement, given deviation 9's search range `(0, π]`:

```
   E_split = ∅   ⟹   Θ_max = min( min_e β_e , π ) .
```

Verified to 8.9e−16 on the 4 split-free corpus patterns (**T4-b**); triangles is precisely the
case where `min β = 4π/3 > π`.

### D4 — T5.2b is not proved to be an inner approximation of `U(ε)`

T5.2b drops the interval and overlap tests and asks only that no orientation harmonic has a
root in `(0,T)`. That is a *stricter* root condition, so the region shrinks — but shrinking is
not the same as being contained in `U(ε)`. `Θ_max ≥ ε` requires the interiors to be disjoint
on `(0,ε)`; the root condition constrains only *changes* of overlap status, and says nothing
about the status just after `θ = 0`. A configuration already penetrating at `θ = 0⁺` has its
contact **at** `θ = 0`, invisible to any root scan on `(0,T)` — which is exactly why
`contact.hpp` carries a separate `penetrates_immediately()`, and exactly the failure mode
T5.1's own counterexamples exhibit.

The description therefore needs one more atom: **no face–face interior overlap at `θ = 0⁺`**
(a finite conjunction of semialgebraic disjointness conditions, so `U(ε)` stays semialgebraic).
Without it, `T5.1 ∧ (no root in (0,T)) ⊄ U(ε)` is unproved.

Search result (honest): 343 positively-oriented samples, 3 with no root in `(0, 0.02)`, 261
with `Θ_max < 0.02` — and **every** one of those 261 did have a root in `(0, 0.02)`. So I found
**no counterexample**; this is a proof gap, not a demonstrated falsehood.

### D5 — T5.2b's stated atom list is not degree ≤ 4, though a degree-≤ 4 list exists

The printed disjunct `g(τ_v)·g(0) > 0` with `τ_v = −r/(p−q)` clears to

```
   g(τ_v) = ( (p+q)(p−q) − r² ) / (p−q) = −(disc/4)/(p−q) ,
   g(τ_v)·g(0) > 0   ⟺   −(disc/4)·(p+q)·(p−q) > 0 ,
```

which is degree **8** in `t`, not 4. The headline claim survives, because an equivalent
degree-≤ 4 atom list exists. With `A = p−q`, `B = 2r`, `C = p+q` (all quadratic in `t`), "both
roots in `(0,T)`" is

```
   disc = B² − 4AC ≥ 0  ∧  A·g(0) > 0  ∧  A·g(T) > 0  ∧  A·B < 0  ∧  −AB − 2A²T < 0 ,
```

every atom of degree ≤ 4, and "no root in `(0,T)`" is `[g(0)g(T) > 0] ∧ ¬[both roots in
`(0,T)`]`. Replace the printed combination with this one.

### D6 — T1.H.2: split-cut cycle ⟹ `Γ` disconnected, but not conversely

Measured over 398 random-`σ` pairs: 210 with a split-cut cycle, **all 210** with `c(Γ) > 1`;
but **124** pairs had `c(Γ) > 1` with a split *forest*, and 0 had a split cycle with `Γ`
connected. So "iff" is wrong in one direction. This matters because T1.H.2 uses the
equivalence to fold case 2 into case 1; the fold is valid, but the converse claim is not.

### D7 — T1 Step 6 is far more robust than the derivation's own risk assessment

The Deriver ranks T1 Step 6 among "the three weakest steps" and proposes the attack: find an
`(M, σ)` where the face-cycle basis of `Γ` and `holes_partition` disagree, e.g. with a
split-cut cycle. **Task 1(c): I ran that attack and it fails.**

I built the closure rows symbolically — BFS `potvec[g] = potvec[f] + σ_g·e_{src(e)}`, one row
`potvec[g] − potvec[f] − σ_g e_{src(e)}` per non-tree hinge edge — giving a `b₁(Γ) × N` matrix
`G`, and compared its **row space** with `L` from `holes_partition` via
`rank(G) = rank(L) = rank([G;L])`. Result: identical on **16/16** corpus cases and **398/398**
random-`σ` cases, with `b₁(Γ) = H` on all of them. Among those 398: 334 with `c(Γ) > 1`, 334
with `M′` disconnected, 210 with split-cut cycles. `holes_seed_growing` and `holes_partition`
agreed 398/398, and the preimages partitioned `E_hinge ⊔ E_split` 398/398 (F11 holds under
fully random `σ`).

**What the theorem then asserts with a split-cut cycle:** `Γ` is disconnected, so the potential
`u` is determined only per component and `M′` falls apart into `c(Γ)` independently placed
pieces (T1.H.1). Eq. (2) and `L` remain correct — they are still exactly the closure conditions
on `Γ`'s cycle space, which is what my row-space test shows. The `2(c(Γ)−1)`-dimensional
relative placement of the components is not in `𝕏` and is not a mechanism of an assembled
structure; the deployment theory applies per component. I verified this directly: on the
projected `X₀` for random `σ`, (T1.5) built from a **different** spanning forest reproduces
`deploy()` to 7.5e−15 on 398/398 (**CE-d**), and `deploy()`'s `max_mismatch` is 0 on all of
them. So T7.H's "all of T7 is then void" should be softened to "the `H`-row interpretation as
*geometric* holes fails (F16), while `L = RD`, the rank identity and the `H` formula survive".

### D8 — `contact.hpp`'s flat-centroid pruning is unsound in principle (the Deriver's (T4.3) is not)

`core.md` T4.5a's static test (T4.3) minimises the *harmonic* `‖γ_f(θ) − γ_g(θ)‖²`, which is
sound. `contact.hpp::candidate_pairs(..., use_static = true)` instead uses the **flat**
centroids `‖c_f − c_g‖`. Because deployment contracts centroid distances while face radii stay
fixed (T1.D), that test can discard a pair whose faces actually approach. Measured: it discards
**19 842** pairs that the sound moving test keeps, with
`max (‖c_f − c_g‖ − min_θ‖γ_f(θ) − γ_g(θ)‖) = 6.59`. It nevertheless returned the correct
`Θ_max` on 96/96 samples (**T4-g**), so this is a soundness defect without an observed
counterexample. Use the moving form (**T4-f**, 0/96 wrong, keeps 42.6 % of pairs).

### D9 — documentation inconsistency in `deploy_basis.hpp`

The header comment says `S_pv = -sigma_f J x_v + u_f`. That `u_f` means `2 J u_f` in `core.md`'s
notation (the header's `u_f` is the whole translation direction, `t_f = sin(θ/2)·u_f`). The code
is right — `S` is computed as `2 dY/dθ|₀`, which equals `J(2u_f − σ_f x_v)` (verified to
7.3e−15, **T1-b**) — but the two files use `u_f` for different objects. Fix one of them before
either becomes a paper formula.

### D10 — the `τ` chart is singular at `θ = π`

(T3.6)'s leading-coefficient branch triggers when `p − q = 0`; near that, the reported root
jumps between "exactly `π`" and "a few `1e−3` past `π`". 54 of 5 404 sampled harmonics were in
the regime `|p−q| ≤ 1e−9·scale`. Since the search range ends at `π`, this decides whether a
contact at the very end of the range is counted. Not an error in T3.4, but a caveat: near
`θ = π` use the amplitude/phase form, not the `τ`-quadratic.

### D11 — the validity certificate: what T5.1 says, and what I would add

**Task 1(f).** T5.1 does **not** claim positive orientation suffices — it states the opposite
explicitly, calls `n_inv = 0` "necessary but not sufficient", and says K1a's `p_valid`
over-counts. I agree and reproduced the counterexamples (51 of 70). The Deriver's proposed
certificate is `Θ_max(X) > 0` by T4.2″.

That is the right *deployment* certificate but it is not the whole story, and the two conditions
are different:

* **Flat-pattern validity (fabrication):** all `A_f > 0` **and** no face–face interior overlap
  at `θ = 0`. This certifies that `M(X)` is an embedded planar graph — a real sheet to cut.
  Positive orientation alone gives only local injectivity.
* **Deployment validity:** `Θ_max(X) > 0`. This is strictly different: a flat state with two
  faces overlapping at `θ = 0` could in principle separate as the cuts open, and conversely.

T6.4's item 2 lists only the second. I recommend the certificate carry **three** items:
`min_f A_f > 0`, no face–face overlap at `θ = 0`, and `Θ_max > 0` by T4.2″ cross-checked
against the bisection. The `θ = 0` overlap test is also the atom missing from T5.2b (D4), so
adding it closes both gaps at once.

---

## 3. Test results

All from `code/tests/derivation_tests.cpp`; 25 test cases, 29 033 assertions, 0 failures.

| id | identity / claim | samples | max error | tolerance | verdict |
|---|---|---|---|---|---|
| T1-a | `Y(θ) = cos(θ/2)C + sin(θ/2)S` vs `deploy()` | 85 806 | 1.42e−14 | 1e−11 | PASS |
| T1-b | my `(C,S)` from (T1.3)/(T1.5) == `deploy_basis()` | 12 712 | 7.32e−15 | 1e−11 | PASS |
| T1-c | `X ∈ 𝕏` ⟹ closure residual 0 on every non-tree edge | 192 | 1.28e−14 | 1e−10 | PASS |
| T1-d | `max_mismatch = 2|sin(θ/2)|·closure defect` (rel) | 192 | 1.90e−15 | 1e−9 | PASS |
| T1-A | `yᵀ(AAᵀ)^{-1}y = 1` (centred ellipse) | 46 890 | 1.98e−09 | 1e−8 | PASS |
| T1-B1 | split duplicates are the same vector | 5 160 | 2.25e−14 | 1e−11 | PASS |
| T1-B2 | offset `= 2 sin(θ/2) J Δu` (T1.8) | 5 160 | 1.87e−14 | 1e−11 | PASS |
| T1-C | `|signed hinge angle| = θ` | 27 990 | 1.48e−13 | 1e−11 | PASS |
| T1-D | face signed area `θ`-independent | 19 800 | 4.97e−14 | 1e−11 | PASS |
| T2-a | pin residual (T2.1) with `ω = −σ/2`, `w` of (T2.3) | 22 392 | 1.26e−14 | 1e−10 | PASS |
| T2-b | `‖A(Y_θ)σ‖∞ / ‖A‖∞` | 480 | 1.02e−14 | 1e−12 | PASS |
| T2-c | pencil `A(θ) = cos(θ/2)A_c + sin(θ/2)A_s` | 480 | 4.00e−15 | 1e−12 | PASS |
| T3-a | (T3.2) det harmonic (relative) | 12 345 | 3.55e−15 | 1e−12 | PASS |
| T3-b | (T3.3) dot harmonic (relative) | 12 345 | 1.55e−15 | 1e−12 | PASS |
| T3-c | (T3.4) squared-length harmonic (relative) | 12 345 | 6.10e−16 | 1e−12 | PASS |
| T3-d | my `(p,q,r)` == `orient_harmonic()` | 12 345 | 4.53e−14 | 1e−10 | PASS |
| T3-e | intra-face `q = 0` | 4 808 | 1.74e−14 | 1e−12 | PASS |
| T3-f | intra-face `r = 0` | 4 808 | 1.78e−14 | 1e−12 | PASS |
| T3-g | intra-face `p =` flat signed area ×2 | 4 808 | 1.74e−14 | 1e−12 | PASS |
| T3-h | intra-face `q″ = r″ = 0` | 4 808 | 1.56e−14 | 1e−12 | PASS |
| T3-i | every `τ`-quadratic root satisfies `h(θ) = 0` (rel) | 5 404 | 8.40e−14 | 1e−9 | PASS |
| T3-j | my roots == `harmonic_roots` on `(0, π−1e−3]` | 5 404 | 3.99e−09 | 1e−7 | PASS |
| T4-a | `β_e = 2π − α_f − α_g` recomputed by hand | 933 | 0.00e+00 | 1e−12 | PASS |
| T4-b | split-free `Θ_max = min(min_e β_e, π)` | 4 | 8.88e−16 | 1e−6 | PASS |
| T4-c | `ρ_max` (face frame) `=` circumradius, exactly | 15 132 | 1.40e−14 | 1e−12 | PASS |
| T4-d | `max_θ‖y_u − γ_f‖ = ‖x_u − x̄_f‖`, exactly | 15 132 | 1.87e−14 | 1e−9 | PASS |
| T4-e | `σ_max([C|S])` is a sound bound (raw frame) | 15 132 | 7.11e−15 | 1e−9 | PASS |
| T4-f | moving-centroid pruning preserves `Θ_max` | 96 | 0.00e+00 | 1e−6 | PASS |
| T4-g | flat-centroid pruning preserves `Θ_max` | 96 | 0.00e+00 | 1e−6 | PASS (see D8) |
| T5-a | `A_f(X₀+Φt)` exactly quadratic in `t` (3rd diff.) | 6 312 | 3.40e−14 | 1e−9 | PASS |
| T5-b | (T5.2) `p = −σ_f det(d, Δu)` | 4 300 | 3.01e−13 | 1e−11 | PASS |
| T5-c | (T5.2) `q = −p`, i.e. `p + q = 0` | 4 300 | 3.01e−13 | 1e−11 | PASS |
| T5-d | (T5.2) `r = ⟨d, Δu⟩ = det(d, JΔu)` | 4 300 | 1.86e−14 | 1e−11 | PASS |
| T5-e | (T5.4) `h(2 arctan(−r/p)) = 0` | 4 300 | 5.40e−14 | 1e−9 | PASS |
| T6-a | (T6.2) `∂θ/∂t_i` vs central difference (relative) | 2 160 | 3.90e−07 | 1e−6 | PASS |
| T7-a | `L = R D` exactly | 16 | 0.00e+00 | 0 | PASS |
| T7-b | every left-null vector of `L` is out-harmonic | 16 | 0.00e+00 | 1e−9 | PASS |
| CE-a | F11 partition under **random** `σ` | 398 | 0 failures | 0 | PASS |
| CE-b | `H = |E_h| − |F| + c(Γ)` under random `σ` | 398 | 0 failures | 0 | PASS |
| CE-c | T1 Step 6 row space under random `σ` | 398 | 0 failures | 0 | PASS |
| CE-d | (T1.5) tree-independent on `X₀`, `c(Γ) > 1` | 398 | 7.54e−15 | 1e−9 | PASS |

Controls that must fail, and do: the sign-flipped `S` misses `deploy()` by 45 (O(1)); off the
shape space the closure defect is ≥ 0.41 on 192/192 samples; the naive row-sum identity holds
0/16.

---

## 4. Counterexample searches

**(a) Grazes on symmetric tilings — the T4.2″ attack (task 1a).** Independent separating-axis
overlap test on all convex-face patterns with split cuts:

| pattern | `|𝒞|` | `θ₁` (first contact) | `Θ_max` (T4.2″) | interior overlap on `(0, Θ_max)` | overlap just after |
|---|---|---|---|---|---|
| hexagons | 3 | 1.047198 | **2.094395** | none at 199/199 sampled angles | yes |
| snub_square | 11 | 0.834111 | 0.834111 | none at 199/199 | yes |
| truncated_square | 4 | 2.356194 | 2.356194 | none at 199/199 | yes |

The hexagon graze is real and `min`-over-roots is wrong there by a factor of 2. The `4.8.8`
arithmetic of T4.Check checks out: `2π − 3π/4 − π/2 = 3π/4 = 2.356194`. **T4.2″'s treatment is
correct**, and I additionally closed the one gap the statement could have had: interior overlap
is an open condition in `θ`, so there are no isolated overlap angles and the interval scan
cannot skip one.

**(b) Random `σ`, including connectivity-violating ones (tasks 1c).** 398 usable `(graph, σ)`
pairs over 6 generators, `σ` i.i.d. uniform:

| observation | count |
|---|---|
| split-cut cycle present | 210 |
| `Γ` disconnected (`c(Γ) > 1`) | 334 |
| `M′` disconnected | 334 |
| split cycle **but** `Γ` connected | **0** |
| `Γ` disconnected **but** split forest | **124** (refutes the "iff" of T1.H.2 — D6) |
| Alg. 1 (seed growing) ≠ partition formulation | 0 |
| preimages fail to partition `E_hinge ⊔ E_split` | 0 |
| `H ≠ |E_h| − |F| + c(Γ)` | 0 |
| T1 Step 6 row-space mismatch | 0 |
| `X₀` with `deploy()` mismatch > 1e−9 | 0 |

**(c) Boundary-touching split components.** All 16 corpus cases have notch-owned hinge edges.
The corrected row-sum identity holds 16/16, the naive one 0/16 — T7(v) confirmed.

**(d) Positive orientation without validity.** 51 of 70 positively-oriented shape-space samples
have `Θ_max = 0`; a witness is saved at
`derivations/check_failures/T5_1_positive_orientation_zero_range_hexagons.json`.

**(e) H-LOC.** See D2 — refuted in the form the `O(n)` argument needs, by the growing-patch
table. No failing sample needed to be saved; the effect is systematic.

**(f) T5.2b's missing `θ = 0` atom.** 343 positively-oriented samples at `ε = 0.02`; 261 with
`Θ_max < ε`, all of which had a root in `(0,ε)`. **No counterexample found** — the gap is a
proof gap (D4), not a demonstrated falsehood.

`derivations/check_failures/` contains the one saved witness; nothing else failed.

---

## 5. What the paper may claim, after this check

1. T1, T2, T3, T4.2/T4.2′/T4.2″, T4.4 (with the `π` cap), T4.5a, T5.1, T5.2a, T5.3, T6.2 and
   T7 are theorems and can be stated as such.
2. T4.5b's `√2` correction must be **withdrawn** for the face frame; K2c's `ρ_f` is exact, and
   K2c must be redesigned because as specified it is a tautology.
3. Nothing may be labelled `O(n)`, "locality theorem" or "certified active set". T4.5 ships as
   the exact broad phase only.
4. T5.2b must gain the `θ = 0⁺` embeddedness atom and the degree-≤ 4 atom list of D5 before
   "quantifier-free description of degree ≤ 4" is claimed.
5. The validity certificate is three items, not two (D11).

---

# Round 2 — check of the revised `core.md` (fresh context, changed material only)

Scope: only what round 2 changed or added. Inputs read: `derivations/core.md` (revised),
this file's round-1 verdicts, `derivations/scratch/check_r2.cpp`, `code/README.md`.
New program work: two test cases appended to `code/tests/derivation_tests.cpp` (that file only;
`CMakeLists.txt` and every other source untouched). Rebuilt with the command in its header:

```
clang++ -std=c++20 -O2 -arch arm64 -I/opt/homebrew/include \
  -I/opt/homebrew/include/eigen3 -Icode/src \
  code/tests/derivation_tests.cpp code/build/libkiri_core.a \
  -o code/build/derivation_tests && ./code/build/derivation_tests
```

Result: **27 test cases, 29 054 assertions, 0 failures.**

## R2.1 Verdicts on the changed material

| item | verdict | reason |
|---|---|---|
| **T4.5b′** exact swept radius | **AGREE** | algebra redone from scratch below; the lemma is not merely true, it is *stronger* than stated |
| **H-LOC refuted**, `O(n)` wording removed | **AGREE** | refutation is sound; grep residue is all negation or explicitly-labelled heuristic — list in R2.3 |
| **D3** `Θ_max = min(min β, π)` | **AGREE on the correction**, one label missing | the `≥` direction is measured, not derived, and carries no `[A]`/`[N]` tag — see R2.4(a) |
| **T5.2b′** `POS ∧ EMB ∧ NOROOT ⟹ Θ_max ≥ ε` | **AGREE** | the openness argument and the compactness argument are both correct as written; `EMB` is the `θ = 0⁺` predicate, which is what the proof needs |
| **T5.1d** `τ = 0` deflation | **AGREE on the algebra**, **DISAGREE with the "0 at every `ε`" claim** | the deflated atom list is exact — but only after a **third** structural class is excluded, which `core.md` does not name. See R2.5 |
| **D6** split cycle ⟹ `c(Γ)>1` only | **AGREE** | one-directional statement is what round 1 established; 124/398 forest counterexamples to the converse stand |
| **D7** Step 6 demoted from weak-step list | **AGREE** | matches this file's CE-c; row-space equality 398/398 including 210 split-cut-cycle cases |
| **D8** flat-centroid pruning as a note | **AGREE** | correctly recorded as a soundness defect with no observed wrong `Θ_max` |
| **D9** `u_f` vs `ν_f` note (§0.10) | **AGREE** | `t_f = 2s·J u_f` and `t_f = s·ν_f` give `ν_f = 2J u_f`, `u_f = ½Jᵀν_f` since `J⁻¹ = Jᵀ`. Correct |
| **D10** `τ` chart singular at `θ = π` | **AGREE** | `h(π) = p − q`, so `θ = π` is a root iff the leading coefficient `A = p − q` vanishes; the `(T3.6)` branch and the amplitude/phase fallback are stated correctly |
| **D11** three validity predicates / five certificate items | **AGREE**, one wording nit | see R2.4(b) |

## R2.2 The two algebras, redone independently

**T4.5b′.** From (T1.5), `C_u = x_u` and `S_u = J(2u_f − σ_f x_u)`. The deployed face centroid is
the image of the flat centroid under the same isometry, so `γ_f(θ) = c·x̄_f + s·J(2u_f − σ_f x̄_f)`,
i.e. `gc_f = x̄_f`, `gs_f = J(2u_f − σ_f x̄_f)`. Subtracting, the `2u_f` term cancels:

```
   x   := C_u − gc_f = x_u − x̄_f ,
   χ   := S_u − gs_f = J( −σ_f x_u + σ_f x̄_f ) = −σ_f · J ( x_u − x̄_f ) = −σ_f J x .
```

Hence `⟨x, χ⟩ = −σ_f⟨x, Jx⟩ = 0` and `‖χ‖ = ‖x‖`, so `[x | χ]` has both singular values equal to
`‖x‖` and `σ_max = max(‖x‖,‖χ‖) = ‖x_u − x̄_f‖`. Stronger than `core.md` states: since the columns
are orthogonal and equinormal,

```
   ‖ y_u(θ) − γ_f(θ) ‖²  =  c²‖x‖² + s²‖χ‖² + 2cs⟨x,χ⟩  =  ‖x‖²   for EVERY θ ,
```

so the distance is **constant in `θ`**, not merely maximised at `‖x_u − x̄_f‖`. That is T1.D in
coordinates, exactly as `core.md` says. The `√2` withdrawal is correct and the K2c gate ratio is
identically 1.

**T5.1d.** With `C := g(0) = 0`, `g(τ) = Cτ⁰ + Bτ + Aτ² = τ(Aτ + B)`. The roots are `τ = 0` and, if
`A ≠ 0`, `τ = −B/A`. Multiplying by `A² > 0`:

```
   −B/A > 0   ⟺   A·B < 0 ,
   −B/A < T   ⟺   −A·B − A²T < 0 .
```

So "a root lies in `(0,T)`" is exactly `A·B < 0 ∧ −A·B − A²T < 0`, and (T5.1d) is its negation.
Both atoms are degree ≤ 4 in `t`. The factor is `A²T`, **not** `2A²T` as in (T5.1b) — (T5.1b) locates
the *vertex* `−B/(2A)`, (T5.1d) the *root* `−B/A`; `core.md` has both factors right.
The `A = 0` case is handled correctly by accident and by design: `A·B = 0` is not `< 0`, so (T5.1d)
returns "no root", which is right because `g = Bτ` has only the root at `0`.

## R2.3 Grep residue for `O(n)` / active-set / locality wording

Every surviving occurrence was read in context. **No residue that asserts the removed claims.**

| line | context | status |
|---|---|---|
| 31, 787, 789–818, 852, 1103, 1192–1193, 1400–1406, 1431–1433, 1461 | the refutation itself, the change log, the dead-end list, the "not to be confused with H-LOC" note | negation — correct |
| 817 | "The candidate list is `O(n²)`" | true and is the claim that replaces `O(n)` |
| 1080, 1295–1296 | "certified-feasible region", "certified valid, with a certified range" | different sense (a certificate at the returned `t̂`), not an active-set claim |
| 1225 | `𝒜 ⊆ 𝒞` "the current **active set**" in (T6.1) | a defined symbol in the optimizer, immediately governed by T6.3 |
| 1257–1261 | T6.3 "**Active set — a heuristic, with no size guarantee**", "`|𝒜|` has no proved `O(n)` bound … never as a 'certified active set'" | correct and explicit |

## R2.4 Two residues that need a label, not a correction

**(a) (T4.1b) has no proof of the `≥` direction.** `E_split = ∅ ⟹ Θ_max = min(min_e β_e, π)` is
printed as "the true statement is", with `[N]` support on 8 split-free patterns (4 in my corpus,
test **T4-b**, 8.88e−16). T4.4 derives `β_e ∈ 𝒞` and that `β_e` is a genuine overlap when the two
edge lengths differ — that gives `Θ_max ≤ min(min β, π)`. It does **not** show that no *earlier*
candidate root binds, which is the `≥` direction. On split-free patterns that is plausible and is
measured, but it is an unproved equality stated without a tag. Label it `[N]`-supported, or prove it.

**(b) `EMB` is a `θ = 0⁺` predicate; T6.4 item 2 reports a `θ = 0` predicate.** T5.2b.0(ii) defines
`EMB` as "pairwise disjoint interiors for all sufficiently small `θ > 0`" and identifies it with
`¬penetrates_immediately()` — which is what the proof of T5.2b′ uses to get `θ* ≥ θ̄ > 0`. But
T5.1 (D11) and T6.4 item 2 both say "no face–face interior overlap **at `θ = 0`**", i.e. flat
embeddedness of `M`. Those are not the same predicate, and flat embeddedness does not immediately
imply the `0⁺` version. The certificate as printed therefore does not literally supply `EMB`.
Fix is one word: item 2 should name `penetrates_immediately()`, or state both.

## R2.5 The one substantive disagreement: a third structural class at `τ = 0`

`core.md` reports (R2-E) **0 mismatches out of 2 363 380 at every `ε`** for (T5.1c)+(T5.1d).
I cannot reproduce that as an unqualified statement. Rebuilding the comparison from my own
`(p,q,r)` and my own `τ`-quadratic root finder, over **1.82 M** candidate harmonics on 16 graphs at
three `ε`, the deflated list mismatches direct root finding on **143** harmonics at every `ε`
(69 spurious roots, 74 missed) — and the count is **`ε`-independent**, which is the tell.

Every one of the 143 has the same shape, e.g. `p = 0.4330127, q = −0.4330127, r = 0`:

```
   C = p + q = 0   AND   B = 2r = 0 ,   A = p − q ≠ 0     ⟹    g(τ) = A τ² ,
   equivalently   h(θ) = p (1 − cos θ) ,   a DOUBLE root at θ = 0 and no other root.
```

This is a **third** structural class, distinct from the two `core.md` names (permanent incidences
and split-edge pairs, both of which have `C = 0` but `B ≠ 0`). Deflating once leaves `Aτ + B` with
`B` at round-off, so the sign of `−B/A` — and hence membership in `(0,T)` — is decided by noise;
direct root finding is equally ill-posed there, reporting a spurious root at `θ ≈ 2.3e−8`. The
clause (T5.1d) is not *wrong* on this class, it is *undecidable* at double precision.

Excluding that class explicitly (`|C| ≤ 1e−11·scale ∧ |B| ≤ 1e−11·scale`, **2 176** of 1.82 M
harmonics, 0.12 %), the deflated list is exact:

| `ε` | harmonics tested | double-root class skipped | (T5.1c) as printed | of those, `g(0)=0` | (T5.1c)+(T5.1d) |
|---|---|---|---|---|---|
| 0.200 | 1 817 732 | 2 176 | **61 201** (all spurious) | 61 201 | **0** |
| 0.020 | 1 817 748 | 2 176 | **62 310** (all spurious) | 62 310 | **0** |
| 0.001 | 1 817 705 | 2 176 | **62 327** (all spurious) | 62 327 | **0** |

So `core.md`'s conclusion survives — the Checker's degree-≤ 4 list is right, round 1's was wrong,
and deflation fixes it — but the claim must read *"0 mismatches once the double-root class
`g(0) = g′(0) = 0` is removed by the identity test of T3.H.1"*, not "0 mismatches". The deflated
atom list needs either a second deflation (`g = Aτ²` has no root in the open interval, so the
correct atom is simply "true") or an explicit membership rule putting this class with the
identically-zero harmonics of T3.H.1. It is decided by the combinatorics like the other two classes,
so this costs nothing in the quantifier-free description.

## R2.6 New tests (appended to `code/tests/derivation_tests.cpp`)

`TEST_CASE("R2 T4.5b' exact swept radius in the face frame (randomized)")` — 16 corpus cases ×
8 shape-space samples (`X₀` plus Gaussian null-space perturbations, σ = 0.25), **20 176 `M′`-copies**,
everything rebuilt from my own BFS potential and my own flat centroid:

| id | claim | samples | max error | tol |
|---|---|---|---|---|
| R2-a | `C_u − gc_f = x_u − x̄_f` (the `2u_f` term cancels) | 20 176 | 0.000e+00 | 1e−12 |
| R2-b | `S_u − gs_f = −σ_f J (x_u − x̄_f)` | 20 176 | 2.234e−14 | 1e−12 |
| R2-c | `⟨x,χ⟩ = 0` and `‖x‖ = ‖χ‖` | 20 176 | 2.063e−14 | 1e−12 |
| R2-d | `σ_max([x|χ]) = max(‖x‖,‖χ‖)` — exact, no `√2` | 20 176 | 5.440e−15 | 1e−12 |
| R2-e | `‖y_u(θ) − γ_f(θ)‖` is **constant** in `θ` | 20 176 | 2.321e−14 | 1e−12 |
| R2-f | that constant `= ‖x_u − x̄_f‖` | 20 176 | 2.149e−14 | 1e−12 |
| R2-g | my `(gc,gs)` `=` `swept_discs()` (library cross-check) | 5 280 | 1.014e−14 | 1e−10 |

`TEST_CASE("R2 tau=0 deflation (T5.1d): spurious roots with and without it")` — table in R2.5;
id **R2-h**, 1.82 M harmonics per `ε`, 0 mismatches, tolerance 0. Truth is `my_roots()`, the
round-1 independent `τ`-quadratic solver, which never reports the endpoint `0`; harmonics whose
root falls within `1e−9` of either endpoint are counted as ambiguous and skipped (≈ 9 500 per `ε`).

## R2.7 Final verdict

**There is ONE unresolved disagreement, plus two labelling residues.** Nothing else remains:
D1, D2, D4, D6, D7, D8, D9, D10, D11 are fully resolved, and the T4.5b′ and T5.2b′ proofs are
correct as written.

1. **Unresolved (R2.5).** `core.md`'s R2-E line — "**0 / 2 363 380**, at every `ε`" — is not
   reproducible as stated. A third systematically degenerate class, `g(0) = g′(0) = 0`
   (`p = −q`, `r = 0`; `h = p(1 − cos θ)`, a double root at `θ = 0`), is undecidable at double
   precision by both the atom list and direct root finding: 143 of 1.82 M harmonics, `ε`-independent.
   After excluding it, 0 / 1 817 732. **Fix:** name the class in T5.2b.2 and route it to T3.H.1's
   identity test (or add the trivial "no root" atom for it), and qualify the R2-E number.
2. **Labelling residue (R2.4a).** (T4.1b)'s `≥` direction is measured on 8 patterns, not derived,
   and is printed untagged as "the true statement is".
3. **Labelling residue (R2.4b).** `EMB` (T5.2b.0(ii), a `θ = 0⁺` predicate) and T6.4 item 2
   ("no overlap at `θ = 0`") are different predicates; the certificate as printed does not supply
   the hypothesis T5.2b′ needs.

None of the three changes a theorem's truth value; (1) changes a printed number and adds one
required clause, (2) and (3) are one-line edits.

---

# Round 3 (final) — confirmation of the three round-3 resolutions

Fresh context. Inputs read: `core.md` round-3 change log and §T3.H.5, T4.4, T5.1, T5.2b.0,
T5.2b′, T5.2b.1, T5.2b.2, T6.4; `derivations/scratch/check_r3.cpp`. New tests appended to
`code/tests/derivation_tests.cpp` (that file only): `R3 Lemma T5.1e and the three-class atom
list (randomized, chart-free truth)` and `R3 certificate POS ^ NOOVERLAP(eps/2) ^ NOROOT
implies Theta_max >= eps (randomized)`. Both rebuilt with the header build line and run;
**2 test cases, 59 982 assertions, 0 failures**.

## R3.1 — Lemma T5.1e / T5.1e′ and the three-class atom list: **AGREE**, with one missing side-condition

**Re-derived independently.** With `τ = tan(θ/2)`, `cos θ = (1−τ²)/(1+τ²)`, `sin θ = 2τ/(1+τ²)`:
`(1+τ²) h(θ) = (p+q) + 2rτ + (p−q)τ² =: g(τ)`. Since `1+τ² > 0`, the roots of `h` in `(0,π)`
correspond exactly to the roots of `g` in `(0,∞)`, and `(0,ε) ↔ (0,T)`. Verified numerically,
**200 000 random `(p,q,r,θ)`, max relative error 1.444e−13** (id **R3-a**).

The order of vanishing of `g` at `τ = 0` is fixed by `C = g(0) = p+q` and `B = g′(0) = 2r`, so
the split `C ≠ 0` / `C = 0, B ≠ 0` / `C = B = 0` is **exhaustive and mutually exclusive by
construction** — there is no fourth class. Class 3 forces `q = −p`, `r = 0`, hence
`h(θ) = p − p cos θ = p(1 − cos θ)`, and `A = p − q = 2p`, i.e. `p = A/2` as `core.md` states.
Verified: the identity holds to **4.44e−16** over 19 989 random class-3 harmonics (**R3-b1**).

**Lemma T5.1e is correct and its proof is correct.** `1 − cos θ > 0` on `(0, 2π)` because
`cos θ = 1` only on `2πZ`; one line, no gap. Confirmed by direct dense sampling rather than by
the argument: **19 989 class-3 harmonics × 3 999 angles in `(0,π)`, 0 sign changes, 0 exact
zeros, and `h/p > 0` throughout** (**R3-b2**).

I also re-derived the other two branches and they are right. (T5.1d): with `C = 0`, `g = τ(Aτ+B)`,
so a root lies in `(0,T)` iff `0 < −B/A < T`; multiplying by `A² > 0` gives exactly
`A·B < 0 ∧ −A·B − A²T < 0`, and if `A = 0` the atoms correctly report no root. (T5.1b)/(T5.1c):
`BOTH` is the standard "two real roots, both in `(0,T)`" condition — `disc ≥ 0`, `A·g(0) > 0`,
`A·g(T) > 0`, vertex `−B/(2A) ∈ (0,T)` cleared by `2A² > 0` — and it correctly catches the
interior tangency `disc = 0` that `g(0)g(T) > 0` alone would miss.

**`check_r3.cpp` rebuilt and rerun.** Its numbers reproduce exactly, at all three `ε`:

| `ε` | decided | ambiguous | class 1 / 2 / 3 | round-2 list | round-3 list |
|---|---|---|---|---|---|
| 0.200 | 2 145 387 | 47 | 2 070 850 / 73 446 / 1 138 | **153** (all class 3) | **0** |
| 0.020 | 2 145 387 | 47 | 2 070 850 / 73 446 / 1 138 | **153** (all class 3) | **0** |
| 0.001 | 2 145 387 | 47 | 2 070 850 / 73 446 / 1 138 | **153** (all class 3) | **0** |

**Its reference is NOT circular.** `truth_crossing()` never forms `τ`, never forms `g`, and never
solves a quadratic: it evaluates `h(θ) = p + q cos θ + r sin θ` at the break points
`{0⁺} ∪ {critical points in (0,ε)} ∪ {ε}` and looks for a sign change. The only place it touches
the class split is the sign at `0⁺`, taken as the first non-vanishing Taylor coefficient
(`C`, else `B`, else `A`) — and that is an independently correct fact about `h`, not an import
from the atom list: `h(0) = C`, `h′(0) = r = B/2`, `h″(0) = −q = A/2` when `C = B = 0`. Crucially,
if Lemma T5.1e were false the oracle would still see the sign change at a break point, so it is
**not blind to class 3** in the way `check_r2.cpp` was. R2.5's circularity objection is answered.

**Independently reconfirmed with a second, different oracle.** Test **R3-c** decides the truth
through the *amplitude/phase* form `h = p + Ρcos(θ−φ)` in `long double` — no `τ` chart, no
quadratic, no Taylor rule, no class split anywhere — by a 4 001-point dense sign scan on
`[θ_lo, ε)` with `θ_lo = min(1e−5, ε/10)`, plus a closed-form check that no root hides in the
un-scanned sliver `(0, θ_lo]` (those harmonics are reported ambiguous, not decided). On my own
corpus (16 cases × 4 shape-space samples):

| `ε` | decided | ambiguous | class 1 / 2 / 3 | three-class list | two-class list |
|---|---|---|---|---|---|
| 0.200 | 1 829 226 | 142 | 1 745 368 / 81 500 / 2 500 | **0** | **393**, all class 3 |
| 0.020 | 1 829 226 | 142 | 1 745 368 / 81 500 / 2 500 | **0** | **393**, all class 3 |

So on a different corpus, with a different oracle, the conclusion is the same: the class-3 clause
is **necessary** (393 failures without it, all in class 3, `ε`-independent) and **sufficient**
(0 failures with it). No mismatch of the conservative kind (atoms say root, truth says no
crossing) occurred either — 0 of 1 829 226 — so the asymmetry `core.md` flags is not exercised here.

**The one missing side-condition (new, R3.1a).** (T5.1e′) is printed as "`true` when `C = B = 0`",
unconditionally, but Lemma T5.1e is stated and proved only for `p ≠ 0`, i.e. `A ≠ 0`. When
`A = C = B = 0` the harmonic is identically zero and **every** `θ` is a root, so the atom `true`
is wrong there. This is not hypothetical: these pairs are systematic and numerous. Measured over
the same corpus, **14 928** candidate harmonics are identically zero, and the numeric class rule
puts **all 14 928 of them in class 3** (**R3-c** printout). On the `squares` case I checked what
they are: all **144** of them have `w` coincident with `a` or with `b` at every `θ` — permanent
hinge incidences (T3.H.2), 144 / 144, none of any other kind.

Both `check_r3.cpp` and my own test exclude them *before* classifying (`h.scale() ≤ tol → skip`,
which is T3.H.1's identity test), so the 0-mismatch numbers are correct **as measured** — but they
are conditional on an exclusion that (T5.1e′) does not state. Two edits close it:

1. (T5.1e′) should read "`true` when `C = B = 0` **and `A ≠ 0`**; when `A = C = B = 0` the pair is
   identically zero and is removed from the candidate list by T3.H.1."
2. T5.2b.0(iii) defines `NOROOT` as "for every candidate pair `π ∈ 𝒞-list`, `h_o,π` has no root in
   `(0,ε)`". Taken literally that is never satisfiable, because every hinge point contributes an
   identically-zero pair with roots everywhere; and taken as *implemented* by the atom list it is
   weaker than the definition on exactly those pairs, which is the unsafe direction. The `𝒞-list`
   must be stated as the candidate list **with identically-zero harmonics removed** (T3.H.1).

Consequence for Proposition T5.2b′: its "O is closed" step concludes from Lemma T4.2 that some
candidate pair has `h_o,π(θ*) = 0`, contradicting `NOROOT`. If the T4.2 witness at `θ*` is a
permanent incidence — and hinge-adjacent faces always have one, at every `θ` — no contradiction
follows, so the proof needs the additional step that a permanent incidence is never the *only*
witness at an angle where the overlap status changes. `core.md` states neither the exclusion nor
this step. See R3.3.

**Verdict R3.1: AGREE** on the lemma, on the exhaustiveness of the three classes, on the atom
list for `A ≠ 0`, on the reproduced numbers, and on non-circularity. **One residue (R3.1a)**: the
`A ≠ 0` side-condition and the T3.H.1 exclusion from the `𝒞-list` are not printed.

## R3.2 — (T4.1b) split into derived `≤` and MEASURED `≥`: **AGREE**

T4.4 now prints `Θ_max ≤ min(min_e β_e, π)` tagged **[A]** with the wedge argument, and the
`≥` direction tagged **[N], MEASURED, not derived**, with the sample size (8 split-free patterns,
8.88e−16; my 4 patterns at 8.9e−16, test **T4-b**, cited), with the explicit sentences "I have
**no proof** of the `≥` direction and I do not assert one" and "That is 8 patterns, not a theorem",
and with the plausible proof route named as open. The change-log row R2.4a matches. This is
exactly what R2.4a asked for; nothing is left. **AGREE, no residue.**

## R3.3 — the certificate `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT` and the proof of T5.2b′: **AGREE**, one gap inherited from R3.1a

*Predicate.* The germ `EMB` is gone from the hypothesis and is now a **conclusion**; the second
ingredient is an exact polygon–polygon interior-disjointness test at the single explicit angle
`θ₁ = ε/2`. T5.1 and T6.4 item 2 were both updated to match, and both now say plainly that the
`θ = 0` test is withdrawn and why it is degenerate (split-edge copies coincide, adjacent faces
touch along whole edges). R2.4b is fully answered.

*Connectedness proof.* Correct as written, modulo R3.1a. `(0,ε)` is connected; `O` open is the
standard `ρ/2` continuity argument and is right; `O` closed uses compactness of the deployed
vertices (T1.A ellipses) to extract `z_k → z ∈ P_f(θ*) ∩ P_g(θ*)` with interiors disjoint at
`θ*`, which is exactly the hypothesis of **Lemma T4.2**. Lemma T4.2 itself I re-read line by line
and it is a correct piece of plane topology: the transversal case forces four local sectors, the
parallel case forces a maximal shared segment whose endpoint is a vertex of one lying on the
other's boundary. It is applied here **verbatim and legitimately** — the hypotheses it needs
(`int ∩ int = ∅`, `P ∩ Q ≠ ∅`, closed simple polygons) are the ones established at `θ*`.

*Completeness of the candidate list for this argument.* Lemma T4.2 yields "a vertex of `P` on
`∂Q` **or** a vertex of `Q` on `∂P`", and T4.1b's list is over **ordered** pairs `(w,(a,b))` with
`w ∈ face g`, `(a,b)` a consecutive pair of face `f ≠ g`, over all ordered distinct face pairs, so
both directions are covered. `NOROOT` is strictly stronger than the contact condition (it drops
the two `E2` interval tests and the overlap selection), which is conservative in the safe
direction. The single-point step is right: `θ₁ = ε/2 ∈ (0,ε)` rules out `O = (0,ε)`, and without
it the argument would only give "the status is constant", which is the failure mode T5.1's 51
counterexamples exhibit. The final `EMB` claim follows because `(0,θ₁) ⊆ (0,ε) \ O`.

*The gap.* As set out in R3.1a, the step "that vertex … has `h_o,π(θ*) = 0`, contradicting
`NOROOT`" is vacuous when the T4.2 witness is a permanent incidence, and hinge-adjacent faces
carry one at every angle. The proof needs one added sentence: identically-zero harmonics are
removed from `𝒞` by T3.H.1, and at an angle where the overlap status changes the T4.2 witness
cannot be *only* such a pair (for a hinge the transition is the `β_e` event of T4.4, which is a
genuine root). I could not find a counterexample, and the measurement below found none, but the
statement is not proved in `core.md`.

*Measured.* Test **R3-e/R3-f**: 10 all-convex corpus cases (convexity is `θ`-independent by
T1.D, so it is decided once in the flat state and separating-axis overlap is then exact), 280
positively-oriented shape-space samples at `ε = 0.006`, `NOROOT` evaluated over the complete
ordered list with the three-class atoms, `NOOVERLAP` at `θ₁ = ε/2` only, and the verdict taken
from an independent 47-point scan of `(0,ε)`:

| id | claim | samples meeting the hypothesis | violations |
|---|---|---|---|
| **R3-e** | `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT ⟹ no interior overlap anywhere on (0,ε)` | 175 / 280 | **0** |
| **R3-f** | `NOROOT ⟹ the overlap status is constant on (0,ε)` (the clopen step alone) | 176 / 280 | **0** |

`NOROOT` holds on 176 samples and the full certificate on 175, so the hypothesis has content and
`NOOVERLAP(ε/2)` is not vacuous — it excludes exactly the one sample that overlaps from the first
instant, which is the case a `θ = 0` test cannot see. **AGREE**, with the R3.1a gap carried over.

## R3.4 Tests added

Appended to `code/tests/derivation_tests.cpp` (only that file; no other file touched), built with
the header build line and run: 2 test cases, **59 982 assertions, 0 failures**.

| id | claim | samples | max error | tol |
|---|---|---|---|---|
| R3-a | `(T3.5)` `(1+τ²)h = (p+q) + 2rτ + (p−q)τ²` | 200 000 | 1.444e−13 | 1e−12 |
| R3-b1 | class 3 (`C=B=0`) is exactly `h = p(1−cos θ)` | 19 989 | 4.439e−16 | 1e−15 |
| R3-b2 | Lemma T5.1e: `sign h = sign p` on `(0,π)`, no zero, no sign change | 19 989 | 0.000e+00 | 0 |
| R3-c1 | 3-class atom list `=` amplitude/phase crossing truth, `ε = 0.2` | 1 829 226 | 0.000e+00 | 0 |
| R3-c2 | same, `ε = 0.02` | 1 829 226 | 0.000e+00 | 0 |
| R3-e | `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT ⟹ Θ_max ≥ ε` | 175 | 0.000e+00 | 0 |
| R3-f | `NOROOT ⟹ overlap status constant on (0,ε)` | 176 | 0.000e+00 | 0 |

The tests also assert the *necessity* of the class-3 clause (`mm2 > 0`, and every round-2 failure
in class 3), so they would fail if (T5.1e′) were dropped again.

## R3.5 Final verdict

R2.5, R2.4a and R2.4b are all resolved, and on all three `core.md`'s round-3 resolution is
correct: Lemma T5.1e is true and proved, the three-class list is exhaustive and exact (on two
corpora, with two independent chart-free oracles), (T4.1b) is now split and tagged honestly, and
the certificate is a pointwise test at `θ₁ = ε/2` whose connectedness proof is valid. Nothing
that was disputed in rounds 1–2 remains disputed.

One **new** residue was found in this round and it is not a labelling nit:

* **R3.1a** — (T5.1e′) is printed as "`true` when `C = B = 0`" with no `A ≠ 0` side-condition,
  while Lemma T5.1e assumes `p ≠ 0`; and T5.2b.0(iii)'s `𝒞-list` is not stated as excluding the
  identically-zero harmonics of T3.H.1. **14 928** measured candidate harmonics are identically
  zero and **all 14 928** fall in class 3 by the numeric rule, so the printed atom asserts "no
  root" for pairs whose every `θ` is a root. Both programs silently exclude them before
  classifying, so no measured number changes; what changes is that Proposition T5.2b′'s
  "`O` is closed" step needs the exclusion stated **and** needs the additional step that a
  permanent incidence is never the sole Lemma-T4.2 witness at an overlap-status change. Fix:
  two sentences in T5.2b.2 and one in the T5.2b′ proof.

**ZERO unresolved disagreements: NO** — one remains: **R3.1a** (the `A ≠ 0` side-condition on
(T5.1e′), the T3.H.1 exclusion from the `𝒞-list` of T5.2b.0(iii), and the corresponding step in
the proof of Proposition T5.2b′).

---

# Round 4 (final) — Sub-lemma T5.2b″ and the round-4 edits

Fresh context. Scope: the round-4 change log of `derivations/core.md`, the (T5.1e′) side-condition,
the `𝒞-list` definition in T5.2b.0(iii), and the new Sub-lemma T5.2b″ inside the proof of
Proposition T5.2b′. New program: the test case
`R4 Sub-lemma T5.2b'' Case A ...` appended to `code/tests/derivation_tests.cpp` (that file only;
build line unchanged, in its header). Whole file re-run: **30 test cases, 89 044 assertions,
0 failures**.

| id | claim | samples | result |
|---|---|---|---|
| **R4-a1** | (T5.2b″-1) `θ > β_e ⟹ int P_f ∩ int P_g ∩ B(p,r) ≠ ∅` | 1 731 | **0 failures** |
| **R4-a2** | (T5.2b″-1) `θ < β_e ⟹ int P_f ∩ int P_g ∩ B(p,r) = ∅` | 3 817 | **0 failures** |
| **R4-a3** | hinges excluded because a face polygon is not simple | 99 | reported, not asserted |
| **R4-b** | no face pair shares two `M′`-vertices (Case B is vacuous) | 933 | **0** |
| **R4-c1** | the `β_e` adjacent-edge collinearity pair is not identically zero | 3 819 | **0** |
| **R4-c2** | its orientation harmonic vanishes at `θ = β_e`, `\|h\|/scale` | 3 819 | 2.466e−13 (tol 1e−9) |

R4-a is deliberately **not** a sector computation: it deploys the two faces from my own `(C,S)`
basis and searches a polar grid inside `B(p,r)` (720 directions × 3 radii) for a point strictly
inside *both* face polygons, with a clearance margin, and compares that search against the predicate
`θ > β_e`. `r` is the round-4 prescription — half the distance from the hinge vertex to the face's
non-incident edges, measured in the flat state, hence `θ`-independent by T1.D. 1 306 positively
oriented null-space samples over the 16-graph corpus, 3 random hinge edges per sample, one angle
tested on each side of `β_e` with a 0.05 rad guard band.

## R4.1 — localization of Lemma T4.2 at the limit point `z`: **AGREE**

The localized form is a correct reading of the round-2 proof, not a new claim. Given
`z ∈ ∂P_f ∩ ∂P_g` with disjoint interiors, exactly three cases arise: `z` is a vertex of one polygon
(then it is a vertex on the other's boundary, since `z ∈ ∂P` of the other); the two edges through `z`
are transversal (locally each boundary is a straight line and each polygon a half-plane, so the four
sectors give an interior point — excluded); or the edges are collinear, and the maximal shared
segment `S ∋ z` stops exactly where an edge of one polygon ends, so each endpoint of `S` is a vertex
of one polygon lying on the other's boundary. The degenerate sub-case where the two collinear edges
meet only at `z` is the first case. Every witness so produced is an ordered (vertex, edge) pair of
two *different* faces, i.e. a member of the complete list of T4.1b. Correct.

## R4.2 — Case A: **AGREE** on the conclusion, with three corrections to the printed argument

* **The `θ`-independent radius is correct.** `d :=` distance from `p` to the union of the face's
  edges not incident to `p` satisfies `d ≤ ` each incident edge length (the edges at the far endpoint
  of an incident edge are non-incident and lie within that length), so `r = d/2` is smaller than both
  incident edges and `P_f(θ) ∩ B(p,r)` really is the full angular sector, reflex corners included.
  `d` is an intra-face distance, so T1.D makes it `θ`-invariant. Verified by construction in R4-a.
* **The equivalence (T5.2b″-1) is confirmed measurably**, 1 731 + 3 817 local tests, 0 failures
  (R4-a1/R4-a2). Sanity: the sector accounting gives `f` on `[θ/2, α_f + θ/2]` and `g` on
  `[−θ/2 − α_g, −θ/2]`, so the far gap is `2π − α_f − α_g − θ = β_e − θ` and the interiors meet iff
  `α_f + α_g + θ > 2π`, i.e. iff `θ > β_e`. This holds for `β_e < 0` too.
* **Correction 1 (hypothesis that is used but not printed): the face polygons must be simple.**
  Lemma T4.2, on which the sub-lemma rests, is stated for closed simple polygons, but Case A's
  sector picture uses simplicity again and independently, and T4–T6's standing hypotheses list
  positive orientation, not simplicity. This is not pedantry: with `POS` alone, my first run produced
  **5 failures of R4-a1 out of 873**, and **all 5** were hinges at which a face polygon was
  self-intersecting (`voronoi2`, `voronoi3`, reflex corners `α ≈ 5.4–6.2`, `β_e < 0`). Excluding
  non-simple faces — 99 of the tested hinges — removes every failure. `core.md` should print
  simplicity of each face polygon among the standing hypotheses of T4.
* **Correction 2 (overstatement): "`h_o,π′` vanishes at that one angle alone" is false.** With
  `a = p`, `h_o,π′(θ) = ± L_f L_g sin(β_e − θ)`, which also vanishes at `θ = β_e ± π` whenever that
  lies in `(0, π]`. The proof needs only `h_o,π′ ≢ 0` and `h_o,π′(β_e) = 0`, both of which hold, so
  the conclusion is untouched; and an extra root only makes `NOROOT` stricter, hence `R(ε)` smaller
  and still inner. The sentence should read "not identically zero", not "vanishes at that angle
  alone". Measured: 3 819 such pairs, none identically zero, `|h(β_e)|/scale ≤ 2.5e−13`.
* **Correction 3 (a step that is used silently): `f` and `g` sharing a hinge point are adjacent
  across exactly one hinge edge `e`, so "the" `β_e` exists.** The proof writes `β_e` for the pair
  `(f, g)` without establishing this. It is true, and the argument is short: by 0.2 a hinge edge
  keeps `v` only at `src(e)`, and its direction is the `σ`-induced half-edge direction; in the
  boundary traversal of one face the two hinge edges at `v` are consecutive, one entering `v` and
  one leaving it, so a face keeps `v` for **at most one** of its incident hinge edges. The
  identification of copies of `v` is therefore a *matching* on the faces around `v`, never a fan of
  three or more, so two faces sharing a hinge point are hinge-adjacent. Without this, a fan
  `f₁—f₂—f₃` all keeping `v` would make Case A's single-`β_e` accounting wrong. Measured
  consequence: R4-b, 0 of 933 face pairs share two or more `M′`-vertices.

## R4.3 — Case B: **AGREE**, and the case is in fact vacuous

The kinematics check the round-4 log asks for: `γ_f` is `R(−σ_f θ/2)` plus a translation (0.7′/T1.5),
so `R(θ) = γ_g^{-1} ∘ γ_f` is a rotation by `(σ_g − σ_f)θ/2`. Two faces sharing a hinge point are
hinge-adjacent (R4.2, correction 3), hence `σ_g = −σ_f`, so `R(θ)` is a rotation by **`±θ`** — this
is exactly T1.C, the statement that every hinge opens by `θ`. A rotation by `±θ` is not constant on
any interval of `θ`, so the printed conclusion "`R(θ)` is constant, hence the overlap predicate is
`θ`-independent, contradiction" is right, and it can be sharpened: **two faces can never share two
distinct material points at all**, so Case B is vacuous rather than merely contradictory. That also
answers the structural question — sharing two hinge points would need two hinge edges between the
same face pair, each keeping a different endpoint, and the rotation argument forbids it outright.
Measured: **0 / 933** face pairs share two `M′`-vertices (R4-b). The proof is correct as printed;
the vacuity is a strictly stronger statement and worth printing in its place.

## R4.4 — the stated uncovered hypothesis and the two interval atoms `s₁, s₂`: **DISAGREE**

`core.md` states openly that permanently collinear but non-coincident triples are not covered, and
proposes `s₁ = ⟨w−a, b−a⟩`, `s₂ = ⟨w−b, a−b⟩` as "the additional candidate atoms needed in that
case". **They do not close the case.** Split the transition at such a witness `z = w(θ*)`:

* **`w` at an endpoint of `[a,b]` at `θ*`.** Then `s₁` or `s₂` has a root at `θ*` and the proposed
  atoms do fire. This sub-case is closed by `s₁, s₂`.
* **`w` in the relative interior of `[a,b]` on a whole neighbourhood of `θ*`.** Then `s₁, s₂ > 0`
  there and neither has a root: **the proposed atoms are silent at exactly the transition they are
  supposed to catch.** The overlap still begins, because near `w` face `f` is a half-plane bounded by
  the line `ab` and face `g` is an angular sector at `w`; the interiors start to meet when a bounding
  ray of `g`'s sector crosses the line `ab`, i.e. when the *neighbour* vertex `w′` of `w` in `g`
  becomes collinear with `a, b`. The atom that fires is the ordinary orientation harmonic of the pair
  `(w′, (a,b))`, already in the complete list, and it is not identically zero whenever the relative
  rotation is non-zero, i.e. whenever `σ_f = −σ_g`; when `σ_f = σ_g` the relative placement is
  constant in `θ` (same computation as R4.3) and no transition can occur at all. This is Case A's
  own mechanism, transplanted, not an interval atom.

So the correct statement is: the case is closed by `s₁, s₂` **together with** the neighbour-vertex
collinearity pairs — and the latter cost nothing, being already in the complete candidate list with
`h ≢ 0`. As printed, the hypothesis is not closed by what `core.md` names. Separately, the
hypothesis is unnecessary for the pairs where the identically-zero harmonics actually arise: for
hinge-adjacent `f, g` the relative motion is a pure rotation by `±θ` about the shared hinge point
`p`, so in `f`'s body frame a material point `w` of `g` traces a **circle** of radius `|w − p|`
centred at `p`; a circle lies on a fixed line only if its radius is zero, i.e. `w = p`, a
coincidence. Permanent collinearity without coincidence is therefore impossible for hinge-adjacent
pairs — which is precisely where all 14 928 measured identically-zero harmonics live.

This is the one item of round 4 on which I disagree with `core.md` as written. It is a gap in a
*stated* hypothesis, not in a claimed theorem, and the fix is one sentence plus one named pair.

## R4.5 — (T5.1e′) side-condition, T5.2b.0(iii) `𝒞-list`, degree ≤ 4: **AGREE**

* **(T5.1e′) with `A ≠ 0`** is exactly the fix R3.1a(a) asked for, and it matches Lemma T5.1e's own
  hypothesis `p ≠ 0`. Confirmed by R3-b1/R3-b2, re-run unchanged.
* **T5.2b.0(iii)'s `𝒞-list = {pairs with f ≠ g} \ {π : h_o,π ≡ 0}`, removed by the identity test of
  T3.H.1 before the three-class split**, is the definition both programs implement, and it makes
  (iii) satisfiable. Correct, and it now says so.
* **Degree ≤ 4 survives.** `s₁, s₂` are the dot-product harmonic of (T3.3), i.e.
  `p′ + q′ cos θ + r′ sin θ` with `p′, q′, r′` bilinear in `(C, S)` and therefore quadratic in `t`
  (T3.2), so the same `g(τ) = C + Bτ + Aτ²` three-class list decides them and the atoms
  `g(0)g(T)`, `disc`, `A·g(0)`, `A·g(T)`, `A·B`, `A·B + 2A²T` stay at degree ≤ 4. The pairs I add in
  R4.4 are ordinary orientation harmonics already in the list and add no new atom shape. The
  standing caveat is unchanged and correctly printed in `core.md`: the degree bound is a statement
  about `POS ∧ NOROOT` only, never about `NOOVERLAP(ε/2)`.

## R4.6 Final verdict

Per item: R4.1 **AGREE** · R4.2 **AGREE** (three corrections to the argument, none to the
conclusion) · R4.3 **AGREE** (and vacuous) · R4.4 **DISAGREE** · R4.5 **AGREE**.

Residues for a round 5, in priority order:

1. **R4.4** — `s₁, s₂` do not close the permanently-collinear non-coincident case; the
   neighbour-vertex collinearity pair must be named alongside them. Substantive.
2. **R4.2 correction 1** — simplicity of each face polygon is used by T4/T5 and is not among the
   standing hypotheses; 5 of my first 873 local tests failed on exactly the self-intersecting faces.
3. **R4.2 correction 2** — "vanishes at that one angle alone" should read "is not identically zero".
4. **R4.2 correction 3 / R4.3** — the matching argument (a face keeps a vertex for at most one of its
   hinge edges) should be printed: it is what makes `β_e` well defined for the pair in Case A, and it
   makes Case B vacuous.

**ZERO unresolved disagreements: NO** — one remains: **R4.4** (the two interval atoms `s₁, s₂` do
not by themselves close the hypothesis Sub-lemma T5.2b″ states openly).

---

# Round 5 (final) — the R4.4 replacement proof, T5.2b‴, and the Case-A corrections

Fresh context. Scope: `derivations/core.md` "Change log — round 5", Sub-lemma T5.2b″ (Case A with
(T5.2b″-1b)/(T5.2b″-1c) and the endpoint/interior sub-cases of [D2]), Sub-lemma T5.2b‴ [A0], the
simplicity hypothesis in the standing hypotheses of T4–T6, and Remark [C] (Case B vacuous). New
program: the test case `R5 same-sigma pairs: pure relative translation, constant local predicate;
and the zero set of h_o,pi' in (0,pi)` appended to `code/tests/derivation_tests.cpp` (that file only;
build line unchanged, in its header). Whole file re-run: **31 test cases, 89 055 assertions,
0 failures.**

| id | claim | samples | result |
|---|---|---|---|
| **R5-a1** | `σ_f = σ_g ⟹` every relative direction is `θ`-independent | 14 496 | 1.954e−13 (tol 1e−9) |
| **R5-a2** | `σ_f = σ_g ⟹` the local predicate (T5.2b″-1c) never flips | 14 496 | **0 flips** |
| **R5-a3** | power check: `σ_f = −σ_g` drifts by exactly `±θ`, and the predicate does flip | 17 208 | 1.901e−13; **5 303 flips** |
| **R5-c1** | `h_o,π′` has **at most one** root in `(0, π)` | 20 000 | **0** with two or more |
| **R5-c2** | `β_e ∈ (0, π) ⟹` that root **is** `β_e` | 5 038 | max err 4.44e−16 |
| **R5-c3** | `θ > β_e` does not change value at `β_e ± π` | 20 000 | **0** changes |
| **R5-d** | (T5.2b″-1b) `h = ± L L′ sin(β_e − θ)` on the corpus | 894 | 2.711e−13 (tol 1e−9) |

R5-a is not a re-run of a `core.md` formula: it deploys both faces from my own `(C, S)` basis, reads
the *directions* `dir(b−a)` of `f` and `dir(w→w_next)`, `dir(w→w_prev)` of `g` off the deployed
polygons on a 13-point `θ`-grid in `(0, π]`, forms the arc `A` and the open half-circle `J` of
(T5.2b″-1c) from them, and evaluates "`A ∩ J ≠ ∅`" as a circle-arc predicate with a `1e−7` guard
band. R5-c is synthetic: 20 000 random `(β_e, L, L′)` with `β_e ∈ (−2π, 2π)`, roots of
`sin(β_e − θ)` found by a 4 000-interval scan plus 200 bisections.

## R5.1 — (T5.2b″-1c), the compactness/continuity argument, and the neighbour pair: **AGREE**

The chain is correct. `A(θ*) ⊆ K(θ*)`; `± u(θ*) ∉ A(θ*)` once neither is an endpoint of `A(θ*)`
(an arc inside a closed half-circle cannot contain that half-circle's endpoint in its relative
interior without leaving it); `A(θ*)` compact inside the open arc `K(θ*) \ {± u(θ*)}` therefore
persists on a neighbourhood, both arcs having continuously moving endpoints; and `z_k ∈ B(w(θ_k), r)`
for large `k` gives the contradiction. Hence a bounding ray is parallel to `ab`, and since
`w ∈ ab` it is collinear with `a, b`.

**Yes to both parts of the question.** The neighbour pair `(w′, (a,b))` vanishes at the transition
by construction (that is what "endpoint of `A(θ*)` equals `± u(θ*)`" says), and it necessarily has
`h ≢ 0`: `h ≡ 0` forces `w′ − w ∥ b − a` for all `θ`, hence `σ_f = σ_g` by (0.7′), hence a pure
relative translation, hence a `θ`-independent local predicate — R5-a1/R5-a2 confirm exactly that
implication measurably, and R5-a3 shows the predicate genuinely varies otherwise (5 303 real flips),
so the contradiction is not vacuous.

Two wording corrections, neither touching the conclusion:

* **(i) `A` should be the open arc, or the equivalence should be argued.** The direction set of
  `int P_g` at `w` is the *open* sector, not the closed arc `A`. The printed (T5.2b″-1c) is
  nonetheless correct, because `J` is open and `α_g > 0`: if an endpoint of `A` lies in `J`, so do
  nearby interior directions. One clause ("`A ∩ J ≠ ∅ ⟺ int A ∩ J ≠ ∅` since `J` is open") repairs it.
* **(ii) the localization condition is under-stated.** "`B(w(θ), r)` meets no vertex of `f` other
  than through the edge `(a,b)`" does not exclude a *non-incident edge* of `f` crossing the ball. The
  condition actually used is `B(w(θ), r) ∩ ∂P_f ⊆ relint(a,b)`, which the same compactness argument
  supplies (`f` simple, `w(θ*)` in `relint(a,b)`, `I` compact).

## R5.2 — `σ_f = σ_g ⟹` pure relative translation ⟹ locally constant predicate: **AGREE**

Derivation: `γ_f = R(−σ_f θ/2) + t_f` (0.7′), so `R(θ) = γ_g^{-1} ∘ γ_f` has rotation part
`R((σ_g − σ_f)θ/2) = I` when `σ_f = σ_g`; a translation preserves directions, so `u`, `A`, `J` and
the predicate `A ∩ J ≠ ∅` are all `θ`-free while only the base point `w` slides along `ab`.
Measured on the corpus: relative directions constant to **1.954e−13** over 14 496 comparisons, and
**0 / 14 496** predicate flips, against **5 303 / 17 208** flips on the opposite-`σ` control.

## R5.3 — Sub-lemma T5.2b‴ (hinge adjacency via matching): **AGREE**

By 0.2 the direction `src(e) → dst(e)` is the common `σ`-induced half-edge direction, so in the
`σ_f`-traversal of `∂f` a hinge edge `e` at `v` keeps `v` for `f` exactly when `e` is the unique edge
*leaving* `v`. Build the graph on faces around `v` whose edges are the hinge edges at `v` with
`src = v`: an edge `e` there is *simultaneously* the leaving edge of both its faces, so every face
has degree ≤ 1. That is a matching, never a fan, and two faces sharing a hinge point are joined by
exactly one hinge edge, so "the" `β_e` is well defined. Confirmed by R4-b (0 / 933). Note the step
"each face visits `v` at most once" is simplicity again — now a standing hypothesis, so this is
covered, but it is worth a half-sentence in [A0]. Remark [C] follows: hinge adjacency gives
`σ_g = −σ_f`, `R(θ)` is a rotation by `±θ`, and for `θ ∈ (0, π]` that fixes exactly one point.

## R5.4 — the zero set of `h_o,π′` and the sector transition: **AGREE on the conclusion, one printed sentence is false**

* **What the proof needs is right.** `θ* = β_e ∈ (0, ε) ⊂ (0, π)` is a root of `h_o,π′` and
  `h_o,π′ ≢ 0`. The closed form (T5.2b″-1b) is confirmed on the corpus to **2.7e−13** (R5-d), and
  `h_o,π′` has **at most one** root in `(0, π)` in all 20 000 synthetic cases (R5-c1), which is that
  root whenever `β_e ∈ (0, π)` (R5-c2, 5 038 / 5 038).
* **Correction: "exactly one of `{β_e, β_e + π, β_e − π}` lies in `(0, π)`" is false in general.**
  It fails whenever `β_e ≤ −π` — measured **5 008 / 20 000** synthetic `β_e ∈ (−2π, 2π)`, where none
  of the three lies in `(0, π)` and the root in range is `β_e + 2π`. The correct general statement is
  "`sin(β_e − θ)` has **at most one** root in the open interval `(0, π)` of length `π`, and exactly
  one unless `β_e ≡ 0 (mod π)`". `β_e ≤ −π` needs `α_f + α_g ≥ 3π`, i.e. both faces strongly reflex
  at the shared hinge, and is outside the case at hand (`β_e = θ* ∈ (0, ε)`), so **no theorem is
  affected** — but the enumerating clause and its three-way case split should be deleted or replaced.
* **"Only `β_e` is a sector-overlap transition" is right**, and trivially so: the predicate is
  `θ > β_e`, whose only transition is `θ = β_e`. Checked at `β_e ± π`: **0 / 20 000** changes (R5-c3).

## R5.5 — nothing else in T5.2b′'s closedness argument was weakened: **AGREE**

Openness, the compactness/subsequence extraction, the appeal to Lemma T4.2 and its localization at
`z`, and the "some witness has `h ≢ 0`" branch are unchanged from round 4, where I checked them
(R4.1). Round 5 only adds hypotheses (simplicity) and replaces the round-4 `s₁, s₂` claim by a proof;
it removes nothing. `s₁, s₂` are retained, correctly, for the endpoint sub-case only.

One presentational gap, not a truth-value gap: [A] assumes "**every** witness at `z` is a permanent
coincidence" and [D2] assumes "**every** witness at `z` is permanently collinear non-coincident", so
the **mixed** case (both kinds present at one `z`) is not addressed by either heading as printed. It
composes without any new argument: Case A needs only `z ∈ H`, and [D2]'s proof uses only the *one*
witness it is applied to, so in the mixed case with `z ∉ H` pick any permanently-collinear witness and
run [D2]. Recommended fix: state [D2] for a single witness rather than for all of them.

## R5.6 Final verdict

Per item: (1) (T5.2b″-1c), compactness, and the neighbour pair — **AGREE** (two wording corrections)
· (2) `σ_f = σ_g ⟹` pure translation ⟹ locally constant predicate — **AGREE**, measured
· (3) Sub-lemma T5.2b‴ — **AGREE** · (4) the zero set and the sector transition — **AGREE on the
conclusion**, with one printed sentence false and to be replaced · (5) nothing else weakened —
**AGREE**.

Non-blocking corrections to print, in priority order:

1. **R5.4** — delete "exactly one of `{β_e, β_e ± π}` lies in `(0, π)`" and its three-way split;
   write "at most one root in `(0, π)`, and it is `β_e` here since `β_e = θ* ∈ (0, ε)`".
2. **R5.1(i)** — `A` in (T5.2b″-1c) is the closed arc; add the one clause that makes
   `A ∩ J ≠ ∅ ⟺ int A ∩ J ≠ ∅` (`J` open, `α_g > 0`).
3. **R5.1(ii)** — state the localization as `B(w(θ), r) ∩ ∂P_f ⊆ relint(a,b)`, not "meets no vertex".
4. **R5.5** — state [D2] for a *single* permanently-collinear witness, closing the mixed case.
5. **R5.3** — note in [A0] that "a face visits `v` once" is the simplicity hypothesis.

**ZERO unresolved disagreements: YES.** The five items above are corrections to printed wording, each
localized to one sentence; none changes a hypothesis, a proof strategy, or the truth value of any
theorem, and the substantive round-4 residue (R4.4) is fully conceded and correctly replaced by the
[D1]/[D2] proof, which I have now verified both logically and numerically.

---

# Round 6 — Checker pass on the F32 `NOROOT` amendment

Fresh context. Scope: the amendment proposed in `results/kill/jitter/cert_diagnosis.md` §4 —
`T5.2b.0(iii)` becomes "no **admissible** root in `(0, ε)`", i.e. `𝒞(X) ∩ (0, ε) = ∅`, with the two
interval tests `E2` of T4.1b included — read against `core.md` T3, T4.1b, T4.2′/T4.2″, T5 (T5.1,
T5.2a, T5.2b, T5.2b′, T5.2b″ with Cases A/B and [D1]/[D2]/[D3], T5.3, T3.H.1) and the fixed code
`code/src/method/contact.cpp` (`validity_certificate`, `contact_angles`) with the F32 doctest in
`code/tests/test_method.cpp`. Two new programs: `derivations/scratch/check_t5_adm.cpp` (standalone,
build line in its header) and the `R6` test case appended to `code/tests/derivation_tests.cpp`
(that file only), plus `derivations/scratch/check_b4_93.cpp` and the `R6-c` case for the B4 add-on.
Whole file re-run: **33 test cases, 89 074 assertions, 0 failures.**

| id | claim | samples | result |
|---|---|---|---|
| **R6-a1** | Case A's substitute root at `β_e` is ADMISSIBLE (`E2` holds) | 9 642 | **0 failures** |
| **R6-a2** | its projection ratio `⟨w−a,b−a⟩/‖b−a‖²` equals `L′/L` | 9 642 | 1.771e−14 (tol 1e−9) |
| **R6-a3** | `h_o,π′(β_e) = 0` (round-5 claim, re-run under jitter) | 9 642 | 3.531e−14 (tol 1e−9) |
| **R6-b** | COMPLETENESS is false: a graze with `Θ_max ≥ ε` is not certified | 1 | **counterexample found** |
| **R6-c1** | `voronoi_93` faces (94, 184): the interiors are DISJOINT (dense probe) | 1 442 401 | **0 common interior points** |
| **R6-c2** | `polygons_overlap` misfires there (pinned referee defect) | 6 | **6 / 6 shrink values wrong** |

R6-a splits the ratio into **9 062** pairs strictly inside the segment and **580** exactly at the far
endpoint (`L′ = L`, the vertex-vertex case), over the whole corpus with null-space jitter `0.12`,
simple faces only. Range of the ratio: `[0.007307, 1.000000]`. The standalone program repeats A on
the 7 tilings only (1 363 pairs, ratio `[0.540396, 1.000000]`, `|s/‖e‖² − L′/L| ≤ 1.554e−15`).

## R6.1 — soundness of `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT_adm`: **VERIFIED**

**Verdict: T5.2b′ survives.** The amendment weakens the hypothesis (admissible roots are a subset of
roots, so `NOROOT_old ⟹ NOROOT_adm`), so soundness has to be re-proved, and it goes through. The
proof needs an admissible root in each of the branches that contradicts `NOROOT`, and there are
three:

1. **The main branch is admissible for free.** A witness at `z` is a vertex of one polygon *on the
   boundary of the other* — Lemma T4.2's localization delivers exactly that, in both the
   vertex-on-edge case and the shared-segment case (each endpoint of the maximal shared segment is a
   vertex of one polygon **on** an edge of the other). "On the boundary" is `E2` at `θ*` with
   non-strict inequalities. The diagnosis's own soundness argument ("the connectedness argument needs
   precisely `𝒞(X) ∩ (0,ε) = ∅`") is correct, and this is why. **VERIFIED.**

2. **Case A: the diagnosis's one flagged gap closes, cleanly.** The substitute pair
   `π′ = (w,(a,b))`, `a = p` the shared hinge point, `w` the far endpoint of the **shorter** far-side
   edge. At `θ = β_e` the far-side gap `β_e − θ` vanishes, so the two far-side edges are collinear
   **and co-directed** — the anti-parallel configuration is `β_e ± π`, which round 5 (R5.4) already
   separated out as a zero of `h_o,π′` that is not a sector transition. Hence with `L = ‖b − a‖`,
   `L′ = ‖w − a‖`, both `θ`-invariant by T1.D,
   ```
        ⟨ w − a , b − a ⟩ = L L′ ,     ‖ b − a ‖² = L² ,     0 < L′ ≤ L ,
   ```
   so `E2` holds, strictly on the lower side and with **equality on the upper side iff `L′ = L`**.
   The segment of T4.1b is **closed**, so equality is admissible; `contact.cpp`'s inclusive
   `-tol ≤ s ≤ l2 + tol` implements the right convention and is load-bearing, not cosmetic —
   **580 of 9 642** measured Case-A pairs sit exactly at that endpoint, and every regular tiling in
   the corpus produces only that case at `X₀` (all 34 unjittered pairs have ratio exactly `1`).
   **VERIFIED**, with the caveat that a *strict* interval test would break Case A on precisely the
   patterns both papers ship.

3. **[D2] interior sub-case: a real gap the diagnosis did not name, and it is repairable.** The
   round-5 neighbour pair `π″ = (w′,(a,b))` has `h_o,π″(θ*) = 0` and `h_o,π″ ≢ 0`, but `w′(θ*)` is
   only known to be on the **line** `ab`; its distance from `w` is the edge length of `g`, which
   nothing bounds. With `a = 0`, `b = L`, `w = c ∈ (0,L)`, `w′ = c + d` on the line at `θ*`:
   `c + d ∈ [0, L]` makes `π″` itself admissible; `c + d > L` puts `b(θ*)` in the relative interior
   of `[w(θ*), w′(θ*)]`, so the **role-swapped** pair `(b, (w,w′))` is admissible and is in the
   complete ordered list; `c + d < 0` is symmetric with `a`. The swapped pair is not identically
   zero: `h_o,(b,(w,w′)) ≡ 0` plus `w` permanently on `ab` forces `line(ab) = line(w w′)` on a
   punctured neighbourhood of `θ*` (`‖y_w − y_b‖²` is a non-zero harmonic, so `w ≠ b` off at most two
   angles), hence `h_o,π″ ≡ 0`, which round 5 refuted. **VERIFIED** after the insertion; without it
   the amended proof has a hole here, not only in Case A.

**Assessment of the diagnosis's own caveat.** `cert_diagnosis.md` §4 flagged exactly one gap (Case A
admissibility) and called its geometric argument "stated here as the one gap the amendment opens, not
as a proved step". The flag was right and the geometric sketch was right; the sentence "the interval
predicate holds with equality at worst — which the inclusive tolerance admits" is exactly the
`L′ ≤ L` fact, now proved and measured. But the diagnosis **missed** item 3 above: the same
admissibility question applies to [D2]'s neighbour pair and is *not* answered by the same argument
there. Neither gap is fatal; both are now closed in `core.md` round 6 (R6.2).

## R6.2 — `[D3]` and the interval atoms `s₁, s₂`: **DISPUTED, and withdrawn in the Deriver's favour**

Round 4/5 left [D2]'s **endpoint** sub-case (`w` at an endpoint of `[a,b]` at `θ*`, with
`h_o,(w,(a,b)) ≡ 0`) to the two extra atoms (T5.2b″-2), and [D3] instructed that they "must be added
to `NOROOT`". `contact.cpp` never implemented them. **The code was right and [D3] was unnecessary:
the endpoint sub-case cannot occur.**

*Proof.* Say `w(θ*) = a(θ*) =: z`. Then the vertex `a` of `f` coincides with the vertex `w` of `g`,
so `a` lies on `∂P_g(θ*)`, on **both** edges of `g` at `w`. Hence `(a,(w,w′))` and `(a,(w,w″))` are
both witnesses at `z`, and [D2]'s hypothesis makes both identically zero: `a` lies permanently on the
line `w w′` and on the line `w w″`. Those two lines are distinct and meet only at `w` provided the
interior angle `α_g` of `g` at `w` is not `0`, `π` or `2π`; `α_g` is `θ`-invariant (T1.D). So
`a ≡ w`, a permanent **coincidence** — contradicting [D2]'s non-coincidence hypothesis, and landing
in T3.H.2 / Case A instead. ∎

This needs one hypothesis that `core.md` never printed: **no straight vertices** (`α ∉ {0, π, 2π}` at
every vertex of every face). It is implied by any reasonable reading of "simple polygon" and is
removable by deleting redundant vertices, but it is now *used*, so it must be stated — the Deriver
has added it as (H6). With it, `E2` appears in the certificate **only** as the admissibility filter
on orientation-harmonic roots, never as atoms of its own, and derivation and code agree for the first
time since round 4.

## R6.3 — completeness (`Θ_max ≥ ε ⟹ certificate`): **DISPUTED — it is FALSE**

The jitter re-measurement (0/5064 exceptions in either direction at `ε = 0.006`) does **not**
generalize, and `core.md`'s own T4.2 table contains the counterexample.

By R6.1 the amended `NOROOT_adm(ε)` *is* `𝒞(X) ∩ (0, ε) = ∅`, i.e.
```
        NOROOT_adm(ε)  ⟺  ε ≤ θ₁(X) = min 𝒞(X) ,          while       U(ε)  asks  ε ≤ Θ_max(X) ,
```
and by T4.2″ `Θ_max = θ_{i*}` with `i* ≥ 1` whenever the range is positive, so `θ₁ ≤ Θ_max` with
**equality iff the first contact is already an overlap transition** — no graze. That is exactly
Proposition T4.3's genericity hypothesis, which T4.3 itself records as failing identically on the
symmetric patterns of both papers.

**Counterexample, measured** (`check_t5_adm.cpp` §B, `derivation_tests` R6-b): `hexagons` at `X₀` has
`θ₁ = 1.047198 = π/3` (the split-duplicate end-to-end contact, a graze) and `Θ_max = 2.094395 = 2π/3`.
At `ε = 1.570796` the exact range clears `ε` while the certificate returns
`(POS, NOOVERLAP, NOROOT) = (1, 1, 0)` with first admissible root `1.047198`. So
`X₀ ∈ U(ε) \ R(ε)`: **still a strict inner approximation.**

**What the jitter number actually says.** `ε = 0.006` sits below `θ₁` on every certified row of that
corpus; on the regular tilings `θ₁ ≥ 0.83`, three orders of magnitude above `ε`. Exactness there is a
property of *(corpus, ε)*, not of the certificate. `cert_diagnosis.md` §6's sentence "the certificate
is now exact on this corpus" is correctly scoped; §7's "on K2a's corpus the certificate is now
*exactly* the valid set, 173 = 173" is also correctly scoped; but the standing claim
"`a*(certified)` now coincides with `a*(Θ_max > 0)`" must **not** be promoted to "the certificate is
complete". Any downstream text that does so is wrong at `ε` above the graze scale.

**Tag: CONJECTURE (T5.2b-complete)** — that `U(ε)` admits a quantifier-free description with an atom
count bounded independently of `t`. A complete certificate needs T4.2″'s overlap *selection*, one
probe per interval of `(0, ε) \ 𝒞(X)`; that is exact and finite but its atom structure depends on
`t`, so semialgebraicity (T5.2a) survives while the bounded description does not.

## R6.5 — B4 `voronoi_93`, scan `0.2484` vs bisection `0`: **REFEREE ARTEFACT**, the scan is right

Replayed through `kill_b4.cpp`'s own pipeline (same graph, same `sigma_def`, same `free` system, same
seeds and weights) by `derivations/scratch/check_b4_93.cpp`. The row reproduces exactly:
`Theta_max` scan `= 0.248400`, referee bisection `= 0`, and `|𝒞(X)| = 3234` with
`min 𝒞(X) = 0.248400`. Pruning is not the cause: the scan gives `0.248400` on the **pruned** 3 625
pairs and on **all** 201 295 pairs, with the same `𝒞`.

**Where the referee trips.** `has_collision(..., 1e-12)` fires from `θ = 1e-7` up to `θ = 0.036255`
on face pair **(94, 184)**, and `referee_theta` returns `0` on its `col(1e-7)` probe before the
bisection even starts. The exact scan reports **no** contact angle below `0.2484` for that pair — its
first is `1.809342`.

**The scan is right.** The two faces share exactly one `M′`-vertex, a hinge point, with interior
angles `α_f = 0.368615` and `α_g = 4.105228`, so `β_e = 2π − α_f − α_g = 1.809342` — *exactly* the
first contact angle the scan reports for the pair, in agreement with (T5.2b″-1). Both faces are
**simple** and positively oriented. A `1 200 × 1 200` grid over the union bounding box at `θ = 0.01`
finds **73 735** points inside face 94, **342 350** inside face 184 and **0 inside both**; 17 280
radial probes around the shared hinge point at six radii find **0** common interior points at every
`θ ∈ {1e−4, 0.01, 0.03, 0.036}`. There is no interior overlap to detect, hence nothing missing from
`𝒞(X)`, and Corollary T4.2′ is not violated.

**Where the referee is wrong.** `polygons_overlap`'s `shrink_poly` displaces the **shared hinge
vertex** by `−s (p − centroid)`, a *different* vector in each of the two copies, so the two edges
meeting exactly at that vertex acquire orientation determinants of size `1e−11 … 1e−17` whose signs
then decide the strict `seg_intersect` crossing test. Every spurious crossing found in the replay is
hinge-incident — `A edge(545,546) × B edge(546,484)` and `A edge(546,543) × B edge(546,484)`, with
determinants `2.1e−11`, `3.7e−12`, `1.5e−17`. The answer flips with `shrink` and with `θ` for no
geometric reason (`shrink = 0` gives `0, 1, 0, 1, 0` across `θ = 1e−4, 0.01, 0.03, 0.036, 0.04`).
This is precisely the degeneracy `shrink` is documented to prevent ("faces that merely touch at a
hinge vertex … do not register as overlapping"); shrinking each polygon toward **its own** centroid
does not achieve that when one of the two sectors at the shared vertex is **reflex** (`α_g > π`
here), because the displacement can carry the reflex apex across the neighbour's thin sector.

**Not an admissibility-filter false negative.** The new `NOROOT` plays no part: `cert_noroot = 0` on
this row anyway, `eps_max = 0.2484` comes from `first_root`, and the pair's genuine on-segment
contact at `β_e = 1.809342` is *kept* by the filter, not removed. The filter's job here is exactly
what it does — it declines to call the permanent hinge incidence a contact, which is right.

**`voronoi_96` agrees** (`0.241884` both ways) because no such pair exists there: its first
`has_collision` is at `θ = 0.2484`, the genuine event.

**Consequence for the bundle.** `theta_bisect` and any statistic built on it (`worst_referee_gap`,
`referee_agree`, and the `col(1e-7)` zero-range shortcut shared by K2a/K2b/K5/K6) can report a
spurious `Θ_max = 0` whenever a design has a hinge whose two sectors are thin-against-reflex. It
biases the referee **downward** only, so it cannot manufacture a PASS; the K5/K6 FAIL verdicts, which
rest on `zero_plus_q ≤ 0` and on the exact scan, are unaffected. But **`theta_bisect` must not be
quoted as ground truth against the exact scan**, and B4's referee-agreement line needs the caveat.
Fixing `polygons_overlap` (offset each polygon along its own inward edge normals rather than toward
its centroid, or exclude edge pairs that share an endpoint) is a change to `code/src/core/collision.cpp`
that would move numbers across the whole kill bundle, so it is reported here, not applied.

**Regression case added.** `derivation_tests` `R6-c` pins the reproduction with the two polygons as
literal coordinates — no cache, no repair, no generator — asserting the ground truth (`0` common
interior points of `1 442 401` probes, `α_g > π`, `β_e = 1.809342`) and pinning the misfire
(`6 / 6` shrink values wrong). It is a pinned reproduction, to be inverted when the predicate is
fixed.

## R6.4 Final verdict

Per item: (1) **soundness of the amended certificate — VERIFIED**, with two insertions the diagnosis
did not supply (Case A admissibility via `L′ ≤ L` and the closed interval; the role-swapped pair in
[D2]'s interior sub-case) · (2) **[D3] withdrawn — VERIFIED** that the endpoint sub-case is vacuous
under the new hypothesis (H6), so `s₁, s₂` are unnecessary and the code was already right ·
(3) **completeness — DISPUTED, refuted** by the hexagon graze, with the sharp statement
`NOROOT_adm(ε) ⟺ ε ≤ min 𝒞(X)` · (4) **B4 `voronoi_93` — REFEREE ARTEFACT**: the exact scan and the
repaired certificate are both correct at `0.2484`; the bisection referee's `0` is a
`polygons_overlap` misfire at a shared hinge vertex with a reflex sector, unrelated to the
admissibility filter.

**Blocking corrections: none.** Required edits, all now present in `core.md` round 6: amend
T5.2b.0(iii) as in `cert_diagnosis.md` §4 with the word *closed* on the interval test; insert the
Case A and [D2] admissibility paragraphs into T5.2b″; print (H6); withdraw [D3]; add R6.6 so that no
sentence claims completeness. **ZERO unresolved disagreements: YES.**
