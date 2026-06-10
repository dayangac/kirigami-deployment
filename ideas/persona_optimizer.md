# Ideas — persona: OPTIMIZER

Numerical optimization / geometry processing. What I want from an idea: a crisp objective,
an exact structure that removes a solver, and a baseline it beats by an explained margin.

## Inputs actually used

- `STATE.md` F1–F14 (F11 partition, F13 body-and-hinge count, F14 `L = R·D`).
- `notes/paper_2026.md` — full (Secs. 3, 6, 7, 8, 10, 12, 13, 15 read line by line).
- `notes/paper_2025.md` — Secs. 3 (Eqs. 1–5), 13 (limitations), 14 (missing supplement).
- `notes/field_kirigami.md` (25 rows, seeds 1–7), `notes/field_rigidity.md` (26 rows, Sec. 0),
  `notes/field_tutte.md` (20 rows, Secs. 0 and 2). All three exist and were read.
- `papers/related/isogami.txt` Secs. 3.1–3.2, 4.1–4.2, 6.
- `code/src/core/*.hpp` — to state kill experiments against routines that exist.

## Two facts I was handed, and my own re-derivation of them

**(i)** *(claimed verified by ideator-geometer / ideator-rigidity; I re-derived it, see below)*
Uniform deployment positions are `Y_θ = cos(θ/2)·C(X) + sin(θ/2)·S(X)` with `C, S` **linear in X**.

Re-derivation from the convention in `code/src/core/kinematics.hpp`: face `f` maps
`y = R(−σ_f θ/2) x + t_f`. Across a hinge edge `(f,g)` at hinge vertex `u`, `σ_g = −σ_f`, so
`t_g − t_f = (R(−σ_f θ/2) − R(σ_f θ/2)) x_u = −2σ_f sin(θ/2) J x_u` with `J` = rotation by 90°.
Seeding `t_root = 0` and propagating by BFS gives `t_f = sin(θ/2) b_f` with `b_f` linear in `X`.
Hence `y_i = cos(θ/2) x_i + sin(θ/2)(−σ_f J x_i + b_f)`, i.e. exactly (i). **I checked this
algebraically by hand; I did not run it.**

**(ii)** *(claimed verified by ideator-rigidity; cited, not re-derived)* Infinitesimal mobility of
the cut structure is `dim ker A − c(Γ)` for a `2·b₁(Γ) × |F|` matrix `A` with one scalar per face,
and `σ ∈ ker A` is Eq. (2).

**Everything I build on top of (i) and (ii) below is my own and is flagged as unverified
numerically.** Nothing here has been run.

## The one observation that generates half my list

Because of (i), write `c = cos(θ/2)`, `s = sin(θ/2)`. *Every* quadratic function of the deployed
positions — orientation determinants, polygon areas, cross products, periodicity determinants — is
a quadratic form in `(c,s)`:

```
   Q(θ) = A c² + B cs + D s²
        = (A+D)/2  +  (A−D)/2 · cos θ  +  (B/2) · sin θ
        = p + q cos θ + r sin θ                                            (★)
```

with `p, q, r` **quadratic in X**. A single harmonic. Its zeros in `(0, π]` are closed form:
`cos(θ − atan2(r,q)) = −p/√(q²+r²)`. So the whole deployment — collisions, closure, area,
macroscopic strain — is governed by first-harmonic algebra, not by simulation. The 2026 paper
computes all of it by forward-kinematics sampling and bisection (Sec. 6, limitation 1:
"such issues can only be detected through explicit deployment simulation via forward kinematics").
That single observation (★) is the engine behind ideas 1, 2, 3, 4 and 10.

---

## Idea 1 — Deployment Jacobian is affine in a fixed matrix

**Type:** theorem + characterization.

**Claim.** For any uniformly deployable periodic kirigami with prescribed periods `P₀ = [t_h, t_v]`,
the macroscopic deployment Jacobian of Eq. (12) is exactly
`J(θ) = cos(θ/2)·I + sin(θ/2)·K`, where `K = P_S P₀⁻¹ ∈ R^{2×2}` is a **constant** matrix and the map
`X ↦ K` is **linear** on the shape space `X` of Eq. (5). Consequently: (a) `J(θ)` is a similarity for
all `θ` iff `K` is a similarity, which is *precisely* the paper's Eq. (13) constraint at `θ = 0` —
so the paper's empirical claim is a two-line theorem; (b) the set of achievable macroscopic
behaviours `{K(X) : X ∈ X}` is an **affine subspace** of `R^{2×2}` of dimension `≤ min(4, dim X)`,
computable by one small linear solve; (c) the Poisson ratio and both principal stretches are closed
form in `θ` for every pattern, conformal or not.

**Gap targeted.** 2026 Sec. 5.1, p. 60:12, verbatim: "Empirically, we observe that for all tilings
with a non-trivial kernel (`E_split ≠ {}`), solving Eq. (13) yields an embedding with conformal
deployment. Moreover, **once the derivative of the conformal distortion vanishes at `θ = 0`, the
deployment remains conformal for all opening angles `θ`.**" My 2026 notes flag this as
"the load-bearing one and it is stated without proof or even a stated experiment count"
(notes §10.1), and list it as unjustified heuristic rows 4 and 5 of §13.2.

**Novelty vs field.**
- *field_kirigami R22 (IsoGami)*: computes expansion profiles (its Fig. 3 polar plots of `U_d`) by
  **sampling** the deployment trajectory with a DAE-avoiding time stepper (its Eqs. 5–9). Numerical
  per-sample; no closed form, no achievable-set characterization, isohedral tilings only.
- *field_kirigami R5/R7 (Konaković 2016/2018)*: assume the triangular linkage is conformal a priori
  and design with conformal maps; they never characterize which linkages are conformal.
- *field_rigidity R17 (Borcea–Streinu geometric auxetics)*: studies periodic auxetic deformations via
  the cone of expansive motions; no per-pattern closed-form Jacobian, and orchestrator's check of
  Borcea–Streinu 2020 (STATE.md) found a single specific 3-DOF framework, not this.

**Why it might be true / sketch.** With (i), the periodic offset between an identified vertex pair is
`p(θ) = y_{i'} − y_i = c(x_{i'} − x_i) + s(S_{i'} − S_i)`, so `P_θ = c P₀ + s P_S` where `P_S` is
linear in `X` and `P₀` is the *prescribed* period matrix of Eq. (3b)/(3c) — constant, hence
`K = P_S P₀⁻¹` is linear in `X`. Then

```
   J(θ) = P_θ P₀⁻¹ = c·I + s·K .
```
Similarity matrices `{aI + bJ}` form a 2D **linear** subspace `Sim ⊂ R^{2×2}`, and `I ∈ Sim`.
Therefore `cI + sK ∈ Sim` for all `θ` ⟺ `K ∈ Sim` ⟺ `K₁₁ = K₂₂` and `K₁₂ = −K₂₁`.
And `dJ/dθ|₀ = ½K`, so the paper's two constraints in Eq. (13) say exactly `K ∈ Sim`. Global
conformality follows. If `K = aI + bJ` then `J(θ) = (c + sa)I + sbJ`, a pure rotation-scale with
scale `√((c+sa)² + (sb)²)` — Poisson ratio `−1` at every `θ`, as the paper asserts. In the general
case `det J(θ) = c² + cs·tr K + s²·det K`, another instance of (★).

**Kill experiment.** `kiri_sweep` over 500 periodic instances: `periodic_squares`,
`periodic_triangles`, `periodic_hexagons`, `periodic_kagome` at sizes 2–6, plus random `σ` flips
kept connected by `describe_orientation`. For each: `assemble_system`/`solve_system` for `X ∈ X`
(random shape-space coefficients), then `deploy()` at 64 angles in `(0, π)`, extract `P_θ` from the
periodic pairs, and compare with `cos(θ/2)I + sin(θ/2)K` where `K` is fitted from `θ = 10⁻³` alone.
**If `max_θ ‖J(θ) − cI − sK‖_∞ > 10⁻⁹` on any instance, the theorem is dead.** Second check: impose
Eq. (13) and verify `‖J₁₁−J₂₂‖ + ‖J₁₂+J₂₁‖ < 10⁻¹⁰` at all 64 angles. Runtime: dominated by 500 ×
64 `deploy()` calls on ≤2,000-face meshes; minutes.

**If it survives, the demo.** (1) A conformality *theorem* replacing the paper's sentence. (2) A plot
over ≥20 patterns of the exact Poisson-ratio curve `ν(θ)` from `K`, overlaid with sampled forward
kinematics — two curves, zero visible gap. (3) The achievable-set figure: for one pattern, the affine
image `K(X) ⊂ R^{2×2}` drawn in the `(K₁₁−K₂₂, K₁₂+K₂₁)` plane, with the origin (= conformal) either
in it or provably not — which **decides** the paper's second empirical claim per pattern. Hero
example: the (3,4,3,12) tiling of Fig. 21, the only named input in the paper.

**Risk.** The periodic pair set `I_h, I_v` may not be uniquely defined for a metatile whose `P₀` the
paper "detects automatically" (Sec. 4.4, procedure not described); if `P₀` itself varies with the
design variables, `K` stops being linear in `X` and part (b) collapses (parts (a), (c) survive).

**Effort.** Derivation S. Code S (all routines exist; needs a `periodicity_jacobian` helper).

---

## Idea 2 — Hole area is a single harmonic; closed-form second closure angle

**Type:** theorem + algorithm.

**Claim.** For a uniformly deployable structure the **total hole area** is
`A(θ) = q(cos θ − 1) + r sin θ` with `q, r` quadratic in `X` (using `A(0) = 0`, the paper's gap-free
assumption). Hence `A(θ) = 0` has exactly one solution in `(0, 2π)`, given in closed form by
`θ_closed = 2·arctan(r/q)`, and for a periodic pattern the same angle is
`tan(θ_closed/2) = tr K / (1 − det K)` in terms of Idea 1's `K` alone. A pattern is
*fully closed at maximal deployment* iff `θ_closed ≤ θ_max` **and** the configuration at `θ_closed`
is overlap-free — replacing the paper's necessary-but-not-sufficient conditions and its
moving-target energy Eq. (14) by one closed-form angle plus one non-overlap test.

**Gap targeted.** 2026 Sec. 5.2, p. 60:12: the paper derives a *necessary* condition ("all hinges
share the same maximal opening angle"), states "However, this condition is not sufficient", adds the
equal-length condition, and offers **no sufficiency proof** (notes §10.2: "no proof is given that the
two conditions together are *sufficient*"). Eq. (14) then optimizes against `ᾱ`, "the average of all
`α_lik`, **updated after each optimization step**" — a block-coordinate scheme with a moving target
and no convergence statement. Also 2025 Observation 3.1 (two maximally open configurations),
`notes/paper_2025.md` §2.10.

**Novelty vs field.**
- *field_kirigami R17 (Dudte, Choi, Mahadevan 2023, additive framework)*: builds compact-to-compact
  quad kirigami by recursive four-bar linkage composition; closure is enforced by construction on
  quads, never as a solvable equation for arbitrary planar graphs.
- *field_kirigami R11 (Choi, Dudte, Mahadevan 2021)*: two compact states linked by a zero-energy
  family — the phenomenon, for quads, without a formula for the second state's angle.
- *field_rigidity R10 (Mitschke et al. 2013)*: finite auxetic deformation of plane tessellations by
  numerical continuation over an archive; no closed-form closure angle.

**Why it might be true / sketch.** Total hole area = (area enclosed by the outer boundary) −
Σ(face areas). Face areas are constant (rigid faces). Each hole's boundary is a fixed cyclic sequence
of duplicated edges — the combinatorics is determined by the hole preimage, not by `θ` — so its
signed area `½ Σ cross(y_i, y_{i+1})` is a quadratic form in `Y_θ`, hence of form (★). `A(0) = 0`
gives `p + q = 0`, so `A(θ) = q(cos θ − 1) + r sin θ = 2 sin(θ/2)[r cos(θ/2) − q sin(θ/2)]`, whose
only root in `(0, 2π)` is `tan(θ/2) = r/q`. Periodic cross-check: hole area per cell is
`det(P₀)·det(cI + sK) − Σ face areas`, so `A(θ) = 0` again ⟺ `c² + cs·tr K + s² det K = 1` ⟺
`s(c·tr K + s(det K − 1)) = 0` ⟺ `tan(θ/2) = tr K/(1 − det K)`. Two independent derivations of the
same number is the cheapest possible self-check. Design: prescribing `θ*` is the single quadratic
equation `r(X) cos(θ*/2) = q(X) sin(θ*/2)` on the affine shape space.

**Kill experiment.** `kiri_sweep` on 300 instances (4 periodic generators × sizes × random `σ`, plus
`voronoi_of_random_points` with fixed boundary): compute `q, r` from the hole cycles returned by
`holes_geometric` at a small `θ`, predict `θ_closed`, then evaluate the true total hole area by
`deploy()` + shoelace on a 200-point grid over `(0, π]`. **If the measured area at the predicted
angle exceeds `10⁻⁸ ·` (total face area) on any instance, or if the sampled area curve is not fit by
`p + q cos θ + r sin θ` to `10⁻¹⁰` relative, dead.** Also cross-check the two formulas agree to
`10⁻¹⁰` on the periodic subset. < 10 min.

**If it survives, the demo.** A gallery of ≥20 patterns annotated with the predicted `θ_closed` and
the measured one (a table with two columns that agree to machine precision), and a designed family:
starting from the (3,4,3,12) tiling, solve the one quadratic for `θ* ∈ {60°, 90°, 120°}` over the
shape space and show three patterns that each close exactly at the requested angle. Baseline: the
paper's Eq. (14) run to convergence, reporting its residual gap area and its iteration count against
our single solve. Export: the `θ*` = 90° pattern as a cut file.

**Risk.** Signed areas cancel under overlap, so `A(θ*) = 0` is only necessary; if for most random
patterns `θ_closed > θ_max`, the closed form is correct but useless as a design tool, and the
contribution shrinks to the harmonic identity. Second risk: hole boundary cycles may change
combinatorially at a collision, invalidating the fixed-cycle premise past `θ_max`.

**Effort.** Derivation S. Code M (needs a hole-area accumulator over `holes_geometric` cycles).

---

## Idea 3 — All orientation predicates are single harmonics: exact contact angles

**Type:** characterization.

**Claim.** For any triple of `M'`-vertices, the orientation determinant along the uniform deployment
is `D(θ) = p + q cos θ + r sin θ`; therefore every vertex–edge and edge–edge contact event during
deployment occurs at a root of an explicit single harmonic, and `θ_max` is
`min` over a finite candidate set of closed-form roots — **no bisection, no time stepping, and exact
to machine precision**. Corollary to be tested and, I predict, **falsified**: the paper's assumption
that "collisions occur only along the split-cut edges" (Sec. 4.5) fails on a measurable fraction of
random graphs.

**Gap targeted.** 2026 Sec. 6 limitation 1 verbatim: "At present, such issues can only be detected
through explicit deployment simulation via forward kinematics and require additional post-processing
to resolve. **Developing geometric characteristics for favorable deployment behavior directly from
the embedding remains an open problem.**" And Sec. 4.5's empirical restriction: "we observe in all
experimental cases that collisions occur along the split-cut edges" (assumption 17 of notes §12).

**Novelty vs field.**
- *field_kirigami R18 (Liu et al. 2024)*: the only closed-form collision inequality I found,
  `length(C1p1) + length(C3p2) < 2a cos θ`, but only for enumerated opposite/adjacent quad pair types
  with curved contours, and they fall back to a manual UI for the general case.
- *field_kirigami R22 (IsoGami)*: handles collisions with IPC log-barriers inside the continuation
  (its Eq. 9) and halves the step size up to 10 times to find the deployment limit — i.e. exactly the
  bisection this idea removes, and only for the one-ring of a unit cell.
- *field_tutte row 13/14 (Streinu parallel redrawing; Spriggs)*: Remark A.4 makes deployment a
  parallel morph; this literature asks *whether* a parallel morph stays planar but supplies no
  per-pair event time.

**Why it might be true / sketch.** With `y_i(θ) = c·C_i + s·S_i`,
`det(y_j−y_i, y_k−y_i) = c²·det(ΔC,ΔC') + cs·[det(ΔC,ΔS') + det(ΔS,ΔC')] + s²·det(ΔS,ΔS')`, which is
(★). Face–face overlap of convex pieces begins exactly when some vertex–edge orientation predicate
changes sign, so the first contact angle is
`θ_max = min { smallest root in (0,π] of D_t : t ∈ T }` over the candidate triples `T`. `T` need not
be all `O(n³)` triples: hinge geometry restricts contact to pairs of faces that are hinge- or
split-adjacent, or become adjacent across a hole, so `|T| = O(n)` after a one-ring filter — and the
filter itself is checkable by comparing against the existing brute-force `theta_max()`.

**Kill experiment.** 500 graphs from `delaunay_of_random_points`, `voronoi_of_random_points`,
`quad_dominant_random` (n = 60…400) with `σ` from `assign_orientation_relaxation`. For each: compute
`theta_max()` (existing grid + bisection, grid = 720) and the closed-form minimum root over the
candidate set. **If the two disagree by more than `10⁻⁶` rad on any instance — in particular if the
closed form ever reports a *larger* angle, meaning a missed candidate pair — the characterization is
dead.** Separately, log for every first-contact event whether the contacting pair involves duplicates
of a split edge; **if that fraction is 100% over all 500 graphs, the corollary is dead** (and the
paper's assumption stands, which is itself a publishable negative). < 20 min.

**If it survives, the demo.** A speed/accuracy table on ≥20 graphs: bisection `θ_max` vs closed form,
error and wall-clock (expect 2–3 orders of magnitude). A counterexample figure: the smallest graph
where the first collision is *not* at a split-edge pair, with the offending pair highlighted, plus the
fraction over the 500-graph sample. The corrected barrier set for Eq. (9) then follows for free.

**Risk.** Contact between two faces meeting only at a hinge vertex is a degenerate (tangential)
predicate root — a double root of `D` — so the "smallest root" rule may need a multiplicity filter,
and shrink tolerances like the `shrink = 1e-6` already used in `polygons_overlap` may be doing real
work. The whole thing is also only as true as fact (i).

**Effort.** Derivation S. Code M (candidate-pair enumeration is the whole job).

---

## Idea 4 — Range-maximizing embedding with a certified upper bound

**Type:** algorithm.

**Claim.** Maximizing the collision-free deployment range over the shape space,
`max_{T} θ_max(X₀ + ΦT)`, is exactly a max–min over roots of the single harmonics of Idea 3, whose
coefficients `(p,q,r)` are **quadratic forms in `T`**. Requiring all predicates to keep their sign on
`[0, Θ]` is therefore a system of quadratic inequalities in `T`; lifting `M = TTᵀ ⪰ 0` gives a
semidefinite relaxation whose optimum is a **certified upper bound** on the achievable range, and the
rank-one rounding of that relaxation beats the paper's Eq. (9) barrier-at-`θ=0` heuristic by a large
margin on ≥20 graphs.

**Gap targeted.** 2026 Eq. (9) is optimized "at `θ = 0`, i.e. on the *undeployed* configuration", and
my notes (§8.5) record: "Optimizing at `θ = 0` only is itself a heuristic: it penalizes the *initial*
velocity direction at split edges, not collisions at the angle where they actually occur. The paper
does not justify why this suffices". Open question 14 of notes §15: "Does optimizing Eq. (9) at
`θ = 0` actually maximize `θ_max`, or merely delay the first collision?" Nobody has ever reported an
upper bound on the attainable range for a given cut structure.

**Novelty vs field.**
- *field_kirigami R22 (IsoGami)*: maximizes expansion by **sampling** tile-shape parameters on a
  coarse grid ("This process is largely heuristic and could benefit from further automation", its
  Sec. 4.2). No objective, no bound.
- *field_kirigami R18 (Liu et al. 2024)*: closed-form range test for enumerated quad pairs, used as a
  feasibility filter, not as an objective with a certificate.
- *field_kirigami R13 (Chen et al. 2021)*: optimizes elastic auxetic cells from a library; different
  model, no rigid-face range guarantee.

**Why it might be true / sketch.** For a fixed candidate pair `t`, sign preservation of
`D_t(θ) = p_t + q_t cos θ + r_t sin θ` on `[0,Θ]` is implied by the three conditions
`D_t(0) > 0`, `D_t(Θ) > 0`, and either `q_t sin θ* − r_t cos θ* ≠ 0` inside or
`D_t(θ*) > 0` at the single interior extremum `θ* = atan2(r_t, q_t)`. Each is affine in `(p,q,r)`,
hence a quadratic form `Tᵀ Q_t(θ) T + linear + const ≥ 0`. Maximizing `Θ` subject to all of them is a
QCQP; the standard Shor lift (`M = TTᵀ`, drop rank) is an SDP in `dim(X)²` variables — and `dim(X)`
is *small*, 6 or 8 in the paper's own examples (Figs. 12, 13), so the SDP is tiny. A bisection on `Θ`
turns it into a feasibility SDP. Since only rank ≤ `dim X` matters, Burer–Monteiro factorization
`M = VVᵀ` with the existing `lbfgs_minimize` is enough; the dual value certifies the bound.

**Kill experiment.** 20 graphs where the Eq. (9) baseline is known to matter (the `Fig. 14`-style
cases: `voronoi_of_random_points` with `σ` giving `|E_split| ≥ 3`). For each: (a) `optimize_collision`
(existing Eq. (9) implementation) → `θ_max^base`; (b) the SDP bisection → `θ_max^ours` and the dual
bound `θ̄`. **If `θ_max^ours ≤ θ_max^base` on more than 5 of the 20, or if the certified bound `θ̄` is
more than 3× the best achieved value on every instance (i.e. the relaxation is vacuous), dead.**
Runtime: 20 × (≈20 bisection steps × a `dim X`-sized Burer–Monteiro solve) — seconds each.

**If it survives, the demo.** A scatter over ≥20 graphs of `θ_max` baseline vs ours vs the certified
upper bound (three columns, one row per graph): the story is "ours is close to the ceiling, the
baseline is not, and now we know where the ceiling is". Hero: one pattern taken from
`θ_max ≈ 20°` to `θ_max ≈ 180°` (the paper's Fig. 14 claims exactly such a jump for its own method,
so it is the right thing to beat). Export: the maximal-range pattern as a cut file.

**Risk.** The relaxation may be uniformly loose because the constraint set is a *conjunction over
many pairs*, and Shor bounds degrade with the number of nonconvex constraints; then the "certified
bound" contribution dies and only a better local optimizer remains, which is a much weaker paper.
Second risk: `dim X = 0` for most random graphs (pure hinge case), leaving nothing to optimize over —
must pre-screen with `solve_system` for `dim_null ≥ 2`.

**Effort.** Derivation M. Code L (Burer–Monteiro + duality certificate; no external SDP solver is
allowed by D3, so this is written from scratch on top of `lbfgs_minimize`).

---

## Idea 5 — Split cuts cancel residuals: max-cut is the wrong objective for σ

**Type:** algorithm + counterexample.

**Claim.** Define the **deployability defect** of an orientation on the *given* geometry,
`D(σ) = Σ_{C ∈ C^M(σ)} ‖ Σ_{e ∈ C ∩ E_hinge} e ‖²`. Then `D` is **not monotone** in `|E_split|`:
because a split edge merges two hole preimages and replaces their residuals `ρ₁, ρ₂` by `ρ₁ + ρ₂`,
adding split cuts can strictly *decrease* the defect, and there exist planar graphs where every
minimum-`|E_split|` orientation has `D > 0` while some orientation with more split cuts has `D = 0`.
Therefore the paper's max-cut objective is the wrong proxy, and a defect-driven local search over `σ`
beats it by a large margin on the paper's own use case (Eq. (6) repair distance).

**Gap targeted.** 2026 Sec. 4.2, p. 60:6: "we empirically favor many balanced, small holes rather
than a few large ones. This motivates choosing face orientations that minimize the number of split
edges". My notes §13.2 record the objection the paper never addresses: "Note the tension: fewer split
edges also means *smaller kernel*, i.e. a smaller design space. The paper does not discuss this
trade-off." And Sec. 4.2 chooses `σ` **before** ever looking at the input geometry `X_ini`, although
the downstream deliverable (Eq. (6)) is precisely "the closest embedding to the initial guess".

**Novelty vs field.**
- *field_kirigami R22 (IsoGami)*: the direct threat — joint topology + shape exploration. But it
  *enumerates* `3^m` joint choices per unit cell with `m ≤ 16` and filters hierarchically; its
  objective is 1-DOF plus expansion, not a geometric defect, and it never treats topology choice as
  an optimization against a target embedding.
- *field_kirigami R9 (Chen, Choi, Mahadevan 2020)*: varies kirigami topology at fixed geometry, by
  deterministic/stochastic control of cut connectivity; no defect objective, quads only.
- *field_kirigami R17 (Dudte et al. 2023)*: varies geometry at fixed topology. Nobody couples the two
  through the residual of Eq. (2).

**Why it might be true / sketch.** By F11/F14, `L = R·D` where `R` is the partition matrix of hole
preimages; flipping a face changes which edges are hinge vs split and thus merges or splits rows of
`R`. Merging rows `K₁, K₂` replaces the pair of constraints `ρ₁ = 0, ρ₂ = 0` by the single weaker
`ρ₁ + ρ₂ = 0`. Hence `D` can only drop when a merge puts two residuals that are near-negatives in the
same group, and it strictly drops iff `ρ₁ ≠ 0` and `ρ₁ + ρ₂ = 0`. Constructing an instance: take any
2-colourable tiling, perturb one interior vertex so two adjacent holes acquire residuals `±ρ`, and
flip the face between them so the connecting edge becomes a split. The minimum-`|E_split|` solution
(the pure 2-colouring) has `D = 2‖ρ‖²`; the flipped one has `D = 0`. This is a *construction*, so the
existence half is nearly certain; what is open is how often it matters on random inputs and how much
of the Eq. (6) repair distance it saves.

**Kill experiment.** (a) Existence: build the perturbed-square construction above by hand
(`periodic_squares(3,3)` + one vertex moved), enumerate all `2^F` orientations with
`brute_force_orientation`-style enumeration, and check that `argmin |E_split| ≠ argmin D`.
**If they always coincide, the whole idea is dead in one minute.** (b) Significance: 40 graphs with
`F ≤ 18` from `quad_dominant_random` and `voronoi_of_random_points`; enumerate all `2^F` `σ`, compute
for each the Eq. (6) repair distance `ρ(σ) = ‖X₀(σ) − X_ini‖` via `solve_system`, and compare
`ρ(σ_GW)` from `assign_orientation_relaxation` against `min_σ ρ(σ)`. **If the GW orientation is within
5% of optimal on every one of the 40 graphs, the algorithmic half is dead.** `2^18 × ` a small dense
SVD is the cost driver — restrict to `F ≤ 16` and use `rank_only_sparse` for screening; < 30 min.

**If it survives, the demo.** Baseline comparison on ≥20 graphs: repair distance and `|E_split|` for
GW-`σ` vs defect-optimized `σ` (greedy face flips with incremental re-partitioning, warm-started).
The headline figure is a scatter of `|E_split|` (x) against repair distance (y) showing they are
uncorrelated or anti-correlated — that single plot refutes the stated principle. Hero: a tiling that
is exactly deployable for a many-split `σ` and irreparably distorted for the minimum-split one.

**Risk.** The counterexample may be non-generic: perhaps for *most* geometries the residuals do not
line up and max-cut is within noise of optimal, leaving a true but unimportant result. Also `σ` must
keep `M'` connected (Sec. 4.2 principle 1), which the enumeration must enforce or the comparison is
unfair.

**Effort.** Derivation S. Code M (enumeration + incremental flip update).

---

## Idea 6 — The dual is planar, so max-cut is not NP-hard here

**Type:** algorithm.

**Claim.** The dual graph `M_d` of a planar `M` is planar, and max-cut on planar graphs is solvable
**exactly in polynomial time** (Orlova–Dorfman 1972; Hadlock 1975, via minimum `T`-join / matching).
So the paper's justification — "Since max-cut is NP-hard, we use the following approximation, based
on the algorithm proposed by Goemans and Williamson" — does not apply to its own instances, and the
orientation assignment it solves heuristically has an exact polynomial algorithm. Quantified claim:
the paper's unit-circle relaxation plus diameter rounding leaves a measurable optimality gap in
`|E_split|` on ≥20 graphs.

**Gap targeted.** 2026 Sec. 4.2, p. 60:6, verbatim as above. My notes §6.5 already record that the
method "is not actually the GW algorithm": "Here the vectors are constrained to the **unit circle in
`R²`**, which is a non-convex problem, not an SDP, and carries no approximation guarantee", and
"How the diameter is chosen is not stated."

**Novelty vs field.** This is not a new algorithm in combinatorial optimization — Hadlock is 1975 —
and I say so plainly. The novelty is (a) the observation that the paper's stated complexity barrier
is inapplicable to its own inputs, and (b) the measurement of what the heuristic costs, which nobody
has reported. *field_kirigami R22 (IsoGami)* enumerates `3^m` joint topologies exhaustively for
`m ≤ 16` precisely because it has no such structure to exploit; *field_kirigami R9* searches topology
stochastically; *field_tutte* has no row on max-cut at all (I checked rows 11, 12 — Whiteley matroid
and Servatius–Whiteley directions/lengths — which are about rank, not cuts).

**Why it might be true / sketch.** `M_d` is the planar dual restricted to interior adjacencies, hence
planar (this is also `field_rigidity` Sec. 0 fact 1: `Γ` is planar). For planar graphs, a cut
corresponds to a cycle in the dual of `M_d` — i.e. in `M` itself — and max-cut = maximum-weight cycle
= complement of a minimum-weight `T`-join over odd-degree terminals, computed by shortest paths plus
a minimum-weight perfect matching on the terminals. All weights here are `1`. Note `M_d`'s own dual
is `M`, so the machinery runs on the input graph directly. Two honesty flags: **I did not verify the
Hadlock reference myself in this session** — treat "planar max-cut ∈ P" as a claim to check before
committing — and the connectivity side condition (principle 1 of Sec. 4.2) is *not* part of max-cut,
so the exact optimum may need a repair pass.

**Kill experiment.** 200 graphs with `F ≤ 20` (`voronoi_of_random_points`, `quad_dominant_random`,
`tiling_snub_square` patches): compute `|E_split|` from `assign_orientation_relaxation` (with its
existing multi-restart + best-diameter rounding) and from `brute_force_orientation` (exact by
enumeration over `2^F`, already implemented for `F ≤ 20`). **If the relaxation attains the exact
optimum on ≥ 99% of instances and is never worse than 1 split edge, the algorithmic contribution is
dead** (the heuristic is fine and only the complexity remark survives, which is not a paper).
< 15 min, uses only existing routines.

**If it survives, the demo.** The exact `T`-join solver on ≥20 large graphs (hundreds to thousands of
faces, where brute force is impossible), reporting `|E_split|` heuristic vs exact and the
corresponding `dim X` and repair distance. Paired with Idea 5 this becomes one narrative: *solve the
paper's stated problem exactly, then show the stated problem is the wrong one.*

**Risk.** The kill experiment plausibly *passes for the paper*: with multiple restarts and a best-of
diameter rounding, the circle relaxation may be optimal on almost all planar duals, which are sparse
and nearly bipartite. Second risk: a blossom implementation is a substantial piece of C++ with no
external library allowed.

**Effort.** Derivation S. Code L (minimum-weight perfect matching from scratch).

---

## Idea 7 — The SVD basis is the wrong basis: a sparsest design-space basis

**Type:** characterization + algorithm.

**Claim.** `x ∈ ker L` iff the edge weighting `w = Dx` (the per-hinge-edge coordinate difference)
sums to zero on every hole preimage; consequently `ker L` admits a basis whose vectors are supported
on **explicit local vertex sets** (one generator per split edge, plus `dim Z` global generators), and
a matroid-greedy over these candidate supports returns a **provably sparsest** basis in
`O(|E| α(|V|))`. The paper's own interpretation of Eq. (5) — "Each `φ_i` identifies a set of vertices
that can be translated together" — is **false for the basis it actually computes**, since SVD returns
dense orthonormal vectors, not indicators.

**Gap targeted.** 2026 Sec. 4.4, p. 60:8: "We can perform a singular value decomposition (SVD) of the
system matrix in Eq. (4) to identify the null space, collecting the vectors `{φ_i}` corresponding to
zero singular values. **Each `φ_i` identifies a set of vertices that can be translated
collectively.**" And Sec. 5.4: the UI exposes these as sliders. Notes §7.4 flags the numerical-rank
problem: "no threshold for 'zero singular value' is given ... For an implementer this is a real gap".

**Novelty vs field.**
- *field_tutte row 15 (incidence-matrix rank) + F14*: `L = R·D` gives `rank(L) = H − dim Z`; that is
  the *dimension*. Nobody in that table gives a **basis**, sparse or otherwise.
- *field_tutte row 7 (Colin de Verdière null space)*: the only place a null space of a planar-graph
  matrix is given combinatorial meaning, via nodal domains — but for a symmetric CdV matrix, and `L`
  is neither symmetric nor vertex-indexed.
- *field_kirigami R17 (additive framework)*: its recursion says the constraint matrix is triangular
  in a good order — the same flavour of argument, for quads only, and it produces a construction, not
  a basis of a design space.
- Sparse null-space basis is a classical NP-hard problem (Coleman–Pothen); the content here is that
  *this particular* `L` has enough structure to make it exactly solvable.

**Why it might be true / sketch.** `Lx = R(Dx) = 0` says: for each preimage `K`,
`Σ_{e ∈ C_K ∩ E_hinge} (x_{head} − x_{tail}) = 0`. In the pure case (`E_split = ∅`) this is one
constraint per interior vertex and the kernel is trivial with a fixed boundary. Each split edge
merges two preimages, deleting one constraint, so `dim ker` grows by exactly one per split edge
(Scout-c's (2.2): `dim = 2(|E_split| + dim Z)`). The natural generator for the split edge `s` is the
kernel vector obtained by solving the *remaining* constraints with a unit source at `s`, whose support
should be the set of vertices reachable from `s` before the constraints re-pin the solution. Matroid
structure: kernel vectors of a matrix form the cocircuit space of a linear matroid, and greedy on
sorted candidate supports gives the sparsest basis whenever the candidates contain a basis
(Coleman–Pothen's condition). Falsifiable both ways.

**Kill experiment.** 500 graphs with `dim_null ≥ 2` (screen with `solve_system`). For each: (a) verify
the local-generator construction spans the same space as `Phi` from the SVD (compare projectors,
`‖P_greedy − P_svd‖_F < 10⁻⁸`); (b) compare greedy sparsity `Σ‖φ_i‖₀` against an exhaustive minimum
over all support sets of size ≤ 6 on the subset with `N ≤ 40`. **If the greedy basis fails to span
the kernel on any instance, or is beaten by exhaustive search on more than 5% of the small ones, the
"provably sparsest" claim is dead** (a heuristic sparse basis would still be a minor result).
< 20 min.

**If it survives, the demo.** Side-by-side figure: the paper's SVD basis vectors rendered as
per-vertex heat maps (dense, global, uninterpretable) against our basis (compact, disjoint, one blob
per split edge) — this is the figure that shows the paper's sentence is about our basis, not theirs.
Plus a UI-relevant number over ≥20 graphs: mean support size, SVD vs ours.

**Risk.** The kernel may simply not have a local basis — a single split edge in a periodic pattern can
produce a globally supported mode (the `dim Z = 1` translation direction certainly is global). Then
the claim survives only in the fixed-boundary case, which is the weaker half.

**Effort.** Derivation M. Code S–M.

---

## Idea 8 — The exact non-uniform deployment variety is a linear section of a torus

**Type:** characterization.

**Claim.** Identify `R²` with `C`. Assign to each face `f` the complex number `z_f = e^{iψ_f}`
(`ψ_f` = its rotation angle). Then the **complete, finite (not infinitesimal)** configuration space of
the cut structure is

```
      Def(M') = { z ∈ C^{|F|} :  W z = 0 ,  |z_f| = 1 ∀f } ,
      W ∈ C^{b₁(Γ) × |F|},   W[γ, ·] = Σ_{e ∈ γ} ± x_{u_e}  (hinge vertex positions, complex),
```
i.e. **complex-linear constraints intersected with a torus**. `W·1 = 0` always (global rotation), and
`W·σ = 0` is exactly Eq. (2); the uniform deployment path is the circle
`z = cos(θ/2)·1 − i sin(θ/2)·σ` inside `span_C{1, σ}`. Extra finite deployment modes exist iff
`dim_C ker W > 1` *and* the extra kernel meets the torus, both decidable.

**Gap targeted.** 2026 Sec. 6, limitation 2, verbatim: "Second, patterns that can be uniformly
deployed **might also be deployed with a different set of hinge angles. Analyzing the degrees of
freedom of the deployment space, which could be achieved by sensitivity analysis, is an interesting
direction to explore.**" Also 2025 Sec. 6 / supplement C (non-uniform deployment "not given" in the
main text, `notes/paper_2025.md` §14).

**Novelty vs field.**
- *fact (ii), from ideator-rigidity*: mobility `= dim ker A − c(Γ)` with `A` real `2b₁(Γ)×|F|`. That
  is the **linearization** of `W` at a configuration; my `W` is the exact finite object, and the
  complexification collapses `A`'s `2b₁ × |F|` real system into a `b₁ × |F|` complex one whose kernel
  is a *complex* subspace — a strictly stronger statement, because it says the constraint manifold is
  linear in `z` and all the nonlinearity sits in `|z_f| = 1`.
- *field_kirigami R22 (IsoGami)*: computes mobility as the count of zero eigenvalues of `JᵀJ` at
  each sampled configuration and continues numerically through turning points. Per-configuration
  numerics on isohedral tilings; no algebraic model of the variety.
- *field_rigidity R25 (Li, Zhu & Qu, rigid-origami mobility to 4th order)* and *R20 (Müller)*: local
  higher-order mobility machinery; both are series expansions at a point, not a global variety.
- *field_rigidity R14 (Kapovich–Millson)*: planar linkage moduli can be arbitrary varieties — the
  upper bound on what any theorem can claim, and the reason a *linear-section-of-a-torus* answer for
  this family is surprising.

**Why it might be true / sketch.** A pin joint at hinge vertex `u` between faces `f, g` reads
`R_f x_u + t_f = R_g x_u + t_g`, so `t_g − t_f = (R_f − R_g)x_u`. Consistency of the translations
around any cycle `γ` of the hinge graph `Γ` is `Σ_{e ∈ γ}(R_{f_e} − R_{g_e})x_{u_e} = 0`, and there
are `b₁(Γ)` independent cycles. In complex form `R_f ↦ z_f`, so each cycle gives one **complex-linear**
equation in `z`. Sanity checks that must hold: `z ≡ 1` (all faces unrotated) is always a solution;
`z = σ` gives `Σ(σ_f − σ_g)x_u = 2Σ σ_f x_u`, which telescopes to the directed-edge sum of Eq. (2);
and `z = c·1 − i s·σ` reproduces fact (i) exactly. Deployment paths are then curves in
`ker_C W ∩ (S¹)^{|F|}`, a system of `|F|` real quadratics in `2·dim_C ker W` unknowns — tiny, and
solvable by Newton from many starts with certification for small graphs.

**Kill experiment.** 300 graphs (all generators, `F ≤ 300`). Assemble `W` from a cycle basis of `Γ`
(spanning tree + non-tree hinge edges). Checks: (1) `‖W·1‖ < 10⁻¹²`; (2) `‖W·σ‖` equals the stacked
Eq. (2) residual from `hole_residuals` up to the known factor — **if these two identities fail, the
model is wrong and the idea is dead in five minutes**; (3) for `X` solved to be uniformly deployable,
`W(c·1 − i s·σ) = 0` for 32 angles; (4) reconstruct `Y_θ` from `z` and compare against `deploy()`
vertex-by-vertex to `10⁻¹⁰`. (5) On the subset with `dim_C ker W ≥ 2`, run Newton from 200 random
torus starts and report whether any solution outside `span_C{1,σ}` is found. < 25 min.

**If it survives, the demo.** A table over ≥20 patterns of `dim_C ker W`, the extra-mode count found
by the torus solve, and whether the extra modes are collision-free; a hero pattern with a **second,
non-uniform** deployment path animated beside the uniform one; and a design experiment: pick `X` in
the shape space to force `dim_C ker W = 2` (a rank-drop condition on a small matrix, linear-algebraic
because `W` is linear in `X`) and exhibit the designed extra mode.

**Risk.** The cycle constraints may be equivalent to the hole-preimage constraints in a way that makes
`ker_C W = span{1, σ}` for essentially every pattern, making the variety a single circle and the
result a (nice) rigidity theorem rather than a design tool. Also `W` linear in `X` and `z`-linear is
the whole strength — if a sign or a hinge-vertex convention is wrong, check (2) will expose it.

**Effort.** Derivation M. Code M.

---

## Idea 9 — Certified inverse design by low-rank relaxation on the shape space

**Type:** algorithm.

**Claim.** Because `Y_θ` is affine in `X` (fact (i)) and the shape space is affine (`X = X₀ + ΦT`),
the 2026 inverse-design step and the fully-closed objective are both of the form
"minimize `Σ_e (‖A_e(X)‖² − ℓ_e²)²` over an affine set" — a low-rank QCQP in `T` whose Shor/
Burer–Monteiro relaxation gives a **lower bound and a global-optimality certificate** when the
recovered matrix has rank 2. Claim: the relaxation is tight (rank 2 to `10⁻⁶`) on a majority of
instances, and where it is not, its bound exposes the paper's Newton scheme as landing in
strictly suboptimal local minima on ≥20 instances.

**Gap targeted.** 2026 Sec. 5.3, p. 60:12–13: "We then adjust `Y_θ` so that edge lengths in `M'`
correspond to those in `S`" — my notes §10.3 record: "the edge-length matching energy is not written
down; the choice of `θ` is not discussed". And Sec. 5.4: "We optimized all the non-linear energies
via Newton's method with a per-element Hessian projection onto the PSD cone" — a local method with
no global statement anywhere in either paper. 2025 Eqs. (1)–(2) are the same story with weights
`ω₁..ω₄` that "never appear in the main text" (`notes/paper_2025.md` §3).

**Novelty vs field.**
- *field_kirigami R22 (IsoGami)*: local penalty minimization with IPC barriers; no bound.
- *field_kirigami R5/R7 (Konaković 2016/2018)*: conformal-map initialization then local
  optimization; no certificate.
- *field_kirigami R23 (Jiang & Choi 2026)*: multi-state length-based optimization on quads; also
  local. No one in this literature reports a global-optimality certificate for a kirigami inverse
  design, which is the entire point.

**Why it might be true / sketch.** Let `u = [T; 1]`. Each edge-length residual is
`(uᵀ B_e u − ℓ_e²)`, with `B_e` PSD of rank ≤ 2 because `A_e` maps into `R²`. Lifting `M = uuᵀ ⪰ 0`
with `M_{last,last} = 1` linearizes every residual: `tr(B_e M) − ℓ_e²`. The objective becomes a convex
quadratic in `M` over a spectrahedron; dropping the rank constraint is the relaxation. `dim T` is
small (6–8 in the paper's figures, ≤ ~100 in its largest case per Sec. 5.4 timings), so `M` is at
worst `101 × 101` — trivial. Burer–Monteiro `M = VVᵀ` with `V ∈ R^{(dimT+1)×r}` and the existing
`lbfgs_minimize` finds the relaxation optimum; a rank-2 result plus a dual residual is the
certificate.

**Kill experiment.** 20 instances: `tiling_3_4_3_12` and 4 periodic generators, each with target edge
lengths sampled from a perturbed copy of the pattern (so a solution is known to exist and the true
optimum is 0). Run (a) Newton/L-BFGS from 10 random starts (the baseline stand-in for the paper's
method), (b) Burer–Monteiro at `r = 2, 3, 4`. **If the relaxation's optimum is strictly below every
Newton run's value on fewer than 3 of the 20 instances, and never certifies rank 2, dead.**
Instances are tiny; < 10 min.

**If it survives, the demo.** A table over ≥20 instances: best-of-10 Newton objective, relaxation
bound, certified-optimal flag, wall clock. The hero is an instance where 10/10 Newton restarts land
above the certified optimum, i.e. the paper's solver provably misses. Export: the certified-optimal
inverse-designed pattern.

**Risk.** The lift squares the problem size in the *shape-space* dimension only if `dim X` stays
small; the paper's Sec. 5.3 pipeline has boundary vertices fixed to a circle, which may leave
`dim X = 0` and nothing to relax over — in which case the idea applies to the fully-closed objective
(Eq. 14) but not to inverse design, halving its scope. Also the "known solution exists" instance
design makes the test easy; harder targets may be uniformly loose.

**Effort.** Derivation M. Code L.

---

## Idea 10 — "Small balanced holes" is not the right proxy, and here is one that is

**Type:** characterization (with an algorithmic payoff).

**Claim.** Over a large random sample, hole size and balance have **no significant correlation** with
any operational measure of auxetic quality, whereas two computable quantities derived from Idea 1 do:
the isotropy defect `d_iso(X) = |K₁₁−K₂₂| + |K₁₂+K₂₁|` and the dimension of the achievable set
`dim{K(X) : X ∈ X}`. Replacing "minimize `|E_split|`, keep holes small and balanced" by
"minimize `d_iso` over `σ` and `X` jointly" changes the selected orientation on a substantial
fraction of graphs and improves the delivered conformality by a large margin.

**Gap targeted.** 2026 Sec. 4.2, p. 60:6: "For auxetic behaviors, **we empirically favor many
balanced, small holes rather than a few large ones.**" My notes §13.2 record the support offered:
"Nothing. No definition of 'balanced', no experiment." This is the only place in either paper where a
design principle is asserted with zero evidence.

**Novelty vs field.**
- *field_kirigami R2 (Mitschke et al. 2013)*: screens an archive of 1- and 2-uniform tessellations for
  auxetic behaviour — the closest thing to a correlation study, but over a fixed archive with a
  bar-and-joint model, and no per-pattern predictor.
- *field_kirigami R22 (IsoGami)*: characterizes topology diversity in its dataset and reports
  expansion profiles, but its filters are thresholds on expansion, not tested predictors.
- *field_kirigami R24 (Dang & Paulino 2025 survey)*: no such correlation reported.

**Why it might be true / sketch.** From Idea 1, `K` is linear in `X` and captures the *entire*
macroscopic behaviour; hole size enters `K` only through the particular geometry, not through the
count of split edges. Two patterns with the same `|E_split|` can have `d_iso` differing by orders of
magnitude, and `dim{K(X)}` — the count of macroscopic behaviours reachable without changing the cut —
is bounded by `min(4, dim X)` and is what actually decides whether Eq. (13) can succeed. So the
hypothesis is that the paper's proxy is correlated with the wrong thing (kernel *size*), and Idea 5
already predicts the sign of the error.

**Kill experiment.** 500 periodic instances with `σ` sampled uniformly among connected assignments.
For each record: `|E_split|`, mean and variance of hole sizes (from `holes_partition`), `dim X`
(`solve_system`), `d_iso`, `dim{K(X)}`, and three outcome measures — `θ_max` (`theta_max`), the
achieved conformal residual after Eq. (13), and the area expansion at `θ_max`. Compute Spearman
correlations. **If hole-size variance correlates with the outcome measures at `|ρ| > 0.5` while
`d_iso` and `dim{K}` do not, the claim is dead and the paper's heuristic is vindicated.**
< 30 min; correlation plots in matplotlib from the C++ CSV dump.

**If it survives, the demo.** A correlation matrix figure (the paper's proxy vs ours vs outcomes) and
a head-to-head on ≥20 graphs: `σ` chosen by max-cut vs `σ` chosen by minimizing `d_iso`, with the
resulting conformal residual per graph.

**Risk.** This is the weakest idea on the list: a correlation study is not a theorem, and a null
result ("nothing predicts anything") is likely and unpublishable on its own. It earns its place only
as a section inside Idea 1 or Idea 5, not as a standalone contribution.

**Effort.** Derivation S. Code S (all measurements exist; this is a sweep driver).

---

## Self-attack — hostile SIGGRAPH/TOG reviewer on my top 3

**Against Idea 1 (`J(θ) = cos(θ/2)I + sin(θ/2)K`). Reject.**
The paper's entire mathematical content is two lines of algebra that follow immediately from the
kinematic model the authors already state. Segall et al. write that faces rotate rigidly by `±θ/2`;
anyone who writes that sentence and then differentiates has this result. That the original authors
called it an "empirical observation" is a matter of exposition, not of knowledge — they clearly
computed `∂J/∂θ` at `θ = 0` to *derive* Eq. (13), so they had the linearity in front of them.
Formalizing a one-line observation as a theorem is a note, not a TOG paper. The "achievable set is an
affine subspace" claim is likewise definitional once you accept that the shape space is affine and
`P₀` is fixed; and it is fragile, since the authors explicitly allow `P₀` to be detected from the
input, in which case it is not fixed. Finally the deliverable is a plot showing two curves that
coincide — a validation figure, not a contribution. Novelty gate: a reader of Sec. 5.1 who wanted the
proof could produce it in an afternoon.

**Against Idea 8 (complex torus variety). Reject.**
This is a change of notation dressed as a theorem. Writing planar rotations as unit complex numbers
and observing that pin-joint loop closures become linear in those complex numbers is standard in the
linkage literature — it is the classical Denavit–Hartenberg / complex-vector loop-equation method
that mechanism design has used since the 1950s (see the mechanism-design texts cited by IsoGami
itself), and it is exactly how planar four-bar analysis is taught. The "surprising" statement that
the configuration space is a linear section of a torus is the *definition* of a planar linkage's
loop equations, and Kapovich–Millson tells us such sections realize arbitrary varieties, so the
structure buys no tractability. The author then proposes to solve `|F|` real quadratics by Newton
"from many starts" — that is precisely the numerical mobility computation IsoGami already performs
(zero eigenvalues of `JᵀJ`, continuation through turning points), with worse robustness and no
contact handling. And the concrete deliverable, "design a pattern with a second finite mode", is a
rank-drop condition whose feasibility is not established anywhere in the submission.

**Against Idea 2 (hole area is a single harmonic). Reject.**
The result as stated is nearly vacuous: `A(θ) = 0` is a *necessary* condition for full closure, and
the submission concedes that signed areas cancel under overlap, so the closed form certifies nothing
without the very collision test the authors claim to have eliminated. Worse, the second root
`θ_closed = 2 arctan(r/q)` is meaningful only when it precedes the collision angle, and the paper
gives no argument, and no measurement, that this happens for any interesting fraction of patterns —
the risk section admits as much. The relationship to Segall Sec. 5.2 is also overstated: their Eq.
(14) is not an attempt to *find* the closure angle, it is an attempt to *make the pattern gap-free
there*, i.e. to enforce equal `β` and equal contact lengths, which this submission does not address
at all. So the claimed replacement of Eq. (14) is not a replacement. Two derivations of the same
formula agreeing is a debugging practice, not evidence of significance. For a graphics venue the
demo — three patterns that close at 60/90/120 degrees — needs the fabricated artifact to be
convincing, and none is shown.

*(What I take from my own attacks: Idea 1 must ship together with the achievable-set decision
procedure and the exact `ν(θ)` curves for non-conformal patterns, or it is a note. Idea 8 must be
positioned against the loop-equation literature explicitly — `notes/field_rigidity.md` has no row for
complex loop equations, which is a gap in my prior-art coverage, not evidence of novelty. Idea 2 must
be merged with Idea 3, so that closure and collision are decided by the same harmonic machinery.)*

---

## Ranked list

| # | Idea | One-line justification |
|---|---|---|
| 1 | **1 — `J(θ) = cos(θ/2)I + sin(θ/2)K`** | Highest confidence, kills a load-bearing empirical claim with a proof, and hands every later idea the matrix `K`. |
| 2 | **3 — exact contact angles from single harmonics** | Directly answers the paper's own named open problem (Sec. 6 limitation 1) and can falsify the split-edge assumption in the same run. |
| 3 | **8 — the complex torus deployment variety** | The most surprising claim on the list and the only exact finite answer to Sec. 6 limitation 2; also the one most likely to be old, so it must be prior-art-checked first. |
| 4 | **2 — hole area harmonic, closed-form closure angle** | Clean, cheap, and turns Sec. 5.2's discovery problem into one equation; weak alone, strong merged with 3. |
| 5 | **5 — deployability defect beats max-cut for `σ`** | A constructive counterexample to a stated design principle, with a real baseline to beat; the existence half is nearly certain. |
| 6 | **4 — range maximization with a certified bound** | Largest practical payoff and the only certificate in the field, but the relaxation may be vacuous and the code is the heaviest. |
| 7 | **7 — sparsest design-space basis** | Fixes a sentence the paper gets wrong about its own output; modest ceiling. |
| 8 | **6 — exact planar max-cut for orientations** | The complexity observation is sharp and the measurement is easy, but the likely outcome is that the heuristic is already optimal. |
| 9 | **9 — certified inverse design by low-rank relaxation** | High value if tight, but `dim X` may be 0 exactly where inverse design lives. |
| 10 | **10 — refuting "small balanced holes"** | A correlation study; belongs inside another idea, not on its own. |

---

## Overlap with `ideas/persona_geometer.md` and `ideas/persona_rigidity.md`

Read **after** my ten were drafted; nothing below was copied, and nothing above was changed in
response to them. Both files are excellent and both beat me to my two highest-ranked ideas.

### Fully duplicated — I withdraw the novelty claim

| My idea | Theirs | Verdict |
|---|---|---|
| **1** — `J(θ) = cos(θ/2)I + sin(θ/2)K`, conformality theorem | **geometer G2** (`J = cI + sW`, `W = QP₀⁻¹`) and **rigidity Idea 4** (`J = cI + sM`, `M = P₁'P₀⁻¹`) | Independently derived three times, identically, including the "Eq. (13) is exact, not a linearization" punchline. This is now a *confirmed* result, not a candidate idea. Priority is theirs; my only additions are the exact `det J(θ) = c² + cs·tr K + s² det K` and the framing of `{K(X) : X ∈ X}` as an **affine subspace**, which turns the paper's second empirical claim ("Eq. (13) always succeeds when the kernel is nontrivial") into a *decidable* rank test. Neither file states that. |
| **3** — orientation predicates are single harmonics, exact contact angles | **geometer G3** (identical, same root formula) and **rigidity Idea 3** (same, in `t = tan(θ/2)` form, with `Y_θ = cC + sS` already verified numerically on 21 graphs) | Fully duplicated, down to the split-edge-assumption test. Theirs is further along: rigidity reports V1, a measured `6×10⁻¹⁵` residual. Drop mine. |
| **4** — range-maximizing embedding | **geometer G10** and **rigidity Idea 9** | Same objective, same baseline. My only distinct content is the **Shor/Burer–Monteiro certified upper bound**; both of theirs claim local optima only and say so explicitly. If this line is pursued, the certificate is the piece worth keeping from my version. |
| **6** — exact planar max-cut for `σ` | **geometer G6** | Theirs is strictly stronger: it proves `deg_split(v) = deg(v) mod 2`, identifies `E_split` as a **T-join**, and derives a *lower* bound `dim(shape space) ≥ |T| + 2 dim Z`. Withdraw mine entirely. |

### Partial overlap — a real difference remains

- **My 8 (complex torus variety)** vs **geometer G9** and **rigidity Idea 6.** G9 has the same
  variety in real form, `Σ_γ (R(ω_{f₁}) − R(ω_{f₂}))x_u = 0`, and computes DOF by linearizing;
  rigidity Idea 6 has the infinitesimal version `Θ = D(ker A)`. What neither states is that under
  `z_f = e^{iψ_f}` those equations are **complex-linear in `z`**, so the configuration space is
  `ker_C W ∩ (S¹)^{|F|}` — a linear section of a torus. Consequences only my form gives: the global
  rotation is scalar multiplication `z ↦ λz`; the uniform path is the circle in `span_C{1, σ}`;
  "design an extra finite mode" becomes a **complex rank-drop condition on a `b₁ × |F|` matrix
  linear in `X`**; and the finite (not infinitesimal) question reduces to `|F|` real quadratics in
  `2 dim_C ker W` unknowns. My own self-attack notes this is close to the classical complex
  loop-equation method of mechanism design, which is a gap in all three of our prior-art tables.
- **My 7 (sparsest design-space basis)** vs **geometer G7** (`L = A_S Δ_h`, `null(L) = Δ_h⁻¹(B)`).
  G7 characterizes the null space as a preimage and bounds its dimension; it does not ask for a
  **basis**, and does not make the point that the paper's sentence "each `φ_i` identifies a set of
  vertices that can be translated together" is false for the SVD basis it actually computes. The
  matroid-greedy sparsest-basis claim is mine and is compatible with G7 as its algebraic engine.

### Genuinely different — not present in either file

1. **My 5 — the deployability defect `D(σ)` as the objective for face orientation.** Neither file
   proposes optimizing `σ` against the *input geometry*. G6's risk paragraph notices that
   `|E_split|` may be the wrong objective but proposes nothing in its place; my claim is sharper and
   constructive: split cuts merge holes and therefore let residuals **cancel**, so more split cuts
   can strictly reduce the defect, and there is a one-minute enumeration that decides it.
2. **My 2 — total hole area is a single harmonic; the second gap-free angle is
   `2 arctan(r/q) = 2 arctan(tr K/(1 − det K))` in closed form.** Neither file mentions hole area,
   Eq. (14), or fully-closed configurations at all (grep: no hits for "hole area", "fully closed",
   "Eq. (14)" in either). This is the whole of 2026 Sec. 5.2 and 2025 Observation 3.1, unclaimed.
3. **My 9 — certified global inverse design by low-rank relaxation.** No convex relaxation, duality
   bound, or global-optimality certificate appears anywhere in either file; rigidity's Idea 10 uses
   the word "certificate" for a prestress test, which is a different object.
4. **My 10 — measuring the "many balanced, small holes" principle.** Both files quote the principle
   as unsupported; neither proposes an experiment that would confirm or refute it.

### Consolidated recommendation

Ideas 1, 3, 4 and 6 on my list should be treated as **already claimed** by the other two personas and
merged into their versions, with three transplants: the affine-achievable-set decision test into G2,
the certified upper bound into G10/rigidity-9, and nothing into G6 (it dominates my Idea 6). My
remaining distinct contributions, ranked, are **5** (defect vs max-cut), **2** (closure angle in
closed form), **8** in its complex-linear form only, **7**, **9**, **10**.

---

## What I could not check

- **Nothing in this file was executed.** No code was built or run; every number quoted is a
  prediction. Fact (i) I re-derived by hand from `code/src/core/kinematics.hpp`'s stated convention;
  fact (ii) I took on trust from ideator-rigidity.
- **"Max-cut on planar graphs is in P" (Orlova–Dorfman 1972 / Hadlock 1975)** — I did not verify this
  reference in this session, and neither did `notes/field_kirigami.md` or `notes/field_tutte.md`
  (no row covers max-cut). `ideas/persona_geometer.md` cites Hadlock 1975 for the same fact, which
  is corroboration by another agent, not verification. It must be checked before use.
- **Classical complex loop equations in mechanism design.** My Idea 8 is very likely adjacent to the
  standard complex-vector loop-closure method for planar linkages. None of the three field files has
  a row for it (`field_rigidity` R1–R26 covers rigidity theory and mobility, not kinematic synthesis
  by complex numbers). This is an unchecked prior-art hole in my strongest distinct idea.
- **Coleman–Pothen sparse-null-basis matroid condition** (Idea 7) is cited from memory; not verified.
- **Whether `dim X > 0` at all** for the random generators — every optimization idea here (4, 7, 9,
  10) is vacuous when the shape space is trivial, and no one has yet measured how often that happens.
  This is the single cheapest measurement that should precede committing to any of them.
- The 2025 supplement (Secs. A–F) is unavailable, so the baseline weights `ω₁..ω₄` for any comparison
  against 2025 must be reconstructed, and any "we beat their optimizer" claim inherits that
  uncertainty.
