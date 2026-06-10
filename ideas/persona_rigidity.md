# Ideas — persona: RIGIDITY (mechanisms, configuration spaces, body-and-hinge frameworks)

Ideator subagent. Inputs actually read (2026-09-03):
- `specs/ideator.md`, `specs/common_preamble.md`
- `STATE.md` facts F1–F13 (F13 = Maxwell count ⇒ uniform deployability is a special-position mechanism)
- `notes/paper_2026.md` §§2, 3 (Eq. 2–13), 7.3 (rank claims), 8 (collision), 12 (assumptions), 13 (limitations), 15 (open questions)
- `notes/paper_2025.md` §2.10–2.11 (θ_max), §13 (limitations)
- `notes/field_rigidity.md` — full (R1–R26, §2 "known vs open")
- `notes/field_kirigami.md` — table R1–R24, §2 seed verdicts, orchestrator IsoGami addendum
- `notes/field_tutte.md` — **partially written when I read it** (it was a 30-line skeleton, then grew to §0–§1
  while I worked). I used its §0.1–0.5 (the `L = R·D` factorisation and the periodic rank-deficiency
  prediction) and credit Scout-c for it below. I did **not** see its §2–§3.
- `papers/related/isogami.txt` + the orchestrator's verified IsoGami addendum in `notes/field_kirigami.md`.

## 0. What I verified numerically before writing (this matters — read this first)

I compiled a standalone C++ program against the existing core
(`code/src/core/{mesh,cut,holes,kinematics,generators,orientation,tutte_auxetic}.cpp`;
the repo's own `CMakeLists.txt` does not configure yet because `collision.cpp`/`kiri_sweep`
are still missing, so I compiled the sources directly with clang++ `-arch arm64`).
Source: `scratchpad/rig_check.cpp`. 21 graphs: squares 3×3 and 5×5, triangles, hexagons,
kagome, periodic squares, 8 random Delaunay (~50 faces), 8 random Voronoi (25 faces).
σ from `brute_force_orientation` (|F| ≤ 18) or `assign_orientation_relaxation`.

| Claim | Result |
|---|---|
| **(V1)** `Y_θ = cos(θ/2)·C + sin(θ/2)·S` with `C, S` independent of θ | max error **6×10⁻¹⁵** at scale ~6 on **21/21** graphs, at θ ∈ {0.05, 0.7, 1.9, 2.5, 3.0} with `C,S` fitted from θ ∈ {0.37, 1.13} |
| **(V2)** `m = dim ker A − c(Γ)` equals `3\|F\| − rank(rigidity matrix) − 3c(Γ)` | agreed on **21/21** graphs, at both `X_ini` and the solved `X_0` |
| **(V3)** `A σ = 0 ⟺ Eq. (2) residuals vanish` | on `X_ini` (not deployable): `‖Aσ‖ ∈ [3.2, 11]`, max hole residual ∈ [0.65, 3.1]. On solved `X_0`: both `< 2×10⁻¹⁴`. **21/21** |
| **(V4)** infinitesimal mobility along the deployment path | **drops discontinuously at θ = 0 on 20/21 graphs** and is then constant for θ > 0. Examples: squares 5×5 `10 → 1`; kagome `18 → 7`; periodic squares `4 → 1`; a random Delaunay `26 → 19`. Only the triangle tiling kept `m = 9` throughout. |
| `rank(L) = H` (paper's Sec. 4.4 claim 2) | held in **21/21** with a real (non-periodic) boundary. Scout-c's prediction that it fails in the *periodic* case is **not tested by my run** — I used `BoundaryMode::Fixed` throughout. |

Two consequences that shape every idea below and that neither Segall paper contains:

1. **Deployment is a rational curve of degree 2.** With `t = tan(θ/2)`, `Y(t) ∝ C + t·S`, so every
   polynomial event (collinearity, contact, area, lattice conformality) is an explicit low-degree
   polynomial in `t`. "Explicit deployment simulation via forward kinematics" (2026 Sec. 6,
   limitation 1) is not needed.
2. **A uniformly deployable Segall-2026 pattern is essentially never a 1-DOF mechanism.** Measured
   mobilities on solved embeddings: 1, 2, 3 (Voronoi, many split cuts) up to 19–26 (Delaunay).
   2026 Sec. 6 limitation 2 asks for exactly this number and offers "sensitivity analysis"; the
   answer is a rank of an `2H × |F|` matrix, not a sensitivity study.

Everything below marked **[verified here]** is one of the rows above. Everything else is a
conjecture with a kill test.

---

## The common algebraic object (used by ideas 1–8)

Body-and-pin model (F13, `notes/field_rigidity.md` §0): faces `F` are rigid bodies (3 DOF each),
each hinge edge `e ∈ E_hinge` is a pin joint at its **source** vertex (F1), split edges impose
nothing. Hinge graph `Γ = (F, E_hinge)`. **Γ is bipartite by σ**, because a hinge edge by
definition joins faces of opposite orientation (F1) — so the alternating vector σ ∈ {±1}^F is
well defined, which is *why* a common opening angle can exist at all.

Infinitesimal motion of face `f`: angular velocity `ω_f`, translation `v_f`. Pin at `p_e`:

```
        (v_f − v_g) + (ω_f − ω_g) J p_e = 0 ,          J = [[0,−1],[1,0]].          (R1)
```

Given ω, a compatible v exists iff the edge vector `e ↦ (ω_f − ω_g) J p_e` is a gradient, i.e. iff
its circulation vanishes on a cycle basis `{z}` of Γ. Dropping the invertible `J`:

```
   A(ω)_z  =  Σ_{e ∈ z} z_e (ω_{h(e)} − ω_{t(e)}) p_e  =  Σ_i ω_{f_i} (p_{i−1} − p_i)  ∈ R².   (R2)
```

`A` is a **2·b₁(Γ) × |F|** matrix, linear in ω and linear in the hinge-point positions `p`. With
`W = ker A` and `c = c(Γ)`:

```
   dim ker(rigidity matrix) = dim W + 2c ,      m = dim W − c ,
   s = 2|E_hinge| − 3|F| + dim W + 2c ,    hence   m − s = 3|F| − 2|E_hinge| − 3c.       (R3)
```

The last identity is exactly Maxwell/Calladine (F13, R1), so (R3) is Maxwell with the unknown `s`
supplied. Three checks I ran by hand: one body ⇒ m = 0; two bodies one pin ⇒ m = 1; three bodies
pinned pairwise ⇒ m = 0 generically and **m = 1 exactly when the three pins are collinear**, which
is the Aronhold–Kennedy theorem falling out of (R2). Not in R1–R26 in this form.

And the punchline connecting the two papers' worlds:

```
   B(ω, X)_z := Σ_i ω_{f_i} (p_{i−1}(X) − p_i(X))   is bilinear;
   L  =  B(σ, ·)   (the paper's design matrix, Eq. 4)  and   A  =  B(·, X)   (the mobility matrix).
   Uniform deployability (Eq. 2)  ⟺  σ ∈ ker A  ⟺  X ∈ ker L.                              (R4)
```

**[verified here, V3]**. The paper's entire Sec. 4.4 is one slot of a bilinear form whose other
slot is the mechanism theory it says (Sec. 6) it does not have.

---

# The ten ideas

---

## Idea 1 — The mobility operator of a hinged kirigami
**Type: characterization (+ theorem).**

**Claim.** For a kirigami structure `M'` with hinge graph `Γ` having `c` components and cycle rank
`b₁(Γ) = |E_hinge| − |F| + c`, the mobility is exactly `m = dim ker A − c`, where `A` is the
`2b₁(Γ) × |F|` matrix of (R2); `A` depends on the embedding only through the hinge points `p_e`, is
affine-covariant, and `σ ∈ ker A` iff the structure is uniformly deployable. Equivalently the
`2|E_hinge| × 3|F|` rigidity matrix is replaced, with no loss, by an operator with **one scalar
unknown per face**.

**Gap targeted.** 2026 Sec. 6, limitation 2 verbatim: "patterns that can be uniformly deployed might
also be deployed with a different set of hinge angles. Analyzing the degrees of freedom of the
deployment space, which could be achieved by sensitivity analysis, is an interesting direction to
explore." Also 2025 open question 19 in `notes/paper_2025.md` ("the paper never performs a
Maxwell/Calladine count"), and F13.

**Novelty vs field.** R26 (Chen–Choi–Mahadevan 2020) computes `DoF = 8L² − rank(A)` numerically for
an `L×L` quad array with a *specific* geometry — the same bookkeeping, but no reduction, no
arbitrary graph, no link to a deployability condition. R4/R6 (Tay–Whiteley, Katoh–Tanigawa) give
generic *rigidity* criteria for body-and-hinge multigraphs; `notes/field_rigidity.md` R7 records
that they return "rigid" for every working pattern, so they cannot produce this number. IsoGami
(F12) counts zero eigenvalues of `JᵀJ` at one configuration for isohedral tilings only, with no
reduction and no combinatorial statement.

**Why it might be true / sketch.** (R1)–(R3) above. The reduction is the observation that (R1) is a
*gradient* condition on Γ: `v` exists iff the edge 1-form `(Dω)_e J p_e` has zero circulation, which
is `b₁(Γ)` vector equations. The `+2c` are the per-component translations, the `−3c` in Maxwell are
translations plus the constant ω. Recovering Maxwell exactly (R3) is the dimensional check; the
Aronhold–Kennedy collinear-pin case is the limiting-case check.

**Kill experiment.** `code/src/core` + a new `mobility.cpp`. Generators: `delaunay_of_random_points`,
`voronoi_of_random_points`, `quad_dominant_random`, plus the 7 tilings, 500 graphs, 100–5000 faces.
For each: σ from `assign_orientation_relaxation`, `make_cut`, build `A` from a spanning-tree cycle
basis and the full body-pin rigidity matrix, compare `dim ker A − c` with `3|F| − rank(R) − 3c` by
sparse QR / SVD with a fixed relative tolerance. **If they disagree on any graph, dead.**
Already run on 21 graphs at both `X_ini` and solved `X_0`: **21/21 agree** [verified here, V2].
Under 30 min for 500 graphs because `A` has `|F|` columns.

**If it survives, the demo.** A table of measured `m` for all figure patterns of both papers,
against the generic (pebble-game / Tay–Whiteley) prediction, which is "rigid" everywhere — the
special-position excess `m_actual − m_generic` plotted over ≥20 graphs. Hero: the Fig. 11 pattern
of 2026, with its 6-dimensional *design* space and its (measured) mobility side by side, showing
they are different numbers with different meanings. Fabricable export: the pattern annotated with
its non-uniform modes.

**Risk.** That some rigidity paper already writes the body-and-pin rigidity matrix in this reduced
`ω`-only form. It is close to classical instant-centre kinematics, and Aronhold–Kennedy is 19th
century. Mitigation: the *combination* with (R4) — that the same bilinear form gives both Segall's
design matrix and the mobility matrix — is what has to be new, and that requires the 2026 paper.

**Effort.** Derivation S (done). Code M.

---

## Idea 2 — The flat state is a branch point: mobility drops the instant you deploy
**Type: theorem.**

**Claim.** For a uniformly deployable kirigami, `m(θ) := dim ker A(Y_θ) − c` satisfies
`m(0) > m(θ) = const` for all `θ` in the open deployment interval, except on a measure-zero set of
embeddings. In particular the flat configuration is a **singular point of the configuration
variety**: infinitesimal flexes counted at `θ = 0` overcount the actual mechanisms by
`m(0) − m(0⁺)`, and the uniform mode is one branch of several meeting there.

**Gap targeted.** 2026 Sec. 4.5 / Eq. (9): the collision optimisation is solved **at `θ = 0`**
("Optimized at `θ = 0`, i.e. on the *undeployed* configuration", `notes/paper_2026.md` §8.4) — that
is, at precisely the configuration where the linearisation is not representative. Also 2026 Sec. 6
limitation 2. `notes/field_rigidity.md` §2(b) lists "whether uniform-θ deployment is a smooth point
or a branch point" as open with no theory answering it.

**Novelty vs field.** R19 (Kumar–Pellegrino 2000) detects bifurcation points numerically along a
traversed path for pin-jointed bars, but must be *told* where to look and studies particular
structures. R9 (Connelly–Servatius) is the warning that cusps exist, not a statement about kirigami.
R25 (PNAS rigid-origami mobility, unverified second-hand in the field notes) classifies origami as
regular/bifurcated/shaky but is 3D edge-hinged. None of them says the *flat* state of a corner-pinned
tiling is systematically the singular one.

**Why it might be true / sketch.** At `θ = 0` every pair of duplicates of a hinge edge coincides and
all hinge points sit at the original vertices `X`; entire families of pins become collinear (whole
edges of `M` are lines of pins), which by the Aronhold–Kennedy case of (R2) manufactures extra
kernel. For `θ > 0` the faces rotate by `∓θ/2` and those collinearities break. Concretely,
`A(θ) = cos(θ/2) A₀ + sin(θ/2) A₁` (idea 5), so `rank A(θ)` is the rank of a pencil: it is constant
off a finite set of `t = tan(θ/2)` and takes its *minimum* exactly on that set. `t = 0` is in the set.

**Kill experiment.** Same driver. For 500 solved embeddings `X_0`, compute `m(θ)` at
`θ ∈ {0, 10⁻⁶, 10⁻³, 0.2, 0.6, 1.0, 1.6, 2.2, 2.8}`. **If `m(0⁺) = m(0)` on more than a negligible
fraction, or if `m(θ)` is not constant on `(0, θ_max)`, the sharp form is dead.** Already run on 21:
drop on **20/21** (only the triangle tiling kept `m = 9`), constant thereafter [verified here, V4].
The triangle tiling is a genuine exception and must be explained, not hidden — that is the first
piece of real work.

**If it survives, the demo.** `m(θ)` staircase plots for ≥20 graphs; a hero example where the two
branches at `θ = 0` are both followed by forward kinematics and shown to be geometrically distinct
patterns from the same flat sheet. Fabricable export: one cut sheet, two deployed states.

**Risk.** That the drop is an artefact of numerical rank at `θ = 10⁻⁶` rather than a genuine
stratification. Mitigated by using the pencil (idea 5), where the drop is an exact statement about
`rank(cA₀ + sA₁)` and needs no small-θ probing.

**Effort.** Derivation M. Code S (the machinery is idea 1's).

---

## Idea 3 — The deployment path is a conic; every contact event has a closed form
**Type: characterization.**

**Claim.** For any `X` and σ, forward kinematics satisfies exactly
`Y_θ = cos(θ/2)·C(X) + sin(θ/2)·S(X)` with `C, S` **linear** in `X`; hence with `t = tan(θ/2)` the
deployed configuration is `(C + tS)/√(1+t²)`, every pairwise cross product
`(y_i − y_j) × (y_k − y_l)` is a **homogeneous quadratic in (cos θ/2, sin θ/2)** whose coefficients
are quadratic forms in `X`, and the first collinearity/contact angle of any edge pair is the
smallest positive root of an explicit quadratic in `t`. Therefore
`θ_max = 2 arctan( min over candidate pairs of the smallest positive root )` in closed form.

**Gap targeted.** 2026 Sec. 6, limitation 1, verbatim: "self-intersections may occur at small
opening angles ... At present, such issues can only be detected through explicit deployment
simulation via forward kinematics and require additional post-processing to resolve. Developing
geometric characteristics for favorable deployment behavior directly from the embedding remains an
open problem." Also 2026 Sec. 4.5's empirical restriction of collisions to `E_split`
("we observe in all experimental cases"), which becomes checkable rather than assumed. 2025
Sec. 4.2's `θ_max = min_v (2π − α_i − α_j)` is only the *local back-side* bound at a hinge and
ignores non-adjacent face pairs entirely.

**Novelty vs field.** Kirigami R18 (Liu et al. 2024) gives a closed-form inequality
`length(C₁p₁) + length(C₃p₂) < 2a cos θ` but only for enumerated adjacent/opposite quad pairs, and
falls back to a manual UI otherwise. IsoGami (F12) handles contact with an IPC log-barrier, i.e.
numerically, and reports contact only "blocking further expansion". R22 (Tachi) and R19 are
numerical path-followers. Nobody states that the deployment path of a corner-pinned rigid tiling is
a degree-2 rational curve.

**Why it might be true / sketch.** Faces rotate by `R(∓θ/2)`; across a hinge at `x_p`,
`t_f − t_g = (R(θ/2) − R(−θ/2)) x_p = 2 sin(θ/2) J x_p`, so every face translation is
`sin(θ/2) · w_f` with `w_f` a signed sum of `2J x_p` along a BFS path — **linear in X**. Then
`y = R(∓θ/2)x + sin(θ/2)w_f = cos(θ/2)·x + sin(θ/2)·(∓Jx + w_f)`. Path-independence of `w_f` around
cycles of Γ is precisely Eq. (2). Products of two such expressions are homogeneous quadratics in
`(cos, sin)`; dividing by `cos²` gives a quadratic in `t`.

**Kill experiment.** Fit `C, S` from `deploy()` at two angles and predict `Y_θ` at five others on
500 graphs. **If the max relative error exceeds 10⁻⁹ on any graph, dead.** Already run on 21 graphs:
max error `6×10⁻¹⁵` at scale ~6 [verified here, V1]. Second half: implement the quadratic-root
`θ_max` and compare against the existing bisection-based collision detector (once
`code/src/core/collision.cpp` lands) on 500 graphs; **if the closed-form and bisection `θ_max`
differ by more than 10⁻⁶ rad on any graph, dead.**

**If it survives, the demo.** `θ_max` computed in closed form vs. bisection on ≥20 graphs — same
answer, orders of magnitude faster and differentiable. A plot of "which pair of edges is the binding
constraint", testing 2026's `E_split`-only assumption directly. Hero: a pattern where the binding
pair is *not* a split-edge duplicate pair, refuting Sec. 4.5's empirical assumption if it exists.

**Risk.** The half-angle convention. The code uses `R(−σ(f)θ/2)`, the paper's prose says "each face
undergoes a rigid rotation by θ" (`notes/paper_2026.md` §8.3 flags this as loose). The result holds
in either convention with `θ` or `θ/2`; only the constant differs. Second risk: the closed form is
per *pair*, and the number of candidate pairs is quadratic unless a broad phase is added.

**Effort.** Derivation S (done). Code M.

---

## Idea 4 — Conformal at θ = 0 implies conformal for every θ (a two-line theorem)
**Type: theorem.**

**Claim.** For a periodic pattern with metatile periodicity matrix `P_θ`, the deployment Jacobian of
2026 Eq. (12) is exactly `J(θ) = cos(θ/2)·I + sin(θ/2)·M` for a fixed `M = P₁'P₀⁻¹`. Consequently
`J(θ)` is a similarity for all θ **iff** `M` is a similarity **iff** the two linear conditions of
Eq. (13) at `θ = 0` hold. The paper's empirical claim is a theorem, and Eq. (13) is not a heuristic
linearisation but the exact condition.

**Gap targeted.** 2026 Sec. 5.1, quoted in `notes/paper_2026.md` §13.2: "Empirically, we observe...
once the derivative of the conformal distortion vanishes at `θ = 0`, the deployment remains conformal
for all opening angles." No proof, no experiment count. Listed as implicit assumption 19 in §12.

**Novelty vs field.** R17 (Borcea–Streinu, geometric auxetics) characterises which one-parameter
paths are auxetic via a positive-semidefinite Gram velocity, *assuming a path exists*, and does not
produce this normal form. R11 (Guest–Hutchinson) is a determinacy dichotomy. R24 (Schenk–Guest)
assumes the 1-DOF kinematics. Kirigami R1 (Grima–Evans) asserts `ν = −1` for rotating squares as a
modelling statement. Nobody writes `J(θ) = cI + sM`.

**Why it might be true / sketch.** By idea 3 every vertex is `cos(θ/2)·(·) + sin(θ/2)·(·)`, so each
periodicity vector is `p^x_θ = c a₁ + s b₁`, `p^y_θ = c a₂ + s b₂`, hence `P_θ = cA + sB` with
`A = P₀`. Then `J = P_θP₀⁻¹ = cI + sM`, `M = BA⁻¹`. `cI + sM` is a similarity for one `(c,s)` with
`s ≠ 0` iff `M` is, and `∂J/∂θ|₀ = ½M`, whose `(11−22)` and `(12+21)` parts are what Eq. (13) sets
to zero. Both directions are immediate. Poisson ratio `−1` at every θ follows.

**Kill experiment.** On all periodic generators (`periodic_squares/triangles/hexagons/kagome`) plus
≥100 random periodic patterns: solve Eq. (13), then measure the conformal distortion
`|J₁₁−J₂₂| + |J₁₂+J₂₁|` at 50 angles across the range. **If it is nonzero beyond 10⁻¹⁰ at any θ on
any pattern, the theorem is false and dead.** Separately: on patterns where Eq. (13) is *not*
imposed, check that `J(θ) − cos(θ/2)I` has rank-1-in-θ structure; **if `J(θ)` is not affine in
(cos θ/2, sin θ/2), dead** (this is implied by V1, already verified).

**If it survives, the demo.** A closed-form Poisson-ratio curve `ν(θ)` for every periodic pattern in
both papers, plus the exact statement of *which* metatiles admit a conformal solution (an algebraic
condition on `M`). Hero: a pattern where Eq. (13) is infeasible, which the paper's empirical claim
("for all tilings with a non-trivial kernel, solving Eq. (13) yields an embedding with conformal
deployment") predicts cannot exist.

**Risk.** It is a small theorem — a referee may call it a remark. Its value is as a lemma inside a
larger paper (ideas 1–3) and as a correction of a stated empirical claim. Second risk: the metatile
`P_θ` may be defined with a re-chosen metatile per θ, in which case `A`, `B` are not fixed.

**Effort.** Derivation S (done). Code S.

---

## Idea 5 — Bifurcation angles are generalised eigenvalues of a matrix pencil
**Type: theorem.**

**Claim.** Along the uniform deployment path the mobility operator is a **linear matrix pencil**
`A(θ) = cos(θ/2)·A₀ + sin(θ/2)·A₁` with `A₀, A₁` fixed and linear in `X`. Hence `rank A(θ)` is
constant except at the finitely many `θ` with `tan(θ/2)` a generalised eigenvalue of `(A₀, A₁)` (in
the Kronecker sense, including `∞`), and every configuration at which the kirigami gains an
infinitesimal mechanism — every candidate bifurcation of the deployment — is a root of an explicit
determinantal polynomial in `t = tan(θ/2)`. No path-following and no probing is needed.

**Gap targeted.** 2026 Sec. 6 limitation 2 again, and `notes/field_rigidity.md` §2(b): "Whether
uniform-θ deployment is a *smooth point* or a *branch point* of the configuration space. R19's
algorithm would settle it computationally for a given pattern; no theory answers it." And "Whether
the mobility is constant along the deployment path. Nothing in the literature addresses this."

**Novelty vs field.** R19 (Kumar–Pellegrino) locates bifurcations numerically by traversing a path.
R20 (Müller) computes higher-order local mobility by nested Lie brackets with no a-priori order
bound. R12/R13 (Fowler–Guest, Schulze–Guest–Fowler) detect hidden modes only via point-group
symmetry and are useless for the arbitrary asymmetric graphs the 2026 paper advertises. A pencil
gives all of them at once, globally, with no symmetry hypothesis and no traversal.

**Why it might be true / sketch.** By idea 3, `p_e(θ) = cos(θ/2)p_e⁰ + sin(θ/2)p_e¹`. `A` in (R2) is
linear in `p`, so `A(θ) = cA₀ + sA₁`. For a rectangular pencil the rank is the *normal rank* except
where the Kronecker structure degenerates; the exceptional `t` are the finite eigenvalues of the
regular part. `t = 0` (the flat state) is generically among them — which is exactly idea 2. Extra
flexes are the eigenvectors.

**Kill experiment.** 500 graphs. Compute the exceptional `t` from a generalised SVD / staircase
algorithm on `(A₀, A₁)`, then independently sample `rank A(θ)` on a dense grid of 2000 angles.
**If any sampled rank drop occurs at a `θ` that is not a computed eigenvalue, or the rank is not
constant between consecutive eigenvalues, dead.** Cheap: `A` has `|F|` columns.

**If it survives, the demo.** A "bifurcation spectrum" per pattern: the deployment interval marked
with its exceptional angles and, at each, the branch directions from the eigenvectors, followed by
forward kinematics into visibly different states. Comparison against R19-style path following on
≥20 graphs: same bifurcations, found without traversal. Hero: a pattern with an interior bifurcation
where one flat sheet reaches two distinct deployed shapes.

**Risk.** That in practice `A₀` and `A₁` share a large common kernel so the pencil is *singular*
(not regular) for every kirigami, in which case the Kronecker minimal indices dominate and the
"finitely many bad angles" statement degenerates. This is a real possibility and it is the first
thing the kill experiment will show. Also: rank drop is necessary but not sufficient for a genuine
branch — see idea 10.

**Effort.** Derivation M. Code M (needs a numerically careful pencil staircase; Eigen has GSVD only
via QZ on square pencils, so this may need a Kronecker-form implementation).

---

## Idea 6 — The space of admissible non-uniform opening-angle fields
**Type: characterization.**

**Claim.** The admissible *infinitesimal* opening-rate fields on the hinge set are exactly
`Θ = D(ker A) ⊂ R^{E_hinge}` where `D` is the incidence map of Γ, i.e. `θ̇_e = ω_{h(e)} − ω_{t(e)}`.
`dim Θ = dim ker A − c = m`. The uniform field `θ̇ ≡ const` lies in `Θ` iff `σ ∈ ker A` iff Eq. (2),
so the paper's condition is the statement "the constant field is admissible", and `Θ` is the *whole*
answer to "which other hinge-angle assignments deploy this pattern".

**Gap targeted.** 2026 Sec. 6 limitation 2, verbatim and in full. Also 2026 Sec. 2.11's unsupported
aside "Note that there can be multiple sets of hinge angles that allow for a rigid deployment",
which offers no example and no proof (`notes/paper_2026.md` §2.11), and 2025's open question list.

**Novelty vs field.** `notes/field_kirigami.md` seed 5 verdict: "Nobody I found treats θ as a field
with a compatibility condition of its own. The obvious mathematical question — what is the space of
admissible θ fields on a given cut structure — appears unasked." R25 (Dorn–Lang–Pellegrino) has many
DOF via sub-folds, a different mechanism. IsoGami's future work names "graded patterns" as open
(F12).

**Why it might be true / sketch.** Immediate from (R1): a flex is `(ω, v)`, the opening rate at a
hinge is the relative angular velocity `ω_{h} − ω_{t}`, and `ω` ranges over `ker A`. `ker D` is the
constants per component, so `dim Θ = dim ker A − c = m`. The measured `m` (idea 1) is therefore the
*dimension of the space of graded deployment fields* — for a random Delaunay pattern in my run that
is 19 to 26 dimensions, not 1. The finite question is whether a given `θ̇ ∈ Θ` integrates; by idea 3
the uniform one does, exactly and in closed form.

**Kill experiment.** 500 graphs: compute `Θ` and verify `dim Θ = m` and `1 ∈ Θ ⟺ hole residuals
vanish`. **If `dim Θ ≠ m` or the equivalence fails on any graph, dead.** Then integrability: pick a
random `θ̇ ∈ Θ` orthogonal to the uniform mode, integrate the flex numerically with re-projection
onto the pin constraints, and measure whether the pin residual can be driven below 10⁻¹⁰ at finite
step size. **If every non-uniform direction is infinitesimal-only on every graph, the "graded
deployment" half is dead** (the dimension count survives).

**If it survives, the demo.** A designed graded pattern: prescribe a target `θ̇` field (say, a
gradient across the sheet), project onto `Θ`, deploy, and show a spatially varying opening. ≥20
graphs comparing prescribed vs achieved fields. Fabricable export: a sheet whose holes open
progressively.

**Risk.** That `Θ` is dominated by boundary-face free rotations (idea 7), so the "graded field" is
just floppy edges rather than a designable interior field. The fix and the real content is to report
`Θ` restricted to the 2-core.

**Effort.** Derivation S. Code M.

---

## Idea 7 — Mobility splits into a combinatorial boundary part and a geometric core part
**Type: theorem.**

**Claim.** `ker A = R^{F ∖ core₂(Γ)} ⊕ ker A|_{core₂(Γ)}`, hence
`m = |F ∖ core₂(Γ)| + (dim ker A_core − c)`, where `core₂(Γ)` is the 2-core of the hinge graph
(iteratively delete degree-≤1 faces). The first term is purely combinatorial and is exactly the
"floppy boundary" degrees of freedom; the second is the only geometrically interesting number.
Corollary: the folklore that rotating squares is a 1-DOF mechanism is **false for every finite
patch** and true only after quotienting the non-2-core, or in the periodic setting.

**Gap targeted.** 2026 assumption 7 / Sec. 4.2 principle (1): full connectivity of `M'` is required
by the definition of a uniformly deployable embedded graph but "never enforced as a hard constraint
anywhere in the optimization" (`notes/paper_2026.md` §2.12). Connectivity is 1-connectivity; this
idea says the right invariant is the 2-core. `notes/field_rigidity.md` §2(b): "Boundary effects are
a near-certain source of extra DOF for finite patches ... I could not find this analysed for hinged
kirigami specifically", and §2(b) again: "The claim that rotating squares have exactly 1 DOF: I
found it asserted repeatedly and proved nowhere" (R15).

**Novelty vs field.** R26 (Chen–Choi–Mahadevan) sees DOF percolation transitions in cut density but
works at fixed geometry and gives no decomposition. R7 (pebble game) computes rigid clusters
generically and, per the field notes, "will confidently give the wrong number" here. R15
(Grima–Evans) asserts 1-DOF without proof. No 2-core statement exists for this model.

**Why it might be true / sketch.** `A(ω)_z` in (R2) is supported on faces lying on the cycle `z`.
Faces outside the 2-core lie on no cycle, so their `ω` is unconstrained: `A` has zero columns there.
That gives the direct sum immediately, and the count follows. The rotating-squares corollary is then
arithmetic: in an `n×n` patch the corner and edge squares that carry fewer than two pins fall out of
the 2-core.

**Kill experiment.** 500 graphs: compute the 2-core of Γ and check
`dim ker A = |F ∖ core₂| + dim ker A_core` exactly (integer equality of ranks). **If it fails on any
graph, dead.** Plus the rotating-squares family `n = 2..40`: check `m(n)` against the closed formula
the decomposition predicts, and separately against `3|F| − rank(R) − 3`.

**If it survives, the demo.** `m` vs patch size for the rotating-squares family with the
combinatorial and geometric parts separated, showing the folklore claim is only the core term.
≥20 graphs with the split reported. Hero: a pattern whose entire apparent mobility is boundary flop,
next to one whose mobility is genuinely interior.

**Risk.** It may be too easy — a referee could call the 2-core observation obvious once (R2) is
written down. Its defence is the rotating-squares corollary, which contradicts a claim repeated
throughout the auxetics literature (R15) and was never proved either way.

**Effort.** Derivation S. Code S.

---

## Idea 8 — Generic embeddings in the shape space have minimal mobility
**Type: theorem (+ design algorithm).**

**Claim.** On the shape space `X = ker[L; B]` (2026 Eq. 5), the function `X ↦ dim ker A(X)` is upper
semicontinuous and its minimum `m*` is attained off a proper algebraic subvariety `V ⊂ X` defined by
the vanishing of the `m*`-sized minors of `A(X)`. Hence: **almost every uniformly deployable
embedding of a given cut structure has the same mobility `m*(σ, G)`, an invariant of the
combinatorics alone**, and the extra modes seen at special `X` (symmetric tilings) are confined to
`V`. Designing a minimal-mobility pattern is then: pick `X ∈ X ∖ V`.

**Gap targeted.** 2026 Sec. 4.4 and Fig. 12/13 report design-space dimensions as measured numbers
with no statement of what is generic (`notes/paper_2026.md` §7.5, and §7.3 claim 4 "Otherwise, the
system becomes rank-deficient" asserted without argument). 2025 Sec. 6 limitation 3: "Our current
investigation of the realizable shape space associated with each pattern is largely empirical;
developing a mathematically rigorous framework for analyzing these shape spaces remains an open
problem."

**Novelty vs field.** R5 (Jackson–Jordán, pin-collinear body-and-pin) is the structural template —
"a *specific* class of special positions can be shown to behave generically" — but for a different
degeneracy class, as the field notes say explicitly. R4/R6 are genericity statements over *all*
placements, not over a constrained linear family, which is what `X` is. This is genericity relative
to a variety, which is exactly the gap `notes/field_rigidity.md` §2(a) identifies: "any correct DOF
statement for kirigami must be a statement about the actual embedding".

**Why it might be true / sketch.** `A(X)` has entries linear in `X` and `X` is a linear subspace, so
`X ↦ A(X)` is a linear map into matrices; `rank` is lower semicontinuous on any irreducible variety,
so the maximum rank is attained on a Zariski-open dense subset, and the complement is cut out by
minors. The content is not the abstract statement but (i) computing `m*` combinatorially and
(ii) showing `V` is non-empty and contains the symmetric tilings everyone uses.

**Kill experiment.** 500 cut structures. For each, sample 30 random points of `X` (random null-space
coefficients `t_{φ_i}` from `solve_system`) and compute `m`. **If `m` is not constant across the 30
samples for at least 95% of structures, the semicontinuity claim is fine but the practical form is
dead.** Then check that the classical tilings (squares, kagome, hexagons) sit at `m > m*`, i.e. in
`V`. **If the symmetric tilings do not have strictly larger `m` than random points of their own
shape space, the "extra modes are a symmetry artefact" story is dead.**

**If it survives, the demo.** A histogram of `m` over the shape space for ≥20 cut structures showing
a sharp mode at `m*` with a thin tail at symmetric points; a table of `m*` vs combinatorics. Hero:
take a symmetric tiling, perturb inside its own shape space, and watch the parasitic modes vanish
while uniform deployability is exactly preserved.

**Risk.** `m*` may be so large (my run: 19–26 for random Delaunay patterns) that "minimal mobility"
is not a useful design target, and the interesting content collapses into idea 7's boundary term.

**Effort.** Derivation M. Code M.

---

## Idea 9 — Range-maximising Tutte auxetic embeddings from the closed-form θ_max
**Type: algorithm.**

**Claim.** Using idea 3, `θ_max(X)` is the minimum over candidate pairs of the smallest positive
root of an explicit quadratic whose coefficients are quadratic forms in `X`. Maximising `θ_max` over
the shape space `X` (i.e. over the null-space coefficients `t_{φ_i}` of Eq. 5) is then a smooth
max-min problem with exact gradients, solvable by a smooth-min surrogate plus Newton, and it beats
the paper's Eq. (9) barrier heuristic — which is a proxy evaluated **at `θ = 0` only** — by a large
margin on the achieved collision-free range.

**Gap targeted.** 2026 Sec. 4.5, Eq. (9), and `notes/paper_2026.md` §8.5: "Optimizing at `θ = 0`
only is itself a heuristic: it penalizes the *initial* velocity direction at split edges, not
collisions at the angle where they actually occur. The paper does not justify why this suffices."
Plus §15 open question 14: "Does optimizing Eq. (9) at `θ = 0` actually maximize `θ_max`, or merely
delay the first collision?" And Sec. 6 limitation 1.

**Novelty vs field.** IsoGami (F12) does collision-aware continuation with IPC but only screens a
sampled catalogue of isohedral topologies with default tile shapes ("We perform topology screening
using default tile shapes"), so it never optimises geometry for range. R18 (Liu et al. 2024) has
closed-form tests for enumerated quad pair types only. Nobody optimises deployment range over a
characterised linear design space, because until 2026 Eq. (5) there was no such space.

**Why it might be true / sketch.** With `X = X₀ + Σ φ_i t_i`, each `y_i(θ)` is bilinear in
`(X, (cos θ/2, sin θ/2))`, so each contact polynomial is `q(t; X) = α(X)t² + β(X)t + γ(X)` with
`α, β, γ` quadratic in the coefficients `t_i`. The smallest positive root is a smooth function of
`(α, β, γ)` away from the discriminant; a log-sum-exp softmin over pairs gives a differentiable
objective; Newton with PSD projection (the paper's own solver, Sec. 5.4) applies.

**Kill experiment.** Implement both. Baseline: Eq. (9) with a standard log-barrier `B` and a swept
`γ` (the paper gives neither, `notes/paper_2026.md` §8.5). On ≥20 graphs from the tiling and random
generators, measure achieved `θ_max` by the *same* independent bisection collision routine for both
methods. **If the new objective does not beat the Eq. (9) baseline by ≥ 25% median relative
improvement in `θ_max`, dead.** Runtime budget under 30 min because the closed form removes the
simulation loop.

**If it survives, the demo.** `θ_max` scatter, new vs baseline, on ≥20 graphs, with the binding
contact pair annotated. Hero: a pattern where Eq. (9) yields a nearly useless range and the new
objective opens it fully. Fabricable export: two cut sheets from the same graph, one per method.

**Risk.** The Segall baseline may already be near-optimal because the shape space is small for the
patterns that matter (few split cuts ⇒ small kernel — a tension the paper itself never discusses,
per `notes/paper_2026.md` §13.2). Then the margin will not be there. Second risk: the true binding
constraint may be a non-adjacent pair, requiring a broad phase that eats the speed advantage.

**Effort.** Derivation S (rests on idea 3). Code L.

---

## Idea 10 — Which flat-state flexes are real: a prestress certificate at θ = 0
**Type: theorem (+ algorithm).**

**Claim.** At the flat configuration a kirigami has `m(0)` infinitesimal flexes but only
`m(0⁺) < m(0)` of them extend to motions (idea 2). The obstruction is second order and computable
from the same data: a flex `ω ∈ ker A(0)` extends to second order iff the quadratic form
`Q(ω) = Σ_{self-stresses ρ} ρᵀ · (second-order pin residual of ω)` vanishes, and `Q` is an explicit
quadratic form on `ker A(0)` whose coefficients are the `A₁` block of idea 5's pencil. Concretely:
`ω` is a genuine branch direction iff `ω ∈ ker A₀ ∩ ker A₁` (up to the pencil's Kronecker structure).

**Gap targeted.** 2026 Sec. 6 limitation 2 (which flexes are real), and 2025 Sec. 6 limitation 2
("stability analysis of the deployed structure and other mechanical factors were not considered").
`notes/field_rigidity.md` §2(c) ranks second-order rigidity (R8) as the most directly transferable
machinery and notes it "certifies *rigidity*", never producing flexes — this idea uses it the other
way, as a filter on candidate branches.

**Novelty vs field.** R8 (Connelly–Whiteley) supplies the second-order/prestress framework at a
given non-generic configuration but says nothing about kirigami. R20 (Müller) is the general
higher-order machinery with no a-priori order bound and a bracket recursion that blows up with loop
count — a kirigami tiling has `|E_hinge| − |F| + 1` loops. R9 (Connelly–Servatius) warns that
order-`n` definitions break for `n ≥ 3`. The claim here is that for *this* model the pencil makes the
second-order test exact and closed form, because the constraint variety restricted to the deployment
direction is degree 2 (idea 3), so no order beyond 2 is needed along that direction.

**Why it might be true / sketch.** The pin constraint (R1) is quadratic in the configuration, and by
idea 3 the deployment path is degree 2 in `t`; so the Taylor expansion of the constraint along a
candidate branch terminates. Writing the pencil `A(t) = A₀ + tA₁` (after dividing by `cos(θ/2)`),
`ω` continues to first order in `t` iff `A₀ω = 0` and `A₁ω ∈ im A₀`; iterating gives the Kronecker
chain, and the chain length bounds the order needed. That is precisely the structure R9 says is
missing in general and that this model supplies.

**Kill experiment.** On 500 graphs: compute `ker A₀`, and for each basis flex test whether it lies in
the pencil's "deflating subspace". Independently, integrate that flex numerically from `θ = 0` with
constraint re-projection. **If a flex certified as extending fails to integrate (pin residual cannot
be driven below 10⁻¹⁰), or a flex certified as blocked does integrate, on any graph, dead.**
Cross-check: the count of certified-extending flexes must equal the measured `m(0⁺)` from V4.
**If those two numbers differ on any graph, dead.**

**If it survives, the demo.** For ≥20 graphs, the table `m(0)` / `m(0⁺)` / certified count, all
three agreeing. Hero: a pattern with two certified branches at `θ = 0` deployed side by side, and a
third "shaky" flex shown to go nowhere.

**Risk.** The single most likely failure: the pencil is singular (non-regular) for kirigami, so the
deflating-subspace computation is numerically delicate and the clean statement dissolves into
Kronecker minimal indices. This is the same risk as idea 5 and they stand or fall together.

**Effort.** Derivation L. Code M.

---

# Self-attack: hostile SIGGRAPH/TOG reviewer on my top three

**Against idea 1 (mobility operator).** *Reject.* The authors rediscover, in coordinates, the
elementary fact that a planar body-and-pin framework's compatibility condition is a circulation
condition on the joint graph — this is instant-centre kinematics, taught in every machine-theory
course, and the "three collinear pins ⇒ mobile" corollary the authors offer as a validation is
Aronhold–Kennedy (1872). Reducing the rigidity matrix from `3|F|` to `|F|` columns is a linear
algebra convenience, not a contribution; the resulting number is still computed by a numerical rank,
exactly as Chen, Choi and Mahadevan did in PNAS in 2020, and exactly as IsoGami did in 2026 with
eigenvalues of `JᵀJ`. There is no theorem here: no combinatorial formula for `dim ker A`, no
characterisation of when it exceeds 2, nothing a designer can act on before running an SVD. The one
genuinely new sentence — that Segall's design matrix and the mobility matrix are two slots of one
bilinear form — is an observation about notation with no consequence drawn from it. The paper is a
reformulation in search of a result.

**Against idea 3 (conic deployment path).** *Reject.* That a rigid mechanism driven by a single
angle has trajectories rational in `tan(θ/2)` is the Weierstrass substitution and is the standard
device in every algebraic treatment of planar linkages; the authors' "theorem" is that a composition
of rotations by `±θ/2` is linear in `(cos θ/2, sin θ/2)`, which is true by inspection of the forward
kinematics the previous paper already wrote down. Segall et al. did not state it because it is not
worth stating. The claimed payoff, a closed-form `θ_max`, is a quadratic root per candidate pair —
but the authors must still enumerate `O(n²)` pairs, so they have replaced a bisection whose cost
they never measured with an enumeration whose cost they also do not measure, and the paper reports
no case in which the closed form finds a collision the simulation missed. Meanwhile the actual open
problem named in Segall et al. Sec. 6 is to predict *favourable deployment behaviour from the
embedding* — a structural characterisation, not a faster evaluator of the same quantity. This
submission answers an easier question and claims the harder one.

**Against idea 2 (flat state is a branch point).** *Reject.* The central empirical claim rests on a
numerical rank taken at `θ = 10⁻⁶` and at `θ = 0`, with a fixed relative tolerance the authors never
justify and never vary; a rank that "drops" between two nearby configurations is the classic
signature of a badly conditioned SVD, not of a stratification, and the authors' own data contain a
counterexample they cannot explain (the triangular tiling, where the drop does not occur). Even
granting the drop, "the flat state is singular" is unsurprising: it is the configuration where every
pair of hinge-edge duplicates coincides, so of course the linearisation degenerates — this is why
the mechanisms literature has warned since Connelly and Servatius (1994) that special positions must
be handled with higher-order tools, which the authors invoke (their idea 10) but do not deliver
here. The consequence drawn for Segall et al.'s Eq. (9) — that optimising collisions at `θ = 0` is
evaluated at a singular configuration — is asserted, never demonstrated to change any result. Show
me one pattern where Eq. (9) demonstrably fails *because* of this, or the observation is a curiosity.

---

# Final ranking

| # | Idea | One-line justification |
|---|---|---|
| 1 | **3 — conic deployment path + closed-form contact angles** | The one claim already verified to 10⁻¹⁵ on 21/21 graphs, it converts the paper's own named open problem from simulation to algebra, and every other idea of mine depends on it. |
| 2 | **1 — mobility operator `m = dim ker A − c`** | Verified against the full rigidity matrix on 21/21 graphs, answers Sec. 6 limitation 2 with a number instead of a "sensitivity analysis", and reproduces Maxwell and Aronhold–Kennedy as checks. |
| 3 | **2 — flat state is a branch point, mobility drops at `θ = 0`** | The most surprising result I found (20/21 graphs, e.g. `10 → 1`), and it lands directly on Eq. (9), which the paper solves at exactly that singular configuration. |
| 4 | **5 — bifurcation angles as pencil eigenvalues** | Turns "where does the mechanism bifurcate" into a determinantal root-finding with no path traversal; highest upside, real risk that the pencil is singular. |
| 5 | **7 — 2-core decomposition of mobility** | Cheap, exact, and its corollary contradicts the universally repeated and never-proved claim that rotating squares is 1-DOF. |
| 6 | **6 — space of admissible non-uniform θ fields** | The exact object Sec. 6 limitation 2 asks for, with dimension `= m`; downgraded only because its content is mostly a corollary of ideas 1 and 7. |
| 7 | **9 — range-maximising embeddings** | The strongest *algorithm* deliverable with a clean baseline, but its margin depends on a shape space that is small precisely when the pattern is good. |
| 8 | **4 — conformality theorem** | Certainly true and it converts a stated empirical claim into a two-line proof, but it is a lemma, not a paper. |
| 9 | **10 — prestress certificate for flat-state flexes** | The right way to close idea 2, but it stands or falls with idea 5's pencil regularity and is the largest derivation effort. |
| 10 | **8 — genericity on the shape space** | Structurally sound but the abstract half is routine algebraic geometry, and my measured `m*` values (19–26) suggest the practical form may be vacuous. |

---

# What I could not check

- **`rank(L) = H` in the periodic case.** Scout-c's `notes/field_tutte.md` §0.4 predicts
  `rank(L) ≤ H − 1` when there is no boundary; my 21 runs all used `BoundaryMode::Fixed`, so they
  neither confirm nor refute it. In the bounded case `rank(L) = H` held 21/21. I did **not** claim
  the rank result as my own idea because Scout-c derived it independently and first. The one
  increment I would add is a consequence they do not state: with a *fixed* boundary the sum of all
  hole rows gives the necessary feasibility condition `Σ_{v ∈ ∂} (indeg_h(v) − outdeg_h(v)) x_v = 0`
  on the prescribed boundary data, which the inverse-design pipeline (2026 Sec. 5.3, boundary fixed
  to a circle) does not check. Translation-invariant because `Σ_{v∈∂}(indeg−outdeg) = 0`. Untested.
- **`notes/field_tutte.md` was incomplete while I worked.** I read its §0 and §1; §2 (proof routes)
  and §3 (search log) were not yet written, so my novelty claims against the Tutte/discrete-geometry
  literature are weaker than those against the rigidity and kirigami tables.
- **Whether `A(θ)` is a regular or singular pencil.** Ideas 5 and 10 both assume this can be handled;
  I did not compute a Kronecker form. This is the single largest unverified assumption in the file.
- **The triangular-tiling exception to idea 2** (mobility 9, constant through `θ = 0`). I have no
  explanation. It may be that its hinge graph has no cycles through collinear pins at the flat state.
- **Borcea–Streinu full texts** (R16, R17), flagged as abstract-only in `notes/field_rigidity.md`.
  *Periodic tilings and auxetic deployments* (2020) is the nearest topical neighbour to idea 4 and
  I could not read it. The orchestrator has checked the 2020 paper (STATE.md) and judged it not a
  threat; I relied on that.
- **R25** (PNAS rigid-origami mobility framework) is second-hand in the field notes (HTTP 403).
  If its fourth-order compatibility pipeline is as described, it is the nearest prior art to idea 10
  and must be read before that idea is committed to.
- No claim in this file about `θ_max` has been compared against an independent collision routine,
  because `code/src/core/collision.cpp` does not exist yet (the repo's `CMakeLists.txt` references it
  and therefore does not configure).
