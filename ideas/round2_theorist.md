# Round 2 — Ideator, persona: THEORIST (discrete geometer / rigidity theorist)

**What I want.** A theorem. Specifically an *existence / non-existence criterion* for usable
(deployable, self-intersection-free) embeddings in the Tutte auxetic shape space `𝕏(G, σ)`, stated in
terms of the combinatorics of `(G, σ)` and not in terms of "we sampled and found nothing". Round 1
measured emptiness (K1a 0/200, K5 0/400, K6 0/31 certified); round 2 must *prove* it, and prove the
complementary statement for the graph classes where the space is not empty.

**Files read.** `specs/common_preamble.md`, `specs/ideator.md`, `specs/ideator_round2.md`, `STATE.md`
(F1–F30, D5–D9, dead ends), `derivations/core.md` §0, T1, T2, T3 headers, T4.4, T5.1, T5.2b, T5.3, T7,
`results/kill/KILL_REPORT.md` (K1a, K5, K2a, K2b, K1c, K3a, F23), `results/kill/k6/k6_final.csv` and
`k6_full.csv` headers + first rows, `code/src/method/{deploy_basis.hpp, zero_plus.hpp}`,
`notes/field_tutte.md` (rows 1–20 + Routes A/B/C), `notes/field_kirigami.md` (R1–R24),
`notes/field_rigidity.md` (R1–R26). I did **not** read the 2025/2026 paper text in this session; every
paper citation below is taken from the F-facts and the reader notes as recorded in `STATE.md`, and I
mark the two places where that matters.

---

## 0. The three algebraic facts everything below is built on

I restate them because each of my ideas is a theorem *about* one of them, and because two of them are
mine rather than quotations.

**(0.a) The split-gap coefficient, in closed form** (`derivations/core.md` T5.3, T1.B). For a split
edge `e = {a, b}` between faces `f, g` (`σ_f = σ_g`), with `d = x_b − x_a` and `Δu = u_g − u_f` the
difference of face potentials, the two duplicates of `e` stay *the same vector* and differ by the pure
translation `2 sin(θ/2) J Δu`. The signed gap is the harmonic `h(θ) = r sin θ + p (1 − cos θ)` with

```
        r = ⟨d, Δu⟩ ,      p = − σ_f det(d, Δu) ,      p + q = 0 .
```

`code/src/method/zero_plus.hpp` computes `q_e = det(dS_e, d_e)`, a positive multiple of `r`. **The
entire measured emptiness is the statement `min_e q_e < 0`** (K6 `min_q0 = −21.5` on graph 0,
`inward_frac0` 0.24–0.41).

**(0.b) The vertex parity identity** (mine; one line). At an interior vertex `v` the faces around `v`
form a closed cycle in the dual; an edge of that cycle is a hinge edge exactly at a sign change of
`σ`, and the number of sign changes around a closed cycle is **even**. Hence

```
        deg_hinge(v)  is even ,      deg_split(v) ≡ deg(v)  (mod 2)      at every interior v.
```

This is `STATE.md` U9 with a proof, and it is the seed of Ideas 1 and 10.

**(0.c) The barycentric row is a *balance* condition, not a Tutte condition** (mine). In the pure
case Eq. (2) at interior `v` reads `k_v x_v = Σ_{u ∈ In(v)} x_u`, i.e.

```
        Σ_{u ∈ In(v)} ( x_u − x_v )  =  0 .                                   (BAL)
```

F14 reads this as "`x_v` is the mean of half its star, so Tutte does not apply", and stops. The
*kinematic* reading is the useful one: by T1.5 the `k_v` copies of `v` separate at `0⁺` with
velocities proportional to `J(x_u − x_v)`, one per in-edge. (BAL) says those velocity vectors **sum to
zero**, hence (for `k_v ≥ 3`) positively span the plane, hence no copy can be pushed into its
neighbours' material: *the barycentric row is precisely what prevents the `0⁺` collapse.* A split edge
at `v` deletes a term from (BAL) without deleting the corresponding motion, and the balance is broken.
Ideas 9 and 2 are built on this.

---

## 1. The ten ideas

### Idea 1 — Parity forces split cuts: the odd-vertex deficiency theorem

**1. Title / type.** *Odd-degree deficiency: a parity obstruction to deployability* — **theorem**.

**2. Claim.** For every `σ` on a disk-topology planar `M`, `E_split(σ)` is a **T-join** of
`T = {v interior : deg(v) odd}`; hence `|E_split(σ)| ≥ |T| / 2` for **every** `σ`, and
`E_split(σ) = ∅` for some `σ` **iff** every interior vertex of `M` has even degree (equivalently: the
interior dual is bipartite; equivalently, for a triangulation, `M` is Heawood 3-vertex-colourable).
Combined with the measured law `Θ_max = min(min_v β_v, π)` **iff** `E_split = ∅`, this gives a graph
class with a *proved* deployment guarantee and, in the complementary direction, a lower bound on the
number of independent split-gap sign conditions that scales with the odd-degree count.

**3. Gap targeted.** 2026 sells the method as working on *arbitrary* planar graphs (its whole framing
against 2025, which per `notes/field_kirigami.md` R21 assumes "2-colorable planar tilings (even
valency)"). 2026 §4.4 discusses `E_split = ∅` only as "2-colorable with alternating orientations"
(F9) — a restatement, not a criterion on `M`. Nowhere in either paper is there a statement of the form
"this graph *cannot* be cut without split edges", nor a bound on how many split edges are forced.
2026 §7 limitation (i) (collision range uncharacterised) is precisely where the consequence bites.

**4. Novelty vs field.** (i) `field_kirigami` R21 (Segall 2025) *assumes* even valency; it does not
prove the equivalence, does not bound `|E_split|` when valency is odd, and has no notion of split cut
at all. (ii) R12 Dang 2021 is quads-only (all interior degrees are 4, so `T = ∅` automatically — their
whole class is the trivial case of my theorem, which is *why* they get an iff). (iii) `field_tutte`
row 11 (Whiteley's hypergraph matroid) is the standard "rows indexed by sets" machinery and produces
generic rank counts, not a parity obstruction on `σ`. The Heawood-3-colourability link (a triangulation
is 3-colourable iff Eulerian) appears in neither paper nor in any row of the three field tables.

**5. Why it might be true / sketch.** (0.b) gives `deg_split(v) ≡ deg(v) (mod 2)` at interior `v`, so
the subgraph `E_split` has odd degree exactly on `T` — that is the definition of a T-join, and every
T-join has `≥ |T|/2` edges. For the iff: `T = ∅ ⇒` every vertex cycle of the interior dual is even; a
cycle of the interior dual encloses a set `S` of interior vertices and has length
`Σ_{v ∈ S} deg(v) − 2|E(S)| ≡ Σ_{v∈S} deg(v) (mod 2)`, so *all* cycles are even, the interior dual is
bipartite, and its proper 2-colouring is a `σ` with `E_split = ∅`. Conversely `E_split = ∅` forces
`deg_split(v) = 0`, hence `deg(v)` even. For triangulations the interior dual is 3-regular and
bipartite iff the primal is Eulerian iff (Heawood) it is 3-colourable.

**6. Kill experiment.** Driver: extend `code/apps/kill_k1a.cpp`'s population loop (200 graphs,
`F ∈ [101, 793]`) plus the 8 authored tilings of `kill_common.hpp::deployable_population()`. For each
graph compute `|T|`, then `|E_split|` under `σ_mc`, `σ_def` (the 200 saved
`results/kill/k5/sigma/*.json`) and 200 uniformly random `σ`. **Dead if any `(graph, σ)` has
`|E_split| < ⌈|T|/2⌉`, or if any graph with `T = ∅` yields `E_split ≠ ∅` under the interior-dual
2-colouring, or if any graph with `T = ∅` fails
`POS ∧ NOOVERLAP ∧ NOROOT` with `Θ_max = min(min β, π)`.** For the `T = ∅` half I need even-degree
generators: a quad grid, a hexagonal-dual (triangular) lattice patch, and Eulerian triangulations
obtained by taking any triangulation's barycentric subdivision. Counting only: seconds. The
certificate half: `kiri_analyze` at ~10 s per graph, ≈ 10 min for 40 graphs.

**7. If it survives, the demo.** A single scatter: `|T|/|V_int|` on the x-axis, certified `Θ_max` on
the y-axis, over ≥ 200 graphs — random Delaunay/Voronoi in a cloud at `Θ_max = 0` with `|T|/|V_int| ≈
0.5`, authored tilings and the constructed Eulerian family at `Θ_max = min β` with `|T| = 0`. Hero:
the *same* point set, triangulated two ways (one Eulerian by subdivision, one not), one deploying to
`min β` with a laser-cuttable export, the other certified dead.

**8. Risk.** The equivalence half is close to folklore (planar `G` Eulerian ⟺ `G*` bipartite is a
textbook fact) and 2025 already restricts to even valency — a referee can call this "the known
hypothesis, restated". The defence has to be the *quantitative* half (`|E_split| ≥ |T|/2` plus the
fatality curve of Idea 2/7), which is not folklore. A second risk: my interior-dual cycle argument
assumes disk topology; with holes in `M` there are dual cycles not generated by vertex cycles.

**9. Effort.** Derivation **S** (the proof above is essentially complete). Code **S**.

---

### Idea 2 — A Farkas certificate of emptiness for the `0⁺` feasibility system

**1. Title / type.** *Certified emptiness of the usable region by aggregated quadratics* —
**theorem** (with a computable certificate).

**2. Claim.** Usability requires `q_e(t) > 0` on every split edge and `μ_c(t) > 0` on every corner
incidence, both **quadratic** polynomials in the design coordinates `t ∈ R^{2|E_split|}`
(`zero_plus.hpp`). Emptiness of `{t : q_e(t) > 0 ∀e, μ_c(t) > 0 ∀c}` is *certified* by any
`λ, ν ≥ 0`, not all zero, with `Σ λ_e q_e(t) + Σ ν_c μ_c(t) ≤ 0` for all `t` — i.e. with the
aggregated quadratic form negative semidefinite and its aggregated affine part non-positive. **Claim:
such a certificate exists on random planar graphs, with `λ` supported on a set of size `O(1)`
(a "frustrated cluster"), and does not exist on the authored tilings.**

**3. Gap targeted.** 2026 §7 limitation (i): "the collision-free range is not characterized from the
embedding, only via FK simulation — open problem". A Farkas/S-procedure certificate is the exact
object that closes it in the negative direction: a *proof*, checkable in `O(1)` per cluster, that no
member of the shape space deploys. 2026 has no negative certificate of any kind; Eq. (9) (F6) is a
penalty, and K1c shows it re-closes on 9.47 % of designs.

**4. Novelty vs field.** (i) `field_tutte` row 14 (Spriggs, parallel-morph counterexamples) warns that
staying in the linear solution space does not imply staying embedded, but gives no certificate. (ii)
IsoGami (F12) detects mobility numerically via zero eigenvalues of `JᵀJ` and resolves contact with IPC
— a simulation, so it can only ever report "this one failed". (iii) `notes/screen_r1.md` records that
"validity certificate in the null space" was searched for and **NOT FOUND**. The S-procedure itself is
textbook; what is new is that the kirigami `0⁺` constraints are exactly quadratic (T1.5 ⇒ `q` and `μ`
quadratic), which is what makes the certificate finite and cheap.

**5. Why it might be true / sketch.** Write `q_e(t) = tᵀ A_e t + b_eᵀ t + c_e`. If
`Σ λ_e A_e ⪯ 0` and the aggregate is non-positive at its unconstrained maximiser, then
`max_t Σ λ_e q_e ≤ 0`, so some `q_e ≤ 0` everywhere: emptiness, proved. The reason to expect a *small*
support: by (0.c) the `0⁺` opening velocities at a vertex sum to zero when all its edges are hinges;
when `k` split edges meet at one vertex they carry `k` velocity terms that were deleted from (BAL), and
the remaining ones must absorb the imbalance — so a **single vertex star** should already supply a
`λ`. That is the local certificate to look for first: `λ` = indicator of the split edges in one star.

**6. Kill experiment.** Driver: new `kill_k8.cpp` on top of `zero_plus_form()` (already returns `GS`,
`DD`, hence `A_e, b_e, c_e` by three evaluations per monomial, or directly by differentiating the
bilinear forms). On the 60-graph K6 subset: (a) test uniform `λ ≡ 1` — is `Σ_e A_e ⪯ 0`? (b) for every
interior vertex with `deg_split ≥ 2`, test the star-supported `λ`. (c) a greedy search over supports of
size ≤ 6. **Dead if no certificate is found on any of the 60 graphs** (then emptiness stays a
measurement, and Idea 2 collapses into K6's negative result). Also **dead if a certificate IS found on
an authored tiling** (it would then be proving something false, i.e. my forms are wrong). Cost:
`|E_split| ≈ 400` eigen-decompositions of `2|E_split| × 2|E_split|` matrices is too much — restrict to
star supports first, dimension ≤ 12 after projecting onto the star's active columns; seconds per graph.

**7. If it survives, the demo.** A table over ≥ 20 random graphs: certificate support size, the
certifying vertex, the aggregated form's largest eigenvalue, wall time — against K1a's `10⁴` samples
per graph that prove nothing. Hero figure: one Voronoi patch with the certifying star highlighted, and
the sentence "no point of this 736-dimensional shape space deploys, and here is the 8×8 matrix that
proves it".

**8. Risk.** The most likely failure is that no small-support certificate exists because the
obstruction is genuinely global; then only a full SDP over `λ ∈ R^{400}` would decide it, which is
neither cheap nor a nice theorem. Second risk: `μ_c` is only *piecewise* quadratic (the max/min branch
in `zero_plus.hpp`), so the aggregate over corner incidences is not a single quadratic — the
certificate may have to be built from `q` alone, which may be too weak (K6 pass 3 says `q > 0` is
achievable and the binding contact is then `vertex-edge`).

**9. Effort.** Derivation **M**. Code **M**.

---

### Idea 3 — The first-order area identity, and the expansion floor

**1. Title / type.** *Areal expansion of a uniform deployment is a closed form in the flat design* —
**characterization**.

**2. Claim.** For every `X ∈ 𝕏`, the area enclosed by the deployed outer contour satisfies

```
      d/dθ |_{θ=0}  Area(θ)   =   ½ Σ_{e ∈ E_hinge} |x_{dst(e)} − x_{src(e)}|²   +   Σ_{e ∈ E_split} ⟨d_e, Δu_e⟩ ,
```

exactly (not asymptotically). Consequently **every usable design obeys an expansion floor**: since
usability forces `⟨d_e, Δu_e⟩ > 0` on every split edge (0.a), the areal expansion rate is at least
`½ Σ_{E_hinge} |e|² / Area(0)`, a quantity computable from `(G, σ, X)` alone; and any target requiring
a smaller expansion rate is **unreachable by that pattern**. Periodic version: the right-hand side
equals `Area(cell) · d(det J)/dθ|_0`, a closed-form link between K7's periodic Jacobian and the hinge
edge lengths.

**3. Gap targeted.** 2025's open problem, `notes/paper_2025.md` Sec 13–15 as recorded in
`specs/ideator_round2.md` §What-we-want (c): "the realizable shape space per pattern (triangles cannot
close a hemisphere, Fig F.4)" — stated there as an empirical observation with no formula. 2026 §5.1's
periodic Jacobian analysis (F7) computes `J(θ)` numerically per pattern and never relates it to the cut
combinatorics.

**4. Novelty vs field.** (i) `field_kirigami` R17 (Dudte–Choi–Mahadevan 2023) reduce quad kirigami to
four-bar recursions and compute achievable shapes numerically, with no expansion bound. (ii) R5/R7
(Konaković et al.) design triangular auxetic linkages to hit a target surface and enforce the
expansion limit as a *per-triangle* numeric bound from the fixed linkage family, not as a formula in
the cut graph. (iii) `notes/screen_bundle.md` verdict S7 records that "hole count / area harmonic /
second closed angle" was searched and **NOT FOUND**; the optimizer persona has the per-hole harmonic
`A(θ) = q(cos θ − 1) + r sin θ`, but not the *global* identity nor the floor.

**5. Why it might be true / sketch.** Faces are rigid (T1.D), so the enclosed area is
`Σ_f area(f) + Σ_{holes and notches} area`. At `θ = 0` every hole is a degenerate slit, so to first
order its area is the sum of the openings of the cut segments on its boundary. A hinge edge opens a
wedge of apex angle `θ` and legs `|d|` (T1.C), area `½|d|² sin θ`. A split edge opens a parallelogram
with sides `R_f d` and `2 sin(θ/2) J Δu` (T1.B), area `2 sin(θ/2) det(R_f d, J Δu) → θ ⟨d, Δu⟩`.
Differentiate at 0 and sum. Every cut edge is owned by exactly one preimage (F11), so no term is
double counted — which is exactly where F11 earns its keep.

**6. Kill experiment.** Driver: a 40-line addition to `kiri_analyze`. On the 8 authored tilings + 30
random designs from the K1a population, compute the polygon area enclosed by the deployed outer
contour at `θ = 10⁻⁵` and `θ = 2·10⁻⁵` by `deploy()`, Richardson-extrapolate the derivative, and
compare with the right-hand side assembled from `cut.hpp` edge lists and `zero_plus`'s `Δu`.
**Dead if the relative deviation exceeds `10⁻⁶` on any design.** Then the second half: on the K1a 200,
check that `Σ_{E_split} ⟨d, Δu⟩ < 0` on every graph (it must be, since some `q_e < 0` on all of them) —
and report how often the *aggregate* is negative, since an aggregate-negative graph is non-deployable
by a **one-number** test. Under 5 minutes.

**7. If it survives, the demo.** Reproduce 2025 Fig F.4's failure ("triangles cannot close a
hemisphere") as an arithmetic inequality: expansion floor of the triangle pattern versus the areal
ratio demanded by the hemisphere, both computed in closed form, on the 16 released patterns of
`baseline/kirigami_tessellations` against their 4 target meshes — a 16×4 table of predicted
feasible/infeasible against the paper's reported successes and failures.

**8. Risk.** The first-order identity is almost certainly true and almost certainly *easy* — a referee
may call it a computation, not a theorem. The floor is the interesting half, and its risk is that the
2025 pipeline's deployment is **non-uniform**, so a single `θ` and hence a single expansion rate may
not be the right object; the floor would then bound only the uniform sub-family.

**9. Effort.** Derivation **S**. Code **S**; the 2025 comparison **M**.

---

### Idea 4 — Injectivity of the Eq. (6) projection is an index condition on half-stars

**1. Title / type.** *A local index criterion for the validity of the Tutte auxetic projection* —
**theorem / characterization**.

**2. Claim.** In the pure-hinge case with a fixed convex boundary, `X₀` is the unique solution of
`x_v = (1/k_v) Σ_{u ∈ In(v)} x_u`; the in- and out-edges **alternate** around every interior star, so
`X₀` is a convex-combination map with weights that are non-reciprocal and zero on half of each star.
`X₀` is injective **iff** every interior vertex has Gortler–Gotsman–Thurston index 1, and the index can
drop below 1 only at vertices where an *out*-neighbour's image crosses the cone spanned by the
in-neighbours' images. This yields an `O(deg)` per-vertex certificate of validity for Eq. (6), and a
counterexample family showing that no hypothesis on `M` alone (3-connectivity, nodal 3-connectivity)
suffices.

**3. Gap targeted.** 2026 §6/§7: the shape space has no validity certificate; F17 quantifies the
consequence (Eq. (6) `X₀` self-intersecting on 57/57, inverted faces on 142/200, face overlap on
200/200, K1a). The paper simply solves Eq. (6) and proceeds. F14 already observes "Tutte/Floater do not
directly apply" and stops there; nobody has said *what does* apply.

**4. Novelty vs field.** (i) `field_tutte` row 2 (Floater 2003 Thm 4.1) requires **positive** weights
on the *whole* neighbourhood — exactly the hypothesis that fails. (ii) row 4 (GGT 2006) is the right
machinery and explicitly handles asymmetric weights `w_ij ≠ w_ji` (§6), but their Cor. 2.8 needs
positivity, and their §6 conclusion says the constrained case "remains to be seen"; zero half-edge
weights are outside their notion of a *valid* one-form. (iii) row 16 (Ó Dúnlaing) gives nodal
3-connectivity as the right hypothesis for the all-neighbour barycentric map — I expect to show, by
counterexample, that it is *not* sufficient here, which is itself the result.

**5. Why it might be true / sketch.** Alternation: around an interior `v` the faces alternate in `σ`
(pure case), and by 0.2 the hinge direction is fixed by the `σ` pair, so consecutive edges alternate
`in`/`out`; hence `In(v)` and `Out(v)` interleave, `k_v = deg(v)/2`, and `x_v ∈ int conv(In(v))` by
(BAL). Local injectivity at `v` is the statement that the star's image winds once; with only `In(v)`
constrained, the `Out(v)` images are free, and the winding number drops exactly when some
`x_w, w ∈ Out(v)`, falls on the wrong side of the fan. The GGT index theorem
`Σ_v ind(v) + Σ_f ind(f) = 2 − 2g` then converts "every vertex is a wheel" into global injectivity, as
in their Thm 3.13. A one-parameter counterexample: keep `In(v)` fixed (so (BAL) holds and `X₀` is
unchanged) and slide one out-neighbour across the fan — the linear system does not see the move, so the
same `X₀` is both injective and not, proving no criterion on `X₀` alone can work and forcing the
statement to be about the *pair* (in-star, out-star).

**6. Kill experiment.** Driver: extend `kill_k1a.cpp`. On its 200 graphs compute, for `X₀`, the
per-vertex winding of the star's image (sum of signed corner angles / 2π) and the inverted-face and
overlap sets already computed there. **Dead if any graph has all interior indices equal to 1 and still
has an inverted face or a face–face overlap** (then index-1 is not sufficient and the criterion is
wrong), or if the 58 graphs with `n_inv(X₀) = 0` but overlapping `X₀` all have all-index-1 vertices
(then the index misses the failure mode that actually matters, and the theorem is about the wrong
predicate). Under 2 minutes on the cached population.

**7. If it survives, the demo.** A validity oracle for Eq. (6) that is `O(|E|)` and exact, replacing
K1a's `10⁴`-sample Monte Carlo; on ≥ 200 graphs, a confusion table of index-certificate vs the exact
overlap test, expected to be `0` false positives. Hero: the sliding-out-neighbour counterexample as a
two-panel figure with identical `X₀` linear systems and opposite validity.

**8. Risk.** The index condition is necessary and local; **global** injectivity needs the boundary
hypothesis of GGT Thm 3.13 (turning numbers, reflex extrema are wheels), and `M'` has many boundary
components (one per hole), so the hypothesis may simply fail and leave me with a necessary condition
only. That is still publishable as a *certificate of invalidity*, which is what K1a needs, but it is
not the iff I am claiming.

**9. Effort.** Derivation **M/L**. Code **S**.

---

### Idea 5 — The split-gap sign is a 1-ring quantity, so `σ` selection is a bounded-arity CSP

**1. Title / type.** *Locality of the split-gap sign and a local orientation rule* —
**theorem + algorithm**.

**2. Claim.** For a split edge `e = {a, b}` between `f` and `g`, the potential difference `Δu_e` is
computable along **any** `Γ`-path from `f` to `g` (path-independence is `X ∈ 𝕏`); I claim there is
always such a path inside `star(a) ∪ star(b)` unless the split forest has a component spanning both
rings, and that on the measured populations the required radius is `≤ 2` on a large majority of split
edges. Consequently `sign(q_e)` depends only on `σ` restricted to a bounded neighbourhood, and "every
split cut opens" is a **constraint-satisfaction problem of bounded arity** on `σ` — solvable by local
search with a guarantee on bounded-treewidth neighbourhoods, and reducible to 2-SAT when both endpoints
have degree 3.

**3. Gap targeted.** 2026 Eq. (1) / §4.2: `σ` is chosen by a Goemans–Williamson max-cut relaxation on
the dual, "then split by a diameter" (F3) — an acknowledged heuristic with no relation to deployment.
K5's own verdict, verbatim in `KILL_REPORT.md`: "no purely combinatorial objective on `σ` will produce
a deployable random graph. The objective has to carry a deployment term."

**4. Novelty vs field.** (i) 2026 Eq. (1) is max-cut; my objective is a different problem class
entirely (CSP over local sign predicates), not a reweighting of theirs. (ii) `field_kirigami` R9
(Chen–Choi–Mahadevan 2020) search cut *topology* at fixed geometry by simulated annealing on a
different objective (target shape), with no locality theorem. (iii) The optimizer persona's Idea 5
("split cuts cancel residuals: max-cut is the wrong objective") proposes a geometry-aware `σ` objective
but has no locality statement, so its objective is global and its search unguided.

**5. Why it might be true / sketch.** `u` is a potential on `Γ` with increments
`u_g − u_f = σ_g x_{src(e)}` (T1 Step 5). Along a `Γ`-path `f = h_0, …, h_m = g` of even length, the
signs alternate (`Γ` is bipartite, 0.5), so
`Δu = σ_g Σ_i (−1)^i x_{src(e_i)}` — a telescoping sum of **hinge-source positions**. For the shortest
possible case `m = 2` (a common `Γ`-neighbour `h`) this collapses to
`Δu = σ_g ( x_{src(h→g)} − x_{src(f→h)} )`: the vector between two hinge points, and
`q_e ∝ ⟨d_e, Δu_e⟩` becomes a statement about the *local* geometry only. A common `Γ`-neighbour exists
iff `f, g` lie in a triangle of the dual, i.e. at a degree-3 vertex — so the clean case is exactly
triangulation-dual patterns; for higher degree the path must wrap the star of `a` or of `b`, giving
arity `deg(a) + deg(b)`.

**6. Kill experiment.** Driver: a 60-line addition to `kill_k5.cpp` (it already builds `σ_mc`,
`σ_def`, `Γ`, and the split lists). For every split edge on the 200 K1a graphs compute (i) the
`Γ`-distance between `f` and `g`, (ii) the smallest `k` such that a connecting `Γ`-path exists inside
the `k`-ring of `{a, b}`. **Dead if the median over graphs of the 95th percentile of `k` exceeds 2**
(unbounded arity ⇒ no CSP, and the algorithm half dies with it). Then verify the `m = 2` closed form
against `zero_plus_q()` to `10⁻¹²` on every split edge where the distance is 2. Under 5 minutes.

**7. If it survives, the demo.** A `σ`-chooser: local search on the CSP, on the same 200 graphs, with
the metric `fraction of split edges with q_e > 0` — baselines `σ_mc` (measured `inward_frac0 ≈ 0.41`)
and `σ_def` (`≈ 0.24`). The paper claim is not "we now deploy random graphs" (Idea 2 says we cannot)
but "the sign structure is local, here is the exact arity, and here is how far a local rule gets".

**8. Risk.** Most likely false as stated: the split forest is a *forest*, and a long split component
separates `f` from `g` in `Γ` by exactly the distance the component spans, so on graphs where
`σ_def` doubles the split count (K5: 126 → 321) the arity should blow up. If so, the honest result is
the negative one — "the sign is **not** local, hence no local orientation rule can exist", which is
still a theorem and directly explains K5's failure.

**9. Effort.** Derivation **S**. Code **S** for the measurement, **M** for the CSP solver.

---

### Idea 6 — The tangent cone of the non-uniform deployment variety at the flat state

**1. Title / type.** *Which angle fields close every hole: the flat state is a singular point, and its
tangent cone is `ker A_c`* — **theorem**.

**2. Claim.** Let each face carry its own rotation `a_f` and set `z_f = e^{i a_f}`. Hole closure for a
**general** (non-uniform) angle field is `Re/Im` of a fixed complex-linear system `W z = 0` with
`|z_f| = 1`; the deployment variety is `V = {W z = 0} ∩ T^{|F|}`, the flat state `z ≡ 1` lies in `V`,
and the tangent cone of `V` at `z ≡ 1` is `ker A_c` (the `cos` half of the mobility pencil
`A(θ) = cos(θ/2) A_c + sin(θ/2) A_s`, T2.3). A tangent direction `ω ∈ ker A_c` extends to a real branch
iff it satisfies the second-order condition `ωᵀ H_z ω ∈ range A_c` for the (explicit) Hessian; the
uniform branch `ω = σ` always does (T2.2). **The number of real branches through the flat state is the
number of solutions of that quadratic on the projectivised `ker A_c`,** which is the precise form of
"the flat state is a branch point" (U5) and of 2026 §7 limitation (ii).

**3. Gap targeted.** 2026 §7 limitation (ii) verbatim (F8): "non-uniform hinge-angle deployment / DOF
of the deployment space — could be achieved by sensitivity analysis". Sensitivity analysis at the flat
state is exactly the *wrong* thing to do, because T2.3 shows `rank A(θ)` attains its minimum at
`τ = 0`. 2025 shows a non-uniform deployment (Fig. 10) with no theory.

**4. Novelty vs field.** (i) The optimizer persona's Idea 8 already has "the exact non-uniform
deployment variety is a linear section of a torus" — **I am not claiming that part**; my contribution
is the singularity analysis at `z ≡ 1` (tangent cone, second-order extension, branch count). (ii)
`field_rigidity` R8 (Connelly–Whiteley 1996, second-order rigidity) and R9 (Connelly–Servatius, cusps)
are the correct general tools and say nothing about this variety. (iii) R20 (Müller, higher-order local
mobility) is the systematic machinery; applying it here is new because `W` is explicit and sparse.

**5. Why it might be true / sketch.** Closure across a hinge at `v` between `f, g` gives
`t_g − t_f = (R(a_f) − R(a_g)) x_v`, which is **linear in `(cos a_f, sin a_f)`**; summing around a
cycle of `Γ` kills `t` and leaves `Σ_{e ∈ z} ± [(c_f − c_g) x_v + (s_f − s_g) J x_v] = 0`, i.e. a fixed
real-linear system in `(c, s) ∈ R^{2|F|}` intersected with `c_f² + s_f² = 1`. Complexify:
`z_f = c_f + i s_f` and the two real equations per cycle combine into one complex equation, hence `W`.
Linearising at `z = 1` (`c = 1, s = 0`) leaves only the `s`-terms: `Σ ± (s_f − s_g) J x_v = 0`, which is
`A_c` up to the factor `J`. Second order: `c_f = 1 − a_f²/2`, feeding the `c`-terms back gives the
quadratic obstruction.

**6. Kill experiment.** Driver: new `kill_k9.cpp` reusing `mobility.cpp`'s cycle-basis assembly. On the
8 authored tilings + 12 K3a Delaunay patches, assemble `A_c`, compute `dim ker A_c`, then for 50 random
`ω ∈ ker A_c` test the second-order condition by projecting the explicit quadratic onto
`coker A_c`. **Dead if the obstruction vanishes identically on every graph** (then every tangent
direction extends, the flat state is a smooth point of the right dimension, there is no branch
structure, and Idea 6 is just optimizer Idea 8 again). Also dead if `σ ∉ ker A_c` on any graph
(a contradiction with T2.2 that would mean my `A_c` is wrong). 20 graphs × dense kernel at
`|F| ≤ 500`: under 10 minutes.

**7. If it survives, the demo.** For each of the 16 released 2025 patterns, the branch count at the
flat state and the identification of which branch 2025's Fig. 10 non-uniform deployment follows. Hero:
a deployment-unfriendly tiling (`Θ_max = 0` for uniform deployment) shown to have a *second* real
branch with positive range — a pattern that cannot be deployed uniformly but can be deployed at all.

**8. Risk.** `dim ker A_c` is huge — K3a measured mobility median 305 after 2-core on Delaunay patches
(F26/F29) — so the "branch count" is a question about a quadratic on a 300-dimensional projective
space, which is not a number, it is a variety. The idea may only be meaningful on the low-mobility
periodic patterns (Voronoi measured `m = 1–3`, U4), which shrinks it to the 2025 pattern library. Also,
`SparseQR` rank is unreliable above ~700 columns (F29) — kernels must be dense-computed, capping size.

**9. Effort.** Derivation **M**. Code **M**.

---

### Idea 7 — Orbit collapse: symmetry is why tilings are feasible and random graphs are not

**1. Title / type.** *The feasibility system has one inequality per symmetry orbit* —
**theorem**, with a quantitative corollary.

**2. Claim.** If `(M, σ, X)` is invariant under a group `G` and `t` is restricted to the `G`-invariant
subspace `𝕏^G`, then `q_e(t)` is constant on `G`-orbits of split edges and `μ_c(t)` on orbits of corner
incidences; the `0⁺` feasibility system therefore reduces to `#orbits` inequalities in
`dim 𝕏^G` unknowns. For a wallpaper-symmetric tiling `#orbits = O(1)` independent of patch size, so
feasibility is **size-independent** — an `n × n` patch of a working tiling works for every `n`. For a
random graph `#orbits = |E_split| + #corners`, growing linearly, and the measured per-edge inward
frequency `ρ ≈ 0.24–0.41` (K6 `inward_frac0`) makes the all-open event exponentially rare.

**3. Gap targeted.** The dichotomy 2026 never explains and `specs/ideator_round2.md` asks for
explicitly: "why authored tilings work and random graphs fail" (K1a: `0/8` tilings fail, `0/200`
random graphs succeed). 2026 §4.2's max-cut `σ` is symmetry-blind by construction.

**4. Novelty vs field.** (i) `field_rigidity` R12 (Fowler–Guest 2000) and R13 (Schulze–Guest–Fowler
2014, symmetric body-hinge) are the symmetry-adapted mobility machinery — they compute how *mobility*
decomposes by irreducible representation, not how a *collision-feasibility system* collapses. (ii)
IsoGami (F12) restricts to isohedral tilings and enumerates `3^m` joint assignments with filters — it
*uses* symmetry as a search restriction, never as an explanation for feasibility. (iii) R22/R2
(Mitschke 2013) search tessellation archives; no orbit-counting statement.

**5. Why it might be true / sketch.** `q_e` and `μ_c` are polynomial functions of `X` that are
equivariant for the action on cut structures (they are built from `d_e` and `Δu_e`, both natural). If
`γ ∈ G` maps split edge `e` to `e'` and fixes `X`, then `q_{e'}(X) = ± q_e(X)`, and the sign is `+`
because `G` acts by orientation-preserving symmetries of the cut (a reflection would swap `f, g` and,
per `zero_plus.hpp`'s convention note, flip both `d_e` and `dS_e`, leaving `q` unchanged). Restricting
`t` to `𝕏^G` keeps the equivariance. The counting corollary is then immediate.

**6. Kill experiment.** Driver: `kiri_analyze` + a 30-line orbit check. On the 8 authored tilings
compute `q_e` for all split edges at `X_ini` and count **distinct values to `10⁻¹²`**; on the 200 K1a
graphs do the same at `X₀`. **Dead if any authored tiling shows a number of distinct `q` values
comparable to `|E_split|`** (then symmetry does not collapse the system and the mechanism is wrong).
Second half: from the K1a sampler's cached per-sample `q` vectors, estimate `ρ` and check whether
`log P(all q_e > 0)` tracks `|E_split| · log(1 − ρ)` across graphs of different sizes. Under 5 minutes.

**7. If it survives, the demo.** The dichotomy as one plot: `#distinct q values` versus `|E_split|`,
tilings on the horizontal axis near `y = O(1)`, random graphs on the diagonal. Plus the size-independence
statement demonstrated on `n × n` patches of a working tiling for `n = 4 … 20`, certificate constant.

**8. Risk.** It is close to a tautology once stated ("symmetric things have few orbits") and a referee
will say so. The defence is that it makes a *falsifiable prediction* about the exponential rate, and
that the exponential rate is the only quantitative account of `0/200` anyone has. Second risk: the
authored tilings in `deployable_population()` are finite patches with boundary, where `G` is only a
subgroup of the wallpaper group and the orbit count is not `O(1)` but `O(perimeter)`.

**9. Effort.** Derivation **S**. Code **S**.

---

### Idea 8 — The vertex balance defect: a per-vertex non-deployability criterion

**1. Title / type.** *Split cuts break the barycentric balance; the defect predicts where deployment
dies* — **theorem / characterization**.

**2. Claim.** Define, at each interior vertex `v`, the **balance defect**
`δ_v := ‖ Σ_{u ∈ In(v)} (x_u − x_v) ‖` restricted to hinge in-edges. In the pure case
`δ_v = 0` for all `v` by (BAL), the `k_v` copies of `v` have `0⁺` separation velocities
`∝ J(x_u − x_v)` summing to zero, so for `k_v ≥ 3` they positively span the plane and **no copy of `v`
can enter its neighbour's material at `0⁺`**. With split cuts, `δ_v > 0` in general, and I claim the
corner margin obeys `min_c μ_c(v) ≤ −C(v) · δ_v + O(local geometry)` — i.e. the vertices where
deployment dies are exactly the high-defect ones. Corollary: a graph is `0⁺`-infeasible if some vertex
has `δ_v` exceeding an explicit local bound.

**3. Gap targeted.** F30 / K6 pass 3: "a design with every `q_e > 0` and every `a_f > 0` still has
`Θ_max = 0` on essentially every graph, and the exact scan reports the binding contact as vertex-edge"
(`zero_plus.hpp`). The corner obstruction is the binding one and neither paper has any account of it —
2026 Eq. (7)–(9) penalises split-edge separation only (F6), so it is blind to it *by construction*.

**4. Novelty vs field.** (i) 2025's `β = 2π − α_i − α_j` bound (F18) is a vertex-local *angle* bound
for the pure case; it says nothing when split edges are present, which is exactly the regime here. (ii)
`field_kirigami` R18 (Liu et al. 2024) give a closed-form collision test for Escher quad tessellations,
a curve-intersection test, not a criterion in the cut combinatorics. (iii) `notes/screen_bundle.md` S3:
"contact harmonic PARTIAL (Grima 2012 locking angle, hinge-adjacent only)" — the hinge-adjacent case is
precisely the pure case my theorem says is safe.

**5. Why it might be true / sketch.** The `0⁺` separation of two copies of `v` is
`sin(θ/2) (S_p − S_a)` with `S_p − S_a = J[2Δu − (σ_p − σ_f) x_v]` (`zero_plus.hpp` corner section). For
copies separated by a hinge in-edge `u → v`, this is `2 σ J(x_u − x_v)`: the copy moves **perpendicular
to the in-edge**, i.e. the fan boundary sweeps outward. The `k_v` such vectors are `J` applied to the
vectors of (BAL), so they sum to zero; a zero-sum set of `k ≥ 3` nonzero vectors is not contained in
any open half-plane, so the fans move apart in a genuine pinwheel. Deleting a term (a split edge at
`v`) leaves a set with nonzero sum `J·(residual)`, which *can* lie in a half-plane, and then some fan
moves into another. `δ_v` is the norm of that residual.

**6. Kill experiment.** Driver: extend `kill_k6.cpp` (it already computes `corner_incidences()` and
`zero_plus_corner_margin()`). On the 60-graph K6 subset, for both `σ_mc` and `σ_def`, compute `δ_v` per
interior vertex and correlate with `min_c μ_c` over the incidences at `v`. **Dead if the rank
correlation between `δ_v` and `−min_c μ_c(v)` is below 0.3, or if the argmin vertex of `μ` is not among
the top decile of `δ_v` on at least half the graphs.** Also check the positive half directly: on the 8
authored tilings and on split-free designs, verify `δ_v = 0` and `min μ_c > 0` at every vertex with
`k_v ≥ 3`, and find the `k_v = 2` exceptions (where the pinwheel degenerates to two antipodal vectors
and the argument gives nothing). Under 10 minutes on the cached K6 population.

**7. If it survives, the demo.** A heat map of `δ_v` over a Voronoi patch with the certified first
contact marked, on ≥ 20 graphs; and the `k_v = 2` degenerate case worked out as the sharp boundary of
the theorem (which is exactly the rotating-squares vertex, so the theorem has to be tight there).

**8. Risk.** The inequality direction is the weak point: I can prove the *pure* case cleanly (zero-sum
⇒ positive span ⇒ no ingress), but "defect large ⇒ ingress" needs a quantitative local bound with a
constant `C(v)` depending on the star's angles, and that constant may be so pessimistic that the
criterion never fires. Then it degrades to a correlation, not a theorem.

**9. Effort.** Derivation **M**. Code **S**.

---

### Idea 9 — Deployability-preserving graph surgery: Eulerian repair with a certificate

**1. Title / type.** *Any planar graph is `O(|T|)` local edits away from a certified-deployable one* —
**algorithm** (with a guarantee from Idea 1).

**2. Claim.** Given any planar `M`, let `T` be its odd-degree interior vertices. Inserting a minimum
T-join's worth of local edits — each edit an edge insertion or a vertex split that flips the parity of
exactly two vertices — produces `M*` with all interior degrees even, hence (Idea 1) admitting `σ*` with
`E_split = ∅`, hence `dim 𝕏 = 0`, `X = X₀`, and `Θ_max = min(min_v β_v, π) > 0` **certified**. The
cost is `|T|/2` edits and is computable by minimum-weight T-join (Edmonds), polynomial.

**3. Gap targeted.** Both papers take the graph as given and optimise `σ` and `X`. 2026 §4.2's `σ`
pipeline (F3) is the only design freedom exercised, and K1a/K5 show it is not enough. The move "change
the graph, not the orientation" is nowhere in either paper.

**4. Novelty vs field.** (i) `field_kirigami` R9 (Chen–Choi–Mahadevan 2020) vary *cut topology* at fixed
tile geometry, not the underlying planar graph, and have no guarantee. (ii) R17 (Dudte et al. 2023)
build patterns additively from four-bar linkages — a constructive family, not a repair of an arbitrary
input. (iii) IsoGami (F12) enumerates joints on isohedral tilings; nothing about non-tiling inputs.

**5. Why it might be true / sketch.** Idea 1 supplies the guarantee; the only content here is that the
parity repair can be done **planarly and locally**. Splitting a vertex `v` of odd degree into two
vertices joined by an edge changes the two new degrees to `a+1` and `b+1` with `a + b = deg(v)`; choose
`a` odd so both are even. Alternatively, insert a degree-2 vertex on an edge between two members of
`T`: that changes both endpoints' degrees by... nothing — so the correct primitive is inserting a chord
inside a face joining two odd vertices, which flips both parities and preserves planarity because both
lie on the same face. Minimum-weight T-join on the face-adjacency structure gives the optimal set.

**6. Kill experiment.** Driver: new `kill_k10.cpp` using `generators.cpp` for input and `mesh.cpp` for
the surgery. On 20 K1a graphs: compute `T`, run the greedy chord insertion, verify all interior degrees
even, 2-colour the interior dual, and run the full certificate. **Dead if fewer than 18/20 repaired
graphs certify `POS ∧ NOOVERLAP ∧ NOROOT` with `Θ_max > 0`.** The likely failure is not parity but
injectivity of `X₀` (Idea 4), which is a *different* obstruction: if `E_split = ∅` gives `dim 𝕏 = 0`
then `X₀` is forced and may self-intersect with no freedom to repair it. Under 15 minutes.

**7. If it survives, the demo.** ≥ 20 pairs (input graph, repaired graph) with certified `Θ_max`,
`0` → `min β`; a hero fabricable export of a repaired Voronoi patch; and the honest counter-table of
how many repairs died on `X₀` injectivity instead of on parity.

**8. Risk.** The obvious one, and it is serious: even-degree removes the split-cut collision source but
does **not** make `X₀` injective, and with `dim 𝕏 = 0` there is nothing left to optimise. The idea may
prove only that the two obstructions are *independent*, which is a weaker (but still real) result.

**9. Effort.** Derivation **S** (inherits Idea 1). Code **M**.

---

### Idea 10 — Certified deployment as a hard constraint in the 2025 inverse-design pipeline

**1. Title / type.** *Replacing simulation checks with the exact certificate* — **algorithm**.

**2. Claim.** The exact `Θ_max` of T4.2″ and the certificate `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT` (validated
187/187 against bisection, K2a) can replace the simulation-based deployability check in **both**
papers' pipelines, and the closed-form second root `θ* = 2 arctan(−r/p)` (T5.4) supplies a
differentiable margin `max(0, −p/r)` to optimise against — a drop-in replacement for the Eq. (9)
penalty that is blind to `p` by construction. On the 16 released patterns of
`baseline/kirigami_tessellations` against their 4 target meshes, this produces designs with *certified*
deployment where the released pipeline's check is a false negative or a false positive.

**3. Gap targeted.** 2026 §7 limitation (i) (range only via FK simulation) and F6/Eq. (9)'s
first-order-only penalty; K1c measured the authors' native `prevent_intersections` silently re-closing
on **9.47 %** of designs. 2025's pipeline checks deployment by simulation as well
(`specs/ideator_round2.md` §what-we-want (e)).

**4. Novelty vs field.** (i) IsoGami (F12) uses IPC contact in a continuation — numerical, per-trajectory.
(ii) `notes/screen_r1.md`: "certified range-optimal projection NOT FOUND (IsoGami / PyKirigami /
arXiv:2509.22002 all numerical)". (iii) `field_kirigami` R18 (Liu 2024) has a closed-form collision test
but for one specific curved-contour quad family.

**5. Why it might be true / sketch.** T5.3 gives the exact false-negative condition (T5.5): Eq. (9)
certifies `r > 0` and is blind to `p`; the second root exists in `(0, π)` iff `p < 0`, and
`τ* = −r/p` is a ratio of two quadratic forms in `X`, hence smooth and cheap. Substituting
`max(0, −p/r)` for `sign(r)` costs the same `O(|E_split|)` per evaluation and removes the entire
false-negative class, subject to the two extra clauses (interval test at `θ*`, wedge/entering test).

**6. Kill experiment.** Driver: `kiri_analyze` on the 16 released patterns as deployed by
`baseline/kirigami_tessellations`. **Dead if none of the 16 shows a certificate violation and none
shows an Eq. (9) false negative** — i.e. if the released library is already certificate-clean, there is
no headroom and the contribution is a re-implementation. Under 20 minutes, given baseline parity is
already established (8/8, `≤ 4e−8`).

**7. If it survives, the demo.** The 16×4 table (pattern × target) of released-pipeline verdict vs
certified verdict, with every disagreement rendered; plus the substituted penalty run head-to-head with
the authors' native `prevent` on ≥ 20 designs, refereed by our bisection.

**8. Risk.** K2b already killed the *range-margin* claim (median gain 0 against native `prevent`). This
idea has to be scoped strictly to **certification**, not to beating them on range, or it dies the same
death. It is also the least theorem-like item on my list and I rank it accordingly.

**9. Effort.** Derivation **S** (all of it is already in T4/T5). Code **M/L** (2025 pipeline
integration).

---

## 2. Self-attack — the strongest rejection I can write against my top three

### Against Idea 8 (vertex balance defect)

*Reject.* The authors present as a "theorem" a statement whose proved half is vacuous and whose useful
half is a correlation. The proved half says: in the pure-hinge case, where the paper's own measurements
already establish that `Θ_max = min(min β, π) > 0` unconditionally (their F18), the copies of a vertex
do not collide. That is not news; it is a restatement of a fact they list as *given*. The half that
would matter — that a large balance defect *forces* an ingress — is stated with an unnamed constant
"`C(v)` depending on the local geometry", and the experimental section then reports a rank correlation,
which is what one reports when one does not have a theorem. Worse, the mechanism is presented as an
explanation of the paper's central negative result, but the paper's own K6 data show the binding contact
is vertex-into-edge on graphs where *every* `q_e > 0` — that is, in a regime where the split edges have
already been repaired, so "split cuts break the balance" cannot be the operative mechanism there. The
`k_v = 2` case, which is the single most important vertex type in the entire kirigami literature
(rotating squares), is explicitly excluded by the hypothesis `k_v ≥ 3`. A theorem about kirigami that
does not cover rotating squares is a theorem about nothing.

### Against Idea 2 (Farkas certificate of emptiness)

*Reject.* The S-procedure is 60 years old and the observation that the constraints are quadratic is two
lines of the authors' own T1. What is offered as the contribution is therefore entirely the *empirical*
claim that a small-support certificate exists — and that claim is unproved, unmotivated by anything
except a hope about vertex stars, and, on the authors' own account, may require an SDP in 400 variables
that they do not run. If the certificate is found only by search, then the paper reports "we searched
and found a matrix", which is not a characterization of anything; a reader cannot look at `(G, σ)` and
tell whether a certificate exists. Meanwhile the corner constraints `μ_c`, which the authors themselves
identify as the binding ones, are only piecewise quadratic, so the aggregation is not even well posed
over the full constraint set; the certificate is built from `q` alone, which K6 pass 3 already showed is
the non-binding half. So the paper proves emptiness of a *relaxation* that the authors' own data say is
non-empty. Finally: certifying that a family of designs *does not work* is a negative result about a
population the authors chose (random Delaunay/Voronoi). No practitioner cuts a Voronoi diagram.

### Against Idea 1 (odd-vertex deficiency)

*Reject.* The equivalence "some `σ` has `E_split = ∅` iff every interior vertex has even degree" is the
planar-duality fact that `G*` is bipartite iff `G` is Eulerian, taught in a first course, composed with
the paper's own definition of a split edge. The prior work the authors themselves cite — Segall et al.
2025 — *assumes* even valency, which is to say the community already knew that even valency is the
regime where the construction is clean; the present paper supplies the converse direction for a
condition nobody wanted to weaken. The Heawood 3-colourability remark is decoration. That leaves
`|E_split| ≥ |T|/2`, a T-join bound whose proof is the handshake lemma, and whose consequence — that
graphs with many odd vertices have many split cuts — is neither surprising nor, by itself, an obstruction:
the paper must still import Ideas 2 and 7 to convert "many split cuts" into "no deployable embedding",
and those are the ideas doing the work. Strip the parity theorem out and the paper is unchanged. There
is also a hypothesis problem: the cycle-space argument is stated for disk topology, and the paper's
headline population (Voronoi patches) is a disk only by construction, while the interesting periodic
case has a torus, where the argument as written fails.

---

## 3. Final ranked list

1. **Idea 8 — vertex balance defect.** The only idea on my list that *explains* the measured mechanism
   rather than restating it: (BAL) is a zero-sum condition on the `0⁺` separation velocities, so the
   barycentric row is what prevents collapse and split cuts are the deleted terms. Proves the pure case
   outright and predicts where the corner contact fires.
2. **Idea 2 — Farkas certificate of emptiness.** Converts K1a/K5/K6's `0/200`, `0/400`, `0/31` from a
   measurement into a proof, and is the only idea that answers 2026 §7 limitation (i) in the negative
   direction. Highest payoff, real chance of failing.
3. **Idea 1 — odd-vertex deficiency.** Cheapest to prove, cheapest to kill, supplies the *positive*
   half of the dichotomy (a graph class with a guarantee) and the lower bound that Ideas 2 and 7 need
   as input. Weakest on novelty; strongest on certainty.
4. **Idea 3 — area identity and expansion floor.** An exact closed form nobody has, an immediate
   one-number non-deployability test, and the only route I found to 2025's "realizable shape space"
   question. Ranked below the top three only because the identity itself is easy.
5. **Idea 4 — index criterion for Eq. (6) injectivity.** Directly closes F17 (the baseline's real
   practical failure) with the right machinery (`field_tutte` row 4), and the sliding-out-neighbour
   counterexample is a genuinely non-guessable outcome. Risk: may deliver only a necessary condition.
6. **Idea 7 — orbit collapse.** The cleanest available account of the tilings-vs-random dichotomy and
   testable in five minutes, but close to tautological once stated.
7. **Idea 5 — locality of the split-gap sign.** Genuinely two-sided: if locality holds we get a new
   `σ` rule; if it fails we get a theorem saying no local `σ` rule can exist, which explains K5.
8. **Idea 6 — tangent cone at the flat state.** The right way to attack 2026 §7 limitation (ii), but
   the measured mobility (median 305) means the branch count is a variety, not a number, so the result
   may only exist for the low-mobility periodic library.
9. **Idea 9 — Eulerian repair.** Attractive as a "we fix your graph" story, but likely to be blocked by
   `X₀` injectivity, which `dim 𝕏 = 0` leaves no freedom to repair.
10. **Idea 10 — certified 2025 pipeline.** Least theorem-like, and adjacent to the already-dead K2b
    range-margin claim. Kept because the certificate exists and the integration cost is the only
    obstacle.

---

## 4. What I could not verify

* I did **not** re-read the 2026 or 2025 paper text this session. Every paper citation above is taken
  from `STATE.md` F1–F30 and the reader notes; the two places where a precise quotation would matter are
  Idea 1's claim that 2025 assumes even valency (taken from `notes/field_kirigami.md` row R21, not from
  the paper) and Idea 3's claim about 2025 Fig. F.4 (taken from `specs/ideator_round2.md` §12).
* The area identity of Idea 3 is derived here but **not** numerically checked; its kill test is the
  check.
* The closed form `Δu = σ_g (x_{src(h→g)} − x_{src(f→h)})` of Idea 5 is derived only for the case where
  `f` and `g` have a common `Γ`-neighbour; I have not established how often that holds, which is why
  the measurement is the kill test rather than a footnote.
* Idea 7's claim that a reflection leaves `q_e` invariant rests on reading the sign-convention paragraph
  of `code/src/method/zero_plus.hpp`, not on a proof.
* Sign conventions: everything above is in the code convention `a_f = −σ_f θ/2` (`derivations/core.md`
  0.7). I have not re-derived T5.3's `(p, q, r)` independently; I use them as given, and Idea 3's
  identity inherits any sign error there. The identity's kill test would catch it.
