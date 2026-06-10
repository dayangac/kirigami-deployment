# Round 2 — adversary persona: the rejection, then ten ideas that survive it

Persona: a hostile, well-read TOG reviewer who has refereed four "we generalize X to arbitrary
graphs" submissions this cycle and rejected three of them.

Inputs actually read in this run: `specs/common_preamble.md`, `specs/ideator.md`,
`specs/ideator_round2.md`, `STATE.md` (all of it: F1-F30, D5-D9, U1-U9, dead ends, session log),
`ideas/ranking.md` Sec. 0 and Sec. 3 (the five round-1 hostile paragraphs, in full),
`results/kill/KILL_REPORT.md` (corrections list, certificate table, verdict table),
`derivations/check.md` Sec. 0 verdict summary and the per-step tables for T1-T7,
`notes/repo_2025.md` (section skeleton; body still being written),
`notes/field_kirigami.md` / `notes/field_rigidity.md` / `notes/field_tutte.md` (row tables and the
rows named below), `notes/screen_bundle.md` + `notes/screen_r1.md` (verdict lines),
`code/src/method/*.hpp` in full (`zero_plus.hpp`, `contact.hpp`, `deploy_basis.hpp`,
`mobility.hpp`, `periodic_jacobian.hpp`), `code/apps/` listing,
`results/kill/k7/*.csv` headers and the `c3_*.log` tails.

**Honesty statement.** I ran no code and no web searches in this session. Every number below is
quoted from `STATE.md` or `results/kill/KILL_REPORT.md`. Three citations I use for positioning
(Connelly-Demaine-Rote expansive motions; Streinu pointed pseudo-triangulations; Maxwell-Cremona
reciprocal diagrams) are **not** in the field tables and I did **not** verify them in this run:
they are flagged `[SCREEN-ME]` wherever used and must be screened before any of them appears in a
paper. Nothing else is asserted from memory.

---

# Part 1 — The rejection of the current plan (D7 + D8)

*Manuscript under review: "Exact deployment calculus and a validity certificate for uniformly
deployable kirigami on arbitrary planar graphs", contributing (i) closed-form `Theta_max` and a
certificate `POS AND NOOVERLAP(eps/2) AND NOROOT`, (ii) an emptiness result on random planar
graphs, (iii) a 9.47% re-closure rate for the prior work's Eq. (9), (iv) rank/hole-count
corrections, and (v) a designable periodic Poisson family.*

**Recommendation: reject.** Five reasons, in decreasing order of how much they hurt.

**R-1. The emptiness result is a statement about two heuristics on one random generator, and the
submission knows it.** The claim on offer is that the "full design space" of Segall et al. 2026 is
generically empty. What was measured is: with `sigma` from their Eq. (1) max-cut *or* from the
submission's own defect-minimizing local search, and with `X0` from their Eq. (6) least-norm
projection, plus Gaussian samples and trust-region probes in the null space, 0 of 400 designs on
Voronoi and Delaunay patches had `Theta_max > 0` (F25, F30). Zero of four hundred *samples* from a
space of dimension up to 1017 (F19) is not a theorem; it is a Monte-Carlo estimate of the measure of
a set whose measure nobody claims is large. Emptiness of a semialgebraic set is a *decidable*
question, the submission has assembled exactly the polynomial description needed to decide it
(T5.2a-d), and then decides it by sampling. Meanwhile the submission's own control says authored
tilings pass 8 out of 8. So the honest reading of the data is "our orientation heuristic and our
projection do not work on our random graphs", and the reviewer cannot distinguish that from "the
design space is empty" on the evidence given. The confound is not exotic: random Voronoi and
Delaunay patches have generic edge lengths and generic vertex figures, and their split-forest
components are long relative to face size, whereas every pattern in the original paper's figures has
equal edge lengths and symmetric vertex figures. The submission never varies that ratio. One
interpolation experiment - jitter a tiling continuously toward a random patch - separates the two
explanations, and it is not in the paper.

**R-2. The mathematical core is the Weierstrass substitution, and the one theorem that was not
folklore has been withdrawn.** `Y_theta = cos(theta/2) C(X) + sin(theta/2) S(X)` says that a rigid
rotation through `+-theta/2` is trigonometric in `theta`. That faces counter-rotate by `+-theta/2`
on a 2-colourable rigid-unit network is Grima-Evans rotating-rigid-unit folklore (field_kirigami R1)
and, for arbitrary bipartite polygon networks, is Acuna et al. 2022 (F21). Once you have that, every
polynomial predicate along the path is a first harmonic `p + q cos theta + r sin theta` because
determinants of trigonometric linear forms are, and the roots come out of the tangent half-angle
substitution, which is nineteenth-century. `T5.2a`'s "the usable region is semialgebraic" is
Tarski-Seidenberg applied to a formula the submission wrote down; the submission's own Checker says
so in as many words. The single claim that would have been a theorem - `H-LOC`, an `O(n)` certified
active set - is listed in `check.md` as **refuted**, with the swept-disc drift growing from 4.6 to
20 with patch diameter. What remains after the withdrawal is: enumerate a candidate list of quadratic
harmonics, deflate the ones with a root at zero, sort the roots, and take the first gap. That is a
correct and careful *implementation*, and the graze handling is a genuine catch (`hexagons_auto`
1.047 vs 2.094), but it is a bug in a naive implementation, not a characterization. And by the
submission's own measurement the whole apparatus is **inert** on split-free patterns, where
`Theta_max = min(min beta_e, pi)` exactly - which is Segall et al. 2025 Sec. 4.2's published local
formula with a cap.

**R-3. The negative results are errata about an under-specified baseline, and the headline number
is convention-dependent.** Eq. (9) is printed in the original with an unspecified barrier `B`, an
unspecified weight `gamma`, and an unspecified diameter-selection rule; the submission's own Reader
records all three as absent from the paper. The submission then reports that this optimizer
"silently re-closes on 9.47% of designs". But `KILL_REPORT.md` states that the rule as literally
written gives **65.09%**, and the same rule applied at its own reference range gives **exactly 0**,
because it is self-contradictory. A quantity that takes the values 0%, 9.47% and 65% depending on
which of three defensible conventions the authors adopt is not a measurement of the prior work; it
is a measurement of the submission's convention. Publishing the middle number as the headline, with
the other two in a footnote, is not acceptable. Similarly for the corrections: `rank(L) = H - 1` on
a boundary-free pattern is the observation that a boundary-free incidence matrix has the all-ones
vector in its left kernel - Kirchhoff, 1847 - and the submission concedes it is "harmless in
practice" because the boundary rows pin the translation. `H = |E_hinge| - |F| + c(Gamma)` is Euler's
formula on the hinge graph. That the authors' code adds a constraint row for boundary-touching
split components is a two-line over-constraint in someone else's repository. These belong in a
technical report or an author correction, not in the contributions list of a TOG paper.

**R-4. The constructive half is the competitor's paper with the interesting part removed.** IsoGami
(F12, field_kirigami R22) computes mobility, runs contact-aware continuation with IPC, and designs
periodic isohedral tilings, concurrently and with fabricated artefacts. Against that, "the
achievable set of periodic Jacobians is an affine subspace" is the statement that a linear map has a
linear image; `J(theta) = cos(theta/2) I + sin(theta/2) K` is the same half-angle identity as R-2
transported to the unit cell; and the resulting `nu(theta)` is Grima-Evans with the unit cell left
symbolic (field_kirigami R1). Worse, the submission's own interim K7 data show the *design* half at
6 of 22 on target-driven design, and the `c3_*.log` lines show the second target failing the
certificate outright (`cert=0`, `thmax=0`) on both `snub_square_2x2` and `voronoi_torus_3_n45`. So
the theorem is trivial and the algorithm does not work yet.

**R-5. Novelty positioning against the nearest neighbours is missing where it matters most.** Segall
et al. 2025 Sec. 4.2 already gives closed-form `theta_max` for the hinge-adjacent case; Liu et al.
2024 (field_kirigami R18) already give closed-form non-adjacent conditions for quad tessellations;
Dang et al. 2021 (F20, R12) already prove an iff for rigid deployability of quad kirigami. The
submission's delta is the non-adjacent complement on non-quad graphs. That is a completion. It may
be worth a section. It is not worth an abstract.

**What would change my mind.** A result whose answer I cannot guess before reading the experiment,
about the design space itself rather than about the authors' two heuristics, and stated so that it
survives being told "your random generator is doing the work".

---

# Part 2 — Ten ideas that survive that reviewer

Design rules I imposed on myself, from the rejection above:

- **No idea may depend on the random generator for its statement.** Every claim is either a
  per-instance *certificate* (a proof object checkable on the graph in front of you) or a *predictive
  scalar* validated across generators including the authored tilings.
- **No idea may be a one-line corollary of the half-angle identity.** The half-angle identity is
  infrastructure; it appears in the *proofs*, never in the *claims*.
- **Every idea must have an outcome I cannot guess.** If I can predict the number before the run,
  the idea is in Part 3.

## The enabling identity used by ideas 1, 2, 3, 5 and 9 (already verified, not claimed as new)

From `derivations/check.md` (all AGREE, with the numerical residuals quoted there):

- T1.3: across a hinge edge `f -> g` with hinge vertex `v_src`, the auxiliary face variables satisfy
  `u_g - u_f = sigma_g * x_{src}` , **independent of `theta`**.
- T1.5: `y_{(v,f)}(theta) = cos(theta/2) x_v + sin(theta/2) J (2 u_f - sigma_f x_v)`.

Read together, `u : F -> R^2` is a **discrete potential on the hinge graph** whose prescribed
increments are the hinge-vertex positions signed by `sigma`, and **Eq. (2) (Prop. 4.1) is exactly the
integrability condition for that potential**: the increments integrate iff their circulation
vanishes on every cycle of `Gamma`, which is F4/F11's hole-preimage sum. Consequences that the
project has not yet used:

1. For a **split** edge `e = {a,b}` between faces `f, g` with `sigma_f = sigma_g`, T1.B gives
   `dS_e = S_{(v,g)} - S_{(v,f)} = 2 J (u_g - u_f)` (the `sigma_f x_v` terms cancel), so with
   `d_e = x_b - x_a`,

   ```
   q_e  =  det(dS_e, d_e)  =  2 det(J Delta u_e, d_e)  =  -+ 2 <Delta u_e , d_e> ,
                                             Delta u_e := u_g - u_f
   ```

   the global sign being fixed once and for all by the stored-CCW convention of `zero_plus.hpp`.
   This is consistent with `check.md` T5.3 (`r = det(d, J Delta u) = <d, Delta u>`, verified to
   1.9e-14 over 4300 split-edge samples) and with `zero_plus.hpp`'s "`q_e` is a positive multiple of
   `r`". **The `0+` obstruction is an inner product between a potential jump and an edge vector.**
2. For a **corner** incidence (copy of `v` in face `h` entering the corner of face `f` at `v`), the
   same substitution gives `dS = 2 J (u_f - u_h) - (sigma_f - sigma_h) J x_v`, so the corner margins
   `mu` of `zero_plus.hpp` are also bilinear in `(u, X)`.

So the entire `0+` collision theory - the mechanism behind F30, the emptiness result - is a system of
sign conditions on `<Delta u, d>`. That reformulation is what makes ideas 1-3 and 9 possible; I
claim it as *machinery*, and I say in every idea below which part is new.

---

## Idea 1 — The expansive cone: non-uniform deployability is an LP, uniform is one ray in it

**Type:** theorem + algorithm.

**Claim.** Fix a valid flat embedding `X` and its cut structure. The set of infinitesimal flexes of
the body-and-pin framework that separate every cut at first order,

```
P(X) = { V in ker A(X, 0) :  det( V_g(x_v) - V_f(x_v), d_e ) > 0   for every split edge e,
                             and every corner copy leaves its host corner cone },
```

is an **open cone that is polyhedral on each corner-branch chart**, `P(X)` is decided by one linear
program per chart, and the uniform deployment is the single ray `R_{>0} * sigma` inside
`ker A(X,0)`. Consequently `Theta_max(X) > 0` at first order **iff** `sigma in P(X)`, and the
measured emptiness F25/F30 is the statement that *one specific ray of a cone of dimension up to 305
(F26) points the wrong way* - not that the cone is empty.

**Gap targeted.** 2026 Sec. 7 limitation (ii), verbatim in F8: non-uniform hinge-angle deployment and
the DOF of the deployment space, which the authors say "could be achieved by sensitivity analysis"
and do not do. Also 2026 Sec. 6/7 limitation (i): they characterize the collision-free range only by
FK simulation.

**Novelty vs field.** (a) field_rigidity **R17** Borcea-Streinu geometric auxetics defines auxetic
paths by a PSD condition on the *lattice Gram matrix* of a periodic framework; the cone here is over
*face copies of a cut structure on a finite patch* and its inequalities are non-penetration of
specific duplicate vertices, not lattice expansion - different variables, different inequalities,
and theirs has no cuts. (b) field_rigidity **R4/R5** Tay-Whiteley and Jackson-Jordan give
combinatorial rank of body-and-hinge frameworks in *generic* position; F13 records that generic
criteria predict rigidity for every working pattern, so they cannot see this cone at all. (c)
`[SCREEN-ME]` Connelly-Demaine-Rote expansive motions and Streinu's pointed pseudo-triangulations
concern *expansive* motions of polygonal linkages (all pairwise distances non-decreasing), a
strictly stronger and differently-indexed condition than "each cut opens"; they must be screened
before this idea is written up, and if a screen finds the exact inequality system there, the idea
degrades to idea 2 alone.

**Why it might be true / sketch.** The flex variables are one angular velocity `omega_f` and one
translation `w_f` per face; `A(X,0)` is `mobility.hpp::build_A` and `ker A` has dimension `m`
(median 305 after 2-core, F26/F29). Every velocity `V_f(y) = omega_f J y + w_f` is **linear** in the
flex. The relative velocity of two copies of the same source vertex is therefore linear in the flex,
and, `X` being **fixed**, `d_e` and the corner edge directions `e1, e2` are **constants**. So each
split-edge condition `det(V_g - V_f, d_e) > 0` is a *strict linear inequality in the flex*. Reflex
corners give `mu = min(-g1, -g2) > 0`, a conjunction of two linear inequalities - convex. Convex
corners give `mu = max(-g1, -g2) > 0`, a **disjunction** - so `P(X)` is a union of polyhedra indexed
by a branch choice per convex incidence. This is the honest statement and it is what makes the idea
falsifiable in both directions: **feasibility of any single branch LP is a sound existence proof**
(the LP witness is a flex that separates everything), whereas **emptiness requires all branches**,
which is exponential in general. The escape is that the branch is not free: at a convex corner only
the two half-planes adjacent to the actual incoming direction are geometrically reachable, so a
default branch assignment (the side the copy is already on at `X`) plus one repair pass covers the
cases that matter, and any success is a proof regardless.

The quadratic-to-linear step is the whole trick, and it is worth stating why the project has not
seen it: `zero_plus.hpp` searches over the **embedding** `X` with the flex slaved to it
(`V = S(X)/2`, hence `q_e` quadratic in `t`, hence non-convex). Freeing the flex from the embedding
turns the same geometry into an LP at the cost of leaving the uniform branch. That is exactly the
trade the 2026 paper refuses to make and the 2025 paper makes only pictorially (Fig. 10).

**Kill experiment.** Driver `kill_k8a.cpp`, reusing `mobility.hpp::build_A`,
`zero_plus.hpp::split_copies` / `corner_incidences`, and the K1a population generator in
`kill_common.hpp`. For each of the **100** K1a graphs (101-793 faces) at `X = X_ini` (which is
injective 200/200 by F25, so the LP data are valid): build `A`, get a kernel basis `N` by dense
`ColPivHouseholderQR` on the 2-core (F29 warns sparse rank is unreliable above ~700 columns - use
dense on the 2-core only), assemble the linear inequality rows, and solve
`max_z min_i l_i(N z)` s.t. `||z||_inf <= 1` by the multiplicative-weights dual
`min_{lambda in simplex} || N^T l(lambda) ||_1` (about 60 lines, no LP library; the dual value being
`<= 0` is a **Farkas certificate** of infeasibility for that branch, and a primal `z` with positive
margin is an existence proof).
**Kills the interesting claim if:** the number of graphs with a strictly positive certified margin is
**below 20 of 100** - i.e. non-uniform flexes do *not* rescue the graphs the uniform ray fails on.
**Kills the whole idea if:** on the 8 authored tilings (which have `Theta_max > 0`) the LP reports
`sigma` outside `P(X)`, since that contradicts the equivalence. Runtime: `m ~ 300`, `~2000` rows,
about 1 s per graph - **< 10 min for the whole population**.

**If it survives, the demo.** Two figures. (i) A scatter over 100 random graphs of certified LP
margin against uniform `Theta_max` (which is 0 on all of them), with the 8 tilings overlaid in the
upper-right quadrant - one picture that says "the space is not empty, the ray is wrong". (ii) The
hero: a Voronoi patch with `Theta_max = 0` under Eq. (6), opened to a visible non-uniform state by
integrating the LP flex (idea 8), exported through `kiri_export` as a laser SVG so the fabricated
piece deploys where the published pipeline says it cannot.

**Risk.** The convex-corner disjunction defeats the LP in practice: every branch assignment I can
choose cheaply is infeasible while some exotic one is feasible, so I get neither an existence proof
nor an emptiness proof and the result is "inconclusive on 80 of 100". Second risk: the screen finds
the identical cone in the expansive-motion literature.

**Effort.** Derivation S (the algebra is T1.3 + T1.5, already checked). Code M.

---

## Idea 2 — A per-graph emptiness certificate for the usable set, by convex duality

**Type:** theorem + algorithm.

**Claim.** For the `0+` feasible set of the *uniform* branch,
`U_0 = { t in R^{2k} : q_e(t) > 0 for all split e, mu_i(t) > 0 for all corner incidences i,
a_f(t) > 0 for all f }`, with every constraint an explicit **quadratic** form in the shape
coordinates `t` (`zero_plus.hpp` states this and `check.md` T3.2/T5-a verify the exact quadratic
dependence to 3.4e-14), the following is a *proof* of `U_0 = {}` on a given graph:

```
there exist weights  lambda >= 0, sum lambda = 1,  with  Q(lambda) := sum_i lambda_i M_i   psd-negative,
i.e.  lambda_max( Q(lambda) ) <= 0   as a form on the homogenised variable w = (1, t).
```

Moreover `lambda -> lambda_max(Q(lambda))` is **convex** on the simplex, so the best certificate is
computable by a subgradient / multiplicative-weights loop whose *success is self-verifying* (one
eigendecomposition) and whose failure gives, by duality, a near-feasible direction.

**Gap targeted.** This is aimed at reviewer point R-1 above, which is a gap in *our own* plan rather
than in the Segall papers; the Segall-paper gap it closes is 2026 Sec. 6/7 limitation (i) - "developing
geometric characteristics for favorable deployment behavior directly from the embedding remains an
open problem" - by supplying the negative half of exactly that characteristic.

**Novelty vs field.** (a) field_tutte rows on Tutte-embedding variants (Groiss-Juttler-Mokris 2021,
screened in `screen_r1.md` as requiring positive weights on the whole star) give injectivity
certificates for the *embedding*, not emptiness certificates for a *deployment* constraint set.
(b) field_kirigami **R22** IsoGami filters candidate joint assignments numerically and reports
failures as failures; it has no infeasibility proof. (c) `[SCREEN-ME]` degree-2 Positivstellensatz /
S-procedure duality is textbook convex optimization; the novelty is not the tool but that the object
being certified - a kirigami design space - has never been given one, and that the certificate is
*sparse and readable*: `lambda` is supported on a handful of constraints and names **which** split
edges and corners are jointly responsible.

**Why it might be true / sketch.** Each `q_e(t) = det(dS_e(t), d_e(t))` is a determinant of two
*affine* functions of `t`, so its homogenised matrix `M_e` has **rank at most 4** and a completely
explicit factorisation `M_e = (1/2)(G_e^T E D_e + D_e^T E^T G_e)` with `E` the 2x2 rotation by 90
degrees and `G_e, D_e` the `2 x (m+1)` blocks already stored in `ZeroPlusForm::GS` and
`ZeroPlusForm::DD`. So `Q(lambda)` is a sum of explicitly low-rank matrices, computable without
forming anything dense in `m` beyond the final `(m+1) x (m+1)` matrix. If `Q(lambda) <= 0` then for
every `t`, `sum lambda_i (constraint_i)(t) <= 0`, so some constraint is non-positive and `U_0` is
empty; the conclusion is a genuine implication with no sampling in it. Existence of such a `lambda`
whenever `U_0` is empty is *not* guaranteed (that is the gap between the S-procedure and exact
Positivstellensatz for more than one quadratic) - which is precisely why this is an experiment and
not a lemma.

**Kill experiment.** Driver `kill_k8b.cpp`. Population: the **50** smallest K1a graphs (so
`m = 2 * dim_null` stays under about 400 and the eigendecomposition is milliseconds). For each:
build `ZeroPlusForm` (exists), assemble `M_e` for split edges, `M_i` for corner incidences (via
`corner_incidences`, taking the reflex/convex branch that is active at `t = 0`), `M_f` for face
areas; run 2000 multiplicative-weights iterations on `lambda`; report `min_lambda lambda_max(Q)`.
**Kill rule:** a certificate (`lambda_max <= -1e-8` after scaling) is found on **fewer than 25 of
50** graphs. Below that the certificate is too weak to replace the sampling argument and the
emptiness claim stays a Monte-Carlo statement. Also a hard soundness check: run the same procedure on
the **8 authored tilings**, where `U_0` is non-empty; **any** certificate found there is a bug and
kills the driver, not the idea. Runtime: dominated by 2000 eigendecompositions of a 400x400 matrix
per graph, about 20 s each - **under 20 min for 50 graphs**, and trivially shardable.

**If it survives, the demo.** A table replacing "0 of 200 samples" with "**proved empty** on N of 50
graphs, with a certificate of median support size k naming the responsible split edges", plus a
figure drawing the certificate's support on the mesh: the reader sees the obstruction as a subgraph,
not as a p-value. This is the single change that converts the project's negative result from a
statistic into a theorem-with-witnesses.

**Risk.** The S-procedure is lossy with many quadratics and `lambda_max` never goes negative, so the
answer is "no certificate found" everywhere - which proves nothing and reads, in a paper, exactly
like a failed idea. Second risk: numerical scaling - the constraints have wildly different units
(areas vs. inner products) and the simplex weights degenerate; mitigate by normalising each `M_i` to
unit Frobenius norm and reporting that the normalisation is part of the certificate.

**Effort.** Derivation S (the low-rank factorisation is two lines). Code M.

---

## Idea 3 — The orientation is the design variable: choosing sigma by potential geometry

**Type:** algorithm + characterization.

**Claim.** Using the potential identity, the `0+` sign of a split edge is
`sign <Delta u_e, d_e>` where `Delta u_e = u_g - u_f` is a **path sum of `sigma`-signed hinge-vertex
positions along `Gamma`**. Therefore `sigma` does not merely control *how many* split edges exist
(the 2026 Sec. 4.2 max-cut heuristic) - it controls, through the potential, *which side each one
opens on*. The claim: an orientation objective that maximises the number of split edges with
`<Delta u_e(sigma), d_e> ` of the correct sign at `X_ini` achieves a strictly positive certified
validity rate on the K1a population where max-cut `sigma` and defect-minimizing `sigma` both achieve
**0 of 200** (F25, F30).

**Gap targeted.** 2026 Sec. 4.2, quoted through `ideas/ranking.md`: "we empirically favor many
balanced, small holes rather than a few large ones. This motivates choosing face orientations that
minimize the number of split edges", for which the Reader records the supporting evidence as
"Nothing. No definition of 'balanced', no experiment."

**Novelty vs field.** (a) 2026 Eq. (1) is a max-cut relaxation on the dual graph, a *purely
combinatorial* objective; it cannot see `X` at all, and this one is purely geometric given the
combinatorics. (b) `ideas/ranking.md` R5 (the round-1 "defect beats max-cut" idea) minimises the
*deployability defect* at `X_ini`, a different and, as K5 measured, insufficient objective - it fixed
the flat sheet (defect down 122x, POS 58/200 -> 186/200) and still gave 0 of 200 range because it
*doubled* the split cuts (126 -> 321) and never looked at their sign. This idea targets the sign
directly and predicts *fewer* but *correctly-oriented* cuts. (c) field_kirigami **R9** Chen-Choi-
Mahadevan 2020 searches cut *topology* at fixed geometry on square tiles; no potential, no sign
condition, quads only.

**Why it might be true / sketch.** Flipping `sigma` on a single face `f` does three things at once:
it re-classifies each edge of `f` between hinge and split, it reverses the increment `sigma_g x_src`
on every hinge edge at `f`, and hence it changes `u` on an entire *side* of `Gamma`. The first effect
is what max-cut optimizes; the second and third are invisible to max-cut and are what determine the
`0+` signs. Concretely, since `du = sigma * x_src` along hinge edges, `Delta u_e` for a split edge is
the signed sum of hinge-vertex positions around any `Gamma`-path from `f` to `g`, so
`<Delta u_e, d_e> = sum_{hinge h on the path} +- <x_{src(h)}, d_e>` - a **linear function of the
orientation-induced signs** along that path. That makes a local-search or a linear relaxation
natural, and it makes the objective computable in closed form without ever forming `X0` or the null
space.

**Kill experiment.** Driver `kill_k8c.cpp`. Population: the same 200 K1a graphs. Compute
`sigma_mc` (existing, `orientation.hpp`), then run a face-flip local search on
`N_bad(sigma) = #{ split e : <Delta u_e, d_e> <= 0 at X_ini }` with `u` recomputed incrementally.
Then for the resulting `sigma_sign`: build the cut, solve Eq. (6) for `X0`, and evaluate
`validity_certificate` at `eps = 0.006` (the K2a value). **Kill rule:** certified valid on **0 of
200** - i.e. the same result as both existing orientation rules - or the local search fails to reduce
`N_bad` below 50% of its initial value on the median graph. Either outcome kills it. Runtime:
the local search is combinatorial and cheap; the 200 Eq. (6) solves dominate at a few seconds each
(F19: worst dense solve 3.4 s at ~3k faces, and this population is smaller) - **about 20 min**.

**If it survives, the demo.** The Pareto plot the round-1 Critic asked for and nobody has run:
`|E_split|` against certified `Theta_max` across three orientation rules on 200 graphs, showing that
the paper's stated principle ("minimize split edges") is optimizing the wrong axis. Hero: the same
graph cut three ways, two of which have `Theta_max = 0`.

**Risk.** The signs are not independently controllable - flipping one face fixes one edge and breaks
two - so `N_bad` plateaus at a constant fraction and the certified rate stays 0. Given K5's history
(a 122x improvement in the defect that bought exactly zero range), this is the most likely outcome
of the ten, and I rank the idea accordingly.

**Effort.** Derivation S. Code S.

---

## Idea 4 — The deployability transition: which scalar predicts that a graph is not deployable

**Type:** characterization.

**Claim.** There is a **generator-independent** scalar computable from `(G, sigma, X_ini)` alone -
before any null space is formed - that predicts certified deployability, and the transition from the
authored-tiling regime (8/8 deployable) to the random regime (0/200) is a **sharp function of one
geometric ratio**, not of "randomness". The candidate: `rho = median over split-forest components of
(component diameter) / (median face inradius)`, with the transition near a critical `rho*`.

**Gap targeted.** Directly the hole in *our* emptiness claim (reviewer R-1), and 2026 Sec. 6/7
limitation (i): "developing geometric characteristics for favorable deployment behavior directly from
the embedding remains an open problem". This is that geometric characteristic, on the negative side.

**Novelty vs field.** (a) field_kirigami **R3** Shan et al. 2015 and **R6** Tang-Yin study specific
perforated patterns experimentally, with no predictor across pattern families. (b) **R22** IsoGami
enumerates and filters `3^m` joint assignments per tiling; a filter is not a predictor and it does
not transfer between tilings. (c) 2026 Sec. 4.2's `|E_split|` is the paper's only proxy and F30 shows
it is not the operative one, since defect-minimizing `sigma` doubled `|E_split|` and changed nothing.

**Why it might be true / sketch.** By the potential identity, `Delta u_e` for a split edge is a path
sum of hinge-vertex positions over a `Gamma`-path whose length grows with the split-forest component
diameter, while `d_e` is a single edge vector. A long path sum has a direction essentially
uncorrelated with any one edge, so `<Delta u_e, d_e>` has near-random sign; with `n_s` split edges
needing the same sign simultaneously, the probability of a good flat configuration decays like
`2^{-n_s}` unless the component is short enough for the path sum to be dominated by one or two terms.
Tilings have split-forest components of diameter 1 to 2 edges (`Remark A.4` forest structure plus
symmetry); random Voronoi patches, per K6's first row, had 162 of 391 split edges pointing inward,
which is the signature of near-random signs. So the prediction is not "randomness hurts" but "path
length hurts", and it is testable by *making a tiling worse continuously without randomising its
combinatorics*.

**Kill experiment.** Driver `kill_k8d.cpp`. Two populations. **(i) Jitter ladder:** each of the 8
authored tilings, vertex positions perturbed by `alpha * (mean edge length) * N(0, I)` for
`alpha` in 20 steps from 0 to 0.5, 20 seeds each = **3200 designs**; combinatorics and `sigma` held
fixed, so only geometry moves. For each, `validity_certificate` at `eps = 0.006` and the exact
`Theta_max` by `exact_theta_max_overlap`. **(ii) Cross-generator:** the 200 K1a graphs plus the 16
pattern JSONs in `baseline/kirigami_tessellations`. Fit nothing; just report, for `rho` and for three
competing scalars (`|E_split|/|F|`, `N_bad/|E_split|` from idea 3, split-forest component diameter),
the ROC AUC for predicting "certified deployable".
**Kill rule:** every candidate scalar has **AUC < 0.80** on the combined population, or the jitter
ladder shows **no monotone collapse** (deployable fraction not decreasing in `alpha`, or dropping to
zero already at `alpha = 0.02`, which would mean the tilings are knife-edge and the whole framing is
wrong). Runtime: 3200 certificates at tens of milliseconds plus 216 `Theta_max` evaluations -
**well under 30 min**.

**If it survives, the demo.** One figure that defends the entire negative result: certified
deployable fraction against `alpha`, eight curves collapsing onto one when the x-axis is changed to
`rho`. Plus the sentence a referee cannot object to: "the design space is empty as a function of a
measured geometric ratio, on patterns whose combinatorics we never randomised."

**Risk.** The collapse is not sharp and `rho` explains 60% of the variance, giving a plot instead of
a characterization. Or the tilings turn out to be knife-edge (deployability dies at `alpha = 0.01`),
in which case the honest conclusion is that authored patterns work only in exact symmetry, which is a
*different* and possibly better paper but kills this one as stated.

**Effort.** Derivation S. Code S.

---

## Idea 5 — The hole-area budget: a per-hole necessary condition that names the culprit

**Type:** theorem (necessary condition) + characterization.

**Claim.** For each hole preimage `C` (F11's partition), the hole's area at deployment angle
`theta` is a first harmonic `A_C(theta) = q_C (cos theta - 1) + r_C sin theta` with
`A_C(0) = 0`, and its opening rate decomposes **exactly** into a combinatorially-signed sum:

```
A_C'(0)  =  (1/2) [ sum over hinge edges e in C  of  w_e   +   sum over split edges e in C  of  q_e ],
w_e  =  |det( x_{tgt(e)} - x_{src(e)} , (the other hinge arm at src) )|  >  0 ,
```

so **the hinge edges of a hole always contribute positively and only its split edges can be
negative**. `A_C'(0) > 0` for every hole is *necessary* for `Theta_max > 0`, it is cheaper than the
full certificate by an order of magnitude, and its violation **names the responsible hole**.

**Gap targeted.** 2026 Sec. 6/7 limitation (i) again, and 2026 Sec. 4.2's unsupported "balanced,
small holes" principle, which this makes precise: the relevant balance is hinge-arm area against
split-cut inner products *within one hole*, not hole size.

**Novelty vs field.** (a) `ideas/persona_optimizer.md` derived the hole-area harmonic
`A(theta) = q(cos theta - 1) + r sin theta` and used it for the *second closed angle*; it did not
decompose `r` by edge type and did not use it as an obstruction. (b) field_kirigami **R1**
Grima-Evans compute hole area for rotating squares as a closed form in one pattern; no per-hole
decomposition, no arbitrary graphs. (c) **R12** Dang et al. 2021 require all voids to be
parallelograms in the deployed state (F20) - a *global* structural hypothesis for quads, whereas this
is a per-hole first-order inequality valid on any planar graph.

**Why it might be true / sketch.** The hole boundary at angle `theta` is a closed polygon assembled
from: for each hinge edge in `C`, a wedge of opening `theta` at the shared source vertex, contributing
area `(1/2) |arm1 x arm2| sin theta + O(theta^2)`; for each split edge, a parallelogram with sides
`d_e` and `sin(theta/2) dS_e`, contributing signed area `sin(theta/2) det(dS_e, d_e) = sin(theta/2)
q_e`. Both are first order in `theta` with the stated coefficients, and the shoelace sum of the two
families is the hole area because F11 says the preimages partition `E_hinge` union `E_split`.
Necessity is immediate: `A_C(0) = 0` and `A_C(theta) < 0` for small `theta` is a self-overlapping
hole boundary, hence an overlap.

**Kill experiment.** Driver `kill_k8e.cpp`. For each of the 200 K1a graphs at `X0` from Eq. (6):
(a) verify the identity by comparing the closed-form `A_C'(0)` against a central difference of the
shoelace area of the deployed hole boundary at `theta = 1e-4` (routine: `deploy_basis` +
`holes.hpp`); (b) for graphs with `Theta_max = 0`, check whether the hole minimising `A_C'(0)` is the
hole containing the contact pair reported by `exact_theta_max_overlap`. **Kill rules:** (a) the
identity fails to relative `1e-8` on **any** graph - dead, the decomposition is wrong; (b) the
predicted culprit hole is the actual binding hole on **fewer than 40% of graphs**, in which case the
condition is true but useless as a diagnostic and demotes to a remark. Runtime: 200 graphs, closed
form plus one deployed evaluation each - **under 10 min**.

**If it survives, the demo.** A per-hole heat map of `A_C'(0)` on a Voronoi patch with the binding
contact circled, next to the same map on `snub_square` where every hole is positive. Then the
combinatorial corollary as a table: holes ranked by (number of split edges)/(number of hinge edges),
against `A_C'(0)`, over 200 graphs.

**Risk.** The identity is right but vacuous: `A_C'(0) > 0` holds on essentially every graph while
`Theta_max = 0` anyway, because the binding contact is a **vertex-into-edge** event (K6's measured
mechanism) that does not close a hole. That is the single most likely failure and part (b) of the
kill test is designed to detect it in one run.

**Effort.** Derivation M (the wedge coefficient needs care at reflex corners). Code S.

---

## Idea 6 — The curvature budget: what a pattern can and cannot wrap, in closed form

**Type:** characterization.

**Claim.** For a periodic pattern with unit-cell area `a`, hole-area harmonic `A(theta)` and certified
`Theta_max`, the achievable areal expansion is the closed interval
`[1, 1 + A(Theta_max)/a]`, and a target surface is realizable by that pattern **only if** the areal
distortion required by its best conformal (or authalic) flattening lies inside that interval
pointwise. This is a *closed-form, per-pattern* necessary condition on the realizable shape space,
and it correctly predicts the 2025 paper's published failure that a triangle pattern cannot close a
hemisphere (2025 Fig. F.4) and its published successes on the other target meshes.

**Gap targeted.** 2025 Sec. 13-15 open problems as recorded in `specs/ideator_round2.md` item 7 and
`notes/paper_2025.md`: "realizable shape space per pattern (triangles cannot close a hemisphere, Fig
F.4)" and "automatic tiling selection from curvature". The 2025 paper reports the failure; it does
not predict it.

**Novelty vs field.** (a) field_kirigami **R5/R7** Konakovic et al. use a *fixed* triangular linkage
and measure its expansion range experimentally per-instance, with no per-pattern closed form and no
prediction across pattern families. (b) **R8/R11** Choi-Dudte-Mahadevan design quad patterns for
target surfaces by optimization; the constraint is enforced numerically, not characterized. (c)
**R22** IsoGami is isohedral and periodic like this, but its expansion range comes out of a
continuation solver per structure.

**Why it might be true / sketch.** For a periodic pattern the deployed metric is
`g(theta) = J(theta)^T J(theta)` with `J(theta) = cos(theta/2) I + sin(theta/2) K` (F7/U8,
`periodic_jacobian.hpp`), so `det g(theta) = det J(theta)^2` is an explicit trigonometric polynomial
in `theta` and the areal expansion is `det J(theta)`. Independently, conservation of face area gives
`det J(theta) = 1 + A(theta)/a`, which is a **cross-check the code can run**, and it ties the
macroscopic Jacobian (`periodic_jacobian.hpp`) to the microscopic hole harmonic (idea 5). A surface
of positive Gaussian curvature demands areal *contraction* relative to its flattening in the outer
region and expansion nowhere above the pattern's `det J(Theta_max)`; a pattern whose maximum
expansion is small therefore cannot supply the curvature integral, and `int K dA` over the wrapped
region is bounded by the integrated excess `int (det J(theta(u)) - 1) dA` available.

**Kill experiment.** Driver `kill_k8f.cpp`, using `baseline/kirigami_tessellations` read-only.
For each of the **16** pattern JSONs: build the quotient (`periodic_jacobian.hpp::build_quotient`),
compute `K`, the certified `Theta_max` (`exact_theta_max_overlap`), and `det J(Theta_max)`; separately
compute `1 + A(Theta_max)/a` from the hole harmonic and require the two to agree to `1e-10` (this is
the internal consistency gate). Then for the **4 target meshes** (`hemisphere2`, `pringles2`,
`bumps_plane`, +1), compute the required areal distortion range from a discrete conformal flattening
and form the predicted 16x4 feasibility matrix. **Kill rule:** the predicted matrix disagrees with
the 2025 paper's reported successes/failures on **2 or more of the known entries**, or the internal
consistency gate fails on any pattern. Runtime: 16 quotients and 16 range computations, seconds each;
the flattening is the only new numerical work - **under 30 min**, and it can be cut to the 4 meshes'
areal distortion computed by a simple mass-ratio argument if a full flattening is too slow.

**If it survives, the demo.** The 16x4 predicted-vs-actual table, with the triangle/hemisphere cell
as the hero, and a nomogram: pattern expansion range on one axis, target areal distortion on the
other, every 2025 experiment plotted as a point on the correct side of the line.

**Risk.** The 2025 pipeline uses a *non-uniform* angle field, so the per-cell expansion is not a
single `theta` and the interval bound is too weak to separate the known cases - everything is
predicted feasible. Mitigation is honest: state the condition pointwise in `theta(u)`, which is still
necessary. Second risk: the paper's reported successes/failures are too few and too soft to
constitute a test (only one crisp documented failure), reducing the validation to `n = 1`.

**Effort.** Derivation M. Code M (the flattening is the only unfamiliar piece).

---

## Idea 7 — The usable set is disconnected, and that is why local repair silently re-closes

**Type:** theorem (a negative structural result) with an explicit witness.

**Claim.** The usable region `U(eps)` inside the shape space is **not connected** on named
Archimedean tilings: there exist two certified points `t_A, t_B` in `U(eps)` such that the straight
segment between them leaves `U(eps)`, and more strongly they lie in different connected components.
Consequently any descent method on a smooth penalty - including 2026 Eq. (9) - is **provably
incomplete**, and the measured 9.47% silent re-closure rate (F28/K1c) is a *structural* consequence
rather than a tuning failure of `gamma`.

**Gap targeted.** 2026 Eq. (7)-(9) and Sec. 4.5: the collision energy is presented as a repair step
with no statement about the geometry of the set it searches. `derivations/check.md` T5.2c records
"not basic" as an **unproved** negative and "not by quadrics" as agreed; this supplies the missing
witness for the stronger structural statement.

**Novelty vs field.** (a) field_tutte / `screen_r1.md`: "semialgebraic configuration spaces" is
textbook (Kapovich-Millson, King) and says nothing about connectivity of *this* set. (b)
field_rigidity **R14** Kapovich-Millson prove planar linkage moduli spaces can be arbitrary
manifolds - a universality theorem *about linkages*, not about the collision-free part of a Tutte
auxetic null space, and it gives no witness on a named pattern. (c) **R22** IsoGami's continuation
solver assumes it can walk to a better configuration; a disconnection witness is a direct statement
about when it cannot.

**Why it might be true / sketch.** `U(eps)` is cut out by the requirement that a *first-order-in-eps*
neighbourhood of `theta = 0` be root-free, plus positivity of `|F|` quadratic forms. Both `q_e` and
the corner margins are quadratic and **indefinite**: `q_e(t) = det(dS_e(t), d_e(t))` changes sign when
either factor rotates through the other, and a quadratic inequality `q > 0` with an indefinite form
already defines a **two-component** region in its own right (the interior of a hyperbola's two
branches). Intersecting many such regions generically yields several components. The symmetry of the
Archimedean tilings makes this checkable by hand: on `snub_square` (`dim` shape space 13 per F18's
discussion; `k7_main` records `dim_null = 8`, `dimK = 4` for the 2x2 cell) reflect a certified point
through the pattern's symmetry to obtain a second certified point in a different branch, then test
the segment.

**Kill experiment.** Driver `kill_k8g.cpp`. On `snub_square`, `trunc_square` (4.8.8) and
`hexagons_auto`: sample 2000 shape-space points by Gaussian and by symmetry reflection, keep those
passing `validity_certificate(eps = 0.006)`, and for every certified pair test 64 equally spaced
points on the segment. **Kill rule:** on all three patterns, **every** certified pair is joined by a
fully-certified segment (i.e. `U(eps)` is empirically convex-connected). Then the idea is dead and
`U(eps)` should simply be reported as well-behaved, which is itself worth one sentence. A partial
outcome - segments fail but a two-step polyline always works - demotes the claim from "disconnected"
to "non-convex", which is much weaker and I will say so. Runtime: 3 patterns x 2000 certificates plus
about `10^4` segment tests at milliseconds each - **under 15 min**.

**If it survives, the demo.** A 2D slice of the shape space of `snub_square` with `U(eps)` shaded,
two components visible, the Eq. (9) descent trajectory drawn from a start in the bad component,
terminating at a boundary it cannot cross - the picture that explains the 9.47%.

**Risk.** Finding two certified points in different components requires `U(eps)` to be non-empty in
more than one place, and on the tilings the certified set may be a single blob around the symmetric
configuration. Then the honest answer is "non-convex, connected", a remark not a theorem.

**Effort.** Derivation S (the statement is the experiment). Code S.

---

## Idea 8 — From the LP ray to a finite non-uniform deployment, with a certificate at every step

**Type:** algorithm.

**Claim.** Given a flex `V in P(X)` from idea 1, the predictor-corrector continuation on the pin
constraint variety, with a step accepted only when the exact contact scan certifies no overlap on the
step, produces a **finite** non-uniform deployment of a random planar graph whose *uniform*
deployment has `Theta_max = 0`. Metric: `theta_eff`, the mean hinge opening angle reached, on graphs
where the published pipeline reaches exactly 0.

**Gap targeted.** 2026 Sec. 7 limitation (ii) verbatim (F8): non-uniform hinge-angle deployment,
"could be achieved by sensitivity analysis". 2025 Fig. 10 shows a non-uniform deployment *picture*
for one tiling; it gives no algorithm, no certificate, and does not treat non-2-colourable graphs.

**Novelty vs field.** (a) 2025 Fig. 10 is the closest and it is an illustration; this is a certified
continuation on arbitrary planar graphs. (b) field_kirigami **R22** IsoGami runs contact-aware
continuation with IPC on isohedral tilings - the nearest method, and the delta is (i) arbitrary
planar graphs, (ii) an *exact* contact predicate rather than a barrier, (iii) a per-step certificate
rather than a converged energy. (c) field_rigidity **R22** Tachi 2009 and **R25** Li-Zhu-Qu do rigid
origami continuation; different constraint variety, no cuts, no contacts of duplicate copies.

**Why it might be true / sketch.** The pin constraints are quadratic in the face configurations, so
Newton correction onto the variety is standard; the novelty is only in the acceptance test. Away from
`theta = 0` the mobility is constant on `(0, theta_max)` (U5, 20 of 21 graphs), so the tangent space
does not jump and the continuation is well-posed; and by the no-locking lemma (T2.4, AGREE) every
termination is a contact rather than a kinematic dead centre, so the acceptance test is *complete*:
if no contact fires, the step is legal. The contact predicates along a **non-uniform** step are no
longer first harmonics - that is the honest cost - so the step must be certified by a conservative
swept-volume bound (`contact.hpp::swept_discs` with the exact radius `max(||x||,||chi||)` from
correction 3 of `KILL_REPORT.md`) plus an exact test at the step endpoint.

**Kill experiment.** Driver `kill_k8h.cpp` on **20** K1a graphs with `Theta_max = 0`. Take the LP
flex, continue with adaptive steps, and record `theta_eff` at termination. **Kill rule:** median
`theta_eff < 0.05` rad, i.e. the graphs open by a visually invisible amount and the "non-uniform
rescues random graphs" story is false. Runtime: 20 graphs, a few hundred Newton steps each with a
swept-disc test per step - **20 to 30 min**, and it should be run only after idea 1 passes.

**If it survives, the demo.** Side by side: the same Voronoi patch under the published pipeline
(flat, `Theta_max = 0`) and under the certified non-uniform continuation (visibly open), both exported
through `kiri_export` to laser SVG. Plus the `theta_eff` distribution over 20 graphs.

**Risk.** The continuation walks a few degrees and jams on the same vertex-into-edge contacts that
kill the uniform branch, because the LP flex is only *first-order* separating and the second-order
behaviour turns back. This is the standard first-order-versus-finite gap and it is why the kill
threshold is set at a modest 0.05 rad.

**Effort.** Derivation M. Code L. Gated on idea 1.

---

## Idea 9 — Certified inverse design: replacing the 2025 simulation check, and auditing its results

**Type:** algorithm.

**Claim.** The 2025 pipeline decides deployment-friendliness and `theta_max` by simulation with a
vertex-fusing preprocessing step that F24 shows is unsound (their `merge_close_verts` at 0.1 x average
edge length fuses hinge duplicates and reports 0.067 on `snub_square` where the true first contact is
about 1.65). Replacing that check with `validity_certificate` plus `exact_theta_max_overlap` inside
their optimizer changes the reported feasible pattern/target pairs on a measurable fraction of the
16 patterns x 4 meshes grid, in **both** directions: some pairs they reject are feasible, some they
accept re-close.

**Gap targeted.** 2025 Sec. 4 deployment-friendliness check and its `theta_max`; `notes/repo_2025.md`
Secs. 3.2 and 3.4. This is the 2025-centric direction D9 names.

**Novelty vs field.** (a) 2026's own Eq. (9) is a *penalty* used at `theta = 0`; this is an exact
predicate used as a *constraint oracle* at every angle. (b) **R22** IsoGami uses IPC, i.e. a barrier,
which is a different guarantee (no interpenetration up to solver tolerance) from an exact certificate
with a proved candidate list. (c) `screen_r1.md` records "certified range-optimal projection NOT
FOUND" and that PyKirigami and arXiv:2509.22002 are numerical.

**Why it might be true / sketch.** Their check is a threshold on a fused mesh; ours is the T4.2''
interval scan, which K2a validates against bisection on 187 of 187 configurations with worst error
2e-10 while the naive min-over-roots errs by up to 1.047 rad on 5 of them. Any pipeline whose feasible
set is defined by the weaker predicate must differ on the boundary cases, and F24 shows their
predicate is not merely weaker but wrong by a factor of 25 on a named pattern.

**Kill experiment.** Driver `kill_k8i.cpp` plus a read-only harness over
`baseline/kirigami_tessellations`. For the 16 patterns: compute their `theta_max` (their code, already
built per the baseline parity work) and ours; count disagreements beyond 1e-3. **Kill rule:** fewer
than **3 of 16** patterns disagree materially - then the fusing bug is a curiosity confined to
`snub_square` and there is no audit to write. Runtime: 16 patterns, both pipelines - **under 20 min**,
most of it already-built code.

**If it survives, the demo.** A 16-row table, their `theta_max` against ours, with the disagreements
highlighted and the fused-duplicate mechanism drawn for one case; then the downstream effect on which
pattern the 2025 selection step would choose for each target.

**Risk.** This is an errata paper by the reviewer's own standard R-3, and it must ship as a *section*
inside a paper whose contribution is elsewhere. It also risks being read as hostile to the authors
rather than as a contribution; the framing must be "the exact predicate changes design outcomes",
with the bug as evidence, not as the point.

**Effort.** Derivation S. Code S (mostly harness).

---

## Idea 10 — Permanent incidences: the identically-zero harmonics as an exact symmetry invariant

**Type:** characterization.

**Claim.** The set of candidate pairs whose contact harmonic is **identically zero** - struck by the
identity test before classification, per correction 2 of `KILL_REPORT.md` - is not numerical noise
but a **combinatorial invariant** of `(G, sigma)`: a pair is permanently incident for *all* `theta`
and *all* `X` in the shape space iff the two copies share a source vertex or lie in the same face
(the trivial cases) **or** the pair's `Gamma`-cycle carries a potential relation `Delta u = 0` forced
by `Eq. (2)`. The count of non-trivial permanent incidences is a computable integer, it is **zero on
random graphs and positive exactly on the symmetric tilings**, and it explains why `T4.3-generic` is
false on the papers' own figures.

**Gap targeted.** `derivations/core.md` T4.3-generic, labelled CONJECTURE and stated to be false on
symmetric patterns; nothing in either Segall paper addresses degenerate contact families at all.
`check.md` records the third harmonic class `g(0) = g'(0) = 0` firing on 143 of 1.82M cases.

**Novelty vs field.** (a) field_rigidity **R12/R13** Fowler-Guest and Schulze-Guest-Fowler give
symmetry-adapted mobility counts for symmetric frameworks - the right tool family, and the delta is
that this is an invariant of the *contact* system, not of the rigidity matrix. (b) **R22** IsoGami
works only on isohedral tilings, where such degeneracies are ubiquitous, and never isolates them.

**Why it might be true / sketch.** A harmonic `p + q cos theta + r sin theta` vanishes identically iff
`p = q = r = 0`, three quadratic equations in `X`; on a shape space of dimension `2 |E_split|` that is
codimension 3 and generically empty - which is T4.3-generic. It fails exactly when the equations are
*not independent*, and the potential identity says when: if `Delta u = 0` on the relevant pair then
`dS = 0` identically and the pair never separates, for every `X` in the space. So the invariant is
`#{pairs with Delta u == 0 identically on the shape space}`, computable by a rank test on the
`ZeroPlusForm` blocks rather than by root-finding.

**Kill experiment.** Driver `kill_k8j.cpp`. On the 8 authored tilings and 100 K1a graphs, count
identically-zero harmonics found by `classify_harmonic`, and independently count pairs with
`Delta u == 0` as a *symbolic* condition (all columns of the `GS` block zero to `1e-12`). **Kill rule:**
the two counts differ on any graph beyond the trivially-incident pairs, or the count is positive on
random graphs (which would mean it is not a symmetry invariant at all). Runtime: **under 10 min**.

**If it survives, the demo.** A table: pattern, symmetry group order, permanent-incidence count,
whether `T4.3-generic` holds - resolving a CONJECTURE in `derivations/core.md` with a combinatorial
criterion rather than an appeal to genericity.

**Risk.** The count is positive only on the trivial pairs, i.e. the invariant is `0` everywhere
interesting and the "explanation" of T4.3's failure is elsewhere (for instance in exact edge-length
equalities rather than in `Delta u`). This is a one-section result at best, which is why it is last.

**Effort.** Derivation S. Code S.

---

# Part 3 — Five "obvious" round-2 ideas, and why they are too obvious

**O1. Prove the 2026 Sec. 5.1 conformality claim ("empirically, we observe that for all tilings with
a non-trivial kernel, solving Eq. (13) yields an embedding with conformal deployment").** This is U8
and it is one line from `J(theta) = cos(theta/2) I + sin(theta/2) K`: if the `theta = 0` conditions
of Eq. (13) hold then `J11 - J22` and `J12 + J21` vanish identically in `theta` because they are
linear in `(cos, sin)` with both coefficients killed. The authors wrote "empirically" because they
had not spent the ten minutes; that is exactly the future-work paragraph MISSION Sec. 8 forbids. Keep
it as a two-sentence proposition inside whatever paper ships, never as a contribution.

**O2. The second fully-closed angle in closed form, `theta = 2 arctan(r/q)` from the hole-area
harmonic.** B20, already in the D8 plan. It is a root of a first harmonic that the project already
computes; Grima's locking-angle literature (screened as S3 PARTIAL in `screen_bundle.md`) is adjacent
enough that a referee will ask what is new, and the answer is "the unit cell is symbolic". A corollary,
not a contribution.

**O3. A better orientation heuristic: SDP-rounded max-cut, annealing, or exhaustive search on small
graphs.** The natural next move after K5, and the wrong one, because K5 already showed that a 122x
improvement in the deployability defect bought exactly **zero** range. Optimizing harder against an
objective that has been measured not to correlate with the outcome is the definition of an obvious
idea. Idea 3 above is admissible only because it changes the *objective* to the `0+` sign, which is
the quantity F30 identifies as binding.

**O4. Replace Eq. (9) with a proper IPC / log-barrier solver and report a better range.** K2b already
failed at this against the authors' native optimizer (median gain 0.0000 on 30 designs; theirs better
on 15). IsoGami already ships IPC on this exact object. Running a third barrier method is an
engineering exercise whose best case is a tie.

**O5. Extend the exact contact calculus to per-hinge angles, to 3D, or to curved cuts.** The
half-angle structure that makes everything a first harmonic dies the moment the angles differ, so
"extend it" is either false or requires the genuinely different machinery of idea 8. Stated as an
extension it is a paragraph; stated as idea 8 with a certificate and a kill threshold it is work.
Curved cuts are Liu et al. 2024 (field_kirigami R18).

**Also-obvious, listed for completeness:** publishing `rank(L) = H - 1` and
`H = |E_hinge| - |F| + c(Gamma)` as results (Kirchhoff and Euler, per reviewer point R-3); and
"characterize the shape space of one more tiling family".

---

# Part 4 — Motivating fact and killing fact, per idea

| idea | motivated by (measured) | could be killed in < 30 min by (measured/driver) |
|---|---|---|
| 1 expansive cone LP | F26/F29 mobility median 305 after 2-core; F30 emptiness is a `0+` split/corner event; T2.4 no locking | `kill_k8a`: positive certified LP margin on < 20/100 K1a graphs, or `sigma` reported outside `P(X)` on the 8 tilings |
| 2 emptiness certificate | F25/F30 "0 of 400 samples", which is not a proof; `zero_plus.hpp` exact quadratic forms | `kill_k8b`: certificate on < 25/50 graphs, or **any** certificate on the 8 tilings (soundness bug) |
| 3 sign-aware orientation | F30: `sigma_def` doubled split cuts 126 -> 321 and gained nothing; K6 row 1: 162/391 split edges inward | `kill_k8c`: certified valid 0/200 again, or `N_bad` not halved by local search |
| 4 deployability transition | F25 0/200 random vs 0/8 tiling failures - the generator confound | `kill_k8d`: all candidate scalars AUC < 0.80, or no monotone collapse on the jitter ladder |
| 5 hole-area budget | F11 preimages partition `E_hinge ∪ E_split`; optimizer's `A(theta) = q(cos-1) + r sin` | `kill_k8e`: identity fails at 1e-8 anywhere, or predicted culprit hole right on < 40% |
| 6 curvature budget | 2025 Fig. F.4 documented failure; `periodic_jacobian.hpp` `J(theta) = cI + sK` | `kill_k8f`: >= 2 mispredictions on the 16x4 grid, or `det J(Theta_max) != 1 + A/a` |
| 7 disconnected `U(eps)` | F28 9.47% silent re-closure; `check.md` T5.2c "not basic" unproved | `kill_k8g`: every certified pair joined by a certified segment on all three patterns |
| 8 certified non-uniform continuation | U5 mobility constant on `(0, theta_max)`; T2.4 every stop is a contact | `kill_k8h`: median `theta_eff < 0.05` rad on 20 graphs |
| 9 certified inverse-design oracle | F24 their `merge_close_verts` reports 0.067 vs true ~1.65 on `snub_square`; K2a exact on 187/187 | `kill_k8i`: < 3 of 16 patterns disagree materially |
| 10 permanent incidences | `KILL_REPORT.md` correction 2: 13 630 spurious roots deflated; `check.md` 143/1.82M third-class | `kill_k8j`: counts disagree, or the count is positive on random graphs |

---

# Part 5 — Self-attack: the strongest rejection of my own top three

**Against idea 1 (the expansive cone).**
> The submission's Theorem 1 says that a set defined by strict linear inequalities on a linear space
> is an open polyhedral cone. Its Theorem 2 says that the uniform deployment is a particular element
> of that space. Neither is a theorem. The actual content is an experiment reporting that a linear
> program is feasible on some graphs, and even that is qualified out of existence: the authors concede
> that at convex corners the feasible set is a *union* of polyhedra, that they cannot enumerate the
> branches, and therefore that infeasibility of their chosen branch proves nothing. So the paper can
> report successes and cannot report failures - a method with no negative half. Meanwhile Borcea and
> Streinu have studied cones of expansive motions for periodic frameworks since 2015, and the
> carpenter's-rule literature has studied expansive motions of planar linkages since 2003; the
> submission cites neither in its statement of Theorem 1, and its own screen note admits it has not
> checked them.
>
> Most damning: the flexes the LP returns are *infinitesimal*. The paper's own kill threshold for the
> finite version is 0.05 radians - under three degrees - which is an admission that nobody expects
> these motions to integrate. A first-order rescue of a design space whose failure is measured at
> `theta = 0` is not a rescue; it is the tangent line to the failure.

*Fair rebuttal:* the disjunction is real and I stated it rather than hiding it, and the asymmetry is
the standard one for existence proofs by certificate. The finite-versus-infinitesimal charge is why
idea 8 exists and is gated on idea 1 rather than assumed. The Borcea-Streinu and carpenter's-rule
positioning must be screened before a word is written; if the identical inequality system is there,
idea 1 dies and idea 2 carries the section.

**Against idea 2 (the emptiness certificate).**
> This is the S-procedure. It is in every convex optimization textbook, it is known to be lossy for
> more than one quadratic constraint, and the submission's own kill rule concedes it expects to fail
> on up to half the instances. A "theorem" that holds on 25 of 50 graphs, with no characterization of
> which 25, is a heuristic with a verification step. Furthermore the certificate certifies emptiness of
> the `0+` set only, which is a first-order proxy for the real question; the authors have not shown
> that `0+` emptiness implies emptiness of the usable set, only that it did on the graphs they sampled.
> Finally, the enterprise is self-referential: the submission is building a proof system to shore up
> a negative result that the same submission produced, about a heuristic in someone else's paper.

*Fair rebuttal:* `0+` emptiness *does* imply `Theta_max = 0`, which is the claim - a design with a
non-positive `0+` margin overlaps immediately, and that implication is one line from the `sin(theta/2)`
expansion in `zero_plus.hpp`, not an empirical observation. The "25 of 50" objection is right that
partial coverage is not a theorem, and the correct framing is a *decision procedure that is sound
always and complete sometimes*, reported with its coverage - which is how every practical
Positivstellensatz result is reported. The self-referential charge is fair only if the emptiness
result is the paper; if the paper's positive content is idea 1 or idea 6, this becomes the section
that makes the negative half rigorous.

**Against idea 6 (the curvature budget).**
> That a surface cannot be wrapped by a material that expands by at most a factor `c` unless its
> area distortion is below `c` is not a theorem about kirigami; it is the definition of area
> distortion. The submission's `det J(theta) = 1 + A(theta)/a` is conservation of area. The
> validation is against a single documented failure in a prior paper, and the prior paper's pipeline
> uses a spatially varying angle field, so the submission's single-`theta` interval does not even
> apply to the experiments it claims to predict. Sixteen patterns and four meshes give at most sixty
> four cells, of which the authors can name one that is known to fail.

*Fair rebuttal:* the reviewer is right that the inequality is elementary and the contribution must be
the *closed form for `A(Theta_max)` per pattern*, which requires the certified `Theta_max` that only
this project has, plus the pointwise-in-`theta(u)` statement. If the validation really reduces to one
documented failure, the idea must be demoted to a section of idea 8's paper or dropped. This is the
idea in my top ten most likely to be true and least likely to be interesting, and I rank it fifth
rather than first for exactly that reason.

---

# Part 6 — Ranked list

1. **Idea 1 — the expansive cone LP.** The only idea whose outcome I genuinely cannot predict, and
   both outcomes are results: either non-uniform flexes rescue the graphs the uniform ray fails
   (reframing the emptiness result as a statement about one wrong ray) or they do not (making the
   emptiness structural, about the cut structure rather than the heuristics). Cheap, < 10 min, and it
   answers 2026 limitation (ii) with an algorithm rather than a picture.
2. **Idea 2 — the per-graph emptiness certificate.** The single change that upgrades "0 of 400
   samples" into proofs with named witnesses, killing reviewer point R-1 outright. Sound whenever it
   succeeds; its risk is coverage, not correctness.
3. **Idea 4 — the deployability transition under jitter.** The cheapest idea here and the one that
   most directly defends the negative result against "your generator did the work". Should be run
   first because it costs an hour and everything else's framing depends on its answer.
4. **Idea 8 — certified non-uniform continuation.** The hero demo of the whole project if idea 1
   passes: a graph that the published pipeline cannot open, opened and exported for laser cutting.
   Gated, expensive, high variance.
5. **Idea 6 — the curvature budget.** Most likely to be true, least likely to surprise; the natural
   2025-facing section and the one that connects to D9. Rank depends entirely on whether the 2025
   successes/failures grid is rich enough to be a test.
6. **Idea 5 — the hole-area budget.** A clean necessary condition with a diagnostic use, at real risk
   of being vacuous because K6 measured the binding contact to be vertex-into-edge rather than
   hole-closing. One run settles it.
7. **Idea 7 — disconnected usable set.** A genuinely surprising structural negative with a witness,
   but its success needs `U(eps)` to be non-empty in two places on a tiling, which may simply be
   false.
8. **Idea 9 — certified inverse-design oracle.** Solid, cheap, and the most 2025-centric, but it is an
   errata by the reviewer's own standard and cannot lead a paper.
9. **Idea 3 — sign-aware orientation.** The most likely of the ten to return "0 of 200 again", based
   on K5's precedent; worth one run because it is a day of work and would be a genuine algorithm if it
   fired.
10. **Idea 10 — permanent incidences.** Resolves an internal CONJECTURE and would tidy `core.md`
    T4.3, but its most probable outcome is "count is zero except trivially", which is a footnote.

**If I could run exactly three things tomorrow:** `kill_k8d` (idea 4, one hour, decides the framing of
everything already measured), then `kill_k8a` (idea 1, ten minutes, decides whether the project has a
positive half), then `kill_k8b` (idea 2, twenty minutes, decides whether the negative half is a
theorem or a statistic).
