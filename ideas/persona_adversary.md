# Ideator — persona: ADVERSARY

A hostile, well-read TOG/SIGGRAPH reviewer. My job is (i) to break the "closed-form theory of
uniform deployment" that two other personas converged on, and (ii) to produce ten ideas that are
*surprising and defensible* — ideas whose kill experiment could genuinely fail, and which Segall
et al. would not obviously do next.

**Read (in this order):** `specs/common_preamble.md`, `specs/ideator.md`, `STATE.md` (F1–F14, U1–U5),
`notes/paper_2026.md` (Secs. 2, 8, 10, 12, 13, 15 — the limitations block quoted verbatim at
`notes/paper_2026.md:943–965`), `notes/paper_2025.md` (Sec. 2.11 `θ_max`, Sec. 4), all three field
notes (`field_kirigami.md` 433 lines / 24 rows, `field_rigidity.md` 473 lines / 26 rows,
`field_tutte.md` 288 lines), `papers/related/isogami.txt` (via the R22 row + STATE F12),
`ideas/persona_rigidity.md` (678 lines, in full), `ideas/persona_geometer.md` (Sec. 0 in full,
idea headings + claims for G1–G10, self-attack section).

**Citation hygiene.** Two field tables independently use the labels `R1…R26`. I always write
`kiri-R<n>` for `notes/field_kirigami.md` and `rig-R<n>` for `notes/field_rigidity.md`. Nothing
below is a citation I have not seen in one of those files or retrieved myself in Part 0.

---

## Part 0 — Verbatim web-search log (novelty check on the convergent bundle)

Eight queries and two fetches, run by me in this session. Queries reproduced exactly as issued.

| # | Query / URL (verbatim) | What came back | Bearing on the bundle |
|---|---|---|---|
| Q1 | `rotating rigid units auxetic "rotating squares" general tiling one degree of freedom proof face rotates half the opening angle` | Grima & Evans 2000 and successors; "Two-dimensional systems with only one degree of freedom can be attained for artificially designed interconnected rotating rigid units that are pin-jointed at their corners." **No** statement of the ±θ/2 rotation as a theorem; no general-tiling proof. | (a) is *asserted* 1-DOF folklore for squares/triangles/rhombi; the general statement is not in these results |
| Q2 | `Mitschke Schroder-Turk "packing and self-assembly" OR "auxetic mechanisms" tessellations finite deformation one-parameter family closed form` | Mitschke, Robins, Mecke, Schröder-Turk, "Finite auxetic deformations of plane tessellations", Proc. R. Soc. A (= kiri-R2, rig-R10). Systematic search over a tessellation *archive*; bar-and-joint skeletal model; ν_ss = −1 for two new hexagonal mechanisms. | Nearest *finite* deformation analysis. Bar-joint, not rigid-face; no closed-form deployment map |
| Q3 | `Choi Dudte Mahadevan kirigami tessellations deployment one degree of freedom closed form trajectory opening angle rigid quads` | kiri-R8/R11/R12/R17 + Dang et al. arXiv:2106.15891 + **arXiv:2608.30032 (new to me)**. "quadrilateral tiles connected at hinges with a single global degree of freedom"; "a linear design strategy … allows control of the degrees of freedom in the **deployment angle field**". | The phrase "deployment angle field" in kiri-R17 is a live threat to any non-uniform-angle idea. Quads only |
| Q4 | `Borcea Streinu auxetic deformation periodic framework "one-parameter" deployment cosine sine explicit parametrization rotating units` | rig-R16/R17. "purely geometric notions of auxetic one-parameter deformations"; "a specific planar periodic framework with three degrees of freedom". | Confirms STATE's note that Borcea–Streinu 2020 is one specific 3-DOF framework. **No** cos/sin parametrization found |
| Q5 | `hinged tessellation rigid polygons corner pin joints "alternating" rotation "half the" angle theorem arbitrary planar tiling mechanism trajectory ellipse` | MathWorld "hinged tessellation"; a 2024 *Front. Archit. Res.*-style paper on generating hinged tessellations by adding hinges; **and the Segall 2025 paper itself** (dl.acm.org/doi/10.1145/3757377.3763895). No half-angle theorem. | Nothing |
| Q6 | `kirigami deployment self-contact limit "vertex" into edge closed form maximum opening angle collision rigid tiles analytic` | Top hits are **Segall 2025 itself** ("The maximum opening angle allowed at a vertex is 2π minus the sum of two face angles"), IsoGami, and **PyKirigami (arXiv:2508.15753)**. | **Decisive:** the only closed-form θ_max in the literature is Segall 2025's *local* per-hinge bound. Everything else is numerical |
| Q7 | `"rigid unit modes" OR "rotating rigid units" framework tilings alternating units rotate "equal and opposite" angles arbitrary polygon tiling general theorem Grima Evans review` | RUM literature (Dove et al.); Grima & Evans 2005 "Auxetic behaviour from rotating rigid units": "**rigid polygons** connected together through hinges at their vertices"; "the infinitesimal flex involves rigid units … rotating infinitesimally in **alternating senses**"; Grima's generic non-equilateral triangles. | **The most damaging hit for (a).** The alternating-sense rotation of corner-hinged *rigid polygons* is explicitly the RUM/auxetics model, stated for arbitrary polygon shapes |
| Q8 | `PyKirigami simulator arXiv 2508.15753 deployment kinematics rotating squares generalized tiles degrees of freedom` | Jiang & Choi, PyKirigami, submitted 2025-08-21, rev. 2026-02-16. "models tessellations as **articulated rigid-body networks**"; "collision detection … identify geometric **locking states**". | A 2026 simulator by the nearest competing group still integrates deployment numerically. Nobody has the closed form |
| F1 | `https://arxiv.org/pdf/2106.15891` (Dang, Feng, Duan, Wang 2021, = kiri-R12) | **Fetch failed usefully**: the PDF's text streams did not extract; the fetcher returned "I cannot locate a clearly stated main theorem". | **I could not verify kiri-R12's theorem statement myself.** I rely on the kiri-R12 row, which itself relies on the Dang & Paulino 2025 survey (kiri-R24 Sec. 2.2). Flagged as unverified |
| F2 | `https://arxiv.org/html/2608.30032` (Jiang & Choi 2026, = kiri-R23) | Retrieved. "focus[es] on kirigami structures consisting of quadrilateral tiles with the standard rotating-squares topology"; length-based constrained optimization via IPOPT; **no** cos(θ/2)/sin(θ/2) parametrization, **no** half-angle theorem, **no** θ_max in closed form, **no** conformality analysis. Contributions are an "inertia transposition law" and an "aspect-ratio law". | The most recent competitor (30 Aug 2026) does **not** contain the bundle |

**What Part 0 establishes.** The ±θ/2 alternating rotation of corner-pinned rigid polygons is
standard in the rotating-rigid-units literature (Q1, Q7) and is not a defensible novelty claim.
The *closed form as an explicit linear-in-X map*, and everything downstream of it, is not in any
of the eleven retrievals. So the bundle's novelty lives entirely in "linear in X", not in "±θ/2".

---

## Part 1 — Attack on the convergent bundle (a)–(e)

I take the five claims in the order I was handed them and answer three questions for each:
**is it true? is it new? is it a paper?**

### (a) `ω_f = σ(f)·θ/2`, hence `Y_θ = cos(θ/2)·C(X) + sin(θ/2)·S(X)` with `C, S` linear in `X`

**True.** Both personas verified it numerically (geometer `check.cpp`, 2.0e-15 on a 5×5 patch;
rigidity `rig_check.cpp` V1, 6e-15 on 21/21 graphs including Voronoi patches with split cuts). I
also accept the derivation: it is three lines. Uniform deployability *defines* every hinge to open
by θ; a hinge joins faces of opposite σ (F1); Γ is connected and bipartite with σ as the
2-colouring; therefore `ω_{f1} − ω_{f2} = σ(f1)θ` on every edge of Γ pins `ω` up to a global
rotation. Nothing here is a theorem. **It is the definition of uniform deployability, unpacked.**

**Not new, in the form it is usually stated.** Q7: Grima & Evans 2005 model auxetics as "rigid
polygons connected together through hinges at their vertices" whose units "rotate infinitesimally
in alternating senses". Every rotating-squares/rectangles/triangles/rhombi paper (kiri-R1) writes
the deployed geometry as a function of a single rotation angle applied with alternating sign. A
reviewer from that community will not be impressed by ±θ/2 for one second. If the paper leads
with it, the paper is dead on the first page.

**Two things in (a) *are* new, and both are smaller than they look.**

1. **The split-cut case.** In every rotating-rigid-unit paper I found, the alternating structure
   lives on the tiling: adjacent tiles rotate oppositely. With split cuts, `M` is *not*
   2-colourable (2026 Sec. 5 advertises exactly this: "non-2-colorable (Fig. 1, 9 and 11)"), and
   the geometer's point (0.1) is that the 2-colouring exists on **Γ**, not on `M`, because split
   edges are *deleted* from Γ. That is a genuine observation and it answers
   `notes/paper_2026.md` open question 19. It is worth one paragraph and a figure. It is not
   worth a section.
2. **Linearity in `X`.** `C` and `S` are fixed matrices built from Γ, σ and a spanning tree of Γ.
   Nobody in kiri-R1/R2/R8/R12/R17/R23 states this, because nobody else has a *linear design
   space* to be linear over — Segall 2026 is the first to produce `X = ker[L;B]`. So the leverage
   is parasitic on Segall's own contribution. That is not a flaw; it is the honest framing. But it
   means the sentence "we give a closed form for deployment" is a *lemma of Segall 2026*, not an
   independent result.

**Verdict on (a): a paragraph.** Half a page in a preliminaries section, titled "Lemma 1", with
the split-cut/Γ-bipartiteness observation as the only claim of originality. Any submission whose
abstract sells (a) will be rejected by the auxetics reviewer.

**One correction to both personas.** The geometer writes `Y_θ = cos(θ/2)X' + sin(θ/2)J(CX)`; the
rigidity persona writes `Y_θ = cos(θ/2)C(X) + sin(θ/2)S(X)`. These agree, but note that the
`cos(θ/2)` coefficient is the *duplication map* `X'`, i.e. the same matrix for every design — so
one of the two "linear maps" is trivial and the whole content is in `S`. Stating it as two general
linear maps overstates the structure by a factor of two.

### (b) collision/collinearity events are roots of `P + Q cos θ + R sin θ = 0`, so `θ_max` is closed form

This is the claim I attack hardest, because it is the one aimed at 2026 Sec. 6 limitation 1
(`notes/paper_2026.md:951`, verbatim: "such issues can only be detected through explicit deployment
simulation via forward kinematics … Developing geometric characteristics for favorable deployment
behavior directly from the embedding remains an open problem").

**The algebra is right.** Every deployed point is `cos(θ/2)a + sin(θ/2)Jb`, so every deployed
*difference* is too, and the cross product of two such differences expands with
`cos² = (1+cos θ)/2`, `sin² = (1−cos θ)/2`, `cos·sin = sin(θ)/2` into `P + Q cos θ + R sin θ`
(geometer (0.8)). Setting `t = tan(θ/2)`, this is a quadratic in `t`, so **at most two roots per
ordered pair**, in closed form. Same for the dot product. I have no objection to the algebra.

**Objection 1 — collinearity is not contact, and the personas' claim (b) as stated is false.**
Both personas write "every collision/collinearity event is a root of `P + Q cos θ + R sin θ = 0`".
The vanishing of a cross product says two segments are **parallel**, or that a point lies on a
**line**. It says nothing about whether the point lies on the *segment*. The 2026 paper's own
collision picture (Fig. 15(b), `notes/paper_2026.md` figure table) is the *split-edge* case where
the two duplicates `e′, e″` of one split edge become collinear — for that special pair,
collinearity and contact do coincide, which is exactly why Eq. (7)'s `z_i` construction works.
Generalising from that one figure to "every collision event" is the error. The generic first
contact between two faces that do **not** share a hinge or a split edge is a **vertex into edge
interior**, and its condition is a *system*:

```
   cross( y_b(θ) − y_a(θ) ,  y_v(θ) − y_a(θ) )  =  0            (one trig equation, ≤ 2 roots)
   0  ≤  dot( y_b − y_a , y_v − y_a )  ≤  ‖y_b − y_a‖²          (two trig INEQUALITIES)
```

**Objection 1, part 2 — but the inequalities are also closed form, so the personas are lucky.**
`dot(y_b − y_a, y_v − y_a)` is *also* of the form `P′ + Q′ cos θ + R′ sin θ` (same expansion, dot
instead of cross), and `‖y_b − y_a‖²` likewise. So each root of the equation can be *tested*
against the window in closed form, in O(1). The claim survives — but only after this repair, and
the repair is exactly what neither persona wrote. **A referee who knows computational geometry
will notice within thirty seconds that "collinearity ⇒ collision" is stated and not true, and will
distrust the rest of the file.** The claim must be restated as: *the first-contact angle of an
ordered (vertex, edge) pair is the smallest positive root of a quadratic in `t = tan(θ/2)` that
satisfies two closed-form interval tests.*

**Objection 2 — "closed form for θ_max" is a category error.** Even after the repair,

```
   θ_max(X)  =  min over O(n²) candidate (vertex, edge) pairs of a filtered smallest positive root.
```

A minimum over a quadratically large, combinatorially varying index set is **not a closed form**.
It is an *exact event enumeration*, which is a different and much more modest contribution than
"θ_max is closed form (2026 Sec. 6 limitation 1)". Segall's Sec. 6 asks for "geometric
characteristics … **directly from the embedding**", i.e. a *predicate on `X`* that certifies a
large range without enumerating contacts. An O(n²) exact event sweep does not answer that
question; it replaces bisection with a sweep. That is a *better simulator*, which is precisely
what PyKirigami (Q8) already is, numerically.

**Objection 3 — prior art on the exact-collision half.** kiri-R18 (Liu, Lu, Cao, Deussen, Tu 2024,
*Graphical Models* 133, 101215) — cited by Segall 2026 itself at Sec. 4.6 — already contains, per
the kiri-R18 row, "a closed-form non-adjacent collision condition (Sec. 4.3): no collision occurs
on opposite edges if `length(C1 p1) + length(C3 p2) < 2a cos(θ)`". That is a *conservative
sufficient* bound for named edge pairs in quad Escher tessellations, not an exact root for
arbitrary graphs, and their own conclusion admits the method "still requires manual optimization
of non-adjacent edges". So (b) is a real advance over kiri-R18 — but the delta must be stated as
*exact vs conservative, arbitrary graph vs quads*, and any paper that does not cite kiri-R18 in
the same sentence gets desk-rejected for missing the one reference Segall themselves point at.

**Objection 4 — the local bound is already in the predecessor paper.** `notes/paper_2025.md:187`:
"maximum opening angle allowed at v = 2π − α_i − α_j", and `:194` "θ_max = min over all hinge
vertices v". A closed-form min-over-a-set θ_max already exists in Segall 2025 for the
hinge-adjacent (backside) case. The reader's own flag at `notes/paper_2025.md:958` is the crux:
"**Collision model is local and backside-only.**" So the contribution of (b) is *exactly* the
non-adjacent complement of a formula the predecessor already has. That is defensible — 2026 Fig.
14(b) shows self-intersection at *small* angles, which the local formula cannot see — but it must
be positioned as "we complete Segall 2025 Sec. 4.2", not as "we solve an open problem".

**Objection 5 — the actual payoff is a gradient, and nobody said so.** The reason to want exact
roots is not to evaluate θ_max faster; bisection with a BVH is already fast enough and the paper
already ships it. The reason is that an exact root has an **exact derivative** `∂θ_max/∂X`, which
turns 2026 Eq. (9) (a soft barrier on first-order separation at θ = 0, with an unstated barrier
`B`, an unstated γ, and an unstated diameter rule — see the Reader's flags in STATE) into a
*directly optimisable objective on the shape space*. That, and not the closed form, is the
contribution. Neither persona's write-up leads with it. Both bury it in "Idea 9"/"G10".

**Verdict on (b): a section, and only if reframed.** As "θ_max is closed form", a reviewer writes
*"this is an O(n²) exact collision sweep for a one-parameter rigid-body motion, i.e. a textbook
kinetic-data-structure exercise; the interesting question — a certificate on `X` — is untouched."*
As "an exactly differentiable deployment-range objective that replaces the paper's undocumented
first-order barrier, plus a proof that the first-order barrier has systematic false negatives", it
is publishable. See my idea A1 and A4.

### (c) conformality of periodic deployment is exact iff the θ = 0 condition of Eq. (13) holds

**True, and it is one line.** Given `J(θ) = cos(θ/2)I + sin(θ/2)W` (geometer (0.8)): the conformal
2×2 matrices are the linear subspace `{aI + bJ_rot}`. `cos(θ/2)I` lies in it for every θ. `cos` and
`sin` are linearly independent as functions of θ. Therefore `J(θ)` is conformal for all θ **iff**
`W` is conformal. And `dJ/dθ|_{θ=0} = ½W` because `d/dθ cos(θ/2)|_0 = 0`. So Eq. (13)'s
infinitesimal condition **is** the condition "W conformal", verbatim, not merely a necessary
consequence of it.

**So (c) is a proof of a stated empirical claim** — `notes/paper_2026.md:820` quotes it verbatim:
"once the derivative of the conformal distortion vanishes at θ = 0, the deployment remains
conformal for all opening angles θ", with the Reader's note "**stated without proof or even a
stated experiment count**". Closing that is legitimate. It is also, unavoidably, **a remark**.

**The killing observation.** (c) is *strictly weaker than it sounds*, in a way both personas miss.
It says: given that `W` is achievable, conformality is preserved. It says **nothing about which
`W` are achievable** — i.e. about the image of the design map `X ↦ W(X)` restricted to
`X ∈ ker[L;B]`. That image is the real object. It controls the achievable finite Poisson ratio
curve `ν(θ)`, it decides whether a given periodic graph can be made isotropic at all (2026 Sec.
5.1's Eq. (13) may have no nondegenerate solution and the paper says only "empirically, we
observe … for all tilings with a non-trivial kernel"), and it is a *linear-algebraic* question
because `W = Q P_0^{-1}` with `Q`, `P_0` both linear in `X` — so `X ↦ W` is a **matrix-valued
rational map of bidegree (1,−1)**, whose image is a constructible set. Nobody has looked at it.
That is my idea A5, and it is worth more than (c) by an order of magnitude.

**Verdict on (c): a paragraph.** Two sentences in a "we explain the paper's empirical observation"
remark. A reviewer will write: *"Section 5 proves that a line segment in a linear subspace stays
in the subspace."* Correct, and fatal if it is a headline.

### (d) `H = |E_hinge| − |F| + c(Γ)`

**True and useful; the easy half of the right question.** It is Euler's formula on Γ once you know
that holes ↔ bounded faces of Γ. The geometer got the `c(Γ) − 1` correction only after the
numerics caught the error (400/400 failures without it, 0/400 with), which is a point in the
file's favour, not against it. It proves 2026 Sec. 4.4's third claim (`H ≤ V_int`, equality iff
`E_split = ∅`) as an identity and it replaces Alg. 1's flood fill with a formula. Also it corrects
`notes/field_tutte.md` (2.1) and STATE U1, both of which drop the `c(Γ)` term.

**Why it is the easy half.** 2026 Sec. 4.4's *load-bearing* claim is not `H ≤ V_int`. It is
`#independent equations = #holes`, i.e. `rank(L) = H`, because that is what fixes
`dim(shape space) = 2(N − rank L)` and therefore the entire "we characterise the full space of
embeddings" claim in the abstract. STATE F14 already shows this claim is **false as stated for
periodic patterns** (`1ᵀL = 0` ⇒ `rank(L) ≤ H − 1`). Counting `H` does not touch `rank(L)`. Both
personas know this — geometer G5 goes at `rank(L) = H − dim Z` — but (d) as handed to me is the
part that follows from Euler, and Euler is not a contribution.

**Verdict on (d): a paragraph**, inside a lemma block, subordinate to a real rank theorem.

### (e) infinitesimal mobility `= dim ker A − c(Γ)`, `A ∈ R^{2b₁(Γ)×|F|}`, measured mobilities ≫ 1

Three attacks, in increasing severity.

**Attack 1 — the headline number is a boundary artefact, and I am fairly confident of it.**
Reported (STATE U4): mobility 1–3 on Voronoi patches, **19–26 on random Delaunay patches with a
free boundary**. A Delaunay patch of ~50 faces has a large boundary; a face attached to Γ by a
*single* hinge edge is a free pendulum contributing exactly 1 to `dim ker A` and nothing to the
mechanism's interior behaviour. Chains of such faces contribute one each. The correct statement
is almost certainly

```
   dim ker A(Γ)  =  |F ∖ core₂(Γ)|  +  dim ker A(core₂(Γ)),
```

where `core₂` is the 2-core (iteratively delete degree-≤1 nodes of Γ). The rigidity persona's own
Idea 7 proposes exactly this decomposition — which means the file simultaneously advertises
`m = 19–26` as "answering 2026 Sec. 6 limitation 2 numerically" (STATE U4, verbatim in the file's
Sec. 0) and proposes the theorem that explains that number away. A reviewer will catch the
tension. **Prediction I would stake the idea on: `dim ker A(core₂)` is 1–3 for every one of those
Delaunay patches.** If so, the headline "measured mobilities ≫ 1" is a statement about dangling
tiles, not about kirigami. This is trivially testable and must be tested *before* the number is
put in an abstract.

**Attack 2 — it is the mobility at the wrong configuration, by the file's own evidence.**
2026 Sec. 6 limitation 2 (`notes/paper_2026.md:955`) asks for "the degrees of freedom of the
deployment space", i.e. of the configuration-space branch the structure actually travels on. The
rigidity persona's own V4 measured that mobility **drops discontinuously at θ = 0** on 20 of 21
graphs (squares 5×5: 10 → 1; kagome: 18 → 7; periodic squares: 4 → 1). So `dim ker A` at the flat
state is *not* the DOF of the deployment space; it is the dimension of the tangent cone at a
singular point, and it over-counts by up to a factor of ten. The number that answers Sec. 6
limitation 2 is `m(θ)` for `θ ∈ (0, θ_max)`, which the file also measured — and which is 1 for
5×5 squares and periodic squares, i.e. **exactly the classical rotating-squares answer**. So the
honest headline is the opposite of the one given: on the deployed branch these are nearly always
low-DOF mechanisms, and the interesting cases are the ones where they are not.

**Attack 3 — the "constant on (0, θ_max)" claim is unfalsifiable by the method used.** U5 asserts
`m(θ)` is constant on the open interval, from five sampled angles. Rank drops on a **measure-zero**
set. Sampling five angles is guaranteed to miss every bifurcation with probability 1, so the
observation is exactly what a bifurcating structure would also produce. U5 is not evidence; it is
the null result of a test with zero power. The only correct method is to solve for the parameter
values where the rank drops, which is a generalised-eigenvalue problem in `t = tan(θ/2)` because
`A(θ) = cos(θ/2)A_c + sin(θ/2)A_s` is a **linear matrix pencil** — the rigidity persona's Idea 5.
That is the best idea in either file and it is ranked fifth in its own file.

**A fourth point, in the bundle's favour, that neither persona states.** `σ ∈ ker A(Y_θ)` for
*every* θ along the uniform path, because the uniform motion is a valid velocity at every
configuration on the path. Therefore **the uniform branch never locks**: there is no dead-centre,
no jamming, no configuration where the structure cannot continue to open. Every stopping event is
a *collision*, never a kinematic singularity. This is a genuine theorem, it is two lines, it is
the correct rebuttal to anyone (PyKirigami, Q8, advertises "identify geometric locking states")
who expects locking, and it is what makes (b)'s enumeration *complete* rather than merely
*necessary*. Neither persona noticed it. It should be Lemma 2 of any paper built on this bundle.

**Verdict on (e): a section with a bug.** The algebra is right and the bilinear-form observation
(`L = B(σ,·)`, `A = B(·,X)`, rigidity persona's (R4)) is genuinely elegant and genuinely absent
from both Segall papers and from IsoGami (kiri-R22, which computes mobility as zero eigenvalues of
`JᵀJ` — a number, not a structure). But the reported mobilities are contaminated by dangling
faces and are measured at a singular configuration, and the stability claim rests on a test that
cannot fail.

### Which of these is "a paragraph, not a paper"?

| Claim | Status | Verdict |
|---|---|---|
| (a) ±θ/2 and the cos/sin closed form | true; ±θ/2 is RUM folklore (Q1, Q7); linearity in `X` is new but parasitic on Segall | **paragraph** (Lemma 1) |
| (b) closed-form contact events / θ_max | algebra right, statement false as written (collinearity ≠ contact), "closed form" oversold, kiri-R18 + Segall 2025 Sec. 4.2 are prior art for the two halves | **section**, only as a differentiable objective |
| (c) conformal iff conformal at θ = 0 | true, one line, proves a quoted unproved claim | **paragraph** (Remark) |
| (d) `H = \|E_hinge\| − \|F\| + c(Γ)` | true, Euler; corrects field_tutte and STATE U1; but `rank(L)`, not `H`, is the load-bearing quantity | **paragraph** (Lemma 2) |
| (e) mobility `= dim ker A − c(Γ)` | true; headline numbers are boundary artefacts at a singular configuration; stability claim untested | **section**, after a 2-core reduction and a re-measurement at θ > 0 |

Individually: five paragraphs and two damaged sections. **The bundle as handed to me is not a
paper.** A TOG reviewer reading it as an abstract writes: *"The authors observe that a
one-parameter rigid-body motion of a corner-hinged tiling is trigonometric in the parameter, and
draw the standard consequences. The main claims are either folklore in the auxetics literature
(the alternating half-angle rotation), immediate from Euler's formula (the hole count), or
one-line linear algebra (the conformality preservation). Reject."*

### The smallest bundle that IS a paper

One title, one theorem, one algorithm, one measured margin. My proposal:

> **"The deployment path is a conic: exact contact analysis and range-optimal Tutte auxetic
> embeddings."**

with exactly this content:

1. **§3 Preliminaries (2 pages).** Lemma 1 = (a), with the split-cut/Γ-bipartiteness observation
   and explicit credit to the RUM literature (kiri-R1, Q7) for the ±θ/2 half. Lemma 2 = (d).
   Lemma 3 = **the no-locking lemma** (`σ ∈ ker A(Y_θ) ∀θ`), which is what licenses everything
   after. This is where (c) goes, as a two-sentence remark.
2. **§4 The exact contact calculus (the theorem).** Each ordered (vertex, edge) pair contributes at
   most 2 candidate angles, roots of a quadratic in `t = tan(θ/2)`, filtered by two closed-form
   interval tests. Then the part that makes it a theorem rather than an enumeration: **a
   certified active-set bound** — a proof that only pairs within a computable geometric
   neighbourhood can be first-contacts, reducing O(n²) to O(n) with a certificate (my idea A4).
   Without this, §4 is an exercise.
3. **§5 The indictment (the surprise).** 2026 Eq. (9) tests the *first-order* separation rate at
   θ = 0. With the closed form, the exact separation of a pair is `P + Q cos θ + R sin θ`, whose
   derivative at θ = 0 is `R` — **independent of `P` and `Q`**. So a design can pass Eq. (9) with a
   large positive rate and still close at a moderate θ. Measure the false-negative rate of the
   paper's own heuristic on 500 graphs. This number is currently unknown to anyone and could come
   back at 2% (idea dies quietly) or at 40% (paper writes itself). **This is the only part of the
   bundle whose answer is not guessable in advance.** (My idea A1.)
4. **§6 The algorithm.** Maximise the exact `θ_max` over `X = ker[L;B]` using the exact
   `∂θ_max/∂X` from the closed-form roots, with the active set from §4. Compare against Segall's
   Eq. (6) + Eq. (9) pipeline on ≥ 20 graphs. The claim to beat is a large margin in achieved
   collision-free range at equal geometric-distortion budget.
5. **§7 The mechanism companion.** The bifurcation pencil (rigidity Idea 5) as the answer to
   "is the uniform branch the one you get?", *after* the 2-core reduction fixes the mobility
   numbers. One page, honest, with the U5 methodology corrected.

(e)'s mobility theory is **a second paper**, not part of this one, and it needs the 2-core
experiment run before anyone writes its abstract.

---

## Part 2 — Five ideas that are TOO OBVIOUS (and exactly why)

The test I apply: *could a competent reader of both Segall papers have written this as a
future-work paragraph, with the result already known?* If yes, it is not an idea, it is a task.

**O1. "Extend the framework to non-uniform hinge angles via sensitivity analysis."**
2026 Sec. 6 limitation 2, `notes/paper_2026.md:955`, **names the method**: "Analyzing the degrees
of freedom of the deployment space, which could be achieved by sensitivity analysis, is an
interesting direction to explore." Writing the paragraph the paper already wrote is not novelty.
Worse, kiri-R17 (Dudte, Choi, Becker, Mahadevan 2023) already has a "deployment angle field" with
controlled DOF, obtained by "elementary linear algebra" (Q3), for quads. Any non-uniform-angle
idea must beat kiri-R17 on a class kiri-R17 cannot touch, and must say so in its first sentence.
Both the rigidity persona (Idea 6) and the geometer (G9) propose non-uniform-angle spaces without
citing kiri-R17 at all. That alone would sink both.

**O2. "Replace the Goemans–Williamson relaxation for σ with an exact or better max-cut solver."**
2026 Sec. 4.2 / Eq. (1) is openly a heuristic, so it looks like free points. It is not, for two
reasons. First, the paper says (`notes/paper_2026.md:619`) that σ "may simply be specified
manually" and the UI lets the user click faces — so σ is not treated as an optimisation the
authors care about winning. Second and fatally, **the objective is unmotivated**: minimising split
cuts is a proxy for connectivity, not for anything measured in the results. Beating a proxy by a
better solver produces a table nobody can interpret. The only interesting version of this — which
graphs admit *no* good σ — is a different idea (my A10), and it is not about the solver.

**O3. "Relax the rigid-planar-face assumption / add thickness / add stiffness."**
2026 Sec. 6 limitations 3 and 4, both quoted verbatim in `notes/paper_2026.md:959` and `:963`,
including the suggested references (Warisaya et al. 2022, Zaman et al. 2025). Also inadmissible
under `specs/ideator.md`: it is a mechanics/fabrication direction, not a statement about the
design space. It fails the "computational at the core" constraint outright.

**O4. "Apply the framework to aperiodic / quasicrystalline / Penrose tilings."**
2026 Sec. 5's own opening sentence (`notes/paper_2026.md:767`) advertises this as *already done*:
"Our analysis applies to arbitrary planar graphs, including those that are non-periodic (Figures 1
and 3), non-2-colorable (Fig. 1, 9 and 11)". Running the existing pipeline on a hat monotile is a
figure, not a result. The one non-obvious version — whether aperiodicity *forces* the shape space
to behave differently, e.g. whether `rank(L)` has a growth law different from the periodic case —
is a special case of my A7 and should be posed there, not as its own paper.

**O5. "Optimise curved cuts for collision-free deployment."**
2026 Sec. 4.6 states the problem *and names the solution*: "unlike straight-line cuts, curved cuts
are more prone to collision and may require dedicated optimization to ensure collision-free with
large range of deployment **[Liu et al. 2024]**" (`notes/paper_2026.md:965`). kiri-R18 is that
paper and it already contains a closed-form non-adjacent collision condition for the curved-quad
case. Doing it again on general graphs is engineering. Additionally the constraint set is trivial
— Eq. (10) is a monotone-radius condition, so the admissible curve space is a convex cone in each
cut — which means the optimisation is easy and the paper would be a UI demonstration.

**A note on what "obvious" does not mean.** The paper's *unproved* claims (Sec. 4.4's
`rank(L) = H`, Sec. 5.1's "remains conformal for all θ") are obvious *targets* but their answers
are not obvious — indeed STATE F14 shows Sec. 4.4 is false as stated for periodic patterns. An
obvious question with a surprising answer is a fine paper. An obvious question with the expected
answer is O1–O5.

---

## Part 3 — Ten ideas

Throughout: `M = (X, F)` the uncut planar graph; `Γ = (F, E_hinge)` the hinge graph;
`X = ker[L; B]` the shape space (2026 Eq. (5)); `Y_θ = cos(θ/2)·X' + sin(θ/2)·J·S X` the closed
form (geometer (0.7), rigidity V1); `t = tan(θ/2)`; `A(θ)` the `2b₁(Γ) × |F|` mobility matrix
(rigidity (R2)). All ideas assume the reference pipeline in `code/src/core/` exists as specified in
`specs/builder_core.md` (mesh, cut, hole preimages, `[L;B]` + null space, forward kinematics,
θ_max by collision bisection).

---

### A1. The collision barrier of Eq. (9) has a measurable false-negative rate

**Type:** characterization + algorithm.

**Claim.** For a split-edge pair `(i,j)`, the exact signed separation along the deployment path is
`s_{ij}(θ) = P + Q cos θ + R sin θ` with `P + Q = 0` (contact at θ = 0) and `s′_{ij}(0) = R`.
2026 Eq. (9) penalises only `sign(R)`. Therefore **the set of designs that Eq. (9) certifies but
which self-intersect at some `θ < θ_max` is exactly `{X : R > 0 and the second root
`t = −R/(2Q)`… lies in (0, tan(θ_max/2))}`, a nonempty semialgebraic subset of `X`**, and its
relative volume — the heuristic's false-negative rate — is a computable number that is currently
unknown.

**Gap targeted.** 2026 Sec. 4.5, Eqs. (7)–(9), the collision-avoidance energy. The paper's own
Sec. 6 limitation 1 (`notes/paper_2026.md:951`, verbatim): "self-intersections may occur at small
opening angles (e.g., Fig. 14). At present, such issues can only be detected through explicit
deployment simulation via forward kinematics." Plus the Reader's flags in STATE: the barrier `B`
has "**no formula, no parameters, no smoothness statement**", γ and the diameter-selection rule are
absent from the paper, and the regularizer is printed unsquared.

**Novelty vs field.** (i) kiri-R18 (Liu et al. 2024) has a *conservative sufficient* non-adjacent
collision inequality for quad Escher tiles — it errs on the safe side, so it has false *positives*,
never false negatives, and it says nothing about Eq. (9). (ii) kiri-R22 IsoGami handles collisions
by IPC barriers inside a nonlinear continuation — a numerical remedy with no characterization of
when a first-order test fails. (iii) rig-R19 (Kumar & Pellegrino 2000) detects bifurcation
numerically for pin-jointed bars, unrelated to contact. None of the three asks "how often is the
published first-order test wrong", which is the question.

**Sketch.** For the two duplicates `e′, e″` of a split edge, both endpoints are of the form
`cos(θ/2)a + sin(θ/2)Jb`, so their cross product is `P + Q cos θ + R sin θ` by geometer (0.8), and
at θ = 0 the duplicates coincide, forcing `P + Q = 0`. Write `u = P + Q cos θ + R sin θ` in
`t = tan(θ/2)`: with `cos θ = (1−t²)/(1+t²)`, `sin θ = 2t/(1+t²)`,

```
   (1 + t²) u(t) = (P + Q) + 2R t + (P − Q) t²  =  2R t + 2P t²  =  2t (R + P t).
```

So (using `P + Q = 0`) the separation has a **double structure**: the root at `t = 0` is the flat
state, and the *only other* root is `t* = −R/P`. Eq. (9) enforces `R > 0`. A re-closure at a
positive angle therefore happens **iff `P < 0`**, and then at exactly `θ* = 2 arctan(−R/P)`.
Eq. (9) never looks at `P`. This is not a heuristic weakness that might exist; it is a
one-parameter family the energy is blind to by construction. The unknown is only *how often*
`P < 0` and `θ* < θ_max` for embeddings that Eq. (9) actually returns.

**Kill experiment.** `gen_voronoi` / `gen_delaunay` / `gen_perturbed_grid`, 500 graphs, 100–800
faces. For each: assign σ, solve Eq. (4)/(6) for `X_0`, run the *reference implementation of
Eq. (9)* to convergence, then for every split-edge pair evaluate `P, R` in closed form and count
pairs with `P < 0, R > 0, θ* < θ_max`. Cross-check every reported `θ*` against the existing
`forward_kinematics` + `collision_bisection` routine to ≤ 1e-9. **Kills the idea if the fraction of
Eq.-(9)-certified designs with at least one re-closure below `θ_max` is < 3%.** Runtime: one
`[L;B]` solve plus one Newton solve plus `O(|E_split|)` closed-form evaluations per graph; well
under 30 min.

**If it survives, the demo.** A histogram of `θ*/θ_max` over the 500 graphs; a scatter of "range
predicted by Eq. (9)" vs "true range" with the diagonal, showing the systematic bias; the drop-in
replacement (penalise `min(1, −P/R)` instead of `sign(R)`, which is smooth and has the same cost)
on ≥ 20 graphs with a range-vs-distortion Pareto front. Hero: reproduce 2026 Fig. 14(b) — the
paper's own picture of collisions at small angles — and predict the collision angle to machine
precision without simulating.

**Risk.** The authors' unpublished implementation may already include a second-order or multi-angle
sampling term that the paper does not describe (the Reader flagged that `B`, γ and the diameter
rule are all absent from the text). Then the "false-negative rate" is a property of the *paper*,
not of the *code*, and a reviewer with access to `github.com/segaviv/tuttekiri` will say so. This
must be checked against that repository before committing.

**Effort.** Derivation S. Code M.

---

### A2. The uniform branch never locks — but it bifurcates, and the angles are a pencil

**Type:** theorem.

**Claim.** (i) **No-locking lemma:** `σ ∈ ker A(Y_θ)` for every `θ`, so the uniform deployment
never reaches a kinematic dead-centre; every terminating event is a contact. (ii)
`A(θ) = cos(θ/2) A_c + sin(θ/2) A_s` is a **linear matrix pencil**, so the set of angles at which
`dim ker A` jumps is the set of real `t = tan(θ/2)` with `rank(A_c + t A_s) < max_t rank(...)`,
i.e. the real points of the determinantal variety of the pencil `(A_c, A_s)` — **finitely many, at
most `min(2b₁(Γ), |F|)`, computable by one generalized eigenvalue / Kronecker canonical form
computation.** (iii) For a positive fraction of Segall-2026 designs at least one such `θ_b` lies
strictly inside `(0, θ_max)`.

**Gap targeted.** 2026 Sec. 6 limitation 2 (`notes/paper_2026.md:955`), and — more sharply —
the *unstated assumption* running through both papers and through IsoGami that the deployment
path is unique. STATE U5 (rigidity persona) asserts mobility is constant on `(0, θ_max)` from five
sampled angles; that test has zero power against a measure-zero bifurcation set, so the question
is genuinely open in the file that raised it.

**Novelty vs field.** (i) rig-R19 Kumar & Pellegrino 2000 detect bifurcation of pin-jointed bar
structures **numerically**, by path-following with a perturbation; they have no pencil because
their configuration is not a one-parameter trigonometric curve. (ii) rig-R9 Connelly & Servatius
1994 give higher-order rigidity and configuration-space cusps — the right warning, no computation.
(iii) kiri-R22 IsoGami runs compliant-mode continuation and *filters for exactly one DOF*; a
continuation solver silently follows whichever branch its predictor lands on, so IsoGami would
cross a bifurcation without reporting it. Nobody has an exact bifurcation spectrum for a kirigami
deployment, and Q8 shows the state of the art (PyKirigami, Feb 2026) still detects "locking
states" numerically.

**Sketch.** `A(ω)_z = Σ_i ω_{f_i}(p_{i−1} − p_i)` (rigidity (R2)), linear in the hinge points `p`.
Each `p_i(θ) = cos(θ/2) a_i + sin(θ/2) J b_i`, so `A(θ) = cos(θ/2)A_c + sin(θ/2)A_s` exactly, with
`A_c, A_s` assembled once from `X`, σ and Γ. Rank of `cos·A_c + sin·A_s` is homogeneous of degree
0, so factor out `cos(θ/2)` and study `A_c + t A_s`. For (i): the uniform motion is a legitimate
velocity field at *every* `Y_θ` (differentiate the closed form), and `A(Y_θ)σ = 0` by (R4)/V3
evaluated at `Y_θ` rather than at `X`. So `σ` is always in the kernel; the branch cannot terminate.
For (ii): generic rank of a pencil is attained off a proper algebraic subset; the exceptional `t`
are the roots of the gcd of the maximal minors, obtainable from the Kronecker canonical form of
`(A_c, A_s)` or, when `A_c` has full column rank, as the finite generalized eigenvalues of the pair.

**Kill experiment.** 500 graphs from the same three generators plus the six named tilings
(squares, triangles, hexagons, kagome, (3,4,3,12), the Fig. 11 pentagon/hexagon graph). For each:
build `A_c, A_s` from the solved `X_0`; compute the real `t` where rank drops, by SVD of
`A_c + tA_s` on a fine sweep **cross-validated against** the generalized eigenvalues of
`(A_c, −A_s)` restricted to a maximal-rank square subblock; compare each `θ_b` against `θ_max` from
`collision_bisection`. Two ways to die: **(α)** if `rank(A_c + tA_s)` is constant for all real `t`
on every one of the 500 graphs, the bifurcation set is empty and part (iii) is dead — the theorem
in (i)/(ii) survives but reduces to a remark; **(β)** if the sweep and the eigenvalue computation
disagree on any graph, the pencil formulation is wrong and everything dies. Cheap: `2b₁ × |F|`
matrices, ~500–3000 columns, one SVD per sample.

**If it survives, the demo.** For ≥ 20 graphs, a one-dimensional "deployment spectrum" plot: the
interval `[0, θ_max]` with collision events (A1/A4) in red and bifurcation angles in blue, showing
which designs bifurcate before they collide. Hero: a hand-authored pattern with `θ_b ≪ θ_max`,
fabricated, where the physical sheet visibly leaves the uniform mode — a claim the geometry makes
in advance and the object either confirms or refutes.

**Risk.** The most likely failure is outcome (α): the rank could be constant because `σ` spans the
kernel generically and the special positions that create extra flexes may be exactly the flat state
`t = 0` and nothing else. STATE U5 is weak evidence *for* this risk. That is why this is the right
experiment: it is genuinely 50/50, it runs in an afternoon, and either answer is worth stating.

**Effort.** Derivation M. Code M.

---

### A3. Fully-closed patterns are exactly root-coincidence, not angle-averaging

**Type:** characterization + algorithm.

**Claim.** A pattern is fully closed at `θ_max` **iff all first-contact roots coincide**:
`t*_1 = t*_2 = … = t*_k` over the active contact set, where each `t*_i` is the closed-form root of
A1/A4. 2026 Eq. (14) instead penalises `(α_lik − ᾱ)²` toward a *running average* plus an
edge-length mismatch — a surrogate that is neither necessary nor sufficient for the exact
condition, and I claim it is provably **not sufficient**: there exist embeddings with all
`α_lik` equal and all named edge lengths equal that still have a gap at `θ_max`, because the
closing pairs that matter are non-adjacent and Eq. (14) only sees hinge-adjacent pairs.

**Gap targeted.** 2026 Sec. 5.2 / Eq. (14). The paper's own text (`notes/paper_2026.md`, Sec. 10.2)
says the `β_i = β` condition is necessary, then adds a second condition, and — quoting the Reader —
"**no proof is given that the two conditions together are sufficient, only that the first alone is
not.**" I claim the pair is in fact *not* sufficient and give the exact condition.

**Novelty vs field.** (i) kiri-R11 Choi, Dudte, Mahadevan 2021 "Compact reconfigurable kirigami"
constructs patterns with two compact states, for **quads**, by a per-linkage construction — not a
condition on an arbitrary graph. (ii) kiri-R23 Jiang & Choi 2026 solve for multiple compact states
by length-based constrained optimization with IPOPT (F2) — again quads/rotating-squares only, and
numerically, with "angle inequality constraints to prevent self-intersection", not an exact closure
condition. (iii) kiri-R17's additive construction grows patterns that are compact by construction
and never characterises the closure condition on a given graph.

**Sketch.** "Fully closed at `θ_max`" means: at `θ = θ_max` the deployed structure has zero total
hole area. Hole area is a signed polygon area of deployed vertices, hence
`Area_h(θ) = quadratic in (cos(θ/2), sin(θ/2)) = P_h + Q_h cos θ + R_h sin θ`. So

```
   fully closed at θ*   ⟺   P_h + Q_h cos θ* + R_h sin θ* = 0  for every hole h,
   with θ* = min over the active contact set of the first-contact roots.
```

Two exact systems, both trigonometric-quadratic in `t*` and both **quadratic in `X`** (because
`P, Q, R` are quadratic in `X` by geometer (0.8)). Eq. (14) is a first-order proxy for the first
system evaluated only on hinge-adjacent pairs. The counterexample to sufficiency should be
constructible by hand: take a pattern with all `α_lik` equal (a regular tiling, so Eq. (14) is at
its global minimum) and perturb *inside* the shape space in a direction that leaves every
hinge-adjacent angle fixed — which exists whenever `dim X > #{α constraints}` — and check whether
a hole area stays positive at `θ_max`.

**Kill experiment.** Two parts. **(1) Sufficiency counterexample:** on the Fig. 11 graph (shape
space known to be 6-dimensional, STATE F9) and on (3,4,3,12) (2026 Fig. 21, a named reconstructible
tiling), enumerate 10⁴ random points of `X` with Eq. (14)'s residual below 1e-10 and evaluate total
hole area at `θ_max` in closed form. **If every such point has hole area < 1e-8, Eq. (14) is
sufficient after all and the "not sufficient" half is dead.** **(2) Algorithm:** solve the exact
root-coincidence system by Gauss–Newton on 20 graphs and compare the residual hole area at `θ_max`
against Eq. (14)'s solutions. **If Eq. (14) achieves hole area within 5% of the exact solve on all
20, the algorithmic half is dead too.**

**If it survives, the demo.** A table of residual hole area at `θ_max`, Eq. (14) vs exact, on
20 graphs; the counterexample rendered as a pair of deployed states that look identical at θ = 0
and differ visibly at `θ_max`; a fabricable export of a fully-closed-at-both-ends pattern on a
graph Segall's Fig. 2 does not contain — ideally a non-2-colourable one, which no prior method
(kiri-R11/R17/R23, all quads) can produce.

**Risk.** Eq. (14) may be sufficient *in practice on the tilings the authors tried*, because
regular tilings have so much symmetry that the hinge-adjacent conditions already pin everything.
Then the counterexample exists only on irregular graphs and the result reads as pedantic. Mitigate
by running part (1) on Voronoi patches first, where symmetry is absent.

**Effort.** Derivation M. Code M.

---

### A4. A certified O(n) active contact set (the locality theorem)

**Type:** theorem + algorithm.

**Claim.** Let `r_f` be the circumradius of face `f` and `d(f,g)` the flat-state distance between
face centroids. Then no pair `(f,g)` with `d(f,g) > κ·(r_f + r_g)` can produce the *first* contact
of the structure, for an explicit constant `κ` depending only on `θ_max ≤ π`; consequently the
active contact set has size `O(n)` with a **certificate**, and exact `θ_max` is computable in
`O(n log n)` rather than `O(n²)`, with a proof that no event is missed.

**Gap targeted.** 2026 Sec. 6 limitation 1 verbatim: "Developing geometric characteristics for
favorable deployment behavior **directly from the embedding** remains an open problem." An O(n²)
sweep does not answer that; a *certificate on the embedding* that bounds which pairs can ever meet
is the first honest step toward it. Also the exact form of the objection I raised against the
convergent bundle's claim (b) in Part 1.

**Novelty vs field.** (i) kiri-R18 Liu et al. 2024 have a *sufficient* non-adjacent collision
inequality but only for named opposite-edge pairs of quad Escher tiles, and their own conclusion
says the method "still requires manual optimization of non-adjacent edges" — i.e. no certified
global set. (ii) Kinetic data structures / rotating-rigid-body collision are standard, but the
standard results are for *arbitrary* motions and give amortised bounds, not a geometric predicate
on a *design*; here the motion is a known conic, which is what makes a static certificate possible.
(iii) kiri-R22 IsoGami uses IPC with a spatial hash: a numerical broad phase, no bound.

**Sketch.** Every vertex traces `y(θ) = cos(θ/2)x + sin(θ/2)Jχ`, an **ellipse** with semi-axes
`‖x‖` and `‖χ‖` (in the frame where the face's fixed point is the origin). So the entire
trajectory of face `f` over `θ ∈ [0, π]` lies in a disc of radius `ρ_f := max_{u ∈ f} max(‖x_u‖,
‖χ_u‖)` about a computable centre. `ρ_f` is a **linear** functional of `X` (a max of norms of
linear images), computable once per face. Two faces can contact only if their swept discs overlap:

```
       ‖c_f − c_g‖  ≤  ρ_f + ρ_g .                                            (A4.1)
```

Since the flat tiling is gap-free and planar, `Σ_f area(f)` is the patch area, so the number of
faces with `ρ_f` bounded away from `r_f` by a fixed factor is `O(1)` per unit area — a packing
argument gives `O(n)` overlapping pairs whenever the `ρ_f/r_f` ratios are uniformly bounded. The
theorem to prove is exactly that this ratio is bounded, i.e. `‖χ_u‖ ≤ κ‖x_u‖` with `κ` controlled
by the diameter of Γ along the spanning tree used to build `S` — which is where the argument can
genuinely fail, since `χ = σ_f x_u + 2w_f` and `w_f` is a *path sum* that can grow linearly with
the patch diameter. **If `w_f` grows linearly, the swept discs are global and the theorem is
false** — that is the real content and the real risk.

**Kill experiment.** 500 graphs, 100–5000 faces, three generators. Compute `ρ_f/r_f` for every
face and plot its maximum against patch diameter. Then compute the exact `θ_max` by (i) the
`O(n²)` filtered-root enumeration and (ii) the (A4.1)-pruned enumeration, and check they agree to
1e-12 on every graph. **Two kills:** if `max_f ρ_f/r_f` grows like the patch diameter (fit slope
> 0.3 on a log-log against `n`), the locality theorem is false as stated and the idea reduces to a
heuristic broad phase; if the pruned and unpruned `θ_max` ever disagree, the pruning is unsound.

**If it survives, the demo.** Wall-clock `θ_max` vs `n` for `n ∈ [100, 5000]` on a log-log plot,
pruned vs unpruned vs the reference bisection, with the exact-agreement column; the swept-disc
picture on a hero pattern; the ≥ 20-graph range-optimisation of A1 running at 5000 faces, which is
only possible with the pruning.

**Risk.** `w_f` really is a path sum over a spanning tree of Γ, and for a large patch with a free
boundary the far faces translate a long way. The saving grace is that `w_f` is a sum of
`σ(·)x_{u_·}` with **alternating signs** along the tree, so it may be `O(1)` by cancellation rather
than `O(diameter)` — but that is precisely the unproven step. Genuine 50/50.

**Effort.** Derivation L. Code M.

---

### A5. Which periodic deployment Jacobians are achievable?

**Type:** characterization.

**Claim.** For a periodic pattern, `J(θ) = cos(θ/2)I + sin(θ/2)W` with `W = Q P_0^{-1}` where `Q`
and `P_0` are both **linear** in `X`. Hence the design map `Φ : X ∩ {P_0 invertible} → R^{2×2}`,
`X ↦ W(X)`, is a **rational map of bidegree (1,−1)** whose image is a constructible set of
dimension `≤ min(4, dim X)`. I claim: (i) `Φ` is surjective onto an open cone whose dimension is
computable as the rank of a fixed `4 × dim X` matrix at a generic point; (ii) the finite Poisson
ratio along the deployment is the closed-form rational function
`ν(θ) = −(λ₂ cot(θ/2) + μ₂)/(λ₁ cot(θ/2) + μ₁)` determined by exactly **two invariants of `W`**
(its symmetric part's eigenvalue ratio and its skew part), so the achievable `ν(θ)` curves of a
given periodic graph form a 2-parameter family; and (iii) `ν ≡ −1` for all θ **iff** `W` is
conformal, recovering claim (c) as the degenerate case.

**Gap targeted.** 2026 Sec. 5.1. The paper solves Eq. (13) and reports (`notes/paper_2026.md:820`,
verbatim) "Empirically, we observe that for all tilings with a non-trivial kernel (`E_split != {}`),
solving Eq. (13) yields an embedding with conformal deployment" — an *existence* claim with no
proof and no experiment count, immediately followed by the "remains conformal for all θ" claim.
This idea answers the existence half, which claim (c) does not touch.

**Novelty vs field.** (i) kiri-R1 Grima & Evans give `ν(θ)` in closed form for rotating
squares/rectangles/triangles — a handful of *specific* structures, one formula each, never as the
image of a design map. (ii) kiri-R18 Liu et al. 2024 compute Poisson's ratio in closed form (their
Eq. (12)) "and controlled to some extent" for dihedral Escher quads — one family. (iii) rig-R17
Borcea & Streinu 2015/2018/2020 develop geometric auxetics for periodic frameworks and characterise
*auxetic* one-parameter deformations by a cone condition on the periodicity lattice — the closest
conceptual relative, but for bar-and-joint frameworks and, per STATE, for one specific 3-DOF
framework in the 2020 paper, never as an image-of-a-linear-design-space question. (iv) kiri-R2
Mitschke et al. search an archive for `ν = −1` mechanisms; searching is the alternative to
characterising.

**Sketch.** Let `p_x(θ), p_y(θ)` be the two period vectors, each `cos(θ/2)p° + sin(θ/2)q` with
`p°, q` linear in `X` (geometer (0.8)). Then `P_θ = cos(θ/2)P_0 + sin(θ/2)Q` and
`J(θ) = P_θ P_0^{-1} = cos(θ/2)I + sin(θ/2)W`. The Cauchy–Green tensor is
`J^T J = cos²(θ/2) I + cos(θ/2)sin(θ/2)(W + W^T) + sin²(θ/2) W^T W`; dividing by `cos²(θ/2)` and
writing `τ = tan(θ/2)`, `C(τ) = I + τ(W + W^T) + τ² W^T W` — a **matrix polynomial pencil of degree
2 in τ**, so principal stretches are the square roots of the eigenvalues of `C(τ)`, algebraic of
degree 2 in `τ`. The Poisson ratio in the principal frame is a ratio of logarithmic stretch rates,
hence a rational function of `τ` with coefficients that are polynomial invariants of `W` — of
which there are exactly two (`tr W`, `det W`) up to the orthogonal conjugation that a global
rotation supplies. Conformal `W = aI + bJ_rot` gives `C(τ) = (1 + aτ)²I + …` proportional to the
identity, hence equal stretches, hence `ν = −1` at all θ: claim (c), as a corollary.

**Kill experiment.** Build a library of ≥ 30 periodic graphs (square, triangular, hexagonal,
kagome, (3,4,3,12), the Fig. 11 graph, and ~24 random periodic Voronoi metatiles with 6–20 faces).
For each: solve for `X = ker[L;B]` with periodic boundary; sample 1000 random points of `X`;
compute `W(X)` in closed form and, independently, by finite-differencing the *simulated* period
vectors from `forward_kinematics` at two angles. **Kill 1:** any disagreement > 1e-10 kills the
`W = QP_0^{-1}` closed form. **Kill 2:** measure `dim image(Φ)` as the numerical rank of the
Jacobian `∂vec(W)/∂X` at 100 random points; **if that rank is 0 or 1 on more than half the library**
(i.e. `W` is essentially forced by the graph and not designable), the characterization is vacuous
and the idea dies. **Kill 3:** check `ν(θ)` from the formula against `ν` measured from simulated
period vectors at 10 angles; disagreement > 1e-8 kills (ii).

**If it survives, the demo.** For each of ≥ 20 periodic graphs, plot the achievable region in the
`(tr W, det W)` plane with the conformal locus `tr²  = 4det` marked, and the corresponding family
of `ν(θ)` curves — the first picture of "what Poisson ratio behaviours a given tiling can be
*designed* to have", as opposed to "what one specific structure does". Hero: a graph where the
conformal locus is *empty*, disproving the paper's empirical existence claim, or one where it is
2-dimensional, giving a whole family of isotropic designs the paper's single Eq. (13) solve
returns only one point of.

**Risk.** `dim X` for periodic metatiles is small (2026 Fig. 11's graph: 6; Fig. 13: 8) and
`P_0`-invertibility plus the boundary rows may leave `Φ` with rank ≤ 1, making the "achievable set"
a curve and the characterization thin. Kill 2 is designed to find this fast. Second risk: the
"exactly two invariants" claim needs care about whether the reference frame is fixed by the
lattice or free — if free, `ν(θ)` is direction-dependent and the clean 2-parameter story needs a
direction argument I have not written.

**Effort.** Derivation M. Code M.

---

### A6. Mobility is mostly a boundary artefact: the 2-core theorem

**Type:** theorem.

**Claim.** `dim ker A(Γ) = |F ∖ core₂(Γ)| + dim ker A(core₂(Γ))` exactly, where `core₂(Γ)` is the
2-core of the hinge graph (repeatedly delete nodes of degree ≤ 1). Consequently the mobility
`m = dim ker A − c(Γ)` reported for finite free-boundary patches is dominated by dangling faces,
and **the intrinsic mobility `m_core := dim ker A(core₂) − c(core₂)` is small — I predict 1–3 — for
every planar hinge cut satisfying the connectivity principle of 2026 Sec. 4.2**, on Delaunay and
Voronoi patches alike. Second half: `m_core` evaluated at `θ > 0` (not at the flat state) is the
number that answers 2026 Sec. 6 limitation 2, and it is smaller still.

**Gap targeted.** 2026 Sec. 6 limitation 2 verbatim (`notes/paper_2026.md:955`): "Analyzing the
degrees of freedom of the deployment space … is an interesting direction to explore." And directly:
STATE U4, which reports "19–26 on random Delaunay (finite patch, free boundary)" and calls it the
answer. I claim that number is an artefact and the real answer is an order of magnitude smaller.

**Novelty vs field.** (i) rig-R1 Maxwell/Calladine gives `m − s`, never `m`; the 2-core reduction is
exactly the step that makes the count useful and it is not in the Maxwell tradition. (ii) rig-R7
Jacobs & Hendrickson's pebble game does a combinatorial DOF count but under **genericity**, which
STATE F13 establishes fails here — every working kirigami pattern is in special position and the
pebble game predicts rigidity. (iii) kiri-R22 IsoGami counts zero eigenvalues of `JᵀJ` and filters
for exactly 1 DOF on isohedral tilings; because their tilings are periodic there is no free
boundary and no dangling faces, so IsoGami never encounters this artefact and never states the
reduction. (iv) kiri-R9 Chen, Choi, Mahadevan 2020 derive DOF for specific square-tile
connectivities — named connectivities, fixed geometry, no general reduction.

**Sketch.** A face `f` of degree 1 in Γ appears in **no** cycle of Γ, so its column of `A` is
identically zero (rows of `A` are indexed by a cycle basis, and `A(ω)_z = Σ_{i ∈ z} ω_{f_i}(p_{i−1}
− p_i)`). A zero column contributes exactly 1 to the kernel and nothing else; deleting `f` also
leaves the cycle space of Γ unchanged, so `b₁(Γ) = b₁(Γ ∖ f)` and the remaining rows are unchanged.
Induct. This gives the decomposition as an *identity*, not an approximation — which is why it is a
theorem and not a heuristic. The empirical half (`m_core` small) is the falsifiable part: it is a
claim that the 2-core of a planar hinge graph in the special position enforced by Eq. (2) has
almost no flexes beyond the uniform one, which is *not* implied by anything above.

**Kill experiment.** Rerun exactly the rigidity persona's `rig_check.cpp` configuration — the same
21 graphs, plus 500 more from the three generators — but report four numbers per graph: `m`,
`m_core` at the flat state, `m` at `θ = 0.7`, and `m_core` at `θ = 0.7`. The identity is checked by
`dim ker A(Γ) − dim ker A(core₂) − |F ∖ core₂| = 0` to machine rank tolerance. **Kills:** any graph
where that residual is nonzero kills the theorem (it should not happen; it is a proof). The
*interesting* kill is the prediction: **if `m_core` at `θ = 0.7` exceeds 5 on more than 10% of the
500 graphs, the "intrinsic mobility is small" claim is dead**, and — importantly — that would be
*good news* for the convergent bundle's claim (e), which I am attacking. Either outcome resolves a
live disagreement.

**If it survives, the demo.** A four-column table over 500 graphs and a scatter of `m` vs `m_core`
with the line `m = |F ∖ core₂| + m_core`; a figure showing a Delaunay patch with the dangling faces
shaded, making the artefact visible; and the corrected statement of the answer to Sec. 6
limitation 2. Hero: a patch where `m = 26` and `m_core = 1`, i.e. the entire reported mobility is
twenty-five pendulums.

**Risk.** The 2-core may not remove enough: near-dangling structures (degree-2 chains, whose columns
are not zero but whose cycles are long and thin) may still inflate the count without being
"boundary artefacts" in any clean sense. Then the theorem is right, the prediction is wrong, and
the honest answer is that mobility genuinely is large on irregular patches — which would make the
convergent bundle's (e) stronger than I claim. I am willing to lose this one.

**Effort.** Derivation S. Code S. (This is the cheapest experiment in the file and it should be run
first, because it decides whether anyone should write (e)'s abstract at all.)

---

### A7. The corank of `L` is the number of closed classes of the hole digraph

**Type:** characterization (theorem + validation vs brute force).

**Claim.** `rank(L) = H − dim Z` where `Z = {y ∈ R^H : y_{K(·)} is out-harmonic on the hole digraph}`
(STATE F14), and `dim Z` equals the number of **closed classes** (terminal strongly connected
components) of the hole digraph `D` whose nodes are hole preimages and whose arcs are induced by
hinge edges. Corollary predictions: (i) `dim Z = 0` when some hole reaches the boundary along
reversed hinge edges from every other hole (fixed-boundary case ⇒ `rank L = H`, matching the
paper); (ii) `dim Z ≥ 1` for every boundary-free/periodic pattern, so **2026 Sec. 4.4's
"#independent equations = #holes" is off by exactly the number of closed classes**; (iii)
`dim(shape space) = 2(N − H + dim Z)`.

**Gap targeted.** 2026 Sec. 4.4 verbatim prose: "#independent equations = #holes", stated without
proof, and restated in Sec. 5.1 as the empirical "for all tilings with a non-trivial kernel". This
is the single load-bearing claim under the abstract's "we characterize the full space of
embeddings". STATE F14 already shows it is false for periodic patterns; nobody has the correct
formula. Note this is one of my own "obvious target, non-obvious answer" cases (Part 2 closing
note), and the geometer's G5 is the same target — I include it because it is the claim the paper
most needs and because my formulation (closed classes of a digraph, a Markov-chain object) is
sharper and more testable than "dim ker(S − W)".

**Novelty vs field.** (i) `notes/field_tutte.md` (0.1)/(0.3) gives `leftnull(L) = Z` and the
out-harmonic characterization but stops there — it does not identify `dim Z` combinatorially, and
its `H = V_int − |E_split|` (STATE U1) is missing the `c(Γ)` term. (ii) Tutte/Floater weighted-
average theory needs each vertex to average over its **whole** star; here row `v` says `x_v` is the
mean of the hinge **in**-neighbours, half the star (STATE F14), so no Tutte theorem applies and no
existing rank result transfers. (iii) M-matrix / irreducible-diagonal-dominance theory (the
standard tool) gives nonsingularity under a reachability hypothesis — STATE U2 — but does **not**
give the corank when the hypothesis fails; the closed-class count is exactly the statement that
does.

**Sketch.** Row `K` of `L` is `Σ_{e ∈ K ∩ E_hinge} (x_{head} − x_{tail})`. Left-multiplying by
`y ∈ R^H` gives a weighted signed edge sum which vanishes for all `X` iff `w(e) := y_{K(e)}` is a
**circulation** on the hinge digraph (F14). By Remark A.1 every interior vertex has
`indeg_h = outdeg_h`, so the circulation condition at `v` reads `k_v g(v) = Σ_{v→w} g(w)` — `g` is
harmonic for the **out**-transition matrix, i.e. `g` is a right eigenvector of a row-stochastic
matrix `S` with eigenvalue 1. By Perron–Frobenius theory for stochastic matrices, the eigenvalue-1
right-eigenspace has dimension exactly the number of closed (recurrent) classes of the chain,
with a basis of indicator-like functions supported on each closed class's basin. Fixed boundary
rows make the chain absorbing at the boundary, killing all closed classes among the interior holes
⇒ `dim Z = 0`. A torus has no boundary ⇒ the all-ones vector survives ⇒ `dim Z ≥ 1`.

**Kill experiment.** 500 random graphs (Voronoi, Delaunay, perturbed grid, plus 30 periodic
metatiles), random and GW-assigned σ. For each: build `L` from `hole_preimages`, compute
`rank(L)` by rank-revealing QR with a tolerance swept over 1e-8…1e-12 (to check tolerance
robustness); build the hole digraph and count its closed classes by Tarjan SCC + a "no outgoing
arc" test; compare. **If `H − rank(L) ≠ #closed classes` on any single graph, the characterization
is dead.** This is exactly the ≥500-graph-vs-brute-force validation `specs/ideator.md` asks for.
Also cross-check `dim(shape space)` against the SVD null-space dimension already computed by the
pipeline for 2026 Fig. 11 (expected 6) and Fig. 13 (expected 8).

**If it survives, the demo.** A table of `H`, `rank(L)`, `#closed classes` and the shape-space
dimension over 500 graphs with zero mismatches; the corrected statement of Sec. 4.4; a periodic
example where the paper's formula predicts a shape space of dimension `2(N − H)` and the true one
is `2(N − H + 1)`, i.e. **an extra design degree of freedom the paper's own tool does not expose**,
rendered as a one-parameter family of embeddings.

**Risk.** The rank of `L` may be numerically ambiguous on large graphs (`L` is `H × N` with tiny
integer entries but the graphs are irregular), making the "zero mismatches" claim tolerance-
dependent. Mitigate by computing `rank(L)` over `Q` exactly with a fraction-free elimination on
graphs up to ~500 faces, since `L`'s entries are `±1` integers. Second risk: the geometer's G5 has
the same target, so this is a shared result, not a differentiator.

**Effort.** Derivation M. Code M.

---

### A8. Eq. (2) is sufficient but not necessary: the missing design space

**Type:** theorem (existence) + counterexample search.

**Claim.** There exist planar graphs `M` with orientation σ such that (i) the uniform system
`LX = 0` with a nondegenerate boundary admits **only degenerate** solutions (self-overlapping or
zero-area faces), yet (ii) the cut structure `M'` admits a **finite non-uniform** rigid deployment
of positive extent from the flat state. That is, **the class of deployable kirigami is strictly
larger than the class Segall 2026 characterises**, and the gap is exactly the set of designs whose
hinge-angle vector leaves the diagonal `θ·1`.

**Gap targeted.** The framing of 2026 Sec. 3 itself: uniform deployability (`θ_e = θ` for all `e`)
is declared as the object of study and 2026 Sec. 6 limitation 2 concedes that "patterns that can be
uniformly deployed might also be deployed with a different set of hinge angles". Nobody asks the
converse — patterns that **cannot** be uniformly deployed but can be deployed. If such patterns
exist and are common, the paper's abstract claim to "characterize the full space of embeddings that
admit uniform deployment" is true but the *interesting* space is larger, and the whole linear-system
framework is a restriction whose cost has never been measured.

**Novelty vs field.** (i) kiri-R12 Dang et al. 2021's theorem (per kiri-R24 Sec. 2.2, since I could
**not** retrieve the paper myself — see Part 0 F1) states a per-cut-loop compatibility condition
`g(cos β₁) = cos β₁`, with the note "When all the cuts are parallelograms… compatibility holds for
any β₁" — i.e. their condition is already non-uniform-angle in spirit, but for **quads only** and
as a loop condition, not as a comparison of two design classes. (ii) kiri-R17's "deployment angle
field" is again quads and is a *construction*, not an existence/non-existence statement. (iii)
rig-R14 Kapovich & Millson 2002 (planar linkage moduli can be arbitrary varieties) is the correct
warning that the non-uniform configuration space can be wild — which is why the claim must be
*existence of a positive-extent path*, not a full characterization.

**Sketch.** The non-uniform deployment variety is cut out by the cycle closure conditions on Γ with
a per-hinge angle `θ_e`. Around a cycle `z` of Γ with faces `f_1 … f_k`, composing the rigid motions
gives `Π_i R(ω_i) = I` and a translation closure; the rotation part forces `Σ_i ±θ_{e_i} = 0 mod 2π`
per cycle (a **linear** condition on the angle vector, independent of `X`), and the translation part
gives two **trigonometric** conditions per cycle. The uniform ansatz `θ_e ≡ θ` satisfies the
rotation part automatically (alternating signs, even cycle length by bipartiteness of Γ) and reduces
the translation part to `L X = 0`. **Any other integer solution of the rotation conditions gives a
different ansatz `θ_e = n_e θ` with `n_e ∈ Z`, and a different linear system `L_n X = 0`.** So the
first concrete question is whether `ker L_n` can be nondegenerate while `ker L_1` is not — a finite
search over small integer vectors `n` in the kernel of the cycle-sign matrix. This is a clean,
computable, non-obvious question and it may well have a yes.

**Kill experiment.** Build the cycle-sign matrix `Σ ∈ {0,±1}^{b₁(Γ) × |E_hinge|}` and enumerate
integer vectors `n ∈ ker Σ ∩ {−3…3}^{|E_hinge|}` for small graphs (`|E_hinge| ≤ 24`, so use the
hexagonal, kagome, (3,4,3,12), Fig. 11 and ~200 random 12–20-face patches). For each `n`, assemble
`L_n` and compute `dim ker[L_n; B]` and whether a random point of it is an embedding (no inverted
faces, checked by `signed_area`). **Kills the idea if, on all ~200 graphs, every `n ≠ ±1` gives a
kernel that is either empty of embeddings or a subset of `ker L_1`'s embeddings.** Additionally the
existence half must be certified dynamically: for any candidate, run `forward_kinematics` with the
non-uniform angle vector `nθ` and verify face rigidity to 1e-12 over `θ ∈ (0, 0.5]`.

**If it survives, the demo.** A graph that is *not* uniformly deployable — for which 2026's Eq. (6)
"repair" would move the vertices — but which deploys rigidly with angle ratios `(1, 2, 1, 2, …)`,
shown side by side with Segall's repaired version, with the repaired version's geometric distortion
quantified. Fabricable export of the non-uniform pattern. This is the only idea in my file that
would make Segall et al. change their framing rather than extend it.

**Risk.** The integer-ratio ansatz `θ_e = n_e θ` may be the *only* accessible non-uniform family
(the general variety is not a linear space and enumerating it is hopeless), and it may be that every
`n ≠ ±1` produces a degenerate `L_n` for parity reasons — the rotation-closure condition on an even
cycle with alternating signs is quite rigid. Then the answer is "no", which is still worth one
paragraph but not a paper. Second risk: kiri-R12's theorem may already contain the quad case of
this, and I could not verify its statement (Part 0, F1). **This idea must not be committed to
before someone reads Dang, Feng, Duan & Wang 2021 in full.**

**Effort.** Derivation L. Code M.

---

### A9. How many analytic branches pass through the flat state?

**Type:** theorem.

**Claim.** The flat configuration `Y_0` is a singular point of the configuration variety (STATE U5:
mobility drops from `m(0)` to a smaller constant the instant `θ > 0`, on 20 of 21 graphs). I claim
the number of distinct analytic branches through `Y_0` is computable exactly, because **one branch
is known in closed form** — the uniform path `Y_θ` — so the others can be found by deflation: write
a candidate motion as `Y_0 + ε v + ε²w + …`, impose the hinge constraints order by order, and
quotient by the known branch. Specifically: the second-order (prestress) condition
`v^T H_λ v = 0 ∀λ ∈ coker(rigidity matrix)` selects which of the `m(0)` first-order flat flexes
extend, and I claim **the number that extend equals `m(θ>0) + 1`, the drop measured in U5 being
exactly the number of first-order flexes killed at second order.**

**Gap targeted.** 2026 Sec. 6 limitation 2, again, but from the branch side rather than the DOF
side. Also, sharply: **2026 Eq. (9) is minimised at exactly this singular configuration** (STATE
U5's closing note). If the flat state is a branch point, a first-order energy evaluated there is
being evaluated at the one configuration where first-order information is least reliable — which
is the mechanism-theoretic sibling of my A1.

**Novelty vs field.** (i) rig-R8 Connelly & Whiteley 1996 second-order rigidity and prestress
stability is precisely the right tool and is *not applied* anywhere in the kirigami literature —
STATE F13 records that no published proof exists that even rotating squares is exactly 1-DOF.
(ii) rig-R9 Connelly & Servatius 1994 supply the cusp warning: a configuration can be second-order
rigid yet flexible, so the count needs care and the theorem must be stated as an upper bound plus a
certificate. (iii) rig-R25 Li, Zhu & Qu's 4th-order rigid-origami mobility framework is the closest
published finite-DOF pipeline — for origami, numerically, without a known closed-form branch to
deflate against. **Having one branch in closed form is what makes the deflation tractable and is
unavailable to every prior work.**

**Sketch.** Let `R(Y)` be the pin-constraint map, `R(Y_θ) = 0` along the known branch. First-order
flexes are `ker DR(Y_0)`, of dimension `m(0) + 3`. Differentiate the closed form twice:
`dY/dθ|_0 = ½ J S X` and `d²Y/dθ²|_0 = −¼ X'`, giving the known branch's 2-jet for free. For any
other candidate `v ∈ ker DR(Y_0)`, second order requires `D²R(Y_0)[v,v] ∈ image(DR(Y_0))`,
equivalently `λ^T D²R(Y_0)[v,v] = 0` for every self-stress `λ`. That is a **quadratic form on
`ker DR`**, i.e. a quadric in `P(ker DR) ≅ P^{m(0)+2}`; its real points, modulo the trivial motions
and the known branch, bound the branch count. The claim `#extending = m(θ>0) + 1` is the sharp
prediction and is exactly the thing that could be false.

**Kill experiment.** On the 21-graph set plus 200 random graphs: compute `m(0)` and `m(0.7)` (with
the A6 2-core correction applied); compute the self-stress space `coker DR(Y_0)`; assemble the
quadratic forms `λ^T D²R[v,v]` on `ker DR(Y_0)`; compute the dimension of the real variety of their
common zeros by random sampling plus rank of the Jacobian at sampled points. **Kills the idea if
the predicted `#extending` differs from the measured `m(0.7) + 1` on any graph**, or if the common
zero set is not a manifold near the known branch (a cusp, rig-R9's warning) on more than a couple of
graphs — in which case the honest statement is an inequality, not an equality, and the idea
downgrades to a bound.

**If it survives, the demo.** For ≥ 20 graphs, a table of `m(0)`, `#second-order survivors`,
`m(0.7)`, showing the identity; the triangle tiling — STATE U5's unexplained exception where
mobility does *not* drop — explained as the case where the quadric is identically zero. Hero: a
pattern with a *second* branch that is not the uniform one, deployed both ways from the same flat
sheet, fabricated.

**Risk.** The quadric's real variety can be computed reliably only in low dimension; `m(0)` reached
26 in U4, and a real variety in `P^28` cut by many quadrics is not something to be sampled
casually. A6 must land first (reducing `m(0)` to the core value) or this idea is not computable.
Explicit dependency: **A9 requires A6.** Second risk: rig-R9's cusps could make the branch count
ill-posed on a positive fraction of graphs.

**Effort.** Derivation L. Code L.

---

### A10. The shape space is bigger than the design space: emptiness of the usable region

**Type:** characterization (with an emptiness theorem as the sharp form).

**Claim.** The advertised design space is the linear space `X = ker[L; B]`, of dimension
`2(N − rank L)` (2026 Eq. (5), abstract: "we characterize the full space of embeddings that admit
uniform deployment"). The **usable** design space is
`U(ε) = {X ∈ X : X is an embedding (all faces positively oriented) and θ_max(X) ≥ ε}`, a
semialgebraic subset. I claim (i) `dim U` can be **strictly less** than `dim X` — because the
embedding constraint is open but the `θ_max ≥ ε` constraint is not, and the exact `θ_max` of A1/A4
is a min of algebraic functions with a nonempty active set on a set of positive codimension; and
(ii) **there exist planar graphs and σ for which `U(ε)` is empty for every `ε` bounded away from 0
while `dim X > 0`** — i.e. patterns the paper's tool reports as "uniformly deployable with a
6-dimensional shape space" that have no usable design at all.

**Gap targeted.** The abstract's word "**full**", and 2026 Sec. 6 limitation 1 verbatim: "although
all embeddings in the solution space `X` are theoretically uniformly deployable, some may be
geometrically undesirable (e.g., exhibiting extremely short edges or sharp angles). Moreover,
uniform deployability does not guarantee a large collision-free deployment range." The paper states
the phenomenon and offers "post-processing". Nobody has measured how much of `X` survives, or
whether `U` can be empty.

**Novelty vs field.** (i) kiri-R2 Mitschke et al. 2013 search an archive for auxetic mechanisms and
report which tessellations work — a *discrete* yes/no over a catalogue, no notion of a design space
with a usable subregion. (ii) kiri-R18 Liu et al. 2024 explicitly report "a **smaller solution
space** compared to traditional Escher dihedral tessellations" as a limitation of their
deployability constraints — the same phenomenon, observed once, for one family, never quantified.
(iii) kiri-R22 IsoGami filters its catalogue by simulation with collisions, i.e. it *empirically*
discards unusable designs, which is precisely the procedure a characterization would replace. No
one has a volume, a dimension, or an emptiness criterion.

**Sketch.** `X` is linear. The embedding condition is `sign(area(f)) > 0` for all `f` — an open
polyhedral-ish condition (each `area(f)` is a quadratic form in `X`). `θ_max(X) ≥ ε` is, by A1/A4,
`min_i t*_i(X) ≥ tan(ε/2)` where each `t*_i = −R_i(X)/P_i(X)` is a **ratio of quadratics in `X`**
(geometer (0.8): `P, Q, R` are quadratic in `X`). So

```
   U(ε) = {X : area_f(X) > 0 ∀f}  ∩  {X : R_i(X) + tan(ε/2)·P_i(X) ≤ 0  ∀ active i} ,
```

a basic semialgebraic set defined by **quadrics**. Its dimension is computable numerically as the
rank of the active constraint gradients at a boundary point; its emptiness is a quadratic
feasibility problem, decidable in practice by an SDP relaxation (Shor / Lasserre level 1) whose
*infeasibility certificate* is a rigorous proof of emptiness. That is the theorem-grade half:
**an SDP certificate that a given graph + σ admits no usable design.**

**Kill experiment.** 500 graphs, three generators. For each: compute `dim X`; sample 10⁴ points of
`X` from a Gaussian in the null-space basis (this is exactly what 2026's UI sliders expose); record
the fraction that are embeddings and the distribution of exact `θ_max`. Then run the Shor SDP
relaxation for `U(π/6) = ∅` on every graph with a positive sampling failure rate. **Kills the idea
if (a) the embedding+range failure rate is below 5% on essentially all 500 graphs (then `X` is
usable and the distinction is academic), or (b) the SDP certificate never fires on any graph with
`dim X > 0` and no sampled usable point (then emptiness is real but uncertifiable, and the sharp
form dies while the measurement survives).**

**If it survives, the demo.** For ≥ 20 graphs, a "usable fraction" bar chart against `dim X`,
showing that shape-space dimension is a poor proxy for design freedom; the hero is a graph the
paper's pipeline accepts, whose sliders produce nothing but self-intersecting or immediately
colliding patterns, together with the SDP certificate proving no slider setting can work. Then the
constructive companion: project onto `U` by minimising the A1 exact-range objective from a
feasible start, i.e. A1's algorithm used as the repair the paper leaves to "post-processing".

**Risk.** The dominant risk is (a): Segall's Eq. (6) least-norm solve starts from a sensible
`X_ini`, so the neighbourhood of the returned `X_0` is probably fine, and the failure rate may only
be large far out along the sliders — where no user goes. Then the result is "the shape space is
locally usable and globally not", which is true, mildly interesting, and not a paper on its own. It
would still be the right §7 of the A1 paper.

**Effort.** Derivation M. Code L (an SDP solver is not in the C++ stack; Shor's relaxation for
these sizes can be done with an Eigen-based projected-gradient / spectral bundle, which is real
work — this is the one idea in my file with a genuine implementation risk under the C++-only
directive D3).

---

## Part 4 — Self-attack: the strongest rejection I can write against my own top three

### Against A1 (Eq. (9)'s false-negative rate) — **Reject, borderline.**

> The paper's contribution is that a published collision heuristic is a first-order test and that
> first-order tests miss second-order events. This is not a discovery; it is the definition of a
> first-order test, and the original authors plainly know it — they say so themselves in
> Section 6, in the sentence the submission quotes. The submission's algebra reduces to observing
> that a function `s(θ) = P + Q cos θ + R sin θ` with `s(0) = 0` has one further root at
> `t = −R/P` and that the published energy constrains only `R`. That is one line of trigonometry.
>
> Worse, the empirical claim is not against the method as implemented. Equation (9) is stated in
> the paper with an *unspecified* barrier `B`, an *unspecified* weight `γ`, and an unspecified
> diameter-selection rule; the authors ship a WebAssembly implementation at a public repository.
> A false-negative rate measured against the submission's own reconstruction of an underspecified
> energy measures the submission's reading, not the method. Until the authors' code is run, the
> headline number is not evidence about anything.
>
> Finally, the proposed remedy — penalise `min(1, −P/R)` — is a one-line change to an existing
> energy. Even at a 40% false-negative rate, this is a fix to somebody else's Section 4.5. It is a
> good bug report and a poor paper.

**My honest reply, as the author.** The reviewer is right that the mathematics is one line and
right that the number must be measured against `github.com/segaviv/tuttekiri`, not against my
reconstruction; that check moves to step zero of the project. The reviewer is wrong that a
quantified systematic failure of the state of the art is worth nothing — but only if the number is
large *and* the exact-range optimiser of the same paper converts it into range gained on ≥ 20
graphs. **A1 is not a paper by itself. It is the §5 of the A4+A1 paper.** If A4's locality theorem
fails, A1 should be a two-page note, not a submission.

### Against A2 (the bifurcation pencil) — **Reject.**

> The submission proves two things. The first is that the uniform deployment path never locks,
> because the uniform velocity field lies in the kernel of the mobility operator at every
> configuration on the path. This is immediate: the path exists, therefore it has a tangent,
> therefore the tangent is in the kernel. The submission presents this as a lemma.
>
> The second is that the mobility matrix along the path is a linear pencil `A_c + t A_s` in
> `t = tan(θ/2)`, so rank drops occur at the roots of the pencil's determinantal ideal. Given the
> closed-form deployment map — which the submission concedes is a restatement of the definition of
> uniform deployability, standard in the rotating-rigid-units literature — this too is immediate:
> a matrix whose entries are affine in `t` is an affine pencil.
>
> What would make this a paper is the third claim, that bifurcations actually occur inside the
> collision-free range and change what the physical structure does. Here the submission offers a
> 500-graph count and no mechanics. A rank drop of the first-order mobility operator is not a
> bifurcation of the configuration variety; it is a necessary condition for one. Connelly and
> Servatius showed thirty years ago that first-order degeneracy can be an artefact of the
> linearisation. Without a second-order analysis at each candidate `θ_b` — which the submission
> defers to a companion idea — the reported "bifurcation spectrum" is a spectrum of *candidates*.
> And with no energy, no actuation model and no stability argument, the claim that a physical
> sheet leaves the uniform branch is unsupported by anything in the paper.

**My honest reply.** Devastating and correct on the third point. A2 must ship with the
second-order test at each candidate `θ_b` — which is A9 restricted to `θ_b` instead of to `θ = 0`,
a much smaller computation because `m(θ_b)` is small. **A2 and A9 are one idea, not two**, and the
combined idea is "the exact branch structure of the uniform path". The no-locking lemma stays a
lemma. If the 500-graph sweep returns an empty bifurcation set (outcome α), the whole thing
collapses to two remarks and must be abandoned within the first afternoon.

### Against A6 (the 2-core theorem) — **Reject, decisively.**

> The submission observes that a rigid body attached to a structure by a single pin joint can
> rotate about that pin, and that such bodies therefore contribute one degree of freedom each. It
> formalises this as a decomposition of the kernel of the mobility operator over the 2-core of the
> hinge graph, and verifies it numerically on 521 graphs.
>
> This is correct and it is a textbook exercise. The 2-core reduction is the first step of every
> rigidity-percolation and pebble-game implementation in the literature; Jacobs and Hendrickson's
> algorithm performs it implicitly. That the reduction has not previously been written down *for
> kirigami* reflects that no previous kirigami paper reported a raw mobility count on a
> free-boundary finite patch, because everyone else works with periodic or fixed-boundary
> patterns where dangling tiles do not arise.
>
> The submission's real content is therefore the empirical claim that the reduced mobility is
> small. But "small" is measured on three synthetic graph generators, and no mechanism is offered
> for why a planar hinge cut in the special position enforced by Equation (2) should have a nearly
> trivial 2-core mobility. An unexplained empirical regularity on synthetic data is a figure.

**My honest reply.** Fair, and I do not have the mechanism. A6 is not a contribution; **it is a
correction I owe the rest of the project before anybody publishes a mobility number**, and its
right home is a paragraph in someone else's §7 plus a line in `STATE.md` retracting U4's headline.
I rank it high for *project value* and low for *paper value*, and those are different axes. Any
submission that leads with A6 deserves the rejection above.

---

## Part 5 — Final ranked list

Ranked by *expected value to this project*, which is (probability the kill test passes) × (value if
it does) ÷ (cost). Where that disagrees with paper-worthiness, I say so.

| # | Idea | One-line justification |
|---|---|---|
| 1 | **A1 — Eq. (9)'s false-negative rate** | The only claim in this file whose answer is genuinely unknown, cheap to get, and directly indicts the state of the art with an exact drop-in fix. Not a paper alone; the core of one with A4. |
| 2 | **A2 (+A9) — no-locking lemma and the bifurcation pencil** | Highest ceiling of anything I can defend: nobody, including the concurrent IsoGami and PyKirigami, knows whether kirigami deployments bifurcate. Genuinely 50/50, resolvable in an afternoon, and either answer is publishable *if* the second-order test ships with it. |
| 3 | **A6 — the 2-core theorem** | Lowest cost in the file and it decides whether the convergent bundle's claim (e) survives at all. Run it first regardless of what gets written. Low paper value, high project value. |
| 4 | **A4 — the certified O(n) active contact set** | The one place a real theorem is required rather than convenient, with a real chance of being false (`w_f` may grow with patch diameter). It is also what makes A1 scale to 5,000 faces, which is what `specs/ideator.md` demands. |
| 5 | **A5 — achievable `W` and the Poisson-ratio family** | The right question hiding behind convergent claim (c): not "is conformality preserved" but "which deployment Jacobians can a given tiling be designed to have". Speaks directly to the auxetics community (kiri-R1, kiri-R2, rig-R17) that the CG papers have been talking past. |
| 6 | **A8 — Eq. (2) is sufficient but not necessary** | The most exciting idea here and the only one that would change Segall's *framing* rather than extend it. Ranked sixth only because it is gated on reading Dang, Feng, Duan & Wang 2021 in full, which I could not retrieve (Part 0, F1). **Promote to rank 2 if kiri-R12 turns out not to contain it.** |
| 7 | **A7 — corank of `L` = number of closed classes** | The paper's load-bearing unproved claim, with a sharp Perron–Frobenius answer and the cleanest 500-graph-vs-brute-force validation available. Ranked here only because the geometer's G5 targets the same result, so it is shared, not differentiating. |
| 8 | **A3 — fully closed ⟺ root coincidence** | Correct reformulation of a surrogate objective and a probable disproof of Eq. (14)'s sufficiency, but it improves one short subsection (2026 Sec. 5.2) and the counterexample may only exist on irregular graphs. |
| 9 | **A10 — emptiness of the usable design region** | Attacks the word "full" in the abstract, and the SDP infeasibility certificate is the only theorem-grade emptiness result available. Ranked ninth because the most likely outcome is "locally usable, globally not", and because an SDP in pure C++ (directive D3) is the largest implementation risk in the file. |
| 10 | **A9 — branch count at the flat state** | The deepest question and the least tractable: it depends on A6 landing first, it lives in `P^{m(0)+2}`, and rig-R9's cusps may make the count ill-posed. Its useful fragment is the second-order test at `θ_b`, which I have already folded into A2. |

**If I were forced to commit to one project today:** A6 in the morning (two hours, decides whether
U4/claim (e) is real), then A4 → A1 as the paper, with A2 run in parallel as a cheap lottery ticket
whose downside is one afternoon and whose upside is the more interesting of the two papers.

---

## Part 6 — What I could not verify

- **kiri-R12 (Dang, Feng, Duan, Wang 2021, Phys. Rev. E 104, 055006).** My fetch of
  `https://arxiv.org/pdf/2106.15891` returned unextractable text (Part 0, F1). Every statement I
  make about its theorem comes from the kiri-R12 row, which itself derives from the Dang & Paulino
  2025 survey. **A8 must not be committed to before someone reads this paper in full.**
- **The authors' implementation.** 2026 Eq. (9)'s barrier `B`, its weight γ and the diameter rule
  are absent from the paper (Reader's flags in STATE). A1's headline number is only meaningful
  against `github.com/segaviv/tuttekiri`, which I did not retrieve.
- **kiri-R22 IsoGami.** The kiri-R22 row is explicitly marked "[search-result summaries only; ACM
  full text returned HTTP 403]". STATE F12 says the orchestrator later read the PDF in full and
  confirms no theorem / no closed form / no rank formula / isohedral only. I relied on F12, not on
  my own reading of `papers/related/isogami.txt`, which I skimmed only for the collision and
  mobility sections.
- **Every numeric claim attributed to another persona** (U3, U4, U5, geometer (0.6)/(0.9)) is
  their measurement, not mine. I ran no code in this session; my contribution is adversarial
  analysis and the design of kill tests. In particular my central prediction — that `m_core` is
  1–3 where `m` is 19–26 — is a **prediction**, not a result, and A6 exists to test it.
- **The rotating-rigid-units priority claim.** Q1 and Q7 returned abstracts and summaries, not the
  full texts of Grima & Evans 2000/2005. I am confident the ±θ/2 alternating rotation is folklore
  there, but I did not read those papers, and the exact form in which they state it (infinitesimal
  flex vs finite motion; specific unit shapes vs arbitrary polygons) matters for how sharply a
  novelty claim can be phrased. Someone must read them before Lemma 1 is written.

*End of file.*
