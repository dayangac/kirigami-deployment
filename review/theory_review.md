# Hostile review — theory bundle (T1–T7, K7, F35)

Reviewer stance: TOG/rigidity-theory referee. No praise offered; only what a referee would write in a rejection
letter, and the honest counter.

---

## T1 — Trig-linear deployment (`Y_θ = cos(θ/2)C + sin(θ/2)S`)

1. **Statement.** Uniform deployment is closed-form linear in `(cos θ/2, sin θ/2)`, with `C, S` linear in
   the flat vertex positions `X`.
2. **Classification: TRIVIAL-CONSEQUENCE-OF-KNOWN.** This is the 2D specialization of the Tay–Whiteley
   body-and-hinge motion-assignment correspondence (scalars on hinges, `Σ ω D = 0` per cycle) — `notes/screen_bundle.md`
   already concedes as much for the closely-related T2 mobility matrix (S6), and T1's derivation is the same
   linear-algebra move one step earlier: solving the hinge-point-matching condition for a potential on `Γ`.
   The half-angle magnitude itself is anticipated (with two extra hypotheses) by **Acuña, Gutiérrez, Silva,
   Palza, Núñez, Düring, *Auxetic behavior on demand*, Commun. Phys. 5 (2022), arXiv:2101.12352** — counter-
   rotation on arbitrary bipartite planar polygon networks.
3. **Dismissal sentence.** "This is Tay–Whiteley's motion assignment plus the observation that a rigid-body
   isometry in the plane linearizes trigonometrically in the half-angle — undergraduate linear algebra once
   the correspondence is known."
4. **Rebuttal.** The bundle does not claim T1 as a novel *fact*; it claims the derivation *from the paper's own
   hole-preimage formalism*, tying Eq. (2) to the graph closure condition with an explicit, checked sign
   convention (`derivations/core.md:190-390`, Check C2 vs `kinematics.jl::deploy()`, ⟨JULIA:check_t1_t2:1.42e−14⟩). Acuña et al.
   require collinearity and a fixed 2-colouring; T1 removes both and handles non-bipartite duals via split
   cuts — a real generalization, just not one that needs new mathematics.
5. **To bullet-proof.** State explicitly, in the paper, that T1 is "the Tay–Whiteley motion assignment made
   trig-linear and specialized to corner-hinged tilings", cite Tay–Whiteley (1984) directly (the bundle's own
   citation for it is second-hand — `notes/screen_bundle.md:306` flags the Handbook excerpts as truncated and
   unverified against the primary source), and cite Acuña et al. for the closest published half-angle result.

---

## T2 — No-locking (`σ ∈ ker A(Y_θ)` for all `θ`)

1. **Statement.** The uniform-deployment velocity field is an infinitesimal flex of the pin-jointed structure
   at every `θ`, so the branch never hits a kinematic dead centre; the only stopping mechanism is contact.
2. **Classification: TRIVIAL-CONSEQUENCE-OF-KNOWN.** Exhibiting an explicit tangent vector to an explicit
   curve of realizations and checking it against the rigidity matrix is the standard "first-order flex along a
   known finite motion" argument; it follows immediately once T1's closed form exists. No literature claims
   otherwise for this specific mechanism, but the *proof technique* is textbook (Tay–Whiteley / classical
   infinitesimal rigidity).
3. **Dismissal sentence.** "Once you have written down an explicit one-parameter family of embeddings, of
   course its velocity lies in the kernel of the linearized constraint operator along the family — that's the
   definition of a smooth path in the constraint variety, not a theorem."
4. **Rebuttal.** The content is not "a flex exists" but the *corollary*: mobility `m ≥ 1` at every `θ ∈ R`
   with no kinematic termination, which is what licenses T4's claim that the candidate-contact list is
   *complete* (T2.4, `derivations/core.md:452-460`). That corollary is doing real work downstream and is not
   itself in any screened source.
5. **To bullet-proof.** T2.H.2 already concedes the real gap: `σ ∈ ker A` is a first-order statement about one
   exhibited curve, not a proof the configuration space is a smooth 1-manifold there, and other branches may
   cross it (`dim ker A(Y_θ)` may exceed 1). A referee will ask for either a transversality argument ruling out
   coincident branches, or an explicit disclaimer that only the *uniform* branch's non-degeneracy is claimed.

---

## T3 — Harmonic predicates (`h(θ) = p + q cos θ + r sin θ`, coefficients quadratic in `X`)

1. **Statement.** Every orientation/dot-product/squared-length predicate on the deployed vertices is a
   single first harmonic in `θ`, coefficients quadratic forms in the flat coordinates.
2. **Classification: KNOWN (form) / NEW (the specific closed-form coefficients and the systematic-degeneracy
   catalogue).** The half-angle substitution `c² = (1+cos θ)/2` etc. producing a harmonic from a bilinear form
   in `(c,s)` is Weierstrass/Chebyshev-level algebra — nothing here is a discovery. `notes/screen_bundle.md`'s
   S3 row records the closest prior *use* of this fact: **Grima, Chetcuti, Manicaro, Attard, Camilleri, Gatt,
   Evans, "On the auxetic properties of generic rotating rigid triangles", Proc. R. Soc. A 468 (2012) 810–830**
   — read in full — give a closed-form "locking angle" for the hinge-adjacent two-triangle case, conceding
   "this will not be possible in the general case".
3. **Dismissal sentence.** "A quadratic form composed with the Weierstrass substitution is a harmonic —
   this is not a predicate theorem, it's trigonometric identity-pushing dressed as Proposition T3.2."
4. **Rebuttal.** The load-bearing fact is explicitly *not* the harmonic form (T3.2, `derivations/core.md:520`
   says so itself: "not the harmonic form itself... This is the load-bearing fact for T5") but that the
   coefficients are quadratic in the *design* coordinates `t`, which is what turns contact detection into a
   closed-form root-find and the certificate/gradient machinery (T5, T6) into algebra rather than numerics.
   Grima et al. solve one degenerate sub-case (hinge-adjacent, collinearity assumed); T3 is uniform over all
   face pairs and all planar graphs, with the phase term (`r sin θ`) that non-adjacency forces and that Grima's
   case doesn't need.
5. **To bullet-proof.** The systematic-degeneracy catalogue (T3.H.5, `p+q=0` at every permanent incidence and
   split-edge pair, 73 446/2 145 387 measured occurrences; the round-3 sub-class `p=−q, r=0`) is genuinely new
   bookkeeping but currently lives only as a derivation note plus a measured count. It needs to be stated as a
   lemma ("`h(0)=0` and `h'(0)=0` occur iff …") with a proof, not a measurement, or a referee will ask why a
   supposedly exact predicate needs an empirically-discovered deflation rule.

---

## T4 — Exact `Θ_max` (graze handling, moving-centroid broad phase)

1. **Statement.** `Θ_max` is computable exactly via a complete, closed-form candidate list of contact angles
   (T4.2″), correctly handling grazing (vertex-vertex) contacts that are *not* overlaps; the sound broad-phase
   pruning must use the moving centroid, not the flat one (T4.5a/b).
2. **Classification: NEW for the exactness/completeness machinery; TRIVIAL-CONSEQUENCE-OF-KNOWN for the
   underlying geometric fact.** Lemma T4.2 ("if two simple polygons' interiors are disjoint but they touch, a
   vertex of one meets the boundary of the other") is a standard computational-geometry fact — the honest
   citation is general polygon-intersection theory (a broad-phase/narrow-phase argument any CG textbook has),
   not a discovery. The 2026 paper's own `β_e = 2π − α_f − α_g` bound (T4.4) is **KNOWN** — the bundle itself
   shows `β_e` is just one member of the candidate root list, and that the paper's claimed equality
   `Θ_max = min_e β_e` is **false on the paper's own triangle example** (`Θ_max = π` vs `min β = 4π/3`).
3. **Dismissal sentence.** "`Θ_max` as a min over closed-form roots is a candidate-event sweep, the same
   pattern every collision-detection broad/narrow-phase system uses; nothing about polygon simplicity or
   grazing contact is new to computational geometry."
4. **Rebuttal.** Two contributions survive scrutiny that are not in the screened literature: (a) the
   *counterexample-driven correction itself* — round 1 asserted `Θ_max = min_e β_e` on split-free patterns,
   which is refuted by the bundle's own hexagon case (a genuine grazing contact at `θ=π/3` that is not an
   overlap, `derivations/core.md:665-693`) and by the triangle row of the T4.2 measurement table; the
   *direction* that survives (`≤`, not `=`) is proved, the reverse only measured on 8 patterns and explicitly
   **not proved** (T4.4, `derivations/core.md:785-800`) — this honesty is itself worth stating as the paper's
   claim, since it directly falsifies 2025/2026's implicit assumption; (b) Lemma T4.5b′, that the swept radius
   about the *moving* centroid equals the flat circumradius exactly (not merely bounded), which shows the
   published K2c locality gate is *identically* `1` and therefore vacuous — a real bug-in-the-test finding, not
   present in Konaković-Luković or Grima.
5. **To bullet-proof.** T4.3 (genericity: first contact is vertex-into-edge-interior) is stated as a
   **CONJECTURE** and the file itself shows it fails on the paper's own hexagon figure — "any implementation
   that assumes T4.3 will be wrong on the paper's own figures" (`derivations/core.md:717-724`). That needs
   either a real proof restricted to a stated symmetry class, or it must be dropped from any claim of generality
   and flagged as a per-pattern check. The `≥` direction of `Θ_max ≤ min(min β, π)` needs either a proof or
   must never be printed as equality without the `[N]` tag and sample size.

---

## T5 — The certificate (sound, not complete)

1. **Statement.** `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT` at one probed angle plus the deflated three-class root list
   proves `Θ_max(X̂) ≥ ε` (T5.2b′, T6.4 item 2); it is a sound inner approximation `R(ε) ⊆ U(ε)` of the true
   usable region, not the region itself.
2. **Classification: KNOWN (the general pattern) / NEW (the explicit inner region for this problem).**
   That configuration spaces of linkages are semialgebraic is textbook (Kapovich–Millson, King — the bundle
   concedes this itself at T5.2a: "the content here is the explicit description, not the fact"). Certifying a
   nonconvex feasible region by a sound-but-incomplete inner witness is standard practice in verified
   optimization / SOS relaxations generally, and in the linkage-cone literature specifically: **Rote, Santos,
   Streinu, *Expansive Motions and the Polytope of Pointed Pseudo-Triangulations*, 2003** already gives an
   exact polyhedral (hence exactly-describable, sound) cone of first-order expansive velocities for pointed
   pseudo-triangulations — the closest prior "certificate" for a cognate problem.
3. **Dismissal sentence.** "A sound-but-incomplete semialgebraic certificate for a nonconvex feasible set is
   exactly what you'd expect from Tarski–Seidenberg plus a convex relaxation of the interval constraints — this
   is off-the-shelf verified-computation methodology applied to a new domain, not a theory contribution."
4. **Rebuttal.** The round-2→round-3 history is itself the strongest asset here, not a weakness: round 1's
   candidate certificate was neither sound nor complete (D4: not proved to be an inner approximation at all;
   D5: the atom list was wrong-degree). What's now proved (T5.2b′) is a genuine containment theorem specific to
   the harmonic-predicate structure of this problem (deflating systematic `g(0)=0` roots, T3.H.5), and the five-
   item T6.4 certificate is exercised end-to-end on 33 K7 patterns with agreement to the tolerance the paper's
   own deviation-9 convention sets. Rote–Santos–Streinu's cone is for *unconstrained* expansive motions of
   pseudo-triangulations with no cut structure; nothing there handles split-edge parallelism or a hinge-graph
   with holes.
5. **To bullet-proof.** State plainly, as the paper's claim, that the certificate is a **sufficient, not
   necessary**, condition — round 2's language ("no face–face interior overlap at θ=0") is explicitly withdrawn
   as *wrong*, not just imprecise (T6.4 item 2, `derivations/core.md:1802-1810`), so any earlier paper draft
   using that phrasing must be corrected. Report completeness gap quantitatively: how much of `U(ε)` does `R(ε)`
   miss, at least on the corpus, or concede it is unknown.

---

## T6 — Exact gradients (implicit differentiation of the root map)

1. **Statement.** `∂θ_π/∂t_i` is closed-form via implicit differentiation of `h_π(θ_π(t); t) = 0`, so the
   deployment-range objective has exact analytic gradients with no finite differences or adjoint solves.
2. **Classification: TRIVIAL-CONSEQUENCE-OF-KNOWN.** Implicit function theorem applied to a scalar equation
   with an explicit closed form for the equation's coefficients. This is first-year calculus once T3's
   quadratic-in-`t` coefficients are established; no citation is needed or claimed, and none was found.
3. **Dismissal sentence.** "Once you have `h(θ;t)=0` in closed form, `∂θ/∂t = −(∂h/∂t)/(∂h/∂θ)` is the
   implicit function theorem, not a contribution."
4. **Rebuttal.** Correctly filed by the bundle itself as routine algebra (T6.2 carries no `[A]`-vs-conjecture
   tension); its only real content is practical — avoiding finite-difference noise in an optimizer — and the
   genuine open problem sits one level up, in T6.H.1: `θ_π(t)` is **not differentiable** at double roots
   (`disc=0`) or where the active pair switches, which is precisely where a range-maximizing optimizer wants to
   sit, and the file states this is an unproved **CONJECTURE (T6-smooth)** about measure-zero non-smoothness
   along the optimizer path.
5. **To bullet-proof.** Nothing to bullet-proof in T6.2 itself. T6.3's active set is explicitly *not*
   certified — "no proved `O(n)` bound", downstream of H-LOC's refutation — and must never be called a
   "certified active set" in the paper (the file already enforces this ban on itself). T6-smooth needs either a
   subgradient/bundle argument or an empirical frequency count of non-smooth events actually hit during
   optimization runs.

---

## T7 — Rank / periodic corrections (`L = R·D`, out-harmonic characterization)

1. **Statement.** The paper's hole-constraint matrix `L` factors as `R·D` (partition times hinge-digraph
   incidence); its left null space is exactly the out-harmonic circulations on the hinge digraph; for
   boundary-free (periodic) patterns `rank(L) ≤ H−1`, refuting the paper's "#independent equations = #holes"
   (2026 §4.4) by exactly one.
2. **Classification: TRIVIAL-CONSEQUENCE-OF-KNOWN, with a genuine off-by-one correction to the published
   paper.** `L = R·D` and the harmonic/circulation characterization of `leftnull(D)` is standard discrete
   potential theory (graph Laplacian factorization; out-harmonic functions on digraphs are classical, related
   to random-walk/Markov-chain theory) — no author claims otherwise. Nothing in the screened literature states
   the periodic-boundary defect specifically for this problem.
3. **Dismissal sentence.** "`L = R·D` with `D` an incidence matrix and `R` a 0/1 partition is graph-Laplacian
   factorization from any algebraic graph theory textbook; the periodic rank deficiency is just Euler's formula
   applied to a torus quotient — not new mathematics."
4. **Rebuttal.** The mathematics is textbook, but the *finding* is a specific, checked correction to a
   published claim (2026 §4.4), with the correction shown to be **harmless in practice** only because the
   paper's own boundary rows happen to pin the missing translational degree of freedom — a fact the published
   paper does not state and would need to, since it is exactly the kind of silent cancellation a careful referee
   would want spelled out. It is also the piece of the theory bundle with the cleanest, fully closed
   verification loop: F22 measured `1ᵀL=0`, `rank(L)=H−1`, `dim Z=1` on boundary-free tori to `~1e−16`.
5. **To bullet-proof.** Nothing mathematically remains; this is solid, checked, and low-risk. The only
   presentational fix needed is to state it as a **correction to the 2026 paper's §4.4**, explicitly, rather
   than folding it silently into "background" — otherwise a referee who checks Eq. (4)'s row count against a
   periodic example will find the discrepancy independently and read it as a bug the authors missed.

---

## K7 — periodic Jacobian claims

**C1 — affine Jacobian `J(θ)=cos(θ/2)I+sin(θ/2)K`.**
Classification: **TRIVIAL-CONSEQUENCE-OF-KNOWN** — a direct corollary of T1's `Y_θ = cC+sS` restricted to the
period vectors of a torus quotient; nothing about the periodic case requires new theory beyond correctly
building the quotient (which K7 does do correctly, unlike the deviation-11 finite-patch shortcut it replaces).
Dismissal: "affine in `(cos θ/2, sin θ/2)` is inherited verbatim from T1 — this is a restriction, not a new
result." Rebuttal: the *quantitative* finding is new — `dim 𝒦 = 2·rank(D)`, not the spec's guessed
`min(4, 2·dim_null)`, falsified on `squares_3x3` where `dim_null=4` but `dim 𝒦=0` (`K` frozen at a pure
rotation). That the shape-space dimension does not bound designability is a genuine, checked (33/33, 1.04e−13)
correction to the project's own prior hypothesis, not the literature's. To bullet-proof: prove `dim 𝒦 = 2 rank(D)`
in closed form (it is currently `[N]`-only, verified but not derived, per `results/kill/KILL_REPORT.md:1216-1226`).

**C2 — exact conformality for all `θ`.**
Classification: **NEW**, closest prior art screened and ruled out. `notes/screen_bundle.md` (S4) reads
**Czajkowski, Coulais, van Hecke, Rocklin, *Conformal Elasticity of Mechanism-Based Metamaterials*, Nat.
Commun. 13 (2022), arXiv:2103.12683** in full, incl. appendices: their conformality is *spatial* (nonuniform
fields over an already-dilational continuum), not the *temporal* statement that first-order conformality at
`θ=0` propagates to every `θ`. Also checked and ruled out: **Konaković-Luković et al. 2018**, abstract-level
only (full text unobtained — flagged as a residual risk). Dismissal: "once `J(θ)` is affine in `(cos θ/2, sin
θ/2)` and similarities form a 2D linear subspace of `2×2` matrices, `K` a similarity forces `J(θ)` a similarity
for every `θ` — one line, not worth a claim number." Rebuttal: it is one line *given* C1, but it directly
settles 2026 §5.1's own hedge — "empirically ... for all θ" — turning an admitted empirical gap in the
published paper into a proof. That is a legitimate, citable contribution regardless of proof length.
To bullet-proof: get the Konaković-Luković 2018 full text (currently abstract-only) before printing "no source
supplies this" as a novelty claim.

**C4 — closed-form `θ_c`.**
Classification: **TRIVIAL-CONSEQUENCE-OF-KNOWN**, direct algebra from C1 (`A(θ)=det P₀(det J(θ)−1)`, half-angle
substitution). No prior art found or needed. Dismissal: "closed-form second root of a harmonic given the
harmonic's own closed form is not a separate theorem." Rebuttal: it is the periodic form the paper's spec asked
for and the 2026 paper does not state; matches forward kinematics to `5.5e−14` on 32/33 periodic + 4.4e−16 on
7 bounded patterns — cheap, correct, useful, but modest.

**Overall K7 caveat that undercuts C1/C2/C4's usefulness:** all three are statements about the achievable
Jacobian *set*; the constructive companion claim, C3 (hit a target `K*` with a certified, collision-free
design), **FAILS**, 48% vs an 80% bar, for the same `θ=0⁺` split-cut-folds-inward mechanism that kills K5/K6.
A referee will correctly note that C1/C2/C4 characterize a set most of which is *not reachable by any valid
design* — the paper's headline conformal/anisotropic-Poisson claims need the K5/K6/K7-C3 caveat attached every
time C1/C2/C4 are invoked, or they overclaim.

---

## F35 — hole-area harmonic and border functional

1. **Statement.** Per-hole area is a first harmonic, `A_C(θ) = a_C sin θ − b_C(1−cos θ)`, and the total border
   functional `Σ_C A_C` equals a closed-form quantity `2B(X)` including notches; periodically,
   `W + R = det P₀ · tr K` per representative preimage.
2. **Classification: TRIVIAL-CONSEQUENCE-OF-KNOWN.** Direct corollary of T3 (every predicate composed of
   `M′`-vertex positions is a first harmonic) applied to the shoelace-formula area of a hole polygon — no new
   argument beyond T3 plus signed-area bilinearity. The corollary "max-area angle is half of `θ_c(C)`" (B10) is
   again one line from the harmonic's own closed form (T3.4-style root/argmax algebra).
3. **Dismissal sentence.** "Area is a quadratic/bilinear functional of vertex positions and vertex positions
   are harmonic in θ by T3 — area is a harmonic by composition; nothing here required a new idea."
4. **Rebuttal.** The bundle concedes this itself — F35's own STATE.md entry frames it as verification of an
   identity, not a discovery, and explicitly notes the resulting design principle (a budget threshold `τ*` on
   `tr K`) **fails as an algorithm** (B3: the threshold co-varies with the very budget it's supposed to gate,
   AUC 0.95 vs an 0.8 bar but constrained re-solve is noise-level, 12→14/25). So F35 is honest, checked
   (8.7e−14 on 471 holes; ≤2.7e−16 border functional on 10 patches; 1e−15 on 30/33 periodic), and low-risk, but
   it is bookkeeping in service of an algorithm (budget-threshold design) that the same document kills.
5. **To bullet-proof.** Nothing mathematically; the only fix is presentational — do not let the paper cite F35
   as support for a design *algorithm* since B3 already shows the natural algorithmic use of it fails.

---

## Ranked

**Carries the paper (top 5):**
1. **T4 — exact `Θ_max` with graze handling.** Falsifies both the project's own round-1 claim and an implicit
   assumption in the 2026 paper (`Θ_max = min β`), on the paper's own figures. This is the theorem a referee
   will actually test by hand (the 4.8.8 arithmetic check is trivial to redo), and it survives.
2. **T3 — harmonic predicates with the systematic-degeneracy catalogue.** Everything downstream (T4, T5, T6,
   F35, K7) is dead without it; it is the one piece that is simultaneously exact, checked to `1e-12`, and load-
   bearing for five other claims.
3. **K7-C3's FAIL, read as a positive result about the *whole* bundle.** Not a claim in the classification
   table above, but the single most defensible scientific statement in the project: it independently reproduces
   the `θ=0⁺` collision obstruction (F30/K5/K6) on a population sharing *none* of their construction (periodic
   tori vs. random fixed-boundary patches). A referee respects a negative result reproduced across
   independent generators far more than another PASS.
4. **T7 — the periodic rank off-by-one.** A concrete, checked correction to a specific published equation
   (2026 §4.4), not a reproof of known linear algebra dressed up.
5. **T5/T6.4 — the sound certificate, with its own withdrawal history intact.** The round 1→3 correction
   trail (D4, D5, and the explicit withdrawal of "no overlap at θ=0" as *wrong*, not imprecise) is evidence the
   project's own error-catching process works, which matters more for credibility than the certificate's
   modest mathematical content.

**Should be demoted to remarks (bottom 3):**
1. **T6 — exact gradients.** Textbook implicit-function-theorem algebra; state it as an implementation detail
   inside T3/T5's discussion, not as a numbered theorem.
2. **T2 — no-locking.** Once T1's closed form exists, exhibiting its tangent vector in the flex kernel is
   routine; keep only the T2.4 corollary (no kinematic termination ⇒ contact list is complete) as a one-line
   remark feeding T4, and drop the theorem-with-corollaries framing.
3. **K7-C4 — closed-form `θ_c` (periodic).** One-line algebra from C1; fold it into the C1 discussion as "and
   the second root is...", not a separate lettered claim next to C1–C3.
