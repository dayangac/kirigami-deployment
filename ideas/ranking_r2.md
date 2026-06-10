# Critic round 2 — dedupe, adjudication, ranking, and the recommended bundle

Written 2026-09-04 by the Critic subagent (round 2). Output of `specs/critic.md` applied to the
round-2 idea files. Format follows `ideas/ranking.md` (round 1).

**Inputs actually read.** `specs/common_preamble.md`, `specs/critic.md`, `STATE.md` in full
(F1–F33, U1–U9, D1–D10, dead ends, session log), `notes/screen_r2.md` in full (580 lines),
`ideas/round2_theorist.md` (§0, Idea 1/2/3/8 in full, §2 self-attack, §3 ranked list, §4),
`ideas/round2_theorist_b.md` (B2/B3/B4 in full, §2 disagreements, §3 self-attack, §4 ranked list, §5),
`ideas/round2_inverse.md` (I9, I10 in full, self-attack, ranked list, caveats),
`ideas/round2_adversary.md` (Part 2 enabling identity, Ideas 1/2/8 in full, Parts 3–6),
`results/kill/KILL_REPORT.md` §A3 in full and its section index, `MISSION.md` §1 and §8,
`ideas/ranking.md` §0 (for format and for the round-1 facts P1–P3).

**Honesty policy.** I ran no code and no web searches. Every number is from `STATE.md`,
`notes/screen_r2.md` or `results/kill/KILL_REPORT.md`, and I say which. Two experiments are
**running and undecided** at the time of writing — B4 (`experimenter-b4`, no `results/kill/b4/`
directory exists yet) and the certificate repair (`checker-cert`, F32) — and the ranking below is
given **conditionally on both**, because both change what the other ideas mean rather than merely
how likely they are to work.

---

## 0. Three facts that reorder every round-2 file

**(Q1) The certificate is broken, and every "certified" count in every idea file is currently
uninterpretable.** F32 (A3, measured): `NOROOT(ε=0.006)` rejects **1622 of 1859** split-bearing
designs whose *exact* `Θ_max` is between 1.0 and 2.4 rad — an **87.3 % false-negative rate**. Every
kill rule in the four files that reads "certified on ≥ N of M" is, as written, a test of the
certificate and not of the idea. This includes B3's `12/50`, B4's own PASS bar, A2's soundness
control on the 8 tilings, and K7 C3's `12/25`. **No round-2 idea may be committed before
`checker-cert` lands.** K5/K6/K7-C3's FAIL verdicts survive F32 because they rest on *exact*
`Θ_max = 0` and on `min q < 0`, not on the certificate — but their *rates* do not.

**(Q2) Half the authored library is not a test of anything.** A3 measured that four of the eight
reference tilings (`rotating_squares`, `triangles_alternating`, `kagome_3636`,
`periodic_squares_4x4`) are **split-free**, so `dim_null = 0`, and Eq. (6) returns the authored
tiling exactly at every jitter amplitude up to `a = 1.0` (`‖X0 − X_base‖_∞ = 0.0000`). Every idea
whose validation set is "the 8 authored tilings" therefore has an effective *n* of 4. This kills
A2's soundness control as stated (it must run on the four split-bearing tilings), and it is why
A1's "`σ` outside `P(X)` on the 8 tilings kills the idea" is a weaker check than it sounds.

**(Q3) The project's only PASSED constructive machinery is periodic, and it is already at
theorem strength.** F33/K7: `J(θ) = cos(θ/2) I + sin(θ/2) K` exact to 1.4e-13 against independent
FK, `K` affine in the shape coordinate to 4.9e-14, `dim 𝒦 = 2·rank(D)` on 33/33 (C1 PASS);
conformal designs conformal at **every** `θ`, residual ≤ 1e-16 over 200 angles on 25/25 — which
settles 2026 §5.1's "empirically" (C2 PASS); `θ_c = 2·atan2(tr K, 1 − det K)` to 5.5e-14 on 32
periodic + 7 bounded designs (C4 PASS). The **only** thing that failed was C3, and it failed for a
reason that is not about the theory: reachability is perfect (target `K` hit to 8.9e-16 on the
designable patterns), and the obstruction is the same `0⁺` split-cut collision as K5/K6. So the
constructive half of the project needs exactly one thing: **a way to steer the design away from the
`0⁺` collision.** That is the criterion I rank against.

---

## 1. Adjudication: T-1 versus B4 on what causes the `0⁺` collapse

The team lead's framing is that K6 refutes T-1, because "all 400 satisfy Eq. (2) and still collapse".
**I do not accept that as a refutation, and the reason is F11.**

Eq. (2) is a constraint **per hole preimage**, not per vertex. F11 (measured 100/100 by Builder-Core,
and the basis of the verified `L = R·D` factorisation in F14/F22) says the preimages are
`C_K = E(K) ∪ {hinge edges whose target ∈ V(K)}`, one per connected component `K` of the split
subgraph. Only in the **pure** case (`E_split = ∅`) does this reduce to `C_v = in-edges of v`. So:

- On a design with split cuts, `L X = 0` asserts that the hinge in-edge vectors sum to zero **over
  the union of all vertices in each split component**. It asserts nothing about any individual
  vertex in a component of size ≥ 2.
- T-1's hypothesis (BAL) is the *per-vertex* statement `Σ_{u ∈ In(v)}(x_u − x_v) = 0`. On the K6
  population — random Voronoi/Delaunay patches, which by F15 have `dim_null = |E_split|` large
  (up to 1017 by F19) — split components are everywhere, so **(BAL) is exactly what fails at the
  split-touched vertices** while Eq. (2) holds globally.

**Therefore T-1 and K6 are consistent, and the theorist-a/theorist-b disagreement is not the
disagreement the lead states.** The real content of theorist-b's objection (§2a) is different and it
is correct: T-1's second "hence" — that zero-sum velocities positively spanning the plane implies no
copy enters a *particular* neighbouring corner cone — does not follow, because the corner predicate
is `μ = max(−g₁, −g₂)` (convex corner, a disjunction), not a half-plane condition.

**The measurement that decides it, and it is five minutes with no new code** (theorist-b's own
proposal, which I endorse verbatim): over the K6 population, restrict `corner_incidences` to
interior vertices touched by **no** split edge, and report how many have `μ ≤ 0`.
- **Count = 0** → T-1's safe half is true as stated, the collapse is *localised on split-touched
  vertices*, and T-1 becomes a real lemma (and the strongest NOT-FOUND item in the screen).
- **Count > 0** → T-1 is dead as stated, and the corner disjunction is the operative obstruction,
  which promotes A1 (whose whole technical content is that disjunction).

**T-1 and B4 are not competing hypotheses.** T-1 answers *which vertices collapse*; B4 answers
*whether the shape space has enough freedom to avoid it*. Both can be true, and under B4 PASS they
compose: the freed `|V_∂|` boundary dimensions (I3's count) are precisely the ones that can pay the
split components' balance defect. I record this as an adjudication, not a tie-break.

---

## 2. Master dedupe table

Six clusters. Personas: **A** = adversary, **T** = theorist-a, **B** = theorist-b, **I** = inverse.
"Screen" is `notes/screen_r2.md`'s verdict. C1–C7 are `specs/critic.md`'s constraints.

| id | cluster / merged idea | members | screen | C1 | C2 | C3 | C4 | C5 | C6 | C7 | verdict |
|---|---|---|--:|--:|--:|--:|--:|--:|--:|--:|---|
| **X1** | Expansive cone LP: non-uniform deployability is an LP, uniform is one ray | A1, A8(=A4 continuation), I10, T-6 (tangent cone), B9 (budget for every flex) | PARTIAL (Rote–Santos–Streinu 2003; must cite) | 1 | 1 | 1 | 1 | 1 | 1 | 1 | **PURSUE — run first** |
| **X2** | `tr K` budget threshold + budget-constrained periodic inverse design | B3, B2 (the `O(\|E\|)` scalar), B5 | not screened (B3 threshold); budget *principle* KNOWN via A5/I2 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | **PURSUE-IF-KILL-PASSES** |
| **X3** | Emptiness certificate by convex duality (S-procedure / Farkas) | A2, T-2, part of B5 | KNOWN as method (C-IRIS, Dai et al. IJRR); NOT FOUND as application | 1 | 1 | 0.5 | 1 | 1 | 0.5 | 1 | **PURSUE-IF-B4-FAILS** |
| **X4** | Boundary is the obstruction + free-boundary dimension count | B4, I3, B6 (notch leakage) | I3 NOT FOUND (two lines of Euler); B4 internal | 1 | 1 | 1 | 1 | 1 | 0.5 | 1 | **RUNNING — decides framing** |
| **X5** | Vertex balance defect / parity of the cut set | T-1, T-3, T-5, T-9 | T-1 **NOT FOUND** (strongest unfound); T-3 **KNOWN** = 2025 Remark 4.1 | 1 | 1 | 1 | 1 | 1 | 0.5 | 1 | T-1 PURSUE-IF-KILL-PASSES; T-3 **DROP** |
| **X6** | Curvature / area budget: hole area is a boundary functional, expansion floor, `det J = 1 + A/a` | A5, T-4, I2, B1, B10, I4, I7, I8 | **KNOWN** (Konaković-Luković, Panetta, Crane, Pauly, SIGGRAPH 2018) except T-4's *floor* and B3's threshold | 1 | 1 | 0 | 1 | 1 | 0 | 1 | **DROP as headline; keep T-4's floor as one lemma** |
| **X7** | 2025-facing audits and certified oracle | I9, I6, A9, A10 | errata-grade | 1 | 1 | 1 | 1 | 1 | 0 | 1 | **DROP as contribution; keep as errata** |
| **X8** | Already run or already known | A3/A-idea-4 (**done**, F32), A6(=A5), A7, T-7, B7, B8, I1, I5, O1–O5 | — | — | — | — | — | — | — | — | **CLOSED** |

**Convergence, recorded as evidence and not as a tie-break.** Four independent personas produced the
curvature/area budget (X6) and three produced an emptiness certificate (X3). The screen says X6 is
the *published* one and X3 is the *known-method* one. Convergence here measures obviousness, not
merit: the ideas four personas all reached are the ideas a referee will also reach.

---

## 3. Top 5 — hostile rejection, fair rebuttal, verdict

### #1 — X1: the expansive cone LP (A1), with A4 as its gated hero

**What it adds to the verified bundle.** T1–T7 characterise the *uniform* branch and prove
(T2.4) that the uniform branch never locks, so every stop is a contact. K7 C1/C2/C4 characterise
the *uniform* periodic Jacobian. A3 measured that the uniform branch dies sharply under jitter.
Nothing in the verified bundle says anything about a **non-uniform** flex, and 2026 §7 limitation
(ii) — verbatim in F8, "could be achieved by sensitivity analysis" — is the one stated limitation
of the source paper the project has not touched. X1 turns it into a decidable linear program, and it
does so by the one move nobody in this project has made: **freeing the flex from the embedding**, so
that `q_e` stops being quadratic in `t` and becomes linear in the flex.

**Hostile TOG reviewer.**
> Theorem 1 states that a set defined by strict linear inequalities on a linear subspace is an open
> polyhedral cone; Theorem 2 states that the uniform deployment is an element of that subspace.
> Neither is a theorem, and the cone itself is Rote–Santos–Streinu 2003, whose Lemma 3.2 the authors
> do not cite in their statement. What remains is an experiment reporting that an LP is feasible on
> some graphs — and the authors concede that at convex corners the feasible set is a *union* of
> polyhedra they cannot enumerate, so infeasibility proves nothing and the method has no negative
> half. Most damning: the flexes are infinitesimal, and the authors' own kill threshold for the
> finite version is 0.05 radians — under three degrees. A first-order rescue of a design space whose
> failure is measured at `θ = 0` is not a rescue; it is the tangent line to the failure.

**Fair rebuttal.** The screen (`screen_r2.md`, A1) read Rote–Santos–Streinu in full and found the
*paradigm* but not the *system*: RSS index velocities by point pairs, `O(n²)` non-strict
inequalities, no rigid bodies, no cuts, and pointedness is a conclusion rather than a branching
hypothesis. The cut-separation system on face-flex variables inside `ker A` is not there, and the
disjunction over convex corners is genuinely absent from that literature. The asymmetry
(feasibility proves, infeasibility does not) is the standard asymmetry of existence certificates and
is stated up front. The infinitesimal charge is exactly why A4 exists and is *gated* on A1 rather
than assumed — and the no-locking lemma T2.4 (AGREE, checked) makes A4's acceptance test complete:
every termination of a continuation is a contact, never a kinematic dead centre.

**Verdict: PURSUE — run `kill_k8a` first.** It is under ten minutes, its outcome is genuinely
unpredictable, and *both* outcomes are publishable: if non-uniform flexes rescue the graphs the
uniform ray fails on, F30 is reframed as "one ray of a 305-dimensional cone points the wrong way";
if they do not, the emptiness result becomes structural — a statement about the cut structure rather
than about two orientation heuristics — which is the single strongest defence against reviewer point
R-1 that the project can buy.

### #2 — X2: `tr K` is the budget, and budget-constrained periodic inverse design (B3)

**What it adds.** This is the only round-2 idea that is a **constructive deliverable** built on
machinery that has already PASSED. K7 C3 failed at 48 % against an 80 % bar, and F33 records exactly
why: reachability is never the problem (`dim 𝒦 = 4`, target hit to 8.9e-16), the `0⁺` split-cut
collision is. B3 supplies the missing steering term — a *linear* functional `tr K` of the design
variable that (B.3)+(B.5) claim controls `Σ_split r_e` — plus a threshold `τ* = W/det P₀`, plus an
algorithm that is one extra linear constraint inside a driver that already runs. It retrodicts C3's
own split (`diag(1,−0.5)`, `tr K = 0.5`, certified 0/25 with `min q < 0` on 25/25, against
`random_in_K` at 12/25). If it holds it converts K7 C3 from a dead end into the paper's algorithm,
with a ≥ 20-graph figure that MISSION §8(c) accepts.

**Hostile TOG reviewer.**
> The authors observe that a linear functional of their design variable controls a first-order area
> rate and christen it a budget. Stripped of vocabulary: `det J(θ) − 1` is the areal expansion, its
> derivative at zero is `½ tr K`, and asking a pattern to expand less makes it open its holes less.
> Equations (11)–(13) of Segall et al. 2026 already parametrise `J(θ)`; that `tr K` is its area
> derivative is one line of calculus any reader supplies. The contribution reduces to a threshold
> `τ*` which the authors' own risk paragraph admits is a linearisation of unknown quality, because
> `W` is not constant along the achievable set. What is measured is that one hand-picked anisotropic
> target fails 25 times and random targets succeed 12 times — a correlation with `n = 50`, no
> controlled variation of `tr K`, and no evidence that the deployment budget rather than the
> anisotropy is doing the work.

**Fair rebuttal.** The reviewer is right about the identity and right that the retrodiction is a
correlation. The paper stands or falls on **kill step (ii)**, the constrained re-solve on the
*identical* targets — a controlled experiment, not a correlation — and I would not let it ship
without it. The reviewer is wrong that the closed form is one line of theirs: 2026 §5.1 constrains
`J` only at `θ = 0` and says nothing about whether the design deploys, and K7 C2 had to *prove* the
all-`θ` statement the paper asserts empirically.

**Verdict: PURSUE-IF-KILL-PASSES**, and the kill must be re-run after `checker-cert`, because the
`12/50` baseline it is measured against is a certificate count and F32 says certificate counts are
currently 87 % false-negative.

### #3 — X3: the emptiness certificate by convex duality (A2 ≡ T-2)

**What it adds.** F25/F30's headline numbers are `0/200` and `0/400` *samples*. X3 replaces them
with per-graph proofs plus a named witness subgraph. It is the only idea that answers 2026 §6/§7
limitation (i) — "geometric characteristics for favourable deployment behaviour directly from the
embedding remains an open problem" — in the negative direction, and its constraint matrices have an
explicit rank-≤4 factorisation from blocks the code already stores (`ZeroPlusForm::GS`, `::DD`).

**Hostile TOG reviewer.**
> This is the S-procedure. It is sixty years old, it is known to be lossy for more than one
> quadratic, and its use for exactly this purpose — certifying emptiness of a collision set by
> Positivstellensatz solved as an SDP — is a published, named robotics method (C-IRIS; Dai, Amice,
> Werner, Zhang, Tedrake, IJRR) that also has the completeness converse these authors concede they
> lack. Their degree-2 relaxation is the lowest rung of a hierarchy the prior work climbs, and the
> literature's answer to their failures is "raise the degree". A decision procedure that succeeds on
> 25 of 50 graphs, with no characterisation of which 25, is a heuristic with a verification step.
> And the enterprise is self-referential: a proof system built to shore up a negative result the
> same submission produced, about a population — random Voronoi patches — that nobody cuts.

**Fair rebuttal.** Soundness is exact: `Q(λ) ⪯ 0` implies emptiness with no sampling, and `0⁺`
emptiness implies `Θ_max = 0` by one line of the `sin(θ/2)` expansion, so this is a *decision
procedure sound always and complete sometimes*, reported with its coverage — which is how every
practical Positivstellensatz result is reported. The application to a kirigami design space is
genuinely unfound (screen, A2). The self-referential charge is the fatal one and it is
**conditional**: it holds if the emptiness result is the paper, and dissolves if X1 or X2 carries the
positive half.

**Verdict: PURSUE-IF-B4-FAILS.** Under B4 PASS this idea certifies the emptiness of a set the 2026
paper *chose* rather than one the problem forces, and that is a much weaker sentence — theorist-b
§2d is right about this and it is the sharpest thing in that file.

### #4 — X5(T-1): the vertex balance defect

**What it adds.** The only NOT-FOUND *mechanism* in the whole round-2 bundle (screen: "strongest
unfound idea"). Everything else in the project explains the collapse by measurement; T-1 explains it
by a per-vertex algebraic defect and predicts *where* the corner contact fires.

**Hostile TOG reviewer.**
> The proved half is vacuous — it says that in the pure-hinge case, where the authors' own F18
> already establishes `Θ_max = min(min β, π) > 0` unconditionally, the copies of a vertex do not
> collide. That is a restatement of a fact they list as given, and it reads as an explanation of why
> Segall et al. 2025's even-valency class works, i.e. as support for someone else's Remark 4.1. The
> half that would matter is stated with an unnamed constant `C(v)`, and the experimental section
> then reports a rank correlation — which is what one reports instead of a theorem. And the
> hypothesis `k_v ≥ 3` excludes `k_v = 2`, the single most important vertex type in the kirigami
> literature. A theorem about kirigami that does not cover rotating squares is a theorem about
> nothing.

**Fair rebuttal.** The pure half is not vacuous *as a proof*: F18's `Θ_max = min β` is a measurement
on 8 patterns and T4.4's `≥` direction is explicitly untagged in `check.md`, so a proof of `0⁺`
safety at pure vertices closes a real gap in the project's own theory. My §1 adjudication removes the
lead's objection (Eq. (2) is per-preimage, not per-vertex, by F11), so the idea is not refuted by K6.
What remains true is that the *quantitative* half is unproved and the fallback is weak.

**Verdict: PURSUE-IF-KILL-PASSES**, and the kill is five minutes with no new code (§4, `kill_t1`).
This is the cheapest idea in the entire round-2 bundle. Run it alongside `kill_k8a`.

### #5 — X6 reduced to one lemma: the expansion **floor** (T-4), plus B2's `O(|E|)` bound

**What it adds.** Almost nothing survives the screen here, and I rank the cluster fifth only after
cutting it down to two items. The screen is unambiguous: the budget design principle, the conformal
scale-factor check, and the Gauss–Bonnet cone-singularity consequence are all
Konaković-Luković, Panetta, Crane and Pauly, SIGGRAPH 2018, with four verbatim quotes; the published
maximum areal expansion of the triangular Kagome linkage is **exactly 4**, which is I2's hemisphere
threshold, so **I2's stated kill rule fires on a published number before any code runs**. What is not
published is (a) T-4's identity as a **floor** — a lower bound on the expansion rate, where every
published bound is a ceiling, obstructing *low*-expansion targets — and (b) B2's one-pass scalar
certificate `min_e r_e ≤ (2B − W)/|E_split|`.

**Hostile TOG reviewer.**
> That a surface cannot be wrapped by a material expanding by at most `c` unless its area distortion
> is below `c` is not a theorem about kirigami; it is the definition of area distortion. `det J =
> 1 + A/a` is conservation of area, and the "budget identity" is Green's theorem with a face
> potential glued on. The validation set is a single documented failure in a prior paper's
> supplement that the authors state they cannot access.

**Fair rebuttal.** Accepted in full for A5, I2, B1, B10, I4 and I7 — those are one section of
related work, not a contribution. The floor is different in kind and I found no counterpart to it,
and B2 is a microsecond rejection test with a hard mathematical kill rule (any design with
`2B − W < 0` and `min_e q_e > 0` refutes it outright).

**Verdict: DROP as a headline; PURSUE T-4's floor and B2 as two lemmas inside whichever bundle
ships.** Merge A5 + T-4 + I2 + B1 into one section led by the 2018 paper, per the screen's action 1.

---

## 4. Kill specs — exact, no further design decisions

**Gate 0 (blocks every rule below that mentions "certified").** `checker-cert` must land first, and
its own acceptance test is fixed by A3's data: on the four **split-bearing** reference tilings at
jitter amplitude `a = 0.005`, where `POS 20/20`, `NOOVERLAP 20/20` and exact `Θ_max ∈ [1.0, 2.4]` on
20/20 seeds, the repaired certificate must pass on **≥ 18 of 20 seeds per tiling** (currently 0, 0,
5, 10 of 20). Until then, every rule below is stated against **exact `Θ_max`** and `min_e q_e`, which
F32 confirms are unaffected.

### `kill_t1` — T-1's vertex balance defect (5 min, no new code, run first)
- **Population.** The K6 designs, seeds recoverable from `results/kill/k6/shards_v4/*.csv` (400 rows).
- **Routine.** `zero_plus.hpp::corner_incidences` at the K6 `X0`; `make_cut` for the split-edge set.
- **Quantity.** Restrict incidences to interior vertices **not** incident to any split edge. Report
  `N_pure` = number of such vertices, and `N_bad` = number of them with `μ ≤ 0`.
- **Tolerance.** `μ ≤ 0` at `1e-12` relative to the local edge scale.
- **PASS.** `N_bad = 0` over all 400 designs, with `N_pure > 0` on at least 300 of them (otherwise
  the test is vacuous — random Voronoi patches may have almost no pure vertices, and if `N_pure`
  is near zero the correct verdict is INCONCLUSIVE, not PASS).
- **FAIL.** Any `N_bad > 0` → T-1 dead as stated; the corner disjunction is the mechanism, which
  transfers the weight to X1.
- **Wall.** Minutes. **This is the cheapest decisive test in the round-2 bundle.**

### `kill_k8a` — X1's expansive cone LP (< 10 min, run second)
- **Population.** The 100 K1a graphs (101–793 faces), at `X = X_ini`, which F25 certifies injective
  200/200 so the LP data are valid. Plus the **four split-bearing** reference tilings as the
  soundness control (per Q2 — *not* all eight).
- **Routine.** `mobility.hpp::build_A`; dense `ColPivHouseholderQR` kernel basis on the **2-core
  only** (F29: sparse rank is unreliable above ~700 columns — do not use SparseQR here);
  `zero_plus.hpp::split_copies` and `::corner_incidences` for the inequality rows; branch at each
  convex corner = the side the copy is already on at `X_ini`, one repair pass.
- **Quantity.** `max_z min_i l_i(N z)` s.t. `‖z‖_∞ ≤ 1`, by multiplicative weights on the dual
  `min_{λ ∈ simplex} ‖Nᵀ l(λ)‖₁`. Report the certified margin and, on infeasible branches, the Farkas
  `λ`.
- **Tolerance.** Margin `> 1e-8` after normalising each inequality row to unit norm.
- **PASS.** Strictly positive certified margin on **≥ 20 of 100** graphs. **FAIL** below that: the
  emptiness is structural, not a wrong ray.
- **Hard soundness kill.** If the LP reports `σ ∉ P(X)` on any of the four split-bearing tilings
  (which have exact `Θ_max` between 1.6 and 2.4 rad), the driver is wrong and must be fixed before
  the number is quoted.
- **Wall.** `m ≈ 300`, ≈ 2000 rows, ~1 s per graph.

### `kill_b3` — X2's `tr K` threshold and constrained re-solve (after Gate 0)
- **Population.** K7's 50 C3 rows verbatim, `AchievableSet` and the pattern list reused unchanged.
- **Step (i), the theorem.** Record `tr K` at each solved design and `τ* = W/det P₀`; also record
  the **variation of `W` over `𝒦`** by evaluating `W(t)` at 20 points of the achievable set — this
  is the risk B3 names and it must be measured before the threshold is quoted.
- **Step (ii), the algorithm.** Re-solve with `tr K ≥ τ* + δ|τ*|` for `δ ∈ {0, 0.1, 0.5}`, re-certify.
- **PASS (theorem).** AUC of `tr K − τ*` for `min_e q_e > 0` ≥ 0.8 over the 50 rows, **and** no row
  with `min_e q_e > 0` at `tr K < τ*`. **PASS (algorithm).** The constrained re-solve raises the
  certified rate on the *same* targets above the re-measured post-Gate-0 baseline by a margin that
  survives the 50-row binomial error bar.
- **FAIL (soft).** Theorem holds, algorithm does not → B3 demotes to a corollary of B1 and X2 dies
  as a deliverable.
- **FAIL (hard).** `W` varies over `𝒦` by an amount comparable to `det P₀ · tr K` → the threshold is
  not a threshold, it is a penalty, and it is not a theorem.
- **Wall.** K7's C3 shards ran in minutes; this adds one linear constraint.

### `kill_b2` — the `O(|E|)` budget bound (5 min, shares B1's app)
- **Population.** The 400 K6 designs, the 8 reference tilings, the 33 K7 periodic designs.
- **Quantity.** `B(X)`, `W(X)`, `min_e q_e` per design.
- **HARD KILL.** Any design with `2B − W < 0` and `min_e q_e > 0` — the inequality is violated.
- **SOFT KILL.** `2B − W > 0` on all 400 K6 designs — true but never binding, demote to a footnote.
- **Note.** B1's global sum (B.4) is **unverified** (theorist-b §5: the geometric tracer returns
  holes but not notches), so `B(X)` itself is unvalidated on bounded patches. This test validates
  the identity and the bound in one pass and must report the identity residual separately.

### `kill_k8b` — X3's emptiness certificate (< 20 min, only if B4 FAILS)
- **Population.** The 50 smallest K1a graphs (so `m = 2·dim_null ≲ 400`).
- **Routine.** `ZeroPlusForm`; `M_e` per split edge via the rank-≤4 factorisation
  `M_e = ½(G_eᵀ E D_e + D_eᵀ Eᵀ G_e)` from `GS`/`DD`; `M_i` per corner incidence at the branch
  active at `t = 0`; `M_f` per face area. 2000 multiplicative-weights iterations; normalise every
  `M_i` to unit Frobenius norm and report the normalisation as part of the certificate.
- **PASS.** `λ_max(Q(λ)) ≤ −1e-8` on **≥ 25 of 50** graphs.
- **Hard soundness kill.** Any certificate found on the **four split-bearing** tilings (Q2 — the
  split-free four have `dim_null = 0` and are not a control) is a bug in the driver.
- **Wall.** ~20 s per graph, shardable.

### `kill_k8h` — A4's certified non-uniform continuation (30 min, **gated on `kill_k8a` PASS**)
- **Population.** 20 K1a graphs with exact `Θ_max = 0`.
- **Routine.** Integrate the LP flex by predictor–corrector on the pin-constraint variety; accept a
  step only when `contact.hpp::swept_discs` with the exact radius `max(‖x‖,‖χ‖)` (KILL_REPORT
  correction 3) plus an exact endpoint test report no overlap.
- **PASS.** Median `θ_eff` ≥ 0.05 rad over the 20 graphs. **FAIL** below → the first-order rescue
  does not integrate, and X1 keeps only its LP half.

### Shared kill experiments
`kill_t1` and `kill_k8b` read the same K6/K1a populations and can share a driver.
`kill_b2` and B1's identity check are one app. `kill_b3` is inside `kill_k7.cpp`.

---

## 5. The ranking, conditional on B4 and on the certificate repair

B4 is not really a candidate idea; it is a **framing switch**, and it flips the meaning of X3 and of
the project's headline negative result. I give the ranking under both outcomes.

### Under **B4 FAIL** (still `0/400` with the boundary rows dropped)

Emptiness is structural: it is about the cut structure of random planar graphs, not about 2026's
Eq. (4). The negative result hardens and deserves a proof.

| # | idea | why |
|--:|---|---|
| 1 | **X1** expansive cone LP (+ A4 gated) | the only positive half left, and 2026 §7(ii) is verbatim unaddressed |
| 2 | **X3** emptiness certificate | the negative result becomes theorem-with-witnesses, and the self-referential objection is answered because X1 carries the positive half |
| 3 | **X2** `tr K` budget + constrained design | the constructive deliverable, now the *periodic* escape hatch from a structurally empty patch space |
| 4 | **T-1** vertex balance defect | explains *which* vertices, feeds X3's witness reading |
| 5 | **T-4 floor + B2** | two lemmas; B2 is the `O(\|E\|)` front end to X3's SDP |

### Under **B4 PASS** (free boundary gives a strictly positive feasible rate)

F30 stops being "the random-graph shape space is empty" and becomes **"2026 Eq. (4)'s fixed boundary
empties it"** — a correction to the source paper rather than an extension of it, and the single
most quotable sentence available to this project. X3 drops hard, because certifying emptiness of a
set the paper chose is a much weaker claim, and I3's dimension count becomes load-bearing rather
than decorative because it counts exactly the freed dimensions.

| # | idea | why |
|--:|---|---|
| 1 | **X4** = B4 + I3 | the headline: identical graphs, three boundary conditions, certified rate; I3 supplies the closed-form count `k = \|V_∂\| + \|E_split\| + 1 − c(Γ)` of the dimensions that were being suppressed |
| 2 | **X2** `tr K` budget + constrained design | unchanged in value, and now consistent with X4 (the periodic case is the free-boundary limit) |
| 3 | **X1** expansive cone LP | still the answer to 2026 §7(ii), but less urgent: the space is not empty, so the uniform ray is no longer the only suspect |
| 4 | **T-1** vertex balance defect | explains the residual failures under a free boundary |
| 5 | **T-4 floor + B2** | unchanged |

### Interim B4 data as of writing (read by me from disk; the verdict is `experimenter-b4`'s to give)

`results/kill/b4/shards/shard_*.csv` had 586 completed rows when I read them. Aggregated by variant:

| variant | rows | `0⁺` feasible | exact `Θ_max > 0` |
|---|--:|--:|--:|
| `fixed` (K6 control) | 198 | 1 | 0 |
| `free` (boundary rows dropped, 1 vertex pinned) | 195 | 0 | 0 |
| `split_free` | 193 | 0 | 0 |

This is **partial and unofficial** — the run is unfinished and I did not verify the driver — but it
points hard at **B4 FAIL**, and the free-boundary column is *worse* than the fixed control, not
better, which is the direction theorist-b's own risk paragraph predicted (F17's inverted faces get
worse without `B`). Plan on the B4-FAIL ranking and treat the B4-PASS column as the contingency.

**B4's own PASS bar must be the full certificate, not `0⁺` feasibility.** Theorist-b's self-attack is
right that a free-boundary design that is `0⁺`-feasible and self-intersecting in the flat state is
not a design, and F17 says inverted faces get *worse* without `B`. And its second self-criticism is
right and cheap to fix: the `11/12` torus evidence has no fixed-boundary control on the *same*
graphs, so B4 must cut a disk-topology patch from the same point set and run both.

---

## 6. Novelty holes — what the Scout must still screen

Every claim below currently rests on the round-2 files' internal reasoning, on `notes/field_*.md`
rows, or on my own reading. Exact query strings:

1. **B3's trace threshold.** Never screened. `periodic auxetic linkage "areal expansion" linear
   functional of Jacobian threshold deployable "unit cell" inverse design closed form`; and
   `"Poisson ratio" target periodic mechanism design "achievable set" affine Jacobian kirigami
   constraint collision`.
2. **B2's one-pass certificate.** `kirigami OR "hinged tiling" non-deployability certificate
   "single pass" edge sum bound "opening angle" necessary condition closed form`.
3. **B7 (universality of `K`), B8 (`θ_c(P) + θ_c(P*) = 2π`).** Both marked `[SCREEN-ME]` by their
   own author. B8 in particular: `rotating squares rotating triangles dual tiling "closing angle"
   sum 2 pi auxetic reciprocal pattern duality`.
4. **A4 / X1's continuation.** `Borcea Streinu "expansive periodic mechanisms" arXiv:1507.03132`
   full text — the screen surfaced but did not fetch it, and it is the nearest periodic analogue of
   the cone.
5. **T-1's proof pattern.** `Connelly Demaine Rote "Straightening polygonal arcs"` DCG 30 (2003)
   full text — the screen could not access it, and T-1 is a local version of its non-collision
   argument. The exact wording must be checked before T-1 is written up.
6. **Konaković et al. 2016 "Beyond developable"** — the 2018 paper supersedes it for X6, but its
   numeric bounds are unverified and it is where the bounded-conformal-factor initialisation
   originates.
7. **The 2025 supplement (Fig. F.4).** Not available locally. Every claim in X6 that "the 16×4 table
   validates against the paper's reported successes and failures" rests on a validation set nobody
   on this project has seen. This is an *evidential* hole, not a novelty hole, and it is on its own
   sufficient reason not to make X6 the headline.

**Two housekeeping actions the screen already prescribed and that are not yet done.** Add
Konaković-Luković, Panetta, Crane, Pauly (SIGGRAPH 2018) as its own row in
`notes/field_kirigami.md`. Rewrite T-3's framing to the T-join bound before any effort is spent on
it, because the equivalence it claims is 2025 Remark 4.1 verbatim.

---

## 7. Recommended bundle and paper skeleton

**Title of the bundle.** *Closed-form deployment of corner-hinged kirigami on arbitrary planar
graphs: exact ranges, achievable periodic Jacobians, and where the design space is empty.*

**The theorem** (already proved and independently checked; MISSION §8(a)). The trig-linear
deployment `Y_θ = cos(θ/2) C(X) + sin(θ/2) S(X)` (T1) and the exact contact calculus T4.2″ with the
crossing and interval clauses, giving `Θ_max` in closed form with graze handling; plus the periodic
half from K7: `J(θ) = cos(θ/2) I + sin(θ/2) K` with `K` affine in the shape coordinate and
`dim 𝒦 = 2·rank(D)` (C1), conformality at **every** `θ` whenever it holds at `θ = 0` — which settles
2026 §5.1's "empirically" (C2) — and the second closed angle `θ_c = 2·atan2(tr K, 1 − det K)` (C4).
The load-bearing lemma is T2.4: the uniform branch never locks, so **every** stop is a contact,
which is what makes the contact calculus a complete characterisation rather than a bound.

**The characterisation, validated at MISSION §8(b) scale.** Exact `Θ_max` against brute-force
bisection. K2a has it on **187/187** to 2e-10, and the naive min-over-roots rule errs by up to
1.047 rad on 5 of them. Scale the same driver to **≥ 500 random graphs** — this is a re-run, not new
work, and it is the single cheapest way to satisfy §8(b) outright.

**The algorithm.** X2 (`tr K` budget-constrained periodic design) if `kill_b3` passes; X1's LP flex
if `kill_k8a` passes and `kill_b3` does not; both if both pass, in which case X1 is the patch-space
algorithm and X2 the periodic one.

**The ≥ 20-graph figure.** ≥ 20 random **Voronoi tori** with prescribed Poisson-ratio / conformality
targets, certified deployable under the budget constraint, plotted against the same graphs designed
by the unconstrained 2026 §5.1 route, with each design's `tr K − τ*` on the x-axis. This is the
MISSION §8(c) figure and it exists only because K7 C1 established that the target is always reachable
and C3 established that reachability is never the failure mode.

**The hero.** An arbitrary random Voronoi periodic graph deployed to a **prescribed anisotropic**
Jacobian, certified collision-free to a closed-form `Θ_max`, exported through `kiri_export` to laser
SVG in the closed and open states and to 3MF. Neither Segall paper produces this object: 2026's
demonstrated periodic examples are authored tilings, and this project has measured its random-graph
patch space to be empty under four strategies.

**The negative half, and how to write it.** The emptiness result (F25/F30, `0/200` and `0/400` under
four strategies) plus A3's sharp jitter transition (`a*` from 0.16 to 0.52 median-edge units on the
four split-bearing tilings, transition width 0.36–0.61 decades). A3 is what defends the negative
against "your generator did the work", and it must lead that section — with the honest disclosure
that four of the eight authored tilings are split-free and therefore trivially invariant. Add X3's
certificates if B4 fails; drop the section to a remark about Eq. (4) if B4 passes.

**Errata section** (F14 `rank(L) ≤ H − 1` for periodic patterns against 2026 §4.4's prose claim;
F23 the authors' extra boundary-component row is an over-constraint; F24 `merge_close_verts` fuses
hinge duplicates so their `θ_max` reads 0.067 where the true first contact is ≈ 1.65 on
`snub_square`; F28 Eq. (9)'s output re-closes on 9.47 % of designs; the 2026 Eq. (18) `2π`/`2φ`
typo). Errata do not carry a paper and must not be presented as contributions.

**What I would cut.** The whole curvature-budget cluster as a claim (A5, I2, B1, B10, I4, I7, I8) —
one related-work paragraph led by Konaković-Luković et al. 2018, plus T-4's floor as a lemma. T-3
(2025 Remark 4.1). I1 (`m ≡ 0` on all three regular families it is motivated by). I9 and the audit
family (errata). Everything in `round2_adversary.md` Part 3.

---

## 8. Final recommendation

**Run, in this order, and nothing else until they report:**
1. `checker-cert` (already running) — Gate 0; no "certified" number means anything until it lands.
2. `kill_t1` — five minutes, no new code, settles the §1 adjudication and the fate of T-1.
3. `kill_k8a` — ten minutes, settles whether the project has a positive half on patches.
4. B4 (already running) — settles the framing of the negative half.
5. `kill_b3` — settles whether K7 C3 becomes the paper's algorithm or stays a dead end.

**Commit to X2 as the constructive core** unless `kill_b3` fails hard, because it is the only
round-2 idea that builds on machinery that has already PASSED (K7 C1/C2/C4), it repairs the only K7
criterion that failed, its driver already exists, and it produces the ≥ 20-graph figure and the hero
that MISSION §1 requires on disk. **Fallback order:** X1 (+A4) → X4 if B4 passes → X3 if B4 fails.

**The single sentence a TOG reviewer would reject the top candidate with.**
For X2: *"The authors observe that the derivative of `det J` at zero is `½ tr K`, add that number as
a linear constraint to a design problem already parametrised by Equations (11)–(13) of the paper they
extend, and report that constraining it helps — on fifty rows, with the threshold's own key quantity
`W` admitted to be non-constant along the set it is defined over."*
For X1, if it is promoted: *"A first-order rescue of a design space whose failure is measured at
`θ = 0`, on a cone whose paradigm is Rote–Santos–Streinu 2003, with an acceptance threshold of 0.05
radians that concedes nobody expects these flexes to integrate."*
