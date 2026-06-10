# IDEA — The usable shape space of Tutte auxetic kirigami

**Exact deployment range, a sound certificate, certified non-deployability of the published
embedding, and a range-maximising constrained embedding that deploys where it does not.**

Every number below cites the file it was measured in. Nothing is quoted from memory, and
nothing from `STATE.md`'s "Unverified claims" appears unless `STATE.md` records it as later
verified.

---

## 1. The claim, in one paragraph

Segall, Ren and Sorkine-Hornung (TOG 2026) give a linear shape space `X` of uniformly
deployable kirigami on arbitrary planar graphs, and state as an open problem that membership
in `X` says nothing about whether the design actually opens
(`notes/paper_2026.md:949-951`, quoted in full in §6 below). We close that problem in three
parts. First, the **deployment range is exactly computable in closed form**: every geometric
predicate along the uniform branch is a single first harmonic in the opening angle with
coefficients quadratic in the flat design coordinates, so the collision-free range `Θ_max` is
a minimum over closed-form roots, with grazing contacts handled explicitly
(`derivations/core.md` T3, T4). Validated against a tolerance-independent bisection referee
on **3,113 designs, agreement 3,113/3,113 at 1e−5, worst gap 3.822e−09**
(`results/final/e1/E1.md`). Second, a sound sufficient **certificate** `Θ_max ≥ ε`, exact
below the first graze (`derivations/check.md` Round 6), holds on 2,150 of those 3,113 designs
with **0 soundness violations** and is **exact on that population**
(`results/final/e1/E1.md`). Third, and this is the finding that changes the framework's
practical reach: the paper's own Eq. (6) projection is **certifiably non-deployable on random
planar graphs** — 0 of 400 designs under four repair strategies, with per-graph dual
certificates covering even non-uniform flexes (`results/kill/KILL_REPORT.md` §K1a, §K5, §K6,
§B4, §K8a) — and the obstruction is not the graph but the projection, which manufactures
reflex corners that jam the mechanism at `θ = 0⁺` (`results/kill/KILL_REPORT.md` §T-1).
Replacing the projection's proximity objective by the **0⁺ deployment margin**, inside
convexity and split-outward barriers in the same null space, makes **307 of the same 400
designs (76.8 %)** deployable and certified, against **0** for every baseline
(`results/final/figures/summary_table.md`, `results/kill/k9c/`).

---

## 2. Theorems — what is actually new

The hostile theory review (`review/theory_review.md`) classified every item. Its demotions
are binding here: the following are the claims that survive as theorems, and the rest are
listed as lemmas because they are trivial consequences of known results.

* **T4 — exact `Θ_max` with graze handling.** A complete closed-form candidate list of contact
  angles, correctly separating a grazing vertex-vertex touch from an overlap. This falsifies
  the implicit `Θ_max = min_e β_e` reading of the published bound *on the paper's own hexagon
  figure*: `θ₁ = π/3` is a graze, `Θ_max = 2π/3` (`derivations/core.md:665-693`). The proved
  direction is `Θ_max ≤ min(min β, π)`; the reverse is **measured on 8 patterns and not
  proved** (`derivations/core.md:785-800`) and must never be printed as an equality.
  Ranked first by the theory review.
* **T3 — harmonic predicates with the systematic-degeneracy catalogue.** The harmonic *form*
  is Weierstrass-level algebra and is not claimed; what is claimed is that the coefficients
  are quadratic in the *design* coordinates, plus the catalogue of systematic degeneracies
  (`p+q=0` at every permanent incidence and split-edge pair, 73,446 of 2,145,387 measured
  occurrences, `review/theory_review.md` T3). The review requires this catalogue be stated as
  a lemma with a proof rather than a measurement; that is an open item.
* **T5 — the certificate: sound, not complete.** `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT(ε)` implies
  `Θ_max ≥ ε` (`derivations/core.md` T5.2b′, verified `derivations/check.md` Round 6).
  Proved scope: **exact for `ε` below the first graze, a strict inner approximation above it**
  — the hexagon at `ε = π/2` has range but fails the certificate (`STATE.md` F32).
  Completeness is **false**, and the earlier phrasing "no face–face interior overlap at
  `θ = 0`" is **withdrawn as wrong**, not merely imprecise (`derivations/core.md:1802-1810`).
* **T7 — the periodic rank correction.** `L = R·D`, left null space = out-harmonic
  circulations, and for boundary-free patterns `rank(L) ≤ H − 1`, so the published
  "#independent equations = #holes" (2026 §4.4) is **off by one** on periodic patterns.
  Measured `1ᵀL = 0`, `rank(L) = H − 1`, `dim Z = 1` to ~1e−16 on boundary-free tori
  (`STATE.md` F22). The mathematics is textbook; the *correction to a published equation* is
  the contribution, and it must be stated as such. The correction is harmless in practice only
  because the paper's boundary rows pin the missing translation — a silent cancellation the
  paper does not state.
* **K7 C2 — exact conformality for all `θ`.** With `J(θ) = cos(θ/2)I + sin(θ/2)K` and `K`
  affine in the shape coordinate, first-order conformality at `θ = 0` propagates to every `θ`:
  residual `≤ 1e−16` over 200 angles on 25/25 designable patterns
  (`results/kill/KILL_REPORT.md` §K7). This settles the 2026 paper's own hedge in §5.1, where
  the property is asserted "empirically ... for all θ". Screened NEW; the nearest prior art,
  Czajkowski, Coulais, van Hecke and Rocklin (*Nat. Commun.* 13, 2022, arXiv:2103.12683), is
  a *spatial* conformality of a continuum, read in full and ruled out (`notes/screen_bundle.md`
  S4).

**Lemmas, not theorems** (demoted per `review/theory_review.md`, and stated as consequences):

* **T1** trig-linear deployment `Y_θ = cos(θ/2)C(X) + sin(θ/2)S(X)` — the Tay–Whiteley
  motion assignment made trig-linear and specialised to corner-hinged tilings; the half-angle
  magnitude is anticipated, under two extra hypotheses, by Acuña et al., *Auxetic behavior on
  demand*, Commun. Phys. 5 (2022), arXiv:2101.12352.
* **T2** no-locking (`σ ∈ ker A(Y_θ)` for all `θ`) — routine once T1's closed form exists.
  Only the corollary survives as a remark: no kinematic termination, hence T4's candidate list
  is complete.
* **T6** exact gradients by implicit differentiation of the root map — the implicit function
  theorem. Never call its active set "certified": `H-LOC` is **refuted**, not unproved, and
  there is no proved `O(n)` bound (`derivations/check.md` D2).
* **K7 C1/C4** — `J(θ)` affine and `θ_c = 2·atan2(tr K, 1 − det K)` follow from T1; the new
  quantitative part is `dim 𝒦 = 2·rank(D)` on 33/33 (1.04e−13), which falsifies the project's
  own prior guess `min(4, 2·dim_null)` on `squares_3x3`. `dim 𝒦 = 2 rank D` is verified, not
  derived.
* **F35** — hole area is a first harmonic `A_C(θ) = a_C sin θ − b_C(1 − cos θ)` (471 holes,
  8.7e−14) and sums to a border functional (≤2.7e−16 on 10 patches). This is T3 composed with
  the shoelace formula. It must **not** be cited as support for a design algorithm: the natural
  algorithmic use of it, a `tr K` budget threshold, is killed by B3 in the same corpus.
* **T-3** — the T-join framing survives; the equivalence itself is 2025 Remark 4.1 verbatim
  (`notes/screen_r2.md`).

---

## 3. The characterization

`results/final/e1/E1.md`, N = 3,113 designs (authored 1,660; jittered 1,322; random 131),
spanning 13 families and both orientation rules:

| quantity | value |
|---|---|
| exact `Θ_max` (T4.2″) vs bisection referee, agreement at 1e−5 | 3,113 / 3,113 |
| worst gap | 3.822e−09 |
| certificate holds (`ε = 0.006`) | 2,150 |
| soundness violations | 0 |
| designs with `Θ_max ≥ ε` | 2,150 — certificate **exact** on this population |

MISSION §8b asks for a characterization validated on ≥ 500 graphs against brute force; this
exceeds it six-fold. The naive rule `Θ_max = min` over roots errs by up to **1.047 rad** on
5 of 187 designs (`results/kill/KILL_REPORT.md` §K2a) — the graze clause is load-bearing, not
bookkeeping.

---

## 4. The certified impossibility, and its mechanism

On 200 random Voronoi/Delaunay/quad graphs (`F ∈ [101, 793]`, median 332.5) under both
orientation rules:

| strategy | certified deployable |
|---|--:|
| Eq. (6) projection alone, `σ_mc` (K1a) | 0 / 200 |
| both `σ` rules, sampled and projected (K5) | 0 / 400 |
| four 0⁺-repair strategies (K6) | 0 / 400 |
| boundary rows dropped (B4) | 2 / 400 exact, **1 / 400** refereed |
| non-uniform expansive-cone LP (K8a) | 0 / 100 |
| authors' native `prevent`, capped side-run (K6) | 0 / 8 |
| authors' full native pipeline (Native200, stopped at 573/600 runs on 193 graphs) | 0 / 172 completed by our exact scan; 5 / 573 by their own collision test (false negatives, K1c mechanism); 362 timed out at 600 s, 39 crashed |

K8a is the strong form: 84/100 designs at `X_ini` and 68/100 at `X0` carry near-Farkas dual
certificates below 1e−6, i.e. **no non-uniform flex in `ker A` avoids the 0⁺ contacts**
either. The certificates are verified sound — `λ ≥ 0` exactly, `|Σλ − 1| ≤ 2.5e−13`, 96/100
and 99/100 below 1e−6 at 1e5 dual steps (`STATE.md` F38). **Caveat, mandatory:** these
certificates cover the **4 visited branch charts of 2^2338**, never all of `P(X)`
(`STATE.md` F38).

**The mechanism.** The obstruction is manufactured by the projection, not by the graph. The
input embeddings have **zero** non-convex face corners; `X0` has **22,103 of 156,220**
(`results/kill/KILL_REPORT.md` §T-1). Of 1,968 balanced pure vertices with a non-positive 0⁺
margin, **1,966** sit at a reflex corner. At a convex corner the margin is a disjunction
`max(−g₁, −g₂)`; at a reflex corner it is a conjunction `min(−g₁, −g₂)`, and that is the one
that fails. The reflex-corner mechanism is screened **NOT FOUND** in the literature
(`notes/screen_k9.md`, item (b)).

A second negative result, from K9b, is reported because it is informative rather than
because it is favourable: under the proximity objective, **feasibility of the constraint pair
does not determine deployability.** Nine designs come back from a four-times-larger search
exactly feasible at the full margin and still `Θ_max = 0`
(`results/kill/KILL_REPORT.md` §K9b). That anomaly is what motivated changing the objective,
and it does not survive the change: in the completed 400-design run under the range objective
no margin-feasible-but-jammed design remains
(`results/final/figures/fig_feasible_vs_deployable.png`).

---

## 5. The algorithm

**Two-stage range-maximising constrained embedding inside the Eq. (4) null space**
(`code/src/method/range_embed.{hpp,cpp}`). Same population, same shape space
`X(t) = X0 + Φt`, same barriers (`cross_i ≥ δ`, `q_e ≥ δ′`), same certificate as K9; only the
objective changes. **Stage A** maximises the 0⁺ margin
`m(X) = min(min_e q_e, min_j μ_j)/med²` through a log-sum-exp softmin that is a **proved
lower bound on `m` for every `t`**, with analytic gradients checked against central
differences; the reported `m` is always the exact minimum, never the surrogate. **Stage B**
maximises the exact `Θ_max` in a shrinking trust region with exact rejection.

Yield on the identical population, all refereed at both shrinks:

| method | deployable | `ε_max ≥ 0.1` rad | `ε_max` median |
|---|--:|--:|--:|
| Eq. (6) alone / `σ_def` / K6 repairs / K8a / convexity alone | 0 | 0 | — |
| B4 (boundary dropped) | 1 / 400 refereed | — | — |
| authors' native `prevent` | 0 / 8 (capped side-run) · 0 / 172 completed in Native200 (5 / 573 by their own test) | — | — |
| K9 (proximity objective) | 36 / 400 = 9.0 % | 33 | 0.25 rad |
| K9b (4× search, proximity) | 27 / 400 | 14 | 0.115 rad, max 1.828 |
| **K9c (range-maximising)** | **307 / 400 = 76.8 %** | **207** | **0.254 rad**, max π |
| best of the three arms, per design | 307 / 400 | 214 | 0.278 rad, max π |

Every positive is refereed (307/307) and `ε_max ≤ referee` holds **400/400**. Per
orientation, `σ_mc` 155/200 and `σ_def` 152/200; per family, Voronoi 122/134, Delaunay
108/134, quad-random 77/132; taking the better `σ` per graph, **185 of the 200 graphs**
deploy. The two proximity arms measured on this same population give 36 and 27, which is the
control that says the jump is the objective and not a change of population or referee.
**Yield is not scale-free:** it falls with face count, Delaunay from about 1.0 to 0.8 and
quad-random from 0.75 to 0.43 across `|F|` bins with Wilson 95 % intervals
(`results/final/figures/fig_yield.png`). In the completed run **no margin-feasible-but-jammed
designs remain** (`results/final/figures/fig_feasible_vs_deployable.png`): the nine such
designs that motivated K9c were a property of the proximity arms, not of the constrained set.

Convexity and split-inward feasibility is reached on 256 of the 400 designs, 244 of them
with a strictly positive 0⁺ margin; stage B supplies the winning point on 157 of the 307
positives, and the `t = 0` start rather than a proximity warm start supplies 278 of them.
Of the 93 failures, `split-inward` binds on 57, `inverted` on 27 and `vertex-edge` on 9.

**Row count used: 400 of 400**, from `results/kill/k9c/{k9c.csv, summary.txt}` (re-merged
from all 12 shards) and `results/final/figures/summary_table.md`, which agree. The 182-row
interim snapshot was concordant at 77.5 %.

**Hero:** K9 population graph `delaunay` id 148, `σ_mc`, F = 101, exact
`Θ_max = 1.9967778150149833` rad, `ε_max` equal to it, bisection referee
1.9967778152171005, certified, 0 inverted faces (`export/hero/README.md`, locked as a
regression in `code/tests/test_design.cpp`). Exported closed, at `Θ_max/2` and at `0.9 Θ_max`
as SVG, STL and 3MF, plus a 6-frame deployment sequence; all 9 SVGs pass `xmllint`, all 3
3MFs pass `unzip -t`, and every solid reports closed, consistently oriented, 0 boundary and 0
non-manifold edges.

---

## 6. Relation to prior work

* **Segall, Ren, Sorkine-Hornung 2026, Limitation 1**, quoted in full because clipping it
  would misrepresent it (`notes/paper_2026.md:949-951`): *"First, although all embeddings in
  the solution space `X` are theoretically uniformly deployable, some may be geometrically
  undesirable (e.g., exhibiting extremely short edges or sharp angles). Moreover, uniform
  deployability does not guarantee a large collision-free deployment range; self-intersections
  may occur at small opening angles (e.g., Fig. 14). At present, such issues can only be
  detected through explicit deployment simulation via forward kinematics and require
  additional post-processing to resolve. Developing geometric characteristics for favorable
  deployment behavior directly from the embedding remains an open problem."* This is exactly
  the problem addressed. It is agnostic about frequency: it does **not** concede that the space
  is generically empty, and the paper never claims a success rate over any population.
* **Choi, Dudte, Mahadevan 2019** (*Nat. Mater.*; arXiv:1812.08644) already use the same
  convexity primitive — the cross-product sign at a corner, Eq. (4) `⟨(b−a)×(c−a), n⟩ ≥ 0` —
  solved with `fmincon`. Screened PARTIAL: no null space, no Tutte embedding, no barrier/SQP
  machinery, and no split-inward condition on duplicated vertices (`notes/screen_k9.md`).
* **Rote, Santos, Streinu 2003** give an exact polyhedral, LP-decidable cone of first-order
  expansive velocities for pointed pseudo-triangulations — the closest prior "certificate" for
  a cognate problem. Screened PARTIAL: their cone is unconstrained by cut structure
  (`notes/screen_r2.md`). Must be cited, not competed with.
* **Grima et al. 2012** (*Proc. R. Soc. A* 468, 810–830) give a closed-form locking angle for
  the hinge-adjacent two-triangle case and concede "this will not be possible in the general
  case" — the closest prior *use* of the harmonic contact algebra (`notes/screen_bundle.md` S3).
* **Acuña et al. 2022** (*Commun. Phys.* 5, arXiv:2101.12352) give a three-step sufficient
  recipe for `ν = −1` counter-rotation on bipartite polygon networks, requiring collinearity
  and a fixed 2-colouring; T1 removes both (`STATE.md` F21).
* **Tay–Whiteley** body-and-hinge motion assignments: T1/T2 are the 2D specialization. The
  project's citation is currently second-hand through Handbook excerpts and must be replaced
  by the primary source (`notes/screen_bundle.md:306`).
* **Konaković-Luković et al. 2018** (SIGGRAPH): expansion bound = conformal scale-factor
  bound, with cone singularities where it fails; this subsumes the curvature-budget family,
  which was consequently dropped as a headline (`notes/screen_r2.md` A5).
* **arXiv:2608.30032** (Jiang and Choi), "A unified geometric design framework for kirigami
  structures" — a concurrent 2026 preprint, **now cleared**. The full text was retrieved and
  screened claim by claim against C1–C10 (`notes/screen_jiang_choi.md`, Tables A and B): it
  has zero collision, self-intersection, deployment-range, null-space and certificate content,
  and works on `N × N` rotating-squares quad patterns rather than arbitrary planar graphs. The
  only overlap is its Eq. (2), the cross-product corner-convexity inequality, so **C9 remains
  PARTIAL on that primitive alone**, already attributed to Choi–Dudte–Mahadevan 2019. This
  supersedes the abstract-only screen in `notes/screen_k9.md`.

---

## 7. Mandatory caveats

These attach to the claims above every time they are made. They come from the three hostile
reviews and are binding.

1. **Population scale.** `F` median 332.5, range [101, 793]; every reference case drawn from
   the 2026 paper's own figures has `F ∈ [4, 97]`. The random population is 3–8× larger by
   face count than anything the paper demonstrates. It is a representative sampling procedure
   at a scale and split-density regime the paper never tests — both halves of that sentence
   must be said (`review/negative_review.md` §3).
2. **The baseline is not the authors' full pipeline.** The `0/400` figures are Eq. (6) plus
   our own repairs. Native `prevent` at published defaults was actually invoked on **8**
   designs. Until Native200 reports, every such sentence reads "all four of *our* repair
   strategies are 0/400 certified; the authors' native `prevent`, tested separately on the 8
   of these 400 designs that reached a valid flat state, is also 0/8"
   (`review/negative_review.md` §5). Native200 (F40; `results/kill/native200/`) ran their full pipeline unconditionally on 193 of the 200 graphs, three variants, and was stopped by the user at 573/600 runs: 172 completed, 362 timed out at 600 s, 39 crashed. The crash cause was reproduced (F42): an empty-vertex-list access in `convert::to_eig_mat`, reached from `Hmesh::merge_close_verts` after a deployed configuration collapses to zero faces — not the `UnitPattern::get_holes()` segfault the earlier F24 note assumed. `results/kill/native200/crashfix.patch` guards it and also adds the missing `prev()->twin()` null check in `get_holes()`, a second latent defect. By our exact scan none of the 172 completed outputs opens; their own forward-kinematics collision test reports 5 positives, which our scan shows are its false negatives (the K1c mechanism). The baseline is therefore zero at scale, with the caveat that our median |F| = 332 exceeds their largest example (|F| ≤ 97) and that 63 % of their runs did not finish.
3. **K7's achievable-Jacobian caveat.** C1/C2/C4 characterize the achievable Jacobian *set*;
   the constructive companion C3 **fails**, 48 % against an 80 % bar, for the same 0⁺
   split-cut mechanism. Most of that set is not reachable by any valid design, and the caveat
   attaches every time C1/C2/C4 are invoked (`review/theory_review.md`, K7 overall).
4. **Slivers.** The K9c gallery was inspected: top designs open to `Θ_max = π`, and slivers
   remain a geometric-quality caveat (`STATE.md` F37). Short edges and sharp angles are
   precisely what the paper's Limitation 1 also names.
5. **The certificate is an inner approximation above the first graze.** It is exact only
   below it (`STATE.md` F32). Completeness is false and is stated as false.
6. **`ε_max ≤ referee`, always.** Every reported certified range is below the refereed exact
   range: 3,113/3,113 in E1, 400/400 in K9c. Where the two referees disagree (the `π/4000`
   bisection grid stepping over a 6.7e−01-wide thin collision interval on `delaunay 34
   σ_def`), the exact scan is the smaller and the conservative number, and the certificate is
   below both (`results/kill/KILL_REPORT.md` §K9c).
7. **"Deployable" needs its qualifiers.** Under an *added* convexity + split-inward
   constrained embedding, not the paper's own Eq. (6)/(9) pipeline; on a specific 200-graph
   Delaunay/Voronoi/quad population, with a yield that falls with face count; with median
   certified opening 0.254 rad (≈15°), not full deployment; certified sound but produced by a non-convex heuristic with no completeness
   guarantee, with feasibility ≠ deployability directly observed
   (`review/constructive_review.md` §4).

---

## 8. What died, and the number that killed it

| idea | verdict | killing number |
|---|---|---|
| R1 range-margin claim vs the authors' code (K2b) | FAIL | median relative gain **0.0000** on 30 designs; ours better 7, theirs 15, tie 16 |
| R3 bifurcation spectrum (K3a) | FAIL | `m_core ≤ 5` on **0 of 167** Delaunay patches, median 305 |
| A8 "Eq. (2) sufficient not necessary" | dead — already published | 2025 Fig. 10; Dang 2021 proves iff for quads |
| Our own Eq. (9) claim (P3) | withdrawn | authors' native `prevent` reaches `Θ_max = π` on snub square |
| K7 C3 target-driven periodic design | FAIL | 12/25 = 48 % vs an 80 % bar |
| T-1 vertex balance as a predictor | FAIL | balance holds exactly (2e−14) yet 1,968 pure vertices still collapse |
| B4 boundary drop | FAIL | 2/400 exact, 1/400 refereed vs a 10 % bar |
| K8a expansive-cone LP | FAIL | 0/100 positive margin, dual-certified structural on 84/100 |
| B3 `tr K` budget as a design constraint | FAIL | spread(W)/spread(det P₀ tr K) median 0.73, max 1.40 — the threshold co-varies with the budget |
| K9 / K9b constrained embedding, proximity objective | FAIL by 3–4 designs | 36/400 = 9.0 % against a 10 % bar |

The last row is the one that mattered: its failure diagnosis — that the objective, not the
feasible set, was the cap — is what K9c then confirmed.
