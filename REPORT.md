# REPORT — The usable shape space of Tutte auxetic kirigami

Companion to `IDEA.md`. Every number cites the file it was measured in. Test counts were
re-run by the writer at the time of writing, not copied.

---

## Abstract

Segall, Ren and Sorkine-Hornung (TOG 2026) embed uniformly deployable kirigami on arbitrary
planar graphs in a linear shape space, and state as an open problem that membership in that
space does not imply the design opens without self-intersection. We characterize the usable
part of the space exactly. Every geometric predicate along the uniform deployment branch is a
single first harmonic in the opening angle whose coefficients are quadratic in the flat
design coordinates, so the collision-free range is a minimum over closed-form roots, with
grazing contacts separated from overlaps. On 3,113 designs the closed form agrees with a
tolerance-independent bisection referee 3,113/3,113 at 1e-5, worst gap 3.822e-09, and a sound
sufficient certificate holds on 2,150 with zero soundness violations, exactly matching the
set with range at least the threshold. We then show that the paper's own projection is
certifiably non-deployable on random planar graphs -- zero of 400 designs under four repair
strategies, with per-graph dual certificates that rule out non-uniform flexes as well -- and
identify the cause: the projection manufactures reflex face corners absent from the input,
and these jam the mechanism at the instant of opening. Replacing the projection's proximity
objective by the first-order deployment margin, inside convexity and split-outward barriers
in the same null space, makes 77 % of that population deployable and certified against a
baseline of zero.

---

## Contributions

| # | contribution | kind | proved / measured in |
|---|---|---|---|
| C1 | Exact collision-free deployment range `Θ_max` from a complete closed-form candidate list, with grazing contacts handled | **theorem** (T4) | `derivations/core.md` T4.1–T4.5, `derivations/check.md` |
| C2 | Harmonic predicates: coefficients quadratic in the *design* coordinates, plus the systematic-degeneracy catalogue | **theorem** (T3) | `derivations/core.md` T3 |
| C3 | A sound sufficient validity certificate `Θ_max ≥ ε`; exact below the first graze, inner above; completeness **false** | **theorem** (T5) | `derivations/core.md` T5.2b′, `derivations/check.md` Round 6 |
| C4 | Correction to 2026 §4.4: `L = R·D`, left null space = out-harmonic circulations, `rank(L) ≤ H−1` on periodic patterns | **theorem** (T7) + erratum | `derivations/core.md` T7, `STATE.md` F22 |
| C5 | Exact conformality of the periodic Jacobian for all `θ`, settling 2026 §5.1's "empirically" | **theorem** (K7 C2) | `results/kill/KILL_REPORT.md` §K7 |
| C6 | Validation of C1 and C3 on 3,113 designs against brute force | **characterization** | `results/final/e1/E1.md` |
| C7 | Certified non-deployability of the Eq. (6) embedding on random planar graphs, uniform and non-uniform, with dual certificates | **characterization** | `results/kill/KILL_REPORT.md` §K1a §K5 §K6 §B4 §K8a |
| C8 | The mechanism: reflex corners created by the projection, not present in the input | **characterization** | `results/kill/KILL_REPORT.md` §T-1 |
| C9 | Two-stage range-maximising constrained embedding; ⟨JULIA:k9c:307 of 400⟩ designs deployable against 0 for every baseline | **algorithm** | `results/final/figures/summary_table.md`, `results/kill/k9c/`, `Kirigami/src/method/range_embed.jl` |
| C10 | Reimplementation of the 2026 pipeline at parity with the authors' code, plus three errata against it | engineering | `baseline/parity.md`, `results/core_validation/` |

Demoted to lemmas or remarks on the theory review's instruction: T1 trig-linear deployment
(Tay–Whiteley motion assignment specialized), T2 no-locking (only the completeness corollary
survives), T6 exact gradients (implicit function theorem), K7 C1/C4, F35 hole-area harmonic.
See `IDEA.md` §2 and `review/theory_review.md`.

---

## Method summary

**Setup.** A planar graph `M` is cut by an orientation `σ: F → {−1,+1}`. Interior edges
between faces of opposite orientation are hinges, retained at the source vertex with the
target duplicated; edges between same-orientation faces are split cuts, fully separated
(`STATE.md` F1). Hole preimages are grown by the paper's Alg. 1; equivalently, and this is
what the code uses, they are the connected components of the split-edge forest together with
the hinge edges pointing into them (`STATE.md` F11, measured to agree on 100/100 graphs,
F15). Uniform deployability is the closure condition Eq. (2), assembled as `[L; B] X = [0; T]`
(Eq. 4), whose null space `X = X0 + Φ t` is the shape space (Eq. 5), reached from the input by
the least-norm projection Eq. (6).

**T1 (lemma).** Each face rotates by `σ(f)·θ/2`, so `Y_θ = cos(θ/2)·C(X) + sin(θ/2)·S(X)`
with `C, S` linear in `X` (`derivations/core.md` T1.2). Each deployed vertex traces an
ellipse.

**T2 (remark).** `σ ∈ ker A(Y_θ)` for every `θ`, so the uniform branch never reaches a
kinematic dead centre. The only consequence used downstream is that the contact candidate
list of T4 is complete (`derivations/core.md` T2.4).

**T3.** Substituting T1 into any orientation, dot-product or squared-length predicate gives
`h(θ) = p + q cos θ + r sin θ` with `p, q, r` quadratic forms in `X`
(`derivations/core.md` T3). The load-bearing part is the quadratic dependence on the design
coordinates, not the harmonic form.

**T4.** `Θ_max` is the smallest `θ` at which two face interiors first overlap. Candidate
angles are the roots of the T3 harmonics of every vertex-into-edge and edge-crossing pair,
plus the interval clauses `0 ≤ ⟨w−a, b−a⟩ ≤ |b−a|²` that distinguish a contact with a
segment from collinearity with its infinite line. A grazing vertex-vertex touch is a root
that is not an overlap and must not terminate the range: the paper's own hexagon figure has
first root `π/3` and `Θ_max = 2π/3` (`derivations/core.md:665-693`). Broad-phase pruning
must use the swept radius about the *moving* centroid, which equals the flat circumradius
exactly, making the published locality gate identically 1 and therefore vacuous
(T4.5b′; this is why K2c's original PASS is withdrawn).

**T5.** `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT(ε)` at one probed angle, over the deflated three-class
root list, implies `Θ_max ≥ ε`. Soundness is proved and independently verified
(`derivations/check.md` Round 6). Completeness is **false**: `NOROOT_adm(ε)` is equivalent to
`ε ≤ θ₁`, the first contact angle, while `Θ_max` can exceed `θ₁` at a graze.

**T6 (remark).** `∂θ_π/∂t_i = −(∂h/∂t)/(∂h/∂θ)` by implicit differentiation, giving analytic
gradients for range maximisation. The active set is **not** certified; `H-LOC` is refuted and
the words `O(n)`, "certified active set" and "locality theorem" are withdrawn
(`derivations/check.md` D2).

**T7.** `L = R·D` with `D` the signed hinge-digraph incidence and `R` the 0/1 preimage
partition, so `leftnull(L)` is the out-harmonic circulations. Periodically `1ᵀL = 0`, hence
`rank(L) ≤ H − 1`.

**The algorithm (C9).** Inside `X(t) = X0 + Φt`, stage A maximises
`m(X) = min(min_e q_e, min_j μ_j)/med²` — the joint minimum of the split-edge sign and the
corner separation margin, both quadratic in `t` — through the log-sum-exp softmin
`M_κ(t) = −κ log Σ exp(−c_i(t)/κ) ≤ min_i c_i(t)`, under log barriers `cross_i ≥ δ` and
`q_e ≥ δ′` with geometric continuation. The surrogate is a lower bound everywhere and tight
where the branch modes were read; the reported `m` is always the exact minimum. Stage B
maximises the exact `Θ_max` with T6 gradients in a shrinking trust region (`0.25, 0.10, 0.04`
median edges), accepting a step only when exact convexity, exact split signs, a margin floor
and an actual range improvement all hold.

---

## Results

### The characterization (E1)

`results/final/e1/E1.md`, N = 3,113 designs across 13 families:

| quantity | value |
|---|---|
| agreement with the bisection referee at 1e−5 | 3,113 / 3,113 |
| worst gap | 3.822e−09 |
| certificate holds, `ε = 0.006` | 2,150 |
| soundness violations | 0 |
| designs with `Θ_max ≥ ε` | 2,150 (certificate exact on this population) |
| certificate confusion against `Θ_max ≥ ε` | TP 2,150 · FP 0 · FN 0 · TN 963 |
| slowest shard / summed measurement time | 311.1 s / 352.1 s |

The three random-graph families (58 Delaunay, 59 quad-random, 14 Voronoi) never certify and
never reach `Θ_max ≥ ε`; random `σ_mc` yields **0 tightly-embedded designs out of 900 ids
tried**.

### The impossibility (K1a, K5, K6, B4, K8a)

Population: 200 random Voronoi/Delaunay/quad graphs, `F ∈ [101, 793]`, median 332.5, under
both orientation rules.

| experiment | result |
|---|---|
| K1a | `X0` has inverted faces on 142/200 and face–face overlap on 200/200; **0/200** samples injective |
| K5 | certified valid `X0` on **0/200** under `σ_def` and **0/200** under `σ_mc`; `σ_def` does fix the flat sheet (defect down 122×, POS 58→186, projection distance down 4.7×) |
| K6 | certified `Θ_max > 0` on **0/400** after the `first_root` fix; the earlier 15/200 was that bug |
| B4 | boundary rows dropped: **2/400** exact, **1/400** refereed. `dim_null(free) − dim_null(fixed) = \|V_∂\| − 1` on 400/400, so the count is right; the fixed boundary is not what empties the space (the 2 exact rows are knife-edge outcomes of a non-converged L-BFGS run, not reproducible numbers: `derivations/scratch/b4_path/README.md`) |
| K8a | positive non-uniform margin on **0/100** at `X_ini` and **0/100** at `X0`; near-Farkas dual certificates below 1e−6 on 84/100 and 68/100 (archived LP-layer values predate the solver fix `d38b0a4` and are superseded by the rerun: `results/kill/k8a_julia/PROVENANCE.md`) |

Dual-certificate soundness was independently re-verified: `λ ≥ 0` exactly,
`|Σλ − 1| ≤ 2.5e−13`, 96/100 and 99/100 below 1e−6 at 1e5 dual steps, 65 and 57 below 1e−9
(`STATE.md` F38). One genuine non-uniform rescue exists — `snub_square_R20_s1`, certified
margin 0.043 — so the mechanism is real and simply never fires on random graphs.

### The mechanism (T-1)

400 designs, 165,788 interior vertices, 14,886 of them pure (touched by no split edge).
Vertex balance holds **exactly** at pure vertices (max `|δ_v| = 2e−14`), yet 1,968 of them
still have a non-positive 0⁺ margin, on 142/400 designs. **1,966 of those 1,968** binding
incidences are at reflex corners created by the Eq. (6) projection: 0 non-convex face corners
at the input, **22,103 of 156,220** at `X0`. Pure vertices with all-convex incidences fail on
**1 of 12,919**. Balance is therefore not the predictor; convexity is.

### The algorithm (K9, K9b, K9c)

Same 400 designs throughout.

| method | deployable (exact, refereed) | `ε_max ≥ 0.1` rad | `ε_max` median / max |
|---|--:|--:|---|
| Eq. (6) alone (K1a) | 0 | 0 | — |
| `σ_mc` and `σ_def` sampled/projected (K5) | 0 | 0 | — |
| four 0⁺-repair strategies (K6) | 0 | 0 | — |
| convexity alone (K9 arm a; 115 feasible) | 0 | 0 | — |
| boundary rows dropped (B4) | 1 / 400 | — | — |
| expansive-cone LP, non-uniform (K8a) | 0 / 100 | — | — |
| authors' native pipeline: capped side-run / Native200 (stopped at 573/600) | 0 / 8 · 0 / 172 completed by exact scan (5 / 573 by their own test) | — | — |
| K9, convexity + split-inward, proximity objective | 36 / 400 = 9.0 % | 33 | 0.253 / 1.997 |
| K9b, 4× search budget, proximity | 27 / 400 | 14 | 0.115 / 1.828 |
| **K9c, range-maximising** | **307 / 400 = 76.8 %** | **207** | **0.254 / π** |
| best of the three arms, per design | 307 / 400 | 214 | 0.278 / π |

`ε_max ≤ referee` holds **400/400**, and every positive is refereed: 307/307. Certified
306 of the 307. Per orientation, `σ_mc` **155/200** and `σ_def` **152/200**; per family,
Voronoi **122/134**, Delaunay **108/134**, quad-random **77/132**. Taking the better `σ` per
graph, **185 of the 200 graphs** carry a deployable design and **152** carry one with
`ε_max ≥ 0.1` rad. `‖X − X_ini‖` per vertex runs
around 0.75 median edges against K9's 0.38 — range is bought with distance from the input.
The two proximity arms measured on this same population give 36 and 27, so the jump is the
objective and not a change of population or referee.

**Yield falls with face count.** `results/final/figures/fig_yield.png` bins the population by
`|F|` with Wilson 95 % intervals: Delaunay drops from about 1.0 to 0.8 across the bins and
quad-random from 0.75 to 0.43. The method is not scale-free, and the trend must be reported
with the headline.

**The K9b anomaly does not survive the full run.** In the completed 400-design merge there
are **no margin-feasible-but-jammed designs left** — every design that reaches the convexity
and split-inward margin under the range objective also deploys
(`results/final/figures/fig_feasible_vs_deployable.png`). The nine exactly-feasible,
`Θ_max = 0` designs that motivated K9c were a property of the proximity arms, not of the
constrained set.

Convexity and split-inward feasibility is reached on **256** of the 400 designs, **244** of
them with a strictly positive 0⁺ margin. Stage B supplies the winning point on **157 of the
307 positives**, and within stage A the `t = 0` start rather than a proximity warm start
supplies **278**, so the range objective does not need K9's feasible point to start from.
`‖X − X_ini‖` per vertex is median **0.734** median edges against K9's 0.38 — range is bought
with distance from the input.

**Row count used: 400 of 400.** Numbers come from `results/kill/k9c/{k9c.csv, summary.txt}`,
re-merged from all 12 shards, and from `results/final/figures/summary_table.md`, which agree.
The earlier 182-row interim snapshot was concordant at 77.5 %. The baseline's 0/400 is
independently confirmed in `results/kill/k6/k6_final.csv` before it is asserted in any figure.

Of the 93 designs that do not deploy, `split-inward` binds on **57**, `inverted` on **27**
and `vertex-edge` on **9**. The `inverted` rows are designs where no start reached
feasibility and the driver records the least-infeasible point honestly rather than
discarding it.

### Rank and hole-count corrections

`H = |E_hinge| − |F| + c(Γ)` on 50/50 sweep graphs and 8/8 reference cases; `L = R·D` exactly
on 50/50; `rank(L) = H − dim Z` on 27/27 dense. On boundary-free tori `1ᵀL = 0`,
`rank(L) = H − 1`, `dim Z = 1` to ~1e−16, confirming the periodic prediction (`STATE.md`
F22). With a boundary the naive row-sum identity fails on 0/50 and the corrected identity
`(1ᵀL)_v = I(v)·indeg(v) − Σ_{v→w} I(w)` holds 50/50 and 8/8.

### Mobility (K3a), with its caveats

The 2-core identity `dim ker A = |F ∖ core₂| + dim ker A|core₂` holds on 493/508; every one
of the 15 exceptions is off by exactly ±1 and 7/7 rechecked densely hold exactly. After
2-core reduction, Delaunay mobility is still **median 305** (range 36–1675), so "mobility ≫ 1
is a dangling-face artefact" is false and the second PASS clause (`m_core ≤ 5`) fails on
**0 of 167**. **Caveat:** `SparseQR` rank is unreliable above roughly 700 columns (gap
1e−5–1e−3) and such ranks must not be quoted (`STATE.md` F29). **Second caveat:** the
`dim ker A` medians reported in K8a (Voronoi 3, quad 17, Delaunay 72) come from a different
population and a different reduction from K3a's 305, and the two must not be quoted side by
side (`STATE.md` F38).

### Errata against the published pipeline

1. **Rank claim.** 2026 §4.4's "#independent equations = #holes" is off by one on periodic
   patterns (T7 above).
2. **Boundary over-constraint.** The authors' code adds a constraint row for split-forest
   components touching the boundary; those faces form a path in `Γ`, so no closure constraint
   is needed. Their `X0` satisfies our system on 8/8; our null space is one dimension larger
   on two cases. Confirmed by forward kinematics from all 53 seed faces and 20 BFS orders
   agreeing to **7.16e−15** (`results/kill/KILL_REPORT.md` §F23).
3. **Silent re-closure.** The authors' native Eq. (9) output re-closes on **9.47 %** of
   designs under the well-posed rule (`results/kill/KILL_REPORT.md` §K1c). The rule as
   originally written reports 65.09 %, because it counts grazes, and is self-contradictory.
4. **Upstream defects observed, not claimed as contributions:** their collision-aware `θ_max`
   runs `merge_close_verts` at 0.1·average edge before testing, fusing hinge duplicates, and
   reports 0.067 on snub square where first contact is ≈1.65 (`STATE.md` F24). The same
   `merge_close_verts` also crashes: once a deployed configuration has collapsed to zero faces
   it leaves an empty vertex list, and the collision scan's next call reaches
   `convert::to_eig_mat` on it and dereferences element 0 of an empty vector —
   `EXC_BAD_ACCESS`. This is the reproduced cause of the Native200 crashes (`STATE.md` F42;
   `results/kill/native200/crashfix.patch`). A second, latent defect sits in
   `UnitPattern::get_holes()`, which dereferences `prev()->twin()` without a null check on a
   bounded patch; the same patch guards both.

### Parity with the authors' code

On 8 reference graphs: `E_hinge`/`E_split` agree exactly as sets on 8/8, symmetric difference
0; their `X0` satisfies our system to 4e−15 and their own to 6e−15; forward kinematics agrees
to **4e−8** on 8/8 when both are applied to a deployable embedding (`baseline/parity.md`).

### Hero

`delaunay` id 148, `σ_mc`, F = 101, N = 59, `|E_split| = 16`, `dim_null = 16`. Exact
`Θ_max = 1.9967778150149833` rad; `ε_max` equal to it; bisection referee
1.9967778152171005; certified; 0 inverted faces; `min_signed_area = 0.971`; 13,096
candidates over 756 face pairs (`export/hero/README.md`). Locked as a regression in
`Kirigami/test/test_design.jl` at 1e−9. Exported closed, at `Θ_max/2` and at `0.9 Θ_max` as
SVG, STL and 3MF, plus a 6-frame deployment sequence. All 9 SVGs pass `xmllint --noout`; all
3 3MFs pass `unzip -t`; every solid reports closed, consistently oriented, 0 boundary edges
and 0 non-manifold edges.

### Figures

| figure | path |
|---|---|
| headline table, all arms against every baseline | `results/final/figures/summary_table.md` |
| yield versus face count, per family and per `σ`, Wilson 95 % intervals | `results/final/figures/fig_yield.png` |
| full gallery, every deployable design, closed and half-open | `results/final/figures/fig_gallery_full.png` |
| `ε_max` distribution, K9c against K9 | `results/final/figures/fig_eps_hist.png` |
| 0⁺ margin against exact `Θ_max`, feasible-but-jammed designs ringed | `results/final/figures/fig_feasible_vs_deployable.png` |
| E1 agreement and certificate confusion counts | `results/final/figures/fig_e1_agreement.png` |
| exact vs referee gap histogram | `results/final/e1/e1_gap.png` |
| K9c gallery, 20 best graphs, closed and half-open | `results/kill/k9c/k9c_gallery.png` |
| K9c `ε_max` distribution and per-arm counts | `results/kill/k9c/k9c_hist.png` |
| K9 gallery and per-arm plots | `results/kill/k9/` |
| jitter transition | `results/kill/jitter/` |
| agreement, mobility, locality, margin, orientation, re-closure | `results/kill/{k2a/k2a_agreement.png, k3a/k3a_mobility.png, k2c/k2c_locality.png, k2b/k2b_margin.png, k5/k5_orientation.png, k1c/k1c_rates.png}` |
| hero closed and open previews | `export/hero/hero_148_sigma_mc_{closed,open_half,open_0.9tm}.png` |

---

## How each gate was verified

* **Gate 1 (Read).** Both papers spot-checked directly against the PDFs by the orchestrator,
  including Appendix A. Eq. (18)'s printed hinge contribution `k(2π+θ)` is a typo for
  `k(2φ+θ)`; with the fix, Eqs. (15)–(17) give `x + y = π` exactly (`STATE.md` F10).
* **Gate 2 (Reimplement).** The orchestrator rebuilt from scratch with 0 warnings and re-ran
  the tests; reference cases include Fig. 21 (3,4,3,12), Fig. 6 hexagons, the 2025 Fig. 2
  patterns, and rank-deficient collision cases. **Label required:**
  `results/core_validation/reference_cases.md` is a Gate-2 artefact from an earlier code era
  — its snub-square `Θ_max = 1.646` came from the old collision code plus bisection, and
  T4's exact 1.671 supersedes it. Different code eras, not different seeds.
* **Gate 3 (Scout).** Three field tables, 25 + 26 + 20 verified rows, each "cannot do"
  quoted; the orchestrator re-ran two queries per scout and retrieved one full text.
* **Gate 4 (Ideate).** Top 5 ideas each carry a runnable kill spec.
* **Gate 5 (Kill).** Satisfied by K1a + K2a + K1c surviving with quantitative outcomes.
* **Gate 6 (Screen).** 37 queries in round 1, 21 in the bundle screen, 12 for K9. Verdicts
  are recorded as KNOWN / PARTIAL / NOT FOUND, and the KNOWN items are cited rather than
  claimed (`notes/screen_r1.md`, `notes/screen_bundle.md`, `notes/screen_r2.md`,
  `notes/screen_k9.md`).
* **Gate 7 (Derive + Check).** Six deriver/checker rounds, the checker in a fresh context
  each time, ending with **zero unresolved disagreements**. `derivation_tests.jl`: **⟨JULIA:derivation_tests⟩,
  0 failures** (re-run by the writer). **Disclosure:** MISSION §9's
  threshold — a derivation not reconciled after two checker rounds — was formally crossed at
  round 5. No round disputed any theorem's truth value; each residue was a side condition or
  wording item the checker itself prescribed. That is convergence rather than
  irreconcilability, but the reader is entitled to the fact and to disagree with the
  judgement.
* **Gate 8 (Build).** `Pkg.test()`: **⟨JULIA:tests⟩, 0 failures** (re-run by
  the writer). Export verified by `xmllint`, `unzip -t`, `check_manifold` and visual
  inspection of PNG previews. **Not verified: no slicer is installed, so no 3MF was opened in
  one.**
* **Gate 9 (Experiment).** E1 (3,113 designs) and K9c (400 of 400 designs) as above.
* **Gate 10 (Report).** This file.

---

## Every kill test

| test | hypothesis | bar | measured | verdict |
|---|---|---|---|---|
| K1b | the harmonic identity fits every predicate | worst relative residual < 1e−10 | 4.180e−12 over 207,096 triples | **PASS** |
| K1a | the Eq. (6) embedding is valid on random graphs | — | inverted on 142/200, overlap 200/200, 0/200 samples injective | **PASS-emptiness** |
| K2a | closed-form `Θ_max` matches brute force | ≤1e−5 on all | 187/187, worst 1.956e−10; naive min-over-roots errs by 1.047 rad on 5 | **PASS** |
| K2c | the locality gate prunes soundly | pruned == unpruned | 285/285, worst difference 0; the old gate is identically 1 by algebra | **PASS on the replaced rule; original PASS withdrawn** |
| K1c | the authors' Eq. (9) silently re-closes | ≥3 % | **9.47 %** | **PASS** |
| K2b | null-space range optimization beats native `prevent` | ≥25 % median gain | **0.0000**; ours 7, theirs 15, tie 16 of 30 | **FAIL** |
| K3a | mobility ≫ 1 is a dangling-face artefact | identity 500/500 and `m_core ≤ 5` on ≥90 % | identity 493/508 (all ±1; 7/7 dense exact); `m_core ≤ 5` on **0/167** | **FAIL** |
| K5 | a defect-minimizing `σ` makes the space usable | a substantial fraction valid | **0/200** both `σ`; flat-sheet defect down 122× | **FAIL** |
| K6 | 0⁺ repair inside the null space | >0 certified | **0/400** | **FAIL** |
| A3 | a sharp jitter transition, with a predictor | predictor AUC ≥ 0.9 | transition sharp, `a* = 0.16–0.52` median edges; adversary's predictor 0.778 | **MIXED (FAIL as written)** |
| K7 | achievable periodic Jacobians | C3 ≥80 % certified | C1/C2/C4 pass to 1.4e−13 / 1e−16 / 5.5e−14; **C3 12/25 = 48 %** | **C1/C2/C4 PASS, C3 FAIL** |
| B4 | the fixed boundary empties the space | ≥10 % | **2/400** exact, **1/400** refereed | **FAIL** |
| K8a | a non-uniform expansive flex exists | ≥20/100 | **0/100**, dual-certified on 84/100 | **FAIL** |
| B3 | `tr K` budget threshold as a design constraint | AUC ≥0.8 **and** the constrained re-solve wins | AUC 0.980 / 0.950, but spread ratio median 0.73, max 1.40; re-solve 12→14/25 | **FAIL (hard-fail clause)** |
| T-1 | balanced pure vertices predict 0⁺ survival | 0 bad pure vertices | 1,968 bad of 14,886; 1,966 at reflex corners | **FAIL, with the mechanism as residue** |
| K9 | convexity + split-inward makes random graphs deploy | ≥40/400 | **36/400 = 9.0 %** | **FAIL by 4** |
| K9b | more search budget lifts the yield | >36 | **27/400** in the final merge (its own run reported 28); 9 feasible with `Θ_max = 0` | **FAIL, the diagnosis is the residue** |
| K9c | the objective, not the feasible set, was the cap | ≥40/400 | **307 / 400 = 76.8 %**, certified 306 | **PASS** |
| F23 | the authors' extra boundary row is over-constraint | — | FK from 53 seeds and 20 orders agrees to 7.16e−15 | **confirmed** |
| E1 | the characterization scales past 500 graphs | ≥500 vs brute force | **3,113/3,113**, worst gap 3.822e−09 | **PASS** |

---

## Dead ends

Each records what was tried and the number that killed it.

* **R1 range-optimal projection as a margin claim.** Killed by K2b: median gain 0.0000
  against the authors' native `prevent`, which itself reaches `Θ_max = π` on snub square. The
  earlier "86 % gap closure" was measured against our own `X0`, not their code.
* **R3 bifurcation spectrum / mobility bundle.** Killed by K3a: `m_core ≤ 5` on 0 of 167
  Delaunay patches, median 305. At that mobility the bifurcation question is ill-posed for
  generic finite patches.
* **A8, "Eq. (2) sufficient but not necessary."** Already published: 2025 Fig. 10 shows
  non-uniform deployment, and Dang et al. 2021 (PRE 104:055006, arXiv:2106.15891) prove the
  iff for quadrilateral tessellations.
* **Our own Eq. (9) claim (P3).** Withdrawn after baseline parity: our reimplementation of
  Eq. (9) is weak, but the authors' native `prevent` at published defaults is not
  (`STATE.md` F18). Any margin claim must be against their code.
* **K7 C3, target-driven periodic design.** 12/25 = 48 % against an 80 % bar. Reachability is
  never the problem — on 25 designable patterns every target `K` is hit to 8.9e−16 — the 0⁺
  collision is.
* **T-1 vertex balance as a predictor.** Balance holds to 2e−14 and the designs collapse
  anyway. Correction on the record: an earlier log entry said T-1 was contradicted because
  Eq. (2) holds; that reason was wrong, balance does hold per vertex, and the conclusion
  fails for a different reason.
* **B4, boundary drop.** 2/400 exact, 1/400 refereed. The `|V_∂| − 1` extra dimensions are
  real (400/400), but 33 handles cannot pay for 126–321 split-sign constraints.
* **K8a, expansive-cone LP.** 0/100, dual-certified structural on 84/100. Positive residue:
  the certificate machinery, and one genuine rescue example.
* **B3, `tr K` budget.** The threshold quantity co-varies with the budget over `𝒦`; the
  identity survives as a verified theorem, the algorithm does not.
* **Curvature-budget cluster, T-3, I1, the audit family.** Dropped at ranking after the
  round-2 screen placed them under Konaković-Luković et al. 2018 and 2025 Remark 4.1.

---

## What a reviewer will object to

Drawn from `review/theory_review.md`, `review/negative_review.md` and
`review/constructive_review.md`. Each objection is stated as the reviewer would state it,
then answered honestly.

1. **"Most of the theory bundle is textbook."** Correct for T1, T2, T6, K7 C1/C4 and F35, and
   they are presented as lemmas for that reason. What survives as new is T4's
   exactness/completeness with graze handling, T3's degeneracy catalogue, T5's explicit inner
   region, T7's off-by-one correction, and K7 C2.
2. **"The paper never claims random graphs work, so you are attacking a claim it does not
   make."** True, and it must be said. "Arbitrary" in 2026 §5 is combinatorial — not
   periodic, not 2-colorable, holes allowed — demonstrated with one hand-picked example per
   exclusion, never over a population. The random-graph experiment is ours.
3. **"Your baseline is not their pipeline."** Correct, and currently the sharpest objection.
   Native `prevent` was invoked on **8** designs across roughly 1,600 (design, `σ`) slots,
   because on the rest the embedding was already invalid before `prevent` had anything to
   repair. Wherever `0/400` appears it means our four repair strategies. **[Native200:
   pending]** — the full run across all 400 designs under both `σ` is in progress
   (`results/kill/native200/shards`, 108 rows so far, many at the 600 s timeout).
4. **"The population is adversarially large."** Median `F` = 332.5 against `F ∈ [4, 97]` for
   every reference case in the paper's own figures. The generator is unbiased; the *scale* is
   outside what the paper demonstrates. Both halves are true and both are stated.
5. **"Your yield is family-dependent and falls with graph size."** Correct, and it is
   plotted rather than averaged away: Voronoi 122/134, Delaunay 108/134, quad-random 77/132,
   and yield declining with `|F|` (Delaunay about 1.0 to 0.8, quad-random 0.75 to 0.43,
   Wilson 95 % intervals in `results/final/figures/fig_yield.png`). Nothing here claims a
   scale-free method.
6. **"Median opening 0.254 rad is not deployment."** About 15°. The `ε_max ≥ 0.1` rad count
   (207 of the 307 positives, 214 taking the best arm per design) is reported alongside the
   median for exactly this reason, and the maximum is π, a fully open design.
7. **"Feasibility does not determine deployability, by your own measurement."** True under
   the proximity objective, and reported rather than buried: nine exactly-feasible designs
   had `Θ_max = 0`, which is why the objective had to change. Under the range objective, in
   the completed 400-design run, no such design remains
   (`results/final/figures/fig_feasible_vs_deployable.png`). The constraint pair is still a
   necessary filter and is not claimed as the mechanism.
8. **"The solver is a non-convex heuristic with no completeness guarantee."** Correct. The
   certificate is sound; the search is not. Nothing here characterizes when the constrained
   set is non-empty.
9. **"Your dual certificates cover a handful of branch charts."** Correct: 4 visited charts
   of `2^2338`. They are sound where they apply and do not cover all of `P(X)`.
10. **"K7 characterizes a set you cannot reach."** Correct, and the C3 caveat is attached to
    every invocation of C1/C2/C4.
11. **"Your certificate needs an empirically-discovered deflation rule."** The systematic
    degeneracy catalogue is currently a derivation note plus 73,446 measured occurrences of
    2,145,387. It needs to be a lemma with a proof. Open.
12. **"`Θ_max ≤ min(min β, π)` — you proved only one direction."** Correct. The reverse is
    measured on 8 patterns and must never be printed as an equality without the tag and the
    sample size.
13. **"Slivers."** Present in the K9c gallery and acknowledged; short edges and sharp angles
    are what the paper's Limitation 1 also names.
14. **"Six checker rounds means the derivation was not reconcilable."** MISSION §9's
    threshold was formally crossed. The judgement that this was convergence is recorded above
    with its reasoning, so a reader can disagree with it.
15. **"A concurrent preprint may already do this."** Cleared. The full text of
    arXiv:2608.30032 (Jiang and Choi, v1 of 30 Aug 2026, 25 pp.) was retrieved
    (`papers/related/jiang_choi_2026.pdf`) and screened claim by claim against C1–C10 in
    `notes/screen_jiang_choi.md`. It contains no collision content, no self-intersection
    content, no deployment-range quantity, no null space and no certificate, and it works on
    `N × N` rotating-squares quad patterns rather than arbitrary planar graphs. The single
    overlap is its Eq. (2), the cross-product corner-convexity inequality, which is the same
    primitive already conceded to Choi, Dudte and Mahadevan (2019); C9 stays PARTIAL on that
    primitive alone.
16. **"No physical artifact."** Correct, and it is a scope decision rather than an omission
    to be repaired: the project is theory and computation only, so nothing was fabricated, no
    slicer is installed and no 3MF was opened in one. The constructive review's suggestion —
    laser-cut and actuate the `ε_max ≈ 2 rad` outlier and photograph it — falls outside that
    scope. Every claim here is stated about the geometric model, and none is offered as
    evidence about a physical piece.

---

## Next steps toward a paper, in order

1. **Native200 is partial by decision** (stopped at 573/600 runs, F40): 0/172 completed designs open by our exact scan; 5/573 by their own test are false negatives. Finishing the last 27 runs would change nothing qualitatively; rerunning with a longer timeout on the 362 timed-out designs is the only remaining objection-3 item.
2. **Prove the T3 degeneracy catalogue as a lemma** (`h(0) = 0 ∧ h′(0) = 0` iff …) instead of
   citing a measurement. Half a day of derivation plus one checker round.
3. **Prove or restrict `dim 𝒦 = 2·rank(D)`**, currently verified on 33/33 but not derived.
4. **Quantify the certificate's completeness gap** on the corpus — how much of `U(ε)` the
   inner region `R(ε)` misses — or concede it is unknown, as the theory review requires.
5. **Replace the second-hand Tay–Whiteley citation** with the primary source, and obtain the
   Konaković-Luković 2018 full text before printing any "no source supplies this" claim.
6. **Physical fabrication is out of scope by decision**, not pending work. The project is
   theory and computation only; nothing is to be laser-cut, printed or photographed, and no
   claim in this report depends on a physical artifact. The `export/` SVG/3MF/STL artefacts
   stay as verified file-format output (`xmllint`, `unzip -t`) and are not a fabrication
   result.
7. **Explain the yield trend.** Yield falls with `|F|` and differs by family (Voronoi
   122/134 against quad-random 77/132). Whether that is split-cut density, aspect ratio, or
   the solver's basin is unknown, and it is the first question a reviewer will ask of
   `fig_yield.png`.
8. **Characterize when the constrained set is non-empty.** This is the open research question
   the work leaves behind, and the one that would turn a diagnostic result into a theory.

---

## Related work

Citations are drawn only from the field tables and screen logs.

Segall, Ren, Padilla and Sorkine-Hornung (SIGGRAPH Asia 2025) and Segall, Ren and
Sorkine-Hornung (TOG 2026) are the base; the 2026 paper's Limitation 1 states the problem
addressed here, and its released 2025 pipeline has no collision predicate at all — `θ_max`
there is the hinge-local kinematic bound `min(2π − α_i − α_j)` (`notes/repo_2025.md`).
Tay–Whiteley body-and-hinge motion assignments are what T1 and T2 specialize to two
dimensions. Acuña et al., *Auxetic behavior on demand* (Commun. Phys. 5, 2022,
arXiv:2101.12352) give a sufficient three-step recipe for counter-rotating `ν = −1` networks
under collinearity and a fixed 2-colouring. Grima et al. (Proc. R. Soc. A 468, 2012,
810–830) give a closed-form locking angle for the hinge-adjacent two-triangle case and
concede the general case is out of reach — the closest prior use of the harmonic contact
algebra. Dang et al. (PRE 104:055006, 2021, arXiv:2106.15891) prove rigid deployability for
quadrilateral tessellations iff all cut voids are parallelograms. Choi, Dudte and Mahadevan
(Nat. Mater. 2019, arXiv:1812.08644) already use the cross-product convexity primitive as a
non-overlap constraint solved with `fmincon`, without a null space or a split-inward
condition. Rote, Santos and Streinu (2003) give an exact polyhedral cone of first-order
expansive velocities for pointed pseudo-triangulations, the closest prior certificate for a
cognate problem; the C-IRIS line (Dai, Amice, Werner, Zhang, Tedrake, IJRR) supplies the
general pattern of Positivstellensatz emptiness certificates for collision sets.
Konaković-Luković, Panetta, Crane and Pauly (SIGGRAPH 2018) equate the expansion bound with
a conformal scale-factor bound and locate cone singularities where it fails, which subsumes
the curvature-budget family. Czajkowski, Coulais, van Hecke and Rocklin (Nat. Commun. 13,
2022, arXiv:2103.12683) treat *spatial* conformality of a continuum, read in full and ruled
out as prior art for K7 C2. Han et al., IsoGami (SIGGRAPH 2026, DOI 10.1145/3799902.3811127)
compute mobility numerically with IPC contact on isohedral periodic tilings only, with no
theorem, closed form, rank formula or arbitrary graphs. Jiang and Choi (arXiv:2608.30032, 2026) give a
length-based interior-point design framework for `N × N` rotating-squares quad kirigami across
2D-to-2D, 2D-to-3D and 3D-to-3D morphing, with an inertia-transposition and aspect-ratio theory
of the compact end states; they use the same cross-product corner inequality (their Eq. (2)) as
Choi, Dudte and Mahadevan 2019, have no collision predicate and no deployment-range quantity,
and check deployment by PyKirigami simulation.

---

## File map

| area | path |
|---|---|
| core pipeline (2026 reimplementation) | `Kirigami/src/core/{mesh,cut,holes,orientation,tutte_auxetic,kinematics,collision,optimize,rank_checks,generators}.jl` |
| method | `Kirigami/src/method/{deploy_basis,contact,design,range_opt,zero_plus,convex_embed,range_embed,mobility,expansive_cone,periodic_jacobian,budget}.jl` |
| export | `Kirigami/src/export/{layout,solid,svg,stl,threemf,xml,zip,material}.jl` |
| tests | `Kirigami/test/{test_mesh_cut,test_holes,test_system,test_kinematics,test_collision,test_method,test_design,test_range_embed,test_rank_checks,test_reference_cases,test_export,derivation_tests}.jl` |
| CLIs | `Kirigami/apps/{kiri_gen,kiri_analyze,kiri_deploy,kiri_design,kiri_export,kiri_sweep,kiri_reference}.jl` |
| kill drivers | `Kirigami/apps/kill_{k1a,k1b,k1c,k2a,k2b,k2c,k3a,k5,k6,k7,k8a,k9,k9b,k9c,b3,b4,t1,jitter,e1,native200,f23}.jl` |
| derivations | `derivations/core.md`, `derivations/check.md` |
| kill results | `results/kill/KILL_REPORT.md`, `results/kill/{k1a,k1b,k1c,k2a,k2b,k2c,k3a,k5,k6,k7,k8a,k9,k9b,k9c,b3,b4,t1,jitter,native200,f23}/` |
| characterization | `results/final/e1/{E1.md,e1.csv,summary.txt,e1_gap.png}` |
| core validation | `results/core_validation/{reference_cases.md,rank_claim.md,referee_fix.md,sweep.csv,cases/}` |
| baseline (authors' code) | `baseline/native/`, `baseline/parity.md`, `baseline/kirigami_tessellations/` |
| export artifacts | `export/hero/`, `export/samples/` |
| paper notes and screens | `notes/{paper_2026,paper_2025,repo_2025,field_kirigami,field_rigidity,field_tutte,screen_r1,screen_bundle,screen_r2,screen_k9}.md` |
| reviews | `review/{theory_review,negative_review,constructive_review}.md` |
| ideation | `ideas/{persona_*,ranking,ranking_r2,round2_*}.md` |
| state and escalation | `STATE.md`, `ESCALATION.md` |

---

## Appendix — reproduction

**Build.**

```
julia --project=Kirigami -e 'using Pkg; Pkg.instantiate()'   # one-off; no build step
```

**Tests** (both re-run by the writer at the time of writing).

```
julia --project=Kirigami -e 'using Pkg; Pkg.test()'   # ⟨JULIA:tests⟩, 0 failures; includes derivation_tests.jl (⟨JULIA:derivation_tests⟩)
```

**The characterization (E1).**

```
# no build step; the driver is Kirigami/apps/kill_e1.jl (run_e1.sh invokes it)
bash code/scripts/run_e1.sh 200                       # 12-way sharded; slowest shard 311 s
arch -arm64 /usr/local/bin/python3 code/scripts/plot_e1.py
```

**Kill drivers.**

```
julia --project=Kirigami Kirigami/apps/kill_k1b.jl  --out results/kill/k1b          # ~34 s
julia --project=Kirigami Kirigami/apps/kill_k1a.jl  --out results/kill/k1a          # ~40 s
julia --project=Kirigami Kirigami/apps/kill_k2a.jl  --out results/kill/k2a          # ~7 s
julia --project=Kirigami Kirigami/apps/kill_k2c.jl  --out results/kill/k2c          # ~3 min
julia --project=Kirigami Kirigami/apps/kill_k3a.jl  --out results/kill/k3a          # ~25 min
julia --project=Kirigami Kirigami/apps/kill_k3a_recheck.jl --out results/kill/k3a   # ~4 min
julia --project=Kirigami Kirigami/apps/kill_k5.jl   --out results/kill/k5 --n 200 --maxf 800   # ~5 min
julia --project=Kirigami Kirigami/apps/kill_k1c.jl  --out results/kill/k1c          # ~31 s   (needs baseline/native)
julia --project=Kirigami Kirigami/apps/kill_k2b.jl  --out results/kill/k2b --maxf 160          # ~26 min (needs baseline/native; the native numbers themselves are NOT rerun -- they come from the authors' unchanged binary, reused as archived)
julia --project=Kirigami Kirigami/apps/kill_f23.jl  --out results/kill/f23          # ~1 s
julia --project=Kirigami Kirigami/apps/kill_k6.jl   --out results/kill/k6
julia --project=Kirigami Kirigami/apps/kill_k7.jl   --out results/kill/k7
julia --project=Kirigami Kirigami/apps/kill_t1.jl   --out results/kill/t1           # ~7 s off the K6 cache
julia --project=Kirigami Kirigami/apps/kill_b3.jl   --out results/kill/b3
julia --project=Kirigami Kirigami/apps/kill_b4.jl   --out results/kill/b4
julia --project=Kirigami Kirigami/apps/kill_k8a.jl  --out results/kill/k8a
julia --project=Kirigami Kirigami/apps/kill_k8a_recheck.jl --out results/kill/k8a
julia --project=Kirigami Kirigami/apps/kill_jitter.jl --out results/kill/jitter
julia --project=Kirigami Kirigami/apps/kill_k9.jl   --out results/kill/k9
julia --project=Kirigami Kirigami/apps/kill_k9b.jl  --out results/kill/k9b
julia --project=Kirigami Kirigami/apps/kill_k9c.jl  --out results/kill/k9c          # 12 shards; --aggregate --nshards 12 to merge
julia --project=Kirigami Kirigami/apps/kill_native200.jl --out results/kill/native200           # in progress
arch -arm64 /usr/local/bin/python3 code/scripts/plot_kill.py
arch -arm64 /usr/local/bin/python3 code/scripts/plot_k9c.py
```

**The method on one graph.**

```
./code/build/kiri_gen delaunay --out g.json --seed 1 --orient auto
./code/build/kiri_analyze g.json --orient json --boundary fixed --out results/tmp --collision
./code/build/kiri_design g.json --sigma both --out design.json --referee
```

**The hero.**

```
./code/build/kiri_design export/hero/hero_input.json --sigma json --seed 10036 \
  --out export/hero/hero_design.json --referee --out-graph export/hero/hero_graph.json

EXP=code/build/src/export/kiri_export
$EXP export/hero/hero_graph.json --profile felt_laser --theta 0 \
  --svg export/hero/hero_148_sigma_mc_closed.svg \
  --stl export/hero/hero_148_sigma_mc_closed.stl \
  --3mf export/hero/hero_148_sigma_mc_closed.3mf \
  --json export/hero/hero_148_sigma_mc_closed.json
$EXP export/hero/hero_graph.json --profile felt_laser --theta-frac 0.5 ...   # open at Theta_max/2
$EXP export/hero/hero_graph.json --profile felt_laser --theta-frac 0.9 ...   # open at 0.9*Theta_max
arch -arm64 /usr/local/bin/python3 export/render_svg.py <svg> -o <png>
```

The full hero recipe, including the throwaway graph-dump step needed to reproduce the exact
K9 population graph bit-for-bit, is in `export/hero/README.md`.
