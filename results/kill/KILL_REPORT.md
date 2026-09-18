# Phase 5 — kill experiments (Experimenter)

Written incrementally, one experiment at a time. Every number below is measured by the
Julia drivers in `Kirigami/apps/kill_*.jl` against `Kirigami/src/core` and `Kirigami/src/method`;
nothing is quoted from memory. Machine: Apple silicon (arm64), native Julia as described in
`docs/NUMERICS.md`. Populations are deterministic functions of the graph id
(`data/corpus/ (frozen populations, see data/corpus/README.md)`), so every experiment sees the same graphs.

Order run: K3a, K1b, K1a, K2c, K2b, K1c, K2a — as instructed. This report was completed
by a second Experimenter after the first was killed mid-run; K1b, K1a, K2c and K3a are
the first Experimenter's runs, re-read and re-verified here, and K2a, K2b, K1c, the K3a
recheck and the F23 check are new.

## Corrections applied, and which of them changed a verdict

These arrived in four rounds, from `derivations/core.md` (rounds 1-3), `derivations/check.md`
(the Checker, 29 033 + 59 982 assertions) and the orchestrator. Every one is implemented in
`Kirigami/src/method/` and covered by a unit test in `Kirigami/test/test_method.jl`.

1. **`Θ_max` is not the minimum over harmonic roots** (T4.2″). A root can be a *graze* — a
   tangency without a crossing. `hexagons_auto` grazes at `π/3 = 1.047198` and does not
   overlap until `2π/3 = 2.094395`, confirmed here with the graze handling in place.
   `Θ_max` is computed by the T4.2″ interval scan (`exact_theta_max_overlap`): take the
   complete contact set `C(X)`, walk the gaps, return the left endpoint of the first gap
   whose midpoint has an interior overlap. `i* = 0` is the zero-range case and needs no
   special probe. **Changed K2a from FAIL to PASS.**

2. **The `τ = 0` deflation, in three classes** (T5.2b.2, `check.md` R2.5). With
   `A = p − q`, `B = 2r`, `C = p + q` and `g(τ) = C + Bτ + Aτ²`, `τ = tan(θ/2)`:
   `C ≠ 0` is class 1; `C = 0, B ≠ 0` is class 2, a **simple** root at `θ = 0`;
   `C = 0, B = 0, A ≠ 0` is class 3, `h = p(1 − cos θ)`, a **double** root at `θ = 0` and,
   by Lemma T5.1e, constant sign on `(0, π)` — never a contact, atom `true`;
   `A = B = C = 0` is a permanent incidence, struck by an identity test on the
   coefficients before any classification. `C = 0` is not a degeneracy: it holds
   *identically on the whole shape space* for every permanent incidence and every
   split-edge duplicate pair. Without the deflation the flat-state root leaks in at `~1e−8`
   and every pattern with split cuts reports `Θ_max = 0`. **Measured here: the deflation
   removes 13 630 spurious roots in `(0, ε)` across K2a's 187 configurations**, and the
   deflated root list is exact on 60 000 synthetic harmonics of all three classes.

3. **The swept radius: the `√2` correction is withdrawn** (`check.md` D1, `core.md` round 2).
   In the face's own frame the per-face translation cancels and
   `S_u − ḡ_f = −σ_f J (C_u − ḡ_f)`, so the two columns are orthogonal and of **equal
   norm** and the trajectory is a **circle**: `max(‖x‖, ‖χ‖) = σ_max = ‖x_u − x̄_f‖`
   **exactly**. Applying `√2` there makes the pruning looser, not sounder. The broad phase
   now uses that exact radius with T4.5a's **moving-centroid** distance test. Effect:
   surviving pairs per face drop from a median of 22.23 to **12.69**, with `Θ_max`
   unchanged. The `√2` factor is real only for the *absolute* trajectory (a genuine
   ellipse; 72 of 98 copies violate `max(‖C‖,‖S‖)`, worst ratio 1.2137), which no pruning
   uses. **My earlier report of this as a two-frame subtlety was right in substance, but
   the prescription for the pruning frame is withdrawn outright, not qualified.**

4. **The flat-centroid pruning is unsound in principle** (`check.md` D8). Deployment
   contracts centroid distances while face radii stay fixed, so a flat-centroid test can
   discard a pair whose faces approach (19 842 such pairs measured). The moving test is now
   the default; the flat variant is kept only to reproduce that measurement.

5. **`Θ_max = min(min_e β_e, π)` on split-free patterns**, not `min_e β_e` (`check.md` D3).
   The triangle tiling has `min β = 4π/3 > π`. **Measured here: 16/16 split-free
   configurations satisfy the capped identity to `8.882e−16`, and the cap bites on 5 of
   them.**

6. **The `τ` chart is singular at `θ = π`** (`check.md` D10). `harmonic_roots` uses the
   amplitude/phase form (`atan2` + `acos`), not the `τ`-quadratic, so this is already
   satisfied; the deflated class-2 branch handles `A = 0` explicitly as a root at `θ = π`
   (`h = r sin θ` there). All searches are capped at `π`.

7. **The validity certificate** (T5.2b.0/T5.2b′, `check.md` R3.3):

   > `VALID(ε)` = `POS` ∧ `NOOVERLAP(θ₁ = ε/2)` ∧ `NOROOT(ε)`

   `POS` is all face signed areas `> 0` at `θ = 0`; `NOOVERLAP` is one exact
   polygon–polygon test at the **single interior angle** `θ₁ = ε/2`; `NOROOT` is no
   admissible **deflated** root in `(0, ε)` over the candidate list minus the
   identically-zero harmonics. **The `θ = 0` overlap test is withdrawn**: at `θ = 0` the two
   copies of every split edge coincide and adjacent faces touch along whole shared edges, so
   it is degenerate, not merely different. **Measured here: the certificate implies
   `Θ_max ≥ ε` on 173/173 configurations where it holds, 0 violations** — and it is exact on
   this population: 173/187 configurations actually have `Θ_max ≥ ε`, and 173 are certified.
   (Before the F32 fix — `NOROOT` lacked the two interval tests — it certified only 65/187;
   see §A3 and `results/kill/jitter/cert_diagnosis.md`.)

8. **H-LOC is refuted, not unproved** (`check.md` D2). The words `O(n)`, "certified active
   set" and "locality theorem" are withdrawn. See K2c.

9. **The referee's grid must be fine** (found here, not in a derivation). At 180 samples the
   bisection steps clean over a **genuine overlap window of width `2.9e−3` rad** at
   `θ = 0.9467` on `trunc_square_R20_s2` and reports `2.025` where the closed form says
   `0.9467`; a `1e−4` fine scan (`Kirigami/apps/dbg_k2a_out.jl`) confirms the closed form. All
   referee bisections here use 4000 grid points.

### Which certificate each number uses

| number | predicate |
|---|---|
| K2a's certificate columns, and the T5.2b′ soundness check | `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT`, `ε = 0.006` |
| K2a's `Θ_max`, and every `Θ_max` reported anywhere | T4.2″ interval scan, refereed by bisection on the exact overlap predicate, 4000-point grid. The predicate was rewritten after F34 (`results/core_validation/referee_fix.md`): it no longer shrinks the faces, and its tolerance argument is a relative degeneracy tolerance, not a shrink, so the referee now agrees at `1e−12`, `1e−9` and `1e−6` alike |
| K1a's `p_valid` and "X0 injective" | `n_inv = 0` **and** no overlap at `θ = 0` — the **withdrawn** predicate. K1a is quoted as it was measured and its numbers are not restated under the new certificate |
| K2b's "valid ladder point" and "certified validity" | no overlap at `θ = 0` — the **withdrawn** predicate (K2b ran before the certificate was defined) |
| K1c's "certified design" | no overlap at `θ = 0` — the **withdrawn** predicate |
| K5's "overlap-free X0" | no overlap at `θ = 0`, because that is the predicate the K5 PASS rule names and the one K1a used, so the two are comparable. K5 **also** reports the proper certificate, and the difference is the whole finding |

## Verdict table

| experiment | verdict | key number |
|---|---|---|
| K3a | **FAIL** (both clauses, but see the recheck) | identity 505/508 as measured; every exception is a sparse-rank estimate on a graph above the dense limit, and 7/7 of the rechecked ids hold **exactly** with a dense rank. Second clause fails outright: `m_core ≤ 5` on **0 of 167** Delaunay patches (median `m_core = 305`). `σ ∈ ker A` on **296/296 at X0**, 0/212 at `X_ini` |
| K1b | **PASS** | worst relative fit residual `1.450e−12` over 207 096 triples (rule: `< 1e−10`) |
| K1a | **PASS-emptiness** | `n_inv(X0) > 0` on 142/200; `p_valid = 0` on 134/200; an injective sample on **1/200** graphs, `Θ_max > 0` on **0/200** |
| K2c | **PASS on the replaced rule; the original PASS is WITHDRAWN** | pruned == unpruned exact `Θ_max` on **285/285**, worst difference **0**. The old gate (`slope ≤ 0.3`) is vacuous: the statistic is identically 1 by algebra. H-LOC refuted |
| K2b | **FAIL** | median relative gain over the authors' native baseline **0.0000** on 30 live designs (rule: `≥ 25 %`); ours better on 7, native better on 15, tie on 16 |
| K1c | **PASS** | corrected re-closure rate **9.47 %** on the authors' native Eq. (9) (rule: `≥ 3 %`). The rule as written reports **65.09 %**; with the crossing clause applied at its own reference range it reports **exactly 0**, because the rule is self-contradictory |
| K2a | **PASS** (with T4.2″) | `\|Θ_max − bisection\| ≤ 1e−5` on **187/187**, worst `1.817e−10`; min over roots gets **182/187**, worst error `1.047` rad. Certificate implies `Θ_max ≥ ε` on **173/173** after the F32 fix (65/65 before) |
| K5 | **FAIL** on the final rule | `σ_def` gives a valid `X0` by the certificate (`ε = 0.3`) on **0 / 200** — as does `σ_mc`. It does fix the flat sheet: defect down **122×**, `POS` from 58/200 to **186/200**, projection distance down **4.7×**, and 62.5 % pass the **withdrawn** `θ = 0` test |
| K6 | **FAIL** | certified `Θ_max > 0` on **2/400** under the primary repair and **15/200** `σ_mc` (7.5 %) under the split-only secondary — **17/400** taking the best of primary / secondary / ladder / stage 2 per design, 15 of them with `ε_max ≥ 0.1` rad — every one by leaving the projection's neighbourhood; 0⁺ split-duplicate collision is what the exact scan names on every other graph |
| A3 | **MIXED** (FAIL as literally written) | exact `Θ_max > 0` transition sharp at `a* = 0.16–0.52` median edges on the 4 split-bearing tilings (4 of 8 authored tilings are split-free: shape space is a point); as first run the certificate `NOROOT(ε=0.006)` rejected **1622/1859** rows with `Θ_max ≥ 1` rad (87.3 % false negatives → F32), after the F32 fix it certifies **1856/1859**; adversary's predictor ρ AUC 0.778 (rule ≥ 0.9) |
| K7 | **C1/C2/C4 PASS, C3 FAIL** | `J(θ) = cos(θ/2)I + sin(θ/2)K` exact to 7.3e−14, `K` affine to 4.8e−14, `dim 𝒦 = 2·rank D` on **33/33**; conformal at every θ (≤1e−15, 25/25); `θ_c = 2·atan2(tr K, 1 − det K)` to 5.6e−14 on 39/39; target-driven certified design **13/25 = 52 %** (rule 80 %) — reachability is total, obstruction is the 0⁺ collision |
| K9 | **FAIL (by 4 designs), first non-zero near the projection** | convexity + split-inward constrained embedding gives exact `Θ_max > 0` on **36/400 = 9.0 %** (refereed **36 / 36** after the F34 predicate fix; it was 32 at the old shrink `1e−12`), certified `ε_max` median 0.25 rad, max 2.0, `≥ 0.1` on 33/36; rule was ≥ 10 % (≥ 40). Every baseline on the same 400 designs is **0** or near it: Eq. (6) alone 0, `σ_mc` 0, `σ_def` 0, K6 repairs 2 (primary) / 15 (split-only, far from the projection), B4 2, K8a 0, the authors' `prevent` 0. Convexity alone: 118 feasible, **0** deployable |
| K9b | **FAIL** | the solver push gives **28/400** refereed (7.0 %), BELOW K9's 36; best over the two solver configurations **37/400**, 34 with `ε_max ≥ 0.1` rad, on 34/200 graphs. Only the `δ` sweep paid (`δ = 1e−2·med²` wins 15 of 28). The 9 K9 positives it loses are all exactly margin-feasible and still `Θ_max = 0` ⇒ the cap is basin selection, not search budget |
| K9c | **PASS** | replacing the PROXIMITY objective by the 0⁺ MARGIN inside the same barriers gives refereed exact `Θ_max > 0` on **312 of 400 designs (78.0 %)**, against a bar of ≥ 40; **310** certified, **192** with `ε_max ≥ 0.1` rad (**201** best-of-3), `ε_max` median **0.162** (best-of-3 **0.172**), max **π** (optimiser-path dependent at the ±2 % level; an earlier pass gave 307 / 306 / 214). The two proximity arms on the same 400 give **36** (K9 settings) and **27** (`δ = 1e−2`). All **9** of K9b's lost-but-feasible designs are recovered. Soundness `ε_max ≤ bisection` **400/400**. Run COMPLETE: 400 of 400 |
| Native200 | **FAIL (partial coverage)** | the authors' full native pipeline (their coloring, Eq. (6), Eq. (9), FK collision test), run unconditionally on the K9 population, merged over the 600 s pass and the 3 600 s reruns (`native200_final.csv`): exact `Θ_max > 0` on **0 / 195** completed cells; their own test says `Θ > 0` on **27 / 195**, all 27 false negatives of the F24/K1c mechanism. Coverage 574/600 dispatched (26 never dispatched): 195 completed, 364 timed out, 15 crashed. Caveat: median `\|F\|` 202 (completed) vs 332.5 (population) vs ≤ 97 (paper's own examples) |
| F23 | **confirmed** | forward kinematics from all 53 seed faces and 20 BFS orders agree to `6.24e−15` |

---

## K3a — the 2-core mobility identity

**PASS rule, copied from `ideas/ranking.md` Sec. 4:** "Check `dim ker A = |F \ core2(Gamma)| +
dim ker A|core2` as an integer identity, ranks by `ColPivHouseholderQR` with threshold
`1e-10 * ||A||_2`. Report `m_core = dim ker A|core2 - c(core2)` on the random Delaunay
patches where U4 reported `m = 19-26`. PASS if the identity holds 500/500 **and**
`m_core <= 5` on `>= 90%` of those patches (confirming the adversary's prediction that
U4's headline is a dangling-face artefact). If the identity fails on any graph, `A` is
misassembled and B11–B16 all stop until it is fixed."

**Measured** (`results/kill/k3a/{k3a.csv,summary.txt}`, 508 graphs, `F ∈ [10, 4918]`,
wall 44.5 s):

| quantity | value |
|---|---|
| identity holds | 505 / 508 |
| identity at `θ = 0.3·Θ_max` (the extra probe asked for) | 8 / 8 |
| `σ ∈ ker A` (Eq. (2)) | 296 / 508 |
| Delaunay patches with `m_core ≤ 5` | **0 / 167 (0.00 %)** |
| Delaunay median `m_full` / median `m_core` | 322 / 305 |
| Delaunay `m_full` range / `m_core` range | [39, 1875] / [36, 1675] |

**The exceptions are a rank-estimator artifact, not a misassembled `A`.** The three
violators of the current run (`delaunay_88`, `F = 1644`; `delaunay_247`, `F = 3103`;
`quad_random_266`, `F = 1079`) all sit above `mobility.jl`'s `dense_limit = 700`, where
`matrix_rank` switches from the dense pivoted QR to a sparse QR whose rank estimate is
unreliable near the threshold; which of the large graphs trip it depends on the sparse
factorisation and is not reproducible across linear-algebra backends (an earlier pass of the
same driver reported 15 such exceptions, every one off by exactly ±1).
`Kirigami/apps/kill_k3a_recheck.jl` regenerates the ids of that earlier list with K3a's exact
configuration (the Eq. (6) projection where `N ≤ 1400`, else `X_ini` — `rank(A)` depends on
the pin positions, so this had to match) and recomputes both ranks densely
(`results/kill/k3a/{k3a_recheck.csv,recheck_summary.txt}`):

| id | kind | F | X | sparse-QR defect (recheck) | dense defect | smallest kept `σ_i / σ_0` |
|---|---|--:|---|--:|--:|--:|
| 257 | quad_random | 1396 | X0 | 0 | **0** | 5.913e−05 |
| 266 | quad_random | 1079 | X0 | −2 | **0** | 1.886e−04 |
| 320 | quad_random | 1358 | X0 | 0 | **0** | 6.678e−04 |
| 425 | quad_random | 1428 | X0 | 0 | **0** | 5.869e−05 |
| 445 | delaunay | 1138 | X0 | 0 | **0** | 1.698e−03 |
| 467 | quad_random | 911 | X0 | 0 | **0** | 6.131e−05 |
| 494 | quad_random | 879 | X0 | 0 | **0** | 2.632e−04 |

7 of 7 rechecked: the identity holds **exactly** with a dense rank, including
`quad_random_266`, the one current violator small enough to recheck. The singular gap at the
cut is `5.9e−5` to `1.7e−3` — three to eight orders above the `1e−10` threshold, so the rank
is not a borderline call and the sparse estimate is simply wrong there. The other
violators (`F` from 1644 to 4902) are too large for a dense QR and were not rechecked; I
cannot claim them, only that they show the same signature.

**Verdict: FAIL on the rule as written, on both clauses, and the two failures are of
completely different kinds.**

* Clause 1 (identity 500/500) fails at 505/508 **as measured**, but the evidence says the
  identity is true and the rank estimator is not. `A` is not misassembled, so B11–B16 are
  not blocked — provided every mobility number above 700 faces is computed densely. That
  is a build item.
* Clause 2 (`m_core ≤ 5` on ≥ 90 % of Delaunay patches) fails **outright, at 0 %**, and
  this is a real result: **the adversary's prediction is refuted.** Deleting the dangling
  faces removes almost nothing — median mobility goes from 322 to 305, about 5 %. The
  finite mobility of a free-boundary Delaunay patch is *not* a dangling-face artefact; it
  is intrinsic to the 2-core.
* **`σ ∈ ker A` is 296/508 only because two populations were mixed, and the orchestrator
  is right that only `X0` counts.** Broken out by configuration: `σ ∈ ker A` on
  **296 of 296** graphs where `A` was evaluated at the solved Eq. (6) projection `X0`
  (worst residual `3.595e−12`), and on **0 of 212** where the solve was skipped (`N > 1400`)
  and `A` was built at `X_ini` instead, with residuals up to `6.373e+02`. That is the
  expected answer, not a defect: `σ ∈ ker A` **is** Eq. (2), so it holds exactly when the
  configuration is uniformly deployable, and `X_ini` is not. All 8 reference cases are in
  the `X0` group and all 8 hold. The correct statement is **296/296 at `X0`**; the
  `296/508` figure should not be quoted.
* Note also that the measured mobility is nowhere near `STATE.md` U4's headline of
  `m = 19–26`. On this population (100–5000 faces) the median is 322 and the range is
  [39, 1875]. U4's numbers were measured on a different and much smaller population, and
  should not be quoted as a general figure.

Artifacts: `results/kill/k3a/k3a.csv`, `summary.txt`, `k3a_recheck.csv`,
`recheck_summary.txt`. Driver: `Kirigami/apps/kill_k3a.jl`, `Kirigami/apps/kill_k3a_recheck.jl`.

---

## K1b — the harmonic identity (the gate)

**PASS rule, copied from `ideas/ranking.md`:** "200 graphs from K1a's set, 20 random `X in X`
each. For 50 random ordered (vertex, edge) triples per graph, sample `deploy()` at 200
angles on `(0, min(pi, theta_max))`, least-squares fit the orientation determinant to
`p + q cos(theta) + r sin(theta)`. PASS if the relative residual is `< 1e-10` on every
triple of every graph. FAIL on any exceedance — and then B1/B3/B17/B20 all die together, so
nothing else may be started before this returns. Also assert face signed areas are constant
in `theta` to `1e-12`."

**Measured** (208 graphs, `F ∈ [10, 793]`, 207 096 triples, wall 42.2 s):

| quantity | value | rule |
|---|---|---|
| worst relative LS residual of `p + q cos + r sin` | **1.450e−12** (`voronoi_129`) | `< 1e−10` ✓ |
| triples with residual / max\|det\| ≥ 1e−10 | 0 | ✓ |
| worst residual / (\|AB\| \|AP\|) | 1.425e−12 (`voronoi_129`) | ✓ |
| worst \|closed form − `deploy()`\| / (\|AB\| \|AP\|) | 1.442e−12 | |
| worst \|fitted − closed-form\| coefficients | 3.141e−13 | |
| worst relative \|`deploy` − (cos C + sin S)\| | 1.854e−15 | |
| worst relative face signed-area drift | **1.002e−13** (`quad_random_5`) | `1e−12` ✓ |
| triples skipped as numerically coincident | 14 (both \|AB\| and \|AP\| below `1e−6` of the diameter) | |

**Verdict: PASS.** The trig-linear form `Y(θ) = cos(θ/2) C + sin(θ/2) S` and the harmonic
`p + q cos θ + r sin θ` structure of every orientation determinant are confirmed to
`1.5e−12` over 207 096 triples, and faces are rigid to `1.0e−13`. Everything downstream
(B1, B3, B17, B20, the whole closed-form contact calculus) is cleared to proceed.

Artifacts: `results/kill/k1b/{k1b.csv,summary.txt}`. Driver: `Kirigami/apps/kill_k1b.jl`.

---

## K1a — is the baseline really invalid, and is the shape space really usable?

**PASS rule, copied from `ideas/ranking.md`:** "**PASS-algorithm** if `n_inv(X0) > 0` on
`>= 150/200` graphs and `p_valid > 0` on `>= 20` graphs. **PASS-emptiness** if
`n_inv(X0) > 0` on `>= 150/200` and `p_valid = 0` on `>= 20` graphs. **FAIL** (premise dead,
F17 was a sweep artefact) if `n_inv(X0) = 0` on `> 100/200`. Both PASS branches keep R1
alive; they select which paper it is."

**Measured** (200 graphs, `F ∈ [101, 793]`, `10^4` Gaussian samples each plus a
trust-region ladder, wall 60.6 s):

| quantity | value |
|---|---|
| `n_inv(X_ini) > 0` (asserted 0) | **0** ✓ |
| `n_inv(X0) > 0` | 142 / 200 |
| `n_inv(X0) = 0` | 58 / 200 |
| mean inverted fraction of `X0` | 0.0410 (max 0.1655) |
| `p_valid > 0` (fraction of samples with `n_inv = 0`) | 66 / 200 |
| `p_valid = 0` | **134 / 200** |
| any trust-region radius valid | 69 / 200 |
| `X_ini` injective (control) | **200 / 200** |
| **`X0` injective (no face overlap)** | **0 / 200** |
| **some sample injective** | **1 / 200** |
| **some sample with `Θ_max > 0`** | **0 / 200** |

**Verdict: PASS-emptiness, and by a wider margin than the rule contemplates.** The rule is
written in terms of `n_inv` — inverted faces — and on that measure the answer is already
`142/200` invalid and `p_valid = 0` on `134/200`, which is PASS-emptiness (`≥ 150` is missed
by 8, `≥ 20` is passed by 114). But `n_inv` is the wrong measure (correction 3): under the
**real** validity criterion, no face–face overlap in the flat state, the Eq. (6) projection
is invalid on **200/200**, not 142/200, and **of the two million samples drawn, a single one
on a single graph is a valid embedding — and it has `Θ_max = 0`** (which sample lands there
depends on the sampler's floating-point path; the 0/200 with `Θ_max > 0` does not). 58 graphs have every face positively oriented and *still* self-intersect
— which is exactly T5.1's point that positive orientation certifies nothing.

I am reporting this branch as it came out and did not retune the sampler. The consequence
for R1 is the one `ranking.md` anticipated: on this population the usable region `U(ε)` is
not reachable by isotropic sampling of the Tutte auxetic null space around `X0`, so R1 is
"the emptiness paper", not "the algorithm paper" — **for random planar graphs.** It is not
empty in general: the authored tilings supply 187 valid deployable configurations
(`deployable_population()` in `kill_common.jl (populations frozen in data/corpus/)`), which is exactly the population K2a and
K1c run on. The honest statement is that validity is a *strong* constraint that random
graphs miss, not that `U(ε)` is empty.

Artifacts: `results/kill/k1a/{k1a.csv,summary.txt}`. Driver: `Kirigami/apps/kill_k1a.jl`.

---

## K2c — the locality gate: the original PASS is withdrawn

**The rule as `ranking.md` writes it:** "PASS if the log-log slope [of `max_f rho_f / r_f`
against `log n`] is `<= 0.3` **and** pruned/unpruned `theta_max` agree to `1e-12` on every
graph."

**The first clause is vacuous and its PASS is withdrawn.** In the face's own frame — the
frame the spec names and `swept_discs()` implements — the algebra gives
`χ_u = −σ_f J x_u`, hence `‖χ_u‖ = ‖x_u‖` and `ρ_f = r_f` for **every face of every patch at
every size**. The statistic is identically 1, so its log-log slope is 0 for any input
whatsoever: the gate cannot fail and cannot test H-LOC. Measured `1.000000000001`. The
"K2c PASS, slope 0.0000" I reported earlier is an artefact and is retracted.

**The replaced rule, per the orchestrator:** (1) PASS iff pruned and unpruned exact `Θ_max`
agree to `1e−12` on every graph; (2) report the empirical growth of the candidate set with
`n`, claiming no `O(n)` theorem.

**Measured** (462 graphs, 285 with `F ≤ 800` for the `Θ_max` half; wall 105 s; broad phase
now the exact face-frame radius with the moving-centroid test):

| quantity | value |
|---|---|
| **pruned == unpruned exact `Θ_max` to `1e−12`** | **285 / 285**, worst difference **0.000e+00** |
| the same for the (unsound in principle) flat-centroid variant | 285 / 285, worst difference `0.000e+00` |
| worst \| ‖Y_pv − m_f(θ)‖ − r_f \| — exactness of the moving disc | 6.729e−13 |
| surviving pairs per face | median **12.69**, max 18.40 |
| log-log slope of pairs per face vs `log F` | 0.1695 (R² 0.32, n = 285) |

The tighter radius is worth a factor of 1.75: the median surviving pairs per face falls from
22.23 (with the loose `√2` radius) to **12.69**, with `Θ_max` unchanged on every graph.

**Verdict: PASS on the replaced rule. The broad phase is sound and it is measurably tight.
No `O(n)` claim is made or supported.**

### H-LOC is refuted

H-LOC would need `max_f ‖γ_f(θ) − x̄_f‖ ≤ κ · r_f` with `κ` independent of the patch — that
is what the packing argument for an `O(n)` active set needs. Measured on a growing square
tiling (the random population lives in a fixed 40×40 box, so its diameter barely varies and
a regression on it is meaningless):

| F | patch diameter | `max_f ‖γ_f − x̄_f‖ / r_f` |
|--:|--:|--:|
| 4 | 2.828 | 1.414 |
| 16 | 5.657 | 2.316 |
| 36 | 8.485 | 3.924 |
| 64 | 11.314 | 5.560 |
| 100 | 14.142 | 7.205 |
| 144 | 16.971 | 8.854 |
| 196 | 19.799 | 10.505 |
| 256 | 22.627 | 12.157 |
| 289 | 24.042 | 12.984 |
| 361 | 26.870 | 14.637 |

A straight line: `drift/r = a + 0.5656 · diameter`. On the fixed-box random population the
drift ranges from 25.8 to 105.0 circumradii. So the swept region of a face is not contained
in any fixed multiple of its flat circumdisc, the packing count does not close, and **H-LOC
is refuted rather than merely unproved.** The words "`O(n)`", "certified active set" and
"locality theorem" must not appear in the paper. What survives is an *exact* broad phase
(sound with no hypothesis) that is empirically tight — a good algorithm, not a theorem.

Artifacts: `results/kill/k2c/{k2c.csv,k2c_drift.csv,summary.txt,k2c_locality.png}`. Driver:
`Kirigami/apps/kill_k2c.jl`.

---

## K2a — exact `Θ_max` versus bisection

**PASS rule, copied from `ideas/ranking.md`:** "Compute `theta_max` two ways: (i) closed-form
smallest positive root over all ordered (vertex, edge) pairs, each root filtered by the two
closed-form interval tests `0 <= dot(y_b - y_a, y_v - y_a) <= ||y_b - y_a||^2` evaluated at
that root; (ii) `theta_max` from `collision.jl` (grid scan + bisection, faces shrunk by
relative `1e-6` per deviation 9). PASS if `|theta_max_exact - theta_max_bisect| <= 1e-5` rad
on every graph. FAIL on any exceedance. Secondary output, reported either way: whether the
binding pair is a split-edge pair."

**Population deviation, recorded.** Not 292 random graphs: K1a shows no random graph in this
project's population has an embedded Eq. (6) projection, and on a non-deployable embedding
`deploy()` disagrees with itself at the shared hinge pins, so the two methods would be
comparing different families rather than two ways of measuring one. The population is
`deployable_population()` — seven authored tilings at five clip radii, each with `X0` and up
to eight still-embedded shape-space samples, **187 configurations**, every one both uniformly
deployable and embedded, with the basis error asserted per member.

**Measured** (187 configurations, wall 40.6 s):

| method | agreement with bisection (`1e−12`, 4000-point grid) | worst error |
|---|--:|--:|
| **T4.2″ interval scan (the corrected closed form)** | **187 / 187 = 100.00 %** | **1.817e−10** |
| min over roots (the rule as written) | 182 / 187 | **1.047 rad** (`hexagons_R40`) |
| min over roots, shipped `1e−6` shrink | 182 / 187 | 1.047 rad |

The referee is tolerance-independent since the F34 predicate fix
(`results/core_validation/referee_fix.md`): the min-over-roots row reads **182 / 187** at
`1e−6`, `1e−9` and `1e−12` alike (before the fix the shipped `1e−6` shrink gave 149 / 187).
The certificate columns do not depend on the referee.

Secondary numbers:

| quantity | value |
|---|---|
| configurations with `Θ_max = 0` (`i* = 0`, immediate penetration) | 14 |
| configurations where min-over-roots < `Θ_max` (a **graze**) | 5 |
| mean \|C(X)\| / mean intervals probed | 40.1 / 1.9 |
| flat state embedded (overlap predicate at `1e−12`) | 187 / 187 |
| binding edge is a split duplicate / hinge / no contact | 130 / 42 / 10 |
| binding faces not adjacent in `M` | **1** |

**Verdict: PASS with the T4.2″ scan; FAIL with the rule exactly as `ranking.md` writes it.**
The correction is worth 5 configurations and, on `hexagons`, `1.047` rad — a third of the
entire deployment range. The mechanism is the one T4.2″ predicts: the two duplicates of a
split edge are congruent translates, so when they become collinear the contact is
automatically vertex-to-**vertex**, the faces touch at a point and separate. This is not a
rare degeneracy; it holds *identically* on the whole shape space of a symmetric tiling, and
symmetric tilings are what both papers ship.

The scan is cheap: the mean candidate set has 40.1 angles and the scan stops after 1.9
overlap probes. Reproducing the Deriver's table exactly (`Kirigami/apps/dbg_t422.jl`):

| pattern | \|C\| | `θ₁` (first contact) | `Θ_max` (T4.2″) | bisection | `min β` |
|---|--:|--:|--:|--:|--:|
| rotating_squares | 1 | — | 3.141593 | 3.141593 | 3.141593 |
| triangles_alternating | 0 | — | 3.141593 | 3.141593 | 4.188790 |
| kagome_3636 | 1 | 3.141593 | 3.141593 | 3.141593 | 3.141593 |
| **hexagons_auto** | 3 | **1.047198** | **2.094395** | 2.094395 | 2.094395 |
| truncated_square_488 | 12 | 2.356194 | 2.356194 | 2.356194 | 2.356194 |
| snub_square_33434 | 18 | 1.646136 | 1.646136 | 1.646136 | 3.383898 |
| tiling_3_4_3_12 | 10 | 2.380636 | 2.380636 | 2.380636 | 2.380636 |

**The corrections, measured on this population.**

| quantity | value |
|---|---|
| split-free configurations with `Θ_max = min(min β, π)` to `1e−12` | **16 / 16**, worst `8.882e−16` |
| of those, patterns where the `π` cap bites (`min β > π`) | **5** |
| candidate harmonics by class 1 / 2 / 3 + identically zero | 728 783 / 98 376 / 631 + 0 |
| **spurious roots in `(0, ε)` removed by the `τ = 0` deflation** | **13 630** |
| certificate `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT` holds (`ε = 0.006`) | 173 / 187 (65 / 187 before the F32 fix) |
| ... of which `POS` / `NOOVERLAP` / `NOROOT` alone | 187 / 173 / 187 |
| **T5.2b′: certificate holds ⟹ bisection `Θ_max ≥ ε`** | **173 / 173, 0 violations** |
| configurations that actually have `Θ_max ≥ ε` | 173 / 187 |

Three things to read off. First, the class-2 population is large (98 376 of 827 790
candidates, 12 %) and every one of them has its `θ = 0` root deflated away; without that,
13 630 of them would have been reported as contacts inside `(0, ε)` and the usable region
would be empty. Second, T5.2b′ is confirmed: wherever the certificate holds, the independent
bisection agrees that `Θ_max ≥ ε`, with no exceptions. Third, as originally run the certificate was a strict
**inner** approximation — 65 certified against 173 actually valid — because `NOROOT` dropped
the two interval tests (collinearity with the edge's infinite line counted as a contact).
That was F32; with the interval tests restored (`contact.jl`, §A3) the certificate holds on
173 against 173 valid, 0 violations. `NOOVERLAP` (173/187) is the binding
clause, `POS` is 187/187.

**One binding pair is not adjacent in `M`.** 2026 §4.5's stated empirical assumption — that
the first contact is between faces sharing a split edge — is false on at least one of these
187 configurations. It is also true 130/187 of the time, so the assumption is a good
heuristic and a bad theorem.

Artifacts: `results/kill/k2a/{k2a.csv,summary.txt}`. Drivers: `Kirigami/apps/kill_k2a.jl`,
`Kirigami/apps/dbg_t422.jl`, `Kirigami/apps/dbg_k2a_out.jl`.

---

## K1c — Eq. (9)'s false-negative rate

**PASS rule, copied from `ideas/ranking.md`:** "Run `optimize_collision_sweep` (the documented
gamma ladder `{1e-4, ..., 10}`, best kept) to get `Y_9`. For every split-edge pair evaluate
`P, R` in closed form and count pairs with `P < 0`, `R > 0`, and
`theta* = 2 arctan(-R/P) < theta_max(Y_9)`. Cross-check every reported `theta*` against
`collision.jl`'s bisection to `<= 1e-6` rad. PASS if the fraction of gamma-ladder-certified
designs with at least one re-closure below `theta_max` is `>= 3%`. FAIL below 3% — B17 dies,
R1 survives without its §5."

**Correction applied** (T5.3, and the orchestrator's instruction). A root of the separation
harmonic is a re-closure only if the duplicates meet **as segments** and the faces actually
**cross** there. Both extra clauses are implemented:

```
   r = <d, du> > 0                        Eq. (9) certifies this, and nothing else
   p = -sigma_f det(d, du) < 0            so tau* = -r/p > 0 exists
   theta* = 2 arctan(-r/p) in range
   INTERVAL clause  0 <= dot <= |.|^2 at theta*      (segments, not lines)
   CROSSING clause  the two faces' interiors overlap just after theta* and not just before
```

The T5.2 closed forms were verified against the generic harmonic on every split edge:
`p + q = 0` to **0.000e+00 exactly** (structural, as T5.2 predicts) and
`|closed − generic| / scale` to **1.899e−08** worst. Recovering the derivation's `du` from
this basis needs `du = −½ J (S(a″) − S(a′))`; with `du` taken naively as `S(a″) − S(a′)` the
check reads `5.0`, which is a rotation, not a disagreement.

**Measured** (171 deployable embedded designs with split cuts, wall 61.1 s):

| `Y_9` source | certified designs | rule **as written** (`θ* < Θ_max`, no clauses) | clauses applied, **same** range `Θ_max` | corrected: `θ* < min β` | corrected: `θ*` **is** the binding contact |
|---|--:|--:|--:|--:|--:|
| authors' native `prevent` (published defaults) | 169 | **110 (65.09 %)** | **0 (0.00 %)** | **16 (9.47 %)** | 12 (7.10 %) |
| our `optimize_collision_sweep` (γ ladder) | 163 | 51 (31.29 %) | **0 (0.00 %)** | 46 (28.22 %) | 12 (7.36 %) |

Worst `|θ* − per-pair bisection|` over qualifying re-closures: **9.390e−11** (rule: `1e−6`).

The native column is a fixed function of the authors' binary; the `optimize_collision_sweep`
column is the endpoint of a gradient descent with a nine-value `γ` ladder and is
optimiser-path dependent (an earlier pass of the same driver ended at 115 certified designs,
40 / 9 / 32 in the last three columns); the verdict does not depend on which endpoint is used.

**The rule as literally written in `ranking.md` is unsatisfiable, and that is the finding.**
Column 4 is exactly **zero on both sources, and it has to be**: the crossing clause says the
faces' interiors overlap just after `θ*`, and `Θ_max` is *by definition* the first angle at
which any two interiors overlap. So `θ* < Θ_max` and "θ* is a genuine crossing" are
contradictory. The `65.09 %` in column 3 is entirely made of roots that are **not**
collisions — grazes, or collinearities where the duplicates never meet as segments. K1c as
specified does not over-count by a factor; it counts a set that is empty once the clauses
are applied, and 110/169 of the authors' designs would have been reported as having a
false negative on the strength of nothing.

Two references therefore remain, and both are reported:

* **`θ* < min_e β_e`** — did Eq. (9) leave a re-closure inside the range the pattern would
  otherwise reach? This is the question B17 is actually about, and it is well posed because
  `min β` is a kinematic bound, computed without reference to collisions. **9.47 %** on the
  authors' code, **28.22 %** on ours.
* **`θ* == Θ_max`** — is the re-closure the binding contact, i.e. is a split-edge
  re-closure what actually stops the deployment? **7.10 %** and **7.36 %**.

**Verdict: PASS on both sources and on both well-posed references (7.10 %, 9.47 %, 7.36 %,
28.22 %; rule `≥ 3 %`).** B17 lives: 2026 Eq. (9) certifies `r = h′(0) > 0` and is blind to
`p` by construction, and on about one design in ten of the authors' own output that blindness
costs range. The `65.09 %` headline the uncorrected experiment would have produced is not a
result and must not be quoted.

**The closed forms.** T5.2 was verified against the generic harmonic on every split edge:
`p + q = 0` to **0.000e+00 exactly** (structural, as T5.2 predicts) and
`|closed − generic| / scale` to **1.899e−08** worst. Recovering the derivation's `du` from
this basis needs `du = −½ J (S(a″) − S(a′))`; taken naively as `S(a″) − S(a′)` the check
reads `5.0`, which is a rotation, not a disagreement. The residual `1.9e−08` is carried
entirely by split edges that are nearly collapsed; the median over all edges is `2.0e−15`.

**A side finding about our reimplementation.** Our `optimize_collision_sweep` collapses **96
split edges to zero length** across this population; the authors' native `prevent` collapses
**0**. That is why only 163 of 171 of our outputs are valid flat embeddings against 169 of
171 of theirs, and it is an independent reason (beyond F18) not to quote our Eq. (9)
reimplementation as a baseline for anything.

Artifacts: `results/kill/k1c/{k1c.csv,summary.txt}`. Driver: `Kirigami/apps/kill_k1c.jl`.

---


## K5 — is the ORIENTATION why the shape space is empty?

Added after K1a, at the orchestrator's request, as the Critic's R5 kill.

**The question.** K1a found that with `σ` from Eq. (1) — a max-cut relaxation on the dual
graph — not one of two million shape-space samples over 200 random graphs is a valid
embedding. Eq. (1) picks `σ` to maximise the number of hinge cuts and knows nothing about
whether the hole-closure system then has a solution near `X_ini`. K5 asks whether choosing
`σ` to minimise the **deployability defect** instead fixes that.

**Definitions.** `D(σ) = Σ_holes ‖Σ_{hinge e ∈ C} (x_target − x_source)‖²` evaluated at
`X_ini` — exactly Eq. (2)'s residual, so `D = 0` iff `X_ini` is already uniformly
deployable. `σ_def` is a greedy local search over face flips minimising `D`: single-face
flips and adjacent-pair flips, strictly improving moves only, any flip rejected that makes
`c(Γ) > 1` or detaches a face; started from `σ_mc` and 3 random `σ`, best kept, at most
`20·|F|` attempts per start, deterministic seed `7000 + id`.

**PASS rule (final):** `σ_def` gives a **valid** `X0` on `≥ 20 %` of graphs, where valid is
the final certificate `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT` with `ε = 0.3`.

**Measured** (the same 200 K1a graphs, `F ∈ [100, 800]`, wall 152 s; `c(Γ) = 1` held on every
accepted `σ`, and the search never failed to keep `Γ` connected):

| quantity | `σ_mc` (Eq. 1) | `σ_def` |
|---|--:|--:|
| **valid `X0` by the certificate, `ε = 0.3`** | **0 / 200** | **0 / 200** |
| ... of which `POS` alone | 58 | **186** |
| ... `NOOVERLAP(0.15)` alone | **0** | **0** |
| ... `NOROOT(0.3)` alone | **0** | **0** |
| **certified `Θ_max > 0`** | **0 / 200** | **0 / 200** |
| median `D(σ)` | 2.245e+03 | **1.840e+01** |
| median inverted faces of `X0` | 9.5 | **0** |
| median `‖X0 − X_ini‖_∞` / median edge length | 1.587 | **0.334** |
| median `\|E_split\|` | 126.5 | **321.5** |
| median `dim_null` | 126.5 | 321.5 |
| graphs with `D = 0` exactly (`X_ini` already deployable) | 0 | **16** |
| overlap-free `X0` at `θ = 0` (the **withdrawn** predicate) | 0 / 200 | **125 / 200** |

Correlation of `D` with `\|E_split\|`: **−0.6446** under `σ_mc`, **−0.0346** under `σ_def`.
Median flips accepted: 138.5.

**Verdict: FAIL on the final rule — 0 %, against a threshold of 20 %.** Under the earlier
form of the rule, which used the `θ = 0` overlap test and asked only about the graphs where
`σ_mc` fails, `σ_def` scores **62.5 %** and would PASS. Both are reported because the gap
between them *is* the result.

### What `σ_def` fixes, and what it does not

**It fixes the flat sheet, decisively.** The defect falls by a factor of **122** in the
median, the Eq. (6) projection has to move the mesh **4.7× less far** (0.334 median edge
lengths against 1.587), the median count of inverted faces goes from 9.5 to **zero**, and
`POS` goes from 58/200 to **186/200**. On 16 graphs the search finds a `σ` for which `X_ini`
is *already* uniformly deployable, so Eq. (6) has nothing to do. Eq. (1) really is choosing
a bad `σ`, and choosing `σ` against the defect instead is cheap and works.

**It does not touch the deployment.** `NOOVERLAP(0.15)` is false on **all 400** designs and
`Θ_max = 0` on all 400. That is not a single-probe artefact: a log-spaced ladder
(`1e−12 … 1.0`) finds an interior overlap at the **first or second rung on every graph that
overlaps at all** (`1e−10` at worst over the 125 `σ_def` graphs; with the pre-F34 predicate
it was `1e−12` on every one — `results/core_validation/referee_fix.md`), and a direct scan (`Kirigami/apps/dbg_k5.jl`) confirms overlap at
`1e−2`, `0.1` and `0.5` as well, where the shrink tolerance is not in question. The flat
sheet is a valid planar cut pattern; it self-intersects the instant the cuts open.

This is exactly the degeneracy that forced the certificate change. At `θ = 0` the two copies
of every split edge coincide and adjacent faces touch along whole shared edges, so a
`θ = 0` test cannot see the failure that matters here — which is why the earlier rule's 62.5 %
and the final rule's 0 % disagree so completely, and why the `θ = 0` predicate was withdrawn.

**Why, mechanically.** `σ_def` buys its low defect by more than **doubling** the split cuts
(median 126.5 → 321.5): a split edge imposes no closure constraint, so cutting more edges
fully apart makes Eq. (2) easy to satisfy. The hinge graph is then sparse, the structure is
barely a mechanism, and it collapses into itself immediately. The correlation of `D` with
`\|E_split\|` bears this out and is strongly non-monotone across the two families — `−0.64`
under `σ_mc`, `−0.03` under `σ_def` — confirming the optimizer persona's prediction that the
relationship is not monotone.

**The lesson.** A defect-only objective is wrong for the same reason max-cut is: **neither
one knows about collisions.** `σ` is a real and badly-chosen design variable — that much R5
gets right and it is worth saying in the paper — but no purely combinatorial objective on
`σ` will produce a deployable random graph. The objective has to carry a deployment term,
and the certificate of T5.2b′ is now cheap enough to be that term.

**The follow-on comparison was not run**, per the instruction gating it on K5 passing. It
would also have been vacuous: the range optimiser and the authors' native `prevent` would
both have started from designs with `Θ_max = 0`, so there is no range to improve and no
referee number to report.

`σ_def` orientations are saved as JSON graphs in `results/kill/k5/sigma/<kind>_<id>.json`
(our contract: `vertices`, `faces`, `orientation`), 200 files, for Builder-Method to reuse.

Artifacts: `results/kill/k5/{k5.csv,summary.txt,sigma/,k5_orientation.png}`. Driver:
`Kirigami/apps/kill_k5.jl`, diagnostic `Kirigami/apps/dbg_k5.jl`.

---

## K2b — the range-optimisation margin, against the authors' native code

**PASS rule as amended by the orchestrator:** "Baseline for K2b is the AUTHORS' native code:
`baseline/native/build/tuttekiri_cli`, subcommand `prevent` (their Eq. 9) with published
defaults plus their parameter ladder, best kept; referee both by OUR collision bisection
(their own `theta_max` routine is buggy, `STATE.md` F24). Keep our `collision.jl` Eq.(9) as
a secondary column. PASS if the median relative gain is `>= 25%` over the best baseline with
certified validity on `>= 20` graphs; report FAIL honestly otherwise."

**`ranking.md`'s original rule, for the record:** "PASS if the median relative gain over the
best-gamma baseline across the 30 graphs is `>= 25%` **and** snub square closes `>= 30%` of
its measured gap (from 1.646 toward `min beta = 3.384`)."

### How it was run

* **Baseline.** `tuttekiri_cli prevent --sweep`, which starts from *their* Eq. (6) projection
  and runs *their* `opt::prevent_intersections`. Every ladder point's optimised pattern is
  dumped (a `--dumpdir` option added to our CLI in `baseline/native/src/main.cpp`; nothing in
  `baseline/tuttekiri/` was touched), plus their `X0` as the `T = 0` fallback. Ladder: the
  full 84-point grid `{0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0} × {1, 10, 50, 100} × {0, 0.1, 1.0}`
  on the four reference cases; the 36-point sub-ladder `{0, 0.1, 0.3, 1.0} × {1, 10, 100} ×
  {0, 0.1, 1.0}` elsewhere, which contains the published default `(0.1, 10, 0.1)` and every
  corner of the full grid. Their own `theta_max` is never used (F24: `merge_close_verts`
  fuses hinge duplicates before their overlap test).
* **Referee.** Our bisection on the true polygon overlap (`1e-12`), 4000-point grid,
  applied identically to every ladder point, to our Eq. (9) output and to ours. A ladder
  point counts only if it is a **valid flat embedding** (no face–face overlap at `θ = 0`).
* **Ours.** `maximize_range` — softmin of the closed-form first-contact roots over the null
  space with a signed-area log barrier, active set refreshed every 10 iterations — with a
  `T = 0` fallback so it is never reported worse than the embedding it started from.
* **Population.** 38 designs: the 4 reference cases with split cuts, 8 random graphs, and 26
  members of `deployable_population()` with `dim_null ≥ 2`. Wall 1533 s.

### Two deviations, both recorded

1. **Cost.** The native `prevent` is the whole cost and it is steep in the null-space
   dimension: measured on voronoi patches, one ladder point costs `0.95 s` at 60 faces
   (`k = 73`), `2.3 s` at 100 (`k = 118`), `8.8 s` at 150 (`k = 168`), `27 s` at 200
   (`k = 222`), `63 s` at 300 (`k = 329`). `ranking.md`'s "26 random graphs, 100–800 faces"
   with a ladder is a ~10 h run. Random graphs are capped at 160 faces.
2. **The random graphs make the experiment vacuous, so most slots went elsewhere.** On
   `voronoi_0/3`, `delaunay_1/4`, `quad_random_2/5` and two more, **0 of 37** native ladder
   points is a valid flat embedding (one graph had 4 valid points, all of which penetrate
   immediately). Baseline, our Eq. (9) and our optimiser all return `Θ_max = 0` and the
   relative gain is `0/0`. This is K1a's emptiness result reappearing, now with the authors'
   own code: it is not that our projection is bad, it is that **nothing in the shape space of
   a random planar graph is a valid embedding**. Eight such graphs are kept to document it;
   the remaining slots went to `deployable_population()`, the only population on which "the
   margin over the baseline" has a numerical answer.

### Measured

| quantity | value | rule |
|---|--:|---|
| designs | 38 (30 live, 8 vacuous) | |
| **median relative gain over the native baseline** | **0.0000** | `≥ 0.25` ✗ |
| median gain restricted to live designs (baseline > 0) | **0.0000** (n = 30) | ✗ |
| designs with certified validity and range > 0 | **30 / 38** | `≥ 20` ✓ |
| ours better / native better / tie | **7 / 15 / 16** | |
| gain quartiles over the 30 live designs | −0.305, −0.018, 0.000, 0.000, +0.106 | |
| snub square (reference case) | native `3.1416`, ours `3.1416`, `min β = 3.3839` | gap closed **0.0 %**, rule `≥ 30 %` ✗ |

**Verdict: FAIL, on both clauses, and not narrowly.** The median gain is zero, not 25 %.

### What the numbers actually say

* **The authors' code is strong, and F18 was right to withdraw the Critic's headline.** On
  `snub_square_33434` their ladder reaches `Θ_max = π`, the top of the search range, certified
  by our own bisection. There is no headroom left to win. The `86.1 %` gap closure quoted by
  the first Experimenter was measured against **our** `X0` (1.6461 → π), not against their
  code; against their code the gap is already closed and our number is a tie.
* **Our optimiser mostly does not move.** On 16 of 30 live designs it returns exactly the
  refereed range of the embedding it started from, i.e. the `T = 0` fallback wins. The
  closed-form softmin objective is correct (K2a certifies the roots to `2e−10`) but the
  optimiser is not finding anything with it on these patterns.
* **Where we lose, we lose badly.** `hexagons_R20_s5` and `_s6`: the native reaches `3.0148`
  and `2.9530`, ours stays at `2.0944` — a 30 % deficit. `2.0944 = 2π/3` is exactly the
  graze-limited hexagon value of K2a, so our optimiser is sitting still at the local barrier
  while their penalty walks the design out of the hexagon geometry entirely.
* **Where we win, we win little.** Best gain `+10.6 %` on `tiling_3_4_3_12`; 7 wins in total.
* **Our Eq. (9) reimplementation is the weakest of the three** and should not be cited: it
  returns `Θ_max = 0` on every `snub_square_R20*` design (invalid flat embedding) and on two
  `trunc_square` designs, consistent with the 96 collapsed split edges K1c found.

The consequence for R1 is direct: R1's algorithmic half — "range-optimal embeddings beating
Eq. (6) + Eq. (9)" — is **not supported by this experiment**. R1's *diagnostic* half survives
intact and is stronger than before (K1a, K2a, K1c, and the vacuity finding above), but the
margin claim must be dropped or the optimiser rebuilt.

Artifacts: `results/kill/k2b/{k2b.csv,summary.txt,k2b_margin.png}`. Driver:
`Kirigami/apps/kill_k2b.jl`; CLI change in `baseline/native/src/main.cpp` (`--dumpdir`,
`--short`).

---

## F23 — is the authors' extra boundary-component row an over-constraint?

**The question** (`STATE.md` F23): the authors' code adds a constraint row for split-forest
components that touch the mesh boundary; we drop them as notches. Their `X0` satisfies our
system 8/8; our `X0` has residual 1.73 (`hexagons_auto`) / 1.00 (`snub_square`) in theirs;
our null space is one dimension larger on those two. The orchestrator's reasoning is that
the faces around a boundary-touching component form a *path* in `Γ` — the exterior breaks
the cycle — so no closure constraint is needed, and their extra row only shrinks their space.
The check asked for: deploy **our** `X0` by forward kinematics from two different seed faces
and confirm all `M′` vertex positions agree to `< 1e−9` and all hinges open by `θ`.

**Measured** (`Kirigami/apps/kill_f23.jl`), stronger than asked: **every** seed face, not two,
plus 20 randomised BFS visit orders, at four angles.

| pattern | F | n_split | H | rank(L) | dim_null | Eq. (2) residual of our `X0` |
|---|--:|--:|--:|--:|--:|--:|
| `hexagons_auto` | 10 | 5 | 4 | 4 | 5 | 2.483e−16 |
| `snub_square_33434` | 53 | 13 | 12 | 12 | 13 | 1.351e−15 |

| pattern | θ | seeds | orders | worst seed disagreement | worst order disagreement | worst hinge-angle error | worst split-duplicate parallel error |
|---|--:|--:|--:|--:|--:|--:|--:|
| hexagons | 0.30 | 10 | 20 | 2.245e−15 | 1.778e−15 | 3.886e−16 | 3.377e−16 |
| hexagons | 0.70 | 10 | 20 | 2.667e−15 | 8.882e−16 | 4.441e−16 | 3.331e−16 |
| hexagons | 1.20 | 10 | 20 | 1.807e−15 | 8.882e−16 | 6.661e−16 | 2.483e−16 |
| hexagons | 2.00 | 10 | 20 | 1.691e−15 | 8.889e−16 | 6.661e−16 | 4.578e−16 |
| snub_square | 0.30 | 53 | 20 | 3.560e−15 | 1.986e−15 | 8.327e−16 | 5.097e−16 |
| snub_square | 0.70 | 53 | 20 | 4.572e−15 | 8.882e−16 | 1.554e−15 | 6.184e−16 |
| snub_square | 1.20 | 53 | 20 | 5.348e−15 | 2.220e−15 | 2.442e−15 | 4.647e−16 |
| snub_square | 2.00 | 53 | 20 | 6.242e−15 | 1.332e−15 | 3.109e−15 | 1.768e−15 |

Worst over both patterns: seed/order disagreement **6.242e−15**, hinge-angle error
**3.109e−15**, split-duplicate parallel error **1.768e−15** — all fourteen orders inside the
`1e−9` the check asks for. (Different seeds place the structure in different rigid frames,
which is not a disagreement about shape; the comparison is after the best rigid alignment,
`rigid_align_residual`.)

**Answer: yes, the authors' extra boundary-component row is an over-constraint.** Our `X0`
already closes every hole (`Eq. (2)` residual `1e−15`), and its forward kinematics is
path-independent from every one of the 53 seed faces and every BFS order tried, with every
hinge opening by exactly `θ` and every split-edge duplicate pair staying parallel and equal
(Remark A.4). There is nothing left for an extra closure equation to enforce. Their row is
harmless in the sense that their own `X0` still satisfies our system — it only removes
admissible designs, one dimension of shape space on each of these two patterns.

Artifacts: `results/kill/f23/{f23.csv,summary.txt}`. Driver: `Kirigami/apps/kill_f23.jl`.

---

## The three orchestrator rank checks, restated

Folded in from `results/core_validation/rank_claim.md` as instructed. Final numbers:

1. **`H == |E_hinge| − |F| + c(Γ)`** — holds on **50 / 50** sweep graphs and **8 / 8**
   reference cases, with `H` counting only all-interior preimages. Counting the notches in
   breaks it (`H_all == pred` on **0 / 50**), so notches are *not* bounded faces of `Γ`.
   This confirms `STATE.md` U1. On a boundary-free patch `Γ` lives on a surface of Euler
   characteristic 0 and the count overshoots by exactly one:
   `H == |E_hinge| − |F| + c(Γ) − 1`, verified on the 4×4 and 6×4 square tori and the 4×4
   triangle torus. Eq. (1) auto-orientation gives `c(Γ) == 1` on 50 / 50.
2. **`L = R·D`** — exact, max absolute difference **0**, on **50 / 50** sweep graphs and all
   8 reference cases, with `D` the `|E_hinge| × N` signed hinge incidence (+1 target, −1
   source) and `R` the `H × |E_hinge|` 0/1 preimage-ownership matrix. Also
   `rank(L) == H − dim Z` on **27 / 27** (dense only).
3. **Periodic `rank(L) = H − 1`** — on the three boundary-free tori (`torus_squares_4x4`,
   `torus_squares_6x4`, `torus_triangles_4x4`): `1ᵀL = 0` exactly, `rank(L) = H − 1`,
   `dim Z = 1`, out-harmonic residual `3.3e−16 / 7.2e−16 / 5.6e−16`. So 2026 §4.4's "the
   number of independent equations equals the number of holes" is **false as stated for
   periodic patterns**. With a boundary the naive row-sum identity fails on 50 / 50; the
   corrected form `(1ᵀL)_v = I(v)·indeg(v) − Σ_{v→w} I(w)` holds on 50 / 50 and 8 / 8.

**Prose/table inconsistency resolved.** The Critic flagged that `rank_claim.md`'s prose and
its tables disagreed about check 2c (every `y ∈ Z` gives an out-harmonic `g`). The table
credited it `27 / 27` while the prose said it was vacuous. The prose was right and the table
was misleading: `dim Z == 0` on all 27 dense sweep graphs **and** on all 8 reference cases,
so zero residuals were evaluated — the sweep's own line "out-harmonic residuals actually
evaluated: 0" says so. I have edited `results/core_validation/rank_claim.md` to mark row 2c
`27 / 27 vacuously`, to note that the reference cases are vacuous too (the prose had claimed
they were not), and to rewrite the closing paragraph so it no longer says checks 2a, 2b and
2c all "held on every graph measured". Check 2c has real content on three patches only —
the tori, where `dim Z == 1` and the residuals are `≤ 7.3e−16` — and that is the honest
strength of the evidence for the out-harmonic characterization.

---

## What this bundle now looks like

**Supported, and stronger than the specs asked for.**
* The closed-form apparatus is real. The harmonic identity holds to `1.5e−12` over 207 096
  triples (K1b); exact `Θ_max` agrees with an independent bisection on **187/187**
  configurations to `2e−10` (K2a), after the graze correction and the `τ = 0` deflation.
* The validity certificate works and is proved to work: `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT`
  implies `Θ_max ≥ ε` on **173/173** configurations where it holds, with 0 violations (K2a,
  re-run after the F32 `NOROOT` fix; 65/65 before it).
* The broad phase is exact and sound and empirically tight: pruned and unpruned `Θ_max`
  agree to `1e−12` on **285/285** graphs, median 12.7 surviving pairs per face (K2c).
* 2026 Eq. (6) has no validity certificate, on **200/200** random graphs, and the authors'
  own code confirms it: **0 of 37** ladder points on random graphs is valid (K1a, K2b).
* Eq. (9) is blind to the second-order coefficient by construction, and on about one design
  in ten of the authors' own output that blindness costs range (K1c).
* The orientation is a real and unexploited design variable: choosing `σ` against the
  deployability defect instead of max-cut cuts the defect by **122×**, the projection
  distance by **4.7×**, and produces a valid flat sheet on **62.5 %** of the graphs where
  Eq. (1) fails (K5).
* Their extra boundary-component constraint row is an over-constraint (F23).

**Not supported.**
* **The margin.** Our range optimiser does not beat the authors' native `prevent` — median
  gain 0.0000 over 30 designs, 7 wins to 15 losses (K2b).
* **`O(n)` anything.** H-LOC is refuted: the centroid drift grows linearly with the patch
  diameter, from 1.4 to 14.6 circumradii over diameters 2.8 to 26.9 (K2c). "Certified active
  set", "locality theorem" and every `O(n)` claim are withdrawn. The broad phase is a good
  algorithm, not a theorem.
* **The K2c PASS as originally reported.** Its gate statistic is identically 1 by algebra
  and cannot fail. Retracted.
* **The dangling-face story.** `m_core ≤ 5` on **0 of 167** Delaunay patches; `m_core` grows
  roughly linearly with the patch (K3a).
* **A deployable random graph.** Neither orientation produces one: `Θ_max = 0` on all 400
  designs in K5, and on 200/200 in K1a. The obstruction is not the projection and not the
  orientation objective alone — it is that nothing in the pipeline optimises against
  collisions before the deployment starts.

**Blocked on a build item (orchestrator decision, recorded).** Every mobility rank above 700
faces goes through a sparse QR, whose rank estimate is wrong on 3 of 508 graphs in the
current run (15 in an earlier pass) where the singular gap is `1e−5` to `1e−3` (K3a recheck).
**No rank above 700 faces is to be quoted** until `mobility.jl`'s `matrix_rank` is fixed or
the dense path is forced; a dense recheck on a subsample is accepted as sufficient evidence
for now, and that is what `Kirigami/apps/kill_k3a_recheck.jl` provides (7 of 7 rechecked ids
small enough for a dense QR have the identity restored exactly). The larger violators, `F`
from 1644 to 4902, are not claimed either way. Every rank-derived number in K3a's table above `F = 700` inherits this
caveat, including the median `m_full` and `m_core` and their ranges.

## K6 — 0⁺ repair of the shape space: FAIL

Added after K5, as the last attempt to save the random-graph design space before pivoting
to periodic patterns.

**The question.** K1a and K5 both end at `Θ_max = 0`, and the exact scan always names the
same culprit: at `θ = 0⁺` the two copies of a split cut translate **into** each other
instead of apart. K5's `σ_def` fixes the flat sheet and does not touch this. K6 asks
whether the jam is *repairable inside the shape space itself* — the null space is large
(median `dim_null` 126 under `σ_mc`, 321 under `σ_def`), so is there a point in it at which
every cut opens outward?

**The 0⁺ calculus** (`Kirigami/src/method/zero_plus.jl`, derived there in full). By T1.B
the two copies of a split edge `e` differ by a pure translation `sin(θ/2)·dS_e`, so the sign
`q_e := det(dS_e, d_e)` decides the side: `q_e > 0` the cut opens, `q_e ≤ 0` the copies move
into each other and `Θ_max = 0`. The same first-order argument for a **corner** — a copy of
vertex `v` entering the corner cone of a face at `v` — gives a margin `μ`. Both `q_e` and
`μ` are exact **quadratics** in the shape-space coefficients `t`, because `deploy()` is
linear in `X`. Nothing here is sampled.

**PASS rule:** certified `Θ_max > 0` on **≥ 20 %** of designs for at least one `σ` family,
under the same final certificate as K5 (`POS ∧ NOOVERLAP(ε/2) ∧ NOROOT`, `ε = 0.3`).

### How it was run

* **Population.** Exactly the 200 K1a/K5 graphs (`voronoi` 67, `delaunay` 67,
  `quad_random` 66; `F` from 101 to 793, median 332.5), each under both orientation rules —
  400 designs. `σ_mc` is Eq. (1); `σ_def` is read back from K5's saved
  `results/kill/k5/sigma/*.json`, so K6 studies the designs K5 measured rather than
  re-running the greedy search. 0 graphs skipped for a missing `σ_def`.
* **The repair.** `F(t) = Σ_split softplus(−q_e/s) + Σ_faces softplus(−a_f/s)
  [+ w_corner Σ softplus(−μ/s)] + w_prox‖X(t) − X0‖²/(N s) + λ‖t‖²`, `s = (median edge)²`,
  minimised by the project's own L-BFGS from `t = 0` and 3 Gaussian starts of scale
  `0.1·med`, `max_iter = 1200`, seed `6000 + 7·id`. **Feasible** is always decided by the
  exact constraint values, never by the objective.
* **Two variants, same designs and seeds.** *Primary* = split + vertex-edge with the
  proximity term, `w_corner = 0.05`, `w_prox = 1e3` (weights fixed on a 3-graph pilot:
  `w_corner = 1` lets the few thousand corner terms swamp the few hundred area terms and the
  minimiser walks ~100 median edges away with inverted faces). *Secondary* = split-only,
  `w_corner = w_prox = 0`. Extensions reported but not part of the rule: a `λ` ladder
  `{1e−6, 1e−3, 1e−1, 1}` and a stage 2 that hands a feasible point to `maximize_range`.
* **Sharding.** 12 disjoint shards by `gidx % 12`, run concurrently; the merged CSV is
  bit-identical to the concatenated shards. Wall 1778 s of solver time across the shards.
* **Native cap.** The authors' `prevent` from the same `X0` was run only on designs that are
  0⁺-feasible with `F ≤ 550`, at most 12 per family, 420 s timeout (180 s in the earlier
  pass — `run.log` records 5 designs killed by the alarm). No design qualified inside any
  12-graph shard, so the comparison lives in a separate **40-graph capped run**
  (`results/kill/k6/native/`).

### Measured — all from `results/kill/k6/k6.csv` (400 rows)

| quantity | `σ_mc` | `σ_def` |
|---|--:|--:|
| designs | 200 | 200 |
| median `\|E_split\|` = median `dim_null` | 126.5 | 321.5 |
| **inward split edges at `t = 0`** (q1 / median / q3 / max) | 31.8 / **65.5** / 119.8 / 265 | 58 / **103.5** / 192.5 / 340 |
| median inward *fraction* of split edges | **0.500** | **0.342** |
| designs with **zero** inward split edges at `t = 0` | **0 / 200** | **0 / 200** |
| inward corner incidences at `t = 0` (q1 / median / q3 / max) | 120.5 / **215.5** / 406.8 / 670 | 149 / **252.5** / 468.3 / 798 |
| designs with zero inward corners | **0 / 200** | **0 / 200** |
| median `min q` at `t = 0` | −77.81 | −130.43 |
| **0⁺-feasible, primary variant** | **0 / 200** | **2 / 200** (1.0 %) |
| **0⁺-feasible, secondary (split-only)** | **59 / 200** | **67 / 200** |
| 0⁺-feasible under the `λ` ladder | **0 / 200** | **2 / 200** |
| **certified `Θ_max > 0` (`eps_max > 0`), primary** | **0 / 200** | **2 / 200** (1.0 %) |
| certified `Θ_max > 0`, secondary (split-only) | **15 / 200** (7.5 %) | **0 / 200** |
| certified `Θ_max > 0`, ladder / stage 2 | **0 / 200** each | **2 / 200** each (the same two designs) |
| certificate at `ε = 0.3` (`POS ∧ NOOVERLAP ∧ NOROOT`) | **0 / 200** | **0 / 200** |
| ... of which `POS` alone | 0 | 2 | 
| bisection `Θ_max > 0` on the secondary points | 15 / 200 | 0 / 200 |
| secondary certified `ε_max`, median / max | 2.466 / π | — |
| binding contact at the primary failures | `split-inward` 200 | `split-inward` 188, `vertex-edge` 9, `inverted` 1 |

**The best of our four repair strategies is `15 / 200 = 7.5 %` on one `σ` family, against a
bar of 20 %; the authors' native `prevent`, run separately on the 8 of these designs that
reached a valid flat state, is `0 / 8` (the unconditional 200-graph run is §Native200).**
Which designs the L-BFGS repair lands on is optimiser-path dependent (an earlier pass of the
same driver reported the primary variant at 0 / 400 and the secondary at 0 / 400 exact, with
the bisection referee already positive on 13 of these `σ_mc` designs); the tallies below are
those of the canonical `k6.csv`.

*Primary variant.* The two designs the primary repair makes 0⁺-feasible are `voronoi_12`
(`F = 521`, 916 split edges) and `voronoi_93` (`F = 635`, 1124 split edges). Both are `POS`,
both have every `q_e > 0` and every `μ > 0`, and — since the F34 predicate fix
(`results/core_validation/referee_fix.md`; the old overlap predicate was rejecting the first
interval of the exact scan at a hinge vertex) — both have a positive exact range:
`Θ_max = 0.002354` (`voronoi_12`) and `Θ_max = 0.239603` (`voronoi_93`), agreed to `2.6e−8` by
the refereeing bisection, so the primary exact-scan row reads **2 / 400**. Neither is certified
at `ε = 0.3` (`NOROOT` fails), and 0.5 % is far below K6's 20 % bar. They also sit **10.6 and
12.8 median edge lengths** away from `X0` — the repair, where it succeeds at all, succeeds by
leaving the neighbourhood of the projection entirely.

*Secondary variant.* The split-only objective (no corner term, no proximity term) reaches a
0⁺-feasible point on 59 `σ_mc` designs, and on **15** of them — all Delaunay, `F` from 101 to
268 — the certificate at that point is positive with a large range: `ε_max` median 2.466 rad,
max π, and the bisection referee returns the same `Θ_max` on 14 of the 15 (on `delaunay_127`
the certificate gives 0.058 rad against a bisection of 3.023). These points are not near
`X0`: with `w_prox = 0` the minimiser is free to leave the projection's neighbourhood, and the
same 15 designs are the ones on which the split-only relaxation happens to land in the
deployable region. 7.5 % is the best family rate of the experiment and is still well under
the 20 % bar.

**The capped native comparison** (`results/kill/k6/native/{k6.csv,summary.txt}`, 40 graphs,
80 designs): 8 designs qualified and were handed to the authors' `prevent` from the same
`X0`, refereed by our bisection. Ours better / theirs better / tie = **1 / 0 / 3** under
`σ_mc` and **0 / 0 / 4** under `σ_def`; median `Θ_max` **0.0000 for both sides**;
`native_theta > 0` on **0 / 80**. Their code does not escape the jam either.

### The `first_root` bug, and the history of the secondary 15

An early pass of this experiment reported the **exact** `Θ_max` (T4.2″) as positive on
**15 / 200** `σ_mc` designs under the secondary variant — all `delaunay` — while the
bisection referee of that time agreed on only 13 and disagreed by up to 2.466 rad on
`delaunay_7`, `delaunay_82` and `delaunay_127`. Two defects were entangled in that pass.
`validity_certificate` set `cert.first_root` to the **first** admissible deflated root the
`(face-pair, vertex, edge)` loop happened to meet, not the **smallest** one over all
candidates; callers use it as the supremum of the `ε` for which `NOROOT(ε)` holds, so any
candidate scanned after the true first contact could raise the reported range. That was fixed
by the three-line guard now at `Kirigami/src/method/contact.jl:526` — keep a root only if it
is smaller than the incumbent — and with it alone the secondary exact scan fell to 0 / 400.
The second defect was the referee's: the pre-F34 overlap predicate rejected the first
interval at hinge vertices (`results/core_validation/referee_fix.md`). With both fixed, the
exact scan and the bisection agree on the same 15 designs (14 exactly; `delaunay_127` is the
one certificate-below-referee case), so the 15 are a measurement, not an artefact, and the
table above reports them. The 2.466 worst-case gap printed in
`results/kill/k6/native/summary.txt` (`referee: … 11 / 12, worst 2.466e+00`) is from the
40-graph native run, which used the early referee and is affected in that column only.

### Verdict: FAIL — 7.5 % at best against a 20 % bar, on both `σ` and all four repair variants

**What it means for F30.** The emptiness of the random-graph design space is **not** an
artefact of the `θ = 0⁺` jam being repairable and nobody having tried. The jam is
first-order, exactly characterised, and the repair problem is solvable on a third of the
designs (every `q_e > 0` and every face area positive on 126 of 400 under the split-only
objective, which does not constrain `μ`) — and solving it buys deployment range on only 17
of the 400 taking the best of primary / secondary / ladder / stage 2 per design (2 primary,
15 secondary; 15 of the 17 with `ε_max ≥ 0.1` rad, median 2.41 rad; all by leaving the
projection's neighbourhood), because on the rest a fresh contact is waiting
behind the one that was repaired. Fixing `q` moves the binding event from `split-inward` to
`vertex-edge` and then to `inverted`; the obstruction is not one bad sign but the whole
first-order cone being empty. F30 stands as a statement about the projection's
neighbourhood, now against four strategies rather than two: Eq. (6) alone, `σ_def`, 0⁺
repair, and the authors' native `prevent`.

**Consequence.** This is the result that ended the random-graph constructive programme
(orchestrator decision D8): the characterisation half of the work survives and is
strengthened, and the constructive half moves to periodic patterns, where the shape space is
not generically empty (K7).

Artifacts: `results/kill/k6/{k6.csv, summary.txt, cache/, native/}`. Driver
`Kirigami/apps/kill_k6.jl`, method `Kirigami/src/method/zero_plus.jl`, diagnostic
`Kirigami/apps/dbg_k6.jl`.

---

## Reproducing

```
julia --project=Kirigami -e 'using Pkg; Pkg.test()'   # 187 test sets, 176,692 assertions, 23 broken, 0 failures
julia --project=Kirigami Kirigami/apps/kill_k3a.jl  --out results/kill/k3a          # ~25 min
julia --project=Kirigami Kirigami/apps/kill_k3a_recheck.jl --out results/kill/k3a   # ~4 min
julia --project=Kirigami Kirigami/apps/kill_k1b.jl  --out results/kill/k1b          # ~34 s
julia --project=Kirigami Kirigami/apps/kill_k1a.jl  --out results/kill/k1a          # ~40 s
julia --project=Kirigami Kirigami/apps/kill_k5.jl   --out results/kill/k5 --n 200 --maxf 800   # ~5 min
julia --project=Kirigami Kirigami/apps/kill_k2c.jl  --out results/kill/k2c          # ~3 min
julia --project=Kirigami Kirigami/apps/kill_k2a.jl  --out results/kill/k2a          # ~7 s
julia --project=Kirigami Kirigami/apps/kill_k1c.jl  --out results/kill/k1c          # ~31 s   (needs baseline/native)
julia --project=Kirigami Kirigami/apps/kill_k2b.jl  --out results/kill/k2b --maxf 160   # ~26 min (needs baseline/native)
julia --project=Kirigami Kirigami/apps/kill_f23.jl  --out results/kill/f23          # ~1 s
julia --project=Kirigami/scripts Kirigami/scripts/plot_kill.jl
```

Figures: `k2a/k2a_agreement.png`, `k1c/k1c_rates.png`, `k3a/k3a_mobility.png`,
`k2c/k2c_locality.png`, `k2c/k2c_hloc.png`, `k2b/k2b_margin.png`, `k5/k5_orientation.png`.

Method code changed for the corrections: `Kirigami/src/method/deploy_basis.jl`
(`classify_harmonic`, `harmonic_roots_deflated`, `harmonic_no_root_in`) and
`Kirigami/src/method/contact.jl` (`exact_theta_max_overlap`, `validity_certificate`, the
exact face-frame radius in the broad phase). New unit test coverage in
`Kirigami/test/test_method.jl`: the hexagon graze; `Θ_max` vs bisection and pruning
invariance; completeness of `C(X)`; the swept radius in both frames; the three harmonic
classes and Lemma T5.1e; deflated roots vs a `τ`-chart-free sign-change scan on 60 000
harmonics; and the certificate implying no overlap on `(0, ε)`.

## A3 — jitter transition

Run after round-2 ideation, from `ideas/round2_adversary.md` Idea 4. It exists to answer the
one objection that the adversary persona puts first (rejection **R-1**): the emptiness result
is "0/200 random graphs, 8/8 authored tilings", and a referee cannot separate *"the design
space is empty"* from *"our generator and our orientation heuristic fail on our random
graphs"*, because random Voronoi/Delaunay patches differ from authored tilings in every
respect at once. The separating experiment is an interpolation with the **combinatorics held
fixed**.

**Hypothesis.** Deployability dies as a **sharp function of one geometric ratio**, not of
"randomness". Jittering an authored tiling's interior vertex positions continuously should
therefore reproduce the random-graph regime, with a per-tiling transition amplitude `a*`, and
some scalar computable from `(G, σ, X_ini)` alone should predict the outcome.

**PASS bar** (restated verbatim in intent from the idea file's kill rule): a sharp `a*` exists
per tiling **and** some predictor reaches **ROC AUC ≥ 0.9** for "certified `Θ_max > 0`".
**FAIL** if the certified fraction is already ≈ 0 at the smallest amplitude (the tilings would
be knife-edge, so emptiness is generic and not a generator artefact), if there is no monotone
collapse, or if every predictor has AUC < 0.9.

### How it was run

* **Driver** `Kirigami/apps/kill_jitter.jl`. No `core/` or `method/` source was modified.
* **Population.** The eight `reference_cases()` authored tilings — the same patches K2a/K5
  use — with `X[v] ← X[v] + a·(median edge)·N(0, I)` on **interior vertices only**, for `a` on
  **20 log-spaced steps in [0.005, 1.0]**, **20 seeds** each, plus one `a = 0` baseline per
  tiling. **3208 designs**, each measured under **two σ rules** = 6416 rows. Boundary vertices
  are not jittered because the pipeline solves Eq. (4)/(6) with `BoundaryMode::Fixed`.
* **σ rules, as in K5.** `σ_mc` = the reference case's own orientation (checkerboard where the
  tiling is 2-colourable, else `assign_orientation_relaxation`, Eq. (1)); it depends on the
  dual graph only, so it is constant along the ladder — the combinatorics are never
  randomised. `σ_def` = K5's greedy defect search (single-face and adjacent-pair flips,
  strictly improving, `c(Γ) = 1` enforced) re-run on each jittered geometry, from `σ_mc`, cap
  `3|F|` attempts (K5 used 4 starts and `20|F|`; the ladder is 3208 designs).
* **Measurement, unchanged pipeline.** `make_cut → holes_partition → assemble_system(Fixed) →
  solve_system → X0`, then `validity_certificate` (`POS ∧ NOOVERLAP(ε/2) ∧ NOROOT`) at
  `ε = 0.006` (the idea file's value) and `ε = 0.3` (K5's), and the **exact** `Θ_max` by
  `exact_theta_max_overlap` — the T4.2″ interval scan K2a validated.
* **Sharding.** 12 shards by design index mod 12, `xargs -P 12`. **Wall 2.2 s.** Analysis
  (`--analyze`) is also Julia; python draws the figure only.
* **Predictors scored.** From `(G, σ, X_ini)` before any null space: `rho_geom` (the idea's
  own ratio — median over split-forest components of geometric diameter / median face
  inradius), `rho_hop`, `rho_max`, `psplit = |E_split|/|F|`, `badq_ini` (fraction of split
  edges with `q_e ≤ 0`), `aniso`, `D(σ)`, and `amp` as a positive control. Two `X0`-derived
  controls are reported separately: `badq_X0`, `margin_X0`.

### Measured — `results/kill/jitter/{jitter.csv, ladder.csv, summary.txt, transition.png}`

**The population splits in two, and the split is `|E_split|`.** Four of the eight reference
cases are **split-free** under their own σ: `rotating_squares`, `triangles_alternating`,
`kagome_3636`, `periodic_squares_4x4`. For these `dim_null = 0`, so the fixed-boundary Eq. (4)
system has a **unique** solution — and the authored positions already satisfy it. The Eq. (6)
projection therefore sends **any** jittered interior geometry back to the authored tiling
**exactly**: measured `‖X0 − X_base‖_∞ = 0.0000` at every amplitude up to `a = 1.0`, while
`‖X0 − X_ini‖_∞` reaches 2.8 median edges. Those four curves are flat at 1.0 not because they
tolerate jitter but because **they never left the start point**. New unit test in
`Kirigami/test/test_method.jl` ("split-free tiling: Eq. (6) projection undoes an arbitrary
interior jitter"): `dim_null = 0`, `‖X0 − X_base‖ < 1e−9` at `a ∈ {0.05, 0.5, 1.0}`. This is a
result in its own right and it means **the ladder is only informative on the four
split-bearing tilings** (`hexagons_auto` 5, `truncated_square_488` 12, `snub_square_33434` 13,
`tiling_3_4_3_12` 4 split cuts under `σ_mc`), where `‖X0 − X_base‖ ≈ ‖X0 − X_ini‖`.

**On those four there is a monotone, reasonably sharp transition in the exact `Θ_max`.**

| tiling (σ_mc) | `a*` (Θ_max > 0) | `a` at frac ≥ 0.9 | `a` at frac ≤ 0.1 | width |
|---|--:|--:|--:|--:|
| `snub_square_33434` | 0.163 | 0.107 | 0.248 | 0.36 dec |
| `hexagons_auto` | 0.277 | 0.142 | 0.573 | 0.61 dec |
| `truncated_square_488` | 0.320 | 0.248 | 0.573 | 0.36 dec |
| `tiling_3_4_3_12` | 0.522 | 0.248 | 0.757 | 0.48 dec |

`a*` under `σ_def` is uniformly **smaller** (0.071 / 0.188 / 0.096 / already dead at 0.005):
consistent with K5 and K6, the defect-minimising orientation buys nothing here and roughly
doubles-to-quintuples `|E_split|` (`snub` 13 → 25, `3_4_3_12` 4 → 20). Over all split-bearing
rows `Θ_max > 0` holds on **1225/1604 (76.4 %)** under `σ_mc` and **634/1604 (39.5 %)** under
`σ_def`.

**As first run, the certificate did not track this, and that finding decided the verdict
and produced F32.** In that pass, at the smallest amplitude `a = 0.005`, every split-bearing
tiling had `POS 20/20` and `NOOVERLAP 20/20` and exact `Θ_max` between 1.0 and 2.4 rad on
**20/20 seeds** — yet the certificate passed on 0/20 (`hexagons_auto`,
`truncated_square_488`), 5/20 (`snub_square`), 10/20 (`tiling_3_4_3_12`); `NOROOT` was the
only clause that failed, and the unjittered `a = 0` baselines certified 6/8 under `σ_mc`.
Over the split-bearing ladder, **1622 of the 1859 rows with `Θ_max > 0` were not certified —
an 87.3 % false-negative rate**. The diagnosis (`results/kill/jitter/cert_diagnosis.md`): an
infinitesimal perturbation of an exactly symmetric tiling splits a degenerate contact and
pushes a **collinearity** root — the duplicate edge on the *line* of its partner, not on the
segment — into `(0, 0.006)` while the first real overlap stays at ≈ 2 rad, and `NOROOT` was
counting it because it lacked the two interval tests. With the interval tests restored
(F32, `contact.jl`) the canonical `jitter.csv` reads: certificate holds on **20/20** seeds of
all four split-bearing tilings at `a = 0.005`, the `a = 0` baselines certify **8/8**, and over
the split-bearing ladder **1856 of the 1859 rows with `Θ_max > 0` are certified** — the three
exceptions (`hexagons_auto` at `a = 0.33` and `0.43`, `snub_square` at `a = 0.027` under
`σ_def`) have `Θ_max` between 0.0018 and 0.0060 rad, i.e. below `ε = 0.006`, and are correct
rejections. `a*(certified)` now coincides with `a*(Θ_max > 0)` on every tiling (0.272 against
0.277 on `hexagons_auto`, identical elsewhere). F30's note that the certificate is a strict
inner approximation stands in principle (exact below the first graze, inner above), but on
this population the inner gap is empty.

**Predictors, ROC AUC (directed = max(AUC, 1−AUC)).** The right population is the
split-bearing rows; the all-rows column is contaminated by the split-free/split-bearing
dichotomy, on which every `E_split`-derived scalar is trivially perfect.

| predictor | AUC(Θ_max > 0), all rows | AUC(Θ_max > 0), split-bearing | AUC(certified), split-bearing |
|---|--:|--:|--:|
| `rho_geom` (the idea's own ratio) | 0.918 | **0.778** | 0.777 |
| `rho_hop` | 0.892 | 0.706 | 0.706 |
| `rho_max` | 0.952 | 0.870 | 0.870 |
| `psplit` = split cuts per face | 0.874 | 0.657 | 0.656 |
| `badq_ini` | 0.987 | **0.976** | 0.976 |
| `aniso` | 0.745 | 0.816 | 0.816 |
| `D(σ)` | 0.585 | 0.504 | 0.504 |
| `amp` (positive control) | 0.743 | 0.831 | 0.831 |
| `margin_X0` (X0 control) | 1.000 | 1.000 | 1.000 |
| `badq_X0` (X0 control) | 1.000 | 1.000 | 0.999 |

(The certified column equals the `Θ_max > 0` column to three decimals because, after F32,
the certificate labels differ from the exact-range labels on 3 of 3208 split-bearing rows.)

The AUC routine was cross-checked against an independent brute-force Mann-Whitney count on
`badq_ini` (0.012789 vs 0.0128 — identical).

### Verdict

**MIXED, and FAIL as the rule is literally written.** Taken clause by clause:

* The FAIL clause *"certified fraction already ≈ 0 at the smallest amplitude"* **fired on the
  first pass** — but its stated meaning ("the tilings are knife-edge, so emptiness is
  generic") is **not supported**. Exact `Θ_max` is 1.0–2.4 rad on 20/20 seeds at `a = 0.005`
  for all four split-bearing tilings, and after F32 the certificate agrees on 20/20. What was
  knife-edge was the **certificate's `NOROOT` clause**, not the geometry.
* On the exact `Θ_max` reading the transition **exists and is monotone and sharp** (0.36–0.61
  decades), with `a* ∈ [0.16, 0.52]` under `σ_mc`. Deployability survives perturbations of
  order **one sixth of an edge length** and dies around **one third to one half**.
* The predictor bar is **met only by a restatement of the mechanism**. `badq_ini` reaches
  0.976 on split-bearing rows, but "some split cut already opens the wrong way at `X_ini`"
  *is* the K5/K6 failure mode, not an independent geometric characteristic. **The adversary's
  own proposed ratio `rho` reaches only 0.778** and fails the 0.9 bar. `rho_max` (worst
  component instead of median) does better at 0.870 and still fails.

**What this does and does not license.** It **does** answer R-1 in part: emptiness is a
function of measured geometry, reached continuously from patterns whose combinatorics were
never randomised, with a transition amplitude that is neither 0 nor 1. It **does not** deliver
the promised generator-independent predictor, and it **does not** support the sentence the
idea file wanted ("eight curves collapsing onto one when the x-axis is `rho`"): four of the
eight tilings carry no signal at all, and on the other four `rho` is a weak separator. Two
consequences that were not anticipated and should be carried forward:

1. **Split-free authored patterns are outside the whole problem.** Their Tutte auxetic shape
   space is a point, and it is the authored tiling. Any claim of the form "authored tilings
   pass 8/8" is really "4 trivially, 4 substantively".
2. **The certificate needed the interval tests** before it could be used as a label on
   authored geometry: as first run, `NOROOT(ε)` at `ε = 0.006` rejected 87 % of designs whose
   true range exceeds 1 rad. That was F32; with the fix the rejection rate on this ladder is
   3 of 1859, all correct.

Nothing here contradicts `STATE.md`. F30's "the certificate is a strict inner approximation"
holds in principle and is not exercised on this population after F32; F25's "authored tilings
are fine, 0/8 fail" survives at `a = 0` under both the `Θ_max > 0` reading and the certificate
(8/8), and the "8" should carry the split-free caveat above.

### Reproducing

```
julia --project=Kirigami -e 'using Pkg; Pkg.test()'   # 187 test sets, 176,692 assertions, 23 broken, 0 failures
mkdir -p results/kill/jitter
seq 0 11 | xargs -P 12 -I{} julia --project=Kirigami Kirigami/apps/kill_jitter.jl --shard {} --nshards 12 \
    --out results/kill/jitter                        # ~2 s wall
cd results/kill/jitter && head -1 jitter_0.csv > jitter.csv \
  && for f in jitter_[0-9]*.csv; do tail -n +2 $f; done | sort -t, -k1,1 -k2,2n -k4,4n -k8,8 >> jitter.csv
julia --project=Kirigami Kirigami/apps/kill_jitter.jl --analyze results/kill/jitter/jitter.csv --out results/kill/jitter
julia --project=Kirigami/scripts Kirigami/scripts/plot_jitter.jl results/kill/jitter
```

Figure: `results/kill/jitter/transition.png` (deployable fraction, certified fraction, and
median `Θ_max` against amplitude; grey = the four split-free tilings whose projection undoes
the jitter). Per-rung table: `results/kill/jitter/ladder.csv`. Every row is a deterministic
function of `(tiling, amplitude index, seed)`, so shards may be run in any order.

## K7 — achievable periodic Jacobians (C1–C4)

Run under D8, after K6's FAIL, to test the constructive half of the deliverable on the
population where the shape space is *not* generically empty: **periodic** patterns on a torus
rather than random graphs with a fixed boundary. Spec `specs/experimenter_k7.md`; the four
claims are R4 (C1–C3) and B20 (C4) from `ideas/ranking.md`.

**Hypothesis.** On a torus the deployment of a periodic pattern is an **affine** family. With
`P_θ = [p_x(θ) p_y(θ)]` the deployed period matrix and `J(θ) = P_θ P_0⁻¹`,

* **C1** `J(θ) = cos(θ/2) I + sin(θ/2) K(X)` **exactly**, with `K` constant in `θ` and
  **affine** in the shape-space coordinate `X = X_0 + Φ t`. The achievable set
  `𝒦 = {K(X) : X ∈ 𝒳}` is then an affine subspace of the 2×2 matrices, and inverse design of
  the whole `θ`-history reduces to **one linear solve**.
* **C2** the conformal designs are the solutions of two linear equations in `t`
  (`K₁₁ = K₂₂`, `K₁₂ = −K₂₁`), and such a design is conformal at **every** `θ`, not just
  infinitesimally — this is the sentence 2026 Sec. 5.1 asserts "empirically for all θ".
* **C3** the Poisson-ratio curves `ν(θ, d)` are parameterised by `K ∈ 𝒦`, so a **target** `K*`
  can be hit by least squares and the leftover freedom spent on maximising certified `Θ_max`.
* **C4** the per-cell hole area is a first harmonic, so the **second fully-closed angle** has a
  closed form.

**PASS bars** (verbatim from the spec): **C1** affine and exact, all tolerances (`1e-9`) met on
≥ 30 patterns. **C3** target hit `< 1e-8` with a certified, `Θ_max > 0` design on **≥ 80 %** of
the patterns with `dim 𝒦 ≥ 1`. **C4** `|θ_c − measured| < 1e-6` on every pattern that has a
second closed state. C2 has no numeric bar in the spec beyond "verify `< 1e-9` at 200 angles".

### How it was run

* **Driver** `Kirigami/apps/kill_k7.jl`, stages `main`, `c3`, `c4bounded`, `nu`, `c4sweep`,
  `perdbg`. The periodic layer it uses — lattice detection, the torus quotient, the super
  patch, `periodic_jacobian`, the achievable set — is `Kirigami/src/method/periodic_jacobian.jl`,
  new in this experiment. `core/` was not modified.
* **Why a quotient and not a finite patch.** `Kirigami/README.md` deviation 11 imposes Eqs. (3b)–(3c)
  on a patch whose boundary edges are *not* topologically identified, so every hole preimage
  straddling the seam contributes no row of `L` (deviation 2 drops preimages touching the mesh
  boundary) and exactly the wrap-around deployability conditions go missing. K7 therefore
  builds the genuine quotient: one face per translation class, edges keyed by
  (class pair, lattice offset), hole rows carrying the integer offset sum on the right. The
  finite-patch `dim_null` is reported alongside as `dim_null_patch` and is systematically
  larger (e.g. `voronoi_torus_3_n45`: 45 on the quotient, 63 on the patch).
* **Population — 33 patterns, the full spec population, no row missing.** Seven periodic
  families (`squares`, `triangles`, `hexagons`, `kagome`, `snub_square`, `trunc_square_488`,
  `t3_4_3_12`) at 2×2, 3×2, 3×3 cells = 21, plus **12 periodic Voronoi tori** of 20–200 faces
  per cell (`voronoi_torus_0_n20` … `_11_n200`), σ from Eq. (1) solved as max-cut **on the
  quotient dual** rather than on a boundary-truncated patch dual.
* **Measurement.** `K` comes from the closed-form deploy basis; every check against it uses an
  **independent** `deploy()` forward-kinematics call, so C1 is basis-vs-FK, not basis-vs-basis.
  Certification is the same `POS ∧ NOOVERLAP(ε/2) ∧ NOROOT` certificate K2a/K6 use, with the
  `first_root` fix.
* **Sharding.** 12 shards by pattern index mod 12. Wall for the `main` stage is set by the
  largest tori (`voronoi_torus_11_n200`, 1824 s alone); the `c3` stage is seconds per pattern.
  The two previous K7 agents were killed by session limits with the `c3` stage 40/50 rows in;
  those ten rows were finished here by addressing single patterns (`--nshard 33 --shard N`),
  **no completed row was recomputed.**

### Measured — `results/kill/k7/{k7_main.csv (33), k7_c3_all.csv (50), k7_c4_bounded.csv (7), k7_periodicity.csv (33), k7_nu_*.csv, k7_c4_sweep.csv}`

**C1 — PASS, and the dimension formula is not the one the spec guessed.**

| quantity | bar | worst over 33 patterns |
|---|--:|--:|
| `‖P_θ^fit − P_θ^FK‖ / ‖P_θ^FK‖`, 50 angles | `1e-9` | **7.29e−14** |
| `‖K_fit − K_basis‖_∞` | `1e-9` | **5.99e−14** |
| second difference of `t ↦ K` (affinity), 31 triples/pattern | `1e-9` | **4.79e−14** |
| `‖P_0 − T‖_∞` (measured period vs the lattice) | — | 6.39e−14 |

`K` is constant in `θ` and affine in `t` to machine precision on **33/33**. The spec asked us
to compare `dim 𝒦` with `min(4, 2·dim_null)`:

* **`dim 𝒦 = 2·rank(D)` on 33/33**, where `D` is the `k × 2` matrix of period-potential
  increments `w_t` — this is the formula, and it is exact.
* **`min(4, 2·dim_null)` fails on 1/33.** `squares_3x3` has `dim_null = 4` (6 split cuts) yet
  `dim 𝒦 = 0`: its `K` is frozen at a **rotation** (`tr K = 0`, `det K = 1`) over the entire
  4-dimensional shape space, so `J(θ)` is a rigid rotation by `θ/2` at every design and the
  cell never opens. Shape-space dimension therefore does **not** bound designability; the rank
  of the period-potential map does.

Where `dim 𝒦 = 4` (25 of 33 patterns) the achievable set is **all** of the 2×2 matrices — the
map `t ↦ K` is surjective. Figure `k7_dimK_hist.png`.

**C2 — PASS, and it settles 2026 Sec. 5.1.** On all **25** patterns with `dim 𝒦 ≥ 1` the two
linear conformality equations are consistent (residual ≤ **9.19e−16**), and at the solution the
conformal distortion of `J(θ)` over **200 angles spanning `[0, π]`** is ≤ **2.24e−14** (best:
1.92e−16). This is not an empirical coincidence: since `J(θ)` is an
affine path through `I` and the similarities are the linear span of `{I, Jrot}`, `J(θ)` is a
similarity for **every** `θ` exactly when `K` is one. The paper's "once the derivative of the
conformal distortion vanishes at `θ = 0`, the deployment remains conformal for all opening
angles" is therefore a one-line consequence of C1, now with a numerical witness.

**C3 — FAIL, 52 % against an 80 % bar. The achievable set is not the obstruction.**

| target | patterns | max hit error | `0⁺` feasible | **certified `Θ_max > 0`** |
|---|--:|--:|--:|--:|
| random `K* ∈ 𝒦` | 25 | **1.44e−15** | 21/25 | **13/25 (52 %)** |
| anisotropic `K* = diag(1, −0.5)` | 25 | 5.97e−15 (24) / **1.0** (1) | **0/25** | **0/25** |

Certification is under the **current** certificate, which acquired the T4.1b interval
(segment) test on admissible NOROOT roots while K7 was being written (F32,
`results/kill/jitter/cert_diagnosis.md`); the certified `ε` agrees with the exact `Θ_max` to
`2.2e−9` rad on all 13 certified designs. The random target `K*` is drawn inside `𝒦` from a
null-space basis, and which basis is returned — hence which target is drawn and which
design realises it — depends on the linear-algebra backend, so the certified count of this
row is not reproducible across floating-point environments (earlier passes of the same driver
gave 11/25 and 12/25); the anisotropic row, whose target is fixed, is. The `c3_cert` /
`c3_eps_cert` columns of `k7_main.csv` are from the `main` stage's own draw and should be read
from `k7_c3_all.csv` instead. No pass comes close to the 80 % bar, and none certifies a single
anisotropic design.

The linear solve does everything C1 promises: on 24 of 25 patterns even the *anisotropic*
target is hit to `8.9e−16`, because `dim 𝒦 = 4` makes `𝒦` the whole matrix space. The single
target-unreachable case is `squares_3x2`, the only pattern with `rank(D) = 1` (2 split cuts,
`dim 𝒦 = 2`); its hit error is exactly `1.0`, and that is a genuine *achievable-set* failure.
Every other failure is the **same `0⁺` obstruction as K5/K6/F30**: the design that realises
`K*` places at least one split cut so that its two duplicates translate **into** each other at
`θ = 0⁺`, i.e. `min_e q_e < 0`, and `Θ_max = 0` before any certificate is consulted.

Per family, and this is where the mechanism differs:

| family | `nsplit` per cell | `dim 𝒦` | `rank D` | certified (random target) | `min q` at the anisotropic design |
|---|--:|--:|--:|--:|--:|
| `squares` 2×2, `triangles`, `kagome` (8 patterns) | **0** | **0** | 0 | — (excluded, `𝒦` is a point) | — |
| `squares_3x3` | 6 | **0** | 0 | — (excluded, `K` is a rotation) | — |
| `squares_3x2` | 2 | 2 | 1 | 0/1 | −33.5, **and the target is unreachable** |
| `hexagons` | 4–9 | 4 | 2 | 2/3 | −5.91 … −0.30 |
| `snub_square` | 8–18 | 4 | 2 | 2/3 | −41.5 … −14.9 |
| `trunc_square_488` (4.8.8) | 8–18 | 4 | 2 | 2/3 | −2.66 … −0.35 |
| `t3_4_3_12` | 8–24 | 4 | 2 | 3/3 | −49.0 … −38.7 |
| `voronoi_torus` | 20–204 | 4 | 2 | **4**/12 | −15.0 … −0.29 |

Read the table as follows. **`kagome`, `triangles` and `squares_2x2` are outside the problem**:
under the max-cut σ on the quotient dual they are **split-free** (`nsplit = 0`), so `dim_null = 0`,
the shape space is a single point and `K` cannot be moved at all — exactly the split-free
dichotomy A3 found on authored patches, now on the torus. **`squares_3x3` is the opposite
degeneracy**: 6 split cuts and a 4-dimensional shape space that moves no vertex in a way the
period Jacobian can see. **`dim 𝒦 = 4` with `rank D = 2` on the remaining 25 means the
achievable set is everything** — so "which `K` can this pattern realise" is the wrong question
for them, and 2026's designable-Poisson-family framing does not need any pattern-dependent
reachability caveat. What it does need is a validity caveat: on `snub_square` and `t3_4_3_12`
the anisotropic design overshoots by `min q ≈ −39` to `−49`, two orders of magnitude past the
boundary, i.e. hitting `diag(1, −0.5)` demands vertex motions that fold split cuts deep into
their neighbours. On `hexagons`, 4.8.8 and Voronoi the overshoot is `O(1)` and sometimes only
`−0.29` — those are the patterns where a margin-aware search might still recover a design, and
`voronoi_torus_0_n20` (`min q = −0.29`) is the closest near-miss in the population.

Where a design **is** certified, the closed-form `ν(θ)` is right: over the 13 certified designs
the largest discrepancy between the predicted `ν` from `K` and forward kinematics at 20 angles
and 2 directions is **1.15e−09**, and both the certified `ε` and the exact `Θ_max` reach
**2.18 rad** (`t3_4_3_12_2x2`). Figures `k7_nu_curves.png` (predicted lines vs FK circles, four
patterns, plus the conformal design as the flat dashed curve) and `k7_c3_certified.png`.

**C4 — PASS. The closed form, derived rather than assumed.** `A(θ) = det(P_0)(det J(θ) − 1)`;
substituting `J(θ) = cos(θ/2) I + sin(θ/2) K` and the half-angle identities gives the first
harmonic `A(θ) = q(cos θ − 1) + r sin θ` with

```
q = -det(P_0) (det K - 1) / 2 ,    r = det(P_0) tr K / 2 ,
theta_c = 2 atan2(r, q) = 2 atan2(tr K, 1 - det K) .
```

The right-hand form is the **periodic version the spec asked for and could not assume**; it is
derived here and holds identically, since `det P_0` cancels. Note `p + q = 0` automatically, so
`θ = 0` is always a root — the second root is the content.

| population | rows | worst `|θ_c − θ_measured|` | bar |
|---|--:|--:|--:|
| periodic (bisection on FK hole area) | **32** (of 33) | **5.62e−14** | 1e−6 |
| bounded patches, `k7_c4_bounded.csv` | 7 | **4.44e−16** | 1e−6 |

The excluded periodic row is `squares_3x3`, which has `q = r = 0`: its `A(θ) ≡ 0`, so it has
**no** second closed state, consistent with `K` being a rotation. On the bounded patches the
closed state is overlap-free (`θ_c ≤ Θ_max`, the 2026 equal-length condition) on **4/7** —
`rotating_squares` (π), `kagome_3636` (π), `hexagons_auto` (2π/3), `truncated_square_488`
(3π/4) — and unreachable on `triangles_alternating` (`θ_c = 4π/3 > π`), `snub_square_33434`
(`θ_c = 3.768` vs `Θ_max = 1.646`) and `tiling_3_4_3_12` (`2.818` vs `2.381`). Predicting the
angle is exact; reaching it is a separate, and often negative, question.

**One anomaly, chased down and resolved.** The area-harmonic fit error over a full
`θ ∈ (0, 2π)` sweep is `1e−15` on the 21 tilings but reaches **2.3e−2** on the Voronoi tori.
It is **not** a failure of the harmonic. `k7_c4_sweep.csv` / `k7_c4_area_fit.png` localise the
deviation to exactly the interval where **`det P_θ < 0`** — the deployed cell has inverted —
and it is an artefact of the measurement using `|det P_θ|` where the harmonic is the signed
`det(P_0)(det J − 1)`. Off that interval the agreement is `1e−15` at every sampled angle. On
`voronoi_torus_0_n20` the inverted interval is `θ ∈ (3.97, 4.44)`, far past both `θ_c = 2.100`
and any deployable range. The `squares_3x2` fit error of 0.207 is the same effect.

**A second anomaly, which is a real caveat.** `periodic_jacobian`'s consistency statistic —
the spread of the per-corner period about its mean — is `≤ 1.1e−13` on **32/33** patterns and
**27.5** on `snub_square_3x3`. `k7_periodicity.csv` separates the two possibilities: the
quotient residual there is `2.4e−15` and the *mean* period still follows
`cos(θ/2) P_0 + sin(θ/2) Q` to `3.0e−14`, but the cut structure `M'` of its super patch has
**9 connected components**, each threading all 9 cell copies. `deploy()` gauges each component
independently, so the period read off a corner in one component differs from another by a
gauge offset that grows with `θ` (0.69 at `θ = 0.05`, 14.4 at `θ = 1.1`, 27.4 at `θ = 3.0`).
**`snub_square_3x3` is the one pattern whose `K` is gauge-dependent and should not be quoted**;
it contributes to the C1/C2/C4 tables above, but removing it changes none of the worst-case
figures there, because it is never the maximum in any of them. Seven other patterns have `2 ≤ ncomp ≤ 12` with no such effect.

### Verdict

| claim | bar | measured | verdict |
|---|---|---|---|
| **C1** affine Jacobian | all tolerances `1e-9`, ≥ 30 patterns | 7.29e−14 / 4.79e−14 on **33** | **PASS** |
| **C2** exact conformality | `< 1e-9` at 200 angles | **2.24e−14**, 25/25 patterns | **PASS** |
| **C3** designable Poisson family | ≥ 80 % certified, `Θ_max > 0` | **13/25 = 52 %** (basis-dependent draw; 11/25 and 12/25 in earlier passes) | **FAIL** |
| **C4** closed-form `θ_c` | `< 1e-6` on every pattern with one | **5.62e−14** on 32 + 7 | **PASS** |

**The split is clean and it is the same split as everywhere else in this bundle.** The
*algebra* of R4 is right and better than claimed: `J(θ)` is exactly affine, the achievable set
is an affine subspace and on 25 of 33 patterns it is the **whole** 2×2 matrix space, conformal
design is two linear equations, and the second closed angle is a two-argument arctangent of
`tr K` and `det K`. Every one of those is a closed form we did not have and that the 2026 paper
either asserts without proof (C2) or does not state (C4's periodic form, and `dim 𝒦 = 2 rank D`).
The *geometry* is where it dies, in the identical way as K5, K6 and F30: hitting a prescribed
`K` is a linear solve with a large solution set, and the solution set is dominated by
embeddings whose split cuts open inward at `θ = 0⁺`. **B20 (C4) passes outright. R4's
characterisation half (C1, C2) passes; R4's algorithm half (C3) fails.**

**Contradiction with `STATE.md`: one, and it is a correction to an in-flight note.** The S2 04:14
session-limit entry records K7 as "31/60 rows done". That undercounts and mis-scopes: the `main`
stage was in fact **complete at 33/33** when the agent was killed, and the missing work was ten
`c3` rows out of 50, now finished. There is no 60-row population. Otherwise K7 is consistent
with the rest: F15's `dim_null = |E_split|` holds on the quotient for every split-bearing
pattern here (e.g. `hexagons_3x3` 9/9, `voronoi_torus_11_n200` 204/204); F22's periodic rank
result is what makes the quotient system well-posed; F30's `0⁺` mechanism is confirmed on a
population that shares **none** of K5/K6's construction — periodic instead of fixed-boundary,
tori instead of random Delaunay patches, tilings as well as random Voronoi. That K7 reproduces
it independently is the strongest evidence in the bundle that the `0⁺` obstruction is a
property of the problem and not of one generator.

**New unit test.** `Kirigami/test/test_method.jl`, "C4: theta_c = 2 atan2(tr K, 1 − det K) is where
forward kinematics recloses" — builds `hexagons_2x2` and `t3_4_3_12_2x2` through the same
`make_tiling_pattern` the experiment uses, checks the two forms of the closed form against each
other and against the recorded values (2.09439510239 and 2.72771226181), checks the harmonic
against forward kinematics on `(0, θ_c]`, checks the hole area is positive at `θ_c/2` and zero
at `θ_c`, and brackets the sign change by 80 bisection steps to the `1e-6` K7 bar. It passes: **1 case, 24
assertions**. Whole suite after the concurrent F32 certificate change and its own test fixes:
**66 cases, 16283 assertions, all passing.**

### Reproducing

```
julia --project=Kirigami -e 'using Pkg; Pkg.test()'   # 187 test sets, 176,692 assertions, 23 broken, 0 failures
julia --project=Kirigami -e 'using Pkg; Pkg.test(test_args=["test_method_3"])'   # the file holding the K7 test (test_args filters by file name); its "C4: theta_c" test set alone: 24 assertions, 0 failures
mkdir -p results/kill/k7

# main stage: 33 rows, 12 shards. Wall is set by voronoi_torus_11_n200 (~30 min alone).
for i in $(seq 0 11); do julia --project=Kirigami Kirigami/apps/kill_k7.jl --stage main --out results/kill/k7 \
    --shard $i --nshard 12 > results/kill/k7/shard_$i.log 2>&1 & done; wait

# C3 (50 rows), the bounded-patch C4 table, the nu curves, and the two diagnostics
for i in $(seq 0 11); do julia --project=Kirigami Kirigami/apps/kill_k7.jl --stage c3 --out results/kill/k7 \
    --shard $i --nshard 12 > results/kill/k7/c3_$i.log 2>&1 & done; wait
julia --project=Kirigami Kirigami/apps/kill_k7.jl --stage c4bounded --out results/kill/k7
julia --project=Kirigami Kirigami/apps/kill_k7.jl --stage nu        --out results/kill/k7
julia --project=Kirigami Kirigami/apps/kill_k7.jl --stage c4sweep   --out results/kill/k7   # the |det P| artefact
for i in $(seq 0 11); do julia --project=Kirigami Kirigami/apps/kill_k7.jl --stage perdbg --out results/kill/k7 \
    --shard $i --nshard 12 > results/kill/k7/perdbg_$i.log 2>&1 & done; wait

# merge (k7_main.csv, k7_c3_all.csv, k7_periodicity.csv), then figures
cd results/kill/k7 && head -1 k7_main_0.csv > k7_main.csv \
  && for f in k7_main_[0-9]*.csv; do tail -n +2 $f; done >> k7_main.csv
head -1 k7_c3_0.csv > k7_c3_all.csv \
  && for f in k7_c3_[0-9]*.csv; do tail -n +2 $f; done | sort -u >> k7_c3_all.csv
head -1 k7_periodicity_0.csv > k7_periodicity.csv \
  && for f in k7_periodicity_[0-9]*.csv; do tail -n +2 $f; done >> k7_periodicity.csv
cd - && julia --project=Kirigami/scripts Kirigami/scripts/plot_k7.jl results/kill/k7
```

Figures: `k7_nu_curves.png` (C3, predicted vs measured `ν(θ)` on four patterns),
`k7_dimK_hist.png` (C1, `dim 𝒦` against both candidate formulas),
`k7_c3_certified.png` (C3, per family: patterns, `0⁺`-feasible, certified),
`k7_c4_area_fit.png` (C4, the `det P_θ < 0` artefact). Every pattern is a deterministic
function of its index in `population()` — the tiling seeds are `20260904 + 7919·i`, the Voronoi
seeds `9100001 + 104729·i` — so shards may be run in any order, and a single pattern `i` can be
re-run alone with `--nshard 33 --shard i`.

---

## B4 — does the fixed boundary empty the space?

Run after round-2 ideation, from `ideas/round2_theorist_b.md` idea **B4**. It exists to
separate two readings of F30/K6 that the K6 run itself cannot tell apart. K6 measured
`2/400` certified under its primary repair on 200 random Voronoi/Delaunay/quad patches; K7
ran the *same* Voronoi generator on tori and got `11/12` zero-plus feasible and `4/12`
certified (`results/kill/k7/k7_c3_all.csv`, `random_in_K` rows). The two populations differ in
exactly one thing: the patch carries 2026 Eq. (4)'s **fixed boundary rows** `B`, the torus
does not. B4's claim is that `B`, not graph randomness, is what empties the space — by the
budget identity (B.4) of that file, the total first-order hole-opening rate is a functional
of the **border** alone, so freezing the border freezes the budget.

**The three boundary variants**, all on K6's population verbatim (the same 200 graphs from
`kill_common.jl::make_graph` (frozen in `data/corpus/`), the same `σ_mc` from Eq. (1) and `σ_def` read back from
`results/kill/k5/sigma/`, the same `X_ini`, the same repair objective, weights
`w_corner = 0.05`, `w_prox = 1e3`, `λ = 1e−6`, `max_iter = 1200`, 3 Gaussian starts and the
same seeds `6000 + 7·id + which`):

| variant | system solved | why |
|---|---|---|
| `fixed` | `[L; B] X = [0; T]` | the K6 system, run as a **control**; the shape space is read straight out of K6's own cache, so this is literally the design K6 measured |
| `free` | `[L; e_0ᵀ] X = [0; x_0]` | `B` dropped entirely. The only gauge `L` cannot see is the global translation (every `L` row sums to zero, T7 (iv)), killed by pinning one vertex — exactly what `BoundaryMode::Periodic` already does for the same reason. Global rotation and scaling are left in: `q_e` is rotation-invariant and positively homogeneous under scaling, so they change no label |
| `split_free` | `B` restricted to boundary vertices touching **no** split edge | the F23 question one step further: free exactly the border where the cuts, and hence the budget, live |

Downstream everything is K6's pipeline unchanged: Eq. (6) projection of the same `X_ini`,
`zero_plus_repair` (secondary split-only first, then the primary split + vertex-edge +
proximity objective warm-started from it), the **exact `Θ_max` of T4.2″ as the primary
label** (the certificate is under repair, `STATE.md` F32) and the certificate as the
secondary label, with the independent overlap bisection as referee. 12 disjoint shards by
`gidx % 12`.

**PASS rule (B4):** exact `Θ_max > 0` on **≥ 10 %** of the 400 `free` designs, against
the `fixed` control.

### Measured — `results/kill/b4/b4.csv`, 1200 rows (200 graphs × 2 σ × 3 variants)

| variant \| σ | designs | median `dim_null` | median `min q` at `t=0` | 0⁺-feasible (split-only) | 0⁺-feasible (primary) | **exact `Θ_max > 0`** | certified `eps_max > 0` |
|---|--:|--:|--:|--:|--:|--:|--:|
| `fixed` \| `σ_mc` | 200 | 126.5 | −77.81 | 59 | 0 | **0** | 0 |
| `split_free` \| `σ_mc` | 200 | 131 | −66.06 | 85 | 0 | **0** | 0 |
| `free` \| `σ_mc` | 200 | 157 | −51.49 | **176** | 0 | **0** | 0 |
| `fixed` \| `σ_def` | 200 | 321.5 | −130.43 | 67 | 2 | **2** | 2 |
| `split_free` \| `σ_def` | 200 | 330 | −130.43 | 180 | 1 | **1** | 1 |
| `free` \| `σ_def` | 200 | 339.5 | −130.43 | **200** | 2 | **2** | 2 |

Which one or two designs the primary repair lands on under each variant is a knife-edge
outcome of a 1 200-iteration L-BFGS run that has not converged
(`derivations/scratch/b4_path/README.md`): the primary-feasible set is optimiser-path
dependent (an earlier pass of the same driver had 3 feasible and 2 positive under `free` and
1 under `split_free`, with a different pair of designs), while every `t = 0` column, the
split-only counts and the dimension identity below are not.

**Verdict: FAIL — 2 / 400 = 0.50 % under `free`, the same rate as the `fixed` control,
against a 10 % bar. B4 is dead as a rescue of the design space.**

### Three things the run nevertheless establishes

**1. The control reproduces K6 exactly.** Under `fixed` the primary repair is 0⁺-feasible
on exactly two designs, `voronoi_12` and `voronoi_93` under `σ_def` — the same two named in
the K6 section above, with the same `Θ_max` (0.002354 and 0.238891 here against 0.002354
and 0.239603 in K6; the primary rows of the two drivers share `X0` and the shape space, and
differ only by the repair's floating-point path). The B4 harness and K6 agree design for
design.

**2. B4's dimension count is exactly right.** `dim_null(free) − dim_null(fixed) = |V_∂| − 1`
on **400 / 400** rows (mean gain 32.8; the `−1` is the pin). Freeing the border really does
buy the `|V_∂|` scalar handles inverse's I3 counts, and `split_free` buys a mean of 12.4 of
them. So the *mechanism* B4 proposes is present and measurable; it is its *consequence* that
fails.

**3. The boundary really does constrain the 0⁺ cone — and that buys almost no range.** The
split-only repair, which K6 used as its relaxation, goes from `59/200` to `176/200` under
`σ_mc` and from `67/200` to `200/200` under `σ_def` when `B` is dropped, monotonically
through `split_free` in between. That is a factor of 3 on `σ_mc`, and on `σ_def` the
relaxed problem becomes **universally** solvable. The median `min q` at `t = 0` also
improves under `σ_mc` (−77.8 → −51.5), i.e. the projection itself lands in a less jammed
place. None of it survives contact with the vertex-edge margin: the primary objective is
feasible on 2 of 200 at best, and the binding contact simply moves, `split-inward` → 185,
`vertex-edge` → 12, `inverted` → 1 under `free`/`σ_def` — the same migration K6 reported.
**The first-order cone is not emptied by the boundary rows; it is emptied by the number of
simultaneous sign conditions, and `|V_∂| ≈ 33` extra dimensions against a median of 126 (`σ_mc`) to 321 (`σ_def`) split-sign
constraints is not enough to change that.**

### The positives, all refereed

| design | variant | `F` | `n_split` | `dim_null` | `min q` | exact `Θ_max` | certified `eps_max` | bisection referee |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| `voronoi_96`, `σ_def` | `free` | 497 | 866 | 947 | 0.1640 | **0.24882** | 0.24882 | **0.24882** |
| `voronoi_12`, `σ_def` | `free` | 521 | 916 | 995 | 0.1058 | 0.07594 | 0.07594 | 0.07594 |
| `voronoi_12`, `σ_def` | `split_free` | 521 | 916 | 960 | 0.0302 | 0.02128 | 0.02128 | 0.02128 |
| `voronoi_93`, `σ_def` | `fixed` | 635 | 1124 | 1124 | 0.1957 | 0.23889 | 0.23889 | 0.23889 |
| `voronoi_12`, `σ_def` | `fixed` | 521 | 916 | 916 | 0.0035 | 0.00235 | 0.00235 | 0.00235 |

Every primary-feasible design is agreed by all three independent labels — the exact T4.2″
scan, the certificate, and the polygon-overlap bisection: `|Θ_exact − bisection| ≤ 1e−5` on
**5 / 5**, worst gap `2.0e−8`. (An earlier pass of this experiment, refereed with the pre-F34
overlap predicate, had one `free` positive contradicted by its referee; with the predicate
fixed, `results/core_validation/referee_fix.md`, no such contradiction remains.) Freeing the
border buys `voronoi_96` outright and lifts `voronoi_12` from 0.002 to 0.076 rad; it loses
`voronoi_93`, which the primary repair reaches only under the fixed boundary.

The `free` positives sit **10.3 and 10.7 median edge lengths** from `X0`, with `min_area > 0`
and `POS` true, so they are embedded flat states — but they are as far from the Eq. (6)
projection as K6's were, which is B4's own risk (8) coming true: where the free boundary
helps at all, it helps by leaving the neighbourhood of the projection entirely, and there is
nothing holding the border in place while it goes.

### What this does and does not change

F30 needs a qualifier and nothing more. "The random-graph design space is empty" was
measured under 2026 Eq. (4)'s fixed boundary; with the boundary rows dropped the rate is
`0.50 %` by the exact scan and refereed, against `0.50 %` with them. So the sentence "the
published boundary condition empties the space" is **not** licensed — the space is
essentially empty either way, and the boundary condition accounts for a factor that is
invisible next to the gap to K7's `4/12` on tori. The torus advantage is therefore **not** explained by the
absence of `B` alone, which is the substantive negative result here and is a constraint on
B1/B2/B3: whatever makes the periodic case work, the `|V_∂|` border handles are not it.

**Consequence.** D8 (the constructive programme moves to periodic patterns) stands. B5,
which `ideas/round2_theorist_b.md` explicitly gates on B4, should be run on the periodic
population rather than on the fixed-boundary one.

No core or method source was modified for this run; `Kirigami/apps/kill_b4.jl` builds its
variant systems from `assemble_system(..., BoundaryMode::None)` and appends the pin rows
itself.

Artifacts: `results/kill/b4/{b4.csv, shard_*.csv, summary_shard_*.txt, cache/}`. Driver
`Kirigami/apps/kill_b4.jl`.

```
julia --project=Kirigami Kirigami/apps/kill_b4.jl --n 200 --nshards 12 --shard I --out results/kill/b4 \
    --cache results/kill/b4/cache        # ~75 min wall with 12 shards concurrent
```

---

## K8a — expansive-cone LP

**Verdict: FAIL on the PASS rule (0 of 100, bar was ≥ 20 of 100), PASS on the hard
soundness control (4 of 4), PASS on the Euler-step sanity (72 of 72, bar was ≥ 10).**
Driver `Kirigami/apps/kill_k8a.jl`, method `Kirigami/src/method/expansive_cone.jl`,
artifacts `results/kill/k8a/` (`k8a.csv`, `summary.txt`, `summary_*.txt`), 5 new unit tests in
`Kirigami/test/test_method.jl` (suite now 187 test sets / 176,692 assertions (23 marked broken)).

*Reproducibility note.* The combinatorial columns of `k8a.csv` (row counts, `dim ker A`,
`dim flex`, `σ ∈ P(X)`, the margin **signs** and the Euler-step outcomes) are exact
functions of the input and reproduce bit for bit. The LP-layer values — margin magnitudes,
dual bounds, the `dual_max < 1e−6` tallies — are the endpoints of an iterative solver
(away-step Frank–Wolfe, 100 000 dual steps; see the recheck below) and move at the last
digit or two between floating-point environments; the σ-chart bound on a non-deployable
`X_ini` and the margin magnitudes on the truncated-square family additionally depend on a
degenerate-row normalisation and are not reproducible across environments at all. The
numbers below are those of the canonical `k8a.csv`.

### What was solved

A flex of the flat structure gives every face an angular velocity and a translation,
`V_f(y) = ω_f J y + w_f`, and the hinge edges impose `V_f(p_i) = V_g(p_i)` at the pin —
which is exactly `mobility.jl::build_rigidity`. Because every copy of a source vertex
sits at `X_v` at `θ = 0` and `X` is held fixed, every `0⁺` separation condition of
`zero_plus.jl` becomes a **strict linear inequality in the flex**:

```
split edge e :  row_e = det( V_{f1}(X_v) - V_{f0}(X_v), d_e )   = q_e / 2
corner (p into the corner of f at v) :  -g1 = det(dV, e1),  -g2 = -det(dV, e2),
                                        mu = max(-g1,-g2) convex, min reflex,  = mu_0+/2
```

The two identities `row_e = q_e/2` and `max/min(-g1,-g2) = mu/2` are asserted exactly (to
`1e-9` relative) against `zero_plus_q` and `zero_plus_corner_margin` in the unit test *"the
cone rows reproduce zero_plus's q_e and corner margins at the sigma ray"*, on three
tilings. The flex basis is built as `ker A` (dense SVD of `build_A`) plus the integrated
translations plus one free translation pair per component of `Γ`; the unit test checks
`dim ker R = dim ker A + 2 c(Γ)` against an independent dense column-pivoted Householder
QR rank of `build_rigidity`, checks orthonormality, and checks that the uniform ray lies in
the computed span. Over all 277 configurations run, `max |R N| / max|R| = 1.2e-13`.

**Prior work, per `notes/screen_r2.md` A1 (verdict PARTIAL).** The expansion cone as *a
polyhedral cone of infinitesimal motions decided by an LP* is Rote–Santos–Streinu 2003
(Lemma 3.2) and Connelly–Demaine–Rote 2003. Nothing here is a new *kind* of object; what
is instantiated is that construction on a **body-and-pin framework with a cut structure**,
where the rows are per split edge and per corner incidence rather than per point pair.
The header of `expansive_cone.jl` carries this citation.

### The LP, and the two-sided bracket

With rows normalised to unit norm in the basis coordinates the quantity solved for is the
matrix-game value

```
val = max_{||z||_2 <= 1} min_i a_i . z  =  min_{lambda in simplex} || A^T lambda ||_2 .
```

Both sides are `>= 0` (take `z = 0`), and `val > 0` iff the strict system is feasible. The
solver is local code, no LP library: a log-sum-exp smoothing annealed from `mu = 1` to
`1e-5` with projected gradient on the unit ball (primal), plus Frank–Wolfe on
`min_lambda ||A^T lambda||_2` over the simplex (dual, 100 000 steps). **Every iterate of
either loop is feasible for its side, so the reported `margin_l2` is a rigorous lower
bound and `dual_bound` a rigorous upper bound on the LP value** — the solver never has to
be trusted. A small dual bound is an approximate Farkas certificate. The three hand-solved
LPs of the unit test *"cone_lp solves three hand-solved linear programs"* pin the solver:
`{z1 > 0, -z1 > 0}` (value 0, `lambda = (1/2,1/2)` exactly), two orthogonal rows (value
`1/sqrt 2`), and `{z1 > 0, z2 - z1 > 0}` (bracket tight to `1e-2`).

**Deviation from the spec.** `ideas/ranking_r2.md` asks for `||z||_inf <= 1` and
multiplicative weights on the `l1` dual. The cone is homogeneous, so the two
normalisations are positive iff each other; the `l2` game is better conditioned and its
dual bound is a by-product of the same gradients. The spec's quantity is reported anyway
as `margin_inf` for the returned witness.

**The convex-corner disjunction, and the honest asymmetry.** Every one of the 244 931
corner incidences on the K1a population at `X_ini` is at a **convex** corner, so `P(X)` is
a union of `2^{n_convex}` polyhedra and no enumeration is possible. Four passes are run:
pass 0 is the conjunctive relaxation (a *subset* of `P(X)`, needing no branch choice at
all), pass 1 is the branch the uniform ray itself selects, passes 2–3 repair the branch
from the best witness so far. Fixed points were reached at pass 2 or 3 on every graph.
**Feasibility of any pass proves `P(X) != {}`. Infeasibility of the passes tried proves
nothing** — what supports the negative reading is the dual bound, not the failure to find
a witness.

### Hard soundness control: 4 of 4

| tiling | F | split | dim ker A | dim flex | σ ∈ P(X) | min `q_e` | min `mu` | LP margin |
|---|---|---|---|---|---|---|---|---|
| `hexagons_auto` | 10 | 5 | 3 | 5 | yes | 1.000 | 1.000 | 6.21e−1 |
| `truncated_square_488` | 21 | 12 | 2 | 4 | yes | 1.414 | 1.414 | 9.25e−1 |
| `snub_square_33434` | 53 | 13 | 29 | 31 | yes | 0.688 | 0.509 | 6.58e−2 |
| `tiling_3_4_3_12` | 25 | 4 | 5 | 7 | yes | 1.742 | 1.216 | 3.12e−1 |

σ is reported **inside** `P(X)` on all four, so the hard soundness kill did not fire, and
the LP independently finds a strictly positive margin on all four. (`tiling_3_4_3_12` is
evaluated at its Eq. (6) projection, the others at `X_ini`, following the rule "X_ini if
Eq. (2) holds there, else the projection" — the same rule the method tests use.)

Extended control, `kill_common.jl::deployable_population` (frozen in `data/corpus/`) (authored tilings at five clip
radii plus shape-space samples, 73 configurations): σ ∈ `P(X)` on **71 of 73**, LP margin
> 0 on **72 of 73**.

### The result on the K1a population

| | at `X_ini` | at `X0` |
|---|---|---|
| LP margin > 0 | **0 / 100** | **0 / 100** |
| `dual_max < 1e-6` (every branch chart tried certified) | 95 / 100 | 99 / 100 |
| `dual_max < 1e-5` | 99 / 100 | 100 / 100 |
| `dual_max` median / max | 6.7e−10 / 1.3e−5 | 6.7e−10 / 1.7e−6 |
| σ is a flex of the framework | 0 / 100 | 100 / 100 |
| σ ∈ `P(X)` | 0 / 100 | 0 / 100 |
| dim ker A, median | 17 | 65 |
| rows, median | 4 911 (142 split, 2 338 corner incidences) | same |

Per family at `X_ini`: Voronoi `dim ker A` median 3 (34/34 dual-certified), quad-random
median 17 (33/33 certified), Delaunay median 72 (28/33 certified, `dual_max` median
4.9e−9). The Delaunay residue is Frank–Wolfe convergence, not a different phenomenon: at
20 000 dual steps the same graphs sat at `1e−4`, at 100 000 they sit at `1e−6` or below,
and a per-graph test at `dim_flex = 156` walked `1.2e−4 → 4.5e−5 → 1.8e−6` over
600 / 4 000 / 100 000 steps. No Delaunay graph has a positive primal margin at any budget.

**Reading.** On this population the emptiness is **not** "one specific ray of a
305-dimensional cone points the wrong way". The whole cone, over every branch chart the
driver reached, is empty: on 95 of 100 graphs there is a non-negative combination of the
split-edge and corner rows summing to within `1e-6` of zero, which is an (approximate)
Farkas certificate that no flex at all — uniform or not — separates every cut and every
corner at first order. That is the **structural** reading of F25/F30 and is, per
`ideas/ranking_r2.md` §3 #1, the stronger defence against reviewer point R-1: the failure
is a property of the cut structure, not of the two orientation heuristics.

`σ ∈ ker A` on 100/100 at `X0` and 0/100 at `X_ini` reproduces K3a's split (`296/296` at
`X0`, `0/212` at `X_ini`) exactly, from independent code. At `X_ini` the uniform
deployment is not even a motion of the framework, so "σ ∉ P(X_ini)" there is a statement
about a vector, not about a deployment; the meaningful σ row is the `X0` one, where σ *is*
a flex and is nevertheless outside `P(X0)` on 100/100 (median `min q_e = -79.2`, median 66.5
split edges with `q_e <= 0`).

### The one non-uniform rescue that does exist

`snub_square_R20_s1` (a shape-space sample of the snub-square tiling) has σ **outside**
`P(X)` — `min q_e = -6.05e-2`, so the uniform deployment collides at `0⁺` and `Θ_max = 0` —
and the LP nevertheless certifies a **strictly positive non-uniform margin of 5.93e-2**,
with a collision-free Euler step at the largest tested amplitude. So the mechanism A1
proposes is real and does fire; it simply does not fire on any of the 100 random graphs.
`trunc_square_R40_s0` is the control for that: σ outside `P(X)` and the LP finds nothing.

### Euler-step sanity (task step 2)

On every configuration with a positive margin the witness flex is integrated one explicit
Euler step, `Y_p = X_v + h V_{face(p)}(X_v)`, at `h = eps * (median edge) / max_p |V_p|`
for `eps` descending through `{0.5, 0.2, 0.1, 0.03, 0.01, 0.003, 0.001}`, and the stepped
structure is tested for face–face overlap. **72 of 72 are collision-free**, 65 of them
already at `eps = 0.5` (half a median edge of vertex travel), 4 at 0.2, 2 at 0.1, 1 at
0.03. All four reference tilings reach `eps = 0.5`.

*Deviation, stated:* `contact.jl`'s T4.2″ scan is parameterised by the **uniform** angle
`θ` and has no meaning along a non-uniform flex, so it cannot be used here. The predicate
used is `has_collision` / `polygons_overlap` — the same exact overlap primitive T4.2″
itself calls at every interval midpoint — and a unit test checks that along the uniform ray
the Euler step agrees with the exact kinematics `deploy(c, X, h)` to `1e-3` median edges.

### Reproducing

```bash
julia --project=Kirigami -e 'using Pkg; Pkg.test()'
for i in $(seq 0 11); do julia --project=Kirigami Kirigami/apps/kill_k8a.jl --n 100 --x0 --dual-iters 100000 \
    --out results/kill/k8a --cache results/kill/cache --shard $i --nshard 12 \
    > results/kill/k8a/log_$i.txt 2>&1 & done; wait
cd results/kill/k8a && head -1 k8a_0.csv > k8a.csv \
  && for f in k8a_[0-9]*.csv; do tail -n +2 $f; done >> k8a.csv && cd -
julia --project=Kirigami/scripts Kirigami/scripts/summarise_k8a.jl results/kill/k8a \
  > results/kill/k8a/summary.txt
```

Every graph is a deterministic function of its id (`kill_common.jl::make_graph(id, 100, 800,
1400)`), so shards may be run in any order and a single graph re-run with
`--nshard 100 --shard i`. Wall time about 11 min of solver time over 12 shards.

### Caveats, stated

1. **Infeasibility is not proved, it is certified to `1e-6`.** The dual bound is a
   numerical Farkas residual, not an exact rational certificate, and it is computed per
   branch chart. 5 of 100 graphs at `X_ini` have `dual_max >= 1e-6` (all Delaunay, all
   below `1.3e-5`); on those the FAIL rests on "no witness found", which is weaker.
2. **The branch enumeration is not exhaustive** and cannot be: `2^{2338}` charts per
   graph. Four charts are visited. A positive result would have been a proof; the negative
   result is evidence.
3. **`dim ker A` medians here (Voronoi 3, quad-random 17, Delaunay 72 at `X_ini`) are not
   the F26/F29 number.** F29's median 305 is `dim ker A` on the **2-core** of a different
   K3a population; these are the full graph at a different size distribution. The two are
   not comparable and neither contradicts the other, but the gap is large enough to be
   worth one explicit check before either is quoted next to the other.
4. The margin is scale-dependent through the row normalisation and the orthonormal flex
   basis; only its **sign** is invariant. No absolute margin value should be compared
   across graphs of different size.

### Recheck (K8a-recheck agent) — the solver was wrong; the certificates were not

Full write-up `results/kill/k8a/recheck.md`; artifacts `results/kill/k8a/recheck_k9.csv`,
`recheck_k9_summary.txt`, `results/kill/k8a/recheck/` (20 000 dual steps, controls
included) and `results/kill/k8a/recheck100k/` (100 000 steps, K1a only). Driver
`Kirigami/apps/kill_k8a_recheck.jl`, 4 new unit tests. **Nothing above is deleted; every number
in this subsection is a re-measurement.**

**The bug.** `cone_lp` computed the two ends of its bracket with two unrelated algorithms.
The dual end (Frank–Wolfe on `min_λ ||A^T λ||`) was sound. The primal end came from a
separate log-sum-exp smoothing loop — about 1 000 projected-gradient steps of size
`μ/σ_max² ≈ 1e−5`, started at `z = 0`, in a space of dimension up to 300 against 5 000 rows
— and was the **only** source of the returned witness. The Frank–Wolfe iterate `w = A^T λ`
was never evaluated as a primal direction, although it *is* the optimum: at the minimum-norm
point of `conv{a_i}` one has `a_i · w ≥ ||w||²` for every `i`, so `z = w/||w||` attains
`min_i a_i · z ≥ ||w||` and the bracket closes at `margin = val = ||w*||`.

**The fix.** Away-step Frank–Wolfe (Lacoste-Julien & Jaggi 2015 — linearly convergent on a
polytope) on the min-norm point, with the primal witness read off the same iterate; the
smoothing loop retained only as an extra witness source, so nothing previously found is
lost. New `ConeLPResult::gap`, new `ExpansiveConeReport::sigma_chart_margin` (the uniform
ray on its own branch chart, computed with no solver in the loop — a rigorous lower bound
the solver must beat), and `farkas_residual()`, which recomputes `||A^T λ||` row by row from
the multipliers alone and checks `λ ≥ 0`, `Σλ = 1`.

**Was it wrong? Yes — decisively.** On K9's 30 variant-(b) embeddings, regenerated
bit-identically:

| | old solver | corrected |
|---|---|---|
| `σ ∈ P(X)` by direct `zero_plus` measurement | 23 | 23 |
| ... LP margin `> 0` | **0 / 23** | **23 / 23** |
| ... margin ≥ the solver-free `sigma_chart_margin` | — | **23 / 23** |
| ... a non-uniform flex strictly beats `σ` | — | **23 / 23** |
| `σ` **outside** `P(X)`, LP margin `> 0` anyway | **0 / 7** | **6 / 7** |
| feasible overall | **0 / 30** | **29 / 30** |

`margin / sigma_chart_margin`: median **8.73**, max 56.9. K9's own summary now carries the
corrected solver (`expansive-cone LP: run 30, feasible 29`); the `feasible 0` lines of the
original K9 pass were **solver artefacts and must not be quoted**.

**Were the certificates real? Yes.** Every dual multiplier vector recomputes as a genuine
simplex point: over all 277 configurations the worst `min_i λ_i` is exactly `0` and the
worst `|Σλ − 1|` is `2.5e−13`. Independently verified counts, 100 000 dual steps, "every
branch chart the driver solved":

| | at `X_ini` | at `X0` |
|---|---|---|
| reported `dual_max < 1e−6` (original run) | 84 / 100 | 68 / 100 |
| **verified** `||A^T λ|| < 1e−6` | **96 / 100** | **99 / 100** |
| **verified** `< 1e−9` (the tolerance asked for) | **65 / 100** | **57 / 100** |
| verified `< 1e−5` | 98 / 100 | 100 / 100 |
| worst verified residual | `3.9e−5` | `1.8e−6` |

Per family at `X_ini`, verified `< 1e−9`: Voronoi **34 / 34**, quad-random **31 / 33**,
Delaunay **0 / 33** — 29 / 33 at `1e−6` (`dim_flex` up to 156; residuals `3e−9 … 4e−5`, two to three orders
better than the original vanilla Frank–Wolfe at the same budget, but short of `1e−9`).

**New margin counts.** LP margin `> 0` on the K1a population: **0 / 100 at `X_ini`, 0 / 100
at `X0`** — identical to the original, so **the FAIL verdict and the 0 / 100 headline
stand**. The corrected solver finds no witness the old one missed on this population. The
controls improve, which is the second signature of the old primal's under-convergence:
`hexagons_auto` `5.40e−1 → 6.21e−1`, `truncated_square_488` `7.35e−1 → 9.25e−1` (this
family's magnitude is the one that depends on the degenerate-row normalisation; the recheck
pass itself printed `8.52e−1`), `snub_square_33434` `3.53e−2 → 6.58e−2`, `tiling_3_4_3_12`
unchanged at `3.12e−1`, all with brackets closing to `1e−9` or better. Hard soundness 4 / 4
and deployable population 72 / 73 margin, 72 / 72 Euler — unchanged. The main table of this
section is the canonical `k8a.csv`, produced with the corrected solver at 100 000 dual
steps; the recheck tables above are the archived `recheck100k/` pass.

**What survives, and what a certificate means.** A verified residual `r` with `margin = 0`
pins the chart's LP value to `[0, r]`. That is a numerical bound, not an exact rational
Farkas certificate, and **caveat 2 above is untouched and remains binding**: four charts of
`2^{2338}` are visited, so the certificate is of *those charts*, never of `P(X)`. The
recheck's contribution is that the evidence is now evidence rather than a solver failure —
and that on the one population where `P(X)` is provably non-empty (K9's 23), the same four
charts find a witness immediately.

**Test suite.** 4 new unit tests in `Kirigami/test/test_method.jl`: the bracket closes to `1e−9`
on two unit rows at a known angle (`val = cos(φ/2)` exactly, four angles); a 1 200-row,
60-dimension narrowly feasible system where the old solver returned 0 and the new one beats
the explicit witness `z = e₀`; `farkas_residual` on hand certificates including one off the
simplex; and `margin_l2 ≥ sigma_chart_margin` on the four split-bearing tilings. Suite
passes (103 cases / 17 504 assertions at the time of the run).

## B3 — tr K as the expansion budget

`Kirigami/apps/kill_b3.jl`, library `Kirigami/src/method/budget.jl` (5 new unit tests in
`Kirigami/test/test_method.jl`; suite 187 test sets / 176,692 assertions (23 marked broken), 0 failures). Outputs
`results/kill/b3/{b3_b1_global.csv, b3_retro.csv, b3_constrained.csv, summary.txt}`.

**Verdict — theorem PASS after a repair to its derivation, algorithm FAIL, and B3's own
hard-fail clause TRIGGERED.**

**B1 global half (the part `derivations/scratch/check_b1.jl` could not close).** On the
same 10 fixed-boundary patches, the total void area — the shoelace of the deployed border
half-edges minus the θ-independent face-area sum, so holes *and* notches with no walk
tracing — equals the border functional's first harmonic to `3.41e-15` at 6 angles, and
the edge-wise budget satisfies `W + R = 2 B(X)` to `2.92e-16`. Bar `1e-10`: **PASS**.
Notches are not a correction: `Σ_C a_C` is 28–46 % *below* the hole-only sum on every
pattern, i.e. notches close as holes open.

**B3's identity is wrong as written, and the repair changes what B3 says.** Both copies of
an edge contribute `⟨x_a − x_b, u_{f1} − u_{f0}⟩`. For a split edge `σ_{f0} = σ_{f1}` so
the term is lattice-translation invariant, but for a hinge edge `σ_{f0} = −σ_{f1}` and
`u_{f1+t} − u_{f0+t} = u_{f1} − u_{f0} + σ_{f1} t`, so the term depends on which **lift**
the preimage walk uses. Summing one arbitrary representative per quotient edge misses a
defect measured at `0.38–1.88 × det P₀` on the 33 K7 patterns. Summing one representative
**preimage** per translation class, each at its own lift (`method::periodic_cell_edges`),
gives `W + R = det P₀ · tr K` to `3.9e-14` on 46 of 50 C3 rows — every row of 23 of the 25
patterns with `dim 𝒦 ≥ 1` (`b3_retro.csv`, `b5_rel`)
(exceptions: `squares_3x2`, where the 3×3 super patch carries no full representative set,
and `snub_square_3x3`, whose `K` is gauge-dependent — F33). Consequence:
`tr K − τ* = R / det P₀` **exactly**, so the budget gap *is* the total split budget
rescaled, `τ*` is a function of the design rather than of the pattern, and
`tr K ≥ τ*` is not a linear constraint on the achievable set.

**Stage 1, retrodiction (PASS bar 1).** The baseline reproduces `k7_c3_all.csv` bit for
bit: 0 mismatches on `hit_err/cert/eps_cert/min_q/min_area/zp_feasible` over all 50 rows,
13 certified, 21 with `min q > 0`.

| quantity | value | bar |
|---|---|---|
| AUC(`tr K − τ*`) for `min q > 0` | **0.989** | 0.80 PASS |
| AUC(`tr K − τ*`) for certified | **0.906** | 0.80 PASS |
| rows with `min q > 0` at `tr K < τ*` | **0** | must be 0 PASS |
| same, restricted to the 46 identity-exact rows | 0.987 / 0.895 | — |
| AUC(`R/(det P₀ n_split)`) for `min q > 0` / certified | 0.982 / 0.963 | — |
| **spread(W) / spread(det P₀ tr K) over 𝒦**, 20 points × 25 patterns | **median 0.776, max 1.401** | **HARD FAIL** |

The last row is the clause B3 §8 and the Critic both name as fatal: `W` moves over the
achievable set by ~0.78 of what `det P₀ tr K` moves, so `τ*` is a penalty, not a threshold.
(The random-target rows inherit K7 C3's basis-dependent draw, so the AUCs and the certified
counts of this section move with it between floating-point environments; the hard-fail
ratio and the identity checks do not.)

**Stage 2, constrained re-solve (FAIL bar 2).** Minimising `‖K − K*‖` over 𝒦 subject to
`tr K ≥ τ* + δ|τ*|`; since `g = Aᵀ vec(I)`, at `dim 𝒦 = 4` the minimiser is
`K* + ((τ*(1+δ) − tr K*)/2) I`, an isotropic expansion added to the target, and the
leftover freedom is the same null space the baseline searched, spent identically. `τ*` is
refreshed by a 4-step fixed point because it depends on the design.

| δ | certified | flips | `min q > 0` | max ‖K − K*‖ |
|---|---|---|---|---|
| 0.0 | 13 → **15** | +2, −0 | 21 → 37 | 10.4 |
| 0.1 | 13 → **14** | +2, −1 | 21 → 37 | 19.8 |
| 0.5 | 13 → **2** | +1, −12 | 21 → 25 | 588.9 |

`+2` and `+1` of 50 are inside the 50-row binomial error bar (±4.7 at `p = 0.26`). The
constraint does what the theorem says — it nearly doubles the designs whose split cuts all
open at `0⁺`, 21 → 37 — but that does not convert into certified deployment range, which
is the K5/K6/F30 obstruction again and is not a budget phenomenon. And the target is
abandoned: the isotropic correction is 10–589 in Frobenius norm.

**Contradictions with STATE.md: none.** K7's C3 row is reproduced exactly with the
repaired certificate, and F33's `snub_square_3x3` gauge caveat is the cause of one of the
two identity exceptions.

---

## K9 — convexity-constrained embedding

Added after the round-2 escalation, as the constructive half of decision D9. It is the first
experiment on random graphs that returns a **non-zero** deployment range.

**The hypothesis.** T-1 measured that 1 969 of 1 970 balanced pure vertices with a
non-positive `0⁺` margin sit at a **reflex** corner, and that the input embeddings have
**zero** non-convex faces while `X0` has thousands of non-convex corners. At a convex corner
the margin is `μ = max(−g₁, −g₂)`, a disjunction; at a reflex corner it is `min(−g₁, −g₂)`, a
conjunction, and that is the one that fails. So the reflex corners are **created** by the
unconstrained projection of Eq. (6), not by the graph, and the question is whether putting
them back inside the shape space buys deployment range that K6's `q`-only repair could not.

**The problem solved** (`Kirigami/src/method/convex_embed.jl`, derived there in full).
Inside the affine shape space `X(t) = X0 + Φ t`,

```
minimise ‖X(t) − X_ini‖²_F   subject to   cross_i(t) ≥ δ   for every corner i of every face
                                          [ and q_e(t) ≥ δ′ for every split cut ]
```

with `cross = det(X[v] − X[v_prev], X[v_next] − X[v])`, positive exactly at a strictly convex
corner of a CCW face. Every corner strictly convex implies every face convex **and**
positively oriented, so this constraint **subsumes** the face-area constraint of K6's repair.
Both families are quadratic in `t` and the feasible set is not convex, so the solve is a
heuristic and is reported as one: phase A minimises a softplus penalty under a continuation
over the target `τ` (`0.01, 0.1, 1` × the solve target, warm-started); phase B refines the
proximity behind a log barrier, keeping **only strictly feasible** iterates. FEASIBLE always
means the exact constraint values at the returned point, never the value of a penalty.

**PASS rule:** refereed exact `Θ_max > 0` on **≥ 10 %** of the 400 designs (≥ 40).

### How it was run

* **Population.** Exactly K6's 200 graphs under both orientation rules — 400 designs, the
  same `σ_def` read back from `results/kill/k5/sigma/*.json`. Nothing is re-sampled, so K9,
  K6 and K5 are the same designs measured three ways.
* **Settings.** `δ = δ′ = 1e−3·med²`, 3 Gaussian starts of scale `0.1·med` plus `t = 0`,
  L-BFGS 600 iterations per continuation stage, 6 barrier stages, certificate `ε = 0.3`.
* **Three arms on every design.** `X0` (no solve, the control), **(a)** convexity only,
  **(b)** convexity + split-inward. (a) is a relaxation of (b), so (a)'s minimiser is handed
  to (b) as a warm start.
* **Sharding.** 12 disjoint shards, merged to `results/kill/k9/k9.csv` (400 rows); 1.60 h of
  solver time across the shards.

### Measured — all from `results/kill/k9/{k9.csv, summary.txt}`

| quantity | `X0` (control) | (a) convexity | (b) convexity + split |
|---|--:|--:|--:|
| convex-feasible designs | 79 / 400 (19.8 %) | **118 / 400 (29.5 %)** | 30 / 400 (7.5 %) |
| **exact `Θ_max > 0`** | **0** | **0** | **36 / 400 (9.0 %)** |
| refereed (bisection, fixed predicate) | 0 | 0 | **36** |
| refereed (bisection, pre-F34 predicate at shrink `1e−12`, earlier pass) | 0 | 0 | 32 |
| certified (`ε_max > 0`) | 0 | 0 | **36** |
| all `0⁺` corner margins `μ > 0` at the returned `X` | 0 / 79 | 0 / 118 | **23 / 30** |
| median `min μ / med²` over feasible | −8.79 | −10.06 | **+4.96e−2** |
| `‖X − X_ini‖` per vertex, median / q90 (median edges) | 0 / 0 | 0.275 / 0.399 | 0.376 / 0.567 |
| expansive-cone LP feasible (corrected solver, §K8a recheck) | — | 0 / 118 | **29 / 30** |
| binding contact at `Θ_max = 0` | split-inward 400 | split-inward 400 | split-inward 259, inverted 73, vertex-edge 32 |

(The (a) arm's convex-feasible count and the (b) arm's endpoints are L-BFGS outcomes and
optimiser-path dependent — an earlier pass gave 115 convex-feasible under (a) — while the
36 positives, their `Θ_max` and the `t = 0` columns reproduce.)

Per orientation, arm (b):

| quantity | `σ_mc` | `σ_def` |
|---|--:|--:|
| convex-feasible | 13 / 200 | 17 / 200 |
| **exact `Θ_max > 0`** | **24 / 200 (12.0 %)** | **12 / 200 (6.0 %)** |
| refereed (fixed predicate; pre-F34 at `1e−12` in brackets) | 24 (20) | 12 (12) |
| `ε_max` median / max | 0.278 / **1.997** | 0.152 / 0.381 |
| `‖X − X_ini‖` per vertex, median | 0.523 | 0.297 |
| expansive-cone LP feasible (corrected solver) | 13 / 13 | 16 / 17 |

Over the 36 positives: `Θ_max` median **0.249**, min 0.027, max **1.997** rad; **33 of 36
have `ε_max ≥ 0.1` rad**; by family delaunay 25, voronoi 10, quad_random 1. Taking the better
`σ` per graph, **33 of the 200 graphs** carry at least one deployable design. Note that only
**23 of the 36** positives are `b_feasible` at the full margin `δ′`: the other 13 miss the
margin the solver was aiming at and deploy anyway, so the margin is a sufficient handle on
the `0⁺` cone, not a necessary one.

**Convexity alone is necessary, not sufficient.** Arm (a) raises convex feasibility from
79 to 118 designs and buys **exactly zero** deployment range — the binding contact stays
`split-inward` on all 400. The two constraints together are what moves the count off zero,
which is the same lesson K6 reached from the other side (split-only feasible on 126 / 400,
certified 15). Non-convex corners at `X0`: **35 336 / 656 736 = 5.38 %**.

### The referee audit — resolved by the F34 predicate fix

Four of the 36 exact positives (delaunay 13, 40, 103, 190, all `σ_mc`) refereed as
`Θ_bisect = 0` at the historical shrink `1e−12`. `Kirigami/apps/dbg_k9.jl` traces every one to
the bisection's early-exit `has_collision(1e−7)` firing on a pair of **hinge-adjacent** faces,
which touch at the pin **by construction** — and the firing is **non-monotone in `θ`**: true
at `1e−8, 1e−7, 1e−5`, false at `1e−9, 1e−6, 1e−4, 1e−3, 1e−2, 0.05, 0.1`. It is numerical
noise of the near-flat state, not a contact: the shrink displaced the two copies of the shared
pin by different amounts. `polygons_overlap` has since been rewritten (F34 fix,
`results/core_validation/referee_fix.md`) and **K9 re-run in full against a matched baseline
binary**: the referee now agrees with the exact scan on **400 / 400** designs with worst gap
**0**, the refereed count is **36 / 36**, and the only four cells that changed in the whole
400-row table are these four. Re-refereeing at the old shrink `1e−9` already reproduced the
same values:

| id | exact `Θ_max` | bisection (`1e−9`) |
|---|--:|--:|
| delaunay 13 | 0.219029 | 0.219029 |
| delaunay 40 | 0.177931 | 0.177931 |
| delaunay 103 | 0.301957 | 0.301957 |
| delaunay 190 | 0.103354 | 0.103354 |

This is the F34 artefact, on the same side as K6's `first_root` bug but with the opposite
sign: there the referee caught the scan, here the scan catches the referee. **32** was the
conservative headline used against the PASS bar; **36** is the number the fixed referee gives,
and it is now the number to quote. Neither reaches 40. Soundness (`ε_max ≤ bisection Θ_max`)
held on 396 / 400 before the fix and on **400 / 400** after; arms `X0` and (a) referee
400 / 400 with worst gap **exactly 0**, unchanged by the fix.

### The expansive-cone LP control

As first run, on the 30 (b)-feasible designs the LP of K8a found **no** witness (0 / 30) —
but on **23** of them the uniform ray `σ` is itself in `P(X)` by direct measurement (`min q > 0` **and**
`min μ > 0`, and `expansive_cone.jl`'s rows equal `q/2` and `μ/2` at the `σ` flex, asserted
by unit test). A strictly feasible point therefore **provably exists** on those 23 and the LP
failed to recover it; its dual bounds there run `2.0e−2` to `5.4e−2`, far above the `1e−7`
threshold at which infeasibility would be claimed, so the LP asserts nothing. These `0 / 30`
are **cone_lp non-convergence, not emptiness**, and K8a's caveat 1 applies with force: K8a's
dual "certificates" on 84 / 100 must be re-examined before any of them is quoted.

**Resolved — see §K8a "Recheck".** The diagnosis was right: `cone_lp`'s primal witness came
from an under-converged smoothing loop and never used the Frank–Wolfe iterate, which is the
optimum. With the corrected away-step min-norm-point solver — the one the canonical
`k9.csv` now carries — the same 30 embeddings give **29 feasible**, all 23 with `σ ∈ P(X)`
among them, at a median **8.7×** the uniform ray's margin, plus 6 of the 7 rescues where `σ`
is outside `P(X)`. **The `0 / 30` figure must not be quoted.** K8a's dual certificates were
sound and independently verify at 95 / 100 (`X_ini`) and 99 / 100 (`X0`) at `1e−6`; K8a's
own `0 / 100` FAIL is unchanged.

### Verdict: FAIL — 9.0 % (36 of 400) against a 10 % bar

The bar is missed by **4 designs**, and it is missed while every comparable baseline near
the projection reports zero or nearly so: Eq. (6) alone 0 / 200 (K1a), `σ_def` 0 / 200 and
`σ_mc` 0 / 200 (K5), K6's primary repair 2 / 400 (and its split-only relaxation 15 / 400,
every one by leaving the projection's neighbourhood), B4 2 / 400, K8a 0 / 100, and the
authors' native `prevent` 0 / 80 on the capped comparison. Against that population, 36
designs with a **certified** `ε_max` up to 2.0 rad is the first constructive foothold on
random graphs that stays near the projection, and F30 needs the qualifier it did not need
before: the random-graph design space is empty **under unconstrained Eq. (6)**, and only
sparsely non-empty when convexity and the split signs are imposed together inside the null
space.

Artifacts: `results/kill/k9/{k9.csv, summary.txt, shard_*.csv, shard_*.txt}`. Driver
`Kirigami/apps/kill_k9.jl`, method `Kirigami/src/method/convex_embed.jl`, referee diagnostic
`Kirigami/apps/dbg_k9.jl`, unit tests in `Kirigami/test/test_method.jl` (hand-checked
`corner_crosses`, a hand-known constrained minimiser on the 2×2 grid, the phase-A gradient
against a central difference, the FEASIBLE-is-exact invariant, and barrier-stage
monotonicity).

---

## K9b — pushing the constrained embedding

K9 missed its 10 % bar by 4 designs. K9b holds the population, the shape space and the
certificate **fixed** and pushes only the solver, to decide whether 9 % is a property of the
design space or of the search. The answer is neither of the two the task expected: **more
search made it worse**, and the reason is informative.

**The four levers** (`Kirigami/apps/kill_k9b.jl`), in the order the plan ranked them:

* **(i) per-graph best of `σ`.** Nothing new is computed — both `σ` already run on every
  graph — so this is a reporting convention, always given as "best-of-2 over graphs" beside
  the per-design count and never mixed into it.
* **(ii) more restarts and more barrier stages.** 8 Gaussian starts + `t = 0` (K9: 3 + 0),
  plus the same two warm starts K9 used, variant (a)'s convex minimiser and K6's split-only
  repair point; 10 barrier stages (K9: 6).
* **(iii) a sweep of `δ = δ′` over `{1e−4, 1e−3, 3e−3, 1e−2}·med²`**, tried best-first and
  cut short once a design is certified past `ε = 0.1` or a 200 s per-design budget is spent.
* **(iv) stage 2**, `range_opt.jl`'s softmin-of-first-contact objective (the T6 analytic
  gradients) re-centred at the feasible point, then **gated**: accepted only if the exact
  convexity and split constraints still hold and the exact `Θ_max` actually improved. The
  barriers are enforced by exact rejection, not inside the stage-2 objective — weaker than
  the plan asked for, and stated here as such.

Every design is refereed at **both** shrinks, `1e−12` and the `1e−9` that K9's audit showed
to be free of the hinge-vertex artefact. 12 shards, 7.69 h of solver time, ~62 min wall.

### Measured — `results/kill/k9b/{k9b.csv, summary.txt}` (400 rows)

| quantity | K9 (b) | **K9b sweep** | best over both |
|---|--:|--:|--:|
| exact `Θ_max > 0` | 36 / 400 (9.0 %) | **28 / 400 (7.0 %)** | **37 / 400 (9.25 %)** |
| refereed at shrink `1e−12` | 36 (32 with the pre-F34 predicate) | **28** | 37 |
| refereed at shrink `1e−9` | 36 | **28** | **37** |
| certified (`ε_max > 0`) | 36 | 28 | 37 |
| ... of which `ε_max ≥ 0.1` rad | 33 | 16 | **34** |
| `ε_max` median / max (rad) | 0.253 / 1.997 | 0.147 / 1.828 | 0.253 / 1.997 |
| convexity + split feasible | 30 | 28 | — |
| graphs with ≥ 1 deployable `σ` (of 200) | 33 | 26 | **34** |

Per orientation, K9b: `σ_mc` **22 / 200** positives (11.0 %), `σ_def` **6 / 200** (3.0 %) —
the same asymmetry K9 measured, sharpened. Referee agreement is exact: worst
`|exact − bisection(1e−9)|` **4.1e−06** over all 400, and soundness `ε_max ≤ bisection`
holds **400 / 400**. Sixteen rows (e.g. `delaunay 40 σ_mc`, `bisection(1e−9) = 4.1e−6` rad)
have a bisection value below `5e−6` with the exact scan at 0; that is below the bisection's
own grid resolution and is counted as a zero, not a positive.

### Which lever mattered, and why the count went down

**The `δ` sweep is the only lever that paid.** Among K9b's 28 positives the winning margin is
`δ = 1e−2·med²` on **15**, K9's `1e−3` on **12**, `1e−4` on **1**. A convexity/split margin
ten times K9's is the better setting on more than half the designs it finds, which says the
constraint should be posed with a *generous* margin rather than a tight one.

**Restarts and barrier stages did not pay, and stage 2 barely ran.** Stage 2 supplied the
accepted point on **2** of 28 positives. Worse, it was gated on `Θ_max > 0`, so it never ran
on the feasible-but-jammed designs where it would have mattered — a design error in this
driver, not a measurement, and the single clearest next step.

**The 9 designs K9b lost say what the real obstruction is.** K9b keeps 27 of K9's 36
positives, adds 1 (`delaunay 22 σ_def`, `Θ_max = 0.151`), and loses 9. All nine are
`delaunay`, and on all nine K9b's returned point is **exactly feasible** — every corner cross
and every split sign clears the margin — and still has `Θ_max = 0`, binding on `vertex-edge`.
So the two solvers found *different points inside the same feasible set*, one deployable and
one not. **Margin feasibility does not determine deployability**, and the proximity objective
`‖X − X_ini‖` that selects the point inside the feasible set is not aligned with deployment
range. That, and not the search budget, is what caps the count near 9 %: only 13 of K9b's 28
positives are margin-feasible, and only 13 of its 28 feasible designs are positive.

Because a positive is exact-verified under the same certificate on the same design, the
per-design maximum over the two solver configurations is a legitimate multi-start count, and
it is reported as such and never in place of either run: **37 / 400 refereed at `1e−9`, all
37 certified, 34 with `ε_max ≥ 0.1` rad, on 34 of the 200 graphs.**

### Verdict: FAIL — the push does not close the gap

| target | achieved |
|---|---|
| ≥ 40 / 400 refereed positive | **28** (K9b alone), **37** (best over both solvers) |
| as many as possible with `ε_max ≥ 0.1` rad | **16** (K9b alone), **34** (best over both) |

The bar is still missed, now by 3 designs rather than 4, and the baselines on this exact
population that stay near the projection are still **0** — Eq. (6) alone, `σ_mc`, `σ_def`,
convexity alone, K8a, and the authors' native `prevent` — with K6's repairs and B4 at 2 / 400
near the projection and K6's split-only relaxation at 15 / 400 by leaving it. The
constructive claim that survives is unchanged in kind and firmer in evidence: **37 random
graphs deploy from a constrained embedding near the projection, against 2 without it.**

**Figures.** `results/kill/k9b/k9b_gallery.png` — **24** random graphs, each closed and at
`θ = Θ_max/2`, all from K9b's own run (`F` from 101 to 442, `Θ_max/2` from 0.0002 to 0.914
rad). `results/kill/k9b/k9b_hist.png` — the `ε_max` distribution and the count against every
baseline at zero.

**Contradictions with STATE.md: one.** F36 records K9 as "9.0 %, first non-zero constructive
result", which stands. It does not record that the count is **basin-selected rather than
budget-limited**; K9b's 9 lost-and-still-feasible designs are the evidence, and F36 should be
amended to say that increasing the search budget by roughly 4× moves the count *down*.

Artifacts: `results/kill/k9b/{k9b.csv, summary.txt, shard_*.csv, shard_*.txt, gallery/,
k9b_gallery.png, k9b_hist.png}`. Driver `Kirigami/apps/kill_k9b.jl`, plots
`Kirigami/scripts/plot_k9b.jl`.

---

## K9c — range-maximising embedding

### The finding K9c acts on

K9b's audit is the premise. Of K9's 36 deployable designs, **9** come back from K9b's
four-times-larger search as points that are **exactly feasible at the full margin** — every
corner cross and every split sign clears `δ` — and still have `Θ_max = 0`, binding
`vertex-edge`. Two solvers found *different points inside the same feasible set*, one
deployable and one not. So the feasible set is not what caps the yield near 9 %; the
**objective that selects the point inside it** is, and that objective is proximity,
`‖X − X_ini‖`, which knows nothing about deployment.

K9c keeps the population (K6's 200 graphs × 2 orientation rules = 400 designs), the shape
space (Eq. (6)'s null space `X = X0 + Φt`), the barriers (`cross_i ≥ δ`, `q_e ≥ δ′`) and the
certificate fixed, and replaces only the objective.

### What is maximised

`Kirigami/src/method/range_embed.jl` derives it in full. The scalar is the **0⁺
deployment margin**, in units of `med²`,

```
m(X) = min( min_e q_e(X) , min_j mu_j(X) ) / med^2
```

— the joint minimum of `zero_plus.jl`'s two first-order separation families, the split-edge
sign `q_e = det(dS_e, d_e)` and the corner (vertex-into-edge) margin `mu_j`. Both are
**quadratic** in `t`; `mu` additionally carries the max/min disjunction whose branch is the
sign of `cross = det(e1, e2)` — the *same* cross that `convex_embed.jl` constrains, since
`det(X_v − X_prev, X_next − X_v) = det(e1, e2)`. Inside the convexity barrier every corner is
convex, so the reflex branch cannot occur while feasible.

The margin list is flattened into smooth entries: one `q_e` per split edge, and per corner
incidence one entry `−g_j` at a convex corner (`j` the argmax of `−g1, −g2`) or two entries
`−g1, −g2` at a reflex one. The objective is the log-sum-exp softmin

```
M_kappa(t) = -kappa * log sum_i exp(-c_i(t)/kappa)   <=   min_i c_i(t)
```

which is a **lower bound on `m` for every `t`**, and **tight at the point where the branch
modes were read** — so maximising it maximises a certified lower bound, and the reported `m`
is always the exact minimum, never the surrogate. Two unit tests check exactly those two
claims (`Kirigami/test/test_range_embed.jl`), and two more finite-difference the analytic
gradient of the whole stage-A objective — softmin plus both log barriers, split term and
corner term — against a central difference on two tilings and two barrier weights.

Stage A minimises `−M_kappa + w·(convexity + split log barriers)` with `w` decreased
geometrically over 6 continuation stages and the branch modes refreshed between them; once
feasibility is reached the continuation never leaves the feasible set, and the point kept is
the one with the largest **exact** `m`. Stage B, from a point with `m > 0`, pushes the exact
`Θ_max` with `range_opt.jl`'s softmin-of-first-contact objective (the T6 analytic
gradients, active set refreshed) under a **shrinking trust region** — caps of `0.25, 0.10,
0.04` median edges on `‖X − X_prev‖_∞` — accepting a step only if the exact convexity, the
exact split signs, a margin floor (`m ≥ ½ m₀`) and an actual `Θ_max` improvement all hold.
The margin constraint in stage B is enforced by exact rejection plus the trust region, not by
a barrier inside `range_opt`'s objective; that is weaker than a barrier and is stated as such.

### How it was run

Three arms on every design, each measured under the same certificate, the best recorded
**with its provenance** so the K9c contribution is never confused with a best-of-three count:

* **arm `k9`** — `convex_embed` with K9's settings (`δ = δ′ = 1e−3·med²`, 3 Gaussian starts
  + `t = 0`, 6 barrier stages), proximity objective.
* **arm `k9b`** — `convex_embed` with `δ = δ′ = 1e−2·med²` (K9b's most successful margin),
  8 Gaussian starts, 10 barrier stages, warm started from arm `k9`. K9b's full four-value
  sweep does not fit the wall budget; 27 of K9b's 28 positives won at `1e−3` or `1e−2`, and
  those two values are what the two proximity arms carry. This arm is therefore **K9b-like,
  not K9b**, and is labelled so everywhere.
* **arm `k9c`** — stage A from each of arm `k9`'s point, arm `k9b`'s point and `t = 0`
  (6 continuation stages × 120 L-BFGS iterations), then stage B.

Only the **winning** point of each design is refereed, at both shrinks `1e−12` and `1e−9`;
candidate comparison uses the exact scan alone, which is where the cost is not. 12 disjoint
shards over the graph index, so the designs completed are spread across the whole population
rather than a prefix.

**Coverage.** The machine was shared with eleven other concurrent experiments (load average
110–145 throughout), so the run overran its wall budget and was first reported at an interim
182 of 400 designs. It has since **finished all 400**; the 12 shards were re-merged with
`kill_k9c --aggregate --nshards 12` and `results/kill/k9c/{k9c.csv, summary.txt}` now hold the
full population. Every number below is against those 400 rows, and the interim snapshot is
preserved as a labelled subsection.

### Measured — `results/kill/k9c/{k9c.csv, summary.txt}` (400 of 400 designs)

| quantity | arm `k9` (proximity) | arm `k9b` (proximity, `δ=1e−2`) | **arm `k9c` (range-max)** | best of three |
|---|--:|--:|--:|--:|
| convexity + split feasible | 30 (7.5 %) | 24 (6.0 %) | **258 (64.5 %)** | 212 |
| ... of which 0⁺ margin `m > 0` | 23 | 12 | **239** | 210 |
| **exact `Θ_max > 0`** | 36 (9.0 %) | 27 (6.8 %) | **312 (78.0 %)** | **312** |
| refereed at shrink `1e−12` | — | — | — | **311** |
| refereed at shrink `1e−9` | — | — | — | **312** |
| certified (`ε_max > 0`) | 36 | 27 | **310** | 310 |
| ... of which `ε_max ≥ 0.1` rad | 33 | 14 | **192** | **201** |
| `ε_max` median / q90 / max (rad) | 0.253 / 0.540 / 2.00 | 0.115 / 0.470 / 1.83 | **0.162 / 0.958 / 3.14** | 0.172 / 0.958 / π |
| 0⁺ margin `m` over feasible, median (med²) | 0.050 | 0.004 | **0.529** | 0.503 |

The K9c arm's endpoints are those of a two-stage optimiser (stage A continuation, stage B
trust-region ascent) and are optimiser-path dependent: which of the 400 designs cross
`Θ_max = 0`, how far each opens and where the winning point comes from all move between
floating-point environments (an earlier pass of the same driver on the same population gave
307 positives, 306 certified, 214 at `ε_max ≥ 0.1` rad and a median `ε_max` of 0.278 rad),
while the two proximity arms, the `t = 0` columns and the referee do not. The rate itself is
stable to about ±2 % of the population; the numbers quoted here are those of the canonical
`k9c.csv`.

Per orientation, best-of-three: `σ_mc` **156 / 200 (78.0 %)** with **80** at
`ε_max ≥ 0.1` rad, `σ_def` **156 / 200 (78.0 %)** with **121** — the deployment asymmetry K9
and K9b measured has closed, but `σ_def`'s positives open wider. Per family,
Voronoi **126 / 134**, Delaunay **109 / 134**, quad-random **77 / 132**. Taking the better
`σ` per graph, **187 of the 200 graphs** carry a deployable design and **147** carry one with
`ε_max ≥ 0.1` rad. `‖X − X_ini‖` per vertex is median **0.927** median edges against K9's
0.38 and K9b's 0.21: the range-maximising point is farther from the input embedding, and that
is the price the objective is buying range with.

**Yield is not scale-free.** `results/final/figures/fig_yield.png` bins the population by
`|F|` with Wilson 95 % intervals: Delaunay falls from about 0.94 to 0.65 across the bins and
quad-random from 0.88 to 0.44. That trend must be reported with the headline count.

**Provenance of the winning point** (400 designs): arm `k9` wins on **102**, arm `k9b` on
**2**, a K9c point on **296** (`k9c/x0` 154, `k9c/x0+B` 127, `k9c/k9+B` 9, `k9c/k9b+B` 4,
`k9c/k9` 1, `k9c/k9b` 1). Most of arm `k9`'s wins (88 of 102) are designs where all three
arms return `Θ_max = 0` and the tie-break keeps the first arm. Within K9c, the `t = 0` start
supplies the winner on **281** and the two warm starts on 15, so **the range objective does
not need K9's feasible point to start from**; it reaches better points from the raw Eq. (6)
projection.

**Stage B is where much of the range comes from.** It supplied the winning point on **140 of
the 312 positives**. Stage A alone gets the design into the `m > 0` region — 239 of 400
designs, against 23 for the K9 proximity arm — and stage B then converts that margin into
angular range. Neither stage is sufficient alone.

### The referee

Soundness (`ε_max ≤ bisection at 1e−9`) holds **400 / 400**, and the two shrinks agree on
every positive but one: 311 at `1e−12` and 312 at `1e−9`, the exception being
`delaunay 136 σ_def` with `Θ_max = 1.9e−7` rad, below the bisection's own grid resolution
(`1e−12` reports 0, `1e−9` reports 1.4e−7); the same grid effect puts `delaunay 136 σ_mc` at
2.2e−7 under `1e−9` with the exact scan at 0. Neither is counted as a positive by the
referee. The F34 hinge-vertex artefact that cost K9 four designs does not appear here.

One disagreement is worth naming and runs the *other* way. On `delaunay 103 σ_def` the exact
scan reports `Θ_max = 0.2986` and the bisection reports `0.3294`; that is the source of the
`worst |exact − bisection(1e−9)| = 3.09e−02` line in the summary, and on `σ_mc` the same
statistic is `4.37e−05`. The bisection walks a grid of `π/4000 ≈ 7.9e−4` rad and can step
over a collision interval narrower than that; the exact scan cannot, because it tests the
overlap status at the closed-form contact angles themselves. The exact scan is therefore the
**smaller and the conservative** number here, the certificate is below both, and nothing in
the count depends on which of the two is preferred.

### What still binds on the designs that fail

Of the 88 designs with `Θ_max = 0`: `split-inward` on **58**, `inverted` on **24**,
`vertex-edge` on **6**. The `inverted` count is new — K9 and K9b never returned an inverted face, because their
barrier point was always feasible. Here the returned point is the least-infeasible one when
no start reaches feasibility, and the driver records it honestly (`best_feas = 0`) rather
than discarding it.

### K9b's nine lost designs

Confirmed on the completed run (`results/kill/k9c/summary.txt`, `K9 / K9b cross-check`):
K9 has 36 positives, K9b 28, and the 9 K9 designs K9b lost are recovered by the K9c arm
alone, 9 / 9, and again by the best-of-three row. All **9** are: `delaunay 10, 37, 43, 64, 112, 163, 190 (σ_def)` and
`delaunay 40, 103 (σ_mc)`. This is the direct test of K9b's diagnosis, and it comes out the
way the diagnosis predicted — those designs' feasible sets always contained a deployable
point, and the proximity objective was simply not looking for it.

### Verdict: PASS — the objective, not the feasible set, was the obstruction

| target | achieved |
|---|---|
| ≥ 40 / 400 refereed exact `Θ_max > 0` | **312 / 400** |
| as many as possible with `ε_max ≥ 0.1` rad | **192** (**201** best-of-3) |
| recover K9b's 9 lost designs | **9 / 9** |

The constructive claim changes in kind. K9 and K9b supported "a sparse subset of random
graphs deploys from a constrained embedding". K9c supports **"most random graphs deploy from
a constrained embedding, once the embedding is chosen to maximise the 0⁺ margin rather than
to stay near the input"** — 78.0 % of the 400 designs, against 0 for Eq. (6) alone
(K1a), `σ_mc` and `σ_def` (K5), convexity alone (K9a), K8a and the authors' native `prevent`,
and 2 / 400 for K6's primary repair and for B4 (17 / 400 counting K6's split-only relaxation,
which succeeds only by leaving the projection's neighbourhood), all on this exact population. The two proximity arms measured
here on the same 400 designs land at 9.0 % and 6.8 %, reproducing K9's published 9.0 % and
K9b's 7.0 %, which is the control that says the jump is the objective and not a change of
population or referee.

**Figures.** The full, non-curated figure set is `results/final/figures/` (see its README):
`fig_gallery_full.png` renders every deployable design, `fig_yield.png` the yield-versus-`|F|`
curves with Wilson intervals, `fig_feasible_vs_deployable.png` the margin-versus-range
scatter, `summary_table.md` the headline table. The earlier 20-graph gallery is
`results/kill/k9c/k9c_gallery.png` (`F` from 101 to 340, `ε_max` from 1.40 to π).
`results/kill/k9c/k9c_hist.png` — the `ε_max` distribution and the per-arm counts against
the baselines.

**Contradictions with STATE.md: one, and it is F36's amendment.** F36 records the K9/K9b
finding as "the ~9 % yield is basin-selected, not budget-limited … the selecting objective
must be the range itself". K9c confirms the *diagnosis* exactly — all 9 lost designs come
back — but the amendment's implicit ceiling does not survive: it treats the constrained slice
as sparsely deployable, and with the range objective it is deployable on **312 of the 400
designs**. F36 should be amended a second time to say that the constrained slice is
broadly deployable and that K9's 9 % was a property of the proximity objective alone. The
"first non-zero constructive result" framing in D9 and F36 is superseded, not contradicted.

Artifacts: `results/kill/k9c/{k9c.csv, summary.txt, shard_*.csv, shard_*.txt, gallery/,
k9c_gallery.png, k9c_hist.png}`. Driver `Kirigami/apps/kill_k9c.jl`, method
`Kirigami/src/method/range_embed.jl`, unit tests `Kirigami/test/test_range_embed.jl`, plots
`Kirigami/scripts/plot_k9c.jl`.

---

## Native200 — the authors' full pipeline on the K9 population (574/600 dispatched)


*Native baseline.* Every number from the authors' pipeline (Native200, the K2b ladders, the `native` cells of the regime study) comes from their unchanged binary (`baseline/native/`, `tuttekiri_cli`) and is reused as archived; only our side of each comparison is rerun.
**Hypothesis (fair baseline).** Every prior comparison against the authors' code was capped:
K2b ran their `prevent` on 30 live designs; K6's native subsection on 8 of 400 that reached a
valid flat state, refereed 0/80. Neither is the authors' *full* pipeline — their own
`initialized_two_face_coloring`, their Eq. (6), their Eq. (9), their FK collision test — run
**unconditionally**, with no pre-filter, on the same 200-graph population K5/K6/K9/K9b/K9c
use. Native200 is that run.

**How run.** `Kirigami/apps/kill_native200.jl`, three variants per graph: *native* (their own
coloring, via the CLI's `color` subcommand, fed to `prevent`), *sigma_mc* (our max-cut `σ`
baked into the input JSON), *sigma_def* (K5's defect-minimising `σ`). Each (graph, variant)
cell is `baseline/native/build/tuttekiri_cli prevent … --collisions --dump …`, wrapped in
`perl -e 'alarm shift; exec @ARGV' 600` (no GNU `timeout` on this machine), so a cell is
`completed`, `timed_out` (SIGALRM), or `crashed` (any other nonzero exit — `Segmentation
fault: 11` in every case observed). The optimised embedding, when one comes back, is
re-evaluated with our own instruments: the exact T4.2″ scan (`exact_theta_max_overlap`),
the bisection referee, and the `ε = 0.3` certificate (`POS ∧ NOOVERLAP ∧ NOROOT`) — the
authors' `native_theta_collisions` is recorded only as a diagnostic column (F24: their
`merge_close_verts` fuses hinge duplicates before the FK test runs). Sharded 12-way
(`gidx % 12`), run concurrently.

**Coverage, and why it is partial.** The 600 s pass was stopped by the user at 573/600
(id, variant) cells (`run_all.log`: the `xargs -P12` pipeline received `SIGTERM 15`, not a
crash): 172 completed, 362 timed out, 39 crashed. The 39 crashed cells were then rerun
against the F42-patched binary at a 3 600 s cap (`crashfix3600/`, 26 rows; F42 traced the
`Segmentation fault: 11` to their `merge_close_verts` being called on a deployed
configuration that had collapsed to zero faces, not to F24's `get_holes()` —
`results/kill/native200/crashfix.patch`; on 3 completed cells the patched and original
binaries agree exactly) and 4 of the timed-out cells at 3 600 s (`rerun3600/`).
`Kirigami/apps/native200_merge.jl` merges the three passes into
`results/kill/native200/native200_final.csv`, one row per cell with the 600 s status of a
superseded cell kept in `first_status`; **every number below is from that file.** On the
600-cell grid: **574 dispatched** (26 never dispatched), **195 completed, 364 timed out,
15 crashed.** Timeouts dominate (63 %) and are concentrated where `|F|` is largest — median
`|F|` is 202 for completed cells, 482 for timed-out ones, 431 for crashed ones, against
**332.5 for the full 200-graph population** (`results/kill/k9/k9.csv`) and **≤ 97 for the
2026 paper's own worked examples**. `quad_random` is worst hit: 150/198 timed out, 9/198
crashed, only 29/198 completed (22 of the 195 completions are crash-fixed cells, 16 of them
`quad_random`).

| | `native` | `sigma_mc` | `sigma_def` | total |
|---|--:|--:|--:|--:|
| completed | 58 | 64 | 73 | **195** |
| timed_out | 130 | 125 | 109 | **364** |
| crashed | 1 | 4 | 10 | **15** |
| never dispatched | 11 | 7 | 8 | 26 |
| dispatched | 189 | 193 | 192 | 574 |

| family | completed | timed_out | crashed | never dispatched |
|---|--:|--:|--:|--:|
| delaunay | 120 | 67 | 5 | 9 |
| voronoi | 46 | 147 | 1 | 7 |
| quad_random | 29 | 150 | 9 | 10 |

**Measured, on the 195 completed cells.** Their own FK collision test reports `Θ > 0` on
**27 / 195** completed (27 / 574 of all dispatched cells); our exact scan on the *same* dumped
embedding gives `Θ_max > 0` on **0 / 195**, the bisection referee **0 / 195**, and the
`ε = 0.3` certificate holds on **0 / 195**. By variant the 27 are 25 `sigma_def` and 2
`native`; by family `quad_random` 16, `voronoi` 8, `delaunay` 3. (The 600 s pass alone had
5 / 172 by their test and 0 / 172 by ours; 16 of the 22 additional positives of their test
come from crash-fixed `quad_random` cells.)

**The false-negative note (K1c mechanism).** All 27 designs where the authors' test says
`Θ > 0` but our exact scan on their own output says `0` fail `POS` at `ε = 0.3` — the exact
scan finds an immediate `0⁺` first-contact collision that the authors' `merge_close_verts`
fuses away before their FK test ever sees it, the same mechanism F24/K1c already documented
(there, the corrected re-closure rate on their `prevent` was 9.47 % against their own
reported 65.09 %). The full 27-row table (id, family, variant, `F`, `dim_null`, their `Θ`,
our exact and bisection values, all `0.0000`) is in `NATIVE200_FINAL.md` §3; the five from
the 600 s pass are:

| id | kind | variant | `F` | `dim_null` | native `Θ` | our exact | our bisect |
|---|---|---|--:|--:|--:|--:|--:|
| 160 | delaunay | native | 292 | 59 | 0.0880 | 0.0000 | 0.0000 |
| 36 | voronoi | sigma_def | 131 | 235 | 0.6726 | 0.0000 | 0.0000 |
| 102 | voronoi | sigma_def | 336 | 600 | 0.3508 | 0.0000 | 0.0000 |
| 168 | voronoi | sigma_def | 138 | 238 | 1.0248 | 0.0000 | 0.0000 |
| 156 | voronoi | sigma_def | 308 | 566 | 0.5704 | 0.0000 | 0.0000 |

These 27 are the K1c false-negative rate, now measured unconditionally instead of on a
hand-picked 169.

**Timing** (`secs` of `native200_final.csv`). Completed cells: median 13.5 s, max 570.7 s
(some finish just under the 600 s cap). Timed-out cells cluster at the alarm (median 600.2 s;
the tail to 3 711 s is the four 3 600 s reruns, the rest scheduling jitter under 12-way
contention). Crashed cells (the 15 that still crash under the patched binary): median
264 s, min 23.5 s.

**Verdict: baseline is 0 / 195 completed by exact scan, 27 / 195 by their own test.** The
authors' full pipeline — their coloring, their Eq. (6), their Eq. (9), their FK collision
test — run unconditionally on the K9 population, does not escape the `0⁺` jam any better than
the four repair strategies F30/K6 already killed, or the capped 8-design comparison in
K6's native subsection (0/8). The 27 apparent positives are false negatives of their own
collision predicate, not evidence against F30. **Caveat:** this run is not a clean 0 % — it
is 0 % *of what finished*, on a population whose median `|F|` (332.5, and worse for the
`quad_random` third that mostly never finished) is 3–4× the paper's own worked examples;
the authors' pipeline's practical operating range on this hardware is smaller graphs than
this population supplies, which is itself informative about the published method's scaling,
independent of the `0⁺` obstruction.

Artifacts: `results/kill/native200/{native200_final.csv, NATIVE200_FINAL.md, native200.csv,
summary.txt, shards/, crashfix3600/, rerun3600/, logs/, run_all.log, run_shard_{}.log}`.
Driver `Kirigami/apps/kill_native200.jl`, merge `Kirigami/apps/native200_merge.jl`. K9's
numbers for the same 200 graphs are `results/kill/k9/k9.csv`; nothing there was recomputed.

```
julia --project=Kirigami Kirigami/apps/kill_native200.jl --n 200 --maxf 800 --timeout 600 --shard {} --nshards 12 \
  --out results/kill/native200/shards --work /tmp/kiri_native200/shard_{} \
  --logdir results/kill/native200/logs --sigma results/kill/k5/sigma   # needs baseline/native
```
