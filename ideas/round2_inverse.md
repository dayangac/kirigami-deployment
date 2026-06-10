# Round 2 ideation — persona: inverse design / 2025-paper-centric geometry processing

Author: ideator2-inverse. Written against `specs/ideator_round2.md`.
Everything below is local; no external sources were consulted beyond the files named.

**What I read.** `specs/common_preamble.md`, `specs/ideator.md`, `specs/ideator_round2.md`,
`STATE.md` (F14–F30, D5–D9, dead ends), `notes/paper_2025.md` (§10.4–§10.7, §12–§15 in full;
§1–§9 skimmed for Def. 4.1 / Prop. 4.1 / Obs. 3.1 / θ_max), `notes/repo_2025.md` (a heading
skeleton only — §1–§7 are still empty as of writing, so I used the repo files directly),
`derivations/core.md` (§0, T1, T3, T4.4, T4.5, T7), `results/kill/KILL_REPORT.md` (K5, K2b,
K1a, K2a headline numbers), `code/src/method/periodic_jacobian.hpp` (the derived comment block),
`code/src/method/contact.hpp`, `code/apps/kill_k7.cpp`, and the **measured K7 output**
`results/kill/k7/k7_main_*.csv` + `k7_c4_bounded.csv`, plus the 16 pattern JSONs and 4 target
meshes at `baseline/kirigami_tessellations/code/data/`.
I did **not** read `notes/field_rigidity.md`, `papers/related/isogami.txt`, or the round-1 persona
files in full; novelty claims below are checked against `notes/field_kirigami.md`'s R-table (R1–R24)
and `notes/screen_bundle.md`'s verdicts S1–S7, and I say so per idea.

---

## The measured facts this round is built on

These are numbers already in the repo, not new claims. Everything below cites them by tag.

**M1 (K7, `results/kill/k7/k7_main_*.csv`, `k7_c4_bounded.csv`).** The second gapless angle
`θ_c = 2·atan2(tr K, 1 − det K)` matches forward kinematics to `≤ 1.4e−14` on every pattern in the
K7 population, and takes these values on the regular families:

| family | `n_split` | `dim_null` (quotient) | `dim K` | `θ_c` | `θ_c` as a fraction |
|---|--:|--:|--:|--:|---|
| triangles (2×2, 3×2, 3×3) | 0 | 0 | 0 | 4.18879 | `4π/3` |
| squares 2×2 | 0 | 0 | 0 | 3.14159 | `π` |
| kagome (2×2, 3×2, 3×3) | 0 | 0 | 0 | 3.14159 | `π` |
| hexagons (2×2, 3×2, 3×3) | 4, 6, 9 | 4, 6, 9 | 4 | 2.09440 | `2π/3` |
| snub square | 8–18 | 8–18 | 4 | 3.60–3.77 | — |
| 4.8.8 | 8–18 | 8–18 | 4 | 1.73–2.00 | — |
| 3.4.3.12 | 8–24 | 8–24 | 4 | 2.52–2.74 | — |

**M2 (T4.4 `[N]`, F18).** On split-free patterns `Θ_max = min(min_e β_e, π)` to `8.9e−16` on 8
patterns, `β_e = 2π − α_f − α_g`. For the three regular rotating patterns `min β` is
`4π/3` (triangles), `π` (squares), `2π/3` (hexagons) — **identical to `θ_c` in M1, exactly.**

**M3 (F15 + T7(vi)).** `dim 𝕏 = V_int − H = |E_split|` with the boundary fixed;
`H = |E_hinge| − |F| + c(Γ) = V_int − |E_split| + (c(Γ) − 1)` on a disk.

**M4 (F30 / K5 / K6).** On random planar graphs the 2026 shape space is certifiably empty
(0/400 with `Θ_max > 0`); authored tilings never fail (0/8). Split cuts are the only design
freedom **and** the only collision source.

**M5 (2025 §6, p.10, verbatim).** "a fully closed hemisphere shape can be achieved using a quad
pattern (see Fig. 1) and hexagonal pattern (Fig. 17), but **not with the triangular pattern**
(see Fig F.4, right)"; and "Our current investigation of the realizable shape space associated
with each pattern is largely **empirical** (Fig. 14); developing a mathematically rigorous
framework … remains an open problem."

**M6 (2025 §5, p.9).** "Patterns with hole-free configurations at maximum opening angle … are
more structurally constrained. Achieving high curvature in such cases requires introducing
singularities" (Fig. 15, full sphere, n-RoSy initialisation).

**M7 (2025 §5, p.9, fabrication).** "thin joints … are fragile and prone to breaking during
rotation; thicker hinges offer durability but **limit the achievable rotation angles**."

---

## P0 — shared infrastructure every kill test below needs (≈ 90 lines of C++, one hour, once)

`baseline/kirigami_tessellations/code/data/patterns/*.json` is a **polygon soup**: a JSON array of
`{points: [{x,y}...], color}` objects, one per tile of a unit patch, with no shared-vertex
information. Every kill test below needs the same importer:

```
kiri::Mesh import_kiri_json(path, weld_tol)   // weld coincident corners -> Mesh + topology
```

plus `detect_lattice` (already in `method/periodic_jacobian.hpp`) to recover `t1, t2`, and the
2025 complete cut (`σ` = the dual 2-colouring, which exists iff every interior valency is even —
2025 Remark 4.1; our `orientation.hpp` max-cut already returns it when the dual is bipartite).
This is **not** an idea; it is the prerequisite. Budget it once and every test below is minutes.
If the weld fails on some patterns, the honest fallback is our own generators (`squares`,
`triangles`, `kagome`, `hexagons`, `snub_square`, `trunc_square_488`, `t3_4_3_12` — already in
`kill_k7.cpp`) plus a hand-authored `hexagon_triangle`, and the tests below still run, with the
sample described as "our reconstruction of their families", not "their library".

---
# The ten ideas

---

## I1 — Closure margin: why a triangle tiling cannot close and a quad tiling can
**Type: characterization.**

**Claim.** For a uniformly deployable split-free tiling define the **closure margin**
`m(X) := min_e β_e(X) − θ_c(X)`, with `β_e = 2π − α_f − α_g` (2025 §4.2) and
`θ_c = 2·atan2(tr K, 1 − det K)` the second gapless angle (M1). A gapless deployed state exists
iff `m(X) ≥ 0`. On the three regular rotating families `m ≡ 0` **exactly** (M1 = M2), so closure
is a knife edge; the pattern families separate not by `m` at the regular point but by
`sup_{X ∈ 𝕏} m(X)` over the pattern's own shape space, which is `0` for triangles
(`dim 𝕏 = 0`, M1) and strictly positive for hexagon/quad families with `|E_split| > 0`.

**Gap targeted.** 2025 §6 p.10 (M5): the triangle/quad/hex hemisphere fact is stated as a
**fabrication observation** with evidence in the unavailable Fig. F.4, and "developing a
mathematically rigorous framework … remains an open problem". 2025 open question 18 in
`notes/paper_2025.md` §15 ("when is the max-angle deployed configuration hole-free?") is exactly
`m ≥ 0` and has no answer in either paper.

**Novelty vs field.** R1 (Grima & Evans, rotating squares/triangles) computes the closed states of
the three regular patterns individually but has no shape space and no margin. R5/R7
(Konaković 2016/2018) fix the triangular linkage a priori and never ask which patterns admit a
gapless state. R12 (Dang et al. 2021) proves rigid deployability for quads ⟺ all voids are
parallelograms at some deployed state — a **parallelogram** condition, not a **zero-area**
condition, and quads only. `notes/screen_bundle.md` S7: "hole count / area harmonic / second
closed angle NOT FOUND" in the literature search.

**Why it might be true / sketch.** Rigid faces ⇒ material area per cell is `θ`-independent
(T1.D), so total hole area per fundamental cell is exactly
`A(θ) = det(P₀)·(det J(θ) − 1)` with `J(θ) = cos(θ/2) I + sin(θ/2) K` (`periodic_jacobian.hpp`).
Expanding with `c² = (1+cos θ)/2`, `s² = (1−cos θ)/2`, `cs = (sin θ)/2`:

```
   det J(θ) − 1 = ½(det K − 1)(1 − cos θ) + ½ (tr K) sin θ ,
```
a first harmonic with `p + q = 0`. Its roots are `θ = 0` and `tan(θ/2) = tr K/(1 − det K)`, i.e.
`θ_c`. Because `A ≥ 0` on the physical range, `A` is a single positive bump on `(0, θ_c)`.
Closure is reachable iff no contact happens first, and on split-free patterns the first contact is
the hinge-sector bound (M2). Hence `m ≥ 0`. For a regular `n`-gon rotating pattern the face must
rotate by `θ_c/2 = 2π/n` to map to itself, giving the closed form **`θ_c = 4π/n`**:
`4π/3, π, 2π/3` for `n = 3, 4, 6` — reproducing M1 to `1e−15`, and simultaneously
`β = 2π − 2·(interior angle) = 4π/n` as well, which is why `m ≡ 0` there.

**Kill experiment.** Driver `kill_i1.cpp` on top of P0. (a) Sample `X ∈ 𝕏` on a free-boundary
`m×m` patch (`m = 4,6,8`) of each of `triangles`, `squares`, `kagome`, `hexagon_triangle`,
`snub_square`, `4.8.8`, `3.4.3.12`, 200 samples each: compute `θ_c` per cell from `K`, `min β`
exactly (`hinge_beta()`), and `m`. (b) Compute `sup m` by local ascent in `𝕏`.
**Kills the idea if:** `sup_{𝕏} m > 0` for the triangle family on any patch size, **or**
`sup_{𝕏} m ≤ 0` for both hexagon_triangle and a quad family. Either outcome breaks the claimed
separation. Runtime: seconds per sample, ≪ 30 min.

**If it survives, the demo.** A single scatter, one dot per (pattern, sample): `sup m` on the
`x`-axis, "gapless state reached in forward kinematics" on the `y`-axis, with the three families
of M5 labelled. Hero: side-by-side rendered deployments of the triangle patch (jams with holes
open) and the hexagon_triangle patch (closes), both certified by `exact_theta_max`. Fabricable
export: the two patches as SVG cut files at their respective `min(θ_c, Θ_max)`.

**Risk.** `m ≡ 0` on all three regular patterns means the whole separation rests on the shape
space, not on the margin. If `dim 𝕏 > 0` for triangles on a free-boundary patch (I3 predicts
`|V_∂|` DOF, which is **not** zero) and those boundary DOFs are enough to open a positive margin
in the interior, the criterion collapses and the honest answer becomes "triangles fail for a
fabrication reason (M7: 120° of face rotation vs 90° and 60°), not a geometric one". I regard this
as the single most likely outcome and it is why I1 is ranked below I3.

**Effort.** Derivation S (done above, modulo the `sup` claim). Code M (needs P0).

---

## I2 — The curvature budget of a pattern, and the `(1 − Ω/4π)^{-2}` threshold
**Type: characterization (with an exact closed form and an exactly-computable threshold).**

**Claim.** Under the homogenisation ansatz below, the areal expansion a pattern can deliver is
```
   λ_A(θ, K) = det J(θ) = ½(1 + det K) + ½(1 − det K) cos θ + ½ (tr K) sin θ ,
   B(pattern) := sup_{K ∈ 𝒦}  λ_A(θ†, K) / inf_{K ∈ 𝒦} λ_A(θ†, K)   over θ† ∈ (0, Θ_max],
```
with `𝒦 = {K₀ + Σ_j t_j M_j}` the **affine** achievable set already computed by
`AchievableSet` in `kill_k7.cpp`. A spherical cap of total Gaussian curvature `Ω` is conformally
realizable from a flat sheet only if `B ≥ (1 − Ω/4π)^{-2}`; for a closed hemisphere `Ω = 2π` and
the threshold is exactly **4**; for a full sphere `Ω = 4π` it is `+∞`, so no pattern can close a
sphere without singularities.

**Gap targeted.** 2025 §6 p.10 verbatim (M5): the realizable shape space per pattern is
"largely empirical (Fig. 14)"; and the same paragraph asks for "a continuous model to analyze the
behavior of the hinged kirigami structures. If such a model exists, one could gain a better
understanding of the space of shapes that can be achieved from a given pattern". This is that
model, in closed form. 2026 §5.1 only asserts conformality "empirically for all θ".

**Novelty vs field.** R5/R7 (Konaković 2016/2018) use bounded-conformal-factor maps for a
**single** fixed triangular linkage with an empirically measured maximum expansion; the budget is
an input constant there, not derived from the pattern. R20 (Zaman 2025) varies tiles spatially but
for actuation, with no metric bound. R22 (IsoGami) is isohedral and does not produce a per-pattern
metric budget. `notes/screen_bundle.md` S4 ("exact conformality") returned NOT FOUND, and S7
(hole-area harmonic) NOT FOUND.

**Why it might be true / sketch.** The cap threshold is exact and I derived it here. Stereographic
projection of a cap of polar half-angle `θ₀` from the far pole gives the unit disk with conformal
factor `λ(r) = 2/(1 + r²)`, `r ≤ r₀ = tan(θ₀/2)`, so
`λ_max/λ_min = 1 + r₀² = sec²(θ₀/2)`. The cap's total curvature is `Ω = 2π(1 − cos θ₀)`, hence
`1 − Ω/4π = (1 + cos θ₀)/2 = cos²(θ₀/2)`, giving
```
   λ_max/λ_min = (1 − Ω/4π)^{-1}          ⇒     areal ratio = (1 − Ω/4π)^{-2} .
```
Hemisphere `θ₀ = π/2`: ratio `4`. Quarter cap `θ₀ = π/3`: `(3/4)^{-2} = 16/9 ≈ 1.78`. Full sphere:
divergent — which is precisely M6, the one place 2025 needs singularities.
**Homogenisation ansatz (stated, not proved):** the deployed 3D surface is, per cell and up to a
rigid motion of ℝ³, the *planar* deployed cell at the common angle `θ`, so the induced metric of
the coarse surface is the planar cell metric `J(θ)ᵀJ(θ)` and the local area scaling is `det J(θ)`.
This is the same ansatz Konaković's auxetic pipeline uses implicitly; here it is explicit and it
has its own kill test (part (c) below).

**Kill experiment.** `kill_i2.cpp`. (a) For each of the 7 families in the K7 population plus
`hexagon_triangle`, compute `𝒦` (`achievable()`), then `B` by maximising and minimising
`λ_A(θ†, K)` over `t` at 64 values of `θ†` in `(0, Θ_max]` — `λ_A` is a quadratic in `t`, so each
extremum is a small eigenproblem or a bounded ascent from 32 starts. (b) Check the ordering
`B(triangles) < 4 ≤ B(quad family), B(hexagon_triangle)`. (c) **Ansatz check**: take the target
`hemisphere2.obj`, compute a discrete conformal map to the plane, histogram the conformal factor,
and verify its `(λ_max/λ_min)²` equals `4` to within the mesh's own boundary truncation.
**Kills the idea if:** `B(triangles) ≥ 4`, **or** `B(quad) < 4` and `B(hex) < 4` (then the
threshold does not separate the paper's own three cases), **or** (c) returns a ratio not within
±20 % of 4 (then the threshold is being computed for the wrong object). Runtime < 30 min;
(a) alone is minutes and is the decisive half.

**If it survives, the demo.** The budget table for all 16 of their patterns, with a horizontal
line at `4` and the paper's own quad/hex/triangle verdicts as ticks and crosses; a
"budget vs achievable cap angle" curve `θ₀(B) = 2 arccos(B^{-1/4})`; and a plot of budget against
`dim 𝒦`, testing the paper's own "fewer holes ⇒ harder" claim (Fig. 14 caption) quantitatively.
Hero: the deepest cap each of the 16 patterns can reach.

**Risk.** The ansatz. 2025's deployed state is genuinely 3D with per-hinge angles, and their
`E_reconfig` is a **soft** penalty (2025 §12 assumption 4), so the local cell need not be a rigid
copy of any planar deployed cell. If the per-cell 3D metric drifts from `cos(θ/2)I + sin(θ/2)K`
by more than the reported `E_r ≈ 3e−4`, the budget bounds nothing. Second risk: `𝒦` as coded is
the **fixed-lattice** achievable set, which is an inner approximation of the true spatially-varying
achievable set, so `B` as measured is a lower bound and a `B < 4` for quads would be inconclusive.
Both must be said out loud in any write-up.

**Effort.** Derivation M (the cap threshold is done; the ansatz needs a careful statement).
Code M.

---

## I3 — The design space of a complete-cut tiling is its boundary, and nothing else
**Type: theorem.**

**Claim.** For a disk-topology tiling with `c(Γ) = 1` and `rank(L) = H`, the Tutte auxetic shape
space with a **free** boundary has scalar dimension
```
   k  =  |V|  −  H  =  |V_∂|  +  |E_split|  +  1 − c(Γ)      ( = |V_∂| + |E_split| when c(Γ)=1 ).
```
In particular for a 2025-style **complete cut** (`E_split = ∅`, forced by Remark 4.1's even-valency
condition) the design space has exactly `|V_∂|` scalar degrees of freedom — **independent of the
pattern type, of the number of faces, and of the number of holes** — so it grows as `Θ(√n)` while
the number of shape-approximation residuals grows as `Θ(n)`.

**Gap targeted.** 2026 §4.4 states "the number of independent equations equals the number of
holes" and never counts the resulting dimension with a free boundary; F15 measured only the
**fixed**-boundary case `k = |E_split|`. 2025 never states the dimension of its own design space at
all, yet Fig. 12 (`n_f` = 213 → 7068, `E_r` **non-monotone**, "the paper does not comment on this",
`notes/paper_2025.md` §10.2) and Fig. 14 ("fewer holes to absorb local distortions") are both
empirical observations about exactly this quantity.

**Novelty vs field.** R12 (Dang 2021) and R17 (Dudte 2023) count DOFs for quad kirigami only and by
linkage mobility, not by a Tutte-type linear system with a free boundary. R2 (Mitschke 2013)
searches an archive. Nothing in the R-table gives a pattern-independent dimension formula.
`notes/screen_bundle.md` S5 records "per-face rank claim NOT FOUND".

**Why it might be true / sketch.** One line from T7(vi). `H = |E_hinge| − |F| + c(Γ)` and
`|E_hinge| = |E| − |E_∂| − |E_split|` with `|E_∂| = |V_∂|` on a disk. Then
```
   k = |V| − H = |V| − |E| + |V_∂| + |E_split| + |F| − c(Γ) = (|V| − |E| + |F|) + |V_∂| + |E_split| − c(Γ)
     = 1 + |V_∂| + |E_split| − c(Γ)                                        (Euler, disk: V − E + F = 1)
```
counting `F` as the bounded faces. Specialising to a fixed boundary removes `|V_∂|` unknowns and
`c(Γ) = 1` gives `k = |E_split|`, **which is exactly F15's measurement on 100/100 graphs** — so the
identity is already validated in its fixed-boundary specialisation. The free-boundary statement is
the one that has never been drawn, and its corollary (all 2025 patterns have `|V_∂|` DOF) is the
non-obvious part: it says the interior of a complete-cut tiling is *design-rigid*.

**Kill experiment.** `kill_i3.cpp`, no P0 needed. Generate `m×m` patches, `m = 2 … 20`, for all 7
K7 families plus 200 random even-valency tilings, with a **free** boundary; build `[L]` (no `B`
rows), take `dim null` by dense SVD for `n ≤ 700` columns (F29 forbids trusting SparseQR above
that) and report `|V_∂| + |E_split| + 1 − c(Γ)` alongside. **Kills the idea if:** the two disagree
on any single instance where `rank(L) = H` holds. Also fit `k` against `n` on a log-log plot and
kill the `Θ(√n)` corollary if the fitted exponent for a split-free family is not `0.5 ± 0.05`.
Runtime: minutes.

**If it survives, the demo.** A log-log plot of `k` vs `n_f` for the 7 families, all collapsing on
one line of slope `1/2` for split-free patterns and slope `1` once split cuts are allowed; overlaid
with 2025's own Fig. 12 numbers (`n_f`, `E_r`) to show the non-monotone `E_r` is a DOF-starvation
effect. Hero: two tilings with the same face count and a 3× difference in `|V_∂|`, one of which can
hit a target the other cannot.

**Risk.** It is a corollary of T7(vi) plus Euler, so a hostile reader calls it bookkeeping. The
defence has to be the corollary, not the formula: "the design space of every hinged kirigami tiling
is its boundary" is a statement about the method, and it says 2025's densest results
(`n_f = 7068`) have essentially the same design freedom as their sparsest. Secondary risk:
`rank(L) = H` can fail (T7(ii): `rank L = H − dim Z`), and on a free-boundary patch I have not
checked that `dim Z = 0`; the periodic case has `dim Z = 1` (F22). Third: the formula is **disk-only**. On the torus
quotient the rank is `H − 1` (F22) and the system carries the lattice-offset right-hand side, so the
count is a different one; I checked my formula against the K7 quotient columns and could not make it
reproduce `dim_null` for `squares_3x2` and `hexagons_2x2` simultaneously without knowing `c(Γ)` per
row, so I make **no** periodic claim here.

**Effort.** Derivation S (done). Code S.

---

## I4 — The maximum-area configuration sits at exactly half the second closed angle
**Type: theorem.**

**Claim.** For a periodic uniformly deployable pattern, the total hole area per cell
`A(θ) = det(P₀)(det J(θ) − 1)` attains its maximum at `θ_A = θ_c/2`, exactly, where `θ_c` is the
second gapless angle. Hence 2025's two maximally open configurations (Observation 3.1) are in a
fixed `1 : 2` ratio whenever the max-angle state is the gapless one, and the max-area state is
always reached strictly before it.

**Gap targeted.** 2025 Observation 3.1 (§3.2, p.4) says a pattern "typically" admits two distinct
maximally open configurations and gives no relation between them; `notes/paper_2025.md` §15 open
question 17 asks explicitly "when do they coincide? Can max-area come *after* max-angle …? The
paper's figures always show max. area before max. angle (Fig. 2, Fig. 5) but never argues that this
ordering is general." This answers that question in closed form.

**Novelty vs field.** R1 (Grima & Evans) knows the max-area state of rotating squares by symmetry,
not as a general identity. Nothing in R1–R24 states a relation between the two maximal states for
an arbitrary pattern. S7 in `notes/screen_bundle.md`: NOT FOUND.

**Why it might be true / sketch.** From I1's expansion, `A(θ)/det(P₀) = p + q cos θ + r sin θ` with
`p = −q = (det K − 1)/2` and `r = (tr K)/2`. A first harmonic with `p + q = 0` has its stationary
point at `θ_A = atan2(r, q)` and its non-zero root at `tan(θ_c/2) = r/(−q) = tr K/(1 − det K)`.
So `θ_A = arctan(r/q)`… and `θ_c = 2 arctan(r/(−q))`; with `q = (det K − 1)/2 < 0` for an expanding
pattern, `θ_A = atan2(r, q)` and `θ_c/2 = arctan(r/|q|)` coincide. Sanity check on rotating squares:
`K = I` (isotropic, `k = 1`) gives `θ_A = π/2`, `θ_c = π` — squares are maximally open at 45° of
face rotation and re-tile at 90°, both correct. Regular `n`-gon: `θ_A = 2π/n`, `θ_c = 4π/n`.

**Kill experiment.** `kill_i4.cpp` (or three lines added to `kill_k7.cpp`, which already computes
`θ_c` and the area harmonic): on the whole K7 population plus 200 random torus patterns, find
`argmax A(θ)` by golden section on the forward-kinematic area (not on the harmonic — that would be
circular) and compare to `θ_c/2`. **Kills the idea if** `|θ_A − θ_c/2| > 1e−9` on any pattern where
both are inside `(0, π]`. Runtime: minutes, and it reuses `k7_c4_bounded.csv`'s validated area
measurement path.

**If it survives, the demo.** One figure: `A(θ)` for the 16 patterns of their library normalised to
`θ/θ_c`, all collapsing onto a single `sin`-bump with its peak at `1/2` — the "universal deployment
curve". Hero: their own Fig. 2 columns 3 and 5 annotated with the measured `θ_A` and `θ_c` and the
2:1 ratio.

**Risk.** It is two lines of trigonometry once the harmonic form is known, and the harmonic form is
already ours (B20/K7 C4). A reviewer will call it a corollary. Its value is entirely that it settles
a question the paper asks out loud, and it must be presented that way, packaged with I1/I2, never
as a standalone contribution. Second risk: when the max-angle state is *not* gapless (`θ_c > Θ_max`,
e.g. triangles), the 2:1 ratio is vacuous, and that is the majority of interesting patterns.

**Effort.** Derivation S (done). Code S.

---
## I5 — Automatic tiling selection from a target's conformal-factor spectrum
**Type: algorithm.**

**Claim.** Given a target surface `S`, compute a discrete conformal map to the plane, extract the
conformal factor field `λ`, and form the **demand pair** `ρ(S) = (λ_max/λ_min)²` (areal) and
`κ(S) = ` the 90th-percentile anisotropy of the principal-stretch ratio. Then the pattern selected
as `argmin { B(P) : B(P) ≥ ρ(S), anisotropy budget of P ≥ κ(S) }` — the *cheapest* pattern that
clears the demand — reproduces 2025's own Fig. 13 assignments and their triangle/quad/hex verdict,
without running their optimizer.

**Gap targeted.** 2025 §6 p.10 verbatim: "**automatically identifying optimal tiling types from
geometric features, such as the curvature profile, is an interesting direction for future
research**". 2025 §12 assumption 8: "the tiling type is provided as input rather than derived from
the target geometry".

**Novelty vs field.** R7 (Konaković-Luković 2018) selects nothing — one linkage, and the conformal
factor is *clamped* to what the linkage can do. R22 (IsoGami) chooses among isohedral tilings by
optimization, not by a computed budget. R20 varies tiles, not tiling type. No row in R1–R24
computes a per-pattern realizability certificate to select against.

**Why it might work / sketch.** I2 makes `B(P)` a computable scalar per pattern, and `ρ(S)` is a
computable scalar per target, so selection collapses to a table lookup with two inequalities. The
"cheapest that clears" rule (rather than "largest budget") is the non-obvious half: M6 says
hole-free / structurally constrained patterns are *desirable* (least material waste, most
constrained, best fabricability) so one wants the tightest fit, not the loosest.

**Kill experiment.** `kill_i5.cpp`. On the 4 target meshes shipped with their repo
(`hemisphere2.obj`, `pringles2.obj`, `bumps_plane.obj`, `lilium.obj`) compute `ρ` and `κ`; on the
16 pattern JSONs compute `B` and the anisotropy budget (I2, I8). Score the selector against the
only ground truth available: their Fig. 13's 6×7 grid as transcribed in `notes/paper_2025.md`
§10.1, plus the three verdicts of M5. **Kills the idea if** the selector puts `triangles` above
either the quad or the hexagon family for `hemisphere2.obj`, or if `ρ(hemisphere2)` is not within
±20 % of `4`. Runtime: minutes after P0 and I2.

**If it survives, the demo.** A 16 × 4 heat map of "demand cleared / not cleared", with their own
successful (pattern, shape) pairs from Fig. 13 marked, and the count of Fig. 13 cells the selector
would have rejected. Hero: a target whose demand *no* pattern in their library clears, with the
minimal singularity count from I7 printed next to it.

**Risk.** The ground truth is thin. Fig. 13 shows 6 patterns × 7 shapes all *succeeding*, so it
provides no negatives; the only negative in the whole paper is the triangle/hemisphere one and its
figure is in the unavailable supplement (§14). A selector validated on one negative is not
validated. Second risk: `ρ` depends on how the boundary of the cut-out flat domain is chosen, which
is a free parameter of the discrete conformal map, and could be tuned to fit.

**Effort.** Derivation S. Code M (needs a discrete conformal map — CETM or a simple Yamabe flow, in
C++; this is the real cost and it may not be worth paying).

---

## I6 — Certified-deployable inverse design: hard constraints where 2025 uses penalties
**Type: algorithm.**

**Claim.** Replacing 2025's soft `E_reconfig` (Eq. 3) and `E_planar` (Eq. 4) by (i) an exact
parametrization `X(t) = X₀ + Φ t` of the shape space and (ii) the exact collision certificate
`POS ∧ NOOVERLAP(ε/2) ∧ NOROOT` of T5.2b′ as a hard filter, yields designs with reconfigurability
error `E_r` at machine precision (vs their reported `2.7e−4 … 7.8e−4`, `notes/paper_2025.md` §10.2)
and a **certified** collision-free deployment range, at comparable shape-approximation error.

**Gap targeted.** 2025 §12 assumption 4 (rigidity is relaxed into soft penalties, hence the nonzero
reported `E_r`, `E_p`) and assumption 7 ("Non-adjacent (global) self-collision during deployment is
not modeled"). 2026 §6 lists the absence of a validity certificate as open; F17 quantified it
(Eq. (6)'s projection self-intersects on 57/57 random graphs).

**Novelty vs field.** Both Segall papers *simulate* deployment and check it visually; no row of
R1–R24 carries a certificate. R18 (Liu 2024) has a closed-form collision test for one dihedral
Escher family only. R12/R17's linear-algebraic quad formulations have no collision model at all.

**Why it might work / sketch.** `𝕏` is affine (T1, T7) so hard deployability costs nothing:
optimise over `t` directly, and every iterate is exactly deployable. Every contact predicate is
`p + q cos θ + r sin θ` with `p, q, r` quadratic in `t` (T3.2), so `Θ_max(t)` has a closed-form
gradient wherever the active contact is unique (T6.2), giving a differentiable barrier
`−log(Θ_max(t) − ε)` that is *exact* rather than the heuristic Eq. (9) of 2026 whose false-negative
rate we measured at **9.47 %** (K1c/F28).

**Kill experiment.** `kill_i6.cpp`. On a 2025 pattern patch (a quad family, ~400 faces) with a free
boundary, run: (a) our constrained optimizer against `hemisphere2.obj` with the certificate as a
hard filter; (b) the same shape term with 2025's soft formulation as we can reconstruct it.
Compare `E_r`, `E_p`, one-sided Hausdorff distance to the target, and certified `Θ_max`.
**Kills the idea if** the hard-constrained result's shape error exceeds the soft result's by more
than 2× on ≥ 3 of 5 targets — meaning the `|V_∂|` degrees of freedom of I3 are simply not enough,
and 2025's softness is load-bearing rather than sloppy. Runtime: the optimizer is the cost; keep it
to 5 targets and a 30-minute cap, and report honestly if it does not converge.

**If it survives, the demo.** A ≥ 20-design table: `E_r` (ours `~1e−15`, theirs `~1e−4`), certified
`Θ_max`, and Hausdorff error, with the 2026 Eq. (9) baseline as a third column. Hero: a fabricated
SVG whose deployment is certified rather than simulated. Fabricable export: `code/src/export`
already emits SVG/STL/3MF.

**Risk.** K2b already failed once at exactly this game: against the authors' **native** code our
range optimisation had median gain 0 (7 better / 15 worse / 16 ties, F27). The lesson stands —
their heuristics are strong. The defensible claim here is *certification*, not *margin*, and the
write-up must not slide from one into the other. Second risk: I3 says the design space is only
`|V_∂|`-dimensional, which makes it quite likely that the hard-constrained optimum is much worse in
shape error, i.e. that this idea's kill test fails for the reason I3 predicts.

**Effort.** Derivation S (all pieces exist). Code L.

---

## I7 — A lower bound on the number of singularities, from budget and Gauss–Bonnet
**Type: theorem.**

**Claim.** If a pattern's areal budget is `B` and the target region carries total Gaussian curvature
`Ω`, then any cone-metric decomposition realizing it needs at least
```
   N  ≥  Ω / ( 4π ( 1 − B^{-1/2} ) )
```
cone points, because each cone can absorb at most `4π(1 − B^{-1/2})` of curvature before the
surrounding cap's demand `(1 − Ω_local/4π)^{-2}` exceeds `B`.

**Gap targeted.** 2025 §6 p.10, last sentence: "one could gain a better understanding of the space
of shapes … **and where singularities could be placed to extend it**". Fig. 15 introduces
singularities for a full sphere with no theory for how many or where (2025 §15, open question 16).

**Novelty vs field.** Meekes & Vaxman 2021 (cited by 2025 for the n-RoSy initialisation) places
singularities by field design, with no metric budget. R5/R7 never use cone points. No R-row bounds
a singularity count by a linkage's expansion capacity.

**Why it might be true / sketch.** Invert I2's cap threshold: a cap of budget-limited curvature has
`Ω_max = 4π(1 − B^{-1/2})`. Decompose the target into `N` budget-feasible caps whose curvatures sum
to `Ω` by Gauss–Bonnet; the count follows. Sanity check: for `B = 2` (rotating squares, area
doubling at full opening) `Ω_max = 4π(1 − 0.7071) = 3.68`, and a full sphere `Ω = 4π = 12.57`
needs `N ≥ 3.4`, i.e. **4 cones**.

**Kill experiment.** `kill_i7.cpp`. Take a family of spherical caps of increasing `θ₀` and, for one
pattern with a measured `B`, run I6's optimizer with `N = 0, 1, 2, …` prescribed cone points and
record the first `N` at which the shape error drops below 2025's own threshold (`E_r < 1e−2`,
`E_p < 0.1°`, their Fig. 13 criterion). **Kills the idea if** the measured `N` is below the bound
for any cap — a lower bound that is violated is simply false. Runtime: this one does **not** fit in
30 minutes honestly, because it needs I6's optimizer per cone count; the cheap surrogate is to test
the bound against the *metric* problem only (does a cone metric with `N` cones and factor ratio
`≤ B` exist?), which is a linear feasibility problem and does run in minutes.

**If it survives, the demo.** A staircase plot: required cones vs cap angle, one curve per pattern
budget, with 2025's Fig. 15 sphere marked. Hero: the minimal-singularity sphere for their quad
pattern, next to their n-RoSy result.

**Risk.** The decomposition argument is hand-wavy as written — caps do not tile a surface and
curvature is not additively partitionable into disjoint caps without overlap loss, so the constant
is probably wrong even if the scaling is right. This needs a real proof (a covering/Bishop-type
argument) before it can be called a theorem, and I do not have one. Also depends entirely on I2.

**Effort.** Derivation L (the proof does not exist yet). Code M.

---

## I8 — Anisotropy budget and a pattern-orientation field for saddle targets
**Type: characterization + algorithm.**

**Claim.** The polar stretch of `J(θ) = cos(θ/2) I + sin(θ/2) K` traces, as `θ` and `K ∈ 𝒦` vary, a
2-parameter region in principal-stretch space `(λ₁, λ₂)`; a pattern is *conformal* iff `𝒦` lies in
`span{I, J_rot}` (proved in `periodic_jacobian.hpp`, verified by K7 C2 to `c2_conf_res ≈ 1e−17` on
hexagons and `≈ 1e−17` on snub square). For anisotropic targets the pattern must additionally be
*rotated* to align its principal stretch direction with the target's principal curvature direction,
turning tiling selection into a field-design problem with a computable per-point feasibility test.

**Gap targeted.** 2025 §6 p.10: "different tiling types exhibit distinct **conformal and shearing**
deformation" — asserted, never quantified. 2026 §5.1 asserts conformality for all `θ` empirically.

**Novelty vs field.** R7 uses conformal maps precisely because the triangular linkage is isotropic;
anisotropic linkages are not in that pipeline. R14 (Dang 2022, spherical) and R19 (Jiang 2024) are
different mechanisms. `screen_bundle.md` S4 ("exact conformality") NOT FOUND.

**Why it might work / sketch.** `J(θ)ᵀJ(θ)` has eigenvalues determined by `tr(KᵀK)`, `tr K`,
`det K` and `θ` alone, so the achievable `(λ₁, λ₂)` region is the image of an affine set under a
low-degree polynomial map — cheap to sample and to bound. The alignment requirement is that the
target's metric distortion tensor and the pattern's achievable one share eigenvectors, which is a
2-form matching condition, i.e. a cross-field problem of the kind Meekes & Vaxman 2021 solves.

**Kill experiment.** `kill_i8.cpp`. For the 16 patterns, sample `𝒦` (32 random `t`) × 64 angles and
plot the achievable `(λ₁, λ₂)` cloud. Then compute, for `pringles2.obj` (a saddle) and
`hemisphere2.obj`, the required `(λ₁, λ₂)` field from a discrete metric-matching map.
**Kills the idea if** the achievable clouds of all 16 patterns are indistinguishable up to rotation
(then anisotropy carries no selection signal), or if the isotropic sub-family already covers every
target's demand (then I8 collapses into I2). Runtime: minutes for the clouds; the demand field
needs I5's map.

**If it survives, the demo.** 16 stretch-region plots with each target's demand overlaid; a
"which pattern for which curvature sign" map. Hero: a saddle realized by an anisotropic pattern
that the isotropic budget of I2 says is infeasible.

**Risk.** Likely subsumed by I2: if the areal budget alone already separates the paper's cases,
anisotropy is a refinement nobody needs. Also, `𝒦` is measured at fixed lattice (same caveat as I2).

**Effort.** Derivation M. Code M.

---

## I9 — Audit: is 2025's `θ_max = min(2π − α_i − α_j)` exact on its own library?
**Type: characterization (supporting; a negative result is still a result).**

**Claim.** For every uniformly deployable split-free tiling, the paper's local hinge bound is the
exact first-contact angle, i.e. no global (non-adjacent) collision binds first — or it is not, and
we can say exactly which of their 16 patterns and which optimized variants break it.

**Gap targeted.** 2025 §12 assumption 7, verbatim in `notes/paper_2025.md`: "Non-adjacent (global)
self-collision during deployment is **not modeled**." 2025 §15 open question 10 asks whether
`θ_max` is even recomputed during their optimization, since the optimizer changes the interior
angles it depends on.

**Novelty vs field.** This is an audit of a stated, untested assumption in the source paper, using
machinery (T4.2″, K2a: exact vs bisection on 187/187 to `2e−10`) that no prior work has.
F24 already found their own `θ_max` routine buggy (merge_close_verts fuses hinge duplicates,
reporting 0.067 on snub square where the truth is ≈ 1.65).

**Why it might be true / sketch.** T4.4 proves `Θ_max ≤ min(min β, π)` and measured equality on
8 split-free patterns to `8.9e−16` — but the `≥` direction is explicitly **not proved**
(T4.4, `[N]` tag). Regular patterns are the easy case; the paper's *optimized* tilings have
irregular tiles (Fig. 14 red arrows), which is exactly where a distant vertex-into-edge contact
could bind first.

**Kill experiment.** `kill_i9.cpp`, after P0: run `exact_theta_max` (full candidate list, no
pruning, `theta_hi = π`) against `min_e β_e` on all 16 imported patterns at their as-shipped
geometry, and on 200 random perturbations of each within `𝕏`. **Kills the idea if** equality holds
on all 16 and all 200 perturbations of each — then the paper's assumption is safe and the audit
yields nothing publishable on its own (it becomes a lemma inside I1/I6). Runtime: minutes.

**If it survives, the demo.** The list of their own patterns where the local bound over-reports,
with the witnessing (vertex, edge) pair rendered. Hero: an optimized tiling that self-intersects at
an angle their pipeline declares safe.

**Risk.** Most likely outcome is "no violations on the regular library", i.e. the negative. That is
still worth one paragraph, but it cannot carry a section. Second risk: their shipped JSONs are
*unoptimized* patterns, so the interesting geometry (Fig. 13/14 outputs) is not in the repo at all
and we would have to re-run their optimizer to find it.

**Effort.** Derivation S (done, T4.4). Code S (after P0).

---

## I10 — A rank test for plain (non-uniform) deployability
**Type: characterization.**

**Claim.** For a tiling with deployment-**un**friendly vertices, let `θ ∈ R^{|E_hinge|}` be the
per-hinge angle field and `G(θ) = 0` the closure system at every hole. Then the tiling is deployable
in the sense of 2025 Def. 3.1 iff `G` has a solution branch through `θ = 0`, and a sufficient
first-order test is `rank DG(0) < |E_hinge|` together with a non-degenerate second-order term;
uniform deployability is the special case `1 ∈ ker DG(0)`.

**Gap targeted.** 2025 §6 first limitation and `notes/paper_2025.md` §15 open question 12: Prop. 4.1
characterizes *uniform* deployability only; Fig. 10 exhibits a tiling that is deployable
non-uniformly and the formulation for that case sits in the unavailable Supplement C. STATE.md D5
dropped the round-1 version of this (A8) as already-shown-by-Fig.-10; the difference here is that a
**rank test** is not shown by Fig. 10.

**Novelty vs field.** R12 (Dang 2021, F20) proves for **quads** that rigid deployability ⟺ a
deployed state with all voids parallelograms — so for quads uniform and plain deployability
coincide and this idea can only be about non-quad tilings. R21 is the source paper. Tay–Whiteley
motion assignments (S6) give the linkage mobility but not a closure-branch criterion.

**Why it might be true / sketch.** `DG(0)` is `L` in disguise: differentiating the hole closure
`Σ_{e ∈ C} R(θ_e-partial-sums)(x_dst − x_src) = 0` at `θ = 0` gives, per hole, a linear map from the
angle field to `R²`, whose matrix is `J_rot` times the same incidence-times-ownership product
`R · D` of T7(i) but with columns indexed by hinge **edges** instead of vertices. Uniform
deployability is `L X = 0` (Eq. 2 for the constant field `θ_e ≡ θ`); the non-uniform question is
whether the same operator has any other kernel direction that integrates to a real branch.

**Kill experiment.** `kill_i10.cpp`. Build `DG(0)` on (a) the 16 imported patterns, (b) a
hand-authored copy of 2025 Fig. 8 (their **non**-deployable four-line-family example, which is
reconstructible from the figure description in `notes/paper_2025.md` §7.1), and (c) a hand-authored
Fig.-10-like tiling with a few unfriendly vertices. **Kills the idea if** the rank test declares
Fig. 8 deployable, or declares a Fig.-10-like tiling non-deployable. Runtime: minutes; the cost is
authoring the two test tilings, which is the real risk.

**If it survives, the demo.** A classification of the 16 patterns into uniform / non-uniform /
non-deployable, with the non-uniform ones' angle fields solved by Newton continuation and rendered.
Hero: a deployed non-uniform tiling with its per-hinge angle field colour-coded, i.e. Fig. 10 done
computationally instead of by hand.

**Risk.** A first-order rank test is necessary, not sufficient — the branch can fail to integrate,
and 2025's Fig. 8 example is presumably exactly a case where the linearisation is fine and the
geometry jams. Without the second-order analysis this is a heuristic wearing a theorem's clothes.
Also, we cannot verify Fig. 8's geometry: it is reconstructed from a figure description, not from
data.

**Effort.** Derivation M. Code M.

---
# Self-attack — the strongest rejection I can write against my top 3

## Against I2 (curvature budget)

> The paper's central object, the "budget" `B`, is computed from a **fixed-lattice periodic**
> Jacobian, and is then used to make claims about **spatially varying, aperiodic, three-dimensional**
> designs. The authors acknowledge the gap and call it a "homogenisation ansatz", but an ansatz that
> is never validated against the artefact it is supposed to explain is a hypothesis, not a model.
> The one quantitative check offered (§I2 kill (c)) tests the *target's* conformal factor, not the
> deployed structure's metric, so it cannot detect the ansatz failing. Worse, the paper's own
> derivation shows the deployed cell metric is `J(θ)ᵀJ(θ)` with a **single global `θ`**; the 2025
> pipeline it claims to explain enforces reconfigurability only as a soft penalty and its deployed
> states are genuinely non-planar. The `(1 − Ω/4π)^{-2}` threshold is a classical fact about
> stereographic projection dressed up as a kirigami result — the interesting number, `B`, is the one
> that is not derived but measured, and it is measured on the wrong object. Finally, the paper leans
> on a single empirical negative (the triangular hemisphere) whose evidence the authors themselves
> report they could not see, since it lives in a supplement they do not have. One unseen negative is
> not a validation set.

## Against I3 (design space = boundary)

> This is Euler's formula. `k = |V| − H`, substitute `H = |E_hinge| − |F| + c(Γ)`, substitute
> `|E_∂| = |V_∂|`, collect. The authors say so themselves — "one line from T7(vi)". The claimed
> non-obvious corollary, that the interior of a complete-cut tiling is design-rigid, is *already
> visible* in their own round-1 measurement F15 (`dim = |E_split|`, and `E_split = ∅` for a complete
> cut gives zero interior freedom immediately). The `Θ(√n)` scaling law is then a statement about
> boundary vertex counts of planar patches, which is not a fact about kirigami. The one genuinely
> empirical hook — that this explains why 2025's Fig. 12 shows non-monotone `E_r` as density rises —
> is a post-hoc reading of four data points from a figure the authors transcribed from a PDF, with
> no access to the underlying statistics (they are in the unavailable Supplement F). A TOG paper
> cannot rest a section on four numbers read off someone else's plot.

## Against I1 (closure margin)

> The proposed criterion is `m = min β − θ_c ≥ 0`, and the authors' own headline table shows
> `m ≡ 0` — **exactly zero, to `1e−15`** — for all three regular families they use to motivate it.
> A criterion that is identically degenerate on every example in the motivating figure is not a
> criterion. The rescue, "take the supremum over the shape space", then requires the shape space to
> be non-trivial, and the companion result I3 says it is supported entirely on the boundary — so the
> supremum is being taken over boundary perturbations of a patch whose interior is rigid, and any
> positive margin found will be a boundary artefact that vanishes in the large-patch limit the
> method is actually used in. Meanwhile the honest explanation for the triangular pattern's failure
> is sitting in the source paper's own fabrication section: closure needs `2π/n` of face rotation,
> which is 120° for triangles against 90° and 60°, and the paper states outright that thick hinges
> "limit the achievable rotation angles". The geometric story is a re-derivation of a mechanical
> one.

---

# Ranked list

| # | idea | one-line justification |
|---|---|---|
| 1 | **I2 curvature budget + `(1−Ω/4π)^{-2}` threshold** | The only idea that answers 2025's own headline open problem (M5) with a closed form, gives a number (`4` for a hemisphere) that can be wrong, and reuses `AchievableSet`, which already exists and is validated. |
| 2 | **I3 design space = `|V_∂| + |E_split| + 1 − c(Γ)`** | Cheapest to kill (minutes, no P0), already half-validated by F15, and its corollary — every 2025 tiling has `Θ(√n)` design freedom regardless of pattern — reframes both papers' inverse-design sections. |
| 3 | **I1 closure margin / gapless criterion** | Directly targets the triangle-vs-quad-vs-hex fact the 2025 paper states and cannot show; `θ_c = 4π/n` is exact and matches M1 to `1e−15`; ranked third only because `m ≡ 0` on the regular families makes the separation rest on I3's shape space. |
| 4 | **I4 max-area at `θ_c/2`** | A two-line theorem that closes 2025 open question 17 and gives a universal normalised deployment curve; too small to stand alone, ideal as a corollary inside I1/I2. |
| 5 | **I6 certified-deployable inverse design** | The strongest *algorithmic* deliverable and the one both papers most obviously lack, but K2b (F27) already showed the authors' heuristics beat ours on margin, so the claim must be certification only. |
| 6 | **I9 audit of 2025's `θ_max`** | Minutes to run, tests a stated untested assumption, and gives a real finding if it fails; most likely yields a clean negative that becomes a lemma. |
| 7 | **I5 automatic tiling selection** | Exactly the paper's stated future work and a satisfying demo, but the validation set is one unseen negative and the discrete conformal map is real code cost. |
| 8 | **I8 anisotropy budget** | Well-posed and cheap to sample, but probably subsumed by I2's areal budget; keep as a refinement. |
| 9 | **I10 rank test for non-uniform deployability** | Highest-value target (2025's first-listed limitation) but the test is only first-order and the two decisive test tilings must be reconstructed from figure descriptions. |
| 10 | **I7 singularity lower bound** | The most attractive statement in the list and the least defensible: the cap-decomposition argument is not a proof, and it is entirely downstream of I2. |

---

# What I could not verify

- I did **not** run any code. Every number quoted as measured comes from files already in the repo
  (`results/kill/k7/*.csv`, `results/kill/KILL_REPORT.md`, `STATE.md` F14–F30, `derivations/core.md`);
  every number I derived myself is marked as such (`θ_c = 4π/n`, the free-boundary dimension count,
  the `(1 − Ω/4π)^{-1}` cap factor) and each has a kill test above.
- `notes/repo_2025.md` was a heading skeleton at the time of writing (§0–§7 empty), so my reading of
  their pattern JSON format comes from the files themselves, and my statements about their C++
  pipeline come from `notes/paper_2025.md` plus the header comments in
  `code/src/method/periodic_jacobian.hpp` — **not** from reading their C++.
- The triangle/quad/hex hemisphere fact (M5) is a *quotation* from the 2025 paper. Its evidence
  (Fig. F.4) is in a supplement not available locally, as `notes/paper_2025.md` §14 records. Any
  claim of "explaining" it inherits that limitation and must say so.
- I did not check the round-1 persona files, so I may be re-proposing something already written
  there; the closest known overlap is the optimizer persona's B20 (hole-area harmonic, second closed
  angle) and R4 (achievable periodic Jacobians), which K7 has already tested — I1/I2/I4 build **on
  top of** those, and I say so rather than claiming them.
