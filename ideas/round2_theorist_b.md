# Round 2 — Ideator, persona: THEORIST (second pass, `_b`)

**Files read this session.** `specs/common_preamble.md`, `specs/ideator.md`, `specs/ideator_round2.md`,
`derivations/core.md` (§0, T1 in full, T3, T4.4, T5.3, T7, T8), `STATE.md` (F1–F31, dead ends, session
log), `results/kill/k7/*.csv` (`k7_c3_all.csv` in full, 50 rows; headers of `k7_main.csv`,
`k7_c4_sweep.csv`, `k7_nu_*.csv`), `code/README.md` (build, JSON contract, library layout),
`code/src/method/zero_plus.hpp` and `code/src/method/periodic_jacobian.hpp` in full,
`notes/field_kirigami.md` and `notes/field_rigidity.md` (R-row tables only),
and the three sibling round-2 files `ideas/round2_theorist.md`, `ideas/round2_inverse.md`,
`ideas/round2_adversary.md`. I also **wrote and ran one program**,
`derivations/scratch/check_b1.cpp`, which verifies the central identity of this file (§0.g); it links
against `code/src/core` and modifies nothing under `code/`. I did **not** re-read the 2026/2025 paper text; every paper citation is
via the F-facts, `code/README.md`, or a sibling file, and is marked where that matters.

**My angle, and why it is not one of the other three.** The three round-2 files all attack the
emptiness (F25/F30) from the *shape space*: theorist-a asks which `σ` and which parity class make it
empty, adversary asks for a convex certificate that it is empty, inverse asks what a *pattern* can
wrap. None of them has an **accounting identity**: a conservation law that says where the opening of
every hole is paid for. This file supplies one, derives everything else from it, and uses it to
explain the single most important unexplained number in the project — that the *same* random Voronoi
graphs are `0/400` on a fixed-boundary patch (F30/K6) and `11/12` zero-plus feasible with `4/12`
fully certified on a **torus** (`results/kill/k7/k7_c3_all.csv`, `random_in_K` rows). That contrast is
not in `STATE.md` as a fact yet, and no sibling file mentions it.

---

## 0. The identity everything below is built on: hole area is a boundary functional

Everything in this file is a corollary of one computation, so I do it once, in full, in the code
convention of `derivations/core.md` §0.7 (`R_f = c I − σ_f s J`, `t_f = 2 s J u_f`, `c = cos(θ/2)`,
`s = sin(θ/2)`, `u` the face potential of T1 Step 3 with `u_g − u_f = σ_g x_{src(e)}` across a hinge
edge traversed `f → g`).

**(0.a) The area of any closed walk in `M′` is a first harmonic, computed edge by edge.**
Let `W` be a closed walk in `M′`; every directed edge `(a → b)` of `W` lies in exactly one face `f`,
so `y_a = R_f x_a + t_f` and `y_b = R_f x_b + t_f` with the *same* `R_f, t_f`. Then

```
   det(y_a, y_b) = det(R_f x_a, R_f x_b) + det(R_f x_a, t_f) + det(t_f, R_f x_b)
                 = det(x_a, x_b) + det( R_f (x_a − x_b), t_f ) .
```

Write `d := x_a − x_b`. Using `det(d, J u) = ⟨d, u⟩`, `det(J d, J u) = det(d, u)`, `2sc = sin θ`
and `2s² = 1 − cos θ`:

```
   det( R_f d, 2 s J u_f ) = 2s [ c ⟨d, u_f⟩ − σ_f s det(d, u_f) ]
                           = sin θ · ⟨d, u_f⟩  −  (1 − cos θ) · σ_f det(d, u_f) .
```

Summing over `W` and halving, and using that the walk is **degenerate in the flat state** (a hole
preimage collapses to a doubly-traversed tree at `θ = 0`, so `Σ_W det(x_a, x_b) = 0`):

```
 ┌──────────────────────────────────────────────────────────────────────────────────────┐
 │   A_C(θ)  =  a_C · sin θ  −  b_C · (1 − cos θ) ,      A_C(0) = 0 ,                    │
 │                                                                                      │
 │   a_C = ½ Σ_{(a→b) ∈ ∂C, face f}  ⟨ x_a − x_b , u_f ⟩ ,                               │
 │   b_C = ½ Σ_{(a→b) ∈ ∂C, face f}  σ_f · det( x_a − x_b , u_f ) .            (B.1)     │
 └──────────────────────────────────────────────────────────────────────────────────────┘
```

for every hole preimage `C`, with `∂C` its boundary walk oriented with the hole on the left. In the
`(p, q, r)` notation of T3, `A_C` has `p = −b_C`, `q = +b_C`, `r = a_C`, i.e. `p + q = 0` — the same
degenerate class as the split-gap harmonic of T5.3, and for the same reason (the flat state is a root).
Individual terms of (B.1) depend on the origin and on the root face of `u`; the sums do not, because
`Σ_{∂C} (x_a − x_b) = 0` for a closed walk.

**(0.b) Per-edge form: hinge terms and split terms.** Both copies of an interior edge `e` lie on the
same preimage walk (F11) and are traversed in opposite senses. For a hinge edge `e` with
`d_e = x_{dst} − x_{src}`, between `f` and `g`:

```
   contribution to a_C  =  ½ ⟨ d_e , u_g − u_f ⟩  =  ½ σ_g ⟨ d_e , x_{src(e)} ⟩ ,       (B.2h)
```

and for a split edge `e = {a, b}` between `f, g` with `Δu_e = u_g − u_f`:

```
   contribution to a_C  =  ½ ⟨ d_e , Δu_e ⟩  =  ½ r_e ,                                 (B.2s)
```

where `r_e` is **exactly** the `r` of `derivations/core.md` T5.3, hence a positive multiple of the
`q_e = det(dS_e, d_e)` that `code/src/method/zero_plus.hpp` computes and that K5/K6 measured. So:

```
   a_C  =  ½ [  Σ_{e ∈ C ∩ E_hinge} σ_{g(e)} ⟨d_e, x_{src(e)}⟩   +   Σ_{e ∈ C ∩ E_split} r_e  ] . (B.3)
```

**The whole measured emptiness (`min_e q_e < 0`, F30) is therefore a statement about how much
first-order area each hole is allowed to open.**

**(0.c) The divergence theorem.** Sum (B.1) over *all* preimages (holes **and** notches — they
partition `E_hinge ⊔ E_split` by F11). Each interior edge copy appears once; each face `f` contributes
`Σ_{edges of f} ⟨x_a − x_b, u_f⟩ = ⟨ Σ (x_a − x_b), u_f ⟩ = 0` because a face polygon is closed, and the
same telescoping kills the `b` sum. What is left is the border of `M`:

```
 ┌──────────────────────────────────────────────────────────────────────────────────────┐
 │   Σ_C a_C  =  ± ½ Σ_{border edges (a→b) ⊂ face f}  ⟨ x_a − x_b , u_f ⟩  =:  B(X) ,    │
 │   Σ_C b_C  =  ± ½ Σ_{border edges}  σ_f det( x_a − x_b , u_f ) .            (B.4)     │
 └──────────────────────────────────────────────────────────────────────────────────────┘
```

I call `B(X)` the **expansion budget**. The global sign is fixed by the face-boundary orientation
convention and I do not assert which; the kill test pins it (see B1).

**(0.d) The periodic case, and an independent confirmation that (B.1) is right.** On a torus there is
no border, so (B.4) would give `Σ_C a_C = 0` — false, because tori do open. The resolution is that `u`
is *not single-valued* on the quotient: `code/src/method/periodic_jacobian.hpp` derives
`u_{f+t} = u_f + w_t + σ_f t/2`, so the preimage walks that cross the seam close only up to the
monodromy `w_t`, and the telescoping in (0.c) leaves exactly that defect. The header also states,
independently of anything here, that the hole area per cell is `det(P₀)(det J(θ) − 1)` with
`J(θ) = c I + s K`. Expanding,

```
   det J(θ) − 1 = c² + cs·tr K + s²·det K − 1 = ½ (tr K) sin θ − ½ (1 − det K)(1 − cos θ) ,
```

which is (B.1) with

```
   Σ_C a_C = ½ det(P₀) · tr K ,      Σ_C b_C = ½ det(P₀) · (1 − det K) .          (B.5)
```

**The budget of a periodic pattern is `½ det(P₀) tr K`, and nothing else.** That two independently
derived routes — my walk computation and the K7 monodromy computation, verified numerically in
`kill_k7` — give the same two coefficients is the strongest evidence I have that (B.1) is correct
before running anything.

**(0.e) The one fact this file exists to explain.** In `results/kill/k7/k7_c3_all.csv` the target
`diag_1_-0.5` (a `K` with `tr K = 0.5`) is certified on **0 of 25** rows, with `min_q` negative on
**25 of 25**; the `random_in_K` targets are certified on **12 of 25** with `min_q > 0` on 20 of 25.
Random Voronoi **tori** — the same generator family that is `0/400` on fixed-boundary patches
(F30/K6) — are zero-plus feasible on **11 of 12** and fully certified on **4 of 12**. By (B.5) the
difference between the two target classes is a difference in `tr K`, i.e. in budget; by (B.4) the
difference between the torus and the fixed-boundary patch is that the torus has 4 free parameters in
`K` (`dimK = 4` on every row of `k7_c3_all.csv`) buying budget, while a patch with `B X = T` pinned
has its budget fixed by boundary data it is not allowed to move. **Emptiness is a budget statement,
not a randomness statement.**

**(0.f) Hand checks of (B.1)/(B.5) against `results/kill/k7/k7_main.csv`.** The columns `c4_q`,
`c4_r`, `c4_thetac` are the pattern's own hole-area harmonic in the form
`A(θ) = c4_r · sin θ + c4_q · (cos θ − 1)`, so `a = c4_r` and `b = c4_q` in my notation and
`θ_c = 2 atan2(a, b)`. All values below are read from that file (33 rows, measured, `c4_area_fit_err`
and `c4_geom_err` columns exist and are the fit residuals).

| pattern | `n_split` | `a = c4_r` | `b = c4_q` | `θ_c` | reading via (B.5) |
|---|---|---|---|---|---|
| `squares_2x2` | 0 | 4.000000 | −0.000000 | `π` | `det K = 1`, `tr K = 2` ⟹ `K = I`, `J(θ) = (c+s) I` |
| `kagome_2x2` | 0 | 3.000000 | −0.000000 | `π` | `det K = 1` again |
| `triangles_2x2` | 0 | 6.000000 | −3.464102 | `4π/3` | `det K > 1` |
| `hexagons_2x2` | 4 | 6.000000 | +3.464102 | `2π/3` | `det K < 1` |
| `voronoi_torus_*` (12 rows, `n = 20 … 200`) | 20…204 | 33.8 … 39.2 | 34.8 … 39.2 | 1.86 … 2.10 | see B7 |

* **Square grid, by hand.** `b = 0` ⟹ `det K = 1`; `a = ½ det(P₀) tr K = 4` with `det P₀ = 4` ⟹
  `tr K = 2`; with `det K = 1` and `tr K = 2` the only diagonalisable `K` is `I`, so
  `J(θ) = (cos(θ/2) + sin(θ/2)) I` — an exact isotropic similarity — and the hole area per cell is
  `det(P₀) sin θ = 4 sin θ`, i.e. `sin θ` per hole. At `θ = π/2` each hole has area 1, exactly the
  area of a unit tile: the textbook rotating-squares open state in which void area equals tile area.
  This is a genuine, independent confirmation of (B.5).
* **Triangles and hexagons are exact duals in `b`.** `a = 6` for both, `b = ∓2√3`, and
  `θ_c(triangles) + θ_c(hexagons) = 4π/3 + 2π/3 = 2π` **exactly**. That is a prediction to test on
  every dual pair, not a coincidence I have proved (idea B8).
* **A single split edge.** By (B.2s) it contributes `½ r_e = ½⟨d_e, Δu_e⟩` to its hole and to nothing
  else; by T5.3 and `zero_plus.hpp` this is a positive multiple of the measured `q_e`. It is the only
  term of (B.3) whose sign is free, which is the algebraic content of "split cuts are the only source
  of design freedom **and** the only source of collision loss" (`specs/ideator_round2.md` item 2).
* **A guess of mine that the data killed, recorded because it is the obvious wrong turn.** I first
  argued that a *pure* hole `C_v` (in-edges of one interior vertex, no split edge) is the `2 sin(θ/2)`
  scaled image of the fixed polygon `{J u_f}`, hence `A ∝ (1 − cos θ)` and `a_C = 0`, hence
  "split-free patterns have zero budget". The `squares_2x2` row (`a = 4 ≠ 0`) refutes it. The error:
  the boundary walk of `C_v` also runs along the two **fixed hinge points** `x_{src(e)}` of its
  in-edges, so the hole is not a scaled copy of anything. (B.1) is unaffected — only the corollary was
  wrong — but anyone reusing this file should not re-derive that shortcut.

**(0.g) (B.1) IS NOW MEASURED, not only derived.** I wrote and ran
`derivations/scratch/check_b1.cpp` (build line in its header; links only against `code/src/core/*.cpp`,
modifies nothing under `code/`). On **10 patterns** (`squares`, `triangles`, `hexagons`, `kagome`,
`snub_square`, `truncated_square`, `t3_4_3_12`, at two sizes for three of them), solved to `X₀` in `𝕏`
with a fixed boundary, it (i) rebuilds `u` by BFS over `Γ` and checks its closure on every hinge edge,
(ii) takes the hole cycles from `holes_geometric_cycles` at `θ = 0.05`, (iii) shoelaces each hole at
`θ ∈ {0.05, 0.2, 0.4, 0.7, 1.0, 1.3}` and least-squares fits `A_C = a sin θ − b(1 − cos θ)`, and
(iv) compares the fit with the closed form (B.1) evaluated along the *same* directed walk.

| check | worst over the run | count |
|---|---|---|
| `u` closure residual on every hinge edge | `2.70e−15` | 10 graphs |
| **(B.1) closed form vs fitted `(a_C, b_C)`, relative** | **`8.72e−14`** | **471 holes** |
| first-harmonic fit residual (is `A_C` really `p + q cos + r sin` with `p + q = 0`?) | `1.60e−14` | 471 holes |
| **(B.2s) `Δu_e` recovered from the deployed split offset vs `u_g − u_f`** | **`6.79e−15`** | 58 split edges |
| (B10) numerically located area maximum vs `½ θ_c` | `6.40e−04` | grid-limited: the scan step is `π/2000 = 1.57e−3`, so this is half a grid cell and consistent with exact |

So (B.1), the first-harmonic form with `A_C(0) = 0`, and the split-edge term (B.2s) are **verified**;
(B.4) and (B.3)'s hinge term are **not** — the geometric hole tracer returns holes only, never
notches, so the global sum could not be closed on a patch in this run. That is B1's remaining kill
test.

---

# 1. The ten ideas

## B1 — The expansion-budget identity: hole opening is a boundary functional

**1. Title / type.** *The budget identity `Σ_C a_C = B(X)`* — **theorem**.

**2. Claim.** For every `X ∈ 𝕏` on a disk-topology patch, the total first-order hole-opening rate is
a sum over the **border** edges alone,
`Σ_C A_C′(0) = Σ_C a_C = ± ½ Σ_{border (a→b) ⊂ f} ⟨x_a − x_b, u_f⟩ =: B(X)`, and on a torus it is
`½ det(P₀) tr K` per cell; combined with (B.3),
`Σ_{e ∈ E_split} r_e = 2 B(X) − Σ_{e ∈ E_hinge} σ_{g(e)} ⟨d_e, x_{src(e)}⟩`. Both sides are explicit
quadratic forms in `X`, hence quadratic polynomials in the design coordinates `t`.

**3. Gap targeted.** 2026 §6/§7 limitation (i) as recorded in F8: "developing geometric
characteristics for favorable deployment behavior directly from the embedding remains an open
problem". This is such a characteristic, and it is an identity rather than a heuristic. It also makes
2026 §4.2's unsupported "balanced, small holes" principle precise in a different way from adversary
Idea 5: the constraint is not per-hole size but a **global conservation law** with a boundary source.

**4. Novelty vs field.** (a) `field_rigidity` **R11** Guest–Hutchinson and **R16/R17**
Borcea–Streinu give periodic-framework deformation theory in terms of the lattice Gram matrix; they
have no per-hole decomposition and no cut structure, so the split/hinge split of (B.3) is not
available there. (b) `field_kirigami` **R12/R17** (Dang 2021, Dudte 2023) are quad-only and
parametrise vertex positions linearly; neither writes a hole-area harmonic. (c) `STATE.md` session log
line for the early Screen-Scout records verdict **S7** — "hole count / area harmonic / 2nd closed
angle **NOT FOUND**" in the literature — which is exactly this object. (d) Internally,
`periodic_jacobian.hpp` has the *periodic total* `det(P₀)(det J − 1)`; it does not have the per-hole
form, the boundary functional, or the split/hinge decomposition, and neither does any sibling
round-2 file.

**5. Why it might be true / sketch.** §0 above, in full: the walk computation (0.a), the per-edge
regrouping (0.b) which uses F11 (both copies of an edge lie on one preimage walk) and T1.3
(`u_g − u_f = σ_g x_src`), and the telescoping (0.c) which uses only that a face polygon is closed.
The periodic case (0.d) reproduces `periodic_jacobian.hpp`'s independently derived coefficients,
which is a non-trivial consistency check.

**6. Kill experiment.** *Partly run already — see (0.g): the per-hole half of this test passed at
`8.72e−14` over 471 holes on 10 patterns, and what remains is the global sum (B.4).* The full version
is a new app `kill_b1.cpp` on the existing 16-graph corpus plus 30 random
Delaunay/Voronoi patches (`kiri_gen`), ≈ 150 lines: (i) deploy at 6 angles, shoelace each hole
polygon from the `"holes"` cycles that `kiri_deploy` already dumps, least-squares fit
`A_C = a_C sin θ − b_C(1 − cos θ)`; (ii) compute `a_C, b_C` from (B.1) using `u` rebuilt by BFS over
`Γ` (30 lines, the same BFS as `derivations/scratch/check_t1_t2.cpp` C1); (iii) compute the
right-hand side of (B.4) from border edges only. **Kills the idea:** any relative deviation
`> 1e−10` in (i)-vs-(ii) on any hole, or in (ii)-vs-(iii) summed. Also assert `A_C(0) = 0` and
`p + q = 0` per hole. Runtime: minutes. A sign flip in (B.4) is *not* a kill — it pins the
convention.

**7. If it survives, the demo.** One figure per pattern: the per-hole budget `a_C` drawn as a
colour map on the flat design, with the border edges carrying `B(X)` drawn as arrows — the picture
literally shows where the opening is paid for. Table of `Σ a_C` vs `B(X)` vs `½ det(P₀) tr K` on 46
designs.

**8. Risk.** The per-hole identity is now measured, so the surviving risk is only in (B.4): the
identity is true but *vacuous as a design tool* because `B(X)` is itself a free
quadratic form that the null space can make large — i.e. it constrains nothing on a patch with a free
boundary. B2 is the version of the claim that can actually fail.

**9. Effort.** Derivation **S** (done above). Code **S** (one app, reuses the hole cycles already
dumped).

---

## B2 — The budget bound: an `O(|E|)` upper bound on the deployment margin, before any null space

**1. Title / type.** *Budget bound on `min_e q_e`* — **theorem + screening algorithm**.

**2. Claim.** For any `X ∈ 𝕏` with `E_split ≠ ∅`,

```
   min_{e ∈ E_split} r_e   ≤   ( 2 B(X) − W(X) ) / |E_split| ,
   W(X) := Σ_{e ∈ E_hinge} σ_{g(e)} ⟨ d_e , x_{src(e)} ⟩ ,
```

so `2B(X) < W(X)` is a **certificate of `Θ_max = 0`** costing one pass over the edges, with no
linear system, no null space and no SVD. On the K6 population (400 designs, `0/400` certified) the
prediction is `2B − W < 0`; on the 8 authored tilings and on the 12 certified K7 torus designs the
prediction is `2B − W > 0`.

**3. Gap targeted.** Same 2026 §7 limitation (i); and directly the hole in this project's own
emptiness claim, which is currently the measurement `0/400` (F30) with no accompanying reason. The
sibling files attack that hole with a per-instance convex certificate (adversary Idea 2, theorist-a
Idea 2 — both are `O(k²)` SDP/Farkas computations on the shape space). This one is a closed-form
scalar inequality and is `O(|E|)`; where it applies it is strictly cheaper, and where it does not the
convex certificates are still needed. The two are complementary, not competing.

**4. Novelty vs field.** (a) `field_rigidity` **R1** Maxwell–Calladine counts are counts of degrees
of freedom, not of areas, and F13 records that they predict rigidity for every working pattern —
they cannot see this. (b) **R8** Connelly–Whiteley second-order rigidity is the right general tool
for special-position flexes but gives no scalar of this form. (c) `field_kirigami` **R22** IsoGami
filters `3^m` joint assignments by simulation; a filter is not a bound.

**5. Why it might be true / sketch.** Averaging: `min ≤ mean`, applied to (B.3) summed over all
holes and combined with (B.4). Every quantity is explicit. The bound is *tight* exactly when all
`r_e` are equal, which is the symmetric case (hexagons), so it should be near-tight on tilings and
loose on random graphs — which is the regime where it is used as a *negative* certificate, where
looseness only costs coverage, never correctness.

**6. Kill experiment.** `kill_b2.cpp` reusing `kill_k6.cpp`'s design population verbatim (the CSVs
in `results/kill/k6/shards_v4/` record the graph seeds, so the same designs are reproducible):
compute `B(X)`, `W(X)`, `min_e q_e` for all 400 K6 designs, the 8 tilings, and the 33 K7 periodic
designs. **Kills the idea:** any design with `2B − W < 0` and `min_e q_e > 0` (the inequality is
violated — a hard mathematical kill), or `2B − W > 0` on all 400 K6 designs (the bound is true but
never binding, so it explains nothing — a soft kill that demotes B2 to a footnote of B1). Runtime
≈ 5 min; the certificate itself is microseconds per design.

**7. If it survives, the demo.** A scatter of `min_e q_e` against `(2B − W)/|E_split|` over ≈ 450
designs with the `y = x` bound line and every point below it, coloured by certified / not; plus a
table "graphs rejected in `O(|E|)` before any solve: N of 400". The hero example is a random Voronoi
patch rejected in microseconds that K6 spent an L-BFGS run per random start failing to repair.

**8. Risk.** The bound is true and non-binding: `B(X)` is large and positive on random patches
because `u_f` grows like the path length from the root face, so `B` is `O(n)` and `W` is `O(n)` with
an unfavourable constant that may go the wrong way. I have **not** checked the sign of `2B − W` on a
single real design; it is the whole content of the kill test.

**9. Effort.** Derivation **S**. Code **S**.

---

## B3 — `tr K` is the budget: a feasible region in Jacobian space for periodic inverse design

**1. Title / type.** *The budget half-space in `K`-space, and budget-projected target design* —
**theorem + algorithm**. **This is the pair I would build the constructive half of the paper on.**

**2. Claim.** For a periodic pattern, (B.5) says the *entire* first-order opening budget is
`½ det(P₀) tr K`, so the design variable `K ∈ 𝒦` (the affine achievable set that
`AchievableSet` in `kill_k7.cpp` already computes, `dimK = 4` on every K7 row) enters the
zero-plus feasibility problem only through the linear functional `tr K`. Consequently there is a
threshold `τ*(pattern) = W(X)/det(P₀)` such that `tr K < τ*` ⟹ `min_e q_e < 0` ⟹ `Θ_max = 0`, and
K7's C3 result splits along exactly that line: the target `diag(1, −0.5)` has `tr K = 0.5` and is
certified on **0 of 25** rows with `min_q < 0` on **25 of 25**, while `random_in_K` targets are
certified on **12 of 25**. The algorithm: given a target Jacobian `K_tgt`, solve the **budget-projected**
problem `min ‖K − K_tgt‖` over `K ∈ 𝒦 ∩ {tr K ≥ τ*}` instead of over `𝒦`, then certify.

**3. Gap targeted.** 2026 §5.1 and Eqs. (11)–(13), which design the periodic Jacobian by linear
constraints at `θ = 0` and say nothing about whether the resulting design deploys at all; and 2026 §7
limitation (i). It is also the fix for this project's own K7 C3 verdict (6/22 at the interim,
12/50 in the final CSV): the affinity theorem holds and the target is hit to `1e−16` on **50 of 50**
rows (`hit_err` column), so the design algorithm's only failure mode is that the hit target is
outside the deployable region — precisely what a budget constraint fixes.

**4. Novelty vs field.** (a) `field_rigidity` **R17** Borcea–Streinu periodic auxetics characterise
auxetic paths by a PSD condition on the lattice Gram matrix — a condition on the *lattice*, blind to
whether the cut duplicates collide; the budget threshold is a statement about the cuts. (b)
`field_kirigami` **R8/R11/R17** Choi–Dudte–Mahadevan design quad patterns for targets by
optimisation with the constraint enforced numerically. (c) **R22** IsoGami obtains its expansion range
from IPC-contact continuation per instance, not from a closed-form feasible region. (d) Inverse's I2
and I8 (`ideas/round2_inverse.md`) design *with* `𝒦` and ask what areal/anisotropy budget a pattern
delivers; neither has the trace threshold, and I2's `λ_A(θ, K) = det J(θ)` is exactly (B.5) rewritten
— so I2 and this idea share an identity and differ in what they do with it.

**5. Why it might be true / sketch.** (B.3) summed over the cell: `Σ_split r_e = det(P₀) tr K − W`,
with `W` the hinge sum, which is *fixed* by the flat combinatorics and the flat positions and does
not depend on `K` at first order along the achievable set. Then `min_e r_e ≤ (det(P₀) tr K − W)/n_split`
gives the threshold. Both `tr K` and `W` are already computable in `kill_k7.cpp` (the `c4_r` column
*is* `½ det(P₀) tr K`).

**6. Kill experiment.** `kill_b3.cpp`, reusing `kill_k7.cpp`'s pattern list and `AchievableSet`
verbatim: for each of the 25 `(pattern, target)` C3 rows, (i) record `tr K` at the solved design and
`τ* = W/det(P₀)`; (ii) re-solve with the extra linear constraint `tr K ≥ τ* + δ` for
`δ ∈ {0, 0.1, 0.5}·|τ*|`; (iii) re-run the existing certificate. **Kills the idea:** if certified rows
and uncertified rows are not separated by `tr K − τ*` (any certified row with `tr K < τ*`, or an
AUC below ≈ 0.8 over the 50 rows), the threshold claim is dead; if the constrained re-solve does not
raise the certified rate above 12/50 on the same targets, the *algorithm* is dead even if the theorem
survives. Runtime: K7's C3 shards ran in minutes each; this is one extra linear constraint.

**7. If it survives, the demo.** (a) A scatter of `tr K − τ*` against `min_e q_e` over 50+ rows with
the threshold line. (b) The design deliverable: ≥ 20 **random Voronoi tori** with prescribed
Poisson-ratio/conformality targets, certified deployable, next to the same graphs designed by the
unconstrained (2026 §5.1) route with their certified rate. (c) Hero: an arbitrary periodic Voronoi
graph deployed to a prescribed anisotropic Jacobian — an object neither Segall paper produces, since
2026's demonstrated periodic examples are authored tilings and its random-graph space is empty.

**8. Risk.** `W` is not exactly constant along `𝒦` (the flat positions move with `t`, so `W = W(t)`
is quadratic, not constant), which makes `tr K ≥ τ*` a *linearisation* of the true feasible set
rather than the set itself. If the variation of `W` over `𝒦` is comparable to `det(P₀) tr K`, the
threshold is not a threshold and the idea collapses to "add a penalty", which is not a theorem. The
kill test's step (i) measures exactly that variation.

**9. Effort.** Derivation **S**. Code **M** (one constrained re-solve inside an existing driver).

---

## B4 — The obstruction is the boundary condition, not the graph

**1. Title / type.** *Fixed-boundary emptiness vs free-boundary/periodic feasibility* —
**theorem (conditional) + the decisive experiment**.

**2. Claim.** The `0/400` emptiness of F30/K6 is an artefact of 2026 Eq. (4)'s **fixed** boundary,
not of graph randomness. Precisely: with the boundary rows `B X = T` dropped, the shape space of the
*same* random graph gains `|V_∂|` scalar dimensions (the count is inverse's I3,
`ideas/round2_inverse.md` §I3, `k = |V_∂| + |E_split|`), and by (B.4) those are exactly the
dimensions that move `B(X)` — the budget is a *border* functional and a fixed border is the only
thing that pins it. Prediction: zero-plus feasibility on the K6 population rises from `0/400` to a
strictly positive rate under a free boundary, and the K7 torus rate (`11/12` feasible, `4/12`
certified on random Voronoi tori) is the same phenomenon with the border removed altogether.

**3. Gap targeted.** F30 and the K6 dead-end entry in `STATE.md`, which currently reads as "the
random-graph shape space is empty" without qualification; and 2026 §4.3/Eq. (4), where the fixed
boundary is presented as a normalisation choice with no discussion of what it costs.

**4. Novelty vs field.** No field row is relevant: this is a statement about a specific linear system
in a specific paper. Internally it is new — theorist-a Idea 2, adversary Ideas 2 and 4, and inverse I3
all take the boundary condition as given; I3 counts the free-boundary dimension but never tests
feasibility in it.

**5. Why it might be true / sketch.** (B.4) shows `B(X)` depends on `X` only through border-edge
vectors and border-**face** potentials. Under Eq. (4) the border vertex positions are frozen, so the
border edge vectors `x_a − x_b` are constants and only the `O(|F_∂|)` potentials `u_f` remain free —
and those are path sums determined by the interior, i.e. not directly steerable. Free the border and
`B` acquires `|V_∂|` direct handles. This also retro-explains F23 (the authors' code adds constraint
rows for split-forest components touching the boundary, over-constraining exactly the border where
the budget lives) and F25 (`X_ini` injective 200/200 but `X₀` self-intersecting 142/200: the
projection has to buy opening somewhere and the border denies it).

**6. Kill experiment.** `kill_b4.cpp` = `kill_k6.cpp` with one change: build `[L]` without the
boundary rows `B`, pin only a single vertex and one direction to kill the rigid motions, then run the
identical `zero_plus_repair` with `w_corner`/`w_prox` at K6's final settings, on the identical 400
designs. **Kills the idea:** still `0/400` feasible. That is a clean, decisive, ≈ 20-minute run
(K6's final pass is already sharded, `results/kill/k6/shards_v4/`), and it is the single cheapest
high-information experiment in this file.

**7. If it survives, the demo.** The headline table of the paper: certified-deployable fraction on
identical random graphs under (fixed boundary, free boundary, periodic), with `B(X)` plotted for each.
The sentence it licenses — "the arbitrary-planar-graph design space is not empty; the published
boundary condition empties it" — is a genuine correction to 2026 rather than an extension.

**8. Risk.** Free-boundary designs may be feasible at `0⁺` and still uncertifiable, or may drift into
non-injective flat states with nothing to hold them (F17's inverted faces get *worse* without `B`).
Then the claim is true and useless, and the honest result becomes B3's periodic route instead.

**9. Effort.** Derivation **S**. Code **S** (a one-flag change to an existing driver).

---

## B5 — Spectral budget allocation: designing `t` by a generalised eigenproblem

**1. Title / type.** *Max-min split gap by a generalised Rayleigh quotient* — **algorithm**.

**2. Claim.** Every `q_e(t)` and every face area `a_f(t)` is an explicit quadratic polynomial in
`t` (`zero_plus.hpp` says so and `check.md` T5-a verifies it to `3.4e−14`), so on the homogenised
variable `w = (1, t)` the feasibility problem is "find `w` with `w₀ = 1` and `wᵀ M_e w > 0` for all
`e`". The claim: the relaxation `max_w min_e wᵀ M_e w / wᵀ N w` with `N` the proximity form
`‖X(t) − X₀‖²`, solved as a sequence of generalised eigenproblems (multiplicative-weights over the
simplex `λ`, exact eigen-solve per step), reaches strictly positive `min_e q_e` on designs where
K6's L-BFGS-on-softplus reached `0/400`, because the failure mode K6 reports is a **local** minimum
of a non-convex softplus sum and this relaxation is scale-invariant and start-free.

**3. Gap targeted.** 2026 Eqs. (7)–(9): a smooth penalty minimised from one start, with the
false-negative rate now measured at 9.47% (F28/K1c) and the silent re-closure documented. And this
project's own K6, whose negative result is currently confounded with its optimiser.

**4. Novelty vs field.** (a) `field_kirigami` **R22** IsoGami's filter is enumeration plus
simulation. (b) Adversary Idea 2 uses the *same* quadratic forms to certify **in**feasibility by
convex duality; this uses them to search for feasibility, and the two are the two sides of one
`λ_max(Σ λ_i M_i)` computation — running both gives either a design or a proof, which is a much
stronger paper claim than either alone. I regard B5 + adversary Idea 2 as one deliverable and would
cite it as joint.

**5. Why it might be true / sketch.** `λ ↦ λ_max(Σ_i λ_i M_i)` is convex on the simplex; its minimiser
either certifies emptiness (value `≤ 0`, adversary Idea 2) or produces a direction `w` with positive
margin on the weighted average, which is then rounded by one line search along `w`. Because K6's
objective is a sum of softplus barriers, its failures are start-dependent; this one has no start.

**6. Kill experiment.** `kill_b5.cpp` on the identical K6 400-design population, 200 MW iterations,
comparing (a) K6's reported feasibility, (b) the spectral relaxation's, (c) the certificate.
**Kills the idea:** `0/400` again, or feasible-but-uncertified everywhere. ≈ 20 min (each iteration
is one `4k × 4k` symmetric eigen-solve; `dim_null` up to 1017 on Voronoi (F19) makes the largest cases
the cost driver — cap the population at `|F| ≤ 400` if needed and say so).

**7. If it survives, the demo.** Feasibility rate vs K6 on 400 designs, and the pairing with the
emptiness certificate: every graph gets either a certified design or a proof there is none.

**8. Risk.** The relaxation is exactly the standard SDP relaxation of a non-convex QCQP and its
integrality gap here is unknown; the honest prior, given four independent strategies failed at
`0/400`, is that the set really is empty on fixed-boundary patches, in which case B5 returns
adversary Idea 2's certificate and no designs. B4 must be run **first**: if the fixed-boundary set is
empty for the structural reason B1/B2 give, B5 on that population is wasted compute.

**9. Effort.** Derivation **S**. Code **M**.

---

## B6 — Notch leakage: the boundary's share of the budget, and a size scaling law

**1. Title / type.** *Notches consume budget* — **theorem (accounting) + scaling law**.

**2. Claim.** The sum in (B.4) runs over **all** preimages, holes and notches alike (F11 partitions
`E_hinge ⊔ E_split`, and F16/`code/README.md` deviation 2 record that boundary-touching preimages are
notches carrying no row of `L`). Hence the budget available to genuine holes is
`B(X) − Σ_{notches} a_N`, and since notches are the preimages of the border, their number scales like
`|∂M| ~ √|F|`. Prediction: the certified fraction on a fixed random-graph family is an increasing
function of `|F|/|∂M|`, and the periodic (`0` notches) case is the limit of that family, which is
exactly the ordering observed between K6 (`0/400`, patches with boundary) and K7 C3 (`11/12`
zero-plus feasible on random Voronoi **tori**).

**3. Gap targeted.** F16 and 2026 §4.1 "border edges are excluded" (Alg. 1), which the paper treats as
a bookkeeping detail. It is not: it removes constraint rows *and* diverts budget.

**4. Novelty vs field.** Nothing in `field_kirigami`/`field_rigidity` treats the boundary of a finite
kirigami patch as a resource; **R10** Mitschke et al. and **R11** Guest–Hutchinson are periodic or
infinite. Internally, F23 documents the authors' *extra* boundary rows as over-constraint but never
connects the boundary to feasibility.

**5. Why it might be true / sketch.** Direct from (0.c): the telescoping leaves only border edges, and
the left side counts every preimage. Notches are typically *large* preimages (a split-forest component
that reaches the border swallows many hinge edges), so their `a_N` share is not proportional to their
count. That is the part that could go either way, and it is what the experiment measures.

**6. Kill experiment.** `kill_b6.cpp`: on 60 Delaunay/Voronoi patches with `|F| ∈ [100, 3000]`
(`kiri_gen`), report `B(X)`, `Σ_holes a_C`, `Σ_notches a_N`, the notch share, and zero-plus
feasibility. **Kills the idea:** the notch share does not decrease with `|F|/|∂M|` (Spearman `ρ > −0.3`),
or feasibility does not increase with it. ≈ 15 min.

**7. If it survives, the demo.** A single curve: certified fraction vs `|F|/|∂M|` over ≈ 60 patches,
with the periodic point plotted at the right end, and the notch-share curve beneath it.

**8. Risk.** Feasibility on patches may be `0` at *every* size (K6 tested up to ~800 faces), in which
case the curve is flat at zero and only the notch-share half of the claim survives — an accounting
statement without a consequence.

**9. Effort.** Derivation **S**. Code **S**.

---

## B7 — Universality: the deployment Jacobian of a random periodic graph concentrates

**1. Title / type.** *`θ_c → θ*` for random Voronoi tori* — **characterization** (empirical law with
a homogenisation argument).

**2. Claim.** The pattern-level harmonic `(a, b)` of a random periodic graph concentrates as the cell
grows: on the 12 Voronoi tori of `k7_main.csv`, `θ_c = 2 atan2(a, b)` is `2.100, 2.059, 1.806, 1.747,
1.988` for `|F_cell| = 20 … 55` and `1.952, 1.944, 1.925, 1.940, 1.860, 1.933, 1.931` for
`|F_cell| = 70 … 200` — mean `1.926`, range `0.092` on the large half against `0.353` on the small
half. Moreover `a + b = 88.8 ± 1.3` on **all** rows with `|F_cell| ≥ 55` (values `90.4, 89.8, 89.3,
88.8, 89.1, 87.8, 89.0, 88.8`) while `a` and `b` separately vary by 8%. Claim: `K` for a random
periodic planar graph with max-cut `σ` converges in probability to a deterministic `K*` depending only
on the point-process family, and the `a + b` invariance is an exact identity, not a coincidence.

**3. Gap targeted.** 2026 §5.1 designs `J(θ)` per pattern and reports it per example; it has no notion
of a *generic* value, and its random-graph experiments are all on fixed-boundary patches. 2025 §13–15's
"automatic tiling selection from curvature" presumes patterns differ; a universality law says random
graphs *do not*, which is a genuinely non-guessable outcome and constrains inverse's I5 (pattern
selection) from below.

**4. Novelty vs field.** (a) `field_rigidity` **R11** Guest–Hutchinson analyse infinite repetitive
structures but give no law of large numbers over random tilings. (b) **R16/R17** Borcea–Streinu are
per-framework. (c) `field_kirigami` **R2** Mitschke et al. enumerate 1- and 2-uniform tessellations —
a catalogue, not a limit theorem. I could not find any statement of a homogenised deployment Jacobian
for random kirigami; marked **[SCREEN-ME]** since I ran no searches this session.

**5. Why it might be true / sketch.** `a = ½ det(P₀) tr K` and (B.3) express `a` as a sum of `|E_hinge|`
local terms `σ_g ⟨d_e, x_src⟩` plus `|E_split|` terms `r_e`; for a stationary ergodic point process
these are sums of local functionals of the Voronoi tessellation, so a law of large numbers applies
provided the max-cut `σ` is itself asymptotically local — which is the questionable step, since Eq. (1)
is a global relaxation. `a + b = ½ det(P₀)(1 + tr K − det K)` and the observed invariance says
`tr K − det K` is pinned; I have no derivation of that and would want it.

**6. Kill experiment.** `kill_b7.cpp` = `kill_k7.cpp`'s C4 block over 60 Voronoi tori (5 seeds ×
`|F_cell| ∈ {50, 100, 200, 400, 800, 1600}`) plus 20 Delaunay tori: report `a, b, θ_c, a+b, det P₀`.
**Kills the idea:** the standard deviation of `θ_c` within a size does not shrink with `|F_cell|`
(no `~|F|^{-1/2}` trend over a 32× range), or `a + b` varies by more than a few percent once `det P₀`
is normalised out. ≈ 25 min (K7's per-torus cost at `n = 200` is seconds).

**7. If it survives, the demo.** `θ_c` vs `|F_cell|` with error bars over 80 tori and the fitted
constant; the accompanying sentence is "an arbitrary random periodic graph has an essentially
determined deployment Jacobian, so `σ` and the graph — not the design freedom — fix the mechanics".

**8. Risk.** `n = 12` samples, one per size, no repeats: the concentration may be an artefact of a
single seed sequence, and `det P₀` may not be constant across the 12 rows (I did not verify that it is,
which is why the kill test normalises by it). This is the least theorem-like idea in the file and I
rank it accordingly.

**9. Effort.** Derivation **M** (the LLN step is real work). Code **S**.

---

## B8 — Planar duality flips the sign of `b`, and `θ_c(P) + θ_c(P*) = 2π`

**1. Title / type.** *The dual-pattern relation for the second closed angle* — **theorem**.

**2. Claim.** For the triangle/hexagon dual pair in `k7_main.csv`: `a = 6.000000` for both,
`b = −3.464102` (triangles) and `+3.464102` (hexagons) — equal and opposite to 7 digits — and
`θ_c = 4π/3` and `2π/3`, summing to `2π` **exactly**. Claim: for a pattern `P` and its planar dual
`P*` carried by the dual `σ`, `a(P*) = a(P)` and `b(P*) = −b(P)` after normalising `det P₀`, hence
`θ_c(P) + θ_c(P*) = 2π`; equivalently `tr K* = tr K` and `det K* = 2 − det K`.

**3. Gap targeted.** 2026 §5.1 and §4.4 treat each pattern separately; no relation between patterns
appears in either paper. 2025's 16-pattern library (`baseline/kirigami_tessellations`, F31) contains
several dual pairs, so the claim is immediately testable on published data.

**4. Novelty vs field.** (a) `field_kirigami` **R1** Grima–Evans rotating squares/triangles are the
canonical special cases and are *not* stated as a dual pair with a conserved quantity. (b) **R2**
Mitschke et al. catalogue tessellation families without such a relation. (c) The Screen-Scout verdict
**S7** ("2nd closed angle NOT FOUND") applies a fortiori to a duality law about it. **[SCREEN-ME]** —
this is the idea in the file most likely to already exist in the auxetics literature under
"rotating squares / rotating triangles duality", and it must be screened before writing.

**5. Why it might be true / sketch.** `a` is the `sin θ` coefficient, an *odd* harmonic; `b` is the
`(1 − cos θ)` coefficient, an *even* one. Passing to the dual exchanges the roles of faces and
vertices, which (in the pure-hinge case) exchanges the two families of the bipartite `Γ`, i.e. flips
`σ` globally. By `derivations/core.md` 0.8 a global `σ` flip is `θ → −θ` **plus** a reversal of every
hinge direction; `θ → −θ` alone sends `(a, b) → (−a, b)`, so the composite acting as `(a, b) →
(a, −b)` is exactly the statement that the direction reversal contributes the extra sign. That is the
step I have **not** proved, and it is the whole content of the theorem.

**6. Kill experiment.** `kill_b8.cpp`: the seven tiling families already in `make_tiling_pattern`
plus the dual of each (constructed by the existing mesh dual, or loaded from
`baseline/kirigami_tessellations`), computing `(a, b, θ_c, det P₀)` for each. **Kills the idea:** any
dual pair with `|θ_c(P) + θ_c(P*) − 2π| > 1e−6` after normalisation, or `b(P*) ≠ −b(P)`. ≈ 10 min;
the squares family is self-dual and is the first sanity case (`b = 0`, `θ_c = π`, `π + π = 2π` ✓ — a
non-trivial consistency check the claim already passes).

**7. If it survives, the demo.** A two-column figure of dual pairs from the 2025 library with their
`θ_c` and the conserved sum; and the corollary that a pattern with `θ_c > π` (triangles) is dual to
one with `θ_c < π`, which is a design rule for choosing between a pattern and its dual.

**8. Risk.** The evidence is **one** dual pair at three cell sizes. Kagome is self-dual-ish and has
`b = 0`, `θ_c = π`, consistent; but three families is not a theorem, and the `σ`-flip step in the
sketch may simply be false for patterns with split cuts (hexagons has `n_split = 4`, triangles `0`,
so the pair I am reasoning from is not even split-symmetric).

**9. Effort.** Derivation **M**. Code **S**.

---

## B9 — The budget law holds for every flex, not only the uniform one

**1. Title / type.** *A discrete divergence theorem for the whole mechanism* — **theorem**, and a
bound on adversary Idea 1's expansive cone.

**2. Claim.** Let `V ∈ ker A(X, 0)` be **any** infinitesimal flex of the body-and-pin framework (not
necessarily the uniform ray `σ`; F26/K3a measured `dim ker A ≈ 305` after 2-core on Delaunay patches,
so this cone is huge). Then the total first-order opening rate of all preimages under `V` is again a
**border** sum:

```
   Σ_C  Ȧ_C[V]  =  ½ Σ_{border (a→b) ⊂ f}  det( x_a − x_b , V_f(x_·) ) -type boundary term ,
```

i.e. the same telescoping as (0.c) with `2 s J u_f` replaced by the flex's per-face velocity field.
Corollary: adversary Idea 1's cone `P(X) = { V : every cut opens }` is contained in the half-space
`{ V : boundary term > 0 }`; on a patch whose deployed border is constrained to move rigidly, that
term vanishes and **`P(X) = ∅` — no flex whatsoever opens every cut.**

**3. Gap targeted.** 2026 §7 limitation (ii) verbatim in F8: non-uniform hinge-angle deployment and
the DOF of the deployment space, "could be achieved by sensitivity analysis". Adversary Idea 1
(`ideas/round2_adversary.md`) formulates that space as an LP cone and asks whether `σ` points the
wrong way inside it; this idea supplies the *linear functional that every ray of the cone must
satisfy*, which is a structural fact their LP does not contain, and which turns their "the cone is
large so `σ` is just unlucky" reading into a testable dichotomy.

**4. Novelty vs field.** (a) `field_rigidity` **R8** Connelly–Whiteley and the expansive-motion
literature (Connelly–Demaine–Rote) work with pairwise distance monotonicity; the functional here is an
*area* functional and is an equality, not an inequality. (b) **R11** Guest–Hutchinson give periodic
mode counts, no area law. (c) **R26** Chen–Choi–Mahadevan compute DOF for corner-pinned quads
numerically. **[SCREEN-ME]** on the expansive-motions line specifically, which is also flagged in
adversary Idea 1.

**5. Why it might be true / sketch.** (0.a)–(0.c) never used that the motion is the uniform one until
the substitution `t_f = 2 s J u_f`. For a general flex, each face `f` has velocity
`v_f(x) = ω_f J (x − c_f) + τ_f`, and the same two steps go through with `⟨x_a − x_b, ·⟩` replaced by
the corresponding bilinear pairing; the per-face telescoping still holds because a face polygon is
closed and `v_f` is affine in `x`. Both the derivation and the exact form of the border term need to be
written out — I have *not* done that here, and the coefficient may not be as clean as (B.4).

**6. Kill experiment.** `kill_b9.cpp` on 20 patches: build `A(X, 0)` with the existing
`method/mobility.hpp`, draw 50 random `V ∈ ker A` per patch, compute the interior sum
`Σ_C Ȧ_C[V]` by finite-differencing each hole polygon's area along `V`, and compare with the border
sum. **Kills the idea:** relative deviation `> 1e−8` on any flex. Then, the useful half: for each
patch report `max_{V ∈ ker A} min_e (cut opening rate)` by LP — if that is `> 0` on patches where the
uniform ray fails, adversary Idea 1's reading is right and the mechanism is that `σ` is unlucky; if it
is `≤ 0` on all 20, the emptiness is structural and B4's boundary explanation is the operative one.
≈ 25 min (`dim ker A` up to 1675 on Delaunay, so cap `|F| ≤ 500`).

**7. If it survives, the demo.** One plot: for 20 patches, the uniform ray's `min_e q_e` against the
LP optimum over the whole cone. It settles, in one figure, whether non-uniform deployment can rescue
what uniform deployment cannot — the question 2026 explicitly leaves open.

**8. Risk.** The general-flex telescoping may leave an interior remainder (the flex is not a rigid
motion per face plus a *global* potential, and the hinge-point identification that gives (B.2h) its
clean form uses `u`, which only exists on the uniform branch). If so, the theorem degrades to "there
is a border term **plus** an interior defect", and only the periodic corollary survives.

**9. Effort.** Derivation **M/L** (this is the one derivation in the file I have not completed).
Code **M**.

---

## B10 — Per-hole closed forms: max-area angle is exactly half the re-closing angle

**1. Title / type.** *`θ_max-area(C) = ½ θ_c(C)`, per hole, and the exact condition for a global
second closed state* — **characterization**.

**2. Claim.** From (B.1) in the `τ = tan(θ/2)` chart, `(1 + τ²) A_C = 2τ(a_C − b_C τ)`, so every hole
re-closes at `τ = a_C/b_C`, i.e. `θ_c(C) = 2 arctan(a_C/b_C)` (`= 2 atan2(a_C, b_C)` on the correct
branch), and `A_C′(θ) = a_C cos θ − b_C sin θ = 0` at `tan θ = a_C/b_C`, i.e. at **exactly half** of
it. Consequences: (i) `b_C ≤ 0` ⟹ the hole never re-closes on `(0, π)` and opens monotonically;
(ii) a pattern has a *global* second closed (gapless) state iff the ratio `a_C/b_C` is the **same for
every hole**, which is a codimension-`H − 1` condition satisfied by the symmetric tilings and by
nothing else generically; (iii) the periodic special case `θ_max-area = ½ θ_c` is exactly
`ideas/round2_inverse.md` I4, which is therefore a one-line corollary of (B.1) rather than a
conjecture needing a homogenisation ansatz — and it holds per hole, on arbitrary graphs, with no
periodicity.

**3. Gap targeted.** 2025 Observation 3.1 (§3.2) as quoted in `ideas/round2_inverse.md` I4: a pattern
"typically" admits two distinct maximally-open configurations. "Typically" is exactly the
codimension-`H − 1` condition in (ii), which is the content here.

**4. Novelty vs field.** (a) Screen-Scout **S7**: hole-area harmonic and second closed angle **NOT
FOUND** in the field. (b) `ideas/persona_optimizer.md` (round 1, per `STATE.md` open-task line) already
has `A(θ) = q(cos θ − 1) + r sin θ` and `2 arctan(r/q)` for the *pattern*; `periodic_jacobian.hpp` has
it for the periodic cell. Neither has the per-hole version, the `½` relation, or the constancy
criterion. (c) Inverse I4 claims the `½` relation for the periodic total and calls it a claim to be
tested; here it is a corollary. **I am building on I4, not duplicating it, and I say so.**

**5. Why it might be true / sketch.** Pure algebra of `h(θ) = a sin θ − b(1 − cos θ)`: written as
`R sin(θ + φ)` shifted, its unique interior maximum and its unique nonzero root are in the ratio
`1 : 2` because the harmonic is *pinned at zero at `θ = 0`* — the `p + q = 0` degeneracy of T3.H.5.
So the `½` is not a coincidence of the periodic case; it is the signature of the whole degenerate
class. Checked by hand on `squares_2x2` (`b = 0`: root at `π`, max at `π/2` ✓, and `A = 4 sin θ` does
peak at `π/2`) and on `hexagons_2x2` (`a = 6, b = 2√3`: `θ_c = 2π/3`, max at `π/3`).

**6. Kill experiment.** *Partly run (§0.g): the numerically located area maximum matches `½ θ_c` to
`6.4e−4`, which is half the scan's grid step, i.e. consistent with exact.* Folded into B1's app: for every hole of every corpus pattern, compare the
fitted `θ_c(C)` and the numerically located area maximum. **Kills the idea:** any hole where the ratio
differs from `2` by more than `1e−8` (a fit-precision kill, since the relation is exact if (B.1) is).
For (ii): report `spread of a_C/b_C over holes` per pattern and check it is `~1e−15` on the symmetric
tilings and `O(1)` on random graphs; if a random graph shows a constant ratio, (ii)'s genericity claim
is wrong. ≈ included in B1's minutes.

**7. If it survives, the demo.** A per-hole `θ_c` colour map on a random pattern next to a tiling: the
tiling is uniform, the random graph is a rainbow — one picture explaining why 2025's gapless second
state is a *tiling* phenomenon. Plus the corrected statement of 2025 Observation 3.1.

**8. Risk.** It is nearly free algebra once (B.1) is in hand, so its novelty rests entirely on the
per-hole/arbitrary-graph generality and on I4 not having got there first. If the Critic decides I4
and this are one idea, this collapses into B1 as a corollary — which is where it belongs anyway.

**9. Effort.** Derivation **S** (done). Code **S** (shares B1's app).

---

# 2. Where I disagree with the three sibling round-2 files

I cite them above where I build on them; here is where I think they are wrong, each with the test
that settles it.

**(a) `ideas/round2_theorist.md` §0.c and Idea 8 ("the barycentric row is precisely what prevents the
`0⁺` collapse").** The argument is: `(BAL) Σ_{u ∈ In(v)}(x_u − x_v) = 0` makes the `k_v` separation
velocities sum to zero, hence they positively span the plane, hence no copy can be pushed into its
neighbours' material. The second "hence" does not follow — vectors summing to zero positively span the
plane, but that says nothing about whether a *particular* copy's velocity lies inside a *particular*
neighbouring corner cone, which is exactly the predicate `zero_plus.hpp` calls `mu` and computes as
`max(−g₁, −g₂)`. The claim also has to survive a measured fact: every K6 design satisfies Eq. (2)
exactly, so (BAL) holds at every pure interior vertex of all 400, and all 400 still have `Θ_max = 0`.
**Test that settles it (5 minutes, no new code):** in the K6 population, restrict
`corner_incidences` to vertices touched by **no** split edge and report how many of those have
`mu ≤ 0`. Theorist-a's Idea 8 predicts zero. If the count is positive, Idea 8 is dead as stated.

**(b) `ideas/round2_adversary.md` Idea 5 ("the hinge edges of a hole always contribute positively").**
Their decomposition assigns each hinge edge `w_e = |det(x_tgt − x_src, other hinge arm at src)| > 0`
by construction. My (B.2h) gives the same total as `½ σ_{g(e)} ⟨d_e, x_{src(e)}⟩`, whose individual
terms are **origin-dependent and not sign-definite**; two decompositions of one origin-independent
quantity can differ by terms summing to zero, so both can be correct, but a decomposition into
manifestly positive terms is a strictly stronger claim and it implies that **every split-free hole has
`a_C > 0`**. That is falsifiable: `k7_main.csv` already shows `b` (the other coefficient) taking both
signs across split-free patterns (`squares` and `kagome` at `0`, `triangles` at `−3.46`), so sign
definiteness is not a general feature of these coefficients. **Test:** B1's app reports `a_C` per hole;
one split-free hole with `a_C ≤ 0` kills their positivity claim and with it the "the hinge edges pay,
the split edges spend" reading of their necessary condition.

**(c) `ideas/round2_inverse.md` I4 ("max-area configuration sits at exactly half the second closed
angle") is not a claim to be tested; it is a corollary.** B10 derives it in two lines from the
`p + q = 0` degeneracy, per hole, on arbitrary graphs, with no homogenisation ansatz and no
periodicity. Its I2 areal budget `λ_A(θ, K) = det J(θ)` is (B.5) rewritten. Both should be re-scoped
as consequences of one identity rather than presented as two independent findings, or a referee will
say the paper has one idea dressed as three.

**(d) All three files treat the fixed boundary as given.** If B4's experiment comes back positive, the
central negative result of this project (F30, `0/400`) is a statement about 2026 Eq. (4) and not about
random graphs, and every emptiness-certificate idea in the sibling files (theorist-a Idea 2, adversary
Idea 2) is certifying the emptiness of a set the paper chose rather than one the problem forces. That
does not make those ideas wrong — it changes what they mean, and the change is large enough that
**B4 should be run before either of them is committed to.**

---

# 3. Self-attack on my top three, as a hostile TOG reviewer

**Against B3 (the budget half-space in `K`-space).** "The authors observe that a linear functional of
the design variable controls a first-order area rate and christen it a *budget*. Stripped of the
vocabulary, the content is: `det J(θ) − 1` is the areal expansion, its derivative at zero is
`½ tr K`, and if you ask a pattern to expand less it opens its holes less. Equation (11)–(13) of the
2026 paper already parametrises `J(θ)`; the observation that `tr K` is its area derivative is one line
of calculus that any reader would supply. The paper's actual contribution is then a *threshold*
`τ* = W/det P₀`, and the authors admit in their own risk paragraph that `W` is not constant along the
achievable set, so `tr K ≥ τ*` is a linearisation of unknown quality. What they have measured is that
one hand-picked anisotropic target fails 25 times and random targets succeed 12 times out of 25 — a
correlation with `n = 50` and no controlled variation of `tr K` at all. Reject unless the constrained
re-solve raises the certified rate on the *same* targets."
— My answer: the reviewer is right that the identity is elementary once seen and that the retrodiction
is a correlation. The paper stands or falls on kill step (ii), the constrained re-solve on identical
targets, which is a controlled experiment; I would not submit without it, and if it fails I would
publish B1 + B4 and drop B3.

**Against B1 (the budget identity).** "This is Green's theorem. The deployed sheet is faces of constant
area plus holes; the total enclosed area is a contour integral over the outer boundary; therefore the
total hole area is a boundary functional. The authors have rediscovered the divergence theorem in a
discrete setting and decorated it with the face potential `u`. Worse, it is *vacuous* on a patch,
because the boundary term is itself a free quadratic form on the shape space: the identity constrains
nothing until you fix the boundary, and once you fix the boundary the interesting statement is about
2026's Eq. (4), not about kirigami."
— My answer: the geometric content *is* the divergence theorem, and the file should say so plainly
rather than pretend otherwise. What is not the divergence theorem is (B.3), the split-into-hinge-and-
split-edge decomposition, which is what connects the global law to the exact quantity `q_e` that all
of K5/K6 measures, and (B.5)'s identification of the boundary term with `tr K`. The reviewer's second
point is correct and is why B2 and B4, not B1, carry the weight.

**Against B4 (the boundary is the obstruction).** "The claim is that removing the boundary constraint
makes the design space non-empty. But the boundary constraint is not decoration: without it the Tutte
projection Eq. (6) has nothing to hold it, and the authors' own F17 reports 4–6% inverted faces
*with* the boundary pinned. A free-boundary design that is `0⁺`-feasible and self-intersecting in the
flat state is not a design. The paper would then be claiming a design space that is non-empty only in
a set the authors themselves certify as invalid. Furthermore the periodic evidence offered — 11 of 12
random Voronoi tori zero-plus feasible — has no fixed-boundary control on the *same* graphs, so the
comparison is between two different generators."
— My answer: the second point is fair and cheap to fix, and I would fix it: generate each Voronoi
torus, then cut a disk-topology patch from the *same* point set and run both. The first point is the
real risk and it is why B4's kill criterion must be the full certificate (`POS ∧ NOOVERLAP ∧ NOROOT`),
never `0⁺` feasibility alone.

---

# 4. Ranked list

| # | idea | type | why here | kill test | effort |
|---|---|---|---|---|---|
| 1 | **B3** `tr K` is the budget; budget half-space in `K`-space | thm + alg | the only pair in this file that is a *constructive* deliverable: certified deployable designs on arbitrary **periodic** graphs, where 2026 shows only authored tilings, and it retro-explains K7 C3's 0/25 vs 12/25 split | `kill_b3`: certified rows must separate by `tr K − τ*`; constrained re-solve must beat 12/50 on the same targets | S / M |
| 2 | **B4** the fixed boundary, not the graph, empties the space | thm + exp | cheapest decisive experiment in the project; if positive it rewrites the meaning of F30 and of two sibling ideas | `kill_b4`: `0/400` again ⟹ dead, ≈ 20 min | S / S |
| 3 | **B1** the budget identity `Σ_C a_C = B(X)` | theorem | the identity everything else is derived from; the per-hole half is **already verified** at `8.7e−14` over 471 holes (§0.g) and it reproduces `periodic_jacobian.hpp`'s independent coefficients | `kill_b1`: only the global sum (B.4) is left; fitted vs closed-form per hole already passed | S / S |
| 4 | **B2** `O(\|E\|)` budget bound on `min_e q_e` | thm + alg | turns the identity into a microsecond non-deployability certificate; genuinely two-sided, since it may be true and never binding | `kill_b2`: violation ⟹ hard kill; never binding ⟹ demoted to a corollary | S / S |
| 5 | **B10** per-hole `θ_max-area = ½ θ_c`, and the constant-ratio criterion | charac. | free algebra, sharp statement of 2025 Observation 3.1, and it subsumes a sibling's headline claim | shares B1's app; ratio ≠ 2 ⟹ dead | S / S |
| 6 | **B9** the budget law for every flex; bound on the expansive cone | theorem | the right way to settle whether non-uniform deployment can rescue uniform failure — 2026 §7 limitation (ii) — and the LP half is decisive either way | `kill_b9`: interior sum vs border sum to `1e−8`; then the cone LP | M-L / M |
| 7 | **B6** notch leakage and the `\|F\|/\|∂M\|` scaling law | thm + law | correct accounting and a plausible mechanism for the patch-vs-torus gap, but likely to produce a curve that is flat at zero | `kill_b6`: no monotone trend ⟹ dead | S / S |
| 8 | **B5** spectral budget allocation | algorithm | pairs with adversary Idea 2 into "a design or a proof there is none"; but if B4 shows the fixed-boundary set really is empty, this returns certificates and no designs | `kill_b5`: `0/400` again | S / M |
| 9 | **B8** duality: `θ_c(P) + θ_c(P*) = 2π` | theorem | exact on the one dual pair measured and self-consistent on the self-dual case, but three families is not a theorem and it is the most likely to be known already | `kill_b8`: any pair off by `> 1e−6` | M / S |
| 10 | **B7** universality of `K` for random periodic graphs | charac. | the most surprising thing in the K7 data (`θ_c = 1.926 ± 0.05` over a 3× size range) and the least theorem-like; `n = 12`, one seed per size | `kill_b7`: no `\|F\|^{-1/2}` narrowing over 32× ⟹ dead | M / S |

**If only one thing is run:** B4, then B1, then B3. B4 is twenty minutes and can invalidate the
framing of three files; B1 is the identity they all then need; B3 is the paper.

---

# 5. What I could not verify

* **(B.1) and (B.2s) are measured** (§0.g, `derivations/scratch/check_b1.cpp`, 471 holes / 58 split
  edges, worst `8.7e−14`). **(B.3)'s hinge term and (B.4) are not** — the geometric tracer returns
  holes but not notches, so the global sum stays unclosed, and every claim that rests on `B(X)`
  (B2, B4, B6) rests on an unverified identity. One of my hand corollaries was already refuted by the
  K7 data (§0.f, last bullet); treat the unmeasured parts with the same suspicion.
* **The global sign in (B.4)** is not determined here; it depends on the face-boundary orientation
  convention and on whether hole walks are traversed with the hole on the left. Every downstream
  statement is sign-symmetric except B2's inequality direction, which must be pinned numerically
  before it is quoted.
* **`W(X)` constant along `𝒦`** (B3's threshold) is assumed at first order and is false at second
  order. Unmeasured.
* **B9's general-flex telescoping is not done.** I asserted the shape of the result, not the result.
* **No literature search was run this session.** Three novelty rows are marked **[SCREEN-ME]** (B7's
  homogenised Jacobian, B8's duality law, B9's expansive-motions overlap). B8 in particular smells
  like something the rotating-squares/rotating-triangles literature may already contain.
* **`det P₀` across the 12 Voronoi tori** is not verified to be constant, which B7's `a + b ≈ 88.8`
  observation silently assumes.
* **Paper citations** are via F1–F31, `code/README.md` and the sibling files; I did not open
  `papers/*.txt` or the PDFs. The 2025 Observation 3.1 wording in B10 is quoted from
  `ideas/round2_inverse.md` I4, not from the paper.
