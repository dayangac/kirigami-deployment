# Critic — master ranking, hostile review, and kill-experiment specs

Written 2026-09-03 by the Critic subagent. Inputs actually read in full or in the stated part:
`specs/critic.md`, `specs/common_preamble.md`, `STATE.md` (F1–F19, U1–U9),
`ideas/persona_geometer.md` (806 lines, all ten ideas + Sec. 0 + self-attack),
`ideas/persona_rigidity.md` (678 lines, all ten + Sec. 0 numerics),
`ideas/persona_optimizer.md` (825 lines, ideas 1–5, 8, 10, overlap section, ranked list),
`ideas/persona_adversary.md` (1222 lines, in full: Part 0 search log, Part 1 attack, Parts 2–6),
`notes/paper_2026.md` (Prop. 4.1 block at :508–:529, Sec. 4.4 rank claims :671, Eq. (9) :734,
limitations :951–:965, assumptions, open questions),
`notes/paper_2025.md` (Sec. 7.1/7.2 — Fig. 8 and **Fig. 10**, the implication lattice at :576–:585),
`notes/field_tutte.md` (rows 1–18, routes A/B/C, summary table),
`notes/field_kirigami.md` + `notes/field_rigidity.md` (row indices and Sec. 2 verdicts, as cited by
the personas; I did not re-read every row),
`notes/screen_bundle.md` (S1 and S4 complete; S2, S3, S5–S7 still `_pending_`),
`results/core_validation/rank_claim.md`, `results/core_validation/reference_cases.md`,
`code/README.md` (full: API, CLI tools, and the ten documented deviations).

**Honesty policy.** I ran no code and no web searches. Every number I quote is from `STATE.md`
F15–F19 or from `results/core_validation/`. Where I overturn a persona's verdict I say which
document overturned it. Where a novelty claim rests only on the field tables I say so in Sec. 6.

---

## 0. Three facts from Phase 2 that reorder everything the personas wrote

The four persona files were written **before** the Builder-Core validation landed. Three measured
facts change the ranking, and no persona had them.

**(P1) The design space and the collision problem have the same source: split cuts.**
F15: `dim_null = |E_split|` on 100/100 random graphs (fixed boundary). F18: `theta_max = min_i beta_i`
exactly on **every** split-free pattern, and strictly smaller as soon as split cuts exist
(snub square 1.646 vs 3.384; `reference_cases.md`). Put together:

> A pattern has a nontrivial shape space **iff** it has split cuts, and it has a non-local collision
> deficit **iff** it has split cuts.

So every range-optimisation idea (G10 / rigidity-9 / optimizer-4 / A1) is *automatically* aimed at
exactly the patterns where it has room to act, and every "the shape space may be empty so there is
nothing to optimise" risk paragraph in those four files is wrong in the regime that matters. This is
the single strongest structural fact available to this project and none of the four files states it.
It also bounds the ambition: on the pure hinge case there is **nothing to do** — the 2025 local
formula `min_v (2pi - alpha_i - alpha_j)` is already exact, measured on 3/3 split-free reference
cases. Any paper built on exact contact analysis must lead with split-cut patterns or it is empty.

**(P2) The paper's own pipeline returns a geometrically invalid embedding on 57/57 random graphs.**
F17: `X0` from Eq. (6) has inverted faces on 57/57 graphs where it was computed (0/27 on the
re-run sweep, see P4), mean 5.9% of faces inverted, max vertex displacement 11.3 in a box of side 40. The *input* `X_ini` is a valid planar
embedding by construction. So Eq. (6) takes a valid embedding and returns an invalid one, on every
random graph tried, while the eight authored tilings in `reference_cases.md` are all fine. The honest
one-sentence framing: **the method works on the tilings the paper shows and fails on the arbitrary
planar graphs its title claims.** No persona had this; the geometer *predicted* it (G8) and was right.

**(P3) Eq. (9) is not merely a first-order heuristic — it fails outright on named tilings.**
F18: on 4.8.8, `gamma = 1e-4` raises `theta_max` from 0.059 to 2.027 but `gamma = 1e-3` gives 0; on
the snub square **no gamma in the ladder `{1e-4 ... 10}` helps at all** (`reference_cases.md`:
`theta_max` stays 1.646 against `min beta = 3.384`, with 13 dimensions of shape space available).
The adversary's A1 asks "what is the false-negative rate of Eq. (9)?" — a subtler question than the
one the data already answers, which is "on a named Archimedean tiling with a 13-dimensional design
space, the published optimiser finds nothing." That is a better hero case than any of A1's histograms
and it is already reproducible in `kiri_reference`.

**(P4) Addendum, added after the orchestrator re-ran the sweep and added three checks to
`results/core_validation/rank_claim.md`.** The sweep is now 50 graphs (115–2734 faces), of which 27
had the dense null space computed, and `X0` is invalid on **0/27** — the same conclusion as F17's
0/57 on the earlier 100-graph sweep, at a lower N. `STATE.md` F17 still quotes 57/57; the artefact on
disk now says 0/27. **Both numbers should be reported, or the sweep re-run at n=200, before either is
put in a paper.** Four further measurements change items in my table:

- `L = R·D` exactly, 50/50, and `rank(L) = H - dim Z` with `Z` out-harmonic, 27/27. **F14 is now
  measured, not just derived.** This confirms the easy direction of B5/G5/A7 and removes it from the
  "unverified" column.
- On boundary-free (torus) patches `1^T L = 0` exactly and **`rank(L) = H - 1`, measured** on the 4x4
  and 6x4 square tori and the 4x4 triangle torus. So 2026 §4.4's "number of independent equations
  equals the number of holes" is now *confirmed false as stated* for periodic patterns, not merely
  predicted. **My K3b below is therefore already half-done**; what remains is to widen it from 3 tori
  to a population and to check whether the extra kernel direction is a usable design degree of
  freedom or only the global translation.
- `H = |E_hinge| - |F| + c(Gamma)` holds **50/50** and `c(Gamma) = 1` on 50/50 under the Eq. (1)
  orientation. So U1 and B4 are **confirmed as stated**, not corrected.
- **An internal inconsistency in that file which someone must resolve.** Its prose says "Check 3, as
  first stated, is false whenever notches exist" and that the empirically holding identity is
  `H_all == |E_hinge| - |F| + c(Gamma)` counting holes *plus* notches. Its own tables say the
  opposite: check 3a (`H == pred`) holds 50/50 and check 3b (`H_all == pred`) holds **0/50**, and
  every reference-case row reads `H = pred: yes`, `H_all = pred: no`. The tables are right and the
  prose paragraph is inverted: the bounded faces of a plane graph number exactly `b_1 = |E| - |V| +
  c`, a notch opens through the patch boundary and therefore merges into `Gamma`'s *outer* face, so
  notches are not bounded faces of `Gamma`. The prose should be corrected before anything cites it.

**One correction to the adversary, from the predecessor paper.** A8 ("Eq. (2) is sufficient but not
necessary") claims as its contribution the *existence* of patterns that deploy non-uniformly but not
uniformly. That existence is **already published, with a photograph**: 2025 Fig. 10, caption verbatim
(`notes/paper_2025.md:552`): "A tiling containing deployment-unfriendly vertices (highlighted in red)
**can still be deployable** (b,c). However, as shown in (c), the **opening angles are unequal**...
In (b) and (d) we show a simulated and a fabricated result." The 2025 implication lattice
(`notes/paper_2025.md:576-585`) states it as a row. A8's existence half is therefore dead on arrival;
only the *characterisation* half survives, and it is much harder. This alone moves A8 from the
adversary's "promote to rank 2" to the bottom half of my table.

---

## 1. Master table of all ideas (deduplicated)

Merged across the four files. `Personas` lists every file that proposed the same claim; convergence
is recorded, not rewarded. C1–C7 are the MISSION §3/§8 constraints, scored 0/1.
Effort: S/M/L for derivation + code.

| id | title | type | personas | C1 | C2 | C3 | C4 | C5 | C6 | C7 | kill experiment (one line) | effort | verdict |
|---|---|---|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|---|---|---|
| B1 | Deployment path is a conic: `Y_th = cos(th/2)C(X) + sin(th/2)S(X)`, `C,S` linear in `X` | thm | G1, rig-3, opt-(i), adv-(a) | 1 | 1 | 0 | 1 | 1 | **0** | 1 | already done: fit `C,S` at 2 angles, predict 5 more (rig V1: 6e-15 on 21/21) | S+S | **PARAGRAPH.** Lemma 1. Half-angle is RUM folklore (screen S1: Acuna 2022 has counter-rotation on arbitrary bipartite networks). Linearity in `X` is new but parasitic on Segall's own Eq. (5). |
| B2 | Infinitesimal conformality ⇒ exact conformality; `J(th)=cos(th/2)I+sin(th/2)W` | thm | G2, rig-4, opt-1, adv-(c) | 1 | 1 | 1 | 1 | 1 | 0 | 1 | 30 periodic patterns, solve Eq. (13), measure `\|J11-J22\|+\|J12+J21\|` at 200 angles; fail if > 1e-9 | S+S | **PARAGRAPH.** Screen S4 verdict NOT FOUND, so it is new; it is still one line of linear algebra proving a quoted empirical sentence. Remark, not section. |
| B3 | Exact contact angles: every orientation predicate is `p + q cos th + r sin th` | char | G3, rig-3b, opt-3, adv-(b)/A1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | closed-form `theta_max` vs `collision.hpp` bisection on 300 graphs to 1e-6 | S+M | **PURSUE, as machinery.** Must be restated: collinearity ≠ contact (adv Objection 1); the repaired form is "smallest positive root of a quadratic in `t`, filtered by two closed-form interval tests". |
| B4 | `H = \|E_hinge\| - \|F\| + c(Gamma)`; holes = bounded faces of Gamma | thm | G4, adv-(d), U1 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | Alg. 1 count vs Euler formula on 500 random-sigma graphs | S+S | **PARAGRAPH.** Euler. Now confirmed 50/50 as stated, with `c(Gamma)=1` 50/50 (P4); `field_tutte` (2.1) is the `c(Gamma)=1` special case. |
| B5 | `corank(L)` = number of closed classes of the hole digraph | thm | G5, A7 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | `rank(L)` by ColPivQR vs Tarjan SCC closed-class count, 500 graphs random sigma; plus corank>=2 search | M+M | **PURSUE-IF-KILL-PASSES.** Corank 0 with fixed boundary (F15), and the periodic corank-1 failure is now measured on 3 tori (P4). The theorem's easy direction (`rank = H - dim Z`, `Z` out-harmonic) is confirmed 27/27. Only the `corank >= 2` search is still open. |
| B6 | `E_split` is a T-join; planar max-cut is polynomial; `dim X >= \|T\| + 2 dim Z` | thm+alg | G6, opt-6 | 1 | 1 | 0 | 1 | 1 | 0 | 1 | parity `deg_split(v) = deg(v) mod 2` on 500 odd-degree-heavy graphs; blossom vs Eq. (1) gap | M+M | **DROP the algorithm** (Hadlock 1975; the geometer concedes this). **Keep the bound** as a lemma: odd valency forces design freedom. |
| B7 | `L = A_S · Delta_h` (split-component-averaged hinge Laplacian) | char | G7 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | assemble `L` both ways, require exact integer equality on 500 graphs | S+S | **PARAGRAPH (machinery).** Refines F14's `L = R·D` to vertex level. Engine for B5, B8. |
| B8 | Valid (injective) Tutte auxetic embeddings: theorem/algorithm + certificate | thm+alg | G8, orchestrator candidate, F17 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 200 graphs: is `X0` invalid? does any of 1e4 sampled `X in X` have all signed areas > 0? | M+M | **PURSUE — top candidate.** Baseline fails 57/57 (F17). Both outcomes publishable (algorithm, or emptiness). See Sec. 3. |
| B9 | The non-uniform deployment variety (real / complex-torus / integer-ratio forms) | char | G9, rig-6, opt-8, A8 | 1 | 1 | 0 | 1 | 1 | **0** | 1 | `dim Theta = m`; `ker_C W = span{1,sigma}`? Newton from 200 torus starts on `dim_C ker W >= 2` | M+M | **DROP as framed.** 2025 Fig. 10 already exhibits and fabricates a non-uniformly-deployed pattern, and 2026 §6 **names sensitivity analysis as the method**. kiri-R17 has a "deployment angle field" for quads. Adversary's own O1 rules this out; A8 does not escape it. |
| B10 | Range-maximising embeddings over the shape space | alg | G10, rig-9, opt-4, adv §6 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | exact-`theta_max` optimiser vs `optimize_collision_sweep` (gamma ladder) on 30 graphs, same bisection referee | S+L | **PURSUE.** The only idea with a measured baseline failure to beat (P3). Merges into R1. |
| B11 | Mobility operator `m = dim ker A - c(Gamma)`, `A` is `2b1(Gamma) x \|F\|` | char | rig-1, U4 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | vs `3\|F\| - rank(rigidity matrix) - 3c` on 500 graphs (21/21 already) | S+M | **PARAGRAPH.** Instant-centre kinematics reduced; Aronhold–Kennedy falls out. Elegant, not a result. Its one new sentence (`L = B(sigma,·)`, `A = B(·,X)`) is notation. |
| B12 | The flat state is a branch point: mobility drops at `theta = 0` | thm | rig-2, U5 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | `m(theta)` at 9 angles on 500 graphs; needs the pencil (B13) to be non-numerical | M+S | **PURSUE-IF-KILL-PASSES**, and only after B14. Measured 20/21 but at `theta = 1e-6` with a fixed SVD tolerance — the rigidity persona's own self-attack is right that this is the signature of a bad rank. |
| B13 | Bifurcation angles = generalised eigenvalues of the pencil `A_c + t A_s` | thm | rig-5, A2 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | pencil eigenvalues vs dense rank sweep over 2000 angles, 500 graphs; two ways to die | M+M | **PURSUE-IF-KILL-PASSES.** Genuinely 50/50, one afternoon. Must ship with a second-order test at each candidate (A2 self-attack), i.e. merged with B16. |
| B14 | 2-core theorem: `dim ker A = \|F \ core2\| + dim ker A\|core2` | thm | rig-7, A6 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | exact integer rank identity on 500 graphs + `m_core` on the Delaunay patches | S+S | **RUN FIRST, NOT A PAPER.** Project hygiene: it decides whether U4's headline "mobility 19–26" is real. Adversary's self-attack is correct that it is a textbook reduction. |
| B15 | Generic embeddings in `X` have minimal mobility `m*` | thm | rig-8 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | 30 random points of `X` per structure, is `m` constant on >= 95%? | M+M | **DROP.** Routine upper-semicontinuity; and `m*` is contaminated until B14 lands. |
| B16 | Prestress certificate: which flat-state flexes extend | thm+alg | rig-10, A9 | 1 | 1 | 1 | 0 | 1 | 1 | 1 | certified-extending count vs measured `m(0+)` on 500 graphs | L+L | **DROP as standalone; keep as B13's second half.** A9 concedes it needs B14 first and lives in `P^{m(0)+2}`; C4 fails at 5,000 faces. |
| B17 | Eq. (9) has a measurable false-negative rate; separation is `2t(R + Pt)` | char+alg | A1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 500 graphs: count Eq.(9)-certified designs with a re-closure root `t* = -R/P` below `theta_max` | S+M | **PURSUE, as a section of R1.** The algebra is one line (adversary's own reviewer says so). P3 makes the empirical half stronger than A1 realised: on snub square Eq. (9) does not merely mis-certify, it achieves nothing for any gamma. |
| B18 | Certified `O(n)` active contact set (locality theorem) | thm+alg | A4 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | does `max_f rho_f / r_f` grow with patch diameter (log-log slope > 0.3)? pruned vs unpruned `theta_max` to 1e-12 | L+M | **PURSUE-IF-KILL-PASSES.** The only place a real theorem is needed rather than convenient, and it can genuinely be false (`w_f` is an alternating path sum). It is what makes R1 scale past ~1,000 faces. |
| B19 | Which periodic Jacobians `W` are achievable; the `nu(theta)` family | char | A5, opt-1b, opt-10 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | rank of `d vec(W)/dX` at 100 random points of `X` on 30 periodic graphs; vacuous if <= 1 on half | M+M | **PURSUE-IF-KILL-PASSES**, as the *second* paper. The right question behind B2, and it speaks to the auxetics community the CG papers talk past. Gated on periodic `dim X`, which is 6–8 in the paper's own examples. |
| B20 | Fully closed ⟺ root coincidence; hole area is a single harmonic, `theta_closed = 2 arctan(r/q)` | char+alg | A3, opt-2 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1e4 points of `X` with Eq. (14) residual < 1e-10; is total hole area at `theta_max` ever > 1e-8? | M+M | **PURSUE-IF-KILL-PASSES**, low priority. Improves one short subsection (2026 §5.2). Optimizer-2's two independent derivations of the same angle is the cheapest self-check in any of the four files. |
| B21 | The usable region `U(eps) = {X embedding} ∩ {theta_max >= eps}` can be small or empty | char+thm | A10 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | sampling failure rate over 1e4 points of `X` on 500 graphs; SDP infeasibility certificate | M+L | **PURSUE — merges with B8.** F17 already supplies the headline (0/57). The SDP-in-C++ risk (directive D3) is real; the sampling half is not. |
| B22 | Deployability defect `D(sigma)` is the right orientation objective, not max-cut | alg+ctrex | opt-5, opt-10 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | brute-force all `2^F` sigma on 40 graphs with `F <= 18`: is `argmin \|E_split\| = argmin` repair distance? | S+M | **PURSUE-IF-KILL-PASSES**, third tier. Genuinely nobody's else's idea (confirmed by the optimizer's own overlap audit). But P1 sharpens the tension it exploits: fewer split cuts ⇒ *smaller* shape space ⇒ less room to repair. That is the paper, not the solver. |
| B23 | The SVD basis is the wrong basis; a sparsest design-space basis | alg | opt-7 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | matroid-greedy sparsest basis vs SVD basis support size on 200 graphs | S+M | **DROP.** Fixes a sentence about the paper's own UI. Modest ceiling, and it is a nicer-basis paper. |
| B24 | Certified global inverse design by low-rank relaxation | alg | opt-9 | 1 | 1 | 1 | 0 | 1 | 1 | 1 | is the dual bound within 3x of the best achieved value on any instance? | M+L | **DROP as standalone; keep the certificate.** Same SDP-in-C++ risk as B21, without B21's measured motivation. Transplant the Shor/Burer–Monteiro bound into R1. |

**Constraint failures worth naming.** Only three ideas fail a constraint outright:
B1 and B4 fail **C6** (a reader of both papers, told the result, would write them as a future-work
sentence — indeed B4 *is* Euler's formula); B9 fails **C6** twice over, because 2026 §6 names the
method (sensitivity analysis) and 2025 Fig. 10 already shows and fabricates the phenomenon; B6's
algorithmic half fails **C3** (Hadlock 1975). B16 and B24 fail **C4** at the 5,000-face end.
Nothing in any file fails C7 — no persona proposed fabrication, UI or user studies, which is to their
credit and means C7 does no discriminating work here.

---

## 2. Bundles: which sets of ideas are one paper

### 2.1 The convergent bundle as handed over ("closed-form theory of uniform deployment")
B1 + B2 + B3 + B4 + B11. **Verdict: a collection of paragraphs.** I agree with the adversary's Part 1
and add two independent reasons.

First, the screen now confirms the shape of the problem. `notes/screen_bundle.md` S1 returns
**PARTIAL** and names Acuna et al. 2022 (Comm. Phys. 5, arXiv:2101.12352), which establishes
counter-rotation on *arbitrary bipartite planar networks* — random networks, a Penrose quasicrystal,
a disordered isotropic network, 3D-printed and FEM-validated — and is cited by neither Segall paper.
S6 returns **PARTIAL** for the mobility half. So two of the bundle's five claims have a named living
competitor, and the surviving novelty in each is a hypothesis-removal, not a discovery.

Second, and decisively: **the bundle's load-bearing theorem, if you write it out, is a definition.**
Uniform deployability says every hinge opens by `theta`; `Gamma` is bipartite by construction (F1);
therefore face rotations are pinned up to a global rotation. There is nothing to prove. Everything
downstream (B2, B3, B20, B19) is a corollary of *substituting* that closed form into an existing
question. A corollary family with no theorem at its head is a preliminaries section.

**What the bundle is actually good for.** It is an excellent §3. Lemma 1 = B1 (with explicit credit
to the rotating-rigid-units literature and to Acuna et al. for the counter-rotation half), Lemma 2 =
B4, Lemma 3 = the adversary's **no-locking lemma** (`sigma ∈ ker A(Y_theta)` for every `theta`, so the
uniform branch never reaches a dead centre and every stop is a contact — this is what makes any
contact enumeration *complete* rather than merely necessary, and no persona but the adversary saw
it), Remark = B2. Four pages of machinery that licenses a real result. It is not the result.

### 2.2 The adversary's "smallest paper" ("the deployment path is a conic")
B1 + B4 + no-locking + B18 + B17 + B10 + B13. **Verdict: a real paper, but it leads with the wrong
thing.** Its §4 headline is an exact contact enumeration, and the adversary's own Objection 2 is
correct that a minimum over `O(n^2)` candidate pairs is not "a geometric characteristic directly from
the embedding" — it is a better simulator, and PyKirigami (Aug 2025, rev. Feb 2026) is already the
numerical version of that. A referee who reads §4 first sees a kinetic-data-structure exercise.

### 2.3 My proposed bundle: **the usable shape space** (R1 below)
B8 + B21 + B3 + B17 + B10 + B18, with B1/B4/no-locking/B7 as machinery. **Verdict: this is the
paper.** It differs from 2.2 by one move: it leads with the *design space*, not with the *simulator*.

- **Load-bearing theorem.** The 2026 abstract claims to "characterize the full space of embeddings
  that admit uniform deployment", and Eq. (5) delivers a *linear* space `X` of dimension `|E_split|`
  (F15). The theorem is that the geometrically **usable** subset is not linear and is exactly
  computable:
  > Because the deployment path is a conic, every validity predicate — the signed area of a face at
  > `theta = 0`, and the separation of any (vertex, edge) pair at any `theta` — is a single harmonic
  > `p + q cos(theta) + r sin(theta)` whose coefficients are **quadratic forms in `X`**. Hence
  > `U(eps) = {X in X : all faces positively oriented, theta_max(X) >= eps}` is a **basic
  > semialgebraic set cut out by quadrics**, with a finite, explicitly enumerable constraint list,
  > and membership is decidable in closed form with no forward-kinematics simulation.

  That is a statement about the design space, of the kind 2026 §6 asks for verbatim ("geometric
  characteristics for favorable deployment behavior **directly from the embedding**"), and it is not
  a restatement of a definition. The no-locking lemma is what makes the constraint list *complete*.

- **Headline algorithm.** Project onto `U`: maximise the exact `theta_max` over the null-space
  coefficients subject to positive signed areas, using exact gradients from the closed-form roots and
  the active set from B18. The baseline is the paper's own Eq. (6) + Eq. (9) pipeline as implemented
  in `code/src/core/collision.hpp` with the documented gamma ladder.

- **Headline figure.** One row per graph: `X_ini` (valid), Eq. (6) `X0` (invalid — 5.9% of faces
  inverted, F17), Eq. (9) repair (still invalid, and on snub square still `theta_max = 1.646` against
  a kinematic ceiling of 3.384, F18/P3), ours (valid, certified, near the ceiling). That single strip
  is the whole paper and the numbers for its first three columns already exist in
  `results/core_validation/`.

- **Why both outcomes are publishable.** If a valid large-range point of `X` is usually findable, the
  paper is an algorithm with a certificate. If it usually is not, the paper is A10's emptiness result:
  *the space the abstract calls "full" is generically unusable*, with an infeasibility certificate.
  This is the only idea in all four files where the kill experiment cannot leave you with nothing,
  which is why it is my first commitment.

### 2.4 The mobility bundle (a second paper, not this one)
B11 + B14 + B12 + B13 + B16. Coherent, but every number in it is currently contaminated: U4's
"mobility 19–26" is measured at a singular configuration (B12) on free-boundary patches with dangling
faces (B14). B14 must be run before any mobility number is written anywhere, including in `STATE.md`.
After that, the bundle's real content is B13 + B16: an exact bifurcation spectrum with a second-order
certificate at each candidate. That is a good paper and it is 50/50. It is not the first commitment
because its most likely outcome (an empty bifurcation set) leaves nothing.

### 2.5 The rank/algebra bundle
B4 + B5 + B7 + B6's bound. Coherent as a *correction* to 2026 §4.4, which is the paper's load-bearing
unproved claim. But F15 settled the fixed-boundary case (`rank(L) = H`, corank 0) and the
orchestrator's new checks (P4) have now *measured* the periodic failure `rank(L) = H - 1` on three
tori, together with `L = R·D` (50/50) and `rank(L) = H - dim Z` (27/27). The bundle is therefore
already largely executed, and what it produced is a correction, not a theorem.
**Verdict: one section of somebody else's paper.** Worth finishing because it is nearly free.

---

## 3. The top five, each with the strongest rejection I can write, a fair rebuttal, and a verdict

Ranked by (probability the kill test passes) x (value if it does) / cost, with paper-worthiness and
project-value separated where they disagree.

---

### R1. The usable shape space — certified valid, range-optimal Tutte auxetic embeddings
**= B8 + B21 + B3 + B17 + B10 + B18.** Type: theorem family + algorithm + measured margin.
Targets 2026 §6 limitation 1 verbatim (both halves: "solutions may exhibit extremely short edges or
sharp angles ... self-intersections may occur at small opening angles ... Developing geometric
characteristics for favorable deployment behavior directly from the embedding remains an open
problem") and 2026 open question 12.

**Hostile TOG reviewer.**
> The submission's mathematical content is that a rigid one-parameter motion is trigonometric in its
> parameter, so that polynomial predicates evaluated along it are first harmonics. This is the
> Weierstrass substitution. Everything the authors call a theorem follows from it in one line, and
> the authors concede as much: their "characterization" of the usable set is the observation that a
> conjunction of quadratic inequalities defines a semialgebraic set, which is the definition of a
> semialgebraic set.
>
> The empirical centrepiece is that the original authors' Eq. (6) projection produces
> self-intersecting embeddings on 57 of 57 random graphs. But Equation (6) is a least-norm projection
> onto a linear space; nobody, least of all Segall et al., ever claimed it was injective — they say
> so themselves in Section 6 and offer post-processing. Reporting that an explicitly non-injective
> projection is not injective, on *random Voronoi and Delaunay patches* that appear in none of the
> original paper's figures, is a strawman with a large N. Show me the failure on the patterns the
> paper actually ships.
>
> The algorithmic comparison is worse. Equation (9) is stated in the original with an unspecified
> barrier `B`, an unspecified weight `gamma`, and an unspecified diameter rule. The submission
> reconstructs all three, sweeps `gamma` over a ladder of its own choosing, and then reports beating
> the reconstruction. That is not a baseline; it is a self-portrait. The authors ship a WebAssembly
> implementation at a public repository, and until it is run the headline margin measures the
> submission's reading of a paragraph.
>
> Finally, the honest version of the algorithm is: add inequality constraints to a linear solve and
> run a local optimiser in the null space. Constrained parameterisation with local-injectivity
> barriers is a solved problem in geometry processing with a decade of literature. Reject.

**Fair rebuttal.**
The strawman charge is the only one that lands, and it lands only halfway. The 2026 title is
*Uniformly Deployable Kirigami on **Arbitrary Planar Graphs***, and §5's opening sentence advertises
non-periodic and non-2-colourable graphs as in scope. Voronoi and Delaunay patches are the arbitrary
planar graphs of the title. A failure rate of 57/57 in the advertised regime, against 0/8 on the
authored tilings shown in the figures, is precisely the gap between claim and demonstration, and it
is a finding about scope, not a strawman — provided the paper states both numbers side by side, which
it must. The reconstruction charge is real and unfixable by argument: the only answer is to run
`github.com/segaviv/tuttekiri` and report against it, which becomes step zero of the project.

On triviality: the reviewer is right that the harmonic identity is one line, and any submission that
sells it as the theorem deserves this rejection. The defensible claim is narrower and is not one line
— that the usable region has a **finite, complete, exactly enumerable** constraint list, which
requires the no-locking lemma (every terminating event is a contact, never a kinematic dead centre)
and requires B18's locality theorem to make the list `O(n)` with a certificate rather than `O(n^2)`
by assumption. Neither is trivial and B18 can be false.

On the "solved problem" charge: constrained parameterisation optimises over vertex positions with a
barrier. Here the feasible set is a *linear subspace of dimension `|E_split|`* (F15, up to 1017
measured, F19) inside which every point is exactly uniformly deployable, and the objective is a
collision-free *range* over a whole one-parameter motion, not injectivity of a single map. The prior
art the reviewer is reaching for does not have either object.

**Verdict: PURSUE.** First commitment. It is the only idea in the four files whose kill experiment
cannot return "nothing": either valid range-optimal points exist and there is an algorithm, or they
do not and there is an emptiness result about the word "full" in the abstract.

---

### R2. Exact contact calculus with a certified `O(n)` active set
**= B3 + B18 + B17, the theorem core of R1, judged on its own** in case R1's algorithmic half is
absorbed elsewhere. Targets 2026 §6 limitation 1 and §4.5's empirical "collisions occur only along
split-cut edges".

**Hostile TOG reviewer.**
> Two prior results bracket this submission and it cites neither prominently enough. Segall et al.
> 2025 §4.2 already gives a closed-form `theta_max` — `min over hinge vertices of (2pi - alpha_i -
> alpha_j)` — for the hinge-adjacent case. Liu et al. 2024 (Graphical Models 133, 101215), which
> Segall et al. 2026 cite themselves, already give a closed-form non-adjacent collision condition.
> The submission's contribution is therefore the non-adjacent complement of one formula and the
> exact version of another, for a graph class neither covers. That is a completion, not a discovery,
> and the abstract must say so.
>
> The submission's own validation undermines its motivation. It reports that on every split-free
> pattern the exact `theta_max` equals `min_i beta_i` — i.e. equals the 2025 local formula, exactly.
> So the entire apparatus is inert on the pure hinge case and acts only when split cuts are present.
> Since split cuts are what the original paper introduces, the contribution is a repair of a
> complication the original authors created, on a family of patterns the original paper handles by
> an optimisation it does not fully specify.
>
> The locality theorem is asserted where it matters least. The swept-disc radius `rho_f` is bounded
> by a path sum over a spanning tree of the hinge graph; the authors' own risk section concedes this
> may grow with patch diameter, in which case the `O(n)` claim is false and the contribution reduces
> to a broad phase, which is a spatial hash. A conditional theorem whose hypothesis the authors
> cannot verify is a conjecture with an experiment attached.

**Fair rebuttal.**
The bracketing is accurate and must be stated in the first paragraph, not defended against: this
completes Segall 2025 §4.2 to the non-adjacent case and generalises Liu et al. 2024 from conservative
sufficient conditions on enumerated quad pairs to exact roots on arbitrary planar graphs. The
reviewer's second paragraph is the strongest point and it is *my* finding (P1), not the submission's
weakness: `theta_max = min beta_i` exactly when `E_split = {}` (F18), and `dim X = |E_split|` (F15),
so the exact calculus is non-trivial on precisely the patterns that have a design space to optimise
over. Stated that way it is not a complication-repair, it is the observation that the paper's design
freedom and its collision failures have a single common cause. That reframing is worth the section.
The locality objection is correct and is why B18 is scored as a kill-experiment gate rather than a
result: if `rho_f/r_f` grows with diameter, the theorem is false and R2 must ship as a broad phase.

**Verdict: PURSUE-IF-KILL-PASSES.** Do not commit to it separately from R1 — it is R1's §4.
Commit to B18's locality kill (K2c) before anyone writes the word "theorem".

---

### R3. The bifurcation spectrum of the uniform branch, with a second-order certificate
**= B13 + B16 + B12, gated on B14.** Targets 2026 §6 limitation 2 and the unstated assumption in
both papers and in IsoGami that the deployment path is unique.

**Hostile TOG reviewer.**
> The submission proves two things and both are immediate. That the uniform branch never locks
> follows from the branch existing and therefore having a tangent. That the mobility matrix along a
> path parameterised by `cos(theta/2)` and `sin(theta/2)` is an affine pencil follows from the
> matrix's entries being affine in those quantities. Presenting these as Lemma 1 and Theorem 1 is
> padding.
>
> The substantive claim is that bifurcations occur inside the collision-free range and change the
> physical behaviour. A rank drop of a first-order operator is a *necessary* condition for a
> bifurcation, not a sufficient one; Connelly and Servatius showed thirty years ago that first-order
> degeneracy can be an artefact of the linearisation, and this submission's own supporting evidence
> is a numerical rank taken with an unjustified fixed tolerance, on data that contains an unexplained
> counterexample (the triangular tiling, where the predicted drop does not occur). With no energy, no
> actuation model, and no stability argument, the claim that a physical sheet leaves the uniform
> branch is unsupported by anything in the paper. What is reported is a spectrum of candidates.

**Fair rebuttal.**
The first paragraph is fair and the lemmas should be labelled as such. The second is the reason the
idea must ship as B13 **plus** B16: a second-order prestress test at each candidate angle, which is
cheap there precisely because `m(theta_b)` is small away from the flat state — far cheaper than A9's
original proposal of doing it at `theta = 0` where `m(0)` reached 26. The triangular-tiling exception
is not a hole to be hidden; it is the one graph in the 21 with an explanation to find, and finding it
is the first piece of real work. The genuine defence is that this is the only question in the four
files whose answer nobody knows: PyKirigami (Aug 2025, rev. Feb 2026) still detects locking states
numerically, and IsoGami's continuation solver would cross a bifurcation without reporting it.

**Verdict: PURSUE-IF-KILL-PASSES.** Run K3a (2-core) first, which costs an hour and is mandatory
anyway; then run the pencil kill. If the bifurcation set is empty on 500 graphs, abandon within the
day — the surviving content is two remarks.

---

### R4. Which periodic deployment Jacobians are achievable; the designable Poisson-ratio family
**= B19 (A5 + optimizer's affine-achievable-set test + optimizer-10).** Targets 2026 §5.1's
existence claim, quoted verbatim at `notes/paper_2026.md:820`: "Empirically, we observe that for all
tilings with a non-trivial kernel (`E_split != {}`), solving Eq. (13) yields an embedding with
conformal deployment."

**Hostile TOG reviewer.**
> Closed-form Poisson ratios for rotating rigid units have been published since Grima and Evans in
> 2000, for squares, rectangles, triangles and rhombi, and Grima et al. (Proc. R. Soc. A 468, 2012)
> give them for generic rotating triangles. The submission's `nu(theta)` is those formulas with the
> unit cell left symbolic. Its "design map" is the observation that a matrix built linearly from the
> vertex positions depends linearly on the vertex positions.
>
> The paper's own examples have shape spaces of dimension six and eight. A four-dimensional target
> reached from a six-dimensional source is not a characterisation problem, it is a rank computation
> on a small matrix, and the submission reports it as one. If the rank turns out to be one — which
> the submission's own kill experiment is designed to detect — the entire achievable set is a curve
> and there is nothing to design.

**Fair rebuttal.**
The screen supports the delta: `notes/screen_bundle.md` S4 returns **NOT FOUND** for the
`cos(theta/2) I + sin(theta/2) K` decomposition, having read Czajkowski et al. (Nat. Commun. 13, 211,
2022) in full and established that their conformality is *spatial* and assumed, not *temporal* and
derived. The published `nu(theta)` formulas are per-structure; the object here is the image of a
design map over a characterised linear space, which did not exist before 2026 Eq. (5). The reviewer's
second paragraph is the real risk and the kill experiment is built around it.

**Verdict: PURSUE-IF-KILL-PASSES**, but as the *second* paper, not the first. It is the only idea in
the four files aimed at the auxetics/metamaterials audience rather than at the graphics one, which is
strategically valuable and tactically slower.

---

### R5. The orientation objective is wrong: deployability defect beats max-cut
**= B22 (optimizer-5 + optimizer-10).** Targets 2026 §4.2, p. 60:6 verbatim: "we empirically favor
many balanced, small holes rather than a few large ones. This motivates choosing face orientations
that minimize the number of split edges", which the Reader records as supported by "Nothing. No
definition of 'balanced', no experiment."

**Hostile TOG reviewer.**
> Section 4.2 of the original is a paragraph of design intuition guarding a heuristic that the
> original authors explicitly say the user may override by clicking faces. The submission proves
> that this heuristic is not optimal for an objective the original authors never adopted. Refuting a
> proxy nobody defended, against a target nobody set, on graphs with at most eighteen faces because
> the comparison requires enumerating `2^F` orientations, is not a contribution.
>
> Worse, the submission's own framework contradicts its recommendation. It reports that the shape
> space has dimension `|E_split|`. Minimising split cuts therefore minimises design freedom, and the
> submission elsewhere spends a section arguing that design freedom is what the original paper fails
> to use. Which is it? A paper cannot argue that fewer split cuts are better and that more split cuts
> are better in adjacent sections.

**Fair rebuttal.**
The contradiction the reviewer names is not a flaw, it is the result, and stating it that way is the
whole rescue: `dim X = |E_split|` (F15) and split cuts are also the sole source of non-local collision
loss (F18) — so the paper's stated principle "minimise split edges" trades away the exact resource
that the rest of the pipeline needs to repair itself. That is a coherent, testable, *quantified*
criticism of a design principle stated with zero evidence, and it needs no `2^F` enumeration to make:
it needs a Pareto plot of `|E_split|` against achieved range and against Eq. (6) repair distance over
500 graphs. The enumeration is only needed for the strong existence half, which is a footnote.

**Verdict: PURSUE-IF-KILL-PASSES**, third tier. It is a genuinely unclaimed idea (the optimizer's own
overlap audit confirms neither other persona proposed it) and it costs little, but it belongs inside
R1 as the section that explains *why* the shape space is the size it is.

---

### Just below the line
**B14 (2-core)** is not in the top five because it is not a paper — the adversary's self-attack is
right that it is a textbook reduction. It is nevertheless the **first thing to run**, because U4's
headline mobility number is in `STATE.md` and is probably an artefact of dangling faces.
**B20 (fully-closed = root coincidence)** and **B5 (periodic corank)** are cheap, decisive, and worth
one section each; neither carries a paper. **B9 (non-uniform deployment)** is excluded on C6 by
2025 Fig. 10 (see Sec. 0), and any future proposal in that direction must open by citing it.

---

## 4. Kill-experiment specs for the top three

Written so an Experimenter can run them against `code/src/core` with no further design decisions.
All C++17/20 per directive D3; Python only for matplotlib on the dumped CSV. All routines named below
exist in `code/README.md`'s library table and CLI list. Note deviation 5: the dense SVD path is used
for `N <= 1400` columns; above that only `rank` is available via `SparseQR`, so every experiment that
needs a **null-space basis** is capped at `N <= 1400` vertices until a sparse null-basis routine is
written. That cap is a build item, not a licence to skip the large end.

### K1 — for R1 (usable shape space). Three parts; K1b is shared with K2.

**K1a — is the baseline really invalid, and is the shape space really usable?**
- Generators: `kiri_gen voronoi`, `kiri_gen delaunay`, `kiri_gen quad_random`, seeds `1..200`
  (≈67 each), point counts chosen to give 100–800 faces, `N <= 1400`.
- `sigma`: `assign_orientation_relaxation` (Eq. (1) + the documented diameter/repair rules).
  Boundary: `BoundaryMode::Fixed`.
- Routines: `make_cut`, `holes_partition`, `solve_system` (null basis `Phi`, `X0` = Eq. (6)),
  `signed_area` per face of `M` under a candidate `X`.
- Quantities per graph: `dim_null`; `n_inv(X_ini)` (must be 0 — assert); `n_inv(X0)`; then draw
  `10^4` samples `X = X0 + Phi t`, `t ~ N(0, s^2 I)` with `s` set so the median `||Phi t||_inf`
  equals the median edge length of `X_ini`, and record `p_valid` = fraction with `n_inv = 0`, and the
  best (largest `theta_max`) valid sample.
- Tolerance: a face counts as inverted if its signed area `<= 1e-12 * (bounding box area)`.
- PASS/FAIL. **PASS-algorithm** if `n_inv(X0) > 0` on `>= 150/200` graphs and `p_valid > 0` on
  `>= 20` graphs. **PASS-emptiness** if `n_inv(X0) > 0` on `>= 150/200` and `p_valid = 0` on
  `>= 20` graphs. **FAIL** (premise dead, F17 was a sweep artefact) if `n_inv(X0) = 0` on `> 100/200`.
  Both PASS branches keep R1 alive; they select which paper it is.
- Wall time: 200 x (one dense SVD at `N <= 1400`, worst 3.4 s measured) + `2 x 10^6` cheap area
  evaluations ≈ **12 min**.

**K1b — the harmonic identity (load-bearing for R1 and R2; run this first, it is 5 minutes).**
- 200 graphs from K1a's set, 20 random `X in X` each. For 50 random ordered (vertex, edge) triples
  per graph, sample `deploy()` at 200 angles on `(0, min(pi, theta_max))`, least-squares fit the
  orientation determinant to `p + q cos(theta) + r sin(theta)`.
- PASS if the relative residual is `< 1e-10` on every triple of every graph. FAIL on any exceedance —
  and then B1/B3/B17/B20 all die together, so nothing else may be started before this returns.
- Also assert face signed areas are constant in `theta` to `1e-12` (rigidity of faces).
- Wall time: **≈ 5 min**.

**K1c — Eq. (9)'s false-negative rate (B17), on the same 200 graphs.**
- Run `optimize_collision_sweep` (the documented gamma ladder `{1e-4, ..., 10}`, best kept) to get
  `Y_9`. For every split-edge pair evaluate `P, R` in closed form and count pairs with `P < 0`,
  `R > 0`, and `theta* = 2 arctan(-R/P) < theta_max(Y_9)`.
- Cross-check every reported `theta*` against `collision.hpp`'s bisection to `<= 1e-6` rad.
- PASS if the fraction of gamma-ladder-certified designs with at least one re-closure below
  `theta_max` is `>= 3%`. FAIL below 3% — B17 dies, R1 survives without its §5.
- Wall time: **≈ 8 min** (one Newton solve per graph dominates).

### K2 — for R2 (exact contact calculus). Three parts; K2c is the gate.

**K2a — exact `theta_max` agrees with bisection.**
- 300 graphs: the 8 `kiri_reference` cases + 292 random from the three generators, 100–800 faces.
- Compute `theta_max` two ways: (i) closed-form smallest positive root over all ordered
  (vertex, edge) pairs, each root filtered by the two closed-form interval tests
  `0 <= dot(y_b - y_a, y_v - y_a) <= ||y_b - y_a||^2` evaluated at that root; (ii) `theta_max` from
  `collision.hpp` (grid scan + bisection, faces shrunk by relative `1e-6` per deviation 9).
- PASS if `|theta_max_exact - theta_max_bisect| <= 1e-5` rad on every graph (the tolerance is set by
  deviation 9, not by the closed form). FAIL on any exceedance.
- Secondary output, reported either way: whether the binding pair is a split-edge pair. If a non-split
  pair ever binds first, 2026 §4.5's stated empirical assumption is falsified — a result in itself.
- Wall time: **≈ 15 min** at `O(n^2)` pairs, 800 faces.

**K2b — the range-optimisation margin (this is B10 and the algorithmic half of R1).**
- 30 graphs: `snub_square`, `truncated_square`, `t3_4_3_12`, `hexagons` (the four reference cases
  with `E_split != {}`), plus 26 random with `dim_null >= 2`, `N <= 1400`.
- Baseline: `optimize_collision_sweep`, full gamma ladder, best `theta_max` kept, `T = 0` fallback.
- Ours: maximise a log-sum-exp softmin of the closed-form first-contact roots over the null-space
  coefficients using `lbfgs_minimize`, with a log barrier on every face's signed area, active set
  refreshed every 10 iterations.
- **Both** measured by the same `collision.hpp` bisection, never by their own objectives.
- PASS if the median relative gain over the best-gamma baseline across the 30 graphs is `>= 25%`
  **and** snub square closes `>= 30%` of its measured gap (from 1.646 toward `min beta = 3.384`).
  FAIL otherwise.
- Wall time: **≈ 25 min**.

**K2c — the locality theorem gate (B18). Run before anyone writes "theorem".**
- 500 graphs, 100–5,000 faces, three generators (rank-only path is fine above `N = 1400`; this test
  needs no null basis).
- Per face compute `rho_f = max over its vertices of max(||x_u||, ||chi_u||)` in the face's own frame
  and `r_f` its circumradius; fit `log(max_f rho_f / r_f)` against `log n`.
- Then compute exact `theta_max` with and without the swept-disc pruning `||c_f - c_g|| <= rho_f +
  rho_g`, on the subset with `n <= 800`.
- PASS if the log-log slope is `<= 0.3` **and** pruned/unpruned `theta_max` agree to `1e-12` on every
  graph. FAIL on either — B18 dies and R2 ships as a heuristic broad phase, which is honest but is
  not a theorem, and R1's §4 must be retitled.
- Wall time: **≈ 20 min**.

### K3 — for R3 (mobility and bifurcation). Two independent parts, both cheap.

**K3a — the 2-core identity and the intrinsic mobility (B14). RUN THIS FIRST OF EVERYTHING.**
- 500 graphs: three generators (100–5,000 faces) + the 8 reference cases. `sigma` as in K1a, `X` the
  solved `X0` where available, else `X_ini`.
- Build `Gamma`, its 2-core by iterated deletion of degree-<=1 nodes, and `A` (`2 b_1(Gamma) x |F|`)
  from a spanning-tree cycle basis, per rigidity (R2).
- Check `dim ker A = |F \ core2(Gamma)| + dim ker A|core2` as an integer identity, ranks by
  `ColPivHouseholderQR` with threshold `1e-10 * ||A||_2`.
- Report `m_core = dim ker A|core2 - c(core2)` on the random Delaunay patches where U4 reported
  `m = 19-26`.
- PASS if the identity holds 500/500 **and** `m_core <= 5` on `>= 90%` of those patches (confirming
  the adversary's prediction that U4's headline is a dangling-face artefact). If the identity fails on
  any graph, `A` is misassembled and B11–B16 all stop until it is fixed.
- Wall time: **≈ 10 min** (`A` has only `|F|` columns).

**K3b — the periodic rank claim (B5, and the whole of Sec. 2.5).**
- `kiri_sweep`-style driver with `BoundaryMode::Periodic` on `periodic_squares`, `periodic_triangles`,
  `periodic_hexagons`, `periodic_kagome`, 25 sizes each (100 instances).
- Measure `rank(L)`, `H`, `||1^T L||_inf`, `dim ker[L; B]`.
- PASS/report if `||1^T L||_inf < 1e-12` and `rank(L) = H - 1` on all 100. **Already measured on 3
  tori** (`rank_claim.md` orchestrator checks), so this is a widening, not a discovery: 2026 §4.4's
  "number of independent equations equals the number of holes" is already known false as stated for
  periodic patterns. The remaining question, and the only one worth the 5 minutes, is whether the
  extra kernel direction is a **usable** design degree of freedom or merely the global translation
  that Remark 4.1 already gives away. Report `dim ker[L;B]` with and without the translation
  quotiented, and whether any point along the extra direction is a valid embedding.
  If `rank(L) = H` on any periodic instance, F14 is wrong and must be retracted from `STATE.md`.
- Wall time: **≈ 5 min**.

**K3c — the bifurcation pencil (B13), only after K3a passes.**
- 500 graphs + the 8 reference cases. Assemble `A_c, A_s` at the solved `X0`; compute the real `t`
  where `rank(A_c + t A_s)` drops, by a dense SVD sweep over 2000 values of `t` cross-validated
  against the generalised eigenvalues of a maximal-rank square subblock.
- PASS if the two agree on every graph **and** at least one `theta_b` lies strictly inside
  `(0, theta_max)` on `>= 5%` of graphs. FAIL-alpha if the rank is constant in `t` on all 500 (the
  bifurcation set is empty; R3 collapses to two remarks — abandon the same day). FAIL-beta if the
  sweep and the eigenvalues ever disagree (the pencil formulation is wrong; everything in R3 dies).
- Wall time: **≈ 20 min**.

### Shared kill experiments
- **K1b is shared by R1 and R2** and is a prerequisite for both. Run it before anything else except
  K3a, which is independent.
- **K1c and K2b share the Eq. (9) baseline run** (`optimize_collision_sweep` on the same graphs);
  compute both from one pass.
- **K3a is a prerequisite for B11, B12, B13, B15 and B16** and must be run before any mobility number
  is quoted, including the `m = 19-26` currently in `STATE.md` U4.
- **Step zero, not a kill experiment but a blocker for R1's and R2's headline margins:** obtain and
  run `github.com/segaviv/tuttekiri`. Until then every "beats the baseline" number is against a
  reconstruction, and both the adversary's self-attack and my R1 reviewer paragraph are unanswerable.

---

## 5. Novelty holes: claims resting only on the field tables, and the exact queries to run

`notes/screen_bundle.md` has closed S1–S7 for the convergent bundle. Everything below is **not**
covered by S1–S7 and rests only on `notes/field_kirigami.md`, `notes/field_rigidity.md`,
`notes/field_tutte.md` or on a persona's assertion. Queries are given verbatim for the Phase-6 Scout.

| # | Claim resting only on field tables | Which idea depends on it | Query strings to run |
|---|---|---|---|
| N1 | No injectivity theorem exists for a Tutte-type linear system with **extra face-indexed rows** and **zero weights on part of each star**; and no one has characterised the valid subset of such a solution space. Rests on `field_tutte` rows 1–4, 10, 14, 16 and route C. | **R1 / B8 — the first commitment.** | `Tutte embedding additional linear constraints injectivity theorem convex combination zero weights partial star`; `"convex combination" map planar graph "not injective" counterexample subset of neighbours barycentric`; `constrained parameterization local injectivity linear subspace null space "flip-free" certificate`; `Aigerman Lipman orbifold Tutte bijective extra constraints generalization arbitrary linear rows` |
| N2 | Nobody has an exact semialgebraic description of the collision-free parameter interval of a **one-DOF planar rigid-body linkage** as a function of design variables. Rests on `field_kirigami` R18 and the adversary's Q6/Q8 only. | **R1 / R2 — B3, B21.** | `collision-free interval one degree of freedom planar linkage exact semialgebraic characterization design parameters`; `"kinetic data structure" rotating rigid bodies exact first contact time trigonometric polynomial certificate`; `self-intersection free range rigid body motion "quadratic in tan(theta/2)" exact event enumeration` |
| N3 | The swept-disc / locality bound making the active contact set `O(n)` with a certificate. Rests on nothing but the adversary's A4 sketch. | **R2 / B18** — the gate. | `swept volume bound rigid tiling mechanism deployment broad phase provable O(n) candidate pairs planar`; `alternating path sum spanning tree bounded translation rotating rigid units patch diameter` |
| N4 | The 2-core reduction of mobility for a **body-and-pin** framework has not been written for kirigami. Rests on `field_rigidity` R1, R7, R26 and the adversary's own concession that it is textbook. | **K3a / B14** — needed to state it honestly, not to justify it. | `2-core reduction degrees of freedom body-and-pin framework dangling bodies pebble game preprocessing`; `finite patch free boundary spurious degrees of freedom rotating squares "1 DOF" proof` |
| N5 | Nobody has an exact **bifurcation spectrum** for a kirigami/auxetic deployment; the state of the art detects locking numerically. Rests on `field_rigidity` R19, R20, R25, `field_kirigami` R22 and the adversary's Q8. | **R3 / B13, B16.** | `matrix pencil generalized eigenvalue bifurcation one-parameter mechanism rigidity matrix rank drop exact angles`; `PyKirigami locking states deployment kinematics bifurcation branch detection tessellation`; `Kumar Pellegrino bifurcation pin-jointed structures pencil closed form successor 2020..2026` |
| N6 | Nobody has studied the **image of the design map** `X -> W` for a periodic hinged tiling, i.e. which macroscopic Jacobians / Poisson-ratio curves are achievable. Rests on `field_kirigami` R1, R2, R18 and `field_rigidity` R17. | **R4 / B19.** | `achievable Poisson ratio range design space periodic rotating rigid unit tiling image of design map`; `"Poisson's ratio" as a function of unit cell geometry rotating polygons inverse design achievable set characterization`; `Borcea Streinu auxetic cone periodic framework achievable deformation Jacobian family` |
| N7 | The **complex loop-closure** form `{z in (S^1)^F : W z = 0}`. `screen_bundle` S5 already returns **KNOWN as a formulation** (classical complex loop equations in mechanism design), so the only surviving claim is "extra finite mode = complex rank drop of `W`, linear in `X`". | **B9 (dropped) / optimizer-8.** | `complex loop closure equations planar mechanism synthesis rank deficiency extra assembly mode`; `Kapovich Millson planar linkage configuration space linear section of torus rigid tiling` |
| N8 | The **orientation-choice objective**: nobody optimises the cut topology `sigma` against a target embedding. Rests on `field_kirigami` R9, R17, R22 and the optimizer's own overlap audit. | **R5 / B22.** | `kirigami cut topology optimization against target geometry orientation assignment objective not max-cut`; `choosing hinge versus separated cuts to minimize geometric distortion tessellation deployable` |
| N9 | The **second gap-free (fully closed) angle in closed form**, and the claim that Eq. (14)'s two conditions are not sufficient. `screen_bundle` S7 returns **NOT FOUND** for all three parts, so this is the best-screened item here; what remains unchecked is only the *insufficiency* claim. | **B20.** | `two compact states kirigami closure angle closed form condition necessary sufficient arbitrary polygon tiling`; `"fully closed" deployable tessellation exact closure condition hole area vanishes rotating units` |

**One item the screen has already changed and the orchestrator must act on.** `screen_bundle` S1's
action line stands: **Acuna et al. 2022 (Comm. Phys. 5, arXiv:2101.12352) must be added as a row of
`notes/field_kirigami.md`.** It establishes counter-rotation on arbitrary bipartite planar networks
(random, Penrose, disordered), it is cited by neither Segall paper, and it is the reference a referee
will raise against any statement of B1.

---

## 6. Final recommendation

**Commit first to R1 — the usable shape space** (`ideas/ranking.md` Sec. 2.3 and 3.R1), assembled as:

- §3 preliminaries: B1 (Lemma 1, with credit to the rotating-rigid-units literature and to Acuna et
  al. 2022), B4 (Lemma 2), the no-locking lemma (Lemma 3, the adversary's, and the one that makes the
  contact list complete), B2 as a two-sentence remark, B7 as the algebraic engine.
- §4 the theorem: the usable region `U(eps)` is a basic semialgebraic set cut out by quadrics, with a
  finite complete constraint list, made `O(n)` by B18's locality theorem **if K2c passes**.
- §5 the indictment: F17 (0/57 valid), F18/P3 (Eq. (9) achieves nothing for any gamma on the snub
  square, which has 13 dimensions of unused shape space), B17's false-negative rate.
- §6 the algorithm and the margin: range-optimal certified-valid embeddings vs Eq. (6) + Eq. (9) on
  >= 20 graphs, refereed by an independent bisection.
- §7 the honest limit: B22's Pareto plot showing the paper's "minimise split edges" principle trades
  away the exact resource (`dim X = |E_split|`) the repair needs.

**Why this and not the adversary's bundle.** Same machinery, one different lead. Leading with the
design space attacks the word "full" in the 2026 abstract and answers §6 limitation 1's actual
request ("characteristics **directly from the embedding**"). Leading with the contact calculus, as
the adversary proposes, invites the correct objection that an `O(n^2)` event sweep is a better
simulator, which PyKirigami already is numerically. And R1 is the only candidate whose kill experiment
cannot return nothing: K1a's two PASS branches select between an algorithm paper and an emptiness
paper, and both are publishable.

**Order of operations.**
1. **K3a (2-core), ~10 min, mandatory regardless of what is written.** It decides whether U4's
   `m = 19-26` in `STATE.md` is real. Retract or confirm that line before anything else quotes it.
2. **K1b (harmonic identity), ~5 min.** If it fails, B1, B3, B17, B20 and half of R1 die at once, so
   nothing may start before it returns.
3. **K1a (validity + usability), ~12 min.** Selects which R1 paper this is.
4. **K2c (locality gate), ~20 min.** Decides whether R1 §4 says "theorem" or "broad phase".
5. **K1c + K2a + K2b, ~50 min combined.** The indictment and the margin.
6. In parallel and cheap: **K3b (periodic rank), ~5 min** — a two-line correction to `STATE.md` and
   to 2026 §4.4 either way.
7. **Step zero in the background:** obtain and run `github.com/segaviv/tuttekiri`, without which every
   margin in §5 and §6 is against a reconstruction.

**Fallback order.**
1. **R3 (bifurcation spectrum)** — after K3a, run K3c the same afternoon. Genuinely 50/50, and it is
   the only question in the four files whose answer nobody in the field knows. If K3c returns
   FAIL-alpha, abandon within the day.
2. **R4 (achievable Jacobians / Poisson family)** — the second paper, aimed at the auxetics audience.
   Slower, needs a periodic-metatile library, and is gated on the achievable-set rank being >= 2.
3. **R2 standalone (exact contact calculus)** — only if R1's algorithmic half is absorbed elsewhere,
   and only with B18 proved; otherwise it is a faster evaluator of a quantity that already has one.
4. **R5 (orientation objective)** — cheap, unclaimed, but it is a section, not a submission.
5. **B20 (fully closed = root coincidence)** — the cleanest remaining single section, best-screened
   of the unscreened items (S7 NOT FOUND on all three parts).

**What I would not spend a day on, and why.** B9 in every form (2025 Fig. 10 published the
phenomenon with a photograph and 2026 §6 named the method). B6's algorithm (Hadlock 1975). B15
(routine semicontinuity). B23 (a nicer basis). B16 and B24 as standalones (C4 fails at scale, and the
SDP-in-C++ requirement under directive D3 is the largest implementation risk in any of the files).

---

## 7. What I could not verify

- **I ran no code and no searches.** Every quantity I cite comes from `STATE.md` F15–F19,
  `results/core_validation/rank_claim.md`, `results/core_validation/reference_cases.md`, or a
  persona's own reported measurement, which I have not reproduced.
- **`notes/screen_bundle.md` S1–S7 are complete but the summary table and the query log are still
  `_pending_`.** I used the seven verdicts and the quoted sources; I did not fetch any of them.
- **The three field tables.** I read `notes/field_tutte.md` in the parts bearing on injectivity
  (rows 1–5, 7, 9, 10, 14, 16, 18 and routes A/B/C) and relied on the personas' citations for
  `notes/field_kirigami.md` and `notes/field_rigidity.md` rather than re-reading all 51 rows.
- **The authors' implementation.** `github.com/segaviv/tuttekiri` was not retrieved by anyone in this
  project. Every "beats the baseline" claim in R1, R2 and R5 is currently against the reconstruction
  documented in `code/README.md` deviations 6–9.
- **Dang, Feng, Duan & Wang 2021 (Phys. Rev. E 104, 055006).** The adversary could not extract it and
  neither did I. It matters for B9/A8, which I drop on other grounds (2025 Fig. 10), so the gap no
  longer blocks a decision — but if anyone revives a non-uniform-angle idea, that paper must be read.
- **Whether `p_valid > 0` in K1a.** My entire first recommendation rests on the claim that the usable
  region question has two publishable answers. That the *sampling* distinguishes them is a design
  choice I made and did not test; if `10^4` Gaussian samples are simply the wrong probe of a
  1000-dimensional space (F19: `dim_null` up to 1017), K1a returns `p_valid = 0` everywhere for
  reasons of measure, not of geometry, and the emptiness branch would be a false positive. **The
  Experimenter must therefore also report `p_valid` for a local probe** — a short trust-region walk
  from `X_ini` projected onto `X` — and treat a disagreement between the two probes as inconclusive
  rather than as the emptiness branch. This is the single weakest link in my recommendation.

*End of file.*
