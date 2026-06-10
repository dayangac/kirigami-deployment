# Hostile referee response: does the negative result survive?

Role: co-author of "Uniformly Deployable Kirigami on Arbitrary Planar Graphs" (Segall,
Ren, Sorkine-Hornung, TOG 2026), defending the paper against `results/kill/KILL_REPORT.md`
and `STATE.md`'s claims that its shape space is empty on random planar graphs.

## 1. Does the paper claim it deploys *random* planar graphs?

**No — never, anywhere.** The paper's generality claim (Sec. 5 opener, `notes/paper_2026.md`
lines 763-771) is combinatorial, not statistical:

> "Our analysis applies to *arbitrary* planar graphs, including those that are non-periodic
> (Figures 1 and 3), non-2-colorable (Fig. 1, 9 and 11) or with inner holes (Figures 17 and
> 18)."

"Arbitrary" here means "we don't require 2-colorability / periodicity / disk-without-holes,"
demonstrated by one hand-picked counterexample per exclusion. It is never demonstrated, or
claimed, over a random sample.

- **Fig. 1** (p. 60:1): "Non-two-colorable, aperiodic planar graph... The specific graph is
  not given" (figure table, `notes/paper_2026.md` line 875). An artistic input (panda
  texture), not a random or standard one.
- **Fig. 3** (p. 60:3): a hemisphere 3D mesh fed through the **inverse-design** pipeline
  (Sec. 5.3b) — flatten via Tutte auxetic embedding with boundary fixed to a circle, cut,
  then optimize `Y_theta` to match target edge lengths. This is a *disk-topology mesh*, not
  a random planar graph exercising the Eq.(4)-Eq.(6)-Eq.(9) deployability pipeline the kill
  campaign tests. Conflating the two pipelines would be our own error, not theirs.

The only *named, standard* input anywhere in the paper is the (3,4,3,12) Archimedean tiling
(Fig. 21), and every other reproducible figure (Figs. 4-29) is a small, curated periodic or
symmetric patch, `F` in the range 4-97 in the reference cases (`rank_claim.md`). No figure,
table, or sentence reports a success rate over a *population* of graphs of any kind. The
"random graph" experiment is entirely ours, not theirs, and the paper cannot be shown to
overclaim on a question it never poses.

## 2. Is the "0/400" baseline fair to the published pipeline?

**No — the 0/400 figures are Eq.(6)-projection-plus-our-certificate, not the authors' full
pipeline, and this is not stated clearly enough where the numbers are quoted.**

The authors' full pipeline is color (Eq. 1) → solve (Eq. 4/6) → **prevent** (Eq. 9, their
native `opt::prevent_intersections`) → deploy with **their** collision test. Tallying how
often that full stack, not just our stand-in for step 3, was actually run on random graphs:

| experiment | population | native `prevent` run on | outcome |
|---|---|---|---|
| K1a (`0/200`, `200/200` overlap, `0/200` valid samples) | 200 random graphs | **0** — Eq.(6)+sampling only, no `prevent` | n/a |
| K5 (`0/200` both σ) | same 200 | **0** | n/a |
| K6 (`0/400`, "all four strategies") | same 200 × 2σ | **0** in the headline table; K6's own "native cap" side-run used **8/80** designs (40 graphs × 2σ, capped `F ≤ 550`) | native `theta_max = 0` on 8/8 |
| B4 (`0/400`, boundary-relaxed) | same 200 × 2σ × 3 variants | **0** | n/a |
| K2b (median gain 0, "30 live designs") | 30 designs incl. 8 random graphs | native `prevent --sweep` run on all 30, but the 8 random ones are vacuous (**0/37 ladder points valid** before `prevent` could even act) | n/a for random subset |

So across roughly **1600 (design, σ) slots** derived from the 200-graph random population
(K1a/K5/K6/B4's four strategies × 2σ, plus K2b's random subset), the authors' *native*
Eq.(9) optimizer was actually invoked on **8 designs total** (K6's capped run) — because on
every other slot the pattern was already invalid (self-intersecting or `Θ_max = 0`) before
`prevent` had anything to work with. `prevent` cannot rescue an embedding that is not a
valid flat state to begin with; it optimizes *within* the null space, and K1a already shows
`0/200` samples in that null space are valid embeddings under `σ_mc`.

**What is and is not supported:**

- **Supported, cleanly:** "Eq.(6)'s projection, under the paper's own Eq.(1) orientation, is
  not a valid (non-self-intersecting) embedding on any of 200 random Voronoi/Delaunay/quad
  graphs, and no sample of the null space around it is either" — this is K1a, measured
  directly, no baseline-fairness question attaches to it since Eq.(9)/`prevent` is a
  *downstream* repair of collisions in an already-valid embedding, not a fix for
  self-intersection at `θ=0`.
- **Supported, narrowly:** "on the 8 designs where a valid flat state was reached at all
  (K6's capped native run), the authors' own `prevent`, at published defaults, also failed
  to produce `Θ_max > 0`" — n=8, not 400.
- **Not supported as literally written:** any sentence of the form "the published pipeline
  yields nothing on random graphs" if read as testing the *full* pipeline including native
  `prevent`. I searched `KILL_REPORT.md` and `STATE.md` for that exact framing and did not
  find it asserted outright — the closest is F30/F36's "0/400" language, which is
  Eq.(6)+repair, not Eq.(6)+Eq.(9)-native+repair. It should be labelled "Eq.(6) (and our
  0⁺/convexity repairs) alone" wherever quoted, with the n=8 native spot-check cited
  separately and explicitly as underpowered.

## 3. Is the population representative or adversarial?

**Adversarial in scale, not in generative process.** Voronoi/Delaunay/quad-dominant
triangulations of uniform random points (`code/src/core/generators.cpp:496-580`) are a
standard, unbiased way to produce "arbitrary planar graphs" — nothing is gerrymandered
against the method. But:

- **Size.** The 200-graph population has `F` in `[101, 793]`, median 332.5
  (`KILL_REPORT.md` K5/K6). Every reference case drawn from the paper's own figures has `F`
  in `[4, 97]` (`rank_claim.md`'s reference-case table; `snub_square_33434` at F=53 is the
  largest reproducible one). The random population is **3-8× larger by face count than any
  figure the paper ever demonstrates**, and shape-space dimension (hence the difficulty of
  finding a jointly-feasible embedding) scales with `|E_split|`, which scales with `F`.
- **Boundary/hole geometry.** The authors' own Sec. 4.2 heuristic explicitly wants "many
  balanced, small holes" (`notes/paper_2026.md` §13.2) as a *proxy* for good auxetic
  behavior — an admission that not every `σ` (and by extension not every graph) is expected
  to behave well. Random Voronoi/Delaunay cells are not selected for this property at all;
  no attempt is made to prefer orientations or graphs with small balanced holes beyond the
  max-cut relaxation's blunt edge-count objective.
- **Split-cut density.** A3 (jitter transition) found 4 of the paper's 8 reference tilings
  are **split-free** under their own orientation, i.e. their shape space is a single point —
  the authored tiling itself, trivially deployable. That is not a coincidence: it is easier
  to *author* a 2-colorable tiling than to hit one by chance. The random population's split
  density is whatever Eq.(1)'s max-cut relaxation happens to produce on unstructured graphs,
  which K5 shows correlates only weakly (`-0.03` to `-0.64`) with the deployability defect.

So the honest description is: **a representative sampling procedure, at a scale and
split-density regime the paper never tests**, not a hand-picked adversarial input. Both
things are true and should both be said. "Representative of arbitrary planar graphs" is
fair; "representative of what the paper demonstrates" is not.

## 4. Does Limitation 1 concede the problem statement?

**Yes, substantially, and it should be quoted in full, not clipped to the last sentence.**
Full text (`notes/paper_2026.md` lines 949-951, Sec. 6):

> "First, although all embeddings in the solution space `X` are theoretically uniformly
> deployable, some may be geometrically undesirable (e.g., exhibiting extremely short edges
> or sharp angles). Moreover, uniform deployability does not guarantee a large
> collision-free deployment range; self-intersections may occur at small opening angles
> (e.g., Fig. 14). At present, such issues can only be detected through explicit deployment
> simulation via forward kinematics and require additional post-processing to resolve.
> Developing geometric characteristics for favorable deployment behavior directly from the
> embedding remains an open problem."

This concedes exactly R1's problem statement: (a) membership in `X` does not imply a valid
or useful embedding, (b) the authors have no a-priori certificate, only "explicit deployment
simulation... and post-processing," and (c) finding one is explicitly open. It does **not**
concede that the shape space is *generically empty* of good points, or that `X` is empty on
random graphs specifically — the limitation is agnostic about frequency. The paper should be
quoted with the whole paragraph, and the framing should be "the paper's own Limitation 1
states the open problem this work answers negatively (on this population) and partially
addresses constructively (K9, F36)" — not "the paper admits its method doesn't work."

## 5. Overclaiming audit — sentence and correction

| location | as written | problem | corrected |
|---|---|---|---|
| `KILL_REPORT.md` §K6 "The whole design space, all four strategies, is `0/400` certified." | reads as if the authors' pipeline was exhaustively tried | the "four strategies" are Eq.(6) alone, `σ_def`, and two internal 0⁺-repair variants — **none is the authors' native `prevent`**; that ran on a separate, capped 8-design side-experiment | "All four of our own repair strategies are `0/400` certified; the authors' native `prevent`, tested separately on the 8 of these 400 designs that reached a valid flat state, is also `0/8`." |
| `STATE.md` F30 "on random graphs the Tutte auxetic shape space under BOTH tested orientation rules ... contains no sampled/projected point with positive deployment range (0/400)" | same conflation, plus "the shape space contains no point" over-generalizes from *sampled/projected* points to the space itself | "no *sampled or Eq.(6)/repair-projected* point... has positive range; whether the true infinite-dimensional space contains such a point is not decided by sampling" |
| `STATE.md` line 31, K1a "PASS-emptiness" | terse enough not to overclaim on its own, but downstream quotes of it drop the qualifier "under `σ_mc`" | keep the qualifier every time it is re-quoted: results are conditional on the Eq.(1) max-cut orientation (K5 shows a different `σ` changes the picture on the flat-sheet defect, though not on range) |
| `STATE.md` F27 "RANGE-MARGIN CLAIM DEAD" | flatly correct per K2b, but reads as "our method never wins," while K2b's own table shows 7/30 wins, best +10.6% | "the *median* margin claim is dead (0% vs a 25% bar); the method still wins on 7/30 designs, by up to 10.6%" |
| `STATE.md` F36 "FIRST NON-ZERO CONSTRUCTIVE RESULT ON RANDOM GRAPHS" | accurate and already hedged with the 36/400 vs 40/400 bar miss in the same line — no correction needed | none |

No sentence found that literally claims "the published pipeline achieves nothing on random
graphs" including `prevent`; the risk is *implication* through juxtaposition (a `0/400`
number sitting next to prose about "the authors' pipeline") rather than a false explicit
claim. Every KILL_REPORT verdict I checked (K1a, K1c, K2b, K5, K6, B4, K8a) states its own
population and baseline precisely in its own section; the risk is entirely in *STATE.md's
compressed one-line summaries and in this review's own headline*, not in the underlying
measurements.

## 6. The one experiment that would make this unassailable

**Run the authors' native pipeline — color → Eq.(6)/(4) solve → their `prevent` (Eq. 9,
published defaults + the existing 36-point ladder) → their (or our corrected) collision test
— on all 200 random graphs under both `σ_mc` and `σ_def`, unconditionally, not gated on
0⁺-feasibility first.**

Currently `prevent` is only ever invoked *after* a design has already been filtered to be
0⁺-feasible or otherwise promising (K6's cap, K2b's live-design selection), which is
reasonable for cost but leaves the "what if `prevent` alone, from the raw `X0`, rescues some
of the 142/200 or 200/200 already-invalid projections" question untested at scale. Given
`baseline/parity.md`'s measured per-graph cost of native `prevent` (0.95 s at 60 faces up to
63 s at 300, and 26-160 ms at `F≈50`) and this population's median `F=332.5` with `dim_null`
in the hundreds, a full run at the single published-default parameter point (skip the ladder)
is **200 graphs × 2 σ × ~10-60 s ≈ 1-7 hours** of wall time on one machine, embarrassingly
parallel across the 12-shard infrastructure already built for K1a/K5/K6/B4 — well under an
hour wall-clock with 12-way sharding. This closes the only gap a hostile reviewer can still
exploit: right now the "0/400" headline is honest about what *was* measured but is
vulnerable to "you never actually ran our optimizer on your graphs," and this experiment
removes that objection entirely, in either direction.
