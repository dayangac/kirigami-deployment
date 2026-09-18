# Hero example 2: K9c Delaunay graph id 130, sigma_mc

Replacement/addition hero picked from K9c's top-10 by ε_max (`results/kill/k9c/k9c.csv`),
restricted to the three geometrically-comparable candidates the orchestrator named:
`delaunay_130` (F=130), `delaunay_28` (F=124), `delaunay_55` (F=120), all `sigma_mc`, all
with `k9c_eps = π` (fully open, `eps_max` capped at π). Selection rule: among candidates
with `eps_max >= 2.5 rad`, pick the one with the largest minimum face-corner angle at the
closed state (a geometric-cleanliness tiebreaker on top of the range).

## Candidate comparison (closed state, `characterize` (`Kirigami/src/method/design.jl`))

| id | F | min corner angle | slivers (< 10°) | eps_max |
|---|---|---|---|---|
| delaunay_28 | 124 | 0.582° | 8 | π |
| delaunay_55 | 120 | 1.846° | 4 | π |
| **delaunay_130** | **130** | **5.069°** | **3** | **π** |

`delaunay_130` wins: largest minimum corner angle and fewest slivers among the three, all
tied at the maximal certified range. **This is hero2.**

## Why kiri_design cannot reproduce this point (and what was done instead)

`kiri_design --maximise-eps` runs `design_constrained` (`method/design.jl`), whose stage-1
embedding is `convex_embed`'s **proximity-to-X_ini** objective (the K9 arm), and only then
conditionally runs `range_opt` stage B from *that* point. K9c's winning point for these
three designs (`best_src = k9c/x0+B` in `k9c.csv`) instead comes from `range_embed`'s
**stage-A margin-maximisation** (`method/range_embed.jl`), warm-started from **t = 0**
(`"x0"`), not from the proximity point — a different starting point feeding a different
objective, so `kiri_design` cannot be relied on to land in the same basin.

This was checked empirically, not just argued: a throwaway program
(`dump_k9c_graph.jl`, not committed) that replicates `kill_k9c.jl`'s exact per-graph
pipeline line-for-line (same `convex_embed`/`range_embed`/`maximize_margin_range` calls,
same seed formula `9300 + 7*id + which`) was run standalone on `delaunay_130
sigma_mc`. It did **not** reproduce the archived point: it found `theta_exact = 2.68985`
via `k9c/k9+B`, versus the archived `theta_exact = π` via `k9c/x0+B` — a different local
optimum of the same non-convex stage-A/stage-B search, most likely from non-associative
floating-point reduction order in the L-BFGS solver compounding over hundreds of
iterations. **Bit-identical re-optimisation is not reliable here, unlike hero's exact
constrained solve.**

So, per the fallback instruction, the export uses the **archived winning embedding
directly**, recovered from the gallery dump of the K9c pass this hero was picked from
(the canonical `results/kill/k9c/k9c.csv`, a later pass of the same driver, ends at a
different local optimum on this design, `theta_exact = 0.572` via the same `k9c/x0+B` route,
and its gallery JSON is that point; the hero geometry survives as `hero2_graph.json` and
the `hero2_130_sigma_mc_*.json` dumps here, and is re-characterised from them).
That file is a `deployment_json` dump (`c.prime_faces` topology — one vertex per face
corner, `theta=0`, `max_mismatch=0`), not the `vertices`/`faces`/`orientation` mesh
contract `kiri_export` reads, so a second throwaway program
(`reconstruct_k9c_graph.jl`) rebuilds the *same deterministic, non-random* topology
(`make_graph` + `make_cut`, no optimisation, no RNG) and folds the gallery's per-corner
positions back onto the original N=71-vertex indexing: for every face `f`, corner `k`,
`X_orig[faces[f][k]] := gallery.vertices[prime_faces[f][k]]`. All prime-copies of the
same original vertex agreed to within `5e-15` (float noise), confirming the fold is exact.
A third throwaway program (`characterize_k9c.jl`) then calls `characterize` — the
same pure, deterministic function `kiri_design` uses to fill its `characterization` field —
directly on the reconstructed graph, independently reproducing `k9c.csv`'s archived
`theta_exact = eps_max = π` to full precision. None of the three throwaway programs were
added to the package; they are graph-topology/measurement plumbing, not new method
code.

## Numbers (`hero2_characterization.json`, `characterize`)

| quantity | value |
|---|---|
| F (faces) | 130 |
| N (vertices) | 71 |
| \|E_split\| (`n_split`) | 18 |
| n_hinge (hinge sites) | 172 |
| Theta_max (T4.2'' exact) | 3.141592653589793 rad (= π) |
| eps_max | 3.1415926535887495 rad |
| theta_max (bisection referee, `--referee`) | 3.141592653589793 rad |
| certified | yes |
| n_inverted | 0 |
| first contact | face_e=3, face_v=56 at theta=3.1415926535897496 |
| min face-corner angle (closed) | 5.068533794889171° |
| sliver count (< 10°, closed) | 3 / 390 corners |

Certificate clauses (`POS /\ NOOVERLAP(eps/2) /\ NOROOT(eps)`, eps=0.3 rad):
`pos = true`, `nooverlap = true`, `noroot = true`, `valid = true`,
`min_signed_area = 1.390707882217665`, `n_candidates = 20696`, `n_pairs = 1188`,
`n_roots_deflated = 0`, `n_roots_inadmissible = 913`, `first_root = -1.0` (no root found
in range), `theta_1 = 0.15`. Matches `results/kill/k9c/k9c.csv` row
`id=130,kind=delaunay,sigma=sigma_mc` (`F=130`, `k9c_eps=k9c_theta=3.14159`,
`best_src=k9c/x0+B`) to the CSV's printed precision.

## Comparison: hero (K9, id 148) vs hero2 (K9c, id 130)

| quantity | hero (K9, id 148) | hero2 (K9c, id 130) |
|---|---|---|
| method | K9 (proximity constrained embed) | K9c (range-maximising: stage A margin max + stage B Theta_max max) |
| F (faces) | 101 | 130 |
| N (vertices) | 59 | 71 |
| \|E_split\| | 16 | 18 |
| Theta_max (exact) | 1.9967778150149833 rad | 3.141592653589793 rad (= π, fully open) |
| eps_max (certified) | 1.9967778150139832 rad | 3.1415926535887495 rad |
| min face-corner angle (closed) | 7.047723201854175° | 5.068533794889171° |
| sliver count (< 10°, closed) | 1 / 303 corners | 3 / 390 corners |
| certified | yes | yes |
| n_inverted | 0 | 0 |

hero2 has 57% more range (fully open vs 2.0 rad) at the cost of a somewhat tighter minimum
corner angle (5.1° vs 7.0°) and more sliver corners (3 vs 1) — the geometric-quality vs.
deployment-range trade-off K9c's summary flagged as a caveat for its top designs.

## Reproduction

1. **Dump the exact K9c-population graph** (topology + sigma_mc, bit-identical to what
   `apps/kill_k9c.jl` builds for `id=130`), same throwaway-utility pattern as hero's
   `dump_hero_graph.jl`:
   ```
   julia --project=Kirigami export/hero2/dump_hero2_input.jl 130 800 export/hero2/hero2_input.json
   ```
   (`dump_hero2_input.jl` calls `make_graph(130, 100, 800, 1400)`, K9c's
   population parameters, and `save_mesh_json`.) Output: V=71, F=130, kind=delaunay.

2. **Recover the archived K9c winning embedding** from the gallery dump (kiri_design
   cannot reproduce it — see above):
   ```
   julia --project=Kirigami export/hero2/reconstruct_k9c_graph.jl 130 800 \
     results/kill/k9c/gallery/delaunay_130_sigma_mc_closed.json export/hero2/hero2_graph.json
   ```
   Prints `gallery theta=0 max_mismatch=0` and `max spread across prime-copies ... 5.02e-15`
   (the fold-back sanity check). Output `hero2_graph.json` is a plain graph JSON
   (`vertices`/`faces`/`orientation`) at the certified embedding, the same contract
   `hero_graph.json` used — `kiri_export` reads it directly.

3. **Independently verify the characterization** (does not depend on step 2's provenance,
   only on the recovered X and `characterize`, the same function `kiri_design`
   uses):
   ```
   julia --project=Kirigami export/hero2/characterize_k9c.jl export/hero2/hero2_graph.json export/hero2/hero2_characterization.json
   ```
   Prints `theta_max=3.14159 eps_max=3.14159 theta_bisect=3.14159 certified=1
   n_inverted=0 min_corner_angle_deg=5.06853 n_sliver=3` — the archived winning point
   (the canonical `k9c.csv` ends at 0.572 rad on this design, see above).

4. **Export.**
   ```
   EXP="julia --project=Kirigami Kirigami/apps/kiri_export.jl"
   $EXP export/hero2/hero2_graph.json --profile felt_laser --theta 0 \
     --svg export/hero2/hero2_130_sigma_mc_closed.svg \
     --stl export/hero2/hero2_130_sigma_mc_closed.stl \
     --3mf export/hero2/hero2_130_sigma_mc_closed.3mf \
     --json export/hero2/hero2_130_sigma_mc_closed.json

   $EXP export/hero2/hero2_graph.json --profile felt_laser --theta-frac 0.5 \
     --svg export/hero2/hero2_130_sigma_mc_open_half.svg \
     --stl export/hero2/hero2_130_sigma_mc_open_half.stl \
     --3mf export/hero2/hero2_130_sigma_mc_open_half.3mf \
     --json export/hero2/hero2_130_sigma_mc_open_half.json

   $EXP export/hero2/hero2_graph.json --profile felt_laser --theta-frac 0.9 \
     --svg export/hero2/hero2_130_sigma_mc_open_0.9tm.svg \
     --stl export/hero2/hero2_130_sigma_mc_open_0.9tm.stl \
     --3mf export/hero2/hero2_130_sigma_mc_open_0.9tm.3mf \
     --json export/hero2/hero2_130_sigma_mc_open_0.9tm.json
   ```
   `--theta-frac` reports `theta_max=3.14159` at every fraction (export/solid.jl's
   independent geometric cross-check), agreeing with the certified value.

   Deployment sequence, 6 frames at `theta = k/5 * Theta_max`, `k = 0..5`:
   ```
   for k in 0 1 2 3 4 5; do
     frac=$(echo "scale=1; $k/5" | bc)
     $EXP export/hero2/hero2_graph.json --profile felt_laser --theta-frac $frac \
       --svg export/hero2/sequence/hero2_130_seq_k${k}.svg \
       --json export/hero2/sequence/hero2_130_seq_k${k}.json
   done
   ```

5. **PNG previews:**
   ```
   python3 export/render_svg.py <svg> -o <png>
   ```
   run once per SVG above (closed, open_half, open_0.9tm, and the 6 sequence frames).

## Validation

- **XML well-formedness**: `xmllint --noout` on all 9 SVGs (closed, open_half,
  open_0.9tm, 6 sequence frames) — all pass.
- **3MF zip integrity**: `unzip -t` on all 3 3MF files — `No errors detected in
  compressed data` for all three.
- **kiri_export's own checks** (`check_manifold`, printed per run): all three solids
  report `closed, consistently oriented, boundary edges 0, non-manifold edges 0, flipped
  edges 0`; 1 component at theta=0 (single living-hinge sheet), 130 components once open
  (per-face bodies separated, as documented for the living-hinge writer in
  `Kirigami/README.md`'s `src/export/` section — a rendering/kinematics fact, not a solid-quality
  defect, since each per-face body is individually closed and consistently oriented).
  `degenerate tris` (234-256 depending on theta) are the flat ears documented in
  `Kirigami/README.md`'s `src/export/` section (flat ears from `build_solid`); they carry no area or
  volume. Ear-clipping fell back to a fan on 3-4 faces per state (faces 6, 53, 56, 116) —
  these are the thinnest slivers (`min_corner_angle` 5.07°), consistent with the sliver
  count above; the fallback still closes the solid (confirmed by the manifold check).
- **Visual check**: PNG previews of the closed state and both open states
  (`hero2_130_sigma_mc_closed.png`, `_open_half.png`, `_open_0.9tm.png`) were viewed.
  Closed: 130 faces tile edge-to-edge with no gaps or overlaps (a few very thin sliver
  triangles are visible near vertices 6/7/41/96/97, consistent with the 3-sliver count).
  Open at both 0.5 and 0.9 * Theta_max: faces have visibly separated (holes open at every
  hinge/split site) and no two faces overlap in the drawing, consistent with the
  certified `nooverlap` clause.

## Files

```
export/hero2/
  hero2_input.json                        # step 1: topology + sigma_mc (make_graph(130,...))
  hero2_graph.json                        # step 2: reconstructed certified embedding, plain graph JSON
  hero2_characterization.json             # step 3: characterize full report
  hero2_130_sigma_mc_closed.{svg,stl,3mf,json,png}      # theta = 0
  hero2_130_sigma_mc_open_half.{svg,stl,3mf,json,png}   # theta = Theta_max / 2
  hero2_130_sigma_mc_open_0.9tm.{svg,stl,3mf,json,png}  # theta = 0.9 * Theta_max
  sequence/
    hero2_130_seq_k{0..5}.{svg,json,png}  # theta = k * Theta_max / 5, k = 0..5 (deployment figure)
```

`dump_hero2_input.jl`, `reconstruct_k9c_graph.jl`, `characterize_k9c.jl` (the
throwaway utilities above) were run from the scratch directory, not added to
the repository or the package — they are graph-topology dump / fold-back /
measurement plumbing around existing library code (`make_graph` (`Kirigami/src/core/kill_common.jl`),
`make_cut`, `characterize`), not new method code or a deliverable.
`export/hero/` is unmodified.
