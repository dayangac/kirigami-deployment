# Hero example: K9 Delaunay graph id 148, sigma_mc

The hero design chosen by the orchestrator: K9's population graph id 148
(`delaunay` kind), the Eq. (1) max-cut orientation (`sigma_mc`), constrained
embedding. This is one of K9's 36/400 positives (the smallest of the two
regression-locked in `code/tests/test_design.cpp`), certified valid with the
largest `Theta_max` in that set.

## Numbers (from `hero_design.json`, `characterization`)

| quantity | value |
|---|---|
| F (faces) | 101 |
| N (vertices) | 59 |
| \|E_split\| (`n_split`) | 16 |
| dim_null | 16 |
| n_hinge (hinge sites) | 128 |
| Theta_max (T4.2'' exact) | 1.9967778150149833 rad |
| eps_max | 1.9967778150139832 rad |
| theta_max (bisection referee, `--referee`) | 1.9967778152171005 rad |
| certified | yes |
| n_inverted | 0 |
| first contact | face_e=53, face_v=52 at theta=1.9967778150149833 |

Certificate clauses (`POS /\ NOOVERLAP(eps/2) /\ NOROOT(eps)`, eps=0.3 rad):
`pos = true`, `nooverlap = true`, `noroot = true`, `valid = true`,
`min_signed_area = 0.9709894601265123`, `n_candidates = 13096`,
`n_pairs = 756`, `n_roots_deflated = 0`, `n_roots_inadmissible = 544`,
`first_root = -1.0` (no root found in range), `theta_1 = 0.15`.

These match `k9.csv` row `id=148,kind=delaunay,sigma=sigma_mc` (columns
`F=101`, `a_theta_exact=a_eps_max=1.99678`, `a_feasible=1`) and the
`tests/test_design.cpp` regression lock (`theta_max == eps_max ==
1.9967778150149833`, epsilon 1e-9) to full precision, i.e. well inside the
requested 1e-4 tolerance.

## Reproduction

1. **Dump the exact K9-population graph** (topology + `assign_orientation_relaxation`
   sigma + `X_ini`, bit-identical to what `apps/kill_k9.cpp` and
   `tests/test_design.cpp`'s `k9_design(148)` build — this matters because
   `results/kill/k5/sigma/delaunay_148.json` on disk holds K5's *sigma_def*
   orientation, not sigma_mc, and `kiri_design --sigma mc`'s own
   `orientation_maxcut` uses a different RNG seed than `make_graph`'s internal
   one, so neither reproduces the locked value on its own). A throwaway app
   (not added to the CMake build, to avoid touching the shared `CMakeLists.txt`
   mid-flight) does this:

   ```
   clang++ -std=c++20 -O2 -arch arm64 \
     -I code/src -I code/apps -I /opt/homebrew/include -I /opt/homebrew/include/eigen3 \
     dump_hero_graph.cpp code/build/libkiri_core.a -o dump_hero_graph
   ./dump_hero_graph 148 export/hero/hero_input.json
   ```
   `dump_hero_graph.cpp` just calls `kiri::kill::make_graph(148, 100, 800, 1400)`
   (K9's population parameters) and `save_mesh_json`. Output: V=59, F=101,
   kind=delaunay.

2. **Regenerate the design** with the method API, `--sigma json` so the CLI
   uses the sigma baked into step 1's file instead of recomputing it, and the
   K9 seed formula (`9000 + 7*id`, `which=0` for sigma_mc):

   ```
   ./code/build/kiri_design export/hero/hero_input.json --sigma json --seed 10036 \
     --out export/hero/hero_design.json --referee --out-graph export/hero/hero_graph.json
   ```
   Output: `sigma_json constrained: Theta_max = 1.996778 rad, eps_max = 1.996778 rad,
   certified = yes, feasible = yes ... inverted faces 0` — matches the target
   to the printed precision; full precision confirmed in `hero_design.json`
   above.

3. **Export.** `hero_graph.json` (the `--out-graph` output) is a plain graph
   JSON (`vertices`/`faces`/`orientation`) at the certified embedding, which
   `kiri_export` reads directly:

   ```
   EXP=code/build/src/export/kiri_export
   $EXP export/hero/hero_graph.json --profile felt_laser --theta 0 \
     --svg export/hero/hero_148_sigma_mc_closed.svg \
     --stl export/hero/hero_148_sigma_mc_closed.stl \
     --3mf export/hero/hero_148_sigma_mc_closed.3mf \
     --json export/hero/hero_148_sigma_mc_closed.json

   $EXP export/hero/hero_graph.json --profile felt_laser --theta-frac 0.5 \
     --svg export/hero/hero_148_sigma_mc_open_half.svg \
     --stl export/hero/hero_148_sigma_mc_open_half.stl \
     --3mf export/hero/hero_148_sigma_mc_open_half.3mf \
     --json export/hero/hero_148_sigma_mc_open_half.json

   $EXP export/hero/hero_graph.json --profile felt_laser --theta-frac 0.9 \
     --svg export/hero/hero_148_sigma_mc_open_0.9tm.svg \
     --stl export/hero/hero_148_sigma_mc_open_0.9tm.stl \
     --3mf export/hero/hero_148_sigma_mc_open_0.9tm.3mf \
     --json export/hero/hero_148_sigma_mc_open_0.9tm.json
   ```
   `--theta-frac` computes `theta_max` internally (`export/solid.hpp`'s
   geometric `theta_max`, a cross-check independent of `method::characterize`)
   and reports `theta_max=1.99678` at every fraction — agrees with the design's
   certified value.

   Deployment sequence for the paper figure, 6 frames at `theta = k/5 * Theta_max`,
   `k = 0..5` (k=5 is exactly `Theta_max`, the first-contact configuration):
   ```
   for k in 0 1 2 3 4 5; do
     frac=$(python3 -c "print($k/5)")
     $EXP export/hero/hero_graph.json --profile felt_laser --theta-frac $frac \
       --svg export/hero/sequence/hero_148_seq_k${k}.svg \
       --json export/hero/sequence/hero_148_seq_k${k}.json
   done
   ```

4. **PNG previews:**
   ```
   arch -arm64 /usr/local/bin/python3 export/render_svg.py <svg> -o <png>
   ```
   run once per SVG above (closed, open_half, open_0.9tm, and the 6 sequence
   frames).

## Validation

- **XML well-formedness**: `xmllint --noout` on all 9 SVGs (closed, open_half,
  open_0.9tm, 6 sequence frames) — all pass.
- **3MF zip integrity**: `unzip -t` on all 3 3MF files (closed, open_half,
  open_0.9tm) — `No errors detected in compressed data` for all three.
- **kiri_export's own checks** (`check_manifold`, printed per run): all three
  solids report `closed, consistently oriented, boundary edges 0, non-manifold
  edges 0, flipped edges 0`, 1 component at theta=0 (single living-hinge
  sheet) and 101 components once open (necks visibly separated at theta>0, as
  documented for the living-hinge writer in `code/README.md`'s Export
  section) — this is a rendering/kinematics fact, not a solid-quality defect,
  since each of the 101 per-face bodies is individually closed and
  consistently oriented. `degenerate tris` (174-197 depending on theta) are
  the flat ears documented in `code/README.md`'s "Deviations and known limits
  (export)" item 2; they carry no area or volume.
- **Visual check**: PNG previews of the closed state and both open states
  (`hero_148_sigma_mc_closed.png`, `_open_half.png`, `_open_0.9tm.png`) were
  viewed. Closed: 101 faces tile edge-to-edge with no gaps or overlaps. Open
  at both 0.5 and 0.9 * Theta_max: faces have visibly separated (holes open
  at every hinge/split site) and no two faces overlap in the drawing,
  consistent with the certified `nooverlap` clause.

## Files

```
export/hero/
  hero_input.json                       # step 1: topology + sigma_mc + X_ini (make_graph(148,...))
  hero_design.json                      # step 2: kiri_design full report (X0, characterization, certificate)
  hero_graph.json                       # step 2: --out-graph, the certified embedding as plain graph JSON
  hero_148_sigma_mc_closed.{svg,stl,3mf,json,png}      # theta = 0
  hero_148_sigma_mc_open_half.{svg,stl,3mf,json,png}   # theta = Theta_max / 2
  hero_148_sigma_mc_open_0.9tm.{svg,stl,3mf,json,png}  # theta = 0.9 * Theta_max
  sequence/
    hero_148_seq_k{0..5}.{svg,json,png}  # theta = k * Theta_max / 5, k = 0..5 (deployment figure)
```

`dump_hero_graph.cpp` (the throwaway graph-dump utility from step 1) was
compiled and run from the scratch directory, not added to the repository or
the CMake build — it is a 20-line wrapper around
`kiri::kill::make_graph` + `save_mesh_json`, not a deliverable.
