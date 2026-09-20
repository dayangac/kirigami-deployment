# `baseline/` — the authors' published code this project builds on

This folder holds the reference implementations that the rest of the repository is
measured against. It is third-party code; nothing in it was written here.

| folder | upstream | what it is |
|---|---|---|
| `tuttekiri/` | github.com/segaviv/tuttekiri | Segall, Ren, Sorkine-Hornung, *Uniformly Deployable Kirigami on Arbitrary Planar Graphs*, ACM TOG (SIGGRAPH) 2026. C++ core (`code/cpp/`), the web UI sources (`code/src/`, kept because `state.js` holds the `prevent` defaults quoted in `parity.md`), and the pattern data (`code/data/`). |
| `kirigami_tessellations/` | github.com/segaviv/kirigami_tessellations | Segall, Ren, Padilla, Sorkine-Hornung, *Kirigami Tessellations*, SIGGRAPH Asia 2025. C++ core (`code/cpp/`), JS sources (`code/js/`), pattern data (`code/data/`) and the five fabrication SVGs (`fabrication_patterns/`, read by `apps/kill_regime.jl`). |
| `native/` | ours | a CMake harness that compiles the tuttekiri C++ natively and adds a CLI (below). |
| `parity.md` | ours | the measured 1:1 comparison of our `Kirigami/src/core` against their code. |

The two upstream trees are **copies, trimmed to source and data**: their figures, demo
videos, the compiled `.wasm`/`.js` bundle and web-page boilerplate were removed, so the
image links in their own `README.md` files are dead. No source or data file was
modified; the two upstream defects found while building (listed at the end of this file)
are worked around in `native/`, never patched in place. For the full repositories,
licences and the web demos, go to the upstream links above.

The rest of this file documents the native build.

---

## Native build of the authors' tuttekiri code

A referee baseline: the **authors' own C++**, compiled natively on arm64 macOS with
Apple clang, driven on **our** JSON graphs and writing **our** JSON conventions, so that
every number in `code/src/core` can be checked 1:1 against the reference implementation.

* `baseline/tuttekiri/` — upstream checkout, **read-only**. Nothing in this repository
  writes into it or patches it.
* `baseline/native/` — the CMake project that compiles those sources and adds a CLI.
* `baseline/parity.md` — the measured comparison.

## Build

```
cmake -S baseline/native -B baseline/native/build \
      -DCMAKE_OSX_ARCHITECTURES=arm64 -DCMAKE_BUILD_TYPE=Release
cmake --build baseline/native/build -j8
```

Verified: configure 44 s, build ~4 min, `tuttekiri_cli` links and runs. Dependencies are
cloned by `FetchContent` from the public repositories the authors pin (nothing is uploaded):

| dependency | version | note |
|---|---|---|
| libigl | `v2.5.0` | the authors' pin, with `LIBIGL_RESTRICTED_TRIANGLE=ON` |
| Optiz | `master` (github.com/segaviv/optiz) | the authors' pin |
| Eigen | `3.4.0` via FetchContent | **not** Homebrew's 5.0.1, see below |

### Three things had to be worked around, none of them inside `baseline/tuttekiri/`

1. **Eigen.** Homebrew ships Eigen 5.0.1 in `/opt/homebrew/include/eigen3`. libigl v2.5.0
   and Optiz are written against the Eigen 3.x internal API, so 5.0.1 is not a candidate;
   the build pins the same 3.4.0 that the authors' own `CMakeLists.txt` falls back to when
   `find_package(Eigen3)` fails. This keeps the numerics identical to theirs rather than
   merely close.
2. **Optiz and OpenMP.** `optiz/CMakeLists.txt` appends a bare `-fopenmp` on any Apple
   platform, prefixing `-Xclang` only when `CMAKE_CXX_COMPILER_ID STREQUAL "Clang"`. Apple's
   compiler reports `AppleClang`, so the prefix is skipped and the bare flag is rejected.
   Optiz's own `find_package(OpenMP)` had already failed on this machine (no `libomp`), so
   the intended configuration is no-OpenMP; a `PATCH_COMMAND` in our `FetchContent_Declare`
   deletes those two lines **in the fetched copy under `build/_deps`**. Single-threaded
   Optiz is also bit-reproducible, which is what a baseline wants.
3. **emscripten.** `utils/Hmesh.h`, `utils/conversions.h` and `geometry/unit_pattern.h`
   `#include <emscripten/bind.h>`, and `conversions.h` has non-template inline functions
   that use `emscripten::val`. `baseline/native/emsc_stub/emscripten/bind.h` is a ~40-line
   stub declaring a `val` whose every operation **throws**. The native CLI never calls the
   JS glue, so the stub is only ever parsed, and a silent wrong answer is impossible.

`bind.cpp` itself is excluded from the build and replaced by `src/main.cpp`.

## What is built, and what is not

Built and exercised: `coloring/` (face-orientation assignment), `geometry/kirigami.cpp`
(forward deployment, `max_opening_angle`), `geometry/deployment.cpp` (`make_deployable`:
hole/constraint construction, the KKT projection, the null-space basis),
`geometry/unit_pattern.cpp`, `opt/prevent_intersections.cpp` (their Eq. 9),
`opt/fully_closed.cpp`, `opt/opt.cpp`, `parameterization/`, `utils/`.

Not exercised: the inverse-design path (`opt::init` / `optimize_rigidity` / `param::lift`)
needs a 3D target mesh, which is outside this task; it compiles and links but is not
called. Nothing had to be dropped for a hard emscripten dependency.

## Source files added (all under `baseline/native/src/`)

| file | contents |
|---|---|
| `json_io.{h,cpp}` | our JSON graph contract <-> `kirigami::UnitPattern`, and the sigma <-> colour map |
| `bind_port.{h,cpp}` | `do_polygons_intersect` and `max_opening_angle` copied **verbatim** from `bind.cpp` (lines 279-418), plus one diagnostic that reports *which* face pair collides |
| `their_system.{h,cpp}` | a transcription of steps 1-5 of `make_deployable`, used only to measure residuals of *other* embeddings in the authors' own system; validated at run time (see below) |
| `main.cpp` | the CLI |

`make_deployable` returns only the projected positions and the kernel, not its constraint
matrix, so cross-residuals need a rebuild of it. Every run checks that the rebuild
reproduces the authors' own printed `num_holes` / `num_constraints` **and** that the
authors' own returned solution satisfies it to `1e-6`; `rebuild_validated` was `true` on
all eight reference cases.

## The sigma <-> colour correspondence (measured, not assumed)

Our contract: `sigma = +1` clockwise, `sigma = -1` counter-clockwise, hinge at the source.
Theirs: `face_colors` in `{0,1}`. For a hinge edge their `kirigami::deploy` pivots about
the **source of the colour-0 face's stored half-edge**, while ours pivots about the source
of the `sigma = -1` face's stored half-edge — opposite ends of the same edge. The mapping
that makes the two forward-kinematics maps agree is therefore

```
sigma = +1 (clockwise)          <->  colour 1
sigma = -1 (counter-clockwise)  <->  colour 0
```

This was established by measurement: with the naive map the deployed configurations differ
by 4.5% of the patch diameter; with this map they agree to `2e-8` on every reference case
(the residual is `float` rounding — their `deploy` takes `float angle`). Their final
normalisation rotates the deployed patch by the opposite sign to ours, which rigid
alignment absorbs. Faces stored clockwise are reversed on load and the count is reported.

## Command line

```
tuttekiri_cli <cmd> <input> [options]

color     <g.json>                       their coloring::initialized_two_face_coloring
solve     <g.json>                       their make_deployable: X0 + null-space basis
deploy    <g.json> --theta t             their kirigami::deploy, dumps M'
thetamax  <g.json>                       their max_opening_angle, with and without collisions
analyze   <g.json> [--dump d] [--collisions] [--theta t]
                                         the whole pipeline as one JSON record
prevent   <g.json> [--barrier b --strength s --close w] [--sweep] [--dump f]
                                         their opt::prevent_intersections (Eq. 9)
closed    <g.json>                       their opt::optimize_for_fully_closed
parity    <case_dir> [--ourscut cut.json] [--theta t]
                                         one row of baseline/parity.md
collide   <g.json>                       diagnostic: where their collision scan fires
compare   --a x.json --b y.json          rigid Procrustes of two deployed dumps
```

Defaults for `prevent` are the web UI's (`src/state.js`): `barrier 0.1`,
`barrier_strength 10`, `close_to_original_weight 0.1`. `--sweep` runs the 84-point grid
`barrier in {0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0}` x `strength in {1, 10, 50, 100}` x
`close in {0, 0.1, 1.0}`.

## Two upstream defects found while building this

1. **`UnitPattern::get_holes()` segfaults on any bounded patch.** In the
   `face_colors == 1` branch (`unit_pattern.cpp`) it evaluates
   `cur_edge->prev()->twin()->fi` with no null check, so it crashes as soon as a boundary
   half-edge is reached. It is only safe on a fully periodic unit pattern. The CLI guards
   the call and reports `n_traced_holes = -1` instead; all eight of our reference cases have
   a boundary, so their hole tracer could not be run on any of them.
2. **`max_opening_angle(mesh, true)` reports false collisions at small angles.** It calls
   `merge_close_verts()` (tolerance `0.1 * avg_edge_len`) on the deployed mesh before the
   overlap test, which fuses hinge duplicates that have opened by less than that tolerance.
   On `snub_square_33434` this makes their scan return `theta_max = 0.0367 rad` when a scan
   that steps past the artifact finds first contact at `1.69-1.75 rad`. See `parity.md`.

## Reproducing `parity.md`

```
for d in results/core_validation/cases/*/; do c=$(basename $d)
  ./code/build/kiri_analyze $d/M.json --orient json --out /tmp/ours/$c
  ./baseline/native/build/tuttekiri_cli parity $d --theta 0.5235987755982988 \
      --collisions --ourscut /tmp/ours/$c/cut.json --out /tmp/parity/$c.json
done
```
