# Kirigami -- the usable shape space of Tutte auxetic kirigami

Julia implementation of the kirigami deployment pipeline of Segall, Ren and
Sorkine-Hornung (TOG 2026), together with this project's exact deployment-range
characterization and its range-maximising constrained embedding. The package is
`Kirigami/`; the numerical conventions that keep its results reproducible against the
frozen corpora are in `docs/NUMERICS.md`.

## What the method is

Segall et al. give a linear shape space `X` of uniformly deployable kirigami on any
oriented planar graph and leave open whether a member of `X` actually opens without
self-intersection (`../notes/IDEA.md` §1, `../notes/reading/paper_2026.md`). This project shows that every
geometric predicate along the uniform deployment branch is a single first harmonic in the
opening angle whose coefficients are quadratic in the flat design coordinates, so the
collision-free range `Θ_max` is exactly computable from a complete closed-form candidate
list with grazing contacts handled explicitly, and a sound sufficient certificate
`Θ_max ≥ ε` follows (`../notes/derivations/core.md` T3-T5, `../paper/REPORT.md` C1-C3). It then shows that
the paper's own Eq. (6) projection is certifiably non-deployable on random planar graphs
because it manufactures reflex face corners that jam the mechanism at `θ = 0⁺`, and that
replacing the projection's proximity objective by the `0⁺` deployment margin inside
convexity and split-outward barriers in the same null space makes most of the same
population deployable and certified (`../paper/REPORT.md` C7-C9, `results/kill/KILL_REPORT.md`,
`results/kill/k9c/`; measured yields: 312 of 400 deployable, 310 certified).

## Repository layout

| path | contents |
|---|---|
| `Kirigami/` | the Julia package (`Project.toml`, one module `Kirigami`, nothing exported; see `Kirigami/README.md`) |
| `Kirigami/src/{core,method,export}/` | one `.jl` per unit (mesh, cut, holes, Tutte auxetic system, kinematics, collision, generators, RNG; deploy basis, contact calculus, design API, constrained/range-maximising embeddings, mobility, periodic Jacobian, expansive cone, budget; layout, material, SVG/STL/3MF writers) |
| `Kirigami/test/` | the unit tests, one file per unit (`runtests.jl`, `helpers.jl`, `test_*.jl`, `derivation_tests.jl`) |
| `Kirigami/apps/` | user-facing CLIs `kiri_*.jl` and the experiment drivers `kill_*.jl`, plus the included helper files `common_app.jl`, `kill_common.jl`, `native_common.jl` |
| `Kirigami/gui/` | GLMakie desktop explorer (`app.jl`, window-free `model.jl`, own `Project.toml`, headless tests in `gui/test/`) |
| `Kirigami/scripts/` | CairoMakie plotting / summarising scripts `plot_*.jl` (own `Project.toml`) and the two shell drivers `run_e1.sh`, `run_scaling.sh` |
| `data/corpus/` | the frozen input populations of every experiment (`k1a_200.json`, `native200.json`, `k3a_500.json`, ..., `reference_cases_8.json`, `mt19937_vectors.json`), the unit-test fixtures in `reference_patterns/`, and the method fixtures in `method_fixtures/`; see `data/corpus/README.md` |
| `data/deploy_frames.json`, `data/web_patterns.json` | the GUI's pre-characterised designs and the raw pattern library |
| `docs/NUMERICS.md` | the numerical conventions (fma placement, libm, RNG, BLAS threads, frozen corpora) |
| `results/` | measured outputs: `core_validation/` (reference cases), `kill/` (K1a ... K9c, `KILL_REPORT.md`), `final/` (E1 and the paper figures), `regime/`, `scaling/`, `yield/`, `propagation/` |
| `export/` | fabrication samples: `hero/`, `hero2/` (JSON/SVG/STL/3MF/PNG of the two hero designs, with READMEs), `samples/`, `render_svg.py` |
| `baseline/` | the authors' 2025/2026 code (third-party, untouched) and the thin harness `baseline/native/` used as the parity referee (`baseline/parity.md`) |

Not in the repository: the LaTeX sources (`paper/siggraph-poster/`, `paper/techreport/`, `paper/siggraph-poster/`, with `REPORT.md`) and the working notes (`notes/`: `derivations/` with the theory `core.md`, `lemmas.md`, `check.md` and its `scratch/check_*.jl` programs, `ideas/`, `review/`, `specs/`, `reading/`, `STATE.md`, `IDEA.md`) live in the sibling folders `../paper/` and `../notes/` next to this checkout; the papers read figures from `results/` and `export/` here via `\graphicspath`, and the scratch programs load this package via `--project=code/Kirigami`.

## How to run

Everything runs from source; there is no build step. One-off setup:

```
export PATH=$HOME/.juliaup/bin:$PATH          # native arm64 Julia, see below
julia --project=Kirigami -e 'using Pkg; Pkg.instantiate()'
```

**Tests** (187 test sets / 176,692 assertions (23 marked broken), 0 failures; the 23 broken are optimiser-path-dependent `@test_broken` locks, see `test_design.jl` and `docs/NUMERICS.md`):

```
julia --project=Kirigami -e 'using Pkg; Pkg.test()'
julia --project=Kirigami -e 'using Pkg; Pkg.test(test_args=["derivation_tests"])'   # one file (basename substring)
julia --project=Kirigami -e 'using Pkg; Pkg.test(test_args=["design", "export"])'   # several
```

`derivation_tests.jl` is the Checker's independent verification of the derivations
(38 test sets / 150,191 assertions, 0 failures; 151,801 assertions with `KIRIGAMI_FULL_DERIVATIONS=1`); the scratch programs behind `../notes/derivations/check.md` are run one
by one from the parent folder as `julia --project=code/Kirigami notes/derivations/scratch/check_X.jl`
(list in `../notes/derivations/scratch/README.md`).

**A driver.** The user-facing CLIs and the experiment drivers (`-h` prints the flags of
each):

```
julia --project=Kirigami Kirigami/apps/kiri_gen.jl delaunay 60 40 --out g.json --orient auto --seed 3
julia --project=Kirigami Kirigami/apps/kiri_design.jl g.json --sigma both --out design.json --baseline --out-graph best.json
julia --project=Kirigami Kirigami/apps/kiri_analyze.jl best.json --out out/
julia --project=Kirigami Kirigami/apps/kiri_export.jl best.json --profile felt_laser --theta-frac 0.5 --svg out.svg --stl out.stl --3mf out.3mf
julia --project=Kirigami Kirigami/apps/kiri_reference.jl --out results/core_validation
julia --project=Kirigami Kirigami/apps/kill_k9c.jl --out results/kill/k9c/k9c.csv        # an experiment; -h prints its flags
```

`kiri_design` with the run's seed (`--seed $((9300 + 7*id + which))`) is the command that
reproduces one row of `results/kill/k9c/k9c.csv`; see the reproducibility note in
`Kirigami/README.md` for what "reproduces" means numerically.

**The GUI** (separate environment so the core package stays light):

```
julia --project=Kirigami/gui -e 'using Pkg; Pkg.instantiate()'
julia --project=Kirigami/gui Kirigami/gui/app.jl [data/deploy_frames.json] [--theta rad] [--screenshot out.png]
julia --project=Kirigami/gui Kirigami/gui/test/runtests.jl          # headless model tests
```

**Plotting scripts and figures** (CairoMakie, separate environment). They read the CSVs
and JSON under `results/` and never re-run a solver:

```
julia --project=Kirigami/scripts -e 'using Pkg; Pkg.instantiate()'
julia --project=Kirigami/scripts Kirigami/scripts/plot_final.jl     # results/final/figures/fig_*.png, summary_table.md
julia --project=Kirigami/scripts Kirigami/scripts/plot_k9c.jl       # the K9c figures under results/kill/k9c/
```

To regenerate one figure end to end (data and picture), run the driver that produced its
CSV and then the script that reads it; `results/final/figures/README.md` lists, per
figure, which files it reads. Example for the yield figure: `kill_k9c.jl` (400 designs;
hours), then `plot_final.jl`.

## Frozen corpus policy and `--regenerate`

Every input population was generated once and written to `data/corpus/`
(producer, toolchain and per-file provenance in `data/corpus/README.md`). The Julia apps
and tests read those files by default, so an experiment always runs on byte-identical
inputs whichever machine runs it. Each driver that takes a population also accepts
`--regenerate`, which rebuilds it through the generators and the bit-exact
`MT19937` (`Kirigami/apps/kill_common.jl`: `population(name; regenerate)`,
`reference_cases(; regenerate)`, `deployable_population(; regenerate)`). On the reference
platform `make_graph` regenerates all 200 graphs of `k1a_200.json` vertex for vertex; the
frozen file remains the definition of the population. The corpus is never rewritten by a
Julia run.

## Julia binary: use the native arm64 build

The frozen corpora and fixtures were produced on Apple arm64 and their bits depend on
that platform's libm and on FMA placement (`docs/NUMERICS.md`). Use the native arm64
Julia from juliaup (`export PATH=$HOME/.juliaup/bin:$PATH`; the package declares
`julia = "1.10"`, the reference runs were made on 1.12). The x86_64 Homebrew build under
Rosetta is slower and its libm differs by 1 ulp in trig/hypot on a fraction of arguments,
which shows up as equal-objective tie flips in the orientation relaxation on symmetric
tilings (documented in `Kirigami/src/core/generators.jl`); random graphs are unaffected.
