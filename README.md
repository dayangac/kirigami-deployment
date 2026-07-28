# kirigami-julia -- the usable shape space of Tutte auxetic kirigami

Julia reimplementation of the kirigami deployment pipeline of Segall, Ren and
Sorkine-Hornung (TOG 2026), together with this project's exact deployment-range
characterization and its range-maximising constrained embedding. The reference
implementation it reproduces is the C++/Eigen code of `kirigami-experiments` (see
`PORTING.md` and `docs/MIGRATION.md`); the Julia package is `Kirigami/`.

## What the method is

Segall et al. give a linear shape space `X` of uniformly deployable kirigami on any
oriented planar graph and leave open whether a member of `X` actually opens without
self-intersection (`IDEA.md` §1, `notes/paper_2026.md`). This project shows that every
geometric predicate along the uniform deployment branch is a single first harmonic in the
opening angle whose coefficients are quadratic in the flat design coordinates, so the
collision-free range `Θ_max` is exactly computable from a complete closed-form candidate
list with grazing contacts handled explicitly, and a sound sufficient certificate
`Θ_max ≥ ε` follows (`derivations/core.md` T3-T5, `REPORT.md` C1-C3). It then shows that
the paper's own Eq. (6) projection is certifiably non-deployable on random planar graphs
because it manufactures reflex face corners that jam the mechanism at `θ = 0⁺`, and that
replacing the projection's proximity objective by the `0⁺` deployment margin inside
convexity and split-outward barriers in the same null space makes most of the same
population deployable and certified (`REPORT.md` C7-C9, `results/kill/KILL_REPORT.md`,
`results/kill/k9c/`; measured yields: ⟨JULIA:k9c:307 of 400⟩).

## Repository layout

| path | contents |
|---|---|
| `Kirigami/` | the Julia package (`Project.toml`, one module `Kirigami`, nothing exported; see `Kirigami/README.md`) |
| `Kirigami/src/{core,method,export}/` | one `.jl` per C++ translation unit, same basenames (mesh, cut, holes, Tutte auxetic system, kinematics, collision, generators, RNG; deploy basis, contact calculus, design API, constrained/range-maximising embeddings, mobility, periodic Jacobian, expansive cone, budget; layout, material, SVG/STL/3MF writers) |
| `Kirigami/test/` | the unit tests, one file per C++ test file (`runtests.jl`, `helpers.jl`, `test_*.jl`, `derivation_tests.jl`) |
| `Kirigami/apps/` | user-facing CLIs `kiri_*.jl` and the experiment drivers `kill_*.jl`, plus the included helper files `common_app.jl`, `kill_common.jl`, `native_common.jl` |
| `Kirigami/gui/` | GLMakie desktop explorer (`app.jl`, window-free `model.jl`, own `Project.toml`, headless tests in `gui/test/`); replaces the C++/wasm web explorer |
| `Kirigami/scripts/` | CairoMakie plotting / summarising scripts `plot_*.jl` (own `Project.toml`), their Python originals, and the two shell drivers `run_e1.sh`, `run_scaling.sh` |
| `data/corpus/` | the frozen input populations of every experiment (`k1a_200.json`, `native200.json`, `k3a_500.json`, ..., `reference_cases_8.json`, `mt19937_vectors.json`), the C++ unit-test fixtures in `reference_patterns/`, and the method fixtures in `method_fixtures/`; see `data/corpus/README.md` |
| `data/deploy_frames.json`, `data/web_patterns.json` | the GUI's pre-characterised designs and the raw pattern library |
| `derivations/` | the theory: `core.md` (T1-T7), `lemmas.md`, `check.md` (the Checker's rounds), `scratch/check_*.{cpp,jl}` (the programs the derivations quote, C++ kept as evidence, Julia ports runnable) |
| `docs/paper/`, `docs/techreport/`, `docs/poster/` | LaTeX sources and PDFs of the paper, the technical report and the poster |
| `docs/MIGRATION.md`, `docs/CPP_MENTIONS.md`, `docs/REWORDING_PLAN.md` | what was migrated from the C++ repo, where the prose still says C++, and the planned rewording |
| `results/` | measured outputs: `core_validation/` (reference cases), `kill/` (K1a ... K9c, `KILL_REPORT.md`), `final/` (E1 and the paper figures), `regime/`, `scaling/`, `yield/`, `propagation/` |
| `export/` | fabrication samples: `hero/`, `hero2/` (JSON/SVG/STL/3MF/PNG of the two hero designs, with READMEs), `samples/`, `render_svg.py` |
| `baseline/` | the authors' 2025/2026 code (third-party, untouched) and the thin harness `baseline/native/` used as the parity referee (`baseline/parity.md`) |
| `papers/` | text transcriptions of the source papers and `papers/related/` |
| `IDEA.md`, `REPORT.md`, `STATE.md`, `ESCALATION.md` | the claim, the companion report with every number's provenance, the orchestrator state, the round-2 stop record |
| `notes/`, `ideas/`, `specs/`, `review/` | reading notes, ideation records, the agent specs of each phase, the theory/paper reviews (historical; they describe the C++ era) |
| `PORTING.md` | the porting conventions the Julia code follows |

## How to run

Everything runs from source; there is no build step. One-off setup:

```
export PATH=$HOME/.juliaup/bin:$PATH          # native arm64 Julia, see below
julia --project=Kirigami -e 'using Pkg; Pkg.instantiate()'
```

**Tests** (⟨JULIA:tests⟩):

```
julia --project=Kirigami -e 'using Pkg; Pkg.test()'
julia --project=Kirigami -e 'using Pkg; Pkg.test(test_args=["derivation_tests"])'   # one file (basename substring)
julia --project=Kirigami -e 'using Pkg; Pkg.test(test_args=["design", "export"])'   # several
```

`derivation_tests.jl` is the port of the Checker's `derivation_tests.cpp`
(⟨JULIA:derivation_tests⟩); the scratch programs behind `derivations/check.md` are run one
by one as `julia --project=Kirigami derivations/scratch/check_X.jl` (list in
`derivations/scratch/README.md`).

**A driver.** The user-facing CLIs and the experiment drivers take the same arguments as
their C++ namesakes and write the same file formats:

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
and JSON under `results/` and never re-run a solver; outputs currently carry a `_julia`
suffix next to the matplotlib originals for visual comparison:

```
julia --project=Kirigami/scripts -e 'using Pkg; Pkg.instantiate()'
julia --project=Kirigami/scripts Kirigami/scripts/plot_final.jl     # results/final/figures/fig_*_julia.png, summary_table_julia.md
julia --project=Kirigami/scripts Kirigami/scripts/plot_k9c.jl       # the K9c figures under results/kill/k9c/
```

To regenerate one figure end to end (data and picture), run the driver that produced its
CSV and then the script that reads it; `results/final/figures/README_julia.md` lists, per
figure, which files it reads. Example for the yield figure: `kill_k9c.jl` (400 designs;
hours), then `plot_final.jl`.

## Frozen corpus policy and `--regenerate`

Every input population the C++ built procedurally was written once to `data/corpus/`
(producer, toolchain and per-file provenance in `data/corpus/README.md`). The Julia apps
and tests read those files by default, so an experiment always runs on byte-identical
inputs whichever machine runs it. Each driver that takes a population also accepts
`--regenerate`, which rebuilds it through the ported generators and the bit-exact
`MT19937` (`Kirigami/apps/kill_common.jl`: `population(name; regenerate)`,
`reference_cases(; regenerate)`, `deployable_population(; regenerate)`). On the reference
platform `make_graph` regenerates all 200 graphs of `k1a_200.json` vertex for vertex; the
frozen file remains the definition of the population. The corpus is never rewritten by a
Julia run.

## Julia binary: use the native arm64 build

The C++ reference was built with clang on Apple arm64 and its numbers depend on that
platform's libm and on FMA contraction (`PORTING.md`, "Bit-exactness note"). Use the
native arm64 Julia from juliaup (`export PATH=$HOME/.juliaup/bin:$PATH`; the package
declares `julia = "1.10"`, the port was developed on 1.12). The x86_64 Homebrew build under
Rosetta is slower and its libm differs by 1 ulp in trig/hypot on a fraction of arguments,
which shows up as equal-objective tie flips in the orientation relaxation on symmetric
tilings (documented in `Kirigami/src/core/generators.jl`); random graphs are unaffected.
