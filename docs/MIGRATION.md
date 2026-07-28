# Migration of non-C++ material from `kirigami-experiments`

Date: 2026-09-19. Source: `~/Documents/kirigami-experiments` (read-only, untouched, at commit `c7c0fce`).
Destination: this repo. Everything below was copied with `rsync -a` / `cp`; nothing in the source was modified.

## What was copied

| source | destination | note |
|---|---|---|
| `IDEA.md`, `REPORT.md`, `STATE.md`, `ESCALATION.md` | same names, repo root | verbatim |
| `derivations/` | `derivations/` | includes `derivations/scratch/*.cpp` (the Checker's scratch programs referenced by page in `check.md`/`lemmas.md`; kept as documentary evidence, not built) |
| `ideas/`, `notes/`, `specs/`, `review/` | same | verbatim (`ideas/` includes `rigidity_rig_check.cpp`, same reasoning) |
| `papers/` | `papers/` | 4 `.txt` transcriptions + `papers/related/` (4 PDF + 4 txt) |
| `docs/paper/` | `docs/paper/` | `paper.tex`, `refs.bib`, `paper.pdf`, `SLOTS.md`, `build.sh` only |
| `docs/techreport/` | `docs/techreport/` | `techreport.tex`, `refs.bib`, `techreport.pdf` only |
| `export/` | `export/` | `render_svg.py`, `samples/`, `hero/`, `hero2/` (JSON/SVG/STL/3MF/PNG + READMEs) |
| `results/` | `results/` | see exclusions below |
| `code/scripts/*.py`, `*.sh` | `Kirigami/scripts/` | 13 Python plotting/summarising scripts stay Python (they read CSVs). The two shell drivers got a header `# TODO(julia-port): invoked C++ binary <name>; replace with Kirigami/apps/<name>.jl` (`run_e1.sh` -> `kill_e1`, `run_scaling.sh` -> `kill_scaling`). |
| `code/web/patterns.json` | `data/web_patterns.json` | the web explorer's pattern library |
| `baseline/` | `baseline/` | `README.md`, `parity.md`, `tuttekiri/` (upstream 2026 checkout), `kirigami_tessellations/` (upstream 2025 checkout), `native/` minus `build*/`. `baseline/native/` is **our** thin CMake harness (`CMakeLists.txt`, `src/{main,json_io,bind_port,their_system}.cpp`, `src_fixed/`, `emsc_stub/`, ~100 KB) around the **authors'** C++ — kept because it is the referee baseline, not part of the pipeline being ported. |

Files > 20 MB copied (2): `papers/related/jiang_choi_2026.pdf` (32 MB), `papers/related/isogami.pdf` (21 MB).
Everything else is under 20 MB (largest result file: `results/kill/k9c/` at 60 MB total across many small files).

## What was skipped and why

| skipped | size | why |
|---|---|---|
| `code/{src,apps,tests,CMakeLists.txt,README.md,build*,Kirigami}` | — | the C++ being ported (other agents' task); `code/Kirigami/` is an untracked local dir |
| `code/web/` except `patterns.json` | 3.4 MB | C++/wasm bindings + HTML/JS explorer, replaced by `Kirigami/gui` |
| `results/kill/cache_k2c/`, `results/kill/cache/`, `results/kill/b4/cache/`, `results/kill/k6/cache/`, `results/kill/b4/smoke/cache/` | 3.4 GB (1,870 files) | regenerable intermediates (`.gitignore` already excludes `results/**/cache/`) |
| `results/scaling/shape/` | 306 MB (39 `.bin`) | binary intermediates of the scaling run; `scaling.csv`/`capped.csv`/`SCALING.md` kept |
| other `results/**/*.bin` | 0 | none exist outside the two categories above |
| `baseline/native/{build,build_fixed,build2}/` | 362 MB | build trees incl. FetchContent clones of libigl/Eigen/Optiz/triangle |
| `docs/poster/` | 5.1 MB | not in the requested list (poster.tex/pdf, figs, `plot_poster_yield.py`) — flagging in case it is wanted |
| `docs/run_notes.{tex,pdf,aux,log,out,toc}` | 4.3 MB | not in the requested list — flagging |
| LaTeX by-products `*.aux *.log *.fls *.fdb_latexmk *.out *.blg *.bbl *.toc`, `docs/paper/build.log` | ~0.3 MB | regenerable |
| `g.json` | 9 KB | a graph-JSON (faces/vertices) test fixture, not a graphify graph — skipped as instructed; trivially recoverable |
| `code/scripts/__pycache__/`, `.DS_Store` | — | junk |
| `scripts/` (repo root) | 16K | not in the list; only `continue-claude-after-reset.sh` (session tooling) |

## Verification (`du -sh` / `find -type f | wc -l`, source vs destination)

| directory | source size | source files | dest size | dest files | reconciliation |
|---|---:|---:|---:|---:|---|
| `derivations/` | 564K | 15 | 564K | 15 | identical |
| `ideas/` | 608K | 11 | 608K | 11 | identical |
| `notes/` | 532K | 11 | 532K | 11 | identical |
| `specs/` | 164K | 26 | 164K | 26 | identical |
| `review/` | 52K | 3 | 52K | 3 | identical |
| `papers/` | 75M | 12 | 75M | 12 | identical |
| `export/` | 27M | 96 | 27M | 96 | identical |
| `results/` | 3.8G | 4,635 | 135M | 2,726 | 4,635 − 1,870 cache − 39 shape = 2,726 |
| `baseline/` | 441M | 7,788 | 79M | 185 | difference is `native/build*/` |
| `docs/paper/` | 9.8M | 13 | 9.7M | 5 | 8 LaTeX by-products + build.log dropped |
| `docs/techreport/` | 5.5M | 9 | 5.3M | 3 | 6 LaTeX by-products dropped |
| `code/scripts/` -> `Kirigami/scripts/` | 164K | 16 | 124K | 15 | `__pycache__` dropped |
| `code/web/patterns.json` -> `data/web_patterns.json` | 20K | 1 | 20K | 1 | identical |

Total migrated: **333 MB, 3,151 files** (excluding `.git`).

## C++ mentions

`docs/CPP_MENTIONS.md` lists **942 lines in 61 files** of the migrated `.md/.tex/.bib` that mention the C++ toolchain, C++ file paths, test binaries/counts, or the wasm web explorer. None have been edited. The heaviest are `notes/repo_2025.md` (116, third-party — keep), `results/kill/KILL_REPORT.md` (114), `docs/techreport/techreport.tex` (74), `derivations/core.md` (51), `derivations/check.md` (46), `derivations/lemmas.md` (42), `ideas/round2_adversary.md` (40), `REPORT.md` (36), `ideas/round2_theorist_b.md` (36), `STATE.md` (31).

## Julia port status (2026-09-20)

What of `kirigami-experiments/code/` now exists in Julia, file by file. "Ported" means the
Julia file mirrors the C++ translation unit under the conventions of `PORTING.md` and its
tests pass on the frozen corpora; the numerical caveats are in `Kirigami/README.md`
("Reproducibility of the numbers") and `data/corpus/method_fixtures/README_design.md`.

| C++ | Julia | status |
|---|---|---|
| `src/core/*.{hpp,cpp}` (13 units) | `Kirigami/src/core/*.jl`, same basenames, plus `core/kill_common.jl` (the population definitions of `apps/kill_common.hpp`) | ported |
| `src/method/*.{hpp,cpp}` (11 units) | `Kirigami/src/method/*.jl`, same basenames | ported |
| `src/export/*.{hpp,cpp}` (8 units) | `Kirigami/src/export/*.jl`, same basenames | ported |
| `tests/test_*.cpp`, `derivation_tests.cpp`, `test_main.cpp`, `helpers.hpp` | `Kirigami/test/test_*.jl`, `derivation_tests.jl`, `runtests.jl`, `helpers.jl`; `test_method.cpp` is split into `test_method_{1,2,3}.jl` | ported (`⟨JULIA:tests⟩`, `⟨JULIA:derivation_tests⟩`) |
| `apps/common_app.hpp`, `kill_common.hpp`, `native_common.hpp` | `Kirigami/apps/common_app.jl`, `kill_common.jl`, `native_common.jl` (included files; populations read from `data/corpus/`, `--regenerate` rebuilds them) | ported |
| `apps/kiri_{gen,analyze,deploy,design,export,reference,sweep}.cpp` | `Kirigami/apps/kiri_*.jl`, same CLI arguments and output formats (checked against the C++ binaries on the same inputs) | ported |
| `apps/kill_{e1,f23,jitter,k1a,k1b,k1c,k2a,k2b,k2c,k3a,k3a_recheck,k5,k6,k7,k8a,k8a_recheck,k9,k9b,k9c,t1}.cpp` | `Kirigami/apps/kill_*.jl`, same names | ported |
| `apps/kill_{b3,b4,basin,basin_agg,native200,regime,scaling,yield}.cpp`, `native200_cell.cpp`, `native200_merge.cpp` | not yet in `Kirigami/apps/` at the time of writing | pending |
| `apps/freeze_corpus.cpp` | -- | not needed: it produced `data/corpus/` once; the Julia side only reads it |
| `derivations/scratch/check_*.cpp`, `ideas/rigidity_rig_check.cpp` | `check_*.jl`, `rigidity_rig_check.jl` next to them | ported (C++ kept as documentary evidence) |
| `scripts/*.py`, `*.sh` | `Kirigami/scripts/plot_*.jl`, `summarise_k8a.jl` (CairoMakie, own `Project.toml`); the Python originals kept beside them; `run_e1.sh` / `run_scaling.sh` still name the C++ binaries (`TODO(julia-port)` header) | ported (figures regenerated with a `_julia` suffix) |
| `web/` (wasm bindings, HTML/JS explorer) | `Kirigami/gui/` (GLMakie desktop app, `app.jl` + window-free `model.jl`, own `Project.toml` and tests); `web/patterns.json` -> `data/web_patterns.json` | replaced, deliberately not ported |
| `apps/dbg_*.cpp` (`dbg_cert`, `dbg_k1b`, `dbg_k2a`, `dbg_k2a_out`, `dbg_k2c`, `dbg_k5`, `dbg_k6`, `dbg_k9`, `dbg_t422`) | -- | deliberately not ported: one-off debugging probes of the C++ era whose findings are recorded in `STATE.md` / `results/kill/KILL_REPORT.md` |
| `CMakeLists.txt`, `README.md` (C++ build and graph contract) | `Kirigami/Project.toml`, `Kirigami/README.md`, this repo's `README.md` | replaced |

Deliberately not ported, in one line: the `dbg_*` apps and the web explorer (the latter
replaced by `Kirigami/gui/`). Everything else in `code/` has a Julia counterpart, except
the ten apps marked pending above (eight `kill_*` drivers and the two `native200_*` helpers).
