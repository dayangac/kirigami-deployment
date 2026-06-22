# Rewording plan: C++ → Julia in the migrated prose

Prepared 2026-09-19 from `docs/CPP_MENTIONS.md`. **Nothing here is applied.** Apply only after the Julia port reproduces the C++ numbers; every placeholder `⟨JULIA:…⟩` is filled from the Julia run named in the table at the end.

## Scope and conventions

**In scope** (every non-KEEP, non-HISTORICAL file of `CPP_MENTIONS.md`): `docs/techreport/techreport.tex`, `docs/paper/paper.tex`, `docs/poster/poster.tex` (scanned in the source repo: **0 mentions**, nothing to do), `REPORT.md`, `STATE.md`, `IDEA.md`, `derivations/*.md`, `results/**/*.md`, `export/hero*/README.md`, `review/*.md`.
**Out of scope, left verbatim:** `baseline/**` and `notes/repo_2025.md`, `notes/paper_2026.md` (third-party C++ — KEEP); `specs/*`, `ideas/*`, `notes/screen_r2.md` (HISTORICAL agent specs / ideation records; they describe what was true in the C++ era and are not documentation of the current implementation). `STATE.md`: only its five current-state lines (7, 13, 67, 173, 184) are reworded; the dated log entries stay as history and are listed as *leave*.

**Mechanical rules applied** (R1–R10), so the same C++ token always maps to the same Julia wording:

| rule | old | new |
|---|---|---|
| R1 path | `code/src/{core,method,export}/X.{hpp,cpp}`, `X.hpp`, `X.cpp` (X a module) | `Kirigami/src/{core,method,export}/X.jl`, `X.jl` |
| R2 tests | `code/tests/test_X.cpp`, `tests/test_X.cpp` | `Kirigami/test/test_X.jl`, `test/test_X.jl` |
| R3 apps | `code/apps/kill_X.cpp`, `code/apps/dbg_X.cpp`, `apps/kill_X.cpp` | `Kirigami/apps/kill_X.jl`, … |
| R4 corpus | `kill_common.hpp`, `kill_common::make_graph` | `kill_common.jl` (populations frozen in `data/corpus/`) — see assumption A1 |
| R5 binaries | `./code/build/kill_X …` | `julia --project=Kirigami Kirigami/apps/kill_X.jl …` |
| R6 suites | `./code/build/kiri_tests`, `` `kiri_tests` `` | `julia --project=Kirigami -e 'using Pkg; Pkg.test()'`, `` `Pkg.test()` `` |
| R6b | `./code/build/derivation_tests`, `` `derivation_tests` `` | `julia --project=Kirigami -e 'using Pkg; Pkg.test(test_args=["derivation_tests"])'`, `` `derivation_tests.jl` `` — see A5 |
| R7 build | `cmake -S code -B code/build … && cmake --build code/build -j` | `julia --project=Kirigami -e 'using Pkg; Pkg.instantiate()'   # one-off; no build step` |
| R7b | `cmake --build … --target kill_X` | `# no build step: Kirigami/apps/kill_X.jl runs from source` |
| R8 scratch | `derivations/scratch/check_X.cpp` + its `clang++ …` build line | `derivations/scratch/check_X.jl`, run with `julia --project=Kirigami derivations/scratch/check_X.jl` — see A2 |
| R9 terms | `doctest`/`doctests`, `TEST_CASE`, `Eigen::ColPivHouseholderQR`, `nlohmann/json`, `code/CMakeLists.txt`, `code/README.md` | `unit test(s)`, `@testset`, `` `qr(A, ColumnNorm())` (LinearAlgebra) ``, `JSON.jl`, `Kirigami/Project.toml`, `Kirigami/README.md` |
| R10 counts | any `N cases / M assertions`, doctest output blocks | `⟨JULIA:tests⟩` or `⟨JULIA:derivation_tests⟩` |
| R11 values | measured values with a decimal mantissa (`1.42e−14`) on a line that cites the program that measured them; `0 mismatches`, `0 violations`, `307 of 400`, `7 of 7` | `⟨JULIA:<program>:<old value>⟩` (the old value is kept inside the placeholder so the applier knows what it replaced) |
| — | bare driver names in prose (`` `kill_k5` ``, `kill_k8a: ten minutes`) | **unchanged** — they name the experiment/driver, and the Julia driver keeps the name (`kill_k5.jl`) |
| — | tolerances (`1e-9`, `1e−12`), line numbers, equation numbers | **unchanged** |

**Assumptions the plan depends on** (each is a decision for the port, not a fact yet):

- **A1** The experiment drivers share a `Kirigami/apps/kill_common.jl` (same name as the C++ header) whose populations are loaded from the frozen `data/corpus/*.json`, never regenerated. If the helper ends up with another name, replace `kill_common.jl` throughout.
- **A2** The Checker's scratch programs `derivations/scratch/check_{t1_t2,t4_t5,r2,r3,l1,l2,b1,t5_adm,b4_93}.cpp` get Julia ports with the same basenames (`.jl`), runnable as `julia --project=Kirigami derivations/scratch/check_X.jl [ε]`. Until they exist, the derivations' "[N] Measured" paragraphs cite programs that cannot be re-run; the placeholders make that visible.
- **A3** The web explorer is replaced by a GLMakie desktop app at `Kirigami/gui/app.jl` (per the team lead; `Kirigami/gui/` is empty today). GLMakie is **not** yet in `Kirigami/Project.toml` — it should live in a separate `Kirigami/gui/Project.toml` so the core package stays light; the wording below assumes that.
- **A4** The throwaway graph-dump utilities in `export/hero*/README.md` (`dump_hero_graph.cpp` etc.) are re-created as `export/hero*/dump_*.jl` scripts (they are 20-line wrappers around `make_graph`/`characterize`; with the corpus frozen they may reduce to reading `data/corpus/k1a_200.json`).
- **A5** `Kirigami/test/runtests.jl` accepts `test_args` to run a single file (`Pkg.test(test_args=["derivation_tests"])`). Today it runs every file in `TEST_FILES` unconditionally; a 3-line filter on `ARGS` is needed, otherwise replace with the whole-suite command everywhere.
- **A6** Determinism wording: the C++ text relies on "Eigen built without OpenMP". The Julia equivalent proposed is `BLAS.set_num_threads(1)` at app start plus the single-module, no-global-state design; `test_design.jl` must check bit-identical repeated calls as `test_design.cpp` did.
- **A7** Julia version: `Project.toml` says `julia = "1.10"` (compat), the lead's target is 1.12. The wording uses "Julia ≥ 1.10"; change to "Julia 1.12" if the compat bound is raised.

## Draft "Implementation" paragraph (paper and tech report)

**Paper (`docs/paper/paper.tex`, replaces the toolchain sentence wherever the C++ version stood; today the paper has only the one `\plan{}` path reference, so this is new text for the Methods/Reproducibility paragraph):**

> All algorithms are implemented in a single Julia package, \texttt{Kirigami} (Julia $\geq 1.10$; \texttt{LinearAlgebra}, \texttt{SparseArrays}, \texttt{StaticArrays}, \texttt{JSON}, \texttt{DelaunayTriangulation}; no external optimiser). The package is organised as \texttt{src/core} (planar mesh, cuts, hole detection, the Tutte-auxetic linear system, forward kinematics, collision, orientation, generators), \texttt{src/method} (deploy basis, exact contact angles, the $0^+$ calculus, convexity- and range-constrained embedding, design, mobility, expansive cone, periodic Jacobian) and \texttt{src/export} (SVG/STL/3MF fabrication files). Every experiment in this paper is a script under \texttt{apps/} that reads a frozen graph population from \texttt{data/corpus/} and writes CSV; the plots are Python/matplotlib over those CSVs. The test suite (\texttt{Pkg.test()}, ⟨JULIA:tests⟩) includes \texttt{derivation\_tests.jl}, which re-checks every numbered identity of the derivations on the same corpus (⟨JULIA:derivation\_tests⟩). A GLMakie desktop explorer (\texttt{gui/app.jl}) animates the exact kinematic path $Y(\theta)=\cos(\theta/2)C+\sin(\theta/2)S$ from the cached $(C,S)$ basis.

**Tech report (`docs/techreport/techreport.tex` §Software, replaces l.3463–3466; the same text minus the last sentence is the replacement for l.429–431):**

> Everything is Julia ($\geq 1.10$), one package \file{Kirigami} (\file{Project.toml}). Dependencies are the standard library (\file{LinearAlgebra} for the dense factorisations --- \file{qr(A, ColumnNorm())}, \file{svd}, \file{nullspace} --- \file{SparseArrays}/SuiteSparse for the sparse ones, \file{Random}, \file{Statistics}, \file{Printf}) plus \file{StaticArrays} (\file{SVector\{2,Float64\}} for every planar point), \file{JSON} (the interchange contract) and \file{DelaunayTriangulation} (Delaunay/Voronoi generators). The L-BFGS with backtracking line search in \file{core/optimize.jl} is ours, so there is no external optimiser anywhere in the pipeline. The C++ \file{std::mt19937} and libc++ \file{uniform\_real\_distribution} are ported bit-exactly (\file{core/mt19937.jl}) so every seeded population is reproduced; the populations the experiments read are nevertheless frozen once in \file{data/corpus/} and never regenerated. Layout: \file{src/core}, \file{src/method}, \file{src/export} (one file per former C++ translation unit, same basenames, same function names); \file{test/} (one file per former test unit, \file{Pkg.test()}); \file{apps/} (the \file{kiri\_*} CLIs and the \file{kill\_*} experiment drivers, run from source with \file{julia --project=Kirigami apps/X.jl}); \file{gui/app.jl} (GLMakie desktop explorer, its own \file{Project.toml}); \file{scripts/} (matplotlib plotting of the CSVs). Python appears only for plotting.

## Proposed edits, by file

Format: `L<line>` (or `L<a>–<b>` for a block) then the old text and the proposed new text. *leave* entries are matched lines that need no change, with the reason.

### `docs/techreport/techreport.tex` — 56 edits, 1 leave

- **L82–84**
  ```
  OLD: soundness violations of the certificate. All code is C++20 with CMake; the two test
       suites re-run for this report at $109$ cases / $17{,}854$ assertions and $33$ cases /
       $89{,}074$ assertions, no failures.
  NEW: soundness violations of the certificate. All code is Julia (package \file{Kirigami}); the two test
       suites re-run for this report at ⟨JULIA:tests⟩ (whole suite) and ⟨JULIA:derivation_tests⟩ (\file{derivation\_tests.jl}), no failures.
  ```
- **L429–431**
  ```
  OLD: The reimplementation is C++20 + CMake, header-only dependencies (Eigen, nlohmann/json,
       doctest), no external optimiser --- the L-BFGS with backtracking line search in
       \file{core/optimize.hpp} is ours. The library layout is in Section~\ref{sec:software}.
  NEW: The reimplementation is a Julia package (\file{Kirigami}, Julia $\geq 1.10$) on the standard library
       (\file{LinearAlgebra}, \file{SparseArrays}, \file{Random}, \file{Statistics}) plus \file{StaticArrays}, \file{JSON}
       and \file{DelaunayTriangulation}; no external optimiser --- the L-BFGS with backtracking line search in
       \file{core/optimize.jl} is ours. The library layout is in Section~\ref{sec:software}.
  ```
- **L434** *leave* — third-party `baseline/` reference, KEEP: `(\file{baseline/native}, libigl 2.5, Optiz, Eigen 3.4 via FetchContent, their \file{bind.cpp}`
- **L535**
  ```
  OLD: \textbf{[F, \file{core/kinematics.cpp}] The code's sign convention.} Face $f$ is transformed
  NEW: \textbf{[F, \file{core/kinematics.jl}] The code's sign convention.} Face $f$ is transformed
  ```
- **L634**
  ```
  OLD: \textbf{[N] Path-independence.} \file{derivations/scratch/check\_t1\_t2.cpp}, check C1: build
  NEW: \textbf{[N] Path-independence.} \file{derivations/scratch/check\_t1\_t2.jl}, check C1: build
  ```
- **L668**
  ```
  OLD: \caption{Assembling the deploy basis $(C,S)$ --- \file{method/deploy\_basis.cpp}}
  NEW: \caption{Assembling the deploy basis $(C,S)$ --- \file{method/deploy\_basis.jl}}
  ```
- **L689**
  ```
  OLD: \textbf{A naming clash, recorded.} The header of \file{deploy\_basis.hpp} writes
  NEW: \textbf{A naming clash, recorded.} The header of \file{deploy\_basis.jl} writes
  ```
- **L697**
  ```
  OLD: independent BFS, against \file{kinematics.cpp::deploy()} at
  NEW: independent BFS, against \file{kinematics.jl::deploy()} at
  ```
- **L1181**
  ```
  OLD: \caption{Exact $\Tmax$ (T4.2$''$) --- \file{method/contact.cpp::exact\_theta\_max\_overlap}}
  NEW: \caption{Exact $\Tmax$ (T4.2$''$) --- \file{method/contact.jl::exact\_theta\_max\_overlap}}
  ```
- **L1288**
  ```
  OLD: \caption{Broad phase --- \file{method/contact.cpp::candidate\_pairs}}
  NEW: \caption{Broad phase --- \file{method/contact.jl::candidate\_pairs}}
  ```
- **L1584**
  ```
  OLD: This is exactly what \file{contact.cpp} implements --- before round 6 the derivation demanded
  NEW: This is exactly what \file{contact.jl} implements --- before round 6 the derivation demanded
  ```
- **L1669**
  ```
  OLD: \textbf{[N] The measured counterexample} (\file{derivations/scratch/check\_t5\_adm.cpp},
  NEW: \textbf{[N] The measured counterexample} (\file{derivations/scratch/check\_t5\_adm.jl},
  ```
- **L1725**
  ```
  OLD: \file{contact.cpp}; a doctest fails without it and passes with it.
  NEW: \file{contact.jl}; a unit test fails without it and passes with it.
  ```
- **L1746**
  ```
  OLD: bisection referee (\file{results/final/e1/E1.md}, driver \file{code/apps/kill\_e1.cpp}, which
  NEW: bisection referee (\file{results/final/e1/E1.md}, driver \file{Kirigami/apps/kill\_e1.jl}, which
  ```
- **L2149**
  ```
  OLD: \file{code/src/method/contact.cpp:339} --- keep a root only if it is smaller than the
  NEW: \file{Kirigami/src/method/contact.jl:339} --- keep a root only if it is smaller than the
  ```
- **L2303**
  ```
  OLD: \file{zero\_plus\_q} and \file{zero\_plus\_corner\_margin} at the $\sigma$ ray, in a doctest,
  NEW: \file{zero\_plus\_q} and \file{zero\_plus\_corner\_margin} at the $\sigma$ ray, in a unit test,
  ```
- **L2305**
  ```
  OLD: one free translation pair per component of $\Gamma$; a doctest checks
  NEW: one free translation pair per component of $\Gamma$; a unit test checks
  ```
- **L2447**
  ```
  OLD: A doctest pins this. \emph{Any claim of the form ``authored tilings pass 8/8'' is really
  NEW: A unit test pins this. \emph{Any claim of the form ``authored tilings pass 8/8'' is really
  ```
- **L2717**
  ```
  OLD: surrogate. Two doctests check exactly those two claims
  NEW: surrogate. Two unit tests check exactly those two claims
  ```
- **L2718**
  ```
  OLD: (\file{code/tests/test\_range\_embed.cpp}), and two more finite-difference the analytic
  NEW: (\file{Kirigami/test/test\_range\_embed.jl}), and two more finite-difference the analytic
  ```
- **L2737**
  ```
  OLD: objective of \file{range\_opt.hpp}, whose gradients are the implicit-function-theorem
  NEW: objective of \file{range\_opt.jl}, whose gradients are the implicit-function-theorem
  ```
- **L2774**
  ```
  OLD: \caption{\texttt{design\_range\_max} --- the deliverable (K9c). \file{method/range\_embed.cpp}, \file{method/design.cpp}}
  NEW: \caption{\texttt{design\_range\_max} --- the deliverable (K9c). \file{method/range\_embed.jl}, \file{method/design.jl}}
  ```
- **L2835–2837**
  ```
  OLD: iterates an unordered container or touches global state, and Eigen is built without OpenMP so
       every reduction runs in one fixed order. Repeated calls are bit-identical and
       \file{tests/test\_design.cpp} checks that. Two consequences worth stating:
  NEW: iterates an unordered container or touches global state, and BLAS is pinned to one thread
       (\file{BLAS.set\_num\_threads(1)}) so every reduction runs in one fixed order. Repeated calls are
       bit-identical and \file{test/test\_design.jl} checks that. Two consequences worth stating:
  ```
- **L2843**
  ```
  OLD: Eigen version or assembly order only if $\Phi$ is \textbf{carried along, not recomputed}.
  NEW: Julia/BLAS version or assembly order only if $\Phi$ is \textbf{carried along, not recomputed}.
  ```
- **L2845**
  ```
  OLD: \item \file{apps/kill\_k9c.cpp} is the one place that was genuinely non-deterministic: it
  NEW: \item \file{apps/kill\_k9c.jl} is the one place that was genuinely non-deterministic: it
  ```
- **L2855**
  ```
  OLD: \file{tests/test\_design.cpp} locks four rows to $10^{-9}$: \file{delaunay\_130}
  NEW: \file{test/test\_design.jl} locks four rows to $10^{-9}$: \file{delaunay\_130}
  ```
- **L3047**
  ```
  OLD: \file{method/periodic\_jacobian.cpp} builds the genuine quotient: one face per translation
  NEW: \file{method/periodic\_jacobian.jl} builds the genuine quotient: one face per translation
  ```
- **L3463–3466**
  ```
  OLD: Everything is C++17/20 with CMake. Python appears only for matplotlib plotting of CSV/JSON
       dumped by the C++ programs. Dependencies are header-only (Eigen, nlohmann/json, doctest); the
       L-BFGS with backtracking line search is ours, so there is no external optimiser anywhere in
       the pipeline.
  NEW: Everything is Julia ($\geq 1.10$), one package \file{Kirigami}. Python appears only for matplotlib plotting of CSV/JSON
       dumped by the Julia apps. Dependencies are the standard library (\file{LinearAlgebra}, \file{SparseArrays},
       \file{Random}, \file{Statistics}, \file{Printf}) plus \file{StaticArrays} ($\mathbb{R}^2$ vectors), \file{JSON} (I/O) and
       \file{DelaunayTriangulation} (random-point generators); the L-BFGS with backtracking line search is ours, so there is
       no external optimiser anywhere in the pipeline. Dense factorisations are \file{LinearAlgebra} (\file{qr(A, ColumnNorm())},
       \file{svd}), sparse ones \file{SparseArrays}/SuiteSparse; \file{std::mt19937} is ported bit-exactly (\file{core/mt19937.jl}).
  ```
- **L3471–3473**
  ```
  OLD: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
       ./code/build/kiri_tests          # 109 cases, 17,854 assertions, 0 failures
       ./code/build/derivation_tests    # 33 cases, 89,074 assertions, 0 failures
  NEW: julia --project=Kirigami -e 'using Pkg; Pkg.instantiate()'   # one-off
       julia --project=Kirigami -e 'using Pkg; Pkg.test()'          # ⟨JULIA:tests⟩, 0 failures; includes derivation_tests.jl (⟨JULIA:derivation_tests⟩)
  ```
- **L3489–3495**
  ```
  OLD: core pipeline (2026 reimplementation) & \file{code/src/core/}: \file{mesh} \file{cut} \file{holes} \file{orientation} \file{tutte\_auxetic} \file{kinematics} \file{collision} \file{optimize} \file{rank\_checks} \file{generators} \\
       method layer & \file{code/src/method/}: \file{deploy\_basis} \file{contact} \file{design} \file{range\_opt} \file{zero\_plus} \file{convex\_embed} \file{range\_embed} \file{mobility} \file{expansive\_cone} \file{periodic\_jacobian} \file{budget} \\
       fabrication export & \file{code/src/export/}: \file{layout} \file{solid} \file{svg} \file{stl} \file{threemf} \file{xml} \file{zip} \file{material} \\
       tests & \file{code/tests/}: \file{test\_mesh\_cut} \file{test\_holes} \file{test\_system} \file{test\_kinematics} \file{test\_collision} \file{test\_method} \file{test\_design} \file{test\_range\_embed} \file{test\_rank\_checks} \file{test\_reference\_cases} \file{test\_export} \file{derivation\_tests} \\
       CLIs & \file{code/apps/}: \file{kiri\_gen} \file{kiri\_analyze} \file{kiri\_deploy} \file{kiri\_design} \file{kiri\_export} \file{kiri\_sweep} \file{kiri\_reference} \\
       experiment drivers & \file{code/apps/kill\_*}: \file{k1a} \file{k1b} \file{k1c} \file{k2a} \file{k2b} \file{k2c} \file{k3a} \file{k5} \file{k6} \file{k7} \file{k8a} \file{k9} \file{k9b} \file{k9c} \file{b3} \file{b4} \file{t1} \file{jitter} \file{e1} \file{native200} \file{f23} \\
       web explorer & \file{code/web/}: \file{bindings.cpp} \file{build.sh} \file{index.html} \file{app.html} \file{kiri.js} \\
  NEW: core pipeline (2026 reimplementation) & \file{Kirigami/src/core/}: \file{mesh} \file{cut} \file{holes} \file{orientation} \file{tutte\_auxetic} \file{kinematics} \file{collision} \file{optimize} \file{rank\_checks} \file{generators} \file{import\_soup} \file{mt19937} \\
       method layer & \file{Kirigami/src/method/}: \file{deploy\_basis} \file{contact} \file{design} \file{range\_opt} \file{zero\_plus} \file{convex\_embed} \file{range\_embed} \file{mobility} \file{expansive\_cone} \file{periodic\_jacobian} \file{budget} \\
       fabrication export & \file{Kirigami/src/export/}: \file{layout} \file{solid} \file{svg} \file{stl} \file{threemf} \file{xml} \file{zip} \file{material} \\
       tests & \file{Kirigami/test/}: \file{test\_mesh\_cut} \file{test\_holes} \file{test\_system} \file{test\_kinematics} \file{test\_collision} \file{test\_method} \file{test\_design} \file{test\_range\_embed} \file{test\_rank\_checks} \file{test\_reference\_cases} \file{test\_export} \file{test\_import\_soup} \file{test\_mt19937} \file{derivation\_tests} \\
       CLIs & \file{Kirigami/apps/}: \file{kiri\_gen} \file{kiri\_analyze} \file{kiri\_deploy} \file{kiri\_design} \file{kiri\_export} \file{kiri\_sweep} \file{kiri\_reference} \\
       experiment drivers & \file{Kirigami/apps/kill\_*}: \file{k1a} \file{k1b} \file{k1c} \file{k2a} \file{k2b} \file{k2c} \file{k3a} \file{k5} \file{k6} \file{k7} \file{k8a} \file{k9} \file{k9b} \file{k9c} \file{b3} \file{b4} \file{t1} \file{jitter} \file{e1} \file{native200} \file{f23} \\
       frozen corpora & \file{data/corpus/}: the graph populations every driver reads (produced once, never regenerated) \\
       desktop explorer & \file{Kirigami/gui/}: \file{app.jl} (GLMakie) \\
  ```
- **L3506**
  ```
  OLD: the web explorer and compared byte for byte:
  NEW: the desktop explorer and compared byte for byte:
  ```
- **L3524**
  ```
  OLD: \file{method/design.hpp} is the single entry point, so that ``certified'' means one thing
  NEW: \file{method/design.jl} is the single entry point, so that ``certified'' means one thing
  ```
- **L3632–3650**
  ```
  OLD: \subsection{The web explorer}
       
       \file{code/web/} compiles \file{src/core} $+$ \file{src/method} $+$ \file{src/export} to
       WebAssembly via emscripten (\file{SINGLE\_FILE=1} embeds the wasm as base64 inside the
       \file{.js}, \file{MODULARIZE=1} gives one factory the page and node can await; no threads, no
       filesystem; the native CMake build is untouched). Everything crossing the JS boundary is a
       JSON string in \emph{exactly} the contract above, so the same payloads can be fed to
       \file{kiri\_gen} / \file{kiri\_design} / \file{kiri\_analyze} from a terminal and compared byte
       for byte. No numerics live in the bindings: every function is a thin marshalling wrapper.
       
       The surface is \file{generate}, \file{orient}, \file{design}, \file{characterize},
       \file{parse\_pattern}, \file{null\_space} (which caches $X_0$ and $\Phi$), \file{characterize\_at}
       (characterise at $X_0+\Phi t$ from the cache), \file{export\_svg} and \file{set\_progress}.
       The point of the design is that \textbf{the page animates the true kinematic path}: it holds
       the $(C,S)$ basis and evaluates $Y(\theta)=\cos(\theta/2)C+\sin(\theta/2)S$ in JavaScript, so
       dragging the opening angle re-solves nothing and is exact rather than interpolated.
       
       \clearpage
       %==============================================================================
  NEW: \subsection{The desktop explorer}
       
       \file{Kirigami/gui/app.jl} is a GLMakie desktop application that calls \file{src/core} $+$
       \file{src/method} $+$ \file{src/export} in-process (no bindings, no marshalling layer; the same
       Julia functions the CLIs and tests call). Every payload it loads or saves is the JSON contract
       above, so the same files can be fed to \file{kiri\_gen} / \file{kiri\_design} /
       \file{kiri\_analyze} from a terminal and compared byte for byte. No numerics live in the GUI:
       it only holds state and draws.
       
       The surface is the same as the CLIs': generate, orient, design, characterise,
       parse pattern, null space (which caches $X_0$ and $\Phi$), characterise at $X_0+\Phi t$ from the
       cache, export SVG. The point of the design is unchanged: \textbf{the window animates the true
       kinematic path}. It holds the $(C,S)$ basis and evaluates $Y(\theta)=\cos(\theta/2)C+\sin(\theta/2)S$
       as the opening-angle slider moves, so dragging the angle re-solves nothing and is exact rather
       than interpolated.
  ```
- **L3820**
  ```
  OLD: \file{expansive\_cone.hpp}; the C-IRIS line supplies the general pattern of
  NEW: \file{expansive\_cone.jl}; the C-IRIS line supplies the general pattern of
  ```
- **L3912–3914**
  ```
  OLD: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
       ./code/build/kiri_tests           # 109 cases, 17854 assertions, 0 failures
       ./code/build/derivation_tests     # 33 cases, 89074 assertions, 0 failures
  NEW: julia --project=Kirigami -e 'using Pkg; Pkg.instantiate()'
       julia --project=Kirigami -e 'using Pkg; Pkg.test()'   # ⟨JULIA:tests⟩ incl. derivation_tests.jl ⟨JULIA:derivation_tests⟩, 0 failures
  ```
- **L3919**
  ```
  OLD: cmake --build code/build --target kill_e1
  NEW: # no build step; the driver is Kirigami/apps/kill_e1.jl (run_e1.sh invokes it)
  ```
- **L3926**
  ```
  OLD: ./code/build/kill_k1b  --out results/kill/k1b          # ~34 s
  NEW: julia --project=Kirigami Kirigami/apps/kill_k1b.jl  --out results/kill/k1b          # ~34 s
  ```
- **L3927**
  ```
  OLD: ./code/build/kill_k1a  --out results/kill/k1a          # ~40 s
  NEW: julia --project=Kirigami Kirigami/apps/kill_k1a.jl  --out results/kill/k1a          # ~40 s
  ```
- **L3928**
  ```
  OLD: ./code/build/kill_k2a  --out results/kill/k2a          # ~7 s
  NEW: julia --project=Kirigami Kirigami/apps/kill_k2a.jl  --out results/kill/k2a          # ~7 s
  ```
- **L3929**
  ```
  OLD: ./code/build/kill_k2c  --out results/kill/k2c          # ~3 min
  NEW: julia --project=Kirigami Kirigami/apps/kill_k2c.jl  --out results/kill/k2c          # ~3 min
  ```
- **L3930**
  ```
  OLD: ./code/build/kill_k3a  --out results/kill/k3a          # ~25 min
  NEW: julia --project=Kirigami Kirigami/apps/kill_k3a.jl  --out results/kill/k3a          # ~25 min
  ```
- **L3931**
  ```
  OLD: ./code/build/kill_k3a_recheck --out results/kill/k3a   # ~4 min
  NEW: julia --project=Kirigami Kirigami/apps/kill_k3a_recheck.jl --out results/kill/k3a   # ~4 min
  ```
- **L3932**
  ```
  OLD: ./code/build/kill_k5   --out results/kill/k5 --n 200 --maxf 800    # ~5 min
  NEW: julia --project=Kirigami Kirigami/apps/kill_k5.jl   --out results/kill/k5 --n 200 --maxf 800    # ~5 min
  ```
- **L3933**
  ```
  OLD: ./code/build/kill_k1c  --out results/kill/k1c          # ~31 s   (needs baseline/native)
  NEW: julia --project=Kirigami Kirigami/apps/kill_k1c.jl  --out results/kill/k1c          # ~31 s   (needs baseline/native)
  ```
- **L3934**
  ```
  OLD: ./code/build/kill_k2b  --out results/kill/k2b --maxf 160   # ~26 min (needs native)
  NEW: julia --project=Kirigami Kirigami/apps/kill_k2b.jl  --out results/kill/k2b --maxf 160   # ~26 min (needs native)
  ```
- **L3935**
  ```
  OLD: ./code/build/kill_f23  --out results/kill/f23          # ~1 s
  NEW: julia --project=Kirigami Kirigami/apps/kill_f23.jl  --out results/kill/f23          # ~1 s
  ```
- **L3936**
  ```
  OLD: ./code/build/kill_k6   --out results/kill/k6
  NEW: julia --project=Kirigami Kirigami/apps/kill_k6.jl   --out results/kill/k6
  ```
- **L3937**
  ```
  OLD: ./code/build/kill_t1   --out results/kill/t1           # ~7 s off the K6 cache
  NEW: julia --project=Kirigami Kirigami/apps/kill_t1.jl   --out results/kill/t1           # ~7 s off the K6 cache
  ```
- **L3938**
  ```
  OLD: ./code/build/kill_jitter --shard I --nshards 12 --out results/kill/jitter  # ~2 s
  NEW: julia --project=Kirigami Kirigami/apps/kill_jitter.jl --shard I --nshards 12 --out results/kill/jitter  # ~2 s
  ```
- **L3939**
  ```
  OLD: ./code/build/kill_k7   --stage main --shard I --nshard 12 --out results/kill/k7
  NEW: julia --project=Kirigami Kirigami/apps/kill_k7.jl   --stage main --shard I --nshard 12 --out results/kill/k7
  ```
- **L3940**
  ```
  OLD: ./code/build/kill_b3   --out results/kill/b3
  NEW: julia --project=Kirigami Kirigami/apps/kill_b3.jl   --out results/kill/b3
  ```
- **L3941**
  ```
  OLD: ./code/build/kill_b4   --n 200 --nshards 12 --shard I \
  NEW: julia --project=Kirigami Kirigami/apps/kill_b4.jl   --n 200 --nshards 12 --shard I \
  ```
- **L3943**
  ```
  OLD: ./code/build/kill_k8a  --n 100 --x0 --dual-iters 100000 --shard I --nshard 8 \
  NEW: julia --project=Kirigami Kirigami/apps/kill_k8a.jl  --n 100 --x0 --dual-iters 100000 --shard I --nshard 8 \
  ```
- **L3945**
  ```
  OLD: ./code/build/kill_k9   --out results/kill/k9
  NEW: julia --project=Kirigami Kirigami/apps/kill_k9.jl   --out results/kill/k9
  ```
- **L3946**
  ```
  OLD: ./code/build/kill_k9b  --out results/kill/k9b
  NEW: julia --project=Kirigami Kirigami/apps/kill_k9b.jl  --out results/kill/k9b
  ```
- **L3947**
  ```
  OLD: ./code/build/kill_k9c  --out results/kill/k9c     # 12 shards; --aggregate to merge
  NEW: julia --project=Kirigami Kirigami/apps/kill_k9c.jl  --out results/kill/k9c     # 12 shards; --aggregate to merge
  ```
- **L3948**
  ```
  OLD: ./code/build/kill_native200 --n 200 --maxf 800 --timeout 600 \
  NEW: julia --project=Kirigami Kirigami/apps/kill_native200.jl --n 200 --maxf 800 --timeout 600 \
  ```

### `docs/paper/paper.tex` — 1 edits, 0 leave

- **L345**
  ```
  OLD: \plan{The algorithm, from \texttt{code/src/method/range\_embed.\{hpp,cpp\}},
  NEW: \plan{The algorithm, from \texttt{Kirigami/src/method/range\_embed.jl},
  ```

### `REPORT.md` — 35 edits, 0 leave

- **L42**
  ```
  OLD: | C9 | Two-stage range-maximising constrained embedding; 307 of 400 designs deployable against 0 for every baseline | **algorithm** | `results/final/figures/summary_table.md`, `results/kill/k9c/`, `code/src/method/range_embed.{hpp,cpp}` |
  NEW: | C9 | Two-stage range-maximising constrained embedding; ⟨JULIA:k9c:307 of 400⟩ designs deployable against 0 for every baseline | **algorithm** | `results/final/figures/summary_table.md`, `results/kill/k9c/`, `Kirigami/src/method/range_embed.jl` |
  ```
- **L271**
  ```
  OLD: `code/tests/test_design.cpp` at 1e−9. Exported closed, at `Θ_max/2` and at `0.9 Θ_max` as
  NEW: `Kirigami/test/test_design.jl` at 1e−9. Exported closed, at `Θ_max/2` and at `0.9 Θ_max` as
  ```
- **L316–317**
  ```
  OLD:   each time, ending with **zero unresolved disagreements**. `derivation_tests`: **33 cases /
         89,074 assertions, 0 failures** (re-run by the writer). **Disclosure:** MISSION §9's
  NEW:   each time, ending with **zero unresolved disagreements**. `derivation_tests.jl`: **⟨JULIA:derivation_tests⟩,
         0 failures** (re-run by the writer). **Disclosure:** MISSION §9's
  ```
- **L323**
  ```
  OLD: * **Gate 8 (Build).** `kiri_tests`: **103 cases / 17,504 assertions, 0 failures** (re-run by
  NEW: * **Gate 8 (Build).** `Pkg.test()`: **⟨JULIA:tests⟩, 0 failures** (re-run by
  ```
- **L531**
  ```
  OLD: | core pipeline (2026 reimplementation) | `code/src/core/{mesh,cut,holes,orientation,tutte_auxetic,kinematics,collision,optimize,rank_checks,generators}.{hpp,cpp}` |
  NEW: | core pipeline (2026 reimplementation) | `Kirigami/src/core/{mesh,cut,holes,orientation,tutte_auxetic,kinematics,collision,optimize,rank_checks,generators}.jl` |
  ```
- **L532**
  ```
  OLD: | method | `code/src/method/{deploy_basis,contact,design,range_opt,zero_plus,convex_embed,range_embed,mobility,expansive_cone,periodic_jacobian,budget}.{hpp,cpp}` |
  NEW: | method | `Kirigami/src/method/{deploy_basis,contact,design,range_opt,zero_plus,convex_embed,range_embed,mobility,expansive_cone,periodic_jacobian,budget}.jl` |
  ```
- **L533**
  ```
  OLD: | export | `code/src/export/{layout,solid,svg,stl,threemf,xml,zip,material}.{hpp,cpp}` |
  NEW: | export | `Kirigami/src/export/{layout,solid,svg,stl,threemf,xml,zip,material}.jl` |
  ```
- **L534**
  ```
  OLD: | tests | `code/tests/{test_mesh_cut,test_holes,test_system,test_kinematics,test_collision,test_method,test_design,test_range_embed,test_rank_checks,test_reference_cases,test_export,derivation_tests}.cpp` |
  NEW: | tests | `Kirigami/test/{test_mesh_cut,test_holes,test_system,test_kinematics,test_collision,test_method,test_design,test_range_embed,test_rank_checks,test_reference_cases,test_export,derivation_tests}.jl` |
  ```
- **L535**
  ```
  OLD: | CLIs | `code/apps/{kiri_gen,kiri_analyze,kiri_deploy,kiri_design,kiri_export,kiri_sweep,kiri_reference}.cpp` |
  NEW: | CLIs | `Kirigami/apps/{kiri_gen,kiri_analyze,kiri_deploy,kiri_design,kiri_export,kiri_sweep,kiri_reference}.jl` |
  ```
- **L536**
  ```
  OLD: | kill drivers | `code/apps/kill_{k1a,k1b,k1c,k2a,k2b,k2c,k3a,k5,k6,k7,k8a,k9,k9b,k9c,b3,b4,t1,jitter,e1,native200,f23}.cpp` |
  NEW: | kill drivers | `Kirigami/apps/kill_{k1a,k1b,k1c,k2a,k2b,k2c,k3a,k5,k6,k7,k8a,k9,k9b,k9c,b3,b4,t1,jitter,e1,native200,f23}.jl` |
  ```
- **L555**
  ```
  OLD: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
  NEW: julia --project=Kirigami -e 'using Pkg; Pkg.instantiate()'   # one-off; no build step
  ```
- **L561–562**
  ```
  OLD: ./code/build/kiri_tests           # 103 cases, 17504 assertions, 0 failures
       ./code/build/derivation_tests     # 33 cases, 89074 assertions, 0 failures
  NEW: julia --project=Kirigami -e 'using Pkg; Pkg.test()'   # ⟨JULIA:tests⟩, 0 failures; includes derivation_tests.jl (⟨JULIA:derivation_tests⟩)
  ```
- **L568**
  ```
  OLD: cmake --build code/build --target kill_e1
  NEW: # no build step; the driver is Kirigami/apps/kill_e1.jl (run_e1.sh invokes it)
  ```
- **L576**
  ```
  OLD: ./code/build/kill_k1b  --out results/kill/k1b          # ~34 s
  NEW: julia --project=Kirigami Kirigami/apps/kill_k1b.jl  --out results/kill/k1b          # ~34 s
  ```
- **L577**
  ```
  OLD: ./code/build/kill_k1a  --out results/kill/k1a          # ~40 s
  NEW: julia --project=Kirigami Kirigami/apps/kill_k1a.jl  --out results/kill/k1a          # ~40 s
  ```
- **L578**
  ```
  OLD: ./code/build/kill_k2a  --out results/kill/k2a          # ~7 s
  NEW: julia --project=Kirigami Kirigami/apps/kill_k2a.jl  --out results/kill/k2a          # ~7 s
  ```
- **L579**
  ```
  OLD: ./code/build/kill_k2c  --out results/kill/k2c          # ~3 min
  NEW: julia --project=Kirigami Kirigami/apps/kill_k2c.jl  --out results/kill/k2c          # ~3 min
  ```
- **L580**
  ```
  OLD: ./code/build/kill_k3a  --out results/kill/k3a          # ~25 min
  NEW: julia --project=Kirigami Kirigami/apps/kill_k3a.jl  --out results/kill/k3a          # ~25 min
  ```
- **L581**
  ```
  OLD: ./code/build/kill_k3a_recheck --out results/kill/k3a   # ~4 min
  NEW: julia --project=Kirigami Kirigami/apps/kill_k3a_recheck.jl --out results/kill/k3a   # ~4 min
  ```
- **L582**
  ```
  OLD: ./code/build/kill_k5   --out results/kill/k5 --n 200 --maxf 800   # ~5 min
  NEW: julia --project=Kirigami Kirigami/apps/kill_k5.jl   --out results/kill/k5 --n 200 --maxf 800   # ~5 min
  ```
- **L583**
  ```
  OLD: ./code/build/kill_k1c  --out results/kill/k1c          # ~31 s   (needs baseline/native)
  NEW: julia --project=Kirigami Kirigami/apps/kill_k1c.jl  --out results/kill/k1c          # ~31 s   (needs baseline/native)
  ```
- **L584**
  ```
  OLD: ./code/build/kill_k2b  --out results/kill/k2b --maxf 160          # ~26 min (needs baseline/native)
  NEW: julia --project=Kirigami Kirigami/apps/kill_k2b.jl  --out results/kill/k2b --maxf 160          # ~26 min (needs baseline/native)
  ```
- **L585**
  ```
  OLD: ./code/build/kill_f23  --out results/kill/f23          # ~1 s
  NEW: julia --project=Kirigami Kirigami/apps/kill_f23.jl  --out results/kill/f23          # ~1 s
  ```
- **L586**
  ```
  OLD: ./code/build/kill_k6   --out results/kill/k6
  NEW: julia --project=Kirigami Kirigami/apps/kill_k6.jl   --out results/kill/k6
  ```
- **L587**
  ```
  OLD: ./code/build/kill_k7   --out results/kill/k7
  NEW: julia --project=Kirigami Kirigami/apps/kill_k7.jl   --out results/kill/k7
  ```
- **L588**
  ```
  OLD: ./code/build/kill_t1   --out results/kill/t1           # ~7 s off the K6 cache
  NEW: julia --project=Kirigami Kirigami/apps/kill_t1.jl   --out results/kill/t1           # ~7 s off the K6 cache
  ```
- **L589**
  ```
  OLD: ./code/build/kill_b3   --out results/kill/b3
  NEW: julia --project=Kirigami Kirigami/apps/kill_b3.jl   --out results/kill/b3
  ```
- **L590**
  ```
  OLD: ./code/build/kill_b4   --out results/kill/b4
  NEW: julia --project=Kirigami Kirigami/apps/kill_b4.jl   --out results/kill/b4
  ```
- **L591**
  ```
  OLD: ./code/build/kill_k8a  --out results/kill/k8a
  NEW: julia --project=Kirigami Kirigami/apps/kill_k8a.jl  --out results/kill/k8a
  ```
- **L592**
  ```
  OLD: ./code/build/kill_k8a_recheck --out results/kill/k8a
  NEW: julia --project=Kirigami Kirigami/apps/kill_k8a_recheck.jl --out results/kill/k8a
  ```
- **L593**
  ```
  OLD: ./code/build/kill_jitter --out results/kill/jitter
  NEW: julia --project=Kirigami Kirigami/apps/kill_jitter.jl --out results/kill/jitter
  ```
- **L594**
  ```
  OLD: ./code/build/kill_k9   --out results/kill/k9
  NEW: julia --project=Kirigami Kirigami/apps/kill_k9.jl   --out results/kill/k9
  ```
- **L595**
  ```
  OLD: ./code/build/kill_k9b  --out results/kill/k9b
  NEW: julia --project=Kirigami Kirigami/apps/kill_k9b.jl  --out results/kill/k9b
  ```
- **L596**
  ```
  OLD: ./code/build/kill_k9c  --out results/kill/k9c          # 12 shards; --aggregate --nshards 12 to merge
  NEW: julia --project=Kirigami Kirigami/apps/kill_k9c.jl  --out results/kill/k9c          # 12 shards; --aggregate --nshards 12 to merge
  ```
- **L597**
  ```
  OLD: ./code/build/kill_native200 --out results/kill/native200           # in progress
  NEW: julia --project=Kirigami Kirigami/apps/kill_native200.jl --out results/kill/native200           # in progress
  ```

### `STATE.md` — 5 edits, 26 leave

- **L7**
  ```
  OLD: - 2026-09-05 00:30: ALL MISSION DELIVERABLES ON DISK AND COMMITTED. IDEA.md, derivations/ (6 rounds), code/ (109 tests / 17,854; derivation 33 / 89,074), results/ (400-design study, figures, E1 3,113), export/hero + hero2 (SVG/3MF/STL), REPORT.md, web explorer code/web. Open: none required by MISSION. Optional follow-ups: rerun Native200 timeouts with longer budget; fabricate hero2; manuscript.
  NEW: - 2026-09-05 00:30: ALL MISSION DELIVERABLES ON DISK AND COMMITTED. IDEA.md, derivations/ (6 rounds), code/ (⟨JULIA:tests⟩; derivation 33 / 89,074), results/ (400-design study, figures, E1 3,113), export/hero + hero2 (SVG/3MF/STL), REPORT.md, web explorer code/web. Open: none required by MISSION. Optional follow-ups: rerun Native200 timeouts with longer budget; fabricate hero2; manuscript.
  ```
- **L13**
  ```
  OLD: - Toolchain: clang++ (Apple), cmake 3.29.2, Eigen 5.0.1, nlohmann-json 3.12, doctest 2.5.3 in /opt/homebrew/include.
  NEW: - Toolchain: Julia ≥ 1.10 (`Kirigami/Project.toml`): LinearAlgebra, SparseArrays, StaticArrays, JSON, DelaunayTriangulation, Random, Statistics, Printf. No compiler.
  ```
- **L20** *leave* — dated log entry (history): `- [x] P2 Builder-Core → code/ (≈4k lines C++, 28 doctest cases / 6436 assertions). Gate 2 PASSED: orchestrator rebuilt from scratch (0 warni`
- **L25** *leave* — dated log entry (history): `- [x] P4 Ideator-rigidity → ideas/persona_rigidity.md (678 lines) + ideas/rigidity_rig_check.cpp (its own numeric checks on 21 graphs, compi`
- **L35** *leave* — dated log entry (history): `- [~] P4b Ideation round 2 (human directive): ideator2-theorist / ideator2-inverse / ideator2-theorist-fable pending; ideator2-theorist-b DO`
- **L36** *leave* — dated log entry (history): `- [x] experimenter-b4 → B4 FAIL (results/kill/b4/, KILL_REPORT §B4, kill_b4 app, 71 tests/16,488 assertions). Exact Θ_max>0 with boundary ro`
- **L37** *leave* — dated log entry (history): `- [x] critic-r2 → ideas/ranking_r2.md (539 lines). Top-5: X1 = A1 expansive-cone LP (+A4 hero, kill_k8a: LP over ker A at X_ini on 100 K1a g`
- **L38** *leave* — dated log entry (history): `- [x] checker-t5 → round 6 complete (derivation_tests 33 cases / 89,074 assertions; kiri_tests 16,622 — both re-run by orchestrator). Add-on`
- **L40** *leave* — dated log entry (history): `- [x] experimenter-t1 → kill_t1 FAIL (results/kill/t1/, KILL_REPORT §T-1, 7 s off K6 cache). 400 designs, 165,788 interior vertices, 14,886 `
- **L45** *leave* — dated log entry (history): `- [x] experimenter-jitter → A3 DONE (results/kill/jitter/, KILL_REPORT §A3, kill_jitter app, 65 tests/16,229 assertions). MIXED. Exact Θ_max`
- **L50** *leave* — dated log entry (history): `- [x] P7 Checker → derivations/check.md + code/tests/derivation_tests.cpp (standalone, 25 cases / 29,033 assertions pass). VERIFIED: T1 (+St`
- **L51** *leave* — dated log entry (history): `- [x] P8 Builder-Export → code/src/export (SVG laser layers cut/score/engrave; living-hinge neck & pin-pad solids; binary STL; 3MF with own `
- **L56** *leave* — dated log entry (history): `- [x] P7 Deriver round 5 + Checker round 5 (fresh, checker-6): ZERO unresolved disagreements: YES. derivation_tests.cpp: 31 cases / 89,055 a`
- **L58** *leave* — dated log entry (history): `- [x] Baseline parity: baseline/native builds the authors' code natively (libigl 2.5, Optiz, Eigen 3.4 via FetchContent; bind.cpp replaced b`
- **L67**
  ```
  OLD: - D3 (human directive 2026-09-03): ALL coding tasks in C++ (core, method, export, tests, experiment drivers, data generation). Python is allowed only for matplotlib plotting of CSV/JSON dumped by C++. Every Builder/Experimenter spec must state this.
  NEW: - D3 (human directive 2026-09-03, superseded 2026-09-19 by the Julia port): ALL coding tasks in Julia (core, method, export, tests, experiment drivers, data generation) in the `Kirigami` package. Python is allowed ONLY for matplotlib plotting of CSV/JSON dumped by the Julia apps.
  ```
- **L74** *leave* — dated log entry (history): `- D9 (orchestrator, 2026-09-04 17:15, resolves ESCALATION.md): constructive half = K9's convexity + split-inward constrained embedding (F36)`
- **L111** *leave* — dated log entry (history): `- F38 (checker-k8a, results/kill/k8a/recheck.md): cone_lp's primal witness came from an under-converged log-sum-exp loop while the Frank–Wol`
- **L112** *leave* — dated log entry (history): `- F37 FINAL (K9c, 400/400, orchestrator's own aggregation 20:35): deployable 307/400 (76.8%), refereed 307/307, certified 306, ε_max ≥ 0.1 o`
- **L113** *leave* — dated log entry (history): `- F36 (K9, convexity-constrained embedding, results/kill/k9/, K6's 400 designs; agent killed by session limit before writing §K9): FIRST NON`
- **L115** *leave* — dated log entry (history): `- F34 (checker-t5, B4 voronoi_93): the bisection referee theta_bisect (collision.cpp has_collision/polygons_overlap) is NOT ground truth. On`
- **L117** *leave* — dated log entry (history): `- F32 (A3 → checker-cert, FIXED): the certificate's NOROOT clause omitted the two interval tests (0 ≤ ⟨w−a,b−a⟩ ≤ |b−a|²) that the exact sca`
- **L123** *leave* — dated log entry (history): `- U10 RESOLVED (builder-method-2): the code is deterministic (no clocks/threads/unordered containers; Eigen single-threaded). The id-130 rep`
- **L163** *leave* — dated log entry (history): `- [~] WP1 Deriver-L DONE → derivations/lemmas.md (869 lines; L1.1–L1.4, L2.1–L2.3 proved; converse of L1.1 under H-L1 which is measurably fa`
- **L164** *leave* — dated log entry (history): `- [x] WP2 Experimenter-Y → results/yield/{YIELD.md, features.csv 400, fresh.csv 200, analyse.py, 3 figs}, code/apps/kill_yield.cpp. FAIL vs `
- **L165** *leave* — dated log entry (history): `- [x] WP3 CANCELLED (D15, human directive: no fabrication, theory only). Builder-fab stopped, export/fab and its apps removed, kiri_export.c`
- **L167** *leave* — dated log entry (history): `- [x] WP6b Propagator-A DONE: derivation_tests CMake target (code/CMakeLists.txt:68); test_reference_cases.cpp 204 lines, kiri_tests 109/17,`
- **L173**
  ```
  OLD: - E1 (environment, 2026-09-08): `cmake --build` / `make` / `lldb` fail (xcode-select shim error 34304). Use `/Applications/Xcode.app/Contents/Developer/usr/bin/make -C code/build <target>` directly. Load average 17–27 with 8 native processes + 2 agents.
  NEW: - E1 (environment, 2026-09-08): obsolete after the Julia port — no cmake/make/Xcode dependency remains.
  ```
- **L184**
  ```
  OLD: - Working tree clean at a18c000. baseline/native/build/tuttekiri_cli present. code/build/{kiri_tests,kill_native200} present. em++ at /opt/homebrew/bin. Bambu Studio installed; PrusaSlicer not. Load average 1.6 (machine idle). python3 matplotlib: use `arch -arm64 /usr/local/bin/python3` as before.
  NEW: - Working tree clean. baseline/native/build/tuttekiri_cli must be rebuilt locally if native comparisons are re-run (build tree not migrated). Julia apps run from source; no build products.
  ```
- **L192** *leave* — dated log entry (history): `- D15 (human directive 2026-09-08 21:40): NO fabrication. Mission 2 stays theory/computation only. WP3 removed from MISSION2.md; builder-fab`
- **L201** *leave* — dated log entry (history): `- 2026-09-08 22:22: load spike to 84 observed (regime native shards + basin + rerun); by 22:22 heavy processes ≈ 18 cores (17 tuttekiri_cli,`
- **L203** *leave* — dated log entry (history): `- 2026-09-08 22:44: USER STOPPED deriver-l (round 2 in progress, 5 R2 corrections already in lemmas.md), checker-l (idle), exp-regime (7a do`

### `IDEA.md` — 2 edits, 0 leave

- **L171**
  ```
  OLD: (`code/src/method/range_embed.{hpp,cpp}`). Same population, same shape space
  NEW: (`Kirigami/src/method/range_embed.jl`). Same population, same shape space
  ```
- **L214**
  ```
  OLD: regression in `code/tests/test_design.cpp`). Exported closed, at `Θ_max/2` and at `0.9 Θ_max`
  NEW: regression in `Kirigami/test/test_design.jl`). Exported closed, at `Θ_max/2` and at `0.9 Θ_max`
  ```

### `derivations/core.md` — 51 edits, 0 leave

- **L22–23**
  ```
  OLD: The Checker read this file with fresh context and wrote `derivations/check.md` (25 test cases,
       29 033 assertions, its own program `code/tests/derivation_tests.cpp`), raising eleven
  NEW: The Checker read this file with fresh context and wrote `derivations/check.md` (⟨JULIA:derivation_tests⟩,
       its own test file `Kirigami/test/derivation_tests.jl`), raising eleven
  ```
- **L25**
  ```
  OLD: `derivations/scratch/check_r2.cpp` (build line in its header; `graphs used: 16`, 197 shape-space
  NEW: `derivations/scratch/check_r2.jl` (run line in its header; `graphs used: 16`, 197 shape-space
  ```
- **L37**
  ```
  OLD: | **D8** | `contact.hpp`'s flat-centroid pruning is unsound in principle | **ADOPTED as a note.** My (T4.3) is the sound static test; the flat-centroid variant discards 19 842 pairs the sound test keeps (no observed wrong `Θ_max`, 96/96) | T4.5a |
  NEW: | **D8** | `contact.jl`'s flat-centroid pruning is unsound in principle | **ADOPTED as a note.** My (T4.3) is the sound static test; the flat-centroid variant discards 19 842 pairs the sound test keeps (no observed wrong `Θ_max`, 96/96) | T4.5a |
  ```
- **L38**
  ```
  OLD: | **D9** | `u_f` means different objects here and in `deploy_basis.hpp` | **CORRECTED** by a notation note fixing this file's `u_f` and naming the header's object explicitly. No code was touched (other agents own `code/`) | §0.10 |
  NEW: | **D9** | `u_f` means different objects here and in `deploy_basis.jl` | **CORRECTED** by a notation note fixing this file's `u_f` and naming the header's object explicitly. No code was touched (other agents own `code/`) | §0.10 |
  ```
- **L53**
  ```
  OLD: Round 3's program is `derivations/scratch/check_r3.cpp` (same 16-graph corpus, same seed, same 197
  NEW: Round 3's program is `derivations/scratch/check_r3.jl` (same 16-graph corpus, same seed, same 197
  ```
- **L58**
  ```
  OLD: | **R2.5** | a **third** structural class `g(0) = g′(0) = 0` (`p = −q`, `r = 0`, `h = p(1 − cos θ)`) is unnamed, and round 2's "0 mismatches" is not reproducible | **CONCEDED, and round 2's number is WITHDRAWN as circular** — `check_r2.cpp` compared the deflated atom list against a reference using the *same* deflation, so it was structurally blind to class 3. New: **Lemma T5.1e** proves `h = p(1 − cos θ)` has constant sign on `(0, π)`, hence is *never* a contact event, so the class deflates by `(1 − cos θ) ∝ τ²` and its atom is the constant `true` (T5.1e′). Re-measured against a `τ`-chart-free crossing test: **0 mismatches of 2 145 387 at every `ε ∈ {0.2, 0.02, 0.001}`**, versus **153** (all class 3, `ε`-independent) without the new clause | T5.2b.2, T3.H.5 |
  NEW: | **R2.5** | a **third** structural class `g(0) = g′(0) = 0` (`p = −q`, `r = 0`, `h = p(1 − cos θ)`) is unnamed, and round 2's "⟨JULIA:check_r3:mismatches⟩ mismatches" is not reproducible | **CONCEDED, and round 2's number is WITHDRAWN as circular** — `check_r2.jl` compared the deflated atom list against a reference using the *same* deflation, so it was structurally blind to class 3. New: **Lemma T5.1e** proves `h = p(1 − cos θ)` has constant sign on `(0, π)`, hence is *never* a contact event, so the class deflates by `(1 − cos θ) ∝ τ²` and its atom is the constant `true` (T5.1e′). Re-measured against a `τ`-chart-free crossing test: **⟨JULIA:check_r3:mismatches⟩ mismatches of 2 145 387 at every `ε ∈ {0.2, 0.02, 0.001}`**, versus **153** (all class 3, `ε`-independent) without the new clause | T5.2b.2, T3.H.5 |
  ```
- **L91**
  ```
  OLD: Deciding numbers for round 4: **unchanged** — `check_r3.cpp` still reports 0 mismatches of
  NEW: Deciding numbers for round 4: **unchanged** — `check_r3.jl` still reports ⟨JULIA:check_r3:mismatches⟩ mismatches of
  ```
- **L92**
  ```
  OLD: 2 145 387 at every `ε ∈ {0.2, 0.02, 0.001}`, and the standalone `code/tests/derivation_tests.cpp`
  NEW: 2 145 387 at every `ε ∈ {0.2, 0.02, 0.001}`, and the standalone `Kirigami/test/derivation_tests.jl`
  ```
- **L109**
  ```
  OLD: Deciding numbers for round 5: **unchanged**. `check_r3.cpp` still reports 0 mismatches of
  NEW: Deciding numbers for round 5: **unchanged**. `check_r3.jl` still reports ⟨JULIA:check_r3:mismatches⟩ mismatches of
  ```
- **L111**
  ```
  OLD: `code/tests/derivation_tests.cpp` re-runs at **30 test cases, 89 044 assertions, 0 failures**.
  NEW: `Kirigami/test/derivation_tests.jl` re-runs at **⟨JULIA:derivation_tests⟩, 0 failures**.
  ```
- **L129**
  ```
  OLD: (`cut.hpp`: `corner_to_prime`); write `(v, f)` for the copy of `v ∈ V` carried by face `f`, and
  NEW: (`cut.jl`: `corner_to_prime`); write `(v, f)` for the copy of `v ∈ V` carried by face `f`, and
  ```
- **L142**
  ```
  OLD: **0.7 [F, `code/src/core/kinematics.cpp`] The code's sign convention.** Face `f` is transformed by
  NEW: **0.7 [F, `Kirigami/src/core/kinematics.jl`] The code's sign convention.** Face `f` is transformed by
  ```
- **L171**
  ```
  OLD: **0.10 [D] The symbol `u_f`, and the clash with `deploy_basis.hpp` (D9).** In this file `u_f ∈ R²`
  NEW: **0.10 [D] The symbol `u_f`, and the clash with `deploy_basis.jl` (D9).** In this file `u_f ∈ R²`
  ```
- **L173**
  ```
  OLD: satisfies the increment law (T1.3). The header comment of `code/src/method/deploy_basis.hpp` writes
  NEW: satisfies the increment law (T1.3). The header comment of `Kirigami/src/method/deploy_basis.jl` writes
  ```
- **L185**
  ```
  OLD: `deploy_basis.hpp` should rename its comment variable to `ν_f` before either symbol becomes a
  NEW: `deploy_basis.jl` should rename its comment variable to `ν_f` before either symbol becomes a
  ```
- **L272**
  ```
  OLD: **[N] Path-independence check.** `derivations/scratch/check_t1_t2.cpp`, check **C1**: build `u` by BFS
  NEW: **[N] Path-independence check.** `derivations/scratch/check_t1_t2.jl`, check **C1**: build `u` by BFS
  ```
- **L307**
  ```
  OLD: BFS, against `kinematics.cpp::deploy()` at `θ ∈ {0, 0.05, 0.3, 0.77, 1.2, 1.9, 2.6, 3.0}` on the same
  NEW: BFS, against `kinematics.jl::deploy()` at `θ ∈ {0, 0.05, 0.3, 0.77, 1.2, 1.9, 2.6, 3.0}` on the same
  ```
- **L390**
  ```
  OLD: Re-run `derivations/scratch/check_t1_t2.cpp`. Tolerance `1e−12` on C1–C6, `1e−9` on C7. Additionally:
  NEW: Re-run `derivations/scratch/check_t1_t2.jl`. Tolerance `1e−12` on C1–C6, `1e−9` on C7. Additionally:
  ```
- **L497**
  ```
  OLD: round-off). Independent of `check_t1_t2.cpp`'s C6, which tests the pin equations rather than the
  NEW: round-off). Independent of `check_t1_t2.jl`'s C6, which tests the pin equations rather than the
  ```
- **L566**
  ```
  OLD: matching `deploy_basis.hpp::harmonic_roots`. Each of `p, q, r` is quadratic in `t`, so the
  NEW: matching `deploy_basis.jl::harmonic_roots`. Each of `p, q, r` is quadratic in `t`, so the
  ```
- **L629**
  ```
  OLD: quantity `collision.hpp::theta_max` estimates by grid scan + bisection with faces shrunk by a
  NEW: quantity `collision.jl::theta_max` estimates by grid scan + bisection with faces shrunk by a
  ```
- **L687**
  ```
  OLD: **[N] Measured (`derivations/scratch/check_t4_t5.cpp`, checks D1/D2).** 16 valid embeddings
  NEW: **[N] Measured (`derivations/scratch/check_t4_t5.jl`, checks D1/D2).** 16 valid embeddings
  ```
- **L782**
  ```
  OLD: **[N] Measured (`derivations/scratch/check_r2.cpp`, R2-C).** On all **8 split-free** patterns
  NEW: **[N] Measured (`derivations/scratch/check_r2.jl`, R2-C).** On all **8 split-free** patterns
  ```
- **L822**
  ```
  OLD: faces actually approach can be discarded. `code/src/method/contact.hpp::candidate_pairs(...,
  NEW: faces actually approach can be discarded. `Kirigami/src/method/contact.jl::candidate_pairs(...,
  ```
- **L841**
  ```
  OLD: which no pruning uses. K2c's spec and `contact.hpp::swept_discs` both work **in the face's own
  NEW: which no pruning uses. K2c's spec and `contact.jl::swept_discs` both work **in the face's own
  ```
- **L862**
  ```
  OLD: **[N] Measured (`derivations/scratch/check_r2.cpp`, R2-A).** Over **1 628** `M′`-copies on 16
  NEW: **[N] Measured (`derivations/scratch/check_r2.jl`, R2-A).** Over **1 628** `M′`-copies on 16
  ```
- **L870**
  ```
  OLD: 1. Of the two radii in `contact.hpp`, `rho` (the `√(‖x‖²+‖χ‖²)` form, commented "sound bound") is
  NEW: 1. Of the two radii in `contact.jl`, `rho` (the `√(‖x‖²+‖χ‖²)` form, commented "sound bound") is
  ```
- **L933**
  ```
  OLD: * Re-run `derivations/scratch/check_t4_t5.cpp` (rules: D1 `≤ 1e−5`, D2 `≤ 1e−12`) **and**
  NEW: * Re-run `derivations/scratch/check_t4_t5.jl` (rules: D1 `≤ 1e−5`, D2 `≤ 1e−12`) **and**
  ```
- **L934**
  ```
  OLD:   `derivations/scratch/check_r2.cpp` (rules: R2-A `≤ 1e−12`, R2-C `≤ 1e−12`, R2-D violations `= 0`;
  NEW:   `derivations/scratch/check_r2.jl` (rules: R2-A `≤ 1e−12`, R2-C `≤ 1e−12`, R2-D violations `= 0`;
  ```
- **L935**
  ```
  OLD:   its R2-E line is **withdrawn**, see T5.2b.2) **and** `derivations/scratch/check_r3.cpp`
  NEW:   its R2-E line is **withdrawn**, see T5.2b.2) **and** `derivations/scratch/check_r3.jl`
  ```
- **L1082**
  ```
  OLD: the negation of `contact.hpp::penetrates_immediately()`. `EMB` is what the proof needs, but it is a
  NEW: the negation of `contact.jl::penetrates_immediately()`. `EMB` is what the proof needs, but it is a
  ```
- **L1501**
  ```
  OLD: for it. Both `derivations/scratch/check_r3.cpp` and `code/tests/derivation_tests.cpp` implement
  NEW: for it. Both `derivations/scratch/check_r3.jl` and `Kirigami/test/derivation_tests.jl` implement
  ```
- **L1523**
  ```
  OLD: **[N] Measured, round 3 (`derivations/scratch/check_r3.cpp`).** Round 2's number
  NEW: **[N] Measured, round 3 (`derivations/scratch/check_r3.jl`).** Round 2's number
  ```
- **L1524**
  ```
  OLD: ("**0 / 2 363 380** at every `ε`") was **not a valid measurement and is withdrawn**: `check_r2.cpp`
  NEW: ("**0 / 2 363 380** at every `ε`") was **not a valid measurement and is withdrawn**: `check_r2.jl`
  ```
- **L1526**
  ```
  OLD: comparison could not see class 3 at all. `check_r3.cpp` fixes this by computing the truth **without
  NEW: comparison could not see class 3 at all. `check_r3.jl` fixes this by computing the truth **without
  ```
- **L1706**
  ```
  OLD: * **Round 3:** re-run `derivations/scratch/check_r3.cpp` with several `ε`. The **three-class**
  NEW: * **Round 3:** re-run `derivations/scratch/check_r3.jl` with several `ε`. The **three-class**
  ```
- **L1711**
  ```
  OLD:   Re-run `derivations/scratch/check_r2.cpp` for R2-A/R2-C/R2-D; R2-D must report **0** violations;
  NEW:   Re-run `derivations/scratch/check_r2.jl` for R2-A/R2-C/R2-D; R2-D must report **0** violations;
  ```
- **L1793**
  ```
  OLD: 5. **Independence** — the same `Θ_max` cross-checked against `collision.hpp`'s bisection, which uses
  NEW: 5. **Independence** — the same `Θ_max` cross-checked against `collision.jl`'s bisection, which uses
  ```
- **L1984**
  ```
  OLD: | `derivations/scratch/check_t1_t2.cpp` | `clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 -I code/src …` linking `code/src/core/*.cpp` | C1 1.60e−14, C2 1.42e−14, C3 9.73e−14, C4 2.08e−14, C5 4.46e−14, C6 1.83e−14, C7 1.68e−10 on 30 graphs × 8 angles |
  NEW: | `derivations/scratch/check_t1_t2.jl` | `julia --project=Kirigami derivations/scratch/check_t1_t2.jl` linking `Kirigami/src/core/*.jl` | C1 ⟨JULIA:check_t1_t2:1.60e−14⟩, C2 ⟨JULIA:check_t1_t2:1.42e−14⟩, C3 ⟨JULIA:check_t1_t2:9.73e−14⟩, C4 ⟨JULIA:check_t1_t2:2.08e−14⟩, C5 ⟨JULIA:check_t1_t2:4.46e−14⟩, C6 ⟨JULIA:check_t1_t2:1.83e−14⟩, C7 ⟨JULIA:check_t1_t2:1.68e−10⟩ on 30 graphs × 8 angles |
  ```
- **L1985**
  ```
  OLD: | `derivations/scratch/check_t4_t5.cpp` | same | D1 2.09e−6 (rule 1e−5), D2 8.88e−16, D3 4.37e−15 over 72 split edges, D4 719/1061 violations in the **raw** frame, ratio 1.4141 — see T4.5b: this says nothing about any pruning |
  NEW: | `derivations/scratch/check_t4_t5.jl` | same | D1 ⟨JULIA:check_t4_t5:2.09e−6⟩ (rule 1e−5), D2 ⟨JULIA:check_t4_t5:8.88e−16⟩, D3 ⟨JULIA:check_t4_t5:4.37e−15⟩ over 72 split edges, D4 719/1061 violations in the **raw** frame, ratio 1.4141 — see T4.5b: this says nothing about any pruning |
  ```
- **L1986**
  ```
  OLD: | `derivations/scratch/check_r2.cpp` (round 2) | same, plus `code/src/core/*.cpp`; build line in the file header, takes `ε` as `argv[1]` | R2-A 3.11e−15 (1 628 copies), R2-B drift/`r_f` 4.62→20.03 over diameters 5.66→19.80, R2-C 8.88e−16 (8 split-free patterns), R2-D 0 violations at every `ε ∈ {0.2, 0.05, 0.02, 0.005, 0.001}` on 197 samples. **R2-E is withdrawn** — its deflated reference applied the same deflation as the atom list under test, so the comparison was circular (T5.2b.2) |
  NEW: | `derivations/scratch/check_r2.jl` (round 2) | same, plus `Kirigami/src/core/*.jl`; run line in the file header, takes `ε` as `ARGS[1]` | R2-A ⟨JULIA:check_r2:3.11e−15⟩ (1 628 copies), R2-B drift/`r_f` 4.62→20.03 over diameters 5.66→19.80, R2-C ⟨JULIA:check_r2:8.88e−16⟩ (8 split-free patterns), R2-D ⟨JULIA:check_r2:violations⟩ violations at every `ε ∈ {0.2, 0.05, 0.02, 0.005, 0.001}` on 197 samples. **R2-E is withdrawn** — its deflated reference applied the same deflation as the atom list under test, so the comparison was circular (T5.2b.2) |
  ```
- **L1987**
  ```
  OLD: | `derivations/scratch/check_r3.cpp` (round 3) | same corpus and seed; build line in the file header, takes `ε` as `argv[1]` | three-class deflated atom list vs a `τ`-chart-free crossing test: **0** mismatches of 2 145 387 at `ε ∈ {0.2, 0.02, 0.001}` (47 ambiguous, skipped); the round-2 two-class list gets **153** wrong, all in class 3, `ε`-independent; class sizes 2 070 850 / 73 446 / 1 138 |
  NEW: | `derivations/scratch/check_r3.jl` (round 3) | same corpus and seed; run line in the file header, takes `ε` as `ARGS[1]` | three-class deflated atom list vs a `τ`-chart-free crossing test: **0** mismatches of 2 145 387 at `ε ∈ {0.2, 0.02, 0.001}` (47 ambiguous, skipped); the round-2 two-class list gets **153** wrong, all in class 3, `ε`-independent; class sizes 2 070 850 / 73 446 / 1 138 |
  ```
- **L1991**
  ```
  OLD: `code/tests/derivation_tests.cpp` (25 cases, 29 033 assertions, 0 failures), was written
  NEW: `Kirigami/test/derivation_tests.jl` (⟨JULIA:derivation_tests⟩, 0 failures), was written
  ```
- **L2037**
  ```
  OLD: exactly at the far endpoint, `derivation_tests` R6-a). `contact.cpp` implements `E2` inclusively at
  NEW: exactly at the far endpoint, `derivation_tests.jl` R6-a). `contact.jl` implements `E2` inclusively at
  ```
- **L2075**
  ```
  OLD: > Measured (`derivation_tests` R6-a1/R6-a2, 9 642 Case-A substitute pairs over the corpus with
  NEW: > Measured (`derivation_tests.jl` R6-a1/R6-a2, 9 642 Case-A substitute pairs over the corpus with
  ```
- **L2148**
  ```
  OLD: not as separate atoms whose own roots are hunted. This is what `contact.cpp` implements, so the
  NEW: not as separate atoms whose own roots are hunted. This is what `contact.jl` implements, so the
  ```
- **L2187**
  ```
  OLD: **[N] Measured counterexample (`derivations/scratch/check_t5_adm.cpp`, and `derivation_tests` R6-b).**
  NEW: **[N] Measured counterexample (`derivations/scratch/check_t5_adm.jl`, and `derivation_tests.jl` R6-b).**
  ```
- **L2231**
  ```
  OLD: `derivations/scratch/check_b4_93.cpp`: **the scan is right.**
  NEW: `derivations/scratch/check_b4_93.jl`: **the scan is right.**
  ```
- **L2261**
  ```
  OLD: | `derivations/scratch/check_t5_adm.cpp` (new; build line in its header) | A: 1 363 Case-A pairs on the 7 tilings with jitter, ratio `∈ [0.5404, 1.0]`, `|s/‖e‖² − L′/L| ≤ 1.55e−15`, **0** admissibility failures. B: 1 completeness counterexample (`hexagons`, `ε = 1.570796`) |
  NEW: | `derivations/scratch/check_t5_adm.jl` (new; run line in its header) | A: 1 363 Case-A pairs on the 7 tilings with jitter, ratio `∈ [0.5404, 1.0]`, `|s/‖e‖² − L′/L| ≤ ⟨JULIA:check_t5_adm:1.55e−15⟩`, **0** admissibility failures. B: 1 completeness counterexample (`hexagons`, `ε = 1.570796`) |
  ```
- **L2262**
  ```
  OLD: | `derivations/scratch/check_b4_93.cpp` (new; build line in its header) | replays B4 ids 93 and 96 from the K6/B4 caches; 93 reproduces scan `0.248400` vs referee `0`, and localizes the referee's misfire to the shared hinge vertex of faces (94, 184) |
  NEW: | `derivations/scratch/check_b4_93.jl` (new; run line in its header) | replays B4 ids 93 and 96 from the K6/B4 caches; 93 reproduces scan `0.248400` vs referee `0`, and localizes the referee's misfire to the shared hinge vertex of faces (94, 184) |
  ```
- **L2263**
  ```
  OLD: | `code/tests/derivation_tests.cpp` R6 + R6-c cases (new) | R6-a1 **0** failures / 9 642 · R6-a2 `1.77e−14` · R6-a3 `3.53e−14` · R6-b counterexample found · R6-c1 **0** common interior points of 1 442 401 · R6-c2 misfire pinned. Whole file: **33 test cases, 89 074 assertions, 0 failures** |
  NEW: | `Kirigami/test/derivation_tests.jl` R6 + R6-c cases (new) | R6-a1 ⟨JULIA:derivation_tests:R6-a1⟩ failures · R6-a2 `⟨JULIA:derivation_tests:1.77e−14⟩` · R6-a3 `⟨JULIA:derivation_tests:3.53e−14⟩` · R6-b counterexample found · R6-c1 ⟨JULIA:derivation_tests:R6-c1⟩ common interior points of 1 442 401 · R6-c2 misfire pinned. Whole file: **⟨JULIA:derivation_tests⟩, 0 failures** |
  ```

### `derivations/check.md` — 41 edits, 0 leave

- **L5**
  ```
  OLD: re-derived by hand from the code's own sign convention (`code/src/core/kinematics.cpp`), and
  NEW: re-derived by hand from the code's own sign convention (`Kirigami/src/core/kinematics.jl`), and
  ```
- **L7**
  ```
  OLD: `code/tests/derivation_tests.cpp`, which recomputes the face potential `u`, the harmonic
  NEW: `Kirigami/test/derivation_tests.jl`, which recomputes the face potential `u`, the harmonic
  ```
- **L9**
  ```
  OLD: `Γ`-cycle closure rows **without using the Deriver's scratch programs** (`check_t1_t2.cpp`,
  NEW: `Γ`-cycle closure rows **without using the Deriver's scratch programs** (`check_t1_t2.jl`,
  ```
- **L10**
  ```
  OLD: `check_t4_t5.cpp` were never compiled or run by me).
  NEW: `check_t4_t5.jl` were never compiled or run by me).
  ```
- **L15–17**
  ```
  OLD: clang++ -std=c++20 -O2 -arch arm64 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
         -Icode/src code/tests/derivation_tests.cpp code/build/libkiri_core.a \
         -o code/build/derivation_tests && ./code/build/derivation_tests
  NEW: julia --project=Kirigami -e 'using Pkg; Pkg.test(test_args=["derivation_tests"])'
  ```
- **L20**
  ```
  OLD: `code/CMakeLists.txt` and every existing file under `code/` are untouched. Result:
  NEW: `Kirigami/Project.toml` and every existing file under `code/` are untouched. Result:
  ```
- **L65**
  ```
  OLD: | 0.2, 0.4–0.6 | AGREE | definitions; match `cut.hpp` / `mesh.hpp` verbatim |
  NEW: | 0.2, 0.4–0.6 | AGREE | definitions; match `cut.jl` / `mesh.jl` verbatim |
  ```
- **L66**
  ```
  OLD: | 0.7, 0.7′ | AGREE | `kinematics.cpp` sets `a = -sigma[f]*theta*0.5`; `R(a) = cI + sJ` expands to `cI − σ_f s J`. Confirmed numerically: the closed form reproduces `deploy()` to 1.4e−14, the sign-flipped form fails by 45 (test **T1-a**) |
  NEW: | 0.7, 0.7′ | AGREE | `kinematics.jl` sets `a = -sigma[f]*theta*0.5`; `R(a) = cI + sJ` expands to `cI − σ_f s J`. Confirmed numerically: the closed form reproduces `deploy()` to ⟨JULIA:check_t1_t2:1.4e−14⟩, the sign-flipped form fails by 45 (test **T1-a**) |
  ```
- **L93**
  ```
  OLD: | T2.1 (T2.2) | AGREE | the pin constraint and the circulation formulation are the standard body-and-pin model; matches `mobility.hpp::build_A` |
  NEW: | T2.1 (T2.2) | AGREE | the pin constraint and the circulation formulation are the standard body-and-pin model; matches `mobility.jl::build_A` |
  ```
- **L117**
  ```
  OLD: | Thm T4.2″ (T4.1) | **AGREE** | I checked the one gap the statement could have: an *isolated* overlap angle. Interior overlap is an **open** condition in `θ` (an open ball inside both interiors persists under a continuous motion), so the overlap set has no isolated points and `Θ_max = θ_{i*}` exactly. **Task 1(a): the graze treatment is right.** Independently reproduced with my own separating-axis overlap test (not `collision.hpp`): hexagons `θ₁ = 1.047198`, **zero interior overlap at 199 sample angles across `(0, 2.094395)`**, overlap immediately after. `min`-over-roots would report `π/3`, the true range is `2π/3` |
  NEW: | Thm T4.2″ (T4.1) | **AGREE** | I checked the one gap the statement could have: an *isolated* overlap angle. Interior overlap is an **open** condition in `θ` (an open ball inside both interiors persists under a continuous motion), so the overlap set has no isolated points and `Θ_max = θ_{i*}` exactly. **Task 1(a): the graze treatment is right.** Independently reproduced with my own separating-axis overlap test (not `collision.jl`): hexagons `θ₁ = 1.047198`, **zero interior overlap at 199 sample angles across `(0, 2.094395)`**, overlap immediately after. `min`-over-roots would report `π/3`, the true range is `2π/3` |
  ```
- **L166**
  ```
  OLD: **in the face's own frame**", and `contact.hpp::swept_discs` implements exactly that, with
  NEW: **in the face's own frame**", and `contact.jl::swept_discs` implements exactly that, with
  ```
- **L191**
  ```
  OLD: 1. `contact.hpp`'s `rho` (the `√2` form, commented "sound bound") is the one that is loose by
  NEW: 1. `contact.jl`'s `rho` (the `√2` form, commented "sound bound") is the one that is loose by
  ```
- **L258**
  ```
  OLD: `contact.hpp` carries a separate `penetrates_immediately()`, and exactly the failure mode
  NEW: `contact.jl` carries a separate `penetrates_immediately()`, and exactly the failure mode
  ```
- **L322**
  ```
  OLD: ### D8 — `contact.hpp`'s flat-centroid pruning is unsound in principle (the Deriver's (T4.3) is not)
  NEW: ### D8 — `contact.jl`'s flat-centroid pruning is unsound in principle (the Deriver's (T4.3) is not)
  ```
- **L325**
  ```
  OLD: sound. `contact.hpp::candidate_pairs(..., use_static = true)` instead uses the **flat**
  NEW: sound. `contact.jl::candidate_pairs(..., use_static = true)` instead uses the **flat**
  ```
- **L333**
  ```
  OLD: ### D9 — documentation inconsistency in `deploy_basis.hpp`
  NEW: ### D9 — documentation inconsistency in `deploy_basis.jl`
  ```
- **L374**
  ```
  OLD: All from `code/tests/derivation_tests.cpp`; 25 test cases, 29 033 assertions, 0 failures.
  NEW: All from `Kirigami/test/derivation_tests.jl`; ⟨JULIA:derivation_tests⟩, 0 failures.
  ```
- **L494**
  ```
  OLD: this file's round-1 verdicts, `derivations/scratch/check_r2.cpp`, `code/README.md`.
  NEW: this file's round-1 verdicts, `derivations/scratch/check_r2.jl`, `Kirigami/README.md`.
  ```
- **L495**
  ```
  OLD: New program work: two test cases appended to `code/tests/derivation_tests.cpp` (that file only;
  NEW: New program work: two test cases appended to `Kirigami/test/derivation_tests.jl` (that file only;
  ```
- **L496**
  ```
  OLD: `CMakeLists.txt` and every other source untouched). Rebuilt with the command in its header:
  NEW: `Project.toml` and every other source untouched). Rebuilt with the command in its header:
  ```
- **L499–502**
  ```
  OLD: clang++ -std=c++20 -O2 -arch arm64 -I/opt/homebrew/include \
         -I/opt/homebrew/include/eigen3 -Icode/src \
         code/tests/derivation_tests.cpp code/build/libkiri_core.a \
         -o code/build/derivation_tests && ./code/build/derivation_tests
  NEW: julia --project=Kirigami -e 'using Pkg; Pkg.test(test_args=["derivation_tests"])'
  ```
- **L627**
  ```
  OLD: ## R2.6 New tests (appended to `code/tests/derivation_tests.cpp`)
  NEW: ## R2.6 New tests (appended to `Kirigami/test/derivation_tests.jl`)
  ```
- **L674**
  ```
  OLD: T5.2b′, T5.2b.1, T5.2b.2, T6.4; `derivations/scratch/check_r3.cpp`. New tests appended to
  NEW: T5.2b′, T5.2b.1, T5.2b.2, T6.4; `derivations/scratch/check_r3.jl`. New tests appended to
  ```
- **L675**
  ```
  OLD: `code/tests/derivation_tests.cpp` (that file only): `R3 Lemma T5.1e and the three-class atom
  NEW: `Kirigami/test/derivation_tests.jl` (that file only): `R3 Lemma T5.1e and the three-class atom
  ```
- **L705**
  ```
  OLD: **`check_r3.cpp` rebuilt and rerun.** Its numbers reproduce exactly, at all three `ε`:
  NEW: **`check_r3.jl` rebuilt and rerun.** Its numbers reproduce exactly, at all three `ε`:
  ```
- **L720**
  ```
  OLD: **not blind to class 3** in the way `check_r2.cpp` was. R2.5's circularity objection is answered.
  NEW: **not blind to class 3** in the way `check_r2.jl` was. R2.5's circularity objection is answered.
  ```
- **L748**
  ```
  OLD: Both `check_r3.cpp` and my own test exclude them *before* classifying (`h.scale() ≤ tol → skip`,
  NEW: Both `check_r3.jl` and my own test exclude them *before* classifying (`h.scale() ≤ tol → skip`,
  ```
- **L831**
  ```
  OLD: Appended to `code/tests/derivation_tests.cpp` (only that file; no other file touched), built with
  NEW: Appended to `Kirigami/test/derivation_tests.jl` (only that file; no other file touched), built with
  ```
- **L878**
  ```
  OLD: `R4 Sub-lemma T5.2b'' Case A ...` appended to `code/tests/derivation_tests.cpp` (that file only;
  NEW: `R4 Sub-lemma T5.2b'' Case A ...` appended to `Kirigami/test/derivation_tests.jl` (that file only;
  ```
- **L1035**
  ```
  OLD: and the zero set of h_o,pi' in (0,pi)` appended to `code/tests/derivation_tests.cpp` (that file only;
  NEW: and the zero set of h_o,pi' in (0,pi)` appended to `Kirigami/test/derivation_tests.jl` (that file only;
  ```
- **L1164**
  ```
  OLD: `code/src/method/contact.cpp` (`validity_certificate`, `contact_angles`) with the F32 doctest in
  NEW: `Kirigami/src/method/contact.jl` (`validity_certificate`, `contact_angles`) with the F32 unit test in
  ```
- **L1165**
  ```
  OLD: `code/tests/test_method.cpp`. Two new programs: `derivations/scratch/check_t5_adm.cpp` (standalone,
  NEW: `Kirigami/test/test_method.jl`. Two new programs: `derivations/scratch/check_t5_adm.jl` (standalone,
  ```
- **L1166**
  ```
  OLD: build line in its header) and the `R6` test case appended to `code/tests/derivation_tests.cpp`
  NEW: run line in its header) and the `R6` test case appended to `Kirigami/test/derivation_tests.jl`
  ```
- **L1167**
  ```
  OLD: (that file only), plus `derivations/scratch/check_b4_93.cpp` and the `R6-c` case for the B4 add-on.
  NEW: (that file only), plus `derivations/scratch/check_b4_93.jl` and the `R6-c` case for the B4 add-on.
  ```
- **L1208**
  ```
  OLD:    The segment of T4.1b is **closed**, so equality is admissible; `contact.cpp`'s inclusive
  NEW:    The segment of T4.1b is **closed**, so equality is admissible; `contact.jl`'s inclusive
  ```
- **L1239**
  ```
  OLD: to `NOROOT`". `contact.cpp` never implemented them. **The code was right and [D3] was unnecessary:
  NEW: to `NOROOT`". `contact.jl` never implemented them. **The code was right and [D3] was unnecessary:
  ```
- **L1271**
  ```
  OLD: **Counterexample, measured** (`check_t5_adm.cpp` §B, `derivation_tests` R6-b): `hexagons` at `X₀` has
  NEW: **Counterexample, measured** (`check_t5_adm.jl` §B, `derivation_tests.jl` R6-b): `hexagons` at `X₀` has
  ```
- **L1292**
  ```
  OLD: Replayed through `kill_b4.cpp`'s own pipeline (same graph, same `sigma_def`, same `free` system, same
  NEW: Replayed through `kill_b4.jl`'s own pipeline (same graph, same `sigma_def`, same `free` system, same
  ```
- **L1293**
  ```
  OLD: seeds and weights) by `derivations/scratch/check_b4_93.cpp`. The row reproduces exactly:
  NEW: seeds and weights) by `derivations/scratch/check_b4_93.jl`. The row reproduces exactly:
  ```
- **L1339**
  ```
  OLD: its centroid, or exclude edge pairs that share an endpoint) is a change to `code/src/core/collision.cpp`
  NEW: its centroid, or exclude edge pairs that share an endpoint) is a change to `Kirigami/src/core/collision.jl`
  ```
- **L1342**
  ```
  OLD: **Regression case added.** `derivation_tests` `R6-c` pins the reproduction with the two polygons as
  NEW: **Regression case added.** `derivation_tests.jl` `R6-c` pins the reproduction with the two polygons as
  ```

### `derivations/check_lemmas.md` — 24 edits, 0 leave

- **L10**
  ```
  OLD:    and the code (`deploy_basis.hpp`, `contact.hpp`, `zero_plus.hpp`, `periodic_jacobian.hpp`,
  NEW:    and the code (`deploy_basis.jl`, `contact.jl`, `zero_plus.jl`, `periodic_jacobian.jl`,
  ```
- **L11**
  ```
  OLD:    `apps/kill_b3.cpp::achievable/shape_point`, `tests/derivation_tests.cpp`).
  NEW:    `apps/kill_b3.jl::achievable/shape_point`, `test/derivation_tests.jl`).
  ```
- **L83**
  ```
  OLD: identity `zero_plus.hpp` states as `dC_e = 0` and it is here proved for hinge- and vertex-only
  NEW: identity `zero_plus.jl` states as `dC_e = 0` and it is here proved for hinge- and vertex-only
  ```
- **L127**
  ```
  OLD: `A ≠ 0`). The two interval predicates of `contact.hpp` are only ever evaluated *at* a root, so no
  NEW: `A ≠ 0`). The two interval predicates of `contact.jl` are only ever evaluated *at* a root, so no
  ```
- **L181**
  ```
  OLD:   `q_e = 0` degenerate case `zero_plus.hpp` is written to detect — so it is a real, describable
  NEW:   `q_e = 0` degenerate case `zero_plus.jl` is written to detect — so it is a real, describable
  ```
- **L218**
  ```
  OLD: This also matches the code: `kill_b3.cpp::achievable` reads `D(i,0..1) = ½ (M_{2i} T)(row 2)`,
  NEW: This also matches the code: `kill_b3.jl::achievable` reads `D(i,0..1) = ½ (M_{2i} T)(row 2)`,
  ```
- **L383**
  ```
  OLD: | 1 (`dC = 0`, `Y_w − Y_{a*} = s·dS`) | **AGREE** | matches `zero_plus.hpp`'s `dC_e = 0`, correctly generalised |
  NEW: | 1 (`dC = 0`, `Y_w − Y_{a*} = s·dS`) | **AGREE** | matches `zero_plus.jl`'s `dC_e = 0`, correctly generalised |
  ```
- **L433**
  ```
  OLD:    `zero_plus.hpp`/K5/F30 report as *failing* on measured designs. So the hypothesis is not free and
  NEW:    `zero_plus.jl`/K5/F30 report as *failing* on measured designs. So the hypothesis is not free and
  ```
- **L554**
  ```
  OLD: `achievable()` in `kill_k7.cpp` / `kill_b3.cpp` **never reads `consistency`**. And note why the
  NEW: `achievable()` in `kill_k7.jl` / `kill_b3.jl` **never reads `consistency`**. And note why the
  ```
- **L614**
  ```
  OLD: Four `TEST_CASE`s appended to `code/tests/derivation_tests.cpp` — that file only. Nothing under
  NEW: Four `@testset`s appended to `Kirigami/test/derivation_tests.jl` — that file only. Nothing under
  ```
- **L615**
  ```
  OLD: `code/src`, `code/apps`, `code/CMakeLists.txt`, `derivations/lemmas.md` or `derivations/core.md`
  NEW: `Kirigami/src`, `Kirigami/apps`, `Kirigami/Project.toml`, `derivations/lemmas.md` or `derivations/core.md`
  ```
- **L616**
  ```
  OLD: was touched. The file stays standalone doctest.
  NEW: was touched. The file stays standalone unit test.
  ```
- **L623**
  ```
  OLD: Two additions to the file's includes: `method/periodic_jacobian.hpp` and `apps/kill_common.hpp`
  NEW: Two additions to the file's includes: `method/periodic_jacobian.jl` and `apps/kill_common.jl`
  ```
- **L626**
  ```
  OLD: and `achievable()` are copied from `code/apps/kill_k7.cpp` / `kill_b3.cpp` so the numbers compare
  NEW: and `achievable()` are copied from `Kirigami/apps/kill_k7.jl` / `kill_b3.jl` so the numbers compare
  ```
- **L632–633**
  ```
  OLD: The build line in `derivations/check.md`'s header, plus **`-Icode/apps`** (the only change), and
       with the Xcode toolchain compiler by its full path, because a bare `clang++` aborts on this machine
  NEW: `derivation_tests.jl` is part of the package test suite; it is run with (no `-I` flags, no compiler --
  ```
- **L637–642**
  ```
  OLD: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang++ \
         -std=c++20 -O2 -arch arm64 \
         -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk \
         -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 -Icode/src -Icode/apps \
         code/tests/derivation_tests.cpp code/build/libkiri_core.a \
         -o code/build/derivation_tests && ./code/build/derivation_tests
  NEW: julia --project=Kirigami -e 'using Pkg; Pkg.test(test_args=["derivation_tests"])'
  ```
- **L648**
  ```
  OLD: [doctest] test cases:     37 |     37 passed | 0 failed | 0 skipped
  NEW: ⟨JULIA:derivation_tests⟩  (test cases)
  ```
- **L649**
  ```
  OLD: [doctest] assertions: 149115 | 149115 passed | 0 failed |
  NEW: ⟨JULIA:derivation_tests⟩  (assertions)
  ```
- **L650**
  ```
  OLD: [doctest] Status: SUCCESS!
  NEW: Test Summary: | Pass  Total  -- 0 failures
  ```
- **L694**
  ```
  OLD: `kill_common.hpp::reference_cases()` plus `make_graph(id, 18, 46, 220)` for `id = 0 …` until 60 are
  NEW: `kill_common.jl::reference_cases (population frozen in data/corpus/)()` plus `make_graph(id, 18, 46, 220)` for `id = 0 …` until 60 are
  ```
- **L756**
  ```
  OLD: **Corpus.** The **33 K7 patterns** (`l2_k7_population()`, byte-identical to `kill_k7.cpp`'s
  NEW: **Corpus.** The **33 K7 patterns** (`l2_k7_population()`, byte-identical to `kill_k7.jl`'s
  ```
- **L802**
  ```
  OLD: exactly against `core.md` (T5.2)/(T5.4) for split edges and against `zero_plus.hpp`'s `dS_e`. The
  NEW: exactly against `core.md` (T5.2)/(T5.4) for split edges and against `zero_plus.jl`'s `dS_e`. The
  ```
- **L817**
  ```
  OLD: `zero_plus.hpp`/K5/F30 report as the *observed* failure mode. Relatedly, the box's "decided by
  NEW: `zero_plus.jl`/K5/F30 report as the *observed* failure mode. Relatedly, the box's "decided by
  ```
- **L821**
  ```
  OLD: is what `contact.cpp` already does. So the engineering claim is safe and the theorem needs its
  NEW: is what `contact.jl` already does. So the engineering claim is safe and the theorem needs its
  ```

### `derivations/lemmas.md` — 25 edits, 0 leave

- **L4**
  ```
  OLD: `code/src` or `code/tests`. Notation is `core.md` §0 throughout (`c = cos(θ/2)`, `s = sin(θ/2)`,
  NEW: `Kirigami/src` or `Kirigami/test`. Notation is `core.md` §0 throughout (`c = cos(θ/2)`, `s = sin(θ/2)`,
  ```
- **L41**
  ```
  OLD:   `zero_plus.hpp`; `[F, K5 / F30]` reports non-opening split cuts as the observed failure mode, so
  NEW:   `zero_plus.jl`; `[F, K5 / F30]` reports non-opening split cuts as the observed failure mode, so
  ```
- **L55**
  ```
  OLD: face `g ≠ f`. Write `ρ : V′ → V` for the source map (`cut.hpp::prime_to_original`) and
  NEW: face `g ≠ f`. Write `ρ : V′ → V` for the source map (`cut.jl::prime_to_original`) and
  ```
- **L148**
  ```
  OLD: count of distinct pairs; `check_l1.cpp` increments once per sample. The Checker's independent run
  NEW: count of distinct pairs; `check_l1.jl` increments once per sample. The Checker's independent run
  ```
- **L253**
  ```
  OLD:    is the same computation as `zero_plus.hpp`'s header (`dC_e = 0`, `dS_e`) generalised from split
  NEW:    is the same computation as `zero_plus.jl`'s header (`dC_e = 0`, `dS_e`) generalised from split
  ```
- **L353**
  ```
  OLD: **H-L3 is not free.** Its first clause is exactly the predicate `zero_plus.hpp` is written to
  NEW: **H-L3 is not free.** Its first clause is exactly the predicate `zero_plus.jl` is written to
  ```
- **L399**
  ```
  OLD:    copies of the same `M`-vertex — the hinge point is welded by `cut.hpp::corner_to_prime` into a
  NEW:    copies of the same `M`-vertex — the hinge point is welded by `cut.jl::corner_to_prime` into a
  ```
- **L402**
  ```
  OLD:    `contact.cpp`'s `if (p == a || p == b) continue;` already does, *before* the numeric identity
  NEW:    `contact.jl`'s `if (p == a || p == b) continue;` already does, *before* the numeric identity
  ```
- **L449**
  ```
  OLD: know that the answer is structural.** That is what `contact.cpp` already does.
  NEW: know that the answer is structural.** That is what `contact.jl` already does.
  ```
- **L456–457**
  ```
  OLD: **Program** `derivations/scratch/check_l1.cpp`. **Build line** (in the file header; the repo's
       `clang++` needs `DEVELOPER_DIR` pointed at Xcode on this machine, see the note at the end of this
  NEW: **Program** `derivations/scratch/check_l1.jl`. **Run line** (in the file header):
  ```
- **L461–465**
  ```
  OLD: clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
         -I code/src -I code/apps derivations/scratch/check_l1.cpp \
         code/src/core/mesh.cpp code/src/core/cut.cpp code/src/core/holes.cpp \
         code/src/core/kinematics.cpp code/src/core/orientation.cpp code/src/core/generators.cpp \
         code/src/core/tutte_auxetic.cpp code/src/core/collision.cpp code/src/core/optimize.cpp \
  NEW: julia --project=Kirigami derivations/scratch/check_l1.jl
  ```
- **L469**
  ```
  OLD: **Corpus.** The 8 Phase-2 reference cases (`kill_common.hpp::reference_cases()`) plus
  NEW: **Corpus.** The 8 Phase-2 reference cases (`kill_common.jl::reference_cases (population frozen in data/corpus/)()`) plus
  ```
- **L545–547**
  ```
  OLD: `xcodebuild` fails to load, so a bare `clang++` aborts. Both scratch programs were built with
       `DEVELOPER_DIR=/Applications/Xcode.app/Contents/Developer` and the Xcode toolchain compiler
       `.../Toolchains/XcodeDefault.xctoolchain/usr/bin/clang++` with
  NEW: (Toolchain note obsolete after the Julia port: no compiler is involved; both scratch programs are run with `julia --project=Kirigami derivations/scratch/check_l{1,2}.jl`.)
  ```
- **L558**
  ```
  OLD: quotient system of `periodic_jacobian.hpp` supplies `X(t) = X₀ + Φ t` with `Φ ∈ R^{n_q × k}`,
  NEW: quotient system of `periodic_jacobian.jl` supplies `X(t) = X₀ + Φ t` with `Φ ∈ R^{n_q × k}`,
  ```
- **L561**
  ```
  OLD: orders `t` as `(t_{2i}, t_{2i+1}) = ` (`φ_i` on `x`, `φ_i` on `y`) — `kill_k7.cpp::shape_point`.
  NEW: orders `t` as `(t_{2i}, t_{2i+1}) = ` (`φ_i` on `x`, `φ_i` on `y`) — `kill_k7.jl::shape_point`.
  ```
- **L576**
  ```
  OLD: independent of `f`. ∎ This is the header derivation of `periodic_jacobian.hpp`, written out.
  NEW: independent of `f`. ∎ This is the header derivation of `periodic_jacobian.jl`, written out.
  ```
- **L775**
  ```
  OLD: *empirically*; and the corresponding one-line note in `periodic_jacobian.hpp`'s header, which states
  NEW: *empirically*; and the corresponding one-line note in `periodic_jacobian.jl`'s header, which states
  ```
- **L781**
  ```
  OLD: **Program** `derivations/scratch/check_l2.cpp`. Its population, quotient/super-patch pipeline and
  NEW: **Program** `derivations/scratch/check_l2.jl`. Its population, quotient/super-patch pipeline and
  ```
- **L782**
  ```
  OLD: `achievable()` are copied **verbatim** from `code/apps/kill_k7.cpp` (lines 36–243) so the numbers
  NEW: `achievable()` are copied **verbatim** from `Kirigami/apps/kill_k7.jl` (lines 36–243) so the numbers
  ```
- **L788–790**
  ```
  OLD: clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
         -I code/src -I code/apps derivations/scratch/check_l2.cpp \
         code/src/core/*.cpp code/src/method/*.cpp -o /tmp/check_l2 && /tmp/check_l2 100
  NEW: julia --project=Kirigami derivations/scratch/check_l2.jl 100
  ```
- **L880**
  ```
  OLD: **How to attack it.** I read `code/src/method/contact.cpp`'s `validity_certificate` and its
  NEW: **How to attack it.** I read `Kirigami/src/method/contact.jl`'s `validity_certificate` and its
  ```
- **L900**
  ```
  OLD: `kill_k7.cpp` — which reads it off the **even** generators only (`if (j % 2 == 0)`) — is not the
  NEW: `kill_k7.jl` — which reads it off the **even** generators only (`if (j % 2 == 0)`) — is not the
  ```
- **L904**
  ```
  OLD: `M_{2i+1}` against row `i` of `D`, which `kill_k7.cpp` never touches; Check L2 reports
  NEW: `M_{2i+1}` against row `i` of `D`, which `kill_k7.jl` never touches; Check L2 reports
  ```
- **L930–940**
  ```
  OLD: clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
         -I code/src -I code/apps derivations/scratch/check_l1.cpp \
         code/src/core/mesh.cpp code/src/core/cut.cpp code/src/core/holes.cpp \
         code/src/core/kinematics.cpp code/src/core/orientation.cpp code/src/core/generators.cpp \
         code/src/core/tutte_auxetic.cpp code/src/core/collision.cpp code/src/core/optimize.cpp \
         -o /tmp/check_l1 && /tmp/check_l1 140 20
       
       # Check L2  (133 periodic patterns, ~1.0 s)
       clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
         -I code/src -I code/apps derivations/scratch/check_l2.cpp \
         code/src/core/*.cpp code/src/method/*.cpp -o /tmp/check_l2 && /tmp/check_l2 100
  NEW: julia --project=Kirigami derivations/scratch/check_l1.jl
       julia --project=Kirigami derivations/scratch/check_l2.jl 100
  ```
- **L943–948**
  ```
  OLD: On this machine `xcode-select` points at an Xcode whose `xcodebuild` aborts, so a bare `clang++`
       fails before it reaches the source. Prefix both lines with
       `DEVELOPER_DIR=/Applications/Xcode.app/Contents/Developer` and use
       `/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang++`
       with `-isysroot .../Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk`. This is an environment
       defect, not a property of the code; `code/tests` were not run because they go through the same
  NEW: (The Xcode/`clang++` toolchain note that stood here is obsolete: the Julia port needs no compiler. The package test suite was run with `Pkg.test()`.)
  ```

### `derivations/lemmas_statements.md` — 4 edits, 0 leave

- **L39**
  ```
  OLD: face `g ≠ f`. Write `ρ : V′ → V` for the source map (`cut.hpp::prime_to_original`) and
  NEW: face `g ≠ f`. Write `ρ : V′ → V` for the source map (`cut.jl::prime_to_original`) and
  ```
- **L189**
  ```
  OLD: know that the answer is structural.** That is what `contact.cpp` already does.
  NEW: know that the answer is structural.** That is what `contact.jl` already does.
  ```
- **L195**
  ```
  OLD: quotient system of `periodic_jacobian.hpp` supplies `X(t) = X₀ + Φ t` with `Φ ∈ R^{n_q × k}`,
  NEW: quotient system of `periodic_jacobian.jl` supplies `X(t) = X₀ + Φ t` with `Φ ∈ R^{n_q × k}`,
  ```
- **L198**
  ```
  OLD: orders `t` as `(t_{2i}, t_{2i+1}) = ` (`φ_i` on `x`, `φ_i` on `y`) — `kill_k7.cpp::shape_point`.
  NEW: orders `t` as `(t_{2i}, t_{2i+1}) = ` (`φ_i` on `x`, `φ_i` on `y`) — `kill_k7.jl::shape_point`.
  ```

### `results/kill/KILL_REPORT.md` — 106 edits, 3 leave

- **L4**
  ```
  OLD: C++ drivers in `code/apps/kill_*.cpp` against `code/src/core` and `code/src/method`
  NEW: Julia drivers in `Kirigami/apps/kill_*.jl` against `Kirigami/src/core` and `Kirigami/src/method`
  ```
- **L7**
  ```
  OLD: (`code/apps/kill_common.hpp`), so every experiment sees the same graphs.
  NEW: (`data/corpus/ (frozen populations, see data/corpus/README.md)`), so every experiment sees the same graphs.
  ```
- **L18**
  ```
  OLD: `code/src/method/` and covered by a doctest in `code/tests/test_method.cpp`.
  NEW: `Kirigami/src/method/` and covered by a unit test in `Kirigami/test/test_method.jl`.
  ```
- **L89**
  ```
  OLD:    `0.9467`; a `1e−4` fine scan (`code/apps/dbg_k2a_out.cpp`) confirms the closed form. All
  NEW:    `0.9467`; a `1e−4` fine scan (`Kirigami/apps/dbg_k2a_out.jl`) confirms the closed form. All
  ```
- **L149**
  ```
  OLD: off by **exactly ±1**, and every one has `F ≥ 879` — above `mobility.hpp`'s
  NEW: off by **exactly ±1**, and every one has `F ≥ 879` — above `mobility.jl`'s
  ```
- **L151**
  ```
  OLD: `code/apps/kill_k3a_recheck.cpp` regenerates each violator with K3a's exact configuration
  NEW: `Kirigami/apps/kill_k3a_recheck.jl` regenerates each violator with K3a's exact configuration
  ```
- **L199**
  ```
  OLD: `recheck_summary.txt`. Driver: `code/apps/kill_k3a.cpp`, `code/apps/kill_k3a_recheck.cpp`.
  NEW: `recheck_summary.txt`. Driver: `Kirigami/apps/kill_k3a.jl`, `Kirigami/apps/kill_k3a_recheck.jl`.
  ```
- **L231**
  ```
  OLD: Artifacts: `results/kill/k1b/{k1b.csv,summary.txt}`. Driver: `code/apps/kill_k1b.cpp`.
  NEW: Artifacts: `results/kill/k1b/{k1b.csv,summary.txt}`. Driver: `Kirigami/apps/kill_k1b.jl`.
  ```
- **L274**
  ```
  OLD: (`deployable_population()` in `kill_common.hpp`), which is exactly the population K2a and
  NEW: (`deployable_population()` in `kill_common.jl (populations frozen in data/corpus/)`), which is exactly the population K2a and
  ```
- **L278**
  ```
  OLD: Artifacts: `results/kill/k1a/{k1a.csv,summary.txt}`. Driver: `code/apps/kill_k1a.cpp`.
  NEW: Artifacts: `results/kill/k1a/{k1a.csv,summary.txt}`. Driver: `Kirigami/apps/kill_k1a.jl`.
  ```
- **L344**
  ```
  OLD: `code/apps/kill_k2c.cpp`.
  NEW: `Kirigami/apps/kill_k2c.jl`.
  ```
- **L353**
  ```
  OLD: that root; (ii) `theta_max` from `collision.hpp` (grid scan + bisection, faces shrunk by
  NEW: that root; (ii) `theta_max` from `collision.jl` (grid scan + bisection, faces shrunk by
  ```
- **L399**
  ```
  OLD: overlap probes. Reproducing the Deriver's table exactly (`code/apps/dbg_t422.cpp`):
  NEW: overlap probes. Reproducing the Deriver's table exactly (`Kirigami/apps/dbg_t422.jl`):
  ```
- **L431**
  ```
  OLD: That was F32; with the interval tests restored (`contact.cpp`, §A3) the re-run gives
  NEW: That was F32; with the interval tests restored (`contact.jl`, §A3) the re-run gives
  ```
- **L440**
  ```
  OLD: Artifacts: `results/kill/k2a/{k2a.csv,summary.txt}`. Drivers: `code/apps/kill_k2a.cpp`,
  NEW: Artifacts: `results/kill/k2a/{k2a.csv,summary.txt}`. Drivers: `Kirigami/apps/kill_k2a.jl`,
  ```
- **L441**
  ```
  OLD: `code/apps/dbg_t422.cpp`, `code/apps/dbg_k2a_out.cpp`.
  NEW: `Kirigami/apps/dbg_t422.jl`, `Kirigami/apps/dbg_k2a_out.jl`.
  ```
- **L451**
  ```
  OLD: `collision.hpp`'s bisection to `<= 1e-6` rad. PASS if the fraction of gamma-ladder-certified
  NEW: `collision.jl`'s bisection to `<= 1e-6` rad. PASS if the fraction of gamma-ladder-certified
  ```
- **L520**
  ```
  OLD: Artifacts: `results/kill/k1c/{k1c.csv,summary.txt}`. Driver: `code/apps/kill_k1c.cpp`.
  NEW: Artifacts: `results/kill/k1c/{k1c.csv,summary.txt}`. Driver: `Kirigami/apps/kill_k1c.jl`.
  ```
- **L584**
  ```
  OLD: after the fix, `1e−12` on 36 and `1e−10` on the other 89 — `results/core_validation/referee_fix.md`), and a direct scan (`code/apps/dbg_k5.cpp`) confirms overlap at
  NEW: after the fix, `1e−12` on 36 and `1e−10` on the other 89 — `results/core_validation/referee_fix.md`), and a direct scan (`Kirigami/apps/dbg_k5.jl`) confirms overlap at
  ```
- **L616**
  ```
  OLD: `code/apps/kill_k5.cpp`, diagnostic `code/apps/dbg_k5.cpp`.
  NEW: `Kirigami/apps/kill_k5.jl`, diagnostic `Kirigami/apps/dbg_k5.jl`.
  ```
- **L625**
  ```
  OLD: (their own `theta_max` routine is buggy, `STATE.md` F24). Keep our `collision.hpp` Eq.(9) as
  NEW: (their own `theta_max` routine is buggy, `STATE.md` F24). Keep our `collision.jl` Eq.(9) as
  ```
- **L637** *leave* — third-party `baseline/` reference, KEEP: `dumped (a `--dumpdir` option added to our CLI in `baseline/native/src/main.cpp`; nothing in`
- **L710** *leave* — bare driver name, unchanged: ``code/apps/kill_k2b.cpp`; CLI change in `baseline/native/src/main.cpp` (`--dumpdir`,`
- **L726**
  ```
  OLD: **Measured** (`code/apps/kill_f23.cpp`), stronger than asked: **every** seed face, not two,
  NEW: **Measured** (`Kirigami/apps/kill_f23.jl`), stronger than asked: **every** seed face, not two,
  ```
- **L759**
  ```
  OLD: Artifacts: `results/kill/f23/{f23.csv,summary.txt}`. Driver: `code/apps/kill_f23.cpp`.
  NEW: Artifacts: `results/kill/f23/{f23.csv,summary.txt}`. Driver: `Kirigami/apps/kill_f23.jl`.
  ```
- **L839**
  ```
  OLD: quoted** until `mobility.hpp`'s `matrix_rank` is fixed or the dense path is forced; a dense
  NEW: quoted** until `mobility.jl`'s `matrix_rank` is fixed or the dense path is forced; a dense
  ```
- **L841**
  ```
  OLD: `code/apps/kill_k3a_recheck.cpp` provides (7 of 7 violators small enough for a dense QR have
  NEW: `Kirigami/apps/kill_k3a_recheck.jl` provides (⟨JULIA:kill_k3a_recheck:7 of 7⟩ violators small enough for a dense QR have
  ```
- **L858**
  ```
  OLD: **The 0⁺ calculus** (`code/src/method/zero_plus.{hpp,cpp}`, derived there in full). By T1.B
  NEW: **The 0⁺ calculus** (`Kirigami/src/method/zero_plus.jl`, derived there in full). By T1.B
  ```
- **L952**
  ```
  OLD: The fix is the three-line guard now at `code/src/method/contact.cpp:339` — keep a root only
  NEW: The fix is the three-line guard now at `Kirigami/src/method/contact.jl:339` — keep a root only
  ```
- **L978**
  ```
  OLD: native/, smoke/}`. Driver `code/apps/kill_k6.cpp`, method `code/src/method/zero_plus.{hpp,cpp}`,
  NEW: native/, smoke/}`. Driver `Kirigami/apps/kill_k6.jl`, method `Kirigami/src/method/zero_plus.jl`,
  ```
- **L979**
  ```
  OLD: diagnostic `code/apps/dbg_k6.cpp`.
  NEW: diagnostic `Kirigami/apps/dbg_k6.jl`.
  ```
- **L986–987**
  ```
  OLD: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
       ./code/build/kiri_tests                    # 53 cases, 16004 assertions, all passing
  NEW: julia --project=Kirigami -e 'using Pkg; Pkg.test()'   # ⟨JULIA:tests⟩, all passing
  ```
- **L988**
  ```
  OLD: ./code/build/kill_k3a  --out results/kill/k3a          # ~25 min
  NEW: julia --project=Kirigami Kirigami/apps/kill_k3a.jl  --out results/kill/k3a          # ~25 min
  ```
- **L989**
  ```
  OLD: ./code/build/kill_k3a_recheck --out results/kill/k3a   # ~4 min
  NEW: julia --project=Kirigami Kirigami/apps/kill_k3a_recheck.jl --out results/kill/k3a   # ~4 min
  ```
- **L990**
  ```
  OLD: ./code/build/kill_k1b  --out results/kill/k1b          # ~34 s
  NEW: julia --project=Kirigami Kirigami/apps/kill_k1b.jl  --out results/kill/k1b          # ~34 s
  ```
- **L991**
  ```
  OLD: ./code/build/kill_k1a  --out results/kill/k1a          # ~40 s
  NEW: julia --project=Kirigami Kirigami/apps/kill_k1a.jl  --out results/kill/k1a          # ~40 s
  ```
- **L992**
  ```
  OLD: ./code/build/kill_k5   --out results/kill/k5 --n 200 --maxf 800   # ~5 min
  NEW: julia --project=Kirigami Kirigami/apps/kill_k5.jl   --out results/kill/k5 --n 200 --maxf 800   # ~5 min
  ```
- **L993**
  ```
  OLD: ./code/build/kill_k2c  --out results/kill/k2c          # ~3 min
  NEW: julia --project=Kirigami Kirigami/apps/kill_k2c.jl  --out results/kill/k2c          # ~3 min
  ```
- **L994**
  ```
  OLD: ./code/build/kill_k2a  --out results/kill/k2a          # ~7 s
  NEW: julia --project=Kirigami Kirigami/apps/kill_k2a.jl  --out results/kill/k2a          # ~7 s
  ```
- **L995**
  ```
  OLD: ./code/build/kill_k1c  --out results/kill/k1c          # ~31 s   (needs baseline/native)
  NEW: julia --project=Kirigami Kirigami/apps/kill_k1c.jl  --out results/kill/k1c          # ~31 s   (needs baseline/native)
  ```
- **L996**
  ```
  OLD: ./code/build/kill_k2b  --out results/kill/k2b --maxf 160   # ~26 min (needs baseline/native)
  NEW: julia --project=Kirigami Kirigami/apps/kill_k2b.jl  --out results/kill/k2b --maxf 160   # ~26 min (needs baseline/native)
  ```
- **L997**
  ```
  OLD: ./code/build/kill_f23  --out results/kill/f23          # ~1 s
  NEW: julia --project=Kirigami Kirigami/apps/kill_f23.jl  --out results/kill/f23          # ~1 s
  ```
- **L1004**
  ```
  OLD: Method code changed for the corrections: `code/src/method/deploy_basis.{hpp,cpp}`
  NEW: Method code changed for the corrections: `Kirigami/src/method/deploy_basis.jl`
  ```
- **L1006**
  ```
  OLD: `code/src/method/contact.{hpp,cpp}` (`exact_theta_max_overlap`, `validity_certificate`, the
  NEW: `Kirigami/src/method/contact.jl` (`exact_theta_max_overlap`, `validity_certificate`, the
  ```
- **L1007**
  ```
  OLD: exact face-frame radius in the broad phase). New doctest coverage in
  NEW: exact face-frame radius in the broad phase). New unit test coverage in
  ```
- **L1008**
  ```
  OLD: `code/tests/test_method.cpp`: the hexagon graze; `Θ_max` vs bisection and pruning
  NEW: `Kirigami/test/test_method.jl`: the hexagon graze; `Θ_max` vs bisection and pruning
  ```
- **L1036**
  ```
  OLD: * **Driver** `code/apps/kill_jitter.cpp`. No `core/` or `method/` source was modified.
  NEW: * **Driver** `Kirigami/apps/kill_jitter.jl`. No `core/` or `method/` source was modified.
  ```
- **L1053**
  ```
  OLD:   (`--analyze`) is also C++; python draws the figure only.
  NEW:   (`--analyze`) is also Julia; python draws the figure only.
  ```
- **L1069**
  ```
  OLD: tolerate jitter but because **they never left the start point**. New doctest in
  NEW: tolerate jitter but because **they never left the start point**. New unit test in
  ```
- **L1070**
  ```
  OLD: `code/tests/test_method.cpp` ("split-free tiling: Eq. (6) projection undoes an arbitrary
  NEW: `Kirigami/test/test_method.jl` ("split-free tiling: Eq. (6) projection undoes an arbitrary
  ```
- **L1164–1165**
  ```
  OLD: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
       ./code/build/kiri_tests                              # 65 cases, 16229 assertions, all passing
  NEW: julia --project=Kirigami -e 'using Pkg; Pkg.test()'   # ⟨JULIA:tests⟩, all passing
  ```
- **L1167**
  ```
  OLD: seq 0 11 | xargs -P 12 -I{} ./code/build/kill_jitter --shard {} --nshards 12 \
  NEW: seq 0 11 | xargs -P 12 -I{} julia --project=Kirigami Kirigami/apps/kill_jitter.jl --shard {} --nshards 12 \
  ```
- **L1171**
  ```
  OLD: ./code/build/kill_jitter --analyze results/kill/jitter/jitter.csv --out results/kill/jitter
  NEW: julia --project=Kirigami Kirigami/apps/kill_jitter.jl --analyze results/kill/jitter/jitter.csv --out results/kill/jitter
  ```
- **L1209**
  ```
  OLD: * **Driver** `code/apps/kill_k7.cpp`, stages `main`, `c3`, `c4bounded`, `nu`, `c4sweep`,
  NEW: * **Driver** `Kirigami/apps/kill_k7.jl`, stages `main`, `c3`, `c4bounded`, `nu`, `c4sweep`,
  ```
- **L1211**
  ```
  OLD:   patch, `periodic_jacobian`, the achievable set — is `code/src/method/periodic_jacobian.{hpp,cpp}`,
  NEW:   patch, `periodic_jacobian`, the achievable set — is `Kirigami/src/method/periodic_jacobian.jl`,
  ```
- **L1410**
  ```
  OLD: **New doctest.** `code/tests/test_method.cpp`, "C4: theta_c = 2 atan2(tr K, 1 − det K) is where
  NEW: **New unit test.** `Kirigami/test/test_method.jl`, "C4: theta_c = 2 atan2(tr K, 1 − det K) is where
  ```
- **L1422–1424**
  ```
  OLD: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
       ./code/build/kiri_tests                              # 66 cases, 16283 assertions, all passing
       ./code/build/kiri_tests -tc="C4: theta_c*"           # the K7 doctest alone: 24 assertions
  NEW: julia --project=Kirigami -e 'using Pkg; Pkg.test()'   # ⟨JULIA:tests⟩, all passing
       julia --project=Kirigami -e 'using Pkg; Pkg.test(test_args=["C4: theta_c"])'   # the K7 test alone: ⟨JULIA:tests:C4⟩
  ```
- **L1428**
  ```
  OLD: for i in $(seq 0 11); do ./code/build/kill_k7 --stage main --out results/kill/k7 \
  NEW: for i in $(seq 0 11); do julia --project=Kirigami Kirigami/apps/kill_k7.jl --stage main --out results/kill/k7 \
  ```
- **L1432**
  ```
  OLD: for i in $(seq 0 11); do ./code/build/kill_k7 --stage c3 --out results/kill/k7 \
  NEW: for i in $(seq 0 11); do julia --project=Kirigami Kirigami/apps/kill_k7.jl --stage c3 --out results/kill/k7 \
  ```
- **L1435**
  ```
  OLD: for i in $(seq 0 11); do ./code/build/kill_k7 --stage c3 --out results/kill/k7/c3_v2 \
  NEW: for i in $(seq 0 11); do julia --project=Kirigami Kirigami/apps/kill_k7.jl --stage c3 --out results/kill/k7/c3_v2 \
  ```
- **L1437**
  ```
  OLD: ./code/build/kill_k7 --stage c4bounded --out results/kill/k7
  NEW: julia --project=Kirigami Kirigami/apps/kill_k7.jl --stage c4bounded --out results/kill/k7
  ```
- **L1438**
  ```
  OLD: ./code/build/kill_k7 --stage nu        --out results/kill/k7
  NEW: julia --project=Kirigami Kirigami/apps/kill_k7.jl --stage nu        --out results/kill/k7
  ```
- **L1439**
  ```
  OLD: ./code/build/kill_k7 --stage c4sweep   --out results/kill/k7   # the |det P| artefact
  NEW: julia --project=Kirigami Kirigami/apps/kill_k7.jl --stage c4sweep   --out results/kill/k7   # the |det P| artefact
  ```
- **L1440**
  ```
  OLD: for i in $(seq 0 11); do ./code/build/kill_k7 --stage perdbg --out results/kill/k7 \
  NEW: for i in $(seq 0 11); do julia --project=Kirigami Kirigami/apps/kill_k7.jl --stage perdbg --out results/kill/k7 \
  ```
- **L1478**
  ```
  OLD: `kill_common::make_graph`, the same `σ_mc` from Eq. (1) and `σ_def` read back from
  NEW: `kill_common.jl::make_graph` (frozen in `data/corpus/`), the same `σ_mc` from Eq. (1) and `σ_def` read back from
  ```
- **L1581**
  ```
  OLD: No core or method source was modified for this run; `code/apps/kill_b4.cpp` builds its
  NEW: No core or method source was modified for this run; `Kirigami/apps/kill_b4.jl` builds its
  ```
- **L1586**
  ```
  OLD: `code/apps/kill_b4.cpp`.
  NEW: `Kirigami/apps/kill_b4.jl`.
  ```
- **L1589**
  ```
  OLD: ./code/build/kill_b4 --n 200 --nshards 12 --shard I --out results/kill/b4/shards \
  NEW: julia --project=Kirigami Kirigami/apps/kill_b4.jl --n 200 --nshards 12 --shard I --out results/kill/b4/shards \
  ```
- **L1599**
  ```
  OLD: Driver `code/apps/kill_k8a.cpp`, method `code/src/method/expansive_cone.{hpp,cpp}`,
  NEW: Driver `Kirigami/apps/kill_k8a.jl`, method `Kirigami/src/method/expansive_cone.jl`,
  ```
- **L1600**
  ```
  OLD: artifacts `results/kill/k8a/` (`k8a.csv`, `summary.txt`, `log_*.txt`), 5 new doctests in
  NEW: artifacts `results/kill/k8a/` (`k8a.csv`, `summary.txt`, `log_*.txt`), 5 new unit tests in
  ```
- **L1601**
  ```
  OLD: `code/tests/test_method.cpp` (suite now 71 cases / 16 488 assertions, all passing).
  NEW: `Kirigami/test/test_method.jl` (suite now ⟨JULIA:tests⟩).
  ```
- **L1607**
  ```
  OLD: which is exactly `mobility.hpp::build_rigidity`. Because every copy of a source vertex
  NEW: which is exactly `mobility.jl::build_rigidity`. Because every copy of a source vertex
  ```
- **L1609**
  ```
  OLD: `zero_plus.hpp` becomes a **strict linear inequality in the flex**:
  NEW: `zero_plus.jl` becomes a **strict linear inequality in the flex**:
  ```
- **L1618**
  ```
  OLD: `1e-9` relative) against `zero_plus_q` and `zero_plus_corner_margin` in the doctest *"the
  NEW: `1e-9` relative) against `zero_plus_q` and `zero_plus_corner_margin` in the unit test *"the
  ```
- **L1621**
  ```
  OLD: translations plus one free translation pair per component of `Γ`; the doctest checks
  NEW: translations plus one free translation pair per component of `Γ`; the unit test checks
  ```
- **L1631**
  ```
  OLD: The header of `expansive_cone.hpp` carries this citation.
  NEW: The header of `expansive_cone.jl` carries this citation.
  ```
- **L1649**
  ```
  OLD: LPs of the doctest *"cone_lp solves three hand-solved linear programs"* pin the solver:
  NEW: LPs of the unit test *"cone_lp solves three hand-solved linear programs"* pin the solver:
  ```
- **L1683**
  ```
  OLD: Extended control, `kill_common::deployable_population` (authored tilings at five clip
  NEW: Extended control, `kill_common.jl::deployable_population` (frozen in `data/corpus/`) (authored tilings at five clip
  ```
- **L1740**
  ```
  OLD: *Deviation, stated:* `contact.hpp`'s T4.2″ scan is parameterised by the **uniform** angle
  NEW: *Deviation, stated:* `contact.jl`'s T4.2″ scan is parameterised by the **uniform** angle
  ```
- **L1743**
  ```
  OLD: itself calls at every interval midpoint — and a doctest checks that along the uniform ray
  NEW: itself calls at every interval midpoint — and a unit test checks that along the uniform ray
  ```
- **L1749–1750**
  ```
  OLD: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
       ./code/build/kiri_tests
  NEW: julia --project=Kirigami -e 'using Pkg; Pkg.test()'
  ```
- **L1751**
  ```
  OLD: for i in $(seq 0 7); do ./code/build/kill_k8a --n 100 --x0 --dual-iters 100000 \
  NEW: for i in $(seq 0 7); do julia --project=Kirigami Kirigami/apps/kill_k8a.jl --n 100 --x0 --dual-iters 100000 \
  ```
- **L1760**
  ```
  OLD: Every graph is a deterministic function of its id (`kill_common::make_graph(id, 100, 800,
  NEW: Every graph is a deterministic function of its id (`kill_common.jl::make_graph(id, 100, 800,
  ```
- **L1787**
  ```
  OLD: `code/apps/kill_k8a_recheck.cpp`, 4 new doctests. **Nothing above is deleted; every number
  NEW: `Kirigami/apps/kill_k8a_recheck.jl`, 4 new unit tests. **Nothing above is deleted; every number
  ```
- **L1857**
  ```
  OLD: **Test suite.** 4 new doctests in `code/tests/test_method.cpp`: the bracket closes to `1e−9`
  NEW: **Test suite.** 4 new unit tests in `Kirigami/test/test_method.jl`: the bracket closes to `1e−9`
  ```
- **L1866**
  ```
  OLD: `code/apps/kill_b3.cpp`, library `code/src/method/budget.{hpp,cpp}` (5 new doctests in
  NEW: `Kirigami/apps/kill_b3.jl`, library `Kirigami/src/method/budget.jl` (5 new unit tests in
  ```
- **L1867**
  ```
  OLD: `code/tests/test_method.cpp`; suite 80 cases / 16,785 assertions pass). Outputs
  NEW: `Kirigami/test/test_method.jl`; suite ⟨JULIA:tests⟩ pass). Outputs
  ```
- **L1873**
  ```
  OLD: **B1 global half (the part `derivations/scratch/check_b1.cpp` could not close).** On the
  NEW: **B1 global half (the part `derivations/scratch/check_b1.jl` could not close).** On the
  ```
- **L1948**
  ```
  OLD: **The problem solved** (`code/src/method/convex_embed.{hpp,cpp}`, derived there in full).
  NEW: **The problem solved** (`Kirigami/src/method/convex_embed.jl`, derived there in full).
  ```
- **L2021**
  ```
  OLD: `Θ_bisect = 0` at the historical shrink `1e−12`. `code/apps/dbg_k9.cpp` traces every one to
  NEW: `Θ_bisect = 0` at the historical shrink `1e−12`. `Kirigami/apps/dbg_k9.jl` traces every one to
  ```
- **L2051**
  ```
  OLD: `min μ > 0`, and `expansive_cone.hpp`'s rows equal `q/2` and `μ/2` at the `σ` flex, asserted
  NEW: `min μ > 0`, and `expansive_cone.jl`'s rows equal `q/2` and `μ/2` at the `σ` flex, asserted
  ```
- **L2052**
  ```
  OLD: by doctest). A strictly feasible point therefore **provably exists** on those 23 and the LP
  NEW: by unit test). A strictly feasible point therefore **provably exists** on those 23 and the LP
  ```
- **L2078**
  ```
  OLD: `code/apps/kill_k9.cpp`, method `code/src/method/convex_embed.{hpp,cpp}`, referee diagnostic
  NEW: `Kirigami/apps/kill_k9.jl`, method `Kirigami/src/method/convex_embed.jl`, referee diagnostic
  ```
- **L2079**
  ```
  OLD: `code/apps/dbg_k9.cpp`, doctests in `code/tests/test_method.cpp` (hand-checked
  NEW: `Kirigami/apps/dbg_k9.jl`, unit tests in `Kirigami/test/test_method.jl` (hand-checked
  ```
- **L2093**
  ```
  OLD: **The four levers** (`code/apps/kill_k9b.cpp`), in the order the plan ranked them:
  NEW: **The four levers** (`Kirigami/apps/kill_k9b.jl`), in the order the plan ranked them:
  ```
- **L2103**
  ```
  OLD: * **(iv) stage 2**, `range_opt.hpp`'s softmin-of-first-contact objective (the T6 analytic
  NEW: * **(iv) stage 2**, `range_opt.jl`'s softmin-of-first-contact objective (the T6 analytic
  ```
- **L2183**
  ```
  OLD: logs/, k9b_gallery.png, k9b_hist.png}`. Driver `code/apps/kill_k9b.cpp`, plots
  NEW: logs/, k9b_gallery.png, k9b_hist.png}`. Driver `Kirigami/apps/kill_k9b.jl`, plots
  ```
- **L2206**
  ```
  OLD: `code/src/method/range_embed.{hpp,cpp}` derives it in full. The scalar is the **0⁺
  NEW: `Kirigami/src/method/range_embed.jl` derives it in full. The scalar is the **0⁺
  ```
- **L2213**
  ```
  OLD: — the joint minimum of `zero_plus.hpp`'s two first-order separation families, the split-edge
  NEW: — the joint minimum of `zero_plus.jl`'s two first-order separation families, the split-edge
  ```
- **L2216**
  ```
  OLD: sign of `cross = det(e1, e2)` — the *same* cross that `convex_embed.hpp` constrains, since
  NEW: sign of `cross = det(e1, e2)` — the *same* cross that `convex_embed.jl` constrains, since
  ```
- **L2230**
  ```
  OLD: is always the exact minimum, never the surrogate. Two doctests check exactly those two
  NEW: is always the exact minimum, never the surrogate. Two unit tests check exactly those two
  ```
- **L2231**
  ```
  OLD: claims (`code/tests/test_range_embed.cpp`), and two more finite-difference the analytic
  NEW: claims (`Kirigami/test/test_range_embed.jl`), and two more finite-difference the analytic
  ```
- **L2239**
  ```
  OLD: `Θ_max` with `range_opt.hpp`'s softmin-of-first-contact objective (the T6 analytic
  NEW: `Θ_max` with `range_opt.jl`'s softmin-of-first-contact objective (the T6 analytic
  ```
- **L2269** *leave* — bare driver name, unchanged: ``kill_k9c --aggregate --nshards 12` and `results/kill/k9c/{k9c.csv, summary.txt}` now hold the`
- **L2423**
  ```
  OLD: k9c_gallery.png, k9c_hist.png}`. Driver `code/apps/kill_k9c.cpp`, method
  NEW: k9c_gallery.png, k9c_hist.png}`. Driver `Kirigami/apps/kill_k9c.jl`, method
  ```
- **L2424**
  ```
  OLD: `code/src/method/range_embed.{hpp,cpp}`, doctests `code/tests/test_range_embed.cpp`, plots
  NEW: `Kirigami/src/method/range_embed.jl`, unit tests `Kirigami/test/test_range_embed.jl`, plots
  ```
- **L2438**
  ```
  OLD: **How run.** `code/apps/kill_native200.cpp`, three variants per graph: *native* (their own
  NEW: **How run.** `Kirigami/apps/kill_native200.jl`, three variants per graph: *native* (their own
  ```
- **L2515**
  ```
  OLD: run_shard_{}.log}`. Driver `code/apps/kill_native200.cpp`. K9's numbers for the same 200
  NEW: run_shard_{}.log}`. Driver `Kirigami/apps/kill_native200.jl`. K9's numbers for the same 200
  ```
- **L2519**
  ```
  OLD: ./code/build/kill_native200 --n 200 --maxf 800 --timeout 600 --shard {} --nshards 12 \
  NEW: julia --project=Kirigami Kirigami/apps/kill_native200.jl --n 200 --maxf 800 --timeout 600 --shard {} --nshards 12 \
  ```

### `results/kill/jitter/cert_diagnosis.md` — 11 edits, 4 leave

- **L4**
  ```
  OLD: C++ code in `code/`; no number is quoted from memory.
  NEW: Julia code in `Kirigami/`; no number is quoted from memory.
  ```
- **L38**
  ```
  OLD: `contact.cpp::contact_angles()` (the exact scan) applies exactly that test, with
  NEW: `contact.jl::contact_angles()` (the exact scan) applies exactly that test, with
  ```
- **L40**
  ```
  OLD: `contact.cpp::validity_certificate()`'s NOROOT scan did not: it called
  NEW: `contact.jl::validity_certificate()`'s NOROOT scan did not: it called
  ```
- **L52**
  ```
  OLD: `code/apps/dbg_cert.cpp` (new) re-runs one jitter row through the identical `measure()`
  NEW: `Kirigami/apps/dbg_cert.jl` (new) re-runs one jitter row through the identical `measure()`
  ```
- **L95**
  ```
  OLD: `code/src/method/contact.cpp`, in `validity_certificate`'s NOROOT scan: build the edge's
  NEW: `Kirigami/src/method/contact.jl`, in `validity_certificate`'s NOROOT scan: build the edge's
  ```
- **L150**
  ```
  OLD: `derivation_tests`' R5-a…R5-d block (which checks (T5.2b''-1b) on the corpus) passes
  NEW: `derivation_tests.jl`' R5-a…R5-d block (which checks (T5.2b''-1b) on the corpus) passes
  ```
- **L155**
  ```
  OLD: **Regression test** (`code/tests/test_method.cpp`, "F32: NOROOT counts only roots with the
  NEW: **Regression test** (`Kirigami/test/test_method.jl`, "F32: NOROOT counts only roots with the
  ```
- **L161**
  ```
  OLD: * Without the fix (`git stash` of `contact.cpp` only, everything else identical):
  NEW: * Without the fix (`git stash` of `contact.jl` only, everything else identical):
  ```
- **L175**
  ```
  OLD: **Suites** (`cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j`):
  NEW: **Suites** (`julia --project=Kirigami -e 'using Pkg; Pkg.test()'`):
  ```
- **L179–180**
  ```
  OLD: | `kiri_tests` | 66/66 test cases, **16283/16283 assertions**, SUCCESS |
       | `derivation_tests` | 31/31 test cases, **89055/89055 assertions**, SUCCESS |
  NEW: | `Pkg.test()` (whole suite) | ⟨JULIA:tests⟩, SUCCESS |
       | `derivation_tests.jl` | ⟨JULIA:derivation_tests⟩, SUCCESS |
  ```
- **L182** *leave* — bare driver name, unchanged: `## 6. Re-measurement of A3 (`kill_jitter`, 12 shards, ~2 s)`
- **L212** *leave* — bare driver name, unchanged: `**K2a — the numbers change, the claims survive and get stronger.** Re-ran `kill_k2a``
- **L241** *leave* — bare driver name, unchanged: `**K5** (`kill_k5`, full 200 graphs, both sigma rules, `eps = 0.3`) — identical to the`
- **L259** *leave* — bare driver name, unchanged: `**K6** (`kill_k6 --n 20 --no-native`, a smoke re-run under the same code): every clause`
- **L278**
  ```
  OLD: `./code/build/kill_k6 --n 20 --no-native` (40 rows = 20 graphs x 2 sigma):
  NEW: `julia --project=Kirigami Kirigami/apps/kill_k6.jl --n 20 --no-native` (40 rows = 20 graphs x 2 sigma):
  ```

### `results/kill/k8a/recheck.md` — 7 edits, 0 leave

- **L22**
  ```
  OLD:   (`expansive_cone.cpp` lines 208–243 of the original), started from `z = 0` with step
  NEW:   (`expansive_cone.jl` lines 208–243 of the original), started from `z = 0` with step
  ```
- **L44**
  ```
  OLD: `code/src/method/expansive_cone.cpp`, `cone_lp`:
  NEW: `Kirigami/src/method/expansive_cone.jl`, `cone_lp`:
  ```
- **L74**
  ```
  OLD: `code/apps/kill_k8a_recheck.cpp` regenerates K9's variant-(b) embeddings bit-identically
  NEW: `Kirigami/apps/kill_k8a_recheck.jl` regenerates K9's variant-(b) embeddings bit-identically
  ```
- **L139**
  ```
  OLD: `code/apps/kill_k8a.cpp`, 100 graphs, both embeddings, 8 shards, corrected solver:
  NEW: `Kirigami/apps/kill_k8a.jl`, 100 graphs, both embeddings, 8 shards, corrected solver:
  ```
- **L221–222**
  ```
  OLD: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
       ./code/build/kiri_tests                      # 4 new K8a-recheck cases
  NEW: julia --project=Kirigami -e 'using Pkg; Pkg.test()'   # 4 K8a-recheck cases included
  ```
- **L225**
  ```
  OLD: for i in $(seq 0 11); do ./code/build/kill_k8a_recheck --out results/kill/k8a \
  NEW: for i in $(seq 0 11); do julia --project=Kirigami Kirigami/apps/kill_k8a_recheck.jl --out results/kill/k8a \
  ```
- **L229**
  ```
  OLD: for i in $(seq 0 7); do ./code/build/kill_k8a --n 100 --x0 --dual-iters 100000 \
  NEW: for i in $(seq 0 7); do julia --project=Kirigami Kirigami/apps/kill_k8a.jl --n 100 --x0 --dual-iters 100000 \
  ```

### `results/kill/native200/NATIVE200_FINAL.md` — 1 edits, 0 leave

- **L5**
  ```
  OLD: Every number below is computed by `code/apps/native200_merge.cpp` from `results/kill/native200/native200_final.csv`, which is itself merged from `native200.csv` (the 600 s run), `rerun3600/shard_*.csv` (the 3600 s rerun of the timed-out and never-dispatched cells) and `crashfix3600/shard_*.csv` (the 39 crashed cells against the F24-patched CLI). One row per (graph id, variant); the 600 s status of a superseded cell is kept in `first_status`.
  NEW: Every number below is computed by `Kirigami/apps/native200_merge.jl` from `results/kill/native200/native200_final.csv`, which is itself merged from `native200.csv` (the 600 s run), `rerun3600/shard_*.csv` (the 3600 s rerun of the timed-out and never-dispatched cells) and `crashfix3600/shard_*.csv` (the 39 crashed cells against the F24-patched CLI). One row per (graph id, variant); the 600 s status of a superseded cell is kept in `first_status`.
  ```

### `results/core_validation/rank_claim.md` — 4 edits, 0 leave

- **L49**
  ```
  OLD: instances by `kiri_tests` (see `tests/test_holes.cpp`).
  NEW: instances by `Pkg.test()` (see `test/test_holes.jl`).
  ```
- **L120**
  ```
  OLD: patches: the three torus rows above and `tests/test_rank_checks.cpp`, where `dim Z == 1`
  NEW: patches: the three torus rows above and `test/test_rank_checks.jl`, where `dim Z == 1`
  ```
- **L140**
  ```
  OLD: notches -- the torus patches in `tests/test_rank_checks.cpp` -- it collapses to the
  NEW: notches -- the torus patches in `test/test_rank_checks.jl` -- it collapses to the
  ```
- **L157**
  ```
  OLD: `tests/test_rank_checks.cpp` on the 4x4 and 6x4 square tori and the 4x4 triangle torus.
  NEW: `test/test_rank_checks.jl` on the 4x4 and 6x4 square tori and the 4x4 triangle torus.
  ```

### `results/core_validation/referee_fix.md` — 10 edits, 0 leave

- **L3**
  ```
  OLD: Scope: `code/src/core/collision.{hpp,cpp}` -- `polygons_overlap` / `has_collision`, the
  NEW: Scope: `Kirigami/src/core/collision.jl` -- `polygons_overlap` / `has_collision`, the
  ```
- **L74**
  ```
  OLD: New, in `code/tests/test_collision.cpp` (written first, failing against the old predicate):
  NEW: New, in `Kirigami/test/test_collision.jl` (written first, failing against the old predicate):
  ```
- **L89**
  ```
  OLD: `derivation_tests.cpp` R6-c, the pinned reproduction from checker round 6, is INVERTED:
  NEW: `derivation_tests.jl` R6-c, the pinned reproduction from checker round 6, is INVERTED:
  ```
- **L93**
  ```
  OLD: Suites (both re-run after the fix, from the repo root; the `kiri_tests` totals include test
  NEW: Suites (both re-run after the fix, from the repo root; the `Pkg.test()` totals include test
  ```
- **L98–99**
  ```
  OLD: | `kiri_tests` | 103 | 17 504 | 0 failures |
       | `derivation_tests` | 33 | 89 074 | 0 failures |
  NEW: | `Pkg.test()` (whole suite) | ⟨JULIA:tests⟩ | 0 failures |
       | `derivation_tests.jl` | ⟨JULIA:derivation_tests⟩ | 0 failures |
  ```
- **L104**
  ```
  OLD: tree with only `code/src/core/collision.{hpp,cpp}` reverted to `HEAD` (built out of tree in
  NEW: tree with only `Kirigami/src/core/collision.jl` reverted to `HEAD` (built out of tree in
  ```
- **L106**
  ```
  OLD: agents made to `contact.cpp` / `zero_plus.cpp` in the same working tree. The baseline K5 run
  NEW: agents made to `contact.jl` / `zero_plus.jl` in the same working tree. The baseline K5 run
  ```
- **L128**
  ```
  OLD: ### B4 -- the F34 case itself (`derivations/scratch/check_b4_93.cpp`, ids 93 and 96)
  NEW: ### B4 -- the F34 case itself (`derivations/scratch/check_b4_93.jl`, ids 93 and 96)
  ```
- **L171**
  ```
  OLD: probe, `contact.cpp`), not only in the referee: the same hinge artefact was rejecting the
  NEW: probe, `contact.jl`), not only in the referee: the same hinge artefact was rejecting the
  ```
- **L178**
  ```
  OLD: `method/zero_plus.cpp` in the same working tree, not this fix.
  NEW: `method/zero_plus.jl` in the same working tree, not this fix.
  ```

### `results/final/e1/E1.md` — 4 edits, 0 leave

- **L11**
  ```
  OLD: Code: `code/apps/kill_e1.cpp` (new; does not modify any core/ or method/ source). Three
  NEW: Code: `Kirigami/apps/kill_e1.jl` (new; does not modify any core/ or method/ source). Three
  ```
- **L16**
  ```
  OLD: - **authored**: `kill_common::deployable_population()` called with (samples_per_base=20,
  NEW: - **authored**: `kill_common.jl::deployable_population()` (frozen in `data/corpus/`) called with (samples_per_base=20,
  ```
- **L21**
  ```
  OLD: - **random**: `kill_common::make_graph(id, 100, 800, 1400)` -- Voronoi / Delaunay /
  NEW: - **random**: `kill_common.jl::make_graph(id, 100, 800, 1400)` (frozen in `data/corpus/`) -- Voronoi / Delaunay /
  ```
- **L101**
  ```
  OLD: cmake --build . --target kill_e1
  NEW: # no build step; Kirigami/apps/kill_e1.jl runs from source (via Kirigami/scripts/run_e1.sh)
  ```

### `results/final/figures/README.md` — 1 edits, 0 leave

- **L5**
  ```
  OLD: re-run any C++ solver).
  NEW: re-run any Julia solver).
  ```

### `results/regime/REGIME.md` — 1 edits, 0 leave

- **L6**
  ```
  OLD: `results/regime/summary.txt`, which `code/build/kill_regime --mode aggregate` derives from
  NEW: `results/regime/summary.txt`, which `julia --project=Kirigami Kirigami/apps/kill_regime.jl --mode aggregate` derives from
  ```

### `results/yield/BASIN.md` — 3 edits, 0 leave

- **L19**
  ```
  OLD: - Driver `code/apps/kill_basin.cpp`; aggregator `code/apps/kill_basin_agg.cpp`; figure
  NEW: - Driver `Kirigami/apps/kill_basin.jl`; aggregator `Kirigami/apps/kill_basin_agg.jl`; figure
  ```
- **L20**
  ```
  OLD:   `results/yield/plot_basin.py`. Nothing under `code/src` and nothing under `results/kill`
  NEW:   `results/yield/plot_basin.py`. Nothing under `Kirigami/src` and nothing under `results/kill`
  ```
- **L37**
  ```
  OLD:   `code/tests/test_design.cpp` ("design_range_max reproduces the K9c CSV row …") pass in
  NEW:   `Kirigami/test/test_design.jl` ("design_range_max reproduces the K9c CSV row …") pass in
  ```

### `results/yield/YIELD.md` — 3 edits, 0 leave

- **L14**
  ```
  OLD: `code/apps/kill_yield.cpp` (added to `code/CMakeLists.txt`) computes, per design, 45 raw
  NEW: `Kirigami/apps/kill_yield.jl` (added to `Kirigami/Project.toml`) computes, per design, 45 raw
  ```
- **L22**
  ```
  OLD:   same graphs and same σ as `code/apps/kill_k9c.cpp`, the Eq. (6) projection read from the
  NEW:   same graphs and same σ as `Kirigami/apps/kill_k9c.jl`, the Eq. (6) projection read from the
  ```
- **L252**
  ```
  OLD: | `code/apps/kill_yield.cpp` | feature extractor, three modes (`k9c`, `fresh`, `check`) |
  NEW: | `Kirigami/apps/kill_yield.jl` | feature extractor, three modes (`k9c`, `fresh`, `check`) |
  ```

### `export/hero/README.md` — 10 edits, 0 leave

- **L6**
  ```
  OLD: regression-locked in `code/tests/test_design.cpp`), certified valid with the
  NEW: regression-locked in `Kirigami/test/test_design.jl`), certified valid with the
  ```
- **L33**
  ```
  OLD: `tests/test_design.cpp` regression lock (`theta_max == eps_max ==
  NEW: `test/test_design.jl` regression lock (`theta_max == eps_max ==
  ```
- **L40**
  ```
  OLD:    sigma + `X_ini`, bit-identical to what `apps/kill_k9.cpp` and
  NEW:    sigma + `X_ini`, bit-identical to what `apps/kill_k9.jl` and
  ```
- **L41**
  ```
  OLD:    `tests/test_design.cpp`'s `k9_design(148)` build — this matters because
  NEW:    `test/test_design.jl`'s `k9_design(148)` build — this matters because
  ```
- **L46**
  ```
  OLD:    (not added to the CMake build, to avoid touching the shared `CMakeLists.txt`
  NEW:    (not added to the package, to avoid touching the shared `Project.toml`
  ```
- **L50–52**
  ```
  OLD:    clang++ -std=c++20 -O2 -arch arm64 \
            -I code/src -I code/apps -I /opt/homebrew/include -I /opt/homebrew/include/eigen3 \
            dump_hero_graph.cpp code/build/libkiri_core.a -o dump_hero_graph
  NEW: julia --project=Kirigami export/hero/dump_hero_graph.jl
  ```
- **L55**
  ```
  OLD:    `dump_hero_graph.cpp` just calls `kiri::kill::make_graph(148, 100, 800, 1400)`
  NEW:    `dump_hero_graph.jl` just calls `kiri::kill::make_graph(148, 100, 800, 1400)`
  ```
- **L96**
  ```
  OLD:    `--theta-frac` computes `theta_max` internally (`export/solid.hpp`'s
  NEW:    `--theta-frac` computes `theta_max` internally (`export/solid.jl`'s
  ```
- **L156**
  ```
  OLD: `dump_hero_graph.cpp` (the throwaway graph-dump utility from step 1) was
  NEW: `dump_hero_graph.jl` (the throwaway graph-dump utility from step 1) was
  ```
- **L158**
  ```
  OLD: the CMake build — it is a 20-line wrapper around
  NEW: the package — it is a 20-line wrapper around
  ```

### `export/hero2/README.md` — 15 edits, 0 leave

- **L23**
  ```
  OLD: `kiri_design --maximise-eps` runs `design_constrained` (`method/design.cpp`), whose stage-1
  NEW: `kiri_design --maximise-eps` runs `design_constrained` (`method/design.jl`), whose stage-1
  ```
- **L27**
  ```
  OLD: **stage-A margin-maximisation** (`method/range_embed.hpp`), warm-started from **t = 0**
  NEW: **stage-A margin-maximisation** (`method/range_embed.jl`), warm-started from **t = 0**
  ```
- **L32**
  ```
  OLD: (`dump_k9c_graph.cpp`, not committed) that replicates `kill_k9c.cpp`'s exact per-graph
  NEW: (`dump_k9c_graph.jl`, not committed) that replicates `kill_k9c.jl`'s exact per-graph
  ```
- **L47**
  ```
  OLD: (`reconstruct_k9c_graph.cpp`) rebuilds the *same deterministic, non-random* topology
  NEW: (`reconstruct_k9c_graph.jl`) rebuilds the *same deterministic, non-random* topology
  ```
- **L52**
  ```
  OLD: A third throwaway program (`characterize_k9c.cpp`) then calls `method::characterize` — the
  NEW: A third throwaway program (`characterize_k9c.jl`) then calls `method::characterize` — the
  ```
- **L56**
  ```
  OLD: added to the CMake build; they are graph-topology/measurement plumbing, not new method
  NEW: added to the package; they are graph-topology/measurement plumbing, not new method
  ```
- **L106**
  ```
  OLD:    `apps/kill_k9c.cpp` builds for `id=130`), same throwaway-utility pattern as hero's
  NEW:    `apps/kill_k9c.jl` builds for `id=130`), same throwaway-utility pattern as hero's
  ```
- **L107**
  ```
  OLD:    `dump_hero_graph.cpp`:
  NEW:    `dump_hero_graph.jl`:
  ```
- **L109–111**
  ```
  OLD:    clang++ -std=c++20 -O2 -arch arm64 \
            -I code/src -I code/apps -I /opt/homebrew/include -I /opt/homebrew/include/eigen3 \
            dump_hero2_input.cpp code/build/libkiri_core.a -o dump_hero2_input
  NEW: julia --project=Kirigami export/hero2/dump_hero2_input.jl
  ```
- **L114**
  ```
  OLD:    (`dump_hero2_input.cpp` calls `kiri::kill::make_graph(130, 100, 800, 1400)`, K9c's
  NEW:    (`dump_hero2_input.jl` calls `kiri::kill::make_graph(130, 100, 800, 1400)`, K9c's
  ```
- **L120–122**
  ```
  OLD:    clang++ -std=c++20 -O2 -arch arm64 \
            -I code/src -I code/apps -I /opt/homebrew/include -I /opt/homebrew/include/eigen3 \
            reconstruct_k9c_graph.cpp code/build/libkiri_core.a -o reconstruct_k9c_graph
  NEW: julia --project=Kirigami export/hero2/reconstruct_k9c_graph.jl
  ```
- **L135–137**
  ```
  OLD:    clang++ -std=c++20 -O2 -arch arm64 \
            -I code/src -I code/apps -I /opt/homebrew/include -I /opt/homebrew/include/eigen3 \
            characterize_k9c.cpp code/build/libkiri_core.a -o characterize_k9c
  NEW: julia --project=Kirigami export/hero2/characterize_k9c.jl
  ```
- **L165**
  ```
  OLD:    `--theta-frac` reports `theta_max=3.14159` at every fraction (export/solid.hpp's
  NEW:    `--theta-frac` reports `theta_max=3.14159` at every fraction (export/solid.jl's
  ```
- **L223**
  ```
  OLD: `dump_hero2_input.cpp`, `reconstruct_k9c_graph.cpp`, `characterize_k9c.cpp` (the
  NEW: `dump_hero2_input.jl`, `reconstruct_k9c_graph.jl`, `characterize_k9c.jl` (the
  ```
- **L225**
  ```
  OLD: the repository or the CMake build — they are graph-topology dump / fold-back /
  NEW: the repository or the package — they are graph-topology dump / fold-back /
  ```

### `review/constructive_review.md` — 1 edits, 0 leave

- **L17**
  ```
  OLD: **Mechanically: yes, it is a standard L-BFGS/softplus/log-barrier feasibility-then-proximity solve over two quadratic-in-t constraint families** (`code/src/method/convex_embed.hpp`): corner convexity (`cross_i(t) ≥ δ`) and split-inward sign (`q_e(t) ≥ δ'`). Nothing about the optimizer is novel — continuation on target, Gaussian restarts, geometric barrier decay are textbook.
  NEW: **Mechanically: yes, it is a standard L-BFGS/softplus/log-barrier feasibility-then-proximity solve over two quadratic-in-t constraint families** (`Kirigami/src/method/convex_embed.jl`): corner convexity (`cross_i(t) ≥ δ`) and split-inward sign (`q_e(t) ≥ δ'`). Nothing about the optimizer is novel — continuation on target, Gaussian restarts, geometric barrier decay are textbook.
  ```

### `review/negative_review.md` — 1 edits, 0 leave

- **L83**
  ```
  OLD: triangulations of uniform random points (`code/src/core/generators.cpp:496-580`) are a
  NEW: triangulations of uniform random points (`Kirigami/src/core/generators.jl:496-580`) are a
  ```

### `review/theory_review.md` — 1 edits, 0 leave

- **L24**
  ```
  OLD:    convention (`derivations/core.md:190-390`, Check C2 vs `kinematics.cpp::deploy()`, 1.42e−14). Acuña et al.
  NEW:    convention (`derivations/core.md:190-390`, Check C2 vs `kinematics.jl::deploy()`, ⟨JULIA:check_t1_t2:1.42e−14⟩). Acuña et al.
  ```

## Placeholders

Every `⟨JULIA:…⟩` in the edits above, and where its value comes from once the Julia suite and experiments have run. `<program>:<old value>` placeholders keep the C++ value inside the marker so the applier can compare.

| placeholder | source of the number |
|---|---|
| `⟨JULIA:tests⟩` | `Pkg.test()` output, whole suite |
| `⟨JULIA:derivation_tests⟩` | `Pkg.test()` output, `derivation_tests.jl` testset |
| `⟨JULIA:tests:C4⟩` | `Pkg.test()` filtered to the K7 C4 testset |
| `⟨JULIA:check_r3:mismatches⟩` | `check_r3.jl` mismatch count (was 0) |
| `⟨JULIA:check_r2:violations⟩` | `check_r2.jl` R2-D violation count (was 0) |
| `⟨JULIA:kill_k3a_recheck:7 of 7⟩` | K3a recheck violator count |
| `⟨JULIA:k9c:307 of 400⟩` | `results/kill/k9c/summary.txt` after the Julia K9c run |
| `⟨JULIA:derivation_tests:R6-a1⟩` | `derivation_tests.jl` R6-a1 failure count / pair count |
| `⟨JULIA:derivation_tests:R6-c1⟩` | `derivation_tests.jl` R6-c1 |
| `⟨JULIA:check_t1_t2:1.60e−14⟩` | measured by `check_t1_t2` on the C++; re-measure with the Julia port (derivations/core.md:1984) |
| `⟨JULIA:check_t1_t2:1.42e−14⟩` | measured by `check_t1_t2` on the C++; re-measure with the Julia port (derivations/core.md:1984) |
| `⟨JULIA:check_t1_t2:9.73e−14⟩` | measured by `check_t1_t2` on the C++; re-measure with the Julia port (derivations/core.md:1984) |
| `⟨JULIA:check_t1_t2:2.08e−14⟩` | measured by `check_t1_t2` on the C++; re-measure with the Julia port (derivations/core.md:1984) |
| `⟨JULIA:check_t1_t2:4.46e−14⟩` | measured by `check_t1_t2` on the C++; re-measure with the Julia port (derivations/core.md:1984) |
| `⟨JULIA:check_t1_t2:1.83e−14⟩` | measured by `check_t1_t2` on the C++; re-measure with the Julia port (derivations/core.md:1984) |
| `⟨JULIA:check_t1_t2:1.68e−10⟩` | measured by `check_t1_t2` on the C++; re-measure with the Julia port (derivations/core.md:1984) |
| `⟨JULIA:check_t4_t5:2.09e−6⟩` | measured by `check_t4_t5` on the C++; re-measure with the Julia port (derivations/core.md:1985) |
| `⟨JULIA:check_t4_t5:8.88e−16⟩` | measured by `check_t4_t5` on the C++; re-measure with the Julia port (derivations/core.md:1985) |
| `⟨JULIA:check_t4_t5:4.37e−15⟩` | measured by `check_t4_t5` on the C++; re-measure with the Julia port (derivations/core.md:1985) |
| `⟨JULIA:check_r2:3.11e−15⟩` | measured by `check_r2` on the C++; re-measure with the Julia port (derivations/core.md:1986) |
| `⟨JULIA:check_r2:8.88e−16⟩` | measured by `check_r2` on the C++; re-measure with the Julia port (derivations/core.md:1986) |
| `⟨JULIA:check_t5_adm:1.55e−15⟩` | measured by `check_t5_adm` on the C++; re-measure with the Julia port (derivations/core.md:2261) |
| `⟨JULIA:derivation_tests:1.77e−14⟩` | measured by `derivation_tests` on the C++; re-measure with the Julia port (derivations/core.md:2263) |
| `⟨JULIA:derivation_tests:3.53e−14⟩` | measured by `derivation_tests` on the C++; re-measure with the Julia port (derivations/core.md:2263) |
| `⟨JULIA:check_t1_t2:1.4e−14⟩` | measured by `check_t1_t2` on the C++; re-measure with the Julia port (derivations/check.md:66) |

**25 distinct placeholders.**

## Summary

| file | proposed edits | leave |
|---|---:|---:|
| `docs/techreport/techreport.tex` | 56 | 1 |
| `docs/paper/paper.tex` | 1 | 0 |
| `REPORT.md` | 35 | 0 |
| `STATE.md` | 5 | 26 |
| `IDEA.md` | 2 | 0 |
| `derivations/core.md` | 51 | 0 |
| `derivations/check.md` | 41 | 0 |
| `derivations/check_lemmas.md` | 24 | 0 |
| `derivations/lemmas.md` | 25 | 0 |
| `derivations/lemmas_statements.md` | 4 | 0 |
| `results/kill/KILL_REPORT.md` | 106 | 3 |
| `results/kill/jitter/cert_diagnosis.md` | 11 | 4 |
| `results/kill/k8a/recheck.md` | 7 | 0 |
| `results/kill/native200/NATIVE200_FINAL.md` | 1 | 0 |
| `results/core_validation/rank_claim.md` | 4 | 0 |
| `results/core_validation/referee_fix.md` | 10 | 0 |
| `results/final/e1/E1.md` | 4 | 0 |
| `results/final/figures/README.md` | 1 | 0 |
| `results/regime/REGIME.md` | 1 | 0 |
| `results/yield/BASIN.md` | 3 | 0 |
| `results/yield/YIELD.md` | 3 | 0 |
| `export/hero/README.md` | 10 | 0 |
| `export/hero2/README.md` | 15 | 0 |
| `review/constructive_review.md` | 1 | 0 |
| `review/negative_review.md` | 1 | 0 |
| `review/theory_review.md` | 1 | 0 |
| `docs/poster/poster.tex` | 0 | 0 (no C++ mention) |
| **total** | **423** | **34** |