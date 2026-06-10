# C++ mentions in the migrated prose (Markdown / TeX / bib)

Generated 2026-09-19 by scanning every `.md`, `.tex`, `.bib` copied from `kirigami-experiments` into this repo
(`PORTING.md` excluded; it is the Julia repo's own). Pattern (regex, case-sensitive except where noted):

```
C\+\+ | \bEigen\b | Eigen:: | Eigen3 | doctest | nlohmann | CMake|cmake | \bclang\b | \.cpp\b | \.hpp\b | kiri_tests | derivation_tests
| kill_[A-Za-z0-9_]+ | code/src | code/apps | code/tests | code/web | wasm|WASM | [Ee]mscripten | [Ww]eb explorer | WebAssembly
```

**942 matching lines in 61 files.** Nothing here has been edited; rewording to Julia happens only after the port reproduces the numbers.

Legend for the per-file note: **KEEP** = describes third-party C++ (the authors' code, the `baseline/` harness) and stays as is except for `code/...` path references; **HISTORICAL** = agent spec / ideation record from the C++ era, factual at the time — reword only the forward-looking statements or add a one-line 'now Julia' header, do not rewrite history.

## Summary by file

| file | lines | what will need rewording |
|---|---:|---|
| `docs/paper/paper.tex` | 1 | one file path reference in a \plan{} slot: `code/src/method/range_embed.{hpp,cpp}`. |
| `docs/techreport/techreport.tex` | 74 | toolchain paragraphs (l.82, 429–434, 3463–3473: 'C++20 with CMake', Eigen/nlohmann/doctest); test counts 109 tests / 17,854 assertions and 33 / 89,074 (l.3472–3473, 3913–3914); software-layout table (l.3489–3495); reproduce-everything command block `./code/build/kill_*` (l.3912–3948); web explorer / WebAssembly / emscripten section (l.3495, 3506, 3632–3637); many \file{...cpp/hpp} algorithm captions and refs; 'a doctest checks' phrasing. |
| `REPORT.md` | 36 | test counts (`kiri_tests` 103 / 17,504; `derivation_tests` 33 / 89,074, l.316–323, 561–562); software-layout table `code/src/{core,method,export}`, `code/tests`, `code/apps/kiri_*.cpp`, `code/apps/kill_*.cpp` (l.531–536); cmake build + `./code/build/kill_*` reproduce block (l.555–597); `code/tests/test_design.cpp` ref. |
| `STATE.md` | 31 | status log: headline '109 tests / 17,854' (l.7); toolchain line (l.13); per-phase test counts (28/6,436; 71/16,488; 65/16,229; 25/29,033; 31/89,0xx; 33/89,074; kiri_tests 16,622); directive D3 'ALL coding tasks in C++' (l.67); E1 xcode-select/cmake environment note (l.173); build-tree presence (l.184); baseline parity (l.58, KEEP). |
| `IDEA.md` | 2 | two file path refs: `code/src/method/range_embed.{hpp,cpp}` (l.171), `code/tests/test_design.cpp` (l.214). |
| `derivations/check.md` | 46 | clang++ build lines for `derivation_tests` (x2); doctest test counts (25 cases / 29,033 assertions); many `contact.hpp/.cpp`, `kinematics.cpp`, `deploy_basis.hpp`, `collision.cpp` refs; scratch `check_*.cpp` programs (Programs-run table). |
| `derivations/check_lemmas.md` | 28 | file path refs (`*.hpp`, `apps/kill_b3.cpp`, `kill_k7.cpp`); clang++ build line + doctest output block (37 cases / 149,115 assertions) for `derivation_tests`; `code/src|apps|tests|CMakeLists` non-modification statements. |
| `derivations/core.md` | 51 | file refs (`kinematics.cpp`, `deploy_basis.hpp`, `contact.hpp`, `collision.hpp`, `cut.hpp`); test counts (25/29,033; 30/89,044); scratch `check_*.cpp` programs and their clang++ build lines (Programs-run table l.1984–1991, 2261–2263). |
| `derivations/lemmas.md` | 42 | full clang++ build lines listing `code/src/core/*.cpp` (l.461–465, 788–790, 930–940); Xcode toolchain workaround note; `zero_plus.hpp`/`cut.hpp`/`contact.cpp`/`periodic_jacobian.hpp`/`kill_k7.cpp` refs; `code/tests were not run` remark. |
| `derivations/lemmas_statements.md` | 4 | file refs only (`cut.hpp`, `contact.cpp`, `periodic_jacobian.hpp`, `kill_k7.cpp::shape_point`). |
| `results/core_validation/rank_claim.md` | 4 | `kiri_tests` and `tests/test_holes.cpp`, `tests/test_rank_checks.cpp` refs. |
| `results/core_validation/referee_fix.md` | 11 | `code/src/core/collision.{hpp,cpp}` scope; `code/tests/test_collision.cpp` failing-first test; test-count table `kiri_tests` 103 / 17,504 and `derivation_tests` 33 / 89,074 (l.98–99); `check_b4_93.cpp`. |
| `results/final/e1/E1.md` | 4 | driver `code/apps/kill_e1.cpp`; `kill_common::` population refs; `cmake --build . --target kill_e1`. |
| `results/final/figures/README.md` | 1 | one phrase: 're-run any C++ solver'. |
| `results/kill/KILL_REPORT.md` | 114 | the largest: 'C++ drivers in `code/apps/kill_*.cpp`' (l.4); per-experiment Driver: lines naming `code/apps/kill_*.cpp`/`dbg_*.cpp`; method file refs `code/src/method/*.{hpp,cpp}`; 'new doctest in `code/tests/test_method.cpp`' + running suite counts (53/16,004; 65/16,229; 66/16,283; 71/16,488; 80/16,785); four cmake + `./code/build/kill_*` reproduce blocks (l.986–997, 1164–1171, 1422–1440, 1749–1751); `baseline/native/src/main.cpp --dumpdir` (KEEP, baseline). |
| `results/kill/jitter/cert_diagnosis.md` | 16 | 'C++ code in `code/`'; `contact.cpp` refs, `code/apps/dbg_cert.cpp`; regression test in `code/tests/test_method.cpp`; cmake line + test counts 66 / 16,283 and 31 / 89,055 (l.175–180); `kill_jitter`/`kill_k2a`/`kill_k5`/`kill_k6` binaries. |
| `results/kill/k8a/recheck.md` | 8 | `expansive_cone.cpp` line refs; `code/apps/kill_k8a{,_recheck}.cpp`; cmake + `./code/build/kiri_tests` + `kill_k8a*` command block (l.221–229). |
| `results/kill/native200/NATIVE200_FINAL.md` | 1 | one ref: `code/apps/native200_merge.cpp` produces the report. |
| `results/regime/REGIME.md` | 1 | one ref: `code/build/kill_regime --mode aggregate` binary. |
| `results/yield/BASIN.md` | 3 | `code/apps/kill_basin{,_agg}.cpp` drivers; `code/src` non-modification; `code/tests/test_design.cpp`. |
| `results/yield/YIELD.md` | 3 | `code/apps/kill_yield.cpp` (added to `code/CMakeLists.txt`); `code/apps/kill_k9c.cpp`; artifact table row. |
| `export/hero/README.md` | 12 | clang++ build line for throwaway `dump_hero_graph.cpp` (l.50–52); `tests/test_design.cpp` regression-lock refs; `apps/kill_k9.cpp`, `export/solid.hpp` refs; 'not added to the CMake build' remarks. |
| `export/hero2/README.md` | 21 | three clang++ build lines for throwaway `dump_hero2_input.cpp`/`reconstruct_k9c_graph.cpp`/`characterize_k9c.cpp` (l.109–137); `method/design.cpp`, `range_embed.hpp`, `kill_k9c.cpp`, `export/solid.hpp` refs; CMake remarks. |
| `review/constructive_review.md` | 1 | one file path ref `code/src/method/convex_embed.hpp` (l.17). |
| `review/negative_review.md` | 1 | one file:line ref `code/src/core/generators.cpp:496-580`. |
| `review/theory_review.md` | 1 | one ref `kinematics.cpp::deploy()`. |
| `notes/paper_2026.md` | 1 | KEEP — quotes the 2026 paper's own 'WebAssembly backend compiled from C++' sentence (third-party). |
| `notes/repo_2025.md` | 116 | KEEP — analysis of the AUTHORS' 2025 C++/emscripten repo (`kiri_mesh.cpp`, `multi_grid.cpp`, `bind.cpp`, their CMakeLists, Eigen); third-party, not our code. |
| `notes/screen_r2.md` | 3 | `periodic_jacobian.hpp` ref; `kill_k8f`/`kill_i2` proposed driver names. |
| `ideas/persona_adversary.md` | 9 | HISTORICAL (ideation) — scratch `check.cpp`/`rig_check.cpp` refs, `code/src/core/` assumption, 'C++-only stack'/Eigen effort estimates, authors' WebAssembly mention. |
| `ideas/persona_geometer.md` | 13 | HISTORICAL — scratch `check*.cpp` programs, hypothetical `core/*.cpp` routine names, `Eigen::ColPivHouseholderQR`, 'code/src/core does not exist yet'. |
| `ideas/persona_optimizer.md` | 5 | HISTORICAL — `code/src/core/*.hpp`, `kinematics.hpp` convention refs; 'substantial piece of C++'; 'C++ CSV dump'. |
| `ideas/persona_rigidity.md` | 9 | HISTORICAL — standalone C++ compile against `code/src/core/*.cpp` with clang++; `CMakeLists.txt`/`collision.cpp` not-yet-existing remarks; Eigen GSVD remark. |
| `ideas/ranking.md` | 11 | 'All C++17/20 per directive D3' (l.462) and 'SDP-in-C++' (l.677); `collision.hpp` refs; `code/src/core` refs; authors' WebAssembly mention. |
| `ideas/ranking_r2.md` | 21 | `kill_k8a`/`kill_t1`/`kill_b3`/`kill_b2`/`kill_k8b`/`kill_k8h` driver names (binary names -> apps/*.jl); `zero_plus.hpp`, `mobility.hpp`, `contact.hpp` refs; `kill_k7.cpp`. |
| `ideas/round2_adversary.md` | 40 | HISTORICAL — `code/src/method/*.hpp` inventory, proposed `kill_k8a..k8j.cpp` drivers, `kill_common.hpp`; kill-rule table with `kill_k8*` names. |
| `ideas/round2_inverse.md` | 23 | HISTORICAL — `code/src/method/*.hpp`, `code/apps/kill_k7.cpp`; proposed `kill_i1..i10.cpp` drivers; '≈ 90 lines of C++' effort notes; `code/src/export`. |
| `ideas/round2_theorist.md` | 17 | HISTORICAL — `code/src/method/*.hpp` refs; proposed extensions of `kill_k1a/k5/k6.cpp`, new `kill_k8/k9/k10.cpp`; `mobility.cpp`, `generators.cpp`, `mesh.cpp`. |
| `ideas/round2_theorist_b.md` | 36 | HISTORICAL — `code/src/method/{zero_plus,periodic_jacobian}.hpp`, scratch `check_b1.cpp`, proposed `kill_b1..b9.cpp` drivers, `kill_k6.cpp`/`kill_k7.cpp`. |
| `specs/builder_baseline.md` | 6 | HISTORICAL spec — KEEP (baseline/native CMake project, clang, bind.cpp); refers to `code/README.md` JSON contract. |
| `specs/builder_core.md` | 15 | HISTORICAL spec — title 'C++ reimplementation'; CMake targets `kiri_core`/`kiri_tests`; Eigen/nlohmann/doctest deps; `*.hpp` module names; build command. |
| `specs/builder_export.md` | 3 | HISTORICAL spec — `code/src/export/`, `code/apps/kiri_export.cpp`, `code/tests/test_export.cpp`, doctest, CMakeLists, build command. |
| `specs/builder_method.md` | 4 | HISTORICAL spec — `code/src/method/`, `collision.hpp`, doctest/CMakeLists. |
| `specs/common_preamble.md` | 2 | HISTORICAL — directive D3 'ALL code is C++17/20 with CMake' (l.8) and toolchain line clang++/cmake/Eigen/nlohmann/doctest versions (l.9). |
| `specs/critic.md` | 2 | HISTORICAL — 'demonstrable in C++', `code/src/core`. |
| `specs/deriver.md` | 1 | HISTORICAL — 'throwaway numeric checks in C++ under derivations/scratch/' (l.3). |
| `specs/experimenter_final.md` | 2 | HISTORICAL — 'all C++', `collision.hpp`. |
| `specs/experimenter_k7.md` | 3 | HISTORICAL — 'All C++ (D3)', `code/src/core/generators`. |
| `specs/experimenter_kill.md` | 3 | HISTORICAL — 'small C++ apps under code/apps/' (l.3), `collision.hpp` bisection (l.14), 'add doctest cases' (l.17). |
| `specs/ideator.md` | 1 | HISTORICAL — 'Demonstrable in C++'. |
| `specs/ideator_round2.md` | 2 | HISTORICAL — 'AVAILABLE in code/src/method' (l.6); 'existing drivers (kill_k*.cpp, kiri_analyze)' (l.15). |
| `specs/m2_checker_lemmas.md` | 5 | HISTORICAL — `code/src/method/*.hpp`, `code/tests/derivation_tests.cpp` (33 cases / ~89k assertions), `kill_common.hpp`, `code/apps/kill_k7.cpp`. |
| `specs/m2_deriver_lemmas.md` | 7 | HISTORICAL — 'scratch C++ checks', `code/src/method/*.hpp`, `derivations/scratch/check_l1/l2.cpp` compiled against `code/src`; `code/src`/`code/tests` non-modification. |
| `specs/m2_experimenter_basin.md` | 5 | HISTORICAL — Xcode make workaround, `code/src/method/design.hpp`, `code/apps/kill_yield.cpp`/`kill_basin.cpp`. |
| `specs/m2_experimenter_native200.md` | 6 | HISTORICAL — `code/apps/kill_native200.cpp`, `kill_common.hpp`, `native200_merge.cpp`, baseline/native CLI (KEEP). |
| `specs/m2_experimenter_regime.md` | 9 | HISTORICAL — xcode-select/cmake environment note, `code/apps/kill_k9c.cpp`, `kill_common.hpp`, `code/src/method/design.hpp`, `code/apps/kill_regime.cpp`, `kill_scaling.cpp`, CMakeLists. |
| `specs/m2_experimenter_yield.md` | 6 | HISTORICAL — `kill_common.hpp`, `code/src/method/design.hpp`, `code/src/core/*.hpp`, `code/apps/kill_yield.cpp` + CMakeLists. |
| `specs/m2_propagator_a.md` | 6 | HISTORICAL — cmake/xcode environment, `derivation_tests` CMake target, `test_reference_cases.cpp`, `kiri_tests` counts, `code/build` hygiene. |
| `baseline/README.md` | 24 | KEEP — third-party C++ build notes (CMake/Eigen/emscripten stub); only the `code/src/core` path reference (l.5) needs rewording. |
| `baseline/parity.md` | 4 | KEEP — describes the authors' C++ baseline; only `code/src/core` path reference (title, l.1/5) needs -> `Kirigami/src/core`. |
| `baseline/tuttekiri/code/README.md` | 5 | KEEP verbatim — upstream authors' README (WebAssembly build); never edit. |

## Full listing (file:line — line text)

### `docs/paper/paper.tex` (1)

_one file path reference in a \plan{} slot: `code/src/method/range_embed.{hpp,cpp}`._

```
docs/paper/paper.tex:345: \plan{The algorithm, from \texttt{code/src/method/range\_embed.\{hpp,cpp\}},
```

### `docs/techreport/techreport.tex` (74)

_toolchain paragraphs (l.82, 429–434, 3463–3473: 'C++20 with CMake', Eigen/nlohmann/doctest); test counts 109 tests / 17,854 assertions and 33 / 89,074 (l.3472–3473, 3913–3914); software-layout table (l.3489–3495); reproduce-everything command block `./code/build/kill_*` (l.3912–3948); web explorer / WebAssembly / emscripten section (l.3495, 3506, 3632–3637); many \file{...cpp/hpp} algorithm captions and refs; 'a doctest checks' phrasing._

```
docs/techreport/techreport.tex:82: soundness violations of the certificate. All code is C++20 with CMake; the two test
docs/techreport/techreport.tex:429: The reimplementation is C++20 + CMake, header-only dependencies (Eigen, nlohmann/json,
docs/techreport/techreport.tex:430: doctest), no external optimiser --- the L-BFGS with backtracking line search in
docs/techreport/techreport.tex:431: \file{core/optimize.hpp} is ours. The library layout is in Section~\ref{sec:software}.
docs/techreport/techreport.tex:434: (\file{baseline/native}, libigl 2.5, Optiz, Eigen 3.4 via FetchContent, their \file{bind.cpp}
docs/techreport/techreport.tex:535: \textbf{[F, \file{core/kinematics.cpp}] The code's sign convention.} Face $f$ is transformed
docs/techreport/techreport.tex:634: \textbf{[N] Path-independence.} \file{derivations/scratch/check\_t1\_t2.cpp}, check C1: build
docs/techreport/techreport.tex:668: \caption{Assembling the deploy basis $(C,S)$ --- \file{method/deploy\_basis.cpp}}
docs/techreport/techreport.tex:689: \textbf{A naming clash, recorded.} The header of \file{deploy\_basis.hpp} writes
docs/techreport/techreport.tex:697: independent BFS, against \file{kinematics.cpp::deploy()} at
docs/techreport/techreport.tex:1181: \caption{Exact $\Tmax$ (T4.2$''$) --- \file{method/contact.cpp::exact\_theta\_max\_overlap}}
docs/techreport/techreport.tex:1288: \caption{Broad phase --- \file{method/contact.cpp::candidate\_pairs}}
docs/techreport/techreport.tex:1584: This is exactly what \file{contact.cpp} implements --- before round 6 the derivation demanded
docs/techreport/techreport.tex:1669: \textbf{[N] The measured counterexample} (\file{derivations/scratch/check\_t5\_adm.cpp},
docs/techreport/techreport.tex:1725: \file{contact.cpp}; a doctest fails without it and passes with it.
docs/techreport/techreport.tex:1746: bisection referee (\file{results/final/e1/E1.md}, driver \file{code/apps/kill\_e1.cpp}, which
docs/techreport/techreport.tex:2149: \file{code/src/method/contact.cpp:339} --- keep a root only if it is smaller than the
docs/techreport/techreport.tex:2303: \file{zero\_plus\_q} and \file{zero\_plus\_corner\_margin} at the $\sigma$ ray, in a doctest,
docs/techreport/techreport.tex:2305: one free translation pair per component of $\Gamma$; a doctest checks
docs/techreport/techreport.tex:2447: A doctest pins this. \emph{Any claim of the form ``authored tilings pass 8/8'' is really
docs/techreport/techreport.tex:2717: surrogate. Two doctests check exactly those two claims
docs/techreport/techreport.tex:2718: (\file{code/tests/test\_range\_embed.cpp}), and two more finite-difference the analytic
docs/techreport/techreport.tex:2737: objective of \file{range\_opt.hpp}, whose gradients are the implicit-function-theorem
docs/techreport/techreport.tex:2774: \caption{\texttt{design\_range\_max} --- the deliverable (K9c). \file{method/range\_embed.cpp}, \file{method/design.cpp}}
docs/techreport/techreport.tex:2835: iterates an unordered container or touches global state, and Eigen is built without OpenMP so
docs/techreport/techreport.tex:2837: \file{tests/test\_design.cpp} checks that. Two consequences worth stating:
docs/techreport/techreport.tex:2843: Eigen version or assembly order only if $\Phi$ is \textbf{carried along, not recomputed}.
docs/techreport/techreport.tex:2845: \item \file{apps/kill\_k9c.cpp} is the one place that was genuinely non-deterministic: it
docs/techreport/techreport.tex:2855: \file{tests/test\_design.cpp} locks four rows to $10^{-9}$: \file{delaunay\_130}
docs/techreport/techreport.tex:3047: \file{method/periodic\_jacobian.cpp} builds the genuine quotient: one face per translation
docs/techreport/techreport.tex:3463: Everything is C++17/20 with CMake. Python appears only for matplotlib plotting of CSV/JSON
docs/techreport/techreport.tex:3464: dumped by the C++ programs. Dependencies are header-only (Eigen, nlohmann/json, doctest); the
docs/techreport/techreport.tex:3471: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
docs/techreport/techreport.tex:3472: ./code/build/kiri_tests          # 109 cases, 17,854 assertions, 0 failures
docs/techreport/techreport.tex:3473: ./code/build/derivation_tests    # 33 cases, 89,074 assertions, 0 failures
docs/techreport/techreport.tex:3489: core pipeline (2026 reimplementation) & \file{code/src/core/}: \file{mesh} \file{cut} \file{holes} \file{orientation} \file{tutte\_auxetic} \file{kinematics} \file{collision} \file{optimize} \file{rank\_checks} \file{generators} \\
docs/techreport/techreport.tex:3490: method layer & \file{code/src/method/}: \file{deploy\_basis} \file{contact} \file{design} \file{range\_opt} \file{zero\_plus} \file{convex\_embed} \file{range\_embed} \file{mobility} \file{expansive\_cone} \file{periodic\_jacobian} \file...
docs/techreport/techreport.tex:3491: fabrication export & \file{code/src/export/}: \file{layout} \file{solid} \file{svg} \file{stl} \file{threemf} \file{xml} \file{zip} \file{material} \\
docs/techreport/techreport.tex:3492: tests & \file{code/tests/}: \file{test\_mesh\_cut} \file{test\_holes} \file{test\_system} \file{test\_kinematics} \file{test\_collision} \file{test\_method} \file{test\_design} \file{test\_range\_embed} \file{test\_rank\_checks} \file{te...
docs/techreport/techreport.tex:3493: CLIs & \file{code/apps/}: \file{kiri\_gen} \file{kiri\_analyze} \file{kiri\_deploy} \file{kiri\_design} \file{kiri\_export} \file{kiri\_sweep} \file{kiri\_reference} \\
docs/techreport/techreport.tex:3494: experiment drivers & \file{code/apps/kill\_*}: \file{k1a} \file{k1b} \file{k1c} \file{k2a} \file{k2b} \file{k2c} \file{k3a} \file{k5} \file{k6} \file{k7} \file{k8a} \file{k9} \file{k9b} \file{k9c} \file{b3} \file{b4} \file{t1} \file{jitt...
docs/techreport/techreport.tex:3495: web explorer & \file{code/web/}: \file{bindings.cpp} \file{build.sh} \file{index.html} \file{app.html} \file{kiri.js} \\
docs/techreport/techreport.tex:3506: the web explorer and compared byte for byte:
docs/techreport/techreport.tex:3524: \file{method/design.hpp} is the single entry point, so that ``certified'' means one thing
docs/techreport/techreport.tex:3632: \subsection{The web explorer}
docs/techreport/techreport.tex:3634: \file{code/web/} compiles \file{src/core} $+$ \file{src/method} $+$ \file{src/export} to
docs/techreport/techreport.tex:3635: WebAssembly via emscripten (\file{SINGLE\_FILE=1} embeds the wasm as base64 inside the
docs/techreport/techreport.tex:3637: filesystem; the native CMake build is untouched). Everything crossing the JS boundary is a
docs/techreport/techreport.tex:3820: \file{expansive\_cone.hpp}; the C-IRIS line supplies the general pattern of
docs/techreport/techreport.tex:3912: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
docs/techreport/techreport.tex:3913: ./code/build/kiri_tests           # 109 cases, 17854 assertions, 0 failures
docs/techreport/techreport.tex:3914: ./code/build/derivation_tests     # 33 cases, 89074 assertions, 0 failures
docs/techreport/techreport.tex:3919: cmake --build code/build --target kill_e1
docs/techreport/techreport.tex:3926: ./code/build/kill_k1b  --out results/kill/k1b          # ~34 s
docs/techreport/techreport.tex:3927: ./code/build/kill_k1a  --out results/kill/k1a          # ~40 s
docs/techreport/techreport.tex:3928: ./code/build/kill_k2a  --out results/kill/k2a          # ~7 s
docs/techreport/techreport.tex:3929: ./code/build/kill_k2c  --out results/kill/k2c          # ~3 min
docs/techreport/techreport.tex:3930: ./code/build/kill_k3a  --out results/kill/k3a          # ~25 min
docs/techreport/techreport.tex:3931: ./code/build/kill_k3a_recheck --out results/kill/k3a   # ~4 min
docs/techreport/techreport.tex:3932: ./code/build/kill_k5   --out results/kill/k5 --n 200 --maxf 800    # ~5 min
docs/techreport/techreport.tex:3933: ./code/build/kill_k1c  --out results/kill/k1c          # ~31 s   (needs baseline/native)
docs/techreport/techreport.tex:3934: ./code/build/kill_k2b  --out results/kill/k2b --maxf 160   # ~26 min (needs native)
docs/techreport/techreport.tex:3935: ./code/build/kill_f23  --out results/kill/f23          # ~1 s
docs/techreport/techreport.tex:3936: ./code/build/kill_k6   --out results/kill/k6
docs/techreport/techreport.tex:3937: ./code/build/kill_t1   --out results/kill/t1           # ~7 s off the K6 cache
docs/techreport/techreport.tex:3938: ./code/build/kill_jitter --shard I --nshards 12 --out results/kill/jitter  # ~2 s
docs/techreport/techreport.tex:3939: ./code/build/kill_k7   --stage main --shard I --nshard 12 --out results/kill/k7
docs/techreport/techreport.tex:3940: ./code/build/kill_b3   --out results/kill/b3
docs/techreport/techreport.tex:3941: ./code/build/kill_b4   --n 200 --nshards 12 --shard I \
docs/techreport/techreport.tex:3943: ./code/build/kill_k8a  --n 100 --x0 --dual-iters 100000 --shard I --nshard 8 \
docs/techreport/techreport.tex:3945: ./code/build/kill_k9   --out results/kill/k9
docs/techreport/techreport.tex:3946: ./code/build/kill_k9b  --out results/kill/k9b
docs/techreport/techreport.tex:3947: ./code/build/kill_k9c  --out results/kill/k9c     # 12 shards; --aggregate to merge
docs/techreport/techreport.tex:3948: ./code/build/kill_native200 --n 200 --maxf 800 --timeout 600 \
```

### `REPORT.md` (36)

_test counts (`kiri_tests` 103 / 17,504; `derivation_tests` 33 / 89,074, l.316–323, 561–562); software-layout table `code/src/{core,method,export}`, `code/tests`, `code/apps/kiri_*.cpp`, `code/apps/kill_*.cpp` (l.531–536); cmake build + `./code/build/kill_*` reproduce block (l.555–597); `code/tests/test_design.cpp` ref._

```
REPORT.md:42: | C9 | Two-stage range-maximising constrained embedding; 307 of 400 designs deployable against 0 for every baseline | **algorithm** | `results/final/figures/summary_table.md`, `results/kill/k9c/`, `code/src/method/range_embed.{hpp,cpp}` |
REPORT.md:271: `code/tests/test_design.cpp` at 1e−9. Exported closed, at `Θ_max/2` and at `0.9 Θ_max` as
REPORT.md:316:   each time, ending with **zero unresolved disagreements**. `derivation_tests`: **33 cases /
REPORT.md:323: * **Gate 8 (Build).** `kiri_tests`: **103 cases / 17,504 assertions, 0 failures** (re-run by
REPORT.md:531: | core pipeline (2026 reimplementation) | `code/src/core/{mesh,cut,holes,orientation,tutte_auxetic,kinematics,collision,optimize,rank_checks,generators}.{hpp,cpp}` |
REPORT.md:532: | method | `code/src/method/{deploy_basis,contact,design,range_opt,zero_plus,convex_embed,range_embed,mobility,expansive_cone,periodic_jacobian,budget}.{hpp,cpp}` |
REPORT.md:533: | export | `code/src/export/{layout,solid,svg,stl,threemf,xml,zip,material}.{hpp,cpp}` |
REPORT.md:534: | tests | `code/tests/{test_mesh_cut,test_holes,test_system,test_kinematics,test_collision,test_method,test_design,test_range_embed,test_rank_checks,test_reference_cases,test_export,derivation_tests}.cpp` |
REPORT.md:535: | CLIs | `code/apps/{kiri_gen,kiri_analyze,kiri_deploy,kiri_design,kiri_export,kiri_sweep,kiri_reference}.cpp` |
REPORT.md:536: | kill drivers | `code/apps/kill_{k1a,k1b,k1c,k2a,k2b,k2c,k3a,k5,k6,k7,k8a,k9,k9b,k9c,b3,b4,t1,jitter,e1,native200,f23}.cpp` |
REPORT.md:555: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
REPORT.md:561: ./code/build/kiri_tests           # 103 cases, 17504 assertions, 0 failures
REPORT.md:562: ./code/build/derivation_tests     # 33 cases, 89074 assertions, 0 failures
REPORT.md:568: cmake --build code/build --target kill_e1
REPORT.md:576: ./code/build/kill_k1b  --out results/kill/k1b          # ~34 s
REPORT.md:577: ./code/build/kill_k1a  --out results/kill/k1a          # ~40 s
REPORT.md:578: ./code/build/kill_k2a  --out results/kill/k2a          # ~7 s
REPORT.md:579: ./code/build/kill_k2c  --out results/kill/k2c          # ~3 min
REPORT.md:580: ./code/build/kill_k3a  --out results/kill/k3a          # ~25 min
REPORT.md:581: ./code/build/kill_k3a_recheck --out results/kill/k3a   # ~4 min
REPORT.md:582: ./code/build/kill_k5   --out results/kill/k5 --n 200 --maxf 800   # ~5 min
REPORT.md:583: ./code/build/kill_k1c  --out results/kill/k1c          # ~31 s   (needs baseline/native)
REPORT.md:584: ./code/build/kill_k2b  --out results/kill/k2b --maxf 160          # ~26 min (needs baseline/native)
REPORT.md:585: ./code/build/kill_f23  --out results/kill/f23          # ~1 s
REPORT.md:586: ./code/build/kill_k6   --out results/kill/k6
REPORT.md:587: ./code/build/kill_k7   --out results/kill/k7
REPORT.md:588: ./code/build/kill_t1   --out results/kill/t1           # ~7 s off the K6 cache
REPORT.md:589: ./code/build/kill_b3   --out results/kill/b3
REPORT.md:590: ./code/build/kill_b4   --out results/kill/b4
REPORT.md:591: ./code/build/kill_k8a  --out results/kill/k8a
REPORT.md:592: ./code/build/kill_k8a_recheck --out results/kill/k8a
REPORT.md:593: ./code/build/kill_jitter --out results/kill/jitter
REPORT.md:594: ./code/build/kill_k9   --out results/kill/k9
REPORT.md:595: ./code/build/kill_k9b  --out results/kill/k9b
REPORT.md:596: ./code/build/kill_k9c  --out results/kill/k9c          # 12 shards; --aggregate --nshards 12 to merge
REPORT.md:597: ./code/build/kill_native200 --out results/kill/native200           # in progress
```

### `STATE.md` (31)

_status log: headline '109 tests / 17,854' (l.7); toolchain line (l.13); per-phase test counts (28/6,436; 71/16,488; 65/16,229; 25/29,033; 31/89,0xx; 33/89,074; kiri_tests 16,622); directive D3 'ALL coding tasks in C++' (l.67); E1 xcode-select/cmake environment note (l.173); build-tree presence (l.184); baseline parity (l.58, KEEP)._

```
STATE.md:7: - 2026-09-05 00:30: ALL MISSION DELIVERABLES ON DISK AND COMMITTED. IDEA.md, derivations/ (6 rounds), code/ (109 tests / 17,854; derivation 33 / 89,074), results/ (400-design study, figures, E1 3,113), export/hero + hero2 (SVG/3MF/STL), ...
STATE.md:13: - Toolchain: clang++ (Apple), cmake 3.29.2, Eigen 5.0.1, nlohmann-json 3.12, doctest 2.5.3 in /opt/homebrew/include.
STATE.md:20: - [x] P2 Builder-Core → code/ (≈4k lines C++, 28 doctest cases / 6436 assertions). Gate 2 PASSED: orchestrator rebuilt from scratch (0 warnings) and re-ran tests; reference cases incl. Fig 21 (3,4,3,12), Fig 6 hex, 2025 Fig 2 patterns, r...
STATE.md:25: - [x] P4 Ideator-rigidity → ideas/persona_rigidity.md (678 lines) + ideas/rigidity_rig_check.cpp (its own numeric checks on 21 graphs, compiled against core sources directly).
STATE.md:35: - [~] P4b Ideation round 2 (human directive): ideator2-theorist / ideator2-inverse / ideator2-theorist-fable pending; ideator2-theorist-b DONE → ideas/round2_theorist_b.md (877 lines, B1–B10). B1 hole area is a first harmonic A_C(θ)=a_C ...
STATE.md:36: - [x] experimenter-b4 → B4 FAIL (results/kill/b4/, KILL_REPORT §B4, kill_b4 app, 71 tests/16,488 assertions). Exact Θ_max>0 with boundary rows dropped: 0/200 (σ_mc), 2/200 (σ_def) = 0.5% vs 10% bar; refereed by bisection 1/400 (voronoi_9...
STATE.md:37: - [x] critic-r2 → ideas/ranking_r2.md (539 lines). Top-5: X1 = A1 expansive-cone LP (+A4 hero, kill_k8a: LP over ker A at X_ini on 100 K1a graphs, PASS ≥20/100 positive certified margin; soundness kill if σ ∉ P(X) on the 4 split-bearing ...
STATE.md:38: - [x] checker-t5 → round 6 complete (derivation_tests 33 cases / 89,074 assertions; kiri_tests 16,622 — both re-run by orchestrator). Add-ons: (a) voronoi_93 = REFEREE ARTEFACT (F34); (b) five round-5 corrections applied to core.md body,...
STATE.md:40: - [x] experimenter-t1 → kill_t1 FAIL (results/kill/t1/, KILL_REPORT §T-1, 7 s off K6 cache). 400 designs, 165,788 interior vertices, 14,886 pure (no split edge); balance (BAL) holds EXACTLY at pure vertices (max |δ_v| 2e-14) yet 1,968 of...
STATE.md:45: - [x] experimenter-jitter → A3 DONE (results/kill/jitter/, KILL_REPORT §A3, kill_jitter app, 65 tests/16,229 assertions). MIXED. Exact Θ_max transition exists and is sharp (0.36–0.61 decades): a* = 0.16 (snub), 0.28 (hex), 0.32 (4.8.8), ...
STATE.md:50: - [x] P7 Checker → derivations/check.md + code/tests/derivation_tests.cpp (standalone, 25 cases / 29,033 assertions pass). VERIFIED: T1 (+Step 6 robust), T2, T3, T4.2/2′/2″ (graze reproduced), T4.4 (with min(minβ,π)), T4.5a (moving centr...
STATE.md:51: - [x] P8 Builder-Export → code/src/export (SVG laser layers cut/score/engrave; living-hinge neck & pin-pad solids; binary STL; 3MF with own zip writer), kiri_export CLI, 12 tests / 1018 assertions; export/samples (snub square, 3.4.3.12: ...
STATE.md:56: - [x] P7 Deriver round 5 + Checker round 5 (fresh, checker-6): ZERO unresolved disagreements: YES. derivation_tests.cpp: 31 cases / 89,055 assertions, 0 failures. Gate 7 PASSED. Three presentational notes left in check.md Round 5 (A writ...
STATE.md:58: - [x] Baseline parity: baseline/native builds the authors' code natively (libigl 2.5, Optiz, Eigen 3.4 via FetchContent; bind.cpp replaced by CLI). Parity 8/8 on edge sets, Eq.(6) residual, FK@30° (≤4e-8). See baseline/parity.md.
STATE.md:67: - D3 (human directive 2026-09-03): ALL coding tasks in C++ (core, method, export, tests, experiment drivers, data generation). Python is allowed only for matplotlib plotting of CSV/JSON dumped by C++. Every Builder/Experimenter spec must...
STATE.md:74: - D9 (orchestrator, 2026-09-04 17:15, resolves ESCALATION.md): constructive half = K9's convexity + split-inward constrained embedding (F36). Round-2 outcome: of 8 kill tests, theorems PASS (B1/F35, K7 C1/C2/C4, certificate round 6), alg...
STATE.md:111: - F38 (checker-k8a, results/kill/k8a/recheck.md): cone_lp's primal witness came from an under-converged log-sum-exp loop while the Frank–Wolfe dual was sound; the FW iterate w = Aᵀλ at the min-norm point IS the primal optimum (a_i·w ≥ ‖w...
STATE.md:112: - F37 FINAL (K9c, 400/400, orchestrator's own aggregation 20:35): deployable 307/400 (76.8%), refereed 307/307, certified 306, ε_max ≥ 0.1 on 214, median 0.278, max π; σ_mc 155/200, σ_def 152/200; Delaunay 108/134, Voronoi 122/134, quad ...
STATE.md:113: - F36 (K9, convexity-constrained embedding, results/kill/k9/, K6's 400 designs; agent killed by session limit before writing §K9): FIRST NON-ZERO CONSTRUCTIVE RESULT ON RANDOM GRAPHS. Variant (b) = convexity margin + split-inward q_e ≥ δ...
STATE.md:115: - F34 (checker-t5, B4 voronoi_93): the bisection referee theta_bisect (collision.cpp has_collision/polygons_overlap) is NOT ground truth. On voronoi_93 it reports a collision of faces (94,184) from θ=1e-7 that a 1200×1200 grid and 17,280...
STATE.md:117: - F32 (A3 → checker-cert, FIXED): the certificate's NOROOT clause omitted the two interval tests (0 ≤ ⟨w−a,b−a⟩ ≤ |b−a|²) that the exact scan applies, so collinearity with an edge's INFINITE LINE counted as a contact (deliberate weakenin...
STATE.md:123: - U10 RESOLVED (builder-method-2): the code is deterministic (no clocks/threads/unordered containers; Eigen single-threaded). The id-130 replay differed only because it used range_embed's 8×200 stage defaults instead of the run's --stage...
STATE.md:163: - [~] WP1 Deriver-L DONE → derivations/lemmas.md (869 lines; L1.1–L1.4, L2.1–L2.3 proved; converse of L1.1 under H-L1 which is measurably false on 0.58 % of pairs; L2.2(e) rank(D)=min(2,k) CONJECTURE 132/133; checks: 84 graphs / 30.1M ca...
STATE.md:164: - [x] WP2 Experimenter-Y → results/yield/{YIELD.md, features.csv 400, fresh.csv 200, analyse.py, 3 figs}, code/apps/kill_yield.cpp. FAIL vs bar (F43); mechanism measured. Solver-basin hypothesis → WP2b.
STATE.md:165: - [x] WP3 CANCELLED (D15, human directive: no fabrication, theory only). Builder-fab stopped, export/fab and its apps removed, kiri_export.cpp reverted.
STATE.md:167: - [x] WP6b Propagator-A DONE: derivation_tests CMake target (code/CMakeLists.txt:68); test_reference_cases.cpp 204 lines, kiri_tests 109/17,854 → 111/17,975 from its change; build duplicates + orphan dbg_k1c removed; 12 text items (REPOR...
STATE.md:173: - E1 (environment, 2026-09-08): `cmake --build` / `make` / `lldb` fail (xcode-select shim error 34304). Use `/Applications/Xcode.app/Contents/Developer/usr/bin/make -C code/build <target>` directly. Load average 17–27 with 8 native proce...
STATE.md:184: - Working tree clean at a18c000. baseline/native/build/tuttekiri_cli present. code/build/{kiri_tests,kill_native200} present. em++ at /opt/homebrew/bin. Bambu Studio installed; PrusaSlicer not. Load average 1.6 (machine idle). python3 ma...
STATE.md:192: - D15 (human directive 2026-09-08 21:40): NO fabrication. Mission 2 stays theory/computation only. WP3 removed from MISSION2.md; builder-fab stopped; export/fab/, kiri_export_tiles/kiri_predict/kiri_fab_scale deleted; kiri_export.cpp and...
STATE.md:201: - 2026-09-08 22:22: load spike to 84 observed (regime native shards + basin + rerun); by 22:22 heavy processes ≈ 18 cores (17 tuttekiri_cli, 3 kill_basin, 1 kill_scaling, 1 clang), load decaying. No throttling. Caveat for WP7b: scaling t...
STATE.md:203: - 2026-09-08 22:44: USER STOPPED deriver-l (round 2 in progress, 5 R2 corrections already in lemmas.md), checker-l (idle), exp-regime (7a done: results/regime/REGIME.md; 7b partial: results/scaling/SCALING.md, scaling.csv, capped.csv), e...
```

### `IDEA.md` (2)

_two file path refs: `code/src/method/range_embed.{hpp,cpp}` (l.171), `code/tests/test_design.cpp` (l.214)._

```
IDEA.md:171: (`code/src/method/range_embed.{hpp,cpp}`). Same population, same shape space
IDEA.md:214: regression in `code/tests/test_design.cpp`). Exported closed, at `Θ_max/2` and at `0.9 Θ_max`
```

### `derivations/check.md` (46)

_clang++ build lines for `derivation_tests` (x2); doctest test counts (25 cases / 29,033 assertions); many `contact.hpp/.cpp`, `kinematics.cpp`, `deploy_basis.hpp`, `collision.cpp` refs; scratch `check_*.cpp` programs (Programs-run table)._

```
derivations/check.md:5: re-derived by hand from the code's own sign convention (`code/src/core/kinematics.cpp`), and
derivations/check.md:7: `code/tests/derivation_tests.cpp`, which recomputes the face potential `u`, the harmonic
derivations/check.md:9: `Γ`-cycle closure rows **without using the Deriver's scratch programs** (`check_t1_t2.cpp`,
derivations/check.md:10: `check_t4_t5.cpp` were never compiled or run by me).
derivations/check.md:15: clang++ -std=c++20 -O2 -arch arm64 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
derivations/check.md:16:   -Icode/src code/tests/derivation_tests.cpp code/build/libkiri_core.a \
derivations/check.md:17:   -o code/build/derivation_tests && ./code/build/derivation_tests
derivations/check.md:20: `code/CMakeLists.txt` and every existing file under `code/` are untouched. Result:
derivations/check.md:65: | 0.2, 0.4–0.6 | AGREE | definitions; match `cut.hpp` / `mesh.hpp` verbatim |
derivations/check.md:66: | 0.7, 0.7′ | AGREE | `kinematics.cpp` sets `a = -sigma[f]*theta*0.5`; `R(a) = cI + sJ` expands to `cI − σ_f s J`. Confirmed numerically: the closed form reproduces `deploy()` to 1.4e−14, the sign-flipped form fails by 45 (test **T1-a**) |
derivations/check.md:93: | T2.1 (T2.2) | AGREE | the pin constraint and the circulation formulation are the standard body-and-pin model; matches `mobility.hpp::build_A` |
derivations/check.md:117: | Thm T4.2″ (T4.1) | **AGREE** | I checked the one gap the statement could have: an *isolated* overlap angle. Interior overlap is an **open** condition in `θ` (an open ball inside both interiors persists under a continuous motion), so th...
derivations/check.md:166: **in the face's own frame**", and `contact.hpp::swept_discs` implements exactly that, with
derivations/check.md:191: 1. `contact.hpp`'s `rho` (the `√2` form, commented "sound bound") is the one that is loose by
derivations/check.md:258: `contact.hpp` carries a separate `penetrates_immediately()`, and exactly the failure mode
derivations/check.md:322: ### D8 — `contact.hpp`'s flat-centroid pruning is unsound in principle (the Deriver's (T4.3) is not)
derivations/check.md:325: sound. `contact.hpp::candidate_pairs(..., use_static = true)` instead uses the **flat**
derivations/check.md:333: ### D9 — documentation inconsistency in `deploy_basis.hpp`
derivations/check.md:374: All from `code/tests/derivation_tests.cpp`; 25 test cases, 29 033 assertions, 0 failures.
derivations/check.md:494: this file's round-1 verdicts, `derivations/scratch/check_r2.cpp`, `code/README.md`.
derivations/check.md:495: New program work: two test cases appended to `code/tests/derivation_tests.cpp` (that file only;
derivations/check.md:496: `CMakeLists.txt` and every other source untouched). Rebuilt with the command in its header:
derivations/check.md:499: clang++ -std=c++20 -O2 -arch arm64 -I/opt/homebrew/include \
derivations/check.md:500:   -I/opt/homebrew/include/eigen3 -Icode/src \
derivations/check.md:501:   code/tests/derivation_tests.cpp code/build/libkiri_core.a \
derivations/check.md:502:   -o code/build/derivation_tests && ./code/build/derivation_tests
derivations/check.md:627: ## R2.6 New tests (appended to `code/tests/derivation_tests.cpp`)
derivations/check.md:674: T5.2b′, T5.2b.1, T5.2b.2, T6.4; `derivations/scratch/check_r3.cpp`. New tests appended to
derivations/check.md:675: `code/tests/derivation_tests.cpp` (that file only): `R3 Lemma T5.1e and the three-class atom
derivations/check.md:705: **`check_r3.cpp` rebuilt and rerun.** Its numbers reproduce exactly, at all three `ε`:
derivations/check.md:720: **not blind to class 3** in the way `check_r2.cpp` was. R2.5's circularity objection is answered.
derivations/check.md:748: Both `check_r3.cpp` and my own test exclude them *before* classifying (`h.scale() ≤ tol → skip`,
derivations/check.md:831: Appended to `code/tests/derivation_tests.cpp` (only that file; no other file touched), built with
derivations/check.md:878: `R4 Sub-lemma T5.2b'' Case A ...` appended to `code/tests/derivation_tests.cpp` (that file only;
derivations/check.md:1035: and the zero set of h_o,pi' in (0,pi)` appended to `code/tests/derivation_tests.cpp` (that file only;
derivations/check.md:1164: `code/src/method/contact.cpp` (`validity_certificate`, `contact_angles`) with the F32 doctest in
derivations/check.md:1165: `code/tests/test_method.cpp`. Two new programs: `derivations/scratch/check_t5_adm.cpp` (standalone,
derivations/check.md:1166: build line in its header) and the `R6` test case appended to `code/tests/derivation_tests.cpp`
derivations/check.md:1167: (that file only), plus `derivations/scratch/check_b4_93.cpp` and the `R6-c` case for the B4 add-on.
derivations/check.md:1208:    The segment of T4.1b is **closed**, so equality is admissible; `contact.cpp`'s inclusive
derivations/check.md:1239: to `NOROOT`". `contact.cpp` never implemented them. **The code was right and [D3] was unnecessary:
derivations/check.md:1271: **Counterexample, measured** (`check_t5_adm.cpp` §B, `derivation_tests` R6-b): `hexagons` at `X₀` has
derivations/check.md:1292: Replayed through `kill_b4.cpp`'s own pipeline (same graph, same `sigma_def`, same `free` system, same
derivations/check.md:1293: seeds and weights) by `derivations/scratch/check_b4_93.cpp`. The row reproduces exactly:
derivations/check.md:1339: its centroid, or exclude edge pairs that share an endpoint) is a change to `code/src/core/collision.cpp`
derivations/check.md:1342: **Regression case added.** `derivation_tests` `R6-c` pins the reproduction with the two polygons as
```

### `derivations/check_lemmas.md` (28)

_file path refs (`*.hpp`, `apps/kill_b3.cpp`, `kill_k7.cpp`); clang++ build line + doctest output block (37 cases / 149,115 assertions) for `derivation_tests`; `code/src|apps|tests|CMakeLists` non-modification statements._

```
derivations/check_lemmas.md:10:    and the code (`deploy_basis.hpp`, `contact.hpp`, `zero_plus.hpp`, `periodic_jacobian.hpp`,
derivations/check_lemmas.md:11:    `apps/kill_b3.cpp::achievable/shape_point`, `tests/derivation_tests.cpp`).
derivations/check_lemmas.md:83: identity `zero_plus.hpp` states as `dC_e = 0` and it is here proved for hinge- and vertex-only
derivations/check_lemmas.md:127: `A ≠ 0`). The two interval predicates of `contact.hpp` are only ever evaluated *at* a root, so no
derivations/check_lemmas.md:181:   `q_e = 0` degenerate case `zero_plus.hpp` is written to detect — so it is a real, describable
derivations/check_lemmas.md:218: This also matches the code: `kill_b3.cpp::achievable` reads `D(i,0..1) = ½ (M_{2i} T)(row 2)`,
derivations/check_lemmas.md:383: | 1 (`dC = 0`, `Y_w − Y_{a*} = s·dS`) | **AGREE** | matches `zero_plus.hpp`'s `dC_e = 0`, correctly generalised |
derivations/check_lemmas.md:433:    `zero_plus.hpp`/K5/F30 report as *failing* on measured designs. So the hypothesis is not free and
derivations/check_lemmas.md:554: `achievable()` in `kill_k7.cpp` / `kill_b3.cpp` **never reads `consistency`**. And note why the
derivations/check_lemmas.md:614: Four `TEST_CASE`s appended to `code/tests/derivation_tests.cpp` — that file only. Nothing under
derivations/check_lemmas.md:615: `code/src`, `code/apps`, `code/CMakeLists.txt`, `derivations/lemmas.md` or `derivations/core.md`
derivations/check_lemmas.md:616: was touched. The file stays standalone doctest.
derivations/check_lemmas.md:623: Two additions to the file's includes: `method/periodic_jacobian.hpp` and `apps/kill_common.hpp`
derivations/check_lemmas.md:626: and `achievable()` are copied from `code/apps/kill_k7.cpp` / `kill_b3.cpp` so the numbers compare
derivations/check_lemmas.md:632: The build line in `derivations/check.md`'s header, plus **`-Icode/apps`** (the only change), and
derivations/check_lemmas.md:633: with the Xcode toolchain compiler by its full path, because a bare `clang++` aborts on this machine
derivations/check_lemmas.md:637: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang++ \
derivations/check_lemmas.md:640:   -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 -Icode/src -Icode/apps \
derivations/check_lemmas.md:641:   code/tests/derivation_tests.cpp code/build/libkiri_core.a \
derivations/check_lemmas.md:642:   -o code/build/derivation_tests && ./code/build/derivation_tests
derivations/check_lemmas.md:648: [doctest] test cases:     37 |     37 passed | 0 failed | 0 skipped
derivations/check_lemmas.md:649: [doctest] assertions: 149115 | 149115 passed | 0 failed |
derivations/check_lemmas.md:650: [doctest] Status: SUCCESS!
derivations/check_lemmas.md:694: `kill_common.hpp::reference_cases()` plus `make_graph(id, 18, 46, 220)` for `id = 0 …` until 60 are
derivations/check_lemmas.md:756: **Corpus.** The **33 K7 patterns** (`l2_k7_population()`, byte-identical to `kill_k7.cpp`'s
derivations/check_lemmas.md:802: exactly against `core.md` (T5.2)/(T5.4) for split edges and against `zero_plus.hpp`'s `dS_e`. The
derivations/check_lemmas.md:817: `zero_plus.hpp`/K5/F30 report as the *observed* failure mode. Relatedly, the box's "decided by
derivations/check_lemmas.md:821: is what `contact.cpp` already does. So the engineering claim is safe and the theorem needs its
```

### `derivations/core.md` (51)

_file refs (`kinematics.cpp`, `deploy_basis.hpp`, `contact.hpp`, `collision.hpp`, `cut.hpp`); test counts (25/29,033; 30/89,044); scratch `check_*.cpp` programs and their clang++ build lines (Programs-run table l.1984–1991, 2261–2263)._

```
derivations/core.md:23: 29 033 assertions, its own program `code/tests/derivation_tests.cpp`), raising eleven
derivations/core.md:25: `derivations/scratch/check_r2.cpp` (build line in its header; `graphs used: 16`, 197 shape-space
derivations/core.md:37: | **D8** | `contact.hpp`'s flat-centroid pruning is unsound in principle | **ADOPTED as a note.** My (T4.3) is the sound static test; the flat-centroid variant discards 19 842 pairs the sound test keeps (no observed wrong `Θ_max`, 96/96)...
derivations/core.md:38: | **D9** | `u_f` means different objects here and in `deploy_basis.hpp` | **CORRECTED** by a notation note fixing this file's `u_f` and naming the header's object explicitly. No code was touched (other agents own `code/`) | §0.10 |
derivations/core.md:53: Round 3's program is `derivations/scratch/check_r3.cpp` (same 16-graph corpus, same seed, same 197
derivations/core.md:58: | **R2.5** | a **third** structural class `g(0) = g′(0) = 0` (`p = −q`, `r = 0`, `h = p(1 − cos θ)`) is unnamed, and round 2's "0 mismatches" is not reproducible | **CONCEDED, and round 2's number is WITHDRAWN as circular** — `check_r2.c...
derivations/core.md:91: Deciding numbers for round 4: **unchanged** — `check_r3.cpp` still reports 0 mismatches of
derivations/core.md:92: 2 145 387 at every `ε ∈ {0.2, 0.02, 0.001}`, and the standalone `code/tests/derivation_tests.cpp`
derivations/core.md:109: Deciding numbers for round 5: **unchanged**. `check_r3.cpp` still reports 0 mismatches of
derivations/core.md:111: `code/tests/derivation_tests.cpp` re-runs at **30 test cases, 89 044 assertions, 0 failures**.
derivations/core.md:129: (`cut.hpp`: `corner_to_prime`); write `(v, f)` for the copy of `v ∈ V` carried by face `f`, and
derivations/core.md:142: **0.7 [F, `code/src/core/kinematics.cpp`] The code's sign convention.** Face `f` is transformed by
derivations/core.md:171: **0.10 [D] The symbol `u_f`, and the clash with `deploy_basis.hpp` (D9).** In this file `u_f ∈ R²`
derivations/core.md:173: satisfies the increment law (T1.3). The header comment of `code/src/method/deploy_basis.hpp` writes
derivations/core.md:185: `deploy_basis.hpp` should rename its comment variable to `ν_f` before either symbol becomes a
derivations/core.md:272: **[N] Path-independence check.** `derivations/scratch/check_t1_t2.cpp`, check **C1**: build `u` by BFS
derivations/core.md:307: BFS, against `kinematics.cpp::deploy()` at `θ ∈ {0, 0.05, 0.3, 0.77, 1.2, 1.9, 2.6, 3.0}` on the same
derivations/core.md:390: Re-run `derivations/scratch/check_t1_t2.cpp`. Tolerance `1e−12` on C1–C6, `1e−9` on C7. Additionally:
derivations/core.md:497: round-off). Independent of `check_t1_t2.cpp`'s C6, which tests the pin equations rather than the
derivations/core.md:566: matching `deploy_basis.hpp::harmonic_roots`. Each of `p, q, r` is quadratic in `t`, so the
derivations/core.md:629: quantity `collision.hpp::theta_max` estimates by grid scan + bisection with faces shrunk by a
derivations/core.md:687: **[N] Measured (`derivations/scratch/check_t4_t5.cpp`, checks D1/D2).** 16 valid embeddings
derivations/core.md:782: **[N] Measured (`derivations/scratch/check_r2.cpp`, R2-C).** On all **8 split-free** patterns
derivations/core.md:822: faces actually approach can be discarded. `code/src/method/contact.hpp::candidate_pairs(...,
derivations/core.md:841: which no pruning uses. K2c's spec and `contact.hpp::swept_discs` both work **in the face's own
derivations/core.md:862: **[N] Measured (`derivations/scratch/check_r2.cpp`, R2-A).** Over **1 628** `M′`-copies on 16
derivations/core.md:870: 1. Of the two radii in `contact.hpp`, `rho` (the `√(‖x‖²+‖χ‖²)` form, commented "sound bound") is
derivations/core.md:933: * Re-run `derivations/scratch/check_t4_t5.cpp` (rules: D1 `≤ 1e−5`, D2 `≤ 1e−12`) **and**
derivations/core.md:934:   `derivations/scratch/check_r2.cpp` (rules: R2-A `≤ 1e−12`, R2-C `≤ 1e−12`, R2-D violations `= 0`;
derivations/core.md:935:   its R2-E line is **withdrawn**, see T5.2b.2) **and** `derivations/scratch/check_r3.cpp`
derivations/core.md:1082: the negation of `contact.hpp::penetrates_immediately()`. `EMB` is what the proof needs, but it is a
derivations/core.md:1501: for it. Both `derivations/scratch/check_r3.cpp` and `code/tests/derivation_tests.cpp` implement
derivations/core.md:1523: **[N] Measured, round 3 (`derivations/scratch/check_r3.cpp`).** Round 2's number
derivations/core.md:1524: ("**0 / 2 363 380** at every `ε`") was **not a valid measurement and is withdrawn**: `check_r2.cpp`
derivations/core.md:1526: comparison could not see class 3 at all. `check_r3.cpp` fixes this by computing the truth **without
derivations/core.md:1706: * **Round 3:** re-run `derivations/scratch/check_r3.cpp` with several `ε`. The **three-class**
derivations/core.md:1711:   Re-run `derivations/scratch/check_r2.cpp` for R2-A/R2-C/R2-D; R2-D must report **0** violations;
derivations/core.md:1793: 5. **Independence** — the same `Θ_max` cross-checked against `collision.hpp`'s bisection, which uses
derivations/core.md:1984: | `derivations/scratch/check_t1_t2.cpp` | `clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 -I code/src …` linking `code/src/core/*.cpp` | C1 1.60e−14, C2 1.42e−14, C3 9.73e−14, C4 2.08e−14, C5 4....
derivations/core.md:1985: | `derivations/scratch/check_t4_t5.cpp` | same | D1 2.09e−6 (rule 1e−5), D2 8.88e−16, D3 4.37e−15 over 72 split edges, D4 719/1061 violations in the **raw** frame, ratio 1.4141 — see T4.5b: this says nothing about any pruning |
derivations/core.md:1986: | `derivations/scratch/check_r2.cpp` (round 2) | same, plus `code/src/core/*.cpp`; build line in the file header, takes `ε` as `argv[1]` | R2-A 3.11e−15 (1 628 copies), R2-B drift/`r_f` 4.62→20.03 over diameters 5.66→19.80, R2-C 8.88e−16...
derivations/core.md:1987: | `derivations/scratch/check_r3.cpp` (round 3) | same corpus and seed; build line in the file header, takes `ε` as `argv[1]` | three-class deflated atom list vs a `τ`-chart-free crossing test: **0** mismatches of 2 145 387 at `ε ∈ {0.2, ...
derivations/core.md:1991: `code/tests/derivation_tests.cpp` (25 cases, 29 033 assertions, 0 failures), was written
derivations/core.md:2037: exactly at the far endpoint, `derivation_tests` R6-a). `contact.cpp` implements `E2` inclusively at
derivations/core.md:2075: > Measured (`derivation_tests` R6-a1/R6-a2, 9 642 Case-A substitute pairs over the corpus with
derivations/core.md:2148: not as separate atoms whose own roots are hunted. This is what `contact.cpp` implements, so the
derivations/core.md:2187: **[N] Measured counterexample (`derivations/scratch/check_t5_adm.cpp`, and `derivation_tests` R6-b).**
derivations/core.md:2231: `derivations/scratch/check_b4_93.cpp`: **the scan is right.**
derivations/core.md:2261: | `derivations/scratch/check_t5_adm.cpp` (new; build line in its header) | A: 1 363 Case-A pairs on the 7 tilings with jitter, ratio `∈ [0.5404, 1.0]`, `|s/‖e‖² − L′/L| ≤ 1.55e−15`, **0** admissibility failures. B: 1 completeness counter...
derivations/core.md:2262: | `derivations/scratch/check_b4_93.cpp` (new; build line in its header) | replays B4 ids 93 and 96 from the K6/B4 caches; 93 reproduces scan `0.248400` vs referee `0`, and localizes the referee's misfire to the shared hinge vertex of fac...
derivations/core.md:2263: | `code/tests/derivation_tests.cpp` R6 + R6-c cases (new) | R6-a1 **0** failures / 9 642 · R6-a2 `1.77e−14` · R6-a3 `3.53e−14` · R6-b counterexample found · R6-c1 **0** common interior points of 1 442 401 · R6-c2 misfire pinned. Whole fi...
```

### `derivations/lemmas.md` (42)

_full clang++ build lines listing `code/src/core/*.cpp` (l.461–465, 788–790, 930–940); Xcode toolchain workaround note; `zero_plus.hpp`/`cut.hpp`/`contact.cpp`/`periodic_jacobian.hpp`/`kill_k7.cpp` refs; `code/tests were not run` remark._

```
derivations/lemmas.md:4: `code/src` or `code/tests`. Notation is `core.md` §0 throughout (`c = cos(θ/2)`, `s = sin(θ/2)`,
derivations/lemmas.md:41:   `zero_plus.hpp`; `[F, K5 / F30]` reports non-opening split cuts as the observed failure mode, so
derivations/lemmas.md:55: face `g ≠ f`. Write `ρ : V′ → V` for the source map (`cut.hpp::prime_to_original`) and
derivations/lemmas.md:148: count of distinct pairs; `check_l1.cpp` increments once per sample. The Checker's independent run
derivations/lemmas.md:253:    is the same computation as `zero_plus.hpp`'s header (`dC_e = 0`, `dS_e`) generalised from split
derivations/lemmas.md:353: **H-L3 is not free.** Its first clause is exactly the predicate `zero_plus.hpp` is written to
derivations/lemmas.md:399:    copies of the same `M`-vertex — the hinge point is welded by `cut.hpp::corner_to_prime` into a
derivations/lemmas.md:402:    `contact.cpp`'s `if (p == a || p == b) continue;` already does, *before* the numeric identity
derivations/lemmas.md:449: know that the answer is structural.** That is what `contact.cpp` already does.
derivations/lemmas.md:456: **Program** `derivations/scratch/check_l1.cpp`. **Build line** (in the file header; the repo's
derivations/lemmas.md:457: `clang++` needs `DEVELOPER_DIR` pointed at Xcode on this machine, see the note at the end of this
derivations/lemmas.md:461: clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
derivations/lemmas.md:462:   -I code/src -I code/apps derivations/scratch/check_l1.cpp \
derivations/lemmas.md:463:   code/src/core/mesh.cpp code/src/core/cut.cpp code/src/core/holes.cpp \
derivations/lemmas.md:464:   code/src/core/kinematics.cpp code/src/core/orientation.cpp code/src/core/generators.cpp \
derivations/lemmas.md:465:   code/src/core/tutte_auxetic.cpp code/src/core/collision.cpp code/src/core/optimize.cpp \
derivations/lemmas.md:469: **Corpus.** The 8 Phase-2 reference cases (`kill_common.hpp::reference_cases()`) plus
derivations/lemmas.md:545: `xcodebuild` fails to load, so a bare `clang++` aborts. Both scratch programs were built with
derivations/lemmas.md:547: `.../Toolchains/XcodeDefault.xctoolchain/usr/bin/clang++` with
derivations/lemmas.md:558: quotient system of `periodic_jacobian.hpp` supplies `X(t) = X₀ + Φ t` with `Φ ∈ R^{n_q × k}`,
derivations/lemmas.md:561: orders `t` as `(t_{2i}, t_{2i+1}) = ` (`φ_i` on `x`, `φ_i` on `y`) — `kill_k7.cpp::shape_point`.
derivations/lemmas.md:576: independent of `f`. ∎ This is the header derivation of `periodic_jacobian.hpp`, written out.
derivations/lemmas.md:775: *empirically*; and the corresponding one-line note in `periodic_jacobian.hpp`'s header, which states
derivations/lemmas.md:781: **Program** `derivations/scratch/check_l2.cpp`. Its population, quotient/super-patch pipeline and
derivations/lemmas.md:782: `achievable()` are copied **verbatim** from `code/apps/kill_k7.cpp` (lines 36–243) so the numbers
derivations/lemmas.md:788: clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
derivations/lemmas.md:789:   -I code/src -I code/apps derivations/scratch/check_l2.cpp \
derivations/lemmas.md:790:   code/src/core/*.cpp code/src/method/*.cpp -o /tmp/check_l2 && /tmp/check_l2 100
derivations/lemmas.md:880: **How to attack it.** I read `code/src/method/contact.cpp`'s `validity_certificate` and its
derivations/lemmas.md:900: `kill_k7.cpp` — which reads it off the **even** generators only (`if (j % 2 == 0)`) — is not the
derivations/lemmas.md:904: `M_{2i+1}` against row `i` of `D`, which `kill_k7.cpp` never touches; Check L2 reports
derivations/lemmas.md:930: clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
derivations/lemmas.md:931:   -I code/src -I code/apps derivations/scratch/check_l1.cpp \
derivations/lemmas.md:932:   code/src/core/mesh.cpp code/src/core/cut.cpp code/src/core/holes.cpp \
derivations/lemmas.md:933:   code/src/core/kinematics.cpp code/src/core/orientation.cpp code/src/core/generators.cpp \
derivations/lemmas.md:934:   code/src/core/tutte_auxetic.cpp code/src/core/collision.cpp code/src/core/optimize.cpp \
derivations/lemmas.md:938: clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
derivations/lemmas.md:939:   -I code/src -I code/apps derivations/scratch/check_l2.cpp \
derivations/lemmas.md:940:   code/src/core/*.cpp code/src/method/*.cpp -o /tmp/check_l2 && /tmp/check_l2 100
derivations/lemmas.md:943: On this machine `xcode-select` points at an Xcode whose `xcodebuild` aborts, so a bare `clang++`
derivations/lemmas.md:946: `/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang++`
derivations/lemmas.md:948: defect, not a property of the code; `code/tests` were not run because they go through the same
```

### `derivations/lemmas_statements.md` (4)

_file refs only (`cut.hpp`, `contact.cpp`, `periodic_jacobian.hpp`, `kill_k7.cpp::shape_point`)._

```
derivations/lemmas_statements.md:39: face `g ≠ f`. Write `ρ : V′ → V` for the source map (`cut.hpp::prime_to_original`) and
derivations/lemmas_statements.md:189: know that the answer is structural.** That is what `contact.cpp` already does.
derivations/lemmas_statements.md:195: quotient system of `periodic_jacobian.hpp` supplies `X(t) = X₀ + Φ t` with `Φ ∈ R^{n_q × k}`,
derivations/lemmas_statements.md:198: orders `t` as `(t_{2i}, t_{2i+1}) = ` (`φ_i` on `x`, `φ_i` on `y`) — `kill_k7.cpp::shape_point`.
```

### `results/core_validation/rank_claim.md` (4)

_`kiri_tests` and `tests/test_holes.cpp`, `tests/test_rank_checks.cpp` refs._

```
results/core_validation/rank_claim.md:49: instances by `kiri_tests` (see `tests/test_holes.cpp`).
results/core_validation/rank_claim.md:120: patches: the three torus rows above and `tests/test_rank_checks.cpp`, where `dim Z == 1`
results/core_validation/rank_claim.md:140: notches -- the torus patches in `tests/test_rank_checks.cpp` -- it collapses to the
results/core_validation/rank_claim.md:157: `tests/test_rank_checks.cpp` on the 4x4 and 6x4 square tori and the 4x4 triangle torus.
```

### `results/core_validation/referee_fix.md` (11)

_`code/src/core/collision.{hpp,cpp}` scope; `code/tests/test_collision.cpp` failing-first test; test-count table `kiri_tests` 103 / 17,504 and `derivation_tests` 33 / 89,074 (l.98–99); `check_b4_93.cpp`._

```
results/core_validation/referee_fix.md:3: Scope: `code/src/core/collision.{hpp,cpp}` -- `polygons_overlap` / `has_collision`, the
results/core_validation/referee_fix.md:74: New, in `code/tests/test_collision.cpp` (written first, failing against the old predicate):
results/core_validation/referee_fix.md:89: `derivation_tests.cpp` R6-c, the pinned reproduction from checker round 6, is INVERTED:
results/core_validation/referee_fix.md:93: Suites (both re-run after the fix, from the repo root; the `kiri_tests` totals include test
results/core_validation/referee_fix.md:98: | `kiri_tests` | 103 | 17 504 | 0 failures |
results/core_validation/referee_fix.md:99: | `derivation_tests` | 33 | 89 074 | 0 failures |
results/core_validation/referee_fix.md:104: tree with only `code/src/core/collision.{hpp,cpp}` reverted to `HEAD` (built out of tree in
results/core_validation/referee_fix.md:106: agents made to `contact.cpp` / `zero_plus.cpp` in the same working tree. The baseline K5 run
results/core_validation/referee_fix.md:128: ### B4 -- the F34 case itself (`derivations/scratch/check_b4_93.cpp`, ids 93 and 96)
results/core_validation/referee_fix.md:171: probe, `contact.cpp`), not only in the referee: the same hinge artefact was rejecting the
results/core_validation/referee_fix.md:178: `method/zero_plus.cpp` in the same working tree, not this fix.
```

### `results/final/e1/E1.md` (4)

_driver `code/apps/kill_e1.cpp`; `kill_common::` population refs; `cmake --build . --target kill_e1`._

```
results/final/e1/E1.md:11: Code: `code/apps/kill_e1.cpp` (new; does not modify any core/ or method/ source). Three
results/final/e1/E1.md:16: - **authored**: `kill_common::deployable_population()` called with (samples_per_base=20,
results/final/e1/E1.md:21: - **random**: `kill_common::make_graph(id, 100, 800, 1400)` -- Voronoi / Delaunay /
results/final/e1/E1.md:101: cmake --build . --target kill_e1
```

### `results/final/figures/README.md` (1)

_one phrase: 're-run any C++ solver'._

```
results/final/figures/README.md:5: re-run any C++ solver).
```

### `results/kill/KILL_REPORT.md` (114)

_the largest: 'C++ drivers in `code/apps/kill_*.cpp`' (l.4); per-experiment Driver: lines naming `code/apps/kill_*.cpp`/`dbg_*.cpp`; method file refs `code/src/method/*.{hpp,cpp}`; 'new doctest in `code/tests/test_method.cpp`' + running suite counts (53/16,004; 65/16,229; 66/16,283; 71/16,488; 80/16,785); four cmake + `./code/build/kill_*` reproduce blocks (l.986–997, 1164–1171, 1422–1440, 1749–1751); `baseline/native/src/main.cpp --dumpdir` (KEEP, baseline)._

```
results/kill/KILL_REPORT.md:4: C++ drivers in `code/apps/kill_*.cpp` against `code/src/core` and `code/src/method`
results/kill/KILL_REPORT.md:7: (`code/apps/kill_common.hpp`), so every experiment sees the same graphs.
results/kill/KILL_REPORT.md:18: `code/src/method/` and covered by a doctest in `code/tests/test_method.cpp`.
results/kill/KILL_REPORT.md:89:    `0.9467`; a `1e−4` fine scan (`code/apps/dbg_k2a_out.cpp`) confirms the closed form. All
results/kill/KILL_REPORT.md:149: off by **exactly ±1**, and every one has `F ≥ 879` — above `mobility.hpp`'s
results/kill/KILL_REPORT.md:151: `code/apps/kill_k3a_recheck.cpp` regenerates each violator with K3a's exact configuration
results/kill/KILL_REPORT.md:199: `recheck_summary.txt`. Driver: `code/apps/kill_k3a.cpp`, `code/apps/kill_k3a_recheck.cpp`.
results/kill/KILL_REPORT.md:231: Artifacts: `results/kill/k1b/{k1b.csv,summary.txt}`. Driver: `code/apps/kill_k1b.cpp`.
results/kill/KILL_REPORT.md:274: (`deployable_population()` in `kill_common.hpp`), which is exactly the population K2a and
results/kill/KILL_REPORT.md:278: Artifacts: `results/kill/k1a/{k1a.csv,summary.txt}`. Driver: `code/apps/kill_k1a.cpp`.
results/kill/KILL_REPORT.md:344: `code/apps/kill_k2c.cpp`.
results/kill/KILL_REPORT.md:353: that root; (ii) `theta_max` from `collision.hpp` (grid scan + bisection, faces shrunk by
results/kill/KILL_REPORT.md:399: overlap probes. Reproducing the Deriver's table exactly (`code/apps/dbg_t422.cpp`):
results/kill/KILL_REPORT.md:431: That was F32; with the interval tests restored (`contact.cpp`, §A3) the re-run gives
results/kill/KILL_REPORT.md:440: Artifacts: `results/kill/k2a/{k2a.csv,summary.txt}`. Drivers: `code/apps/kill_k2a.cpp`,
results/kill/KILL_REPORT.md:441: `code/apps/dbg_t422.cpp`, `code/apps/dbg_k2a_out.cpp`.
results/kill/KILL_REPORT.md:451: `collision.hpp`'s bisection to `<= 1e-6` rad. PASS if the fraction of gamma-ladder-certified
results/kill/KILL_REPORT.md:520: Artifacts: `results/kill/k1c/{k1c.csv,summary.txt}`. Driver: `code/apps/kill_k1c.cpp`.
results/kill/KILL_REPORT.md:584: after the fix, `1e−12` on 36 and `1e−10` on the other 89 — `results/core_validation/referee_fix.md`), and a direct scan (`code/apps/dbg_k5.cpp`) confirms overlap at
results/kill/KILL_REPORT.md:616: `code/apps/kill_k5.cpp`, diagnostic `code/apps/dbg_k5.cpp`.
results/kill/KILL_REPORT.md:625: (their own `theta_max` routine is buggy, `STATE.md` F24). Keep our `collision.hpp` Eq.(9) as
results/kill/KILL_REPORT.md:637:   dumped (a `--dumpdir` option added to our CLI in `baseline/native/src/main.cpp`; nothing in
results/kill/KILL_REPORT.md:710: `code/apps/kill_k2b.cpp`; CLI change in `baseline/native/src/main.cpp` (`--dumpdir`,
results/kill/KILL_REPORT.md:726: **Measured** (`code/apps/kill_f23.cpp`), stronger than asked: **every** seed face, not two,
results/kill/KILL_REPORT.md:759: Artifacts: `results/kill/f23/{f23.csv,summary.txt}`. Driver: `code/apps/kill_f23.cpp`.
results/kill/KILL_REPORT.md:839: quoted** until `mobility.hpp`'s `matrix_rank` is fixed or the dense path is forced; a dense
results/kill/KILL_REPORT.md:841: `code/apps/kill_k3a_recheck.cpp` provides (7 of 7 violators small enough for a dense QR have
results/kill/KILL_REPORT.md:858: **The 0⁺ calculus** (`code/src/method/zero_plus.{hpp,cpp}`, derived there in full). By T1.B
results/kill/KILL_REPORT.md:952: The fix is the three-line guard now at `code/src/method/contact.cpp:339` — keep a root only
results/kill/KILL_REPORT.md:978: native/, smoke/}`. Driver `code/apps/kill_k6.cpp`, method `code/src/method/zero_plus.{hpp,cpp}`,
results/kill/KILL_REPORT.md:979: diagnostic `code/apps/dbg_k6.cpp`.
results/kill/KILL_REPORT.md:986: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
results/kill/KILL_REPORT.md:987: ./code/build/kiri_tests                    # 53 cases, 16004 assertions, all passing
results/kill/KILL_REPORT.md:988: ./code/build/kill_k3a  --out results/kill/k3a          # ~25 min
results/kill/KILL_REPORT.md:989: ./code/build/kill_k3a_recheck --out results/kill/k3a   # ~4 min
results/kill/KILL_REPORT.md:990: ./code/build/kill_k1b  --out results/kill/k1b          # ~34 s
results/kill/KILL_REPORT.md:991: ./code/build/kill_k1a  --out results/kill/k1a          # ~40 s
results/kill/KILL_REPORT.md:992: ./code/build/kill_k5   --out results/kill/k5 --n 200 --maxf 800   # ~5 min
results/kill/KILL_REPORT.md:993: ./code/build/kill_k2c  --out results/kill/k2c          # ~3 min
results/kill/KILL_REPORT.md:994: ./code/build/kill_k2a  --out results/kill/k2a          # ~7 s
results/kill/KILL_REPORT.md:995: ./code/build/kill_k1c  --out results/kill/k1c          # ~31 s   (needs baseline/native)
results/kill/KILL_REPORT.md:996: ./code/build/kill_k2b  --out results/kill/k2b --maxf 160   # ~26 min (needs baseline/native)
results/kill/KILL_REPORT.md:997: ./code/build/kill_f23  --out results/kill/f23          # ~1 s
results/kill/KILL_REPORT.md:1004: Method code changed for the corrections: `code/src/method/deploy_basis.{hpp,cpp}`
results/kill/KILL_REPORT.md:1006: `code/src/method/contact.{hpp,cpp}` (`exact_theta_max_overlap`, `validity_certificate`, the
results/kill/KILL_REPORT.md:1007: exact face-frame radius in the broad phase). New doctest coverage in
results/kill/KILL_REPORT.md:1008: `code/tests/test_method.cpp`: the hexagon graze; `Θ_max` vs bisection and pruning
results/kill/KILL_REPORT.md:1036: * **Driver** `code/apps/kill_jitter.cpp`. No `core/` or `method/` source was modified.
results/kill/KILL_REPORT.md:1053:   (`--analyze`) is also C++; python draws the figure only.
results/kill/KILL_REPORT.md:1069: tolerate jitter but because **they never left the start point**. New doctest in
results/kill/KILL_REPORT.md:1070: `code/tests/test_method.cpp` ("split-free tiling: Eq. (6) projection undoes an arbitrary
results/kill/KILL_REPORT.md:1164: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
results/kill/KILL_REPORT.md:1165: ./code/build/kiri_tests                              # 65 cases, 16229 assertions, all passing
results/kill/KILL_REPORT.md:1167: seq 0 11 | xargs -P 12 -I{} ./code/build/kill_jitter --shard {} --nshards 12 \
results/kill/KILL_REPORT.md:1171: ./code/build/kill_jitter --analyze results/kill/jitter/jitter.csv --out results/kill/jitter
results/kill/KILL_REPORT.md:1209: * **Driver** `code/apps/kill_k7.cpp`, stages `main`, `c3`, `c4bounded`, `nu`, `c4sweep`,
results/kill/KILL_REPORT.md:1211:   patch, `periodic_jacobian`, the achievable set — is `code/src/method/periodic_jacobian.{hpp,cpp}`,
results/kill/KILL_REPORT.md:1410: **New doctest.** `code/tests/test_method.cpp`, "C4: theta_c = 2 atan2(tr K, 1 − det K) is where
results/kill/KILL_REPORT.md:1422: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
results/kill/KILL_REPORT.md:1423: ./code/build/kiri_tests                              # 66 cases, 16283 assertions, all passing
results/kill/KILL_REPORT.md:1424: ./code/build/kiri_tests -tc="C4: theta_c*"           # the K7 doctest alone: 24 assertions
results/kill/KILL_REPORT.md:1428: for i in $(seq 0 11); do ./code/build/kill_k7 --stage main --out results/kill/k7 \
results/kill/KILL_REPORT.md:1432: for i in $(seq 0 11); do ./code/build/kill_k7 --stage c3 --out results/kill/k7 \
results/kill/KILL_REPORT.md:1435: for i in $(seq 0 11); do ./code/build/kill_k7 --stage c3 --out results/kill/k7/c3_v2 \
results/kill/KILL_REPORT.md:1437: ./code/build/kill_k7 --stage c4bounded --out results/kill/k7
results/kill/KILL_REPORT.md:1438: ./code/build/kill_k7 --stage nu        --out results/kill/k7
results/kill/KILL_REPORT.md:1439: ./code/build/kill_k7 --stage c4sweep   --out results/kill/k7   # the |det P| artefact
results/kill/KILL_REPORT.md:1440: for i in $(seq 0 11); do ./code/build/kill_k7 --stage perdbg --out results/kill/k7 \
results/kill/KILL_REPORT.md:1478: `kill_common::make_graph`, the same `σ_mc` from Eq. (1) and `σ_def` read back from
results/kill/KILL_REPORT.md:1581: No core or method source was modified for this run; `code/apps/kill_b4.cpp` builds its
results/kill/KILL_REPORT.md:1586: `code/apps/kill_b4.cpp`.
results/kill/KILL_REPORT.md:1589: ./code/build/kill_b4 --n 200 --nshards 12 --shard I --out results/kill/b4/shards \
results/kill/KILL_REPORT.md:1599: Driver `code/apps/kill_k8a.cpp`, method `code/src/method/expansive_cone.{hpp,cpp}`,
results/kill/KILL_REPORT.md:1600: artifacts `results/kill/k8a/` (`k8a.csv`, `summary.txt`, `log_*.txt`), 5 new doctests in
results/kill/KILL_REPORT.md:1601: `code/tests/test_method.cpp` (suite now 71 cases / 16 488 assertions, all passing).
results/kill/KILL_REPORT.md:1607: which is exactly `mobility.hpp::build_rigidity`. Because every copy of a source vertex
results/kill/KILL_REPORT.md:1609: `zero_plus.hpp` becomes a **strict linear inequality in the flex**:
results/kill/KILL_REPORT.md:1618: `1e-9` relative) against `zero_plus_q` and `zero_plus_corner_margin` in the doctest *"the
results/kill/KILL_REPORT.md:1621: translations plus one free translation pair per component of `Γ`; the doctest checks
results/kill/KILL_REPORT.md:1631: The header of `expansive_cone.hpp` carries this citation.
results/kill/KILL_REPORT.md:1649: LPs of the doctest *"cone_lp solves three hand-solved linear programs"* pin the solver:
results/kill/KILL_REPORT.md:1683: Extended control, `kill_common::deployable_population` (authored tilings at five clip
results/kill/KILL_REPORT.md:1740: *Deviation, stated:* `contact.hpp`'s T4.2″ scan is parameterised by the **uniform** angle
results/kill/KILL_REPORT.md:1743: itself calls at every interval midpoint — and a doctest checks that along the uniform ray
results/kill/KILL_REPORT.md:1749: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
results/kill/KILL_REPORT.md:1750: ./code/build/kiri_tests
results/kill/KILL_REPORT.md:1751: for i in $(seq 0 7); do ./code/build/kill_k8a --n 100 --x0 --dual-iters 100000 \
results/kill/KILL_REPORT.md:1760: Every graph is a deterministic function of its id (`kill_common::make_graph(id, 100, 800,
results/kill/KILL_REPORT.md:1787: `code/apps/kill_k8a_recheck.cpp`, 4 new doctests. **Nothing above is deleted; every number
results/kill/KILL_REPORT.md:1857: **Test suite.** 4 new doctests in `code/tests/test_method.cpp`: the bracket closes to `1e−9`
results/kill/KILL_REPORT.md:1866: `code/apps/kill_b3.cpp`, library `code/src/method/budget.{hpp,cpp}` (5 new doctests in
results/kill/KILL_REPORT.md:1867: `code/tests/test_method.cpp`; suite 80 cases / 16,785 assertions pass). Outputs
results/kill/KILL_REPORT.md:1873: **B1 global half (the part `derivations/scratch/check_b1.cpp` could not close).** On the
results/kill/KILL_REPORT.md:1948: **The problem solved** (`code/src/method/convex_embed.{hpp,cpp}`, derived there in full).
results/kill/KILL_REPORT.md:2021: `Θ_bisect = 0` at the historical shrink `1e−12`. `code/apps/dbg_k9.cpp` traces every one to
results/kill/KILL_REPORT.md:2051: `min μ > 0`, and `expansive_cone.hpp`'s rows equal `q/2` and `μ/2` at the `σ` flex, asserted
results/kill/KILL_REPORT.md:2052: by doctest). A strictly feasible point therefore **provably exists** on those 23 and the LP
results/kill/KILL_REPORT.md:2078: `code/apps/kill_k9.cpp`, method `code/src/method/convex_embed.{hpp,cpp}`, referee diagnostic
results/kill/KILL_REPORT.md:2079: `code/apps/dbg_k9.cpp`, doctests in `code/tests/test_method.cpp` (hand-checked
results/kill/KILL_REPORT.md:2093: **The four levers** (`code/apps/kill_k9b.cpp`), in the order the plan ranked them:
results/kill/KILL_REPORT.md:2103: * **(iv) stage 2**, `range_opt.hpp`'s softmin-of-first-contact objective (the T6 analytic
results/kill/KILL_REPORT.md:2183: logs/, k9b_gallery.png, k9b_hist.png}`. Driver `code/apps/kill_k9b.cpp`, plots
results/kill/KILL_REPORT.md:2206: `code/src/method/range_embed.{hpp,cpp}` derives it in full. The scalar is the **0⁺
results/kill/KILL_REPORT.md:2213: — the joint minimum of `zero_plus.hpp`'s two first-order separation families, the split-edge
results/kill/KILL_REPORT.md:2216: sign of `cross = det(e1, e2)` — the *same* cross that `convex_embed.hpp` constrains, since
results/kill/KILL_REPORT.md:2230: is always the exact minimum, never the surrogate. Two doctests check exactly those two
results/kill/KILL_REPORT.md:2231: claims (`code/tests/test_range_embed.cpp`), and two more finite-difference the analytic
results/kill/KILL_REPORT.md:2239: `Θ_max` with `range_opt.hpp`'s softmin-of-first-contact objective (the T6 analytic
results/kill/KILL_REPORT.md:2269: `kill_k9c --aggregate --nshards 12` and `results/kill/k9c/{k9c.csv, summary.txt}` now hold the
results/kill/KILL_REPORT.md:2423: k9c_gallery.png, k9c_hist.png}`. Driver `code/apps/kill_k9c.cpp`, method
results/kill/KILL_REPORT.md:2424: `code/src/method/range_embed.{hpp,cpp}`, doctests `code/tests/test_range_embed.cpp`, plots
results/kill/KILL_REPORT.md:2438: **How run.** `code/apps/kill_native200.cpp`, three variants per graph: *native* (their own
results/kill/KILL_REPORT.md:2515: run_shard_{}.log}`. Driver `code/apps/kill_native200.cpp`. K9's numbers for the same 200
results/kill/KILL_REPORT.md:2519: ./code/build/kill_native200 --n 200 --maxf 800 --timeout 600 --shard {} --nshards 12 \
```

### `results/kill/jitter/cert_diagnosis.md` (16)

_'C++ code in `code/`'; `contact.cpp` refs, `code/apps/dbg_cert.cpp`; regression test in `code/tests/test_method.cpp`; cmake line + test counts 66 / 16,283 and 31 / 89,055 (l.175–180); `kill_jitter`/`kill_k2a`/`kill_k5`/`kill_k6` binaries._

```
results/kill/jitter/cert_diagnosis.md:4: C++ code in `code/`; no number is quoted from memory.
results/kill/jitter/cert_diagnosis.md:38: `contact.cpp::contact_angles()` (the exact scan) applies exactly that test, with
results/kill/jitter/cert_diagnosis.md:40: `contact.cpp::validity_certificate()`'s NOROOT scan did not: it called
results/kill/jitter/cert_diagnosis.md:52: `code/apps/dbg_cert.cpp` (new) re-runs one jitter row through the identical `measure()`
results/kill/jitter/cert_diagnosis.md:95: `code/src/method/contact.cpp`, in `validity_certificate`'s NOROOT scan: build the edge's
results/kill/jitter/cert_diagnosis.md:150: `derivation_tests`' R5-a…R5-d block (which checks (T5.2b''-1b) on the corpus) passes
results/kill/jitter/cert_diagnosis.md:155: **Regression test** (`code/tests/test_method.cpp`, "F32: NOROOT counts only roots with the
results/kill/jitter/cert_diagnosis.md:161: * Without the fix (`git stash` of `contact.cpp` only, everything else identical):
results/kill/jitter/cert_diagnosis.md:175: **Suites** (`cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j`):
results/kill/jitter/cert_diagnosis.md:179: | `kiri_tests` | 66/66 test cases, **16283/16283 assertions**, SUCCESS |
results/kill/jitter/cert_diagnosis.md:180: | `derivation_tests` | 31/31 test cases, **89055/89055 assertions**, SUCCESS |
results/kill/jitter/cert_diagnosis.md:182: ## 6. Re-measurement of A3 (`kill_jitter`, 12 shards, ~2 s)
results/kill/jitter/cert_diagnosis.md:212: **K2a — the numbers change, the claims survive and get stronger.** Re-ran `kill_k2a`
results/kill/jitter/cert_diagnosis.md:241: **K5** (`kill_k5`, full 200 graphs, both sigma rules, `eps = 0.3`) — identical to the
results/kill/jitter/cert_diagnosis.md:259: **K6** (`kill_k6 --n 20 --no-native`, a smoke re-run under the same code): every clause
results/kill/jitter/cert_diagnosis.md:278: `./code/build/kill_k6 --n 20 --no-native` (40 rows = 20 graphs x 2 sigma):
```

### `results/kill/k8a/recheck.md` (8)

_`expansive_cone.cpp` line refs; `code/apps/kill_k8a{,_recheck}.cpp`; cmake + `./code/build/kiri_tests` + `kill_k8a*` command block (l.221–229)._

```
results/kill/k8a/recheck.md:22:   (`expansive_cone.cpp` lines 208–243 of the original), started from `z = 0` with step
results/kill/k8a/recheck.md:44: `code/src/method/expansive_cone.cpp`, `cone_lp`:
results/kill/k8a/recheck.md:74: `code/apps/kill_k8a_recheck.cpp` regenerates K9's variant-(b) embeddings bit-identically
results/kill/k8a/recheck.md:139: `code/apps/kill_k8a.cpp`, 100 graphs, both embeddings, 8 shards, corrected solver:
results/kill/k8a/recheck.md:221: cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j
results/kill/k8a/recheck.md:222: ./code/build/kiri_tests                      # 4 new K8a-recheck cases
results/kill/k8a/recheck.md:225: for i in $(seq 0 11); do ./code/build/kill_k8a_recheck --out results/kill/k8a \
results/kill/k8a/recheck.md:229: for i in $(seq 0 7); do ./code/build/kill_k8a --n 100 --x0 --dual-iters 100000 \
```

### `results/kill/native200/NATIVE200_FINAL.md` (1)

_one ref: `code/apps/native200_merge.cpp` produces the report._

```
results/kill/native200/NATIVE200_FINAL.md:5: Every number below is computed by `code/apps/native200_merge.cpp` from `results/kill/native200/native200_final.csv`, which is itself merged from `native200.csv` (the 600 s run), `rerun3600/shard_*.csv` (the 3600 s rerun of the timed-out ...
```

### `results/regime/REGIME.md` (1)

_one ref: `code/build/kill_regime --mode aggregate` binary._

```
results/regime/REGIME.md:6: `results/regime/summary.txt`, which `code/build/kill_regime --mode aggregate` derives from
```

### `results/yield/BASIN.md` (3)

_`code/apps/kill_basin{,_agg}.cpp` drivers; `code/src` non-modification; `code/tests/test_design.cpp`._

```
results/yield/BASIN.md:19: - Driver `code/apps/kill_basin.cpp`; aggregator `code/apps/kill_basin_agg.cpp`; figure
results/yield/BASIN.md:20:   `results/yield/plot_basin.py`. Nothing under `code/src` and nothing under `results/kill`
results/yield/BASIN.md:37:   `code/tests/test_design.cpp` ("design_range_max reproduces the K9c CSV row …") pass in
```

### `results/yield/YIELD.md` (3)

_`code/apps/kill_yield.cpp` (added to `code/CMakeLists.txt`); `code/apps/kill_k9c.cpp`; artifact table row._

```
results/yield/YIELD.md:14: `code/apps/kill_yield.cpp` (added to `code/CMakeLists.txt`) computes, per design, 45 raw
results/yield/YIELD.md:22:   same graphs and same σ as `code/apps/kill_k9c.cpp`, the Eq. (6) projection read from the
results/yield/YIELD.md:252: | `code/apps/kill_yield.cpp` | feature extractor, three modes (`k9c`, `fresh`, `check`) |
```

### `export/hero/README.md` (12)

_clang++ build line for throwaway `dump_hero_graph.cpp` (l.50–52); `tests/test_design.cpp` regression-lock refs; `apps/kill_k9.cpp`, `export/solid.hpp` refs; 'not added to the CMake build' remarks._

```
export/hero/README.md:6: regression-locked in `code/tests/test_design.cpp`), certified valid with the
export/hero/README.md:33: `tests/test_design.cpp` regression lock (`theta_max == eps_max ==
export/hero/README.md:40:    sigma + `X_ini`, bit-identical to what `apps/kill_k9.cpp` and
export/hero/README.md:41:    `tests/test_design.cpp`'s `k9_design(148)` build — this matters because
export/hero/README.md:46:    (not added to the CMake build, to avoid touching the shared `CMakeLists.txt`
export/hero/README.md:50:    clang++ -std=c++20 -O2 -arch arm64 \
export/hero/README.md:51:      -I code/src -I code/apps -I /opt/homebrew/include -I /opt/homebrew/include/eigen3 \
export/hero/README.md:52:      dump_hero_graph.cpp code/build/libkiri_core.a -o dump_hero_graph
export/hero/README.md:55:    `dump_hero_graph.cpp` just calls `kiri::kill::make_graph(148, 100, 800, 1400)`
export/hero/README.md:96:    `--theta-frac` computes `theta_max` internally (`export/solid.hpp`'s
export/hero/README.md:156: `dump_hero_graph.cpp` (the throwaway graph-dump utility from step 1) was
export/hero/README.md:158: the CMake build — it is a 20-line wrapper around
```

### `export/hero2/README.md` (21)

_three clang++ build lines for throwaway `dump_hero2_input.cpp`/`reconstruct_k9c_graph.cpp`/`characterize_k9c.cpp` (l.109–137); `method/design.cpp`, `range_embed.hpp`, `kill_k9c.cpp`, `export/solid.hpp` refs; CMake remarks._

```
export/hero2/README.md:23: `kiri_design --maximise-eps` runs `design_constrained` (`method/design.cpp`), whose stage-1
export/hero2/README.md:27: **stage-A margin-maximisation** (`method/range_embed.hpp`), warm-started from **t = 0**
export/hero2/README.md:32: (`dump_k9c_graph.cpp`, not committed) that replicates `kill_k9c.cpp`'s exact per-graph
export/hero2/README.md:47: (`reconstruct_k9c_graph.cpp`) rebuilds the *same deterministic, non-random* topology
export/hero2/README.md:52: A third throwaway program (`characterize_k9c.cpp`) then calls `method::characterize` — the
export/hero2/README.md:56: added to the CMake build; they are graph-topology/measurement plumbing, not new method
export/hero2/README.md:106:    `apps/kill_k9c.cpp` builds for `id=130`), same throwaway-utility pattern as hero's
export/hero2/README.md:107:    `dump_hero_graph.cpp`:
export/hero2/README.md:109:    clang++ -std=c++20 -O2 -arch arm64 \
export/hero2/README.md:110:      -I code/src -I code/apps -I /opt/homebrew/include -I /opt/homebrew/include/eigen3 \
export/hero2/README.md:111:      dump_hero2_input.cpp code/build/libkiri_core.a -o dump_hero2_input
export/hero2/README.md:114:    (`dump_hero2_input.cpp` calls `kiri::kill::make_graph(130, 100, 800, 1400)`, K9c's
export/hero2/README.md:120:    clang++ -std=c++20 -O2 -arch arm64 \
export/hero2/README.md:121:      -I code/src -I code/apps -I /opt/homebrew/include -I /opt/homebrew/include/eigen3 \
export/hero2/README.md:122:      reconstruct_k9c_graph.cpp code/build/libkiri_core.a -o reconstruct_k9c_graph
export/hero2/README.md:135:    clang++ -std=c++20 -O2 -arch arm64 \
export/hero2/README.md:136:      -I code/src -I code/apps -I /opt/homebrew/include -I /opt/homebrew/include/eigen3 \
export/hero2/README.md:137:      characterize_k9c.cpp code/build/libkiri_core.a -o characterize_k9c
export/hero2/README.md:165:    `--theta-frac` reports `theta_max=3.14159` at every fraction (export/solid.hpp's
export/hero2/README.md:223: `dump_hero2_input.cpp`, `reconstruct_k9c_graph.cpp`, `characterize_k9c.cpp` (the
export/hero2/README.md:225: the repository or the CMake build — they are graph-topology dump / fold-back /
```

### `review/constructive_review.md` (1)

_one file path ref `code/src/method/convex_embed.hpp` (l.17)._

```
review/constructive_review.md:17: **Mechanically: yes, it is a standard L-BFGS/softplus/log-barrier feasibility-then-proximity solve over two quadratic-in-t constraint families** (`code/src/method/convex_embed.hpp`): corner convexity (`cross_i(t) ≥ δ`) and split-inward s...
```

### `review/negative_review.md` (1)

_one file:line ref `code/src/core/generators.cpp:496-580`._

```
review/negative_review.md:83: triangulations of uniform random points (`code/src/core/generators.cpp:496-580`) are a
```

### `review/theory_review.md` (1)

_one ref `kinematics.cpp::deploy()`._

```
review/theory_review.md:24:    convention (`derivations/core.md:190-390`, Check C2 vs `kinematics.cpp::deploy()`, 1.42e−14). Acuña et al.
```

### `notes/paper_2026.md` (1)

_KEEP — quotes the 2026 paper's own 'WebAssembly backend compiled from C++' sentence (third-party)._

```
notes/paper_2026.md:853: > "The core algorithms are implemented in a WebAssembly backend, compiled from C++ code. We optimized all the non-linear energies (discussed in Sec. 4.5, 5.1 and 5.2) via Newton's method with a per-element Hessian projection onto the PSD...
```

### `notes/repo_2025.md` (116)

_KEEP — analysis of the AUTHORS' 2025 C++/emscripten repo (`kiri_mesh.cpp`, `multi_grid.cpp`, `bind.cpp`, their CMakeLists, Eigen); third-party, not our code._

```
notes/repo_2025.md:12: Total hand-written source: ~2 600 lines of C++ (excluding the vendored
notes/repo_2025.md:13: `code/cpp/json.hpp`) plus ~2 500 lines of JS. It is small; this note covers all of it.
notes/repo_2025.md:30: - **Does the C++ build natively? Not as shipped — it is emscripten-only, but only
notes/repo_2025.md:32:   1. `code/CMakeLists.txt:43` sets `target_link_options(... -sMODULARIZE
notes/repo_2025.md:33:      -sEXPORT_NAME="createModule" --bind -sWASM=1 -sALLOW_MEMORY_GROWTH=1
notes/repo_2025.md:35:      `code/CMakeLists.txt:45-53` post-build steps copy `kiri_patterns.js/.wasm`.
notes/repo_2025.md:36:   2. `code/cpp/bind.cpp` includes `<emscripten/bind.h>` (`bind.cpp:10`) and is the only
notes/repo_2025.md:38:      `file(GLOB SRC_FILES cpp/*.cpp ...)` (`CMakeLists.txt:39`), so a native build must
notes/repo_2025.md:43:      *not pinned at all* — a reproducibility hazard), Eigen 3.4.0 as a fallback.
notes/repo_2025.md:44: - Actual dependency surface is small: Eigen everywhere; **libigl only for
notes/repo_2025.md:45:   `igl/AABB.h`** (`opt/opt.cpp:9`, `parameterization/lift.cpp:7`); Optiz for
notes/repo_2025.md:47:   CGAL, no solver beyond Eigen's sparse. So a native CLI is straightforward, and the
notes/repo_2025.md:49: - A **prebuilt `code/kiri_patterns.wasm` + `code/kiri_patterns.js` are committed**, so
notes/repo_2025.md:69: `multi_grid.cpp:167-194` (`MultiGrid::from_json`) and `bind.cpp:56-67`
notes/repo_2025.md:71: **`p /= 100` on every coordinate** (`multi_grid.cpp:177-179`, `bind.cpp:63-65`).
notes/repo_2025.md:73: **Units.** JSON coordinates are **HTML-canvas pixels**; the C++ works in units of
notes/repo_2025.md:76: scale-free: the target mesh is normalized to unit area (`bind.cpp:155-156`) and the
notes/repo_2025.md:84: | **bare array** | array of the same 3-point group objects | `fig13_col1..fig13_col6` | `MultiGrid::from_json` (`multi_grid.cpp:167`) — the **offline C++ path**; these are **not in the UI dropdown** (`index.html:479-489`) |
notes/repo_2025.md:87: (`js/design-module.js:40-50, 144-175`) and are never passed to C++. `color` is
notes/repo_2025.md:98: (verbatim port, `multi_grid.cpp:25-46` and `utils.cpp:31-65`) on all 16 files.
notes/repo_2025.md:123:   to the deployed appearance, not to any origami waterbomb crease pattern. The C++
notes/repo_2025.md:202: ## 3. The C++ pipeline
notes/repo_2025.md:208:   → MultiGrid              (multi_grid.cpp)      arrangement + unit-cell extraction
notes/repo_2025.md:209:   → SymPattern             (sym_pattern.cpp)     half-edge pattern, faces, cut labels
notes/repo_2025.md:210:   → SymPattern::mesh(r,r)  (sym_pattern.cpp:28)  tile the unit cell r×r → Hmesh
notes/repo_2025.md:219: ### 3.1 The complete cut (paper Supp. A) — `kiri_mesh.cpp:17-54`
notes/repo_2025.md:228: if (e->ti != -1) q.push({e->twin(), !is_cut});      // kiri_mesh.cpp:46-50
notes/repo_2025.md:234: relaxation, nothing. `SymPattern::determine_cut_edges()` (`sym_pattern.cpp:109-127`) is
notes/repo_2025.md:241: > (`kiri_mesh.cpp:29-31`)
notes/repo_2025.md:259:                                                           // kiri_mesh.cpp:38-43
notes/repo_2025.md:270: ### 3.2 Deployment-friendliness (paper Def. 4.1, Prop. 4.1, Supp. B) — `kiri_mesh.cpp:169-239`
notes/repo_2025.md:274: of the uncut mesh (boundary vertices skipped, `kiri_mesh.cpp:192-193`) it walks the star
notes/repo_2025.md:279:   angles** at that step (`kiri_mesh.cpp:208-211`),
notes/repo_2025.md:286: return point.norm() < 1e-5;                        // kiri_mesh.cpp:180-189
notes/repo_2025.md:291: `res = 1` on the corners of their *non-cut* incident faces (`kiri_mesh.cpp:227-236`),
notes/repo_2025.md:302:   (`kiri_mesh.cpp:213`) but then `if (iter > 10) { ...; exit(1); }`
notes/repo_2025.md:303:   (`kiri_mesh.cpp:215-219`). Any vertex whose cut-star cycle has more than ~10 steps
notes/repo_2025.md:307: ### 3.3 Unit tile and periodicity (paper Sec. 4.2, Supp. D) — `multi_grid.cpp:25-156`
notes/repo_2025.md:309: **`MultiGrid::calculate_periodicity()` (`multi_grid.cpp:25-46`) is the algorithm the
notes/repo_2025.md:317:    `approximate_ratio(|t0/t|) → (n, m)` by **continued fractions** (`utils.cpp:31-65`,
notes/repo_2025.md:321:    (`multi_grid.cpp:33-36`). A returned `(100,100)` is the sentinel for
notes/repo_2025.md:323:    **"No unit pattern found."** (`bind.cpp:69-72`).
notes/repo_2025.md:330: `get_sym_pattern()` (`multi_grid.cpp:52-156`) then builds the actual arrangement inside
notes/repo_2025.md:335: Vertex identification is an **O(n²) linear scan at 1e-4** (`multi_grid.cpp:134-142`).
notes/repo_2025.md:342: Face extraction: `SymPattern::build_data_structure()` (`sym_pattern.cpp:129-157`) adds
notes/repo_2025.md:344: (`find_next_edge`, `sym_pattern.cpp:64-81`), then `build_faces()` (`sym_pattern.cpp:170-199`)
notes/repo_2025.md:347: translation by the `periodicity` rows (`sym_pattern.cpp:28-57`), and the seams are then
notes/repo_2025.md:348: welded by `merge_close_verts()` at the call sites (`bind.cpp:75,95,103`).
notes/repo_2025.md:350: `SymPattern::get_bounding_parallelogram(opening_angle)` (`sym_pattern.cpp:83-107`)
notes/repo_2025.md:354: only for visual alignment (`bind.cpp:120-125`) and for the lattice basis in `lift`.
notes/repo_2025.md:356: ### 3.4 Deployment and θ_max — `kiri_mesh.cpp:60-86, 88-151`
notes/repo_2025.md:362: struct RotatingInfo { Eigen::Vector2d center; bool clockwise; int prev_face; };
notes/repo_2025.md:369: (`kiri_mesh.cpp:131-135`). `open_t` then composes the chain from the face back to the
notes/repo_2025.md:386: stepping.** `KiriMesh::max_opening_angle()` (`kiri_mesh.cpp:60-86`):
notes/repo_2025.md:410:   `rep_x=rep_y=2`, `bind.cpp:135-138, 168-171`), so it is a per-unit-tile quantity.
notes/repo_2025.md:424: ### 3.5 Deployed-state vertex fusing — `Hmesh.cpp:270-294`, `lift.cpp:152-161`
notes/repo_2025.md:432:                                                           // Hmesh.cpp:273-284
notes/repo_2025.md:436: (`lift.cpp:152-161`) to three different vertex sets:
notes/repo_2025.md:459: ### 3.6 Inverse design — `opt/opt.cpp`, `parameterization/{param,lift}.cpp`
notes/repo_2025.md:461: **Initialization by parameterization** (`bind.cpp:151-181`):
notes/repo_2025.md:464:    `bind.cpp:155-156`), then `param::isometric_param`:
notes/repo_2025.md:465:    - `init_vf` (`param.cpp:24-47`): a smooth **unit vector field** on faces by a
notes/repo_2025.md:471:    - `init_uv` (`param.cpp:49-77`): a Poisson solve aligning `∇u, ∇v` to that field
notes/repo_2025.md:474:      (`param.cpp:92-99`), i.e. an as-isometric-as-possible flattening.
notes/repo_2025.md:475:    - The result is centred (`bind.cpp:159`).
notes/repo_2025.md:478:    (`state.h:8`) and defined (`state.cpp:5`) and **never read or written anywhere**.
notes/repo_2025.md:483: 2. `liftPattern` (`bind.cpp:163-181`) fixes the two free parameters automatically:
notes/repo_2025.md:488:      `// Scale to achieve an area_ratio of 5.` (`bind.cpp:173-175`).
notes/repo_2025.md:490:      (`bind.cpp:183-198`, `js/algorithm-module.js:693-703`) with `scale=1, rot=0`
notes/repo_2025.md:493: 3. `param::lift` (`lift.cpp:36-164`) tiles the deployed unit pattern over the UV domain:
notes/repo_2025.md:495:    deployed lattice `basis` is read off the 2×2 block (`lift.cpp:63-69`); the integer
notes/repo_2025.md:497:    (`lift.cpp:71-78`); each replica vertex is located in the UV triangulation by an
notes/repo_2025.md:499:    (`lift.cpp:90-104`); a face is kept **only if all its vertices land inside the UV**
notes/repo_2025.md:500:    (`lift.cpp:106-118`). Optional `remove_dangling_faces` iteratively drops faces with
notes/repo_2025.md:501:    fewer than two non-boundary edges (`lift.cpp:127-150`) — but `bind.cpp:178,188`
notes/repo_2025.md:504: **The energy — this is Supplement Sec. F's numbers.** Declared at `opt/opt.cpp:15-17`:
notes/repo_2025.md:516: | `ω1 E_reconfig` | `rigid_weight * get_rigid_error` (polygon edges) and `* get_face_rigid_error` (all other in-face pairs) | `opt.cpp:134, 178` | **ω1 = 10.0** |
notes/repo_2025.md:517: | `ω2 E_planar` | `planarity_weight * err` | `opt.cpp:175` | **ω2 = 10** |
notes/repo_2025.md:518: | `ω3 E_shape` | `closeness_weight * sqr((x_p − cp)·n)` | `opt.cpp:147` | **ω3 = 0.1** |
notes/repo_2025.md:519: | `ω4 E_fairness` | `close_to_init_weight * ‖x_p − Y⁰‖²` | `opt.cpp:154` | **ω4 = 0.5 static, but `reset_optimization()` sets it to `orig_close_to_init_weight = 10`, then anneals `×0.8` every 20 iterations** (`opt.cpp:186, 203-204`) |
notes/repo_2025.md:520: | — | `smoothness_weight = 0.1` | `opt.cpp:17` | **declared and never used anywhere.** Dead. |
notes/repo_2025.md:523: (`bind.cpp:179, 189`), so the **effective ω4 at iteration 0 is 10**, not 0.5, and it
notes/repo_2025.md:532: return sqr((v1-v2).norm() - (v3-v4).norm()) / (l1*l1);       // opt.cpp:48-52
notes/repo_2025.md:539: (`Optiz::Problem prob(state::opt_lifted, state::opt_ground)`, `opt.cpp:127`), i.e. `X`
notes/repo_2025.md:541: (`opt.cpp:57-85`) covers the non-edge in-face pairs `j ≥ i+2` (skipping the wrap-around
notes/repo_2025.md:544: `max_face_degree() > 3`** (`opt.cpp:158`) — for a pure triangle tiling the edge term
notes/repo_2025.md:547: **Planarity** (`opt.cpp:159-176`): the face normal `n` is computed by PCA (smallest
notes/repo_2025.md:548: eigenvector of the centred covariance, `opt.cpp:96-102`) **from the current iterate and
notes/repo_2025.md:553: **Shape/closeness** (`opt.cpp:138-148`): `igl::AABB` closest-point query against the
notes/repo_2025.md:558: own library. Configuration (`opt.cpp:128-130`):
notes/repo_2025.md:592: | reported metrics | `rigid_avg/max`, `close_avg/max`, `planarity_avg/max` | `bind.cpp:223-233`, `js/algorithm-module.js:667-670` |
notes/repo_2025.md:596: emscripten bindings. Note `sw.js:43` calls
notes/repo_2025.md:598: (`merged_mesh`) is **ignored by the C++** (`bind.cpp:91` takes it and never reads it).
notes/repo_2025.md:606: | **A** — complete cut construction | the algorithm producing `E_c` | **Yes, fully.** `kiri_mesh.cpp:17-54` (and `sym_pattern.cpp:109-127`): alternating BFS 2-colouring, one boolean flipped across twins, corner duplication paired across ...
notes/repo_2025.md:607: | **B** — proof of Prop. 4.1, `P₂ ~ P₁` | the proof by construction | **No.** Only the *predicate*: the polygon-closure test at tolerance `1e-5` (`kiri_mesh.cpp:180-189`). A proof cannot be recovered from it. |
notes/repo_2025.md:609: | **D** — periodicity / unit tile (`λ₁, λ₂`) | the computation | **Yes, fully.** `multi_grid.cpp:25-46`: continued-fraction commensurability of family spacings, with the hard `n,m < 10` rationality cut-off and the `(100,100)` "no unit pa...
notes/repo_2025.md:614: - **the ω4 annealing schedule** (`10 · 0.8^{floor(k/20)}`, `opt.cpp:203-204`);
notes/repo_2025.md:616:   area ratio of 5 (`bind.cpp:172-175`).
notes/repo_2025.md:630:   (`kiri_mesh.cpp:181-227`) and `sum_angles` in face orientation classification
notes/repo_2025.md:631:   (`sym_pattern.cpp:179-193`).
notes/repo_2025.md:633:   (`opt.cpp:127`) — **positions only**. No angle is ever a degree of freedom.
notes/repo_2025.md:636:   (`kiri_mesh.cpp:145-167`, `lift.cpp:52-53`, `sym_pattern.cpp:90,99`).
notes/repo_2025.md:645: exists): `SymPattern::make_friendly()` (`sym_pattern.cpp:59-62`, computes a mesh,
notes/repo_2025.md:646: discards it, returns `*this`), `SymPattern::from_json_graph` (`sym_pattern.cpp:201-225`,
notes/repo_2025.md:648: `param::flip_combinatorics` (`lift.cpp:231-251`, not even declared in `lift.h`),
notes/repo_2025.md:659: `opt::get_errors` (`opt.cpp:207-257`), surfaced to JS by `bind.cpp:223-233`, is the
notes/repo_2025.md:664: | `rigid_avg`, `rigid_max` | `\|‖v_i−v_j‖ − ‖g_i−g_j‖\| / l₁` with `l₁ = ½(‖v_i−v_j‖+‖g_i−g_j‖)`, over all polygon edges plus (if `max_face_degree > 3`) all in-face pairs (`opt.cpp:214-224`, `relative=true, sqrd_err=false`) | dimensionle...
notes/repo_2025.md:665: | `close_avg`, `close_max` | `\|(p − cp)·n\|`, point-to-plane distance to the target (`opt.cpp:87-94`) | model units (target normalized to unit area) | the shape-approximation error; **the paper reports no threshold for it** |
notes/repo_2025.md:666: | `planarity_avg`, `planarity_max` | `Σ_j \|π/2 − acos(n̂·v̂_j)\| / π · 180` over consecutive edges of each face with >3 vertices, `n̂` the PCA normal (`opt.cpp:104-116`) | **degrees** | **average planarity error `E_p`**; fabricability c...
notes/repo_2025.md:673: 2. The planarity sum is **not normalized by face degree** (`opt.cpp:110-114` sums over
notes/repo_2025.md:693:    `KiriMesh::open_separate_faces(θ)` (`kiri_mesh.cpp:153-167`) returns each face as its
notes/repo_2025.md:695:    `get_bounding_parallelogram(θ)` (`sym_pattern.cpp:83-107`) minus the sum of face
notes/repo_2025.md:696:    areas (`Hmesh::Face::area()`, `Hmesh.cpp:64-82`, exact for planar polygons) gives the
notes/repo_2025.md:721:    `planarity_*`; a native build (§0) plus a small CLI driver replacing `bind.cpp` would
notes/repo_2025.md:735: injectivity test, and of hole-area logging; the license statement; the emscripten-only
notes/repo_2025.md:750: - whether the C++ compiles natively — **assessed, not attempted**, per instructions;
```

### `notes/screen_r2.md` (3)

_`periodic_jacobian.hpp` ref; `kill_k8f`/`kill_i2` proposed driver names._

```
notes/screen_r2.md:253:    deployment angle, valid for **any** periodic corner-hinged pattern via `periodic_jacobian.hpp`, is
notes/screen_r2.md:263: `kill_k8f` / `kill_i2` must be run knowing that **the published maximum area expansion of the
notes/screen_r2.md:516: 3. **Do not run `kill_i2`'s stated kill rule as written.** "`B(triangles) >= 4` kills the idea" fires
```

### `ideas/persona_adversary.md` (9)

_HISTORICAL (ideation) — scratch `check.cpp`/`rig_check.cpp` refs, `code/src/core/` assumption, 'C++-only stack'/Eigen effort estimates, authors' WebAssembly mention._

```
ideas/persona_adversary.md:53: **True.** Both personas verified it numerically (geometer `check.cpp`, 2.0e-15 on a 5×5 patch;
ideas/persona_adversary.md:54: rigidity `rig_check.cpp` V1, 6e-15 on 21/21 graphs including Voronoi patches with split cuts). I
ideas/persona_adversary.md:405: (rigidity (R2)). All ideas assume the reference pipeline in `code/src/core/` exists as specified in
ideas/persona_adversary.md:785: **Kill experiment.** Rerun exactly the rigidity persona's `rig_check.cpp` configuration — the same
ideas/persona_adversary.md:1080: **Effort.** Derivation M. Code L (an SDP solver is not in the C++ stack; Shor's relaxation for
ideas/persona_adversary.md:1081: these sizes can be done with an Eigen-based projected-gradient / spectral bundle, which is real
ideas/persona_adversary.md:1082: work — this is the one idea in my file with a genuine implementation risk under the C++-only
ideas/persona_adversary.md:1100: > diameter-selection rule; the authors ship a WebAssembly implementation at a public repository.
ideas/persona_adversary.md:1189: | 9 | **A10 — emptiness of the usable design region** | Attacks the word "full" in the abstract, and the SDP infeasibility certificate is the only theorem-grade emptiness result available. Ranked ninth because the most likely outcome is ...
```

### `ideas/persona_geometer.md` (13)

_HISTORICAL — scratch `check*.cpp` programs, hypothetical `core/*.cpp` routine names, `Eigen::ColPivHouseholderQR`, 'code/src/core does not exist yet'._

```
ideas/persona_geometer.md:104: **Verified numerically** (`scratchpad/check3.cpp`, 400 random-`sigma` trials on a 6x6 patch): the
ideas/persona_geometer.md:162: Three throwaway C++ programs in the session scratchpad
ideas/persona_geometer.md:163: (`.../scratchpad/check.cpp`, `check2.cpp`, `check3.cpp`; C++17, `clang++ -arch arm64`, no
ideas/persona_geometer.md:164: dependencies) on square-grid patches and tori, since `code/src/core/` does not exist yet.
ideas/persona_geometer.md:219: `core/forward_kinematics.cpp` at 50 values of `theta in (0, theta_max)` on **200 graphs**
ideas/persona_geometer.md:266: periodic Voronoi metatiles), solve Eq. (13) with `core/conformal_solve.cpp`, then sample `theta` on a
ideas/persona_geometer.md:317: `core/theta_max_bisect.cpp` at tolerance `1e-8`. **If `|theta_max^(a) - theta_max^(b)| > 1e-6` on any
ideas/persona_geometer.md:362: (`scratchpad/check3.cpp`): 0/400 failures of the corrected identity, plus 0 failures of Remark A.1.
ideas/persona_geometer.md:365: split-heavy regime), run `core/hole_preimage.cpp` (Algorithm 1 with the three prose rules) to get
ideas/persona_geometer.md:413: compute `rank(L)` by rank-revealing QR (`Eigen::ColPivHouseholderQR`, tol `1e-10 * ||L||`) and
ideas/persona_geometer.md:576: polygon, solve `Delta_h X = 0` with `core/tutte_auxetic_solve.cpp`, and compute the signed area of
ideas/persona_geometer.md:628: pin-joint constraint Jacobian (`3|F|` DOF, `2|E_hinge|` constraints) from `core/rigidity_jacobian.cpp`.
ideas/persona_geometer.md:777:   `Gamma`; I could not verify that against Algorithm 1 itself, because `code/src/core/` does not yet
```

### `ideas/persona_optimizer.md` (5)

_HISTORICAL — `code/src/core/*.hpp`, `kinematics.hpp` convention refs; 'substantial piece of C++'; 'C++ CSV dump'._

```
ideas/persona_optimizer.md:14: - `code/src/core/*.hpp` — to state kill experiments against routines that exist.
ideas/persona_optimizer.md:21: Re-derivation from the convention in `code/src/core/kinematics.hpp`: face `f` maps
ideas/persona_optimizer.md:429: and nearly bipartite. Second risk: a blossom implementation is a substantial piece of C++ with no
ideas/persona_optimizer.md:665: < 30 min; correlation plots in matplotlib from the C++ CSV dump.
ideas/persona_optimizer.md:809:   prediction. Fact (i) I re-derived by hand from `code/src/core/kinematics.hpp`'s stated convention;
```

### `ideas/persona_rigidity.md` (9)

_HISTORICAL — standalone C++ compile against `code/src/core/*.cpp` with clang++; `CMakeLists.txt`/`collision.cpp` not-yet-existing remarks; Eigen GSVD remark._

```
ideas/persona_rigidity.md:17: I compiled a standalone C++ program against the existing core
ideas/persona_rigidity.md:18: (`code/src/core/{mesh,cut,holes,kinematics,generators,orientation,tutte_auxetic}.cpp`;
ideas/persona_rigidity.md:19: the repo's own `CMakeLists.txt` does not configure yet because `collision.cpp`/`kiri_sweep`
ideas/persona_rigidity.md:20: are still missing, so I compiled the sources directly with clang++ `-arch arm64`).
ideas/persona_rigidity.md:21: Source: `scratchpad/rig_check.cpp`. 21 graphs: squares 3×3 and 5×5, triangles, hexagons,
ideas/persona_rigidity.md:130: **Kill experiment.** `code/src/core` + a new `mobility.cpp`. Generators: `delaunay_of_random_points`,
ideas/persona_rigidity.md:240: `code/src/core/collision.cpp` lands) on 500 graphs; **if the closed-form and bisection `θ_max`
ideas/persona_rigidity.md:347: **Effort.** Derivation M. Code M (needs a numerically careful pencil staircase; Eigen has GSVD only
ideas/persona_rigidity.md:677:   because `code/src/core/collision.cpp` does not exist yet (the repo's `CMakeLists.txt` references it
```

### `ideas/ranking.md` (11)

_'All C++17/20 per directive D3' (l.462) and 'SDP-in-C++' (l.677); `collision.hpp` refs; `code/src/core` refs; authors' WebAssembly mention._

```
ideas/ranking.md:112: | B3 | Exact contact angles: every orientation predicate is `p + q cos th + r sin th` | char | G3, rig-3b, opt-3, adv-(b)/A1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | closed-form `theta_max` vs `collision.hpp` bisection on 300 graphs to 1e-6 | S+M |...
ideas/ranking.md:130: | B21 | The usable region `U(eps) = {X embedding} ∩ {theta_max >= eps}` can be small or empty | char+thm | A10 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | sampling failure rate over 1e4 points of `X` on 500 graphs; SDP infeasibility certificate | M+L ...
ideas/ranking.md:133: | B24 | Certified global inverse design by low-rank relaxation | alg | opt-9 | 1 | 1 | 1 | 0 | 1 | 1 | 1 | is the dual bound within 3x of the best achieved value on any instance? | M+L | **DROP as standalone; keep the certificate.** Same...
ideas/ranking.md:200:   in `code/src/core/collision.hpp` with the documented gamma ladder.
ideas/ranking.md:265: > the reconstruction. That is not a baseline; it is a self-portrait. The authors ship a WebAssembly
ideas/ranking.md:461: Written so an Experimenter can run them against `code/src/core` with no further design decisions.
ideas/ranking.md:462: All C++17/20 per directive D3; Python only for matplotlib on the dumped CSV. All routines named below
ideas/ranking.md:502: - Cross-check every reported `theta*` against `collision.hpp`'s bisection to `<= 1e-6` rad.
ideas/ranking.md:514:   `collision.hpp` (grid scan + bisection, faces shrunk by relative `1e-6` per deviation 9).
ideas/ranking.md:528: - **Both** measured by the same `collision.hpp` bisection, never by their own objectives.
ideas/ranking.md:677: SDP-in-C++ requirement under directive D3 is the largest implementation risk in any of the files).
```

### `ideas/ranking_r2.md` (21)

_`kill_k8a`/`kill_t1`/`kill_b3`/`kill_b2`/`kill_k8b`/`kill_k8h` driver names (binary names -> apps/*.jl); `zero_plus.hpp`, `mobility.hpp`, `contact.hpp` refs; `kill_k7.cpp`._

```
ideas/ranking_r2.md:152: **Verdict: PURSUE — run `kill_k8a` first.** It is under ten minutes, its outcome is genuinely
ideas/ranking_r2.md:248: **Verdict: PURSUE-IF-KILL-PASSES**, and the kill is five minutes with no new code (§4, `kill_t1`).
ideas/ranking_r2.md:249: This is the cheapest idea in the entire round-2 bundle. Run it alongside `kill_k8a`.
ideas/ranking_r2.md:289: ### `kill_t1` — T-1's vertex balance defect (5 min, no new code, run first)
ideas/ranking_r2.md:291: - **Routine.** `zero_plus.hpp::corner_incidences` at the K6 `X0`; `make_cut` for the split-edge set.
ideas/ranking_r2.md:302: ### `kill_k8a` — X1's expansive cone LP (< 10 min, run second)
ideas/ranking_r2.md:306: - **Routine.** `mobility.hpp::build_A`; dense `ColPivHouseholderQR` kernel basis on the **2-core
ideas/ranking_r2.md:308:   `zero_plus.hpp::split_copies` and `::corner_incidences` for the inequality rows; branch at each
ideas/ranking_r2.md:321: ### `kill_b3` — X2's `tr K` threshold and constrained re-solve (after Gate 0)
ideas/ranking_r2.md:337: ### `kill_b2` — the `O(|E|)` budget bound (5 min, shares B1's app)
ideas/ranking_r2.md:346: ### `kill_k8b` — X3's emptiness certificate (< 20 min, only if B4 FAILS)
ideas/ranking_r2.md:357: ### `kill_k8h` — A4's certified non-uniform continuation (30 min, **gated on `kill_k8a` PASS**)
ideas/ranking_r2.md:360:   step only when `contact.hpp::swept_discs` with the exact radius `max(‖x‖,‖χ‖)` (KILL_REPORT
ideas/ranking_r2.md:366: `kill_t1` and `kill_k8b` read the same K6/K1a populations and can share a driver.
ideas/ranking_r2.md:367: `kill_b2` and B1's identity check are one app. `kill_b3` is inside `kill_k7.cpp`.
ideas/ranking_r2.md:482: **The algorithm.** X2 (`tr K` budget-constrained periodic design) if `kill_b3` passes; X1's LP flex
ideas/ranking_r2.md:483: if `kill_k8a` passes and `kill_b3` does not; both if both pass, in which case X1 is the patch-space
ideas/ranking_r2.md:522: 2. `kill_t1` — five minutes, no new code, settles the §1 adjudication and the fate of T-1.
ideas/ranking_r2.md:523: 3. `kill_k8a` — ten minutes, settles whether the project has a positive half on patches.
ideas/ranking_r2.md:525: 5. `kill_b3` — settles whether K7 C3 becomes the paper's algorithm or stays a dead end.
ideas/ranking_r2.md:527: **Commit to X2 as the constructive core** unless `kill_b3` fails hard, because it is the only
```

### `ideas/round2_adversary.md` (40)

_HISTORICAL — `code/src/method/*.hpp` inventory, proposed `kill_k8a..k8j.cpp` drivers, `kill_common.hpp`; kill-rule table with `kill_k8*` names._

```
ideas/round2_adversary.md:14: `code/src/method/*.hpp` in full (`zero_plus.hpp`, `contact.hpp`, `deploy_basis.hpp`,
ideas/round2_adversary.md:15: `mobility.hpp`, `periodic_jacobian.hpp`), `code/apps/` listing,
ideas/round2_adversary.md:151:    the global sign being fixed once and for all by the stored-CCW convention of `zero_plus.hpp`.
ideas/round2_adversary.md:153:    1.9e-14 over 4300 split-edge samples) and with `zero_plus.hpp`'s "`q_e` is a positive multiple of
ideas/round2_adversary.md:157:    `mu` of `zero_plus.hpp` are also bilinear in `(u, X)`.
ideas/round2_adversary.md:202: translation `w_f` per face; `A(X,0)` is `mobility.hpp::build_A` and `ker A` has dimension `m`
ideas/round2_adversary.md:218: seen it: `zero_plus.hpp` searches over the **embedding** `X` with the flex slaved to it
ideas/round2_adversary.md:223: **Kill experiment.** Driver `kill_k8a.cpp`, reusing `mobility.hpp::build_A`,
ideas/round2_adversary.md:224: `zero_plus.hpp::split_copies` / `corner_incidences`, and the K1a population generator in
ideas/round2_adversary.md:225: `kill_common.hpp`. For each of the **100** K1a graphs (101-793 faces) at `X = X_ini` (which is
ideas/round2_adversary.md:262: coordinates `t` (`zero_plus.hpp` states this and `check.md` T3.2/T5-a verify the exact quadratic
ideas/round2_adversary.md:301: **Kill experiment.** Driver `kill_k8b.cpp`. Population: the **50** smallest K1a graphs (so
ideas/round2_adversary.md:368: **Kill experiment.** Driver `kill_k8c.cpp`. Population: the same 200 K1a graphs. Compute
ideas/round2_adversary.md:369: `sigma_mc` (existing, `orientation.hpp`), then run a face-flip local search on
ideas/round2_adversary.md:424: **Kill experiment.** Driver `kill_k8d.cpp`. Two populations. **(i) Jitter ladder:** each of the 8
ideas/round2_adversary.md:490: **Kill experiment.** Driver `kill_k8e.cpp`. For each of the 200 K1a graphs at `X0` from Eq. (6):
ideas/round2_adversary.md:493: `holes.hpp`); (b) for graphs with `Theta_max = 0`, check whether the hole minimising `A_C'(0)` is the
ideas/round2_adversary.md:540: `periodic_jacobian.hpp`), so `det g(theta) = det J(theta)^2` is an explicit trigonometric polynomial
ideas/round2_adversary.md:543: macroscopic Jacobian (`periodic_jacobian.hpp`) to the microscopic hole harmonic (idea 5). A surface
ideas/round2_adversary.md:549: **Kill experiment.** Driver `kill_k8f.cpp`, using `baseline/kirigami_tessellations` read-only.
ideas/round2_adversary.md:550: For each of the **16** pattern JSONs: build the quotient (`periodic_jacobian.hpp::build_quotient`),
ideas/round2_adversary.md:610: **Kill experiment.** Driver `kill_k8g.cpp`. On `snub_square`, `trunc_square` (4.8.8) and
ideas/round2_adversary.md:660: swept-volume bound (`contact.hpp::swept_discs` with the exact radius `max(||x||,||chi||)` from
ideas/round2_adversary.md:663: **Kill experiment.** Driver `kill_k8h.cpp` on **20** K1a graphs with `Theta_max = 0`. Take the LP
ideas/round2_adversary.md:709: **Kill experiment.** Driver `kill_k8i.cpp` plus a read-only harness over
ideas/round2_adversary.md:759: **Kill experiment.** Driver `kill_k8j.cpp`. On the 8 authored tilings and 100 K1a graphs, count
ideas/round2_adversary.md:821: | 1 expansive cone LP | F26/F29 mobility median 305 after 2-core; F30 emptiness is a `0+` split/corner event; T2.4 no locking | `kill_k8a`: positive certified LP margin on < 20/100 K1a graphs, or `sigma` reported outside `P(X)` on the 8 ...
ideas/round2_adversary.md:822: | 2 emptiness certificate | F25/F30 "0 of 400 samples", which is not a proof; `zero_plus.hpp` exact quadratic forms | `kill_k8b`: certificate on < 25/50 graphs, or **any** certificate on the 8 tilings (soundness bug) |
ideas/round2_adversary.md:823: | 3 sign-aware orientation | F30: `sigma_def` doubled split cuts 126 -> 321 and gained nothing; K6 row 1: 162/391 split edges inward | `kill_k8c`: certified valid 0/200 again, or `N_bad` not halved by local search |
ideas/round2_adversary.md:824: | 4 deployability transition | F25 0/200 random vs 0/8 tiling failures - the generator confound | `kill_k8d`: all candidate scalars AUC < 0.80, or no monotone collapse on the jitter ladder |
ideas/round2_adversary.md:825: | 5 hole-area budget | F11 preimages partition `E_hinge ∪ E_split`; optimizer's `A(theta) = q(cos-1) + r sin` | `kill_k8e`: identity fails at 1e-8 anywhere, or predicted culprit hole right on < 40% |
ideas/round2_adversary.md:826: | 6 curvature budget | 2025 Fig. F.4 documented failure; `periodic_jacobian.hpp` `J(theta) = cI + sK` | `kill_k8f`: >= 2 mispredictions on the 16x4 grid, or `det J(Theta_max) != 1 + A/a` |
ideas/round2_adversary.md:827: | 7 disconnected `U(eps)` | F28 9.47% silent re-closure; `check.md` T5.2c "not basic" unproved | `kill_k8g`: every certified pair joined by a certified segment on all three patterns |
ideas/round2_adversary.md:828: | 8 certified non-uniform continuation | U5 mobility constant on `(0, theta_max)`; T2.4 every stop is a contact | `kill_k8h`: median `theta_eff < 0.05` rad on 20 graphs |
ideas/round2_adversary.md:829: | 9 certified inverse-design oracle | F24 their `merge_close_verts` reports 0.067 vs true ~1.65 on `snub_square`; K2a exact on 187/187 | `kill_k8i`: < 3 of 16 patterns disagree materially |
ideas/round2_adversary.md:830: | 10 permanent incidences | `KILL_REPORT.md` correction 2: 13 630 spurious roots deflated; `check.md` 143/1.82M third-class | `kill_k8j`: counts disagree, or the count is positive on random graphs |
ideas/round2_adversary.md:872: expansion in `zero_plus.hpp`, not an empirical observation. The "25 of 50" objection is right that
ideas/round2_adversary.md:930: **If I could run exactly three things tomorrow:** `kill_k8d` (idea 4, one hour, decides the framing of
ideas/round2_adversary.md:931: everything already measured), then `kill_k8a` (idea 1, ten minutes, decides whether the project has a
ideas/round2_adversary.md:932: positive half), then `kill_k8b` (idea 2, twenty minutes, decides whether the negative half is a
```

### `ideas/round2_inverse.md` (23)

_HISTORICAL — `code/src/method/*.hpp`, `code/apps/kill_k7.cpp`; proposed `kill_i1..i10.cpp` drivers; '≈ 90 lines of C++' effort notes; `code/src/export`._

```
ideas/round2_inverse.md:11: K1a, K2a headline numbers), `code/src/method/periodic_jacobian.hpp` (the derived comment block),
ideas/round2_inverse.md:12: `code/src/method/contact.hpp`, `code/apps/kill_k7.cpp`, and the **measured K7 output**
ideas/round2_inverse.md:65: ## P0 — shared infrastructure every kill test below needs (≈ 90 lines of C++, one hour, once)
ideas/round2_inverse.md:75: plus `detect_lattice` (already in `method/periodic_jacobian.hpp`) to recover `t1, t2`, and the
ideas/round2_inverse.md:77: 2025 Remark 4.1; our `orientation.hpp` max-cut already returns it when the dual is bipartite).
ideas/round2_inverse.md:81: `kill_k7.cpp`) plus a hand-authored `hexagon_triangle`, and the tests below still run, with the
ideas/round2_inverse.md:116: `A(θ) = det(P₀)·(det J(θ) − 1)` with `J(θ) = cos(θ/2) I + sin(θ/2) K` (`periodic_jacobian.hpp`).
ideas/round2_inverse.md:130: **Kill experiment.** Driver `kill_i1.cpp` on top of P0. (a) Sample `X ∈ 𝕏` on a free-boundary
ideas/round2_inverse.md:164: `AchievableSet` in `kill_k7.cpp`. A spherical cap of total Gaussian curvature `Ω` is conformally
ideas/round2_inverse.md:198: **Kill experiment.** `kill_i2.cpp`. (a) For each of the 7 families in the K7 population plus
ideas/round2_inverse.md:266: **Kill experiment.** `kill_i3.cpp`, no P0 needed. Generate `m×m` patches, `m = 2 … 20`, for all 7
ideas/round2_inverse.md:322: **Kill experiment.** `kill_i4.cpp` (or three lines added to `kill_k7.cpp`, which already computes
ideas/round2_inverse.md:369: **Kill experiment.** `kill_i5.cpp`. On the 4 target meshes shipped with their repo
ideas/round2_inverse.md:389: C++; this is the real cost and it may not be worth paying).
ideas/round2_inverse.md:418: **Kill experiment.** `kill_i6.cpp`. On a 2025 pattern patch (a quad family, ~400 faces) with a free
ideas/round2_inverse.md:429: SVG whose deployment is certified rather than simulated. Fabricable export: `code/src/export`
ideas/round2_inverse.md:468: **Kill experiment.** `kill_i7.cpp`. Take a family of spherical caps of increasing `θ₀` and, for one
ideas/round2_inverse.md:495: `span{I, J_rot}` (proved in `periodic_jacobian.hpp`, verified by K7 C2 to `c2_conf_res ≈ 1e−17` on
ideas/round2_inverse.md:513: **Kill experiment.** `kill_i8.cpp`. For the 16 patterns, sample `𝒦` (32 random `t`) × 64 angles and
ideas/round2_inverse.md:555: **Kill experiment.** `kill_i9.cpp`, after P0: run `exact_theta_max` (full candidate list, no
ideas/round2_inverse.md:601: **Kill experiment.** `kill_i10.cpp`. Build `DG(0)` on (a) the 16 imported patterns, (b) a
ideas/round2_inverse.md:695:   their pattern JSON format comes from the files themselves, and my statements about their C++
ideas/round2_inverse.md:697:   `code/src/method/periodic_jacobian.hpp` — **not** from reading their C++.
```

### `ideas/round2_theorist.md` (17)

_HISTORICAL — `code/src/method/*.hpp` refs; proposed extensions of `kill_k1a/k5/k6.cpp`, new `kill_k8/k9/k10.cpp`; `mobility.cpp`, `generators.cpp`, `mesh.cpp`._

```
ideas/round2_theorist.md:12: `k6_full.csv` headers + first rows, `code/src/method/{deploy_basis.hpp, zero_plus.hpp}`,
ideas/round2_theorist.md:34: `code/src/method/zero_plus.hpp` computes `q_e = det(dS_e, d_e)`, a positive multiple of `r`. **The
ideas/round2_theorist.md:103: **6. Kill experiment.** Driver: extend `code/apps/kill_k1a.cpp`'s population loop (200 graphs,
ideas/round2_theorist.md:104: `F ∈ [101, 793]`) plus the 8 authored tilings of `kill_common.hpp::deployable_population()`. For each
ideas/round2_theorist.md:137: (`zero_plus.hpp`). Emptiness of `{t : q_e(t) > 0 ∀e, μ_c(t) > 0 ∀c}` is *certified* by any
ideas/round2_theorist.md:165: **6. Kill experiment.** Driver: new `kill_k8.cpp` on top of `zero_plus_form()` (already returns `GS`,
ideas/round2_theorist.md:184: in `zero_plus.hpp`), so the aggregate over corner incidences is not a single quadratic — the
ideas/round2_theorist.md:235: compare with the right-hand side assembled from `cut.hpp` edge lists and `zero_plus`'s `Δu`.
ideas/round2_theorist.md:295: **6. Kill experiment.** Driver: extend `kill_k1a.cpp`. On its 200 graphs compute, for `X₀`, the
ideas/round2_theorist.md:355: **6. Kill experiment.** Driver: a 60-line addition to `kill_k5.cpp` (it already builds `σ_mc`,
ideas/round2_theorist.md:413: **6. Kill experiment.** Driver: new `kill_k9.cpp` reusing `mobility.cpp`'s cycle-basis assembly. On the
ideas/round2_theorist.md:465: per `zero_plus.hpp`'s convention note, flip both `d_e` and `dS_e`, leaving `q` unchanged). Restricting
ideas/round2_theorist.md:505: (`zero_plus.hpp`). The corner obstruction is the binding one and neither paper has any account of it —
ideas/round2_theorist.md:516: `sin(θ/2) (S_p − S_a)` with `S_p − S_a = J[2Δu − (σ_p − σ_f) x_v]` (`zero_plus.hpp` corner section). For
ideas/round2_theorist.md:524: **6. Kill experiment.** Driver: extend `kill_k6.cpp` (it already computes `corner_incidences()` and
ideas/round2_theorist.md:574: **6. Kill experiment.** Driver: new `kill_k10.cpp` using `generators.cpp` for input and `mesh.cpp` for
ideas/round2_theorist.md:736:   of `code/src/method/zero_plus.hpp`, not on a proof.
```

### `ideas/round2_theorist_b.md` (36)

_HISTORICAL — `code/src/method/{zero_plus,periodic_jacobian}.hpp`, scratch `check_b1.cpp`, proposed `kill_b1..b9.cpp` drivers, `kill_k6.cpp`/`kill_k7.cpp`._

```
ideas/round2_theorist_b.md:7: `code/src/method/zero_plus.hpp` and `code/src/method/periodic_jacobian.hpp` in full,
ideas/round2_theorist_b.md:11: `derivations/scratch/check_b1.cpp`, which verifies the central identity of this file (§0.g); it links
ideas/round2_theorist_b.md:12: against `code/src/core` and modifies nothing under `code/`. I did **not** re-read the 2026/2025 paper text; every paper citation is
ideas/round2_theorist_b.md:84: `q_e = det(dS_e, d_e)` that `code/src/method/zero_plus.hpp` computes and that K5/K6 measured. So:
ideas/round2_theorist_b.md:110: is *not single-valued* on the quotient: `code/src/method/periodic_jacobian.hpp` derives
ideas/round2_theorist_b.md:128: `kill_k7` — give the same two coefficients is the strongest evidence I have that (B.1) is correct
ideas/round2_theorist_b.md:166:   else; by T5.3 and `zero_plus.hpp` this is a positive multiple of the measured `q_e`. It is the only
ideas/round2_theorist_b.md:178: `derivations/scratch/check_b1.cpp` (build line in its header; links only against `code/src/core/*.cpp`,
ideas/round2_theorist_b.md:227: `periodic_jacobian.hpp` has the *periodic total* `det(P₀)(det J − 1)`; it does not have the per-hole
ideas/round2_theorist_b.md:234: The periodic case (0.d) reproduces `periodic_jacobian.hpp`'s independently derived coefficients,
ideas/round2_theorist_b.md:239: is a new app `kill_b1.cpp` on the existing 16-graph corpus plus 30 random
ideas/round2_theorist_b.md:243: `Γ` (30 lines, the same BFS as `derivations/scratch/check_t1_t2.cpp` C1); (iii) compute the
ideas/round2_theorist_b.md:299: **6. Kill experiment.** `kill_b2.cpp` reusing `kill_k6.cpp`'s design population verbatim (the CSVs
ideas/round2_theorist_b.md:328: `AchievableSet` in `kill_k7.cpp` already computes, `dimK = 4` on every K7 row) enters the
ideas/round2_theorist_b.md:356: gives the threshold. Both `tr K` and `W` are already computable in `kill_k7.cpp` (the `c4_r` column
ideas/round2_theorist_b.md:359: **6. Kill experiment.** `kill_b3.cpp`, reusing `kill_k7.cpp`'s pattern list and `AchievableSet`
ideas/round2_theorist_b.md:417: **6. Kill experiment.** `kill_b4.cpp` = `kill_k6.cpp` with one change: build `[L]` without the
ideas/round2_theorist_b.md:442: `t` (`zero_plus.hpp` says so and `check.md` T5-a verifies it to `3.4e−14`), so on the homogenised
ideas/round2_theorist_b.md:466: **6. Kill experiment.** `kill_b5.cpp` on the identical K6 400-design population, 200 MW iterations,
ideas/round2_theorist_b.md:511: **6. Kill experiment.** `kill_b6.cpp`: on 60 Delaunay/Voronoi patches with `|F| ∈ [100, 3000]`
ideas/round2_theorist_b.md:560: **6. Kill experiment.** `kill_b7.cpp` = `kill_k7.cpp`'s C4 block over 60 Voronoi tori (5 seeds ×
ideas/round2_theorist_b.md:608: **6. Kill experiment.** `kill_b8.cpp`: the seven tiling families already in `make_tiling_pattern`
ideas/round2_theorist_b.md:668: **6. Kill experiment.** `kill_b9.cpp` on 20 patches: build `A(X, 0)` with the existing
ideas/round2_theorist_b.md:669: `method/mobility.hpp`, draw 50 random `V ∈ ker A` per patch, compute the interior sum
ideas/round2_theorist_b.md:713: has `A(θ) = q(cos θ − 1) + r sin θ` and `2 arctan(r/q)` for the *pattern*; `periodic_jacobian.hpp` has
ideas/round2_theorist_b.md:755: neighbouring corner cone, which is exactly the predicate `zero_plus.hpp` calls `mu` and computes as
ideas/round2_theorist_b.md:840: | 1 | **B3** `tr K` is the budget; budget half-space in `K`-space | thm + alg | the only pair in this file that is a *constructive* deliverable: certified deployable designs on arbitrary **periodic** graphs, where 2026 shows only authore...
ideas/round2_theorist_b.md:841: | 2 | **B4** the fixed boundary, not the graph, empties the space | thm + exp | cheapest decisive experiment in the project; if positive it rewrites the meaning of F30 and of two sibling ideas | `kill_b4`: `0/400` again ⟹ dead, ≈ 20 min ...
ideas/round2_theorist_b.md:842: | 3 | **B1** the budget identity `Σ_C a_C = B(X)` | theorem | the identity everything else is derived from; the per-hole half is **already verified** at `8.7e−14` over 471 holes (§0.g) and it reproduces `periodic_jacobian.hpp`'s independ...
ideas/round2_theorist_b.md:843: | 4 | **B2** `O(\|E\|)` budget bound on `min_e q_e` | thm + alg | turns the identity into a microsecond non-deployability certificate; genuinely two-sided, since it may be true and never binding | `kill_b2`: violation ⟹ hard kill; never ...
ideas/round2_theorist_b.md:845: | 6 | **B9** the budget law for every flex; bound on the expansive cone | theorem | the right way to settle whether non-uniform deployment can rescue uniform failure — 2026 §7 limitation (ii) — and the LP half is decisive either way | `k...
ideas/round2_theorist_b.md:846: | 7 | **B6** notch leakage and the `\|F\|/\|∂M\|` scaling law | thm + law | correct accounting and a plausible mechanism for the patch-vs-torus gap, but likely to produce a curve that is flat at zero | `kill_b6`: no monotone trend ⟹ dead...
ideas/round2_theorist_b.md:847: | 8 | **B5** spectral budget allocation | algorithm | pairs with adversary Idea 2 into "a design or a proof there is none"; but if B4 shows the fixed-boundary set really is empty, this returns certificates and no designs | `kill_b5`: `0/...
ideas/round2_theorist_b.md:848: | 9 | **B8** duality: `θ_c(P) + θ_c(P*) = 2π` | theorem | exact on the one dual pair measured and self-consistent on the self-dual case, but three families is not a theorem and it is the most likely to be known already | `kill_b8`: any p...
ideas/round2_theorist_b.md:849: | 10 | **B7** universality of `K` for random periodic graphs | charac. | the most surprising thing in the K7 data (`θ_c = 1.926 ± 0.05` over a 3× size range) and the least theorem-like; `n = 12`, one seed per size | `kill_b7`: no `\|F\|^...
ideas/round2_theorist_b.md:858: * **(B.1) and (B.2s) are measured** (§0.g, `derivations/scratch/check_b1.cpp`, 471 holes / 58 split
```

### `specs/builder_baseline.md` (6)

_HISTORICAL spec — KEEP (baseline/native CMake project, clang, bind.cpp); refers to `code/README.md` JSON contract._

```
specs/builder_baseline.md:3: Output: baseline/native/ (CMake project wrapping baseline/tuttekiri/code/cpp), a CLI `tuttekiri_cli`, and baseline/README.md. Do NOT modify anything inside baseline/tuttekiri/ (treat it as read-only upstream); copy or reference its sourc...
specs/builder_baseline.md:5: Read first: specs/common_preamble.md; code/README.md (our JSON graph contract); baseline/tuttekiri/README.md; baseline/tuttekiri/code/cpp/ (all files: utils/Hmesh, coloring/*, geometry/{kirigami,deployment,unit_pattern}, opt/{opt,deploym...
specs/builder_baseline.md:7: Goal: a native (arm64 macOS, clang) executable that runs the AUTHORS' implementations of: face-orientation assignment (coloring/*), hole/constraint construction + Tutte auxetic embedding solve (their Eq. 4/6 equivalents), shape-space bas...
specs/builder_baseline.md:10: 1. Dependencies: libigl v2.5.0 and Optiz (github.com/segaviv/optiz) are pulled by their CMake via FetchContent — cloning public repos is allowed (nothing leaves the machine). Eigen: prefer /opt/homebrew/include/eigen3 (Eigen 5.0.1) if it...
specs/builder_baseline.md:11: 2. Mirror bind.cpp's exported functions: read what the web UI calls (sequence of calls for load → color → solve → basis → deploy → prevent intersections → fully closed) and replicate that sequence in the CLI with subcommands: `color`, `s...
specs/builder_baseline.md:13: 4. Document exact build commands. If something cannot be built natively (e.g. a hard emscripten dependency), isolate the smallest subset that can (at minimum: coloring + constraint solve + deployment) and say precisely what is missing.
```

### `specs/builder_core.md` (15)

_HISTORICAL spec — title 'C++ reimplementation'; CMake targets `kiri_core`/`kiri_tests`; Eigen/nlohmann/doctest deps; `*.hpp` module names; build command._

```
specs/builder_core.md:1: # Builder-Core: C++ reimplementation of the Segall 2026 pipeline
specs/builder_core.md:3: Output: code/ (CMake project), all tests passing, plus code/README.md documenting the API and the CLI tools.
specs/builder_core.md:10: - code/CMakeLists.txt (C++20, `-O2`, `CMAKE_OSX_ARCHITECTURES=arm64`, warnings on). Targets: `kiri_core` (static lib), `kiri_tests` (doctest), CLI apps under code/apps/.
specs/builder_core.md:11: - Headers under code/src/core/, one concern per file. Namespace `kiri`.
specs/builder_core.md:12: - Dependencies: Eigen (dense+sparse), nlohmann/json, doctest — all header-only in /opt/homebrew/include.
specs/builder_core.md:13: - Build command must be documented and verified: `cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j && ./code/build/kiri_tests`.
specs/builder_core.md:15: ### 1. Graph I/O and data model (`mesh.hpp`)
specs/builder_core.md:27: ### 2. Edge classification and cutting (`cut.hpp`)
specs/builder_core.md:32: ### 3. Hole preimage detection (Algorithm 1, 2026 Sec 4.1) (`holes.hpp`)
specs/builder_core.md:38: ### 4. Deployability test and the linear system (Eq. 2–6) (`tutte_auxetic.hpp`)
specs/builder_core.md:42: - Solve: dense for N ≤ ~5000 unknowns (Eigen ColPivHouseholderQR / BDCSVD), sparse SparseQR otherwise. Provide: rank, null-space basis {φ_i} (columns, N×k) via SVD of the constraint matrix, particular solution X0 = argmin ‖X − X_ini‖ s.t...
specs/builder_core.md:45: ### 5. Forward kinematics (`kinematics.hpp`)
specs/builder_core.md:58: ### 8. Test-graph generators (`generators.hpp`) and CLI apps
specs/builder_core.md:63: ### 9. Tests (doctest, code/tests/) — all must pass; write failing-first where possible
specs/builder_core.md:82: 1. `cmake --build` succeeds with zero warnings-as-errors disabled but no errors; `kiri_tests` reports all passed.
```

### `specs/builder_export.md` (3)

_HISTORICAL spec — `code/src/export/`, `code/apps/kiri_export.cpp`, `code/tests/test_export.cpp`, doctest, CMakeLists, build command._

```
specs/builder_export.md:3: Output: code/src/export/ (namespace kiri::export_), app code/apps/kiri_export.cpp, tests in code/tests/test_export.cpp, docs appended to code/README.md (section "Export"). All C++ (D3). Another agent is concurrently editing code/apps and...
specs/builder_export.md:12: 5. Tests (doctest, standalone-buildable like the others; register in your CMakeLists): SVG parses as XML (write a tiny well-formedness check), path count equals face count, neck present at every hinge vertex and absent at split edges, ST...
specs/builder_export.md:14: Build: `cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build code/build -j && ./code/build/kiri_tests` must still pass, plus your own test binary. Do not commit. Reply with the file list, test results, and the sli...
```

### `specs/builder_method.md` (4)

_HISTORICAL spec — `code/src/method/`, `collision.hpp`, doctest/CMakeLists._

```
specs/builder_method.md:10: Output: code/src/method/ (namespace kiri::method; extend what the Experimenter started: deploy_basis, mobility, contact, range_opt), app code/apps/kiri_usable.cpp, tests code/tests/test_method.cpp (extend), README section "Method". All C...
specs/builder_method.md:18: 4. CERTIFICATE: (i) all face signed areas > 0 at θ=0; (ii) no face–face overlap at θ=0 (collision.hpp exact test); (iii) exact θ_max recomputed with the COMPLETE candidate list (no pruning) and cross-checked against bisection (tol 1e-5);...
specs/builder_method.md:19: 5. Baselines for comparison, callable from the same CLI: (a) Eq.(1) σ + Eq.(6) only; (b) (a) + our collision.hpp Eq.(9); (c) (a) + the authors' native `prevent` via baseline/native/build/tuttekiri_cli (subprocess; parse its JSON; then re...
specs/builder_method.md:21: Tests (doctest, register in CMakeLists): certificate agrees with bisection on 100 random valid embeddings; the σ search never returns c(Γ) > 1; on the 8 reference cases the method's θ_max ≥ the authors' prevent result (or document where ...
```

### `specs/common_preamble.md` (2)

_HISTORICAL — directive D3 'ALL code is C++17/20 with CMake' (l.8) and toolchain line clang++/cmake/Eigen/nlohmann/doctest versions (l.9)._

```
specs/common_preamble.md:8: - Coding directive (D3): ALL code is C++17/20 with CMake. Python is allowed ONLY for matplotlib plotting of CSV/JSON dumped by C++ programs, and must be run as `arch -arm64 /usr/local/bin/python3` (the shell runs under Rosetta; the numpy...
specs/common_preamble.md:9: - Toolchain (verified): clang++ (Apple), cmake 3.29, Eigen 5.0.1, nlohmann-json 3.12, doctest 2.5.3, all headers in /opt/homebrew/include. Build with `-arch arm64` / `CMAKE_OSX_ARCHITECTURES=arm64`.
```

### `specs/critic.md` (2)

_HISTORICAL — 'demonstrable in C++', `code/src/core`._

```
specs/critic.md:11: C4 demonstrable in C++ within the run on random + authored planar graphs, 100–5,000 faces, on top of code/src/core.
specs/critic.md:20: 4. For each of the top 5, restate the kill experiment as an exact spec that an Experimenter can run against code/src/core with no further design decisions: generator(s), number of graphs, σ choice, routine to call, quantity to compute, t...
```

### `specs/deriver.md` (1)

_HISTORICAL — 'throwaway numeric checks in C++ under derivations/scratch/' (l.3)._

```
specs/deriver.md:3: Output: derivations/core.md (step-numbered; every step either a definition, a cited fact, or an algebraic step that a Checker can verify by hand or by a numeric test). No code, except that you MAY write small throwaway numeric checks in ...
```

### `specs/experimenter_final.md` (2)

_HISTORICAL — 'all C++', `collision.hpp`._

```
specs/experimenter_final.md:6: Experiments (all C++; ≥ 20 graphs everywhere, ≥ 500 for the characterization validation):
specs/experimenter_final.md:7: E1 Characterization validation (MISSION §8b): exact Θ_max (T4.2″) vs brute-force bisection (collision.hpp, tol 1e-5) on ≥ 500 random graphs (Voronoi/Delaunay/quad-random, 100–5000 faces, σ from Eq.(1) AND σ_def if K5 passed) at ≥ 3 embed...
```

### `specs/experimenter_k7.md` (3)

_HISTORICAL — 'All C++ (D3)', `code/src/core/generators`._

```
specs/experimenter_k7.md:3: Output: results/kill/k7/ + a K7 section in results/kill/KILL_REPORT.md. All C++ (D3); sharded across cores from the start (12 processes) if any step exceeds ~5 min serial.
specs/experimenter_k7.md:6: Population: ≥ 30 periodic patterns: the periodic generators in code/src/core/generators (torus/periodic squares, triangles, hexagons, kagome, snub square, 3.4.3.12, 4.8.8 at 2–3 cell sizes) plus periodic random patterns (Voronoi of rando...
specs/experimenter_k7.md:11: C3 (designable Poisson family): Poisson ratio ν(θ) along direction d from U = polar stretch of J(θ) (as in IsoGami Sec 4.3 / 2026 Sec 5.1); the achievable ν-curves are parameterized by K ∈ 𝒦. Algorithm: given a target K* (random in 𝒦, an...
```

### `specs/experimenter_kill.md` (3)

_HISTORICAL — 'small C++ apps under code/apps/' (l.3), `collision.hpp` bisection (l.14), 'add doctest cases' (l.17)._

```
specs/experimenter_kill.md:3: Output: results/kill/ (one subdirectory per experiment with CSV/JSON dumps, PNG figures via `arch -arm64 /usr/local/bin/python3` + matplotlib, and a KILL_REPORT.md at results/kill/KILL_REPORT.md), plus any small C++ apps under code/apps/...
specs/experimenter_kill.md:14: For K2b implement the closed-form first-contact roots properly: vertex-into-edge-interior contact = harmonic root PLUS the two interval (projection) inequalities, each also a harmonic; enumerate candidate (vertex, edge) pairs with the sw...
specs/experimenter_kill.md:17: All existing tests must still pass after your additions; add doctest cases for any new numerical routine (closed-form roots vs bisection, harmonic fit). Do not commit. Reply with the verdict table (experiment, PASS/FAIL, key number) and ...
```

### `specs/ideator.md` (1)

_HISTORICAL — 'Demonstrable in C++'._

```
specs/ideator.md:16: - Demonstrable in C++ within the run on random and hand-authored planar graphs with 100–5,000 faces, on top of the reference pipeline (mesh, cut, hole preimages, linear system + null space, forward kinematics, θ_max by collision bisectio...
```

### `specs/ideator_round2.md` (2)

_HISTORICAL — 'AVAILABLE in code/src/method' (l.6); 'existing drivers (kill_k*.cpp, kiri_analyze)' (l.15)._

```
specs/ideator_round2.md:6: 1. Uniform deployment is trig-linear: Y_θ = cos(θ/2)C(X) + sin(θ/2)S(X); every contact predicate is p + q cos θ + r sin θ with p,q,r quadratic in X; exact Θ_max and a certificate POS ∧ NOOVERLAP(ε/2) ∧ NOROOT are implemented and validate...
specs/ideator_round2.md:15: Ideas that turn the above into a THEOREM, CHARACTERIZATION, or ALGORITHM (MISSION §8) that could be demonstrated within hours on top of the existing code. Especially: (a) an existence / non-existence criterion for usable embeddings (why ...
```

### `specs/m2_checker_lemmas.md` (5)

_HISTORICAL — `code/src/method/*.hpp`, `code/tests/derivation_tests.cpp` (33 cases / ~89k assertions), `kill_common.hpp`, `code/apps/kill_k7.cpp`._

```
specs/m2_checker_lemmas.md:12: 4. code/src/method/deploy_basis.hpp, contact.hpp (`harmonic_roots`, the class split, the identity test), zero_plus.hpp (header derivation of `dS_e`, `dC_e = 0`), periodic_jacobian.hpp (header derivation of `J(θ)`, `K = Q P₀⁻¹`, `Achievab...
specs/m2_checker_lemmas.md:13: 5. code/tests/derivation_tests.cpp — the existing standalone test file (33 cases / ~89k assertions); its build line is documented in derivations/check.md. Your new tests go here.
specs/m2_checker_lemmas.md:22: **Tests.** Add to code/tests/derivation_tests.cpp (keep the file standalone as it is; follow its existing style and the build line in check.md):
specs/m2_checker_lemmas.md:23: - L1: on ≥ 1,000 random `(graph, σ, t)` triples across kill_common.hpp's population (≥ 50 graphs) and the reference tilings, for every candidate pair: predicted class from the combinatorial type per L1.1/L1.3 vs the numeric class from th...
specs/m2_checker_lemmas.md:24: - L2: `rank(A) == 2·rank(D)` on the 33 K7 patterns and ≥ 100 random Voronoi tori with random σ (find the generator in code/apps/kill_k7.cpp), dense SVD, tolerance 1e-10 relative; the reconstruction of each `M_j` from row `j` of `D` by th...
```

### `specs/m2_deriver_lemmas.md` (7)

_HISTORICAL — 'scratch C++ checks', `code/src/method/*.hpp`, `derivations/scratch/check_l1/l2.cpp` compiled against `code/src`; `code/src`/`code/tests` non-modification._

```
specs/m2_deriver_lemmas.md:9: You write `derivations/lemmas.md` with full statements and proofs, plus scratch C++ checks. A separate Checker, in a fresh context, will see only your statements and the code and will try to break them.
specs/m2_deriver_lemmas.md:14: 3. code/src/method/deploy_basis.hpp (the C, S basis), code/src/method/contact.hpp (`harmonic_roots`, the three-class deflation, the identity test) and code/src/method/zero_plus.hpp (the header comment derives `dS_e`, `dC_e = 0` for split...
specs/m2_deriver_lemmas.md:15: 4. results/kill/KILL_REPORT.md §K7 (lines ~1180–1465) and code/src/method/periodic_jacobian.hpp (header comment derives `J(θ) = cos(θ/2)I + sin(θ/2)K`, `K = Q P₀⁻¹`, `Q = 2 Jrot [w_h w_v]`, the `AchievableSet` struct with `A` (4 × 2k), `...
specs/m2_deriver_lemmas.md:28: **Scratch checks you run yourself** (in derivations/scratch/check_l1.cpp, compiled against code/src as the existing scratch files are — see the header of derivations/scratch/check_r3.cpp for the build line): on ≥ 50 random graphs of code...
specs/m2_deriver_lemmas.md:35: - Set up from periodic_jacobian.hpp's derivation: `u_{f+t} = u_f + w_t + σ_f t/2`, `Q = 2 Jrot [w_h w_v]`, `K = Q P₀⁻¹`. The null space `Φ` acts on the x and y coordinates separately (t ∈ ℝ^{2k}). Show that the map `t ↦ K(t) − K0` is lin...
specs/m2_deriver_lemmas.md:40: **Scratch checks** (derivations/scratch/check_l2.cpp): rebuild `A` and `D` via the existing `periodic_jacobian` API on the 33 K7 patterns plus ≥ 100 random Voronoi tori (use `make_tiling_pattern` and the Voronoi torus generator used by k...
specs/m2_deriver_lemmas.md:50: - Do not modify derivations/core.md, code/src, or code/tests. Scratch programs go in derivations/scratch/.
```

### `specs/m2_experimenter_basin.md` (5)

_HISTORICAL — Xcode make workaround, `code/src/method/design.hpp`, `code/apps/kill_yield.cpp`/`kill_basin.cpp`._

```
specs/m2_experimenter_basin.md:4: Environment: build with `/Applications/Xcode.app/Contents/Developer/usr/bin/make -C code/build -j2 <target>` after `cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64` (bare `make`/`cmake --build` are broken by an xcode-select s...
specs/m2_experimenter_basin.md:11: - code/src/method/design.hpp `design_range_max` / `RangeMaxOptions` (seed = 9300 + 7·id + which in the run; arms; stage settings 6 × 120) and code/apps/kill_k9c.cpp (how graph, σ_mc, σ_def and X_ini are produced per design — reuse exactl...
specs/m2_experimenter_basin.md:12: - code/apps/kill_yield.cpp (WP2's driver regenerates the same designs; copy its design-regeneration code path).
specs/m2_experimenter_basin.md:15: 1. Driver code/apps/kill_basin.cpp. Designs: all 93 K9c rows with `theta_exact == 0` (or ≤ 1e-9), plus 30 successes drawn deterministically (every 10th success by row order). For each design run `design_range_max` with the run defaults b...
specs/m2_experimenter_basin.md:20: - Do not modify code/src. Do not touch results/kill/*. Do not commit; do not edit REPORT.md/IDEA.md/STATE.md.
```

### `specs/m2_experimenter_native200.md` (6)

_HISTORICAL — `code/apps/kill_native200.cpp`, `kill_common.hpp`, `native200_merge.cpp`, baseline/native CLI (KEEP)._

```
specs/m2_experimenter_native200.md:10: - code/apps/kill_native200.cpp in full (driver; options --n --maxf --out --sigma --cli --work --timeout --shard --nshards --logdir).
specs/m2_experimenter_native200.md:11: - code/apps/kill_common.hpp (population is a deterministic function of graph id; kind = voronoi/delaunay/quad_random by id % 3).
specs/m2_experimenter_native200.md:14: - baseline/native/ (the vendored authors' code and CMakeLists; the CLI is baseline/native/build/tuttekiri_cli).
specs/m2_experimenter_native200.md:17: 1. **Make the driver resumable and cell-selectable.** Add to kill_native200.cpp: `--cells <file>` (a list of `id,variant` lines to run, ignoring the shard split), `--resume <csv>` (skip any cell whose status in that CSV is `completed` or...
specs/m2_experimenter_native200.md:21: 5. **Merge + report tooling.** Write code/apps/native200_merge.cpp (or extend kill_native200 with `--merge`) that combines native200.csv and every shard CSV under rerun3600/ and crashfix3600/ into results/kill/native200/native200_final.c...
specs/m2_experimenter_native200.md:26:    The report must be regenerable: NATIVE200_FINAL.md is produced by a script (C++ writer or a ≤ 100-line Python that only reads the CSV) so the orchestrator can rerun it as the detached job progresses. Generate a first version now, labe...
```

### `specs/m2_experimenter_regime.md` (9)

_HISTORICAL — xcode-select/cmake environment note, `code/apps/kill_k9c.cpp`, `kill_common.hpp`, `code/src/method/design.hpp`, `code/apps/kill_regime.cpp`, `kill_scaling.cpp`, CMakeLists._

```
specs/m2_experimenter_regime.md:5: Environment (verified 2026-09-08): `cmake --build` and bare `make` are broken by an xcode-select shim; build with `/Applications/Xcode.app/Contents/Developer/usr/bin/make -C code/build -j4 <target>` after `cmake -S code -B code/build -DC...
specs/m2_experimenter_regime.md:12: - code/apps/kill_k9c.cpp (population `make_graph(id, 100, maxf, 1400)` from code/apps/kill_common.hpp; σ_mc from `orientation_maxcut`, σ_def read from results/kill/k5/sigma — for NEW graphs you compute σ_def with `method::orientation_def...
specs/m2_experimenter_regime.md:13: - code/src/method/design.hpp (`characterize`, `design_range_max`, `design_baseline`, `orientation_maxcut`, `orientation_defect`), code/src/core/tutte_auxetic.hpp (`solve_system` dense SVD; `rank_only_sparse`), code/src/core/generators.hp...
specs/m2_experimenter_regime.md:14: - baseline/kirigami_tessellations/README.md and fabrication_patterns/ (the 2025 paper's 16 pattern files; STATE.md P4b notes they are "polygon soup" needing a ~90-line welding importer — check whether one exists in code/ already: grep fo...
specs/m2_experimenter_regime.md:17: 1. New driver code/apps/kill_regime.cpp (add to CMakeLists app list). Population: ids 2000..2099 through `make_graph(id, 20, 100, 1400)` (three families by id % 3; confirm `F ∈ [20, 100]` and report the distribution). Both σ rules per gr...
specs/m2_experimenter_regime.md:20:    - **native** = the authors' full pipeline via the crash-fixed CLI, their own colouring AND our σ (two cells), 600 s cap, scored by our exact scan + referee + ε = 0.3 certificate exactly as kill_native200 does (reuse its per-cell code ...
specs/m2_experimenter_regime.md:21:    Also run the same three arms on the 8 reference tilings of `kiri_reference` (results/core_validation/reference_cases.json lists them) and on every pattern of baseline/kirigami_tessellations that an importer can load (if no importer ex...
specs/m2_experimenter_regime.md:26: 1. Driver code/apps/kill_scaling.cpp. Graphs: Delaunay and Voronoi, target face counts 50, 100, 200, 500, 1,000, 2,000, 5,000 (use `generate` directly with the site counts needed; 3 seeds each; report actual `F`, `N`). For each: wall tim...
specs/m2_experimenter_regime.md:32: - Do not modify code/src/method or code/src/core except to add the importer. Do not touch results/kill/*. Do not commit; do not edit REPORT.md/IDEA.md/STATE.md.
```

### `specs/m2_experimenter_yield.md` (6)

_HISTORICAL — `kill_common.hpp`, `code/src/method/design.hpp`, `code/src/core/*.hpp`, `code/apps/kill_yield.cpp` + CMakeLists._

```
specs/m2_experimenter_yield.md:11: - code/apps/kill_common.hpp (`make_graph(id, min_faces, max_faces, n_cap)`: population is a deterministic function of id; kind by id % 3; the K9 population is ids 0..199 with min 100, max 800, n_cap per kill_k9c.cpp — read that file's ca...
specs/m2_experimenter_yield.md:12: - code/src/method/design.hpp — `design_range_max` (the K9c pipeline as a library call; defaults reproduce the run), `characterize`, `orientation_maxcut`, `orientation_defect`.
specs/m2_experimenter_yield.md:13: - code/src/core/{cut.hpp, holes.hpp, tutte_auxetic.hpp} for the input-side quantities; code/src/method/zero_plus.hpp (`zero_plus_q`, `zero_plus_form`) and convex_embed.hpp for the 0⁺ constraint counts.
specs/m2_experimenter_yield.md:17: 1. **Input-side features, per design** (computed at `X_ini` and from the cut structure only — nothing that the solver produces). Write code/apps/kill_yield.cpp that, for each of the 400 K9c designs (same ids, same σ_mc / σ_def as kill_k9...
specs/m2_experimenter_yield.md:26: 2. **Fresh graphs for held-out validation.** Generate 100 new graphs with `make_graph(id, 100, 800, n_cap)` for ids 1000..1099 (outside the K9 population; confirm they are not equal to any id 0..199 graph), both σ rules, run `design_rang...
specs/m2_experimenter_yield.md:40: - Do not modify code/src/method or code/src/core. New code goes in code/apps/kill_yield.cpp (add to CMakeLists.txt's app list) and results/yield/analyse.py.
```

### `specs/m2_propagator_a.md` (6)

_HISTORICAL — cmake/xcode environment, `derivation_tests` CMake target, `test_reference_cases.cpp`, `kiri_tests` counts, `code/build` hygiene._

```
specs/m2_propagator_a.md:5: Environment (verified 2026-09-08): `cmake --build` and bare `make` are broken by an xcode-select shim; use `/Applications/Xcode.app/Contents/Developer/usr/bin/make -C code/build -j4 <target>` after `cmake -S code -B code/build -DCMAKE_OS...
specs/m2_propagator_a.md:14: - code/CMakeLists.txt; code/tests/derivation_tests.cpp header and the build line for it in derivations/check.md ("Programs run" / round sections — grep for `derivation_tests`); code/tests/test_reference_cases.cpp (a one-line stub); resul...
specs/m2_propagator_a.md:19: 1. **derivation_tests CMake target.** Add `derivation_tests` as an executable target in code/CMakeLists.txt reproducing the documented hand build line (same sources, defines, include paths), so `make -C code/build derivation_tests` build...
specs/m2_propagator_a.md:20: 2. **test_reference_cases.cpp.** Replace the stub with real doctests: load results/core_validation/reference_cases.json, and for each of the 8 reference cases regenerate the pattern with the same generator call `kiri_reference` uses (rea...
specs/m2_propagator_a.md:21: 3. Build tree hygiene: delete the stale duplicates `code/build/libkiri_core 2.a`, `libkiri_core 3.a`, `CMakeFiles 2` if present, and the orphan `code/build/dbg_k1c` (no source). Do not touch anything else in code/build.
specs/m2_propagator_a.md:37: ≤15 lines: derivation_tests target result (cases/assertions); test_reference_cases result and kiri_tests counts before/after; files cleaned; number of sentences changed per document; vault articles added/updated and link-check result; an...
```

### `baseline/README.md` (24)

_KEEP — third-party C++ build notes (CMake/Eigen/emscripten stub); only the `code/src/core` path reference (l.5) needs rewording._

```
baseline/README.md:3: A referee baseline: the **authors' own C++**, compiled natively on arm64 macOS with
baseline/README.md:4: Apple clang, driven on **our** JSON graphs and writing **our** JSON conventions, so that
baseline/README.md:5: every number in `code/src/core` can be checked 1:1 against the reference implementation.
baseline/README.md:9: * `baseline/native/` — the CMake project that compiles those sources and adds a CLI.
baseline/README.md:15: cmake -S baseline/native -B baseline/native/build \
baseline/README.md:17: cmake --build baseline/native/build -j8
baseline/README.md:27: | Eigen | `3.4.0` via FetchContent | **not** Homebrew's 5.0.1, see below |
baseline/README.md:31: 1. **Eigen.** Homebrew ships Eigen 5.0.1 in `/opt/homebrew/include/eigen3`. libigl v2.5.0
baseline/README.md:32:    and Optiz are written against the Eigen 3.x internal API, so 5.0.1 is not a candidate;
baseline/README.md:33:    the build pins the same 3.4.0 that the authors' own `CMakeLists.txt` falls back to when
baseline/README.md:34:    `find_package(Eigen3)` fails. This keeps the numerics identical to theirs rather than
baseline/README.md:36: 2. **Optiz and OpenMP.** `optiz/CMakeLists.txt` appends a bare `-fopenmp` on any Apple
baseline/README.md:43: 3. **emscripten.** `utils/Hmesh.h`, `utils/conversions.h` and `geometry/unit_pattern.h`
baseline/README.md:44:    `#include <emscripten/bind.h>`, and `conversions.h` has non-template inline functions
baseline/README.md:45:    that use `emscripten::val`. `baseline/native/emsc_stub/emscripten/bind.h` is a ~40-line
baseline/README.md:49: `bind.cpp` itself is excluded from the build and replaced by `src/main.cpp`.
baseline/README.md:53: Built and exercised: `coloring/` (face-orientation assignment), `geometry/kirigami.cpp`
baseline/README.md:54: (forward deployment, `max_opening_angle`), `geometry/deployment.cpp` (`make_deployable`:
baseline/README.md:56: `geometry/unit_pattern.cpp`, `opt/prevent_intersections.cpp` (their Eq. 9),
baseline/README.md:57: `opt/fully_closed.cpp`, `opt/opt.cpp`, `parameterization/`, `utils/`.
baseline/README.md:61: called. Nothing had to be dropped for a hard emscripten dependency.
baseline/README.md:68: | `bind_port.{h,cpp}` | `do_polygons_intersect` and `max_opening_angle` copied **verbatim** from `bind.cpp` (lines 279-418), plus one diagnostic that reports *which* face pair collides |
baseline/README.md:70: | `main.cpp` | the CLI |
baseline/README.md:125:    `face_colors == 1` branch (`unit_pattern.cpp`) it evaluates
```

### `baseline/parity.md` (4)

_KEEP — describes the authors' C++ baseline; only `code/src/core` path reference (title, l.1/5) needs -> `Kirigami/src/core`._

```
baseline/parity.md:1: # Parity: our `code/src/core` vs. the authors' published implementation
baseline/parity.md:4: authors' own C++ on the eight reference graphs in `results/core_validation/cases/*/M.json`,
baseline/parity.md:102: not worse than theirs on any case measured. Note their routine seeds from `Eigen::Random`
baseline/parity.md:108: 1%-step forward collision scan, i.e. `theta_max` exactly as `bind.cpp` computes it for the
```

### `baseline/tuttekiri/code/README.md` (5)

_KEEP verbatim — upstream authors' README (WebAssembly build); never edit._

```
baseline/tuttekiri/code/README.md:14: ## Compile the WebAssembly part
baseline/tuttekiri/code/README.md:16: If you change anything in the cpp files, you'll have to recompile the c++ code into WebAssembly.
baseline/tuttekiri/code/README.md:17: To build the `cpp/` WebAssembly component, install Emscripten and run:
baseline/tuttekiri/code/README.md:23: emcmake cmake -DCMAKE_BUILD_TYPE=Release ..
baseline/tuttekiri/code/README.md:27: The compilation will generate kirigami_cpp.js & kirigami_cpp.wasm in the build directory and copy these files into the src/ direcory.
```
