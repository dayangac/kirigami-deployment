# Mission 2 / WP6b — Propagator-A: engineering fixes, vault refresh, and the already-decided text changes

Paste of specs/common_preamble.md applies (read it first).

Environment (verified 2026-09-08): `cmake --build` and bare `make` are broken by an xcode-select shim; use `/Applications/Xcode.app/Contents/Developer/usr/bin/make -C code/build -j4 <target>` after `cmake -S code -B code/build -DCMAKE_OSX_ARCHITECTURES=arm64`. The machine runs ~14 native baseline processes plus three agents; keep your own concurrency ≤ 2. Theory and computation only (D15): nothing about fabrication.

## Why this exists
Three things can be done now without waiting for the pending work packages: (1) two engineering defects an artifact reviewer will trip on; (2) the Obsidian vault is stale at K9's 36/400 and lacks everything since; (3) two text changes are already decided by verified facts F41 (Jiang & Choi cleared) and F42 (the Native200 crash cause) and by D15 (no fabrication).

## Read first
- STATE.md: "Mission 2" section in full (F41, F42, D12–D16, E1), plus F24, F40.
- notes/screen_jiang_choi.md §6 (exact replacement sentences).
- results/kill/native200/crashfix.patch and NATIVE200_FINAL.md (partial) for the crash-cause wording.
- code/CMakeLists.txt; code/tests/derivation_tests.cpp header and the build line for it in derivations/check.md ("Programs run" / round sections — grep for `derivation_tests`); code/tests/test_reference_cases.cpp (a one-line stub); results/core_validation/reference_cases.json.
- ~/Documents/ObsidianPT/kirigami-exp/wiki/index.md, sources.md, lint.md, and the article style of e.g. wiki/k9-constrained-embedding.md; the global vault rules in ~/.claude/CLAUDE.md "Knowledge bases" section (one .md per concept, [[wikilinks]], one-line definition + why it matters + links + source citations with page/line refs, index.md and sources.md maintained, incremental, never clobber).

## Task
### A. Engineering
1. **derivation_tests CMake target.** Add `derivation_tests` as an executable target in code/CMakeLists.txt reproducing the documented hand build line (same sources, defines, include paths), so `make -C code/build derivation_tests` builds it. Build it and run it; it must report the same 33 cases / 89,074 assertions / 0 failures (or the current count if the Checker-L has already added cases — report what you see). Do NOT alter the test file.
2. **test_reference_cases.cpp.** Replace the stub with real doctests: load results/core_validation/reference_cases.json, and for each of the 8 reference cases regenerate the pattern with the same generator call `kiri_reference` uses (read code/apps/kiri_reference.cpp), run the pipeline, and assert `F`, `n_split`, `dim_null`, and the exact `Θ_max` from `method::characterize` to 1e-9 against the JSON — **except** snub_square, whose JSON `Θ_max = 1.646` is a Gate-2 artefact from the old collision code (REPORT.md "How each gate was verified", Gate 2); for that case assert the T4 value 1.671 from results/kill/k2a (find the exact number there) and put a comment saying why. Rebuild and run kiri_tests; report counts before and after.
3. Build tree hygiene: delete the stale duplicates `code/build/libkiri_core 2.a`, `libkiri_core 3.a`, `CMakeFiles 2` if present, and the orphan `code/build/dbg_k1c` (no source). Do not touch anything else in code/build.

### B. Text changes already decided
4. In REPORT.md, IDEA.md and docs/techreport/techreport.tex: replace the Jiang & Choi sentences exactly as notes/screen_jiang_choi.md §6 prescribes (REPORT.md objection 15, Related work last sentence; IDEA.md §6 bullet; find the techreport's equivalent by grepping `2608.30032`). Add the Jiang & Choi BibTeX entry to docs/techreport/refs.bib.
5. Correct the crash attribution (F42): REPORT.md "Errata against the published pipeline" item 4 and the Native200 paragraph in "Mandatory caveats"/objection 3, IDEA.md §4 table row, and the techreport §2 upstream-defects sentence and §9 Native200 paragraph currently attribute the 39 crashes to `UnitPattern::get_holes()`; the reproduced cause is an empty-vertex-list access in `convert::to_eig_mat` from `Hmesh::merge_close_verts` after a deployed configuration collapses to zero faces (crashfix.patch). Reword; keep the `get_holes()` null-check as a second latent defect fixed by the same patch. Do NOT yet change any Native200 *count* — those wait for WP5's final report; add nothing that says the rerun is complete.
6. REPORT.md "Next steps toward a paper" item 6 ("Fabricate one hero physically") and objection 16 ("No physical artifact"): reword as out of scope by decision (D15), not as pending work. Same in the techreport §11 open-items list. Leave the export/ artefacts and the sentences that describe the existing SVG/3MF exports as they are (they exist and were verified by xmllint/unzip); remove only forward-looking fabrication promises.
7. List every changed sentence (file, old → new) in results/propagation/CHANGES_A.md.

### C. Vault refresh (incremental, never rebuild)
8. Add articles to ~/Documents/ObsidianPT/kirigami-exp/wiki/ for: k9c-range-maximising (307/400, objective argument, provenance, yield trend), e1-characterization (3,113 designs), native200-baseline (state as of today: 172 completed at 600 s, F42 crash cause, rerun in progress — label as in progress), hero-examples (hero id 148 and hero2 id 130 as *rendered* exemplars), jiang-choi-screen (F41), mission-2 (the work packages and their status from STATE.md). Update k9-constrained-embedding.md and index.md so the "where the project stands" paragraph reflects K9c, not K9; update sources.md, facts.md (F37–F42), decisions.md (D12–D16), dead-ends.md, open-questions.md. Every article: one-line definition, why it matters, [[links]], citations to repo files with line/section refs. Run a link check (every [[target]] resolves to a file) and report orphans/stubs; do not run graphify unless graphify-out/ update is cheap (skip and say so if not).

## Rules
- One logical change per file group; report exact counts and diffs. Do not commit.
- Do not touch results/kill/*, derivations/*, results/yield, results/regime, results/scaling (other agents own them).

## Reply
≤15 lines: derivation_tests target result (cases/assertions); test_reference_cases result and kiri_tests counts before/after; files cleaned; number of sentences changed per document; vault articles added/updated and link-check result; anything not done.
