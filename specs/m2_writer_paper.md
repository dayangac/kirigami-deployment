# Mission 2 / WP6a — Writer-P: the TOG paper skeleton

Paste of specs/common_preamble.md applies (read it first). No code. No fabrication content (D15: theory and computation only).

## Why this exists
MISSION2.md §1 item 6: a `docs/paper/` skeleton in ACM TOG format with the abstract, contribution list and figure list frozen, so that when WP1 (lemmas), WP2 (yield), WP5 (Native200) and WP7 (regime, scaling) land, their numbers drop into named slots rather than being written from scratch. You write the skeleton now from what is verified; you mark every slot that a pending work package will fill with a `\todo{WP-n: ...}` macro so nothing pending is ever stated as done.

## Read first
- REPORT.md in full; IDEA.md §1, §2, §7 (mandatory caveats), §8.
- review/theory_review.md, review/negative_review.md, review/constructive_review.md — their demands are binding on wording.
- docs/poster/poster.tex (the 2-page version: reuse its abstract sentences where they are still true) and docs/techreport/techreport.tex §1 reading guide, §11 (open items) — the techreport's [D]/[F]/[A]/[N] tagging is a working-document device and does NOT go in the paper.
- notes/screen_jiang_choi.md §6 (the exact replacement wording for the Jiang & Choi sentences; F41) and STATE.md "Mission 2 verified facts" (F41, F42) and decisions D12–D16.
- results/final/figures/README.md and summary_table.md; results/kill/KILL_REPORT.md verdict table (lines ~103–125).
- docs/techreport/refs.bib and docs/poster/refs.bib (the only bibliography sources you may cite from; add the Jiang & Choi entry from screen_jiang_choi.md).

## Task
Create docs/paper/ with:
1. `paper.tex` — `\documentclass[acmtog]{acmart}` (check the installed acmart version with `kpsewhich acmart.cls`; the poster already builds with acmart sigconf so the class is present). Title (working): *Exact Deployment Range, Sound Certificates and Range-Maximising Embeddings for Uniformly Deployable Kirigami on Arbitrary Planar Graphs*. Author block: Emre Dayangaç, affiliation `\todo{affiliation}`; a second author slot commented out with `% mentor/co-author TBD`.
   - **Abstract** ≤ 250 words, every number from REPORT.md's abstract or summary_table.md, with `\todo{WP5}` where the Native200 baseline sentence goes and `\todo{WP7}` where the paper-regime sentence goes.
   - **Contributions**: the C1–C10 list with the theory review's demotions applied (T1, T2, T6, K7 C1/C4, F35 as lemmas, not contributions), with the two lemmas of WP1 listed as `\todo{WP1: proved / measured}`.
   - Section headings with a one-paragraph plan each (what goes there, which file supplies it), no body text: 1 Introduction; 2 Related work (with the Jiang & Choi sentence per F41, Choi–Dudte–Mahadevan, Rote–Santos–Streinu, Grima, Acuña, Dang, Konaković-Luković, Tay–Whiteley primary source `\todo{primary citation}`, IsoGami); 3 Background: the 2026 framework (σ, cuts, Eq. (2)–(6)); 4 Exact kinematics and contact calculus (T1–T4, with the lemma slots); 5 The certificate (T5, sound, incomplete, F32 history in one sentence); 6 Why the published projection jams (K1a, T-1 reflex corners, K8a dual certificates with the branch-chart caveat); 7 Range-maximising constrained embedding (K9→K9b→K9c, the objective argument, 307/400, yield trend `\todo{WP2 mechanism}`); 8 Periodic patterns (K7 C1/C2/C4, C3 caveat, dim formula `\todo{WP1}`); 9 Evaluation (E1 3,113; baselines incl. `\todo{WP5}`; paper-regime `\todo{WP7a}`; scaling `\todo{WP7b}`); 10 Limitations (every caveat of IDEA.md §7, minus fabrication); 11 Conclusion.
   - **Figure list** as a LaTeX table in a comment block and as `\begin{figure}` stubs with `\includegraphics` paths to the existing PNGs under results/final/figures, results/kill/k9c, results/final/e1, export/hero2 (the teaser: hero2 closed / half / 0.9 — as *rendered images*, no fabrication), and `\todo` stubs for fig_yield_predictor (WP2), fig_regime (WP7a), fig_scaling (WP7b).
   - A `\todo` macro that prints in red in draft and can be switched off.
2. `refs.bib` — merged from the two existing .bib files plus Jiang & Choi; no entry that is not in one of those sources.
3. `Makefile` or `build.sh` running `latexmk -pdf` (check `which latexmk`; fall back to pdflatex/bibtex ×2). **Build it**; the PDF must compile with zero errors; warnings listed in docs/paper/build.log summary.
4. `docs/paper/SLOTS.md`: a table of every `\todo{WP-n}` in paper.tex with the file that will fill it, so the Propagator can close them mechanically.

## Rules
- No sentence claims anything not in REPORT.md/IDEA.md/STATE.md verified facts. No fabrication, slicer, photo, or physical-artefact wording anywhere (D15).
- Do not edit REPORT.md/IDEA.md/STATE.md or the techreport. Do not commit.

## Reply
≤15 lines: files created; build result (errors/warnings count); number of `\todo` slots and which WP each waits on; abstract word count; anything you could not source.
