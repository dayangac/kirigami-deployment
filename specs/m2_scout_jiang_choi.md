# Mission 2 / WP4 — Scout-JC: full-text screen of arXiv:2608.30032 (Jiang & Choi 2026)

Paste of specs/common_preamble.md applies (read it first). Web search and web fetch ARE allowed for this role.

## Why this exists
REPORT.md §"What a reviewer will object to" item 15 and IDEA.md §6 flag arXiv:2608.30032, "A unified geometric design framework for kirigami structures" (Jiang & Choi, 2026), as an unresolved novelty risk: only the abstract was retrievable in the earlier screen (notes/screen_k9.md items 9–10, notes/screen_r1.md). A TOG submission cannot go out with that open. Your job is to obtain the full text and screen it claim by claim.

## Read first
- REPORT.md §Contributions (the table C1–C10) and §Related work.
- IDEA.md §2 (theorems), §5 (the algorithm), §6 (relation to prior work).
- notes/screen_k9.md and notes/screen_r1.md (what was tried before, and the vocabulary used).
- derivations/core.md headings only (T1–T7 statements), to know what our theory claims.

## Task
1. **Retrieve the full text.** Try, in order, logging every URL and the HTTP result / what came back:
   a. https://arxiv.org/pdf/2608.30032 (and /pdf/2608.30032v1, v2)
   b. https://arxiv.org/html/2608.30032
   c. https://arxiv.org/abs/2608.30032 → any "other formats" / source link
   d. Semantic Scholar / Google Scholar search on the exact title; any PDF mirror
   e. The authors' homepages or group pages (search "Jiang Choi kirigami unified geometric design framework")
   f. Any GitHub repository, talk slides, or citing papers.
   If a PDF is obtained, save it to papers/related/jiang_choi_2026.pdf and convert with `pdftotext -layout` to papers/related/jiang_choi_2026.txt. If only HTML, save the text. If nothing beyond the abstract is obtainable after all six routes, say so explicitly with the log; then screen the abstract plus whatever secondary material you found.
2. **Read it in full** and write notes/screen_jiang_choi.md with:
   - Bibliographic record: authors with affiliations, date of v1 and latest version, venue if any, abstract verbatim.
   - A one-page summary of what the paper actually does: the objects (which graphs / tilings; whether arbitrary planar graphs; whether cuts are hinge/split as in Segall 2026), the design space it parametrises, whether there is a linear or null-space formulation, whether there is any collision / self-intersection treatment, whether any closed form in the opening angle appears, whether there are certificates or proofs, what experiments and scales.
   - **Table A**: one row per contribution C1–C10 of REPORT.md. Columns: our claim (one line); what Jiang & Choi state on the same object (quote the deciding sentence, with section/page); verdict KNOWN / PARTIAL / NOT FOUND; one-line justification. A KNOWN verdict means their paper already contains the result or an equivalent; PARTIAL means the same primitive or a special case; NOT FOUND means nothing on that object.
   - **Table B**: the constrained-embedding algorithm specifically (C9): convexity constraints at face corners, split-cut outward/inward sign constraint on duplicated vertices, a null-space (Tutte auxetic) parametrisation, a range or 0⁺-margin objective, a sound certificate of the resulting range. One row per ingredient, same columns as Table A.
   - **Verdict paragraph**: is any of C1–C10 KNOWN? If yes, which, and what exact sentence of REPORT.md / IDEA.md must change. Do not soften a KNOWN into a PARTIAL; do not inflate a PARTIAL into a KNOWN.
   - **Citation block**: a correct BibTeX entry for the preprint as retrieved.
3. **Freshness re-screen.** Run these searches restricted to 2026 and record hits with one line each and a verdict: "kirigami deployment range closed form"; "kirigami self-intersection certificate"; "Tutte embedding kirigami null space"; "auxetic tessellation collision-free opening angle"; "hinged kirigami convexity constraint". Anything published after 2026-09-04 that touches C1–C10 goes in a **Table C** with the same columns. Also check whether Segall/Ren/Sorkine-Hornung have a newer paper or preprint (2026) that follows up on Limitation 1 of the TOG paper.

## Rules
- Every verdict cites a quoted sentence with a section number. No verdict from the abstract alone if the full text was obtained.
- Log every query and URL in a final "Query log" section, with what came back.
- Write notes/screen_jiang_choi.md incrementally: skeleton with headings first.
- Do not edit REPORT.md, IDEA.md or STATE.md; the orchestrator propagates.

## Reply
≤15 lines: whether full text was obtained and how; the verdict per C1–C10 as a single line (e.g. "C1 NF, C2 NF, C3 PARTIAL, …"); the Table B verdict; any Table C hit; anything you could not do.
