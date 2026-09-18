# CHANGES_A.md — text changes propagated by Propagator-A (Mission 2 / WP6b)

Date: 2026-09-08. Agent: Propagator-A. Nothing here is committed; the orchestrator commits.

Three decided changes were propagated, each traceable to a verified fact or a directive:

| # | Driver | What changes |
|---|---|---|
| B4 | **F41** (scout-jc; `notes/screen_jiang_choi.md` §6) | Jiang & Choi arXiv:2608.30032 is cleared, not an open risk |
| B5 | **F42** (exp-native200; `results/kill/native200/crashfix.patch`) | The 39 Native200 crashes are an empty-vertex-list access in `convert::to_eig_mat`, not the `get_holes()` segfault |
| B6 | **D15** (human directive) | Physical fabrication is out of scope by decision, not pending work |

**Counts.** 9 sentences/items changed in total: REPORT.md 4, IDEA.md 2, techreport.tex 5
(one of them the bibliography entry's sibling change in refs.bib, listed separately below).
Per document: REPORT.md 4 items, IDEA.md 2 items, techreport.tex 5 items, refs.bib 1 entry.
**No Native200 count was touched** anywhere — the counts stay at 573/600 runs, 172 completed,
362 timed out, 39 crashed, pending WP5's final report — and nothing was added that says the
rerun is complete.

---

## B4 — Jiang & Choi cleared (F41)

### 1. `REPORT.md`, §"What a reviewer will object to", item 15

**Old**

> 15. **"A concurrent preprint may already do this."** arXiv:2608.30032 was retrievable only as
>     an abstract; overlap with the constrained embedding could not be ruled out. Flagged as an
>     open risk, not a clearance.

**New**

> 15. **"A concurrent preprint may already do this."** Cleared. The full text of
>     arXiv:2608.30032 (Jiang and Choi, v1 of 30 Aug 2026, 25 pp.) was retrieved
>     (`papers/related/jiang_choi_2026.pdf`) and screened claim by claim against C1–C10 in
>     `notes/screen_jiang_choi.md`. It contains no collision content, no self-intersection
>     content, no deployment-range quantity, no null space and no certificate, and it works on
>     `N × N` rotating-squares quad patterns rather than arbitrary planar graphs. The single
>     overlap is its Eq. (2), the cross-product corner-convexity inequality, which is the same
>     primitive already conceded to Choi, Dudte and Mahadevan (2019); C9 stays PARTIAL on that
>     primitive alone.

### 2. `REPORT.md`, §Related work, last sentence

**Old**

> Jiang and Choi (arXiv:2608.30032) is a concurrent 2026 preprint whose full text could not be
> retrieved.

**New**

> Jiang and Choi (arXiv:2608.30032, 2026) give a length-based interior-point design framework
> for `N × N` rotating-squares quad kirigami across 2D-to-2D, 2D-to-3D and 3D-to-3D morphing,
> with an inertia-transposition and aspect-ratio theory of the compact end states; they use the
> same cross-product corner inequality (their Eq. (2)) as Choi, Dudte and Mahadevan 2019, have
> no collision predicate and no deployment-range quantity, and check deployment by PyKirigami
> simulation.

### 3. `IDEA.md` §6, the arXiv:2608.30032 bullet

**Old**

> * **arXiv:2608.30032**, "A unified geometric design framework for kirigami structures" — a
>   concurrent 2026 preprint. Only the abstract was retrievable; overlap with the constrained
>   embedding **could not be ruled out** and this is flagged as an open risk, not a clearance
>   (`notes/screen_k9.md`, inaccessible sources). A separate screen found it carries zero
>   collision content (`notes/screen_r1.md`).

**New**

> * **arXiv:2608.30032** (Jiang and Choi), "A unified geometric design framework for kirigami
>   structures" — a concurrent 2026 preprint, **now cleared**. The full text was retrieved and
>   screened claim by claim against C1–C10 (`notes/screen_jiang_choi.md`, Tables A and B): it
>   has zero collision, self-intersection, deployment-range, null-space and certificate
>   content, and works on `N × N` rotating-squares quad patterns rather than arbitrary planar
>   graphs. The only overlap is its Eq. (2), the cross-product corner-convexity inequality, so
>   **C9 remains PARTIAL on that primitive alone**, already attributed to
>   Choi–Dudte–Mahadevan 2019. This supersedes the abstract-only screen in
>   `notes/screen_k9.md`.

### 4. `paper/techreport/techreport.tex`, objections list (was line 3724)

**Old**

> \item \emph{``A concurrent preprint may already do this.''} arXiv:2608.30032 (Jiang and Choi)
> was retrievable only as an abstract; overlap with the constrained embedding could not be ruled
> out. Flagged as an open risk, not a clearance.

**New** — cleared, with the C1–C10 screen cited, the Eq. (2) overlap named, C9 kept
\textsc{partial}, and `\cite{Jiang2026}` / `\cite{Choi2019}` wired in.

### 5. `paper/techreport/refs.bib`, entry `Jiang2026`

A placeholder entry existed (`author = {Jiang and Choi}`, no eprint fields). Replaced with the
verified record from `notes/screen_jiang_choi.md` §2:

```bibtex
@misc{Jiang2026,
  author        = {Jiang, Qinghai and Choi, Gary P. T.},
  title         = {A unified geometric design framework for kirigami structures},
  year          = {2026},
  month         = aug,
  eprint        = {2608.30032},
  archivePrefix = {arXiv},
  primaryClass  = {cond-mat.soft},
  doi           = {10.48550/arXiv.2608.30032},
  note          = {arXiv:2608.30032v1, submitted 30 Aug 2026, 25 pp.}
}
```

---

## B5 — Native200 crash attribution corrected (F42)

The reproduced cause is an empty-vertex-list access in `convert::to_eig_mat`, reached from
`Hmesh::merge_close_verts` after a deployed configuration has collapsed to zero faces. The
`UnitPattern::get_holes()` null dereference is real but **latent**, and the same patch
(`results/kill/native200/crashfix.patch`) fixes both. No count was changed.

### 6. `REPORT.md`, §"Errata against the published pipeline", item 4

**Old**

> 4. **Upstream defects observed, not claimed as contributions:** `UnitPattern::get_holes()`
>    segfaults on bounded patches; their collision-aware `θ_max` runs `merge_close_verts` at
>    0.1·average edge before testing, fusing hinge duplicates, and reports 0.067 on snub square
>    where first contact is ≈1.65 (`STATE.md` F24).

**New**

> 4. **Upstream defects observed, not claimed as contributions:** their collision-aware `θ_max`
>    runs `merge_close_verts` at 0.1·average edge before testing, fusing hinge duplicates, and
>    reports 0.067 on snub square where first contact is ≈1.65 (`STATE.md` F24). The same
>    `merge_close_verts` also crashes: once a deployed configuration has collapsed to zero
>    faces it leaves an empty vertex list, and the collision scan's next call reaches
>    `convert::to_eig_mat` on it and dereferences element 0 of an empty vector —
>    `EXC_BAD_ACCESS`. This is the reproduced cause of the Native200 crashes (`STATE.md` F42;
>    `results/kill/native200/crashfix.patch`). A second, latent defect sits in
>    `UnitPattern::get_holes()`, which dereferences `prev()->twin()` without a null check on a
>    bounded patch; the same patch guards both.

### 7. `IDEA.md` §7, mandatory-caveats Native200 paragraph

**Old** (clause only)

> 172 completed, 362 timed out at 600 s, 39 crashed on the known upstream bugs (F24).

**New**

> 172 completed, 362 timed out at 600 s, 39 crashed. The crash cause was reproduced (F42): an
> empty-vertex-list access in `convert::to_eig_mat`, reached from `Hmesh::merge_close_verts`
> after a deployed configuration collapses to zero faces — not the `UnitPattern::get_holes()`
> segfault the earlier F24 note assumed. `results/kill/native200/crashfix.patch` guards it and
> also adds the missing `prev()->twin()` null check in `get_holes()`, a second latent defect.

**Not changed:** the `IDEA.md` §4 table row (line 140) prints "39 crashed" with **no** cause
attribution, so it needed no edit and none was made. Its count is untouched.

### 8. `paper/techreport/techreport.tex` §2, "Two upstream defects we observed"

**Old** (first clause)

> \file{UnitPattern::get\_holes()} segfaults on bounded patches through an unchecked
> \file{prev()->twin()}; and their collision-aware $\theta_{\max}$ runs …

**New** — the `merge_close_verts` fusing defect is stated first and kept as the false-negative
mechanism; the crash is then attributed to the empty-vertex-list access in
`convert::to_eig_mat` reached from `merge_close_verts`, citing F42 and the patch; and
`get_holes()` is kept as the **second, latent** defect fixed by the same patch. (The
subsubsection title "Two upstream defects we observed" still holds: there are two.)

### 9. `paper/techreport/techreport.tex` §9, Native200 caveat paragraph (was line 3401)

**Old**

> \file{quad\_random} is worst hit: 149/188 timed out and 26/188 crashed, consistent with the
> unchecked \file{prev()->twin()} in \file{get\_holes()} on bounded patches.

**New**

> \file{quad\_random} is worst hit: 149/188 timed out and 26/188 crashed. The crash cause was
> reproduced afterwards (\file{STATE.md} F42): an empty-vertex-list access in
> \file{convert::to\_eig\_mat}, reached from \file{Hmesh::merge\_close\_verts} once a deployed
> configuration has collapsed to zero faces, not the unchecked \file{prev()->twin()} in
> \file{get\_holes()} that the original F24 note assumed --- that null dereference is real but
> latent, and the same patch (\file{results/kill/native200/crashfix.patch}) fixes both.

Counts in this paragraph (573/600, 172/362/39, 149/188, 26/188, median $|F|$) are untouched,
as is the "rerunning the 362 timeouts with a longer budget is the one remaining item" sentence.

---

## B6 — Fabrication out of scope by decision (D15)

The `export/` artefacts and every sentence describing the existing SVG/3MF/STL output are left
exactly as they are; they exist and were verified with `xmllint --noout` and `unzip -t`. Only
forward-looking fabrication promises were removed.

### 10. `REPORT.md`, §"Next steps toward a paper", item 6

**Old**

> 6. **Fabricate one hero physically** and photograph it closed and open. The highest-value
>    remaining item for a graphics venue, and the only one needing hardware.

**New**

> 6. **Physical fabrication is out of scope by decision**, not pending work. The project is
>    theory and computation only; nothing is to be laser-cut, printed or photographed, and no
>    claim in this report depends on a physical artifact. The `export/` SVG/3MF/STL artefacts
>    stay as verified file-format output (`xmllint`, `unzip -t`) and are not a fabrication
>    result.

### 11. `REPORT.md`, objection 16

**Old**

> 16. **"No physical artifact."** None was fabricated; no slicer is installed and no 3MF was
>     opened in one. The constructive review's strongest suggestion — laser-cut and actuate the
>     `ε_max ≈ 2 rad` outlier and photograph it — was not done.

**New**

> 16. **"No physical artifact."** Correct, and it is a scope decision rather than an omission
>     to be repaired: the project is theory and computation only, so nothing was fabricated, no
>     slicer is installed and no 3MF was opened in one. The constructive review's suggestion —
>     laser-cut and actuate the `ε_max ≈ 2 rad` outlier and photograph it — falls outside that
>     scope. Every claim here is stated about the geometric model, and none is offered as
>     evidence about a physical piece.

### 12. `paper/techreport/techreport.tex` §11, open-items list (was line 3775)

**Old**

> \item \textbf{Fabricate one design physically} and photograph it closed and open. The highest
> value remaining item, and the only one needing hardware.

**New**

> \item \textbf{Physical fabrication is out of scope by decision}, not an open item. This work
> is theory and computation only: nothing is to be cut, printed or photographed, and no claim
> made here rests on a physical artifact. The \file{export/} SVG/3MF/STL output stays as
> verified file-format generation, not as a fabrication result.

### 13. `paper/techreport/techreport.tex`, "Review 3" demands paragraph

**Old**

> Two have not: no design has been physically fabricated, and there is still \textbf{no
> characterisation of when the constrained set is non-empty}.

**New**

> One has not: there is still \textbf{no characterisation of when the constrained set is
> non-empty}. The remaining demand, physical fabrication, is out of scope by decision rather
> than outstanding --- this work is theory and computation only.

**Not changed:** the §"Fabrication export" subsection and its "Not verified: no slicer is
installed on this machine, so no 3MF has been opened in one, and nothing has been physically
fabricated" paragraph. That is a factual statement about artefacts that exist and was left
standing, per the spec.

---

## Rebuild check

`paper/techreport/techreport.tex` was rebuilt (`pdflatex`, `bibtex`, `pdflatex` ×2): **63 pages,
0 errors, 0 undefined citations**, with the new `Jiang2026` and existing `Choi2019` keys
resolving. `techreport.pdf`, `.aux`, `.bbl`, `.blg`, `.log`, `.toc` are regenerated artefacts
of that build.
