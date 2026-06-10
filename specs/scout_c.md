# Scout-c: Tutte embeddings, discrete geometry, convex/spectral embedding theory

Output: notes/field_tutte.md

Goal: related-work table on the mathematics of Tutte-type linear embeddings, because Segall et al. 2026 define a "Tutte auxetic embedding" as the solution of a Laplacian-like linear system [L; B] X = [0; T] where rows of L are vector-sum constraints over hinge-cut edges in each hole preimage (one row per hole), and they state WITHOUT PROOF that (i) the system is full rank iff there are no split-cut edges (one hole per interior vertex), and (ii) otherwise the null space dimension equals #interior vertices − #holes. They also do not address injectivity (whether the solution is a valid straight-line planar embedding). Read notes/paper_2026.md Sec 4.4 if present, else papers/paper_2026_tutte_flow.txt around "Tutte auxetic embedding".

Use WebSearch / WebFetch. Cover AT LEAST:
- Tutte 1963 "How to draw a graph" (spring embedding theorem: 3-connected planar + convex boundary ⇒ planar straight-line embedding); Floater 2003 (One-to-one piecewise linear mappings over triangulations; convex combination maps); Gortler, Gotsman, Thurston 2006 (Discrete one-forms on meshes and applications to 3D mesh parameterization) — the one-form / index argument for injectivity; Colin de Verdière; Linial–Lovász–Wigderson rubber-band embeddings.
- Generalizations: non-symmetric weights, negative weights (Gortler et al.), Tutte embeddings with non-convex boundary, Tutte on non-triangulations / polygonal meshes, Tutte-type embeddings with additional linear constraints, weighted Laplacians whose rows are not per-vertex (e.g. constraints per face or per cycle).
- Rank of incidence/Laplacian-type matrices whose rows are sums over edge sets: connections to cycle/cut spaces, matroid rank, planar duality (rows indexed by holes ~ faces of a derived graph?). Look for a known counting formula.
- Discrete differential geometry of cut-and-rotate structures: Grima rotating rigid units as a linear map; "auxetic" Laplacians; harmonic maps on non-manifold meshes.
- Spectral / convex embedding theory that might give injectivity for solutions of constraint systems (e.g. maximum principle arguments, discrete harmonic functions with mixed constraints).
- Any work on "parallel-redrawing" / "parallel morphs" of planar graphs (Whiteley's parallel redrawing, Ross/Streinu): drawings with prescribed edge directions — highly relevant since Eq.(2) is a vector-sum (closure) constraint and Remark A.4 gives parallel duplicated edges.

For EACH row: citation, statement of the main theorem/algorithm (precise), hypotheses, what it does NOT cover (quoted/paraphrased), and relevance to (a) proving the rank formula, (b) proving injectivity of Tutte auxetic embeddings, (c) parametrizing the null space combinatorially. Aim for ≥ 12 rows.

Then a section "Candidate proof routes" (≤ 1 page): sketch, with references, how one might prove (i) rank(L) = #holes for generic/any orientation, (ii) a formula dim(null) = f(#interior vertices, #holes, #split edges, #components of split forest), (iii) injectivity under convex fixed boundary. Mark each route as 'known technique', 'plausible', or 'speculative'.

Finally list EVERY search query you ran, verbatim, and which returned nothing useful.
