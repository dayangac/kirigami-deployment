# Scout-a: kirigami / auxetic geometry literature screen

Output: notes/field_kirigami.md

Goal: a related-work table for the geometry of rigid hinged kirigami / rotating-rigid-unit auxetics, so a later ideation phase knows exactly what already exists. The base papers are Segall et al. 2025 (Reconfigurable Hinged Kirigami Tessellations, SA'25) and Segall et al. 2026 (Uniformly Deployable Kirigami on Arbitrary Planar Graphs, TOG 45(4)). Read notes/paper_2026.md and notes/paper_2025.md first if present, else papers/paper_2026_tutte_flow.txt Sec 1–2.

Use WebSearch / WebFetch. Cover AT LEAST these works and anything they lead to:
- Konaković et al. 2016 (Beyond Developable), 2018 (Rapid Deployment via Programmable Auxetics)
- Choi, Dudte, Mahadevan 2019 (Programming shape using kirigami tessellations), 2021 (Compact reconfigurable kirigami); Dudte et al. 2023 (Additive framework for kirigami design); Chen, Choi, Mahadevan 2020 (Deterministic and stochastic control of kirigami topology)
- Jiang, Rist, Pottmann, Wallner 2020 (Freeform quad-based kirigami), 2022 (Shape-morphing mechanical metamaterials), Jiang et al. 2024 (Quad mesh mechanisms)
- Dang, Feng, Duan, Wang 2021 (PRE theorem deployable kirigami tessellations different topologies), 2022 (PRL spherical kirigami compatibility)
- Warisaya, Sato, Tachi 2022 (Freeform auxetic mechanisms based on corner-connected tiles)
- Zaman et al. 2025 (One string to pull them all, TOG)
- Liu, Lu, Cao, Deussen, Tu 2024 (Auxetic dihedral Escher tessellations)
- Mitschke, Robins, Mecke, Schröder-Turk 2013 (Finite auxetic deformations of plane tessellations) — geometry side only; Scout-b covers rigidity
- Grima & Evans rotating rigid units (rotating squares/triangles/rectangles), Rafsanjani & Pasini 2016, Shan et al. 2015 (planar isotropic negative Poisson ratio), Tang & Yin 2017
- Any 2024–2026 work on: collision-free deployment range of kirigami, non-uniform hinge angles / multi-DOF kirigami mechanisms, design-space / null-space characterizations of auxetic tilings, joint topology+geometry optimization of cut patterns, bistable/multistable kirigami as an inverse problem, spatially varying deployment.

For EACH row of the table record: citation (authors, year, venue), pattern class handled (triangles / quads / 2-colorable / arbitrary), what it ASSUMES, what it ACHIEVES (theorem / algorithm / fabrication), what it CANNOT do (quote or precise paraphrase, with a page/section if you can see it), and relation to the Segall 2026 framework (subsumed by / orthogonal / extends). Aim for ≥ 15 rows.

Then a section "Gaps": for each of the seven idea seeds below, list the closest existing work and say in one sentence whether the seed appears open:
1. closed-form collision / deployment range per split-cut pair; range-maximizing embeddings
2. full configuration space via rigidity Jacobian; designed DOF
3. existence / injectivity theorem for Tutte auxetic embeddings; combinatorial rank formula for the null space
4. joint orientation–geometry optimization replacing max-cut-then-solve
5. spatially varying deployment fields (non-uniform θ as a designed field)
6. multi-target reconfiguration (one flat sheet, several deployed configurations)
7. stability / bistability from hinge stiffness as an inverse problem on the null space

Finally list EVERY search query you ran, verbatim, and which returned nothing useful. Honesty over coverage.
