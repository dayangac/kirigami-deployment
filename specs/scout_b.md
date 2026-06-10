# Scout-b: rigidity theory, linkages, mechanism configuration spaces

Output: notes/field_rigidity.md

Goal: related-work table on the mathematics that could support a theorem about the configuration space of rigid-face hinged kirigami (faces = rigid bodies, hinges = pin joints in 2D). Base framework: Segall et al. 2026 (Tutte auxetic embeddings; uniform opening angle θ is one path in configuration space; paper's Sec 6 says analyzing the DOF of the deployment space "could be achieved by sensitivity analysis" and leaves it open). Read notes/paper_2026.md if present.

Use WebSearch / WebFetch. Cover AT LEAST:
- Bar-and-joint and body-and-hinge (body-bar) rigidity: Maxwell counting, Laman's theorem, Tay's theorem for body-bar frameworks, Whiteley; generic vs. special-position rigidity; infinitesimal vs finite flexes; pebble game algorithms (Jacobs & Hendrickson).
- Mitschke, Robins, Mecke, Schröder-Turk 2013 (Finite auxetic deformations of plane tessellations) and related Schröder-Turk works on skeletal / tessellation mechanisms and Poisson ratio.
- Guest & Hutchinson 2003 (On the determinacy of repetitive structures), Guest & Fowler on symmetry-extended mobility, Kangwai/Guest.
- Connelly, Connelly–Whiteley (second-order rigidity, prestress stability), Connelly & Servatius; tensegrity stability framework.
- Configuration spaces of planar linkages (Kapovich–Millson, Thurston), rotating-rigid-unit (Grima) mechanism DOF counts, Wunderlich / Kokotsakis / Stachel flexible polyhedra if relevant.
- Periodic frameworks: Borcea & Streinu (periodic rigidity, deformations of crystal frameworks), Ross, Schulze & Whiteley periodic/symmetric rigidity.
- Any computational work computing the full nonlinear configuration space of a rigid-body mechanism from the Jacobian (numerical continuation, bifurcation detection).
- Origami rigidity analogues (Tachi rigid-foldability, Demaine, Schenk & Guest 2013) only insofar as the DOF-counting machinery transfers.

For EACH row record: citation, object studied (bar-joint / body-hinge / periodic / linkage), main result (theorem or algorithm), assumptions (generic position? planar? finite vs infinitesimal?), what it CANNOT do or does not cover (quoted or precisely paraphrased), and applicability to hinged kirigami (does the counting apply directly to faces-as-bodies with shared pin joints at hinge vertices?). Aim for ≥ 12 rows.

Then a section "What is known vs open" about: (a) DOF of a 2D body-hinge framework whose graph is a planar 'hinge graph' (nodes = faces, edges = hinge vertices); (b) whether uniform-θ deployability implies a 1-DOF mechanism or can coexist with extra finite flexes; (c) tools to compute finite (not just infinitesimal) mobility. Say explicitly what you could not verify.

Finally list EVERY search query you ran, verbatim, and which returned nothing useful.
