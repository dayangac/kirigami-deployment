# Experimenter: Phase-5 kill experiments for the top-ranked ideas

Output: results/kill/ (one subdirectory per experiment with CSV/JSON dumps, PNG figures via `arch -arm64 /usr/local/bin/python3` + matplotlib, and a KILL_REPORT.md at results/kill/KILL_REPORT.md), plus any small C++ apps under code/apps/ and helper code under code/src/method/ (namespace kiri::method) needed to run them. All computation in C++ (D3).

Read first: specs/common_preamble.md; code/README.md (API); STATE.md (F1–F21, U1–U9, D5); ideas/ranking.md Sec 4 (exact kill specs K1a, K1b, K1c, K2a, K2b, K2c, K3a, K3b) and Sec 0 (P1–P3); results/core_validation/{rank_claim.md,reference_cases.md}.

Run, in this order, exactly as specified in ideas/ranking.md Sec 4 (if a spec detail is missing, choose the simplest faithful option and record it):
1. K3a (2-core mobility identity) — also re-measure mobility at θ = 0.3·θ_max, not only at θ = 0, and report both.
2. K1b (harmonic identity: every (vertex, edge) orientation determinant and every face signed area along deployment is p + q cos θ + r sin θ; face areas constant) — this gates everything else.
3. K1a (validity/usability of the shape space: inverted faces of X_ini, of Eq.(6) X0, and fraction p_valid of samples of X = X0 + Φt that are valid; both Gaussian sampling and a local trust-region probe from the projection of X_ini).
4. K2c (locality gate for the active set).
5. K2b (the margin: closed-form-root softmin range maximization with a signed-area log barrier, over null coefficients, vs the paper's Eq.(6)+Eq.(9) pipeline with its γ ladder; both refereed by the same independent collision bisection).
6. K1c (Eq.(9) false-negative rate) and K2a (closed-form θ_max vs bisection) if time permits (< 30 min each).
For K2b implement the closed-form first-contact roots properly: vertex-into-edge-interior contact = harmonic root PLUS the two interval (projection) inequalities, each also a harmonic; enumerate candidate (vertex, edge) pairs with the swept-disc broad phase from K2c; validate every closed-form θ_max against collision.hpp bisection (tolerance 1e-5) and report the agreement rate.

KILL_REPORT.md must contain, per experiment: the exact PASS/FAIL rule copied from ranking.md, the measured numbers, the verdict, wall time, and the paths of the artifacts. If K1b fails, stop and report. If K1a returns the emptiness branch, say so plainly; do not tune the sampler to get a different answer. Also fold in the three orchestrator rank checks already in rank_claim.md: restate their final numbers (H == |E_hinge| − |F| + c(Γ); L = R·D; periodic rank(L) = H − 1) and resolve the prose/table inconsistency the Critic flagged in rank_claim.md by fixing the prose to match the tables.
All existing tests must still pass after your additions; add doctest cases for any new numerical routine (closed-form roots vs bisection, harmonic fit). Do not commit. Reply with the verdict table (experiment, PASS/FAIL, key number) and the two or three most surprising numbers.
