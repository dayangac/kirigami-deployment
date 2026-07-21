# method/expansive_cone.jl -- the expansive cone of a cut structure and the linear program
# that decides it (ideas/round2_adversary.md A1, ideas/ranking_r2.md X1 / kill_k8a).
#
# Port of code/src/method/expansive_cone.{hpp,cpp}. 1-based indices throughout; a flex
# vector stacks (omega_f, w_fx, w_fy) at positions 3(f-1)+1..3(f-1)+3. C++ out-pointers
# become extra return values (`farkas_residual` returns `(resid, lambda_min, sum_err)`,
# `normalise_rows!` returns `n_degenerate`). The LP solver is the C++'s own (smoothing
# loop + away-step Frank-Wolfe on the min-norm point); no LP library is used, so it is
# ported literally.
#
# WHAT IS NEW HERE, AND WHAT IS NOT. The object "polyhedral cone of first-order motions
# cut out by homogeneous linear inequalities on the velocities, decided by one linear
# program" is classical: G. Rote, F. Santos, I. Streinu, "Expansive Motions and the
# Polytope of Pointed Pseudo-Triangulations", Discrete and Computational Geometry -- The
# Goodman-Pollack Festschrift, Algorithms and Combinatorics 25, Springer 2003, Lemma 3.2,
# and R. Connelly, E. D. Demaine, G. Rote, Discrete Comput. Geom. 30 (2003) 205-239.
# See notes/screen_r2.md, section A1: the verdict there is PARTIAL. Their variables are
# point velocities with one inequality per PAIR of points; here the variables are face
# velocities of a body-and-pin framework with a CUT structure and the inequalities are
# one per split edge plus a branch per corner incidence. Nothing in this file may be
# described as a new kind of object; only its instantiation on a cut structure is.
#
# THE CONSTRUCTION. A flex of the flat structure assigns every face f an angular
# velocity omega_f and a translation w_f, so the velocity of the material point of f
# sitting at y is
#
#     V_f(y) = omega_f * J y + w_f ,        J = rot(+pi/2) .
#
# The hinge edges impose V_f(p_i) = V_g(p_i) at the pin p_i (the hinge SOURCE vertex),
# which is exactly `mobility.jl::build_rigidity`. Its kernel is the space of flexes,
# and it factors through `build_A`: given omega with A omega = 0 the translations are
# determined along the hinge forest up to one free translation per component of Gamma,
# so
#     dim ker R  =  dim ker A  +  2 c(Gamma) .
# `flex_basis` builds an orthonormal basis of ker R this way and checks the residual.
#
# THE INEQUALITIES. At theta = 0 every copy of a source vertex v sits at X_v, so the
# relative velocity of two copies is a linear function of the flex, and X being FIXED
# the edge vectors are constants. With the sign convention of `zero_plus.jl`:
#
#   split edge e = {v_from, v_to} between f0 (which stores the direction) and f1,
#   d_e = X[v_to] - X[v_from], dV = V_{f1}(X[v_from]) - V_{f0}(X[v_from]) :
#
#       row_e(flex) = det( dV, d_e )   > 0     <=>   the cut opens.
#
#   Since zero_plus.jl's dS_e = 2 dV, row_e = q_e / 2 exactly; the doctests assert it.
#
#   corner incidence (copy p of v, in face f_p, against the corner of face f_a at v),
#   e1 = X[v_next] - X[v], e2 = X[v_prev] - X[v], dV = V_{f_p}(X_v) - V_{f_a}(X_v) :
#
#       g1 = det(e1, dV),  g2 = det(dV, e2),   cross = det(e1, e2),
#       mu = max(-g1, -g2)  at a convex corner (cross > 0),
#       mu = min(-g1, -g2)  at a reflex corner.
#
#   So the two rows carried per incidence are  row_G1 = -g1 = det(dV, e1)  and
#   row_G2 = -g2 = -det(dV, e2), and mu = mu_zero_plus / 2 exactly.
#
# THE BRANCHES. At a reflex corner the margin is a MIN, so both rows are required and
# the constraint is convex. At a convex corner it is a MAX, so the feasible set is a
# UNION of two half-spaces and P(X) is a union of polyhedra indexed by one branch per
# convex incidence. This file does not enumerate them. It solves, in order:
#
#   pass 0  ("and"):   require BOTH rows at every incidence, convex ones included.
#                      This is a SUBSET of P(X), so feasibility here is already a sound
#                      existence proof and needs no branch selection at all.
#   pass 1  (sigma):   at each convex incidence keep only the row the UNIFORM ray itself
#                      satisfies better -- "the side the copy is already on", read off the
#                      sigma velocities, the only geometrically distinguished branch at a
#                      flat state where all copies of a vertex coincide.
#   pass k > 1:        at each convex incidence keep only the row that the best witness so
#                      far satisfies better (argmax of the two), then re-solve.
#
# Feasibility of any pass is a proof that P(X) is non-empty. Infeasibility of every pass
# proves nothing, because the branch enumeration is not exhaustive. That asymmetry is
# stated here so that no caller can mistake a FAIL for an emptiness theorem.
#
# THE LINEAR PROGRAM. With rows normalised to unit Euclidean norm in the basis
# coordinates z, the quantity solved for is the matrix-game value
#
#     val = max_{||z||_2 <= 1} min_i a_i . z  =  min_{lambda in simplex} || A^T lambda ||_2
#
# (Sion; both sides are >= 0 because z = 0 is feasible, and val > 0 iff the strict system
# is feasible). Any primal z gives a rigorous LOWER bound min_i a_i.z / ||z||_2 and any
# lambda in the simplex a rigorous UPPER bound ||A^T lambda||_2, so the solver returns a
# two-sided bracket and never has to be trusted. A lambda with a small upper bound is an
# approximate Farkas certificate: sum_i lambda_i a_i ~ 0 with lambda >= 0, sum = 1.
#
# THE SOLVER (rewritten; see results/kill/k8a/recheck.md). The dual problem is the
# MINIMUM-NORM POINT of the convex hull of the rows,  w* = argmin { ||w|| : w in
# conv{a_i} },  and the primal witness is NOT a separate object: at the minimum-norm
# point the supporting-hyperplane inequality  a_i . w* >= ||w*||^2  holds for every i
# (else moving from w* towards a_i decreases the norm), so
#
#     z* = w* / ||w*||   satisfies   min_i a_i . z*  >=  ||w*||  =  dual bound,
#
# and with margin <= val <= dual the bracket CLOSES: margin = val = ||w*||. The original
# solver computed the dual by Frank-Wolfe but took its primal witness from an unrelated
# log-sum-exp smoothing loop, which is very slowly convergent and returned margin = 0 on
# systems whose feasibility is known independently (K9's 23 designs with sigma in P(X)).
# This version runs AWAY-STEP Frank-Wolfe (Lacoste-Julien and Jaggi, NeurIPS 2015 --
# linearly convergent on a polytope, unlike vanilla FW's O(1/t)) on the min-norm point
# and reads the primal witness off the same iterate w, so both ends of the bracket
# improve together. The smoothing loop is kept only as an extra source of witnesses; the
# reported margin is the best of the two. No external LP library is used.
#
# DEVIATION FROM THE SPEC. ideas/ranking_r2.md's kill_k8a asks for the normalisation
# ||z||_inf <= 1 and multiplicative weights on the l1 dual. The cone is homogeneous, so
# the two normalisations are positive iff each other; the l2 game is better conditioned
# and its dual bound comes free from the same softmax weights. `margin_inf` reports the
# spec's quantity for the returned witness.

# ---------------------------------------------------------------------------
# Flexes of the body-and-pin framework.

"""Velocity of the material point of face f at y under the flex `u` (length 3F)."""
function flex_velocity(u::AbstractVector{Float64}, f::Int, y::Vec2)
    om = u[3 * (f - 1) + 1]
    return Vec2(-om * y[2] + u[3 * (f - 1) + 2], om * y[1] + u[3 * (f - 1) + 3])
end

# Coefficients of  flex |-> det( V_f(y), cvec )  on the three unknowns of face f.
# det(V_f(y), c) = omega_f * det(J y, c) + det(w_f, c)
#               = -omega_f * (y . c) + w_fx * c.y - w_fy * c.x .
function _det_row!(I::Vector{Int}, J::Vector{Int}, V::Vector{Float64}, f::Int, y::Vec2,
                   cvec::Vec2, s::Float64, row::Int)
    push!(I, row); push!(J, 3 * (f - 1) + 1); push!(V, -s * dot(y, cvec))
    push!(I, row); push!(J, 3 * (f - 1) + 2); push!(V, s * cvec[2])
    push!(I, row); push!(J, 3 * (f - 1) + 3); push!(V, -s * cvec[1])
    return nothing
end

# M'-vertex -> some face containing it (0 if none). Two faces sharing an M'-vertex are
# hinged at that vertex, so their velocities there agree and the choice is immaterial.
function _prime_face_map(c::CutStructure)
    pv_face = zeros(Int, c.n_prime_vertices)
    for f in eachindex(c.prime_faces), pv in c.prime_faces[f]
        (pv >= 1 && pv_face[pv] == 0) && (pv_face[pv] = f)
    end
    return pv_face
end

mutable struct FlexBasis
    N::Matrix{Float64}       # 3F x n, orthonormal columns
    F::Int
    dim_ker_A::Int           # dim of the angular-velocity kernel (mobility.jl::build_A)
    components::Int          # c(Gamma)
    residual::Float64        # max |R N| / max|R|, the self-check that N really is a flex
end
FlexBasis() = FlexBasis(zeros(0, 0), 0, 0, 0, 0.0)
dim(fb::FlexBasis) = size(fb.N, 2)

"""
    flex_basis(g, pins, rel_tol = 1e-10) -> FlexBasis

Orthonormal basis of ker(build_rigidity(g, pins)). `rel_tol` is the relative singular
value threshold used for ker A.
"""
function flex_basis(g::HingeGraph, pins::Vector{Vec2}, rel_tol::Float64 = 1e-10)
    fb = FlexBasis()
    fb.F = g.F
    fb.components = g.components
    g.F == 0 && return fb

    # ker A : the admissible angular velocities.
    local Omega::Matrix{Float64}
    if n_cycles(g) == 0
        Omega = Matrix{Float64}(I, g.F, g.F)
    else
        Ad = Matrix(build_A(g, pins))
        nrm = maximum(abs, Ad)
        nrm > 0 && (Ad ./= nrm)
        S = svd(Ad; full = true)
        sv = S.S
        smax = isempty(sv) ? 0.0 : sv[1]
        thr = rel_tol * max(1.0, smax)
        r = count(s -> s > thr, sv)
        Omega = S.V[:, r+1:g.F]
    end
    fb.dim_ker_A = size(Omega, 2)

    # Integrate the translations along the BFS forest, and add the per-component free
    # translations. Order: roots first so that a child's parent is already done
    # (std::stable_sort by depth; Julia's default sort is stable).
    order = sort(collect(1:g.F); by = f -> g.depth[f])
    # Component label per face (root of its BFS tree).
    comp = zeros(Int, g.F)
    nc = 0
    for f in order
        comp[f] = g.parent[f] == 0 ? (nc += 1) : comp[g.parent[f]]
    end

    n = fb.dim_ker_A + 2 * nc
    B = zeros(3 * g.F, n)
    for j in 1:fb.dim_ker_A
        for f in 1:g.F
            B[3 * (f - 1) + 1, j] = Omega[f, j]
        end
        for f in order
            g.parent[f] == 0 && continue  # w = 0 at the root
            p = g.parent[f]
            i = g.pedge[f]
            Jp = Vec2(-pins[i][2], pins[i][1])
            dom = Omega[p, j] - Omega[f, j]
            B[3 * (f - 1) + 2, j] = B[3 * (p - 1) + 2, j] + dom * Jp[1]
            B[3 * (f - 1) + 3, j] = B[3 * (p - 1) + 3, j] + dom * Jp[2]
        end
    end
    for f in 1:g.F
        B[3 * (f - 1) + 2, fb.dim_ker_A + 2 * (comp[f] - 1) + 1] = 1.0
        B[3 * (f - 1) + 3, fb.dim_ker_A + 2 * (comp[f] - 1) + 2] = 1.0
    end

    # Householder QR, thin Q (the C++ householderQ() * Identity(3F, n)).
    fb.N = Matrix(qr(B).Q)[:, 1:n]

    # Self-check: N really is a flex.
    if n_edges(g) > 0
        R = build_rigidity(g, pins)
        rn = maximum(abs, R)
        rn > 0 && (fb.residual = maximum(abs, R * fb.N) / rn)
    end
    return fb
end

"""
    sigma_flex(c, X) -> Vector{Float64}

The uniform deployment ray sigma, as a flex vector of length 3F, read off the exact
kinematics: omega_f = -sigma_f / 2 and w_f = (S_{(v,f)} + sigma_f J X_v) / 2 for any
vertex v of f, where S = 2 dY/dtheta|_0. Whether it IS a flex of the framework at this
X is not assumed: `expansive_cone` reports ||R u|| / (||R|| ||u||), which is ~0 only when
X satisfies the deployability condition Eq. (2).
"""
function sigma_flex(c::CutStructure, X::Vector{Vec2})
    m = c.mesh
    F = n_faces(m)
    d = deploy(c, X, 0.0)
    u = zeros(3 * F)
    for f in 1:F
        isempty(m.faces[f]) && continue
        sg = Float64(m.sigma[f])
        v = m.faces[f][1]
        pv = prime_vertex(c, f, v)
        pv < 1 && continue
        Jx = Vec2(-X[v][2], X[v][1])
        w = d.dY_dtheta[pv] + 0.5 * sg * Jx
        u[3 * (f - 1) + 1] = -0.5 * sg
        u[3 * (f - 1) + 2] = w[1]
        u[3 * (f - 1) + 3] = w[2]
    end
    return u
end

# ---------------------------------------------------------------------------
# The linearised 0+ non-collision system.

# The C++ enum class ConeRowKind {Split, CornerG1, CornerG2}; the members carry a `Cone`
# prefix because Julia enum members are module-level names and `Split` is an EdgeType.
@enum ConeRowKind::UInt8 ConeSplit = 0 ConeCornerG1 = 1 ConeCornerG2 = 2

struct ConeRow
    kind::ConeRowKind
    index::Int      # into split_copies() / corner_incidences()
    convex::Bool    # corner rows only: the incidence is at a convex corner
end

mutable struct ConeSystem
    A::SparseMatrixCSC{Float64,Int}   # n_rows x 3F, unnormalised
    rows::Vector{ConeRow}
    n_split::Int                      # number of split edges / corner incidences
    n_corner::Int
    n_convex::Int                     # convex incidences (the disjunctive ones)
end
ConeSystem() = ConeSystem(spzeros(0, 0), ConeRow[], 0, 0, 0)

"""
    cone_system(c, X) -> ConeSystem

All rows, at the flat embedding X. Row values are exactly q_e/2 and -g1, -g2 of
zero_plus.jl when evaluated at the sigma flex.
"""
function cone_system(c::CutStructure, X::Vector{Vec2})
    cs = ConeSystem()
    m = c.mesh
    F = n_faces(m)
    sp = split_copies(c)
    ci = corner_incidences(c)
    pv_face = _prime_face_map(c)
    cs.n_split = length(sp)
    cs.n_corner = length(ci)

    I = Int[]; J = Int[]; V = Float64[]
    row = 0
    for i in 1:cs.n_split
        s = sp[i]
        y = X[s.v_from]
        de = X[s.v_to] - X[s.v_from]
        row += 1
        _det_row!(I, J, V, s.f1, y, de, +1.0, row)
        _det_row!(I, J, V, s.f0, y, de, -1.0, row)
        push!(cs.rows, ConeRow(ConeSplit, i, false))
    end
    for i in 1:cs.n_corner
        z = ci[i]
        fp = pv_face[z.pv_other]
        (fp < 1 || z.face < 1) && continue
        y = X[z.v]
        e1 = X[z.v_next] - X[z.v]
        e2 = X[z.v_prev] - X[z.v]
        convex = _det2(e1, e2) > 0
        convex && (cs.n_convex += 1)
        # row_G1 = -g1 = det(dV, e1)
        row += 1
        _det_row!(I, J, V, fp, y, e1, +1.0, row)
        _det_row!(I, J, V, z.face, y, e1, -1.0, row)
        push!(cs.rows, ConeRow(ConeCornerG1, i, convex))
        # row_G2 = -g2 = -det(dV, e2)
        row += 1
        _det_row!(I, J, V, fp, y, e2, -1.0, row)
        _det_row!(I, J, V, z.face, y, e2, +1.0, row)
        push!(cs.rows, ConeRow(ConeCornerG2, i, convex))
    end
    cs.A = sparse(I, J, V, row, 3 * F)   # duplicates summed, as setFromTriplets
    return cs
end

# ---------------------------------------------------------------------------
# The LP.

Base.@kwdef mutable struct ConeLPOptions
    mu0::Float64 = 1.0          # initial smoothing, in units of the (unit) row norms
    mu_min::Float64 = 1e-5
    mu_factor::Float64 = 0.5
    iters_per_stage::Int = 60
    tol::Float64 = 1e-8         # margin above which the witness is declared feasible
    dual_tol::Float64 = 1e-7    # dual bound below which we stop and call it infeasible
    dual_iters::Int = 20000     # away-step Frank-Wolfe steps on the min-norm point
    gap_tol::Float64 = 1e-10    # stop when the FW duality gap ||w||^2 - min_i a_i.w is small
end

mutable struct ConeLPResult
    z::Vector{Float64}        # witness, ||z||_2 = 1 (empty if there were no rows)
    lambda::Vector{Float64}   # dual weights on the simplex
    margin_l2::Float64        # min_i a_i.z with ||z||_2 = 1   -- rigorous LOWER bound
    margin_inf::Float64       # min_i a_i.z / ||z||_inf        -- the spec's quantity
    dual_bound::Float64       # ||A^T lambda||_2               -- rigorous UPPER bound
    n_active::Int             # rows within 1e-6 of the witness minimum
    n_dual_support::Int       # rows with lambda_i > 1e-6 / m
    iterations::Int
    feasible::Bool            # margin_l2 > tol
    gap::Float64              # dual_bound - margin_l2, the width of the rigorous bracket
end
ConeLPResult() = ConeLPResult(Float64[], Float64[], 0.0, 0.0, 0.0, 0, 0, 0, false, 0.0)

"""
    normalise_rows!(A) -> n_degenerate

Scales every row of `A` to unit Euclidean norm in place; rows of norm <= 0 (or not
finite) are left alone and counted in the return value.

Replicated C++ behaviour, flagged: a row whose exact value on the flex space is zero
(a corner incidence whose two copies never separate) arrives as ~1e-16 rounding noise,
passes the `n > 0` test and is scaled to a unit row of noise. On tilings with such rows
(squares checkerboard: 80 of 160; truncated square: 16 of 55) the LP value depends on
the rounding of the flex basis and is not reproducible across linear-algebra backends
(see test_method_3.jl, frozen reference). Not fixed here because the C++ numbers are
the acceptance criterion; a threshold relative to max |A N| would be the fix.
"""
function normalise_rows!(A::Matrix{Float64})
    bad = 0
    for i in 1:size(A, 1)
        n = norm(@view A[i, :])
        if n > 0 && isfinite(n)
            A[i, :] ./= n
        else
            bad += 1
        end
    end
    return bad
end

"""
    farkas_residual(A, lambda) -> (resid, lambda_min, sum_err)

Independent verification of a Farkas / Gordan certificate. Recomputes ||A^T lambda||_2
from scratch (no solver state) and checks that lambda is in the simplex. `A` must be
the SAME matrix the certificate was computed on.
"""
function farkas_residual(A::Matrix{Float64}, lambda::Vector{Float64})
    (length(lambda) != size(A, 1) || size(A, 1) == 0) && return (Inf, -1.0, 1.0)
    # Recomputed from scratch, row by row, with no reuse of any solver iterate.
    w = zeros(size(A, 2))
    lmin = Inf
    lsum = 0.0
    for i in 1:size(A, 1)
        li = lambda[i]
        lmin = min(lmin, li)
        lsum += li
        li != 0 && (w .+= li .* @view A[i, :])
    end
    return norm(w), lmin, abs(lsum - 1.0)
end

"""
    cone_lp(A, opt = ConeLPOptions()) -> ConeLPResult

`A` must already have unit-norm rows (use `normalise_rows!`).
"""
function cone_lp(A::Matrix{Float64}, opt::ConeLPOptions = ConeLPOptions())
    r = ConeLPResult()
    m, n = size(A)
    (m == 0 || n == 0) && return r

    # Largest singular value of A, by power iteration on A^T A.
    v = fill(1.0 / sqrt(n), n)
    smax = 1.0
    for _ in 1:30
        w = A' * (A * v)
        nw = norm(w)
        nw <= 0 && break
        v = w / nw
        smax = sqrt(nw)
    end
    (!(smax > 0) || !isfinite(smax)) && (smax = 1.0)

    z = zeros(n)
    lam = fill(1.0 / m, m)
    best_primal = 0.0
    best_dual = Inf
    best_z = copy(z)
    best_lam = copy(lam)

    iters = 0
    stop = false
    mu = opt.mu0
    while mu >= opt.mu_min * 0.999 && !stop
        step = mu / (smax * smax)
        for _ in 1:opt.iters_per_stage
            a = A * z                         # m
            amin = minimum(a)
            # Primal bound, free: the witness z rescaled to the unit sphere.
            zn = norm(z)
            if zn > 0 && amin / zn > best_primal
                best_primal = amin / zn
                best_z = z / zn
            end
            # lambda = softmax(-a / mu)
            e = exp.(-(a .- amin) ./ mu)
            se = sum(e)
            if se > 0
                lam = e ./ se
            else
                fill!(lam, 1.0 / m)
            end
            gr = A' * lam                     # n
            # Dual bound, free: ||A^T lambda||_2 for a lambda in the simplex.
            d = norm(gr)
            if d < best_dual
                best_dual = d
                best_lam = copy(lam)
            end
            z .+= step .* gr
            zn2 = norm(z)
            zn2 > 1.0 && (z ./= zn2)
            iters += 1
            if best_dual < opt.dual_tol
                stop = true
                break
            end
        end
        mu *= opt.mu_factor
    end
    # ---- the min-norm point of conv{rows}, by AWAY-STEP Frank-Wolfe ----------------
    #
    # min_{l in simplex} 0.5 ||A^T l||^2.  Write w = A^T l.  Every iterate l is in the
    # simplex, so ||w|| is a rigorous UPPER bound on the LP value; and w/||w|| is a
    # feasible unit vector, so min_i a_i . w / ||w|| is a rigorous LOWER bound. At the
    # minimiser a_i . w >= ||w||^2 for every i (otherwise the segment from w to a_i would
    # contain a shorter point), so the two bounds meet and the bracket closes. Vanilla FW
    # converges at O(1/t) and stalls at the "zig-zag" boundary; the away step restores
    # linear convergence (Lacoste-Julien and Jaggi 2015).
    #
    # The ORIGINAL solver never evaluated w/||w|| as a primal witness -- that omission,
    # not the dual, is what produced margin_l2 = 0 on strictly feasible systems.
    let
        l = fill(1.0 / m, m)
        w = A' * l                            # n
        function try_primal(wv)
            wn = norm(wv)
            wn > 0 || return
            mn = minimum(A * wv) / wn
            if mn > best_primal
                best_primal = mn
                best_z = wv / wn
            end
        end
        try_primal(w)
        for t in 0:opt.dual_iters-1
            wn = norm(w)
            if wn < best_dual
                best_dual = wn
                best_lam = copy(l)
            end
            wn <= 0 && break
            cg = A * w                        # m: c_i = a_i . w
            is = argmin(cg)                   # Frank-Wolfe (toward) vertex, first minimum
            ia = 0
            cmax = -Inf
            for i in 1:m
                if l[i] > 0 && cg[i] > cmax
                    cmax = cg[i]
                    ia = i
                end
            end
            w2 = wn * wn
            gap_fw = w2 - cg[is]                        # >= 0
            gap_aw = ia > 0 ? cmax - w2 : -1.0          # >= 0
            # The FW gap also bounds the suboptimality of ||w||^2 / 2, and cg(is)/wn is
            # exactly the primal margin of the current witness, so track it here.
            if wn > 0 && cg[is] / wn > best_primal
                best_primal = cg[is] / wn
                best_z = w / wn
            end
            (gap_fw <= opt.gap_tol && (gap_aw <= opt.gap_tol || ia < 1)) && break
            (best_dual < opt.dual_tol && gap_fw <= opt.dual_tol * opt.dual_tol) && break

            local pick::Int     # vertex moved to/from
            local gmax::Float64 # step ceiling
            local sgn::Float64  # +1 toward a_pick, -1 away from a_pick
            if ia < 1 || gap_fw >= gap_aw
                pick = is; gmax = 1.0; sgn = +1.0
            else
                pick = ia; sgn = -1.0
                lv = l[ia]
                gmax = lv < 1.0 ? lv / (1.0 - lv) : Inf
                if !(gmax > 0)
                    pick = is; gmax = 1.0; sgn = +1.0
                end
            end
            dv = sgn .* (vec(A[pick, :]) .- w)
            dd = dot(dv, dv)
            dd > 0 || break
            gamma = -dot(w, dv) / dd
            gamma = min(gmax, max(0.0, gamma))
            gamma > 0 || break
            if sgn > 0
                l .*= (1.0 - gamma)
                l[pick] += gamma
            else
                l .*= (1.0 + gamma)
                l[pick] -= gamma
                l[pick] < 0 && (l[pick] = 0.0)   # a drop step: the vertex leaves the active set
                ls = sum(l)                      # keep l EXACTLY on the simplex after the drop
                ls > 0 && (l ./= ls)
            end
            w .+= gamma .* dv
            iters += 1
            (t & 255) == 255 && (w = A' * l)     # periodic refresh against drift
        end
        w = A' * l
        wn = norm(w)
        if wn < best_dual
            best_dual = wn
            best_lam = copy(l)
        end
        try_primal(w)
    end

    r.z = best_z
    r.lambda = best_lam
    r.margin_l2 = best_primal
    r.dual_bound = best_dual
    r.iterations = iters
    r.feasible = best_primal > opt.tol
    r.gap = best_dual - best_primal
    if !isempty(best_z)
        zi = maximum(abs, best_z)
        a = A * best_z
        amin = minimum(a)
        r.margin_inf = zi > 0 ? amin / zi : 0.0
        lim = amin + 1e-6 * max(1.0, abs(amin))
        r.n_active = count(x -> x <= lim, a)
    end
    lthr = 1e-6 / m
    r.n_dual_support = count(x -> x > lthr, best_lam)
    return r
end

# ---------------------------------------------------------------------------
# The whole experiment for one (cut structure, embedding).

Base.@kwdef mutable struct ExpansiveConeOptions
    max_passes::Int = 4          # pass 0 = "and", pass 1 = the sigma branch, then repairs
    lp::ConeLPOptions = ConeLPOptions()
    rel_tol::Float64 = 1e-10
end

Base.@kwdef mutable struct ExpansiveConeReport
    ok::Bool = false
    F::Int = 0
    n_rows::Int = 0
    n_dropped::Int = 0
    dim_ker_A::Int = 0
    dim_flex::Int = 0
    components::Int = 0
    flex_residual::Float64 = 0.0

    passes::Int = 0                 # passes actually solved
    pass_feasible::Int = -1         # first pass (0-based, as the C++) with a strictly
                                    # positive margin, -1 if none
    margin_l2::Float64 = 0.0        # best over passes
    margin_inf::Float64 = 0.0
    dual_bound::Float64 = 0.0       # dual bound of the best pass
    dual_bound_max::Float64 = 0.0   # WORST dual bound over the passes solved: every branch
                                    # chart tried has LP value <= this
    n_active::Int = 0
    n_dual_support::Int = 0
    flex::Vector{Float64} = Float64[]   # the witness, as a 3F flex vector (empty if none)

    # The uniform ray.
    sigma_residual::Float64 = 1.0   # ||R u_sigma|| / (||R|| ||u_sigma||)
    sigma_is_flex::Bool = false     # sigma_residual < 1e-8
    sigma_min_q::Float64 = 0.0      # min over split edges of q_e  (zero_plus_q)
    sigma_min_mu::Float64 = 0.0     # min over corner incidences of mu
    sigma_in_cone::Bool = false     # both minima > -tol_geom, i.e. sigma in P(X)
    sigma_n_bad_q::Int = 0
    sigma_n_bad_mu::Int = 0

    # The uniform ray evaluated directly in the LP's own coordinates, on the branch chart
    # it itself selects (pass 1). This is a rigorous LOWER bound on that chart's LP value
    # that the solver plays no part in: whenever sigma_in_cone is true it must be > 0, and
    # the solver must return margin_l2 >= it. The recheck of K8a uses exactly this as the
    # ground truth against which the solver is tested.
    sigma_chart_margin::Float64 = 0.0
    sigma_in_flex_span::Bool = false   # ||u_sigma - N N^T u_sigma|| / ||u_sigma|| < 1e-8

    # Independent verification of the dual certificates, recomputed by `farkas_residual`
    # from the returned multipliers on each pass's own row subset.
    farkas_resid_max::Float64 = 0.0    # worst recomputed ||A^T lambda|| over the passes solved
    farkas_lambda_min::Float64 = 0.0   # worst (most negative) multiplier over the passes
    farkas_sum_err::Float64 = 0.0      # worst |sum lambda - 1| over the passes
    passes_certified_1e9::Int = 0      # passes whose RECOMPUTED residual is < 1e-9
    passes_certified_1e6::Int = 0      # ... < 1e-6
    lp_gap_max::Float64 = 0.0          # worst bracket width dual_bound - margin_l2
end

"""
    expansive_cone(c, X, opt = ExpansiveConeOptions()) -> ExpansiveConeReport

The branch schedule of the header: pass 0 the conjunctive relaxation, pass 1 the branch
the uniform ray selects, later passes repaired from the best witness so far. Every
feasible pass is a sound existence proof; infeasibility of the passes tried proves nothing.
"""
function expansive_cone(c::CutStructure, X::Vector{Vec2},
                        opt::ExpansiveConeOptions = ExpansiveConeOptions())
    rep = ExpansiveConeReport()
    m = c.mesh
    rep.F = n_faces(m)

    g = build_hinge_graph(c)
    pins = pins_flat(c, g, X)
    fb = flex_basis(g, pins, opt.rel_tol)
    rep.dim_ker_A = fb.dim_ker_A
    rep.dim_flex = dim(fb)
    rep.components = g.components
    rep.flex_residual = fb.residual
    dim(fb) == 0 && return rep

    cs = cone_system(c, X)
    rep.n_rows = size(cs.A, 1)
    if rep.n_rows == 0
        # No split edges and no multi-copy vertices: nothing to separate, the cone is
        # everything. Report the whole flex space as feasible with an infinite margin.
        rep.ok = true
        rep.pass_feasible = 0
        rep.margin_l2 = rep.margin_inf = Inf
        return rep
    end

    M = Matrix(cs.A * fb.N)   # n_rows x dim_flex (sparse * dense)
    rep.n_dropped = normalise_rows!(M)

    # The branch schedule (see the header). Pass 0 is the conjunctive relaxation, a
    # SUBSET of P(X) that needs no branch choice at all. Pass 1 is the branch the uniform
    # ray itself selects at every convex incidence -- "the side the copy is already on",
    # read off the sigma velocities, which is the only geometrically distinguished branch
    # available at a flat state where all copies of a vertex coincide. Later passes repair
    # the branch from the best witness so far. Every pass is a SOUND existence proof when
    # feasible; infeasibility of the passes tried proves nothing.
    u_sigma = sigma_flex(c, X)
    z_sigma = fb.N' * u_sigma
    a_sigma = M * z_sigma
    let un = norm(u_sigma)
        rep.sigma_in_flex_span = (un > 0) && (norm(u_sigma - fb.N * z_sigma) / un < 1e-8)
    end

    # Index of the G2 partner of a convex G1 row, or 0.
    function partner(i::Int)
        (cs.rows[i].kind != ConeCornerG1 || !cs.rows[i].convex) && return 0
        j = i + 1
        j > rep.n_rows && return 0
        (cs.rows[j].kind != ConeCornerG2 || cs.rows[j].index != cs.rows[i].index) && return 0
        return j
    end
    function branch_from(a::Vector{Float64})
        k = trues(rep.n_rows)
        for i in 1:rep.n_rows-1
            j = partner(i)
            j == 0 && continue
            if a[i] >= a[j]
                k[j] = false
            else
                k[i] = false
            end
        end
        return k
    end

    # The uniform ray evaluated on ITS OWN branch chart, with no solver in the loop: a
    # rigorous lower bound on that chart's LP value. If sigma is inside P(X) this must be
    # positive, and the solver's margin on pass 1 must be at least this. Used by the K8a
    # recheck as ground truth (results/kill/k8a/recheck.md).
    if length(a_sigma) == rep.n_rows && norm(z_sigma) > 0
        ks = branch_from(a_sigma)
        mn = Inf
        for i in 1:rep.n_rows
            ks[i] && (mn = min(mn, a_sigma[i]))
        end
        isfinite(mn) && (rep.sigma_chart_margin = mn / norm(z_sigma))
    end

    keep = trues(rep.n_rows)
    best = ConeLPResult()
    have_best = false
    rep.farkas_lambda_min = Inf
    for pass in 0:opt.max_passes-1
        idx = findall(keep)
        isempty(idx) && break
        Mk = M[idx, :]
        r = cone_lp(Mk, opt.lp)
        rep.dual_bound_max = max(rep.dual_bound_max, r.dual_bound)
        rep.lp_gap_max = max(rep.lp_gap_max, r.gap)
        rep.passes += 1
        # Independent verification of this pass's Farkas / Gordan certificate: recompute
        # ||Mk^T lambda|| from the multipliers alone, and check lambda is on the simplex.
        let
            fres, lmin, serr = farkas_residual(Mk, r.lambda)
            rep.farkas_resid_max = max(rep.farkas_resid_max, fres)
            rep.farkas_lambda_min = min(rep.farkas_lambda_min, lmin)
            rep.farkas_sum_err = max(rep.farkas_sum_err, serr)
            if lmin >= -1e-12 && serr <= 1e-9
                fres < 1e-9 && (rep.passes_certified_1e9 += 1)
                fres < 1e-6 && (rep.passes_certified_1e6 += 1)
            end
        end
        if !have_best || r.margin_l2 > best.margin_l2
            have_best = true
            best = r
            r.feasible && (rep.pass_feasible = pass)
            rep.margin_l2 = r.margin_l2
            rep.margin_inf = r.margin_inf
            rep.dual_bound = r.dual_bound
            rep.n_active = r.n_active
            rep.n_dual_support = r.n_dual_support
            rep.flex = fb.N * r.z
        end
        r.feasible && break
        cs.n_convex == 0 && break  # no disjunction to repair
        local next::BitVector
        if pass == 0 && length(a_sigma) == rep.n_rows
            next = branch_from(a_sigma)
        elseif !isempty(r.z)
            next = branch_from(M * r.z)
        else
            break
        end
        next == keep && break
        keep = next
    end
    rep.ok = true
    isfinite(rep.farkas_lambda_min) || (rep.farkas_lambda_min = 0.0)

    # The uniform ray.
    us = u_sigma
    if n_edges(g) > 0
        R = build_rigidity(g, pins)
        rn = maximum(abs, R)
        un = norm(us)
        (rn > 0 && un > 0) && (rep.sigma_residual = norm(R * us) / (rn * un))
    else
        rep.sigma_residual = 0.0
    end
    rep.sigma_is_flex = rep.sigma_residual < 1e-8
    q = zero_plus_q(c, X)
    mu = zero_plus_corner_margin(c, X)
    scale = 0.0
    for v in q
        scale = max(scale, abs(v))
    end
    for v in mu
        scale = max(scale, abs(v))
    end
    gtol = 1e-9 * max(1.0, scale)
    rep.sigma_min_q = isempty(q) ? Inf : minimum(q)
    rep.sigma_min_mu = isempty(mu) ? Inf : minimum(mu)
    rep.sigma_n_bad_q = count(v -> v < -gtol, q)
    rep.sigma_n_bad_mu = count(v -> v < -gtol, mu)
    rep.sigma_in_cone = (rep.sigma_n_bad_q == 0 && rep.sigma_n_bad_mu == 0)
    return rep
end

"""
    euler_step(c, X, flex, h) -> Vector{Vec2}

Euler step along a flex: Y_p = X[orig(p)] + h * V_{face(p)}(X[orig(p)]) for every
M'-vertex p. Well defined because two faces sharing an M'-vertex are hinged at it.
"""
function euler_step(c::CutStructure, X::Vector{Vec2}, flex::Vector{Float64}, h::Float64)
    pv_face = _prime_face_map(c)
    Y = fill(Vec2(0, 0), c.n_prime_vertices)
    for p in 1:c.n_prime_vertices
        v = c.prime_to_original[p]
        f = pv_face[p]
        Y[p] = f >= 1 ? X[v] + h * flex_velocity(flex, f, X[v]) : X[v]
    end
    return Y
end

"""max_p |V_p| over the M'-vertices, the step-size scale."""
function flex_speed(c::CutStructure, X::Vector{Vec2}, flex::Vector{Float64})
    pv_face = _prime_face_map(c)
    s = 0.0
    for p in 1:c.n_prime_vertices
        f = pv_face[p]
        f < 1 && continue
        s = max(s, norm(flex_velocity(flex, f, X[c.prime_to_original[p]])))
    end
    return s
end
