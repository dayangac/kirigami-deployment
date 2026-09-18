# method/zero_plus.jl -- the 0+ (immediate) separation calculus of a split cut, and the
# feasibility problem "make every split cut open OUTWARD and every face positively
# oriented, somewhere in the Tutte auxetic null space".
#
# By T1.B the two copies of a split edge e = {a, b} between faces f, g (sigma_f = sigma_g)
# are the SAME vector at every theta and differ by a pure translation, so with
# Y = cos(theta/2) C + sin(theta/2) S,
#
#     y_{(v,g)}(theta) - y_{(v,f)}(theta) = sin(theta/2) * dS_e ,
#     dS_e := S_{(v,g)} - S_{(v,f)}   (the same vector for v = a and v = b),
#
# with dC_e = 0 identically. So the relative motion at 0+ is one fixed direction and the
# only question is which SIDE of e it points to.
#
# THE SIGN CONVENTION. Take the half-edge of e inside f in the face's STORED (geometrically
# CCW) order, from -> to, and let d_e := x_to - x_from. Face f is then to the LEFT of d_e and
# face g to the RIGHT. The copies separate iff g's copy moves to g's own side:
#
#     q_e := det(dS_e, d_e) = -det(d_e, dS_e) > 0    <=>   the cut OPENS
#     q_e <= 0                                       <=>   the copies move INTO each
#                                                          other and Theta_max = 0.
#
# Swapping the roles of f and g reverses BOTH d_e and dS_e, and det(-u, -v) = det(u, v), so
# q_e does not depend on which face is called f.
#
# Relation to 2026 Eq. (9) and T5.3: with the T5.3 coefficients (p, q, r) of the separation
# harmonic, r = h'(0) is exactly Eq. (9)'s first-order separation term, and q_e is a
# positive multiple of r.
#
# Everything here is exact rather than sampled: deploy() is linear and homogeneous in X, so
# dS_e(X0 + Phi t) = dS_e(X0) + sum_j t_j dS_e(Phi_j) and d_e is affine in t; hence q_e is a
# QUADRATIC polynomial in t, as is every face signed area.
#
# Shape-space coefficient layout: Phi is N x k and acts on each coordinate
# separately, so t has length m = 2k ordered (column 1 x, column 1 y, column 2 x, ...):
# t[2j-1] is the x-coefficient of column j, t[2j] its y-coefficient.

# _det2 lives in core/mesh.jl (shared).

# softplus(z) = log(1 + e^z), stable; sigma(z) = softplus'(z).
function _softplus(z::Float64)
    z > 30.0 && return z
    z < -30.0 && return exp(z)
    return log1p(exp(z))
end
function _sigmoid(z::Float64)
    z >= 0 && return 1.0 / (1.0 + exp(-z))
    e = exp(z)
    return e / (1.0 + e)
end

# S = 2 dY/dtheta at theta = 0, as a flat 2*n_prime vector (x1, y1, x2, y2, ...).
function _s_vector(c::CutStructure, X::Vector{Vec2})
    d = deploy(c, X, 0.0)
    s = Vector{Float64}(undef, 2 * length(d.dY_dtheta))
    for i in eachindex(d.dY_dtheta)
        s[2i-1] = 2.0 * d.dY_dtheta[i][1]
        s[2i] = 2.0 * d.dY_dtheta[i][2]
    end
    return s
end

# One split edge, resolved to the two M'-vertex copies whose separation q_e measures.
mutable struct SplitCopies
    edge::Int              # index into Mesh.edges
    f0::Int                # the two incident faces (f0 supplies the stored direction)
    f1::Int
    v_from::Int            # d_e = X[v_to] - X[v_from], in f0's stored CCW order
    v_to::Int
    pv0::Int               # the copies of v_from in f0 and in f1
    pv1::Int
end

function split_copies(c::CutStructure)
    m = c.mesh
    out = SplitCopies[]
    for e in c.split_edges
        ed = m.edges[e]
        ed.n_faces != 2 && continue  # a split cut always has two faces
        h0 = m.half_edges[ed.he[1]]
        f0 = h0.face
        f1 = m.half_edges[ed.he[2]].face
        push!(out, SplitCopies(e, f0, f1, h0.from, h0.to,
                               prime_vertex(c, f0, h0.from), prime_vertex(c, f1, h0.from)))
    end
    return out
end

"""q_e at a single embedding X (no shape space). Same order as `split_copies`."""
function zero_plus_q(c::CutStructure, X::Vector{Vec2})
    sc = split_copies(c)
    d = deploy(c, X, 0.0)
    q = Vector{Float64}(undef, length(sc))
    for (i, s) in enumerate(sc)
        dS = 2.0 * (d.dY_dtheta[s.pv1] - d.dY_dtheta[s.pv0])
        de = X[s.v_to] - X[s.v_from]
        q[i] = _det2(dS, de)
    end
    return q
end

# ---------------------------------------------------------------------------
# The 0+ CORNER (vertex-into-edge) separation, K6 pass 3.
#
# At theta = 0 every copy of a source vertex v sits at X_v, so for two copies p, a of v
# (a in face f, p in another face)
#
#     y_p(theta) - y_a(theta) = sin(theta/2) dS ,        dS := S_p - S_a
#
# exactly, and the two edges of f at v have directions e1, e2 (to the next and the previous
# vertex of f in stored CCW order) to leading order in theta. So p enters the interior of f
# at 0+ iff dS lies strictly inside the corner cone of f at v. Writing
#     g1 = det(e1, dS),  g2 = det(dS, e2),  cross = det(e1, e2),
# the cone is {g1 > 0 and g2 > 0} at a convex corner (cross > 0) and its complement at a
# reflex corner, so the SEPARATION MARGIN, positive exactly when p stays out of f, is
#
#     mu = max(-g1, -g2)   (convex corner) ,     mu = min(-g1, -g2)   (reflex corner).
#
# Like q_e, mu is QUADRATIC in the shape-space coefficients t; unlike q_e it is only
# piecewise smooth. FEASIBILITY is always decided by the exact margin, never by the
# smoothed objective.
mutable struct CornerIncidence
    v::Int             # source vertex of M
    face::Int          # the face whose corner at v is being entered
    pv_corner::Int     # M'-vertex: the corner of `face` at v
    pv_other::Int      # M'-vertex: the other copy of v that might enter
    v_next::Int        # e1 = X[v_next] - X[v]
    v_prev::Int        # e2 = X[v_prev] - X[v]
end

function corner_incidences(c::CutStructure)
    m = c.mesh
    out = CornerIncidence[]
    for v in 1:n_vertices(m)
        fs = m.vertex_faces[v]
        length(fs) < 2 && continue
        # The distinct copies of v over the faces around it. Hinge-joined faces share the
        # same M'-vertex and so contribute one copy, not two.
        copies = Int[]
        for f in fs
            pv = prime_vertex(c, f, v)
            (pv >= 1 && !(pv in copies)) && push!(copies, pv)
        end
        length(copies) < 2 && continue
        for f in fs
            vs = m.faces[f]
            n = length(vs)
            i = findfirst(==(v), vs)
            i === nothing && continue
            pv_a = prime_vertex(c, f, v)
            for pv_p in copies
                pv_p == pv_a && continue
                push!(out, CornerIncidence(v, f, pv_a, pv_p, vs[mod1(i + 1, n)], vs[mod1(i - 1, n)]))
            end
        end
    end
    return out
end

# mu from the three vectors, with the convexity branch. Shared by the direct evaluator and
# the objective so the two can never disagree about the branch.
function _corner_margin(e1::Vec2, e2::Vec2, dS::Vec2)
    g1 = _det2(e1, dS)
    g2 = _det2(dS, e2)
    convex = _det2(e1, e2) > 0
    return convex ? max(-g1, -g2) : min(-g1, -g2)
end

"""mu for every incidence at a single embedding X, in `corner_incidences` order."""
function zero_plus_corner_margin(c::CutStructure, X::Vector{Vec2})
    ci = corner_incidences(c)
    d = deploy(c, X, 0.0)
    mu = Vector{Float64}(undef, length(ci))
    for (i, z) in enumerate(ci)
        e1 = X[z.v_next] - X[z.v]
        e2 = X[z.v_prev] - X[z.v]
        dS = 2.0 * (d.dY_dtheta[z.pv_other] - d.dY_dtheta[z.pv_corner])
        mu[i] = _corner_margin(e1, e2, dS)
    end
    return mu
end

# The quadratic form of q over the shape space X(t) = X0 + Phi_2 t:
#
#   dS_e(t) = GS[2e-1, :], GS[2e, :]  dotted with  w = (1, t)
#   d_e(t)  = DD[2e-1, :], DD[2e, :]  dotted with  w
#   q_e(t)  = det( dS_e(t), d_e(t) )
mutable struct ZeroPlusForm
    m::Int                    # = 2 * size(Phi, 2)
    ns::Int                   # number of split edges
    GS::Matrix{Float64}       # 2*ns x (m+1)
    DD::Matrix{Float64}       # 2*ns x (m+1)
    split::Vector{SplitCopies}
end

# The unit shape-space direction (column j of Phi on coordinate `comp` = 1 (x) or 2 (y)),
# as a point set; shared by the S-column builders below.
function _phi_direction(Phi::Matrix{Float64}, j::Int, comp::Int)
    N = size(Phi, 1)
    P = Vector{Vec2}(undef, N)
    for v in 1:N
        P[v] = comp == 1 ? Vec2(Phi[v, j], 0.0) : Vec2(0.0, Phi[v, j])
    end
    return P
end

function zero_plus_form(c::CutStructure, X0::Vector{Vec2}, Phi::Matrix{Float64})
    split = split_copies(c)
    ns = length(split)
    k = size(Phi, 2)
    m = 2k
    GS = zeros(2ns, m + 1)
    DD = zeros(2ns, m + 1)

    function fill_col!(col::Int, P::Vector{Vec2})
        s = _s_vector(c, P)
        for i in 1:ns
            sp = split[i]
            GS[2i-1, col] = s[2sp.pv1-1] - s[2sp.pv0-1]
            GS[2i, col] = s[2sp.pv1] - s[2sp.pv0]
            de = P[sp.v_to] - P[sp.v_from]
            DD[2i-1, col] = de[1]
            DD[2i, col] = de[2]
        end
    end

    fill_col!(1, X0)
    for j in 1:k, comp in 1:2
        # column 2j is the x-direction of Phi column j, column 2j+1 its y-direction
        fill_col!(2j + comp - 1, _phi_direction(Phi, j, comp))
    end
    return ZeroPlusForm(m, ns, GS, DD, split)
end

"""q_e(t) for every split edge, from the form. `t` has length form.m."""
function zero_plus_eval(form::ZeroPlusForm, t::Vector{Float64})
    w = vcat(1.0, t)
    g = form.GS * w
    d = form.DD * w
    q = Vector{Float64}(undef, form.ns)
    for i in 1:form.ns
        q[i] = g[2i-1] * d[2i] - g[2i] * d[2i-1]
    end
    return q
end

"""X(t) = X0 + Phi_2 t."""
function shape_point(X0::Vector{Vec2}, Phi::Matrix{Float64}, t::Vector{Float64})
    X = copy(X0)
    k = size(Phi, 2)
    length(t) != 2k && return X
    for j in 1:k
        tx = t[2j-1]
        ty = t[2j]
        (tx == 0.0 && ty == 0.0) && continue
        for v in eachindex(X)
            X[v] += Vec2(tx * Phi[v, j], ty * Phi[v, j])
        end
    end
    return X
end

# ---------------------------------------------------------------------------
# The feasibility problem.
#
#   F(t) = sum_split softplus(-q_e(t)/s) + sum_faces softplus(-a_f(t)/s) + lambda |t|^2
#
# with s = (median edge length)^2 (both q and a are areas) and lambda = 1e-6 / s, so the
# regulariser is ~1e-6 when |t| is one median edge. Minimised by the project's own L-BFGS
# from t = 0 and `n_random` Gaussian starts of scale 0.1 * median edge.
#
# FEASIBLE means all q_e > 0 AND all a_f > 0 STRICTLY at the returned minimiser. The
# softplus never reaches zero, so feasibility is decided by the constraint values, never
# by the objective value.
Base.@kwdef mutable struct ZeroPlusRepairOptions
    n_random::Int = 3
    max_iter::Int = 300
    lambda_rel::Float64 = 1e-6
    start_scale::Float64 = 0.1   # in median edge lengths
    seed::UInt32 = 6000
    # K6 pass 3. w_corner > 0 adds sum_incidences softplus(-mu/s) to the objective (the
    # "split + vertex-edge" variant) and makes FEASIBLE additionally require every mu > 0.
    # w_prox > 0 adds w_prox * |X(t) - X0|^2 / (N * med^2), the PROXIMITY term: unlike the
    # |t|^2 regulariser it is the mesh-metric distance actually moved, so it is comparable
    # between graphs and is what the reported t_norm measures.
    w_corner::Float64 = 0.0
    w_prox::Float64 = 0.0
    # Warm start. Empty means t = 0.
    t_init::Vector{Float64} = Float64[]
end

Base.@kwdef mutable struct ZeroPlusRepairResult
    t::Vector{Float64} = Float64[]
    X::Vector{Vec2} = Vec2[]
    feasible::Bool = false
    min_q::Float64 = 0.0             # at the returned t
    min_area::Float64 = 0.0
    min_q_start::Float64 = 0.0       # at t = 0
    n_bad_q_start::Int = 0           # split edges with q <= 0 at t = 0
    n_bad_area_start::Int = 0
    n_bad_q::Int = 0                 # at the returned t
    n_bad_area::Int = 0
    min_margin::Float64 = 0.0        # min corner margin mu at the returned t (w_corner > 0)
    min_margin_start::Float64 = 0.0  # ... at t = 0
    n_bad_margin::Int = 0
    n_bad_margin_start::Int = 0
    n_incidences::Int = 0
    best_start::Int = -1             # 0 = t = 0, >0 = random start index
    iterations::Int = 0
    f_start::Float64 = 0.0
    f_end::Float64 = 0.0
end

# The exact constraint values one `_repair_eval` produces alongside F.
mutable struct _RepairMeasure
    min_q::Float64
    min_a::Float64
    n_bad_q::Int
    n_bad_a::Int
    min_mu::Float64
    n_bad_mu::Int
end
_RepairMeasure() = _RepairMeasure(0.0, 0.0, 0, 0, 0.0, 0)

# The 0+ repair objective, as a reusable context: the expensive setup (the split form,
# and the S-columns of the corner term) is built once and every evaluation reuses the
# same scratch. zero_plus_repair minimises `_repair_eval`; zero_plus_objective exposes exactly
# the same `_repair_eval` so a test can finite-difference the gradient the optimiser actually uses.
mutable struct _RepairCtx
    c::CutStructure
    X0::Vector{Vec2}
    Phi::Matrix{Float64}
    opt::ZeroPlusRepairOptions
    mesh::Mesh
    N::Int
    k::Int
    m::Int
    form::ZeroPlusForm
    s::Float64
    lambda::Float64
    use_corner::Bool
    inc::Vector{CornerIncidence}
    mode::Vector{Int}   # 0 = enforce -g1 > 0, 1 = enforce -g2 > 0, 2 = both
    n_inc::Int
    np::Int
    Smat::Matrix{Float64}
    uvec::Vector{Float64}
    vvec::Vector{Float64}
    gradS::Vector{Float64}
    gradX::Matrix{Float64}
end

function _RepairCtx(c::CutStructure, X0::Vector{Vec2}, Phi::Matrix{Float64},
                    med_edge::Float64, opt::ZeroPlusRepairOptions)
    N = length(X0)
    k = size(Phi, 2)
    m = 2k
    form = zero_plus_form(c, X0, Phi)
    s = max(1e-300, med_edge * med_edge)
    lambda = opt.lambda_rel / s

    # --- the corner (vertex-into-edge) term, only when asked for -------------------
    # S(t) is LINEAR in X and X is affine in t, so one column of S per shape-space
    # direction is the whole dependence: Smat * (1, t) is S at the point X(t).
    use_corner = opt.w_corner > 0
    inc = use_corner ? corner_incidences(c) : CornerIncidence[]
    n_inc = length(inc)
    np = c.n_prime_vertices
    Smat = zeros(0, 0)
    if use_corner
        Smat = zeros(2np, m + 1)
        Smat[:, 1] = _s_vector(c, X0)
        for j in 1:k, comp in 1:2
            Smat[:, 2j+comp-1] = _s_vector(c, _phi_direction(Phi, j, comp))
        end
    end

    # --- fix the corner ACTIVE SET at the start point ---------------------------
    # "dS outside the corner cone" is a DISJUNCTION at a convex corner
    # (max(-g1, -g2) > 0) and a CONJUNCTION at a reflex one (min(-g1, -g2) > 0), and
    # the formula switches at det(e1, e2) = 0 -- so the exact margin is discontinuous
    # in t wherever a corner passes through straightness, and a line search cannot
    # cross that. The objective therefore fixes, once, at the start point:
    #   convex corner -> the single half-plane -g_j > 0 that is closest to holding
    #                    (sufficient: it implies max(-g1, -g2) > 0),
    #   reflex corner -> both half-planes (which is the exact condition there).
    # Every term is then a plain softplus of one smooth quadratic: no max, no branch,
    # no discontinuity anywhere on the path. The set is a SUBSET of the true feasible
    # set, so a feasible answer is still a true one; and feasibility is decided by the
    # exact margin afterwards, never by this surrogate.
    mode = Int[]
    if use_corner
        mode = fill(2, n_inc)
        t0 = length(opt.t_init) == m ? opt.t_init : zeros(m)
        Xs = shape_point(X0, Phi, t0)
        s0 = Smat * vcat(1.0, t0)
        for i in 1:n_inc
            z = inc[i]
            e1 = Xs[z.v_next] - Xs[z.v]
            e2 = Xs[z.v_prev] - Xs[z.v]
            dS = Vec2(s0[2z.pv_other-1] - s0[2z.pv_corner-1], s0[2z.pv_other] - s0[2z.pv_corner])
            if _det2(e1, e2) <= 0
                mode[i] = 2  # reflex: both
                continue
            end
            mode[i] = (-_det2(e1, dS) >= -_det2(dS, e2)) ? 0 : 1
        end
    end

    return _RepairCtx(c, X0, Phi, opt, c.mesh, N, k, m, form, s, lambda, use_corner, inc, mode,
                      n_inc, np, Smat, zeros(2form.ns), zeros(2form.ns), zeros(2np), zeros(N, 2))
end

# F(t) with its gradient written into `grad` (length m); the exact constraint values go
# into `meas` when given.
function _repair_eval(ctx::_RepairCtx, t::Vector{Float64}, grad::Vector{Float64},
               meas::Union{_RepairMeasure,Nothing} = nothing)
    m = ctx.m
    form = ctx.form
    s = ctx.s
    opt = ctx.opt
    fill!(grad, 0.0)
    F = 0.0
    # --- split-edge term -----------------------------------------------------
    w = vcat(1.0, t)
    gvec = form.GS * w
    dvec = form.DD * w
    mq = Inf
    nbq = 0
    uvec = ctx.uvec
    vvec = ctx.vvec
    for i in 1:form.ns
        gx = gvec[2i-1]; gy = gvec[2i]
        dx = dvec[2i-1]; dy = dvec[2i]
        q = gx * dy - gy * dx
        mq = min(mq, q)
        q <= 0 && (nbq += 1)
        z = -q / s
        F += _softplus(z)
        wq = -_sigmoid(z) / s  # dF/dq
        # dq/dt_j = det(GS_j, d) + det(g, DD_j)
        uvec[2i-1] = wq * dy
        uvec[2i] = -wq * dx
        vvec[2i-1] = -wq * gy
        vvec[2i] = wq * gx
    end
    if form.ns > 0 && m > 0
        grad .+= transpose(@view(form.GS[:, 2:end])) * uvec
        grad .+= transpose(@view(form.DD[:, 2:end])) * vvec
    end
    # --- face-area term ------------------------------------------------------
    X = shape_point(ctx.X0, ctx.Phi, t)
    gradX = ctx.gradX
    fill!(gradX, 0.0)
    ma = Inf
    nba = 0
    mesh = ctx.mesh
    for f in 1:n_faces(mesh)
        vs = mesh.faces[f]
        n = length(vs)
        a = 0.0
        for i in 1:n
            a += _det2(X[vs[i]], X[vs[mod1(i + 1, n)]])
        end
        a *= 0.5
        ma = min(ma, a)
        a <= 0 && (nba += 1)
        z = -a / s
        F += _softplus(z)
        wa = -_sigmoid(z) / s
        for i in 1:n
            u = X[vs[mod1(i - 1, n)]] - X[vs[mod1(i + 1, n)]]
            gradX[vs[i], 1] += wa * 0.5 * (-u[2])
            gradX[vs[i], 2] += wa * 0.5 * u[1]
        end
    end
    # --- corner (vertex-into-edge) term ---------------------------------------
    # The convexity branch is read from the CURRENT point and held fixed inside the
    # evaluation: d(cross)/dt is not differentiated, so the objective is piecewise
    # smooth across a corner that changes convexity. Feasibility below is decided by
    # the EXACT margin, so a mis-taken branch can cost convergence, never soundness.
    mmu = Inf
    nbm = 0
    if ctx.use_corner
        svec = ctx.Smat * w
        gradS = ctx.gradS
        fill!(gradS, 0.0)
        for i in 1:ctx.n_inc
            z = ctx.inc[i]
            e1 = X[z.v_next] - X[z.v]
            e2 = X[z.v_prev] - X[z.v]
            dS = Vec2(svec[2z.pv_other-1] - svec[2z.pv_corner-1], svec[2z.pv_other] - svec[2z.pv_corner])
            g1 = _det2(e1, dS)
            g2 = _det2(dS, e2)
            # The EXACT margin, for reporting and for the feasibility verdict.
            convex = _det2(e1, e2) > 0
            mu = convex ? max(-g1, -g2) : min(-g1, -g2)
            mmu = min(mmu, mu)
            mu <= 0 && (nbm += 1)
            # The SURROGATE actually minimised: the fixed half-plane(s) of mode[i].
            dF_dg1 = 0.0
            dF_dg2 = 0.0
            if ctx.mode[i] != 1  # enforce -g1 > 0
                z_ = g1 / s
                F += opt.w_corner * _softplus(z_)
                dF_dg1 += opt.w_corner * _sigmoid(z_) / s
            end
            if ctx.mode[i] != 0  # enforce -g2 > 0
                z_ = g2 / s
                F += opt.w_corner * _softplus(z_)
                dF_dg2 += opt.w_corner * _sigmoid(z_) / s
            end
            # g1 = det(e1, dS), g2 = det(dS, e2)
            de1 = Vec2(dF_dg1 * dS[2], -dF_dg1 * dS[1])
            de2 = Vec2(-dF_dg2 * dS[2], dF_dg2 * dS[1])
            ddS = Vec2(dF_dg1 * (-e1[2]) + dF_dg2 * e2[2], dF_dg1 * e1[1] + dF_dg2 * (-e2[1]))
            gradX[z.v_next, 1] += de1[1]
            gradX[z.v_next, 2] += de1[2]
            gradX[z.v_prev, 1] += de2[1]
            gradX[z.v_prev, 2] += de2[2]
            gradX[z.v, 1] -= de1[1] + de2[1]
            gradX[z.v, 2] -= de1[2] + de2[2]
            gradS[2z.pv_other-1] += ddS[1]
            gradS[2z.pv_other] += ddS[2]
            gradS[2z.pv_corner-1] -= ddS[1]
            gradS[2z.pv_corner] -= ddS[2]
        end
        m > 0 && (grad .+= transpose(@view(ctx.Smat[:, 2:end])) * gradS)
    end
    # --- proximity term -------------------------------------------------------
    # |X(t) - X0|^2 / (N med^2): the mesh-metric distance actually moved, per vertex.
    if opt.w_prox > 0
        cprox = opt.w_prox / (max(1, ctx.N) * s)
        d2 = 0.0
        for v in 1:ctx.N
            dv = X[v] - ctx.X0[v]
            d2 += dv[1]^2 + dv[2]^2
            gradX[v, 1] += 2.0 * cprox * dv[1]
            gradX[v, 2] += 2.0 * cprox * dv[2]
        end
        F += cprox * d2
    end
    # The X-gradient of the area, corner and proximity terms, pulled back to t at once.
    if m > 0
        gt = transpose(ctx.Phi) * gradX  # k x 2
        for j in 1:ctx.k
            grad[2j-1] += gt[j, 1]
            grad[2j] += gt[j, 2]
        end
    end
    # --- regulariser ---------------------------------------------------------
    F += ctx.lambda * dot(t, t)
    grad .+= (2.0 * ctx.lambda) .* t
    if meas !== nothing
        meas.min_q = form.ns > 0 ? mq : 1.0
        meas.min_a = ma
        meas.n_bad_q = nbq
        meas.n_bad_a = nba
        meas.min_mu = ctx.n_inc > 0 ? mmu : 1.0
        meas.n_bad_mu = nbm
    end
    return F
end

"""
    zero_plus_objective(c, X0, Phi, med_edge, opt, t) -> (F, grad)

The repair objective F(t) and its gradient, exposed so the gradient the optimiser uses can
be finite-differenced directly. Same options, same objective as `zero_plus_repair`.
"""
function zero_plus_objective(c::CutStructure, X0::Vector{Vec2}, Phi::Matrix{Float64},
                             med_edge::Float64, opt::ZeroPlusRepairOptions, t::Vector{Float64})
    ctx = _RepairCtx(c, X0, Phi, med_edge, opt)
    g = zeros(ctx.m)
    F = _repair_eval(ctx, t, g)
    return F, g
end

"""
    zero_plus_repair(c, X0, Phi, med_edge, opt = ZeroPlusRepairOptions()) -> ZeroPlusRepairResult

Minimises the 0+ repair objective from t = 0 (or `opt.t_init`) and `opt.n_random` Gaussian
starts (`MT19937(opt.seed)`, libc++ normal_distribution draw order); stops at the first
FEASIBLE minimiser (all q_e > 0, all a_f > 0, and all mu > 0 when `opt.w_corner > 0`),
otherwise keeps the smallest F.
"""
function zero_plus_repair(c::CutStructure, X0::Vector{Vec2}, Phi::Matrix{Float64},
                          med_edge::Float64, opt::ZeroPlusRepairOptions = ZeroPlusRepairOptions())
    ctx = _RepairCtx(c, X0, Phi, med_edge, opt)
    m = ctx.m
    n_inc = ctx.n_inc
    use_corner = ctx.use_corner

    fg = (t, g) -> _repair_eval(ctx, t, g)

    best = ZeroPlusRepairResult()
    best.n_incidences = n_inc
    best.t = length(opt.t_init) == m ? copy(opt.t_init) : zeros(m)
    t_start = copy(best.t)
    let g0 = zeros(m), ms = _RepairMeasure()
        best.f_start = _repair_eval(ctx, best.t, g0, ms)
        best.min_q_start = ms.min_q
        best.min_area = ms.min_a
        best.n_bad_q_start = ms.n_bad_q
        best.n_bad_area_start = ms.n_bad_a
        best.min_margin_start = ms.min_mu
        best.n_bad_margin_start = ms.n_bad_mu
    end
    best.X = shape_point(X0, Phi, best.t)
    best.min_q = best.min_q_start
    best.n_bad_q = best.n_bad_q_start
    best.n_bad_area = best.n_bad_area_start
    best.min_margin = best.min_margin_start
    best.n_bad_margin = best.n_bad_margin_start
    best.f_end = best.f_start
    m == 0 && return best

    lo = LbfgsOptions(max_iter = opt.max_iter)
    rng = MT19937(opt.seed)
    G = NormalDist(0.0, opt.start_scale * med_edge)

    best_f = Inf
    have = false
    for start in 0:opt.n_random
        t0 = copy(t_start)
        if start > 0
            for i in 1:m
                t0[i] = t_start[i] + normal(G, rng)
            end
        end
        r = lbfgs_minimize(fg, t0, lo)
        g = zeros(m)
        ms = _RepairMeasure()
        F = _repair_eval(ctx, r.x, g, ms)
        # With the corner term on, FEASIBLE means the vertex-edge margins are strictly
        # positive too: the split signs alone were measured (K6 pass 2) not to imply it.
        feas = ms.n_bad_q == 0 && ms.n_bad_a == 0 && (!use_corner || ms.n_bad_mu == 0)
        # A feasible point always beats an infeasible one; among equals, the smaller F.
        better = !have || (feas && !best.feasible) || ((feas == best.feasible) && F < best_f)
        if better
            have = true
            best_f = F
            best.t = r.x
            best.X = shape_point(X0, Phi, r.x)
            best.feasible = feas
            best.min_q = ms.min_q
            best.min_area = ms.min_a
            best.n_bad_q = ms.n_bad_q
            best.n_bad_area = ms.n_bad_a
            best.min_margin = ms.min_mu
            best.n_bad_margin = ms.n_bad_mu
            best.best_start = start
            best.iterations = r.iterations
            best.f_end = F
        end
        best.feasible && break  # feasibility is the goal; stop at the first one found
    end
    return best
end
