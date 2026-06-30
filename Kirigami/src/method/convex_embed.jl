# method/convex_embed.jl -- the CONVEXITY-CONSTRAINED point of the Tutte auxetic shape space
# (ESCALATION.md option 2, experiment K9).
#
# Port of code/src/method/convex_embed.{hpp,cpp}.
#
# THE PROBLEM SOLVED HERE. Inside the affine shape space X(t) = X0 + Phi t (Phi acting on
# each coordinate separately, so t has dimension m = 2k as in zero_plus.jl),
#
#     minimise  || X(t) - X_ini ||_F^2      subject to   cross_i(t) >= delta  for every
#                                                        corner i of every face,
#
# where, for the corner of face f at v with the stored (geometrically CCW) neighbours
# v_prev, v_next,
#
#     a = X[v] - X[v_prev],   b = X[v_next] - X[v],   cross = det(a, b) ,
#
# so cross > 0 exactly at a strictly convex corner of a CCW face. Every corner strictly
# convex implies every face convex AND positively oriented, so this constraint SUBSUMES
# the face-area constraint of zero_plus.jl's repair.
#
# Optionally (variant (b) of the K9 task) the split-cut signs of zero_plus.jl are added
# as further constraints  q_e(t) >= delta' , with the identical sign convention.
#
# A NOTE ON THE OBJECTIVE. Eq. (6)'s X0 is the minimum-correction solution, i.e. the
# ORTHOGONAL projection of X_ini onto the affine solution set, so X_ini - X0 is orthogonal
# to range(Phi) and, with Phi orthonormal (it comes from an SVD),
#     || X(t) - X_ini ||_F^2 = || X0 - X_ini ||_F^2 + || t ||^2 .
# The code nevertheless evaluates the true || X(t) - X_ini ||_F^2 so that nothing depends
# on Phi being exactly orthonormal.
#
# HOW IT IS SOLVED. Both constraint families are QUADRATIC in t and the feasible set is
# not convex, so this is a heuristic and is reported as one:
#   phase A (feasibility): minimise  sum_i softplus((tau - cross_i)/s) [+ the split terms]
#       + lambda |t|^2  by the project's own L-BFGS, with a CONTINUATION over the target
#       tau (tau = 0.01, 0.1, 1 times the solve target, warm-started in sequence) from
#       t = 0 and n_random Gaussian starts. s = med_edge^2, the natural scale of a cross.
#   phase B (proximity): from a feasible point, minimise
#       || X(t) - X_ini ||^2 / (N s)  +  w * sum_i B((cross_i - delta)/s)  [+ split]
#       with B(u) = -log u (linearly extended below u0, so the line search never sees an
#       infinity) and w decreased geometrically. Only STRICTLY FEASIBLE iterates are kept.
#
# FEASIBLE always means the exact constraint values at the returned point, never the
# value of a penalty or a barrier. An infeasible answer is a statement about this solver,
# not a proof that the convex slice of the shape space is empty.

# The log barrier -log(u), extended linearly (C^1) below u0 so that the backtracking line
# search never evaluates an infinity and can walk back into the feasible region on its
# own. Strict feasibility is checked separately, on the exact constraint values.
# Returns (B, dB/du). Shared with range_embed.jl (the C++ duplicates it per TU).
const kBarrierU0 = 1e-8
function _log_barrier(u::Float64)
    u >= kBarrierU0 && return -log(u), -1.0 / u
    return -log(kBarrierU0) + (kBarrierU0 - u) / kBarrierU0, -1.0 / kBarrierU0
end

# One corner of one face, in the face's stored (CCW) order.
mutable struct FaceCorner
    face::Int
    idx::Int          # position of v inside Mesh.faces[face]
    v_prev::Int
    v::Int
    v_next::Int
end

function face_corners(m::Mesh)
    out = FaceCorner[]
    for f in 1:n_faces(m)
        vs = m.faces[f]
        n = length(vs)
        n < 3 && continue
        for i in 1:n
            push!(out, FaceCorner(f, i, vs[mod1(i - 1, n)], vs[i], vs[mod1(i + 1, n)]))
        end
    end
    return out
end

"""cross_i = det(X[v] - X[v_prev], X[v_next] - X[v]) for every corner, in `face_corners`
order. Strictly positive at a strictly convex corner of a CCW-stored face."""
function corner_crosses(m::Mesh, X::Vector{Vec2})
    fc = face_corners(m)
    out = Vector{Float64}(undef, length(fc))
    for (i, z) in enumerate(fc)
        a = X[z.v] - X[z.v_prev]
        b = X[z.v_next] - X[z.v]
        out[i] = _det2(a, b)
    end
    return out
end

Base.@kwdef mutable struct ConvexEmbedOptions
    # Convexity margin delta = delta_rel * med_edge^2 (a cross has units of area).
    delta_rel::Float64 = 1e-3
    # The penalty phase aims at solve_factor * delta so that feasibility at delta is not
    # decided on the boundary of what the smooth penalty can reach.
    solve_factor::Float64 = 2.0
    # > 0 turns on the split-cut constraints q_e >= split_delta_rel * med_edge^2.
    split_delta_rel::Float64 = 0.0
    n_random::Int = 3            # Gaussian restarts, in addition to t = 0
    start_scale::Float64 = 0.1   # restart scale, in median edge lengths
    seed::UInt32 = 9000
    max_iter::Int = 600          # L-BFGS iterations per continuation stage
    lambda_rel::Float64 = 1e-6   # |t|^2 regulariser of phase A, relative to med_edge^2
    # Phase B (proximity). barrier_stages = 0 disables it.
    barrier_stages::Int = 6
    barrier_w0::Float64 = 1.0
    barrier_factor::Float64 = 0.25
    barrier_iter::Int = 300
    # Warm start: the origin of phase A's continuation (and of its Gaussian restarts).
    # Empty means t = 0. Variant (a) of K9 is a RELAXATION of variant (b), so handing (a)'s
    # minimiser to (b) is a strictly better start than the origin whenever (a) was solved.
    t_init::Vector{Float64} = Float64[]
end

Base.@kwdef mutable struct ConvexEmbedResult
    t::Vector{Float64} = Float64[]
    X::Vector{Vec2} = Vec2[]
    feasible::Bool = false         # EXACT: every cross >= delta (and every q >= delta')
    feasible_strict::Bool = false  # every cross > 0 (and every q > 0), margin ignored

    n_corners::Int = 0
    n_split::Int = 0
    delta::Float64 = 0.0
    delta_split::Float64 = 0.0

    min_cross_start::Float64 = 0.0  # at t = 0, i.e. at X0
    n_bad_start::Int = 0            # corners with cross < delta at t = 0
    n_nonconvex_start::Int = 0      # corners with cross <= 0 at t = 0
    min_cross::Float64 = 0.0
    n_bad::Int = 0
    n_nonconvex::Int = 0

    min_q_start::Float64 = 0.0
    min_q::Float64 = 0.0
    n_bad_q_start::Int = 0
    n_bad_q::Int = 0

    # || X - X_ini ||_F / (sqrt(N) * med_edge): the per-vertex RMS move, in median edges.
    dist_ini::Float64 = 0.0         # at the returned point
    dist_ini_x0::Float64 = 0.0      # at X0 (t = 0)
    dist_phase_a::Float64 = 0.0     # after phase A, before the barrier refinement
    t_norm_rel::Float64 = 0.0       # |t| / med_edge

    best_start::Int = -1            # 0 = t = 0, > 0 = restart index
    iterations::Int = 0
    barrier_stages_kept::Int = 0    # barrier stages that ended strictly feasible
end

# Shared context: the corner list, the split form and all scratch, built once.
mutable struct _ConvexCtx
    c::CutStructure
    mesh::Mesh
    X0::Vector{Vec2}
    Phi::Matrix{Float64}
    X_ini::Vector{Vec2}
    opt::ConvexEmbedOptions
    N::Int
    k::Int
    m::Int
    s::Float64
    lambda::Float64
    delta::Float64
    delta_q::Float64
    use_split::Bool
    corners::Vector{FaceCorner}
    n_corner::Int
    form::ZeroPlusForm
    uvec::Vector{Float64}
    vvec::Vector{Float64}
    gradX::Matrix{Float64}
end

function _ConvexCtx(c::CutStructure, X0::Vector{Vec2}, Phi::Matrix{Float64},
                    X_ini::Vector{Vec2}, med_edge::Float64, opt::ConvexEmbedOptions)
    N = length(X0)
    k = size(Phi, 2)
    m = 2k
    s = max(1e-300, med_edge * med_edge)
    lambda = opt.lambda_rel / s
    delta = opt.delta_rel * s
    delta_q = opt.split_delta_rel * s
    use_split = opt.split_delta_rel > 0
    corners = face_corners(c.mesh)
    form = use_split ? zero_plus_form(c, X0, Phi) :
           ZeroPlusForm(m, 0, zeros(0, m + 1), zeros(0, m + 1), SplitCopies[])
    return _ConvexCtx(c, c.mesh, X0, Phi, X_ini, opt, N, k, m, s, lambda, delta, delta_q,
                      use_split, corners, length(corners), form,
                      zeros(2form.ns), zeros(2form.ns), zeros(N, 2))
end

_n_split(ctx::_ConvexCtx) = ctx.use_split ? ctx.form.ns : 0

# Exact constraint values at X(t): the ONLY thing feasibility is ever decided on.
# Returns (min_cross, n_bad, n_nonconvex, min_q, n_bad_q).
function _measure(ctx::_ConvexCtx, t::Vector{Float64})
    Xt = shape_point(ctx.X0, ctx.Phi, t)
    mc = Inf
    nb = 0
    nnc = 0
    for z in ctx.corners
        a = Xt[z.v] - Xt[z.v_prev]
        b = Xt[z.v_next] - Xt[z.v]
        cr = _det2(a, b)
        mc = min(mc, cr)
        cr < ctx.delta && (nb += 1)
        cr <= 0 && (nnc += 1)
    end
    min_cross = ctx.n_corner > 0 ? mc : 1.0
    mq = Inf
    nbq = 0
    if ctx.use_split
        q = zero_plus_q(ctx.c, Xt)
        for v in q
            mq = min(mq, v)
            v < ctx.delta_q && (nbq += 1)
        end
        isempty(q) && (mq = 1.0)
    else
        mq = 1.0
    end
    return min_cross, nb, nnc, mq, nbq
end

# ||X(t) - X_ini||_F / (sqrt(N) med): the per-vertex RMS move, in median edges.
function _dist_to_ini(ctx::_ConvexCtx, t::Vector{Float64})
    Xt = shape_point(ctx.X0, ctx.Phi, t)
    d2 = 0.0
    for v in 1:ctx.N
        dv = Xt[v] - ctx.X_ini[v]
        d2 += dv[1]^2 + dv[2]^2
    end
    return sqrt(d2 / max(1, ctx.N) / ctx.s)
end

# Accumulates the X-gradient of a corner term with dF/dcross = wc.
# cross = det(a, b), a = X[v] - X[v_prev], b = X[v_next] - X[v].
@inline function _accum_corner!(gradX::Matrix{Float64}, z::FaceCorner, a::Vec2, b::Vec2, wc::Float64)
    gradX[z.v, 1] += wc * (b[2] + a[2])
    gradX[z.v, 2] += wc * (-b[1] - a[1])
    gradX[z.v_prev, 1] += wc * (-b[2])
    gradX[z.v_prev, 2] += wc * b[1]
    gradX[z.v_next, 1] += wc * (-a[2])
    gradX[z.v_next, 2] += wc * a[1]
    return nothing
end

# grad += Phi^T gradX, interleaved (x, y) per column.
function _pullback!(grad::Vector{Float64}, Phi::Matrix{Float64}, gradX::Matrix{Float64}, k::Int)
    k == 0 && return nothing
    gt = transpose(Phi) * gradX  # k x 2
    for j in 1:k
        grad[2j-1] += gt[j, 1]
        grad[2j] += gt[j, 2]
    end
    return nothing
end

# The split-term gradient in t from per-edge weights wq_i = dF/dq_i, via the form:
# dq/dt_j = det(GS_j, d) + det(g, DD_j).
function _split_grad!(grad::Vector{Float64}, ctx::_ConvexCtx, gvec::Vector{Float64},
                      dvec::Vector{Float64}, wq_of_q)
    form = ctx.form
    F = 0.0
    for i in 1:form.ns
        gx = gvec[2i-1]; gy = gvec[2i]
        dx = dvec[2i-1]; dy = dvec[2i]
        q = gx * dy - gy * dx
        Fi, wq = wq_of_q(q)
        F += Fi
        ctx.uvec[2i-1] = wq * dy
        ctx.uvec[2i] = -wq * dx
        ctx.vvec[2i-1] = -wq * gy
        ctx.vvec[2i] = wq * gx
    end
    if form.ns > 0 && ctx.m > 0
        grad .+= transpose(@view(form.GS[:, 2:end])) * ctx.uvec
        grad .+= transpose(@view(form.DD[:, 2:end])) * ctx.vvec
    end
    return F
end

# ---- phase A: the penalty, at target fraction `frac` of the solve targets ----------
function _eval_penalty(ctx::_ConvexCtx, t::Vector{Float64}, grad::Vector{Float64}, frac::Float64)
    opt = ctx.opt
    s = ctx.s
    fill!(grad, 0.0)
    gradX = ctx.gradX
    fill!(gradX, 0.0)
    X = shape_point(ctx.X0, ctx.Phi, t)
    tau = frac * opt.solve_factor * ctx.delta
    F = 0.0
    for z in ctx.corners
        a = X[z.v] - X[z.v_prev]
        b = X[z.v_next] - X[z.v]
        cr = _det2(a, b)
        zz = (tau - cr) / s
        F += _softplus(zz)
        _accum_corner!(gradX, z, a, b, -_sigmoid(zz) / s)
    end
    if ctx.use_split
        tau_q = frac * opt.solve_factor * ctx.delta_q
        w = vcat(1.0, t)
        gvec = ctx.form.GS * w
        dvec = ctx.form.DD * w
        F += _split_grad!(grad, ctx, gvec, dvec, q -> begin
            zz = (tau_q - q) / s
            (_softplus(zz), -_sigmoid(zz) / s)
        end)
    end
    _pullback!(grad, ctx.Phi, gradX, ctx.k)
    F += ctx.lambda * dot(t, t)
    grad .+= (2.0 * ctx.lambda) .* t
    return F
end

# ---- phase B: proximity + log barrier, at barrier weight `bw` ---------------------
function _eval_barrier(ctx::_ConvexCtx, t::Vector{Float64}, grad::Vector{Float64}, bw::Float64)
    s = ctx.s
    fill!(grad, 0.0)
    gradX = ctx.gradX
    fill!(gradX, 0.0)
    X = shape_point(ctx.X0, ctx.Phi, t)
    F = 0.0
    cprox = 1.0 / (max(1, ctx.N) * s)
    for v in 1:ctx.N
        dv = X[v] - ctx.X_ini[v]
        F += cprox * (dv[1]^2 + dv[2]^2)
        gradX[v, 1] += 2.0 * cprox * dv[1]
        gradX[v, 2] += 2.0 * cprox * dv[2]
    end
    for z in ctx.corners
        a = X[z.v] - X[z.v_prev]
        b = X[z.v_next] - X[z.v]
        B, dB = _log_barrier((_det2(a, b) - ctx.delta) / s)
        F += bw * B
        _accum_corner!(gradX, z, a, b, bw * dB / s)
    end
    if ctx.use_split
        w = vcat(1.0, t)
        gvec = ctx.form.GS * w
        dvec = ctx.form.DD * w
        F += _split_grad!(grad, ctx, gvec, dvec, q -> begin
            B, dB = _log_barrier((q - ctx.delta_q) / s)
            (bw * B, bw * dB / s)
        end)
    end
    _pullback!(grad, ctx.Phi, gradX, ctx.k)
    return F
end

"""
    convex_embed_objective(c, X0, Phi, med_edge, opt, tau, t) -> (F, grad)

Phase A's objective and its gradient at `t`, with the penalty target fraction `tau`,
exposed so that the gradient the optimiser uses can be finite-differenced directly.
"""
function convex_embed_objective(c::CutStructure, X0::Vector{Vec2}, Phi::Matrix{Float64},
                                med_edge::Float64, opt::ConvexEmbedOptions, tau::Float64,
                                t::Vector{Float64})
    ctx = _ConvexCtx(c, X0, Phi, X0, med_edge, opt)
    g = zeros(ctx.m)
    F = _eval_penalty(ctx, t, g, tau)
    return F, g
end

"""
    convex_embed(c, X0, Phi, X_ini, med_edge, opt = ConvexEmbedOptions()) -> ConvexEmbedResult

The whole solve (phase A continuation from t = 0 and `opt.n_random` Gaussian restarts,
then the phase B barrier refinement). `X_ini` is the input embedding (the proximity
target); pass X0 itself to make the objective the plain |t|^2.
"""
function convex_embed(c::CutStructure, X0::Vector{Vec2}, Phi::Matrix{Float64},
                      X_ini::Vector{Vec2}, med_edge::Float64,
                      opt::ConvexEmbedOptions = ConvexEmbedOptions())
    ctx = _ConvexCtx(c, X0, Phi, X_ini, med_edge, opt)
    m = ctx.m

    r = ConvexEmbedResult()
    r.n_corners = ctx.n_corner
    r.n_split = _n_split(ctx)
    r.delta = ctx.delta
    r.delta_split = ctx.delta_q
    r.t = zeros(max(0, m))
    r.min_cross_start, r.n_bad_start, r.n_nonconvex_start, r.min_q_start, r.n_bad_q_start =
        _measure(ctx, r.t)
    r.min_cross = r.min_cross_start
    r.n_bad = r.n_bad_start
    r.n_nonconvex = r.n_nonconvex_start
    r.min_q = r.min_q_start
    r.n_bad_q = r.n_bad_q_start
    r.dist_ini_x0 = _dist_to_ini(ctx, r.t)
    r.dist_ini = r.dist_ini_x0
    r.dist_phase_a = r.dist_ini_x0
    r.X = copy(X0)
    r.feasible = (r.n_bad == 0 && r.n_bad_q == 0)
    r.feasible_strict = (r.n_nonconvex == 0 && (!ctx.use_split || r.min_q > 0))
    m == 0 && return r

    # ---- phase A: continuation over the penalty target, from t = 0 and n_random starts --
    fracs = (0.01, 0.1, 1.0)
    lo = LbfgsOptions(max_iter = opt.max_iter)
    rng = MT19937(opt.seed)
    G = NormalDist(0.0, opt.start_scale * med_edge)

    t_start = length(opt.t_init) == m ? copy(opt.t_init) : zeros(m)
    best_t = r.t
    have_feasible = false
    best_dist = Inf
    best_bad = r.n_bad + r.n_bad_q
    iterations = 0
    best_start = -1

    for start in 0:opt.n_random
        have_feasible && break
        t0 = copy(t_start)
        if start > 0
            for i in 1:m
                t0[i] = t_start[i] + normal(G, rng)
            end
        end
        for frac in fracs
            fg = (t, g) -> _eval_penalty(ctx, t, g, frac)
            lr = lbfgs_minimize(fg, t0, lo)
            t0 = lr.x
            iterations += lr.iterations
        end
        mc, nb, nnc, mq, nbq = _measure(ctx, t0)
        feas = (nb == 0 && nbq == 0)
        dist = _dist_to_ini(ctx, t0)
        better = feas ? (!have_feasible || dist < best_dist) :
                        (!have_feasible && nb + nbq < best_bad)
        if better
            best_t = t0
            best_dist = dist
            best_bad = nb + nbq
            best_start = start
            r.min_cross = mc
            r.n_bad = nb
            r.n_nonconvex = nnc
            r.min_q = mq
            r.n_bad_q = nbq
        end
        feas && (have_feasible = true)
    end
    r.best_start = best_start
    r.iterations = iterations
    r.feasible = have_feasible
    r.dist_phase_a = _dist_to_ini(ctx, best_t)
    r.dist_ini = r.dist_phase_a

    # ---- phase B: barrier refinement of the proximity, only from a feasible point -------
    if have_feasible && opt.barrier_stages > 0
        lb = LbfgsOptions(max_iter = opt.barrier_iter)
        bw = opt.barrier_w0
        t_cur = best_t
        for _ in 1:opt.barrier_stages
            fg = (t, g) -> _eval_barrier(ctx, t, g, bw)
            lr = lbfgs_minimize(fg, t_cur, lb)
            mc, nb, nnc, mq, nbq = _measure(ctx, lr.x)
            (nb != 0 || nbq != 0) && break  # left the feasible set: keep the last good point
            dist = _dist_to_ini(ctx, lr.x)
            r.barrier_stages_kept += 1
            t_cur = lr.x
            if dist < best_dist
                best_dist = dist
                best_t = lr.x
                r.min_cross = mc
                r.n_bad = nb
                r.n_nonconvex = nnc
                r.min_q = mq
                r.n_bad_q = nbq
            end
            bw *= opt.barrier_factor
        end
        r.dist_ini = best_dist
    end

    r.t = best_t
    r.X = shape_point(X0, Phi, best_t)
    r.t_norm_rel = norm(best_t) / med_edge
    r.feasible_strict = (r.n_nonconvex == 0 && (!ctx.use_split || r.min_q > 0))
    return r
end
