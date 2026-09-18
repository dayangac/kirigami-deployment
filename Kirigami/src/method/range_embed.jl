# method/range_embed.jl -- the RANGE-MAXIMISING point of the convexity-constrained slice of
# the Tutte auxetic shape space (experiment K9c).
#
# The paper's headline algorithm.
#
# WHY THIS EXISTS. K9 and K9b solve the same feasibility problem inside the null space
# X(t) = X0 + Phi t: every face corner strictly convex (cross_i >= delta) and every split
# cut opening outward (q_e >= delta'). They then select a point inside that feasible set by
# PROXIMITY, minimising ||X(t) - X_ini||. K9b's audit measured what that costs: two solvers
# found different points inside the same feasible set, one deployable and one not. Margin
# feasibility therefore does not determine deployability, and the objective that picks the
# point is what caps the yield -- not the search budget and not the feasible set.
#
# K9c replaces the objective. The barriers stay exactly as K9 posed them; what is
# maximised inside them is the DEPLOYMENT MARGIN itself.
#
# THE 0+ MARGIN. Write s = med_edge^2 (every quantity below is an area). Two families of
# first-order separation are already derived in zero_plus.jl:
#
#   split edges           q_e(X) = det(dS_e, d_e)          > 0  <=>  the cut opens at 0+
#   corner incidences     mu_j(X) = max(-g1, -g2)  (convex corner)
#                         mu_j(X) = min(-g1, -g2)  (reflex corner)
#                         g1 = det(e1, dS), g2 = det(dS, e2),
#                         e1 = X[v_next] - X[v], e2 = X[v_prev] - X[v]
#
# and the design deploys at 0+ exactly when both families are strictly positive. The
# scalar this file maximises is their joint minimum, in units of s,
#
#     m(X) := min( min_e q_e(X), min_j mu_j(X) ) / s .
#
# SMOOTHING, AND WHY THE SURROGATE IS A LOWER BOUND. m is a min of quantities that are
# QUADRATIC in t, except that mu carries the max/min disjunction whose branch is chosen by
# cross = det(e1, e2) -- the SAME cross that convex_embed.jl constrains (with
# a = X[v] - X[v_prev], b = X[v_next] - X[v], det(a, b) = det(e1, e2)). So inside the
# convexity barrier every corner is convex and the reflex branch cannot occur; the code
# still handles it, for start points that have not yet reached feasibility.
#
# The margin list is flattened into plain smooth entries c_i(t), one per split edge and,
# per corner incidence,
#
#   convex corner  ->  ONE entry -g_j, j = argmax(-g1, -g2) at the current point
#                      (a lower bound on mu = max(-g1, -g2)),
#   reflex corner  ->  TWO entries -g1 and -g2 (their min is exactly mu there),
#
# and the objective is the log-sum-exp softmin
#
#     M_kappa(t) = -kappa * log sum_i exp(-c_i(t) / (s kappa))   <=   min_i c_i(t) / s ,
#
# a LOWER BOUND on m(X(t)) for every t and every choice of modes. The reported m is always
# the EXACT min over the exact margins, never the surrogate. The modes are refreshed at the
# start of every continuation stage, and at the refresh point the convex-corner surrogate
# equals mu exactly.
#
# THE SOLVE. Stage A minimises, over t in the null space,
#
#     F(t) = -M_kappa(t) + w * [ sum_corners B((cross_i - delta)/s)
#                                + sum_split B((q_e - delta')/s) ] ,
#
# B(u) = -log u linearly extended below u0 (the same barrier convex_embed.jl uses), with w
# decreased geometrically over `stages` continuation stages and the modes refreshed between
# them. Only iterates that are EXACTLY feasible are kept; the returned point is the one with
# the largest EXACT m over all kept iterates.
#
# Stage B is `maximize_margin_range`: from a point with m > 0, push the certified range with
# range_opt.jl's softmin-of-first-contact objective under a TRUST REGION -- a decreasing
# sequence of caps on ||X - X_prev||_inf -- and accept a step only if the exact convexity
# and split constraints still hold, the exact margin has not fallen below
# `margin_floor_frac` of where it started, and the exact Theta_max improved.
#
# NOTHING HERE IS A PROOF. An infeasible or zero-margin answer is a statement about this
# solver, never about the constrained slice being empty.

"""
    zero_plus_margin(c, X, med_edge) -> (margin, min_q, min_mu)

The EXACT 0+ margin m(X) = min(min_e q_e, min_j mu_j) / med_edge^2, together with the two
families separately (all in units of med_edge^2). A design with no split edge reports
min_q = 1.0 (the +inf substitute zero_plus.jl uses); likewise for no corner incidence.
"""
function zero_plus_margin(c::CutStructure, X::Vector{Vec2}, med_edge::Float64)
    s = max(1e-300, med_edge * med_edge)
    q = zero_plus_q(c, X)
    mq = isempty(q) ? s : minimum(q)
    mu = zero_plus_corner_margin(c, X)
    mmu = isempty(mu) ? s : minimum(mu)
    return min(mq, mmu) / s, mq / s, mmu / s
end

"""
    range_embed_modes(c, X) -> Vector{Int}

The per-incidence surrogate mode at X, in `corner_incidences(c)` order:
  0 = convex corner, enforce -g1;  1 = convex corner, enforce -g2;  2 = reflex, both.
At a convex corner the mode is the argmax of (-g1, -g2), so the surrogate entry EQUALS mu
there; anywhere else it is a lower bound.
"""
function range_embed_modes(c::CutStructure, X::Vector{Vec2})
    inc = corner_incidences(c)
    d = deploy(c, X, 0.0)
    mode = fill(2, length(inc))
    for (i, z) in enumerate(inc)
        e1 = X[z.v_next] - X[z.v]
        e2 = X[z.v_prev] - X[z.v]
        dS = 2.0 * (d.dY_dtheta[z.pv_other] - d.dY_dtheta[z.pv_corner])
        if _det2(e1, e2) <= 0
            mode[i] = 2  # reflex: both half-planes
            continue
        end
        # Convex: the single half-plane closest to holding. -g_j <= max(-g1, -g2) = mu, and
        # the argmax attains it, so the surrogate entry equals mu at this very point.
        mode[i] = (-_det2(e1, dS) >= -_det2(dS, e2)) ? 0 : 1
    end
    return mode
end

Base.@kwdef mutable struct RangeEmbedOptions
    # Barriers, identical in form and units to convex_embed.jl.
    delta_rel::Float64 = 1e-3        # convexity margin delta  = delta_rel * med^2
    split_delta_rel::Float64 = 1e-3  # split margin      delta' = split_delta_rel * med^2
    # Softmin temperature, in units of med^2 (the margins are measured in med^2). Small
    # kappa tracks the true min more tightly and is harder to optimise.
    kappa::Float64 = 5e-3
    stages::Int = 8                  # barrier continuation stages (mode refresh between)
    iter_per_stage::Int = 200        # L-BFGS iterations per stage
    barrier_w0::Float64 = 1e-2       # w of stage 0
    barrier_factor::Float64 = 0.5    # w *= factor per stage
    n_random::Int = 0                # Gaussian restarts around t_init, in addition to it
    start_scale::Float64 = 0.05      # restart scale, in median edges
    seed::UInt32 = 9300
    t_init::Vector{Float64} = Float64[]  # warm start; empty means t = 0
end

Base.@kwdef mutable struct RangeEmbedResult
    t::Vector{Float64} = Float64[]
    X::Vector{Vec2} = Vec2[]
    feasible::Bool = false       # EXACT: every cross >= delta and every q >= delta'
    margin::Float64 = 0.0        # EXACT m(X) at the returned point, in med^2
    margin_start::Float64 = 0.0  # ... at the start point
    min_cross::Float64 = 0.0     # in med^2
    min_q::Float64 = 0.0
    min_mu::Float64 = 0.0
    n_bad::Int = 0
    n_bad_q::Int = 0
    n_bad_mu::Int = 0
    stages_kept::Int = 0
    iterations::Int = 0
    best_start::Int = -1
    n_entries::Int = 0           # size of the flattened margin list at the start modes
end

# Exact constraint and margin values at X(t). The ONLY thing feasibility and the reported
# margin are ever decided on.
mutable struct _RangeExact
    min_cross::Float64
    min_q::Float64
    min_mu::Float64
    margin::Float64
    n_bad::Int
    n_bad_q::Int
    n_bad_mu::Int
    feasible::Bool
end

# Everything the stage-A objective needs, built once per solve.
mutable struct _RangeCtx
    c::CutStructure
    mesh::Mesh
    X0::Vector{Vec2}
    Phi::Matrix{Float64}
    opt::RangeEmbedOptions
    N::Int
    k::Int
    m::Int
    np::Int
    s::Float64
    delta::Float64
    delta_q::Float64
    corners::Vector{FaceCorner}
    inc::Vector{CornerIncidence}
    n_corner::Int
    n_inc::Int
    form::ZeroPlusForm
    Smat::Matrix{Float64}   # 2*np x (m+1): S at X(t) is Smat * (1, t)
    # scratch
    uvec::Vector{Float64}
    vvec::Vector{Float64}
    gradS::Vector{Float64}
    gradX::Matrix{Float64}
end

function _RangeCtx(c::CutStructure, X0::Vector{Vec2}, Phi::Matrix{Float64},
                   med_edge::Float64, opt::RangeEmbedOptions)
    N = length(X0)
    k = size(Phi, 2)
    m = 2k
    np = c.n_prime_vertices
    s = max(1e-300, med_edge * med_edge)
    delta = opt.delta_rel * s
    delta_q = opt.split_delta_rel * s
    corners = face_corners(c.mesh)
    inc = corner_incidences(c)
    form = zero_plus_form(c, X0, Phi)

    # S is LINEAR in X and X is affine in t, so one column per shape-space direction is
    # the whole dependence.
    Smat = zeros(2np, m + 1)
    Smat[:, 1] = _s_vector(c, X0)
    for j in 1:k, comp in 1:2
        Smat[:, 2j+comp-1] = _s_vector(c, _phi_direction(Phi, j, comp))
    end

    return _RangeCtx(c, c.mesh, X0, Phi, opt, N, k, m, np, s, delta, delta_q, corners, inc,
                     length(corners), length(inc), form, Smat,
                     zeros(2 * max(1, form.ns)), zeros(2 * max(1, form.ns)), zeros(2np), zeros(N, 2))
end

function _n_entries(ctx::_RangeCtx, mode::Vector{Int})
    n = ctx.form.ns
    for md in mode
        n += (md == 2) ? 2 : 1
    end
    return n
end

function _measure(ctx::_RangeCtx, t::Vector{Float64})
    Xt = shape_point(ctx.X0, ctx.Phi, t)
    mc = Inf
    n_bad = 0
    for z in ctx.corners
        a = Xt[z.v] - Xt[z.v_prev]
        b = Xt[z.v_next] - Xt[z.v]
        cr = _det2(a, b)
        mc = min(mc, cr)
        cr < ctx.delta && (n_bad += 1)
    end
    min_cross = (ctx.n_corner > 0 ? mc : ctx.s) / ctx.s
    margin, mq, mmu = zero_plus_margin(ctx.c, Xt, sqrt(ctx.s))
    n_bad_q = count(v -> v < ctx.delta_q, zero_plus_q(ctx.c, Xt))
    n_bad_mu = count(v -> v <= 0, zero_plus_corner_margin(ctx.c, Xt))
    return _RangeExact(min_cross, mq, mmu, margin, n_bad, n_bad_q, n_bad_mu,
                       n_bad == 0 && n_bad_q == 0)
end

# F(t) = -softmin_kappa(margin entries / s) + bw * sum barriers, with its gradient.
function _range_eval(ctx::_RangeCtx, t::Vector{Float64}, grad::Vector{Float64},
                     mode::Vector{Int}, bw::Float64)
    m = ctx.m
    s = ctx.s
    form = ctx.form
    fill!(grad, 0.0)
    gradX = ctx.gradX
    fill!(gradX, 0.0)
    X = shape_point(ctx.X0, ctx.Phi, t)
    w = vcat(1.0, t)
    gvec = form.GS * w
    dvec = form.DD * w
    svec = ctx.Smat * w

    # ---- pass 1: collect every margin entry, in units of s ---------------------------
    ne = _n_entries(ctx, mode)
    cvals = zeros(max(1, ne))
    idx = 0
    for i in 1:form.ns
        q = gvec[2i-1] * dvec[2i] - gvec[2i] * dvec[2i-1]
        cvals[idx += 1] = q / s
    end
    for i in 1:ctx.n_inc
        z = ctx.inc[i]
        e1 = X[z.v_next] - X[z.v]
        e2 = X[z.v_prev] - X[z.v]
        dS = Vec2(svec[2z.pv_other-1] - svec[2z.pv_corner-1], svec[2z.pv_other] - svec[2z.pv_corner])
        g1 = _det2(e1, dS)
        g2 = _det2(dS, e2)
        mode[i] != 1 && (cvals[idx += 1] = -g1 / s)
        mode[i] != 0 && (cvals[idx += 1] = -g2 / s)
    end

    # ---- the log-sum-exp softmin, stably ---------------------------------------------
    kappa = max(1e-12, ctx.opt.kappa)
    M = 0.0
    wts = zeros(max(1, ne))
    if ne > 0
        cmin = minimum(@view cvals[1:ne])
        Z = 0.0
        for i in 1:ne
            e = exp(-(cvals[i] - cmin) / kappa)
            wts[i] = e
            Z += e
        end
        M = cmin - kappa * log(Z)
        wts[1:ne] ./= Z   # dM/dc_i, summing to 1
    end
    F = -M

    # ---- pass 2: the gradient of -M, plus the barriers --------------------------------
    # dF/dc_i = -wts(i); c_i = raw_i / s, so dF/draw_i = -wts(i)/s.
    uvec = ctx.uvec
    vvec = ctx.vvec
    if form.ns > 0
        fill!(uvec, 0.0)
        fill!(vvec, 0.0)
    end
    gradS = ctx.gradS
    fill!(gradS, 0.0)
    idx = 0
    for i in 1:form.ns
        gx = gvec[2i-1]; gy = gvec[2i]
        dx = dvec[2i-1]; dy = dvec[2i]
        q = gx * dy - gy * dx
        wq = -wts[idx += 1] / s                    # objective part
        B, dB = _log_barrier((q - ctx.delta_q) / s)  # split-outward barrier
        F += bw * B
        wq += bw * dB / s
        # dq/dt = det(GS_j, d) + det(g, DD_j)
        uvec[2i-1] = wq * dy
        uvec[2i] = -wq * dx
        vvec[2i-1] = -wq * gy
        vvec[2i] = wq * gx
    end
    if form.ns > 0 && m > 0
        grad .+= transpose(@view(form.GS[:, 2:end])) * @view(uvec[1:2form.ns])
        grad .+= transpose(@view(form.DD[:, 2:end])) * @view(vvec[1:2form.ns])
    end
    for i in 1:ctx.n_inc
        z = ctx.inc[i]
        e1 = X[z.v_next] - X[z.v]
        e2 = X[z.v_prev] - X[z.v]
        dS = Vec2(svec[2z.pv_other-1] - svec[2z.pv_corner-1], svec[2z.pv_other] - svec[2z.pv_corner])
        dF_dg1 = 0.0
        dF_dg2 = 0.0
        mode[i] != 1 && (dF_dg1 += wts[idx += 1] / s)   # entry is -g1, dF/dc = -w
        mode[i] != 0 && (dF_dg2 += wts[idx += 1] / s)
        # g1 = det(e1, dS): dg1/de1 = (dS.y, -dS.x), dg1/ddS = (-e1.y, e1.x)
        # g2 = det(dS, e2): dg2/de2 = (-dS.y, dS.x), dg2/ddS = ( e2.y, -e2.x)
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

    # ---- the convexity barrier --------------------------------------------------------
    for z in ctx.corners
        a = X[z.v] - X[z.v_prev]
        b = X[z.v_next] - X[z.v]
        B, dB = _log_barrier((_det2(a, b) - ctx.delta) / s)
        F += bw * B
        _accum_corner!(gradX, z, a, b, bw * dB / s)
    end

    _pullback!(grad, ctx.Phi, gradX, ctx.k)
    return F
end

"""
    range_embed_objective(c, X0, Phi, med_edge, opt, mode, barrier_w, t) -> (F, grad)

Stage A's objective and its analytic gradient at fixed modes and fixed barrier weight,
exposed so that a test can finite-difference the gradient the optimiser actually uses.
"""
function range_embed_objective(c::CutStructure, X0::Vector{Vec2}, Phi::Matrix{Float64},
                               med_edge::Float64, opt::RangeEmbedOptions, mode::Vector{Int},
                               barrier_w::Float64, t::Vector{Float64})
    ctx = _RangeCtx(c, X0, Phi, med_edge, opt)
    g = zeros(ctx.m)
    F = _range_eval(ctx, t, g, mode, barrier_w)
    return F, g
end

function _copy_exact!(r::RangeEmbedResult, e::_RangeExact)
    r.margin = e.margin
    r.min_cross = e.min_cross
    r.min_q = e.min_q
    r.min_mu = e.min_mu
    r.n_bad = e.n_bad
    r.n_bad_q = e.n_bad_q
    r.n_bad_mu = e.n_bad_mu
    r.feasible = e.feasible
    return r
end

"""
    range_embed(c, X0, Phi, med_edge, opt = RangeEmbedOptions()) -> RangeEmbedResult

Stage A: the barrier continuation that maximises the softmin surrogate of the exact 0+
margin inside the convexity and split-sign barriers, from `opt.t_init` (or t = 0) and
`opt.n_random` Gaussian restarts (`MT19937(opt.seed)`). The returned point is the FEASIBLE
kept iterate with the largest EXACT margin; if no start ever becomes feasible, the
least-infeasible one, reported as infeasible.
"""
function range_embed(c::CutStructure, X0::Vector{Vec2}, Phi::Matrix{Float64},
                     med_edge::Float64, opt::RangeEmbedOptions = RangeEmbedOptions())
    ctx = _RangeCtx(c, X0, Phi, med_edge, opt)
    m = ctx.m

    r = RangeEmbedResult()
    t_start = length(opt.t_init) == m ? copy(opt.t_init) : zeros(max(0, m))
    r.t = t_start
    let e0 = _measure(ctx, t_start)
        r.margin_start = e0.margin
        _copy_exact!(r, e0)
    end
    r.X = shape_point(X0, Phi, r.t)
    m == 0 && return r
    r.n_entries = _n_entries(ctx, range_embed_modes(c, r.X))

    lo = LbfgsOptions(max_iter = opt.iter_per_stage)
    rng = MT19937(opt.seed)
    G = NormalDist(0.0, opt.start_scale * med_edge)

    # The point kept is the FEASIBLE one with the largest EXACT margin; if no start ever
    # becomes feasible, the least-infeasible one is kept and reported as infeasible.
    have_feasible = r.feasible
    best_margin = r.feasible ? r.margin : -Inf
    best_bad = r.n_bad + r.n_bad_q

    for start in 0:opt.n_random
        t_cur = copy(t_start)
        if start > 0
            for i in 1:m
                t_cur[i] = t_start[i] + normal(G, rng)
            end
        end
        cur = _measure(ctx, t_cur)
        cur_feasible = cur.feasible
        bw = opt.barrier_w0
        for _ in 1:opt.stages
            # Refresh the surrogate branch at the current point: there the convex-corner
            # entry equals mu exactly, so the lower bound is tight where the stage begins.
            mode = range_embed_modes(c, shape_point(X0, Phi, t_cur))
            fg = (t, g) -> _range_eval(ctx, t, g, mode, bw)
            lr = lbfgs_minimize(fg, t_cur, lo)
            r.iterations += lr.iterations
            e = _measure(ctx, lr.x)
            # Once feasible, never leave: an infeasible iterate ends this start's continuation.
            (cur_feasible && !e.feasible) && break
            t_cur = lr.x
            cur = e
            cur_feasible = e.feasible
            r.stages_kept += 1
            better = e.feasible ? (!have_feasible || e.margin > best_margin) :
                                  (!have_feasible && e.n_bad + e.n_bad_q < best_bad)
            if better
                if e.feasible
                    have_feasible = true
                    best_margin = e.margin
                end
                best_bad = e.n_bad + e.n_bad_q
                r.t = lr.x
                _copy_exact!(r, e)
                r.best_start = start
            end
            bw *= opt.barrier_factor
        end
    end
    r.X = shape_point(X0, Phi, r.t)
    return r
end

# ---------------------------------------------------------------------------
# Stage B: maximise the exact Theta_max from a point with m > 0, under a trust region and
# exact rejection of any step that breaks convexity, the split signs, or the margin floor.

Base.@kwdef mutable struct MarginRangeOptions
    delta_rel::Float64 = 1e-3
    split_delta_rel::Float64 = 1e-3
    margin_floor_frac::Float64 = 0.5   # accept only if m >= frac * m at the point handed in
    step_caps::Vector{Float64} = [0.25, 0.10, 0.04]  # in median edges
    rounds::Int = 8
    iters_per_round::Int = 10
    kappa::Float64 = 0.05              # range_opt softmin temperature (radians)
    barrier_eps::Float64 = 0.2         # range_opt's face-area floor
end

Base.@kwdef mutable struct MarginRangeResult
    X::Vector{Vec2} = Vec2[]
    improved::Bool = false
    theta_before::Float64 = 0.0
    theta_after::Float64 = 0.0
    margin_after::Float64 = 0.0
    caps_tried::Int = 0
    caps_accepted::Int = 0
end

function _exact_theta(c::CutStructure, X::Vector{Vec2})
    B = deploy_basis(c, X)
    sd = swept_discs(c, B)
    pairs = candidate_pairs(c, sd, Float64(pi), true)
    return exact_theta_max_overlap(c, B, pairs, 1e-9, Float64(pi), 1e-9).theta_max
end

"""
    maximize_margin_range(c, X, Phi, med_edge, o = MarginRangeOptions()) -> MarginRangeResult

Stage B: for each cap in `o.step_caps` runs `maximize_range` under that trust region and
accepts the step only if, measured exactly at the returned point, every cross >= delta,
every q >= delta', the 0+ margin is >= `o.margin_floor_frac` of its starting value, and
Theta_max improved.
"""
function maximize_margin_range(c::CutStructure, X::Vector{Vec2}, Phi::Matrix{Float64},
                               med_edge::Float64, o::MarginRangeOptions = MarginRangeOptions())
    mesh = c.mesh
    s = max(1e-300, med_edge * med_edge)
    delta = o.delta_rel * s
    delta_q = o.split_delta_rel * s

    res = MarginRangeResult()
    res.X = copy(X)
    res.theta_before = _exact_theta(c, X)
    res.theta_after = res.theta_before
    margin0, _, _ = zero_plus_margin(c, X, med_edge)
    res.margin_after = margin0
    floor_m = o.margin_floor_frac * margin0

    for cap in o.step_caps
        res.caps_tried += 1
        rop = RangeOptOptions()
        rop.kappa = o.kappa
        rop.rounds = o.rounds
        rop.iters_per_round = o.iters_per_round
        rop.barrier_eps = o.barrier_eps
        rop.step_cap = cap * med_edge
        rr = maximize_range(c, res.X, Phi, rop)
        # Exact rejection: convexity, the split signs, and the margin floor, all measured at
        # the returned point -- range_opt knows about none of them.
        mc = minimum(corner_crosses(mesh, rr.X_opt); init = Inf)
        q = zero_plus_q(c, rr.X_opt)
        mq = isempty(q) ? s : minimum(q)
        mg, _, _ = zero_plus_margin(c, rr.X_opt, med_edge)
        (mc < delta || mq < delta_q || mg < floor_m) && continue
        th = _exact_theta(c, rr.X_opt)
        th <= res.theta_after + 1e-9 && continue
        res.X = rr.X_opt
        res.theta_after = th
        res.margin_after = mg
        res.improved = true
        res.caps_accepted += 1
    end
    return res
end
