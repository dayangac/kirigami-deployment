# method/range_opt.jl -- K2b's alternative to Eq. (6) + Eq. (9): maximise a softmin of the
# CLOSED-FORM first-contact angles over the Tutte auxetic null space, under a smooth
# barrier that keeps every face positively oriented.
#
#   maximise   softmin_kappa { theta*_i(T) : i in the active set }
#   over       X = X0 + Phi T
#   subject to (as a barrier)  area_f(X) >= eps * area_f(X0)
#
# theta*_i is the smallest root in (0, theta_hi] of the orientation harmonic of a
# (vertex, edge) pair that also satisfies the two interval harmonics, so the objective
# is differentiable wherever the active root is simple; the gradient is analytic,
# obtained by implicit differentiation of  p + q cos + r sin = 0  and pulled back
# through the two linear operators C and S of deploy_basis.jl.
#
# The active set is rebuilt from the swept-disc broad phase every `iters_per_round`
# L-BFGS iterations. Nothing here is refereed by its own objective: the caller
# measures the result with collision.jl's bisection.
#
# Port of code/src/method/range_opt.{hpp,cpp}. The C++ member functions of
# `RangeObjective` are `setup!(ob, ...)`, `basis(ob, t)`, `rebuild_active!(ob, t)` and
# `value_and_grad(ob, t, g)` (g written in place, the `lbfgs_minimize` contract).
# The flattened coefficient vector t stores T row-major: t[2i-1] = T[i,1], t[2i] = T[i,2].

Base.@kwdef mutable struct RangeOptOptions
    kappa::Float64 = 0.05          # softmin temperature (radians)
    rounds::Int = 20               # active-set refreshes
    iters_per_round::Int = 10      # L-BFGS iterations between refreshes
    active_max::Int = 250          # candidates kept per refresh
    barrier_w::Float64 = 1.0
    barrier_eps::Float64 = 0.2     # floor as a fraction of the initial signed area
    barrier_kappa::Float64 = 0.05
    theta_hi::Float64 = Float64(pi)
    theta_lo::Float64 = 1e-9
    step_cap::Float64 = 0.0        # if > 0, ||Phi T||_inf is kept below this
end

mutable struct RangeOptResult
    T::Matrix{Float64}             # k x 2
    X_opt::Vector{Vec2}
    theta_ref_before::Float64      # refereed by collision.jl bisection
    theta_ref_after::Float64
    theta_closed_before::Float64; theta_closed_after::Float64
    rounds_run::Int
    n_active::Int
    trace::Vector{Float64}         # refereed theta_max after each round
end

# Exposed for testing: the objective and its analytic gradient at a fixed active set.
mutable struct RangeObjective
    c::Union{CutStructure,Nothing}
    X0m::Matrix{Float64}        # N x 2
    Phi::Matrix{Float64}        # N x k
    U0::Matrix{Float64}; Uphi::Matrix{Float64}   # n_prime x 2, n_prime x k
    area0::Vector{Float64}      # per face of M
    active::Vector{NTuple{3,Int}}  # (pa, pb, pv) M'-vertices
    opt::RangeOptOptions
end
RangeObjective() = RangeObjective(nothing, zeros(0, 2), zeros(0, 0), zeros(0, 2), zeros(0, 0),
                                  Float64[], NTuple{3,Int}[], RangeOptOptions())

_unflatten_T(t::AbstractVector{Float64}, k::Int) = (T = Matrix{Float64}(undef, k, 2);
    for i in 1:k; T[i, 1] = t[2i - 1]; T[i, 2] = t[2i]; end; T)

function setup!(ob::RangeObjective, cc::CutStructure, X0::Vector{Vec2}, P::AbstractMatrix{Float64},
                o::RangeOptOptions)
    ob.c = cc
    ob.opt = o
    ob.Phi = Matrix{Float64}(P)
    ob.X0m = points_to_matrix(X0)
    N = length(X0); k = size(ob.Phi, 2)
    cols = Matrix{Float64}(undef, N, 2 + k)
    cols[:, 1:2] = ob.X0m
    k > 0 && (cols[:, 3:end] = ob.Phi)
    U = velocity_operator_apply(cc, cols)
    ob.U0 = U[:, 1:2]
    ob.Uphi = U[:, 3:end]
    m = cc.mesh
    ob.area0 = zeros(n_faces(m))
    for f in 1:n_faces(m)
        vs = m.faces[f]
        s = 0.0
        for i in eachindex(vs)
            a = X0[vs[i]]
            b = X0[vs[mod1(i + 1, length(vs))]]
            s += a[1] * b[2] - a[2] * b[1]
        end
        ob.area0[f] = 0.5 * s
    end
    return ob
end

"""The deploy basis at X = X0 + Phi T, from the linear operators (no deploy() call)."""
function basis(ob::RangeObjective, t::AbstractVector{Float64})
    k = size(ob.Phi, 2)
    T = _unflatten_T(t, k)
    X = k > 0 ? ob.X0m + ob.Phi * T : copy(ob.X0m)
    Um = k > 0 ? ob.U0 + ob.Uphi * T : copy(ob.U0)
    np = ob.c.n_prime_vertices
    C = Matrix{Float64}(undef, np, 2)
    S = Matrix{Float64}(undef, np, 2)
    for i in 1:np
        v = ob.c.prime_to_original[i]
        C[i, 1] = X[v, 1]
        C[i, 2] = X[v, 2]
        S[i, 1] = -2.0 * Um[i, 2]   # S = 2 * rot90(U)
        S[i, 2] = 2.0 * Um[i, 1]
    end
    return DeployBasis(C, S)
end

function rebuild_active!(ob::RangeObjective, t::AbstractVector{Float64})
    c = ob.c
    B = basis(ob, t)
    sd = swept_discs(c, B)
    pairs = candidate_pairs(c, sd, ob.opt.theta_hi, true)
    found = Tuple{Float64,NTuple{3,Int}}[]
    PF = c.prime_faces
    function scan(fe::Int, fv::Int)
        E = PF[fe]
        V = PF[fv]
        ne = length(E)
        for i in 1:ne
            a = E[i]; b = E[mod1(i + 1, ne)]
            U = basis_c(B, b) - basis_c(B, a); Vv = basis_s(B, b) - basis_s(B, a)
            L2 = dot_from_vectors(U, Vv, U, Vv)
            lscale = abs(L2.p) + amp(L2)
            lscale <= 0 && continue
            for p in V
                (p == a || p == b) && continue
                P = basis_c(B, p) - basis_c(B, a); Q = basis_s(B, p) - basis_s(B, a)
                det = orient_from_vectors(U, Vv, P, Q)
                abs(det.p) + amp(det) <= 0 && continue
                for th in harmonic_roots(det, ob.opt.theta_lo, ob.opt.theta_hi)
                    D = dot_from_vectors(U, Vv, P, Q)
                    s = harmonic_eval(D, th); l2 = harmonic_eval(L2, th)
                    tol = 1e-12 * lscale
                    (s < -tol || s > l2 + tol) && continue
                    push!(found, (th, (a, b, p)))
                    break
                end
            end
        end
    end
    for (f, g) in pairs
        scan(f, g)
        scan(g, f)
    end
    # std::sort by the angle only: ties keep no particular order in the C++ either
    sort!(found; by = first)
    empty!(ob.active)
    for i in 1:min(length(found), ob.opt.active_max)
        push!(ob.active, found[i][2])
    end
    return ob
end

"""Objective (-softmin of the active first-contact angles + barrier) and its analytic
gradient, written into `g` (length 2k). Returns the objective value."""
function value_and_grad(ob::RangeObjective, t::AbstractVector{Float64}, g::AbstractVector{Float64})
    c = ob.c
    k = size(ob.Phi, 2)
    N = size(ob.X0m, 1)
    np = c.n_prime_vertices
    opt = ob.opt
    fill!(g, 0.0)
    B = basis(ob, t)
    T = _unflatten_T(t, k)
    X = k > 0 ? ob.X0m + ob.Phi * T : copy(ob.X0m)

    # ---- softmin of the closed-form first-contact angles -----------------------
    na = length(ob.active)
    th = fill(opt.theta_hi, na)
    dFdtheta = zeros(na)
    has = falses(na)
    for i in 1:na
        a, b, p = ob.active[i]
        U = basis_c(B, b) - basis_c(B, a); Vv = basis_s(B, b) - basis_s(B, a)
        P = basis_c(B, p) - basis_c(B, a); Q = basis_s(B, p) - basis_s(B, a)
        det = orient_from_vectors(U, Vv, P, Q)
        L2 = dot_from_vectors(U, Vv, U, Vv)
        D = dot_from_vectors(U, Vv, P, Q)
        lscale = abs(L2.p) + amp(L2)
        for r in harmonic_roots(det, opt.theta_lo, opt.theta_hi)
            s = harmonic_eval(D, r); l2 = harmonic_eval(L2, r)
            tol = 1e-12 * max(1e-300, lscale)
            (s < -tol || s > l2 + tol) && continue
            th[i] = r
            dFdtheta[i] = -det.q * sin(r) + det.r * cos(r)
            has[i] = true
            break
        end
    end
    mn = opt.theta_hi
    for v in th
        mn = min(mn, v)
    end
    Z = 0.0
    for v in th
        Z += exp(-(v - mn) / opt.kappa)
    end
    Z <= 0 && (Z = 1.0)
    softmin = mn - opt.kappa * log(Z)
    obj = -softmin

    gC = zeros(np, 2); gS = zeros(np, 2)
    for i in 1:na
        has[i] || continue
        abs(dFdtheta[i]) < 1e-14 && continue
        alpha = exp(-(th[i] - mn) / opt.kappa) / Z
        alpha < 1e-12 && continue
        a, b, p = ob.active[i]
        cc = cos(0.5 * th[i]); ss = sin(0.5 * th[i])
        Ya = cc * basis_c(B, a) + ss * basis_s(B, a)
        Yb = cc * basis_c(B, b) + ss * basis_s(B, b)
        Yp = cc * basis_c(B, p) + ss * basis_s(B, p)
        A = Yb - Ya; Bv = Yp - Ya
        w = alpha / dFdtheta[i]
        gb = Vec2(Bv[2], -Bv[1])
        gp = Vec2(-A[2], A[1])
        ga = -(gb + gp)
        for (idx, gv) in ((a, ga), (b, gb), (p, gp))
            gC[idx, 1] += w * cc * gv[1]
            gC[idx, 2] += w * cc * gv[2]
            gS[idx, 1] += w * ss * gv[1]
            gS[idx, 2] += w * ss * gv[2]
        end
    end

    # ---- signed-area barrier ---------------------------------------------------
    m = c.mesh
    gX = zeros(N, 2)
    for f in 1:n_faces(m)
        vs = m.faces[f]
        nf = length(vs)
        s = 0.0
        for i in 1:nf
            u = vs[i]; v = vs[mod1(i + 1, nf)]
            s += X[u, 1] * X[v, 2] - X[u, 2] * X[v, 1]
        end
        area = 0.5 * s
        ob.area0[f] == 0 && continue
        u = area / ob.area0[f] - opt.barrier_eps
        kb = opt.barrier_kappa
        e = u / kb
        local Bv, dB
        if e > 30
            Bv = kb * exp(-e); dB = -exp(-e)
        elseif e < -30
            Bv = -u; dB = -1.0
        else
            Bv = kb * log1p(exp(-e)); dB = -1.0 / (1.0 + exp(e))
        end
        obj += opt.barrier_w * Bv
        coef = opt.barrier_w * dB / ob.area0[f] * 0.5
        for i in 1:nf
            u0 = vs[i]
            nx = vs[mod1(i + 1, nf)]; pv = vs[mod1(i - 1, nf)]
            gX[u0, 1] += coef * (X[nx, 2] - X[pv, 2])
            gX[u0, 2] += coef * (X[pv, 1] - X[nx, 1])
        end
    end

    # ---- pull back to T --------------------------------------------------------
    for i in 1:np
        v = c.prime_to_original[i]
        gX[v, 1] += gC[i, 1]
        gX[v, 2] += gC[i, 2]
    end
    if k > 0
        gT = ob.Phi' * gX          # k x 2
        gU = ob.Uphi' * gS         # k x 2
        for j in 1:k
            g[2j - 1] = gT[j, 1] + 2.0 * gU[j, 2]
            g[2j] = gT[j, 2] - 2.0 * gU[j, 1]
        end
    end
    return obj
end

"""
    maximize_range(c, X0, Phi, o = RangeOptOptions()) -> RangeOptResult

Softmin range maximisation over X = X0 + Phi T with active-set refreshes; the returned
T is the best REFEREED (collision.jl bisection) iterate.
"""
function maximize_range(c::CutStructure, X0::Vector{Vec2}, Phi::AbstractMatrix{Float64},
                        o::RangeOptOptions = RangeOptOptions())
    k = size(Phi, 2)
    res = RangeOptResult(zeros(max(k, 0), 2), copy(X0), 0.0, 0.0, 0.0, 0.0, 0, 0, Float64[])
    res.theta_ref_before = theta_max(c, X0, 90, 40).theta_max_geometric
    res.theta_ref_after = res.theta_ref_before
    k == 0 && return res

    ob = RangeObjective()
    setup!(ob, c, X0, Phi, o)
    t = zeros(2 * k)
    best_t = copy(t)
    rebuild_active!(ob, t)
    res.n_active = length(ob.active)
    let B = basis(ob, t), sd = swept_discs(c, B)
        res.theta_closed_before = exact_theta_max(c, B, candidate_pairs(c, sd, o.theta_hi, true)).theta_max
    end
    best_ref = res.theta_ref_before

    fg = (x, g) -> value_and_grad(ob, x, g)
    for round in 1:o.rounds
        rebuild_active!(ob, t)
        isempty(ob.active) && break
        lo = LbfgsOptions()
        lo.max_iter = o.iters_per_round
        lr = lbfgs_minimize(fg, t, lo)
        t = lr.x
        if o.step_cap > 0
            mx = maximum(abs, ob.Phi * _unflatten_T(t, k))
            mx > o.step_cap && (t *= o.step_cap / mx)
        end
        Xc = matrix_to_points(ob.X0m + ob.Phi * _unflatten_T(t, k))
        ref = theta_max(c, Xc, 90, 40).theta_max_geometric
        push!(res.trace, ref)
        if ref > best_ref + 1e-12
            best_ref = ref
            best_t = copy(t)
        end
        res.rounds_run += 1
    end
    T = _unflatten_T(best_t, k)
    res.T = T
    res.X_opt = matrix_to_points(ob.X0m + ob.Phi * T)
    res.theta_ref_after = best_ref
    let B = basis(ob, best_t), sd = swept_discs(c, B)
        res.theta_closed_after = exact_theta_max(c, B, candidate_pairs(c, sd, o.theta_hi, true)).theta_max
    end
    return res
end
