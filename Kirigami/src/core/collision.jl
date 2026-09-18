# core/collision.jl -- polygon overlap, the deployment range, and the collision-aware
# optimization of Eq. (7)-(9), 2026 Sec. 4.5 / Sec. 5.2.
#
# The orientation determinants are accumulated in plain Float64 (the reference platform,
# Apple arm64, has no wider hardware type), and every reduction below is written in a
# fixed contraction form so the predicate's tie decisions are platform-independent.

rot90(v::Vec2) = Vec2(-v[2], v[1])

# Dot products and squared norms are UNFUSED sums of products, whereas StaticArrays' `dot`
# is a muladd that becomes an fma on aarch64.  Inline cross products `ux*vy - uy*vx` and
# `ux*ux + uy*uy` ARE contracted (fma on the first product).  Both forms are written out
# explicitly so the predicate's tie decisions are the same on every platform
# (docs/NUMERICS.md).
_dot(a::Vec2, b::Vec2) = a[1] * b[1] + a[2] * b[2]
_sq(a::Vec2) = a[1] * a[1] + a[2] * a[2]

# ---------------------------------------------------------------------------
# Robust primitives for the overlap predicate (F34 fix).
#
# The predicate decides whether the OPEN interiors of two simple polygons meet.
# Touching at a shared vertex (the hinge pin of two faces incident to the same
# hinge edge) or along a shared edge is NOT an overlap.  The two copies of a
# shared pin are bit-for-bit equal (both read the same Y[pv]), so every
# orientation determinant that involves the pin twice is exactly zero; the old
# `shrink_poly` heuristic destroyed that cancellation by moving each copy toward
# a different centroid, which is the whole of F34.  Nothing is shrunk here.
# ---------------------------------------------------------------------------

# Orientation determinant of (b - a, c - a).
function orient_raw(a::Vec2, b::Vec2, c::Vec2)
    ux = b[1] - a[1]
    uy = b[2] - a[2]
    vx = c[1] - a[1]
    vy = c[2] - a[2]
    return fma(ux, vy, -(uy * vx))
end

# Sign of the orientation, with a RELATIVE tolerance: the point c counts as ON the
# line ab when its distance from that line is at most `eps` times the scale of the
# configuration, i.e. |det| <= eps * |u| * S with S = max(|v|, |a|, |b|, |c|).
#
# The |a|,|b|,|c| terms matter: the points fed to this predicate include midpoints
# a + t (b - a), whose own rounding error is an absolute ulp of the WORLD coordinate
# (~4e-15 on a sheet 30 units across), not of the edge length.  A purely
# |u|*|v|-relative threshold rejects those points from their own edge, which produced
# a false overlap between two faces meeting along a split edge whose two copies
# coincide exactly at theta = 0 (K5 voronoi_75, faces 106/390).  The predicate stays
# exact (eps plays no role) whenever two of the three points coincide.
function orient_sign(a::Vec2, b::Vec2, c::Vec2, eps::Float64)
    d = orient_raw(a, b, c)
    d == 0 && return 0
    ux = b[1] - a[1]; uy = b[2] - a[2]
    vx = c[1] - a[1]; vy = c[2] - a[2]
    s2 = fma(vx, vx, vy * vy)
    s2 = max(s2, _sq(a))
    s2 = max(s2, _sq(b))
    s2 = max(s2, _sq(c))
    # |d| <= eps * |u| * S, compared squared so no sqrt is needed.
    thr2 = eps * eps * fma(ux, ux, uy * uy) * s2
    d * d <= thr2 && return 0
    return d > 0 ? 1 : -1
end

# Is p on the CLOSED segment [a,b]?  (T4.1b's convention: closed segments.)
function on_segment(p::Vec2, a::Vec2, b::Vec2, eps::Float64)
    orient_sign(a, b, p, eps) != 0 && return false
    e = b - a
    l2 = _sq(e)
    l2 == 0.0 && return _sq(p - a) == 0.0
    s = _dot(p - a, e)
    return s >= -eps * l2 && s <= l2 * (1.0 + eps)
end

# +1 strictly inside, 0 on the boundary, -1 strictly outside.
function point_in_poly3(p::Vec2, P::Vector{Vec2}, eps::Float64)
    n = length(P)
    j = n
    for i in 1:n
        on_segment(p, P[j], P[i], eps) && return 0
        j = i
    end
    inside = false
    j = n
    for i in 1:n
        if ((P[i][2] > p[2]) != (P[j][2] > p[2])) &&
           (p[1] < (P[j][1] - P[i][1]) * (p[2] - P[i][2]) / (P[j][2] - P[i][2]) + P[i][1])
            inside = !inside
        end
        j = i
    end
    return inside ? 1 : -1
end

# Parameters t in [0,1] along [a,b] at which [a,b] meets the CLOSED segment [c,d].
# Collinear overlaps contribute the projections of c and d.  Signs are taken from
# orient_sign (tolerant); the parameter itself from the raw determinants.
function segment_hits!(ts::Vector{Float64}, a::Vec2, b::Vec2, c::Vec2, d::Vec2, eps::Float64)
    s1 = orient_sign(a, b, c, eps); s2 = orient_sign(a, b, d, eps)
    s3 = orient_sign(c, d, a, eps); s4 = orient_sign(c, d, b, eps)
    e = b - a
    l2 = _sq(e)
    l2 == 0.0 && return
    push_t(t) = (t >= 0.0 && t <= 1.0) && push!(ts, t)
    if s1 == 0 && s2 == 0  # collinear
        push_t(_dot(c - a, e) / l2)
        push_t(_dot(d - a, e) / l2)
        return
    end
    (s1 * s2 > 0 || s3 * s4 > 0) && return  # one segment entirely on one side
    s3 == 0 && push_t(0.0)  # a lies on [c,d]
    s4 == 0 && push_t(1.0)  # b lies on [c,d]
    (s3 == 0 || s4 == 0) && return
    d3 = orient_raw(c, d, a); d4 = orient_raw(c, d, b)
    den = d3 - d4
    den == 0 && return
    push_t(d3 / den)
    return
end

# A point strictly inside a simple polygon (ear construction; no shrinking).
# Returns the point, or `nothing`.
function interior_point(P::Vector{Vec2}, eps::Float64)
    n = length(P)
    n < 3 && return nothing
    area2 = 0.0
    for i in 1:n
        u = P[i]
        v = P[mod1(i + 1, n)]
        area2 += fma(u[1], v[2], -(u[2] * v[1]))
    end
    orient = area2 >= 0 ? 1 : -1
    for i in 1:n
        ip = mod1(i - 1, n); inx = mod1(i + 1, n)
        u = P[ip]; v = P[i]; w = P[inx]
        orient_sign(u, v, w, eps) != orient && continue  # not a strictly convex corner
        # the vertex inside triangle (u,v,w) that is farthest from the line uw, if any
        best = 0
        bestd = 0.0
        for j in 1:n
            (j == i || j == ip || j == inx) && continue
            q = P[j]
            if orient_sign(u, v, q, eps) != orient || orient_sign(v, w, q, eps) != orient ||
               orient_sign(w, u, q, eps) != orient
                continue
            end
            dd = abs(orient_raw(u, w, q))
            if dd > bestd
                bestd = dd
                best = j
            end
        end
        cand = best == 0 ? Vec2((u + v + w) / 3.0) : Vec2(0.5 * (v + P[best]))
        point_in_poly3(cand, P, eps) == 1 && return cand
    end
    return nothing
end

function bbox(P::Vector{Vec2})
    lo = P[1]; hi = P[1]
    for p in P
        lo = min.(lo, p)
        hi = max.(hi, p)
    end
    return lo, hi
end

"""
    polygons_overlap(A, B, tol_rel = 1e-12) -> Bool

Overlap test for two simple polygons: true iff their OPEN interiors meet. Touching
at a shared vertex (the hinge pin two incident faces share) or along a shared edge
is NOT an overlap, and is decided exactly -- the shared pin has identical coordinates
in both faces, so the orientation determinants that involve it vanish exactly.

`tol_rel` is a RELATIVE degeneracy tolerance on the orientation predicate: a triple
with |det| <= tol_rel * |u| * |v| (i.e. |sin angle| <= tol_rel) counts as collinear.
Values are clamped to [1e-14, 1e-9], so the historical call sites (1e-12, 1e-9, 1e-6)
all agree; the floor because the midpoints the routine constructs carry an ulp of the
world coordinate.
"""
function polygons_overlap(A::Vector{Vec2}, B::Vector{Vec2}, tol_rel::Float64 = 1e-12)
    (length(A) < 3 || length(B) < 3) && return false
    eps = min(max(tol_rel, 1e-14), 1e-9)
    la, ha = bbox(A)
    lb, hb = bbox(B)
    (ha[1] < lb[1] || hb[1] < la[1] || ha[2] < lb[2] || hb[2] < la[2]) && return false
    # Interiors meet iff some sub-arc of one boundary runs strictly inside the other, or
    # the two polygons have identical boundaries.  Split each edge at every crossing with
    # the other boundary and classify the midpoint of each piece.
    ts = Float64[]
    function sweep(P::Vector{Vec2}, Q::Vector{Vec2})
        np = length(P); nq = length(Q)
        for i in 1:np
            a = P[i]; b = P[mod1(i + 1, np)]
            empty!(ts)
            push!(ts, 0.0)
            push!(ts, 1.0)
            for j in 1:nq
                segment_hits!(ts, a, b, Q[j], Q[mod1(j + 1, nq)], eps)
            end
            sort!(ts)
            for k in 1:length(ts)-1
                (ts[k+1] > ts[k]) || continue
                mid = a + 0.5 * (ts[k] + ts[k+1]) * (b - a)
                point_in_poly3(mid, Q, eps) == 1 && return true
            end
        end
        return false
    end
    sweep(A, B) && return true
    sweep(B, A) && return true

    # Boundaries never enter the other's interior: the remaining overlap is containment
    # with coincident boundaries (A == B, or one inside the other touching everywhere).
    pa = interior_point(A, eps)
    pa !== nothing && point_in_poly3(pa, B, eps) == 1 && return true
    pb = interior_point(B, eps)
    pb !== nothing && point_in_poly3(pb, A, eps) == 1 && return true
    return false
end

"""True iff any two faces of the deployed structure overlap."""
function has_collision(c::CutStructure, Y::Vector{Vec2}, tol_rel::Float64 = 1e-12)
    F = n_faces(c.mesh)
    polys = [Vec2[Y[pv] for pv in c.prime_faces[f]] for f in 1:F]
    lo = Vector{Vec2}(undef, F)
    hi = Vector{Vec2}(undef, F)
    for f in 1:F
        lo[f], hi[f] = bbox(polys[f])
    end
    for i in 1:F, j in i+1:F
        (hi[i][1] < lo[j][1] || hi[j][1] < lo[i][1] || hi[i][2] < lo[j][2] ||
         hi[j][2] < lo[i][2]) && continue
        polygons_overlap(polys[i], polys[j], tol_rel) && return true
    end
    return false
end

mutable struct ThetaMaxReport
    theta_max_geometric::Float64  # first overlap, by grid scan + bisection
    min_beta::Float64             # min_i (2*pi - alpha_lik), 2026 Sec. 5.2
    argmin_beta_edge::Int         # edge index of M (0 = no hinge edges)
    beta::Vector{Float64}         # one entry per hinge edge, same order as cut.hinge_edges
    collided::Bool                # whether an overlap was found at all in (0, pi]
end
ThetaMaxReport() = ThetaMaxReport(0.0, 0.0, 0, Float64[], false)

"""
    hinge_beta(c, X) -> Vector{Float64}

Analytic per-hinge bound beta = 2*pi - alpha_1 - alpha_2 where alpha_1, alpha_2 are
the interior angles of the two incident faces at the hinge (source) vertex.
"""
function hinge_beta(c::CutStructure, X::Vector{Vec2})
    m = c.mesh
    tmp = Mesh(X, m.faces)  # face_angle_at only reads X and faces
    beta = Float64[]
    for e in c.hinge_edges
        v = c.hinge_dir[e].src
        ed = m.edges[e]
        f1 = m.half_edges[ed.he[1]].face
        f2 = m.half_edges[ed.he[2]].face
        a1 = face_angle_at(tmp, f1, v)
        a2 = face_angle_at(tmp, f2, v)
        push!(beta, 2pi - a1 - a2)
    end
    return beta
end

"""
    theta_max(c, X, grid = 90, bisect_iters = 40, theta_hi = pi) -> ThetaMaxReport

Grid scan over (0, theta_hi] for the first colliding deployment, refined by bisection.
"""
function theta_max(c::CutStructure, X::Vector{Vec2}, grid::Int = 90, bisect_iters::Int = 40,
                   theta_hi::Float64 = Float64(pi))
    r = ThetaMaxReport()
    r.beta = hinge_beta(c, X)
    if !isempty(r.beta)
        r.min_beta, imin = findmin(r.beta)
        r.argmin_beta_edge = c.hinge_edges[imin]
    end

    collides(th) = has_collision(c, deploy(c, X, th).Y)

    lo = 0.0
    hi = -1.0
    for i in 1:grid
        th = theta_hi * i / grid
        if collides(th)
            hi = th
            break
        end
        lo = th
    end
    if hi < 0
        r.theta_max_geometric = theta_hi
        r.collided = false
        return r
    end
    r.collided = true
    for _ in 1:bisect_iters
        mid = 0.5 * (lo + hi)
        if collides(mid)
            hi = mid
        else
            lo = mid
        end
    end
    r.theta_max_geometric = lo
    return r
end

"""
    velocity_operator_apply(c, cols) -> Matrix

Applies the (linear) velocity operator column-wise: input N x C, output n_prime x C,
returning u with z = rot90(u). Used to make Eq. (9) affine in the shape-space
coefficients.
"""
function velocity_operator_apply(c::CutStructure, cols::AbstractMatrix{Float64})
    m = c.mesh
    F = n_faces(m)
    C = size(cols, 2)
    cf = [-0.5 * m.sigma[f] for f in 1:F]

    adj = _hinge_adjacency(c)
    W = zeros(F, C)
    seen = falses(F)
    for s in 1:F
        seen[s] && continue
        seen[s] = true
        W[s, :] .= 0.0
        q = Int[s]
        while !isempty(q)
            f = popfirst!(q)
            for (e, g) in adj[f]
                seen[g] && continue
                seen[g] = true
                v = c.hinge_dir[e].src
                W[g, :] = W[f, :] + (cf[f] - cf[g]) * cols[v, :]
                push!(q, g)
            end
        end
    end
    U = zeros(c.n_prime_vertices, C)
    for f in 1:F
        vs = m.faces[f]
        for k in eachindex(vs)
            U[c.prime_faces[f][k], :] = cf[f] * cols[vs[k], :] + W[f, :]
        end
    end
    return U
end

"""Rows z_i = d Y_theta / d theta at theta = 0 (Eq. 7), for the M' vertices."""
function deployment_velocity(c::CutStructure, X::Vector{Vec2})
    U = velocity_operator_apply(c, points_to_matrix(X))
    return [rot90(Vec2(U[i, 1], U[i, 2])) for i in 1:size(U, 1)]
end

Base.@kwdef mutable struct CollisionOptOptions
    gamma::Float64 = 1e-2   # regularizer weight (the paper gives no value)
    kappa::Float64 = 0.05   # barrier smoothing (the paper gives no barrier formula)
    max_iter::Int = 400
end

mutable struct CollisionOptResult
    T::Matrix{Float64}          # k x 2 shape-space coefficients
    X_opt::Vector{Vec2}         # resulting embedding
    f0::Float64                 # objective before / after
    f1::Float64
    theta_max_before::Float64
    theta_max_after::Float64
    iterations::Int
end
CollisionOptResult() = CollisionOptResult(zeros(0, 2), Vec2[], 0.0, 0.0, 0.0, 0.0, 0)

# B(s) = kappa * log(1 + exp(-s/kappa)), B'(s) = -1/(1+exp(s/kappa)); returns (B, dB)
function _barrier(s::Float64, kappa::Float64)
    t = s / kappa
    if t > 30
        return kappa * exp(-t), -exp(-t)
    elseif t < -30
        return -s, -1.0
    else
        return kappa * log1p(exp(-t)), -1.0 / (1.0 + exp(t))
    end
end

"""
    optimize_collision(c, X0, Phi, opt = CollisionOptOptions()) -> CollisionOptResult

Eq. (9): minimize sum_{split edges} B(z_i . n_hat) + gamma ||Y - Y0||^2 over the
shape space Y = X0 + Phi T. Barrier B(s) = kappa * log(1 + exp(-s/kappa)).
"""
function optimize_collision(c::CutStructure, X0::Vector{Vec2}, Phi::AbstractMatrix{Float64},
                            opt::CollisionOptOptions = CollisionOptOptions())
    m = c.mesh
    N = n_vertices(m)
    k = size(Phi, 2)
    res = CollisionOptResult()
    res.T = zeros(max(k, 0), 2)
    res.X_opt = X0

    # Per split edge and per incident face: the M'-vertex pair (pa, pb) in the face's
    # CCW-stored order, plus the corresponding original vertices (va, vb).
    items = NTuple{4,Int}[]
    for e in c.split_edges
        ed = m.edges[e]
        for s in 1:2
            he = m.half_edges[ed.he[s]]
            f = he.face
            nf = length(m.faces[f])
            kk = he.corner
            push!(items, (c.prime_faces[f][kk], c.prime_faces[f][mod1(kk + 1, nf)], he.from, he.to))
        end
    end

    # Everything is affine in T: u = U0 + Urho * T, x = X0 + Phi * T.
    cols = zeros(N, 2 + k)
    cols[:, 1:2] = points_to_matrix(X0)
    k > 0 && (cols[:, 3:end] = Phi)
    U = velocity_operator_apply(c, cols)
    U0 = U[:, 1:2]           # n_prime x 2
    Urho = U[:, 3:end]       # n_prime x k

    Xref = points_to_matrix(X0)
    kappa = opt.kappa
    gamma = opt.gamma

    function eval(tvec::Vector{Float64}, grad::Vector{Float64})
        T = zeros(k, 2)
        for i in 1:k
            T[i, 1] = tvec[2i-1]
            T[i, 2] = tvec[2i]
        end
        Xm = Xref + (k > 0 ? Phi * T : zeros(N, 2))
        Um = U0 + (k > 0 ? Urho * T : zeros(size(U0, 1), 2))
        grad .= 0.0
        f = 0.0
        for (pa, pb, va, vb) in items
            d = Vec2(Xm[vb, 1] - Xm[va, 1], Xm[vb, 2] - Xm[va, 2])
            dn = norm(d)
            dn < 1e-12 && continue
            for side in 1:2
                pv = side == 1 ? pa : pb
                u = Vec2(Um[pv, 1], Um[pv, 2])
                # z . n_hat with z = rot90(u), n_hat = rot90(d)/|d|  ==>  (u . d)/|d|
                s = dot(u, d) / dn
                B, dB = _barrier(s, kappa)
                f += B
                for mIdx in 1:k
                    rho = Urho[pv, mIdx]
                    del = Phi[vb, mIdx] - Phi[va, mIdx]
                    # ds/dT_m = (rho*d + del*u)/|d| - s*del*d/|d|^2
                    ds = (rho * d + del * u) / dn - s * del * d / (dn * dn)
                    grad[2mIdx-1] += dB * ds[1]
                    grad[2mIdx] += dB * ds[2]
                end
            end
        end
        diff = Xm - Xref
        f += gamma * sum(abs2, diff)
        if k > 0
            gp = 2.0 * gamma * (Phi' * diff)  # k x 2
            for mIdx in 1:k
                grad[2mIdx-1] += gp[mIdx, 1]
                grad[2mIdx] += gp[mIdx, 2]
            end
        end
        return f
    end

    res.theta_max_before = theta_max(c, X0).theta_max_geometric
    if k == 0 || isempty(items)
        res.f0 = res.f1 = 0.0
        res.theta_max_after = res.theta_max_before
        return res
    end
    t0 = zeros(2k)
    g = zeros(2k)
    res.f0 = eval(t0, g)
    lo = LbfgsOptions(max_iter = opt.max_iter)
    lr = lbfgs_minimize(eval, t0, lo)
    res.iterations = lr.iterations
    res.f1 = lr.f
    for i in 1:k
        res.T[i, 1] = lr.x[2i-1]
        res.T[i, 2] = lr.x[2i]
    end
    Xm = Xref + Phi * res.T
    res.X_opt = matrix_to_points(Xm)
    res.theta_max_after = theta_max(c, res.X_opt).theta_max_geometric
    return res
end

# CollisionSweepResult extends CollisionOptResult; the fields are flattened into one struct.
mutable struct CollisionSweepResult
    T::Matrix{Float64}
    X_opt::Vector{Vec2}
    f0::Float64
    f1::Float64
    theta_max_before::Float64
    theta_max_after::Float64
    iterations::Int
    gamma_used::Float64
    ladder::Vector{Tuple{Float64,Float64}}  # (gamma, theta_max_after)
end

const DEFAULT_COLLISION_GAMMAS = [1e-4, 1e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0, 3.0, 10.0]

"""
    optimize_collision_sweep(c, X0, Phi, gammas = DEFAULT_COLLISION_GAMMAS, kappa = 0.05)

The paper gives no value for gamma, and a single gamma can make the deployment
range WORSE (measured; see results/core_validation/reference_cases.md): the
objective of Eq. (9) contains nothing that keeps the flat embedding embedded.
This sweeps a ladder of gammas and keeps the run with the largest theta_max,
falling back to T = 0 (i.e. X0 unchanged), so the result is never worse than
the input. `gamma_used` reports the winner (Inf for the fallback).
"""
function optimize_collision_sweep(c::CutStructure, X0::Vector{Vec2}, Phi::AbstractMatrix{Float64},
                                  gammas::Vector{Float64} = DEFAULT_COLLISION_GAMMAS,
                                  kappa::Float64 = 0.05)
    best = CollisionSweepResult(zeros(max(size(Phi, 2), 0), 2), X0, 0.0, 0.0, 0.0, 0.0, 0, Inf,
                                Tuple{Float64,Float64}[])
    best.theta_max_before = theta_max(c, X0).theta_max_geometric
    best.theta_max_after = best.theta_max_before  # T = 0 fallback
    for g in gammas
        o = CollisionOptOptions(gamma = g, kappa = kappa)
        r = optimize_collision(c, X0, Phi, o)
        push!(best.ladder, (g, r.theta_max_after))
        if r.theta_max_after > best.theta_max_after + 1e-9
            best.T = r.T
            best.X_opt = r.X_opt
            best.f0 = r.f0
            best.f1 = r.f1
            best.iterations = r.iterations
            best.theta_max_after = r.theta_max_after
            best.gamma_used = g
        end
    end
    return best
end
