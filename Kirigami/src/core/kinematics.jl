# core/kinematics.jl -- forward kinematics of the kirigami structure M' (2026 Sec. 4.5).
#
# Convention (derived, see code/README.md): every face f is transformed rigidly by
#     y = R(-sigma(f) * theta / 2) * x + t_f,
# so that across a hinge edge the two duplicates open by exactly |theta|, and across
# a split edge both faces carry the same rotation (hence the duplicates of a split
# edge stay parallel -- Remark A.4, which this convention reproduces automatically).
# The translations t_f are propagated by BFS over the hinge adjacency of M'.
#
# Port of code/src/core/kinematics.{hpp,cpp}. 1-based faces; `seed_face` defaults to 1
# (the C++ default 0 = the first face).

const Mat2 = SMatrix{2,2,Float64,4}

mutable struct Deployment
    theta::Float64
    Y::Vector{Vec2}           # positions of the M' vertices
    dY_dtheta::Vector{Vec2}   # derivative of Y w.r.t. theta
    reached::Vector{Bool}     # per face: reached by the BFS
    max_mismatch::Float64     # largest disagreement over non-tree hinge adjacencies
    bfs_order::Vector{Int}    # faces in BFS order
end

# Trig via the project's libm shim (mesh.jl) so the rotation entries carry the same bits
# as the C++ (which linked Apple's libm) wherever the platform allows.
_rot(a::Float64) = Mat2(libm_cos(a), libm_sin(a), -libm_sin(a), libm_cos(a))     # column-major
_drot(a::Float64) = Mat2(-libm_sin(a), libm_cos(a), -libm_cos(a), -libm_sin(a))  # d/da rot(a)

# Eigen's fixed-size 2x2 * vector product as clang -O2 (-ffp-contract=on) compiled it:
# y_i = fma(R(i,1), x_1, R(i,0) * x_0). Verified bit-exact against the C++ deployed
# positions of four frozen patches (382/382 vertices); a plain `R * x` differs by an ulp.
_mv(R::Mat2, x::Vec2) = Vec2(fma(R[1, 2], x[2], R[1, 1] * x[1]), fma(R[2, 2], x[2], R[2, 1] * x[1]))

# face -> list of (hinge edge, neighbour face)
function _hinge_adjacency(c::CutStructure)
    m = c.mesh
    adj = [Tuple{Int,Int}[] for _ in 1:n_faces(m)]
    for e in c.hinge_edges
        ed = m.edges[e]
        f1 = m.half_edges[ed.he[1]].face
        f2 = m.half_edges[ed.he[2]].face
        push!(adj[f1], (e, f2))
        push!(adj[f2], (e, f1))
    end
    return adj
end

function _run_kinematics(c::CutStructure, X::Vector{Vec2}, theta::Float64,
                         order_hint::Vector{Int}, seed_face::Int)
    m = c.mesh
    F = n_faces(m)
    d = Deployment(theta, fill(Vec2(0, 0), c.n_prime_vertices),
                   fill(Vec2(0, 0), c.n_prime_vertices), falses(F), 0.0, Int[])

    R = Vector{Mat2}(undef, F)
    dR = Vector{Mat2}(undef, F)
    for f in 1:F
        a = -m.sigma[f] * theta * 0.5
        R[f] = _rot(a)
        dR[f] = _drot(a) * (-m.sigma[f] * 0.5)
    end
    t = fill(Vec2(0, 0), F)
    dt = fill(Vec2(0, 0), F)

    adj = _hinge_adjacency(c)

    if !isempty(order_hint)
        seeds = order_hint
    else
        seeds = collect(1:F)
        if seed_face > 1 && seed_face <= F
            seeds[1], seeds[seed_face] = seeds[seed_face], seeds[1]
        end
    end

    queued = falses(F)
    for s in seeds
        queued[s] && continue
        queued[s] = true
        t[s] = Vec2(0, 0)
        dt[s] = Vec2(0, 0)
        q = Int[s]
        while !isempty(q)
            f = popfirst!(q)
            d.reached[f] = true
            push!(d.bfs_order, f)
            for (e, g) in adj[f]
                v = c.hinge_dir[e].src  # the hinge stays at the source vertex
                xv = X[v]
                p = _mv(R[f], xv) + t[f]
                dp = _mv(dR[f], xv) + dt[f]
                if !queued[g]
                    queued[g] = true
                    t[g] = p - _mv(R[g], xv)
                    dt[g] = dp - _mv(dR[g], xv)
                    push!(q, g)
                else
                    p2 = _mv(R[g], xv) + t[g]
                    d.max_mismatch = max(d.max_mismatch, norm(p - p2))
                end
            end
        end
    end

    for f in 1:F
        vs = m.faces[f]
        for k in eachindex(vs)
            pv = c.prime_faces[f][k]
            d.Y[pv] = _mv(R[f], X[vs[k]]) + t[f]
            d.dY_dtheta[pv] = _mv(dR[f], X[vs[k]]) + dt[f]
        end
    end
    return d
end

"""
    deploy(c, X, theta, seed_face = 1) -> Deployment

Deploys M' at angle theta from the flat embedding X of M, seeding at `seed_face`.
If `theta` is 0 the result is the flat structure (all duplicates coincident).
"""
deploy(c::CutStructure, X::Vector{Vec2}, theta::Real, seed_face::Int = 1) =
    _run_kinematics(c, X, Float64(theta), Int[], seed_face)

"""Deploys with a caller-supplied face visit order (used to check order independence)."""
deploy_with_order(c::CutStructure, X::Vector{Vec2}, theta::Real, face_order::Vector{Int}) =
    _run_kinematics(c, X, Float64(theta), face_order, 1)

"""Best rigid alignment (rotation + translation) of `a` onto `b`; returns max residual."""
function rigid_align_residual(a::Vector{Vec2}, b::Vector{Vec2})
    (length(a) != length(b) || isempty(a)) && return Inf
    ca = sum(a) / length(a)
    cb = sum(b) / length(b)
    H = zeros(Mat2)
    for i in eachindex(a)
        H += (a[i] - ca) * (b[i] - cb)'
    end
    S = svd(H)
    R = S.V * S.U'
    if det(R) < 0
        V = Matrix(S.V)
        V[:, 2] .*= -1
        R = Mat2(V * S.U')
    end
    worst = 0.0
    for i in eachindex(a)
        worst = max(worst, norm(R * (a[i] - ca) - (b[i] - cb)))
    end
    return worst
end
