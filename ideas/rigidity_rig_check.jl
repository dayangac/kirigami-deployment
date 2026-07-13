# Standalone check of two claims made by the rigidity ideator.
# Port of rigidity_rig_check.cpp (same graphs, same seeds, same printed lines).
#
# Run (from repo root):
#   julia --project=Kirigami ideas/rigidity_rig_check.jl
#
#  V1: Y_theta = cos(theta/2)*C + sin(theta/2)*S exactly (forward kinematics is
#      linear in (cos(theta/2), sin(theta/2))).
#  V2: mobility m = dim ker A - 1, where A(omega)_z = sum_e z_e (omega_{f1}-omega_{f2}) p_e
#      over a cycle basis {z} of the hinge graph Gamma, compared with
#      m = 3|F| - rank(RigidityMatrix) - 3.
#  V3: sigma in ker A  <=>  Eq.(2) holds (uniform deployability).
#  V4: dim ker A(theta) along the deployment path.

using Kirigami
using LinearAlgebra
using Printf
using StaticArrays
const K = Kirigami
const Vec2 = K.Vec2

const Jrot = SMatrix{2,2,Float64}(0, 1, -1, 0)   # [0 -1; 1 0]

mutable struct Gamma
    F::Int
    h::Vector{Int}                 # per hinge edge: head face, tail face
    t::Vector{Int}
    eid::Vector{Int}               # mesh edge id
    cyc::Vector{Vector{Float64}}   # cycle basis vectors in R^{|E_h|}
    comps::Int
end

function build_gamma(c::K.CutStructure)
    m = c.mesh
    g = Gamma(K.n_faces(m), Int[], Int[], Int[], Vector{Float64}[], 0)
    for e in c.hinge_edges
        ed = m.edges[e]
        push!(g.h, m.half_edges[ed.he[1]].face)
        push!(g.t, m.half_edges[ed.he[2]].face)
        push!(g.eid, e)
    end
    E = length(g.h)
    adj = [Tuple{Int,Int}[] for _ in 1:g.F]   # (edge, other face)
    for i in 1:E
        push!(adj[g.h[i]], (i, g.t[i])); push!(adj[g.t[i]], (i, g.h[i]))
    end
    parent = fill(-2, g.F); pedge = fill(-1, g.F)
    gpath = [zeros(E) for _ in 1:g.F]
    intree = falses(E)
    for s in 1:g.F
        parent[s] != -2 && continue
        g.comps += 1
        parent[s] = -1; stack = [s]
        while !isempty(stack)
            f = pop!(stack)
            for (i, o) in adj[f]
                parent[o] != -2 && continue
                parent[o] = f; pedge[o] = i; intree[i] = true
                gpath[o] = copy(gpath[f])
                gpath[o][i] += (g.h[i] == o ? 1.0 : -1.0)   # traversing f -> o
                push!(stack, o)
            end
        end
    end
    for i in 1:E
        intree[i] && continue
        z = copy(gpath[g.t[i]])
        z[i] += 1.0
        z .-= gpath[g.h[i]]
        push!(g.cyc, z)
    end
    return g
end

# A: (2 * n_cycles) x F
function build_A(g::Gamma, p::Vector{Vec2})
    nz = length(g.cyc)
    A = zeros(2 * nz, g.F)
    for k in 1:nz, i in 1:length(g.h)
        z = g.cyc[k][i]
        z == 0 && continue
        A[2k-1, g.h[i]] += z * p[i][1]; A[2k-1, g.t[i]] -= z * p[i][1]
        A[2k,   g.h[i]] += z * p[i][2]; A[2k,   g.t[i]] -= z * p[i][2]
    end
    return A
end

# full body-pin rigidity matrix: 2|E_h| x 3|F|, unknowns (w_f, vx_f, vy_f)
function build_R(g::Gamma, p::Vector{Vec2})
    E = length(g.h)
    R = zeros(2 * E, 3 * g.F)
    for i in 1:E
        Jp = Jrot * p[i]
        hf = 3 * (g.h[i] - 1); tf = 3 * (g.t[i] - 1)
        R[2i-1, hf+1] += Jp[1]; R[2i-1, tf+1] -= Jp[1]
        R[2i,   hf+1] += Jp[2]; R[2i,   tf+1] -= Jp[2]
        R[2i-1, hf+2] += 1;     R[2i-1, tf+2] -= 1
        R[2i,   hf+3] += 1;     R[2i,   tf+3] -= 1
    end
    return R
end

function rank_of(M::Matrix{Float64}, tol::Float64 = 1e-9)
    isempty(M) && return 0
    s = svdvals(M)
    t = tol * (isempty(s) ? 1.0 : s[1])
    return count(x -> x > t, s)
end

function run(name::String, m::K.Mesh, rng::K.MT19937)
    try
        run1(name, m, rng)
    catch e
        @printf("%-20s EXCEPTION: %s\n", name, sprint(showerror, e))
    end
end

function run1(name::String, m::K.Mesh, rng::K.MT19937)
    K.build_topology!(m)
    if isempty(m.sigma)
        orep = K.n_faces(m) <= 18 ? K.brute_force_orientation(m) :
                                    K.assign_orientation_relaxation(m, rng)
        m.sigma = orep.sigma
    end
    if isempty(m.sigma) || length(m.sigma) != K.n_faces(m)
        @printf("%-20s sigma assignment failed, skipped\n", name); return
    end
    K.build_topology!(m)
    c = K.make_cut(m)
    g = build_gamma(c)
    if isempty(g.h)
        @printf("%-22s no hinge edges, skipped\n", name); return
    end

    # hinge points at theta = 0 are just X[src]
    p0 = [m.X[c.hinge_dir[g.eid[i]].src] for i in 1:length(g.h)]

    # ---- V1: linearity of forward kinematics in (cos(th/2), sin(th/2)) ----
    Yof(th) = K.deploy(c, m.X, th, 1).Y
    ta = 0.37; tb = 1.13
    Ya = Yof(ta); Yb = Yof(tb)
    ca = cos(ta / 2); sa = sin(ta / 2); cb = cos(tb / 2); sb = sin(tb / 2)
    dt = ca * sb - sa * cb
    NP = c.n_prime_vertices
    C = [(sb * Ya[i] - sa * Yb[i]) / dt for i in 1:NP]
    S = [(-cb * Ya[i] + ca * Yb[i]) / dt for i in 1:NP]
    v1err = 0.0; scale = 0.0
    for th in (0.05, 0.7, 1.9, 2.5, 3.0)
        Y = Yof(th)
        cc = cos(th / 2); ss = sin(th / 2)
        for i in 1:NP
            v1err = max(v1err, norm(Y[i] - (cc * C[i] + ss * S[i])))
            scale = max(scale, norm(Y[i]))
        end
    end

    # ---- V2: mobility ----
    A0 = build_A(g, p0)
    R0 = build_R(g, p0)
    dimW = g.F - rank_of(A0)
    m_A = dimW - g.comps                        # one trivial global rotation per component
    m_R = 3 * g.F - rank_of(R0) - 3 * g.comps   # Maxwell/Calladine, trivial motions removed

    # ---- V3: sigma in ker A  <=>  Eq.(2) ----
    sg = Float64.(m.sigma)
    sigres = norm(A0 * sg)
    hs = K.holes_partition(c)
    rr = K.hole_residuals(c, m.X, hs)

    # ---- V3b: repeat on a SOLVED (deployable) embedding X0 ----
    sys = K.assemble_system(c, hs, m.X, K.Fixed)
    sr = K.solve_system(sys, m.X)
    X0 = K.matrix_to_points(sr.X0)
    p0s = [X0[c.hinge_dir[g.eid[i]].src] for i in 1:length(g.h)]
    As = build_A(g, p0s)
    dimWs = g.F - rank_of(As)
    mAs = dimWs - g.comps
    mRs = 3 * g.F - rank_of(build_R(g, p0s)) - 3 * g.comps
    sigres_s = norm(As * sg)
    rrs = K.hole_residuals(c, X0, hs)

    # ---- V4: dim ker A along the path ----
    nchange = 0; dim0 = -1; dims = ""
    for th in (0.0, 0.2, 0.6, 1.0, 1.6, 2.2, 2.8)
        Y = K.deploy(c, X0, th, 1).Y
        pth = [Y[K.prime_vertex(c, g.h[i], c.hinge_dir[g.eid[i]].src)] for i in 1:length(g.h)]
        d = g.F - rank_of(build_A(g, pth))
        dims *= @sprintf("%d ", d - g.comps)
        if dim0 < 0
            dim0 = d
        elseif d != dim0
            nchange += 1
        end
    end

    @printf("%-20s |F|=%4d |Eh|=%4d |Es|=%3d H=%4d c=%d | V1=%.2e(sc %.1f) | ini: m_A=%3d m_R=%3d %s |As|=%.1e res=%.1e | sol: m_A=%3d m_R=%3d %s |As|=%.1e res=%.1e | rkL=%d/%d m(th)= %s\n",
            name, g.F, length(g.h), K.n_split(c), K.n_interior_holes(hs), g.comps,
            v1err, scale, m_A, m_R, (m_A == m_R ? "OK" : "BAD"), sigres, rr.max_norm,
            mAs, mRs, (mAs == mRs ? "OK" : "BAD"), sigres_s, rrs.max_norm,
            sr.rank_L, sys.n_hole_rows, dims)
end

function main()
    rng = K.MT19937(12345)
    run("squares 3x3", K.tiling_squares(K.rect(Vec2(0, 0), 1.6, 1.6)), rng)
    run("squares 5x5", K.tiling_squares(K.rect(Vec2(0, 0), 2.6, 2.6)), rng)
    run("triangles", K.tiling_triangles(K.disk(Vec2(0, 0), 2.0)), rng)
    run("hexagons", K.tiling_hexagons(K.disk(Vec2(0, 0), 2.5)), rng)
    run("kagome", K.tiling_kagome(K.disk(Vec2(0, 0), 2.5)), rng)
    run("periodic squares", K.periodic_squares(3, 3), rng)
    for k in 0:7
        r2 = K.MT19937(100 + k)
        run("delaunay random", K.delaunay_of_random_points(30, 3.0, r2), rng)
    end
    for k in 0:7
        r2 = K.MT19937(200 + k)
        run("voronoi random", K.voronoi_of_random_points(25, 3.0, r2), rng)
    end
end

main()
