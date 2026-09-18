# method/mobility.jl -- the hinge graph Gamma = (F, E_hinge), its cycle space, the
# angular-velocity matrix A (2 b_1(Gamma) x |F|), its 2-core, and the mobility
# identity of kill experiment K3a:
#
#     dim ker A  ==  |F \ core2(Gamma)|  +  dim ker A|core2 .
#
# A is the closure operator of the rigidity check of ideas/ (K3a): for each cycle z of a
# spanning-forest cycle basis of Gamma,
#     sum_i z_i (omega_{head(i)} - omega_{tail(i)}) p_i = 0     in R^2,
# with p_i the pin position (the hinge SOURCE vertex) in the current configuration.
# Here the cycle basis is stored implicitly as a BFS forest plus the non-tree edges,
# so A is assembled sparsely and the whole thing scales to a few thousand faces.
#
# Faces, hinge edges and mesh edge ids are 1-based; a root's parent is 0. `subgraph`
# returns (HingeGraph, sub_pins).

mutable struct HingeGraph
    F::Int
    head::Vector{Int}   # per hinge edge: the two incident faces
    tail::Vector{Int}
    eid::Vector{Int}    # mesh edge index
    components::Int     # c(Gamma), isolated faces counted
    # BFS spanning forest
    parent::Vector{Int}   # 0 at a root
    pedge::Vector{Int}    # hinge edge to the parent, 0 at a root
    depth::Vector{Int}
    nontree::Vector{Int}  # non-tree hinge edge indices
end
HingeGraph() = HingeGraph(0, Int[], Int[], Int[], 0, Int[], Int[], Int[], Int[])
n_edges(g::HingeGraph) = length(g.head)
n_cycles(g::HingeGraph) = length(g.nontree)

function build_forest!(g::HingeGraph)
    E = n_edges(g)
    adj = [Tuple{Int,Int}[] for _ in 1:g.F]
    for i in 1:E
        push!(adj[g.head[i]], (i, g.tail[i]))
        push!(adj[g.tail[i]], (i, g.head[i]))
    end
    g.parent = zeros(Int, g.F)
    g.pedge = zeros(Int, g.F)
    g.depth = zeros(Int, g.F)
    seen = falses(g.F)
    intree = falses(E)
    g.components = 0
    for s in 1:g.F
        seen[s] && continue
        g.components += 1
        seen[s] = true
        q = Int[s]
        while !isempty(q)
            f = popfirst!(q)
            for (i, o) in adj[f]
                seen[o] && continue
                seen[o] = true
                g.parent[o] = f
                g.pedge[o] = i
                g.depth[o] = g.depth[f] + 1
                intree[i] = true
                push!(q, o)
            end
        end
    end
    g.nontree = Int[i for i in 1:E if !intree[i]]
    return g
end

# sign of the tree edge to f as traversed parent(f) -> f
tree_sign(g::HingeGraph, f::Int) = g.head[g.pedge[f]] == f ? 1.0 : -1.0

# Fundamental cycle of non-tree edge i, as (edge, coefficient) pairs.
# z = gpath[tail] - gpath[head] + e_i, gpath the signed root->f tree path.
function fundamental_cycle!(g::HingeGraph, i::Int, out::Vector{Tuple{Int,Float64}})
    empty!(out)
    push!(out, (i, 1.0))
    a = g.tail[i]; b = g.head[i]
    while g.depth[a] > g.depth[b]
        push!(out, (g.pedge[a], +tree_sign(g, a)))
        a = g.parent[a]
    end
    while g.depth[b] > g.depth[a]
        push!(out, (g.pedge[b], -tree_sign(g, b)))
        b = g.parent[b]
    end
    while a != b
        push!(out, (g.pedge[a], +tree_sign(g, a)))
        a = g.parent[a]
        push!(out, (g.pedge[b], -tree_sign(g, b)))
        b = g.parent[b]
    end
    return out
end

"""The hinge graph of the cut structure with its BFS spanning forest."""
function build_hinge_graph(c::CutStructure)
    m = c.mesh
    g = HingeGraph()
    g.F = n_faces(m)
    for e in c.hinge_edges
        ed = m.edges[e]
        push!(g.head, m.half_edges[ed.he[1]].face)
        push!(g.tail, m.half_edges[ed.he[2]].face)
        push!(g.eid, e)
    end
    build_forest!(g)
    return g
end

"""
    subgraph(g, keep_face, keep_edge, pins) -> (HingeGraph, sub_pins)

Same construction on an arbitrary edge subset (used for the 2-core), with the
faces relabelled 1..n; `keep_face` / `keep_edge` are boolean masks.
"""
function subgraph(g::HingeGraph, keep_face::AbstractVector{Bool}, keep_edge::AbstractVector{Bool},
                  pins::Vector{Vec2})
    relabel = zeros(Int, g.F)
    n = 0
    for f in 1:g.F
        keep_face[f] && (n += 1; relabel[f] = n)
    end
    h = HingeGraph()
    h.F = n
    sub_pins = Vec2[]
    for i in 1:n_edges(g)
        keep_edge[i] || continue
        push!(h.head, relabel[g.head[i]])
        push!(h.tail, relabel[g.tail[i]])
        push!(h.eid, g.eid[i])
        push!(sub_pins, pins[i])
    end
    build_forest!(h)
    return h, sub_pins
end

pins_flat(c::CutStructure, g::HingeGraph, X::Vector{Vec2}) =
    Vec2[X[c.hinge_dir[g.eid[i]].src] for i in 1:n_edges(g)]

pins_deployed(c::CutStructure, g::HingeGraph, Y::Vector{Vec2}) =
    Vec2[Y[prime_vertex(c, g.head[i], c.hinge_dir[g.eid[i]].src)] for i in 1:n_edges(g)]

"""The angular-velocity matrix A (2 n_cycles x F), sparse, from the fundamental cycles."""
function build_A(g::HingeGraph, pins::Vector{Vec2})
    nz = n_cycles(g)
    I = Int[]; J = Int[]; V = Float64[]
    cyc = Tuple{Int,Float64}[]
    for k in 1:nz
        fundamental_cycle!(g, g.nontree[k], cyc)
        for (i, z) in cyc
            z == 0 && continue
            push!(I, 2k - 1); push!(J, g.head[i]); push!(V, z * pins[i][1])
            push!(I, 2k - 1); push!(J, g.tail[i]); push!(V, -z * pins[i][1])
            push!(I, 2k);     push!(J, g.head[i]); push!(V, z * pins[i][2])
            push!(I, 2k);     push!(J, g.tail[i]); push!(V, -z * pins[i][2])
        end
    end
    return sparse(I, J, V, 2 * nz, g.F)   # duplicates summed
end

"""
    matrix_rank(A, rel_tol = 1e-10, dense_limit = 700) -> (rank, used_sparse)

Rank: dense column-pivoted Householder QR below `dense_limit` columns (count
|R_ii| > rel_tol * |R_11|), sparse QR (SPQR) above.
"""
function matrix_rank(A::SparseMatrixCSC{Float64,Int}, rel_tol::Float64 = 1e-10,
                     dense_limit::Int = 700)
    (size(A, 1) == 0 || size(A, 2) == 0) && return 0, false
    nrm = nnz(A) == 0 ? 0.0 : maximum(abs, nonzeros(A))
    nrm == 0 && return 0, false
    if size(A, 2) <= dense_limit && size(A, 1) <= 4 * dense_limit
        D = Matrix(A) / nrm
        F = qr(D, ColumnNorm())
        R = F.R
        d = [abs(R[i, i]) for i in 1:min(size(R)...)]
        isempty(d) && return 0, false
        thr = rel_tol * maximum(d)
        return count(>(thr), d), false
    end
    S = A / nrm
    F = qr(S; tol = rel_tol)
    return rank(F), true
end

"""The 2|E_hinge| x 3|F| body-and-pin rigidity matrix, unknowns (omega_f, v_f)."""
function build_rigidity(g::HingeGraph, pins::Vector{Vec2})
    E = n_edges(g)
    R = zeros(2 * E, 3 * g.F)
    for i in 1:E
        Jp = Vec2(-pins[i][2], pins[i][1])
        h = 3 * (g.head[i] - 1); t = 3 * (g.tail[i] - 1)
        R[2i - 1, h + 1] += Jp[1]
        R[2i - 1, t + 1] -= Jp[1]
        R[2i, h + 1] += Jp[2]
        R[2i, t + 1] -= Jp[2]
        R[2i - 1, h + 2] += 1
        R[2i - 1, t + 2] -= 1
        R[2i, h + 3] += 1
        R[2i, t + 3] -= 1
    end
    return R
end

mutable struct TwoCore
    in_core::Vector{Bool}       # per face
    edge_in_core::Vector{Bool}  # per hinge edge
    n_core_faces::Int
    n_dangling::Int  # |F \ core2|, isolated faces included
    c_core::Int      # components of the 2-core (0 if empty)
end

"""The 2-core of the hinge graph by iterated leaf stripping."""
function two_core(g::HingeGraph)
    tc = TwoCore(trues(g.F), trues(n_edges(g)), 0, 0, 0)
    deg = zeros(Int, g.F)
    adj = [Tuple{Int,Int}[] for _ in 1:g.F]
    for i in 1:n_edges(g)
        deg[g.head[i]] += 1
        deg[g.tail[i]] += 1
        push!(adj[g.head[i]], (i, g.tail[i]))
        push!(adj[g.tail[i]], (i, g.head[i]))
    end
    q = Int[f for f in 1:g.F if deg[f] <= 1]
    while !isempty(q)
        f = pop!(q)
        tc.in_core[f] || continue
        tc.in_core[f] = false
        for (i, o) in adj[f]
            tc.edge_in_core[i] || continue
            tc.edge_in_core[i] = false
            if tc.in_core[o]
                deg[o] -= 1
                deg[o] <= 1 && push!(q, o)
            end
        end
    end
    tc.n_core_faces = count(tc.in_core)
    tc.n_dangling = g.F - tc.n_core_faces
    seen = falses(g.F)
    for f in 1:g.F
        (!tc.in_core[f] || seen[f]) && continue
        tc.c_core += 1
        st = Int[f]
        seen[f] = true
        while !isempty(st)
            u = pop!(st)
            for (i, o) in adj[u]
                (!tc.edge_in_core[i] || seen[o]) && continue
                seen[o] = true
                push!(st, o)
            end
        end
    end
    return tc
end

mutable struct MobilityReport
    F::Int; n_hinge::Int; n_cycles::Int; components::Int
    dim_ker_A::Int
    n_dangling::Int
    core_faces::Int; core_edges::Int; core_cycles::Int; c_core::Int
    dim_ker_A_core::Int
    identity_holds::Bool
    m_full::Int  # dim ker A - c(Gamma)
    m_core::Int  # dim ker A|core2 - c(core2)
    sigma_in_ker::Bool
    sigma_residual::Float64
    used_sparse::Bool
end
MobilityReport() = MobilityReport(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false, 0, 0, false, 0.0, false)

"""
    mobility_at(c, g, pins, rel_tol = 1e-10) -> MobilityReport

dim ker A at the given pin positions, the 2-core identity, and whether sigma is in ker A.
"""
function mobility_at(c::CutStructure, g::HingeGraph, pins::Vector{Vec2}, rel_tol::Float64 = 1e-10)
    r = MobilityReport()
    r.F = g.F
    r.n_hinge = n_edges(g)
    r.n_cycles = n_cycles(g)
    r.components = g.components
    A = build_A(g, pins)
    rk, sp = matrix_rank(A, rel_tol)
    r.dim_ker_A = g.F - rk
    r.used_sparse = sp
    r.m_full = r.dim_ker_A - g.components

    tc = two_core(g)
    r.n_dangling = tc.n_dangling
    r.core_faces = tc.n_core_faces
    r.c_core = tc.c_core
    h, cp = subgraph(g, tc.in_core, tc.edge_in_core, pins)
    r.core_edges = n_edges(h)
    r.core_cycles = n_cycles(h)
    if h.F == 0
        r.dim_ker_A_core = 0
        r.m_core = 0
    else
        Ac = build_A(h, cp)
        rkc, spc = matrix_rank(Ac, rel_tol)
        r.dim_ker_A_core = h.F - rkc
        r.used_sparse = r.used_sparse || spc
        r.m_core = r.dim_ker_A_core - h.components
    end
    r.identity_holds = (r.dim_ker_A == r.n_dangling + r.dim_ker_A_core)

    m = c.mesh
    sg = Float64[m.sigma[f] for f in 1:g.F]
    nrm = nnz(A) == 0 ? 0.0 : maximum(abs, nonzeros(A))
    r.sigma_residual = size(A, 1) > 0 ? norm(A * sg) : 0.0
    r.sigma_in_ker = r.sigma_residual <= 1e-8 * max(1.0, nrm) * sqrt(Float64(g.F))
    return r
end
