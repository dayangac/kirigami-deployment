# core/cut.jl -- edge classification into E_hinge / E_split / E_border and the
# construction of the kirigami structure M' = (X', F') (2026 Sec. 3).
#
# 1-based indices throughout.

@enum EdgeType::UInt8 Border Hinge Split

struct DirectedEdge
    src::Int  # hinge stays at src, dst is duplicated
    dst::Int
end
DirectedEdge() = DirectedEdge(0, 0)

mutable struct CutStructure
    mesh::Mesh

    edge_type::Vector{EdgeType}      # per edge of M
    hinge_dir::Vector{DirectedEdge}  # valid only where edge_type == Hinge
    hinge_edges::Vector{Int}
    split_edges::Vector{Int}
    border_edges::Vector{Int}

    # M' vertices: one per equivalence class of face-corners.
    n_prime_vertices::Int
    corner_to_prime::Vector{Int}    # index by half-edge id (== corner id)
    face_corner_base::Vector{Int}   # half-edge id preceding the first corner of each face
    prime_to_original::Vector{Int}  # M'-vertex -> original vertex of M
    prime_faces::Vector{Vector{Int}}  # |F| faces, M'-vertex indices, same cyclic order

    # Per interior vertex: count of hinge edges directed in / out (Remark A.1).
    hinge_in::Vector{Int}
    hinge_out::Vector{Int}
end

n_hinge(c::CutStructure) = length(c.hinge_edges)
n_split(c::CutStructure) = length(c.split_edges)

# Union-find with path compression; `unite!` returns false if already joined.
struct DSU
    p::Vector{Int}
end
DSU(n::Int) = DSU(collect(1:n))
function find!(d::DSU, x::Int)
    while d.p[x] != x
        d.p[x] = d.p[d.p[x]]
        x = d.p[x]
    end
    return x
end
function unite!(d::DSU, a::Int, b::Int)
    a = find!(d, a)
    b = find!(d, b)
    a == b && return false
    d.p[a] = b
    return true
end

"""The M'-vertex carrying the corner of face `face` at original vertex `v`."""
function prime_vertex(c::CutStructure, face::Int, v::Int)
    vs = c.mesh.faces[face]
    for k in eachindex(vs)
        vs[k] == v && return c.corner_to_prime[c.face_corner_base[face] + k]
    end
    error("prime_vertex: vertex not on face")
end

"""
    make_cut(m::Mesh) -> CutStructure

Classifies edges and builds M'. Requires `m.sigma` to be set and topology built.
"""
function make_cut(m::Mesh)
    isempty(m.sigma) && error("make_cut: sigma not set")
    E = n_edges(m)
    edge_type = fill(Border, E)
    hinge_dir = fill(DirectedEdge(), E)
    hinge_edges = Int[]
    split_edges = Int[]
    border_edges = Int[]

    for e in 1:E
        ed = m.edges[e]
        if boundary(ed)
            push!(border_edges, e)
            continue
        end
        f1 = m.half_edges[ed.he[1]].face
        f2 = m.half_edges[ed.he[2]].face
        d1 = half_edge_direction(m, ed.he[1])
        d2 = half_edge_direction(m, ed.he[2])
        if m.sigma[f1] != m.sigma[f2]
            # Opposite orientations => the two induced half-edges have the SAME
            # direction => hinge cut. Hinge stays at the source, target is duplicated.
            (d1[1] != d2[1] || d1[2] != d2[2]) &&
                error("internal: opposite sigma but half-edges disagree")
            edge_type[e] = Hinge
            hinge_dir[e] = DirectedEdge(d1[1], d1[2])
            push!(hinge_edges, e)
        else
            (d1[1] != d2[2] || d1[2] != d2[1]) &&
                error("internal: equal sigma but half-edges are not antiparallel")
            edge_type[e] = Split
            push!(split_edges, e)
        end
    end

    # Corner ids == half-edge ids (half-edge h of face f at corner c represents the
    # corner at vertex half_edges[h].from).
    n_corners = length(m.half_edges)
    dsu = DSU(n_corners)
    # corner lookup: (face, vertex) -> corner/half-edge id
    at_vertex = [Tuple{Int,Int}[] for _ in 1:n_vertices(m)]
    for h in 1:n_corners
        push!(at_vertex[m.half_edges[h].from], (m.half_edges[h].face, h))
    end
    function corner_of(face::Int, v::Int)
        for (f, h) in at_vertex[v]
            f == face && return h
        end
        error("corner_of: vertex not on face")
    end

    # Two face-corners at v are identified iff their shared edge is a hinge cut
    # whose SOURCE is v.
    for e in hinge_edges
        ed = m.edges[e]
        v = hinge_dir[e].src
        f1 = m.half_edges[ed.he[1]].face
        f2 = m.half_edges[ed.he[2]].face
        unite!(dsu, corner_of(f1, v), corner_of(f2, v))
    end

    corner_to_prime = zeros(Int, n_corners)
    prime_to_original = Int[]
    next = 0
    for h in 1:n_corners
        r = find!(dsu, h)
        if corner_to_prime[r] == 0
            next += 1
            corner_to_prime[r] = next
            push!(prime_to_original, m.half_edges[r].from)
        end
        corner_to_prime[h] = corner_to_prime[r]
    end

    prime_faces = [Int[] for _ in 1:n_faces(m)]
    face_corner_base = zeros(Int, n_faces(m))
    base = 0
    for f in 1:n_faces(m)
        face_corner_base[f] = base
        nf = length(m.faces[f])
        for k in 1:nf
            push!(prime_faces[f], corner_to_prime[base + k])
        end
        base += nf
    end

    hinge_in = zeros(Int, n_vertices(m))
    hinge_out = zeros(Int, n_vertices(m))
    for e in hinge_edges
        hinge_out[hinge_dir[e].src] += 1
        hinge_in[hinge_dir[e].dst] += 1
    end
    return CutStructure(m, edge_type, hinge_dir, hinge_edges, split_edges, border_edges,
                        next, corner_to_prime, face_corner_base, prime_to_original,
                        prime_faces, hinge_in, hinge_out)
end

"""
    check_remark_A1(c) -> (ok, bad)

Remark A.1: at every interior vertex the number of hinge edges directed in
equals the number directed out. `ok` is true if it holds everywhere;
`bad` lists the offending vertices.
"""
function check_remark_A1(c::CutStructure)
    bad = Int[]
    m = c.mesh
    for v in 1:n_vertices(m)
        m.vertex_is_boundary[v] && continue
        c.hinge_in[v] != c.hinge_out[v] && push!(bad, v)
    end
    return isempty(bad), bad
end

# faces joined through shared M'-vertices: node f for face f, node F+pv for M'-vertex pv
function face_prime_dsu(c::CutStructure)
    F = n_faces(c.mesh)
    dsu = DSU(F + c.n_prime_vertices)
    for f in 1:F, pv in c.prime_faces[f]
        unite!(dsu, f, F + pv)
    end
    return dsu
end

"""Per-face component label of M' (1..count_components)."""
function face_components(c::CutStructure)
    F = n_faces(c.mesh)
    dsu = face_prime_dsu(c)
    label = zeros(Int, F)
    root_id = zeros(Int, F + c.n_prime_vertices)
    next = 0
    for f in 1:F
        r = find!(dsu, f)
        if root_id[r] == 0
            next += 1
            root_id[r] = next
        end
        label[f] = root_id[r]
    end
    return label
end

"""Number of connected components of M' (faces joined through shared M'-vertices)."""
function count_components(c::CutStructure)
    F = n_faces(c.mesh)
    dsu = face_prime_dsu(c)
    seen = falses(F + c.n_prime_vertices)
    comps = 0
    for f in 1:F
        r = find!(dsu, f)
        if !seen[r]
            seen[r] = true
            comps += 1
        end
    end
    return comps
end
