# core/holes.jl -- hole preimage detection (2026 Sec. 4.1, Algorithm 1), an
# independent partition formulation, and a geometric cross-check on M'.
#
# Port of code/src/core/holes.{hpp,cpp}. 1-based indices throughout.

mutable struct HolePreimage
    edges::Vector{Int}     # edge indices of M, sorted
    vertices::Vector{Int}  # the "hole vertices": V(K), sorted
    all_interior::Bool     # true iff every vertex of `vertices` is interior
end
HolePreimage() = HolePreimage(Int[], Int[], false)

mutable struct HoleSet
    all::Vector{HolePreimage}      # every group found, boundary-touching included
    interior_indices::Vector{Int}  # indices into `all` with all_interior == true
end
HoleSet() = HoleSet(HolePreimage[], Int[])
n_interior_holes(hs::HoleSet) = length(hs.interior_indices)

function finalize!(hs::HoleSet, c::CutStructure)
    for h in hs.all
        sort!(h.edges)
        unique!(h.edges)
        sort!(h.vertices)
        unique!(h.vertices)
        h.all_interior = !any(v -> c.mesh.vertex_is_boundary[v], h.vertices)
    end
    sort!(hs.all; by = h -> h.vertices)
    hs.interior_indices = [i for i in eachindex(hs.all) if hs.all[i].all_interior]
    return hs
end

"""
    holes_seed_growing(c::CutStructure) -> HoleSet

Algorithm 1 of the paper, seeded from every hinge edge, implementing the three
PROSE rules (for a hinge edge only its target endpoint is a growth vertex).
"""
function holes_seed_growing(c::CutStructure)
    m = c.mesh
    hs = HoleSet()
    hinge_assigned = falses(n_edges(m))

    # Growth vertices of an edge: for a hinge edge only its target; for a split
    # edge both endpoints (the prose rules are conditions on the shared vertex).
    function growth_vertices(e::Int)
        if c.edge_type[e] == Hinge
            return (c.hinge_dir[e].dst,)
        elseif c.edge_type[e] == Split
            return (m.edges[e].key.a, m.edges[e].key.b)
        end
        return ()
    end

    for seed in c.hinge_edges
        hinge_assigned[seed] && continue
        in_C = Set{Int}(seed)
        checked = Set{Int}()
        queue = Int[seed]
        verts = Set{Int}()
        hinge_assigned[seed] = true

        while !isempty(queue)
            ep = popfirst!(queue)
            ep in checked && continue
            push!(checked, ep)
            for w in growth_vertices(ep)
                push!(verts, w)
                for e in m.vertex_edges[w]
                    e == ep && continue
                    c.edge_type[e] == Border && continue  # line 4: \ E_border
                    if c.edge_type[e] == Split
                        # rule 1 (split-split) and rule 3 with the hinge ep pointing at w
                        if !(e in in_C)
                            push!(in_C, e)
                            push!(queue, e)
                        end
                    else
                        # rule 2 / rule 3: the candidate hinge must point INTO w
                        if c.hinge_dir[e].dst == w
                            push!(in_C, e)
                            hinge_assigned[e] = true
                            # hinge edges are terminal (not enqueued), but their target is w,
                            # already a growth vertex of this preimage.
                        end
                    end
                end
            end
        end
        C = HolePreimage()
        C.edges = sort!(collect(in_C))
        # The hole vertices are the growth vertices reached.
        for e in C.edges, w in growth_vertices(e)
            push!(verts, w)
        end
        C.vertices = sort!(collect(verts))
        push!(hs.all, C)
    end

    # Split edges never reached from a hinge seed still form their own component.
    split_seen = falses(n_edges(m))
    for h in hs.all, e in h.edges
        c.edge_type[e] == Split && (split_seen[e] = true)
    end
    for seed in c.split_edges
        split_seen[seed] && continue
        C = HolePreimage()
        in_C = Set{Int}(seed)
        verts = Set{Int}()
        q = Int[seed]
        while !isempty(q)
            ep = popfirst!(q)
            split_seen[ep] = true
            for w in (m.edges[ep].key.a, m.edges[ep].key.b)
                push!(verts, w)
                for e in m.vertex_edges[w]
                    if c.edge_type[e] == Split && !(e in in_C)
                        push!(in_C, e)
                        push!(q, e)
                    end
                    if c.edge_type[e] == Hinge && c.hinge_dir[e].dst == w
                        push!(in_C, e)
                    end
                end
            end
        end
        C.edges = sort!(collect(in_C))
        C.vertices = sort!(collect(verts))
        push!(hs.all, C)
    end

    finalize!(hs, c)
    return hs
end

"""
    holes_partition(c::CutStructure) -> HoleSet

Independent formulation: components K of the split-edge subgraph (isolated
vertices carrying incoming hinge edges are trivial components);
C_K = E(K) union { hinge edges whose target lies in V(K) }.
"""
function holes_partition(c::CutStructure)
    m = c.mesh
    N = n_vertices(m)
    dsu = DSU(N)
    for e in c.split_edges
        unite!(dsu, m.edges[e].key.a, m.edges[e].key.b)
    end

    groups = Dict{Int,HolePreimage}()
    function touch(v::Int)
        g = get!(HolePreimage, groups, find!(dsu, v))
        push!(g.vertices, v)
        return g
    end
    # split components
    for e in c.split_edges
        g = touch(m.edges[e].key.a)
        push!(g.edges, e)
        touch(m.edges[e].key.b)
    end
    # hinge edges attach to the component of their TARGET; a target with no split
    # edge forms a trivial component {v}.
    for e in c.hinge_edges
        g = touch(c.hinge_dir[e].dst)
        push!(g.edges, e)
    end

    hs = HoleSet()
    for k in sort!(collect(keys(groups)))
        push!(hs.all, groups[k])
    end
    finalize!(hs, c)
    return hs
end

"""
    same_hole_sets(a, b) -> (same, msg)

True iff the two constructions give the same collection of preimages (compared as
sorted lists of edge sets). `msg` is empty when they agree.
"""
function same_hole_sets(a::HoleSet, b::HoleSet)
    key(h) = sort([p.edges for p in h.all])
    ka = key(a)
    kb = key(b)
    ka == kb && return true, ""
    return false, "hole sets differ: $(length(ka)) vs $(length(kb)) preimages"
end

"""
    holes_partition_edges(c, hs) -> (ok, multiply_covered, uncovered)

True iff the preimages partition E_hinge union E_split.
`multiply_covered` / `uncovered` list the offending edges.
"""
function holes_partition_edges(c::CutStructure, hs::HoleSet)
    cnt = zeros(Int, n_edges(c.mesh))
    for h in hs.all, e in h.edges
        cnt[e] += 1
    end
    multiply_covered = Int[]
    uncovered = Int[]
    for e in Iterators.flatten((c.hinge_edges, c.split_edges))
        if cnt[e] == 0
            push!(uncovered, e)
        elseif cnt[e] > 1
            push!(multiply_covered, e)
        end
    end
    return isempty(multiply_covered) && isempty(uncovered), multiply_covered, uncovered
end

"""
    split_subgraph_is_forest(c) -> (ok, cycle_vertices)

Split-edge subgraph acyclicity (Remark A.4's tree claim). `ok` is true if the
split subgraph is a forest; `cycle_vertices` receives the endpoints of the
cycle-closing split edges.
"""
function split_subgraph_is_forest(c::CutStructure)
    dsu = DSU(n_vertices(c.mesh))
    cycle_vertices = Int[]
    for e in c.split_edges
        k = c.mesh.edges[e].key
        if !unite!(dsu, k.a, k.b)
            push!(cycle_vertices, k.a)
            push!(cycle_vertices, k.b)
        end
    end
    return isempty(cycle_vertices), cycle_vertices
end

# Darts: two per M'-edge. Dart 2k-1 is (u->v) as stored by the face, 2k the reverse.
struct Dart
    from::Int
    to::Int
    orig_edge::Int
    face_side::Bool  # true if this dart is the face's own boundary direction
end

function holes_trace(c::CutStructure, Y::Vector{Vec2}, want_cycles::Bool)
    m = c.mesh
    darts = Dart[]
    out_darts = [Int[] for _ in 1:c.n_prime_vertices]

    for f in 1:n_faces(m)
        pv = c.prime_faces[f]
        nf = length(pv)
        for k in 1:nf
            u = pv[k]
            v = pv[mod1(k + 1, nf)]
            oe = m.half_edges[c.face_corner_base[f] + k].edge
            push!(darts, Dart(u, v, oe, true))
            push!(darts, Dart(v, u, oe, false))
            d0 = length(darts) - 1
            push!(out_darts[u], d0)
            push!(out_darts[v], d0 + 1)
        end
    end
    nd = length(darts)
    is_face_dart = [isodd(d) for d in 1:nd]

    # Angular order of outgoing darts at each M'-vertex.
    ang = [(t = Y[d.to] - Y[d.from]; atan(t[2], t[1])) for d in darts]
    for v in 1:c.n_prime_vertices
        sort!(out_darts[v]; by = d -> ang[d])
    end
    # position of a dart inside its vertex's sorted list
    pos = zeros(Int, nd)
    for v in 1:c.n_prime_vertices, (i, d) in enumerate(out_darts[v])
        pos[d] = i
    end

    twin(d) = isodd(d) ? d + 1 : d - 1
    # next dart of the face traversal keeping the face on the left
    function next_dart(d::Int)
        t = twin(d)
        lst = out_darts[darts[t].from]
        n = length(lst)
        return lst[mod1(pos[t] - 1, n)]  # immediately clockwise from the reversed dart
    end

    visited = falses(nd)
    holes = Vector{Int}[]
    for start in 1:nd
        visited[start] && continue
        cycle = Int[]
        d = start
        while !visited[d]
            visited[d] = true
            push!(cycle, d)
            d = next_dart(d)
        end
        d != start && continue  # not a clean cycle (should not happen)
        # signed area
        area = 0.0
        for e in cycle
            p = Y[darts[e].from]
            q = Y[darts[e].to]
            area += p[1] * q[2] - q[1] * p[2]
        end
        area <= 0 && continue  # outer face
        # Is this one of the M' faces? Then all its darts are face darts of one face.
        all_face = all(e -> is_face_dart[e], cycle)
        has_border = any(e -> c.edge_type[darts[e].orig_edge] == Border, cycle)
        all_face && continue
        # Definition 4.1: a hole's boundary is composed of DUPLICATED interior edges.
        # A bounded region bounded partly by the mesh boundary is not a hole.
        has_border && continue
        if want_cycles
            push!(holes, [darts[e].from for e in cycle])
            continue
        end
        oes = [darts[e].orig_edge for e in cycle]
        sort!(oes)
        unique!(oes)
        push!(holes, oes)
    end
    want_cycles || sort!(holes)
    return holes
end

"""
    holes_geometric(c, Y) -> Vector{Vector{Int}}

Geometric hole detection on the deployed kirigami M' with M'-vertex positions `Y`.
Returns, for every bounded complement component bounded ONLY by duplicated interior
edges (Definition 4.1), the sorted set of ORIGINAL M edges whose duplicates bound it.
PRECONDITION: the deployed configuration must be free of face-face overlaps --
the routine traces the planar arrangement using M' vertices as its only nodes.
"""
holes_geometric(c::CutStructure, Y::Vector{Vec2}) = holes_trace(c, Y, false)

"""Same traversal as `holes_geometric`, but returning the M'-vertex cycles bounding each hole (for plotting)."""
holes_geometric_cycles(c::CutStructure, Y::Vector{Vec2}) = holes_trace(c, Y, true)
