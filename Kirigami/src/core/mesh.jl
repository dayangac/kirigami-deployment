# core/mesh.jl -- planar graph data model, JSON I/O, half-edge structure and the
# face-orientation (sigma) induced half-edge direction of Segall et al. 2026 Sec. 3.
#
# All indices are 1-based internally; the JSON file format is 0-based and is converted at
# the I/O boundary. Functions that mutate their argument carry the Julia `!` suffix
# (`build_topology!`, `normalize_face_ccw!`).

# 2-D cross product (signed parallelogram area); shared by several method modules.
_det2(u::Vec2, v::Vec2) = u[1] * v[2] - u[2] * v[1]

# Platform libm shims. The frozen corpora and fixtures were produced with Apple's arm64
# libm; Julia's Base trig differs from it by 1 ulp on some arguments, which is visible in
# bit-exact reproductions (tiling generators, float32 STL normals). On macOS call the
# system libm; elsewhere use Base (docs/NUMERICS.md).
# True only when the system libm is the arm64 Apple libm (an x86_64 Julia under Rosetta
# gets Apple's x86_64 libm, which differs by 1 ulp on some arguments); bit-exact
# trig-dependent tests key on this.
const _USE_SYSTEM_LIBM = Sys.isapple() && Sys.ARCH === :aarch64
if Sys.isapple()
    const _LIBM = "libSystem.B.dylib"
    libm_cos(x::Float64) = ccall((:cos, _LIBM), Float64, (Float64,), x)
    libm_sin(x::Float64) = ccall((:sin, _LIBM), Float64, (Float64,), x)
    libm_tan(x::Float64) = ccall((:tan, _LIBM), Float64, (Float64,), x)
    libm_pow(x::Float64, y::Float64) = ccall((:pow, _LIBM), Float64, (Float64, Float64), x, y)
    libm_atan2(y::Float64, x::Float64) = ccall((:atan2, _LIBM), Float64, (Float64, Float64), y, x)
else
    libm_cos(x::Float64) = cos(x)
    libm_sin(x::Float64) = sin(x)
    libm_tan(x::Float64) = tan(x)
    libm_pow(x::Float64, y::Float64) = x^y
    libm_atan2(y::Float64, x::Float64) = atan(y, x)
end

# Where BOTH sin(x) and cos(x) of one argument are needed in one function (the rotation of
# kinematics.jl, the harmonic p + q cos + r sin, ...), the reference values come from
# Apple libm's combined `__sincos_stret`, whose sine differs from the standalone `sin` by
# 1 ulp on ~4% of arguments.  `libm_sincos(x) -> (sin, cos)` calls that same routine.
struct _SinCosRet
    s::Float64
    c::Float64
end
if Sys.isapple()
    function libm_sincos(x::Float64)
        r = ccall((:__sincos_stret, _LIBM), _SinCosRet, (Float64,), x)
        return r.s, r.c
    end
else
    libm_sincos(x::Float64) = sincos(x)
end

# A pair of vertex indices in canonical (sorted) order -- the undirected edge key.
struct EdgeKey
    a::Int
    b::Int
    EdgeKey(u::Integer, v::Integer) = u < v ? new(u, v) : new(v, u)
end
Base.isless(x::EdgeKey, y::EdgeKey) = x.a != y.a ? x.a < y.a : x.b < y.b
# == and hash are the default structural ones (immutable struct), so EdgeKey works as
# a Dict key and `edge_index[EdgeKey(u, v)]` finds the edge regardless of order.

# One half-edge: the occurrence of an undirected edge inside one face.
# `from`/`to` are the vertices in the face's *stored* (given, geometrically CCW)
# cyclic order; the sigma-induced direction is computed separately.
mutable struct HalfEdge
    face::Int
    corner::Int   # index within the face's vertex list of `from`
    from::Int
    to::Int
    edge::Int     # index into Mesh.edges
end
HalfEdge() = HalfEdge(0, 0, 0, 0, 0)

mutable struct Edge
    key::EdgeKey
    # half-edge indices of the (1 or 2) faces incident to this edge (0 = unset)
    he::Vector{Int}
    n_faces::Int
end
Edge(key::EdgeKey) = Edge(key, [0, 0], 0)
boundary(e::Edge) = e.n_faces == 1

mutable struct PeriodicSpec
    present::Bool
    th::Vec2
    tv::Vec2
    pairs_h::Vector{Tuple{Int,Int}}
    pairs_v::Vector{Tuple{Int,Int}}
end
PeriodicSpec() = PeriodicSpec(false, Vec2(0, 0), Vec2(0, 0), Tuple{Int,Int}[], Tuple{Int,Int}[])

# The uncut embedded planar graph M = (X, F).
mutable struct Mesh
    X::Vector{Vec2}                    # vertex positions
    faces::Vector{Vector{Int}}         # each face CCW in the given geometry
    sigma::Vector{Int}                 # +1 clockwise, -1 counter-clockwise; may be empty
    periodic::PeriodicSpec

    # derived (filled by build_topology!)
    half_edges::Vector{HalfEdge}
    edges::Vector{Edge}
    edge_index::Dict{EdgeKey,Int}
    vertex_half_edges::Vector{Vector{Int}}  # outgoing half-edges at each vertex
    vertex_edges::Vector{Vector{Int}}       # incident edge indices
    vertex_faces::Vector{Vector{Int}}       # incident face indices
    vertex_is_boundary::Vector{Bool}
end
Mesh() = Mesh(Vec2[], Vector{Int}[], Int[], PeriodicSpec(), HalfEdge[], Edge[],
              Dict{EdgeKey,Int}(), Vector{Int}[], Vector{Int}[], Vector{Int}[], Bool[])
Mesh(X::Vector{Vec2}, faces::Vector{Vector{Int}}) =
    (m = Mesh(); m.X = X; m.faces = faces; m)

n_vertices(m::Mesh) = length(m.X)
n_faces(m::Mesh) = length(m.faces)
n_edges(m::Mesh) = length(m.edges)
n_interior_vertices(m::Mesh) = count(!, m.vertex_is_boundary)

"""
    build_topology!(m::Mesh)

Builds half_edges/edges/vertex_* and the boundary flags.
Throws on non-manifold input (edge with >2 faces).
"""
function build_topology!(m::Mesh)
    empty!(m.half_edges)
    empty!(m.edges)
    empty!(m.edge_index)
    N = n_vertices(m)
    m.vertex_half_edges = [Int[] for _ in 1:N]
    m.vertex_edges = [Int[] for _ in 1:N]
    m.vertex_faces = [Int[] for _ in 1:N]
    m.vertex_is_boundary = falses(N)

    for f in 1:n_faces(m)
        vs = m.faces[f]
        length(vs) < 3 && error("face with fewer than 3 vertices")
        nv = length(vs)
        for c in 1:nv
            u = vs[c]
            v = vs[mod1(c + 1, nv)]
            (u < 1 || u > N || v < 1 || v > N) && error("face references out-of-range vertex")
            u == v && error("degenerate face edge (u == v)")
            k = EdgeKey(u, v)
            ei = get(m.edge_index, k, 0)
            if ei == 0
                push!(m.edges, Edge(k))
                ei = length(m.edges)
                m.edge_index[k] = ei
            end
            push!(m.half_edges, HalfEdge(f, c, u, v, ei))
            hi = length(m.half_edges)
            e = m.edges[ei]
            if e.n_faces >= 2
                # error message reports 0-based vertex ids (the JSON file's numbering)
                error("non-manifold input: edge ($(k.a - 1),$(k.b - 1)) is shared by more than two faces")
            end
            e.n_faces += 1
            e.he[e.n_faces] = hi
            push!(m.vertex_half_edges[u], hi)
        end
        for v in vs
            push!(m.vertex_faces[v], f)
        end
    end
    for fs in m.vertex_faces
        sort!(fs)
        unique!(fs)
    end
    for ei in 1:n_edges(m)
        k = m.edges[ei].key
        push!(m.vertex_edges[k.a], ei)
        push!(m.vertex_edges[k.b], ei)
        if boundary(m.edges[ei])
            m.vertex_is_boundary[k.a] = true
            m.vertex_is_boundary[k.b] = true
        end
    end
    # Isolated vertices (no incident edge) count as boundary -- they carry no constraint.
    for v in 1:N
        isempty(m.vertex_edges[v]) && (m.vertex_is_boundary[v] = true)
    end
    return m
end

"""Signed area of face `f` using the stored vertex order (>0 iff CCW)."""
function face_signed_area(m::Mesh, f::Int)
    vs = m.faces[f]
    n = length(vs)
    a = 0.0
    for i in 1:n
        p = m.X[vs[i]]
        q = m.X[vs[mod1(i + 1, n)]]
        a += p[1] * q[2] - q[1] * p[2]
    end
    return 0.5 * a
end

"""Interior angle of face `f` at its corner `c` (index into faces[f])."""
function corner_angle(m::Mesh, f::Int, c::Int)
    vs = m.faces[f]
    n = length(vs)
    p = m.X[vs[mod1(c - 1, n)]]
    o = m.X[vs[c]]
    q = m.X[vs[mod1(c + 1, n)]]
    a = p - o
    b = q - o
    # Interior angle measured inside the polygon. The stored order is CCW, so the
    # interior lies to the left of o->q; the interior angle is the CCW turn from
    # (o->q) to (o->p).
    ang = atan(b[1] * a[2] - b[2] * a[1], dot(b, a))
    ang < 0 && (ang += 2pi)
    return ang
end

"""Interior angle of face `f` at vertex `v` (v must be a vertex of f)."""
function face_angle_at(m::Mesh, f::Int, v::Int)
    vs = m.faces[f]
    for c in eachindex(vs)
        vs[c] == v && return corner_angle(m, f, c)
    end
    error("face_angle_at: vertex not on face")
end

"""
    half_edge_direction(m, h) -> (src, dst)

The sigma-induced direction of the half-edge `h`.
sigma(f) = -1 (counter-clockwise) traverses the stored CCW order;
sigma(f) = +1 (clockwise) traverses it reversed.
"""
function half_edge_direction(m::Mesh, h::Int)
    he = m.half_edges[h]
    isempty(m.sigma) && error("half_edge_direction: sigma not set")
    s = m.sigma[he.face]
    (s != 1 && s != -1) && error("sigma must be +1 or -1")
    # sigma = -1 (CCW) keeps the stored order; sigma = +1 (CW) reverses it.
    s == -1 && return (he.from, he.to)
    return (he.to, he.from)
end

"""Same as `half_edge_direction(m, h)`, addressed by (face, edge)."""
function half_edge_direction(m::Mesh, face::Int, edge::Int)
    e = m.edges[edge]
    for i in 1:e.n_faces
        m.half_edges[e.he[i]].face == face && return half_edge_direction(m, e.he[i])
    end
    error("half_edge_direction: face not incident to edge")
end

"""
    normalize_face_ccw!(m::Mesh)

Ensures every face's stored vertex order is counter-clockwise in the geometry,
which the sigma-direction rule assumes.
"""
function normalize_face_ccw!(m::Mesh)
    for f in eachindex(m.faces)
        vs = m.faces[f]
        n = length(vs)
        a = 0.0
        for i in 1:n
            p = m.X[vs[i]]
            q = m.X[vs[mod1(i + 1, n)]]
            a += p[1] * q[2] - q[1] * p[2]
        end
        a < 0 && reverse!(m.faces[f])
    end
    return m
end

# ---------------------------------------------------------------------------- JSON I/O

function mesh_from_json(j::AbstractDict)
    m = Mesh()
    for v in j["vertices"]
        length(v) != 2 && error("vertex must have 2 coordinates")
        push!(m.X, Vec2(Float64(v[1]), Float64(v[2])))
    end
    for f in j["faces"]
        # file is 0-based -> internal 1-based
        push!(m.faces, [Int(i) + 1 for i in f])
    end
    if haskey(j, "orientation") && j["orientation"] !== nothing
        m.sigma = [Int(s) for s in j["orientation"]]
        length(m.sigma) != length(m.faces) && error("orientation array length != number of faces")
    end
    if haskey(j, "periodic") && j["periodic"] !== nothing
        p = j["periodic"]
        m.periodic.present = true
        m.periodic.th = Vec2(Float64(p["th"][1]), Float64(p["th"][2]))
        m.periodic.tv = Vec2(Float64(p["tv"][1]), Float64(p["tv"][2]))
        # vertex pairs are 0-based in the file -> 1-based internally
        if haskey(p, "pairs_h")
            for q in p["pairs_h"]
                push!(m.periodic.pairs_h, (Int(q[1]) + 1, Int(q[2]) + 1))
            end
        end
        if haskey(p, "pairs_v")
            for q in p["pairs_v"]
                push!(m.periodic.pairs_v, (Int(q[1]) + 1, Int(q[2]) + 1))
            end
        end
    end
    normalize_face_ccw!(m)
    build_topology!(m)
    return m
end

function mesh_to_json(m::Mesh)
    j = Dict{String,Any}()
    j["vertices"] = [[x[1], x[2]] for x in m.X]
    # internal 1-based -> 0-based in the file
    j["faces"] = [[i - 1 for i in f] for f in m.faces]
    isempty(m.sigma) || (j["orientation"] = copy(m.sigma))
    if m.periodic.present
        p = Dict{String,Any}()
        p["th"] = [m.periodic.th[1], m.periodic.th[2]]
        p["tv"] = [m.periodic.tv[1], m.periodic.tv[2]]
        # 1-based -> 0-based
        p["pairs_h"] = [[q[1] - 1, q[2] - 1] for q in m.periodic.pairs_h]
        p["pairs_v"] = [[q[1] - 1, q[2] - 1] for q in m.periodic.pairs_v]
        j["periodic"] = p
    end
    return j
end

"""Loads a mesh JSON file (0-based indices in the file), normalises faces CCW and builds topology."""
function load_mesh_json(path::AbstractString)
    isfile(path) || error("cannot open $path")
    return mesh_from_json(JSON.parsefile(path))
end

"""Writes the mesh as JSON in the project file format (0-based indices)."""
function save_mesh_json(m::Mesh, path::AbstractString)
    open(path, "w") do io
        JSON.print(io, mesh_to_json(m), 2)
        write(io, "\n")
    end
    return nothing
end

mesh_from_json_string(text::AbstractString) = mesh_from_json(JSON.parse(text))
mesh_to_json_string(m::Mesh) = JSON.json(mesh_to_json(m), 2)
