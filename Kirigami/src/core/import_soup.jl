# core/import_soup.jl -- a welding importer for polygon soup.
#
# The 2025 paper's fabrication SVGs draw every
# tile edge as an independent segment with each tile's corners repeated once per incident
# tile; turning that into a Mesh is three steps:
#   1. WELD   endpoints within `weld_tol` become one vertex (uniform grid hash).
#   2. PRUNE  vertices of degree <= 1 are removed with their edge, repeatedly.
#   3. WALK   the standard planar face traversal: neighbours sorted by angle, successor of
#             (u -> v) is (v -> w) with w immediately clockwise from u at v. Walks of
#             positive signed area are the bounded faces; the outer face is dropped.
# Curved `<path>` elements of an SVG are ignored (counted).

mutable struct SoupReport
    mesh::Mesh
    ok::Bool
    status::String
    n_segments::Int         # input segments, after zero-length and duplicate removal
    n_vertices_welded::Int  # vertices after welding
    n_pruned::Int           # vertices removed as degree <= 1
    n_faces::Int            # bounded faces kept
    n_faces_dropped::Int    # walks dropped as outer face or below the area floor
    min_face_area::Float64
    max_face_area::Float64
    total_area::Float64
end
SoupReport() = SoupReport(Mesh(), false, "not_run", 0, 0, 0, 0, 0, 0.0, 0.0, 0.0)

# Uniform grid hash: welds points that lie within `tol` (inclusive) of one another.
mutable struct Welder
    tol::Float64
    cells::Dict{Int64,Vector{Int}}
    pts::Vector{Vec2}
end
Welder(t::Float64) = Welder(t > 0 ? t : 1e-12, Dict{Int64,Vector{Int}}(), Vec2[])
_welder_key(i::Int64, j::Int64) = i * Int64(2147483647) + j

function add!(w::Welder, p::Vec2)
    ci = floor(Int64, p[1] / w.tol)
    cj = floor(Int64, p[2] / w.tol)
    for di in -1:1, dj in -1:1
        ids = get(w.cells, _welder_key(ci + di, cj + dj), nothing)
        ids === nothing && continue
        for k in ids
            _norm2(w.pts[k] - p) <= w.tol && return k
        end
    end
    push!(w.pts, p)
    id = length(w.pts)
    push!(get!(w.cells, _welder_key(ci, cj), Int[]), id)
    return id
end

# polygon_area (shoelace) is shared with export/layout.jl.

# `std::stod`: optional leading whitespace, then the longest valid floating-point prefix.
const _STOD_RE = r"^\s*[+-]?(?:(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?|inf(?:inity)?|nan)"i
function _stod(s::AbstractString)
    m = match(_STOD_RE, s)
    m === nothing && throw(ArgumentError("stod: no conversion for \"$s\""))
    return parse(Float64, strip(m.match))
end

# `std::stringstream >> double` loop: reads numbers until the first token that is not one.
function parse_points(s::AbstractString)
    nums = Float64[]
    for tok in split(replace(s, ',' => ' '))
        m = match(_STOD_RE, tok)
        m === nothing && break
        push!(nums, parse(Float64, m.match))
    end
    return [Vec2(nums[i], nums[i + 1]) for i in 1:2:(length(nums) - 1)]
end

"""
    weld_segments(segs, weld_tol, area_floor_frac=0.0) -> SoupReport

Welds a set of undirected segments `(a::Vec2, b::Vec2)` into a planar Mesh by the three
steps above. `weld_tol` is absolute; `area_floor_frac` drops faces smaller than that
fraction of the total bounded area.
"""
function weld_segments(segs_in::AbstractVector, weld_tol::Float64, area_floor_frac::Float64 = 0.0)
    rep = SoupReport()
    if isempty(segs_in)
        rep.status = "no_segments"
        return rep
    end

    # ---- 1. weld -------------------------------------------------------------
    w = Welder(weld_tol)
    undirected = Set{Tuple{Int,Int}}()
    for s in segs_in
        a = add!(w, s[1])
        b = add!(w, s[2])
        a == b && continue  # zero length after welding
        push!(undirected, (min(a, b), max(a, b)))
    end
    rep.n_segments = length(undirected)
    rep.n_vertices_welded = length(w.pts)
    if isempty(undirected)
        rep.status = "all_segments_degenerate"
        return rep
    end

    # ---- 2. prune degree <= 1 -------------------------------------------------
    # adjacency as sorted vectors (std::set semantics: `first` = smallest neighbour)
    adj = [Int[] for _ in 1:length(w.pts)]
    for e in sort!(collect(undirected))
        push!(adj[e[1]], e[2])
        push!(adj[e[2]], e[1])
    end
    foreach(sort!, adj)
    stack = [v for v in 1:length(adj) if length(adj[v]) <= 1]
    while !isempty(stack)
        v = pop!(stack)
        isempty(adj[v]) && continue
        u = adj[v][1]
        empty!(adj[v])
        deleteat!(adj[u], searchsortedfirst(adj[u], v))
        rep.n_pruned += 1
        length(adj[u]) <= 1 && push!(stack, u)
    end

    # ---- 3. face walk ---------------------------------------------------------
    # Neighbours of each vertex, sorted CCW by angle; `pos` gives a neighbour's rank.
    nv = length(w.pts)
    nb = [Int[] for _ in 1:nv]
    pos = [Dict{Int,Int}() for _ in 1:nv]
    n_edges = 0
    for v in 1:nv
        isempty(adj[v]) && continue
        c = w.pts[v]
        nb[v] = sort(adj[v]; by = a -> libm_atan2(w.pts[a][2] - c[2], w.pts[a][1] - c[1]))
        for (i, a) in enumerate(nb[v])
            pos[v][a] = i
        end
        n_edges += length(nb[v])
    end
    n_edges ÷= 2
    if n_edges == 0
        rep.status = "everything_pruned"
        return rep
    end

    seen = Set{Tuple{Int,Int}}()
    walks = Vector{Int}[]
    guard = 4 * n_edges + 8
    for u in 1:nv, v in nb[u]
        de = (u, v)
        de in seen && continue
        cyc = Int[]
        cur = de
        while !(cur in seen) && length(cyc) < guard
            push!(seen, cur)
            push!(cyc, cur[1])
            # successor: at cur[2], the neighbour immediately CLOCKWISE from cur[1]
            ring = nb[cur[2]]
            k = pos[cur[2]][cur[1]]
            nk = mod1(k - 1, length(ring))
            cur = (cur[2], ring[nk])
        end
        length(cyc) >= 3 && push!(walks, cyc)
    end

    # ---- keep the bounded faces ----------------------------------------------
    faces = Vector{Int}[]
    areas = Float64[]
    total = 0.0
    for cyc in walks
        a = polygon_area([w.pts[v] for v in cyc])
        if a <= 0
            rep.n_faces_dropped += 1  # outer face / degenerate walk
            continue
        end
        push!(faces, cyc)
        push!(areas, a)
        total += a
    end
    if isempty(faces)
        rep.status = "no_bounded_faces"
        return rep
    end
    if area_floor_frac > 0
        keep = [i for i in eachindex(faces) if areas[i] >= area_floor_frac * total]
        rep.n_faces_dropped += length(faces) - length(keep)
        faces = faces[keep]
        areas = areas[keep]
        if isempty(faces)
            rep.status = "no_faces_above_area_floor"
            return rep
        end
    end

    # ---- compact to the used vertices ----------------------------------------
    remap = zeros(Int, nv)
    m = Mesh()
    for f in faces
        for (k, v) in enumerate(f)
            if remap[v] == 0
                push!(m.X, w.pts[v])
                remap[v] = length(m.X)
            end
            f[k] = remap[v]
        end
    end
    m.faces = faces
    build_topology!(m)

    rep.mesh = m
    rep.n_faces = length(faces)
    rep.min_face_area = minimum(areas)
    rep.max_face_area = maximum(areas)
    rep.total_area = foldl(+, areas; init = 0.0)   # std::accumulate order
    rep.ok = true
    rep.status = "ok"
    return rep
end

"""Same as `weld_segments`, from a closed-polygon soup: each polygon contributes its boundary segments."""
function weld_polygons(polys::AbstractVector, weld_tol::Float64, area_floor_frac::Float64 = 0.0)
    segs = Tuple{Vec2,Vec2}[]
    for p in polys
        n = length(p)
        for i in 1:n
            push!(segs, (p[i], p[mod1(i + 1, n)]))
        end
    end
    return weld_segments(segs, weld_tol, area_floor_frac)
end

"""
    read_svg_segments(path) -> (segments, n_path_ignored)

Reads the drawable straight segments of an SVG: `<line x1 y1 x2 y2>`, `<polyline points>`
and `<polygon points>`. Curved `<path>` elements are IGNORED and counted.
"""
function read_svg_segments(path::AbstractString)
    out = Tuple{Vec2,Vec2}[]
    n_path_ignored = 0
    isfile(path) || return (out, n_path_ignored)
    s = read(path, String)
    tag_re = r"<(line|polyline|polygon|path)\b[^>]*>"
    attr_re = r"(x1|y1|x2|y2|points|d)\s*=\s*\"([^\"]*)\""
    for t in eachmatch(tag_re, s)
        tag = t.captures[1]
        body = t.match
        if tag == "path"
            n_path_ignored += 1
            continue
        end
        at = Dict{String,String}()
        for a in eachmatch(attr_re, body)
            at[a.captures[1]] = a.captures[2]
        end
        if tag == "line"
            if all(k -> haskey(at, k), ("x1", "y1", "x2", "y2"))
                push!(out, (Vec2(_stod(at["x1"]), _stod(at["y1"])),
                            Vec2(_stod(at["x2"]), _stod(at["y2"]))))
            end
        elseif haskey(at, "points")
            P = parse_points(at["points"])
            n = length(P)
            n < 2 && continue
            last = tag == "polygon" ? n : n - 1
            for i in 1:last
                push!(out, (P[i], P[mod1(i + 1, n)]))
            end
        end
    end
    return (out, n_path_ignored)
end

"""Convenience: read an SVG and weld it. `weld_tol_frac` is a fraction of the bounding-box
diagonal of the drawing."""
function import_svg_soup(path::AbstractString, weld_tol_frac::Float64 = 1e-4,
                         area_floor_frac::Float64 = 0.0)
    segs, _ = read_svg_segments(path)
    if isempty(segs)
        r = SoupReport()
        r.status = "no_segments_in_svg"
        return r
    end
    lo = segs[1][1]
    hi = segs[1][1]
    for s in segs
        lo = min.(lo, s[1], s[2])
        hi = max.(hi, s[1], s[2])
    end
    diag = _norm2(hi - lo)
    return weld_segments(segs, max(1e-12, weld_tol_frac * diag), area_floor_frac)
end
