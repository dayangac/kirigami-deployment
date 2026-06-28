# core/generators.jl -- test-graph generators: regular / Archimedean tilings clipped
# to a region, and random planar graphs (own Bowyer-Watson Delaunay, half-plane Voronoi,
# greedy quad merge).
#
# Port of code/src/core/generators.{hpp,cpp}. The random generators are meant to reproduce
# the C++ meshes BIT-EXACTLY for a given `MT19937` (same draw order, same welding order,
# same face order), so the frozen populations in data/corpus can be regenerated. Two
# things are needed for that beyond a literal translation:
#   * clang -O2 on arm64 contracts `a*b + c` written in ONE plain-double expression into
#     an fma; Eigen expressions (`v.dot(w)`, `v.norm()`, `A + (B-A)*t`) are NOT contracted
#     because the multiply and the add live in different (inlined) functions. `fma` is
#     used below exactly where the C++ has a plain-double `a*b ± c*d` / `a*b + c`.
#   * `std::cos/sin/tan/pow/atan2` are Apple's arm64 libm, which is not correctly rounded
#     (~8% of arguments) and is not what an x86_64 (Rosetta) Julia reaches through ccall:
#     the x86_64 slice of libsystem_m differs from the arm64 slice on ~33% of arguments,
#     Julia's Base trig on ~8%. So `libm_*` below call the system libm only on an arm64
#     Apple build (then tilings and the relaxation are bit-exact with the C++) and fall
#     back to Base elsewhere (then tiling vertices can differ from the C++ by 1 ulp).

const _LIBM = "/usr/lib/libSystem.B.dylib"
const _USE_SYSTEM_LIBM = Sys.isapple() && Sys.ARCH === :aarch64
if _USE_SYSTEM_LIBM
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

# Eigen `v.dot(w)`, `v.squaredNorm()`, `v.norm()` for Vector2d: x*x' + y*y', not contracted.
_dot2(a::Vec2, b::Vec2) = a[1] * b[1] + a[2] * b[2]
_sqnorm2(a::Vec2) = a[1] * a[1] + a[2] * a[2]
_norm2(a::Vec2) = sqrt(_sqnorm2(a))

const kSqrt3 = 1.7320508075688772

# Welds points closer than `tol` (grid hash of cell size tol, 3x3 neighbourhood). The
# grid stores ONE id per cell and `add!` overwrites it, exactly like the C++ std::map
# assignment; replicated on purpose.
mutable struct VertexWelder
    tol::Float64
    grid::Dict{Tuple{Int64,Int64},Int}
    pts::Vector{Vec2}
end
VertexWelder(tol::Float64) = VertexWelder(tol, Dict{Tuple{Int64,Int64},Int}(), Vec2[])

"""Returns the (1-based) id of the welded point equal to `p` within tol, adding it if new."""
function add!(w::VertexWelder, p::Vec2)
    # std::llround: nearest, ties away from zero
    gx = round(Int64, p[1] / w.tol, RoundNearestTiesAway)
    gy = round(Int64, p[2] / w.tol, RoundNearestTiesAway)
    for dx in -1:1, dy in -1:1
        id = get(w.grid, (gx + dx, gy + dy), 0)
        if id != 0 && _norm2(w.pts[id] - p) < w.tol
            return id
        end
    end
    push!(w.pts, p)
    id = length(w.pts)
    w.grid[(gx, gy)] = id
    return id
end

function centroid(p::Vector{Vec2})
    c = Vec2(0.0, 0.0)
    for q in p
        c += q
    end
    return c / Float64(length(p))
end

function clip_polys(polys::Vector{Vector{Vec2}}, r)
    return [p for p in polys if contains(r, centroid(p))]
end

function polar(r::Float64, deg::Float64)
    a = deg * pi / 180.0
    return Vec2(r * libm_cos(a), r * libm_sin(a))
end

@enum ClipKind Disk Rect

struct ClipRegion
    kind::ClipKind
    center::Vec2
    radius::Float64
    half_w::Float64
    half_h::Float64
end
disk(c::Vec2, r::Real) = ClipRegion(Disk, c, Float64(r), 1.0, 1.0)
rect(c::Vec2, hw::Real, hh::Real) = ClipRegion(Rect, c, 1.0, Float64(hw), Float64(hh))

function Base.contains(r::ClipRegion, p::Vec2)
    r.kind == Disk && return _norm2(p - r.center) <= r.radius
    return abs(p[1] - r.center[1]) <= r.half_w && abs(p[2] - r.center[2]) <= r.half_h
end

# Range of lattice indices needed to cover the clip region generously.
function lattice_extent(r::ClipRegion, step::Float64)
    reach = r.kind == Disk ? r.radius : hypot(r.half_w, r.half_h)
    return Int(ceil(reach / max(step, 1e-9))) + 3
end

"""
    mesh_from_polygons(polys; weld_tol=1e-6) -> Mesh

Welds a polygon soup into a Mesh (faces normalized CCW, topology built). Vertex order is
the order of first appearance; duplicate faces (same vertex set) are dropped.
"""
function mesh_from_polygons(polys::Vector{Vector{Vec2}}, weld_tol::Float64 = 1e-6)
    w = VertexWelder(weld_tol)
    m = Mesh()
    seen = Set{Vector{Int}}()
    for p in polys
        f = Int[]
        for q in p
            id = add!(w, q)
            (isempty(f) || f[end] != id) && push!(f, id)
        end
        while length(f) > 1 && f[1] == f[end]
            pop!(f)
        end
        length(f) < 3 && continue
        key = sort(f)
        key in seen && continue  # duplicate face
        push!(seen, key)
        push!(m.faces, f)
    end
    m.X = w.pts
    normalize_face_ccw!(m)
    build_topology!(m)
    return m
end

"""Keeps only the largest edge-connected component of faces (first largest on ties)."""
function largest_component(m::Mesh)
    F = n_faces(m)
    comp = fill(0, F)
    nc = 0
    adj = [Int[] for _ in 1:F]
    for e in m.edges
        if e.n_faces == 2
            a = m.half_edges[e.he[1]].face
            b = m.half_edges[e.he[2]].face
            push!(adj[a], b)
            push!(adj[b], a)
        end
    end
    for f in 1:F
        comp[f] != 0 && continue
        nc += 1
        stack = [f]
        comp[f] = nc
        while !isempty(stack)
            g = pop!(stack)
            for h in adj[g]
                if comp[h] == 0
                    comp[h] = nc
                    push!(stack, h)
                end
            end
        end
    end
    size = zeros(Int, nc)
    for f in 1:F
        size[comp[f]] += 1
    end
    best = argmax(size)
    polys = Vector{Vec2}[]
    for f in 1:F
        comp[f] != best && continue
        push!(polys, [m.X[v] for v in m.faces[f]])
    end
    return mesh_from_polygons(polys)
end

# ---------------------------------------------------------------- tilings

# Triangular-lattice point (i, j); `i + 0.5*j` is contracted in the C++ but 0.5*j is exact.
_tri_P(i::Int, j::Int) = Vec2(i + 0.5 * j, kSqrt3 * 0.5 * j)
_tri_mid(i1, j1, i2, j2) = 0.5 * (_tri_P(i1, j1) + _tri_P(i2, j2))

function tiling_triangles(r::ClipRegion)
    K = lattice_extent(r, 1.0)
    polys = Vector{Vec2}[]
    for i in -K:K, j in -K:K
        push!(polys, [_tri_P(i, j), _tri_P(i + 1, j), _tri_P(i, j + 1)])
        push!(polys, [_tri_P(i + 1, j), _tri_P(i + 1, j + 1), _tri_P(i, j + 1)])
    end
    return largest_component(mesh_from_polygons(clip_polys(polys, r)))
end

function tiling_squares(r::ClipRegion)
    K = lattice_extent(r, 1.0)
    polys = Vector{Vec2}[]
    for i in -K:K, j in -K:K
        push!(polys, [Vec2(i, j), Vec2(i + 1, j), Vec2(i + 1, j + 1), Vec2(i, j + 1)])
    end
    return largest_component(mesh_from_polygons(clip_polys(polys, r)))
end

function tiling_hexagons(r::ClipRegion)
    K = lattice_extent(r, kSqrt3)
    polys = Vector{Vec2}[]
    for i in -K:K, j in -K:K
        c = Vec2(1.5 * i, kSqrt3 * 0.5 * i) + Vec2(0.0, kSqrt3 * j)
        push!(polys, [c + polar(1.0, 60.0 * k) for k in 0:5])
    end
    return largest_component(mesh_from_polygons(clip_polys(polys, r)))
end

const _KAGOME_DIRS = ((1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1))

# Medial tiling of the triangular lattice: hexagon per lattice vertex, triangle per
# lattice triangle; vertices are the midpoints of triangular-lattice edges.
function _kagome_polys!(polys::Vector{Vector{Vec2}}, i::Int, j::Int)
    push!(polys, [_tri_mid(i, j, i + d[1], j + d[2]) for d in _KAGOME_DIRS])
    # up triangle (i,j),(i+1,j),(i,j+1)
    push!(polys, [_tri_mid(i, j, i + 1, j), _tri_mid(i + 1, j, i, j + 1), _tri_mid(i, j + 1, i, j)])
    # down triangle (i+1,j),(i+1,j+1),(i,j+1)
    push!(polys, [_tri_mid(i + 1, j, i + 1, j + 1), _tri_mid(i + 1, j + 1, i, j + 1),
                  _tri_mid(i, j + 1, i + 1, j)])
    return polys
end

"""3.6.3.6 trihexagonal (kagome) tiling."""
function tiling_kagome(r::ClipRegion)
    K = lattice_extent(r, 1.0)
    polys = Vector{Vec2}[]
    for i in -K:K, j in -K:K
        _kagome_polys!(polys, i, j)
    end
    return largest_component(mesh_from_polygons(clip_polys(polys, r)))
end

"""2-uniform tiling [3.4.3.12; 3.12.12], unit edge length. Dodecagons on a square lattice of
spacing 2*apothem, vertices at angles -15 + 30k; the 4 axis-aligned edges are shared with
the neighbouring dodecagons, the other 8 carry triangles, one square per diagonal gap."""
function tiling_3_4_3_12(r::ClipRegion)
    a12 = 0.5 / libm_tan(pi / 12.0)   # apothem = 1.8660254
    R12 = 0.5 / libm_sin(pi / 12.0)   # circumradius = 1.9318517
    s = 2 * a12
    K = lattice_extent(r, s)
    polys = Vector{Vec2}[]
    for i in -K:K, j in -K:K
        c = Vec2(i * s, j * s)
        dod = [c + polar(R12, -15.0 + 30.0 * k) for k in 0:11]
        push!(polys, dod)
        for k in 0:11
            k % 3 == 0 && continue  # axis edges are shared with dodecagons
            mid_deg = 30.0 * k
            apex = c + polar(a12 + kSqrt3 * 0.5, mid_deg)
            push!(polys, [dod[k + 1], dod[(k + 1) % 12 + 1], apex])
        end
        sc = c + Vec2(a12, a12)
        push!(polys, [sc + Vec2(-0.5, -0.5), sc + Vec2(0.5, -0.5), sc + Vec2(0.5, 0.5),
                      sc + Vec2(-0.5, 0.5)])
    end
    return largest_component(mesh_from_polygons(clip_polys(polys, r)))
end

"""3.3.4.3.4 snub square: squares on a square lattice with orientation alternating 0/30
degrees by parity; the triangles are the 3-cliques of the unit-distance graph."""
function tiling_snub_square(r::ClipRegion)
    A = (1 + kSqrt3) / 4
    B = (3 + kSqrt3) / 4
    u = Vec2(A, B)
    v = Vec2(B, -A)
    K = lattice_extent(r, _norm2(u))
    squares = Vector{Vec2}[]
    for i in -K:K, j in -K:K
        # Eigen `i*u + j*v`: coefficient-wise i*u.x + j*v.x, not contracted
        c = Vec2(Float64(i) * u[1] + Float64(j) * v[1], Float64(i) * u[2] + Float64(j) * v[2])
        rot = ((i + j) & 1) != 0 ? 30.0 : 0.0
        push!(squares, [c + polar(sqrt(2.0) / 2, rot + 45.0 + 90.0 * k) for k in 0:3])
    end
    w = VertexWelder(1e-6)
    sq_idx = [[add!(w, p) for p in q] for q in squares]
    pts = w.pts
    n = length(pts)
    nbr = [Int[] for _ in 1:n]
    for i in 1:n, j in (i + 1):n
        if abs(_norm2(pts[i] - pts[j]) - 1.0) < 1e-6
            push!(nbr[i], j)
            push!(nbr[j], i)
        end
    end
    foreach(sort!, nbr)  # std::set iteration order
    polys = Vector{Vec2}[]
    for f in sq_idx
        push!(polys, [pts[i] for i in f])
    end
    for i in 1:n, j in nbr[i]
        j <= i && continue
        for k in nbr[j]
            k <= j && continue
            k in nbr[i] && push!(polys, [pts[i], pts[j], pts[k]])
        end
    end
    return largest_component(mesh_from_polygons(clip_polys(polys, r)))
end

"""4.8.8 truncated square: octagons on a square lattice of spacing 1 + sqrt(2)."""
function tiling_truncated_square(r::ClipRegion)
    R8 = 0.5 / libm_sin(pi / 8.0)
    s = 1 + sqrt(2.0)
    K = lattice_extent(r, s)
    polys = Vector{Vec2}[]
    for i in -K:K, j in -K:K
        c = Vec2(i * s, j * s)
        push!(polys, [c + polar(R8, 22.5 + 45.0 * k) for k in 0:7])
        sc = c + Vec2(s / 2, s / 2)
        push!(polys, [sc + polar(sqrt(2.0) / 2, 90.0 * k) for k in 0:3])
    end
    return largest_component(mesh_from_polygons(clip_polys(polys, r)))
end

# ---------------------------------------------------------------- periodic

"""Detects periodic vertex pairs (j, i) with x_j - x_i = th (pairs_h) / tv (pairs_v)."""
function add_periodic_pairs!(m::Mesh, th::Vec2, tv::Vec2, tol::Float64 = 1e-6)
    w = VertexWelder(tol)
    for p in m.X
        add!(w, p)
    end
    m.periodic.present = true
    m.periodic.th = th
    m.periodic.tv = tv
    empty!(m.periodic.pairs_h)
    empty!(m.periodic.pairs_v)
    N = n_vertices(m)
    for (t, out) in ((th, m.periodic.pairs_h), (tv, m.periodic.pairs_v))
        for i in 1:N
            j = add!(w, m.X[i] + t)   # may append a new (unmatched) point, as in the C++
            (j <= N && j != i) && push!(out, (j, i))  # x_j - x_i = t
        end
    end
    return m
end

function periodic_squares(n::Int, m_::Int)
    polys = Vector{Vec2}[]
    for i in 0:(n - 1), j in 0:(m_ - 1)
        push!(polys, [Vec2(i, j), Vec2(i + 1, j), Vec2(i + 1, j + 1), Vec2(i, j + 1)])
    end
    g = mesh_from_polygons(polys)
    add_periodic_pairs!(g, Vec2(n, 0), Vec2(0, m_))
    return g
end

function periodic_triangles(n::Int, m_::Int)
    polys = Vector{Vec2}[]
    for i in 0:(n - 1), j in 0:(m_ - 1)
        push!(polys, [_tri_P(i, j), _tri_P(i + 1, j), _tri_P(i, j + 1)])
        push!(polys, [_tri_P(i + 1, j), _tri_P(i + 1, j + 1), _tri_P(i, j + 1)])
    end
    g = mesh_from_polygons(polys)
    add_periodic_pairs!(g, Vec2(n, 0), Vec2(0.5 * m_, kSqrt3 * 0.5 * m_))
    return g
end

function periodic_hexagons(n::Int, m_::Int)
    polys = Vector{Vec2}[]
    for i in 0:(n - 1), j in 0:(m_ - 1)
        c = Vec2(1.5 * i, kSqrt3 * 0.5 * i) + Vec2(0.0, kSqrt3 * j)
        push!(polys, [c + polar(1.0, 60.0 * k) for k in 0:5])
    end
    g = mesh_from_polygons(polys)
    add_periodic_pairs!(g, Vec2(1.5 * n, kSqrt3 * 0.5 * n), Vec2(0, kSqrt3 * m_))
    return g
end

function periodic_kagome(n::Int, m_::Int)
    polys = Vector{Vec2}[]
    for i in 0:(n - 1), j in 0:(m_ - 1)
        _kagome_polys!(polys, i, j)
    end
    g = mesh_from_polygons(polys)
    add_periodic_pairs!(g, Vec2(n, 0), Vec2(0.5 * m_, kSqrt3 * 0.5 * m_))
    return g
end

# Boundary-free (torus) patches: the n x m square / triangle grid with the wrap-around
# faces indexing the same vertices as their opposite side. The wrapping faces are
# geometrically degenerate, so these meshes serve ONLY the combinatorial checks.
function _torus_mesh(n::Int, m_::Int, P, tri::Bool)
    g = Mesh()
    V(i, j) = mod(i, n) * m_ + mod(j, m_) + 1   # 1-based
    for i in 0:(n - 1), j in 0:(m_ - 1)
        push!(g.X, P(i, j))
    end
    for i in 0:(n - 1), j in 0:(m_ - 1)
        if tri
            push!(g.faces, [V(i, j), V(i + 1, j), V(i, j + 1)])
            push!(g.faces, [V(i + 1, j), V(i + 1, j + 1), V(i, j + 1)])
        else
            push!(g.faces, [V(i, j), V(i + 1, j), V(i + 1, j + 1), V(i, j + 1)])
        end
    end
    build_topology!(g)
    return g
end
torus_squares(n::Int, m_::Int) = _torus_mesh(n, m_, (i, j) -> Vec2(i, j), false)
torus_triangles(n::Int, m_::Int) = _torus_mesh(n, m_, _tri_P, true)

# ---------------------------------------------------------------- random

mutable struct _Tri
    v::NTuple{3,Int}
    cc::Vec2
    r2::Float64
    alive::Bool
end

"""
    delaunay_triangles(pts) -> Vector{NTuple{3,Int}}

Delaunay triangulation of a point set (Bowyer-Watson with a super-triangle); returns
triangle index triples (1-based) in the C++ triangle order.
"""
function delaunay_triangles(pts::Vector{Vec2})
    n = length(pts)
    n < 3 && return NTuple{3,Int}[]
    minx = maxx = pts[1][1]
    miny = maxy = pts[1][2]
    for p in pts
        minx = min(minx, p[1]); maxx = max(maxx, p[1])
        miny = min(miny, p[2]); maxy = max(maxy, p[2])
    end
    d = max(maxx - minx, maxy - miny) * 20 + 10
    mid = Vec2((minx + maxx) / 2, (miny + maxy) / 2)
    P = copy(pts)
    push!(P, mid + Vec2(-d, -d))
    push!(P, mid + Vec2(d, -d))
    push!(P, mid + Vec2(0.0, d))
    a, b, c = n + 1, n + 2, n + 3

    # circumcircle; plain-double products/sums are fma-contracted as in the C++ build
    function circum(i, j, k)
        A = P[i]; B = P[j]; C = P[k]
        ax = B[1] - A[1]; ay = B[2] - A[2]
        bx = C[1] - A[1]; by = C[2] - A[2]
        det = 2 * fma(ax, by, -(ay * bx))
        abs(det) < 1e-18 && return (false, Vec2(0.0, 0.0), 0.0)
        a2 = fma(ax, ax, ay * ay)
        b2 = fma(bx, bx, by * by)
        cc = A + Vec2(fma(by, a2, -(ay * b2)) / det, fma(ax, b2, -(bx * a2)) / det)
        r2 = _sqnorm2(cc - A)
        return (true, cc, r2)
    end

    tris = _Tri[]
    let (_, cc, r2) = circum(a, b, c)
        push!(tris, _Tri((a, b, c), cc, r2, true))
    end
    edge_count = Dict{Tuple{Int,Int},Int}()
    for i in 1:n
        empty!(edge_count)
        for t in tris
            t.alive || continue
            if _sqnorm2(P[i] - t.cc) <= t.r2 * (1 + 1e-12)
                t.alive = false
                for e in 1:3
                    u = t.v[e]; w = t.v[mod1(e + 1, 3)]
                    key = (min(u, w), max(u, w))
                    edge_count[key] = get(edge_count, key, 0) + 1
                end
            end
        end
        # std::map iteration: keys in ascending (min, max) order
        for e in sort!(collect(keys(edge_count)))
            edge_count[e] != 1 && continue
            ok, cc, r2 = circum(e[1], e[2], i)
            ok && push!(tris, _Tri((e[1], e[2], i), cc, r2, true))
        end
        filter!(t -> t.alive, tris)
    end
    out = NTuple{3,Int}[]
    for t in tris
        t.alive || continue
        (t.v[1] > n || t.v[2] > n || t.v[3] > n) && continue
        push!(out, t.v)
    end
    return out
end

# n distinct points U(0,box)^2 (x drawn before y), rejecting points within 1e-4 of an
# earlier one (the rejected draw is consumed).
function random_points(n::Int, box::Float64, rng::MT19937)
    pts = Vec2[]
    w = VertexWelder(1e-4)
    while length(pts) < n
        x = uniform_real(rng, 0.0, box)
        y = uniform_real(rng, 0.0, box)
        p = Vec2(x, y)
        before = length(w.pts)
        add!(w, p) == before + 1 && push!(pts, p)
    end
    return pts
end

function delaunay_of_random_points(n::Int, box::Float64, rng::MT19937)
    pts = random_points(n, box, rng)
    tris = delaunay_triangles(pts)
    polys = [[pts[t[1]], pts[t[2]], pts[t[3]]] for t in tris]
    return largest_component(mesh_from_polygons(polys))
end

function voronoi_of_random_points(n::Int, box::Float64, rng::MT19937)
    pts = random_points(n, box, rng)
    polys = Vector{Vec2}[]
    # order sites by distance for early termination of the half-plane clipping
    order = collect(1:n)
    for i in 1:n
        cell = [Vec2(0.0, 0.0), Vec2(box, 0.0), Vec2(box, box), Vec2(0.0, box)]
        pi_ = pts[i]
        sort!(order; by = a -> _sqnorm2(pts[a] - pi_))
        for j in order
            j == i && continue
            dij = _norm2(pts[j] - pi_)
            rmax = 0.0
            for q in cell
                rmax = max(rmax, _norm2(q - pi_))
            end
            dij > 2 * rmax && break  # no further site can cut this cell
            # clip by the half-plane { x : |x - p_i| <= |x - p_j| }
            nvec = pts[j] - pi_
            off = _dot2(nvec, 0.5 * (pi_ + pts[j]))
            out = Vec2[]
            nc = length(cell)
            for k in 1:nc
                A = cell[k]
                B = cell[mod1(k + 1, nc)]
                da = _dot2(nvec, A) - off
                db = _dot2(nvec, B) - off
                da <= 0 && push!(out, A)
                if (da < 0 && db > 0) || (da > 0 && db < 0)
                    push!(out, A + (B - A) * (da / (da - db)))
                end
            end
            cell = out
            length(cell) < 3 && break
        end
        length(cell) >= 3 && push!(polys, cell)
    end
    g = mesh_from_polygons(polys, 1e-7)
    return largest_component(g)
end

# strict convexity of the quad (a, b, c, d) with the C++ contracted cross product
function _quad_convex(X::Vector{Vec2}, q::Vector{Int})
    convex = true
    for k in 1:4
        A = X[q[k]]; B = X[q[mod1(k + 1, 4)]]; C = X[q[mod1(k + 2, 4)]]
        u = B - A; v = C - B
        fma(u[1], v[2], -(u[2] * v[1])) <= 1e-12 && (convex = false)
    end
    return convex
end

"""Delaunay triangulation with pairs of triangles greedily merged across shuffled
interior edges whenever the union is convex."""
function quad_dominant_random(n::Int, box::Float64, rng::MT19937)
    tri = delaunay_of_random_points(n, box, rng)
    used = falses(n_faces(tri))
    eidx = collect(1:n_edges(tri))
    shuffle!(eidx, rng)
    polys = Vector{Vec2}[]
    for e in eidx
        ed = tri.edges[e]
        ed.n_faces != 2 && continue
        f1 = tri.half_edges[ed.he[1]].face
        f2 = tri.half_edges[ed.he[2]].face
        (used[f1] || used[f2]) && continue
        # build the quad: the shared edge's endpoints plus the two apexes
        a1 = a2 = 0
        for v in tri.faces[f1]
            (v != ed.key.a && v != ed.key.b) && (a1 = v)
        end
        for v in tri.faces[f2]
            (v != ed.key.a && v != ed.key.b) && (a2 = v)
        end
        quad = [ed.key.a, a1, ed.key.b, a2]
        if !_quad_convex(tri.X, quad)
            quad = [ed.key.a, a2, ed.key.b, a1]
            _quad_convex(tri.X, quad) || continue
        end
        used[f1] = used[f2] = true
        push!(polys, [tri.X[v] for v in quad])
    end
    for f in 1:n_faces(tri)
        used[f] && continue
        push!(polys, [tri.X[v] for v in tri.faces[f]])
    end
    return largest_component(mesh_from_polygons(polys))
end

"""
    generate(kind, params, rng) -> Mesh

Named factory (kiri_gen / sweep driver): `params[i]` defaults as in the C++
(`static_cast<int>` = truncation toward zero).
"""
function generate(kind::AbstractString, p::AbstractVector{<:Real}, rng::MT19937)
    par(i, d) = i <= length(p) ? Float64(p[i]) : d
    R = disk(Vec2(0.13, 0.07), par(1, 5.0))
    kind == "triangles" && return tiling_triangles(R)
    kind == "squares" && return tiling_squares(R)
    kind == "hexagons" && return tiling_hexagons(R)
    kind == "kagome" && return tiling_kagome(R)
    kind == "t3_4_3_12" && return tiling_3_4_3_12(R)
    kind == "snub_square" && return tiling_snub_square(R)
    kind == "truncated_square" && return tiling_truncated_square(R)
    kind == "squares_rect" &&
        return tiling_squares(rect(Vec2(par(1, 3.0) / 2.0, par(2, 3.0) / 2.0),
                                   par(1, 3.0) / 2.0 + 0.01, par(2, 3.0) / 2.0 + 0.01))
    kind == "periodic_squares" && return periodic_squares(trunc(Int, par(1, 3.0)), trunc(Int, par(2, 3.0)))
    kind == "periodic_triangles" && return periodic_triangles(trunc(Int, par(1, 3.0)), trunc(Int, par(2, 3.0)))
    kind == "periodic_hexagons" && return periodic_hexagons(trunc(Int, par(1, 3.0)), trunc(Int, par(2, 3.0)))
    kind == "torus_squares" && return torus_squares(trunc(Int, par(1, 4.0)), trunc(Int, par(2, 4.0)))
    kind == "torus_triangles" && return torus_triangles(trunc(Int, par(1, 4.0)), trunc(Int, par(2, 4.0)))
    kind == "periodic_kagome" && return periodic_kagome(trunc(Int, par(1, 3.0)), trunc(Int, par(2, 3.0)))
    kind == "delaunay" && return delaunay_of_random_points(trunc(Int, par(1, 50.0)), par(2, 10.0), rng)
    kind == "voronoi" && return voronoi_of_random_points(trunc(Int, par(1, 50.0)), par(2, 10.0), rng)
    kind == "quad_random" && return quad_dominant_random(trunc(Int, par(1, 50.0)), par(2, 10.0), rng)
    error("unknown generator kind: $kind")
end
