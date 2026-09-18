# export/solid.jl -- turning a Layout into a closed, consistently oriented triangle
# soup: every face is extruded by the sheet thickness, and
#
#   * living-hinge-neck: all faces share one z slab, and the wall that two faces
#     would both build over their common neck bridge is dropped, so the whole
#     sheet comes out as ONE closed manifold solid;
#   * pin-pad: each face is its own closed solid with a through-hole at every
#     hinge site, and the two faces of a hinge sit in different z slabs
#     (sigma 2-colours the hinge adjacency), so they overlap only over the pad.
#
# Triangle indices are 1-based into `V`.

const Vec3 = SVector{3,Float64}

mutable struct TriMesh
    V::Vector{Vec3}
    T::Vector{NTuple{3,Int}}
    # Optional: which solid body each triangle belongs to. Pin-pad plates are
    # separate bodies that may overlap in plan view (they are stacked in z), so
    # they must be welded and checked per body, not globally.
    tri_group::Vector{Int}
end
TriMesh() = TriMesh(Vec3[], NTuple{3,Int}[], Int[])
n_tris(M::TriMesh) = length(M.T)

# Every `x*y - z*w` is evaluated as fma(x, y, -(z*w)) (docs/NUMERICS.md). The
# ear-clipping tie-breaks and the STL float32 normals depend on those sub-ulp
# residuals, so the 2D/3D cross products spell the contraction out explicitly
# (verified byte-for-byte against the archived hero exports).
fms(x, y, z, w) = fma(x, y, -(z * w))  # x*y - z*w, contracted
cross_fma(a::Vec3, b::Vec3) = Vec3(fms(a[2], b[3], a[3], b[2]),
                                   fms(a[3], b[1], a[1], b[3]),
                                   fms(a[1], b[2], a[2], b[1]))

function normal(M::TriMesh, t::Int)
    a = M.V[M.T[t][1]]
    b = M.V[M.T[t][2]]
    c = M.V[M.T[t][3]]
    return cross_fma(b - a, c - a)
end

"""Axis-aligned bounding box `(lo, hi)` of the vertices (zeros when empty)."""
function bbox(M::TriMesh)
    isempty(M.V) && return (Vec3(0, 0, 0), Vec3(0, 0, 0))
    lo = hi = M.V[1]
    for p in M.V
        lo = min.(lo, p)
        hi = max.(hi, p)
    end
    return lo, hi
end

mutable struct ManifoldReport
    closed::Bool                 # every undirected edge used by exactly 2 triangles
    consistently_oriented::Bool  # every directed edge used exactly once
    n_boundary_edges::Int        # edges used once
    n_nonmanifold_edges::Int     # edges used more than twice
    n_flipped_edges::Int         # directed edges used more than once
    n_degenerate::Int            # zero-area triangles
    n_components::Int
    volume::Float64              # signed volume via the divergence theorem (mm^3)
end
ManifoldReport() = ManifoldReport(false, false, 0, 0, 0, 0, 0, 0.0)
ok(r::ManifoldReport) = r.closed && r.consistently_oriented && r.n_degenerate == 0

function summary(r::ManifoldReport)
    return (r.closed ? "closed" : "NOT closed") * ", " *
           (r.consistently_oriented ? "consistently oriented" : "INCONSISTENT orientation") *
           ", boundary edges $(r.n_boundary_edges), non-manifold edges $(r.n_nonmanifold_edges)" *
           ", flipped edges $(r.n_flipped_edges), degenerate tris $(r.n_degenerate)" *
           ", components $(r.n_components), volume $(fmt_num(r.volume)) mm^3"
end

cross2(a::Vec2, b::Vec2) = fms(a[1], b[2], a[2], b[1])

function point_in_tri(p::Vec2, a::Vec2, b::Vec2, c::Vec2, eps::Float64)
    d1 = cross2(b - a, p - a)
    d2 = cross2(c - b, p - b)
    d3 = cross2(a - c, p - c)
    return d1 > eps && d2 > eps && d3 > eps
end

# Quantized 3D key, so vertices that should coincide actually do (std::llround
# rounds ties away from zero).
const Key3 = NTuple{3,Int64}
llround(x::Float64) = round(Int64, x, RoundNearestTiesAway)
key_of(p::Vec3, tol::Float64) = (llround(p[1] / tol), llround(p[2] / tol), llround(p[3] / tol))

# deliberate: `join` does not path-compress the second argument (the component
# representatives, and with them the 3MF object order, depend on it)
function dsu_join!(d::DSU, a::Int, b::Int)
    d.p[find!(d, a)] = find!(d, b)
end

# per-body split of a grouped mesh: one TriMesh per tri_group value, in order of
# first appearance (std::map<int, TriMesh> in check_manifold iterates by key; the
# aggregate there is order-independent)
function split_by_group(M::TriMesh)
    parts = Dict{Int,TriMesh}()
    for t in eachindex(M.T)
        p = get!(TriMesh, parts, M.tri_group[t])
        tri = ntuple(k -> (push!(p.V, M.V[M.T[t][k]]); length(p.V)), 3)
        push!(p.T, tri)
    end
    return parts
end

"""Manifold / orientation / volume report of a triangle soup (welded on `weld_tol`)."""
function check_manifold(M::TriMesh, weld_tol::Float64 = 1e-6)
    if !isempty(M.tri_group)
        # Analyse each body on its own, then aggregate.
        parts = split_by_group(M)
        all_ = ManifoldReport()
        all_.closed = true
        all_.consistently_oriented = true
        for k in sort!(collect(keys(parts)))
            r = check_manifold(parts[k], weld_tol)
            all_.closed = all_.closed && r.closed
            all_.consistently_oriented = all_.consistently_oriented && r.consistently_oriented
            all_.n_boundary_edges += r.n_boundary_edges
            all_.n_nonmanifold_edges += r.n_nonmanifold_edges
            all_.n_flipped_edges += r.n_flipped_edges
            all_.n_degenerate += r.n_degenerate
            all_.n_components += r.n_components
            all_.volume += r.volume
        end
        return all_
    end
    r = ManifoldReport()
    # Weld vertices first: STL round-trips lose the index structure.
    # Grid hash with a 27-cell probe, so two vertices that differ by less than the
    # tolerance still weld even when they straddle a cell boundary.
    weld = Dict{Key3,Int}()
    id = zeros(Int, length(M.V))
    nv = 0
    for i in eachindex(M.V)
        k = key_of(M.V[i], weld_tol)
        found = 0
        for dx in -1:1, dy in -1:1, dz in -1:1
            found != 0 && break
            f = get(weld, (k[1] + dx, k[2] + dy, k[3] + dz), 0)
            f != 0 && (found = f)
        end
        if found == 0
            nv += 1
            found = nv
        end
        haskey(weld, k) || (weld[k] = found)
        id[i] = found
    end
    dsu = DSU(max(nv, 1))
    directed = Dict{Tuple{Int,Int},Int}()
    undirected = Dict{Tuple{Int,Int},Int}()
    for (ti, t) in enumerate(M.T)
        a = id[t[1]]
        b = id[t[2]]
        c = id[t[3]]
        if a == b || b == c || a == c || norm(normal(M, ti)) < 1e-12
            r.n_degenerate += 1
        end
        dsu_join!(dsu, a, b)
        dsu_join!(dsu, b, c)
        for e in ((a, b), (b, c), (c, a))
            directed[e] = get(directed, e, 0) + 1
            ue = (min(e[1], e[2]), max(e[1], e[2]))
            undirected[ue] = get(undirected, ue, 0) + 1
        end
    end
    for v in values(undirected)
        v == 1 && (r.n_boundary_edges += 1)
        v > 2 && (r.n_nonmanifold_edges += 1)
    end
    for v in values(directed)
        v > 1 && (r.n_flipped_edges += 1)
    end
    r.closed = (r.n_boundary_edges == 0 && r.n_nonmanifold_edges == 0)
    r.consistently_oriented = (r.n_flipped_edges == 0)
    seen = falses(max(nv, 1))
    for t in M.T
        seen[find!(dsu, id[t[1]])] = true
    end
    r.n_components = 0
    for i in 1:nv
        (seen[i] && find!(dsu, i) == i) && (r.n_components += 1)
    end
    vol = 0.0
    for t in M.T
        a = M.V[t[1]]
        b = M.V[t[2]]
        c = M.V[t[3]]
        vol += dot(a, cross_fma(b, c))
    end
    r.volume = vol / 6.0
    return r
end

"""
    triangulate_with_holes(outer, holes) -> (clean, poly, tris)

Ear-clipping triangulation of a simple polygon with holes (outer CCW, holes CW).
`tris` are index triples into `poly`, the concatenated vertex list (outer, then
each hole spliced in with a bridge). `clean` is false when the clipping stalled
and the remainder was fanned.
"""
function triangulate_with_holes(outer::Vector{Vec2}, holes::Vector{Vector{Vec2}})
    tris = NTuple{3,Int}[]
    length(outer) < 3 && return false, Vec2[], tris

    # Work on an index list into `pts`; holes are spliced in with a bridge.
    pts = copy(outer)
    loop = collect(1:length(outer))

    # Sort holes by descending max-x so the bridges never cross a later hole.
    hs = [h for h in holes if length(h) >= 3]
    sort!(hs; by = h -> -maximum(p -> p[1], h), alg = MergeSort)

    for h in hs
        mi = 1
        for i in 2:length(h)
            h[i][1] > h[mi][1] && (mi = i)
        end
        Mp = h[mi]

        # Ray M -> +x: find the loop edge it first hits.
        best_t = 1e30
        best_edge = 0
        I = Vec2(0, 0)
        n = length(loop)
        for i in 1:n
            a = pts[loop[i]]
            b = pts[loop[mod1(i + 1, n)]]
            (a[2] > Mp[2]) == (b[2] > Mp[2]) && continue
            s = (Mp[2] - a[2]) / (b[2] - a[2])
            x = a[1] + s * (b[1] - a[1])
            if x >= Mp[1] && x - Mp[1] < best_t
                best_t = x - Mp[1]
                best_edge = i
                I = Vec2(x, Mp[2])
            end
        end
        best_edge == 0 && return false, Vec2[], tris
        # Bridge vertex: the endpoint of that edge with the larger x, then any
        # reflex loop vertex inside triangle (M, I, P) that is more "visible".
        P = loop[best_edge]
        if pts[loop[mod1(best_edge + 1, n)]][1] > pts[P][1]
            P = loop[mod1(best_edge + 1, n)]
        end
        let
            bestP = P
            best_ang = 1e30
            for i in 1:n
                vi = loop[i]
                vi == P && continue
                p = pts[vi]
                p[1] < Mp[1] && continue
                (!point_in_tri(p, Mp, I, pts[P], -1e-12) && vi != P) && continue
                ang = libm_atan2(abs(p[2] - Mp[2]), p[1] - Mp[1])
                if ang < best_ang
                    best_ang = ang
                    bestP = vi
                end
            end
            P = bestP
        end
        # Splice: ... P, M, hole..., M, P ...
        base = length(pts)
        append!(pts, h)
        push!(pts, Mp)
        dupM = length(pts)
        push!(pts, pts[P])
        dupP = length(pts)

        merged = Int[]
        hn = length(h)
        at = findfirst(==(P), loop)
        at === nothing && return false, Vec2[], tris
        for i in 1:at
            push!(merged, loop[i])
        end
        for k in 0:hn-1
            push!(merged, base + mod(mi - 1 + k, hn) + 1)
        end
        push!(merged, dupM)
        push!(merged, dupP)
        for i in at+1:n
            push!(merged, loop[i])
        end
        loop = merged
    end

    # ---- ear clipping --------------------------------------------------------
    scale = 0.0
    for p in pts
        scale = max(scale, norm(p))
    end
    eps = max(1e-12, 1e-12 * scale * scale)
    idx = copy(loop)
    guard = 4 * length(idx) * length(idx) + 64
    while length(idx) > 3 && guard > 0
        guard -= 1
        n = length(idx)
        clipped = false
        for i in 1:n
            ia = idx[mod1(i - 1, n)]
            ib = idx[i]
            ic = idx[mod1(i + 1, n)]
            a = pts[ia]
            b = pts[ib]
            c = pts[ic]
            cross2(b - a, c - a) <= eps && continue  # reflex or collinear
            ok_ = true
            for j in 1:n
                ok_ || break
                iv = idx[j]
                (iv == ia || iv == ib || iv == ic) && continue
                point_in_tri(pts[iv], a, b, c, eps) && (ok_ = false)
            end
            ok_ || continue
            push!(tris, (ia, ib, ic))
            deleteat!(idx, i)
            clipped = true
            break
        end
        clipped && continue
        # No ear. Clipping valid ears can leave a zero-width needle (two boundary
        # edges doubling back along each other), which happens wherever a neck tab
        # is attached to a face that is also inset along the adjacent side. Clip the
        # needle tip: the triangle is degenerate, but the boundary edges still
        # cancel, so the extruded solid stays closed.
        needle = 0
        best_flat = 1e30
        for i in 1:n
            a = pts[idx[mod1(i - 1, n)]]
            b = pts[idx[i]]
            c = pts[idx[mod1(i + 1, n)]]
            cr = abs(cross2(b - a, c - a))
            if cr <= eps && dot(a - b, c - b) > 0 && cr < best_flat
                best_flat = cr
                needle = i
            end
        end
        if needle > 0
            push!(tris, (idx[mod1(needle - 1, n)], idx[needle], idx[mod1(needle + 1, n)]))
            deleteat!(idx, needle)
            continue
        end
        # Last resort: clip the most convex corner regardless of containment. The
        # triangulation may then self-overlap, but the solid is still closed.
        best = 0
        best_cr = eps
        for i in 1:n
            a = pts[idx[mod1(i - 1, n)]]
            b = pts[idx[i]]
            c = pts[idx[mod1(i + 1, n)]]
            cr = cross2(b - a, c - a)
            if cr > best_cr
                best_cr = cr
                best = i
            end
        end
        best == 0 && break
        push!(tris, (idx[mod1(best - 1, n)], idx[best], idx[mod1(best + 1, n)]))
        deleteat!(idx, best)
    end
    clean = (length(idx) == 3)
    if clean
        push!(tris, (idx[1], idx[2], idx[3]))
    elseif length(idx) > 3
        # Terminal fallback: fan the remainder. The triangles may be degenerate or
        # overlap, but every boundary edge is still used exactly once, so the
        # extruded solid stays closed and consistently oriented.
        for i in 2:length(idx)-1
            push!(tris, (idx[1], idx[i], idx[i + 1]))
        end
    end
    return clean, pts, tris
end

struct SolidLoop
    pts::Vector{Vec2}
    z0::Float64
    z1::Float64
end

"""
    build_solid(L::Layout) -> (TriMesh, warnings)

Builds the 3D solid described at the top of this file.
"""
function build_solid(L::Layout)
    warnings = String[]
    M = TriMesh()
    R = pin_hole_radius(L.profile)
    pad_mode = (L.profile.hinge == PinPad)

    # Every wall segment we would build, so the pair over a neck bridge cancels.
    # The two faces of a hinge compute the tab tip from their own copies of the
    # edge, so the coordinates agree only to rounding: weld all outline points to
    # canonical ids first, or the two walls fail to cancel and the edge ends up
    # carried by four triangles.
    tol = 1e-6
    weld = Dict{Key3,Int}()
    n_weld = 0
    function weld_id(p::Vec2, z::Float64)
        k = key_of(Vec3(p[1], p[2], z), tol)
        for dx in -1:1, dy in -1:1, dz in -1:1
            f = get(weld, (k[1] + dx, k[2] + dy, k[3] + dz), 0)
            if f != 0
                haskey(weld, k) || (weld[k] = f)
                return f
            end
        end
        n_weld += 1
        weld[k] = n_weld
        return n_weld
    end
    seg_count = Dict{Tuple{Int,Int},Int}()

    nf = length(L.pieces)
    face_loops = [SolidLoop[] for _ in 1:nf]
    face_poly = [Vec2[] for _ in 1:nf]
    face_tris = [NTuple{3,Int}[] for _ in 1:nf]

    for fi in 1:nf
        fp = L.pieces[fi]
        length(fp.outline) < 3 && continue
        holes = Vector{Vec2}[]
        if pad_mode
            for ctr in fp.pin_holes
                circle = Vec2[]
                nseg = 24
                for k in 0:nseg-1  # clockwise: a hole
                    a = -2 * pi * k / nseg
                    push!(circle, ctr + R * Vec2(libm_cos(a), libm_sin(a)))
                end
                push!(holes, circle)
            end
        end
        clean, poly, tris = triangulate_with_holes(fp.outline, holes)
        face_poly[fi] = poly
        face_tris[fi] = tris
        if !clean
            push!(warnings, "face $(fp.face - 1): ear clipping stalled; the remaining ring was closed with a fan (possibly degenerate triangles, but still a closed solid)")
        end
        push!(face_loops[fi], SolidLoop(fp.outline, fp.z0, fp.z1))
        for h in holes
            push!(face_loops[fi], SolidLoop(h, fp.z0, fp.z1))
        end
        for lp in face_loops[fi]
            n = length(lp.pts)
            for i in 1:n
                key = (weld_id(lp.pts[i], lp.z0), weld_id(lp.pts[mod1(i + 1, n)], lp.z0))
                seg_count[key] = get(seg_count, key, 0) + 1
            end
        end
    end

    # tri_group values are 0-based face ids so the 3MF object order and the per-body
    # checks match the archived exports. Living-hinge sheets are one body (0).
    tri_group_of(face) = pad_mode ? face - 1 : 0
    function add_v(p::Vec2, z::Float64)
        push!(M.V, Vec3(p[1], p[2], z))
        return length(M.V)
    end

    for fi in 1:nf
        fp = L.pieces[fi]
        isempty(face_poly[fi]) && continue
        z0 = fp.z0
        z1 = fp.z1
        # caps
        np = length(face_poly[fi])
        bot = zeros(Int, np)
        top = zeros(Int, np)
        for i in 1:np
            bot[i] = add_v(face_poly[fi][i], z0)
            top[i] = add_v(face_poly[fi][i], z1)
        end
        grp = tri_group_of(fp.face)
        for t in face_tris[fi]
            push!(M.T, (bot[t[1]], bot[t[3]], bot[t[2]]))  # downward normal
            push!(M.T, (top[t[1]], top[t[2]], top[t[3]]))  # upward normal
            push!(M.tri_group, grp)
            push!(M.tri_group, grp)
        end
        # walls, skipping a segment whose reverse another face also builds
        for lp in face_loops[fi]
            n = length(lp.pts)
            for i in 1:n
                a = lp.pts[i]
                b = lp.pts[mod1(i + 1, n)]
                if get(seg_count, (weld_id(b, lp.z0), weld_id(a, lp.z0)), 0) > 0
                    continue  # shared neck bridge
                end
                a0 = add_v(a, z0)
                b0 = add_v(b, z0)
                a1 = add_v(a, z1)
                b1 = add_v(b, z1)
                push!(M.T, (a0, b0, b1))
                push!(M.T, (a0, b1, a1))
                push!(M.tri_group, grp)
                push!(M.tri_group, grp)
            end
        end
    end
    return M, warnings
end
