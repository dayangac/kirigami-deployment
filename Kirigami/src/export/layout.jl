# export/layout.jl -- the shared 2D stage of every writer: deploy the kirigami at an
# angle, scale it to millimetres, and offset every face away from the edges that
# are actually cut, keeping a neck of material at each hinge vertex.
#
# Both the SVG writer and the 3D writers consume `Layout`; nothing downstream
# touches `Mesh` or `CutStructure` directly, so the fabrication geometry is
# defined in exactly one place.
#
# Geometry decisions:
#
#  * A face's cut outline is its deployed polygon offset inward by
#    `face_inset(profile)` along every side that is actually cut (hinge or
#    split) and by 0 along border sides, so the sheet perimeter is the true
#    mesh boundary and two faces sharing a cut edge end up `gap` apart once the
#    kerf has burned away.
#
#  * Living-hinge neck. The hinge of edge e = (src -> dst) lives at `src`, which
#    is a *corner* of both incident faces, so a neck "centred" on it can only be
#    a bridge straddling e. Each incident face therefore keeps a rectangular tab
#    of material against e, anchored at `src` and running `neck_width` along e.
#    The union of the two tabs is one bridge of length `neck_width` centred on
#    the cut line, located at the hinge vertex. It is clipped to 0.45*|e| so the
#    tab can never swallow the edge.
#
#  * Pin-pad. The tab is replaced by a circular lobe of radius `pad_radius`
#    centred on the hinge vertex, plus a through-hole of radius
#    `pin_hole_radius()`. The lobe deliberately overhangs the face outline: the
#    two faces of a hinge overlap exactly over the lobe and are stacked in z by
#    sigma (hinge edges are precisely the dual edges with opposite sigma, so
#    sigma is a proper 2-colouring of the hinge adjacency), which is what makes
#    "alternating up/down, overlapping only at the pad" well defined.
#
# Indices are 1-based; "none" markers are 0.

# One face of M', ready to cut or extrude. All coordinates are millimetres.
mutable struct FacePiece
    face::Int
    raw::Vector{Vec2}      # deployed face polygon, CCW, un-inset
    outline::Vector{Vec2}  # the cut outline: `raw` inset by the profile, with
                           # a tab of material at every hinge source vertex
    # per side i (raw[i] -> raw[i+1]) of the face
    side_type::Vector{EdgeType}
    side_edge::Vector{Int}         # edge index in the original mesh (0 = none)
    side_neck_start::Vector{Bool}  # hinge side whose source is corner i
    side_neck_end::Vector{Bool}    # hinge side whose source is corner i+1
    # Indices into `outline` of the two endpoints of each neck tab, or 0.
    side_tab_begin::Vector{Int}
    side_tab_end::Vector{Int}
    # Pin-pad mode only: centres of the through-holes drilled at the hinge sites
    # of this face (radius pin_hole_radius(profile)).
    pin_holes::Vector{Vec2}
    # Corners of the face that carry a neck tab / pad lobe (a hinge source).
    corner_is_hinge_source::Vector{Bool}
    # z-range of the extruded plate (mm); hinge-adjacent faces alternate layers.
    z0::Float64
    z1::Float64
end
FacePiece() = FacePiece(0, Vec2[], Vec2[], EdgeType[], Int[], Bool[], Bool[], Int[], Int[],
                        Vec2[], Bool[], 0.0, 0.0)

# A place where two faces stay connected: the source vertex of a hinge edge.
mutable struct HingeSite
    edge::Int          # edge index in the original mesh
    face_a::Int
    face_b::Int
    src_vertex::Int
    dst_vertex::Int
    p::Vec2            # the hinge point, in mm (shared by both faces)
    dir_a::Vec2        # unit direction src->dst inside face_a, in mm space
    dir_b::Vec2        # same inside face_b (equal to dir_a only when theta == 0)
    neck::Float64      # neck length actually used on this edge (mm)
    setback::Float64   # distance from the hinge vertex at which the neck starts (mm)
end
HingeSite() = HingeSite(0, 0, 0, 0, 0, Vec2(0, 0), Vec2(0, 0), Vec2(0, 0), 0.0, 0.0)

mutable struct Layout
    mesh::Mesh
    cut::CutStructure
    profile::MaterialProfile
    theta::Float64
    scale::Float64  # mm per input unit

    Y::Vector{Vec2}  # M' vertices in mm at this theta
    pieces::Vector{FacePiece}
    hinges::Vector{HingeSite}

    bb_min::Vec2  # bounding box of the raw (un-inset) geometry
    bb_max::Vec2
    warnings::Vector{String}
end

n_faces(L::Layout) = length(L.pieces)
width(L::Layout) = L.bb_max[1] - L.bb_min[1]
height(L::Layout) = L.bb_max[2] - L.bb_min[2]

const kLayoutEps = 1e-9

# libm_* shims live in core/mesh.jl (shared with generators).

# rot90 comes from core/collision.jl

# Intersection of the lines p0 + s*d0 and p1 + s*d1. Falls back to the midpoint
# of the two base points when the directions are (numerically) parallel.
function line_intersect(p0::Vec2, d0::Vec2, p1::Vec2, d1::Vec2)
    den = d0[1] * d1[2] - d0[2] * d1[1]
    abs(den) < 1e-12 && return 0.5 * (p0 + p1)
    w = p1 - p0
    s = (w[1] * d1[2] - w[2] * d1[1]) / den
    return p0 + s * d0
end

# The intersection of circle(c, R) with the line p + s*d that lies furthest along
# `pref`. Returns `nothing` when the line misses the circle.
function circle_line(c::Vec2, R::Float64, p::Vec2, d::Vec2, pref::Vec2)
    f = p - c
    a = dot(d, d)
    a < 1e-18 && return nothing
    b = 2.0 * dot(f, d)
    cc = dot(f, f) - R * R
    disc = b * b - 4 * a * cc
    disc < 0 && return nothing
    sq = sqrt(disc)
    h0 = p + ((-b - sq) / (2 * a)) * d
    h1 = p + ((-b + sq) / (2 * a)) * d
    return dot(h0 - c, pref) >= dot(h1 - c, pref) ? h0 : h1
end

function push_unique!(poly::Vector{Vec2}, p::Vec2)
    (!isempty(poly) && norm(poly[end] - p) < 5e-6) && return  # matches the 1e-6 mm weld
    push!(poly, p)
end

"""Signed area of a polygon (> 0 iff counter-clockwise)."""
function polygon_area(p::Vector{Vec2})
    a = 0.0
    n = length(p)
    for i in 1:n
        u = p[i]
        v = p[mod1(i + 1, n)]
        a += u[1] * v[2] - v[1] * u[2]
    end
    return 0.5 * a
end

"""
    build_layout(m, c, X, theta, scale, profile) -> Layout

`X` is the flat embedding in input units (typically mesh.X or an optimized X0).
"""
function build_layout(m::Mesh, c::CutStructure, X::Vector{Vec2}, theta::Float64,
                      scale::Float64, profile::MaterialProfile)
    warnings = String[]
    for w in check(profile)
        push!(warnings, "profile: " * w)
    end

    dep = deploy(c, X, theta)
    Y = [scale * y for y in dep.Y]

    d0 = face_inset(profile)
    pad_mode = (profile.hinge == PinPad)
    R = profile.pad_radius

    # ---- per-face pieces -----------------------------------------------------
    pieces = [FacePiece() for _ in 1:n_faces(m)]
    bb_min = Vec2(0, 0)
    bb_max = Vec2(0, 0)
    bb_init = false
    for f in 1:n_faces(m)
        fp = pieces[f]
        fp.face = f
        fv = m.faces[f]
        pv = c.prime_faces[f]
        n = length(fv)

        fp.raw = [Y[pv[i]] for i in 1:n]
        for q in fp.raw
            if !bb_init
                bb_min = bb_max = q
                bb_init = true
            else
                bb_min = min.(bb_min, q)
                bb_max = max.(bb_max, q)
            end
        end

        fp.side_type = fill(Border, n)
        fp.side_edge = zeros(Int, n)
        fp.side_neck_start = falses(n)
        fp.side_neck_end = falses(n)
        fp.side_tab_begin = zeros(Int, n)
        fp.side_tab_end = zeros(Int, n)
        fp.corner_is_hinge_source = falses(n)

        # side i runs raw[i] -> raw[i+1]
        u = Vector{Vec2}(undef, n)
        nrm = Vector{Vec2}(undef, n)
        len = zeros(n)
        ins = zeros(n)
        for i in 1:n
            j = mod1(i + 1, n)
            e = get(m.edge_index, EdgeKey(fv[i], fv[j]), 0)
            fp.side_edge[i] = e
            fp.side_type[i] = (e == 0) ? Border : c.edge_type[e]
            dvec = fp.raw[j] - fp.raw[i]
            len[i] = norm(dvec)
            u[i] = (len[i] > kLayoutEps) ? dvec / len[i] : Vec2(1, 0)
            nrm[i] = rot90(u[i])  # interior side of a CCW polygon
            ins[i] = (fp.side_type[i] == Border) ? 0.0 : d0
            if fp.side_type[i] == Hinge && e > 0
                hd = c.hinge_dir[e]
                hd.src == fv[i] && (fp.side_neck_start[i] = true)
                hd.src == fv[j] && (fp.side_neck_end[i] = true)
            end
            if len[i] > kLayoutEps && ins[i] > 0.4 * len[i]
                # face / side ids are reported 0-based (the JSON file's numbering)
                push!(warnings, "face $(f - 1) side $(i - 1): inset $(fmt_num(ins[i])) mm exceeds 40% of the side length $(fmt_num(len[i])) mm")
            end
        end

        # corner i is a hinge source iff one of its two sides is a hinge edge whose
        # source vertex is fv[i].
        for i in 1:n
            pi_ = mod1(i - 1, n)
            hp = fp.side_neck_end[pi_]  # prev side ends at a hinge source
            hn = fp.side_neck_start[i]
            fp.corner_is_hinge_source[i] = hp || hn
        end

        # ---- outline -----------------------------------------------------------
        # Corners are the plain inset corner. A hinge side carries a rectangular
        # TAB against the uncut edge, set back by `setback` from the hinge (source)
        # vertex and `neck_width` long. Both incident faces build the same tab, so
        # the union is one bridge of material straddling the cut. The setback is
        # what keeps two independent hinges that share a vertex from meeting at a
        # point -- a pinch that is neither manufacturable nor 2-manifold.
        empty!(fp.outline)
        for i in 1:n
            pi_ = mod1(i - 1, n)
            p = fp.raw[i]
            bp = p + ins[pi_] * nrm[pi_]
            bn = p + ins[i] * nrm[i]

            if pad_mode && fp.corner_is_hinge_source[i]
                A = circle_line(p, R, bp, u[pi_], -u[pi_])
                B = circle_line(p, R, bn, u[i], u[i])
                if A !== nothing && B !== nothing
                    a0 = libm_atan2((A - p)[2], (A - p)[1])
                    a1 = libm_atan2((B - p)[2], (B - p)[1])
                    # The union with the disc bulges OUTWARD, so the boundary sweeps the
                    # long way round p: counter-clockwise for a CCW outline.
                    sweep = a1 - a0
                    while sweep <= 1e-9
                        sweep += 2pi
                    end
                    steps = max(6, Int(ceil(sweep / 0.12)))
                    for st in 0:steps
                        ang = a0 + sweep * (st / steps)
                        push_unique!(fp.outline, p + R * Vec2(libm_cos(ang), libm_sin(ang)))
                    end
                    push!(fp.pin_holes, p)
                    continue
                end
                push!(warnings, "face $(f - 1) corner $(i - 1): pad radius does not reach the offset outline; no pad emitted")
            end
            push_unique!(fp.outline, line_intersect(bp, u[pi_], bn, u[i]))

            # tab on side i, if that side is a hinge
            (pad_mode || fp.side_type[i] != Hinge) && continue
            (!fp.side_neck_start[i] && !fp.side_neck_end[i]) && continue
            s_back = min(max(ins[i], 1e-3), 0.2 * len[i])
            Ln = min(profile.neck_width, len[i] - 2 * s_back - 1e-6)
            if Ln <= 0
                push!(warnings, "face $(f - 1) side $(i - 1): edge too short for a neck; the hinge is left uncut-free")
                continue
            end
            t0 = fp.side_neck_start[i] ? s_back : (len[i] - s_back - Ln)
            t1 = t0 + Ln
            T0 = fp.raw[i] + t0 * u[i]
            T1 = fp.raw[i] + t1 * u[i]
            push_unique!(fp.outline, T0 + ins[i] * nrm[i])
            push_unique!(fp.outline, T0)
            fp.side_tab_begin[i] = length(fp.outline)
            push_unique!(fp.outline, T1)
            fp.side_tab_end[i] = length(fp.outline)
            push_unique!(fp.outline, T1 + ins[i] * nrm[i])
        end
        if length(fp.outline) > 1 && norm(fp.outline[1] - fp.outline[end]) < 1e-7
            pop!(fp.outline)
        end
        # Drop straight-through vertices. They are harmless in 2D but make the
        # extruded caps carry zero-area triangles, and removing them here keeps the
        # caps and the side walls using exactly the same ring.
        changed = true
        while changed && length(fp.outline) > 3
            changed = false
            no = length(fp.outline)
            for k in 1:no
                a = fp.outline[mod1(k - 1, no)]
                b = fp.outline[k]
                cc = fp.outline[mod1(k + 1, no)]
                e0 = b - a
                e1 = cc - b
                n0 = norm(e0)
                n1 = norm(e1)
                if n0 < 1e-9 || n1 < 1e-9 ||
                   (abs(e0[1] * e1[2] - e0[2] * e1[1]) < 1e-9 * n0 * n1 && dot(e0, e1) > 0)
                    deleteat!(fp.outline, k)
                    changed = true
                    break
                end
            end
        end
        # Index bookkeeping for the tabs is invalidated by the cleanup above; recompute.
        for k in 1:n
            fp.side_tab_begin[k] == 0 && continue
            e = fp.side_edge[k]
            hd = c.hinge_dir[e]
            src_corner = (hd.src == fv[k]) ? k : mod1(k + 1, n)
            sb = min(max(ins[k], 1e-3), 0.2 * len[k])
            Ln = min(profile.neck_width, len[k] - 2 * sb - 1e-6)
            dir = (src_corner == k) ? u[k] : -u[k]
            T0 = fp.raw[src_corner] + sb * dir
            T1 = fp.raw[src_corner] + (sb + Ln) * dir
            fp.side_tab_begin[k] = fp.side_tab_end[k] = 0
            for q in eachindex(fp.outline)
                norm(fp.outline[q] - T0) < 1e-7 && (fp.side_tab_begin[k] = q)
                norm(fp.outline[q] - T1) < 1e-7 && (fp.side_tab_end[k] = q)
            end
            if fp.side_tab_begin[k] > fp.side_tab_end[k]
                fp.side_tab_begin[k], fp.side_tab_end[k] = fp.side_tab_end[k], fp.side_tab_begin[k]
            end
        end
        if polygon_area(fp.outline) <= 0
            push!(warnings, "face $(f - 1): outline is degenerate or inverted (area $(fmt_num(polygon_area(fp.outline))) mm^2) -- inset/neck too large for this face")
        end

        # z layer: sigma is a proper 2-colouring of the hinge adjacency.
        # The two layers are separated by `clearance` so the stacked plates are two
        # disjoint solids that can actually rotate against each other.
        sg = (f <= length(m.sigma)) ? m.sigma[f] : -1
        fp.z0 = (pad_mode && sg > 0) ? (profile.thickness + profile.clearance) : 0.0
        fp.z1 = fp.z0 + profile.thickness
    end

    L = Layout(m, c, profile, theta, scale, Y, pieces, HingeSite[], bb_min, bb_max, warnings)

    # ---- hinge sites ---------------------------------------------------------
    for e in c.hinge_edges
        ed = m.edges[e]
        ed.n_faces != 2 && continue
        hs = HingeSite()
        hs.edge = e
        hs.src_vertex = c.hinge_dir[e].src
        hs.dst_vertex = c.hinge_dir[e].dst
        hs.face_a = m.half_edges[ed.he[1]].face
        hs.face_b = m.half_edges[ed.he[2]].face
        pa = prime_vertex(c, hs.face_a, hs.src_vertex)
        pb = prime_vertex(c, hs.face_b, hs.src_vertex)
        hs.p = Y[pa]
        if norm(Y[pa] - Y[pb]) > 1e-6 * max(1.0, width(L))
            push!(warnings, "hinge edge $(e - 1): the two copies of the source vertex are $(fmt_num(norm(Y[pa] - Y[pb]))) mm apart")
        end
        da = Y[prime_vertex(c, hs.face_a, hs.dst_vertex)] - Y[pa]
        db = Y[prime_vertex(c, hs.face_b, hs.dst_vertex)] - Y[pb]
        la = norm(da)
        lb = norm(db)
        hs.dir_a = (la > kLayoutEps) ? da / la : Vec2(1, 0)
        hs.dir_b = (lb > kLayoutEps) ? db / lb : Vec2(1, 0)
        le = min(la, lb)
        hs.setback = min(max(face_inset(profile), 1e-3), 0.2 * le)
        hs.neck = min(profile.neck_width, le - 2 * hs.setback - 1e-6)
        # Design rule for pin-pad: the pad lobe is a full disc about the hinge
        # vertex, so it reaches into every face meeting there. Faces in the SAME z
        # layer (same sigma) would then interfere. Report it rather than silently
        # shipping a colliding assembly.
        if pad_mode && hs.src_vertex <= length(m.vertex_faces)
            for g in m.vertex_faces[hs.src_vertex]
                (g == hs.face_a || g == hs.face_b) && continue
                sg = (g <= length(m.sigma)) ? m.sigma[g] : -1
                sa = (hs.face_a <= length(m.sigma)) ? m.sigma[hs.face_a] : -1
                sb = (hs.face_b <= length(m.sigma)) ? m.sigma[hs.face_b] : -1
                if sg == sa || sg == sb
                    push!(warnings, "pin-pad: the pad of hinge edge $(e - 1) at vertex $(hs.src_vertex - 1) (radius $(fmt_num(profile.pad_radius)) mm) reaches face $(g - 1), which prints in the same layer -- reduce pad_radius or increase the scale")
                    break
                end
            end
        end
        push!(L.hinges, hs)
    end
    return L
end

# Default number formatting of the warning texts (%g, precision 6).
fmt_num(v::Float64) = @sprintf("%g", v)
