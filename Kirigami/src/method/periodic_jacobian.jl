# method/periodic_jacobian.jl -- the periodic (torus) quotient of a planar kirigami
# pattern, its Tutte auxetic shape space, and the deployment Jacobian of the
# fundamental parallelogram (2026 Sec. 5.1, Eqs. (11)-(13)).
#
# Port of code/src/method/periodic_jacobian.{hpp,cpp}. 1-based indices throughout;
# lattice offsets stay integer 2-vectors (`SVector{2,Int}`). C++ out-pointers become
# extra return values (`fk_period_matrix` returns `(P, spread)`).
#
# WHY A QUOTIENT AND NOT THE FINITE PATCH.  code/README.md deviation 11 imposes
# Eqs. (3b)-(3c) as linear constraints on a finite patch "whose boundary edges are
# not topologically identified".  The hole preimages that straddle the seam of that
# patch therefore contribute NO row of L (deviation 2 drops every preimage touching
# the mesh boundary), so the patch system is missing exactly the deployability
# conditions of the holes that wrap around.  Those conditions are what makes the
# face potential u of derivations/core.md T1 exist on the infinite tiling, which is
# what makes the deployed structure periodic at all.  This file therefore builds the
# genuine quotient: one face per translation class, edges keyed by (class pair,
# lattice offset), and hole rows carrying the integer offset sum on the right-hand
# side.
#
# THE RESULT (derived here, verified by kill_k7):
#   Let u be the face potential of T1 (u_g - u_f = sigma_g x_src(e) across a hinge
#   edge).  Translating the whole picture by a period t gives
#        u_{f+t} = u_f + w_t + sigma_f t / 2 ,        w_t constant over faces,
#   because the increment (sigma_g - sigma_f)/2 * t = sigma_g t absorbs the shift of
#   x_src.  Substituting into T1.5, y_{(v,f)} = cos(th/2) x_v + sin(th/2) Jrot (2 u_f
#   - sigma_f x_v), the sigma_f t terms cancel and
#        y_{(v+t, f+t)}(th) - y_{(v,f)}(th)  =  cos(th/2) t + 2 sin(th/2) Jrot w_t .
#   Hence with P_0 = [t_h t_v] and Q = 2 Jrot [w_h w_v],
#        P_theta = cos(theta/2) P_0 + sin(theta/2) Q ,
#        J(theta) = P_theta P_0^{-1} = cos(theta/2) I + sin(theta/2) K ,  K = Q P_0^{-1},
#   with K CONSTANT in theta and LINEAR in the vertex positions (P_0 is fixed on the
#   shape space because the periods are fixed data).  Two corollaries:
#     * J(theta) is a similarity for every theta  <=>  K is a similarity, because the
#       similarities are the linear span of {I, Jrot} and J(theta) is an affine path
#       through I.  This PROVES 2026 Sec. 5.1's "once the derivative of the conformal
#       distortion vanishes at theta = 0, the deployment remains conformal for all
#       opening angles", which the paper states without proof.
#     * the hole area per cell is det(P_0) (det J(theta) - 1), a first harmonic with
#       p + q = 0, whose second root is theta_c = 2 atan2(tr K, 1 - det K).

const Vec2i = SVector{2,Int}

"""The +90 degree rotation used throughout (same Jrot as deploy_basis)."""
rot90() = Mat2(0.0, 1.0, -1.0, 0.0)   # column-major: [0 -1; 1 0]

# ---------------------------------------------------------------- lattice

mutable struct Lattice
    t1::Vec2
    t2::Vec2
    ok::Bool
    err::String
end
Lattice() = Lattice(Vec2(1, 0), Vec2(0, 1), false, "")

_face_centroid(m::Mesh, f::Int) = sum(m.X[v] for v in m.faces[f]) / length(m.faces[f])

# std::llround: nearest integer, ties away from zero.
_llround(x::Float64) = round(Int64, x, RoundNearestTiesAway)

# side count + sorted edge lengths, quantised, as a face fingerprint.
function _face_signature(m::Mesh, f::Int, tol::Float64)
    F = m.faces[f]
    L = Int64[]
    for k in eachindex(F)
        d = norm(m.X[F[mod1(k + 1, length(F))]] - m.X[F[k]])
        push!(L, _llround(d / tol))
    end
    sort!(L)
    s = string(length(F))
    for v in L
        s *= "," * string(v)
    end
    return s
end

struct CentroidIndex
    bins::Dict{Tuple{Int64,Int64},Vector{Int}}
    h::Float64
    pts::Vector{Vec2}
end
function CentroidIndex(p::Vector{Vec2}, cell::Float64)
    bins = Dict{Tuple{Int64,Int64},Vector{Int}}()
    for i in eachindex(p)
        k = (floor(Int64, p[i][1] / cell), floor(Int64, p[i][2] / cell))
        push!(get!(bins, k, Int[]), i)
    end
    return CentroidIndex(bins, cell, p)
end
# Index of a point within `tol` of q, 0 if none (the C++ `find`).
function _find(idx::CentroidIndex, q::Vec2, tol::Float64)
    bx = floor(Int64, q[1] / idx.h)
    by = floor(Int64, q[2] / idx.h)
    for dx in -1:1, dy in -1:1
        lst = get(idx.bins, (bx + dx, by + dy), nothing)
        lst === nothing && continue
        for i in lst
            norm(idx.pts[i] - q) < tol && return i
        end
    end
    return 0
end

"""
    detect_lattice(big, tol = 1e-6) -> Lattice

Finds the two shortest independent translations that map the tiling to itself,
by matching face centroids of equal "signature" (side count + sorted edge lengths).
"""
function detect_lattice(big::Mesh, tol::Float64 = 1e-6)
    out = Lattice()
    NV = n_vertices(big)
    if n_faces(big) < 4 || NV < 8
        out.err = "too few faces"
        return out
    end
    # Vertex fingerprint: sorted sizes of the incident faces plus the sorted incident
    # edge lengths. A translation symmetry preserves it.
    vsig = Vector{String}(undef, NV)
    for v in 1:NV
        fs = sort!([length(big.faces[f]) for f in big.vertex_faces[v]])
        el = Int64[]
        for e in big.vertex_edges[v]
            k = big.edges[e].key
            push!(el, _llround(norm(big.X[k.a] - big.X[k.b]) / tol))
        end
        sort!(el)
        s = ""
        for a in fs
            s *= string(a) * "."
        end
        s *= "|"
        for a in el
            s *= string(a) * "."
        end
        vsig[v] = s
    end
    lo = big.X[1]
    hi = big.X[1]
    for p in big.X
        lo = min.(lo, p)
        hi = max.(hi, p)
    end
    ctr = 0.5 * (lo + hi)
    R = 0.0
    medge = 0.0
    for p in big.X
        R = max(R, norm(p - ctr))
    end
    for e in big.edges
        medge += norm(big.X[e.key.a] - big.X[e.key.b])
    end
    medge /= max(1, length(big.edges))
    margin = 4.0 * medge
    idx = CentroidIndex(big.X, max(1e-6, medge))

    v0 = 1
    for v in 2:NV
        norm(big.X[v] - ctr) < norm(big.X[v0] - ctr) && (v0 = v)
    end
    maxlen = max(3.0 * medge, 0.40 * R)

    good = Vec2[]
    for g in 1:NV
        (g == v0 || vsig[g] != vsig[v0]) && continue
        v = big.X[g] - big.X[v0]
        L = norm(v)
        (L < tol || L > maxlen) && continue
        tested = 0
        hit = 0
        for i in 1:NV
            qp = big.X[i] + v
            # both the source and its image must sit well inside the clipped patch, so
            # that a mismatch means a genuine failure of the symmetry and not the clip.
            (norm(qp - ctr) > R - margin || norm(big.X[i] - ctr) > R - margin) && continue
            tested += 1
            j = _find(idx, qp, 1e-6)
            (j > 0 && vsig[j] == vsig[i]) && (hit += 1)
        end
        (tested >= max(8, NV ÷ 8) && hit == tested) && push!(good, v)
    end
    if length(good) < 2
        out.err = "no lattice found ($(length(good)) candidates)"
        return out
    end
    # The C++ std::sort is unstable (libc++ introsort), so among candidates of equal
    # length the order it produced is not reproducible here; Julia's default sort is
    # stable and keeps vertex order among ties. data/corpus/method_fixtures freezes the
    # C++ lattice of every tiling family so the test can compare.
    sort!(good; by = v -> dot(v, v))
    t1 = good[1]
    t2 = Vec2(0, 0)
    found = false
    best_area = 0.0
    for v in good
        d = abs(t1[1] * v[2] - t1[2] * v[1])
        d < 1e-9 && continue
        if !found || d < best_area - 1e-9
            found = true
            best_area = d
            t2 = v
        end
    end
    if !found
        out.err = "all lattice candidates parallel"
        return out
    end
    for _ in 1:8
        k = round(dot(t2, t1) / dot(t1, t1))
        k == 0 && break
        t2 -= k * t1
    end
    out.t1 = t1
    out.t2 = t2
    # keep the basis right-handed, so det P_0 is the (positive) cell area
    if out.t1[1] * out.t2[2] - out.t1[2] * out.t2[1] < 0
        out.t1, out.t2 = out.t2, out.t1
    end
    out.ok = true
    return out
end

"""
    fundamental_domain(big, t1, t2, origin, tol = 1e-9) -> Mesh

The faces of `big` whose centroid has lattice coordinates in [0,1)^2 relative to
`origin`, extracted as a standalone Mesh (topology built, face order preserved).
"""
function fundamental_domain(big::Mesh, t1::Vec2, t2::Vec2, origin::Vec2, tol::Float64 = 1e-9)
    T = Mat2(t1[1], t1[2], t2[1], t2[2])
    Ti = inv(T)
    cell = Mesh()
    vmap = zeros(Int, n_vertices(big))
    for f in 1:n_faces(big)
        r = Ti * (_face_centroid(big, f) - origin)
        (floor(r[1] + 1e-7) != 0 || floor(r[2] + 1e-7) != 0) && continue
        face = Int[]
        for v in big.faces[f]
            if vmap[v] == 0
                push!(cell.X, big.X[v])
                vmap[v] = length(cell.X)
            end
            push!(face, vmap[v])
        end
        push!(cell.faces, face)
        isempty(big.sigma) || push!(cell.sigma, big.sigma[f])
    end
    build_topology!(cell)
    return cell
end

# ---------------------------------------------------------------- quotient

mutable struct QuotientEdge
    ca::Int                 # quotient vertex classes, canonical order
    cb::Int
    d::Vec2i                # lattice offset of cb relative to ca
    he::Vector{Int}         # cell half-edges realising the two sides
    f::Vector{Int}          # the two quotient faces
    hinge::Bool
    src_class::Int
    dst_class::Int
    src_off::Vec2i
    dst_off::Vec2i
end
QuotientEdge() = QuotientEdge(0, 0, Vec2i(0, 0), [0, 0], [0, 0], false, 0, 0, Vec2i(0, 0), Vec2i(0, 0))

mutable struct Quotient
    ok::Bool
    err::String
    cell::Mesh                          # one face per translation class
    T::Mat2                             # columns t_h, t_v
    nq::Int                             # quotient vertices
    vclass::Vector{Int}                 # per cell vertex -> class
    voff::Vector{Vec2i}                 # per cell vertex -> lattice offset
    class_rep::Vector{Int}              # per class -> a cell vertex
    Xq::Vector{Vec2}                    # per class -> representative position
    edges::Vector{QuotientEdge}
    n_hinge::Int
    n_split::Int
    H::Int
    hole_edges::Vector{Vector{Int}}     # hinge quotient-edge indices per hole
    L::Matrix{Float64}                  # H x nq
    Noff::Matrix{Float64}               # H x 2, the integer lattice-offset sums
end
Quotient() = Quotient(false, "", Mesh(), Mat2(1.0, 0.0, 0.0, 1.0), 0, Int[], Vec2i[], Int[],
                      Vec2[], QuotientEdge[], 0, 0, 0, Vector{Int}[], zeros(0, 0), zeros(0, 2))
n_faces(q::Quotient) = n_faces(q.cell)

"""
    build_quotient(cell, T, tol = 1e-6) -> Quotient

`cell` must carry sigma and be a fundamental domain for T (one face per class).
"""
function build_quotient(cell_in::Mesh, T::Mat2, tol::Float64 = 1e-6)
    q = Quotient()
    q.cell = deepcopy(cell_in)   # the C++ copies the mesh by value
    q.T = T
    cell = q.cell
    isempty(cell.half_edges) && build_topology!(cell)
    if length(cell.sigma) != n_faces(cell)
        q.err = "sigma missing"
        return q
    end
    Ti = inv(T)
    NV = n_vertices(cell)
    q.vclass = fill(0, NV)
    q.voff = fill(Vec2i(0, 0), NV)

    # Class = position modulo the lattice. Reduce, then match canonical points.
    canon = Vec2[]
    for i in 1:NV
        r = Ti * cell.X[i]
        o = Vec2i(floor(Int, r[1] + tol), floor(Int, r[2] + tol))
        p = cell.X[i] - T * o
        q.voff[i] = o
        cls = 0
        for c in eachindex(canon)
            if norm(canon[c] - p) < 1e-6
                cls = c
                break
            end
        end
        if cls == 0
            push!(canon, p)
            push!(q.class_rep, i)
            cls = length(canon)
        end
        q.vclass[i] = cls
    end
    q.nq = length(canon)
    q.Xq = canon

    # Quotient edges, keyed by (class pair, lattice offset), canonicalised. The C++
    # std::map iterates in key order, so the edge list is built from the sorted keys.
    bykey = Dict{NTuple{4,Int},Vector{Int}}()
    for h in eachindex(cell.half_edges)
        he = cell.half_edges[h]
        ca = q.vclass[he.from]
        cb = q.vclass[he.to]
        d = q.voff[he.to] - q.voff[he.from]
        if ca < cb || (ca == cb && (d[1] > 0 || (d[1] == 0 && d[2] > 0)))
            k = (ca, cb, d[1], d[2])
        else
            k = (cb, ca, -d[1], -d[2])
        end
        if ca == cb && d[1] == 0 && d[2] == 0
            q.err = "degenerate self-loop edge"
            return q
        end
        push!(get!(bykey, k, Int[]), h)
    end
    n_unpaired = 0
    for k in sort!(collect(keys(bykey)))
        hs = bykey[k]
        if length(hs) == 1
            n_unpaired += 1
            continue
        end
        if length(hs) != 2
            q.err = "quotient edge with $(length(hs)) sides"
            return q
        end
        e = QuotientEdge()
        e.ca = k[1]
        e.cb = k[2]
        e.d = Vec2i(k[3], k[4])
        e.he = [hs[1], hs[2]]
        e.f = [cell.half_edges[hs[1]].face, cell.half_edges[hs[2]].face]
        if e.f[1] == e.f[2]
            q.err = "quotient edge with both sides on one face"
            return q
        end
        e.hinge = cell.sigma[e.f[1]] != cell.sigma[e.f[2]]
        dir = half_edge_direction(cell, e.he[1])
        e.src_class = q.vclass[dir[1]]
        e.dst_class = q.vclass[dir[2]]
        e.src_off = q.voff[dir[1]]
        e.dst_off = q.voff[dir[2]]
        if e.hinge
            dir2 = half_edge_direction(cell, e.he[2])
            rel1 = e.dst_off - e.src_off
            rel2 = q.voff[dir2[2]] - q.voff[dir2[1]]
            if q.vclass[dir2[1]] != e.src_class || q.vclass[dir2[2]] != e.dst_class || rel1 != rel2
                q.err = "hinge direction disagrees across the two sides"
                return q
            end
            q.n_hinge += 1
        else
            q.n_split += 1
        end
        push!(q.edges, e)
    end
    if n_unpaired > 0
        q.err = "not a fundamental domain: $n_unpaired unpaired edges"
        return q
    end

    # Hole preimages: components of the split subgraph on the quotient vertices.
    dsu = DSU(q.nq)
    for e in q.edges
        e.hinge || (dsu.p[find!(dsu, e.ca)] = find!(dsu, e.cb))
    end
    comp = zeros(Int, q.nq)
    for v in 1:q.nq
        r = find!(dsu, v)
        if comp[r] == 0
            q.H += 1
            comp[r] = q.H
        end
        comp[v] = comp[r]
    end
    for v in 1:q.nq
        comp[v] = comp[find!(dsu, v)]
    end

    q.L = zeros(q.H, q.nq)
    q.Noff = zeros(q.H, 2)
    q.hole_edges = [Int[] for _ in 1:q.H]
    for i in eachindex(q.edges)
        e = q.edges[i]
        e.hinge || continue
        k = comp[e.dst_class]
        push!(q.hole_edges[k], i)
        q.L[k, e.dst_class] += 1.0
        q.L[k, e.src_class] -= 1.0
        q.Noff[k, 1] += e.dst_off[1] - e.src_off[1]
        q.Noff[k, 2] += e.dst_off[2] - e.src_off[2]
    end
    q.ok = true
    return q
end

"""
    quotient_system(q, Xq_ini) -> LinearSystem

[L; e_pin] X = [-(T n_k)^T ; x_pin], the periodic Tutte auxetic system on the
quotient. Rows H+1, columns nq.
"""
function quotient_system(q::Quotient, Xq_ini::Vector{Vec2})
    sys = LinearSystem()
    sys.N = q.nq
    sys.n_hole_rows = q.H
    sys.n_boundary_rows = 1
    sys.A = zeros(q.H + 1, q.nq)
    sys.rhs = zeros(q.H + 1, 2)
    sys.A[1:q.H, :] = q.L
    for k in 1:q.H
        n = Vec2(q.Noff[k, 1], q.Noff[k, 2])
        r = -(q.T * n)
        sys.rhs[k, 1] = r[1]
        sys.rhs[k, 2] = r[2]
    end
    sys.A[q.H + 1, 1] = 1.0
    sys.rhs[q.H + 1, 1] = Xq_ini[1][1]
    sys.rhs[q.H + 1, 2] = Xq_ini[1][2]
    return sys
end

"""Residual of Eq. (2) per hole for quotient positions Xq (max norm)."""
function quotient_residual(q::Quotient, Xq::Vector{Vec2})
    worst = 0.0
    for k in 1:q.H
        s = q.T * Vec2(q.Noff[k, 1], q.Noff[k, 2])
        for v in 1:q.nq
            q.L[k, v] != 0.0 && (s += q.L[k, v] * Xq[v])
        end
        worst = max(worst, norm(s))
    end
    return worst
end

# ---------------------------------------------------------------- super patch

# (2*half+1)^2 copies of the fundamental domain, used as the forward-kinematics
# object: the translation vectors of the deployed tiling are read off it directly.
mutable struct SuperPatch
    mesh::Mesh
    half::Int
    nfc::Int                          # cell faces
    vert_class::Vector{Int}           # per super vertex
    vert_off::Vector{Vec2i}
    face_cell::Vector{Int}            # per super face -> cell face
    face_off::Vector{Vec2i}
end
SuperPatch() = SuperPatch(Mesh(), 1, 0, Int[], Vec2i[], Int[], Vec2i[])
# 1-based: f is a 1-based cell face, the result a 1-based super face.
face_index(sp::SuperPatch, i::Int, j::Int, f::Int) =
    ((i + sp.half) * (2 * sp.half + 1) + (j + sp.half)) * sp.nfc + f

function build_super(q::Quotient, half::Int = 1)
    sp = SuperPatch()
    sp.half = half
    sp.nfc = n_faces(q.cell)
    R = 2 * half + 1
    vid = Dict{Tuple{Int,Int,Int},Int}()
    function getid(cls::Int, o::Vec2i)
        key = (cls, o[1], o[2])
        id = get(vid, key, 0)
        id != 0 && return id
        push!(sp.vert_class, cls)
        push!(sp.vert_off, o)
        id = length(sp.vert_class)
        vid[key] = id
        return id
    end
    nsf = R * R * sp.nfc
    sp.mesh.faces = [Int[] for _ in 1:nsf]
    sp.mesh.sigma = zeros(Int, nsf)
    sp.face_cell = zeros(Int, nsf)
    sp.face_off = fill(Vec2i(0, 0), nsf)
    for i in -half:half, j in -half:half, f in 1:sp.nfc
        fi = face_index(sp, i, j, f)
        face = Int[]
        for v in q.cell.faces[f]
            push!(face, getid(q.vclass[v], q.voff[v] + Vec2i(i, j)))
        end
        sp.mesh.faces[fi] = face
        sp.mesh.sigma[fi] = q.cell.sigma[f]
        sp.face_cell[fi] = f
        sp.face_off[fi] = Vec2i(i, j)
    end
    sp.mesh.X = fill(Vec2(0, 0), length(sp.vert_class))
    set_super_positions!(sp, q, q.Xq)
    build_topology!(sp.mesh)
    return sp
end

"""Writes the super-patch vertex positions from quotient positions Xq (C++ `set_super_positions`)."""
function set_super_positions!(sp::SuperPatch, q::Quotient, Xq::Vector{Vec2})
    for i in eachindex(sp.vert_class)
        sp.mesh.X[i] = Xq[sp.vert_class[i]] + q.T * sp.vert_off[i]
    end
    return sp
end

# ---------------------------------------------------------------- Jacobian

mutable struct PeriodicJac
    ok::Bool
    P0::Mat2
    Q::Mat2
    K::Mat2
    consistency::Float64   # max spread of the measured period over (face, corner)
    p0_err::Float64        # max |measured P_0 column - the lattice column|
end
PeriodicJac() = PeriodicJac(false, Mat2(1.0, 0.0, 0.0, 1.0), zeros(Mat2), zeros(Mat2), 0.0, 0.0)

"""
    periodic_jacobian(sp, q, c, B) -> PeriodicJac

P_theta = cos(th/2) P0 + sin(th/2) Q read off the deploy basis of the super patch.
"""
# `B` is a DeployBasis (deploy_basis.jl is included after this file, so no annotation).
function periodic_jacobian(sp::SuperPatch, q::Quotient, c::CutStructure, B)
    out = PeriodicJac()
    sp.half < 1 && return out
    n = 0
    cs = (Vec2[], Vec2[])
    ss = (Vec2[], Vec2[])
    Bc(i) = Vec2(B.C[i, 1], B.C[i, 2])
    Bs(i) = Vec2(B.S[i, 1], B.S[i, 2])
    for f in 1:sp.nfc
        nk = length(q.cell.faces[f])
        for k in 1:nk
            a = c.corner_to_prime[c.face_corner_base[face_index(sp, 0, 0, f)] + k]
            bh = c.corner_to_prime[c.face_corner_base[face_index(sp, 1, 0, f)] + k]
            bv = c.corner_to_prime[c.face_corner_base[face_index(sp, 0, 1, f)] + k]
            push!(cs[1], Bc(bh) - Bc(a))
            push!(ss[1], Bs(bh) - Bs(a))
            push!(cs[2], Bc(bv) - Bc(a))
            push!(ss[2], Bs(bv) - Bs(a))
            n += 1
        end
    end
    n == 0 && return out
    Csum = zeros(2, 2)
    Ssum = zeros(2, 2)
    for d in 1:2
        mc = Vec2(0, 0)
        ms = Vec2(0, 0)
        for i in 1:n
            mc += cs[d][i]
            ms += ss[d][i]
        end
        mc /= n
        ms /= n
        for i in 1:n
            out.consistency = max(out.consistency, norm(cs[d][i] - mc))
            out.consistency = max(out.consistency, norm(ss[d][i] - ms))
        end
        Csum[:, d] = mc
        Ssum[:, d] = ms
    end
    out.P0 = Mat2(Csum)
    out.Q = Mat2(Ssum)
    out.p0_err = maximum(abs, Csum - q.T)
    abs(det(out.P0)) < 1e-12 && return out
    out.K = out.Q * inv(out.P0)
    out.ok = true
    return out
end

"""J(theta) from K."""
jac_at(K::Mat2, theta::Float64) = cos(0.5 * theta) * Mat2(1.0, 0.0, 0.0, 1.0) + sin(0.5 * theta) * K

"""Relative conformal distortion of a 2x2 matrix: (smax - smin) / (smax + smin)."""
function conformal_distortion(A::Mat2)
    s = svdvals(A)
    s0, s1 = s[1], s[2]
    s0 + s1 <= 0 && return 0.0
    return (s0 - s1) / (s0 + s1)
end

"""
    poisson_ratio(K, theta, d)

Secant Poisson ratio along the reference direction d at opening angle theta:
  nu = -(|J d_perp| - 1) / (|J d| - 1)      (d, d_perp unit, J = J(theta)).
Rotation invariant, hence a function of the polar stretch U of J alone.
"""
function poisson_ratio(K::Mat2, theta::Float64, d::Vec2)
    Jt = jac_at(K, theta)
    u = d / norm(d)
    up = Vec2(-u[2], u[1])
    lp = norm(Jt * u) - 1.0
    lt = norm(Jt * up) - 1.0
    return -lt / lp
end

"""
    quotient_sigma(nf, dual, rng) -> Vector{Int}

Max-cut sigma on the QUOTIENT dual graph (2026 Eq. (1)'s purpose, solved directly
on the torus adjacency rather than on a patch's boundary-truncated dual).  Returns the
bipartite 2-colouring when the dual is bipartite and Gamma is connected.
"""
function quotient_sigma(nf::Int, dual::Vector{Tuple{Int,Int}}, rng::MT19937)
    nf <= 0 && return Int[]
    adj = [Int[] for _ in 1:nf]
    for (a, b) in dual
        push!(adj[a], b)
        push!(adj[b], a)
    end
    cut(s) = count(((a, b),) -> s[a] != s[b], dual)
    function gamma_connected(s)
        seen = falses(nf)
        st = Int[1]
        seen[1] = true
        n = 1
        while !isempty(st)
            f = pop!(st)
            for g in adj[f]
                if !seen[g] && s[g] != s[f]
                    seen[g] = true
                    n += 1
                    push!(st, g)
                end
            end
        end
        return n == nf
    end
    # bipartite?
    col = zeros(Int, nf)
    bip = true
    for s in 1:nf
        bip || break
        col[s] != 0 && continue
        col[s] = -1
        st = Int[s]
        while !isempty(st) && bip
            f = pop!(st)
            for g in adj[f]
                if col[g] == 0
                    col[g] = -col[f]
                    push!(st, g)
                elseif col[g] == col[f]
                    bip = false
                end
            end
        end
    end
    bip && gamma_connected(col) && return col
    best = Int[]
    best_cut = -1
    for _ in 1:40
        s = [uniform_int(rng, 0, 1) != 0 ? 1 : -1 for _ in 1:nf]
        improved = true
        while improved
            improved = false
            for f in 1:nf
                same = 0
                diff = 0
                for g in adj[f]
                    s[g] == s[f] ? (same += 1) : (diff += 1)
                end
                if same > diff
                    s[f] = -s[f]
                    improved = true
                end
            end
        end
        c = cut(s) + (gamma_connected(s) ? 1000000 : 0)
        if c > best_cut
            best_cut = c
            best = s
        end
    end
    return best
end

"""The quotient dual graph as face-index pairs, one per quotient edge."""
quotient_dual(q::Quotient) = Tuple{Int,Int}[(e.f[1], e.f[2]) for e in q.edges]

"""
    fk_period_matrix(sp, q, cut, theta) -> (P, spread)

P_theta measured by an INDEPENDENT forward-kinematics call: deploy() the super
patch to opening angle theta and average, over every (cell face, corner), the
displacement between the copy at lattice offset (0,0) and the copy at (1,0) resp.
(0,1).  `spread` receives the largest deviation of any single such displacement
from that average, together with the deployment mismatch -- it is the numerical
witness that the deployed tiling really is a translate of itself.
"""
function fk_period_matrix(sp::SuperPatch, q::Quotient, cut::CutStructure, theta::Float64)
    D = deploy(cut, sp.mesh.X, theta)
    M = zeros(2, 2)
    sp2 = 0.0
    for d in 1:2
        v = Vec2[]
        for f in 1:sp.nfc, kk in 1:length(q.cell.faces[f])
            a = cut.corner_to_prime[cut.face_corner_base[face_index(sp, 0, 0, f)] + kk]
            b = cut.corner_to_prime[cut.face_corner_base[face_index(sp, d == 1 ? 1 : 0, d == 1 ? 0 : 1, f)] + kk]
            push!(v, D.Y[b] - D.Y[a])
        end
        m = sum(v) / length(v)
        for p in v
            sp2 = max(sp2, norm(p - m))
        end
        M[:, d] = m
    end
    return Mat2(M), max(sp2, D.max_mismatch)
end

"""
    cell_face_area_sum(sp)

Sum of |signed area| over the faces of one cell of `sp` at its current positions.
Constant along the deployment, so the hole area per cell is |det P_theta| minus it.
"""
function cell_face_area_sum(sp::SuperPatch)
    a = 0.0
    for f in 1:sp.nfc
        a += abs(face_signed_area(sp.mesh, face_index(sp, 0, 0, f)))
    end
    return a
end

# ---------------------------------------------------------------- periodic patterns

# A periodic kirigami pattern: one fundamental domain `cell` (sigma assigned) plus
# the period matrix T = [t_h t_v].  Shared by the K7 experiment and its tests so both
# exercise the same object.
mutable struct PeriodicPattern
    name::String
    family::String
    cell::Mesh
    T::Mat2
    ok::Bool
    err::String
end
PeriodicPattern() = PeriodicPattern("", "", Mesh(), Mat2(1.0, 0.0, 0.0, 1.0), false, "")

"""
    make_tiling_pattern(family, n, m, rng) -> PeriodicPattern

An n x m super-cell of one of the periodic tiling families ("squares",
"triangles", "hexagons", "kagome", "snub_square", "trunc_square_488",
"t3_4_3_12"), with sigma from quotient_sigma.  `rng` seeds the max-cut restarts.
"""
function make_tiling_pattern(family::AbstractString, n::Int, m::Int, rng::MT19937)
    P = PeriodicPattern()
    P.family = family
    P.name = family * "_" * string(n) * "x" * string(m)
    local big::Mesh
    try
        if family == "squares"
            big = tiling_squares(rect(Vec2(0, 0), 6.5, 6.5))
        elseif family == "triangles"
            big = tiling_triangles(disk(Vec2(0.05, 0.03), 7.0))
        elseif family == "hexagons"
            big = tiling_hexagons(disk(Vec2(0.05, 0.03), 9.0))
        elseif family == "kagome"
            big = tiling_kagome(disk(Vec2(0.05, 0.03), 7.0))
        elseif family == "snub_square"
            big = tiling_snub_square(disk(Vec2(0.05, 0.03), 9.0))
        elseif family == "trunc_square_488"
            big = tiling_truncated_square(disk(Vec2(0.05, 0.03), 13.0))
        elseif family == "t3_4_3_12"
            big = tiling_3_4_3_12(disk(Vec2(0.05, 0.03), 15.0))
        else
            P.err = "unknown family"
            return P
        end
    catch e
        P.err = "generate: " * sprint(showerror, e)
        return P
    end
    lat = detect_lattice(big)
    if !lat.ok
        P.err = "lattice: " * lat.err
        return P
    end
    t1 = Float64(n) * lat.t1
    t2 = Float64(m) * lat.t2
    # origin: the centroid of the face nearest the patch centre
    lo = big.X[1]
    hi = big.X[1]
    for p in big.X
        lo = min.(lo, p)
        hi = max.(hi, p)
    end
    mid = 0.5 * (lo + hi)
    f0 = 1
    bd = 1e300
    for f in 1:n_faces(big)
        c = _face_centroid(big, f)
        if norm(c - mid) < bd
            bd = norm(c - mid)
            f0 = f
        end
    end
    org = _face_centroid(big, f0)
    # shift the origin off the face centroid lattice so no centroid sits on the cut
    local cell::Mesh
    try
        cell = fundamental_domain(big, t1, t2, org - 0.001 * (t1 + t2))
    catch e
        P.err = "fundamental_domain: " * sprint(showerror, e)
        return P
    end
    if n_faces(cell) < 2
        P.err = "cell has $(n_faces(cell)) faces"
        return P
    end
    T = Mat2(t1[1], t1[2], t2[1], t2[2])
    cell.sigma = fill(-1, n_faces(cell))
    q0 = build_quotient(cell, T)
    if !q0.ok
        P.err = "quotient(probe): " * q0.err
        return P
    end
    cell.sigma = quotient_sigma(n_faces(cell), quotient_dual(q0), rng)
    P.cell = cell
    P.T = T
    P.ok = true
    return P
end

# ---------------------------------------------------------------- achievable set

# The affine map t (in R^{2k}, k = dim null) -> K, as K0 + sum_j t_j M_j.
mutable struct AchievableSet
    K0::Mat2
    M::Vector{Mat2}          # 2k generators
    A::Matrix{Float64}       # 4 x 2k, columns = vec(M_j)
    dimK::Int                # rank(A)
    D::Matrix{Float64}       # k x 2 period-potential increments
    rankD::Int
end
AchievableSet() = AchievableSet(zeros(Mat2), Mat2[], zeros(4, 0), 0, zeros(0, 2), 0)
