# check_b4_93.jl -- Checker round 6, add-on: the B4 `voronoi_93` disagreement.
# Port of check_b4_93.cpp (same graph, same sigma_def, same repair seeds, same printed lines).
#
# results/kill/b4/b4.csv, row (id 93, voronoi, sigma_def, free):
#     cert_pos=1  cert_noovl=1  cert_noroot=0  eps_max=0.2484
#     theta_exact=0.2484   theta_bisect=0
# while (id 96, ..., free) has theta_exact == theta_bisect == 0.241884.
#
# This program replays both rows through kill_b4's own pipeline (same graph, same
# sigma_def from results/kill/k5/sigma, same `free` variant system, same repair seeds and
# weights) and then asks WHICH predicate separates the two answers:
#
#   1. the BROAD PHASE.  exact_theta_max_overlap is run on candidate_pairs(..., prune=true)
#      while the referee's has_collision() tests EVERY face pair.  Re-run the scan with
#      prune=false: if theta_exact collapses, a pair the grid dropped is the culprit.
#   2. the SHRINK.  the scan is called with shrink = 1e-9, the referee with 1e-12.
#   3. the theta = 0+ PROBE.  referee_theta returns 0 as soon as has_collision fires at
#      theta = 1e-7; the scan's interval sweep starts at theta_lo = 1e-9 and decides the
#      first slab by an overlap probe at its MIDPOINT.
#
# Run (from repo root):
#   julia --project=Kirigami derivations/scratch/check_b4_93.jl [id ...]     (default: 93 96)
#
# CAVEAT (port): the C++ read the `free`-variant shape (X0, Phi) from the B4 cache
# results/kill/b4/cache/free_sigma_def/shape_<id>.bin, which was NOT migrated (caches are
# regenerable intermediates).  Phi is only determined up to an orthogonal transform inside the
# null space and the repair is NOT invariant under one, so a recomputed shape does not replay
# b4.csv (id 93 comes out infeasible, id 96 at a different Theta_max).  Set the environment
# variable KIRI_B4_CACHE to a directory holding the C++ `shape_<id>.bin` files (e.g. the source
# repo's results/kill/b4/cache/free_sigma_def) to replay the archived rows exactly; without it
# the shape is recomputed and the printout is a fresh design, not the B4 row.
#
# Depends on Kirigami/src/method/{contact,zero_plus,deploy_basis,design}.jl (C++ names).

using Kirigami
using LinearAlgebra
using Printf
using StaticArrays
const K = Kirigami
const Vec2 = K.Vec2
include(joinpath(@__DIR__, "corpus_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

Bc(B, i::Int) = Vec2(B.C[i, 1], B.C[i, 2])
Bs(B, i::Int) = Vec2(B.S[i, 1], B.S[i, 2])
heval(h, th::Float64) = h.p + h.q * cos(th) + h.r * sin(th)
hamp(h) = hypot(h.q, h.r)
cross2(a::Vec2, b::Vec2) = a[1] * b[2] - a[2] * b[1]

# kill_b4's referee_theta, verbatim.
function referee_theta(c::K.CutStructure, X::Vector{Vec2}, grid::Int = 4000, iters::Int = 50)
    col(th) = K.has_collision(c, K.deploy(c, X, th).Y, 1e-12)
    col(1e-7) && return 0.0
    lo = 0.0; hi = -1.0
    for i in 1:grid
        th = pi * i / grid
        if col(th)
            hi = th
            break
        end
        lo = th
    end
    hi < 0 && return Float64(pi)
    for i in 1:iters
        mid = 0.5 * (lo + hi)
        if col(mid)
            hi = mid
        else
            lo = mid
        end
    end
    return lo
end

# kill_common's Shape (the part this program reads).
mutable struct Shape
    N::Int
    k::Int
    X0::Matrix{Float64}   # N x 2
    Phi::Matrix{Float64}  # N x k
    ok::Bool
end

# kill_common's load_shape: int32[5] header {N, k, rank_L, H, n_interior}, then X0 (N x 2)
# and Phi (N x k) as column-major little-endian doubles.
function load_shape(dir::String, id::Int)
    path = joinpath(dir, "shape_" * string(id) * ".bin")
    isfile(path) || return nothing
    open(path, "r") do io
        hdr = [read(io, Int32) for _ in 1:5]
        N = Int(hdr[1]); k = Int(hdr[2])
        X0 = Matrix{Float64}(undef, N, 2); read!(io, X0)
        Phi = Matrix{Float64}(undef, N, k)
        k > 0 && read!(io, Phi)
        return Shape(N, k, X0, Phi, true)
    end
end

# kill_b4's `free` variant: drop every Eq. (4) boundary row, pin vertex 1.
# (The C++ consulted the B4 cache first; see the caveat in the header.)
function shape_free(m::K.Mesh, c::K.CutStructure, hs::K.HoleSet, id::Int)
    cache = get(ENV, "KIRI_B4_CACHE", "")
    if !isempty(cache)
        s = load_shape(cache, id)
        if s !== nothing && s.N == K.n_vertices(m)
            @printf("  shape (X0, Phi) loaded from %s\n", cache)
            return s
        end
    end
    sys = K.assemble_system(c, hs, m.X, K.None)
    nh = size(sys.A, 1)
    N = sys.N
    A = zeros(nh + 1, N); rhs = zeros(nh + 1, 2)
    A[1:nh, :] = sys.A
    rhs[1:nh, :] = sys.rhs
    A[nh + 1, 1] = 1.0
    rhs[nh + 1, 1] = m.X[1][1]
    rhs[nh + 1, 2] = m.X[1][2]
    sys.A = A
    sys.rhs = rhs
    sys.n_boundary_rows = 1
    sr = K.solve_system(sys, m.X)
    return Shape(N, sr.dim_null, sr.X0, sr.Phi, sr.projection_ok)
end

# proper (open) segment crossing test, for the simplicity screen
function seg_cross(p1::Vec2, p2::Vec2, q1::Vec2, q2::Vec2)
    d1 = cross2(p2 - p1, q1 - p1); d2 = cross2(p2 - p1, q2 - p1)
    d3 = cross2(q2 - q1, p1 - q1); d4 = cross2(q2 - q1, p2 - q1)
    return ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0))
end

# which face pair collides at theta, over ALL pairs, at the given shrink (0-based ids
# like the C++ printout, (-1,-1) if none)
function colliding_pair(c::K.CutStructure, Y::Vector{Vec2}, shrink::Float64)
    F = K.n_faces(c.mesh)
    for f in 1:F
        A = [Y[pv] for pv in c.prime_faces[f]]
        for g in f+1:F
            B = [Y[pv] for pv in c.prime_faces[g]]
            K.polygons_overlap(A, B, shrink) && return (f - 1, g - 1)
        end
    end
    return (-1, -1)
end

poly_of(c::K.CutStructure, Y::Vector{Vec2}, f::Int) = [Y[pv] for pv in c.prime_faces[f]]

function pip(P::Vector{Vec2}, z::Vec2)
    inside = false
    n = length(P)
    j = n
    for i in 1:n
        if ((P[i][2] > z[2]) != (P[j][2] > z[2])) &&
           (z[1] < (P[j][1] - P[i][1]) * (z[2] - P[i][2]) / (P[j][2] - P[i][2]) + P[i][1])
            inside = !inside
        end
        j = i
    end
    return inside
end

function edge_dist(P::Vector{Vec2}, z::Vec2)
    d = 1e300
    n = length(P)
    for i in 1:n
        a2 = P[i]; b2 = P[i % n + 1]
        e2 = b2 - a2
        t2 = max(0.0, min(1.0, dot(z - a2, e2) / max(dot(e2, e2), 1e-300)))
        d = min(d, norm(z - (a2 + t2 * e2)))
    end
    return d
end

function replay(id::Int)
    @printf("================ id %d ================\n", id)
    g = make_graph(id, 100, 800, 1400)
    if !g.ok
        @printf("  graph not ok\n"); return
    end
    m0 = g.mesh
    K.build_topology!(m0)
    p = joinpath(REPO, "results", "kill", "k5", "sigma", g.kind * "_" * string(id) * ".json")
    if !isfile(p)
        @printf("  no sigma_def at %s\n", p); return
    end
    sm = K.load_mesh_json(p)
    m = m0
    m.sigma = sm.sigma
    K.build_topology!(m)
    c = K.make_cut(m)
    hs = K.holes_partition(c)
    sh = shape_free(m, c, hs, id)
    if !sh.ok || sh.k < 1
        @printf("  shape not ok\n"); return
    end
    X0 = K.matrix_to_points(sh.X0)
    med = K.median_edge_length(m)
    seed0 = 6000 + 7 * id + 1   # which = 1 (sigma_def)

    ro = K.ZeroPlusRepairOptions()
    ro.n_random = 3; ro.max_iter = 1200; ro.lambda_rel = 1e-6; ro.seed = seed0
    sr = K.zero_plus_repair(c, X0, sh.Phi, med, ro)
    ro2 = deepcopy(ro)
    ro2.w_corner = 0.05; ro2.w_prox = 1e3
    sr.feasible && (ro2.t_init = sr.t)
    rr = K.zero_plus_repair(c, X0, sh.Phi, med, ro2)
    @printf("  repair feasible=%d min_q=%.6g  (kind=%s F=%d N=%d)\n", Int(rr.feasible), rr.min_q,
            g.kind, K.n_faces(m), K.n_vertices(m))
    rr.feasible || return
    X = rr.X

    B = K.deploy_basis(c, X)
    sd = K.swept_discs(c, B)
    pruned = K.candidate_pairs(c, sd, Float64(pi), true)
    allp = K.candidate_pairs(c, sd, Float64(pi), false)
    r_pruned = K.exact_theta_max_overlap(c, B, pruned, 1e-9, Float64(pi), 1e-9)
    r_all = K.exact_theta_max_overlap(c, B, allp, 1e-9, Float64(pi), 1e-9)
    tb = referee_theta(c, X)
    @printf("  pairs: pruned %d of %d\n", length(pruned), length(allp))
    @printf("  Theta_max  scan(pruned, shrink 1e-9) = %.6f   |C| = %d\n", r_pruned.theta_max,
            length(r_pruned.candidates))
    @printf("  Theta_max  scan(ALL   , shrink 1e-9) = %.6f   |C| = %d\n", r_all.theta_max,
            length(r_all.candidates))
    @printf("  Theta_max  referee bisection          = %.6f\n", tb)

    # where does has_collision fire, and on which pair?
    probes = [1e-7, 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.2484, 0.3]
    for th in probes
        Y = K.deploy(c, X, th).Y
        hc = K.has_collision(c, Y, 1e-12)
        p12 = colliding_pair(c, Y, 1e-12)
        p9 = colliding_pair(c, Y, 1e-9)
        in_pruned = false
        if p12[1] >= 0
            for q in pruned
                (q[1] - 1 == p12[1] && q[2] - 1 == p12[2]) && (in_pruned = true)
            end
        end
        @printf("    th=%-9.6g has_collision(1e-12)=%d  pair@1e-12=(%d,%d) in_pruned=%d  pair@1e-9=(%d,%d)\n",
                th, Int(hc), p12[1], p12[2], Int(in_pruned), p9[1], p9[2])
    end

    # If a pair collides at some theta below theta_exact, report its scan status.
    Yc = K.deploy(c, X, 1e-7).Y
    bad = colliding_pair(c, Yc, 1e-12)
    bad[1] >= 0 || return @printf("\n")
    f = bad[1] + 1; gg = bad[2] + 1     # back to 1-based for indexing
    one = [(f, gg)]
    ca = K.contact_angles(c, B, one, 1e-9, Float64(pi), 1e-9)
    @printf("  offending pair (%d,%d) at theta=1e-7: |C(pair)| = %d, first = %s\n", bad[1],
            bad[2], length(ca), isempty(ca) ? "none" : @sprintf("%f", ca[1]))
    # is it a hinge-adjacent / shared-vertex pair?  (a permanent incidence)
    shared = 0
    for u in c.prime_faces[f], v in c.prime_faces[gg]
        u == v && (shared += 1)
    end
    @printf("  offending pair shares %d M'-vertices (hinge point => boundary contact, not interior overlap)\n", shared)
    for s in (1e-12, 1e-9, 1e-7, 1e-6)
        @printf("    polygons_overlap(shrink=%.0e) = %d\n", s,
                Int(K.polygons_overlap(poly_of(c, Yc, f), poly_of(c, Yc, gg), s)))
    end
    # ---- deep dive on the offending pair -----------------------------------------
    ov(th) = K.polygons_overlap(poly_of(c, K.deploy(c, X, th).Y, f),
                                poly_of(c, K.deploy(c, X, th).Y, gg), 1e-12)
    # bisect the end of the overlap
    lo = 1e-9; hi = 0.1
    if ov(lo) && !ov(hi)
        for i in 1:60
            mid = 0.5 * (lo + hi)
            if ov(mid)
                lo = mid
            else
                hi = mid
            end
        end
        @printf("  pair (%d,%d): interiors overlap on (0, %.9f), clear above\n", bad[1], bad[2], hi)
        # is that angle in C(X) for this pair?  and for ANY pair?
        ca1 = K.contact_angles(c, B, one, 1e-9, 0.5, 1e-9)
        @printf("  C(pair) on (0,0.5): %d angles ->", length(ca1))
        for v in ca1
            @printf(" %.6f", v)
        end
        @printf("\n")
        n_near = count(v -> abs(v - hi) < 1e-4, r_all.candidates)
        @printf("  C(X) angles within 1e-4 of %.9f (over ALL pairs): %d;  min C(X) = %.9f\n",
                hi, n_near, isempty(r_all.candidates) ? -1.0 : r_all.candidates[1])
        # simplicity of the two faces at theta = 0 and at the transition
        function simple_at(ff::Int, th::Float64)
            P = poly_of(c, K.deploy(c, X, th).Y, ff)
            n = length(P)
            for i in 1:n, j in i+1:n
                (j == i + 1 || (i == 1 && j == n)) && continue
                seg_cross(P[i], P[i % n + 1], P[j], P[j % n + 1]) && return false
            end
            return true
        end
        @printf("  face %d simple at 1e-7/mid/0.1: %d %d %d ; face %d: %d %d %d\n",
                bad[1], Int(simple_at(f, 1e-7)), Int(simple_at(f, 0.5 * hi)), Int(simple_at(f, 0.1)),
                bad[2], Int(simple_at(gg, 1e-7)), Int(simple_at(gg, 0.5 * hi)), Int(simple_at(gg, 0.1)))
        # signed areas (POS) of the two faces
        t = deepcopy(m); t.X = X
        @printf("  flat signed areas: face %d = %.6g, face %d = %.6g\n",
                bad[1], K.face_signed_area(t, f), bad[2], K.face_signed_area(t, gg))
        # does the exact scan's own candidate list contain ANY angle below hi?
        below = count(v -> v < hi, r_all.candidates)
        @printf("  C(X) angles strictly below %.9f: %d\n", hi, below)

        # ---- WHO is the witness at the transition, and why did contact_angles miss it?
        thc = hi
        Yt = K.deploy(c, X, thc).Y
        gsc = 0.0
        for y in Yt
            gsc = max(gsc, norm(y))
        end
        function scan_pair(fe::Int, fv::Int)
            E = c.prime_faces[fe]
            V = c.prime_faces[fv]
            ne = length(E)
            for i in 1:ne
                a = E[i]; b = E[i % ne + 1]
                U = Bc(B, b) - Bc(B, a); Vv = Bs(B, b) - Bs(B, a)
                L2 = K.dot_from_vectors(U, Vv, U, Vv)
                for pv in V
                    P = Bc(B, pv) - Bc(B, a); Q = Bs(B, pv) - Bs(B, a)
                    det_ = K.orient_from_vectors(U, Vv, P, Q)
                    D = K.dot_from_vectors(U, Vv, P, Q)
                    l2 = heval(L2, thc); sdot = heval(D, thc)
                    geo = sqrt(max(l2, 0.0)) * norm(Yt[pv] - Yt[a])
                    hval = heval(det_, thc)
                    dist = geo > 0 ? abs(hval) / sqrt(max(l2, 1e-300)) : 0.0
                    on_seg = (sdot >= -1e-9 * l2 && sdot <= l2 * (1 + 1e-9))
                    if dist < 1e-6 * gsc && on_seg
                        skipped_endpoint = (pv == a || pv == b)
                        ident_zero = (abs(det_.p) + hamp(det_) <= 0)
                        idsc = abs(det_.p) + abs(det_.q) + abs(det_.r)
                        @printf("    witness edge(%d,%d) of face %d  vertex %d of face %d : dist=%.3e s/|e|^2=%.6f  skipped(p==a||p==b)=%d ident_zero=%d |p|+|q|+|r|=%.3e  h(0)=%.3e\n",
                                a - 1, b - 1, fe - 1, pv - 1, fv - 1, dist, sdot / l2, Int(skipped_endpoint),
                                Int(ident_zero), idsc, heval(det_, 0.0))
                        if !skipped_endpoint && !ident_zero
                            rts = K.harmonic_roots_deflated(det_, 1e-9, Float64(pi))
                            @printf("      deflated roots in (0,pi):")
                            for rv in rts
                                @printf(" %.6f", rv)
                            end
                            @printf("   (undeflated:")
                            for rv in K.harmonic_roots(det_, 1e-9, Float64(pi))
                                @printf(" %.6f", rv)
                            end
                            @printf(")\n")
                        end
                    end
                end
            end
        end
        scan_pair(f, gg)
        scan_pair(gg, f)

        # ---- WHERE is the overlap?  sample points strictly inside both.
        shared_pv = 0
        for u in c.prime_faces[f], v2 in c.prime_faces[gg]
            u == v2 && (shared_pv = u)
        end
        for th in (1e-4, 0.01, 0.03, 0.036, 0.04)
            Y2 = K.deploy(c, X, th).Y
            A2 = poly_of(c, Y2, f); B2 = poly_of(c, Y2, gg)
            lo2 = Vec2(1e300, 1e300); hi2 = Vec2(-1e300, -1e300)
            for z in A2
                lo2 = min.(lo2, z); hi2 = max.(hi2, z)
            end
            best = -1.0; bz = Vec2(0, 0)
            G = 400
            for i in 0:G, j in 0:G
                z = Vec2(lo2[1] + (hi2[1] - lo2[1]) * i / G, lo2[2] + (hi2[2] - lo2[2]) * j / G)
                (!pip(A2, z) || !pip(B2, z)) && continue
                d = min(edge_dist(A2, z), edge_dist(B2, z))
                if d > best
                    best = d; bz = z
                end
            end
            nA = 0; nB = 0; nAB = 0
            for i in 0:G, j in 0:G
                z = Vec2(lo2[1] + (hi2[1] - lo2[1]) * i / G, lo2[2] + (hi2[2] - lo2[2]) * j / G)
                ia = pip(A2, z); ib = pip(B2, z)
                ia && (nA += 1); ib && (nB += 1); (ia && ib) && (nAB += 1)
            end
            # radial probe around the shared hinge point, like the R4-a local test
            nloc = 0
            if shared_pv >= 1
                pz = Y2[shared_pv]
                r0 = 1e300
                for z in A2
                    norm(z - pz) > 1e-14 && (r0 = min(r0, norm(z - pz)))
                end
                for k in 1:6
                    rad = 0.4 * r0 * k / 7.0
                    for aa in 0:2879
                        ph = 2 * pi * aa / 2880
                        z = pz + rad * Vec2(cos(ph), sin(ph))
                        (pip(A2, z) && pip(B2, z)) && (nloc += 1)
                    end
                end
            end
            @printf("      grid inside A=%d B=%d BOTH=%d ; radial probes at the hinge inside BOTH=%d\n",
                    nA, nB, nAB, nloc)
            # which edge pair does the shrunk-polygon test call a crossing?
            for sk in (0.0, 1e-9)
                ca_ = sum(A2) / length(A2)
                As = [ca_ + (z - ca_) * (1.0 - sk) for z in A2]
                cb_ = sum(B2) / length(B2)
                Bs_ = [cb_ + (z - cb_) * (1.0 - sk) for z in B2]
                na = length(As); nb = length(Bs_)
                for i in 1:na, j in 1:nb
                    a3 = As[i]; b3 = As[i % na + 1]; c3 = Bs_[j]; d3 = Bs_[j % nb + 1]
                    e1 = cross2(b3 - a3, c3 - a3); e2 = cross2(b3 - a3, d3 - a3)
                    e3 = cross2(d3 - c3, a3 - c3); e4 = cross2(d3 - c3, b3 - c3)
                    if ((e1 > 0 && e2 < 0) || (e1 < 0 && e2 > 0)) && ((e3 > 0 && e4 < 0) || (e3 < 0 && e4 > 0))
                        pa = c.prime_faces[f][i]; pb = c.prime_faces[f][i % na + 1]
                        pc = c.prime_faces[gg][j]; pd = c.prime_faces[gg][j % nb + 1]
                        @printf("      shrink=%.0e SPURIOUS crossing: A edge(%d,%d) x B edge(%d,%d) crosses=%.2e %.2e %.2e %.2e  hinge-incident=%d\n",
                                sk, pa - 1, pb - 1, pc - 1, pd - 1, e1, e2, e3, e4,
                                Int(pa == shared_pv || pb == shared_pv || pc == shared_pv || pd == shared_pv))
                    end
                end
            end
            dh = shared_pv >= 1 ? norm(bz - Y2[shared_pv]) : -1.0
            @printf("  th=%-8.5g deepest common interior point: depth=%.4e  dist to hinge=%.4e (hinge pv=%d)  polygons_overlap(1e-12)=%d\n",
                    th, best, best > 0 ? dh : -1.0, shared_pv - 1, Int(K.polygons_overlap(A2, B2, 1e-12)))
            @printf("      polygons_overlap shrink 0 / 1e-9 / 1e-6 / 1e-4 = %d %d %d %d\n",
                    Int(K.polygons_overlap(A2, B2, 0.0)), Int(K.polygons_overlap(A2, B2, 1e-9)),
                    Int(K.polygons_overlap(A2, B2, 1e-6)), Int(K.polygons_overlap(A2, B2, 1e-4)))
        end
        # ---- dump the two polygons at theta = 0.01 as literals, for a deterministic
        #      regression case in derivation_tests (no cache, no repair, no generator).
        Yd = K.deploy(c, X, 0.01).Y
        for ff in (f, gg)
            @printf("  # face %d (%d vertices), theta = 0.01\n  P%d = Vec2[", ff - 1, length(c.prime_faces[ff]), ff == f ? 1 : 2)
            k2 = 0
            for pv in c.prime_faces[ff]
                @printf("%sVec2(%.17g, %.17g)", (k2 > 0 ? ", " : ""), Yd[pv][1], Yd[pv][2])
                k2 += 1
            end
            @printf("]\n")
        end

        # hinge angles and beta_e
        if shared_pv >= 1
            function ang_at(ff::Int, th::Float64)
                Y2 = K.deploy(c, X, th).Y
                PF = c.prime_faces[ff]
                n = length(PF)
                i = 0
                for j in 1:n
                    PF[j] == shared_pv && (i = j)
                end
                i < 1 && return -1.0
                a2 = Y2[PF[mod1(i - 1, n)]] - Y2[shared_pv]; b2 = Y2[PF[mod1(i + 1, n)]] - Y2[shared_pv]
                t2 = atan(cross2(b2, a2), dot(a2, b2))
                t2 < 0 && (t2 += 2 * pi)
                return t2
            end
            af = ang_at(f, 0.01); ag2 = ang_at(gg, 0.01)
            @printf("  hinge angles at th=0.01: alpha_%d=%.6f alpha_%d=%.6f  beta_e=%.6f\n",
                    bad[1], af, bad[2], ag2, 2 * pi - af - ag2)
        end
    else
        @printf("  pair (%d,%d): ov(1e-9)=%d ov(0.1)=%d (not the simple on/off shape)\n",
                bad[1], bad[2], Int(ov(lo)), Int(ov(hi)))
    end
    @printf("\n")
end

function main()
    if !isempty(ARGS)
        for a in ARGS
            replay(parse(Int, a))
        end
    else
        replay(93)
        replay(96)
    end
end

main()
