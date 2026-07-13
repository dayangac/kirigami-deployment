# check_t5_adm.jl -- Checker round 6.  Two questions about the AMENDED certificate
#   VALID(eps) = POS /\ NOOVERLAP(eps/2) /\ NOROOT_adm(eps)
# where NOROOT_adm asks for no ADMISSIBLE root (vertex on the edge SEGMENT), i.e.
# C(X) ^ (0,eps) = {} exactly (results/kill/jitter/cert_diagnosis.md Sec. 4).
# Port of check_t5_adm.cpp (same corpus, same seeds, same printed lines).
#
# Run (from repo root):
#   julia --project=Kirigami derivations/scratch/check_t5_adm.jl
#
# A  SOUNDNESS side, Sub-lemma T5.2b'' Case A.  The substitute pair pi' of core.md
#    (far endpoint w of the SHORTER far-side edge, against the other far-side edge
#    (a=p, b)) has h_o,pi'(beta_e) = 0.  The amendment needs more: that root must be
#    ADMISSIBLE.  Claim (proved in check.md Round 6, R6.1): at theta = beta_e the two
#    far-side edges are collinear and CO-DIRECTED, so
#        <w - a, b - a> = L' L   and   |b - a|^2 = L^2 ,   with  0 < L' <= L,
#    hence the projection ratio s/|e|^2 = L'/L lies in (0, 1].  Measured here.
#
# B  COMPLETENESS side.  Is "exact Theta_max >= eps  =>  certificate" true?  NO.
#    A GRAZE -- an admissible contact that is not an overlap transition -- sits in
#    (0, Theta_max) on the hexagon pattern (core.md T4.2 table: theta_1 = pi/3,
#    Theta_max = 2pi/3).  For any eps strictly between the two, Theta_max >= eps holds
#    while NOROOT_adm fails.  The 0/5064 exceptions of the jitter run are an artefact of
#    eps = 0.006 being far below every contact angle on that corpus.
#
# Depends on Kirigami/src/method/{deploy_basis,contact}.jl: deploy_basis, orient_from_vectors,
# dot_from_vectors, swept_discs, candidate_pairs, contact_angles, exact_theta_max_overlap,
# validity_certificate (C++ names, per PORTING.md).

using Kirigami
using LinearAlgebra
using Printf
using StaticArrays
const K = Kirigami
const Vec2 = K.Vec2

# row i of a DeployBasis matrix as a Vec2 (the C++ DeployBasis::c(i) / s(i))
Bc(B, i::Int) = Vec2(B.C[i, 1], B.C[i, 2])
Bs(B, i::Int) = Vec2(B.S[i, 1], B.S[i, 2])
heval(h, th::Float64) = h.p + h.q * cos(th) + h.r * sin(th)
hamp(h) = hypot(h.q, h.r)

mutable struct Case
    name::String
    m::K.Mesh
    c::K.CutStructure
    hs::K.HoleSet
    sys::K.LinearSystem
    rep::K.SolveReport
    X0::Vector{Vec2}
    k::Int
    ok::Bool
end

function sample(cs::Case, rng::K.MT19937, scale::Float64)
    X = copy(cs.rep.X0)
    if cs.k > 0 && scale != 0.0
        g = K.NormalDist(0.0, scale)
        T = Matrix{Float64}(undef, cs.k, 2)
        for i in 1:cs.k
            T[i, 1] = K.normal(g, rng); T[i, 2] = K.normal(g, rng)
        end
        X += cs.rep.Phi * T
    end
    return K.matrix_to_points(X)
end

function corpus()
    out = Case[]
    rng = K.MT19937(20260903)
    function add(m::K.Mesh, nm::String)
        K.n_faces(m) < 3 && return
        K.build_topology!(m)
        m.sigma = K.assign_orientation_relaxation(m, rng, 6, 400, 90).sigma
        K.build_topology!(m)
        c = K.make_cut(m)
        hs = K.holes_partition(c)
        sys = K.assemble_system(c, hs, m.X, K.Fixed)
        rep = K.solve_system(sys, m.X)
        X0 = K.matrix_to_points(rep.X0)
        r = K.hole_residuals(c, X0, hs)
        ok = rep.projection_ok && r.max_norm < 1e-7 && K.n_hinge(c) > 0
        ok && push!(out, Case(nm, m, c, hs, sys, rep, X0, rep.dim_null, ok))
    end
    add(K.tiling_squares(K.disk(Vec2(0, 0), 2.6)), "squares")
    add(K.tiling_triangles(K.disk(Vec2(0, 0), 2.2)), "triangles")
    add(K.tiling_hexagons(K.disk(Vec2(0, 0), 3.0)), "hexagons")
    add(K.tiling_kagome(K.disk(Vec2(0, 0), 2.4)), "kagome")
    add(K.tiling_snub_square(K.disk(Vec2(0, 0), 2.2)), "snub_square")
    add(K.tiling_truncated_square(K.disk(Vec2(0, 0), 2.8)), "truncated_square")
    add(K.tiling_3_4_3_12(K.disk(Vec2(0, 0), 3.2)), "t3_4_3_12")
    return out
end

function interior_angle(vs::Vector{Int}, X::Vector{Vec2}, v::Int)
    n = length(vs)
    i = 0
    for j in 1:n
        vs[j] == v && (i = j)
    end
    i < 1 && return -1.0
    a = X[vs[mod1(i - 1, n)]] - X[v]; b = X[vs[mod1(i + 1, n)]] - X[v]
    (norm(a) < 1e-14 || norm(b) < 1e-14) && return -1.0
    t = atan(b[1] * a[2] - b[2] * a[1], dot(a, b))
    t < 0 && (t += 2 * pi)
    return t   # interior angle at v (same convention as derivation_tests)
end

function main()
    C = corpus()
    @printf("corpus: %d cases\n\n", length(C))

    # ==================================================================== A
    # Case A: admissibility of the substitute pair's root at theta = beta_e.
    @printf("A  Case A substitute pair pi' = (w, (a,b)) at theta = beta_e\n")
    @printf("   admissible <=> 0 <= <w-a,b-a> <= |b-a|^2 at beta_e\n\n")
    nA = 0; nA_bad = 0; nA_degen = 0; n_beta_seen = 0; n_beta_out = 0; nA_strict = 0
    worst_ratio_lo = 1e300; worst_ratio_hi = -1e300; e_ratio = 0.0; e_root = 0.0
    rngA = K.MT19937(60061)
    for cs in C
        for rep in 0:39
            Xs = sample(cs, rngA, rep == 0 ? 0.0 : 0.12)
            B = K.deploy_basis(cs.c, Xs)
            for e in cs.c.hinge_edges
                ed = cs.m.edges[e]
                f = cs.m.half_edges[ed.he[1]].face
                g = cs.m.half_edges[ed.he[2]].face
                (f < 1 || g < 1) && continue
                v = cs.c.hinge_dir[e].src; dst = cs.c.hinge_dir[e].dst
                vf = cs.m.faces[f]
                vg = cs.m.faces[g]
                af = interior_angle(vf, Xs, v); ag = interior_angle(vg, Xs, v)
                (af < 0 || ag < 0) && continue
                beta = 2 * pi - af - ag
                n_beta_seen += 1
                if !(beta > 1e-3 && beta < pi - 1e-3)   # the case the proof uses
                    n_beta_out += 1
                    continue
                end

                # far-side neighbour of v in each face: the one that is NOT dst
                function far_of(vs::Vector{Int})
                    n = length(vs)
                    i = 0
                    for j in 1:n
                        vs[j] == v && (i = j)
                    end
                    pv = vs[mod1(i - 1, n)]; nx = vs[mod1(i + 1, n)]
                    far = (nx == dst) ? pv : nx
                    fi = 0
                    for j in 1:n
                        vs[j] == far && (fi = j)
                    end
                    return (i, fi)
                end
                (if_v, if_far) = far_of(vf)
                (ig_v, ig_far) = far_of(vg)
                (if_far < 1 || ig_far < 1 || if_v < 1 || ig_v < 1) && continue
                lf = norm(Xs[vf[if_far]] - Xs[v])
                lg = norm(Xs[vg[ig_far]] - Xs[v])
                if lf <= lg   # shorter far edge belongs to f: its far endpoint is w
                    w = cs.c.prime_faces[f][if_far]
                    a = cs.c.prime_faces[g][ig_v]
                    b = cs.c.prime_faces[g][ig_far]
                    Lp = lf; L = lg
                else
                    w = cs.c.prime_faces[g][ig_far]
                    a = cs.c.prime_faces[f][if_v]
                    b = cs.c.prime_faces[f][if_far]
                    Lp = lg; L = lf
                end
                U = Bc(B, b) - Bc(B, a); Vv = Bs(B, b) - Bs(B, a)
                P = Bc(B, w) - Bc(B, a); Q = Bs(B, w) - Bs(B, a)
                det = K.orient_from_vectors(U, Vv, P, Q)
                D = K.dot_from_vectors(U, Vv, P, Q)
                L2 = K.dot_from_vectors(U, Vv, U, Vv)
                sc = abs(det.p) + abs(det.q) + abs(det.r)
                gsc = norm(U) * norm(P)
                if sc <= 1e-11 * max(gsc, 1e-300)
                    nA_degen += 1
                    continue
                end
                nA += 1
                e_root = max(e_root, abs(heval(det, beta)) / max(sc, 1e-300))
                s = heval(D, beta); l2 = heval(L2, beta)
                ratio = s / l2
                worst_ratio_lo = min(worst_ratio_lo, ratio)
                worst_ratio_hi = max(worst_ratio_hi, ratio)
                ratio < 1.0 - 1e-9 && (nA_strict += 1)
                e_ratio = max(e_ratio, abs(ratio - Lp / L))
                tol = 1e-12 * (abs(L2.p) + hamp(L2))
                if s < -tol || s > l2 + tol
                    nA_bad += 1
                    @printf("   [ADM-FAIL] %-18s beta=%.6f  s/|e|^2=%.6f  L'/L=%.6f\n", cs.name,
                            beta, ratio, Lp / L)
                end
            end
        end
    end
    @printf("   hinge edges with both angles    : %d  (beta_e outside (0,pi): %d)\n", n_beta_seen, n_beta_out)
    @printf("   pairs tested                    : %d  (identically-zero, skipped: %d)\n", nA, nA_degen)
    @printf("   max |h(beta_e)|/scale           : %.3e   (root of the substitute pair)\n", e_root)
    @printf("   projection ratio s/|e|^2 range  : [%.6f, %.6f]  (claim: (0, 1])\n",
            worst_ratio_lo, worst_ratio_hi)
    @printf("   pairs with ratio < 1 - 1e-9    : %d  (L' < L: vertex STRICTLY inside)\n", nA_strict)
    @printf("   max |s/|e|^2 - L'/L|            : %.3e   (claim: equal)\n", e_ratio)
    @printf("   ADMISSIBILITY FAILURES          : %d\n\n", nA_bad)

    # ==================================================================== B
    # Completeness: a graze inside (0, Theta_max) breaks "Theta_max >= eps => cert".
    @printf("B  completeness: Theta_max >= eps  =>  certificate ?\n\n")
    @printf("   %-18s %8s %10s %10s  %10s  %s\n", "pattern", "|C|", "theta_1", "Theta_max", "eps",
            "cert(pos,noov,noroot)")
    n_counterex = 0
    for cs in C
        B = K.deploy_basis(cs.c, cs.X0)
        sd = K.swept_discs(cs.c, B)
        pairs = K.candidate_pairs(cs.c, sd, Float64(pi), false)
        cand = K.contact_angles(cs.c, B, pairs, 1e-9, Float64(pi), 1e-9)
        rep = K.exact_theta_max_overlap(cs.c, B, pairs, 1e-9, Float64(pi), 1e-9)
        if isempty(cand)
            @printf("   %-18s %8d %10s %10.6f  %10s  %s\n", cs.name, 0, "-", rep.theta_max,
                    "-", "(no contact in range)")
            continue
        end
        th1 = cand[1]; Tm = rep.theta_max
        if !(th1 < Tm - 1e-6)
            @printf("   %-18s %8d %10.6f %10.6f  %10s  %s\n", cs.name, length(cand), th1, Tm,
                    "-", "(no graze: theta_1 == Theta_max)")
            continue
        end
        eps = 0.5 * (th1 + Tm)   # strictly between the graze and the true range
        cert = K.validity_certificate(cs.c, B, cs.X0, pairs, eps, 1e-9)
        valid = cert.pos && cert.nooverlap && cert.noroot
        @printf("   %-18s %8d %10.6f %10.6f  %10.6f  (%d,%d,%d)%s\n", cs.name,
                length(cand), th1, Tm, eps, Int(cert.pos), Int(cert.nooverlap), Int(cert.noroot),
                (Tm >= eps && !valid) ? "   <-- COUNTEREXAMPLE to completeness" : "")
        if Tm >= eps && !valid
            n_counterex += 1
            @printf("        admissible root at %.6f in (0,eps) though no interior overlap until %.6f\n",
                    cert.first_root, Tm)
        end
    end
    @printf("\n   completeness counterexamples: %d\n", n_counterex)
    @printf("\nVERDICT  A: %s   B: %s\n", nA_bad == 0 ? "Case A root is ADMISSIBLE (sound)" : "FAIL",
            n_counterex > 0 ? "completeness is FALSE" : "no counterexample found")
end

main()
