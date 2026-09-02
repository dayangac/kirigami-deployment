# test_method_1.jl -- port of the tests/test_method.cpp cases that exercise
# method/{deploy_basis, mobility, contact, range_opt}, case by case.
#
# Inputs: the C++ `make_case` tilings (checkerboard or relaxation sigma, Eq. (6)
# projection), the delaunay meshes of the mobility cases and the range-objective case were
# frozen by replaying test_method.cpp exactly into
# CORPUS/method_fixtures/test_method_1.json (producer: freeze_method_1.cpp next to it,
# linked against the C++ libkiri_core.a). Each block carries a `provenance` string with
# the generator call, seed and sigma method. Besides the inputs the fixture holds every
# number the C++ computed (theta_max variants, contact-angle sets, certificate counters,
# mobility reports, objective value and gradient), which the tests below check against in
# addition to the C++ CHECKs themselves. The synthetic-harmonic cases draw their
# coefficients from the bit-exact MT19937 with the C++ seeds.
#
# Bit-exactness: on the arm64 Julia (PORTING.md) the deploy basis, every harmonic
# coefficient, every root, C(X), the witnesses and all certificate counters reproduce the C++
# bit for bit (fma placed as clang contracts it, __sincos_stret for sin/cos pairs, unfused
# Eigen dots, libm atan2/acos/hypot/atan), and the tie-sensitive fields are asserted exactly.
# On a non-arm64 Julia (x86_64 under Rosetta: different libm ulps) the tie-sensitive
# DIAGNOSTIC counters -- ExactRangeReport.n_roots, ValidityCertificate.n_roots_undeflated,
# the 1e-9 dedup length of C(X) -- are asserted loosely (M1_EXACT_TIES = false).
include("helpers.jl")
import JSON
using SparseArrays: spzeros

const K = Kirigami
const M1_EXACT_TIES = K._USE_SYSTEM_LIBM
const M1 = JSON.parsefile(joinpath(CORPUS, "method_fixtures", "test_method_1.json"))

struct M1Case
    m::K.Mesh
    c::K.CutStructure
    hs::K.HoleSet
    X::Vector{K.Vec2}
end

m1_points(P) = [K.Vec2(Float64(p[1]), Float64(p[2])) for p in P]
m1_matrix(A) = isempty(A) ? zeros(0, 0) : Float64[Float64(A[i][j]) for i in eachindex(A), j in eachindex(A[1])]
# fixture face pairs are 0-based (C++) -> 1-based
m1_pairs(P) = Tuple{Int,Int}[(Int(p[1]) + 1, Int(p[2]) + 1) for p in P]

# The frozen mesh carries its sigma under "orientation"; X is the C++ make_case's X.
function m1_case(rec)
    m = K.mesh_from_json_string(JSON.json(rec["mesh"]))
    m.sigma = Int[s for s in rec["sigma"]]
    c = K.make_cut(m)
    hs = K.holes_partition(c)
    return M1Case(m, c, hs, m1_points(rec["X"]))
end
m1_tiling(name) = m1_case(M1["tilings"][name])

# The C++ make_case, on the Julia generators (used only to check the frozen inputs).
function m1_make_case(m::K.Mesh, checker::Bool)
    rng = K.MT19937(2026)
    K.build_topology!(m)
    m.sigma = checker ? K.checkerboard(m) : K.assign_orientation_relaxation(m, rng, 8, 500, 180).sigma
    K.build_topology!(m)
    c = K.make_cut(m)
    hs = K.holes_partition(c)
    X = m.X
    if !K.deployable(K.hole_residuals(c, X, hs), 1e-9)
        sys = K.assemble_system(c, hs, m.X, K.Fixed)
        sr = K.solve_system(sys, m.X)
        X = K.matrix_to_points(sr.X0)
    end
    return M1Case(m, c, hs, X)
end

# Bisection on the true polygon overlap. collision.hpp's default shrink of 1e-6
# (README deviation 9) reports first contact up to ~5e-5 PAST the true angle on the
# 4.8.8 pattern -- measured, see results/kill/KILL_REPORT.md K2a -- so the reference
# used here shrinks by 1e-12 instead.
function bisect_theta_max(c::K.CutStructure, X::Vector{K.Vec2}, shrink::Float64)
    col(th) = K.has_collision(c, K.deploy(c, X, th).Y, shrink)
    lo = 0.0
    hi = -1.0
    for i in 1:360
        th = pi * i / 360
        if col(th)
            hi = th
            break
        end
        lo = th
    end
    hi < 0 && return Float64(pi)
    for _ in 1:60
        mid = 0.5 * (lo + hi)
        if col(mid)
            hi = mid
        else
            lo = mid
        end
    end
    return lo
end

const M1_TILING_GEN = Dict(
    "squares_checker" => (() -> K.tiling_squares(K.rect(K.Vec2(2.5, 2.5), 2.51, 2.51)), true),
    "triangles_checker" => (() -> K.tiling_triangles(K.disk(K.Vec2(0.13, 0.07), 2.5)), true),
    "kagome_checker" => (() -> K.tiling_kagome(K.disk(K.Vec2(0.13, 0.07), 2.5)), true),
    "snub_2_6" => (() -> K.tiling_snub_square(K.disk(K.Vec2(0.13, 0.07), 2.6)), false),
    "snub_3_0" => (() -> K.tiling_snub_square(K.disk(K.Vec2(0.13, 0.07), 3.0)), false),
    "snub_3_2" => (() -> K.tiling_snub_square(K.disk(K.Vec2(0.13, 0.07), 3.2)), false),
    "trunc_3_0" => (() -> K.tiling_truncated_square(K.disk(K.Vec2(0.13, 0.07), 3.0)), false),
    "trunc_4_0" => (() -> K.tiling_truncated_square(K.disk(K.Vec2(0.13, 0.07), 4.0)), false),
    "hexagons_3_0" => (() -> K.tiling_hexagons(K.disk(K.Vec2(0.13, 0.07), 3.0)), false),
    "t3_4_3_12_4_2" => (() -> K.tiling_3_4_3_12(K.disk(K.Vec2(0.13, 0.07), 4.2)), false),
)

@testset "frozen make_case inputs reproduce from the Julia generators" begin
    # Not a C++ case: pins that the fixture inputs are what generators.jl / orientation.jl /
    # tutte_auxetic.jl produce for the same calls, so the fixture is a cross-check and not
    # a fork of the inputs.
    for (name, (gen, checker)) in M1_TILING_GEN
        fx = M1["tilings"][name]
        cs = m1_make_case(gen(), checker)
        @test K.n_vertices(cs.m) == fx["N"]
        @test K.n_faces(cs.m) == fx["F"]
        @test cs.m.sigma == Int[s for s in fx["sigma"]]
        @test cs.c.n_prime_vertices == fx["n_prime"]
        Xf = m1_points(fx["X"])
        @test length(cs.X) == length(Xf)
        length(cs.X) == length(Xf) && @test maximum(norm.(cs.X .- Xf)) < 1e-10
    end
end

@testset "deploy_basis reproduces deploy() exactly (U3/U6)" begin
    for name in ("squares_checker", "kagome_checker", "snub_2_6", "trunc_3_0")
        cs = m1_tiling(name)
        B = K.deploy_basis(cs.c, cs.X)
        scale = 0.0
        for i in 1:K.n_prime(B)
            scale = max(scale, norm(K.basis_c(B, i)))
        end
        @test scale > 0
        for th in (0.05, 0.4, 1.1, 2.0, 3.0)
            Y = K.deploy(cs.c, cs.X, th).Y
            Yb = K.basis_eval(B, th)
            err = maximum(norm(Y[i] - Yb[i]) for i in eachindex(Y))
            @test err < 1e-12 * scale
        end
        # and equals the C++ basis
        fx = M1["tilings"][name]
        @test maximum(abs, B.C - m1_matrix(fx["C"])) < 1e-12
        @test maximum(abs, B.S - m1_matrix(fx["S"])) < 1e-12
    end
end

@testset "face signed areas are constant along the deployment" begin
    cs = m1_tiling("snub_2_6")
    function area(Y, f)
        pf = cs.c.prime_faces[f]
        s = 0.0
        for i in eachindex(pf)
            u = Y[pf[i]]
            v = Y[pf[mod1(i + 1, length(pf))]]
            s += u[1] * v[2] - u[2] * v[1]
        end
        return 0.5 * s
    end
    Y0 = K.deploy(cs.c, cs.X, 0.0).Y
    for th in (0.3, 1.0, 2.2)
        Y = K.deploy(cs.c, cs.X, th).Y
        for f in 1:K.n_faces(cs.m)
            a0 = area(Y0, f)
            @test abs(area(Y, f) - a0) < 1e-11 * max(1.0, abs(a0))
        end
    end
end

@testset "orientation determinant is the closed-form harmonic" begin
    cs = m1_tiling("trunc_3_0")
    B = K.deploy_basis(cs.c, cs.X)
    rng = K.MT19937(99)
    checked = 0
    for t in 1:300
        checked >= 200 && break
        pf = cs.c.prime_faces[K.uniform_int(rng, 0, K.n_faces(cs.m) - 1) + 1]
        i = K.uniform_int(rng, 0, length(pf) - 1) + 1
        a = pf[i]; b = pf[mod1(i + 1, length(pf))]
        p = K.uniform_int(rng, 0, cs.c.n_prime_vertices - 1) + 1
        (p == a || p == b) && continue
        h = K.orient_harmonic(B, a, b, p)
        d = K.dot_harmonic(B, a, b, p)
        l = K.len2_harmonic(B, a, b)
        geo = 0.0
        for th in (0.13, 0.9, 1.7, 2.6)
            Y = K.deploy(cs.c, cs.X, th).Y
            u = Y[b] - Y[a]; v = Y[p] - Y[a]
            geo = max(geo, norm(u) * norm(v))
            @test abs(K.harmonic_eval(h, th) - (u[1] * v[2] - u[2] * v[1])) < 1e-10 * max(1.0, geo)
            @test abs(K.harmonic_eval(d, th) - dot(u, v)) < 1e-10 * max(1.0, geo)
            @test abs(K.harmonic_eval(l, th) - dot(u, u)) < 1e-10 * max(1.0, geo)
        end
        checked += 1
    end
    @test checked > 100
end

@testset "harmonic_roots finds every root and no spurious ones" begin
    rng = K.MT19937(7)
    for t in 1:500
        h = K.Harmonic(K.uniform_real(rng, -2.0, 2.0), K.uniform_real(rng, -2.0, 2.0),
                       K.uniform_real(rng, -2.0, 2.0))
        roots = K.harmonic_roots(h, 0.0, 2pi)
        for r in roots
            @test abs(K.harmonic_eval(h, r)) < 1e-9 * max(1.0, K.scale(h))
        end
        # grid scan for sign changes must find the same count
        changes = 0
        G = 20000
        prev = K.harmonic_eval(h, 0.0)
        for i in 1:G
            v = K.harmonic_eval(h, 2pi * i / G)
            ((prev > 0 && v < 0) || (prev < 0 && v > 0)) && (changes += 1)
            prev = v
        end
        @test length(roots) >= changes
        @test length(roots) <= changes + 1
    end
end

@testset "harmonic_fit recovers exact coefficients" begin
    th = Float64[]
    y = Float64[]
    truth = K.Harmonic(0.7, -1.3, 2.1)
    for i in 1:50
        push!(th, pi * i / 50.0)
        push!(y, K.harmonic_eval(truth, th[end]))
    end
    f, mr, sc = K.harmonic_fit(th, y)
    @test mr < 1e-12
    @test isapprox(f.p, truth.p; rtol = 1e-12)
    @test isapprox(f.q, truth.q; rtol = 1e-12)
    @test isapprox(f.r, truth.r; rtol = 1e-12)
end

@testset "closed-form theta_max agrees with the collision bisection" begin
    for name in ("squares_checker", "triangles_checker", "kagome_checker", "trunc_3_0")
        cs = m1_tiling(name)
        fx = M1["tilings"][name]
        flat = K.has_collision(cs.c, K.deploy(cs.c, cs.X, 0.0).Y)
        @test flat == fx["flat_collides_default"]
        flat && continue  # flat state not embedded
        B = K.deploy_basis(cs.c, cs.X)
        sd = K.swept_discs(cs.c, B)
        all = K.candidate_pairs(cs.c, sd, pi, false)
        pruned = K.candidate_pairs(cs.c, sd, pi, true)
        te = K.exact_theta_max(cs.c, B, all).theta_max
        tp = K.exact_theta_max(cs.c, B, pruned).theta_max
        @test abs(te - tp) < 1e-12
        @test length(pruned) <= length(all)
        tb = bisect_theta_max(cs.c, cs.X, 1e-12)
        @test abs(te - tb) < 1e-5
        # and the shipped bisection is within its own documented slack
        @test K.theta_max(cs.c, cs.X, 90, 40).theta_max_geometric >= te - 1e-9
        # C++ values
        @test length(all) == fx["n_pairs_all"]
        @test length(pruned) == fx["n_pairs_pruned"]
        @test pruned == m1_pairs(fx["pairs_pruned"])
        @test isapprox(te, fx["theta_exact_all"]; rtol = 1e-12, atol = 1e-13)
        @test isapprox(tb, fx["theta_bisect_1e12"]; rtol = 1e-12, atol = 1e-13)
    end
end

@testset "contact: every C++ number of the exact scan, per frozen tiling" begin
    # Not a C++ case: the C++ ExactRangeReport / OverlapRangeReport / swept radii on each
    # of the 10 make_case tilings, field by field.
    for (name, fx) in M1["tilings"]
        cs = m1_tiling(name)
        B = K.deploy_basis(cs.c, cs.X)
        sd = K.swept_discs(cs.c, B)
        @test maximum(abs.(sd.rho_max .- Float64.(fx["rho_max"]))) < 1e-12
        @test maximum(abs.(sd.circum .- Float64.(fx["circum"]))) < 1e-12
        @test maximum(abs.(sd.rho .- Float64.(fx["rho"]))) < 1e-12
        pruned = K.candidate_pairs(cs.c, sd, pi, true)
        @test pruned == m1_pairs(fx["pairs_pruned"])
        ep = K.exact_theta_max(cs.c, B, pruned)
        e = fx["exact_pruned"]
        @test ep.n_pairs == e["n_pairs"]
        @test ep.n_candidates == e["n_candidates"]
        # n_roots counts roots at or below the running best and so depends on the order in
        # which exact ties (regular tilings have many) resolve at the 1e-15 slack.
        if M1_EXACT_TIES
            @test ep.n_roots == e["n_roots"]
        else
            ep.n_roots == e["n_roots"] || @info "$name: n_roots $(ep.n_roots) (C++ $(e["n_roots"])), tie-order dependent"
            @test ep.n_roots >= (e["found"] ? 1 : 0)
        end
        @test ep.n_zero_contacts == e["n_zero_contacts"]
        @test ep.first.found == e["found"]
        @test isapprox(ep.theta_max, fx["theta_exact_pruned"]; rtol = 1e-12, atol = 1e-13)
        if e["found"]
            # C++ M'-vertex / face ids are 0-based
            @test ep.first.pv == e["pv"] + 1
            @test ep.first.pa == e["pa"] + 1
            @test ep.first.pb == e["pb"] + 1
            @test ep.first.face_v == e["face_v"] + 1
            @test ep.first.face_e == e["face_e"] + 1
            @test ep.first.corner_e == e["corner_e"] + 1
        end
        ov = K.exact_theta_max_overlap(cs.c, B, pruned, 1e-9, pi, 1e-9)
        o = fx["overlap_pruned"]
        @test isapprox(ov.theta_max, fx["theta_overlap_pruned"]; rtol = 1e-12, atol = 1e-13)
        cc = Float64.(o["candidates"])
        if M1_EXACT_TIES
            @test ov.candidates == cc
        else
            # C(X) as a set at 1e-6 resolution: a tangential (double) root has |p| = amp, so
            # its existence and position (acos near 1: error ~ sqrt(ulp) ~ 1e-8) are decided
            # by the last ulp, and the 1e-9 dedup chain can then split or merge one
            # near-coincident pair (trunc_4_0: an extra root 2.6e-8 from 3pi/4).
            @test all(any(abs(v - w) < 1e-6 for w in cc) for v in ov.candidates)
            @test all(any(abs(v - w) < 1e-6 for v in ov.candidates) for w in cc)
            @test abs(length(ov.candidates) - length(cc)) <= 1
        end
        @test ov.i_star == o["i_star"]
        @test ov.n_intervals_tested == o["n_intervals_tested"]
        @test ov.zero_range == o["zero_range"]
        @test ov.n_overlap_tests == o["n_overlap_tests"]
        @test ov.first_contact.found == o["first_contact_found"]
        @test isapprox(ov.first_contact.theta, o["first_contact_theta"]; rtol = 1e-12, atol = 1e-13)
        @test isapprox(K.theta_max(cs.c, cs.X, 90, 40).theta_max_geometric, fx["theta_shipped"];
                       rtol = 1e-12, atol = 1e-13)
    end
end

@testset "A: sigma is in its kernel exactly when Eq. (2) holds" begin
    # rng mt19937(31337): delaunay_of_random_points(40 + 5t, 3.0, rng), relaxation(4, 300, 90)
    for rec in M1["mobility_sigma_kernel"]["cases"]
        rec["skipped"] && continue
        m = K.mesh_from_json_string(JSON.json(rec["mesh"]))
        m.sigma = Int[s for s in rec["sigma"]]
        c = K.make_cut(m)
        hs = K.holes_partition(c)
        g = K.build_hinge_graph(c)
        @test K.n_edges(g) == rec["n_hinge_edges"]
        K.n_edges(g) == 0 && continue
        hg = rec["hinge_graph"]
        @test g.F == hg["F"]
        @test g.components == hg["components"]
        @test K.n_cycles(g) == hg["n_cycles"]
        @test g.head == Int[h + 1 for h in hg["head"]]   # 0-based face ids in the fixture
        @test g.tail == Int[t + 1 for t in hg["tail"]]
        @test g.eid == Int[e + 1 for e in hg["eid"]]
        @test g.nontree == Int[i + 1 for i in hg["nontree"]]
        # X_ini is generically NOT deployable; the Eq. (6) projection is.
        r_ini = K.mobility_at(c, g, K.pins_flat(c, g, m.X))
        X0 = m1_points(rec["X0"])
        r_0 = K.mobility_at(c, g, K.pins_flat(c, g, X0))
        @test r_ini.sigma_in_ker == K.deployable(K.hole_residuals(c, m.X, hs), 1e-9)
        @test r_0.sigma_in_ker == K.deployable(K.hole_residuals(c, X0, hs), 1e-9)
        @test r_0.sigma_in_ker
        # the 2-core identity
        @test r_0.identity_holds
        @test r_ini.identity_holds
        # the C++ reports, field by field
        for (r, fr) in ((r_ini, rec["r_ini"]), (r_0, rec["r_0"]))
            for fld in (:F, :n_hinge, :n_cycles, :components, :dim_ker_A, :n_dangling,
                        :core_faces, :core_edges, :core_cycles, :c_core, :dim_ker_A_core,
                        :identity_holds, :m_full, :m_core, :sigma_in_ker, :used_sparse)
                @test getfield(r, fld) == fr[string(fld)]
            end
            @test isapprox(r.sigma_residual, fr["sigma_residual"]; atol = 1e-9)
        end
        # the Julia Eq. (6) projection lands on the same X0
        sys = K.assemble_system(c, hs, m.X, K.Fixed)
        sr = K.solve_system(sys, m.X)
        @test maximum(norm.(K.matrix_to_points(sr.X0) .- X0)) < 1e-9
    end
end

@testset "A and the body-and-pin rigidity matrix give the same mobility" begin
    # rng mt19937(4711): delaunay_of_random_points(25 + 4t, 3.0, rng), relaxation(4, 300, 90)
    for rec in M1["mobility_rigidity"]["cases"]
        rec["skipped"] && continue
        m = K.mesh_from_json_string(JSON.json(rec["mesh"]))
        m.sigma = Int[s for s in rec["sigma"]]
        c = K.make_cut(m)
        g = K.build_hinge_graph(c)
        K.n_edges(g) == 0 && continue
        pins = K.pins_flat(c, g, m.X)
        r = K.mobility_at(c, g, pins)
        R = K.build_rigidity(g, pins)
        # Eigen ColPivHouseholderQR rank with threshold 1e-10 on R / max|R|
        F = qr(R / maximum(abs, R), ColumnNorm())
        d = abs.(diag(F.R))
        rank_R = count(>(1e-10 * maximum(d)), d)
        m_R = 3 * g.F - rank_R - 3 * g.components
        @test r.m_full == m_R
        @test rank_R == rec["rank_R"]
        @test m_R == rec["m_R"]
        @test r.m_full == rec["r"]["m_full"]
        @test r.dim_ker_A == rec["r"]["dim_ker_A"]
    end
end

@testset "matrix_rank: the SparseQR path agrees with the dense one" begin
    # Not a C++ case. The C++ tests never exceed `dense_limit` = 700 columns, so the sparse
    # branch (SPQR here, Eigen::SparseQR there) is exercised by forcing it on a small A.
    for rec in M1["mobility_sigma_kernel"]["cases"]
        rec["skipped"] && continue
        m = K.mesh_from_json_string(JSON.json(rec["mesh"]))
        m.sigma = Int[s for s in rec["sigma"]]
        c = K.make_cut(m)
        g = K.build_hinge_graph(c)
        A = K.build_A(g, K.pins_flat(c, g, m.X))
        rd, sd = K.matrix_rank(A, 1e-10)
        rs, ss = K.matrix_rank(A, 1e-10, 1)
        @test !sd
        @test ss
        @test rd == rs
        @test g.F - rd == rec["r_ini"]["dim_ker_A"]
    end
    @test K.matrix_rank(spzeros(0, 5)) == (0, false)
    @test K.matrix_rank(spzeros(4, 5)) == (0, false)
end

@testset "range objective: analytic gradient matches central differences" begin
    rec = M1["range_objective"]
    cs = m1_tiling("trunc_3_0")
    sys = K.assemble_system(cs.c, cs.hs, cs.m.X, K.Fixed)
    sr = K.solve_system(sys, cs.m.X)
    @test sr.dim_null > 0
    @test sr.dim_null == rec["dim_null"]
    # The null-space basis is only defined up to an orthogonal change of basis, so the
    # C++ Phi is used from here on (the C++ t, f0 and g refer to it).
    X0 = m1_points(rec["X0"])
    @test maximum(norm.(K.matrix_to_points(sr.X0) .- X0)) < 1e-9
    Phi = m1_matrix(rec["Phi"])
    o = K.RangeOptOptions()
    ob = K.RangeObjective()
    K.setup!(ob, cs.c, X0, Phi, o)
    k = size(Phi, 2)
    t = zeros(2k)
    K.rebuild_active!(ob, t)
    @test !isempty(ob.active)
    @test length(ob.active) == length(rec["active_at_0"])
    # same active triples (0-based in the fixture), as a set: std::sort ties are unordered
    @test Set(ob.active) == Set(NTuple{3,Int}((a[1] + 1, a[2] + 1, a[3] + 1)) for a in rec["active_at_0"])
    # the basis reconstructed from (X0, Phi, T) must equal deploy()'s at the same X
    rng = K.MT19937(5)
    G = K.NormalDist(0.0, 1e-3)
    for i in 1:2k
        t[i] = K.normal(G, rng)
    end
    @test maximum(abs.(t .- Float64.(rec["t"]))) < 1e-15   # bit-exact normal stream
    t = Float64.(rec["t"])
    let T = K._unflatten_T(t, k)
        X = K.matrix_to_points(K.points_to_matrix(X0) + Phi * T)
        Bd = K.deploy_basis(cs.c, X)
        Bo = K.basis(ob, t)
        @test maximum(abs, Bd.C - Bo.C) < 1e-10
        @test maximum(abs, Bd.S - Bo.S) < 1e-9
        @test maximum(abs, Bo.C - m1_matrix(rec["basis_C"])) < 1e-12
        @test maximum(abs, Bo.S - m1_matrix(rec["basis_S"])) < 1e-12
    end
    g = zeros(2k)
    f0 = K.value_and_grad(ob, t, g)
    @test isfinite(f0)
    @test isapprox(f0, rec["f0"]; rtol = 1e-9)
    @test maximum(abs.(g .- Float64.(rec["g"]))) < 1e-9 * max(1.0, maximum(abs, g))
    checked = 0
    worst = 0.0
    for trial in 1:12
        i = K.uniform_int(rng, 0, 2k - 1) + 1
        h = 1e-6
        tp = copy(t); tm = copy(t); dummy = zeros(2k)
        tp[i] += h
        tm[i] -= h
        fd = (K.value_and_grad(ob, tp, dummy) - K.value_and_grad(ob, tm, dummy)) / (2h)
        err = abs(fd - g[i]) / max(1.0, abs(fd))
        worst = max(worst, err)
        checked += 1
        # the C++ drew the same index and got the same difference quotient
        fc = rec["fd_checks"][trial]
        @test i == fc["i"] + 1
        @test isapprox(fd, fc["fd"]; rtol = 1e-6, atol = 1e-9)
    end
    @test checked == 12
    @test worst < 1e-4

    # a short maximize_range run against the C++ end-to-end numbers
    o2 = K.RangeOptOptions(rounds = 2, iters_per_round = 5)
    rr = K.maximize_range(cs.c, X0, Phi, o2)
    mr = rec["maximize_range_2x5"]
    @test rr.rounds_run == mr["rounds_run"]
    @test rr.n_active == mr["n_active"]
    @test isapprox(rr.theta_ref_before, mr["theta_ref_before"]; rtol = 1e-10)
    @test isapprox(rr.theta_closed_before, mr["theta_closed_before"]; rtol = 1e-10)
    @test isapprox(rr.theta_ref_after, mr["theta_ref_after"]; rtol = 1e-6)
    @test isapprox(rr.theta_closed_after, mr["theta_closed_after"]; rtol = 1e-6)
    @test length(rr.trace) == length(mr["trace"])
end

# ---------------------------------------------------------------------------
# T4.2" (derivations/core.md): min-over-roots is not Theta_max.

@testset "T4.2\": the hexagon graze -- first contact is not the first overlap" begin
    cs = m1_tiling("hexagons_3_0")
    @test !K.has_collision(cs.c, K.deploy(cs.c, cs.X, 0.0).Y, 1e-12)
    B = K.deploy_basis(cs.c, cs.X)
    sd = K.swept_discs(cs.c, B)
    pairs = K.candidate_pairs(cs.c, sd, pi, true)

    theta_1 = K.exact_theta_max(cs.c, B, pairs).theta_max   # min over roots
    ov = K.exact_theta_max_overlap(cs.c, B, pairs, 1e-9, pi, 1e-9)

    # The two disagree by exactly pi/3: the vertex lands on an edge ENDPOINT, the faces
    # touch at a point and separate again.
    @test isapprox(theta_1, pi / 3; rtol = 1e-6)
    @test isapprox(ov.theta_max, 2pi / 3; rtol = 1e-6)
    @test ov.theta_max > theta_1 + 1.0

    # No two faces overlap anywhere strictly between the graze and Theta_max.
    for th in (1.05, 1.2, 1.5, 1.8, 2.0, 2.09)
        @test !K.has_collision(cs.c, K.deploy(cs.c, cs.X, th).Y, 1e-12)
    end
    # and they do just past it.
    @test K.has_collision(cs.c, K.deploy(cs.c, cs.X, ov.theta_max + 1e-3).Y, 1e-12)
    @info "hexagons: |C| = $(length(ov.candidates)), theta_1 = $theta_1, Theta_max = $(ov.theta_max)"
end

@testset "T4.2\": Theta_max agrees with bisection and is pruning-invariant" begin
    for name in ("squares_checker", "triangles_checker", "kagome_checker", "trunc_3_0", "snub_3_2")
        cs = m1_tiling(name)
        fx = M1["tilings"][name]
        K.has_collision(cs.c, K.deploy(cs.c, cs.X, 0.0).Y, 1e-12) && continue
        B = K.deploy_basis(cs.c, cs.X)
        sd = K.swept_discs(cs.c, B)
        all = K.candidate_pairs(cs.c, sd, pi, false)
        pruned = K.candidate_pairs(cs.c, sd, pi, true)
        ta = K.exact_theta_max_overlap(cs.c, B, all, 1e-9, pi, 1e-9).theta_max
        tp = K.exact_theta_max_overlap(cs.c, B, pruned, 1e-9, pi, 1e-9).theta_max
        @test isapprox(ta, tp; rtol = 1e-12)  # the broad phase is sound
        @test abs(ta - bisect_theta_max(cs.c, cs.X, 1e-12)) < 1e-5
        # Theta_max is never earlier than the first contact.
        @test ta >= K.exact_theta_max(cs.c, B, pruned).theta_max - 1e-12
        @test isapprox(ta, fx["theta_overlap_all"]; rtol = 1e-12, atol = 1e-13)
    end
end

@testset "T4.2\": the candidate set is complete -- no overlap change off C(X)" begin
    cs = m1_tiling("snub_3_2")
    if !K.has_collision(cs.c, K.deploy(cs.c, cs.X, 0.0).Y, 1e-12)
        B = K.deploy_basis(cs.c, cs.X)
        sd = K.swept_discs(cs.c, B)
        ov = K.exact_theta_max_overlap(cs.c, B, K.candidate_pairs(cs.c, sd, pi, true), 1e-9, pi, 1e-9)
        # Corollary T4.2': the overlap status is locally constant off C(X). Sample each gap
        # at three interior points and require one verdict per gap.
        lo = 0.0
        for i in 0:length(ov.candidates)
            hi = i < length(ov.candidates) ? ov.candidates[i + 1] : Float64(pi)
            if hi - lo < 1e-6
                lo = hi
                continue
            end
            first = K.has_collision(cs.c, K.deploy(cs.c, cs.X, lo + 0.25 * (hi - lo)).Y, 1e-12)
            for f in (0.5, 0.75)
                @test K.has_collision(cs.c, K.deploy(cs.c, cs.X, lo + f * (hi - lo)).Y, 1e-12) == first
            end
            lo = hi
        end
        @info "snub square: |C| = $(length(ov.candidates)) gaps checked for constancy"
    end
end

@testset "T4.5b: the swept radius -- absolute frame needs sqrt(2), face frame does not" begin
    cs = m1_tiling("snub_3_0")
    B = K.deploy_basis(cs.c, cs.X)
    sd = K.swept_discs(cs.c, B)

    @testset "absolute frame: max(|C|,|S|) is violated, by at most sqrt(2)" begin
        # y(theta) = cos(theta/2) C + sin(theta/2) S is a genuine ellipse about the origin,
        # so max_theta |y| = sigma_max([C|S]), which max(|C|,|S|) under-estimates. This is
        # the 68% / ratio 1.4141 measurement of derivations/core.md T4.5b.
        violations = 0
        total = 0
        worst_ratio = 0.0
        for pv in 1:K.n_prime(B)
            C = K.basis_c(B, pv); S = K.basis_s(B, pv)
            rmax = max(norm(C), norm(S))
            rsound = sqrt(dot(C, C) + dot(S, S))
            peak = 0.0
            for i in 0:400
                t = 0.5 * pi * i / 400
                peak = max(peak, norm(cos(t) * C + sin(t) * S))
            end
            total += 1
            @test peak <= rsound * (1 + 1e-12)  # sqrt(|C|^2+|S|^2) always bounds it (T4.4)
            peak > rmax * (1 + 1e-9) && (violations += 1)
            rmax > 0 && (worst_ratio = max(worst_ratio, peak / rmax))
        end
        @test violations > 0                        # the spec's bound really is violated
        @test worst_ratio <= sqrt(2.0) + 1e-9       # and never by more than sqrt(2) (T4.4)
        @info "T4.5b absolute frame: $violations/$total copies exceed max(|C|,|S|), worst ratio $worst_ratio"
    end

    @testset "face-local frame: the swept region is exactly the flat circumdisc" begin
        # Relative to the moving face centroid, C - gc = x_v - xbar_f and
        # S - gs = -sigma_f J (x_v - xbar_f) (the per-face translation u_f cancels), so the
        # trajectory is cos t (x_v - xbar) - sigma sin t J(x_v - xbar) = R(-sigma t)(x_v - xbar):
        # a CIRCLE of radius |x_v - xbar_f|. Hence in this frame max(|C|,|S|) = sigma_max
        # exactly and the sqrt(2) correction does not bite. swept_discs() uses this frame,
        # which is why K2c measures max_f rho_f / r_f = 1 to 1e-12.
        worst_ratio = 0.0
        worst_iso = 0.0
        for f in 1:K.n_faces(cs.m)
            gc = K.Vec2(sd.gc[f, 1], sd.gc[f, 2]); gs = K.Vec2(sd.gs[f, 1], sd.gs[f, 2])
            for pv in cs.c.prime_faces[f]
                C = K.basis_c(B, pv) - gc; S = K.basis_s(B, pv) - gs
                worst_iso = max(worst_iso, abs(norm(C) - norm(S)) / max(1e-300, norm(C)))
                peak = 0.0
                for i in 0:400
                    t = 0.5 * pi * i / 400
                    peak = max(peak, norm(cos(t) * C + sin(t) * S))
                end
                norm(C) > 0 && (worst_ratio = max(worst_ratio, peak / norm(C)))
            end
            @test isapprox(sd.rho_max[f], sd.circum[f]; rtol = 1e-12)
            @test isapprox(sd.rho[f], sqrt(2.0) * sd.circum[f]; rtol = 1e-12)
        end
        @test worst_iso < 1e-12            # |C| == |S| in the face frame
        @test worst_ratio < 1 + 1e-9       # the trajectory never leaves the circumdisc
        @info "T4.5b face frame: worst | |C|-|S| | / |C| = $worst_iso, worst peak/r_f = $worst_ratio"
    end
end

# ---------------------------------------------------------------------------
# The tau = 0 deflation and the validity certificate (core.md T5.2b.0/T5.2b.2,
# check.md R2.5 / R3.3).

@testset "T5.2b.2: the three harmonic classes are what the algebra says" begin
    # class 1: C = p + q != 0
    @test K.classify_harmonic(K.Harmonic(1.0, 0.3, 0.2)) == K.Class1
    # class 2: C = 0, B = 2r != 0  -> g = tau (A tau + B), simple root at theta = 0
    @test K.classify_harmonic(K.Harmonic(1.0, -1.0, 0.5)) == K.Class2
    # class 3: C = 0, B = 0, A != 0 -> h = p (1 - cos theta)
    @test K.classify_harmonic(K.Harmonic(0.4330127, -0.4330127, 0.0)) == K.Class3
    # identically zero
    @test K.classify_harmonic(K.Harmonic(0.0, 0.0, 0.0)) == K.Zero

    @testset "class 3 has constant sign on (0, pi) and is never a contact (Lemma T5.1e)" begin
        h = K.Harmonic(0.4330127, -0.4330127, 0.0)
        for i in 1:399
            th = pi * i / 400.0
            @test K.harmonic_eval(h, th) > 0
            # doctest Approx.epsilon(1e-14) is |a-b| <= 1e-14 (1 + max|a|,|b|): absolute near 0
            @test isapprox(K.harmonic_eval(h, th), h.p * (1 - cos(th)); rtol = 1e-14, atol = 1e-14)
        end
        @test isempty(K.harmonic_roots_deflated(h, 0.0, pi))
        @test K.harmonic_no_root_in(h, 0.2)
    end

    @testset "class 2: the deflated root is T5.3's theta* = 2 atan(-r/p)" begin
        # p = -q so A = 2p, B = 2r, tau* = -B/A = -r/p
        for p in (-1.0, -0.4, -2.5), r in (0.2, 1.0, 3.0)
            h = K.Harmonic(p, -p, r)
            rts = K.harmonic_roots_deflated(h, 0.0, pi)
            @test length(rts) == 1
            length(rts) == 1 || continue
            @test isapprox(rts[1], 2.0 * atan(-r / p); rtol = 1e-12)
            @test abs(K.harmonic_eval(h, rts[1])) < 1e-12 * K.scale(h)
            @test isapprox(K.harmonic_eval(h, 0.0), 0.0; atol = 1e-14)  # theta = 0 is a root
        end
    end

    @testset "class 1 falls through to the amplitude/phase roots unchanged" begin
        rng = K.MT19937(4242)
        for i in 1:2000
            h = K.Harmonic(K.uniform_real(rng, -2.0, 2.0), K.uniform_real(rng, -2.0, 2.0),
                           K.uniform_real(rng, -2.0, 2.0))
            K.classify_harmonic(h) != K.Class1 && continue
            a = K.harmonic_roots(h, 0.0, pi)
            b = K.harmonic_roots_deflated(h, 0.0, pi)
            @test length(a) == length(b)
            length(a) == length(b) || continue
            for k in eachindex(a)
                @test isapprox(a[k], b[k]; rtol = 1e-12)
            end
        end
    end
end

@testset "T5.2b.2: deflated roots agree with a direct sign-change scan" begin
    # Truth: a genuine CROSSING of h on (lo, hi), found by scanning, with no reference to
    # the tau chart. Skip harmonics whose root sits within 1e-7 of an endpoint (ambiguous).
    function crossings(h, lo, hi, n)
        # Strict sign CHANGE between two nonzero samples. h(1e-9) is exactly 0.0 in double
        # for a class-3 harmonic (1 - cos(1e-9) underflows), so "prev == 0" must not count
        # as a crossing -- class 3 has constant sign and no crossing at all.
        out = Float64[]
        prev = K.harmonic_eval(h, lo)
        for i in 1:n
            th = lo + (hi - lo) * i / n
            cur = K.harmonic_eval(h, th)
            if prev != 0.0 && cur != 0.0 && (prev < 0) != (cur < 0)
                a = lo + (hi - lo) * (i - 1) / n; b = th
                for _ in 1:60
                    mid = 0.5 * (a + b)
                    if (K.harmonic_eval(h, a) < 0) != (K.harmonic_eval(h, mid) < 0)
                        b = mid
                    else
                        a = mid
                    end
                end
                push!(out, 0.5 * (a + b))
            end
            prev = cur
        end
        return out
    end
    rng = K.MT19937(99991)
    U() = K.uniform_real(rng, -1.0, 1.0)
    tested = 0; mism_deflated = 0; mism_raw = 0; c2 = 0; c3 = 0
    for i in 0:59999
        kind = i % 3
        # The three classes exactly as the algebra defines them. A class-3 harmonic with
        # p + q perturbed off zero is NOT tested: there the remaining root sits at
        # theta ~ sqrt(noise/p) and its existence is decided by rounding, so direct root
        # finding is equally ill-posed and there is no truth to compare against. That is
        # the class derivations/check.md R2.5 excludes explicitly, and for the same reason.
        local h
        if kind == 0
            h = K.Harmonic(U(), U(), U())                       # generic
        elseif kind == 1
            p = U(); h = K.Harmonic(p, -p, U())                 # class 2
        else
            p = U(); h = K.Harmonic(p, -p, 0.0)                 # class 3
        end
        K.scale(h) < 1e-6 && continue
        truth = crossings(h, 1e-9, pi, 6000)
        ambiguous = any(t -> t < 1e-7 || t > pi - 1e-7, truth)
        ambiguous && continue
        tested += 1
        k = K.classify_harmonic(h)
        c2 += (k == K.Class2)
        c3 += (k == K.Class3)
        length(K.harmonic_roots_deflated(h, 1e-9, pi)) != length(truth) && (mism_deflated += 1)
        length(K.harmonic_roots(h, 1e-9, pi)) != length(truth) && (mism_raw += 1)
    end
    @test tested > 10000
    @test c2 > 1000
    @test c3 > 1000
    @test mism_deflated == 0  # the deflated list is exact on all three classes
    # On synthetic coefficients of size O(1) the raw list happens to agree too: the
    # tau = 0 root rounds to ~1e-16, below the 1e-9 lower limit. It is on REAL geometry,
    # where the coefficients come from differences with cancellation, that the raw list
    # reports it at ~1e-8 and inside the interval -- measured in the certificate test
    # below (255 spurious roots on 4 certified reference patterns).
    @info "deflation: $tested harmonics ($c2 class 2, $c3 class 3), mismatches deflated $mism_deflated vs raw $mism_raw"
end

const M1_CERT_FIELDS = (:pos, :nooverlap, :noroot, :n_inverted, :n_pairs, :n_candidates,
                        :n_identically_zero, :n_class1, :n_class2, :n_class3,
                        :n_roots_deflated, :n_roots_inadmissible)
function m1_check_cert(cert, fx)
    for fld in M1_CERT_FIELDS
        @test getfield(cert, fld) == fx[string(fld)]
    end
    # n_roots_undeflated counts the tau = 0 artefact, a root at ~1e-16 whose side of 0 is
    # decided by the last ulp of atan2/acos/hypot: exact on the arm64 libm only.
    if M1_EXACT_TIES
        @test cert.n_roots_undeflated == fx["n_roots_undeflated"]
    else
        @test abs(cert.n_roots_undeflated - fx["n_roots_undeflated"]) <= 0.25 * fx["n_roots_undeflated"] + 5
    end
    @test isapprox(cert.eps, fx["eps"]; rtol = 1e-15)
    @test isapprox(cert.theta_1, fx["theta_1"]; rtol = 1e-15)
    @test isapprox(cert.min_signed_area, fx["min_signed_area"]; rtol = 1e-12)
    @test isapprox(cert.first_root, fx["first_root"]; rtol = 1e-12, atol = 1e-13)
    # 0-based M'-vertex ids, -1 = none
    @test cert.bad_pv == (fx["bad_pv"] < 0 ? 0 : fx["bad_pv"] + 1)
    @test cert.bad_a == (fx["bad_a"] < 0 ? 0 : fx["bad_a"] + 1)
    @test cert.bad_b == (fx["bad_b"] < 0 ? 0 : fx["bad_b"] + 1)
end

@testset "certificate: POS and NOOVERLAP(eps/2) and NOROOT imply Theta_max >= eps" begin
    eps = 0.006
    held = 0; checked = 0
    undeflated_total = 0
    for name in ("squares_checker", "triangles_checker", "kagome_checker", "hexagons_3_0", "snub_3_2")
        cs = m1_tiling(name)
        fx = M1["tilings"][name]
        cert = K.validity_certificate(cs.c, cs.X, eps)
        m1_check_cert(cert, fx["cert_0006_conv"])
        # Permanent incidences must not reach the root test. On these patterns the hinge
        # cases are removed earlier still, by the `p == a || p == b` guard: a hinge keeps the
        # SOURCE vertex, so both faces carry the same M'-vertex there and the pair is never
        # formed. The identity test catches whatever is left (accidental permanent
        # collinearity); it may legitimately fire zero times.
        @test cert.n_class1 + cert.n_class2 + cert.n_class3 + cert.n_identically_zero == cert.n_candidates
        checked += 1
        K.valid(cert) || continue
        held += 1
        # With the deflation there is no root in (0, eps); without it the same scan reports
        # the tau = 0 artefact as a root. This is the measurement, on real geometry.
        @test cert.n_roots_deflated == 0
        undeflated_total += cert.n_roots_undeflated
        # The conclusion: no interior overlap anywhere on (0, eps), by an independent scan.
        for i in 1:46
            @test !K.has_collision(cs.c, K.deploy(cs.c, cs.X, eps * i / 47.0).Y, 1e-12)
        end
        # and therefore Theta_max >= eps
        B = K.deploy_basis(cs.c, cs.X)
        sd = K.swept_discs(cs.c, B)
        @test K.exact_theta_max_overlap(cs.c, B, K.candidate_pairs(cs.c, sd, pi, true), 1e-9, pi, 1e-9).theta_max >= eps
        @info "certificate holds on $name: $(cert.n_identically_zero) permanent incidences struck, classes $(cert.n_class1)/$(cert.n_class2)/$(cert.n_class3), undeflated would report $(cert.n_roots_undeflated) roots in (0,eps)"
    end
    @test checked == 5
    @test held > 0              # the hypothesis has content
    @test undeflated_total > 0  # and the deflation is doing real work on real geometry
    @info "without the tau = 0 deflation these certified cases would report $undeflated_total spurious roots in (0, eps)"
end

# ---------------------------------------------------------------------------
# Regression (F32): NOROOT must apply the two interval (projection) inequalities of
# T4.1b, not merely find a root of the orientation harmonic.
#
# A root of h_o,pi says the vertex is collinear with the INFINITE LINE through the edge.
# The C-list of Corollary T4.2' also demands 0 <= <w-a, b-a> <= |b-a|^2 -- the vertex on
# the SEGMENT -- and only such a root can change the overlap status. Without that filter
# the certificate rejected 1622 of 1859 jittered authored designs whose exact Theta_max
# was 1.0-2.4 rad (results/kill/jitter/cert_diagnosis.md); on the five rows diagnosed
# there, ALL 168 roots reported in (0, eps) had the vertex off the segment.
#
# The three cases below are the authored reference tilings that fail: 4.8.8, snub square
# and 3.4.3.12, at their unjittered positions. Each has exact Theta_max > 1.5 rad, so
# NOROOT(0.006) must hold; before the fix it did not.
@testset "F32: NOROOT counts only roots with the vertex ON the edge segment" begin
    eps = 0.006
    inadmissible_total = 0
    checked = 0
    for name in ("trunc_4_0", "snub_3_2", "t3_4_3_12_4_2")
        cs = m1_tiling(name)
        fx = M1["tilings"][name]
        B = K.deploy_basis(cs.c, cs.X)
        sd = K.swept_discs(cs.c, B)
        pairs = K.candidate_pairs(cs.c, sd, pi, true)
        theta_max = K.exact_theta_max_overlap(cs.c, B, pairs, 1e-9, pi, 1e-9).theta_max
        @test (theta_max < 1.0) == (fx["theta_overlap_pruned"] < 1.0)
        theta_max < 1.0 && continue  # not one of the rows this regression is about
        checked += 1
        cert = K.validity_certificate(cs.c, B, cs.X, pairs, eps)
        m1_check_cert(cert, fx["cert_0006"])
        # The geometry says the range is huge, so every clause of the certificate must hold.
        @test cert.pos
        @test cert.nooverlap
        @test cert.noroot
        @test cert.n_roots_deflated == 0
        # and the artefact the filter removes is real, not hypothetical
        inadmissible_total += cert.n_roots_inadmissible
        # An admissible root would have to be a contact angle of the exact scan; there is
        # none below eps, since the first one is theta_max itself.
        ca = K.contact_angles(cs.c, B, pairs, 0.0, eps)
        @test isempty(ca)
        @test length(ca) == length(fx["contact_angles_0_0006"])
        @info "F32 $name: Theta_max = $theta_max, NOROOT($eps) = $(cert.noroot), off-segment roots removed = $(cert.n_roots_inadmissible)"
    end
    # The 4.8.8 case gets its sigma from this file's own relaxation stream, which differs
    # from kill_common's, and does not always clear 1 rad; the other two always do.
    @test checked >= 2
    @test inadmissible_total > 0  # without the interval test these would be rejections
end

# ---------------------------------------------------------------------------
# Regression (K6): ValidityCertificate::first_root must be the SMALLEST admissible
# deflated root over EVERY candidate, not the first one the pair loop happens to hit.
# Callers (kill_k6) use it as the supremum of the eps for which NOROOT(eps) holds --
# taking a later root there certifies a range the structure does not have.

@testset "harmonic roots are returned in ascending order" begin
    # rts[1] in the certificate scan is only the smallest root of that candidate
    # because harmonic_roots sorts. Pin that contract: p + q cos + r sin with two roots
    # in (0, eps), whose natural arctan branches (phi - psi, phi + psi) come out in the
    # opposite order unless the routine sorts.
    t1 = 0.4; t2 = 1.1
    # Build h with roots exactly at t1 and t2: h = A * (cos(t - phi) - cos(psi)) with
    # phi = (t1 + t2)/2, psi = (t2 - t1)/2.
    phi = 0.5 * (t1 + t2); psi = 0.5 * (t2 - t1)
    h = K.Harmonic(-cos(psi), cos(phi), sin(phi))
    @test abs(K.harmonic_eval(h, t1)) < 1e-12
    @test abs(K.harmonic_eval(h, t2)) < 1e-12
    rts = K.harmonic_roots(h, 0.0, 2.0)
    @test length(rts) >= 2
    @test issorted(rts)
    @test isapprox(rts[1], t1; rtol = 1e-9)
    dfl = K.harmonic_roots_deflated(h, 0.0, 2.0)
    @test issorted(dfl)
end

const M1_SIX = ("squares_checker", "triangles_checker", "kagome_checker", "hexagons_3_0",
                "snub_3_2", "t3_4_3_12_4_2")

@testset "certificate: NOROOT holds just below first_root, by more than the root tolerance" begin
    # The "largest certified eps" is the supremum of the eps for which NOROOT(eps) holds,
    # i.e. first_root probed from below. harmonic_roots accepts a root up to eps + 1e-15,
    # so the step down must CLEAR that tolerance: nextafter(first_root, 0) is one ulp
    # (~2e-16 at first_root ~ 1) and does NOT, which silently makes the largest certified
    # eps identically 0 on every design that has a root at all.
    # eps must be wide enough that at least two of the four patterns carry an ADMISSIBLE
    # root (F32: off-segment collinearities no longer count), hence 3.0 and not 1.2.
    eps = 3.0
    checked = 0
    # The last two carry genuine contacts below eps; the regular ones mostly do not, now
    # that off-segment collinearities have stopped counting (F32).
    for name in M1_SIX
        cs = m1_tiling(name)
        fx = M1["tilings"][name]
        B = K.deploy_basis(cs.c, cs.X)
        sd = K.swept_discs(cs.c, B)
        pairs = K.candidate_pairs(cs.c, sd, pi, true)
        cert = K.validity_certificate(cs.c, B, cs.X, pairs, eps)
        m1_check_cert(cert, fx["cert_3"])
        cert.noroot && continue
        @test cert.first_root > 1e-6
        checked += 1
        # At the root itself, NOROOT must fail (the interval is half-open).
        @test !K.validity_certificate(cs.c, B, cs.X, pairs, cert.first_root).noroot
        @test fx["cert_3_at_root_noroot"] == false
        # One ulp below is NOT enough to clear the 1e-15 acceptance slack.
        @test !K.validity_certificate(cs.c, B, cs.X, pairs, prevfloat(cert.first_root)).noroot
        @test fx["cert_3_ulp_below_noroot"] == false
        # 1e-12 below is, and that is the step the K6 driver uses.
        @test K.validity_certificate(cs.c, B, cs.X, pairs, cert.first_root - 1e-12).noroot
        @test fx["cert_3_1e12_below_noroot"] == true
    end
    @test checked >= 2
end

@testset "certificate: first_root is the minimum over all candidates, not the first seen" begin
    # eps large enough that many candidates have admissible roots inside (0, eps), so the
    # pair loop meets a LARGER root before the smallest one on at least one pattern.
    eps = 3.0
    discriminating = 0; checked = 0
    for name in M1_SIX
        cs = m1_tiling(name)
        B = K.deploy_basis(cs.c, cs.X)
        sd = K.swept_discs(cs.c, B)
        pairs = K.candidate_pairs(cs.c, sd, pi, true)
        cert = K.validity_certificate(cs.c, B, cs.X, pairs, eps)
        cert.noroot && continue  # no root at all: nothing to compare
        checked += 1

        # Independent rescan in the SAME loop order, collecting every admissible deflated
        # root in (0, eps): the minimum, and the one a first-seen rule would have reported.
        min_root = -1.0; first_seen = -1.0
        PF = cs.c.prime_faces
        function rescan(fe, fv)
            E = PF[fe]; V = PF[fv]
            ne = length(E)
            for i in 1:ne
                a = E[i]; b = E[mod1(i + 1, ne)]
                U = K.basis_c(B, b) - K.basis_c(B, a); Vv = K.basis_s(B, b) - K.basis_s(B, a)
                for p in V
                    (p == a || p == b) && continue
                    P = K.basis_c(B, p) - K.basis_c(B, a); Q = K.basis_s(B, p) - K.basis_s(B, a)
                    det = K.orient_from_vectors(U, Vv, P, Q)
                    sc = max(K.scale(det), norm(U) * norm(P))
                    abs(det.p) + abs(det.q) + abs(det.r) <= 1e-11 * max(1e-300, sc) && continue
                    # Same admissibility filter the certificate applies (F32): a root counts
                    # only when the vertex is on the edge SEGMENT, not merely on its line.
                    L2 = K.dot_from_vectors(U, Vv, U, Vv)
                    D = K.dot_from_vectors(U, Vv, P, Q)
                    lscale = abs(L2.p) + K.amp(L2)
                    lscale <= 0 && continue
                    r0 = -1.0
                    for th in K.harmonic_roots_deflated(det, 0.0, eps)
                        sp = K.harmonic_eval(D, th); l2 = K.harmonic_eval(L2, th); tol = 1e-12 * lscale
                        (sp < -tol || sp > l2 + tol) && continue
                        r0 = th
                        break  # roots are ascending
                    end
                    r0 < 0 && continue
                    first_seen < 0 && (first_seen = r0)
                    (min_root < 0 || r0 < min_root) && (min_root = r0)
                end
            end
        end
        for (f, g) in pairs
            rescan(f, g)
            rescan(g, f)
        end
        @test min_root > 0
        @test isapprox(cert.first_root, min_root; rtol = 1e-12)
        first_seen > min_root * (1 + 1e-9) && (discriminating += 1)
        @info "$name: min root $min_root, first-seen root $first_seen"
    end
    @test checked >= 2
    # The bug is observable: on at least one pattern the first candidate with a root is
    # NOT the one carrying the smallest root, so the old code returned the wrong eps.
    @test discriminating >= 1
end
