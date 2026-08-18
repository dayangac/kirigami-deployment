# test_kinematics.jl -- port of code/tests/test_kinematics.cpp, case by case.
#
# Inputs (meshes, sigmas, the shuffled face orders and the theta samples drawn from the
# case-local std::mt19937) and the C++ inline results (Eq. (6) embeddings X, mismatches,
# the finite-difference tally) were frozen by replaying the C++ test bodies
# (data/corpus/reference_patterns/freeze_fixtures_2a.cpp) into
# CORPUS/reference_patterns/test_fixtures_kinematics.json.
include("helpers.jl")

const K = Kirigami
# TODO(generators): regenerate each block via the call in its "provenance" field
# (generate(kind, {R}, rng) + assign_orientation_relaxation / random_sigma with the stated
# seed) once generators.jl + orientation.jl reproduce the C++ streams.
const KIN_FIXTURES = JSON.parsefile(joinpath(CORPUS, "reference_patterns", "test_fixtures_kinematics.json"))

# A deployable embedding of `g` with the given sigma (projection of Eq. 6).
function deployable_embedding(g::K.Mesh, c::K.CutStructure, hs::K.HoleSet)
    sys = K.assemble_system(c, hs, g.X, K.Fixed)
    rep = K.solve_system(sys, g.X)
    return K.matrix_to_points(rep.X0)
end

signed_angle(a::K.Vec2, b::K.Vec2) = atan(a[1] * b[2] - a[2] * b[1], dot(a, b))

@testset "forward kinematics: every hinge opens by exactly theta" begin
    for fam in KIN_FIXTURES["hinge_opens"]["families"]
        g = fixture_mesh_raw(fam["mesh"])
        c = K.make_cut(g)
        hs = K.holes_partition(c)
        X = deployable_embedding(g, c, hs)
        @test K.hole_residuals(c, X, hs).max_norm < 1e-9
        # the port's Eq. (6) projection equals the C++ one
        @test maximum(norm.(X .- fixture_points(fam["X"]))) < 1e-9
        for (ti, th) in enumerate((0.1, 0.5, 1.0))
            d = K.deploy(c, X, th)
            @test d.max_mismatch < 1e-9
            @test isapprox(d.max_mismatch, fam["max_mismatch_0.1_0.5_1.0"][ti]; atol = 1e-12)
            for e in c.hinge_edges
                ed = g.edges[e]
                f1 = g.half_edges[ed.he[1]].face
                f2 = g.half_edges[ed.he[2]].face
                src = c.hinge_dir[e].src
                dst = c.hinge_dir[e].dst
                p = K.prime_vertex(c, f1, src)
                @test p == K.prime_vertex(c, f2, src)  # hinge stays at the SOURCE
                d1 = d.Y[K.prime_vertex(c, f1, dst)] - d.Y[p]
                d2 = d.Y[K.prime_vertex(c, f2, dst)] - d.Y[p]
                @test abs(abs(signed_angle(d1, d2)) - th) < 1e-9
            end
        end
    end
end

@testset "forward kinematics is independent of the BFS order on a deployable structure" begin
    for fam in KIN_FIXTURES["bfs_order"]["families"]
        g = fixture_mesh_raw(fam["mesh"])
        c = K.make_cut(g)
        hs = K.holes_partition(c)
        X = deployable_embedding(g, c, hs)
        @test K.hole_residuals(c, X, hs).max_norm < 1e-9
        ref = K.deploy(c, X, 0.7, 1)
        @test maximum(norm.(ref.Y .- fixture_points(fam["Y_ref"]))) < 1e-9
        # from the C++ X the port's deploy() is BIT-exact (FMA-contracted 2x2 mat-vec, Apple
        # __sincos_stret for the rotation entries)
        @test K.deploy(c, fixture_points(fam["X"]), 0.7, 1).Y == fixture_points(fam["Y_ref"])
        for (trial, order0) in enumerate(fam["orders"])
            order = [Int(f) + 1 for f in order0]  # 0-based C++ face order -> 1-based
            alt = K.deploy_with_order(c, X, 0.7, order)
            @test alt.max_mismatch < 1e-9
            @test K.rigid_align_residual(alt.Y, ref.Y) < 1e-9
            mm, al = fam["mismatch_and_align"][trial]
            @test isapprox(alt.max_mismatch, mm; atol = 1e-12)
            @test isapprox(K.rigid_align_residual(alt.Y, ref.Y), al; atol = 1e-12)
        end
    end
end

@testset "non-deployable (3,4,3,12) has a strictly positive kinematic mismatch" begin
    fx = KIN_FIXTURES["t3_4_3_12_r3_5"]
    g = fixture_mesh_raw(fx["mesh"])
    c = K.make_cut(g)
    hs = K.holes_partition(c)
    r = K.hole_residuals(c, g.X, hs)
    @test r.max_norm > 1e-3  # Eq. (2) is violated (2026 Fig. 21b)
    @test isapprox(r.max_norm, fx["max_norm"]; rtol = 1e-12)
    d = K.deploy(c, g.X, 0.3)
    @test d.max_mismatch > 1e-6
    @test isapprox(d.max_mismatch, fx["deploy_0_3_mismatch"]; rtol = 1e-9)
    # After the Eq. (6) projection it becomes deployable and the mismatch vanishes.
    X0 = deployable_embedding(g, c, hs)
    @test maximum(norm.(X0 .- fixture_points(fx["X0"]))) < 1e-9
    @test K.hole_residuals(c, X0, hs).max_norm < 1e-9
    @test K.deploy(c, X0, 0.3).max_mismatch < 1e-9
end

@testset "dY/dtheta matches central finite differences (>= 1000 samples)" begin
    fx = KIN_FIXTURES["fd"]
    samples = 0
    worst_rel = 0.0
    for fam in fx["families"]
        g = fixture_mesh_raw(fam["mesh"])
        for trial in fam["trials"]
            g.sigma = Int[s for s in trial["sigma"]]
            c = K.make_cut(g)
            for th in trial["thetas"]
                h = 1e-6
                d = K.deploy(c, g.X, th)
                dp = K.deploy(c, g.X, th + h)
                dm = K.deploy(c, g.X, th - h)
                for i in 1:c.n_prime_vertices
                    fd = (dp.Y[i] - dm.Y[i]) / (2h)
                    denom = max(1.0, norm(fd))
                    worst_rel = max(worst_rel, norm(fd - d.dY_dtheta[i]) / denom)
                    samples += 1
                end
            end
        end
    end
    @test samples >= 1000
    @test worst_rel < 1e-6
    # C++: 20532 samples, worst relative error 1.08751e-09 (finite-difference noise level)
    @test samples == fx["samples"]
    @test worst_rel < 10 * fx["worst_rel"]
end

@testset "Remark A.4: split-cut duplicates stay parallel throughout deployment" begin
    fx = KIN_FIXTURES["remark_A4"]
    edges_checked = 0
    for fam in fx["families"]
        g = fixture_mesh_raw(fam["mesh"])
        c = K.make_cut(g)
        @test K.n_split(c) == fam["n_split"]
        K.n_split(c) == 0 && continue
        hs = K.holes_partition(c)
        X = deployable_embedding(g, c, hs)
        @test K.hole_residuals(c, X, hs).max_norm < 1e-9
        @test maximum(norm.(X .- fixture_points(fam["X"]))) < 1e-9
        tm_geo = K.theta_max(c, X).theta_max_geometric
        @test isapprox(tm_geo, fam["theta_max_geometric"]; atol = 1e-9)
        tmax = max(0.4, tm_geo)
        for th in (0.1, 0.3, 0.6, 0.9 * tmax)
            d = K.deploy(c, X, th)
            for e in c.split_edges
                ed = g.edges[e]
                f1 = g.half_edges[ed.he[1]].face
                f2 = g.half_edges[ed.he[2]].face
                u = d.Y[K.prime_vertex(c, f1, ed.key.b)] - d.Y[K.prime_vertex(c, f1, ed.key.a)]
                v = d.Y[K.prime_vertex(c, f2, ed.key.b)] - d.Y[K.prime_vertex(c, f2, ed.key.a)]
                @test abs(u[1] * v[2] - u[2] * v[1]) < 1e-9
                edges_checked += 1
            end
        end
    end
    @test edges_checked > 0
    @test edges_checked == fx["edges_checked"]  # C++: 76 (split edge, theta) pairs
end
