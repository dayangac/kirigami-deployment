# test_system.jl -- port of code/tests/test_system.cpp, case by case.
#
# The C++ cases build their meshes with generators.hpp and (for some) sigma from
# assign_orientation_relaxation with a case-local std::mt19937. Those inputs, plus every
# number the C++ computes inline (ranks, residuals, X0, the random affine maps / null-space
# offsets and the resulting "worst" values), were frozen by replaying the C++ test bodies
# (data/corpus/reference_patterns/freeze_fixtures_2a.cpp, linked against the C++ library)
# into CORPUS/reference_patterns/test_fixtures_system.json. The checks below are the C++
# CHECKs, and where the C++ printed a MESSAGE its number is asserted against the port.
include("helpers.jl")

const K = Kirigami
# TODO(generators): the RNG-dependent blocks (t3_4_3_12, nullspace, orientation_small_patch)
# are still read from the fixture; regenerate via the call recorded in each block's
# "provenance" field once assign_orientation_relaxation reproduces the C++ stream bit-exactly.
const SYSTEM_FIXTURES = JSON.parsefile(joinpath(CORPUS, "reference_patterns", "test_fixtures_system.json"))

# Deterministic (RNG-free) cases call the Julia generators directly and assert that the
# result is bit-identical to the frozen C++ mesh; the checkerboard sigma is recomputed and
# compared to the frozen one too. The RNG / relaxation-dependent cases stay on the fixtures.
function checkerboard_case(g::K.Mesh, fx_mesh)
    @test same_mesh_as_fixture(g, fx_mesh)
    g.sigma = checkerboard_sigma(g)
    @test g.sigma == Int[s for s in fx_mesh["orientation"]]
    return g
end

struct Analysis
    c::K.CutStructure
    hs::K.HoleSet
    sys::K.LinearSystem
    rep::K.SolveReport
end
function analyze(g::K.Mesh, mode::K.BoundaryMode = K.Fixed)
    c = K.make_cut(g)
    hs = K.holes_partition(c)
    sys = K.assemble_system(c, hs, g.X, mode)
    rep = K.solve_system(sys, g.X)
    return Analysis(c, hs, sys, rep)
end
# the C++ solve numbers frozen next to each mesh
function check_solve_matches(a::Analysis, s)
    @test a.rep.rank_full == s["rank_full"]
    @test a.rep.rank_L == s["rank_L"]
    @test a.rep.dim_null == s["dim_null"]
    @test a.sys.n_hole_rows == s["n_hole_rows"]
    @test a.sys.n_boundary_rows == s["n_boundary_rows"]
    @test K.n_split(a.c) == s["n_split"]
    @test K.n_hinge(a.c) == s["n_hinge"]
    @test K.n_interior_holes(a.hs) == s["n_interior_holes"]
    @test a.rep.projection_ok == s["projection_ok"]
    @test isapprox(a.rep.sv_tol, s["sv_tol"]; rtol = 1e-9)
    X0_cpp = fixture_matrix(s["X0"])
    @test maximum(abs, a.rep.X0 - X0_cpp) < 1e-9
end

@testset "rotating squares: deployable, unique solution, one hole per interior vertex" begin
    fx = SYSTEM_FIXTURES["rotating_squares"]
    g = checkerboard_case(K.tiling_squares(K.rect(K.Vec2(2.5, 2.5), 2.51, 2.51)), fx["mesh"])
    @test K.n_faces(g) == 25
    a = analyze(g)
    @test K.n_split(a.c) == 0
    @test K.count_components(a.c) == 1
    @test K.count_components(a.c) == fx["components"]
    @test K.n_interior_holes(a.hs) == K.n_interior_vertices(g)
    @test K.deployable(K.hole_residuals(a.c, g.X, a.hs))
    @test a.rep.rank_L == K.n_interior_holes(a.hs)  # rows of L are independent
    @test a.rep.dim_null == 0                       # unique solution
    d = K.deploy(a.c, g.X, 0.6)
    @test d.max_mismatch < 1e-12
    @test isapprox(d.max_mismatch, fx["deploy_0_6_max_mismatch"]; atol = 1e-12)
    check_solve_matches(a, fx["solve"])
end

@testset "equilateral triangle tiling with alternating sigma is deployable with dim_null 0" begin
    fx = SYSTEM_FIXTURES["triangles_alternating"]
    g = checkerboard_case(K.tiling_triangles(K.disk(K.Vec2(0.13, 0.07), 3.5)), fx["mesh"])
    a = analyze(g)
    @test K.n_split(a.c) == 0  # the dual of a triangular tiling is 2-colourable
    @test K.deployable(K.hole_residuals(a.c, g.X, a.hs))
    @test a.rep.dim_null == 0
    @test K.n_interior_holes(a.hs) == K.n_interior_vertices(g)
    check_solve_matches(a, fx["solve"])
end

@testset "kagome (3.6.3.6) with alternating sigma is deployable with dim_null 0" begin
    fx = SYSTEM_FIXTURES["kagome_alternating"]
    g = checkerboard_case(K.tiling_kagome(K.disk(K.Vec2(0.13, 0.07), 3.0)), fx["mesh"])
    a = analyze(g)
    @test K.n_split(a.c) == 0
    @test K.deployable(K.hole_residuals(a.c, g.X, a.hs))
    @test a.rep.dim_null == 0
    check_solve_matches(a, fx["solve"])
end

@testset "hexagonal tiling: every connected sigma satisfies dim_null == #interior - H" begin
    fx = SYSTEM_FIXTURES["hexagons_exhaustive"]
    # TODO(generators): tiling_hexagons(disk(Vec2(0.13, 0.07), 2.4)) -- same faces, but the
    # Julia vertex coordinates differ from the C++ by 1 ulp (FMA contraction), so the
    # frozen mesh is used to keep the numbers identical.
    g = fixture_mesh_raw(fx["mesh"])
    @test K.n_faces(g) <= 20
    @test K.n_faces(g) >= 6
    F = K.n_faces(g)
    @test F == fx["F"]
    connected = 0; violations = 0; with_split = 0
    for mask in 0:(UInt64(1) << F)-1
        (mask & 1) != 0 && continue  # global flip symmetry
        g.sigma = [((mask >> (i - 1)) & 1) != 0 ? 1 : -1 for i in 1:F]
        c = K.make_cut(g)
        K.count_components(c) != 1 && continue
        connected += 1
        K.n_split(c) > 0 && (with_split += 1)
        hs = K.holes_partition(c)
        @test K.same_hole_sets(hs, K.holes_seed_growing(c))[1]
        @test K.holes_partition_edges(c, hs)[1]
        sys = K.assemble_system(c, hs, g.X, K.Fixed)
        rep = K.solve_system(sys, g.X)
        rep.dim_null != K.n_interior_vertices(g) - K.n_interior_holes(hs) && (violations += 1)
        # The paper's claim (Sec. 4.4): full rank exactly when E_split is empty.
        K.n_split(c) == 0 && @test rep.dim_null == 0
    end
    @test connected > 0
    @test with_split > 0
    @test violations == 0
    # the C++ tallies: F=7: 39 connected sigma assignments (39 with split cuts), 0 violations
    @test connected == fx["connected"]
    @test with_split == fx["with_split"]
    @test violations == fx["violations"]
end

@testset "(3,4,3,12): not deployable as the Euclidean tiling, deployable after Eq. (6)" begin
    fx = SYSTEM_FIXTURES["t3_4_3_12"]
    g = fixture_mesh_raw(fx["mesh"])
    a = analyze(g)
    r = K.hole_residuals(a.c, g.X, a.hs)
    @test !K.deployable(r)
    violating = count(v -> norm(v) > 1e-6, r.per_hole)
    @test violating > 0
    @test a.rep.projection_ok
    X0 = K.matrix_to_points(a.rep.X0)
    @test K.hole_residuals(a.c, X0, a.hs).max_norm < 1e-9
    @test K.deploy(a.c, X0, 0.4).max_mismatch < 1e-9
    # C++: 12 of 16 hole preimages violate Eq. (2), max residual 0.517638
    @test violating == fx["violating"]
    @test length(r.per_hole) == fx["n_per_hole"]
    @test isapprox(r.max_norm, fx["max_norm"]; rtol = 1e-12)
    @test isapprox(r.l2_norm, fx["l2_norm"]; rtol = 1e-12)
    per_hole_cpp = fixture_points(fx["per_hole"])
    @test maximum(norm.(r.per_hole .- per_hole_cpp)) < 1e-12
    check_solve_matches(a, fx["solve"])
end

@testset "Remark 4.1: uniform deployability is affine invariant" begin
    fx = SYSTEM_FIXTURES["affine"]
    g = checkerboard_case(K.tiling_kagome(K.disk(K.Vec2(0.13, 0.07), 3.0)),
                          SYSTEM_FIXTURES["kagome_alternating"]["mesh"])
    a = analyze(g)
    X0 = K.matrix_to_points(a.rep.X0)
    @test K.hole_residuals(a.c, X0, a.hs).max_norm < 1e-12
    worst = 0.0
    # the 100 maps [a00, a01, a10, a11, bx, by] the C++ drew from mt19937(19), U(-1.5, 1.5)
    for mp in fx["maps"]
        A = K.Mat2(mp[1], mp[3], mp[2], mp[4])  # column-major
        b = K.Vec2(mp[5], mp[6])
        Y = [A * x + b for x in X0]
        scale = maximum(norm, Y)
        worst = max(worst, K.hole_residuals(a.c, Y, a.hs).max_norm / max(scale, 1e-12))
    end
    @test worst < 1e-12
    @test worst < 1e-14   # C++: 4.11271e-16 (rounding-noise level); fx["worst"] holds it
end

@testset "null space: X0 + phi_i t^T stays deployable" begin
    fx = SYSTEM_FIXTURES["nullspace"]
    g = fixture_mesh_raw(fx["mesh"])
    a = analyze(g)
    @test a.rep.dim_null > 0
    @test a.rep.dim_null == fx["solve"]["dim_null"]   # C++: 5 basis vectors
    X0 = K.matrix_to_points(a.rep.X0)
    @test K.hole_residuals(a.c, X0, a.hs).max_norm < 1e-9
    worst = 0.0
    # offsets[i][t] = [ux, uy], the C++ draws from mt19937(23), U(-4, 4). The null-space
    # basis of the port need not equal Eigen's column for column, so the offsets are
    # applied to the port's Phi; the property is basis independent.
    for i in 1:a.rep.dim_null, (ux, uy) in fx["offsets"][i]
        X = copy(a.rep.X0)
        X[:, 1] += a.rep.Phi[:, i] * ux
        X[:, 2] += a.rep.Phi[:, i] * uy
        Y = K.matrix_to_points(X)
        worst = max(worst, K.hole_residuals(a.c, Y, a.hs).max_norm)
        # the boundary conditions of Eq. (4) must also survive
        for v in 1:K.n_vertices(g)
            g.vertex_is_boundary[v] && (worst = max(worst, norm(Y[v] - X0[v])))
        end
    end
    @test worst < 1e-9
    @test worst < 1e-13  # C++: 3.49368e-15 over 500 random offsets (fx["worst"])
    # the two null-space bases span the same subspace (C++ Phi vs port Phi)
    Phi_cpp = fixture_matrix(fx["Phi"])
    @test size(Phi_cpp) == size(a.rep.Phi)
    @test opnorm(Phi_cpp * Phi_cpp' - a.rep.Phi * a.rep.Phi') < 1e-9
    check_solve_matches(a, fx["solve"])
end

@testset "periodic boundary conditions (Eq. 3b-3c)" begin
    fx = SYSTEM_FIXTURES["periodic_squares_4x4"]
    g = checkerboard_case(K.periodic_squares(4, 4), fx["mesh"])
    @test g.periodic.present
    @test !isempty(g.periodic.pairs_h)
    @test !isempty(g.periodic.pairs_v)
    a = analyze(g, K.Periodic)
    @test a.rep.projection_ok
    X0 = K.matrix_to_points(a.rep.X0)
    @test K.hole_residuals(a.c, X0, a.hs).max_norm < 1e-9
    for p in g.periodic.pairs_h
        @test norm(X0[p[1]] - X0[p[2]] - g.periodic.th) < 1e-9
    end
    for p in g.periodic.pairs_v
        @test norm(X0[p[1]] - X0[p[2]] - g.periodic.tv) < 1e-9
    end
    check_solve_matches(a, fx["solve"])
end

@testset "orientation assignment (Eq. 1) and brute force agree on a small patch" begin
    # TODO(generators): bf = brute_force_orientation(g); re = assign_orientation_relaxation(g,
    # MT19937(101), 10, 600, 180) -- both belong to orientation.jl; until then the two
    # C++ reports are read from the fixture and their statistics re-derived here.
    fx = SYSTEM_FIXTURES["orientation_small_patch"]
    g = fixture_mesh_raw(fx["mesh"])
    @test K.n_faces(g) <= 20
    bf = fx["brute_force"]
    re = fx["relaxation"]
    @test bf["components"] == 1
    @test re["components"] == 1
    # brute force minimizes the number of split cuts among connected assignments
    @test bf["n_split"] <= re["n_split"]
    # the reported statistics must match a direct recomputation
    g.sigma = Int[s for s in bf["sigma"]]
    c = K.make_cut(g)
    @test K.n_split(c) == bf["n_split"]
    @test K.n_hinge(c) == bf["n_hinge"]
    @test K.n_interior_holes(K.holes_partition(c)) == bf["n_holes"]
    @test K.count_components(c) == bf["components"]
    g.sigma = Int[s for s in re["sigma"]]
    c = K.make_cut(g)
    @test K.n_split(c) == re["n_split"]
    @test K.n_hinge(c) == re["n_hinge"]
    @test K.n_interior_holes(K.holes_partition(c)) == re["n_holes"]
    @test K.count_components(c) == re["components"]
end

@testset "JSON round trip preserves the graph, sigma and periodic data" begin
    fx = SYSTEM_FIXTURES["periodic_kagome_3x3"]
    g = checkerboard_case(K.periodic_kagome(3, 3), fx["mesh"])
    @test K.n_vertices(g) == fx["N"]
    @test K.n_faces(g) == fx["F"]
    @test K.n_edges(g) == fx["E"]
    h = K.mesh_from_json_string(K.mesh_to_json_string(g))
    @test K.n_vertices(h) == K.n_vertices(g)
    @test K.n_faces(h) == K.n_faces(g)
    @test K.n_edges(h) == K.n_edges(g)
    @test h.sigma == g.sigma
    @test h.periodic.present
    @test length(h.periodic.pairs_h) == length(g.periodic.pairs_h)
    @test length(g.periodic.pairs_h) == fx["n_pairs_h"]
    for i in 1:K.n_vertices(g)
        @test norm(h.X[i] - g.X[i]) < 1e-12
    end
end
