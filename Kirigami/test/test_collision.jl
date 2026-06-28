# test_collision.jl -- port of code/tests/test_collision.cpp, case by case.
#
# The polygon-predicate cases are literal. The tiling cases read the meshes/sigmas the C++
# drew from its case-local std::mt19937 and the C++ results (theta_max, min beta, the
# collision-sweep ladder) from CORPUS/reference_patterns/test_fixtures_collision.json
# (frozen by data/corpus/reference_patterns/freeze_fixtures_2a.cpp).
include("helpers.jl")

const K = Kirigami
# TODO(generators): the relaxation-dependent blocks (theta_max, velocity, collision_opt) are
# read from the fixture; regenerate via the call in each block's "provenance" field once
# assign_orientation_relaxation reproduces the C++ stream bit-exactly.
const COL_FIXTURES = JSON.parsefile(joinpath(CORPUS, "reference_patterns", "test_fixtures_collision.json"))

poly(pts...) = [K.Vec2(p[1], p[2]) for p in pts]

@testset "polygon overlap primitive" begin
    A = poly((0, 0), (1, 0), (1, 1), (0, 1))
    B = poly((0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5))
    C = poly((1, 0), (2, 0), (2, 1), (1, 1))        # edge-adjacent
    D = poly((2, 2), (3, 2), (3, 3), (2, 3))        # disjoint
    E = poly((0.2, 0.2), (0.8, 0.2), (0.5, 0.8))    # contained
    @test K.polygons_overlap(A, B)
    @test !K.polygons_overlap(A, C)  # touching along an edge is not an overlap
    @test !K.polygons_overlap(A, D)
    @test K.polygons_overlap(A, E)
    @test K.polygons_overlap(E, A)
end

@testset "theta_max: geometric first collision never exceeds min_i beta_i" begin
    cases = 0
    for fam in COL_FIXTURES["theta_max"]["families"]
        g = fixture_mesh_raw(fam["mesh"])
        c = K.make_cut(g)
        hs = K.holes_partition(c)
        sys = K.assemble_system(c, hs, g.X, K.Fixed)
        rep = K.solve_system(sys, g.X)
        X0 = K.matrix_to_points(rep.X0)
        @test maximum(norm.(X0 .- fixture_points(fam["X0"]))) < 1e-9
        K.hole_residuals(c, X0, hs).max_norm > 1e-9 && continue
        tm = K.theta_max(c, X0, 90, 40)
        @test tm.theta_max_geometric <= tm.min_beta + 1e-4
        cases += 1
        # the C++ numbers for this kind (e.g. hexagons: 2.0944 / 2.0944 / n_split 3)
        @test K.n_split(c) == fam["n_split"]
        @test isapprox(tm.theta_max_geometric, fam["theta_max_geometric"]; atol = 1e-9)
        @test isapprox(tm.min_beta, fam["min_beta"]; atol = 1e-12)
        @test tm.collided == fam["collided"]
        @test maximum(abs.(tm.beta .- Float64.(fam["beta"]))) < 1e-12
    end
    @test cases >= 5
    @test cases == 7  # the C++ ran all seven kinds
end

@testset "rotating squares: theta_max equals min_i beta_i = pi" begin
    fx = COL_FIXTURES["rotating_squares"]
    g = K.tiling_squares(K.rect(K.Vec2(2.5, 2.5), 2.51, 2.51))
    @test same_mesh_as_fixture(g, fx["mesh"])
    g.sigma = checkerboard_sigma(g)
    @test g.sigma == Int[s for s in fx["mesh"]["orientation"]]
    c = K.make_cut(g)
    tm = K.theta_max(c, g.X, 180, 45)
    @test isapprox(tm.min_beta, pi; rtol = 1e-12)
    @test isapprox(tm.theta_max_geometric, pi; rtol = 1e-6)
    @test !tm.collided  # no overlap anywhere in (0, pi]
    @test isapprox(tm.min_beta, fx["min_beta"]; atol = 1e-14)
    @test tm.theta_max_geometric == fx["theta_max_geometric"]
end

@testset "Eq. (7): deployment velocity equals dY/dtheta at theta = 0" begin
    fx = COL_FIXTURES["velocity"]
    g = fixture_mesh_raw(fx["mesh"])
    c = K.make_cut(g)
    z = K.deployment_velocity(c, g.X)
    d = K.deploy(c, g.X, 0.0)
    @test length(z) == length(d.dY_dtheta)
    worst = maximum(norm.(z .- d.dY_dtheta))
    @test worst < 1e-12
    @test maximum(norm.(z .- fixture_points(fx["z"]))) < 1e-12
end

@testset "Eq. (9): collision optimization enlarges the deployment range" begin
    improved = 0
    tried = 0
    for fam in COL_FIXTURES["collision_opt"]["families"]
        g = fixture_mesh_raw(fam["mesh"])
        c = K.make_cut(g)
        @test K.n_split(c) == fam["n_split"]
        K.n_split(c) == 0 && continue
        hs = K.holes_partition(c)
        sys = K.assemble_system(c, hs, g.X, K.Fixed)
        rep = K.solve_system(sys, g.X)
        @test rep.dim_null == fam["dim_null"]
        rep.dim_null == 0 && continue
        X0 = K.matrix_to_points(rep.X0)
        K.hole_residuals(c, X0, hs).max_norm > 1e-9 && continue
        r = K.optimize_collision_sweep(c, X0, rep.Phi)
        tried += 1
        # the optimized embedding stays in the shape space
        @test K.hole_residuals(c, r.X_opt, hs).max_norm < 1e-8
        # the gamma ladder includes the T = 0 fallback, so it can never do worse
        @test r.theta_max_after >= r.theta_max_before - 1e-9
        sw = fam["sweep"]
        @test isapprox(r.theta_max_before, sw["theta_max_before"]; atol = 1e-9)
        # The optimizer's path (L-BFGS over a different null-space basis) is not expected
        # to reproduce Eigen's iterates; the C++ outcome is reported for comparison only.
        ladder = join(["($(g)->$(round(t; digits = 6)))" for (g, t) in r.ladder], " ")
        println(fam["kind"], ": theta_max ", r.theta_max_before, " -> ", r.theta_max_after,
                " at gamma = ", r.gamma_used, "; ladder: ", ladder,
                "\n    C++: -> ", sw["theta_max_after"], " at gamma = ", sw["gamma_used"],
                "; ladder: ", join(["($(g)->$(round(t; digits = 6)))" for (g, t) in sw["ladder"]], " "))
        r.theta_max_after > r.theta_max_before + 1e-6 && (improved += 1)
    end
    @test tried > 0
    @test tried == 3
    @test improved > 0
end

# ============================================================================
# F34 -- the hinge-vertex referee artefact (results/kill/b4 voronoi_93, faces 94/184;
# K9 delaunay ids 13, 40, 103, 190).  Two deployed prime faces that share exactly one
# M'-vertex (the hinge pin) and whose interiors are disjoint.  The predicate must report
# NO overlap here, at every tolerance, and must not depend on the tolerance value at all
# for this configuration.
# ============================================================================
@testset "polygons_overlap: faces sharing a hinge vertex with a reflex sector do not overlap" begin
    P1 = poly((35.522612034574969, 2.8976170473256473),
              (35.919893509036022, 3.1484087335745494),
              (36.013918812893856, 3.6292236769982966),
              (35.070313922385914, 3.0868107146716834),
              (34.512614462367168, 2.8219234379478291))
    P2 = poly((34.233106280612688, 1.0182292700260629),
              (34.512614462367168, 2.8219234379478291),
              (35.067637209025698, 3.0923743720689583),
              (34.742880186886559, 3.4379946577381681),
              (33.800478898505453, 3.2058972270991357),
              (32.960043232689365, 2.3447653598242049),
              (33.184898605051494, 2.3139255252053297),
              (34.239945838106621, 1.3526146778850854))
    # the shared pin is P1[5] == P2[2], bit-for-bit (both copies read the same Y[pv])
    @test norm(P1[5] - P2[2]) == 0.0
    for tol in (0.0, 1e-15, 1e-12, 1e-9, 1e-6)
        @test !K.polygons_overlap(P1, P2, tol)
        @test !K.polygons_overlap(P2, P1, tol)
    end
end

# A minimal synthetic version of the same configuration: two triangles meeting at one
# shared vertex, the second with a REFLEX corner there, interiors strictly disjoint.
# Also the positive control: nudge one of them across the pin and the overlap must fire.
@testset "polygons_overlap: shared-vertex touching vs. genuine overlap through the pin" begin
    pin = K.Vec2(1.0, 1.0)
    A = [pin, K.Vec2(2.0, 0.0), K.Vec2(3.0, 1.0), K.Vec2(2.0, 1.5)]
    B = [pin, K.Vec2(0.0, 1.4), K.Vec2(-1.0, 0.0), K.Vec2(0.2, 0.2)]  # reflex at (0.2,0.2)
    for tol in (0.0, 1e-12, 1e-9, 1e-6)
        @test !K.polygons_overlap(A, B, tol)
        @test !K.polygons_overlap(B, A, tol)
    end
    Bx = [p + K.Vec2(0.6, 0.0) for p in B]  # now the two interiors genuinely intersect
    @test K.polygons_overlap(A, Bx, 1e-12)
    @test K.polygons_overlap(Bx, A, 1e-12)
end

# Collinear / shared-edge degeneracies: an edge shared exactly (no overlap), a partially
# shared edge (no overlap), and a polygon whose vertices all lie on the other's boundary
# while the interiors do intersect (must fire).
@testset "polygons_overlap: collinear and boundary-only degeneracies" begin
    S = poly((0, 0), (1, 0), (1, 1), (0, 1))
    right = poly((1, 0), (2, 0), (2, 1), (1, 1))                    # exact shared edge
    partial = poly((1, 0.25), (2, 0.25), (2, 0.75), (1, 0.75))      # sub-edge
    diag = poly((0, 0), (1, 1), (0, 1))                             # half of S
    for tol in (0.0, 1e-12, 1e-9)
        @test !K.polygons_overlap(S, right, tol)
        @test !K.polygons_overlap(S, partial, tol)
        @test K.polygons_overlap(S, diag, tol)   # interiors intersect, all vertices on dS
        @test K.polygons_overlap(diag, S, tol)
    end
end

# Two faces adjacent across a SPLIT edge, at theta = 0, where the two copies of that edge
# coincide exactly, ~43 units from the origin (K5 voronoi_75, faces 106 and 390).
@testset "polygons_overlap: faces meeting along a coincident split edge far from the origin" begin
    A = poly((25.157593658847794, 34.013688249821989),
             (23.525630356123127, 33.741539519199222),
             (23.435629113495715, 32.737836500697981),
             (23.438280868613827, 32.732383750554867),
             (23.647595416965657, 32.664869500090575),
             (24.478156409438640, 32.756084473098767))
    B = poly((23.435629113495715, 32.737836500697981),
             (22.077319787908216, 33.409498398804359),
             (21.901162010051728, 33.169673228066642),
             (21.913041261020524, 32.789008667540571),
             (22.819656910734363, 31.980928766465244),
             (23.438280868613827, 32.732383750554867))
    # A[3] == B[1] and A[4] == B[6]: the same segment, traversed in opposite directions
    @test norm(A[3] - B[1]) == 0.0
    @test norm(A[4] - B[6]) == 0.0
    for tol in (0.0, 1e-15, 1e-12, 1e-9, 1e-6)
        @test !K.polygons_overlap(A, B, tol)
        @test !K.polygons_overlap(B, A, tol)
    end
    # translating the whole configuration must not change the answer
    A2 = [p + K.Vec2(1000.0, -500.0) for p in A]
    B2 = [p + K.Vec2(1000.0, -500.0) for p in B]
    @test !K.polygons_overlap(A2, B2, 1e-12)
end
