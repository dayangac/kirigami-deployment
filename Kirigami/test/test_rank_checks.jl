# test_rank_checks.jl -- core/rank_checks, case by case.
# Three structural measurements on L that the paper never states:
#   1. the row sum 1^T L,
#   2. the factorization L = R D and the out-harmonic left null space,
#   3. the hinge graph Gamma and the Euler-type hole count.
#
# The 12 bounded and 3 torus cases (meshes + sigma) and every report field of the
# reference run are frozen in CORPUS/reference_patterns/test_fixtures_rank_checks.json
# (provenance in data/corpus/README.md). Torus faces are loaded AS
# STORED (fixture_mesh_raw): their wrap-around faces are geometrically degenerate.
include("helpers.jl")

const K = Kirigami
# The 4 relaxation cases and the 4 delaunay + random_sigma cases are read from the
# fixture; they regenerate from one MT19937(20260903) shared over the relaxation + delaunay
# calls in file order (see the fixture's "provenance").
const RANK_FIXTURES = JSON.parsefile(joinpath(CORPUS, "reference_patterns", "test_fixtures_rank_checks.json"))

struct RankCase
    name::String
    mesh::K.Mesh
    mode::K.BoundaryMode
    boundary_free::Bool  # torus: every vertex interior
    fx::Any              # the frozen reference reports
end
# The RNG-free cases are built by the generators (checked bit-identical to the frozen
# mesh, checkerboard sigma recomputed and compared); the relaxation / random-sigma
# cases come from the fixture.
const GENERATED_CASES = Dict(
    "rotating_squares" => () -> K.tiling_squares(K.rect(K.Vec2(2.5, 2.5), 2.51, 2.51)),
    "triangles_alternating" => () -> K.tiling_triangles(K.disk(K.Vec2(0.13, 0.07), 3.5)),
    "kagome_3636" => () -> K.tiling_kagome(K.disk(K.Vec2(0.13, 0.07), 3.0)),
    "periodic_squares_4x4" => () -> K.periodic_squares(4, 4),
    "torus_squares_4x4" => () -> K.torus_squares(4, 4),
    "torus_squares_6x4" => () -> K.torus_squares(6, 4),
    "torus_triangles_4x4" => () -> K.torus_triangles(4, 4),
)
function case_mesh(j)
    haskey(GENERATED_CASES, j["name"]) || return fixture_mesh_raw(j["mesh"])
    g = GENERATED_CASES[j["name"]]()
    same_mesh_as_fixture(g, j["mesh"]) || error("generator mesh differs from the frozen mesh: $(j["name"])")
    g.sigma = checkerboard_sigma(g)
    g.sigma == Int[s for s in j["mesh"]["orientation"]] || error("checkerboard sigma differs from the frozen one: $(j["name"])")
    return g
end
function load_cases(key::String)
    modes = Dict("Fixed" => K.Fixed, "Periodic" => K.Periodic, "None" => K.None)
    return [RankCase(j["name"], case_mesh(j), modes[j["mode"]], j["boundary_free"], j)
            for j in RANK_FIXTURES[key]]
end
bounded_cases() = load_cases("bounded_cases")
torus_cases() = load_cases("torus_cases")

struct Bundle
    c::K.CutStructure
    hs::K.HoleSet
    sys::K.LinearSystem
end
function build(g::K.Mesh, mode::K.BoundaryMode)
    c = K.make_cut(g)
    hs = K.holes_partition(c)
    return Bundle(c, hs, K.assemble_system(c, hs, g.X, mode))
end

# every scalar field of the frozen reference report equals the fresh one
function check_fields(rep, fx, fields; float_atol = 1e-12)
    for f in fields
        v = getfield(rep, Symbol(f))
        if v isa AbstractFloat
            @test isapprox(v, fx[f]; atol = float_atol)
        else
            @test v == fx[f]
        end
    end
end
const ROW_SUM_FIELDS = ["N", "H", "all_vertices_interior", "degree_identity", "n_degree_mismatch",
                        "restricted_degree_identity", "n_restricted_mismatch", "n_hinge",
                        "n_hinge_in_L", "no_notch_hinges", "support_on_boundary",
                        "n_interior_nonzero", "max_abs_r", "r_is_zero"]
const FACT_FIELDS = ["H", "n_hinge", "N", "L_equals_RD", "max_abs_diff", "Z_computed", "rank_L",
                     "dim_Z", "rank_identity", "n_harmonic_checks", "n_harmonic_violations",
                     "out_harmonic"]
const HINGE_FIELDS = ["n_faces", "n_hinge", "c_gamma", "H", "H_all", "n_notches", "predicted",
                      "identity_holds", "identity_with_notches", "defect"]

# ------------------------------------------------------------------ check 1 --

@testset "check 1: the row sum of L is the NOTCH-RESTRICTED degree difference" begin
    # The identity r_v == indeg_hinge(v) - outdeg_hinge(v) is FALSE on every patch
    # with a boundary: L omits the boundary-touching preimages (notches), so the
    # hinge edges those preimages own contribute nothing to the row sum. The
    # corrected identity restricts both degrees to hinge edges owned by hole rows.
    naive_ok = 0
    total = 0
    for cs in bounded_cases()
        @testset "$(cs.name)" begin
            b = build(cs.mesh, cs.mode)
            r = K.check_row_sum(b.c, b.hs, b.sys)
            @test length(r.r) == K.n_vertices(cs.mesh)
            @test r.restricted_degree_identity
            @test r.n_restricted_mismatch == 0
            # the naive form holds exactly when no hinge edge belongs to a notch
            @test r.degree_identity == r.no_notch_hinges
            # and the support of r is confined to boundary vertices only in that case
            @test r.support_on_boundary == (r.n_interior_nonzero == 0)
            naive_ok += r.degree_identity
            total += 1
            check_fields(r, cs.fx["row_sum"], ROW_SUM_FIELDS)
            @test r.r == Float64.(cs.fx["row_sum"]["r"])
        end
    end
    @test total == 12
    @test naive_ok == 0  # reference: naive identity holds on 0 / 12 bounded patches
end

@testset "check 1: the torus patches are boundary free, so r == 0 and rank(L) == H - 1" begin
    for cs in torus_cases()
        @testset "$(cs.name)" begin
            b = build(cs.mesh, cs.mode)
            @test K.n_split(b.c) == 0
            r = K.check_row_sum(b.c, b.hs, b.sys)
            @test r.all_vertices_interior
            @test r.r_is_zero
            @test r.max_abs_r == 0.0
            @test r.degree_identity       # no notches, so the naive form holds too
            @test r.no_notch_hinges
            @test r.restricted_degree_identity
            # every vertex carries a hole, and the rows sum to zero, so L drops one rank
            H = K.n_interior_holes(b.hs)
            @test H == K.n_vertices(cs.mesh)
            @test b.sys.n_boundary_rows == 0
            rep = K.solve_system(b.sys, cs.mesh.X)
            @test rep.rank_L == H - 1
            @test rep.rank_L == cs.fx["solve_rank_L"]
            check_fields(r, cs.fx["row_sum"], ROW_SUM_FIELDS)
        end
    end
end

@testset "check 1: the generator 'periodic' patches still have boundary vertices" begin
    # Eqs. (3b)-(3c) are imposed as constraints on a finite patch whose boundary
    # edges are NOT topologically identified, so r is not zero there. This records
    # the distinction the rank statement of Sec. 4.4 turns on.
    cs = only(filter(c -> c.name == "periodic_squares_4x4", bounded_cases()))
    b = build(cs.mesh, K.Periodic)
    r = K.check_row_sum(b.c, b.hs, b.sys)
    @test !r.all_vertices_interior
    @test !r.r_is_zero
    @test r.restricted_degree_identity
    # notches exist here, so the naive identity and the boundary-support claim fail
    @test !r.no_notch_hinges
    @test !r.degree_identity
    @test !r.support_on_boundary
end

# ------------------------------------------------------------------ check 2 --

@testset "check 2: L == R * D exactly" begin
    for cs in bounded_cases()
        @testset "$(cs.name)" begin
            b = build(cs.mesh, cs.mode)
            f = K.check_factorization(b.c, b.hs, b.sys, false)
            @test f.max_abs_diff == 0.0
            @test f.L_equals_RD
            @test size(f.R, 1) == b.sys.n_hole_rows
            @test size(f.R, 2) == K.n_hinge(b.c)
            @test size(f.D, 1) == K.n_hinge(b.c)
            @test size(f.D, 2) == K.n_vertices(cs.mesh)
        end
    end
    for cs in torus_cases()
        @testset "$(cs.name)" begin
            b = build(cs.mesh, cs.mode)
            f = K.check_factorization(b.c, b.hs, b.sys, false)
            @test f.L_equals_RD
        end
    end
end

@testset "check 2: rank(L) == H - dim Z and every y in Z is out-harmonic" begin
    function run(cases)
        for cs in cases
            @testset "$(cs.name)" begin
                b = build(cs.mesh, cs.mode)
                f = K.check_factorization(b.c, b.hs, b.sys, true)
                @test f.Z_computed
                @test f.rank_identity
                @test f.rank_L + f.dim_Z == f.H
                @test f.max_harmonic_residual < 1e-9
                @test f.n_harmonic_violations == 0
                @test f.out_harmonic
                # cross-check the rank against the solver's own SVD of L
                rep = K.solve_system(b.sys, cs.mesh.X)
                @test rep.rank_L == f.rank_L
                check_fields(f, cs.fx["factorization"], FACT_FIELDS)
                @test f.max_harmonic_residual < 10 * cs.fx["factorization"]["max_harmonic_residual"] + 1e-14
            end
        end
    end
    run(bounded_cases())
    run(torus_cases())
end

@testset "check 2: the torus left null space is exactly the constant vector" begin
    for cs in torus_cases()
        @testset "$(cs.name)" begin
            b = build(cs.mesh, cs.mode)
            f = K.check_factorization(b.c, b.hs, b.sys)
            @test f.dim_Z == 1
            # a unit constant vector: every entry has the same magnitude
            y = f.Z[:, 1]
            for i in eachindex(y)
                @test abs(abs(y[i]) - abs(y[1])) < 1e-12
            end
            # constant g is trivially out-harmonic because indeg == outdeg (Remark A.1)
            @test f.out_harmonic
        end
    end
end

# ------------------------------------------------------------------ check 3 --

@testset "check 3: hinge graph Euler count H == |E_hinge| - |F| + c(Gamma)" begin
    for cs in bounded_cases()
        @testset "$(cs.name)" begin
            b = build(cs.mesh, cs.mode)
            h = K.check_hinge_graph(b.c, b.hs)
            @test h.c_gamma >= 1
            # The identity holds as stated, with H counting only the all-interior
            # preimages: the boundary-touching notches are NOT bounded faces of Gamma.
            @test h.identity_holds
            @test h.H == h.n_hinge - h.n_faces + h.c_gamma
            h.n_notches > 0 && @test !h.identity_with_notches
            check_fields(h, cs.fx["hinge_graph"], HINGE_FIELDS)
        end
    end
end

@testset "check 3: on the torus the identity carries a genus correction" begin
    # Gamma is embedded on a torus (chi = 0) rather than the sphere (chi = 2), so
    # the plane count overshoots by exactly one.
    for cs in torus_cases()
        @testset "$(cs.name)" begin
            b = build(cs.mesh, cs.mode)
            h = K.check_hinge_graph(b.c, b.hs)
            @test h.c_gamma == 1
            @test h.n_notches == 0
            @test h.H == h.predicted - 1
            check_fields(h, cs.fx["hinge_graph"], HINGE_FIELDS)
        end
    end
end

@testset "check 3: split-free bounded patches -- Gamma is connected and H matches" begin
    # With no split cuts every preimage is a single vertex star, so the notches are
    # exactly the boundary vertices carrying an incoming hinge edge.
    cs = only(filter(c -> c.name == "rotating_squares", bounded_cases()))
    b = build(cs.mesh, K.Fixed)
    @test K.n_split(b.c) == 0
    h = K.check_hinge_graph(b.c, b.hs)
    @test h.c_gamma == 1
    @test h.H == K.n_interior_vertices(cs.mesh)
    @test h.H_all == length(b.hs.all)
    @test h.predicted == h.n_hinge - h.n_faces + 1
end
