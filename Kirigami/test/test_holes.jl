# test_holes.jl -- port of code/tests/test_holes.cpp, case by case.
#
# The C++ cases build their meshes with generators.hpp + std::mt19937 random sigmas and,
# for the geometric cross-check, deploy them through the Tutte auxetic system / FK /
# collision units. Those sequences were frozen by replaying the C++ tests exactly into
# CORPUS/reference_patterns/test_fixtures_holes.json: per graph the mesh with its sigma,
# and -- where the C++ reached the geometric stage -- the deployed M'-vertex positions Yd
# and the C++ holes_geometric output `geo` (0-based edge ids). The checks below are the
# C++ CHECKs applied to the frozen data; where the C++ skipped a graph (no deployable,
# overlap-free embedding) Yd is null and the Julia test skips it too.
include("helpers.jl")
import JSON

const K = Kirigami
# TODO(generators): regenerate via generators.jl + tutte_auxetic/kinematics/collision once
# those are ported; each fixture block carries a "provenance" record (test case, generator
# call + args, seed, sigma method, how X/Yd were produced). See data/corpus/README.md.
const HOLES_FIXTURES = JSON.parsefile(joinpath(CORPUS, "reference_patterns", "test_fixtures_holes.json"))

holes_fixture_mesh(j) = K.mesh_from_json_string(JSON.json(j))
holes_fixture_sigma(s) = Int[x for x in s]
holes_fixture_points(P) = [K.Vec2(Float64(p[1]), Float64(p[2])) for p in P]
# C++ edge ids are 0-based in the fixture -> 1-based
holes_fixture_geo(g) = sort!([Int[e + 1 for e in cyc] for cyc in g])

# The two independent hole-preimage constructions must agree, must partition
# E_hinge union E_split, and (on a deployable embedding) must reproduce the
# bounded complement components of M' traced geometrically.
@testset "hole preimages: seed-growing == partition formulation == geometry" begin
    # TODO(generators): rng = MT19937(2024); see fx["provenance"] for the exact sequence.
    fx = HOLES_FIXTURES["case1"]
    graphs = 0
    geometric_checks = 0
    split_cycles = 0
    disconnected = 0
    cycle_mismatch = 0
    for rec in fx["graphs"]
        g = holes_fixture_mesh(rec["mesh"])
        c = K.make_cut(g)
        a = K.holes_seed_growing(c)
        b = K.holes_partition(c)
        same, msg = K.same_hole_sets(a, b)
        @test same
        ok, multi, uncov = K.holes_partition_edges(c, b)
        @test ok
        @test isempty(multi)
        @test isempty(uncov)
        graphs += 1
        # Remark A.4's tree claim is conditional on the connectivity principle of
        # Sec. 4.2; a random sigma can violate it, and then Def. 4.2 stops matching
        # the geometry (measured below).
        forest = K.split_subgraph_is_forest(c)[1]
        connected = K.count_components(c) == 1
        @test forest == rec["forest"]
        @test connected == rec["connected"]
        forest || (split_cycles += 1)
        connected || (disconnected += 1)

        # geometric cross-check on the frozen overlap-free deployment
        rec["Yd"] === nothing && continue
        Yd = holes_fixture_points(rec["Yd"])
        @test length(Yd) == c.n_prime_vertices
        geo = K.holes_geometric(c, Yd)
        @test geo == holes_fixture_geo(rec["geo"])
        comb = sort([b.all[i].edges for i in b.interior_indices])
        if forest && connected
            @test geo == comb
            geometric_checks += 1
        elseif geo != comb
            cycle_mismatch += 1
        end
    end
    @test graphs >= 50
    # Every geometric mismatch we saw came with a split cycle or a disconnected M'.
    @test cycle_mismatch <= split_cycles + disconnected
    # the frozen C++ run's own tallies
    cnt = fx["counts"]
    @test graphs == cnt["graphs"]
    @test geometric_checks == cnt["geometric_checks"]
    @test split_cycles == cnt["split_cycles"]
    @test disconnected == cnt["disconnected"]
    @test cycle_mismatch == cnt["cycle_mismatch"]
    @info "combinatorial hole checks: $graphs graphs, $geometric_checks geometric agreements; " *
          "$split_cycles with split cycles, $disconnected disconnected, " *
          "$cycle_mismatch geometric mismatches (all of them in that set)"
end

@testset "hole preimages of the rotating-squares pattern are the interior vertices" begin
    # TODO(generators): tiling_squares(rect(Vec2(2, 2), 2.01, 2.01)) + checkerboard_sigma.
    m = holes_fixture_mesh(HOLES_FIXTURES["rotating_squares"])  # checkerboard sigma frozen
    c = K.make_cut(m)
    @test K.n_split(c) == 0
    hs = K.holes_partition(c)
    @test K.n_interior_holes(hs) == K.n_interior_vertices(m)
    for i in hs.interior_indices
        @test length(hs.all[i].vertices) == 1
        # each interior vertex of a square grid has 2 hinge edges pointing into it
        @test length(hs.all[i].edges) == 2
    end
end

@testset "split cuts merge holes: H = #interior - #interior split edges (forest case)" begin
    # TODO(generators): rng = MT19937(99); generate(kind, [3.2], rng) + 8 x random_sigma per kind
    # (see HOLES_FIXTURES["case3_provenance"]).
    for fam in HOLES_FIXTURES["case3"]
        g = holes_fixture_mesh(fam["mesh"])
        for s in fam["sigmas"]
            g.sigma = holes_fixture_sigma(s)
            c = K.make_cut(g)
            K.split_subgraph_is_forest(c)[1] || continue
            K.count_components(c) == 1 || continue
            interior_split = count(c.split_edges) do e
                k = g.edges[e].key
                !g.vertex_is_boundary[k.a] && !g.vertex_is_boundary[k.b]
            end
            # Only exact when no split component mixes interior and boundary vertices.
            hs = K.holes_partition(c)
            mixed = any(hs.all) do h
                has_b = any(v -> g.vertex_is_boundary[v], h.vertices)
                has_i = any(v -> !g.vertex_is_boundary[v], h.vertices)
                has_i && has_b
            end
            mixed && continue
            @test K.n_interior_holes(hs) == K.n_interior_vertices(g) - interior_split
        end
    end
end

# The spec's cross-check: combinatorial vs geometric hole detection on >= 50 random
# graphs. Two populations are used:
#  (A) random topology (Delaunay / Voronoi / quad-dominant) with the paper's own
#      orientation assignment (Sec. 4.2). For these the Eq. (6) projection almost
#      always lands on a SELF-INTERSECTING embedding (measured below), on which the
#      geometric trace is not defined -- this is the paper's own open limitation.
#  (B) random embeddings inside the shape space of deployable tilings: a random
#      affine map (Remark 4.1) plus a random null-space offset (Eq. 5).
# The embeddings X and the overlap-free deployments Yd come frozen from the C++ run.
@testset "hole preimages: geometric agreement on >= 50 random instances" begin
    # TODO(generators): rng = MT19937(31337); see fx["provenance"] for populations (A) and (B).
    fx = HOLES_FIXTURES["case4"]
    checked = 0
    for rec in fx["instances"]
        g = holes_fixture_mesh(rec["mesh"])
        c = K.make_cut(g)
        @test K.check_remark_A1(c)[1]
        a = K.holes_seed_growing(c)
        b = K.holes_partition(c)
        @test K.same_hole_sets(a, b)[1]
        @test K.holes_partition_edges(c, b)[1]
        rec["Yd"] === nothing && continue
        # the C++ only deployed forest && connected instances
        @test K.split_subgraph_is_forest(c)[1]
        @test K.count_components(c) == 1
        Yd = holes_fixture_points(rec["Yd"])
        geo = K.holes_geometric(c, Yd)
        @test geo == holes_fixture_geo(rec["geo"])
        comb = sort([b.all[i].edges for i in b.interior_indices])
        @test geo == comb
        # the cycle variant bounds the same holes: one M'-vertex loop per hole, whose
        # original edges are exactly the hole's edge set
        cycles = K.holes_geometric_cycles(c, Yd)
        @test length(cycles) == length(geo)
        checked += 1
    end
    @test checked >= 50
    @test checked == fx["counts"]["checked"]
    @info "geometric hole agreement verified on $checked frozen instances; " *
          "$(fx["counts"]["random_topology_checked"]) of $(fx["counts"]["random_topology_total"]) " *
          "random-topology graphs had a non-self-intersecting Tutte auxetic embedding; " *
          "$(fx["counts"]["skipped_overlap"]) instances skipped (self-overlapping deployment)"
end
