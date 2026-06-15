# test_mesh_cut.jl -- port of code/tests/test_mesh_cut.cpp, case by case.
#
# The C++ cases build their meshes with generators.hpp + std::mt19937 random sigmas.
# Those sequences were frozen by replaying the C++ tests exactly (same seeds, same
# call order) into CORPUS/reference_patterns/test_fixtures_mesh_cut.json; the checks
# below are the C++ CHECKs applied to the frozen meshes.
include("helpers.jl")
import JSON

const K = Kirigami
const MESH_CUT_FIXTURES = JSON.parsefile(joinpath(CORPUS, "reference_patterns", "test_fixtures_mesh_cut.json"))

# a mesh JSON object (0-based file format) -> Mesh, and a frozen sigma array
fixture_mesh(j) = K.mesh_from_json_string(JSON.json(j))
fixture_sigma(s) = Int[x for x in s]

@testset "half-edge direction follows sigma (2026 Fig. 4)" begin
    # Two unit squares sharing the edge (1,0)-(1,1); vertex/face layout as welded by
    # the C++ mesh_from_polygons (vertices in order of first appearance).
    function two_squares()
        m = K.Mesh([K.Vec2(0, 0), K.Vec2(1, 0), K.Vec2(1, 1), K.Vec2(0, 1), K.Vec2(2, 0), K.Vec2(2, 1)],
                   [[1, 2, 3, 4], [2, 5, 6, 3]])
        K.normalize_face_ccw!(m)
        K.build_topology!(m)
        return m
    end
    m = two_squares()
    @test K.n_faces(m) == 2
    shared = 0
    for e in 1:K.n_edges(m)
        m.edges[e].n_faces == 2 && (shared = e)
    end
    @test shared >= 1
    f1 = m.half_edges[m.edges[shared].he[1]].face
    f2 = m.half_edges[m.edges[shared].he[2]].face

    @testset "opposite sigma -> half-edges have the SAME direction (hinge cut)" begin
        m.sigma = fill(-1, 2)
        m.sigma[f2] = +1
        d1 = K.half_edge_direction(m, f1, shared)
        d2 = K.half_edge_direction(m, f2, shared)
        @test d1[1] == d2[1]
        @test d1[2] == d2[2]
        c = K.make_cut(m)
        @test c.edge_type[shared] == K.Hinge
        @test K.n_hinge(c) == 1
        @test K.n_split(c) == 0
        # Hinge stays at the source; the target is duplicated.
        src = c.hinge_dir[shared].src
        dst = c.hinge_dir[shared].dst
        @test K.prime_vertex(c, f1, src) == K.prime_vertex(c, f2, src)
        @test K.prime_vertex(c, f1, dst) != K.prime_vertex(c, f2, dst)
    end
    @testset "equal sigma -> half-edges are antiparallel (split cut)" begin
        m.sigma = fill(-1, 2)
        d1 = K.half_edge_direction(m, f1, shared)
        d2 = K.half_edge_direction(m, f2, shared)
        @test d1[1] == d2[2]
        @test d1[2] == d2[1]
        c = K.make_cut(m)
        @test c.edge_type[shared] == K.Split
        # both endpoints duplicated
        ka = m.edges[shared].key.a
        kb = m.edges[shared].key.b
        @test K.prime_vertex(c, f1, ka) != K.prime_vertex(c, f2, ka)
        @test K.prime_vertex(c, f1, kb) != K.prime_vertex(c, f2, kb)
    end
end

@testset "non-manifold input is rejected" begin
    m = K.Mesh([K.Vec2(0, 0), K.Vec2(1, 0), K.Vec2(1, 1), K.Vec2(0, 1), K.Vec2(2, 0)],
               [[1, 2, 3, 4], [1, 2, 5], [2, 1, 3]])
    @test_throws ErrorException K.build_topology!(m)
end

@testset "Remark A.1: in-degree == out-degree at every interior vertex" begin
    fx = MESH_CUT_FIXTURES["remark_A1"]
    checked = 0
    for fam in fx["families"]
        base = fixture_mesh(fam["mesh"])
        for (trial, s) in enumerate(fam["sigmas"])
            base.sigma = fixture_sigma(s)
            c = K.make_cut(base)
            ok, bad = K.check_remark_A1(c)
            @test ok
            ok || @info "failed on $(fam["kind"]) trial $(trial - 1)"
            checked += 1
        end
    end
    for j in fx["delaunay"]
        g = fixture_mesh(j)  # sigma stored in the fixture's "orientation"
        @test !isempty(g.sigma)
        c = K.make_cut(g)
        @test K.check_remark_A1(c)[1]
        checked += 1
    end
    @test checked >= 200
end

@testset "M' vertex count matches the Remark A.2 / A.3 bookkeeping" begin
    # At an interior vertex with k hinge edges in (and k out) and s split edges,
    # the number of duplicated copies is k + s.
    fx = MESH_CUT_FIXTURES["A2A3_kagome"]
    g = fixture_mesh(fx["mesh"])
    for s in fx["sigmas"]
        g.sigma = fixture_sigma(s)
        c = K.make_cut(g)
        copies = zeros(Int, K.n_vertices(g))
        for pv in 1:c.n_prime_vertices
            copies[c.prime_to_original[pv]] += 1
        end
        for v in 1:K.n_vertices(g)
            g.vertex_is_boundary[v] && continue
            nsplit = count(e -> c.edge_type[e] == K.Split, g.vertex_edges[v])
            @test copies[v] == c.hinge_in[v] + nsplit
        end
    end
end

# Julia addition (no C++ counterpart): I/O boundary and geometry helpers are exercised
# implicitly by the C++ suite through the generators; here they need a direct check.
@testset "JSON round trip keeps the 0-based file format (Julia addition)" begin
    fx = MESH_CUT_FIXTURES["A2A3_kagome"]
    g = fixture_mesh(fx["mesh"])
    g.sigma = fixture_sigma(fx["sigmas"][1])
    j = JSON.parse(K.mesh_to_json_string(g))
    @test minimum(minimum.(j["faces"])) == 0
    @test j["faces"] == fx["mesh"]["faces"]
    g2 = K.mesh_from_json_string(K.mesh_to_json_string(g))
    @test g2.faces == g.faces
    @test g2.sigma == g.sigma
    @test K.n_edges(g2) == K.n_edges(g)
    @test all(isapprox.(g2.X, g.X; atol=0.0))
    # save/load through a file
    p = joinpath(mktempdir(), "m.json")
    K.save_mesh_json(g, p)
    g3 = K.load_mesh_json(p)
    @test g3.faces == g.faces && g3.sigma == g.sigma
    # face geometry helpers on a unit square
    sq = K.Mesh([K.Vec2(0, 0), K.Vec2(1, 0), K.Vec2(1, 1), K.Vec2(0, 1)], [[1, 2, 3, 4]])
    K.build_topology!(sq)
    @test isapprox(K.face_signed_area(sq, 1), 1.0; rtol=1e-12)
    @test isapprox(K.corner_angle(sq, 1, 1), pi / 2; rtol=1e-12)
    @test isapprox(K.face_angle_at(sq, 1, 3), pi / 2; rtol=1e-12)
    @test_throws ErrorException K.face_angle_at(sq, 1, 7)
    @test K.n_interior_vertices(sq) == 0
end
