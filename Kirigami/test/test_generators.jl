# test_generators.jl -- the generators regenerate the
# frozen corpora. For the first 10 ids of every `make_graph` population, the reference
# cases and the scaling cells, the Julia mesh must equal the frozen one: vertex coordinates
# with `==` (bit-exact) for the random generators, which involve no libm transcendental,
# and `isapprox(atol = 1e-12)` for the tilings (`==` as well when the system libm is the
# arm64 Apple libm the corpora were produced with, see generators.jl); face lists and the
# checkerboard sigma exactly. On the first divergence the message names the id, the
# vertex and both values.
include("helpers.jl")
import JSON

const K = Kirigami

# Reads the first `n` rows of a top-level JSON array without parsing the whole (up to
# 165 MB) file: a byte scanner that is string- and escape-aware.
function corpus_first_rows(path::AbstractString, n::Int)
    s = read(path)
    depth = 0
    instr = false
    esc = false
    rows = 0
    stop = 0
    for (i, b) in enumerate(s)
        if instr
            if esc
                esc = false
            elseif b == UInt8('\\')
                esc = true
            elseif b == UInt8('"')
                instr = false
            end
            continue
        end
        if b == UInt8('"')
            instr = true
        elseif b == UInt8('[') || b == UInt8('{')
            depth += 1
        elseif b == UInt8(']') || b == UInt8('}')
            depth -= 1
            if depth == 1
                rows += 1
                if rows == n
                    stop = i
                    break
                end
            end
        end
    end
    stop == 0 && return JSON.parse(String(s))
    return JSON.parse(String(s[1:stop]) * "]")
end

corpus_vertices(row) = [K.Vec2(Float64(v[1]), Float64(v[2])) for v in row["mesh"]["vertices"]]
# file is 0-based -> 1-based
corpus_faces(row) = [Int[Int(i) + 1 for i in f] for f in row["mesh"]["faces"]]

# Reports the first differing vertex, or `nothing` if the two coordinate lists are equal.
function first_vertex_divergence(V::Vector{K.Vec2}, X::Vector{K.Vec2}; exact::Bool)
    length(V) != length(X) && return "vertex count $(length(V)) (corpus) vs $(length(X)) (julia)"
    for i in eachindex(V)
        same = exact ? V[i] == X[i] : isapprox(V[i], X[i]; atol = 1e-12)
        same || return "vertex $i (1-based): corpus $(V[i]) vs julia $(X[i]) (diff $(X[i] - V[i]))"
    end
    return nothing
end

function check_population(file, args; n = 10, exact = true)
    rows = corpus_first_rows(joinpath(CORPUS, file), n)
    @testset "$file first $(length(rows)) ids" begin
        for row in rows
            id = Int(row["id"])
            g = K.make_graph(id, args...)
            @test g.ok == row["ok"]
            g.ok || continue
            @test g.kind == row["kind"]
            m = g.mesh
            @test K.n_vertices(m) == row["N"]
            @test K.n_faces(m) == row["F"]
            @test K.n_edges(m) == row["E"]
            @test K.n_interior_vertices(m) == row["n_interior"]
            d = first_vertex_divergence(corpus_vertices(row), m.X; exact = exact)
            @test d === nothing
            d === nothing || @info "divergence" file id d
            @test corpus_faces(row) == m.faces
            @test Int.(row["sigma_mc"]) == m.sigma
            @test Int.(row["sigma_checker"]) == K.checkerboard(m)
        end
    end
end

@testset "make_graph populations are bit-exact" begin
    check_population("k1a_200.json", (100, 800, 1400))
    check_population("native200.json", (100, 800, 1400))
    check_population("regime_100.json", (20, 100, 1400))
    check_population("yield_fresh_100.json", (100, 800, 1400))
    check_population("e1_900.json", (100, 800, 1400))
    check_population("k2b_random_8.json", (100, 160, 1400))
    check_population("derivation_l1_small.json", (18, 46, 220))
    check_population("k2c_500.json", (100, 5000, 1400))
    check_population("k3a_500.json", (100, 5000, 1 << 30))
end

@testset "scaling cells are bit-exact" begin
    rows = JSON.parsefile(joinpath(CORPUS, "scaling_42.json"))
    for row in rows
        row["sites"] > 1000 && continue   # the 2000/5000-site cells take minutes in Bowyer-Watson
        seed = UInt32(1000003) * UInt32(row["seed"]) + UInt32(20260908) + UInt32(7919) * UInt32(row["sites"])
        @test seed == UInt32(row["rng_seed"])
        rng = K.MT19937(seed)
        m = K.generate(row["kind"], [Float64(row["sites"]), 40.0], rng)
        m.sigma = K.assign_orientation_relaxation(m, rng, 4, 300, 90).sigma
        @test K.n_vertices(m) == row["N"]
        @test K.n_faces(m) == row["F"]
        d = first_vertex_divergence(corpus_vertices(row), m.X; exact = true)
        @test d === nothing
        d === nothing || @info "divergence" row["kind"] row["sites"] row["seed"] d
        @test corpus_faces(row) == m.faces
        @test Int.(row["sigma_mc"]) == m.sigma
    end
end

@testset "reference cases (tilings) reproduce the frozen meshes" begin
    rows = JSON.parsefile(joinpath(CORPUS, "reference_cases_8.json"))
    cases = K.reference_cases()
    @test length(cases) == length(rows)
    for (rc, row) in zip(cases, rows)
        @test rc.name == row["name"]
        @test rc.periodic == row["periodic"]
        @test K.n_vertices(rc.mesh) == row["N"]
        @test K.n_faces(rc.mesh) == row["F"]
        d = first_vertex_divergence(corpus_vertices(row), rc.mesh.X; exact = K._USE_SYSTEM_LIBM)
        @test d === nothing
        d === nothing || @info "divergence" rc.name d
        @test corpus_faces(row) == rc.mesh.faces
        # sigma: exact for the checkerboard cases; the relaxation cases are only expected
        # to be exact with the arm64 libm (their symmetric patches have rounding ties)
        sigma = Int.(row["sigma"])
        if rc.name in ("rotating_squares", "triangles_alternating", "kagome_3636", "periodic_squares_4x4") ||
           K._USE_SYSTEM_LIBM
            @test sigma == rc.mesh.sigma
        else
            # must at least be an equivalent orientation: same split / component / hole counts
            a = K.describe_orientation(rc.mesh, rc.mesh.sigma)
            b = K.describe_orientation(rc.mesh, sigma)
            @test (a.n_split, a.components, a.n_holes) == (b.n_split, b.components, b.n_holes)
        end
        if rc.periodic
            p = row["mesh"]["periodic"]
            @test rc.mesh.periodic.th == K.Vec2(p["th"][1], p["th"][2])
            @test rc.mesh.periodic.tv == K.Vec2(p["tv"][1], p["tv"][2])
            # file is 0-based -> 1-based
            @test rc.mesh.periodic.pairs_h == [(Int(q[1]) + 1, Int(q[2]) + 1) for q in p["pairs_h"]]
            @test rc.mesh.periodic.pairs_v == [(Int(q[1]) + 1, Int(q[2]) + 1) for q in p["pairs_v"]]
        end
    end
end

@testset "generate: named kinds and errors" begin
    rng = K.MT19937(1)
    @test_throws ErrorException K.generate("no_such_kind", Float64[], rng)
    for kind in ("triangles", "squares", "hexagons", "kagome", "t3_4_3_12", "snub_square", "truncated_square")
        m = K.generate(kind, [2.0], rng)
        @test K.n_faces(m) > 0
        @test all(f -> K.face_signed_area(m, f) > 0, 1:K.n_faces(m))
    end
    m = K.generate("squares_rect", [3.0, 2.0], rng)
    @test K.n_faces(m) == 6
    @test K.n_vertices(m) == 12
    for kind in ("periodic_squares", "periodic_triangles", "periodic_hexagons", "periodic_kagome")
        m = K.generate(kind, [3.0, 2.0], rng)
        @test m.periodic.present
        @test !isempty(m.periodic.pairs_h) && !isempty(m.periodic.pairs_v)
    end
    m = K.periodic_squares(4, 4)
    @test K.n_faces(m) == 16 && K.n_vertices(m) == 25
    @test length(m.periodic.pairs_h) == 5 && length(m.periodic.pairs_v) == 5
    t = K.torus_squares(4, 3)
    @test K.n_vertices(t) == 12 && K.n_faces(t) == 12
    @test K.n_interior_vertices(t) == 12   # no boundary vertex on the torus
    t = K.torus_triangles(4, 3)
    @test K.n_faces(t) == 24 && K.n_interior_vertices(t) == 12
    # random kinds: seed-deterministic
    for kind in ("delaunay", "voronoi", "quad_random")
        m1 = K.generate(kind, [30.0, 10.0], K.MT19937(42))
        m2 = K.generate(kind, [30.0, 10.0], K.MT19937(42))
        @test m1.X == m2.X && m1.faces == m2.faces
        @test K.n_faces(m1) > 10
    end
    @test isempty(K.delaunay_triangles(K.Vec2[]))
    @test K.delaunay_triangles([K.Vec2(0, 0), K.Vec2(1, 0), K.Vec2(0, 1)]) == [(1, 2, 3)]
end

@testset "brute_force_orientation agrees with the relaxation on a small patch" begin
    m = K.tiling_squares(K.rect(K.Vec2(1.5, 1.5), 1.51, 1.51))   # 3x3 squares, F = 9
    @test K.n_faces(m) == 9
    bf = K.brute_force_orientation(m)
    @test bf.n_split == 0 && bf.components == 1
    rl = K.assign_orientation_relaxation(m, K.MT19937(3), 4, 300, 90)
    @test rl.n_split == 0 && rl.components == 1
    @test abs.(rl.sigma) == ones(Int, 9)
    @test K.checkerboard(m) == bf.sigma || K.checkerboard(m) == -bf.sigma
    @test_throws ErrorException K.brute_force_orientation(K.tiling_squares(K.rect(K.Vec2(2.5, 2.5), 2.51, 2.51)))
end
