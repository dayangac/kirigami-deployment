# helpers.jl -- shared test fixtures.
using Kirigami, Test, LinearAlgebra
const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const CORPUS = joinpath(REPO, "data", "corpus")

import JSON

# A mesh JSON object (0-based file format) -> Mesh with faces taken AS STORED: no CCW
# normalisation (the torus fixtures have geometrically degenerate wrap faces that
# `mesh_from_json` would flip), only `build_topology!`. `orientation` and `periodic`
# are read like `mesh_from_json`. Used by the test_fixtures_*.json loaders.
# (helpers.jl is included once per test file; the guard avoids redefinition warnings.)
if !@isdefined(fixture_mesh_raw)
function fixture_mesh_raw(j)
    X = [Kirigami.Vec2(Float64(v[1]), Float64(v[2])) for v in j["vertices"]]
    faces = [[Int(i) + 1 for i in f] for f in j["faces"]]  # 0-based file -> 1-based
    m = Kirigami.Mesh(X, faces)
    if haskey(j, "orientation") && j["orientation"] !== nothing
        m.sigma = [Int(s) for s in j["orientation"]]
    end
    if haskey(j, "periodic") && j["periodic"] !== nothing
        p = j["periodic"]
        m.periodic.present = true
        m.periodic.th = Kirigami.Vec2(Float64(p["th"][1]), Float64(p["th"][2]))
        m.periodic.tv = Kirigami.Vec2(Float64(p["tv"][1]), Float64(p["tv"][2]))
        m.periodic.pairs_h = [(Int(q[1]) + 1, Int(q[2]) + 1) for q in get(p, "pairs_h", [])]
        m.periodic.pairs_v = [(Int(q[1]) + 1, Int(q[2]) + 1) for q in get(p, "pairs_v", [])]
    end
    Kirigami.build_topology!(m)
    return m
end

# [[x, y], ...] -> Vector{Vec2}
fixture_points(a) = [Kirigami.Vec2(Float64(p[1]), Float64(p[2])) for p in a]
# [[...], ...] -> Matrix{Float64}
fixture_matrix(a) = isempty(a) ? zeros(0, 0) : Float64[a[i][j] for i in eachindex(a), j in eachindex(a[1])]

# true iff `m` has bit-identical vertices and identical face lists to the fixture mesh `j`
function same_mesh_as_fixture(m::Kirigami.Mesh, j)
    length(m.X) == length(j["vertices"]) || return false
    all(m.X[i] == Kirigami.Vec2(Float64(j["vertices"][i][1]), Float64(j["vertices"][i][2]))
        for i in eachindex(m.X)) || return false
    return m.faces == [[Int(i) + 1 for i in f] for f in j["faces"]]
end

# --- sigma helpers ---------------------------------------------------------------------

# Checkerboard sigma on a mesh whose dual is bipartite (else falls back to DFS
# 2-colouring, leaving some monochromatic dual edges = split cuts).
function checkerboard_sigma(m::Kirigami.Mesh)
    adj = Kirigami.dual_graph(m)
    sig = zeros(Int, Kirigami.n_faces(m))
    for s in 1:Kirigami.n_faces(m)
        sig[s] != 0 && continue
        sig[s] = -1
        stack = Int[s]
        while !isempty(stack)
            f = pop!(stack)
            for g in adj[f]
                if sig[g] == 0
                    sig[g] = -sig[f]
                    push!(stack, g)
                end
            end
        end
    end
    return sig
end

# bernoulli(0.5) per face, true -> +1 (libc++ semantics: uniform_real(0,1) < p)
random_sigma(m::Kirigami.Mesh, rng::Kirigami.MT19937) =
    [Kirigami.uniform_real(rng, 0.0, 1.0) < 0.5 ? 1 : -1 for _ in 1:Kirigami.n_faces(m)]
end
