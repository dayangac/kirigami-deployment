# core/kill_common.jl -- the graph populations shared by the Phase-5 kill experiments and
# the eight Phase-2 reference cases (the parts that define INPUTS; the shape-space cache
# helpers live app-side in apps/kill_common.jl).
#
# Everything is a deterministic function of the graph id through the bit-exact `MT19937`,
# so `make_graph(id, ...)` regenerates the frozen populations in data/corpus.

mutable struct Graph
    mesh::Mesh
    kind::String
    id::Int
    ok::Bool
end
Graph() = Graph(Mesh(), "", -1, false)

"""2-colouring of the dual graph by DFS from face 1 (K2c's alternative sigma)."""
function checkerboard(m::Mesh)
    adj = dual_graph(m)
    sig = zeros(Int, n_faces(m))
    for s in 1:n_faces(m)
        sig[s] != 0 && continue
        sig[s] = -1
        st = [s]
        while !isempty(st)
            f = pop!(st)
            for g in adj[f]
                if sig[g] == 0
                    sig[g] = -sig[f]
                    push!(st, g)
                end
            end
        end
    end
    return sig
end

const _KILL_KINDS = ("voronoi", "delaunay", "quad_random")

"""
    make_graph(id, min_faces, max_faces, n_cap) -> Graph

One graph of the population: kind = voronoi / delaunay / quad_random by id % 3,
`rng = MT19937(1000003*id + 20260903)`, site count drawn so that the face count lands in
[min_faces, max_faces] and the vertex count stays <= n_cap (up to 8 attempts on the same
rng); sigma from `assign_orientation_relaxation(m, rng, 4, 300, 90)`.
"""
function make_graph(id::Int, min_faces::Int, max_faces::Int, n_cap::Int)
    g = Graph()
    g.id = id
    g.kind = _KILL_KINDS[id % 3 + 1]
    rng = MT19937((UInt32(1000003) * UInt32(id) + UInt32(20260903)) % UInt32)
    # geometric-ish spread over the face range
    target = min_faces * libm_pow(Float64(max_faces) / min_faces, uniform_real(rng, 0.0, 1.0))
    for _ in 1:8
        sites = if g.kind == "delaunay"
            max(20, trunc(Int, target / 2.0))
        elseif g.kind == "voronoi"
            max(20, trunc(Int, target * 1.15))
        else
            max(20, trunc(Int, target * 1.6))
        end
        m = try
            generate(g.kind, [Float64(sites), 40.0], rng)
        catch
            target *= 1.2
            continue
        end
        if n_faces(m) < min_faces
            target *= max(1.1, Float64(min_faces) / max(1, n_faces(m)))
            continue
        end
        if n_faces(m) > max_faces || n_vertices(m) > n_cap
            target *= 0.8
            continue
        end
        m.sigma = assign_orientation_relaxation(m, rng, 4, 300, 90).sigma
        isempty(m.sigma) && continue
        g.mesh = m
        g.ok = true
        return g
    end
    return g
end

struct RefCase
    name::String
    mesh::Mesh
    periodic::Bool
end

"""The eight Phase-2 reference cases, rebuilt exactly as apps/kiri_reference.jl does: ONE shared
`MT19937(20260903)` feeds the four relaxation calls in file order."""
function reference_cases()
    rng = MT19937(20260903)
    cs = RefCase[]
    let g = tiling_squares(rect(Vec2(2.5, 2.5), 2.51, 2.51))
        g.sigma = checkerboard(g)
        push!(cs, RefCase("rotating_squares", g, false))
    end
    let g = tiling_triangles(disk(Vec2(0.13, 0.07), 3.5))
        g.sigma = checkerboard(g)
        push!(cs, RefCase("triangles_alternating", g, false))
    end
    let g = tiling_kagome(disk(Vec2(0.13, 0.07), 3.0))
        g.sigma = checkerboard(g)
        push!(cs, RefCase("kagome_3636", g, false))
    end
    let g = tiling_hexagons(disk(Vec2(0.13, 0.07), 3.0))
        g.sigma = assign_orientation_relaxation(g, rng, 8, 500, 180).sigma
        push!(cs, RefCase("hexagons_auto", g, false))
    end
    let g = tiling_truncated_square(disk(Vec2(0.13, 0.07), 4.0))
        g.sigma = assign_orientation_relaxation(g, rng, 8, 500, 180).sigma
        push!(cs, RefCase("truncated_square_488", g, false))
    end
    let g = tiling_snub_square(disk(Vec2(0.13, 0.07), 3.2))
        g.sigma = assign_orientation_relaxation(g, rng, 8, 500, 180).sigma
        push!(cs, RefCase("snub_square_33434", g, false))
    end
    let g = tiling_3_4_3_12(disk(Vec2(0.13, 0.07), 4.2))
        g.sigma = assign_orientation_relaxation(g, rng, 8, 500, 180).sigma
        push!(cs, RefCase("tiling_3_4_3_12", g, false))
    end
    let g = periodic_squares(4, 4)
        g.sigma = checkerboard(g)
        push!(cs, RefCase("periodic_squares_4x4", g, true))
    end
    return cs
end
