# kiri_gen <kind> [p0 p1 ...] --out g.json [--seed s] [--orient auto|checker|brute|none]
#
# Run as
#     julia --project=Kirigami Kirigami/apps/kiri_gen.jl <kind> [params...] --out g.json
include(joinpath(@__DIR__, "common_app.jl"))

function main(args::Vector{String})
    if length(args) < 1
        println(stderr, "usage: kiri_gen <kind> [params...] --out <g.json> [--seed s] ",
                "[--orient auto|checker|brute|none]\n",
                "kinds: triangles squares squares_rect hexagons kagome t3_4_3_12 snub_square\n",
                "       truncated_square periodic_squares periodic_triangles periodic_hexagons\n",
                "       periodic_kagome delaunay voronoi quad_random")
        return 1
    end
    kind = args[1]
    out = "g.json"
    orient = "none"
    seed = UInt32(1)
    params = Float64[]
    i = 2
    while i <= length(args)
        a = args[i]
        if a == "--out" && i + 1 <= length(args)
            out = args[i += 1]
        elseif a == "--seed" && i + 1 <= length(args)
            seed = arg_u(args[i += 1])
        elseif a == "--orient" && i + 1 <= length(args)
            orient = args[i += 1]
        else
            push!(params, arg_f(a))
        end
        i += 1
    end
    rng = K.MT19937(seed)
    m = K.generate(kind, params, rng)
    if orient == "auto"
        m.sigma = K.assign_orientation_relaxation(m, rng).sigma
    elseif orient == "brute"
        m.sigma = K.brute_force_orientation(m).sigma
    elseif orient == "checker"
        m.sigma = checkerboard(m)
    end
    d = dirname(out)
    isempty(d) || mkpath(d)
    K.save_mesh_json(m, out)
    println(kind, ": N=", K.n_vertices(m), " F=", K.n_faces(m), " E=", K.n_edges(m),
            " interior=", K.n_interior_vertices(m), " -> ", out)
    return 0
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
