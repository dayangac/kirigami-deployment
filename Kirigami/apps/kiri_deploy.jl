# kiri_deploy <graph.json> --theta t --out y.json
#
include(joinpath(@__DIR__, "common_app.jl"))

function main(args::Vector{String})
    if length(args) < 1
        println(stderr, "usage: kiri_deploy <graph.json> --theta <t> --out <y.json>")
        return 1
    end
    in_ = args[1]
    out = "y.json"
    theta = 0.5
    i = 2
    while i <= length(args)
        a = args[i]
        if a == "--theta" && i + 1 <= length(args)
            theta = arg_f(args[i += 1])
        elseif a == "--out" && i + 1 <= length(args)
            out = args[i += 1]
        end
        i += 1
    end
    m = K.load_mesh_json(in_)
    if isempty(m.sigma)
        println(stderr, "error: input has no \"orientation\" array; run kiri_gen --orient or ",
                "kiri_analyze --orient auto first")
        return 2
    end
    c = K.make_cut(m)
    d = K.deploy(c, m.X, theta)
    write_json(deployment_json(c, m.X, theta), out)
    println("theta=", fmt_g(theta), " max_mismatch=", fmt_g(d.max_mismatch), " -> ", out)
    return 0
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
