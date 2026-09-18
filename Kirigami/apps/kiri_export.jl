# kiri_export <graph.json> --profile felt_laser --theta 0|<rad> --svg out.svg
#             --stl out.stl --3mf out.3mf [--embedding X.json] [--scale mm/unit]
#             [--theta-frac f] [--json report.json]
#
# Writes fabricable geometry for an oriented planar graph: a laser SVG in mm with
# cut / score / engrave layers, and a 3D solid as binary STL and 3MF.
include(joinpath(@__DIR__, "common_app.jl"))

function usage()
    println(stderr, "usage: kiri_export <graph.json> [--profile felt_laser|paper_laser|pla_print]\n",
            "                   [--theta <rad>] [--theta-frac <f of theta_max>]\n",
            "                   [--svg out.svg] [--stl out.stl] [--3mf out.3mf]\n",
            "                   [--embedding X.json] [--scale <mm per input unit>]\n",
            "                   [--json report.json]")
end

function ensure_dir(path::AbstractString)
    d = dirname(path)
    isempty(d) || mkpath(d)
end

function main(args::Vector{String})
    if length(args) < 1
        usage()
        return 1
    end
    in_ = args[1]
    profile_name = "felt_laser"
    svg_path = stl_path = mf_path = emb_path = json_path = ""
    theta = 0.0
    theta_frac = -1.0
    scale = 10.0
    i = 2
    while i <= length(args)
        a = args[i]
        function next()
            i + 1 > length(args) && error("missing value after " * a)
            return args[i += 1]
        end
        try
            if a == "--profile"
                profile_name = next()
            elseif a == "--theta"
                theta = arg_f(next())
            elseif a == "--theta-frac"
                theta_frac = arg_f(next())
            elseif a == "--svg"
                svg_path = next()
            elseif a == "--stl"
                stl_path = next()
            elseif a == "--3mf"
                mf_path = next()
            elseif a == "--embedding"
                emb_path = next()
            elseif a == "--scale"
                scale = arg_f(next())
            elseif a == "--json"
                json_path = next()
            else
                println(stderr, "unknown argument: ", a)
                usage()
                return 1
            end
        catch e
            println(stderr, "error: ", sprint(showerror, e))
            return 1
        end
        i += 1
    end

    try
        m = K.load_mesh_json(in_)
        if isempty(m.sigma)
            println(stderr, "error: ", in_, " has no \"orientation\" array; run kiri_gen --orient ",
                    "or kiri_analyze --orient auto first")
            return 2
        end
        X = m.X
        if !isempty(emb_path)
            e = K.load_mesh_json(emb_path)
            if length(e.X) != length(m.X)
                println(stderr, "error: embedding ", emb_path, " has ", length(e.X),
                        " vertices, graph has ", length(m.X))
                return 2
            end
            X = e.X
        end
        profile = K.profile_by_name(profile_name)
        c = K.make_cut(m)

        theta_max_used = 0.0
        if theta_frac >= 0
            tm = K.theta_max(c, X)
            theta_max_used = tm.theta_max_geometric
            theta = theta_frac * theta_max_used
        end

        L = K.build_layout(m, c, X, theta, scale, profile)
        for w in L.warnings
            println(stderr, "warning: ", w)
        end

        print("profile=", profile.name, " hinge=", K.to_string(profile.hinge), " theta=", fmt_g(theta))
        theta_frac >= 0 && print(" (", fmt_g(theta_frac), " x theta_max=", fmt_g(theta_max_used), ")")
        println(" faces=", K.n_faces(L), " hinges=", length(L.hinges), " bbox=", fmt_g(K.width(L)),
                "x", fmt_g(K.height(L)), " mm")

        rep = Dict{String,Any}()
        rep["input"] = in_
        rep["profile"] = profile.name
        rep["hinge"] = K.to_string(profile.hinge)
        rep["theta"] = theta
        rep["scale_mm_per_unit"] = scale
        rep["faces"] = K.n_faces(L)
        rep["hinge_sites"] = length(L.hinges)
        rep["split_edges"] = K.n_split(c)
        rep["bbox_mm"] = [K.width(L), K.height(L)]
        rep["warnings"] = copy(L.warnings)

        if !isempty(svg_path)
            ensure_dir(svg_path)
            K.write_svg(L, svg_path)
            println("  svg -> ", svg_path)
            rep["svg"] = svg_path
        end
        if !isempty(stl_path) || !isempty(mf_path) || !isempty(json_path)
            solid, w3 = K.build_solid(L)
            for w in w3
                println(stderr, "warning: ", w)
            end
            mr = K.check_manifold(solid, 1e-6)
            println("  solid: ", K.n_tris(solid), " triangles, ", K.summary(mr))
            rep["solid"] = Dict{String,Any}("triangles" => K.n_tris(solid),
                                            "closed" => mr.closed,
                                            "consistently_oriented" => mr.consistently_oriented,
                                            "components" => mr.n_components,
                                            "volume_mm3" => mr.volume,
                                            "boundary_edges" => mr.n_boundary_edges,
                                            "nonmanifold_edges" => mr.n_nonmanifold_edges,
                                            "warnings" => w3)
            if !isempty(stl_path)
                ensure_dir(stl_path)
                K.write_stl_binary(solid, stl_path, "kiri " * profile.name)
                println("  stl -> ", stl_path)
                rep["stl"] = stl_path
            end
            if !isempty(mf_path)
                ensure_dir(mf_path)
                K.write_3mf(solid, mf_path, "kiri " * profile.name)
                n_obj = length(K.split_components(solid))
                println("  3mf -> ", mf_path, " (", n_obj, " objects)")
                rep["3mf"] = mf_path
                rep["solid"]["objects"] = n_obj
            end
        end
        if !isempty(json_path)
            ensure_dir(json_path)
            write_json(rep, json_path)
            println("  report -> ", json_path)
        end
    catch e
        println(stderr, "error: ", sprint(showerror, e))
        return 3
    end
    return 0
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
