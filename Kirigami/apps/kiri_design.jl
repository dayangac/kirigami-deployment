# kiri_design -- the deliverable method as a command-line tool.
#
#   kiri_design in.json --sigma mc|def|both|json --out design.json [--proximity] [--baseline]
#
# For each requested orientation the tool runs, BY
# DEFAULT, design_range_max (K9c: stage A maximises the exact 0+ margin from each of the
# starts {K9 point, K9b point, t = 0}, stage B pushes the exact Theta_max, and the answer
# is the best of {K9, K9b, stage-A/B winner}; the winning arm is printed and stored as
# `provenance`), and on request design_constrained (--proximity, the K9 variant-(b)
# embedding) and design_baseline (--baseline, the 2026 Eq. (6) projection alone). All
# three are characterized with the same exact T4.2'' range and the same repaired validity
# certificate, and written to one JSON.
#
# DETERMINISM. Every solve here is a pure function of (mesh, sigma, X_ini, options,
# --seed): no clock, no thread pool and no global state enters it, so the same command
# line gives bit-identical output on the same build.
include(joinpath(@__DIR__, "exp_common.jl"))

function usage()
    println(stderr,
        "usage: kiri_design <in.json> [--sigma mc|def|both|json] --out <design.json>\n",
        "                   [--proximity] [--baseline] [--eps r] [--delta d]\n",
        "                   [--split-delta d] [--delta-wide d] [--stages n]\n",
        "                   [--stage-iters n] [--no-stage-b] [--no-arm-k9]\n",
        "                   [--no-arm-k9b] [--restarts n] [--barrier n] [--iters n]\n",
        "                   [--seed s] [--maximise-eps] [--referee]\n",
        "                   [--no-warm-starts] [--defect-cap m] [--out-graph g.json]\n",
        "  (default)     K9c range-maximising design (design_range_max)\n",
        "  --proximity   the K9 proximity-constrained design (design_constrained)\n",
        "  --baseline    additionally report the Eq. (6) projection alone\n",
        "  --sigma json  uses the `orientation` field of the input (an error if absent)\n",
        "  --sigma mc    Eq. (1) max-cut orientation (2026 Sec. 4.2)\n",
        "  --sigma def   K5's defect-minimizing orientation, started from mc\n",
        "  --sigma both  runs mc and def and reports both\n",
        "  --stages/--stage-iters are stage A's continuation stages and L-BFGS\n",
        "  iterations per stage; the defaults 6 / 120 are the K9c run's settings, so\n",
        "  --seed \$((9300 + 7 * id + which)) reproduces results/experiments/k9c/k9c.csv.")
end

function certificate_json(v::K.ValidityCertificate)
    return Dict{String,Any}("valid" => K.valid(v), "pos" => v.pos, "nooverlap" => v.nooverlap,
                            "noroot" => v.noroot, "eps" => v.eps, "theta_1" => v.theta_1,
                            "min_signed_area" => v.min_signed_area, "n_inverted" => v.n_inverted,
                            "n_pairs" => v.n_pairs, "n_candidates" => v.n_candidates,
                            "n_identically_zero" => v.n_identically_zero,
                            "n_roots_deflated" => v.n_roots_deflated,
                            "n_roots_inadmissible" => v.n_roots_inadmissible,
                            "first_root" => v.first_root)
end

# M'-vertex / face ids of the witness are 0-based in the file.
_idx0(i::Int) = i < 1 ? -1 : i - 1

function characterization_json(ch::K.Characterization)
    j = Dict{String,Any}()
    j["theta_max"] = ch.theta_max
    j["zero_range"] = ch.zero_range
    j["binding"] = ch.binding
    j["n_contacts"] = length(ch.contacts)
    j["contacts"] = copy(ch.contacts)
    fc = ch.first_contact
    j["first_contact"] = Dict{String,Any}("found" => fc.found, "theta" => fc.theta,
                                          "pv" => _idx0(fc.pv), "pa" => _idx0(fc.pa),
                                          "pb" => _idx0(fc.pb), "face_v" => _idx0(fc.face_v),
                                          "face_e" => _idx0(fc.face_e))
    j["certified"] = ch.certified
    j["eps_max"] = ch.eps_max
    j["certificate"] = certificate_json(ch.certificate)
    j["min_signed_area"] = ch.min_signed_area
    j["n_inverted"] = ch.n_inverted
    j["n_split"] = ch.n_split
    j["n_hinge"] = ch.n_hinge
    j["med_edge"] = ch.med_edge
    ch.theta_bisect >= 0 && (j["theta_bisect"] = ch.theta_bisect)
    return j
end

function design_json(m::K.Mesh, d::K.DesignResult, rm::Union{K.RangeMaxResult,Nothing} = nothing)
    j = Dict{String,Any}()
    j["method"] = d.method
    j["ok"] = d.ok
    j["status"] = d.status
    # Every design entry carries the README contract fields (vertices, faces,
    # orientation) so it can be lifted out as a graph on its own; --out-graph writes the
    # best one as a plain graph JSON that kiri_analyze / kiri_export / kiri_deploy read.
    j["vertices"] = [[x[1], x[2]] for x in d.X]
    j["faces"] = [[i - 1 for i in f] for f in m.faces]   # 1-based -> 0-based in the file
    j["orientation"] = copy(m.sigma)
    j["X0"] = [[x[1], x[2]] for x in d.X0]
    j["dim_null"] = d.dim_null
    j["feasible"] = d.feasible
    j["min_cross_rel"] = d.min_cross
    j["min_q_rel"] = d.min_q
    j["min_mu_rel"] = d.min_mu
    j["n_bad"] = d.n_bad
    j["n_bad_q"] = d.n_bad_q
    j["n_bad_mu"] = d.n_bad_mu
    j["n_corners"] = d.n_corners
    j["n_nonconvex_x0"] = d.n_nonconvex_x0
    j["n_inverted"] = d.n_inverted
    j["dist_ini"] = d.dist_ini
    j["dist_ini_x0"] = d.dist_ini_x0
    j["t_norm_rel"] = d.t_norm_rel
    j["barrier_stages_kept"] = d.barrier_stages_kept
    j["used_range_opt"] = d.used_range_opt
    j["med_edge"] = d.med_edge
    j["characterization"] = characterization_json(d.ch)
    # The K9c provenance: which arm supplied the point that was returned, whether stage B
    # moved it, and every candidate that was scored on the way.
    if rm !== nothing
        j["provenance"] = rm.provenance
        j["stage_b_used"] = rm.stage_b_used
        j["margin_zero_plus"] = rm.margin
        j["arms"] = [Dict{String,Any}("tag" => a.tag, "feasible" => a.feasible,
                                      "margin" => a.margin, "theta_max" => a.theta_max,
                                      "eps_max" => a.eps_max,
                                      "stage_a_winner" => a.stage_a_winner) for a in rm.arms]
    end
    return j
end

function summary_line(tag::AbstractString, d::K.DesignResult,
                      rm::Union{K.RangeMaxResult,Nothing} = nothing)
    print(tag, ": Theta_max = ", fx(d.ch.theta_max, 6), " rad, eps_max = ", fx(d.ch.eps_max, 6),
          " rad, certified = ", (d.ch.certified ? "yes" : "no"),
          ", feasible = ", (d.feasible ? "yes" : "no"), " (min cross ", sci(d.min_cross),
          ", min q ", sci(d.min_q), " med^2)", ", inverted faces ", d.n_inverted)
    d.ch.theta_max <= 1e-9 && print(", binding ", d.ch.binding)
    if rm !== nothing
        print(", winning arm = ", rm.provenance,
              (rm.stage_b_used ? " (stage B improved it)" : " (stage A)"),
              ", 0+ margin = ", sci(rm.margin), " med^2")
    end
    d.ok || print("  [", d.status, "]")
    println()
    if rm !== nothing
        for a in rm.arms
            println("    arm ", a.tag, ": Theta_max = ", fx(a.theta_max, 6), ", eps_max = ",
                    fx(a.eps_max, 6), ", m = ", sci(a.margin), ", feasible = ",
                    (a.feasible ? "yes" : "no"), (a.stage_a_winner ? "   <- stage A winner" : ""))
        end
    end
end

function main(args::Vector{String})
    in_ = out = out_graph = ""
    sigma_mode = "mc"
    baseline = false
    proximity = false
    defect_cap_mult = 20
    opt = K.DesignOptions()
    rmo = K.RangeMaxOptions()
    i = 1
    n = length(args)
    while i <= n
        a = args[i]
        if a == "--sigma" && i + 1 <= n
            sigma_mode = args[i += 1]
        elseif a == "--out" && i + 1 <= n
            out = args[i += 1]
        elseif a == "--out-graph" && i + 1 <= n
            out_graph = args[i += 1]
        elseif a == "--baseline"
            baseline = true
        elseif a == "--proximity"
            proximity = true
        elseif a == "--delta-wide" && i + 1 <= n
            rmo.delta_wide = arg_f(args[i += 1])
        elseif a == "--stages" && i + 1 <= n
            rmo.stages = arg_i(args[i += 1])
        elseif a == "--stage-iters" && i + 1 <= n
            rmo.iter_per_stage = arg_i(args[i += 1])
        elseif a == "--no-stage-b"
            rmo.stage_b = false
        elseif a == "--no-arm-k9"
            rmo.arm_proximity = false
        elseif a == "--no-arm-k9b"
            rmo.arm_proximity_wide = false
        elseif a == "--eps" && i + 1 <= n
            opt.characterize.eps = arg_f(args[i += 1])
        elseif a == "--delta" && i + 1 <= n
            opt.delta_convex = arg_f(args[i += 1])
        elseif a == "--split-delta" && i + 1 <= n
            opt.delta_split = arg_f(args[i += 1])
        elseif a == "--restarts" && i + 1 <= n
            opt.restarts = arg_i(args[i += 1])
        elseif a == "--barrier" && i + 1 <= n
            opt.barrier_stages = arg_i(args[i += 1])
        elseif a == "--iters" && i + 1 <= n
            opt.max_iter = arg_i(args[i += 1])
        elseif a == "--seed" && i + 1 <= n
            opt.seed = arg_u(args[i += 1])
        elseif a == "--maximise-eps"
            opt.maximise_eps = true
        elseif a == "--referee"
            opt.characterize.referee = true
        elseif a == "--no-warm-starts"
            opt.warm_starts = false
        elseif a == "--defect-cap" && i + 1 <= n
            defect_cap_mult = arg_i(args[i += 1])
        elseif a == "-h" || a == "--help"
            usage()
            return 0
        elseif !isempty(a) && a[1] == '-'
            usage()
            return 2
        elseif isempty(in_)
            in_ = a
        else
            usage()
            return 2
        end
        i += 1
    end
    if isempty(in_) || isempty(out)
        usage()
        return 2
    end
    # The knobs the two methods share are parsed once, into DesignOptions, and mirrored
    # here, so that --eps / --delta / --split-delta / --iters / --seed / --referee mean
    # the same thing whichever method runs.
    rmo.delta_convex = opt.delta_convex
    rmo.delta_split = opt.delta_split
    rmo.max_iter = opt.max_iter
    rmo.seed = opt.seed
    rmo.characterize = opt.characterize
    if opt.maximise_eps && !proximity
        println(stderr, "kiri_design: --maximise-eps applies to --proximity only; the default ",
                "method already maximises the range")
    end
    if !(sigma_mode in ("mc", "def", "both", "json"))
        usage()
        return 2
    end

    m = try
        K.load_mesh_json(in_)
    catch e
        println(stderr, "kiri_design: cannot read ", in_, ": ", sprint(showerror, e))
        return 1
    end
    K.build_topology!(m)
    if K.n_faces(m) == 0
        println(stderr, "kiri_design: ", in_, " has no faces")
        return 1
    end

    # --- the orientations to run ----------------------------------------------
    sigmas = Tuple{String,Vector{Int}}[]
    if sigma_mode == "json"
        if length(m.sigma) != K.n_faces(m)
            println(stderr, "kiri_design: --sigma json but the input has no `orientation` field")
            return 1
        end
        push!(sigmas, ("json", m.sigma))
    else
        mc = K.orientation_maxcut(m, opt.seed)
        if length(mc) != K.n_faces(m)
            println(stderr, "kiri_design: Eq. (1) orientation failed")
            return 1
        end
        (sigma_mode == "mc" || sigma_mode == "both") && push!(sigmas, ("mc", mc))
        if sigma_mode == "def" || sigma_mode == "both"
            dr = K.orientation_defect(m, mc, defect_cap_mult * K.n_faces(m), opt.seed)
            if !dr.ok
                println(stderr, "kiri_design: the defect search found no sigma with c(Gamma) = 1")
                sigma_mode == "def" && return 1
            else
                push!(sigmas, ("def", dr.sigma))
            end
        end
    end

    j = Dict{String,Any}()
    j["input"] = in_
    j["n_vertices"] = K.n_vertices(m)
    j["n_faces"] = K.n_faces(m)
    j["options"] = Dict{String,Any}("sigma" => sigma_mode, "eps" => opt.characterize.eps,
                                    "delta_convex" => opt.delta_convex,
                                    "delta_split" => opt.delta_split,
                                    "restarts" => opt.restarts,
                                    "barrier_stages" => opt.barrier_stages,
                                    "max_iter" => opt.max_iter, "seed" => Int(opt.seed),
                                    "warm_starts" => opt.warm_starts,
                                    "maximise_eps" => opt.maximise_eps,
                                    "referee" => opt.characterize.referee,
                                    "method" => (proximity ? "constrained" : "range_max"),
                                    "delta_wide" => rmo.delta_wide, "stages" => rmo.stages,
                                    "stage_iters" => rmo.iter_per_stage,
                                    "stage_b" => rmo.stage_b, "arm_k9" => rmo.arm_proximity,
                                    "arm_k9b" => rmo.arm_proximity_wide)
    j["designs"] = Any[]

    rc = 0
    # --out-graph gets the constrained design with the largest exact Theta_max, ties
    # broken by the certified eps_max -- the same ranking K9b uses.
    have_best = false
    best_mesh = K.Mesh()
    best = K.DesignResult()
    for (src, sig) in sigmas
        ms = K.Mesh(m.X, m.faces)
        ms.sigma = sig
        ms.periodic = m.periodic
        K.build_topology!(ms)
        try
            rm = nothing
            d = K.DesignResult()
            if proximity
                d = K.design_constrained(ms, sig, ms.X, opt)
            else
                rm = K.design_range_max(ms, sig, ms.X, rmo)
                d = rm.design
            end
            if !have_best || d.ch.theta_max > best.ch.theta_max + 1e-12 ||
               (abs(d.ch.theta_max - best.ch.theta_max) <= 1e-12 && d.ch.eps_max > best.ch.eps_max)
                have_best = true
                best = d
                best_mesh = ms
            end
            dj = design_json(ms, d, proximity ? nothing : rm)
            dj["sigma_source"] = src
            push!(j["designs"], dj)
            summary_line("sigma_" * src * (proximity ? " constrained" : " range_max  "), d,
                         proximity ? nothing : rm)
            if baseline
                b = K.design_baseline(ms, sig, ms.X, opt)
                bj = design_json(ms, b)
                bj["sigma_source"] = src
                push!(j["designs"], bj)
                summary_line("sigma_" * src * " baseline   ", b)
            end
        catch e
            println(stderr, "kiri_design: sigma_", src, ": ", sprint(showerror, e))
            rc = 1
        end
    end

    write_json(j, out)
    println("wrote ", out)
    if !isempty(out_graph)
        if !have_best
            println(stderr, "kiri_design: nothing to write to ", out_graph)
            return 1
        end
        best_mesh.X = best.X
        write_json(mesh_json(best_mesh), out_graph)
        println("wrote ", out_graph, " (the best constrained design, Theta_max = ",
                fx(best.ch.theta_max, 6), " rad)")
    end
    return rc
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
