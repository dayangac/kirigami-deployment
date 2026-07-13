# kiri_analyze <graph.json> [--orient auto|json|brute] [--boundary fixed|periodic|none]
#              [--out <dir>] [--collision]
#
# Port of code/apps/kiri_analyze.cpp. Edge / vertex ids in the written JSON are 0-based,
# as the C++ wrote them.
include(joinpath(@__DIR__, "common_app.jl"))

zero_based(v::Vector{Int}) = [i - 1 for i in v]   # internal 1-based -> 0-based in the file

function main(args::Vector{String})
    if length(args) < 1
        println(stderr, "usage: kiri_analyze <graph.json> [--orient auto|json|brute] ",
                "[--boundary fixed|periodic|none] [--out <dir>] [--collision] [--seed s]")
        return 1
    end
    in_ = args[1]
    orient = "json"
    bnd = "fixed"
    outdir = "out"
    do_collision = false
    seed = UInt32(1)
    i = 2
    while i <= length(args)
        a = args[i]
        if a == "--orient" && i + 1 <= length(args)
            orient = args[i += 1]
        elseif a == "--boundary" && i + 1 <= length(args)
            bnd = args[i += 1]
        elseif a == "--out" && i + 1 <= length(args)
            outdir = args[i += 1]
        elseif a == "--collision"
            do_collision = true
        elseif a == "--seed" && i + 1 <= length(args)
            seed = arg_u(args[i += 1])
        end
        i += 1
    end
    rng = K.MT19937(seed)
    m = K.load_mesh_json(in_)
    if orient == "auto"
        m.sigma = K.assign_orientation_relaxation(m, rng).sigma
    elseif orient == "brute"
        m.sigma = K.brute_force_orientation(m).sigma
    end
    if isempty(m.sigma)
        println(stderr, "error: no orientation in the JSON and --orient json requested")
        return 2
    end
    mkpath(outdir)

    c = K.make_cut(m)
    a1, bad = K.check_remark_A1(c)

    hs_seed = K.holes_seed_growing(c)
    hs_part = K.holes_partition(c)
    holes_agree, diff = K.same_hole_sets(hs_seed, hs_part)
    partitions, multi, uncov = K.holes_partition_edges(c, hs_part)
    forest, cyc = K.split_subgraph_is_forest(c)

    res = K.hole_residuals(c, m.X, hs_part)
    mode = K.Fixed
    if bnd == "periodic"
        mode = K.Periodic
    elseif bnd == "none"
        mode = K.None
    end
    sys = K.assemble_system(c, hs_part, m.X, mode)
    rep = K.solve_system(sys, m.X)
    X0 = K.matrix_to_points(rep.X0)
    res0 = K.hole_residuals(c, X0, hs_part)
    tm = K.theta_max(c, X0)

    # Geometric cross-check runs on the projected (deployable) embedding: on a
    # non-deployable embedding the forward kinematics is inconsistent and the traced
    # complement components are meaningless.
    geo = K.holes_geometric(c, K.deploy(c, X0, 0.25).Y)
    comb = sort([hs_part.all[i].edges for i in hs_part.interior_indices])
    geo_agree = (geo == comb)

    j = Dict{String,Any}()
    j["input"] = in_
    j["N"] = K.n_vertices(m)
    j["F"] = K.n_faces(m)
    j["E"] = K.n_edges(m)
    j["n_interior"] = K.n_interior_vertices(m)
    j["n_hinge"] = K.n_hinge(c)
    j["n_split"] = K.n_split(c)
    j["n_border"] = length(c.border_edges)
    j["n_prime_vertices"] = c.n_prime_vertices
    j["components_of_Mprime"] = K.count_components(c)
    j["remark_A1_ok"] = a1
    j["remark_A1_bad_vertices"] = zero_based(bad)
    j["H_all_preimages"] = length(hs_part.all)
    j["H"] = K.n_interior_holes(hs_part)
    j["holes_seed_vs_partition_agree"] = holes_agree
    j["holes_diff"] = diff
    j["preimages_partition_edges"] = partitions
    j["multiply_covered_edges"] = zero_based(multi)
    j["uncovered_edges"] = zero_based(uncov)
    j["split_subgraph_is_forest"] = forest
    j["H_geometric"] = length(geo)
    j["holes_geometric_agree"] = geo_agree
    j["residual_max_initial"] = res.max_norm
    j["deployable_initial"] = K.deployable(res)
    j["rank_full"] = rep.rank_full
    j["rank_L"] = rep.rank_L
    j["dim_null"] = rep.dim_null
    j["n_interior_minus_H"] = K.n_interior_vertices(m) - K.n_interior_holes(hs_part)
    j["rank_claim_holds"] = (rep.dim_null == K.n_interior_vertices(m) - K.n_interior_holes(hs_part))
    j["projection_ok"] = rep.projection_ok
    j["residual_max_after_projection"] = res0.max_norm
    j["theta_max_geometric"] = tm.theta_max_geometric
    j["min_beta"] = tm.min_beta
    write_json(j, outdir * "/summary.json")

    holes_j = [Dict("edges" => zero_based(h.edges), "vertices" => zero_based(h.vertices),
                    "all_interior" => h.all_interior) for h in hs_part.all]
    write_json(holes_j, outdir * "/hole_preimages.json")

    cut_j = Dict{String,Any}()
    cut_j["edge_type"] = [c.edge_type[e] == K.Hinge ? "hinge" :
                          (c.edge_type[e] == K.Split ? "split" : "border") for e in 1:K.n_edges(m)]
    cut_j["hinge"] = [Dict("edge" => e - 1, "src" => c.hinge_dir[e].src - 1,
                           "dst" => c.hinge_dir[e].dst - 1) for e in c.hinge_edges]
    cut_j["prime_faces"] = [zero_based(f) for f in c.prime_faces]
    cut_j["prime_to_original"] = zero_based(c.prime_to_original)
    write_json(cut_j, outdir * "/cut.json")

    let m0 = K.Mesh(X0, m.faces)
        m0.sigma = m.sigma
        write_json(mesh_json(m0), outdir * "/X0.json")
        write_json(mesh_json(m), outdir * "/M.json")
    end
    for th in (0.0, 30.0 * pi / 180.0, 60.0 * pi / 180.0, max(1e-3, tm.theta_max_geometric))
        write_json(deployment_json(c, X0, th), outdir * @sprintf("/deploy_%.4f.json", th))
    end

    if do_collision && rep.dim_null > 0
        cr = K.optimize_collision(c, X0, rep.Phi)
        cj = Dict{String,Any}("f_before" => cr.f0, "f_after" => cr.f1,
                              "theta_max_before" => cr.theta_max_before,
                              "theta_max_after" => cr.theta_max_after,
                              "iterations" => cr.iterations)
        write_json(cj, outdir * "/collision.json")
        mc = K.Mesh(cr.X_opt, m.faces)
        mc.sigma = m.sigma
        write_json(mesh_json(mc), outdir * "/X_collision_opt.json")
        println("collision opt: theta_max ", cpp_g(cr.theta_max_before), " -> ",
                cpp_g(cr.theta_max_after))
    end

    open(outdir * "/summary.csv", "w") do csv
        print(csv, "N,F,E,n_interior,n_hinge,n_split,H,H_geometric,rank_full,rank_L,dim_null,",
              "n_interior_minus_H,deployable,residual_max,theta_max,min_beta\n")
        print(csv, K.n_vertices(m), ",", K.n_faces(m), ",", K.n_edges(m), ",",
              K.n_interior_vertices(m), ",", K.n_hinge(c), ",", K.n_split(c), ",",
              K.n_interior_holes(hs_part), ",", length(geo), ",", rep.rank_full, ",",
              rep.rank_L, ",", rep.dim_null, ",",
              K.n_interior_vertices(m) - K.n_interior_holes(hs_part), ",",
              cpp_b(K.deployable(res)), ",", cpp_g(res.max_norm), ",",
              cpp_g(tm.theta_max_geometric), ",", cpp_g(tm.min_beta), "\n")
    end

    println(JSON.json(_json_clean(j), 2))
    return 0
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
