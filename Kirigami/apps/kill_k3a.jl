# K3a -- the 2-core mobility identity (ideas/ranking.md Sec. 4, K3a).
#
#   dim ker A == |F \ core2(Gamma)| + dim ker A|core2       on every graph
#   m_core = dim ker A|core2 - c(core2) <= 5                on >= 90% of the
#                                                           random Delaunay patches
#
# Extension requested by the orchestrator: the same measurement at the deployed
# configuration theta = 0.3 * theta_max, not only at theta = 0.
#
#   julia --project=Kirigami Kirigami/apps/kill_k3a.jl [--n 500] [--out DIR]
#         [--max-faces 5000] [--limit K] [--regenerate]
#
# Population: `make_graph(id, 100, max_faces, 1 << 30)` for id < n, i.e. the frozen
# data/corpus/k3a_500.json when max_faces == 5000 (the K3a run), regenerated through
# the bit-exact generators otherwise or with --regenerate. Outputs go to
# results/kill/k3a/ by default.
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n = 500
    outdir = joinpath(REPO, "results", "kill", "k3a")
    max_faces = 5000
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--max-faces" && i < length(args); max_faces = arg_i(args[i+1]); i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    csv = open(joinpath(outdir, "k3a.csv"), "w")
    print(csv, "id,kind,N,F,n_hinge,n_split,c_gamma,n_cycles,dim_kerA,n_dangling,core_faces,c_core,",
          "dim_kerA_core,identity,m_full,m_core,sigma_in_ker,sigma_res,X_source,",
          "theta_max,theta_probe,dim_kerA_th,m_full_th,m_core_th,identity_th,secs\n")
    wall = Timer()

    acc = (n = Ref(0), ident = Ref(0), ident_th = Ref(0), n_th = Ref(0), sigma_ok = Ref(0),
           del = Ref(0), del_core_le5 = Ref(0), minF = Ref(1 << 30), maxF = Ref(0))
    bad = String[]
    del_mcore = Int[]
    del_mfull = Int[]

    function run_one(name::String, kind::String, m::K.Mesh, id::Int)
        t = Timer()
        K.build_topology!(m)
        (isempty(m.sigma) || length(m.sigma) != K.n_faces(m)) && return
        c = K.make_cut(m)
        g = K.build_hinge_graph(c)
        K.n_edges(g) == 0 && return
        hs = K.holes_partition(c)
        # X: the solved Eq. (6) projection where the dense SVD is affordable, else X_ini.
        X = copy(m.X)
        src = "X_ini"
        solved = false
        if K.n_vertices(m) <= 1400
            sys = K.assemble_system(c, hs, m.X, K.Fixed)
            sr = K.solve_system(sys, m.X)
            if sr.projection_ok
                X = K.matrix_to_points(sr.X0)
                src = "X0"
                solved = true
            end
        end
        pins0 = K.pins_flat(c, g, X)
        r0 = K.mobility_at(c, g, pins0)

        # deployed probe at 0.3 * theta_max -- only meaningful where the flat state is
        # embedded and the embedding is actually deployable (else the BFS is inconsistent).
        tmax = -1.0; tprobe = -1.0
        rth = K.MobilityReport()
        did_th = false
        if solved && K.n_faces(m) <= 900
            res = K.hole_residuals(c, X, hs)
            if K.deployable(res, 1e-7) && !K.has_collision(c, K.deploy(c, X, 0.0).Y)
                tmax = K.theta_max(c, X, 40, 25).theta_max_geometric
                if tmax > 1e-6
                    tprobe = 0.3 * tmax
                    Y = K.deploy(c, X, tprobe).Y
                    rth = K.mobility_at(c, g, K.pins_deployed(c, g, Y))
                    did_th = true
                end
            end
        end
        print(csv, id, ",", kind, ",", K.n_vertices(m), ",", K.n_faces(m), ",",
              K.n_hinge(c), ",", K.n_split(c), ",", r0.components, ",", r0.n_cycles, ",",
              r0.dim_ker_A, ",", r0.n_dangling, ",", r0.core_faces, ",", r0.c_core, ",",
              r0.dim_ker_A_core, ",", fmt_b(r0.identity_holds), ",", r0.m_full, ",", r0.m_core,
              ",", fmt_b(r0.sigma_in_ker), ",", fmt_g(r0.sigma_residual), ",", src, ",", fmt_g(tmax), ",",
              fmt_g(tprobe), ",", (did_th ? rth.dim_ker_A : -1), ",", (did_th ? rth.m_full : -1),
              ",", (did_th ? rth.m_core : -1), ",", (did_th ? Int(rth.identity_holds) : -1),
              ",", fmt_g(s(t)), "\n")
        flush(csv)
        acc.n[] += 1
        acc.ident[] += r0.identity_holds
        acc.sigma_ok[] += r0.sigma_in_ker
        acc.minF[] = min(acc.minF[], K.n_faces(m))
        acc.maxF[] = max(acc.maxF[], K.n_faces(m))
        r0.identity_holds || push!(bad, name)
        if did_th
            acc.n_th[] += 1
            acc.ident_th[] += rth.identity_holds
            rth.identity_holds || push!(bad, name * "@theta")
        end
        if kind == "delaunay"
            acc.del[] += 1
            acc.del_core_le5[] += (r0.m_core <= 5)
            push!(del_mcore, r0.m_core)
            push!(del_mfull, r0.m_full)
        end
        if acc.n[] % 25 == 0
            println("  ", acc.n[], " graphs, identity ", acc.ident[], "/", acc.n[], ", ", fmt_g(s(wall)), " s")
        end
    end

    nrun = 0
    for rc in reference_cases(regenerate = regenerate)
        nrun >= limit && break
        run_one(rc.name, "reference", rc.mesh, -1)
        nrun += 1
    end
    # the population: frozen k3a_500 is exactly make_graph(id, 100, 5000, 1 << 30)
    frozen = !regenerate && max_faces == 5000
    graphs = frozen ? population("k3a_500"; first = min(n, limit)) : K.Graph[]
    for id in 0:(n - 1)
        nrun >= limit && break
        g = if frozen
            id + 1 <= length(graphs) ? graphs[id + 1] : K.make_graph(id, 100, max_faces, 1 << 30)
        else
            K.make_graph(id, 100, max_faces, 1 << 30)
        end
        g.ok || continue
        run_one(g.kind * "_" * string(id), g.kind, g.mesh, id)
        nrun += 1
    end

    med(v::Vector{Int}) = isempty(v) ? -1 : sort(v)[length(v) ÷ 2 + 1]
    sm = open(joinpath(outdir, "summary.txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("K3a: n=$(acc.n[]) graphs, F in [$(acc.minF[]),$(acc.maxF[])]\n")
    both("identity dim ker A == |F\\core2| + dim ker A|core2 : $(acc.ident[])/$(acc.n[])\n")
    both("identity at theta = 0.3*theta_max            : $(acc.ident_th[])/$(acc.n_th[])\n")
    both("sigma in ker A                               : $(acc.sigma_ok[])/$(acc.n[])\n")
    both("delaunay patches: $(acc.del[]), m_core <= 5 on $(acc.del_core_le5[]) (" *
         fx(acc.del[] > 0 ? 100.0 * acc.del_core_le5[] / acc.del[] : 0.0, 1) * "%)\n")
    both("delaunay median m_full = $(med(del_mfull)), median m_core = $(med(del_mcore))\n")
    if !isempty(del_mfull)
        both("delaunay m_full range = [$(minimum(del_mfull)),$(maximum(del_mfull))], m_core range = [" *
             "$(minimum(del_mcore)),$(maximum(del_mcore))]\n")
    end
    for b in bad
        both("VIOLATION: " * b * "\n")
    end
    both("wall " * fx(s(wall), 1) * " s\n")
    close(sm)
    close(csv)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
