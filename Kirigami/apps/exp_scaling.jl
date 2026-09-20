# WP7b -- the SCALING BENCHMARK.
#
# MISSION.md Sec. 3 asked for 100-5,000 faces; nothing in this project above |F| = 800 has
# ever been run, so REPORT.md's objection 4 ("how does the exact range, the certificate
# and the design routine scale?") has no data behind it. This driver produces that data.
#
# Graphs: `generate("delaunay"|"voronoi", {sites, 40.0}, rng)` at the site counts that land
# near target |F| in {50, 100, 200, 500, 1000, 2000, 5000}, three seeds each, with the same
# orientation rule the experiment population uses (assign_orientation_relaxation, 4 restarts /
# 300 iters / 90 diameters). The ACTUAL |F| and |N| are what the CSV records.
#
# One process per (cell, routine), so the peak resident set at exit (Sys.maxrss, i.e.
# getrusage(RUSAGE_SELF).ru_maxrss) is that routine's own high-water mark and not one
# inherited from an earlier stage. The caller wraps each invocation in the repo's
# `perl -e alarm` cap; a routine that is killed simply leaves no row, and the merge step
# reports it as "capped".
#
#   --stage solve        assemble_system + solve_system (DENSE SVD): rank, dim_null, and
#                        the Eq. (6) projection X0. Caches X0/Phi through exp_common's
#                        save_shape so the later stages do not repeat it.
#   --stage sparse       rank_only_sparse: the rank alone, by sparse QR. This is what is
#                        left once the dense SVD is out of budget -- it gives no
#                        null-space basis, so neither characterize nor design_range_max
#                        can run from it. The CSV says so in the `note` column.
#   --stage characterize method::characterize at X0 (exact Theta_max + the eps = 0.3
#                        certificate), reporting the candidate-pair count of the exact
#                        scan.
#   --stage rangemax     method::design_range_max at run defaults. NOTE: this entry point
#                        re-solves Eq. (4) itself, so its wall time INCLUDES one dense
#                        solve; the `solve` row is what to subtract.
#
# The 1-minute load average at the start of every timed section is recorded, because this
# machine carries two detached baseline reruns throughout (MISSION note): a reader must be
# able to discount contention.
#
#   julia --project=Kirigami Kirigami/apps/exp_scaling.jl [--kind delaunay] [--sites 100]
#         [--seed 0] [--stage solve] [--out DIR] [--csv PATH] [--regenerate]
#
# The graph comes from the frozen data/corpus/scaling_42.json (the (kind, sites, seed)
# cell, sigma_mc included); `--regenerate` rebuilds it through `build_graph` instead.
# Outputs go to results/scaling/ by default.
include(joinpath(@__DIR__, "exp_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# Peak resident set in MB (Sys.maxrss reads ru_maxrss and already accounts for the
# bytes-vs-kilobytes difference between macOS and Linux).
peak_rss_mb() = Sys.maxrss() / (1024.0 * 1024.0)

load1() = Sys.loadavg()[1]

# `std::setprecision(10) << v`.
num(v::Real) = @sprintf("%.10g", Float64(v))

# The graph of one cell, a deterministic function of (kind, sites, seed).
function build_graph(kind::AbstractString, sites::Int, seed::Int, with_sigma::Bool)
    rng = K.MT19937(UInt32(1000003) * UInt32(seed) + UInt32(20260908) + UInt32(7919) * UInt32(sites))
    m = K.generate(kind, [Float64(sites), 40.0], rng)
    K.build_topology!(m)
    with_sigma && (m.sigma = K.assign_orientation_relaxation(m, rng, 4, 300, 90).sigma)
    return m
end

# The frozen cell of data/corpus/scaling_42.json, or `nothing` when the grid has no
# such (kind, sites, seed).
function frozen_graph(kind::AbstractString, sites::Int, seed::Int)
    path = joinpath(CORPUS_DIR, "scaling_42.json")
    isfile(path) || error("exp_scaling: no corpus file $path")
    for r in JSON.parsefile(path)
        (String(r["kind"]) == kind && Int(r["sites"]) == sites && Int(r["seed"]) == seed) || continue
        m = K.mesh_from_json(r["mesh"])
        m.sigma = _intvec(r["sigma_mc"])
        K.build_topology!(m)
        return m
    end
    return nothing
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    kind = "delaunay"; stage = "solve"
    outdir = joinpath(REPO, "results", "scaling")
    csv_path = ""
    sites = 100; seed = 0
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--kind" && i < length(args); kind = args[i+1]; i += 2
        elseif a == "--sites" && i < length(args); sites = arg_i(args[i+1]); i += 2
        elseif a == "--seed" && i < length(args); seed = arg_i(args[i+1]); i += 2
        elseif a == "--stage" && i < length(args); stage = args[i+1]; i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--csv" && i < length(args); csv_path = args[i+1]; i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    cache = joinpath(outdir, "shape")
    mkpath(cache)
    isempty(csv_path) && (csv_path = joinpath(outdir, "scaling.csv"))

    exists = isfile(csv_path) && filesize(csv_path) > 0
    csv = open(csv_path, exists ? "a" : "w")
    if !exists
        print(csv, "kind,sites,seed,N,F,n_split,n_hinge,routine,secs,rss_mb,load1,rank_L,",
              "dim_null,n_pairs,theta_max,eps_max,status,note\n")
        flush(csv)
    end

    tag = kind * "_" * string(sites) * "_" * string(seed)

    tgen = Timer()
    m = regenerate ? build_graph(kind, sites, seed, true) : frozen_graph(kind, sites, seed)
    if m === nothing
        println(stderr, "[scaling] no frozen cell for ", tag, " (use --regenerate)")
        return 2
    end
    if K.n_faces(m) == 0 || length(m.sigma) != K.n_faces(m)
        println(stderr, "[scaling] graph build failed for ", tag)
        return 2
    end
    K.build_topology!(m)
    c = K.make_cut(m)
    hs = K.holes_partition(c)
    gen_s = s(tgen)

    function row(routine, secs, rss, load, rank_L, dim_null, n_pairs, th, eps, status, note)
        print(csv, kind, ",", sites, ",", seed, ",", K.n_vertices(m), ",", K.n_faces(m), ",",
              K.n_split(c), ",", K.n_hinge(c), ",", routine, ",", num(secs), ",", num(rss), ",",
              num(load), ",", rank_L, ",", dim_null, ",", n_pairs, ",", num(th), ",", num(eps),
              ",", status, ",", note, "\n")
        flush(csv)
    end

    if stage == "graph"
        row("generate", gen_s, peak_rss_mb(), load1(), -1, -1, -1, 0, 0, "ok", "-")
        println(tag, " N=", K.n_vertices(m), " F=", K.n_faces(m), " split=", K.n_split(c),
                " hinge=", K.n_hinge(c), " gen ", fmt_g(gen_s), " s")
        return 0
    end

    if stage == "solve"
        l0 = load1()
        t = Timer()
        sys = K.assemble_system(c, hs, m.X, K.Fixed)
        sr = K.solve_system(sys, m.X)
        secs = s(t)
        sh = Shape(K.n_vertices(m), sr.dim_null, sr.X0, sr.Phi, sr.rank_L, sr.H,
                   K.n_interior_vertices(m), sr.projection_ok)
        if sh.ok
            mkpath(joinpath(cache, tag))
            save_shape(joinpath(cache, tag), 0, sh)
        end
        row("solve_dense", secs, peak_rss_mb(), l0, sr.rank_L, sr.dim_null, -1, 0, 0,
            sr.projection_ok ? "ok" : "projection_failed", "assemble+dense_svd")
        println(tag, " solve_dense ", fmt_g(secs), " s, rank ", sr.rank_L, ", k ", sr.dim_null,
                ", rss ", fmt_g(peak_rss_mb()), " MB")
        return 0
    end

    if stage == "sparse"
        l0 = load1()
        t = Timer()
        sys = K.assemble_system(c, hs, m.X, K.Fixed)
        sr = K.rank_only_sparse(sys)
        secs = s(t)
        row("rank_sparse", secs, peak_rss_mb(), l0, sr.rank_L, sr.dim_null, -1, 0, 0, "ok",
            "rank_only;no_null_space_basis")
        println(tag, " rank_sparse ", fmt_g(secs), " s, rank ", sr.rank_L)
        return 0
    end

    if stage == "characterize"
        sh = load_shape(joinpath(cache, tag), 0)
        if sh === nothing || sh.N != K.n_vertices(m)
            row("characterize", 0, peak_rss_mb(), load1(), -1, -1, -1, 0, 0, "no_shape",
                "dense_solve_unavailable")
            println(stderr, tag, " characterize: no cached shape")
            return 3
        end
        X0 = K.matrix_to_points(sh.X0)
        l0 = load1()
        t = Timer()
        co = K.CharacterizeOptions()   # eps = 0.3, no bisection referee (O(n^2) already)
        ch = K.characterize(m, m.sigma, X0, co)
        secs = s(t)
        row("characterize", secs, peak_rss_mb(), l0, sh.rank_L, sh.k, ch.n_pairs, ch.theta_max,
            ch.eps_max, "ok", "at_eq6_point;" * ch.binding)
        println(tag, " characterize ", fmt_g(secs), " s, pairs ", ch.n_pairs, ", theta ",
                fmt_g(ch.theta_max))
        return 0
    end

    if stage == "rangemax"
        l0 = load1()
        t = Timer()
        o = K.RangeMaxOptions()   # run defaults
        o.characterize.referee = false
        r = K.design_range_max(m, m.sigma, m.X, o)
        secs = s(t)
        row("design_range_max", secs, peak_rss_mb(), l0, -1, r.design.dim_null, -1,
            r.design.ch.theta_max, r.design.ch.eps_max, r.design.ok ? "ok" : r.design.status,
            "includes_own_solve;" * (isempty(r.provenance) ? "-" : r.provenance))
        println(tag, " design_range_max ", fmt_g(secs), " s, theta ", fmt_g(r.design.ch.theta_max),
                " (", r.provenance, ")")
        return 0
    end

    println(stderr, "unknown --stage ", stage)
    return 2
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(copy(ARGS)))
end
