# WP2b -- is the K9c failure set a solver basin or a geometric obstruction?
#
# results/yield/YIELD.md (WP2) could not predict the 93 K9c failures from input-side
# features (best held-out AUC 0.739) and left open the referee's question: "did you just
# not search hard enough?" K9b answered the analogous question for the PROXIMITY objective
# (a 4x budget made things worse), but nobody has re-seeded the RANGE objective on the
# failures. This driver does exactly that and nothing else.
#
# Population: every row of results/experiments/k9c/k9c.csv with theta_exact <= 1e-9 (93 designs),
# plus 30 successes drawn deterministically as every 10th success in row order (rows 10,
# 20, ..., 300 of the 307 successes, 1-based).
#
# Per design, method::design_range_max is run at the K9c run defaults with EIGHT new
# seeds, seed = 9300 + 7*id + which + 1000*k for k = 1..8. The graph, sigma_mc, sigma_def
# and X_ini are regenerated exactly as apps/exp_k9c_range_embedding.jl and apps/exp_yield_features.jl produce
# them: make_graph(id, 100, 800, 1400), sigma_mc = mesh.sigma, sigma_def read from
# results/experiments/k5/sigma, X_ini = mesh.X. The K9c run's own answer is copied straight
# out of k9c.csv and emitted as k = 0, so the nine seeds sit in one file and the k = 0 row
# is never a recomputation.
#
# The instrument is method::characterize, the same one K9c/WP7a use. The bisection referee
# (shrink 1e-12 through characterize, and 1e-9 by the local copy of K9c's routine) is run
# ONLY on points with exact Theta_max > 0, which is what the spec asks for and what keeps
# the cost of the 93 failure designs down.
#
#   julia --project=Kirigami Kirigami/apps/exp_basin.jl [--k9c PATH] [--sigma DIR]
#         [--out DIR] [--shard S] [--nshards M] [--seeds 8] [--limit K] [--regenerate]
#
# The graph and sigma_def come from the frozen data/corpus/native200.json (the archived
# K5 sigma_def; `--sigma DIR` reads the archived files instead), `--regenerate` rebuilds
# the graph through make_graph. The design list is read from this repo's migrated
# results/experiments/k9c/k9c.csv. Outputs go to results/yield/.
include(joinpath(@__DIR__, "exp_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# K9c's referee at a caller-chosen shrink (apps/exp_k9c_range_embedding.jl `referee_theta`, verbatim).
function referee_theta_shrink(c::K.CutStructure, X::Vector{Vec2}, shrink::Float64,
                              grid::Int = 4000, iters::Int = 50)
    col(th) = K.has_collision(c, K.deploy(c, X, th).Y, shrink)
    col(1e-7) && return 0.0
    lo = 0.0; hi = -1.0
    for i in 1:grid
        th = pi * i / grid
        if col(th)
            hi = th
            break
        end
        lo = th
    end
    hi < 0 && return Float64(pi)
    for _ in 1:iters
        mid = 0.5 * (lo + hi)
        if col(mid); hi = mid; else; lo = mid; end
    end
    return lo
end

# `std::setprecision(12) << v`.
num(v::Real) = isfinite(v) ? @sprintf("%.12g", Float64(v)) : _fmt_nonfinite(Float64(v))

# `std::getline(ss, x, ',')` over a line: a trailing empty field is dropped.
function split_csv(line::AbstractString)
    f = String.(split(line, ','))
    (!isempty(f) && isempty(f[end]) && endswith(line, ',')) && pop!(f)
    return f
end

# One selected design, with the archived K9c answer that becomes its k = 0 row.
Base.@kwdef mutable struct Design
    id::Int = -1
    which::Int = 0                 # 0 = sigma_mc, 1 = sigma_def
    kind::String = ""; sigma_name::String = ""
    cls::String = ""               # "fail" or "success"
    # the k = 0 row, copied from k9c.csv
    o_feas::String = ""; o_margin::String = ""; o_src::String = ""; o_stage_b::String = ""
    o_theta::String = ""; o_eps::String = ""; o_binding::String = ""
    o_ref12::String = ""; o_ref9::String = ""; o_secs::String = ""
end

const kHeader =
    "id,kind,sigma,which,orig_class,k,seed,status,feasible,margin,provenance,stage_b," *
    "theta_exact,eps_max,certified,binding,theta_ref12,theta_ref9,secs"

# sigma_def of one graph: the archived K5 file if `sigmadir` is given, else the frozen row
function sigma_def_of(row::PopRow, sigmadir::AbstractString)
    if !isempty(sigmadir)
        p = joinpath(sigmadir, row.kind * "_" * string(row.id) * ".json")
        isfile(p) || return Int[]
        sm = K.load_mesh_json(p)
        (length(sm.sigma) == K.n_faces(row.mesh) && length(sm.X) == length(row.mesh.X)) || return Int[]
        return sm.sigma
    end
    return row.sigma_def === nothing ? Int[] : row.sigma_def
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    k9c_csv = joinpath(REPO, "results", "experiments", "k9c", "k9c.csv")
    sigmadir = ""
    outdir = joinpath(REPO, "results", "yield")
    shard = 0; nshards = 1; nseeds = 8; limit = -1
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--k9c" && i < length(args); k9c_csv = args[i+1]; i += 2
        elseif a == "--sigma" && i < length(args); sigmadir = args[i+1]; i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshards" && i < length(args); nshards = arg_i(args[i+1]); i += 2
        elseif a == "--seeds" && i < length(args); nseeds = arg_i(args[i+1]); i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else
            println(stderr, "unknown argument ", a)
            return 2
        end
    end
    mkpath(outdir)

    # ---- the population, read from the archived K9c CSV -------------------------------
    designs = Design[]
    let
        if !isfile(k9c_csv)
            println(stderr, "cannot open ", k9c_csv)
            return 1
        end
        lines = readlines(k9c_csv)
        if isempty(lines)
            println(stderr, "empty ", k9c_csv)
            return 1
        end
        head = split_csv(lines[1])
        function col(nm)
            j = findfirst(==(nm), head)
            if j === nothing
                println(stderr, "missing column ", nm)
                exit(1)
            end
            return j
        end
        c_id = col("id"); c_kind = col("kind"); c_sig = col("sigma")
        c_th = col("theta_exact"); c_eps = col("eps_max"); c_bind = col("binding")
        c_src = col("best_src"); c_feas = col("best_feas"); c_marg = col("best_margin")
        c_sb = col("k9c_stage_b"); c_r12 = col("theta_ref12"); c_r9 = col("theta_ref9")
        c_secs = col("secs")

        fails = Design[]; succs = Design[]
        for line in lines[2:end]
            isempty(line) && continue
            v = split_csv(line)
            length(v) != length(head) && continue
            d = Design()
            d.id = parse(Int, v[c_id])
            d.kind = v[c_kind]
            d.sigma_name = v[c_sig]
            d.which = d.sigma_name == "sigma_def" ? 1 : 0
            d.o_theta = v[c_th]
            d.o_eps = v[c_eps]
            d.o_binding = v[c_bind]
            d.o_src = v[c_src]
            d.o_feas = v[c_feas]
            d.o_margin = v[c_marg]
            d.o_stage_b = v[c_sb]
            d.o_ref12 = v[c_r12]
            d.o_ref9 = v[c_r9]
            d.o_secs = v[c_secs]
            if parse(Float64, v[c_th]) <= 1e-9
                d.cls = "fail"; push!(fails, d)
            else
                d.cls = "success"; push!(succs, d)
            end
        end
        designs = copy(fails)
        # every 10th success in row order, 1-based: rows 10, 20, ... of the success list
        for j in 10:10:length(succs)
            push!(designs, succs[j])
        end
        println(stderr, "population: ", length(fails), " failures + ", length(designs) - length(fails),
                " sampled successes = ", length(designs), " designs")
    end
    # Deterministic run order, grouped by graph so one make_graph serves both sigma rows.
    sort!(designs; by = d -> (d.id, d.which), alg = MergeSort)   # stable, as std::stable_sort

    csv_path = nshards > 1 ? joinpath(outdir, "basin_shard_" * string(shard) * ".csv") :
                             joinpath(outdir, "basin.csv")
    csv = open(csv_path, "w")
    print(csv, kHeader, "\n")
    flush(csv)

    wall = Timer()
    rows = 0; done = 0
    cached_id = -1
    m0 = K.Mesh()
    sigma_mc = Int[]; sigma_def = Int[]
    pop = load_population("native200")

    for (di, d) in enumerate(designs)
        di0 = di - 1   # 0-based design index the shards are dealt by
        (nshards > 1 && di0 % nshards != shard) && continue
        (limit >= 0 && done >= limit) && break

        if d.id != cached_id
            ri = findfirst(r -> r.id == d.id, pop)
            if ri === nothing
                println(stderr, "id ", d.id, " generate FAILED")
                continue
            end
            row = pop[ri]
            if regenerate
                g = K.make_graph(d.id, 100, 800, 1400)
                if !g.ok
                    println(stderr, "id ", d.id, " generate FAILED")
                    continue
                end
                row.mesh = g.mesh
            elseif !row.ok
                println(stderr, "id ", d.id, " generate FAILED")
                continue
            end
            m0 = row.mesh
            K.build_topology!(m0)
            sigma_mc = copy(m0.sigma)
            sigma_def = sigma_def_of(row, sigmadir)
            cached_id = d.id
        end
        sigma = d.which == 1 ? sigma_def : sigma_mc
        if isempty(sigma)
            println(stderr, "id ", d.id, " ", d.sigma_name, " sigma MISSING")
            continue
        end
        done += 1

        m = K.Mesh(m0.X, m0.faces)
        m.sigma = copy(sigma)
        K.build_topology!(m)
        X_ini = copy(m0.X)
        c = K.make_cut(m)

        function row_out(k, seed, status, feas, margin, prov, stage_b, theta, eps, cert, binding, r12, r9, secs)
            print(csv, d.id, ",", d.kind, ",", d.sigma_name, ",", d.which, ",", d.cls, ",", k, ",",
                  seed, ",", status, ",", feas, ",", margin, ",", prov, ",", stage_b, ",", theta, ",",
                  eps, ",", cert, ",", binding, ",", r12, ",", r9, ",", secs, "\n")
            flush(csv)
            rows += 1
        end

        # k = 0: the archived K9c run, copied, never recomputed.
        base_seed = UInt32(9300) + UInt32(7) * UInt32(d.id) + UInt32(d.which)
        row_out(0, base_seed, "archived", d.o_feas, d.o_margin, d.o_src, d.o_stage_b, d.o_theta,
                d.o_eps, parse(Float64, d.o_eps) > 0 ? "1" : "0", d.o_binding, d.o_ref12, d.o_ref9,
                d.o_secs)

        for k in 1:nseeds
            t = Timer()
            o = K.RangeMaxOptions()           # run defaults == the archived K9c settings
            o.characterize.referee = false    # positives are refereed below, negatives never
            o.seed = base_seed + UInt32(1000) * UInt32(k)
            r = K.design_range_max(m, sigma, X_ini, o)
            if isempty(r.design.X)
                row_out(k, o.seed, isempty(r.design.status) ? "no_solution" : r.design.status,
                        "0", "0", "-", "0", "0", "0", "0", "n/a", "-1", "-1", num(s(t)))
                println(stderr, "id ", d.id, " ", d.sigma_name, " k=", k, " no solution: ", r.design.status)
                continue
            end
            ch = r.design.ch
            r12 = -1.0; r9 = -1.0
            if ch.theta_max > 1e-9
                co = K.CharacterizeOptions()
                co.referee = true             # shrink 1e-12, grid 4000, 50 bisections
                cr = K.characterize(m, sigma, r.design.X, co)
                r12 = cr.theta_bisect
                r9 = referee_theta_shrink(c, r.design.X, 1e-9)
            end
            row_out(k, o.seed, "ok", r.design.feasible ? "1" : "0", num(r.margin),
                    isempty(r.provenance) ? "-" : r.provenance, r.stage_b_used ? "1" : "0",
                    num(ch.theta_max), num(ch.eps_max), ch.eps_max > 0 ? "1" : "0", ch.binding,
                    num(r12), num(r9), num(s(t)))

            println("id ", d.id, " ", d.sigma_name, " (", d.cls, ") k=", k, " F=", K.n_faces(m),
                    " theta=", fmt_g(ch.theta_max), " eps=", fmt_g(ch.eps_max), " src=", r.provenance,
                    " [", fmt_g(s(t)), " s, wall ", fmt_g(s(wall)), "]")
            flush(stdout)
        end
    end
    close(csv)
    println("wrote ", csv_path, ": ", rows, " rows from ", done, " designs in ", fmt_g(s(wall)), " s")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
