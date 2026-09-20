# K9c -- the RANGE-MAXIMISING constrained embedding.
#
# K9 and K9b both solve the same feasibility problem inside the Eq. (6) null space (every
# face corner convex, every split cut opening outward) and then pick the point inside that
# feasible set by PROXIMITY to X_ini. K9b's audit measured the price: 9 of K9's 36
# deployable designs come back from K9b's four-times-larger search as points that are
# exactly feasible at the full margin and still have Theta_max = 0, binding vertex-edge.
# The feasible set is not what caps the yield near 9 %; the objective that selects the
# point inside it is.
#
# K9c keeps the population, the shape space, the barriers and the certificate fixed and
# replaces the objective (method/range_embed.jl):
#
#   stage A  maximise the exact 0+ margin  m(X) = min(min_e q_e, min_j mu_j) / med^2
#            through its log-sum-exp softmin surrogate (a LOWER bound on m, tight at the
#            point where the surrogate's branch modes are read), inside the same convexity
#            and split-outward log barriers, warm started from K9's point, from K9b's
#            point, and from t = 0;
#   stage B  once m > 0, push the exact Theta_max with range_opt.jl's softmin-of-first-
#            contact objective (the T6 analytic gradients, active set refreshed) under a
#            shrinking trust region, accepting a step only if the exact convexity, the
#            exact split signs and a margin floor all still hold and Theta_max improved.
#
# DETERMINISM. The per-design wall-clock budget (--budget) gates which stage-A starts
# run, so this driver's answer depends on machine load; that gate is why a replay of id
# 130 sigma_mc could land on a different local optimum (STATE.md U10). --deterministic
# turns the gate off. The library entry point design_range_max never had it, and
# reproduces this run's rows exactly at the run's settings (--stages 6 --stage-iters 120,
# seed 9300 + 7 * id + which); test_design.jl locks four of them.
#
# Every design is measured on THREE points -- K9's, K9b's and K9c's -- and the best is
# recorded WITH ITS PROVENANCE, so the K9c contribution is never confused with the
# best-of-three count. Feasibility and the margin are always the exact values at the
# returned X, never a penalty or a barrier value; an infeasible or zero-range answer is a
# statement about this solver and not about the constrained slice.
#
#   julia --project=Kirigami Kirigami/apps/exp_k9c_range_embedding.jl [--n 200] [--maxf 800] [--out DIR]
#         [--sigma DIR] [--cache DIR] [--iters 600] [--stages 8] [--stage-iters 200]
#         [--budget 120] [--shard S] [--nshards M] [--no-stage-b] [--deterministic]
#         [--no-gallery] [--aggregate] [--limit K] [--regenerate]
#
# Population as exp_k9_convex_embedding.jl. Outputs go to results/experiments/k9c/; `--aggregate` also reads
# results/experiments/k9/k9.csv and results/experiments/k9b/k9b.csv for the cross-check.
include(joinpath(@__DIR__, "exp_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# The referee, at a caller-chosen shrink; 1e-9 is the one K9's audit showed to be free of
# the hinge-vertex artefact (F34), 1e-12 the historical setting.
function referee_theta(c::K.CutStructure, X::Vector{Vec2}, shrink::Float64,
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

Base.@kwdef mutable struct CertReport
    eps_max::Float64 = 0.0
    theta_exact::Float64 = 0.0
    theta_ref12::Float64 = 0.0; theta_ref9::Float64 = 0.0
end

# K6's / K9's / K9b's `certify`, unchanged in substance. `referee` is off while candidate
# points are being compared and on once the winner is chosen: the bisection is the
# expensive half and comparing candidates does not need it.
function certify(c::K.CutStructure, X::Vector{Vec2}, referee::Bool)
    r = CertReport()
    B = K.deploy_basis(c, X)
    sd = K.swept_discs(c, B)
    pairs = K.candidate_pairs(c, sd, Float64(pi), true)

    c_pi = K.validity_certificate(c, B, X, pairs, Float64(pi))
    kRootTol = 1e-12
    eps_star = (c_pi.first_root > kRootTol) ? (c_pi.first_root - kRootTol) : 0.0
    if c_pi.first_root <= 0
        c_star = K.validity_certificate(c, B, X, pairs, Float64(pi))
        K.valid(c_star) && (r.eps_max = Float64(pi))
    end
    if c_pi.pos && eps_star > 0
        c_star = K.validity_certificate(c, B, X, pairs, eps_star)
        K.valid(c_star) && (r.eps_max = eps_star)
    end
    orr = K.exact_theta_max_overlap(c, B, pairs, 1e-9, Float64(pi), 1e-9)
    r.theta_exact = orr.theta_max
    if referee
        r.theta_ref12 = referee_theta(c, X, 1e-12)
        r.theta_ref9 = referee_theta(c, X, 1e-9)
    end
    return r
end

function binding_type(c::K.CutStructure, X::Vector{Vec2})
    any(v -> v <= 0, K.zero_plus_q(c, X)) && return "split-inward"
    m = c.mesh
    for f in 1:K.n_faces(m)
        K.polygon_area(X[m.faces[f]]) <= 0 && return "inverted"
    end
    B = K.deploy_basis(c, X)
    sd = K.swept_discs(c, B)
    pairs = K.candidate_pairs(c, sd, Float64(pi), true)
    orr = K.exact_theta_max_overlap(c, B, pairs, 1e-9, Float64(pi), 1e-9)
    return orr.zero_range ? "vertex-edge" : "hinge-wedge"
end

function median_of(v::Vector{Float64})
    isempty(v) && return NaN
    sv = sort(v)
    return sv[length(sv) ÷ 2 + 1]
end
function quantile_of(v::Vector{Float64}, p::Float64)
    isempty(v) && return NaN
    sv = sort(v)
    i = min(length(sv) - 1, trunc(Int, p * (length(sv) - 1) + 0.5))
    return sv[i + 1]
end

# One candidate point, with everything that decides whether it wins.
Base.@kwdef mutable struct Cand
    have::Bool = false
    src::String = "none"
    X::Vector{Vec2} = Vec2[]
    feasible::Bool = false
    margin::Float64 = 0.0; min_cross::Float64 = 0.0; min_q::Float64 = 0.0; min_mu::Float64 = 0.0; dist::Float64 = 0.0
    cert::CertReport = CertReport()
end

function dist_to_ini(X::Vector{Vec2}, X_ini::Vector{Vec2}, med::Float64)
    d2 = 0.0
    N = length(X)
    for v in 1:N
        d2 += K._sqnorm2(X[v] - X_ini[v])
    end
    return sqrt(d2 / max(1, N)) / med
end

Base.@kwdef mutable struct VariantStats
    n::Int = 0; feasible::Int = 0; margin_pos::Int = 0; theta_pos::Int = 0; ref12::Int = 0; ref9::Int = 0
    certified::Int = 0; eps_ge_01::Int = 0; stage_b::Int = 0
    eps_max::Vector{Float64} = Float64[]; dist::Vector{Float64} = Float64[]; margin::Vector{Float64} = Float64[]
    binding::Dict{String,Int} = Dict{String,Int}(); provenance::Dict{String,Int} = Dict{String,Int}()
    sound::Int = 0; sound_n::Int = 0
    worst_gap9::Float64 = 0.0
end

function report_variant(o::IO, name::String, S::VariantStats)
    print(o, name, ": designs ", S.n, ", convex+split feasible ", S.feasible, " (",
          fx(100.0 * S.feasible / max(1, S.n), 1), " %)",
          ", of which 0+ margin m > 0: ", S.margin_pos, "\n")
    print(o, "  exact Theta_max > 0: ", S.theta_pos, " (",
          fx(100.0 * S.theta_pos / max(1, S.n), 1), " %)",
          ", refereed 1e-12: ", S.ref12, ", refereed 1e-9: ", S.ref9, "\n")
    print(o, "  certified (eps_max > 0): ", S.certified, ", of which eps_max >= 0.1 rad: ",
          S.eps_ge_01, "\n")
    print(o, "  positives whose WINNING point came from stage B (range_opt): ", S.stage_b, "\n")
    isempty(S.eps_max) || print(o, "  eps_max: median ", sci(median_of(S.eps_max)), ", q90 ",
                                sci(quantile_of(S.eps_max, 0.9)), ", max ",
                                sci(quantile_of(S.eps_max, 1.0)), "\n")
    isempty(S.margin) || print(o, "  0+ margin m (med^2) over feasible: median ", sci(median_of(S.margin)),
                               ", q90 ", sci(quantile_of(S.margin, 0.9)), "\n")
    isempty(S.dist) || print(o, "  |X - X_ini| per vertex, median edges: median ", fx(median_of(S.dist)),
                             ", q90 ", fx(quantile_of(S.dist, 0.9)), "\n")
    if !isempty(S.provenance)
        print(o, "  provenance of the winning point:")
        for k in sort!(collect(keys(S.provenance)))
            print(o, " ", k, "=", S.provenance[k])
        end
        print(o, "\n")
    end
    print(o, "  referee: worst |exact - bisection(1e-9)| ", sci(S.worst_gap9),
          "; sound (eps_max <= bisection 1e-9) ", S.sound, " / ", S.sound_n, "\n")
    print(o, "  binding at Theta_max = 0:")
    for k in sort!(collect(keys(S.binding)))
        print(o, " ", k, "=", S.binding[k])
    end
    print(o, "\n")
end

const kHeader =
    "id,kind,sigma,N,F,med_edge,n_split,dim_null,n_corner,nonconvex0," *
    "k9_feas,k9_margin,k9_theta,k9_eps," *
    "k9b_feas,k9b_margin,k9b_theta,k9b_eps," *
    "k9c_feas,k9c_margin,k9c_theta,k9c_eps,k9c_stage_b," *
    "best_src,best_feas,best_margin,best_min_cross,best_min_q,best_min_mu," *
    "eps_max,theta_exact,theta_ref12,theta_ref9,binding,dist,secs"

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

# "<kind>_<id>_<sigma>" -> value of `col_name` from an earlier run's CSV (empty if absent)
function load_theta(path::String, col_name::String)
    out = Dict{String,Float64}()
    isfile(path) || return out
    lines = readlines(path)
    isempty(lines) && return out
    h = split(lines[1], ',')
    ci = findfirst(==(col_name), h); cid = findfirst(==("id"), h)
    ck = findfirst(==("kind"), h); cs = findfirst(==("sigma"), h)
    (ci === nothing || cid === nothing || ck === nothing || cs === nothing) && return out
    for ln in lines[2:end]
        isempty(ln) && continue
        v = split(ln, ',')
        length(v) < ci && continue
        out[v[ck] * "_" * v[cid] * "_" * v[cs]] = parse(Float64, v[ci])
    end
    return out
end

function aggregate(outdir::String, nshards::Int)
    out = open(joinpath(outdir, "k9c.csv"), "w")
    header = false
    for sh in 0:nshards-1
        p = joinpath(outdir, "shard_" * string(sh) * ".csv")
        if !isfile(p)
            println(stderr, "missing ", p)
            continue
        end
        for (i, line) in enumerate(readlines(p))
            if i == 1
                if !header
                    print(out, line, "\n")
                    header = true
                end
                continue
            end
            isempty(line) || print(out, line, "\n")
        end
    end
    close(out)

    # Re-derive every statistic FROM THE MERGED CSV, so the summary is a function of the
    # file on disk and of nothing the shards kept in memory.
    lines = readlines(joinpath(outdir, "k9c.csv"))
    head = isempty(lines) ? String[] : split(lines[1], ',')
    col(nm) = something(findfirst(==(nm), head), -1)
    ag = Dict{String,VariantStats}()
    best_by_graph = Dict{String,Float64}(); beps_by_graph = Dict{String,Float64}()
    key_theta = Dict{String,Float64}(); key_theta_k9c = Dict{String,Float64}()
    rows = 0
    for line in lines[2:end]
        isempty(line) && continue
        f = split(line, ',')
        length(f) != length(head) && continue
        rows += 1
        sig = String(f[col("sigma")])
        # --- the best-of-three row, and the same row restricted to each single solver ----
        arm(pre) = (theta = parse(Float64, f[col(pre * "_theta")]),
                    eps = parse(Float64, f[col(pre * "_eps")]),
                    margin = parse(Float64, f[col(pre * "_margin")]),
                    feas = f[col(pre * "_feas")] == "1")
        th_best = parse(Float64, f[col("theta_exact")])
        eps_best = parse(Float64, f[col("eps_max")])
        r9 = parse(Float64, f[col("theta_ref9")])
        r12 = parse(Float64, f[col("theta_ref12")])
        binding = String(f[col("binding")])
        src = String(f[col("best_src")])

        function add(name, r, is_best)
            S = get!(VariantStats, ag, name)
            S.n += 1
            if r.feas
                S.feasible += 1
                push!(S.margin, r.margin)
                r.margin > 0 && (S.margin_pos += 1)
            end
            if r.theta > 1e-9
                S.theta_pos += 1
                if is_best
                    r12 > 1e-9 && (S.ref12 += 1)
                    r9 > 1e-9 && (S.ref9 += 1)
                else  # per-arm rows are refereed only when they win
                    S.ref12 += 1
                    S.ref9 += 1
                end
            end
            if r.eps > 0
                S.certified += 1
                push!(S.eps_max, r.eps)
                r.eps >= 0.1 && (S.eps_ge_01 += 1)
            end
            if is_best
                S.binding[binding] = get(S.binding, binding, 0) + 1
                S.provenance[src] = get(S.provenance, src, 0) + 1
                push!(S.dist, parse(Float64, f[col("dist")]))
                (length(src) > 2 && endswith(src, "+B") && r.theta > 1e-9) && (S.stage_b += 1)
                S.sound_n += 1
                r.eps <= r9 + 1e-9 && (S.sound += 1)
                S.worst_gap9 = max(S.worst_gap9, abs(r.theta - r9))
            end
        end
        bestrow = (theta = th_best, eps = eps_best, margin = parse(Float64, f[col("best_margin")]),
                   feas = f[col("best_feas")] == "1")
        add("k9c:best-of-3", bestrow, true)
        add("k9c:best-of-3:" * sig, bestrow, true)
        add("arm k9  (proximity, K9 settings)", arm("k9"), false)
        add("arm k9b (proximity, delta = 1e-2)", arm("k9b"), false)
        add("arm k9c (range-maximising)", arm("k9c"), false)
        add("arm k9c (range-maximising):" * sig, arm("k9c"), false)

        gk = f[col("kind")] * "_" * f[col("id")]
        best_by_graph[gk] = max(get(best_by_graph, gk, 0.0), th_best)
        beps_by_graph[gk] = max(get(beps_by_graph, gk, 0.0), eps_best)
        key_theta[gk * "_" * sig] = th_best
        key_theta_k9c[gk * "_" * sig] = parse(Float64, f[col("k9c_theta")])
    end
    # ---- the K9 / K9b cross-check: which designs each earlier run had, and which of
    # K9b's nine LOST-but-feasible designs this run gets back. Read straight from the two
    # earlier CSVs so the claim is a function of files on disk.
    th_k9 = load_theta(joinpath(REPO, "results", "experiments", "k9", "k9.csv"), "b_theta_exact")
    th_k9b = load_theta(joinpath(REPO, "results", "experiments", "k9b", "k9b.csv"), "theta_exact")

    g_pos = count(v -> v > 1e-9, values(best_by_graph))
    g_eps = count(v -> v >= 0.1, values(beps_by_graph))

    so = open(joinpath(outdir, "summary.txt"), "w")
    print(so, "K9c -- the range-maximising constrained embedding\n")
    print(so, "merged from ", nshards, " shards, ", rows, " designs\n")
    print(so, "NOTE: the per-arm rows are counted from the same CSV; only the best-of-3 row is\n",
          "refereed by bisection, so the per-arm 'refereed' columns repeat the exact scan.\n\n")
    for k in sort!(collect(keys(ag)))
        report_variant(so, k, ag[k])
    end
    print(so, "\nbest-of-2 sigma, per GRAPH (", length(best_by_graph), " graphs):\n")
    print(so, "  graphs with exact Theta_max > 0 under at least one sigma: ", g_pos, "\n")
    print(so, "  graphs with certified eps_max >= 0.1 rad: ", g_eps, "\n")

    if !isempty(th_k9) && !isempty(th_k9b)
        lost = 0; recovered_c = 0; recovered_best = 0
        lost_names = String[]; rec_names = String[]
        for k in sort!(collect(keys(th_k9)))   # std::map order
            v = th_k9[k]
            (v <= 1e-9 || !haskey(th_k9b, k) || th_k9b[k] > 1e-9) && continue
            lost += 1                                    # K9 positive, K9b zero
            push!(lost_names, k)
            if haskey(key_theta_k9c, k) && key_theta_k9c[k] > 1e-9
                recovered_c += 1
                push!(rec_names, k)
            end
            (haskey(key_theta, k) && key_theta[k] > 1e-9) && (recovered_best += 1)
        end
        print(so, "\nK9 / K9b cross-check (from results/experiments/k9/k9.csv and k9b/k9b.csv):\n")
        print(so, "  K9 positives ", count(v -> v > 1e-9, values(th_k9)),
              ", K9b positives ", count(v -> v > 1e-9, values(th_k9b)), "\n")
        print(so, "  designs K9 had and K9b LOST: ", lost, "\n")
        print(so, "  ... recovered by the K9c arm alone: ", recovered_c, "\n")
        print(so, "  ... recovered by the best-of-3 row: ", recovered_best, "\n")
        print(so, "  lost designs:")
        for s2 in lost_names
            print(so, " ", s2)
        end
        print(so, "\n  recovered by K9c:")
        for s2 in rec_names
            print(so, " ", s2)
        end
        print(so, "\n")
    end
    close(so)
    println("merged into ", outdir, "/k9c.csv and wrote ", outdir, "/summary.txt")
    return 0
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n = 200; maxf = 800; shard = 0; nshards = 1
    outdir = joinpath(REPO, "results", "experiments", "k9c")
    sigmadir = ""
    cache = joinpath(REPO, "results", "experiments", "k6", "cache")
    max_iter = 600; stages = 8; iter_per_stage = 200
    budget_s = 120.0
    do_stage_b = true; do_aggregate = false; dump_gallery = true
    # The per-design wall-clock budget below decides WHICH stage-A starts run, so the
    # answer depends on machine load: this is the one source of non-determinism in the K9c
    # pipeline (STATE.md U10). --deterministic disables the three gates and always runs
    # every start, which is what design_range_max does.
    deterministic = false
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--maxf" && i < length(args); maxf = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--sigma" && i < length(args); sigmadir = args[i+1]; i += 2
        elseif a == "--cache" && i < length(args); cache = args[i+1]; i += 2
        elseif a == "--iters" && i < length(args); max_iter = arg_i(args[i+1]); i += 2
        elseif a == "--stages" && i < length(args); stages = arg_i(args[i+1]); i += 2
        elseif a == "--stage-iters" && i < length(args); iter_per_stage = arg_i(args[i+1]); i += 2
        elseif a == "--budget" && i < length(args); budget_s = arg_f(args[i+1]); i += 2
        elseif a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshards" && i < length(args); nshards = arg_i(args[i+1]); i += 2
        elseif a == "--no-stage-b"; do_stage_b = false; i += 1
        elseif a == "--deterministic"; deterministic = true; i += 1
        elseif a == "--no-gallery"; dump_gallery = false; i += 1
        elseif a == "--aggregate"; do_aggregate = true; i += 1
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end

    mkpath(outdir)
    dump_gallery && mkpath(joinpath(outdir, "gallery"))

    do_aggregate && return aggregate(outdir, nshards)

    csv_path = nshards > 1 ? joinpath(outdir, "shard_" * string(shard) * ".csv") :
                             joinpath(outdir, "k9c.csv")
    csv = open(csv_path, "w")
    print(csv, kHeader, "\n")
    wall = Timer()

    graphs = 0; missing_sigma = 0; designs = 0; gallery_n = 0
    gallery = Any[]

    maxf == 800 || regenerate || @warn "frozen k1a_200 was built with --maxf 800; use --regenerate for maxf = $maxf"
    gidx = -1
    nrun = 0
    for row in load_population("k1a_200")
        gidx + 1 < n || break
        row.ok || continue
        if regenerate
            g = K.make_graph(row.id, 100, maxf, 1400)
            g.ok || continue
            row.mesh = g.mesh
        end
        m0 = row.mesh
        K.build_topology!(m0)
        sigma_mc = copy(m0.sigma)
        X_ini = copy(m0.X)
        id = row.id

        sigma_def = sigma_def_of(row, sigmadir)
        if isempty(sigma_def)
            missing_sigma += 1
            continue
        end
        gidx += 1
        gidx >= n && break
        (nshards > 1 && gidx % nshards != shard) && continue
        nrun >= limit && break
        nrun += 1

        med = median_edge_length(m0)
        s2 = med * med
        counted = false

        for which in 0:1
            fname = which == 1 ? "sigma_def" : "sigma_mc"
            t = Timer()
            m = K._with_sigma(m0, which == 1 ? sigma_def : sigma_mc)
            K.build_topology!(m)
            c = K.make_cut(m)
            K.n_split(c) == 0 && continue
            hs = K.holes_partition(c)
            sh = shape_space(m, c, hs, joinpath(cache, fname), id)
            (!sh.ok || sh.k < 1) && continue
            X0 = K.matrix_to_points(sh.X0)
            if !counted
                graphs += 1
                counted = true
            end
            designs += 1

            cr0 = K.corner_crosses(m, X0)
            n_corner = length(cr0)
            nonconvex0 = count(v -> v <= 0, cr0)

            seed = (UInt32(9300) + UInt32(7) * UInt32(id) + UInt32(which)) % UInt32

            function make_cand(src::String, X::Vector{Vec2}, delta_rel::Float64)
                cd = Cand(have = true, src = src, X = X)
                mc = Inf
                nbad = 0
                for v in K.corner_crosses(m, X)
                    mc = min(mc, v)
                    v < delta_rel * s2 && (nbad += 1)
                end
                cd.min_cross = (n_corner != 0 ? mc : s2) / s2
                nbq = count(v -> v < delta_rel * s2, K.zero_plus_q(c, X))
                cd.feasible = (nbad == 0 && nbq == 0)
                cd.margin, cd.min_q, cd.min_mu = K.zero_plus_margin(c, X, med)
                cd.dist = dist_to_ini(X, X_ini, med)
                cd.cert = certify(c, X, false)
                return cd
            end

            # ---- arm 1: K9's point (proximity, delta = 1e-3, K9's search settings) ---------
            o9 = K.ConvexEmbedOptions()
            o9.delta_rel = 1e-3
            o9.split_delta_rel = 1e-3
            o9.n_random = 3
            o9.max_iter = max_iter
            o9.barrier_stages = 6
            o9.seed = seed
            r9 = K.convex_embed(c, X0, sh.Phi, X_ini, med, o9)
            c9 = make_cand("k9", r9.X, 1e-3)
            t_k9 = s(t)

            # ---- arm 2: K9b's point (proximity, delta = 1e-2, K9b's larger search) ---------
            # K9b's four-value sweep is cut to its two winning values: 27 of its 28 positives
            # won at 1e-3 or 1e-2, and the full sweep does not fit the wall budget.
            o9b = K.ConvexEmbedOptions()
            o9b.delta_rel = 1e-2
            o9b.split_delta_rel = 1e-2
            o9b.n_random = 8
            o9b.max_iter = max_iter
            o9b.barrier_stages = 10
            o9b.seed = (seed + UInt32(101)) % UInt32
            r9.feasible && (o9b.t_init = r9.t)
            r9b = K.convex_embed(c, X0, sh.Phi, X_ini, med, o9b)
            c9b = make_cand("k9b", r9b.X, 1e-2)
            t_k9b = s(t)

            # ---- arm 3: K9c, stage A -- maximise the 0+ margin from every available start ---
            c9c = Cand()
            stage_b_used = false
            let
                starts = Tuple{String,Vector{Float64},Float64}[]
                r9.feasible && push!(starts, ("k9", r9.t, 1e-3))
                r9b.feasible && push!(starts, ("k9b", r9b.t, 1e-2))
                push!(starts, ("x0", zeros(2 * sh.k), 1e-3))
                for (tag, tinit, drel) in starts
                    (!deterministic && tag == "x0" && c9c.have && c9c.feasible &&
                     s(t) > 0.4 * budget_s) && continue
                    (!deterministic && s(t) > 1.2 * budget_s && c9c.have) && break
                    ro = K.RangeEmbedOptions()
                    ro.delta_rel = drel
                    ro.split_delta_rel = drel
                    ro.stages = stages
                    ro.iter_per_stage = iter_per_stage
                    ro.seed = (seed + UInt32(211)) % UInt32
                    ro.t_init = tinit
                    rr = K.range_embed(c, X0, sh.Phi, med, ro)
                    cd = make_cand("k9c/" * tag, rr.X, drel)
                    better = !c9c.have ||
                             cd.cert.theta_exact > c9c.cert.theta_exact + 1e-12 ||
                             (abs(cd.cert.theta_exact - c9c.cert.theta_exact) <= 1e-12 &&
                              cd.margin > c9c.margin)
                    better && (c9c = cd)
                    println(stderr, "    start ", tag, " done at ", fx(s(t), 1), " s")
                    (!deterministic && c9c.cert.theta_exact > 1e-9 && s(t) > 0.5 * budget_s) && break
                end
                # ---- stage B: push Theta_max under the margin floor and the trust region ------
                if do_stage_b && c9c.have && c9c.margin > 0 && c9c.feasible
                    mo = K.MarginRangeOptions()
                    mo.delta_rel = (c9c.src == "k9c/k9b") ? 1e-2 : 1e-3
                    mo.split_delta_rel = mo.delta_rel
                    mr = K.maximize_margin_range(c, c9c.X, sh.Phi, med, mo)
                    if mr.improved
                        cd = make_cand(c9c.src * "+B", mr.X, mo.delta_rel)
                        if cd.cert.theta_exact > c9c.cert.theta_exact + 1e-9
                            c9c = cd
                            stage_b_used = true
                        end
                    end
                end
            end

            # ---- the best of the three, with provenance ------------------------------------
            t_k9c = s(t)
            best = c9
            for cd in (c9b, c9c)
                cd.have || continue
                if cd.cert.theta_exact > best.cert.theta_exact + 1e-12 ||
                   (abs(cd.cert.theta_exact - best.cert.theta_exact) <= 1e-12 &&
                    cd.cert.eps_max > best.cert.eps_max)
                    best = cd
                end
            end
            # The winner, and only the winner, is refereed at both shrinks.
            bc = certify(c, best.X, true)
            binding = "n/a"
            bc.theta_exact <= 1e-9 && (binding = binding_type(c, best.X))

            if dump_gallery && bc.theta_exact > 1e-9 && bc.theta_ref9 > 1e-9
                tag = row.kind * "_" * string(id) * "_" * fname
                write_json(deployment_json(c, best.X, 0.0), joinpath(outdir, "gallery", tag * "_closed.json"))
                write_json(deployment_json(c, best.X, 0.5 * bc.theta_exact),
                           joinpath(outdir, "gallery", tag * "_half.json"))
                push!(gallery, Dict("tag" => tag, "id" => id, "kind" => row.kind, "sigma" => fname,
                                    "theta_max" => bc.theta_exact, "eps_max" => bc.eps_max,
                                    "src" => best.src, "F" => K.n_faces(m)))
                gallery_n += 1
            end

            secs = s(t)
            g_ = fmt_g
            print(csv, id, ",", row.kind, ",", fname, ",", K.n_vertices(m), ",",
                  K.n_faces(m), ",", g_(med), ",", K.n_split(c), ",", sh.k, ",",
                  n_corner, ",", nonconvex0, ",",
                  (c9.feasible ? 1 : 0), ",", g_(c9.margin), ",", g_(c9.cert.theta_exact),
                  ",", g_(c9.cert.eps_max), ",",
                  (c9b.feasible ? 1 : 0), ",", g_(c9b.margin), ",", g_(c9b.cert.theta_exact),
                  ",", g_(c9b.cert.eps_max), ",",
                  (c9c.feasible ? 1 : 0), ",", g_(c9c.margin), ",", g_(c9c.cert.theta_exact),
                  ",", g_(c9c.cert.eps_max), ",", (stage_b_used ? 1 : 0), ",",
                  best.src, ",", (best.feasible ? 1 : 0), ",", g_(best.margin), ",",
                  g_(best.min_cross), ",", g_(best.min_q), ",", g_(best.min_mu), ",",
                  g_(bc.eps_max), ",", g_(bc.theta_exact), ",", g_(bc.theta_ref12), ",",
                  g_(bc.theta_ref9), ",", binding, ",", g_(best.dist), ",", g_(secs), "\n")
            flush(csv)
            println("id ", id, " ", fname, " F=", K.n_faces(m), " k=", sh.k,
                    " | k9 th=", fx(c9.cert.theta_exact, 3),
                    " m=", fx(c9.margin, 4),
                    " | k9b th=", fx(c9b.cert.theta_exact, 3),
                    " m=", fx(c9b.margin, 4),
                    " | k9c(", c9c.src, ") th=", fx(c9c.cert.theta_exact, 3),
                    " m=", fx(c9c.margin, 4),
                    " | best=", best.src, " (", fx(secs, 1), " s: k9 ",
                    fx(t_k9, 1), ", k9b ", fx(t_k9b - t_k9, 1), ", k9c ",
                    fx(t_k9c - t_k9b, 1), ")")
        end
    end
    close(csv)

    dump_gallery && write_json(gallery, joinpath(outdir, "gallery", "index_" * string(shard) * ".json"))

    sum_path = nshards > 1 ? joinpath(outdir, "shard_" * string(shard) * ".txt") :
                             joinpath(outdir, "summary.txt")
    o = open(sum_path, "w")
    print(o, "K9c -- the range-maximising constrained embedding\n")
    print(o, "population: ", graphs, " graphs (K9's = K6's), ", designs,
          " designs (2 sigma), missing sigma_def ", missing_sigma, "\n")
    print(o, "arms: k9 (proximity, delta=1e-3), k9b (proximity, delta=1e-2, 8 starts), ",
          "k9c (max 0+ margin from k9/k9b/x0, ", stages, " stages x ",
          iter_per_stage, " iters, stage B ", do_stage_b ? "ON" : "OFF", ")\n")
    print(o, "per-design budget ", fx(budget_s, 0), " s",
          deterministic ? " (IGNORED: --deterministic, every start always runs)" : "", "\n")
    print(o, "gallery designs dumped: ", gallery_n, "\n")
    print(o, "\nwall ", fx(s(wall), 1), " s\n")
    close(o)
    println("wrote ", csv_path, " and ", sum_path)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
