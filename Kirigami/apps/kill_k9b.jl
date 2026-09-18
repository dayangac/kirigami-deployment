# K9b -- pushing the convexity-constrained embedding of K9.
#
# K9 reached exact Theta_max > 0 on 36 / 400 designs (9.0 %) against a 10 % bar, with every
# comparable baseline at exactly 0. K9b keeps the population, the shape space and the
# certificate FIXED and pushes only the solver, to see whether the 9 % is a property of the
# design space or of the search. Four levers, in the order the task ranks them:
#
#   (i)   per-graph best of sigma_mc / sigma_def. Nothing new is computed for this -- both
#         sigma are already run on every graph -- but the count is reported as an honest
#         "best-of-2 over graphs" alongside the per-design count, never mixed with it.
#   (ii)  more restarts (8 Gaussian starts + t = 0, plus the two warm starts K9 already
#         used: variant (a)'s convex minimiser and K6's split-only repair point) and more
#         barrier stages (10 vs 6).
#   (iii) a sweep of the convexity margin delta and the split margin delta' over
#         {1e-4, 1e-3, 3e-3, 1e-2} * med^2, tried in order of K9's measured success and cut
#         short as soon as a design is certified past the eps target.
#   (iv)  after feasibility, a second stage that MAXIMISES the certified range along the
#         null space: range_opt.jl's softmin-of-first-contact objective (T6's analytic
#         gradients), re-centred at the feasible point, then GATED -- the stage-2 point is
#         accepted only if the exact convexity and split constraints still hold there and
#         the exact Theta_max actually improved. The barriers are therefore enforced by
#         exact rejection, not inside the stage-2 objective; that is weaker than the task's
#         wording and is reported as such.
#
# REFEREE. K9's post-hoc audit showed the bisection's early-exit
# has_collision(1e-7) firing NON-MONOTONELY at the historical shrink 1e-12 on hinge-adjacent
# faces that touch by construction. Every design here is refereed at BOTH shrinks and both
# counts are reported; the headline count is the 1e-9 one the task names.
#
# SOUNDNESS. As in K9: "feasible" is always the exact constraint values at the returned
# point, so a positive is a true positive; an infeasible answer is a statement about this
# solver and never a proof that the convex slice is empty.
#
#   julia --project=Kirigami Kirigami/apps/kill_k9b.jl [--n 200] [--maxf 800] [--out DIR]
#         [--sigma DIR] [--cache DIR] [--eps 0.3] [--eps-target 0.1] [--iters 600]
#         [--starts 8] [--barrier 10] [--budget 150] [--shard S] [--nshards M]
#         [--no-stage2] [--no-gallery] [--aggregate] [--limit K] [--regenerate]
#
# Population as kill_k9.jl (frozen k1a_200 + archived sigma_def, shared shape-cache
# layout).
# NOTE the per-design wall budget (`--budget`) cuts the delta sweep on elapsed time, so a
# slower or faster machine can legitimately stop the sweep at a different delta. Outputs go
# to results/kill/k9b/.
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# The referee, at a caller-chosen shrink. shrink = 1e-12 is the historical setting every
# earlier kill experiment used; 1e-9 is the one K9's audit showed to be free of the
# hinge-vertex artefact.
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
    theta_ref12::Float64 = 0.0   # bisection at the historical shrink 1e-12
    theta_ref9::Float64 = 0.0    # bisection at shrink 1e-9
    zero_range::Bool = false
end

# K6's / K9's `certify`, unchanged in substance: the certificate predicate is monotone in
# eps, so the largest certified eps is the first admissible deflated root.
function certify(c::K.CutStructure, X::Vector{Vec2}, eps::Float64, referee::Bool = true)
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
    r.zero_range = orr.zero_range
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

# The exact constraint values at an arbitrary X, so that a stage-2 point produced by a
# solver that knows nothing about them can still be gated on them.
Base.@kwdef mutable struct ExactFeas
    ok::Bool = false
    min_cross::Float64 = 0.0; min_q::Float64 = 0.0
    n_bad::Int = 0; n_bad_q::Int = 0
end
function exact_feas(m::K.Mesh, c::K.CutStructure, X::Vector{Vec2}, delta::Float64, delta_q::Float64)
    f = ExactFeas(min_cross = Inf, min_q = Inf)
    for v in K.corner_crosses(m, X)
        f.min_cross = min(f.min_cross, v)
        v < delta && (f.n_bad += 1)
    end
    q = K.zero_plus_q(c, X)
    for v in q
        f.min_q = min(f.min_q, v)
        v < delta_q && (f.n_bad_q += 1)
    end
    isempty(q) && (f.min_q = 0.0)
    f.ok = (f.n_bad == 0 && f.n_bad_q == 0)
    return f
end

# One design's best answer over the whole delta sweep.
Base.@kwdef mutable struct Best
    have::Bool = false
    delta_rel::Float64 = 0.0
    feasible::Bool = false
    stage2::Int = 0              # 1 if the accepted point came from range_opt
    min_cross_rel::Float64 = 0.0; min_q_rel::Float64 = 0.0; min_mu_rel::Float64 = 0.0
    n_bad::Int = 0; n_bad_q::Int = 0; n_bad_mu::Int = 0; n_inv::Int = 0
    dist::Float64 = 0.0; tnorm::Float64 = 0.0
    barrier::Int = 0; start::Int = -1
    cert::CertReport = CertReport()
    binding::String = "n/a"
    X::Vector{Vec2} = Vec2[]
end

Base.@kwdef mutable struct VariantStats
    n::Int = 0; feasible::Int = 0; theta_pos::Int = 0; ref12::Int = 0; ref9::Int = 0; certified::Int = 0
    eps_ge_01::Int = 0; mu_ok::Int = 0; stage2::Int = 0
    theta::Vector{Float64} = Float64[]; eps_max::Vector{Float64} = Float64[]; dist::Vector{Float64} = Float64[]
    binding::Dict{String,Int} = Dict{String,Int}()
    best_delta::Dict{Int,Int} = Dict{Int,Int}()   # delta_rel * 1e6 -> wins
    sound::Int = 0; sound_n::Int = 0
    worst_gap9::Float64 = 0.0
end

function accumulate!(S::VariantStats, b::Best)
    S.n += 1
    if b.feasible
        S.feasible += 1
        push!(S.dist, b.dist)
        b.n_bad_mu == 0 && (S.mu_ok += 1)
    end
    push!(S.theta, b.cert.theta_exact)
    if b.cert.theta_exact > 1e-9
        S.theta_pos += 1
        b.cert.theta_ref12 > 1e-9 && (S.ref12 += 1)
        b.cert.theta_ref9 > 1e-9 && (S.ref9 += 1)
        k = trunc(Int, b.delta_rel * 1e6 + 0.5)
        S.best_delta[k] = get(S.best_delta, k, 0) + 1
        b.stage2 != 0 && (S.stage2 += 1)
    end
    if b.cert.eps_max > 0
        S.certified += 1
        push!(S.eps_max, b.cert.eps_max)
        b.cert.eps_max >= 0.1 && (S.eps_ge_01 += 1)
    end
    S.binding[b.binding] = get(S.binding, b.binding, 0) + 1
    S.sound_n += 1
    b.cert.eps_max <= b.cert.theta_ref9 + 1e-9 && (S.sound += 1)
    S.worst_gap9 = max(S.worst_gap9, abs(b.cert.theta_exact - b.cert.theta_ref9))
end

function report_variant(o::IO, name::String, S::VariantStats)
    print(o, name, ": designs ", S.n, ", convex+split feasible ", S.feasible, " (",
          fx(100.0 * S.feasible / max(1, S.n), 1), " %)\n")
    print(o, "  exact Theta_max > 0: ", S.theta_pos, " (",
          fx(100.0 * S.theta_pos / max(1, S.n), 1), " %)",
          ", refereed 1e-12: ", S.ref12, ", refereed 1e-9: ", S.ref9, "\n")
    print(o, "  certified (eps_max > 0): ", S.certified, ", of which eps_max >= 0.1 rad: ",
          S.eps_ge_01, "\n")
    print(o, "  all 0+ corner margins > 0 at the returned X: ", S.mu_ok, " / ", S.feasible,
          " feasible\n")
    print(o, "  positives whose accepted point came from stage 2 (range_opt): ", S.stage2, "\n")
    isempty(S.eps_max) || print(o, "  eps_max: median ", sci(median_of(S.eps_max)), ", q90 ",
                                sci(quantile_of(S.eps_max, 0.9)), ", max ",
                                sci(quantile_of(S.eps_max, 1.0)), "\n")
    isempty(S.dist) || print(o, "  |X - X_ini| per vertex, median edges: median ", fx(median_of(S.dist)),
                             ", q90 ", fx(quantile_of(S.dist, 0.9)), "\n")
    print(o, "  winning delta_rel among positives:")
    for k in sort!(collect(keys(S.best_delta)))
        print(o, " ", sci(k / 1e6, 1), "=", S.best_delta[k])
    end
    print(o, "\n")
    print(o, "  referee: worst |exact - bisection(1e-9)| ", sci(S.worst_gap9),
          "; sound (eps_max <= bisection 1e-9) ", S.sound, " / ", S.sound_n, "\n")
    print(o, "  binding at Theta_max = 0:")
    for k in sort!(collect(keys(S.binding)))
        print(o, " ", k, "=", S.binding[k])
    end
    print(o, "\n")
end

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

function fill_best_from_point!(best::Best, c::K.CutStructure, X::Vector{Vec2}, s2::Float64)
    best.X = X
    best.n_bad_mu = 0
    mmu = Inf
    mu = K.zero_plus_corner_margin(c, X)
    for v in mu
        mmu = min(mmu, v)
        v <= 0 && (best.n_bad_mu += 1)
    end
    best.min_mu_rel = isempty(mu) ? NaN : mmu / s2
end

function aggregate(outdir::String, nshards::Int)
    # Merge shard CSVs into k9b.csv, then re-derive every statistic from the merged file
    # so that the summary is a function of the CSV and of nothing else.
    out = open(joinpath(outdir, "k9b.csv"), "w")
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
    lines = readlines(joinpath(outdir, "k9b.csv"))
    head = isempty(lines) ? String[] : split(lines[1], ',')
    col(nm) = something(findfirst(==(nm), head), -1)
    ag = Dict{String,VariantStats}()
    best_by_graph = Dict{String,Float64}()   # "kind_id" -> best exact Theta_max
    beps_by_graph = Dict{String,Float64}()
    rows = 0
    for line in lines[2:end]
        isempty(line) && continue
        f = split(line, ',')
        length(f) != length(head) && continue
        rows += 1
        b = Best()
        b.have = true
        b.delta_rel = parse(Float64, f[col("best_delta")])
        b.feasible = f[col("feasible")] == "1"
        b.stage2 = parse(Int, f[col("stage2")])
        b.n_bad_mu = parse(Int, f[col("nbadmu")])
        b.dist = parse(Float64, f[col("dist")])
        b.cert.eps_max = parse(Float64, f[col("eps_max")])
        b.cert.theta_exact = parse(Float64, f[col("theta_exact")])
        b.cert.theta_ref12 = parse(Float64, f[col("theta_ref12")])
        b.cert.theta_ref9 = parse(Float64, f[col("theta_ref9")])
        b.binding = String(f[col("binding")])
        accumulate!(get!(VariantStats, ag, "k9b"), b)
        accumulate!(get!(VariantStats, ag, "k9b:" * f[col("sigma")]), b)
        gk = f[col("kind")] * "_" * f[col("id")]
        best_by_graph[gk] = max(get(best_by_graph, gk, 0.0), b.cert.theta_exact)
        beps_by_graph[gk] = max(get(beps_by_graph, gk, 0.0), b.cert.eps_max)
    end
    g_pos = count(v -> v > 1e-9, values(best_by_graph))
    g_eps = count(v -> v >= 0.1, values(beps_by_graph))

    so = open(joinpath(outdir, "summary.txt"), "w")
    print(so, "K9b -- pushing the convexity-constrained embedding\n")
    print(so, "merged from ", nshards, " shards, ", rows, " designs\n\n")
    for k in sort!(collect(keys(ag)))
        report_variant(so, k, ag[k])
    end
    print(so, "\nbest-of-2 sigma, per GRAPH (", length(best_by_graph), " graphs):\n")
    print(so, "  graphs with exact Theta_max > 0 under at least one sigma: ", g_pos, "\n")
    print(so, "  graphs with certified eps_max >= 0.1 rad: ", g_eps, "\n")
    close(so)
    println("merged into ", outdir, "/k9b.csv and wrote ", outdir, "/summary.txt")
    return 0
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n = 200; maxf = 800; shard = 0; nshards = 1
    outdir = joinpath(REPO, "results", "kill", "k9b")
    sigmadir = ""
    cache = joinpath(REPO, "results", "kill", "k6", "cache")
    eps = 0.3; eps_target = 0.1
    max_iter = 600; n_random = 8; barrier_stages = 10
    budget_s = 150.0   # per-design wall budget for the delta sweep
    do_stage2 = true; do_aggregate = false; dump_gallery = true
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--maxf" && i < length(args); maxf = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--sigma" && i < length(args); sigmadir = args[i+1]; i += 2
        elseif a == "--cache" && i < length(args); cache = args[i+1]; i += 2
        elseif a == "--eps" && i < length(args); eps = arg_f(args[i+1]); i += 2
        elseif a == "--eps-target" && i < length(args); eps_target = arg_f(args[i+1]); i += 2
        elseif a == "--iters" && i < length(args); max_iter = arg_i(args[i+1]); i += 2
        elseif a == "--starts" && i < length(args); n_random = arg_i(args[i+1]); i += 2
        elseif a == "--barrier" && i < length(args); barrier_stages = arg_i(args[i+1]); i += 2
        elseif a == "--budget" && i < length(args); budget_s = arg_f(args[i+1]); i += 2
        elseif a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshards" && i < length(args); nshards = arg_i(args[i+1]); i += 2
        elseif a == "--no-stage2"; do_stage2 = false; i += 1
        elseif a == "--no-gallery"; dump_gallery = false; i += 1
        elseif a == "--aggregate"; do_aggregate = true; i += 1
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end

    # The delta sweep, in the order K9's measurement makes most promising first.
    deltas = [1e-3, 3e-3, 1e-2, 1e-4]

    mkpath(outdir)
    dump_gallery && mkpath(joinpath(outdir, "gallery"))

    do_aggregate && return aggregate(outdir, nshards)

    csv_path = nshards > 1 ? joinpath(outdir, "shard_" * string(shard) * ".csv") :
                             joinpath(outdir, "k9b.csv")
    csv = open(csv_path, "w")
    print(csv, "id,kind,sigma,N,F,med_edge,n_split,dim_null,n_corner,nonconvex0,",
          "best_delta,feasible,stage2,min_cross,min_q,min_mu,nbad,nbadq,nbadmu,ninv,",
          "dist,tnorm,barrier,eps_max,theta_exact,theta_ref12,theta_ref9,binding,",
          "n_delta_tried,secs\n")
    wall = Timer()

    stat = Dict{String,VariantStats}()
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

            nonconvex0 = count(v -> v <= 0, K.corner_crosses(m, X0))

            seed = (UInt32(9000) + UInt32(7) * UInt32(id) + UInt32(which)) % UInt32

            # ---- the two warm starts, computed once and reused across the delta sweep -------
            oa = K.ConvexEmbedOptions()
            oa.delta_rel = 1e-3
            oa.n_random = 4
            oa.max_iter = max_iter
            oa.barrier_stages = barrier_stages
            oa.seed = seed
            ra = K.convex_embed(c, X0, sh.Phi, X_ini, med, oa)

            ro = K.ZeroPlusRepairOptions()
            ro.n_random = 4
            ro.max_iter = max_iter
            ro.seed = seed
            sp = K.zero_plus_repair(c, X0, sh.Phi, med, ro)

            # ---- lever (iii): the delta sweep --------------------------------------------
            best = Best()
            n_delta = 0
            for drel in deltas
                n_delta += 1
                ob = K.ConvexEmbedOptions()
                ob.delta_rel = drel
                ob.split_delta_rel = drel
                ob.n_random = n_random          # lever (ii): 8 Gaussian starts + t = 0
                ob.max_iter = max_iter
                ob.barrier_stages = barrier_stages
                ob.seed = (seed + UInt32(101) * UInt32(n_delta)) % UInt32
                rb = K.convex_embed(c, X0, sh.Phi, X_ini, med, ob)
                if !rb.feasible && ra.feasible
                    o2 = deepcopy(ob)
                    o2.t_init = ra.t
                    o2.n_random = 2
                    r2 = K.convex_embed(c, X0, sh.Phi, X_ini, med, o2)
                    (r2.feasible || r2.n_bad + r2.n_bad_q < rb.n_bad + rb.n_bad_q) && (rb = r2)
                end
                if !rb.feasible && sp.feasible
                    o3 = deepcopy(ob)
                    o3.t_init = sp.t
                    o3.n_random = 2
                    r3 = K.convex_embed(c, X0, sh.Phi, X_ini, med, o3)
                    (r3.feasible || r3.n_bad + r3.n_bad_q < rb.n_bad + rb.n_bad_q) && (rb = r3)
                end

                cr = certify(c, rb.X, eps)
                # The delta that wins is the one with the larger exact Theta_max, ties broken by
                # the certified eps_max. Feasibility is recorded but does not decide: K9 measured
                # 13 of its 36 positives to miss the margin and deploy anyway.
                better = !best.have ||
                         cr.theta_exact > best.cert.theta_exact + 1e-12 ||
                         (abs(cr.theta_exact - best.cert.theta_exact) <= 1e-12 &&
                          cr.eps_max > best.cert.eps_max)
                if better
                    best.have = true
                    best.delta_rel = drel
                    best.feasible = rb.feasible
                    best.stage2 = 0
                    best.min_cross_rel = rb.min_cross / s2
                    best.min_q_rel = rb.min_q / s2
                    best.n_bad = rb.n_bad
                    best.n_bad_q = rb.n_bad_q
                    best.n_inv = count_inverted(m, rb.X)
                    best.dist = rb.dist_ini
                    best.tnorm = rb.t_norm_rel
                    best.barrier = rb.barrier_stages_kept
                    best.cert = cr
                    fill_best_from_point!(best, c, rb.X, s2)
                end
                # Cut the sweep short once the design is already past the eps target, or once the
                # per-design budget is spent -- the sweep is a search, not a measurement.
                best.cert.eps_max >= eps_target && break
                s(t) > budget_s && break
            end

            # ---- lever (iv): stage 2, maximise the range, then GATE on the exact margins ----
            if do_stage2 && best.have && best.cert.theta_exact > 1e-9
                delta = best.delta_rel * s2
                rop = K.RangeOptOptions()
                rop.rounds = 8
                rop.iters_per_round = 10
                rop.barrier_eps = 0.2
                rop.step_cap = 0.25 * med
                rr = K.maximize_range(c, best.X, sh.Phi, rop)
                ef = exact_feas(m, c, rr.X_opt, delta, delta)
                if ef.ok
                    cr2 = certify(c, rr.X_opt, eps)
                    if cr2.theta_exact > best.cert.theta_exact + 1e-9
                        best.stage2 = 1
                        best.feasible = true
                        best.min_cross_rel = ef.min_cross / s2
                        best.min_q_rel = ef.min_q / s2
                        best.n_bad = ef.n_bad
                        best.n_bad_q = ef.n_bad_q
                        best.n_inv = count_inverted(m, rr.X_opt)
                        best.cert = cr2
                        fill_best_from_point!(best, c, rr.X_opt, s2)
                    end
                end
            end

            best.cert.theta_exact <= 1e-9 && (best.binding = binding_type(c, best.X))

            accumulate!(get!(VariantStats, stat, "k9b"), best)
            accumulate!(get!(VariantStats, stat, "k9b:" * fname), best)

            # ---- the gallery: closed and half-deployed, for the >= 20-graph figure ----------
            if dump_gallery && best.cert.theta_exact > 1e-9 && best.cert.theta_ref9 > 1e-9
                tag = row.kind * "_" * string(id) * "_" * fname
                write_json(deployment_json(c, best.X, 0.0), joinpath(outdir, "gallery", tag * "_closed.json"))
                write_json(deployment_json(c, best.X, 0.5 * best.cert.theta_exact),
                           joinpath(outdir, "gallery", tag * "_half.json"))
                push!(gallery, Dict("tag" => tag, "id" => id, "kind" => row.kind, "sigma" => fname,
                                    "theta_max" => best.cert.theta_exact,
                                    "eps_max" => best.cert.eps_max, "F" => K.n_faces(m)))
                gallery_n += 1
            end

            secs = s(t)
            g_ = fmt_g
            print(csv, id, ",", row.kind, ",", fname, ",", K.n_vertices(m), ",",
                  K.n_faces(m), ",", g_(med), ",", K.n_split(c), ",", sh.k, ",",
                  length(K.corner_crosses(m, X0)), ",", nonconvex0, ",",
                  g_(best.delta_rel), ",", (best.feasible ? 1 : 0), ",", best.stage2, ",",
                  g_(best.min_cross_rel), ",", g_(best.min_q_rel), ",", g_(best.min_mu_rel), ",",
                  best.n_bad, ",", best.n_bad_q, ",", best.n_bad_mu, ",",
                  best.n_inv, ",", g_(best.dist), ",", g_(best.tnorm), ",", best.barrier,
                  ",", g_(best.cert.eps_max), ",", g_(best.cert.theta_exact), ",",
                  g_(best.cert.theta_ref12), ",", g_(best.cert.theta_ref9), ",", best.binding,
                  ",", n_delta, ",", g_(secs), "\n")
            flush(csv)
            println("id ", id, " ", fname, " F=", K.n_faces(m),
                    " k=", sh.k, " delta=", g_(best.delta_rel),
                    " feas=", fmt_b(best.feasible), " s2=", best.stage2,
                    " theta=", g_(best.cert.theta_exact),
                    " eps=", g_(best.cert.eps_max), " (", fx(secs, 1), " s)")
        end
    end
    close(csv)

    dump_gallery && write_json(gallery, joinpath(outdir, "gallery", "index_" * string(shard) * ".json"))

    sum_path = nshards > 1 ? joinpath(outdir, "shard_" * string(shard) * ".txt") :
                             joinpath(outdir, "summary.txt")
    o = open(sum_path, "w")
    print(o, "K9b -- pushing the convexity-constrained embedding\n")
    print(o, "population: ", graphs, " graphs (K9's = K6's), ", designs,
          " designs (2 sigma), missing sigma_def ", missing_sigma, "\n")
    print(o, "levers: delta sweep {1e-4,1e-3,3e-3,1e-2}*med^2 (delta' = delta), ", n_random,
          " Gaussian starts + t=0 + (a)'s minimiser + the split-only repair point, ",
          barrier_stages, " barrier stages, stage-2 range_opt ",
          do_stage2 ? "ON (gated on the exact margins)" : "OFF", "\n")
    print(o, "per-design budget ", fx(budget_s, 0), " s, eps target ", fx(eps_target, 2),
          ", certificate eps = ", fx(eps, 2), "\n")
    print(o, "gallery designs dumped: ", gallery_n, "\n\n")
    for k in sort!(collect(keys(stat)))
        report_variant(o, k, stat[k])
    end
    print(o, "\nwall ", fx(s(wall), 1), " s\n")
    close(o)
    println("wrote ", csv_path, " and ", sum_path)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
