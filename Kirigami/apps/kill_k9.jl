# K9 -- convexity-constrained embedding in the Tutte auxetic null space.
# Port of code/apps/kill_k9.cpp.
#
# THE HYPOTHESIS (ESCALATION.md option 2). Every kill experiment so far dies on the same
# 0+ contact, and T-1 localised it: 1,966 of the 1,968 balanced pure vertices whose 0+
# margin is non-positive sit at a REFLEX corner, and the input embeddings have 0 non-convex
# faces while the Eq. (6) projection X0 has 22,103 non-convex corners out of 156,220. At a
# convex corner the 0+ margin is mu = max(-g1, -g2), a disjunction that is easy to satisfy;
# at a reflex corner it is min(-g1, -g2), a conjunction that is not. So the obstruction may
# be MANUFACTURED by Eq. (6)'s unconstrained projection rather than inherent to the graph.
#
# WHAT IS MEASURED. Inside the same affine shape space X(t) = X0 + Phi t, solve
#
#     minimise || X(t) - X_ini ||^2   s.t.  every face corner has cross_i >= delta,
#
# (variant a: `method/convex_embed.jl`), and the same with the split-cut signs added as
# constraints q_e >= delta' (variant b). Then, at the resulting X:
#   - the exact Theta_max (T4.2''), refereed by the same bisection K2a/K5/K6 used;
#   - the repaired validity certificate (POS /\ NOOVERLAP(eps/2) /\ NOROOT(eps)) and the
#     largest certified eps;
#   - the 0+ corner margins mu, so that "did convexity buy the margin the hypothesis says
#     it should" is answered directly and not only through Theta_max;
#   - K8a's expansive-cone LP re-run at the new X (the non-uniform question): is there ANY
#     first-order expansive motion there, uniform or not.
#
# The population is K6's: the same 200 random graphs (voronoi / delaunay / quad_random)
# and the same two sigma families (sigma_mc from Eq. (1), sigma_def from K5), reusing K6's
# shape-space cache so that X0 and Phi are bit-identical to the ones K6 measured.
#
# PASS bar (task): exact Theta_max > 0 on >= 10 % of the 400 designs (>= 40), refereed by
# bisection. HARD FAIL: convex-feasible points exist on >= 50 % of the designs yet
# Theta_max > 0 stays below 2 % -- then convexity is not the lever.
#
# SOUNDNESS NOTE, stated once. `convex_embed` is a penalty/barrier heuristic on a NON-convex
# feasible set. "Feasible" is always the exact constraint values at the returned point, so a
# positive is a true positive; an infeasible answer is a statement about this solver, never a
# proof that the convex slice of the shape space is empty.
#
#   julia --project=Kirigami Kirigami/apps/kill_k9.jl [--n 200] [--maxf 800] [--out DIR]
#         [--sigma DIR] [--cache DIR] [--eps 0.3] [--delta 1e-3] [--split-delta 1e-3]
#         [--iters 600] [--starts 3] [--barrier 6] [--shard S] [--nshards M]
#         [--cone-maxf 800] [--no-cone] [--aggregate] [--limit K] [--regenerate]
#
# Population: the K1a graphs (frozen k1a_200 with the archived K5 `sigma_def`; `--sigma DIR`
# reads results/kill/k5/sigma/<kind>_<id>.json instead, `--regenerate` rebuilds the graphs).
# The shape cache is the C++ layout (results/kill/k6/cache/<sigma>/). Outputs go to
# results/kill/k9_julia/ (the C++ wrote results/kill/k9/).
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# The referee: bisection on the true polygon overlap, shrink 1e-12, 4000-point grid.
# Identical to K2a/K2b/K5/K6 so the numbers are comparable across experiments.
function referee_theta(c::K.CutStructure, X::Vector{Vec2}, grid::Int = 4000, iters::Int = 50)
    col(th) = K.has_collision(c, K.deploy(c, X, th).Y, 1e-12)
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
    valid_at_eps::Bool = false
    pos::Bool = false; noovl::Bool = false; noroot::Bool = false
    eps_max::Float64 = 0.0
    theta_exact::Float64 = 0.0
    theta_bisect::Float64 = 0.0
    zero_range::Bool = false
end

# Exactly K6's `certify`, so that the certified numbers of the two experiments are the
# same quantity. The certificate predicate is monotone in eps, so the largest certified
# eps is the first admissible deflated root (minus the root routine's own tolerance).
function certify(c::K.CutStructure, X::Vector{Vec2}, eps::Float64)
    r = CertReport()
    B = K.deploy_basis(c, X)
    sd = K.swept_discs(c, B)
    pairs = K.candidate_pairs(c, sd, Float64(pi), true)

    c_eps = K.validity_certificate(c, B, X, pairs, eps)
    r.valid_at_eps = K.valid(c_eps)
    r.pos = c_eps.pos
    r.noovl = c_eps.nooverlap
    r.noroot = c_eps.noroot

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
    r.theta_bisect = referee_theta(c, X)
    return r
end

# K6's classification of what binds at a design with Theta_max = 0.
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

# Everything measured at one returned shape-space point.
Base.@kwdef mutable struct PointReport
    feasible::Bool = false
    min_cross_rel::Float64 = 0.0   # in med^2
    min_q_rel::Float64 = 0.0
    n_bad::Int = 0; n_bad_q::Int = 0; n_nonconvex::Int = 0
    dist::Float64 = 0.0; tnorm::Float64 = 0.0
    iters::Int = 0; barrier::Int = 0
    min_mu_rel::Float64 = 0.0      # 0+ corner margin at the returned X, in med^2
    n_bad_mu::Int = 0
    n_inv::Int = 0
    cert::CertReport = CertReport()
    binding::String = "n/a"
    cone_run::Int = 0; cone_feasible::Int = 0
    cone_margin::Float64 = 0.0; cone_dual::Float64 = 0.0
    cone_sigma_in::Int = 0
end

function evaluate(c::K.CutStructure, m::K.Mesh, r::K.ConvexEmbedResult, med::Float64,
                  eps::Float64, do_cone::Bool, cone_opt::K.ExpansiveConeOptions)
    p = PointReport()
    s2 = med * med
    p.feasible = r.feasible
    p.min_cross_rel = r.min_cross / s2
    p.min_q_rel = r.min_q / s2
    p.n_bad = r.n_bad
    p.n_bad_q = r.n_bad_q
    p.n_nonconvex = r.n_nonconvex
    p.dist = r.dist_ini
    p.tnorm = r.t_norm_rel
    p.iters = r.iterations
    p.barrier = r.barrier_stages_kept
    p.n_inv = count_inverted(m, r.X)

    mu = K.zero_plus_corner_margin(c, r.X)
    mmu = Inf
    for v in mu
        mmu = min(mmu, v)
        v <= 0 && (p.n_bad_mu += 1)
    end
    p.min_mu_rel = isempty(mu) ? NaN : mmu / s2

    p.cert = certify(c, r.X, eps)
    p.cert.theta_exact <= 1e-9 && (p.binding = binding_type(c, r.X))

    if do_cone
        er = K.expansive_cone(c, r.X, cone_opt)
        if er.ok
            p.cone_run = 1
            p.cone_feasible = (er.pass_feasible >= 0) ? 1 : 0
            p.cone_margin = er.margin_l2
            p.cone_dual = er.dual_bound_max
            p.cone_sigma_in = er.sigma_in_cone ? 1 : 0
        end
    end
    return p
end

Base.@kwdef mutable struct VariantStats
    n::Int = 0
    feasible::Int = 0
    theta_pos::Int = 0          # exact Theta_max > 0
    theta_pos_ref::Int = 0      # AND the referee agrees Theta_max > 0
    certified::Int = 0          # eps_max > 0
    mu_ok::Int = 0              # every 0+ corner margin strictly positive
    cone_run::Int = 0; cone_feasible::Int = 0
    theta::Vector{Float64} = Float64[]; eps_max::Vector{Float64} = Float64[]
    dist::Vector{Float64} = Float64[]; tnorm::Vector{Float64} = Float64[]; min_mu::Vector{Float64} = Float64[]
    binding::Dict{String,Int} = Dict{String,Int}()
    referee_agree::Int = 0; referee_n::Int = 0
    worst_gap::Float64 = 0.0
    sound::Int = 0; sound_n::Int = 0
end

function accumulate!(S::VariantStats, p::PointReport)
    S.n += 1
    if p.feasible
        S.feasible += 1
        push!(S.dist, p.dist)
        push!(S.tnorm, p.tnorm)
        isfinite(p.min_mu_rel) && push!(S.min_mu, p.min_mu_rel)
        p.n_bad_mu == 0 && (S.mu_ok += 1)
    end
    if p.cone_run != 0
        S.cone_run += 1
        S.cone_feasible += p.cone_feasible
    end
    push!(S.theta, p.cert.theta_exact)
    if p.cert.theta_exact > 1e-9
        S.theta_pos += 1
        p.cert.theta_bisect > 1e-9 && (S.theta_pos_ref += 1)
    end
    if p.cert.eps_max > 0
        S.certified += 1
        push!(S.eps_max, p.cert.eps_max)
    end
    S.binding[p.binding] = get(S.binding, p.binding, 0) + 1
    S.referee_n += 1
    gap = abs(p.cert.theta_exact - p.cert.theta_bisect)
    S.worst_gap = max(S.worst_gap, gap)
    gap <= 1e-5 && (S.referee_agree += 1)
    S.sound_n += 1
    p.cert.eps_max <= p.cert.theta_bisect + 1e-9 && (S.sound += 1)
end

function report_variant(o::IO, name::String, S::VariantStats)
    print(o, name, ": designs ", S.n, ", convex-feasible ", S.feasible, " (",
          fx(100.0 * S.feasible / max(1, S.n), 1), " %)\n")
    print(o, "  exact Theta_max > 0: ", S.theta_pos, " (",
          fx(100.0 * S.theta_pos / max(1, S.n), 1), " %), refereed ", S.theta_pos_ref, "\n")
    print(o, "  certified (eps_max > 0): ", S.certified, "\n")
    print(o, "  all 0+ corner margins > 0 at the returned X: ", S.mu_ok, " / ", S.feasible,
          " feasible\n")
    isempty(S.min_mu) || print(o, "  min mu / med^2 over feasible: median ", sci(median_of(S.min_mu)),
                               ", q10 ", sci(quantile_of(S.min_mu, 0.1)), ", q90 ",
                               sci(quantile_of(S.min_mu, 0.9)), "\n")
    isempty(S.dist) || print(o, "  |X - X_ini| per vertex, in median edges: median ", fx(median_of(S.dist)),
                             ", q90 ", fx(quantile_of(S.dist, 0.9)), "  (|t|/med median ",
                             fx(median_of(S.tnorm)), ")\n")
    isempty(S.eps_max) || print(o, "  eps_max distribution: median ", sci(median_of(S.eps_max)), ", max ",
                                sci(quantile_of(S.eps_max, 1.0)), "\n")
    S.cone_run != 0 && print(o, "  expansive-cone LP: run ", S.cone_run, ", feasible ", S.cone_feasible, "\n")
    print(o, "  referee: agree ", S.referee_agree, " / ", S.referee_n, ", worst gap ",
          sci(S.worst_gap), "; sound (eps_max <= bisection) ", S.sound, " / ", S.sound_n, "\n")
    print(o, "  binding at Theta_max = 0:")
    for k in sort!(collect(keys(S.binding)))   # std::map order
        print(o, " ", k, "=", S.binding[k])
    end
    print(o, "\n")
end

const VARIANT_KEYS = ["a", "a:sigma_mc", "a:sigma_def", "b", "b:sigma_mc", "b:sigma_def"]

function variant_label(key::String)
    base = key[1] == 'a' ? "(a) convexity only" : "(b) convexity + split"
    return occursin(':', key) ? base * " [" * split(key, ':')[2] * "]" : base
end

function write_verdict(o::IO, stat::Dict{String,VariantStats})
    A = get!(VariantStats, stat, "a")
    B = get!(VariantStats, stat, "b")
    best_theta = max(A.theta_pos_ref, B.theta_pos_ref)
    best_feas = max(A.feasible, B.feasible)
    N = max(1, A.n)
    pass = best_theta >= (10 * N + 99) ÷ 100
    hard_fail = (best_feas >= N ÷ 2) && (best_theta < (2 * N + 99) ÷ 100)
    print(o, "PASS bar: refereed Theta_max > 0 on >= 10 % of ", N, " designs (>= ",
          (10 * N + 99) ÷ 100, ")\n")
    print(o, "best variant: ", best_theta, " refereed positives, ", best_feas, " convex-feasible\n")
    print(o, "VERDICT: ", pass ? "PASS" : "FAIL", "\n")
    print(o, "HARD FAIL (feasible >= 50 % and Theta_max > 0 < 2 %): ",
          hard_fail ? "YES -- convexity is not the lever" : "no", "\n")
    print(o, "baselines: Eq. (6) alone 0/200 (K1a); K5 both sigma 0/200; K6 0/400\n")
end

# ---------------------------------------------------------------------------
# Shard aggregation. The run is sharded 12 ways, so the merged CSV and the merged
# summary are rebuilt from the shard CSVs by the SAME reporting code that a
# single-process run uses -- there is no second implementation of the statistics.

function aggregate(outdir::String, delta_rel::Float64, split_delta_rel::Float64,
                   n_random::Int, max_iter::Int, barrier_stages::Int, eps::Float64)
    files = sort!([joinpath(outdir, f) for f in readdir(outdir)
                   if startswith(f, "shard_") && endswith(f, ".csv")])
    if isempty(files)
        println(stderr, "aggregate: no shard CSVs in ", outdir)
        return 1
    end
    stat = Dict{String,VariantStats}()
    header = String[]
    col = Dict{String,Int}()
    merged = open(joinpath(outdir, "k9.csv"), "w")
    designs = 0; convex_at_x0 = 0; convex_at_x0_theta_pos = 0
    split_only_feasible = 0
    tot_corner = 0; tot_nonconvex0 = 0
    seen_graph = Set{String}()

    for f in files
        lines = readlines(f)
        isempty(lines) && continue
        if isempty(header)
            header = split(lines[1], ',')
            for (i, h) in enumerate(header)
                col[h] = i
            end
            print(merged, lines[1], "\n")
        end
        for line in lines[2:end]
            isempty(line) && continue
            print(merged, line, "\n")
            r = split(line, ',')
            if length(r) != length(header)
                println(stderr, "aggregate: bad row width in ", f)
                return 1
            end
            D(k) = parse(Float64, r[col[k]])
            I(k) = parse(Int, r[col[k]])
            S(k) = String(r[col[k]])
            designs += 1
            push!(seen_graph, S("id"))
            tot_corner += I("n_corner")
            tot_nonconvex0 += I("nonconvex0")
            I("split_only_feasible") != 0 && (split_only_feasible += 1)

            p0 = PointReport()
            p0.feasible = (I("nonconvex0") == 0)
            p0.cert.eps_max = D("x0_eps_max")
            p0.cert.theta_exact = D("x0_theta_exact")
            p0.cert.theta_bisect = D("x0_theta_bisect")
            p0.binding = S("x0_binding")
            p0.min_mu_rel = D("min_mu0")
            p0.n_bad_mu = I("n_bad_mu0")
            accumulate!(get!(VariantStats, stat, "x0"), p0)
            if p0.feasible
                convex_at_x0 += 1
                p0.cert.theta_exact > 1e-9 && (convex_at_x0_theta_pos += 1)
            end

            fname = S("sigma")
            for v in 0:1
                pre = v == 1 ? "b_" : "a_"
                p = PointReport()
                p.feasible = I(pre * "feasible") != 0
                p.min_cross_rel = D(pre * "min_cross")
                p.n_bad = I(pre * "nbad")
                p.dist = D(pre * "dist")
                p.tnorm = D(pre * "tnorm")
                p.min_mu_rel = D(pre * "min_mu")
                p.n_bad_mu = I(pre * "n_bad_mu")
                p.cert.eps_max = D(pre * "eps_max")
                p.cert.theta_exact = D(pre * "theta_exact")
                p.cert.theta_bisect = D(pre * "theta_bisect")
                p.binding = S(pre * "binding")
                p.cone_run = I(pre * "cone_run")
                p.cone_feasible = I(pre * "cone_feasible")
                tag = v == 1 ? "b" : "a"
                accumulate!(get!(VariantStats, stat, tag), p)
                accumulate!(get!(VariantStats, stat, tag * ":" * fname), p)
            end
        end
    end
    close(merged)

    o = open(joinpath(outdir, "summary.txt"), "w")
    print(o, "K9 -- convexity-constrained embedding in the Tutte auxetic null space\n")
    print(o, "population: ", length(seen_graph), " graphs (K6's), ", designs,
          " designs (2 sigma), merged from ", length(files), " shards\n")
    print(o, "delta = ", cpp_g(delta_rel), " * med^2, split delta' = ", cpp_g(split_delta_rel),
          " * med^2, starts ", n_random, " + t=0, L-BFGS iters ", max_iter,
          ", barrier stages ", barrier_stages, ", eps = ", cpp_g(eps), "\n")
    print(o, "non-convex corners at X0: ", tot_nonconvex0, " / ", tot_corner, " (",
          fx(100.0 * tot_nonconvex0 / max(1, tot_corner), 2), " %)\n")
    print(o, "split-only repair (K6's secondary objective, the (b) fallback start) feasible on ",
          split_only_feasible, " / ", designs, " designs\n")
    print(o, "designs whose X0 is ALREADY fully convex: ", convex_at_x0, " / ", designs,
          ", of which exact Theta_max > 0: ", convex_at_x0_theta_pos, "\n\n")
    if haskey(stat, "x0")
        report_variant(o, "X0 baseline (no solve)", stat["x0"])
        print(o, "\n")
    end
    for key in VARIANT_KEYS
        haskey(stat, key) || continue
        report_variant(o, variant_label(key), stat[key])
        print(o, "\n")
    end
    write_verdict(o, stat)
    close(o)
    println("aggregated ", designs, " designs from ", length(files), " shards into ",
            outdir, "/k9.csv and ", outdir, "/summary.txt")
    return 0
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

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n = 200; maxf = 800; max_iter = 600; n_random = 3
    eps = 0.3; delta_rel = 1e-3; split_delta_rel = 1e-3
    barrier_stages = 6
    outdir = joinpath(REPO, "results", "kill", "k9_julia")
    sigmadir = ""
    cache = joinpath(REPO, "results", "kill", "k6", "cache")
    shard = 0; nshards = 1; cone_maxf = 800
    do_cone = true; do_aggregate = false
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
        elseif a == "--delta" && i < length(args); delta_rel = arg_f(args[i+1]); i += 2
        elseif a == "--split-delta" && i < length(args); split_delta_rel = arg_f(args[i+1]); i += 2
        elseif a == "--iters" && i < length(args); max_iter = arg_i(args[i+1]); i += 2
        elseif a == "--starts" && i < length(args); n_random = arg_i(args[i+1]); i += 2
        elseif a == "--barrier" && i < length(args); barrier_stages = arg_i(args[i+1]); i += 2
        elseif a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshards" && i < length(args); nshards = arg_i(args[i+1]); i += 2
        elseif a == "--cone-maxf" && i < length(args); cone_maxf = arg_i(args[i+1]); i += 2
        elseif a == "--no-cone"; do_cone = false; i += 1
        elseif a == "--aggregate"; do_aggregate = true; i += 1
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    do_aggregate && return aggregate(outdir, delta_rel, split_delta_rel, n_random, max_iter,
                                     barrier_stages, eps)
    mkpath(outdir)
    csv_path = nshards > 1 ? joinpath(outdir, "shard_" * string(shard) * ".csv") :
                             joinpath(outdir, "k9.csv")
    csv = open(csv_path, "w")
    print(csv, "id,kind,sigma,N,F,med_edge,n_split,dim_null,m,n_corner,",
          "nonconvex0,min_cross0,n_bad_mu0,min_mu0,",
          "x0_eps_max,x0_theta_exact,x0_theta_bisect,x0_binding,",
          "a_feasible,a_min_cross,a_nbad,a_nonconvex,a_ninv,a_dist,a_dist_x0,a_tnorm,",
          "a_iters,a_barrier,a_min_mu,a_n_bad_mu,a_eps_max,a_theta_exact,a_theta_bisect,",
          "a_binding,a_cone_run,a_cone_feasible,a_cone_margin,a_cone_dual,a_cone_sigma,",
          "b_feasible,b_min_cross,b_min_q,b_nbad,b_nbadq,b_ninv,b_dist,b_tnorm,b_barrier,",
          "b_min_mu,b_n_bad_mu,b_eps_max,b_theta_exact,b_theta_bisect,b_binding,",
          "b_cone_run,b_cone_feasible,b_cone_margin,b_cone_dual,b_cone_sigma,",
          "split_only_feasible,split_only_min_q,secs\n")
    wall = Timer()

    stat = Dict{String,VariantStats}()  # "a" / "b", and per-sigma "a:sigma_mc" etc
    graphs = 0; missing_sigma = 0; designs = 0
    tot_corner = 0; tot_nonconvex0 = 0
    convex_at_x0 = 0; convex_at_x0_theta_pos = 0; split_only_feasible = 0

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

            # ---- baseline at X0 -------------------------------------------------
            cr0 = K.corner_crosses(m, X0)
            nonconvex0 = 0
            min_cross0 = Inf
            for v in cr0
                v <= 0 && (nonconvex0 += 1)
                min_cross0 = min(min_cross0, v)
            end
            tot_corner += length(cr0)
            tot_nonconvex0 += nonconvex0
            mu0 = K.zero_plus_corner_margin(c, X0)
            n_bad_mu0 = 0
            min_mu0 = Inf
            for v in mu0
                v <= 0 && (n_bad_mu0 += 1)
                min_mu0 = min(min_mu0, v)
            end

            # The X0 baseline in the SAME columns: K6 reports 0/400 there, and the designs
            # whose X0 is ALREADY fully convex are the cleanest test of the hypothesis --
            # convexity costs them nothing, so if it were the lever they would deploy.
            x0c = certify(c, X0, eps)
            x0_bind = "n/a"
            x0c.theta_exact <= 1e-9 && (x0_bind = binding_type(c, X0))
            let S = get!(VariantStats, stat, "x0")
                p0 = PointReport()
                p0.feasible = (nonconvex0 == 0)
                p0.cert = x0c
                p0.binding = x0_bind
                p0.min_mu_rel = isempty(mu0) ? NaN : min_mu0 / s2
                p0.n_bad_mu = n_bad_mu0
                accumulate!(S, p0)
                if nonconvex0 == 0
                    convex_at_x0 += 1
                    x0c.theta_exact > 1e-9 && (convex_at_x0_theta_pos += 1)
                end
            end

            cone_opt = K.ExpansiveConeOptions()
            cone_here = do_cone && K.n_faces(m) <= cone_maxf

            # ---- variant (a): convexity only -------------------------------------
            oa = K.ConvexEmbedOptions()
            oa.delta_rel = delta_rel
            oa.n_random = n_random
            oa.max_iter = max_iter
            oa.barrier_stages = barrier_stages
            oa.seed = (UInt32(9000) + UInt32(7) * UInt32(id) + UInt32(which)) % UInt32
            ra = K.convex_embed(c, X0, sh.Phi, X_ini, med, oa)
            pa = evaluate(c, m, ra, med, eps, cone_here && ra.feasible, cone_opt)

            # ---- variant (b): convexity + the split-cut signs ---------------------
            # Both constraint families are relaxations of (b), so (b) is started from THREE
            # points and the best answer kept: the origin, (a)'s convex minimiser, and the
            # split-only minimiser of K6's secondary objective (zero_plus_repair with
            # w_corner = w_prox = 0). K6 measured the split-only problem to be solvable on
            # most designs, so if (b) still fails it is not for want of a good start.
            ro = K.ZeroPlusRepairOptions()
            ro.n_random = n_random
            ro.max_iter = max_iter
            ro.seed = (UInt32(9000) + UInt32(7) * UInt32(id) + UInt32(which)) % UInt32
            sp = K.zero_plus_repair(c, X0, sh.Phi, med, ro)

            ob = deepcopy(oa)
            ob.split_delta_rel = split_delta_rel
            rb = K.convex_embed(c, X0, sh.Phi, X_ini, med, ob)
            if !rb.feasible && ra.feasible
                o2 = deepcopy(ob)
                o2.t_init = ra.t
                r2 = K.convex_embed(c, X0, sh.Phi, X_ini, med, o2)
                (r2.feasible || r2.n_bad + r2.n_bad_q < rb.n_bad + rb.n_bad_q) && (rb = r2)
            end
            if !rb.feasible && sp.feasible
                o3 = deepcopy(ob)
                o3.t_init = sp.t
                r3 = K.convex_embed(c, X0, sh.Phi, X_ini, med, o3)
                (r3.feasible || r3.n_bad + r3.n_bad_q < rb.n_bad + rb.n_bad_q) && (rb = r3)
            end
            pb = evaluate(c, m, rb, med, eps, cone_here && rb.feasible, cone_opt)

            sp.feasible && (split_only_feasible += 1)
            accumulate!(get!(VariantStats, stat, "a"), pa)
            accumulate!(get!(VariantStats, stat, "b"), pb)
            accumulate!(get!(VariantStats, stat, "a:" * fname), pa)
            accumulate!(get!(VariantStats, stat, "b:" * fname), pb)

            secs = s(t)
            g_ = cpp_g
            print(csv, id, ",", row.kind, ",", fname, ",", K.n_vertices(m), ",",
                  K.n_faces(m), ",", g_(med), ",", K.n_split(c), ",", sh.k, ",",
                  2 * sh.k, ",", length(cr0), ",", nonconvex0, ",",
                  g_(min_cross0 / s2), ",", n_bad_mu0, ",", g_(min_mu0 / s2), ",",
                  g_(x0c.eps_max), ",", g_(x0c.theta_exact), ",", g_(x0c.theta_bisect), ",",
                  x0_bind, ",",
                  (pa.feasible ? 1 : 0), ",", g_(pa.min_cross_rel), ",", pa.n_bad, ",",
                  pa.n_nonconvex, ",", pa.n_inv, ",", g_(pa.dist), ",", g_(ra.dist_ini_x0),
                  ",", g_(pa.tnorm), ",", pa.iters, ",", pa.barrier, ",",
                  g_(pa.min_mu_rel), ",", pa.n_bad_mu, ",", g_(pa.cert.eps_max), ",",
                  g_(pa.cert.theta_exact), ",", g_(pa.cert.theta_bisect), ",", pa.binding, ",",
                  pa.cone_run, ",", pa.cone_feasible, ",", g_(pa.cone_margin), ",",
                  g_(pa.cone_dual), ",", pa.cone_sigma_in, ",",
                  (pb.feasible ? 1 : 0), ",", g_(pb.min_cross_rel), ",", g_(pb.min_q_rel), ",",
                  pb.n_bad, ",", pb.n_bad_q, ",", pb.n_inv, ",", g_(pb.dist), ",",
                  g_(pb.tnorm), ",", pb.barrier, ",", g_(pb.min_mu_rel), ",", pb.n_bad_mu,
                  ",", g_(pb.cert.eps_max), ",", g_(pb.cert.theta_exact), ",",
                  g_(pb.cert.theta_bisect), ",", pb.binding, ",", pb.cone_run, ",",
                  pb.cone_feasible, ",", g_(pb.cone_margin), ",", g_(pb.cone_dual), ",",
                  pb.cone_sigma_in, ",", (sp.feasible ? 1 : 0), ",", g_(sp.min_q / s2),
                  ",", g_(secs), "\n")
            flush(csv)
            println("id ", id, " ", fname, " N=", K.n_vertices(m),
                    " k=", sh.k, " nonconvex0=", nonconvex0, "/", length(cr0),
                    " a_feas=", cpp_b(pa.feasible), " a_theta=", g_(pa.cert.theta_exact),
                    " b_feas=", cpp_b(pb.feasible), " b_theta=", g_(pb.cert.theta_exact),
                    " (", fx(secs, 1), " s)")
        end
    end
    close(csv)

    sum_path = nshards > 1 ? joinpath(outdir, "shard_" * string(shard) * ".txt") :
                             joinpath(outdir, "summary.txt")
    o = open(sum_path, "w")
    print(o, "K9 -- convexity-constrained embedding in the Tutte auxetic null space\n")
    print(o, "population: ", graphs, " graphs (K6's), ", designs,
          " designs (2 sigma), missing sigma_def ", missing_sigma, "\n")
    print(o, "delta = ", cpp_g(delta_rel), " * med^2, split delta' = ", cpp_g(split_delta_rel),
          " * med^2, starts ", n_random, " + t=0, L-BFGS iters ", max_iter,
          ", barrier stages ", barrier_stages, ", eps = ", cpp_g(eps), "\n")
    print(o, "non-convex corners at X0: ", tot_nonconvex0, " / ", tot_corner, " (",
          fx(100.0 * tot_nonconvex0 / max(1, tot_corner), 2), " %)\n\n")
    print(o, "split-only repair (K6's secondary objective, the (b) fallback start) feasible on ",
          split_only_feasible, " / ", designs, " designs\n")
    print(o, "designs whose X0 is ALREADY fully convex: ", convex_at_x0, " / ", designs,
          ", of which exact Theta_max > 0: ", convex_at_x0_theta_pos, "\n\n")
    if haskey(stat, "x0")
        report_variant(o, "X0 baseline (no solve)", stat["x0"])
        print(o, "\n")
    end
    for key in VARIANT_KEYS
        haskey(stat, key) || continue
        report_variant(o, variant_label(key), stat[key])
        print(o, "\n")
    end
    write_verdict(o, stat)
    print(o, "wall ", fx(s(wall), 1), " s\n")
    close(o)
    println("wrote ", csv_path, " and ", sum_path)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
