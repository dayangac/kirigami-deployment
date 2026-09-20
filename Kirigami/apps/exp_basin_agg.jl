# WP2b aggregator -- merges the exp_basin shards and derives every number BASIN.md quotes.
#
# Kept in its own file so that re-running the analysis never touches the exp_basin
# driver while its shards are still running. Every statistic is a function of
# results/yield/basin.csv on disk (and, for the last section, of results/yield/features.csv),
# so BASIN.md can be regenerated from the archived CSVs alone.
#
#   --mode merge      concatenate basin_shard_*.csv into basin.csv
#   --mode stats      read basin.csv and write basin_stats.txt (implies merge unless
#                     --no-merge)
#
# Definitions used throughout, all matching K9c's own conventions:
#   deploys        exact Theta_max > 1e-9 (the threshold K9c's summary uses)
#   flip           an originally failing design for which at least one of the eight NEW
#                  seeds (k = 1..8) deploys; k = 0 is the archived K9c answer and is zero
#                  by construction on the failure set
#
#   julia --project=Kirigami Kirigami/apps/exp_basin_agg.jl [--out DIR] [--mode stats]
#         [--nshards 3] [--no-merge]
#
# Reads and writes under results/yield/ by default.
include(joinpath(@__DIR__, "common_app.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# `std::getline(ss, x, ',')` over a line: a trailing empty field is dropped.
function split_csv(line::AbstractString)
    f = String.(split(line, ','))
    (!isempty(f) && isempty(f[end]) && endswith(line, ',')) && pop!(f)
    return f
end

# `o << std::setprecision(10) << double`.
p10(v::Real) = isfinite(v) ? @sprintf("%.10g", Float64(v)) : _fmt_nonfinite(Float64(v))

# Type-7 (linear interpolation) quantile, the convention numpy.percentile uses.
function quantile7(v_in::Vector{Float64}, p::Float64)
    isempty(v_in) && return NaN
    v = sort(v_in)
    length(v) == 1 && return v[1]
    h = p * (length(v) - 1.0)
    lo = Int(floor(h))
    hi = min(length(v) - 1, lo + 1)
    return v[lo + 1] + (h - lo) * (v[hi + 1] - v[lo + 1])
end

# Rank-based AUC (Mann-Whitney), ties averaged. > 0.5 means larger x is more often a
# positive. Returns NaN when either class is empty.
function auc(x::Vector{Float64}, y::Vector{Int})
    idx = sortperm(x; alg = MergeSort)   # std::sort with the `<` comparator; ties are grouped below
    rank = zeros(length(x))
    i = 1
    while i <= length(idx)
        j = i
        while j + 1 <= length(idx) && x[idx[j + 1]] == x[idx[i]]
            j += 1
        end
        r = 0.5 * (i + j)
        for t in i:j
            rank[idx[t]] = r
        end
        i = j + 1
    end
    sum_pos = 0.0
    n_pos = 0; n_neg = 0
    for k in eachindex(x)
        if y[k] != 0
            n_pos += 1
            sum_pos += rank[k]
        else
            n_neg += 1
        end
    end
    (n_pos == 0 || n_neg == 0) && return NaN
    return (sum_pos - 0.5 * n_pos * (n_pos + 1)) / (Float64(n_pos) * Float64(n_neg))
end

# Wilson 95 % interval for a binomial proportion.
function wilson(k::Int, n::Int)
    n == 0 && return (NaN, NaN)
    z = 1.959963984540054
    p = k / n
    d = 1.0 + z * z / n
    c = p + z * z / (2.0 * n)
    s = z * sqrt(p * (1 - p) / n + z * z / (4.0 * n * n))
    return ((c - s) / d, (c + s) / d)
end

Base.@kwdef mutable struct Seed
    k::Int = -1
    theta::Float64 = 0.0; eps::Float64 = 0.0; margin::Float64 = 0.0; secs::Float64 = 0.0
    ref9::Float64 = -1.0; ref12::Float64 = -1.0
    prov::String = ""; binding::String = ""; status::String = ""
    feasible::Int = 0
end

Base.@kwdef mutable struct DesignRows
    id::Int = -1; which::Int = 0
    kind::String = ""; sigma::String = ""; cls::String = ""
    orig_binding::String = ""
    seeds::Vector{Seed} = Seed[]   # k = 0..8, in file order
end

key_of(id::Int, sigma::AbstractString) = string(id) * "|" * sigma

# std::map iteration order: keys sorted.
sorted_keys(d::AbstractDict) = sort!(collect(keys(d)))

function main(args::Vector{String})
    outdir = joinpath(REPO, "results", "yield")
    mode = "stats"
    nshards = 3
    do_merge = true
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--mode" && i < length(args); mode = args[i+1]; i += 2
        elseif a == "--nshards" && i < length(args); nshards = arg_i(args[i+1]); i += 2
        elseif a == "--no-merge"; do_merge = false; i += 1
        else
            println(stderr, "unknown argument ", a)
            return 2
        end
    end

    merged = joinpath(outdir, "basin.csv")
    if do_merge
        out = open(merged, "w")
        header = false
        n = 0
        for sh in 0:nshards-1
            p = joinpath(outdir, "basin_shard_" * string(sh) * ".csv")
            if !isfile(p)
                println(stderr, "missing ", p)
                continue
            end
            first = true
            for line in readlines(p)
                if first
                    first = false
                    if !header
                        print(out, line, "\n")
                        header = true
                    end
                    continue
                end
                if !isempty(line)
                    print(out, line, "\n")
                    n += 1
                end
            end
        end
        close(out)
        println("merged ", n, " rows into ", merged)
    end
    mode == "merge" && return 0

    # ---- read the merged CSV back ------------------------------------------------------
    D = Dict{String,DesignRows}()
    order = String[]
    let
        if !isfile(merged)
            println(stderr, "cannot open ", merged)
            return 1
        end
        lines = readlines(merged)
        head = split_csv(isempty(lines) ? "" : lines[1])
        function col(nm)
            j = findfirst(==(nm), head)
            if j === nothing
                println(stderr, "missing column ", nm)
                exit(1)
            end
            return j
        end
        c_id = col("id"); c_kind = col("kind"); c_sig = col("sigma")
        c_which = col("which"); c_cls = col("orig_class"); c_k = col("k")
        c_st = col("status"); c_feas = col("feasible"); c_marg = col("margin")
        c_prov = col("provenance"); c_th = col("theta_exact"); c_eps = col("eps_max")
        c_bind = col("binding"); c_r12 = col("theta_ref12"); c_r9 = col("theta_ref9")
        c_secs = col("secs")
        for line in lines[2:end]
            isempty(line) && continue
            v = split_csv(line)
            length(v) != length(head) && continue
            id = parse(Int, v[c_id])
            key = key_of(id, v[c_sig])
            if !haskey(D, key)
                d = DesignRows()
                d.id = id
                d.kind = v[c_kind]
                d.sigma = v[c_sig]
                d.which = parse(Int, v[c_which])
                d.cls = v[c_cls]
                D[key] = d
                push!(order, key)
            end
            sd = Seed()
            sd.k = parse(Int, v[c_k])
            sd.status = v[c_st]
            sd.feasible = parse(Int, v[c_feas])
            sd.margin = parse(Float64, v[c_marg])
            sd.prov = v[c_prov]
            sd.theta = parse(Float64, v[c_th])
            sd.eps = parse(Float64, v[c_eps])
            sd.binding = v[c_bind]
            sd.ref12 = parse(Float64, v[c_r12])
            sd.ref9 = parse(Float64, v[c_r9])
            sd.secs = parse(Float64, v[c_secs])
            sd.k == 0 && (D[key].orig_binding = sd.binding)
            push!(D[key].seeds, sd)
        end
    end

    o = open(joinpath(outdir, "basin_stats.txt"), "w")
    print(o, "WP2b -- solver basin vs geometry, derived from ", merged, "\n")
    print(o, "deploy threshold: exact Theta_max > 1e-9\n\n")

    kEps = 1e-9
    n_fail = 0; n_succ = 0; incomplete = 0
    for key in order
        d = D[key]
        d.cls == "fail" ? (n_fail += 1) : (n_succ += 1)
        length(d.seeds) != 9 && (incomplete += 1)
    end
    print(o, "designs: ", length(D), " (", n_fail, " originally failing, ", n_succ,
          " sampled successes); designs without all 9 rows: ", incomplete, "\n\n")

    # ---- 1. flip rate ------------------------------------------------------------------
    flipped = 0; fail_complete = 0
    by_kind = Dict{String,Vector{Int}}(); by_binding = Dict{String,Vector{Int}}()   # [flipped, total]
    by_sigma = Dict{String,Vector{Int}}()
    ndeploy_hist = zeros(Int, 9)
    flip_names = String[]
    flip_theta = Float64[]
    flip_eps01 = 0; flip_ref9 = 0
    bump!(M, k, fl) = (v = get!(M, k, [0, 0]); v[1] += fl ? 1 : 0; v[2] += 1)
    for key in order
        d = D[key]
        d.cls != "fail" && continue
        nd = 0
        best = 0.0; best_eps = 0.0; best_ref9 = -1.0
        for sd in d.seeds
            sd.k == 0 && continue
            if sd.theta > kEps
                nd += 1
                if sd.theta > best
                    best = sd.theta; best_eps = sd.eps; best_ref9 = sd.ref9
                end
            end
        end
        fail_complete += 1
        nd < 9 && (ndeploy_hist[nd + 1] += 1)
        fl = nd >= 1
        if fl
            flipped += 1
            # std::to_string(double): fixed, 6 decimals
            push!(flip_names, d.kind * "_" * string(d.id) * "_" * d.sigma * " (n=" * string(nd) *
                              ", theta=" * fx(best, 6) * ")")
            push!(flip_theta, best)
            best_eps >= 0.1 && (flip_eps01 += 1)
            best_ref9 > kEps && (flip_ref9 += 1)
        end
        bump!(by_kind, d.kind, fl)
        bump!(by_binding, d.orig_binding, fl)
        bump!(by_sigma, d.sigma, fl)
    end
    lo, hi = wilson(flipped, fail_complete)
    print(o, "== 1. FLIP RATE OF THE FAILURES (>= 1 of 8 new seeds deploying)\n")
    print(o, "overall: ", flipped, " / ", fail_complete, " = ",
          p10(fail_complete != 0 ? 100.0 * flipped / fail_complete : 0.0), " %, Wilson 95 % [",
          p10(100 * lo), ", ", p10(100 * hi), "] %\n")
    print(o, "of the flipped, eps_max >= 0.1 rad at the best seed: ", flip_eps01,
          "; confirmed by the 1e-9 bisection referee: ", flip_ref9, "\n")
    function dump(title, M)
        print(o, title, "\n")
        for k in sorted_keys(M)
            a_, b_ = wilson(M[k][1], M[k][2])
            print(o, "  ", k, ": ", M[k][1], " / ", M[k][2], " = ",
                  p10(M[k][2] != 0 ? 100.0 * M[k][1] / M[k][2] : 0.0), " %  [", p10(100 * a_), ", ",
                  p10(100 * b_), "]\n")
        end
    end
    dump("by family:", by_kind)
    dump("by original binding:", by_binding)
    dump("by orientation rule:", by_sigma)
    print(o, "distribution of the number of deploying seeds (out of 8), failures only:\n")
    for j in 0:8
        print(o, "  n=", j, ": ", ndeploy_hist[j + 1], "\n")
    end
    if !isempty(flip_theta)
        print(o, "best Theta_max over the flipped designs: min ", p10(quantile7(flip_theta, 0.0)),
              ", median ", p10(quantile7(flip_theta, 0.5)), ", max ", p10(quantile7(flip_theta, 1.0)), "\n")
    end
    print(o, "flipped designs:\n")
    for sname in flip_names
        print(o, "  ", sname, "\n")
    end
    print(o, "\n")

    # ---- 2. seed-to-seed spread on the successes ---------------------------------------
    print(o, "== 2. SEED-TO-SEED SPREAD ON THE 30 SAMPLED SUCCESSES (9 seeds each)\n")
    print(o, "id,sigma,n_seeds,min,q25,median,q75,max,iqr,rel_iqr,theta_k0,best_k,orig_is_best,",
          "n_zero\n")
    orig_best = 0; succ_n = 0; succ_any_zero = 0; succ_all_deploy = 0
    rel_iqr = Float64[]; rel_gain = Float64[]
    succ_eps01_k0 = 0; succ_eps01_best = 0
    for key in order
        d = D[key]
        d.cls != "success" && continue
        succ_n += 1
        th = Float64[]
        th0 = 0.0; best = -1.0; best_eps = 0.0; eps0 = 0.0
        best_k = -1; n_zero = 0
        for sd in d.seeds
            push!(th, sd.theta)
            sd.theta <= kEps && (n_zero += 1)
            if sd.k == 0
                th0 = sd.theta; eps0 = sd.eps
            end
            if sd.theta > best
                best = sd.theta; best_k = sd.k; best_eps = sd.eps
            end
        end
        n_zero > 0 ? (succ_any_zero += 1) : (succ_all_deploy += 1)
        q25 = quantile7(th, 0.25); q50 = quantile7(th, 0.5); q75 = quantile7(th, 0.75)
        mn = quantile7(th, 0.0); mx = quantile7(th, 1.0)
        ob = th0 >= best - 1e-12
        ob && (orig_best += 1)
        q50 > 0 && push!(rel_iqr, (q75 - q25) / q50)
        th0 > kEps && push!(rel_gain, best / th0)
        eps0 >= 0.1 && (succ_eps01_k0 += 1)
        best_eps >= 0.1 && (succ_eps01_best += 1)
        print(o, d.id, ",", d.sigma, ",", length(th), ",", p10(mn), ",", p10(q25), ",", p10(q50), ",",
              p10(q75), ",", p10(mx), ",", p10(q75 - q25), ",",
              p10(q50 > 0 ? (q75 - q25) / q50 : NaN), ",", p10(th0), ",", best_k, ",",
              ob ? 1 : 0, ",", n_zero, "\n")
    end
    print(o, "successes: ", succ_n, "; original seed was the best of the 9 on ", orig_best, " (",
          p10(succ_n != 0 ? 100.0 * orig_best / succ_n : 0.0), " %)\n")
    print(o, "successes with at least one NON-deploying seed: ", succ_any_zero, "; all 9 deploy: ",
          succ_all_deploy, "\n")
    isempty(rel_iqr) || print(o, "relative IQR (IQR / median Theta_max): median ", p10(quantile7(rel_iqr, 0.5)),
                              ", q90 ", p10(quantile7(rel_iqr, 0.9)), ", max ", p10(quantile7(rel_iqr, 1.0)), "\n")
    isempty(rel_gain) || print(o, "best / original Theta_max: median ", p10(quantile7(rel_gain, 0.5)), ", q90 ",
                               p10(quantile7(rel_gain, 0.9)), ", max ", p10(quantile7(rel_gain, 1.0)), "\n")
    print(o, "sampled successes with eps_max >= 0.1: k = 0: ", succ_eps01_k0, ", best of 9: ",
          succ_eps01_best, "\n\n")

    # ---- 3. best-of-9 yield on the 400 --------------------------------------------------
    print(o, "== 3. BEST-OF-9 YIELD ON THE 400 K9c DESIGNS\n")
    base_pos = 307   # K9c single-seed positives (results/experiments/k9c/k9c.csv)
    total = 400
    print(o, "single-seed (archived K9c): ", base_pos, " / ", total, " = ", p10(100.0 * base_pos / total), " %\n")
    best9 = base_pos + flipped
    lo, hi = wilson(best9, total)
    print(o, "best-of-9 (EXACT on the deploy count: every one of the 93 failures was rerun, ",
          "and a success stays a success because k = 0 is one of the nine): ", best9, " / ", total,
          " = ", p10(100.0 * best9 / total), " %, Wilson 95 % [", p10(100 * lo), ", ", p10(100 * hi), "] %\n")
    eps01_succ_rate_k0 = succ_n != 0 ? succ_eps01_k0 / succ_n : 0.0
    eps01_succ_rate_b9 = succ_n != 0 ? succ_eps01_best / succ_n : 0.0
    print(o, "eps_max >= 0.1 rad: EXTRAPOLATED from the ", succ_n, " sampled successes only ",
          "(the other ", base_pos - succ_n, " successes were not rerun).\n")
    print(o, "  sampled-success rate at k = 0: ", p10(eps01_succ_rate_k0), ", best of 9: ",
          p10(eps01_succ_rate_b9), "\n")
    print(o, "  implied count over the 307 successes: ", p10(eps01_succ_rate_k0 * base_pos), " -> ",
          p10(eps01_succ_rate_b9 * base_pos), ", plus ", flip_eps01, " from the flipped failures\n")
    print(o, "  implied 400-design rate: ", p10(100.0 * (eps01_succ_rate_b9 * base_pos + flip_eps01) / total),
          " % (best of 9); treat as an estimate, not a measurement.\n\n")

    # ---- 4. verdict ---------------------------------------------------------------------
    flip_pct = fail_complete != 0 ? 100.0 * flipped / fail_complete : 0.0
    print(o, "== 4. PRE-REGISTERED RULE\n")
    print(o, "flip rate ", p10(flip_pct), " % ", flip_pct >= 20.0 ? ">=" : "<", " 20 % -> ",
          flip_pct >= 20.0 ? "SEARCH FLOOR: the paper must report best-of-9 alongside single-seed." :
                             "PREDOMINANTLY GEOMETRIC: the failure set is not a solver basin.",
          "\n\n")

    # ---- 5. what distinguishes flipped from never-deploying (only if the rule fires) -----
    if flip_pct >= 20.0
        print(o, "== 5. FLIPPED vs NEVER-DEPLOYING, univariate AUC on WP2's features.csv\n")
        print(o, "AUC > 0.5 means a LARGER feature value makes a flip more likely. No fitting, no\n",
              "multivariate model, no held-out split: this is a description of the 93 failures.\n")
        label = Dict{String,Int}()   # key -> 1 flipped, 0 never
        for key in order
            d = D[key]
            d.cls != "fail" && continue
            nd = count(sd -> sd.k != 0 && sd.theta > kEps, d.seeds)
            label[key_of(d.id, d.sigma)] = nd >= 1 ? 1 : 0
        end
        fpath = joinpath(outdir, "features.csv")
        if !isfile(fpath)
            print(o, "features.csv not found; section skipped.\n")
        else
            flines = readlines(fpath)
            head = split_csv(flines[1])
            c_id = 0; c_sig = 0; c_stop = length(head)
            for (j, h) in enumerate(head)
                h == "id" && (c_id = j)
                h == "sigma" && (c_sig = j)
                h == "theta_exact" && (c_stop = j - 1)
            end
            cols = [Float64[] for _ in head]
            y = Int[]
            matched = 0
            for line in flines[2:end]
                isempty(line) && continue
                v = split_csv(line)
                length(v) != length(head) && continue
                key = key_of(parse(Int, v[c_id]), v[c_sig])
                haskey(label, key) || continue   # not one of the 93 failures
                matched += 1
                push!(y, label[key])
                for j in 1:c_stop
                    if j == c_id || j == c_sig || head[j] == "kind"
                        push!(cols[j], 0.0)
                        continue
                    end
                    x = tryparse(Float64, v[j])
                    push!(cols[j], x === nothing ? NaN : x)
                end
            end
            print(o, "joined ", matched, " of the ", length(label), " failing designs on (id, sigma)\n")
            ranked = Tuple{Float64,String}[]
            for j in 1:c_stop
                (j == c_id || j == c_sig || head[j] == "kind") && continue
                (any(isnan, cols[j]) || length(cols[j]) != length(y)) && continue
                a_ = auc(cols[j], y)
                isnan(a_) || push!(ranked, (abs(a_ - 0.5), head[j]))
            end
            sort!(ranked; rev = true)   # std::sort(rbegin, rend): descending on (|AUC-0.5|, name)
            print(o, "feature,AUC (sorted by |AUC - 0.5|)\n")
            for r in ranked
                j = findlast(==(r[2]), head)
                print(o, "  ", r[2], ",", p10(auc(cols[j], y)), "\n")
            end
        end
    end
    close(o)
    println("wrote ", outdir, "/basin_stats.txt")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
