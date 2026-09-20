# WP7a -- the PAPER-REGIME population.
#
# REPORT.md's objection 2: every worked example in the 2026 paper has F in [4, 97], while
# the population K1a/K5/K6/K9/K9b/K9c/Native200 all share has median |F| = 332.5 and
# range [101, 793]. A referee is entitled to ask whether the paper's own Eq. (6)
# projection fails on random graphs AT THE PAPER'S OWN SCALE, and whether either the
# authors' full pipeline or ours works there. This driver answers that with the same
# instruments, on a population built the same way but with the face window moved down:
#
#   make_graph(id, 20, 100, 1400) for id = 2000 .. 2099  (three families by id % 3)
#
# Two orientation rules per graph -- sigma_mc = Eq. (1) max-cut (what make_graph already
# puts in mesh.sigma) and sigma_def = K5's defect-minimising search (method::
# orientation_defect at exp_k5_orientation.jl's cap of 20 |F| attempts, seed 7000 + id) -- so 200
# designs, and FOUR arms per design:
#
#   baseline    method::design_baseline           -- Eq. (6) alone, t = 0
#   k9c         method::design_range_max          -- the deliverable method, run defaults
#   natour      the authors' CLI on OUR sigma     -- their solve + their Eq. (9)
#   natcol      the authors' CLI on THEIR colour  -- one cell per GRAPH, shared by the
#                                                    graph's two sigma rows
#
# Every arm is scored by exactly one instrument: method::characterize at the arm's point,
# giving the exact T4.2'' Theta_max, the binding classification, and eps_max (the largest
# eps the repaired POS /\ NOOVERLAP(eps/2) /\ NOROOT(eps) certificate admits), plus the
# independent bisection referee at shrink 1e-12 (characterize's own) and at 1e-9 (the
# shrink K9c's headline "refereed" counts use). The native arms are RUN exactly as
# apps/exp_native200.jl runs them -- the shared code path is native_common.jl, which
# that driver also uses -- and Native200's own eps = 0.3 certificate and
# theta_max_with_collisions columns are kept alongside so the two runs compare cell for
# cell.
#
# The two constructive arms cost tenths of a second per design; the authors' CLI can sit
# on a 34-face graph for the full 600 s cap. They are therefore run as SEPARATE PASSES
# over the same population, and joined in --mode aggregate on (set, id, sigma), so that a
# slow native pass can never hold up the answer the constructive arms give.
#
# Modes: --mode pop (the 200 random designs, constructive arms, shardable by graph),
# --mode native (the native cells of the same population, shardable), --mode ref (the
# eight Phase-2 reference tilings, both passes), --mode aggregate (join into regime.csv).
#
#   julia --project=Kirigami Kirigami/apps/exp_regime.jl [--mode pop] [--out DIR]
#         [--cli PATH] [--work DIR] [--logdir DIR] [--shard S] [--nshards M]
#         [--merge-shards 4] [--timeout 600] [--patdir DIR] [--weld 1e-4]
#         [--id-lo 2000] [--id-hi 2099] [--minf 20] [--maxf 100] [--limit K] [--regenerate]
#
# Population from the frozen data/corpus/regime_100.json, whose sigma_def IS the
# orientation_defect(m0, sigma_mc, 20 F, 7000 + id) result this driver recomputes
# (data/corpus/README.md: "recomputed exactly as exp_regime.jl does"); `--regenerate`
# rebuilds both the graph (make_graph) and sigma_def (orientation_defect) instead. Outputs
# go to results/regime/ by default. `--limit K`
# stops
# after K graphs of the shard (pop / native modes).
include(joinpath(@__DIR__, "native_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# K9c's referee at a caller-chosen shrink (apps/exp_k9c_range_embedding.jl, verbatim). 1e-9 is the
# shrink K9c's headline "refereed" counts use; characterize's own referee runs at 1e-12.
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

# One arm's answer, as it goes into the CSV.
Base.@kwdef mutable struct Arm
    status::String = "not_run"   # completed / timed_out / crashed / no_solution
    theta::Float64 = 0.0; eps::Float64 = 0.0; ref12::Float64 = 0.0; ref9::Float64 = 0.0
    certified::Bool = false
    cert03::Int = 0               # Native200's eps = 0.3 certificate column
    native_theta::Float64 = -1.0  # their own theta_max_with_collisions (diagnostic)
    binding::String = "n/a"
    secs::Float64 = 0.0
    extra::String = "-"           # k9c provenance
end

# The six shared columns, then the arm-specific tail.
core(a::Arm) = a.status * "," * num(a.theta) * "," * num(a.eps) * "," * (a.certified ? "1" : "0") *
               "," * num(a.ref12) * "," * num(a.ref9) * "," * a.binding * "," * num(a.secs)

# A fresh mesh with `sigma` and its topology (a fresh object so m0's topology arrays are
# never rebuilt in place).
function with_sigma(m0::K.Mesh, sigma::Vector{Int})
    m = K.Mesh(m0.X, m0.faces)
    m.sigma = copy(sigma)
    K.build_topology!(m)
    return m
end

# THE instrument. Everything in this file that reports a Theta_max or an eps goes here.
function score_point(m0::K.Mesh, sigma::Vector{Int}, X::Vector{Vec2})
    a = Arm()
    a.status = "completed"
    co = K.CharacterizeOptions()
    co.referee = true   # shrink 1e-12, grid 4000, 50 bisections -- characterize's default
    ch = K.characterize(m0, sigma, X, co)
    a.theta = ch.theta_max
    a.eps = ch.eps_max
    a.certified = ch.eps_max > 0
    a.cert03 = ch.eps_max >= 0.3 ? 1 : 0
    a.binding = ch.binding
    a.ref12 = ch.theta_bisect
    c = K.make_cut(with_sigma(m0, sigma))
    a.ref9 = referee_theta_shrink(c, X, 1e-9)
    return a
end

const kNativeHeader =
    "set,id,kind,cell,status,theta,eps,cert,ref12,ref9,binding,secs,cert03,nativetheta"

const kHeader =
    "set,id,kind,sigma,N,F,n_split,dim_null," *
    "base_status,base_theta,base_eps,base_cert,base_ref12,base_ref9,base_binding,base_secs," *
    "k9c_status,k9c_theta,k9c_eps,k9c_cert,k9c_ref12,k9c_ref9,k9c_binding,k9c_secs,k9c_prov," *
    "natour_status,natour_theta,natour_eps,natour_cert,natour_ref12,natour_ref9," *
    "natour_binding,natour_secs,natour_cert03,natour_nativetheta," *
    "natcol_status,natcol_theta,natcol_eps,natcol_cert,natcol_ref12,natcol_ref9," *
    "natcol_binding,natcol_secs,natcol_cert03,natcol_nativetheta"

# Runs the two constructive arms on one (mesh, sigma) and writes one CSV row.
function run_design(csv::IO, set_name::AbstractString, id::Int, kind::AbstractString,
                    m0::K.Mesh, sigma::Vector{Int}, sigma_name::AbstractString,
                    nat_our::Arm, nat_col::Arm)
    m = with_sigma(m0, sigma)
    c = K.make_cut(m)

    base = Arm(); k9c = Arm()
    dim_null = -1
    let t = Timer()
        o = K.DesignOptions()
        o.characterize.referee = false   # the point itself is refereed by score_point
        d = K.design_baseline(m, sigma, m.X, o)
        dim_null = d.dim_null
        if !isempty(d.X)
            base = score_point(m, sigma, d.X)
        else
            base.status = isempty(d.status) ? "no_solution" : d.status
        end
        base.secs = s(t)
    end
    let t = Timer()
        o = K.RangeMaxOptions()   # run defaults == the archived K9c run's settings
        o.characterize.referee = false
        r = K.design_range_max(m, sigma, m.X, o)
        if !isempty(r.design.X)
            k9c = score_point(m, sigma, r.design.X)
            k9c.extra = isempty(r.provenance) ? "-" : r.provenance
        else
            k9c.status = isempty(r.design.status) ? "no_solution" : r.design.status
        end
        k9c.secs = s(t)
    end

    print(csv, set_name, ",", id, ",", kind, ",", sigma_name, ",", K.n_vertices(m), ",",
          K.n_faces(m), ",", K.n_split(c), ",", dim_null, ",", core(base), ",", core(k9c), ",",
          k9c.extra, ",", core(nat_our), ",", nat_our.cert03, ",", num(nat_our.native_theta),
          ",", core(nat_col), ",", nat_col.cert03, ",", num(nat_col.native_theta), "\n")
    flush(csv)
    println(set_name, " ", id, " ", kind, " ", sigma_name, " F=", K.n_faces(m), " k=", dim_null,
            " | base th=", fmt_g(base.theta), " | k9c th=", fmt_g(k9c.theta), " eps=",
            fmt_g(k9c.eps), " (", k9c.extra, ")", " | natour ", nat_our.status, " th=",
            fmt_g(nat_our.theta), " | natcol ", nat_col.status, " th=", fmt_g(nat_col.theta),
            " [", fmt_g(base.secs + k9c.secs), " s]")
    flush(stdout)
    return nothing
end

# One native cell, run through the shared exp_native200 code path and then scored by the
# SAME instrument the two constructive arms use.
function native_arm(id::Int, kind::AbstractString, m0::K.Mesh,
                    sigma::Union{Vector{Int},Nothing}, variant::AbstractString,
                    cli::AbstractString, work::AbstractString, timeout_s::Int,
                    logdir::AbstractString)
    a = Arm()
    r, X, used_sigma = process_cell(id, kind, m0, sigma, variant, cli, work, timeout_s, logdir)
    a.status = r.status
    a.secs = r.secs
    a.native_theta = r.native_theta_collisions
    if r.status == "completed" && r.embedding_ok && X !== nothing && !isempty(X) &&
       used_sigma !== nothing && length(used_sigma) == K.n_faces(m0)
        secs = a.secs
        nth = a.native_theta
        a = score_point(m0, used_sigma, X)
        a.secs = secs
        a.native_theta = nth
    elseif r.status == "completed"
        a.status = "no_embedding"
    end
    return a
end

# ---------------------------------------------------------------------------
# Everything REGIME.md quotes is derived HERE, from results/regime/regime.csv on disk and
# from nothing this process kept in memory, and written to results/regime/summary.txt.

Base.@kwdef mutable struct ArmStats
    n::Int = 0; deployable::Int = 0; certified::Int = 0; eps01::Int = 0; refereed9::Int = 0
    eps::Vector{Float64} = Float64[]
end

function add_arm(S::ArmStats, theta::Float64, eps::Float64, ref9::Float64)
    S.n += 1
    theta > 1e-9 && (S.deployable += 1)
    (theta > 1e-9 && ref9 > 1e-9) && (S.refereed9 += 1)
    if eps > 0
        S.certified += 1
        push!(S.eps, eps)
        eps >= 0.1 && (S.eps01 += 1)
    end
    return nothing
end

med_of(v::Vector{Float64}) = isempty(v) ? NaN : sort(v)[length(v) ÷ 2 + 1]
max_of(v::Vector{Float64}) = isempty(v) ? NaN : maximum(v)
cell(v::Float64) = isnan(v) ? "--" : fx(v, 3)

# Wilson score interval for a binomial proportion at 95 %.
function wilson(k::Int, n::Int)
    n == 0 && return (NaN, NaN)
    z = 1.959963984540054
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    r = z * sqrt(p * (1 - p) / n + z * z / (4.0 * n * n))
    return ((c - r) / d, (c + r) / d)
end

# std::map iteration order: keys sorted.
sorted_keys(d::AbstractDict) = sort!(collect(keys(d)))
getstats!(d::Dict{String,ArmStats}, k::String) = get!(d, k) do; ArmStats(); end

function write_summary(outdir::AbstractString)
    path = joinpath(outdir, "regime.csv")
    if !isfile(path)
        println(stderr, "write_summary: no regime.csv")
        return
    end
    lines = readlines(path)
    head = split_csv(lines[1])
    col = Dict{String,Int}(h => i for (i, h) in enumerate(head))

    arms = ["base", "natcol", "natour", "k9c"]
    arm_label = ["Eq. (6) baseline (design_baseline)", "native pipeline (their colouring)",
                 "native pipeline (our sigma)", "k9c range-max (design_range_max)"]
    by_set = Dict{String,Dict{String,ArmStats}}()     # set -> arm -> stats
    best_by_set = Dict{String,ArmStats}()
    # native status counts: set|kind|cell -> status -> n
    nat_status = Dict{String,Dict{String,Int}}()
    # yield by |F| bin: bin|arm -> stats
    by_bin = Dict{String,ArmStats}()
    pop_F = Int[]
    pop_F_kind = Dict{String,Vector{Int}}()
    pop_seen = Set{String}()

    nrows = 0
    for line in lines[2:end]
        isempty(line) && continue
        f = split_csv(line)
        length(f) != length(head) && continue
        nrows += 1
        set_name = f[col["set"]]; kind = f[col["kind"]]
        F = parse(Int, f[col["F"]])
        if set_name == "pop" && !(f[col["id"]] in pop_seen)
            push!(pop_seen, f[col["id"]])
            push!(pop_F, F)
            push!(get!(pop_F_kind, kind, Int[]), F)
        end
        best_theta = 0.0; best_eps = 0.0; best_ref9 = 0.0
        for a in arms
            st = f[col[a * "_status"]]
            th = parse(Float64, f[col[a * "_theta"]])
            ep = parse(Float64, f[col[a * "_eps"]])
            r9 = parse(Float64, f[col[a * "_ref9"]])
            st == "not_run" && continue   # a cell that was never dispatched is not a zero
            add_arm(getstats!(get!(by_set, set_name, Dict{String,ArmStats}()), a), th, ep, r9)
            if set_name == "pop"
                bin = F < 50 ? "[20,50)" : "[50,100]"
                add_arm(getstats!(by_bin, bin * "|" * a), th, ep, r9)
            end
            if th > best_theta || (th == best_theta && ep > best_eps)
                best_theta = th; best_eps = ep; best_ref9 = r9
            end
        end
        add_arm(getstats!(best_by_set, set_name), best_theta, best_eps, best_ref9)
        for a in ("natour", "natcol")
            cellname = a == "natcol" ? "native_col" : f[col["sigma"]]
            d = get!(nat_status, set_name * "|" * kind * "|" * cellname, Dict{String,Int}())
            d[f[col[a * "_status"]]] = get(d, f[col[a * "_status"]], 0) + 1
        end
    end

    o = open(joinpath(outdir, "summary.txt"), "w")
    print(o, "WP7a regime summary -- every number below is a function of ", outdir, "/regime.csv\n")
    print(o, "design rows read: ", nrows, "\n\n")

    print(o, "== population: |F| distribution over the ", length(pop_F),
          " random graphs (ids 2000..2099, make_graph(id, 20, 100, 1400)) ==\n")
    if !isempty(pop_F)
        s2 = sort(pop_F)
        print(o, "  all: min ", s2[1], ", median ", s2[length(s2) ÷ 2 + 1], ", max ", s2[end],
              ", in [20,100] ", (s2[1] >= 20 && s2[end] <= 100 ? "YES" : "NO"), "\n")
        for k in sorted_keys(pop_F_kind)
            v = sort(pop_F_kind[k])
            print(o, "  ", k, ": n ", length(v), ", min ", v[1], ", median ", v[length(v) ÷ 2 + 1],
                  ", max ", v[end], "\n")
        end
        lo = count(F -> F < 50, pop_F)
        print(o, "  bin [20,50): ", lo, " graphs; bin [50,100]: ", length(pop_F) - lo, " graphs\n")
    end

    for sname in sorted_keys(by_set)
        sv = by_set[sname]
        print(o, "\n== arm table, set = ", sname, " ==\n")
        print(o, "| arm | N | deployable (Theta_max>0) | refereed (bisect 1e-9 > 0) | certified ",
              "(eps_max>0) | eps_max >= 0.1 rad | median eps_max | max eps_max |\n")
        print(o, "|---|---|---|---|---|---|---|---|\n")
        for (i, a) in enumerate(arms)
            if !haskey(sv, a)
                print(o, "| ", arm_label[i], " | 0 (none dispatched) | -- | -- | -- | -- | -- | -- |\n")
                continue
            end
            S = sv[a]
            print(o, "| ", arm_label[i], " | ", S.n, " | ", S.deployable, " | ", S.refereed9, " | ",
                  S.certified, " | ", S.eps01, " | ", cell(med_of(S.eps)), " | ", cell(max_of(S.eps)), " |\n")
        end
        B = getstats!(best_by_set, sname)
        print(o, "| best per design (max over the arms run) | ", B.n, " | ", B.deployable, " | ",
              B.refereed9, " | ", B.certified, " | ", B.eps01, " | ", cell(med_of(B.eps)), " | ",
              cell(max_of(B.eps)), " |\n")
    end

    print(o, "\n== native cell status, per set / family / cell ==\n")
    print(o, "| set | family | cell | ")
    sts = ["completed", "timed_out", "crashed", "no_embedding", "not_run"]
    for st in sts
        print(o, st, " | ")
    end
    print(o, "\n|---|---|---|---|---|---|---|---|\n")
    for key in sorted_keys(nat_status)
        parts = split(key, '|')
        a = parts[1]; b = length(parts) >= 2 ? parts[2] : ""; c2 = length(parts) >= 3 ? parts[3] : ""
        print(o, "| ", a, " | ", b, " | ", c2, " | ")
        for st in sts
            print(o, get(nat_status[key], st, 0), " | ")
        end
        print(o, "\n")
    end

    print(o, "\n== yield vs |F|, random population, Wilson 95 % intervals ==\n")
    print(o, "| bin | arm | N | deployable | yield | Wilson 95 % |\n|---|---|---|---|---|---|\n")
    for bin in ("[20,50)", "[50,100]"), (i, a) in enumerate(arms)
        haskey(by_bin, bin * "|" * a) || continue
        S = by_bin[bin * "|" * a]
        lo, hi = wilson(S.deployable, S.n)
        print(o, "| ", bin, " | ", arm_label[i], " | ", S.n, " | ", S.deployable, " | ",
              fx(100.0 * S.deployable / max(1, S.n), 1), " % | [", fx(100 * lo, 1), ", ",
              fx(100 * hi, 1), "] % |\n")
    end
    close(o)
    println(stderr, "wrote ", outdir, "/summary.txt")
    return nothing
end

# sigma_def of one graph: the frozen row's recorded orientation_defect result, or (with
# `--regenerate`, or when the row carries none) recomputed.
function regime_sigma_def(id::Int, m0::K.Mesh, sigma_mc::Vector{Int},
                          frozen::Union{Vector{Int},Nothing})
    (frozen !== nothing && length(frozen) == K.n_faces(m0)) && return frozen
    # K5's defect-minimising sigma, recomputed here (this population has no archived
    # results/experiments/k5/sigma directory): exp_k5_orientation.jl's cap of 20 |F| attempts, seed
    # 7000 + id, exactly the promotion documented in method/design.jl.
    dor = K.orientation_defect(m0, sigma_mc, 20 * K.n_faces(m0), 7000 + id)
    def_ok = dor.ok && length(dor.sigma) == K.n_faces(m0)
    return def_ok ? dor.sigma : sigma_mc
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    mode = "pop"
    outdir = joinpath(REPO, "results", "regime")
    cli = joinpath(REPO, "baseline", "native", "build_fixed", "tuttekiri_cli")   # F24-patched build, baseline/README.md
    work = "/tmp/kiri_regime"
    logdir = ""
    patdir = joinpath(REPO, "baseline", "kirigami_tessellations")
    weld_frac = 1e-4
    shard = 0; nshards = 1; timeout_s = 600; nshards_in = 4
    id_lo = 2000; id_hi = 2099; minf = 20; maxf = 100; ncap = 1400
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--mode" && i < length(args); mode = args[i+1]; i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--cli" && i < length(args); cli = args[i+1]; i += 2
        elseif a == "--work" && i < length(args); work = args[i+1]; i += 2
        elseif a == "--logdir" && i < length(args); logdir = args[i+1]; i += 2
        elseif a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshards" && i < length(args); nshards = arg_i(args[i+1]); i += 2
        elseif a == "--merge-shards" && i < length(args); nshards_in = arg_i(args[i+1]); i += 2
        elseif a == "--timeout" && i < length(args); timeout_s = arg_i(args[i+1]); i += 2
        elseif a == "--patdir" && i < length(args); patdir = args[i+1]; i += 2
        elseif a == "--weld" && i < length(args); weld_frac = arg_f(args[i+1]); i += 2
        elseif a == "--id-lo" && i < length(args); id_lo = arg_i(args[i+1]); i += 2
        elseif a == "--id-hi" && i < length(args); id_hi = arg_i(args[i+1]); i += 2
        elseif a == "--minf" && i < length(args); minf = arg_i(args[i+1]); i += 2
        elseif a == "--maxf" && i < length(args); maxf = arg_i(args[i+1]); i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    isempty(logdir) && (logdir = joinpath(outdir, "logs"))
    mkpath(outdir)
    mkpath(logdir)
    mkpath(work)

    # ---- join the constructive and native passes into regime.csv ---------------
    # The native pass writes one row per CELL, keyed (set, id, cell) with cell in
    # {sigma_mc, sigma_def, native_col}; a design row takes its `natour_*` columns from the
    # cell named by its own sigma and its `natcol_*` columns from the graph's native_col
    # cell. Cells that were never run stay "not_run" and are counted as such in REGIME.md.
    if mode == "aggregate"
        ncell = Dict{String,Vector{String}}()   # key -> 10 native fields
        nparts = [joinpath(outdir, "native_shard_" * string(sh) * ".csv") for sh in 0:nshards_in-1]
        push!(nparts, joinpath(outdir, "native_ref.csv"))
        push!(nparts, joinpath(outdir, "native_patterns.csv"))
        ncells = 0
        for p in nparts
            if !isfile(p)
                println(stderr, "no native pass file ", p, " (skipped)")
                continue
            end
            for (li, line) in enumerate(readlines(p))
                li == 1 && continue
                isempty(line) && continue
                fl = split_csv(line)
                length(fl) < 14 && continue
                # status,theta,eps,cert,ref12,ref9,binding,secs then cert03,nativetheta
                ncell[fl[1] * "|" * fl[2] * "|" * fl[4]] = fl[5:14]
                ncells += 1
            end
        end
        out = open(joinpath(outdir, "regime.csv"), "w")
        print(out, kHeader, "\n")
        rows = 0; patched_our = 0; patched_col = 0
        parts = [joinpath(outdir, "shard_" * string(sh) * ".csv") for sh in 0:nshards_in-1]
        push!(parts, joinpath(outdir, "ref.csv"))
        push!(parts, joinpath(outdir, "patterns.csv"))
        for p in parts
            if !isfile(p)
                println(stderr, "missing ", p, " (skipped)")
                continue
            end
            for (li, line) in enumerate(readlines(p))
                li == 1 && continue
                isempty(line) && continue
                fl = split_csv(line)
                if length(fl) != 45
                    println(stderr, "bad row width ", length(fl))
                    continue
                end
                base = fl[1] * "|" * fl[2] * "|"
                io_ = get(ncell, base * fl[4], nothing)
                if io_ !== nothing
                    fl[26:35] = io_   # 0-based columns 25 + k, k = 0..9
                    patched_our += 1
                end
                ic = get(ncell, base * "native_col", nothing)
                if ic !== nothing
                    fl[36:45] = ic    # 0-based columns 35 + k
                    patched_col += 1
                end
                print(out, join(fl, ","), "\n")
                rows += 1
            end
        end
        println(stderr, "regime.csv: ", rows, " design rows, ", ncells, " native cells read, natour patched ",
                patched_our, ", natcol patched ", patched_col)
        close(out)
        write_summary(outdir)
        return 0
    end

    wall = Timer()

    # ---- helper: write one native-cell row -------------------------------------
    function native_row(o::IO, set_name, id, kind, cellname, a::Arm)
        print(o, set_name, ",", id, ",", kind, ",", cellname, ",", core(a), ",", a.cert03, ",",
              num(a.native_theta), "\n")
        flush(o)
        println("[native] ", set_name, " ", id, " ", cellname, " -> ", a.status, " theta=",
                fmt_g(a.theta), " (", fmt_g(a.secs), " s)")
        flush(stdout)
    end

    # ---- the eight Phase-2 reference tilings -----------------------------------
    if mode == "ref" || mode == "native_ref"
        nat_pass = mode == "native_ref"
        path = joinpath(outdir, nat_pass ? "native_ref.csv" : "ref.csv")
        csv = open(path, "w")
        print(csv, nat_pass ? kNativeHeader : kHeader, "\n")
        flush(csv)
        cs = reference_cases(; regenerate = regenerate)
        idx = 0   # 0-based case index (the CSV numbering)
        for rc in cs
            m0 = K.Mesh(rc.mesh.X, rc.mesh.faces)
            m0.sigma = copy(rc.mesh.sigma)
            m0.periodic = rc.mesh.periodic
            K.build_topology!(m0)
            sig = copy(m0.sigma)
            if length(sig) != K.n_faces(m0)
                idx += 1
                continue
            end
            if nat_pass
                native_row(csv, "ref", idx, rc.name, "sigma_given",
                           native_arm(idx, rc.name, m0, sig, "sigma_given", cli, work, timeout_s, logdir))
                native_row(csv, "ref", idx, rc.name, "native_col",
                           native_arm(idx, rc.name, m0, nothing, "native", cli, work, timeout_s, logdir))
            else
                run_design(csv, "ref", idx, rc.name, m0, sig, "sigma_given", Arm(), Arm())
            end
            idx += 1
        end
        close(csv)
        open(joinpath(outdir, mode * "_done.txt"), "w") do o
            print(o, mode, " done, wall ", fmt_g(s(wall)), " s\n")
        end
        return 0
    end

    # ---- the 2025 paper's own patterns -----------------------------------------
    # baseline/kirigami_tessellations ships two kinds of "pattern file". The 16 JSONs under
    # its code/data/patterns are NOT geometry: each is a list of parallel-line GROUPS plus a
    # gridType/gridSize, i.e. the input to the 2025 web demo's "cut tiling into kirigami"
    # step, and turning one into a planar graph means reimplementing that step, not
    # importing a file. The five SVGs under fabrication_patterns ARE polygon soup -- every
    # tile edge is an independent <line> element -- and those are what core/import_soup.jl
    # welds. Both facts are recorded in patterns_load.csv, which lists every candidate file
    # and whether it loaded.
    if mode == "patterns" || mode == "native_patterns"
        nat_p = mode == "native_patterns"
        path = joinpath(outdir, nat_p ? "native_patterns.csv" : "patterns.csv")
        csv = open(path, "w")
        print(csv, nat_p ? kNativeHeader : kHeader, "\n")
        flush(csv)
        ld = open(joinpath(outdir, "patterns_load.csv"), nat_p ? "a" : "w")
        nat_p || print(ld, "file,kind,loaded,status,segments,welded_vertices,pruned,faces,faces_dropped,",
                       "min_area,max_area,N_kept,F_kept,note\n")

        svgs = String[]; jsons = String[]
        let d = joinpath(patdir, "fabrication_patterns")
            isdir(d) && for e in readdir(d; join = true)
                endswith(e, ".svg") && push!(svgs, e)
            end
        end
        let d = joinpath(patdir, "code", "data", "patterns")
            isdir(d) && for e in readdir(d; join = true)
                endswith(e, ".json") && push!(jsons, e)
            end
        end
        sort!(svgs); sort!(jsons)

        if !nat_p
            for j in jsons
                print(ld, basename(j), ",line_group_spec,0,not_polygon_soup,0,0,0,0,0,0,0,0,0,",
                      "groups+gridType: input to the 2025 web demo cutting step; no geometry\n")
            end
        end

        idx = 0
        for sv in svgs
            name = first(splitext(basename(sv)))
            sr = K.import_svg_soup(sv, weld_frac, 0.0)
            m0 = K.Mesh()
            usable = false
            if sr.ok
                m0 = K.largest_component(sr.mesh)
                K.build_topology!(m0)
                usable = K.n_faces(m0) >= 4
            end
            nat_p || print(ld, name, ".svg,fabrication_svg,", usable ? 1 : 0, ",", sr.status, ",",
                           sr.n_segments, ",", sr.n_vertices_welded, ",", sr.n_pruned, ",", sr.n_faces, ",",
                           sr.n_faces_dropped, ",", num(sr.min_face_area), ",", num(sr.max_face_area), ",",
                           usable ? K.n_vertices(m0) : 0, ",", usable ? K.n_faces(m0) : 0,
                           ",largest_edge_connected_component\n")
            if !usable
                idx += 1
                continue
            end
            rng = K.MT19937(UInt32(20260908) + UInt32(idx))
            m0.sigma = K.assign_orientation_relaxation(m0, rng, 4, 300, 90).sigma
            if length(m0.sigma) != K.n_faces(m0)
                idx += 1
                continue
            end
            K.build_topology!(m0)
            if nat_p
                native_row(csv, "pat2025", idx, name, "sigma_mc",
                           native_arm(idx, name, m0, m0.sigma, "sigma_mc", cli, work, timeout_s, logdir))
                native_row(csv, "pat2025", idx, name, "native_col",
                           native_arm(idx, name, m0, nothing, "native", cli, work, timeout_s, logdir))
            else
                run_design(csv, "pat2025", idx, name, m0, m0.sigma, "sigma_mc", Arm(), Arm())
            end
            idx += 1
        end
        close(csv); close(ld)
        open(joinpath(outdir, mode * "_done.txt"), "w") do o
            print(o, mode, " done, wall ", fmt_g(s(wall)), " s\n")
        end
        return 0
    end

    # ---- the random paper-regime population ------------------------------------
    nat_pass = mode == "native"
    csv_path = joinpath(outdir, (nat_pass ? "native_shard_" : "shard_") * string(shard) * ".csv")
    csv = open(csv_path, "w")
    print(csv, nat_pass ? kNativeHeader : kHeader, "\n")
    flush(csv)

    fdist = nothing
    if !nat_pass
        fdist = open(joinpath(outdir, "faces_shard_" * string(shard) * ".csv"), "w")
        print(fdist, "id,kind,N,F,n_split_mc,n_split_def,dim_null_mc\n")
    end

    (regenerate || (minf == 20 && maxf == 100 && ncap == 1400)) ||
        @warn "frozen regime_100 was built with make_graph(id, 20, 100, 1400); use --regenerate for ($minf, $maxf, $ncap)"
    rows = regenerate ? PopRow[] : load_population("regime_100")
    gidx = -1
    nrun = 0
    for id in id_lo:id_hi
        m0 = K.Mesh(); kind = ""; frozen_def = nothing
        if regenerate
            g = K.make_graph(id, minf, maxf, ncap)
            if !g.ok
                println(stderr, "make_graph failed for id ", id)
                continue
            end
            m0 = g.mesh; kind = g.kind
        else
            ri = findfirst(r -> r.id == id, rows)
            if ri === nothing || !rows[ri].ok
                println(stderr, "make_graph failed for id ", id)
                continue
            end
            m0 = rows[ri].mesh; kind = rows[ri].kind; frozen_def = rows[ri].sigma_def
        end
        gidx += 1
        (nshards > 1 && gidx % nshards != shard) && continue
        nrun >= limit && break
        nrun += 1
        K.build_topology!(m0)
        sigma_mc = copy(m0.sigma)
        sigma_def = regime_sigma_def(id, m0, sigma_mc, frozen_def)

        if nat_pass
            native_row(csv, "pop", id, kind, "native_col",
                       native_arm(id, kind, m0, nothing, "native", cli, work, timeout_s, logdir))
            native_row(csv, "pop", id, kind, "sigma_mc",
                       native_arm(id, kind, m0, sigma_mc, "sigma_mc", cli, work, timeout_s, logdir))
            native_row(csv, "pop", id, kind, "sigma_def",
                       native_arm(id, kind, m0, sigma_def, "sigma_def", cli, work, timeout_s, logdir))
            continue
        end

        let a = with_sigma(m0, sigma_mc), b = with_sigma(m0, sigma_def)
            print(fdist, id, ",", kind, ",", K.n_vertices(m0), ",", K.n_faces(m0), ",",
                  K.n_split(K.make_cut(a)), ",", K.n_split(K.make_cut(b)), ",-1\n")
            flush(fdist)
        end
        run_design(csv, "pop", id, kind, m0, sigma_mc, "sigma_mc", Arm(), Arm())
        run_design(csv, "pop", id, kind, m0, sigma_def, "sigma_def", Arm(), Arm())
    end
    close(csv)
    fdist === nothing || close(fdist)
    open(joinpath(outdir, (nat_pass ? "native_shard_" : "shard_") * string(shard) * "_done.txt"), "w") do o
        print(o, "wall ", fmt_g(s(wall)), " s\n")
    end
    println(stderr, mode, " shard ", shard, " done in ", fmt_g(s(wall)), " s -> ", csv_path)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
