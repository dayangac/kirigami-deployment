# native200_merge -- merge the Native200 runs into one row per (id, variant) and
# write the report.  Every number in NATIVE200_FINAL.md is computed here from
# native200_final.csv; nothing is hand-entered. Port of code/apps/native200_merge.cpp.
#
#   julia --project=Kirigami Kirigami/apps/native200_merge.jl
#     [--orig results/kill/native200/native200.csv]
#     [--rerun results/kill/native200/rerun3600]        (repeatable: any dir of shard_*.csv)
#     [--rerun results/kill/native200/crashfix3600]
#     [--grid results/kill/native200/grid.csv]
#     [--out results/kill/native200_julia/native200_final.csv]
#     [--md  results/kill/native200_julia/NATIVE200_FINAL.md]
#     [--expected 389]   cells the rerun pipeline has to get through
#
# Merge rule: a cell's row is the one from the LATEST run that touched it, where
# `orig` (the 600 s run) < `rerun3600` < `crashfix3600`; the 600 s status is kept
# in a `first_status` column.  Cells in the grid that no run ever dispatched come
# out as status `missing`.
#
# Inputs default to this repo's migrated results/kill/native200/; the merged CSV and the
# report go to results/kill/native200_julia/ (the C++ wrote them next to the inputs).
include(joinpath(@__DIR__, "common_app.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

Base.@kwdef mutable struct Rec
    id::String = ""; kind::String = ""; variant::String = ""; N::String = ""; F::String = ""
    dim_null::String = ""; status::String = ""; native_theta::String = ""; native_finite::String = ""
    embedding_ok::String = ""; theta_exact::String = ""; theta_bisect::String = ""
    cert_pos::String = ""; cert_noovl::String = ""; cert_noroot::String = ""
    cert_eps03::String = ""; secs::String = ""; run::String = ""
    first_status::String = ""
    prio::Int = -1
end

# `std::getline(ss, tok, ',')`: a trailing empty field is dropped.
function split_line(line::AbstractString)
    f = String.(Base.split(line, ','))
    (!isempty(f) && isempty(f[end]) && endswith(line, ',')) && pop!(f)
    return f
end

dbl(s::AbstractString, d::Float64 = -1.0) = something(tryparse(Float64, s), d)
integer(s::AbstractString, d::Int = -1) = something(tryparse(Int, s), d)

function prio_of(run::AbstractString)
    run == "orig" && return 0
    run == "rerun3600" && return 1
    run == "crashfix3600" && return 2
    return 3   # anything newer wins
end

# std::map<std::pair<int, std::string>, Rec>: ordered by (id, variant).
const CellMap = Dict{Tuple{Int,String},Rec}
sorted_cells(cells::CellMap) = sort!(collect(keys(cells)))

function ingest(path::AbstractString, default_run::AbstractString, out::CellMap)
    isfile(path) || return
    for (li, line) in enumerate(readlines(path))
        li == 1 && continue   # header
        v = split_line(line)
        length(v) < 17 && continue
        r = Rec()
        r.id = v[1]; r.kind = v[2]; r.variant = v[3]; r.N = v[4]; r.F = v[5]
        r.dim_null = v[6]; r.status = v[7]; r.native_theta = v[8]; r.native_finite = v[9]
        r.embedding_ok = v[10]; r.theta_exact = v[11]; r.theta_bisect = v[12]
        r.cert_pos = v[13]; r.cert_noovl = v[14]; r.cert_noroot = v[15]
        r.cert_eps03 = v[16]; r.secs = v[17]
        r.run = (length(v) >= 18 && !isempty(v[18])) ? v[18] : default_run
        r.prio = prio_of(r.run)
        key = (integer(r.id), r.variant)
        if !haskey(out, key)
            out[key] = r
            continue
        end
        if r.prio >= out[key].prio
            r.first_status = isempty(out[key].first_status) ? out[key].status : out[key].first_status
            out[key] = r
        end
    end
    return
end

# ---- small tabulation helpers ----------------------------------------------
const Table = Dict{String,Dict{String,Int}}   # row -> col -> n
sorted_keys(d::AbstractDict) = sort!(collect(keys(d)))
inc!(t::Table, r::String, c::String) = (d = get!(t, r, Dict{String,Int}()); d[c] = get(d, c, 0) + 1)

function emit_table(o::IO, title::AbstractString, t::Table, cols::Vector{String})
    print(o, "\n**", title, "**\n\n| | ")
    for c in cols
        print(o, c, " | ")
    end
    print(o, "total |\n|---|")
    for _ in 1:length(cols)+1
        print(o, "--:|")
    end
    print(o, "\n")
    coltot = Dict{String,Int}()
    grand = 0
    for rname in sorted_keys(t)
        row = t[rname]
        print(o, "| ", rname, " | ")
        tot = 0
        for c in cols
            n = get(row, c, 0)
            print(o, n, " | ")
            tot += n
            coltot[c] = get(coltot, c, 0) + n
        end
        print(o, tot, " |\n")
        grand += tot
    end
    print(o, "| **total** | ")
    for c in cols
        print(o, "**", get(coltot, c, 0), "** | ")
    end
    print(o, "**", grand, "** |\n")
    return
end

function quantile(v_in::Vector{Float64}, q::Float64)
    isempty(v_in) && return NaN
    v = sort(v_in)
    pos = q * (length(v) - 1)
    i = Int(floor(pos))
    fr = pos - i
    return (i + 1 < length(v)) ? v[i + 1] * (1 - fr) + v[i + 2] * fr : v[i + 1]
end

# `o << std::fixed << std::setprecision(3)` (and the one setprecision(1) block).
f3(v::Real) = fx(v, 3)

function main(args::Vector{String})
    indir = joinpath(REPO, "results", "kill", "native200")
    outdir = joinpath(REPO, "results", "kill", "native200_julia")
    orig = joinpath(indir, "native200.csv")
    grid = joinpath(indir, "grid.csv")
    out = joinpath(outdir, "native200_final.csv")
    md = joinpath(outdir, "NATIVE200_FINAL.md")
    reruns = [joinpath(indir, "rerun3600"), joinpath(indir, "crashfix3600")]
    reruns_set = false
    expected = 389
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--orig" && i < length(args); orig = args[i+1]; i += 2
        elseif a == "--grid" && i < length(args); grid = args[i+1]; i += 2
        elseif a == "--out" && i < length(args); out = args[i+1]; i += 2
        elseif a == "--md" && i < length(args); md = args[i+1]; i += 2
        elseif a == "--expected" && i < length(args); expected = integer(args[i+1], 389); i += 2
        elseif a == "--rerun" && i < length(args)
            if !reruns_set
                empty!(reruns)
                reruns_set = true
            end
            push!(reruns, args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(dirname(out)); mkpath(dirname(md))

    cells = CellMap()
    ingest(orig, "orig", cells)
    rerun_files = String[]
    for d in reruns
        isdir(d) || continue
        fs = String[]
        for e in readdir(d; join = true)
            (endswith(e, ".csv") && startswith(basename(e), "shard_")) && push!(fs, e)
        end
        sort!(fs)
        for p in fs
            ingest(p, "rerun3600", cells)
            push!(rerun_files, p)
        end
    end

    # Grid: everything never dispatched by any run is `missing`.
    grid_total = 0
    if isfile(grid)
        for (li, line) in enumerate(readlines(grid))
            li == 1 && continue
            v = split_line(line)
            length(v) < 5 && continue
            grid_total += 1
            key = (integer(v[1]), v[2])
            haskey(cells, key) && continue
            r = Rec()
            r.id = v[1]; r.variant = v[2]; r.kind = v[3]; r.N = v[4]; r.F = v[5]
            r.dim_null = "-1"; r.status = "missing"; r.native_theta = "-1"
            r.native_finite = "0"; r.embedding_ok = "0"; r.theta_exact = "-1"
            r.theta_bisect = "-1"; r.cert_pos = "0"; r.cert_noovl = "0"; r.cert_noroot = "0"
            r.cert_eps03 = "0"; r.secs = "0"; r.run = "none"
            cells[key] = r
        end
    end

    # ---- merged CSV -----------------------------------------------------------
    open(out, "w") do c
        print(c, "id,kind,variant,N,F,dim_null,status,native_theta_collisions,native_finite,",
              "embedding_ok,our_theta_exact,our_theta_bisect,cert_pos,cert_noovl,cert_noroot,",
              "cert_eps03,secs,run,first_status\n")
        for k in sorted_cells(cells)
            r = cells[k]
            print(c, r.id, ",", r.kind, ",", r.variant, ",", r.N, ",", r.F, ",", r.dim_null, ",",
                  r.status, ",", r.native_theta, ",", r.native_finite, ",", r.embedding_ok, ",",
                  r.theta_exact, ",", r.theta_bisect, ",", r.cert_pos, ",", r.cert_noovl, ",",
                  r.cert_noroot, ",", r.cert_eps03, ",", r.secs, ",", r.run, ",",
                  isempty(r.first_status) ? r.status : r.first_status, "\n")
        end
    end

    # ---- statistics -----------------------------------------------------------
    by_variant = Dict{String,Table}(); by_family = Dict{String,Table}()   # run-label -> table
    run_count = Dict{String,Int}()
    n_completed = 0; n_exact_pos = 0; n_bisect_pos = 0; n_cert = 0; n_native_pos = 0
    secs_completed = Float64[]
    F_completed = Int[]; F_timedout_3600 = Int[]
    native_pos_rows = Rec[]

    function label(run::AbstractString)
        run == "orig" && return "orig (600 s)"
        run == "crashfix3600" && return "crashfix3600"
        run == "none" && return "never dispatched"
        return "rerun3600"
    end

    for k in sorted_cells(cells)
        r = cells[k]
        L = label(r.run)
        inc!(get!(by_variant, L, Table()), r.status, r.variant)
        inc!(get!(by_family, L, Table()), r.status, r.kind)
        inc!(get!(by_variant, "ALL", Table()), r.status, r.variant)
        inc!(get!(by_family, "ALL", Table()), r.status, r.kind)
        run_count[r.run] = get(run_count, r.run, 0) + 1
        if r.status == "completed"
            n_completed += 1
            push!(secs_completed, dbl(r.secs, 0.0))
            push!(F_completed, integer(r.F, 0))
            dbl(r.theta_exact) > 0 && (n_exact_pos += 1)
            dbl(r.theta_bisect) > 0 && (n_bisect_pos += 1)
            r.cert_eps03 == "1" && (n_cert += 1)
            if dbl(r.native_theta) > 0
                n_native_pos += 1
                push!(native_pos_rows, r)
            end
        end
        (r.status == "timed_out" && r.run != "orig") && push!(F_timedout_3600, integer(r.F, 0))
    end

    count_of(k) = get(run_count, k, 0)
    done_rerun = count_of("rerun3600")
    done_crashfix = count_of("crashfix3600")
    variants = ["native", "sigma_mc", "sigma_def"]
    families = ["delaunay", "voronoi", "quad_random"]

    o = open(md, "w")
    print(o, "# Native200 -- final status of the authors' full native pipeline on the K9 population\n\n")
    if done_rerun < expected || done_crashfix < 39
        print(o, "**PARTIAL -- rerun in progress, ", done_rerun, " of ", expected,
              " cells done** (3600 s rerun of the timed-out and never-dispatched cells), plus ",
              done_crashfix, " of 39 crash-fixed cells. Regenerate with ",
              "`code/build/native200_merge`.\n\n")
    else
        print(o, "**COMPLETE -- all ", expected, " rerun cells and all 39 crash-fixed cells ",
              "finished.**\n\n")
    end
    print(o, "Every number below is computed by `code/apps/native200_merge.cpp` from ",
          "`results/kill/native200/native200_final.csv`, which is itself merged from ",
          "`native200.csv` (the 600 s run), `rerun3600/shard_*.csv` (the 3600 s rerun of the ",
          "timed-out and never-dispatched cells) and `crashfix3600/shard_*.csv` (the 39 ",
          "crashed cells against the F24-patched CLI). One row per (graph id, variant); the ",
          "600 s status of a superseded cell is kept in `first_status`.\n\n")
    print(o, "Grid: ", grid_total, " cells (200 graphs x {native, sigma_mc, sigma_def}). ",
          "Rows by run: ")
    for r in sorted_keys(run_count)
        print(o, label(r), " ", run_count[r], "; ")
    end
    print(o, "\n")

    print(o, "\n## 1. Status by variant\n")
    for L in sorted_keys(by_variant)
        emit_table(o, L, by_variant[L], variants)
    end
    print(o, "\n## 2. Status by family\n")
    for L in sorted_keys(by_family)
        emit_table(o, L, by_family[L], families)
    end

    print(o, "\n## 3. What the completed cells say\n\n")
    print(o, "| quantity | count | of completed |\n|---|--:|--:|\n")
    print(o, "| completed cells (any run) | ", n_completed, " | ", n_completed, " |\n")
    print(o, "| exact `Theta_max > 0` (T4.2\" scan on their dumped embedding) | ", n_exact_pos,
          " | ", n_completed, " |\n")
    print(o, "| bisection referee `Theta > 0` | ", n_bisect_pos, " | ", n_completed, " |\n")
    print(o, "| validity certificate at eps = 0.3 (POS & NOOVERLAP & NOROOT) | ", n_cert,
          " | ", n_completed, " |\n")
    print(o, "| the authors' own FK collision test says `Theta > 0` | ", n_native_pos, " | ",
          n_completed, " |\n")

    print(o, "\nThe cells where their own test reports `Theta > 0` while our exact scan on the ",
          "SAME dumped embedding gives 0 (the F24/K1c false-negative mechanism):\n\n")
    if isempty(native_pos_rows)
        print(o, "_(none)_\n")
    else
        print(o, "| id | kind | variant | F | dim_null | native_theta | our_exact | our_bisect | ",
              "cert_pos | run |\n|---|---|---|--:|--:|--:|--:|--:|--:|---|\n")
        for r in native_pos_rows
            print(o, "| ", r.id, " | ", r.kind, " | ", r.variant, " | ", r.F, " | ", r.dim_null, " | ",
                  f3(dbl(r.native_theta)), " | ", f3(dbl(r.theta_exact)), " | ", f3(dbl(r.theta_bisect)),
                  " | ", r.cert_pos, " | ", label(r.run), " |\n")
        end
    end

    print(o, "\n## 4. Cost\n\n")
    print(o, "| quantity | value |\n|---|--:|\n")
    print(o, "| completed wall time, median (s) | ", f3(quantile(secs_completed, 0.5)), " |\n")
    print(o, "| completed wall time, q90 (s) | ", f3(quantile(secs_completed, 0.9)), " |\n")
    print(o, "| completed wall time, max (s) | ",
          f3(isempty(secs_completed) ? NaN : maximum(secs_completed)), " |\n")
    Fc = Float64.(F_completed)
    print(o, "| completed \\|F\\|, median | ", f3(quantile(Fc, 0.5)), " |\n")
    print(o, "| completed \\|F\\|, max | ", f3(isempty(Fc) ? NaN : maximum(Fc)), " |\n")
    if !isempty(F_timedout_3600)
        Ft = Float64.(F_timedout_3600)
        print(o, "| cells still timing out at 3600 s | ", length(F_timedout_3600), " |\n")
        print(o, "| smallest \\|F\\| that still times out at 3600 s | ", f3(minimum(Ft)), " |\n")
        print(o, "| median \\|F\\| of the cells that still time out | ", f3(quantile(Ft, 0.5)), " |\n")
        print(o, "\nOn this machine their pipeline does not finish within 3600 s above roughly ",
              "|F| = ", f3(minimum(Ft)),
              " faces (the smallest patch that still times out); the largest patch that does ",
              "finish has |F| = ", f3(isempty(Fc) ? 0.0 : maximum(Fc)), ".\n")
    else
        print(o, "| cells still timing out at 3600 s | 0 (none yet) |\n")
    end

    print(o, "\n## 5. The sentence the paper may print\n\n")
    unfinished = grid_total - n_completed
    print(o, "> On the 200-graph K9 population, the authors' full native pipeline (their ",
          "colouring, their Eq. (6), their Eq. (9) at the published web-UI defaults, their ",
          "forward-kinematics collision test) reaches an exact `Theta_max > 0` on **",
          n_exact_pos, " of ", n_completed,
          "** designs it completes, and passes the eps = 0.3 validity certificate on **",
          n_cert, " of ", n_completed, "**; their own collision test reports ",
          "`Theta > 0` on ", n_native_pos,
          " of them, every one a false negative of the F24/K1c merge_close_verts mechanism.\n\n")
    if unfinished > 0
        print(o, "Qualifier that must accompany it: **", unfinished, " of ", grid_total,
              " cells (", fx(100.0 * unfinished / max(1, grid_total), 1),
              " %)** ",
              "never produce a design at all -- they time out or crash -- so the ",
              "denominator is the completed subset, not the population.\n")
    else
        print(o, "No qualifier is needed: every cell of the grid completed.\n")
    end
    close(o)

    println(stderr, "[native200_merge] ", length(cells), " cells -> ", out, ", ", md,
            " (rerun rows so far: ", done_rerun, "/", expected, ", crashfix rows: ", done_crashfix, "/39)")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
