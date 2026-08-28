# Native200 -- the authors' FULL native pipeline on ALL 200 K6/K9 graphs.
# Port of code/apps/kill_native200.cpp.
#
# A referee baseline requested by the team lead: not just Eq.(9) (K2b) but the
# authors' whole chain -- THEIR coloring (coloring::initialized_two_face_coloring),
# THEIR make_deployable (Eq. 6), THEIR opt::prevent_intersections (Eq. 9) at the
# web UI's published defaults, and THEIR forward-kinematics collision test for
# theta_max -- run natively via baseline/native/build/tuttekiri_cli on every one
# of the 200 random graphs of kill_common.hpp's population (make_graph(id, 100,
# maxf, 1400), same ids as K1a/K5/K6/K9/K9b).
#
# Three variants per graph, each a separate native run:
#   native   : no `orientation` field in the input JSON for `prevent`; the CLI's
#              `color` subcommand is run first to get THEIR OWN coloring
#              (coloring::initialized_two_face_coloring), whose output (with an
#              `orientation` field written back by io::pattern_to_json) is then fed
#              to `prevent`. This is the fully-native pipeline, their sigma included.
#   sigma_mc : our max-cut sigma (Eq. 1 relaxation) baked into the input JSON's
#              `orientation` field; native solve + native Eq. (9) on top of it.
#   sigma_def: K5's defect-minimising sigma (results/kill/k5/sigma/*.json), same
#              convention as K6.
#
# For every (graph, variant) that completes, THEIR theta_max_with_collisions is
# recorded (known buggy, F24: merge_close_verts fuses hinge duplicates -- kept as a
# diagnostic column only, never as the referee), and the optimised embedding is
# read back (`--dump`) and re-evaluated with OUR OWN instruments: the exact T4.2"
# theta_max (exact_theta_max_overlap), the bisection referee (identical to
# K2a/K2b/K5/K6), and the repaired validity certificate (POS /\ NOOVERLAP /\
# NOROOT) at eps = 0.3, method/contact.jl's validity_certificate.
#
# A per-(graph,variant) wall-clock cap of 600 s is enforced with `perl -e alarm`
# (no GNU coreutils `timeout` on this machine, same trick as kill_k6.cpp). Runs are
# classified completed / timed_out / crashed (nonzero exit, not a timeout) --
# crashes are candidate F24-class upstream bugs and are logged verbatim.
#
# K9's numbers for the SAME 200 graphs (already computed, results/kill/k9/k9.csv:
# columns a_* = phase-A convex embedding, b_* = phase-B solver push, both sigma
# families) are read back for a side-by-side table; nothing there is recomputed.
#
#   julia --project=Kirigami Kirigami/apps/kill_native200.jl [--n 200] [--maxf 800]
#         [--out DIR] [--sigma DIR] [--cli PATH] [--work DIR] [--timeout 600] [--shard S]
#         [--nshards M] [--logdir DIR] [--cells FILE] [--resume CSV]... [--tag orig]
#         [--csv PATH] [--listcells FILE] [--limit K] [--regenerate]
#
# The grid comes from the frozen data/corpus/native200.json (the population index walk
# of make_graph with the archived K5 sigma_def); `--sigma DIR` reads the archived files
# instead of the frozen sigma_def, `--regenerate` rebuilds the graphs through make_graph.
# The CLI defaults to the reference repo's F24-patched build
# (~/Documents/kirigami-experiments/baseline/native/build_fixed/tuttekiri_cli). Outputs
# go to results/kill/native200_julia/ (the C++ wrote results/kill/native200/). `--limit K`
# stops after K cells.
include(joinpath(@__DIR__, "native_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const CPP_REPO = expanduser("~/Documents/kirigami-experiments")   # PORTING.md: the reference repo

# ---- cell model ------------------------------------------------------------
# A "cell" is one (graph id, variant) pair of the 200 x {sigma_mc, sigma_def,
# native} grid.  The grid is a deterministic function of (n, maxf, sigmadir):
# make_graph(id, 100, maxf, 1400) is seeded from id alone, and sigma_def exists
# only where results/kill/k5/sigma/<kind>_<id>.json is present and matches.
Base.@kwdef mutable struct Cell
    id::Int = -1
    kind::String = ""; variant::String = ""
    N::Int = 0; F::Int = 0
end

# sigma_def of one graph: the archived K5 file if `sigmadir` is given, else the frozen row
function sigma_def_of(row::PopRow, sigmadir::AbstractString, m0::K.Mesh)
    isempty(sigmadir) || return load_sigma_def(sigmadir, row.kind, row.id, m0)
    return row.sigma_def === nothing ? Int[] : row.sigma_def
end

# The frozen population as (row, m0 with topology, sigma_mc, sigma_def), walked in the
# order of the C++ `for (id = 0; id < 400 && gidx + 1 < n; ++id)` loop.
function walk_population(n::Int, maxf::Int, sigmadir::AbstractString, regenerate::Bool)
    (maxf == 800 || regenerate) || @warn "frozen native200 was built with --maxf 800; use --regenerate for maxf = $maxf"
    out = Tuple{PopRow,K.Mesh,Vector{Int},Vector{Int}}[]
    gidx = -1
    for row in load_population("native200")
        gidx + 1 < n || break
        if regenerate
            g = K.make_graph(row.id, 100, maxf, 1400)
            g.ok || continue
            row.mesh = g.mesh
        elseif !row.ok
            continue
        end
        m0 = row.mesh
        K.build_topology!(m0)
        gidx += 1
        gidx >= n && break
        push!(out, (row, m0, copy(m0.sigma), sigma_def_of(row, sigmadir, m0)))
    end
    return out
end

# The full grid, in the same order the original 600-cell run dispatched it.
function enumerate_grid(pop)
    out = Cell[]
    for (row, m0, _, sigma_def) in pop
        c = Cell(id = row.id, kind = row.kind, N = K.n_vertices(m0), F = K.n_faces(m0))
        push!(out, Cell(id = c.id, kind = c.kind, N = c.N, F = c.F, variant = "sigma_mc"))
        isempty(sigma_def) || push!(out, Cell(id = c.id, kind = c.kind, N = c.N, F = c.F, variant = "sigma_def"))
        push!(out, Cell(id = c.id, kind = c.kind, N = c.N, F = c.F, variant = "native"))
    end
    return out
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n = 200; maxf = 800; timeout_s = 600; shard = 0; nshards = 1
    outdir = joinpath(REPO, "results", "kill", "native200_julia")
    sigmadir = ""
    cli = joinpath(CPP_REPO, "baseline", "native", "build_fixed", "tuttekiri_cli")
    isfile(cli) || (cli = joinpath(REPO, "baseline", "native", "build_fixed", "tuttekiri_cli"))
    work = "/tmp/kiri_native200"
    logdir = ""
    cells_file = ""; csv_override = ""; tag = "orig"; listcells = ""
    resume_csvs = String[]
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--maxf" && i < length(args); maxf = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--sigma" && i < length(args); sigmadir = args[i+1]; i += 2
        elseif a == "--cli" && i < length(args); cli = args[i+1]; i += 2
        elseif a == "--work" && i < length(args); work = args[i+1]; i += 2
        elseif a == "--timeout" && i < length(args); timeout_s = arg_i(args[i+1]); i += 2
        elseif a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshards" && i < length(args); nshards = arg_i(args[i+1]); i += 2
        elseif a == "--logdir" && i < length(args); logdir = args[i+1]; i += 2
        # --- new in WP5 ---------------------------------------------------------
        elseif a == "--cells" && i < length(args); cells_file = args[i+1]; i += 2
        elseif a == "--resume" && i < length(args); push!(resume_csvs, args[i+1]); i += 2
        elseif a == "--tag" && i < length(args); tag = args[i+1]; i += 2
        elseif a == "--csv" && i < length(args); csv_override = args[i+1]; i += 2
        elseif a == "--listcells" && i < length(args); listcells = args[i+1]; i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    isempty(logdir) && (logdir = joinpath(outdir, "logs"))

    pop = walk_population(n, maxf, sigmadir, regenerate)
    by_id = Dict(row.id => (m0, smc, sdef) for (row, m0, smc, sdef) in pop)
    kind_of = Dict(row.id => row.kind for (row, _, _, _) in pop)

    # --listcells: dump the whole (id, variant) grid with |F|, then exit.  This is
    # the authoritative definition of the 600-cell grid; the cell lists in
    # results/kill/native200/cells_*.txt are derived from it.
    if !isempty(listcells)
        grid = enumerate_grid(pop)
        open(listcells, "w") do f
            print(f, "id,variant,kind,N,F\n")
            for c in grid
                print(f, c.id, ",", c.variant, ",", c.kind, ",", c.N, ",", c.F, "\n")
            end
        end
        println(stderr, "[native200] grid: ", length(grid), " cells -> ", listcells)
        return 0
    end

    mkpath(outdir)
    mkpath(work)
    mkpath(logdir)

    # --resume: any cell already `completed` or `crashed` in a previous CSV is done.
    done = Set{Tuple{Int,String}}()
    for rc in resume_csvs
        if !isfile(rc)
            println(stderr, "[native200] --resume: cannot read ", rc)
            continue
        end
        for (li, line) in enumerate(readlines(rc))
            li == 1 && continue   # header
            fl = split_csv(line)
            length(fl) < 7 && continue
            (fl[7] == "completed" || fl[7] == "crashed") && push!(done, (parse(Int, fl[1]), fl[3]))
        end
    end

    csv_path = !isempty(csv_override) ? csv_override :
               (nshards > 1 ? joinpath(outdir, "shard_" * string(shard) * ".csv") :
                              joinpath(outdir, "native200.csv"))
    # Append (and skip the header) when the target CSV already exists, so that a
    # killed shard can be restarted with `--resume <its own csv>` without losing the
    # rows it had already produced.
    csv_exists = isfile(csv_path) && filesize(csv_path) > 0
    csv = open(csv_path, csv_exists ? "a" : "w")
    if !csv_exists
        print(csv, "id,kind,variant,N,F,dim_null,status,native_theta_collisions,native_finite,",
              "embedding_ok,our_theta_exact,our_theta_bisect,cert_pos,cert_noovl,cert_noroot,",
              "cert_eps03,secs,run\n")
        flush(csv)
    end
    wall = Timer()
    nrun = 0

    # ---- explicit cell list: run exactly these (id, variant), in file order ----
    if !isempty(cells_file)
        want = Tuple{Int,String}[]
        if !isfile(cells_file)
            println(stderr, "[native200] cannot read --cells ", cells_file)
            return 2
        end
        for line in readlines(cells_file)
            (isempty(line) || line[1] == '#') && continue
            fl = split_csv(line)
            length(fl) < 2 && continue
            fl[1] == "id" && continue   # header
            push!(want, (parse(Int, fl[1]), fl[2]))
        end
        ncell = 0
        for (id, variant) in want
            ncell += 1
            if (id, variant) in done
                println(stderr, "[native200] skip (resume) id ", id, " ", variant)
                continue
            end
            if !haskey(by_id, id)
                println(stderr, "[native200] make_graph failed for id ", id)
                continue
            end
            m0, sigma_mc, sigma_def = by_id[id]
            sp = nothing
            if variant == "sigma_mc"
                sp = sigma_mc
            elseif variant == "sigma_def"
                if isempty(sigma_def)
                    println(stderr, "[native200] no sigma_def for id ", id, ", skipping")
                    continue
                end
                sp = sigma_def
            elseif variant != "native"
                println(stderr, "[native200] unknown variant ", variant)
                continue
            end
            nrun >= limit && break
            nrun += 1
            r, _, _ = process_cell(id, kind_of[id], m0, sp, variant, cli, work, timeout_s, logdir)
            write_row(csv, r, tag)
            println(stderr, "[native200] ", tag, " cell ", ncell, "/", length(want), " id ", id, " ",
                    variant, " F=", r.F, " -> ", r.status, " (", cpp_g(r.secs), " s), wall ",
                    cpp_g(s(wall)), " s")
        end
        close(csv)
        open(csv_path * ".done", "w") do sm
            print(sm, "wall ", cpp_g(s(wall)), " s, cells ", length(want), "\n")
        end
        return 0
    end

    # ---- original mode: the whole grid, split by graph index across shards -----
    grid = enumerate_grid(pop)
    gidx = -1
    gidx_of_id = Dict{Int,Int}()
    for c in grid
        haskey(gidx_of_id, c.id) || (gidx_of_id[c.id] = (gidx += 1))
    end
    for c in grid
        (nshards > 1 && gidx_of_id[c.id] % nshards != shard) && continue
        (c.id, c.variant) in done && continue
        m0, sigma_mc, sigma_def = by_id[c.id]
        sp = c.variant == "sigma_mc" ? sigma_mc : (c.variant == "sigma_def" ? sigma_def : nothing)
        nrun >= limit && break
        nrun += 1
        r, _, _ = process_cell(c.id, c.kind, m0, sp, c.variant, cli, work, timeout_s, logdir)
        write_row(csv, r, tag)
        println(stderr, "[native200] shard ", shard, " id ", c.id, " ", c.variant, " -> ", r.status,
                ", wall ", cpp_g(s(wall)), " s")
    end
    close(csv)
    open(joinpath(outdir, "shard_" * string(shard) * "_done.txt"), "w") do sm
        print(sm, "wall ", cpp_g(s(wall)), " s\n")
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
