# native200_cell -- write the exact input JSON one Native200 cell feeds to the
# authors' `tuttekiri_cli prevent`, so a crashing cell can be reproduced by hand
# (under lldb, or against a patched binary) without running the whole driver.
# Port of code/apps/native200_cell.cpp.
#
#   julia --project=Kirigami Kirigami/apps/native200_cell.jl --id 14 --variant sigma_def
#         --out /tmp/cell14 [--maxf 800] [--sigma DIR] [--regenerate]
#
# writes <out>/input.json (for sigma_mc / sigma_def: the mesh with that
# orientation field; for `native`: <out>/raw.json without an orientation field --
# the caller must then run `tuttekiri_cli color raw.json --out input.json`
# exactly as kill_native200.jl does) and prints id/kind/N/F. The graph and sigma_def
# come from the frozen data/corpus/native200.json (`--sigma DIR` reads the archived K5
# file instead; `--regenerate` rebuilds the graph through make_graph).
include(joinpath(@__DIR__, "kill_common.jl"))

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    id = -1; maxf = 800
    variant = "sigma_mc"; out = "/tmp/native200_cell"
    sigmadir = ""
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--id" && i < length(args); id = arg_i(args[i+1]); i += 2
        elseif a == "--variant" && i < length(args); variant = args[i+1]; i += 2
        elseif a == "--out" && i < length(args); out = args[i+1]; i += 2
        elseif a == "--maxf" && i < length(args); maxf = arg_i(args[i+1]); i += 2
        elseif a == "--sigma" && i < length(args); sigmadir = args[i+1]; i += 2
        else; i += 1
        end
    end
    if id < 0
        println(stderr, "need --id")
        return 2
    end
    rows = load_population("native200")
    ri = findfirst(r -> r.id == id, rows)
    m0 = K.Mesh(); kind = ""; frozen_def = nothing
    if regenerate || ri === nothing
        g = K.make_graph(id, 100, maxf, 1400)
        if !g.ok
            println(stderr, "make_graph failed for id ", id)
            return 2
        end
        m0 = g.mesh; kind = g.kind
    else
        row = rows[ri]
        if !row.ok
            println(stderr, "make_graph failed for id ", id)
            return 2
        end
        m0 = row.mesh; kind = row.kind; frozen_def = row.sigma_def
    end
    K.build_topology!(m0)
    mkpath(out)

    if variant == "native"
        raw = K.Mesh(m0.X, m0.faces)
        K.save_mesh_json(raw, joinpath(out, "raw.json"))
    else
        m = K.Mesh(m0.X, m0.faces)
        m.sigma = copy(m0.sigma)
        if variant == "sigma_def"
            sd = Int[]
            if !isempty(sigmadir)
                p = joinpath(sigmadir, kind * "_" * string(id) * ".json")
                isfile(p) && (sd = K.load_mesh_json(p).sigma)
            elseif frozen_def !== nothing
                sd = frozen_def
            end
            if length(sd) != K.n_faces(m0)
                println(stderr, "no usable sigma_def for id ", id, " (sigmadir '", sigmadir, "')")
                return 2
            end
            m.sigma = copy(sd)
        end
        K.save_mesh_json(m, joinpath(out, "input.json"))
    end
    println("id=", id, " kind=", kind, " variant=", variant, " N=", K.n_vertices(m0), " F=",
            K.n_faces(m0), " out=", out)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(copy(ARGS)))
end
