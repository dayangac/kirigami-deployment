# apps/kill_common.jl -- the graph populations shared by every Phase-5 kill experiment and
# the app-side helpers around them. Port of code/apps/kill_common.hpp.
#
# The population DEFINITIONS (`Graph`, `checkerboard`, `make_graph`, `RefCase`,
# `reference_cases`) live in the package (src/core/kill_common.jl) because the tests need
# them; this file adds what only the apps use: the frozen-corpus loaders, the shape-space
# cache, the deployable population, formatting and a timer.
#
# FROZEN CORPUS FIRST. Every population the C++ built procedurally was written once to
# data/corpus/<name>.json (see data/corpus/README.md). `population(name)` loads that file;
# `population(name; regenerate = true)` rebuilds it through the bit-exact generators
# instead (`make_graph` reproduces all 200 K1a graphs vertex-for-vertex, verified
# 2026-09-19). Apps accept `--regenerate` and default to the frozen file, so every
# experiment runs on identical inputs whichever machine it runs on.

include(joinpath(@__DIR__, "common_app.jl"))

const CORPUS_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "corpus"))

# ---- frozen populations ---------------------------------------------------------------

"""One row of a frozen population file, with the mesh rebuilt (topology included) and
`mesh.sigma = sigma_mc`. Rows with `ok == false` carry an empty mesh."""
mutable struct PopRow
    id::Int
    kind::String
    ok::Bool
    gidx::Int                      # population index of the sharded apps (-1 if absent)
    mesh::K.Mesh                   # with sigma = sigma_mc
    sigma_mc::Vector{Int}
    sigma_checker::Vector{Int}
    sigma_def::Union{Vector{Int},Nothing}
    make_graph_args::Vector{Int}   # [min_faces, max_faces, n_cap]
    seed::UInt32
end

_intvec(x) = x === nothing ? nothing : Int[Int(v) for v in x]

function _pop_row(r::AbstractDict)
    ok = Bool(r["ok"])
    mesh = ok ? K.mesh_from_json(r["mesh"]) : K.Mesh()
    sigma_mc = _intvec(get(r, "sigma_mc", Int[]))
    ok && (mesh.sigma = copy(sigma_mc))
    args = haskey(r, "make_graph_args") ? _intvec(r["make_graph_args"]) : Int[]
    return PopRow(Int(r["id"]), String(r["kind"]), ok, Int(get(r, "gidx", -1)), mesh,
                  sigma_mc, _intvec(get(r, "sigma_checker", Int[])),
                  _intvec(get(r, "sigma_def", nothing)), args,
                  UInt32(get(r, "seed", 0) % UInt32))
end

"""
    load_population(name) -> Vector{PopRow}

Loads `data/corpus/<name>.json` (`name` may also be a path). Rows are in file order, i.e.
ascending id; for the k1a / native200 / k3a / k2c / regime / yield / e1 / k2b /
derivation_l1_small files that is `rows[id - first_id + 1]`.
"""
function load_population(name::AbstractString; first::Union{Nothing,Int} = nothing)
    path = isfile(name) ? String(name) : joinpath(CORPUS_DIR, name * ".json")
    isfile(path) || error("load_population: no corpus file $path")
    rows = first === nothing ? JSON.parsefile(path) : corpus_first_rows(path, first)
    return PopRow[_pop_row(r) for r in rows]
end

# Reads the first `n` rows of a top-level JSON array without parsing the whole (up to
# 165 MB) file: a byte scanner that is string- and escape-aware.
function corpus_first_rows(path::AbstractString, n::Int)
    s = read(path)
    depth = 0
    instr = false
    esc = false
    rows = 0
    stop = 0
    for (i, b) in enumerate(s)
        if instr
            if esc
                esc = false
            elseif b == UInt8('\\')
                esc = true
            elseif b == UInt8('"')
                instr = false
            end
            continue
        end
        if b == UInt8('"')
            instr = true
        elseif b == UInt8('[') || b == UInt8('{')
            depth += 1
        elseif b == UInt8(']') || b == UInt8('}')
            depth -= 1
            if depth == 1
                rows += 1
                if rows == n
                    stop = i
                    break
                end
            end
        end
    end
    stop == 0 && return JSON.parse(String(s))
    return JSON.parse(String(s[1:stop]) * "]")
end

"""`population_row(rows, id)`: the row with graph id `id` (error if absent)."""
function population_row(rows::Vector{PopRow}, id::Int)
    i = findfirst(r -> r.id == id, rows)
    i === nothing && error("population_row: id $id is not in the population")
    return rows[i]
end

"""
    population(name; ids = nothing, regenerate = false) -> Vector{K.Graph}

The kill population `name` as `K.Graph`s (mesh with sigma_mc, kind, id, ok). By default the
frozen corpus file is read; with `regenerate = true` every graph is rebuilt through
`K.make_graph(id, make_graph_args...)`, the args being taken from the corpus file so the
two paths cannot drift apart. `ids` restricts to those graph ids (default: all rows).
"""
function population(name::AbstractString; ids = nothing, regenerate::Bool = false,
                    first::Union{Nothing,Int} = nothing)
    rows = load_population(name; first = first)
    ids === nothing || (rows = PopRow[population_row(rows, id) for id in ids])
    out = K.Graph[]
    for r in rows
        if regenerate
            isempty(r.make_graph_args) && error("population: $name carries no make_graph_args, cannot regenerate")
            push!(out, K.make_graph(r.id, r.make_graph_args...))
        else
            push!(out, K.Graph(r.mesh, r.kind, r.id, r.ok))
        end
    end
    return out
end

"""
    reference_cases(; regenerate = false) -> Vector{K.RefCase}

The eight Phase-2 reference cases, from `data/corpus/reference_cases_8.json` or rebuilt by
`K.reference_cases()`.
"""
function reference_cases(; regenerate::Bool = false)
    regenerate && return K.reference_cases()
    path = joinpath(CORPUS_DIR, "reference_cases_8.json")
    isfile(path) || error("reference_cases: no corpus file $path")
    cs = K.RefCase[]
    for r in JSON.parsefile(path)
        m = K.mesh_from_json(r["mesh"])
        m.sigma = _intvec(r["sigma"])
        push!(cs, K.RefCase(String(r["name"]), m, Bool(r["periodic"])))
    end
    return cs
end

"""Removes `--regenerate` from `args` and returns (args_without_it, regenerate::Bool)."""
function take_regenerate_flag(args::Vector{String})
    keep = String[a for a in args if a != "--regenerate"]
    return keep, length(keep) != length(args)
end

# ---- per-graph helpers -----------------------------------------------------------------

"""Number of faces of `m` (under positions X) with signed area <= tol * bbox area."""
count_inverted(m::K.Mesh, X::Vector{Vec2}, rel_tol::Float64 = 1e-12) =
    K.count_inverted_faces(m, X, rel_tol)

median_edge_length(m::K.Mesh) = K.median_edge_length(m)

# ---- shape-space cache -----------------------------------------------------------------
# The dense SVD of [L;B] costs seconds at N ~ 1400 and four kill experiments need the same
# X0 and Phi, so it is cached to disk (raw doubles, column-major; the C++ layout, so files
# written by either implementation are read by the other).
mutable struct Shape
    N::Int
    k::Int
    X0::Matrix{Float64}    # N x 2
    Phi::Matrix{Float64}   # N x k
    rank_L::Int
    H::Int
    n_interior::Int
    ok::Bool
end
Shape() = Shape(0, 0, zeros(0, 2), zeros(0, 0), 0, 0, 0, false)

shape_path(dir::AbstractString, id::Int) = joinpath(dir, "shape_$(id).bin")

function load_shape(dir::AbstractString, id::Int)
    p = shape_path(dir, id)
    isfile(p) || return nothing
    s = Shape()
    open(p, "r") do f
        hdr = Vector{Int32}(undef, 5)
        read!(f, hdr)
        s.N, s.k, s.rank_L, s.H, s.n_interior = Int.(hdr)
        s.X0 = Matrix{Float64}(undef, s.N, 2)
        read!(f, s.X0)
        s.Phi = Matrix{Float64}(undef, s.N, s.k)
        s.k > 0 && read!(f, s.Phi)
        s.ok = true
    end
    return s
end

function save_shape(dir::AbstractString, id::Int, s::Shape)
    mkpath(dir)
    open(shape_path(dir, id), "w") do f
        write(f, Int32[s.N, s.k, s.rank_L, s.H, s.n_interior])
        write(f, s.X0)
        s.k > 0 && write(f, s.Phi)
    end
    return nothing
end

"""Solves (or loads from `cache_dir`) the Tutte auxetic system for a graph, fixed boundary."""
function shape_space(m::K.Mesh, c::K.CutStructure, hs::K.HoleSet, cache_dir::AbstractString, id::Int)
    if id >= 0 && !isempty(cache_dir)
        s = load_shape(cache_dir, id)
        s !== nothing && s.N == K.n_vertices(m) && return s
    end
    sys = K.assemble_system(c, hs, m.X, K.Fixed)
    sr = K.solve_system(sys, m.X)
    s = Shape(K.n_vertices(m), sr.dim_null, sr.X0, sr.Phi, sr.rank_L, sr.H,
              K.n_interior_vertices(m), sr.projection_ok)
    id >= 0 && !isempty(cache_dir) && s.ok && save_shape(cache_dir, id, s)
    return s
end

# ---- the deployable population -----------------------------------------------------------
# Configurations that are BOTH uniformly deployable (Eq. (2) holds) AND embedded in the flat
# state: the authored tilings at several clip radii, each with the Eq. (6) projection plus
# a few shape-space samples that are still embedded. `deployable_population()` reads the
# frozen file (solver OUTPUT `X`, matched to ~1e-9 rather than bit-exactly);
# `deployable_population(regenerate = true)` rebuilds it as kill_common.hpp does.
mutable struct DeployableConfig
    name::String
    family::String
    mesh::K.Mesh
    X::Vector{Vec2}
    dim_null::Int
    sample::Int   # -1 = X0 / X_ini, else the sample index
end

function deployable_population(; regenerate::Bool = false, samples_per_base::Int = 8,
                               radius_frac::Float64 = 0.2)
    regenerate && return build_deployable_population(samples_per_base, radius_frac)
    # the (8, 0.2) default is deployable_population.json; the other variants the C++ apps
    # used are frozen as deployable_population_<samples>_<radius>.json (k8a: (2, 0.2);
    # e1: (20, 0.2) and (20, 0.35)), with the C++ Eigen-basis samples
    path = (samples_per_base == 8 && radius_frac == 0.2) ?
           joinpath(CORPUS_DIR, "deployable_population.json") :
           joinpath(CORPUS_DIR, "deployable_population_$(samples_per_base)_$(radius_frac).json")
    isfile(path) || error("deployable_population: no corpus file $path (use regenerate = true)")
    out = DeployableConfig[]
    for r in JSON.parsefile(path)
        m = K.mesh_from_json(r["mesh"])
        X = Vec2[Vec2(Float64(p[1]), Float64(p[2])) for p in r["X"]]
        push!(out, DeployableConfig(String(r["name"]), String(r["family"]), m, X,
                                    Int(r["dim_null"]), Int(r["sample"])))
    end
    return out
end

function build_deployable_population(samples_per_base::Int, radius_frac::Float64)
    out = DeployableConfig[]
    radii = [2.0, 2.5, 3.0, 3.5, 4.0]
    rng0 = K.MT19937(20260903)
    for (ri, R) in enumerate(radii)
        ri0 = ri - 1   # the C++ index
        fams = [
            ("squares", K.tiling_squares(K.rect(Vec2(0.5 * R, 0.5 * R), R / 2 + 0.01, R / 2 + 0.01)), true),
            ("triangles", K.tiling_triangles(K.disk(Vec2(0.13, 0.07), R)), true),
            ("kagome", K.tiling_kagome(K.disk(Vec2(0.13, 0.07), R)), true),
            ("hexagons", K.tiling_hexagons(K.disk(Vec2(0.13, 0.07), R)), false),
            ("trunc_square", K.tiling_truncated_square(K.disk(Vec2(0.13, 0.07), R + 1.0)), false),
            ("snub_square", K.tiling_snub_square(K.disk(Vec2(0.13, 0.07), R)), false),
            ("t3_4_3_12", K.tiling_3_4_3_12(K.disk(Vec2(0.13, 0.07), R + 1.2)), false),
        ]
        for (fname, m, checker) in fams
            K.n_faces(m) < 4 && continue
            K.build_topology!(m)
            rng = K.MT19937(UInt32(777) + UInt32(13) * UInt32(ri0) + UInt32(101))
            m.sigma = checker ? K.checkerboard(m) : K.assign_orientation_relaxation(m, rng, 8, 500, 180).sigma
            isempty(m.sigma) && continue
            K.build_topology!(m)
            c = K.make_cut(m)
            hs = K.holes_partition(c)
            X = copy(m.X)
            Phi = zeros(0, 0)
            k = 0
            if !K.deployable(K.hole_residuals(c, X, hs), 1e-9)
                sys = K.assemble_system(c, hs, m.X, K.Fixed)
                sr = K.solve_system(sys, m.X)
                sr.projection_ok || continue
                X = K.matrix_to_points(sr.X0)
                Phi = sr.Phi
                k = sr.dim_null
            else
                sys = K.assemble_system(c, hs, m.X, K.Fixed)
                sr = K.solve_system(sys, m.X)
                if sr.projection_ok
                    Phi = sr.Phi
                    k = sr.dim_null
                end
            end
            K.deployable(K.hole_residuals(c, X, hs), 1e-8) || continue
            K.has_collision(c, K.deploy(c, X, 0.0).Y) && continue
            base = fname * "_R" * string(trunc(Int, 10 * R))
            push!(out, DeployableConfig(base, fname, m, X, k, -1))
            k <= 0 && continue
            med = median_edge_length(m)
            G = K.NormalDist(0.0, 1.0)   # fresh per base, as the C++ declares it (cached second draw)
            kept = 0
            for _ in 1:60
                kept >= samples_per_base && break
                T = Matrix{Float64}(undef, k, 2)
                for i in eachindex(T)   # column-major, as Eigen's data() order
                    T[i] = K.normal(G, rng0)
                end
                D = Phi * T
                mx = maximum(abs, D)
                mx <= 0 && break
                D .*= radius_frac * med / mx
                Xs = K.matrix_to_points(K.points_to_matrix(X) + D)
                count_inverted(m, Xs) != 0 && continue
                K.has_collision(c, K.deploy(c, Xs, 0.0).Y) && continue
                push!(out, DeployableConfig(base * "_s" * string(kept), fname, m, Xs, k, kept))
                kept += 1
            end
        end
    end
    return out
end

# ---- timer -------------------------------------------------------------------------------
struct Timer
    t0::UInt64
end
Timer() = Timer(time_ns())
"""Elapsed seconds, at millisecond resolution like the C++."""
s(t::Timer) = floor((time_ns() - t.t0) / 1e6) / 1000.0
