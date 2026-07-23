# corpus_common.jl -- the part of code/apps/kill_common.hpp the scratch checks need, read
# from the frozen populations in data/corpus/ (PORTING.md: apps never regenerate them).
#
# `make_graph(id, min_faces, max_faces, n_cap)` returns the frozen graph the C++ call of the
# same arguments produced (kind, ok, mesh with sigma_mc as `mesh.sigma`, plus sigma_def when
# the corpus has it); `reference_cases()` the 8 Phase-2 tilings with their sigma.
# TODO(julia-port): replace by Kirigami/apps/kill_common.jl once the apps land.

import JSON

const _CORPUS_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "corpus"))

# (min_faces, max_faces, n_cap) -> corpus file, from data/corpus/README.md
const _POPULATION_FILES = Dict(
    (18, 46, 220)       => "derivation_l1_140.json",   # ids 0..139 (0..79 == derivation_l1_small.json)
    (100, 800, 1400)    => "k1a_200.json",
    (100, 5000, 1 << 30) => "k3a_500.json",
)

const _population_cache = Dict{String,Vector{Any}}()

function _population(file::String)
    get!(_population_cache, file) do
        JSON.parsefile(joinpath(_CORPUS_DIR, file))
    end
end

struct Graph
    mesh::Kirigami.Mesh
    kind::String
    id::Int
    ok::Bool
    sigma_def::Union{Nothing,Vector{Int}}
end

"""Frozen `kill::make_graph(id, min_faces, max_faces, n_cap)`; `ok == false` (empty mesh)
when the id is outside the frozen range or the C++ call failed."""
function make_graph(id::Int, min_faces::Int, max_faces::Int, n_cap::Int)
    file = get(_POPULATION_FILES, (min_faces, max_faces, n_cap), nothing)
    file === nothing && error("make_graph($id, $min_faces, $max_faces, $n_cap) is not frozen in data/corpus/")
    pop = _population(file)
    (id < 0 || id >= length(pop)) && return Graph(Kirigami.Mesh(), "", id, false, nothing)
    row = pop[id + 1]                      # ids are 0-based in the corpus, ascending
    row["id"] == id || error("corpus $file is not id-ordered at $id")
    row["ok"] || return Graph(Kirigami.Mesh(), String(row["kind"]), id, false, nothing)
    m = Kirigami.mesh_from_json(row["mesh"])   # `orientation` == sigma_mc
    sd = get(row, "sigma_def", nothing)
    sd = sd === nothing ? nothing : [Int(s) for s in sd]
    return Graph(m, String(row["kind"]), id, true, sd)
end

struct RefCase
    name::String
    mesh::Kirigami.Mesh
    periodic::Bool
end

"""Frozen `kill::reference_cases()`: the 8 Phase-2 tilings, sigma already in `mesh.sigma`."""
function reference_cases()
    out = RefCase[]
    for row in _population("reference_cases_8.json")
        m = Kirigami.mesh_from_json(row["mesh"])
        if isempty(m.sigma) && haskey(row, "sigma") && row["sigma"] !== nothing
            m.sigma = [Int(s) for s in row["sigma"]]
            Kirigami.build_topology!(m)
        end
        push!(out, RefCase(String(row["name"]), m, Bool(row["periodic"])))
    end
    return out
end
