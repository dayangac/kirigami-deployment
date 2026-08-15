# gui/backend.jl -- the live compute backend of the explorer: pattern -> orientation ->
# design -> characterize -> a GuiModel.Design with the exact deployment basis.
#
# Everything here is window-free. `ComputeBackend` is the phase-2 implementation of
# GuiModel.AbstractBackend; `run_design` is a pure function of (mesh, sigma, options), and
# `DesignJob` wraps it in a task with a progress observable and a cancel flag for app.jl.
# The same calls as apps/kiri_design.jl are made (design_range_max / design_baseline /
# design_constrained, orientation_maxcut / orientation_defect, characterize), so a design
# computed here reproduces the CLI's numbers for the same seed.
module GuiBackend

using Kirigami
using LinearAlgebra: norm
import JSON
using ..GuiModel: GuiModel, Design, face_sense
const K = Kirigami
const Vec2 = K.Vec2

# ---- pattern files ---------------------------------------------------------------------

"""One entry of data/web_patterns.json (the web app's UKP pattern library)."""
struct PatternEntry
    name::String
    kind::String
    text::String
    has_sigma::Bool
    periodic::Bool
end

function load_patterns(path::AbstractString)
    return [PatternEntry(string(j["name"]), string(get(j, "kind", "ukp")), string(j["text"]),
                         Bool(get(j, "has_sigma", false)), Bool(get(j, "periodic", false)))
            for j in JSON.parsefile(path)]
end

"""The .UKP format: an .OBJ of the 2D unit pattern plus `px x y` / `py x y` (the unit
parallelogram; [0 0] = not periodic), `fc 0|1` per face and a `rep2x2` hint. A plain .OBJ
is the same parser with none of the extras. Port of web/bindings.cpp parse_pattern_text."""
mutable struct ParsedPattern
    V::Vector{Vec2}
    F::Vector{Vector{Int}}     # 1-based into V
    fc::Vector{Int}            # 0/1 per face, empty when the file has none
    px::Vec2
    py::Vec2
    periodic::Bool
    rep2x2::Bool
end

function parse_pattern_text(text::AbstractString)
    p = ParsedPattern(Vec2[], Vector{Int}[], Int[], Vec2(0, 0), Vec2(0, 0), false, false)
    for line in split(text, '\n')
        line = rstrip(line, '\r')
        toks = split(line)
        isempty(toks) && continue
        tag = toks[1]
        num(i) = length(toks) >= i ? parse(Float64, toks[i]) : 0.0
        if tag == "v"
            push!(p.V, Vec2(num(2), num(3)))          # a z, when present, is ignored
        elseif tag == "f"
            f = Int[]
            for tok in toks[2:end]
                tok = first(split(tok, '/'))          # v/vt/vn slashes
                isempty(tok) && continue
                idx = parse(Int, tok)
                idx = idx < 0 ? length(p.V) + idx + 1 : idx   # OBJ relative index -> 1-based
                push!(f, idx)
            end
            length(f) >= 3 && push!(p.F, f)
        elseif tag == "px" || tag == "py"
            v = Vec2(num(2), num(3))
            tag == "px" ? (p.px = v) : (p.py = v)
        elseif tag == "fc"
            push!(p.fc, length(toks) >= 2 ? parse(Int, toks[2]) : 0)
        elseif tag == "rep2x2"
            p.rep2x2 = true
        end
    end
    p.periodic = sum(abs2, p.px) > 1e-18 && sum(abs2, p.py) > 1e-18
    for f in p.F, i in f
        (1 <= i <= length(p.V)) || throw(ArgumentError("face index out of range"))
    end
    (isempty(p.V) || isempty(p.F)) && throw(ArgumentError("no vertices or faces in this file"))
    return p
end

# fc 0/1 is the app's two-colouring; sigma is +-1 with the same meaning (fc 0 -> +1).
sigma_of_fc(c::Int) = c != 0 ? -1 : 1

"""Tile the unit pattern rep_x x rep_y times (periodic patterns only), weld the polygon
soup into a Mesh at 1e-6 of the pattern scale and return (mesh, sigma_from_file) where
sigma is `nothing` unless the file colours every face."""
function pattern_mesh(p::ParsedPattern, rep_x::Int = 2, rep_y::Int = 2)
    polys = Vector{Vec2}[]
    fc = Int[]
    per = p.periodic && (rep_x > 1 || rep_y > 1)
    nx = per ? rep_x : 1; ny = per ? rep_y : 1
    for a in 0:nx-1, b in 0:ny-1
        off = a * p.px + b * p.py
        for (fi, f) in enumerate(p.F)
            push!(polys, [p.V[i] + off for i in f])
            push!(fc, fi <= length(p.fc) ? p.fc[fi] : -1)
        end
    end
    m = K.mesh_from_polygons(polys, 1e-6)
    K.build_topology!(m)
    sigma = nothing
    if K.n_faces(m) == length(fc) && all(c -> c >= 0, fc)
        sigma = [sigma_of_fc(c) for c in fc]
    end
    return m, sigma
end

"""A graph of a frozen corpus file (data/corpus/<name>.json), with its `sigma_mc`."""
function load_corpus_graph(path::AbstractString, id::Int)
    for r in JSON.parsefile(path)
        Int(r["id"]) == id || continue
        Bool(r["ok"]) || error("corpus graph $id is marked ok = false")
        m = K.mesh_from_json(r["mesh"])
        m.sigma = Int[Int(v) for v in r["sigma_mc"]]
        K.build_topology!(m)
        return m
    end
    error("corpus graph $id not found in $path")
end

# ---- orientation -------------------------------------------------------------------------

"""Orientation by rule: "file" (the pattern's own fc colours, else an error), "mc" (Eq. (1)
max-cut, as kiri_design --sigma mc) or "def" (K5's defect-minimising search from mc)."""
function choose_sigma(m::K.Mesh, rule::AbstractString, file_sigma; seed::Integer = 9000,
                      defect_cap_mult::Int = 20)
    if rule == "file"
        file_sigma === nothing && throw(ArgumentError("this pattern carries no orientation; choose mc or def"))
        return Vector{Int}(file_sigma)
    end
    mc = K.orientation_maxcut(m, seed)
    length(mc) == K.n_faces(m) || error("Eq. (1) orientation failed")
    rule == "mc" && return mc
    if rule == "def"
        dr = K.orientation_defect(m, mc, defect_cap_mult * K.n_faces(m), seed)
        dr.ok || error("the defect search found no sigma with c(Gamma) = 1")
        return dr.sigma
    end
    throw(ArgumentError("unknown sigma rule $rule"))
end

# ---- design --------------------------------------------------------------------------------

const METHODS = ("range_max", "baseline", "constrained")

"""Run the method on (mesh, sigma) from X_ini = mesh.X, exactly as kiri_design.jl does:
range_max = design_range_max (the paper's method, K9c), baseline = design_baseline (the
2026 Eq. (6) projection alone), constrained = design_constrained (K9 variant (b)).
`progress(msg)` is called at stage boundaries. Returns the DesignResult."""
function run_design(m::K.Mesh, sigma::Vector{Int}; method::AbstractString = "range_max",
                    seed::Integer = 9000, progress::Function = _ -> nothing)
    method in METHODS || throw(ArgumentError("unknown method $method"))
    ms = K.Mesh(copy(m.X), deepcopy(m.faces))
    ms.sigma = copy(sigma)
    ms.periodic = m.periodic
    K.build_topology!(ms)
    progress("designing ($(method), F = $(K.n_faces(ms)))")
    if method == "range_max"
        rmo = K.RangeMaxOptions()
        rmo.seed = UInt32(seed)
        return K.design_range_max(ms, sigma, ms.X, rmo).design
    end
    opt = K.DesignOptions()
    opt.seed = UInt32(seed)
    return method == "baseline" ? K.design_baseline(ms, sigma, ms.X, opt) :
                                  K.design_constrained(ms, sigma, ms.X, opt)
end

"""Where a computed design came from, kept on the Design so exports can rebuild the
real cut layout (mesh in M coordinates, sigma, flat embedding X)."""
struct DesignSource
    mesh::K.Mesh
    sigma::Vector{Int}
    X::Vector{Vec2}
    method::String
    ch::K.Characterization
end

"""Build a GuiModel.Design from a flat embedding X of (mesh, sigma). The frames live on the
M' (prime) vertices: v0 = basis_eval(B, 0) is the deployed flat state, v1 the frame at
theta_half = pi/2 from the exact basis Y(theta) = cos(theta/2) C + sin(theta/2) S, and
faces are c.prime_faces. theta_max/eps_max come from `characterize` (opt: K9's eps 0.3)."""
function design_from_embedding(name::AbstractString, kind::AbstractString, sigma_label::AbstractString,
                               m::K.Mesh, sigma::Vector{Int}, X::Vector{Vec2}, method::AbstractString;
                               ch::Union{K.Characterization,Nothing} = nothing)
    ms = K.Mesh(copy(X), deepcopy(m.faces))
    ms.sigma = copy(sigma)
    ms.periodic = m.periodic
    K.build_topology!(ms)
    c = K.make_cut(ms)
    ch === nothing && (ch = K.characterize(ms, sigma, X))
    B = K.deploy_basis(c, X)
    theta_half = pi / 2            # any angle with sin(theta_half/2) != 0; the basis is exact
    C = Matrix{Float64}(B.C); S = Matrix{Float64}(B.S)
    faces = [copy(f) for f in c.prime_faces]
    return Design(String(name), String(kind), String(sigma_label), face_sense(C, S, faces),
                  K.n_faces(ms), K.n_split(c), ch.theta_max, ch.eps_max, theta_half, faces, C, S,
                  DesignSource(ms, copy(sigma), copy(X), String(method), ch))
end

# ---- the backend --------------------------------------------------------------------------

"""Phase-2 backend: designs are computed live from patterns and corpus graphs. `designs`
returns what has been computed so far (the app adds to it as jobs finish)."""
mutable struct ComputeBackend <: GuiModel.AbstractBackend
    patterns::Vector{PatternEntry}
    corpus_path::String
    designs::Vector{Design}
end
ComputeBackend(patterns_path::AbstractString, corpus_path::AbstractString = "") =
    ComputeBackend(isfile(patterns_path) ? load_patterns(patterns_path) : PatternEntry[],
                   String(corpus_path), Design[])
GuiModel.designs(b::ComputeBackend) = (b.designs, 0)

"""Full pipeline for one request; returns the Design. `progress(msg)` at each stage."""
function design_pattern(b::ComputeBackend, p::PatternEntry; rep_x::Int = 2, rep_y::Int = 2,
                        sigma_rule::AbstractString = "file", method::AbstractString = "range_max",
                        seed::Integer = 9000, progress::Function = _ -> nothing)
    progress("parsing " * p.name)
    parsed = parse_pattern_text(p.text)
    m, file_sigma = pattern_mesh(parsed, rep_x, rep_y)
    progress("orientation ($sigma_rule)")
    sigma = choose_sigma(m, sigma_rule, file_sigma; seed = seed)
    return design_mesh(p.name * "_$(rep_x)x$(rep_y)", "ukp", m, sigma, sigma_rule, method, seed, progress)
end

function design_corpus(b::ComputeBackend, id::Int; sigma_rule::AbstractString = "mc",
                       method::AbstractString = "range_max", seed::Integer = 9000,
                       progress::Function = _ -> nothing)
    progress("loading corpus graph $id")
    m = load_corpus_graph(b.corpus_path, id)
    # "file" here means the frozen sigma_mc of the corpus row
    sigma = sigma_rule == "file" ? copy(m.sigma) : choose_sigma(m, sigma_rule, m.sigma; seed = seed)
    return design_mesh("k1a_$id", "corpus", m, sigma, sigma_rule, method, seed, progress)
end

function design_mesh(name, kind, m::K.Mesh, sigma::Vector{Int}, sigma_rule, method, seed, progress)
    d = run_design(m, sigma; method = method, seed = seed, progress = progress)
    d.ok || error("design failed: " * d.status)
    progress("characterizing")
    return design_from_embedding(name, kind, sigma_rule, m, sigma, d.X, method; ch = d.ch)
end

"""Re-characterize a computed design at its own embedding (same call as the CLI)."""
function recharacterize(d::Design)
    d.source isa DesignSource || throw(ArgumentError("this design has no mesh; it was loaded from a frame file"))
    s = d.source
    return K.characterize(s.mesh, s.sigma, s.X)
end

# ---- exports -----------------------------------------------------------------------------

"""Laser-cut SVG of a computed design at angle theta via the export layer (felt/laser
profile, average edge = `avg_edge_mm`, hinge neck = `hinge_mm`), as the web app did."""
function export_svg_string(d::Design, theta::Real; avg_edge_mm::Real = 40.0, hinge_mm::Real = 1.5)
    L = _layout(d, theta, avg_edge_mm, hinge_mm)
    return K.svg_string(L)
end

function _layout(d::Design, theta, avg_edge_mm, hinge_mm)
    d.source isa DesignSource || throw(ArgumentError("this design has no mesh; export needs a computed design"))
    s = d.source
    c = K.make_cut(s.mesh)
    L = [norm(s.X[e.key.a] - s.X[e.key.b]) for e in s.mesh.edges]
    avg = isempty(L) ? 1.0 : sum(L) / length(L)
    scale = avg > 1e-12 ? avg_edge_mm / avg : 1.0
    prof = K.profile_felt_laser()
    hinge_mm > 0 && (prof.neck_width = Float64(hinge_mm))
    return K.build_layout(s.mesh, c, s.X, Float64(theta), scale, prof)
end

"""3MF of the printed solid at angle theta (build_layout -> build_solid -> write_3mf)."""
function export_3mf(d::Design, theta::Real, path::AbstractString; avg_edge_mm::Real = 40.0, hinge_mm::Real = 1.5)
    L = _layout(d, theta, avg_edge_mm, hinge_mm)
    solid, warnings = K.build_solid(L)
    K.write_3mf(solid, path, "kiri " * L.profile.name)
    return warnings
end

# ---- async job -----------------------------------------------------------------------------

"""A design request running off the UI. State is plain and lock-protected -- never an
Observable, because GLMakie must only be touched from the render thread; app.jl polls
`snapshot(job)` on a timer. The solver has no interruption points, so `cancel!` lets the
current stage finish and then discards the result. Runs on a thread when Julia has more
than one (`julia -t auto`), else as a task that yields only at stage boundaries."""
mutable struct DesignJob
    lock::ReentrantLock
    progress::String
    result::Any                   # nothing | Design | Exception
    done::Bool
    cancelled::Bool
    task::Union{Task,Nothing}
end

function start_job(f::Function)
    job = DesignJob(ReentrantLock(), "starting", nothing, false, false, nothing)
    prog(msg) = (lock(job.lock) do; job.progress = msg; end; yield())
    body = () -> begin
        r = try
            d = f(prog)
            lock(job.lock) do; job.cancelled; end ? ErrorException("cancelled") : d
        catch e
            e
        end
        lock(job.lock) do
            job.result = r; job.done = true
        end
    end
    job.task = Threads.nthreads() > 1 ? Threads.@spawn(body()) : @async(body())
    return job
end

"""(progress, done, result) as of now."""
snapshot(job::DesignJob) = lock(job.lock) do
    (job.progress, job.done, job.result)
end

function cancel!(job::DesignJob)
    lock(job.lock) do
        job.cancelled = true
        job.progress = "cancelling after the current stage"
    end
    return nothing
end

end # module
