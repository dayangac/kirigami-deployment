# gui/model.jl -- pure (window-free) state and computations behind the desktop explorer.
#
# Phase 1 reads pre-characterised designs from data/deploy_frames.json (data/web_patterns.json
# is the raw UKP pattern library for phase 2's design flow).  Everything the UI
# shows is a function of a `Design` and an angle θ, so `app.jl` only wires Observables to
# the functions here.  Phase 2 will swap the `load_designs` file backend for live
# Kirigami.method calls without touching the UI (see `AbstractBackend`).
module GuiModel

import JSON

# ---- palette (web explorer --face-a / --face-b / --accent / --bad / --good, light theme)
const FACE_A = "#5F82A8"
const FACE_B = "#B9C9DA"
const FACE_EDGE = "#2C3E52"
const ACCENT = "#C25A21"
const ACCENT_SOFT = "#F1D9CC"
const BAD = "#B23A3A"
const BAD_SOFT = "#F4D9D9"
const GOOD = "#2E7D5B"
const BG = "#F3F4F1"
const PANEL = "#FFFFFF"
const INK = "#1B1F24"
const INK2 = "#5A6270"
const LINE = "#D9DCD6"

"""One pre-characterised design.  `faces` are 1-based vertex lists (converted from the
0-based file), `C`/`S` the deployment basis of Y(θ) = cos(θ/2)·C + sin(θ/2)·S, `sigma` the
orientation label of the file ("def"/"mc"), `sense[i] = ±1` the rotation sense of face i
(+θ/2 or −θ/2), derived from the basis by `face_sense`."""
struct Design
    id::String
    kind::String
    sigma::String
    sense::Vector{Int}
    F::Int
    nsplit::Int
    theta_max::Float64
    eps_max::Float64
    theta_half::Float64
    faces::Vector{Vector{Int}}
    C::Matrix{Float64}          # n × 2, closed frame v0
    S::Matrix{Float64}          # n × 2
    source::Any                 # nothing for a frame file; GuiBackend.DesignSource when computed
end

"""Deployment basis from two frames: Y(0)=v0=C and Y(θh)=v1 give
S = (v1 − cos(θh/2)·C) / sin(θh/2).  θh must be in (0, 2π)."""
function deploy_basis(v0::AbstractMatrix, v1::AbstractMatrix, theta_half::Real)
    size(v0) == size(v1) || throw(ArgumentError("v0 and v1 differ in size"))
    sh = sin(theta_half / 2)
    abs(sh) > 1e-12 || throw(ArgumentError("theta_half must have sin(theta_half/2) != 0"))
    C = Matrix{Float64}(v0)
    S = (Matrix{Float64}(v1) .- cos(theta_half / 2) .* C) ./ sh
    return C, S
end

"""Rotation sense of every face.  A rigid face turning by s·θ/2 has, for each edge,
eS = s·rot90(eC), so s = sign(eC × eS).  The majority over the face's edges is taken because
the file stores coordinates rounded to 4 decimals."""
function face_sense(C::AbstractMatrix, S::AbstractMatrix, faces)
    sense = Vector{Int}(undef, length(faces))
    for (i, f) in enumerate(faces)
        acc = 0.0
        for k in eachindex(f)
            a, b = f[k], f[mod1(k + 1, length(f))]
            ecx, ecy = C[b, 1] - C[a, 1], C[b, 2] - C[a, 2]
            esx, esy = S[b, 1] - S[a, 1], S[b, 2] - S[a, 2]
            acc += sign(ecx * esy - ecy * esx)
        end
        sense[i] = acc >= 0 ? 1 : -1
    end
    return sense
end

"""Vertex positions at opening angle θ: Y(θ) = cos(θ/2)·C + sin(θ/2)·S (exact path)."""
frame(d::Design, theta::Real) = cos(theta / 2) .* d.C .+ sin(theta / 2) .* d.S

_tomatrix(v) = length(v) == 0 ? zeros(0, 2) : reduce(vcat, (Float64.(p)' for p in v))

"""Build a `Design` from one JSON entry (faces in the file are 0-based; +1 here)."""
function design_from_json(j::AbstractDict)
    faces = [Int.(f) .+ 1 for f in j["faces"]]      # 0-based file -> 1-based Julia
    C, S = deploy_basis(_tomatrix(j["v0"]), _tomatrix(j["v1"]), Float64(j["theta_half"]))
    Design(string(j["id"]), string(get(j, "kind", "")), string(j["sigma"]), face_sense(C, S, faces), Int(j["F"]),
           Int(j["nsplit"]), Float64(j["theta_max"]), Float64(j["eps_max"]),
           Float64(j["theta_half"]), faces, C, S, nothing)
end

has_frame_schema(j) = j isa AbstractDict && all(haskey(j, k) for k in ("id", "faces", "v0", "v1", "theta_half", "theta_max", "eps_max", "sigma"))

"""Load every design in a deploy_frames.json file that carries the frame schema.  Entries of
another shape (e.g. raw UKP text) are skipped and counted in the second return value."""
function load_designs(path::AbstractString)
    raw = JSON.parsefile(path)
    entries = raw isa AbstractDict ? get(raw, "designs", Any[]) : raw
    designs = Design[]
    skipped = 0
    for j in entries
        has_frame_schema(j) ? push!(designs, design_from_json(j)) : (skipped += 1)
    end
    return designs, skipped
end

# ---- backend seam for phase 2: anything that can produce designs
abstract type AbstractBackend end
struct FileBackend <: AbstractBackend
    path::String
end
designs(b::FileBackend) = load_designs(b.path)

# ---- per-θ readouts
@enum DeployState closed certified usable over

"""State badge: closed at θ=0, certified for θ ≤ εmax, usable for θ ≤ Θmax, else over."""
function deploy_state(d::Design, theta::Real)
    theta <= 0 && return closed
    theta <= d.eps_max + 1e-9 && return certified
    theta <= d.theta_max + 1e-9 && return usable
    return over
end

badge_text(s::DeployState) = s == closed ? "closed" : s == certified ? "certified" :
                             s == usable ? "usable (exact scan)" : "past Θmax — collisions"

bbox(P::AbstractMatrix) = size(P, 1) == 0 ? (0.0, 0.0, 0.0, 0.0) :
    (minimum(view(P, :, 1)), minimum(view(P, :, 2)), maximum(view(P, :, 1)), maximum(view(P, :, 2)))

bbox_area(P) = (b = bbox(P); (b[3] - b[1]) * (b[4] - b[2]))

"""Areal expansion at θ: bounding-box area relative to the closed frame (as the web app)."""
function areal_expansion(d::Design, theta::Real)
    a0 = max(bbox_area(d.C), 1e-12)
    return bbox_area(frame(d, theta)) / a0
end

"""Track fractions of the θ ∈ [0, π] slider: (usable end, certified end), clamped to [0,1]."""
track_fractions(d::Design) = (clamp(d.theta_max / pi, 0, 1), clamp(d.eps_max / pi, 0, 1))

"""Bounding box covering both the closed and the fully-open (min(Θmax, π)) frames, so the
camera does not jump while θ moves (web fitView)."""
function view_box(d::Design)
    b0 = bbox(d.C)
    bt = bbox(frame(d, min(d.theta_max, pi)))
    return (min(b0[1], bt[1]), min(b0[2], bt[2]), max(b0[3], bt[3]), max(b0[4], bt[4]))
end

"""Ease-in-out θ for the play animation: u ∈ [0,1] → start + (end−start)(½ − ½cos πu)."""
ease_theta(start, stop, u) = start + (stop - start) * (0.5 - 0.5 * cos(pi * clamp(u, 0, 1)))

"""Play target: max(Θmax, 0.05); restart from 0 when already at the end (web behaviour)."""
function play_range(d::Design, theta::Real)
    stop = max(d.theta_max, 0.05)
    start = theta >= stop - 1e-6 ? 0.0 : Float64(theta)
    return start, stop
end

# ---- exports
face_color(sigma::Int) = sigma > 0 ? FACE_A : FACE_B

# Preview SVG (faces filled by sense, no kerf/hinge geometry) for designs loaded from a
# frame file, which carry no mesh. Computed designs use the real cut layout through
# GuiBackend.export_svg_string (Kirigami build_layout + svg_string).
"""Preview SVG of the frame at θ, 1 unit = `scale` px, y up.  Returns the SVG text."""
function svg_string(d::Design, theta::Real; scale::Real = 40.0, margin::Real = 10.0)
    P = frame(d, theta)
    x0, y0, x1, y1 = bbox(P)
    w = (x1 - x0) * scale + 2margin
    h = (y1 - y0) * scale + 2margin
    io = IOBuffer()
    print(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"", round(w; digits=2),
          "\" height=\"", round(h; digits=2), "\" viewBox=\"0 0 ", round(w; digits=2), " ",
          round(h; digits=2), "\">\n<!-- ", d.id, " theta=", round(theta; digits=6), " rad -->\n")
    for (i, f) in enumerate(d.faces)
        pts = join((string(round((P[v, 1] - x0) * scale + margin; digits=3), ",",
                           round((y1 - P[v, 2]) * scale + margin; digits=3)) for v in f), " ")
        print(io, "<polygon points=\"", pts, "\" fill=\"", face_color(d.sense[i]),
              "\" stroke=\"", FACE_EDGE, "\" stroke-width=\"0.7\" stroke-linejoin=\"round\"/>\n")
    end
    print(io, "</svg>\n")
    return String(take!(io))
end

"""JSON-serialisable state of the viewer (faces written back 0-based, like the file)."""
function state_dict(d::Design, theta::Real)
    P = frame(d, theta)
    Dict("id" => d.id, "kind" => d.kind, "theta" => Float64(theta),
         "sigma" => d.sigma, "sense" => d.sense, "faces" => [f .- 1 for f in d.faces],     # 1-based -> 0-based
         "vertices" => [P[i, :] for i in 1:size(P, 1)],
         "state" => string(deploy_state(d, theta)),
         "characterization" => Dict("theta_max" => d.theta_max, "eps_max" => d.eps_max,
                                    "F" => d.F, "nsplit" => d.nsplit))
end

state_json(d::Design, theta::Real) = JSON.json(state_dict(d, theta), 1)

end # module
