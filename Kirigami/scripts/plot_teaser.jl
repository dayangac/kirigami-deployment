#!/usr/bin/env julia
# Render the poster teaser panels (paper/poster/poster.tex, a sibling folder outside the repo) (closed / half-open / 0.9*Theta_max) directly from the
# kiri_export SVGs for hero2 (export/hero2/*.svg, "cut" layer = one closed straight-line
# polygon path per face). Each state on its own square canvas,
# filled faces in mid-blue, dark thin edges, no axes, tight centred square crop.
# The living-hinge neck notches make each path self-intersecting; a small Ramer-Douglas-
# Peucker pass (RDP_EPS_MM, in the SVG's mm units) recovers the macro polygon first.
# Run: julia --project=Kirigami/scripts Kirigami/scripts/plot_teaser.jl
include(joinpath(@__DIR__, "plot_common.jl"))

const FACE_COLOR = "#5F82A8"
const EDGE_COLOR = "#1a1a1a"
const EDGE_LW = 0.6
const SIZE_PX = 1400
const MARGIN_FRAC = 0.04  # padding around the shape's own bounding box
const RDP_EPS_MM = 3.0

const STATES = [("closed", "hero2_130_sigma_mc_closed.svg", "teaser_closed.png"),
                ("half", "hero2_130_sigma_mc_open_half.svg", "teaser_half.png"),
                ("open90", "hero2_130_sigma_mc_open_0.9tm.svg", "teaser_open90.png")]
const SRC_DIR = "export/hero2"
const OUT_DIR = "results/final/figures"  # read by paper/poster/poster.tex via \graphicspath

function perp_dist(pt, a, b)
    (x, y), (x1, y1), (x2, y2) = pt, a, b
    dx, dy = x2 - x1, y2 - y1
    dx == 0 && dy == 0 && return hypot(x - x1, y - y1)
    t = ((x - x1) * dx + (y - y1) * dy) / (dx * dx + dy * dy)
    return hypot(x - (x1 + t * dx), y - (y1 + t * dy))
end

function rdp(points, eps)
    length(points) < 3 && return points
    dmax, idx = 0.0, 1
    for i in 2:length(points)-1
        d = perp_dist(points[i], points[1], points[end])
        d > dmax && ((dmax, idx) = (d, i))
    end
    if dmax > eps
        left = rdp(points[1:idx], eps)
        right = rdp(points[idx:end], eps)
        return vcat(left[1:end-1], right)
    end
    return [points[1], points[end]]
end

const NUM = r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"

# The "cut" <g> block, then every path's d attribute inside it.
function face_polygons(svg_path)
    txt = read(svg_path, String)
    m = match(r"<g[^>]*id=\"cut\"[^>]*>(.*?)</g>"s, txt)
    m === nothing && error("no <g id=\"cut\"> in $svg_path")
    polys = Vector{Tuple{Float64,Float64}}[]
    for pm in eachmatch(r"<path[^>]*\sd=\"([^\"]*)\"", m.captures[1])
        tokens = [t.match for t in eachmatch(Regex("[MLZmlz]|" * NUM.pattern), pm.captures[1])]
        pts = Tuple{Float64,Float64}[]
        i = 1
        while i <= length(tokens)
            t = tokens[i]
            if t in ("M", "L", "m", "l")
                push!(pts, (parse(Float64, tokens[i+1]), parse(Float64, tokens[i+2])))
                i += 3
            else
                i += 1
            end
        end
        if length(pts) >= 3
            simplified = rdp(vcat(pts, [pts[1]]), RDP_EPS_MM)[1:end-1]
            push!(polys, length(simplified) >= 3 ? simplified : pts)
        end
    end
    return polys
end

function render(name, svg_name, out_name)
    polys = face_polygons(joinpath(SRC_DIR, svg_name))
    xs = [x for poly in polys for (x, y) in poly]
    ys = [y for poly in polys for (x, y) in poly]
    xmin, xmax, ymin, ymax = minimum(xs), maximum(xs), minimum(ys), maximum(ys)
    w, h = xmax - xmin, ymax - ymin
    half = max(w, h) * (1.0 + MARGIN_FRAC) / 2.0
    cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2

    fig = Figure(size = (SIZE_PX, SIZE_PX), backgroundcolor = :white, figure_padding = 0)
    ax = Axis(fig[1, 1]; aspect = DataAspect(), backgroundcolor = :white, yreversed = true)  # SVG y points down
    poly!(ax, [Point2f.(p) for p in polys]; color = FACE_COLOR, strokecolor = EDGE_COLOR,
          strokewidth = EDGE_LW)
    limits!(ax, cx - half, cx + half, cy - half, cy + half)
    hidedecorations!(ax); hidespines!(ax)
    out_path = outpath(joinpath(OUT_DIR, out_name))
    save(out_path, fig; px_per_unit = 1)
    @printf("%s: %d faces, bbox %.1f x %.1f mm\n", out_path, length(polys), w, h)
end

for (name, svg_name, out_name) in STATES
    render(name, svg_name, out_name)
end
