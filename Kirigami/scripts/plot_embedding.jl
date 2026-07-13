#!/usr/bin/env julia
# Plot the uncut mesh M and the kirigami structure M' at several opening angles.
# Reads the JSON produced by the C++ tools (kiri_analyze / kiri_deploy / kiri_reference).
# Port of plot_embedding.py.
#
# usage:
#   plot_embedding.jl case <case_dir> [-o out.png]     # M.json + deploy_*.json panels
#   plot_embedding.jl hist <sweep.csv> [-o out.png]    # dim_null vs (#interior - H)
include(joinpath(@__DIR__, "plot_common.jl"))

const CW = "#5B8FF9"    # sigma = +1  (clockwise)
const CCW = "#F6BD16"   # sigma = -1  (counter-clockwise)
const HOLE = "#E8684A"

function draw_mesh!(ax, path, title)
    d = JSON.parsefile(path)
    V = d["vertices"]
    sig = get(d, "orientation", nothing)
    faces = d["faces"]
    cols = [sig !== nothing && sig[i] == 1 ? CW : CCW for i in eachindex(faces)]
    draw_faces!(ax, V, faces; color = parse.(Makie.Colors.Colorant, cols),
                strokecolor = GREY(0.25), strokewidth = 0.4)
    frame_axis!(ax, V)
    ax.title = title
end

function draw_deployment!(ax, path, title)
    d = JSON.parsefile(path)
    V = d["vertices"]
    draw_faces!(ax, V, d["faces"]; color = "#BDD2FD", strokecolor = GREY(0.25),
                strokewidth = 0.4)
    holes = [h for h in get(d, "holes", []) if length(h) >= 3]
    isempty(holes) || draw_faces!(ax, V, holes; color = (HOLE, 0.45), strokewidth = 0)
    frame_axis!(ax, V)
    ax.title = "$title  ($(length(get(d, "holes", []))) holes)"
end

function plot_case(case_dir, out)
    panels = [("M.json", "M (input, sigma coloured)"),
              ("X0.json", "X0 (Eq. 6 projection)"),
              ("deploy_30deg.json", "M' at theta = 30 deg"),
              ("deploy_60deg.json", "M' at theta = 60 deg"),
              ("deploy_thetamax.json", "M' at theta_max")]
    have = [(f, t) for (f, t) in panels if isfile(joinpath(case_dir, f))]
    fig = Figure(size = (310 * length(have), 330))
    for (i, (f, t)) in enumerate(have)
        ax = Axis(fig[1, i]; titlesize = 11)
        p = joinpath(case_dir, f)
        startswith(f, "deploy") ? draw_deployment!(ax, p, t) : draw_mesh!(ax, p, t)
    end
    Label(fig[0, 1:length(have)], basename(normpath(case_dir)); fontsize = 13)
    savefig(fig, out)
end

function plot_hist(csv_path, out)
    R = rows(csv_path)
    dim = [inum(r, "dim_null") for r in R]
    claim = [inum(r, "n_interior_minus_H") for r in R]
    diff = dim .- claim
    fig = Figure(size = (1000, 400))
    ax1 = Axis(fig[1, 1]; xlabel = "#interior vertices - H", ylabel = "dim null(Eq. 4)",
               title = "Rank claim, $(length(dim)) sweep graphs")
    scatter!(ax1, claim, dim; markersize = 7, color = (CW, 0.7))
    lim = max(maximum(claim; init = 1), maximum(dim; init = 1)) * 1.05
    lines!(ax1, [0, lim], [0, lim]; color = :black, linestyle = :dash, linewidth = 0.8)
    ax2 = Axis(fig[1, 2]; xlabel = "dim_null - (#interior - H)", ylabel = "count",
               title = "Deviation (0 = claim holds)")
    lo, hi = minimum(diff; init = 0) - 1, maximum(diff; init = 0) + 2
    hist_bars!(ax2, diff, collect(lo:hi); color = CCW, strokecolor = GREY(0.3), strokewidth = 0.5)
    savefig(fig, out)
    println("wrote ", out, " | violations: ", count(!=(0), diff), " of ", length(diff))
end

function main()
    length(ARGS) >= 2 || error("usage: plot_embedding.jl case|hist <path> [-o out.png]")
    mode, path = ARGS[1], ARGS[2]
    out = nothing
    for i in 3:length(ARGS)
        ARGS[i] in ("-o", "--out") && (out = ARGS[i+1])
    end
    if mode == "case"
        plot_case(path, outpath(something(out, joinpath(path, "figure.png"))))
    elseif mode == "hist"
        plot_hist(path, outpath(something(out, "rank_claim.png")))
    else
        error("mode must be case or hist")
    end
end

main()
