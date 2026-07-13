#!/usr/bin/env julia
# K6 figure: inward split edges at t = 0, and the certified Theta_max distributions.
# Port of plot_k6.py.
# Run: julia --project=Kirigami/scripts Kirigami/scripts/plot_k6.jl results/kill/k6/k6.csv \
#          -o results/kill/k6/k6_zero_plus.png
include(joinpath(@__DIR__, "plot_common.jl"))
using Random

const FAMS = ["sigma_mc", "sigma_def"]
const COL = Dict("sigma_mc" => "#3b6ea5", "sigma_def" => "#c2622d")

function parse_args(args)
    csvp, out = nothing, "k6.png"
    i = 1
    while i <= length(args)
        if args[i] in ("-o", "--out")
            out = args[i+1]; i += 2
        else
            csvp = args[i]; i += 1
        end
    end
    csvp === nothing && error("usage: plot_k6.jl <csv> [-o out.png]")
    return csvp, out
end

function main()
    csvp, out = parse_args(ARGS)
    byfam = Dict(f => Row[] for f in FAMS)
    for r in rows(csvp)
        haskey(byfam, r["sigma"]) && push!(byfam[r["sigma"]], r)
    end
    all(isempty, values(byfam)) && error("no rows")
    col(rs, key) = [fnum(r, key) for r in rs]

    fig = Figure(size = (1500, 440))

    # (a) inward split edges at t = 0, as a fraction of |E_split|
    ax1 = Axis(fig[1, 1]; xlabel = L"inward split edges $q_e\leq 0$ / $|E_{split}|$ at $t=0$",
               ylabel = "graphs", title = L"(a) why K5 saw $\Theta_{max}=0$")
    for f in FAMS
        v = col(byfam[f], "inward_frac0")
        isempty(v) || hist_bars!(ax1, v, collect(range(0, 1, length = 26));
                                 color = COL[f], alpha = 0.6, label = f, strokewidth = 0)
    end
    axislegend(ax1)

    # (b) absolute counts vs |E_split|
    ax2 = Axis(fig[1, 2]; xlabel = L"|E_{split}|", ylabel = L"inward split edges at $t=0$",
               title = "(b) count vs split-cut count")
    lim = 1.0
    for f in FAMS
        x, y = col(byfam[f], "n_split"), col(byfam[f], "inward0")
        isempty(x) && continue
        scatter!(ax2, x, y; markersize = 7, color = (COL[f], 0.6), label = f)
        lim = max(lim, maximum(x))
    end
    lines!(ax2, [0, lim], [0, lim]; color = :black, linestyle = :dash, linewidth = 0.8,
           label = "all inward")
    axislegend(ax2)

    # (c) certified Theta_max after the repair, plus the refereeing bisection
    labels, data, colors = String[], Vector{Float64}[], String[]
    for f in FAMS
        rs = [r for r in byfam[f] if inum(r, "feasible") != 0]
        for (key, tag) in (("eps_max", "certified"), ("theta_bisect", "bisection"))
            v = [fnum(r, key) for r in rs]
            push!(labels, "$f\n$tag\n(n=$(length(v)))")
            push!(data, isempty(v) ? [0.0] : v)
            push!(colors, COL[f])
        end
    end
    ax3 = Axis(fig[1, 3]; ylabel = L"\Theta_{max}\ \mathrm{(rad)}",
               title = L"(c) after the 0+ repair, at the $0^+$-feasible points",
               xticks = (0:length(data)-1, labels), xticklabelsize = 9)
    rng = MersenneTwister(0)
    for (p, d, c) in zip(0:length(data)-1, data, colors)
        jitter = (rand(rng, length(d)) .- 0.5) .* 0.3
        scatter!(ax3, p .+ jitter, d; markersize = 7, color = (c, 0.7))
    end
    ylims!(ax3, -0.05, max(0.35, maximum(maximum.(data)) * 1.1))

    Label(fig[0, 1:3], "K6 -- 0+ repair in the Tutte auxetic null space, 200 random graphs";
          fontsize = 14)
    savefig(fig, outpath(out))
end

main()
