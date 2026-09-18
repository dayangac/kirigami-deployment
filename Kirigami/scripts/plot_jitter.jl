#!/usr/bin/env julia
# A3 jitter-transition figure: deployable / certified fraction against jitter amplitude,
# one curve per authored tiling per sigma rule, drawn from ladder.csv written by
# `kill_jitter --analyze`.
# Run: julia --project=Kirigami/scripts Kirigami/scripts/plot_jitter.jl [results/kill/jitter]
include(joinpath(@__DIR__, "plot_common.jl"))

D = length(ARGS) >= 1 ? ARGS[1] : "results/kill/jitter"

R = rows(joinpath(D, "ladder.csv"))
# (tiling, rule) -> (a, f_cert, f_theta, med)
series = Dict{Tuple{String,String},NTuple{4,Vector{Float64}}}()
nsplit = Dict{String,Int}()
for r in R
    inum(r, "amp_idx") < 0 && continue  # the a = 0 baseline is drawn as a marker, not a point
    k = (r["tiling"], r["rule"])
    s = get!(series, k) do
        (Float64[], Float64[], Float64[], Float64[])
    end
    push!(s[1], fnum(r, "amp"))
    push!(s[2], fnum(r, "frac_cert"))
    push!(s[3], fnum(r, "frac_theta_pos"))
    push!(s[4], r["median_theta"] in ("", "nan") ? NaN : fnum(r, "median_theta"))
    if r["rule"] == "mc"  # |E_split| is sigma-dependent; label with the max-cut value
        nsplit[r["tiling"]] = inum(r, "n_split")
    end
end

tilings = sort(unique(first.(collect(keys(series)))))
# split-free tilings first, in a muted grey; split-bearing ones get the colour cycle
free = [t for t in tilings if nsplit[t] == 0]
bear = [t for t in tilings if nsplit[t] > 0]
colour = Dict{String,Any}(t => GREY(0.65) for t in free)
for (i, t) in enumerate(bear)
    colour[t] = MPL_CYCLE[(i-1) % length(MPL_CYCLE) + 1]
end

fig = Figure(size = (1450, 520))
titles = [L"fraction with exact $\Theta_{max} > 0$",
          L"fraction CERTIFIED (POS $\wedge$ NOOVERLAP $\wedge$ NOROOT, $\epsilon$ = 0.006)"]
axes = Axis[]
for (col, idx) in ((1, 3), (2, 2))
    ax = Axis(fig[1, col]; xscale = log10, title = titles[col], titlesize = 11,
              xlabel = L"jitter amplitude $a$  (median edge lengths, interior vertices)",
              ylabel = col == 1 ? "fraction of 20 seeds" : "")
    push!(axes, ax)
    hlines!(ax, [0.5]; color = GREY(0.8), linewidth = 0.8)
    for t in vcat(free, bear), (rule, ls) in (("mc", :solid), ("def", :dash))
        k = (t, rule)
        haskey(series, k) || continue
        a = series[k][1]
        y = series[k][idx]
        lines!(ax, a, y; color = colour[t], linestyle = ls,
               linewidth = rule == "mc" ? 1.8 : 1.2)
        rule == "mc" && scatter!(ax, a, y; color = colour[t], markersize = 6)
    end
    ylims!(ax, -0.04, 1.04)
end

ax3 = Axis(fig[1, 3]; xscale = log10, xlabel = L"jitter amplitude $a$",
           title = L"median exact $\Theta_{max}$ (rad), $\sigma_{mc}$", titlesize = 11)
push!(axes, ax3)
handles = []
labels = String[]
for t in vcat(free, bear)
    a, fc, ft, med = series[(t, "mc")]
    keep = .!isnan.(med)
    l = lines!(ax3, a[keep], med[keep]; color = colour[t], linewidth = 1.8)
    scatter!(ax3, a[keep], med[keep]; color = colour[t], markersize = 6)
    push!(handles, l)
    push!(labels, "$t  (|E_split| = $(nsplit[t]))")
end
Legend(fig[2, 1:3], handles, labels; orientation = :horizontal, nbanks = 2,
       framevisible = false, labelsize = 10)
Label(fig[0, 1:3],
      "A3 jitter transition: authored tilings, combinatorics fixed, interior geometry perturbed\n" *
      "(solid = σ max-cut, dashed = σ defect-minimising; grey = split-free patterns, whose Eq. (6) projection undoes the jitter exactly)";
      fontsize = 12)
out = outpath(joinpath(D, "transition.png"))
savefig(fig, out)
