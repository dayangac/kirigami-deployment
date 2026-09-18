#!/usr/bin/env julia
# E1 -- histogram of |Theta_max(T4.2") - bisect| over the extended population.
# Run: julia --project=Kirigami/scripts Kirigami/scripts/plot_e1.jl
include(joinpath(@__DIR__, "plot_common.jl"))

const CSVP = "results/final/e1/e1.csv"
const OUT = outpath("results/final/e1/e1_gap.png")

gaps = Float64[]
for r in rows(CSVP)
    g = tryparse(Float64, r["gap"])
    g === nothing || push!(gaps, g)
end

fig = Figure(size = (7 * 100, 4.5 * 100))
ax = Axis(fig[1, 1];
          xlabel = L"|\Theta_{max}(\mathrm{T4.2''}) - \Theta_{max}(\mathrm{bisection})|\ \mathrm{(rad)}",
          ylabel = "count (log)", yscale = log10,
          title = "E1: exact vs. bisection referee, N=$(length(gaps))")
gap_hist_log!(ax, gaps; color = "#3b6ea5")
vlines!(ax, [1e-5]; color = :crimson, linestyle = :dash, linewidth = 1,
        label = "1e-5 agreement bar")
axislegend(ax)
savefig(fig, OUT)
@printf("wrote %s, N=%d, worst=%g\n", OUT, length(gaps), isempty(gaps) ? NaN : maximum(gaps))
