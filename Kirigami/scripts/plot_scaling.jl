#!/usr/bin/env julia
# fig_scaling.png + the fitted exponents of SCALING.md, from results/scaling/scaling.csv.
# Log-log wall time against |F| per routine with an OLS fit of
# log10(t) on log10(|F|) (95 % interval from a small t table), and peak RSS against |F|.
# Rows below the driver's 1 ms clock resolution are excluded from the fit and drawn hollow.
# Run: julia --project=Kirigami/scripts Kirigami/scripts/plot_scaling.jl
include(joinpath(@__DIR__, "plot_common.jl"))

const SC = joinpath(ROOT, "results", "scaling")
const TCRIT = Dict(1 => 12.706, 2 => 4.303, 3 => 3.182, 4 => 2.776, 5 => 2.571, 6 => 2.447,
                   7 => 2.365, 8 => 2.306, 9 => 2.262, 10 => 2.228, 11 => 2.201, 12 => 2.179,
                   13 => 2.160, 14 => 2.145, 15 => 2.131, 16 => 2.120, 17 => 2.110, 18 => 2.101,
                   19 => 2.093, 20 => 2.086, 21 => 2.080, 22 => 2.074, 23 => 2.069, 24 => 2.064,
                   25 => 2.060, 26 => 2.056, 27 => 2.052, 28 => 2.048, 29 => 2.045, 30 => 2.042)
const ROUTINES = [("solve_dense", "assemble + solve_system (dense SVD)", "#b03030"),
                  ("rank_sparse", "rank_only_sparse (sparse QR)", "#c07a20"),
                  ("characterize", "characterize (exact scan + certificate)", "#1f6f4a"),
                  ("design_range_max", "design_range_max (incl. its own solve)", "#3050b0")]
const MIN_T = 0.002   # the driver's clock is milliseconds; below this the datum is quantised

"OLS slope of log10(t) on log10(F); returns (slope, lo, hi, n, r2)"
function fit(F, t)
    x, y = log10.(F), log10.(t)
    n = length(x)
    n < 3 && return (NaN, NaN, NaN, n, NaN)
    A = hcat(x, ones(n))
    beta = A \ y
    resid = y - A * beta
    dof = n - 2
    s2 = (resid' * resid) / dof
    cov = s2 * inv(A' * A)
    se = sqrt(cov[1, 1])
    tc = get(TCRIT, dof, 1.96)
    ss_tot = sum((y .- mean(y)) .^ 2)
    r2 = ss_tot > 0 ? 1 - (resid' * resid) / ss_tot : NaN
    return beta[1], beta[1] - tc * se, beta[1] + tc * se, n, r2
end

function main()
    R = rows(joinpath(SC, "scaling.csv"))
    fig = Figure(size = (1260, 500))
    ax0 = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = L"faces $|F|$",
               ylabel = "wall time [s]", titlesize = 12,
               title = L"wall time vs. $|F|$, fitted $t \propto |F|^{\alpha}$")
    ax1 = Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = L"faces $|F|$",
               ylabel = "peak RSS [MB]  (ru_maxrss at exit)", titlesize = 12,
               title = L"peak resident set vs. $|F|$")
    lines = []
    for (key, label, color) in ROUTINES
        sel = [r for r in R if r["routine"] == key && r["status"] == "ok"]
        isempty(sel) && continue
        F = [fnum(r, "F") for r in sel]
        t = [fnum(r, "secs") for r in sel]
        rss = [fnum(r, "rss_mb") for r in sel]
        good = t .>= MIN_T
        any(.!good) && scatter!(ax0, F[.!good], max.(t[.!good], 1e-4); markersize = 8,
                                color = :transparent, strokecolor = (color, 0.6), strokewidth = 1)
        scatter!(ax0, F[good], t[good]; markersize = 10, color = color)
        s, lo, hi, n, r2 = fit(F[good], t[good])
        if !isnan(s)
            xx = [minimum(F[good]), maximum(F[good])]
            c = mean(log10.(t[good])) - s * mean(log10.(F[good]))
            lines!(ax0, xx, 10 .^ (c .+ s .* log10.(xx)); linewidth = 1.6, color = color,
                   label = @sprintf("%s: α = %.2f [%.2f, %.2f], n=%d, R²=%.3f", label, s, lo, hi, n, r2))
        end
        push!(lines, (key, label, s, lo, hi, n, r2, Int(maximum(F)), maximum(t), maximum(rss)))
        ord = sortperm(F)
        scatterlines!(ax1, F[ord], rss[ord]; markersize = 8, linewidth = 1.2, color = color, label = label)
    end
    axislegend(ax0; position = :lt, labelsize = 9)
    axislegend(ax1; position = :lt, labelsize = 9)
    savefig(fig, outpath(joinpath(SC, "fig_scaling.png")))

    fits = outpath(joinpath(SC, "fits.csv"))
    open(fits, "w") do f
        write(f, "routine,label,exponent,ci_lo,ci_hi,n_fit,r2,max_F_completed,max_secs,max_rss_mb\n")
        for r in lines
            @printf(f, "%s,%s,%.4f,%.4f,%.4f,%d,%.4f,%d,%.4f,%.2f\n", r...)
        end
    end
    println("wrote ", fits)
    foreach(println, lines)
end

main()
