#!/usr/bin/env julia
# fig_regime.png -- deployment yield against |F| from 20 to 800, three arms.
# Joins the paper-regime population (results/regime/regime.csv, |F| in [20,100]) with the
# population every earlier kill experiment used (|F| in [101,793]), so the two halves of the
# face-count axis are read off one figure:
#
#   Eq. (6) baseline   regime.csv base_theta                 | results/kill/k1a/k1a.csv
#                                                              theta_max_X0 (the Eq. (6)
#                                                              projection, one row per graph)
#   native pipeline    regime.csv natour_theta / natcol_theta | results/kill/native200/
#                                                              native200.csv our_theta_exact,
#                                                              completed cells only
#   k9c range-max      regime.csv k9c_theta                   | results/kill/k9c/k9c.csv
#                                                              k9c_theta
#
# Every column above is the same exact scan; a design counts as deployable when the exact
# Theta_max > 1e-9. Bars are Wilson 95 % intervals.
# Run: julia --project=Kirigami/scripts Kirigami/scripts/plot_regime.jl
include(joinpath(@__DIR__, "plot_common.jl"))

const OUT = outpath(joinpath(ROOT, "results", "regime", "fig_regime.png"))
const EDGES = [20, 35, 50, 71, 100, 150, 220, 330, 480, 800]

# -> Dict(arm => [(F, deployable), ...])
function collect_arms()
    arms = Dict("Eq. (6) baseline" => Tuple{Int,Bool}[],
                "authors' native pipeline" => Tuple{Int,Bool}[],
                "k9c range-max" => Tuple{Int,Bool}[])
    for r in rows(joinpath(ROOT, "results", "regime", "regime.csv"))
        r["set"] == "pop" || continue
        F = inum(r, "F")
        r["base_status"] != "not_run" && push!(arms["Eq. (6) baseline"], (F, fnum(r, "base_theta") > 1e-9))
        r["k9c_status"] != "not_run" && push!(arms["k9c range-max"], (F, fnum(r, "k9c_theta") > 1e-9))
        for pre in ("natour", "natcol")
            r[pre * "_status"] == "completed" &&
                push!(arms["authors' native pipeline"], (F, fnum(r, pre * "_theta") > 1e-9))
        end
    end
    k1a = joinpath(ROOT, "results", "kill", "k1a", "k1a.csv")
    isfile(k1a) && for r in rows(k1a)
        push!(arms["Eq. (6) baseline"], (inum(r, "F"), fnum(r, "theta_max_X0") > 1e-9))
    end
    k9c = joinpath(ROOT, "results", "kill", "k9c", "k9c.csv")
    isfile(k9c) && for r in rows(k9c)
        push!(arms["k9c range-max"], (inum(r, "F"), fnum(r, "k9c_theta") > 1e-9))
    end
    nat = joinpath(ROOT, "results", "kill", "native200", "native200.csv")
    isfile(nat) && for r in rows(nat)
        r["status"] == "completed" && r["embedding_ok"] == "1" &&
            push!(arms["authors' native pipeline"], (inum(r, "F"), fnum(r, "our_theta_exact") > 1e-9))
    end
    return arms
end

function main()
    arms = collect_arms()
    colors = Dict("Eq. (6) baseline" => "#b03030", "authors' native pipeline" => "#7a5cc0",
                  "k9c range-max" => "#1f6f4a")
    fig = Figure(size = (840, 500))
    ax = Axis(fig[1, 1]; xscale = log10, xlabel = L"faces $|F|$  (log scale)",
              ylabel = L"designs with exact $\Theta_{\max} > 0$   [%]",
              xticks = ([20, 50, 100, 200, 400, 800], ["20", "50", "100", "200", "400", "800"]),
              title = "Deployment yield vs. face count, three arms\n" *
                      "bars: Wilson 95 % intervals; left of the dashed line is the 2026 paper's own scale",
              titlesize = 12)
    for name in ["Eq. (6) baseline", "authors' native pipeline", "k9c range-max"]
        data = arms[name]
        xs, ys, lo, hi, ns = Float64[], Float64[], Float64[], Float64[], Int[]
        for (a, b) in zip(EDGES[1:end-1], EDGES[2:end])
            sel = [d for (F, d) in data if a <= F < b]
            length(sel) < 3 && continue
            p, l, h = wilson(count(sel), length(sel))
            push!(xs, sqrt(a * b)); push!(ys, 100p); push!(lo, 100(p - l)); push!(hi, 100(h - p))
            push!(ns, length(sel))
        end
        isempty(xs) && continue
        errorbars!(ax, xs, ys, lo, hi; color = colors[name], whiskerwidth = 6, linewidth = 1.2)
        scatterlines!(ax, xs, ys; color = colors[name], markersize = 10, linewidth = 1.8,
                      label = "$name (n=$(sum(ns)))")
    end
    vlines!(ax, [100]; color = GREY(0.6), linestyle = :dash, linewidth = 1)
    text!(ax, 103, 52; text = "WP7a population  |  earlier kill population", fontsize = 10,
          color = GREY(0.35), rotation = pi / 2, align = (:center, :bottom))
    ylims!(ax, -4, 104)
    xlims!(ax, 18, 900)
    axislegend(ax; position = :rc, labelsize = 10)
    savefig(fig, OUT)
    for (name, data) in arms
        println(name, " n = ", length(data))
    end
end

main()
