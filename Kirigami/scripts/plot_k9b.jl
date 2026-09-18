#!/usr/bin/env julia
# K9b figures.
#   k9b_gallery.png -- >= 20 deployed random graphs, each shown closed (theta = 0) and at
#                      theta = Theta_max / 2, from the JSON dumped by kill_k9b.
#   k9b_hist.png    -- the certified eps_max distribution of K9b against the baselines.
# Run: julia --project=Kirigami/scripts Kirigami/scripts/plot_k9b.jl [results/kill/k9b]
include(joinpath(@__DIR__, "plot_common.jl"))

D = length(ARGS) >= 1 ? ARGS[1] : "results/kill/k9b"
G = joinpath(D, "gallery")

const CLOSED = "#BDD2FD"
const OPEN = "#5B8FF9"

function load_index()
    all = Any[]
    for p in sort(filter(x -> startswith(basename(x), "index_") && endswith(x, ".json"),
                         readdir(G; join = true)))
        append!(all, JSON.parsefile(p))
    end
    # one entry per graph: the sigma with the larger Theta_max (the honest best-of-2)
    best = Dict{Tuple{String,Int},Any}()
    for r in all
        k = (r["kind"], Int(r["id"]))
        if !haskey(best, k) || r["theta_max"] > best[k]["theta_max"]
            best[k] = r
        end
    end
    return sort(collect(values(best)); by = r -> -r["theta_max"]), all
end

function draw!(ax, path, color, title)
    d = JSON.parsefile(path)
    V = d["vertices"]
    draw_faces!(ax, V, d["faces"]; color, strokecolor = GREY(0.3), strokewidth = 0.15)
    frame_axis!(ax, V; margin = 0.04)
    ax.title = title
end

function gallery(n_want = 20)
    entries, _ = load_index()
    entries = [e for e in entries if isfile(joinpath(G, e["tag"] * "_closed.json")) &&
                                     isfile(joinpath(G, e["tag"] * "_half.json"))]
    if isempty(entries)
        println("gallery: nothing dumped")
        return 0
    end
    entries = entries[1:min(length(entries), max(n_want, min(length(entries), 24)))]
    ncol = 8  # 4 designs per row, each as a (closed, half) pair
    nrow = (length(entries) + 3) ÷ 4
    fig = Figure(size = (155 * ncol, 175 * nrow + 60))
    for (i, e) in enumerate(entries)
        r, cpair = divrem(i - 1, 4)
        a0 = Axis(fig[r + 1, 2cpair + 1]; titlesize = 8)
        a1 = Axis(fig[r + 1, 2cpair + 2]; titlesize = 8)
        draw!(a0, joinpath(G, e["tag"] * "_closed.json"), CLOSED,
              @sprintf("%s %d %s\nF=%d  closed", first(e["kind"], 4), e["id"],
                       endswith(e["sigma"], "mc") ? "mc" : "def", e["F"]))
        draw!(a1, joinpath(G, e["tag"] * "_half.json"), OPEN,
              @sprintf("θ = Θmax/2\n%.3f rad", 0.5 * e["theta_max"]))
    end
    Label(fig[0, 1:ncol],
          "K9b: $(length(entries)) random graphs that deploy under the convexity + split-inward " *
          "constrained embedding\n(every baseline on the same population is Θmax = 0 on all 400 designs)";
          fontsize = 11)
    out = outpath(joinpath(D, "k9b_gallery.png"))
    savefig(fig, out)
    println("wrote $out with $(length(entries)) designs")
    return length(entries)
end

# Per design, the better of K9's solver configuration and K9b's sweep. Both are exact-verified
# positives under the same certificate on the same 400 designs, so the per-design maximum is a
# legitimate multi-start count -- reported next to the two runs, never in place of either.
function union_positives(R)
    path = "results/kill/k9/k9.csv"
    isfile(path) || return count(r -> fnum(r, "theta_ref9") > 1e-9, R)
    K = Dict((r["kind"], r["id"], r["sigma"]) => r for r in rows(path))
    n = 0
    for r in R
        k = (r["kind"], r["id"], r["sigma"])
        best = fnum(r, "theta_exact")
        haskey(K, k) && (best = max(best, fnum(K[k], "b_theta_exact")))
        best > 1e-9 && (n += 1)
    end
    return n
end

function hist()
    path = joinpath(D, "k9b.csv")
    isfile(path) || (println("no k9b.csv"); return)
    R = rows(path)
    eps = [fnum(r, "eps_max") for r in R if fnum(r, "eps_max") > 0]
    fig = Figure(size = (1000, 380))
    ax0 = Axis(fig[1, 1]; xlabel = L"certified $\varepsilon_{max}$ (rad)", ylabel = "designs",
               title = "K9b: $(length(eps)) certified designs of $(length(R))", titlesize = 11)
    hist_bars!(ax0, eps, linbins(eps, 24); color = "#5B8FF9", strokecolor = GREY(0.25),
               strokewidth = 0.5)
    vlines!(ax0, [0.1]; color = "#E8684A", linestyle = :dash, linewidth = 1.2)
    ymax = maximum(bincounts(eps, linbins(eps, 24)); init = 1)
    text!(ax0, 0.105, ymax * 0.92; text = "0.1 rad target", color = "#E8684A", fontsize = 10)

    # Every baseline on this exact population is zero except the K6 repair, read from k6.csv
    # under its best-of-four rule.
    names = ["Eq. (6)\n(K1a)", "σ_mc\n(K5)", "σ_def\n(K5)", "0+ repair\n(K6, best of 4)",
             "convexity\nonly (K9a)", "K9 (b)", "K9b\nsweep", "best of\nboth"]
    k9 = 36
    # A positive needs BOTH the exact scan and the referee: delaunay 40 sigma_mc has
    # bisection(1e-9) = 2.07e-7 rad, below the bisection's own grid resolution, with the exact
    # scan at 0. It is bisection noise, not a design that deploys.
    k9b = count(r -> fnum(r, "theta_exact") > 1e-9 && fnum(r, "theta_ref9") > 1e-9, R)
    union = union_positives(R)
    vals = [0, 0, 0, k6_baseline_deployable(), 0, k9, k9b, union]
    cols = vcat(fill(GREY(0.7), 5), parse.(Makie.Colors.Colorant, ["#9BB7F0", "#9BB7F0", "#5B8FF9"]))
    ax1 = Axis(fig[1, 2]; xticks = (0:length(vals)-1, names), xticklabelsize = 9,
               ylabel = L"designs with refereed $\Theta_{max} > 0$",
               title = "of 400 designs (200 graphs x 2 sigma)", titlesize = 11)
    barplot!(ax1, 0:length(vals)-1, vals; color = cols, strokecolor = GREY(0.25), strokewidth = 0.5)
    hlines!(ax1, [40]; color = "#E8684A", linestyle = :dash, linewidth = 1.2)
    text!(ax1, 0.05, 41; text = "PASS bar 40 / 400", color = "#E8684A", fontsize = 10)
    for (i, v) in enumerate(vals)
        text!(ax1, i - 1, v + 0.8; text = string(v), align = (:center, :bottom), fontsize = 10)
    end
    out = outpath(joinpath(D, "k9b_hist.png"))
    savefig(fig, out)
end

gallery()
hist()
