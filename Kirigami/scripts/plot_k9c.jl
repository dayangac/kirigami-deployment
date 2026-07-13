#!/usr/bin/env julia
# K9c figures (port of plot_k9c.py).
#   k9c_gallery.png -- the 20 best designs by certified eps_max, closed and at Theta_max/2.
#   k9c_hist.png    -- certified eps_max distribution + per-solver counts.
# Run: julia --project=Kirigami/scripts Kirigami/scripts/plot_k9c.jl [results/kill/k9c]
include(joinpath(@__DIR__, "plot_common.jl"))

D = length(ARGS) >= 1 ? ARGS[1] : "results/kill/k9c"
G = joinpath(D, "gallery")

const CLOSED = "#BDD2FD"
const OPEN = "#5B8FF9"

function load_index()
    all = Any[]
    for p in sort(filter(x -> startswith(basename(x), "index_") && endswith(x, ".json"),
                         readdir(G; join = true)))
        append!(all, JSON.parsefile(p))
    end
    if isempty(all)
        # Rebuild the same records from the dumped geometry plus k9c.csv, so a partially
        # aggregated run still produces the gallery.
        csvp = joinpath(D, "k9c.csv")
        if isfile(csvp)
            by_tag = Dict("$(r["kind"])_$(r["id"])_$(r["sigma"])" => r for r in rows(csvp))
            for q in sort(filter(x -> endswith(x, "_closed.json"), readdir(G; join = true)))
                tag = basename(q)[1:end-length("_closed.json")]
                r = get(by_tag, tag, nothing)
                r === nothing && continue
                push!(all, Dict("tag" => tag, "id" => inum(r, "id"), "kind" => r["kind"],
                                "sigma" => r["sigma"], "F" => inum(r, "F"),
                                "theta_max" => fnum(r, "theta_exact"),
                                "eps_max" => fnum(r, "eps_max"), "src" => r["best_src"]))
            end
        end
    end
    # one entry per graph: the sigma with the larger certified eps_max (best-of-2)
    best = Dict{Tuple{String,Int},Any}()
    for r in all
        k = (r["kind"], Int(r["id"]))
        if !haskey(best, k) || r["eps_max"] > best[k]["eps_max"]
            best[k] = r
        end
    end
    return sort(collect(values(best)); by = r -> -r["eps_max"])
end

function draw!(ax, path, color, title)
    d = JSON.parsefile(path)
    V = d["vertices"]
    draw_faces!(ax, V, d["faces"]; color, strokecolor = GREY(0.3), strokewidth = 0.15)
    frame_axis!(ax, V; margin = 0.04)
    ax.title = title
end

function gallery(n_want = 20)
    entries = [e for e in load_index() if isfile(joinpath(G, e["tag"] * "_closed.json")) &&
                                          isfile(joinpath(G, e["tag"] * "_half.json"))]
    if isempty(entries)
        println("gallery: nothing dumped")
        return 0
    end
    entries = entries[1:min(n_want, length(entries))]
    ncol = 8
    nrow = (length(entries) + 3) ÷ 4
    fig = Figure(size = (155 * ncol, 175 * nrow + 60))
    for (i, e) in enumerate(entries)
        r, cpair = divrem(i - 1, 4)
        draw!(Axis(fig[r + 1, 2cpair + 1]; titlesize = 8),
              joinpath(G, e["tag"] * "_closed.json"), CLOSED,
              @sprintf("%s %d %s\nF=%d  closed", first(e["kind"], 4), e["id"],
                       endswith(e["sigma"], "mc") ? "mc" : "def", e["F"]))
        draw!(Axis(fig[r + 1, 2cpair + 2]; titlesize = 8),
              joinpath(G, e["tag"] * "_half.json"), OPEN,
              @sprintf("θ = Θmax/2\n%.3f rad, ε=%.2f", 0.5 * e["theta_max"], e["eps_max"]))
    end
    Label(fig[0, 1:ncol],
          "K9c: the $(length(entries)) best random graphs by certified εmax under the RANGE-MAXIMISING " *
          "constrained embedding\n(closed and half deployed; every baseline on the same 400 designs is Θmax = 0)";
          fontsize = 11)
    out = outpath(joinpath(D, "k9c_gallery.png"))
    savefig(fig, out)
    println("wrote $out with $(length(entries)) designs")
    return length(entries)
end

function hist()
    path = joinpath(D, "k9c.csv")
    isfile(path) || (println("no k9c.csv"); return)
    R = rows(path)
    eps = [fnum(r, "eps_max") for r in R if fnum(r, "eps_max") > 0]
    fig = Figure(size = (1040, 380))
    ax0 = Axis(fig[1, 1]; xlabel = L"certified $\varepsilon_{max}$ (rad)", ylabel = "designs",
               title = "K9c best-of-3: $(length(eps)) certified designs of $(length(R))",
               titlesize = 11)
    edges = linbins(eps, 24)
    hist_bars!(ax0, eps, edges; color = "#5B8FF9", strokecolor = GREY(0.25), strokewidth = 0.5)
    vlines!(ax0, [0.1]; color = "#E8684A", linestyle = :dash, linewidth = 1.2)
    text!(ax0, 0.105, maximum(bincounts(eps, edges); init = 1) * 0.92;
          text = "0.1 rad target", color = "#E8684A", fontsize = 10)

    k9 = count(r -> fnum(r, "k9_theta") > 1e-9, R)
    k9b = count(r -> fnum(r, "k9b_theta") > 1e-9, R)
    k9c = count(r -> fnum(r, "k9c_theta") > 1e-9, R)
    best = count(r -> fnum(r, "theta_exact") > 1e-9 && fnum(r, "theta_ref9") > 1e-9, R)
    names = ["Eq. (6)\n(K1a)", "σ_mc\n(K5)", "σ_def\n(K5)", "0+ repair\n(K6)",
             "convexity\nonly (K9a)", "prox. arm\n(K9)", "prox. arm\n(K9b-like)",
             "K9c\nrange-max", "best of\nthree"]
    vals = [0, 0, 0, 0, 0, k9, k9b, k9c, best]
    cols = vcat(fill(GREY(0.7), 5),
                parse.(Makie.Colors.Colorant, ["#9BB7F0", "#9BB7F0", "#5B8FF9", "#3363C7"]))
    ax1 = Axis(fig[1, 2]; xticks = (0:length(vals)-1, names), xticklabelsize = 8,
               ylabel = L"designs with exact $\Theta_{max} > 0$",
               title = "of the $(length(R)) designs measured (population: 200 graphs x 2 sigma = 400)",
               titlesize = 11)
    barplot!(ax1, 0:length(vals)-1, vals; color = cols, strokecolor = GREY(0.25), strokewidth = 0.5)
    hlines!(ax1, [40]; color = "#E8684A", linestyle = :dash, linewidth = 1.2)
    text!(ax1, 0.05, 41; text = "PASS bar 40 / 400", color = "#E8684A", fontsize = 10)
    for (i, v) in enumerate(vals)
        text!(ax1, i - 1, v + 0.8; text = string(v), align = (:center, :bottom), fontsize = 10)
    end
    savefig(fig, outpath(joinpath(D, "k9c_hist.png")))
end

gallery()
hist()
