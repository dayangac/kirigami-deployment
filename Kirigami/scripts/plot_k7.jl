#!/usr/bin/env julia
# K7 figures: Poisson-ratio curves (predicted vs measured), the
# achievable-K dimension histogram, the C4 area-fit deviation and the C3 outcome bars.
# Run: julia --project=Kirigami/scripts Kirigami/scripts/plot_k7.jl [results/experiments/k7]
include(joinpath(@__DIR__, "plot_common.jl"))

D = length(ARGS) >= 1 ? ARGS[1] : "results/experiments/k7"

function nu_figure()
    names = ["hexagons_2x2", "snub_square_2x2", "t3_4_3_12_2x2", "voronoi_torus_0_n20"]
    files = [(n, joinpath(D, "k7_nu_$n.csv")) for n in names]
    files = [(n, p) for (n, p) in files if isfile(p)]
    isempty(files) && return
    fig = Figure(size = (420 * length(files), 360))
    axes = Axis[]
    for (i, (name, path)) in enumerate(files)
        ax = Axis(fig[1, i]; title = name, titlesize = 11, xlabel = L"\theta",
                  ylabel = i == 1 ? L"Poisson ratio $\nu(\theta)$" : "")
        push!(axes, ax)
        R = rows(path)
        series = Dict{Int,NTuple{3,Vector{Float64}}}()
        for r in R
            r["design"] == "random_in_K" || continue
            th, p, m = fnum(r, "theta"), fnum(r, "nu_pred"), fnum(r, "nu_meas")
            abs(p) > 20 && continue
            s = get!(series, inum(r, "dir_deg")) do
                (Float64[], Float64[], Float64[])
            end
            push!(s[1], th); push!(s[2], p); push!(s[3], m)
        end
        for (j, k) in enumerate(sort(collect(keys(series))))
            th, p, m = series[k]
            c = MPL_CYCLE[(j-1) % length(MPL_CYCLE) + 1]
            lines!(ax, th, p; color = c, linewidth = 1.6, label = "predicted, d = $(k)°")
            scatter!(ax, th[1:4:end], m[1:4:end]; color = :transparent, strokecolor = c,
                     strokewidth = 1, markersize = 7)
        end
        # the conformal design of the same pattern
        cth = [fnum(r, "theta") for r in R if r["design"] == "conformal" && inum(r, "dir_deg") == 0]
        cnu = [fnum(r, "nu_pred") for r in R if r["design"] == "conformal" && inum(r, "dir_deg") == 0]
        isempty(cth) || lines!(ax, cth, cnu; color = :black, linestyle = :dash, linewidth = 1.2,
                               label = "conformal design")
    end
    linkyaxes!(axes...)
    axislegend(axes[1]; labelsize = 8)
    Label(fig[0, 1:length(files)],
          "K7 C3: Poisson ratio from the closed-form K (lines) vs forward kinematics (circles)";
          fontsize = 12)
    savefig(fig, outpath(joinpath(D, "k7_nu_curves.png")))
end

function dim_figure()
    path = joinpath(D, "k7_main.csv")
    isfile(path) || return
    R = [r for r in rows(path) if get(r, "dimK", "") != ""]
    c = Dict{Int,Int}(); ub = Dict{Int,Int}()
    for r in R
        c[inum(r, "dimK")] = get(c, inum(r, "dimK"), 0) + 1
        ub[inum(r, "dimK_min4_2k")] = get(ub, inum(r, "dimK_min4_2k"), 0) + 1
    end
    fig = Figure(size = (900, 340))
    ks = 0:4
    ax0 = Axis(fig[1, 1]; xticks = ks, xlabel = L"\dim\mathcal{K}", ylabel = "patterns",
               title = "achievable-set dimension ($(length(R)) periodic patterns)", titlesize = 11)
    barplot!(ax0, ks .- 0.18, [get(c, k, 0) for k in ks]; width = 0.36, color = MPL_CYCLE[1],
             label = L"measured $\dim\mathcal{K}$")
    barplot!(ax0, ks .+ 0.18, [get(ub, k, 0) for k in ks]; width = 0.36, color = MPL_CYCLE[2],
             label = L"bound $\min(4, 2\,\dim\,\mathrm{null})$")
    axislegend(ax0; labelsize = 9)
    ax1 = Axis(fig[1, 2]; xlabel = L"2\,\mathrm{rank}(D)", ylabel = L"\dim\mathcal{K}",
               title = L"$\dim\mathcal{K} = 2\,\mathrm{rank}(D)$ on every pattern", titlesize = 11)
    scatter!(ax1, [2 * inum(r, "rankD") for r in R], [inum(r, "dimK") for r in R];
             markersize = 10, color = (MPL_CYCLE[1], 0.6))
    lines!(ax1, [0, 4], [0, 4]; color = :black, linestyle = :dash, linewidth = 1)
    savefig(fig, outpath(joinpath(D, "k7_dimK_hist.png")))
end

# Where the hole-area harmonic and forward kinematics part company: only on the theta-interval
# where the deployed cell has inverted (det P_theta < 0), because the measured area uses |det P|
# while the harmonic is the signed det P_0 (det J - 1).
function c4_figure()
    path = joinpath(D, "k7_c4_sweep.csv")
    isfile(path) || return
    R = rows(path)
    names = unique([r["name"] for r in R])
    fig = Figure(size = (360 * length(names), 320))
    for (i, name) in enumerate(names)
        rr = [r for r in R if r["name"] == name]
        th = [fnum(r, "theta") for r in rr]
        dev = [max(fnum(r, "rel_dev"), 1e-17) for r in rr]
        det = [fnum(r, "detP") for r in rr]
        ax = Axis(fig[1, i]; yscale = log10, title = name, titlesize = 11, xlabel = L"\theta",
                  ylabel = i == 1 ? L"relative $|A_{\mathrm{fk}} - A_{\mathrm{harmonic}}|$" : "")
        # shade where det P_theta < 0
        for j in 1:length(th)-1
            det[j] < 0 && vspan!(ax, th[j], th[j+1]; color = (MPL_CYCLE[1], 0.18))
        end
        lines!(ax, th, dev; linewidth = 1.4, color = MPL_CYCLE[4])
    end
    Label(fig[0, 1:length(names)],
          L"K7 C4: the harmonic is exact except where $\det P_\theta < 0$ (shaded), an artefact of measuring $|\det P_\theta|$";
          fontsize = 11)
    savefig(fig, outpath(joinpath(D, "k7_c4_area_fit.png")))
end

# C3 outcome per pattern: the target is always hit to machine precision; what fails is
# certifying a positive deployment range at the hitting design. Prefers the post-F32
# certificate run.
function c3_figure()
    path = joinpath(D, "k7_c3_all_v2.csv")
    isfile(path) || (path = joinpath(D, "k7_c3_all.csv"))
    isfile(path) || return
    R = rows(path)
    fams = unique([r["family"] for r in R])
    fig = Figure(size = (940, 380))
    for (i, tgt) in enumerate(["random_in_K", "diag_1_-0.5"])
        tot = [count(r -> r["family"] == f && r["target"] == tgt, R) for f in fams]
        cer = [count(r -> r["family"] == f && r["target"] == tgt && r["cert"] == "1", R) for f in fams]
        fea = [count(r -> r["family"] == f && r["target"] == tgt && r["zp_feasible"] == "1", R) for f in fams]
        x = 0:length(fams)-1
        ax = Axis(fig[1, i]; xticks = (x, fams), xticklabelrotation = 35 * pi / 180,
                  xticklabelsize = 8, title = "target = $tgt", titlesize = 11, ylabel = "patterns")
        barplot!(ax, x .- 0.26, tot; width = 0.26, color = GREY(0.75), label = "patterns")
        barplot!(ax, x, fea; width = 0.26, color = MPL_CYCLE[1], label = L"$0^+$ feasible")
        barplot!(ax, x .+ 0.26, cer; width = 0.26, color = MPL_CYCLE[3],
                 label = L"certified $\Theta_{\max}>0$")
        i == 1 && axislegend(ax; labelsize = 8)
    end
    Label(fig[0, 1:2], "K7 C3: every target is hit exactly; certification is what fails"; fontsize = 12)
    savefig(fig, outpath(joinpath(D, "k7_c3_certified.png")))
end

nu_figure()
dim_figure()
c4_figure()
c3_figure()
println("wrote ", D)
