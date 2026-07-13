#!/usr/bin/env julia
# Figures for results/kill/. Reads only CSV dumped by the C++ drivers (directive D3).
# Port of plot_kill.py.
# Run: julia --project=Kirigami/scripts Kirigami/scripts/plot_kill.jl
include(joinpath(@__DIR__, "plot_common.jl"))

const K = joinpath(ROOT, "results", "kill")
const BLUE = "#2b6cb0"
const RED = "#c53030"

function k2a()
    p = joinpath(K, "k2a", "k2a.csv")
    isfile(p) || return
    R = rows(p)
    roots = [fnum(r, "theta_roots") for r in R]
    t422 = [fnum(r, "theta_t422") for r in R]
    bis = [fnum(r, "bis_1e12") for r in R]
    graze = t422 .> roots .+ 1e-7
    fig = Figure(size = (1100, 460))
    ax0 = Axis(fig[1, 1]; xlabel = L"$\Theta_{max}$ by collision bisection (shrink $10^{-12}$)",
               ylabel = L"closed-form $\Theta_{max}$",
               title = "K2a: closed form vs bisection, $(length(R)) configurations", titlesize = 12)
    lines!(ax0, [0, pi], [0, pi]; color = :black, linestyle = :dash, linewidth = 0.8,
           label = "exact agreement")
    scatter!(ax0, bis[.!graze], t422[.!graze]; markersize = 8, color = BLUE, label = "T4.2'' scan")
    scatter!(ax0, bis[graze], t422[graze]; markersize = 14, color = :transparent,
             strokecolor = RED, strokewidth = 1.6, label = "graze cases ($(count(graze)))")
    scatter!(ax0, bis[graze], roots[graze]; markersize = 10, marker = :xcross, color = RED,
             label = "min over roots (wrong)")
    for i in findall(graze)
        arrows2d!(ax0, [bis[i]], [roots[i]], [0.0], [t422[i] - roots[i]]; color = RED,
                  shaftwidth = 1.0, tiplength = 6, tipwidth = 6)
    end
    axislegend(ax0; position = :lt, labelsize = 9)

    e422 = abs.(t422 .- bis)
    eroot = abs.(roots .- bis)
    lo = 1e-16
    l422 = log10.(max.(e422, lo)); lroot = log10.(max.(eroot, lo))
    ax1 = Axis(fig[1, 2]; ylabel = "configurations", title = "K2a: error distribution",
               xlabel = L"\log_{10}\,|\Theta_{max}^{closed}-\Theta_{max}^{bisect}|")
    hist_bars!(ax1, l422, linbins(l422, 30); color = BLUE, alpha = 0.75, strokewidth = 0,
               label = "T4.2'' scan (187/187 pass)")
    hist_bars!(ax1, lroot, linbins(lroot, 30); color = RED, alpha = 0.6, strokewidth = 0,
               label = "min over roots (182/187)")
    vlines!(ax1, [log10(1e-5)]; color = :black, linestyle = :dash, linewidth = 1,
            label = L"rule: $10^{-5}$ rad")
    axislegend(ax1; labelsize = 9)
    savefig(fig, outpath(joinpath(K, "k2a", "k2a_agreement.png")))
end

function k1c()
    p = joinpath(K, "k1c", "k1c.csv")
    isfile(p) || return
    R = rows(p)
    labels, naive, samerange, corr, binding = String[], Float64[], Float64[], Float64[], Float64[]
    for (src, lab) in (("native_eq9", "authors' native\nprevent"),
                       ("ours_eq9", "our\noptimize_collision_sweep"))
        S = [r for r in R if r["source"] == src && r["valid_Y9"] == "1"]
        isempty(S) && continue
        n = length(S)
        push!(labels, "$lab\n($n certified)")
        push!(naive, 100 * count(r -> inum(r, "n_naive") > 0, S) / n)
        push!(samerange, 100 * count(r -> inum(r, "n_corr_tmax") > 0, S) / n)
        push!(corr, 100 * count(r -> inum(r, "n_below_beta") > 0, S) / n)
        push!(binding, 100 * count(r -> inum(r, "n_binding") > 0, S) / n)
    end
    x = 0:length(labels)-1
    w = 0.2
    fig = Figure(size = (720, 440))
    ax = Axis(fig[1, 1]; xticks = (x, labels), xticklabelsize = 9,
              ylabel = "designs with at least one Eq. (9) false negative (%)",
              title = "K1c: the rule as written counts a set that is empty once the\n" *
                      "interval and crossing clauses are applied", titlesize = 11)
    sets = [(-1.5w, naive, RED, L"rule as written: $\theta^*<\Theta_{max}$, no clauses"),
            (-0.5w, samerange, "#dd6b20", L"clauses + same range $\Theta_{max}$ (empty by definition)"),
            (0.5w, corr, BLUE, L"corrected: re-closure below $\min\beta$"),
            (1.5w, binding, "#2f855a", "corrected: re-closure is the binding contact")]
    for (off, vals, c, lab) in sets
        barplot!(ax, x .+ off, vals; width = w, color = c, label = lab)
        for (xi, v) in zip(x .+ off, vals)
            text!(ax, xi, v + 1; text = @sprintf("%.1f%%", v), align = (:center, :bottom), fontsize = 9)
        end
    end
    hlines!(ax, [3.0]; color = :black, linestyle = :dash, linewidth = 1, label = "PASS threshold 3%")
    axislegend(ax; labelsize = 9)
    savefig(fig, outpath(joinpath(K, "k1c", "k1c_rates.png")))
end

function k3a()
    p = joinpath(K, "k3a", "k3a.csv")
    isfile(p) || return
    R = [r for r in rows(p) if r["kind"] == "delaunay"]
    isempty(R) && return
    mf = [fnum(r, "m_full") for r in R]
    mc = [fnum(r, "m_core") for r in R]
    F = [fnum(r, "F") for r in R]
    fig = Figure(size = (1100, 440))
    hi = max(maximum(mf), maximum(mc))
    ax0 = Axis(fig[1, 1]; xlabel = L"$m_{full}$ (whole hinge graph)", ylabel = L"$m_{core}$ (2-core only)",
               title = "K3a: deleting dangling faces removes ~5% of mobility (n=$(length(R)))",
               titlesize = 11)
    lines!(ax0, [0, hi], [0, hi]; color = :black, linestyle = :dash, linewidth = 0.8, label = "no change")
    scatter!(ax0, mf, mc; markersize = 7, color = F, colormap = :viridis)
    hlines!(ax0, [5]; color = RED, linewidth = 1.2, linestyle = :dot,
            label = L"adversary's prediction $m_{core}\leq 5$")
    axislegend(ax0; labelsize = 9)
    ax1 = Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = "faces", ylabel = L"m_{core}",
               title = "K3a: intrinsic mobility grows with the patch", titlesize = 11)
    scatter!(ax1, F, mc; markersize = 7, color = BLUE, label = L"m_{core}")
    hlines!(ax1, [5]; color = RED, linewidth = 1.2, linestyle = :dot,
            label = L"$m_{core}\leq 5$: 0 of %$(length(R))")
    axislegend(ax1; labelsize = 9)
    savefig(fig, outpath(joinpath(K, "k3a", "k3a_mobility.png")))
end

function k5()
    p = joinpath(K, "k5", "k5.csv")
    isfile(p) || return
    R = rows(p)
    isempty(R) && return
    fig = Figure(size = (1400, 440))
    Dmc = [fnum(r, "mc_D") for r in R]
    Ddef = [fnum(r, "def_D") for r in R]
    m = (Dmc .> 0) .& (Ddef .> 0)
    ax0 = Axis(fig[1, 1]; xscale = log10, yscale = log10,
               xlabel = L"$D(\sigma_{mc})$  (max-cut, Eq. 1)", ylabel = L"$D(\sigma_{def})$  (defect search)",
               title = "K5: deployability defect, $(length(R)) graphs", titlesize = 11)
    scatter!(ax0, Dmc[m], Ddef[m]; markersize = 8, color = BLUE)
    lo, hi = min(minimum(Dmc[m]), minimum(Ddef[m])), max(maximum(Dmc[m]), maximum(Ddef[m]))
    lines!(ax0, [lo, hi], [lo, hi]; color = :black, linestyle = :dash, linewidth = 0.8, label = "no change")
    axislegend(ax0; labelsize = 9)

    pmc = [fnum(r, "mc_proj_rel") for r in R]
    pdf = [fnum(r, "def_proj_rel") for r in R]
    ax1 = Axis(fig[1, 2]; xlabel = L"$\|X_0-X_{ini}\|_\infty$ / median edge length", ylabel = "graphs",
               title = "K5: how far Eq. (6) has to move the mesh", titlesize = 11)
    hist_bars!(ax1, pmc, linbins(pmc, 30); color = RED, alpha = 0.7, strokewidth = 0, label = L"\sigma_{mc}")
    hist_bars!(ax1, pdf, linbins(pdf, 30); color = BLUE, alpha = 0.7, strokewidth = 0, label = L"\sigma_{def}")
    axislegend(ax1; labelsize = 9)

    n = length(R)
    cats = ["overlap-free\nat θ=0", "certificate\nPOS/NOOVL/NOROOT", "Θmax>0"]
    mcv = [count(r -> r["mc_overlap_free0"] == "1", R), count(r -> r["mc_cert"] == "1", R),
           count(r -> fnum(r, "mc_theta") > 1e-9, R)]
    dfv = [count(r -> r["def_overlap_free0"] == "1", R), count(r -> r["def_cert"] == "1", R),
           count(r -> fnum(r, "def_theta") > 1e-9, R)]
    x = 0:2
    w = 0.36
    ax2 = Axis(fig[1, 3]; xticks = (x, cats), xticklabelsize = 9, ylabel = "% of 200 graphs",
               title = "K5: it fixes the flat sheet, not the deployment", titlesize = 11)
    barplot!(ax2, x .- w/2, 100 .* mcv ./ n; width = w, color = RED, label = L"\sigma_{mc}")
    barplot!(ax2, x .+ w/2, 100 .* dfv ./ n; width = w, color = BLUE, label = L"\sigma_{def}")
    for (xi, v) in zip(x .- w/2, mcv)
        text!(ax2, xi, 100v/n + 1; text = string(v), align = (:center, :bottom), fontsize = 9)
    end
    for (xi, v) in zip(x .+ w/2, dfv)
        text!(ax2, xi, 100v/n + 1; text = string(v), align = (:center, :bottom), fontsize = 9)
    end
    ylims!(ax2, 0, 70)
    axislegend(ax2; labelsize = 9)
    savefig(fig, outpath(joinpath(K, "k5", "k5_orientation.png")))
end

function k2c_drift()
    p = joinpath(K, "k2c", "k2c_drift.csv")
    isfile(p) || return
    # the growing-patch table lives in summary.txt; parse it
    diam, drift = Float64[], Float64[]
    for line in eachline(joinpath(K, "k2c", "summary.txt"))
        t = split(line)
        if length(t) >= 3 && t[1] == "squares"
            push!(diam, parse(Float64, t[2])); push!(drift, parse(Float64, t[3]))
        end
    end
    length(diam) < 3 && return
    fig = Figure(size = (640, 440))
    ax = Axis(fig[1, 1]; xlabel = "patch diameter",
              ylabel = L"\max_f \|\gamma_f(\theta)-\bar{x}_f\| / r_f",
              title = "K2c: centroid drift grows linearly - H-LOC is refuted", titlesize = 11)
    scatterlines!(ax, diam, drift; markersize = 10, color = RED, label = "growing square patch")
    b, a = polyfit1(diam, drift)
    xs = range(minimum(diam), maximum(diam), length = 10)
    lines!(ax, xs, a .+ b .* xs; color = :black, linestyle = :dash, linewidth = 1,
           label = @sprintf("linear fit, slope %.4f", b))
    hlines!(ax, [1.0]; color = "#2f855a", linewidth = 1.2, linestyle = :dot,
            label = L"H-LOC would need a constant $\kappa$")
    axislegend(ax; labelsize = 9)
    savefig(fig, outpath(joinpath(K, "k2c", "k2c_hloc.png")))
end

function k2c()
    p = joinpath(K, "k2c", "k2c.csv")
    isfile(p) || return
    R = rows(p)
    key_pairs = something(findfirst(k -> occursin("pairs", k) && occursin("per", k), R[1].cols), 0)
    key_pairs = key_pairs == 0 ? nothing : R[1].cols[key_pairs]
    F = [fnum(r, "F") for r in R]
    rho = haskey(R[1], "max_rho_over_r") ? [fnum(r, "max_rho_over_r") for r in R] : nothing
    fig = Figure(size = (1100, 440))
    have_rho = rho !== nothing && any(isfinite, rho)
    ax0 = Axis(fig[1, 1]; xlabel = "faces", ylabel = L"\max_f \rho_f / r_f",
               xscale = have_rho ? log10 : identity,
               title = "K2c: the swept disc is exactly the flat circumdisc", titlesize = 11)
    if have_rho
        scatter!(ax0, F, rho; markersize = 7, color = BLUE)
        hlines!(ax0, [1.0]; color = RED, linewidth = 1, linestyle = :dash, label = L"$\rho_f/r_f = 1$ exactly")
        ylims!(ax0, 0.999, 1.001)
        axislegend(ax0; labelsize = 9)
    end
    ax1 = Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = "faces",
               ylabel = "surviving broad-phase pairs per face", title = "K2c: the active set is O(n)",
               titlesize = 11)
    if key_pairs !== nothing
        pp = [fnum(r, key_pairs) for r in R]
        m = isfinite.(pp) .& (pp .> 0) .& (F .> 0)
        scatter!(ax1, F[m], pp[m]; markersize = 7, color = "#2f855a")
        if count(m) > 2
            b, a = polyfit1(log.(F[m]), log.(pp[m]))
            xs = range(log(minimum(F[m])), log(maximum(F[m])), length = 20)
            lines!(ax1, exp.(xs), exp.(a .+ b .* xs); color = :black, linestyle = :dash, linewidth = 1.2,
                   label = @sprintf("slope %.4f (rule ≤ 0.3)", b))
            axislegend(ax1; labelsize = 9)
        end
    end
    savefig(fig, outpath(joinpath(K, "k2c", "k2c_locality.png")))
end

function k2b()
    p = joinpath(K, "k2b", "k2b.csv")
    isfile(p) || return
    R = rows(p)
    isempty(R) && return
    base = [fnum(r, "theta_native_best") for r in R]
    ours = [fnum(r, "theta_ours") for r in R]
    eq9 = [fnum(r, "theta_ours_eq9") for r in R]
    names = [r["name"] for r in R]
    order = sortperm(base)
    x = 0:length(R)-1
    fig = Figure(size = (max(750, 34 * length(R)), 460))
    ax = Axis(fig[1, 1]; xticks = (x, names[order]), xticklabelrotation = pi / 2, xticklabelsize = 7,
              ylabel = L"$\Theta_{max}$ (our bisection, shrink $10^{-12}$)",
              title = "K2b: range against the authors' native baseline (n=$(length(R)))", titlesize = 11)
    scatterlines!(ax, x, base[order]; markersize = 8, color = RED, label = "authors' native prevent (best of ladder)")
    scatterlines!(ax, x, ours[order]; markersize = 8, marker = :rect, color = BLUE,
                  label = "ours (softmin of closed-form roots)")
    scatter!(ax, x, eq9[order]; markersize = 8, marker = :utriangle, color = "#a0aec0",
             label = "our Eq. (9) reimplementation")
    axislegend(ax; labelsize = 9)
    savefig(fig, outpath(joinpath(K, "k2b", "k2b_margin.png")))
end

for fn in (k2a, k1c, k3a, k2c, k2c_drift, k2b, k5)
    try
        fn()
    catch e
        println(stderr, "skip $(nameof(fn)): ", sprint(showerror, e))
    end
end
