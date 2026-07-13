#!/usr/bin/env julia
# Final paper figures (port of plot_final.py): baseline (Eq. (6), K6) vs. the three
# constrained-embedding arms (K9 proximity, K9b, K9c range-maximising), plus E1's
# exact-vs-referee agreement. See plot_final.py's header for the source table.
# Figures -> results/final/figures/ (with the _julia suffix while both versions coexist):
#   fig_yield, fig_eps_hist, fig_feasible_vs_deployable, fig_gallery_full, fig_e1_agreement,
#   summary_table.md, README.md
# Run: julia --project=Kirigami/scripts Kirigami/scripts/plot_final.jl
include(joinpath(@__DIR__, "plot_common.jl"))

const K9C_DIR = "results/kill/k9c"
const GALLERY = joinpath(K9C_DIR, "gallery")
const K6_CSV = "results/kill/k6/k6_final.csv"
const E1_CSV = "results/final/e1/e1.csv"
const OUT_DIR = "results/final/figures"
mkpath(OUT_DIR)

const FAMILIES = ["delaunay", "voronoi", "quad_random"]
const FAMILY_LABEL = Dict("delaunay" => "Delaunay", "voronoi" => "Voronoi", "quad_random" => "quad")
const SIGMAS = ["sigma_mc", "sigma_def"]
const ARMS = ["baseline", "k9", "k9b", "k9c"]
const ARM_LABEL = Dict("baseline" => "Eq.(6)+repairs\n(K6)", "k9" => "K9 proximity",
                       "k9b" => "K9b", "k9c" => "K9c range-max")
const ARM_COLOR = Dict("baseline" => GREY(0.6), "k9" => parse(Makie.Colors.Colorant, "#9BB7F0"),
                       "k9b" => parse(Makie.Colors.Colorant, "#5B8FF9"),
                       "k9c" => parse(Makie.Colors.Colorant, "#3363C7"))

# Merge the 12 live shards by (kind, id, sigma), last-write-wins.
function load_k9c_merged()
    out = Dict{Tuple{String,String,String},Row}()
    for f in sort(filter(x -> startswith(basename(x), "shard_") && endswith(x, ".csv"),
                         readdir(K9C_DIR; join = true)))
        for r in rows(f)
            out[(r["kind"], r["id"], r["sigma"])] = r
        end
    end
    return collect(values(out))
end

load_k6() = Dict((r["kind"], r["id"], r["sigma"]) => r for r in rows(K6_CSV))
f(x) = parse(Float64, x)

# ---------------------------------------------------------------------------
# fig_yield: yield vs |F|, per family x sigma, 4 arms, Wilson 95% CI
# ---------------------------------------------------------------------------
function fig_yield(R, k6)
    edges = [100, 230, 420, 800]
    bin_labels = ["F<230\n(F=$(edges[1])-229)", "230<=F<420", "F>=420\n(<=$(edges[end]))"]
    function bin_of(F)
        for i in 1:length(edges)-1
            edges[i] <= F < edges[i+1] && return i
        end
        return length(edges) - 1
    end
    get_(row, arm, field) = arm == "baseline" ? 0.0 : f(row["$(arm)_$(field)"])

    fig = Figure(size = (1350, 760))
    metrics = [("theta", 1e-9, "yield: exact Θmax > 0"),
               ("eps", 0.1, "yield: certified εmax ≥ 0.1 rad")]
    handles, labels = [], String[]
    for (mrow, (field, thresh, mtitle)) in enumerate(metrics), (fcol, fam) in enumerate(FAMILIES)
        ax = Axis(fig[mrow, fcol]; xticks = (0:length(bin_labels)-1, bin_labels), xticklabelsize = 9,
                  title = "$(FAMILY_LABEL[fam]) -- $mtitle", titlesize = 11,
                  ylabel = fcol == 1 ? "fraction of designs (Wilson 95% CI)" : "")
        fam_rows = [r for r in R if r["kind"] == fam]
        n_series = length(ARMS) * length(SIGMAS)
        offsets = [(-0.5 + (j - 0.5) / n_series) * 0.75 for j in 1:n_series]
        j = 0
        for arm in ARMS, sig in SIGMAS
            j += 1
            sub = [r for r in fam_rows if r["sigma"] == sig]
            binned = [Row[] for _ in 1:length(edges)-1]
            for r in sub
                push!(binned[bin_of(inum(r, "F"))], r)
            end
            xs, ys, lo, hi = Float64[], Float64[], Float64[], Float64[]
            for (b, rs) in enumerate(binned)
                isempty(rs) && continue
                k = field == "theta" ? count(r -> get_(r, arm, field) > thresh, rs) :
                                       count(r -> get_(r, arm, field) >= thresh, rs)
                p, lo_, hi_ = wilson(k, length(rs))
                push!(xs, b - 1 + offsets[j]); push!(ys, p)
                push!(lo, max(0.0, p - lo_)); push!(hi, max(0.0, hi_ - p))
            end
            ls = sig == "sigma_mc" ? :solid : :dash
            marker = sig == "sigma_mc" ? :circle : :rect
            errorbars!(ax, xs, ys, lo, hi; color = ARM_COLOR[arm], whiskerwidth = 4, linewidth = 0.8)
            h = scatterlines!(ax, xs, ys; marker, markersize = 7, color = ARM_COLOR[arm],
                              linestyle = ls, linewidth = 0.9)
            if fcol == 1 && mrow == 1
                push!(handles, h)
                push!(labels, "$(split(ARM_LABEL[arm], "\n")[1]) $(sig == "sigma_mc" ? "mc" : "def")")
            end
        end
        ylims!(ax, -0.03, 1.03)
    end
    Legend(fig[3, 1:3], handles, labels; orientation = :horizontal, nbanks = 1,
           framevisible = false, labelsize = 9)
    Label(fig[0, 1:3], "Yield vs. face count |F|, by graph family and orientation rule " *
                       "(N=$(length(R)) designs, $(length(R) ÷ 2) graphs x 2 sigma)"; fontsize = 13)
    savefig(fig, outpath(joinpath(OUT_DIR, "fig_yield.png")))
end

# ---------------------------------------------------------------------------
# fig_eps_hist: eps_max distribution, K9c vs K9, log-x, zero-count bars per arm
# ---------------------------------------------------------------------------
function fig_eps_hist(R)
    fig = Figure(size = (1120, 420))
    k9c_pos = [f(r["k9c_eps"]) for r in R if f(r["k9c_eps"]) > 0]
    k9_pos = [f(r["k9_eps"]) for r in R if f(r["k9_eps"]) > 0]
    lo = log10(min(minimum(k9c_pos), minimum(k9_pos)))
    bins = [10.0^(lo + i * (0 - lo) / 28) for i in 0:28]
    ax0 = Axis(fig[1, 1]; xscale = log10, ylabel = "designs",
               xlabel = L"certified $\varepsilon_{max}$ (rad), designs with $\varepsilon_{max}>0$",
               title = L"$\varepsilon_{max}$ distribution (log-x), K9c vs. K9", titlesize = 12)
    hist_bars!(ax0, k9c_pos, bins; color = ARM_COLOR["k9c"], alpha = 0.75, strokecolor = GREY(0.25),
               strokewidth = 0.4, label = "K9c (n=$(length(k9c_pos)))")
    hist_bars!(ax0, k9_pos, bins; color = ARM_COLOR["k9"], alpha = 0.75, strokecolor = GREY(0.25),
               strokewidth = 0.4, label = "K9 proximity (n=$(length(k9_pos)))")
    vlines!(ax0, [0.1]; color = "#E8684A", linestyle = :dash, linewidth = 1.1)
    ymax = max(maximum(bincounts(k9c_pos, bins)), maximum(bincounts(k9_pos, bins)))
    text!(ax0, 0.11, ymax * 0.9; text = "0.1 rad", color = "#E8684A", fontsize = 9)
    axislegend(ax0; labelsize = 9, framevisible = false)

    n = length(R)
    zero_counts = [arm == "baseline" ? n : n - count(r -> f(r["$(arm)_eps"]) > 0, R) for arm in ARMS]
    ax1 = Axis(fig[1, 2]; xticks = (0:3, [ARM_LABEL[a] for a in ARMS]), xticklabelsize = 9,
               ylabel = L"designs with $\varepsilon_{max}=0$ (of %$n)", title = "count at zero, per arm",
               titlesize = 12)
    barplot!(ax1, 0:3, zero_counts; color = [ARM_COLOR[a] for a in ARMS], strokecolor = GREY(0.25),
             strokewidth = 0.5)
    for (i, v) in enumerate(zero_counts)
        text!(ax1, i - 1, v + n * 0.01; text = string(v), align = (:center, :bottom), fontsize = 9)
    end
    savefig(fig, outpath(joinpath(OUT_DIR, "fig_eps_hist.png")))
end

# ---------------------------------------------------------------------------
# fig_feasible_vs_deployable: best 0+ margin vs exact Theta_max, all K9c rows
# ---------------------------------------------------------------------------
function fig_feasible_vs_deployable(R)
    fig = Figure(size = (740, 600))
    ax = Axis(fig[1, 1];
              xlabel = L"best 0$^+$ margin $m$ (convex + split-feasibility slack; $m>0$ = feasible)",
              ylabel = L"exact $\Theta_{max}$ (rad)",
              title = "Feasible vs. deployable: K9c's best-of-3 margin against exact Θmax\n(N=$(length(R)) designs)",
              titlesize = 12)
    fam_color = Dict("delaunay" => "#3363C7", "voronoi" => "#E8684A", "quad_random" => "#2CA58D")
    for fam in FAMILIES
        xs = [f(r["best_margin"]) for r in R if r["kind"] == fam]
        ys = [f(r["theta_exact"]) for r in R if r["kind"] == fam]
        scatter!(ax, xs, ys; markersize = 7, color = (fam_color[fam], 0.65), label = FAMILY_LABEL[fam])
    end
    stuck = [(f(r["best_margin"]), f(r["theta_exact"])) for r in R
             if f(r["best_margin"]) > 0 && f(r["theta_exact"]) == 0.0]
    if !isempty(stuck)
        scatter!(ax, first.(stuck), last.(stuck); markersize = 14, color = :transparent,
                 strokecolor = "#E8684A", strokewidth = 1.4,
                 label = L"margin-feasible, $\Theta_{max}=0$ (n=%$(length(stuck)))")
    end
    vlines!(ax, [0]; color = GREY(0.5), linewidth = 0.7, linestyle = :dot)
    axislegend(ax; position = :lt, labelsize = 9, framevisible = false)
    out = outpath(joinpath(OUT_DIR, "fig_feasible_vs_deployable.png"))
    savefig(fig, out)
    println("wrote $out ($(length(stuck)) margin-feasible-but-Theta_max=0 designs)")
end

# ---------------------------------------------------------------------------
# fig_gallery_full: ALL K9c-deployable designs, sorted by eps_max descending
# ---------------------------------------------------------------------------
function draw!(ax, path, color, title)
    d = JSON.parsefile(path)
    V = d["vertices"]
    draw_faces!(ax, V, d["faces"]; color, strokecolor = GREY(0.3), strokewidth = 0.12)
    frame_axis!(ax, V; margin = 0.04)
    ax.title = title
end

function fig_gallery_full(R)
    CLOSED, OPEN = "#BDD2FD", "#5B8FF9"
    tagof(r) = "$(r["kind"])_$(r["id"])_$(r["sigma"])"
    deployable = [r for r in R if f(r["theta_exact"]) > 0 &&
                  isfile(joinpath(GALLERY, tagof(r) * "_closed.json")) &&
                  isfile(joinpath(GALLERY, tagof(r) * "_half.json"))]
    sort!(deployable; by = r -> -f(r["eps_max"]))
    n_missing = count(r -> f(r["theta_exact"]) > 0, R) - length(deployable)

    per_row = 6
    ncol = 2 * per_row
    nrow = (length(deployable) + per_row - 1) ÷ per_row
    fig = Figure(size = (82 * ncol, 95 * nrow + 80))
    for (i, r) in enumerate(deployable)
        row, cpair = divrem(i - 1, per_row)
        tag = tagof(r)
        draw!(Axis(fig[row + 1, 2cpair + 1]; titlesize = 6),
              joinpath(GALLERY, tag * "_closed.json"), CLOSED,
              "$(first(r["kind"], 3))$(r["id"]) $(endswith(r["sigma"], "mc") ? "mc" : "def")\nF=$(r["F"])")
        draw!(Axis(fig[row + 1, 2cpair + 2]; titlesize = 6),
              joinpath(GALLERY, tag * "_half.json"), OPEN, @sprintf("ε=%.2f", f(r["eps_max"])))
    end
    Label(fig[0, 1:ncol],
          "All $(length(deployable)) K9c-deployable designs (of $(length(R)) measured), sorted by certified εmax descending\n" *
          "left of each pair = closed (θ=0); right = deployed to Θmax/2\n" *
          "baseline (Eq.(6)+repairs, K6): closed only, Θmax=0 on all $(length(R)) designs in this population -- no deployed panel exists to show";
          fontsize = 10)
    out = outpath(joinpath(OUT_DIR, "fig_gallery_full.png"))
    savefig(fig, out; px_per_unit = 1.5)
    println("wrote $out with $(length(deployable)) designs ($n_missing deployable rows missing gallery json, excluded)")
end

# ---------------------------------------------------------------------------
# fig_e1_agreement: |exact - referee| histogram + certificate confusion counts
# ---------------------------------------------------------------------------
function fig_e1_agreement()
    R = rows(E1_CSV)
    gaps = [abs(f(r["gap"])) for r in R]
    cert = [r["cert_valid"] == "1" for r in R]
    actual = [r["actually_ge_eps"] == "1" for r in R]
    tp = count(cert .& actual); fp = count(cert .& .!actual)
    fn = count(.!cert .& actual); tn = count(.!cert .& .!actual)

    fig = Figure(size = (1060, 420))
    ax0 = Axis(fig[1, 1]; yscale = log10, ylabel = "count (log)",
               xlabel = L"|\Theta_{max}(\mathrm{T4.2''}) - \Theta_{max}(\mathrm{bisection})|\ \mathrm{(rad)}",
               title = @sprintf("E1: exact vs. bisection referee, N=%d\nworst gap %.3e", length(R), maximum(gaps)),
               titlesize = 12)
    hist_bars!(ax0, gaps, linbins(gaps, 60); color = "#3b6ea5", logy = true)
    vlines!(ax0, [1e-5]; color = :crimson, linestyle = :dash, linewidth = 1, label = "1e-5 agreement bar")
    axislegend(ax0; labelsize = 9)

    labels = ["certified &\nΘmax ≥ ε\n(TP)", "certified but\nΘmax < ε\n(FP, unsound)",
              "not certified,\nΘmax ≥ ε\n(FN, over-strict)", "not certified,\nΘmax < ε\n(TN)"]
    vals = [tp, fp, fn, tn]
    cols = [parse(Makie.Colors.Colorant, "#2CA58D"), parse(Makie.Colors.Colorant, "#E8684A"),
            parse(Makie.Colors.Colorant, "#E8684A"), GREY(0.65)]
    ax1 = Axis(fig[1, 2]; xticks = (0:3, labels), xticklabelsize = 9, ylabel = "designs",
               title = "certificate confusion counts (eps=0.006), N=$(length(R))\n0 unsound / 0 over-strict on this population",
               titlesize = 11)
    barplot!(ax1, 0:3, vals; color = cols, strokecolor = GREY(0.25), strokewidth = 0.5)
    for (i, v) in enumerate(vals)
        text!(ax1, i - 1, v + length(R) * 0.01; text = string(v), align = (:center, :bottom), fontsize = 9)
    end
    out = outpath(joinpath(OUT_DIR, "fig_e1_agreement.png"))
    savefig(fig, out)
    @printf("wrote %s (N=%d, worst gap %.3e, TP=%d FP=%d FN=%d TN=%d)\n", out, length(R), maximum(gaps), tp, fp, fn, tn)
    return length(R), maximum(gaps), tp, fp, fn, tn
end

# ---------------------------------------------------------------------------
# summary_table.md
# ---------------------------------------------------------------------------
function summary_table(R)
    n = length(R)
    complete = n >= 400 ? " (full target)" : " (short of the 400-design target; re-run once more shards land)"
    lines = ["# Headline table: baseline vs. method arms\n",
             "K9c live-merged population: N=$n designs, $(n ÷ 2) unique (kind,id) graphs x 2 sigma rules$complete.", "",
             "| arm | N | deployable (\$\\Theta_{max}>0\$) | certified (\$\\varepsilon_{max}>0\$) " *
             "| \$\\varepsilon_{max}\\geq0.1\$ rad | median \$\\varepsilon_{max}\$ | max \$\\varepsilon_{max}\$ |",
             "|---|---|---|---|---|---|---|"]
    for arm in ARMS
        lab = replace(ARM_LABEL[arm], "\n" => " ")
        if arm == "baseline"
            push!(lines, "| $lab | $n | 0 | 0 | 0 | -- | -- |")
            continue
        end
        dep = count(r -> f(r["$(arm)_theta"]) > 1e-9, R)
        eps_pos = [f(r["$(arm)_eps"]) for r in R if f(r["$(arm)_eps"]) > 0]
        ge01 = count(>=(0.1), eps_pos)
        med = isempty(eps_pos) ? 0.0 : sorted_mid(eps_pos)
        mx = isempty(eps_pos) ? 0.0 : maximum(eps_pos)
        push!(lines, @sprintf("| %s | %d | %d | %d | %d | %.3f | %.3f |", lab, n, dep, length(eps_pos), ge01, med, mx))
    end
    eps_pos = [f(r["eps_max"]) for r in R if f(r["eps_max"]) > 0]
    dep = count(r -> f(r["theta_exact"]) > 1e-9, R)
    ge01 = count(>=(0.1), eps_pos)
    med = isempty(eps_pos) ? 0.0 : sorted_mid(eps_pos)
    mx = isempty(eps_pos) ? 0.0 : maximum(eps_pos)
    push!(lines, @sprintf("| best-of-3 (per-design max over K9/K9b/K9c) | %d | %d | %d | %d | %.3f | %.3f |",
                          n, dep, length(eps_pos), ge01, med, mx))
    push!(lines, "")
    out = outpath(joinpath(OUT_DIR, "summary_table.md"))
    write(out, join(lines, "\n") * "\n")
    println("wrote $out")
end

function main()
    k9c_rows = load_k9c_merged()
    k6 = load_k6()
    println("K9c merged: $(length(k9c_rows)) rows (target 400)")
    fig_yield(k9c_rows, k6)
    fig_eps_hist(k9c_rows)
    fig_feasible_vs_deployable(k9c_rows)
    fig_gallery_full(k9c_rows)
    e1_n, e1_worst, tp, fp, fn, tn = fig_e1_agreement()
    summary_table(k9c_rows)

    n = length(k9c_rows)
    readme = """
# results/final/figures -- baseline vs. method figure set

Data sources and the exact command that produced everything here:
`julia --project=Kirigami/scripts Kirigami/scripts/plot_final.jl` (reads only; does not
re-run any solver).

| figure | what it shows | data source(s) |
|---|---|---|
| `fig_yield.png` | Yield vs. face count \\|F\\| in 3 bins, per family (Delaunay/Voronoi/quad) and per sigma rule, 4 arms (baseline Eq.(6)+repairs, K9 proximity, K9b, K9c range-maximising), Wilson 95% CIs, two rows: exact Theta_max>0 and certified eps_max>=0.1 rad. | `results/kill/k9c/shard_*.csv` merged by (kind,id,sigma), N=$n rows at merge time; baseline confirmed identically 0 from `results/kill/k6/k6_final.csv`. |
| `fig_eps_hist.png` | eps_max distribution (log-x) for K9c vs. K9 proximity, plus a bar of the zero-eps_max count per arm. | same K9c merge. |
| `fig_feasible_vs_deployable.png` | Scatter of best 0+ margin (convex+split feasibility slack) vs. exact Theta_max for every row, colored by family, with margin-feasible-but-Theta_max=0 designs ringed. | same K9c merge. |
| `fig_gallery_full.png` | ALL K9c-deployable designs (Theta_max>0), sorted by certified eps_max descending, shown closed and at Theta_max/2, with a caption noting the baseline has no deployed panel to show (Theta_max=0 throughout). | `results/kill/k9c/gallery/*.json` (per-design geometry dumps) + the same merged CSV for sort order. |
| `fig_e1_agreement.png` | E1's \\|exact - bisection referee\\| histogram (log-y) and the certificate's confusion counts (TP=$tp, FP=$fp, FN=$fn, TN=$tn against Theta_max>=eps=0.006). | `results/final/e1/e1.csv`, N=$e1_n rows. |
| `summary_table.md` | Headline table: arms x {N, deployable, certified, eps_max>=0.1 rad, median eps_max, max eps_max}, computed directly from the CSVs above. | same. |

## Notes

- K9c's `results/kill/k9c/summary.txt` and `k9c.csv` are a stale 182-row snapshot; this
  script always merges the live `shard_*.csv` files instead. At merge time the 12 shards
  covered N=$n of the 400-design target (200 graphs x 2 sigma). If the run had still been
  short of 400, this note would say so and the script would need re-running once more
  shards land -- it reads the shards live every time, never the stale `k9c.csv`.
- The baseline (Eq. (6) + the K6 0+ repair) is independently confirmed 0/400 designs
  deployable and 0/400 certified in `results/kill/k6/k6_final.csv` before it is asserted
  in any figure here.
- `fig_gallery_full.png` intentionally renders every deployable design, not a curated
  subset (review/constructive_review.md section 5, item 1) -- it is a large image, one
  row per 6 designs.
"""
    write(outpath(joinpath(OUT_DIR, "README.md")), readme)
    println("wrote ", outpath(joinpath(OUT_DIR, "README.md")))
end

main()
