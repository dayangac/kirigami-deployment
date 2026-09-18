# plot_common.jl -- helpers shared by the plot_*.jl CairoMakie figure scripts.
#
# CSVs are read as String cells (csv.DictReader semantics) and parsed on use, so every
# per-script tally is explicit about which column it counts and how. `outpath` is the
# single hook through which every script names its figure.
using CairoMakie
using CSV
using DataFrames
using JSON
using Printf
using Statistics

CairoMakie.activate!(type = "png")

const ROOT = dirname(dirname(@__DIR__))  # <repo>/Kirigami/scripts -> <repo>

# ---- CSV as a vector of String=>String dicts (csv.DictReader semantics) --------------
struct Row
    d::Dict{String,String}
    cols::Vector{String}
end
Base.getindex(r::Row, k::AbstractString) = r.d[k]
Base.get(r::Row, k::AbstractString, dflt) = get(r.d, k, dflt)
Base.haskey(r::Row, k::AbstractString) = haskey(r.d, k)

function rows(path::AbstractString)
    df = CSV.read(path, DataFrame; types = String, pool = false)
    cols = String.(names(df))
    out = Row[]
    for i in 1:nrow(df)
        d = Dict{String,String}()
        for c in cols
            v = df[i, c]
            d[c] = v === missing ? "" : String(v)
        end
        push!(out, Row(d, cols))
    end
    return out
end

fnum(r::Row, k, d = NaN) = (v = tryparse(Float64, get(r, k, "")); v === nothing ? d : v)
inum(r::Row, k, d = 0) = (v = tryparse(Int, get(r, k, "")); v === nothing ? d : v)
fnum(s::AbstractString, d = NaN) = (v = tryparse(Float64, s); v === nothing ? d : v)

# ---- output naming ----------------------------------------------------------------------
"canonical output path (identity; kept so every script names its figure through one hook)"
outpath(p::AbstractString) = String(p)

function savefig(fig, p::AbstractString; px_per_unit = 2)
    mkpath(dirname(p))
    save(p, fig; px_per_unit)
    println("wrote ", p)
end

# ---- K6 baseline (Eq. (6) + repairs) --------------------------------------------------
# The headline K6 rule: a design counts as deployable when the best of the four repair variants
# (primary, split-only secondary, lambda ladder, stage 2) has a refereed range > 0.
const K6_BASELINE_CSV = "results/kill/k6/k6.csv"
const K6_THETA_COLS = ["theta_exact", "sec_theta_bisect", "lad_theta_bisect", "s2_theta_bisect"]
k6_best_theta(r::Row) = maximum(fnum(r, c, 0.0) for c in K6_THETA_COLS)
"number of K6 designs deployable under the best-of-four rule (17 of 400 in the committed run)"
k6_baseline_deployable(path = K6_BASELINE_CSV) = isfile(path) ? count(r -> k6_best_theta(r) > 1e-9, rows(path)) : 0

# ---- statistics -----------------------------------------------------------------------
"Wilson 95 % interval -> (p, lo, hi)"
function wilson(k, n; z = 1.959963984540054)
    n == 0 && return (NaN, NaN, NaN)
    p = k / n
    d = 1 + z^2 / n
    c = p + z^2 / (2n)
    r = z * sqrt(p * (1 - p) / n + z^2 / (4n^2))
    return p, max(0.0, min(p, (c - r) / d)), min(1.0, max(p, (c + r) / d))
end

median_sorted(v) = sorted_mid(v)
sorted_mid(v) = (s = sort(v); s[length(s) ÷ 2 + 1])  # python sorted(v)[len(v)//2]

"OLS slope and intercept of y on x (np.polyfit(x, y, 1))"
function polyfit1(x, y)
    mx, my = mean(x), mean(y)
    b = sum((x .- mx) .* (y .- my)) / sum((x .- mx) .^ 2)
    return b, my - b * mx
end

# ---- histogram helpers ----------------------------------------------------------------
"bin counts for explicit edges (numpy.histogram semantics, last bin closed)"
function bincounts(v, edges)
    c = zeros(Int, length(edges) - 1)
    for x in v
        (x < edges[1] || x > edges[end]) && continue
        i = searchsortedlast(edges, x)
        i = min(i, length(edges) - 1)
        c[i] += 1
    end
    return c
end

"n equal-width bins over the data range (matplotlib hist(bins=n) semantics)"
function linbins(v, n)
    lo, hi = isempty(v) ? (0.0, 1.0) : (minimum(v), maximum(v))
    hi == lo && return collect(range(lo - 0.5, hi + 0.5, length = n + 1))
    return collect(range(lo, hi, length = n + 1))
end

"histogram as bars; bars with zero count are skipped so log-y axes work"
function hist_bars!(ax, v, edges; color, alpha = 1.0, strokecolor = :white,
                    strokewidth = 0.3, label = nothing, logy = false)
    c = bincounts(v, edges)
    ctr = [(edges[i] + edges[i+1]) / 2 for i in 1:length(edges)-1]
    w = [edges[i+1] - edges[i] for i in 1:length(edges)-1]
    keep = c .> 0
    any(keep) || return
    barplot!(ax, ctr[keep], c[keep]; width = w[keep], gap = 0,
             color = (color, alpha), strokecolor, strokewidth, label,
             fillto = logy ? 0.5 : 0.0)
end

"""log-x histogram of non-negative gaps: quarter-decade bins from `floor` to `top`; exact zeros
are floored to `floor` so they land in the leftmost bin (the x tick there reads `0 / <=floor`)."""
function gap_hist_log!(ax, gaps; floor = 1e-14, top = 1e-4, color, logy = true)
    g = [max(x, floor) for x in gaps]
    edges = 10 .^ collect(log10(floor):0.25:log10(top))
    hist_bars!(ax, g, edges; color, logy)
    ax.xscale = log10
    ax.xticks = (10.0 .^ (log10(floor):2:log10(top)),
                 [i == 0 ? "0 / ≤1e$(Int(round(log10(floor))))" : "1e$(Int(round(log10(floor) + 2i)))"
                  for i in 0:length(log10(floor):2:log10(top))-1])
    xlims!(ax, floor / 1.5, top * 1.5)
end

# ---- polygon drawing (matplotlib PatchCollection equivalents) ------------------------
"draw faces of a pattern JSON {vertices, faces[, orientation, holes]}"
function draw_faces!(ax, V, faces; color, strokecolor = RGBf(0.3, 0.3, 0.3),
                     strokewidth = 0.4)
    polys = [Point2f[Point2f(V[k+1][1], V[k+1][2]) for k in f] for f in faces if length(f) >= 3]
    isempty(polys) && return
    if color isa AbstractVector
        poly!(ax, polys; color, strokecolor, strokewidth)
    else
        poly!(ax, polys; color, strokecolor, strokewidth)
    end
end

function frame_axis!(ax, V; margin = 0.05)
    xs = [p[1] for p in V]
    ys = [p[2] for p in V]
    m = margin * max(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys), 1e-9)
    limits!(ax, minimum(xs) - m, maximum(xs) + m, minimum(ys) - m, maximum(ys) + m)
    ax.aspect = DataAspect()
    hidedecorations!(ax)
    hidespines!(ax)
end

const GREY(v) = RGBf(v, v, v)
# matplotlib default colour cycle
const MPL_CYCLE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
                   "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
