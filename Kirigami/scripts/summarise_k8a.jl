#!/usr/bin/env julia
# Summarise results/kill/k8a/k8a.csv (written by Kirigami/apps/kill_k8a.jl).
# A text reduction of the merged CSV, printed to stdout so the
# caller can tee it into summary.txt.
# Run: julia --project=Kirigami/scripts Kirigami/scripts/summarise_k8a.jl results/kill/k8a
using CSV, DataFrames, Printf
include(joinpath(@__DIR__, "plot_common.jl"))

d = length(ARGS) >= 1 ? ARGS[1] : "results/kill/k8a"
R = rows(d * "/k8a.csv")
f(r, k) = fnum(r, k)
med(v) = sorted_mid(v)
counter(v) = (c = Dict{eltype(v),Int}(); for x in v; c[x] = get(c, x, 0) + 1; end; c)
fmtcounter(c) = "{" * join(["$(repr(k)): $(v)" for (k, v) in sort(collect(c); by = first)], ", ") * "}"

println("K8a -- expansive-cone LP")
println("configurations in k8a.csv: ", length(R))

ref = [r for r in R if r["kind"] == "reference"]
dep = [r for r in R if startswith(r["kind"], "deployable_")]

println("\n-- HARD SOUNDNESS CONTROL: the four split-bearing authored tilings")
for r in ref
    @printf("   %-22s F=%3s split=%3s dim_ker_A=%3s dim_flex=%3s  sigma_in_cone=%s  min q=%s  min mu=%s  LP margin=%s  euler_eps=%s\n",
            r["id"], r["F"], r["n_split"], r["dim_ker_A"], r["dim_flex"], r["sigma_in_cone"],
            r["sigma_min_q"], r["sigma_min_mu"], r["margin_l2"], r["euler_eps"])
end
@printf("   sigma in P(X): %d/%d ; LP margin > 0: %d/%d\n",
        sum(inum(r, "sigma_in_cone") for r in ref; init = 0), length(ref),
        count(r -> f(r, "margin_l2") > 1e-8, ref), length(ref))

println("\n-- DEPLOYABLE POPULATION (kill_common::deployable_population)")
@printf("   configs %d ; sigma in P(X) %d ; LP margin > 0 %d ; Euler step tested %d, collision-free %d\n",
        length(dep), sum(inum(r, "sigma_in_cone") for r in dep; init = 0),
        count(r -> f(r, "margin_l2") > 1e-8, dep),
        sum(inum(r, "euler_tested") for r in dep; init = 0),
        count(r -> f(r, "euler_eps") > 0, dep))
println("   euler_eps histogram (median-edge units): ",
        fmtcounter(counter([r["euler_eps"] for r in dep if inum(r, "euler_tested") != 0])))
for r in dep
    if inum(r, "sigma_in_cone") == 0
        @printf("   sigma OUTSIDE P(X): %s  min q=%s  min mu=%s  LP margin=%s\n",
                r["id"], r["sigma_min_q"], r["sigma_min_mu"], r["margin_l2"])
    end
end

for w in ("X_ini", "X0")
    S = [r for r in R if r["where"] == w]
    isempty(S) && continue
    println("\n-- K1a POPULATION at $w  (n=$(length(S)))")
    @printf("   LP margin > 0                                     : %d/%d   (PASS bar: >= 20/100)\n",
            count(r -> f(r, "margin_l2") > 1e-8, S), length(S))
    for t in (1e-6, 1e-5)
        @printf("   dual_max < %g (every branch chart tried certified) : %d/%d\n",
                t, count(r -> f(r, "dual_max") < t, S), length(S))
    end
    @printf("   dual_max: median %.2e  max %.2e\n", med([f(r, "dual_max") for r in S]),
            maximum(f(r, "dual_max") for r in S))
    @printf("   sigma is a flex of the framework : %d/%d\n",
            sum(inum(r, "sigma_is_flex") for r in S), length(S))
    @printf("   sigma in P(X)                    : %d/%d\n",
            sum(inum(r, "sigma_in_cone") for r in S), length(S))
    @printf("   dim ker A median %d, dim flex median %d, c(Gamma) values %s\n",
            med([inum(r, "dim_ker_A") for r in S]), med([inum(r, "dim_flex") for r in S]),
            string(sort(unique(inum(r, "components") for r in S))))
    @printf("   rows median %d (split %d, corner incidences %d, all convex: %s)\n",
            med([inum(r, "n_rows") for r in S]), med([inum(r, "n_split") for r in S]),
            med([inum(r, "n_corner") for r in S]),
            sum(inum(r, "n_convex") for r in S) == sum(inum(r, "n_corner") for r in S) ? "True" : "False")
    for kind in ("voronoi", "delaunay", "quad_random")
        T = [r for r in S if r["kind"] == kind]
        isempty(T) && continue
        @printf("     %-12s n=%3d dim_ker_A med %4d  margin>0 %3d  dual_max<1e-6 %3d  dual_max med %.1e\n",
                kind, length(T), med([inum(r, "dim_ker_A") for r in T]),
                count(r -> f(r, "margin_l2") > 1e-8, T), count(r -> f(r, "dual_max") < 1e-6, T),
                med([f(r, "dual_max") for r in T]))
    end
end

println("\n-- SELF-CHECKS")
@printf("   max |R N| / max|R| over every configuration : %.2e\n",
        maximum(f(r, "flex_resid") for r in R))
println("   passes solved per configuration            : ",
        fmtcounter(counter([inum(r, "passes") for r in R])))
