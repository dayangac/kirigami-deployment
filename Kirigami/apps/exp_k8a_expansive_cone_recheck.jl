# K8a RECHECK -- is the expansive-cone LP's negative answer a property of the cone, or
# of the solver?  (results/experiments/k8a/recheck.md; EXPERIMENTS.md section K8a, "Recheck".)
#
# THE CONTRADICTION THIS RESOLVES. K9 (F36) ran K8a's `cone_lp` at its 30 convexity +
# split-inward constrained embeddings and got 0 feasible, with dual bounds 2e-2 to 5e-2 --
# three to five orders ABOVE the 1e-7 infeasibility threshold, so the LP claimed nothing
# either way. But on 23 of those 30 the UNIFORM ray sigma is inside P(X) BY DIRECT
# MEASUREMENT (min q_e > 0 and min mu > 0 from `zero_plus.jl`), and sigma is a flex
# there because the embedding lies in the Eq. (4) null space. So the cone is non-empty on
# 23 of 30 and a correct solver MUST return a positive margin on every one of them.
#
# THE TEST. For each K9 design, regenerate variant (b) exactly as `exp_k9_convex_embedding` does, then
# compare three numbers that a correct solver has to order:
#
#     min(min_e q_e, min mu) / 2      the geometric margin of the sigma ray
#  <= sigma_chart_margin              the same ray in the LP's own coordinates, on the
#                                     branch chart it selects -- computed with NO solver
#  <= margin_l2                       what the LP returns
#  <= dual_bound                      the rigorous upper bound
#
# The middle two differ from the first only by the row normalisation and the projection
# onto the orthonormal flex basis, both of which rescale but never change a sign. PASS if
# margin_l2 >= sigma_chart_margin on every design with sigma in P(X); the spec's weaker
# bar is margin > 0 with margin >= min(q, mu)/2 in the unnormalised reading.
#
# STEP 4 OF THE TASK. The same run answers "does a NON-UNIFORM flex with a larger margin
# than the uniform ray exist at K9's constrained embeddings?" -- that is exactly
# margin_l2 > sigma_chart_margin, and the gap between them is reported per design.
#
#   julia --project=Kirigami Kirigami/apps/exp_k8a_expansive_cone_recheck.jl [--n 200] [--maxf 800]
#         [--out DIR] [--sigma DIR] [--cache DIR] [--delta 1e-3] [--split-delta 1e-3]
#         [--iters 600] [--starts 3] [--barrier 6] [--dual-iters 20000] [--cone-maxf 800]
#         [--shard S] [--nshards M] [--limit K] [--regenerate]
#
# Population: the K1a graphs (frozen k1a_200 with its archived K5 `sigma_def`; `--sigma DIR`
# reads results/experiments/k5/sigma/<kind>_<id>.json instead, `--regenerate` rebuilds the graphs
# with make_graph). The shape cache is results/experiments/k6/cache/<sigma>/, so `--cache`
# reuses K6's X0 / Phi. Outputs go to results/experiments/k8a/.
include(joinpath(@__DIR__, "exp_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

Base.@kwdef mutable struct Row
    id::Int = 0
    kind::String = ""; sigma::String = ""; variant::String = ""
    N::Int = 0; F::Int = 0; n_split::Int = 0; n_rows::Int = 0; dim_flex::Int = 0
    embed_feasible::Int = 0
    min_q::Float64 = 0.0; min_mu::Float64 = 0.0; geom_margin::Float64 = 0.0  # min(min q, min mu) / 2
    sigma_in_cone::Int = 0; sigma_span::Int = 0
    sigma_chart::Float64 = 0.0; margin::Float64 = 0.0; dual::Float64 = 0.0; gap::Float64 = 0.0
    pass_feasible::Int = -1; passes::Int = 0
    farkas_resid::Float64 = 0.0; farkas_lmin::Float64 = 0.0; farkas_sumerr::Float64 = 0.0
    cert_1e9::Int = 0
    secs::Float64 = 0.0
end

const kHeader =
    "id,kind,sigma,variant,N,F,n_split,n_rows,dim_flex,embed_feasible,min_q,min_mu," *
    "geom_margin,sigma_in_cone,sigma_span,sigma_chart,margin,dual,gap,pass_feasible," *
    "passes,farkas_resid,farkas_lmin,farkas_sumerr,cert_1e9,secs\n"

function write_row(o::IO, r::Row)
    print(o, r.id, ",", r.kind, ",", r.sigma, ",", r.variant, ",", r.N, ",",
          r.F, ",", r.n_split, ",", r.n_rows, ",", r.dim_flex, ",",
          r.embed_feasible, ",", sci(r.min_q), ",", sci(r.min_mu), ",",
          sci(r.geom_margin), ",", r.sigma_in_cone, ",", r.sigma_span, ",",
          sci(r.sigma_chart), ",", sci(r.margin), ",", sci(r.dual), ",",
          sci(r.gap), ",", r.pass_feasible, ",", r.passes, ",",
          sci(r.farkas_resid), ",", sci(r.farkas_lmin), ",", sci(r.farkas_sumerr),
          ",", r.cert_1e9, ",", fx(r.secs, 2), "\n")
end

function evaluate(id::Int, kind::String, sname::String, variant::String, m::K.Mesh,
                  c::K.CutStructure, X::Vector{Vec2}, embed_feasible::Bool,
                  opt::K.ExpansiveConeOptions)
    t = Timer()
    r = Row(id = id, kind = kind, sigma = sname, variant = variant,
            N = K.n_vertices(m), F = K.n_faces(m), embed_feasible = embed_feasible ? 1 : 0)

    q = K.zero_plus_q(c, X)
    mu = K.zero_plus_corner_margin(c, X)
    r.min_q = isempty(q) ? 0.0 : minimum(q)
    r.min_mu = isempty(mu) ? 0.0 : minimum(mu)
    r.geom_margin = 0.5 * min(r.min_q, r.min_mu)

    rep = K.expansive_cone(c, X, opt)
    r.n_split = K.n_split(c)
    r.n_rows = rep.n_rows
    r.dim_flex = rep.dim_flex
    r.sigma_in_cone = rep.sigma_in_cone ? 1 : 0
    r.sigma_span = rep.sigma_in_flex_span ? 1 : 0
    r.sigma_chart = rep.sigma_chart_margin
    r.margin = rep.margin_l2
    r.dual = rep.dual_bound_max
    r.gap = rep.lp_gap_max
    r.pass_feasible = rep.pass_feasible
    r.passes = rep.passes
    r.farkas_resid = rep.farkas_resid_max
    r.farkas_lmin = rep.farkas_lambda_min
    r.farkas_sumerr = rep.farkas_sum_err
    r.cert_1e9 = (rep.passes > 0 && rep.passes_certified_1e9 == rep.passes) ? 1 : 0
    r.secs = s(t)
    return r
end

function median_of(v::Vector{Float64})
    isempty(v) && return NaN
    sv = sort(v)
    return sv[length(sv) ÷ 2 + 1]
end

# sigma_def of one graph: the archived K5 file if `sigmadir` is given, else the frozen row
function sigma_def_of(row::PopRow, sigmadir::AbstractString)
    if !isempty(sigmadir)
        p = joinpath(sigmadir, row.kind * "_" * string(row.id) * ".json")
        isfile(p) || return Int[]
        sm = K.load_mesh_json(p)
        (length(sm.sigma) == K.n_faces(row.mesh) && length(sm.X) == length(row.mesh.X)) || return Int[]
        return sm.sigma
    end
    return row.sigma_def === nothing ? Int[] : row.sigma_def
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n = 200; maxf = 800; max_iter = 600; n_random = 3; barrier_stages = 6
    delta_rel = 1e-3; split_delta_rel = 1e-3
    dual_iters = 20000; cone_maxf = 800
    shard = 0; nshards = 1
    limit = typemax(Int)
    outdir = joinpath(REPO, "results", "experiments", "k8a")
    sigmadir = ""
    cache = joinpath(REPO, "results", "experiments", "k6", "cache")
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--maxf" && i < length(args); maxf = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--sigma" && i < length(args); sigmadir = args[i+1]; i += 2
        elseif a == "--cache" && i < length(args); cache = args[i+1]; i += 2
        elseif a == "--delta" && i < length(args); delta_rel = arg_f(args[i+1]); i += 2
        elseif a == "--split-delta" && i < length(args); split_delta_rel = arg_f(args[i+1]); i += 2
        elseif a == "--iters" && i < length(args); max_iter = arg_i(args[i+1]); i += 2
        elseif a == "--starts" && i < length(args); n_random = arg_i(args[i+1]); i += 2
        elseif a == "--barrier" && i < length(args); barrier_stages = arg_i(args[i+1]); i += 2
        elseif a == "--dual-iters" && i < length(args); dual_iters = arg_i(args[i+1]); i += 2
        elseif a == "--cone-maxf" && i < length(args); cone_maxf = arg_i(args[i+1]); i += 2
        elseif a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshards" && i < length(args); nshards = arg_i(args[i+1]); i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    sfx = nshards > 1 ? "_" * string(shard) : ""
    csv = open(joinpath(outdir, "recheck_k9" * sfx * ".csv"), "w")
    print(csv, kHeader)
    flush(csv)

    opt = K.ExpansiveConeOptions()
    opt.lp.dual_iters = dual_iters

    rows = Row[]
    wall = Timer()
    gidx = -1; missing_sigma = 0
    nrun = 0
    # the K1a population = the first graphs of make_graph(id, 100, maxf, 1400) that build;
    # the frozen k1a_200 carries maxf == 800 only
    maxf == 800 || regenerate || @warn "frozen k1a_200 was built with --maxf 800; use --regenerate for maxf = $maxf"
    for row in load_population("k1a_200")
        gidx + 1 < n || break
        row.ok || continue
        if regenerate
            g = K.make_graph(row.id, 100, maxf, 1400)
            g.ok || continue
            row.mesh = g.mesh
        end
        m0 = row.mesh
        K.build_topology!(m0)
        sigma_mc = copy(m0.sigma)
        X_ini = copy(m0.X)
        id = row.id

        sigma_def = sigma_def_of(row, sigmadir)
        if isempty(sigma_def)
            missing_sigma += 1
            continue
        end
        gidx += 1
        gidx >= n && break
        (nshards > 1 && gidx % nshards != shard) && continue
        nrun >= limit && break
        nrun += 1

        med = median_edge_length(m0)

        for which in 0:1
            fname = which == 1 ? "sigma_def" : "sigma_mc"
            m = K._with_sigma(m0, which == 1 ? sigma_def : sigma_mc)
            K.build_topology!(m)
            c = K.make_cut(m)
            K.n_split(c) == 0 && continue
            hs = K.holes_partition(c)
            sh = shape_space(m, c, hs, joinpath(cache, fname), id)
            (!sh.ok || sh.k < 1) && continue
            K.n_faces(m) > cone_maxf && continue
            X0 = K.matrix_to_points(sh.X0)

            # Variant (b), regenerated exactly as exp_k9_convex_embedding does it (same seeds, same three
            # starts, same fallbacks) so that the embedding is bit-identical to K9's.
            oa = K.ConvexEmbedOptions()
            oa.delta_rel = delta_rel
            oa.n_random = n_random
            oa.max_iter = max_iter
            oa.barrier_stages = barrier_stages
            oa.seed = (UInt32(9000) + UInt32(7) * UInt32(id) + UInt32(which)) % UInt32
            ra = K.convex_embed(c, X0, sh.Phi, X_ini, med, oa)

            ro = K.ZeroPlusRepairOptions()
            ro.n_random = n_random
            ro.max_iter = max_iter
            ro.seed = (UInt32(9000) + UInt32(7) * UInt32(id) + UInt32(which)) % UInt32
            sp = K.zero_plus_repair(c, X0, sh.Phi, med, ro)

            ob = deepcopy(oa)
            ob.split_delta_rel = split_delta_rel
            rb = K.convex_embed(c, X0, sh.Phi, X_ini, med, ob)
            if !rb.feasible && ra.feasible
                o2 = deepcopy(ob)
                o2.t_init = ra.t
                r2 = K.convex_embed(c, X0, sh.Phi, X_ini, med, o2)
                (r2.feasible || r2.n_bad + r2.n_bad_q < rb.n_bad + rb.n_bad_q) && (rb = r2)
            end
            if !rb.feasible && sp.feasible
                o3 = deepcopy(ob)
                o3.t_init = sp.t
                r3 = K.convex_embed(c, X0, sh.Phi, X_ini, med, o3)
                (r3.feasible || r3.n_bad + r3.n_bad_q < rb.n_bad + rb.n_bad_q) && (rb = r3)
            end
            rb.feasible || continue  # K9 only ran the cone on the feasible embeddings

            r = evaluate(id, row.kind, fname, "b", m, c, rb.X, rb.feasible, opt)
            write_row(csv, r)
            flush(csv)
            push!(rows, r)
            println("  ", id, " ", fname, " F=", r.F, " sigma_in_cone=",
                    r.sigma_in_cone, " geom=", sci(r.geom_margin), " chart=",
                    sci(r.sigma_chart), " margin=", sci(r.margin), " dual=",
                    sci(r.dual), " gap=", sci(r.gap), " [", fx(r.secs, 1), " s]")
        end
    end
    close(csv)

    # ---- summary -------------------------------------------------------------------
    sm = open(joinpath(outdir, "recheck_k9_summary" * sfx * ".txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("K8a recheck -- the corrected cone LP at K9's variant (b) embeddings\n")
    both("designs with a feasible (b) embedding and F <= $cone_maxf : $(length(rows))\n")

    n_sig = 0; n_pos = 0; n_ge_chart = 0; n_ge_geom = 0; n_strictly_better = 0
    n_nosig = 0; n_nosig_pos = 0; n_cert = 0
    gaps = Float64[]; ratios = Float64[]; margins = Float64[]
    for r in rows
        if r.sigma_in_cone != 0
            n_sig += 1
            n_pos += (r.margin > 1e-8)
            n_ge_chart += (r.margin >= r.sigma_chart - 1e-12)
            n_ge_geom += (r.margin > 0 && r.geom_margin > 0)
            r.margin > r.sigma_chart * (1.0 + 1e-6) && (n_strictly_better += 1)
            r.sigma_chart > 0 && push!(ratios, r.margin / r.sigma_chart)
            push!(margins, r.margin)
        else
            n_nosig += 1
            n_nosig_pos += (r.margin > 1e-8)
            n_cert += r.cert_1e9
        end
        push!(gaps, r.gap)
    end
    both("sigma in P(X) by direct zero_plus measurement       : $n_sig\n")
    both("  LP margin > 0 (was 0 with the old solver)         : $n_pos/$n_sig\n")
    both("  LP margin >= the solver-free sigma chart bound    : $n_ge_chart/$n_sig\n")
    both("  LP margin > 0 where min(q, mu)/2 > 0              : $n_ge_geom/$n_sig\n")
    both("  a NON-UNIFORM flex strictly beats the sigma ray   : $n_strictly_better/$n_sig\n")
    isempty(ratios) || both("  margin / sigma_chart: median " * sci(median_of(ratios)) * ", max " *
                            sci(maximum(ratios)) * "\n")
    isempty(margins) || both("  margin: median " * sci(median_of(margins)) * ", max " *
                             sci(maximum(margins)) * "\n")
    both("sigma OUTSIDE P(X)                                  : $n_nosig\n")
    both("  LP margin > 0 anyway (a genuine non-uniform rescue): $n_nosig_pos/$n_nosig\n")
    both("  every pass Farkas-certified < 1e-9 (verified)     : $n_cert/$n_nosig\n")
    isempty(gaps) || both("bracket width dual - margin: median " * sci(median_of(gaps)) * ", max " *
                          sci(maximum(gaps)) * "\n")
    both("wall " * fx(s(wall), 1) * " s\n")
    close(sm)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
