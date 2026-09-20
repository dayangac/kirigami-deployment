# K8a -- the expansive-cone LP (ideas/ranking_r2.md Sec. 4, `exp_k8a_expansive_cone`; A1 of
# ideas/round2_adversary.md; screened PARTIAL in notes/screen_r2.md against
# Rote-Santos-Streinu 2003, whose expansion cone is the classical antecedent).
#
# QUESTION. F25/F30 measured that on 200 random planar graphs the UNIFORM deployment
# has Theta_max = 0, because split-cut duplicates and vertex copies move inward at
# 0+. Is that a statement about the design space, or about one ray of it? Freeing the
# flex from the embedding turns every 0+ separation condition into a strict LINEAR
# inequality on the face velocities inside ker A(X), so the question
#
#     is there ANY infinitesimal flex that separates every cut and every corner?
#
# is one linear program per branch chart (method/expansive_cone.jl).
#
# PASS rule (ranking_r2.md): a strictly positive certified margin on >= 20 of 100
# graphs. FAIL below that -- the emptiness would then be structural rather than a
# wrong ray. HARD SOUNDNESS KILL: if the LP machinery reports sigma OUTSIDE P(X) on
# any of the four split-bearing authored tilings (snub_square_33434, hexagons_auto,
# truncated_square_488, tiling_3_4_3_12), which deploy with exact Theta_max in
# [1.6, 2.4], the driver is wrong and no number here may be quoted.
#
# ASYMMETRY, STATED UP FRONT. The convex-corner margin is a MAX, so P(X) is a UNION
# of polyhedra. This driver solves the conjunctive relaxation first (a SUBSET of
# P(X), so any feasible witness is a sound existence proof), then the branch the
# uniform ray itself selects, then two repair passes from the best witness so far.
# Feasibility proves non-emptiness; infeasibility of the passes tried proves nothing.
# What DOES support a negative reading is the Frank-Wolfe dual bound: a lambda >= 0 on
# the simplex with ||sum_i lambda_i a_i|| ~ 0 is an approximate Farkas certificate for
# the branch chart it was computed on.
#
# SANITY (task step 2). On the graphs with a positive margin, take a short Euler step
# along the LP flex and check with the exact polygon-overlap predicate that the
# stepped structure is embedded. `contact.jl`'s T4.2" scan is parameterised by the
# UNIFORM angle theta and does not apply to a non-uniform flex; the predicate used
# here, `has_collision` / `polygons_overlap`, is the same exact primitive T4.2" calls
# at every interval midpoint.
#
#   julia --project=Kirigami Kirigami/apps/exp_k8a_expansive_cone.jl [--n 100] [--shard S] [--nshard M]
#         [--out DIR] [--cache DIR] [--x0] [--refs-only] [--no-deployable]
#         [--dual-iters 600] [--limit K] [--regenerate]
#
# Population: the K1a graphs (frozen k1a_200, ids < n; `--regenerate` rebuilds them with
# make_graph). The deployable population here is the frozen (2, 0.2) variant
# (data/corpus/deployable_population_2_0.2.json; `--regenerate` rebuilds it). Outputs go to
# results/experiments/k8a/.
include(joinpath(@__DIR__, "exp_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# The displacement ladder for the Euler step, in median-edge units.
const kEps = [0.5, 0.2, 0.1, 0.03, 0.01, 0.003, 0.001]

Base.@kwdef mutable struct Row
    id::String = ""; kind::String = ""; where::String = ""
    N::Int = 0; F::Int = 0; n_split::Int = 0; n_corner::Int = 0; n_convex::Int = 0; n_rows::Int = 0
    dim_ker_A::Int = 0; dim_flex::Int = 0; components::Int = 0
    flex_resid::Float64 = 0.0
    passes::Int = 0; pass_feasible::Int = -1
    margin_l2::Float64 = 0.0; margin_inf::Float64 = 0.0; dual_bound::Float64 = 0.0; dual_max::Float64 = 0.0
    n_active::Int = 0; n_dual_support::Int = 0
    sigma_resid::Float64 = 0.0
    sigma_is_flex::Int = 0; sigma_bad_q::Int = 0; sigma_bad_mu::Int = 0; sigma_in_cone::Int = 0
    sigma_min_q::Float64 = 0.0; sigma_min_mu::Float64 = 0.0
    euler_eps::Float64 = 0.0   # largest collision-free step, in median-edge units (0 = none)
    euler_tested::Int = 0
    # --- recheck columns (results/experiments/k8a/recheck.md) ---------------------------
    sigma_chart::Float64 = 0.0   # rigorous lower bound on the sigma chart's LP value
    sigma_span::Int = 0          # sigma lies in span(flex basis)
    farkas_resid::Float64 = 0.0  # worst INDEPENDENTLY recomputed ||A^T lambda|| over passes
    farkas_lmin::Float64 = 0.0; farkas_sumerr::Float64 = 0.0
    cert_1e9::Int = 0; cert_1e6::Int = 0  # passes whose recomputed certificate survives the tol
    lp_gap::Float64 = 0.0        # worst bracket width dual - margin over the passes
    secs::Float64 = 0.0
end

const kHeader =
    "id,kind,where,N,F,n_split,n_corner,n_convex,n_rows,dim_ker_A,dim_flex,components," *
    "flex_resid,passes,pass_feasible,margin_l2,margin_inf,dual_bound,dual_max,n_active," *
    "n_dual_support,sigma_resid,sigma_is_flex,sigma_min_q,sigma_min_mu,sigma_bad_q," *
    "sigma_bad_mu,sigma_in_cone,euler_eps,euler_tested,sigma_chart,sigma_span," *
    "farkas_resid,farkas_lmin,farkas_sumerr,cert_1e9,cert_1e6,lp_gap,secs\n"

function write_row(o::IO, r::Row)
    print(o, r.id, ",", r.kind, ",", r.where, ",", r.N, ",", r.F, ",",
          r.n_split, ",", r.n_corner, ",", r.n_convex, ",", r.n_rows, ",",
          r.dim_ker_A, ",", r.dim_flex, ",", r.components, ",", sci(r.flex_resid),
          ",", r.passes, ",", r.pass_feasible, ",", sci(r.margin_l2), ",",
          sci(r.margin_inf), ",", sci(r.dual_bound), ",", sci(r.dual_max), ",",
          r.n_active, ",",
          r.n_dual_support, ",", sci(r.sigma_resid), ",", r.sigma_is_flex, ",",
          sci(r.sigma_min_q), ",", sci(r.sigma_min_mu), ",", r.sigma_bad_q, ",",
          r.sigma_bad_mu, ",", r.sigma_in_cone, ",", fmt_g(r.euler_eps), ",",
          r.euler_tested, ",", sci(r.sigma_chart), ",", r.sigma_span, ",",
          sci(r.farkas_resid), ",", sci(r.farkas_lmin), ",", sci(r.farkas_sumerr),
          ",", r.cert_1e9, ",", r.cert_1e6, ",", sci(r.lp_gap), ",",
          fx(r.secs, 2), "\n")
end

# One (cut structure, embedding): the LP, the sigma ray, and the Euler-step sanity.
function run_one(id::String, kind::String, where::String, m::K.Mesh, c::K.CutStructure,
                 X::Vector{Vec2}, opt::K.ExpansiveConeOptions, do_euler::Bool)
    t = Timer()
    r = Row(id = id, kind = kind, where = where, N = K.n_vertices(m), F = K.n_faces(m))
    cs = K.cone_system(c, X)
    r.n_split = cs.n_split
    r.n_corner = cs.n_corner
    r.n_convex = cs.n_convex

    rep = K.expansive_cone(c, X, opt)
    r.n_rows = rep.n_rows
    r.dim_ker_A = rep.dim_ker_A
    r.dim_flex = rep.dim_flex
    r.components = rep.components
    r.flex_resid = rep.flex_residual
    r.passes = rep.passes
    r.pass_feasible = rep.pass_feasible
    r.margin_l2 = rep.margin_l2
    r.margin_inf = rep.margin_inf
    r.dual_bound = rep.dual_bound
    r.dual_max = rep.dual_bound_max
    r.n_active = rep.n_active
    r.n_dual_support = rep.n_dual_support
    r.sigma_resid = rep.sigma_residual
    r.sigma_is_flex = Int(rep.sigma_is_flex)
    r.sigma_min_q = rep.sigma_min_q
    r.sigma_min_mu = rep.sigma_min_mu
    r.sigma_bad_q = rep.sigma_n_bad_q
    r.sigma_bad_mu = rep.sigma_n_bad_mu
    r.sigma_in_cone = Int(rep.sigma_in_cone)
    r.sigma_chart = rep.sigma_chart_margin
    r.sigma_span = Int(rep.sigma_in_flex_span)
    r.farkas_resid = rep.farkas_resid_max
    r.farkas_lmin = rep.farkas_lambda_min
    r.farkas_sumerr = rep.farkas_sum_err
    r.cert_1e9 = (rep.passes > 0 && rep.passes_certified_1e9 == rep.passes) ? 1 : 0
    r.cert_1e6 = (rep.passes > 0 && rep.passes_certified_1e6 == rep.passes) ? 1 : 0
    r.lp_gap = rep.lp_gap_max

    if do_euler && rep.margin_l2 > opt.lp.tol && !isempty(rep.flex)
        med = median_edge_length(m)
        sp = K.flex_speed(c, X, rep.flex)
        if sp > 0 && med > 0
            r.euler_tested = 1
            for eps in kEps
                Y = K.euler_step(c, X, rep.flex, eps * med / sp)
                if !K.has_collision(c, Y)
                    r.euler_eps = eps
                    break
                end
            end
        end
    end
    r.secs = s(t)
    return r
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n = 100; shard = 0; nshard = 1
    do_x0 = false; refs_only = false; do_deployable = true
    dual_iters = 600
    limit = typemax(Int)
    outdir = joinpath(REPO, "results", "experiments", "k8a")
    cache = joinpath(REPO, "results", "experiments", "cache")
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshard" && i < length(args); nshard = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--cache" && i < length(args); cache = args[i+1]; i += 2
        elseif a == "--x0"; do_x0 = true; i += 1
        elseif a == "--refs-only"; refs_only = true; i += 1
        elseif a == "--no-deployable"; do_deployable = false; i += 1
        elseif a == "--dual-iters" && i < length(args); dual_iters = arg_i(args[i+1]); i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    suffix = nshard > 1 ? "_" * string(shard) : ""
    csv = open(joinpath(outdir, "k8a" * suffix * ".csv"), "w")
    print(csv, kHeader)

    opt = K.ExpansiveConeOptions()
    opt.lp.dual_iters = dual_iters

    # ---- the four split-bearing authored tilings: the hard soundness control -------
    refs_done = 0; refs_sound = 0; refs_lp_pos = 0
    if shard == 0
        want = ["snub_square_33434", "hexagons_auto", "truncated_square_488", "tiling_3_4_3_12"]
        for rc in reference_cases(regenerate = regenerate)
            rc.name in want || continue
            m = rc.mesh
            K.build_topology!(m)
            c = K.make_cut(m)
            hs = K.holes_partition(c)
            X = copy(m.X)
            if !K.deployable(K.hole_residuals(c, X, hs), 1e-9)
                sys = K.assemble_system(c, hs, m.X, K.Fixed)
                sr = K.solve_system(sys, m.X)
                if !sr.projection_ok
                    println("  reference ", rc.name, ": Eq. (6) projection FAILED, skipped")
                    continue
                end
                X = K.matrix_to_points(sr.X0)
            end
            r = run_one(rc.name, "reference", "deployable", m, c, X, opt, true)
            write_row(csv, r)
            flush(csv)
            refs_done += 1
            refs_sound += r.sigma_in_cone
            refs_lp_pos += (r.margin_l2 > opt.lp.tol)
            println("  ref ", rc.name, ": F=", r.F, " split=", r.n_split,
                    " dim_flex=", r.dim_flex, " sigma_in_cone=", r.sigma_in_cone,
                    " (min q ", sci(r.sigma_min_q), ", min mu ", sci(r.sigma_min_mu),
                    ") LP margin ", sci(r.margin_l2), " [", fmt_g(r.secs), " s]")
        end
    end

    # ---- the deployable population: the Euler-step sanity, and a wider soundness
    # control than the four tilings (exp_common::deployable_population, authored
    # tilings at five clip radii plus shape-space samples, all uniformly deployable and
    # embedded at theta = 0). These are the configurations where the LP is EXPECTED to
    # be feasible, so they are where a first-order flex can actually be integrated.
    dep_done = 0; dep_sound = 0; dep_pos = 0; dep_euler = 0; dep_euler_ok = 0
    if shard == 0 && do_deployable
        # the frozen (2, 0.2) variant (deployable_population_2_0.2.json); --regenerate rebuilds
        for d in deployable_population(regenerate = regenerate, samples_per_base = 2, radius_frac = 0.2)
            m = d.mesh
            K.build_topology!(m)
            c = K.make_cut(m)
            r = run_one(d.name, "deployable_" * d.family, "deployable", m, c, d.X, opt, true)
            write_row(csv, r)
            flush(csv)
            dep_done += 1
            dep_sound += r.sigma_in_cone
            dep_pos += (r.margin_l2 > opt.lp.tol)
            dep_euler += r.euler_tested
            dep_euler_ok += (r.euler_eps > 0)
        end
        println("  deployable population: ", dep_done, " configs, sigma in P(X) ",
                dep_sound, ", LP margin > 0 ", dep_pos, ", Euler collision-free ",
                dep_euler_ok, "/", dep_euler)
    end

    # ---- the K1a population -------------------------------------------------------
    graphs = 0; pos = 0; pos_x0 = 0; sigma_flex = 0; sigma_cone = 0
    euler_tested = 0; euler_ok = 0
    dual_small = 0  # dual bound < 1e-6: an approximate Farkas certificate
    cert9 = 0; cert6 = 0; cert9_x0 = 0  # INDEPENDENTLY verified, every pass
    margins = Float64[]; dims = Float64[]
    wall = Timer()
    if !refs_only
        nrun = 0
        for g in population("k1a_200"; regenerate = regenerate)
            g.id < n || continue
            g.id % nshard == shard || continue
            nrun >= limit && break
            g.ok || continue
            m = g.mesh
            K.build_topology!(m)
            isempty(m.sigma) && continue
            c = K.make_cut(m)
            r = run_one(g.kind * "_" * string(g.id), g.kind, "X_ini", m, c, m.X, opt, true)
            write_row(csv, r)
            flush(csv)
            nrun += 1
            graphs += 1
            pos += (r.margin_l2 > opt.lp.tol)
            sigma_flex += r.sigma_is_flex
            sigma_cone += r.sigma_in_cone
            euler_tested += r.euler_tested
            euler_ok += (r.euler_eps > 0)
            dual_small += (r.dual_max < 1e-6)
            cert9 += r.cert_1e9
            cert6 += r.cert_1e6
            push!(margins, r.margin_l2)
            push!(dims, r.dim_flex)

            if do_x0
                hs = K.holes_partition(c)
                sh = shape_space(m, c, hs, cache, g.id)
                if sh.ok
                    X0 = K.matrix_to_points(sh.X0)
                    r0 = run_one(g.kind * "_" * string(g.id), g.kind, "X0", m, c, X0, opt, false)
                    write_row(csv, r0)
                    flush(csv)
                    pos_x0 += (r0.margin_l2 > opt.lp.tol)
                    cert9_x0 += r0.cert_1e9
                end
            end
            println("  ", r.id, " F=", r.F, " rows=", r.n_rows,
                    " dim_flex=", r.dim_flex, " margin=", sci(r.margin_l2),
                    " dual=", sci(r.dual_bound), " eps=", fmt_g(r.euler_eps), " [",
                    fx(r.secs, 1), " s]")
        end
    end
    close(csv)

    sm = open(joinpath(outdir, "summary" * suffix * ".txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("K8a -- expansive-cone LP\n")
    if shard == 0
        both("reference tilings (split-bearing) tested : $refs_done\n")
        both("  sigma in P(X) (HARD soundness)         : $refs_sound/$refs_done\n")
        both("  LP margin > 0 on the same tilings      : $refs_lp_pos/$refs_done\n")
    end
    if shard == 0 && do_deployable
        both("deployable population configs           : $dep_done\n")
        both("  sigma in P(X)                          : $dep_sound/$dep_done\n")
        both("  LP margin > 0                          : $dep_pos/$dep_done\n")
        both("  Euler step collision-free              : $dep_euler_ok/$dep_euler\n")
    end
    if !refs_only
        both("K1a graphs                               : $graphs\n")
        both("LP margin > 0 at X_ini (PASS >= 20/100)  : $pos/$graphs\n")
        do_x0 && both("LP margin > 0 at X0                      : $pos_x0/$graphs\n")
        both("dual_max < 1e-6 (approx. Farkas)         : $dual_small/$graphs\n")
        both("Farkas VERIFIED < 1e-9, every pass       : $cert9/$graphs\n")
        both("Farkas VERIFIED < 1e-6, every pass       : $cert6/$graphs\n")
        do_x0 && both("Farkas VERIFIED < 1e-9 at X0             : $cert9_x0/$graphs\n")
        both("sigma is a flex at X_ini                 : $sigma_flex/$graphs\n")
        both("sigma in P(X_ini)                        : $sigma_cone/$graphs\n")
        both("Euler step tested / collision-free       : $euler_ok/$euler_tested\n")
        if !isempty(margins)
            sv = sort(margins)
            both("margin: median " * sci(sv[length(sv) ÷ 2 + 1]) * ", max " * sci(sv[end]) * "\n")
            dv = sort(dims)
            both("dim flex: median " * fx(dv[length(dv) ÷ 2 + 1], 0) * ", range [" * fx(dv[1], 0) *
                 ", " * fx(dv[end], 0) * "]\n")
        end
    end
    both("wall " * fx(s(wall), 1) * " s\n")
    close(sm)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
