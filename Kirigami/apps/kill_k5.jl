# K5 -- is the ORIENTATION the reason the shape space is empty? (ideas/ranking.md R5,
# the Critic's "the orientation objective is wrong: deployability defect beats max-cut".)
# Port of code/apps/kill_k5.cpp, including its app-local `defect_search` (the corpus only
# carries that search's ARCHIVED output as `sigma_def`; the search is reproduced here with
# the bit-exact MT19937 so the archive can be regenerated and checked).
#
# Definitions.
#   D(sigma) = sum over hole preimages C of || sum_{hinge e in C} (x_target - x_source) ||^2,
#   evaluated at X_ini: Eq. (2)'s residual, D = 0 iff X_ini is already uniformly deployable.
#   sigma_mc  = assign_orientation_relaxation (Eq. (1) + the documented repair rules).
#   sigma_def = greedy local search over face flips minimising D: single-face flips and
#               adjacent-PAIR flips, accept strictly improving moves, reject any flip that
#               makes c(Gamma) > 1 or detaches a face. Started from sigma_mc and from 3
#               random sigma, best kept, at most 20*|F| flip attempts per start.
#
# Measured for both sigma: |E_split|, D(sigma), the Eq. (6) projection distance
# ||X0 - X_ini||_inf relative to the median edge length, inverted faces, whether X0 is
# overlap-free at theta = 0, and Theta_max by bisection.
#
# PASS rule (orchestrator, FINAL): sigma_def gives a VALID X0 on >= 20% of graphs, where
# valid means the final certificate POS /\ NOOVERLAP(eps/2) /\ NOROOT with eps = 0.3.
# The earlier "overlap-free X0 at theta = 0" rule is computed and reported as well.
#
#   julia --project=Kirigami Kirigami/apps/kill_k5.jl [--n 200] [--maxf 800] [--cap 20]
#         [--out DIR] [--eps 0.3] [--no-range] [--limit K] [--regenerate]
#
# Population: make_graph(id, 100, maxf, 1400) for id in 0..399 until n graphs = the frozen
# k1a_200 / e1_900 files when maxf == 800. Outputs go to results/kill/k5_julia/ (the C++
# wrote results/kill/k5/); sigma_def is saved to <out>/sigma/<kind>_<id>.json for K6.
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

mutable struct Eval
    ok::Bool        # topology usable: c(Gamma) == 1, no detached face
    D::Float64      # the deployability defect at X_ini
    n_split::Int
    H::Int
    c_gamma::Int
end
Eval() = Eval(false, 0.0, 0, 0, 0)

# One evaluation of D(sigma) plus the topological admissibility test.
function eval_sigma(m0::K.Mesh, sigma::Vector{Int})
    e = Eval()
    m = K._with_sigma(m0, sigma)
    c = K.make_cut(m)
    hg = K.build_hinge_graph(c)
    e.c_gamma = hg.components
    hg.components != 1 && return e  # rejects disconnected Gamma AND isolated faces
    hs = K.holes_partition(c)
    r = K.hole_residuals(c, m.X, hs)
    d = 0.0
    for v in r.per_hole
        d += K._sqnorm2(v)
    end
    e.D = d
    e.n_split = K.n_split(c)
    e.H = length(r.per_hole)
    e.ok = true
    return e
end

mutable struct SearchResult
    sigma::Vector{Int}
    D::Float64
    attempts::Int
    accepted::Int
    restarts_used::Int
    ok::Bool
end
SearchResult() = SearchResult(Int[], Inf, 0, 0, 0, false)

function random_sigma!(sig::Vector{Int}, F::Int, rng::K.MT19937)
    for f in 1:F
        sig[f] = K.uniform_int(rng, 0, 1) != 0 ? 1 : -1
    end
    return sig
end

improves(e::Eval, cur::Eval) = e.ok && e.D < cur.D - 1e-12 * max(1.0, cur.D)

"""The greedy flip search of kill_k5.cpp: `cap` flip attempts per start, 4 starts
(sigma_mc + 3 random), MT19937(seed) for the random sigma and the visiting orders."""
function defect_search(m0::K.Mesh, sigma_mc::Vector{Int}, cap::Int, seed::Integer)
    best = SearchResult()
    rng = K.MT19937(seed)
    F = K.n_faces(m0)
    adj = K.dual_graph(m0)
    dual_edges = Tuple{Int,Int}[]
    for f in 1:F, g in adj[f]
        g > f && push!(dual_edges, (f, g))
    end

    for start in 0:3
        sig = start == 0 ? copy(sigma_mc) : random_sigma!(zeros(Int, F), F, rng)
        cur = eval_sigma(m0, sig)
        if !cur.ok
            # a random sigma can disconnect Gamma; repair by re-randomising a few times
            fixed = false
            for _ in 1:5
                fixed && break
                random_sigma!(sig, F, rng)
                cur = eval_sigma(m0, sig)
                fixed = cur.ok
            end
            cur.ok || continue
        end
        attempts = 0; accepted = 0
        improved = true
        order = collect(1:F)
        while improved && attempts < cap
            improved = false
            K.shuffle!(order, rng)
            for f in order
                attempts >= cap && break
                attempts += 1
                sig[f] = -sig[f]
                e = eval_sigma(m0, sig)
                if improves(e, cur)
                    cur = e
                    accepted += 1
                    improved = true
                else
                    sig[f] = -sig[f]
                end
            end
            K.shuffle!(dual_edges, rng)
            for (f, g) in dual_edges
                attempts >= cap && break
                attempts += 1
                sig[f] = -sig[f]
                sig[g] = -sig[g]
                e = eval_sigma(m0, sig)
                if improves(e, cur)
                    cur = e
                    accepted += 1
                    improved = true
                else
                    sig[f] = -sig[f]
                    sig[g] = -sig[g]
                end
            end
        end
        if cur.ok && cur.D < best.D
            best.sigma = copy(sig)
            best.D = cur.D
            best.attempts = attempts
            best.accepted = accepted
            best.restarts_used = start
            best.ok = true
        end
    end
    return best
end

mutable struct Outcome
    solved::Bool
    D::Float64; proj_rel::Float64; theta::Float64
    n_split::Int; n_inv::Int; dim_null::Int; c_gamma::Int
    overlap_free0::Bool   # the withdrawn theta = 0 test (K1a's predicate)
    cert::Bool            # POS /\ NOOVERLAP(eps/2) /\ NOROOT at eps
    cert_pos::Bool; cert_noovl::Bool; cert_noroot::Bool
    first_overlap::Float64    # smallest ladder angle with an interior overlap, -1 = none
end
Outcome() = Outcome(false, 0.0, 0.0, 0.0, 0, 0, 0, 0, false, false, false, false, false, -1.0)

# The referee: bisection on the true polygon overlap, shrink 1e-12, `grid`-point scan.
function referee_theta(c::K.CutStructure, X::Vector{Vec2}, grid::Int = 2000)
    col(th) = K.has_collision(c, K.deploy(c, X, th).Y, 1e-12)
    col(1e-7) && return 0.0
    lo = 0.0; hi = -1.0
    for i in 1:grid
        th = pi * i / grid
        if col(th)
            hi = th
            break
        end
        lo = th
    end
    hi < 0 && return Float64(pi)
    for _ in 1:50
        mid = 0.5 * (lo + hi)
        col(mid) ? (hi = mid) : (lo = mid)
    end
    return lo
end

function measure(m0::K.Mesh, sigma::Vector{Int}, med_edge::Float64, eps::Float64)
    o = Outcome()
    m = K._with_sigma(m0, sigma)
    c = K.make_cut(m)
    o.c_gamma = K.build_hinge_graph(c).components
    o.n_split = K.n_split(c)
    hs = K.holes_partition(c)
    for v in K.hole_residuals(c, m.X, hs).per_hole
        o.D += K._sqnorm2(v)
    end
    sr = K.solve_system(K.assemble_system(c, hs, m.X, K.Fixed), m.X)
    sr.projection_ok || return o, Vec2[], zeros(0, 0)
    o.solved = true
    o.dim_null = sr.dim_null
    X0 = K.matrix_to_points(sr.X0)
    mx = 0.0
    for i in eachindex(X0)
        mx = max(mx, maximum(abs, X0[i] - m.X[i]))
    end
    o.proj_rel = med_edge > 0 ? mx / med_edge : mx
    o.n_inv = count_inverted(m, X0)
    o.overlap_free0 = !K.has_collision(c, K.deploy(c, X0, 0.0).Y, 1e-12)
    cert = K.validity_certificate(c, X0, eps)
    o.cert = K.valid(cert)
    o.cert_pos = cert.pos
    o.cert_noovl = cert.nooverlap
    o.cert_noroot = cert.noroot
    o.theta = referee_theta(c, X0)
    # Where does overlap actually begin? A log-spaced ladder, so "Theta_max = 0" never
    # rests on one probe at 1e-7.
    for th in (1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 1.0)
        if K.has_collision(c, K.deploy(c, X0, th).Y, 1e-12)
            o.first_overlap = th
            break
        end
    end
    return o, X0, sr.Phi
end

function pearson(a::Vector{Float64}, b::Vector{Float64})
    (length(a) < 3 || length(a) != length(b)) && return NaN
    ma = sum(a) / length(a); mb = sum(b) / length(b)
    sab = saa = sbb = 0.0
    for i in eachindex(a)
        sab += (a[i] - ma) * (b[i] - mb)
        saa += (a[i] - ma)^2
        sbb += (b[i] - mb)^2
    end
    return (saa > 0 && sbb > 0) ? sab / sqrt(saa * sbb) : NaN
end

median_or_nan(v::Vector{Float64}) = isempty(v) ? NaN : sort(v)[length(v) ÷ 2 + 1]

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n = 200; maxf = 800; cap_mult = 20
    eps = 0.3
    outdir = joinpath(REPO, "results", "kill", "k5_julia")
    do_range = true
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--maxf" && i < length(args); maxf = arg_i(args[i+1]); i += 2
        elseif a == "--cap" && i < length(args); cap_mult = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--eps" && i < length(args); eps = arg_f(args[i+1]); i += 2
        elseif a == "--no-range"; do_range = false; i += 1
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(joinpath(outdir, "sigma"))
    csv = open(joinpath(outdir, "k5.csv"), "w")
    print(csv, "id,kind,N,F,med_edge,",
          "mc_c_gamma,mc_n_split,mc_D,mc_dim_null,mc_proj_rel,mc_n_inv,mc_overlap_free0,",
          "mc_cert,mc_cert_pos,mc_cert_noovl,mc_cert_noroot,mc_theta,mc_first_overlap,",
          "def_c_gamma,def_n_split,def_D,def_dim_null,def_proj_rel,def_n_inv,def_overlap_free0,",
          "def_cert,def_cert_pos,def_cert_noovl,def_cert_noroot,def_theta,def_first_overlap,",
          "flips_attempted,flips_accepted,best_start,D_ratio,",
          "range_mc_before,range_mc_after,range_def_before,range_def_after,secs\n")
    wall = Timer()

    graphs = 0; mc_fail = 0; def_fix = 0; mc_ok = 0; def_ok = 0; def_break = 0
    mc_cert_n = 0; def_cert_n = 0; def_disconnected = 0
    mc_pos = 0; mc_noovl = 0; mc_noroot = 0; def_pos = 0; def_noovl = 0; def_noroot = 0
    mc_theta_pos = 0; def_theta_pos = 0
    first_ov_def = Float64[]
    Dmc = Float64[]; Ddef = Float64[]; split_mc = Float64[]; split_def = Float64[]
    projmc = Float64[]; projdef = Float64[]
    gain_mc = Float64[]; gain_def = Float64[]

    # the population: ids 0..399 through make_graph(id, 100, maxf, 1400); frozen when maxf == 800
    frozen = !regenerate && maxf == 800
    want = min(n, limit)
    pop = frozen ? population(want <= 200 ? "k1a_200" : "e1_900"; first = want) : K.Graph[]
    graph_of(id) = (frozen && id + 1 <= length(pop)) ? pop[id + 1] : K.make_graph(id, 100, maxf, 1400)

    for id in 0:399
        graphs >= n && break
        graphs >= limit && break
        g = graph_of(id)
        g.ok || continue
        t = Timer()
        m = g.mesh
        K.build_topology!(m)
        sigma_mc = copy(m.sigma)
        med = median_edge_length(m)

        sr = defect_search(m, sigma_mc, cap_mult * K.n_faces(m), UInt32(7000) + UInt32(id))
        if !sr.ok
            def_disconnected += 1
            continue
        end

        omc, X0mc, Phimc = measure(m, sigma_mc, med, eps)
        odef, X0def, Phidef = measure(m, sr.sigma, med, eps)
        (!omc.solved || !odef.solved) && continue
        graphs += 1

        # Save sigma_def next to the graph so Builder-Method can reuse it.
        K.save_mesh_json(K._with_sigma(m, sr.sigma), joinpath(outdir, "sigma", g.kind * "_" * string(id) * ".json"))

        # K2b's range optimisation, started from each sigma's own X0, refereed by the same
        # bisection. Only run where the start is overlap-free and the shape space is nontrivial.
        rmb = rma = rdb = rda = 0.0
        if do_range
            ro = K.RangeOptOptions(rounds = 10)
            function run_range(sig, X0, Phi, start_ok)
                (!start_ok || size(Phi, 2) < 2) && return 0.0, 0.0
                cc = K.make_cut(K._with_sigma(m, sig))
                before = referee_theta(cc, X0)
                r = K.maximize_range(cc, X0, Phi, ro)
                vo = !K.has_collision(cc, K.deploy(cc, r.X_opt, 0.0).Y, 1e-12)
                after = max(before, vo ? referee_theta(cc, r.X_opt) : 0.0)
                return before, after
            end
            rmb, rma = run_range(sigma_mc, X0mc, Phimc, omc.overlap_free0)
            rdb, rda = run_range(sr.sigma, X0def, Phidef, odef.overlap_free0)
            rmb > 1e-9 && push!(gain_mc, (rma - rmb) / rmb)
            rdb > 1e-9 && push!(gain_def, (rda - rdb) / rdb)
        end

        outcome_cols(o) = join([o.c_gamma, o.n_split, cpp_g(o.D), o.dim_null, cpp_g(o.proj_rel), o.n_inv,
                                cpp_b(o.overlap_free0), cpp_b(o.cert), cpp_b(o.cert_pos), cpp_b(o.cert_noovl),
                                cpp_b(o.cert_noroot), cpp_g(o.theta), cpp_g(o.first_overlap)], ",")
        print(csv, id, ",", g.kind, ",", K.n_vertices(m), ",", K.n_faces(m), ",", cpp_g(med), ",",
              outcome_cols(omc), ",", outcome_cols(odef), ",",
              sr.attempts, ",", sr.accepted, ",", sr.restarts_used, ",",
              cpp_g(omc.D > 0 ? odef.D / omc.D : 0.0), ",",
              cpp_g(rmb), ",", cpp_g(rma), ",", cpp_g(rdb), ",", cpp_g(rda), ",", cpp_g(s(t)), "\n")
        flush(csv)

        push!(Dmc, omc.D); push!(Ddef, odef.D)
        push!(split_mc, omc.n_split); push!(split_def, odef.n_split)
        push!(projmc, omc.proj_rel); push!(projdef, odef.proj_rel)
        mc_ok += omc.overlap_free0
        def_ok += odef.overlap_free0
        mc_cert_n += omc.cert
        def_cert_n += odef.cert
        mc_pos += omc.cert_pos; mc_noovl += omc.cert_noovl; mc_noroot += omc.cert_noroot
        def_pos += odef.cert_pos; def_noovl += odef.cert_noovl; def_noroot += odef.cert_noroot
        mc_theta_pos += (omc.theta > 1e-9)
        def_theta_pos += (odef.theta > 1e-9)
        odef.first_overlap > 0 && push!(first_ov_def, odef.first_overlap)
        if !omc.overlap_free0
            mc_fail += 1
            odef.overlap_free0 && (def_fix += 1)
        elseif !odef.overlap_free0
            def_break += 1
        end
        if graphs % 10 == 0
            println("  ", graphs, " graphs, mc_ok ", mc_ok, ", def_ok ", def_ok, ", ", fx(s(wall), 1), " s")
        end
    end
    close(csv)

    sm = open(joinpath(outdir, "summary.txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("K5: $graphs graphs (R5, the orientation objective)\n")
    both("sigma_def = greedy flip search on the deployability defect D(sigma), 4 starts,\n")
    both("            cap $(cap_mult)*|F| attempts, c(Gamma)=1 enforced\n")
    both("graphs where the search could not keep Gamma connected : $def_disconnected\n")
    both("--- the headline: overlap-free X0 at theta = 0 ---\n")
    both("sigma_mc  gives an overlap-free X0 : $mc_ok/$graphs\n")
    both("sigma_def gives an overlap-free X0 : $def_ok/$graphs\n")
    both("graphs where sigma_mc FAILS        : $mc_fail\n")
    both("  ... and sigma_def SUCCEEDS       : $def_fix (" *
         fx(mc_fail > 0 ? 100.0 * def_fix / mc_fail : 0.0, 2) * " %)   [PASS needs >= 20 %]\n")
    both("graphs where sigma_def breaks one sigma_mc solved : $def_break\n")
    both("--- THE PASS PREDICATE: certificate POS /\\ NOOVERLAP(eps/2) /\\ NOROOT, eps = " *
         fx(eps, 3) * " ---\n")
    both("sigma_mc  valid X0 : $mc_cert_n/$graphs   (POS $mc_pos, NOOVERLAP $mc_noovl, NOROOT $mc_noroot)\n")
    both("sigma_def valid X0 : $def_cert_n/$graphs   (POS $def_pos, NOOVERLAP $def_noovl, NOROOT $def_noroot)\n")
    both("certified Theta_max > 0, sigma_mc / sigma_def : $mc_theta_pos / $def_theta_pos  (of $graphs)\n")
    both("sigma_def: smallest ladder angle with an interior overlap, worst over graphs : " *
         (isempty(first_ov_def) ? "none overlapped" : sci(maximum(first_ov_def))) *
         "\n    (so 'Theta_max = 0' does not rest on a single probe)\n")
    both("--- the defect and the split count ---\n")
    both("median D(sigma_mc) / D(sigma_def)   : " * sci(median_or_nan(Dmc)) * " / " * sci(median_or_nan(Ddef)) * "\n")
    both("correlation of D with |E_split|, sigma_mc  : " * fx(pearson(Dmc, split_mc)) * "\n")
    both("correlation of D with |E_split|, sigma_def : " * fx(pearson(Ddef, split_def)) * "\n")
    both("--- the Eq. (6) projection distance ||X0 - X_ini||_inf / median edge ---\n")
    both("median, sigma_mc  : " * fx(median_or_nan(projmc)) * "\n")
    both("median, sigma_def : " * fx(median_or_nan(projdef)) * "\n")
    rate = mc_fail > 0 ? 100.0 * def_fix / mc_fail : 0.0
    if !isempty(gain_mc) || !isempty(gain_def)
        both("--- range optimisation started from each sigma's X0 ---\n")
        both("designs with a usable start, sigma_mc / sigma_def : $(length(gain_mc)) / $(length(gain_def))\n")
        both("median relative range gain, sigma_mc / sigma_def  : " * fx(median_or_nan(gain_mc)) * " / " *
             fx(median_or_nan(gain_def)) * "\n")
    end
    cert_rate = graphs > 0 ? 100.0 * def_cert_n / graphs : 0.0
    both("--- verdicts ---\n")
    both("FINAL RULE  : sigma_def gives a VALID X0 (certificate POS /\\ NOOVERLAP(eps/2)\n")
    both("              /\\ NOROOT, eps = " * fx(eps, 3) * ") on " * fx(cert_rate, 2) * " % of graphs; needs >= 20 %\n")
    both("VERDICT (FINAL RULE): " * (cert_rate >= 20.0 ? "PASS" : "FAIL") * "\n")
    both("EARLIER RULE: sigma_def gives an overlap-free X0 (the WITHDRAWN theta = 0 test)\n")
    both("              on " * fx(rate, 2) * " % of the graphs where sigma_mc fails\n")
    both("VERDICT (earlier rule, withdrawn predicate): " * (rate >= 20.0 ? "PASS" : "FAIL") * "\n")
    both("Range optimiser vs the authors' native prevent from sigma_def's X0: " *
         (cert_rate >= 20.0 ? "run below" : "NOT RUN -- gated on the final rule passing,\n" *
          "  and there is nothing to optimise: every start has Theta_max = 0.") * "\n")
    both("wall " * fx(s(wall), 1) * " s\n")
    close(sm)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
