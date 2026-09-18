# K2a -- exact theta_max agrees with bisection (ideas/ranking.md Sec. 4, K2a).
# The recorded corrections and deviations apply verbatim (T4.2" interval scan, tau = 0 deflation, exact moving-centroid broad phase,
# split-free cap Theta_max = min(min beta, pi), population = deployable_population(),
# shrinks 1e-6 / 1e-9 / 1e-12 with the last as the reference).
#
#   julia --project=Kirigami Kirigami/apps/kill_k2a.jl [--n 292] [--out DIR] [--cache DIR]
#         [--limit K] [--regenerate]
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# Referee grid: 4000 samples (step 7.9e-4) resolve the 2.9e-3-wide overlap window of
# trunc_square_R20_s2 that a 180-point scan stepped over.
function bisect_theta_max(c::K.CutStructure, X::Vector{Vec2}, shrink::Float64,
                          grid::Int = 4000, iters::Int = 50)
    col(th) = K.has_collision(c, K.deploy(c, X, th).Y, shrink)
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
    for _ in 1:iters
        mid = 0.5 * (lo + hi)
        if col(mid); hi = mid; else; lo = mid; end
    end
    return lo
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n = 292; limit = typemax(Int)
    outdir = joinpath(REPO, "results", "kill", "k2a")
    cache = ""
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--cache" && i < length(args); cache = args[i+1]; i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    csv = open(joinpath(outdir, "k2a.csv"), "w")
    print(csv, "name,kind,N,F,n_split,X_source,theta_roots,theta_t422,n_candidates,i_star,zero_range,",
          "split_free,min_beta_capped,d_beta_capped,cert_valid,cert_pos,cert_nooverlap,cert_noroot,",
          "cert_class1,cert_class2,cert_class3,cert_zero,cert_undeflated,",
          "bis_1e6,bis_1e9,bis_1e12,theta_hybrid,immediate,d_1e6,d_1e9,d_1e12,dT_1e9,dT_1e12,d12_pure,",
          "tight_flat_embedded,overlap_just_after,binding_edge_type,",
          "binding_faces_adjacency,min_beta,secs\n")
    wall = Timer()

    graphs = 0; ok6 = 0; ok9 = 0; ok12 = 0; split_bind = 0; hinge_bind = 0; none_bind = 0
    nonadj_bind = 0; tight = 0; tight_ok = 0; tangential = 0; n_immediate = 0; ok12_pure = 0
    worst6 = 0.0; worst9 = 0.0; worst12 = 0.0
    worst12_graph = ""
    okT9 = 0; okT12 = 0; n_zero_range = 0; n_graze = 0
    worstT9 = 0.0; worstT12 = 0.0
    worstT12_graph = ""
    tot_candidates = 0; tot_intervals = 0
    n_split_free = 0; ok_beta_cap = 0; n_beta_over_pi = 0
    worst_beta_cap = 0.0
    n_cert_valid = 0; n_cert_pos = 0; n_cert_nooverlap = 0; n_cert_noroot = 0
    tot_undeflated = 0; tot_c1 = 0; tot_c2 = 0; tot_c3 = 0; tot_zero = 0
    n_cert_valid_checked = 0; n_cert_sound = 0; n_cert_violation = 0; n_range_ge_eps = 0

    function run_one(name::String, kind::String, m::K.Mesh, X::Vector{Vec2})
        t = Timer()
        K.build_topology!(m)
        c = K.make_cut(m)
        src = "deployable"
        B = K.deploy_basis(c, X)
        # guard: the trig-linear form must reproduce deploy() on this configuration
        let err = 0.0, sc = 1e-300
            for th in (0.4, 1.5, 2.8)
                Y = K.deploy(c, X, th).Y
                Yb = K.basis_eval(B, th)
                for i in eachindex(Y)
                    err = max(err, norm(Y[i] - Yb[i]))
                    sc = max(sc, norm(Y[i]))
                end
            end
            if err > 1e-9 * sc
                println(stderr, "SKIP (basis mismatch) ", name)
                return
            end
        end
        sd = K.swept_discs(c, B)
        pairs = K.candidate_pairs(c, sd, Float64(pi), true)
        ex = K.exact_theta_max(c, B, pairs)
        ov = K.exact_theta_max_overlap(c, B, pairs, 1e-9, Float64(pi), 1e-9)
        immediate = K.penetrates_immediately(c, X)
        theta_hybrid = immediate ? 0.0 : ex.theta_max
        b6 = bisect_theta_max(c, X, 1e-6)
        b9 = bisect_theta_max(c, X, 1e-9)
        b12 = bisect_theta_max(c, X, 1e-12)
        d6 = abs(theta_hybrid - b6); d9 = abs(theta_hybrid - b9); d12 = abs(theta_hybrid - b12)
        d12_pure = abs(ex.theta_max - b12)
        dT9 = abs(ov.theta_max - b9); dT12 = abs(ov.theta_max - b12)
        tot_candidates += length(ov.candidates)
        tot_intervals += ov.n_intervals_tested

        etype = "none"; adj = "none"
        if ex.first.found
            # corner_e is 0-based in the witness (offset within prime_faces[face_e])
            he = c.face_corner_base[ex.first.face_e] + ex.first.corner_e
            e = m.half_edges[he].edge
            etype = c.edge_type[e] == K.Split ? "split" : c.edge_type[e] == K.Hinge ? "hinge" : "border"
            # adjacency of the two faces in M
            adj = "not_adjacent"
            for ed in m.edges
                ed.n_faces == 2 || continue
                f1 = m.half_edges[ed.he[1]].face; f2 = m.half_edges[ed.he[2]].face
                if (f1 == ex.first.face_e && f2 == ex.first.face_v) ||
                   (f2 == ex.first.face_e && f1 == ex.first.face_v)
                    ei = m.edge_index[ed.key]
                    adj = c.edge_type[ei] == K.Split ? "split_adjacent" :
                          c.edge_type[ei] == K.Hinge ? "hinge_adjacent" : "border_adjacent"
                    break
                end
            end
        end
        # Diagnostics for a disagreement: is the flat state embedded at TIGHT tolerance,
        # and does an actual overlap start just past the closed-form contact angle?
        tight_flat = !K.has_collision(c, K.deploy(c, X, 0.0).Y, 1e-12)
        overlap_after = false
        if ex.first.found
            thp = ex.theta_max * (1.0 + 1e-4) + 1e-9
            thp <= pi && (overlap_after = K.has_collision(c, K.deploy(c, X, thp).Y, 1e-12))
        end
        beta = K.hinge_beta(c, X)
        mb = isempty(beta) ? 0.0 : minimum(beta)
        # check.md D3: on a split-free pattern Theta_max = min(min beta, pi).
        split_free = K.n_split(c) == 0
        mb_cap = min(mb, Float64(pi))
        d_beta = split_free ? abs(ov.theta_max - mb_cap) : -1.0
        if split_free
            n_split_free += 1
            worst_beta_cap = max(worst_beta_cap, d_beta)
            d_beta <= 1e-12 && (ok_beta_cap += 1)
            mb > pi && (n_beta_over_pi += 1)
        end
        # The validity certificate POS /\ NOOVERLAP(eps/2) /\ NOROOT at eps = 0.006.
        cert = K.validity_certificate(c, B, X, pairs, 0.006)
        cv = K.valid(cert)
        n_cert_valid += cv
        n_cert_pos += cert.pos
        n_cert_nooverlap += cert.nooverlap
        n_cert_noroot += cert.noroot
        tot_undeflated += cert.n_roots_undeflated
        tot_c1 += cert.n_class1; tot_c2 += cert.n_class2; tot_c3 += cert.n_class3
        tot_zero += cert.n_identically_zero
        # T5.2b': the certificate must IMPLY Theta_max >= eps. Refereed by the bisection.
        if cv
            n_cert_valid_checked += 1
            b12 >= 0.006 - 1e-9 && (n_cert_sound += 1)
            ov.theta_max < 0.006 - 1e-9 && (n_cert_violation += 1)
        end
        b12 >= 0.006 - 1e-9 && (n_range_ge_eps += 1)
        print(csv, name, ",", kind, ",", K.n_vertices(m), ",", K.n_faces(m), ",",
              K.n_split(c), ",", src, ",", fmt_g(ex.theta_max), ",", fmt_g(ov.theta_max), ",",
              length(ov.candidates), ",", ov.i_star, ",", fmt_b(ov.zero_range), ",",
              fmt_b(split_free), ",", fmt_g(mb_cap), ",", fmt_g(d_beta), ",",
              fmt_b(cv), ",", fmt_b(cert.pos), ",", fmt_b(cert.nooverlap), ",",
              fmt_b(cert.noroot), ",", cert.n_class1, ",", cert.n_class2, ",",
              cert.n_class3, ",", cert.n_identically_zero, ",", cert.n_roots_undeflated,
              ",", fmt_g(b6), ",", fmt_g(b9), ",",
              fmt_g(b12), ",", fmt_g(theta_hybrid), ",", fmt_b(immediate), ",", fmt_g(d6), ",", fmt_g(d9), ",", fmt_g(d12),
              ",", fmt_g(dT9), ",", fmt_g(dT12), ",", fmt_g(d12_pure), ",", fmt_b(tight_flat), ",",
              fmt_b(overlap_after), ",", etype, ",", adj, ",", fmt_g(mb),
              ",", fmt_g(s(t)), "\n")
        flush(csv)
        graphs += 1
        ok6 += (d6 <= 1e-5)
        ok9 += (d9 <= 1e-5)
        ok12 += (d12 <= 1e-5)
        worst6 = max(worst6, d6)
        worst9 = max(worst9, d9)
        if d12 > worst12
            worst12 = d12; worst12_graph = name
        end
        n_immediate += immediate
        okT9 += (dT9 <= 1e-5)
        okT12 += (dT12 <= 1e-5)
        if dT12 > worstT12
            worstT12 = dT12; worstT12_graph = name
        end
        worstT9 = max(worstT9, dT9)
        n_zero_range += ov.zero_range
        (ex.first.found && ov.theta_max > ex.theta_max + 1e-7) && (n_graze += 1)
        ok12_pure += (d12_pure <= 1e-5)
        tight += tight_flat
        tight_flat && (tight_ok += (d12 <= 1e-5))
        (tight_flat && d12 > 1e-5 && !overlap_after) && (tangential += 1)
        if etype == "split"; split_bind += 1
        elseif etype == "hinge"; hinge_bind += 1
        elseif etype == "none"; none_bind += 1
        end
        adj == "not_adjacent" && (nonadj_bind += 1)
        if graphs % 20 == 0
            println("  ", graphs, " graphs, T4.2\" agree(1e-12) ", okT12, ", ", fmt_g(s(wall)), " s")
        end
    end

    nrun = 0
    for d in deployable_population(regenerate = regenerate)
        nrun >= limit && break
        run_one(d.name, d.family, d.mesh, d.X)
        nrun += 1
    end
    close(csv)

    pct(a, b) = b > 0 ? 100.0 * a / b : 0.0
    sm = open(joinpath(outdir, "summary.txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("K2a: $graphs graphs with an embedded flat state\n")
    both("--- T4.2\" interval scan (the corrected closed form) ---\n")
    both("|Theta_max(T4.2\") - bisect(1e-9)|  <= 1e-5 : $okT9/$graphs, worst " * sci(worstT9) * "\n")
    both("|Theta_max(T4.2\") - bisect(1e-12)| <= 1e-5 : $okT12/$graphs, worst " * sci(worstT12) * " (" * worstT12_graph * ")\n")
    both("agreement rate (T4.2\", 1e-12 reference)     : " * fx(pct(okT12, graphs), 2) * " %\n")
    both("configurations with Theta_max = 0 (i* = 0)   : $n_zero_range\n")
    both("configurations where min-over-roots < Theta_max (GRAZE) : $n_graze\n")
    both("mean |C(X)| / mean intervals probed          : " * fx(graphs > 0 ? tot_candidates / graphs : 0.0, 1) *
         " / " * fx(graphs > 0 ? tot_intervals / graphs : 0.0, 1) * "\n")
    both("split-free patterns: Theta_max == min(min beta, pi) to 1e-12 : $ok_beta_cap/$n_split_free, worst " *
         sci(worst_beta_cap) * "\n")
    both("  of those, patterns with min beta > pi (the cap bites)      : $n_beta_over_pi\n")
    both("--- validity certificate POS /\\ NOOVERLAP(eps/2) /\\ NOROOT, eps = 0.006 ---\n")
    both("certificate holds / POS / NOOVERLAP / NOROOT : $n_cert_valid / $n_cert_pos / $n_cert_nooverlap / $n_cert_noroot  (of $graphs)\n")
    both("candidate harmonics by class 1/2/3 + identically zero : $tot_c1 / $tot_c2 / $tot_c3 + $tot_zero\n")
    both("spurious roots in (0,eps) the tau=0 deflation removes : $tot_undeflated\n")
    both("T5.2b': certificate holds => bisection Theta_max >= eps : $n_cert_sound/$n_cert_valid_checked   (violations by the closed form: $n_cert_violation)\n")
    both("configurations that actually have Theta_max >= eps      : $n_range_ge_eps/$graphs   (so the certificate is a strict INNER approximation)\n")
    both("--- min-over-roots (the rule as written in ranking.md) ---\n")
    both("|exact - bisect| <= 1e-5, shrink 1e-6 (shipped) : $ok6/$graphs, worst " * sci(worst6) * "\n")
    both("|exact - bisect| <= 1e-5, shrink 1e-9           : $ok9/$graphs, worst " * sci(worst9) * "\n")
    both("|exact - bisect| <= 1e-5, shrink 1e-12          : $ok12/$graphs, worst " * sci(worst12) * " (" * worst12_graph * ")\n")
    both("configurations that penetrate immediately (range 0) : $n_immediate\n")
    both("pure root enumeration agrees (no theta=0 test)  : $ok12_pure/$graphs\n")
    both("flat state embedded at shrink 1e-12            : $tight/$graphs\n")
    both("  of those, |exact - bisect(1e-12)| <= 1e-5     : $tight_ok/$tight\n")
    both("  disagreements that are tangential contact     : $tangential  (no overlap just past the closed-form angle)\n")
    both("binding edge is a split duplicate / hinge / no contact : $split_bind / $hinge_bind / $none_bind\n")
    both("binding faces not adjacent in M                       : $nonadj_bind\n")
    both("VERDICT (corrected rule, T4.2\" vs bisect 1e-12, all configs): " *
         ((okT12 == graphs && graphs > 0) ? "PASS" : "FAIL") * "\n")
    both("VERDICT (rule as written, min-over-roots): " * ((ok12 == graphs && graphs > 0) ? "PASS" : "FAIL") * "\n")
    both("VERDICT (tightly embedded configurations, tangential contacts allowed): " *
         ((tight > 0 && tight_ok + tangential == tight) ? "PASS" : "FAIL") * "\n")
    both("wall " * fx(s(wall), 1) * " s\n")
    close(sm)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
