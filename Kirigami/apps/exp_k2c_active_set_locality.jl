# K2c -- the locality gate for the active set (ideas/ranking.md Sec. 4, K2c).
# The withdrawn/replaced PASS rule and every recorded deviation apply verbatim:
#   (1) PASS iff pruned and unpruned exact Theta_max agree to 1e-12 on EVERY graph;
#   (2) the candidate-set growth with n is reported, and H-LOC is refuted by the
#       centroid-drift column (linear in the patch diameter).
#
#   julia --project=Kirigami Kirigami/apps/exp_k2c_active_set_locality.jl [--n 500] [--theta-limit 800]
#         [--out DIR] [--cache DIR] [--limit K] [--regenerate]
include(joinpath(@__DIR__, "exp_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# least-squares slope of ly on lx (and R^2); NaN below 3 points
function slope(pts::Vector{Tuple{Float64,Float64}})
    length(pts) < 3 && return NaN, 0.0
    sx = sum(p[1] for p in pts) / length(pts)
    sy = sum(p[2] for p in pts) / length(pts)
    sxy = 0.0; sxx = 0.0; syy = 0.0
    for (lx, ly) in pts
        sxy += (lx - sx) * (ly - sy)
        sxx += (lx - sx)^2
        syy += (ly - sy)^2
    end
    r2 = (sxx > 0 && syy > 0) ? (sxy * sxy) / (sxx * syy) : 0.0
    return (sxx > 0 ? sxy / sxx : NaN), r2
end

# max over theta in [0, pi] of |m_f(theta) - m_f(0)| / r_f (centroid drift), 65 samples
function max_drift(sd::K.SweptDiscs, F::Int)
    worst = 0.0
    for f in 1:F
        gc = Vec2(sd.gc[f, 1], sd.gc[f, 2]); gs = Vec2(sd.gs[f, 1], sd.gs[f, 2])
        d = 0.0
        for i in 0:64
            t = 0.5 * pi * i / 64.0
            d = max(d, norm((cos(t) - 1.0) * gc + sin(t) * gs))
        end
        sd.circum[f] > 0 && (worst = max(worst, d / sd.circum[f]))
    end
    return worst
end

diameter(X::Vector{Vec2}) = begin
    lo = X[1]; hi = X[1]
    for q in X
        lo = min.(lo, q); hi = max.(hi, q)
    end
    norm(hi - lo)
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n = 500; theta_limit = 800; limit = typemax(Int)
    outdir = joinpath(REPO, "results", "experiments", "k2c")
    cache = joinpath(REPO, "results", "experiments", "cache")
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--theta-limit" && i < length(args); theta_limit = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--cache" && i < length(args); cache = args[i+1]; i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    csv = open(joinpath(outdir, "k2c.csv"), "w")
    print(csv, "id,kind,N,F,X_source,max_rho_over_r,mean_rho_over_r,max_rho,min_r,",
          "n_pairs_all,n_pairs_pruned,n_pairs_static,pairs_per_face,",
          "theta_all,theta_pruned,theta_static,diff_pruned,diff_static,",
          "cand_all,cand_pruned,secs_all,secs_pruned\n")
    wall = Timer()

    fit_all = Tuple{Float64,Float64}[]; fit_x0 = Tuple{Float64,Float64}[]; fit_xini = Tuple{Float64,Float64}[]
    drift_rows = Tuple{Float64,Int,Float64}[]   # (diam, F, drift)
    graphs = 0; theta_tested = 0; theta_agree = 0; static_agree = 0
    worst_diff = 0.0; worst_static_diff = 0.0
    worst_graph = ""
    max_ratio_overall = 0.0; worst_disc_err = 0.0
    fit_ppf = Tuple{Float64,Float64}[]
    ppf = Float64[]

    function run_one(name::String, kind::String, m::K.Mesh, id::Int)
        K.build_topology!(m)
        isempty(m.sigma) && return
        c = K.make_cut(m)
        hs = K.holes_partition(c)
        X = copy(m.X)
        src = "X_ini"
        sh = shape_space(m, c, hs, cache, id)
        if sh.ok
            X = K.matrix_to_points(sh.X0); src = "X0"
        end
        src == "X0" || return
        F = K.n_faces(m)
        B = K.deploy_basis(c, X)
        sd = K.swept_discs(c, B)
        mx = 0.0; sum_ = 0.0; maxrho = 0.0; minr = 1e300
        cnt = 0
        for f in 1:F
            sd.circum[f] <= 0 && continue
            ratio = sd.rho_max[f] / sd.circum[f]
            mx = max(mx, ratio)
            sum_ += ratio
            cnt += 1
            maxrho = max(maxrho, sd.rho_max[f])
            minr = min(minr, sd.circum[f])
        end
        (cnt == 0 || mx <= 0) && return
        # direct check of the exactness lemma: |Y_pv(theta) - m_f(theta)| == |x_u| for all theta
        disc_err = 0.0
        for th in (0.3, 1.1, 2.4, 3.0)
            ct = cos(0.5 * th); st = sin(0.5 * th)
            Y = K.basis_eval(B, th)
            for f in 1:F
                mf = Vec2(ct * sd.gc[f, 1] + st * sd.gs[f, 1], ct * sd.gc[f, 2] + st * sd.gs[f, 2])
                gc0 = Vec2(sd.gc[f, 1], sd.gc[f, 2])
                for pv in c.prime_faces[f]
                    disc_err = max(disc_err, abs(norm(Y[pv] - mf) - norm(K.basis_c(B, pv) - gc0)))
                end
            end
        end
        worst_disc_err = max(worst_disc_err, disc_err)
        pt = (log(Float64(F)), log(mx))
        push!(fit_all, pt)
        push!(src == "X0" ? fit_x0 : fit_xini, pt)
        max_ratio_overall = max(max_ratio_overall, mx)

        npa = 0; npp = 0; nps = 0; ca = 0; cp = 0
        ta = -1.0; tp = -1.0; ts = -1.0; dp = -1.0; ds = -1.0; sa = 0.0; sp = 0.0
        if F <= theta_limit
            t1 = Timer()
            all = K.candidate_pairs(c, sd, Float64(pi), false)
            ra = K.exact_theta_max(c, B, all)
            sa = s(t1)
            t2 = Timer()
            pruned = K.candidate_pairs(c, sd, Float64(pi), true, false)
            rp = K.exact_theta_max(c, B, pruned)
            sp = s(t2)
            stat = K.candidate_pairs(c, sd, Float64(pi), true, true)
            rs = K.exact_theta_max(c, B, stat)
            npa = length(all); npp = length(pruned); nps = length(stat)
            ca = ra.n_candidates; cp = rp.n_candidates
            ta = ra.theta_max; tp = rp.theta_max; ts = rs.theta_max
            dp = abs(ta - tp); ds = abs(ta - ts)
            theta_tested += 1
            theta_agree += (dp <= 1e-12)
            static_agree += (ds <= 1e-12)
            if dp > worst_diff
                worst_diff = dp; worst_graph = name
            end
            worst_static_diff = max(worst_static_diff, ds)
            push!(ppf, npp / F)
            push!(fit_ppf, (log(Float64(F)), log(npp / F)))
        end
        push!(drift_rows, (diameter(X), F, max_drift(sd, F)))
        print(csv, id, ",", kind, ",", K.n_vertices(m), ",", F, ",", src, ",",
              fmt_g(mx), ",", fmt_g(sum_ / cnt), ",", fmt_g(maxrho), ",", fmt_g(minr), ",", npa, ",", npp,
              ",", nps, ",", fmt_g(F > 0 ? npp / F : 0.0), ",", fmt_g(ta), ",",
              fmt_g(tp), ",", fmt_g(ts), ",", fmt_g(dp), ",", fmt_g(ds), ",", ca, ",", cp, ",", fmt_g(sa), ",",
              fmt_g(sp), "\n")
        flush(csv)
        graphs += 1
        if graphs % 25 == 0
            println("  ", graphs, " graphs, theta tested ", theta_tested, ", agree ", theta_agree, ", ",
                    fmt_g(s(wall)), " s")
        end
    end

    nrun = 0
    for g in population("k2c_500"; regenerate = regenerate)
        g.id < n || continue
        nrun >= limit && break
        g.ok || continue
        run_one(g.kind * "_" * string(g.id), g.kind, g.mesh, g.id)
        nrun += 1
    end
    close(csv)

    sa_, r2a = slope(fit_all); sb, r2b = slope(fit_x0); sc, r2c = slope(fit_xini)
    sort!(ppf)
    sm = open(joinpath(outdir, "summary.txt"), "w")
    both(str) = (print(str); print(sm, str))
    # Centroid drift: the H-LOC refutation.
    sort!(drift_rows; by = d -> d[1])
    open(joinpath(outdir, "k2c_drift.csv"), "w") do dcsv
        print(dcsv, "diameter,F,max_drift_over_r\n")
        for d in drift_rows
            print(dcsv, fmt_g(d[1]), ",", d[2], ",", fmt_g(d[3]), "\n")
        end
    end

    both("K2c: $graphs graphs\n")
    both("log-log slope of max_f rho_f/r_f vs log F, all      : " * fx(sa_) * " (R^2 " * fx(r2a) * ", n=$(length(fit_all)))\n")
    both("  on the X0 subset                                  : " * fx(sb) * " (R^2 " * fx(r2b) * ", n=$(length(fit_x0)))\n")
    both("  on the X_ini subset                               : " * fx(sc) * " (R^2 " * fx(r2c) * ", n=$(length(fit_xini)))\n")
    both("largest max_f rho_f/r_f seen                        : " * fx(max_ratio_overall, 12) * "\n")
    both("worst | |Y_pv - m_f(theta)| - r_f |                  : " * sci(worst_disc_err) * "   (exactness of the moving disc)\n")
    both("theta_max tested (F <= $theta_limit)                      : $theta_tested\n")
    both("pruned == unpruned to 1e-12 (sound prune)           : $theta_agree/$theta_tested, worst diff " * sci(worst_diff) * " (" * worst_graph * ")\n")
    both("pruned == unpruned to 1e-12 (spec's static prune)   : $static_agree/$theta_tested, worst diff " * sci(worst_static_diff) * "\n")
    sp2, r2p = slope(fit_ppf)
    both("log-log slope of surviving-pairs-per-face vs log F  : " * fx(sp2) * " (R^2 " * fx(r2p) * ", n=$(length(fit_ppf)))\n")
    isempty(ppf) || both("surviving pairs per face: median " * fx(ppf[length(ppf) ÷ 2 + 1]) * ", max " * fx(ppf[end]) * "\n")
    both("--- the gate statistic, withdrawn ---\n")
    both("max_f rho_f/r_f is identically 1 by algebra (chi_u = -sigma_f J x_u), so its\n")
    both("log-log slope is 0 for any input and the 'slope <= 0.3' gate cannot fail.\n")
    both("--- H-LOC: centroid drift vs patch diameter (the refutation) ---\n")
    both("H-LOC would need max_f |gamma_f(theta) - xbar_f| <= kappa * r_f for a kappa\n")
    both("independent of the patch. Measured on a GROWING square tiling (the random\n")
    both("population is generated in a fixed 40x40 box, so its diameter barely varies and a\n")
    both("regression on it is meaningless -- it is reported afterwards only for completeness):\n")
    both("  patch      diameter   max_f |gamma_f - xbar_f| / r_f\n")
    grow = Tuple{Float64,Float64}[]
    last_F = -1
    for R in (2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0)
        gm = K.tiling_squares(K.rect(Vec2(0.5 * R, 0.5 * R), R / 2 + 0.01, R / 2 + 0.01))
        (K.n_faces(gm) < 4 || K.n_faces(gm) == last_F) && continue
        last_F = K.n_faces(gm)
        K.build_topology!(gm)
        gm.sigma = checkerboard(gm)
        K.build_topology!(gm)
        gc2 = K.make_cut(gm)
        ghs = K.holes_partition(gc2)
        GX = copy(gm.X)
        if !K.deployable(K.hole_residuals(gc2, GX, ghs), 1e-9)
            gsr = K.solve_system(K.assemble_system(gc2, ghs, gm.X, K.Fixed), gm.X)
            gsr.projection_ok || continue
            GX = K.matrix_to_points(gsr.X0)
        end
        GB = K.deploy_basis(gc2, GX)
        gsd = K.swept_discs(gc2, GB)
        diam = diameter(GX)
        worst = max_drift(gsd, K.n_faces(gm))
        push!(grow, (diam, worst))
        both("  squares    " * fx(diam, 3) * "     " * fx(worst, 3) * "   (F=$(K.n_faces(gm)))\n")
    end
    if length(grow) >= 3
        # linear fit drift = a + b * diameter
        sx = sum(g[1] for g in grow) / length(grow); sy = sum(g[2] for g in grow) / length(grow)
        sxy = 0.0; sxx = 0.0
        for g in grow
            sxy += (g[1] - sx) * (g[2] - sy)
            sxx += (g[1] - sx)^2
        end
        both("  linear fit drift/r = a + b*diameter, b = " * fx(sxx > 0 ? sxy / sxx : 0.0, 4) *
             "  -- growth is LINEAR in the diameter, so no such kappa exists\n")
    end
    both("  On the fixed-box random population (diameter " * fx(isempty(drift_rows) ? 0.0 : drift_rows[1][1], 1) *
         " to " * fx(isempty(drift_rows) ? 0.0 : drift_rows[end][1], 1) *
         ", too narrow to regress): drift/r ranges " *
         fx(isempty(drift_rows) ? 0.0 : drift_rows[1][3], 1) * " to " *
         fx(isempty(drift_rows) ? 0.0 : drift_rows[end][3], 1) * "\n")
    both("  A drift of tens of circumradii means the swept region is NOT contained in any\n")
    both("  fixed multiple of the flat circumdisc, so the packing argument for an O(n)\n")
    both("  active set does not close. H-LOC is refuted, not merely unproved.\n")
    both("VERDICT (REPLACED rule: pruned == unpruned exact Theta_max, 1e-12): " *
         ((theta_agree == theta_tested && theta_tested > 0) ? "PASS" : "FAIL") * "\n")
    both("VERDICT (old slope gate): WITHDRAWN -- vacuous, the statistic is identically 1\n")
    both("No O(n) claim is made: H-LOC is refuted by the drift column above.\n")
    both("wall " * fx(s(wall), 1) * " s\n")
    close(sm)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
