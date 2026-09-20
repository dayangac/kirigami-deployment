# K1a -- is the baseline really invalid, and is the shape space really usable?
# (ideas/ranking.md Sec. 4, K1a.)
#
# PASS rule, copied from ranking.md:
#   "PASS-algorithm if n_inv(X0) > 0 on >= 150/200 graphs and p_valid > 0 on >= 20
#    graphs. PASS-emptiness if n_inv(X0) > 0 on >= 150/200 and p_valid = 0 on >= 20
#    graphs. FAIL (premise dead, F17 was a sweep artefact) if n_inv(X0) = 0 on
#    > 100/200. Both PASS branches keep R1 alive; they select which paper it is."
# Tolerance: a face counts as inverted if its signed area <= 1e-12 * (bbox area).
#
# Samples: X = X0 + Phi t, t ~ N(0, s^2 I) with s calibrated so the median
# ||Phi t||_inf equals the median edge length of X_ini. In addition to the spec's
# Gaussian cloud, a local trust-region probe re-runs the same count at a ladder of
# radii r * (median edge length), r in {2, 1, 0.5, 0.25, 0.1, 0.03, 0.01}, which
# answers whether the invalidity of X0 is local or global in the shape space.
#
#   julia --project=Kirigami Kirigami/apps/exp_k1a_projection_validity.jl [--n 200] [--samples 10000]
#         [--out DIR] [--cache DIR] [--limit K] [--regenerate]
#
# Outputs go to results/experiments/k1a/ by default; `--cache <dir>` reuses a shape cache (same
# X0, same Phi, hence the same samples).
include(joinpath(@__DIR__, "exp_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# signed area of every face under the flattened positions in `X` (N x 2)
function count_inverted_matrix(faces::Vector{Vector{Int}}, X::AbstractMatrix{Float64}, tol::Float64)
    n = 0
    for f in faces
        s = 0.0
        k = length(f)
        for i in 1:k
            a = f[i]; b = f[mod1(i + 1, k)]
            s += X[a, 1] * X[b, 2] - X[a, 2] * X[b, 1]
        end
        0.5 * s <= tol && (n += 1)
    end
    return n
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n = 200; n_samples = 10000; batch = 250; n_theta = 25; limit = typemax(Int)
    outdir = joinpath(REPO, "results", "experiments", "k1a")
    cache = joinpath(REPO, "results", "experiments", "cache")
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--samples" && i < length(args); n_samples = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--cache" && i < length(args); cache = args[i+1]; i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    radii = [2.0, 1.0, 0.5, 0.25, 0.1, 0.03, 0.01]
    csv = open(joinpath(outdir, "k1a.csv"), "w")
    print(csv, "id,kind,N,F,n_interior,n_split,H,dim_null,n_inv_Xini,n_inv_X0,frac_inv_X0,max_move,",
          "xini_injective,sigma_cal,n_samples,n_valid,p_valid,min_n_inv,x0_injective,n_inj_tested,n_injective,",
          "theta_max_X0,theta_max_best,n_theta_tried")
    for r in radii
        print(csv, ",p_valid_r", fmt_g(r))
    end
    print(csv, ",secs\n")
    wall = Timer()

    acc = (graphs = Ref(0), x0_invalid = Ref(0), x0_valid = Ref(0), pvalid_pos = Ref(0), pvalid_zero = Ref(0),
           xini_invalid = Ref(0), xini_inj_n = Ref(0), x0_inj = Ref(0), inj_pos = Ref(0), tm_pos = Ref(0),
           minF = Ref(1 << 30), maxF = Ref(0), sum_frac = Ref(0.0), max_frac = Ref(0.0),
           any_radius_valid = Ref(0))
    gains = Float64[]

    function run_one(name::String, kind::String, m::K.Mesh, id::Int)
        t = Timer()
        K.build_topology!(m)
        isempty(m.sigma) && return
        c = K.make_cut(m)
        hs = K.holes_partition(c)
        sh = shape_space(m, c, hs, cache, id)
        sh.ok || return
        N = K.n_vertices(m); F = K.n_faces(m); k = sh.k
        lo = m.X[1]; hi = m.X[1]
        for p in m.X
            lo = min.(lo, p); hi = max.(hi, p)
        end
        bbox = max(1e-300, (hi[1] - lo[1]) * (hi[2] - lo[2]))
        tol = 1e-12 * bbox
        Xini = K.points_to_matrix(m.X)
        n_inv_ini = count_inverted_matrix(m.faces, Xini, tol)
        n_inv_X0 = count_inverted_matrix(m.faces, sh.X0, tol)
        max_move = 0.0
        for v in 1:N
            max_move = max(max_move, hypot(sh.X0[v, 1] - Xini[v, 1], sh.X0[v, 2] - Xini[v, 2]))
        end

        med_edge = median_edge_length(m)
        rng = K.MT19937(UInt32(4242) + UInt32(97) * UInt32(id + 1))
        G = K.NormalDist(0.0, 1.0)

        # calibrate s so that median ||Phi t||_inf == median edge length
        s_cal = 0.0
        if k > 0
            pilot = 120
            T = Matrix{Float64}(undef, k, 2 * pilot)
            for i in eachindex(T)   # column-major draw order
                T[i] = K.normal(G, rng)
            end
            D = sh.Phi * T
            inf = [maximum(abs, view(D, :, 2b+1:2b+2)) for b in 0:pilot-1]
            med = sort(inf)[pilot ÷ 2 + 1]   # nth_element at pilot/2 (0-based)
            s_cal = med > 0 ? med_edge / med : 0.0
        end

        n_valid = 0; min_inv = 1 << 30
        # "no inverted face" is weaker than "embedded": the first `inj_budget` valid samples
        # are also tested for injectivity with has_collision at theta = 0.
        inj_budget = 200
        n_inj_tested = 0; n_injective = 0
        valid_samples = Matrix{Float64}[]
        if k > 0 && s_cal > 0
            done = 0
            while done < n_samples
                B = min(batch, n_samples - done)
                T = Matrix{Float64}(undef, k, 2 * B)
                for i in eachindex(T)
                    T[i] = K.normal(G, rng) * s_cal
                end
                D = sh.Phi * T
                for b in 0:B-1
                    X = sh.X0 + D[:, 2b+1:2b+2]
                    ni = count_inverted_matrix(m.faces, X, tol)
                    min_inv = min(min_inv, ni)
                    if ni == 0
                        n_valid += 1
                        if n_inj_tested < inj_budget
                            n_inj_tested += 1
                            inj = !K.has_collision(c, K.deploy(c, K.matrix_to_points(X), 0.0).Y)
                            if inj
                                n_injective += 1
                                length(valid_samples) < n_theta && push!(valid_samples, X)
                            end
                        end
                    end
                end
                done += batch
            end
        else
            min_inv = n_inv_X0
            n_inv_X0 == 0 && (n_valid = n_samples)
        end
        p_valid = n_valid / n_samples

        # trust-region ladder
        pr = zeros(length(radii))
        if k > 0 && s_cal > 0
            R = 200
            for (ri, r) in enumerate(radii)
                T = Matrix{Float64}(undef, k, 2 * R)
                for i in eachindex(T)
                    T[i] = K.normal(G, rng) * s_cal * r
                end
                D = sh.Phi * T
                ok = 0
                for b in 0:R-1
                    X = sh.X0 + D[:, 2b+1:2b+2]
                    count_inverted_matrix(m.faces, X, tol) == 0 && (ok += 1)
                end
                pr[ri] = ok / R
            end
        end

        # theta_max of X0 and of the best valid sample
        tm0 = -1.0; tmb = -1.0
        xini_inj = !K.has_collision(c, K.deploy(c, m.X, 0.0).Y)
        x0_inj = (n_inv_X0 == 0) && !K.has_collision(c, K.deploy(c, K.matrix_to_points(sh.X0), 0.0).Y)
        x0_inj && (tm0 = K.theta_max(c, K.matrix_to_points(sh.X0), 40, 25).theta_max_geometric)
        for X in valid_samples
            tm = K.theta_max(c, K.matrix_to_points(X), 40, 25).theta_max_geometric
            tmb = max(tmb, tm)
        end
        print(csv, id, ",", kind, ",", N, ",", F, ",", K.n_interior_vertices(m), ",",
              K.n_split(c), ",", K.n_interior_holes(hs), ",", k, ",", n_inv_ini, ",",
              n_inv_X0, ",", fmt_g(n_inv_X0 / F), ",", fmt_g(max_move), ",", fmt_b(xini_inj), ",", fmt_g(s_cal), ",",
              n_samples, ",", n_valid, ",", fmt_g(p_valid), ",",
              (min_inv == (1 << 30) ? -1 : min_inv), ",", fmt_b(x0_inj), ",", n_inj_tested,
              ",", n_injective, ",", fmt_g(tm0), ",", fmt_g(tmb), ",", length(valid_samples))
        for v in pr
            print(csv, ",", fmt_g(v))
        end
        print(csv, ",", fmt_g(s(t)), "\n")
        flush(csv)
        acc.graphs[] += 1
        acc.xini_invalid[] += (n_inv_ini > 0)
        acc.x0_invalid[] += (n_inv_X0 > 0)
        acc.x0_valid[] += (n_inv_X0 == 0)
        acc.pvalid_pos[] += (p_valid > 0)
        acc.pvalid_zero[] += (p_valid == 0)
        acc.minF[] = min(acc.minF[], F)
        acc.maxF[] = max(acc.maxF[], F)
        acc.sum_frac[] += n_inv_X0 / F
        acc.max_frac[] = max(acc.max_frac[], n_inv_X0 / F)
        any(v -> v > 0, pr) && (acc.any_radius_valid[] += 1)
        acc.xini_inj_n[] += xini_inj
        acc.x0_inj[] += x0_inj
        acc.inj_pos[] += (n_injective > 0)
        acc.tm_pos[] += (tmb > 0)
        (tm0 > 0 && tmb > 0) && push!(gains, tmb / tm0)
        if acc.graphs[] % 10 == 0
            println("  ", acc.graphs[], " graphs, X0 invalid ", acc.x0_invalid[], ", p_valid>0 ",
                    acc.pvalid_pos[], ", ", fmt_g(s(wall)), " s")
        end
    end

    pop = population("k1a_200"; regenerate = regenerate)
    nrun = 0
    for g in pop
        g.id < n || continue
        nrun >= limit && break
        g.ok || continue
        run_one(g.kind * "_" * string(g.id), g.kind, g.mesh, g.id)
        nrun += 1
    end
    close(csv)

    sm = open(joinpath(outdir, "summary.txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("K1a: $(acc.graphs[]) graphs, F in [$(acc.minF[]),$(acc.maxF[])]\n")
    both("n_inv(X_ini) > 0 (must be 0)      : $(acc.xini_invalid[])\n")
    both("n_inv(X0) > 0                     : $(acc.x0_invalid[])/$(acc.graphs[])\n")
    both("n_inv(X0) == 0                    : $(acc.x0_valid[])/$(acc.graphs[])\n")
    both("mean inverted fraction of X0      : " * fx(acc.graphs[] > 0 ? acc.sum_frac[] / acc.graphs[] : 0.0) *
         ", max " * fx(acc.max_frac[]) * "\n")
    both("p_valid > 0                       : $(acc.pvalid_pos[])/$(acc.graphs[])\n")
    both("p_valid == 0                      : $(acc.pvalid_zero[])/$(acc.graphs[])\n")
    both("any trust-region radius valid     : $(acc.any_radius_valid[])/$(acc.graphs[])\n")
    both("X_ini injective (control)         : $(acc.xini_inj_n[])/$(acc.graphs[])\n")
    both("X0 injective (no face overlap)    : $(acc.x0_inj[])/$(acc.graphs[])\n")
    both("some sample injective             : $(acc.inj_pos[])/$(acc.graphs[])\n")
    both("some sample with theta_max > 0    : $(acc.tm_pos[])/$(acc.graphs[])\n")
    if !isempty(gains)
        sort!(gains)
        both("theta_max(best sample)/theta_max(X0), median = " * fx(gains[length(gains) ÷ 2 + 1]) *
             " over $(length(gains)) graphs\n")
    end
    both("wall " * fx(s(wall), 1) * " s\n")
    close(sm)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
