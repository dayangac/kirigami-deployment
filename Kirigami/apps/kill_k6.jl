# K6 -- 0+ repair in the null space.
#
# Background (K5 / STATE.md F30). On the 200 random graphs of K1a, every projected
# design -- with sigma from Eq. (1) (sigma_mc) and with the defect-minimising sigma of
# K5 (sigma_def) -- has Theta_max = 0, because some split-cut duplicate translates INTO
# its neighbour at theta = 0+. By T1.B this is a SIGN condition, one scalar per split edge:
#     q_e(t) = det( dS_e(t), d_e(t) ) > 0   <=>  the cut e opens at 0+
# (method/zero_plus.jl). q_e is QUADRATIC in the shape-space coefficients t, as is every
# face signed area a_f. This experiment asks, in order:
#   (1) how many split edges are "inward" (q_e <= 0) at t = 0, for both sigma;
#   (1b) how many CORNER incidences are inward (mu <= 0) at t = 0;
#   (2) can a point of the null space be found with all q_e > 0, all a_f > 0 and all
#       mu > 0 -- PRIMARY objective (pass 3) with the corner and proximity terms, SECONDARY
#       (pass 2) split-only, by L-BFGS from t = 0 and 3 Gaussian starts, same seeds;
#   (3) at every such point, the full validity certificate at eps = 0.3, the largest
#       certified eps, the exact Theta_max (T4.2") and the refereeing bisection;
#   (4) the authors' native `prevent` from the same X0, where the CLI is available.
#
#   julia --project=Kirigami Kirigami/apps/kill_k6.jl [--n 200] [--maxf 800] [--out DIR]
#         [--sigma DIR] [--cache DIR] [--eps 0.3] [--iters 1200] [--starts 3] [--cli PATH]
#         [--work DIR] [--native-maxf 550] [--native-max 12] [--native-timeout 420]
#         [--shard i] [--nshards k] [--w-corner 0.05] [--w-prox 1e3] [--no-native]
#         [--limit K] [--regenerate]
#
# sigma_def comes from K5's saved orientations: `--sigma DIR` reads <DIR>/<kind>_<id>.json
# (K5's layout); without it the archived K5 sigma_def of the frozen population
# (data/corpus/k1a_200.json, verified bit-identical to results/kill/k5/sigma) is used.
# Outputs go to results/kill/k6/; the shape cache lives under <out>/cache.
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# The referee: bisection on the true polygon overlap, shrink 1e-12, 4000-point grid.
# Identical to K2a/K2b/K5 so the numbers are comparable.
function referee_theta(c::K.CutStructure, X::Vector{Vec2}, grid::Int = 4000, iters::Int = 50)
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
    for _ in 1:iters
        mid = 0.5 * (lo + hi)
        col(mid) ? (hi = mid) : (lo = mid)
    end
    return lo
end

Base.@kwdef mutable struct CertReport
    valid_at_eps::Bool = false   # POS /\ NOOVERLAP(eps/2) /\ NOROOT(eps) at eps = 0.3
    pos::Bool = false; noovl::Bool = false; noroot::Bool = false
    eps_max::Float64 = 0.0       # largest eps for which the certificate holds
    theta_exact::Float64 = 0.0   # T4.2" Theta_max
    theta_bisect::Float64 = 0.0  # the referee
    zero_range::Bool = false
    i_star::Int = -1
end

# The certificate predicate is MONOTONE in eps, so the largest certified eps is the first
# deflated root (or pi if there is none), probed 1e-12 below it (the root routine accepts a
# root up to eps + 1e-15, so a one-ulp step down would not clear it).
function certify(c::K.CutStructure, X::Vector{Vec2}, eps::Float64)
    r = CertReport()
    B = K.deploy_basis(c, X)
    sd = K.swept_discs(c, B)
    pairs = K.candidate_pairs(c, sd, Float64(pi), true)

    c_eps = K.validity_certificate(c, B, X, pairs, eps)
    r.valid_at_eps = K.valid(c_eps)
    r.pos = c_eps.pos
    r.noovl = c_eps.nooverlap
    r.noroot = c_eps.noroot

    c_pi = K.validity_certificate(c, B, X, pairs, Float64(pi))
    kRootTol = 1e-12
    eps_star = c_pi.first_root > kRootTol ? c_pi.first_root - kRootTol : 0.0
    if c_pi.first_root <= 0
        # no admissible root anywhere in (0, pi]: the certificate can be probed at pi.
        c_star = K.validity_certificate(c, B, X, pairs, Float64(pi))
        K.valid(c_star) && (r.eps_max = Float64(pi))
    end
    if c_pi.pos && eps_star > 0
        c_star = K.validity_certificate(c, B, X, pairs, eps_star)
        K.valid(c_star) && (r.eps_max = eps_star)
    end

    orr = K.exact_theta_max_overlap(c, B, pairs, 1e-9, Float64(pi), 1e-9)
    r.theta_exact = orr.theta_max
    r.zero_range = orr.zero_range
    r.i_star = orr.i_star
    r.theta_bisect = referee_theta(c, X)
    return r
end

# Which contact type binds at a design with Theta_max = 0? Reported for the FAIL case.
function binding_type(c::K.CutStructure, X::Vector{Vec2})
    for v in K.zero_plus_q(c, X)
        v <= 0 && return "split-inward"
    end
    m = c.mesh
    t = K.Mesh(X, m.faces)
    for f in 1:K.n_faces(m)
        K.face_signed_area(t, f) <= 0 && return "inverted"
    end
    B = K.deploy_basis(c, X)
    sd = K.swept_discs(c, B)
    pairs = K.candidate_pairs(c, sd, Float64(pi), true)
    orr = K.exact_theta_max_overlap(c, B, pairs, 1e-9, Float64(pi), 1e-9)
    orr.zero_range && return "vertex-edge"
    return "hinge-wedge"
end

function read_vertices(path::AbstractString, expect_n::Int)
    isfile(path) || return Vec2[]
    j = try
        JSON.parsefile(path)
    catch
        return Vec2[]
    end
    (j isa AbstractDict && haskey(j, "vertices")) || return Vec2[]
    X = Vec2[]
    for v in j["vertices"]
        length(v) < 2 && return Vec2[]
        a = Float64(v[1]); b = Float64(v[2])
        (isfinite(a) && isfinite(b)) || return Vec2[]
        push!(X, Vec2(a, b))
    end
    return length(X) == expect_n ? X : Vec2[]
end

median_of(v::Vector{Float64}) = isempty(v) ? NaN : sort(v)[length(v) ÷ 2 + 1]
function quantile_of(v::Vector{Float64}, p::Float64)
    isempty(v) && return NaN
    w = sort(v)
    i = min(length(w) - 1, trunc(Int, p * (length(w) - 1) + 0.5))   # 0-based rank
    return w[i + 1]
end
max_or(v::Vector{Float64}, d::Float64) = isempty(v) ? d : maximum(v)

Base.@kwdef mutable struct FamilyStats
    name::String = ""
    graphs::Int = 0
    inward::Vector{Float64} = Float64[]        # # split edges with q <= 0 at t = 0
    inward_frac::Vector{Float64} = Float64[]
    n_split::Vector{Float64} = Float64[]
    feasible::Int = 0
    cert_valid::Int = 0                # certificate at eps = 0.3
    certified_positive::Int = 0        # eps_max > 0
    theta_cert::Vector{Float64} = Float64[]    # eps_max, over graphs where it is > 0
    theta_exact::Vector{Float64} = Float64[]
    referee_agree::Int = 0; referee_n::Int = 0
    worst_referee_gap::Float64 = 0.0
    sound::Int = 0; sound_n::Int = 0   # eps_max <= bisection Theta_max
    inward_corner::Vector{Float64} = Float64[]      # # corner incidences with mu <= 0 at t = 0
    inward_corner_frac::Vector{Float64} = Float64[]
    n_corner::Vector{Float64} = Float64[]
    corner_feasible::Int = 0           # primary: all q > 0, all a > 0 AND all mu > 0
    sec_feasible::Int = 0; sec_certified::Int = 0  # secondary: split-only
    sec_theta_cert::Vector{Float64} = Float64[]
    lad_feasible::Int = 0; lad_certified::Int = 0
    s2_n::Int = 0; s2_feasible::Int = 0; s2_certified::Int = 0; s2_theta_pos::Int = 0
    s2_theta_cert::Vector{Float64} = Float64[]
    lad_theta_cert::Vector{Float64} = Float64[]
    binding::Dict{String,Int} = Dict{String,Int}()
    secs::Float64 = 0.0
    # native comparison
    nat_n::Int = 0; nat_tried::Int = 0; nat_failed::Int = 0
    ours_better::Int = 0; theirs_better::Int = 0; ties::Int = 0
    nat_theta::Vector{Float64} = Float64[]; our_theta::Vector{Float64} = Float64[]
    nat_secs::Float64 = 0.0
end

# `std::system("perl -e 'alarm N; exec @ARGV' ...")`: the native CLI under a wall-clock
# guard; returns the exit code (non-zero on timeout or failure).
function run_native(cli::String, timeout::Int, gjson::String, outjson::String)
    cmd = `/usr/bin/perl -e 'alarm shift; exec @ARGV' $timeout $cli prevent $gjson --dump $outjson`
    p = run(pipeline(ignorestatus(cmd); stdout = devnull, stderr = devnull))
    return p.exitcode
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n = 200; maxf = 800; max_iter = 1200; n_random = 3
    eps = 0.3
    outdir = joinpath(REPO, "results", "kill", "k6")
    sigmadir = ""
    cache = ""
    cli = joinpath(REPO, "baseline", "native", "build", "tuttekiri_cli")
    work = "/tmp/kiri_k6"
    native_maxf = 550; native_max = 12; native_timeout = 420
    shard = 0; nshards = 1
    lambdas = [1e-6, 1e-3, 1e-1, 1.0]
    do_native = true
    # K6 pass 3 weights, chosen on a 3-graph pilot (see KILL_REPORT K6).
    w_corner = 0.05; w_prox = 1e3
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--maxf" && i < length(args); maxf = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--sigma" && i < length(args); sigmadir = args[i+1]; i += 2
        elseif a == "--cache" && i < length(args); cache = args[i+1]; i += 2
        elseif a == "--eps" && i < length(args); eps = arg_f(args[i+1]); i += 2
        elseif a == "--iters" && i < length(args); max_iter = arg_i(args[i+1]); i += 2
        elseif a == "--starts" && i < length(args); n_random = arg_i(args[i+1]); i += 2
        elseif a == "--cli" && i < length(args); cli = args[i+1]; i += 2
        elseif a == "--work" && i < length(args); work = args[i+1]; i += 2
        elseif a == "--native-maxf" && i < length(args); native_maxf = arg_i(args[i+1]); i += 2
        elseif a == "--native-max" && i < length(args); native_max = arg_i(args[i+1]); i += 2
        elseif a == "--native-timeout" && i < length(args); native_timeout = arg_i(args[i+1]); i += 2
        elseif a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshards" && i < length(args); nshards = arg_i(args[i+1]); i += 2
        elseif a == "--w-corner" && i < length(args); w_corner = arg_f(args[i+1]); i += 2
        elseif a == "--w-prox" && i < length(args); w_prox = arg_f(args[i+1]); i += 2
        elseif a == "--no-native"; do_native = false; i += 1
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    isempty(cache) && (cache = joinpath(outdir, "cache"))
    mkpath(outdir)
    mkpath(work)
    csv_path = nshards > 1 ? joinpath(outdir, "shard_$(shard).csv") : joinpath(outdir, "k6.csv")
    csv = open(csv_path, "w")
    print(csv, "id,kind,sigma,N,F,med_edge,n_split,dim_null,m,",
          "inward0,inward_frac0,min_q0,n_inv0,n_corner,inward_corner0,min_mu0,",
          "feasible,best_start,iters,min_q,min_area,min_mu,n_bad_q,n_bad_area,n_bad_mu,",
          "t_norm_rel,",
          "cert_eps03,cert_pos,cert_noovl,cert_noroot,eps_max,theta_exact,theta_bisect,",
          "binding,sec_feasible,sec_eps_max,sec_theta_bisect,sec_binding,",
          "lad_feasible,lad_eps_max,lad_lambda,lad_tnorm,lad_theta_bisect,",
          "s2_feasible,s2_eps_max,s2_theta_bisect,",
          "native_theta,native_ok,secs,native_secs\n")
    wall = Timer()

    fam = Dict("sigma_mc" => FamilyStats(name = "sigma_mc"), "sigma_def" => FamilyStats(name = "sigma_def"))

    graphs = 0; missing_sigma = 0
    # Deterministic population index: incremented for every graph that exists AND has a
    # K5 sigma_def, i.e. exactly the 200-graph population K1a/K5 measured, in id order.
    # Shard i of k processes the graphs with gidx % k == i.
    frozen = !regenerate && maxf == 800
    rows = frozen ? load_population(n <= 200 ? "k1a_200" : "e1_900"; first = min(n, limit)) : PopRow[]
    gidx = -1
    nrun = 0
    for id in 0:399
        gidx + 1 < n || break
        nrun >= limit && break
        row = frozen && id + 1 <= length(rows) ? rows[id + 1] : nothing
        g = row === nothing ? K.make_graph(id, 100, maxf, 1400) : K.Graph(row.mesh, row.kind, row.id, row.ok)
        g.ok || continue
        m0 = g.mesh
        K.build_topology!(m0)
        sigma_mc = copy(m0.sigma)

        # sigma_def comes from K5's saved orientations, so K6 studies exactly the designs
        # K5 measured rather than re-running the greedy search.
        sigma_def = Int[]
        if !isempty(sigmadir)
            p = joinpath(sigmadir, g.kind * "_" * string(id) * ".json")
            if isfile(p)
                sm_ = K.load_mesh_json(p)
                (length(sm_.sigma) == K.n_faces(m0) && length(sm_.X) == length(m0.X)) && (sigma_def = sm_.sigma)
            end
        elseif row !== nothing && row.sigma_def !== nothing && length(row.sigma_def) == K.n_faces(m0)
            sigma_def = copy(row.sigma_def)
        end
        if isempty(sigma_def)
            missing_sigma += 1
            continue
        end
        gidx += 1
        gidx >= n && break
        (nshards > 1 && gidx % nshards != shard) && continue
        nrun += 1

        med = median_edge_length(m0)
        counted = false

        for which in 0:1
            fname = which == 1 ? "sigma_def" : "sigma_mc"
            t = Timer()
            m = K._with_sigma(m0, which == 1 ? sigma_def : sigma_mc)
            c = K.make_cut(m)
            K.n_split(c) == 0 && continue
            hs = K.holes_partition(c)
            sh = shape_space(m, c, hs, joinpath(cache, fname), id)
            (!sh.ok || sh.k < 1) && continue
            X0 = K.matrix_to_points(sh.X0)

            F = fam[fname]
            F.graphs += 1
            if !counted
                graphs += 1
                counted = true
            end

            # ---- (1) inward split edges at t = 0 ---------------------------------
            q0 = K.zero_plus_q(c, X0)
            inward0 = count(v -> v <= 0, q0)
            min_q0 = isempty(q0) ? Inf : minimum(q0)
            push!(F.inward, inward0)
            push!(F.inward_frac, isempty(q0) ? 0.0 : inward0 / length(q0))
            push!(F.n_split, K.n_split(c))
            n_inv0 = count_inverted(m, X0)

            # ---- (1b) inward CORNERS at t = 0 (K6 pass 3) -------------------------
            mu0 = K.zero_plus_corner_margin(c, X0)
            inward_c0 = count(v -> v <= 0, mu0)
            min_mu0 = isempty(mu0) ? NaN : minimum(mu0)
            push!(F.inward_corner, inward_c0)
            push!(F.inward_corner_frac, isempty(mu0) ? 0.0 : inward_c0 / length(mu0))
            push!(F.n_corner, length(mu0))

            # ---- (2) the feasibility solve, and (3) the certificate at its minimiser.
            # SECONDARY first: the split-only objective of pass 2, same lambda and seed;
            # its minimiser warm-starts the primary.
            seed0 = UInt32(6000) + UInt32(7) * UInt32(id) + UInt32(which)
            sec_bind = "n/a"
            sc = CertReport()
            sr = let ro = K.ZeroPlusRepairOptions(n_random = n_random, max_iter = max_iter,
                                                  lambda_rel = lambdas[1], seed = seed0)
                K.zero_plus_repair(c, X0, sh.Phi, med, ro)   # w_corner = w_prox = 0
            end
            if sr.feasible
                F.sec_feasible += 1
                sc = certify(c, sr.X, eps)
                if sc.eps_max > 0
                    F.sec_certified += 1
                    push!(F.sec_theta_cert, sc.eps_max)
                end
                sc.theta_exact <= 1e-9 && (sec_bind = binding_type(c, sr.X))
            else
                sec_bind = binding_type(c, sr.X)
            end

            rr = K.ZeroPlusRepairResult()   # the lambda = 1e-6 run
            cr = CertReport()               # its certificate
            bind = "n/a"
            lad_feasible = false
            lad_eps = 0.0; lad_lambda = 0.0; lad_theta = 0.0; lad_tnorm = 0.0
            for li in eachindex(lambdas)
                ro = K.ZeroPlusRepairOptions(n_random = n_random, max_iter = max_iter,
                                             lambda_rel = lambdas[li],
                                             w_corner = w_corner,   # PRIMARY: split + vertex-edge
                                             w_prox = w_prox,       #          with the proximity term
                                             seed = seed0)
                sr.feasible && (ro.t_init = sr.t)   # warm start from the relaxation
                r = K.zero_plus_repair(c, X0, sh.Phi, med, ro)
                c1 = CertReport()
                r.feasible && (c1 = certify(c, r.X, eps))
                if li == 1
                    rr = r
                    cr = c1
                    if r.feasible
                        F.feasible += 1
                        c1.valid_at_eps && (F.cert_valid += 1)
                        if c1.eps_max > 0
                            F.certified_positive += 1
                            push!(F.theta_cert, c1.eps_max)
                        end
                        push!(F.theta_exact, c1.theta_exact)
                        F.referee_n += 1
                        gap = abs(c1.theta_exact - c1.theta_bisect)
                        F.worst_referee_gap = max(F.worst_referee_gap, gap)
                        gap <= 1e-5 && (F.referee_agree += 1)
                        F.sound_n += 1
                        c1.eps_max <= c1.theta_bisect + 1e-9 && (F.sound += 1)
                        c1.theta_exact <= 1e-9 && (bind = binding_type(c, r.X))
                    else
                        bind = binding_type(c, r.X)
                    end
                end
                if r.feasible
                    lad_feasible = true
                    if c1.eps_max > lad_eps
                        lad_eps = c1.eps_max
                        lad_lambda = lambdas[li]
                        lad_theta = c1.theta_bisect
                        lad_tnorm = norm(r.t) / med
                    end
                end
            end
            tnorm = isempty(rr.t) ? 0.0 : norm(rr.t) / med
            rr.feasible && (F.corner_feasible += 1)
            lad_feasible && (F.lad_feasible += 1)
            if lad_eps > 0
                F.lad_certified += 1
                push!(F.lad_theta_cert, lad_eps)
            end
            F.binding[bind] = get(F.binding, bind, 0) + 1

            # ---- stage 2: the existing range optimiser, started from the repaired point.
            s2_eps = 0.0; s2_theta = 0.0
            s2_feasible = 0
            if rr.feasible && sh.k >= 2
                o2 = K.RangeOptOptions(rounds = 10)
                r2 = K.maximize_range(c, rr.X, sh.Phi, o2)
                q2 = K.zero_plus_q(c, r2.X_opt)
                ok = count_inverted(m, r2.X_opt) == 0 && all(v -> v > 0, q2)
                if ok
                    s2_feasible = 1
                    c2 = certify(c, r2.X_opt, eps)
                    s2_eps = c2.eps_max
                    s2_theta = c2.theta_bisect
                end
                F.s2_n += 1
                s2_feasible == 1 && (F.s2_feasible += 1)
                if s2_eps > 0
                    F.s2_certified += 1
                    push!(F.s2_theta_cert, s2_eps)
                end
                s2_theta > 1e-9 && (F.s2_theta_pos += 1)
            end

            # ---- (4) the authors' native `prevent`, from the same X0 -------------
            nat_theta = -1.0
            nat_ok = 0
            nat_secs = 0.0
            if do_native && rr.feasible && K.n_faces(m) <= native_maxf && F.nat_tried < native_max
                gdir = joinpath(work, fname * "_" * string(id))
                rm(gdir; force = true, recursive = true)
                mkpath(gdir)
                K.save_mesh_json(K.Mesh(X0, m.faces, m.sigma, m.periodic, m.half_edges, m.edges, m.edge_index,
                                        m.vertex_half_edges, m.vertex_edges, m.vertex_faces, m.vertex_is_boundary),
                                 joinpath(gdir, "g.json"))   # the SAME X0 our repair started from
                F.nat_tried += 1
                tn = Timer()
                rc = run_native(cli, native_timeout, joinpath(gdir, "g.json"), joinpath(gdir, "out.json"))
                nat_secs = s(tn)
                if rc == 0
                    Xb = read_vertices(joinpath(gdir, "out.json"), K.n_vertices(m))
                    if !isempty(Xb)
                        nat_ok = 1
                        nat_theta = referee_theta(c, Xb)
                    end
                end
                nat_ok == 0 && (F.nat_failed += 1)
                if nat_ok == 1
                    F.nat_n += 1
                    F.nat_secs += nat_secs
                    push!(F.nat_theta, nat_theta)
                    push!(F.our_theta, cr.theta_bisect)
                    if cr.theta_bisect > nat_theta + 1e-6
                        F.ours_better += 1
                    elseif nat_theta > cr.theta_bisect + 1e-6
                        F.theirs_better += 1
                    else
                        F.ties += 1
                    end
                end
                rm(gdir; force = true, recursive = true)
            end

            secs = s(t)
            F.secs += secs
            print(csv, id, ",", g.kind, ",", fname, ",", K.n_vertices(m), ",", K.n_faces(m), ",",
                  fmt_g(med), ",", K.n_split(c), ",", sh.k, ",", 2 * sh.k, ",", inward0, ",",
                  fmt_g(isempty(q0) ? 0.0 : inward0 / length(q0)), ",", fmt_g(min_q0), ",",
                  n_inv0, ",", length(mu0), ",", inward_c0, ",", fmt_g(min_mu0), ",",
                  fmt_b(rr.feasible), ",", rr.best_start, ",",
                  rr.iterations, ",", fmt_g(rr.min_q), ",", fmt_g(rr.min_area), ",", fmt_g(rr.min_margin),
                  ",", rr.n_bad_q, ",", rr.n_bad_area, ",", rr.n_bad_margin, ",",
                  fmt_g(tnorm), ",", fmt_b(cr.valid_at_eps), ",",
                  fmt_b(cr.pos), ",", fmt_b(cr.noovl), ",", fmt_b(cr.noroot), ",",
                  fmt_g(cr.eps_max), ",", fmt_g(cr.theta_exact), ",", fmt_g(cr.theta_bisect), ",",
                  bind, ",", fmt_b(sr.feasible), ",", fmt_g(sc.eps_max), ",",
                  fmt_g(sc.theta_bisect), ",", sec_bind, ",",
                  fmt_b(lad_feasible), ",", fmt_g(lad_eps), ",", fmt_g(lad_lambda), ",",
                  fmt_g(lad_tnorm), ",", fmt_g(lad_theta), ",", s2_feasible, ",",
                  fmt_g(s2_eps), ",", fmt_g(s2_theta), ",", fmt_g(nat_theta), ",", nat_ok, ",",
                  fmt_g(secs), ",", fmt_g(nat_secs), "\n")
            flush(csv)
        end
        if graphs % 10 == 0 && counted
            println("  ", graphs, " graphs, feasible mc/def ", fam["sigma_mc"].feasible, "/",
                    fam["sigma_def"].feasible, ", certified ", fam["sigma_mc"].certified_positive, "/",
                    fam["sigma_def"].certified_positive, ", ", fx(s(wall), 1), " s")
        end
    end
    close(csv)

    sm = open(joinpath(outdir, nshards > 1 ? "summary_shard_$(shard).txt" : "summary.txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("K6 -- 0+ repair in the null space. $graphs graphs, eps = " * fx(eps, 3) * "\n")
    both("graphs skipped for a missing K5 sigma_def : $missing_sigma\n")
    both("t lives in the FULL null space, dim m = 2 * dim_null (Phi is N x k and acts on\n")
    both("each coordinate separately). Feasible = all q_e > 0 AND all a_f > 0 strictly.\n\n")

    best_rate = 0.0; best_lad_rate = 0.0
    best_fam = ""; best_lad_fam = ""
    for key in ("sigma_mc", "sigma_def")
        F = fam[key]
        pct(a, b) = fx(b > 0 ? 100.0 * a / b : 0.0, 2)
        both("=== $(F.name) ($(F.graphs) designs) ===\n")
        both("  (1) inward split edges at t = 0 (q_e <= 0)\n")
        both("      median / q1 / q3 / max      : " * fx(median_of(F.inward), 1) * " / " *
             fx(quantile_of(F.inward, 0.25), 1) * " / " * fx(quantile_of(F.inward, 0.75), 1) *
             " / " * fx(max_or(F.inward, 0.0), 1) * "\n")
        both("      median fraction of split edges : " * fx(median_of(F.inward_frac), 4) * "\n")
        both("      designs with ZERO inward edges : $(count(==(0.0), F.inward)) / $(F.graphs)\n")
        both("      median |E_split|            : " * fx(median_of(F.n_split), 1) * "\n")
        both("  (1b) inward CORNER incidences at t = 0 (mu <= 0)\n")
        both("      median / q1 / q3 / max      : " * fx(median_of(F.inward_corner), 1) * " / " *
             fx(quantile_of(F.inward_corner, 0.25), 1) * " / " *
             fx(quantile_of(F.inward_corner, 0.75), 1) * " / " *
             fx(max_or(F.inward_corner, 0.0), 1) * "\n")
        both("      median fraction of incidences  : " * fx(median_of(F.inward_corner_frac), 4) * "\n")
        both("      designs with ZERO inward corners : $(count(==(0.0), F.inward_corner)) / $(F.graphs)\n")
        both("      median # incidences         : " * fx(median_of(F.n_corner), 1) * "\n")
        both("  (2) 0+-feasible after the L-BFGS repair : $(F.feasible) / $(F.graphs)  (" *
             pct(F.feasible, F.graphs) * " %)\n")
        both("  (3) certificate POS /\\ NOOVERLAP(eps/2) /\\ NOROOT at eps = " * fx(eps, 2) *
             " : $(F.cert_valid) / $(F.graphs)\n")
        both("      certified Theta_max > 0 (largest certified eps) : $(F.certified_positive) / $(F.graphs)  (" *
             pct(F.certified_positive, F.graphs) * " %)\n")
        if !isempty(F.theta_cert)
            both("      certified Theta_max  median / q1 / q3 / max : " *
                 fx(median_of(F.theta_cert)) * " / " * fx(quantile_of(F.theta_cert, 0.25)) * " / " *
                 fx(quantile_of(F.theta_cert, 0.75)) * " / " * fx(maximum(F.theta_cert)) * "\n")
        end
        if !isempty(F.theta_exact)
            both("      exact Theta_max (T4.2\")  median / max : " * fx(median_of(F.theta_exact)) *
                 " / " * fx(maximum(F.theta_exact)) * "\n")
        end
        both("      referee: |Theta_exact - bisection| <= 1e-5 on $(F.referee_agree) / $(F.referee_n), worst " *
             sci(F.worst_referee_gap) * "\n")
        both("      soundness: certified eps <= bisection Theta_max on $(F.sound) / $(F.sound_n)\n")
        both("  (3b) SECONDARY (split-only objective, pass 2), same designs and seeds\n")
        both("      0+-feasible : $(F.sec_feasible) / $(F.graphs), certified Theta_max > 0 on " *
             "$(F.sec_certified) / $(F.graphs)  (" * pct(F.sec_certified, F.graphs) * " %)\n")
        if !isempty(F.sec_theta_cert)
            both("      secondary certified Theta_max  median / max : " *
                 fx(median_of(F.sec_theta_cert)) * " / " * fx(maximum(F.sec_theta_cert)) * "\n")
        end
        both("      lambda ladder {1e-6,1e-3,1e-1,1} (extension): feasible $(F.lad_feasible) / $(F.graphs)" *
             ", certified Theta_max > 0 on $(F.lad_certified) / $(F.graphs)  (" *
             pct(F.lad_certified, F.graphs) * " %)\n")
        if !isempty(F.lad_theta_cert)
            both("      ladder certified Theta_max  median / q1 / q3 / max : " *
                 fx(median_of(F.lad_theta_cert)) * " / " * fx(quantile_of(F.lad_theta_cert, 0.25)) *
                 " / " * fx(quantile_of(F.lad_theta_cert, 0.75)) * " / " * fx(maximum(F.lad_theta_cert)) * "\n")
        end
        both("      stage 2 (repair, then range_opt): run on $(F.s2_n), still 0+-feasible $(F.s2_feasible)" *
             ", certified Theta_max > 0 on $(F.s2_certified) / $(F.graphs)  (" *
             pct(F.s2_certified, F.graphs) * " %), bisection Theta_max > 0 on $(F.s2_theta_pos)\n")
        if !isempty(F.s2_theta_cert)
            both("      stage-2 certified Theta_max  median / max : " *
                 fx(median_of(F.s2_theta_cert)) * " / " * fx(maximum(F.s2_theta_cert)) * "\n")
        end
        both("      binding contact type at the failures:")
        for k in sort!(collect(keys(F.binding)))   # std::map order
            both("  " * k * " " * string(F.binding[k]))
        end
        both("\n")
        if F.nat_n > 0
            both("  (4) vs the authors' native `prevent` from the same X0, refereed by OUR bisection\n")
            both("      designs compared : $(F.nat_n) of $(F.nat_tried) attempted (timed out or failed: $(F.nat_failed))\n")
            both("      ours better / theirs better / tie : $(F.ours_better) / $(F.theirs_better) / $(F.ties)\n")
            both("      median Theta_max ours / native : " * fx(median_of(F.our_theta)) * " / " *
                 fx(median_of(F.nat_theta)) * "\n")
            both("      native wall : " * fx(F.nat_secs, 1) * " s\n")
        else
            both("  (4) native comparison: no design qualified (feasible and |F| <= $native_maxf)\n")
        end
        both("      wall : " * fx(F.secs, 1) * " s\n\n")
        # The headline rate is the better of the primary and the secondary objective.
        rate = F.graphs > 0 ? 100.0 * max(F.certified_positive, F.sec_certified) / F.graphs : 0.0
        if rate > best_rate
            best_rate = rate; best_fam = F.name
        end
        lrate = F.graphs > 0 ? 100.0 * max(F.lad_certified, F.s2_certified) / F.graphs : 0.0
        if lrate > best_lad_rate
            best_lad_rate = lrate; best_lad_fam = F.name
        end
    end
    both("--- PASS rule: certified Theta_max > 0 on >= 20 % for at least one sigma ---\n")
    both("(primary = split + vertex-edge objective with proximity term; secondary =\n")
    both(" split-only. The rate below is the better of the two, per sigma family.)\n")
    both("best family : " * (isempty(best_fam) ? "none" : best_fam) * " at " * fx(best_rate, 2) * " %\n")
    both("VERDICT (lambda = 1e-6, as specified): " * (best_rate >= 20.0 ? "PASS" : "FAIL") * "\n")
    both("with the lambda ladder and stage 2 (extensions), best family " *
         (isempty(best_lad_fam) ? "none" : best_lad_fam) * " at " * fx(best_lad_rate, 2) * " %\n")
    both("VERDICT (extensions): " * (best_lad_rate >= 20.0 ? "PASS" : "FAIL") * "\n")
    both("baselines: Eq. (6) alone 0/200 (K1a); K5 both sigma 0/200\n")
    both("wall " * fx(s(wall), 1) * " s\n")
    close(sm)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
