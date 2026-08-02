# B4 -- is the FIXED BOUNDARY what empties the shape space, or is it the graph?
# Port of code/apps/kill_b4.cpp.
#
# Background. K6 (results/kill/KILL_REPORT.md) repaired the 0+ jam inside the Tutte
# auxetic null space on 400 designs (200 random graphs x 2 sigma) and certified
# Theta_max > 0 on 0 of 400. K7 ran the SAME Voronoi generator on tori and found
# 11/12 zero-plus feasible and 4/12 fully certified (results/kill/k7/k7_c3_all.csv).
# The two populations differ in exactly one thing: the patch carries the fixed
# boundary rows B of 2026 Eq. (4), the torus does not.
#
# Idea B4 (ideas/round2_theorist_b.md) says the emptiness is an artefact of B. The
# budget identity (B.4) of that file makes the total first-order hole-opening rate a
# functional of the BORDER alone, so freezing the border freezes the budget.
#
# The experiment. Rerun K6's population verbatim -- same graphs (make_graph, ids 0..399
# filtered to the 200 that K5 saved a sigma_def for), same sigma_mc / sigma_def, same
# X_ini, same repair objective, weights, seeds and iteration cap -- under three boundary
# variants:
#
#   fixed      : [L; B] X = [0; T]           the K6 system, run as a CONTROL. Must
#                                            reproduce K6's 0/400.
#   free       : [L; e_0^T] X = [0; x_0]     the boundary rows dropped entirely. Only
#                                            gauge left is the global translation, which
#                                            L cannot see (every L row sums to 0, T7(iv)),
#                                            so it is killed by pinning ONE vertex --
#                                            exactly what BoundaryMode Periodic does
#                                            (core/tutte_auxetic.jl) for the same reason.
#                                            The remaining gauge directions (global
#                                            rotation and scaling: both columns of X0 lie
#                                            in null(L)) are left in, because q_e is
#                                            invariant under rotation and positively
#                                            homogeneous under scaling, so they change no
#                                            label -- they only add 2 harmless dimensions.
#   split_free : B restricted to the boundary vertices that touch NO split edge.
#                F23 found the authors' extra boundary-component row to be an
#                over-constraint; this variant asks the same question one step further --
#                free exactly the border where the cuts (and hence the budget) live, and
#                keep the rest of Eq. (4).
#
# Everything downstream is K6's pipeline unchanged: Eq. (6) projection of the same
# X_ini, zero_plus_repair (secondary split-only, then primary split + vertex-edge +
# proximity warm-started from it), the exact Theta_max of T4.2" as the PRIMARY label
# (the certificate is under repair, STATE.md F32) and the certificate as SECONDARY.
#
# PASS rule (B4): exact Theta_max > 0 on >= 10 % of the 400 designs under `free`,
# against 0/400 under `fixed`. Still 0/400 => B4 is dead.
#
#   julia --project=Kirigami Kirigami/apps/kill_b4.jl [--n 200] [--maxf 800] [--out DIR]
#         [--sigma DIR] [--cache DIR] [--k6cache DIR] [--eps 0.3] [--iters 1200]
#         [--starts 3] [--shard S] [--nshards M] [--w-corner 0.05] [--w-prox 1e3]
#         [--lambda 1e-6] [--no-fixed] [--no-free] [--no-splitfree] [--limit K]
#         [--regenerate]
#
# Population as kill_k9.jl (frozen k1a_200 + archived sigma_def). The shape caches use the
# C++ layout: `--k6cache` for the fixed control, `--cache` for the two variants. Outputs go
# to results/kill/b4_julia/ (the C++ wrote results/kill/b4/).
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# The referee: bisection on the true polygon overlap. Identical to K2a/K2b/K5/K6 so the
# numbers are comparable across the bundle.
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
        if col(mid); hi = mid; else; lo = mid; end
    end
    return lo
end

Base.@kwdef mutable struct CertReport
    valid_at_eps::Bool = false
    pos::Bool = false; noovl::Bool = false; noroot::Bool = false
    eps_max::Float64 = 0.0
    theta_exact::Float64 = 0.0
    theta_bisect::Float64 = 0.0
    zero_range::Bool = false
end

# Byte-for-byte K6's certify(), so `fixed` reproduces K6 and the three variants are
# scored by one rule.
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
    eps_star = (c_pi.first_root > kRootTol) ? (c_pi.first_root - kRootTol) : 0.0
    if c_pi.first_root <= 0
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
    r.theta_bisect = referee_theta(c, X)
    return r
end

# K6's binding_type, verbatim.
function binding_type(c::K.CutStructure, X::Vector{Vec2})
    any(v -> v <= 0, K.zero_plus_q(c, X)) && return "split-inward"
    m = c.mesh
    for f in 1:K.n_faces(m)
        K.polygon_area(X[m.faces[f]]) <= 0 && return "inverted"
    end
    B = K.deploy_basis(c, X)
    sd = K.swept_discs(c, B)
    pairs = K.candidate_pairs(c, sd, Float64(pi), true)
    orr = K.exact_theta_max_overlap(c, B, pairs, 1e-9, Float64(pi), 1e-9)
    return orr.zero_range ? "vertex-edge" : "hinge-wedge"
end

# ---- the three boundary variants -------------------------------------------------
# Which boundary vertices keep their Eq. (4) row. `fixed`: all of them. `free`: none
# (one pin instead). `split_free`: those with no incident split edge.
@enum Variant VFixed VFree VSplitFree
variant_name(v::Variant) = v == VFixed ? "fixed" : (v == VFree ? "free" : "split_free")

function pinned_vertices(m::K.Mesh, c::K.CutStructure, v::Variant)
    N = K.n_vertices(m)
    v == VFree && return [1]
    touches_split = falses(N)
    if v == VSplitFree
        for e in c.split_edges
            touches_split[m.edges[e].key.a] = true
            touches_split[m.edges[e].key.b] = true
        end
    end
    out = [i for i in 1:N if m.vertex_is_boundary[i] && !touches_split[i]]
    isempty(out) && push!(out, 1)  # still need the translation gauge
    return out
end

# [L; rows pinning `pin`] X = [0; X_ini(pin)], then the Eq. (6) projection. Core is not
# touched: assemble_system(..., None) gives L, the rows are appended here.
# Returns (Shape, n_pin).
function shape_variant(m::K.Mesh, c::K.CutStructure, hs::K.HoleSet, v::Variant,
                       cache_dir::String, id::Int)
    pin = pinned_vertices(m, c, v)
    if id >= 0 && !isempty(cache_dir)
        s = load_shape(cache_dir, id)
        s !== nothing && s.N == K.n_vertices(m) && return s, length(pin)
    end
    sys = K.assemble_system(c, hs, m.X, K.None)
    nh = size(sys.A, 1)
    N = sys.N
    A = zeros(nh + length(pin), N)
    rhs = zeros(nh + length(pin), 2)
    A[1:nh, :] = sys.A
    rhs[1:nh, :] = sys.rhs
    for (i, p) in enumerate(pin)
        A[nh + i, p] = 1.0
        rhs[nh + i, 1] = m.X[p][1]
        rhs[nh + i, 2] = m.X[p][2]
    end
    sys.A = A
    sys.rhs = rhs
    sys.n_boundary_rows = length(pin)
    sr = K.solve_system(sys, m.X)
    s = Shape(N, sr.dim_null, sr.X0, sr.Phi, sr.rank_L, sr.H, K.n_interior_vertices(m), sr.projection_ok)
    id >= 0 && !isempty(cache_dir) && s.ok && save_shape(cache_dir, id, s)
    return s, length(pin)
end

function median_of(v::Vector{Float64})
    isempty(v) && return NaN
    sv = sort(v)
    return sv[length(sv) ÷ 2 + 1]
end
function quantile_of(v::Vector{Float64}, p::Float64)
    isempty(v) && return NaN
    sv = sort(v)
    i = min(length(sv) - 1, trunc(Int, p * (length(sv) - 1) + 0.5))
    return sv[i + 1]
end

Base.@kwdef mutable struct Cell  # one (variant, sigma) bucket
    designs::Int = 0
    dim_null::Vector{Float64} = Float64[]; n_split::Vector{Float64} = Float64[]
    min_q0::Vector{Float64} = Float64[]; inward_frac0::Vector{Float64} = Float64[]
    feasible::Int = 0; sec_feasible::Int = 0
    theta_pos::Int = 0          # exact Theta_max > 0  -- the PRIMARY label
    certified::Int = 0          # eps_max > 0          -- the SECONDARY label
    cert_eps03::Int = 0
    theta_pos_vals::Vector{Float64} = Float64[]; theta_cert_vals::Vector{Float64} = Float64[]
    min_q_feas::Vector{Float64} = Float64[]
    referee_agree::Int = 0; referee_n::Int = 0
    worst_referee_gap::Float64 = 0.0
    binding::Dict{String,Int} = Dict{String,Int}()
    secs::Float64 = 0.0
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
    n = 200; maxf = 800; max_iter = 1200; n_random = 3
    eps = 0.3
    outdir = joinpath(REPO, "results", "kill", "b4_julia")
    sigmadir = ""
    cache = joinpath(REPO, "results", "kill", "b4_julia", "cache")
    k6cache = joinpath(REPO, "results", "kill", "k6", "cache")
    shard = 0; nshards = 1
    w_corner = 0.05; w_prox = 1e3; lambda = 1e-6
    do_fixed = true; do_free = true; do_splitfree = true
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--maxf" && i < length(args); maxf = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--sigma" && i < length(args); sigmadir = args[i+1]; i += 2
        elseif a == "--cache" && i < length(args); cache = args[i+1]; i += 2
        elseif a == "--k6cache" && i < length(args); k6cache = args[i+1]; i += 2
        elseif a == "--eps" && i < length(args); eps = arg_f(args[i+1]); i += 2
        elseif a == "--iters" && i < length(args); max_iter = arg_i(args[i+1]); i += 2
        elseif a == "--starts" && i < length(args); n_random = arg_i(args[i+1]); i += 2
        elseif a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshards" && i < length(args); nshards = arg_i(args[i+1]); i += 2
        elseif a == "--w-corner" && i < length(args); w_corner = arg_f(args[i+1]); i += 2
        elseif a == "--w-prox" && i < length(args); w_prox = arg_f(args[i+1]); i += 2
        elseif a == "--lambda" && i < length(args); lambda = arg_f(args[i+1]); i += 2
        elseif a == "--no-fixed"; do_fixed = false; i += 1
        elseif a == "--no-free"; do_free = false; i += 1
        elseif a == "--no-splitfree"; do_splitfree = false; i += 1
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    csv_path = nshards > 1 ? joinpath(outdir, "shard_" * string(shard) * ".csv") :
                             joinpath(outdir, "b4.csv")
    csv = open(csv_path, "w")
    print(csv, "id,kind,sigma,variant,N,F,med_edge,n_bdry,n_pin,n_split,dim_null,m,",
          "inward0,inward_frac0,min_q0,n_inv0,n_corner,inward_corner0,min_mu0,",
          "feasible,best_start,iters,min_q,min_area,min_mu,n_bad_q,n_bad_area,n_bad_mu,",
          "t_norm_rel,sec_feasible,sec_min_q,",
          "cert_eps03,cert_pos,cert_noovl,cert_noroot,eps_max,theta_exact,theta_bisect,",
          "binding,secs\n")
    wall = Timer()

    variants = Variant[]
    do_fixed && push!(variants, VFixed)
    do_free && push!(variants, VFree)
    do_splitfree && push!(variants, VSplitFree)

    cell = Dict{String,Cell}()   # key = "<variant>|<sigma>"

    graphs = 0; missing_sigma = 0
    maxf == 800 || regenerate || @warn "frozen k1a_200 was built with --maxf 800; use --regenerate for maxf = $maxf"
    gidx = -1
    nrun = 0
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
        counted = false

        for which in 0:1
            fname = which == 1 ? "sigma_def" : "sigma_mc"
            m = K._with_sigma(m0, which == 1 ? sigma_def : sigma_mc)
            K.build_topology!(m)
            c = K.make_cut(m)
            K.n_split(c) == 0 && continue
            hs = K.holes_partition(c)
            n_bdry = count(m.vertex_is_boundary)

            for var in variants
                t = Timer()
                n_pin = 0
                # The `fixed` control reads K6's own shape cache, so it is literally the same
                # X0 and Phi K6 measured, not a re-solve.
                sh = if var == VFixed
                    n_pin = n_bdry
                    shape_space(m, c, hs, joinpath(k6cache, fname), id)
                else
                    s_, n_pin = shape_variant(m, c, hs, var,
                                              joinpath(cache, variant_name(var) * "_" * fname), id)
                    s_
                end
                (!sh.ok || sh.k < 1) && continue
                X0 = K.matrix_to_points(sh.X0)
                if !counted
                    graphs += 1
                    counted = true
                end

                C = get!(Cell, cell, variant_name(var) * "|" * fname)
                C.designs += 1
                push!(C.dim_null, sh.k)
                push!(C.n_split, K.n_split(c))

                q0 = K.zero_plus_q(c, X0)
                inward0 = 0
                min_q0 = Inf
                for v in q0
                    v <= 0 && (inward0 += 1)
                    min_q0 = min(min_q0, v)
                end
                push!(C.min_q0, min_q0)
                push!(C.inward_frac0, isempty(q0) ? 0.0 : inward0 / length(q0))
                n_inv0 = count_inverted(m, X0)

                mu0 = K.zero_plus_corner_margin(c, X0)
                inward_c0 = 0
                min_mu0 = Inf
                for v in mu0
                    v <= 0 && (inward_c0 += 1)
                    min_mu0 = min(min_mu0, v)
                end
                isempty(mu0) && (min_mu0 = NaN)

                # K6's seeds exactly: the repair sees the same random starts it saw there.
                seed0 = (UInt32(6000) + UInt32(7) * UInt32(id) + UInt32(which)) % UInt32

                # SECONDARY: split-only, K6 pass 2
                ro2 = K.ZeroPlusRepairOptions()
                ro2.n_random = n_random
                ro2.max_iter = max_iter
                ro2.lambda_rel = lambda
                ro2.seed = seed0
                sr = K.zero_plus_repair(c, X0, sh.Phi, med, ro2)
                sr.feasible && (C.sec_feasible += 1)

                ro = K.ZeroPlusRepairOptions()  # PRIMARY: split + vertex-edge + proximity
                ro.n_random = n_random
                ro.max_iter = max_iter
                ro.lambda_rel = lambda
                ro.w_corner = w_corner
                ro.w_prox = w_prox
                sr.feasible && (ro.t_init = sr.t)
                ro.seed = seed0
                rr = K.zero_plus_repair(c, X0, sh.Phi, med, ro)

                cr = CertReport()
                bind = "n/a"
                if rr.feasible
                    C.feasible += 1
                    cr = certify(c, rr.X, eps)
                    cr.valid_at_eps && (C.cert_eps03 += 1)
                    if cr.eps_max > 0
                        C.certified += 1
                        push!(C.theta_cert_vals, cr.eps_max)
                    end
                    if cr.theta_exact > 1e-9
                        C.theta_pos += 1
                        push!(C.theta_pos_vals, cr.theta_exact)
                    else
                        bind = binding_type(c, rr.X)
                    end
                    push!(C.min_q_feas, rr.min_q)
                    C.referee_n += 1
                    gap = abs(cr.theta_exact - cr.theta_bisect)
                    C.worst_referee_gap = max(C.worst_referee_gap, gap)
                    gap <= 1e-5 && (C.referee_agree += 1)
                else
                    # Some q_e <= 0 or some face inverted at the minimiser. By T1.B a non-positive
                    # q_e means the two copies of that cut translate INTO each other at 0+, so
                    # Theta_max = 0 without a scan; K6 scores these the same way.
                    bind = binding_type(c, rr.X)
                end
                C.binding[bind] = get(C.binding, bind, 0) + 1

                tnorm = isempty(rr.t) ? 0.0 : norm(rr.t) / med
                secs = s(t)
                C.secs += secs
                g_ = cpp_g
                print(csv, id, ",", row.kind, ",", fname, ",", variant_name(var), ",",
                      K.n_vertices(m), ",", K.n_faces(m), ",", g_(med), ",", n_bdry, ",",
                      n_pin, ",", K.n_split(c), ",", sh.k, ",", 2 * sh.k, ",",
                      inward0, ",", g_(isempty(q0) ? 0.0 : inward0 / length(q0)), ",",
                      g_(min_q0), ",", n_inv0, ",", length(mu0), ",", inward_c0, ",",
                      g_(min_mu0), ",", Int(rr.feasible), ",", rr.best_start, ",",
                      rr.iterations, ",", g_(rr.min_q), ",", g_(rr.min_area), ",", g_(rr.min_margin),
                      ",", rr.n_bad_q, ",", rr.n_bad_area, ",", rr.n_bad_margin, ",",
                      g_(tnorm), ",", Int(sr.feasible), ",", g_(sr.min_q), ",",
                      Int(cr.valid_at_eps), ",", Int(cr.pos), ",", Int(cr.noovl), ",",
                      Int(cr.noroot), ",", g_(cr.eps_max), ",", g_(cr.theta_exact), ",",
                      g_(cr.theta_bisect), ",", bind, ",", g_(secs), "\n")
                flush(csv)
            end
        end
        counted && println("  id ", id, " (", graphs, " graphs) ", fx(s(wall), 1), " s")
    end
    close(csv)

    sm = open(joinpath(outdir, nshards > 1 ? "summary_shard_" * string(shard) * ".txt" : "summary.txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("B4 -- does the FIXED BOUNDARY empty the shape space? $graphs graphs, eps = " * fx(eps, 3) * "\n")
    both("graphs skipped for a missing K5 sigma_def : $missing_sigma\n")
    both("PRIMARY label = exact Theta_max (T4.2\") > 0. SECONDARY = certificate eps_max > 0.\n\n")
    for var in variants, sname in ("sigma_mc", "sigma_def")
        key = variant_name(var) * "|" * sname
        haskey(cell, key) || continue
        C = cell[key]
        both("=== " * key * " (" * string(C.designs) * " designs) ===\n")
        both("    median dim_null / |E_split| : " * fx(median_of(C.dim_null), 1) * " / " *
             fx(median_of(C.n_split), 1) * "\n")
        both("    median min q at t = 0       : " * fx(median_of(C.min_q0), 3) *
             ",  median inward fraction " * fx(median_of(C.inward_frac0), 4) * "\n")
        both("    0+-feasible primary / secondary : " * string(C.feasible) * " / " *
             string(C.sec_feasible) * "  of " * string(C.designs) * "  (" *
             fx(C.designs != 0 ? 100.0 * C.feasible / C.designs : 0.0, 2) * " %)\n")
        both("    exact Theta_max > 0 (PRIMARY)   : " * string(C.theta_pos) * " / " *
             string(C.designs) * "  (" *
             fx(C.designs != 0 ? 100.0 * C.theta_pos / C.designs : 0.0, 2) * " %)\n")
        isempty(C.theta_pos_vals) ||
            both("      Theta_max  min / q1 / median / q3 / max : " *
                 fx(minimum(C.theta_pos_vals)) * " / " *
                 fx(quantile_of(C.theta_pos_vals, 0.25)) * " / " * fx(median_of(C.theta_pos_vals)) *
                 " / " * fx(quantile_of(C.theta_pos_vals, 0.75)) * " / " *
                 fx(maximum(C.theta_pos_vals)) * "\n")
        both("    certified eps_max > 0 (SECONDARY): " * string(C.certified) * " / " *
             string(C.designs) * ", certificate at eps = 0.3 on " * string(C.cert_eps03) * "\n")
        isempty(C.min_q_feas) ||
            both("    min q at the feasible minimisers, median : " * sci(median_of(C.min_q_feas)) * "\n")
        both("    referee |Theta_exact - bisection| <= 1e-5 on " * string(C.referee_agree) *
             " / " * string(C.referee_n) * ", worst " * sci(C.worst_referee_gap) * "\n")
        both("    binding contact at the Theta_max = 0 designs:")
        for k in sort!(collect(keys(C.binding)))
            both("  " * k * " " * string(C.binding[k]))
        end
        both("\n    wall : " * fx(C.secs, 1) * " s\n\n")
    end
    free_pos = 0; free_n = 0; fixed_pos = 0; fixed_n = 0
    for (k, C) in cell
        if startswith(k, "free|")
            free_pos += C.theta_pos; free_n += C.designs
        end
        if startswith(k, "fixed|")
            fixed_pos += C.theta_pos; fixed_n += C.designs
        end
    end
    both("--- PASS rule (B4): exact Theta_max > 0 on >= 10 % of the population under `free` ---\n")
    both("free  : $free_pos / $free_n  (" * fx(free_n != 0 ? 100.0 * free_pos / free_n : 0.0, 2) * " %)\n")
    both("fixed : $fixed_pos / $fixed_n   (the K6 control)\n")
    both("VERDICT: " * ((free_n != 0 && 100.0 * free_pos / free_n >= 10.0) ?
                        "PASS -- the boundary is the culprit" : "FAIL -- B4 is dead") * "\n")
    both("wall " * fx(s(wall), 1) * " s\n")
    close(sm)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
