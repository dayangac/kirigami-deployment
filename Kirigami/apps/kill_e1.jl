# E1 -- exact Theta_max (T4.2") vs the bisection referee, extended to >= 500 designs
# (specs/experimenter_final.md E1). Port of code/apps/kill_e1.cpp. Extends K2a
# (results/kill/KILL_REPORT.md, section K2a; 187 tightly embedded designs, all
# authored-tiling-derived) with three more populations, so the >=500 designs actually span
# random Voronoi/Delaunay/quad graphs and both orientation rules, not only the authored
# family K2a was restricted to:
#
#   (A) authored: kill_common's deployable_population(), called twice with a larger
#       sample count and two different jitter radii than K2a used (K2a: (8, 0.2);
#       here: (20, 0.2) and (20, 0.35)) -- "more seeds", same generator, same filter
#       (embedded at theta=0, uniformly deployable). sigma is whatever
#       deployable_population assigned (checkerboard on 2-colourable families,
#       Eq. (1) relaxation otherwise) -- this is the "mc" rule. A "def" row is added
#       per design by re-orienting via the K5 defect search (defect_search_at below,
#       copied from kill_k5 / kill_jitter -- not exposed by any header) at the
#       SAME X and re-solving; kept only if the def projection is also tightly
#       embedded.
#   (B) random Voronoi / Delaunay / quad-random graphs, make_graph(id, 100, 800, 1400)
#       -- K5's exact population call, the one measured to give sigma_def a
#       61.5% embedding rate (results/core_validation/referee_fix.md). Both sigma_mc
#       (Eq. (1), already assigned by make_graph) and sigma_def (K5's greedy defect
#       search) are solved and kept if the projection is tightly embedded.
#   (C) the jitter ladder (ideas/round2_adversary.md Idea 4 / kill_jitter) on the
#       four SPLIT-BEARING authored tilings only (hexagons_auto, truncated_square_488,
#       snub_square_33434, tiling_3_4_3_12 -- reference_cases() indices 3,4,5,6; STATE.md
#       F30/A3's a* = 0.28, 0.32, 0.16, 0.52 median edges resp.), amplitudes log-spaced
#       in [0.01, 0.6] (below and around a*, where designs are still embeddable often
#       enough to be worth the compute), both sigma rules.
#
# All three populations are filtered the same way K2a's `deployable_population` filters
# its own: uniformly deployable (Eq. (2) residual < 1e-8), no inverted face, and NOT
# colliding at theta = 0 at the tight tolerance 1e-12 -- "tightly embedded", matching
# K2a's own description of its population.
#
# Per accepted design: exact Theta_max (T4.2", exact_theta_max_overlap), the bisection
# referee (4000-point grid + 50 bisection steps, collision.jl's F34-fixed predicate,
# tol_rel = 1e-12 -- referee_fix.md: this predicate is now tolerance-independent over
# 1e-12..1e-6, so a single tolerance is reported, not three), the certificate columns
# (POS / NOOVERLAP(eps/2) / NOROOT at eps = 0.006, matching K2a) and whether the design
# ACTUALLY has Theta_max >= eps (by the exact scan) -- the certificate-exactness column.
#
# Sharding: --shard k --nshards n, exactly kill_jitter's convention (every task whose
# running counter is == k mod n is processed; the counter always advances, so `cat`ting
# the shards in any order reproduces the whole run).
#
#   julia --project=Kirigami Kirigami/apps/kill_e1.jl [--shard S] [--nshards M] [--out DIR]
#         [--n-random 900] [--limit K] [--regenerate]
#
# Population: (A) is always rebuilt (the frozen deployable_population.json is the (8, 0.2)
# variant); (B) reads the frozen e1_900 rows (ids 0..899, `--regenerate` rebuilds them
# with make_graph); (C) reads the frozen reference cases. E1's own sigma_def (cap 20*F,
# seed 5000000 + id) is recomputed here -- it is NOT the archived K5 sigma_def of
# e1_900.json (4 starts, K5's seed). `--limit K` stops after K tasks owned by this shard.
# Outputs go to results/final/e1_julia/ (the C++ wrote results/final/e1/).
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# ---------------------------------------------------------------------------
# K5's defect objective and greedy flip search (kill_k5 / kill_jitter). Semantics
# unchanged: D(sigma) = sum of squared per-hole Eq. (2) residuals at a given X; accept
# single-face and adjacent-pair flips that strictly improve D, reject any flip that
# disconnects the hinge graph.
Base.@kwdef mutable struct SigmaEval
    ok::Bool = false
    D::Float64 = 0.0
end

# the C++ `Mesh m = m0; m.X = X; m.sigma = sigma; m.build_topology()`
function mesh_with(m0::K.Mesh, X::Vector{Vec2}, sigma::Vector{Int})
    m = K.Mesh(X, m0.faces)
    m.sigma = sigma
    m.periodic = m0.periodic
    K.build_topology!(m)
    return m
end

function eval_sigma_at(m0::K.Mesh, X::Vector{Vec2}, sigma::Vector{Int})
    e = SigmaEval()
    m = mesh_with(m0, X, sigma)
    c = K.make_cut(m)
    K.build_hinge_graph(c).components != 1 && return e
    hs = K.holes_partition(c)
    r = K.hole_residuals(c, m.X, hs)
    for v in r.per_hole
        e.D += K._sqnorm2(v)
    end
    e.ok = true
    return e
end

function defect_search_at(m0::K.Mesh, X::Vector{Vec2}, sigma_mc::Vector{Int}, cap::Int, seed::UInt32)
    rng = K.MT19937(seed)
    F = K.n_faces(m0)
    sig = copy(sigma_mc)
    cur = eval_sigma_at(m0, X, sig)
    cur.ok || return sigma_mc
    adj = K.dual_graph(m0)
    dual_edges = Tuple{Int,Int}[]
    for f in 1:F, g in adj[f]
        g > f && push!(dual_edges, (f, g))
    end
    order = collect(1:F)

    attempts = 0
    improved = true
    while improved && attempts < cap
        improved = false
        K.shuffle!(order, rng)
        for f in order
            attempts >= cap && break
            attempts += 1
            sig[f] = -sig[f]
            e = eval_sigma_at(m0, X, sig)
            if e.ok && e.D < cur.D - 1e-12 * max(1.0, cur.D)
                cur = e
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
            e = eval_sigma_at(m0, X, sig)
            if e.ok && e.D < cur.D - 1e-12 * max(1.0, cur.D)
                cur = e
                improved = true
            else
                sig[f] = -sig[f]
                sig[g] = -sig[g]
            end
        end
    end
    return sig
end

# ---------------------------------------------------------------------------
# kill_jitter's amplitude ladder and interior-vertex Gaussian jitter.
function jitter_interior(m::K.Mesh, a::Float64, med::Float64, seed::UInt32)
    X = copy(m.X)
    a <= 0 && return X
    rng = K.MT19937(seed)
    G = K.NormalDist(0.0, 1.0)
    for v in 1:K.n_vertices(m)
        m.vertex_is_boundary[v] && continue
        g1 = K.normal(G, rng)
        g2 = K.normal(G, rng)
        X[v] += Vec2(g1, g2) * (a * med)
    end
    return X
end

# ---------------------------------------------------------------------------
# One CSV row: everything E1.md needs, per design.
Base.@kwdef mutable struct Row
    name::String = ""; family::String = ""; source::String = ""; sigma_rule::String = ""
    N::Int = 0; F::Int = 0; n_split::Int = 0
    theta_exact::Float64 = -1.0; theta_bisect::Float64 = -1.0; gap::Float64 = -1.0
    cert_valid::Bool = false; cert_pos::Bool = false; cert_nooverlap::Bool = false; cert_noroot::Bool = false
    actually_ge_eps::Bool = false
    secs::Float64 = 0.0
end

function write_header(csv::IO)
    print(csv, "name,family,source,sigma_rule,N,F,n_split,theta_exact,theta_bisect,gap,",
          "cert_valid,cert_pos,cert_nooverlap,cert_noroot,actually_ge_eps,secs\n")
end

function write_row(csv::IO, r::Row)
    print(csv, r.name, ",", r.family, ",", r.source, ",", r.sigma_rule, ",", r.N,
          ",", r.F, ",", r.n_split, ",", cpp_g(r.theta_exact), ",", cpp_g(r.theta_bisect),
          ",", cpp_g(r.gap), ",", Int(r.cert_valid), ",", Int(r.cert_pos), ",",
          Int(r.cert_nooverlap), ",", Int(r.cert_noroot), ",", Int(r.actually_ge_eps),
          ",", cpp_g(r.secs), "\n")
    flush(csv)
end

# Referee: 4000-point grid + 50 bisection steps, F34-fixed predicate (tolerance-
# independent; a single tol_rel = 1e-12 is used, matching K2a's "reference" column).
function bisect_theta_max(c::K.CutStructure, X::Vector{Vec2}, tol_rel::Float64 = 1e-12,
                          grid::Int = 4000, iters::Int = 50)
    col(th) = K.has_collision(c, K.deploy(c, X, th).Y, tol_rel)
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

const kEps = 0.006  # K2a's certificate eps

# Runs the full E1 measurement (exact scan, referee, certificate) on an already-embedded
# (mesh, X) pair and writes one row. Returns 1 if a row was written.
function measure_and_emit(csv::IO, name::String, family::String, source::String,
                          sigma_rule::String, m::K.Mesh, X::Vector{Vec2})
    t = Timer()
    K.build_topology!(m)
    c = K.make_cut(m)
    B = K.deploy_basis(c, X)
    let err = 0.0, sc = 1e-300  # guard: trig-linear form must reproduce deploy()
        for th in (0.4, 1.5, 2.8)
            Y = K.deploy(c, X, th).Y
            Yb = K.basis_eval(B, th)
            for i in eachindex(Y)
                err = max(err, K._norm2(Y[i] - Yb[i]))
                sc = max(sc, K._norm2(Y[i]))
            end
        end
        err > 1e-9 * sc && return 0  # SKIP: basis mismatch
    end
    sd = K.swept_discs(c, B)
    pairs = K.candidate_pairs(c, sd, Float64(pi), true)
    ov = K.exact_theta_max_overlap(c, B, pairs, 1e-9, Float64(pi), 1e-9)
    b = bisect_theta_max(c, X)
    cert = K.validity_certificate(c, B, X, pairs, kEps)

    r = Row(name = name, family = family, source = source, sigma_rule = sigma_rule,
            N = K.n_vertices(m), F = K.n_faces(m), n_split = K.n_split(c),
            theta_exact = ov.theta_max, theta_bisect = b, gap = abs(ov.theta_max - b),
            cert_valid = K.valid(cert), cert_pos = cert.pos, cert_nooverlap = cert.nooverlap,
            cert_noroot = cert.noroot, actually_ge_eps = (ov.theta_max >= kEps - 1e-9))
    r.secs = s(t)
    write_row(csv, r)
    return 1
end

# Tightly-embedded test matching deployable_population's own filter: uniformly
# deployable at X (residual < 1e-8), no inverted face, no collision at theta = 0,
# tol_rel = 1e-12 ("tightly embedded", K2a's phrase).
function tightly_embedded(c::K.CutStructure, m::K.Mesh, hs::K.HoleSet, X::Vector{Vec2})
    K.deployable(K.hole_residuals(c, X, hs), 1e-8) || return false
    count_inverted(m, X) != 0 && return false
    K.has_collision(c, K.deploy(c, X, 0.0).Y, 1e-12) && return false
    return true
end

# the sigma_def re-orientation of a (mesh, X): re-solve and emit if tightly embedded
function emit_def(csv::IO, m::K.Mesh, X::Vector{Vec2}, sigma_def::Vector{Int}, name::String,
                  family::String, source::String)
    md = K._with_sigma(m, sigma_def)
    K.build_topology!(md)
    cd = K.make_cut(md)
    hsd = K.holes_partition(cd)
    sr = K.solve_system(K.assemble_system(cd, hsd, X, K.Fixed), X)
    sr.projection_ok || return 0
    X0d = K.matrix_to_points(sr.X0)
    tightly_embedded(cd, md, hsd, X0d) || return 0
    return measure_and_emit(csv, name, family, source, "def", md, X0d)
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    shard = 0; nshards = 1
    n_random_ids = 900      # ids fed to make_graph, id % 3 selects the kind
    outdir = joinpath(REPO, "results", "final", "e1_julia")
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshards" && i < length(args); nshards = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--n-random" && i < length(args); n_random_ids = arg_i(args[i+1]); i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    csv = open(joinpath(outdir, "e1_" * string(shard) * ".csv"), "w")
    write_header(csv)
    wall = Timer()
    idx = 0; kept = 0; seen = 0

    owns() = (idx += 1; (idx - 1) % nshards == shard)
    done() = seen >= limit

    # (A) authored population, twice, with more samples/seeds than K2a's default (8, 0.2).
    for (rf, tag) in ((0.2, "a1"), (0.35, "a2"))
        done() && break
        for d in deployable_population(regenerate = true, samples_per_base = 20, radius_frac = rf)
            done() && break
            name = tag * "_" * d.name
            if owns()
                seen += 1
                kept += measure_and_emit(csv, name, d.family, "authored", "mc", d.mesh, d.X)
            end
            # sigma_def variant at the SAME X: re-orient, re-check tight embedding.
            if owns()
                seen += 1
                m = d.mesh
                K.build_topology!(m)
                sigma_mc = copy(m.sigma)
                sigma_def = defect_search_at(m, d.X, sigma_mc, 3 * K.n_faces(m),
                                             (UInt32(900001) + UInt32(seen)) % UInt32)
                kept += emit_def(csv, m, d.X, sigma_def, name * "_def", d.family, "authored")
            end
        end
    end

    # (B) random Voronoi / Delaunay / quad-random graphs, K5's exact population call.
    kinds = ("voronoi", "delaunay", "quad_random")
    rows = load_population("e1_900")
    for id in 0:n_random_ids-1
        done() && break
        g = if regenerate
            K.make_graph(id, 100, 800, 1400)
        else
            r = population_row(rows, id)
            K.Graph(r.mesh, r.kind, r.id, r.ok)
        end
        if !g.ok
            owns(); owns()   # keep idx in step whether or not it built
            continue
        end
        m = g.mesh
        K.build_topology!(m)
        sigma_mc = copy(m.sigma)
        name = "rand_" * kinds[id % 3 + 1] * "_" * string(id)

        if owns()
            seen += 1
            c = K.make_cut(m)
            hs = K.holes_partition(c)
            sr = K.solve_system(K.assemble_system(c, hs, m.X, K.Fixed), m.X)
            if sr.projection_ok
                X0 = K.matrix_to_points(sr.X0)
                tightly_embedded(c, m, hs, X0) && (kept += measure_and_emit(csv, name, g.kind, "random", "mc", m, X0))
            end
        end
        if owns()
            seen += 1
            sigma_def = defect_search_at(m, m.X, sigma_mc, 20 * K.n_faces(m),
                                         (UInt32(5000000) + UInt32(id)) % UInt32)
            kept += emit_def(csv, m, m.X, sigma_def, name, g.kind, "random")
        end
        (seen % 200 == 0 && seen > 0) && println(stderr, "shard ", shard, ": ", seen,
                                                 " tasks seen, ", kept, " kept, ", cpp_g(s(wall)), " s")
    end

    # (C) jitter ladder on the 4 split-bearing authored tilings (reference_cases indices
    # 3,4,5,6), amplitudes log-spaced below/around the known a* (STATE.md F30/A3).
    let cases = reference_cases(regenerate = regenerate)
        split_bearing = (3, 4, 5, 6)   # C++ 0-based indices
        kNAmp = 12; kNSeed = 15
        kAmpLo = 0.01; kAmpHi = 0.6
        amp(i) = kAmpLo * K.libm_pow(kAmpHi / kAmpLo, i / (kNAmp - 1))
        for si in 0:3
            done() && break
            rc = cases[split_bearing[si + 1] + 1]
            base = rc.mesh
            K.build_topology!(base)
            sigma_mc = copy(base.sigma)
            med = median_edge_length(base)
            tname = rc.name

            for ai in -1:kNAmp-1
                done() && break
                a = ai < 0 ? 0.0 : amp(ai)
                nseeds = ai < 0 ? 1 : kNSeed
                for seed in 0:nseeds-1
                    done() && break
                    sd = (UInt32(7300000) + UInt32(100000) * UInt32(si) + UInt32(1000) * UInt32(ai + 1) +
                          UInt32(seed)) % UInt32
                    X = jitter_interior(base, a, med, sd)
                    name = tname * "_a" * string(ai) * "_s" * string(seed)

                    if owns()
                        seen += 1
                        c = K.make_cut(base)
                        hs = K.holes_partition(c)
                        sr = K.solve_system(K.assemble_system(c, hs, X, K.Fixed), X)
                        if sr.projection_ok
                            X0 = K.matrix_to_points(sr.X0)
                            tightly_embedded(c, base, hs, X0) && (kept += measure_and_emit(csv, name, tname, "jitter", "mc", base, X0))
                        end
                    end
                    if owns()
                        seen += 1
                        sigma_def = defect_search_at(base, X, sigma_mc, 3 * K.n_faces(base), (sd + UInt32(77)) % UInt32)
                        kept += emit_def(csv, base, X, sigma_def, name * "_def", tname, "jitter")
                    end
                end
            end
        end
    end
    close(csv)

    println(stderr, "shard ", shard, " done: ", seen, " tasks seen, ", kept, " kept, ", cpp_g(s(wall)), " s")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
