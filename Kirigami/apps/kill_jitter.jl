# A3 -- the JITTER TRANSITION (ideas/round2_adversary.md, Idea 4). Port of
# code/apps/kill_jitter.cpp.
#
# The hole in F25/F30 that the adversary names (rejection R-1): the emptiness result is
# "0/200 random graphs, 8/8 authored tilings". A referee cannot tell "the design space is
# empty" from "our generator and our orientation heuristic fail on our random graphs",
# because random Voronoi / Delaunay patches differ from authored tilings in EVERY respect
# at once -- edge lengths, vertex figures, face sizes, split-forest component length.
#
# This driver removes the confound by interpolating with the combinatorics HELD FIXED.
# Each of the eight authored reference cases has its INTERIOR vertex positions perturbed,
#
#     X[v] <- X[v] + a * (median edge length) * N(0, I),      v interior,
#
# for `a` on 20 log-spaced steps in [0.005, 1.0], 20 seeds each = 3200 designs per sigma
# rule, plus one a = 0 baseline per tiling. Faces, sigma rule and boundary are untouched,
# so the only thing that moves is geometry. Boundary vertices are NOT jittered because the
# pipeline solves Eq. (4)/(6) with Fixed boundary, i.e. they are the data the projection
# is anchored to.
#
# Per design and per sigma rule the existing pipeline is run unchanged:
#   sigma -> make_cut -> holes_partition -> assemble_system(Fixed) -> solve_system
#         -> X0 (Eq. (6) projection)
#         -> validity_certificate(POS /\ NOOVERLAP(eps/2) /\ NOROOT(eps)) at TWO eps
#         -> exact Theta_max by exact_theta_max_overlap (the T4.2" interval scan K2a used).
#
# Two sigma rules, as in K5:
#   sigma_mc  : the reference case's own orientation (checkerboard where the tiling is
#               2-colourable, else assign_orientation_relaxation = Eq. (1) max-cut). This
#               is a function of the DUAL GRAPH only, so it is constant along the ladder --
#               which is the point: the combinatorics are not randomised.
#   sigma_def : greedy flip search minimising the deployability defect D(sigma) at the
#               JITTERED X_ini, exactly K5's objective and moves, started from sigma_mc,
#               with a smaller flip cap because 3200 designs x 8 tilings is the budget.
#
# PASS BAR (restated from ideas/round2_adversary.md Idea 4, "Kill rule"):
#   PASS iff a sharp transition amplitude a* exists per tiling (certified fraction
#        decreasing in a, crossing 1/2 at a finite a*) AND some candidate predictor
#        reaches ROC AUC >= 0.9 for "certified Theta_max > 0".
#   FAIL iff the certified fraction is already ~0 at the smallest amplitude (the tilings
#        are knife-edge, so emptiness is generic and not a generator artefact), or no
#        monotone collapse, or every predictor has AUC < 0.9.
# Either outcome is informative; the adversary's own text says the knife-edge outcome is
# "a different and possibly better paper".
#
# The predictors scored (all computable from (G, sigma, X_ini) BEFORE any null space, as
# the idea demands; the two X0 ones are marked and reported separately as a control):
#   rho_geom  median over split-forest components of (geometric diameter of the component's
#             vertex set) / (median face inradius, 2A/P)
#   rho_hop   median over split-forest components of the graph-hop diameter
#   psplit    |E_split| / |F|                               (2026 Sec. 4.2's only proxy)
#   badq_ini  #{split edges with q_e <= 0 at X_ini} / |E_split|   (idea 3's N_bad)
#   defect    D(sigma) at X_ini, normalised by |holes| * med^4
#   aniso     median over faces of (longest edge / shortest edge) -- the "generic edge
#             lengths" the referee names, made into a number
#   amp       the jitter amplitude itself (a positive control: it MUST score well, and it
#             is not a predictor of anything since it is the independent variable)
#   badq_X0   [X0-derived control] same as badq_ini but at the projected X0
#   margin_X0 [X0-derived control] min corner margin mu at X0 / med^2
#
# Sharding: --shard k --nshards n takes every design whose index is == k (mod n), and
# writes jitter_<k>.csv. The rows are self-describing so `cat` in any order is valid.
#
#   julia --project=Kirigami Kirigami/apps/kill_jitter.jl [--out DIR] [--shard S]
#         [--nshards M] [--cap K] [--tiling T] [--eps 0.006] [--analyze merged.csv]
#         [--limit K] [--regenerate]
#
# Population: the 8 reference cases (frozen reference_cases_8; `--regenerate` rebuilds
# them). Outputs go to results/kill/jitter_julia/ (the C++ wrote results/kill/jitter/).
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# ---------------------------------------------------------------------------
# The ladder.
const kNAmp = 20
const kNSeed = 20
const kAmpLo = 0.005
const kAmpHi = 1.0

function amplitude(i::Int)  # i in [0, kNAmp), log-spaced; i < 0 means the a = 0 baseline
    i < 0 && return 0.0
    return kAmpLo * K.libm_pow(kAmpHi / kAmpLo, i / (kNAmp - 1))
end

# Gaussian displacement of the INTERIOR vertices only. Deterministic in (tiling, i, seed).
function jitter(m::K.Mesh, a::Float64, med::Float64, seed::UInt32)
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
# Predictors from (G, sigma, X) alone.

# 2 * area / perimeter: the inradius of a tangential polygon, and the standard proxy
# otherwise. Uses |signed area| so an inverted face still contributes a positive length.
function face_inradius(m::K.Mesh, X::Vector{Vec2}, f::Int)
    fv = m.faces[f]
    A = 0.0; P = 0.0
    nf = length(fv)
    for i in 1:nf
        p = X[fv[i]]
        q = X[fv[mod1(i + 1, nf)]]
        A += p[1] * q[2] - q[1] * p[2]
        P += K._norm2(q - p)
    end
    A = abs(0.5 * A)
    return P > 0 ? 2 * A / P : 0.0
end

function median_of(v::Vector{Float64})
    isempty(v) && return NaN
    sv = sort(v)
    return sv[length(sv) ÷ 2 + 1]
end

Base.@kwdef mutable struct ForestStats
    rho_geom::Float64 = 0.0   # median over components of geom diameter / median inradius
    rho_hop::Float64 = 0.0    # median over components of hop diameter
    max_geom::Float64 = 0.0   # worst component, in inradius units
    n_components::Int = 0     # components with at least one split edge
end

# Components of the subgraph of M spanned by E_split. Isolated vertices are ignored:
# the idea's rho is about the components a split cut actually creates.
function forest_stats(m::K.Mesh, c::K.CutStructure, X::Vector{Vec2})
    fs = ForestStats()
    N = K.n_vertices(m)
    adj = [Int[] for _ in 1:N]
    for e in c.split_edges
        k = m.edges[e].key
        push!(adj[k.a], k.b)
        push!(adj[k.b], k.a)
    end
    inr = [face_inradius(m, X, f) for f in 1:K.n_faces(m)]
    med_in = median_of(inr)

    seen = falses(N)
    geom = Float64[]; hop = Float64[]
    for s0 in 1:N
        (seen[s0] || isempty(adj[s0])) && continue
        comp = Int[]
        q = Int[s0]
        seen[s0] = true
        while !isempty(q)
            v = popfirst!(q)
            push!(comp, v)
            for w in adj[v]
                if !seen[w]
                    seen[w] = true
                    push!(q, w)
                end
            end
        end
        dmax = 0.0
        for i in 1:length(comp), j in i+1:length(comp)
            dmax = max(dmax, K._norm2(X[comp[i]] - X[comp[j]]))
        end
        # hop diameter by BFS from every member (components are small)
        hmax = 0
        for r in comp
            dist = fill(-1, N)
            bq = Int[r]
            dist[r] = 0
            while !isempty(bq)
                v = popfirst!(bq)
                for w in adj[v]
                    if dist[w] < 0
                        dist[w] = dist[v] + 1
                        push!(bq, w)
                    end
                end
            end
            for v in comp
                hmax = max(hmax, dist[v])
            end
        end
        push!(geom, med_in > 0 ? dmax / med_in : 0.0)
        push!(hop, Float64(hmax))
        fs.n_components += 1
    end
    fs.rho_geom = isempty(geom) ? 0.0 : median_of(geom)
    fs.rho_hop = isempty(hop) ? 0.0 : median_of(hop)
    fs.max_geom = isempty(geom) ? 0.0 : maximum(geom)
    return fs
end

# median over faces of (longest edge / shortest edge)
function face_anisotropy(m::K.Mesh, X::Vector{Vec2})
    r = Float64[]
    for f in 1:K.n_faces(m)
        fv = m.faces[f]
        lo = 1e300; hi = 0.0
        nf = length(fv)
        for i in 1:nf
            L = K._norm2(X[fv[mod1(i + 1, nf)]] - X[fv[i]])
            lo = min(lo, L)
            hi = max(hi, L)
        end
        lo > 0 && push!(r, hi / lo)
    end
    return median_of(r)
end

# ---------------------------------------------------------------------------
# K5's defect objective and flip search, trimmed for the ladder budget.
Base.@kwdef mutable struct Eval
    ok::Bool = false
    D::Float64 = 0.0
    n_split::Int = 0
end

# the C++ `Mesh m = m0; m.X = X; m.sigma = sigma; m.build_topology()`
function mesh_with(m0::K.Mesh, X::Vector{Vec2}, sigma::Vector{Int})
    m = K.Mesh(X, m0.faces)
    m.sigma = sigma
    m.periodic = m0.periodic
    K.build_topology!(m)
    return m
end

function eval_sigma(m0::K.Mesh, X::Vector{Vec2}, sigma::Vector{Int})
    e = Eval()
    m = mesh_with(m0, X, sigma)
    c = K.make_cut(m)
    K.build_hinge_graph(c).components != 1 && return e
    hs = K.holes_partition(c)
    r = K.hole_residuals(c, m.X, hs)
    d = 0.0
    for v in r.per_hole
        d += K._sqnorm2(v)
    end
    e.D = d
    e.n_split = K.n_split(c)
    e.ok = true
    return e
end

# Greedy single-face + adjacent-pair flips, accept strictly improving, start = sigma_mc.
# K5 used 4 starts and a 20*|F| cap; here 1 start and `cap` attempts, because the ladder
# is 3200 designs. The comparison that matters is sigma_def vs sigma_mc on the SAME design.
# Returns (sigma, D).
function defect_search(m0::K.Mesh, X::Vector{Vec2}, sigma_mc::Vector{Int}, cap::Int, seed::UInt32)
    rng = K.MT19937(seed)
    F = K.n_faces(m0)
    sig = copy(sigma_mc)
    cur = eval_sigma(m0, X, sig)
    cur.ok || return sigma_mc, NaN
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
            e = eval_sigma(m0, X, sig)
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
            e = eval_sigma(m0, X, sig)
            if e.ok && e.D < cur.D - 1e-12 * max(1.0, cur.D)
                cur = e
                improved = true
            else
                sig[f] = -sig[f]
                sig[g] = -sig[g]
            end
        end
    end
    return sig, cur.D
end

# ---------------------------------------------------------------------------
# One (design, sigma) measurement.
Base.@kwdef mutable struct Row
    solved::Bool = false
    n_split::Int = 0; H::Int = 0; dim_null::Int = 0; n_inv::Int = 0; c_gamma::Int = 0
    D::Float64 = 0.0                 # defect at X_ini
    proj_rel::Float64 = 0.0          # ||X0 - X_ini||_inf / med
    back_rel::Float64 = 0.0          # ||X0 - X_base||_inf / med, X_base = the UNJITTERED tiling
    theta_max::Float64 = -1.0        # exact, T4.2" interval scan
    zero_range::Bool = false
    cert_small::Bool = false; cert_large::Bool = false
    cert_pos::Bool = false; cert_noovl::Bool = false; cert_noroot::Bool = false
    rho_geom::Float64 = 0.0; rho_hop::Float64 = 0.0; rho_max::Float64 = 0.0; psplit::Float64 = 0.0
    badq_ini::Float64 = 0.0; aniso::Float64 = 0.0
    badq_X0::Float64 = 0.0; margin_X0::Float64 = 0.0
    n_forest_comp::Int = 0
end

linf(v::Vec2) = max(abs(v[1]), abs(v[2]))

function measure(m0::K.Mesh, X_ini::Vector{Vec2}, X_base::Vector{Vec2}, sigma::Vector{Int},
                 med::Float64, eps_small::Float64, eps_large::Float64)
    r = Row()
    m = mesh_with(m0, X_ini, sigma)
    c = K.make_cut(m)
    r.c_gamma = K.build_hinge_graph(c).components
    r.n_split = K.n_split(c)
    hs = K.holes_partition(c)
    let res = K.hole_residuals(c, m.X, hs)
        r.H = length(res.per_hole)
        for v in res.per_hole
            r.D += K._sqnorm2(v)
        end
    end
    # predictors from (G, sigma, X_ini)
    fs = forest_stats(m, c, X_ini)
    r.rho_geom = fs.rho_geom
    r.rho_hop = fs.rho_hop
    r.rho_max = fs.max_geom
    r.n_forest_comp = fs.n_components
    r.psplit = K.n_faces(m) != 0 ? r.n_split / K.n_faces(m) : 0.0
    r.aniso = face_anisotropy(m, X_ini)
    let q = K.zero_plus_q(c, X_ini)
        bad = count(v -> v <= 0, q)
        r.badq_ini = isempty(q) ? 0.0 : bad / length(q)
    end

    sr = K.solve_system(K.assemble_system(c, hs, m.X, K.Fixed), m.X)
    sr.projection_ok || return r
    r.solved = true
    r.dim_null = sr.dim_null
    X0 = K.matrix_to_points(sr.X0)
    mx = 0.0
    for i in eachindex(X0)
        mx = max(mx, linf(X0[i] - X_ini[i]))
    end
    r.proj_rel = med > 0 ? mx / med : mx
    # How far the Eq. (6) projection lands from the ORIGINAL tiling. If this is ~0 while
    # proj_rel is ~a, the projection has simply UNDONE the jitter and the design is not on
    # the ladder at all -- the single most important control in this experiment.
    bk = 0.0
    for i in eachindex(X0)
        bk = max(bk, linf(X0[i] - X_base[i]))
    end
    r.back_rel = med > 0 ? bk / med : bk
    r.n_inv = count_inverted(m, X0)
    let q = K.zero_plus_q(c, X0)
        bad = count(v -> v <= 0, q)
        r.badq_X0 = isempty(q) ? 0.0 : bad / length(q)
        mu = K.zero_plus_corner_margin(c, X0)
        mn = 1e300
        for v in mu
            mn = min(mn, v)
        end
        r.margin_X0 = isempty(mu) ? 0.0 : mn / (med * med)
    end

    # The certificate and the exact range, on the same basis and the same broad phase
    # K2a used. A basis mismatch means deploy() disagrees with the trig-linear form, i.e.
    # the configuration is not uniformly deployable; the row is then left unsolved.
    B = K.deploy_basis(c, X0)
    let err = 0.0, sc = 1e-300
        for th in (0.4, 1.5, 2.8)
            Y = K.deploy(c, X0, th).Y
            Yb = K.basis_eval(B, th)
            for i in eachindex(Y)
                err = max(err, K._norm2(Y[i] - Yb[i]))
                sc = max(sc, K._norm2(Y[i]))
            end
        end
        if err > 1e-9 * sc
            r.solved = false
            return r
        end
    end
    sd = K.swept_discs(c, B)
    pairs = K.candidate_pairs(c, sd, Float64(pi), true)
    cs = K.validity_certificate(c, B, X0, pairs, eps_small)
    cl = K.validity_certificate(c, B, X0, pairs, eps_large)
    r.cert_small = K.valid(cs)
    r.cert_large = K.valid(cl)
    r.cert_pos = cs.pos
    r.cert_noovl = cs.nooverlap
    r.cert_noroot = cs.noroot
    ov = K.exact_theta_max_overlap(c, B, pairs, 1e-9, Float64(pi), 1e-9)
    r.theta_max = ov.theta_max
    r.zero_range = ov.zero_range
    return r
end

# ---------------------------------------------------------------------------
# ANALYSIS (--analyze <merged csv>).

struct Table
    header::Vector{String}
    rows::Vector{Vector{String}}
end
function tcol(t::Table, name::String)
    i = findfirst(==(name), t.header)
    i === nothing && error("no column " * name)
    return i
end

function read_csv(path::String)
    isfile(path) || error("cannot open " * path)
    header = String[]
    rows = Vector{String}[]
    first = true
    for line in readlines(path)
        isempty(line) && continue
        cells = String.(split(line, ','))
        if first
            header = cells
            first = false
        elseif cells[1] != "tiling"
            push!(rows, cells)
        end
    end
    return Table(header, rows)
end

# Mann-Whitney ROC AUC with mid-ranks for ties. Returns NaN if either class is empty.
function auc(score::Vector{Float64}, label::Vector{Int})
    ord = sortperm(score; alg = MergeSort)   # std::sort tie order is irrelevant: mid-ranks
    rank = zeros(length(score))
    i = 1
    while i <= length(ord)
        j = i
        while j + 1 <= length(ord) && score[ord[j + 1]] == score[ord[i]]
            j += 1
        end
        r = 0.5 * ((i - 1) + (j - 1)) + 1.0
        for k in i:j
            rank[ord[k]] = r
        end
        i = j + 1
    end
    sum_pos = 0.0
    np = 0; nn = 0
    for k in eachindex(label)
        if label[k] != 0
            sum_pos += rank[k]
            np += 1
        else
            nn += 1
        end
    end
    (np == 0 || nn == 0) && return NaN
    return (sum_pos - np * (np + 1) / 2.0) / (Float64(np) * nn)
end

# The transition amplitude: the largest ladder amplitude at which the fraction is still
# >= 1/2, linearly interpolated in log(a) towards the next step. NaN if the fraction is
# never >= 1/2 (already dead at the smallest amplitude) or never < 1/2 (never dies).
# Returns (a_star, never_dies, dead_at_start).
function a_star(amps::Vector{Float64}, frac::Vector{Float64})
    isempty(amps) && return NaN, false, false
    frac[1] < 0.5 && return NaN, false, true
    for i in 1:length(amps)-1
        if frac[i] >= 0.5 && frac[i + 1] < 0.5
            l0 = log(amps[i]); l1 = log(amps[i + 1])
            w = (frac[i] - 0.5) / max(1e-12, frac[i] - frac[i + 1])
            return exp(l0 + w * (l1 - l0)), false, false
        end
    end
    return NaN, true, false
end

pad_to(s::String, n::Int) = s * " "^max(1, n - length(s))

function run_analysis(csv_path::String, outdir::String)
    t = read_csv(csv_path)
    c_til = tcol(t, "tiling"); c_ai = tcol(t, "amp_idx"); c_rule = tcol(t, "rule")
    c_solved = tcol(t, "solved"); c_nsplit = tcol(t, "n_split")
    c_cs = tcol(t, "cert_small"); c_th = tcol(t, "theta_max")
    c_pos = tcol(t, "cert_pos"); c_ovl = tcol(t, "cert_noovl"); c_root = tcol(t, "cert_noroot")

    preds = ["rho_geom", "rho_hop", "rho_max", "psplit", "badq_ini", "aniso", "D", "amp",
             "badq_X0", "margin_X0", "back_rel", "proj_rel"]
    pcol = [tcol(t, p) for p in preds]

    tilings = sort!(unique([r[c_til] for r in t.rows]))

    sm = open(joinpath(outdir, "summary.txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("A3 -- jitter transition. $(length(t.rows)) (design, sigma) rows from " * csv_path * "\n\n")

    # ---- per tiling x rule: the ladder ----
    lad = open(joinpath(outdir, "ladder.csv"), "w")
    print(lad, "tiling,rule,amp_idx,amp,n,frac_cert,frac_theta_pos,median_theta,n_split\n")
    both("=== the ladder: fraction with a POSITIVE exact Theta_max, and fraction CERTIFIED " *
         "(eps = 0.006) ===\n")
    both("(a = 0 is the unjittered authored tiling; amplitudes are in median edge lengths)\n\n")

    stars = NamedTuple[]
    for tl in tilings, rule in ("mc", "def")
        amps = Float64[]; fc = Float64[]; ft = Float64[]
        nsplit0 = -1
        for ai in -1:kNAmp-1
            n = 0; nc = 0; nt = 0
            ths = Float64[]
            for r in t.rows
                (r[c_til] != tl || r[c_rule] != rule || parse(Int, r[c_ai]) != ai) && continue
                n += 1
                if parse(Int, r[c_solved]) != 0
                    nc += parse(Int, r[c_cs])
                    th = parse(Float64, r[c_th])
                    nt += (th > 1e-9)
                    push!(ths, th)
                end
                nsplit0 < 0 && (nsplit0 = parse(Int, r[c_nsplit]))
            end
            n == 0 && continue
            f_c = nc / n; f_t = nt / n
            print(lad, tl, ",", rule, ",", ai, ",", cpp_g(amplitude(ai)), ",", n, ",",
                  cpp_g(f_c), ",", cpp_g(f_t), ",", cpp_g(median_of(ths)), ",", nsplit0, "\n")
            if ai >= 0
                push!(amps, amplitude(ai)); push!(fc, f_c); push!(ft, f_t)
            end
        end
        a_cert, nd_c, dead_c = a_star(amps, fc)
        a_theta, nd_t, dead_t = a_star(amps, ft)
        push!(stars, (tiling = tl, rule = rule, a_cert = a_cert, a_theta = a_theta,
                      nd_c = nd_c, dead_c = dead_c, nd_t = nd_t, dead_t = dead_t))
    end
    close(lad)

    star_str(v, nd, dead) = dead ? "DEAD at a=0.005" : (nd ? "survives a=1.0" : fx(v, 4))
    both("tiling                 rule   a*(Theta_max > 0)      a*(certified)\n")
    for st in stars
        both("  " * pad_to(st.tiling, 23) * st.rule * "    " *
             pad_to(star_str(st.a_theta, st.nd_t, st.dead_t), 23) *
             star_str(st.a_cert, st.nd_c, st.dead_c) * "\n")
    end
    both("\n")

    # ---- AUC ----
    function score_auc(subset::String, split_only::Bool)
        both("=== ROC AUC, population = " * subset * " ===\n")
        both("predictor      AUC(Theta_max>0)  AUC(certified)   [directed = max(AUC, 1-AUC)]\n")
        for p in eachindex(preds)
            sc_t = Float64[]; sc_c = Float64[]
            lb_t = Int[]; lb_c = Int[]
            for r in t.rows
                parse(Int, r[c_solved]) == 0 && continue
                (split_only && parse(Int, r[c_nsplit]) == 0) && continue
                v = parse(Float64, r[pcol[p]])
                isfinite(v) || continue
                push!(sc_t, v)
                push!(lb_t, Int(parse(Float64, r[c_th]) > 1e-9))
                push!(sc_c, v)
                push!(lb_c, parse(Int, r[c_cs]))
            end
            at = auc(sc_t, lb_t); ac = auc(sc_c, lb_c)
            d(x) = isnan(x) ? x : max(x, 1 - x)
            both("  " * pad_to(preds[p], 14) * fx(at, 4) * " [" * fx(d(at), 4) * "]     " *
                 fx(ac, 4) * " [" * fx(d(ac), 4) * "]\n")
        end
        both("\n")
    end
    score_auc("all solved rows (both sigma rules)", false)
    score_auc("split-bearing tilings only (n_split > 0)", true)

    # ---- which clause binds ----
    let n = 0, npos = 0, novl = 0, nroot = 0, nzr = 0
        for r in t.rows
            parse(Int, r[c_solved]) == 0 && continue
            n += 1
            npos += parse(Int, r[c_pos])
            novl += parse(Int, r[c_ovl])
            nroot += parse(Int, r[c_root])
            nzr += (parse(Float64, r[c_th]) <= 1e-9)
        end
        both("=== which certificate clause binds (all solved rows, n = $n) ===\n")
        both("  POS        $npos\n")
        both("  NOOVERLAP  $novl\n")
        both("  NOROOT     $nroot\n")
        both("  exact Theta_max == 0 : $nzr\n\n")
    end
    close(sm)
    println("wrote ", outdir, "/summary.txt and ", outdir, "/ladder.csv")
    return 0
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    outdir = joinpath(REPO, "results", "kill", "jitter_julia")
    shard = 0; nshards = 1; cap = 0; only_tiling = -1
    analyze = ""
    eps_small = 0.006; eps_large = 0.3
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshards" && i < length(args); nshards = arg_i(args[i+1]); i += 2
        elseif a == "--cap" && i < length(args); cap = arg_i(args[i+1]); i += 2
        elseif a == "--tiling" && i < length(args); only_tiling = arg_i(args[i+1]); i += 2
        elseif a == "--eps" && i < length(args); eps_small = arg_f(args[i+1]); i += 2
        elseif a == "--analyze" && i < length(args); analyze = args[i+1]; i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    isempty(analyze) || return run_analysis(analyze, outdir)
    csv = open(joinpath(outdir, "jitter_" * string(shard) * ".csv"), "w")
    print(csv, "tiling,amp_idx,amp,seed,N,F,med_edge,rule,",
          "solved,c_gamma,n_split,H,dim_null,n_inv,D,proj_rel,",
          "cert_small,cert_large,cert_pos,cert_noovl,cert_noroot,theta_max,zero_range,back_rel,",
          "rho_geom,rho_hop,rho_max,n_forest_comp,psplit,badq_ini,aniso,badq_X0,margin_X0,secs\n")

    cases = reference_cases(regenerate = regenerate)
    wall = Timer()
    idx = 0
    done = 0

    g_ = cpp_g
    function emit(tiling::String, ai::Int, a::Float64, seed::Int, m::K.Mesh, med::Float64,
                  rule::String, r::Row, secs::Float64)
        print(csv, tiling, ",", ai, ",", g_(a), ",", seed, ",", K.n_vertices(m), ",",
              K.n_faces(m), ",", g_(med), ",", rule, ",", Int(r.solved), ",",
              r.c_gamma, ",", r.n_split, ",", r.H, ",", r.dim_null, ",",
              r.n_inv, ",", g_(r.D), ",", g_(r.proj_rel), ",", Int(r.cert_small), ",",
              Int(r.cert_large), ",", Int(r.cert_pos), ",", Int(r.cert_noovl), ",",
              Int(r.cert_noroot), ",", g_(r.theta_max), ",", Int(r.zero_range), ",",
              g_(r.back_rel), ",",
              g_(r.rho_geom), ",", g_(r.rho_hop), ",", g_(r.rho_max), ",", r.n_forest_comp,
              ",", g_(r.psplit), ",", g_(r.badq_ini), ",", g_(r.aniso), ",", g_(r.badq_X0),
              ",", g_(r.margin_X0), ",", g_(secs), "\n")
        flush(csv)
    end

    for (ti, rc) in enumerate(cases)
        t = ti - 1   # the C++ tiling index (0-based, also in the seed)
        (only_tiling >= 0 && t != only_tiling) && continue
        base = rc.mesh
        K.build_topology!(base)
        sigma_mc = copy(base.sigma)
        med = median_edge_length(base)
        flip_cap = cap > 0 ? cap : 3 * K.n_faces(base)

        for ai in -1:kNAmp-1
            a = amplitude(ai)
            nseeds = ai < 0 ? 1 : kNSeed
            for seed in 0:nseeds-1
                idx += 1
                ((idx - 1) % nshards) == shard || continue
                done >= limit && break
                tm = Timer()
                sd = (UInt32(4200000) + UInt32(100000) * UInt32(t) +
                      UInt32(1000) * UInt32(ai + 1) + UInt32(seed)) % UInt32
                X = jitter(base, a, med, sd)

                rmc = measure(base, X, base.X, sigma_mc, med, eps_small, eps_large)
                emit(rc.name, ai, a, seed, base, med, "mc", rmc, s(tm))

                td = Timer()
                sigma_def, Ddef = defect_search(base, X, sigma_mc, flip_cap, (sd + UInt32(77)) % UInt32)
                rdef = measure(base, X, base.X, sigma_def, med, eps_small, eps_large)
                emit(rc.name, ai, a, seed, base, med, "def", rdef, s(td))

                done += 1
                done % 25 == 0 && println("shard ", shard, ": ", done, " designs, ", fx(s(wall), 1), " s")
            end
            done >= limit && break
        end
        done >= limit && break
    end
    close(csv)
    println("shard ", shard, " done: ", done, " designs, ", fx(s(wall), 1), " s")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
