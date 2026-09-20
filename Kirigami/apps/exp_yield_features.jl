# exp_yield_features -- INPUT-SIDE features of the K9c population, for WP2 (why the K9c yield is
# not scale-free).
#
# REPORT.md Sec. Results records that the K9c method reaches exact Theta_max > 0 on
# 307/400 designs but that the rate falls with |F| (Delaunay ~1.0 -> ~0.8, quad-random
# 0.75 -> 0.43 across the bins [100, 230, 420, 800]) and differs by family. The open
# question ("Next steps" 7) is WHICH property of the input decides it: split-cut density,
# aspect ratio, or the solver's basin.
#
# This driver answers the measurement half. It computes, for every design, a vector of
# features that depend ONLY on the input -- the graph, the orientation sigma, the input
# embedding X_ini and the cut structure derived from them -- plus a clearly separated
# group of PROJECTION-side quantities (things the Eq. (6) projection produces before any
# optimiser runs). Nothing the K9c solver returns is a feature; those live in the
# outcome columns and are excluded from the pre-registered bar by specs/m2_experimenter_yield.md.
#
# Two modes:
#   --mode k9c    ids 0..199 of the K9 population (make_graph(id, 100, 800, 1400)), both
#                 sigma rules, sigma_def read from the results/experiments/k5/sigma cache and the
#                 Eq. (6) projection read from the results/experiments/k6/cache shape cache, i.e.
#                 EXACTLY the objects exp_k9c_range_embedding.jl uses. No solver is run: the outcomes are
#                 joined in from results/experiments/k9c/k9c.csv on (id, kind, sigma).
#   --mode fresh  ids 1000..1099, held out. sigma_mc from make_graph, sigma_def from
#                 method::orientation_defect (the library promotion of K5's search, cap
#                 20|F|, seed 7000 + id), then method::design_range_max at the K9c run's
#                 defaults with seed 9300 + 7 id + which. Features and outcomes both new.
#   --mode check  confirms that the 100 fresh graphs are not any of the 200 K9 graphs.
#
# Rows are flushed as they are produced, so a killed shard leaves usable output.
#
#   julia --project=Kirigami Kirigami/apps/exp_yield_features.jl [--mode k9c] [--out DIR]
#         [--sigma DIR] [--cache DIR] [--k9c PATH] [--shard S] [--nshards M] [--n 100]
#         [--id0 1000] [--limit K] [--regenerate]
#
# Populations: --mode k9c walks data/corpus/native200.json (the K9 population with the
# archived K5 sigma_def), --mode fresh data/corpus/yield_fresh_100.json (whose sigma_def
# is the archived results/yield/fresh_sigma orientation); `--regenerate` rebuilds the
# graphs through make_graph and recomputes the fresh sigma_def with orientation_defect.
# One deliberate deviation in --mode k9c: the Eq. (6) projection comes from the shape
# cache when present (`--cache`) and is otherwise SOLVED and cached, where the original
# K9c run skipped the design (this repo carries no
# results/experiments/k6/cache). Outputs go to results/yield/ by default. `--limit K` stops
# after K graphs of the shard.
include(joinpath(@__DIR__, "exp_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

median_of(v::Vector{Float64}) = isempty(v) ? NaN : sort(v)[length(v) ÷ 2 + 1]

# ---------------------------------------------------------------------------
# The feature vector. Column order here is the column order of the CSV; both modes
# emit the identical block so features.csv and fresh.csv are directly comparable.
Base.@kwdef mutable struct Feat
    # size / topology
    N::Int = 0; F::Int = 0; E::Int = 0; E_int::Int = 0; n_border::Int = 0
    n_split::Int = 0; n_hinge::Int = 0; dim_null::Int = 0
    split_density::Float64 = 0.0      # n_split / E_interior
    split_per_face::Float64 = 0.0     # n_split / F
    dim_null_per_N::Float64 = 0.0

    # the 0+ constraint count that range_embed actually enforces
    n_corner::Int = 0                 # face corners (convexity constraints)
    n_incid::Int = 0                  # corner incidences (mu constraints)
    n_constraints_0plus::Int = 0      # n_split + n_incid
    dof_ratio::Float64 = 0.0          # 2 dim_null / n_constraints_0plus
    dof_ratio_all::Float64 = 0.0      # 2 dim_null / (n_split + n_incid + n_corner)

    # face geometry at X_ini
    ang_min::Float64 = 0.0; ang_med::Float64 = 0.0             # corner angles, degrees
    asp_min::Float64 = 0.0; asp_med::Float64 = 0.0; asp_max::Float64 = 0.0   # longest/shortest face edge
    frac_face_sharp::Float64 = 0.0    # faces with a corner < 20 deg
    face_deg_med::Float64 = 0.0

    # cut structure
    n_split_comp::Int = 0; max_split_comp::Int = 0; split_depth::Int = 0
    mean_split_comp::Float64 = 0.0
    hinge_diam::Int = 0
    hinge_frac_largest::Float64 = 0.0
    pure_frac::Float64 = 0.0          # interior vertices touched by no split edge

    # hole structure
    H::Int = 0; n_notch::Int = 0
    hole_size_mean::Float64 = 0.0
    hole_size_max::Int = 0

    # the 0+ signs AT X_ini
    frac_q_pos_ini::Float64 = 0.0; min_q_ini::Float64 = 0.0
    frac_mu_pos_ini::Float64 = 0.0; min_mu_ini::Float64 = 0.0
    frac_convex_ini::Float64 = 0.0

    # --- PROJECTION side (produced by the Eq. (6) solve, before any optimiser) ---
    proj_dist_x0::Float64 = 0.0       # ||X0 - X_ini||_F / (sqrt(N) med), per-vertex RMS
    proj_frac_convex_x0::Float64 = 0.0
    proj_nonconvex0::Int = 0
    proj_frac_q_pos_x0::Float64 = 0.0
    proj_min_q_x0::Float64 = 0.0
end

const kFeatHeader =
    "id,kind,sigma,N,F,E,E_int,n_border,n_split,n_hinge,dim_null,split_density," *
    "split_per_face,dim_null_per_N,n_corner,n_incid,n_constraints_0plus,dof_ratio," *
    "dof_ratio_all,ang_min,ang_med,asp_min,asp_med,asp_max,frac_face_sharp,face_deg_med," *
    "n_split_comp,max_split_comp,mean_split_comp,split_depth,hinge_diam," *
    "hinge_frac_largest,pure_frac,H,n_notch,hole_size_mean,hole_size_max," *
    "frac_q_pos_ini,min_q_ini,frac_mu_pos_ini,min_mu_ini,frac_convex_ini," *
    "proj_dist_x0,proj_frac_convex_x0,proj_nonconvex0,proj_frac_q_pos_x0,proj_min_q_x0," *
    "med_edge"

function write_feat(o::IO, id::Int, kind::AbstractString, sig::AbstractString, f::Feat, med::Float64)
    g = fmt_g
    print(o, id, ",", kind, ",", sig, ",", f.N, ",", f.F, ",", f.E, ",", f.E_int, ",", f.n_border,
          ",", f.n_split, ",", f.n_hinge, ",", f.dim_null, ",", g(f.split_density), ",",
          g(f.split_per_face), ",", g(f.dim_null_per_N), ",", f.n_corner, ",", f.n_incid, ",",
          f.n_constraints_0plus, ",", g(f.dof_ratio), ",", g(f.dof_ratio_all), ",", g(f.ang_min),
          ",", g(f.ang_med), ",", g(f.asp_min), ",", g(f.asp_med), ",", g(f.asp_max), ",",
          g(f.frac_face_sharp), ",", g(f.face_deg_med), ",", f.n_split_comp, ",", f.max_split_comp,
          ",", g(f.mean_split_comp), ",", f.split_depth, ",", f.hinge_diam, ",",
          g(f.hinge_frac_largest), ",", g(f.pure_frac), ",", f.H, ",", f.n_notch, ",",
          g(f.hole_size_mean), ",", f.hole_size_max, ",", g(f.frac_q_pos_ini), ",", g(f.min_q_ini),
          ",", g(f.frac_mu_pos_ini), ",", g(f.min_mu_ini), ",", g(f.frac_convex_ini), ",",
          g(f.proj_dist_x0), ",", g(f.proj_frac_convex_x0), ",", f.proj_nonconvex0, ",",
          g(f.proj_frac_q_pos_x0), ",", g(f.proj_min_q_x0), ",", g(med))
    return nothing
end

# Longest path (in edges) of a tree component, by double BFS.
function tree_diameter(adj::Vector{Vector{Int}}, comp::Vector{Int})
    length(comp) < 2 && return 0
    function bfs(s::Int)
        d = Dict{Int,Int}(s => 0)
        q = Int[s]
        best = s; bd = 0
        while !isempty(q)
            v = popfirst!(q)
            if d[v] > bd
                bd = d[v]; best = v
            end
            for w in adj[v]
                if !haskey(d, w) && w != s
                    d[w] = d[v] + 1
                    push!(q, w)
                end
            end
        end
        return bd, best
    end
    _, a = bfs(comp[1])
    bd, _ = bfs(a)
    return bd
end

# Exact diameter of the largest connected component of the hinge subgraph, by BFS from
# every vertex of that component.
function hinge_stats(m::K.Mesh, c::K.CutStructure)
    N = K.n_vertices(m)
    adj = [Int[] for _ in 1:N]
    for e in c.hinge_edges
        k = m.edges[e].key
        push!(adj[k.a], k.b)
        push!(adj[k.b], k.a)
    end
    comp = fill(-1, N)
    nc = 0; best = 0
    comps = Vector{Int}[]
    for v in 1:N
        (comp[v] >= 0 || isempty(adj[v])) && continue
        cur = Int[]
        q = Int[v]
        comp[v] = nc
        while !isempty(q)
            u = popfirst!(q)
            push!(cur, u)
            for w in adj[u]
                if comp[w] < 0
                    comp[w] = nc
                    push!(q, w)
                end
            end
        end
        (best == 0 || length(cur) > length(comps[best])) && (best = nc + 1)
        push!(comps, cur)
        nc += 1
    end
    best == 0 && return 0, 0.0
    L = comps[best]
    frac_largest = length(L) / max(1, N)
    d = fill(-1, N)
    dm = 0
    for s in L
        fill!(d, -1)
        q = Int[s]
        d[s] = 0
        while !isempty(q)
            u = popfirst!(q)
            dm = max(dm, d[u])
            for w in adj[u]
                if d[w] < 0
                    d[w] = d[u] + 1
                    push!(q, w)
                end
            end
        end
    end
    return dm, frac_largest
end

function features(m::K.Mesh, c::K.CutStructure, hs::K.HoleSet, X_ini::Vector{Vec2}, med::Float64)
    f = Feat()
    s = med * med
    f.N = K.n_vertices(m)
    f.F = K.n_faces(m)
    f.E = K.n_edges(m)
    f.n_split = K.n_split(c)
    f.n_hinge = K.n_hinge(c)
    f.n_border = length(c.border_edges)
    f.E_int = f.n_split + f.n_hinge
    f.split_density = f.n_split / max(1, f.E_int)
    f.split_per_face = f.n_split / max(1, f.F)

    inc = K.corner_incidences(c)
    fc = K.face_corners(m)
    f.n_corner = length(fc)
    f.n_incid = length(inc)
    f.n_constraints_0plus = f.n_split + f.n_incid
    f.dof_ratio = 0.0   # dim_null filled by the caller; see finish_dof

    # --- face geometry at X_ini -------------------------------------------------
    mi = K.Mesh(X_ini, m.faces)
    angs = Float64[]; asps = Float64[]; degs = Float64[]
    sharp = 0
    for fi in 1:K.n_faces(m)
        vs = m.faces[fi]
        push!(degs, Float64(length(vs)))
        amin = 1e9; lmin = 1e300; lmax = 0.0
        for i in 1:length(vs)
            a = K.corner_angle(mi, fi, i) * 180.0 / pi
            push!(angs, a)
            amin = min(amin, a)
            L = norm(X_ini[vs[mod1(i + 1, length(vs))]] - X_ini[vs[i]])
            lmin = min(lmin, L)
            lmax = max(lmax, L)
        end
        amin < 20.0 && (sharp += 1)
        push!(asps, lmax / max(1e-300, lmin))
    end
    f.ang_min = isempty(angs) ? 0.0 : minimum(angs)
    f.ang_med = median_of(angs)
    f.asp_min = isempty(asps) ? 0.0 : minimum(asps)
    f.asp_med = median_of(asps)
    f.asp_max = isempty(asps) ? 0.0 : maximum(asps)
    f.frac_face_sharp = sharp / max(1, f.F)
    f.face_deg_med = median_of(degs)

    # --- split forest ------------------------------------------------------------
    N = K.n_vertices(m)
    sadj = [Int[] for _ in 1:N]
    touched = falses(N)
    for e in c.split_edges
        k = m.edges[e].key
        push!(sadj[k.a], k.b)
        push!(sadj[k.b], k.a)
        touched[k.a] = true
        touched[k.b] = true
    end
    let seen = falses(N), sizes = Int[]
        for v in 1:N
            (seen[v] || !touched[v]) && continue
            cur = Int[]
            q = Int[v]
            seen[v] = true
            while !isempty(q)
                u = popfirst!(q)
                push!(cur, u)
                for w in sadj[u]
                    if !seen[w]
                        seen[w] = true
                        push!(q, w)
                    end
                end
            end
            push!(sizes, length(cur))
            f.split_depth = max(f.split_depth, tree_diameter(sadj, cur))
        end
        f.n_split_comp = length(sizes)
        f.max_split_comp = isempty(sizes) ? 0 : maximum(sizes)
        f.mean_split_comp = isempty(sizes) ? 0.0 : sum(Float64, sizes) / length(sizes)
    end
    let n_int = 0, pure = 0
        for v in 1:N
            if !m.vertex_is_boundary[v]
                n_int += 1
                touched[v] || (pure += 1)
            end
        end
        f.pure_frac = pure / max(1, n_int)
    end
    f.hinge_diam, f.hinge_frac_largest = hinge_stats(m, c)

    # --- holes -------------------------------------------------------------------
    f.H = K.n_interior_holes(hs)
    f.n_notch = 0
    let sz = Float64[]
        for h in hs.all
            if !h.all_interior
                f.n_notch += 1
                continue
            end
            push!(sz, Float64(length(h.edges)))
            f.hole_size_max = max(f.hole_size_max, length(h.edges))
        end
        f.hole_size_mean = isempty(sz) ? 0.0 : sum(sz) / length(sz)
    end

    # --- the 0+ signs AT X_ini ----------------------------------------------------
    let q = K.zero_plus_q(c, X_ini)
        pos = count(v -> v > 0, q)
        mn = isempty(q) ? 1e300 : minimum(q)
        f.frac_q_pos_ini = isempty(q) ? 1.0 : pos / length(q)
        f.min_q_ini = isempty(q) ? 1.0 : mn / s
    end
    let mu = K.zero_plus_corner_margin(c, X_ini)
        pos = count(v -> v > 0, mu)
        mn = isempty(mu) ? 1e300 : minimum(mu)
        f.frac_mu_pos_ini = isempty(mu) ? 1.0 : pos / length(mu)
        f.min_mu_ini = isempty(mu) ? 1.0 : mn / s
    end
    let cr = K.corner_crosses(m, X_ini)
        pos = count(v -> v > 0, cr)
        f.frac_convex_ini = isempty(cr) ? 1.0 : pos / length(cr)
    end
    return f
end

# The projection-side group, from X0.
function projection_features!(m::K.Mesh, c::K.CutStructure, X0::Vector{Vec2},
                              X_ini::Vector{Vec2}, med::Float64, f::Feat)
    s = med * med
    d2 = 0.0
    for v in eachindex(X0)
        d2 += dot(X0[v] - X_ini[v], X0[v] - X_ini[v])
    end
    f.proj_dist_x0 = sqrt(d2 / max(1, length(X0))) / med
    cr = K.corner_crosses(m, X0)
    pos = count(v -> v > 0, cr)
    bad = length(cr) - pos
    f.proj_frac_convex_x0 = isempty(cr) ? 1.0 : pos / length(cr)
    f.proj_nonconvex0 = bad
    q = K.zero_plus_q(c, X0)
    qpos = count(v -> v > 0, q)
    mn = isempty(q) ? 1e300 : minimum(q)
    f.proj_frac_q_pos_x0 = isempty(q) ? 1.0 : qpos / length(q)
    f.proj_min_q_x0 = isempty(q) ? 1.0 : mn / s
    return nothing
end

function finish_dof!(f::Feat)
    f.dim_null_per_N = f.dim_null / max(1, f.N)
    f.dof_ratio = 2.0 * f.dim_null / max(1, f.n_constraints_0plus)
    f.dof_ratio_all = 2.0 * f.dim_null / max(1, f.n_constraints_0plus + f.n_corner)
    return nothing
end

# ---------------------------------------------------------------------------
# k9c.csv, indexed by "id|kind|sigma".
Base.@kwdef mutable struct Outcome
    theta_exact::Float64 = 0.0; eps_max::Float64 = 0.0; best_margin::Float64 = 0.0
    dist::Float64 = 0.0; secs::Float64 = 0.0
    best_feas::Int = 0
    binding::String = ""; best_src::String = ""
    ok::Bool = false
end

function load_k9c(path::AbstractString)
    out = Dict{String,Outcome}()
    isfile(path) || return out
    lines = readlines(path)
    isempty(lines) && return out
    h = String.(split(lines[1], ','))
    col(nm) = something(findfirst(==(nm), h), 0)
    ci = col("id"); ck = col("kind"); cs = col("sigma")
    for ln in lines[2:end]
        isempty(ln) && continue
        v = String.(split(ln, ','))
        length(v) != length(h) && continue
        o = Outcome()
        o.ok = true
        o.theta_exact = parse(Float64, v[col("theta_exact")])
        o.eps_max = parse(Float64, v[col("eps_max")])
        o.best_margin = parse(Float64, v[col("best_margin")])
        o.dist = parse(Float64, v[col("dist")])
        o.secs = parse(Float64, v[col("secs")])
        o.best_feas = parse(Int, v[col("best_feas")])
        o.binding = v[col("binding")]
        o.best_src = v[col("best_src")]
        out[v[ci] * "|" * v[ck] * "|" * v[cs]] = o
    end
    return out
end

const kOutHeader = ",theta_exact,eps_max,best_feas,best_margin,best_src,binding,dist,secs"

function write_outcome(o::IO, u::Outcome)
    print(o, ",", fmt_g(u.theta_exact), ",", fmt_g(u.eps_max), ",", u.best_feas, ",",
          fmt_g(u.best_margin), ",", u.best_src, ",", u.binding, ",", fmt_g(u.dist), ",",
          fmt_g(u.secs))
    return nothing
end

# A cheap fingerprint of a graph, for the fresh-vs-K9 disjointness check.
function fingerprint(m::K.Mesh)
    sx = 0.0; sy = 0.0; sq = 0.0
    for p in m.X
        sx += p[1]; sy += p[2]; sq += dot(p, p)
    end
    return string(K.n_vertices(m), ":", K.n_faces(m), ":", sci(sx, 9), ":", sci(sy, 9), ":", sci(sq, 9))
end

# The mesh of a frozen row, rebuilt through make_graph on request.
function row_mesh(row::PopRow, regenerate::Bool)
    if regenerate
        g = K.make_graph(row.id, 100, 800, 1400)
        return g.ok ? g.mesh : nothing
    end
    return row.ok ? row.mesh : nothing
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    mode = "k9c"
    outdir = joinpath(REPO, "results", "yield")
    sigmadir = ""
    cache = joinpath(REPO, "results", "experiments", "k6", "cache")
    k9c_csv = joinpath(REPO, "results", "experiments", "k9c", "k9c.csv")
    shard = 0; nshards = 1; n_fresh = 100; id0 = 1000
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--mode" && i < length(args); mode = args[i+1]; i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--sigma" && i < length(args); sigmadir = args[i+1]; i += 2
        elseif a == "--cache" && i < length(args); cache = args[i+1]; i += 2
        elseif a == "--k9c" && i < length(args); k9c_csv = args[i+1]; i += 2
        elseif a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshards" && i < length(args); nshards = arg_i(args[i+1]); i += 2
        elseif a == "--n" && i < length(args); n_fresh = arg_i(args[i+1]); i += 2
        elseif a == "--id0" && i < length(args); id0 = arg_i(args[i+1]); i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)

    # ---------------------------------------------------------------------------
    if mode == "check"
        old = Dict{String,String}()
        for row in load_population("k1a_200")
            row.id < 200 || continue
            m = row_mesh(row, regenerate)
            m === nothing && continue
            K.build_topology!(m)
            old[fingerprint(m)] = row.kind * "_" * string(row.id)
        end
        o = open(joinpath(outdir, "fresh_check.txt"), "w")
        print(o, "fingerprint = N:F:sum_x:sum_y:sum_|x|^2 of the generated graph\n")
        print(o, "K9 population graphs fingerprinted: ", length(old), " (of 200 ids)\n")
        coll = 0; made = 0
        fresh = load_population("yield_fresh_100")
        for id in id0:id0+n_fresh-1
            ri = findfirst(r -> r.id == id, fresh)
            m = (ri === nothing || regenerate) ? (g = K.make_graph(id, 100, 800, 1400); g.ok ? g.mesh : nothing) :
                                                 row_mesh(fresh[ri], false)
            if m === nothing
                print(o, "id ", id, " FAILED to generate\n")
                continue
            end
            made += 1
            K.build_topology!(m)
            fp = fingerprint(m)
            if haskey(old, fp)
                coll += 1
                print(o, "id ", id, " COLLIDES with ", old[fp], "\n")
            end
        end
        print(o, "fresh graphs generated: ", made, " / ", n_fresh, "\n")
        print(o, "collisions with the K9 population: ", coll, "\n")
        close(o)
        println("wrote ", outdir, "/fresh_check.txt: ", made, " fresh graphs, ", coll, " collisions")
        return 0
    end

    # ---------------------------------------------------------------------------
    if mode == "k9c"
        k9c = load_k9c(k9c_csv)
        if isempty(k9c)
            println(stderr, "cannot read ", k9c_csv)
            return 1
        end
        csv = open(joinpath(outdir, "features.csv"), "w")
        print(csv, kFeatHeader, kOutHeader, "\n")
        rows = 0; missing = 0
        gidx = -1
        nrun = 0
        for row in load_population("native200")
            gidx + 1 < 200 || break
            m0 = row_mesh(row, regenerate)
            m0 === nothing && continue
            K.build_topology!(m0)
            sigma_mc = copy(m0.sigma)
            X_ini = copy(m0.X)
            id = row.id
            sigma_def = Int[]
            if !isempty(sigmadir)
                p = joinpath(sigmadir, row.kind * "_" * string(id) * ".json")
                if isfile(p)
                    sm = K.load_mesh_json(p)
                    (length(sm.sigma) == K.n_faces(m0) && length(sm.X) == length(m0.X)) && (sigma_def = sm.sigma)
                end
            elseif row.sigma_def !== nothing && length(row.sigma_def) == K.n_faces(m0)
                sigma_def = row.sigma_def
            end
            isempty(sigma_def) && continue
            gidx += 1
            gidx >= 200 && break
            nrun >= limit && break
            nrun += 1
            med = median_edge_length(m0)
            for which in 0:1
                fname = which == 1 ? "sigma_def" : "sigma_mc"
                m = K.Mesh(m0.X, m0.faces)
                m.sigma = copy(which == 1 ? sigma_def : sigma_mc)
                K.build_topology!(m)
                c = K.make_cut(m)
                K.n_split(c) == 0 && continue
                hs = K.holes_partition(c)
                sh = shape_space(m, c, hs, joinpath(cache, fname), id)
                if !sh.ok || sh.N != K.n_vertices(m) || sh.k < 1
                    missing += 1
                    println(stderr, "no shape for id ", id, " ", fname)
                    continue
                end
                u = get(k9c, string(id) * "|" * row.kind * "|" * fname, nothing)
                if u === nothing
                    missing += 1
                    continue
                end
                f = features(m, c, hs, X_ini, med)
                f.dim_null = sh.k
                projection_features!(m, c, K.matrix_to_points(sh.X0), X_ini, med, f)
                finish_dof!(f)
                write_feat(csv, id, row.kind, fname, f, med)
                write_outcome(csv, u)
                print(csv, "\n")
                flush(csv)
                rows += 1
            end
        end
        close(csv)
        println("wrote ", outdir, "/features.csv: ", rows, " rows, ", missing, " skipped")
        return 0
    end

    # ---------------------------------------------------------------------------
    if mode == "fresh"
        csv_path = nshards > 1 ? joinpath(outdir, "fresh_shard_" * string(shard) * ".csv") :
                                 joinpath(outdir, "fresh.csv")
        csv = open(csv_path, "w")
        print(csv, kFeatHeader, kOutHeader, "\n")
        mkpath(joinpath(outdir, "fresh_sigma"))
        rows = 0
        wall = Timer()
        fresh = load_population("yield_fresh_100")
        nrun = 0
        for k in 0:n_fresh-1
            (nshards > 1 && k % nshards != shard) && continue
            nrun >= limit && break
            id = id0 + k
            ri = findfirst(r -> r.id == id, fresh)
            row = ri === nothing ? nothing : fresh[ri]
            m0 = (row === nothing || regenerate) ? (g = K.make_graph(id, 100, 800, 1400); g.ok ? g.mesh : nothing) :
                                                    row_mesh(row, false)
            if m0 === nothing
                println(stderr, "id ", id, " generate FAILED")
                continue
            end
            nrun += 1
            K.build_topology!(m0)
            sigma_mc = copy(m0.sigma)
            X_ini = copy(m0.X)
            med = median_edge_length(m0)

            # sigma_def: K5's greedy defect search, cached so a restart does not redo it.
            # Between the run's own cache file and the search sits the frozen row's
            # archived sigma_def (skipped with --regenerate).
            sigma_def = Int[]
            kind = row === nothing ? K.make_graph(id, 100, 800, 1400).kind : row.kind
            sp = joinpath(outdir, "fresh_sigma", kind * "_" * string(id) * ".json")
            if isfile(sp)
                sm = K.load_mesh_json(sp)
                length(sm.sigma) == K.n_faces(m0) && (sigma_def = sm.sigma)
            end
            if isempty(sigma_def) && !regenerate && row !== nothing && row.sigma_def !== nothing &&
               length(row.sigma_def) == K.n_faces(m0)
                sigma_def = row.sigma_def
            end
            if isempty(sigma_def)
                dr = K.orientation_defect(m0, sigma_mc, 20 * K.n_faces(m0), 7000 + id)
                if dr.ok && !isempty(dr.sigma)
                    sigma_def = dr.sigma
                    sm = K.Mesh(m0.X, m0.faces)
                    sm.sigma = copy(sigma_def)
                    K.save_mesh_json(sm, sp)
                end
            end
            if isempty(sigma_def)
                println(stderr, "id ", id, " sigma_def FAILED")
                continue
            end

            for which in 0:1
                fname = which == 1 ? "sigma_def" : "sigma_mc"
                t = Timer()
                sig = which == 1 ? sigma_def : sigma_mc
                m = K.Mesh(m0.X, m0.faces)
                m.sigma = copy(sig)
                K.build_topology!(m)
                c = K.make_cut(m)
                K.n_split(c) == 0 && continue
                hs = K.holes_partition(c)

                ro = K.RangeMaxOptions()
                ro.seed = UInt32(9300) + UInt32(7) * UInt32(id) + UInt32(which)
                rr = K.design_range_max(m, sig, X_ini, ro)
                if !rr.design.ok
                    println(stderr, "id ", id, " ", fname, " design not ok: ", rr.design.status)
                    continue
                end
                f = features(m, c, hs, X_ini, med)
                f.dim_null = rr.design.dim_null
                projection_features!(m, c, rr.design.X0, X_ini, med, f)
                finish_dof!(f)

                u = Outcome()
                u.ok = true
                u.theta_exact = rr.design.ch.theta_max
                u.eps_max = rr.design.ch.eps_max
                u.best_feas = rr.design.feasible ? 1 : 0
                u.best_margin = rr.margin
                u.best_src = rr.provenance
                u.binding = rr.design.ch.binding
                u.dist = rr.design.dist_ini
                u.secs = s(t)

                write_feat(csv, id, kind, fname, f, med)
                write_outcome(csv, u)
                print(csv, "\n")
                flush(csv)
                rows += 1
                println("id ", id, " ", fname, " F=", f.F, " k=", f.dim_null, " theta=",
                        fx(u.theta_exact, 4), " eps=", fx(u.eps_max, 3), " src=", u.best_src, " (",
                        fx(u.secs, 1), " s, wall ", fx(s(wall), 0), ")")
                flush(stdout)
            end
        end
        close(csv)
        println("wrote ", csv_path, ": ", rows, " rows in ", fx(s(wall), 1), " s")
        return 0
    end

    println(stderr, "unknown --mode ", mode)
    return 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
