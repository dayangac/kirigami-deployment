# T-1 -- the vertex balance defect. Port of code/apps/kill_t1.cpp.
#
# The claim under test (ideas/round2_theorist.md Idea 1 / §0.c, "T-1"). At an interior
# vertex v of the PURE case the barycentric row of Eq. (2) reads
#
#     (BAL)   delta_v := sum_{u in In(v)} ( x_u - x_v )  =  0 ,
#
# where In(v) = { src(e) : e a hinge edge with dst(e) = v } (core/tutte_auxetic.jl
# assembles exactly this row per hole preimage, and by F11 a vertex touched by NO split
# edge is its own preimage, so its row IS (BAL)). T-1 reads (BAL) kinematically: by T1.5
# the k_v copies of v separate at 0+ with velocities proportional to J(x_u - x_v), one
# per in-edge; (BAL) says those velocities sum to zero, hence for k_v >= 3 they
# positively span the plane, hence -- T-1's second "hence" -- no copy of v can be pushed
# into a neighbouring face's material. A split edge at v deletes a term from (BAL)
# without deleting the corresponding motion, so the balance is broken there.
#
# theorist-b (ideas/round2_theorist_b.md, the logged disagreement) objects that the
# second "hence" does not follow: the 0+ corner predicate of method/zero_plus.jl is
#     mu = max(-g1, -g2)   at a convex corner,
# a DISJUNCTION, not a half-plane condition, so positively spanning velocities do not by
# themselves keep every copy out of every corner cone. ideas/ranking_r2.md §1 adjudicates
# the two and specifies this measurement:
#
#   over the K6 population (200 random graphs x 2 sigma = 400 designs, at the projected
#   design X0 -- the SAME point at which K6 reports inward0 / min_mu0), classify every
#   interior vertex by whether any incident edge is a split edge, and report
#     N_pure = # interior vertices touched by NO split edge,
#     N_bad  = # of those with mu_v <= 0,
#   where mu_v = min of method/zero_plus.jl's corner margin over the incidences at v.
#
#   PASS  N_bad = 0 (T-1's safe half is true as stated; the 0+ collapse is localised on
#         split-touched vertices and T-1 becomes a real lemma).
#   FAIL  N_bad > 0 (T-1 dead as stated; the corner disjunction is the mechanism).
#   INCONCLUSIVE if N_pure is too small to decide (< 50 pure vertices in total).
#
# Also reported, because they are free at the same X0:
#   * corr(|delta_v|, mu_v) over the SPLIT-TOUCHED vertices -- T-1's quantitative half
#     predicts that the larger the balance defect, the smaller the separation margin;
#   * per design, whether the culprit split edge (argmin_e q_e at X0, i.e. the first 0+
#     collision K6 reports through min_q0) is incident to the vertex of largest |delta_v|.
#
# Nothing is recomputed: X0 and Phi come from K6's shape-space cache and sigma_def from
# K5's saved orientations, so this driver reads exactly the designs K6 measured.
#
#   julia --project=Kirigami Kirigami/apps/kill_t1.jl [--n 200] [--maxf 800] [--out DIR]
#         [--sigma DIR] [--cache DIR] [--mu-tol 1e-12] [--delta-tol 1e-9] [--limit K]
#         [--regenerate]
#
# Population as kill_k9.jl (frozen k1a_200 + archived sigma_def; C++ shape-cache layout).
# Outputs go to results/kill/t1_julia/ (the C++ wrote results/kill/t1/).
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# Pearson correlation of two equal-length samples; NaN if either is constant.
function pearson(a::Vector{Float64}, b::Vector{Float64})
    n = length(a)
    (n != length(b) || n < 2) && return NaN
    ma = 0.0; mb = 0.0
    for i in 1:n
        ma += a[i]; mb += b[i]
    end
    ma /= n; mb /= n
    saa = 0.0; sbb = 0.0; sab = 0.0
    for i in 1:n
        da = a[i] - ma; db = b[i] - mb
        saa += da * da; sbb += db * db; sab += da * db
    end
    (saa <= 0 || sbb <= 0) && return NaN
    return sab / sqrt(saa * sbb)
end

# Fractional ranks (average rank for ties), for Spearman.
function ranks_of(v::Vector{Float64})
    n = length(v)
    idx = sortperm(v; alg = MergeSort)
    r = zeros(n)
    i = 1
    while i <= n
        j = i
        while j + 1 <= n && v[idx[j + 1]] == v[idx[i]]
            j += 1
        end
        avg = 0.5 * ((i - 1) + (j - 1)) + 1.0
        for k in i:j
            r[idx[k]] = avg
        end
        i = j + 1
    end
    return r
end

function spearman(a::Vector{Float64}, b::Vector{Float64})
    (length(a) != length(b) || length(a) < 2) && return NaN
    return pearson(ranks_of(a), ranks_of(b))
end

function median_of(v::Vector{Float64})
    isempty(v) && return NaN
    sv = sort(v)
    return sv[length(sv) ÷ 2 + 1]
end

Base.@kwdef mutable struct Totals
    designs::Int = 0
    n_int::Int = 0; n_pure::Int = 0; n_touched::Int = 0
    n_bad::Int = 0                       # pure vertices with mu_v <= 0
    n_pure_k3::Int = 0; n_bad_k3::Int = 0  # ... restricted to k_v >= 3 (T-1's hypothesis)
    n_pure_k2::Int = 0; n_bad_k2::Int = 0  # ... and to k_v = 2 (outside it)
    n_touched_bad::Int = 0               # touched vertices with mu_v <= 0
    n_pure_unbalanced::Int = 0           # pure vertices with |delta_v| above tolerance
    # Of the binding (minimising) incidence at a BAD pure vertex: is the corner convex?
    n_bad_convex::Int = 0; n_bad_reflex::Int = 0
    # pure vertices ALL of whose incidences sit at convex corners (T-1's best case)
    n_pure_allcvx::Int = 0; n_bad_allcvx::Int = 0
    n_faces_tot::Int = 0; n_faces_nonconvex_X0::Int = 0; n_faces_nonconvex_in::Int = 0
    designs_with_pure::Int = 0           # designs having at least one pure vertex
    designs_with_bad::Int = 0
    culprit_hits::Int = 0; culprit_n::Int = 0   # culprit edge incident to the argmax-|delta| vertex
    culprit_top10::Int = 0                       # ... or to a vertex in the top 10 % of |delta|
    corr_per_design::Vector{Float64} = Float64[]   # Spearman(|delta|, mu) on touched, per design
    dv_touched::Vector{Float64} = Float64[]; mu_touched::Vector{Float64} = Float64[]  # pooled
    pure_frac::Vector{Float64} = Float64[]
    max_pure_delta_rel::Float64 = 0.0
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
    n = 200; maxf = 800
    outdir = joinpath(REPO, "results", "kill", "t1_julia")
    sigmadir = ""
    cache = joinpath(REPO, "results", "kill", "k6", "cache")
    # mu and q are AREAS (a 2x2 determinant of two vectors linear in X), so the zero test
    # is relative to med^2; delta is a LENGTH, relative to med.
    mu_tol_rel = 1e-12; delta_tol_rel = 1e-9
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--maxf" && i < length(args); maxf = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--sigma" && i < length(args); sigmadir = args[i+1]; i += 2
        elseif a == "--cache" && i < length(args); cache = args[i+1]; i += 2
        elseif a == "--mu-tol" && i < length(args); mu_tol_rel = arg_f(args[i+1]); i += 2
        elseif a == "--delta-tol" && i < length(args); delta_tol_rel = arg_f(args[i+1]); i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)

    vcsv = open(joinpath(outdir, "t1.csv"), "w")
    print(vcsv, "id,kind,sigma,v,deg,k_in,n_split_inc,touched,delta_rel,mu_rel,n_inc,",
          "is_culprit_endpoint\n")
    gcsv = open(joinpath(outdir, "t1_graph.csv"), "w")
    print(gcsv, "id,kind,sigma,N,F,med_edge,n_split,n_int,n_pure,n_bad,n_touched,n_touched_bad,",
          "n_pure_unbalanced,n_pure_k3,n_bad_k3,n_pure_k2,n_bad_k2,max_pure_delta_rel,min_mu_pure,min_mu_touched,",
          "spearman_touched,pearson_touched,culprit_edge,culprit_q,culprit_is_argmax_delta,",
          "culprit_delta_pctile\n")

    T = Totals()
    graphs = 0; missing_sigma = 0
    wall = Timer()

    # The identical deterministic population index as kill_k6: a graph counts only if
    # it exists AND has a K5 sigma_def, and the ids are walked in order.
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
        nrun >= limit && break
        nrun += 1

        med = median_edge_length(m0)
        mu_tol = mu_tol_rel * med * med
        delta_tol = delta_tol_rel * med
        counted = false

        for which in 0:1
            fname = which == 1 ? "sigma_def" : "sigma_mc"
            m = K._with_sigma(m0, which == 1 ? sigma_def : sigma_mc)
            K.build_topology!(m)
            c = K.make_cut(m)
            K.n_split(c) == 0 && continue
            hs = K.holes_partition(c)
            sh = shape_space(m, c, hs, joinpath(cache, fname), id)
            (!sh.ok || sh.k < 1) && continue
            X0 = K.matrix_to_points(sh.X0)
            if !counted
                graphs += 1
                counted = true
            end
            T.designs += 1

            N = K.n_vertices(m)

            # ---- the balance defect delta_v, per interior vertex -------------------
            # delta_v = sum over hinge edges with dst == v of ( X0[src] - X0[v] ), i.e. the
            # residual of the barycentric row (BAL). Equal to the hole-preimage residual of
            # core/tutte_auxetic.jl with the sign flipped, and identically zero at a vertex
            # that is its own preimage (F11: not touched by any split edge).
            delta = fill(Vec2(0, 0), N)
            k_in = zeros(Int, N)
            for e in c.hinge_edges
                d = c.hinge_dir[e]
                delta[d.dst] += X0[d.src] - X0[d.dst]
                k_in[d.dst] += 1
            end

            # ---- split incidence per vertex ---------------------------------------
            n_split_inc = zeros(Int, N)
            for e in c.split_edges
                n_split_inc[m.edges[e].key.a] += 1
                n_split_inc[m.edges[e].key.b] += 1
            end

            # ---- mu_v = min corner margin over the incidences at v ------------------
            inc = K.corner_incidences(c)
            mu = K.zero_plus_corner_margin(c, X0)
            mu_v = fill(Inf, N)
            n_inc = zeros(Int, N)
            mu_v_convex = falses(N)   # convexity of the MINIMISING incidence
            all_convex = trues(N)     # every incidence at v sits at a convex corner
            for i in eachindex(inc)
                v = inc[i].v
                e1 = X0[inc[i].v_next] - X0[v]
                e2 = X0[inc[i].v_prev] - X0[v]
                cr = K._det2(e1, e2)
                if mu[i] < mu_v[v]
                    mu_v[v] = mu[i]
                    mu_v_convex[v] = cr > 0
                end
                (cr > 0) || (all_convex[v] = false)
                n_inc[v] += 1
            end

            # How much of the reflex population is projection damage: count faces that are
            # non-convex at X0 against the same count at the INPUT embedding m.X.
            for f in 1:K.n_faces(m)
                vs = m.faces[f]
                nf = length(vs)
                bad0 = false; badi = false
                for j in 1:nf
                    a_ = vs[j]; b_ = vs[mod1(j + 1, nf)]; cc = vs[mod1(j + 2, nf)]
                    p1 = X0[b_] - X0[a_]; p2 = X0[cc] - X0[b_]
                    K._det2(p1, p2) <= 0 && (bad0 = true)
                    r1 = m.X[b_] - m.X[a_]; r2 = m.X[cc] - m.X[b_]
                    K._det2(r1, r2) <= 0 && (badi = true)
                end
                T.n_faces_tot += 1
                bad0 && (T.n_faces_nonconvex_X0 += 1)
                badi && (T.n_faces_nonconvex_in += 1)
            end

            # ---- the culprit split edge: argmin_e q_e at X0 -------------------------
            q0 = K.zero_plus_q(c, X0)
            sc = K.split_copies(c)
            culprit = 0; cu = 0; cv = 0
            culprit_q = NaN
            for i in eachindex(q0)
                if culprit == 0 || q0[i] < culprit_q
                    culprit = sc[i].edge
                    culprit_q = q0[i]
                    cu = sc[i].v_from
                    cv = sc[i].v_to
                end
            end

            # ---- per-design accumulation -------------------------------------------
            n_int = 0; n_pure = 0; n_bad = 0; n_touched = 0; n_touched_bad = 0
            n_pure_unbal = 0
            n_pure_k3 = 0; n_bad_k3 = 0; n_pure_k2 = 0; n_bad_k2 = 0
            min_mu_pure = Inf
            min_mu_touched = Inf
            max_pure_delta = 0.0
            best_delta = -1.0
            argmax_delta = 0
            dv_t = Float64[]; mu_t = Float64[]
            all_delta = Float64[]

            for v in 1:N
                m.vertex_is_boundary[v] && continue
                n_inc[v] == 0 && continue  # fewer than two copies: no 0+ corner predicate
                n_int += 1
                dn = K._norm2(delta[v])
                drel = dn / med
                mrel = mu_v[v] / (med * med)
                touched = n_split_inc[v] > 0
                push!(all_delta, drel)
                if dn > best_delta
                    best_delta = dn
                    argmax_delta = v
                end
                if touched
                    n_touched += 1
                    mu_v[v] <= mu_tol && (n_touched_bad += 1)
                    min_mu_touched = min(min_mu_touched, mrel)
                    push!(dv_t, drel)
                    push!(mu_t, mrel)
                    push!(T.dv_touched, drel)
                    push!(T.mu_touched, mrel)
                else
                    n_pure += 1
                    bad = mu_v[v] <= mu_tol
                    if bad
                        n_bad += 1
                        if mu_v_convex[v]; T.n_bad_convex += 1; else; T.n_bad_reflex += 1; end
                    end
                    if all_convex[v]
                        T.n_pure_allcvx += 1
                        bad && (T.n_bad_allcvx += 1)
                    end
                    # T-1 assumes k_v >= 3 (three or more zero-sum velocities positively span the
                    # plane); k_v = 2 is outside the hypothesis and is counted separately.
                    if k_in[v] >= 3
                        n_pure_k3 += 1
                        bad && (n_bad_k3 += 1)
                    elseif k_in[v] == 2
                        n_pure_k2 += 1
                        bad && (n_bad_k2 += 1)
                    end
                    min_mu_pure = min(min_mu_pure, mrel)
                    dn > delta_tol && (n_pure_unbal += 1)
                    max_pure_delta = max(max_pure_delta, drel)
                end
                is_cul = (v == cu || v == cv)
                # vertex ids are written 0-based, as the C++ did
                print(vcsv, id, ',', row.kind, ',', fname, ',', v - 1, ',',
                      length(m.vertex_edges[v]), ',', k_in[v], ',', n_split_inc[v], ',',
                      touched ? 1 : 0, ',', sci(drel, 6), ',', sci(mrel, 6), ',',
                      n_inc[v], ',', is_cul ? 1 : 0, '\n')
            end

            sp = spearman(dv_t, mu_t)
            pe = pearson(dv_t, mu_t)

            # Is the culprit split edge incident to the largest-|delta| interior vertex, and
            # where do its endpoints sit in the |delta| distribution?
            cul_argmax = (culprit > 0 && argmax_delta > 0 &&
                          (cu == argmax_delta || cv == argmax_delta)) ? 1 : 0
            cul_pct = NaN
            if culprit > 0 && !isempty(all_delta)
                dc = max(K._norm2(delta[cu]), K._norm2(delta[cv])) / med
                below = count(x -> x <= dc, all_delta)
                cul_pct = below / length(all_delta)
            end

            # edge id written 0-based (-1 = none), as the C++ did
            print(gcsv, id, ',', row.kind, ',', fname, ',', N, ',', K.n_faces(m), ',',
                  fx(med, 6), ',', K.n_split(c), ',', n_int, ',', n_pure, ',',
                  n_bad, ',', n_touched, ',', n_touched_bad, ',',
                  n_pure_unbal, ',', n_pure_k3, ',', n_bad_k3, ',',
                  n_pure_k2, ',', n_bad_k2, ',', sci(max_pure_delta, 4), ',',
                  sci(min_mu_pure, 4), ',', sci(min_mu_touched, 4), ',',
                  fx(sp, 6), ',', fx(pe, 6), ',', culprit - 1, ',',
                  sci(culprit_q, 4), ',', cul_argmax, ',', fx(cul_pct, 4), '\n')

            T.n_int += n_int; T.n_pure += n_pure; T.n_touched += n_touched
            T.n_bad += n_bad; T.n_touched_bad += n_touched_bad
            T.n_pure_unbalanced += n_pure_unbal
            T.n_pure_k3 += n_pure_k3; T.n_bad_k3 += n_bad_k3
            T.n_pure_k2 += n_pure_k2; T.n_bad_k2 += n_bad_k2
            n_pure > 0 && (T.designs_with_pure += 1)
            n_bad > 0 && (T.designs_with_bad += 1)
            T.max_pure_delta_rel = max(T.max_pure_delta_rel, max_pure_delta)
            push!(T.pure_frac, n_int != 0 ? n_pure / n_int : 0.0)
            isfinite(sp) && push!(T.corr_per_design, sp)
            if culprit > 0
                T.culprit_n += 1
                T.culprit_hits += cul_argmax
                (isfinite(cul_pct) && cul_pct >= 0.9) && (T.culprit_top10 += 1)
            end
        end
        (gidx + 1) % 25 == 0 && println(stderr, "  ", gidx + 1, " graphs, ", T.designs,
                                        " designs, ", cpp_g(s(wall)), " s")
    end
    close(vcsv)
    close(gcsv)

    pooled_sp = spearman(T.dv_touched, T.mu_touched)
    pooled_pe = pearson(T.dv_touched, T.mu_touched)

    # The verdict is decided on the vertices T-1's hypothesis actually covers (k_v >= 3);
    # the k_v = 2 counts are reported alongside because a k_v = 2 vertex has only two
    # opposite velocities, which span a LINE, not the plane, so T-1 says nothing there.
    verdict = T.n_pure_k3 < 50 ? "INCONCLUSIVE" : (T.n_bad_k3 == 0 ? "PASS" : "FAIL")

    o = IOBuffer()
    print(o, "\n=== T-1: vertex balance defect (K6 population, at X0) ===\n")
    print(o, "graphs ", graphs, ", designs ", T.designs, " (missing sigma_def: ", missing_sigma, ")\n")
    print(o, "interior vertices with >=2 copies : ", T.n_int, "\n")
    print(o, "  pure    (no incident split edge): ", T.n_pure,
          "   on ", T.designs_with_pure, " / ", T.designs, " designs",
          ", median pure fraction ", fx(median_of(T.pure_frac), 4), "\n")
    print(o, "  touched (>=1 incident split edge): ", T.n_touched, "\n")
    print(o, "N_bad  (pure, mu_v <= 0)          : ", T.n_bad,
          "   on ", T.designs_with_bad, " / ", T.designs, " designs\n")
    print(o, "        (touched, mu_v <= 0)      : ", T.n_touched_bad, "\n")
    print(o, "  restricted to k_v >= 3 (T-1's hypothesis): pure ", T.n_pure_k3, ", bad ", T.n_bad_k3, "\n")
    print(o, "  restricted to k_v == 2 (outside it)      : pure ", T.n_pure_k2, ", bad ", T.n_bad_k2, "\n")
    print(o, "  binding incidence at a bad pure vertex: convex corner ", T.n_bad_convex,
          ", reflex corner ", T.n_bad_reflex, "\n")
    print(o, "  pure vertices with EVERY incidence at a convex corner: ", T.n_pure_allcvx,
          ", bad ", T.n_bad_allcvx, "\n")
    print(o, "  faces non-convex at X0: ", T.n_faces_nonconvex_X0, " / ", T.n_faces_tot,
          " ; at the INPUT embedding: ", T.n_faces_nonconvex_in, "\n")
    print(o, "collapse RATE: pure ", fx(T.n_pure != 0 ? T.n_bad / T.n_pure : 0.0, 4),
          "  vs touched ", fx(T.n_touched != 0 ? T.n_touched_bad / T.n_touched : 0.0, 4), "\n")
    print(o, "pure vertices with |delta_v| > tol : ", T.n_pure_unbalanced,
          "   (max |delta_v|/med = ", sci(T.max_pure_delta_rel, 3),
          "; must be ~0, it is the F11 self-check)\n")
    print(o, "corr(|delta_v|, mu_v) on touched   : Spearman ", fx(pooled_sp, 4),
          ", Pearson ", fx(pooled_pe, 4),
          ", median per-design Spearman ", fx(median_of(T.corr_per_design), 4),
          "  (n = ", length(T.dv_touched), ")\n")
    print(o, "culprit split edge incident to the argmax-|delta| vertex : ",
          T.culprit_hits, " / ", T.culprit_n,
          " ; to a top-10% |delta| vertex : ", T.culprit_top10, " / ", T.culprit_n, "\n")
    print(o, "VERDICT: ", verdict, "\n")
    print(o, "wall ", fx(s(wall), 1), " s\n")

    txt = String(take!(o))
    print(txt)
    open(joinpath(outdir, "t1_summary.txt"), "w") do f
        print(f, txt)
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
