# K1b -- the harmonic identity (ideas/ranking.md Sec. 4, K1b). Gates R1 and R2.
# Port of code/apps/kill_k1b.cpp.
#
# PASS rule, copied from ranking.md:
#   "200 graphs from K1a's set, 20 random X in X each. For 50 random ordered
#    (vertex, edge) triples per graph, sample deploy() at 200 angles on
#    (0, min(pi, theta_max)), least-squares fit the orientation determinant to
#    p + q cos(theta) + r sin(theta). PASS if the relative residual is < 1e-10 on
#    every triple of every graph. FAIL on any exceedance. Also assert face signed
#    areas are constant in theta to 1e-12."
#
# Deviation, recorded: the sampling interval is (0, pi] rather than
# (0, min(pi, theta_max)) -- theta_max is a collision quantity and the identity under
# test is algebraic, so the wider interval is a strictly harder test.
#
# Beyond the spec: the closed-form coefficients
#   p = (det(U,P)+det(V,Q))/2, q = (det(U,P)-det(V,Q))/2, r = (det(U,Q)+det(V,P))/2
# are compared with the fitted ones, and Y(theta) = cos(theta/2) C + sin(theta/2) S is
# compared with deploy() directly (claim U3/U6).
#
#   julia --project=Kirigami Kirigami/apps/kill_k1b.jl [--n 200] [--shapes 20] [--out DIR]
#         [--cache DIR] [--limit K] [--regenerate]
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# C++ std::to_string(double): fixed, 6 decimals
to_string_d(v::Real) = fx(v, 6)

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n = 200; n_shapes = 20; n_triples = 50; n_angles = 200; limit = typemax(Int)
    outdir = joinpath(REPO, "results", "kill", "k1b_julia")
    cache = joinpath(REPO, "results", "kill", "cache")
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i < length(args); n = arg_i(args[i+1]); i += 2
        elseif a == "--shapes" && i < length(args); n_shapes = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--cache" && i < length(args); cache = args[i+1]; i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    csv = open(joinpath(outdir, "k1b.csv"), "w")
    print(csv, "id,kind,N,F,dim_null,n_triples,max_rel_res,max_rel_res_geo,max_closedform_err,",
          "max_coeff_rel_err,max_basis_err,basis_scale,max_area_drift,area_scale,min_degeneracy,",
          "secs\n")
    wall = Timer()

    worst_res = 0.0; worst_coef = 0.0; worst_basis = 0.0; worst_area = 0.0; worst_res_geo = 0.0
    worst_cf = 0.0
    worst_res_graph = ""; worst_area_graph = ""; worst_res_geo_graph = ""
    n_triples_total = 0; n_fail = 0; n_fail_geo = 0; n_degen_total = 0
    worst_degeneracy = 1.0
    graphs = 0; minF = 1 << 30; maxF = 0

    function run_one(name::String, kind::String, m::K.Mesh, id::Int)
        t = Timer()
        K.build_topology!(m)
        isempty(m.sigma) && return
        c = K.make_cut(m)
        hs = K.holes_partition(c)
        sh = shape_space(m, c, hs, cache, id)
        sh.ok || return
        rng = K.MT19937(UInt32(7777) + UInt32(31) * UInt32(graphs))
        G = K.NormalDist(0.0, 1.0)
        s_scale = median_edge_length(m)

        g_res = 0.0; g_coef = 0.0; g_basis = 0.0; g_scale = 0.0; g_area = 0.0; g_ascale = 0.0
        g_res_geo = 0.0; g_cf = 0.0; g_deg = 1.0
        g_degen = 0; g_tr = 0
        lo = m.X[1]; hi = m.X[1]
        for p in m.X
            lo = min.(lo, p); hi = max.(hi, p)
        end
        g_scale_diam = norm(hi - lo)
        ths = [pi * (a + 1.0) / n_angles for a in 0:n_angles-1]
        F = K.n_faces(m)

        for _ in 1:n_shapes
            T = Matrix{Float64}(undef, sh.k, 2)
            for i in 1:sh.k   # row by row: T(i,0) then T(i,1), as the C++ draws them
                T[i, 1] = K.normal(G, rng); T[i, 2] = K.normal(G, rng)
            end
            Xm = copy(sh.X0)
            if sh.k > 0
                D = sh.Phi * T
                mx = maximum(abs, D)
                mx > 0 && (D .*= s_scale / mx)
                Xm .+= D
            end
            X = K.matrix_to_points(Xm)
            B = K.deploy_basis(c, X)

            # all deployments once, shared by every triple
            Ys = Vector{Vector{Vec2}}(undef, n_angles)
            for a in 1:n_angles
                Ys[a] = K.deploy(c, X, ths[a]).Y
                Yb = K.basis_eval(B, ths[a])
                for i in eachindex(Yb)
                    g_basis = max(g_basis, norm(Ys[a][i] - Yb[i]))
                    g_scale = max(g_scale, norm(Ys[a][i]))
                end
            end
            # face signed areas constant in theta
            Y0 = K.basis_eval(B, 0.0)
            for f in 1:F
                pf = c.prime_faces[f]
                area(Y) = begin
                    s = 0.0
                    for i in eachindex(pf)
                        u = Y[pf[i]]; v = Y[pf[mod1(i + 1, length(pf))]]
                        s += u[1] * v[2] - u[2] * v[1]
                    end
                    0.5 * s
                end
                a0 = area(Y0)
                g_ascale = max(g_ascale, abs(a0))
                for a in 1:17:n_angles
                    g_area = max(g_area, abs(area(Ys[a]) - a0))
                end
            end
            # random ordered (vertex, edge) triples; the distribution objects are stateless
            for _ in 1:n_triples
                f = K.uniform_int(rng, 0, F - 1) + 1
                pf = c.prime_faces[f]
                k0 = K.uniform_int(rng, 0, length(pf) - 1)
                A = pf[k0 + 1]; Bb = pf[mod1(k0 + 2, length(pf))]
                P = K.uniform_int(rng, 0, c.n_prime_vertices - 1) + 1
                (P == A || P == Bb) && continue
                y = Vector{Float64}(undef, n_angles)
                # geometric scale |Yb-Ya| * |Yp-Ya|: the largest |det| the triple can reach
                geo = 0.0; umax = 0.0; vmax = 0.0
                for a in 1:n_angles
                    u = Ys[a][Bb] - Ys[a][A]; v = Ys[a][P] - Ys[a][A]
                    y[a] = u[1] * v[2] - u[2] * v[1]
                    geo = max(geo, norm(u) * norm(v))
                    umax = max(umax, norm(u))
                    vmax = max(vmax, norm(v))
                end
                # numerically coincident triples are counted separately, never dropped silently
                if umax < 1e-6 * g_scale_diam || vmax < 1e-6 * g_scale_diam
                    g_degen += 1
                    continue
                end
                geo <= 0 && continue
                fit, mr, sc = K.harmonic_fit(ths, y)
                rel = sc > 0 ? mr / sc : mr
                rel_geo = mr / geo
                g_res = max(g_res, rel)
                g_res_geo = max(g_res_geo, rel_geo)
                g_tr += 1
                if rel >= 1e-10
                    n_fail += 1
                    g_deg = min(g_deg, sc / geo)
                end
                rel_geo >= 1e-10 && (n_fail_geo += 1)
                cf = K.orient_harmonic(B, A, Bb, P)
                cfe = 0.0
                for a in 1:n_angles
                    cfe = max(cfe, abs(K.harmonic_eval(cf, ths[a]) - y[a]))
                end
                g_cf = max(g_cf, cfe / geo)
                g_coef = max(g_coef, (abs(fit.p - cf.p) + abs(fit.q - cf.q) + abs(fit.r - cf.r)) / geo)
            end
        end
        print(csv, id, ",", kind, ",", K.n_vertices(m), ",", F, ",", sh.k, ",",
              g_tr, ",", cpp_g(g_res), ",", cpp_g(g_res_geo), ",", cpp_g(g_cf), ",", cpp_g(g_coef), ",",
              cpp_g(g_basis), ",", cpp_g(g_scale), ",", cpp_g(g_area), ",", cpp_g(g_ascale), ",", cpp_g(g_deg), ",",
              g_degen, ",", cpp_g(s(t)), "\n")
        flush(csv)
        graphs += 1
        n_triples_total += g_tr
        n_degen_total += g_degen
        minF = min(minF, F)
        maxF = max(maxF, F)
        if g_res > worst_res
            worst_res = g_res; worst_res_graph = name
        end
        if g_res_geo > worst_res_geo
            worst_res_geo = g_res_geo; worst_res_geo_graph = name
        end
        worst_cf = max(worst_cf, g_cf)
        g_deg < worst_degeneracy && (worst_degeneracy = g_deg)
        if g_area / max(1e-300, g_ascale) > worst_area
            worst_area = g_area / max(1e-300, g_ascale)
            worst_area_graph = name
        end
        worst_coef = max(worst_coef, g_coef)
        worst_basis = max(worst_basis, g_basis / max(1e-300, g_scale))
        if graphs % 20 == 0
            println("  ", graphs, " graphs, worst rel res ", cpp_g(worst_res), ", ", cpp_g(s(wall)), " s")
        end
    end

    for rc in reference_cases(regenerate = regenerate)
        run_one(rc.name, "reference", rc.mesh, -1)
    end
    nrun = 0
    for g in population("k1a_200"; regenerate = regenerate)
        g.id < n || continue
        nrun >= limit && break
        g.ok || continue
        run_one(g.kind * "_" * string(g.id), g.kind, g.mesh, g.id)
        nrun += 1
    end
    close(csv)

    sm = open(joinpath(outdir, "summary.txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("K1b: $graphs graphs, F in [$minF,$maxF], $n_triples_total triples\n")
    both("worst relative LS residual of p+q cos+r sin : " * sci(worst_res) * " (" * worst_res_graph * ")\n")
    both("triples with residual/max|det| >= 1e-10      : $n_fail  (smallest max|det|/(|AB||AP|) among them: " *
         sci(worst_degeneracy) * ")\n")
    both("worst residual / (|AB| |AP|)                : " * sci(worst_res_geo) * " (" * worst_res_geo_graph * ")\n")
    both("triples with residual/(|AB||AP|) >= 1e-10   : $n_fail_geo\n")
    both("worst |closed form - deploy| / (|AB| |AP|)  : " * sci(worst_cf) * "\n")
    both("triples skipped as numerically coincident   : $n_degen_total (both |AB| and |AP| < 1e-6 of the structure diameter)\n")
    both("worst |fitted - closed-form| coefficients   : " * sci(worst_coef) * "\n")
    both("worst relative |deploy - (cos C + sin S)|   : " * sci(worst_basis) * "\n")
    both("worst relative face signed-area drift       : " * sci(worst_area) * " (" * worst_area_graph * ")\n")
    both("VERDICT (spec rule, residual/max|det|)     : " * (n_fail == 0 ? "PASS" : "FAIL") * "\n")
    both("VERDICT (residual/(|AB||AP|), see note)    : " * (n_fail_geo == 0 ? "PASS" : "FAIL") * "\n")
    both("wall " * to_string_d(s(wall)) * " s\n")
    close(sm)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
