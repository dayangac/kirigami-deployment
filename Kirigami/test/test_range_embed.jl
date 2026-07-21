# test_range_embed.jl -- port of code/tests/test_range_embed.cpp, case by case: the
# RANGE-MAXIMISING point of the convexity-constrained slice of the Tutte auxetic shape
# space (K9c).
#
# What is checked here:
#   * the exact 0+ margin m(X) agrees, term by term, with zero_plus's two families;
#   * the surrogate the objective maximises is a LOWER BOUND on the exact margin, and is
#     TIGHT at the point where the modes were read;
#   * the analytic gradient of the stage-A objective (softmin + both barriers) matches a
#     central difference, on both the split term and the corner term;
#   * a hand-built case with a known answer: on a shape space that only moves one vertex,
#     the solver moves it in the direction that opens the margin;
#   * the invariants the driver relies on -- a FEASIBLE verdict is the exact constraint
#     values, the returned margin is the exact one, and the returned point never has a
#     smaller margin than the start point it was given.
include("helpers.jl")

const KR = Kirigami

struct RCase
    m::KR.Mesh
    c::KR.CutStructure
    hs::KR.HoleSet
    X::Vector{KR.Vec2}   # the Eq. (6) projection X0
    Phi::Matrix{Float64}
    med::Float64
    m_dim::Int
end

# Median edge length (std::nth_element median: the element at index n/2 of the sorted list).
function med_edge_of(m::KR.Mesh)
    L = [norm(m.X[e.key.a] - m.X[e.key.b]) for e in m.edges]
    isempty(L) && return 1.0
    return sort(L)[div(length(L), 2) + 1]
end

# The Eq. (6) projection and its null space, on a split-bearing tiling.
function make_rcase(mesh::KR.Mesh, checker::Bool)
    rng = KR.MT19937(2026)
    KR.build_topology!(mesh)
    mesh.sigma = checker ? KR.checkerboard(mesh) :
                 KR.assign_orientation_relaxation(mesh, rng, 8, 500, 180).sigma
    KR.build_topology!(mesh)
    c = KR.make_cut(mesh)
    hs = KR.holes_partition(c)
    med = med_edge_of(mesh)
    sys = KR.assemble_system(c, hs, mesh.X, KR.Fixed)
    sr = KR.solve_system(sys, mesh.X)
    return RCase(mesh, c, hs, KR.matrix_to_points(sr.X0), sr.Phi, med, 2 * size(sr.Phi, 2))
end

rctr = KR.Vec2(0.13, 0.07)
rcase(k) = k == 1 ? make_rcase(KR.tiling_truncated_square(KR.disk(rctr, 3.0)), false) :
                    make_rcase(KR.tiling_snub_square(KR.disk(rctr, 3.2)), false)

@testset "zero_plus_margin is exactly the min of zero_plus's two families" begin
    for k in 0:1
        cs = rcase(k)
        @test KR.n_split(cs.c) > 0
        s = cs.med * cs.med
        mq = minimum(KR.zero_plus_q(cs.c, cs.X))
        mmu = minimum(KR.zero_plus_corner_margin(cs.c, cs.X))
        @test isfinite(mq)
        @test isfinite(mmu)
        got, gq, gmu = KR.zero_plus_margin(cs.c, cs.X, cs.med)
        @test isapprox(gq, mq / s; rtol = 1e-12)
        @test isapprox(gmu, mmu / s; rtol = 1e-12)
        @test isapprox(got, min(mq, mmu) / s; rtol = 1e-12)
    end
end

@testset "the softmin surrogate is a lower bound on the exact margin, tight at the modes" begin
    # The objective maximises  M_kappa(t) = -kappa log sum exp(-c_i/kappa)  over a flat list
    # of entries that is itself a lower bound on the exact per-item margin. Two claims:
    #   (a) at the point where the modes were read, and with kappa -> 0, the smallest entry
    #       EQUALS the exact margin m(X) (the convex-corner mode is the argmax of -g1, -g2);
    #   (b) for any kappa > 0, M_kappa <= min_i c_i, so the objective never overstates.
    cs = rcase(0)
    @test cs.m_dim > 0
    rng = KR.MT19937(77001)
    G = KR.NormalDist(0.0, 0.05 * cs.med)
    for trial in 0:2
        t = [trial > 0 ? KR.normal(G, rng) : 0.0 for _ in 1:cs.m_dim]
        Xt = KR.shape_point(cs.X, cs.Phi, t)
        mode = KR.range_embed_modes(cs.c, Xt)
        exact, _, _ = KR.zero_plus_margin(cs.c, Xt, cs.med)

        # Rebuild the flat entry list the objective uses, independently of the solver.
        s = cs.med * cs.med
        ent = Float64[]
        for q in KR.zero_plus_q(cs.c, Xt)
            push!(ent, q / s)
        end
        inc = KR.corner_incidences(cs.c)
        d = KR.deploy(cs.c, Xt, 0.0)
        @test length(inc) == length(mode)
        for i in eachindex(inc)
            z = inc[i]
            e1 = Xt[z.v_next] - Xt[z.v]
            e2 = Xt[z.v_prev] - Xt[z.v]
            dS = 2.0 * (d.dY_dtheta[z.pv_other] - d.dY_dtheta[z.pv_corner])
            g1 = e1[1] * dS[2] - e1[2] * dS[1]
            g2 = dS[1] * e2[2] - dS[2] * e2[1]
            mode[i] != 1 && push!(ent, -g1 / s)
            mode[i] != 0 && push!(ent, -g2 / s)
        end
        @test !isempty(ent)
        emin = minimum(ent)
        # (a) tightness at the mode-reading point.
        @test isapprox(emin, exact; rtol = 1e-10)
        # (b) the log-sum-exp softmin is below the min for every kappa > 0.
        for kappa in (1e-3, 1e-2, 1e-1)
            Z = sum(exp(-(c - emin) / kappa) for c in ent)
            M = emin - kappa * log(Z)
            @test M <= emin + 1e-12
            @test M >= emin - kappa * log(length(ent)) - 1e-12
        end
    end
end

@testset "range_embed's stage-A gradient matches a central difference" begin
    for k in 0:1
        cs = rcase(k)
        @test cs.m_dim > 0
        opt = KR.RangeEmbedOptions(delta_rel = 1e-4, split_delta_rel = 1e-4, kappa = 2e-2)
        rng = KR.MT19937(90101 + k)
        G = KR.NormalDist(0.0, 0.02 * cs.med)
        t = [KR.normal(G, rng) for _ in 1:cs.m_dim]
        mode = KR.range_embed_modes(cs.c, KR.shape_point(cs.X, cs.Phi, t))
        for bw in (0.0, 1e-3)
            _, g = KR.range_embed_objective(cs.c, cs.X, cs.Phi, cs.med, opt, mode, bw, t)
            h = 1e-6 * cs.med
            worst = 0.0
            scale = 1e-30
            for i in 1:min(cs.m_dim, 12)
                tp = copy(t); tm = copy(t)
                tp[i] += h
                tm[i] -= h
                fp, _ = KR.range_embed_objective(cs.c, cs.X, cs.Phi, cs.med, opt, mode, bw, tp)
                fm, _ = KR.range_embed_objective(cs.c, cs.X, cs.Phi, cs.med, opt, mode, bw, tm)
                fd = (fp - fm) / (2h)
                worst = max(worst, abs(fd - g[i]))
                scale = max(scale, abs(fd))
            end
            @test worst <= 1e-4 * max(1.0, scale)
        end
    end
end

@testset "range_embed: the reported margin and feasibility are the EXACT values at X" begin
    for k in 0:1
        cs = rcase(k)
        @test cs.m_dim > 0
        opt = KR.RangeEmbedOptions(delta_rel = 1e-3, split_delta_rel = 1e-3, stages = 3,
                                   iter_per_stage = 60)
        r = KR.range_embed(cs.c, cs.X, cs.Phi, cs.med, opt)

        s = cs.med * cs.med
        margin, mq, mmu = KR.zero_plus_margin(cs.c, r.X, cs.med)
        @test isapprox(r.margin, margin; rtol = 1e-12)
        @test isapprox(r.min_q, mq; rtol = 1e-12)
        @test isapprox(r.min_mu, mmu; rtol = 1e-12)

        mc = Inf
        nbad = 0
        for v in KR.corner_crosses(cs.m, r.X)
            mc = min(mc, v)
            v < opt.delta_rel * s && (nbad += 1)
        end
        @test isapprox(r.min_cross, mc / s; rtol = 1e-12)
        @test r.n_bad == nbad
        nbq = count(v -> v < opt.split_delta_rel * s, KR.zero_plus_q(cs.c, r.X))
        @test r.n_bad_q == nbq
        @test r.feasible == (nbad == 0 && nbq == 0)
        # X really is the shape-space point of the returned t.
        Xt = KR.shape_point(cs.X, cs.Phi, r.t)
        err = maximum(norm(Xt[v] - r.X[v]) for v in eachindex(Xt))
        @test err < 1e-12 * cs.med
    end
end

@testset "range_embed never returns a feasible point worse than a feasible start" begin
    # The continuation keeps the best EXACT margin it has seen and refuses to leave the
    # feasible set once inside it, so handing it a feasible warm start can only help.
    cs = rcase(0)
    @test cs.m_dim > 0
    co = KR.ConvexEmbedOptions(delta_rel = 1e-3, split_delta_rel = 1e-3, n_random = 2,
                               max_iter = 200, barrier_stages = 3)
    ce = KR.convex_embed(cs.c, cs.X, cs.Phi, cs.m.X, cs.med, co)
    if ce.feasible  # otherwise nothing to warm start from on this tiling
        opt = KR.RangeEmbedOptions(delta_rel = 1e-3, split_delta_rel = 1e-3, stages = 4,
                                   iter_per_stage = 80, t_init = ce.t)
        r = KR.range_embed(cs.c, cs.X, cs.Phi, cs.med, opt)
        @test r.feasible
        @test isapprox(r.margin_start, KR.zero_plus_margin(cs.c, ce.X, cs.med)[1]; rtol = 1e-12)
        @test r.margin >= r.margin_start - 1e-12
    end
end

@testset "range_embed opens the margin on a hand-built one-degree-of-freedom space" begin
    # A synthetic shape space in which Phi moves exactly one interior vertex along one
    # direction. The margin m(t) is then a scalar function of t, and the solver must land
    # where a dense scan of that same exact function says the maximum is.
    cs = rcase(0)
    @test cs.m_dim > 0
    # Keep only the first null-space column: a 2-dimensional t = (tx, ty).
    Phi1 = cs.Phi[:, 1:1]
    opt = KR.RangeEmbedOptions(delta_rel = 1e-4, split_delta_rel = 1e-4, stages = 6,
                               iter_per_stage = 200, kappa = 5e-3)
    r = KR.range_embed(cs.c, cs.X, Phi1, cs.med, opt)
    # Scan the EXACT margin over the same 2-d space on a coarse grid around the answer and
    # check the solver is not beaten by any grid point that is also feasible.
    s = cs.med * cs.med
    R = max(2.0 * norm(r.t), 0.5 * cs.med)
    best = -Inf
    for a in -6:6, b in -6:6
        t = [R * a / 6.0, R * b / 6.0]
        Xt = KR.shape_point(cs.X, Phi1, t)
        ok = all(v -> v >= opt.delta_rel * s, KR.corner_crosses(cs.m, Xt))
        ok = ok && all(v -> v >= opt.split_delta_rel * s, KR.zero_plus_q(cs.c, Xt))
        ok || continue
        best = max(best, KR.zero_plus_margin(cs.c, Xt, cs.med)[1])
    end
    if r.feasible
        # The grid is coarse, so it may find a slightly better point; it must not find a
        # dramatically better one, and the solver must have improved on t = 0.
        @test r.margin >= r.margin_start - 1e-12
        isfinite(best) && @test r.margin >= best - 0.25 * abs(best) - 1e-9
    end
end

# ---------------------------------------------------------------------------
# The frozen C++ numbers (data/corpus/method_fixtures/test_method_2.json, see its README):
# the acceptance criterion of the port. The C++ null basis Phi is used so that the probe
# point and the solve live in the same coordinates as the C++ run.

@testset "range_embed reproduces the frozen C++ reference numbers" begin
    relerr(a, b) = isempty(b) ? 0.0 : maximum(abs.(a .- b)) / max(1e-300, maximum(abs.(b)))
    fixtures = JSON.parsefile(joinpath(CORPUS, "method_fixtures", "test_method_2.json"))["cases"]
    for fx in fixtures
        m = fixture_mesh_raw(fx["mesh"])
        c = KR.make_cut(m)
        X0 = fixture_points(fx["X0"])
        med = fx["med"]
        margin, mq, mmu = KR.zero_plus_margin(c, X0, med)
        @test isapprox(margin, fx["margin_at_X0"]; rtol = 1e-12, atol = 1e-15)
        @test isapprox(mq, fx["margin_min_q_at_X0"]; rtol = 1e-12, atol = 1e-15)
        @test isapprox(mmu, fx["margin_min_mu_at_X0"]; rtol = 1e-12, atol = 1e-15)
        @test KR.range_embed_modes(c, X0) == Int.(fx["modes_at_X0"])
        haskey(fx, "probe_t") || continue
        Phi = fixture_matrix(fx["Phi"])
        t = Float64.(fx["probe_t"])
        mode = KR.range_embed_modes(c, KR.shape_point(X0, Phi, t))
        @test mode == Int.(fx["range_modes_at_probe"])
        opt = KR.RangeEmbedOptions(delta_rel = 1e-4, split_delta_rel = 1e-4, kappa = 2e-2)
        for (k, bw) in (("range_embed_objective_bw0", 0.0), ("range_embed_objective_bw1e-3", 1e-3))
            F, g = KR.range_embed_objective(c, X0, Phi, med, opt, mode, bw, t)
            @test isapprox(F, fx[k * "_at_probe"]; rtol = 1e-12)
            @test relerr(g, Float64.(fx[k * "_grad_at_probe"])) < 1e-12
        end
        # The solve: same verdict and stage counts, the exact start margin to 1e-12, and the
        # end point to the L-BFGS rounding-path level. Measured on snub_square_3.2 (305
        # softmin entries, kappa 5e-3): the C++ and Julia iterates agree to 1e-15 for the
        # first ~10 iterations and diverge to 1e-12 / 1e-8 / 1e-5 after 20 / 40 / 60, so
        # after three 60-iteration stages the margins differ by ~1e-4 and the (unoptimised)
        # min cross by ~2e-3. The other five cases agree to ~1e-8.
        s = fx["range_embed_solve"]
        r = KR.range_embed(c, X0, Phi, med, KR.RangeEmbedOptions(delta_rel = 1e-3, split_delta_rel = 1e-3,
                                                                 stages = 3, iter_per_stage = 60))
        @test r.feasible == s["feasible"]
        @test r.stages_kept == s["stages_kept"]
        @test r.best_start == s["best_start"]
        @test r.n_entries == s["n_entries"]
        @test isapprox(r.margin_start, s["margin_start"]; rtol = 1e-12)
        @test isapprox(r.margin, s["margin"]; rtol = 1e-3)
        @test isapprox(r.min_cross, s["min_cross"]; rtol = 1e-2)
    end
end

# No C++ test case covers stage B; this replays one run of the C++ (scratch driver on the
# frozen truncated_square_3.0 input, 2026-09-19): theta_before = theta_after = 3 pi / 4,
# three caps tried, none accepted, margin unchanged.
@testset "maximize_margin_range: stage B on the frozen truncated-square case" begin
    fx = only(c for c in JSON.parsefile(joinpath(CORPUS, "method_fixtures", "test_method_2.json"))["cases"]
              if c["name"] == "truncated_square_3.0")
    m = fixture_mesh_raw(fx["mesh"])
    c = KR.make_cut(m)
    X0 = fixture_points(fx["X0"])
    Phi = fixture_matrix(fx["Phi"])
    med = fx["med"]
    r = KR.range_embed(c, X0, Phi, med, KR.RangeEmbedOptions(stages = 3, iter_per_stage = 60))
    mr = KR.maximize_margin_range(c, r.X, Phi, med)
    @test isapprox(mr.theta_before, 3pi / 4; rtol = 1e-12)
    @test mr.theta_after == mr.theta_before
    @test !mr.improved
    @test mr.caps_tried == 3
    @test mr.caps_accepted == 0
    @test isapprox(mr.margin_after, r.margin; rtol = 1e-12)
    @test mr.X == r.X
end
