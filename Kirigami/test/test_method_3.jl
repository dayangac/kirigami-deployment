# test_method_3.jl -- port of the tests/test_method.cpp cases that exercise
# method/periodic_jacobian, method/budget and method/expansive_cone, case by case.
#
# Inputs are built the way the C++ builds them (generators + MT19937 + orientation /
# quotient_sigma). data/corpus/method_fixtures/method_3_reference.json holds the C++
# numbers of the same calls (produced by the scratchpad freezer described in its
# "provenance" field); the last testset compares against them, so a mismatch there
# separates "the port disagrees with the C++" from "the self-consistency identity fails".
include("helpers.jl")
import JSON

const K = Kirigami
const M3_REF = JSON.parsefile(joinpath(CORPUS, "method_fixtures", "method_3_reference.json"))

# ---------------------------------------------------------------- helpers (C++ test fixtures)

# tests/helpers.hpp checkerboard_sigma: BFS 2-colouring of the dual graph.
function m3_checkerboard_sigma(m::K.Mesh)
    adj = K.dual_graph(m)
    sig = zeros(Int, K.n_faces(m))
    for s in 1:K.n_faces(m)
        sig[s] != 0 && continue
        sig[s] = -1
        st = Int[s]
        while !isempty(st)
            f = pop!(st)
            for g in adj[f]
                if sig[g] == 0
                    sig[g] = -sig[f]
                    push!(st, g)
                end
            end
        end
    end
    return sig
end

struct M3Case
    m::K.Mesh
    c::K.CutStructure
    hs::K.HoleSet
    X::Vector{K.Vec2}
end

# A deployable, embedded configuration: X_ini if it already satisfies Eq. (2),
# otherwise the Eq. (6) projection.
function m3_make_case(m::K.Mesh, checker::Bool)
    rng = K.MT19937(2026)
    K.build_topology!(m)
    m.sigma = checker ? m3_checkerboard_sigma(m) :
              K.assign_orientation_relaxation(m, rng, 8, 500, 180).sigma
    K.build_topology!(m)
    c = K.make_cut(m)
    hs = K.holes_partition(c)
    X = copy(m.X)
    if !K.deployable(K.hole_residuals(c, X, hs), 1e-9)
        sys = K.assemble_system(c, hs, m.X, K.Fixed)
        sr = K.solve_system(sys, m.X)
        X = K.matrix_to_points(sr.X0)
    end
    return M3Case(m, c, hs, X)
end

m3_disk(r) = K.disk(K.Vec2(0.13, 0.07), r)
m3_case_hexagons() = m3_make_case(K.tiling_hexagons(m3_disk(2.0)), false)
m3_case_trunc() = m3_make_case(K.tiling_truncated_square(m3_disk(2.5)), false)
m3_case_snub() = m3_make_case(K.tiling_snub_square(m3_disk(2.0)), false)
m3_case_t34312() = m3_make_case(K.tiling_3_4_3_12(m3_disk(2.5)), false)
m3_case_squares() = m3_make_case(K.tiling_squares(K.rect(K.Vec2(2.5, 2.5), 2.01, 2.01)), true)

# The closed-form second fully-closed angle of a periodic pattern (K7 / C4).
struct M3PeriodicCase
    q::K.Quotient
    sp::K.SuperPatch
    cut::K.CutStructure
    J::K.PeriodicJac
    face_area::Float64
end

function m3_make_periodic_case(family, n, m, seed)
    rng = K.MT19937(seed)
    P = K.make_tiling_pattern(family, n, m, rng)
    @test P.ok
    P.ok || error("make_tiling_pattern: " * P.err)
    q = K.build_quotient(P.cell, P.T)
    @test q.ok
    q.ok || error("build_quotient: " * q.err)
    sys = K.quotient_system(q, q.Xq)
    sr = K.solve_system(sys, q.Xq)
    @test sr.projection_ok
    sp = K.build_super(q, 1)
    K.set_super_positions!(sp, q, K.matrix_to_points(sr.X0))
    cut = K.make_cut(sp.mesh)
    K.set_super_positions!(sp, q, K.matrix_to_points(sr.X0))
    face_area = K.cell_face_area_sum(sp)
    J = K.periodic_jacobian(sp, q, cut, K.deploy_basis(cut, sp.mesh.X))
    @test J.ok
    return M3PeriodicCase(q, sp, cut, J, face_area)
end

# Signed per-cell hole area at theta, from an independent forward-kinematics deploy.
function m3_fk_hole_area(C::M3PeriodicCase, theta)
    P, _ = K.fk_period_matrix(C.sp, C.q, C.cut, Float64(theta))
    s = det(C.J.P0) >= 0 ? 1.0 : -1.0
    return s * det(P) - C.face_area
end

# A deployable fixed-boundary design on one of the generator families, with the
# relaxation sigma (so split cuts exist).
struct M3BudgetCase
    m::K.Mesh
    c::K.CutStructure
    X::Vector{K.Vec2}
    ok::Bool
end

function m3_make_budget_case(kind, par, seed)
    rng = K.MT19937(seed)
    m = K.generate(kind, [par], rng)
    K.build_topology!(m)
    m.sigma = K.assign_orientation_relaxation(m, rng).sigma
    K.build_topology!(m)
    c = K.make_cut(m)
    hs = K.holes_partition(c)
    sys = K.assemble_system(c, hs, m.X, K.Fixed)
    sr = K.solve_system(sys, m.X)
    sr.projection_ok || return M3BudgetCase(m, c, K.Vec2[], false)
    X = K.matrix_to_points(sr.X0)
    return M3BudgetCase(m, c, X, K.deployable(K.hole_residuals(c, X, hs), 1e-7))
end

const M3_BUDGET_KINDS = ["squares", "triangles", "hexagons", "kagome", "snub_square"]

# ---------------------------------------------------------------- periodic_jacobian

@testset "C4: theta_c = 2 atan2(tr K, 1 - det K) is where forward kinematics recloses" begin
    # (family, n, m, seed) reproduce the K7 population rows exactly: the seed is
    # 20260904 + 7919 * (index in specs/experimenter_k7.md's family x size grid).
    rows = [("hexagons", 2, 2, 20260904 + 7919 * 6, 2.09439510239),
            ("t3_4_3_12", 2, 2, 20260904 + 7919 * 18, 2.72771226181)]
    for (fam, n, m, seed, expect) in rows
        C = m3_make_periodic_case(fam, n, m, seed)
        Km = C.J.K
        dP0 = abs(det(C.J.P0))
        q = -dP0 * (det(Km) - 1) / 2
        rr = dP0 * tr(Km) / 2

        # the two equivalent forms of the closed form agree
        thc = 2.0 * atan(rr, q)
        thc <= 0 && (thc += 2pi)
        thc_K = 2.0 * atan(tr(Km), 1.0 - det(Km))
        thc_K <= 0 && (thc_K += 2pi)
        @test isapprox(thc, thc_K; atol = 1e-12)
        @test isapprox(thc, expect; rtol = 1e-9)

        # the harmonic itself is exact against forward kinematics on (0, theta_c]
        fit = 0.0
        for i in 1:40
            th = thc * i / 40.0
            fit = max(fit, abs(m3_fk_hole_area(C, th) - (q * (cos(th) - 1) + rr * sin(th))))
        end
        @test fit < 1e-9 * (abs(q) + abs(rr) + dP0)

        # the hole area is strictly positive inside (0, theta_c) and closes AT theta_c
        @test m3_fk_hole_area(C, 0.5 * thc) > 1e-6 * dP0
        @test abs(m3_fk_hole_area(C, thc)) < 1e-9 * dP0

        # and the sign change bracketed by bisection lands on theta_c to 1e-6 (the K7 bar)
        lo = 0.5 * thc
        hi = min(2pi, thc + 0.5)
        @test m3_fk_hole_area(C, lo) > 0
        @test m3_fk_hole_area(C, hi) < 0
        for _ in 1:80
            mid = 0.5 * (lo + hi)
            if m3_fk_hole_area(C, mid) > 0
                lo = mid
            else
                hi = mid
            end
        end
        @test abs(0.5 * (lo + hi) - thc) < 1e-6
    end
end

# ---------------------------------------------------------------- expansive_cone

@testset "cone_lp solves three hand-solved linear programs" begin
    # The value solved for is  max_{||z||_2 <= 1} min_i a_i.z  =  min_{l in simplex}
    # ||A^T l||_2, and the rows below are already unit norm.
    @testset "infeasible: z1 > 0 and -z1 > 0" begin
        A = reshape([1.0, -1.0], 2, 1)
        r = K.cone_lp(A)
        @test r.feasible == false
        @test r.margin_l2 <= 1e-8
        # Farkas: l = (1/2, 1/2) gives A^T l = 0 exactly, so the dual bound must reach 0.
        @test r.dual_bound < 1e-6
        @test abs(r.lambda[1] - 0.5) < 1e-3
    end
    @testset "orthogonal rows: value = 1/sqrt(2)" begin
        A = Matrix{Float64}(I, 2, 2)
        r = K.cone_lp(A)
        @test r.feasible
        @test isapprox(r.margin_l2, 1 / sqrt(2); rtol = 1e-3)
        @test isapprox(r.dual_bound, 1 / sqrt(2); rtol = 1e-3)
        # ||z||_inf = 1/sqrt(2) at the optimum, so the l_inf margin is 1.
        @test isapprox(r.margin_inf, 1.0; rtol = 1e-3)
    end
    @testset "z1 > 0 and z2 - z1 > 0: value = 1/sqrt(5)" begin
        # Rows (1,0) and (-1,1)/sqrt(2); the optimum is bracketed by the two bounds.
        A = [1.0 0.0; -1 / sqrt(2) 1 / sqrt(2)]
        r = K.cone_lp(A)
        @test r.feasible
        @test r.margin_l2 > 0.3
        @test r.dual_bound >= r.margin_l2 - 1e-9       # the bracket is two-sided
        @test r.dual_bound - r.margin_l2 < 1e-2        # and it is tight
    end
end

@testset "the cone rows reproduce zero_plus's q_e and corner margins at the sigma ray" begin
    for (k, cs) in enumerate([m3_case_hexagons(), m3_case_trunc(), m3_case_snub()])
        sys = K.cone_system(cs.c, cs.X)
        u = K.sigma_flex(cs.c, cs.X)
        val = sys.A * u
        q = K.zero_plus_q(cs.c, cs.X)
        mu = K.zero_plus_corner_margin(cs.c, cs.X)
        @test sys.n_split == length(q)
        scale = 1e-12
        for v in q
            scale = max(scale, abs(v))
        end
        for v in mu
            scale = max(scale, abs(v))
        end

        # Split rows: row_e = q_e / 2 exactly (dS_e = 2 dV).
        for i in 1:sys.n_split
            @test abs(val[i] - 0.5 * q[i]) < 1e-9 * scale
        end

        # Corner rows: mu = max/min(-g1, -g2) with the branch read off the corner.
        for r in eachindex(sys.rows)
            sys.rows[r].kind != K.ConeCornerG1 && continue
            idx = sys.rows[r].index
            m1 = val[r]
            m2 = val[r + 1]
            got = sys.rows[r].convex ? max(m1, m2) : min(m1, m2)
            @test abs(got - 0.5 * mu[idx]) < 1e-9 * scale
        end
    end
end

@testset "flex_basis: dimension identity, residual, and the sigma ray inside it" begin
    for (k, cs) in enumerate([m3_case_squares(), m3_case_hexagons(), m3_case_t34312()])
        g = K.build_hinge_graph(cs.c)
        pins = K.pins_flat(cs.c, g, cs.X)
        fb = K.flex_basis(g, pins)
        # dim ker R = dim ker A + 2 c(Gamma), and N really is a flex.
        @test K.dim(fb) == fb.dim_ker_A + 2 * g.components
        @test fb.residual < 1e-9
        @test size(fb.N, 2) > 0
        # Orthonormal columns.
        G = fb.N' * fb.N
        @test maximum(abs, G - I) < 1e-10
        # Against an independent dense rank of the rigidity matrix (the C++ uses
        # ColPivHouseholderQR with threshold 1e-10 on R / max|R|).
        R = K.build_rigidity(g, pins)
        Rs = R / maximum(abs, R)
        qrR = qr(Rs, ColumnNorm())
        dR = abs.(diag(qrR.R))
        rk = count(x -> x > 1e-10 * maximum(dR), dR)
        @test K.dim(fb) == 3 * g.F - rk
        # The uniform ray is a flex of a deployable embedding, and lies in span(N).
        u = K.sigma_flex(cs.c, cs.X)
        resid = u - fb.N * (fb.N' * u)
        @test norm(resid) / norm(u) < 1e-8
        @test norm(R * u) / (maximum(abs, R) * norm(u)) < 1e-9
    end
end

@testset "expansive_cone: sigma is inside P(X) on the split-bearing tilings" begin
    for cs in [m3_case_trunc(), m3_case_snub()]
        opt = K.ExpansiveConeOptions()
        opt.lp.mu_min = 1e-4
        rep = K.expansive_cone(cs.c, cs.X, opt)
        @test rep.ok
        @test rep.sigma_is_flex
        @test rep.sigma_n_bad_q == 0
        @test rep.sigma_n_bad_mu == 0
        @test rep.sigma_in_cone
        @test rep.flex_residual < 1e-9
    end
end

@testset "euler_step along a feasible flex leaves the structure collision free" begin
    cs = m3_case_trunc()
    u = K.sigma_flex(cs.c, cs.X)
    sp = K.flex_speed(cs.c, cs.X, u)
    @test sp > 0
    L = [norm(cs.m.X[e.key.a] - cs.m.X[e.key.b]) for e in cs.m.edges]
    med = sort(L)[length(L) ÷ 2 + 1]   # nth_element at index size/2 (0-based)
    h = 0.01 * med / sp
    Y = K.euler_step(cs.c, cs.X, u, h)
    @test !K.has_collision(cs.c, Y)
    # The Euler step along the uniform ray agrees with the exact kinematics to O(h^2).
    Yex = K.deploy(cs.c, cs.X, h).Y
    dmax = maximum(norm(Y[i] - Yex[i]) for i in eachindex(Y))
    @test dmax < 1e-3 * med
end

# --- K8a RECHECK: the corrected min-norm-point solver ------------------------------

@testset "cone_lp closes the bracket: two unit rows at a known angle" begin
    # For two UNIT rows at angle phi (phi < pi) the min-norm point of the segment is the
    # bisector scaled by cos(phi/2), so val = cos(phi/2) exactly.
    for phi in (0.4, 1.0, 2pi / 3, 3.0)
        A = [1.0 0.0; cos(phi) sin(phi)]
        r = K.cone_lp(A)
        val = cos(0.5 * phi)
        @test isapprox(r.margin_l2, val; rtol = 1e-9)
        @test isapprox(r.dual_bound, val; rtol = 1e-9)
        @test r.gap < 1e-9
        @test r.feasible
    end
end

@testset "cone_lp finds a small positive margin in a many-row, high-dimension system" begin
    # 1200 rows in 60 dimensions, all inside a narrow cone about e_1: strictly feasible
    # but only by a small margin; z = e_1 is an explicit witness. Deterministic seed and
    # the libc++ uniform_real stream (mt19937.jl), so A is the C++ matrix.
    m, n = 1200, 60
    rng = K.MT19937(20260904)
    A = zeros(m, n)
    for i in 1:m
        row = zeros(n)
        for j in 2:n
            row[j] = K.uniform_real(rng, -1.0, 1.0)
        end
        row ./= norm(row)
        A[i, :] = 0.30 .* [1.0; zeros(n - 1)] .+ 0.05 .* row
    end
    K.normalise_rows!(A)
    witness = minimum(A[:, 1])   # margin of z = e_1, a feasible unit vector
    @test witness > 0
    r = K.cone_lp(A)
    @test r.feasible
    @test r.margin_l2 >= witness - 1e-9          # never worse than the explicit witness
    @test r.dual_bound >= r.margin_l2 - 1e-9     # the bracket is two-sided ...
    @test r.gap < 1e-6                           # ... and it closes
    ref = M3_REF["cone_lp_1200x60"]
    @test isapprox(A[1, 1:5], Float64.(ref["A_first_row_first5"]); rtol = 1e-12)
    @test isapprox(witness, ref["witness"]; rtol = 1e-12)
end

@testset "farkas_residual verifies a certificate independently of the solver" begin
    A = reshape([1.0, -1.0], 2, 1)
    res, lmin, serr = K.farkas_residual(A, [0.5, 0.5])
    @test res < 1e-15
    @test isapprox(lmin, 0.5)
    @test serr < 1e-15
    res, lmin, serr = K.farkas_residual(A, [1.0, 0.0])
    @test isapprox(res, 1.0)
    @test isapprox(lmin, 0.0; atol = 1e-15)
    # A lambda off the simplex must be caught by sum_err, not silently accepted.
    res, lmin, serr = K.farkas_residual(A, [0.25, 0.25])
    @test res < 1e-15
    @test isapprox(serr, 0.5)
    # And the solver's own certificate on this system survives the independent recompute.
    r = K.cone_lp(A)
    res, lmin, serr = K.farkas_residual(A, r.lambda)
    @test res < 1e-9
    @test lmin >= -1e-12
    @test serr < 1e-9
end

@testset "expansive_cone: the LP beats the sigma ray whenever sigma is inside P(X)" begin
    # sigma_chart_margin is computed with no solver in the loop, so it is a rigorous lower
    # bound on that chart's LP value; a correct solver can never return less.
    for cs in [m3_case_hexagons(), m3_case_trunc(), m3_case_snub(), m3_case_t34312()]
        rep = K.expansive_cone(cs.c, cs.X)
        @test rep.ok
        @test rep.sigma_in_cone
        @test rep.sigma_in_flex_span
        @test rep.sigma_chart_margin > 0
        @test rep.margin_l2 >= rep.sigma_chart_margin - 1e-9
        @test rep.dual_bound >= rep.margin_l2 - 1e-9
        @test rep.lp_gap_max < 1e-6
        @test rep.pass_feasible >= 0
        @test rep.farkas_lambda_min >= -1e-12
        @test rep.farkas_sum_err < 1e-9
    end
end

# ---------------------------------------------------------------- budget (B1/B3)

@testset "face_potential closes on every hinge edge (T1 potential exists)" begin
    for k in M3_BUDGET_KINDS
        bc = m3_make_budget_case(k, 3.5, 20260904)
        @test bc.ok
        bc.ok || continue
        u, clo = K.face_potential(bc.c, bc.X)
        @test clo >= 0.0
        @test clo < 1e-10
        @test length(u) == K.n_faces(bc.m)
    end
end

@testset "the split budget term is exactly half the 0+ margin: q_e = 2 r_e" begin
    n = 0
    for k in M3_BUDGET_KINDS
        bc = m3_make_budget_case(k, 3.5, 20260904)
        @test bc.ok
        bc.ok || continue
        u, _ = K.face_potential(bc.c, bc.X)
        bt = K.budget_terms(bc.c, bc.X, u)
        q = K.zero_plus_q(bc.c, bc.X)
        sc = K.split_copies(bc.c)
        @test length(q) == length(sc)
        # split_copies drops non-manifold split edges; match by edge index.
        for i in eachindex(sc)
            pos = findfirst(==(sc[i].edge), bc.c.split_edges)
            @test pos !== nothing
            r = bt.split[pos]
            @test isapprox(q[i], 2.0 * r; rtol = 1e-11)
            n += 1
        end
    end
    @test n > 0
end

@testset "budget identity on a patch: W + R = 2 B(X) over the border alone" begin
    for k in M3_BUDGET_KINDS
        bc = m3_make_budget_case(k, 3.5, 20260904)
        @test bc.ok
        bc.ok || continue
        u, _ = K.face_potential(bc.c, bc.X)
        bt = K.budget_terms(bc.c, bc.X, u)
        B = K.border_functional(bc.c, bc.X, u)
        sa = abs(bt.W) + abs(bt.R) + abs(2 * B.a)
        sb = abs(bt.bW) + abs(bt.bR) + abs(2 * B.b)
        @test abs(bt.W + bt.R - 2 * B.a) < 1e-10 * max(1.0, sa)
        @test abs(bt.bW + bt.bR - 2 * B.b) < 1e-10 * max(1.0, sb)
    end
end

@testset "the measured void area is the border functional's first harmonic" begin
    th = [0.05, 0.2, 0.4, 0.7, 1.0, 1.3]
    for k in M3_BUDGET_KINDS
        bc = m3_make_budget_case(k, 3.5, 20260904)
        @test bc.ok
        bc.ok || continue
        u, _ = K.face_potential(bc.c, bc.X)
        B = K.border_functional(bc.c, bc.X, u)
        for t in th
            meas = K.measured_void_area(bc.c, K.deploy(bc.c, bc.X, t).Y)
            pred = B.a * sin(t) - B.b * (1 - cos(t))
            @test abs(meas - pred) < 1e-10 * max(1.0, abs(B.a) + abs(B.b))
        end
        @test abs(K.measured_void_area(bc.c, K.deploy(bc.c, bc.X, 0.0).Y)) < 1e-9
    end
end

@testset "periodic budget identity: W + R = det(P0) tr K over one cell" begin
    fams = ["triangles", "hexagons", "kagome", "trunc_square_488", "t3_4_3_12"]
    n = 0
    for f in fams
        rng = K.MT19937(20260904)
        P = K.make_tiling_pattern(f, 2, 2, rng)
        @test P.ok
        P.ok || continue
        q = K.build_quotient(P.cell, P.T)
        @test q.ok
        q.ok || continue
        sys = K.quotient_system(q, q.Xq)
        sr = K.solve_system(sys, q.Xq)
        @test sr.projection_ok
        sp = K.build_super(q, 1)
        c = K.make_cut(sp.mesh)
        K.set_super_positions!(sp, q, K.matrix_to_points(sr.X0))
        keep, npre = K.periodic_cell_edges(c, sp)
        nh = count(e -> keep[e], c.hinge_edges)
        ns = count(e -> keep[e], c.split_edges)
        @test nh == q.n_hinge
        @test ns == q.n_split
        @test npre == q.H
        u, clo = K.face_potential(c, sp.mesh.X)
        @test clo < 1e-10
        bt = K.budget_terms(c, sp.mesh.X, u, keep)
        B = K.deploy_basis(c, sp.mesh.X)
        pj = K.periodic_jacobian(sp, q, c, B)
        budget = det(pj.P0) * tr(pj.K)
        @test abs(bt.W + bt.R - budget) < 1e-10 * max(1.0, abs(bt.W) + abs(bt.R))
        # and the split half is exactly half the 0+ margins of the same cell
        qv = K.zero_plus_q(c, sp.mesh.X)
        sc = K.split_copies(c)
        sumq = 0.0
        for i in eachindex(sc)
            keep[sc[i].edge] && (sumq += qv[i])
        end
        @test abs(sumq - 2 * bt.R) < 1e-9 * max(1.0, abs(sumq))
        n += 1
    end
    @test n == 5
end

# ---------------------------------------------------------------- frozen C++ reference
#
# The same calls as above compared with the numbers the C++ produced. detect_lattice's
# candidate order comes from an unstable libc++ std::sort, so the lattice basis is compared
# up to the symmetries that leave the downstream cell equivalent (same |det|, same lengths);
# everything after it (sigma from quotient_sigma, K, theta_c, budgets) must agree.

@testset "frozen C++ reference: lattices, C4, periodic budget, patch budget, cone" begin
    bigs = Dict(
        "squares" => K.tiling_squares(K.rect(K.Vec2(0, 0), 6.5, 6.5)),
        "triangles" => K.tiling_triangles(K.disk(K.Vec2(0.05, 0.03), 7.0)),
        "hexagons" => K.tiling_hexagons(K.disk(K.Vec2(0.05, 0.03), 9.0)),
        "kagome" => K.tiling_kagome(K.disk(K.Vec2(0.05, 0.03), 7.0)),
        "snub_square" => K.tiling_snub_square(K.disk(K.Vec2(0.05, 0.03), 9.0)),
        "trunc_square_488" => K.tiling_truncated_square(K.disk(K.Vec2(0.05, 0.03), 13.0)),
        "t3_4_3_12" => K.tiling_3_4_3_12(K.disk(K.Vec2(0.05, 0.03), 15.0)))
    for (name, ref) in M3_REF["lattice"]
        big = bigs[name]
        @test K.n_vertices(big) == ref["n_vertices"]
        @test K.n_faces(big) == ref["n_faces"]
        lat = K.detect_lattice(big)
        @test lat.ok == ref["ok"]
        t1r = fixture_points([ref["t1"]])[1]
        t2r = fixture_points([ref["t2"]])[1]
        @test isapprox(norm(lat.t1), norm(t1r); rtol = 1e-12)
        @test isapprox(norm(lat.t2), norm(t2r); rtol = 1e-12)
        @test isapprox(K._det2(lat.t1, lat.t2), K._det2(t1r, t2r); rtol = 1e-12)
        @test K._det2(lat.t1, lat.t2) > 0
    end
    for ref in M3_REF["c4"]
        C = m3_make_periodic_case(ref["family"], ref["n"], ref["m"], ref["seed"])
        @test C.q.nq == ref["nq"]
        @test C.q.H == ref["H"]
        @test C.q.n_hinge == ref["n_hinge"]
        @test C.q.n_split == ref["n_split"]
        Kr = fixture_matrix(ref["K"])
        thc = 2.0 * atan(tr(C.J.K), 1.0 - det(C.J.K))
        thc <= 0 && (thc += 2pi)
        @test isapprox(thc, ref["theta_c"]; rtol = 1e-9)
        @test isapprox(tr(C.J.K), tr(Kr); rtol = 1e-9)
        @test isapprox(det(C.J.K), det(Kr); rtol = 1e-9)
        @test isapprox(abs(det(C.J.P0)), abs(det(fixture_matrix(ref["P0"]))); rtol = 1e-12)
        @test isapprox(C.face_area, ref["face_area"]; rtol = 1e-12)
        # sigma of the cell, when the lattice representative agreed with the C++
        cell_ref = fixture_mesh_raw(ref["cell"])
        if K.n_vertices(C.q.cell) == K.n_vertices(cell_ref) &&
           isapprox(C.q.cell.X, cell_ref.X; atol = 1e-12)
            @test C.q.cell.sigma == cell_ref.sigma
        else
            @info "C4 $(ref["family"]): lattice representative differs from the C++ (libc++ sort tie); cell compared through K only"
        end
    end
    for ref in M3_REF["periodic_budget"]
        rng = K.MT19937(20260904)
        P = K.make_tiling_pattern(ref["family"], 2, 2, rng)
        @test P.ok == ref["ok"]
        P.ok || continue
        @test K.n_faces(P.cell) == ref["n_cell_faces"]
        q = K.build_quotient(P.cell, P.T)
        @test q.n_hinge == ref["n_hinge"]
        @test q.n_split == ref["n_split"]
        @test q.H == ref["H"]
        sr = K.solve_system(K.quotient_system(q, q.Xq), q.Xq)
        sp = K.build_super(q, 1)
        c = K.make_cut(sp.mesh)
        K.set_super_positions!(sp, q, K.matrix_to_points(sr.X0))
        keep, npre = K.periodic_cell_edges(c, sp)
        @test npre == ref["n_preimage"]
        u, _ = K.face_potential(c, sp.mesh.X)
        bt = K.budget_terms(c, sp.mesh.X, u, keep)
        pj = K.periodic_jacobian(sp, q, c, K.deploy_basis(c, sp.mesh.X))
        @test isapprox(bt.W + bt.R, ref["W"] + ref["R"]; rtol = 1e-9)
        @test isapprox(det(pj.P0) * tr(pj.K), ref["det_P0_tr_K"]; rtol = 1e-9)
    end
    for ref in M3_REF["patch_budget"]
        bc = m3_make_budget_case(ref["kind"], 3.5, 20260904)
        @test bc.ok == ref["ok"]
        bc.ok || continue
        mref = fixture_mesh_raw(ref["mesh"])
        same_input = bc.m.sigma == mref.sigma && isapprox(bc.m.X, mref.X; atol = 1e-12)
        @test same_input
        # With the C++ inputs (frozen mesh + sigma + X) the budget numbers must agree.
        cr = K.make_cut(mref)
        Xr = fixture_points(ref["X"])
        u, clo = K.face_potential(cr, Xr)
        @test clo < 1e-10
        bt = K.budget_terms(cr, Xr, u)
        B = K.border_functional(cr, Xr, u)
        @test isapprox(bt.W, ref["W"]; rtol = 1e-10, atol = 1e-12)
        @test isapprox(bt.R, ref["R"]; rtol = 1e-10, atol = 1e-12)
        @test isapprox(bt.bW, ref["bW"]; rtol = 1e-10, atol = 1e-12)
        @test isapprox(bt.bR, ref["bR"]; rtol = 1e-10, atol = 1e-12)
        @test isapprox(B.a, ref["border_a"]; rtol = 1e-10, atol = 1e-12)
        @test isapprox(B.b, ref["border_b"]; rtol = 1e-10, atol = 1e-12)
        @test isapprox(K.measured_void_area(cr, K.deploy(cr, Xr, 0.7).Y), ref["void_area_0_7"];
                       rtol = 1e-10, atol = 1e-12)
    end
    for ref in M3_REF["cone"]
        # Run on the frozen C++ input (mesh + sigma + X), so the comparison isolates
        # flex_basis / cone_system / cone_lp from the generator and orientation ports.
        mref = fixture_mesh_raw(ref["mesh"])
        cr = K.make_cut(mref)
        Xr = fixture_points(ref["X"])
        g = K.build_hinge_graph(cr)
        pins = K.pins_flat(cr, g, Xr)
        fb = K.flex_basis(g, pins)
        @test K.dim(fb) == ref["dim_flex"]
        @test fb.dim_ker_A == ref["dim_ker_A"]
        @test g.components == ref["components"]
        sys = K.cone_system(cr, Xr)
        @test size(sys.A, 1) == ref["n_rows"]
        @test sys.n_split == ref["n_split"]
        @test sys.n_corner == ref["n_corner"]
        @test sys.n_convex == ref["n_convex"]
        rep = K.expansive_cone(cr, Xr)
        @test rep.pass_feasible == ref["pass_feasible"]
        @test rep.passes == ref["passes"]
        @test rep.sigma_in_cone == ref["sigma_in_cone"]
        @test isapprox(rep.sigma_chart_margin, ref["sigma_chart_margin"]; rtol = 1e-8)
        # nlohmann-json writes +inf (no split edges) as null.
        m3_inf(x) = x === nothing ? Inf : Float64(x)
        @test isapprox(rep.sigma_min_q, m3_inf(ref["sigma_min_q"]); rtol = 1e-8, atol = 1e-12)
        @test isapprox(rep.sigma_min_mu, m3_inf(ref["sigma_min_mu"]); rtol = 1e-8, atol = 1e-12)
        # C++ FRAGILITY (replicated, not fixed): rows of A N whose exact value is zero on
        # the flex space (a corner incidence whose two copies never separate) come out as
        # ~1e-16 rounding noise, and `normalise_rows` only skips rows of norm <= 0, so they
        # are scaled to unit rows of pure noise. The LP value then depends on the rounding
        # of the flex basis (Eigen Householder vs LAPACK): squares_checker (80 of 160 rows)
        # and truncated_square (16 of 55) are not reproducible and their margin is not
        # compared; the counts were measured on the C++ N itself. Where no such row exists
        # the LP value is a well-defined number and the closed brackets must agree.
        M0 = Matrix(sys.A * fb.N)
        n_noise = count(i -> norm(@view M0[i, :]) < 1e-12, 1:size(M0, 1))
        expected_noise = Dict("squares_checker" => 80, "truncated_square" => 16)
        @test n_noise == get(expected_noise, ref["name"], 0)
        if n_noise == 0
            @test isapprox(rep.margin_l2, ref["margin_l2"]; atol = 2e-6)
            @test isapprox(rep.dual_bound, ref["dual_bound"]; atol = 2e-6)
        end
        # And the generator + orientation ports reproduce the C++ input itself.
        cs = ref["name"] == "squares_checker" ? m3_case_squares() :
             ref["name"] == "hexagons" ? m3_case_hexagons() :
             ref["name"] == "truncated_square" ? m3_case_trunc() :
             ref["name"] == "snub_square" ? m3_case_snub() : m3_case_t34312()
        @test cs.m.sigma == mref.sigma
        @test isapprox(cs.X, Xr; atol = 1e-9)
    end
end
