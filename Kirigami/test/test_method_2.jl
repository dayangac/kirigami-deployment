# test_method_2.jl -- method/zero_plus and method/convex_embed, case by case (the
# related method cases live in test_method_1.jl and test_method_3.jl).
#
# The cases build their meshes with generators + assign_orientation_relaxation on
# MT19937(2026), which reproduces the frozen meshes.
include("helpers.jl")

const K2 = Kirigami

struct MCase
    m::K2.Mesh
    c::K2.CutStructure
    hs::K2.HoleSet
    X::Vector{K2.Vec2}
end

# A deployable, embedded configuration: X_ini if it already satisfies Eq. (2),
# otherwise the Eq. (6) projection.
function make_case(m::K2.Mesh, checker::Bool)
    rng = K2.MT19937(2026)
    K2.build_topology!(m)
    m.sigma = checker ? K2.checkerboard(m) : K2.assign_orientation_relaxation(m, rng, 8, 500, 180).sigma
    K2.build_topology!(m)
    c = K2.make_cut(m)
    hs = K2.holes_partition(c)
    X = copy(m.X)
    if !K2.deployable(K2.hole_residuals(c, X, hs), 1e-9)
        sys = K2.assemble_system(c, hs, m.X, K2.Fixed)
        sr = K2.solve_system(sys, m.X)
        X = K2.matrix_to_points(sr.X0)
    end
    return MCase(m, c, hs, X)
end

# The null-space basis Phi of the Eq. (6) system, as the K6 driver builds it.
null_basis(cs::MCase) =
    K2.solve_system(K2.assemble_system(cs.c, cs.hs, cs.m.X, K2.Fixed), cs.m.X).Phi

# Median edge length (std::nth_element median: the element at index n/2 of the sorted list).
function med_edge(m::K2.Mesh)
    L = [norm(m.X[e.key.a] - m.X[e.key.b]) for e in m.edges]
    isempty(L) && return 1.0
    return sort(L)[div(length(L), 2) + 1]
end

function n_inverted(m::K2.Mesh, X::Vector{K2.Vec2})
    t = deepcopy(m)
    t.X = X
    return count(f -> K2.face_signed_area(t, f) <= 0, 1:K2.n_faces(m))
end

# The natural scale of the corner margin mu = det(e, dS): |e| |dS|.
function corner_scale(c::K2.CutStructure, X::Vector{K2.Vec2}, z::K2.CornerIncidence)
    d = K2.deploy(c, X, 0.0)
    dS = 2.0 * (d.dY_dtheta[z.pv_other] - d.dY_dtheta[z.pv_corner])
    e = max(norm(X[z.v_next] - X[z.v]), norm(X[z.v_prev] - X[z.v]))
    return max(1e-300, e * norm(dS))
end

# Winding-number point-in-polygon, strict interior.
function strictly_inside(P::Vector{K2.Vec2}, p::K2.Vec2)
    wn = 0
    n = length(P)
    for i in 1:n
        a = P[i]; b = P[mod1(i + 1, n)]
        d = (b[1] - a[1]) * (p[2] - a[2]) - (b[2] - a[2]) * (p[1] - a[1])
        if a[2] <= p[2]
            (b[2] > p[2] && d > 0) && (wn += 1)
        elseif b[2] <= p[2] && d < 0
            wn -= 1
        end
    end
    return wn != 0
end

det2(u, v) = u[1] * v[2] - u[2] * v[1]
ctr = K2.Vec2(0.13, 0.07)

# The reference numbers (mesh, X0, Phi, exact quantities, objective values, solves) are
# frozen in data/corpus/method_fixtures/test_method_2.json (provenance in the README there).
const METHOD2_FIXTURES = let d = JSON.parsefile(joinpath(CORPUS, "method_fixtures", "test_method_2.json"))
    Dict(c["name"] => c for c in d["cases"])
end
# The MCase of a frozen entry: the frozen mesh (with its sigma) and make_case's X.
function fixture_case(fx)
    m = fixture_mesh_raw(fx["mesh"])
    c = K2.make_cut(m)
    return MCase(m, c, K2.holes_partition(c), fixture_points(fx["X_case"]))
end

# ---------------------------------------------------------------------------
# K6: the 0+ split-cut separation q_e (method/zero_plus).

@testset "q_e sign convention: positive on every split edge of a deployable pattern" begin
    # hexagons_auto has Theta_max = 2.0944 > 0 (K2a), so no split-cut duplicate may be
    # moving inward at theta = 0+: every q_e must be strictly positive.
    cs = make_case(K2.tiling_hexagons(K2.disk(ctr, 3.0)), false)
    @test K2.n_split(cs.c) > 0
    q = K2.zero_plus_q(cs.c, cs.X)
    @test length(q) == K2.n_split(cs.c)
    for v in q
        @test v > 0
    end

    # The same on every other deployable reference pattern that has split cuts.
    checked = 0
    for k in 0:2
        c2 = k == 0 ? make_case(K2.tiling_snub_square(K2.disk(ctr, 3.2)), false) :
             k == 1 ? make_case(K2.tiling_truncated_square(K2.disk(ctr, 4.0)), false) :
                      make_case(K2.tiling_3_4_3_12(K2.disk(ctr, 4.2)), false)
        K2.n_split(c2.c) == 0 && continue
        K2.has_collision(c2.c, K2.deploy(c2.c, c2.X, 0.0).Y, 1e-12) && continue
        B = K2.deploy_basis(c2.c, c2.X)
        sd = K2.swept_discs(c2.c, B)
        ov = K2.exact_theta_max_overlap(c2.c, B, K2.candidate_pairs(c2.c, sd, Float64(pi), true),
                                        1e-9, Float64(pi), 1e-9)
        ov.theta_max <= 1e-9 && continue  # a zero-range pattern is allowed inward edges
        for v in K2.zero_plus_q(c2.c, c2.X)
            @test v > 0
        end
        checked += 1
    end
    @test checked >= 2
end

# ---------------------------------------------------------------------------
# K6 pass 3: the 0+ CORNER (vertex-into-edge) margin and the objective that uses it.

@testset "corner margin decides the 0+ vertex-into-edge status, against the geometry" begin
    # mu > 0  <=>  the other copy of the vertex stays OUT of the face at theta = 0+.
    # Checked against an actual point-in-polygon test on the deployed faces at a small
    # theta, on every incidence of every reference pattern that has more than one copy.
    tested = 0; in_agree = 0; out_agree = 0; n_in = 0; n_out = 0
    for k in 0:3
        cs = k == 0 ? make_case(K2.tiling_hexagons(K2.disk(ctr, 3.0)), false) :
             k == 1 ? make_case(K2.tiling_snub_square(K2.disk(ctr, 3.2)), false) :
             k == 2 ? make_case(K2.tiling_truncated_square(K2.disk(ctr, 4.0)), false) :
                      make_case(K2.tiling_kagome(K2.disk(ctr, 2.5)), true)
        inc = K2.corner_incidences(cs.c)
        isempty(inc) && continue
        mu = K2.zero_plus_corner_margin(cs.c, cs.X)
        @test length(mu) == length(inc)
        th = 1e-7
        Y = K2.deploy(cs.c, cs.X, th).Y
        for i in eachindex(inc)
            # Only decide the ones that are not on the knife edge: mu must be well clear of
            # zero RELATIVE to its own scale |e| |dS|.
            abs(mu[i]) < 1e-3 * corner_scale(cs.c, cs.X, inc[i]) && continue
            poly = [Y[pv] for pv in cs.c.prime_faces[inc[i].face]]
            inside = strictly_inside(poly, Y[inc[i].pv_other])
            tested += 1
            if mu[i] < 0
                n_in += 1; inside && (in_agree += 1)
            else
                n_out += 1; inside || (out_agree += 1)
            end
        end
    end
    # Every reference pattern above is deployable, so all its margins separate. The
    # ENTERING half of the criterion is exercised by walking off X0 inside the null space
    # until some copy does move into a neighbouring face: the point-in-polygon test on
    # WELL-CONDITIONED corners, and "some margin < 0  =>  the structure collides at 0+",
    # which uses collision.jl rather than any part of the corner calculus.
    # The Gaussian walk is in null-space COORDINATES, and a fresh SVD null basis spans the
    # same space as the frozen one with different columns, so this block replays the walk
    # on the frozen Phi (data/corpus/method_fixtures/test_method_2.json) and reproduces
    # its tallies exactly: 60 states, 3 with an inward corner, 6/6 entering, 3 collisions.
    let
        fx = METHOD2_FIXTURES["hexagons_3.0"]
        cs = fixture_case(fx)
        Phi = fixture_matrix(fx["Phi"])
        @test size(Phi, 2) > 0
        mdim = 2 * size(Phi, 2)
        inc = K2.corner_incidences(cs.c)
        @test !isempty(inc)
        rng = K2.MT19937(99)
        G = K2.NormalDist(0.0, 0.15 * med_edge(cs.m))
        n_states = 0; n_neg_states = 0; n_coll = 0
        for trial in 1:60
            t = [K2.normal(G, rng) for _ in 1:mdim]
            Xt = K2.shape_point(cs.X, Phi, t)
            n_inverted(cs.m, Xt) != 0 && continue  # a folded flat state is a different bug
            n_states += 1
            mu = K2.zero_plus_corner_margin(cs.c, Xt)
            Y = K2.deploy(cs.c, Xt, 1e-7).Y
            any_neg = false
            for i in eachindex(inc)
                mu[i] >= -1e-2 * corner_scale(cs.c, Xt, inc[i]) && continue
                any_neg = true
                e1 = Xt[inc[i].v_next] - Xt[inc[i].v]
                e2 = Xt[inc[i].v_prev] - Xt[inc[i].v]
                sinang = det2(e1, e2) / (norm(e1) * norm(e2))
                abs(sinang) < 0.2 && continue  # sliver corner: winding number unusable
                poly = [Y[pv] for pv in cs.c.prime_faces[inc[i].face]]
                n_in += 1
                tested += 1
                strictly_inside(poly, Y[inc[i].pv_other]) && (in_agree += 1)
            end
            if any_neg
                n_neg_states += 1
                K2.has_collision(cs.c, K2.deploy(cs.c, Xt, 1e-6).Y, 1e-12) && (n_coll += 1)
            end
        end
        @test n_states > 20
        @test n_in > 5               # the negative branch is actually exercised
        @test n_neg_states > 0
        @test n_coll == n_neg_states  # a negative margin always means a 0+ collision
        @test n_states == 60          # the reference tallies on the same Phi
        @test n_neg_states == 3
        @test n_coll == 3
    end
    @test tested > 50
    @test in_agree == n_in
    @test out_agree == n_out
end

@testset "0+ repair objective: analytic gradient matches central differences" begin
    # The gradient the L-BFGS actually uses, on the full objective: split signs, face
    # areas, the corner term and the proximity term all switched on at once.
    cs = make_case(K2.tiling_hexagons(K2.disk(ctr, 3.0)), false)
    Phi = null_basis(cs)
    @test size(Phi, 2) > 0
    med = med_edge(cs.m)
    m = 2 * size(Phi, 2)

    opt = K2.ZeroPlusRepairOptions(lambda_rel = 1e-3, w_corner = 1.0, w_prox = 0.5)

    rng = K2.MT19937(4242)
    G = K2.NormalDist(0.0, 0.05 * med)
    for trial in 1:3
        t = [K2.normal(G, rng) for _ in 1:m]
        _, g = K2.zero_plus_objective(cs.c, cs.X, Phi, med, opt, t)
        @test length(g) == m
        h = 1e-6 * med
        for i in 1:m
            tp = copy(t); tm = copy(t)
            tp[i] += h
            tm[i] -= h
            fp, _ = K2.zero_plus_objective(cs.c, cs.X, Phi, med, opt, tp)
            fm, _ = K2.zero_plus_objective(cs.c, cs.X, Phi, med, opt, tm)
            fd = (fp - fm) / (2h)
            # |g - fd| <= 1e-5 * (1e-6 + |fd|)
            @test abs(g[i] - fd) <= 1e-5 * (1e-6 + abs(fd))
        end
    end
end

@testset "0+ repair with the corner term: FEASIBLE means every exact margin is positive" begin
    # Feasibility is decided by the exact q_e, a_f and mu, never by the smoothed
    # objective, so re-measuring them at the returned point must reproduce the verdict.
    ran = 0
    for k in 0:2
        cs = k == 0 ? make_case(K2.tiling_hexagons(K2.disk(ctr, 3.0)), false) :
             k == 1 ? make_case(K2.tiling_snub_square(K2.disk(ctr, 3.2)), false) :
                      make_case(K2.tiling_truncated_square(K2.disk(ctr, 4.0)), false)
        Phi = null_basis(cs)
        size(Phi, 2) == 0 && continue
        med = med_edge(cs.m)
        opt = K2.ZeroPlusRepairOptions(n_random = 2, max_iter = 200, w_corner = 1.0, w_prox = 1e-2)
        r = K2.zero_plus_repair(cs.c, cs.X, Phi, med, opt)
        ran += 1
        @test r.n_incidences == length(K2.corner_incidences(cs.c))
        mu = K2.zero_plus_corner_margin(cs.c, r.X)
        mmin = Inf
        nbad = 0
        for v in mu
            mmin = min(mmin, v)
            v <= 0 && (nbad += 1)
        end
        if !isempty(mu)
            @test isapprox(r.min_margin, mmin; rtol = 1e-9)
            @test r.n_bad_margin == nbad
        end
        if r.feasible
            @test nbad == 0
            for v in K2.zero_plus_q(cs.c, r.X)
                @test v > 0
            end
            @test n_inverted(cs.m, r.X) == 0
        end
        @info "case $k: $(length(mu)) incidences, feasible $(r.feasible), min margin $mmin, inward $nbad"
    end
    @test ran >= 2
end

@testset "q_e is independent of which incident face is called f0" begin
    # Swapping the two faces reverses BOTH d_e and dS_e, and det(-u,-v) = det(u,v).
    cs = make_case(K2.tiling_snub_square(K2.disk(ctr, 3.2)), false)
    @test K2.n_split(cs.c) > 0
    sc = K2.split_copies(cs.c)
    d = K2.deploy(cs.c, cs.X, 0.0)
    q = K2.zero_plus_q(cs.c, cs.X)
    for i in eachindex(sc)
        # Rebuild q with the roles of the two faces exchanged, using f1's stored direction.
        ed = cs.m.edges[sc[i].edge]
        h1 = cs.m.half_edges[ed.he[2]]
        pv_a = K2.prime_vertex(cs.c, h1.face, h1.from)
        pv_b = K2.prime_vertex(cs.c, sc[i].f0, h1.from)
        dS = 2.0 * (d.dY_dtheta[pv_b] - d.dY_dtheta[pv_a])
        de = cs.X[h1.to] - cs.X[h1.from]
        q_swapped = dS[1] * de[2] - dS[2] * de[1]
        @test isapprox(q_swapped, q[i]; rtol = 1e-10)
    end
end

@testset "q_e(t) is exactly the quadratic the shape-space form predicts" begin
    cs = make_case(K2.tiling_snub_square(K2.disk(ctr, 3.2)), false)
    sys = K2.assemble_system(cs.c, cs.hs, cs.m.X, K2.Fixed)
    sr = K2.solve_system(sys, cs.m.X)
    @test sr.projection_ok
    @test sr.dim_null >= 1
    X0 = K2.matrix_to_points(sr.X0)
    form = K2.zero_plus_form(cs.c, X0, sr.Phi)
    @test form.m == 2 * sr.dim_null

    rng = K2.MT19937(4242)
    G = K2.NormalDist(0.0, 0.05)
    for trial in 0:4
        t = [trial > 0 ? K2.normal(G, rng) : 0.0 for _ in 1:form.m]
        qf = K2.zero_plus_eval(form, t)
        qd = K2.zero_plus_q(cs.c, K2.shape_point(X0, sr.Phi, t))
        @test length(qd) == length(qf)
        worst = 0.0
        scale = 1e-30
        for i in eachindex(qf)
            worst = max(worst, abs(qf[i] - qd[i]))
            scale = max(scale, abs(qd[i]))
        end
        @test worst <= 1e-9 * scale
    end
end

@testset "q_e <= 0 forces immediate penetration: Theta_max = 0" begin
    # The mechanism claim behind K6: an inward split edge is not merely correlated with
    # Theta_max = 0, it causes it. Walk the null space until some q_e goes negative while
    # the flat state is still an embedding, then check the structure collides at 0+.
    cs = make_case(K2.tiling_snub_square(K2.disk(ctr, 3.2)), false)
    sys = K2.assemble_system(cs.c, cs.hs, cs.m.X, K2.Fixed)
    sr = K2.solve_system(sys, cs.m.X)
    @test sr.projection_ok
    @test sr.dim_null >= 1
    X0 = K2.matrix_to_points(sr.X0)
    for v in K2.zero_plus_q(cs.c, X0)
        @test v > 0
    end

    rng = K2.MT19937(31337)
    G = K2.NormalDist(0.0, 1.0)
    found = 0
    tested = 0
    for trial in 0:399
        found >= 3 && break
        t = [K2.normal(G, rng) for _ in 1:2sr.dim_null]
        t .*= 0.02 * (1 + trial % 40)
        X = K2.shape_point(X0, sr.Phi, t)
        q = K2.zero_plus_q(cs.c, X)
        minimum(q) >= 0 && continue
        tested += 1
        # Only judge embeddings that are still valid flat states, so the collision at 0+
        # cannot be blamed on the flat state itself.
        K2.has_collision(cs.c, K2.deploy(cs.c, X, 0.0).Y, 1e-12) && continue
        @test K2.has_collision(cs.c, K2.deploy(cs.c, X, 1e-6).Y, 1e-12)
        found += 1
    end
    @test tested > 0
    @info "snub square: $tested null-space points with an inward split edge, $found of them valid flat states"
end

# ---------------------------------------------------------------------------
# convex_embed -- the convexity-constrained shape-space point of K9.

# The 2x2 unit-square grid, built by hand so that both the convex reference embedding
# and the reflex one are exact rational data a reader can check with a pencil.
#   7---8---9        faces (CCW):  {1,2,5,4} {2,3,6,5} {4,5,8,7} {5,6,9,8}
#   | 3 | 4 |        vertex 5 is the only interior vertex.
#   4---5---6
#   | 1 | 2 |
#   1---2---3
function grid2x2()
    X = [K2.Vec2(i, j) for j in 0:2 for i in 0:2]
    m = K2.Mesh(X, [[1, 2, 5, 4], [2, 3, 6, 5], [4, 5, 8, 7], [5, 6, 9, 8]])
    m.sigma = [-1, 1, 1, -1]
    K2.build_topology!(m)
    return m
end

@testset "corner_crosses is det(X_v - X_prev, X_next - X_v), hand-checked" begin
    @testset "unit square: every cross is +1" begin
        m = K2.Mesh([K2.Vec2(0, 0), K2.Vec2(1, 0), K2.Vec2(1, 1), K2.Vec2(0, 1)], [[1, 2, 3, 4]])
        m.sigma = [-1]
        K2.build_topology!(m)
        cr = K2.corner_crosses(m, m.X)
        @test length(cr) == 4
        for v in cr
            @test isapprox(v, 1.0; rtol = 1e-12)
        end
    end
    @testset "reflex quad: the dart's tip is the one negative corner" begin
        # (0,0) (2,0) (0.5,0.5) (0,2) in CCW order: corner 3 is the reflex tip.
        m = K2.Mesh([K2.Vec2(0, 0), K2.Vec2(2, 0), K2.Vec2(0.5, 0.5), K2.Vec2(0, 2)], [[1, 2, 3, 4]])
        m.sigma = [-1]
        K2.build_topology!(m)
        cr = K2.corner_crosses(m, m.X)
        @test length(cr) == 4
        @test isapprox(cr[1], 4.0; rtol = 1e-12)   # a=(0,-2)   b=(2,0)
        @test isapprox(cr[2], 1.0; rtol = 1e-12)   # a=(2,0)    b=(-1.5,0.5)
        @test isapprox(cr[3], -2.0; rtol = 1e-12)  # a=(-1.5,0.5) b=(-0.5,1.5)
        @test isapprox(cr[4], 1.0; rtol = 1e-12)   # a=(-0.5,1.5) b=(0,-2)
    end
    @testset "the corner list matches the faces, in stored order" begin
        m = grid2x2()
        fc = K2.face_corners(m)
        @test length(fc) == 16
        for z in fc
            vs = m.faces[z.face]
            n = length(vs)
            @test vs[z.idx] == z.v
            @test vs[mod1(z.idx + 1, n)] == z.v_next
            @test vs[mod1(z.idx - 1, n)] == z.v_prev
        end
    end
end

@testset "convex_embed recovers a hand-known constrained minimiser" begin
    # One free vertex, one shape-space direction per coordinate: X(t) moves ONLY the
    # interior vertex 5, from the reflex start (1.7, 1.7) back into the square. The
    # constrained minimiser is exactly X_ini and the solver must land on it.
    m = grid2x2()
    X_ini = copy(m.X)
    X0 = copy(m.X)
    X0[5] = K2.Vec2(1.7, 1.7)
    m.X = X0
    K2.build_topology!(m)
    c = K2.make_cut(m)

    Phi = zeros(K2.n_vertices(m), 1)
    Phi[5, 1] = 1.0  # orthonormal: |X(t) - X_ini|^2 = |t|^2 + const

    # The start really is non-convex, and the hand computation says where.
    cr0 = K2.corner_crosses(m, X0)
    n_neg = count(<=(0), cr0)
    @test n_neg > 0

    opt = K2.ConvexEmbedOptions(delta_rel = 1e-3, n_random = 2, max_iter = 400)
    r = K2.convex_embed(c, X0, Phi, X_ini, 1.0, opt)

    @test r.n_corners == 16
    @test r.n_nonconvex_start == n_neg
    @test isapprox(r.min_cross_start, minimum(cr0); rtol = 1e-12)
    @test r.feasible
    # The exact constraint values at the returned point really do clear the margin.
    for v in K2.corner_crosses(m, r.X)
        @test v >= r.delta
    end
    # and the returned point is X_ini to well under the margin's own length scale.
    @test norm(r.X[5] - K2.Vec2(1, 1)) < 1e-2
    @test r.dist_ini < 1e-2
    @test r.dist_ini < r.dist_ini_x0
    for v in 1:K2.n_vertices(m)
        v != 5 && @test norm(r.X[v] - X_ini[v]) < 1e-12  # Phi moves nothing else
    end
end

@testset "convex_embed's phase-A gradient matches a central difference" begin
    cs = make_case(K2.tiling_snub_square(K2.disk(ctr, 3.2)), false)
    sys = K2.assemble_system(cs.c, cs.hs, cs.m.X, K2.Fixed)
    sr = K2.solve_system(sys, cs.m.X)
    @test sr.projection_ok
    @test sr.dim_null >= 1
    X0 = K2.matrix_to_points(sr.X0)
    med = med_edge(cs.m)
    m = 2 * sr.dim_null

    rng = K2.MT19937(90001)
    G = K2.NormalDist(0.0, 0.05 * med)
    for variant in 0:1
        # (a) convexity only, (b) + split signs
        opt = K2.ConvexEmbedOptions(delta_rel = 1e-3, split_delta_rel = variant == 1 ? 1e-3 : 0.0)
        t = [K2.normal(G, rng) for _ in 1:m]
        _, g = K2.convex_embed_objective(cs.c, X0, sr.Phi, med, opt, 1.0, t)
        h = 1e-6 * med
        worst = 0.0
        scale = 1e-30
        for i in 1:min(m, 12)
            tp = copy(t); tm = copy(t)
            tp[i] += h
            tm[i] -= h
            fp, _ = K2.convex_embed_objective(cs.c, X0, sr.Phi, med, opt, 1.0, tp)
            fm, _ = K2.convex_embed_objective(cs.c, X0, sr.Phi, med, opt, 1.0, tm)
            fd = (fp - fm) / (2h)
            worst = max(worst, abs(fd - g[i]))
            scale = max(scale, abs(fd))
        end
        @test worst <= 1e-4 * max(1.0, scale)
    end
end

@testset "convex_embed: a FEASIBLE verdict is always the exact constraint values" begin
    # On the split-bearing tilings, whatever the solver returns, the reported minima and
    # counts are the ones a direct evaluation at r.X gives -- feasibility is never decided
    # by a penalty or a barrier value.
    ran = 0
    for k in 0:2
        cs = k == 0 ? make_case(K2.tiling_snub_square(K2.disk(ctr, 3.2)), false) :
             k == 1 ? make_case(K2.tiling_truncated_square(K2.disk(ctr, 4.0)), false) :
                      make_case(K2.tiling_hexagons(K2.disk(ctr, 3.0)), false)
        Phi = null_basis(cs)
        size(Phi, 2) == 0 && continue
        med = med_edge(cs.m)
        opt = K2.ConvexEmbedOptions(delta_rel = 1e-3, split_delta_rel = 1e-3, n_random = 1,
                                    max_iter = 150, barrier_stages = 2)
        r = K2.convex_embed(cs.c, cs.X, Phi, cs.X, med, opt)
        ran += 1
        cr = K2.corner_crosses(cs.m, r.X)
        q = K2.zero_plus_q(cs.c, r.X)
        @test length(cr) == r.n_corners
        mc = Inf; nb = 0
        for v in cr
            mc = min(mc, v)
            v < r.delta && (nb += 1)
        end
        mq = Inf; nbq = 0
        for v in q
            mq = min(mq, v)
            v < r.delta_split && (nbq += 1)
        end
        @test isapprox(r.min_cross, mc; rtol = 1e-9)
        @test r.n_bad == nb
        if !isempty(q)
            @test isapprox(r.min_q, mq; rtol = 1e-9)
            @test r.n_bad_q == nbq
        end
        @test r.feasible == (nb == 0 && nbq == 0)
        if r.feasible
            # Every corner strictly convex implies every face convex and positively oriented.
            @test n_inverted(cs.m, r.X) == 0
        end
        @info "case $k: $(r.n_corners) corners, nonconvex at X0 $(r.n_nonconvex_start), feasible $(r.feasible)"
    end
    @test ran >= 2
end

@testset "convex_embed: the barrier stages only ever move CLOSER, and stay feasible" begin
    # Phase B is a refinement, not a re-solve: it starts at phase A's feasible point and
    # keeps an iterate only if the EXACT constraints still hold there and the proximity
    # has improved. So, with phase A held bit-identical (same seed, same starts):
    #   (i) dist_ini <= dist_phase_a always, with equality when no stage was kept;
    #   (ii) running with the barrier can never end further from X_ini than running
    #        without it;
    #   (iii) every returned point of a FEASIBLE solve clears the exact margins.
    ran = 0
    for k in 0:2
        cs = k == 0 ? make_case(K2.tiling_snub_square(K2.disk(ctr, 3.2)), false) :
             k == 1 ? make_case(K2.tiling_truncated_square(K2.disk(ctr, 4.0)), false) :
                      make_case(K2.tiling_hexagons(K2.disk(ctr, 3.0)), false)
        Phi = null_basis(cs)
        size(Phi, 2) == 0 && continue
        med = med_edge(cs.m)

        # convexity only: variant (a), the relaxation
        no_b = K2.ConvexEmbedOptions(delta_rel = 1e-3, split_delta_rel = 0.0, n_random = 1,
                                     max_iter = 200, seed = 4242, barrier_stages = 0)
        with_b = K2.ConvexEmbedOptions(delta_rel = 1e-3, split_delta_rel = 0.0, n_random = 1,
                                       max_iter = 200, seed = 4242, barrier_stages = 4,
                                       barrier_iter = 150)

        # Phase A is the same computation in both, so its answer must agree exactly.
        r0 = K2.convex_embed(cs.c, cs.X, Phi, cs.X, med, no_b)
        rb = K2.convex_embed(cs.c, cs.X, Phi, cs.X, med, with_b)
        ran += 1
        @test r0.feasible == rb.feasible
        @test r0.barrier_stages_kept == 0
        @test isapprox(r0.dist_ini, r0.dist_phase_a; rtol = 1e-12)
        @test isapprox(rb.dist_phase_a, r0.dist_phase_a; rtol = 1e-9)

        @test rb.barrier_stages_kept >= 0
        @test rb.barrier_stages_kept <= with_b.barrier_stages
        @test rb.dist_ini <= rb.dist_phase_a * (1 + 1e-12)
        @test rb.dist_ini <= r0.dist_ini * (1 + 1e-12)
        if rb.barrier_stages_kept == 0
            @test isapprox(rb.dist_ini, rb.dist_phase_a; rtol = 1e-12)
        end

        if rb.feasible
            # The refined point is still exactly feasible, measured directly.
            for v in K2.corner_crosses(cs.m, rb.X)
                @test v >= rb.delta
            end
            @test n_inverted(cs.m, rb.X) == 0
        end
        @info "case $k: stages kept $(rb.barrier_stages_kept), dist $(r0.dist_ini) -> $(rb.dist_ini)"
    end
    @test ran >= 2
end

# ---------------------------------------------------------------------------
# The frozen reference numbers themselves: the acceptance criterion.

@testset "zero_plus / convex_embed reproduce the frozen reference numbers" begin
    relerr(a, b) = isempty(b) ? 0.0 : maximum(abs.(a .- b)) / max(1e-300, maximum(abs.(b)))
    for name in ("hexagons_3.0", "snub_square_3.2", "truncated_square_4.0", "truncated_square_3.0",
                 "3_4_3_12_4.2", "kagome_2.5_checker")
        fx = METHOD2_FIXTURES[name]
        cs = fixture_case(fx)
        # The generator + orientation reproduce the frozen mesh and sigma.
        mg = name == "hexagons_3.0" ? K2.tiling_hexagons(K2.disk(ctr, 3.0)) :
             name == "snub_square_3.2" ? K2.tiling_snub_square(K2.disk(ctr, 3.2)) :
             name == "truncated_square_4.0" ? K2.tiling_truncated_square(K2.disk(ctr, 4.0)) :
             name == "truncated_square_3.0" ? K2.tiling_truncated_square(K2.disk(ctr, 3.0)) :
             name == "3_4_3_12_4.2" ? K2.tiling_3_4_3_12(K2.disk(ctr, 4.2)) :
                                       K2.tiling_kagome(K2.disk(ctr, 2.5))
        gc = make_case(mg, name == "kagome_2.5_checker")
        @test gc.m.faces == cs.m.faces
        @test maximum(norm.(gc.m.X .- cs.m.X)) < 1e-12
        @test gc.m.sigma == cs.m.sigma
        @test maximum(norm.(gc.X .- cs.X)) < 1e-12
        @test isapprox(med_edge(cs.m), fx["med"]; rtol = 1e-15)
        # The Eq. (6) projection and the null space (as a subspace) agree.
        sr = K2.solve_system(K2.assemble_system(cs.c, cs.hs, cs.m.X, K2.Fixed), cs.m.X)
        X0c = fixture_points(fx["X0"])
        Phic = fixture_matrix(fx["Phi"])
        @test sr.dim_null == fx["dim_null"]
        @test maximum(norm.(K2.matrix_to_points(sr.X0) .- X0c)) < 1e-12
        if sr.dim_null > 0
            @test maximum(abs.(sr.Phi * sr.Phi' .- Phic * Phic')) < 1e-12
        end
        # Exact quantities at X_case / X0.
        @test K2.n_split(cs.c) == fx["n_split"]
        @test length(K2.corner_incidences(cs.c)) == fx["n_incidences"]
        @test relerr(K2.zero_plus_q(cs.c, cs.X), Float64.(fx["q_at_case"])) < 1e-12
        @test relerr(K2.zero_plus_corner_margin(cs.c, cs.X), Float64.(fx["mu_at_case"])) < 1e-12
        @test relerr(K2.corner_crosses(cs.m, cs.X), Float64.(fx["crosses_at_case"])) < 1e-12
        haskey(fx, "probe_t") || continue
        med = fx["med"]
        t = Float64.(fx["probe_t"])
        @test relerr(K2.zero_plus_eval(K2.zero_plus_form(cs.c, X0c, Phic), t),
                     Float64.(fx["q_form_at_probe"])) < 1e-12
        F, g = K2.zero_plus_objective(cs.c, X0c, Phic, med,
                                      K2.ZeroPlusRepairOptions(lambda_rel = 1e-3, w_corner = 1.0, w_prox = 0.5), t)
        @test isapprox(F, fx["zero_plus_objective_at_probe"]; rtol = 1e-12)
        @test relerr(g, Float64.(fx["zero_plus_grad_at_probe"])) < 1e-12
        for (k, sd) in (("convex_embed_a", 0.0), ("convex_embed_b", 1e-3))
            F, g = K2.convex_embed_objective(cs.c, X0c, Phic, med,
                                             K2.ConvexEmbedOptions(delta_rel = 1e-3, split_delta_rel = sd), 1.0, t)
            @test isapprox(F, fx[k * "_objective_at_probe"]; rtol = 1e-12)
            @test relerr(g, Float64.(fx[k * "_grad_at_probe"])) < 1e-12
        end
        # The solves: same verdicts and stage counts; the numbers agree to the L-BFGS
        # rounding-path level (the reference and fresh paths diverge at ~1e-8 relative).
        s = fx["convex_embed_solve"]
        r = K2.convex_embed(cs.c, cs.X, Phic, cs.X, med,
                            K2.ConvexEmbedOptions(delta_rel = 1e-3, split_delta_rel = 1e-3, n_random = 1,
                                                  max_iter = 150, barrier_stages = 2))
        @test r.feasible == s["feasible"]
        @test r.n_nonconvex_start == s["n_nonconvex_start"]
        @test r.barrier_stages_kept == s["barrier_stages_kept"]
        @test isapprox(r.dist_ini, s["dist_ini"]; rtol = 1e-6, atol = 1e-12)
        @test isapprox(r.min_cross, s["min_cross"]; rtol = 1e-6)
        s = fx["zero_plus_repair_solve"]
        r = K2.zero_plus_repair(cs.c, cs.X, Phic, med,
                                K2.ZeroPlusRepairOptions(n_random = 2, max_iter = 200, w_corner = 1.0, w_prox = 1e-2))
        @test r.feasible == s["feasible"]
        @test isapprox(r.f_start, s["f_start"]; rtol = 1e-12)
        @test isapprox(r.min_margin_start, s["min_margin_start"]; rtol = 1e-12)
        @test r.n_bad_margin_start == s["n_bad_margin_start"]
        @test isapprox(r.f_end, s["f_end"]; rtol = 1e-6)
        @test isapprox(r.min_margin, s["min_margin"]; rtol = 1e-6)
    end
end
