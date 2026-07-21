# derivation_tests.jl -- port of code/tests/derivation_tests.cpp (Phase 7 Checker): the
# INDEPENDENT adversarial verification of derivations/core.md and lemmas.md.
#
# Everything the Deriver asserts algebraically is RE-DERIVED here from the code's own sign
# convention (kinematics.jl) and compared against the library, never against the Deriver's
# scratch programs: the face potential u_f, the harmonic coefficients (p,q,r), the
# tau-quadratic, the T5.3 closed forms, the T6.2 gradient and the L2 covector d_tau are all
# recomputed from scratch in this file.
#
# Inputs.  Every population the C++ builds (corpus(), the H-LOC patches, the T7(iv) tori,
# l1_corpus(), the L2 K7 population, the 110 fresh and the 2000 search Voronoi tori) is
# frozen in CORPUS/reference_patterns/derivation_inputs{,_l2}.json by
# freeze_derivation_inputs.cpp.  The Julia generators reproduce them (see the
# "generators reproduce the frozen inputs" testset for the exact tally: every sigma but
# hexagons_auto, every random graph and all 2122 Voronoi tori are bit-identical; four
# tilings differ by trig ulps and make_tiling_pattern differs for squares/kagome), so the
# frozen copies are what the cases below run on -- the C++ numbers are then exactly the
# reference.  Shape-space SAMPLES cannot be identical: Gaussian coefficients are drawn on
# the port's null-space basis Phi, which spans the same space as Eigen's but with different
# columns; the cases therefore carry the C++ X0 and Phi (frozen alongside the inputs, and
# checked against the port's: same X0 to 1e-9, same projector Phi Phi^T to 1e-9), so every
# Gaussian shape-space sample below IS the C++ sample and the sample tallies are asserted
# equal to the C++ ones.  The only counts that are not asserted exactly are the
# tie-sensitive ones (roots within 1e-9 of an interval endpoint, a rounding-noise
# classification), which are checked to within a few units per million.
#
# Runtime.  C++ wall times (clang -O2, arm64) are quoted per case; the whole C++ run is
# 42 s, of which R3 (the chart-free crossing truth) is 37.6 s.  The largest cases are
# gated behind ENV["KIRIGAMI_FULL_DERIVATIONS"] == "1"; the default runs a documented
# subset (fewer corpus cases / draws) that keeps every property check meaningful.
include("helpers.jl")
using Printf

const K = Kirigami
const FULL = get(ENV, "KIRIGAMI_FULL_DERIVATIONS", "0") == "1"
const DERIV_INPUTS = JSON.parsefile(joinpath(CORPUS, "reference_patterns", "derivation_inputs.json"))

# ---------------------------------------------------------------- reporting ----
struct DRow
    id::String
    what::String
    samples::Int
    max_err::Float64
    tol::Float64
end
const D_ROWS = DRow[]
record(id, what, n, err, tol) = push!(D_ROWS, DRow(id, what, Int(n), Float64(err), Float64(tol)))
function print_report()
    @printf("\n%-8s %-62s %10s %12s %12s %6s\n", "id", "identity / claim", "samples", "max error", "tolerance", "verdict")
    println("-"^114)
    for r in D_ROWS
        @printf("%-8s %-62s %10d %12.3e %12.3e %6s\n", r.id, r.what, r.samples, r.max_err, r.tol,
                r.max_err <= r.tol ? "PASS" : "FAIL")
    end
    println("-"^114)
end

# ------------------------------------------------------------------ algebra ----
Jm(v::K.Vec2) = K.Vec2(-v[2], v[1])   # rotation by +pi/2
const det2 = K.det2                    # the fma-contracted form the C++ compiled to
sq(v::K.Vec2) = dot(v, v)

# ------------------------------------------------------------------- corpus ----
mutable struct Case
    name::String
    m::K.Mesh
    c::K.CutStructure
    hs::K.HoleSet
    sys::K.LinearSystem
    rep::K.SolveReport
    X0::Vector{K.Vec2}
    k::Int          # dim of the null space
    ok::Bool        # usable (deployable X0, connected Gamma)
end

# The C++ `sample` constructs a FRESH std::normal_distribution per call (no cached
# second variate carries over), so a fresh NormalDist is made here too.
function sample(cs::Case, rng::K.MT19937, scale::Float64)
    X = copy(cs.rep.X0)
    if cs.k > 0 && scale != 0.0
        g = K.NormalDist(0.0, scale)
        T = zeros(cs.k, 2)
        for i in 1:cs.k
            T[i, 1] = K.normal(g, rng)
            T[i, 2] = K.normal(g, rng)
        end
        X += cs.rep.Phi * T
    end
    return K.matrix_to_points(X)
end
function from_T(cs::Case, T::Matrix{Float64})
    X = copy(cs.rep.X0)
    cs.k > 0 && (X += cs.rep.Phi * T)
    return K.matrix_to_points(X)
end

# `cpp` (a frozen corpus row) supplies the C++ solve_system X0 and null-space basis Phi:
# the port's own X0 must agree to 1e-9 and its Phi span the same subspace, after which the
# case carries the C++ pair so that every Gaussian shape-space sample below is the C++ one.
function build_case(m::K.Mesh, name::String, cpp = nothing)
    c = K.make_cut(m)
    hs = K.holes_partition(c)
    sys = K.assemble_system(c, hs, m.X, K.Fixed)
    rep = K.solve_system(sys, m.X)
    if cpp !== nothing
        X0c = fixture_matrix(cpp["X0"])
        Phic = fixture_matrix(cpp["Phi"])
        @test rep.dim_null == cpp["k"]
        @test maximum(abs, rep.X0 - X0c) < 1e-9
        if rep.dim_null > 0
            @test opnorm(Phic * Phic' - rep.Phi * rep.Phi') < 1e-9
            rep.Phi = Phic
        end
        rep.X0 = X0c
    end
    X0 = K.matrix_to_points(rep.X0)
    r = K.hole_residuals(c, X0, hs)
    ok = rep.projection_ok && r.max_norm < 1e-7 && K.n_hinge(c) > 0
    return Case(name, m, c, hs, sys, rep, X0, rep.dim_null, ok)
end

# The whole corpus, built once, from the frozen C++ inputs (16 of 16 candidates are ok).
const CORPUS_CASES = Case[]
function corpus()
    isempty(CORPUS_CASES) || return CORPUS_CASES
    @testset "corpus X0/Phi agree with the C++" begin
        for j in DERIV_INPUTS["corpus"]
            cs = build_case(fixture_mesh_raw(j["mesh"]), j["name"], j)
            cs.ok && push!(CORPUS_CASES, cs)
        end
    end
    return CORPUS_CASES
end

# --------------------- INDEPENDENT re-implementation of T1 (the face potential) --
# u_g - u_f = sigma_g * x_{src(e)} for the hinge edge e traversed f -> g   (T1.3),
# built by the same BFS order the code uses (seeds 1..F, u = 0 at each seed).
mutable struct Potential
    u::Vector{K.Vec2}
    reached::Vector{Bool}
    comp::Vector{Int}
    n_comp::Int
    closure_residual::Float64  # worst violation of (T1.3) on NON-tree hinge edges
end

function potential(cs::Case, X::Vector{K.Vec2})
    m = cs.m
    F = K.n_faces(m)
    P = Potential(fill(K.Vec2(0, 0), F), falses(F), zeros(Int, F), 0, 0.0)
    adj = K._hinge_adjacency(cs.c)
    for s in 1:F
        P.reached[s] && continue
        P.reached[s] = true
        P.u[s] = K.Vec2(0, 0)
        P.n_comp += 1
        P.comp[s] = P.n_comp
        q = Int[s]
        while !isempty(q)
            f = popfirst!(q)
            for (e, g) in adj[f]
                inc = m.sigma[g] * X[cs.c.hinge_dir[e].src]  # sigma_g * x_src
                if !P.reached[g]
                    P.reached[g] = true
                    P.comp[g] = P.n_comp
                    P.u[g] = P.u[f] + inc
                    push!(q, g)
                else
                    P.closure_residual = max(P.closure_residual, norm(P.u[g] - P.u[f] - inc))
                end
            end
        end
    end
    return P
end

# C, S of (T1.5): C_(v,f) = x_v, S_(v,f) = J(2 u_f - sigma_f x_v).  Built from MY u.
struct MyBasis
    C::Matrix{Float64}
    S::Matrix{Float64}
end
function my_basis(cs::Case, X::Vector{K.Vec2}, P::Potential)
    n = cs.c.n_prime_vertices
    C = zeros(n, 2)
    S = zeros(n, 2)
    for f in 1:K.n_faces(cs.m)
        vs = cs.m.faces[f]
        for kk in eachindex(vs)
            pv = cs.c.prime_faces[f][kk]
            xv = X[vs[kk]]
            s = Jm(2.0 * P.u[f] - Float64(cs.m.sigma[f]) * xv)
            C[pv, 1] = xv[1]; C[pv, 2] = xv[2]
            S[pv, 1] = s[1];  S[pv, 2] = s[2]
        end
    end
    return MyBasis(C, S)
end
row(M::Matrix{Float64}, i::Int) = K.Vec2(M[i, 1], M[i, 2])

# --------------------- INDEPENDENT harmonic coefficients (T3.2)-(T3.4) ----------
my_orient(Cab, Sab, Caw, Saw) = K.Harmonic(0.5 * (det2(Cab, Caw) + det2(Sab, Saw)),
                                           0.5 * (det2(Cab, Caw) - det2(Sab, Saw)),
                                           0.5 * (det2(Cab, Saw) + det2(Sab, Caw)))
my_dot(Cab, Sab, Caw, Saw) = K.Harmonic(0.5 * (dot(Cab, Caw) + dot(Sab, Saw)),
                                        0.5 * (dot(Cab, Caw) - dot(Sab, Saw)),
                                        0.5 * (dot(Cab, Saw) + dot(Sab, Caw)))
my_len2(Cab, Sab) = K.Harmonic(0.5 * (sq(Cab) + sq(Sab)), 0.5 * (sq(Cab) - sq(Sab)), dot(Cab, Sab))
heval(h::K.Harmonic, th) = K.harmonic_eval(h, th)
hscale(h::K.Harmonic) = K.scale(h)

# Roots of h in (lo, hi] via the tau-quadratic of (T3.5)/(T3.6), re-derived here.
function my_roots(h::K.Harmonic, lo::Float64, hi::Float64, tol::Float64)
    out = Float64[]
    scl = abs(h.p) + abs(h.q) + abs(h.r)
    scl <= tol && return out                          # identically zero (T3.H.1)
    A = h.p - h.q; B = 2 * h.r; C = h.p + h.q
    lead = max(abs(A), abs(B), abs(C))
    taus = Float64[]
    if abs(A) <= 1e-14 * lead
        abs(B) > 1e-14 * lead && push!(taus, -C / B)
        push!(out, Float64(pi))                        # root at theta = pi
    else
        disc = fma(B, B, -(4 * A * C))                # = 4 (q^2 + r^2 - p^2); contracted as the C++
        if disc >= 0
            sqd = sqrt(disc)
            push!(taus, (-B + sqd) / (2 * A))
            push!(taus, (-B - sqd) / (2 * A))
        end
    end
    for t in taus
        push!(out, 2.0 * atan(t))
    end
    keep = [th for th in out if th > lo && th <= hi + 1e-15]
    sort!(keep)
    return keep
end

# -------------------- INDEPENDENT polygon overlap (SAT, convex polygons) --------
function sat_overlap_convex(A::Vector{K.Vec2}, B::Vector{K.Vec2}, eps::Float64)
    function axes_separate(P, Q)
        n = length(P)
        for i in 1:n
            e = P[mod1(i + 1, n)] - P[i]
            norm(e) < 1e-14 && continue
            nrm = normalize(Jm(e))
            a0 = 1e300; a1 = -1e300; b0 = 1e300; b1 = -1e300
            for p in P
                d = dot(nrm, p); a0 = min(a0, d); a1 = max(a1, d)
            end
            for q in Q
                d = dot(nrm, q); b0 = min(b0, d); b1 = max(b1, d)
            end
            (a1 < b0 + eps || b1 < a0 + eps) && return true   # separated (with margin eps)
        end
        return false
    end
    return !axes_separate(A, B) && !axes_separate(B, A)
end

face_poly(cs::Case, Y::Vector{K.Vec2}, f::Int) = K.Vec2[Y[pv] for pv in cs.c.prime_faces[f]]

# signed area x2 of face f of mesh m under positions X
function area2(m::K.Mesh, X::Vector{K.Vec2}, f::Int)
    vs = m.faces[f]
    a = 0.0
    for j in eachindex(vs)
        a += det2(X[vs[j]], X[vs[mod1(j + 1, length(vs))]])
    end
    return a
end
all_faces_positive(m::K.Mesh, X::Vector{K.Vec2}) = all(f -> area2(m, X, f) > 0, 1:K.n_faces(m))

# =============================================================================
# Generators vs the frozen inputs (documents what the Julia pipeline reproduces)
# =============================================================================
@testset "generators reproduce the frozen derivation inputs" begin
    sig(j) = Int[s for s in j["orientation"]]
    faces_of(j) = [[Int(i) + 1 for i in f] for f in j["faces"]]
    coord_diff(m, j) = maximum(norm(m.X[i] - K.Vec2(j["vertices"][i]...)) for i in eachindex(m.X))
    rng = K.MT19937(20260903)
    gens = Any[("squares", () -> K.tiling_squares(K.disk(K.Vec2(0, 0), 2.6))),
               ("triangles", () -> K.tiling_triangles(K.disk(K.Vec2(0, 0), 2.2))),
               ("hexagons", () -> K.tiling_hexagons(K.disk(K.Vec2(0, 0), 3.0))),
               ("kagome", () -> K.tiling_kagome(K.disk(K.Vec2(0, 0), 2.4))),
               ("snub_square", () -> K.tiling_snub_square(K.disk(K.Vec2(0, 0), 2.2))),
               ("truncated_square", () -> K.tiling_truncated_square(K.disk(K.Vec2(0, 0), 2.8))),
               ("t3_4_3_12", () -> K.tiling_3_4_3_12(K.disk(K.Vec2(0, 0), 3.2)))]
    for s in 1:3
        push!(gens, ("delaunay$s", () -> K.largest_component(K.delaunay_of_random_points(45, 10.0, K.MT19937(1000 + s)))))
        push!(gens, ("voronoi$s", () -> K.largest_component(K.voronoi_of_random_points(40, 10.0, K.MT19937(2000 + s)))))
        push!(gens, ("quad$s", () -> K.largest_component(K.quad_dominant_random(40, 10.0, K.MT19937(3000 + s)))))
    end
    n_exact = 0
    for (i, (nm, g)) in enumerate(gens)
        m = g()
        j = DERIV_INPUTS["corpus"][i]["mesh"]
        @test DERIV_INPUTS["corpus"][i]["name"] == nm
        m.sigma = K.assign_orientation_relaxation(m, rng, 6, 400, 90).sigma
        @test m.faces == faces_of(j)
        @test coord_diff(m, j) < 1e-12      # tiling trig ulps on this (Rosetta) libm
        @test m.sigma == sig(j)             # the relaxation IS reproduced on all 16
        same_mesh_as_fixture(m, j) && (n_exact += 1)
    end
    println("  corpus: $n_exact / 16 meshes bit-identical to the C++, 16 / 16 sigmas identical")
    for (i, R) in enumerate((1.6, 2.6, 3.6, 4.6, 5.6, 7.0))
        m = K.tiling_squares(K.disk(K.Vec2(0, 0), R))
        m.sigma = K.assign_orientation_relaxation(m, K.MT19937(5), 4, 300, 60).sigma
        j = DERIV_INPUTS["hloc"][i]["mesh"]
        @test same_mesh_as_fixture(m, j) && m.sigma == sig(j)
    end
    for (mk, j) in enumerate(DERIV_INPUTS["tori"])
        m = mk == 1 ? K.torus_squares(4, 4) : mk == 2 ? K.torus_squares(6, 4) : K.torus_triangles(4, 4)
        m.sigma = K.assign_orientation_relaxation(m, K.MT19937(9), 4, 300, 60).sigma
        @test same_mesh_as_fixture(m, j["mesh"]) && m.sigma == sig(j["mesh"])
    end
    # L1: the non-periodic reference cases and 53 make_graph populations.  On the native
    # arm64 Julia (system libm = the C++'s) every sigma is reproduced; on the x86_64
    # (Rosetta) build hexagons_auto's relaxation lands on a different sigma (libm ulps).
    l1 = DERIV_INPUTS["l1_corpus"]
    i = 0
    for r in K.reference_cases()
        r.periodic && continue
        i += 1
        j = l1[i]["mesh"]
        @test r.mesh.faces == faces_of(j)
        @test coord_diff(r.mesh, j) < 1e-12
        if r.name == "hexagons_auto" && Sys.ARCH != :aarch64
            @test r.mesh.sigma != sig(j)
        else
            @test r.mesh.sigma == sig(j)
        end
    end
    for j in l1[i+1:end]
        g = K.make_graph(j["id"], 18, 46, 220)
        @test same_mesh_as_fixture(g.mesh, j["mesh"]) && g.mesh.sigma == sig(j["mesh"])
    end
end

# =============================================================================
# T1 -- trig-linear deployment
# =============================================================================

# C++: 0.030 s
@testset "T1 corpus is non-degenerate" begin
    C = corpus()
    println("corpus: $(length(C)) cases")
    for (cs, j) in zip(C, DERIV_INPUTS["corpus"])
        @printf("  %-18s F=%4d N=%4d hinge=%4d split=%3d H=%3d k=%3d\n", cs.name, K.n_faces(cs.m),
                K.n_vertices(cs.m), K.n_hinge(cs.c), K.n_split(cs.c), K.n_interior_holes(cs.hs), cs.k)
        # the C++ per-case statistics
        @test j["ok"]
        @test cs.k == j["k"]
        @test K.n_interior_holes(cs.hs) == j["H"]
        @test K.n_hinge(cs.c) == j["n_hinge"]
        @test K.n_split(cs.c) == j["n_split"]
    end
    @test length(C) >= 10
    @test length(C) == 16
end

# T1 Steps 2-8: Y_theta = cos(theta/2) C + sin(theta/2) S with C_(v,f) = x_v and
# S_(v,f) = J(2 u_f - sigma_f x_v), u from MY independent BFS over (T1.3).
# C++: 0.003 s
@testset "T1.5 closed form vs deploy() -- independent potential" begin
    rng = K.MT19937(7)
    ths = [0.0, 0.05, 0.3, 0.77, 1.2, 1.9, 2.6, 3.0, Float64(pi)]
    worst = 0.0; worst_flipped = 0.0; n = 0
    for cs in corpus(), rep in 0:5
        X = sample(cs, rng, rep == 0 ? 0.0 : 0.15)
        P = potential(cs, X)
        B = my_basis(cs, X, P)
        for th in ths
            D = K.deploy(cs.c, X, th, 1)
            c = cos(th / 2); s = sin(th / 2)
            for i in 1:cs.c.n_prime_vertices
                y = c * row(B.C, i) + s * row(B.S, i)
                yf = c * row(B.C, i) - s * row(B.S, i)   # the flipped convention
                worst = max(worst, norm(y - D.Y[i]))
                worst_flipped = max(worst_flipped, norm(yf - D.Y[i]))
                n += 1
            end
        end
    end
    record("T1-a", "Y(th) = cos(th/2) C + sin(th/2) S vs deploy()", n, worst, 1e-11)
    @test worst <= 1e-11
    @printf("  sign-flip control (S -> -S): max deviation %.3e (must be O(1))\n", worst_flipped)
    @test worst_flipped > 1e-3
end

# The code's DeployBasis (C = Y(0), S = 2 dY/dtheta|_0) must equal my algebraic C, S.
# C++: 0.001 s
@testset "T1.5 my (C,S) == deploy_basis()" begin
    rng = K.MT19937(11)
    worst = 0.0; n = 0
    for cs in corpus(), rep in 0:7
        X = sample(cs, rng, rep == 0 ? 0.0 : 0.2)
        P = potential(cs, X)
        B = my_basis(cs, X, P)
        DB = K.deploy_basis(cs.c, X)
        for i in 1:cs.c.n_prime_vertices
            worst = max(worst, norm(row(B.C, i) - K.basis_c(DB, i)))
            worst = max(worst, norm(row(B.S, i) - K.basis_s(DB, i)))
            n += 1
        end
    end
    record("T1-b", "my (C,S) from (T1.3)/(T1.5) == library deploy_basis()", n, worst, 1e-11)
    @test worst <= 1e-11
end

# T1 Step 5/6: closure of the 1-form on Gamma <=> Eq. (2).  Tested BOTH ways.
# C++: 0.002 s
@testset "T1 Step 5/6 closure of the potential <=> Eq. (2)" begin
    rng = K.MT19937(13)
    worst_in = 0.0; min_out = 1e300; worst_link = 0.0; n_in = 0; n_out = 0
    for cs in corpus()
        K.n_hinge(cs.c) == 0 && continue
        for rep in 0:11
            X = sample(cs, rng, rep == 0 ? 0.0 : 0.25)
            worst_in = max(worst_in, potential(cs, X).closure_residual)
            n_in += 1
        end
        # deliberately leave the shape space
        g = K.NormalDist(0.0, 0.3)
        for rep in 0:11
            X = [x + K.Vec2(K.normal(g, rng), K.normal(g, rng)) for x in cs.X0]
            R = K.hole_residuals(cs.c, X, cs.hs)
            R.max_norm < 1e-6 && continue
            P = potential(cs, X)
            min_out = min(min_out, P.closure_residual)
            th = 1.1
            mm = K.deploy(cs.c, X, th, 1).max_mismatch
            # T1 Step 3: t_f = 2 sin(th/2) J u_f, so a potential defect du maps to
            # a position mismatch of 2|sin(th/2)| |du|.
            worst_link = max(worst_link, abs(mm - 2 * abs(sin(th / 2)) * P.closure_residual) / max(1.0, mm))
            n_out += 1
        end
    end
    record("T1-c", "X in shape space => closure residual 0 on every non-tree edge", n_in, worst_in, 1e-10)
    @test worst_in <= 1e-10
    @printf("  off-shape-space controls: %d samples, smallest closure defect %.3e\n", n_out, min_out)
    @test min_out > 1e-6
    @test n_out == 192   # C++: 192 samples
    record("T1-d", "deploy() max_mismatch == 2|sin(th/2)| * closure defect (rel)", n_out, worst_link, 1e-9)
    @test worst_link <= 1e-9
end

# C++: 0.072 s
@testset "T1.A ellipse, T1.B split translate, T1.C hinge angle, T1.D rigid faces" begin
    rng = K.MT19937(17)
    ths = [0.13, 0.55, 1.0, 1.77, 2.4, 3.0]
    e_ell = 0.0; e_split = 0.0; e_off = 0.0; e_hinge = 0.0; e_area = 0.0
    n_ell = 0; n_split = 0; n_hinge = 0; n_area = 0
    for cs in corpus(), rep in 0:4
        X = sample(cs, rng, rep == 0 ? 0.0 : 0.15)
        P = potential(cs, X)
        B = my_basis(cs, X, P)
        for th in ths
            D = K.deploy(cs.c, X, th, 1)
            s = sin(th / 2)
            # T1.A -- every copy lies on the ellipse y^T (A A^T)^-1 y = 1, A = [C|S].
            for i in 1:cs.c.n_prime_vertices
                A = [row(B.C, i) row(B.S, i)]
                abs(det(A)) < 1e-6 && continue
                M = inv(A * A')
                e_ell = max(e_ell, abs(dot(D.Y[i], M * D.Y[i]) - 1.0))
                n_ell += 1
            end
            # T1.B -- both duplicates of a split edge are the SAME vector R_f d,
            # offset by exactly 2 sin(th/2) J (u_g - u_f).
            for e in cs.c.split_edges
                ed = cs.m.edges[e]
                f = cs.m.half_edges[ed.he[1]].face; g = cs.m.half_edges[ed.he[2]].face
                a = ed.key.a; b = ed.key.b
                d1 = D.Y[K.prime_vertex(cs.c, f, b)] - D.Y[K.prime_vertex(cs.c, f, a)]
                d2 = D.Y[K.prime_vertex(cs.c, g, b)] - D.Y[K.prime_vertex(cs.c, g, a)]
                e_split = max(e_split, norm(d1 - d2))
                off = D.Y[K.prime_vertex(cs.c, g, a)] - D.Y[K.prime_vertex(cs.c, f, a)]
                e_off = max(e_off, norm(off - 2 * s * Jm(P.u[g] - P.u[f])))
                n_split += 1
            end
            # T1.C -- every hinge opens by exactly theta.
            for e in cs.c.hinge_edges
                ed = cs.m.edges[e]
                f = cs.m.half_edges[ed.he[1]].face; g = cs.m.half_edges[ed.he[2]].face
                src = cs.c.hinge_dir[e].src; dst = cs.c.hinge_dir[e].dst
                pf = D.Y[K.prime_vertex(cs.c, f, src)]; pg = D.Y[K.prime_vertex(cs.c, g, src)]
                @test norm(pf - pg) < 1e-9
                v1 = D.Y[K.prime_vertex(cs.c, f, dst)] - pf
                v2 = D.Y[K.prime_vertex(cs.c, g, dst)] - pg
                (norm(v1) < 1e-9 || norm(v2) < 1e-9) && continue
                ang = atan(det2(v1, v2), dot(v1, v2))
                e_hinge = max(e_hinge, abs(abs(ang) - th))
                n_hinge += 1
            end
            # T1.D -- signed face areas are theta-independent.
            Y0 = K.deploy(cs.c, X, 0.0, 1).Y
            for f in 1:K.n_faces(cs.m)
                pol = face_poly(cs, D.Y, f)
                pol0 = face_poly(cs, Y0, f)
                A1 = 0.0; A0 = 0.0
                for i in eachindex(pol)
                    A1 += det2(pol[i], pol[mod1(i + 1, length(pol))])
                    A0 += det2(pol0[i], pol0[mod1(i + 1, length(pol0))])
                end
                e_area = max(e_area, abs(0.5 * (A1 - A0)))
                n_area += 1
            end
        end
    end
    record("T1-A", "T1.A  y^T (A A^T)^-1 y = 1 (centred ellipse)", n_ell, e_ell, 1e-8)
    record("T1-B1", "T1.B  split duplicates are the same vector", n_split, e_split, 1e-11)
    record("T1-B2", "T1.8  offset = 2 sin(th/2) J (u_g - u_f)", n_split, e_off, 1e-11)
    record("T1-C", "T1.C  |signed hinge angle| = theta", n_hinge, e_hinge, 1e-11)
    record("T1-D", "T1.D  face signed area is theta-independent", n_area, e_area, 1e-11)
    @test e_ell <= 1e-8; @test e_split <= 1e-11; @test e_off <= 1e-11
    @test e_hinge <= 1e-11; @test e_area <= 1e-11
end

# =============================================================================
# T2 -- no locking:  A(Y_theta) sigma = 0 for every theta
# =============================================================================

# C++: 0.0003 s
@testset "T2.2 pin equation (T2.1) with the explicit w of (T2.3)" begin
    rng = K.MT19937(23)
    worst = 0.0; scale_seen = 0.0; n = 0
    for cs in corpus(), rep in 0:3
        X = sample(cs, rng, rep == 0 ? 0.0 : 0.2)
        P = potential(cs, X)
        for k in 1:6
            th = K.uniform_real(rng, -4.0, 4.0)
            c = cos(th / 2); s = sin(th / 2)
            # omega_f = -sigma_f/2 ; w_f = c J u_f - sigma_f s u_f   (T2.3)
            for e in cs.c.hinge_edges
                ed = cs.m.edges[e]
                f = cs.m.half_edges[ed.he[1]].face; g = cs.m.half_edges[ed.he[2]].face
                v = cs.c.hinge_dir[e].src
                sf = Float64(cs.m.sigma[f]); sg = Float64(cs.m.sigma[g])
                wf = c * Jm(P.u[f]) - sf * s * P.u[f]
                wg = c * Jm(P.u[g]) - sg * s * P.u[g]
                of = -sf / 2; og = -sg / 2
                pe = c * X[v] + s * Jm(2.0 * P.u[g] - sg * X[v])
                res = (wf - wg) + (of - og) * Jm(pe)
                worst = max(worst, norm(res))
                scale_seen = max(scale_seen, norm(pe) + norm(P.u[f]))
                n += 1
            end
        end
    end
    record("T2-a", "(T2.1) residual with omega=-sigma/2, w of (T2.3), all theta", n, worst, 1e-10)
    @printf("  (T2-a scale of the terms involved: %.3e)\n", scale_seen)
    @test worst <= 1e-10
end

# C++: 0.004 s
@testset "T2 assembled A(theta) sigma = 0, and the pencil (T2.5)" begin
    rng = K.MT19937(29)
    worst_sigma = 0.0; worst_pencil = 0.0; n_sig = 0; n_pen = 0
    for cs in corpus()
        g = K.build_hinge_graph(cs.c)
        K.n_cycles(g) == 0 && continue
        for rep in 0:2
            X = sample(cs, rng, rep == 0 ? 0.0 : 0.2)
            sig = Float64.(cs.m.sigma)
            # A_c and A_s from the two extreme configurations
            Y0 = K.deploy(cs.c, X, 0.0, 1).Y
            Ypi = K.deploy(cs.c, X, Float64(pi), 1).Y
            Ac = Matrix(K.build_A(g, K.pins_deployed(cs.c, g, Y0)))
            As = Matrix(K.build_A(g, K.pins_deployed(cs.c, g, Ypi)))
            for k in 0:9
                th = k == 0 ? 0.0 : K.uniform_real(rng, 0.0, Float64(pi))
                Y = K.deploy(cs.c, X, th, 1).Y
                A = Matrix(K.build_A(g, K.pins_deployed(cs.c, g, Y)))
                nrm = max(1.0, maximum(abs, A))
                worst_sigma = max(worst_sigma, maximum(abs, A * sig) / nrm)
                n_sig += 1
                Ap = cos(th / 2) * Ac + sin(th / 2) * As
                worst_pencil = max(worst_pencil, maximum(abs, A - Ap) / nrm)
                n_pen += 1
            end
        end
    end
    record("T2-b", "||A(Y_theta) sigma||_inf / ||A||_inf", n_sig, worst_sigma, 1e-12)
    record("T2-c", "(T2.5)  A(th) = cos(th/2) A_c + sin(th/2) A_s", n_pen, worst_pencil, 1e-12)
    @test worst_sigma <= 1e-12
    @test worst_pencil <= 1e-12
end

# =============================================================================
# T3 -- harmonic predicates
# =============================================================================

# C++: 0.001 s
@testset "T3.2-T3.4 coefficient formulas vs direct evaluation" begin
    rng = K.MT19937(31)
    e_or = 0.0; e_dt = 0.0; e_l2 = 0.0; e_lib = 0.0; n = 0
    for cs in corpus()
        nv = cs.c.n_prime_vertices
        for rep in 0:3
            X = sample(cs, rng, rep == 0 ? 0.0 : 0.2)
            P = potential(cs, X)
            B = my_basis(cs, X, P)
            DB = K.deploy_basis(cs.c, X)
            for k in 1:40
                a = K.uniform_int(rng, 0, nv - 1) + 1
                b = K.uniform_int(rng, 0, nv - 1) + 1
                w = K.uniform_int(rng, 0, nv - 1) + 1
                (a == b || a == w || b == w) && continue
                Cab = row(B.C, b) - row(B.C, a); Sab = row(B.S, b) - row(B.S, a)
                Caw = row(B.C, w) - row(B.C, a); Saw = row(B.S, w) - row(B.S, a)
                ho = my_orient(Cab, Sab, Caw, Saw)
                hd = my_dot(Cab, Sab, Caw, Saw)
                hl = my_len2(Cab, Sab)
                lo = K.orient_harmonic(DB, a, b, w)
                e_lib = max(e_lib, abs(ho.p - lo.p), abs(ho.q - lo.q), abs(ho.r - lo.r))
                for j in 1:5
                    th = K.uniform_real(rng, -Float64(pi), Float64(pi))
                    c = cos(th / 2); s = sin(th / 2)
                    dab = c * Cab + s * Sab; daw = c * Caw + s * Saw
                    sc = max(1.0, hscale(ho))
                    e_or = max(e_or, abs(det2(dab, daw) - heval(ho, th)) / sc)
                    e_dt = max(e_dt, abs(dot(dab, daw) - heval(hd, th)) / max(1.0, hscale(hd)))
                    e_l2 = max(e_l2, abs(sq(dab) - heval(hl, th)) / max(1.0, hscale(hl)))
                    n += 1
                end
            end
        end
    end
    record("T3-a", "(T3.2) det harmonic p+q cos+r sin (relative)", n, e_or, 1e-12)
    record("T3-b", "(T3.3) dot harmonic (relative)", n, e_dt, 1e-12)
    record("T3-c", "(T3.4) squared-length harmonic (relative)", n, e_l2, 1e-12)
    record("T3-d", "my (p,q,r) == library orient_harmonic()", n, e_lib, 1e-10)
    @test e_or <= 1e-12; @test e_dt <= 1e-12; @test e_l2 <= 1e-12; @test e_lib <= 1e-10
end

# C++: 0.0003 s
@testset "T3.3 intra-face triples have q = r = 0 and p = flat area" begin
    rng = K.MT19937(37)
    e_q = 0.0; e_r = 0.0; e_p = 0.0; e_q2 = 0.0; e_r2 = 0.0; n = 0
    for cs in corpus(), rep in 0:3
        X = sample(cs, rng, rep == 0 ? 0.0 : 0.2)
        P = potential(cs, X)
        B = my_basis(cs, X, P)
        for f in 1:K.n_faces(cs.m)
            pf = cs.c.prime_faces[f]
            length(pf) < 3 && continue
            for i in 1:length(pf)-2
                a = pf[i]; b = pf[i+1]; w = pf[i+2]
                Cab = row(B.C, b) - row(B.C, a); Sab = row(B.S, b) - row(B.S, a)
                Caw = row(B.C, w) - row(B.C, a); Saw = row(B.S, w) - row(B.S, a)
                h = my_orient(Cab, Sab, Caw, Saw)
                l = my_len2(Cab, Sab)
                sc = max(1.0, hscale(h))
                e_q = max(e_q, abs(h.q) / sc)
                e_r = max(e_r, abs(h.r) / sc)
                e_p = max(e_p, abs(h.p - det2(Cab, Caw)) / sc)
                e_q2 = max(e_q2, abs(l.q) / max(1.0, hscale(l)))
                e_r2 = max(e_r2, abs(l.r) / max(1.0, hscale(l)))
                n += 1
            end
        end
    end
    record("T3-e", "(T3.3) intra-face det harmonic: q = 0", n, e_q, 1e-12)
    record("T3-f", "(T3.3) intra-face det harmonic: r = 0", n, e_r, 1e-12)
    record("T3-g", "(T3.3) intra-face p = flat signed area x2", n, e_p, 1e-12)
    record("T3-h", "(T3.3') intra-face edge length constant: q'' = r'' = 0", n, max(e_q2, e_r2), 1e-12)
    @test e_q <= 1e-12; @test e_r <= 1e-12; @test e_p <= 1e-12
    @test max(e_q2, e_r2) <= 1e-12
end

# C++: 0.001 s
@testset "T3.4 tau-quadratic roots and the discriminant condition" begin
    rng = K.MT19937(41)
    e_root = 0.0; e_lib = 0.0; n = 0; n_disc_ok = 0; n_disc_bad = 0; n_near_pi = 0
    function dedup(v)
        o = Float64[]
        for x in v
            (isempty(o) || abs(x - o[end]) > 1e-9) && push!(o, x)
        end
        return o
    end
    for cs in corpus()
        nv = cs.c.n_prime_vertices
        for rep in 0:3
            X = sample(cs, rng, rep == 0 ? 0.0 : 0.2)
            P = potential(cs, X)
            B = my_basis(cs, X, P)
            for k in 1:60
                a = K.uniform_int(rng, 0, nv - 1) + 1
                b = K.uniform_int(rng, 0, nv - 1) + 1
                w = K.uniform_int(rng, 0, nv - 1) + 1
                (a == b || a == w || b == w) && continue
                Cab = row(B.C, b) - row(B.C, a); Sab = row(B.S, b) - row(B.S, a)
                Caw = row(B.C, w) - row(B.C, a); Saw = row(B.S, w) - row(B.S, a)
                h = my_orient(Cab, Sab, Caw, Saw)
                sc = max(1e-12, hscale(h))
                # real-root criterion p^2 <= q^2 + r^2   <=>   |p| <= amplitude
                # (the C++ `q*q + r*r` is contracted to fma(q, q, r*r); it decides the
                # exactly-degenerate double roots of the flat tilings)
                has_real = h.p * h.p <= fma(h.q, h.q, h.r * h.r)
                roots = my_roots(h, -Float64(pi), Float64(pi), 1e-14 * sc)
                has_real ? (n_disc_ok += 1) : (n_disc_bad += 1)
                has_real || @test isempty(roots)
                for th in roots
                    e_root = max(e_root, abs(heval(h, th)) / sc)
                    n += 1
                end
                # Library agreement on (0, pi], comparing root SETS (a double root is listed
                # twice by the amplitude/phase route, once by the tau-quadratic); the tau
                # chart is singular at theta = pi, so a neighbourhood of pi is excluded and
                # the |p - q| ~ 0 cases are counted separately.
                HI = Float64(pi) - 1e-3
                if abs(h.p - h.q) <= 1e-9 * sc
                    n_near_pi += 1
                    continue
                end
                mine = dedup(my_roots(h, 1e-9, HI, 1e-14 * sc))
                lib = dedup(K.harmonic_roots(h, 1e-9, HI, 1e-14 * sc))
                if length(mine) == length(lib)
                    for i in eachindex(mine)
                        e_lib = max(e_lib, abs(mine[i] - lib[i]))
                    end
                else
                    e_lib < 1.0 && @printf("  [T3-j disagreement] p=%.17g q=%.17g r=%.17g  mine=%d lib=%d\n",
                                           h.p, h.q, h.r, length(mine), length(lib))
                    e_lib = max(e_lib, 1.0)
                end
            end
        end
    end
    record("T3-i", "(T3.6) every tau-quadratic root satisfies h(theta)=0 (rel)", n, e_root, 1e-9)
    record("T3-j", "my tau-quadratic roots == library harmonic_roots on (0,pi-1e-3]", n, e_lib, 1e-7)
    @printf("  (real-root criterion: %d with p^2<=q^2+r^2, %d without; %d cases with\n   |p-q| <= 1e-9*scale excluded)\n",
            n_disc_ok, n_disc_bad, n_near_pi)
    @test (n_disc_ok, n_disc_bad, n_near_pi) == (2711, 971, 54)   # C++
    @test e_root <= 1e-9; @test e_lib <= 1e-7
end

# =============================================================================
# T4 -- contact calculus
# =============================================================================

# Independent confirmation of the graze: the hexagon pattern has a contact at
# theta_1 but NO interior overlap until Theta_max.  The overlap probe here is my own
# separating-axis test, not collision.jl's.
# C++: 0.011 s
@testset "T4.2'' the graze: contact at theta_1 is not an overlap" begin
    checked = 0
    grazes = String[]
    for cs in corpus()
        K.n_split(cs.c) == 0 && continue
        # convexity screen -- SAT needs convex polygons
        convex = true
        for f in 1:K.n_faces(cs.m)
            vs = cs.m.faces[f]
            n = length(vs)
            for i in 1:n
                a = cs.X0[vs[i]]; b = cs.X0[vs[mod1(i + 1, n)]]; c2 = cs.X0[vs[mod1(i + 2, n)]]
                if det2(b - a, c2 - b) < -1e-12
                    convex = false; break
                end
            end
            convex || break
        end
        convex || continue
        B = K.deploy_basis(cs.c, cs.X0)
        sd = K.swept_discs(cs.c, B)
        pairs = K.candidate_pairs(cs.c, sd, Float64(pi), false)
        cand = K.contact_angles(cs.c, B, pairs, 1e-9, Float64(pi), 1e-9)
        rep = K.exact_theta_max_overlap(cs.c, B, pairs, 1e-9, Float64(pi), 1e-9)
        isempty(cand) && continue
        th1 = cand[1]
        # my own overlap scan on (0, Theta_max) and just past it
        function my_overlap_at(th)
            Y = K.deploy(cs.c, cs.X0, th, 1).Y
            sc = maximum(norm, Y)
            for (f, g) in pairs
                sat_overlap_convex(face_poly(cs, Y, f), face_poly(cs, Y, g), 1e-9 * max(1.0, sc)) && return true
            end
            return false
        end
        Tm = rep.theta_max
        n_clear = 0
        any_overlap_before = false
        if Tm > 1e-6
            for i in 1:199
                th = Tm * i / 200.0
                if my_overlap_at(th)
                    any_overlap_before = true; break
                end
                n_clear += 1
            end
        end
        overlap_after = (Tm < pi - 1e-3) ? my_overlap_at(min(Float64(pi), Tm + 1e-3)) : false
        @printf("  %-18s |C|=%3d theta_1=%.6f Theta_max=%.6f  clear-below=%d/%d  overlap just after=%d\n",
                cs.name, length(cand), th1, Tm, n_clear, 199, Int(overlap_after))
        @test !any_overlap_before
        th1 < Tm - 1e-6 && push!(grazes, cs.name)
        # the C++ numbers (deterministic: X0 only)
        if cs.name == "hexagons"
            @test length(cand) == 3
            @test isapprox(th1, 1.047198; atol = 1e-6)
            @test isapprox(Tm, 2.094395; atol = 1e-6)
            @test overlap_after
        elseif cs.name == "snub_square"
            @test length(cand) == 11
            @test isapprox(th1, 0.834111; atol = 1e-6)
            @test isapprox(Tm, 0.834111; atol = 1e-6)
        elseif cs.name == "truncated_square"
            @test length(cand) == 4
            @test isapprox(Tm, 2.356194; atol = 1e-6)
        end
        checked += 1
    end
    println("  (independent SAT overlap scan on $checked convex-face patterns with split cuts)")
    @test checked >= 1
    @test checked == 3
    @test grazes == ["hexagons"]   # GRAZE CONFIRMED on the hexagons only
end

# T4.4 -- beta_e = 2 pi - alpha_f - alpha_g, recomputed by hand, and Theta_max = min beta
# on split-free patterns.
# C++: 0.007 s
@testset "T4.4 the beta bound, recomputed independently" begin
    e_beta = 0.0; e_split_free = 0.0; n_beta = 0; n_sf = 0
    for cs in corpus()
        lib_beta = K.hinge_beta(cs.c, cs.X0)
        minb = 1e300
        for (i, e) in enumerate(cs.c.hinge_edges)
            ed = cs.m.edges[e]
            f = cs.m.half_edges[ed.he[1]].face; g = cs.m.half_edges[ed.he[2]].face
            v = cs.c.hinge_dir[e].src
            # interior angle of a face at v, from the flat embedding, computed here
            function ang(face)
                vs = cs.m.faces[face]
                n = length(vs)
                idx = findlast(==(v), vs)
                # faces are stored CCW, so the INTERIOR angle is the CCW angle from the
                # outgoing edge to the incoming one, in [0, 2 pi) -- reflex corners included.
                a = cs.X0[vs[mod1(idx - 1, n)]] - cs.X0[v]
                b = cs.X0[vs[mod1(idx + 1, n)]] - cs.X0[v]
                t = atan(det2(b, a), dot(b, a))
                t < 0 && (t += 2pi)
                return t
            end
            beta = 2pi - ang(f) - ang(g)
            e_beta = max(e_beta, abs(beta - lib_beta[i]))
            minb = min(minb, beta)
            n_beta += 1
        end
        if K.n_split(cs.c) == 0 && minb < 1e299
            B = K.deploy_basis(cs.c, cs.X0)
            sd = K.swept_discs(cs.c, B)
            pairs = K.candidate_pairs(cs.c, sd, Float64(pi), false)
            rep = K.exact_theta_max_overlap(cs.c, B, pairs, 1e-9, Float64(pi), 1e-9)
            expect = min(minb, Float64(pi))
            @printf("  split-free %-18s Theta_max=%.9f  min beta=%.9f  (capped at pi: %.9f)\n",
                    cs.name, rep.theta_max, minb, expect)
            e_split_free = max(e_split_free, abs(rep.theta_max - expect))
            n_sf += 1
        end
    end
    record("T4-a", "beta_e = 2pi - alpha_f - alpha_g recomputed by hand", n_beta, e_beta, 1e-12)
    record("T4-b", "split-free: Theta_max = min(min_e beta_e, pi)", n_sf, e_split_free, 1e-6)
    @test e_beta <= 1e-12
    @test e_split_free <= 1e-6
    @test n_sf == 4   # squares, triangles, kagome, t3_4_3_12
    # the 4.8.8 arithmetic of T4.Check, done by hand: 2pi - 3pi/4 - pi/2 = 3pi/4
    @test abs((2pi - 3pi / 4 - pi / 2) - 3pi / 4) < 1e-15
end

# T4.5b -- the swept radius.  THREE separate statements, only one of which is what
# K2c and contact.jl actually use.
# C++: 0.022 s
@testset "T4.5b swept radius: raw frame vs the face's own frame" begin
    e_exact_centroid = 0.0; e_rhomax_circum = 0.0; worst_ratio_raw = 0.0; e_sound = 0.0
    n_raw = 0; n_viol_raw = 0; n_c = 0
    worst_centroid_drift = 0.0
    rng = K.MT19937(131)
    for cs in corpus(), rep in 0:5
        Xs = sample(cs, rng, rep == 0 ? 0.0 : 0.2)
        B = K.deploy_basis(cs.c, Xs)
        sd = K.swept_discs(cs.c, B)
        for f in 1:K.n_faces(cs.m)
            gcf = K.Vec2(sd.gc[f, 1], sd.gc[f, 2]); gsf = K.Vec2(sd.gs[f, 1], sd.gs[f, 2])
            # (1) rho_max (the K2c / contact.jl form, in the FACE'S OWN frame) == circumradius
            e_rhomax_circum = max(e_rhomax_circum, abs(sd.rho_max[f] - sd.circum[f]))
            for pv in cs.c.prime_faces[f]
                # (2) the true swept radius about the MOVING centroid is exactly the circumradius
                mx = 0.0
                for i in 0:100
                    th = pi * i / 100.0
                    c = cos(th / 2); s = sin(th / 2)
                    y = c * K.basis_c(B, pv) + s * K.basis_s(B, pv)
                    g = c * gcf + s * gsf
                    mx = max(mx, norm(y - g))
                end
                rad_u = norm(K.basis_c(B, pv) - gcf)
                e_exact_centroid = max(e_exact_centroid, abs(mx - rad_u))
                n_c += 1
                # (3) the RAW frame (the one T4.5b states): max(|C|,|S|) is NOT a bound
                raw_max = max(norm(K.basis_c(B, pv)), norm(K.basis_s(B, pv)))
                smax = svdvals([K.basis_c(B, pv) K.basis_s(B, pv)])[1]
                traj = 0.0
                for i in 0:200
                    th = pi * i / 200.0
                    traj = max(traj, norm(cos(th / 2) * K.basis_c(B, pv) + sin(th / 2) * K.basis_s(B, pv)))
                end
                e_sound = max(e_sound, traj - smax)                 # must be <= 0
                if traj > raw_max * (1 + 1e-12)
                    n_viol_raw += 1
                    worst_ratio_raw = max(worst_ratio_raw, traj / max(1e-12, raw_max))
                end
                n_raw += 1
            end
            # (4) H-LOC as the O(n) argument needs it: displacement of the moving centroid
            #     from the FLAT centroid, in units of the circumradius.
            for i in 0:40
                th = pi * i / 40.0
                g = cos(th / 2) * gcf + sin(th / 2) * gsf
                worst_centroid_drift = max(worst_centroid_drift, norm(g - gcf) / max(1e-9, sd.circum[f]))
            end
        end
    end
    record("T4-c", "rho_max (face frame) == circumradius, EXACTLY", n_c, e_rhomax_circum, 1e-12)
    record("T4-d", "max_th |y_u(th)-gamma_f(th)| == |x_u - xbar_f|, EXACTLY", n_c, e_exact_centroid, 1e-9)
    record("T4-e", "sigma_max([C|S]) is a sound bound on |y(th)| (raw frame)", n_raw, max(0.0, e_sound), 1e-9)
    @printf("  raw frame: max(|C|,|S|) exceeded on %d / %d copies (%.0f%%), worst ratio %.4f\n",
            n_viol_raw, n_raw, 100.0 * n_viol_raw / max(1, n_raw), worst_ratio_raw)
    @printf("  H-LOC (the version the O(n) count needs): worst |gamma_f(th) - c_f| / r_f = %.2f\n", worst_centroid_drift)
    @test n_raw == 15132                       # C++: 9187 / 15132, ratio 1.4142, drift 51.05
    @test n_viol_raw == 9187
    @test isapprox(worst_ratio_raw, 1.4142; atol = 1e-4)
    @test isapprox(worst_centroid_drift, 51.05; atol = 0.01)
    @test e_rhomax_circum <= 1e-12
    @test e_exact_centroid <= 1e-9
    @test e_sound <= 1e-9
    @test worst_ratio_raw <= sqrt(2.0) + 1e-6
end

# H-LOC, decisively: does the moving centroid stay within O(r_f) of the flat centroid
# as the patch grows?  If it does not, the O(n) active-set argument of T4.5 fails.
# C++: 0.023 s
@testset "H-LOC counterexample search: centroid drift vs patch size" begin
    println("  patch      F     diam   max|gamma_f(th)-c_f|/r_f   max rho_max/r_f")
    prev = 0.0
    grows = false
    cpp_drift = Dict(1.6 => 3.601, 2.6 => 5.748, 3.6 => 8.053, 4.6 => 10.229, 5.6 => 12.404, 7.0 => 15.884)
    cpp_F = Dict(1.6 => 12, 2.6 => 24, 3.6 => 44, 4.6 => 68, 5.6 => 96, 7.0 => 156)
    for j in DERIV_INPUTS["hloc"]
        R = j["R"]
        cs = build_case(fixture_mesh_raw(j["mesh"]), "sq")
        cs.ok || continue
        B = K.deploy_basis(cs.c, cs.X0)
        sd = K.swept_discs(cs.c, B)
        drift = 0.0; rr = 0.0; diam = 0.0
        for f in 1:K.n_faces(cs.m)
            cf = K.Vec2(sd.gc[f, 1], sd.gc[f, 2])
            diam = max(diam, norm(cf))
            rr = max(rr, sd.rho_max[f] / max(1e-9, sd.circum[f]))
            for i in 0:40
                th = pi * i / 40.0
                g = cos(th / 2) * cf + sin(th / 2) * K.Vec2(sd.gs[f, 1], sd.gs[f, 2])
                drift = max(drift, norm(g - cf) / max(1e-9, sd.circum[f]))
            end
        end
        @printf("  R=%-5.1f %5d  %6.2f   %22.3f   %15.6f\n", R, K.n_faces(cs.m), 2 * diam, drift, rr)
        @test K.n_faces(cs.m) == cpp_F[R]
        @test isapprox(drift, cpp_drift[R]; atol = 1e-3)
        @test isapprox(rr, 1.0; atol = 1e-6)
        (drift > prev * 1.3 && prev > 0) && (grows = true)
        prev = drift
    end
    println("  => centroid drift grows with patch diameter: ", grows ? "YES" : "no")
    @test grows  # H-LOC in the form the O(n) argument needs is FALSE
end

# =============================================================================
# T5 -- the usable region
# =============================================================================

# C++: 0.149 s
@testset "T5.1 positive orientation is quadratic in t, and is NOT a validity certificate" begin
    rng = K.MT19937(53)
    e_quad = 0.0; n = 0; n_pos = 0; n_pos_zero_range = 0
    for cs in corpus()
        cs.k == 0 && continue
        g = K.NormalDist(0.0, 0.25)   # persists over the reps (cached second variate)
        for rep in 0:11
            T = zeros(cs.k, 2)
            for i in 1:cs.k
                T[i, 1] = K.normal(g, rng); T[i, 2] = K.normal(g, rng)
            end
            # A_f(X0 + Phi t) quadratic in t  =>  the third central difference vanishes
            i0 = K.uniform_int(rng, 0, cs.k - 1) + 1
            function areas(h)
                T2 = copy(T); T2[i0, 1] += h
                X = from_T(cs, T2)
                return [0.5 * area2(cs.m, X, f) for f in 1:K.n_faces(cs.m)]
            end
            Am2 = areas(-2e-2); Am1 = areas(-1e-2); A0 = areas(0.0); Ap1 = areas(1e-2); Ap2 = areas(2e-2)
            for f in eachindex(A0)
                d3 = Ap2[f] - 2 * Ap1[f] + 2 * Am1[f] - Am2[f]   # ~ 2 h^3 A'''
                e_quad = max(e_quad, abs(d3) / max(1.0, abs(A0[f])))
                n += 1
            end
            # the negative claim: all areas positive yet Theta_max == 0
            X = from_T(cs, T)
            all(a -> a > 0, A0) || continue
            n_pos += 1
            B = K.deploy_basis(cs.c, X)
            sd = K.swept_discs(cs.c, B)
            pairs = K.candidate_pairs(cs.c, sd, Float64(pi), false)
            r = K.exact_theta_max_overlap(cs.c, B, pairs, 1e-9, Float64(pi), 1e-9)
            r.theta_max <= 1e-9 && (n_pos_zero_range += 1)
        end
    end
    record("T5-a", "A_f(X0 + Phi t) is exactly quadratic in t (3rd difference)", n, e_quad, 1e-9)
    println("  positive-orientation samples: $n_pos; of those Theta_max == 0 on $n_pos_zero_range")
    @test (n_pos, n_pos_zero_range) == (70, 51)   # C++
    @test e_quad <= 1e-9
    @test n_pos_zero_range > 0  # T5.1's negative claim reproduced
end

# C++: 0.001 s
@testset "T5.3 closed forms for the split-edge separation harmonic" begin
    rng = K.MT19937(59)
    e_p = 0.0; e_q = 0.0; e_r = 0.0; e_pq = 0.0; e_tau = 0.0; n = 0; n_fn = 0
    for cs in corpus()
        K.n_split(cs.c) == 0 && continue
        for rep in 0:24
            X = sample(cs, rng, rep == 0 ? 0.0 : 0.2)
            P = potential(cs, X)
            B = my_basis(cs, X, P)
            for e in cs.c.split_edges
                ed = cs.m.edges[e]
                f = cs.m.half_edges[ed.he[1]].face; g = cs.m.half_edges[ed.he[2]].face
                a = ed.key.a; b = ed.key.b
                ap = K.prime_vertex(cs.c, f, a); bp = K.prime_vertex(cs.c, f, b)
                aq = K.prime_vertex(cs.c, g, a)
                Cab = row(B.C, bp) - row(B.C, ap); Sab = row(B.S, bp) - row(B.S, ap)
                Caw = row(B.C, aq) - row(B.C, ap); Saw = row(B.S, aq) - row(B.S, ap)
                h = my_orient(Cab, Sab, Caw, Saw)       # generic (T3.2)
                d = X[b] - X[a]; du = P.u[g] - P.u[f]
                sf = Float64(cs.m.sigma[f])
                p = -sf * det2(d, du); q = -p; r = dot(d, du)   # (T5.2)
                sc = max(1e-9, abs(p) + abs(r))
                e_p = max(e_p, abs(h.p - p) / sc)
                e_q = max(e_q, abs(h.q - q) / sc)
                e_r = max(e_r, abs(h.r - r) / sc)
                e_pq = max(e_pq, abs(h.p + h.q) / sc)
                # r = det(d, J du) as well
                e_r = max(e_r, abs(det2(d, Jm(du)) - r) / sc)
                # (T5.3)/(T5.4): the only non-zero root is tau* = -r/p
                if abs(p) > 1e-8 * sc
                    th = 2 * atan(-r / p)
                    e_tau = max(e_tau, abs(heval(h, th)) / sc)
                    (r > 0 && p < 0) && (n_fn += 1)   # Eq.(9) certifies r>0 yet a second root at tau*>0
                end
                n += 1
            end
        end
    end
    record("T5-b", "(T5.2) p = -sigma_f det(d, du)", n, e_p, 1e-11)
    record("T5-c", "(T5.2) q = -p  (so p + q = 0 identically)", n, max(e_q, e_pq), 1e-11)
    record("T5-d", "(T5.2) r = <d, du> = det(d, J du)", n, e_r, 1e-11)
    record("T5-e", "(T5.4) h(2 arctan(-r/p)) = 0", n, e_tau, 1e-9)
    println("  Eq.(9) false-negative ALGEBRAIC precondition (r>0 and p<0): $n_fn / $n split-edge samples  (C++: 871 / 4300)")
    @test n == 4300
    @test e_p <= 1e-11; @test max(e_q, e_pq) <= 1e-11; @test e_r <= 1e-11
    @test e_tau <= 1e-9
    @test n_fn > 0
end

# T5.2b -- is the "conservative" description an INNER approximation of U(eps)?
# C++: 0.709 s
@testset "T5.2b conservative region: inner approximation, and the missing theta=0 atom" begin
    rng = K.MT19937(61)
    n_tested = 0; n_no_root_but_zero_range = 0; n_no_root = 0; n_range_below = 0
    eps = 0.02; T = tan(eps / 2)
    for cs in corpus()
        cs.k == 0 && continue
        g = K.NormalDist(0.0, 0.3)
        for rep in 0:59
            Tm = zeros(cs.k, 2)
            if rep > 0
                for i in 1:cs.k
                    Tm[i, 1] = K.normal(g, rng); Tm[i, 2] = K.normal(g, rng)
                end
            end
            X = from_T(cs, Tm)
            all_faces_positive(cs.m, X) || continue
            B = K.deploy_basis(cs.c, X)
            sd = K.swept_discs(cs.c, B)
            pairs = K.candidate_pairs(cs.c, sd, Float64(pi), false)
            # "no orientation-harmonic root in (0, eps)" over the full candidate list,
            # with NO interval test and NO overlap test -- exactly T5.2b.
            root_in = false
            for (f, gg) in pairs
                for pass in 0:1
                    root_in && break
                    fv = pass == 1 ? gg : f; fe = pass == 1 ? f : gg
                    for pv in cs.c.prime_faces[fv]
                        pe = cs.c.prime_faces[fe]
                        for i in eachindex(pe)
                            root_in && break
                            h = K.orient_harmonic(B, pe[i], pe[mod1(i + 1, length(pe))], pv)
                            hscale(h) < 1e-13 && continue
                            for th in my_roots(h, 1e-9, eps, 1e-14 * hscale(h))
                                if th > 1e-9 && th < eps
                                    root_in = true; break
                                end
                            end
                        end
                        root_in && break
                    end
                end
                root_in && break
            end
            n_tested += 1
            root_in || (n_no_root += 1)
            r = K.exact_theta_max_overlap(cs.c, B, pairs, 1e-9, Float64(pi), 1e-9)
            r.theta_max < eps && (n_range_below += 1)
            (!root_in && r.theta_max < eps) && (n_no_root_but_zero_range += 1)
        end
    end
    @printf("  eps=%.3f (T=%.4f): %d positively-oriented samples; %d with no root in (0,eps);\n  %d with Theta_max < eps; %d satisfy BOTH (no root, yet range < eps)\n",
            eps, T, n_tested, n_no_root, n_range_below, n_no_root_but_zero_range)
    @test (n_tested, n_no_root, n_range_below, n_no_root_but_zero_range) == (343, 3, 261, 0)   # C++
    @test n_no_root > 0
    println("  => T5.2b is an inner approximation only after adding the theta=0+ embeddedness atom: ",
            n_no_root_but_zero_range > 0 ? "CONFIRMED GAP" : "no counterexample found")
end

# =============================================================================
# T6.2 -- exact root gradients by implicit differentiation
# =============================================================================
struct Triple
    fv::Int; fe::Int; pv::Int; pa::Int; pb::Int
end
function harm_at(cs::Case, T::Matrix{Float64}, t::Triple)
    X = from_T(cs, T)
    B = K.deploy_basis(cs.c, X)
    return K.orient_harmonic(B, t.pa, t.pb, t.pv)
end
# the root of h nearest `ref`, or NaN
function root_near(h::K.Harmonic, ref::Float64)
    sc = max(1e-12, hscale(h))
    rs = my_roots(h, -Float64(pi), Float64(pi), 1e-14 * sc)
    best = NaN; bd = 1e300
    for r in rs
        if abs(r - ref) < bd
            bd = abs(r - ref); best = r
        end
    end
    return best
end

# C++: 0.029 s
@testset "T6.2 implicit-differentiation gradient vs central differences" begin
    rng = K.MT19937(67)
    worst_rel = 0.0; n = 0; graphs = 0
    for cs in corpus()
        cs.k < 2 && continue
        graphs += 1
        T0 = zeros(cs.k, 2)
        gg = K.NormalDist(0.0, 0.05)
        for i in 1:cs.k
            T0[i, 1] = K.normal(gg, rng); T0[i, 2] = K.normal(gg, rng)
        end
        X = from_T(cs, T0)
        B0 = K.deploy_basis(cs.c, X)
        F = K.n_faces(cs.m)
        tries = 0
        for attempt in 1:9000
            tries >= 45 && break
            fv = K.uniform_int(rng, 0, F - 1) + 1
            fe = K.uniform_int(rng, 0, F - 1) + 1
            fv == fe && continue
            PV = cs.c.prime_faces[fv]; PE = cs.c.prime_faces[fe]
            pv = PV[K.uniform_int(rng, 0, length(PV) - 1) + 1]
            i0 = K.uniform_int(rng, 0, length(PE) - 1) + 1
            tr = Triple(fv, fe, pv, PE[i0], PE[mod1(i0 + 1, length(PE))])
            h = K.orient_harmonic(B0, tr.pa, tr.pb, tr.pv)
            sc = hscale(h)
            sc < 1e-6 && continue
            rs = my_roots(h, 0.05, Float64(pi) - 0.05, 1e-14 * sc)
            isempty(rs) && continue
            th = rs[1]
            denom = -h.q * sin(th) + h.r * cos(th)   # h'(theta*)
            abs(denom) < 0.05 * sc && continue        # not a simple root, T6.H.1
            tries += 1
            # gradients with respect to a few random design coordinates
            for rep in 1:4
                i = K.uniform_int(rng, 0, cs.k - 1) + 1
                col = K.uniform_int(rng, 0, 1) + 1
                # dp/dt_i, dq/dt_i, dr/dt_i -- exact, because p,q,r are quadratic in t and a
                # central difference is exact on a quadratic.
                hq = 1e-3
                Tp = copy(T0); Tmn = copy(T0)
                Tp[i, col] += hq; Tmn[i, col] -= hq
                hp = harm_at(cs, Tp, tr); hm = harm_at(cs, Tmn, tr)
                dp = (hp.p - hm.p) / (2 * hq)
                dq = (hp.q - hm.q) / (2 * hq)
                dr = (hp.r - hm.r) / (2 * hq)
                grad = -(dp + cos(th) * dq + sin(th) * dr) / denom   # (T6.2)
                # central difference of the ROOT itself, recomputed from scratch
                D = 1e-5
                Ta = copy(T0); Tb = copy(T0)
                Ta[i, col] += D; Tb[i, col] -= D
                ra = root_near(harm_at(cs, Ta, tr), th)
                rb = root_near(harm_at(cs, Tb, tr), th)
                (isfinite(ra) && isfinite(rb)) || continue
                fd = (ra - rb) / (2 * D)
                rel = abs(grad - fd) / max(1e-3, abs(fd))
                worst_rel = max(worst_rel, rel)
                n += 1
            end
        end
    end
    record("T6-a", "(T6.2) dtheta/dt_i vs central difference (relative)", n, worst_rel, 1e-6)
    @test n >= 1000
    @test worst_rel <= 1e-6
end

# =============================================================================
# T7 -- rank and periodic corrections, recomputed independently
# =============================================================================
# JacobiSVD-style rank: singular values above 1e-10 * s_max * max(rows, cols)
function rank_of(A::AbstractMatrix)
    (size(A, 1) == 0 || size(A, 2) == 0) && return 0
    s = svdvals(A)
    t = 1e-10 * s[1] * max(size(A)...)
    return count(x -> x > t, s)
end

# C++: 0.002 s
@testset "T7 (i)-(iii) L = R D, rank(L) = H - dim Z, out-harmonic Z" begin
    e_RD = 0.0; e_harm = 0.0; n = 0; n_rank_ok = 0; n_rank_tot = 0
    for cs in corpus()
        H = K.n_interior_holes(cs.hs)
        H == 0 && continue
        N = K.n_vertices(cs.m); EH = K.n_hinge(cs.c)
        R = zeros(H, EH); D = zeros(EH, N)
        eidx = Dict{Int,Int}()
        for (i, e) in enumerate(cs.c.hinge_edges)
            eidx[e] = i
            D[i, cs.c.hinge_dir[e].dst] += 1.0
            D[i, cs.c.hinge_dir[e].src] -= 1.0
        end
        for r in 1:H, e in cs.hs.all[cs.hs.interior_indices[r]].edges
            haskey(eidx, e) && (R[r, eidx[e]] = 1.0)
        end
        L = cs.sys.A[1:H, :]
        e_RD = max(e_RD, maximum(abs, L - R * D))
        n += 1
        # rank(L) == H - dim Z with Z = left null space
        F = svd(L; full = true)
        tol = 1e-10 * F.S[1] * max(size(L)...)
        rk = count(x -> x > tol, F.S)
        dimZ = H - rk
        n_rank_tot += 1
        rk == H - dimZ && (n_rank_ok += 1)
        # out-harmonic check (T7.1) on each left-null vector
        for z in rk+1:H
            y = F.U[:, z]
            gv = zeros(N); has = falses(N)
            for r in 1:H, v in cs.hs.all[cs.hs.interior_indices[r]].vertices
                gv[v] = y[r]; has[v] = true
            end
            for v in 1:N
                cs.m.vertex_is_boundary[v] && continue
                rhs = 0.0; kin = 0
                for e in cs.m.vertex_edges[v]
                    cs.c.edge_type[e] != K.Hinge && continue
                    hd = cs.c.hinge_dir[e]
                    hd.dst == v && (kin += 1)
                    hd.src == v && (rhs += has[hd.dst] ? gv[hd.dst] : 0.0)
                end
                lhs = kin * (has[v] ? gv[v] : 0.0)
                e_harm = max(e_harm, abs(lhs - rhs))
            end
        end
    end
    record("T7-a", "(i) L == R D exactly", n, e_RD, 0.0)
    record("T7-b", "(iii) every left-null vector of L is out-harmonic", n, e_harm, 1e-9)
    println("  rank(L) == H - dim Z on $n_rank_ok / $n_rank_tot cases (dim Z = 0 on every bounded patch)")
    @test n_rank_tot == 16
    @test e_RD == 0.0
    @test e_harm <= 1e-9
    @test n_rank_ok == n_rank_tot
end

# C++: 0.002 s
@testset "T7 (iv) boundary-free patches: 1^T L = 0 and rank(L) = H - 1" begin
    cases = 0
    cpp_H = [16, 24, 16]
    for (mk, j) in enumerate(DERIV_INPUTS["tori"])
        m = fixture_mesh_raw(j["mesh"])   # torus_* + relaxation(mt19937(9), 4, 300, 60)
        c = K.make_cut(m)
        hs = K.holes_partition(c)
        sys = K.assemble_system(c, hs, m.X, K.None)
        H = K.n_interior_holes(hs)
        H == 0 && continue
        L = sys.A[1:H, :]
        rowsum = L' * ones(H)
        F = svd(L; full = true)
        tol = 1e-10 * F.S[1] * max(size(L)...)
        rk = count(x -> x > tol, F.S)
        @printf("  torus case %d: H=%d rank(L)=%d  ||1^T L||_inf=%.3e  (expect rank = H-1)\n", mk - 1, H, rk, maximum(abs, rowsum))
        @test H == cpp_H[mk]
        @test maximum(abs, rowsum) < 1e-12
        @test rk == H - 1
        # dim Z = 1 with y = 1 (T7 iv), and y = 1 is out-harmonic (T7.1) at every vertex
        y = F.U[:, H]
        spread = maximum(abs, y .- sum(y) / length(y))
        harm = 0.0
        gv = zeros(K.n_vertices(m))
        for r in 1:H, v in hs.all[hs.interior_indices[r]].vertices
            gv[v] = y[r]
        end
        for v in 1:K.n_vertices(m)
            rhs = 0.0; kin = 0
            for e in m.vertex_edges[v]
                c.edge_type[e] != K.Hinge && continue
                c.hinge_dir[e].dst == v && (kin += 1)
                c.hinge_dir[e].src == v && (rhs += gv[c.hinge_dir[e].dst])
            end
            harm = max(harm, abs(kin * gv[v] - rhs))
        end
        @printf("      left-null vector is constant to %.3e; out-harmonic residual %.3e\n", spread, harm)
        @test spread < 1e-9
        @test harm < 1e-9
        cases += 1
    end
    @test cases == 3
end

# C++: 0.0001 s
@testset "T7 (vi) H = |E_hinge| - |F| + c(Gamma)" begin
    ok = 0; tot = 0
    for cs in corpus()
        g = K.build_hinge_graph(cs.c)
        predicted = K.n_hinge(cs.c) - K.n_faces(cs.m) + g.components
        tot += 1
        if predicted == K.n_interior_holes(cs.hs)
            ok += 1
        else
            println("  H mismatch on $(cs.name): H=$(K.n_interior_holes(cs.hs)) predicted=$predicted c(Gamma)=$(g.components)")
        end
    end
    println("  H = |E_h| - |F| + c(Gamma) on $ok / $tot corpus cases")
    @test ok == tot
    @test tot == 16
end

# =============================================================================
# T1 Step 6, attacked directly:  closure of delta on a cycle basis of Gamma
# must have the SAME row space as L (the Eq. (2) rows from holes_partition).
# =============================================================================
struct Step6
    n_cycles::Int; H::Int; rank_G::Int; rank_L::Int; rank_joint::Int
    same_row_space::Bool
end

function step6_check(m::K.Mesh, c::K.CutStructure, hs::K.HoleSet, sys::K.LinearSystem)
    F = K.n_faces(m); N = K.n_vertices(m)
    pot = [zeros(N) for _ in 1:F]
    seen = falses(F)
    adj = K._hinge_adjacency(c)
    tree_edge = falses(K.n_edges(m))
    for s in 1:F
        seen[s] && continue
        seen[s] = true
        q = Int[s]
        while !isempty(q)
            f = popfirst!(q)
            for (e, g) in adj[f]
                if !seen[g]
                    seen[g] = true
                    tree_edge[e] = true
                    pot[g] = copy(pot[f])
                    pot[g][c.hinge_dir[e].src] += m.sigma[g]
                    push!(q, g)
                end
            end
        end
    end
    rows = Vector{Float64}[]
    for e in c.hinge_edges
        tree_edge[e] && continue
        ed = m.edges[e]
        f = m.half_edges[ed.he[1]].face; g = m.half_edges[ed.he[2]].face
        r = pot[g] - pot[f]
        r[c.hinge_dir[e].src] -= m.sigma[g]
        push!(rows, r)
    end
    n_cycles = length(rows)
    H = K.n_interior_holes(hs)
    G = zeros(n_cycles, N)
    for i in 1:n_cycles
        G[i, :] = rows[i]
    end
    L = sys.A[1:H, :]
    Joint = vcat(G, L)
    rG = rank_of(G); rL = rank_of(L); rJ = rank_of(Joint)
    return Step6(n_cycles, H, rG, rL, rJ, rG == rL && rL == rJ)
end

# C++: 0.005 s
@testset "T1 Step 6: Gamma cycle closure and Eq. (2) span the same rows (corpus)" begin
    ok = 0; tot = 0
    for cs in corpus()
        s = step6_check(cs.m, cs.c, cs.hs, cs.sys)
        tot += 1
        if s.same_row_space && s.n_cycles == s.H
            ok += 1
        else
            @printf("  Step-6 mismatch on %-16s b1(Gamma)=%d H=%d rank G=%d L=%d joint=%d\n",
                    cs.name, s.n_cycles, s.H, s.rank_G, s.rank_L, s.rank_joint)
        end
    end
    println("  Step 6 holds on $ok / $tot corpus cases (auto orientation, c(Gamma)=1)")
    @test ok == tot
    @test tot == 16
end

# =============================================================================
# COUNTEREXAMPLE SEARCH: random sigma, including connectivity-violating ones and
# ones that create split-cut cycles.  Attacks F11, T1 Step 6, T7 and Def 4.2.
# The 400 (graph, sigma) pairs are drawn from the Julia generators (bit-exact for the
# random graphs; the tilings are RNG-free), so every tally equals the C++ one.
# =============================================================================
# C++: 0.067 s
@testset "Counterexample search: random sigma, split-cut cycles, disconnected Gamma" begin
    rng = K.MT19937(97)
    n_total = 0; n_split_cycle = 0; n_gamma_disc = 0; n_mprime_disc = 0
    n_partition_fail = 0; n_H_fail = 0; n_step6_fail = 0; n_seed_vs_part_fail = 0
    n_cycle_and_gamma_connected = 0; n_gamma_disc_no_cycle = 0; n_nondeployable_X0 = 0
    n_pot_matches_deploy = 0; n_pot_tested = 0
    worst_pot = 0.0
    gens = [r -> K.tiling_squares(K.disk(K.Vec2(0, 0), 2.2)),
            r -> K.tiling_hexagons(K.disk(K.Vec2(0, 0), 2.4)),
            r -> K.tiling_kagome(K.disk(K.Vec2(0, 0), 2.0)),
            r -> K.tiling_snub_square(K.disk(K.Vec2(0, 0), 2.0)),
            r -> K.largest_component(K.delaunay_of_random_points(30, 8.0, r)),
            r -> K.largest_component(K.voronoi_of_random_points(28, 8.0, r))]
    for trial in 0:399
        m = gens[trial % length(gens) + 1](rng)
        K.n_faces(m) < 4 && continue
        m.sigma = random_sigma(m, rng)   # bernoulli(0.5) per face, true -> +1
        c = try
            K.make_cut(m)
        catch
            continue
        end
        K.n_hinge(c) == 0 && continue
        n_total += 1
        part = K.holes_partition(c)
        seed = K.holes_seed_growing(c)
        K.same_hole_sets(seed, part)[1] || (n_seed_vs_part_fail += 1)
        K.holes_partition_edges(c, part)[1] || (n_partition_fail += 1)
        forest = K.split_subgraph_is_forest(c)[1]
        g = K.build_hinge_graph(c)
        comps_mprime = K.count_components(c)
        forest || (n_split_cycle += 1)
        g.components > 1 && (n_gamma_disc += 1)
        comps_mprime > 1 && (n_mprime_disc += 1)
        (!forest && g.components == 1) && (n_cycle_and_gamma_connected += 1)
        (forest && g.components > 1) && (n_gamma_disc_no_cycle += 1)
        predicted = K.n_hinge(c) - K.n_faces(m) + g.components
        predicted != K.n_interior_holes(part) && (n_H_fail += 1)
        sys = K.assemble_system(c, part, m.X, K.Fixed)
        s6 = step6_check(m, c, part, sys)
        s6ok = s6.same_row_space && s6.n_cycles == s6.H
        if !s6ok
            n_step6_fail += 1
            @printf("  Step-6 FAILS: b1(Gamma)=%d H=%d rank G=%d L=%d joint=%d split-cycle=%d c(Gamma)=%d c(M')=%d\n",
                    s6.n_cycles, s6.H, s6.rank_G, s6.rank_L, s6.rank_joint, Int(!forest), g.components, comps_mprime)
        end
        # T1.H.1: with c(Gamma) > 1 the closed form (T1.5) still holds, per component,
        # PROVIDED X lies in the shape space.  Test on the projected X0, not on the raw X.
        rp = K.solve_system(sys, m.X)
        if rp.projection_ok
            X0 = K.matrix_to_points(rp.X0)
            F = K.n_faces(m)
            u = fill(K.Vec2(0, 0), F)
            seen = falses(F)
            adj = K._hinge_adjacency(c)
            # deliberately a DIFFERENT spanning forest from deploy()'s: depth-first, with
            # the neighbours visited in reverse order.
            for sd in 1:F
                seen[sd] && continue
                seen[sd] = true
                q = Int[sd]
                while !isempty(q)
                    f = pop!(q)
                    for (e, gg) in Iterators.reverse(adj[f])
                        seen[gg] && continue
                        seen[gg] = true
                        u[gg] = u[f] + Float64(m.sigma[gg]) * X0[c.hinge_dir[e].src]
                        push!(q, gg)
                    end
                end
            end
            th = 0.8; cc = cos(th / 2); ss = sin(th / 2)
            D = K.deploy(c, X0, th, 1)
            w = 0.0
            for f in 1:F
                vs = m.faces[f]
                for kk in eachindex(vs)
                    pv = c.prime_faces[f][kk]
                    xv = X0[vs[kk]]
                    y = cc * xv + ss * Jm(2.0 * u[f] - Float64(m.sigma[f]) * xv)
                    w = max(w, norm(y - D.Y[pv]))
                end
            end
            n_pot_tested += 1
            worst_pot = max(worst_pot, w)
            w < 1e-9 && (n_pot_matches_deploy += 1)
            D.max_mismatch > 1e-9 && (n_nondeployable_X0 += 1)
        end
    end
    println("\n  --- random-sigma counterexample search ($n_total usable (graph,sigma) pairs) ---")
    println("  split-cut cycle present                : $n_split_cycle")
    println("  Gamma disconnected  (c(Gamma) > 1)     : $n_gamma_disc")
    println("  M' disconnected                        : $n_mprime_disc")
    println("  split cycle BUT Gamma connected        : $n_cycle_and_gamma_connected   <- refutes 'cycle <=> disconnected'")
    println("  Gamma disconnected BUT split forest    : $n_gamma_disc_no_cycle")
    println("  Alg.1 (seed) != partition formulation  : $n_seed_vs_part_fail")
    println("  preimages do NOT partition E_h u E_s   : $n_partition_fail   (F11)")
    println("  H != |E_h| - |F| + c(Gamma)            : $n_H_fail")
    println("  T1 Step 6 row-space mismatch           : $n_step6_fail")
    @printf("  (T1.5) with a DIFFERENT spanning forest reproduces deploy() on X0 : %d/%d, worst %.3e\n",
            n_pot_matches_deploy, n_pot_tested, worst_pot)
    println("  X0 with deploy() max_mismatch > 1e-9   : $n_nondeployable_X0\n")
    record("CE-a", "F11 partition of E_hinge u E_split under RANDOM sigma", n_total, n_partition_fail, 0.0)
    record("CE-b", "H = |E_h|-|F|+c(Gamma) under RANDOM sigma", n_total, n_H_fail, 0.0)
    record("CE-c", "T1 Step 6 row space under RANDOM sigma", n_total, n_step6_fail, 0.0)
    record("CE-d", "(T1.5) tree-independent on X0, random sigma, c(Gamma)>1", n_pot_tested, worst_pot, 1e-9)
    @test n_total > 100
    @test worst_pot < 1e-9
    @test n_pot_tested > 100
    # the C++ tallies (the whole population is reproduced bit-exactly)
    @test n_total == 398
    @test n_split_cycle == 210
    @test n_gamma_disc == 334
    @test n_mprime_disc == 334
    @test n_cycle_and_gamma_connected == 0
    @test n_gamma_disc_no_cycle == 124
    @test n_seed_vs_part_fail == 0
    @test n_partition_fail == 0
    @test n_H_fail == 0
    @test n_step6_fail == 0
    @test n_pot_matches_deploy == 398
    @test n_pot_tested == 398
    @test n_nondeployable_X0 == 0
end

# =============================================================================
# T4.5a -- is the broad phase actually sound?  Two versions:
#   (i)  the moving-centroid test of (T4.2)/(T4.3), which is what core.md derives;
#   (ii) contact.jl's `use_static` option, which uses the FLAT centroids.
# =============================================================================
# C++: 0.396 s
@testset "T4.5a broad phase: moving-centroid pruning is sound, flat-centroid is not" begin
    n_moving_bad = 0; n_static_bad = 0; cases = 0
    n_flat_unsound = 0
    worst_flat_gap = 0.0
    worst_moving = 0.0; worst_static = 0.0
    rng = K.MT19937(131)
    tot_all = 0; tot_mov = 0; tot_st = 0
    for cs in corpus(), rep in 0:5
        Xs = sample(cs, rng, rep == 0 ? 0.0 : 0.2)
        B = K.deploy_basis(cs.c, Xs)
        sd = K.swept_discs(cs.c, B)
        all_ = K.candidate_pairs(cs.c, sd, Float64(pi), false)
        moving = K.candidate_pairs(cs.c, sd, Float64(pi), true, false)
        stat = K.candidate_pairs(cs.c, sd, Float64(pi), true, true)
        tot_all += length(all_); tot_mov += length(moving); tot_st += length(stat)
        # direct soundness audit: pairs the FLAT test discards although the sound moving
        # test keeps them.  Each is a pair whose faces can actually approach.
        keep = Set(stat)
        for pr in moving
            pr in keep && continue
            n_flat_unsound += 1
            dmin = K.min_center_distance(sd, pr[1], pr[2], Float64(pi))
            c0 = norm(K.Vec2(sd.gc[pr[1], 1], sd.gc[pr[1], 2]) - K.Vec2(sd.gc[pr[2], 1], sd.gc[pr[2], 2]))
            worst_flat_gap = max(worst_flat_gap, c0 - dmin)
        end
        t_all = K.exact_theta_max_overlap(cs.c, B, all_, 1e-9, Float64(pi), 1e-9).theta_max
        t_mov = K.exact_theta_max_overlap(cs.c, B, moving, 1e-9, Float64(pi), 1e-9).theta_max
        t_st = K.exact_theta_max_overlap(cs.c, B, stat, 1e-9, Float64(pi), 1e-9).theta_max
        worst_moving = max(worst_moving, abs(t_mov - t_all))
        worst_static = max(worst_static, abs(t_st - t_all))
        abs(t_mov - t_all) > 1e-6 && (n_moving_bad += 1)
        if abs(t_st - t_all) > 1e-6
            n_static_bad += 1
            @printf("  flat-centroid pruning WRONG on %-16s: %d/%d pairs kept, Theta_max %.6f vs %.6f\n",
                    cs.name, length(stat), length(all_), t_st, t_all)
        end
        cases += 1
    end
    record("T4-f", "(T4.2)/(T4.3) moving-centroid pruning preserves Theta_max", cases, worst_moving, 1e-6)
    record("T4-g", "flat-centroid pruning (K2c's stated form) preserves Theta_max", cases, worst_static, 1e-6)
    @printf("  pairs kept: all %d, moving-centroid %d (%.1f%%), flat-centroid %d (%.1f%%)\n",
            tot_all, tot_mov, 100.0 * tot_mov / max(1, tot_all), tot_st, 100.0 * tot_st / max(1, tot_all))
    @printf("  flat-centroid test discards %d pairs the SOUND moving test keeps; worst gap %.3f\n", n_flat_unsound, worst_flat_gap)
    println("  moving-centroid pruning wrong on $n_moving_bad/$cases samples; flat-centroid on $n_static_bad/$cases")
    @test (tot_all, tot_mov, tot_st, n_flat_unsound) == (102534, 27463, 23817, 3646)   # C++
    @test isapprox(worst_flat_gap, 6.022; atol = 1e-3)
    @test n_static_bad == 0
    @test cases == 96
    @test n_moving_bad == 0
end

# =============================================================================
# T7 (v) -- boundary-touching (notch) preimages contribute no row of L, and the
# corrected row-sum identity holds while the naive one does not.
# =============================================================================
# C++: 0.0001 s
@testset "T7 (v) notches: corrected row-sum identity, naive one fails" begin
    ok_corrected = 0; ok_naive = 0; tot = 0; with_notch = 0
    for cs in corpus()
        rs = K.check_row_sum(cs.c, cs.hs, cs.sys)
        tot += 1
        rs.restricted_degree_identity && (ok_corrected += 1)
        rs.degree_identity && (ok_naive += 1)
        rs.no_notch_hinges || (with_notch += 1)
    end
    println("  corrected identity $ok_corrected/$tot, naive identity $ok_naive/$tot, cases with notch-owned hinges $with_notch")
    @test ok_corrected == tot
    @test with_notch > 0
    @test ok_naive < tot   # the naive identity is false on patches with a boundary
    @test (ok_corrected, ok_naive, with_notch) == (16, 0, 16)
end

# =============================================================================
# ROUND 2 -- checks of the NEW material only
# =============================================================================

# Lemma T4.5b' : in the FACE'S OWN frame the swept radius is exact, and in fact
# the distance to the moving centroid is CONSTANT in theta.
# C++: 0.013 s
@testset "R2 T4.5b' exact swept radius in the face frame (randomized)" begin
    rng = K.MT19937(90210)
    n = 0
    e_x = 0.0; e_chi = 0.0; e_orth = 0.0; e_norm = 0.0; e_smax = 0.0; e_const = 0.0; e_max = 0.0
    e_lib = 0.0; n_lib = 0
    for cs in corpus(), rep in 0:7
        X = sample(cs, rng, rep == 0 ? 0.0 : 0.25)
        P = potential(cs, X)
        B = my_basis(cs, X, P)
        DB = K.deploy_basis(cs.c, X)
        sd = K.swept_discs(cs.c, DB)
        for f in 1:K.n_faces(cs.m)
            vs = cs.m.faces[f]
            sg = Float64(cs.m.sigma[f])
            # MY flat centroid, and MY moving-centroid pair (gc, gs) from (T1.5).
            xbar = sum(X[v] for v in vs) / length(vs)
            gc = xbar
            gs = Jm(2.0 * P.u[f] - sg * xbar)
            # cross-check against the library's centroid pair (not used in the algebra)
            e_lib = max(e_lib, norm(gc - K.Vec2(sd.gc[f, 1], sd.gc[f, 2])))
            e_lib = max(e_lib, norm(gs - K.Vec2(sd.gs[f, 1], sd.gs[f, 2])))
            n_lib += 1
            for kk in eachindex(vs)
                pv = cs.c.prime_faces[f][kk]
                Cu = row(B.C, pv); Su = row(B.S, pv)
                x = Cu - gc; chi = Su - gs
                d = X[vs[kk]] - xbar
                sc = max(1.0, norm(d))
                # (T4.4') the two identities -- the 2u_f term must cancel
                e_x = max(e_x, norm(x - d) / sc)
                e_chi = max(e_chi, norm(chi + sg * Jm(d)) / sc)
                # orthogonal, equal norm  =>  [x|chi] is a similarity
                e_orth = max(e_orth, abs(dot(x, chi)) / (sc * sc))
                e_norm = max(e_norm, abs(norm(x) - norm(chi)) / sc)
                smax = svdvals([x chi])[1]
                rmax = max(norm(x), norm(chi))
                e_smax = max(e_smax, abs(smax - rmax) / sc)
                # the trajectory radius about the MOVING centroid is constant = |x_u - xbar_f|
                lo = 1e300; hi = 0.0
                for i in 0:120
                    th = pi * i / 120.0
                    c = cos(th / 2); s = sin(th / 2)
                    rr = norm((c * Cu + s * Su) - (c * gc + s * gs))
                    lo = min(lo, rr); hi = max(hi, rr)
                end
                e_const = max(e_const, (hi - lo) / sc)
                e_max = max(e_max, abs(hi - norm(d)) / sc)
                n += 1
            end
        end
    end
    println("  R2 T4.5b': $n copies over $(length(corpus())) cases x 8 shape samples")
    record("R2-a", "C_u - gc_f == x_u - xbar_f  (2u_f cancels)", n, e_x, 1e-12)
    record("R2-b", "S_u - gs_f == -sigma_f J (x_u - xbar_f)", n, e_chi, 1e-12)
    record("R2-c", "<x,chi> = 0 and |x| = |chi|", n, max(e_orth, e_norm), 1e-12)
    record("R2-d", "sigma_max([x|chi]) == max(|x|,|chi|)  (EXACT, no sqrt2)", n, e_smax, 1e-12)
    record("R2-e", "|y_u(th)-gamma_f(th)| is CONSTANT in theta", n, e_const, 1e-12)
    record("R2-f", "that constant == |x_u - xbar_f| = flat radius", n, e_max, 1e-12)
    record("R2-g", "my (gc,gs) == swept_discs() (library cross-check)", n_lib, e_lib, 1e-10)
    @test n >= 1000
    @test n == 20176
    @test e_x <= 1e-12; @test e_chi <= 1e-12; @test e_orth <= 1e-12; @test e_norm <= 1e-12
    @test e_smax <= 1e-12; @test e_const <= 1e-12; @test e_max <= 1e-12; @test e_lib <= 1e-10
end

# (T5.1c) / (T5.1d): "no root of g in the OPEN interval (0,T)" as a quantifier-free
# atom list, WITH and WITHOUT the tau = 0 deflation clause.
# (T5.1c) as printed by the round-1 Checker: correct only when g(0) != 0.
function noroot_printed(A, B, C, T)
    g0 = C; gT = C + B * T + A * T * T
    disc = B * B - 4 * A * C
    both = (disc >= 0) && (A * g0 > 0) && (A * gT > 0) && (A * B < 0) && (-A * B - 2 * A * A * T < 0)
    return (g0 * gT > 0) && !both
end
# (T5.1d): the deflated clause, used when g(0) = 0.
noroot_deflated(A, B, C, T) = !((A * B < 0) && (-A * B - A * A * T < 0))

# C++: 0.377 s (1.8M harmonics per eps; the Julia loop runs the same population)
@testset "R2 tau=0 deflation (T5.1d): spurious roots with and without it" begin
    rng = K.MT19937(777)
    cpp = Dict(0.2 => (1817732, 74540, 9460, 2176, 61201), 0.02 => (1817748, 74556, 9444, 2176, 62310),
               0.001 => (1817705, 74513, 9487, 2176, 62327))
    for eps in (0.2, 0.02, 0.001)
        T = tan(eps / 2)
        n = 0; mm_printed = 0; mm_printed_g0 = 0; mm_defl = 0
        spur_printed = 0; spur_defl = 0; miss_printed = 0; miss_defl = 0
        n_g0 = 0; n_dbl = 0; ambiguous = 0
        for cs in corpus(), rep in 0:3
            X = sample(cs, rng, rep == 0 ? 0.0 : 0.25)
            Bs = K.deploy_basis(cs.c, X)
            sd = K.swept_discs(cs.c, Bs)
            pairs = K.candidate_pairs(cs.c, sd, Float64(pi), false)
            for (f0, g0f) in pairs, pass in 0:1
                fv = pass == 1 ? g0f : f0; fe = pass == 1 ? f0 : g0f
                pe = cs.c.prime_faces[fe]
                for pv in cs.c.prime_faces[fv], i in eachindex(pe)
                    h = K.orient_harmonic(Bs, pe[i], pe[mod1(i + 1, length(pe))], pv)
                    sc = hscale(h)
                    sc < 1e-13 && continue                       # identically zero (T3.H.1)
                    A = h.p - h.q; B = 2 * h.r; C = h.p + h.q
                    # TRUTH: roots in the OPEN interval (0, eps), from the tau-quadratic.
                    rs = my_roots(h, 0.0, Float64(pi), 1e-14 * sc)
                    amb = false; truth = false
                    for th in rs
                        if abs(th) < 1e-9 || abs(th - eps) < 1e-9
                            amb = true; break
                        end
                        (th > 0 && th < eps) && (truth = true)
                    end
                    if amb
                        ambiguous += 1; continue
                    end
                    structural_zero = abs(C) <= 1e-11 * sc
                    structural_zero && (n_g0 += 1)
                    # THIRD structural class: C = 0 AND B = 0 (tau = 0 a DOUBLE root);
                    # ill-posed for both sides, counted apart.
                    if structural_zero && abs(B) <= 1e-11 * sc
                        n_dbl += 1; continue
                    end
                    p_noroot = noroot_printed(A, B, C, T)
                    d_noroot = structural_zero ? noroot_deflated(A, B, C, T) : noroot_printed(A, B, C, T)
                    if p_noroot == truth   # p_noroot means "no root"; truth means "root"
                        mm_printed += 1
                        structural_zero && (mm_printed_g0 += 1)
                        truth ? (miss_printed += 1) : (spur_printed += 1)
                    end
                    if d_noroot == truth
                        mm_defl += 1
                        truth ? (miss_defl += 1) : (spur_defl += 1)
                    end
                    n += 1
                end
            end
        end
        @printf("  eps=%-6.3f  harmonics=%d  (g(0)=0 on %d, %.1f%%)  ambiguous skipped=%d  double-root-at-0 class skipped=%d\n",
                eps, n, n_g0, 100.0 * n_g0 / max(1, n), ambiguous, n_dbl)
        @printf("    (T5.1c) as printed : mismatches %d  (spurious roots %d, missed %d); %d of them have g(0)=0\n",
                mm_printed, spur_printed, miss_printed, mm_printed_g0)
        @printf("    (T5.1c)+(T5.1d)    : mismatches %d  (spurious roots %d, missed %d)\n", mm_defl, spur_defl, miss_defl)
        c = cpp[eps]
        println("    C++: harmonics=$(c[1]) g(0)=0 on $(c[2]) ambiguous=$(c[3]) double=$(c[4]) printed mismatches=$(c[5])")
        record("R2-h", "deflated atom list (T5.1c)+(T5.1d) == direct root finding", n, mm_defl, 0.0)
        # the ambiguous/decided split is tie-sensitive (roots within 1e-9 of 0 or eps decided
        # by the last bits of atan/sqrt): the port differs from the C++ by < 20 of 1.8M
        @test abs(n - c[1]) <= 20
        @test abs(ambiguous - c[3]) <= 20
        @test n_dbl == c[4]
        @test abs(mm_printed - c[5]) <= 20
        @test n >= 1000
        @test mm_defl == 0
        # every failure of the undeflated list is a g(0) = 0 pair, and it is SPURIOUS
        @test mm_printed_g0 == mm_printed
        @test miss_printed == 0
    end
end

# =============================================================================
# ROUND 3 (final confirmation).  Two randomized tests, both of which decide the
# truth WITHOUT the tau chart and WITHOUT the class split.
# =============================================================================

# h(theta) = p + q cos + r sin, rewritten in amplitude/phase form (core.md T3.4 /
# T3.H.4).  This is the ONLY root finder used below.  (The C++ `long double` is the
# 64-bit double on the reference platform; libm hypot/atan2/acos/cos as in the C++.)
struct AmpPhase
    p::Float64; R::Float64; phi::Float64
end
amp_phase(h::K.Harmonic) = AmpPhase(h.p, K.libm_hypot(h.q, h.r), K.libm_atan2(h.r, h.q))
ap_eval(a::AmpPhase, th::Float64) = fma(a.R, K.libm_cos(th - a.phi), a.p)   # p + R cos(th - phi)
# All roots of h in (0, hi), from cos(theta - phi) = -p/R.  Closed form.
function amp_phase_roots(a::AmpPhase, hi::Float64)
    out = Float64[]
    a.R <= 0 && return out
    c = -a.p / a.R
    (c < -1.0 || c > 1.0) && return out
    psi = K.libm_acos(c)
    for k in -2:2, sgn in (-1, 1)
        th = a.phi + sgn * psi + 2.0 * Float64(pi) * k
        (th > 0 && th < hi) && push!(out, th)
    end
    sort!(out)
    return out
end

# Independent truth for "does h CROSS zero somewhere in (0, eps)?" -- a dense scan of
# sign(h) on [th_lo, eps); returns 1 = crosses, 0 = does not cross, -1 = undecidable.
function truth_crosses_dense(h::K.Harmonic, eps::Float64)
    scl = abs(h.p) + abs(h.q) + abs(h.r)
    scl <= 0 && return -1
    noise = 1e-13 * scl      # coefficients arrive with ~1e-16 rel. error
    a = amp_phase(h)
    th_lo = min(1e-5, eps / 10.0)
    for r in amp_phase_roots(a, th_lo)
        r > 1e-12 && return -1                   # a root in the un-scanned sliver: undecidable
    end
    N = 4001
    prev = 0
    for i in 0:N-1
        th = th_lo + (eps - th_lo) * i / (N - 1)
        v = ap_eval(a, th)
        abs(v) <= noise && return -1             # a sample sits in the noise band
        s = v > 0 ? 1 : -1
        (i > 0 && s != prev) && return 1
        prev = s
    end
    return 0
end

# The atom list of core.md (T5.1c)/(T5.1d)/(T5.1e'), transcribed literally.
# three_classes = false reproduces round 2 (class 3 routed to (T5.1d)).
@enum R3Klass R3_GENERIC = 0 R3_SIMPLE0 = 1 R3_DOUBLE0 = 2
function r3_classify(h::K.Harmonic, tol::Float64)
    abs(h.p + h.q) > tol && return R3_GENERIC
    abs(2 * h.r) > tol && return R3_SIMPLE0
    return R3_DOUBLE0
end
function r3_atoms_no_root(h::K.Harmonic, T::Float64, k::R3Klass, three_classes::Bool)
    A = h.p - h.q; B = 2 * h.r; C = h.p + h.q
    if k == R3_DOUBLE0
        three_classes && return true             # (T5.1e')
        k = R3_SIMPLE0                           # round 2: only the single deflation
    end
    k == R3_SIMPLE0 && return !((A * B < 0) && (-A * B - A * A * T < 0))   # (T5.1d)
    g0 = C; gT = C + B * T + A * T * T                                     # (T5.1c)
    (g0 * gT > 0) || return false
    disc = B * B - 4 * A * C                                               # (T5.1b)
    both = (disc >= 0) && (A * g0 > 0) && (A * gT > 0) && (A * B < 0) && (-A * B - 2 * A * A * T < 0)
    return !both
end

function face_is_convex(P::Vector{K.Vec2})
    n = length(P)
    n < 3 && return false
    sign = 0
    for i in 1:n
        u = P[mod1(i + 1, n)] - P[i]
        v = P[mod1(i + 2, n)] - P[mod1(i + 1, n)]
        d = det2(u, v)
        abs(d) < 1e-12 && continue
        s = d > 0 ? 1 : -1
        if sign == 0
            sign = s
        elseif s != sign
            return false
        end
    end
    return sign != 0
end

# C++: 37.6 s -- the largest case.  R3-a/R3-b are synthetic and always run; R3-c (the
# corpus scan, ~1.8M harmonics per eps) runs on the first 8 corpus cases by default and
# on all 16 (with the C++ tallies asserted) under KIRIGAMI_FULL_DERIVATIONS=1.
@testset "R3 Lemma T5.1e and the three-class atom list (randomized, chart-free truth)" begin
    rng = K.MT19937(30903)

    # ---- R3-a: the change of variable (T3.5).  (1 + tau^2) h(theta) = C + B tau + A tau^2.
    let
        worst = 0.0; n = 0
        for i in 1:200000
            p = K.uniform_real(rng, -2.0, 2.0); q = K.uniform_real(rng, -2.0, 2.0); r = K.uniform_real(rng, -2.0, 2.0)
            A = p - q; B = 2 * r; C = p + q
            th = K.uniform_real(rng, 1e-4, Float64(pi) - 1e-4)
            tau = tan(th / 2)
            lhs = (1 + tau * tau) * (p + q * cos(th) + r * sin(th))
            rhs = C + B * tau + A * tau * tau
            sc = max(1.0, abs(lhs))
            worst = max(worst, abs(lhs - rhs) / sc)
            n += 1
        end
        record("R3-a", "(T3.5)  (1+tau^2) h = (p+q) + 2r tau + (p-q) tau^2", n, worst, 1e-12)
        @test worst <= 1e-12
    end

    # ---- R3-b: the class conditions, re-derived.  C = 0 and B = 0  <=>  h = p (1 - cos).
    let
        worst_id = 0.0; worst_sign = 0.0; n = 0; sign_changes = 0; zeros_inside = 0
        for i in 1:20000
            p = K.uniform_real(rng, -2.0, 2.0)
            abs(p) < 1e-3 && continue
            q = -p; r = 0.0             # exactly class 3: C = p+q = 0, B = 2r = 0
            @test abs(p + q) == 0.0
            @test abs(2 * r) == 0.0
            A = p - q                   # A = 2p, so p = A/2 as core.md states
            @test abs(A - 2 * p) <= 1e-15
            prev = 0
            mn = 1e300
            for j in 1:3999
                th = pi * j / 4000.0
                h = p + q * cos(th) + r * sin(th)
                worst_id = max(worst_id, abs(h - p * (1 - cos(th))))
                h == 0.0 && (zeros_inside += 1)
                s = h > 0 ? 1 : (h < 0 ? -1 : 0)
                (prev != 0 && s != prev) && (sign_changes += 1)
                s != 0 && (prev = s)
                mn = min(mn, h / p)     # must stay > 0: h has the sign of p
            end
            worst_sign = max(worst_sign, -min(0.0, mn))
            n += 1
        end
        record("R3-b1", "class 3 (C=B=0) is exactly h(th) = p (1 - cos th)", n, worst_id, 1e-15)
        record("R3-b2", "T5.1e: sign h = sign p on (0,pi), no zero, no sign change", n, worst_sign + sign_changes + zeros_inside, 0.0)
        @test n == 19989   # C++
        @test worst_id <= 1e-15
        @test sign_changes == 0
        @test zeros_inside == 0
    end

    # ---- R3-c: the three-class atom list vs the chart-free crossing truth, on the
    #      real corpus.  Also the round-2 two-class list, which must be WRONG somewhere.
    r3_cases = FULL ? corpus() : corpus()[1:8]
    FULL || println("  R3-c: default subset = first 8 corpus cases (KIRIGAMI_FULL_DERIVATIONS=1 runs all 16; C++ 37.6 s)")
    for eps in (0.2, 0.02)
        T = tan(eps / 2)
        tested = 0; ambiguous = 0; ident_zero = 0; ident_zero_in_class3 = 0
        nclass = zeros(Int, 3)
        mm3 = 0; mm2 = 0; mm3_by = zeros(Int, 3); mm2_by = zeros(Int, 3)
        mm3_conservative = 0
        r2 = K.MT19937(70707)
        for cs in r3_cases, rep in 0:3
            X = sample(cs, r2, rep == 0 ? 0.0 : 0.25)
            B = K.deploy_basis(cs.c, X)
            gs = maximum(norm, X)
            tol = 1e-11 * gs * gs
            F = K.n_faces(cs.m)
            for f in 1:F
                pf = cs.c.prime_faces[f]
                nf = length(pf)
                for kk in 1:nf
                    a = pf[kk]; b = pf[mod1(kk + 1, nf)]
                    for g in 1:F
                        g == f && continue
                        for w in cs.c.prime_faces[g]
                            h = K.orient_harmonic(B, a, b, w)
                            kl = r3_classify(h, tol)
                            if hscale(h) <= tol          # h == 0 identically (T3.H.1)
                                ident_zero += 1
                                kl == R3_DOUBLE0 && (ident_zero_in_class3 += 1)
                                continue
                            end
                            nclass[Int(kl)+1] += 1
                            tr = truth_crosses_dense(h, eps)
                            if tr < 0
                                ambiguous += 1; continue
                            end
                            tested += 1
                            truth_noroot = (tr == 0)
                            a3 = r3_atoms_no_root(h, T, kl, true)
                            a2 = r3_atoms_no_root(h, T, kl, false)
                            if a3 != truth_noroot
                                mm3 += 1; mm3_by[Int(kl)+1] += 1
                                (truth_noroot && !a3) && (mm3_conservative += 1)
                            end
                            if a2 != truth_noroot
                                mm2 += 1; mm2_by[Int(kl)+1] += 1
                            end
                        end
                    end
                end
            end
        end
        @printf("  R3-c eps=%-6.3f decided=%d ambiguous=%d  classes %d/%d/%d\n       identically-zero harmonics excluded=%d (of which class 3 by the numeric rule: %d)\n       three-class list: %d mismatches [%d/%d/%d] (of which conservative: %d)\n       two-class  list: %d mismatches [%d/%d/%d]\n",
                eps, tested, ambiguous, nclass[1], nclass[2], nclass[3], ident_zero, ident_zero_in_class3,
                mm3, mm3_by[1], mm3_by[2], mm3_by[3], mm3_conservative, mm2, mm2_by[1], mm2_by[2], mm2_by[3])
        record("R3-c" * (eps > 0.1 ? "1" : "2"), "3-class atom list == chart-free crossing truth", tested, mm3, 0.0)
        if FULL
            # C++ (both eps): decided 1829226, ambiguous 142, classes 1745368/81500/2500,
            # identically zero 14928 (all class 3), two-class mismatches 393 (all class 3)
            @test tested > 100000
            @test abs(tested - 1829226) <= 20
            @test abs(ambiguous - 142) <= 20
            @test nclass == [1745368, 81500, 2500]
            @test ident_zero == 14928 && ident_zero_in_class3 == 14928
            @test abs(mm2 - 393) <= 5
        else
            @test tested > 100000
        end
        @test mm3 == 0
        @test mm2 > 0                       # the class-3 clause is NECESSARY, not decorative
        @test mm2_by[3] == mm2              # and every round-2 failure is in class 3
    end
end

# C++: 0.091 s
@testset "R3 certificate POS ^ NOOVERLAP(eps/2) ^ NOROOT implies Theta_max >= eps (randomized)" begin
    rng = K.MT19937(31003)
    eps = 0.006; T = tan(eps / 2); th1 = eps / 2
    NSCAN = 48
    n_convex_cases = 0; n_samples = 0; n_cert = 0; n_violation = 0
    n_noroot = 0; n_status_nonconstant = 0
    sat_eps = 1e-9
    function overlaps_at(cs::Case, B::K.DeployBasis, th::Float64)
        Y = K.basis_eval(B, th)
        F = K.n_faces(cs.m)
        P = [face_poly(cs, Y, f) for f in 1:F]
        for f in 1:F, g in f+1:F
            sat_overlap_convex(P[f], P[g], sat_eps) && return true
        end
        return false
    end
    for cs in corpus()
        # SAT is exact only for convex polygons; faces are rigid (T1.D) so convexity is
        # theta-independent and can be decided once, in the flat state.
        all_convex = all(f -> face_is_convex(K.Vec2[cs.X0[v] for v in cs.m.faces[f]]), 1:K.n_faces(cs.m))
        all_convex || continue
        n_convex_cases += 1
        for rep in 0:39
            X = sample(cs, rng, rep == 0 ? 0.0 : 0.12)
            # (i) POS
            all_faces_positive(cs.m, X) || continue
            n_samples += 1
            B = K.deploy_basis(cs.c, X)
            gs = maximum(norm, X)
            tol = 1e-11 * gs * gs
            # (iii) NOROOT, over the COMPLETE ordered (vertex, edge) list, three-class atoms.
            noroot = true
            F = K.n_faces(cs.m)
            for f in 1:F
                noroot || break
                pf = cs.c.prime_faces[f]
                nf = length(pf)
                for kk in 1:nf
                    noroot || break
                    a = pf[kk]; b = pf[mod1(kk + 1, nf)]
                    for g in 1:F
                        noroot || break
                        g == f && continue
                        for w in cs.c.prime_faces[g]
                            h = K.orient_harmonic(B, a, b, w)
                            hscale(h) <= tol && continue              # identically zero (T3.H.1)
                            if !r3_atoms_no_root(h, T, r3_classify(h, tol), true)
                                noroot = false; break
                            end
                        end
                    end
                end
            end
            noroot && (n_noroot += 1)
            # (ii) NOOVERLAP(theta_1), one exact test at the single angle theta_1 = eps/2.
            nooverlap = !overlaps_at(cs, B, th1)
            # Independent verdict: scan the whole open interval (0, eps).
            any_overlap = false; n_ov = 0
            for i in 1:NSCAN-1
                if overlaps_at(cs, B, eps * i / NSCAN)
                    any_overlap = true; n_ov += 1
                end
            end
            # Proposition T5.2b': the certificate must imply no overlap anywhere on (0, eps).
            if nooverlap && noroot
                n_cert += 1
                any_overlap && (n_violation += 1)
            end
            # The clopen step on its own: NOROOT => the overlap status is CONSTANT on (0, eps).
            (noroot && n_ov != 0 && n_ov != NSCAN - 1) && (n_status_nonconstant += 1)
        end
    end
    @printf("  R3-e eps=%.4f: %d all-convex cases, %d positively-oriented samples;\n       NOROOT holds on %d; full certificate holds on %d; VIOLATIONS %d\n       R3-f NOROOT but overlap status NOT constant on (0,eps): %d\n",
            eps, n_convex_cases, n_samples, n_noroot, n_cert, n_violation, n_status_nonconstant)
    record("R3-e", "POS ^ NOOVERLAP(eps/2) ^ NOROOT => no overlap on (0,eps)", n_cert, n_violation, 0.0)
    record("R3-f", "NOROOT => overlap status constant on (0,eps) (clopen step)", n_noroot, n_status_nonconstant, 0.0)
    @test (n_convex_cases, n_samples, n_noroot, n_cert) == (10, 280, 176, 175)   # C++
    @test n_cert > 0                 # the hypothesis must have content
    @test n_violation == 0
    @test n_status_nonconstant == 0
end

# =============================================================================
# ROUND 4 -- Sub-lemma T5.2b'' (core.md, inside the proof of Proposition T5.2b')
# =============================================================================

# Crossing-number test plus a strict clearance margin from every edge.
function strictly_inside(P::Vector{K.Vec2}, z::K.Vec2, margin::Float64)
    n = length(P)
    inside = false
    j = n
    for i in 1:n
        a = P[j]; b = P[i]
        if (b[2] > z[2]) != (a[2] > z[2])
            xx = b[1] + (z[2] - b[2]) * (a[1] - b[1]) / (a[2] - b[2])
            z[1] < xx && (inside = !inside)
        end
        j = i
    end
    inside || return false
    j = n
    for i in 1:n
        d = P[i] - P[j]
        L2 = sq(d)
        if L2 >= 1e-30
            t = dot(z - P[j], d) / L2
            t = max(0.0, min(1.0, t))
            norm(z - (P[j] + t * d)) <= margin && return false
        end
        j = i
    end
    return true
end

# Lemma T4.2 and all of T4/T5 are stated for CLOSED SIMPLE polygons.
function segs_cross(a::K.Vec2, b::K.Vec2, c::K.Vec2, d::K.Vec2)
    o(p, q, r) = (v = det2(q - p, r - p); v > 1e-14 ? 1 : (v < -1e-14 ? -1 : 0))
    return o(a, b, c) * o(a, b, d) < 0 && o(c, d, a) * o(c, d, b) < 0
end
function poly_is_simple(P::Vector{K.Vec2})
    n = length(P)
    n < 3 && return false
    for i in 1:n, j in i+1:n
        (mod1(j + 1, n) == i || mod1(i + 1, n) == j) && continue
        segs_cross(P[i], P[mod1(i + 1, n)], P[j], P[mod1(j + 1, n)]) && return false
    end
    return true
end

function dist_point_seg(z::K.Vec2, a::K.Vec2, b::K.Vec2)
    d = b - a
    L2 = sq(d)
    L2 < 1e-30 && return norm(z - a)
    t = dot(z - a, d) / L2
    t = max(0.0, min(1.0, t))
    return norm(z - (a + t * d))
end

# theta-independent radius of Sub-lemma T5.2b'': half the distance from the hinge
# vertex to the edges of the face NOT incident to it.
function sector_radius(vs::Vector{Int}, X::Vector{K.Vec2}, v::Int)
    n = length(vs)
    d = 1e300
    for i in 1:n
        a = vs[i]; b = vs[mod1(i + 1, n)]
        (a == v || b == v) && continue                 # incident edge: skipped
        d = min(d, dist_point_seg(X[v], X[a], X[b]))
    end
    return d
end

# CCW interior angle of face at vertex v, in [0, 2pi) -- reflex corners included.
function interior_angle(vs::Vector{Int}, X::Vector{K.Vec2}, v::Int)
    n = length(vs)
    idx = findlast(==(v), vs)
    idx === nothing && return -1.0
    a = X[vs[mod1(idx - 1, n)]] - X[v]
    b = X[vs[mod1(idx + 1, n)]] - X[v]
    t = atan(det2(b, a), dot(b, a))
    t < 0 && (t += 2pi)
    return t
end

# (idx of v in vs, idx of the FAR endpoint of the far-side edge at v), the C++ far_of
function far_of(vs::Vector{Int}, v::Int, dst::Int)
    n = length(vs)
    idx = findlast(==(v), vs)
    idx === nothing && return (0, 0)
    pv = vs[mod1(idx - 1, n)]; nx = vs[mod1(idx + 1, n)]
    far = (nx == dst) ? pv : nx
    fi = findlast(==(far), vs)
    return (idx, fi === nothing ? 0 : fi)
end

# C++: 0.093 s
@testset "R4 Sub-lemma T5.2b'' Case A: local overlap at the hinge <=> theta > beta_e" begin
    rng = K.MT19937(40041)
    delta = 0.05          # stay clear of the transition angle beta_e
    NPHI = 720; NRAD = 3
    n_samples = 0; n_local_tests = 0
    n_above = 0; n_below = 0; n_fail_above = 0; n_fail_below = 0
    n_pairs_two_shared = 0; n_face_pairs = 0; n_nonsimple = 0
    n_beta_pairs = 0; n_beta_degenerate = 0
    e_beta_root = 0.0

    # ---- R4-b (combinatorial, X-independent): can two faces share two M'-vertices?
    for cs in corpus()
        owners = Dict{Int,Vector{Int}}()                 # prime vertex -> faces
        for f in 1:K.n_faces(cs.m), pv in cs.c.prime_faces[f]
            push!(get!(owners, pv, Int[]), f)
        end
        shared = Dict{Tuple{Int,Int},Int}()              # face pair -> # shared
        for (pv, fs) in owners
            u = sort(unique(fs))
            for i in eachindex(u), j in i+1:length(u)
                shared[(u[i], u[j])] = get(shared, (u[i], u[j]), 0) + 1
            end
        end
        for (pr, cnt) in shared
            n_face_pairs += 1
            cnt >= 2 && (n_pairs_two_shared += 1)
        end
    end

    for cs in corpus()
        isempty(cs.c.hinge_edges) && continue
        he = copy(cs.c.hinge_edges)
        for rep in 0:129
            X = sample(cs, rng, rep == 0 ? 0.0 : 0.10)
            all_faces_positive(cs.m, X) || continue      # faces stay CCW: (T5.1)
            n_samples += 1
            P = potential(cs, X)
            B = my_basis(cs, X, P)
            K.shuffle!(he, rng)
            NE = min(3, length(he))
            for ei in 1:NE
                e = he[ei]
                ed = cs.m.edges[e]
                f = cs.m.half_edges[ed.he[1]].face
                g = cs.m.half_edges[ed.he[2]].face
                v = cs.c.hinge_dir[e].src
                dst = cs.c.hinge_dir[e].dst
                vf = cs.m.faces[f]; vg = cs.m.faces[g]
                af = interior_angle(vf, X, v); ag = interior_angle(vg, X, v)
                (af < 0 || ag < 0) && continue
                Ff = K.Vec2[X[vv] for vv in vf]; Fg = K.Vec2[X[vv] for vv in vg]
                # Lemma T4.2 is stated for SIMPLE polygons; non-simple ones are counted and excluded.
                simple = poly_is_simple(Ff) && poly_is_simple(Fg)
                if !simple
                    n_nonsimple += 1; continue
                end
                beta = 2pi - af - ag
                r = 0.45 * min(sector_radius(vf, X, v), sector_radius(vg, X, v))
                (r > 1e-6) || continue

                # ---- R4-a: the equivalence (T5.2b''-1), tested geometrically.
                function local_overlap(th)
                    c = cos(th / 2); s = sin(th / 2)
                    Pf = K.Vec2[c * row(B.C, pv) + s * row(B.S, pv) for pv in cs.c.prime_faces[f]]
                    Pg = K.Vec2[c * row(B.C, pv) + s * row(B.S, pv) for pv in cs.c.prime_faces[g]]
                    idx = findlast(==(v), vf)
                    p = Pf[idx]                       # the shared hinge point
                    margin = 1e-7 * r
                    for k in 1:NRAD
                        rho = r * k / (NRAD + 1)
                        for a in 0:NPHI-1
                            phi = 2pi * a / NPHI
                            z = p + rho * K.Vec2(cos(phi), sin(phi))
                            (strictly_inside(Pf, z, margin) && strictly_inside(Pg, z, margin)) && return true
                        end
                    end
                    return false
                end

                if beta + delta < pi                     # above the transition
                    th = min(Float64(pi), beta + delta + 0.4 * (pi - beta - delta))
                    n_above += 1; n_local_tests += 1
                    if !local_overlap(th)
                        n_fail_above += 1
                        @printf("    [above-FAIL] %s af=%.6f ag=%.6f beta=%.6f th=%.6f r=%.3e\n", cs.name, af, ag, beta, th, r)
                    end
                elseif beta <= 0
                    th = 0.5 * pi
                    n_above += 1; n_local_tests += 1
                    if !local_overlap(th)
                        n_fail_above += 1
                        @printf("    [above-FAIL neg beta] %s af=%.6f ag=%.6f beta=%.6f th=%.6f r=%.3e\n", cs.name, af, ag, beta, th, r)
                    end
                end
                if beta - delta > 0                      # below the transition
                    th = min(Float64(pi), beta - delta)
                    n_below += 1; n_local_tests += 1
                    if local_overlap(th)
                        n_fail_below += 1
                        @printf("    [below-FAIL] %s af=%.6f ag=%.6f beta=%.6f th=%.6f r=%.3e\n", cs.name, af, ag, beta, th, r)
                    end
                end

                # ---- R4-c: the adjacent-edge collinearity pair of T4.4 at theta = beta_e.
                if beta > 1e-3 && beta < 2pi - 1e-3
                    if_v, if_far = far_of(vf, v, dst)
                    ig_v, ig_far = far_of(vg, v, dst)
                    if if_far > 0 && ig_far > 0
                        lf = norm(X[vf[if_far]] - X[v])
                        lg = norm(X[vg[ig_far]] - X[v])
                        if lf <= lg                       # shorter far edge belongs to f
                            w = cs.c.prime_faces[f][if_far]
                            a = cs.c.prime_faces[g][ig_v]; b = cs.c.prime_faces[g][ig_far]
                        else
                            w = cs.c.prime_faces[g][ig_far]
                            a = cs.c.prime_faces[f][if_v]; b = cs.c.prime_faces[f][if_far]
                        end
                        h = my_orient(row(B.C, b) - row(B.C, a), row(B.S, b) - row(B.S, a),
                                      row(B.C, w) - row(B.C, a), row(B.S, w) - row(B.S, a))
                        sc = abs(h.p) + abs(h.q) + abs(h.r)
                        n_beta_pairs += 1
                        gscale = maximum(norm, X)
                        if sc <= 1e-11 * gscale * gscale
                            n_beta_degenerate += 1        # would be struck by T3.H.1: FATAL
                        else
                            hv = h.p + h.q * cos(beta) + h.r * sin(beta)
                            e_beta_root = max(e_beta_root, abs(hv) / sc)
                        end
                    end
                end
            end
        end
    end

    @printf("  R4-a %d samples, %d local tests at the hinge (%d above beta_e, %d below); failures above %d, below %d\n",
            n_samples, n_local_tests, n_above, n_below, n_fail_above, n_fail_below)
    println("  R4-a hinges EXCLUDED because a face polygon is not simple: $n_nonsimple")
    println("  R4-b $n_face_pairs face pairs sharing an M'-vertex; sharing TWO or more: $n_pairs_two_shared")
    @printf("  R4-c %d beta_e collinearity pairs; identically zero: %d; max |h(beta_e)|/scale = %.3e\n",
            n_beta_pairs, n_beta_degenerate, e_beta_root)
    record("R4-a1", "(T5.2b''-1) theta > beta_e => local overlap in B(p,r)", n_above, n_fail_above, 0.0)
    record("R4-a2", "(T5.2b''-1) theta < beta_e => NO local overlap in B(p,r)", n_below, n_fail_below, 0.0)
    record("R4-a3", "hinges excluded: face polygon not simple (T4.2 hypothesis)", n_nonsimple, 0.0, 0.0)
    record("R4-b", "Case B is vacuous: no face pair shares two M'-vertices", n_face_pairs, n_pairs_two_shared, 0.0)
    record("R4-c1", "beta_e collinearity pair is NOT identically zero", n_beta_pairs, n_beta_degenerate, 0.0)
    record("R4-c2", "its orientation harmonic vanishes at theta = beta_e", n_beta_pairs, e_beta_root, 1e-9)
    # C++: 1306 samples, 5548 local tests (1731 above, 3817 below), 99 non-simple,
    # 933 face pairs, 3819 beta pairs
    @test (n_samples, n_local_tests, n_above, n_below) == (1306, 5548, 1731, 3817)
    @test n_nonsimple == 99
    @test n_face_pairs == 933
    @test n_beta_pairs == 3819
    @test n_samples >= 1000
    @test n_above > 0
    @test n_below > 0
    @test n_fail_above == 0
    @test n_fail_below == 0
    @test n_pairs_two_shared == 0
    @test n_beta_degenerate == 0
    @test e_beta_root <= 1e-9
end

# ============================================================================
# R5 -- round-5 items (derivations/core.md [D2]).
# ============================================================================
ang_of(d::K.Vec2) = atan(d[2], d[1])
function wrap2pi(a::Float64)
    t = rem(a, 2 * Float64(pi))   # C fmod: sign of the dividend
    return t < 0 ? t + 2pi : t
end
# open arc (s1, s1+l1) meets open arc (s2, s2+l2) on the circle?  0 < l < 2pi.
function open_arcs_meet(s1, l1, s2, l2)
    x = wrap2pi(s2 - s1)
    return (x < l1) || (x > 2pi - l2)
end
# how far the configuration is from flipping that boolean (for a guard band)
function arcs_margin(s1, l1, s2, l2)
    x = wrap2pi(s2 - s1)
    return min(abs(x - l1), abs(x - (2pi - l2)))
end

# (4a) of R5 as a function (a top-level loop with closures runs ~50x slower)
function r5_synthetic_roots(r2::K.MT19937)
    n_beta = 0; n_roots_gt1 = 0; n_roots_0 = 0; n_root_is_beta = 0; n_beta_in_range = 0
    n_enum_fail = 0; n_pred_change_at_shift = 0
    e_root = 0.0
    NS = 4000
    for i in 1:20000
        beta = K.uniform_real(r2, -2pi, 2pi); L = K.uniform_real(r2, 0.2, 5.0); Lp = K.uniform_real(r2, 0.2, 5.0)
        h(th) = L * Lp * sin(beta - th)
        roots = Float64[]
        for j in 0:NS-1
            t0 = pi * j / NS; t1 = pi * (j + 1) / NS
            h0 = h(t0); h1 = h(t1)
            (j > 0 && h0 == 0.0) && push!(roots, t0)
            if h0 * h1 < 0
                lo = t0; hi = t1
                for it in 1:200
                    mid = 0.5 * (lo + hi)
                    if h(lo) * h(mid) <= 0
                        hi = mid
                    else
                        lo = mid
                    end
                end
                push!(roots, 0.5 * (lo + hi))
            end
        end
        n_beta += 1
        length(roots) > 1 && (n_roots_gt1 += 1)
        isempty(roots) && (n_roots_0 += 1)
        if beta > 0 && beta < pi
            n_beta_in_range += 1
            if length(roots) == 1
                e_root = max(e_root, abs(roots[1] - beta))
                abs(roots[1] - beta) < 1e-6 && (n_root_is_beta += 1)
            end
        end
        # core.md prints: "exactly one of {beta, beta+pi, beta-pi} lies in (0,pi)"
        cnt = count(k -> (z = beta + k * pi; z > 0 && z < pi), -1:1)
        cnt != 1 && (n_enum_fail += 1)
        # the sector-overlap predicate "theta > beta_e" does not change at beta +- pi
        for k in (-1, 1)
            z = beta + k * pi
            (z <= 1e-6 || z >= pi - 1e-6) && continue
            ((z - 1e-4 > beta) != (z + 1e-4 > beta)) && (n_pred_change_at_shift += 1)
        end
    end
    return (n_beta, n_roots_gt1, n_roots_0, n_root_is_beta, n_beta_in_range, n_enum_fail,
            n_pred_change_at_shift, e_root)
end

# C++: 0.469 s
@testset "R5 same-sigma pairs: pure relative translation, constant local predicate; and the zero set of h_o,pi' in (0,pi)" begin
    rng = K.MT19937(50051)
    NTH = 13                       # theta grid inside (0, pi]
    TH(i) = pi * (i + 1) / NTH

    # ---------------- (2) relative rotation and the local predicate ----------------
    n_same = 0; n_opp = 0
    e_same_angle = 0.0; e_opp_angle = 0.0
    n_pred_same = 0; n_pred_same_flip = 0
    n_pred_opp = 0; n_pred_opp_flip = 0
    for cs in corpus(), rep in 0:24
        X = sample(cs, rng, rep == 0 ? 0.0 : 0.10)
        all_faces_positive(cs.m, X) || continue
        P = potential(cs, X)
        B = my_basis(cs, X, P)
        F = K.n_faces(cs.m)
        F < 2 && continue
        for trial in 1:12
            f = K.uniform_int(rng, 0, F - 1) + 1
            g = K.uniform_int(rng, 0, F - 1) + 1
            f == g && continue
            P.comp[f] != P.comp[g] && continue           # same rigid component only
            vf = cs.m.faces[f]; vg = cs.m.faces[g]
            (length(vf) < 3 || length(vg) < 3) && continue
            nf = length(vf); ng = length(vg)
            ia = K.uniform_int(rng, 0, nf - 1) + 1
            iw = K.uniform_int(rng, 0, ng - 1) + 1
            pa = cs.c.prime_faces[f][ia]
            pb = cs.c.prime_faces[f][mod1(ia + 1, nf)]
            pw = cs.c.prime_faces[g][iw]
            pwn = cs.c.prime_faces[g][mod1(iw + 1, ng)]        # w -> w_next
            pwp = cs.c.prime_faces[g][mod1(iw - 1, ng)]        # w -> w_prev
            same = (cs.m.sigma[f] == cs.m.sigma[g])
            d0 = 0.0; pred0 = false; have0 = false
            for i in 0:NTH-1
                th = TH(i); c = cos(th / 2); s = sin(th / 2)
                Y(pv) = c * row(B.C, pv) + s * row(B.S, pv)
                ab = Y(pb) - Y(pa)
                e1 = Y(pwn) - Y(pw); e2 = Y(pwp) - Y(pw)
                (norm(ab) < 1e-9 || norm(e1) < 1e-9 || norm(e2) < 1e-9) && break
                # J: open half-circle of directions with det(ab, d) > 0
                sJ = ang_of(ab); lJ = Float64(pi)
                # A: open sector of g at w, from dir(w->w_next) CCW to dir(w->w_prev)
                sA = ang_of(e1)
                lA = wrap2pi(ang_of(e2) - sA)
                (lA < 1e-9 || lA > 2pi - 1e-9) && break
                drift = wrap2pi(sA - sJ)
                pr = open_arcs_meet(sA, lA, sJ, lJ)
                arcs_margin(sA, lA, sJ, lJ) < 1e-7 && break
                if !have0
                    d0 = drift; pred0 = pr; have0 = true
                    continue
                end
                # relative direction drift:  delta(theta) = delta0 + (sigma_f - sigma_g) theta/2
                dd = wrap2pi(drift - d0)
                dd > pi && (dd -= 2pi)              # signed drift in (-pi, pi]
                if same
                    n_same += 1
                    e_same_angle = max(e_same_angle, abs(dd))
                    n_pred_same += 1
                    pr != pred0 && (n_pred_same_flip += 1)
                else
                    n_opp += 1
                    expect = -Float64(cs.m.sigma[g]) * (th - TH(0))
                    diff = wrap2pi(dd - expect)
                    diff > pi && (diff -= 2pi)
                    e_opp_angle = max(e_opp_angle, abs(diff))
                    n_pred_opp += 1
                    pr != pred0 && (n_pred_opp_flip += 1)
                end
            end
        end
    end

    # ---------------- (4) zero set of h_o,pi' = +- L L' sin(beta_e - theta) ---------
    # (4a) synthetic: roots of sin(beta - theta) inside the OPEN interval (0, pi).
    (n_beta, n_roots_gt1, n_roots_0, n_root_is_beta, n_beta_in_range, n_enum_fail,
     n_pred_change_at_shift, e_root) = r5_synthetic_roots(K.MT19937(50052))

    # (4b) corpus: the closed form (T5.2b''-1b) h(theta) = +- L L' sin(beta_e - theta)
    n_form = 0; e_form = 0.0
    for cs in corpus()
        isempty(cs.c.hinge_edges) && continue
        he = copy(cs.c.hinge_edges)
        for rep in 0:29
            X = sample(cs, rng, rep == 0 ? 0.0 : 0.10)
            all_faces_positive(cs.m, X) || continue
            P = potential(cs, X)
            B = my_basis(cs, X, P)
            K.shuffle!(he, rng)
            for ei in 1:min(3, length(he))
                e = he[ei]
                ed = cs.m.edges[e]
                f = cs.m.half_edges[ed.he[1]].face
                g = cs.m.half_edges[ed.he[2]].face
                v = cs.c.hinge_dir[e].src; dst = cs.c.hinge_dir[e].dst
                vf = cs.m.faces[f]; vg = cs.m.faces[g]
                af = interior_angle(vf, X, v); ag = interior_angle(vg, X, v)
                (af < 0 || ag < 0) && continue
                beta = 2pi - af - ag
                if_v, if_far = far_of(vf, v, dst)
                ig_v, ig_far = far_of(vg, v, dst)
                (if_v < 1 || ig_v < 1 || if_far < 1 || ig_far < 1) && continue
                lf = norm(X[vf[if_far]] - X[v]); lg = norm(X[vg[ig_far]] - X[v])
                if lf <= lg
                    w = cs.c.prime_faces[f][if_far]
                    a = cs.c.prime_faces[g][ig_v]; b = cs.c.prime_faces[g][ig_far]
                    L = lg; Lp = lf
                else
                    w = cs.c.prime_faces[g][ig_far]
                    a = cs.c.prime_faces[f][if_v]; b = cs.c.prime_faces[f][if_far]
                    L = lf; Lp = lg
                end
                h = my_orient(row(B.C, b) - row(B.C, a), row(B.S, b) - row(B.S, a),
                              row(B.C, w) - row(B.C, a), row(B.S, w) - row(B.S, a))
                sc = abs(h.p) + abs(h.q) + abs(h.r)
                gscale = maximum(norm, X)
                sc <= 1e-11 * gscale * gscale && continue
                # fix the overall sign at one angle, then check the whole grid
                sgn = 0; worst = 0.0
                for i in 0:39
                    th = pi * (i + 0.5) / 40.0
                    hv = h.p + h.q * cos(th) + h.r * sin(th)
                    model = L * Lp * sin(beta - th)
                    if sgn == 0
                        abs(model) > 1e-3 * L * Lp && (sgn = (hv * model >= 0) ? 1 : -1)
                    end
                    sgn != 0 && (worst = max(worst, abs(hv - sgn * model) / (L * Lp)))
                end
                if sgn != 0
                    n_form += 1; e_form = max(e_form, worst)
                end
            end
        end
    end

    @printf("  R5-a same-sigma %d comparisons, max relative-direction drift %.3e; predicate flips %d/%d\n",
            n_same, e_same_angle, n_pred_same_flip, n_pred_same)
    @printf("  R5-a opposite-sigma %d comparisons, max |drift - theta| %.3e; predicate flips %d/%d (power check, must be > 0)\n",
            n_opp, e_opp_angle, n_pred_opp_flip, n_pred_opp)
    @printf("  R5-c %d synthetic beta: >1 root in (0,pi) %d, 0 roots %d; beta in (0,pi): %d of which root==beta %d (max err %.3e)\n",
            n_beta, n_roots_gt1, n_roots_0, n_beta_in_range, n_root_is_beta, e_root)
    @printf("  R5-c core.md's {beta,beta+pi,beta-pi} enumeration fails on %d/%d; predicate 'theta>beta' changes at beta+-pi: %d\n",
            n_enum_fail, n_beta, n_pred_change_at_shift)
    @printf("  R5-d (T5.2b''-1b) closed form on %d corpus pairs, max err %.3e\n", n_form, e_form)
    record("R5-a1", "sigma_f == sigma_g => relative direction is theta-independent", n_same, e_same_angle, 1e-9)
    record("R5-a2", "sigma_f == sigma_g => local predicate (T5.2b''-1c) never flips", n_pred_same, n_pred_same_flip, 0.0)
    record("R5-a3", "power check: opposite sigma drifts by exactly theta", n_opp, e_opp_angle, 1e-9)
    record("R5-c1", "h_o,pi' has AT MOST ONE root in (0,pi)", n_beta, n_roots_gt1, 0.0)
    record("R5-c2", "beta in (0,pi) => that root IS beta", n_beta_in_range, n_beta_in_range - n_root_is_beta, 0.0)
    record("R5-c3", "'theta > beta_e' does not change value at beta_e +- pi", n_beta, n_pred_change_at_shift, 0.0)
    record("R5-d", "(T5.2b''-1b) h = +- L L' sin(beta_e - theta) on the corpus", n_form, e_form, 1e-9)
    # C++: 14496 same / 17208 opposite (5303 flips), 5038 beta in range, 5008 enum
    # failures, 894 corpus pairs
    @test (n_same, n_opp, n_pred_opp_flip) == (14496, 17208, 5303)
    @test (n_beta_in_range, n_enum_fail, n_form) == (5038, 5008, 894)
    @test n_same > 1000
    @test n_opp > 1000
    @test n_form > 100
    @test e_same_angle <= 1e-9
    @test n_pred_same_flip == 0
    @test e_opp_angle <= 1e-9
    @test n_pred_opp_flip > 0           # the predicate test really can flip
    @test n_roots_gt1 == 0
    @test n_beta_in_range - n_root_is_beta == 0
    @test n_pred_change_at_shift == 0
    @test e_form <= 1e-9
end

# ============================================================================
# R6 -- the amended NOROOT (admissible roots only).
# ============================================================================
# C++: 0.073 s
@testset "R6 amended NOROOT: Case A's substitute root is admissible; completeness is false" begin
    rng = K.MT19937(60061)
    n_pairs = 0; n_bad = 0; n_degen = 0; n_strict = 0; n_equal = 0
    lo_ratio = 1e300; hi_ratio = -1e300; e_ratio = 0.0; e_root = 0.0
    for cs in corpus()
        isempty(cs.c.hinge_edges) && continue
        for rep in 0:39
            X = sample(cs, rng, rep == 0 ? 0.0 : 0.12)
            B = K.deploy_basis(cs.c, X)
            for e in cs.c.hinge_edges
                ed = cs.m.edges[e]
                f = cs.m.half_edges[ed.he[1]].face
                g = cs.m.half_edges[ed.he[2]].face
                v = cs.c.hinge_dir[e].src; dst = cs.c.hinge_dir[e].dst
                vf = cs.m.faces[f]; vg = cs.m.faces[g]
                af = interior_angle(vf, X, v); ag = interior_angle(vg, X, v)
                (af < 0 || ag < 0) && continue
                beta = 2pi - af - ag
                (beta > 1e-3 && beta < pi - 1e-3) || continue   # theta* = beta_e in (0, eps)
                Ff = K.Vec2[X[vv] for vv in vf]; Fg = K.Vec2[X[vv] for vv in vg]
                (poly_is_simple(Ff) && poly_is_simple(Fg)) || continue   # T4.2 hypothesis
                if_v, if_far = far_of(vf, v, dst)
                ig_v, ig_far = far_of(vg, v, dst)
                (if_v < 1 || ig_v < 1 || if_far < 1 || ig_far < 1) && continue
                lf = norm(X[vf[if_far]] - X[v]); lg = norm(X[vg[ig_far]] - X[v])
                if lf <= lg
                    w = cs.c.prime_faces[f][if_far]
                    a = cs.c.prime_faces[g][ig_v]; b = cs.c.prime_faces[g][ig_far]
                    Lp = lf; L = lg
                else
                    w = cs.c.prime_faces[g][ig_far]
                    a = cs.c.prime_faces[f][if_v]; b = cs.c.prime_faces[f][if_far]
                    Lp = lg; L = lf
                end
                U = K.basis_c(B, b) - K.basis_c(B, a); Vv = K.basis_s(B, b) - K.basis_s(B, a)
                Pp = K.basis_c(B, w) - K.basis_c(B, a); Q = K.basis_s(B, w) - K.basis_s(B, a)
                det_ = K.orient_from_vectors(U, Vv, Pp, Q)
                D = K.dot_from_vectors(U, Vv, Pp, Q)
                L2 = K.dot_from_vectors(U, Vv, U, Vv)
                sc = abs(det_.p) + abs(det_.q) + abs(det_.r)
                if sc <= 1e-11 * max(norm(U) * norm(Pp), 1e-300)
                    n_degen += 1; continue
                end
                n_pairs += 1
                e_root = max(e_root, abs(heval(det_, beta)) / sc)
                s = heval(D, beta); l2 = heval(L2, beta)
                ratio = s / l2
                lo_ratio = min(lo_ratio, ratio)
                hi_ratio = max(hi_ratio, ratio)
                e_ratio = max(e_ratio, abs(ratio - Lp / L))
                ratio < 1.0 - 1e-9 ? (n_strict += 1) : (n_equal += 1)
                tol = 1e-12 * (abs(L2.p) + K.amp(L2))
                if s < -tol || s > l2 + tol
                    n_bad += 1
                    @printf("    [R6 ADM-FAIL] %-18s beta=%.6f s/|e|^2=%.6f L'/L=%.6f\n", cs.name, beta, ratio, Lp / L)
                end
            end
        end
    end
    @printf("  R6-a %d Case-A substitute pairs (identically zero: %d); s/|e|^2 in [%.6f, %.6f]; strictly inside %d, at the far endpoint %d\n",
            n_pairs, n_degen, lo_ratio, hi_ratio, n_strict, n_equal)

    # ---- R6-b: completeness counterexample (the hexagon graze).
    n_counterex = 0; n_grazes = 0
    for cs in corpus()
        B = K.deploy_basis(cs.c, cs.X0)
        sd = K.swept_discs(cs.c, B)
        pairs = K.candidate_pairs(cs.c, sd, Float64(pi), false)
        cand = K.contact_angles(cs.c, B, pairs, 1e-9, Float64(pi), 1e-9)
        orep = K.exact_theta_max_overlap(cs.c, B, pairs, 1e-9, Float64(pi), 1e-9)
        (isempty(cand) || !(cand[1] < orep.theta_max - 1e-6)) && continue
        n_grazes += 1
        eps = 0.5 * (cand[1] + orep.theta_max)
        cert = K.validity_certificate(cs.c, B, cs.X0, pairs, eps, 1e-9)
        @printf("  R6-b %-18s theta_1=%.6f Theta_max=%.6f eps=%.6f  Theta_max>=eps=%d  cert(pos,noov,noroot)=(%d,%d,%d) first admissible root=%.6f\n",
                cs.name, cand[1], orep.theta_max, eps, Int(orep.theta_max >= eps), Int(cert.pos), Int(cert.nooverlap), Int(cert.noroot), cert.first_root)
        (orep.theta_max >= eps && !K.valid(cert)) && (n_counterex += 1)
        @test cs.name == "hexagons"
        @test isapprox(cert.first_root, 1.047198; atol = 1e-6)
        @test (cert.pos, cert.nooverlap, cert.noroot) == (true, true, false)
    end

    record("R6-a1", "Case A substitute root at beta_e is ADMISSIBLE (on the segment)", n_pairs, n_bad, 0.0)
    record("R6-a2", "its projection ratio <w-a,b-a>/|b-a|^2 equals L'/L", n_pairs, e_ratio, 1e-9)
    record("R6-a3", "h_o,pi' vanishes at beta_e (round-5 claim, re-run here)", n_pairs, e_root, 1e-9)
    record("R6-b", "COMPLETENESS is false: a graze with Theta_max >= eps is not certified", n_grazes, n_counterex == 0, 0.0)
    # C++: 9642 pairs, 0 degenerate, ratio in [0.007307, 1], 9062 strict / 580 equal, 1 graze
    @test (n_pairs, n_degen, n_strict, n_equal, n_grazes) == (9642, 0, 9062, 580, 1)
    @test isapprox(lo_ratio, 0.007307; atol = 1e-6)
    @test n_pairs >= 100
    @test n_bad == 0
    @test e_ratio <= 1e-9
    @test e_root <= 1e-9
    @test lo_ratio > 0.0
    @test hi_ratio <= 1.0 + 1e-12
    @test n_equal > 0          # the closed interval test is load-bearing
    @test n_counterex >= 1     # completeness fails
end

# R6-c -- the B4 `voronoi_93` scan-vs-referee disagreement: the F34 fix makes
# polygons_overlap agree with the dense probe at every tolerance.
# C++: 0.010 s
@testset "R6-c B4 voronoi_93: polygons_overlap misfires at a shared hinge vertex" begin
    P1 = K.Vec2[(35.522612034574969, 2.8976170473256473), (35.919893509036022, 3.1484087335745494),
                (36.013918812893856, 3.6292236769982966), (35.070313922385914, 3.0868107146716834),
                (34.512614462367168, 2.8219234379478291)]
    P2 = K.Vec2[(34.233106280612688, 1.0182292700260629), (34.512614462367168, 2.8219234379478291),
                (35.067637209025698, 3.0923743720689583), (34.742880186886559, 3.4379946577381681),
                (33.800478898505453, 3.2058972270991357), (32.960043232689365, 2.3447653598242049),
                (33.184898605051494, 2.3139255252053297), (34.239945838106621, 1.3526146778850854)]
    poly_area(P) = 0.5 * sum(det2(P[i], P[mod1(i + 1, length(P))]) for i in eachindex(P))
    function pip(P, z)
        inside = false
        n = length(P)
        j = n
        for i in 1:n
            if ((P[i][2] > z[2]) != (P[j][2] > z[2])) &&
               (z[1] < (P[j][1] - P[i][1]) * (z[2] - P[i][2]) / (P[j][2] - P[i][2]) + P[i][1])
                inside = !inside
            end
            j = i
        end
        return inside
    end
    # both simple, both positively oriented, exactly one shared vertex (the hinge point)
    @test poly_is_simple(P1)
    @test poly_is_simple(P2)
    @test poly_area(P1) > 0
    @test poly_area(P2) > 0
    shared = count(norm(a - b) == 0.0 for a in P1, b in P2)
    @test shared == 1
    # GROUND TRUTH: a dense grid finds points inside each face and NONE inside both.
    lo = K.Vec2(1e300, 1e300); hi = K.Vec2(-1e300, -1e300)
    for P in (P1, P2), z in P
        lo = min.(lo, z); hi = max.(hi, z)
    end
    nA = 0; nB = 0; nAB = 0
    G = 1200
    for i in 0:G, j in 0:G
        z = K.Vec2(lo[1] + (hi[1] - lo[1]) * i / G, lo[2] + (hi[2] - lo[2]) * j / G)
        ia = pip(P1, z); ib = pip(P2, z)
        ia && (nA += 1)
        ib && (nB += 1)
        (ia && ib) && (nAB += 1)
    end
    # beta_e at the shared vertex, from the two interior angles (T4.4).
    function ang(P, v)
        n = length(P)
        i = findlast(p -> norm(p - v) == 0.0, P)
        a = P[mod1(i - 1, n)] - v; b = P[mod1(i + 1, n)] - v
        t = atan(det2(b, a), dot(b, a))
        t < 0 && (t += 2pi)
        return t
    end
    hinge = P1[5]
    aF = ang(P1, hinge); aG = ang(P2, hinge)
    beta = 2pi - aF - aG
    misfires = 0
    shrinks = (0.0, 1e-12, 1e-9, 1e-6, 1e-4, 1e-3)
    for s in shrinks
        (K.polygons_overlap(P1, P2, s) || K.polygons_overlap(P2, P1, s)) && (misfires += 1)
    end
    @printf("  R6-c grid inside A=%d B=%d BOTH=%d; alpha_f=%.6f alpha_g=%.6f beta_e=%.6f; polygons_overlap says OVERLAP at %d of %d tolerance values\n",
            nA, nB, nAB, aF, aG, beta, misfires, length(shrinks))
    record("R6-c1", "voronoi_93 faces (94,184): interiors are DISJOINT (dense probe)", (G + 1) * (G + 1), nAB, 0.0)
    record("R6-c2", "polygons_overlap agrees with the dense probe at every tolerance (F34 fixed)", length(shrinks), 0.0, 0.0)
    @test (nA, nB) == (73735, 342350)   # C++
    @test nA > 1000
    @test nB > 1000
    @test nAB == 0                       # ground truth: no interior overlap
    @test aG > pi                        # the hinge corner of face 184 is REFLEX
    @test isapprox(beta, 1.809342; rtol = 1e-5)
    @test misfires == 0                  # F34 fixed: no misfire at any tolerance
end

# =============================================================================
# L1 -- derivations/lemmas.md Part L1 (Checker-L, mission 2 / WP1)
# =============================================================================

# The wider L1 corpus: the 7 non-periodic Phase-2 reference tilings plus the
# kill_common population make_graph(id, 18, 46, 220), id = 0 ... , until >= 60 graphs
# are usable.  Built once, from the frozen C++ inputs (see the generator testset: all but
# hexagons_auto's sigma are reproduced by kill_common.jl), with the C++ X0/Phi.
const L1_CASES = Case[]
function l1_corpus()
    isempty(L1_CASES) || return L1_CASES
    @testset "L1 corpus X0/Phi agree with the C++" begin
        for j in DERIV_INPUTS["l1_corpus"]
            length(L1_CASES) >= 60 && break
            cs = build_case(fixture_mesh_raw(j["mesh"]), j["name"], j)
            cs.ok && push!(L1_CASES, cs)
        end
    end
    return L1_CASES
end

# The three class coefficients of core.md T5.2b.2, and the class, from an ABSOLUTE
# tolerance so that the "predicted" and the "numeric" side use the same rule.
@enum LKlass LK_ZERO = 0 LK_C1 = 1 LK_C2 = 2 LK_C3 = 3
function lclass(A, B, C, tol)
    a0 = abs(A) <= tol; b0 = abs(B) <= tol; c0 = abs(C) <= tol
    (a0 && b0 && c0) && return LK_ZERO
    c0 || return LK_C1
    b0 || return LK_C2
    return LK_C3
end

# C++: 0.336 s (24.5M candidate harmonics)
@testset "L1 (lemmas.md) class prediction from the combinatorics vs the numeric class" begin
    C = l1_corpus()
    println("L1 corpus: $(length(C)) graphs")
    @test length(C) >= 50
    @test length(C) == 60

    rng = K.MT19937(51001)
    NSAMP = 18   # (graph, sigma, t) triples = |corpus| * NSAMP  >= 1000

    n_triples = 0; n_cand = 0
    shared = 0; shared_C_nonzero = 0                     # L1.1(ii): must be 0
    notshared_C_zero = 0                                 # the H-L1 population
    notshared_accidental = 0; notshared_structural = 0
    pers_structural = 0; pers_accidental = 0; pers_rigid = 0
    pred_agree = 0; pred_disagree = 0                    # L1.2 trichotomy on coincident pairs
    zero_cls = 0; zero_and_index_equal = 0; index_equal_not_zero = 0
    zero_not_index_equal = 0
    cc_hinge = 0; cc_split = 0; cc_none = 0; cc_edge_not_at_vw = 0
    c3_hinge = 0; c3_hinge_perp = 0
    nclass = zeros(Int, 4)
    e_p = 0.0; e_q = 0.0; e_r = 0.0; e_hinge_dS = 0.0; e_split_dS = 0.0
    n_hinge_dS = 0; n_split_dS = 0

    for cs in C
        m = cs.m
        F = K.n_faces(m)
        nzero_C = Int[]        # per candidate pair: # samples with |C| <= tol
        pair_shared = Bool[]   # per candidate pair: v_w in {v_a, v_b}

        # face pair -> the shared edges of M (there can be more than one)
        shared_edges = Dict{Tuple{Int,Int},Vector{Int}}()
        for e in 1:K.n_edges(m)
            ed = m.edges[e]
            ed.n_faces != 2 && continue
            f1 = m.half_edges[ed.he[1]].face; f2 = m.half_edges[ed.he[2]].face
            push!(get!(shared_edges, (min(f1, f2), max(f1, f2)), Int[]), e)
        end

        for rep in 0:NSAMP-1
            X = sample(cs, rng, rep == 0 ? 0.0 : 0.25)
            # the 1e-3 probe of L1.1(iii): does a small move of t change the class?
            Tp = zeros(max(cs.k, 1), 2)
            if cs.k > 0
                g = K.NormalDist(0.0, 1e-3)
                for i in 1:cs.k
                    Tp[i, 1] = K.normal(g, rng); Tp[i, 2] = K.normal(g, rng)
                end
            end
            Xp = copy(X)
            if cs.k > 0
                Dm = cs.rep.Phi * Tp
                for v in eachindex(Xp)
                    Xp[v] += K.Vec2(Dm[v, 1], Dm[v, 2])
                end
            end

            pid = 0
            B = K.deploy_basis(cs.c, X)
            Bp = K.deploy_basis(cs.c, Xp)
            gs = maximum(norm, X)
            tol = 1e-11 * gs * gs
            n_triples += 1

            for f in 1:F
                vf = m.faces[f]
                pf = cs.c.prime_faces[f]
                nf = length(pf)
                sf = Float64(m.sigma[f])
                for kk in 1:nf
                    a = pf[kk]; b = pf[mod1(kk + 1, nf)]
                    va = vf[kk]; vb = vf[mod1(kk + 1, nf)]
                    d = X[vb] - X[va]
                    for g in 1:F
                        g == f && continue
                        vg = m.faces[g]
                        pg = cs.c.prime_faces[g]
                        for j in eachindex(pg)
                            w = pg[j]; vw = vg[j]
                            n_cand += 1
                            pid += 1
                            if rep == 0
                                push!(nzero_C, 0); push!(pair_shared, false)
                            end

                            h = K.orient_harmonic(B, a, b, w)
                            A_ = h.p - h.q; B_ = 2 * h.r; C_ = h.p + h.q
                            kl = lclass(A_, B_, C_, tol)
                            nclass[Int(kl)+1] += 1

                            is_shared = (vw == va || vw == vb)
                            index_eq = (w == a || w == b)
                            abs(C_) <= tol && (nzero_C[pid] += 1)
                            pair_shared[pid] = is_shared

                            # ---- L1.3 (corpus form): h == 0  <=>  w IS a or b as an M'-vertex
                            if kl == LK_ZERO
                                zero_cls += 1
                                index_eq ? (zero_and_index_equal += 1) : (zero_not_index_equal += 1)
                            elseif index_eq
                                index_equal_not_zero += 1
                            end

                            # ---- L1.1(ii): coincident copy => C == 0 exactly, for every X
                            if is_shared
                                shared += 1
                                abs(C_) > tol && (shared_C_nonzero += 1)

                                # ---- L1.2 (L1.8): the closed forms, and the class from dS alone
                                astar = (vw == va) ? a : b
                                dS = K.basis_s(B, w) - K.basis_s(B, astar)
                                p_pred = 0.5 * sf * dot(d, dS)
                                r_pred = 0.5 * det2(d, dS)
                                e_p = max(e_p, abs(h.p - p_pred))
                                e_q = max(e_q, abs(h.q + h.p))
                                e_r = max(e_r, abs(h.r - r_pred))
                                kpred = lclass(sf * dot(d, dS), det2(d, dS), 0.0, tol)
                                kpred == kl ? (pred_agree += 1) : (pred_disagree += 1)

                                # ---- the combinatorial type of L1.1(ii)
                                key = (min(f, g), max(f, g))
                                e_at_vw = 0
                                if haskey(shared_edges, key)
                                    for e in shared_edges[key]
                                        if m.edges[e].key.a == vw || m.edges[e].key.b == vw
                                            e_at_vw = e; break
                                        end
                                    end
                                end
                                if e_at_vw == 0
                                    cc_none += 1
                                    haskey(shared_edges, key) && (cc_edge_not_at_vw += 1)  # the P1-L1-a config
                                elseif cs.c.edge_type[e_at_vw] == K.Hinge
                                    cc_hinge += 1
                                    # L1.2b type (b): v_w == dst(e)  =>  dS = 2 sigma_g J (x_src - x_v)
                                    hd = cs.c.hinge_dir[e_at_vw]
                                    if hd.dst == vw
                                        sg = Float64(m.sigma[g])
                                        pred = 2.0 * sg * Jm(X[hd.src] - X[vw])
                                        e_hinge_dS = max(e_hinge_dS, norm(dS - pred))
                                        n_hinge_dS += 1
                                        if kl == LK_C3
                                            c3_hinge += 1
                                            abs(dot(d, X[hd.src] - X[hd.dst])) <= 1e-7 * gs * gs && (c3_hinge_perp += 1)
                                        end
                                    end
                                elseif cs.c.edge_type[e_at_vw] == K.Split
                                    cc_split += 1
                                    # T1.B / L1.2b type (c): dS is the SAME vector for both endpoints
                                    vother = (m.edges[e_at_vw].key.a == vw) ? m.edges[e_at_vw].key.b : m.edges[e_at_vw].key.a
                                    w2 = K.prime_vertex(cs.c, g, vother)
                                    a2 = K.prime_vertex(cs.c, f, vother)
                                    dS2 = K.basis_s(B, w2) - K.basis_s(B, a2)
                                    e_split_dS = max(e_split_dS, norm(dS - dS2))
                                    n_split_dS += 1
                                end
                            elseif abs(C_) <= tol
                                # ---- L1.1(iii)/(iv): C == 0 without being a coincident-copy pair
                                notshared_C_zero += 1
                                hp = K.orient_harmonic(Bp, a, b, w)
                                gp = maximum(norm, Xp)
                                tp = 1e-11 * gp * gp
                                kp = lclass(hp.p - hp.q, 2 * hp.r, hp.p + hp.q, tp)
                                kp != kl ? (notshared_accidental += 1) : (notshared_structural += 1)
                            end
                        end
                    end
                end
            end
        end
        for i in eachindex(nzero_C)
            pair_shared[i] && continue
            nzero_C[i] == 0 && continue
            if cs.k == 0
                pers_rigid += 1; continue
            end
            nzero_C[i] == NSAMP ? (pers_structural += 1) : (pers_accidental += 1)
        end
    end

    @printf("  L1  triples=%d  candidates=%d\n      classes: zero %d  class1 %d  class2 %d  class3 %d\n      coincident-copy pairs: %d   of which C != 0 (must be 0): %d\n        by shared-edge type at v_w:  hinge %d  split %d  none %d  (of the 'none', f,g share an edge NOT at v_w: %d)\n      class predicted from dS == numeric class: %d agree, %d disagree\n      NOT coincident but C == 0: %d  (1e-3 probe: accidental %d / structural %d)\n        per-graph persistence (L1.1(iv) criterion): structural %d  accidental %d  rigid/undecidable %d\n      h == 0 pairs: %d   w IS a or b as an M'-vertex: %d   not: %d\n        index-equal but NOT class zero: %d\n      class-3 hinge pairs: %d   of which <d, x_src - x_dst> = 0: %d\n",
            n_triples, n_cand, nclass[1], nclass[2], nclass[3], nclass[4], shared, shared_C_nonzero,
            cc_hinge, cc_split, cc_none, cc_edge_not_at_vw, pred_agree, pred_disagree, notshared_C_zero,
            notshared_accidental, notshared_structural, pers_structural, pers_accidental, pers_rigid,
            zero_cls, zero_and_index_equal, zero_not_index_equal, index_equal_not_zero, c3_hinge, c3_hinge_perp)

    record("L1-a", "L1.1(ii) coincident copy => h(0) = 0 (0 exceptions)", shared, shared_C_nonzero, 0.0)
    record("L1-b", "L1.2 (L1.8) p = sigma_f <d,dS>/2", shared, e_p, 1e-9)
    record("L1-c", "L1.2 (L1.8) q = -p", shared, e_q, 1e-9)
    record("L1-d", "L1.2 (L1.8) r = det(d,dS)/2", shared, e_r, 1e-9)
    record("L1-e", "L1.2 class from dS alone == numeric class", shared, pred_disagree, 0.0)
    record("L1-f", "L1.2b type (b) dS = 2 sigma_g J (x_src - x_v)", n_hinge_dS, e_hinge_dS, 1e-9)
    record("L1-g", "T1.B/L1.2b type (c) dS same for both split endpoints", n_split_dS, e_split_dS, 1e-9)
    record("L1-h", "L1.2b type (b) class 3 <=> tested edge PERP hinge edge", c3_hinge, c3_hinge - c3_hinge_perp, 0.0)
    record("L1-i", "L1.3 h == 0 => w is a or b as an M'-vertex", zero_cls, zero_not_index_equal, 0.0)
    record("L1-j", "L1.3 w is a or b as an M'-vertex => h == 0", zero_and_index_equal, index_equal_not_zero, 0.0)
    record("L1-k", "H-L1 REFUTED: non-coincident pairs with C == 0 at every sample", pers_structural, 0.0, 0.0)
    record("L1-k2", "L1.1(iii) accidental set (C = 0 at some samples, not all)", pers_accidental, 0.0, 0.0)
    record("L1-k3", "L1.1(ii) type-(d) as printed: f,g share an edge NOT at v_w", shared, cc_edge_not_at_vw, 0.0)

    # C++ tallies: triples 1080, candidates 24460884, classes 223848/23119218/1097596/20222,
    # coincident 835272 (hinge 447696, split 98352, none 289224), NOT coincident C==0 506394
    # (accidental 1120 / structural 505274), persistence 7624 / 1116 / 20447, class-3 hinge 1472.
    # hexagons_auto is the ONE L1 graph whose sigma differs from the C++ (see the
    # generators testset) -- but it is loaded from the frozen C++ input here, so the tallies
    # are exactly comparable; the class-tolerance counts are tie-sensitive at the 1e-11 level.
    @test (n_triples, n_cand) == (1080, 24460884)
    @test shared == 835272
    @test (cc_hinge, cc_split, cc_none, cc_edge_not_at_vw) == (447696, 98352, 289224, 0)
    @test abs(nclass[1] - 223848) <= 50 && abs(nclass[4] - 20222) <= 50
    @test abs(notshared_C_zero - 506394) <= 200
    @test abs(pers_structural - 7624) <= 20 && abs(pers_accidental - 1116) <= 20
    @test c3_hinge == 1472
    @test n_triples >= 1000
    @test shared_C_nonzero == 0
    @test pred_disagree == 0
    @test e_p <= 1e-9
    @test e_q <= 1e-9
    @test e_r <= 1e-9
    @test e_hinge_dS <= 1e-9
    @test e_split_dS <= 1e-9
    @test c3_hinge_perp == c3_hinge
    @test zero_not_index_equal == 0
    @test index_equal_not_zero == 0
    # H-L1 is FALSE and lemmas.md says so (L1.1(iv)): the unexplained instances exist and
    # are STRUCTURAL, i.e. a 1e-3 move of t does not change their class.
    @test notshared_structural > 0
    @test pers_structural > 0
end

# C++: 0.086 s
@testset "L1.2a class-3 harmonics: constant sign on (0, pi], theta = pi included" begin
    rng = K.MT19937(51003)
    n = 0; sign_wrong = 0; zeros_ = 0; pi_zero = 0
    worst_id = 0.0; worst_ratio = 1e300
    for i in 1:10000
        p = K.uniform_real(rng, -3.0, 3.0)
        abs(p) < 1e-3 && (p = (p < 0 ? -1 : 1) * 1e-3)
        h = K.Harmonic(p, -p, 0.0)                # class 3 exactly: C = 0, B = 0, A = 2p
        @test h.p + h.q == 0.0
        @test 2 * h.r == 0.0
        # h''(0) = -q = p  (L1.2a's identification of p with h''(0))
        @test abs(-h.q - p) <= 1e-15
        for j in 1:2000
            th = pi * j / 2000.0                  # j = 2000 is exactly theta = pi
            v = heval(h, th)
            worst_id = max(worst_id, abs(v - p * (1 - cos(th))))
            if v == 0.0
                zeros_ += 1
                j == 2000 && (pi_zero += 1)
            end
            ratio = v / p
            ratio <= 0 && (sign_wrong += 1)
            worst_ratio = min(worst_ratio, ratio)
        end
        # theta = pi explicitly: h(pi) = p - q cos pi = 2p, never 0
        @test abs(heval(h, Float64(pi)) - 2 * p) <= 1e-14 * abs(p)
        # and the deflated root list must be EMPTY on (0, pi]
        @test isempty(K.harmonic_roots_deflated(h, 0.0, Float64(pi)))
        @test K.classify_harmonic(h) == K.Class3
        n += 1
    end
    record("L1-l", "L1.2a class 3 == p(1 - cos theta), identically", n, worst_id, 1e-15)
    record("L1-m", "L1.2a sign h = sign p on (0,pi], no zero (pi included)", n, sign_wrong + zeros_, 0.0)
    @test sign_wrong == 0
    @test zeros_ == 0
    @test pi_zero == 0
    @test worst_ratio > 0.0
end

# =============================================================================
# L2 -- derivations/lemmas.md Part L2 (Checker-L, mission 2 / WP1)
#
# The population, quotient pipeline and achievable() below are COPIED from
# code/apps/kill_k7.cpp / kill_b3.cpp so that the numbers compare directly with
# results/kill/k7/k7_main.csv.  The genuinely INDEPENDENT part is `dtau_covector`.
# =============================================================================
const DERIV_L2 = JSON.parsefile(joinpath(CORPUS, "reference_patterns", "derivation_inputs_l2.json"))

# Eigen's Vector2d::dot is an UNFUSED a0*b0 + a1*b1 (a redux over separate statements, so
# clang's -ffp-contract=on does not touch it); StaticArrays' `dot` is a muladd that becomes
# an fma on aarch64.  The plain form is what reproduces the C++ cells bit for bit.
dot2(a::K.Vec2, b::K.Vec2) = a[1] * b[1] + a[2] * b[2]

# ---- copied from kill_k7.cpp (population); bit-identical to the C++ on all 2122 tori
function l2_voronoi_pattern(inst::Int, nsites::Int, L::Float64, rng::K.MT19937)
    P = K.PeriodicPattern()
    P.family = "voronoi_torus"
    P.name = "voronoi_torus_$(inst)_n$(nsites)"
    site = K.Vec2[]
    guard = 0
    while length(site) < nsites && (guard += 1) <= 100000
        p = K.Vec2(K.uniform_real(rng, 0.0, L), K.uniform_real(rng, 0.0, L))
        ok = true
        for s in site
            d = p - s
            # std::round: half away from zero
            d = K.Vec2(d[1] - L * round(d[1] / L, RoundNearestTiesAway), d[2] - L * round(d[2] / L, RoundNearestTiesAway))
            norm(d) < 0.15 * L / sqrt(nsites) && (ok = false)
        end
        ok && push!(site, p)
    end
    if length(site) < nsites
        P.err = "site rejection failed"; return P
    end
    rep = K.Vec2[]
    for i in -1:1, j in -1:1, s in site
        push!(rep, s + K.Vec2(i * L, j * L))
    end
    polys = Vector{K.Vec2}[]
    for i in 1:nsites
        pi_ = site[i]
        cellp = K.Vec2[pi_ + K.Vec2(-L, -L), pi_ + K.Vec2(L, -L), pi_ + K.Vec2(L, L), pi_ + K.Vec2(-L, L)]
        for pj in rep
            norm(pj - pi_) < 1e-12 && continue
            nvec = pj - pi_
            off = dot2(nvec, 0.5 * (pi_ + pj))
            out = K.Vec2[]
            nc = length(cellp)
            for k in 1:nc
                A = cellp[k]; Bv = cellp[mod1(k + 1, nc)]
                da = dot2(nvec, A) - off; db = dot2(nvec, Bv) - off
                da <= 0 && push!(out, A)
                ((da < 0 && db > 0) || (da > 0 && db < 0)) && push!(out, A + (Bv - A) * (da / (da - db)))
            end
            cellp = out
            length(cellp) < 3 && break
        end
        if length(cellp) < 3
            P.err = "empty voronoi cell"; return P
        end
        push!(polys, cellp)
    end
    cell = try
        K.mesh_from_polygons(polys, 1e-7)
    catch e
        P.err = "weld: " * sprint(showerror, e); return P
    end
    T = K.Mat2(L, 0.0, 0.0, L)
    cell.sigma = fill(-1, K.n_faces(cell))
    q0 = K.build_quotient(cell, T)
    if !q0.ok
        P.err = "quotient(probe): " * q0.err; return P
    end
    cell.sigma = K.quotient_sigma(K.n_faces(cell), K.quotient_dual(q0), rng)
    P.cell = cell; P.T = T; P.ok = true
    return P
end

# a frozen pattern row -> PeriodicPattern (faces AS STORED; T row-major in the file)
function pattern_from_fixture(j)
    P = K.PeriodicPattern()
    P.name = j["name"]; P.family = j["family"]; P.ok = j["ok"]; P.err = j["err"]
    t = j["T"]
    P.T = K.Mat2(t[1], t[3], t[2], t[4])
    P.ok && (P.cell = fixture_mesh_raw(j["cell"]))
    return P
end
same_pattern(P::K.PeriodicPattern, j) =
    P.ok == j["ok"] && (!P.ok || (same_mesh_as_fixture(P.cell, j["cell"]) && P.cell.sigma == Int[s for s in j["cell"]["orientation"]]))

# The K7 population: 7 families x {2x2, 3x2, 3x3} tiling patterns + 12 Voronoi tori.
# The tiling patterns come from the frozen C++ file: make_tiling_pattern reproduces the
# C++ cell for triangles/hexagons/trunc_square_488/t3_4_3_12, differs by trig ulps for
# snub_square and is NOT the same cell for squares (translated) and kagome (different
# vertex/face layout) -- reported to the method port.  The Voronoi tori are rebuilt here
# and checked bit-identical.
function l2_k7_population()
    out = K.PeriodicPattern[]
    for i in 1:21
        push!(out, pattern_from_fixture(DERIV_L2["k7"][i]))
    end
    ns = [20, 28, 36, 45, 55, 70, 85, 100, 120, 140, 170, 200]
    for i in 0:11
        rng = K.MT19937(9100001 + 104729 * i)
        P = l2_voronoi_pattern(i, ns[i+1], 10.0, rng)
        @test same_pattern(P, DERIV_L2["k7"][22+i])
        push!(out, P)
    end
    return out
end

# ---- copied from kill_b3.cpp (PState / shape_point / prepare / K_at / achievable) --
mutable struct L2State
    ok::Bool
    err::String
    q::K.Quotient
    sp::K.SuperPatch
    cut::K.CutStructure
    X0::Matrix{Float64}
    Phi::Matrix{Float64}
    k::Int
    med_edge::Float64
end

function l2_shape_point(S::L2State, t::Vector{Float64})
    T = zeros(max(S.k, 1), 2)
    for i in 1:S.k
        T[i, 1] = t[2i-1]; T[i, 2] = t[2i]
    end
    X = copy(S.X0)
    S.k > 0 && (X += S.Phi * T)
    return K.matrix_to_points(X)
end

function l2_prepare(P::K.PeriodicPattern)
    q = K.build_quotient(P.cell, P.T)
    q.ok || return L2State(false, "quotient: " * q.err, q, K.SuperPatch(), K.make_cut(P.cell), zeros(0, 2), zeros(0, 0), 0, 1.0)
    sys = K.quotient_system(q, q.Xq)
    sr = K.solve_system(sys, q.Xq)
    sr.projection_ok || return L2State(false, "projection failed", q, K.SuperPatch(), K.make_cut(P.cell), zeros(0, 2), zeros(0, 0), 0, 1.0)
    sp = K.build_super(q, 1)
    cut = K.make_cut(sp.mesh)
    return L2State(true, "", q, sp, cut, sr.X0, sr.Phi, sr.dim_null, K.median_edge_length(P.cell))
end

function l2_K_at(S::L2State, Xq::Vector{K.Vec2})
    K.set_super_positions!(S.sp, S.q, Xq)
    B = K.deploy_basis(S.cut, S.sp.mesh.X)
    return K.periodic_jacobian(S.sp, S.q, S.cut, B)
end

# JacobiSVD rank with the kill_b3 rule: sv > max(1e-10 * s_max, 1e-13)
function l2_rank(M::AbstractMatrix)
    s = svdvals(M)
    isempty(s) && return 0
    t = 1e-10 * s[1]
    return count(x -> x > max(t, 1e-13), s)
end

function l2_achievable(S::L2State)
    A = K.AchievableSet()
    z = zeros(2 * S.k)
    A.K0 = l2_K_at(S, l2_shape_point(S, z)).K
    A.A = zeros(4, max(2 * S.k, 1))
    A.D = zeros(max(S.k, 1), 2)
    sc = S.med_edge
    for j in 0:2*S.k-1
        e = zeros(2 * S.k)
        e[j+1] = sc
        M = (l2_K_at(S, l2_shape_point(S, e)).K - A.K0) / sc
        push!(A.M, M)
        A.A[:, j+1] = [M[1, 1], M[2, 1], M[1, 2], M[2, 2]]
        if j % 2 == 0
            MT = M * S.q.T
            A.D[j ÷ 2 + 1, 1] = 0.5 * MT[2, 1]
            A.D[j ÷ 2 + 1, 2] = 0.5 * MT[2, 2]
        end
    end
    if 2 * S.k > 0
        A.dimK = l2_rank(A.A)
        A.rankD = l2_rank(A.D)
    end
    return A
end

# ---- INDEPENDENT: the covector d_tau of (L2.5), from the LIFTED hinge graph ----
# Returns (found, coef).
function dtau_covector(q::K.Quotient, tau::K.Vec2i, reverse_edge_order::Bool = false)
    nf = K.n_faces(q.cell)
    W = 4; R = 2 * W + 1
    sid(f, i, j) = ((f - 1) * R + (i + W)) * R + (j + W) + 1   # 1-based state id
    nstates = nf * R * R
    prev = zeros(Int, nstates); cls = zeros(Int, nstates); sgn = zeros(Float64, nstates)
    seen = falses(nstates)
    order = [i for i in eachindex(q.edges) if q.edges[i].hinge]
    reverse_edge_order && reverse!(order)
    start = sid(1, 0, 0); goal = sid(1, tau[1], tau[2])
    seen[start] = true
    bfs = Int[start]
    while !isempty(bfs)
        s = popfirst!(bfs)
        s == goal && break
        s0 = s - 1
        j0 = s0 % R - W; i0 = (s0 ÷ R) % R - W; f0 = s0 ÷ (R * R) + 1
        for ei in order
            e = q.edges[ei]
            d1 = K.half_edge_direction(q.cell, e.he[2])
            delta = e.src_off - q.voff[d1[1]]
            if e.f[1] == f0
                fto = e.f[2]; oto = K.Vec2i(i0, j0) + delta
            elseif e.f[2] == f0
                fto = e.f[1]; oto = K.Vec2i(i0, j0) - delta
            else
                continue
            end
            (abs(oto[1]) > W || abs(oto[2]) > W) && continue
            t = sid(fto, oto[1], oto[2])
            seen[t] && continue
            seen[t] = true
            prev[t] = s; cls[t] = e.src_class; sgn[t] = Float64(q.cell.sigma[fto])
            push!(bfs, t)
        end
    end
    seen[goal] || return false, zeros(q.nq)
    coef = zeros(q.nq)
    s = goal
    while s != start
        coef[cls[s]] += sgn[s]
        s = prev[s]
    end
    return true, coef
end

# C++: 0.916 s
@testset "L2 (lemmas.md) dim K = 2 rank(D), the (L2.8) generators, and an INDEPENDENT D" begin
    pats = l2_k7_population()
    n_k7 = length(pats)
    # >= 100 fresh random Voronoi tori with random sigma (rebuilt; bit-identical to the C++)
    for i in 0:109
        rng = K.MT19937(4400011 + 7717 * i)
        n = K.uniform_int(rng, 12, 45)
        P = l2_voronoi_pattern(1000 + i, n, 10.0, rng)
        @test same_pattern(P, DERIV_L2["fresh110"][i+1])
        push!(pats, P)
    end

    n_ok = 0; n_k7_ok = 0; eq_fail = 0; dimK_parity_fail = 0
    rankD_hist = zeros(Int, 3)
    rankD_eq_min = 0
    e_fact_even = 0.0; e_fact_odd = 0.0; e_vanish = 0.0; e_Dind = 0.0; e_path = 0.0; e_ones = 0.0
    n_Dind = 0; n_path = 0; path_offbasis_nonzero = 0
    covec_ok = 0; covec_mismatch = 0
    n_inconsistent = 0; n_gamma_split = 0; n_Dind_bad = 0; n_copies_split = 0
    e_Dind_consistent = 0.0
    sq33_maxD = -1.0; sq33_trK = 0.0; sq33_detK = 0.0; sq33_K21 = 0.0
    sq33_k = -1; sq33_dimK = -1

    for (pi_, P) in enumerate(pats)
        P.ok || continue
        S = l2_prepare(P)
        S.ok || continue
        A = l2_achievable(S)
        n_ok += 1
        pi_ <= n_k7 && (n_k7_ok += 1)

        # H2 (Gamma connected) on the object K is measured on: the 3x3 SUPER patch.
        J0 = l2_K_at(S, l2_shape_point(S, zeros(2 * S.k)))
        consistent = J0.consistency <= 1e-9
        consistent || (n_inconsistent += 1)
        let
            F = K.n_faces(S.sp.mesh)
            adj = [Int[] for _ in 1:F]
            for e in S.cut.hinge_edges
                ed = S.sp.mesh.edges[e]
                f1 = S.sp.mesh.half_edges[ed.he[1]].face
                f2 = S.sp.mesh.half_edges[ed.he[2]].face
                push!(adj[f1], f2); push!(adj[f2], f1)
            end
            comp = zeros(Int, F)
            nc = 0
            for s0 in 1:F
                comp[s0] != 0 && continue
                nc += 1
                comp[s0] = nc
                Q = Int[s0]
                while !isempty(Q)
                    f = popfirst!(Q)
                    for g in adj[f]
                        if comp[g] == 0
                            comp[g] = nc; push!(Q, g)
                        end
                    end
                end
            end
            # the sharp predicate: the (0,0) and the (1,0)/(0,1) copies of every cell face
            # lie in the SAME component
            copies_split = false
            for f in 1:S.sp.nfc
                if comp[K.face_index(S.sp, 0, 0, f)] != comp[K.face_index(S.sp, 1, 0, f)] ||
                   comp[K.face_index(S.sp, 0, 0, f)] != comp[K.face_index(S.sp, 0, 1, f)]
                    copies_split = true
                end
            end
            nc > 1 && (n_gamma_split += 1)
            if copies_split
                n_copies_split += 1
                @printf("    [H2 violated where it matters] %-26s c(Gamma_super)=%d  the (0,0) and (1,0) copies of a face are in DIFFERENT components; consistency=%.3e\n",
                        P.name, nc, J0.consistency)
            end
        end

        # ---- L2.1 (L2.9): rank(A) == 2 rank(D)
        A.dimK != 2 * A.rankD && (eq_fail += 1)
        (A.dimK != 0 && A.dimK != 2 && A.dimK != 4) && (dimK_parity_fail += 1)
        rankD_hist[min(A.rankD, 2)+1] += 1
        A.rankD == min(2, S.k) && (rankD_eq_min += 1)

        # ---- L2.1 (L2.8): M_{2i} = 2 e_y d_i^T P0^{-1},  M_{2i+1} = -2 e_x d_i^T P0^{-1}
        P0 = S.q.T; P0i = inv(P0)
        for i in 1:S.k
            di = A.D[i, :]
            Fe = K.Mat2(0.0, 2 * di[1], 0.0, 2 * di[2]) * P0i      # row 2 = 2 d_i^T
            Fo = K.Mat2(-2 * di[1], 0.0, -2 * di[2], 0.0) * P0i    # row 1 = -2 d_i^T
            Me = A.M[2i-1]; Mo = A.M[2i]
            sc = max(1.0, maximum(abs, Me))
            e_fact_even = max(e_fact_even, maximum(abs, Me - Fe) / sc)
            so = max(1.0, maximum(abs, Mo))
            e_fact_odd = max(e_fact_odd, maximum(abs, Mo - Fo) / so)
            # the row of M_j P0 that (L2.8) says vanishes
            MeT = Me * P0; MoT = Mo * P0
            s2 = max(1.0, maximum(abs, MeT))
            e_vanish = max(e_vanish, maximum(abs, MeT[1, :]) / s2)
            s3 = max(1.0, maximum(abs, MoT))
            e_vanish = max(e_vanish, maximum(abs, MoT[2, :]) / s3)
        end

        # ---- INDEPENDENT D from (L2.5), and the two facts L2.2(b)/(c) lean on
        okh, ch = dtau_covector(S.q, K.Vec2i(1, 0))
        okv, cv = dtau_covector(S.q, K.Vec2i(0, 1))
        okh2, ch2 = dtau_covector(S.q, K.Vec2i(1, 0), true)
        if okh && okv && S.k > 0
            Dind = hcat(S.Phi' * ch, S.Phi' * cv)
            sc = max(1.0, maximum(abs, A.D))
            ed = maximum(abs, Dind - A.D) / sc
            if ed > 1e-9
                n_Dind_bad += 1
                @printf("    [Dind mismatch] %-26s k=%3d rankD=%d err=%.3e  consistency=%.3e\n", P.name, S.k, A.rankD, ed, J0.consistency)
            end
            e_Dind = max(e_Dind, ed)
            consistent && (e_Dind_consistent = max(e_Dind_consistent, ed))
            n_Dind += 1
            # L2.2(b)/(c) need d_tau(1) = 0 (bipartite Gamma: every lifted walk has EVEN length)
            e_ones = max(e_ones, max(abs(sum(ch)), abs(sum(cv))))
            # (c) rank(D) = rank([L; d_h; d_v]) - rank(L)
            if S.q.nq <= 80
                Laug = vcat(S.q.L, ch', cv')
                (l2_rank(Laug) - l2_rank(S.q.L) == A.rankD) ? (covec_ok += 1) : (covec_mismatch += 1)
            end
        end
        if okh && okh2 && S.k > 0
            # L2.0c: path-independence holds ON ker[L; e_pin] ...
            diff = ch - ch2
            e_path = max(e_path, maximum(abs, S.Phi' * diff))
            n_path += 1
            # ... but NOT as covectors on all of R^{n_q}
            maximum(abs, diff) > 1e-12 && (path_offbasis_nonzero += 1)
        end

        if P.name == "squares_3x3"
            sq33_maxD = maximum(abs, A.D)
            sq33_trK = tr(A.K0)
            sq33_detK = det(A.K0)
            sq33_K21 = A.K0[2, 1]
            sq33_k = S.k
            sq33_dimK = A.dimK
        end
    end

    @printf("  L2  patterns evaluated: %d  (K7 population: %d of %d)\n      rank(A) == 2 rank(D): %d / %d   failures %d\n      dim K in {0,2,4}: failures %d\n      rank(D) distribution: 0 -> %d   1 -> %d   2 -> %d\n      rank(D) == min(2,k)  [L2.2(e) conjecture]: %d / %d\n      max rel |M_2i   - 2 e_y d_i^T P0^-1| (CIRCULAR, D is read off M_2i) : %.3e\n      max rel |M_2i+1 + 2 e_x d_i^T P0^-1| (INDEPENDENT of how D is built): %.3e\n      max rel vanishing row of M_j P0                                     : %.3e\n      max rel |D from (L2.5) - AchievableSet::D|  (%d patterns)         : %.3e\n        restricted to patterns whose period is CONSISTENT                  : %.3e\n      patterns with PeriodicJac::consistency > 1e-9 : %d   c(Gamma_super) > 1 : %d   face copies in different components : %d   Dind mismatches: %d\n      d_tau(1) == 0 (bipartite parity; needed by L2.2(b)/(c))             : %.3e\n      path-independence on ker (%d patterns)                            : %.3e\n        of those, the two path covectors differ OFF ker: %d\n      L2.2(c) rank([L;d_h;d_v]) - rank(L) == rank(D): %d ok, %d mismatch\n      squares_3x3: k=%d dimK=%d max|D|=%.3e  tr K0=%.6f det K0=%.6f K0(1,0)=%+.6f\n",
            n_ok, n_k7_ok, n_k7, n_ok - eq_fail, n_ok, eq_fail, dimK_parity_fail,
            rankD_hist[1], rankD_hist[2], rankD_hist[3], rankD_eq_min, n_ok,
            e_fact_even, e_fact_odd, e_vanish, n_Dind, e_Dind, e_Dind_consistent,
            n_inconsistent, n_gamma_split, n_copies_split, n_Dind_bad, e_ones, n_path, e_path,
            path_offbasis_nonzero, covec_ok, covec_mismatch,
            sq33_k, sq33_dimK, sq33_maxD, sq33_trK, sq33_detK, sq33_K21)

    record("L2-a", "L2.1 (L2.9) rank(A) == 2 rank(D)", n_ok, eq_fail, 0.0)
    record("L2-b", "L2.1 dim K in {0,2,4} -- never 1, never 3", n_ok, dimK_parity_fail, 0.0)
    record("L2-c", "L2.1 (L2.8) M_2i   = 2 e_y d_i^T P0^-1 (circular)", n_ok, e_fact_even, 1e-10)
    record("L2-d", "L2.1 (L2.8) M_2i+1 = -2 e_x d_i^T P0^-1 (non-circular)", n_ok, e_fact_odd, 1e-10)
    record("L2-e", "L2.1 (L2.8) the vanishing row of M_j P0", n_ok, e_vanish, 1e-10)
    record("L2-f", "L2.0c D rebuilt from (L2.5) == AchievableSet::D (H2 holds)", n_Dind - n_Dind_bad, e_Dind_consistent, 1e-9)
    record("L2-f2", "the D mismatches are EXACTLY the inconsistent-period patterns", n_ok, abs(n_Dind_bad - n_inconsistent), 0.0)
    record("L2-g", "L2.2(b)/(c) need d_tau(1) = 0 (bipartite parity, unprinted)", n_Dind, e_ones, 1e-12)
    record("L2-h", "L2.0c path-independence holds ON ker[L;e_pin]", n_path, e_path, 1e-9)
    record("L2-i", "L2.2(c) rank([L;d_h;d_v]) - rank(L) == rank(D)", covec_ok + covec_mismatch, covec_mismatch, 0.0)

    # C++: 143 patterns (33 K7), rank(D) 8/1/134, 142 with rank(D) == min(2,k), 134 D
    # rebuilt, 2 inconsistent, 26 c(Gamma_super) > 1, 4 copies split, 2 Dind mismatches,
    # 134 path patterns of which 120 differ off ker, 108 ok covector ranks
    @test (n_ok, n_k7_ok) == (143, 33)
    @test rankD_hist == [8, 1, 134]
    @test rankD_eq_min == 142
    @test (n_Dind, n_inconsistent, n_gamma_split, n_copies_split, n_Dind_bad) == (134, 2, 26, 4, 2)
    @test (n_path, path_offbasis_nonzero, covec_ok) == (134, 120, 108)
    @test n_k7_ok == 33
    @test n_ok >= 133
    @test eq_fail == 0
    @test dimK_parity_fail == 0
    @test e_fact_even <= 1e-10
    @test e_fact_odd <= 1e-10
    @test e_vanish <= 1e-10
    @test e_Dind_consistent <= 1e-9
    @test n_Dind_bad == n_inconsistent
    @test n_gamma_split > n_Dind_bad
    @test n_copies_split >= n_Dind_bad
    @test e_ones <= 1e-12
    @test e_path <= 1e-9
    @test covec_mismatch == 0
    # L2.2(d): squares_3x3 has k = 4 and dim K = 0; the sign of K0(1,0) tells J from -J.
    @test sq33_k == 4
    @test sq33_dimK == 0
    @test sq33_maxD >= 0.0
    @test sq33_maxD <= 1e-12
    @test abs(sq33_trK) <= 1e-10
    @test abs(sq33_detK - 1.0) <= 1e-10
    @test sq33_K21 > 0
end

# C++: 1.185 s (2000 draws).  Default: the first 400 draws; the full 2000 under
# KIRIGAMI_FULL_DERIVATIONS=1 (with the C++ tally 1999 built, rank(D) 0/0/1999).
@testset "L2 counterexample search: 2000 random (torus, sigma) draws" begin
    ndraw = FULL ? 2000 : 400
    drawn = 0; built = 0; fail_eq = 0; fail_parity = 0
    hist = zeros(Int, 3)
    for i in 0:ndraw-1
        rng = K.MT19937(880011 + 65537 * i)
        n = K.uniform_int(rng, 6, 16)
        drawn += 1
        P = l2_voronoi_pattern(90000 + i, n, 10.0, rng)
        @test same_pattern(P, DERIV_L2["ce2000"][i+1])
        P.ok || continue
        S = l2_prepare(P)
        S.ok || continue
        A = l2_achievable(S)
        built += 1
        hist[min(A.rankD, 2)+1] += 1
        A.dimK != 2 * A.rankD && (fail_eq += 1)
        (A.dimK != 0 && A.dimK != 2 && A.dimK != 4) && (fail_parity += 1)
    end
    @printf("  L2 counterexample search: %d drawn, %d built;  rank(D) 0/1/2 = %d/%d/%d;  rank(A) != 2 rank(D): %d;  dim K odd: %d\n",
            drawn, built, hist[1], hist[2], hist[3], fail_eq, fail_parity)
    record("L2-j", "$(ndraw)-draw counterexample search for rank(A) == 2 rank(D)", built, fail_eq, 0.0)
    if FULL
        @test built >= 500
        @test (built, hist[3]) == (1999, 1999)
    else
        @test built >= 300
        @test hist[3] == built
    end
    @test fail_eq == 0
    @test fail_parity == 0
end

print_report()
