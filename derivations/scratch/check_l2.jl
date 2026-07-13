# check_l2.jl -- Deriver-L, mission 2 / WP1.  Tests Lemmas L2.1 and L2.2 of
# derivations/lemmas.md: the periodic dimension formula dim K = 2 rank(D), and the
# explicit factorisation of the generators M_j through the period-potential row d_i.
# Port of check_l2.cpp (same population, same seeds, same printed lines).
#
# Run (from repo root):
#   julia --project=Kirigami derivations/scratch/check_l2.jl [NEXTRA=100]
#
# The population, the quotient/super-patch pipeline and `achievable()` are copied
# VERBATIM from code/apps/kill_k7.cpp (lines 36-243) so that the numbers compare
# directly with results/kill/k7/k7_main.csv.  Nothing in Kirigami/ is modified.
#
# What is predicted (lemmas.md L2.1).  Write d_i in R^2 for row i of D (the pair of
# period-potential increments of the i-th null vector phi_i) and P0 = T.  Then
#
#     M_{2i}   =  2 e_y d_i^T P0^{-1}          (phi_i applied to the x coordinate)
#     M_{2i+1} = -2 e_x d_i^T P0^{-1}          (phi_i applied to the y coordinate)
#
# equivalently  M_{2i} T = 2 e_y d_i^T  and  M_{2i+1} T = -2 e_x d_i^T.  The span of
# {M_j} is therefore the image of row(D) (+) row(D) under an injective linear map, so
# dim K = rank(A) = 2 rank(D).
#
# Depends on Kirigami/src/method/periodic_jacobian.jl (build_quotient, quotient_system,
# build_super, set_super_positions!, cell_face_area_sum, periodic_jacobian, fk_period_matrix,
# quotient_dual, quotient_sigma, make_tiling_pattern), deploy_basis.jl (deploy_basis) and
# design.jl (median_edge_length) -- C++ names, per PORTING.md.

using Kirigami
using LinearAlgebra
using Printf
using StaticArrays
const K = Kirigami
const Vec2 = K.Vec2
const Mat2 = SMatrix{2,2,Float64,4}

# Local mirror of method::PeriodicPattern for the Voronoi tori (make_tiling_pattern returns
# the package's own type; `prepare` reads only .cell .T .ok .name .err of either).
mutable struct Pattern
    name::String
    family::String
    cell::K.Mesh
    T::Mat2
    ok::Bool
    err::String
end
Pattern() = Pattern("", "", K.Mesh(), Mat2(1, 0, 0, 1), false, "")

# Periodic Voronoi on the torus [0,L)^2: the cells of n random sites, computed
# against the 3x3 replicated site set, form a fundamental domain.
function make_voronoi_pattern(inst::Int, nsites::Int, L::Float64, rng::K.MT19937)
    P = Pattern()
    P.family = "voronoi_torus"
    P.name = "voronoi_torus_" * string(inst) * "_n" * string(nsites)
    site = Vec2[]
    guard = 0
    while length(site) < nsites && (guard += 1) <= 100000
        p = Vec2(K.uniform_real(rng, 0.0, L), K.uniform_real(rng, 0.0, L))
        ok = true
        for s in site
            d = p - s
            d = Vec2(d[1] - L * round(d[1] / L), d[2] - L * round(d[2] / L))
            norm(d) < 0.15 * L / sqrt(Float64(nsites)) && (ok = false)
        end
        ok && push!(site, p)
    end
    if length(site) < nsites
        P.err = "site rejection failed"
        return P
    end
    rep = Vec2[]
    for i in -1:1, j in -1:1, s in site
        push!(rep, s + Vec2(i * L, j * L))
    end
    polys = Vector{Vec2}[]
    for i in 1:nsites
        pi_ = site[i]
        cellp = [pi_ + Vec2(-L, -L), pi_ + Vec2(L, -L), pi_ + Vec2(L, L), pi_ + Vec2(-L, L)]
        for pj in rep
            norm(pj - pi_) < 1e-12 && continue
            nvec = pj - pi_
            off = dot(nvec, 0.5 * (pi_ + pj))
            out = Vec2[]
            for k in eachindex(cellp)
                A = cellp[k]
                B = cellp[k % length(cellp) + 1]
                da = dot(nvec, A) - off; db = dot(nvec, B) - off
                da <= 0 && push!(out, A)
                ((da < 0 && db > 0) || (da > 0 && db < 0)) && push!(out, A + (B - A) * (da / (da - db)))
            end
            cellp = out
            length(cellp) < 3 && break
        end
        if length(cellp) < 3
            P.err = "empty voronoi cell"
            return P
        end
        push!(polys, cellp)
    end
    cell = K.Mesh()
    try
        cell = K.mesh_from_polygons(polys, 1e-7)
    catch e
        P.err = "weld: " * sprint(showerror, e)
        return P
    end
    T = Mat2(L, 0, 0, L)   # column-major: [L 0; 0 L]
    cell.sigma = fill(-1, K.n_faces(cell))
    q0 = K.build_quotient(cell, T)
    if !q0.ok
        P.err = "quotient(probe): " * q0.err
        return P
    end
    cell.sigma = K.quotient_sigma(K.n_faces(cell), K.quotient_dual(q0), rng)
    P.cell = cell
    P.T = T
    P.ok = true
    return P
end

function population()
    out = Any[]
    fams = ["squares", "triangles", "hexagons", "kagome", "snub_square", "trunc_square_488", "t3_4_3_12"]
    sizes = [(2, 2), (3, 2), (3, 3)]
    for f in fams, s in sizes
        rng = K.MT19937(20260904 + 7919 * length(out))
        push!(out, K.make_tiling_pattern(f, s[1], s[2], rng))
    end
    ns = [20, 28, 36, 45, 55, 70, 85, 100, 120, 140, 170, 200]
    for i in 0:11
        rng = K.MT19937(9100001 + 104729 * i)
        push!(out, make_voronoi_pattern(i, ns[i + 1], 10.0, rng))
    end
    return out
end

# ------------------------------------------------------------------ per-pattern state

mutable struct PState
    ok::Bool
    err::String
    q::Any             # Quotient
    sp::Any            # SuperPatch
    cut::Any           # CutStructure
    X0::Matrix{Float64}   # nq x 2
    Phi::Matrix{Float64}  # nq x k
    k::Int                # dim null
    med_edge::Float64
end
PState() = PState(false, "", nothing, nothing, nothing, zeros(0, 2), zeros(0, 0), 0, 1.0)

function shape_point(S::PState, t::Vector{Float64})
    T = Matrix{Float64}(undef, S.k, 2)
    for i in 1:S.k
        T[i, 1] = t[2i - 1]
        T[i, 2] = t[2i]
    end
    X = copy(S.X0)
    S.k > 0 && (X += S.Phi * T)
    return K.matrix_to_points(X)
end

function prepare(P)
    S = PState()
    S.q = K.build_quotient(P.cell, P.T)
    if !S.q.ok
        S.err = "quotient: " * S.q.err
        return S
    end
    sys = K.quotient_system(S.q, S.q.Xq)
    sr = K.solve_system(sys, S.q.Xq)
    if !sr.projection_ok
        S.err = "projection failed"
        return S
    end
    S.X0 = sr.X0
    S.Phi = sr.Phi
    S.k = sr.dim_null
    S.sp = K.build_super(S.q, 1)
    S.cut = K.make_cut(S.sp.mesh)
    S.med_edge = K.median_edge_length(P.cell)
    S.ok = true
    return S
end

# K at a shape-space point, from the closed-form deploy basis of the super patch.
function K_at(S::PState, Xq::Vector{Vec2})
    K.set_super_positions!(S.sp, S.q, Xq)
    B = K.deploy_basis(S.cut, S.sp.mesh.X)
    return K.periodic_jacobian(S.sp, S.q, S.cut, B)
end

mutable struct AchievableSet
    K0::Matrix{Float64}
    M::Vector{Matrix{Float64}}   # 2k generators
    A::Matrix{Float64}           # 4 x 2k, columns = vec(M_j)
    dimK::Int                    # rank(A)
    D::Matrix{Float64}           # k x 2 period-potential increments
    rankD::Int
end

function svd_rank(Mx::AbstractMatrix{Float64}, rel::Float64)
    isempty(Mx) && return 0
    sv = svdvals(Mx)
    isempty(sv) && return 0
    tol = max(rel * sv[1], 1e-13)
    return count(x -> x > tol, sv)
end

# The affine map t -> K, plus the achievable-set rank and the period-potential matrix D.
function achievable(S::PState)
    z = zeros(2 * S.k)
    K0 = Matrix{Float64}(K_at(S, shape_point(S, z)).K)
    A = AchievableSet(K0, Matrix{Float64}[], zeros(4, 2 * S.k), 0, zeros(max(S.k, 1), 2), 0)
    sc = S.med_edge
    for j in 1:2*S.k
        e = zeros(2 * S.k)
        e[j] = sc
        M = (Matrix{Float64}(K_at(S, shape_point(S, e)).K) - A.K0) / sc
        push!(A.M, M)
        A.A[:, j] = [M[1, 1], M[2, 1], M[1, 2], M[2, 2]]
        if isodd(j)   # the (phi_i, e_1) generator (C++ j even) carries d_i in row 2 of M*T
            MT = M * S.q.T
            A.D[(j + 1) ÷ 2, 1] = 0.5 * MT[2, 1]
            A.D[(j + 1) ÷ 2, 2] = 0.5 * MT[2, 2]
        end
    end
    if 2 * S.k > 0
        A.dimK = svd_rank(A.A, 1e-9)
        A.rankD = svd_rank(A.D, 1e-9)
    end
    return A
end

rank_svd(Mx::AbstractMatrix{Float64}, rel::Float64 = 1e-10) = svd_rank(Mx, rel)

mutable struct Row
    name::String
    k::Int
    dimK::Int
    rankD::Int
    rankD_svd::Int
    err_fact::Float64    # max |M_j - (formula from row j/2 of D)|
    err_xrow::Float64    # max |row 0 of (M_2i T)|  and |row 1 of (M_2i+1 T)|
    rankD_cov::Int
    have_cov::Bool
    consistency::Float64   # R2 (D-L2-5): H2 on the SUPER PATCH
    p0_err::Float64
    eq::Bool
end
Row() = Row("", 0, 0, 0, 0, 0.0, 0.0, -1, false, 0.0, 0.0, false)

function main()
    NEXTRA = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 100
    rows = Row[]
    n_eq = 0; n_fail = 0; n_cov_ok = 0; n_cov_bad = 0; n_inconsistent = 0
    worst_consistency = 0.0
    worst_fact = 0.0; worst_xrow = 0.0
    rankD_hist = Dict{Int,Int}()
    failures = String[]

    pats = population()          # the 33 K7 patterns
    # ... plus NEXTRA random Voronoi tori, fresh sigma seeds
    for i in 0:NEXTRA-1
        ns = [12, 16, 20, 24, 30, 36, 45]
        rng = K.MT19937(7000000 + 65537 * i)
        P = make_voronoi_pattern(1000 + i, ns[i % 7 + 1], 10.0, rng)
        P.ok && push!(pats, P)
    end

    for P in pats
        P.ok || continue
        S = prepare(P)
        S.ok || continue
        A = achievable(S)
        # R2 (D-L2-5): the sharp test of H2 on the patch the period is MEASURED on is
        # PeriodicJac::consistency -- the spread of the measured period over (face, corner).
        # achievable() never reads it, so Check L2 must.
        PJ0 = K_at(S, K.matrix_to_points(S.X0))
        R = Row()
        R.consistency = PJ0.consistency
        R.p0_err = PJ0.p0_err
        R.name = P.name
        R.k = S.k
        R.dimK = A.dimK
        R.rankD = A.rankD
        R.rankD_svd = rank_svd(A.D)
        rA = rank_svd(A.A)
        # the factorisation
        T = Matrix{Float64}(S.q.T)
        Ti = inv(T)
        for i in 1:S.k
            d = A.D[i, :]
            Mx = zeros(2, 2); My = zeros(2, 2)
            Mx[2, :] = 2.0 * d             # M_{2i}  T   = 2 e_y d^T
            My[1, :] = -2.0 * d            # M_{2i+1} T = -2 e_x d^T
            Mx_pred = Mx * Ti; My_pred = My * Ti
            R.err_fact = max(R.err_fact, maximum(abs.(A.M[2i - 1] - Mx_pred)))
            R.err_fact = max(R.err_fact, maximum(abs.(A.M[2i] - My_pred)))
            MxT = A.M[2i - 1] * T; MyT = A.M[2i] * T
            R.err_xrow = max(R.err_xrow, maximum(abs.(MxT[1, :])))
            R.err_xrow = max(R.err_xrow, maximum(abs.(MyT[2, :])))
        end
        # L2.2: rank(D) = rank([L; delta_h; delta_v]) - rank(L), where delta_tau is the
        # period-potential increment read as a covector on ALL of R^nq (not just the null
        # space).  Finite differences in the quotient x-coordinates; only for small cells.
        if S.q.nq <= 60 && !isempty(S.q.L)
            nq = S.q.nq
            G = zeros(2, nq)
            hstep = 1e-4 * S.med_edge
            for v in 1:nq
                Xv = K.matrix_to_points(S.X0)
                Xv[v] = Vec2(Xv[v][1] + hstep, Xv[v][2])
                Mv = (Matrix{Float64}(K_at(S, Xv).K) - A.K0) / hstep
                MvT = Mv * T
                G[1, v] = 0.5 * MvT[2, 1]
                G[2, v] = 0.5 * MvT[2, 2]
            end
            LG = vcat(Matrix{Float64}(S.q.L), G)
            rL = rank_svd(Matrix{Float64}(S.q.L)); rLG = rank_svd(LG)
            R.rankD_cov = rLG - rL
            R.have_cov = true
        end
        if R.rankD_svd == 0 && S.k > 0
            # R2 (D-L2-3): tr and det do NOT separate K0 = +J from -J; print K0(1,0) too.
            @printf("    [rank(D) = 0 with k = %d]  max|D| = %.3e   tr K0 = %.6f  det K0 = %.6f  K0(1,0) = %+.6f\n",
                    S.k, maximum(abs.(A.D)), tr(A.K0), det(A.K0), A.K0[2, 1])
        end
        if R.consistency > 1e-9
            n_inconsistent += 1
            @printf("    [H2 FAILS on the super patch] %-24s consistency = %.3e  p0_err = %.3e   -> EXCLUDED from the L2.1 verification\n",
                    R.name, R.consistency, R.p0_err)
        end
        R.eq = (rA == 2 * R.rankD_svd)
        if R.eq
            n_eq += 1
        else
            n_fail += 1; push!(failures, R.name)
        end
        R.consistency <= 1e-9 && (worst_consistency = max(worst_consistency, R.consistency))
        worst_fact = max(worst_fact, R.err_fact)
        worst_xrow = max(worst_xrow, R.err_xrow)
        rankD_hist[R.rankD_svd] = get(rankD_hist, R.rankD_svd, 0) + 1
        push!(rows, R)
        @printf("%-28s k=%3d  dimK=%d rank(A)=%d rank(D)=%d  2rD=%d %s  fact_err=%.2e  xrow_err=%.2e\n",
                R.name, R.k, R.dimK, rA, R.rankD_svd, 2 * R.rankD_svd,
                R.eq ? "OK " : "FAIL", R.err_fact, R.err_xrow)
        if R.have_cov
            if R.rankD_cov == R.rankD_svd
                n_cov_ok += 1
            else
                n_cov_bad += 1
                @printf("    [covector rank mismatch] %s: rank([L;delta])-rank(L) = %d, rank(D) = %d\n",
                        R.name, R.rankD_cov, R.rankD_svd)
            end
        end
    end

    @printf("\ncheck_l2 -- patterns evaluated: %d\n", length(rows))
    @printf("  rank(A) == 2 rank(D):  %d / %d    failures: %d\n", n_eq, length(rows), n_fail)
    for f in failures
        @printf("    FAIL: %s\n", f)
    end
    @printf("  rank(D) distribution:")
    for k in sort!(collect(keys(rankD_hist)))
        @printf("   %d -> %d", k, rankD_hist[k])
    end
    @printf("\n  max |M_j - formula(D)|            : %.3e\n", worst_fact)
    @printf("  max |vanishing row of M_j T|      : %.3e\n", worst_xrow)
    @printf("  H2 on the super patch (PeriodicJac::consistency > 1e-9): %d of %d patterns  [EXCLUDED]\n",
            n_inconsistent, length(rows))
    @printf("  worst consistency among the retained patterns: %.3e\n", worst_consistency)
    @printf("  L2.2 covector test  rank([L;delta_h;delta_v]) - rank(L) == rank(D):  %d ok, %d mismatch\n",
            n_cov_ok, n_cov_bad)
    @printf("\n")
end

main()
