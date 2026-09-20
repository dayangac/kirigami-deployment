# K7 -- achievable periodic deployment Jacobians (R4) and the closed-form second
# closed angle (B20).
#
# C1  J(theta) = P_theta P_0^{-1} = cos(theta/2) I + sin(theta/2) K, K constant in
#     theta and linear in the vertex positions; dim of the achievable set K.
# C2  K a similarity  <=>  J(theta) conformal for every theta.
# C3  designable Poisson family: hit a target K*, then maximise certified Theta_max
#     over the remaining freedom.
# C4  total hole area A(theta) = q(cos theta - 1) + r sin theta, second closed angle
#     theta_c = 2 atan2(r, q); periodic form theta_c = 2 atan2(tr K, 1 - det K).
#
# Everything periodic is done on the QUOTIENT (method/periodic_jacobian.jl), with a 3x3
# super patch as the forward-kinematics witness.
#
#   julia --project=Kirigami Kirigami/apps/exp_k7_periodic_jacobian.jl [--stage run|pop|bounded|perdbg|
#         c4sweep|nu|c3] [--out DIR] [--shard i] [--nshard k] [--limit K]
#
# The population (21 tiling super-cells + 12 torus Voronoi patterns) is rebuilt from its
# seeds through the bit-exact MT19937 every run (it is not in the corpus).
# Outputs go to results/experiments/k7/. `--limit K`
# processes only the first K patterns of the population (after sharding).
include(joinpath(@__DIR__, "exp_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const Mat2 = K.Mat2

# 12 significant digits (the CSV number format).
g12(v::Real) = @sprintf("%.12g", Float64(v))

# ------------------------------------------------------------------ population

# Periodic Voronoi on the torus [0,L)^2: the cells of n random sites, computed
# against the 3x3 replicated site set, form a fundamental domain.
function make_voronoi_pattern(inst::Int, nsites::Int, L::Float64, rng::K.MT19937)
    P = K.PeriodicPattern()
    P.family = "voronoi_torus"
    P.name = "voronoi_torus_" * string(inst) * "_n" * string(nsites)
    site = Vec2[]
    guard = 0
    while length(site) < nsites && guard < 100000
        guard += 1
        x = K.uniform_real(rng, 0.0, L)
        y = K.uniform_real(rng, 0.0, L)
        p = Vec2(x, y)
        ok = true
        for s_ in site
            d = p - s_
            d = Vec2(d[1] - L * round(d[1] / L, RoundNearestTiesAway),
                     d[2] - L * round(d[2] / L, RoundNearestTiesAway))
            K._norm2(d) < 0.15 * L / sqrt(Float64(nsites)) && (ok = false)
        end
        ok && push!(site, p)
    end
    if length(site) < nsites
        P.err = "site rejection failed"
        return P
    end
    rep = Vec2[]
    for i in -1:1, j in -1:1, s_ in site
        push!(rep, s_ + Vec2(i * L, j * L))
    end
    polys = Vector{Vec2}[]
    for i in 1:nsites
        pi_ = site[i]
        cellp = [pi_ + Vec2(-L, -L), pi_ + Vec2(L, -L), pi_ + Vec2(L, L), pi_ + Vec2(-L, L)]
        for pj in rep
            K._norm2(pj - pi_) < 1e-12 && continue
            nvec = pj - pi_
            off = K._dot2(nvec, 0.5 * (pi_ + pj))
            out = Vec2[]
            nc = length(cellp)
            for k in 1:nc
                A = cellp[k]
                B = cellp[mod1(k + 1, nc)]
                da = K._dot2(nvec, A) - off
                db = K._dot2(nvec, B) - off
                da <= 0 && push!(out, A)
                if (da < 0 && db > 0) || (da > 0 && db < 0)
                    push!(out, A + (B - A) * (da / (da - db)))
                end
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
    cell = try
        K.mesh_from_polygons(polys, 1e-7)
    catch e
        P.err = "weld: " * sprint(showerror, e)
        return P
    end
    T = Mat2(L, 0.0, 0.0, L)
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

function population_k7()
    out = K.PeriodicPattern[]
    fams = ("squares", "triangles", "hexagons", "kagome", "snub_square", "trunc_square_488", "t3_4_3_12")
    sizes = ((2, 2), (3, 2), (3, 3))
    for f in fams, sz in sizes
        rng = K.MT19937(UInt32(20260904) + UInt32(7919) * UInt32(length(out)))
        push!(out, K.make_tiling_pattern(f, sz[1], sz[2], rng))
    end
    ns = (20, 28, 36, 45, 55, 70, 85, 100, 120, 140, 170, 200)
    for i in 0:11
        rng = K.MT19937(UInt32(9100001) + UInt32(104729) * UInt32(i))
        push!(out, make_voronoi_pattern(i, ns[i + 1], 10.0, rng))
    end
    return out
end

# ------------------------------------------------------------------ per-pattern state

mutable struct PState
    ok::Bool
    err::String
    q::K.Quotient
    sp::K.SuperPatch
    cut::K.CutStructure
    X0::Matrix{Float64}    # nq x 2
    Phi::Matrix{Float64}   # nq x k
    k::Int                 # dim null
    med_edge::Float64
    cell_face_area::Float64  # sum of the |F| face areas of one cell (theta-independent)
end

# t (length 2k, row-major layout t(2i) = T(i,0), t(2i+1) = T(i,1)) -> quotient positions
function shape_point(S::PState, t::AbstractVector{Float64})
    X = copy(S.X0)
    if S.k > 0
        T = Matrix{Float64}(undef, S.k, 2)
        for i in 1:S.k
            T[i, 1] = t[2i - 1]
            T[i, 2] = t[2i]
        end
        X += S.Phi * T
    end
    return K.matrix_to_points(X)
end

function prepare(P::K.PeriodicPattern)
    q = K.build_quotient(P.cell, P.T)
    q.ok || return (nothing, "quotient: " * q.err)
    sys = K.quotient_system(q, q.Xq)
    sr = K.solve_system(sys, q.Xq)
    sr.projection_ok || return (nothing, "projection failed")
    sp = K.build_super(q, 1)
    cut = K.make_cut(sp.mesh)
    S = PState(true, "", q, sp, cut, sr.X0, sr.Phi, sr.dim_null, median_edge_length(P.cell), 0.0)
    return (S, "")
end

# The quotient positions are REDUCED into the origin cell, so face geometry must be
# read off the super patch (which restores the lattice offsets), never off `cell`.
function cell_area_sum(S::PState, Xq::Vector{Vec2})
    K.set_super_positions!(S.sp, S.q, Xq)
    inv = 0
    for f in 1:S.sp.nfc
        K.face_signed_area(S.sp.mesh, K.face_index(S.sp, 0, 0, f)) <= 0 && (inv += 1)
    end
    return K.cell_face_area_sum(S.sp), inv
end

# K at a shape-space point, from the closed-form deploy basis of the super patch.
function K_at(S::PState, Xq::Vector{Vec2})
    K.set_super_positions!(S.sp, S.q, Xq)
    B = K.deploy_basis(S.cut, S.sp.mesh.X)
    return K.periodic_jacobian(S.sp, S.q, S.cut, B)
end

# P_theta measured by an INDEPENDENT forward-kinematics call (deploy(), not the basis).
P_fk(S::PState, theta::Float64) = K.fk_period_matrix(S.sp, S.q, S.cut, theta)

# rank at the SVD threshold 1e-9 * s_max (floored at 1e-13)
function sv_rank(A::AbstractMatrix{Float64})
    (size(A, 1) == 0 || size(A, 2) == 0) && return 0
    sv = svdvals(A)
    tol = 1e-9 * (isempty(sv) ? 1.0 : sv[1])
    return count(v -> v > max(tol, 1e-13), sv)
end

# The affine map t -> K, plus the achievable-set rank and the period-potential matrix D.
function achievable(S::PState)
    A = K.AchievableSet()
    z = zeros(2 * S.k)
    A.K0 = K_at(S, shape_point(S, z)).K
    A.A = zeros(4, 2 * S.k)
    A.D = zeros(max(S.k, 1), 2)
    sc = S.med_edge
    for j in 0:(2 * S.k - 1)
        e = zeros(2 * S.k)
        e[j + 1] = sc
        M = (K_at(S, shape_point(S, e)).K - A.K0) / sc
        push!(A.M, M)
        A.A[:, j + 1] = [M[1, 1], M[2, 1], M[1, 2], M[2, 2]]
        if j % 2 == 0   # the (phi_i, e_1) generator carries d_i in row 2 of M*T
            MT = M * S.q.T
            A.D[j ÷ 2 + 1, 1] = 0.5 * MT[2, 1]
            A.D[j ÷ 2 + 1, 2] = 0.5 * MT[2, 2]
        end
    end
    if 2 * S.k > 0
        A.dimK = sv_rank(A.A)
        A.rankD = sv_rank(A.D)
    end
    return A
end

function K_of_t(A::K.AchievableSet, t::AbstractVector{Float64})
    Km = A.K0
    for j in eachindex(A.M)
        Km += t[j] * A.M[j]
    end
    return Km
end

# the minimum-norm least-squares solution
cod_solve(A::AbstractMatrix{Float64}, b::AbstractVector{Float64}) = pinv(A) * b

# the free subspace of the target-hitting affine set: the right singular vectors of A
# past its rank (the last nfree columns of V of the full SVD)
function free_subspace(A::AbstractMatrix{Float64}, nfree::Int)
    n = size(A, 2)
    nfree > 0 || return zeros(n, 0)
    V = svd(A; full = true).V
    return Matrix(V[:, (n - nfree + 1):n])
end

# ------------------------------------------------------------------ certificate

mutable struct Cert
    pos::Bool; noovl::Bool; valid::Bool
    eps_max::Float64      # largest certified eps  (= first admissible deflated root)
    theta_exact::Float64  # T4.2" exact Theta_max on the super patch
end

function certify(c::K.CutStructure, X::Vector{Vec2}, want_exact::Bool)
    r = Cert(false, false, false, 0.0, 0.0)
    # Pass 1 finds the first admissible deflated root, i.e. the largest eps for which
    # NOROOT can hold. NOOVERLAP is tested at eps/2, so it must be re-evaluated INSIDE
    # that window.
    v1 = K.validity_certificate(c, X, Float64(pi), 1e-12)
    cap = v1.first_root > 0 ? v1.first_root : Float64(pi)
    vc = K.validity_certificate(c, X, max(1e-12, cap * (1 - 1e-9)), 1e-12)
    r.pos = vc.pos
    r.noovl = vc.nooverlap
    r.eps_max = (vc.pos && vc.nooverlap && vc.noroot) ? cap * (1 - 1e-9) : 0.0
    r.valid = r.eps_max > 0
    if want_exact
        B = K.deploy_basis(c, X)
        sd = K.swept_discs(c, B)
        pairs = K.candidate_pairs(c, sd, Float64(pi), true)
        r.theta_exact = K.exact_theta_max_overlap(c, B, pairs).theta_max
    end
    return r
end

# ------------------------------------------------------------------ hole harmonics

dt2(a::Vec2, b::Vec2) = a[1] * b[2] - a[2] * b[1]

# A(theta) = p + q cos + r sin for the sum of the shoelace areas of the given
# M'-vertex cycles, from the closed form of derivations/core.md T3.2.
function cycles_area_harmonic(B::K.DeployBasis, cyc::Vector{Vector{Int}})
    p = q = r = 0.0
    for cy in cyc, i in eachindex(cy)
        a = cy[i]; b = cy[mod1(i + 1, length(cy))]
        Ca = K.basis_c(B, a); Cb = K.basis_c(B, b); Sa = K.basis_s(B, a); Sb = K.basis_s(B, b)
        p += 0.25 * (dt2(Ca, Cb) + dt2(Sa, Sb))
        q += 0.25 * (dt2(Ca, Cb) - dt2(Sa, Sb))
        r += 0.25 * (dt2(Ca, Sb) + dt2(Sa, Cb))
    end
    return K.Harmonic(p, q, r)
end

function cycles_area(Y::Vector{Vec2}, cyc::Vector{Vector{Int}})
    a = 0.0
    for cy in cyc, i in eachindex(cy)
        p = Y[cy[i]]
        q = Y[cy[mod1(i + 1, length(cy))]]
        a += 0.5 * (p[1] * q[2] - p[2] * q[1])
    end
    return a
end

# second root of p + q cos + r sin with p = -q:  theta_c = 2 atan2(r, q), in (0, 2pi)
function second_closed_angle(q::Float64, r::Float64)
    t = 2.0 * atan(r, q)
    while t <= 1e-12
        t += 2 * pi
    end
    while t > 2 * pi
        t -= 2 * pi
    end
    return t
end

# The 0+ margin objective of K6 / method/zero_plus.jl, evaluated on the super patch and
# restricted to a linear subspace of the periodic shape space.
softplus(x::Float64) = x > 30 ? x : log1p(exp(x))

function zp_objective(S::PState, Xq::Vector{Vec2})
    K.set_super_positions!(S.sp, S.q, Xq)
    sq = max(1e-12, 0.05 * S.med_edge * S.med_edge)
    f = 0.0
    for v in K.zero_plus_q(S.cut, S.sp.mesh.X)
        f += softplus(-v / sq)
    end
    for fi in 1:K.n_faces(S.sp.mesh)
        f += softplus(-K.face_signed_area(S.sp.mesh, fi) / sq)
    end
    return f
end

# nu_pred / nu_meas for one direction angle
function nu_pair(Kb::Mat2, Jm::Mat2, th::Float64, u::Vec2)
    up = Vec2(-u[2], u[1])
    J = K.jac_at(Kb, th)
    pr = -(norm(J * up) - 1) / (norm(J * u) - 1)
    me = -(norm(Jm * up) - 1) / (norm(Jm * u) - 1)
    return pr, me
end

# the conformal design: t with K(t) a similarity (min-norm least squares)
function conformal_t(AS::K.AchievableSet, k::Int)
    C = zeros(2, 2 * k)
    for j in 1:(2 * k)
        C[1, j] = AS.M[j][1, 1] - AS.M[j][2, 2]
        C[2, j] = AS.M[j][1, 2] + AS.M[j][2, 1]
    end
    rhs = [-(AS.K0[1, 1] - AS.K0[2, 2]), -(AS.K0[1, 2] + AS.K0[2, 1])]
    return C, rhs, cod_solve(C, rhs)
end

k_target_rhs(Kstar::Mat2, K0::Mat2) = [Kstar[1, 1] - K0[1, 1], Kstar[2, 1] - K0[2, 1],
                                       Kstar[1, 2] - K0[1, 2], Kstar[2, 2] - K0[2, 2]]

function gaussian_vec(n::Int, scale::Float64, G::K.NormalDist, rng::K.MT19937)
    v = Vector{Float64}(undef, n)
    for i in 1:n
        v[i] = scale * K.normal(G, rng)
    end
    return v
end

# ------------------------------------------------------------------ stages

function stage_bounded(outdir::String)
    csv = open(joinpath(outdir, "k7_c4_bounded.csv"), "w")
    print(csv, "case,H_traced,p,q,r,p_plus_q,theta_c,theta_c_meas,abs_err,theta_max,",
          "closed_before_contact,area_fit_maxerr\n")
    for rc in reference_cases()
        rc.periodic && continue
        m = rc.mesh
        K.build_topology!(m)
        c = K.make_cut(m)
        hs = K.holes_partition(c)
        X = copy(m.X)
        if !K.deployable(K.hole_residuals(c, X, hs), 1e-9)
            sys = K.assemble_system(c, hs, m.X, K.Fixed)
            sr = K.solve_system(sys, m.X)
            sr.projection_ok || continue
            X = K.matrix_to_points(sr.X0)
        end
        B = K.deploy_basis(c, X)
        sd = K.swept_discs(c, B)
        pairs = K.candidate_pairs(c, sd, Float64(pi), true)
        tmax = K.exact_theta_max_overlap(c, B, pairs).theta_max
        th_tr = 0.5 * max(1e-3, min(tmax, Float64(pi)))
        cyc = try
            K.holes_geometric_cycles(c, K.deploy(c, X, th_tr).Y)
        catch
            continue
        end
        isempty(cyc) && continue
        h = cycles_area_harmonic(B, cyc)
        fit_err = 0.0
        scale = abs(h.p) + hypot(h.q, h.r)
        for i in 1:40
            th = pi * i / 40.0
            meas = cycles_area(K.basis_eval(B, th), cyc)
            fit_err = max(fit_err, abs(meas - K.harmonic_eval(h, th)))
        end
        thc = second_closed_angle(h.q, h.r)
        # measured: bisect the traced-cycle area for its second zero
        A(th) = cycles_area(K.basis_eval(B, th), cyc)
        lo = 1e-6; hi = -1.0
        G = 20000
        prev = A(lo)
        for i in 1:G
            th = 2 * pi * i / G
            v = A(th)
            if prev == 0 || (v < 0) != (prev < 0)
                hi = th
                break
            end
            lo = th
            prev = v
        end
        thm = -1.0
        if hi > 0
            for _ in 1:80
                mid = 0.5 * (lo + hi)
                (A(mid) < 0) != (A(lo) < 0) ? (hi = mid) : (lo = mid)
            end
            thm = 0.5 * (lo + hi)
        end
        print(csv, rc.name, ",", length(cyc), ",", fmt_g(h.p), ",", fmt_g(h.q), ",", fmt_g(h.r), ",",
              fmt_g((h.p + h.q) / max(1e-300, scale)), ",", fmt_g(thc), ",", fmt_g(thm), ",",
              fmt_g(thm > 0 ? abs(thc - thm) : -1), ",", fmt_g(tmax), ",",
              (thc <= tmax + 1e-9 ? 1 : 0), ",", fmt_g(fit_err / max(1e-300, scale)), "\n")
        @printf("%-24s H=%2d thc=%.9f meas=%.9f err=%.2e tmax=%.6f p+q/s=%.2e\n",
                rc.name, length(cyc), thc, thm, thm > 0 ? abs(thc - thm) : -1.0, tmax, (h.p + h.q) / scale)
    end
    close(csv)
end

function stage_perdbg(outdir::String, pats, shard::Int, nshard::Int, limit::Int)
    csv = open(joinpath(outdir, "k7_periodicity_$(shard).csv"), "w")
    print(csv, "name,family,dim_null,n_components,quot_res,cons,p0err,fk_spread_1p1,mean_vs_closed_form,",
          "deploy_mismatch\n")
    nrun = 0
    for (pi_, P) in enumerate(pats)
        (pi_ - 1) % nshard != shard && continue
        P.ok || continue
        nrun >= limit && break
        nrun += 1
        S, _ = prepare(P)
        S === nothing && continue
        X0 = K.matrix_to_points(S.X0)
        qres = K.quotient_residual(S.q, X0)
        K.set_super_positions!(S.sp, S.q, X0)
        J = K.periodic_jacobian(S.sp, S.q, S.cut, K.deploy_basis(S.cut, S.sp.mesh.X))
        th = 1.1
        Pm, spr = P_fk(S, th)
        Pc = cos(th / 2) * J.P0 + sin(th / 2) * J.Q
        Dp = K.deploy(S.cut, S.sp.mesh.X, th)
        ncomp = K.count_components(S.cut)
        print(csv, P.name, ",", P.family, ",", S.k, ",", ncomp, ",", g12(qres), ",", g12(J.consistency),
              ",", g12(J.p0_err), ",", g12(spr), ",", g12(maximum(abs, Pm - Pc)), ",",
              g12(Dp.max_mismatch), "\n")
        flush(csv)
        @printf("%-24s ncomp=%d qres=%.2e cons=%.3g spread=%.3g\n", P.name, ncomp, qres, J.consistency, spr)
    end
    close(csv)
end

function stage_c4sweep(outdir::String, pats)
    want = ("hexagons_2x2", "t3_4_3_12_2x2", "voronoi_torus_0_n20", "voronoi_torus_2_n36")
    f = open(joinpath(outdir, "k7_c4_sweep.csv"), "w")
    print(f, "name,theta,A_meas,A_model,rel_dev,spread,detP\n")
    for P in pats
        (P.ok && P.name in want) || continue
        S, _ = prepare(P)
        S === nothing && continue
        J0 = K_at(S, K.matrix_to_points(S.X0))
        S.cell_face_area, _ = cell_area_sum(S, K.matrix_to_points(S.X0))
        dP0 = abs(det(J0.P0))
        q = -dP0 * (det(J0.K) - 1) / 2
        r = dP0 * tr(J0.K) / 2
        sc = abs(q) + abs(r) + abs(dP0)
        K.set_super_positions!(S.sp, S.q, K.matrix_to_points(S.X0))
        for i in 1:200
            th = 2 * pi * i / 201.0
            Pm, sp = P_fk(S, th)
            am = abs(det(Pm)) - S.cell_face_area
            ad = q * (cos(th) - 1) + r * sin(th)
            print(f, P.name, ",", g12(th), ",", g12(am), ",", g12(ad), ",", g12(abs(am - ad) / sc), ",",
                  g12(sp), ",", g12(det(Pm)), "\n")
        end
        println("c4sweep: ", P.name)
    end
    close(f)
end

function stage_nu(outdir::String, pats)
    want = ("hexagons_2x2", "snub_square_2x2", "t3_4_3_12_2x2", "voronoi_torus_0_n20")
    for P in pats
        (P.ok && P.name in want) || continue
        S, _ = prepare(P)
        S === nothing && continue
        AS = achievable(S)
        # three designs: X0, the conformal solution, and a random achievable target
        designs = Tuple{String,Vector{Float64}}[("X0", zeros(2 * S.k))]
        if AS.dimK >= 1
            _, _, tc = conformal_t(AS, S.k)
            push!(designs, ("conformal", tc))
            rng = K.MT19937(5150)
            G = K.NormalDist(0.0, 1.0)
            push!(designs, ("random_in_K", gaussian_vec(2 * S.k, 0.25 * S.med_edge, G, rng)))
        end
        f = open(joinpath(outdir, "k7_nu_" * P.name * ".csv"), "w")
        print(f, "design,theta,dir_deg,nu_pred,nu_meas,distortion\n")
        for (dn, t) in designs
            Xd = shape_point(S, t)
            Jd = K_at(S, Xd)
            K.set_super_positions!(S.sp, S.q, Xd)
            for i in 1:60
                th = pi * i / 61.0
                Pm, _ = P_fk(S, th)
                Jm = Pm * inv(Jd.P0)
                for a in 0:3
                    ang = a * pi / 4
                    pr, me = nu_pair(Jd.K, Jm, th, Vec2(cos(ang), sin(ang)))
                    print(f, dn, ",", g12(th), ",", a * 45, ",", g12(pr), ",", g12(me), ",",
                          g12(K.conformal_distortion(K.jac_at(Jd.K, th))), "\n")
                end
            end
        end
        close(f)
        println("nu curves: ", P.name)
    end
end

# max |nu_pred - nu_meas| over 20 angles and two directions at the design t_best
function nu_max_error(S::PState, Kb::Mat2, P0::Mat2, Xb::Vector{Vec2})
    numax = 0.0
    dirs = (Vec2(1.0, 0.0), Vec2(cos(0.7), sin(0.7)))
    K.set_super_positions!(S.sp, S.q, Xb)
    for i in 1:20
        th = pi * i / 21.0
        Pm, _ = P_fk(S, th)
        Jm = Pm * inv(P0)
        for d in dirs
            pr, me = nu_pair(Kb, Jm, th, d / norm(d))
            (isfinite(pr) && isfinite(me)) && (numax = max(numax, abs(pr - me)))
        end
    end
    return numax
end

function stage_c3(outdir::String, pats, shard::Int, nshard::Int, limit::Int)
    csv = open(joinpath(outdir, "k7_c3_$(shard).csv"), "w")
    print(csv, "name,family,dim_null,dimK,nfree,target,hit_err,cert,eps_cert,theta_exact,",
          "nu_maxerr,min_q,min_area,zp_feasible\n")
    nrun = 0
    for (pi_, P) in enumerate(pats)
        pi0 = pi_ - 1   # 0-based population index (part of the seeds)
        pi0 % nshard != shard && continue
        P.ok || continue
        nrun >= limit && break
        nrun += 1
        S, _ = prepare(P)
        S === nothing && continue
        AS = achievable(S)
        AS.dimK < 1 && continue
        J0 = K_at(S, K.matrix_to_points(S.X0))
        rng = K.MT19937(UInt32(31337) + UInt32(613) * UInt32(pi0))
        G = K.NormalDist(0.0, 1.0)
        nfree = 2 * S.k - AS.dimK
        Vf = free_subspace(AS.A, nfree)
        for tgt in 0:1
            Kstar = if tgt == 0
                K_of_t(AS, gaussian_vec(2 * S.k, 0.25 * S.med_edge, G, rng))
            else
                Mat2(1.0, 0.0, 0.0, -0.5)
            end
            tstar = cod_solve(AS.A, k_target_rhs(Kstar, AS.K0))
            hit = maximum(abs, K_of_t(AS, tstar) - Kstar)
            # pattern search on the 0+ margin objective over the free subspace
            s_best = zeros(max(nfree, 1))
            tof(sv) = nfree > 0 ? tstar + Vf * sv[1:nfree] : tstar
            fbest = zp_objective(S, shape_point(S, tof(s_best)))
            for _ in 1:4
                nfree > 0 || break
                sv = gaussian_vec(nfree, 0.3 * S.med_edge, G, rng)
                f = zp_objective(S, shape_point(S, tof(sv)))
                if f < fbest
                    fbest = f
                    s_best = sv
                end
            end
            if nfree > 0
                step = 0.4 * S.med_edge
                it = 0
                while it < 60 && step > 1e-3 * S.med_edge
                    it += 1
                    imp = false
                    for i in 1:nfree, sg in (-1, 1)
                        sv = copy(s_best)
                        sv[i] += sg * step
                        f = zp_objective(S, shape_point(S, tof(sv)))
                        if f < fbest - 1e-12
                            fbest = f
                            s_best = sv
                            imp = true
                        end
                    end
                    imp || (step *= 0.5)
                end
            end
            tb = tof(s_best)
            Xb = shape_point(S, tb)
            K.set_super_positions!(S.sp, S.q, Xb)
            minq = 1e300; mina = 1e300
            for v in K.zero_plus_q(S.cut, S.sp.mesh.X)
                minq = min(minq, v)
            end
            for fi in 1:K.n_faces(S.sp.mesh)
                mina = min(mina, K.face_signed_area(S.sp.mesh, fi))
            end
            ct = certify(S.cut, S.sp.mesh.X, K.n_faces(S.sp.mesh) <= 500)
            numax = -1.0
            if ct.eps_max > 0
                Kb = K_at(S, Xb).K
                numax = nu_max_error(S, Kb, J0.P0, Xb)
            end
            print(csv, P.name, ",", P.family, ",", S.k, ",", AS.dimK, ",", nfree, ",",
                  (tgt == 0 ? "random_in_K" : "diag_1_-0.5"), ",", g12(hit), ",",
                  (ct.eps_max > 0 ? 1 : 0), ",", g12(ct.eps_max), ",", g12(ct.theta_exact), ",",
                  g12(numax), ",", g12(minq), ",", g12(mina), ",",
                  ((minq > 0 && mina > 0) ? 1 : 0), "\n")
            flush(csv)
            @printf("[%2d] %-26s tgt=%d hit=%.1e cert=%d eps=%.5f thmax=%.5f minq=%.2e\n", pi0,
                    P.name, tgt, hit, ct.eps_max > 0 ? 1 : 0, ct.eps_max, ct.theta_exact, minq)
        end
    end
    close(csv)
end

const K7_MAIN_EMPTY = ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"   # 34 empty columns before `err`

function stage_main(outdir::String, pats, shard::Int, nshard::Int, limit::Int)
    csv = open(joinpath(outdir, "k7_main_$(shard).csv"), "w")
    print(csv, "name,family,Fcell,nq,nhinge,nsplit,H,dim_null,dim_null_patch,dimK,rankD,dimK_min4_2k,",
          "dimK_2rankD,c1_fit_relerr,c1_K_vs_basis,c1_affine_sec,cons,p0err,",
          "c2_conf_res,c2_maxdist,c2_tnorm,c2_inv,",
          "c3_hit_rand,c3_hit_aniso,c3_cert,c3_eps_cert,c3_theta_exact,c3_nu_maxerr,",
          "c4_q,c4_r,c4_thetac,c4_thetac_meas,c4_abs_err,c4_area_fit_err,c4_geom_err,",
          "reached,nsuperfaces,max_mismatch,err\n")
    nrun = 0
    for (pi_, P) in enumerate(pats)
        pi0 = pi_ - 1
        pi0 % nshard != shard && continue
        nrun >= limit && break
        nrun += 1
        if !P.ok
            print(csv, P.name, K7_MAIN_EMPTY, P.err, "\n")
            continue
        end
        S, err = prepare(P)
        if S === nothing
            print(csv, P.name, K7_MAIN_EMPTY, err, "\n")
            continue
        end
        tm = Timer()

        # the repo's finite-patch periodic system on the same cell, for comparison
        dnull_patch = -1
        let cm = K.Mesh(copy(P.cell.X), deepcopy(P.cell.faces))
            cm.sigma = copy(P.cell.sigma)
            K.add_periodic_pairs!(cm, Vec2(P.T[1, 1], P.T[2, 1]), Vec2(P.T[1, 2], P.T[2, 2]))
            K.build_topology!(cm)
            cc = K.make_cut(cm)
            hh = K.holes_partition(cc)
            sys = K.assemble_system(cc, hh, cm.X, K.Periodic)
            dnull_patch = K.solve_system(sys, cm.X).dim_null
        end

        AS = achievable(S)
        X0pts = K.matrix_to_points(S.X0)
        J0 = K_at(S, X0pts)
        reached = 0
        mism = 0.0
        let
            K.set_super_positions!(S.sp, S.q, X0pts)
            D = K.deploy(S.cut, S.sp.mesh.X, 0.3)
            reached = count(D.reached)
            mism = D.max_mismatch
        end

        # ---- C1: fit K from two angles of FORWARD KINEMATICS, predict 50 angles
        fit_relerr = -1.0; K_vs_basis = -1.0
        let
            K.set_super_positions!(S.sp, S.q, X0pts)
            th1 = 0.4; th2 = 1.1
            P1, _ = P_fk(S, th1)
            P2, _ = P_fk(S, th2)
            Msys = Mat2(cos(th1 / 2), cos(th2 / 2), sin(th1 / 2), sin(th2 / 2))  # column-major
            Mi = inv(Msys)
            P0f = zeros(2, 2); Qf = zeros(2, 2)
            for cix in 1:2
                a = P1[:, cix]; b = P2[:, cix]
                P0f[:, cix] = Mi[1, 1] * a + Mi[1, 2] * b
                Qf[:, cix] = Mi[2, 1] * a + Mi[2, 2] * b
            end
            Kfit = Mat2(Qf * inv(P0f))
            K_vs_basis = maximum(abs, Kfit - J0.K)
            fit_relerr = 0.0
            for i in 1:50
                th = pi * i / 50.0
                Pm, _ = P_fk(S, th)
                Pp = cos(th / 2) * P0f + sin(th / 2) * Qf
                fit_relerr = max(fit_relerr, norm(Pp - Pm) / max(1e-12, norm(Pm)))
            end
        end

        # ---- C1b: K affine in t (second differences)
        affine_sec = 0.0
        if S.k > 0
            rng = K.MT19937(UInt32(4242) + UInt32(977) * UInt32(pi0))
            G = K.NormalDist(0.0, 1.0)
            for _ in 1:31
                t = zeros(2 * S.k); d = zeros(2 * S.k)
                for i in 1:(2 * S.k)   # interleaved draws (the draw order is part of the population's definition)
                    t[i] = 0.3 * S.med_edge * K.normal(G, rng)
                    d[i] = 0.3 * S.med_edge * K.normal(G, rng)
                end
                A0 = K_at(S, shape_point(S, t)).K
                A1 = K_at(S, shape_point(S, t + d)).K
                A2 = K_at(S, shape_point(S, t + 2 * d)).K
                affine_sec = max(affine_sec, maximum(abs, A0 - 2 * A1 + A2))
            end
        end

        # ---- C2: conformal design
        conf_res = -1.0; maxdist = -1.0; tnorm = -1.0
        c2_inv = -1
        if AS.dimK >= 1
            C, rhs, t = conformal_t(AS, S.k)
            conf_res = norm(C * t - rhs) / max(1e-12, norm(rhs) + 1.0)
            tnorm = norm(t)
            Xc = shape_point(S, t)
            Kc = K_at(S, Xc).K
            maxdist = 0.0
            for i in 0:200
                th = pi * i / 200.0
                maxdist = max(maxdist, K.conformal_distortion(K.jac_at(Kc, th)))
            end
            _, c2_inv = cell_area_sum(S, Xc)
        end

        # ---- C3: designable Poisson family
        hit_rand = -1.0; hit_aniso = -1.0; eps_cert = 0.0; theta_exact = 0.0; nu_maxerr = -1.0
        c3_cert = 0
        if AS.dimK >= 1
            rng = K.MT19937(UInt32(999) + UInt32(613) * UInt32(pi0))
            G = K.NormalDist(0.0, 1.0)
            Kstar = K_of_t(AS, gaussian_vec(2 * S.k, 0.25 * S.med_edge, G, rng))
            tstar = cod_solve(AS.A, k_target_rhs(Kstar, AS.K0))
            hit_rand = maximum(abs, K_of_t(AS, tstar) - Kstar)
            Kan = Mat2(1.0, 0.0, 0.0, -0.5)
            hit_aniso = maximum(abs, K_of_t(AS, cod_solve(AS.A, k_target_rhs(Kan, AS.K0))) - Kan)

            # free subspace of the target-hitting affine set
            nfree = 2 * S.k - AS.dimK
            Vf = free_subspace(AS.A, nfree)
            function score(t)
                Xq = shape_point(S, t)
                K.set_super_positions!(S.sp, S.q, Xq)
                return certify(S.cut, S.sp.mesh.X, false).eps_max
            end
            t_best = tstar
            best = score(tstar)
            budget = K.n_faces(S.sp.mesh) > 600 ? 8 : 30
            for _ in 1:budget
                nfree > 0 || break
                t = tstar + Vf * gaussian_vec(nfree, 0.4 * S.med_edge, G, rng)
                v = score(t)
                if v > best
                    best = v
                    t_best = t
                end
            end
            # pattern search on the free coordinates
            if nfree > 0
                step = 0.3 * S.med_edge
                it = 0
                itmax = budget > 8 ? 25 : 8
                while it < itmax && step > 1e-3 * S.med_edge
                    it += 1
                    improved = false
                    for i in 1:nfree
                        improved && break
                        for sgn in (-1, 1)
                            improved && break
                            s_ = zeros(nfree)
                            s_[i] = sgn * step
                            t = t_best + Vf * s_
                            v = score(t)
                            if v > best + 1e-12
                                best = v
                                t_best = t
                                improved = true
                            end
                        end
                    end
                    improved || (step *= 0.5)
                end
            end
            eps_cert = best
            c3_cert = best > 0 ? 1 : 0
            Xb = shape_point(S, t_best)
            K.set_super_positions!(S.sp, S.q, Xb)
            theta_exact = K.n_faces(S.sp.mesh) <= 500 ? certify(S.cut, S.sp.mesh.X, true).theta_exact : -1.0
            Kb = K_at(S, Xb).K
            nu_maxerr = nu_max_error(S, Kb, J0.P0, Xb)
        end

        # ---- C4: second closed angle (periodic form)
        c4q = 0.0; c4r = 0.0; thc = -1.0; thc_meas = -1.0; area_fit = -1.0; geom_err = -1.0
        let
            S.cell_face_area, _ = cell_area_sum(S, X0pts)
            Kj = J0.K
            dP0 = abs(det(J0.P0))
            c4q = -dP0 * (det(Kj) - 1) / 2
            c4r = dP0 * tr(Kj) / 2
            thc = second_closed_angle(c4q, c4r)
            function Ameas(th)
                Pm, _ = P_fk(S, th)
                return abs(det(Pm)) - S.cell_face_area
            end
            area_fit = 0.0
            sc = abs(c4q) + abs(c4r) + abs(dP0)
            for i in 1:50
                th = 2 * pi * i / 51.0
                area_fit = max(area_fit, abs(Ameas(th) - (c4q * (cos(th) - 1) + c4r * sin(th))) / sc)
            end
            lo = 1e-5; hi = -1.0
            prev = Ameas(lo)
            for i in 1:4000
                th = 2 * pi * i / 4000.0
                v = Ameas(th)
                if (v < 0) != (prev < 0)
                    hi = th
                    break
                end
                lo = th; prev = v
            end
            if hi > 0
                for _ in 1:70
                    mid = 0.5 * (lo + hi)
                    (Ameas(mid) < 0) != (Ameas(lo) < 0) ? (hi = mid) : (lo = mid)
                end
                thc_meas = 0.5 * (lo + hi)
            end
            # geometric cross-check: trace the holes of the super patch and sum the ones
            # whose centroid lands in the central cell
            thg = 0.35 * min(thc > 0 ? thc : Float64(pi), Float64(pi))
            try
                D = K.deploy(S.cut, S.sp.mesh.X, thg)
                if !K.has_collision(S.cut, D.Y, 1e-12)
                    cyc = K.holes_geometric_cycles(S.cut, D.Y)
                    # assign a traced hole to the cell of its (min class, min offset) vertex --
                    # a translation-covariant rule, so each quotient hole lands in one cell.
                    acc = 0.0
                    ncell = 0
                    for cy in cyc
                        bc = 1 << 30
                        bo = K.Vec2i(0, 0)
                        for v in cy
                            sv = S.cut.prime_to_original[v]
                            cl = S.sp.vert_class[sv]
                            o = S.sp.vert_off[sv]
                            if cl < bc || (cl == bc && (o[1] < bo[1] || (o[1] == bo[1] && o[2] < bo[2])))
                                bc = cl
                                bo = o
                            end
                        end
                        if bo[1] == 0 && bo[2] == 0
                            acc += cycles_area(D.Y, [cy])
                            ncell += 1
                        end
                    end
                    pred = c4q * (cos(thg) - 1) + c4r * sin(thg)
                    geom_err = (ncell == S.q.H && abs(pred) > 1e-12) ? abs(acc - pred) / abs(pred) : -2.0
                end
            catch
            end
        end

        print(csv, P.name, ",", P.family, ",", K.n_faces(S.q), ",", S.q.nq, ",",
              S.q.n_hinge, ",", S.q.n_split, ",", S.q.H, ",", S.k, ",", dnull_patch,
              ",", AS.dimK, ",", AS.rankD, ",", min(4, 2 * S.k), ",",
              2 * AS.rankD, ",", g12(fit_relerr), ",", g12(K_vs_basis), ",", g12(affine_sec), ",",
              g12(J0.consistency), ",", g12(J0.p0_err), ",", g12(conf_res), ",", g12(maxdist), ",",
              g12(tnorm), ",", c2_inv, ",", g12(hit_rand), ",", g12(hit_aniso), ",", c3_cert,
              ",", g12(eps_cert), ",", g12(theta_exact), ",", g12(nu_maxerr), ",", g12(c4q), ",",
              g12(c4r), ",", g12(thc), ",", g12(thc_meas), ",",
              g12(thc_meas > 0 ? abs(thc - thc_meas) : -1), ",", g12(area_fit), ",", g12(geom_err),
              ",", reached, ",", K.n_faces(S.sp.mesh), ",", g12(mism), ",\n")
        flush(csv)
        @printf("[%2d] %-26s k=%2d dimK=%d rankD=%d c1=%.1e aff=%.1e conf=%.1e dist=%.1e cert=%d eps=%.4f thc=%.6f/%.6f  %.1fs\n",
                pi0, P.name, S.k, AS.dimK, AS.rankD, fit_relerr, affine_sec, conf_res,
                maxdist, c3_cert, eps_cert, thc, thc_meas, s(tm))
    end
    close(csv)
end

function main(args::Vector{String})
    stage = "run"
    outdir = joinpath(REPO, "results", "experiments", "k7")
    shard = 0; nshard = 1
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--stage" && i < length(args); stage = args[i+1]; i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--shard" && i < length(args); shard = arg_i(args[i+1]); i += 2
        elseif a == "--nshard" && i < length(args); nshard = arg_i(args[i+1]); i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    pats = population_k7()

    if stage == "pop"
        for P in pats
            @printf("%-28s ok=%d %s\n", P.name, Int(P.ok), P.err)
        end
        return 0
    end
    stage == "bounded" && (stage_bounded(outdir); return 0)
    stage == "perdbg" && (stage_perdbg(outdir, pats, shard, nshard, limit); return 0)
    stage == "c4sweep" && (stage_c4sweep(outdir, pats); return 0)
    stage == "nu" && (stage_nu(outdir, pats); return 0)
    stage == "c3" && (stage_c3(outdir, pats, shard, nshard, limit); return 0)
    stage_main(outdir, pats, shard, nshard, limit)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
