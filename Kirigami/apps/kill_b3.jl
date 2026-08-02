# kill_b3 -- B3 of ideas/round2_theorist_b.md: "tr K is the expansion budget".
# Port of code/apps/kill_b3.cpp.
#
# STAGE b1     the GLOBAL half of B1 that derivations/scratch/check_b1 could not
#              close: on a fixed-boundary patch the total void area (holes AND notches)
#              must equal the border functional's first harmonic, and the edge-wise
#              budget W + R must equal 2 B(X).  Bar 1e-10.
# STAGE retro  retrodiction.  For every K7 C3 row (reproduced BIT FOR BIT -- the state
#              machinery is kill_k7.jl's, included below exactly as the C++ copied it
#              verbatim from kill_k7.cpp, and the CSV is diffed against
#              results/kill/k7/k7_c3_all_v2.csv) record tr K of the realised design,
#              tau* = W / det P_0, the 0+ margin min_e q_e and the certificate, plus the
#              VARIATION of W over the achievable set, which is the hard-fail clause: if
#              W moves as much as det P_0 tr K does, tau* is not a threshold.
# STAGE con    the algorithm.  Re-solve the target least squares over the achievable
#              set under the affine constraint tr K >= tau* + delta |tau*|, spend the
#              leftover freedom exactly as C3 did, and re-certify.
#
# The budget algebra lives in method/budget.jl; only the drivers are here.
#
#   julia --project=Kirigami Kirigami/apps/kill_b3.jl [--stage run|b1|pdiag|summary]
#         [--out DIR] [--shard S] [--nshard M] [--limit K]
#
# Population: kill_k7.jl's periodic patterns (21 tilings + 12 torus Voronoi), rebuilt from
# the bit-exact generators. Outputs go to results/kill/b3_julia/ (the C++ wrote
# results/kill/b3/). `--limit K` stops after K patterns (or K graphs of stage b1).
include(joinpath(@__DIR__, "kill_k7.jl"))   # population, PState, achievable set, certificate, 0+ objective

# ===========================================================================
# B3 additions
# ===========================================================================

# The periodic budget of one cell at a shape-space point.
# THE PERIODIC TELESCOPING, and why a bare per-cell sum is not enough.  On a patch the
# forward sum of <x_a - x_b, u_f> over each closed face polygon vanishes and leaves the
# border; on a torus there is no border, and what is left instead is the MONODROMY of u.
# Summing over the 2|E_0| half-edge CLASSES, each at its representative whose face lies
# in the reference cell F_0, the side-1 representative of edge e sits at e + s(e) where
# s(e) = -(cell offset of f1), and periodic_jacobian.jl's
#      u_{f+t} = u_f + w_t + sigma_f t/2
# turns the vanishing face sum into
#      W + R = M ,   M := sum_{e in E_0} <x_b - x_a, w_{s(e)} + sigma_{f1} T s(e) / 2> ,
# with [w_h w_v] = 1/2 Jrot^T Q read off the same PeriodicJac.  M = det(P_0) tr K is the
# periodic form of (B.5) and is what `identity_err` tests.
Base.@kwdef mutable struct CellBudget
    W::Float64 = 0.0; R::Float64 = 0.0       # hinge and split sums over one cell
    bW::Float64 = 0.0; bR::Float64 = 0.0     # their (1 - cos) partners
    trK::Float64 = 0.0; detK::Float64 = 0.0; detP0::Float64 = 0.0
    identity_err::Float64 = 0.0   # |W + R - det(P0) tr K| / scale        -- (B.5) sin part
    b5_err::Float64 = 0.0         # |bW + bR - det(P0)(1 - det K)| / scale -- (B.5) cos part
    u_closure::Float64 = 0.0
    tau::Float64 = 0.0            # W / det P0
    n_split_cell::Int = 0; n_hinge_cell::Int = 0; n_preimage::Int = 0
    mask_ok::Bool = false
end

function cell_budget(S::PState, Xq::Vector{Vec2}, keep::Vector{Bool}, nh_expect::Int, ns_expect::Int)
    B = CellBudget()
    K.set_super_positions!(S.sp, S.q, Xq)
    m = S.sp.mesh
    X = m.X
    u, B.u_closure = K.face_potential(S.cut, X)

    # WHY NOT A BARE PER-EDGE SUM.  Both copies of an interior edge contribute
    #     <x_a - x_b, u_{f1} - u_{f0}>
    # to the a-coefficient of the preimage that carries them.  For a SPLIT edge
    # sigma_{f0} = sigma_{f1}, so u_{f1} - u_{f0} is unchanged by a lattice translation and
    # the term is a function of the quotient edge alone.  For a HINGE edge
    # sigma_{f0} = -sigma_{f1} and  u_{f1+t} - u_{f0+t} = u_{f1} - u_{f0} + sigma_{f1} t,
    # so the term depends on WHICH LIFT of the edge the preimage walk actually uses.
    # Summing one arbitrary representative per quotient edge therefore misses a defect that
    # is measured, not small (up to 1.9 det P_0 on the K7 population).  The budget of a cell
    # is instead the sum over the H PREIMAGES anchored in the reference cell, each taken at
    # the lift its own walk uses.
    bt = K.budget_terms(S.cut, X, u, keep)
    B.W = bt.W
    B.R = bt.R
    B.bW = bt.bW
    B.bR = bt.bR
    B.n_hinge_cell = count(e -> keep[e], S.cut.hinge_edges)
    B.n_split_cell = count(e -> keep[e], S.cut.split_edges)
    B.mask_ok = (B.n_hinge_cell == nh_expect && B.n_split_cell == ns_expect)

    Bs = K.deploy_basis(S.cut, X)
    pj = K.periodic_jacobian(S.sp, S.q, S.cut, Bs)
    B.trK = tr(pj.K)
    B.detK = det(pj.K)
    B.detP0 = det(pj.P0)
    scale = max(1.0, abs(B.W) + abs(B.R))
    B.identity_err = abs(B.W + B.R - B.detP0 * B.trK) / scale
    B.b5_err = abs(B.bW + B.bR - B.detP0 * (1 - B.detK)) / max(1.0, abs(B.bW) + abs(B.bR))
    B.tau = B.W / B.detP0
    return B
end

# Rank-based AUC with ties at 0.5.  label true = positive class.
function auc_b3(score::Vector{Float64}, label::Vector{Bool})
    np = count(label); nn = length(label) - np
    (np == 0 || nn == 0) && return -1.0
    s = 0.0
    for i in eachindex(score), j in eachindex(score)
        (!label[i] || label[j]) && continue
        s += score[i] > score[j] ? 1.0 : (score[i] == score[j] ? 0.5 : 0.0)
    end
    return s / (np * nn)
end

# ---------------------------------------------------------------- stage b1

function stage_b1(outdir::String, limit::Int)
    csv = open(joinpath(outdir, "b3_b1_global.csv"), "w")
    print(csv, "kind,par,nfaces,nsplit,nborder,n_preimage,n_notch,B_a,B_b,W,R,",
          "identity_rel,void_rel_worst,hole_only_a,notch_a\n")
    rng = K.MT19937(UInt32(20260904))
    cases = [("squares", 3.5), ("triangles", 3.5), ("hexagons", 3.5), ("kagome", 3.5),
             ("snub_square", 3.5), ("truncated_square", 3.5), ("t3_4_3_12", 4.0),
             ("squares", 5.0), ("hexagons", 5.0), ("kagome", 4.5)]
    th = [0.05, 0.2, 0.4, 0.7, 1.0, 1.3]
    worst_id = 0.0; worst_void = 0.0
    ngraph = 0; n_with_notch = 0
    for (kind, par) in cases
        ngraph >= limit && break
        m = K.generate(kind, [par], rng)
        K.build_topology!(m)
        m.sigma = K.assign_orientation_relaxation(m, rng).sigma
        K.build_topology!(m)
        c = K.make_cut(m)
        hs = K.holes_partition(c)
        sys = K.assemble_system(c, hs, m.X, K.Fixed)
        sr = K.solve_system(sys, m.X)
        sr.projection_ok || continue
        X = K.matrix_to_points(sr.X0)
        K.deployable(K.hole_residuals(c, X, hs), 1e-7) || continue
        ngraph += 1
        u, clo = K.face_potential(c, X)
        bt = K.budget_terms(c, X, u)
        B = K.border_functional(c, X, u)
        sa = max(1.0, abs(bt.W) + abs(bt.R))
        id_rel = abs(bt.W + bt.R - 2 * B.a) / sa
        vrel = 0.0
        for t in th
            meas = K.measured_void_area(c, K.deploy(c, X, t).Y)
            pred = B.a * K.libm_sin(t) - B.b * (1 - K.libm_cos(t))
            vrel = max(vrel, abs(meas - pred) / max(1.0, abs(B.a) + abs(B.b)))
        end
        # holes only, via the geometric tracer -- the half check_b1 already closed
        hole_a = 0.0
        pv_faces = [Int[] for _ in 1:c.n_prime_vertices]
        for f in 1:K.n_faces(m), pv in c.prime_faces[f]
            push!(pv_faces[pv], f)
        end
        cycles = K.holes_geometric_cycles(c, K.deploy(c, X, th[1]).Y)
        for cy in cycles
            length(cy) < 3 && continue
            ok = true
            a = 0.0
            for i in eachindex(cy)
                pa = cy[i]; pb = cy[mod1(i + 1, length(cy))]
                f = 0
                for fa in pv_faces[pa], fb in pv_faces[pb]
                    fa == fb && (f = fa)
                end
                if f == 0
                    ok = false
                    break
                end
                a += 0.5 * K._dot2(X[c.prime_to_original[pa]] - X[c.prime_to_original[pb]], u[f])
            end
            ok && (hole_a += a)
        end
        n_notch = length(hs.all) - K.n_interior_holes(hs)
        n_notch > 0 && (n_with_notch += 1)
        worst_id = max(worst_id, id_rel)
        worst_void = max(worst_void, vrel)
        print(csv, kind, ",", g12(par), ",", K.n_faces(m), ",", K.n_split(c), ",",
              length(c.border_edges), ",", length(hs.all), ",", n_notch, ",", g12(B.a), ",",
              g12(B.b), ",", g12(bt.W), ",", g12(bt.R), ",", g12(id_rel), ",", g12(vrel), ",",
              g12(hole_a), ",", g12(B.a - hole_a), "\n")
        @printf("%-18s par=%.1f  W+R vs 2B rel=%.2e  void rel=%.2e  notches=%d  holes_a=%.4f total_a=%.4f\n",
                kind, par, id_rel, vrel, n_notch, hole_a, B.a)
    end
    close(csv)
    @printf("\nB1 GLOBAL: %d graphs (%d with notches), worst |W+R-2B|/scale = %.3e, worst void-area rel = %.3e  -> %s (bar 1e-10)\n",
            ngraph, n_with_notch, worst_id, worst_void,
            (worst_id < 1e-10 && worst_void < 1e-10) ? "PASS" : "FAIL")
    return 0
end

# The C3 free-subspace pattern search, verbatim from kill_k7's c3 stage, factored
# so the baseline and the constrained re-solve spend the leftover freedom identically.
function free_search(S::PState, tstar::Vector{Float64}, Vf::Matrix{Float64}, nfree::Int, rng::K.MT19937)
    G = K.NormalDist(0.0, 1.0)
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
    return tof(s_best)
end

mutable struct Realised
    X::Vector{Vec2}
    minq::Float64
    mina::Float64
    ct::Cert
end

function realise(S::PState, t::Vector{Float64})
    X = shape_point(S, t)
    K.set_super_positions!(S.sp, S.q, X)
    minq = 1e300; mina = 1e300
    for v in K.zero_plus_q(S.cut, S.sp.mesh.X)
        minq = min(minq, v)
    end
    for fi in 1:K.n_faces(S.sp.mesh)
        mina = min(mina, K.face_signed_area(S.sp.mesh, fi))
    end
    ct = certify(S.cut, S.sp.mesh.X, K.n_faces(S.sp.mesh) <= 500)
    return Realised(X, minq, mina, ct)
end

# Diagnostic: the per-cell void-area harmonic measured by forward kinematics against
# (a) the periodic closed form dP0 tr K / 2 that K7 C4 already verified and (b) the
# edge-wise budget (W + R) / 2.
function stage_pdiag(pats, nmax::Int)
    th = [0.05, 0.2, 0.4, 0.7, 1.0, 1.3]
    done = 0
    for P in pats
        done >= nmax && break
        P.ok || continue
        S, _ = prepare(P)
        S === nothing && continue
        keep, npre = K.periodic_cell_edges(S.cut, S.sp)
        Xq = K.matrix_to_points(S.X0)
        B = cell_budget(S, Xq, keep, S.q.n_hinge, S.q.n_split)
        S.cell_face_area, _ = cell_area_sum(S, Xq)
        Mfit = zeros(6, 2)
        rhs = zeros(6)
        for i in 1:6
            Pm, _ = P_fk(S, th[i])
            rhs[i] = abs(det(Pm)) - S.cell_face_area
            Mfit[i, 1] = K.libm_sin(th[i])
            Mfit[i, 2] = -(1.0 - K.libm_cos(th[i]))
        end
        fit = Mfit \ rhs
        # origin dependence: shift every quotient position by a non-lattice vector
        cshift = Vec2(0.37 * S.med_edge, -0.61 * S.med_edge)
        Xs = [v + cshift for v in Xq]
        Bs2 = cell_budget(S, Xs, keep, S.q.n_hinge, S.q.n_split)
        K.set_super_positions!(S.sp, S.q, Xq)
        q0 = 1e300
        for v in K.zero_plus_q(S.cut, S.sp.mesh.X)
            q0 = min(q0, v)
        end
        K.set_super_positions!(S.sp, S.q, Xs)
        q1 = 1e300
        for v in K.zero_plus_q(S.cut, S.sp.mesh.X)
            q1 = min(q1, v)
        end
        K.set_super_positions!(S.sp, S.q, Xq)
        @printf("   shift: dW=%10.4f dR=%10.4f dtrK=%10.4f ddetP0=%9.2e dminq=%10.4f\n",
                Bs2.W - B.W, Bs2.R - B.R, Bs2.trK - B.trK, Bs2.detP0 - B.detP0, q1 - q0)
        @printf("%-24s H=%2d npre=%2d hinge %d/%d(%d) split %d/%d(%d) | 2a=%10.4f dP0trK=%10.4f W+R=%10.4f (W=%9.4f R=%9.4f) ratio=%7.4f\n",
                P.name, S.q.H, npre, B.n_hinge_cell, S.q.n_hinge, npre,
                B.n_split_cell, S.q.n_split, npre, 2 * fit[1], B.detP0 * B.trK, B.W + B.R, B.W,
                B.R, (B.detP0 * B.trK) / (B.W + B.R))
        done += 1
    end
    return 0
end

function read_csv_b3(path::String)
    isfile(path) || return String[], Vector{String}[]
    lines = readlines(path)
    isempty(lines) && return String[], Vector{String}[]
    hdr = String.(split(lines[1], ','))
    rows = [String.(split(l, ',')) for l in lines[2:end] if !isempty(l)]
    return hdr, rows
end
colof(h::Vector{String}, name::String) = something(findfirst(==(name), h), -1)

function stage_summary(outdir::String)
    H, R = read_csv_b3(joinpath(outdir, "b3_retro.csv"))
    if isempty(R)
        @printf("no b3_retro.csv in %s\n", outdir)
        return 1
    end
    cGap = colof(H, "gap"); cQ = colof(H, "min_q"); cCert = colof(H, "cert")
    cId = colof(H, "identity_rel"); cNs = colof(H, "nsplit_cell"); cdP = colof(H, "detP0")
    cR = colof(H, "R"); cWr = colof(H, "W_ratio"); cName = colof(H, "name")
    gap = Float64[]; rsplit = Float64[]; gapc = Float64[]; rsplitc = Float64[]
    lq = Bool[]; lc = Bool[]; lqc = Bool[]; lcc = Bool[]
    worst_id = 0.0; worst_ratio = 0.0; med_ratio = 0.0
    nq = 0; nc = 0; viol = 0
    ratios = Float64[]
    seen = Set{String}()
    for r in R
        g = parse(Float64, r[cGap]); q = parse(Float64, r[cQ])
        ce = parse(Int, r[cCert])
        parse(Int, r[cNs]) == 0 && continue   # no split cuts: min q undefined
        push!(gap, g)
        push!(rsplit, parse(Float64, r[cR]) / parse(Float64, r[cdP]) / parse(Float64, r[cNs]))
        push!(lq, q > 0)
        push!(lc, ce == 1)
        nq += (q > 0)
        nc += ce
        (q > 0 && g < 0) && (viol += 1)
        idv = parse(Float64, r[cId])
        worst_id = max(worst_id, idv)
        if idv < 1e-9   # rows on which (B.5) is exact
            push!(gapc, g)
            push!(rsplitc, rsplit[end])
            push!(lqc, q > 0)
            push!(lcc, ce == 1)
        end
        if !(r[cName] in seen)
            push!(seen, r[cName])
            wr = parse(Float64, r[cWr])
            if wr >= 0
                push!(ratios, wr)
                worst_ratio = max(worst_ratio, wr)
            end
        end
    end
    sort!(ratios)
    isempty(ratios) || (med_ratio = ratios[length(ratios) ÷ 2 + 1])
    @printf("=== B3 stage 1 (retrodiction) ===\n")
    @printf("rows with split cuts        %d of %d\n", length(gap), length(R))
    @printf("worst |W+R-detP0 trK|/scale %.3e\n", worst_id)
    @printf("min q > 0                   %d      certified %d\n", nq, nc)
    @printf("AUC(gap = trK - tau*)  for min q > 0   %.4f\n", auc_b3(gap, lq))
    @printf("AUC(gap = trK - tau*)  for certified   %.4f\n", auc_b3(gap, lc))
    @printf("AUC(R / (detP0 n_split)) for min q > 0 %.4f\n", auc_b3(rsplit, lq))
    @printf("AUC(R / (detP0 n_split)) for certified %.4f\n", auc_b3(rsplit, lc))
    @printf("rows with min q > 0 at tr K < tau*     %d   (must be 0)\n", viol)
    @printf("--- restricted to the %d rows where (B.5) is exact to 1e-9 ---\n", length(gapc))
    @printf("AUC(gap) min q > 0 %.4f   certified %.4f | AUC(R/(detP0 n_split)) %.4f / %.4f\n",
            auc_b3(gapc, lqc), auc_b3(gapc, lcc), auc_b3(rsplitc, lqc), auc_b3(rsplitc, lcc))
    @printf("W spread / budget spread over K: median %.4f, max %.4f (%d patterns)\n",
            med_ratio, worst_ratio, length(ratios))

    H2, C = read_csv_b3(joinpath(outdir, "b3_constrained.csv"))
    isempty(C) && return 0
    dD = colof(H2, "delta"); dBC = colof(H2, "base_cert"); dCC = colof(H2, "con_cert")
    dBQ = colof(H2, "base_minq"); dCQ = colof(H2, "con_minq"); dN = colof(H2, "name")
    dT = colof(H2, "target"); dK = colof(H2, "K_shift"); dOK = colof(H2, "con_ok")
    @printf("\n=== B3 stage 2 (constrained re-solve) ===\n")
    for dl in (0.0, 0.1, 0.5)
        base = 0; con = 0; up = 0; down = 0; n = 0; qb = 0; qc = 0
        kshift = 0.0
        flips = ""
        for r in C
            (abs(parse(Float64, r[dD]) - dl) > 1e-9 || parse(Int, r[dOK]) == 0) && continue
            n += 1
            b = parse(Int, r[dBC]); c2 = parse(Int, r[dCC])
            base += b
            con += c2
            qb += parse(Float64, r[dBQ]) > 0
            qc += parse(Float64, r[dCQ]) > 0
            kshift = max(kshift, parse(Float64, r[dK]))
            c2 > b && (up += 1; flips *= " +" * r[dN] * "/" * r[dT])
            c2 < b && (down += 1; flips *= " -" * r[dN] * "/" * r[dT])
        end
        @printf("delta=%.1f  n=%3d  certified %2d -> %2d  (up %d, down %d)   min q > 0: %2d -> %2d   max |K - K*| = %.3f\n",
                dl, n, base, con, up, down, qb, qc, kshift)
        isempty(flips) || @printf("           flips:%s\n", flips)
    end
    return 0
end

function main(args::Vector{String})
    stage = "run"
    outdir = joinpath(REPO, "results", "kill", "b3_julia")
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
    stage == "b1" && return stage_b1(outdir, limit)
    stage == "pdiag" && return stage_pdiag(population_k7(), min(33, limit))
    stage == "summary" && return stage_summary(outdir)

    pats = population_k7()
    rcsv = open(joinpath(outdir, "b3_retro_" * string(shard) * ".csv"), "w")
    ccsv = open(joinpath(outdir, "b3_constrained_" * string(shard) * ".csv"), "w")
    print(rcsv, "name,family,dimK,nfree,nsplit_cell,nhinge_cell,mask_ok,u_closure,identity_rel,",
          "b5_rel,npre,target,hit_err,detP0,trK,W,R,tau,gap,min_q,min_area,cert,eps_cert,zp_feasible,",
          "W_lo,W_hi,W_spread,budget_spread,trK_lo,trK_hi,W_ratio\n")
    print(ccsv, "name,family,target,delta,base_trK,base_tau,base_gap,base_cert,base_minq,",
          "con_ok,con_trK,con_tau,con_gap,con_cert,con_minq,con_eps,K_shift,tau_fp_err\n")

    deltas = (0.0, 0.1, 0.5)
    nrun = 0
    for (pi_, P) in enumerate(pats)
        pi0 = pi_ - 1   # the C++ 0-based population index (seeds)
        pi0 % nshard != shard && continue
        P.ok || continue
        nrun >= limit && break
        S, _ = prepare(P)
        S === nothing && continue
        AS = achievable(S)
        AS.dimK < 1 && continue
        nrun += 1

        nh_exp = S.q.n_hinge; ns_exp = S.q.n_split
        keep, npre = K.periodic_cell_edges(S.cut, S.sp)

        rng = K.MT19937(UInt32(31337) + UInt32(613) * UInt32(pi0))
        G = K.NormalDist(0.0, 1.0)   # ONE distribution object for rng and r2, as the C++
        nfree = 2 * S.k - AS.dimK
        Vf = free_subspace(AS.A, nfree)
        vecI = [1.0, 0.0, 0.0, 1.0]
        pu = cod_solve(AS.A, vecI)
        den = dot(vecI, AS.A * pu)

        # W over the achievable set: 20 Gaussian points at the search scale.
        Wlo = 1e300; Whi = -1e300; Blo = 1e300; Bhi = -1e300; Tlo = 1e300; Thi = -1e300
        let r2 = K.MT19937(UInt32(777) + UInt32(9176) * UInt32(pi0))
            for _ in 1:20
                t = gaussian_vec(2 * S.k, 0.3 * S.med_edge, G, r2)
                B = cell_budget(S, shape_point(S, t), keep, nh_exp, ns_exp)
                Wlo = min(Wlo, B.W); Whi = max(Whi, B.W)
                Tlo = min(Tlo, B.trK); Thi = max(Thi, B.trK)
                bud = B.detP0 * B.trK
                Blo = min(Blo, bud); Bhi = max(Bhi, bud)
            end
        end
        Wspread = Whi - Wlo; Bspread = Bhi - Blo
        Wratio = Bspread > 0 ? Wspread / Bspread : -1.0

        for tgt in 0:1
            Kstar = if tgt == 0
                K_of_t(AS, gaussian_vec(2 * S.k, 0.25 * S.med_edge, G, rng))
            else
                Mat2(1.0, 0.0, 0.0, -0.5)
            end
            tstar = cod_solve(AS.A, k_target_rhs(Kstar, AS.K0))
            hit = maximum(abs, K_of_t(AS, tstar) - Kstar)
            trKstar = tr(K_of_t(AS, tstar))

            # ---- baseline: exactly the C3 design
            tb = free_search(S, tstar, Vf, nfree, rng)
            RB = realise(S, tb)
            BB = cell_budget(S, RB.X, keep, nh_exp, ns_exp)
            tname = tgt == 0 ? "random_in_K" : "diag_1_-0.5"
            print(rcsv, P.name, ",", P.family, ",", AS.dimK, ",", nfree, ",", ns_exp, ",",
                  nh_exp, ",", (BB.mask_ok ? 1 : 0), ",", g12(BB.u_closure), ",",
                  g12(BB.identity_err), ",", g12(BB.b5_err), ",", npre, ",", tname, ",",
                  g12(hit), ",", g12(BB.detP0), ",", g12(BB.trK), ",", g12(BB.W), ",", g12(BB.R), ",",
                  g12(BB.tau), ",", g12(BB.trK - BB.tau), ",", g12(RB.minq), ",",
                  g12(RB.mina), ",", (RB.ct.eps_max > 0 ? 1 : 0), ",", g12(RB.ct.eps_max), ",",
                  ((RB.minq > 0 && RB.mina > 0) ? 1 : 0), ",", g12(Wlo), ",", g12(Whi), ",",
                  g12(Wspread), ",", g12(Bspread), ",", g12(Tlo), ",", g12(Thi), ",", g12(Wratio), "\n")
            flush(rcsv)
            @printf("[%2d] %-26s %-12s trK=%9.4f tau=%9.4f gap=%9.4f minq=%9.2e cert=%d id=%.1e b5=%.1e npre=%.0f\n",
                    pi0, P.name, tname, BB.trK, BB.tau, BB.trK - BB.tau, RB.minq,
                    RB.ct.eps_max > 0 ? 1 : 0, BB.identity_err, BB.b5_err, Float64(npre))

            # ---- constrained re-solve, tr K >= tau* + delta |tau*|
            for dl in deltas
                tau = BB.tau; fp_err = 0.0
                tc = tstar
                ok = den > 1e-12 * max(1.0, norm(AS.A))
                if ok
                    for _ in 1:4
                        thresh = tau + dl * abs(tau)
                        cneed = thresh - trKstar
                        tc = cneed > 0 ? tstar + (cneed / den) * pu : tstar
                        Bn = cell_budget(S, shape_point(S, tc), keep, nh_exp, ns_exp)
                        fp_err = abs(Bn.tau - tau) / max(1.0, abs(tau))
                        tau = Bn.tau
                        fp_err < 1e-9 && break
                    end
                end
                r3 = K.MT19937((UInt32(4242) + UInt32(613) * UInt32(pi0) + UInt32(17) * UInt32(tgt) +
                                UInt32(trunc(Int, 1000 * dl))) % UInt32)
                tcf = ok ? free_search(S, tc, Vf, nfree, r3) : tb
                RC = realise(S, tcf)
                BC = cell_budget(S, RC.X, keep, nh_exp, ns_exp)
                kshift = norm(K_of_t(AS, tc) - Kstar)
                print(ccsv, P.name, ",", P.family, ",", tname, ",", g12(dl), ",", g12(BB.trK), ",",
                      g12(BB.tau), ",", g12(BB.trK - BB.tau), ",", (RB.ct.eps_max > 0 ? 1 : 0), ",",
                      g12(RB.minq), ",", (ok ? 1 : 0), ",", g12(BC.trK), ",", g12(BC.tau), ",",
                      g12(BC.trK - BC.tau), ",", (RC.ct.eps_max > 0 ? 1 : 0), ",", g12(RC.minq), ",",
                      g12(RC.ct.eps_max), ",", g12(kshift), ",", g12(fp_err), "\n")
                flush(ccsv)
                @printf("        d=%.1f  trK=%9.4f tau=%9.4f gap=%9.4f minq=%9.2e cert=%d |dK|=%.3f\n",
                        dl, BC.trK, BC.tau, BC.trK - BC.tau, RC.minq, RC.ct.eps_max > 0 ? 1 : 0, kshift)
            end
        end
    end
    close(rcsv)
    close(ccsv)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
