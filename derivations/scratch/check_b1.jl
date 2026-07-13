# check_b1.jl -- numeric check of the budget identity (B.1) of ideas/round2_theorist_b.md.
# Port of check_b1.cpp (same corpus, same seed, same printed lines).
#
# Run (from repo root):
#   julia --project=Kirigami derivations/scratch/check_b1.jl
#
# B1  per hole C:  A_C(theta) = a_C sin(theta) - b_C (1 - cos(theta))  with
#       a_C = 1/2 sum_{(a->b) in dC, face f} <x_a - x_b, u_f>
#       b_C = 1/2 sum_{(a->b) in dC, face f} sigma_f det(x_a - x_b, u_f)
#     measured: shoelace of the hole cycle at 6 angles, least-squares fit of (a, b).
# B1s split-edge term: the contribution of a split edge to a_C must be 1/2 <d_e, du_e>.
# B10 the area maximum sits at exactly half the second root: analytic, checked by
#     locating the max numerically.

using Kirigami
using LinearAlgebra
using Printf
using StaticArrays
const K = Kirigami
const Vec2 = K.Vec2

const J = SMatrix{2,2,Float64}(0, 1, -1, 0)   # [0 -1; 1 0]
det2(a::Vec2, b::Vec2) = a[1] * b[2] - a[2] * b[1]

function potential_u(c::K.CutStructure, X::Vector{Vec2})
    m = c.mesh
    F = K.n_faces(m)
    adj = [Tuple{Int,Int}[] for _ in 1:F]
    for e in c.hinge_edges
        ed = m.edges[e]
        f1 = m.half_edges[ed.he[1]].face; f2 = m.half_edges[ed.he[2]].face
        push!(adj[f1], (e, f2))
        push!(adj[f2], (e, f1))
    end
    u = fill(Vec2(0, 0), F)
    seen = falses(F)
    for s in 1:F
        seen[s] && continue
        seen[s] = true; u[s] = Vec2(0, 0)
        q = [s]; qi = 1
        while qi <= length(q)
            f = q[qi]; qi += 1
            for (e, g) in adj[f]
                target = u[f] + Float64(m.sigma[g]) * X[c.hinge_dir[e].src]
                if !seen[g]
                    seen[g] = true; u[g] = target; push!(q, g)
                end
            end
        end
    end
    w = 0.0
    for e in c.hinge_edges
        ed = m.edges[e]
        f1 = m.half_edges[ed.he[1]].face; f2 = m.half_edges[ed.he[2]].face
        xv = X[c.hinge_dir[e].src]
        w = max(w, norm(u[f2] - u[f1] - Float64(m.sigma[f2]) * xv))
    end
    return u, w
end

function shoelace(P::Vector{Vec2})
    a = 0.0
    n = length(P)
    for i in 1:n
        p = P[i]; q = P[i % n + 1]
        a += p[1] * q[2] - p[2] * q[1]
    end
    return 0.5 * a
end

function main()
    rng = K.MT19937(20260904)
    cases = [("squares", [3.5]), ("triangles", [3.5]), ("hexagons", [3.5]),
             ("kagome", [3.5]), ("snub_square", [3.5]),
             ("truncated_square", [3.5]), ("t3_4_3_12", [4.0]),
             ("squares", [5.0]), ("hexagons", [5.0]), ("kagome", [4.5])]
    th = [0.05, 0.2, 0.4, 0.7, 1.0, 1.3]

    worst_rel = 0.0; worst_abs = 0.0; worst_ratio = 0.0; worst_split = 0.0; worst_clo = 0.0
    where_rel = ""; where_ratio = ""
    n_holes = 0; n_graphs = 0; n_split_terms = 0

    for (kind, par) in cases
        m = K.generate(kind, par, rng)
        K.build_topology!(m)
        orep = K.assign_orientation_relaxation(m, rng)
        m.sigma = orep.sigma
        K.build_topology!(m)
        c = K.make_cut(m)
        hs = K.holes_partition(c)
        sys = K.assemble_system(c, hs, m.X, K.Fixed)
        sr = K.solve_system(sys, m.X)
        sr.projection_ok || continue
        X = K.matrix_to_points(sr.X0)
        res = K.hole_residuals(c, X, hs)
        K.deployable(res, 1e-7) || continue
        u, clo = potential_u(c, X)
        worst_clo = max(worst_clo, clo)

        # pv -> faces containing it
        pv_faces = [Int[] for _ in 1:c.n_prime_vertices]
        for f in 1:K.n_faces(m), pv in c.prime_faces[f]
            push!(pv_faces[pv], f)
        end

        d0 = K.deploy(c, X, th[1])
        cycles = K.holes_geometric_cycles(c, d0.Y)
        if isempty(cycles)
            @printf("%-18s no cycles\n", kind)
            continue
        end
        n_graphs += 1

        # deployed positions at every angle
        Y = [K.deploy(c, X, t).Y for t in th]

        gw = 0.0
        for cyc in cycles
            length(cyc) < 3 && continue
            # closed form (B.1) along the SAME directed walk
            a_cf = 0.0; b_cf = 0.0; ok = true
            for i in eachindex(cyc)
                pa = cyc[i]; pb = cyc[i % length(cyc) + 1]
                f = 0
                for fa in pv_faces[pa], fb in pv_faces[pb]
                    fa == fb && (f = fa)
                end
                if f < 1
                    ok = false
                    break
                end
                dv = X[c.prime_to_original[pa]] - X[c.prime_to_original[pb]]
                a_cf += 0.5 * dot(dv, u[f])
                b_cf += 0.5 * Float64(m.sigma[f]) * det2(dv, u[f])
            end
            ok || continue
            # measured areas + least squares fit of (a, b)
            M = Matrix{Float64}(undef, length(th), 2); rhs = Vector{Float64}(undef, length(th))
            for k in eachindex(th)
                P = [Y[k][pv] for pv in cyc]
                rhs[k] = shoelace(P)
                M[k, 1] = sin(th[k])
                M[k, 2] = -(1.0 - cos(th[k]))
            end
            fit = qr(M, ColumnNorm()) \ rhs
            resid = maximum(abs.(M * fit - rhs))
            scale = max(1e-12, abs(fit[1]) + abs(fit[2]))
            rel = (abs(fit[1] - a_cf) + abs(fit[2] - b_cf)) / scale
            n_holes += 1
            if rel > worst_rel
                worst_rel = rel; where_rel = kind
            end
            worst_abs = max(worst_abs, resid / scale)
            gw = max(gw, rel)
            # B10: numeric argmax vs half the analytic second root
            if abs(b_cf) > 1e-9 * scale
                tc = 2.0 * atan(a_cf, b_cf)
                tmax = atan(a_cf, b_cf)
                best = -1e300; targ = 0.0
                for s in 1:1999
                    t = 3.14159265358979 * s / 2000.0
                    A = a_cf * sin(t) - b_cf * (1 - cos(t))
                    if A > best
                        best = A; targ = t
                    end
                end
                if tc > 0.05 && tc < 3.1
                    r = abs(targ - tmax)
                    if r > worst_ratio
                        worst_ratio = r; where_ratio = kind
                    end
                end
            end
        end
        # (B.2s): split-edge contribution to its hole = 1/2 <d_e, du_e>, compared with the
        # same quantity read off the deployed offset  y_(a,g) - y_(a,f) = 2 sin(th/2) J du.
        for e in c.split_edges
            ed = m.edges[e]
            f = m.half_edges[ed.he[1]].face; g = m.half_edges[ed.he[2]].face
            va = ed.key.a
            du = u[g] - u[f]
            s0 = sin(th[3] * 0.5)
            off = Y[3][K.prime_vertex(c, g, va)] - Y[3][K.prime_vertex(c, f, va)]
            du_meas = (transpose(J) * off) / (2.0 * s0)
            worst_split = max(worst_split, norm(du_meas - du) / max(1.0, norm(du)))
            n_split_terms += 1
        end
        @printf("%-18s holes=%3d  worst rel (B.1) = %.3e\n", kind, length(cycles), gw)
    end

    @printf("\n--- summary over %d graphs, %d holes, %d split edges ---\n",
            n_graphs, n_holes, n_split_terms)
    @printf("u closure residual              %.3e\n", worst_clo)
    @printf("B.1  closed form vs fit (rel)   %.3e   (%s)\n", worst_rel, where_rel)
    @printf("B.1  first-harmonic fit resid   %.3e\n", worst_abs)
    @printf("B.2s du from deployed offset    %.3e\n", worst_split)
    @printf("B10  |argmax - theta_c/2|       %.3e   (%s)\n", worst_ratio, where_ratio)
end

main()
