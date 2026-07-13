# check_t1_t2.jl -- Deriver's own numeric checks for derivations/core.md, T1 and T2.
# Port of check_t1_t2.cpp (same corpus, same seed, same printed lines).
#
# Run (from repo root):
#   julia --project=Kirigami derivations/scratch/check_t1_t2.jl
#
# Checks, all against Kirigami/src/core (which the Deriver does not modify):
#   C1  the potential u on Gamma closes on EVERY hinge edge, not just a spanning tree
#       (path-independence of the path sum; T1 step 6)
#   C2  Y_theta = cos(th/2) C + sin(th/2) S with C_(v,f) = x_v and
#       S_(v,f) = J (2 u_f - sigma_f x_v), against deploy() -- fixes the sign convention
#   C3  hinge opening angle == theta at every hinge edge
#   C4  the two duplicates of a split edge are the SAME vector (parallel, T1 cor. ii)
#   C5  face signed areas are constant in theta (T3)
#   C6  A(Y_theta) sigma = 0 for every theta: the explicit face velocity field
#       (omega_f, w_f) satisfies the pin equation at every hinge, at every theta (T2)
#   C7  the vertex trajectory lies on the centred conic through C, S (T1 cor. i)

using Kirigami
using LinearAlgebra
using Printf
using StaticArrays
const K = Kirigami
const Vec2 = K.Vec2

const J = SMatrix{2,2,Float64}(0, 1, -1, 0)   # column-major: [0 -1; 1 0]

mutable struct Worst
    v::Float64
    where::String
end
Worst() = Worst(0.0, "")
function hit!(w::Worst, x::Float64, s::AbstractString)
    if x > w.v
        w.v = x; w.where = String(s)
    end
end

# BFS on Gamma computing the potential u_f with u_g - u_f = sigma_g * x_src(e).
# Returns u and the worst closure residual over NON-tree hinge edges (C1).
function potential_u(c::K.CutStructure, X::Vector{Vec2})
    m = c.mesh
    F = K.n_faces(m)
    adj = [Tuple{Int,Int}[] for _ in 1:F]   # face -> (edge, other face)
    for e in c.hinge_edges
        ed = m.edges[e]
        f1 = m.half_edges[ed.he[1]].face
        f2 = m.half_edges[ed.he[2]].face
        push!(adj[f1], (e, f2))
        push!(adj[f2], (e, f1))
    end
    u = fill(Vec2(0, 0), F)
    seen = falses(F)
    for s in 1:F
        seen[s] && continue
        seen[s] = true
        u[s] = Vec2(0, 0)
        q = [s]
        qi = 1
        while qi <= length(q)
            f = q[qi]; qi += 1
            for (e, g) in adj[f]
                xv = X[c.hinge_dir[e].src]
                target = u[f] + Float64(m.sigma[g]) * xv
                if !seen[g]
                    seen[g] = true
                    u[g] = target
                    push!(q, g)
                end
            end
        end
    end
    # closure over every hinge edge
    w = 0.0
    for e in c.hinge_edges
        ed = m.edges[e]
        f1 = m.half_edges[ed.he[1]].face
        f2 = m.half_edges[ed.he[2]].face
        xv = X[c.hinge_dir[e].src]
        w = max(w, norm(u[f2] - u[f1] - Float64(m.sigma[f2]) * xv))
        w = max(w, norm(u[f1] - u[f2] - Float64(m.sigma[f1]) * xv))
    end
    return u, w
end

function signed_area(P::Vector{Vec2})
    a = 0.0
    n = length(P)
    for i in 1:n
        p = P[i]; q = P[i % n + 1]
        a += p[1] * q[2] - p[2] * q[1]
    end
    return 0.5 * a
end

function main()
    rng = K.MT19937(20260903)
    cases = [("squares", [3.5]), ("triangles", [3.5]), ("hexagons", [3.5]),
             ("kagome", [3.5]), ("snub_square", [3.5]), ("truncated_square", [3.5]),
             ("t3_4_3_12", [4.0]), ("delaunay", [70.0, 10.0]), ("voronoi", [60.0, 10.0]),
             ("quad_random", [60.0, 10.0])]
    thetas = [0.0, 0.05, 0.3, 0.77, 1.2, 1.9, 2.6, 3.0]

    c1, c2, c3, c4, c5, c6, c7 = (Worst() for _ in 1:7)
    n_graphs = 0; n_deployable = 0

    for (kind, par) in cases
        for rep in 1:3
            m = K.generate(kind, par, rng)
            K.build_topology!(m)
            orep = K.assign_orientation_relaxation(m, rng)
            m.sigma = orep.sigma
            K.build_topology!(m)
            c = K.make_cut(m)
            hs = K.holes_partition(c)
            sys = K.assemble_system(c, hs, m.X, K.Fixed)
            sr = K.solve_system(sys, m.X)
            (!sr.projection_ok || sr.dim_null < 0) && continue
            X = K.matrix_to_points(sr.X0)
            res = K.hole_residuals(c, X, hs)
            n_graphs += 1
            K.deployable(res, 1e-7) || continue   # only test on genuinely deployable X
            n_deployable += 1

            u, closure = potential_u(c, X)
            hit!(c1, closure, kind)

            # C2/C3/C4/C5/C6/C7
            for th in thetas
                cc = cos(th * 0.5); ss = sin(th * 0.5)
                d = K.deploy(c, X, th)
                # closed form
                Y = fill(Vec2(0, 0), c.n_prime_vertices)
                for f in 1:K.n_faces(m)
                    for v in m.faces[f]
                        pv = K.prime_vertex(c, f, v)
                        C = X[v]
                        S = J * (2.0 * u[f] - Float64(m.sigma[f]) * X[v])
                        Y[pv] = cc * C + ss * S
                    end
                end
                # deploy() anchors face 1 at t=0 too, so the two agree exactly (no alignment).
                for i in 1:c.n_prime_vertices
                    hit!(c2, norm(Y[i] - d.Y[i]), kind * " th=" * cpp_to_string(th))
                end

                # C3: hinge opening angle
                for e in c.hinge_edges
                    ed = m.edges[e]
                    f1 = m.half_edges[ed.he[1]].face
                    f2 = m.half_edges[ed.he[2]].face
                    a = c.hinge_dir[e].src; b = c.hinge_dir[e].dst
                    d1 = Y[K.prime_vertex(c, f1, b)] - Y[K.prime_vertex(c, f1, a)]
                    d2 = Y[K.prime_vertex(c, f2, b)] - Y[K.prime_vertex(c, f2, a)]
                    ang = atan(d1[1] * d2[2] - d1[2] * d2[1], dot(d1, d2))
                    hit!(c3, abs(abs(ang) - th), kind)
                end
                # C4: split duplicates identical vectors
                for e in c.split_edges
                    ed = m.edges[e]
                    f1 = m.half_edges[ed.he[1]].face
                    f2 = m.half_edges[ed.he[2]].face
                    a = ed.key.a; b = ed.key.b
                    d1 = Y[K.prime_vertex(c, f1, b)] - Y[K.prime_vertex(c, f1, a)]
                    d2 = Y[K.prime_vertex(c, f2, b)] - Y[K.prime_vertex(c, f2, a)]
                    hit!(c4, norm(d1 - d2), kind)
                end
                # C5: face signed areas constant
                for f in 1:K.n_faces(m)
                    Pf = [Y[K.prime_vertex(c, f, v)] for v in m.faces[f]]
                    Pflat = [X[v] for v in m.faces[f]]
                    hit!(c5, abs(signed_area(Pf) - signed_area(Pflat)), kind)
                end
                # C6: pin equation for the velocity field (omega_f, w_f) at this theta
                #     omega_f = -sigma_f/2,  t_f = 2 s J u_f,  w_f = dt_f/dth - omega_f J t_f
                #                                        = c J u_f - sigma_f s u_f
                for e in c.hinge_edges
                    ed = m.edges[e]
                    f = m.half_edges[ed.he[1]].face
                    g = m.half_edges[ed.he[2]].face
                    p = Y[K.prime_vertex(c, f, c.hinge_dir[e].src)]
                    wf = cc * (J * u[f]) - Float64(m.sigma[f]) * ss * u[f]
                    wg = cc * (J * u[g]) - Float64(m.sigma[g]) * ss * u[g]
                    om_f = -0.5 * m.sigma[f]; om_g = -0.5 * m.sigma[g]
                    hit!(c6, norm((wf - wg) + (om_f - om_g) * (J * p)), kind)
                end
                # C7: trajectory on the centred conic y^T (A A^T)^-1 y = 1, A = [C S]
                if th > 0
                    for f in 1:K.n_faces(m)
                        for v in m.faces[f]
                            pv = K.prime_vertex(c, f, v)
                            Sv = J * (2.0 * u[f] - Float64(m.sigma[f]) * X[v])
                            A = SMatrix{2,2,Float64}(X[v][1], X[v][2], Sv[1], Sv[2])
                            dt = det(A)
                            abs(dt) < 1e-8 && continue   # degenerate ellipse, skip
                            Q = inv(A * transpose(A))
                            hit!(c7, abs(dot(Y[pv], Q * Y[pv]) - 1.0), kind)
                        end
                    end
                end
            end
        end
    end

    @printf("graphs solved: %d, deployable X0 used: %d\n", n_graphs, n_deployable)
    rep(n, w) = @printf("  %-58s max = %.3e   (%s)\n", n, w.v, w.where)
    rep("C1 potential closes on EVERY hinge edge (path-independence)", c1)
    rep("C2 Y = cos(th/2) C + sin(th/2) S vs deploy()", c2)
    rep("C3 |hinge opening angle| - theta", c3)
    rep("C4 split-edge duplicates are the same vector", c4)
    rep("C5 face signed area constant in theta", c5)
    rep("C6 pin equation for (omega,w) = no-locking, A(Y_th) sigma = 0", c6)
    rep("C7 vertex trajectory on the centred conic through C, S", c7)
end

# std::to_string(double) prints with "%f" (6 decimals).
cpp_to_string(x::Float64) = @sprintf("%f", x)

main()
