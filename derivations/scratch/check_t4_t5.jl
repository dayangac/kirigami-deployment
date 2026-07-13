# check_t4_t5.jl -- Deriver's own numeric checks for derivations/core.md, T4 and T5.
# Port of check_t4_t5.cpp (same corpus, same seed, same printed lines).
#
# Deliberately re-implements the harmonic calculus from scratch (it does NOT use
# Kirigami/src/method/deploy_basis.jl, which the Experimenter owns) so that agreement is
# an independent cross-check rather than a tautology.
#
# Run (from repo root):
#   julia --project=Kirigami derivations/scratch/check_t4_t5.jl
#
# D1  exact theta_c by closed-form enumeration over ALL ordered (vertex, edge) pairs
#     with the two interval harmonics, vs collision.jl's grid+bisection theta_max
# D2  on split-free patterns, theta_c == min_i beta_i = min_i (2pi - alpha_i - alpha_j)
#     (2026 Sec. 5.2 / 2025 Sec. 4.2 falls out of the candidate list as a special case)
# D3  split-edge separation: closed-form p = -sigma_f det(d, du), q = -p,
#     r = det(d, J du) against the generic (T3.2) coefficients, and p + q == 0
# D4  is max(|C|,|S|) an upper bound for |Y(theta)| on [0, pi)?  (K2c's rho_f as
#     specified).  Reports the worst overshoot ratio.

using Kirigami
using LinearAlgebra
using Printf
using StaticArrays
const K = Kirigami
const Vec2 = K.Vec2

const J = SMatrix{2,2,Float64}(0, 1, -1, 0)   # [0 -1; 1 0]
cross2(a::Vec2, b::Vec2) = a[1] * b[2] - a[2] * b[1]

mutable struct Harm
    p::Float64
    q::Float64
    r::Float64
end
Harm() = Harm(0.0, 0.0, 0.0)
heval(h::Harm, th::Float64) = h.p + h.q * cos(th) + h.r * sin(th)
hscale(h::Harm) = abs(h.p) + hypot(h.q, h.r)

# (T3.2): det(Yb - Ya, Yw - Ya)
function orient_h(Cab::Vec2, Sab::Vec2, Caw::Vec2, Saw::Vec2)
    dcc = cross2(Cab, Caw); dss = cross2(Sab, Saw)
    Harm(0.5 * (dcc + dss), 0.5 * (dcc - dss), 0.5 * (cross2(Cab, Saw) + cross2(Sab, Caw)))
end
# (T3.3): <Yb - Ya, Yw - Ya>
function dot_h(Cab::Vec2, Sab::Vec2, Caw::Vec2, Saw::Vec2)
    dcc = dot(Cab, Caw); dss = dot(Sab, Saw)
    Harm(0.5 * (dcc + dss), 0.5 * (dcc - dss), 0.5 * (dot(Cab, Saw) + dot(Sab, Caw)))
end

# (T3.5)-(T3.6): roots of p + q cos th + r sin th in (lo, hi], ascending.
function roots_in(h::Harm, lo::Float64, hi::Float64, tol::Float64)
    out = Float64[]
    hscale(h) <= tol && return out   # identically zero to tolerance -> degenerate
    A = h.p - h.q; B = 2 * h.r; Cc = h.p + h.q
    taus = Float64[]
    if abs(A) <= 1e-14 * (abs(B) + abs(Cc) + 1e-300)
        abs(B) > 0 && push!(taus, -Cc / B)   # the other root is at theta = pi
        push!(out, pi)                        # tau = infinity
    else
        disc = B * B - 4 * A * Cc
        disc < 0 && return out
        sq = sqrt(disc)
        push!(taus, (-B + sq) / (2 * A))
        push!(taus, (-B - sq) / (2 * A))
    end
    for t in taus
        push!(out, 2.0 * atan(t))
    end
    keep = [th for th in out if th > lo && th <= hi]
    sort!(keep)
    return keep
end

mutable struct ExactRange
    theta_c::Float64              # first CONTACT angle (smallest admissible root)
    theta_max::Float64            # first INTERIOR-OVERLAP angle
    candidates::Vector{Float64}   # all admissible roots in (0, pi], ascending, deduped
end

# T4: collect every admissible contact root over ordered (M'-vertex w of face g,
# edge (a,b) of face f != g): a root of the orientation harmonic in (0, pi] whose
# two interval harmonics place w on the CLOSED edge.  Lemma T4.1 says the pairwise
# interior-overlap status of the structure can only change at such an angle, so a
# midpoint test between consecutive candidates gives the exact theta_max.
function exact_range(c::K.CutStructure, X::Vector{Vec2}, C::Vector{Vec2}, S::Vector{Vec2},
                     geom_scale::Float64)
    m = c.mesh
    F = K.n_faces(m)
    tol = 1e-11 * geom_scale * geom_scale
    out = ExactRange(pi, pi, Float64[])
    for f in 1:F
        pf = c.prime_faces[f]
        nf = length(pf)
        for k in 1:nf
            a = pf[k]; b = pf[k % nf + 1]
            Cab = C[b] - C[a]; Sab = S[b] - S[a]
            len2 = dot(Cab, Cab)   # constant in theta (T3.3 corollary)
            len2 <= 0 && continue
            for g in 1:F
                g == f && continue
                for w in c.prime_faces[g]
                    Caw = C[w] - C[a]; Saw = S[w] - S[a]
                    ho = orient_h(Cab, Sab, Caw, Saw)
                    hd = dot_h(Cab, Sab, Caw, Saw)
                    for th in roots_in(ho, 1e-9, Float64(pi), tol)
                        d = heval(hd, th)
                        eps = 1e-9 * len2
                        (d < -eps || d > len2 + eps) && continue   # CLOSED edge, endpoints included
                        push!(out.candidates, th)
                    end
                end
            end
        end
    end
    sort!(out.candidates)
    ded = Float64[]
    for t in out.candidates
        (isempty(ded) || t - ded[end] > 1e-9) && push!(ded, t)
    end
    out.candidates = ded
    out.theta_c = isempty(out.candidates) ? Float64(pi) : out.candidates[1]
    # midpoint scan
    out.theta_max = pi
    SHRINK = 1e-7
    # the first interval is (0, cand[1]): an X that is positively oriented but NOT
    # embedded already overlaps there, and theta_max = 0.
    let hi0 = isempty(out.candidates) ? Float64(pi) : out.candidates[1]
        d = K.deploy(c, X, 0.5 * hi0)
        if K.has_collision(c, d.Y, SHRINK)
            out.theta_max = 0.0
            return out
        end
    end
    for i in eachindex(out.candidates)
        hi = i < length(out.candidates) ? out.candidates[i + 1] : Float64(pi)
        hi <= out.candidates[i] + 1e-12 && continue
        mid = 0.5 * (out.candidates[i] + hi)
        d = K.deploy(c, X, mid)
        if K.has_collision(c, d.Y, SHRINK)
            out.theta_max = out.candidates[i]
            break
        end
    end
    return out
end

function main()
    rng = K.MT19937(20260903)
    cases = [("squares", [2.5]), ("triangles", [2.2]), ("hexagons", [2.5]),
             ("kagome", [2.2]), ("snub_square", [2.2]), ("truncated_square", [2.5]),
             ("t3_4_3_12", [3.0]), ("delaunay", [40.0, 8.0]), ("voronoi", [35.0, 8.0]),
             ("quad_random", [40.0, 8.0])]

    @printf("%-18s %5s %5s %5s %6s  %-10s %-10s %-10s %-10s %-8s\n", "graph", "F", "hinge",
            "split", "cands", "theta_c", "theta_max", "bisect", "min_beta", "d1")
    d1_worst = 0.0; d2_worst = 0.0; d3_worst = 0.0; d4_worst_ratio = 1.0
    d3_pairs = 0; d4_violations = 0; d4_total = 0
    n_used = 0; n_split_free = 0

    for (kind, par) in cases
        for rep in 1:2
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
            K.deployable(K.hole_residuals(c, X, hs), 1e-7) || continue
            K.n_faces(m) > 260 && continue
            # D1/D2 are statements about VALID embeddings; X0 is self-intersecting on random
            # graphs (F17), where "first contact" is at theta = 0 and both numbers are 0.
            positively_oriented = true
            for f in 1:K.n_faces(m)
                P = [X[v] for v in m.faces[f]]
                A2 = 0.0
                for i in eachindex(P)
                    a = P[i]; b = P[i % length(P) + 1]
                    A2 += a[1] * b[2] - a[2] * b[1]
                end
                if A2 <= 0
                    positively_oriented = false
                    break
                end
            end
            if !positively_oriented
                @printf("%-18s  SKIPPED (X0 has inverted faces, F17)\n", kind)
                continue
            end
            n_used += 1

            # C and S from T1.5. C = Y(0); S = 2 dY/dtheta|_0 -- taken from the code's
            # analytic derivative so that this program's basis is independent of
            # check_t1_t2.jl's BFS.
            d0 = K.deploy(c, X, 0.0)
            C = d0.Y
            S = [2.0 * d0.dY_dtheta[i] for i in 1:c.n_prime_vertices]

            scale = 0.0
            for p in X
                scale = max(scale, norm(p))
            end

            er = exact_range(c, X, C, S, scale)
            tm = K.theta_max(c, X)
            d1 = abs(er.theta_max - tm.theta_max_geometric)
            d1_worst = max(d1_worst, d1)

            # D2: split-free patterns
            beta_min = tm.min_beta
            if K.n_split(c) == 0
                n_split_free += 1
                d2_worst = max(d2_worst, abs(min(er.theta_max, pi) - min(beta_min, pi)))
            end

            # D3: split-edge separation closed form
            # Recover u from S: S_(v,f) = J(2 u_f - sigma_f x_v) => u_f = (J^-1 S + sigma_f x_v)/2
            u = fill(Vec2(0, 0), K.n_faces(m))
            for f in 1:K.n_faces(m)
                v = m.faces[f][1]
                pv = K.prime_vertex(c, f, v)
                u[f] = 0.5 * (transpose(J) * S[pv] + Float64(m.sigma[f]) * X[v])
            end
            for e in c.split_edges
                ed = m.edges[e]
                f = m.half_edges[ed.he[1]].face
                g = m.half_edges[ed.he[2]].face
                va = ed.key.a; vb = ed.key.b
                dvec = X[vb] - X[va]
                du = u[g] - u[f]
                pred = Harm()
                pred.p = -Float64(m.sigma[f]) * cross2(dvec, du)
                pred.q = -pred.p
                pred.r = cross2(dvec, J * du)
                # generic route: orientation harmonic of (a',b' in f ; a'' in g)
                A = K.prime_vertex(c, f, va); B = K.prime_vertex(c, f, vb)
                W = K.prime_vertex(c, g, va)
                gen = orient_h(C[B] - C[A], S[B] - S[A], C[W] - C[A], S[W] - S[A])
                sc = max(1.0, hscale(gen))
                d3_worst = max(d3_worst, abs(pred.p - gen.p) / sc,
                               abs(pred.q - gen.q) / sc, abs(pred.r - gen.r) / sc,
                               abs(gen.p + gen.q) / sc)
                d3_pairs += 1
            end

            # D4: is max(|C|,|S|) an upper bound for |Y(theta)| on [0, pi)?
            for i in 1:c.n_prime_vertices
                bound = max(norm(C[i]), norm(S[i]))
                bound < 1e-12 && continue
                d4_total += 1
                worst = 0.0
                for k in 0:64
                    th = pi * k / 65.0
                    worst = max(worst, norm(cos(th / 2) * C[i] + sin(th / 2) * S[i]))
                end
                if worst > bound * (1 + 1e-9)
                    d4_violations += 1
                    d4_worst_ratio = max(d4_worst_ratio, worst / bound)
                end
            end

            @printf("%-18s %5d %5d %5d %6d  %-10.6f %-10.6f %-10.6f %-10.6f %-8.1e\n",
                    kind, K.n_faces(m), K.n_hinge(c), K.n_split(c), length(er.candidates),
                    er.theta_c, er.theta_max, tm.theta_max_geometric, beta_min, d1)
        end
    end

    @printf("\ngraphs used: %d (split-free: %d)\n", n_used, n_split_free)
    @printf("D1  |theta_c(exact) - theta_max(bisection)|      worst = %.3e   (rule <= 1e-5)\n",
            d1_worst)
    @printf("D2  |theta_max - min_i beta_i| on split-free     worst = %.3e\n", d2_worst)
    @printf("D3  split-edge closed form vs generic (T3.2), and p+q == 0\n")
    @printf("      pairs = %d   worst relative deviation = %.3e\n", d3_pairs, d3_worst)
    @printf("D4  max(|C|,|S|) as a swept-radius bound: %d / %d copies VIOLATE it, worst ratio %.4f\n",
            d4_violations, d4_total, d4_worst_ratio)
end

main()
