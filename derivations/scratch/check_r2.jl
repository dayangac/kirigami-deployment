# check_r2.jl -- Deriver round 2.  Decides the five numeric questions raised by
# derivations/check.md (D1, D2, D3, D4, D5) with a program written for this purpose.
# Port of check_r2.cpp (same corpus, same seed, same printed lines).
#
# Run (from repo root):
#   julia --project=Kirigami derivations/scratch/check_r2.jl [eps]      (default eps = 0.02)
#
# R2-A (D1)  In the face's OWN moving-centroid frame -- the frame ideas/ranking.md K2c
#            specifies -- is max(|x|,|chi|) an exact swept radius, or a sqrt(2) under-
#            estimate?  Reports max_theta |y_u(theta) - gamma_f(theta)| against both
#            max(|x|,|chi|) and the flat circumradius |x_u - xbar_f|.
# R2-B (D2)  Does the swept region stay within O(r_f) of the FLAT centroid (the form the
#            O(n) packing count needs)?  Growing square patch, reports
#            max_f max_theta |gamma_f(theta) - xbar_f| / r_f against patch diameter.
# R2-C (D3)  Split-free patterns: Theta_max == min(min_e beta_e, pi)?
# R2-D (D4)  Proposition T5.2b' (the inner-approximation direction).  For samples in the
#            shape space: (all A_f > 0) AND (not penetrating at 0+) AND (no orientation
#            harmonic of the COMPLETE candidate list has a root in (0,eps))  ==>
#            Theta_max >= eps.  Counts hypothesis hits and violations.
# R2-E (D5)  The degree-<= 4 atom list for "the tau-quadratic has no root in (0,T)",
#            checked against direct root computation on every harmonic encountered.

using Kirigami
using LinearAlgebra
using Printf
using StaticArrays
const K = Kirigami
const Vec2 = K.Vec2

cross2(a::Vec2, b::Vec2) = a[1] * b[2] - a[2] * b[1]

struct Harm
    p::Float64
    q::Float64
    r::Float64
end
heval(h::Harm, th::Float64) = h.p + h.q * cos(th) + h.r * sin(th)
hscale(h::Harm) = abs(h.p) + hypot(h.q, h.r)

function orient_h(Cab::Vec2, Sab::Vec2, Caw::Vec2, Saw::Vec2)
    dcc = cross2(Cab, Caw); dss = cross2(Sab, Saw)
    Harm(0.5 * (dcc + dss), 0.5 * (dcc - dss), 0.5 * (cross2(Cab, Saw) + cross2(Sab, Caw)))
end
function dot_h(Cab::Vec2, Sab::Vec2, Caw::Vec2, Saw::Vec2)
    dcc = dot(Cab, Caw); dss = dot(Sab, Saw)
    Harm(0.5 * (dcc + dss), 0.5 * (dcc - dss), 0.5 * (dot(Cab, Saw) + dot(Sab, Caw)))
end

# Roots of p + q cos th + r sin th in (lo, hi], ascending.  (T3.5)/(T3.6).
function roots_in(h::Harm, lo::Float64, hi::Float64, tol::Float64)
    out = Float64[]
    hscale(h) <= tol && return out
    A = h.p - h.q; B = 2 * h.r; Cc = h.p + h.q
    taus = Float64[]
    if abs(A) <= 1e-14 * (abs(B) + abs(Cc) + 1e-300)
        abs(B) > 0 && push!(taus, -Cc / B)
        push!(out, pi)
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

# ---------------------------------------------------------------------------
# R2-E: the degree-<= 4 atom list of check.md D5.
# g(tau) = C + B tau + A tau^2 with A = p-q, B = 2r, C = p+q.
# "both roots in (0,T)" :=  disc >= 0 AND A g(0) > 0 AND A g(T) > 0 AND A B < 0
#                           AND -A B - 2 A^2 T < 0      (i.e. vertex -B/(2A) < T)
# "no root in (0,T)"    :=  g(0) g(T) > 0  AND NOT["both roots in (0,T)"]
function no_root_in_atoms(h::Harm, T::Float64)
    A = h.p - h.q; B = 2 * h.r; C = h.p + h.q
    g0 = C; gT = C + B * T + A * T * T
    (g0 * gT > 0) || return false
    disc = B * B - 4 * A * C
    both = (disc >= 0) && (A * g0 > 0) && (A * gT > 0) && (A * B < 0) &&
           (-A * B - 2 * A * A * T < 0)
    return !both
end
# Deflated form: when g(0) = 0 identically (a permanent incidence at theta = 0,
# T3.H.2 / T5.3), divide out the tau factor and test the linear remainder.
function no_root_in_atoms_deflated(h::Harm, T::Float64, czero_tol::Float64)
    A = h.p - h.q; B = 2 * h.r; C = h.p + h.q
    if abs(C) <= czero_tol
        # g = tau (A tau + B); root -B/A in (0,T)?
        inside = (A * B < 0) && (-A * B - A * A * T < 0)
        return !inside
    end
    return no_root_in_atoms(h, T)
end
# Direct reference: does g have a root in the OPEN interval (0,T)?
function no_root_in_direct(h::Harm, T::Float64)
    A = h.p - h.q; B = 2 * h.r; C = h.p + h.q
    if abs(A) <= 1e-300
        abs(B) <= 1e-300 && return true   # constant, no finite root (C != 0 assumed)
        t = -C / B
        return !(t > 0 && t < T)
    end
    disc = B * B - 4 * A * C
    disc < 0 && return true
    sq = sqrt(disc)
    t1 = (-B + sq) / (2 * A); t2 = (-B - sq) / (2 * A)
    in1 = (t1 > 0 && t1 < T); in2 = (t2 > 0 && t2 < T)
    return !(in1 || in2)
end

# ---------------------------------------------------------------------------
mutable struct RangeInfo
    theta_max::Float64
    candidates::Vector{Float64}
    any_orient_root_below::Bool   # any orientation-harmonic root in (0, eps)
end

# The five R2-E tallies the C++ passed by pointer.
mutable struct AtomStats
    tests::Int
    mismatch::Int
    mm_czero::Int
    tests_d::Int
    mm_d::Int
end

# Complete candidate enumeration over ordered (M'-vertex w of g, edge (a,b) of f != g).
# `eps` only drives the extra flag; theta_max is by the T4.2" midpoint scan.
function exact_range(c::K.CutStructure, X::Vector{Vec2}, C::Vector{Vec2}, S::Vector{Vec2},
                     geom_scale::Float64, eps::Float64, st::AtomStats)
    m = c.mesh
    F = K.n_faces(m)
    tol = 1e-11 * geom_scale * geom_scale
    T = tan(0.5 * eps)
    out = RangeInfo(pi, Float64[], false)
    for f in 1:F
        pf = c.prime_faces[f]
        nf = length(pf)
        for k in 1:nf
            a = pf[k]; b = pf[k % nf + 1]
            Cab = C[b] - C[a]; Sab = S[b] - S[a]
            len2 = dot(Cab, Cab)
            len2 <= 0 && continue
            for g in 1:F
                g == f && continue
                for w in c.prime_faces[g]
                    Caw = C[w] - C[a]; Saw = S[w] - S[a]
                    ho = orient_h(Cab, Sab, Caw, Saw)
                    hd = dot_h(Cab, Sab, Caw, Saw)
                    if hscale(ho) > tol
                        # R2-E: atom list vs direct roots, on this very harmonic.
                        czt = 1e-11 * geom_scale * geom_scale
                        a1 = no_root_in_atoms(ho, T); a2 = no_root_in_direct(ho, T)
                        st.tests += 1
                        if a1 != a2
                            st.mismatch += 1
                            abs(ho.p + ho.q) <= czt && (st.mm_czero += 1)
                        end
                        a3 = no_root_in_atoms_deflated(ho, T, czt)
                        st.tests_d += 1
                        # reference for the deflated test: roots of the deflated polynomial
                        ref = let A = ho.p - ho.q, B = 2 * ho.r, Cz = ho.p + ho.q
                            if abs(Cz) <= czt
                                if abs(A) <= 1e-300
                                    true
                                else
                                    t = -B / A
                                    !(t > 0 && t < T)
                                end
                            else
                                a2
                            end
                        end
                        a3 != ref && (st.mm_d += 1)
                        # the D4 hypothesis is about ORIENTATION roots only (interval test dropped),
                        # with the theta = 0 incidence deflated away.
                        ref || (out.any_orient_root_below = true)
                    end
                    for th in roots_in(ho, 1e-9, Float64(pi), tol)
                        d = heval(hd, th)
                        e2 = 1e-9 * len2
                        (d < -e2 || d > len2 + e2) && continue
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

    SHRINK = 1e-7
    let hi0 = isempty(out.candidates) ? Float64(pi) : out.candidates[1]
        d = K.deploy(c, X, 0.5 * hi0)
        if K.has_collision(c, d.Y, SHRINK)
            out.theta_max = 0.0
            return out
        end
    end
    out.theta_max = pi
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

function positively_oriented(m::K.Mesh, X::Vector{Vec2})
    for f in 1:K.n_faces(m)
        A2 = 0.0
        Fv = m.faces[f]
        for i in eachindex(Fv)
            a = X[Fv[i]]; b = X[Fv[i % length(Fv) + 1]]
            A2 += a[1] * b[2] - a[2] * b[1]
        end
        A2 <= 0 && return false
    end
    return true
end

# "no interior overlap at theta = 0+": overlap probe just above zero, below the first
# candidate.  This is exactly the negation of contact.jl::penetrates_immediately().
function embedded_at_zero_plus(c::K.CutStructure, X::Vector{Vec2}, first_candidate::Float64)
    probe = min(1e-6, 0.5 * first_candidate)
    d = K.deploy(c, X, probe)
    return !K.has_collision(c, d.Y, 1e-9)
end

function main()
    rng = K.MT19937(20260904)
    EPS = length(ARGS) > 0 ? parse(Float64, ARGS[1]) : 0.02

    st = AtomStats(0, 0, 0, 0, 0)

    # -----------------------------------------------------------------  R2-A / R2-C
    cases = [("squares", [2.5]), ("triangles", [2.2]), ("hexagons", [2.5]),
             ("kagome", [2.2]), ("snub_square", [2.2]), ("truncated_square", [2.5]),
             ("t3_4_3_12", [3.0]), ("delaunay", [40.0, 8.0]), ("voronoi", [35.0, 8.0]),
             ("quad_random", [40.0, 8.0])]

    a_worst_specform = 0.0   # | max_theta |y-gamma| - max(|x|,|chi|) |
    a_worst_circum = 0.0     # | max_theta |y-gamma| - |x_u - xbar_f| |
    a_copies = 0
    c_worst = 0.0
    n_split_free = 0; n_used = 0

    # -----------------------------------------------------------------  R2-D
    d_hyp = 0; d_viol = 0; d_samples = 0; d_short = 0
    d_worst_viol = 0.0

    @printf("%-18s %5s %6s  %-10s %-10s %-10s\n", "graph", "F", "cands", "Theta_max",
            "min_beta", "min(b,pi)")

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
            K.n_faces(m) > 260 && continue
            X = K.matrix_to_points(sr.X0)
            K.deployable(K.hole_residuals(c, X, hs), 1e-7) || continue
            positively_oriented(m, X) || continue
            n_used += 1

            scale = 0.0
            for p in X
                scale = max(scale, norm(p))
            end

            d0 = K.deploy(c, X, 0.0)
            C = d0.Y
            S = [2.0 * d0.dY_dtheta[i] for i in 1:c.n_prime_vertices]

            # ---- R2-A: swept radius in the face's own moving-centroid frame.
            for f in 1:K.n_faces(m)
                pf = c.prime_faces[f]
                gc = Vec2(0, 0); gs = Vec2(0, 0); xbar = Vec2(0, 0)
                for i in eachindex(pf)
                    gc += C[pf[i]]; gs += S[pf[i]]
                end
                gc /= Float64(length(pf))
                gs /= Float64(length(pf))
                for v in m.faces[f]
                    xbar += X[v]
                end
                xbar /= Float64(length(m.faces[f]))
                for u in pf
                    x = C[u] - gc; chi = S[u] - gs
                    specform = max(norm(x), norm(chi))
                    worst = 0.0
                    for k in 0:128
                        th = pi * k / 128.0
                        worst = max(worst, norm(cos(th / 2) * x + sin(th / 2) * chi))
                    end
                    a_worst_specform = max(a_worst_specform, abs(worst - specform))
                    a_worst_circum = max(a_worst_circum, abs(worst - norm(x)))
                    a_copies += 1
                end
            end

            ri = exact_range(c, X, C, S, scale, EPS, st)
            tm = K.theta_max(c, X)

            # ---- R2-C: split-free  Theta_max == min(min beta, pi)
            if K.n_split(c) == 0
                n_split_free += 1
                c_worst = max(c_worst, abs(ri.theta_max - min(tm.min_beta, pi)))
                @printf("%-18s %5d %6d  %-10.6f %-10.6f %-10.6f\n", kind, K.n_faces(m),
                        length(ri.candidates), ri.theta_max, tm.min_beta,
                        min(tm.min_beta, pi))
            end

            # ---- R2-D: Proposition T5.2b' on X0 and on null-space samples.
            gauss = K.NormalDist(0.0, 1.0)   # fresh per graph, as the C++ constructs it
            n_samp = sr.dim_null > 0 ? 24 : 1
            for s in 0:n_samp-1
                Xs = X
                if s > 0
                    t = Matrix{Float64}(undef, sr.dim_null, 2)
                    for i in 1:sr.dim_null, j in 1:2
                        t[i, j] = 0.02 * scale * K.normal(gauss, rng)
                    end
                    Xs = K.matrix_to_points(sr.X0 + sr.Phi * t)
                end
                positively_oriented(m, Xs) || continue
                d_samples += 1
                ds = K.deploy(c, Xs, 0.0)
                Cs = ds.Y
                Ss = [2.0 * ds.dY_dtheta[i] for i in 1:c.n_prime_vertices]
                sc = 0.0
                for p in Xs
                    sc = max(sc, norm(p))
                end
                rs = exact_range(c, Xs, Cs, Ss, sc, EPS, st)
                firstc = isempty(rs.candidates) ? Float64(pi) : rs.candidates[1]
                emb0 = embedded_at_zero_plus(c, Xs, firstc)
                rs.theta_max < EPS - 1e-9 && (d_short += 1)   # samples the implication must exclude
                (rs.any_orient_root_below || !emb0) && continue   # hypothesis not met
                d_hyp += 1
                if rs.theta_max < EPS - 1e-9
                    d_viol += 1
                    d_worst_viol = max(d_worst_viol, EPS - rs.theta_max)
                end
            end
        end
    end

    # -----------------------------------------------------------------  R2-B
    @printf("\nR2-B  growing square patch: drift of the moving centroid from the FLAT one\n")
    @printf("%8s %5s %10s %14s %14s\n", "clip", "F", "diameter",
            "max|gam-xbar|/r", "max rho_spec/r")
    for R in (1.6, 2.6, 3.6, 4.6, 5.6, 7.0)
        m = K.tiling_squares(K.rect(Vec2(0, 0), R, R))
        K.build_topology!(m)
        r2 = K.MT19937(7)
        orep = K.assign_orientation_relaxation(m, r2)
        m.sigma = orep.sigma
        K.build_topology!(m)
        c = K.make_cut(m)
        hs = K.holes_partition(c)
        sys = K.assemble_system(c, hs, m.X, K.Fixed)
        sr = K.solve_system(sys, m.X)
        sr.projection_ok || continue
        X = K.matrix_to_points(sr.X0)
        d0 = K.deploy(c, X, 0.0)
        C = d0.Y
        S = [2.0 * d0.dY_dtheta[i] for i in 1:c.n_prime_vertices]
        diam = 0.0
        for i in eachindex(X), j in i+1:length(X)
            diam = max(diam, norm(X[i] - X[j]))
        end
        worst_drift = 0.0; worst_spec = 0.0
        for f in 1:K.n_faces(m)
            pf = c.prime_faces[f]
            gc = Vec2(0, 0); gs = Vec2(0, 0); xbar = Vec2(0, 0)
            for u in pf
                gc += C[u]; gs += S[u]
            end
            gc /= Float64(length(pf))
            gs /= Float64(length(pf))
            for v in m.faces[f]
                xbar += X[v]
            end
            xbar /= Float64(length(m.faces[f]))
            rf = 0.0
            for v in m.faces[f]
                rf = max(rf, norm(X[v] - xbar))
            end
            rf <= 0 && continue
            drift = 0.0
            for k in 0:64
                th = pi * k / 64.0
                gam = cos(th / 2) * gc + sin(th / 2) * gs
                drift = max(drift, norm(gam - xbar))
            end
            worst_drift = max(worst_drift, drift / rf)
            spec = 0.0
            for u in pf
                spec = max(spec, max(norm(C[u] - gc), norm(S[u] - gs)))
            end
            worst_spec = max(worst_spec, spec / rf)
        end
        @printf("%8.1f %5d %10.3f %14.3f %14.6f\n", R, K.n_faces(m), diam, worst_drift,
                worst_spec)
    end

    @printf("\n== results ==\n")
    @printf("R2-A  |max_theta|y-gamma| - max(|x|,|chi|)|  worst = %.3e   (copies %d)\n",
            a_worst_specform, a_copies)
    @printf("R2-A  |max_theta|y-gamma| - |x_u - xbar_f||  worst = %.3e\n", a_worst_circum)
    @printf("R2-C  split-free |Theta_max - min(min beta, pi)|  worst = %.3e  (%d patterns)\n",
            c_worst, n_split_free)
    @printf("R2-D  T5.2b' hypothesis met on %d / %d samples; VIOLATIONS = %d (worst shortfall %.3e), eps = %.4f; %d / %d samples have Theta_max < eps\n",
            d_hyp, d_samples, d_viol, d_worst_viol, EPS, d_short, d_samples)
    @printf("R2-E  D5 atom list as printed vs direct roots: %d tests, %d mismatches (%d of them with g(0) = 0)\n",
            st.tests, st.mismatch, st.mm_czero)
    @printf("R2-E  DEFLATED atom list vs deflated roots:  %d tests, %d mismatches\n",
            st.tests_d, st.mm_d)
    @printf("graphs used: %d\n", n_used)
end

main()
