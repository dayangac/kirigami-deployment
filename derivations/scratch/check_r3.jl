# check_r3.jl -- Deriver round 3.  Settles the ONE numeric question left by
# derivations/check.md R2.5: with the THIRD structural class g(0) = g'(0) = 0
# deflated as well, how many harmonics does the atom list get wrong?
# Port of check_r3.cpp (same corpus, same seed, same printed lines).
#
# Run (from repo root):
#   julia --project=Kirigami derivations/scratch/check_r3.jl 0.02
#
# The round-2 program compared the deflated atom list against a reference that used the
# SAME deflation rule, so it could not see the third class; check.md R2.5 is right about
# that.  Here the truth is computed WITHOUT the tau chart at all: h(theta) = p + q cos +
# r sin is a degree-1 trig polynomial, so it is monotone between its critical points
# (tan theta = r/q).  Split (0, eps) at those, and a CROSSING exists iff two consecutive
# break points carry strictly opposite signs.  That is an independent decision procedure
# for "is there a contact event in (0, eps)", and it is well posed on the third class,
# where h = p(1 - cos theta) has constant sign and therefore never crosses.

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

# ---------------------------------------------------------------------------
# The atom list of T5.2b.1 / T5.2b.2, now with THREE structural classes.
#   A = p - q,  B = 2r,  C = p + q,   g(tau) = C + B tau + A tau^2,  T = tan(eps/2).
const K_GENERIC = 1; const K_SIMPLE0 = 2; const K_DOUBLE0 = 3   # 1-based class index

function classify(h::Harm, tol::Float64)
    C = h.p + h.q; B = 2 * h.r
    abs(C) > tol && return K_GENERIC
    abs(B) > tol && return K_SIMPLE0
    return K_DOUBLE0
end

# returns true iff the atom list says "no admissible root in (0,T)".
function atoms_no_root(h::Harm, T::Float64, k::Int)
    A = h.p - h.q; B = 2 * h.r; C = h.p + h.q
    k == K_DOUBLE0 && return true                      # (T5.1e): g = A tau^2, root only at 0
    if k == K_SIMPLE0                                  # (T5.1d): g = tau (A tau + B)
        return !((A * B < 0) && (-A * B - A * A * T < 0))
    end
    g0 = C; gT = C + B * T + A * T * T                 # (T5.1c)
    (g0 * gT > 0) || return false
    disc = B * B - 4 * A * C
    both = (disc >= 0) && (A * g0 > 0) && (A * gT > 0) && (A * B < 0) &&
           (-A * B - 2 * A * A * T < 0)
    return !both
end

# ---------------------------------------------------------------------------
# Independent truth: does h CROSS zero somewhere in the open interval (0, eps)?
# Uses only evaluations of h and its critical points -- no tau chart, no quadratic.
# Returns  0 = no crossing, 1 = crossing, -1 = ambiguous (a break-point value sits in
# the noise band, so double precision cannot decide; those are reported separately).
function truth_crossing(h::Harm, eps::Float64, tol::Float64)
    # break points: 0+, the critical points of h inside (0, eps), and eps-.
    br = Float64[]
    base = atan(h.r, h.q)   # h' = -q sin + r cos = 0  <=>  tan th = r/q
    for k in -2:2
        th = base + k * pi
        (th > 0 && th < eps) && push!(br, th)
    end
    sort!(br)

    sgn = Int[]
    # sign of h just to the RIGHT of 0, from the Taylor coefficients:
    #   h(0) = C,  h'(0) = r = B/2,  h''(0) = -q = (A - C)/2 ... -> class-wise.
    let C = h.p + h.q, B = 2 * h.r, A = h.p - h.q
        s = abs(C) > tol ? C : (abs(B) > tol ? B : A)   # h = p(1 - cos th): sign of A on (0, pi)
        abs(s) <= tol && return -1
        push!(sgn, s > 0 ? 1 : -1)
    end
    for th in br
        v = heval(h, th)
        abs(v) <= tol && return -1
        push!(sgn, v > 0 ? 1 : -1)
    end
    let v = heval(h, eps)
        abs(v) <= tol && return -1     # root at (or within noise of) the endpoint
        push!(sgn, v > 0 ? 1 : -1)
    end
    for i in 2:length(sgn)
        sgn[i] != sgn[i - 1] && return 1
    end
    return 0
end

# ---------------------------------------------------------------------------
mutable struct Tally
    tested::Int
    ambiguous::Int
    n_class::Vector{Int}
    mm3::Int                    # mismatches, three-class atom list
    mm3_by_class::Vector{Int}
    mm2::Int                    # mismatches, round-2 (two-class) atom list
    mm2_by_class::Vector{Int}
end
Tally() = Tally(0, 0, zeros(Int, 3), 0, zeros(Int, 3), 0, zeros(Int, 3))

# round-2 rule: only the C = 0 deflation, no third class.
atoms_no_root_r2(h::Harm, T::Float64, k::Int) = atoms_no_root(h, T, k == K_DOUBLE0 ? K_SIMPLE0 : k)

function scan_harmonics(c::K.CutStructure, C::Vector{Vec2}, S::Vector{Vec2},
                        geom_scale::Float64, eps::Float64, tal::Tally)
    m = c.mesh
    F = K.n_faces(m)
    tol = 1e-11 * geom_scale * geom_scale
    T = tan(0.5 * eps)
    for f in 1:F
        pf = c.prime_faces[f]
        nf = length(pf)
        for k in 1:nf
            a = pf[k]; b = pf[k % nf + 1]
            Cab = C[b] - C[a]; Sab = S[b] - S[a]
            dot(Cab, Cab) <= 0 && continue
            for g in 1:F
                g == f && continue
                for w in c.prime_faces[g]
                    Caw = C[w] - C[a]; Saw = S[w] - S[a]
                    ho = orient_h(Cab, Sab, Caw, Saw)
                    hscale(ho) <= tol && continue
                    tr = truth_crossing(ho, eps, tol)
                    kl = classify(ho, tol)
                    tal.n_class[kl] += 1
                    if tr < 0
                        tal.ambiguous += 1
                        continue
                    end
                    tal.tested += 1
                    truth_noroot = (tr == 0)
                    if atoms_no_root(ho, T, kl) != truth_noroot
                        tal.mm3 += 1; tal.mm3_by_class[kl] += 1
                    end
                    if atoms_no_root_r2(ho, T, kl) != truth_noroot
                        tal.mm2 += 1; tal.mm2_by_class[kl] += 1
                    end
                end
            end
        end
    end
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

function main()
    EPS = length(ARGS) > 0 ? parse(Float64, ARGS[1]) : 0.02
    rng = K.MT19937(20260904)   # same seed / same corpus as check_r2.jl

    cases = [("squares", [2.5]), ("triangles", [2.2]), ("hexagons", [2.5]),
             ("kagome", [2.2]), ("snub_square", [2.2]), ("truncated_square", [2.5]),
             ("t3_4_3_12", [3.0]), ("delaunay", [40.0, 8.0]), ("voronoi", [35.0, 8.0]),
             ("quad_random", [40.0, 8.0])]

    tal = Tally()
    n_used = 0; n_samples = 0

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
                n_samples += 1
                ds = K.deploy(c, Xs, 0.0)
                Cs = ds.Y
                Ss = [2.0 * ds.dY_dtheta[i] for i in 1:c.n_prime_vertices]
                sc = 0.0
                for p in Xs
                    sc = max(sc, norm(p))
                end
                scan_harmonics(c, Cs, Ss, sc, EPS, tal)
            end
        end
    end

    @printf("\nR3  eps = %.4f   graphs used: %d   shape-space samples: %d\n", EPS, n_used,
            n_samples)
    @printf("    harmonics decided: %d    ambiguous (skipped): %d\n", tal.tested, tal.ambiguous)
    @printf("    class sizes:  g(0) != 0 : %d    g(0)=0, g'(0) != 0 : %d    g(0)=g'(0)=0 : %d\n",
            tal.n_class[1], tal.n_class[2], tal.n_class[3])
    @printf("    round-2 atom list (two classes):   %d mismatches  [%d / %d / %d by class]\n",
            tal.mm2, tal.mm2_by_class[1], tal.mm2_by_class[2], tal.mm2_by_class[3])
    @printf("    round-3 atom list (three classes): %d mismatches  [%d / %d / %d by class]\n",
            tal.mm3, tal.mm3_by_class[1], tal.mm3_by_class[2], tal.mm3_by_class[3])
end

main()
