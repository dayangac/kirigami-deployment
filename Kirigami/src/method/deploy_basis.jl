# method/deploy_basis.jl -- the trig-linear form of the forward kinematics and the
# harmonic (p + q cos(theta) + r sin(theta)) calculus built on top of it.
#
#   Y_pv(theta) = cos(theta/2) * C_pv + sin(theta/2) * S_pv,     C, S linear in X.
#
# This is claim U3/U6 of STATE.md, here derived rather than fitted: with
# R_f = rot(-sigma_f theta/2) = cos(theta/2) I - sigma_f sin(theta/2) J and the BFS
# translations t_f satisfying t_g - t_f = (R_f - R_g) x_v = (sigma_g - sigma_f)
# sin(theta/2) J x_v, every translation is sin(theta/2) times a constant vector, so
# C_pv = x_v and S_pv = -sigma_f J x_v + u_f. Numerically C = Y(0) and S = 2 dY/dtheta|_0.
#
# Consequence used everywhere below: for M'-vertices a, b, p every one of
#   det(Y_b - Y_a, Y_p - Y_a),  dot(Y_b - Y_a, Y_p - Y_a),  |Y_b - Y_a|^2
# is exactly p + q cos(theta) + r sin(theta) -- a "harmonic" in this file's language.
#
# Naming: the harmonic operations are the free functions `harmonic_eval(h, th)`, `amp(h)`,
# `scale(h)` (`eval` is reserved in a Julia module) and the basis columns are
# `basis_c(b, i)/basis_s(b, i)`; `classify_harmonic_A` and `harmonic_fit` return their
# secondary outputs (`A_zero`; `max_res`, `scale`) as tuples.

struct DeployBasis
    C::Matrix{Float64}  # n_prime x 2
    S::Matrix{Float64}  # n_prime x 2
end
n_prime(b::DeployBasis) = size(b.C, 1)
basis_c(b::DeployBasis, i::Int) = Vec2(b.C[i, 1], b.C[i, 2])   # row i of C
basis_s(b::DeployBasis, i::Int) = Vec2(b.S[i, 1], b.S[i, 2])   # row i of S

# Bit-faithful scalar kernels (docs/NUMERICS.md "Trig and reductions"): the inline cross
# product `u.x*v.y - u.y*v.x` is an fma, while dot()/squaredNorm()/norm() are UNFUSED
# reductions. The tie-sensitive predicates of contact.jl are decided by the last ulp of
# these, so they are written out.
det2(u::Vec2, v::Vec2) = fma(u[1], v[2], -(u[2] * v[1]))
_dotu(a::Vec2, b::Vec2) = a[1] * b[1] + a[2] * b[2]      # unfused dot
_sqnormu(a::Vec2) = a[1] * a[1] + a[2] * a[2]            # unfused squared norm
_normu(a::Vec2) = sqrt(_sqnormu(a))
if Sys.isapple()
    libm_atan(x::Float64) = ccall((:atan, _LIBM), Float64, (Float64,), x)
else
    libm_atan(x::Float64) = atan(x)
end

"""C = Y(0), S = 2 dY/dtheta|_0, from a single call to deploy()."""
function deploy_basis(c::CutStructure, X::Vector{Vec2})
    d = deploy(c, X, 0.0)
    n = length(d.Y)
    C = Matrix{Float64}(undef, n, 2)
    S = Matrix{Float64}(undef, n, 2)
    for i in 1:n
        C[i, 1] = d.Y[i][1]
        C[i, 2] = d.Y[i][2]
        S[i, 1] = 2.0 * d.dY_dtheta[i][1]
        S[i, 2] = 2.0 * d.dY_dtheta[i][2]
    end
    return DeployBasis(C, S)
end

"""Y(theta) from the basis."""
function basis_eval(b::DeployBasis, theta::Real)
    ss, cc = libm_sincos(0.5 * Float64(theta))   # cos/sin pair -> __sincos_stret
    # `cc * C + ss * S`: unfused per component
    return [Vec2(cc * b.C[i, 1] + ss * b.S[i, 1], cc * b.C[i, 2] + ss * b.S[i, 2]) for i in 1:n_prime(b)]
end

"""h(theta) = p + q cos(theta) + r sin(theta)."""
struct Harmonic
    p::Float64
    q::Float64
    r::Float64
end
Harmonic() = Harmonic(0.0, 0.0, 0.0)
# `p + q*cos + r*sin`: the cos/sin pair is one __sincos_stret call, and the sum is
# contracted left to right: fma(r, sin, fma(q, cos, p)).
function harmonic_eval(h::Harmonic, th::Real)
    s, c = libm_sincos(Float64(th))
    return fma(h.r, s, fma(h.q, c, h.p))
end
# hypot / atan2 / acos are Apple libm (as in generators.jl); the
# sub-ulp differences to Julia's own decide whether the tau = 0 artefact root of a class-2
# harmonic lands at +1e-17 or -1e-17, i.e. inside or outside (0, eps].
if Sys.isapple()
    libm_hypot(x::Float64, y::Float64) = ccall((:hypot, _LIBM), Float64, (Float64, Float64), x, y)
    libm_acos(x::Float64) = ccall((:acos, _LIBM), Float64, (Float64,), x)
else
    libm_hypot(x::Float64, y::Float64) = hypot(x, y)
    libm_acos(x::Float64) = acos(x)
end
amp(h::Harmonic) = libm_hypot(h.q, h.r)
scale(h::Harmonic) = abs(h.p) + amp(h)

# Same three, from explicit (C,S) pairs -- used by the range optimizer, which
# needs the coefficients as functions of the shape-space coefficients.
function orient_from_vectors(U::Vec2, V::Vec2, P::Vec2, Q::Vec2)
    dUP = det2(U, P); dVQ = det2(V, Q); dUQ = det2(U, Q); dVP = det2(V, P)
    return Harmonic(0.5 * (dUP + dVQ), 0.5 * (dUP - dVQ), 0.5 * (dUQ + dVP))
end

function dot_from_vectors(U::Vec2, V::Vec2, P::Vec2, Q::Vec2)
    a = _dotu(U, P); b = _dotu(V, Q); cc = _dotu(U, Q); dd = _dotu(V, P)
    return Harmonic(0.5 * (a + b), 0.5 * (a - b), 0.5 * (cc + dd))
end

"""det(Yb-Ya, Yp-Ya) as a harmonic."""
orient_harmonic(B::DeployBasis, a::Int, b::Int, p::Int) =
    orient_from_vectors(basis_c(B, b) - basis_c(B, a), basis_s(B, b) - basis_s(B, a),
                        basis_c(B, p) - basis_c(B, a), basis_s(B, p) - basis_s(B, a))
"""(Yb-Ya).(Yp-Ya) as a harmonic."""
dot_harmonic(B::DeployBasis, a::Int, b::Int, p::Int) =
    dot_from_vectors(basis_c(B, b) - basis_c(B, a), basis_s(B, b) - basis_s(B, a),
                     basis_c(B, p) - basis_c(B, a), basis_s(B, p) - basis_s(B, a))
"""|Yb-Ya|^2 as a harmonic."""
function len2_harmonic(B::DeployBasis, a::Int, b::Int)
    U = basis_c(B, b) - basis_c(B, a)
    V = basis_s(B, b) - basis_s(B, a)
    return dot_from_vectors(U, V, U, V)
end

"""
    harmonic_roots(h, lo, hi, tol = 0.0) -> Vector{Float64}

All roots of h in (lo, hi], ascending. Returns [] if |p| > amp (no real root)
or if h is identically zero to `tol` (a degenerate triple).

This uses the AMPLITUDE/PHASE form (atan2 + acos), not the tau = tan(theta/2)
quadratic, so it is well conditioned at theta = pi -- the singularity of the tau
chart (derivations/check.md D10). A contact exactly at the end of the range is
therefore decided by the same arithmetic as any other.
"""
function harmonic_roots(h::Harmonic, lo::Real, hi::Real, tol::Real = 0.0)
    out = Float64[]
    A = amp(h)
    (A <= tol && abs(h.p) <= tol) && return out  # identically zero
    A <= 0.0 && return out
    ratio = -h.p / A
    (ratio > 1.0 || ratio < -1.0) && return out
    phi = libm_atan2(h.r, h.q)
    psi = libm_acos(max(-1.0, min(1.0, ratio)))
    twopi = 2.0 * pi
    for base in (phi - psi, phi + psi)
        # shift into (lo, lo + 2pi]
        t = base
        t = fma(-twopi, floor((t - lo) / twopi), t)   # `t -= twopi * floor(...)`, contracted
        t <= lo && (t += twopi)
        while t <= hi + 1e-15
            push!(out, t)
            t += twopi
        end
    end
    sort!(out)
    return out
end

# ---------------------------------------------------------------------------
# The tau = 0 deflation (derivations/core.md T5.2b.2, three classes).
#
# Writing g(tau) = C + B tau + A tau^2 with tau = tan(theta/2) and
#     A = p - q ,   B = 2r ,   C = p + q  =  h(0) ,
# the order of vanishing of g at tau = 0 splits every candidate harmonic into
# three structural classes:
#
#   class 1   C != 0                g has no root at theta = 0
#   class 2   C = 0, B != 0         g = tau (A tau + B): a SIMPLE root at theta = 0
#   class 3   C = 0, B = 0, A != 0  g = A tau^2 = p(1 - cos theta): a DOUBLE root at
#                                   theta = 0 and, by Lemma T5.1e, NO other root on
#                                   (0, pi) -- constant sign, never a contact event
#   zero      A = B = C = 0         h identically zero: a permanent incidence
#                                   (T3.H.1/T3.H.2), struck from the candidate list
#
# C = 0 is not a rare degeneracy. It holds IDENTICALLY on the whole shape space for
# every permanent incidence and for every split-edge duplicate pair (T5.3's p + q = 0),
# so without deflation the theta = 0 root is reported as a root at ~1e-8 inside the
# interval and every pattern with split cuts reports Theta_max = 0. The Deriver
# measured 71 877 spurious mismatches without the deflation and 0 with it; the Checker
# found a further 143 that need the class-3 rule as well.
@enum HarmonicClass Zero Class1 Class2 Class3

"""Classification by |C| and |B| against rel_tol * scale(). Returns (class, A_zero) where
`A_zero` reports whether A = p - q also vanishes (which sends class 3 to Zero and class 2
to a theta = pi root)."""
function classify_harmonic_A(h::Harmonic, rel_tol::Real = 1e-11)
    sc = max(1e-300, scale(h))
    tol = rel_tol * sc
    A = h.p - h.q; B = 2.0 * h.r; C = h.p + h.q
    A_zero = abs(A) <= tol
    (abs(A) <= tol && abs(B) <= tol && abs(C) <= tol) && return (Zero, A_zero)
    abs(C) > tol && return (Class1, A_zero)
    abs(B) > tol && return (Class2, A_zero)
    # C = B = 0. A != 0 here (the all-zero case returned above), so this is class 3.
    return (Class3, A_zero)
end
classify_harmonic(h::Harmonic, rel_tol::Real = 1e-11) = classify_harmonic_A(h, rel_tol)[1]

"""
Roots of h in (lo, hi] with the tau = 0 root deflated away, ascending. Class 3 and
the identically-zero class return []. Class 2 returns the single remaining root
theta* = 2 atan(-B/A) when it lands in range (and theta = pi when A = 0, since
h = r sin theta then). Class 1 falls through to harmonic_roots.
"""
function harmonic_roots_deflated(h::Harmonic, lo::Real, hi::Real, rel_tol::Real = 1e-11)
    k, A_zero = classify_harmonic_A(h, rel_tol)
    if k == Zero
        return Float64[]  # permanent incidence: removed by the identity test, never a candidate
    elseif k == Class3
        # Lemma T5.1e: h = p(1 - cos theta) has constant sign on (0, pi). Never a contact.
        return Float64[]
    elseif k == Class2
        # g(tau) = tau (A tau + B); the tau = 0 root is the flat state and is deflated.
        A = h.p - h.q; B = 2.0 * h.r
        if A_zero
            # A = C = 0 => h = r sin(theta): the only remaining root on (0, pi] is pi.
            # M_PI as a double (an Irrational `pi` would compare exactly and lose the endpoint)
            PI = Float64(pi)
            return (PI > lo && PI <= hi + 1e-15) ? [PI] : Float64[]
        end
        tau = -B / A
        (!(tau > 0) || !isfinite(tau)) && return Float64[]
        th = 2.0 * libm_atan(tau)
        return (th > lo && th <= hi + 1e-15) ? [th] : Float64[]
    else
        return harmonic_roots(h, lo, hi)
    end
end

"""NOROOT for one candidate pair: no admissible deflated root in (0, eps).
Identically-zero harmonics are removed BEFORE the class split (T5.2b.0), so they
return true here; the caller is expected to have dropped them already."""
harmonic_no_root_in(h::Harmonic, eps::Real, rel_tol::Real = 1e-11) =
    isempty(harmonic_roots_deflated(h, 0.0, eps, rel_tol))

"""
    harmonic_fit(th, y) -> (h, max_res, scale)

Least-squares fit of samples (th_i, y_i) to p + q cos + r sin; returns the fit, the
max |residual| and the sample scale.
"""
function harmonic_fit(th::Vector{Float64}, y::Vector{Float64})
    n = length(th)
    M = Matrix{Float64}(undef, n, 3)
    sc = 0.0
    for i in 1:n
        M[i, 1] = 1.0
        M[i, 3], M[i, 2] = libm_sincos(th[i])
        sc = max(sc, abs(y[i]))
    end
    s = qr(M, ColumnNorm()) \ y   # column-pivoted Householder QR solve
    h = Harmonic(s[1], s[2], s[3])
    mr = 0.0
    for i in 1:n
        mr = max(mr, abs(harmonic_eval(h, th[i]) - y[i]))
    end
    return h, mr, sc
end
