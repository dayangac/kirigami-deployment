# method/contact.jl -- the exact first-contact calculus of ideas/ranking.md R2.
#
# Because Y(theta) = cos(theta/2) C + sin(theta/2) S (see deploy_basis.jl), a
# vertex-into-edge-interior contact between M'-vertex p and the M'-edge (a, b) is
# the simultaneous solution of one harmonic EQUATION and two harmonic INEQUALITIES:
#
#   det(Y_b - Y_a, Y_p - Y_a) = 0                                  (collinearity)
#   0 <= dot(Y_b - Y_a, Y_p - Y_a) <= |Y_b - Y_a|^2                (inside the segment)
#
# all three of the form p + q cos(theta) + r sin(theta), so the contact angles are
# closed-form arctangents and the interval tests are evaluations, not searches.
#
# The broad phase is the swept disc of a face: writing m_f(theta) for the moving
# centroid and rho_f for a bound on the distance from it to any vertex of the face,
# two faces can only touch while |m_f - m_g| <= rho_f + rho_g.
#
# Theorems T4/T5 implemented literally. Face pairs are `Tuple{Int,Int}` with 1-based
# faces; M'-vertex ids in the witnesses are 1-based with 0 = none (`first_root = -1.0`
# is the "no root" sentinel).

mutable struct SweptDiscs
    gc::Matrix{Float64}      # |F| x 2: m_f(theta) = cos(theta/2) gc + sin(theta/2) gs
    gs::Matrix{Float64}
    # rho_max is the EXACT swept radius about the moving centroid, and equals the flat
    # circumradius: in the face frame C_u - gc_f = x_u - xbar_f and
    # S_u - gs_f = -sigma_f J (x_u - xbar_f), so the columns are orthogonal and of equal
    # norm and the trajectory is a CIRCLE. rho = sqrt(|x|^2 + |chi|^2) = sqrt(2)*rho_max
    # is sound but loose by exactly sqrt(2); it is kept only to document that fact.
    rho::Vector{Float64}      # loose: max_u sqrt(|x_u|^2 + |chi_u|^2)  ( = sqrt(2) rho_max )
    rho_max::Vector{Float64}  # EXACT: max_u max(|x_u|, |chi_u|)  ( = circum )
    circum::Vector{Float64}   # circumradius of the face at the flat state
end

"""Per-face moving centroid (gc, gs) and the swept radii."""
function swept_discs(c::CutStructure, B::DeployBasis)
    F = n_faces(c.mesh)
    sd = SweptDiscs(zeros(F, 2), zeros(F, 2), zeros(F), zeros(F), zeros(F))
    for f in 1:F
        pf = c.prime_faces[f]
        gc = Vec2(0, 0); gs = Vec2(0, 0)
        for pv in pf
            gc += basis_c(B, pv)
            gs += basis_s(B, pv)
        end
        gc /= length(pf)
        gs /= length(pf)
        sd.gc[f, 1] = gc[1]; sd.gc[f, 2] = gc[2]
        sd.gs[f, 1] = gs[1]; sd.gs[f, 2] = gs[2]
        for pv in pf
            x = basis_c(B, pv) - gc
            chi = basis_s(B, pv) - gs
            # squaredNorm()/norm(): unfused reductions
            sd.rho[f] = max(sd.rho[f], sqrt(_sqnormu(x) + _sqnormu(chi)))
            sd.rho_max[f] = max(sd.rho_max[f], max(_normu(x), _normu(chi)))
            sd.circum[f] = max(sd.circum[f], _normu(x))
        end
    end
    return sd
end

"""min over theta in [0, theta_hi] of |m_f(theta) - m_g(theta)| (exact, closed form)."""
function min_center_distance(sd::SweptDiscs, f::Int, g::Int, theta_hi::Real)
    a = Vec2(sd.gc[f, 1] - sd.gc[g, 1], sd.gc[f, 2] - sd.gc[g, 2])
    b = Vec2(sd.gs[f, 1] - sd.gs[g, 1], sd.gs[f, 2] - sd.gs[g, 2])
    # |cos(t) a + sin(t) b|^2 on t in [0, theta_hi/2]; stationary at
    # tan(2t) = 2 a.b / (|a|^2 - |b|^2).
    aa = _sqnormu(a); bb = _sqnormu(b); ab = _dotu(a, b)
    hi = 0.5 * theta_hi
    # `aa*ct*ct + 2*ab*ct*st + bb*st*st` with (fadd (fmul x y) z) -> fma(x, y, z) applied
    # twice left to right; cos/sin pair -> __sincos_stret (docs/NUMERICS.md)
    function val(t)
        st, ct = libm_sincos(t)
        return sqrt(max(0.0, fma(bb * st, st, fma(aa * ct, ct, ((2 * ab) * ct) * st))))
    end
    best = min(val(0.0), val(hi))
    phi = 0.5 * atan(2 * ab, aa - bb)
    for k in -2:2
        t = phi + 0.5 * pi * k
        (t > 0 && t < hi) && (best = min(best, val(t)))
    end
    return best
end

"""
    candidate_pairs(c, sd, theta_hi, prune, use_static = false) -> Vector{Tuple{Int,Int}}

Face pairs surviving the broad phase. `prune == false` returns every unordered pair.
Both variants use the EXACT face-frame radius rho_max. `use_static == false` (the
default and the only sound choice) is T4.5a's MOVING test
min_theta |gamma_f(theta) - gamma_g(theta)| <= rho_f + rho_g, evaluated in closed form.
`use_static == true` is the flat-centroid test |c_f - c_g| <= rho_f + rho_g, which is
UNSOUND in principle: deployment contracts centroid distances while face radii stay
fixed, so it can discard a pair whose faces actually approach (derivations/check.md D8
measured 19 842 such pairs). It is retained only to reproduce that measurement.
"""
function candidate_pairs(c::CutStructure, sd::SweptDiscs, theta_hi::Real, prune::Bool,
                         use_static::Bool = false)
    F = n_faces(c.mesh)
    out = Tuple{Int,Int}[]
    if !prune
        sizehint!(out, F * (F - 1) ÷ 2)
        for i in 1:F, j in i+1:F
            push!(out, (i, j))
        end
        return out
    end
    # Uniform grid on the flat centroids, cell size = 2 * max radius, so only the
    # 3x3 neighbourhood can survive the radius test.
    # The EXACT swept radius in the face's own frame (derivations/check.md D1, core.md
    # T4.5b as corrected in round 2): relative to the moving centroid,
    #   C_u - gc_f = x_u - xbar_f  and  S_u - gs_f = -sigma_f J (x_u - xbar_f),
    # so the two columns are orthogonal and of EQUAL norm and the trajectory is a circle:
    #   max_theta |y_u(theta) - gamma_f(theta)| = max(|x|,|chi|) = |x_u - xbar_f| exactly.
    # sd.rho (= sqrt(|x|^2+|chi|^2) = sqrt(2) * rho_max) is sound but loose by exactly
    # sqrt(2); the sqrt(2) "correction" once prescribed for this frame is withdrawn.
    rr = sd.rho_max
    rmax = isempty(rr) ? 0.0 : maximum(rr)
    rmax <= 0 && return out
    cell = 2.0 * rmax
    key(f) = (Int(floor(sd.gc[f, 1] / cell)), Int(floor(sd.gc[f, 2] / cell)))
    grid = Dict{Tuple{Int,Int},Vector{Int}}()
    for f in 1:F
        push!(get!(grid, key(f), Int[]), f)
    end
    for f in 1:F
        k = key(f)
        for dx in -1:1, dy in -1:1
            cellv = get(grid, (k[1] + dx, k[2] + dy), nothing)
            cellv === nothing && continue
            for g in cellv
                g <= f && continue
                d = use_static ?
                    _normu(Vec2(sd.gc[f, 1] - sd.gc[g, 1], sd.gc[f, 2] - sd.gc[g, 2])) :
                    min_center_distance(sd, f, g, theta_hi)
                d <= rr[f] + rr[g] && push!(out, (f, g))
            end
        end
    end
    sort!(out)
    return out
end

mutable struct ContactWitness
    theta::Float64
    pv::Int; pa::Int; pb::Int  # M'-vertices: the vertex, the edge endpoints (0 = none)
    face_v::Int; face_e::Int
    corner_e::Int  # index of the edge within prime_faces[face_e]
    found::Bool
end
ContactWitness() = ContactWitness(0.0, 0, 0, 0, 0, 0, 0, false)

mutable struct ExactRangeReport
    first::ContactWitness
    # Candidates already in contact at theta = 0. Every pair of edge-adjacent faces is
    # one, and so is every pair of faces meeting only at a vertex of M, so this is a
    # diagnostic count, not a criterion.
    n_zero_contacts::Int
    theta_max::Float64   # = first.theta if found, else theta_hi
    n_pairs::Int         # face pairs actually tested
    n_candidates::Int    # (vertex, edge) candidates examined
    n_roots::Int         # harmonic roots that passed both interval tests
end

"""Smallest contact angle in (theta_lo, theta_hi] over the given face pairs."""
function exact_theta_max(c::CutStructure, B::DeployBasis, pairs::Vector{Tuple{Int,Int}},
                         theta_lo::Real = 1e-9, theta_hi::Real = Float64(pi),
                         check_zero_contacts::Bool = true)
    rep = ExactRangeReport(ContactWitness(), 0, Float64(theta_hi), 0, 0, 0)
    best = Float64(theta_hi)
    w = ContactWitness()
    PF = c.prime_faces

    function test(fe::Int, fv::Int)
        E = PF[fe]
        V = PF[fv]
        ne = length(E)
        for i in 1:ne
            a = E[i]; b = E[mod1(i + 1, ne)]
            U = basis_c(B, b) - basis_c(B, a); Vv = basis_s(B, b) - basis_s(B, a)
            L2 = dot_from_vectors(U, Vv, U, Vv)
            lscale = abs(L2.p) + amp(L2)
            lscale <= 0 && continue
            for p in V
                (p == a || p == b) && continue
                P = basis_c(B, p) - basis_c(B, a); Q = basis_s(B, p) - basis_s(B, a)
                det = orient_from_vectors(U, Vv, P, Q)
                rep.n_candidates += 1
                dscale = abs(det.p) + amp(det)
                dscale <= 0 && continue
                if check_zero_contacts
                    # Diagnostic only: how many candidates are already in contact at theta = 0.
                    geo0 = _normu(U) * _normu(P)
                    if abs(harmonic_eval(det, 0.0)) <= 1e-9 * geo0
                        s0 = _dotu(U, P); l0 = _sqnormu(U)
                        t0 = 1e-12 * max(1e-300, l0)
                        (s0 >= -t0 && s0 <= l0 + t0) && (rep.n_zero_contacts += 1)
                    end
                end
                for th in harmonic_roots_deflated(det, theta_lo, best)
                    D = dot_from_vectors(U, Vv, P, Q)
                    s = harmonic_eval(D, th); l2 = harmonic_eval(L2, th)
                    tol = 1e-12 * lscale
                    (s < -tol || s > l2 + tol) && continue
                    rep.n_roots += 1
                    if th < best
                        best = th
                        w = ContactWitness(th, p, a, b, fv, fe, i, true)
                    end
                    break  # roots are ascending
                end
            end
        end
    end

    for (f, g) in pairs
        rep.n_pairs += 1
        test(f, g)
        test(g, f)
    end
    rep.first = w
    rep.theta_max = w.found ? w.theta : Float64(theta_hi)
    return rep
end

# ---------------------------------------------------------------------------
# T4.2" (derivations/core.md): the EXACT deployment range.
#
# exact_theta_max() above returns theta_1, the first CONTACT angle. That is not
# Theta_max: a contact can be a graze (tangency without crossing). Measured
# counterexample: hexagons_auto has theta_1 = 1.047198 (two split-edge duplicates
# become collinear end-to-end, vertex lands exactly on an edge ENDPOINT, the faces
# touch at a point and separate) while the first interior overlap is at 2.094395.
#
# T4.2" : let C(X) = {theta_1 < ... < theta_M} be the complete contact-angle set and
# theta_0 := 0, theta_{M+1} := theta_hi. Then
#     Theta_max = theta_{i*},  i* = min{ i >= 0 : interiors overlap at the midpoint
#                                        of (theta_i, theta_{i+1}) },
# and Theta_max = theta_hi if no such i exists. i* = 0 is the zero-range case: the
# structure penetrates immediately (split-cut duplicates moving into each other),
# which no root scan can see because that contact sits at theta = 0.
#
# Completeness of C(X) is Lemma T4.2 / Corollary T4.2': two closed simple polygons
# with disjoint interiors that intersect must have a vertex of one on the boundary
# of the other, so every overlap-status change is a (vertex, edge) contact.
mutable struct OverlapRangeReport
    theta_max::Float64
    candidates::Vector{Float64}  # C(X), ascending, deduplicated
    i_star::Int                  # 0-based index into `candidates`: 0 => zero range, -1 => none
    n_intervals_tested::Int
    first_contact::ContactWitness  # the theta_1 that min-over-roots would have returned
    zero_range::Bool
    n_overlap_tests::Int
end

"""
The complete contact-angle set C(X) over the given face pairs: harmonic roots of the
orientation determinant that also satisfy both interval (projection) inequalities.
"""
function contact_angles(c::CutStructure, B::DeployBasis, pairs::Vector{Tuple{Int,Int}},
                        theta_lo::Real = 1e-9, theta_hi::Real = Float64(pi), dedup::Real = 1e-9)
    PF = c.prime_faces
    out = Float64[]
    function test(fe::Int, fv::Int)
        E = PF[fe]
        V = PF[fv]
        ne = length(E)
        for i in 1:ne
            a = E[i]; b = E[mod1(i + 1, ne)]
            U = basis_c(B, b) - basis_c(B, a); Vv = basis_s(B, b) - basis_s(B, a)
            L2 = dot_from_vectors(U, Vv, U, Vv)
            lscale = abs(L2.p) + amp(L2)
            lscale <= 0 && continue
            for p in V
                (p == a || p == b) && continue
                P = basis_c(B, p) - basis_c(B, a); Q = basis_s(B, p) - basis_s(B, a)
                det = orient_from_vectors(U, Vv, P, Q)
                abs(det.p) + amp(det) <= 0 && continue
                D = dot_from_vectors(U, Vv, P, Q)
                # Every root, not just the first: T4.2" needs the complete set, and a graze
                # early on must not hide a genuine crossing later.
                for th in harmonic_roots_deflated(det, theta_lo, theta_hi)
                    s = harmonic_eval(D, th); l2 = harmonic_eval(L2, th)
                    tol = 1e-12 * lscale
                    (s < -tol || s > l2 + tol) && continue
                    push!(out, th)
                end
            end
        end
    end
    for (f, g) in pairs
        test(f, g)
        test(g, f)
    end
    sort!(out)
    ded = Float64[]
    for v in out
        (isempty(ded) || v - ded[end] > dedup) && push!(ded, v)
    end
    return ded
end

"""
Theta_max by T4.2". `shrink` is passed to polygons_overlap; the overlap probe is
restricted to `pairs`, which is sound whenever `pairs` came from the (sound) broad
phase, since a pair outside it can never touch.
"""
function exact_theta_max_overlap(c::CutStructure, B::DeployBasis, pairs::Vector{Tuple{Int,Int}},
                                 theta_lo::Real = 1e-9, theta_hi::Real = Float64(pi),
                                 shrink::Real = 1e-9)
    rep = OverlapRangeReport(0.0, contact_angles(c, B, pairs, theta_lo, theta_hi), -1, 0,
                             ContactWitness(), false, 0)
    # theta_1, for reporting what the (wrong) min-over-roots rule would have said.
    rep.first_contact = exact_theta_max(c, B, pairs, theta_lo, theta_hi, false).first

    PF = c.prime_faces
    function overlaps_at(th::Float64)
        rep.n_overlap_tests += 1
        Y = basis_eval(B, th)
        for (f, g) in pairs
            poly_a = Vec2[Y[pv] for pv in PF[f]]
            poly_b = Vec2[Y[pv] for pv in PF[g]]
            polygons_overlap(poly_a, poly_b, Float64(shrink)) && return true
        end
        return false
    end

    M = length(rep.candidates)
    lo = 0.0
    for i in 0:M
        hi = i < M ? rep.candidates[i + 1] : Float64(theta_hi)
        if hi <= lo  # duplicate / degenerate slab, nothing to probe
            lo = hi
            continue
        end
        rep.n_intervals_tested += 1
        if overlaps_at(0.5 * (lo + hi))
            rep.i_star = i
            rep.theta_max = lo
            rep.zero_range = (i == 0)
            return rep
        end
        lo = hi
    end
    rep.i_star = -1
    rep.theta_max = Float64(theta_hi)
    return rep
end

# ---------------------------------------------------------------------------
# The validity certificate (derivations/core.md T5.2b.0 / T5.2b', round 3;
# derivations/check.md R3.3, measured 0 violations on 175 samples).
#
#     VALID(eps)  :=  POS  /\  NOOVERLAP(theta_1)  /\  NOROOT(eps)
#
#   POS               every face of M has signed area > 0 at theta = 0.
#   NOOVERLAP(theta_1) the deployed faces have pairwise disjoint interiors at the SINGLE
#                     angle theta_1 = eps/2 -- one exact polygon-polygon test per pair.
#   NOROOT(eps)       for every candidate pair in the C-list, the orientation harmonic
#                     has no admissible DEFLATED root in (0, eps). ADMISSIBLE means the
#                     root also passes the two interval (projection) inequalities
#                     0 <= <w-a, b-a> <= |b-a|^2 of T4.1b -- the vertex is on the edge
#                     SEGMENT, not merely on its infinite line. Without them NOROOT is
#                     not the C-list of Corollary T4.2' but a strictly larger root set,
#                     and it rejected 87.3% of designs with exact Theta_max >= 1 rad
#                     (F32, results/kill/jitter/cert_diagnosis.md).
#
# Together these imply Theta_max >= eps (Proposition T5.2b'): NOROOT makes the overlap
# status constant on (0, eps) by a connectedness argument (the status can change only at
# a candidate root), and the single test at theta_1 fixes which constant it is.
#
# The theta = 0 overlap test is NOT part of this certificate and is withdrawn: at
# theta = 0 the two copies of every split edge COINCIDE and adjacent faces touch along
# whole shared edges, so the test is degenerate rather than merely different. Numbers
# computed with a theta = 0 test are labelled as such wherever they are reported.
#
# The C-list is the ordered candidate list MINUS the pairs whose harmonic is identically
# zero (the permanent incidences of T3.H.1/T3.H.2: hinge-adjacent copies, and the
# coincident copies of a shared vertex). Those are removed by an identity test on the
# coefficients, before the three-class split, never by a root test -- without the removal
# NOROOT is unsatisfiable on any pattern with a hinge.
mutable struct ValidityCertificate
    pos::Bool
    nooverlap::Bool
    noroot::Bool

    eps::Float64; theta_1::Float64
    min_signed_area::Float64
    n_inverted::Int

    n_pairs::Int             # face pairs tested (broad phase)
    n_candidates::Int        # ordered (vertex, edge) candidates seen
    n_identically_zero::Int  # struck by the identity test before classification
    n_class1::Int; n_class2::Int; n_class3::Int
    n_roots_deflated::Int    # admissible roots found in (0, eps)
    # How many roots the UNdeflated test would have reported in (0, eps): the size of the
    # tau = 0 artefact this deflation removes.
    n_roots_undeflated::Int
    # Roots of the orientation harmonic in (0, eps) that are NOT contacts because the
    # vertex lies off the edge SEGMENT (the two interval inequalities of T4.1b). These are
    # not violations of NOROOT; they are counted only to report the size of the artefact
    # that dropping the interval tests used to produce (F32).
    n_roots_inadmissible::Int

    first_root::Float64  # smallest admissible deflated root, if any (-1 = none)
    bad_pv::Int; bad_a::Int; bad_b::Int  # 0 = none
end
ValidityCertificate() = ValidityCertificate(false, false, false, 0.0, 0.0, 0.0, 0,
                                            0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0)
valid(ct::ValidityCertificate) = ct.pos && ct.nooverlap && ct.noroot

"""
    validity_certificate(c, B, X, pairs, eps, shrink = 1e-12) -> ValidityCertificate

`pairs` should come from the sound (moving) broad phase; a pair outside it cannot touch
anywhere in (0, pi], so restricting NOROOT to it is sound.
"""
function validity_certificate(c::CutStructure, B::DeployBasis, X::Vector{Vec2},
                              pairs::Vector{Tuple{Int,Int}}, eps::Real, shrink::Real = 1e-12)
    cert = ValidityCertificate()
    cert.eps = eps
    cert.theta_1 = 0.5 * eps
    m = c.mesh
    PF = c.prime_faces

    # --- POS: every face of M positively oriented at theta = 0 (T5.1) -------------
    # (the same shoelace sum as face_signed_area, evaluated on X directly)
    let lo = 0.0, hi = 0.0, first = true
        for f in 1:n_faces(m)
            vs = m.faces[f]
            n = length(vs)
            s = 0.0
            for i in 1:n
                p = X[vs[i]]; q = X[vs[mod1(i + 1, n)]]
                s += fma(p[1], q[2], -(q[1] * p[2]))   # as face_signed_area: inline a*b - c*d, contracted
            end
            a = 0.5 * s
            if first
                lo = hi = a
                first = false
            end
            lo = min(lo, a)
            hi = max(hi, a)
            a <= 0 && (cert.n_inverted += 1)
        end
        cert.min_signed_area = lo
        cert.pos = (cert.n_inverted == 0)
    end

    # --- NOOVERLAP(theta_1): one exact polygon test per pair, at theta_1 = eps/2 ----
    let Y = basis_eval(B, cert.theta_1), ok = true
        for (f, g) in pairs
            A = Vec2[Y[pv] for pv in PF[f]]
            Bp = Vec2[Y[pv] for pv in PF[g]]
            if polygons_overlap(A, Bp, Float64(shrink))
                ok = false
                break
            end
        end
        cert.nooverlap = ok
    end

    # --- NOROOT(eps): no admissible deflated root in (0, eps) ----------------------
    noroot = true
    function scan(fe::Int, fv::Int)
        E = PF[fe]
        V = PF[fv]
        ne = length(E)
        for i in 1:ne
            a = E[i]; b = E[mod1(i + 1, ne)]
            U = basis_c(B, b) - basis_c(B, a); Vv = basis_s(B, b) - basis_s(B, a)
            # The edge's squared length as a harmonic, for the two interval (projection)
            # inequalities below. Identical to the exact scan's, so the two agree root by root.
            L2 = dot_from_vectors(U, Vv, U, Vv)
            lscale = abs(L2.p) + amp(L2)
            lscale <= 0 && continue
            for p in V
                (p == a || p == b) && continue
                P = basis_c(B, p) - basis_c(B, a); Q = basis_s(B, p) - basis_s(B, a)
                det = orient_from_vectors(U, Vv, P, Q)
                cert.n_candidates += 1
                # The identity test of T3.H.1, on the COEFFICIENTS, before any classification:
                # these are the permanent incidences and are struck from the C-list entirely.
                geo = _normu(U) * _normu(P)
                sc = max(scale(det), geo)
                if abs(det.p) + abs(det.q) + abs(det.r) <= 1e-11 * max(1e-300, sc)
                    cert.n_identically_zero += 1
                    continue
                end
                k = classify_harmonic(det)
                if k == Class1
                    cert.n_class1 += 1
                elseif k == Class2
                    cert.n_class2 += 1
                elseif k == Class3
                    cert.n_class3 += 1
                end
                # The size of the artefact the deflation removes, for reporting.
                isempty(harmonic_roots(det, 0.0, eps)) || (cert.n_roots_undeflated += 1)
                # ADMISSIBILITY (F32). A root of the orientation harmonic only says the vertex
                # is collinear with the INFINITE LINE through the edge. A contact -- and hence a
                # transition of the overlap status, which is all Proposition T5.2b' needs -- also
                # requires the vertex to lie on the SEGMENT. That is the pair of interval
                # inequalities 0 <= <w-a, b-a> <= |b-a|^2 of T4.1b, which the C-list of Corollary
                # T4.2' carries and which contact_angles() applies. Dropping them (round 3's
                # T5.2b.0 says so explicitly) made NOROOT reject 87.3% of designs whose exact
                # Theta_max was 1-2.4 rad: on the five diagnosed rows all 168 reported roots had
                # the vertex outside the segment (|s|/|e|^2 up to 6.8), none inside, none a graze.
                # The test is inclusive at the tolerance, so a borderline root stays admissible
                # and the certificate stays conservative.
                D = dot_from_vectors(U, Vv, P, Q)
                rts = Float64[]
                for th in harmonic_roots_deflated(det, 0.0, eps)
                    s = harmonic_eval(D, th); l2 = harmonic_eval(L2, th)
                    tol = 1e-12 * lscale
                    if s < -tol || s > l2 + tol
                        cert.n_roots_inadmissible += 1
                        continue
                    end
                    push!(rts, th)
                end
                if !isempty(rts)
                    cert.n_roots_deflated += 1
                    # first_root must be the SMALLEST admissible root over every candidate, not
                    # the first one the pair loop happens to meet: it is the supremum of the eps
                    # for which NOROOT(eps) holds, and callers use it as exactly that.
                    if cert.first_root < 0 || rts[1] < cert.first_root
                        cert.first_root = rts[1]
                        cert.bad_pv = p
                        cert.bad_a = a
                        cert.bad_b = b
                    end
                    noroot = false
                end
            end
        end
    end
    for (f, g) in pairs
        cert.n_pairs += 1
        scan(f, g)
        scan(g, f)
    end
    cert.noroot = noroot
    return cert
end

"""Convenience: builds the basis and the sound broad phase itself."""
function validity_certificate(c::CutStructure, X::Vector{Vec2}, eps::Real, shrink::Real = 1e-12)
    B = deploy_basis(c, X)
    sd = swept_discs(c, B)
    return validity_certificate(c, B, X, candidate_pairs(c, sd, Float64(pi), true), eps, shrink)
end

"""
The root enumeration above finds TRANSVERSAL contacts only. A flat state in which
two faces are coincident-adjacent and move into each other -- the failure mode Eq. (9)
is written to prevent -- has its contact at theta = 0 as a tangency, invisible to the
root scan; the deployment range is then 0. This decides that case directly.
"""
penetrates_immediately(c::CutStructure, X::Vector{Vec2}, probe::Real = 1e-6, shrink::Real = 1e-12) =
    has_collision(c, deploy(c, X, probe).Y, Float64(shrink))
