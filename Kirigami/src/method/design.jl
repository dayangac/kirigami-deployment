# method/design.jl -- the deliverable method, as one API (specs/builder_method.md D7 update).
#
# Port of code/src/method/design.{hpp,cpp}. The project has two halves and this file is
# the entry point to both.
#
# (A) CHARACTERIZATION -- `characterize(mesh, sigma, X)`. Given ANY oriented planar graph
#     and ANY flat embedding of it, return the EXACT deployment range Theta_max by T4.2''
#     (the complete contact-angle set C(X), then the first interval whose midpoint has an
#     interior overlap), the contact list itself, and the repaired three-part validity
#     certificate POS /\ NOOVERLAP(eps/2) /\ NOROOT(eps) of T5.2b'. Nothing here is new
#     numerics: it is the thin, named surface over method/contact.jl.
#
# (B) CONSTRUCTION -- `design_constrained(mesh, sigma, X_ini, DesignOptions)`. The K9
#     variant-(b) embedding: inside the Eq. (4) null space X(t) = X0 + Phi t, solve for a
#     point where every face corner is strictly convex (cross_i >= delta) AND every split
#     cut opens outward (q_e >= delta'), then characterize it (STATE.md F36: 36/400, i.e.
#     9.0 %, short of its own 10 % pre-registered bar).
#
# (C) BASELINE -- `design_baseline(mesh, sigma, X_ini, ...)` is the 2026 Eq. (6)
#     projection and nothing else, i.e. t = 0, characterized by exactly the same code.
#
# (D) RANGE-MAXIMISING CONSTRUCTION -- `design_range_max` runs EXACTLY the per-design
#     pipeline of apps/kill_k9c.cpp (arms k9, k9b, stage A from each start, stage B on
#     the stage-A winner, best-of-three by exact Theta_max), with no wall-clock gating.
#
# DEFAULTS ARE K9'S / K9c's, so that with the runs' seeds these functions reproduce
# results/kill/k9/k9.csv and results/kill/k9c/k9c.csv row for row (test/test_design.jl).
#
# Shape-space coefficients `t` are plain `Vector{Float64}` of length 2 * dim_null; an
# empty vector stands for Eigen's empty VectorXd (t = 0 / no coefficients).

# ---------------------------------------------------------------------------
# (A) Characterization.

Base.@kwdef mutable struct CharacterizeOptions
    # The angle the certificate is asked for. K9/K9b use 0.3 rad.
    eps::Float64 = 0.3
    theta_hi::Float64 = pi
    theta_lo::Float64 = 1e-9
    # Shrink passed to polygons_overlap inside the T4.2'' interval probe. K9 uses 1e-9;
    # the historical value in K2a/K5/K6 was 1e-12. F34 removed the shrink from the exact
    # open-interior test, so this only affects the probe's degenerate-touch tolerance.
    overlap_shrink::Float64 = 1e-9
    # Shrink of the certificate's single NOOVERLAP(eps/2) polygon test.
    cert_shrink::Float64 = 1e-12
    # Run the independent bisection cross-check (collision.jl has_collision on a uniform
    # ladder, then bisection). It is a CROSS-CHECK, not ground truth (F34). Off by default
    # because it costs a few thousand overlap tests.
    referee::Bool = false
    referee_grid::Int = 4000
    referee_iters::Int = 50
    referee_shrink::Float64 = 1e-12
end

mutable struct Characterization
    # --- the range -----------------------------------------------------------
    theta_max::Float64            # EXACT, by T4.2''
    zero_range::Bool              # penetrates at 0+ (i* == first candidate)
    contacts::Vector{Float64}     # C(X): every contact angle, ascending, deduplicated
    first_contact::ContactWitness # theta_1, i.e. what a min-over-roots would return
    i_star::Int
    # Why the range is zero, when it is: "split-inward", "inverted", "vertex-edge",
    # "hinge-wedge"; "n/a" when theta_max > 0. K6's classification.
    binding::String

    # --- the certificate -----------------------------------------------------
    certificate::ValidityCertificate  # evaluated at opt.eps
    certified::Bool                   # valid(certificate)
    # The LARGEST eps this X certifies: the first admissible deflated root (minus the
    # root routine's tolerance), or theta_hi when the harmonic has no admissible root.
    eps_max::Float64

    # --- diagnostics ---------------------------------------------------------
    min_signed_area::Float64
    n_inverted::Int
    n_pairs::Int
    theta_bisect::Float64   # < 0 when opt.referee == false
    n_split::Int
    n_hinge::Int
    n_faces::Int
    n_vertices::Int
    med_edge::Float64
end
Characterization() = Characterization(0.0, false, Float64[], ContactWitness(), -1, "n/a",
                                      ValidityCertificate(), false, 0.0,
                                      0.0, 0, 0, -1.0, 0, 0, 0, 0, 0.0)

# ---------------------------------------------------------------------------
# (B)/(C) Construction.

Base.@kwdef mutable struct DesignOptions
    # --- the constrained embedding (K9 variant (b)) ---------------------------
    delta_convex::Float64 = 1e-3  # corner margin, in med_edge^2
    delta_split::Float64 = 1e-3   # split-cut margin q_e, in med_edge^2
    restarts::Int = 3             # Gaussian restarts IN ADDITION to t = 0
    barrier_stages::Int = 6       # proximity (phase B) stages; 0 disables phase B
    max_iter::Int = 600           # L-BFGS iterations per continuation stage
    start_scale::Float64 = 0.1    # restart scale, in median edge lengths
    seed::UInt32 = 9000
    # The remaining convex_embed knobs, at their K9 values.
    solve_factor::Float64 = 2.0   # the penalty phase aims at solve_factor * delta
    lambda_rel::Float64 = 1e-6    # |t|^2 regulariser of phase A, relative to med_edge^2
    barrier_w0::Float64 = 1.0
    barrier_factor::Float64 = 0.25
    barrier_iter::Int = 300
    # K9 hands variant (b) two extra warm starts when the cold solve is infeasible: the
    # convexity-only minimiser (variant (a)) and the split-only repair of K6.
    warm_starts::Bool = true

    # --- optional stage 2 (K9b lever (iv)) ------------------------------------
    # After feasibility, maximise the range with range_opt.jl and KEEP the result only if
    # it is still exactly feasible and its exact Theta_max strictly improved.
    maximise_eps::Bool = false
    range_opt::RangeOptOptions = RangeOptOptions(rounds = 8, iters_per_round = 10,
                                                barrier_eps = 0.2)
    # step_cap of range_opt, in median edge lengths (0 = uncapped). K9b uses 0.25.
    range_step_cap::Float64 = 0.25

    # --- what to measure at the end -------------------------------------------
    characterize::CharacterizeOptions = CharacterizeOptions()
end

mutable struct DesignResult
    ok::Bool                  # the Eq. (4) system solved and dim_null >= 1
    status::String            # "ok", or why not
    method::String            # "constrained", "baseline" or "range_max"

    X::Vector{Vec2}           # the design
    X0::Vector{Vec2}          # the Eq. (6) projection (t = 0)
    t::Vector{Float64}        # shape-space coefficients of X, length 2 * dim_null
    Phi::Matrix{Float64}      # N x dim_null
    dim_null::Int

    # Exact constraint values at X, in med_edge^2. Meaningless for the baseline, where
    # they are still reported (they are what makes the baseline infeasible).
    feasible::Bool            # every cross >= delta AND every q >= delta', EXACT
    min_cross::Float64
    min_q::Float64
    min_mu::Float64
    n_bad::Int
    n_bad_q::Int
    n_bad_mu::Int
    n_corners::Int
    n_nonconvex_x0::Int
    n_inverted::Int

    dist_ini::Float64         # ||X - X_ini||_F / (sqrt(N) med_edge), per-vertex RMS
    dist_ini_x0::Float64      # the same at X0
    t_norm_rel::Float64
    barrier_stages_kept::Int
    used_range_opt::Bool      # stage 2 accepted
    med_edge::Float64

    ch::Characterization      # at X
end
DesignResult() = DesignResult(false, "", "", Vec2[], Vec2[], Float64[], zeros(0, 0), 0,
                              false, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0,
                              0.0, 0.0, 0.0, 0, false, 0.0, Characterization())

# ---------------------------------------------------------------------------
# (D) The range-maximising construction (K9c).

Base.@kwdef mutable struct RangeMaxOptions
    delta_convex::Float64 = 1e-3      # the k9 and t = 0 arms' margin, in med^2
    delta_split::Float64 = 1e-3
    delta_wide::Float64 = 1e-2        # the k9b arm's margin, in med^2
    kappa::Float64 = 5e-3             # range_embed softmin temperature, in med^2
    stages::Int = 6                   # range_embed continuation stages  (K9c run: 6)
    iter_per_stage::Int = 120         # L-BFGS iterations per stage      (K9c run: 120)
    max_iter::Int = 600               # L-BFGS iterations of the proximity arms
    seed::UInt32 = 9300               # base seed; the arms use +0, +101, +211
    # Which starts stage A is given. Turning an arm off also removes it as a candidate.
    arm_proximity::Bool = true        # the K9 point
    arm_proximity_wide::Bool = true   # the K9b point
    start_x0::Bool = true             # t = 0
    stage_b::Bool = true              # maximize_margin_range on the stage-A winner
    margin_range::MarginRangeOptions = MarginRangeOptions()
    characterize::CharacterizeOptions = CharacterizeOptions()
end

# One candidate point, as it was scored. `tag` is K9c's provenance string: "k9", "k9b",
# "k9c/k9", "k9c/k9b", "k9c/x0", and the same with "+B" once stage B improved it.
mutable struct RangeMaxArm
    tag::String
    feasible::Bool      # EXACT, at this arm's own delta
    margin::Float64     # EXACT 0+ margin m(X), in med^2
    theta_max::Float64  # EXACT, T4.2''
    eps_max::Float64    # certified
    stage_a_winner::Bool
end

mutable struct RangeMaxResult
    design::DesignResult        # at the winning point; method = "range_max"
    provenance::String          # the winning arm's tag, i.e. K9c's `best_src`
    stage_b_used::Bool          # the winner came through stage B
    margin::Float64             # EXACT m(X) at the winner, in med^2
    arms::Vector{RangeMaxArm}   # every candidate scored, in the order it was scored
end
RangeMaxResult() = RangeMaxResult(DesignResult(), "", false, 0.0, RangeMaxArm[])

# ---------------------------------------------------------------------------
# Internals.

# The scale every area-valued quantity here is measured in. Copied verbatim from
# apps/kill_common.hpp (nth_element at L.size()/2, i.e. the upper median).
function median_edge_length(m::Mesh)
    L = [norm(m.X[e.key.a] - m.X[e.key.b]) for e in m.edges]
    isempty(L) && return 1.0
    return sort(L)[length(L) ÷ 2 + 1]
end

"""Number of faces of `m` (under positions X) with signed area <= rel_tol * bbox area."""
function count_inverted_faces(m::Mesh, X::Vector{Vec2}, rel_tol::Float64 = 1e-12)
    t = Mesh(X, m.faces)   # face_signed_area only reads X and faces
    lo = X[1]
    hi = X[1]
    for p in X
        lo = min.(lo, p)
        hi = max.(hi, p)
    end
    bbox = max(1e-300, (hi[1] - lo[1]) * (hi[2] - lo[2]))
    n = 0
    for f in 1:n_faces(m)
        face_signed_area(t, f) <= rel_tol * bbox && (n += 1)
    end
    return n
end

# The bisection cross-check of K2a/K5/K6/K9: uniform ladder to the first colliding
# angle, then bisection. NOT ground truth (STATE.md F34); reported as a cross-check.
function referee_theta(c::CutStructure, X::Vector{Vec2}, grid::Int, iters::Int, shrink::Float64)
    col(th) = has_collision(c, deploy(c, X, th).Y, shrink)
    col(1e-7) && return 0.0
    lo = 0.0
    hi = -1.0
    for i in 1:grid
        th = pi * i / grid
        if col(th)
            hi = th
            break
        end
        lo = th
    end
    hi < 0 && return Float64(pi)
    for _ in 1:iters
        mid = 0.5 * (lo + hi)
        if col(mid)
            hi = mid
        else
            lo = mid
        end
    end
    return lo
end

# Mesh + sigma + X, validated once, so every entry point rejects the same bad inputs
# (C++ std::invalid_argument -> ArgumentError).
function prepare(mesh::Mesh, sigma::Vector{Int}, X::Vector{Vec2}, who::String)
    n_vertices(mesh) == 0 && throw(ArgumentError("$who: mesh has no vertices"))
    n_faces(mesh) == 0 && throw(ArgumentError("$who: mesh has no faces"))
    length(sigma) != n_faces(mesh) &&
        throw(ArgumentError("$who: sigma has $(length(sigma)) entries for $(n_faces(mesh)) faces"))
    for s in sigma
        (s != 1 && s != -1) && throw(ArgumentError("$who: sigma entries must be +-1"))
    end
    length(X) != n_vertices(mesh) &&
        throw(ArgumentError("$who: X has $(length(X)) points for $(n_vertices(mesh)) vertices"))
    m = Mesh(copy(X), mesh.faces)
    m.sigma = copy(sigma)
    m.periodic = mesh.periodic
    build_topology!(m)
    median_edge_length(m) > 0 ||
        throw(ArgumentError("$who: degenerate embedding, median edge length is 0"))
    return m
end

# K6's classification of what binds at a design whose range is zero.
function binding_type(c::CutStructure, m::Mesh, X::Vector{Vec2}, zero_range::Bool)
    for v in zero_plus_q(c, X)
        v <= 0 && return "split-inward"
    end
    t = Mesh(X, m.faces)
    for f in 1:n_faces(m)
        face_signed_area(t, f) <= 0 && return "inverted"
    end
    return zero_range ? "vertex-edge" : "hinge-wedge"
end

# The exact constraint values of variant (b) at an arbitrary X, in med^2.
mutable struct ExactMargins
    min_cross::Float64
    min_q::Float64
    min_mu::Float64
    n_bad::Int
    n_bad_q::Int
    n_bad_mu::Int
    n_corners::Int
    n_nonconvex::Int
    feasible::Bool
end

function exact_margins(m::Mesh, c::CutStructure, X::Vector{Vec2}, delta::Float64, delta_q::Float64)
    e = ExactMargins(0.0, 0.0, 0.0, 0, 0, 0, 0, 0, false)
    cr = corner_crosses(m, X)
    e.n_corners = length(cr)
    e.min_cross = Inf
    for v in cr
        e.min_cross = min(e.min_cross, v)
        v < delta && (e.n_bad += 1)
        v <= 0 && (e.n_nonconvex += 1)
    end
    q = zero_plus_q(c, X)
    e.min_q = Inf
    for v in q
        e.min_q = min(e.min_q, v)
        v < delta_q && (e.n_bad_q += 1)
    end
    mu = zero_plus_corner_margin(c, X)
    e.min_mu = Inf
    for v in mu
        e.min_mu = min(e.min_mu, v)
        v <= 0 && (e.n_bad_mu += 1)
    end
    isempty(cr) && (e.min_cross = 0.0)
    isempty(q) && (e.min_q = 0.0)
    isempty(mu) && (e.min_mu = NaN)
    e.feasible = (e.n_bad == 0 && e.n_bad_q == 0)
    return e
end

function rms_move(X::Vector{Vec2}, Y::Vector{Vec2}, med::Float64)
    s = 0.0
    for i in eachindex(X)
        d = X[i] - Y[i]
        s += dot(d, d)
    end
    return sqrt(s / max(1, length(X))) / med
end

function embed_options(o::DesignOptions, delta_split::Float64)
    return ConvexEmbedOptions(delta_rel = o.delta_convex, split_delta_rel = delta_split,
                              solve_factor = o.solve_factor, n_random = o.restarts,
                              start_scale = o.start_scale, seed = o.seed, max_iter = o.max_iter,
                              lambda_rel = o.lambda_rel, barrier_stages = o.barrier_stages,
                              barrier_w0 = o.barrier_w0, barrier_factor = o.barrier_factor,
                              barrier_iter = o.barrier_iter)
end

# The Eq. (4)/(6) shape space of one (mesh, sigma).
mutable struct ShapeSpace
    X0::Vector{Vec2}
    Phi::Matrix{Float64}
    k::Int
    ok::Bool
    status::String
end

function shape_space_of(m::Mesh, c::CutStructure)
    s = ShapeSpace(Vec2[], zeros(0, 0), 0, false, "")
    hs = holes_partition(c)
    sys = assemble_system(c, hs, m.X, Fixed)
    sr = solve_system(sys, m.X)
    if !sr.projection_ok
        s.status = "Eq. (6) projection failed"
        return s
    end
    s.X0 = matrix_to_points(sr.X0)
    s.Phi = sr.Phi
    s.k = sr.dim_null
    s.ok = true
    return s
end

# ---------------------------------------------------------------------------

"""
    characterize(mesh, sigma, X, opt = CharacterizeOptions()) -> Characterization

The exact deployment range Theta_max by T4.2'', the contact list C(X), and the repaired
validity certificate POS /\\ NOOVERLAP(eps/2) /\\ NOROOT(eps) of T5.2b' at `opt.eps`, plus
the largest certified eps. Throws ArgumentError on an empty mesh, a sigma of the wrong
length, an X of the wrong length, or a mesh with no faces.
"""
function characterize(mesh::Mesh, sigma::Vector{Int}, X::Vector{Vec2},
                      opt::CharacterizeOptions = CharacterizeOptions())
    m = prepare(mesh, sigma, X, "characterize")
    c = make_cut(m)
    med = median_edge_length(m)

    r = Characterization()
    r.n_faces = n_faces(m)
    r.n_vertices = n_vertices(m)
    r.n_split = n_split(c)
    r.n_hinge = n_hinge(c)
    r.med_edge = med

    B = deploy_basis(c, X)
    sd = swept_discs(c, B)
    pairs = candidate_pairs(c, sd, opt.theta_hi, true)
    r.n_pairs = length(pairs)

    # --- the exact range, T4.2'' ---------------------------------------------
    orr = exact_theta_max_overlap(c, B, pairs, opt.theta_lo, opt.theta_hi, opt.overlap_shrink)
    r.theta_max = orr.theta_max
    r.zero_range = orr.zero_range
    r.contacts = orr.candidates
    r.first_contact = orr.first_contact
    r.i_star = orr.i_star
    if r.theta_max <= opt.theta_lo
        r.binding = binding_type(c, m, X, orr.zero_range)
    end

    # --- the certificate at eps, and the largest certified eps ----------------
    r.certificate = validity_certificate(c, B, X, pairs, opt.eps, opt.cert_shrink)
    r.certified = valid(r.certificate)
    r.min_signed_area = r.certificate.min_signed_area
    r.n_inverted = r.certificate.n_inverted

    c_pi = validity_certificate(c, B, X, pairs, opt.theta_hi, opt.cert_shrink)
    kRootTol = 1e-12
    eps_star = c_pi.first_root > kRootTol ? c_pi.first_root - kRootTol : 0.0
    if c_pi.first_root <= 0
        # No admissible root anywhere in (0, theta_hi]: the certificate, if it holds at
        # all, holds up to theta_hi.
        valid(c_pi) && (r.eps_max = opt.theta_hi)
    elseif c_pi.pos && eps_star > 0
        c_star = validity_certificate(c, B, X, pairs, eps_star, opt.cert_shrink)
        valid(c_star) && (r.eps_max = eps_star)
    end

    if opt.referee
        r.theta_bisect = referee_theta(c, X, opt.referee_grid, opt.referee_iters, opt.referee_shrink)
    end
    return r
end

# ---------------------------------------------------------------------------

# Everything design_constrained and design_baseline share: solve the shape space, then
# fill the report at whichever X the caller settled on.
function finish(m::Mesh, c::CutStructure, med::Float64, ss::ShapeSpace, X_ini::Vector{Vec2},
                X::Vector{Vec2}, t::Vector{Float64}, opt::DesignOptions, method::String)
    d = DesignResult()
    d.ok = true
    d.status = "ok"
    d.method = method
    d.X = X
    d.X0 = ss.X0
    d.t = t
    d.Phi = ss.Phi
    d.dim_null = ss.k
    d.med_edge = med

    s = med * med
    e = exact_margins(m, c, X, opt.delta_convex * s, opt.delta_split * s)
    d.feasible = e.feasible
    d.min_cross = e.min_cross / s
    d.min_q = e.min_q / s
    d.min_mu = e.min_mu / s
    d.n_bad = e.n_bad
    d.n_bad_q = e.n_bad_q
    d.n_bad_mu = e.n_bad_mu
    d.n_corners = e.n_corners
    d.n_inverted = count_inverted_faces(m, X)

    d.n_nonconvex_x0 = count(v -> v <= 0, corner_crosses(m, ss.X0))

    d.dist_ini = rms_move(X, X_ini, med)
    d.dist_ini_x0 = rms_move(ss.X0, X_ini, med)
    d.t_norm_rel = isempty(t) ? 0.0 : norm(t) / med

    d.ch = characterize(m, m.sigma, X, opt.characterize)
    return d
end

"""
    design_constrained(mesh, sigma, X_ini, opt = DesignOptions()) -> DesignResult

The K9 variant-(b) constrained embedding: the proximity point of the convexity +
split-outward feasible set inside the Eq. (4) null space, then characterized. `X_ini` is
the proximity target and the geometry the Eq. (4) system is assembled from; pass an empty
vector to use `mesh.X`.
"""
function design_constrained(mesh::Mesh, sigma::Vector{Int}, X_ini_in::Vector{Vec2},
                            opt::DesignOptions = DesignOptions())
    X_ini = isempty(X_ini_in) ? mesh.X : X_ini_in
    m = prepare(mesh, sigma, X_ini, "design_constrained")
    c = make_cut(m)
    med = median_edge_length(m)

    d = DesignResult()
    d.method = "constrained"
    if n_split(c) == 0
        d.status = "no split cuts: the shape space is a point, nothing to design"
    end
    ss = shape_space_of(m, c)
    if !ss.ok
        d.status = ss.status
        return d
    end
    if ss.k < 1
        # A split-free (or otherwise rigid) design: X0 is the only member of the shape
        # space, so the constrained problem has a single candidate. Report it as the
        # answer rather than as a failure -- four of the eight authored tilings are like this.
        r = finish(m, c, med, ss, X_ini, ss.X0, Float64[], opt, "constrained")
        r.status = "dim_null = 0: the shape space is the single point X0"
        return r
    end

    # --- K9 variant (b), verbatim ---------------------------------------------
    # Variant (a) (convexity only) and K6's split-only repair are both RELAXATIONS of
    # (b), so their minimisers are strictly better starts than the origin. They are
    # computed only when the cold solve fails, exactly as in apps/kill_k9.cpp.
    ra = nothing
    sp = nothing
    if opt.warm_starts
        oa = embed_options(opt, 0.0)
        ra = convex_embed(c, ss.X0, ss.Phi, X_ini, med, oa)
        ro = ZeroPlusRepairOptions(n_random = opt.restarts, max_iter = opt.max_iter,
                                   start_scale = opt.start_scale, seed = opt.seed)
        sp = zero_plus_repair(c, ss.X0, ss.Phi, med, ro)
    end

    ob = embed_options(opt, opt.delta_split)
    rb = convex_embed(c, ss.X0, ss.Phi, X_ini, med, ob)
    if opt.warm_starts && !rb.feasible && ra.feasible
        o2 = embed_options(opt, opt.delta_split)
        o2.t_init = ra.t
        r2 = convex_embed(c, ss.X0, ss.Phi, X_ini, med, o2)
        if r2.feasible || r2.n_bad + r2.n_bad_q < rb.n_bad + rb.n_bad_q
            rb = r2
        end
    end
    if opt.warm_starts && !rb.feasible && sp.feasible
        o3 = embed_options(opt, opt.delta_split)
        o3.t_init = sp.t
        r3 = convex_embed(c, ss.X0, ss.Phi, X_ini, med, o3)
        if r3.feasible || r3.n_bad + r3.n_bad_q < rb.n_bad + rb.n_bad_q
            rb = r3
        end
    end

    best = finish(m, c, med, ss, X_ini, rb.X, rb.t, opt, "constrained")
    best.barrier_stages_kept = rb.barrier_stages_kept

    # --- optional stage 2: maximise the range, then GATE on the exact margins ---
    if opt.maximise_eps && best.ch.theta_max > opt.characterize.theta_lo
        s = med * med
        rop = deepcopy(opt.range_opt)
        opt.range_step_cap > 0 && (rop.step_cap = opt.range_step_cap * med)
        rr = maximize_range(c, best.X, ss.Phi, rop)
        ef = exact_margins(m, c, rr.X_opt, opt.delta_convex * s, opt.delta_split * s)
        if ef.feasible
            cand = finish(m, c, med, ss, X_ini, rr.X_opt, best.t, opt, "constrained")
            if cand.ch.theta_max > best.ch.theta_max + 1e-9
                cand.barrier_stages_kept = best.barrier_stages_kept
                cand.used_range_opt = true
                best = cand
            end
        end
    end
    return best
end

"""
    design_baseline(mesh, sigma, X_ini, opt = DesignOptions()) -> DesignResult

The 2026 Eq. (6) projection alone (t = 0), characterized identically.
"""
function design_baseline(mesh::Mesh, sigma::Vector{Int}, X_ini_in::Vector{Vec2},
                         opt::DesignOptions = DesignOptions())
    X_ini = isempty(X_ini_in) ? mesh.X : X_ini_in
    m = prepare(mesh, sigma, X_ini, "design_baseline")
    c = make_cut(m)
    med = median_edge_length(m)

    ss = shape_space_of(m, c)
    if !ss.ok
        d = DesignResult()
        d.method = "baseline"
        d.status = ss.status
        return d
    end
    return finish(m, c, med, ss, X_ini, ss.X0, zeros(2 * ss.k), opt, "baseline")
end

# ---------------------------------------------------------------------------
# (D) The range-maximising construction (K9c), as apps/kill_k9c.cpp runs it.

# One scored candidate point of design_range_max. `delta` is the arm's OWN margin: the
# k9b arm is solved (and judged) at 1e-2, every other arm at 1e-3, exactly as in K9c.
mutable struct RangeCand
    have::Bool
    tag::String
    X::Vector{Vec2}
    t::Vector{Float64}
    feasible::Bool
    margin::Float64
    theta_max::Float64
    eps_max::Float64
    ch::Characterization
end
RangeCand() = RangeCand(false, "", Vec2[], Float64[], false, 0.0, 0.0, 0.0, Characterization())

# K9c's ordering of candidates: exact Theta_max first, the certified eps_max as the
# tie-break. Bit-identical to the `best` loop of apps/kill_k9c.cpp.
function beats(a::RangeCand, b::RangeCand)
    !b.have && return true
    a.theta_max > b.theta_max + 1e-12 && return true
    return abs(a.theta_max - b.theta_max) <= 1e-12 && a.eps_max > b.eps_max
end

"""
    design_range_max(mesh, sigma, X_ini, opt = RangeMaxOptions()) -> RangeMaxResult

The K9c range-maximising construction: arms "k9" and "k9b" (proximity points), stage A
(`range_embed` from each available start, in the order k9, k9b, t = 0), stage B
(`maximize_margin_range` from the stage-A winner, kept only if the exact Theta_max
strictly improved), and the best of {k9, k9b, stage-A/B winner} by exact Theta_max, ties
broken by eps_max. With seed 9300 + 7 * id + which this reproduces results/kill/k9c/k9c.csv.
"""
function design_range_max(mesh::Mesh, sigma::Vector{Int}, X_ini_in::Vector{Vec2},
                          opt::RangeMaxOptions = RangeMaxOptions())
    X_ini = isempty(X_ini_in) ? mesh.X : X_ini_in
    m = prepare(mesh, sigma, X_ini, "design_range_max")
    c = make_cut(m)
    med = median_edge_length(m)
    s = med * med

    out = RangeMaxResult()
    ss = shape_space_of(m, c)
    if !ss.ok
        out.design.method = "range_max"
        out.design.status = ss.status
        return out
    end

    # The DesignOptions the shared `finish` reporter measures every point against.
    rep = DesignOptions(delta_convex = opt.delta_convex, delta_split = opt.delta_split,
                        seed = opt.seed, max_iter = opt.max_iter, characterize = opt.characterize)

    if ss.k < 1
        out.design = finish(m, c, med, ss, X_ini, ss.X0, Float64[], rep, "range_max")
        out.design.status = "dim_null = 0: the shape space is the single point X0"
        out.provenance = "x0"
        out.margin = zero_plus_margin(c, ss.X0, med)[1]
        return out
    end

    # Score one point exactly as K9c's `make_cand` does: exact convexity and split margins
    # at the arm's own delta, the exact 0+ margin, and the exact range + certificate.
    function score(tag::String, X::Vector{Vec2}, t::Vector{Float64}, delta_rel::Float64)
        cd = RangeCand()
        cd.have = true
        cd.tag = tag
        cd.X = X
        cd.t = t
        e = exact_margins(m, c, X, delta_rel * s, delta_rel * s)
        cd.feasible = e.feasible
        cd.margin = zero_plus_margin(c, X, med)[1]
        cd.ch = characterize(m, m.sigma, X, opt.characterize)
        cd.theta_max = cd.ch.theta_max
        cd.eps_max = cd.ch.eps_max
        return cd
    end
    function record(cd::RangeCand, stage_a_winner::Bool)
        push!(out.arms, RangeMaxArm(cd.tag, cd.feasible, cd.margin, cd.theta_max, cd.eps_max,
                                    stage_a_winner))
    end

    # ---- arm "k9": the proximity point at delta = delta' = delta_convex ---------------
    c9 = RangeCand()
    r9 = nothing
    if opt.arm_proximity
        o9 = ConvexEmbedOptions(delta_rel = opt.delta_convex, split_delta_rel = opt.delta_split,
                                n_random = 3, max_iter = opt.max_iter, barrier_stages = 6,
                                seed = opt.seed)
        r9 = convex_embed(c, ss.X0, ss.Phi, X_ini, med, o9)
        c9 = score("k9", r9.X, r9.t, opt.delta_convex)
        record(c9, false)
    end
    r9_feasible = r9 !== nothing && r9.feasible

    # ---- arm "k9b": the proximity point at the wide margin, warm started from k9 -------
    c9b = RangeCand()
    r9b = nothing
    if opt.arm_proximity_wide
        o9b = ConvexEmbedOptions(delta_rel = opt.delta_wide, split_delta_rel = opt.delta_wide,
                                 n_random = 8, max_iter = opt.max_iter, barrier_stages = 10,
                                 seed = opt.seed + UInt32(101))
        r9_feasible && (o9b.t_init = r9.t)
        r9b = convex_embed(c, ss.X0, ss.Phi, X_ini, med, o9b)
        c9b = score("k9b", r9b.X, r9b.t, opt.delta_wide)
        record(c9b, false)
    end
    r9b_feasible = r9b !== nothing && r9b.feasible

    # ---- stage A: maximise the exact 0+ margin from every available start --------------
    # The start ORDER is fixed (k9, k9b, t = 0) and every start is always run: K9c's binary
    # skips starts once a wall-clock budget is exceeded, which is the one place its answer
    # depended on machine load. Nothing here reads the clock.
    c9c = RangeCand()
    starts = Tuple{String,Vector{Float64},Float64}[]
    opt.arm_proximity && r9_feasible && push!(starts, ("k9", r9.t, opt.delta_convex))
    opt.arm_proximity_wide && r9b_feasible && push!(starts, ("k9b", r9b.t, opt.delta_wide))
    opt.start_x0 && push!(starts, ("x0", zeros(2 * ss.k), opt.delta_convex))
    for (tag, t0, drel) in starts
        ro = RangeEmbedOptions(delta_rel = drel, split_delta_rel = drel, kappa = opt.kappa,
                               stages = opt.stages, iter_per_stage = opt.iter_per_stage,
                               seed = opt.seed + UInt32(211), t_init = t0)
        rr = range_embed(c, ss.X0, ss.Phi, med, ro)
        cd = score("k9c/" * tag, rr.X, rr.t, drel)
        # Stage A's own ordering: exact Theta_max, the exact margin as the tie-break.
        better = !c9c.have || cd.theta_max > c9c.theta_max + 1e-12 ||
                 (abs(cd.theta_max - c9c.theta_max) <= 1e-12 && cd.margin > c9c.margin)
        record(cd, false)
        better && (c9c = cd)
    end
    if c9c.have
        for a in out.arms
            a.tag == c9c.tag && (a.stage_a_winner = true)
        end
    end

    # ---- stage B: push Theta_max from the stage-A winner, under exact rejection ------
    if opt.stage_b && c9c.have && c9c.feasible && c9c.margin > 0
        mo = deepcopy(opt.margin_range)
        mo.delta_rel = (c9c.tag == "k9c/k9b") ? opt.delta_wide : opt.delta_convex
        mo.split_delta_rel = mo.delta_rel
        mr = maximize_margin_range(c, c9c.X, ss.Phi, med, mo)
        if mr.improved
            cd = score(c9c.tag * "+B", mr.X, c9c.t, mo.delta_rel)
            record(cd, false)
            if cd.theta_max > c9c.theta_max + 1e-9
                c9c = cd
                out.stage_b_used = true
            end
        end
    end

    # ---- the best of the three arms, with provenance -----------------------------------
    best = RangeCand()
    for cd in (c9, c9b, c9c)
        cd.have || continue
        beats(cd, best) && (best = cd)
    end
    if !best.have
        out.design.method = "range_max"
        out.design.status = "every arm was disabled"
        return out
    end
    startswith(best.tag, "k9c/") || (out.stage_b_used = false)

    final_rep = deepcopy(rep)
    final_rep.delta_convex = (best.tag == "k9b" || startswith(best.tag, "k9c/k9b")) ?
                             opt.delta_wide : opt.delta_convex
    final_rep.delta_split = final_rep.delta_convex
    out.design = finish(m, c, med, ss, X_ini, best.X, best.t, final_rep, "range_max")
    out.design.used_range_opt = out.stage_b_used
    out.provenance = best.tag
    out.margin = best.margin
    return out
end

# ---------------------------------------------------------------------------
# Orientations, for the CLI's --sigma switch.

"""Eq. (1) max-cut sigma (2026 Sec. 4.2), i.e. `assign_orientation_relaxation`."""
function orientation_maxcut(mesh::Mesh, seed::Integer = 20260903, restarts::Int = 4,
                            iters::Int = 300, n_diameters::Int = 90)
    rng = MT19937(seed)
    return assign_orientation_relaxation(mesh, rng, restarts, iters, n_diameters).sigma
end

mutable struct DefectOrientationResult
    sigma::Vector{Int}
    defect::Float64   # D(sigma) = sum_holes ||sum_{e in C} (x_t - x_s)||^2 at mesh.X
    attempts::Int
    accepted::Int
    restart_used::Int
    ok::Bool
end
DefectOrientationResult() = DefectOrientationResult(Int[], 0.0, 0, 0, 0, false)

# (ok, D) of one sigma; ok = false rejects a disconnected Gamma AND isolated faces.
function eval_defect(mesh::Mesh, sigma::Vector{Int})
    m = Mesh(mesh.X, mesh.faces)
    m.sigma = copy(sigma)
    m.periodic = mesh.periodic
    build_topology!(m)
    c = make_cut(m)
    hg = build_hinge_graph(c)
    hg.components != 1 && return (false, 0.0)
    hs = holes_partition(c)
    r = hole_residuals(c, m.X, hs)
    d = 0.0
    for v in r.per_hole
        d += dot(v, v)
    end
    return (true, d)
end

"""
    orientation_defect(mesh, sigma_start, cap, seed = 7000) -> DefectOrientationResult

K5's defect-minimizing sigma: greedy single-face and adjacent-pair flips on D(sigma),
4 starts (Eq. (1)'s sigma plus 3 random), rejecting any sigma whose hinge graph is
disconnected. `restart_used` is 0-based like the C++ (0 = sigma_start).
"""
function orientation_defect(mesh::Mesh, sigma_start::Vector{Int}, cap::Int, seed::Integer = 7000)
    length(sigma_start) != n_faces(mesh) &&
        throw(ArgumentError("orientation_defect: sigma_start has the wrong length"))
    best = DefectOrientationResult()
    best.defect = Inf
    rng = MT19937(seed)
    F = n_faces(mesh)
    adj = dual_graph(mesh)
    dual_edges = Tuple{Int,Int}[]
    for f in 1:F, g in adj[f]
        g > f && push!(dual_edges, (f, g))
    end

    for start in 0:3
        sig = start == 0 ? copy(sigma_start) : [uniform_int(rng, 0, 1) != 0 ? 1 : -1 for _ in 1:F]
        (ok, D) = eval_defect(mesh, sig)
        if !ok
            fixed = false
            for _ in 1:5
                fixed && break
                for f in 1:F
                    sig[f] = uniform_int(rng, 0, 1) != 0 ? 1 : -1
                end
                (ok, D) = eval_defect(mesh, sig)
                fixed = ok
            end
            ok || continue
        end
        attempts = 0
        accepted = 0
        improved = true
        order = collect(1:F)
        while improved && attempts < cap
            improved = false
            shuffle!(order, rng)
            for f in order
                attempts >= cap && break
                attempts += 1
                sig[f] = -sig[f]
                (eok, eD) = eval_defect(mesh, sig)
                if eok && eD < D - 1e-12 * max(1.0, D)
                    D = eD
                    accepted += 1
                    improved = true
                else
                    sig[f] = -sig[f]
                end
            end
            shuffle!(dual_edges, rng)
            for (f, g) in dual_edges
                attempts >= cap && break
                attempts += 1
                sig[f] = -sig[f]
                sig[g] = -sig[g]
                (eok, eD) = eval_defect(mesh, sig)
                if eok && eD < D - 1e-12 * max(1.0, D)
                    D = eD
                    accepted += 1
                    improved = true
                else
                    sig[f] = -sig[f]
                    sig[g] = -sig[g]
                end
            end
        end
        if ok && D < best.defect
            best.sigma = copy(sig)
            best.defect = D
            best.attempts = attempts
            best.accepted = accepted
            best.restart_used = start
            best.ok = true
        end
    end
    return best
end
