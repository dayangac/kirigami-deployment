# K1c -- the false-negative rate of 2026 Eq. (9) (ideas/ranking.md Sec. 4, K1c).
# Port of code/apps/kill_k1c.cpp; the CORRECTION recorded there (derivations/core.md T5.3)
# applies verbatim: a root of the separation harmonic is a re-closure only if r > 0,
# p < 0, theta* = 2 arctan(-r/p) is in range, the INTERVAL clause E2 holds at theta*
# and the CROSSING clause holds (the faces' interiors overlap just after theta*).
# The closed forms of T5.2 (p = -sigma_f det(d, du), q = -p, r = <d, du>) are
# cross-checked against the generic orientation harmonic and p + q = 0 is asserted.
# "Range of interest" is reported against Theta_max(Y_9) (`binding`) and min_e beta_e
# (`below_beta`). Y_9 comes from our optimize_collision_sweep and from the authors'
# native `prevent` (baseline/native, published defaults), in separate columns.
#
#   julia --project=Kirigami Kirigami/apps/kill_k1c.jl [--out DIR] [--cli PATH] [--work DIR]
#         [--no-native] [--limit K] [--regenerate]
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const CPP_REPO = expanduser("~/Documents/kirigami-experiments")   # PORTING.md: the reference repo

function referee_theta(c::K.CutStructure, X::Vector{Vec2}, grid::Int = 4000, iters::Int = 50)
    col(th) = K.has_collision(c, K.deploy(c, X, th).Y, 1e-12)
    col(1e-7) && return 0.0
    lo = 0.0; hi = -1.0
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
        if col(mid); hi = mid; else; lo = mid; end
    end
    return lo
end

valid_flat(c::K.CutStructure, X::Vector{Vec2}) = !K.has_collision(c, K.deploy(c, X, 0.0).Y, 1e-12)

mutable struct SplitScan
    n_split::Int
    n_r_pos::Int        # Eq. (9) certifies separation at first order
    n_root::Int         # ... and a second root tau* > 0 exists (p < 0)
    n_interval::Int     # ... and the interval clause E2 holds there
    n_crossing::Int     # ... and the faces actually cross there
    n_below_beta::Int   # ... and theta* < min beta
    n_binding::Int      # ... and theta* is the binding contact (== Theta_max)
    n_degenerate::Int   # split edges collapsed to zero length in this embedding
    n_naive::Int        # what the rule AS WRITTEN counts: r>0, p<0, theta* < Theta_max
    n_corr_tmax::Int    # all clauses, but with the rule's OWN reference range Theta_max
    worst_pq::Float64      # worst |p + q| / scale  (T5.2 asserts p + q == 0)
    worst_closed::Float64  # worst |closed form - generic harmonic| / scale
    first_theta::Float64   # smallest qualifying theta*
    worst_bisect::Float64  # |theta* - bisection of the same face pair| for qualifiers
end
SplitScan() = SplitScan(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, -1.0, 0.0)

const J_ROT = K.Mat2(0.0, 1.0, -1.0, 0.0)   # column-major: [0 -1; 1 0]

function scan_split_edges(m::K.Mesh, c::K.CutStructure, X::Vector{Vec2}, theta_max_Y9::Float64, min_beta::Float64)
    s = SplitScan()
    B = K.deploy_basis(c, X)
    lo = X[1]; hi = X[1]
    for q in X
        lo = min.(lo, q); hi = max.(hi, q)
    end
    diam = norm(hi - lo)
    for ei in c.split_edges
        ed = m.edges[ei]
        ed.n_faces == 2 || continue
        f = m.half_edges[ed.he[1]].face; g = m.half_edges[ed.he[2]].face
        a = ed.key.a; b = ed.key.b
        a1 = K.prime_vertex(c, f, a); b1 = K.prime_vertex(c, f, b); a2 = K.prime_vertex(c, g, a)
        (a1 < 1 || b1 < 1 || a2 < 1) && continue
        s.n_split += 1
        # A split edge collapsed to zero length carries no separation information.
        if norm(X[b] - X[a]) < 1e-9 * diam
            s.n_degenerate += 1
            continue
        end

        # Generic harmonic of (a', b'; a''), and the T5.2 closed form.
        h = K.orient_harmonic(B, a1, b1, a2)
        d = X[b] - X[a]
        # T5.2's offset between the duplicates is 2 sin(theta/2) J du; in this basis it is
        # sin(theta/2) (S(a'') - S(a')), so du = -0.5 J (S(a'') - S(a')).
        off = K.basis_s(B, a2) - K.basis_s(B, a1)
        du = -0.5 * (J_ROT * off)
        p_cf = -Float64(m.sigma[f]) * (d[1] * du[2] - d[2] * du[1])
        r_cf = dot(d, du)
        # Normalise by the GEOMETRIC scale |d| |offset|, not by h.scale().
        sc = max(1e-300, max(K.scale(h), norm(d) * norm(off)))
        s.worst_pq = max(s.worst_pq, abs(h.p + h.q) / sc)
        s.worst_closed = max(s.worst_closed, (abs(h.p - p_cf) + abs(h.q + p_cf) + abs(h.r - r_cf)) / sc)

        p = h.p; r = h.r
        (r > 0) || continue   # Eq. (9) did not certify this edge at first order
        s.n_r_pos += 1
        (p < 0) || continue   # no second root in (0, pi)
        s.n_root += 1
        th = 2.0 * atan(-r / p)
        (th > 1e-9 && th <= pi) || continue
        th < theta_max_Y9 - 1e-9 && (s.n_naive += 1)   # the rule as written

        # Clause E2: the duplicates meet as SEGMENTS at theta*, not merely as lines.
        D = K.dot_harmonic(B, a1, b1, a2); L2 = K.len2_harmonic(B, a1, b1)
        sdot = K.harmonic_eval(D, th); l2 = K.harmonic_eval(L2, th)
        tol = 1e-9 * max(1e-300, K.scale(L2))
        (sdot < -tol || sdot > l2 + tol) && continue
        s.n_interval += 1

        # Crossing clause: the two faces' interiors must overlap just after theta*.
        dth = 1e-5
        function pair_overlaps(t)
            Y = K.basis_eval(B, t)
            A = [Y[pv] for pv in c.prime_faces[f]]
            Bp = [Y[pv] for pv in c.prime_faces[g]]
            return K.polygons_overlap(A, Bp, 1e-12)
        end
        (th + dth <= pi && pair_overlaps(th + dth) && !pair_overlaps(max(1e-9, th - dth))) || continue
        s.n_crossing += 1

        # Cross-check theta* against a bisection on THIS face pair alone.
        let lo = max(1e-9, th - 0.05), hi = min(Float64(pi), th + 0.05)
            if !pair_overlaps(lo) && pair_overlaps(hi)
                for _ in 1:60
                    mid = 0.5 * (lo + hi)
                    if pair_overlaps(mid); hi = mid; else; lo = mid; end
                end
                s.worst_bisect = max(s.worst_bisect, abs(th - lo))
            end
        end

        th < theta_max_Y9 - 1e-9 && (s.n_corr_tmax += 1)
        th < min_beta - 1e-9 && (s.n_below_beta += 1)
        abs(th - theta_max_Y9) <= 1e-5 && (s.n_binding += 1)
        (s.first_theta < 0 || th < s.first_theta) && (s.first_theta = th)
    end
    return s
end

mutable struct Tally
    designs::Int; with_reclosure_beta::Int; with_reclosure_binding::Int
    with_naive::Int; with_corr_tmax::Int
end
Tally() = Tally(0, 0, 0, 0, 0)

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    outdir = joinpath(REPO, "results", "kill", "k1c_julia")
    cli = joinpath(CPP_REPO, "baseline", "native", "build", "tuttekiri_cli")
    work = "/tmp/kiri_k1c_julia"
    use_native = true; limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--cli" && i < length(args); cli = args[i+1]; i += 2
        elseif a == "--work" && i < length(args); work = args[i+1]; i += 2
        elseif a == "--no-native"; use_native = false; i += 1
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir); mkpath(work)
    csv = open(joinpath(outdir, "k1c.csv"), "w")
    print(csv, "name,family,N,F,n_split,dim_null,source,valid_Y9,theta_Y9,min_beta,",
          "n_degenerate,n_r_pos,n_root,n_interval,n_crossing,n_corr_tmax,n_below_beta,n_binding,n_naive,",
          "first_theta,worst_pq,worst_closed,worst_bisect\n")
    wall = Timer()

    t_ours = Tally(); t_native = Tally()
    worst_pq = 0.0; worst_closed = 0.0; worst_bisect = 0.0
    n_seen = 0; n_degen_ours = 0; n_degen_native = 0

    function emit(name, fam, m, c, k, src, Y9, tal)
        ok = valid_flat(c, Y9)
        th9 = ok ? referee_theta(c, Y9) : 0.0
        beta = K.hinge_beta(c, Y9)
        mb = isempty(beta) ? 0.0 : minimum(beta)
        s = scan_split_edges(m, c, Y9, th9, mb)
        print(csv, name, ",", fam, ",", K.n_vertices(m), ",", K.n_faces(m), ",",
              K.n_split(c), ",", k, ",", src, ",", cpp_b(ok), ",", cpp_g(th9), ",", cpp_g(mb),
              ",", s.n_degenerate, ",", s.n_r_pos, ",", s.n_root, ",", s.n_interval, ",", s.n_crossing,
              ",", s.n_corr_tmax, ",", s.n_below_beta, ",", s.n_binding, ",", s.n_naive, ",",
              cpp_g(s.first_theta), ",", cpp_g(s.worst_pq), ",", cpp_g(s.worst_closed), ",", cpp_g(s.worst_bisect), "\n")
        flush(csv)
        worst_pq = max(worst_pq, s.worst_pq)
        worst_closed = max(worst_closed, s.worst_closed)
        worst_bisect = max(worst_bisect, s.worst_bisect)
        if src == "native_eq9"; n_degen_native += s.n_degenerate; else; n_degen_ours += s.n_degenerate; end
        ok || return   # only designs Eq. (9) actually certified count
        tal.designs += 1
        tal.with_reclosure_beta += (s.n_below_beta > 0)
        tal.with_reclosure_binding += (s.n_binding > 0)
        tal.with_naive += (s.n_naive > 0)
        tal.with_corr_tmax += (s.n_corr_tmax > 0)
    end

    nrun = 0
    for d in deployable_population(regenerate = regenerate)
        nrun >= limit && break
        m = d.mesh
        K.build_topology!(m)
        c = K.make_cut(m)
        (K.n_split(c) == 0 || d.dim_null < 1) && continue
        hs = K.holes_partition(c)
        sr = K.solve_system(K.assemble_system(c, hs, m.X, K.Fixed), m.X)
        (!sr.projection_ok || sr.dim_null < 1) && continue
        n_seen += 1
        nrun += 1

        # (a) our Eq. (9) reimplementation, full gamma ladder
        eq9 = K.optimize_collision_sweep(c, d.X, sr.Phi)
        emit(d.name, d.family, m, c, sr.dim_null, "ours_eq9", eq9.X_opt, t_ours)

        # (b) the authors' native prevent_intersections at the published defaults
        if use_native
            gdir = joinpath(work, d.name)
            mkpath(gdir)
            mm = deepcopy(m)
            mm.X = d.X
            K.save_mesh_json(mm, joinpath(gdir, "g.json"))
            cmd = `$cli prevent $(joinpath(gdir, "g.json")) --dump $(joinpath(gdir, "y9.json")) --out $(joinpath(gdir, "rep.json"))`
            rc = success(pipeline(ignorestatus(cmd), stdout = devnull, stderr = devnull))
            if rc && isfile(joinpath(gdir, "y9.json"))
                j = try
                    JSON.parsefile(joinpath(gdir, "y9.json"))
                catch
                    Dict{String,Any}()
                end
                if haskey(j, "vertices") && length(j["vertices"]) == K.n_vertices(m)
                    Y = Vec2[]
                    fin = true
                    for v in j["vertices"]
                        x = Float64(v[1]); y = Float64(v[2])
                        if !isfinite(x) || !isfinite(y)
                            fin = false
                            break
                        end
                        push!(Y, Vec2(x, y))
                    end
                    fin && emit(d.name, d.family, m, c, sr.dim_null, "native_eq9", Y, t_native)
                end
            end
        end
        n_seen % 10 == 0 && println("  ", n_seen, " designs, ", fx(s(wall), 1), " s")
    end
    close(csv)

    pct(a, b) = b > 0 ? 100.0 * a / b : 0.0
    sm = open(joinpath(outdir, "summary.txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("K1c: $n_seen deployable, embedded designs with split cuts\n")
    both("closed forms T5.2 verified: worst |p+q|/scale " * sci(worst_pq) * ", worst |closed - generic|/scale " * sci(worst_closed) * "\n")
    both("worst |theta* - per-pair bisection| for qualifying re-closures : " * sci(worst_bisect) * "\n")
    both("split edges collapsed to zero length, ours / native Eq. (9) : $n_degen_ours / $n_degen_native\n")
    for (T, tag) in ((t_ours, "our optimize_collision_sweep (gamma ladder)"),
                     (t_native, "authors' native prevent (published defaults)"))
        T.designs == 0 && continue
        both("--- Y_9 from " * tag * " ---\n")
        both("  certified designs (valid flat embedding)      : $(T.designs)\n")
        both("  with a re-closure below min beta (CORRECTED)  : $(T.with_reclosure_beta) (" * fx(pct(T.with_reclosure_beta, T.designs), 2) * " %)\n")
        both("  where the re-closure is the binding contact   : $(T.with_reclosure_binding) (" * fx(pct(T.with_reclosure_binding, T.designs), 2) * " %)\n")
        both("  rule AS WRITTEN (no interval/crossing clause) : $(T.with_naive) (" * fx(pct(T.with_naive, T.designs), 2) * " %)\n")
        both("  clauses applied, SAME range Theta_max         : $(T.with_corr_tmax) (" * fx(pct(T.with_corr_tmax, T.designs), 2) * " %)   <- isolates the clauses' effect\n")
        both("  VERDICT (corrected, >= 3% below min beta)    : " * (pct(T.with_reclosure_beta, T.designs) >= 3.0 ? "PASS" : "FAIL") * "\n")
    end
    both("wall " * fx(s(wall), 1) * " s\n")
    close(sm)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
