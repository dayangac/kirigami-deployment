# K2b -- the range-optimisation margin (ideas/ranking.md Sec. 4, K2b).
# Port of code/apps/kill_k2b.cpp; the amended PASS rule, the three corrections and the
# two recorded deviations there apply verbatim:
#   baseline = the authors' native tuttekiri_cli `prevent` (their Eq. 9) with the
#   parameter ladder, best kept, every point refereed by OUR bisection (grid 4000,
#   shrink 1e-12) and required to be a valid flat embedding; ours = maximize_range;
#   PASS if the median relative gain on live graphs is >= 25% with >= 20 certified.
#   Population: the 4 split-carrying reference cases, 8 random graphs (documented as
#   vacuous), and deployable_population() members with dim_null >= 2.
#
#   julia --project=Kirigami Kirigami/apps/kill_k2b.jl [--random 8] [--deployable 26]
#         [--rounds 20] [--maxf 160] [--out DIR] [--cache DIR] [--cli PATH] [--work DIR]
#         [--limit K] [--regenerate]
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const CPP_REPO = expanduser("~/Documents/kirigami-experiments")   # PORTING.md: the reference repo

# The referee. Bisection on the true polygon overlap, shrink 1e-12, fine grid.
function referee_theta(c::K.CutStructure, X::Vector{Vec2}, grid::Int = 4000, iters::Int = 50)
    col(th) = K.has_collision(c, K.deploy(c, X, th).Y, 1e-12)
    col(1e-7) && return 0.0   # penetrates immediately
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

# Validity: no face-face overlap in the FLAT state (T5.1: positive orientation is not a
# certificate).
valid_flat(c::K.CutStructure, X::Vector{Vec2}) = !K.has_collision(c, K.deploy(c, X, 0.0).Y, 1e-12)

function closed_form_theta(c::K.CutStructure, X::Vector{Vec2})
    B = K.deploy_basis(c, X)
    sd = K.swept_discs(c, B)
    return K.exact_theta_max_overlap(c, B, K.candidate_pairs(c, sd, Float64(pi), true), 1e-9, Float64(pi), 1e-9).theta_max
end

function read_vertices(path::AbstractString, expect_n::Int)
    isfile(path) || return Vec2[]
    j = try
        JSON.parsefile(path)
    catch
        return Vec2[]
    end
    (j isa AbstractDict && haskey(j, "vertices")) || return Vec2[]
    X = Vec2[]
    for v in j["vertices"]
        length(v) < 2 && return Vec2[]
        a = Float64(v[1]); b = Float64(v[2])
        (isfinite(a) && isfinite(b)) || return Vec2[]
        push!(X, Vec2(a, b))
    end
    return length(X) == expect_n ? X : Vec2[]
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    n_random = 8; n_deployable = 26; rounds = 20; maxf = 160; limit = typemax(Int)
    outdir = joinpath(REPO, "results", "kill", "k2b_julia")
    cache = joinpath(REPO, "results", "kill", "cache")
    cli = joinpath(CPP_REPO, "baseline", "native", "build", "tuttekiri_cli")
    work = "/tmp/kiri_k2b_julia"
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--random" && i < length(args); n_random = arg_i(args[i+1]); i += 2
        elseif a == "--deployable" && i < length(args); n_deployable = arg_i(args[i+1]); i += 2
        elseif a == "--rounds" && i < length(args); rounds = arg_i(args[i+1]); i += 2
        elseif a == "--maxf" && i < length(args); maxf = arg_i(args[i+1]); i += 2
        elseif a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--cache" && i < length(args); cache = args[i+1]; i += 2
        elseif a == "--cli" && i < length(args); cli = args[i+1]; i += 2
        elseif a == "--work" && i < length(args); work = args[i+1]; i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir); mkpath(work)
    csv = open(joinpath(outdir, "k2b.csv"), "w")
    print(csv, "name,kind,N,F,n_split,dim_null,ladder_points,ladder_valid,",
          "theta_ourX0,valid_ourX0,theta_theirX0,valid_theirX0,",
          "theta_native_best,native_best_idx,theta_ours_eq9,theta_ours,valid_ours,",
          "theta_closed_ours,rel_gain,gap_closed,min_beta,n_active,secs_native,secs_ours\n")
    wall = Timer()

    gains = Float64[]; gains_live = Float64[]
    graphs = 0; certified = 0; ours_better = 0; base_better = 0; tie = 0; vacuous = 0
    snub_gap = -1.0; snub_base = -1.0; snub_ours = -1.0; snub_beta = -1.0

    function run_one(name::String, kind::String, m::K.Mesh, id::Int, full_ladder::Bool,
                     X_given::Union{Vector{Vec2},Nothing} = nothing)
        K.build_topology!(m)
        isempty(m.sigma) && return
        c = K.make_cut(m)
        K.n_split(c) == 0 && return
        hs = K.holes_partition(c)
        sh = shape_space(m, c, hs, cache, id)
        (!sh.ok || sh.k < 2) && return
        # The starting embedding: the member's own X for a deployable_population member,
        # else the Eq. (6) projection. The native baseline is fed the same X.
        X0 = X_given === nothing ? K.matrix_to_points(sh.X0) : X_given
        X_given === nothing || (m.X = X_given)

        # ---- the authors' native baseline -------------------------------------------
        gdir = joinpath(work, name)
        rm(gdir; force = true, recursive = true)
        mkpath(gdir)
        K.save_mesh_json(m, joinpath(gdir, "g.json"))
        tn = Timer()
        cmd = full_ladder ?
            `$cli prevent $(joinpath(gdir, "g.json")) --sweep --dumpdir $(joinpath(gdir, "sw")) --out $(joinpath(gdir, "sw.json"))` :
            `$cli prevent $(joinpath(gdir, "g.json")) --sweep --short --dumpdir $(joinpath(gdir, "sw")) --out $(joinpath(gdir, "sw.json"))`
        rc = success(pipeline(ignorestatus(cmd), stdout = devnull, stderr = devnull))
        sn = s(tn)
        theta_native = -1.0; theta_theirX0 = -1.0
        best_idx = -1; n_points = 0; n_valid = 0
        valid_theirX0 = false
        if rc && isdir(joinpath(gdir, "sw"))
            # std::filesystem::directory_iterator order is unspecified; on APFS it is the
            # directory's hash order. Ties in theta are resolved by that order in the C++
            # (first best wins), so the winner index may differ on exact ties.
            for fn in readdir(joinpath(gdir, "sw"))
                Xb = read_vertices(joinpath(gdir, "sw", fn), K.n_vertices(m))
                isempty(Xb) && continue
                n_points += 1
                # Every candidate must be a VALID flat embedding before its range counts.
                valid_flat(c, Xb) || continue
                n_valid += 1
                th = referee_theta(c, Xb)
                if fn == "run_X0.json"
                    theta_theirX0 = th; valid_theirX0 = true
                end
                if th > theta_native
                    theta_native = th
                    best_idx = fn == "run_X0.json" ? -2 : something(tryparse(Int, match(r"^\D*(\d+)", fn).captures[1]), 0)
                end
            end
        end
        theta_native < 0 && (theta_native = 0.0)   # no ladder point is a valid embedding

        # ---- our Eq. (9) reimplementation, secondary column ---------------------------
        eq9 = K.optimize_collision_sweep(c, X0, sh.Phi)
        theta_eq9 = valid_flat(c, eq9.X_opt) ? referee_theta(c, eq9.X_opt) : 0.0

        # ---- ours ---------------------------------------------------------------------
        to = Timer()
        o = K.RangeOptOptions()
        o.rounds = rounds
        ro = K.maximize_range(c, X0, sh.Phi, o)
        so = s(to)
        vo = valid_flat(c, ro.X_opt)
        theta_ourX0 = valid_flat(c, X0) ? referee_theta(c, X0) : 0.0
        # T = 0 fallback: never report worse than the input we started from.
        theta_ours = max(theta_ourX0, vo ? referee_theta(c, ro.X_opt) : 0.0)
        theta_closed = vo ? closed_form_theta(c, ro.X_opt) : 0.0

        beta = K.hinge_beta(c, X0)
        mb = isempty(beta) ? 0.0 : minimum(beta)
        base = theta_native
        rel = base > 1e-9 ? (theta_ours - base) / base : (theta_ours > 1e-9 ? Inf : 0.0)
        gap = (mb - base > 1e-9) ? (theta_ours - base) / (mb - base) : 0.0

        print(csv, name, ",", kind, ",", K.n_vertices(m), ",", K.n_faces(m), ",",
              K.n_split(c), ",", sh.k, ",", n_points, ",", n_valid, ",",
              cpp_g(theta_ourX0), ",", cpp_b(valid_flat(c, X0)), ",", cpp_g(theta_theirX0), ",",
              cpp_b(valid_theirX0), ",", cpp_g(theta_native), ",", best_idx, ",", cpp_g(theta_eq9), ",",
              cpp_g(theta_ours), ",", cpp_b(vo), ",", cpp_g(theta_closed), ",", cpp_g(rel), ",", cpp_g(gap), ",",
              cpp_g(mb), ",", ro.n_active, ",", cpp_g(sn), ",", cpp_g(so), "\n")
        flush(csv)
        graphs += 1
        (vo && theta_ours > 1e-9) && (certified += 1)
        push!(gains, isfinite(rel) ? rel : 10.0)
        # "Live" = the baseline had a valid ladder point with a positive range.
        if base > 1e-9; push!(gains_live, rel); else; vacuous += 1; end
        if theta_ours > base + 1e-6; ours_better += 1
        elseif base > theta_ours + 1e-6; base_better += 1
        else; tie += 1
        end
        if name == "snub_square_33434"
            snub_base = base; snub_ours = theta_ours; snub_beta = mb; snub_gap = gap
        end
        println("  ", name, " F=", K.n_faces(m), " k=", sh.k, " native=", fx(base), " (", n_valid, "/", n_points,
                " valid) eq9=", fx(theta_eq9), " ours=", fx(theta_ours), " minbeta=", fx(mb),
                "  [", fx(sn, 1), "s + ", fx(so, 1), "s]  (", fx(s(wall), 1), " s)")
    end

    for rc in reference_cases(regenerate = regenerate)
        rc.periodic && continue
        graphs >= limit && break
        run_one(rc.name, "reference", rc.mesh, -1, true)
    end
    taken = 0
    for g in population("k2b_random_8"; regenerate = regenerate)
        taken < n_random || break
        graphs >= limit && break
        g.ok || continue
        before = graphs
        run_one(g.kind * "_" * string(g.id), g.kind, g.mesh, g.id, false)
        graphs > before && (taken += 1)
    end
    dtaken = 0
    for d in deployable_population(regenerate = regenerate)
        dtaken >= n_deployable && break
        graphs >= limit && break
        d.dim_null < 2 && continue
        before = graphs
        run_one(d.name, d.family, d.mesh, -1, false, d.X)
        graphs > before && (dtaken += 1)
    end
    close(csv)

    median(v) = isempty(v) ? NaN : sort(v)[length(v) ÷ 2 + 1]
    sm = open(joinpath(outdir, "summary.txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("K2b: $graphs graphs with split cuts and dim_null >= 2\n")
    both("baseline = authors' native tuttekiri_cli `prevent`, parameter ladder, best kept,\n")
    both("           their X0 as the T=0 fallback; every point refereed by OUR bisection\n")
    both("           (grid 4000, shrink 1e-12) and required to be a valid flat embedding.\n")
    both("median relative gain over the native baseline : " * fx(median(gains)) * "\n")
    both("  restricted to graphs where the baseline > 0  : " * fx(median(gains_live)) * " (n=$(length(gains_live)))\n")
    both("graphs where BOTH sides are 0 (vacuous)        : $vacuous\n")
    both("designs with certified validity and range > 0 : $certified/$graphs\n")
    both("ours better / native better / tie             : $ours_better / $base_better / $tie\n")
    snub_base >= 0 && both("snub square: native " * fx(snub_base) * ", ours " * fx(snub_ours) * ", min beta " *
                           fx(snub_beta) * ", gap closed " * fx(100 * snub_gap, 1) * "%\n")
    pass = (median(gains_live) >= 0.25) && (certified >= 20)
    both("VERDICT (median gain on live graphs >= 25% and >= 20 certified): " * (pass ? "PASS" : "FAIL") * "\n")
    both("wall " * fx(s(wall), 1) * " s\n")
    close(sm)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
