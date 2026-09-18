# apps/native_common.jl -- the per-cell code path of apps/kill_native200.jl, factored out
# so that a second driver (kill_regime, WP7) can run the authors' full native pipeline and
# score it with EXACTLY the same instruments.
#
# The authors' CLI is an external program (`cli`), run through `/usr/bin/perl -e 'alarm ...'`
# so the timeout semantics (and the "timed_out" classification at 0.9 x timeout) are
# identical for every driver.

include(joinpath(@__DIR__, "kill_common.jl"))

# ---- our own referee, identical to K2a/K2b/K5/K6 --------------------------
referee_theta(c::K.CutStructure, X::Vector{Vec2}, grid::Int = 4000, iters::Int = 50) =
    K.referee_theta(c, X, grid, iters, 1e-12)

function closed_form_theta(c::K.CutStructure, X::Vector{Vec2})
    B = K.deploy_basis(c, X)
    sd = K.swept_discs(c, B)
    return K.exact_theta_max_overlap(c, B, K.candidate_pairs(c, sd, Float64(pi), true),
                                     1e-9, Float64(pi), 1e-9).theta_max
end

"""Reads a dumped pattern JSON (io::pattern_to_json format): vertices + orientation.
Returns `nothing` if the file is missing, malformed, or the vertex count doesn't match;
otherwise `(X, sigma)` with `sigma` empty when the file has no orientation."""
function read_pattern(path::AbstractString, expect_n::Int)
    isfile(path) || return nothing
    j = try
        JSON.parsefile(path)
    catch
        return nothing
    end
    (j isa AbstractDict && haskey(j, "vertices")) || return nothing
    X = Vec2[]
    for v in j["vertices"]
        length(v) < 2 && return nothing
        a = Float64(v[1])
        b = Float64(v[2])
        (isfinite(a) && isfinite(b)) || return nothing
        push!(X, Vec2(a, b))
    end
    length(X) != expect_n && return nothing
    sigma = haskey(j, "orientation") ? Int[Int(s) for s in j["orientation"]] : Int[]
    return X, sigma
end

function jget(j, path, dflt::Float64 = -1.0)
    cur = j
    for k in path
        (cur isa AbstractDict && haskey(cur, k)) || return dflt
        cur = cur[k]
    end
    return cur isa Real && !(cur isa Bool) ? Float64(cur) : dflt
end

Base.@kwdef mutable struct Row
    id::Int = -1
    kind::String = ""
    variant::String = ""
    N::Int = 0
    F::Int = 0
    dim_null::Int = -1
    status::String = ""      # completed / timed_out / crashed
    native_theta_collisions::Float64 = -1.0
    native_finite::Bool = false
    embedding_ok::Bool = false
    our_theta_exact::Float64 = -1.0
    our_theta_bisect::Float64 = -1.0
    cert_pos::Bool = false
    cert_noovl::Bool = false
    cert_noroot::Bool = false
    secs::Float64 = 0.0
end

Base.@kwdef mutable struct NativeRun
    status::String = "crashed"
    theta_collisions::Float64 = -1.0
    finite::Bool = false
    dim_null::Int = -1
    has_pattern::Bool = false
end

# Runs a shell line; exit status returned (non-zero on failure).
function _system(cmd::AbstractString)
    p = run(ignorestatus(`/bin/sh -c $cmd`))
    return p.exitcode
end

function run_prevent_native(cli::AbstractString, input::AbstractString, gdir::AbstractString,
                            timeout_s::Int)
    r = NativeRun()
    mkpath(gdir)
    pat = gdir * "/pat.json"
    full = gdir * "/full.json"
    log = gdir * "/log.txt"
    cmd = "/usr/bin/perl -e 'alarm shift; exec @ARGV' $timeout_s $cli prevent $input " *
          "--collisions --dump $pat --out $full > $log 2>&1"
    rc = _system(cmd)
    if rc != 0
        r.status = "crashed"
        isfile(full) || return r
    end
    isfile(full) || (r.status = "crashed"; return r)
    j = try
        JSON.parsefile(full)
    catch
        r.status = "crashed"
        return r
    end
    r.status = "completed"
    r.dim_null = trunc(Int, jget(j, ("dim_null_per_coordinate",), -1.0))
    if haskey(j, "note")
        r.theta_collisions = jget(j, ("theta_max_X0_with_collisions",), -1.0)
        r.finite = true
        r.has_pattern = false
    elseif haskey(j, "default")
        r.theta_collisions = jget(j, ("default", "theta_max_with_collisions"), -1.0)
        r.finite = Bool(get(j["default"], "finite", false))
        r.has_pattern = isfile(pat)
    end
    return r
end

"""Loads K5's defect-minimising sigma for one graph; empty if absent/mismatched."""
function load_sigma_def(sigmadir::AbstractString, kind::AbstractString, id::Int, m0::K.Mesh)
    p = sigmadir * "/" * kind * "_" * string(id) * ".json"
    isfile(p) || return Int[]
    sm = K.load_mesh_json(p)
    (length(sm.sigma) != K.n_faces(m0) || length(sm.X) != length(m0.X)) && return Int[]
    return sm.sigma
end

function _copy_file(src, dst)
    try
        cp(src, dst; force = true)
    catch
    end
end

# ---- one cell --------------------------------------------------------------
# Their `color` pre-step for the native variant, their `prevent` at published defaults,
# their theta_max_with_collisions as a diagnostic column, then OUR exact T4.2" scan, OUR
# bisection referee and the eps = 0.3 validity certificate on the embedding they dumped.
#
# Returns `(row, X, sigma)`: `X` / `sigma` are the embedding the CLI produced and the sigma
# it was scored under, `nothing` when the cell did not
# complete or the dump could not be read.
function process_cell(id::Int, kind::AbstractString, m0::K.Mesh,
                      sigma::Union{Vector{Int},Nothing}, variant::AbstractString,
                      cli::AbstractString, work::AbstractString, timeout_s::Int,
                      logdir::AbstractString)
    t = Timer()
    r = Row(id = id, kind = String(kind), variant = String(variant),
            N = K.n_vertices(m0), F = K.n_faces(m0))

    gdir = work * "/" * kind * "_" * string(id) * "_" * variant
    rm(gdir; force = true, recursive = true)
    mkpath(gdir)
    input_path = gdir * "/g.json"

    if sigma !== nothing
        m = K.Mesh(m0.X, m0.faces)
        m.sigma = copy(sigma)
        K.save_mesh_json(m, input_path)
    else
        # native variant: no orientation field, then their own coloring.
        raw = K.Mesh(m0.X, m0.faces)
        K.save_mesh_json(raw, gdir * "/raw.json")
        colored = gdir * "/colored.json"
        colorcmd = "/usr/bin/perl -e 'alarm shift; exec @ARGV' $timeout_s $cli color " *
                   "$gdir/raw.json --out $colored > $gdir/color_log.txt 2>&1"
        rc = _system(colorcmd)
        if rc != 0 || !isfile(colored)
            r.status = "crashed"
            r.secs = s(t)
            _copy_file(gdir * "/color_log.txt",
                       logdir * "/" * kind * "_" * string(id) * "_native_color.log")
            rm(gdir; force = true, recursive = true)
            return r, nothing, nothing
        end
        input_path = colored
    end

    nr = run_prevent_native(cli, input_path, gdir, timeout_s)
    secs = s(t)

    if nr.status != "completed"
        src = gdir * "/log.txt"
        isfile(src) && _copy_file(src, logdir * "/" * kind * "_" * string(id) * "_" * variant * ".log")
    end

    status = nr.status
    (status == "crashed" && secs >= 0.9 * timeout_s) && (status = "timed_out")

    our_sigma = Int[]
    if sigma !== nothing
        our_sigma = copy(sigma)
    elseif isfile(input_path)
        cj = try
            JSON.parsefile(input_path)
        catch
            Dict{String,Any}()
        end
        if cj isa AbstractDict && haskey(cj, "orientation")
            our_sigma = Int[Int(x) for x in cj["orientation"]]
        end
    end

    r.dim_null = nr.dim_null
    r.status = status
    r.native_theta_collisions = nr.theta_collisions
    r.native_finite = nr.finite
    r.secs = secs

    X_out = nothing
    sigma_out = nothing
    if status == "completed" && length(our_sigma) == K.n_faces(m0)
        mm = K.Mesh(m0.X, m0.faces)
        mm.sigma = our_sigma
        K.build_topology!(mm)
        c = K.make_cut(mm)
        X = Vec2[]
        have_X = false
        if nr.has_pattern
            rp = read_pattern(gdir * "/pat.json", K.n_vertices(m0))
            if rp !== nothing
                X = rp[1]
                have_X = true
            end
        end
        if !have_X && nr.dim_null == 0
            hs = K.holes_partition(c)
            sys = K.assemble_system(c, hs, mm.X, K.Fixed)
            sr = K.solve_system(sys, mm.X)
            if sr.projection_ok
                X = K.matrix_to_points(sr.X0)
                have_X = true
            end
        end
        if have_X && K.n_split(c) + K.n_hinge(c) > 0
            r.embedding_ok = true
            r.our_theta_exact = closed_form_theta(c, X)
            r.our_theta_bisect = referee_theta(c, X)
            vc = K.validity_certificate(c, X, 0.3)
            r.cert_pos = vc.pos
            r.cert_noovl = vc.nooverlap
            r.cert_noroot = vc.noroot
            X_out = X
            sigma_out = our_sigma
        end
    end
    rm(gdir; force = true, recursive = true)
    return r, X_out, sigma_out
end

function write_row(csv::IO, r::Row, tag::AbstractString)
    print(csv, r.id, ",", r.kind, ",", r.variant, ",", r.N, ",", r.F, ",",
          r.dim_null, ",", r.status, ",", fmt_g(r.native_theta_collisions), ",",
          fmt_b(r.native_finite), ",", fmt_b(r.embedding_ok), ",",
          fmt_g(r.our_theta_exact), ",", fmt_g(r.our_theta_bisect), ",", fmt_b(r.cert_pos),
          ",", fmt_b(r.cert_noovl), ",", fmt_b(r.cert_noroot), ",",
          fmt_b(r.cert_pos && r.cert_noovl && r.cert_noroot), ",", fmt_g(r.secs), ",",
          tag, "\n")
    flush(csv)
    return nothing
end

# CSV split: a trailing empty field is dropped (deliberate, the CSV convention).
function split_csv(line::AbstractString)
    f = String.(split(line, ','))
    (!isempty(f) && isempty(f[end]) && endswith(line, ',')) && pop!(f)
    return f
end
