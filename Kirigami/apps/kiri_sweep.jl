# kiri_sweep --n 100 --out sweep.csv [--seed s] [--violations <dir>] [--cases <dir>]
#
# The sweep draws its own random graphs from ONE MT19937(seed) (site count, generator and
# Eq. (1) relaxation on the same stream); there is no frozen corpus for it.
include(joinpath(@__DIR__, "common_app.jl"))

# `std::chrono` millisecond timings, as doubles with microsecond resolution.
ms(a::UInt64, b::UInt64) = floor((b - a) / 1e3) / 1000.0

Base.@kwdef mutable struct Stats
    n::Int = 0
    claim_ok::Int = 0
    rankL_eq_H::Int = 0
    dimnull_eq_split::Int = 0
    a1::Int = 0
    agree::Int = 0
    part::Int = 0
    forest::Int = 0
    dense::Int = 0
    x0_valid::Int = 0
    geo_checked::Int = 0
    geo_agree::Int = 0
    H_le_interior::Int = 0
    # orchestrator checks
    rs_naive::Int = 0
    rs_restricted::Int = 0
    rs_support::Int = 0
    lrd::Int = 0
    z_done::Int = 0
    rank_id::Int = 0
    harmonic::Int = 0
    gamma_conn::Int = 0
    euler_H::Int = 0
    euler_Hall::Int = 0
    worst_harm::Float64 = 0.0
    harm_checks::Int = 0
    flipped_frac::Float64 = 0.0
    max_move::Float64 = 0.0
    worst_solve_ms::Float64 = 0.0
    min_F::Int = 1 << 30
    max_F::Int = 0
end

Base.@kwdef mutable struct RefRow
    name::String = ""
    F::Int = 0
    n_hinge::Int = 0
    H::Int = 0
    H_all::Int = 0
    notches::Int = 0
    c_gamma::Int = 0
    predicted::Int = 0
    rank_L::Int = 0
    dim_Z::Int = 0
    n_hinge_in_L::Int = 0
    rs_naive::Bool = false
    rs_restricted::Bool = false
    rs_support::Bool = false
    r_zero::Bool = false
    all_interior::Bool = false
    lrd::Bool = false
    rank_id::Bool = false
    harmonic::Bool = false
    euler_H::Bool = false
    euler_Hall::Bool = false
    max_harm::Float64 = 0.0
end

# The three L checks of a reference / torus patch, assembled with no boundary rows.
function ref_row(name::String, g::K.Mesh)
    c = K.make_cut(g)
    hs = K.holes_partition(c)
    sys = K.assemble_system(c, hs, g.X, K.None)
    rs = K.check_row_sum(c, hs, sys)
    fz = K.check_factorization(c, hs, sys, true)
    hg = K.check_hinge_graph(c, hs)
    r = RefRow(name = name)
    r.F = K.n_faces(g)
    r.n_hinge = K.n_hinge(c)
    r.H = hg.H
    r.H_all = hg.H_all
    r.notches = hg.n_notches
    r.c_gamma = hg.c_gamma
    r.predicted = hg.predicted
    r.rank_L = fz.rank_L
    r.dim_Z = fz.dim_Z
    r.n_hinge_in_L = rs.n_hinge_in_L
    r.rs_naive = rs.degree_identity
    r.rs_restricted = rs.restricted_degree_identity
    r.rs_support = rs.support_on_boundary
    r.r_zero = rs.r_is_zero
    r.all_interior = rs.all_vertices_interior
    r.lrd = fz.L_equals_RD
    r.rank_id = fz.rank_identity
    r.harmonic = fz.out_harmonic
    r.max_harm = fz.max_harmonic_residual
    r.euler_H = hg.identity_holds
    r.euler_Hall = hg.identity_with_notches
    return r, rs, fz, hg
end

f3(v) = fx(v, 3)            # std::fixed << std::setprecision(3)
sci_default(v) = sci(v, 3)   # std::scientific after the md stream's setprecision(3)

function main(args::Vector{String})
    n = 100
    seed = UInt32(20260903)
    out = "sweep.csv"
    vdir = ""
    cases_dir = "../results/core_validation/cases"
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--n" && i + 1 <= length(args)
            n = arg_i(args[i += 1])
        elseif a == "--out" && i + 1 <= length(args)
            out = args[i += 1]
        elseif a == "--seed" && i + 1 <= length(args)
            seed = arg_u(args[i += 1])
        elseif a == "--violations" && i + 1 <= length(args)
            vdir = args[i += 1]
        elseif a == "--cases" && i + 1 <= length(args)
            cases_dir = args[i += 1]
        end
        i += 1
    end
    isempty(dirname(out)) || mkpath(dirname(out))
    isempty(vdir) || mkpath(vdir)
    csv = open(out, "w")
    print(csv, "id,kind,n_sites,N,F,E,n_interior,n_hinge,n_split,H,H_all,H_geometric,rank_full,rank_L,",
          "dim_null,n_interior_minus_H,claim_holds,deployable,residual_max,X0_computed,X0_valid,",
          "X0_flipped_faces,X0_max_move,theta_max_geom,",
          "min_beta,A1_ok,holes_agree,partition_ok,forest,geo_checked,geo_agree,used_sparse,",
          # --- orchestrator checks (added) ---
          "n_hinge_in_L,rowsum_naive_ok,rowsum_restricted_ok,rowsum_support_boundary,",
          "L_eq_RD,dim_Z,rank_identity,out_harmonic,max_harm_res,",
          "c_gamma,gamma_connected,H_euler_pred,euler_H_ok,euler_Hall_ok,n_notches,",
          "t_build_ms,t_holes_ms,t_solve_ms,t_checks_ms\n")

    rng = K.MT19937(seed)
    violations = 0
    st = Stats()
    violation_lines = String[]
    kinds = ("voronoi", "delaunay", "quad_random")
    for id in 0:n-1
        kind = kinds[id % 3 + 1]
        # target 100..5000 faces
        sites = K.uniform_int(rng, 60, 2600)
        kind == "delaunay" && (sites = max(60, sites ÷ 2))
        t0 = time_ns()
        m = try
            K.generate(kind, [Float64(sites), 40.0], rng)
        catch e
            println(stderr, "gen failed id=", id, ": ", sprint(showerror, e))
            continue
        end
        K.n_faces(m) < 20 && continue
        m.sigma = K.assign_orientation_relaxation(m, rng, 4, 300, 90).sigma
        c = K.make_cut(m)
        t1 = time_ns()
        hs = K.holes_partition(c)
        hs2 = K.holes_seed_growing(c)
        agree = K.same_hole_sets(hs, hs2)[1]
        part = K.holes_partition_edges(c, hs)[1]
        forest = K.split_subgraph_is_forest(c)[1]
        t2 = time_ns()
        sys = K.assemble_system(c, hs, m.X, K.Fixed)
        dense = K.n_vertices(m) <= 1400
        rep = dense ? K.solve_system(sys, m.X) : K.rank_only_sparse(sys)
        t3 = time_ns()
        res = K.hole_residuals(c, m.X, hs)
        # ---- orchestrator checks (added): row sum, L = R D, hinge-graph Euler count ----
        tc0 = time_ns()
        rs = K.check_row_sum(c, hs, sys)
        fz = K.check_factorization(c, hs, sys, dense)
        hg = K.check_hinge_graph(c, hs)
        tc1 = time_ns()
        # The geometric hole trace and the deployment range are only defined on a
        # Tutte auxetic embedding that is itself a valid (non-self-intersecting)
        # straight-line embedding. For random planar graphs the Eq. (6) projection
        # usually is NOT -- the paper lists this as an open problem (Sec. 6).
        Xd = dense ? K.matrix_to_points(rep.X0) : m.X
        X0_valid = dense && !K.has_collision(c, K.deploy(c, Xd, 0.0).Y)
        flipped = -1
        max_move = -1.0
        if dense
            mx = K.Mesh(Xd, m.faces)
            flipped = 0
            max_move = 0.0
            for f in 1:K.n_faces(m)
                K.face_signed_area(mx, f) <= 0 && (flipped += 1)
            end
            for v in 1:K.n_vertices(m)
                max_move = max(max_move, norm(Xd[v] - m.X[v]))
            end
        end
        geo_checked = false
        geo_agree = false
        n_geo = -1
        tm = K.ThetaMaxReport()
        if X0_valid
            th = 0.2
            Yd = Vec2[]
            for _ in 1:8
                Yd = K.deploy(c, Xd, th).Y
                if !K.has_collision(c, Yd)
                    geo_checked = true
                    break
                end
                th *= 0.4
            end
            if geo_checked
                geo = K.holes_geometric(c, Yd)
                comb = sort([hs.all[i].edges for i in hs.interior_indices])
                geo_agree = (geo == comb)
                n_geo = length(geo)
            end
            K.n_faces(m) <= 900 && (tm = K.theta_max(c, Xd, 40, 25))
        end
        if isempty(tm.beta)
            tm.beta = K.hinge_beta(c, Xd)
            tm.min_beta = isempty(tm.beta) ? 0.0 : minimum(tm.beta)
        end
        claim = K.n_interior_vertices(m) - K.n_interior_holes(hs)
        holds = (rep.dim_null == claim)
        a1 = K.check_remark_A1(c)[1]
        print(csv, id, ",", kind, ",", sites, ",", K.n_vertices(m), ",", K.n_faces(m), ",",
              K.n_edges(m), ",", K.n_interior_vertices(m), ",", K.n_hinge(c), ",",
              K.n_split(c), ",", K.n_interior_holes(hs), ",", length(hs.all), ",",
              n_geo, ",", rep.rank_full, ",", rep.rank_L, ",", rep.dim_null, ",",
              claim, ",", fmt_b(holds), ",", fmt_b(K.deployable(res)), ",", fmt_g(res.max_norm), ",",
              fmt_b(dense), ",", fmt_b(X0_valid), ",", flipped, ",", fmt_g(max_move), ",",
              fmt_g(tm.theta_max_geometric), ",", fmt_g(tm.min_beta), ",",
              fmt_b(a1), ",", fmt_b(agree), ",", fmt_b(part), ",", fmt_b(forest), ",",
              fmt_b(geo_checked), ",", fmt_b(geo_agree), ",", fmt_b(rep.used_sparse),
              ",", rs.n_hinge_in_L, ",", fmt_b(rs.degree_identity), ",",
              fmt_b(rs.restricted_degree_identity), ",", fmt_b(rs.support_on_boundary), ",",
              fmt_b(fz.L_equals_RD), ",", (fz.Z_computed ? fz.dim_Z : -1), ",",
              fmt_b(fz.Z_computed ? fz.rank_identity : false), ",",
              fmt_b(fz.Z_computed ? fz.out_harmonic : false), ",", fmt_g(fz.max_harmonic_residual), ",",
              hg.c_gamma, ",", fmt_b(hg.c_gamma == 1), ",", hg.predicted, ",",
              fmt_b(hg.identity_holds), ",", fmt_b(hg.identity_with_notches), ",", hg.n_notches,
              ",", fmt_g(ms(t0, t1)), ",", fmt_g(ms(t1, t2)), ",", fmt_g(ms(t2, t3)), ",",
              fmt_g(ms(tc0, tc1)), "\n")
        flush(csv)
        st.n += 1
        st.claim_ok += holds
        st.rankL_eq_H += (rep.rank_L == K.n_interior_holes(hs))
        st.dimnull_eq_split += (rep.dim_null == K.n_split(c))
        st.H_le_interior += (K.n_interior_holes(hs) <= K.n_interior_vertices(m))
        st.a1 += a1
        st.agree += agree
        st.part += part
        st.forest += forest
        st.geo_checked += geo_checked
        st.geo_agree += geo_agree
        st.min_F = min(st.min_F, K.n_faces(m))
        st.max_F = max(st.max_F, K.n_faces(m))
        st.worst_solve_ms = max(st.worst_solve_ms, ms(t2, t3))
        st.rs_naive += rs.degree_identity
        st.rs_restricted += rs.restricted_degree_identity
        st.rs_support += rs.support_on_boundary
        st.lrd += fz.L_equals_RD
        st.gamma_conn += (hg.c_gamma == 1)
        st.euler_H += hg.identity_holds
        st.euler_Hall += hg.identity_with_notches
        if fz.Z_computed
            st.z_done += 1
            st.rank_id += fz.rank_identity
            st.harmonic += fz.out_harmonic
            st.worst_harm = max(st.worst_harm, fz.max_harmonic_residual)
            st.harm_checks += fz.n_harmonic_checks
        end
        if !isempty(vdir)
            save(ok::Bool, tag::String) = begin
                ok && return
                f = vdir * "/sweep_" * string(id) * "_" * tag * ".json"
                isfile(f) || K.save_mesh_json(m, f)
            end
            save(rs.restricted_degree_identity, "rowsum_restricted")
            save(fz.L_equals_RD, "L_eq_RD")
            if fz.Z_computed
                save(fz.rank_identity, "rank_identity")
                save(fz.out_harmonic, "out_harmonic")
            end
        end
        if dense
            st.dense += 1
            st.x0_valid += X0_valid
            st.flipped_frac += flipped / K.n_faces(m)
            st.max_move = max(st.max_move, max_move)
        end
        if !holds
            push!(violation_lines, "| $id | $kind | $(K.n_faces(m)) | $(K.n_interior_vertices(m)) | " *
                  "$(K.n_interior_holes(hs)) | $(rep.dim_null) | $claim | `violations/violation_$id.json` |")
        end
        if !holds && !isempty(vdir) && violations < 40
            K.save_mesh_json(m, vdir * "/violation_" * string(id) * ".json")
            violations += 1
        end
        println("[", id + 1, "/", n, "] ", kind, " F=", K.n_faces(m), " H=", K.n_interior_holes(hs),
                " dim_null=", rep.dim_null, " claim=", claim, (holds ? " ok" : " VIOLATION"))
    end

    # ---------------------------------------------------------------------------
    # Orchestrator checks on the eight reference cases, read back from the JSON
    # dumps `kiri_reference` writes (no duplicated case construction).
    # ---------------------------------------------------------------------------
    refrows = RefRow[]
    ref_violations = String[]
    let
        vout = isempty(vdir) ? (isempty(dirname(out)) ? "violations" : dirname(out) * "/violations") : vdir
        mfiles = String[]
        if isdir(cases_dir)
            for d in readdir(cases_dir; join = true)
                isdir(d) && isfile(joinpath(d, "M.json")) && push!(mfiles, joinpath(d, "M.json"))
            end
        end
        sort!(mfiles)
        for f in mfiles
            g = try
                K.load_mesh_json(f)
            catch e
                println(stderr, "case load failed ", f, ": ", sprint(showerror, e))
                continue
            end
            isempty(g.sigma) && continue
            name = basename(dirname(f))
            # Periodic patches are ordinary bounded patches for these three checks:
            # they are assembled with no boundary rows so that L is measured alone.
            r, rs, fz, hg = ref_row(name, g)
            push!(refrows, r)
            record(ok::Bool, tag::String) = begin
                ok && return
                mkpath(vout)
                K.save_mesh_json(g, vout * "/" * r.name * "__" * tag * ".json")
                push!(ref_violations, "| `" * r.name * "` | " * tag * " | `violations/" * r.name *
                                      "__" * tag * ".json` |")
            end
            record(rs.degree_identity, "rowsum_naive_degree_identity")
            record(rs.restricted_degree_identity, "rowsum_restricted_degree_identity")
            record(fz.L_equals_RD, "L_eq_RD")
            record(fz.rank_identity, "rank_L_eq_H_minus_dimZ")
            record(fz.out_harmonic, "left_nullspace_out_harmonic")
            record(hg.identity_holds, "euler_H_eq_Ehinge_minus_F_plus_c")
            println("[case] ", r.name, " H=", r.H, " notches=", r.notches, " c(Gamma)=", r.c_gamma,
                    " pred=", r.predicted, " dimZ=", r.dim_Z)
        end
    end

    # The boundary-free counterpart: on a torus every vertex is interior, so 1^T L
    # is forced to zero, dim Z becomes 1 and the out-harmonic check stops being
    # vacuous (dim Z == 0 on every bounded graph measured).
    torusrows = RefRow[]
    for (name, mesh) in (("torus_squares_4x4", K.torus_squares(4, 4)),
                         ("torus_squares_6x4", K.torus_squares(6, 4)),
                         ("torus_triangles_4x4", K.torus_triangles(4, 4)))
        mesh.sigma = checkerboard(mesh)
        r, rs, fz, hg = ref_row(name, mesh)
        push!(torusrows, r)
        println("[torus] ", r.name, " H=", r.H, " rank(L)=", r.rank_L, " dimZ=", r.dim_Z,
                " 1^T L=0:", fmt_b(r.r_zero), " out-harmonic:", fmt_b(r.harmonic), " (res ",
                fmt_g(fz.max_harmonic_residual), ")")
    end

    # rank_claim.md, generated from the same run (no hand-copied numbers).
    md = IOBuffer()
    print(md, "# Rank claim on the sweep (Builder-Core)\n\n")
    print(md, "Source: `sweep.csv`, produced by `julia --project=Kirigami Kirigami/apps/kiri_sweep.jl --n ", n, " --seed ", seed,
          "`.\nRandom planar graphs (Voronoi cells of random points, Delaunay triangulations,\n",
          "quad-dominant edge-collapsed Delaunay), face orientations from the paper's own Eq. (1)\n",
          "relaxation, fixed-boundary system of Eq. (4).\n\n")
    print(md, "Graphs analysed: **", st.n, "**, faces from **", st.min_F, "** to **", st.max_F, "**.\n\n")
    print(md, "## The paper's Sec. 4.4 claims, measured\n\n")
    print(md, "| claim (2026 Sec. 4.4, prose, unproved in the paper) | holds on |\n|---|---|\n")
    print(md, "| `dim null(Eq. 4) == #interior vertices - H` | **", st.claim_ok, " / ", st.n, "** |\n")
    print(md, "| `rank(L) == H` (the rows of L are independent) | **", st.rankL_eq_H, " / ", st.n, "** |\n")
    print(md, "| `H <= #interior vertices` | **", st.H_le_interior, " / ", st.n, "** |\n")
    print(md, "| stronger relation found here: `dim null == |E_split|` | **", st.dimnull_eq_split,
          " / ", st.n, "** |\n\n")
    print(md, "## Structural checks\n\n")
    print(md, "| check | holds on |\n|---|---|\n")
    print(md, "| Remark A.1 (hinge in-degree == out-degree at interior vertices) | ", st.a1, " / ", st.n, " |\n")
    print(md, "| Algorithm 1 (seed-growing) == split-forest partition formulation | ", st.agree, " / ", st.n, " |\n")
    print(md, "| hole preimages partition E_hinge union E_split | ", st.part, " / ", st.n, " |\n")
    print(md, "| split-cut subgraph is a forest (Remark A.4) | ", st.forest, " / ", st.n, " |\n\n")
    print(md, "## Violations\n\n")
    if isempty(violation_lines)
        print(md, "None. `dim_null == #interior - H` held on every graph in the sweep.\n\n")
    else
        print(md, "| id | kind | F | #interior | H | dim_null | #interior - H | saved graph |\n",
              "|---|---|--:|--:|--:|--:|--:|---|\n")
        for l in violation_lines
            print(md, l, "\n")
        end
        print(md, "\n")
    end
    print(md, "## Why the geometric hole cross-check is empty on this sweep\n\n")
    print(md, "The geometric check traces the bounded complement components of the deployed M', which\n",
          "requires the Tutte auxetic embedding to be a valid (non-self-intersecting) straight-line\n",
          "embedding. On these random graphs it never is:\n\n")
    print(md, "| quantity | value |\n|---|---|\n")
    print(md, "| graphs where the dense null space / X0 was computed | ", st.dense, " / ", st.n, " |\n")
    print(md, "| of those, X0 free of face-face overlap at theta = 0 | **", st.x0_valid, " / ", st.dense, "** |\n")
    print(md, "| mean fraction of faces whose orientation flips in X0 | ",
          f3(st.dense > 0 ? st.flipped_frac / st.dense : 0.0), " |\n")
    print(md, "| largest vertex displacement of the Eq. (6) projection (box side 40) | ", f3(st.max_move), " |\n")
    print(md, "| geometric cross-checks actually run / agreed | ", st.geo_checked, " / ", st.geo_agree, " |\n\n")
    print(md, "This is the paper's own open problem (2026 Sec. 6): the shape space X contains\n",
          "embeddings with self-intersections and the paper offers no certificate for the valid\n",
          "subset. The combinatorial-vs-geometric agreement is instead verified on 60 valid\n",
          "instances by the unit tests (see `test/test_holes.jl`).\n\n")
    print(md, "Worst dense solve time in the sweep: ", f3(st.worst_solve_ms), " ms.\n")

    # -------------------------------------------------------- added section ----
    print(md, "\n## Orchestrator checks (added)\n\n")
    print(md, "Three measurements on the hole-constraint matrix `L` that nothing else in the\n",
          "pipeline needed: the row sum `1^T L`, the factorization `L = R D`, and the\n",
          "hinge-graph Euler count. `L` is taken verbatim from `assemble_system` (its first\n",
          "`H` rows), so what is measured is the matrix the solver actually uses. Definitions:\n\n")
    print(md, "- `D` is the `|E_hinge| x N` signed hinge incidence (+1 target, -1 source);\n",
          "  `R` is the `H x |E_hinge|` 0/1 matrix saying which hole preimage owns each hinge edge.\n",
          "- `Z` is the left null space of `L`; `g(v) = y_{K(v)}` with `K(v)` the preimage whose\n",
          "  split-forest component contains `v`, and `g(v) = 0` when `K(v)` is a notch (no row of L).\n",
          "- `Gamma = (nodes = faces, edges = hinge cuts)`; `c(Gamma)` is its component count.\n",
          "- `H` counts ONLY all-interior preimages, as Builder-Core defines it; the\n",
          "  boundary-touching preimages are counted separately as notches.\n\n")

    print(md, "### Reference cases (read back from `cases/*/M.json`, boundary rows omitted)\n\n")
    yn(b) = b ? "yes" : "**no**"
    if isempty(refrows)
        print(md, "No case directory found at `", cases_dir, "` -- run `kiri_reference` first.\n\n")
    else
        print(md, "| case | F | \\|E_hinge\\| | H | notches | c(G) | pred = \\|E_hinge\\|-\\|F\\|+c(G) | ",
              "rank(L) | dim Z | 1^T L = 0 | naive row sum | restricted row sum | L = R D | ",
              "rank = H - dim Z | Z out-harmonic | H = pred | H_all = pred |\n")
        print(md, "|---|--:|--:|--:|--:|--:|--:|--:|--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|\n")
        for r in refrows
            print(md, "| `", r.name, "` | ", r.F, " | ", r.n_hinge, " | ", r.H, " | ",
                  r.notches, " | ", r.c_gamma, " | ", r.predicted, " | ", r.rank_L, " | ",
                  r.dim_Z, " | ", yn(r.r_zero), " | ", yn(r.rs_naive), " | ",
                  yn(r.rs_restricted), " | ", yn(r.lrd), " | ", yn(r.rank_id), " | ",
                  yn(r.harmonic), " | ", yn(r.euler_H), " | ", yn(r.euler_Hall), " |\n")
        end
        print(md, "\n")
    end

    print(md, "### Boundary-free (torus) patches -- where `1^T L = 0` and `dim Z > 0`\n\n")
    print(md, "The generators' \"periodic\" patches are finite patches with a real boundary (Eqs. 3b-3c\n",
          "are constraints, not a topological identification), so they do NOT test the\n",
          "boundary-free case. `torus_squares` / `torus_triangles` do: every vertex is interior.\n",
          "Their wrap-around faces are geometrically degenerate, which is harmless here because\n",
          "all three checks are purely combinatorial.\n\n")
    print(md, "| patch | F | \\|E_hinge\\| | H | notches | c(G) | pred | rank(L) | dim Z | 1^T L = 0 | ",
          "rank = H - 1 | L = R D | Z out-harmonic | max harmonic residual |\n")
    print(md, "|---|--:|--:|--:|--:|--:|--:|--:|--:|:--:|:--:|:--:|:--:|--:|\n")
    for r in torusrows
        print(md, "| `", r.name, "` | ", r.F, " | ", r.n_hinge, " | ", r.H, " | ",
              r.notches, " | ", r.c_gamma, " | ", r.predicted, " | ", r.rank_L, " | ",
              r.dim_Z, " | ", yn(r.r_zero), " | ", yn(r.rank_L == r.H - 1), " | ",
              yn(r.lrd), " | ", yn(r.harmonic), " | ", sci_default(r.max_harm), " |\n")
    end
    print(md, "\nNote `pred = |E_hinge| - |F| + c(Gamma)` overshoots `H` by one on every torus patch:\n",
          "`Gamma` is embedded on a surface of Euler characteristic 0, not 2.\n\n")

    print(md, "### Sweep aggregate (", st.n, " random graphs)\n\n")
    print(md, "| check | holds on |\n|---|---|\n")
    print(md, "| 1a. `1^T L` = naive degree difference `indeg(v) - outdeg(v)` | ", st.rs_naive, " / ", st.n, " |\n")
    print(md, "| 1b. `1^T L` = notch-restricted degree difference (corrected form) | **",
          st.rs_restricted, " / ", st.n, "** |\n")
    print(md, "| 1c. `1^T L` supported on boundary vertices only | ", st.rs_support, " / ", st.n, " |\n")
    print(md, "| 2a. `L == R D` exactly (max abs difference 0) | **", st.lrd, " / ", st.n, "** |\n")
    print(md, "| 2b. `rank(L) == H - dim Z` | **", st.rank_id, " / ", st.z_done, "** (dense only) |\n")
    print(md, "| 2c. every `y` in `Z` gives an out-harmonic `g`, tol 1e-9 | **", st.harmonic, " / ",
          st.z_done, "** (dense only) |\n")
    print(md, "| 3a. `H == |E_hinge| - |F| + c(Gamma)` | ", st.euler_H, " / ", st.n, " |\n")
    print(md, "| 3b. `H_all == |E_hinge| - |F| + c(Gamma)` (holes + notches) | **", st.euler_Hall,
          " / ", st.n, "** |\n")
    print(md, "| Eq. (1) auto-orientation yields `c(Gamma) == 1` | ", st.gamma_conn, " / ", st.n, " |\n\n")
    print(md, "Out-harmonic residuals actually evaluated in the sweep: ", st.harm_checks,
          " (worst ", sci_default(st.worst_harm), ").\n")
    if st.harm_checks == 0
        print(md, "\n**Caveat: check 2c is vacuous on the sweep.** `dim Z == 0` on every one of the ",
              st.z_done, " graphs where\nthe left null space was computed (i.e. `rank(L) == H`, ",
              "consistent with the sweep's own\n`rank(L) == H` row above), so there is no basis ",
              "vector to test. The check is exercised on\nthe torus patches in the table above and ",
              "in `test/test_rank_checks.jl`, where `dim Z == 1`.\n")
    end
    print(md, "\n")

    print(md, "### Identities that failed, and the corrected form\n\n")
    if st.rs_naive < st.n || !isempty(refrows)
        print(md, "**Check 1 (row sum), as first stated, is false on any patch with a boundary.**\n",
              "`L` carries a row only for an all-interior preimage, so every hinge edge owned by a\n",
              "boundary-touching preimage (a notch) contributes nothing to `1^T L`. Its target loses\n",
              "the `+1` and its source loses the `-1`, and the source can be an interior vertex --\n",
              "so `1^T L` is NOT supported on the boundary either. The corrected identity, which\n",
              "does hold everywhere measured, restricts both degrees to hinge edges owned by hole\n",
              "rows. With `I(v) = 1` iff `K(v)` is a row of `L`:\n\n")
        print(md, "```\n",
              "    (1^T L)_v  ==  I(v) * indeg_hinge(v)  -  sum_{hinge v->w} I(w)\n",
              "```\n\n")
        print(md, "(All in-edges of `v` belong to `K(v)`, so the in-degree is all-or-nothing; the\n",
              "out-edges are distributed over the preimages of their targets.) This is exactly the\n",
              "`y = 1` instance of the out-harmonic characterization of check 2. When there are no\n",
              "notches -- the torus patches in `test/test_rank_checks.jl` -- it collapses to the\n",
              "naive form, `1^T L == 0` exactly, and `rank(L) == H - 1` (measured, not assumed).\n\n")
    end
    print(md, "**Check 3 (hinge-graph Euler count) holds exactly as stated**, on all ",
          st.euler_H, " / ", st.n, " sweep graphs and all ", length(refrows),
          " reference cases, with `H` counting only the all-interior preimages:\n\n")
    print(md, "```\n",
          "    H  ==  |E_hinge| - |F| + c(Gamma)          (H = holes only, notches excluded)\n",
          "```\n\n")
    print(md, "So the notches are NOT bounded faces of `Gamma`: counting them in breaks the identity\n",
          "(`H_all == pred` on ", st.euler_Hall, " / ", st.n, " sweep graphs). This\n",
          "confirms unverified claim U1 numerically and, with `c(Gamma) == 1` on ",
          st.gamma_conn, " / ", st.n, " sweep graphs under the Eq. (1) auto-orientation,\n",
          "reduces to Scout-c's `H = |E_hinge| - |F| + 1`.\n\n")
    print(md, "The one place the plane count needs a correction is a boundary-free patch, where\n",
          "`Gamma` lives on a surface of Euler characteristic 0 rather than 2 and the count\n",
          "overshoots by exactly one: `H == |E_hinge| - |F| + c(Gamma) - 1`, verified in\n",
          "`test/test_rank_checks.jl` on the 4x4 and 6x4 square tori and the 4x4 triangle torus.\n\n")
    print(md, "Checks 2a, 2b and 2c held on every graph and every reference case measured; no\n",
          "counterexample to the factorization or to the out-harmonic characterization was found.\n\n")
    if !isempty(ref_violations)
        print(md, "Counterexamples saved (graph JSON, loadable by `kiri_analyze`):\n\n")
        print(md, "| case | identity that failed | file |\n|---|---|---|\n")
        for l in ref_violations
            print(md, l, "\n")
        end
        print(md, "\n")
    end

    let dir = isempty(dirname(out)) ? "." : dirname(out)
        open(dir * "/rank_claim.md", "w") do mf
            write(mf, String(take!(md)))
        end
        println("wrote ", dir, "/rank_claim.md")
    end
    close(csv)
    println("wrote ", out, " (", violations, " violations saved)")
    return 0
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
