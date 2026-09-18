# kiri_reference --out <dir> [--regenerate]
# Runs every Phase-2 reference case and writes results/core_validation artifacts.
#
# The eight cases are read from the frozen
# data/corpus/reference_cases_8.json by default; --regenerate rebuilds them through the
# tiling generators and the bit-exact MT19937 (kill_common.jl `reference_cases`).
include(joinpath(@__DIR__, "kill_common.jl"))

Base.@kwdef mutable struct RefRow
    name::String = ""
    note::String = ""
    N::Int = 0
    F::Int = 0
    interior::Int = 0
    hinge::Int = 0
    split::Int = 0
    H::Int = 0
    H_geo::Int = 0
    comps::Int = 0
    rank_L::Int = 0
    dim_null::Int = 0
    claim::Int = 0
    claim_ok::Bool = false
    claim_applicable::Bool = true
    a1::Bool = false
    holes_agree::Bool = false
    geo_agree::Bool = false
    forest::Bool = false
    partition::Bool = false
    deployable::Bool = false
    res_before::Float64 = 0.0
    res_after::Float64 = 0.0
    theta_max::Float64 = 0.0
    min_beta::Float64 = 0.0
    theta_max_opt::Float64 = -1.0
    gamma_used::Float64 = -1.0
end

const CASE_NOTES = Dict(
    "rotating_squares" => "5x5 square grid, checkerboard sigma",
    "triangles_alternating" => "equilateral triangle tiling, alternating sigma",
    "kagome_3636" => "trihexagonal 3.6.3.6, alternating sigma",
    "hexagons_auto" => "hexagonal tiling, sigma from Eq. (1); dual not 2-colourable",
    "truncated_square_488" => "4.8.8 truncated square, sigma from Eq. (1)",
    "snub_square_33434" => "3.3.4.3.4 snub square, sigma from Eq. (1)",
    "tiling_3_4_3_12" => "2-uniform [3.4.3.12; 3.12.12] (2026 Fig. 21): NOT deployable as the " *
                         "Euclidean tiling",
    "periodic_squares_4x4" => "4x4 periodic square tiling, Eq. (3b)-(3c)",
)

f6(v) = fx(v, 6)   # std::fixed << std::setprecision(6)

function main(args_in::Vector{String})
    args, regenerate = take_regenerate_flag(args_in)
    outdir = "../results/core_validation"
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--out" && i + 1 <= length(args)
            outdir = args[i += 1]
        end
        i += 1
    end
    mkpath(outdir * "/cases")

    cases = reference_cases(regenerate = regenerate)

    rows = RefRow[]
    for cs in cases
        g = cs.mesh
        r = RefRow(name = cs.name, note = get(CASE_NOTES, cs.name, ""))
        c = K.make_cut(g)
        hs = K.holes_partition(c)
        hs2 = K.holes_seed_growing(c)
        r.N = K.n_vertices(g)
        r.F = K.n_faces(g)
        r.interior = K.n_interior_vertices(g)
        r.hinge = K.n_hinge(c)
        r.split = K.n_split(c)
        r.H = K.n_interior_holes(hs)
        r.comps = K.count_components(c)
        r.a1 = K.check_remark_A1(c)[1]
        r.holes_agree = K.same_hole_sets(hs, hs2)[1]
        r.partition = K.holes_partition_edges(c, hs)[1]
        r.forest = K.split_subgraph_is_forest(c)[1]

        rb = K.hole_residuals(c, g.X, hs)
        r.res_before = rb.max_norm
        r.deployable = K.deployable(rb)
        sys = K.assemble_system(c, hs, g.X, cs.periodic ? K.Periodic : K.Fixed)
        rep = K.solve_system(sys, g.X)
        r.rank_L = rep.rank_L
        r.dim_null = rep.dim_null
        r.claim = r.interior - r.H
        r.claim_ok = (r.dim_null == r.claim)
        # The paper states the rank claim for the fixed-boundary system ("the number of
        # unknowns equals the number of interior vertices"). Under periodic boundary
        # conditions the boundary vertices are not pinned, so the comparison does not apply.
        r.claim_applicable = !cs.periodic
        X0 = K.matrix_to_points(rep.X0)
        r.res_after = K.hole_residuals(c, X0, hs).max_norm

        tm = K.theta_max(c, X0, 180, 45)
        r.theta_max = tm.theta_max_geometric
        r.min_beta = tm.min_beta

        # geometric hole cross-check on the deployable embedding
        th = min(0.25, 0.4 * max(r.theta_max, 1e-3))
        Yd = K.deploy(c, X0, th).Y
        k = 0
        while k < 6 && K.has_collision(c, Yd)
            th *= 0.4
            Yd = K.deploy(c, X0, th).Y
            k += 1
        end
        geo = K.holes_geometric(c, Yd)
        r.H_geo = length(geo)
        comb = sort([hs.all[i].edges for i in hs.interior_indices])
        r.geo_agree = (geo == comb)

        if rep.dim_null > 0 && K.n_split(c) > 0
            cr = K.optimize_collision_sweep(c, X0, rep.Phi)
            r.theta_max_opt = cr.theta_max_after
            r.gamma_used = cr.gamma_used
        end

        dir = outdir * "/cases/" * cs.name
        mkpath(dir)
        m0 = K.Mesh(X0, g.faces)
        m0.sigma = g.sigma
        write_json(mesh_json(g), dir * "/M.json")
        write_json(mesh_json(m0), dir * "/X0.json")
        write_json(deployment_json(c, X0, 30.0 * pi / 180.0), dir * "/deploy_30deg.json")
        write_json(deployment_json(c, X0, 60.0 * pi / 180.0), dir * "/deploy_60deg.json")
        write_json(deployment_json(c, X0, max(1e-3, 0.999 * r.theta_max)), dir * "/deploy_thetamax.json")
        push!(rows, r)
        println(cs.name, ": F=", r.F, " H=", r.H, " dim_null=", r.dim_null, " claim=", r.claim,
                (!r.claim_applicable ? " (n/a: periodic)" : (r.claim_ok ? " ok" : " VIOLATION")),
                " theta_max=", fmt_g(r.theta_max))
    end

    md = IOBuffer()
    print(md, "# Reference cases (Builder-Core, Phase 2 gate)\n\n")
    print(md, "Generated by `kiri_reference`. All numbers measured, none copied from the paper.\n\n")
    print(md, "Conventions: sigma = +1 clockwise, sigma = -1 counter-clockwise; a hinge cut keeps the\n",
          "hinge at the SOURCE vertex and duplicates the TARGET (2026 Sec. 3).\n",
          "`H` counts hole preimages all of whose vertices are interior; `H_geo` counts the bounded\n",
          "complement components of the deployed M' traced geometrically (Definition 4.1).\n\n")
    print(md, "| case | N | F | interior | hinge | split | H | H_geo | comps | rank(L) | dim_null | ",
          "interior-H | claim | residual before | residual after | theta_max | min beta |\n")
    print(md, "|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--:|--:|--:|--:|--:|\n")
    for r in rows
        print(md, "| ", r.name, " | ", r.N, " | ", r.F, " | ", r.interior, " | ", r.hinge,
              " | ", r.split, " | ", r.H, " | ", r.H_geo, " | ", r.comps, " | ",
              r.rank_L, " | ", r.dim_null, " | ", r.claim, " | ",
              (!r.claim_applicable ? "n/a" : (r.claim_ok ? "OK" : "**FAIL**")), " | ",
              f6(r.res_before), " | ", f6(r.res_after), " | ", f6(r.theta_max), " | ",
              f6(r.min_beta), " |\n")
    end
    print(md, "\n## Checks\n\n")
    print(md, "| case | Remark A.1 | Alg.1 == partition | preimages partition E | split forest | ",
          "geometric == combinatorial | deployable as given | theta_max after Eq. (9) | gamma |\n")
    print(md, "|---|:--:|:--:|:--:|:--:|:--:|:--:|--:|--:|\n")
    for r in rows
        print(md, "| ", r.name, " | ", (r.a1 ? "OK" : "FAIL"), " | ",
              (r.holes_agree ? "OK" : "FAIL"), " | ", (r.partition ? "OK" : "FAIL"), " | ",
              (r.forest ? "OK" : "FAIL"), " | ", (r.geo_agree ? "OK" : "FAIL"), " | ",
              (r.deployable ? "yes" : "no"), " | ")
        if r.theta_max_opt >= 0
            print(md, f6(r.theta_max_opt), " | ", f6(r.gamma_used), " |\n")
        else
            print(md, "n/a | n/a |\n")
        end
    end
    print(md, "\n## Notes and figures\n\n")
    print(md, "Each figure shows M coloured by sigma, the Eq. (6) projection X0, and M' deployed at\n",
          "30 deg, 60 deg and theta_max with the holes shaded.\n\n")
    for r in rows
        print(md, "- **", r.name, "**: ", r.note, "  \n  ![", r.name, "](cases/", r.name, ".png)\n")
    end
    open(outdir * "/reference_cases.md", "w") do f
        write(f, String(take!(md)))
    end

    jrows = [Dict{String,Any}("name" => r.name, "N" => r.N, "F" => r.F, "interior" => r.interior,
                              "hinge" => r.hinge, "split" => r.split, "H" => r.H, "H_geo" => r.H_geo,
                              "rank_L" => r.rank_L, "dim_null" => r.dim_null, "claim" => r.claim,
                              "claim_ok" => r.claim_ok, "claim_applicable" => r.claim_applicable,
                              "res_before" => r.res_before, "res_after" => r.res_after,
                              "theta_max" => r.theta_max, "min_beta" => r.min_beta,
                              "theta_max_opt" => r.theta_max_opt) for r in rows]
    write_json(jrows, outdir * "/reference_cases.json")
    println("wrote ", outdir, "/reference_cases.md")
    return 0
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
