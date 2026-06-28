# test_reference_cases.jl -- port of code/tests/test_reference_cases.cpp, case by case.
#
# Each of the eight Phase-2 reference cases is regenerated with the SAME generator calls,
# in the SAME order, off the SAME seeded rng as kiri_reference.cpp (`Kirigami.reference_cases`
# = kill_common.hpp::reference_cases, byte-for-byte the C++ test's build_reference_cases),
# then run through make_cut / holes_partition / assemble_system / solve_system and
# compared against the stored results/core_validation/reference_cases.json (frozen copy in
# CORPUS/reference_patterns/).
#
# WHAT IS PINNED AGAINST THE JSON, AND WHAT IS NOT (from the C++ test):
#   * F, n_split and dim_null are combinatorial and are asserted exactly against the JSON.
#   * The JSON's `theta_max` column is a Gate-2 artefact of the pre-F34 collision code plus
#     bisection; it disagrees with today's exact T4.2'' value by ~2e-6 on the split
#     patterns, so it is cross-checked at 3e-6 only.
#   * The exact range is asserted against `characterize` (T4.2''): with no split cuts,
#     Theta_max == min(min_e beta_e, pi) exactly (seven of eight cases, including the three
#     that have split cuts but are hinge-bound). snub_square_33434 is the one split-bound
#     case and is pinned to its measured T4.2'' value.
#
# Julia note: on an x86_64 (Rosetta) Julia the relaxation cases hexagons_auto and
# truncated_square_488 can round to an orientation that differs from the frozen one by a
# symmetry (see generators.jl on libm); the quantities checked here are invariant to that.
include("helpers.jl")
import JSON

const K = Kirigami

const REFERENCE_JSON = joinpath(CORPUS, "reference_patterns", "reference_cases.json")

# The one split-BOUND reference case: its first contact is a split duplicate, so
# Theta_max sits far below min beta (3.3839) and no hinge-local bound predicts it.
# Value measured with the C++ exact T4.2'' code (method::characterize).
const kSnubThetaMaxT4 = 1.646135970500460

@testset "reference cases: reference_cases.json is present and complete" begin
    @test isfile(REFERENCE_JSON)
    j = JSON.parsefile(REFERENCE_JSON)
    @test j isa AbstractVector
    @test length(j) == 8
    for row in j
        for k in ("name", "F", "split", "dim_null", "theta_max", "min_beta")
            @test haskey(row, k)
        end
    end
end

@testset "reference cases: F, n_split, dim_null and the exact T4 range reproduce" begin
    j = JSON.parsefile(REFERENCE_JSON)
    @test length(j) == 8
    cases = K.reference_cases()
    @test length(cases) == length(j)

    have_system = isdefined(K, :assemble_system) && isdefined(K, :solve_system) &&
                  isdefined(K, :matrix_to_points)
    have_characterize = isdefined(K, :characterize)
    have_system || @warn "assemble_system/solve_system not available: dim_null and theta_max checks skipped"
    have_characterize || @warn "characterize not available: theta_max checks skipped"

    for (cs, row) in zip(cases, j)
        @testset "$(cs.name)" begin
            @test row["name"] == cs.name
            g = cs.mesh
            K.build_topology!(g)
            c = K.make_cut(g)
            hs = K.holes_partition(c)

            # --- combinatorics: exact against the stored artefact ---------------------
            @test K.n_faces(g) == row["F"]
            @test K.n_split(c) == row["split"]

            X0 = K.Vec2[]
            if have_system
                sys = K.assemble_system(c, hs, g.X, cs.periodic ? K.Periodic : K.Fixed)
                rep = K.solve_system(sys, g.X)
                X0 = K.matrix_to_points(rep.X0)
                @test rep.dim_null == row["dim_null"]
            end

            if have_system && have_characterize
                # --- the exact range, T4.2'' -----------------------------------------
                ch = K.characterize(g, g.sigma, X0)
                @test ch.n_split == row["split"]
                @test !ch.zero_range

                min_beta = Float64(row["min_beta"])
                hinge_bound = min(min_beta, pi)
                if cs.name == "snub_square_33434"
                    # Split-bound: the first contact is a split duplicate, well below min beta.
                    @test ch.theta_max < hinge_bound - 1.0
                    @test abs(ch.theta_max - kSnubThetaMaxT4) < 1e-9
                else
                    # Hinge-bound: T4.1 says Theta_max == min(min_e beta_e, pi) exactly here.
                    @test abs(ch.theta_max - hinge_bound) < 1e-9
                end

                # --- cross-check against the stored (pre-F34) theta_max column --------
                @test abs(ch.theta_max - Float64(row["theta_max"])) < 3e-6
            end
        end
    end
end
