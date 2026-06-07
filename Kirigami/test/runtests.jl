using Test
using Kirigami

const TEST_FILES = [
    "test_mesh_cut.jl", "test_holes.jl", "test_system.jl", "test_kinematics.jl",
    "test_collision.jl", "test_reference_cases.jl", "test_rank_checks.jl", "test_mt19937.jl",
    "test_import_soup.jl", "test_method.jl", "test_design.jl", "test_range_embed.jl",
    "test_export.jl", "derivation_tests.jl",
]

@testset "Kirigami" begin
    for f in TEST_FILES
        p = joinpath(@__DIR__, f)
        isfile(p) || (@warn "missing test file" f; continue)
        @testset "$f" begin
            include(p)
        end
    end
end
