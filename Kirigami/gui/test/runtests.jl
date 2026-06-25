# gui/test/runtests.jl -- headless tests of gui/model.jl (no window, no GLMakie).
using Test
import JSON
include(joinpath(@__DIR__, "..", "model.jl"))
using .GuiModel: GuiModel, Design, deploy_basis, face_sense, frame, design_from_json, load_designs,
    deploy_state, closed, certified, usable, over, badge_text, areal_expansion,
    track_fractions, view_box, ease_theta, play_range, svg_string, state_dict, state_json,
    FileBackend, designs

const DATA = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "deploy_frames.json"))

# unit square that expands along x with a hinge-like path; a 2-face synthetic design
function synthetic(; theta_half = 1.0, theta_max = 2.0, eps_max = 1.8)
    C = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 2.0 0.0; 2.0 1.0]
    S = [0.0 0.0; 1.5 0.2; 1.5 1.2; 0.0 1.0; 3.0 0.0; 3.0 1.0]
    v1 = cos(theta_half / 2) .* C .+ sin(theta_half / 2) .* S
    j = Dict("id" => "synthetic", "kind" => "test", "sigma" => "def", "F" => 2, "nsplit" => 1,
             "theta_max" => theta_max, "eps_max" => eps_max, "theta_half" => theta_half,
             "faces" => [[0, 1, 2, 3], [1, 4, 5, 2]],
             "v0" => [C[i, :] for i in 1:6], "v1" => [v1[i, :] for i in 1:6])
    return j, C, S
end

@testset "deploy_basis / frame" begin
    j, C, S = synthetic()
    d = design_from_json(j)
    @test d.faces == [[1, 2, 3, 4], [2, 5, 6, 3]]          # 0-based file -> 1-based
    @test d.sigma == "def"
    @test isapprox(d.C, C; atol = 1e-12)
    @test isapprox(d.S, S; atol = 1e-9)
    @test isapprox(frame(d, 0.0), C; atol = 1e-12)
    v1 = reduce(vcat, (Float64.(p)' for p in j["v1"]))
    @test isapprox(frame(d, d.theta_half), v1; atol = 1e-9)   # Y(θ_half) == v1
    @test_throws ArgumentError deploy_basis(C, C, 0.0)          # sin(0) == 0
    @test_throws ArgumentError deploy_basis(C, C[1:2, :], 1.0)
    # the basis is independent of which θ_half was used to sample v1
    j2, _, _ = synthetic(theta_half = 2.5)
    @test isapprox(design_from_json(j2).S, S; atol = 1e-9)
end

@testset "face_sense" begin
    # two unit squares, one rotating +θ/2 (S = rot90 C about its centre), one −θ/2
    sq = [-0.5 -0.5; 0.5 -0.5; 0.5 0.5; -0.5 0.5]
    rot90(P) = hcat(-P[:, 2], P[:, 1])
    C = vcat(sq, sq .+ [2.0 0.0])
    S = vcat(rot90(sq), -rot90(sq) .+ [2.0 0.0])
    faces = [[1, 2, 3, 4], [5, 6, 7, 8]]
    @test face_sense(C, S, faces) == [1, -1]
    # rounding noise of the file (4 decimals) must not flip the majority
    @test face_sense(C, S .+ 1e-4 .* [0.3 -0.7; 0.9 0.1; -0.2 0.5; 0.4 0.4; 0.1 0.6; -0.8 0.2; 0.7 -0.3; 0.5 0.5], faces) == [1, -1]
    @test face_sense(C, S, Vector{Int}[]) == Int[]
end

@testset "state / readouts" begin
    d = design_from_json(synthetic()[1])
    @test deploy_state(d, 0.0) == closed
    @test deploy_state(d, 1.0) == certified
    @test deploy_state(d, 1.8) == certified                  # boundary inclusive
    @test deploy_state(d, 1.9) == usable
    @test deploy_state(d, 2.0) == usable
    @test deploy_state(d, 2.0 + 1e-6) == over
    @test badge_text(over) == "past Θmax — collisions"
    @test areal_expansion(d, 0.0) ≈ 1.0
    @test areal_expansion(d, 1.0) > 1.0
    @test track_fractions(d) == (2.0 / pi, 1.8 / pi)
    dd = design_from_json(synthetic(theta_max = 4.0, eps_max = 3.5)[1])
    @test track_fractions(dd) == (1.0, 1.0)
    vb = view_box(d)
    @test vb[1] <= 0.0 && vb[2] <= 0.0 && vb[3] >= 2.0 && vb[4] >= 1.0
    @test ease_theta(0.0, 2.0, 0.0) == 0.0
    @test ease_theta(0.0, 2.0, 1.0) ≈ 2.0
    @test ease_theta(0.0, 2.0, 0.5) ≈ 1.0
    @test ease_theta(0.0, 2.0, 7.0) ≈ 2.0                    # clamped
    @test play_range(d, 0.5) == (0.5, 2.0)
    @test play_range(d, 2.0) == (0.0, 2.0)                   # at the end: restart
    d0 = design_from_json(synthetic(theta_max = 0.0, eps_max = 0.0)[1])
    @test play_range(d0, 0.0) == (0.0, 0.05)
end

@testset "exports" begin
    d = design_from_json(synthetic()[1])
    s = svg_string(d, 1.0)
    @test startswith(s, "<svg")
    @test count("<polygon", s) == 2
    @test occursin(GuiModel.FACE_A, s) && occursin(GuiModel.FACE_B, s)
    st = state_dict(d, 1.0)
    @test st["faces"] == [[0, 1, 2, 3], [1, 4, 5, 2]]        # written back 0-based
    @test st["sense"] == d.sense
    @test st["state"] == "certified"
    @test length(st["vertices"]) == 6
    back = JSON.parse(state_json(d, 1.0))
    @test back["theta"] == 1.0 && back["id"] == "synthetic"
end

@testset "load_designs" begin
    path = tempname() * ".json"
    j, _, _ = synthetic()
    open(path, "w") do io
        JSON.print(io, [j, Dict("name" => "raw", "text" => "v 0 0\n")])
    end
    ds, skipped = load_designs(path)
    @test length(ds) == 1 && skipped == 1
    @test ds[1].id == "synthetic"
    ds2, _ = designs(FileBackend(path))
    @test length(ds2) == 1
    rm(path)
    if isfile(DATA)
        ds, skipped = load_designs(DATA)
        @info "data/deploy_frames.json" designs = length(ds) skipped
        @test length(ds) == 26 && skipped == 0
        for d in ds
            @test all(all(1 .<= f .<= size(d.C, 1)) for f in d.faces)
            @test 0 <= d.eps_max <= d.theta_max + 1e-12
            @test length(d.sense) == d.F
            @test 0 < count(>(0), d.sense) < d.F               # both senses present
            # rigidity of the frame path: every face edge keeps its length at θ_half
            # (file coordinates are rounded to 4 decimals -> 3e-4 absolute / 1e-3 relative)
            P = frame(d, d.theta_half)
            for f in d.faces, k in eachindex(f)
                a, b = f[k], f[mod1(k + 1, length(f))]
                @test isapprox(hypot(P[b, 1] - P[a, 1], P[b, 2] - P[a, 2]),
                               hypot(d.C[b, 1] - d.C[a, 1], d.C[b, 2] - d.C[a, 2]); atol = 3e-4, rtol = 1e-3)
            end
        end
    end
end
