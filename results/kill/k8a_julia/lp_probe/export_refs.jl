using Kirigami, JSON; const K = Kirigami
include("/Users/emredayangac/Documents/kirigami-julia/Kirigami/apps/kill_common.jl")
out = []
for rc in reference_cases()
    rc.name in ("snub_square_33434", "truncated_square_488", "hexagons_auto", "tiling_3_4_3_12") || continue
    m = rc.mesh; K.build_topology!(m); c = K.make_cut(m); hs = K.holes_partition(c); X = copy(m.X)
    if !K.deployable(K.hole_residuals(c, X, hs), 1e-9)
        sr = K.solve_system(K.assemble_system(c, hs, m.X, K.Fixed), m.X); X = K.matrix_to_points(sr.X0)
    end
    push!(out, Dict("id" => rc.name, "kind" => rc.name, "mesh" => K.mesh_to_json(m), "X" => [[p[1], p[2]] for p in X]))
end
open(joinpath(ARGS[1], "ref_inputs.json"), "w") do io; JSON.print(io, out); end
