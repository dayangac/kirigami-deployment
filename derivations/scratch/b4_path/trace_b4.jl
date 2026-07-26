using Kirigami, LinearAlgebra, Printf
const K = Kirigami
include("/Users/emredayangac/Documents/kirigami-julia/derivations/scratch/corpus_common.jl")
const REPO = "/Users/emredayangac/Documents/kirigami-julia"
function load_shape(dir, id)
    open(joinpath(dir, "shape_$id.bin"), "r") do io
        hdr = [read(io, Int32) for _ in 1:5]; N = Int(hdr[1]); k = Int(hdr[2])
        X0 = Matrix{Float64}(undef, N, 2); read!(io, X0)
        Phi = Matrix{Float64}(undef, N, k); read!(io, Phi)
        return X0, Phi
    end
end
id = parse(Int, ARGS[1])
g = make_graph(id, 100, 800, 1400); m = g.mesh; K.build_topology!(m)
sm = K.load_mesh_json(joinpath(REPO, "results/kill/k5/sigma", g.kind * "_$id.json")); m.sigma = sm.sigma; K.build_topology!(m)
c = K.make_cut(m)
X0m, Phi = load_shape(expanduser("~/Documents/kirigami-experiments/results/kill/b4/cache/free_sigma_def"), id)
X0 = K.matrix_to_points(X0m); med = K.median_edge_length(m)
seed0 = 6000 + 7id + 1
for it in (1, 2, 3, 5, 10, 20, 40, 80, 160, 320, 640)
    r = K.zero_plus_repair(c, X0, Phi, med, K.ZeroPlusRepairOptions(n_random=3, max_iter=it, lambda_rel=1e-6, seed=seed0))
    @printf("Jul iters %d best_start %d f_end %.17g min_q %.17g |t| %.17g\n", it, r.best_start, r.f_end, r.min_q, norm(r.t))
end
