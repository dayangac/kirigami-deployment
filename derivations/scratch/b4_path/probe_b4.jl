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
function dump(tag, r)
    @printf("%s: feasible %d best_start %d iters %d f_start %.17g f_end %.17g min_q_start %.17g n_bad_q_start %d n_bad_area_start %d min_q %.17g min_area %.17g n_bad_q %d n_bad_area %d min_margin_start %.17g n_bad_margin_start %d min_margin %.17g n_bad_margin %d n_inc %d |t| %.17g\n",
        tag, r.feasible, r.best_start, r.iterations, r.f_start, r.f_end, r.min_q_start, r.n_bad_q_start, r.n_bad_area_start, r.min_q, r.min_area, r.n_bad_q, r.n_bad_area, r.min_margin_start, r.n_bad_margin_start, r.min_margin, r.n_bad_margin, r.n_incidences, norm(r.t))
end
id = parse(Int, ARGS[1]); max_iter = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 1200
g = make_graph(id, 100, 800, 1400); m = g.mesh; K.build_topology!(m)
sm = K.load_mesh_json(joinpath(REPO, "results/kill/k5/sigma", g.kind * "_$id.json")); m.sigma = sm.sigma; K.build_topology!(m)
c = K.make_cut(m)
X0m, Phi = load_shape(expanduser("~/Documents/kirigami-experiments/results/kill/b4/cache/free_sigma_def"), id)
X0 = K.matrix_to_points(X0m); med = K.median_edge_length(m)
@printf("N %d k %d med %.17g n_split %d\n", size(Phi,1), size(Phi,2), med, K.n_split(c))
seed0 = 6000 + 7id + 1
ro = K.ZeroPlusRepairOptions(n_random=3, max_iter=max_iter, lambda_rel=1e-6, seed=seed0)
sr = K.zero_plus_repair(c, X0, Phi, med, ro); dump("secondary", sr)
ro2 = K.ZeroPlusRepairOptions(n_random=3, max_iter=max_iter, lambda_rel=1e-6, seed=seed0, w_corner=0.05, w_prox=1e3, t_init = sr.feasible ? sr.t : Float64[])
rr = K.zero_plus_repair(c, X0, Phi, med, ro2); dump("primary", rr)
