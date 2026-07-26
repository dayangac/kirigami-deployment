#include <cstdio>
#include <filesystem>
#include "kill_common.hpp"
#include "method/zero_plus.hpp"
using namespace kiri; using namespace kiri::kill; using namespace kiri::method;
static void dump(const char* tag, const ZeroPlusRepairResult& r) {
  std::printf("%s: feasible %d best_start %d iters %d f_start %.17g f_end %.17g min_q_start %.17g n_bad_q_start %d n_bad_area_start %d "
              "min_q %.17g min_area %.17g n_bad_q %d n_bad_area %d min_margin_start %.17g n_bad_margin_start %d min_margin %.17g n_bad_margin %d n_inc %d |t| %.17g\n",
              tag, r.feasible, r.best_start, r.iterations, r.f_start, r.f_end, r.min_q_start, r.n_bad_q_start, r.n_bad_area_start,
              r.min_q, r.min_area, r.n_bad_q, r.n_bad_area, r.min_margin_start, r.n_bad_margin_start, r.min_margin, r.n_bad_margin, r.n_incidences, r.t.norm());
}
int main(int argc, char** argv) {
  const int id = argc > 1 ? std::atoi(argv[1]) : 93;
  const int max_iter = argc > 2 ? std::atoi(argv[2]) : 1200;
  Graph g = make_graph(id, 100, 800, 1400);
  Mesh m = g.mesh; m.build_topology();
  const Mesh sm = load_mesh_json("results/kill/k5/sigma/" + g.kind + "_" + std::to_string(id) + ".json");
  m.sigma = sm.sigma; m.build_topology();
  const CutStructure c = make_cut(m);
  Shape sh; load_shape("results/kill/b4/cache/free_sigma_def", id, &sh);
  const std::vector<Vec2> X0 = matrix_to_points(sh.X0);
  const double med = median_edge_length(m);
  std::printf("N %d k %d med %.17g n_split %d\n", sh.N, sh.k, med, c.n_split());
  const unsigned seed0 = 6000u + 7u * id + 1u;
  ZeroPlusRepairOptions ro; ro.n_random = 3; ro.max_iter = max_iter; ro.lambda_rel = 1e-6; ro.seed = seed0;
  const ZeroPlusRepairResult sr = zero_plus_repair(c, X0, sh.Phi, med, ro);
  dump("secondary", sr);
  ZeroPlusRepairOptions ro2 = ro; ro2.w_corner = 0.05; ro2.w_prox = 1e3; if (sr.feasible) ro2.t_init = sr.t;
  const ZeroPlusRepairResult rr = zero_plus_repair(c, X0, sh.Phi, med, ro2);
  dump("primary", rr);
}
