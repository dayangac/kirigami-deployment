// freeze_design.cpp -- intermediate values of method/design.cpp for the Julia port
// (data/corpus/method_fixtures/design_intermediates.json). Built against
// code/build/libkiri_core.a, same toolchain as the corpus (Apple clang, -O2, arm64).
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include "kill_common.hpp"
#include "method/design.hpp"
#include "method/convex_embed.hpp"
#include "method/zero_plus.hpp"

using namespace kiri;
using namespace kiri::kill;
using namespace kiri::method;
using json = nlohmann::json;

static json pts(const std::vector<Vec2>& X) {
  json a = json::array();
  for (const auto& p : X) a.push_back({p.x(), p.y()});
  return a;
}
static json vec(const Eigen::VectorXd& v) {
  json a = json::array();
  for (int i = 0; i < v.size(); ++i) a.push_back(v[i]);
  return a;
}
static json mat(const Eigen::MatrixXd& M) {  // row-major nested
  json a = json::array();
  for (int i = 0; i < M.rows(); ++i) {
    json r = json::array();
    for (int j = 0; j < M.cols(); ++j) r.push_back(M(i, j));
    a.push_back(r);
  }
  return a;
}
static json ch_json(const Characterization& ch) {
  json j;
  j["theta_max"] = ch.theta_max;
  j["eps_max"] = ch.eps_max;
  j["zero_range"] = ch.zero_range;
  j["contacts"] = ch.contacts;
  j["i_star"] = ch.i_star;
  j["first_contact_theta"] = ch.first_contact.theta;
  j["certified"] = ch.certified;
  j["first_root_at_pi"] = -1;
  j["binding"] = ch.binding;
  j["n_pairs"] = ch.n_pairs;
  return j;
}

int main(int argc, char** argv) {
  const int id = argc > 1 ? std::atoi(argv[1]) : 148;
  const std::string out = argc > 2 ? argv[2] : "design_intermediates_148.json";
  Graph g = make_graph(id, 100, 800, 1400);
  Mesh m = g.mesh;
  m.build_topology();
  const std::vector<Vec2> X_ini = m.X;
  DesignOptions opt;
  opt.seed = 9000u + 7u * static_cast<unsigned>(id);

  json j;
  j["provenance"] = "freeze_design.cpp (scratchpad), make_graph(id,100,800,1400), K9 seed 9000+7*id, "
                    "DesignOptions defaults; clang -O2 arm64 against code/build/libkiri_core.a";
  j["id"] = id;
  j["seed"] = opt.seed;

  const CutStructure c = make_cut(m);
  const double med = median_edge_length(m);
  j["med_edge"] = med;
  const HoleSet hs = holes_partition(c);
  const LinearSystem sys = assemble_system(c, hs, m.X, BoundaryMode::Fixed);
  const SolveReport sr = solve_system(sys, m.X);
  const std::vector<Vec2> X0 = matrix_to_points(sr.X0);
  j["dim_null"] = sr.dim_null;
  j["X0"] = pts(X0);
  j["Phi"] = mat(sr.Phi);

  // the three convex_embed / zero_plus_repair calls of design_constrained, in order
  ConvexEmbedOptions oa;
  oa.delta_rel = opt.delta_convex; oa.split_delta_rel = 0.0; oa.n_random = opt.restarts;
  oa.seed = opt.seed; oa.max_iter = opt.max_iter; oa.barrier_stages = opt.barrier_stages;
  const ConvexEmbedResult ra = convex_embed(c, X0, sr.Phi, X_ini, med, oa);
  j["ra"] = {{"t", vec(ra.t)}, {"X", pts(ra.X)}, {"feasible", ra.feasible}, {"n_bad", ra.n_bad},
             {"n_bad_q", ra.n_bad_q}, {"iterations", ra.iterations}, {"best_start", ra.best_start},
             {"min_cross", ra.min_cross}, {"barrier_stages_kept", ra.barrier_stages_kept}};
  ZeroPlusRepairOptions ro;
  ro.n_random = opt.restarts; ro.max_iter = opt.max_iter; ro.start_scale = opt.start_scale; ro.seed = opt.seed;
  const ZeroPlusRepairResult sp = zero_plus_repair(c, X0, sr.Phi, med, ro);
  j["sp"] = {{"t", vec(sp.t)}, {"feasible", sp.feasible}, {"min_q", sp.min_q}, {"iterations", sp.iterations},
             {"best_start", sp.best_start}};
  ConvexEmbedOptions ob = oa;
  ob.split_delta_rel = opt.delta_split;
  const ConvexEmbedResult rb = convex_embed(c, X0, sr.Phi, X_ini, med, ob);
  j["rb"] = {{"t", vec(rb.t)}, {"X", pts(rb.X)}, {"feasible", rb.feasible}, {"n_bad", rb.n_bad},
             {"n_bad_q", rb.n_bad_q}, {"iterations", rb.iterations}, {"best_start", rb.best_start},
             {"min_cross", rb.min_cross}, {"min_q", rb.min_q}, {"barrier_stages_kept", rb.barrier_stages_kept},
             {"dist_phase_a", rb.dist_phase_a}};

  const DesignResult d = design_constrained(m, m.sigma, X_ini, opt);
  j["design"] = {{"X", pts(d.X)}, {"t", vec(d.t)}, {"feasible", d.feasible}, {"min_cross", d.min_cross},
                 {"min_q", d.min_q}, {"min_mu", d.min_mu}, {"dist_ini", d.dist_ini},
                 {"n_nonconvex_x0", d.n_nonconvex_x0}, {"barrier_stages_kept", d.barrier_stages_kept},
                 {"ch", ch_json(d.ch)}};
  // characterize at the C++ design point, so the Julia characterize can be tested on the SAME X
  j["characterize_at_design_X"] = ch_json(characterize(m, m.sigma, d.X));
  j["characterize_at_X0"] = ch_json(characterize(m, m.sigma, X0));
  j["characterize_at_X_ini"] = ch_json(characterize(m, m.sigma, X_ini));

  std::ofstream f(out);
  f << j.dump(1) << "\n";
  std::cout << "wrote " << out << " theta_max=" << d.ch.theta_max << " rb.iters=" << rb.iterations << "\n";
  return 0;
}
