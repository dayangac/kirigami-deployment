// freeze_k9c.cpp -- every optimiser call of design_range_max for one K9c row, with its
// inputs (X0, Phi, t_init) and outputs, so the Julia port can be checked call by call on
// IDENTICAL inputs. Built against code/build/libkiri_core.a (Apple clang -O2 arm64).
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include "kill_common.hpp"
#include "method/design.hpp"
#include "method/convex_embed.hpp"
#include "method/range_embed.hpp"

using namespace kiri;
using namespace kiri::kill;
using namespace kiri::method;
using json = nlohmann::json;

static json pts(const std::vector<Vec2>& X) { json a = json::array(); for (const auto& p : X) a.push_back({p.x(), p.y()}); return a; }
static json vec(const Eigen::VectorXd& v) { json a = json::array(); for (int i = 0; i < v.size(); ++i) a.push_back(v[i]); return a; }
static json mat(const Eigen::MatrixXd& M) { json a = json::array(); for (int i = 0; i < M.rows(); ++i) { json r = json::array(); for (int j = 0; j < M.cols(); ++j) r.push_back(M(i, j)); a.push_back(r);} return a; }
static json ce(const ConvexEmbedResult& r) {
  return {{"t", vec(r.t)}, {"X", pts(r.X)}, {"feasible", r.feasible}, {"n_bad", r.n_bad}, {"n_bad_q", r.n_bad_q},
          {"iterations", r.iterations}, {"best_start", r.best_start}, {"min_cross", r.min_cross}, {"min_q", r.min_q},
          {"barrier_stages_kept", r.barrier_stages_kept}, {"dist_phase_a", r.dist_phase_a}};
}
static json re(const RangeEmbedResult& r) {
  return {{"t", vec(r.t)}, {"X", pts(r.X)}, {"feasible", r.feasible}, {"margin", r.margin}, {"margin_start", r.margin_start},
          {"min_cross", r.min_cross}, {"min_q", r.min_q}, {"min_mu", r.min_mu}, {"stages_kept", r.stages_kept},
          {"iterations", r.iterations}, {"best_start", r.best_start}, {"n_entries", r.n_entries}};
}

int main(int argc, char** argv) {
  const int id = argc > 1 ? std::atoi(argv[1]) : 130;
  const std::string out = argc > 2 ? argv[2] : "k9c_calls_130.json";
  Graph g = make_graph(id, 100, 800, 1400);
  Mesh m = g.mesh; m.build_topology();
  const std::vector<Vec2> X_ini = m.X;
  RangeMaxOptions opt; opt.seed = 9300u + 7u * static_cast<unsigned>(id);
  const CutStructure c = make_cut(m);
  const double med = median_edge_length(m);
  const HoleSet hs = holes_partition(c);
  const SolveReport sr = solve_system(assemble_system(c, hs, m.X, BoundaryMode::Fixed), m.X);
  const std::vector<Vec2> X0 = matrix_to_points(sr.X0);
  json j;
  j["provenance"] = "freeze_k9c.cpp (scratchpad): make_graph(id,100,800,1400), RangeMaxOptions defaults, seed 9300+7*id; "
                    "every optimiser call of design_range_max with its inputs; clang -O2 arm64, code/build/libkiri_core.a";
  j["id"] = id; j["seed"] = opt.seed; j["med_edge"] = med; j["dim_null"] = sr.dim_null;
  j["X0"] = pts(X0); j["Phi"] = mat(sr.Phi);

  ConvexEmbedOptions o9; o9.delta_rel = opt.delta_convex; o9.split_delta_rel = opt.delta_split; o9.n_random = 3;
  o9.max_iter = opt.max_iter; o9.barrier_stages = 6; o9.seed = opt.seed;
  const ConvexEmbedResult r9 = convex_embed(c, X0, sr.Phi, X_ini, med, o9);
  j["r9"] = ce(r9);
  ConvexEmbedOptions o9b; o9b.delta_rel = opt.delta_wide; o9b.split_delta_rel = opt.delta_wide; o9b.n_random = 8;
  o9b.max_iter = opt.max_iter; o9b.barrier_stages = 10; o9b.seed = opt.seed + 101u;
  if (r9.feasible) o9b.t_init = r9.t;
  const ConvexEmbedResult r9b = convex_embed(c, X0, sr.Phi, X_ini, med, o9b);
  j["r9b"] = ce(r9b);
  j["stage_a"] = json::array();
  struct Start { std::string tag; Eigen::VectorXd t; double drel; };
  std::vector<Start> starts;
  if (r9.feasible) starts.push_back({"k9", r9.t, opt.delta_convex});
  if (r9b.feasible) starts.push_back({"k9b", r9b.t, opt.delta_wide});
  starts.push_back({"x0", Eigen::VectorXd::Zero(2 * sr.dim_null), opt.delta_convex});
  for (const Start& st : starts) {
    RangeEmbedOptions ro; ro.delta_rel = st.drel; ro.split_delta_rel = st.drel; ro.kappa = opt.kappa;
    ro.stages = opt.stages; ro.iter_per_stage = opt.iter_per_stage; ro.seed = opt.seed + 211u; ro.t_init = st.t;
    const RangeEmbedResult rr = range_embed(c, X0, sr.Phi, med, ro);
    json e = re(rr); e["tag"] = st.tag; e["drel"] = st.drel; e["t_init"] = vec(st.t);
    e["theta_max"] = characterize(m, m.sigma, rr.X).theta_max;
    e["margin_exact"] = zero_plus_margin(c, rr.X, med);
    j["stage_a"].push_back(e);
  }
  const RangeMaxResult rm = design_range_max(m, m.sigma, X_ini, opt);
  j["result"] = {{"provenance", rm.provenance}, {"stage_b_used", rm.stage_b_used}, {"margin", rm.margin},
                 {"theta_max", rm.design.ch.theta_max}, {"eps_max", rm.design.ch.eps_max}, {"X", pts(rm.design.X)}};
  j["arms"] = json::array();
  for (const auto& a : rm.arms) j["arms"].push_back({{"tag", a.tag}, {"feasible", a.feasible}, {"margin", a.margin}, {"theta_max", a.theta_max}, {"eps_max", a.eps_max}, {"stage_a_winner", a.stage_a_winner}});
  // stage B from the stage-A winner, replayed on its own
  {
    const RangeMaxArm* win = nullptr;
    for (const auto& a : rm.arms) if (a.stage_a_winner) win = &a;
    if (win) {
      for (const auto& e : j["stage_a"]) if ("k9c/" + e["tag"].get<std::string>() == win->tag) {
        std::vector<Vec2> Xw; for (const auto& p : e["X"]) Xw.emplace_back(p[0].get<double>(), p[1].get<double>());
        MarginRangeOptions mo; mo.delta_rel = (win->tag == "k9c/k9b") ? opt.delta_wide : opt.delta_convex; mo.split_delta_rel = mo.delta_rel;
        const MarginRangeResult mr = maximize_margin_range(c, Xw, sr.Phi, med, mo);
        j["stage_b"] = {{"from", win->tag}, {"X", pts(mr.X)}, {"improved", mr.improved}, {"theta_before", mr.theta_before},
                        {"theta_after", mr.theta_after}, {"margin_after", mr.margin_after}, {"caps_tried", mr.caps_tried}, {"caps_accepted", mr.caps_accepted}};
      }
    }
  }
  std::ofstream f(out); f << j.dump(1) << "\n";
  std::cout << "wrote " << out << " best=" << rm.provenance << " theta=" << rm.design.ch.theta_max << "\n";
  return 0;
}
