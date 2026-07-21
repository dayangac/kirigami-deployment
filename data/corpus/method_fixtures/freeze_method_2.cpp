// freeze_method2.cpp -- freezes the reference cases of tests/test_method.cpp (zero_plus /
// convex_embed cases) and tests/test_range_embed.cpp, with the C++ numbers the Julia port
// (Kirigami/test/test_method_2.jl, test_range_embed.jl) must reproduce.
//   clang++ ... freeze_method2.cpp code/build/libkiri_core.a -o freeze_method2
//   ./freeze_method2 <out.json>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>
#include <nlohmann/json.hpp>
#include "helpers.hpp"
#include "core/tutte_auxetic.hpp"
#include "method/convex_embed.hpp"
#include "method/range_embed.hpp"
#include "method/zero_plus.hpp"

using namespace kiri;
using namespace kiri::method;
using json = nlohmann::json;

namespace {

double med_edge(const Mesh& m) {
  std::vector<double> L;
  for (const auto& e : m.edges) L.push_back((m.X[e.key.a] - m.X[e.key.b]).norm());
  if (L.empty()) return 1.0;
  std::nth_element(L.begin(), L.begin() + L.size() / 2, L.end());
  return L[L.size() / 2];
}

json pts(const std::vector<Vec2>& X) {
  json a = json::array();
  for (const auto& p : X) a.push_back({p.x(), p.y()});
  return a;
}
json mat(const Eigen::MatrixXd& M) {
  json a = json::array();
  for (int i = 0; i < M.rows(); ++i) {
    json r = json::array();
    for (int j = 0; j < M.cols(); ++j) r.push_back(M(i, j));
    a.push_back(r);
  }
  return a;
}
json vec(const Eigen::VectorXd& v) {
  json a = json::array();
  for (int i = 0; i < v.size(); ++i) a.push_back(v(i));
  return a;
}

// Deterministic test point in the shape space: t_i = scale * med * sin(i + 1).
Eigen::VectorXd probe_t(int m, double med, double scale) {
  Eigen::VectorXd t(m);
  for (int i = 0; i < m; ++i) t(i) = scale * med * std::sin(i + 1.0);
  return t;
}

json freeze_case(const std::string& name, Mesh mesh, bool checker) {
  std::mt19937 rng(2026);
  mesh.build_topology();
  mesh.sigma = checker ? test::checkerboard_sigma(mesh)
                       : assign_orientation_relaxation(mesh, rng, 8, 500, 180).sigma;
  mesh.build_topology();
  CutStructure c = make_cut(mesh);
  c.mesh = &mesh;
  const HoleSet hs = holes_partition(c);
  const double med = med_edge(mesh);
  const LinearSystem sys = assemble_system(c, hs, mesh.X, BoundaryMode::Fixed);
  const SolveReport sr = solve_system(sys, mesh.X);
  const std::vector<Vec2> X0 = matrix_to_points(sr.X0);
  const bool ini_deployable = hole_residuals(c, mesh.X, hs).deployable(1e-9);
  // make_case's X: X_ini if already deployable, else the projection.
  const std::vector<Vec2> Xc = ini_deployable ? mesh.X : X0;
  const int m = 2 * static_cast<int>(sr.Phi.cols());

  json j;
  j["name"] = name;
  j["provenance"] = "tests/test_method.cpp make_case / tests/test_range_embed.cpp make_rcase: "
                    "std::mt19937(2026), assign_orientation_relaxation(m, rng, 8, 500, 180) "
                    "(or checkerboard_sigma), make_cut, holes_partition, assemble_system(Fixed), "
                    "solve_system; X_case = X_ini if hole_residuals deployable(1e-9) else X0";
  j["mesh"] = json::parse(mesh_to_json_string(mesh));
  j["ini_deployable"] = ini_deployable;
  j["X0"] = pts(X0);
  j["X_case"] = pts(Xc);
  j["Phi"] = mat(sr.Phi);
  j["dim_null"] = sr.dim_null;
  j["med"] = med;
  j["n_split"] = c.n_split();
  j["q_at_case"] = zero_plus_q(c, Xc);
  j["mu_at_case"] = zero_plus_corner_margin(c, Xc);
  j["crosses_at_case"] = corner_crosses(mesh, Xc);
  j["n_incidences"] = corner_incidences(c).size();
  {
    double mq = 0, mmu = 0;
    j["margin_at_X0"] = zero_plus_margin(c, X0, med, &mq, &mmu);
    j["margin_min_q_at_X0"] = mq;
    j["margin_min_mu_at_X0"] = mmu;
    j["modes_at_X0"] = range_embed_modes(c, X0);
  }
  if (m > 0) {
    const Eigen::VectorXd t = probe_t(m, med, 0.02);
    j["probe_t"] = vec(t);
    j["q_form_at_probe"] = vec(zero_plus_eval(zero_plus_form(c, X0, sr.Phi), t));
    {
      ZeroPlusRepairOptions opt;
      opt.lambda_rel = 1e-3;
      opt.w_corner = 1.0;
      opt.w_prox = 0.5;
      Eigen::VectorXd g;
      j["zero_plus_objective_at_probe"] = zero_plus_objective(c, X0, sr.Phi, med, opt, t, &g);
      j["zero_plus_grad_at_probe"] = vec(g);
    }
    for (int variant = 0; variant < 2; ++variant) {
      ConvexEmbedOptions opt;
      opt.delta_rel = 1e-3;
      opt.split_delta_rel = variant ? 1e-3 : 0.0;
      Eigen::VectorXd g;
      const std::string k = variant ? "convex_embed_b" : "convex_embed_a";
      j[k + "_objective_at_probe"] = convex_embed_objective(c, X0, sr.Phi, med, opt, 1.0, t, &g);
      j[k + "_grad_at_probe"] = vec(g);
    }
    {
      RangeEmbedOptions opt;
      opt.delta_rel = 1e-4;
      opt.split_delta_rel = 1e-4;
      opt.kappa = 2e-2;
      const std::vector<int> mode = range_embed_modes(c, shape_point(X0, sr.Phi, t));
      j["range_modes_at_probe"] = mode;
      for (double bw : {0.0, 1e-3}) {
        Eigen::VectorXd g;
        const std::string k = bw > 0 ? "range_embed_objective_bw1e-3" : "range_embed_objective_bw0";
        j[k + "_at_probe"] = range_embed_objective(c, X0, sr.Phi, med, opt, mode, bw, t, &g);
        j[k + "_grad_at_probe"] = vec(g);
      }
    }
    {  // the solve of "convex_embed: a FEASIBLE verdict is always the exact constraint values"
      ConvexEmbedOptions opt;
      opt.delta_rel = 1e-3;
      opt.split_delta_rel = 1e-3;
      opt.n_random = 1;
      opt.max_iter = 150;
      opt.barrier_stages = 2;
      const ConvexEmbedResult r = convex_embed(c, Xc, sr.Phi, Xc, med, opt);
      json s;
      s["feasible"] = r.feasible;
      s["min_cross"] = r.min_cross;
      s["n_bad"] = r.n_bad;
      s["min_q"] = r.min_q;
      s["dist_ini"] = r.dist_ini;
      s["dist_phase_a"] = r.dist_phase_a;
      s["best_start"] = r.best_start;
      s["barrier_stages_kept"] = r.barrier_stages_kept;
      s["n_nonconvex_start"] = r.n_nonconvex_start;
      j["convex_embed_solve"] = s;
    }
    {  // the solve of "range_embed: the reported margin and feasibility are the EXACT values"
      RangeEmbedOptions opt;
      opt.delta_rel = 1e-3;
      opt.split_delta_rel = 1e-3;
      opt.stages = 3;
      opt.iter_per_stage = 60;
      const RangeEmbedResult r = range_embed(c, X0, sr.Phi, med, opt);
      json s;
      s["feasible"] = r.feasible;
      s["margin"] = r.margin;
      s["margin_start"] = r.margin_start;
      s["min_cross"] = r.min_cross;
      s["min_q"] = r.min_q;
      s["min_mu"] = r.min_mu;
      s["n_bad"] = r.n_bad;
      s["n_bad_q"] = r.n_bad_q;
      s["stages_kept"] = r.stages_kept;
      s["best_start"] = r.best_start;
      s["n_entries"] = r.n_entries;
      j["range_embed_solve"] = s;
    }
    {  // the solve of "0+ repair with the corner term"
      ZeroPlusRepairOptions opt;
      opt.n_random = 2;
      opt.max_iter = 200;
      opt.w_corner = 1.0;
      opt.w_prox = 1e-2;
      const ZeroPlusRepairResult r = zero_plus_repair(c, Xc, sr.Phi, med, opt);
      json s;
      s["feasible"] = r.feasible;
      s["min_q"] = r.min_q;
      s["min_area"] = r.min_area;
      s["min_margin"] = r.min_margin;
      s["n_bad_margin"] = r.n_bad_margin;
      s["min_margin_start"] = r.min_margin_start;
      s["n_bad_margin_start"] = r.n_bad_margin_start;
      s["best_start"] = r.best_start;
      s["f_start"] = r.f_start;
      s["f_end"] = r.f_end;
      j["zero_plus_repair_solve"] = s;
    }
  }
  return j;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) return 1;
  json out;
  out["provenance"] = "freeze_method2.cpp (scratch freezer, kirigami-julia port) linked against "
                      "code/build/libkiri_core.a, arm64 clang -O2 -ffp-contract=on";
  json cases = json::array();
  const Vec2 ctr(0.13, 0.07);
  cases.push_back(freeze_case("hexagons_3.0", tiling_hexagons(disk(ctr, 3.0)), false));
  cases.push_back(freeze_case("snub_square_3.2", tiling_snub_square(disk(ctr, 3.2)), false));
  cases.push_back(freeze_case("truncated_square_4.0", tiling_truncated_square(disk(ctr, 4.0)), false));
  cases.push_back(freeze_case("truncated_square_3.0", tiling_truncated_square(disk(ctr, 3.0)), false));
  cases.push_back(freeze_case("3_4_3_12_4.2", tiling_3_4_3_12(disk(ctr, 4.2)), false));
  cases.push_back(freeze_case("kagome_2.5_checker", tiling_kagome(disk(ctr, 2.5)), true));
  out["cases"] = cases;
  std::ofstream f(argv[1]);
  f << out.dump(1) << "\n";
  return 0;
}
