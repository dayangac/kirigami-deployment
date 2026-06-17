// freeze_fixtures.cpp -- replays the mesh/sigma (and deployed-Y) sequences of
// tests/test_mesh_cut.cpp and tests/test_holes.cpp and dumps them as JSON so the
// Julia port can run the identical checks without the (unported) generators.
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include "helpers.hpp"

using namespace kiri;
using namespace kiri::test;
using nlohmann::json;

static json mesh_json(const Mesh& m) { return json::parse(mesh_to_json_string(m)); }
static json pts_json(const std::vector<Vec2>& P) {
  json a = json::array();
  for (const auto& p : P) a.push_back({p.x(), p.y()});
  return a;
}
static json edges_json(const std::vector<std::vector<int>>& v) { return json(v); }

// ------------------------------------------------------------ test_mesh_cut.cpp
static json freeze_mesh_cut() {
  json out;
  {
    std::mt19937 rng(12345);
    const std::vector<std::string> kinds = {"squares", "triangles", "hexagons", "kagome",
                                            "snub_square", "truncated_square", "t3_4_3_12"};
    json fam = json::array();
    for (const auto& kind : kinds) {
      Mesh base = generate(kind, {4.0}, rng);
      json sigmas = json::array();
      for (int trial = 0; trial < 30; ++trial) sigmas.push_back(random_sigma(base, rng));
      fam.push_back({{"kind", kind}, {"mesh", mesh_json(base)}, {"sigmas", sigmas}});
    }
    json dl = json::array();
    for (int trial = 0; trial < 40; ++trial) {
      Mesh g = delaunay_of_random_points(30 + trial, 10.0, rng);
      g.sigma = random_sigma(g, rng);
      dl.push_back(mesh_json(g));
    }
    out["remark_A1"] = {{"families", fam}, {"delaunay", dl},
                        {"provenance", {{"test_case", "Remark A.1: in-degree == out-degree at every interior vertex"},
                                        {"source", "tests/test_mesh_cut.cpp"},
                                        {"rng", "std::mt19937 rng(12345), shared across the whole case in this order"},
                                        {"families", "for kind in kinds: base = generate(kind, {4.0}, rng); 30 x random_sigma(base, rng)"},
                                        {"delaunay", "for trial in 0..39: g = delaunay_of_random_points(30 + trial, 10.0, rng); g.sigma = random_sigma(g, rng) (stored as orientation)"},
                                        {"sigma_method", "test::random_sigma = std::bernoulli_distribution(0.5) per face, true -> +1"}}}};
  }
  {
    std::mt19937 rng(7);
    Mesh g = generate("kagome", {4.0}, rng);
    json sigmas = json::array();
    for (int trial = 0; trial < 20; ++trial) sigmas.push_back(random_sigma(g, rng));
    out["A2A3_kagome"] = {{"mesh", mesh_json(g)}, {"sigmas", sigmas},
                          {"provenance", {{"test_case", "M' vertex count matches the Remark A.2 / A.3 bookkeeping"},
                                          {"source", "tests/test_mesh_cut.cpp"},
                                          {"rng", "std::mt19937 rng(7)"},
                                          {"mesh", "generate(\"kagome\", {4.0}, rng)"},
                                          {"sigmas", "20 x random_sigma(g, rng)"},
                                          {"sigma_method", "test::random_sigma = std::bernoulli_distribution(0.5) per face, true -> +1"}}}};
  }
  return out;
}

// ---------------------------------------------------------------- test_holes.cpp
static json freeze_holes() {
  json out;
  // Case 1: seed-growing == partition == geometry
  {
    std::mt19937 rng(2024);
    const std::vector<std::string> kinds = {"squares", "triangles", "hexagons", "kagome",
                                            "snub_square", "truncated_square", "t3_4_3_12"};
    json graphs = json::array();
    int ngraphs = 0, geometric_checks = 0, split_cycles = 0, disconnected = 0, cycle_mismatch = 0;
    auto check_one = [&](Mesh& g, const std::string& tag) {
      json rec;
      rec["tag"] = tag;
      rec["mesh"] = mesh_json(g);
      const CutStructure c = make_cut(g);
      const HoleSet b = holes_partition(c);
      ++ngraphs;
      const bool forest = split_subgraph_is_forest(c);
      const bool connected = count_components(c) == 1;
      if (!forest) ++split_cycles;
      if (!connected) ++disconnected;
      rec["forest"] = forest;
      rec["connected"] = connected;
      rec["Yd"] = nullptr;
      const LinearSystem sys = assemble_system(c, b, g.X, BoundaryMode::Fixed);
      if (sys.A.rows() == 0 || g.n_vertices() > 400) { graphs.push_back(rec); return; }
      const SolveReport rep = solve_system(sys, g.X);
      const auto X0 = matrix_to_points(rep.X0);
      if (hole_residuals(c, X0, b).max_norm > 1e-9) { graphs.push_back(rec); return; }
      double th = 0.25;
      std::vector<Vec2> Yd;
      bool clean = false;
      for (int k = 0; k < 6 && !clean; ++k, th *= 0.35) {
        Yd = deploy(c, X0, th).Y;
        clean = !has_collision(c, Yd);
      }
      if (!clean) { graphs.push_back(rec); return; }
      const auto geo = holes_geometric(c, Yd);
      std::vector<std::vector<int>> comb;
      for (int i : b.interior_indices) comb.push_back(b.all[i].edges);
      std::sort(comb.begin(), comb.end());
      rec["Yd"] = pts_json(Yd);
      rec["geo"] = edges_json(geo);
      if (forest && connected) {
        if (geo != comb) std::cerr << "MISMATCH in forest&&connected case " << tag << "\n";
        ++geometric_checks;
      } else if (geo != comb) {
        ++cycle_mismatch;
      }
      graphs.push_back(rec);
    };
    for (const auto& kind : kinds) {
      Mesh base = generate(kind, {3.2}, rng);
      for (int t = 0; t < 5; ++t) {
        base.sigma = (t == 0) ? checkerboard_sigma(base) : random_sigma(base, rng);
        check_one(base, kind + "/" + std::to_string(t));
      }
    }
    for (int t = 0; t < 30; ++t) {
      Mesh g = (t % 2) ? voronoi_of_random_points(24 + t, 8.0, rng)
                       : delaunay_of_random_points(20 + t, 8.0, rng);
      if (g.n_faces() < 6) continue;
      g.sigma = random_sigma(g, rng);
      check_one(g, std::string(t % 2 ? "voronoi/" : "delaunay/") + std::to_string(t));
    }
    out["case1"] = {{"graphs", graphs},
                    {"provenance", {{"test_case", "hole preimages: seed-growing == partition formulation == geometry"},
                                    {"source", "tests/test_holes.cpp"},
                                    {"rng", "std::mt19937 rng(2024), shared across the whole case in this order"},
                                    {"families", "for kind in kinds: base = generate(kind, {3.2}, rng); t=0: checkerboard_sigma(base), t=1..4: random_sigma(base, rng); tag = kind/t"},
                                    {"random", "for t in 0..29: odd t -> voronoi_of_random_points(24 + t, 8.0, rng), even t -> delaunay_of_random_points(20 + t, 8.0, rng); skip if n_faces < 6; sigma = random_sigma(g, rng); tag = voronoi|delaunay/t"},
                                    {"Yd", "null if the C++ skipped the geometric stage; else deploy(c, X0, th).Y with X0 = solve_system(assemble_system(c, holes_partition(c), g.X, Fixed), g.X).X0 and th the first of 0.25*0.35^k (k<6) with !has_collision"},
                                    {"geo", "C++ holes_geometric(c, Yd), 0-based edge ids"}}},
                    {"counts", {{"graphs", ngraphs}, {"geometric_checks", geometric_checks},
                                {"split_cycles", split_cycles}, {"disconnected", disconnected},
                                {"cycle_mismatch", cycle_mismatch}}}};
  }
  // Case 2: rotating squares
  {
    Mesh m = tiling_squares(rect(Vec2(2, 2), 2.01, 2.01));
    m.sigma = checkerboard_sigma(m);
    out["rotating_squares"] = mesh_json(m);
    out["rotating_squares"]["provenance"] = {{"test_case", "hole preimages of the rotating-squares pattern are the interior vertices"},
                                             {"source", "tests/test_holes.cpp"},
                                             {"mesh", "tiling_squares(rect(Vec2(2, 2), 2.01, 2.01))"},
                                             {"sigma_method", "test::checkerboard_sigma (stored as orientation)"}};
  }
  // Case 3: split cuts merge holes
  {
    std::mt19937 rng(99);
    json fam = json::array();
    for (const std::string& kind : {"hexagons", "truncated_square", "snub_square", "squares"}) {
      Mesh g = generate(kind, {3.2}, rng);
      json sigmas = json::array();
      for (int t = 0; t < 8; ++t) sigmas.push_back(random_sigma(g, rng));
      fam.push_back({{"kind", kind}, {"mesh", mesh_json(g)}, {"sigmas", sigmas}});
    }
    out["case3"] = fam;
    out["case3_provenance"] = {{"test_case", "split cuts merge holes: H = #interior - #interior split edges (forest case)"},
                               {"source", "tests/test_holes.cpp"},
                               {"rng", "std::mt19937 rng(99), shared across the whole case in this order"},
                               {"families", "for kind in {hexagons, truncated_square, snub_square, squares}: g = generate(kind, {3.2}, rng); 8 x random_sigma(g, rng)"},
                               {"sigma_method", "test::random_sigma = std::bernoulli_distribution(0.5) per face, true -> +1"}};
  }
  // Case 4: geometric agreement on >= 50 random instances
  {
    std::mt19937 rng(31337);
    int checked = 0, skipped_overlap = 0, random_topology_checked = 0, random_topology_total = 0;
    json instances = json::array();
    auto try_check = [&](const Mesh& g, const std::vector<Vec2>& X, const std::string& tag) {
      json rec;
      rec["tag"] = tag;
      rec["mesh"] = mesh_json(g);
      rec["X"] = pts_json(X);
      rec["Yd"] = nullptr;
      const CutStructure c = make_cut(g);
      if (!check_remark_A1(c)) std::cerr << "A1 failure " << tag << "\n";
      const HoleSet b = holes_partition(c);
      auto finish = [&](bool r) { instances.push_back(rec); return r; };
      if (!split_subgraph_is_forest(c) || count_components(c) != 1) return finish(false);
      if (hole_residuals(c, X, b).max_norm > 1e-9) return finish(false);
      double th = 0.2;
      std::vector<Vec2> Yd;
      bool clean = false;
      for (int k = 0; k < 8 && !clean; ++k, th *= 0.4) {
        Yd = deploy(c, X, th).Y;
        clean = !has_collision(c, Yd);
      }
      if (!clean) { ++skipped_overlap; return finish(false); }
      const auto geo = holes_geometric(c, Yd);
      std::vector<std::vector<int>> comb;
      for (int i : b.interior_indices) comb.push_back(b.all[i].edges);
      std::sort(comb.begin(), comb.end());
      if (geo != comb) std::cerr << "MISMATCH case4 " << tag << "\n";
      rec["Yd"] = pts_json(Yd);
      rec["geo"] = edges_json(geo);
      ++checked;
      return finish(true);
    };
    for (int t = 0; t < 30; ++t) {
      Mesh g = (t % 3 == 0)   ? delaunay_of_random_points(18 + (t % 11), 8.0, rng)
               : (t % 3 == 1) ? voronoi_of_random_points(22 + (t % 13), 8.0, rng)
                              : quad_dominant_random(20 + (t % 9), 8.0, rng);
      if (g.n_faces() < 6) continue;
      g.sigma = assign_orientation_relaxation(g, rng, 4, 250, 90).sigma;
      const CutStructure c = make_cut(g);
      const HoleSet b = holes_partition(c);
      const LinearSystem sys = assemble_system(c, b, g.X, BoundaryMode::Fixed);
      if (sys.A.rows() == 0) continue;
      const SolveReport rep = solve_system(sys, g.X);
      ++random_topology_total;
      if (try_check(g, matrix_to_points(rep.X0), "A/" + std::to_string(t))) ++random_topology_checked;
    }
    std::uniform_real_distribution<double> U(-0.35, 0.35);
    const std::vector<std::string> kinds = {"squares", "triangles", "kagome",  "hexagons",
                                            "snub_square", "truncated_square", "t3_4_3_12"};
    for (int t = 0; t < 90 && checked < 60; ++t) {
      Mesh g = generate(kinds[t % kinds.size()], {2.6 + 0.2 * (t % 4)}, rng);
      if (g.n_faces() < 6) continue;
      g.sigma = (t % 2) ? checkerboard_sigma(g) : assign_orientation_relaxation(g, rng, 3, 200, 60).sigma;
      const CutStructure c = make_cut(g);
      const HoleSet b = holes_partition(c);
      const LinearSystem sys = assemble_system(c, b, g.X, BoundaryMode::Fixed);
      if (sys.A.rows() == 0) continue;
      const SolveReport rep = solve_system(sys, g.X);
      Eigen::MatrixXd X = rep.X0;
      if (rep.dim_null > 0) {
        Eigen::MatrixXd T(rep.dim_null, 2);
        for (int i = 0; i < rep.dim_null; ++i) {
          T(i, 0) = U(rng);
          T(i, 1) = U(rng);
        }
        X += rep.Phi * T;
      }
      Eigen::Matrix2d A;
      A << 1 + 0.3 * U(rng), 0.4 * U(rng), 0.4 * U(rng), 1 + 0.3 * U(rng);
      if (std::abs(A.determinant()) < 0.2) A.setIdentity();
      X = X * A.transpose();
      try_check(g, matrix_to_points(X), "B/" + std::to_string(t));
    }
    out["case4"] = {{"instances", instances},
                    {"provenance", {{"test_case", "hole preimages: geometric agreement on >= 50 random instances"},
                                    {"source", "tests/test_holes.cpp"},
                                    {"rng", "std::mt19937 rng(31337), shared across the whole case in this order; U = uniform_real_distribution(-0.35, 0.35)"},
                                    {"A", "for t in 0..29: t%3==0 delaunay_of_random_points(18 + t%11, 8.0, rng); t%3==1 voronoi_of_random_points(22 + t%13, 8.0, rng); else quad_dominant_random(20 + t%9, 8.0, rng); skip n_faces<6; sigma = assign_orientation_relaxation(g, rng, 4, 250, 90).sigma; skip if assemble_system rows == 0; X = solve_system(...).X0; tag = A/t"},
                                    {"B", "for t in 0..89 while checked < 60: g = generate(kinds[t%7], {2.6 + 0.2*(t%4)}, rng), kinds = squares,triangles,kagome,hexagons,snub_square,truncated_square,t3_4_3_12; skip n_faces<6; odd t: checkerboard_sigma, even t: assign_orientation_relaxation(g, rng, 3, 200, 60).sigma; X = X0 + Phi*T (T dim_null x 2 of U draws, row-major) then X * A^T with A = [1+0.3U, 0.4U; 0.4U, 1+0.3U] (identity if |det| < 0.2); tag = B/t"},
                                    {"X", "the embedding passed to try_check (M vertex positions)"},
                                    {"Yd", "null if try_check returned before deploying (not forest, disconnected, residual > 1e-9, or all 8 deployments th = 0.2*0.4^k collided); else deploy(c, X, th).Y for the first collision-free th"},
                                    {"geo", "C++ holes_geometric(c, Yd), 0-based edge ids"}}},
                    {"counts", {{"checked", checked}, {"skipped_overlap", skipped_overlap},
                                {"random_topology_checked", random_topology_checked},
                                {"random_topology_total", random_topology_total}}}};
  }
  return out;
}

int main(int argc, char** argv) {
  const std::string dir = argc > 1 ? argv[1] : ".";
  {
    std::ofstream o(dir + "/test_fixtures_mesh_cut.json");
    o << freeze_mesh_cut().dump() << "\n";
  }
  {
    std::ofstream o(dir + "/test_fixtures_holes.json");
    o << freeze_holes().dump() << "\n";
  }
  std::cout << "ok\n";
}
