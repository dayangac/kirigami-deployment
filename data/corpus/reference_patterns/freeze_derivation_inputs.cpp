// freeze_derivation_inputs.cpp -- dumps every INPUT population of
// tests/derivation_tests.cpp (corpus(), the H-LOC patches, the T7(iv) tori, l1_corpus(),
// the L2 K7 population + 110 fresh Voronoi tori + the 2000-draw search) as JSON, so the
// Julia port can check its own generators/relaxation against the C++ and fall back to the
// frozen inputs where the relaxation is not reproduced bit-exactly.
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <nlohmann/json.hpp>
#include "core/collision.hpp"
#include "core/cut.hpp"
#include "core/generators.hpp"
#include "core/holes.hpp"
#include "core/kinematics.hpp"
#include "core/mesh.hpp"
#include "core/orientation.hpp"
#include "core/tutte_auxetic.hpp"
#include "method/periodic_jacobian.hpp"
#include "kill_common.hpp"

using namespace kiri;
using nlohmann::json;
using kiri::method::PeriodicPattern;
using kiri::method::Quotient;

static json mj_(const Mesh& m) { return json::parse(mesh_to_json_string(m)); }

// ------------------------------------------------------------ corpus() (verbatim)
static json mat_json(const Eigen::MatrixXd& M) {
  json a = json::array();
  for (int i = 0; i < M.rows(); ++i) {
    json r = json::array();
    for (int j = 0; j < M.cols(); ++j) r.push_back(M(i, j));
    a.push_back(r);
  }
  return a;
}
struct CaseLite { std::string name; Mesh m; bool ok; int k; int H; int nh; int ns; Eigen::MatrixXd X0, Phi; };
static CaseLite build_lite(Mesh m, const std::string& name) {
  CaseLite c{name, std::move(m), false, 0, 0, 0, 0, {}, {}};
  const CutStructure cut = make_cut(c.m);
  const HoleSet hs = holes_partition(cut);
  const LinearSystem sys = assemble_system(cut, hs, c.m.X, BoundaryMode::Fixed);
  const SolveReport rep = solve_system(sys, c.m.X);
  const auto X0 = matrix_to_points(rep.X0);
  const Residuals r = hole_residuals(cut, X0, hs);
  c.ok = rep.projection_ok && r.max_norm < 1e-7 && cut.n_hinge() > 0;
  c.k = rep.dim_null; c.H = hs.n_interior_holes(); c.nh = cut.n_hinge(); c.ns = cut.n_split();
  c.X0 = rep.X0; c.Phi = rep.Phi;
  return c;
}
static json freeze_corpus() {
  json out = json::array();
  std::mt19937 rng(20260903u);
  auto add = [&](Mesh m, const std::string& nm) {
    if (m.n_faces() < 3) return;
    m.sigma = assign_orientation_relaxation(m, rng, 6, 400, 90).sigma;
    CaseLite c = build_lite(std::move(m), nm);
    out.push_back({{"name", nm}, {"mesh", mj_(c.m)}, {"ok", c.ok}, {"k", c.k}, {"H", c.H},
                   {"n_hinge", c.nh}, {"n_split", c.ns}, {"X0", mat_json(c.X0)}, {"Phi", mat_json(c.Phi)}});
  };
  add(tiling_squares(disk({0, 0}, 2.6)), "squares");
  add(tiling_triangles(disk({0, 0}, 2.2)), "triangles");
  add(tiling_hexagons(disk({0, 0}, 3.0)), "hexagons");
  add(tiling_kagome(disk({0, 0}, 2.4)), "kagome");
  add(tiling_snub_square(disk({0, 0}, 2.2)), "snub_square");
  add(tiling_truncated_square(disk({0, 0}, 2.8)), "truncated_square");
  add(tiling_3_4_3_12(disk({0, 0}, 3.2)), "t3_4_3_12");
  for (int s = 1; s <= 3; ++s) {
    std::mt19937 r2(1000u + s);
    add(largest_component(delaunay_of_random_points(45, 10.0, r2)), "delaunay" + std::to_string(s));
    std::mt19937 r3(2000u + s);
    add(largest_component(voronoi_of_random_points(40, 10.0, r3)), "voronoi" + std::to_string(s));
    std::mt19937 r4(3000u + s);
    add(largest_component(quad_dominant_random(40, 10.0, r4)), "quad" + std::to_string(s));
  }
  return out;
}

static json freeze_hloc() {
  json out = json::array();
  for (double R : {1.6, 2.6, 3.6, 4.6, 5.6, 7.0}) {
    Mesh m = tiling_squares(disk({0, 0}, R));
    if (m.n_faces() < 4) continue;
    std::mt19937 rng(5u);
    m.sigma = assign_orientation_relaxation(m, rng, 4, 300, 60).sigma;
    out.push_back({{"R", R}, {"mesh", mj_(m)}});
  }
  return out;
}

static json freeze_tori() {
  json out = json::array();
  for (int mk : {0, 1, 2}) {
    Mesh m = (mk == 0) ? torus_squares(4, 4) : (mk == 1) ? torus_squares(6, 4) : torus_triangles(4, 4);
    std::mt19937 rng(9u);
    m.sigma = assign_orientation_relaxation(m, rng, 4, 300, 60).sigma;
    out.push_back({{"mk", mk}, {"mesh", mj_(m)}});
  }
  return out;
}

static json freeze_l1() {
  json out = json::array();
  int n = 0;
  for (auto& rc : kill::reference_cases()) {
    if (rc.periodic) continue;
    CaseLite c = build_lite(std::move(rc.mesh), rc.name);
    out.push_back({{"name", rc.name}, {"mesh", mj_(c.m)}, {"ok", c.ok}, {"k", c.k}, {"X0", mat_json(c.X0)}, {"Phi", mat_json(c.Phi)}});
    if (c.ok) ++n;
  }
  for (int id = 0; id < 80 && n < 60; ++id) {
    kill::Graph g = kill::make_graph(id, 18, 46, 220);
    if (!g.ok) continue;
    CaseLite c = build_lite(std::move(g.mesh), g.kind + std::to_string(id));
    out.push_back({{"name", c.name}, {"id", id}, {"mesh", mj_(c.m)}, {"ok", c.ok}, {"k", c.k}, {"X0", mat_json(c.X0)}, {"Phi", mat_json(c.Phi)}});
    if (c.ok) ++n;
  }
  return out;
}

// ---- copied from derivation_tests.cpp / kill_k7.cpp (population) ---------------
static PeriodicPattern l2_voronoi_pattern(int inst, int nsites, double L, std::mt19937& rng) {
  PeriodicPattern P;
  P.family = "voronoi_torus";
  P.name = "voronoi_torus_" + std::to_string(inst) + "_n" + std::to_string(nsites);
  std::uniform_real_distribution<double> U(0, L);
  std::vector<Vec2> site;
  int guard = 0;
  while (static_cast<int>(site.size()) < nsites && guard++ < 100000) {
    const Vec2 p(U(rng), U(rng));
    bool ok = true;
    for (const auto& s : site) {
      Vec2 d = p - s;
      d.x() -= L * std::round(d.x() / L);
      d.y() -= L * std::round(d.y() / L);
      if (d.norm() < 0.15 * L / std::sqrt(static_cast<double>(nsites))) ok = false;
    }
    if (ok) site.push_back(p);
  }
  if (static_cast<int>(site.size()) < nsites) { P.err = "site rejection failed"; return P; }
  std::vector<Vec2> rep;
  for (int i = -1; i <= 1; ++i)
    for (int j = -1; j <= 1; ++j)
      for (const auto& s : site) rep.push_back(s + Vec2(i * L, j * L));
  std::vector<std::vector<Vec2>> polys;
  for (int i = 0; i < nsites; ++i) {
    const Vec2 pi = site[i];
    std::vector<Vec2> cellp{pi + Vec2(-L, -L), pi + Vec2(L, -L), pi + Vec2(L, L), pi + Vec2(-L, L)};
    for (const auto& pj : rep) {
      if ((pj - pi).norm() < 1e-12) continue;
      const Vec2 nvec = pj - pi;
      const double off = nvec.dot(0.5 * (pi + pj));
      std::vector<Vec2> out;
      for (size_t k = 0; k < cellp.size(); ++k) {
        const Vec2& A = cellp[k];
        const Vec2& Bv = cellp[(k + 1) % cellp.size()];
        const double da = nvec.dot(A) - off, db = nvec.dot(Bv) - off;
        if (da <= 0) out.push_back(A);
        if ((da < 0 && db > 0) || (da > 0 && db < 0)) out.push_back(A + (Bv - A) * (da / (da - db)));
      }
      cellp = out;
      if (cellp.size() < 3) break;
    }
    if (cellp.size() < 3) { P.err = "empty voronoi cell"; return P; }
    polys.push_back(cellp);
  }
  Mesh cell;
  try { cell = mesh_from_polygons(polys, 1e-7); }
  catch (const std::exception& e) { P.err = std::string("weld: ") + e.what(); return P; }
  Eigen::Matrix2d T;
  T << L, 0, 0, L;
  cell.sigma.assign(cell.n_faces(), -1);
  Quotient q0 = method::build_quotient(cell, T);
  if (!q0.ok) { P.err = "quotient(probe): " + q0.err; return P; }
  cell.sigma = method::quotient_sigma(cell.n_faces(), method::quotient_dual(q0), rng);
  P.cell = cell;
  P.T = T;
  P.ok = true;
  return P;
}
static json pattern_json(const PeriodicPattern& P) {
  json j = {{"name", P.name}, {"family", P.family}, {"ok", P.ok}, {"err", P.err},
            {"T", {P.T(0, 0), P.T(0, 1), P.T(1, 0), P.T(1, 1)}}};
  if (P.ok) j["cell"] = mj_(P.cell);
  return j;
}
static json freeze_l2() {
  json out;
  json k7 = json::array();
  const char* fams[] = {"squares", "triangles", "hexagons", "kagome", "snub_square", "trunc_square_488", "t3_4_3_12"};
  const int sizes[][2] = {{2, 2}, {3, 2}, {3, 3}};
  int idx = 0;
  for (const char* f : fams)
    for (const auto& s : sizes) {
      std::mt19937 rng(20260904u + 7919u * static_cast<unsigned>(idx++));
      k7.push_back(pattern_json(method::make_tiling_pattern(f, s[0], s[1], rng)));
    }
  const int ns[] = {20, 28, 36, 45, 55, 70, 85, 100, 120, 140, 170, 200};
  for (int i = 0; i < 12; ++i) {
    std::mt19937 rng(9100001u + 104729u * static_cast<unsigned>(i));
    k7.push_back(pattern_json(l2_voronoi_pattern(i, ns[i], 10.0, rng)));
  }
  out["k7"] = k7;
  json fresh = json::array();
  for (int i = 0; i < 110; ++i) {
    std::mt19937 rng(4400011u + 7717u * static_cast<unsigned>(i));
    std::uniform_int_distribution<int> nsd(12, 45);
    fresh.push_back(pattern_json(l2_voronoi_pattern(1000 + i, nsd(rng), 10.0, rng)));
  }
  out["fresh110"] = fresh;
  json ce = json::array();
  for (int i = 0; i < 2000; ++i) {
    std::mt19937 rng(880011u + 65537u * static_cast<unsigned>(i));
    std::uniform_int_distribution<int> nsd(6, 16);
    ce.push_back(pattern_json(l2_voronoi_pattern(90000 + i, nsd(rng), 10.0, rng)));
  }
  out["ce2000"] = ce;
  return out;
}

int main(int argc, char** argv) {
  const std::string dir = argc > 1 ? argv[1] : ".";
  json out;
  out["corpus"] = freeze_corpus();
  out["hloc"] = freeze_hloc();
  out["tori"] = freeze_tori();
  out["l1_corpus"] = freeze_l1();
  out["provenance"] = "freeze_derivation_inputs.cpp: verbatim replay of the input populations of "
                      "tests/derivation_tests.cpp (corpus(): mt19937(20260903) shared over the 16 "
                      "assign_orientation_relaxation(m,rng,6,400,90) calls in file order; hloc: "
                      "tiling_squares(disk(0,R)) + relaxation(mt19937(5),4,300,60); tori: torus_* + "
                      "relaxation(mt19937(9),4,300,60); l1_corpus: kill::reference_cases() non-periodic + "
                      "kill::make_graph(id,18,46,220) until 60 ok). Meshes are 0-based JSON with sigma as "
                      "'orientation'; torus faces must be loaded AS STORED. corpus/l1_corpus rows also carry the C++ "
                      "solve_system X0 (N x 2) and null-space basis Phi (N x k): the Julia port samples the shape "
                      "space with the C++ Phi so that its Gaussian samples are the C++ ones.";
  {
    std::ofstream o(dir + "/derivation_inputs.json");
    o << out.dump() << "\n";
  }
  {
    json l2 = freeze_l2();
    l2["provenance"] = "L2 populations of tests/derivation_tests.cpp: k7 = make_tiling_pattern(fam,n,m, "
                       "mt19937(20260904+7919*i)) for 7 families x {2x2,3x2,3x3} then 12 l2_voronoi_pattern(i, "
                       "ns[i], 10, mt19937(9100001+104729*i)); fresh110 = l2_voronoi_pattern(1000+i, U(12,45), 10, "
                       "mt19937(4400011+7717*i)); ce2000 = l2_voronoi_pattern(90000+i, U(6,16), 10, "
                       "mt19937(880011+65537*i)). Cells carry sigma from quotient_sigma; T = [[a,b],[c,d]] row-major.";
    std::ofstream o(dir + "/derivation_inputs_l2.json");
    o << l2.dump() << "\n";
  }
  std::cout << "done\n";
  return 0;
}
