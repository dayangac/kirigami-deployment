// freeze_method_1.cpp -- replays the inputs of the tests/test_method.cpp cases that
// exercise method/{deploy_basis,mobility,contact,range_opt} (the `make_case` tilings
// with checkerboard / relaxation sigma, the delaunay_of_random_points meshes of the
// mobility cases, and the range-objective case) and dumps them plus every number the
// C++ tests compute, so the Julia port runs the identical checks against C++ values.
// Linked against code/build/libkiri_core.a; same functions, same arguments as the tests.
#include <fstream>
#include <iostream>
#include <random>
#include <nlohmann/json.hpp>
#include "helpers.hpp"
#include "method/contact.hpp"
#include "method/deploy_basis.hpp"
#include "method/mobility.hpp"
#include "method/range_opt.hpp"

using namespace kiri;
using namespace kiri::method;
using nlohmann::json;

static json mesh_json(const Mesh& m) { return json::parse(mesh_to_json_string(m)); }
static json pts_json(const std::vector<Vec2>& P) {
  json a = json::array();
  for (const auto& p : P) a.push_back({p.x(), p.y()});
  return a;
}
static json mat_json(const Eigen::MatrixXd& M) {
  json a = json::array();
  for (int i = 0; i < M.rows(); ++i) {
    json r = json::array();
    for (int j = 0; j < M.cols(); ++j) r.push_back(M(i, j));
    a.push_back(r);
  }
  return a;
}
static json vec_json(const Eigen::VectorXd& v) {
  json a = json::array();
  for (int i = 0; i < v.size(); ++i) a.push_back(v(i));
  return a;
}

struct Case {
  Mesh m;
  CutStructure c;
  HoleSet hs;
  std::vector<Vec2> X;
};

// Verbatim copy of test_method.cpp's make_case.
static Case make_case(Mesh m, bool checker) {
  std::mt19937 rng(2026);
  m.build_topology();
  m.sigma = checker ? test::checkerboard_sigma(m)
                    : assign_orientation_relaxation(m, rng, 8, 500, 180).sigma;
  Case cs;
  cs.m = std::move(m);
  cs.m.build_topology();
  cs.c = make_cut(cs.m);
  cs.c.mesh = &cs.m;
  cs.hs = holes_partition(cs.c);
  cs.X = cs.m.X;
  if (!hole_residuals(cs.c, cs.X, cs.hs).deployable(1e-9)) {
    const LinearSystem sys = assemble_system(cs.c, cs.hs, cs.m.X, BoundaryMode::Fixed);
    const SolveReport sr = solve_system(sys, cs.m.X);
    cs.X = matrix_to_points(sr.X0);
  }
  return cs;
}

// Verbatim copy of test_method.cpp's bisect_theta_max.
static double bisect_theta_max(const CutStructure& c, const std::vector<Vec2>& X, double shrink) {
  auto col = [&](double th) { return has_collision(c, deploy(c, X, th).Y, shrink); };
  double lo = 0, hi = -1;
  for (int i = 1; i <= 360; ++i) {
    const double th = M_PI * i / 360;
    if (col(th)) { hi = th; break; }
    lo = th;
  }
  if (hi < 0) return M_PI;
  for (int i = 0; i < 60; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (col(mid)) hi = mid; else lo = mid;
  }
  return lo;
}

static json cert_json(const ValidityCertificate& ct) {
  return {{"pos", ct.pos}, {"nooverlap", ct.nooverlap}, {"noroot", ct.noroot},
          {"eps", ct.eps}, {"theta_1", ct.theta_1}, {"min_signed_area", ct.min_signed_area},
          {"n_inverted", ct.n_inverted}, {"n_pairs", ct.n_pairs},
          {"n_candidates", ct.n_candidates}, {"n_identically_zero", ct.n_identically_zero},
          {"n_class1", ct.n_class1}, {"n_class2", ct.n_class2}, {"n_class3", ct.n_class3},
          {"n_roots_deflated", ct.n_roots_deflated},
          {"n_roots_undeflated", ct.n_roots_undeflated},
          {"n_roots_inadmissible", ct.n_roots_inadmissible}, {"first_root", ct.first_root},
          // M'-vertex ids, 0-based as the C++ stores them (-1 = none)
          {"bad_pv", ct.bad_pv}, {"bad_a", ct.bad_a}, {"bad_b", ct.bad_b}};
}

static json mob_json(const MobilityReport& r) {
  return {{"F", r.F}, {"n_hinge", r.n_hinge}, {"n_cycles", r.n_cycles},
          {"components", r.components}, {"dim_ker_A", r.dim_ker_A},
          {"n_dangling", r.n_dangling}, {"core_faces", r.core_faces},
          {"core_edges", r.core_edges}, {"core_cycles", r.core_cycles}, {"c_core", r.c_core},
          {"dim_ker_A_core", r.dim_ker_A_core}, {"identity_holds", r.identity_holds},
          {"m_full", r.m_full}, {"m_core", r.m_core}, {"sigma_in_ker", r.sigma_in_ker},
          {"sigma_residual", r.sigma_residual}, {"used_sparse", r.used_sparse}};
}

// One tiling case: the input (mesh with sigma, X) and everything the contact tests compute.
static json freeze_tiling(const std::string& name, Case cs, const std::string& prov) {
  json j;
  j["provenance"] = prov;
  j["mesh"] = mesh_json(cs.m);
  j["sigma"] = cs.m.sigma;
  j["X"] = pts_json(cs.X);
  j["N"] = cs.m.n_vertices();
  j["F"] = cs.m.n_faces();
  j["n_prime"] = cs.c.n_prime_vertices;
  j["n_split"] = cs.c.n_split();
  j["flat_collides_1e12"] = has_collision(cs.c, deploy(cs.c, cs.X, 0.0).Y, 1e-12);
  j["flat_collides_default"] = has_collision(cs.c, deploy(cs.c, cs.X, 0.0).Y);
  const DeployBasis B = deploy_basis(cs.c, cs.X);
  j["C"] = mat_json(B.C);
  j["S"] = mat_json(B.S);
  const SweptDiscs sd = swept_discs(cs.c, B);
  j["rho_max"] = sd.rho_max;
  j["circum"] = sd.circum;
  j["rho"] = sd.rho;
  const auto all = candidate_pairs(cs.c, sd, M_PI, false);
  const auto pruned = candidate_pairs(cs.c, sd, M_PI, true);
  j["n_pairs_all"] = all.size();
  j["n_pairs_pruned"] = pruned.size();
  {
    json pp = json::array();  // 0-based face ids
    for (const auto& [f, g] : pruned) pp.push_back({f, g});
    j["pairs_pruned"] = pp;
  }
  const ExactRangeReport ea = exact_theta_max(cs.c, B, all);
  const ExactRangeReport ep = exact_theta_max(cs.c, B, pruned);
  j["theta_exact_all"] = ea.theta_max;
  j["theta_exact_pruned"] = ep.theta_max;
  j["exact_pruned"] = {{"n_pairs", ep.n_pairs}, {"n_candidates", ep.n_candidates},
                       {"n_roots", ep.n_roots}, {"n_zero_contacts", ep.n_zero_contacts},
                       {"found", ep.first.found}, {"theta", ep.first.theta},
                       {"pv", ep.first.pv}, {"pa", ep.first.pa}, {"pb", ep.first.pb},
                       {"face_v", ep.first.face_v}, {"face_e", ep.first.face_e},
                       {"corner_e", ep.first.corner_e}};
  const auto ova = exact_theta_max_overlap(cs.c, B, all, 1e-9, M_PI, 1e-9);
  const auto ovp = exact_theta_max_overlap(cs.c, B, pruned, 1e-9, M_PI, 1e-9);
  j["theta_overlap_all"] = ova.theta_max;
  j["theta_overlap_pruned"] = ovp.theta_max;
  j["overlap_pruned"] = {{"candidates", ovp.candidates}, {"i_star", ovp.i_star},
                         {"n_intervals_tested", ovp.n_intervals_tested},
                         {"zero_range", ovp.zero_range},
                         {"n_overlap_tests", ovp.n_overlap_tests},
                         {"first_contact_theta", ovp.first_contact.theta},
                         {"first_contact_found", ovp.first_contact.found}};
  j["theta_bisect_1e12"] = bisect_theta_max(cs.c, cs.X, 1e-12);
  j["theta_shipped"] = theta_max(cs.c, cs.X, 90, 40).theta_max_geometric;
  j["cert_0006_conv"] = cert_json(validity_certificate(cs.c, cs.X, 0.006));
  j["cert_0006"] = cert_json(validity_certificate(cs.c, B, cs.X, pruned, 0.006));
  const ValidityCertificate c3 = validity_certificate(cs.c, B, cs.X, pruned, 3.0);
  j["cert_3"] = cert_json(c3);
  if (!c3.noroot) {
    j["cert_3_at_root_noroot"] = validity_certificate(cs.c, B, cs.X, pruned, c3.first_root).noroot;
    j["cert_3_ulp_below_noroot"] =
        validity_certificate(cs.c, B, cs.X, pruned, std::nextafter(c3.first_root, 0.0)).noroot;
    j["cert_3_1e12_below_noroot"] =
        validity_certificate(cs.c, B, cs.X, pruned, c3.first_root - 1e-12).noroot;
  }
  j["contact_angles_0_0006"] = contact_angles(cs.c, B, pruned, 0.0, 0.006);
  return j;
}

int main(int argc, char** argv) {
  const std::string out = argc > 1 ? argv[1] : "test_fixtures_method_1.json";
  json J;
  J["_provenance"] = {
      {"producer", "scratchpad freeze_method_1.cpp (port-method-1), linked against code/build/libkiri_core.a"},
      {"replays", "tests/test_method.cpp make_case(...) inputs for the deploy_basis / contact / mobility / range_opt cases"},
      {"index_base", "0 (C++): face ids, M'-vertex ids, pair lists, bad_* fields"},
      {"toolchain", "Apple clang, arm64, -O2 -std=c++20, default -ffp-contract=on"}};

  // ---- the make_case tilings used by the contact / deploy_basis tests ----------
  struct T { const char* name; Mesh m; bool checker; const char* prov; };
  std::vector<T> tilings;
  tilings.push_back({"squares_checker", tiling_squares(rect(Vec2(2.5, 2.5), 2.51, 2.51)), true,
                     "make_case(tiling_squares(rect(Vec2(2.5,2.5),2.51,2.51)), true)"});
  tilings.push_back({"triangles_checker", tiling_triangles(disk(Vec2(0.13, 0.07), 2.5)), true,
                     "make_case(tiling_triangles(disk(Vec2(0.13,0.07),2.5)), true)"});
  tilings.push_back({"kagome_checker", tiling_kagome(disk(Vec2(0.13, 0.07), 2.5)), true,
                     "make_case(tiling_kagome(disk(Vec2(0.13,0.07),2.5)), true)"});
  tilings.push_back({"snub_2_6", tiling_snub_square(disk(Vec2(0.13, 0.07), 2.6)), false,
                     "make_case(tiling_snub_square(disk(Vec2(0.13,0.07),2.6)), false); rng mt19937(2026), relaxation(8,500,180)"});
  tilings.push_back({"snub_3_0", tiling_snub_square(disk(Vec2(0.13, 0.07), 3.0)), false,
                     "make_case(tiling_snub_square(disk(Vec2(0.13,0.07),3.0)), false); rng mt19937(2026), relaxation(8,500,180)"});
  tilings.push_back({"snub_3_2", tiling_snub_square(disk(Vec2(0.13, 0.07), 3.2)), false,
                     "make_case(tiling_snub_square(disk(Vec2(0.13,0.07),3.2)), false); rng mt19937(2026), relaxation(8,500,180)"});
  tilings.push_back({"trunc_3_0", tiling_truncated_square(disk(Vec2(0.13, 0.07), 3.0)), false,
                     "make_case(tiling_truncated_square(disk(Vec2(0.13,0.07),3.0)), false); rng mt19937(2026), relaxation(8,500,180)"});
  tilings.push_back({"trunc_4_0", tiling_truncated_square(disk(Vec2(0.13, 0.07), 4.0)), false,
                     "make_case(tiling_truncated_square(disk(Vec2(0.13,0.07),4.0)), false); rng mt19937(2026), relaxation(8,500,180)"});
  tilings.push_back({"hexagons_3_0", tiling_hexagons(disk(Vec2(0.13, 0.07), 3.0)), false,
                     "make_case(tiling_hexagons(disk(Vec2(0.13,0.07),3.0)), false); rng mt19937(2026), relaxation(8,500,180)"});
  tilings.push_back({"t3_4_3_12_4_2", tiling_3_4_3_12(disk(Vec2(0.13, 0.07), 4.2)), false,
                     "make_case(tiling_3_4_3_12(disk(Vec2(0.13,0.07),4.2)), false); rng mt19937(2026), relaxation(8,500,180)"});
  for (auto& t : tilings) {
    std::cerr << "tiling " << t.name << "\n";
    J["tilings"][t.name] = freeze_tiling(t.name, make_case(std::move(t.m), t.checker), t.prov);
  }

  // ---- mobility: "A: sigma is in its kernel exactly when Eq. (2) holds" --------
  {
    json arr = json::array();
    std::mt19937 rng(31337);
    for (int t = 0; t < 6; ++t) {
      Mesh m = delaunay_of_random_points(40 + 5 * t, 3.0, rng);
      m.build_topology();
      m.sigma = assign_orientation_relaxation(m, rng, 4, 300, 90).sigma;
      json rec;
      rec["t"] = t;
      rec["skipped"] = m.sigma.empty();
      if (m.sigma.empty()) { arr.push_back(rec); continue; }
      m.build_topology();
      const CutStructure c = make_cut(m);
      const HoleSet hs = holes_partition(c);
      const HingeGraph g = build_hinge_graph(c);
      rec["mesh"] = mesh_json(m);
      rec["sigma"] = m.sigma;
      rec["n_hinge_edges"] = g.n_edges();
      if (g.n_edges() == 0) { rec["skipped"] = true; arr.push_back(rec); continue; }
      const MobilityReport r_ini = mobility_at(c, g, pins_flat(c, g, m.X));
      const LinearSystem sys = assemble_system(c, hs, m.X, BoundaryMode::Fixed);
      const SolveReport sr = solve_system(sys, m.X);
      const auto X0 = matrix_to_points(sr.X0);
      const MobilityReport r_0 = mobility_at(c, g, pins_flat(c, g, X0));
      rec["X0"] = pts_json(X0);
      rec["r_ini"] = mob_json(r_ini);
      rec["r_0"] = mob_json(r_0);
      rec["deployable_ini"] = hole_residuals(c, m.X, hs).deployable(1e-9);
      rec["deployable_0"] = hole_residuals(c, X0, hs).deployable(1e-9);
      rec["hinge_graph"] = {{"F", g.F}, {"components", g.components}, {"n_cycles", g.n_cycles()},
                            {"head", g.head}, {"tail", g.tail}, {"eid", g.eid},
                            {"nontree", g.nontree}};
      arr.push_back(rec);
    }
    J["mobility_sigma_kernel"] = {{"provenance", "std::mt19937 rng(31337); t=0..5: delaunay_of_random_points(40+5t, 3.0, rng); sigma = assign_orientation_relaxation(m, rng, 4, 300, 90)"},
                                  {"cases", arr}};
  }
  // ---- mobility: "A and the body-and-pin rigidity matrix give the same mobility" --
  {
    json arr = json::array();
    std::mt19937 rng(4711);
    for (int t = 0; t < 5; ++t) {
      Mesh m = delaunay_of_random_points(25 + 4 * t, 3.0, rng);
      m.build_topology();
      m.sigma = assign_orientation_relaxation(m, rng, 4, 300, 90).sigma;
      json rec;
      rec["t"] = t;
      rec["skipped"] = m.sigma.empty();
      if (m.sigma.empty()) { arr.push_back(rec); continue; }
      m.build_topology();
      const CutStructure c = make_cut(m);
      const HingeGraph g = build_hinge_graph(c);
      rec["mesh"] = mesh_json(m);
      rec["sigma"] = m.sigma;
      if (g.n_edges() == 0) { rec["skipped"] = true; arr.push_back(rec); continue; }
      const auto pins = pins_flat(c, g, m.X);
      const MobilityReport r = mobility_at(c, g, pins);
      const Eigen::MatrixXd R = build_rigidity(g, pins);
      Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(R / R.cwiseAbs().maxCoeff());
      qr.setThreshold(1e-10);
      const int m_R = 3 * g.F - static_cast<int>(qr.rank()) - 3 * g.components;
      rec["r"] = mob_json(r);
      rec["rank_R"] = static_cast<int>(qr.rank());
      rec["m_R"] = m_R;
      arr.push_back(rec);
    }
    J["mobility_rigidity"] = {{"provenance", "std::mt19937 rng(4711); t=0..4: delaunay_of_random_points(25+4t, 3.0, rng); sigma = assign_orientation_relaxation(m, rng, 4, 300, 90)"},
                              {"cases", arr}};
  }
  // ---- range objective gradient case ------------------------------------------
  {
    Case cs = make_case(tiling_truncated_square(disk(Vec2(0.13, 0.07), 3.0)), false);
    const LinearSystem sys = assemble_system(cs.c, cs.hs, cs.m.X, BoundaryMode::Fixed);
    const SolveReport sr = solve_system(sys, cs.m.X);
    const std::vector<Vec2> X0 = matrix_to_points(sr.X0);
    RangeOptOptions o;
    RangeObjective ob;
    ob.setup(cs.c, X0, sr.Phi, o);
    const int k = static_cast<int>(sr.Phi.cols());
    Eigen::VectorXd t = Eigen::VectorXd::Zero(2 * k);
    ob.rebuild_active(t);
    json rec;
    rec["provenance"] = "make_case(tiling_truncated_square(disk(Vec2(0.13,0.07),3.0)), false); solve_system(assemble_system(c, hs, m.X, Fixed), m.X); RangeObjective::setup(c, X0, Phi, default opts); t ~ normal(0,1e-3) from mt19937(5)";
    rec["mesh"] = mesh_json(cs.m);
    rec["sigma"] = cs.m.sigma;
    rec["X0"] = pts_json(X0);
    rec["Phi"] = mat_json(sr.Phi);
    rec["dim_null"] = sr.dim_null;
    json act0 = json::array();
    for (const auto& a : ob.active) act0.push_back({a[0], a[1], a[2]});
    rec["active_at_0"] = act0;  // 0-based (pa, pb, pv)
    std::mt19937 rng(5);
    std::normal_distribution<double> G(0.0, 1e-3);
    for (int i = 0; i < 2 * k; ++i) t(i) = G(rng);
    rec["t"] = vec_json(t);
    Eigen::VectorXd g;
    const double f0 = ob.value_and_grad(t, g);
    rec["f0"] = f0;
    rec["g"] = vec_json(g);
    const DeployBasis Bo = ob.basis(t);
    rec["basis_C"] = mat_json(Bo.C);
    rec["basis_S"] = mat_json(Bo.S);
    std::uniform_int_distribution<int> RI(0, 2 * k - 1);
    json fds = json::array();
    for (int trial = 0; trial < 12; ++trial) {
      const int i = RI(rng);
      const double h = 1e-6;
      Eigen::VectorXd tp = t, tm = t, dummy;
      tp(i) += h;
      tm(i) -= h;
      const double fd = (ob.value_and_grad(tp, dummy) - ob.value_and_grad(tm, dummy)) / (2 * h);
      fds.push_back({{"i", i}, {"fd", fd}, {"g", g(i)}});
    }
    rec["fd_checks"] = fds;
    // a short maximize_range run, for an end-to-end regression number
    RangeOptOptions o2;
    o2.rounds = 2;
    o2.iters_per_round = 5;
    const RangeOptResult rr = maximize_range(cs.c, X0, sr.Phi, o2);
    rec["maximize_range_2x5"] = {{"theta_ref_before", rr.theta_ref_before},
                                 {"theta_ref_after", rr.theta_ref_after},
                                 {"theta_closed_before", rr.theta_closed_before},
                                 {"theta_closed_after", rr.theta_closed_after},
                                 {"rounds_run", rr.rounds_run}, {"n_active", rr.n_active},
                                 {"trace", rr.trace}, {"T", mat_json(rr.T)}};
    J["range_objective"] = rec;
  }

  std::ofstream(out) << J.dump() << "\n";
  std::cerr << "wrote " << out << "\n";
  return 0;
}
