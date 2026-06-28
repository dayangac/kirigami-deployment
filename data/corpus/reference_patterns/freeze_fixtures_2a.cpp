// freeze_fixtures_2a.cpp -- replays the mesh/sigma/RNG sequences of
// tests/test_system.cpp, tests/test_kinematics.cpp, tests/test_collision.cpp and
// tests/test_rank_checks.cpp and dumps them (plus every number the C++ tests compute
// inline) as JSON, so the Julia port can run the identical checks against the C++
// values without the (concurrently ported) generators / orientation code.
#include <fstream>
#include <iostream>
#include <numeric>
#include <nlohmann/json.hpp>
#include "helpers.hpp"

using namespace kiri;
using namespace kiri::test;
using nlohmann::json;

static json mesh_json(const Mesh& m) {
  json j = json::parse(mesh_to_json_string(m));
  // faces are stored as-is; the Julia loader must NOT re-normalise (torus faces are
  // geometrically degenerate and would flip)
  const Mesh back = mesh_from_json_string(mesh_to_json_string(m));
  j["faces_roundtrip_stable"] = (back.faces == m.faces);
  return j;
}
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

struct Analysis {
  CutStructure c;
  HoleSet hs;
  LinearSystem sys;
  SolveReport rep;
};
static Analysis analyze(const Mesh& g, BoundaryMode mode = BoundaryMode::Fixed) {
  Analysis a;
  a.c = make_cut(g);
  a.hs = holes_partition(a.c);
  a.sys = assemble_system(a.c, a.hs, g.X, mode);
  a.rep = solve_system(a.sys, g.X);
  return a;
}
static json solve_json(const Analysis& a) {
  return {{"rank_full", a.rep.rank_full}, {"rank_L", a.rep.rank_L}, {"dim_null", a.rep.dim_null},
          {"n_hole_rows", a.sys.n_hole_rows}, {"n_boundary_rows", a.sys.n_boundary_rows},
          {"n_split", a.c.n_split()}, {"n_hinge", a.c.n_hinge()},
          {"n_interior_holes", a.hs.n_interior_holes()},
          {"projection_ok", a.rep.projection_ok}, {"projection_res", a.rep.projection_res},
          {"sv_tol", a.rep.sv_tol}, {"X0", mat_json(a.rep.X0)}};
}

// ------------------------------------------------------------------ test_system.cpp
static json freeze_system() {
  json out;
  {
    Mesh g = tiling_squares(rect(Vec2(2.5, 2.5), 2.51, 2.51));
    g.sigma = checkerboard_sigma(g);
    const Analysis a = analyze(g);
    out["rotating_squares"] = {{"mesh", mesh_json(g)}, {"solve", solve_json(a)},
                               {"components", count_components(a.c)},
                               {"deploy_0_6_max_mismatch", deploy(a.c, g.X, 0.6).max_mismatch},
                               {"provenance", "tiling_squares(rect(Vec2(2.5,2.5),2.51,2.51)); sigma = checkerboard_sigma(g)"}};
  }
  {
    Mesh g = tiling_triangles(disk(Vec2(0.13, 0.07), 3.5));
    g.sigma = checkerboard_sigma(g);
    out["triangles_alternating"] = {{"mesh", mesh_json(g)}, {"solve", solve_json(analyze(g))},
                                    {"provenance", "tiling_triangles(disk(Vec2(0.13,0.07),3.5)); sigma = checkerboard_sigma(g)"}};
  }
  {
    Mesh g = tiling_kagome(disk(Vec2(0.13, 0.07), 3.0));
    g.sigma = checkerboard_sigma(g);
    out["kagome_alternating"] = {{"mesh", mesh_json(g)}, {"solve", solve_json(analyze(g))},
                                 {"provenance", "tiling_kagome(disk(Vec2(0.13,0.07),3.0)); sigma = checkerboard_sigma(g)"}};
  }
  {
    Mesh g = tiling_hexagons(disk(Vec2(0.13, 0.07), 2.4));
    const int F = g.n_faces();
    int connected = 0, violations = 0, with_split = 0;
    for (unsigned long long mask = 0; mask < (1ull << F); ++mask) {
      if (mask & 1ull) continue;
      std::vector<int> sig(F);
      for (int i = 0; i < F; ++i) sig[i] = ((mask >> i) & 1ull) ? 1 : -1;
      Mesh h = g;
      h.sigma = sig;
      const CutStructure c = make_cut(h);
      if (count_components(c) != 1) continue;
      ++connected;
      if (c.n_split() > 0) ++with_split;
      const HoleSet hs = holes_partition(c);
      const LinearSystem sys = assemble_system(c, hs, h.X, BoundaryMode::Fixed);
      const SolveReport rep = solve_system(sys, h.X);
      if (rep.dim_null != h.n_interior_vertices() - hs.n_interior_holes()) ++violations;
    }
    out["hexagons_exhaustive"] = {{"mesh", mesh_json(g)}, {"F", F}, {"connected", connected},
                                  {"with_split", with_split}, {"violations", violations},
                                  {"provenance", "tiling_hexagons(disk(Vec2(0.13,0.07),2.4)); all 2^F sigma with bit 0 == -1"}};
  }
  {
    std::mt19937 rng(3);
    Mesh g = tiling_3_4_3_12(disk(Vec2(0.13, 0.07), 4.2));
    g.sigma = assign_orientation_relaxation(g, rng, 6, 400, 120).sigma;
    const Analysis a = analyze(g);
    const Residuals r = hole_residuals(a.c, g.X, a.hs);
    int violating = 0;
    for (const auto& v : r.per_hole)
      if (v.norm() > 1e-6) ++violating;
    const auto X0 = matrix_to_points(a.rep.X0);
    out["t3_4_3_12"] = {{"mesh", mesh_json(g)}, {"solve", solve_json(a)},
                        {"violating", violating}, {"n_per_hole", r.per_hole.size()},
                        {"max_norm", r.max_norm}, {"l2_norm", r.l2_norm},
                        {"per_hole", pts_json(r.per_hole)},
                        {"X0_max_norm", hole_residuals(a.c, X0, a.hs).max_norm},
                        {"X0_deploy_0_4_mismatch", deploy(a.c, X0, 0.4).max_mismatch},
                        {"provenance", "mt19937 rng(3); tiling_3_4_3_12(disk(Vec2(0.13,0.07),4.2)); sigma = assign_orientation_relaxation(g, rng, 6, 400, 120).sigma"}};
  }
  {
    std::mt19937 rng(19);
    std::uniform_real_distribution<double> U(-1.5, 1.5);
    Mesh g = tiling_kagome(disk(Vec2(0.13, 0.07), 3.0));
    g.sigma = checkerboard_sigma(g);
    const Analysis a = analyze(g);
    const auto X0 = matrix_to_points(a.rep.X0);
    json maps = json::array();
    double worst = 0;
    for (int t = 0; t < 100; ++t) {
      Eigen::Matrix2d A;
      A << U(rng), U(rng), U(rng), U(rng);
      const Vec2 b(U(rng), U(rng));
      maps.push_back({A(0, 0), A(0, 1), A(1, 0), A(1, 1), b.x(), b.y()});
      std::vector<Vec2> Y(X0.size());
      double scale = 0;
      for (size_t i = 0; i < X0.size(); ++i) {
        Y[i] = A * X0[i] + b;
        scale = std::max(scale, Y[i].norm());
      }
      worst = std::max(worst, hole_residuals(a.c, Y, a.hs).max_norm / std::max(scale, 1e-12));
    }
    out["affine"] = {{"maps", maps}, {"worst", worst}, {"X0_max_norm", hole_residuals(a.c, X0, a.hs).max_norm},
                     {"provenance", "kagome_alternating mesh; mt19937 rng(19), U(-1.5,1.5): per map A(0,0),A(0,1),A(1,0),A(1,1),b.x,b.y in draw order; map = [a00,a01,a10,a11,bx,by]"}};
  }
  {
    std::mt19937 rng(23);
    std::uniform_real_distribution<double> U(-4.0, 4.0);
    Mesh g = tiling_hexagons(disk(Vec2(0.13, 0.07), 3.0));
    g.sigma = assign_orientation_relaxation(g, rng, 5, 300, 120).sigma;
    const Analysis a = analyze(g);
    const auto X0 = matrix_to_points(a.rep.X0);
    json offsets = json::array();  // [i][t] = [ux, uy]
    double worst = 0;
    for (int i = 0; i < a.rep.dim_null; ++i) {
      json row = json::array();
      for (int t = 0; t < 100; ++t) {
        Eigen::MatrixXd X = a.rep.X0;
        const double ux = U(rng), uy = U(rng);
        row.push_back({ux, uy});
        X.col(0) += a.rep.Phi.col(i) * ux;
        X.col(1) += a.rep.Phi.col(i) * uy;
        const auto Y = matrix_to_points(X);
        worst = std::max(worst, hole_residuals(a.c, Y, a.hs).max_norm);
        for (int v = 0; v < g.n_vertices(); ++v)
          if (g.vertex_is_boundary[v]) worst = std::max(worst, (Y[v] - X0[v]).norm());
      }
      offsets.push_back(row);
    }
    out["nullspace"] = {{"mesh", mesh_json(g)}, {"solve", solve_json(a)}, {"offsets", offsets},
                        {"worst", worst}, {"Phi", mat_json(a.rep.Phi)},
                        {"provenance", "mt19937 rng(23); tiling_hexagons(disk(Vec2(0.13,0.07),3.0)); sigma = assign_orientation_relaxation(g, rng, 5, 300, 120).sigma; then U(-4,4) draws ux,uy per (i,t)"}};
  }
  {
    Mesh g = periodic_squares(4, 4);
    g.sigma = checkerboard_sigma(g);
    const Analysis a = analyze(g, BoundaryMode::Periodic);
    out["periodic_squares_4x4"] = {{"mesh", mesh_json(g)}, {"solve", solve_json(a)},
                                   {"X0_max_norm", hole_residuals(a.c, matrix_to_points(a.rep.X0), a.hs).max_norm},
                                   {"provenance", "periodic_squares(4,4); sigma = checkerboard_sigma(g); BoundaryMode::Periodic"}};
  }
  {
    std::mt19937 rng(101);
    Mesh g = tiling_hexagons(disk(Vec2(0.13, 0.07), 2.4));
    const OrientationReport bf = brute_force_orientation(g);
    const OrientationReport re = assign_orientation_relaxation(g, rng, 10, 600, 180);
    auto rep_json = [](const OrientationReport& r) {
      return json{{"sigma", r.sigma}, {"n_split", r.n_split}, {"n_hinge", r.n_hinge},
                  {"components", r.components}, {"n_holes", r.n_holes}, {"objective", r.objective}};
    };
    out["orientation_small_patch"] = {{"mesh", mesh_json(g)}, {"brute_force", rep_json(bf)},
                                      {"relaxation", rep_json(re)},
                                      {"provenance", "mt19937 rng(101); tiling_hexagons(disk(Vec2(0.13,0.07),2.4)); brute_force_orientation(g); assign_orientation_relaxation(g, rng, 10, 600, 180)"}};
  }
  {
    Mesh g = periodic_kagome(3, 3);
    g.sigma = checkerboard_sigma(g);
    out["periodic_kagome_3x3"] = {{"mesh", mesh_json(g)}, {"N", g.n_vertices()}, {"F", g.n_faces()},
                                  {"E", g.n_edges()}, {"n_pairs_h", g.periodic.pairs_h.size()},
                                  {"provenance", "periodic_kagome(3,3); sigma = checkerboard_sigma(g)"}};
  }
  return out;
}

// -------------------------------------------------------------- test_kinematics.cpp
static std::vector<Vec2> deployable_embedding(const Mesh& g, const CutStructure& c, const HoleSet& hs) {
  const LinearSystem sys = assemble_system(c, hs, g.X, BoundaryMode::Fixed);
  const SolveReport rep = solve_system(sys, g.X);
  return matrix_to_points(rep.X0);
}

static json freeze_kinematics() {
  json out;
  {
    std::mt19937 rng(5);
    json fam = json::array();
    for (const std::string& kind : {"squares", "kagome", "triangles", "hexagons", "snub_square"}) {
      Mesh g = generate(kind, {3.0}, rng);
      g.sigma = assign_orientation_relaxation(g, rng, 4, 250, 90).sigma;
      const CutStructure c = make_cut(g);
      const HoleSet hs = holes_partition(c);
      const auto X = deployable_embedding(g, c, hs);
      json mm = json::array();
      for (double th : {0.1, 0.5, 1.0}) mm.push_back(deploy(c, X, th).max_mismatch);
      fam.push_back({{"kind", kind}, {"mesh", mesh_json(g)}, {"X", pts_json(X)},
                     {"X_max_norm", hole_residuals(c, X, hs).max_norm}, {"max_mismatch_0.1_0.5_1.0", mm}});
    }
    out["hinge_opens"] = {{"families", fam},
                          {"provenance", "mt19937 rng(5); per kind: generate(kind,{3.0},rng); sigma = assign_orientation_relaxation(g,rng,4,250,90).sigma; X = Eq.(6) projection (Fixed)"}};
  }
  {
    std::mt19937 rng(11);
    json fam = json::array();
    for (const std::string& kind : {"squares", "kagome", "hexagons", "truncated_square"}) {
      Mesh g = generate(kind, {3.0}, rng);
      g.sigma = assign_orientation_relaxation(g, rng, 4, 250, 90).sigma;
      const CutStructure c = make_cut(g);
      const HoleSet hs = holes_partition(c);
      const auto X = deployable_embedding(g, c, hs);
      const Deployment ref = deploy(c, X, 0.7, 0);
      json orders = json::array(), resid = json::array();
      for (int trial = 0; trial < 8; ++trial) {
        std::vector<int> order(g.n_faces());
        std::iota(order.begin(), order.end(), 0);
        std::shuffle(order.begin(), order.end(), rng);
        const Deployment alt = deploy_with_order(c, X, 0.7, order);
        orders.push_back(order);
        resid.push_back({alt.max_mismatch, rigid_align_residual(alt.Y, ref.Y)});
      }
      fam.push_back({{"kind", kind}, {"mesh", mesh_json(g)}, {"X", pts_json(X)},
                     {"orders", orders}, {"Y_ref", pts_json(ref.Y)}, {"mismatch_and_align", resid}});
    }
    out["bfs_order"] = {{"families", fam},
                        {"provenance", "mt19937 rng(11); per kind: generate(kind,{3.0},rng); sigma = assign_orientation_relaxation(g,rng,4,250,90).sigma; 8 x std::shuffle(0..F-1, rng) (orders 0-based)"}};
  }
  {
    std::mt19937 rng(3);
    Mesh g = tiling_3_4_3_12(disk(Vec2(0.13, 0.07), 3.5));
    g.sigma = assign_orientation_relaxation(g, rng, 6, 400, 120).sigma;
    const CutStructure c = make_cut(g);
    const HoleSet hs = holes_partition(c);
    const Residuals r = hole_residuals(c, g.X, hs);
    const auto X0 = deployable_embedding(g, c, hs);
    out["t3_4_3_12_r3_5"] = {{"mesh", mesh_json(g)}, {"max_norm", r.max_norm},
                             {"deploy_0_3_mismatch", deploy(c, g.X, 0.3).max_mismatch},
                             {"X0", pts_json(X0)}, {"X0_max_norm", hole_residuals(c, X0, hs).max_norm},
                             {"X0_deploy_0_3_mismatch", deploy(c, X0, 0.3).max_mismatch},
                             {"provenance", "mt19937 rng(3); tiling_3_4_3_12(disk(Vec2(0.13,0.07),3.5)); sigma = assign_orientation_relaxation(g,rng,6,400,120).sigma"}};
  }
  {
    std::mt19937 rng(77);
    std::uniform_real_distribution<double> TH(0.05, 1.4);
    int samples = 0;
    double worst_rel = 0;
    json fam = json::array();
    const std::vector<std::string> kinds = {"squares", "kagome", "hexagons", "triangles", "snub_square", "truncated_square"};
    for (const auto& kind : kinds) {
      Mesh g = generate(kind, {2.6}, rng);
      json trials = json::array();
      for (int trial = 0; trial < 6; ++trial) {
        g.sigma = (trial == 0) ? checkerboard_sigma(g) : random_sigma(g, rng);
        const CutStructure c = make_cut(g);
        json ths = json::array();
        for (int s = 0; s < 6; ++s) {
          const double th = TH(rng);
          ths.push_back(th);
          const double h = 1e-6;
          const Deployment d = deploy(c, g.X, th);
          const Deployment dp = deploy(c, g.X, th + h);
          const Deployment dm = deploy(c, g.X, th - h);
          for (int i = 0; i < c.n_prime_vertices; ++i) {
            const Vec2 fd = (dp.Y[i] - dm.Y[i]) / (2 * h);
            const double denom = std::max(1.0, fd.norm());
            worst_rel = std::max(worst_rel, (fd - d.dY_dtheta[i]).norm() / denom);
            ++samples;
          }
        }
        trials.push_back({{"sigma", g.sigma}, {"thetas", ths}});
      }
      fam.push_back({{"kind", kind}, {"mesh", mesh_json(g)}, {"trials", trials}});
    }
    out["fd"] = {{"families", fam}, {"samples", samples}, {"worst_rel", worst_rel},
                 {"provenance", "mt19937 rng(77), TH = U(0.05,1.4); per kind: generate(kind,{2.6},rng); 6 trials: sigma = checkerboard (trial 0) else random_sigma(g,rng); 6 x TH(rng)"}};
  }
  {
    std::mt19937 rng(13);
    int edges_checked = 0;
    json fam = json::array();
    for (const std::string& kind : {"hexagons", "truncated_square", "snub_square", "t3_4_3_12"}) {
      Mesh g = generate(kind, {3.0}, rng);
      g.sigma = assign_orientation_relaxation(g, rng, 5, 300, 120).sigma;
      const CutStructure c = make_cut(g);
      json rec = {{"kind", kind}, {"mesh", mesh_json(g)}, {"n_split", c.n_split()}};
      if (c.n_split() > 0) {
        const HoleSet hs = holes_partition(c);
        const auto X = deployable_embedding(g, c, hs);
        const double tmax_geo = theta_max(c, X).theta_max_geometric;
        const double tmax = std::max(0.4, tmax_geo);
        for (double th : {0.1, 0.3, 0.6, 0.9 * tmax}) {
          (void)th;
          edges_checked += c.n_split();
        }
        rec["X"] = pts_json(X);
        rec["theta_max_geometric"] = tmax_geo;
      }
      fam.push_back(rec);
    }
    out["remark_A4"] = {{"families", fam}, {"edges_checked", edges_checked},
                        {"provenance", "mt19937 rng(13); per kind: generate(kind,{3.0},rng); sigma = assign_orientation_relaxation(g,rng,5,300,120).sigma"}};
  }
  return out;
}

// --------------------------------------------------------------- test_collision.cpp
static json freeze_collision() {
  json out;
  {
    std::mt19937 rng(41);
    json fam = json::array();
    const std::vector<std::string> kinds = {"squares", "kagome", "triangles", "hexagons", "snub_square", "truncated_square", "t3_4_3_12"};
    for (const auto& kind : kinds) {
      Mesh g = generate(kind, {2.8}, rng);
      g.sigma = assign_orientation_relaxation(g, rng, 5, 300, 120).sigma;
      const CutStructure c = make_cut(g);
      const HoleSet hs = holes_partition(c);
      const LinearSystem sys = assemble_system(c, hs, g.X, BoundaryMode::Fixed);
      const SolveReport rep = solve_system(sys, g.X);
      const auto X0 = matrix_to_points(rep.X0);
      json rec = {{"kind", kind}, {"mesh", mesh_json(g)}, {"X0", pts_json(X0)},
                  {"X0_max_norm", hole_residuals(c, X0, hs).max_norm}, {"n_split", c.n_split()}};
      if (hole_residuals(c, X0, hs).max_norm <= 1e-9) {
        const ThetaMaxReport tm = theta_max(c, X0, 90, 40);
        rec["theta_max_geometric"] = tm.theta_max_geometric;
        rec["min_beta"] = tm.min_beta;
        rec["collided"] = tm.collided;
        rec["beta"] = tm.beta;
      }
      fam.push_back(rec);
    }
    out["theta_max"] = {{"families", fam},
                        {"provenance", "mt19937 rng(41); per kind: generate(kind,{2.8},rng); sigma = assign_orientation_relaxation(g,rng,5,300,120).sigma; X0 = Eq.(6) projection; theta_max(c,X0,90,40)"}};
  }
  {
    Mesh g = tiling_squares(rect(Vec2(2.5, 2.5), 2.51, 2.51));
    g.sigma = checkerboard_sigma(g);
    const CutStructure c = make_cut(g);
    const ThetaMaxReport tm = theta_max(c, g.X, 180, 45);
    out["rotating_squares"] = {{"mesh", mesh_json(g)}, {"min_beta", tm.min_beta},
                               {"theta_max_geometric", tm.theta_max_geometric}, {"collided", tm.collided},
                               {"provenance", "tiling_squares(rect(Vec2(2.5,2.5),2.51,2.51)); checkerboard; theta_max(c,g.X,180,45)"}};
  }
  {
    std::mt19937 rng(53);
    Mesh g = generate("truncated_square", {3.0}, rng);
    g.sigma = assign_orientation_relaxation(g, rng, 4, 250, 90).sigma;
    const CutStructure c = make_cut(g);
    const auto z = deployment_velocity(c, g.X);
    const Deployment d = deploy(c, g.X, 0.0);
    double worst = 0;
    for (size_t i = 0; i < z.size(); ++i) worst = std::max(worst, (z[i] - d.dY_dtheta[i]).norm());
    out["velocity"] = {{"mesh", mesh_json(g)}, {"z", pts_json(z)}, {"worst", worst},
                       {"provenance", "mt19937 rng(53); generate(\"truncated_square\",{3.0},rng); sigma = assign_orientation_relaxation(g,rng,4,250,90).sigma"}};
  }
  {
    std::mt19937 rng(67);
    json fam = json::array();
    for (const std::string& kind : {"truncated_square", "snub_square", "hexagons"}) {
      Mesh g = generate(kind, {3.0}, rng);
      g.sigma = assign_orientation_relaxation(g, rng, 5, 300, 120).sigma;
      const CutStructure c = make_cut(g);
      json rec = {{"kind", kind}, {"mesh", mesh_json(g)}, {"n_split", c.n_split()}};
      if (c.n_split() > 0) {
        const HoleSet hs = holes_partition(c);
        const LinearSystem sys = assemble_system(c, hs, g.X, BoundaryMode::Fixed);
        const SolveReport rep = solve_system(sys, g.X);
        rec["dim_null"] = rep.dim_null;
        if (rep.dim_null > 0) {
          const auto X0 = matrix_to_points(rep.X0);
          rec["X0"] = pts_json(X0);
          rec["Phi"] = mat_json(rep.Phi);
          rec["X0_max_norm"] = hole_residuals(c, X0, hs).max_norm;
          if (hole_residuals(c, X0, hs).max_norm <= 1e-9) {
            const CollisionSweepResult r = optimize_collision_sweep(c, X0, rep.Phi);
            json ladder = json::array();
            for (const auto& [gm, tm] : r.ladder) ladder.push_back({gm, tm});
            rec["sweep"] = {{"theta_max_before", r.theta_max_before}, {"theta_max_after", r.theta_max_after},
                            {"gamma_used", std::isinf(r.gamma_used) ? json(nullptr) : json(r.gamma_used)},
                            {"ladder", ladder}, {"X_opt_max_norm", hole_residuals(c, r.X_opt, hs).max_norm},
                            {"iterations", r.iterations}, {"f0", r.f0}, {"f1", r.f1}};
          }
        }
      }
      fam.push_back(rec);
    }
    out["collision_opt"] = {{"families", fam},
                            {"provenance", "mt19937 rng(67); per kind: generate(kind,{3.0},rng); sigma = assign_orientation_relaxation(g,rng,5,300,120).sigma; optimize_collision_sweep(c,X0,Phi) default gammas, kappa 0.05"}};
  }
  return out;
}

// ------------------------------------------------------------- test_rank_checks.cpp
static json freeze_rank_checks() {
  json out;
  json cases = json::array();
  auto add = [&](const std::string& name, const Mesh& g, BoundaryMode mode, bool boundary_free) {
    const CutStructure c = make_cut(g);
    const HoleSet hs = holes_partition(c);
    const LinearSystem sys = assemble_system(c, hs, g.X, mode);
    const RowSumReport r = check_row_sum(c, hs, sys);
    const FactorizationReport f = check_factorization(c, hs, sys, true);
    const HingeGraphReport h = check_hinge_graph(c, hs);
    const SolveReport rep = solve_system(sys, g.X);
    std::vector<double> rv(r.r.data(), r.r.data() + r.r.size());
    cases.push_back({{"name", name}, {"mesh", mesh_json(g)},
                     {"mode", mode == BoundaryMode::Fixed ? "Fixed" : mode == BoundaryMode::Periodic ? "Periodic" : "None"},
                     {"boundary_free", boundary_free},
                     {"row_sum", {{"r", rv}, {"N", r.N}, {"H", r.H}, {"all_vertices_interior", r.all_vertices_interior},
                                  {"degree_identity", r.degree_identity}, {"n_degree_mismatch", r.n_degree_mismatch},
                                  {"restricted_degree_identity", r.restricted_degree_identity},
                                  {"n_restricted_mismatch", r.n_restricted_mismatch}, {"n_hinge", r.n_hinge},
                                  {"n_hinge_in_L", r.n_hinge_in_L}, {"no_notch_hinges", r.no_notch_hinges},
                                  {"support_on_boundary", r.support_on_boundary}, {"n_interior_nonzero", r.n_interior_nonzero},
                                  {"max_abs_r", r.max_abs_r}, {"r_is_zero", r.r_is_zero}}},
                     {"factorization", {{"H", f.H}, {"n_hinge", f.n_hinge}, {"N", f.N}, {"L_equals_RD", f.L_equals_RD},
                                        {"max_abs_diff", f.max_abs_diff}, {"Z_computed", f.Z_computed}, {"rank_L", f.rank_L},
                                        {"dim_Z", f.dim_Z}, {"rank_identity", f.rank_identity},
                                        {"n_harmonic_checks", f.n_harmonic_checks}, {"n_harmonic_violations", f.n_harmonic_violations},
                                        {"max_harmonic_residual", f.max_harmonic_residual}, {"out_harmonic", f.out_harmonic}}},
                     {"hinge_graph", {{"n_faces", h.n_faces}, {"n_hinge", h.n_hinge}, {"c_gamma", h.c_gamma}, {"H", h.H},
                                      {"H_all", h.H_all}, {"n_notches", h.n_notches}, {"predicted", h.predicted},
                                      {"identity_holds", h.identity_holds}, {"identity_with_notches", h.identity_with_notches},
                                      {"defect", h.defect}}},
                     {"solve_rank_L", rep.rank_L}, {"n_boundary_rows", sys.n_boundary_rows},
                     {"n_interior_holes", hs.n_interior_holes()}, {"n_split", c.n_split()}});
  };
  {
    std::mt19937 rng(20260903);
    { Mesh g = tiling_squares(rect(Vec2(2.5, 2.5), 2.51, 2.51)); g.sigma = checkerboard_sigma(g); add("rotating_squares", g, BoundaryMode::Fixed, false); }
    { Mesh g = tiling_triangles(disk(Vec2(0.13, 0.07), 3.5)); g.sigma = checkerboard_sigma(g); add("triangles_alternating", g, BoundaryMode::Fixed, false); }
    { Mesh g = tiling_kagome(disk(Vec2(0.13, 0.07), 3.0)); g.sigma = checkerboard_sigma(g); add("kagome_3636", g, BoundaryMode::Fixed, false); }
    { Mesh g = tiling_hexagons(disk(Vec2(0.13, 0.07), 3.0)); g.sigma = assign_orientation_relaxation(g, rng, 8, 500, 180).sigma; add("hexagons_auto", g, BoundaryMode::Fixed, false); }
    { Mesh g = tiling_truncated_square(disk(Vec2(0.13, 0.07), 4.0)); g.sigma = assign_orientation_relaxation(g, rng, 8, 500, 180).sigma; add("truncated_square_488", g, BoundaryMode::Fixed, false); }
    { Mesh g = tiling_snub_square(disk(Vec2(0.13, 0.07), 3.2)); g.sigma = assign_orientation_relaxation(g, rng, 8, 500, 180).sigma; add("snub_square_33434", g, BoundaryMode::Fixed, false); }
    { Mesh g = tiling_3_4_3_12(disk(Vec2(0.13, 0.07), 4.2)); g.sigma = assign_orientation_relaxation(g, rng, 8, 500, 180).sigma; add("tiling_3_4_3_12", g, BoundaryMode::Fixed, false); }
    { Mesh g = periodic_squares(4, 4); g.sigma = checkerboard_sigma(g); add("periodic_squares_4x4", g, BoundaryMode::Periodic, false); }
    for (int k = 0; k < 4; ++k) {
      Mesh g = delaunay_of_random_points(70 + 20 * k, 12.0, rng);
      g.sigma = random_sigma(g, rng);
      add("delaunay_random_sigma_" + std::to_string(k), g, BoundaryMode::Fixed, false);
    }
  }
  out["bounded_cases"] = cases;
  cases = json::array();
  { Mesh g = torus_squares(4, 4); g.sigma = checkerboard_sigma(g); add("torus_squares_4x4", g, BoundaryMode::None, true); }
  { Mesh g = torus_squares(6, 4); g.sigma = checkerboard_sigma(g); add("torus_squares_6x4", g, BoundaryMode::None, true); }
  { Mesh g = torus_triangles(4, 4); g.sigma = checkerboard_sigma(g); add("torus_triangles_4x4", g, BoundaryMode::None, true); }
  out["torus_cases"] = cases;
  out["provenance"] = "bounded_cases(): one mt19937 rng(20260903) shared in file order over the four assign_orientation_relaxation(g,rng,8,500,180) calls and the 4 delaunay_of_random_points(70+20k,12.0,rng) + random_sigma(g,rng) pairs; torus_cases(): torus_squares(4,4), torus_squares(6,4), torus_triangles(4,4) with checkerboard_sigma. Torus faces are geometrically degenerate: load faces AS STORED (no CCW normalisation).";
  return out;
}

int main(int argc, char** argv) {
  const std::string dir = argc > 1 ? argv[1] : ".";
  auto dump = [&](const std::string& name, const json& j) {
    std::ofstream o(dir + "/test_fixtures_" + name + ".json");
    o << j.dump() << "\n";
    std::cout << "wrote test_fixtures_" << name << ".json\n";
  };
  dump("system", freeze_system());
  dump("kinematics", freeze_kinematics());
  dump("collision", freeze_collision());
  dump("rank_checks", freeze_rank_checks());
  return 0;
}
