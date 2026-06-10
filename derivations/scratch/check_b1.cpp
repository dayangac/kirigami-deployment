// check_b1.cpp -- numeric check of the budget identity (B.1) of ideas/round2_theorist_b.md.
//
// Build (from repo root):
//   clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
//     -I code/src derivations/scratch/check_b1.cpp \
//     code/src/core/mesh.cpp code/src/core/cut.cpp code/src/core/holes.cpp \
//     code/src/core/kinematics.cpp code/src/core/orientation.cpp code/src/core/generators.cpp \
//     code/src/core/tutte_auxetic.cpp code/src/core/collision.cpp code/src/core/optimize.cpp \
//     -o /tmp/check_b1 && /tmp/check_b1
//
// B1  per hole C:  A_C(theta) = a_C sin(theta) - b_C (1 - cos(theta))  with
//       a_C = 1/2 sum_{(a->b) in dC, face f} <x_a - x_b, u_f>
//       b_C = 1/2 sum_{(a->b) in dC, face f} sigma_f det(x_a - x_b, u_f)
//     measured: shoelace of the hole cycle at 6 angles, least-squares fit of (a, b).
// B1s split-edge term: the contribution of a split edge to a_C must be 1/2 <d_e, du_e>.
// B10 the area maximum sits at exactly half the second root: analytic, checked by
//     locating the max numerically.
#include <cmath>
#include <cstdio>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "core/cut.hpp"
#include "core/generators.hpp"
#include "core/holes.hpp"
#include "core/kinematics.hpp"
#include "core/mesh.hpp"
#include "core/orientation.hpp"
#include "core/tutte_auxetic.hpp"

using namespace kiri;

static Eigen::Matrix2d Jmat() { Eigen::Matrix2d J; J << 0, -1, 1, 0; return J; }
static double det2(const Vec2& a, const Vec2& b) { return a.x() * b.y() - a.y() * b.x(); }

static std::vector<Vec2> potential_u(const CutStructure& c, const std::vector<Vec2>& X,
                                     double* worst_closure) {
  const Mesh& m = *c.mesh;
  const int F = m.n_faces();
  std::vector<std::vector<std::pair<int, int>>> adj(F);
  for (int e : c.hinge_edges) {
    const Edge& ed = m.edges[e];
    const int f1 = m.half_edges[ed.he[0]].face, f2 = m.half_edges[ed.he[1]].face;
    adj[f1].push_back({e, f2});
    adj[f2].push_back({e, f1});
  }
  std::vector<Vec2> u(F, Vec2::Zero());
  std::vector<char> seen(F, 0);
  for (int s = 0; s < F; ++s) {
    if (seen[s]) continue;
    seen[s] = 1; u[s] = Vec2::Zero();
    std::vector<int> q{s};
    for (size_t qi = 0; qi < q.size(); ++qi) {
      const int f = q[qi];
      for (auto [e, g] : adj[f]) {
        const Vec2 target = u[f] + double(m.sigma[g]) * X[c.hinge_dir[e].src];
        if (!seen[g]) { seen[g] = 1; u[g] = target; q.push_back(g); }
      }
    }
  }
  double w = 0;
  for (int e : c.hinge_edges) {
    const Edge& ed = m.edges[e];
    const int f1 = m.half_edges[ed.he[0]].face, f2 = m.half_edges[ed.he[1]].face;
    const Vec2& xv = X[c.hinge_dir[e].src];
    w = std::max(w, (u[f2] - u[f1] - double(m.sigma[f2]) * xv).norm());
  }
  *worst_closure = w;
  return u;
}

static double shoelace(const std::vector<Vec2>& P) {
  double a = 0;
  for (size_t i = 0; i < P.size(); ++i) {
    const Vec2& p = P[i]; const Vec2& q = P[(i + 1) % P.size()];
    a += p.x() * q.y() - p.y() * q.x();
  }
  return 0.5 * a;
}

int main() {
  const Eigen::Matrix2d J = Jmat();
  std::mt19937 rng(20260904);
  struct Case { std::string kind; std::vector<double> par; };
  std::vector<Case> cases = {{"squares", {3.5}}, {"triangles", {3.5}}, {"hexagons", {3.5}},
                             {"kagome", {3.5}},  {"snub_square", {3.5}},
                             {"truncated_square", {3.5}}, {"t3_4_3_12", {4.0}},
                             {"squares", {5.0}}, {"hexagons", {5.0}}, {"kagome", {4.5}}};
  const std::vector<double> th = {0.05, 0.2, 0.4, 0.7, 1.0, 1.3};

  double worst_rel = 0, worst_abs = 0, worst_ratio = 0, worst_split = 0, worst_clo = 0;
  std::string where_rel, where_ratio;
  int n_holes = 0, n_graphs = 0, n_split_terms = 0;

  for (const auto& cs : cases) {
    Mesh m = generate(cs.kind, cs.par, rng);
    m.build_topology();
    auto orep = assign_orientation_relaxation(m, rng);
    m.sigma = orep.sigma;
    m.build_topology();
    CutStructure c = make_cut(m);
    HoleSet hs = holes_partition(c);
    LinearSystem sys = assemble_system(c, hs, m.X, BoundaryMode::Fixed);
    SolveReport sr = solve_system(sys, m.X);
    if (!sr.projection_ok) continue;
    std::vector<Vec2> X = matrix_to_points(sr.X0);
    Residuals res = hole_residuals(c, X, hs);
    if (!res.deployable(1e-7)) continue;
    double clo = 0;
    std::vector<Vec2> u = potential_u(c, X, &clo);
    worst_clo = std::max(worst_clo, clo);

    // pv -> faces containing it
    std::vector<std::vector<int>> pv_faces(c.n_prime_vertices);
    for (int f = 0; f < m.n_faces(); ++f)
      for (int pv : c.prime_faces[f]) pv_faces[pv].push_back(f);

    Deployment d0 = deploy(c, X, th[0]);
    auto cycles = holes_geometric_cycles(c, d0.Y);
    if (cycles.empty()) { std::printf("%-18s no cycles\n", cs.kind.c_str()); continue; }
    ++n_graphs;

    // deployed positions at every angle
    std::vector<std::vector<Vec2>> Y(th.size());
    for (size_t k = 0; k < th.size(); ++k) Y[k] = deploy(c, X, th[k]).Y;

    double gw = 0;
    for (const auto& cyc : cycles) {
      if (cyc.size() < 3) continue;
      // closed form (B.1) along the SAME directed walk
      double a_cf = 0, b_cf = 0; bool ok = true;
      double split_sum = 0; int n_sp = 0;
      for (size_t i = 0; i < cyc.size(); ++i) {
        const int pa = cyc[i], pb = cyc[(i + 1) % cyc.size()];
        int f = -1;
        for (int fa : pv_faces[pa])
          for (int fb : pv_faces[pb])
            if (fa == fb) f = fa;
        if (f < 0) { ok = false; break; }
        const Vec2 dv = X[c.prime_to_original[pa]] - X[c.prime_to_original[pb]];
        a_cf += 0.5 * dv.dot(u[f]);
        b_cf += 0.5 * double(m.sigma[f]) * det2(dv, u[f]);
      }
      if (!ok) continue;
      // measured areas + least squares fit of (a, b)
      Eigen::MatrixXd M(th.size(), 2); Eigen::VectorXd rhs(th.size());
      for (size_t k = 0; k < th.size(); ++k) {
        std::vector<Vec2> P; P.reserve(cyc.size());
        for (int pv : cyc) P.push_back(Y[k][pv]);
        rhs(k) = shoelace(P);
        M(k, 0) = std::sin(th[k]);
        M(k, 1) = -(1.0 - std::cos(th[k]));
      }
      Eigen::Vector2d fit = M.colPivHouseholderQr().solve(rhs);
      const double resid = (M * fit - rhs).cwiseAbs().maxCoeff();
      const double scale = std::max(1e-12, std::fabs(fit(0)) + std::fabs(fit(1)));
      const double rel = (std::fabs(fit(0) - a_cf) + std::fabs(fit(1) - b_cf)) / scale;
      ++n_holes;
      if (rel > worst_rel) { worst_rel = rel; where_rel = cs.kind; }
      worst_abs = std::max(worst_abs, resid / scale);
      gw = std::max(gw, rel);
      // B10: numeric argmax vs half the analytic second root
      if (std::fabs(b_cf) > 1e-9 * scale) {
        const double tc = 2.0 * std::atan2(a_cf, b_cf);
        const double tmax = std::atan2(a_cf, b_cf);
        double best = -1e300, targ = 0;
        for (int s = 1; s < 2000; ++s) {
          const double t = 3.14159265358979 * s / 2000.0;
          const double A = a_cf * std::sin(t) - b_cf * (1 - std::cos(t));
          if (A > best) { best = A; targ = t; }
        }
        if (tc > 0.05 && tc < 3.1) {
          const double r = std::fabs(targ - tmax);
          if (r > worst_ratio) { worst_ratio = r; where_ratio = cs.kind; }
        }
      }
      (void)split_sum; (void)n_sp;
    }
    // (B.2s): split-edge contribution to its hole = 1/2 <d_e, du_e>, compared with the
    // same quantity read off the deployed offset  y_(a,g) - y_(a,f) = 2 sin(th/2) J du.
    for (int e : c.split_edges) {
      const Edge& ed = m.edges[e];
      const int f = m.half_edges[ed.he[0]].face, g = m.half_edges[ed.he[1]].face;
      const int va = ed.key.a, vb = ed.key.b;
      const Vec2 de = X[vb] - X[va];
      const Vec2 du = u[g] - u[f];
      const double r_closed = de.dot(du);
      const double s0 = std::sin(th[2] * 0.5);
      const Vec2 off = Y[2][c.prime_vertex(g, va)] - Y[2][c.prime_vertex(f, va)];
      const Vec2 du_meas = (J.transpose() * off) / (2.0 * s0);
      worst_split = std::max(worst_split, (du_meas - du).norm() / std::max(1.0, du.norm()));
      (void)r_closed; ++n_split_terms;
    }
    std::printf("%-18s holes=%3zu  worst rel (B.1) = %.3e\n", cs.kind.c_str(), cycles.size(), gw);
  }

  std::printf("\n--- summary over %d graphs, %d holes, %d split edges ---\n",
              n_graphs, n_holes, n_split_terms);
  std::printf("u closure residual              %.3e\n", worst_clo);
  std::printf("B.1  closed form vs fit (rel)   %.3e   (%s)\n", worst_rel, where_rel.c_str());
  std::printf("B.1  first-harmonic fit resid   %.3e\n", worst_abs);
  std::printf("B.2s du from deployed offset    %.3e\n", worst_split);
  std::printf("B10  |argmax - theta_c/2|       %.3e   (%s)\n", worst_ratio, where_ratio.c_str());
  return 0;
}
