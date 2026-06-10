// check_t1_t2.cpp -- Deriver's own numeric checks for derivations/core.md, T1 and T2.
//
// Build (from repo root):
//   clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
//     -I code/src derivations/scratch/check_t1_t2.cpp \
//     code/src/core/mesh.cpp code/src/core/cut.cpp code/src/core/holes.cpp \
//     code/src/core/kinematics.cpp code/src/core/orientation.cpp code/src/core/generators.cpp \
//     code/src/core/tutte_auxetic.cpp code/src/core/collision.cpp code/src/core/optimize.cpp \
//     -o /tmp/check_t1_t2 && /tmp/check_t1_t2
//
// Checks, all against code/src/core (which the Deriver does not modify):
//   C1  the potential u on Gamma closes on EVERY hinge edge, not just a spanning tree
//       (path-independence of the path sum; T1 step 6)
//   C2  Y_theta = cos(th/2) C + sin(th/2) S with C_(v,f) = x_v and
//       S_(v,f) = J (2 u_f - sigma_f x_v), against deploy() -- fixes the sign convention
//   C3  hinge opening angle == theta at every hinge edge
//   C4  the two duplicates of a split edge are the SAME vector (parallel, T1 cor. ii)
//   C5  face signed areas are constant in theta (T3)
//   C6  A(Y_theta) sigma = 0 for every theta: the explicit face velocity field
//       (omega_f, w_f) satisfies the pin equation at every hinge, at every theta (T2)
//   C7  the vertex trajectory lies on the centred conic through C, S (T1 cor. i)

#include <cmath>
#include <cstdio>
#include <random>
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

static Eigen::Matrix2d Jmat() {
  Eigen::Matrix2d J;
  J << 0, -1, 1, 0;
  return J;
}

struct Worst {
  double v = 0;
  std::string where;
  void hit(double x, const std::string& w) {
    if (x > v) { v = x; where = w; }
  }
};

// BFS on Gamma computing the potential u_f with u_g - u_f = sigma_g * x_src(e).
// Returns u and the worst closure residual over NON-tree hinge edges (C1).
static std::vector<Vec2> potential_u(const CutStructure& c, const std::vector<Vec2>& X,
                                     double* worst_closure) {
  const Mesh& m = *c.mesh;
  const int F = m.n_faces();
  std::vector<std::vector<std::pair<int, int>>> adj(F);  // face -> (edge, other face)
  for (int e : c.hinge_edges) {
    const Edge& ed = m.edges[e];
    const int f1 = m.half_edges[ed.he[0]].face;
    const int f2 = m.half_edges[ed.he[1]].face;
    adj[f1].push_back({e, f2});
    adj[f2].push_back({e, f1});
  }
  std::vector<Vec2> u(F, Vec2::Zero());
  std::vector<char> seen(F, 0);
  for (int s = 0; s < F; ++s) {
    if (seen[s]) continue;
    seen[s] = 1;
    u[s] = Vec2::Zero();
    std::vector<int> q{s};
    for (size_t qi = 0; qi < q.size(); ++qi) {
      const int f = q[qi];
      for (auto [e, g] : adj[f]) {
        const Vec2& xv = X[c.hinge_dir[e].src];
        const Vec2 target = u[f] + double(m.sigma[g]) * xv;
        if (!seen[g]) {
          seen[g] = 1;
          u[g] = target;
          q.push_back(g);
        }
      }
    }
  }
  // closure over every hinge edge
  double w = 0;
  for (int e : c.hinge_edges) {
    const Edge& ed = m.edges[e];
    const int f1 = m.half_edges[ed.he[0]].face;
    const int f2 = m.half_edges[ed.he[1]].face;
    const Vec2& xv = X[c.hinge_dir[e].src];
    w = std::max(w, (u[f2] - u[f1] - double(m.sigma[f2]) * xv).norm());
    w = std::max(w, (u[f1] - u[f2] - double(m.sigma[f1]) * xv).norm());
  }
  *worst_closure = w;
  return u;
}

static double signed_area(const std::vector<Vec2>& P) {
  double a = 0;
  for (size_t i = 0; i < P.size(); ++i) {
    const Vec2& p = P[i];
    const Vec2& q = P[(i + 1) % P.size()];
    a += p.x() * q.y() - p.y() * q.x();
  }
  return 0.5 * a;
}

int main() {
  const Eigen::Matrix2d J = Jmat();
  std::mt19937 rng(20260903);

  struct Case { std::string kind; std::vector<double> par; };
  std::vector<Case> cases = {
      {"squares", {3.5}},   {"triangles", {3.5}},   {"hexagons", {3.5}},
      {"kagome", {3.5}},    {"snub_square", {3.5}}, {"truncated_square", {3.5}},
      {"t3_4_3_12", {4.0}}, {"delaunay", {70, 10}}, {"voronoi", {60, 10}},
      {"quad_random", {60, 10}},
  };
  const std::vector<double> thetas = {0.0, 0.05, 0.3, 0.77, 1.2, 1.9, 2.6, 3.0};

  Worst c1, c2, c3, c4, c5, c6, c7;
  int n_graphs = 0, n_deployable = 0;

  for (const auto& cs : cases) {
    for (int rep = 0; rep < 3; ++rep) {
      Mesh m = generate(cs.kind, cs.par, rng);
      m.build_topology();
      auto orep = assign_orientation_relaxation(m, rng);
      m.sigma = orep.sigma;
      m.build_topology();
      CutStructure c = make_cut(m);
      HoleSet hs = holes_partition(c);
      LinearSystem sys = assemble_system(c, hs, m.X, BoundaryMode::Fixed);
      SolveReport sr = solve_system(sys, m.X);
      if (!sr.projection_ok || sr.dim_null < 0) continue;
      std::vector<Vec2> X = matrix_to_points(sr.X0);
      Residuals res = hole_residuals(c, X, hs);
      ++n_graphs;
      if (!res.deployable(1e-7)) continue;  // only test on genuinely deployable X
      ++n_deployable;

      double closure = 0;
      std::vector<Vec2> u = potential_u(c, X, &closure);
      c1.hit(closure, cs.kind);

      // C2/C3/C4/C5/C6/C7
      for (double th : thetas) {
        const double cc = std::cos(th * 0.5), ss = std::sin(th * 0.5);
        Deployment d = deploy(c, X, th);
        // closed form
        std::vector<Vec2> Y(c.n_prime_vertices, Vec2::Zero());
        for (int f = 0; f < m.n_faces(); ++f) {
          for (int v : m.faces[f]) {
            const int pv = c.prime_vertex(f, v);
            const Vec2 C = X[v];
            const Vec2 S = J * (2.0 * u[f] - double(m.sigma[f]) * X[v]);
            Y[pv] = cc * C + ss * S;
          }
        }
        // deploy() anchors face 0 at t=0 too, so the two agree exactly (no alignment).
        for (int i = 0; i < c.n_prime_vertices; ++i)
          c2.hit((Y[i] - d.Y[i]).norm(), cs.kind + " th=" + std::to_string(th));

        // C3: hinge opening angle
        for (int e : c.hinge_edges) {
          const Edge& ed = m.edges[e];
          const int f1 = m.half_edges[ed.he[0]].face;
          const int f2 = m.half_edges[ed.he[1]].face;
          const int a = c.hinge_dir[e].src, b = c.hinge_dir[e].dst;
          const Vec2 d1 = Y[c.prime_vertex(f1, b)] - Y[c.prime_vertex(f1, a)];
          const Vec2 d2 = Y[c.prime_vertex(f2, b)] - Y[c.prime_vertex(f2, a)];
          const double ang = std::atan2(d1.x() * d2.y() - d1.y() * d2.x(), d1.dot(d2));
          c3.hit(std::fabs(std::fabs(ang) - th), cs.kind);
        }
        // C4: split duplicates identical vectors
        for (int e : c.split_edges) {
          const Edge& ed = m.edges[e];
          const int f1 = m.half_edges[ed.he[0]].face;
          const int f2 = m.half_edges[ed.he[1]].face;
          const int a = ed.key.a, b = ed.key.b;
          const Vec2 d1 = Y[c.prime_vertex(f1, b)] - Y[c.prime_vertex(f1, a)];
          const Vec2 d2 = Y[c.prime_vertex(f2, b)] - Y[c.prime_vertex(f2, a)];
          c4.hit((d1 - d2).norm(), cs.kind);
        }
        // C5: face signed areas constant
        for (int f = 0; f < m.n_faces(); ++f) {
          std::vector<Vec2> Pf, Pflat;
          for (int v : m.faces[f]) {
            Pf.push_back(Y[c.prime_vertex(f, v)]);
            Pflat.push_back(X[v]);
          }
          c5.hit(std::fabs(signed_area(Pf) - signed_area(Pflat)), cs.kind);
        }
        // C6: pin equation for the velocity field (omega_f, w_f) at this theta
        //     omega_f = -sigma_f/2,  t_f = 2 s J u_f,  w_f = dt_f/dth - omega_f J t_f
        //                                        = c J u_f - sigma_f s u_f
        for (int e : c.hinge_edges) {
          const Edge& ed = m.edges[e];
          const int f = m.half_edges[ed.he[0]].face;
          const int g = m.half_edges[ed.he[1]].face;
          const Vec2 p = Y[c.prime_vertex(f, c.hinge_dir[e].src)];
          const Vec2 wf = cc * (J * u[f]) - double(m.sigma[f]) * ss * u[f];
          const Vec2 wg = cc * (J * u[g]) - double(m.sigma[g]) * ss * u[g];
          const double om_f = -0.5 * m.sigma[f], om_g = -0.5 * m.sigma[g];
          c6.hit(((wf - wg) + (om_f - om_g) * (J * p)).norm(), cs.kind);
        }
        // C7: trajectory on the centred conic y^T (A A^T)^-1 y = 1, A = [C S]
        if (th > 0) {
          for (int f = 0; f < m.n_faces(); ++f) {
            for (int v : m.faces[f]) {
              const int pv = c.prime_vertex(f, v);
              Eigen::Matrix2d A;
              A.col(0) = X[v];
              A.col(1) = J * (2.0 * u[f] - double(m.sigma[f]) * X[v]);
              const double det = A.determinant();
              if (std::fabs(det) < 1e-8) continue;  // degenerate ellipse, skip
              const Eigen::Matrix2d Q = (A * A.transpose()).inverse();
              c7.hit(std::fabs(Y[pv].dot(Q * Y[pv]) - 1.0), cs.kind);
            }
          }
        }
      }
    }
  }

  std::printf("graphs solved: %d, deployable X0 used: %d\n", n_graphs, n_deployable);
  auto rep = [](const char* n, const Worst& w) {
    std::printf("  %-58s max = %.3e   (%s)\n", n, w.v, w.where.c_str());
  };
  rep("C1 potential closes on EVERY hinge edge (path-independence)", c1);
  rep("C2 Y = cos(th/2) C + sin(th/2) S vs deploy()", c2);
  rep("C3 |hinge opening angle| - theta", c3);
  rep("C4 split-edge duplicates are the same vector", c4);
  rep("C5 face signed area constant in theta", c5);
  rep("C6 pin equation for (omega,w) = no-locking, A(Y_th) sigma = 0", c6);
  rep("C7 vertex trajectory on the centred conic through C, S", c7);
  return 0;
}
