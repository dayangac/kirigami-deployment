// check_t4_t5.cpp -- Deriver's own numeric checks for derivations/core.md, T4 and T5.
//
// Deliberately re-implements the harmonic calculus from scratch (it does NOT use
// code/src/method/deploy_basis.*, which the Experimenter owns) so that agreement is
// an independent cross-check rather than a tautology.
//
// Build (from repo root):
//   clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
//     -I code/src derivations/scratch/check_t4_t5.cpp \
//     code/src/core/mesh.cpp code/src/core/cut.cpp code/src/core/holes.cpp \
//     code/src/core/kinematics.cpp code/src/core/orientation.cpp code/src/core/generators.cpp \
//     code/src/core/tutte_auxetic.cpp code/src/core/collision.cpp code/src/core/optimize.cpp \
//     -o /tmp/check_t4_t5 && /tmp/check_t4_t5
//
// D1  exact theta_c by closed-form enumeration over ALL ordered (vertex, edge) pairs
//     with the two interval harmonics, vs collision.hpp's grid+bisection theta_max
// D2  on split-free patterns, theta_c == min_i beta_i = min_i (2pi - alpha_i - alpha_j)
//     (2026 Sec. 5.2 / 2025 Sec. 4.2 falls out of the candidate list as a special case)
// D3  split-edge separation: closed-form p = -sigma_f det(d, du), q = -p,
//     r = det(d, J du) against the generic (T3.2) coefficients, and p + q == 0
// D4  is max(|C|,|S|) an upper bound for |Y(theta)| on [0, pi)?  (K2c's rho_f as
//     specified).  Reports the worst overshoot ratio.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "core/collision.hpp"
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
static double cross2(const Vec2& a, const Vec2& b) { return a.x() * b.y() - a.y() * b.x(); }

struct Harm {
  double p = 0, q = 0, r = 0;
  double eval(double th) const { return p + q * std::cos(th) + r * std::sin(th); }
  double scale() const { return std::fabs(p) + std::hypot(q, r); }
};

// (T3.2): det(Yb - Ya, Yw - Ya)
static Harm orient_h(const Vec2& Cab, const Vec2& Sab, const Vec2& Caw, const Vec2& Saw) {
  Harm h;
  const double dcc = cross2(Cab, Caw), dss = cross2(Sab, Saw);
  h.p = 0.5 * (dcc + dss);
  h.q = 0.5 * (dcc - dss);
  h.r = 0.5 * (cross2(Cab, Saw) + cross2(Sab, Caw));
  return h;
}
// (T3.3): <Yb - Ya, Yw - Ya>
static Harm dot_h(const Vec2& Cab, const Vec2& Sab, const Vec2& Caw, const Vec2& Saw) {
  Harm h;
  const double dcc = Cab.dot(Caw), dss = Sab.dot(Saw);
  h.p = 0.5 * (dcc + dss);
  h.q = 0.5 * (dcc - dss);
  h.r = 0.5 * (Cab.dot(Saw) + Sab.dot(Caw));
  return h;
}

// (T3.5)-(T3.6): roots of p + q cos th + r sin th in (lo, hi], ascending.
static std::vector<double> roots_in(const Harm& h, double lo, double hi, double tol) {
  std::vector<double> out;
  if (h.scale() <= tol) return out;  // identically zero to tolerance -> degenerate
  const double A = h.p - h.q, B = 2 * h.r, Cc = h.p + h.q;
  std::vector<double> taus;
  if (std::fabs(A) <= 1e-14 * (std::fabs(B) + std::fabs(Cc) + 1e-300)) {
    if (std::fabs(B) > 0) taus.push_back(-Cc / B);  // the other root is at theta = pi
    out.push_back(M_PI);                            // tau = infinity
  } else {
    const double disc = B * B - 4 * A * Cc;
    if (disc < 0) return out;
    const double sq = std::sqrt(disc);
    taus.push_back((-B + sq) / (2 * A));
    taus.push_back((-B - sq) / (2 * A));
  }
  for (double t : taus) out.push_back(2.0 * std::atan(t));
  std::vector<double> keep;
  for (double th : out)
    if (th > lo && th <= hi) keep.push_back(th);
  std::sort(keep.begin(), keep.end());
  return keep;
}

struct ExactRange {
  double theta_c = M_PI;              // first CONTACT angle (smallest admissible root)
  double theta_max = M_PI;            // first INTERIOR-OVERLAP angle
  std::vector<double> candidates;     // all admissible roots in (0, pi], ascending, deduped
};

// T4: collect every admissible contact root over ordered (M'-vertex w of face g,
// edge (a,b) of face f != g): a root of the orientation harmonic in (0, pi] whose
// two interval harmonics place w on the CLOSED edge.  Lemma T4.1 says the pairwise
// interior-overlap status of the structure can only change at such an angle, so a
// midpoint test between consecutive candidates gives the exact theta_max.
static ExactRange exact_range(const CutStructure& c, const std::vector<Vec2>& X,
                              const std::vector<Vec2>& C, const std::vector<Vec2>& S,
                              double geom_scale) {
  const Mesh& m = *c.mesh;
  const int F = m.n_faces();
  const double tol = 1e-11 * geom_scale * geom_scale;
  ExactRange out;
  for (int f = 0; f < F; ++f) {
    const auto& pf = c.prime_faces[f];
    const int nf = static_cast<int>(pf.size());
    for (int k = 0; k < nf; ++k) {
      const int a = pf[k], b = pf[(k + 1) % nf];
      const Vec2 Cab = C[b] - C[a], Sab = S[b] - S[a];
      const double len2 = Cab.squaredNorm();  // constant in theta (T3.3 corollary)
      if (len2 <= 0) continue;
      for (int g = 0; g < F; ++g) {
        if (g == f) continue;
        for (int w : c.prime_faces[g]) {
          const Vec2 Caw = C[w] - C[a], Saw = S[w] - S[a];
          const Harm ho = orient_h(Cab, Sab, Caw, Saw);
          const Harm hd = dot_h(Cab, Sab, Caw, Saw);
          for (double th : roots_in(ho, 1e-9, M_PI, tol)) {
            const double d = hd.eval(th);
            const double eps = 1e-9 * len2;
            if (d < -eps || d > len2 + eps) continue;   // CLOSED edge, endpoints included
            out.candidates.push_back(th);
          }
        }
      }
    }
  }
  std::sort(out.candidates.begin(), out.candidates.end());
  std::vector<double> ded;
  for (double t : out.candidates)
    if (ded.empty() || t - ded.back() > 1e-9) ded.push_back(t);
  out.candidates.swap(ded);
  out.theta_c = out.candidates.empty() ? M_PI : out.candidates.front();
  // midpoint scan
  out.theta_max = M_PI;
  const double SHRINK = 1e-7;
  // the first interval is (0, cand[0]): an X that is positively oriented but NOT
  // embedded already overlaps there, and theta_max = 0.
  {
    const double hi0 = out.candidates.empty() ? M_PI : out.candidates.front();
    Deployment d = deploy(c, X, 0.5 * hi0);
    if (has_collision(c, d.Y, SHRINK)) { out.theta_max = 0.0; return out; }
  }
  for (size_t i = 0; i < out.candidates.size(); ++i) {
    const double hi = (i + 1 < out.candidates.size()) ? out.candidates[i + 1] : M_PI;
    if (hi <= out.candidates[i] + 1e-12) continue;
    const double mid = 0.5 * (out.candidates[i] + hi);
    Deployment d = deploy(c, X, mid);
    if (has_collision(c, d.Y, SHRINK)) { out.theta_max = out.candidates[i]; break; }
  }
  return out;
}

int main() {
  const Eigen::Matrix2d J = Jmat();
  std::mt19937 rng(20260903);

  struct Case { std::string kind; std::vector<double> par; };
  std::vector<Case> cases = {
      {"squares", {2.5}},   {"triangles", {2.2}},   {"hexagons", {2.5}},
      {"kagome", {2.2}},    {"snub_square", {2.2}}, {"truncated_square", {2.5}},
      {"t3_4_3_12", {3.0}}, {"delaunay", {40, 8}},  {"voronoi", {35, 8}},
      {"quad_random", {40, 8}},
  };

  std::printf("%-18s %5s %5s %5s %6s  %-10s %-10s %-10s %-10s %-8s\n", "graph", "F", "hinge",
              "split", "cands", "theta_c", "theta_max", "bisect", "min_beta", "d1");
  double d1_worst = 0, d2_worst = 0, d3_worst = 0, d4_worst_ratio = 1.0;
  long d3_pairs = 0, d4_violations = 0, d4_total = 0;
  int n_used = 0, n_split_free = 0;

  for (const auto& cs : cases) {
    for (int rep = 0; rep < 2; ++rep) {
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
      if (!hole_residuals(c, X, hs).deployable(1e-7)) continue;
      if (m.n_faces() > 260) continue;
      // D1/D2 are statements about VALID embeddings; X0 is self-intersecting on random
      // graphs (F17), where "first contact" is at theta = 0 and both numbers are 0.
      bool positively_oriented = true;
      for (int f = 0; f < m.n_faces(); ++f) {
        std::vector<Vec2> P;
        for (int v : m.faces[f]) P.push_back(X[v]);
        double A2 = 0;
        for (size_t i = 0; i < P.size(); ++i) {
          const Vec2& a = P[i]; const Vec2& b = P[(i + 1) % P.size()];
          A2 += a.x() * b.y() - a.y() * b.x();
        }
        if (A2 <= 0) { positively_oriented = false; break; }
      }
      if (!positively_oriented) { std::printf("%-18s  SKIPPED (X0 has inverted faces, F17)\n", cs.kind.c_str()); continue; }
      ++n_used;

      // C and S from T1.5, via my own BFS potential (validated in check_t1_t2.cpp).
      Deployment d0 = deploy(c, X, 0.0);
      Deployment dh = deploy(c, X, 1e-5);
      // C = Y(0); S = 2 dY/dtheta|_0 -- taken from the code's analytic derivative so that
      // this program's basis is independent of check_t1_t2.cpp's BFS.
      std::vector<Vec2> C = d0.Y, S(c.n_prime_vertices);
      for (int i = 0; i < c.n_prime_vertices; ++i) S[i] = 2.0 * d0.dY_dtheta[i];
      (void)dh;

      double scale = 0;
      for (auto& p : X) scale = std::max(scale, p.norm());

      ExactRange er = exact_range(c, X, C, S, scale);
      ThetaMaxReport tm = theta_max(c, X);
      const double d1 = std::fabs(er.theta_max - tm.theta_max_geometric);
      d1_worst = std::max(d1_worst, d1);

      // D2: split-free patterns
      double beta_min = tm.min_beta;
      if (c.n_split() == 0) {
        ++n_split_free;
        d2_worst = std::max(d2_worst, std::fabs(std::min(er.theta_max, M_PI) - std::min(beta_min, M_PI)));
      }

      // D3: split-edge separation closed form
      // Recover u from S: S_(v,f) = J(2 u_f - sigma_f x_v) => u_f = (J^-1 S + sigma_f x_v)/2
      std::vector<Vec2> u(m.n_faces(), Vec2::Zero());
      for (int f = 0; f < m.n_faces(); ++f) {
        const int v = m.faces[f][0];
        const int pv = c.prime_vertex(f, v);
        u[f] = 0.5 * (J.transpose() * S[pv] + double(m.sigma[f]) * X[v]);
      }
      for (int e : c.split_edges) {
        const Edge& ed = m.edges[e];
        const int f = m.half_edges[ed.he[0]].face;
        const int g = m.half_edges[ed.he[1]].face;
        const int va = ed.key.a, vb = ed.key.b;
        const Vec2 dvec = X[vb] - X[va];
        const Vec2 du = u[g] - u[f];
        Harm pred;
        pred.p = -double(m.sigma[f]) * cross2(dvec, du);
        pred.q = -pred.p;
        pred.r = cross2(dvec, J * du);
        // generic route: orientation harmonic of (a',b' in f ; a'' in g)
        const int A = c.prime_vertex(f, va), B = c.prime_vertex(f, vb);
        const int W = c.prime_vertex(g, va);
        Harm gen = orient_h(C[B] - C[A], S[B] - S[A], C[W] - C[A], S[W] - S[A]);
        const double sc = std::max(1.0, gen.scale());
        d3_worst = std::max({d3_worst, std::fabs(pred.p - gen.p) / sc,
                             std::fabs(pred.q - gen.q) / sc, std::fabs(pred.r - gen.r) / sc,
                             std::fabs(gen.p + gen.q) / sc});
        ++d3_pairs;
      }

      // D4: is max(|C|,|S|) an upper bound for |Y(theta)| on [0, pi)?
      for (int i = 0; i < c.n_prime_vertices; ++i) {
        const double bound = std::max(C[i].norm(), S[i].norm());
        if (bound < 1e-12) continue;
        ++d4_total;
        double worst = 0;
        for (int k = 0; k <= 64; ++k) {
          const double th = M_PI * k / 65.0;
          worst = std::max(worst, (std::cos(th / 2) * C[i] + std::sin(th / 2) * S[i]).norm());
        }
        if (worst > bound * (1 + 1e-9)) {
          ++d4_violations;
          d4_worst_ratio = std::max(d4_worst_ratio, worst / bound);
        }
      }

      std::printf("%-18s %5d %5d %5d %6zu  %-10.6f %-10.6f %-10.6f %-10.6f %-8.1e\n",
                  cs.kind.c_str(), m.n_faces(), c.n_hinge(), c.n_split(), er.candidates.size(),
                  er.theta_c, er.theta_max, tm.theta_max_geometric, beta_min, d1);
    }
  }

  std::printf("\ngraphs used: %d (split-free: %d)\n", n_used, n_split_free);
  std::printf("D1  |theta_c(exact) - theta_max(bisection)|      worst = %.3e   (rule <= 1e-5)\n",
              d1_worst);
  std::printf("D2  |theta_max - min_i beta_i| on split-free     worst = %.3e\n", d2_worst);
  std::printf("D3  split-edge closed form vs generic (T3.2), and p+q == 0\n");
  std::printf("      pairs = %ld   worst relative deviation = %.3e\n", d3_pairs, d3_worst);
  std::printf("D4  max(|C|,|S|) as a swept-radius bound: %ld / %ld copies VIOLATE it, "
              "worst ratio %.4f\n", d4_violations, d4_total, d4_worst_ratio);
  return 0;
}
