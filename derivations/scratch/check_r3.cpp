// check_r3.cpp -- Deriver round 3.  Settles the ONE numeric question left by
// derivations/check.md R2.5: with the THIRD structural class g(0) = g'(0) = 0
// deflated as well, how many harmonics does the atom list get wrong?
//
// Build (from repo root):
//   clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
//     -I code/src derivations/scratch/check_r3.cpp \
//     code/src/core/mesh.cpp code/src/core/cut.cpp code/src/core/holes.cpp \
//     code/src/core/kinematics.cpp code/src/core/orientation.cpp code/src/core/generators.cpp \
//     code/src/core/tutte_auxetic.cpp code/src/core/collision.cpp code/src/core/optimize.cpp \
//     -o /tmp/check_r3 && /tmp/check_r3 0.02
//
// The round-2 program compared the deflated atom list against a reference that used the
// SAME deflation rule, so it could not see the third class; check.md R2.5 is right about
// that.  Here the truth is computed WITHOUT the tau chart at all: h(theta) = p + q cos +
// r sin is a degree-1 trig polynomial, so it is monotone between its critical points
// (tan theta = r/q).  Split (0, eps) at those, and a CROSSING exists iff two consecutive
// break points carry strictly opposite signs.  That is an independent decision procedure
// for "is there a contact event in (0, eps)", and it is well posed on the third class,
// where h = p(1 - cos theta) has constant sign and therefore never crosses.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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

static double cross2(const Vec2& a, const Vec2& b) { return a.x() * b.y() - a.y() * b.x(); }

struct Harm {
  double p = 0, q = 0, r = 0;
  double eval(double th) const { return p + q * std::cos(th) + r * std::sin(th); }
  double scale() const { return std::fabs(p) + std::hypot(q, r); }
};

static Harm orient_h(const Vec2& Cab, const Vec2& Sab, const Vec2& Caw, const Vec2& Saw) {
  Harm h;
  const double dcc = cross2(Cab, Caw), dss = cross2(Sab, Saw);
  h.p = 0.5 * (dcc + dss);
  h.q = 0.5 * (dcc - dss);
  h.r = 0.5 * (cross2(Cab, Saw) + cross2(Sab, Caw));
  return h;
}

// ---------------------------------------------------------------------------
// The atom list of T5.2b.1 / T5.2b.2, now with THREE structural classes.
//   A = p - q,  B = 2r,  C = p + q,   g(tau) = C + B tau + A tau^2,  T = tan(eps/2).
enum Klass { K_GENERIC = 0, K_SIMPLE0 = 1, K_DOUBLE0 = 2 };

static Klass classify(const Harm& h, double tol) {
  const double C = h.p + h.q, B = 2 * h.r;
  if (std::fabs(C) > tol) return K_GENERIC;
  if (std::fabs(B) > tol) return K_SIMPLE0;
  return K_DOUBLE0;
}

// returns true iff the atom list says "no admissible root in (0,T)".
static bool atoms_no_root(const Harm& h, double T, Klass k) {
  const double A = h.p - h.q, B = 2 * h.r, C = h.p + h.q;
  if (k == K_DOUBLE0) return true;                      // (T5.1e): g = A tau^2, root only at 0
  if (k == K_SIMPLE0) {                                 // (T5.1d): g = tau (A tau + B)
    return !((A * B < 0) && (-A * B - A * A * T < 0));
  }
  const double g0 = C, gT = C + B * T + A * T * T;      // (T5.1c)
  if (!(g0 * gT > 0)) return false;
  const double disc = B * B - 4 * A * C;
  const bool both = (disc >= 0) && (A * g0 > 0) && (A * gT > 0) && (A * B < 0) &&
                    (-A * B - 2 * A * A * T < 0);
  return !both;
}

// ---------------------------------------------------------------------------
// Independent truth: does h CROSS zero somewhere in the open interval (0, eps)?
// Uses only evaluations of h and its critical points -- no tau chart, no quadratic.
// Returns  0 = no crossing, 1 = crossing, -1 = ambiguous (a break-point value sits in
// the noise band, so double precision cannot decide; those are reported separately).
static int truth_crossing(const Harm& h, double eps, double tol) {
  // break points: 0+, the critical points of h inside (0, eps), and eps-.
  std::vector<double> br;
  const double base = std::atan2(h.r, h.q);   // h' = -q sin + r cos = 0  <=>  tan th = r/q
  for (int k = -2; k <= 2; ++k) {
    const double th = base + k * M_PI;
    if (th > 0 && th < eps) br.push_back(th);
  }
  std::sort(br.begin(), br.end());

  std::vector<int> sgn;
  // sign of h just to the RIGHT of 0, from the Taylor coefficients:
  //   h(0) = C,  h'(0) = r = B/2,  h''(0) = -q = (A - C)/2 ... -> class-wise.
  {
    const double C = h.p + h.q, B = 2 * h.r, A = h.p - h.q;
    double s;
    if (std::fabs(C) > tol) s = C;
    else if (std::fabs(B) > tol) s = B;
    else s = A;                 // h = p(1 - cos th), p = A/2, sign of A on all of (0, pi)
    if (std::fabs(s) <= tol) return -1;
    sgn.push_back(s > 0 ? 1 : -1);
  }
  for (double th : br) {
    const double v = h.eval(th);
    if (std::fabs(v) <= tol) return -1;
    sgn.push_back(v > 0 ? 1 : -1);
  }
  {
    const double v = h.eval(eps);
    if (std::fabs(v) <= tol) return -1;     // root at (or within noise of) the endpoint
    sgn.push_back(v > 0 ? 1 : -1);
  }
  for (size_t i = 1; i < sgn.size(); ++i)
    if (sgn[i] != sgn[i - 1]) return 1;
  return 0;
}

// ---------------------------------------------------------------------------
struct Tally {
  long long tested = 0, ambiguous = 0;
  long long n_class[3] = {0, 0, 0};
  long long mm3 = 0;                 // mismatches, three-class atom list
  long long mm3_by_class[3] = {0, 0, 0};
  long long mm2 = 0;                 // mismatches, round-2 (two-class) atom list
  long long mm2_by_class[3] = {0, 0, 0};
};

static bool atoms_no_root_r2(const Harm& h, double T, Klass k) {
  // round-2 rule: only the C = 0 deflation, no third class.
  return atoms_no_root(h, T, k == K_DOUBLE0 ? K_SIMPLE0 : k);
}

static void scan_harmonics(const CutStructure& c, const std::vector<Vec2>& C,
                           const std::vector<Vec2>& S, double geom_scale, double eps,
                           Tally* tal) {
  const Mesh& m = *c.mesh;
  const int F = m.n_faces();
  const double tol = 1e-11 * geom_scale * geom_scale;
  const double T = std::tan(0.5 * eps);
  for (int f = 0; f < F; ++f) {
    const auto& pf = c.prime_faces[f];
    const int nf = static_cast<int>(pf.size());
    for (int k = 0; k < nf; ++k) {
      const int a = pf[k], b = pf[(k + 1) % nf];
      const Vec2 Cab = C[b] - C[a], Sab = S[b] - S[a];
      if (Cab.squaredNorm() <= 0) continue;
      for (int g = 0; g < F; ++g) {
        if (g == f) continue;
        for (int w : c.prime_faces[g]) {
          const Vec2 Caw = C[w] - C[a], Saw = S[w] - S[a];
          const Harm ho = orient_h(Cab, Sab, Caw, Saw);
          if (ho.scale() <= tol) continue;
          const int tr = truth_crossing(ho, eps, tol);
          const Klass kl = classify(ho, tol);
          ++tal->n_class[kl];
          if (tr < 0) { ++tal->ambiguous; continue; }
          ++tal->tested;
          const bool truth_noroot = (tr == 0);
          if (atoms_no_root(ho, T, kl) != truth_noroot) { ++tal->mm3; ++tal->mm3_by_class[kl]; }
          if (atoms_no_root_r2(ho, T, kl) != truth_noroot) { ++tal->mm2; ++tal->mm2_by_class[kl]; }
        }
      }
    }
  }
}

static bool positively_oriented(const Mesh& m, const std::vector<Vec2>& X) {
  for (int f = 0; f < m.n_faces(); ++f) {
    double A2 = 0;
    const auto& F = m.faces[f];
    for (size_t i = 0; i < F.size(); ++i) {
      const Vec2& a = X[F[i]];
      const Vec2& b = X[F[(i + 1) % F.size()]];
      A2 += a.x() * b.y() - a.y() * b.x();
    }
    if (A2 <= 0) return false;
  }
  return true;
}

int main(int argc, char** argv) {
  const double EPS = (argc > 1) ? std::atof(argv[1]) : 0.02;
  std::mt19937 rng(20260904);   // same seed / same corpus as check_r2.cpp

  struct Case { std::string kind; std::vector<double> par; };
  std::vector<Case> cases = {
      {"squares", {2.5}},   {"triangles", {2.2}},   {"hexagons", {2.5}},
      {"kagome", {2.2}},    {"snub_square", {2.2}}, {"truncated_square", {2.5}},
      {"t3_4_3_12", {3.0}}, {"delaunay", {40, 8}},  {"voronoi", {35, 8}},
      {"quad_random", {40, 8}},
  };

  Tally tal;
  int n_used = 0, n_samples = 0;

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
      if (m.n_faces() > 260) continue;
      std::vector<Vec2> X = matrix_to_points(sr.X0);
      if (!hole_residuals(c, X, hs).deployable(1e-7)) continue;
      if (!positively_oriented(m, X)) continue;
      ++n_used;

      double scale = 0;
      for (auto& p : X) scale = std::max(scale, p.norm());

      std::normal_distribution<double> gauss(0.0, 1.0);
      const int n_samp = (sr.dim_null > 0) ? 24 : 1;
      for (int s = 0; s < n_samp; ++s) {
        std::vector<Vec2> Xs = X;
        if (s > 0) {
          Eigen::MatrixXd t(sr.dim_null, 2);
          for (int i = 0; i < sr.dim_null; ++i)
            for (int j = 0; j < 2; ++j) t(i, j) = 0.02 * scale * gauss(rng);
          Eigen::MatrixXd Xm = sr.X0 + sr.Phi * t;
          Xs = matrix_to_points(Xm);
        }
        if (!positively_oriented(m, Xs)) continue;
        ++n_samples;
        Deployment ds = deploy(c, Xs, 0.0);
        std::vector<Vec2> Cs = ds.Y, Ss(c.n_prime_vertices);
        for (int i = 0; i < c.n_prime_vertices; ++i) Ss[i] = 2.0 * ds.dY_dtheta[i];
        double sc = 0;
        for (auto& p : Xs) sc = std::max(sc, p.norm());
        scan_harmonics(c, Cs, Ss, sc, EPS, &tal);
      }
    }
  }

  std::printf("\nR3  eps = %.4f   graphs used: %d   shape-space samples: %d\n", EPS, n_used,
              n_samples);
  std::printf("    harmonics decided: %lld    ambiguous (skipped): %lld\n", tal.tested,
              tal.ambiguous);
  std::printf("    class sizes:  g(0) != 0 : %lld    g(0)=0, g'(0) != 0 : %lld    "
              "g(0)=g'(0)=0 : %lld\n",
              tal.n_class[0], tal.n_class[1], tal.n_class[2]);
  std::printf("    round-2 atom list (two classes):   %lld mismatches  [%lld / %lld / %lld "
              "by class]\n", tal.mm2, tal.mm2_by_class[0], tal.mm2_by_class[1],
              tal.mm2_by_class[2]);
  std::printf("    round-3 atom list (three classes): %lld mismatches  [%lld / %lld / %lld "
              "by class]\n", tal.mm3, tal.mm3_by_class[0], tal.mm3_by_class[1],
              tal.mm3_by_class[2]);
  return 0;
}
