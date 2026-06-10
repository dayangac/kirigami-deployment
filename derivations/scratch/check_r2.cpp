// check_r2.cpp -- Deriver round 2.  Decides the five numeric questions raised by
// derivations/check.md (D1, D2, D3, D4, D5) with a program written for this purpose.
//
// Build (from repo root):
//   clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
//     -I code/src derivations/scratch/check_r2.cpp \
//     code/src/core/mesh.cpp code/src/core/cut.cpp code/src/core/holes.cpp \
//     code/src/core/kinematics.cpp code/src/core/orientation.cpp code/src/core/generators.cpp \
//     code/src/core/tutte_auxetic.cpp code/src/core/collision.cpp code/src/core/optimize.cpp \
//     -o /tmp/check_r2 && /tmp/check_r2
//
// R2-A (D1)  In the face's OWN moving-centroid frame -- the frame ideas/ranking.md K2c
//            specifies -- is max(|x|,|chi|) an exact swept radius, or a sqrt(2) under-
//            estimate?  Reports max_theta |y_u(theta) - gamma_f(theta)| against both
//            max(|x|,|chi|) and the flat circumradius |x_u - xbar_f|.
// R2-B (D2)  Does the swept region stay within O(r_f) of the FLAT centroid (the form the
//            O(n) packing count needs)?  Growing square patch, reports
//            max_f max_theta |gamma_f(theta) - xbar_f| / r_f against patch diameter.
// R2-C (D3)  Split-free patterns: Theta_max == min(min_e beta_e, pi)?
// R2-D (D4)  Proposition T5.2b' (the inner-approximation direction).  For samples in the
//            shape space: (all A_f > 0) AND (not penetrating at 0+) AND (no orientation
//            harmonic of the COMPLETE candidate list has a root in (0,eps))  ==>
//            Theta_max >= eps.  Counts hypothesis hits and violations.
// R2-E (D5)  The degree-<= 4 atom list for "the tau-quadratic has no root in (0,T)",
//            checked against direct root computation on every harmonic encountered.

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
static Harm dot_h(const Vec2& Cab, const Vec2& Sab, const Vec2& Caw, const Vec2& Saw) {
  Harm h;
  const double dcc = Cab.dot(Caw), dss = Sab.dot(Saw);
  h.p = 0.5 * (dcc + dss);
  h.q = 0.5 * (dcc - dss);
  h.r = 0.5 * (Cab.dot(Saw) + Sab.dot(Caw));
  return h;
}

// Roots of p + q cos th + r sin th in (lo, hi], ascending.  (T3.5)/(T3.6).
static std::vector<double> roots_in(const Harm& h, double lo, double hi, double tol) {
  std::vector<double> out;
  if (h.scale() <= tol) return out;
  const double A = h.p - h.q, B = 2 * h.r, Cc = h.p + h.q;
  std::vector<double> taus;
  if (std::fabs(A) <= 1e-14 * (std::fabs(B) + std::fabs(Cc) + 1e-300)) {
    if (std::fabs(B) > 0) taus.push_back(-Cc / B);
    out.push_back(M_PI);
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

// ---------------------------------------------------------------------------
// R2-E: the degree-<= 4 atom list of check.md D5.
// g(tau) = C + B tau + A tau^2 with A = p-q, B = 2r, C = p+q.
// "both roots in (0,T)" :=  disc >= 0 AND A g(0) > 0 AND A g(T) > 0 AND A B < 0
//                           AND -A B - 2 A^2 T < 0      (i.e. vertex -B/(2A) < T)
// "no root in (0,T)"    :=  g(0) g(T) > 0  AND NOT["both roots in (0,T)"]
static bool no_root_in_atoms(const Harm& h, double T) {
  const double A = h.p - h.q, B = 2 * h.r, C = h.p + h.q;
  const double g0 = C, gT = C + B * T + A * T * T;
  if (!(g0 * gT > 0)) return false;
  const double disc = B * B - 4 * A * C;
  const bool both = (disc >= 0) && (A * g0 > 0) && (A * gT > 0) && (A * B < 0) &&
                    (-A * B - 2 * A * A * T < 0);
  return !both;
}
// Deflated form: when g(0) = 0 identically (a permanent incidence at theta = 0,
// T3.H.2 / T5.3), divide out the tau factor and test the linear remainder.
static bool no_root_in_atoms_deflated(const Harm& h, double T, double czero_tol) {
  const double A = h.p - h.q, B = 2 * h.r, C = h.p + h.q;
  if (std::fabs(C) <= czero_tol) {
    // g = tau (A tau + B); root -B/A in (0,T)?
    const bool inside = (A * B < 0) && (-A * B - A * A * T < 0);
    return !inside;
  }
  return no_root_in_atoms(h, T);
}
// Direct reference: does g have a root in the OPEN interval (0,T)?
static bool no_root_in_direct(const Harm& h, double T) {
  const double A = h.p - h.q, B = 2 * h.r, C = h.p + h.q;
  if (std::fabs(A) <= 1e-300) {
    if (std::fabs(B) <= 1e-300) return true;  // constant, no finite root (C != 0 assumed)
    const double t = -C / B;
    return !(t > 0 && t < T);
  }
  const double disc = B * B - 4 * A * C;
  if (disc < 0) return true;
  const double sq = std::sqrt(disc);
  const double t1 = (-B + sq) / (2 * A), t2 = (-B - sq) / (2 * A);
  const bool in1 = (t1 > 0 && t1 < T), in2 = (t2 > 0 && t2 < T);
  return !(in1 || in2);
}

// ---------------------------------------------------------------------------
struct RangeInfo {
  double theta_max = M_PI;
  std::vector<double> candidates;
  bool any_orient_root_below = false;  // any orientation-harmonic root in (0, eps)
};

// Complete candidate enumeration over ordered (M'-vertex w of g, edge (a,b) of f != g).
// `eps` only drives the extra flag; theta_max is by the T4.2" midpoint scan.
static RangeInfo exact_range(const CutStructure& c, const std::vector<Vec2>& X,
                             const std::vector<Vec2>& C, const std::vector<Vec2>& S,
                             double geom_scale, double eps, long long* n_atom_tests,
                             long long* n_atom_mismatch, long long* n_mm_czero,
                             long long* n_atom_tests_d, long long* n_mm_d) {
  const Mesh& m = *c.mesh;
  const int F = m.n_faces();
  const double tol = 1e-11 * geom_scale * geom_scale;
  const double T = std::tan(0.5 * eps);
  RangeInfo out;
  for (int f = 0; f < F; ++f) {
    const auto& pf = c.prime_faces[f];
    const int nf = static_cast<int>(pf.size());
    for (int k = 0; k < nf; ++k) {
      const int a = pf[k], b = pf[(k + 1) % nf];
      const Vec2 Cab = C[b] - C[a], Sab = S[b] - S[a];
      const double len2 = Cab.squaredNorm();
      if (len2 <= 0) continue;
      for (int g = 0; g < F; ++g) {
        if (g == f) continue;
        for (int w : c.prime_faces[g]) {
          const Vec2 Caw = C[w] - C[a], Saw = S[w] - S[a];
          const Harm ho = orient_h(Cab, Sab, Caw, Saw);
          const Harm hd = dot_h(Cab, Sab, Caw, Saw);
          if (ho.scale() > tol) {
            // R2-E: atom list vs direct roots, on this very harmonic.
            const double czt = 1e-11 * geom_scale * geom_scale;
            const bool a1 = no_root_in_atoms(ho, T), a2 = no_root_in_direct(ho, T);
            ++(*n_atom_tests);
            if (a1 != a2) {
              ++(*n_atom_mismatch);
              if (std::fabs(ho.p + ho.q) <= czt) ++(*n_mm_czero);
            }
            const bool a3 = no_root_in_atoms_deflated(ho, T, czt);
            const bool a4 = (std::fabs(ho.p + ho.q) <= czt)
                                ? no_root_in_direct(Harm{-0.5 * (ho.p - ho.q) * 0 + 0.0,
                                                         0.0, 0.0}, T)
                                : a2;
            (void)a4;
            ++(*n_atom_tests_d);
            // reference for the deflated test: roots of the deflated polynomial
            bool ref;
            {
              const double A = ho.p - ho.q, B = 2 * ho.r, C = ho.p + ho.q;
              if (std::fabs(C) <= czt) {
                if (std::fabs(A) <= 1e-300) ref = true;
                else { const double t = -B / A; ref = !(t > 0 && t < T); }
              } else ref = a2;
            }
            if (a3 != ref) ++(*n_mm_d);
            // the D4 hypothesis is about ORIENTATION roots only (interval test dropped),
            // with the theta = 0 incidence deflated away.
            if (!ref) out.any_orient_root_below = true;
          }
          for (double th : roots_in(ho, 1e-9, M_PI, tol)) {
            const double d = hd.eval(th);
            const double e2 = 1e-9 * len2;
            if (d < -e2 || d > len2 + e2) continue;
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

  const double SHRINK = 1e-7;
  {
    const double hi0 = out.candidates.empty() ? M_PI : out.candidates.front();
    Deployment d = deploy(c, X, 0.5 * hi0);
    if (has_collision(c, d.Y, SHRINK)) { out.theta_max = 0.0; return out; }
  }
  out.theta_max = M_PI;
  for (size_t i = 0; i < out.candidates.size(); ++i) {
    const double hi = (i + 1 < out.candidates.size()) ? out.candidates[i + 1] : M_PI;
    if (hi <= out.candidates[i] + 1e-12) continue;
    const double mid = 0.5 * (out.candidates[i] + hi);
    Deployment d = deploy(c, X, mid);
    if (has_collision(c, d.Y, SHRINK)) { out.theta_max = out.candidates[i]; break; }
  }
  return out;
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

// "no interior overlap at theta = 0+": overlap probe just above zero, below the first
// candidate.  This is exactly the negation of contact.hpp::penetrates_immediately().
static bool embedded_at_zero_plus(const CutStructure& c, const std::vector<Vec2>& X,
                                  double first_candidate) {
  const double probe = std::min(1e-6, 0.5 * first_candidate);
  Deployment d = deploy(c, X, probe);
  return !has_collision(c, d.Y, 1e-9);
}

int main(int argc, char** argv) {
  std::mt19937 rng(20260904);
  const double EPS = (argc > 1) ? std::atof(argv[1]) : 0.02;

  long long atom_tests = 0, atom_mismatch = 0, mm_czero = 0, atom_tests_d = 0, mm_d = 0;

  // -----------------------------------------------------------------  R2-A / R2-C
  struct Case { std::string kind; std::vector<double> par; };
  std::vector<Case> cases = {
      {"squares", {2.5}},   {"triangles", {2.2}},   {"hexagons", {2.5}},
      {"kagome", {2.2}},    {"snub_square", {2.2}}, {"truncated_square", {2.5}},
      {"t3_4_3_12", {3.0}}, {"delaunay", {40, 8}},  {"voronoi", {35, 8}},
      {"quad_random", {40, 8}},
  };

  double a_worst_specform = 0;   // | max_theta |y-gamma| - max(|x|,|chi|) |
  double a_worst_circum = 0;     // | max_theta |y-gamma| - |x_u - xbar_f| |
  long long a_copies = 0;
  double c_worst = 0;
  int n_split_free = 0, n_used = 0;

  // -----------------------------------------------------------------  R2-D
  long long d_hyp = 0, d_viol = 0, d_samples = 0, d_short = 0;
  double d_worst_viol = 0;

  std::printf("%-18s %5s %6s  %-10s %-10s %-10s\n", "graph", "F", "cands", "Theta_max",
              "min_beta", "min(b,pi)");

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

      Deployment d0 = deploy(c, X, 0.0);
      std::vector<Vec2> C = d0.Y, S(c.n_prime_vertices);
      for (int i = 0; i < c.n_prime_vertices; ++i) S[i] = 2.0 * d0.dY_dtheta[i];

      // ---- R2-A: swept radius in the face's own moving-centroid frame.
      for (int f = 0; f < m.n_faces(); ++f) {
        const auto& pf = c.prime_faces[f];
        Vec2 gc = Vec2::Zero(), gs = Vec2::Zero(), xbar = Vec2::Zero();
        for (size_t i = 0; i < pf.size(); ++i) { gc += C[pf[i]]; gs += S[pf[i]]; }
        gc /= double(pf.size());
        gs /= double(pf.size());
        for (int v : m.faces[f]) xbar += X[v];
        xbar /= double(m.faces[f].size());
        for (int u : pf) {
          const Vec2 x = C[u] - gc, chi = S[u] - gs;
          const double specform = std::max(x.norm(), chi.norm());
          double worst = 0;
          for (int k = 0; k <= 128; ++k) {
            const double th = M_PI * k / 128.0;
            worst = std::max(worst, (std::cos(th / 2) * x + std::sin(th / 2) * chi).norm());
          }
          a_worst_specform = std::max(a_worst_specform, std::fabs(worst - specform));
          a_worst_circum = std::max(a_worst_circum, std::fabs(worst - x.norm()));
          ++a_copies;
        }
      }

      RangeInfo ri = exact_range(c, X, C, S, scale, EPS, &atom_tests, &atom_mismatch, &mm_czero, &atom_tests_d, &mm_d);
      ThetaMaxReport tm = theta_max(c, X);

      // ---- R2-C: split-free  Theta_max == min(min beta, pi)
      if (c.n_split() == 0) {
        ++n_split_free;
        c_worst = std::max(c_worst,
                           std::fabs(ri.theta_max - std::min(tm.min_beta, M_PI)));
        std::printf("%-18s %5d %6zu  %-10.6f %-10.6f %-10.6f\n", cs.kind.c_str(), m.n_faces(),
                    ri.candidates.size(), ri.theta_max, tm.min_beta,
                    std::min(tm.min_beta, M_PI));
      }

      // ---- R2-D: Proposition T5.2b' on X0 and on null-space samples.
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
        ++d_samples;
        Deployment ds = deploy(c, Xs, 0.0);
        std::vector<Vec2> Cs = ds.Y, Ss(c.n_prime_vertices);
        for (int i = 0; i < c.n_prime_vertices; ++i) Ss[i] = 2.0 * ds.dY_dtheta[i];
        double sc = 0;
        for (auto& p : Xs) sc = std::max(sc, p.norm());
        RangeInfo rs = exact_range(c, Xs, Cs, Ss, sc, EPS, &atom_tests, &atom_mismatch, &mm_czero, &atom_tests_d, &mm_d);
        const double firstc = rs.candidates.empty() ? M_PI : rs.candidates.front();
        const bool emb0 = embedded_at_zero_plus(c, Xs, firstc);
        if (rs.theta_max < EPS - 1e-9) ++d_short;  // samples the implication must exclude
        if (rs.any_orient_root_below || !emb0) continue;  // hypothesis not met
        ++d_hyp;
        if (rs.theta_max < EPS - 1e-9) {
          ++d_viol;
          d_worst_viol = std::max(d_worst_viol, EPS - rs.theta_max);
        }
      }
    }
  }

  // -----------------------------------------------------------------  R2-B
  std::printf("\nR2-B  growing square patch: drift of the moving centroid from the FLAT one\n");
  std::printf("%8s %5s %10s %14s %14s\n", "clip", "F", "diameter",
              "max|gam-xbar|/r", "max rho_spec/r");
  for (double R : {1.6, 2.6, 3.6, 4.6, 5.6, 7.0}) {
    Mesh m = tiling_squares(rect(Vec2(0, 0), R, R));
    m.build_topology();
    std::mt19937 r2(7);
    auto orep = assign_orientation_relaxation(m, r2);
    m.sigma = orep.sigma;
    m.build_topology();
    CutStructure c = make_cut(m);
    HoleSet hs = holes_partition(c);
    LinearSystem sys = assemble_system(c, hs, m.X, BoundaryMode::Fixed);
    SolveReport sr = solve_system(sys, m.X);
    if (!sr.projection_ok) continue;
    std::vector<Vec2> X = matrix_to_points(sr.X0);
    Deployment d0 = deploy(c, X, 0.0);
    std::vector<Vec2> C = d0.Y, S(c.n_prime_vertices);
    for (int i = 0; i < c.n_prime_vertices; ++i) S[i] = 2.0 * d0.dY_dtheta[i];
    double diam = 0;
    for (size_t i = 0; i < X.size(); ++i)
      for (size_t j = i + 1; j < X.size(); ++j) diam = std::max(diam, (X[i] - X[j]).norm());
    double worst_drift = 0, worst_spec = 0;
    for (int f = 0; f < m.n_faces(); ++f) {
      const auto& pf = c.prime_faces[f];
      Vec2 gc = Vec2::Zero(), gs = Vec2::Zero(), xbar = Vec2::Zero();
      for (int u : pf) { gc += C[u]; gs += S[u]; }
      gc /= double(pf.size());
      gs /= double(pf.size());
      for (int v : m.faces[f]) xbar += X[v];
      xbar /= double(m.faces[f].size());
      double rf = 0;
      for (int v : m.faces[f]) rf = std::max(rf, (X[v] - xbar).norm());
      if (rf <= 0) continue;
      double drift = 0;
      for (int k = 0; k <= 64; ++k) {
        const double th = M_PI * k / 64.0;
        const Vec2 gam = std::cos(th / 2) * gc + std::sin(th / 2) * gs;
        drift = std::max(drift, (gam - xbar).norm());
      }
      worst_drift = std::max(worst_drift, drift / rf);
      double spec = 0;
      for (int u : pf)
        spec = std::max(spec, std::max((C[u] - gc).norm(), (S[u] - gs).norm()));
      worst_spec = std::max(worst_spec, spec / rf);
    }
    std::printf("%8.1f %5d %10.3f %14.3f %14.6f\n", R, m.n_faces(), diam, worst_drift,
                worst_spec);
  }

  std::printf("\n== results ==\n");
  std::printf("R2-A  |max_theta|y-gamma| - max(|x|,|chi|)|  worst = %.3e   (copies %lld)\n",
              a_worst_specform, a_copies);
  std::printf("R2-A  |max_theta|y-gamma| - |x_u - xbar_f||  worst = %.3e\n", a_worst_circum);
  std::printf("R2-C  split-free |Theta_max - min(min beta, pi)|  worst = %.3e  (%d patterns)\n",
              c_worst, n_split_free);
  std::printf("R2-D  T5.2b' hypothesis met on %lld / %lld samples; VIOLATIONS = %lld "
              "(worst shortfall %.3e), eps = %.4f; %lld / %lld samples have "
              "Theta_max < eps\n", d_hyp, d_samples, d_viol, d_worst_viol, EPS, d_short,
              d_samples);
  std::printf("R2-E  D5 atom list as printed vs direct roots: %lld tests, %lld mismatches "
              "(%lld of them with g(0) = 0)\n", atom_tests, atom_mismatch, mm_czero);
  std::printf("R2-E  DEFLATED atom list vs deflated roots:  %lld tests, %lld mismatches\n",
              atom_tests_d, mm_d);
  std::printf("graphs used: %d\n", n_used);
  return 0;
}
