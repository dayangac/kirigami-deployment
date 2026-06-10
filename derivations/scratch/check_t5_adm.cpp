// check_t5_adm.cpp -- Checker round 6.  Two questions about the AMENDED certificate
//   VALID(eps) = POS /\ NOOVERLAP(eps/2) /\ NOROOT_adm(eps)
// where NOROOT_adm asks for no ADMISSIBLE root (vertex on the edge SEGMENT), i.e.
// C(X) ^ (0,eps) = {} exactly (results/kill/jitter/cert_diagnosis.md Sec. 4).
//
// Build (from repo root; libkiri_core.a must exist -- cmake --build code/build):
//   clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
//     -Icode/src derivations/scratch/check_t5_adm.cpp code/build/libkiri_core.a \
//     -o /tmp/check_t5_adm && /tmp/check_t5_adm
//
// A  SOUNDNESS side, Sub-lemma T5.2b'' Case A.  The substitute pair pi' of core.md
//    (far endpoint w of the SHORTER far-side edge, against the other far-side edge
//    (a=p, b)) has h_o,pi'(beta_e) = 0.  The amendment needs more: that root must be
//    ADMISSIBLE.  Claim (proved in check.md Round 6, R6.1): at theta = beta_e the two
//    far-side edges are collinear and CO-DIRECTED, so
//        <w - a, b - a> = L' L   and   |b - a|^2 = L^2 ,   with  0 < L' <= L,
//    hence the projection ratio s/|e|^2 = L'/L lies in (0, 1].  Measured here.
//
// B  COMPLETENESS side.  Is "exact Theta_max >= eps  =>  certificate" true?  NO.
//    A GRAZE -- an admissible contact that is not an overlap transition -- sits in
//    (0, Theta_max) on the hexagon pattern (core.md T4.2 table: theta_1 = pi/3,
//    Theta_max = 2pi/3).  For any eps strictly between the two, Theta_max >= eps holds
//    while NOROOT_adm fails.  The 0/5064 exceptions of the jitter run are an artefact of
//    eps = 0.006 being far below every contact angle on that corpus.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
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
#include "method/contact.hpp"
#include "method/deploy_basis.hpp"

using namespace kiri;
using namespace kiri::method;

namespace {

struct Case {
  std::string name;
  Mesh m;
  CutStructure c;
  HoleSet hs;
  LinearSystem sys;
  SolveReport rep;
  std::vector<Vec2> X0;
  int k = 0;
  bool ok = false;

  std::vector<Vec2> sample(std::mt19937& rng, double scale) const {
    Eigen::MatrixXd X = rep.X0;
    if (k > 0 && scale != 0.0) {
      std::normal_distribution<double> g(0.0, scale);
      Eigen::MatrixXd T(k, 2);
      for (int i = 0; i < k; ++i) { T(i, 0) = g(rng); T(i, 1) = g(rng); }
      X += rep.Phi * T;
    }
    return matrix_to_points(X);
  }
};

std::vector<std::unique_ptr<Case>> corpus() {
  std::vector<std::unique_ptr<Case>> out;
  std::mt19937 rng(20260903u);
  auto add = [&](Mesh m, const std::string& nm) {
    if (m.n_faces() < 3) return;
    m.sigma = assign_orientation_relaxation(m, rng, 6, 400, 90).sigma;
    auto p = std::make_unique<Case>();
    p->name = nm;
    p->m = std::move(m);
    p->c = make_cut(p->m);
    p->c.mesh = &p->m;
    p->hs = holes_partition(p->c);
    p->sys = assemble_system(p->c, p->hs, p->m.X, BoundaryMode::Fixed);
    p->rep = solve_system(p->sys, p->m.X);
    p->X0 = matrix_to_points(p->rep.X0);
    p->k = p->rep.dim_null;
    const Residuals r = hole_residuals(p->c, p->X0, p->hs);
    p->ok = p->rep.projection_ok && r.max_norm < 1e-7 && p->c.n_hinge() > 0;
    if (p->ok) out.push_back(std::move(p));
  };
  add(tiling_squares(disk({0, 0}, 2.6)), "squares");
  add(tiling_triangles(disk({0, 0}, 2.2)), "triangles");
  add(tiling_hexagons(disk({0, 0}, 3.0)), "hexagons");
  add(tiling_kagome(disk({0, 0}, 2.4)), "kagome");
  add(tiling_snub_square(disk({0, 0}, 2.2)), "snub_square");
  add(tiling_truncated_square(disk({0, 0}, 2.8)), "truncated_square");
  add(tiling_3_4_3_12(disk({0, 0}, 3.2)), "t3_4_3_12");
  return out;
}

double interior_angle(const std::vector<int>& vs, const std::vector<Vec2>& X, int v) {
  const int n = static_cast<int>(vs.size());
  int i = -1;
  for (int j = 0; j < n; ++j)
    if (vs[j] == v) i = j;
  if (i < 0) return -1;
  const Vec2 a = X[vs[(i + n - 1) % n]] - X[v], b = X[vs[(i + 1) % n]] - X[v];
  if (a.norm() < 1e-14 || b.norm() < 1e-14) return -1;
  double t = std::atan2(b.x() * a.y() - b.y() * a.x(), a.dot(b));
  if (t < 0) t += 2 * M_PI;
  return t;  // interior angle at v (same convention as derivation_tests.cpp)
}

}  // namespace

int main() {
  auto C = corpus();
  std::printf("corpus: %zu cases\n\n", C.size());

  // ==================================================================== A
  // Case A: admissibility of the substitute pair's root at theta = beta_e.
  std::printf("A  Case A substitute pair pi' = (w, (a,b)) at theta = beta_e\n");
  std::printf("   admissible <=> 0 <= <w-a,b-a> <= |b-a|^2 at beta_e\n\n");
  long long nA = 0, nA_bad = 0, nA_degen = 0, n_beta_seen = 0, n_beta_out = 0, nA_strict = 0;
  double worst_ratio_lo = 1e300, worst_ratio_hi = -1e300, e_ratio = 0, e_root = 0;
  std::mt19937 rngA(60061u);
  for (auto& cs : C) {
   for (int rep = 0; rep < 40; ++rep) {
    const std::vector<Vec2> Xs = cs->sample(rngA, rep == 0 ? 0.0 : 0.12);
    const DeployBasis B = deploy_basis(cs->c, Xs);
    for (size_t ei = 0; ei < cs->c.hinge_edges.size(); ++ei) {
      const int e = cs->c.hinge_edges[ei];
      const Edge& ed = cs->m.edges[e];
      const int f = cs->m.half_edges[ed.he[0]].face;
      const int g = cs->m.half_edges[ed.he[1]].face;
      if (f < 0 || g < 0) continue;
      const int v = cs->c.hinge_dir[e].src, dst = cs->c.hinge_dir[e].dst;
      const auto& vf = cs->m.faces[f];
      const auto& vg = cs->m.faces[g];
      const double af = interior_angle(vf, Xs, v), ag = interior_angle(vg, Xs, v);
      if (af < 0 || ag < 0) continue;
      const double beta = 2 * M_PI - af - ag;
      ++n_beta_seen;
      if (!(beta > 1e-3 && beta < M_PI - 1e-3)) { ++n_beta_out; continue; }  // the case the proof uses

      // far-side neighbour of v in each face: the one that is NOT dst
      auto far_of = [&](const std::vector<int>& vs) {
        const int n = static_cast<int>(vs.size());
        int i = -1;
        for (int j = 0; j < n; ++j)
          if (vs[j] == v) i = j;
        const int pv = vs[(i + n - 1) % n], nx = vs[(i + 1) % n];
        const int far = (nx == dst) ? pv : nx;
        int fi = -1;
        for (int j = 0; j < n; ++j)
          if (vs[j] == far) fi = j;
        return std::pair<int, int>{i, fi};
      };
      const auto [if_v, if_far] = far_of(vf);
      const auto [ig_v, ig_far] = far_of(vg);
      if (if_far < 0 || ig_far < 0 || if_v < 0 || ig_v < 0) continue;
      const double lf = (Xs[vf[if_far]] - Xs[v]).norm();
      const double lg = (Xs[vg[ig_far]] - Xs[v]).norm();
      int w, a, b;
      double Lp, L;
      if (lf <= lg) {  // shorter far edge belongs to f: its far endpoint is w
        w = cs->c.prime_faces[f][if_far];
        a = cs->c.prime_faces[g][ig_v];
        b = cs->c.prime_faces[g][ig_far];
        Lp = lf;
        L = lg;
      } else {
        w = cs->c.prime_faces[g][ig_far];
        a = cs->c.prime_faces[f][if_v];
        b = cs->c.prime_faces[f][if_far];
        Lp = lg;
        L = lf;
      }
      const Vec2 U = B.c(b) - B.c(a), Vv = B.s(b) - B.s(a);
      const Vec2 P = B.c(w) - B.c(a), Q = B.s(w) - B.s(a);
      const Harmonic det = orient_from_vectors(U, Vv, P, Q);
      const Harmonic D = dot_from_vectors(U, Vv, P, Q);
      const Harmonic L2 = dot_from_vectors(U, Vv, U, Vv);
      const double sc = std::abs(det.p) + std::abs(det.q) + std::abs(det.r);
      const double gsc = U.norm() * P.norm();
      if (sc <= 1e-11 * std::max(gsc, 1e-300)) { ++nA_degen; continue; }
      ++nA;
      e_root = std::max(e_root, std::abs(det.eval(beta)) / std::max(sc, 1e-300));
      const double s = D.eval(beta), l2 = L2.eval(beta);
      const double ratio = s / l2;
      worst_ratio_lo = std::min(worst_ratio_lo, ratio);
      worst_ratio_hi = std::max(worst_ratio_hi, ratio);
      if (ratio < 1.0 - 1e-9) ++nA_strict;
      e_ratio = std::max(e_ratio, std::abs(ratio - Lp / L));
      const double tol = 1e-12 * (std::abs(L2.p) + L2.amp());
      if (s < -tol || s > l2 + tol) {
        ++nA_bad;
        std::printf("   [ADM-FAIL] %-18s beta=%.6f  s/|e|^2=%.6f  L'/L=%.6f\n", cs->name.c_str(),
                    beta, ratio, Lp / L);
      }
    }
   }
  }
  std::printf("   hinge edges with both angles    : %lld  (beta_e outside (0,pi): %lld)\n", n_beta_seen, n_beta_out);
  std::printf("   pairs tested                    : %lld  (identically-zero, skipped: %lld)\n", nA,
              nA_degen);
  std::printf("   max |h(beta_e)|/scale           : %.3e   (root of the substitute pair)\n", e_root);
  std::printf("   projection ratio s/|e|^2 range  : [%.6f, %.6f]  (claim: (0, 1])\n",
              worst_ratio_lo, worst_ratio_hi);
  std::printf("   pairs with ratio < 1 - 1e-9    : %lld  (L' < L: vertex STRICTLY inside)\n", nA_strict);
  std::printf("   max |s/|e|^2 - L'/L|            : %.3e   (claim: equal)\n", e_ratio);
  std::printf("   ADMISSIBILITY FAILURES          : %lld\n\n", nA_bad);

  // ==================================================================== B
  // Completeness: a graze inside (0, Theta_max) breaks "Theta_max >= eps => cert".
  std::printf("B  completeness: Theta_max >= eps  =>  certificate ?\n\n");
  std::printf("   %-18s %8s %10s %10s  %10s  %s\n", "pattern", "|C|", "theta_1", "Theta_max", "eps",
              "cert(pos,noov,noroot)");
  int n_counterex = 0;
  for (auto& cs : C) {
    const DeployBasis B = deploy_basis(cs->c, cs->X0);
    const auto sd = swept_discs(cs->c, B);
    const auto pairs = candidate_pairs(cs->c, sd, M_PI, false);
    const auto cand = contact_angles(cs->c, B, pairs, 1e-9, M_PI, 1e-9);
    const auto rep = exact_theta_max_overlap(cs->c, B, pairs, 1e-9, M_PI, 1e-9);
    if (cand.empty()) {
      std::printf("   %-18s %8d %10s %10.6f  %10s  %s\n", cs->name.c_str(), 0, "-", rep.theta_max,
                  "-", "(no contact in range)");
      continue;
    }
    const double th1 = cand.front(), Tm = rep.theta_max;
    if (!(th1 < Tm - 1e-6)) {
      std::printf("   %-18s %8zu %10.6f %10.6f  %10s  %s\n", cs->name.c_str(), cand.size(), th1, Tm,
                  "-", "(no graze: theta_1 == Theta_max)");
      continue;
    }
    const double eps = 0.5 * (th1 + Tm);  // strictly between the graze and the true range
    const auto cert = validity_certificate(cs->c, B, cs->X0, pairs, eps, 1e-9);
    std::printf("   %-18s %8zu %10.6f %10.6f  %10.6f  (%d,%d,%d)%s\n", cs->name.c_str(),
                cand.size(), th1, Tm, eps, (int)cert.pos, (int)cert.nooverlap, (int)cert.noroot,
                (Tm >= eps && !cert.valid()) ? "   <-- COUNTEREXAMPLE to completeness" : "");
    if (Tm >= eps && !cert.valid()) {
      ++n_counterex;
      std::printf("        admissible root at %.6f in (0,eps) though no interior overlap "
                  "until %.6f\n",
                  cert.first_root, Tm);
    }
  }
  std::printf("\n   completeness counterexamples: %d\n", n_counterex);
  std::printf("\nVERDICT  A: %s   B: %s\n", nA_bad == 0 ? "Case A root is ADMISSIBLE (sound)" : "FAIL",
              n_counterex > 0 ? "completeness is FALSE" : "no counterexample found");
  return 0;
}
