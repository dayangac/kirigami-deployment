// check_b4_93.cpp -- Checker round 6, add-on: the B4 `voronoi_93` disagreement.
//
// results/kill/b4/b4.csv, row (id 93, voronoi, sigma_def, free):
//     cert_pos=1  cert_noovl=1  cert_noroot=0  eps_max=0.2484
//     theta_exact=0.2484   theta_bisect=0
// while (id 96, ..., free) has theta_exact == theta_bisect == 0.241884.
//
// This program replays both rows through kill_b4.cpp's own pipeline (same graph, same
// sigma_def from results/kill/k5/sigma, same `free` variant system, same repair seeds and
// weights) and then asks WHICH predicate separates the two answers:
//
//   1. the BROAD PHASE.  exact_theta_max_overlap is run on candidate_pairs(..., prune=true)
//      while the referee's has_collision() tests EVERY face pair.  Re-run the scan with
//      prune=false: if theta_exact collapses, a pair the grid dropped is the culprit.
//   2. the SHRINK.  the scan is called with shrink = 1e-9, the referee with 1e-12.
//   3. the theta = 0+ PROBE.  referee_theta returns 0 as soon as has_collision fires at
//      theta = 1e-7; the scan's interval sweep starts at theta_lo = 1e-9 and decides the
//      first slab by an overlap probe at its MIDPOINT.
//
// Build (from repo root; libkiri_core.a must exist):
//   clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
//     -Icode/src -Icode/apps derivations/scratch/check_b4_93.cpp code/build/libkiri_core.a \
//     -o /tmp/check_b4_93 && /tmp/check_b4_93
#include <cstdio>
#include <string>
#include <vector>

#include "kill_common.hpp"
#include "method/contact.hpp"
#include "method/zero_plus.hpp"

using namespace kiri;
using namespace kiri::kill;
using namespace kiri::method;

namespace {

// kill_b4.cpp's referee_theta, verbatim.
double referee_theta(const CutStructure& c, const std::vector<Vec2>& X, int grid = 4000,
                     int iters = 50) {
  auto col = [&](double th) { return has_collision(c, deploy(c, X, th).Y, 1e-12); };
  if (col(1e-7)) return 0.0;
  double lo = 0, hi = -1;
  for (int i = 1; i <= grid; ++i) {
    const double th = M_PI * i / grid;
    if (col(th)) { hi = th; break; }
    lo = th;
  }
  if (hi < 0) return M_PI;
  for (int i = 0; i < iters; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (col(mid)) hi = mid; else lo = mid;
  }
  return lo;
}

// kill_b4.cpp's `free` variant: drop every Eq. (4) boundary row, pin vertex 0.
Shape shape_free(const Mesh& m, const CutStructure& c, const HoleSet& hs,
                 const std::string& cache_dir, int id) {
  Shape s;
  if (!cache_dir.empty() && load_shape(cache_dir, id, &s)) {
    if (s.N == m.n_vertices()) return s;
    s = Shape();
  }
  LinearSystem sys = assemble_system(c, hs, m.X, BoundaryMode::None);
  const int nh = static_cast<int>(sys.A.rows());
  const int N = sys.N;
  Eigen::MatrixXd A(nh + 1, N), rhs(nh + 1, 2);
  A.setZero();
  A.topRows(nh) = sys.A;
  rhs.topRows(nh) = sys.rhs;
  A(nh, 0) = 1.0;
  rhs(nh, 0) = m.X[0].x();
  rhs(nh, 1) = m.X[0].y();
  sys.A = std::move(A);
  sys.rhs = std::move(rhs);
  sys.n_boundary_rows = 1;
  const SolveReport sr = solve_system(sys, m.X);
  s.N = N; s.k = sr.dim_null; s.X0 = sr.X0; s.Phi = sr.Phi;
  s.rank_L = sr.rank_L; s.H = sr.H; s.n_interior = m.n_interior_vertices();
  s.ok = sr.projection_ok;
  return s;
}

inline double cross2(const Vec2& a, const Vec2& b) { return a.x() * b.y() - a.y() * b.x(); }
// proper (open) segment crossing test, for the simplicity screen
bool seg_cross(const Vec2& p1, const Vec2& p2, const Vec2& q1, const Vec2& q2) {
  const double d1 = cross2(p2 - p1, q1 - p1), d2 = cross2(p2 - p1, q2 - p1);
  const double d3 = cross2(q2 - q1, p1 - q1), d4 = cross2(q2 - q1, p2 - q1);
  return ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0));
}

// which face pair collides at theta, over ALL pairs, at the given shrink
std::pair<int, int> colliding_pair(const CutStructure& c, const std::vector<Vec2>& Y,
                                   double shrink) {
  const int F = c.mesh->n_faces();
  std::vector<Vec2> A, B;
  for (int f = 0; f < F; ++f) {
    A.clear();
    for (int pv : c.prime_faces[f]) A.push_back(Y[pv]);
    for (int g = f + 1; g < F; ++g) {
      B.clear();
      for (int pv : c.prime_faces[g]) B.push_back(Y[pv]);
      if (polygons_overlap(A, B, shrink)) return {f, g};
    }
  }
  return {-1, -1};
}

void replay(int id) {
  std::printf("================ id %d ================\n", id);
  Graph g = make_graph(id, 100, 800, 1400);
  if (!g.ok) { std::printf("  graph not ok\n"); return; }
  Mesh m0 = g.mesh;
  m0.build_topology();
  const std::string p = "results/kill/k5/sigma/" + g.kind + "_" + std::to_string(id) + ".json";
  if (!std::filesystem::exists(p)) { std::printf("  no sigma_def at %s\n", p.c_str()); return; }
  const Mesh sm = load_mesh_json(p);
  Mesh m = m0;
  m.sigma = sm.sigma;
  m.build_topology();
  const CutStructure c = make_cut(m);
  const HoleSet hs = holes_partition(c);
  const Shape sh = shape_free(m, c, hs, "results/kill/b4/cache/free_sigma_def", id);
  if (!sh.ok || sh.k < 1) { std::printf("  shape not ok\n"); return; }
  const std::vector<Vec2> X0 = matrix_to_points(sh.X0);
  const double med = median_edge_length(m);
  const unsigned seed0 = 6000u + 7u * static_cast<unsigned>(id) + 1u;  // which = 1 (sigma_def)

  ZeroPlusRepairOptions ro;
  ro.n_random = 3; ro.max_iter = 1200; ro.lambda_rel = 1e-6; ro.seed = seed0;
  const ZeroPlusRepairResult sr = zero_plus_repair(c, X0, sh.Phi, med, ro);
  ZeroPlusRepairOptions ro2 = ro;
  ro2.w_corner = 0.05; ro2.w_prox = 1e3;
  if (sr.feasible) ro2.t_init = sr.t;
  const ZeroPlusRepairResult rr = zero_plus_repair(c, X0, sh.Phi, med, ro2);
  std::printf("  repair feasible=%d min_q=%.6g  (kind=%s F=%d N=%d)\n", (int)rr.feasible, rr.min_q,
              g.kind.c_str(), m.n_faces(), m.n_vertices());
  if (!rr.feasible) return;
  const std::vector<Vec2>& X = rr.X;

  const DeployBasis B = deploy_basis(c, X);
  const SweptDiscs sd = swept_discs(c, B);
  const auto pruned = candidate_pairs(c, sd, M_PI, true);
  const auto allp = candidate_pairs(c, sd, M_PI, false);
  const auto r_pruned = exact_theta_max_overlap(c, B, pruned, 1e-9, M_PI, 1e-9);
  const auto r_all = exact_theta_max_overlap(c, B, allp, 1e-9, M_PI, 1e-9);
  const double tb = referee_theta(c, X);
  std::printf("  pairs: pruned %zu of %zu\n", pruned.size(), allp.size());
  std::printf("  Theta_max  scan(pruned, shrink 1e-9) = %.6f   |C| = %zu\n", r_pruned.theta_max,
              r_pruned.candidates.size());
  std::printf("  Theta_max  scan(ALL   , shrink 1e-9) = %.6f   |C| = %zu\n", r_all.theta_max,
              r_all.candidates.size());
  std::printf("  Theta_max  referee bisection          = %.6f\n", tb);

  // where does has_collision fire, and on which pair?
  const double probes[] = {1e-7, 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.2484, 0.3};
  for (double th : probes) {
    const auto Y = deploy(c, X, th).Y;
    const bool hc = has_collision(c, Y, 1e-12);
    const auto p12 = colliding_pair(c, Y, 1e-12);
    const auto p9 = colliding_pair(c, Y, 1e-9);
    bool in_pruned = false;
    if (p12.first >= 0)
      for (auto& q : pruned)
        if (q.first == p12.first && q.second == p12.second) in_pruned = true;
    std::printf("    th=%-9.6g has_collision(1e-12)=%d  pair@1e-12=(%d,%d) in_pruned=%d  "
                "pair@1e-9=(%d,%d)\n",
                th, (int)hc, p12.first, p12.second, (int)in_pruned, p9.first, p9.second);
  }

  // If a pair collides at some theta below theta_exact, report its scan status.
  const auto Yc = deploy(c, X, 1e-7).Y;
  const auto bad = colliding_pair(c, Yc, 1e-12);
  if (bad.first >= 0) {
    std::vector<std::pair<int, int>> one{bad};
    const auto ca = contact_angles(c, B, one, 1e-9, M_PI, 1e-9);
    std::printf("  offending pair (%d,%d) at theta=1e-7: |C(pair)| = %zu, first = %s\n", bad.first,
                bad.second, ca.size(), ca.empty() ? "none" : std::to_string(ca.front()).c_str());
    // is it a hinge-adjacent / shared-vertex pair?  (a permanent incidence)
    int shared = 0;
    for (int u : c.prime_faces[bad.first])
      for (int v : c.prime_faces[bad.second])
        if (u == v) ++shared;
    std::printf("  offending pair shares %d M'-vertices (hinge point => boundary contact, "
                "not interior overlap)\n", shared);
    for (double s : {1e-12, 1e-9, 1e-7, 1e-6}) {
      std::vector<Vec2> A2, B2;
      for (int pv : c.prime_faces[bad.first]) A2.push_back(Yc[pv]);
      for (int pv : c.prime_faces[bad.second]) B2.push_back(Yc[pv]);
      std::printf("    polygons_overlap(shrink=%.0e) = %d\n", s,
                  (int)polygons_overlap(A2, B2, s));
    }
  }
  // ---- deep dive on the offending pair -----------------------------------------
  if (bad.first >= 0) {
    const int f = bad.first, gg = bad.second;
    auto ov = [&](double th) {
      const auto Y = deploy(c, X, th).Y;
      std::vector<Vec2> A2, B2;
      for (int pv : c.prime_faces[f]) A2.push_back(Y[pv]);
      for (int pv : c.prime_faces[gg]) B2.push_back(Y[pv]);
      return polygons_overlap(A2, B2, 1e-12);
    };
    // bisect the end of the overlap
    double lo = 1e-9, hi = 0.1;
    if (ov(lo) && !ov(hi)) {
      for (int i = 0; i < 60; ++i) { const double mid = 0.5 * (lo + hi); if (ov(mid)) lo = mid; else hi = mid; }
      std::printf("  pair (%d,%d): interiors overlap on (0, %.9f), clear above\n", f, gg, hi);
      // is that angle in C(X) for this pair?  and for ANY pair?
      std::vector<std::pair<int,int>> one{{f, gg}};
      const auto ca1 = contact_angles(c, B, one, 1e-9, 0.5, 1e-9);
      std::printf("  C(pair) on (0,0.5): %zu angles ->", ca1.size());
      for (double v : ca1) std::printf(" %.6f", v);
      std::printf("\n");
      int n_near = 0;
      for (double v : r_all.candidates) if (std::abs(v - hi) < 1e-4) ++n_near;
      std::printf("  C(X) angles within 1e-4 of %.9f (over ALL pairs): %d;  min C(X) = %.9f\n",
                  hi, n_near, r_all.candidates.empty() ? -1.0 : r_all.candidates.front());
      // simplicity of the two faces at theta = 0 and at the transition
      auto simple_at = [&](int ff, double th) {
        const auto Y = deploy(c, X, th).Y;
        std::vector<Vec2> P;
        for (int pv : c.prime_faces[ff]) P.push_back(Y[pv]);
        const int n = (int)P.size();
        for (int i = 0; i < n; ++i)
          for (int j = i + 1; j < n; ++j) {
            if (j == i + 1 || (i == 0 && j == n - 1)) continue;
            if (seg_cross(P[i], P[(i+1)%n], P[j], P[(j+1)%n])) return false;
          }
        return true;
      };
      std::printf("  face %d simple at 1e-7/mid/0.1: %d %d %d ; face %d: %d %d %d\n",
                  f, (int)simple_at(f,1e-7), (int)simple_at(f,0.5*hi), (int)simple_at(f,0.1),
                  gg, (int)simple_at(gg,1e-7), (int)simple_at(gg,0.5*hi), (int)simple_at(gg,0.1));
      // signed areas (POS) of the two faces
      Mesh t = m; t.X = X;
      std::printf("  flat signed areas: face %d = %.6g, face %d = %.6g\n",
                  f, t.face_signed_area(f), gg, t.face_signed_area(gg));
      // does the exact scan's own candidate list contain ANY angle below hi?
      int below = 0;
      for (double v : r_all.candidates) if (v < hi) ++below;
      std::printf("  C(X) angles strictly below %.9f: %d\n", hi, below);

      // ---- WHO is the witness at the transition, and why did contact_angles miss it?
      const double thc = hi;
      const auto Yt = deploy(c, X, thc).Y;
      double gsc = 0;
      for (const auto& y : Yt) gsc = std::max(gsc, y.norm());
      auto scan_pair = [&](int fe, int fv) {
        const auto& E = c.prime_faces[fe];
        const auto& V = c.prime_faces[fv];
        const int ne = (int)E.size();
        for (int i = 0; i < ne; ++i) {
          const int a = E[i], b = E[(i + 1) % ne];
          const Vec2 U = B.c(b) - B.c(a), Vv = B.s(b) - B.s(a);
          const Harmonic L2 = dot_from_vectors(U, Vv, U, Vv);
          for (int pv : V) {
            const Vec2 P = B.c(pv) - B.c(a), Q = B.s(pv) - B.s(a);
            const Harmonic det = orient_from_vectors(U, Vv, P, Q);
            const Harmonic D = dot_from_vectors(U, Vv, P, Q);
            const double l2 = L2.eval(thc), sdot = D.eval(thc);
            const double geo = std::sqrt(std::max(l2, 0.0)) * (Yt[pv] - Yt[a]).norm();
            const double hval = det.eval(thc);
            const double dist = (geo > 0) ? std::abs(hval) / std::sqrt(std::max(l2, 1e-300)) : 0.0;
            const bool on_seg = (sdot >= -1e-9 * l2 && sdot <= l2 * (1 + 1e-9));
            if (dist < 1e-6 * gsc && on_seg) {
              const bool skipped_endpoint = (pv == a || pv == b);
              const bool ident_zero = (std::abs(det.p) + det.amp() <= 0);
              const double idsc = std::abs(det.p) + std::abs(det.q) + std::abs(det.r);
              std::vector<std::pair<int,int>> one2{{fe, fv}};
              std::printf("    witness edge(%d,%d) of face %d  vertex %d of face %d : "
                          "dist=%.3e s/|e|^2=%.6f  skipped(p==a||p==b)=%d ident_zero=%d "
                          "|p|+|q|+|r|=%.3e  h(0)=%.3e\n",
                          a, b, fe, pv, fv, dist, sdot / l2, (int)skipped_endpoint,
                          (int)ident_zero, idsc, det.eval(0.0));
              if (!skipped_endpoint && !ident_zero) {
                const auto rts = harmonic_roots_deflated(det, 1e-9, M_PI);
                std::printf("      deflated roots in (0,pi):");
                for (double rv : rts) std::printf(" %.6f", rv);
                std::printf("   (undeflated:");
                for (double rv : harmonic_roots(det, 1e-9, M_PI)) std::printf(" %.6f", rv);
                std::printf(")\n");
              }
            }
          }
        }
      };
      scan_pair(f, gg);
      scan_pair(gg, f);

      // ---- WHERE is the overlap?  sample points strictly inside both.
      auto pip = [&](const std::vector<Vec2>& P, const Vec2& z) {
        bool in = false;
        const int n = (int)P.size();
        for (int i = 0, j = n - 1; i < n; j = i++)
          if (((P[i].y() > z.y()) != (P[j].y() > z.y())) &&
              (z.x() < (P[j].x() - P[i].x()) * (z.y() - P[i].y()) / (P[j].y() - P[i].y()) + P[i].x()))
            in = !in;
        return in;
      };
      auto edge_dist = [&](const std::vector<Vec2>& P, const Vec2& z) {
        double d = 1e300;
        const int n = (int)P.size();
        for (int i = 0; i < n; ++i) {
          const Vec2 a2 = P[i], b2 = P[(i + 1) % n];
          const Vec2 e2 = b2 - a2;
          const double t2 = std::max(0.0, std::min(1.0, (z - a2).dot(e2) / std::max(e2.squaredNorm(), 1e-300)));
          d = std::min(d, (z - (a2 + t2 * e2)).norm());
        }
        return d;
      };
      int shared_pv = -1;
      for (int u : c.prime_faces[f]) for (int v2 : c.prime_faces[gg]) if (u == v2) shared_pv = u;
      for (double th : {1e-4, 0.01, 0.03, 0.036, 0.04}) {
        const auto Y2 = deploy(c, X, th).Y;
        std::vector<Vec2> A2, B2;
        for (int pv : c.prime_faces[f]) A2.push_back(Y2[pv]);
        for (int pv : c.prime_faces[gg]) B2.push_back(Y2[pv]);
        Vec2 lo2(1e300, 1e300), hi2(-1e300, -1e300);
        for (auto& z : A2) { lo2 = lo2.cwiseMin(z); hi2 = hi2.cwiseMax(z); }
        double best = -1; Vec2 bz(0, 0);
        const int G = 400;
        for (int i = 0; i <= G; ++i) for (int j = 0; j <= G; ++j) {
          const Vec2 z(lo2.x() + (hi2.x() - lo2.x()) * i / G, lo2.y() + (hi2.y() - lo2.y()) * j / G);
          if (!pip(A2, z) || !pip(B2, z)) continue;
          const double d = std::min(edge_dist(A2, z), edge_dist(B2, z));
          if (d > best) { best = d; bz = z; }
        }
        int nA = 0, nB = 0, nAB = 0;
        for (int i = 0; i <= G; ++i) for (int j = 0; j <= G; ++j) {
          const Vec2 z(lo2.x() + (hi2.x() - lo2.x()) * i / G, lo2.y() + (hi2.y() - lo2.y()) * j / G);
          const bool ia = pip(A2, z), ib = pip(B2, z);
          if (ia) ++nA; if (ib) ++nB; if (ia && ib) ++nAB;
        }
        // radial probe around the shared hinge point, like the R4-a local test
        int nloc = 0;
        if (shared_pv >= 0) {
          const Vec2 pz = Y2[shared_pv];
          double r0 = 1e300;
          for (auto& z : A2) if ((z - pz).norm() > 1e-14) r0 = std::min(r0, (z - pz).norm());
          for (int k = 1; k <= 6; ++k) {
            const double rad = 0.4 * r0 * k / 7.0;
            for (int aa = 0; aa < 2880; ++aa) {
              const double ph = 2 * M_PI * aa / 2880;
              const Vec2 z = pz + rad * Vec2(std::cos(ph), std::sin(ph));
              if (pip(A2, z) && pip(B2, z)) ++nloc;
            }
          }
        }
        std::printf("      grid inside A=%d B=%d BOTH=%d ; radial probes at the hinge inside BOTH=%d\n",
                    nA, nB, nAB, nloc);
        // which edge pair does the shrunk-polygon test call a crossing?
        for (double sk : {0.0, 1e-9}) {
          std::vector<Vec2> As, Bs;
          {
            Vec2 ca = Vec2::Zero(); for (auto& z : A2) ca += z; ca /= (double)A2.size();
            for (auto& z : A2) As.push_back(ca + (z - ca) * (1.0 - sk));
            Vec2 cb = Vec2::Zero(); for (auto& z : B2) cb += z; cb /= (double)B2.size();
            for (auto& z : B2) Bs.push_back(cb + (z - cb) * (1.0 - sk));
          }
          const int na = (int)As.size(), nb = (int)Bs.size();
          for (int i = 0; i < na; ++i) for (int j = 0; j < nb; ++j) {
            const Vec2 &a3 = As[i], &b3 = As[(i+1)%na], &c3 = Bs[j], &d3 = Bs[(j+1)%nb];
            const double e1 = cross2(b3-a3, c3-a3), e2 = cross2(b3-a3, d3-a3);
            const double e3 = cross2(d3-c3, a3-c3), e4 = cross2(d3-c3, b3-c3);
            if (((e1>0&&e2<0)||(e1<0&&e2>0)) && ((e3>0&&e4<0)||(e3<0&&e4>0))) {
              const int pa = c.prime_faces[f][i], pb = c.prime_faces[f][(i+1)%na];
              const int pc = c.prime_faces[gg][j], pd = c.prime_faces[gg][(j+1)%nb];
              std::printf("      shrink=%.0e SPURIOUS crossing: A edge(%d,%d) x B edge(%d,%d) "
                          "crosses=%.2e %.2e %.2e %.2e  hinge-incident=%d\n", sk, pa, pb, pc, pd,
                          e1, e2, e3, e4,
                          (int)(pa==shared_pv||pb==shared_pv||pc==shared_pv||pd==shared_pv));
            }
          }
        }
        const double dh = (shared_pv >= 0) ? (bz - Y2[shared_pv]).norm() : -1;
        std::printf("  th=%-8.5g deepest common interior point: depth=%.4e  dist to hinge=%.4e "
                    "(hinge pv=%d)  polygons_overlap(1e-12)=%d\n",
                    th, best, best > 0 ? dh : -1.0, shared_pv, (int)polygons_overlap(A2, B2, 1e-12));
        std::printf("      polygons_overlap shrink 0 / 1e-9 / 1e-6 / 1e-4 = %d %d %d %d\n",
                    (int)polygons_overlap(A2,B2,0.0), (int)polygons_overlap(A2,B2,1e-9),
                    (int)polygons_overlap(A2,B2,1e-6), (int)polygons_overlap(A2,B2,1e-4));
      }
      // ---- dump the two polygons at theta = 0.01 as C++ literals, for a deterministic
      //      regression case in derivation_tests.cpp (no cache, no repair, no generator).
      {
        const auto Yd = deploy(c, X, 0.01).Y;
        for (int ff : {f, gg}) {
          std::printf("  // face %d (%zu vertices), theta = 0.01\n  const std::vector<Vec2> P%d = {",
                      ff, c.prime_faces[ff].size(), ff == f ? 1 : 2);
          int k2 = 0;
          for (int pv : c.prime_faces[ff]) {
            std::printf("%s{%.17g, %.17g}", (k2 ? ", " : ""), Yd[pv].x(), Yd[pv].y());
            ++k2;
          }
          std::printf("};\n");
        }
      }

      // hinge angles and beta_e
      if (shared_pv >= 0) {
        auto ang_at = [&](int ff, double th) {
          const auto Y2 = deploy(c, X, th).Y;
          const auto& PF = c.prime_faces[ff];
          const int n = (int)PF.size();
          int i = -1;
          for (int j = 0; j < n; ++j) if (PF[j] == shared_pv) i = j;
          if (i < 0) return -1.0;
          const Vec2 a2 = Y2[PF[(i + n - 1) % n]] - Y2[shared_pv], b2 = Y2[PF[(i + 1) % n]] - Y2[shared_pv];
          double t2 = std::atan2(cross2(b2, a2), a2.dot(b2));
          if (t2 < 0) t2 += 2 * M_PI;
          return t2;
        };
        const double af = ang_at(f, 0.01), ag2 = ang_at(gg, 0.01);
        std::printf("  hinge angles at th=0.01: alpha_%d=%.6f alpha_%d=%.6f  beta_e=%.6f\n",
                    f, af, gg, ag2, 2 * M_PI - af - ag2);
      }
    } else {
      std::printf("  pair (%d,%d): ov(1e-9)=%d ov(0.1)=%d (not the simple on/off shape)\n",
                  f, gg, (int)ov(lo), (int)ov(hi));
    }
  }
  std::printf("\n");
}

}  // namespace

int main(int argc, char** argv) {
  if (argc > 1) {
    for (int i = 1; i < argc; ++i) replay(std::stoi(argv[i]));
  } else {
    replay(93);
    replay(96);
  }
  return 0;
}
