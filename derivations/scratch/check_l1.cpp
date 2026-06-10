// check_l1.cpp -- Deriver-L, mission 2 / WP1.  Tests Lemmas L1.1, L1.2, L1.3 of
// derivations/lemmas.md: the COMBINATORIAL classification of contact harmonics.
//
// Build (from repo root):
//   clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
//     -I code/src -I code/apps derivations/scratch/check_l1.cpp \
//     code/src/core/mesh.cpp code/src/core/cut.cpp code/src/core/holes.cpp \
//     code/src/core/kinematics.cpp code/src/core/orientation.cpp code/src/core/generators.cpp \
//     code/src/core/tutte_auxetic.cpp code/src/core/collision.cpp code/src/core/optimize.cpp \
//     -o /tmp/check_l1 && /tmp/check_l1
//
// What is predicted (lemmas.md L1.0-L1.3).  For a candidate (w, (a,b)) with (a,b) an
// edge of face f and w a corner of face g != f, write va, vb, vw for the M-vertices the
// three M'-vertices are copies of, d := x_vb - x_va, and -- when vw == va or vw == vb --
// dS := S_w - S_(the copy of vw inside f).  Then
//
//     C = p + q = det(x_vb - x_va, x_vw - x_va)          (the FLAT determinant)
//
// and, in the coincident-copy case vw in {va, vb},
//
//     C = 0 identically,   p = (sigma_f/2) <d, dS>,   q = -p,   r = (1/2) det(d, dS).
//
// Hence the whole class of such a pair is decided by dS alone:
//     dS = 0                          -> h identically zero      (L1.3)
//     dS != 0, det(d,dS) = 0          -> class 3, never a contact (L1.2)
//     det(d,dS) != 0                  -> class 2                  (L1.1)
// and a pair with vw not in {va, vb} is class 1 unless the three FLAT source vertices
// happen to be collinear -- structurally (frozen triple) or accidentally (isolated t).
//
// The program measures: (i) the combinatorial C=0 predicate against |C| <= tol;
// (ii) the three closed forms above; (iii) the Zero/3/2 split against the numeric class;
// (iv) for every disagreement, whether it persists over all shape-space samples of that
// graph (structural) or occurs at isolated t (accidental).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
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
#include "kill_common.hpp"

using namespace kiri;
using namespace kiri::kill;

static double cross2(const Vec2& a, const Vec2& b) { return a.x() * b.y() - a.y() * b.x(); }
static double dot2(const Vec2& a, const Vec2& b) { return a.x() * b.x() + a.y() * b.y(); }

struct Harm { double p = 0, q = 0, r = 0;
  double scale() const { return std::fabs(p) + std::hypot(q, r); } };

static Harm orient_h(const Vec2& Cab, const Vec2& Sab, const Vec2& Caw, const Vec2& Saw) {
  Harm h;
  const double dcc = cross2(Cab, Caw), dss = cross2(Sab, Saw);
  h.p = 0.5 * (dcc + dss);
  h.q = 0.5 * (dcc - dss);
  h.r = 0.5 * (cross2(Cab, Saw) + cross2(Sab, Caw));
  return h;
}

enum Klass { K_ZERO = 0, K_C1 = 1, K_C2 = 2, K_C3 = 3 };
static const char* kname(int k) {
  return k == K_ZERO ? "zero" : k == K_C1 ? "class1" : k == K_C2 ? "class2" : "class3";
}

// numeric class from (p,q,r) -- the rule the code uses (deploy_basis.hpp)
static Klass numeric_class(const Harm& h, double tol) {
  const double C = h.p + h.q, B = 2 * h.r, A = h.p - h.q;
  const bool c0 = std::fabs(C) <= tol, b0 = std::fabs(B) <= tol, a0 = std::fabs(A) <= tol;
  if (c0 && b0 && a0) return K_ZERO;
  if (!c0) return K_C1;
  return b0 ? K_C3 : K_C2;
}

struct Tally {
  long long pairs = 0;
  // L1.1 : combinatorial C = 0 predicate
  long long shared = 0, czero_numeric = 0;
  long long agree_c = 0;
  long long shared_but_C_nonzero = 0;      // must be 0
  long long cz_not_shared = 0;             // flat-collinear triples (structural or accidental)
  long long cz_not_shared_structural = 0;  // persists over every sample of its graph
  long long cz_not_shared_frozen = 0;      // ... and all three source vertices are frozen
  // L1.2/L1.3 : class prediction on the coincident-copy pairs
  long long agree_class = 0, disagree_class = 0;
  long long n_zero = 0, n_c1 = 0, n_c2 = 0, n_c3 = 0;
  long long c3_persist = 0, c3_transient = 0, c3_single = 0;
  long long cz_ns_single = 0;
  // closed forms
  double max_err_p = 0, max_err_q = 0, max_err_r = 0;
  double max_err_hinge_dS = 0, max_err_split_dS = 0;
  long long n_hinge_dS = 0, n_split_dS = 0;
  // L1.3 structure of the identically-zero pairs
  long long zero_same_prime = 0;
  long long zero_shared = 0, zero_hinge_src = 0, zero_split_zero_du = 0, zero_other = 0;
  long long zero_lambda_interior = 0;      // h == 0 with vw not in {va, vb}
  // shared-edge type of a coincident-copy pair
  long long sh_hinge = 0, sh_split = 0, sh_none = 0;
  // class 3 via the hinge formula:  class3 <=> <d, x_src - x_v> = 0  (d perp the hinge edge)
  long long c3_hinge = 0, c3_hinge_perp_ok = 0;
  long long c3_by_type[3] = {0, 0, 0};
  double max_err_split_dS_const = 0;
  long long n_split_dS_const = 0;
  // persistent non-shared C == 0, by how many of the three source vertices are frozen
  long long cz_ns_frozen_cnt[4] = {0, 0, 0, 0};
  long long ex_printed = 0;
  std::map<std::string, long long> disagreements;
};

static bool positively_oriented(const Mesh& m, const std::vector<Vec2>& X) {
  for (int f = 0; f < m.n_faces(); ++f) {
    double A2 = 0;
    const auto& F = m.faces[f];
    for (size_t i = 0; i < F.size(); ++i)
      A2 += X[F[i]].x() * X[F[(i + 1) % F.size()]].y() -
            X[F[i]].y() * X[F[(i + 1) % F.size()]].x();
    if (A2 <= 0) return false;
  }
  return true;
}

// ---------------------------------------------------------------------------
// One graph, many shape-space samples.  Candidates are enumerated in a FIXED order so
// the per-candidate persistence counters line up across samples.
static void run_graph(const std::string& name, Mesh m, int n_samp, std::mt19937& rng,
                      Tally* T, int* used) {
  m.build_topology();
  if (m.sigma.empty()) return;
  CutStructure c = make_cut(m);
  HoleSet hs = holes_partition(c);
  LinearSystem sys = assemble_system(c, hs, m.X, BoundaryMode::Fixed);
  SolveReport sr = solve_system(sys, m.X);
  if (!sr.projection_ok) return;
  std::vector<Vec2> X0 = matrix_to_points(sr.X0);
  if (!hole_residuals(c, X0, hs).deployable(1e-7)) return;
  if (!positively_oriented(m, X0)) return;
  ++*used;

  const int F = m.n_faces();
  double scale = 0;
  for (auto& p : X0) scale = std::max(scale, p.norm());

  // frozen vertices: zero row of Phi
  std::vector<char> frozen(m.n_vertices(), 1);
  for (int v = 0; v < m.n_vertices(); ++v)
    for (int j = 0; j < sr.dim_null; ++j)
      if (std::fabs(sr.Phi(v, j)) > 1e-12) { frozen[v] = 0; break; }

  // pass 0: count candidates
  long long ncand = 0;
  for (int f = 0; f < F; ++f)
    for (size_t k = 0; k < c.prime_faces[f].size(); ++k)
      for (int g = 0; g < F; ++g)
        if (g != f) ncand += static_cast<long long>(c.prime_faces[g].size());
  std::vector<int> cz_count(ncand, 0), c3_count(ncand, 0), samp_seen(ncand, 0);

  std::normal_distribution<double> gauss(0.0, 1.0);
  int done = 0;
  for (int s = 0; s < n_samp; ++s) {
    std::vector<Vec2> X = X0;
    if (s > 0 && sr.dim_null > 0) {
      Eigen::MatrixXd t(sr.dim_null, 2);
      for (int i = 0; i < sr.dim_null; ++i)
        for (int j = 0; j < 2; ++j) t(i, j) = 0.03 * scale * gauss(rng);
      X = matrix_to_points(Eigen::MatrixXd(sr.X0 + sr.Phi * t));
    } else if (s > 0) {
      break;  // rigid shape space: one sample only
    }
    if (!positively_oriented(m, X)) continue;
    ++done;
    Deployment dp = deploy(c, X, 0.0);
    std::vector<Vec2> Cv = dp.Y, Sv(c.n_prime_vertices);
    for (int i = 0; i < c.n_prime_vertices; ++i) Sv[i] = 2.0 * dp.dY_dtheta[i];
    double sc = 0;
    for (auto& p : X) sc = std::max(sc, p.norm());
    const double tol2 = 1e-11 * sc * sc;   // determinant scale
    const double tol1 = 1e-11 * sc;        // vector scale

    long long idx = -1;
    for (int f = 0; f < F; ++f) {
      const auto& pf = c.prime_faces[f];
      const int nf = static_cast<int>(pf.size());
      for (int kk = 0; kk < nf; ++kk) {
        const int a = pf[kk], b = pf[(kk + 1) % nf];
        const int va = c.prime_to_original[a], vb = c.prime_to_original[b];
        const Vec2 Cab = Cv[b] - Cv[a], Sab = Sv[b] - Sv[a];
        const Vec2 d = X[vb] - X[va];
        for (int g = 0; g < F; ++g) {
          if (g == f) continue;
          for (int w : c.prime_faces[g]) {
            ++idx;
            const int vw = c.prime_to_original[w];
            const Vec2 Caw = Cv[w] - Cv[a], Saw = Sv[w] - Sv[a];
            const Harm h = orient_h(Cab, Sab, Caw, Saw);
            ++T->pairs;
            ++samp_seen[idx];
            const double C = h.p + h.q;
            const bool cz = std::fabs(C) <= tol2;
            const bool sh = (vw == va || vw == vb);
            if (cz) ++cz_count[idx];
            if (sh) ++T->shared;
            if (cz) ++T->czero_numeric;
            if (sh == cz) ++T->agree_c;
            if (sh && !cz) ++T->shared_but_C_nonzero;
            if (!sh && cz) ++T->cz_not_shared;

            const Klass kn = numeric_class(h, tol2);
            switch (kn) {
              case K_ZERO: ++T->n_zero; break;
              case K_C1: ++T->n_c1; break;
              case K_C2: ++T->n_c2; break;
              case K_C3: ++T->n_c3; break;
            }
            if (kn == K_C3) ++c3_count[idx];

            if (!sh) {
              if (kn == K_ZERO) ++T->zero_lambda_interior;
              continue;
            }
            // ---- coincident-copy pair: the closed forms and the predicted class
            const int acopy = (vw == va) ? a : b;
            const Vec2 dS = Sv[w] - Sv[acopy];
            const double sf = m.sigma[f];
            const double pp = 0.5 * sf * dot2(d, dS);
            const double rr = 0.5 * cross2(d, dS);
            T->max_err_p = std::max(T->max_err_p, std::fabs(h.p - pp));
            T->max_err_q = std::max(T->max_err_q, std::fabs(h.q + pp));
            T->max_err_r = std::max(T->max_err_r, std::fabs(h.r - rr));

            Klass kp;
            if (dS.norm() <= tol1) kp = K_ZERO;
            else if (std::fabs(cross2(d, dS)) <= tol2) kp = K_C3;
            else kp = K_C2;
            if (kp == kn) ++T->agree_class;
            else {
              ++T->disagree_class;
              char buf[128];
              std::snprintf(buf, sizeof buf, "%s: pred %s, num %s (shared=%d)",
                            name.c_str(), kname(kp), kname(kn), (int)sh);
              ++T->disagreements[buf];
            }
            {
              int et = -1;
              for (int e : m.vertex_edges[vw]) {
                const auto& E = m.edges[e];
                if (E.n_faces != 2) continue;
                const int fa = m.half_edges[E.he[0]].face, fb = m.half_edges[E.he[1]].face;
                if ((fa == f && fb == g) || (fa == g && fb == f))
                  et = (c.edge_type[e] == EdgeType::Hinge) ? 0
                     : (c.edge_type[e] == EdgeType::Split) ? 1 : 2;
              }
              if (et == 0) ++T->sh_hinge; else if (et == 1) ++T->sh_split; else ++T->sh_none;
              if (kn == K_C3) ++T->c3_by_type[et < 0 ? 2 : et];
            }
            if (kn == K_ZERO) {
              ++T->zero_shared;
              if (w == a || w == b) ++T->zero_same_prime;
              // which structural reason?
              int e_shared = -1;
              for (int e : m.vertex_edges[vw]) {
                const auto& E = m.edges[e];
                if (E.n_faces != 2) continue;
                const int fa = m.half_edges[E.he[0]].face, fb = m.half_edges[E.he[1]].face;
                if ((fa == f && fb == g) || (fa == g && fb == f)) e_shared = e;
              }
              if (e_shared >= 0 && c.edge_type[e_shared] == EdgeType::Hinge &&
                  c.hinge_dir[e_shared].src == vw) ++T->zero_hinge_src;
              else if (e_shared >= 0 && c.edge_type[e_shared] == EdgeType::Split)
                ++T->zero_split_zero_du;
              else ++T->zero_other;
            }
            // ---- dS closed forms across a shared edge
            int e_sh = -1;
            for (int e : m.vertex_edges[vw]) {
              const auto& E = m.edges[e];
              if (E.n_faces != 2) continue;
              const int fa = m.half_edges[E.he[0]].face, fb = m.half_edges[E.he[1]].face;
              if ((fa == f && fb == g) || (fa == g && fb == f)) e_sh = e;
            }
            if (e_sh >= 0 && c.edge_type[e_sh] == EdgeType::Split) {
              // T1.B: dS_e is the SAME vector for both endpoints of the split edge
              const auto& E = m.edges[e_sh];
              const int vo = (E.key.a == vw) ? E.key.b : E.key.a;
              const int wo = c.prime_vertex(g, vo), ao = c.prime_vertex(f, vo);
              if (wo >= 0 && ao >= 0) {
                const Vec2 dS2 = Sv[wo] - Sv[ao];
                T->max_err_split_dS_const =
                    std::max(T->max_err_split_dS_const, (dS - dS2).norm());
                ++T->n_split_dS_const;
              }
            }
            if (e_sh >= 0 && c.edge_type[e_sh] == EdgeType::Hinge) {
              const int src = c.hinge_dir[e_sh].src;
              const double sg = m.sigma[g];
              Vec2 pred(0, 0);
              const Vec2 dd = X[src] - X[vw];
              pred = Vec2(-dd.y(), dd.x()) * (2.0 * sg);
              T->max_err_hinge_dS = std::max(T->max_err_hinge_dS, (dS - pred).norm());
              ++T->n_hinge_dS;
              if (kn == K_C3) {
                ++T->c3_hinge;
                if (std::fabs(dot2(d, dd)) <= tol2) ++T->c3_hinge_perp_ok;
              }
            }
          }
        }
      }
    }
  }
  // persistence bookkeeping
  for (long long i = 0; i < ncand; ++i) {
    if (samp_seen[i] == 0) continue;
    if (cz_count[i] > 0 && cz_count[i] == samp_seen[i]) { /* structural, counted below */ }
  }
  // recount: which of the "C=0 but not shared" are structural, and frozen
  {
    long long idx = -1;
    for (int f = 0; f < F; ++f) {
      const auto& pf = c.prime_faces[f];
      const int nf = static_cast<int>(pf.size());
      for (int kk = 0; kk < nf; ++kk) {
        const int a = pf[kk], b = pf[(kk + 1) % nf];
        const int va = c.prime_to_original[a], vb = c.prime_to_original[b];
        for (int g = 0; g < F; ++g) {
          if (g == f) continue;
          for (int w : c.prime_faces[g]) {
            ++idx;
            const int vw = c.prime_to_original[w];
            const bool sh = (vw == va || vw == vb);
            if (sh || samp_seen[idx] == 0) continue;
            if (samp_seen[idx] < 5) {                 // rigid graph: one t only, cannot decide
              T->cz_ns_single += cz_count[idx];
              T->c3_single += c3_count[idx];
              continue;
            }
            if (cz_count[idx] == samp_seen[idx]) {
              T->cz_not_shared_structural += samp_seen[idx];
              const int nf3 = (frozen[va] ? 1 : 0) + (frozen[vb] ? 1 : 0) + (frozen[vw] ? 1 : 0);
              T->cz_ns_frozen_cnt[nf3] += samp_seen[idx];
              if (nf3 == 3) T->cz_not_shared_frozen += samp_seen[idx];
              if (nf3 < 3 && T->ex_printed < 8) {
                ++T->ex_printed;
                std::printf("    [example] %s  f=%d g=%d  va=%d(%s) vb=%d(%s) vw=%d(%s)  "
                            "bnd=%d%d%d\n", name.c_str(), f, g, va, frozen[va] ? "frz" : "free",
                            vb, frozen[vb] ? "frz" : "free", vw, frozen[vw] ? "frz" : "free",
                            (int)m.vertex_is_boundary[va], (int)m.vertex_is_boundary[vb],
                            (int)m.vertex_is_boundary[vw]);
              }
            }
            if (c3_count[idx] > 0) {
              if (c3_count[idx] == samp_seen[idx]) T->c3_persist += c3_count[idx];
              else T->c3_transient += c3_count[idx];
            }
          }
        }
      }
    }
  }
  // class-3 persistence for shared pairs too
  {
    long long idx = -1;
    for (int f = 0; f < F; ++f) {
      const auto& pf = c.prime_faces[f];
      const int nf = static_cast<int>(pf.size());
      for (int kk = 0; kk < nf; ++kk) {
        const int a = pf[kk], b = pf[(kk + 1) % nf];
        const int va = c.prime_to_original[a], vb = c.prime_to_original[b];
        for (int g = 0; g < F; ++g) {
          if (g == f) continue;
          for (int w : c.prime_faces[g]) {
            ++idx;
            const int vw = c.prime_to_original[w];
            if (!(vw == va || vw == vb) || samp_seen[idx] == 0) continue;
            if (samp_seen[idx] < 5) { T->c3_single += c3_count[idx]; continue; }
            if (c3_count[idx] > 0) {
              if (c3_count[idx] == samp_seen[idx]) T->c3_persist += c3_count[idx];
              else T->c3_transient += c3_count[idx];
            }
          }
        }
      }
    }
  }
  (void)done;
}

int main(int argc, char** argv) {
  const int NRAND = (argc > 1) ? std::atoi(argv[1]) : 52;
  const int NSAMP = (argc > 2) ? std::atoi(argv[2]) : 20;
  std::mt19937 rng(20260908);
  Tally T;
  int used = 0;

  for (auto& rc : reference_cases()) {
    if (rc.mesh.sigma.empty()) continue;
    run_graph("ref:" + rc.name, rc.mesh, NSAMP, rng, &T, &used);
  }
  for (int id = 0; id < NRAND; ++id) {
    Graph g = make_graph(id, 18, 46, 220);
    if (!g.ok) continue;
    run_graph("rand:" + g.kind + std::to_string(id), g.mesh, NSAMP, rng, &T, &used);
  }

  std::printf("\ncheck_l1 -- graphs used: %d  (random requested %d, samples/graph %d)\n",
              used, NRAND, NSAMP);
  std::printf("  candidate harmonics examined            : %lld\n", T.pairs);
  std::printf("\nL1.1  combinatorial predicate  C == 0  <=>  vw in {va, vb}\n");
  std::printf("  coincident-copy pairs (vw in {va,vb})   : %lld\n", T.shared);
  std::printf("  numerically |C| <= tol                  : %lld\n", T.czero_numeric);
  std::printf("  AGREEMENTS                              : %lld / %lld\n", T.agree_c, T.pairs);
  std::printf("  shared but C != 0  (must be 0)          : %lld\n", T.shared_but_C_nonzero);
  std::printf("  C == 0 but NOT shared (flat-collinear)  : %lld\n", T.cz_not_shared);
  std::printf("     ... on rigid graphs (dim_null = 0, one sample, undecidable) : %lld\n",
              T.cz_ns_single);
  std::printf("     ... of these, persistent over every sample of the graph : %lld\n",
              T.cz_not_shared_structural);
  std::printf("     ... of these, all three source vertices frozen          : %lld\n",
              T.cz_not_shared_frozen);
  std::printf("\nL1.2/L1.3  class of a coincident-copy pair from dS alone\n");
  std::printf("  numeric classes over ALL pairs: zero %lld  class1 %lld  class2 %lld  class3 %lld\n",
              T.n_zero, T.n_c1, T.n_c2, T.n_c3);
  std::printf("  predicted == numeric on coincident-copy pairs : %lld agree, %lld disagree\n",
              T.agree_class, T.disagree_class);
  std::printf("  closed forms  max |p - sf<d,dS>/2| = %.3e   max |q + p| = %.3e   "
              "max |r - det(d,dS)/2| = %.3e\n", T.max_err_p, T.max_err_q, T.max_err_r);
  std::printf("  hinge dS = 2 sigma_g J (x_src - x_v):  %lld pairs, max err %.3e\n",
              T.n_hinge_dS, T.max_err_hinge_dS);
  std::printf("  class 3 occurrences: persistent %lld   transient (isolated t) %lld   "
              "undecidable (rigid graph, 1 sample) %lld\n",
              T.c3_persist, T.c3_transient, T.c3_single);
  std::printf("\nL1.3  the identically-zero pairs\n");
  std::printf("  zero & coincident-copy                  : %lld\n", T.zero_shared);
  std::printf("     ... of which w IS the same M'-vertex as a or b : %lld\n",
              T.zero_same_prime);
  std::printf("     shared edge is a HINGE with src == v : %lld\n", T.zero_hinge_src);
  std::printf("     shared edge is a SPLIT (Delta u = 0) : %lld\n", T.zero_split_zero_du);
  std::printf("     no shared edge between f and g       : %lld\n", T.zero_other);
  std::printf("  zero with vw NOT in {va,vb} (lambda interior) : %lld\n",
              T.zero_lambda_interior);
  std::printf("\n  coincident-copy pairs by shared-edge type: hinge %lld  split %lld  none %lld\n",
              T.sh_hinge, T.sh_split, T.sh_none);
  std::printf("  class-3 pairs by shared-edge type: hinge %lld  split %lld  none %lld\n",
              T.c3_by_type[0], T.c3_by_type[1], T.c3_by_type[2]);
  std::printf("  class-3 pairs across a HINGE: %lld, of which <d, x_src - x_v> = 0 : %lld\n",
              T.c3_hinge, T.c3_hinge_perp_ok);
  std::printf("  split dS constant over the two endpoints (T1.B): %lld pairs, max err %.3e\n",
              T.n_split_dS_const, T.max_err_split_dS_const);
  std::printf("  persistent non-shared C==0, by #frozen source vertices: "
              "0:%lld  1:%lld  2:%lld  3:%lld\n",
              T.cz_ns_frozen_cnt[0], T.cz_ns_frozen_cnt[1], T.cz_ns_frozen_cnt[2],
              T.cz_ns_frozen_cnt[3]);
  if (!T.disagreements.empty()) {
    std::printf("\n  DISAGREEMENT TYPES:\n");
    for (auto& kv : T.disagreements)
      std::printf("    %-70s  x %lld\n", kv.first.c_str(), kv.second);
  }
  std::printf("\n");
  return 0;
}
