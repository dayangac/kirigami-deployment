// check_l2.cpp -- Deriver-L, mission 2 / WP1.  Tests Lemmas L2.1 and L2.2 of
// derivations/lemmas.md: the periodic dimension formula dim K = 2 rank(D), and the
// explicit factorisation of the generators M_j through the period-potential row d_i.
//
// Build (from repo root):
//   clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
//     -I code/src -I code/apps derivations/scratch/check_l2.cpp \
//     code/src/core/*.cpp code/src/method/deploy_basis.cpp code/src/method/contact.cpp \
//     code/src/method/periodic_jacobian.cpp code/src/method/zero_plus.cpp -o /tmp/check_l2 \
//     && /tmp/check_l2
//
// The population, the quotient/super-patch pipeline and `achievable()` are copied
// VERBATIM from code/apps/kill_k7.cpp (lines 36-243) so that the numbers compare
// directly with results/kill/k7/k7_main.csv.  Nothing in code/ is modified.
//
// What is predicted (lemmas.md L2.1).  Write d_i in R^2 for row i of D (the pair of
// period-potential increments of the i-th null vector phi_i) and P0 = T.  Then
//
//     M_{2i}   =  2 e_y d_i^T P0^{-1}          (phi_i applied to the x coordinate)
//     M_{2i+1} = -2 e_x d_i^T P0^{-1}          (phi_i applied to the y coordinate)
//
// equivalently  M_{2i} T = 2 e_y d_i^T  and  M_{2i+1} T = -2 e_x d_i^T.  The span of
// {M_j} is therefore the image of row(D) (+) row(D) under an injective linear map, so
// dim K = rank(A) = 2 rank(D).

#include <cstdlib>
#include <functional>
#include <random>
#include <set>
#include <map>
#include "kill_common.hpp"
#include "method/contact.hpp"
#include "method/periodic_jacobian.hpp"
#include "method/zero_plus.hpp"

using namespace kiri;
using namespace kiri::kill;
using namespace kiri::method;

namespace {

using Pattern = kiri::method::PeriodicPattern;

// Periodic Voronoi on the torus [0,L)^2: the cells of n random sites, computed
// against the 3x3 replicated site set, form a fundamental domain.
Pattern make_voronoi_pattern(int inst, int nsites, double L, std::mt19937& rng) {
  Pattern P;
  P.family = "voronoi_torus";
  P.name = "voronoi_torus_" + std::to_string(inst) + "_n" + std::to_string(nsites);
  std::uniform_real_distribution<double> U(0, L);
  std::vector<Vec2> site;
  int guard = 0;
  while (static_cast<int>(site.size()) < nsites && guard++ < 100000) {
    const Vec2 p(U(rng), U(rng));
    bool ok = true;
    for (const auto& s : site) {
      Vec2 d = p - s;
      d.x() -= L * std::round(d.x() / L);
      d.y() -= L * std::round(d.y() / L);
      if (d.norm() < 0.15 * L / std::sqrt(static_cast<double>(nsites))) ok = false;
    }
    if (ok) site.push_back(p);
  }
  if (static_cast<int>(site.size()) < nsites) {
    P.err = "site rejection failed";
    return P;
  }
  std::vector<Vec2> rep;
  for (int i = -1; i <= 1; ++i)
    for (int j = -1; j <= 1; ++j)
      for (const auto& s : site) rep.push_back(s + Vec2(i * L, j * L));
  std::vector<std::vector<Vec2>> polys;
  for (int i = 0; i < nsites; ++i) {
    const Vec2 pi = site[i];
    std::vector<Vec2> cellp{pi + Vec2(-L, -L), pi + Vec2(L, -L), pi + Vec2(L, L),
                            pi + Vec2(-L, L)};
    for (const auto& pj : rep) {
      if ((pj - pi).norm() < 1e-12) continue;
      const Vec2 nvec = pj - pi;
      const double off = nvec.dot(0.5 * (pi + pj));
      std::vector<Vec2> out;
      for (size_t k = 0; k < cellp.size(); ++k) {
        const Vec2& A = cellp[k];
        const Vec2& B = cellp[(k + 1) % cellp.size()];
        const double da = nvec.dot(A) - off, db = nvec.dot(B) - off;
        if (da <= 0) out.push_back(A);
        if ((da < 0 && db > 0) || (da > 0 && db < 0))
          out.push_back(A + (B - A) * (da / (da - db)));
      }
      cellp = out;
      if (cellp.size() < 3) break;
    }
    if (cellp.size() < 3) {
      P.err = "empty voronoi cell";
      return P;
    }
    polys.push_back(cellp);
  }
  Mesh cell;
  try {
    cell = mesh_from_polygons(polys, 1e-7);
  } catch (const std::exception& e) {
    P.err = std::string("weld: ") + e.what();
    return P;
  }
  Eigen::Matrix2d T;
  T << L, 0, 0, L;
  cell.sigma.assign(cell.n_faces(), -1);
  Quotient q0 = build_quotient(cell, T);
  if (!q0.ok) {
    P.err = "quotient(probe): " + q0.err;
    return P;
  }
  cell.sigma = quotient_sigma(cell.n_faces(), quotient_dual(q0), rng);
  P.cell = cell;
  P.T = T;
  P.ok = true;
  return P;
}

std::vector<Pattern> population() {
  std::vector<Pattern> out;
  const char* fams[] = {"squares",     "triangles",        "hexagons", "kagome",
                        "snub_square", "trunc_square_488", "t3_4_3_12"};
  const int sizes[][2] = {{2, 2}, {3, 2}, {3, 3}};
  for (const char* f : fams)
    for (const auto& s : sizes) {
      std::mt19937 rng(20260904u + 7919u * static_cast<unsigned>(out.size()));
      out.push_back(make_tiling_pattern(f, s[0], s[1], rng));
    }
  const int ns[] = {20, 28, 36, 45, 55, 70, 85, 100, 120, 140, 170, 200};
  for (int i = 0; i < 12; ++i) {
    std::mt19937 rng(9100001u + 104729u * static_cast<unsigned>(i));
    out.push_back(make_voronoi_pattern(i, ns[i], 10.0, rng));
  }
  return out;
}


// ------------------------------------------------------------------ per-pattern state

struct PState {
  bool ok = false;
  std::string err;
  Quotient q;
  SuperPatch sp;
  CutStructure cut;
  Eigen::MatrixXd X0, Phi;  // nq x 2, nq x k
  int k = 0;                // dim null
  double med_edge = 1.0;
  double cell_face_area = 0;  // sum of the |F| face areas of one cell (theta-independent)
};

std::vector<Vec2> shape_point(const PState& S, const Eigen::VectorXd& t) {
  Eigen::MatrixXd T(S.k, 2);
  for (int i = 0; i < S.k; ++i) {
    T(i, 0) = t(2 * i);
    T(i, 1) = t(2 * i + 1);
  }
  Eigen::MatrixXd X = S.X0;
  if (S.k) X += S.Phi * T;
  return matrix_to_points(X);
}

PState prepare(const Pattern& P) {
  PState S;
  S.q = build_quotient(P.cell, P.T);
  if (!S.q.ok) {
    S.err = "quotient: " + S.q.err;
    return S;
  }
  const LinearSystem sys = quotient_system(S.q, S.q.Xq);
  const SolveReport sr = solve_system(sys, S.q.Xq);
  if (!sr.projection_ok) {
    S.err = "projection failed";
    return S;
  }
  S.X0 = sr.X0;
  S.Phi = sr.Phi;
  S.k = sr.dim_null;
  S.sp = build_super(S.q, 1);
  S.cut = make_cut(S.sp.mesh);
  S.med_edge = median_edge_length(P.cell);
  S.ok = true;
  return S;
}

// The quotient positions are REDUCED into the origin cell, so face geometry must be
// read off the super patch (which restores the lattice offsets), never off `cell`.
double cell_area_sum(PState& S, const std::vector<Vec2>& Xq, int* inverted) {
  set_super_positions(S.sp, S.q, Xq);
  if (inverted) {
    int inv = 0;
    for (int f = 0; f < S.sp.nfc; ++f)
      if (S.sp.mesh.face_signed_area(S.sp.face_index(0, 0, f)) <= 0) ++inv;
    *inverted = inv;
  }
  return cell_face_area_sum(S.sp);
}

// K at a shape-space point, from the closed-form deploy basis of the super patch.
PeriodicJac K_at(PState& S, const std::vector<Vec2>& Xq) {
  set_super_positions(S.sp, S.q, Xq);
  const DeployBasis B = deploy_basis(S.cut, S.sp.mesh.X);
  return periodic_jacobian(S.sp, S.q, S.cut, B);
}

// P_theta measured by an INDEPENDENT forward-kinematics call (deploy(), not the basis).
bool P_fk(PState& S, double theta, Eigen::Matrix2d* P, double* spread) {
  *P = fk_period_matrix(S.sp, S.q, S.cut, theta, spread);
  return true;
}

// The affine map t -> K, plus the achievable-set rank and the period-potential matrix D.
AchievableSet achievable(PState& S) {
  AchievableSet A;
  const Eigen::VectorXd z = Eigen::VectorXd::Zero(2 * S.k);
  A.K0 = K_at(S, shape_point(S, z)).K;
  A.A = Eigen::MatrixXd::Zero(4, 2 * S.k);
  A.D = Eigen::MatrixXd::Zero(std::max(S.k, 1), 2);
  const double sc = S.med_edge;
  for (int j = 0; j < 2 * S.k; ++j) {
    Eigen::VectorXd e = Eigen::VectorXd::Zero(2 * S.k);
    e(j) = sc;
    const Eigen::Matrix2d M = (K_at(S, shape_point(S, e)).K - A.K0) / sc;
    A.M.push_back(M);
    A.A.col(j) << M(0, 0), M(1, 0), M(0, 1), M(1, 1);
    if (j % 2 == 0) {  // the (phi_i, e_1) generator carries d_i in row 2 of M*T
      const Eigen::Matrix2d MT = M * S.q.T;
      A.D(j / 2, 0) = 0.5 * MT(1, 0);
      A.D(j / 2, 1) = 0.5 * MT(1, 1);
    }
  }
  if (2 * S.k > 0) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A.A);
    const double tol = 1e-9 * (svd.singularValues().size() ? svd.singularValues()(0) : 1.0);
    for (int i = 0; i < svd.singularValues().size(); ++i)
      if (svd.singularValues()(i) > std::max(tol, 1e-13)) ++A.dimK;
    Eigen::JacobiSVD<Eigen::MatrixXd> sd(A.D);
    const double t2 = 1e-9 * (sd.singularValues().size() ? sd.singularValues()(0) : 1.0);
    for (int i = 0; i < sd.singularValues().size(); ++i)
      if (sd.singularValues()(i) > std::max(t2, 1e-13)) ++A.rankD;
  }
  return A;
}

Eigen::Matrix2d K_of_t(const AchievableSet& A, const Eigen::VectorXd& t) {
  Eigen::Matrix2d K = A.K0;
  for (int j = 0; j < static_cast<int>(A.M.size()); ++j) K += t(j) * A.M[j];
  return K;
}
// ---------------------------------------------------------------------------
static int rank_svd(const Eigen::MatrixXd& Mx, double rel = 1e-10) {
  if (Mx.size() == 0) return 0;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Mx);
  const auto& sv = svd.singularValues();
  if (sv.size() == 0) return 0;
  const double tol = std::max(rel * sv(0), 1e-13);
  int r = 0;
  for (int i = 0; i < sv.size(); ++i)
    if (sv(i) > tol) ++r;
  return r;
}

struct Row {
  std::string name;
  int k = 0, dimK = 0, rankD = 0, rankD_svd = 0;
  double err_fact = 0;    // max |M_j - (formula from row j/2 of D)|
  double err_xrow = 0;    // max |row 0 of (M_2i T)|  and |row 1 of (M_2i+1 T)|
  int rankD_cov = -1;
  bool have_cov = false;
  double consistency = 0, p0_err = 0;   // R2 (D-L2-5): H2 on the SUPER PATCH
  bool eq = false;
};

}  // namespace

int main(int argc, char** argv) {
  const int NEXTRA = (argc > 1) ? std::atoi(argv[1]) : 100;
  std::vector<Row> rows;
  long long n_eq = 0, n_fail = 0, n_cov_ok = 0, n_cov_bad = 0, n_inconsistent = 0;
  double worst_consistency = 0;
  double worst_fact = 0, worst_xrow = 0;
  std::map<int, int> rankD_hist;
  std::vector<std::string> failures;

  std::vector<Pattern> pats = population();          // the 33 K7 patterns
  // ... plus NEXTRA random Voronoi tori, fresh sigma seeds
  for (int i = 0; i < NEXTRA; ++i) {
    const int ns[] = {12, 16, 20, 24, 30, 36, 45};
    std::mt19937 rng(7000000u + 65537u * static_cast<unsigned>(i));
    Pattern P = make_voronoi_pattern(1000 + i, ns[i % 7], 10.0, rng);
    if (P.ok) pats.push_back(P);
  }

  for (auto& P : pats) {
    if (!P.ok) continue;
    PState S = prepare(P);
    if (!S.ok) continue;
    AchievableSet A = achievable(S);
    // R2 (D-L2-5): the sharp test of H2 on the patch the period is MEASURED on is
    // PeriodicJac::consistency -- the spread of the measured period over (face, corner).
    // achievable() never reads it, so Check L2 must.
    const PeriodicJac PJ0 = K_at(S, matrix_to_points(S.X0));
    Row R;
    R.consistency = PJ0.consistency;
    R.p0_err = PJ0.p0_err;
    R.name = P.name;
    R.k = S.k;
    R.dimK = A.dimK;
    R.rankD = A.rankD;
    R.rankD_svd = rank_svd(A.D);
    const int rA = rank_svd(A.A);
    // the factorisation
    const Eigen::Matrix2d T = S.q.T;
    for (int i = 0; i < S.k; ++i) {
      const Eigen::Vector2d d = A.D.row(i).transpose();
      Eigen::Matrix2d Mx = Eigen::Matrix2d::Zero(), My = Eigen::Matrix2d::Zero();
      Mx.row(1) = 2.0 * d.transpose();            // M_{2i}  T   = 2 e_y d^T
      My.row(0) = -2.0 * d.transpose();           // M_{2i+1} T = -2 e_x d^T
      const Eigen::Matrix2d Ti = T.inverse();
      const Eigen::Matrix2d Mx_pred = Mx * Ti, My_pred = My * Ti;
      R.err_fact = std::max(R.err_fact, (A.M[2 * i] - Mx_pred).cwiseAbs().maxCoeff());
      R.err_fact = std::max(R.err_fact, (A.M[2 * i + 1] - My_pred).cwiseAbs().maxCoeff());
      const Eigen::Matrix2d MxT = A.M[2 * i] * T, MyT = A.M[2 * i + 1] * T;
      R.err_xrow = std::max(R.err_xrow, MxT.row(0).cwiseAbs().maxCoeff());
      R.err_xrow = std::max(R.err_xrow, MyT.row(1).cwiseAbs().maxCoeff());
    }
    // L2.2: rank(D) = rank([L; delta_h; delta_v]) - rank(L), where delta_tau is the
    // period-potential increment read as a covector on ALL of R^nq (not just the null
    // space).  Finite differences in the quotient x-coordinates; only for small cells.
    if (S.q.nq <= 60 && S.q.L.size() > 0) {
      const int nq = S.q.nq;
      Eigen::MatrixXd G(2, nq);
      const double hstep = 1e-4 * S.med_edge;
      for (int v = 0; v < nq; ++v) {
        std::vector<Vec2> Xv = matrix_to_points(S.X0);
        Xv[v].x() += hstep;
        const Eigen::Matrix2d Mv = (K_at(S, Xv).K - A.K0) / hstep;
        const Eigen::Matrix2d MvT = Mv * T;
        G(0, v) = 0.5 * MvT(1, 0);
        G(1, v) = 0.5 * MvT(1, 1);
      }
      Eigen::MatrixXd LG(S.q.L.rows() + 2, nq);
      LG << S.q.L, G;
      const int rL = rank_svd(S.q.L), rLG = rank_svd(LG);
      R.rankD_cov = rLG - rL;
      R.have_cov = true;
    }
    if (R.rankD_svd == 0 && S.k > 0) {
      // R2 (D-L2-3): tr and det do NOT separate K0 = +J from -J; print K0(1,0) too.
      std::printf("    [rank(D) = 0 with k = %d]  max|D| = %.3e   "
                  "tr K0 = %.6f  det K0 = %.6f  K0(1,0) = %+.6f\n",
                  S.k, A.D.cwiseAbs().maxCoeff(), A.K0.trace(), A.K0.determinant(),
                  A.K0(1, 0));
    }
    if (R.consistency > 1e-9) {
      ++n_inconsistent;
      std::printf("    [H2 FAILS on the super patch] %-24s consistency = %.3e  "
                  "p0_err = %.3e   -> EXCLUDED from the L2.1 verification\n",
                  R.name.c_str(), R.consistency, R.p0_err);
    }
    R.eq = (rA == 2 * R.rankD_svd);
    if (R.eq) ++n_eq; else { ++n_fail; failures.push_back(R.name); }
    if (R.consistency <= 1e-9) worst_consistency = std::max(worst_consistency, R.consistency);
    worst_fact = std::max(worst_fact, R.err_fact);
    worst_xrow = std::max(worst_xrow, R.err_xrow);
    ++rankD_hist[R.rankD_svd];
    rows.push_back(R);
    std::printf("%-28s k=%3d  dimK=%d rank(A)=%d rank(D)=%d  2rD=%d %s  "
                "fact_err=%.2e  xrow_err=%.2e\n",
                R.name.c_str(), R.k, R.dimK, rA, R.rankD_svd, 2 * R.rankD_svd,
                R.eq ? "OK " : "FAIL", R.err_fact, R.err_xrow);
    if (R.have_cov) {
      if (R.rankD_cov == R.rankD_svd) ++n_cov_ok; else { ++n_cov_bad;
        std::printf("    [covector rank mismatch] %s: rank([L;delta])-rank(L) = %d, "
                    "rank(D) = %d\n", R.name.c_str(), R.rankD_cov, R.rankD_svd); }
    }
  }

  std::printf("\ncheck_l2 -- patterns evaluated: %zu\n", rows.size());
  std::printf("  rank(A) == 2 rank(D):  %lld / %zu    failures: %lld\n", n_eq, rows.size(),
              n_fail);
  for (auto& f : failures) std::printf("    FAIL: %s\n", f.c_str());
  std::printf("  rank(D) distribution:");
  for (auto& kv : rankD_hist) std::printf("   %d -> %d", kv.first, kv.second);
  std::printf("\n  max |M_j - formula(D)|            : %.3e\n", worst_fact);
  std::printf("  max |vanishing row of M_j T|      : %.3e\n", worst_xrow);
  std::printf("  H2 on the super patch (PeriodicJac::consistency > 1e-9): %lld of %zu patterns"
              "  [EXCLUDED]\n", n_inconsistent, rows.size());
  std::printf("  worst consistency among the retained patterns: %.3e\n", worst_consistency);
  std::printf("  L2.2 covector test  rank([L;delta_h;delta_v]) - rank(L) == rank(D):"
              "  %lld ok, %lld mismatch\n", n_cov_ok, n_cov_bad);
  std::printf("\n");
  return 0;
}
