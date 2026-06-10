// Standalone check of two claims made by the rigidity ideator:
//  V1: Y_theta = cos(theta/2)*C + sin(theta/2)*S exactly (forward kinematics is
//      linear in (cos(theta/2), sin(theta/2))).
//  V2: mobility m = dim ker A - 1, where A(omega)_z = sum_e z_e (omega_{f1}-omega_{f2}) p_e
//      over a cycle basis {z} of the hinge graph Gamma, compared with
//      m = 3|F| - rank(RigidityMatrix) - 3.
//  V3: sigma in ker A  <=>  Eq.(2) holds (uniform deployability).
//  V4: dim ker A(theta) along the deployment path.
#include <cstdio>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include "core/mesh.hpp"
#include "core/cut.hpp"
#include "core/holes.hpp"
#include "core/kinematics.hpp"
#include "core/tutte_auxetic.hpp"
#include "core/generators.hpp"
#include "core/orientation.hpp"

using namespace kiri;

static Eigen::Matrix2d Jrot() { Eigen::Matrix2d J; J << 0,-1, 1,0; return J; }

struct Gamma {
  int F = 0;
  std::vector<int> h, t;              // per hinge edge: head face, tail face
  std::vector<int> eid;               // mesh edge id
  std::vector<Eigen::VectorXd> cyc;   // cycle basis vectors in R^{|E_h|}
  int comps = 0;
};

static Gamma build_gamma(const CutStructure& c) {
  const Mesh& m = *c.mesh;
  Gamma g; g.F = m.n_faces();
  for (int e : c.hinge_edges) {
    const Edge& ed = m.edges[e];
    g.h.push_back(m.half_edges[ed.he[0]].face);
    g.t.push_back(m.half_edges[ed.he[1]].face);
    g.eid.push_back(e);
  }
  const int E = (int)g.h.size();
  std::vector<std::vector<std::pair<int,int>>> adj(g.F); // (edge, other face)
  for (int i = 0; i < E; ++i) { adj[g.h[i]].push_back({i,g.t[i]}); adj[g.t[i]].push_back({i,g.h[i]}); }
  std::vector<int> parent(g.F,-2), pedge(g.F,-1);
  std::vector<Eigen::VectorXd> gpath(g.F, Eigen::VectorXd::Zero(E));
  std::vector<char> intree(E,0);
  for (int s = 0; s < g.F; ++s) {
    if (parent[s] != -2) continue;
    g.comps++;
    parent[s] = -1; std::vector<int> stack{s};
    while (!stack.empty()) {
      int f = stack.back(); stack.pop_back();
      for (auto [i,o] : adj[f]) {
        if (parent[o] != -2) continue;
        parent[o] = f; pedge[o] = i; intree[i] = 1;
        gpath[o] = gpath[f];
        gpath[o][i] += (g.h[i] == o ? +1.0 : -1.0);   // traversing f -> o
        stack.push_back(o);
      }
    }
  }
  for (int i = 0; i < E; ++i) {
    if (intree[i]) continue;
    Eigen::VectorXd z = gpath[g.t[i]];
    z[i] += 1.0;
    z -= gpath[g.h[i]];
    g.cyc.push_back(z);
  }
  return g;
}

// A: (2 * n_cycles) x F
static Eigen::MatrixXd build_A(const Gamma& g, const std::vector<Vec2>& p) {
  const int nz = (int)g.cyc.size();
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*nz, g.F);
  for (int k = 0; k < nz; ++k)
    for (int i = 0; i < (int)g.h.size(); ++i) {
      double z = g.cyc[k][i];
      if (z == 0) continue;
      A(2*k+0, g.h[i]) += z * p[i].x();  A(2*k+0, g.t[i]) -= z * p[i].x();
      A(2*k+1, g.h[i]) += z * p[i].y();  A(2*k+1, g.t[i]) -= z * p[i].y();
    }
  return A;
}

// full body-pin rigidity matrix: 2|E_h| x 3|F|, unknowns (w_f, vx_f, vy_f)
static Eigen::MatrixXd build_R(const Gamma& g, const std::vector<Vec2>& p) {
  const int E = (int)g.h.size();
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(2*E, 3*g.F);
  Eigen::Matrix2d J = Jrot();
  for (int i = 0; i < E; ++i) {
    Vec2 Jp = J * p[i];
    R(2*i+0, 3*g.h[i]+0) += Jp.x();  R(2*i+0, 3*g.t[i]+0) -= Jp.x();
    R(2*i+1, 3*g.h[i]+0) += Jp.y();  R(2*i+1, 3*g.t[i]+0) -= Jp.y();
    R(2*i+0, 3*g.h[i]+1) += 1;       R(2*i+0, 3*g.t[i]+1) -= 1;
    R(2*i+1, 3*g.h[i]+2) += 1;       R(2*i+1, 3*g.t[i]+2) -= 1;
  }
  return R;
}

static int rank_of(const Eigen::MatrixXd& M, double tol = 1e-9) {
  if (M.size() == 0) return 0;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(M);
  auto s = svd.singularValues();
  double t = tol * (s.size() ? s(0) : 1.0);
  int r = 0; for (int i = 0; i < s.size(); ++i) if (s(i) > t) ++r;
  return r;
}

static void run1(const char* name, Mesh m, std::mt19937& rng);
static void run(const char* name, Mesh m, std::mt19937& rng) {
  try { run1(name, m, rng); } catch (const std::exception& e) { printf("%-20s EXCEPTION: %s\n", name, e.what()); }
}
static void run1(const char* name, Mesh m, std::mt19937& rng) {
  m.build_topology();
  if (m.sigma.empty()) {
    OrientationReport orep = (m.n_faces() <= 18) ? brute_force_orientation(m)
                                                 : assign_orientation_relaxation(m, rng);
    m.sigma = orep.sigma;
  }
  if (m.sigma.empty() || (int)m.sigma.size() != m.n_faces()) {
    printf("%-20s sigma assignment failed, skipped\n", name); return; }
  m.build_topology();
  CutStructure c = make_cut(m);
  Gamma g = build_gamma(c);
  if (g.h.empty()) { printf("%-22s no hinge edges, skipped\n", name); return; }

  // hinge points at theta = 0 are just X[src]
  std::vector<Vec2> p0(g.h.size());
  for (size_t i = 0; i < g.h.size(); ++i) p0[i] = m.X[c.hinge_dir[g.eid[i]].src];

  // ---- V1: linearity of forward kinematics in (cos(th/2), sin(th/2)) ----
  auto Yof = [&](double th){ return deploy(c, m.X, th, 0).Y; };
  const double ta = 0.37, tb = 1.13;
  auto Ya = Yof(ta), Yb = Yof(tb);
  double ca = std::cos(ta/2), sa = std::sin(ta/2), cb = std::cos(tb/2), sb = std::sin(tb/2);
  double det = ca*sb - sa*cb;
  const int NP = c.n_prime_vertices;
  std::vector<Vec2> C(NP), S(NP);
  for (int i = 0; i < NP; ++i) {
    C[i] = ( sb*Ya[i] - sa*Yb[i]) / det;
    S[i] = (-cb*Ya[i] + ca*Yb[i]) / det;
  }
  double v1err = 0, scale = 0;
  for (double th : {0.05, 0.7, 1.9, 2.5, 3.0}) {
    auto Y = Yof(th);
    double cc = std::cos(th/2), ss = std::sin(th/2);
    for (int i = 0; i < NP; ++i) {
      v1err = std::max(v1err, (Y[i] - (cc*C[i] + ss*S[i])).norm());
      scale = std::max(scale, Y[i].norm());
    }
  }

  // ---- V2: mobility ----
  Eigen::MatrixXd A0 = build_A(g, p0);
  Eigen::MatrixXd R0 = build_R(g, p0);
  int dimW = g.F - rank_of(A0);
  int m_A  = dimW - g.comps;                     // one trivial global rotation per component
  int m_R  = 3*g.F - rank_of(R0) - 3*g.comps;    // Maxwell/Calladine, trivial motions removed

  // ---- V3: sigma in ker A  <=>  Eq.(2) ----
  Eigen::VectorXd sg(g.F);
  for (int f = 0; f < g.F; ++f) sg[f] = m.sigma[f];
  double sigres = (A0*sg).norm();
  HoleSet hs = holes_partition(c);
  Residuals rr = hole_residuals(c, m.X, hs);

  // ---- V3b: repeat on a SOLVED (deployable) embedding X0 ----
  LinearSystem sys = assemble_system(c, hs, m.X, BoundaryMode::Fixed);
  SolveReport sr = solve_system(sys, m.X);
  std::vector<Vec2> X0 = matrix_to_points(sr.X0);
  std::vector<Vec2> p0s(g.h.size());
  for (size_t i = 0; i < g.h.size(); ++i) p0s[i] = X0[c.hinge_dir[g.eid[i]].src];
  Eigen::MatrixXd As = build_A(g, p0s);
  int dimWs = g.F - rank_of(As);
  int mAs = dimWs - g.comps;
  int mRs = 3*g.F - rank_of(build_R(g, p0s)) - 3*g.comps;
  double sigres_s = (As*sg).norm();
  Residuals rrs = hole_residuals(c, X0, hs);

  // ---- V4: dim ker A along the path ----
  int nchange = 0; int dim0 = -1; char dims[128] = {0}; int dp = 0;
  for (double th : {0.0, 0.2, 0.6, 1.0, 1.6, 2.2, 2.8}) {
    auto Y = deploy(c, X0, th, 0).Y;
    std::vector<Vec2> pth(g.h.size());
    for (size_t i = 0; i < g.h.size(); ++i)
      pth[i] = Y[c.prime_vertex(g.h[i], c.hinge_dir[g.eid[i]].src)];
    int d = g.F - rank_of(build_A(g, pth));
    dp += snprintf(dims+dp, sizeof(dims)-dp, "%d ", d - g.comps);
    if (dim0 < 0) dim0 = d; else if (d != dim0) ++nchange;
  }

  printf("%-20s |F|=%4d |Eh|=%4d |Es|=%3d H=%4d c=%d | V1=%.2e(sc %.1f)"
         " | ini: m_A=%3d m_R=%3d %s |As|=%.1e res=%.1e"
         " | sol: m_A=%3d m_R=%3d %s |As|=%.1e res=%.1e | rkL=%d/%d m(th)= %s\n",
         name, g.F, (int)g.h.size(), c.n_split(), hs.n_interior_holes(), g.comps,
         v1err, scale, m_A, m_R, (m_A==m_R?"OK":"BAD"), sigres, rr.max_norm,
         mAs, mRs, (mAs==mRs?"OK":"BAD"), sigres_s, rrs.max_norm,
         sr.rank_L, sys.n_hole_rows, dims);
}

int main() {
  std::mt19937 rng(12345);
  run("squares 3x3",        tiling_squares(rect({0,0},1.6,1.6)), rng);
  run("squares 5x5",        tiling_squares(rect({0,0},2.6,2.6)), rng);
  run("triangles",          tiling_triangles(disk({0,0},2.0)), rng);
  run("hexagons",           tiling_hexagons(disk({0,0},2.5)), rng);
  run("kagome",             tiling_kagome(disk({0,0},2.5)), rng);
  run("periodic squares",   periodic_squares(3,3), rng);
  for (int k = 0; k < 8; ++k) {
    std::mt19937 r2(100+k);
    run("delaunay random",  delaunay_of_random_points(30, 3.0, r2), rng);
  }
  for (int k = 0; k < 8; ++k) {
    std::mt19937 r2(200+k);
    run("voronoi random",   voronoi_of_random_points(25, 3.0, r2), rng);
  }
  return 0;
}
