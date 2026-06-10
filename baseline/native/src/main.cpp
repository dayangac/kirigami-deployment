// tuttekiri_cli - native driver for the AUTHORS' implementation.
//
// Every numerical step below is a call into baseline/tuttekiri/code/cpp; this
// file only decodes our JSON, mirrors the call sequence bind.cpp exposes to the
// web UI, and encodes the results back into our JSON conventions.
#include "bind_port.h"
#include "coloring/init.h"
#include "geometry/deployment.h"
#include "geometry/kirigami.h"
#include "geometry/unit_pattern.h"
#include "json_io.h"
#include "opt/fully_close.h"
#include "opt/prevent_intersections.h"
#include <filesystem>
#include "their_system.h"

#include <chrono>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <map>
#include <string>
#include <unistd.h>
#include <vector>

using json = nlohmann::json;
using Clock = std::chrono::steady_clock;

// ---------------------------------------------------------------- utilities

// The authors' code prints progress on stdout ("Solving system, ...",
// "No twin for edge ..."), which would corrupt our JSON.  Capture it.
class StdoutCapture {
public:
  StdoutCapture() {
    std::cout.flush();
    fflush(stdout);
    saved_ = dup(1);
    tmp_ = tmpfile();
    dup2(fileno(tmp_), 1);
  }
  std::string stop() {
    if (saved_ < 0)
      return text_;
    std::cout.flush();
    fflush(stdout);
    dup2(saved_, 1);
    close(saved_);
    saved_ = -1;
    fseek(tmp_, 0, SEEK_SET);
    char buf[4096];
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), tmp_)) > 0)
      text_.append(buf, n);
    fclose(tmp_);
    return text_;
  }
  ~StdoutCapture() { stop(); }

private:
  int saved_ = -1;
  FILE *tmp_ = nullptr;
  std::string text_;
};

struct Args {
  std::string cmd, in, out;
  std::map<std::string, std::string> opt;
  bool has(const std::string &k) const { return opt.count(k) > 0; }
  double num(const std::string &k, double d) const {
    auto it = opt.find(k);
    return it == opt.end() ? d : std::stod(it->second);
  }
};

static Args parse(int argc, char **argv) {
  Args a;
  a.cmd = argv[1];
  std::vector<std::string> pos;
  for (int i = 2; i < argc; i++) {
    std::string s = argv[i];
    if (s.rfind("--", 0) == 0) {
      std::string k = s.substr(2);
      if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0)
        a.opt[k] = argv[++i];
      else
        a.opt[k] = "1";
    } else
      pos.push_back(s);
  }
  if (!pos.empty())
    a.in = pos[0];
  if (a.opt.count("out"))
    a.out = a.opt["out"];
  return a;
}

// Edge classification in the authors' own terms: interior edge, hinge if the
// two incident face colours differ, split if they are equal.
struct Cuts {
  std::vector<std::pair<int, int>> hinge, split; // (v_from, v_to), v_from<v_to
  int n_interior = 0;
};
static Cuts classify(const kirigami::UnitPattern &p) {
  Cuts c;
  auto &mesh = p.mesh;
  for (auto &e : mesh.edges) {
    if (!e.twin() || e.index > e.twin()->index)
      continue;
    c.n_interior++;
    int a = e.vi, b = e.next()->vi;
    if (a > b)
      std::swap(a, b);
    if (p.face_colors[e.fi] != p.face_colors[e.twin()->fi])
      c.hinge.push_back({a, b});
    else
      c.split.push_back({a, b});
  }
  std::sort(c.hinge.begin(), c.hinge.end());
  std::sort(c.split.begin(), c.split.end());
  return c;
}

static int count_interior_vertices(const kirigami::UnitPattern &p) {
  auto &mesh = const_cast<utils::Hmesh &>(p.mesh);
  int n = 0;
  for (int i = 0; i < mesh.nv(); i++)
    if (!mesh.is_boundary_vertex(i))
      n++;
  return n;
}

static json edges_to_json(const std::vector<std::pair<int, int>> &e) {
  json a = json::array();
  for (auto &p : e)
    a.push_back({p.first, p.second});
  return a;
}

// Their make_deployable, with the printed diagnostics recovered.
struct SolveResult {
  kirigami::UnitPattern X0;
  Eigen::MatrixXd kernel; // nv x k, empty if trivial
  int printed_holes = -1, printed_constraints = -1;
  double seconds = 0;
};
static SolveResult run_solve(const kirigami::UnitPattern &p) {
  SolveResult r;
  auto t0 = Clock::now();
  std::string log;
  {
    StdoutCapture cap;
    auto [res, kernel] = kirigami::make_deployable(p, false);
    r.X0 = res;
    r.kernel = kernel;
    log = cap.stop();
  }
  r.seconds = std::chrono::duration<double>(Clock::now() - t0).count();
  auto pos = log.find("num_holes:");
  if (pos != std::string::npos)
    r.printed_holes = std::stoi(log.substr(pos + 10));
  pos = log.find("num_constraints:");
  if (pos != std::string::npos)
    r.printed_constraints = std::stoi(log.substr(pos + 16));
  if (r.kernel.cols() > 0 && r.kernel.norm() < 1e-4)
    r.kernel = Eigen::MatrixXd(p.mesh.nv(), 0); // Eigen's trivial-kernel column
  return r;
}

static double theta_max(const kirigami::UnitPattern &p, bool collisions) {
  StdoutCapture cap;
  double v = max_opening_angle_with_collisions(p, collisions);
  cap.stop();
  return v;
}

// UnitPattern::get_holes() dereferences cur_edge->prev()->twin() without a null
// check (unit_pattern.cpp, the face_colors == 1 branch), so it crashes on any
// mesh that still has boundary half-edges after make_periodic().  It is only
// safe on a fully periodic unit pattern.  Returns -1 when it cannot be run.
static int n_traced_holes(const kirigami::UnitPattern &p) {
  StdoutCapture cap;
  auto q = p;
  q.make_periodic();
  bool has_boundary = false;
  for (auto &e : q.mesh.edges)
    if (!e.twin())
      has_boundary = true;
  int n = has_boundary ? -1 : (int)q.get_holes().size();
  cap.stop();
  return n;
}

// ------------------------------------------------------------------ commands

static json cmd_color(const io::Input &in) {
  auto p = in.pattern;
  auto t0 = Clock::now();
  std::vector<int> colors;
  {
    StdoutCapture cap;
    colors = coloring::initialized_two_face_coloring(p.mesh);
    cap.stop();
  }
  double secs = std::chrono::duration<double>(Clock::now() - t0).count();
  p.face_colors = colors;
  json j = io::pattern_to_json(p);
  auto c = classify(p);
  j["n_hinge"] = c.hinge.size();
  j["n_split"] = c.split.size();
  j["seconds"] = secs;
  return j;
}

static json cmd_deploy(const kirigami::UnitPattern &p, double theta) {
  auto t0 = Clock::now();
  StdoutCapture cap;
  auto [m, per] = kirigami::deploy(p.mesh, p.face_colors, p.periodicity, theta);
  cap.stop();
  json j;
  j["theta"] = theta;
  j["vertices"] = io::verts_to_json(m.V);
  j["faces"] = io::faces_to_json(m.F);
  j["seconds"] = std::chrono::duration<double>(Clock::now() - t0).count();
  return j;
}

// Full pipeline, the parity workhorse.
static json cmd_analyze(const io::Input &in, const std::string &outdir,
                        bool collisions, double deploy_theta) {
  auto p = in.pattern;
  json j;
  j["nv"] = p.mesh.nv();
  j["nf"] = p.mesh.nf();
  j["n_cw_faces_reversed_on_load"] = in.n_cw_faces_reversed;
  j["had_orientation"] = in.had_orientation;
  j["had_periodic"] = in.had_periodic;
  j["n_interior_vertices"] = count_interior_vertices(p);

  auto c = classify(p);
  j["n_interior_edges"] = c.n_interior;
  j["n_hinge"] = c.hinge.size();
  j["n_split"] = c.split.size();
  j["E_hinge"] = edges_to_json(c.hinge);
  j["E_split"] = edges_to_json(c.split);

  j["n_traced_holes"] = n_traced_holes(p);

  auto sol = run_solve(p);
  j["solve_seconds"] = sol.seconds;
  j["their_num_holes_rows"] = sol.printed_holes;
  j["their_num_constraints"] = sol.printed_constraints;
  int k = sol.kernel.cols();
  j["dim_null_per_coordinate"] = k;
  j["dim_null_2d"] = 2 * k;
  j["rank_A"] = p.mesh.nv() - k;

  // Residuals in THEIR system.  The rebuild is validated first.
  auto pp = p;
  if (pp.periodicity.norm() > 1e-6) {
    StdoutCapture cap;
    pp.make_periodic();
    cap.stop();
  }
  auto sys = their::build(pp);
  j["rebuilt_num_holes_rows"] = sys.num_holes;
  j["rebuilt_num_constraints"] = sys.num_constraints;
  bool counts_ok = (sol.printed_holes < 0 || sol.printed_holes == sys.num_holes) &&
                   (sol.printed_constraints < 0 ||
                    sol.printed_constraints == sys.num_constraints);
  double res_after = (sys.A * sol.X0.mesh.V - sys.b).norm();
  bool sol_ok = res_after < 1e-6 * std::max(1.0, sys.b.norm());
  j["rebuild_validated"] = counts_ok && sol_ok;
  j["residual_before"] = (sys.A * p.mesh.V - sys.b).norm();
  j["residual_after"] = res_after;

  j["theta_max_kinematic"] = theta_max(p, false);
  if (collisions)
    j["theta_max_with_collisions"] = theta_max(p, true);

  if (!outdir.empty()) {
    json x0 = io::pattern_to_json(sol.X0);
    if (k > 0) {
      json basis = json::array();
      for (int i = 0; i < k; i++) {
        json col = json::array();
        for (int v = 0; v < sol.kernel.rows(); v++)
          col.push_back(sol.kernel(v, i));
        basis.push_back(col);
      }
      x0["kernel"] = basis;
    }
    io::write(outdir + "/X0.json", x0);
    io::write(outdir + "/deploy.json", cmd_deploy(p, deploy_theta));
    auto x0p = sol.X0;
    io::write(outdir + "/deploy_X0.json", cmd_deploy(x0p, deploy_theta));
    j["theta_max_X0_kinematic"] = theta_max(x0p, false);
    if (collisions)
      j["theta_max_X0_with_collisions"] = theta_max(x0p, true);
  }
  return j;
}

static json run_prevent(const kirigami::UnitPattern &p,
                        const Eigen::MatrixXd &kernel, double barrier,
                        double strength, double close_w, bool collisions,
                        const std::string &dump = "") {
  auto q = p;
  auto t0 = Clock::now();
  {
    StdoutCapture cap;
    q.mesh.V = opt::prevent_intersections(p, kernel, barrier, strength, close_w);
    cap.stop();
  }
  json j;
  j["barrier"] = barrier;
  j["barrier_strength"] = strength;
  j["close_to_original_weight"] = close_w;
  j["seconds"] = std::chrono::duration<double>(Clock::now() - t0).count();
  j["theta_max_kinematic"] = theta_max(q, false);
  if (collisions)
    j["theta_max_with_collisions"] = theta_max(q, true);
  j["finite"] = q.mesh.V.allFinite();
  if (!dump.empty())
    io::write(dump, io::pattern_to_json(q));
  return j;
}

// ---- parity helpers --------------------------------------------------------

// Rigid (rotation + translation) Procrustes alignment of P onto Q.
struct Fit {
  double rms = 0, max_dev = 0;
  bool reflected = false;
};
static Fit rigid_fit(const Eigen::MatrixXd &P, const Eigen::MatrixXd &Q) {
  Eigen::RowVector2d cp = P.colwise().mean(), cq = Q.colwise().mean();
  Eigen::MatrixXd Pc = P.rowwise() - cp, Qc = Q.rowwise() - cq;
  Eigen::Matrix2d H = Pc.transpose() * Qc;
  Eigen::JacobiSVD<Eigen::Matrix2d> svd(H, Eigen::ComputeFullU |
                                               Eigen::ComputeFullV);
  Eigen::Matrix2d R = svd.matrixV() * svd.matrixU().transpose();
  Fit f;
  if (R.determinant() < 0) {
    f.reflected = true;
    Eigen::Matrix2d D = Eigen::Matrix2d::Identity();
    D(1, 1) = -1;
    R = svd.matrixV() * D * svd.matrixU().transpose();
  }
  Eigen::MatrixXd d = Pc * R.transpose() - Qc;
  f.rms = std::sqrt(d.squaredNorm() / d.rows());
  f.max_dev = d.rowwise().norm().maxCoeff();
  return f;
}

// Face-corner point cloud in (face, corner) order: the one ordering both
// implementations agree on, since M' has one vertex per (face, corner).
static Eigen::MatrixXd corner_cloud(const nlohmann::json &j) {
  const auto &V = j.at("vertices");
  const auto &F = j.at("faces");
  int n = 0;
  for (const auto &f : F)
    n += f.size();
  Eigen::MatrixXd P(n, 2);
  int r = 0;
  for (const auto &f : F)
    for (const auto &vi : f) {
      int i = vi.get<int>();
      P(r, 0) = V[i][0].get<double>();
      P(r, 1) = V[i][1].get<double>();
      r++;
    }
  return P;
}

static json cmd_parity(const std::string &dir, double theta, bool collisions,
                       const std::string &ours_cut) {
  io::Input in = io::load(dir + "/M.json");
  json j = cmd_analyze(in, "", collisions, theta);

  // sigma: their coloring heuristic vs. the sigma stored in M.json.
  {
    std::vector<int> theirs;
    {
      StdoutCapture cap;
      theirs = coloring::initialized_two_face_coloring(in.pattern.mesh);
      cap.stop();
    }
    auto q = in.pattern;
    q.face_colors = theirs;
    auto ct = classify(q);
    int same = 0;
    for (size_t i = 0; i < theirs.size(); i++)
      same += (theirs[i] == in.pattern.face_colors[i]);
    int agree = std::max(same, (int)theirs.size() - same); // allow global flip
    j["sigma_agree_up_to_flip"] = agree;
    j["sigma_total"] = (int)theirs.size();
    j["their_coloring_n_hinge"] = (int)ct.hinge.size();
    j["their_coloring_n_split"] = (int)ct.split.size();
  }

  // E_hinge as a SET, against our kiri_analyze cut.json when it is supplied.
  if (!ours_cut.empty()) {
    std::ifstream f(ours_cut);
    if (f) {
      json cj; f >> cj;
      std::vector<std::pair<int, int>> ours_hinge;
      for (const auto &h : cj.at("hinge")) {
        int a = h.at("src").get<int>(), b = h.at("dst").get<int>();
        if (a > b) std::swap(a, b);
        ours_hinge.push_back({a, b});
      }
      std::sort(ours_hinge.begin(), ours_hinge.end());
      auto ct = classify(in.pattern);
      std::vector<std::pair<int, int>> diff;
      std::set_symmetric_difference(ct.hinge.begin(), ct.hinge.end(),
                                    ours_hinge.begin(), ours_hinge.end(),
                                    std::back_inserter(diff));
      j["ours_n_hinge"] = (int)ours_hinge.size();
      j["E_hinge_symdiff"] = (int)diff.size();
      j["E_hinge_sets_equal"] = diff.empty();
    }
  }

  // Forward kinematics at theta, against our dump.
  std::string ours_path =
      dir + "/deploy_" + std::to_string((int)std::lround(theta * 180 / M_PI)) +
      "deg.json";
  std::ifstream test(ours_path);
  if (test) {
    json ours;
    { std::ifstream f(ours_path); f >> ours; }
    json theirs = cmd_deploy(in.pattern, theta);
    Eigen::MatrixXd P = corner_cloud(theirs), Q = corner_cloud(ours);
    if (P.rows() != Q.rows()) {
      j["fk_error"] = "corner count mismatch";
      j["fk_corners_theirs"] = (int)P.rows();
      j["fk_corners_ours"] = (int)Q.rows();
    } else {
      double scale = Q.rowwise().norm().maxCoeff();
      Fit f = rigid_fit(P, Q);
      j["fk_corners"] = (int)P.rows();
      j["fk_rms"] = f.rms;
      j["fk_max_dev"] = f.max_dev;
      j["fk_rel_max_dev"] = f.max_dev / std::max(scale, 1e-12);
      j["fk_reflected"] = f.reflected;
    }
  } else {
    j["fk_error"] = "no " + ours_path;
  }
  return j;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr
        << "usage: tuttekiri_cli <color|solve|deploy|thetamax|analyze|"
           "prevent|closed> <in.json> [--out f] [--theta t] [--collisions]\n";
    return 2;
  }
  Args a = parse(argc, argv);
  try {
    if (a.in.empty() && a.cmd != "compare")
      throw std::runtime_error("no input file");
    io::Input in = (a.cmd == "parity" || a.cmd == "compare")
                       ? io::Input{}
                       : io::load(a.in);
    bool collisions = a.has("collisions");
    double theta = a.num("theta", M_PI / 6);
    json out;

    if (a.cmd == "color") {
      out = cmd_color(in);
    } else if (a.cmd == "deploy") {
      out = cmd_deploy(in.pattern, theta);
    } else if (a.cmd == "thetamax") {
      out["theta_max_kinematic"] = theta_max(in.pattern, false);
      out["theta_max_with_collisions"] = theta_max(in.pattern, true);
    } else if (a.cmd == "solve") {
      auto sol = run_solve(in.pattern);
      out = io::pattern_to_json(sol.X0);
      out["dim_null_per_coordinate"] = sol.kernel.cols();
      out["seconds"] = sol.seconds;
      json basis = json::array();
      for (int i = 0; i < sol.kernel.cols(); i++) {
        json col = json::array();
        for (int v = 0; v < sol.kernel.rows(); v++)
          col.push_back(sol.kernel(v, i));
        basis.push_back(col);
      }
      out["kernel"] = basis;
    } else if (a.cmd == "analyze") {
      out = cmd_analyze(in, a.opt.count("dump") ? a.opt["dump"] : "",
                        collisions, theta);
    } else if (a.cmd == "prevent") {
      auto sol = run_solve(in.pattern);
      auto base = sol.X0;
      out["dim_null_per_coordinate"] = sol.kernel.cols();
      out["theta_max_input_kinematic"] = theta_max(in.pattern, false);
      out["theta_max_X0_kinematic"] = theta_max(base, false);
      if (collisions) {
        out["theta_max_input_with_collisions"] = theta_max(in.pattern, true);
        out["theta_max_X0_with_collisions"] = theta_max(base, true);
      }
      if (sol.kernel.cols() == 0) {
        out["note"] = "trivial shape space: Eq. (9) has nothing to optimise";
      } else if (a.has("sweep")) {
        // --dumpdir d  writes each ladder point's optimised pattern to d/run_<i>.json,
        // so an external referee (our collision bisection) can score every point rather
        // than trusting the authors' own theta_max, which merge_close_verts corrupts
        // (baseline/README.md upstream defect 2).
        // --short     runs the 36-point sub-ladder {0,0.1,0.3,1.0} x {1,10,100} x
        //             {0,0.1,1.0}, which contains the published default and every corner
        //             of the full grid; used where the full 84 points are too slow.
        const std::string dumpdir =
            a.opt.count("dumpdir") ? a.opt["dumpdir"] : std::string();
        if (!dumpdir.empty())
          std::filesystem::create_directories(dumpdir);
        const std::vector<double> bs =
            a.has("short") ? std::vector<double>{0.0, 0.1, 0.3, 1.0}
                           : std::vector<double>{0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0};
        const std::vector<double> ss = a.has("short")
                                           ? std::vector<double>{1.0, 10.0, 100.0}
                                           : std::vector<double>{1.0, 10.0, 50.0, 100.0};
        const std::vector<double> ws{0.0, 0.1, 1.0};
        json runs = json::array();
        int idx = 0;
        for (double b : bs)
          for (double s : ss)
            for (double w : ws) {
              const std::string dmp =
                  dumpdir.empty() ? std::string()
                                  : dumpdir + "/run_" + std::to_string(idx) + ".json";
              json r = run_prevent(base, sol.kernel, b, s, w, collisions, dmp);
              r["index"] = idx;
              if (!dmp.empty()) r["dump"] = dmp;
              runs.push_back(r);
              ++idx;
            }
        out["sweep"] = runs;
        if (!dumpdir.empty()) {
          // X0 itself is the T = 0 fallback of the ladder; dump it too.
          io::write(dumpdir + "/run_X0.json", io::pattern_to_json(base));
          out["dump_X0"] = dumpdir + "/run_X0.json";
        }
      } else {
        out["default"] = run_prevent(
            base, sol.kernel, a.num("barrier", 0.1), a.num("strength", 10.0),
            a.num("close", 0.1), collisions,
            a.opt.count("dump") ? a.opt["dump"] : "");
      }
    } else if (a.cmd == "closed") {
      auto sol = run_solve(in.pattern);
      auto q = sol.X0;
      auto t0 = Clock::now();
      {
        StdoutCapture cap;
        q.mesh.V = opt::optimize_for_fully_closed(sol.X0, sol.kernel);
        cap.stop();
      }
      out = io::pattern_to_json(q);
      out["seconds"] = std::chrono::duration<double>(Clock::now() - t0).count();
      out["dim_null_per_coordinate"] = sol.kernel.cols();
      out["theta_max_before"] = theta_max(sol.X0, false);
      out["theta_max_after"] = theta_max(q, false);
      out["finite"] = q.mesh.V.allFinite();
    } else if (a.cmd == "collide") {
      // Bisect their own collision predicate to locate first contact.
      double hi = a.num("hi", theta_max(in.pattern, false));
      out["theta_max_kinematic"] = hi;
      json scan = json::array();
      double first = hi;
      for (int i = 1; i <= 200; i++) {
        double t = hi * i / 200.0;
        int f0 = -1, f1 = -1;
        bool hit;
        { StdoutCapture cap; hit = first_collision_at(in.pattern, t, &f0, &f1); cap.stop(); }
        if (hit) {
          first = hi * (i - 1) / 200.0;
          out["first_contact"] = first;
          out["first_contact_faces"] = {f0, f1};
          break;
        }
      }
      out["first_contact_or_kinematic"] = out.contains("first_contact")
                                              ? out["first_contact"]
                                              : json(hi);
      // Merge tolerance their scan applies before testing.
      out["merge_tolerance"] = 0.1 * in.pattern.mesh.avg_edge_len();
      // Full occupancy scan: their predicate can fire at small angles because
      // merge_close_verts() fuses hinge duplicates that are still closer than
      // 10% of the average edge length, then release again.
      json occ = json::array();
      for (int i = 1; i <= 60; i++) {
        double t = hi * i / 60.0;
        int f0 = -1, f1 = -1;
        bool hit;
        { StdoutCapture cap; hit = first_collision_at(in.pattern, t, &f0, &f1); cap.stop(); }
        occ.push_back({t, hit, f0, f1});
      }
      out["scan"] = occ;
    } else if (a.cmd == "compare") {
      json A, B;
      { std::ifstream f(a.opt["a"]); f >> A; }
      { std::ifstream f(a.opt["b"]); f >> B; }
      Eigen::MatrixXd P = corner_cloud(A), Q = corner_cloud(B);
      Fit f = rigid_fit(P, Q);
      out["corners"] = (int)P.rows();
      out["rms"] = f.rms;
      out["max_dev"] = f.max_dev;
      out["reflected"] = f.reflected;
    } else if (a.cmd == "parity") {
      out = cmd_parity(a.in, theta, collisions,
                       a.opt.count("ourscut") ? a.opt["ourscut"] : "");
    } else {
      throw std::runtime_error("unknown command: " + a.cmd);
    }

    if (!a.out.empty())
      io::write(a.out, out);
    else
      std::cout << out.dump(2) << "\n";
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
