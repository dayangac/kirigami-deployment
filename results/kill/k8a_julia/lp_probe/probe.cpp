// C++ side of the K8a LP-layer probe: same inputs (JSON), same quantities as probe.jl.
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include "core/cut.hpp"
#include "core/mesh.hpp"
#include "method/expansive_cone.hpp"
#include "method/mobility.hpp"
using namespace kiri; using namespace kiri::method; using json = nlohmann::json;
static std::vector<char> chart(const ConeSystem& cs, const std::vector<int>& idx, const Eigen::VectorXd& a) {
  const int n = static_cast<int>(idx.size());
  std::vector<char> keep(n, 1);
  for (int i = 0; i + 1 < n; ++i) {
    const ConeRow& r0 = cs.rows[idx[i]]; const ConeRow& r1 = cs.rows[idx[i + 1]];
    if (r0.kind != ConeRowKind::CornerG1 || !r0.convex || r1.kind != ConeRowKind::CornerG2 || r1.index != r0.index) continue;
    if (a(i) >= a(i + 1)) keep[i + 1] = 0; else keep[i] = 0;
  }
  return keep;
}
int main(int argc, char** argv) {
  json J = json::parse(std::ifstream(argv[1]));
  json out = json::array();
  for (const auto& r : J) {
    Mesh m = mesh_from_json_string(r["mesh"].dump());
    m.build_topology();
    std::vector<Vec2> X = m.X;
    if (r.contains("X")) { X.clear(); for (const auto& p : r["X"]) X.emplace_back(p[0].get<double>(), p[1].get<double>()); }
    CutStructure c = make_cut(m); c.mesh = &m;
    const HingeGraph g = build_hinge_graph(c);
    const std::vector<Vec2> pins = pins_flat(c, g, X);
    const FlexBasis fb = flex_basis(g, pins);
    const ConeSystem cs = cone_system(c, X);
    Eigen::MatrixXd M0 = cs.A * fb.N;
    const int nr = static_cast<int>(M0.rows());
    std::vector<double> rn(nr); int noise = 0;
    std::vector<int> all, okr;
    for (int i = 0; i < nr; ++i) { rn[i] = M0.row(i).norm(); all.push_back(i); if (rn[i] < 1e-12) ++noise; else okr.push_back(i); }
    Eigen::MatrixXd M = M0; normalise_rows(M);
    const Eigen::VectorXd z = fb.N.transpose() * sigma_flex(c, X);
    const Eigen::VectorXd a = M * z;
    const std::vector<char> keep = chart(cs, all, a);
    double mn = 1e300; int noise_in_chart = 0;
    for (int i = 0; i < nr; ++i) if (keep[i]) { mn = std::min(mn, a(i)); if (rn[i] < 1e-12) ++noise_in_chart; }
    const double sigma_chart = mn / z.norm();
    Eigen::MatrixXd Mc(okr.size(), M0.cols()); for (size_t i = 0; i < okr.size(); ++i) Mc.row(i) = M0.row(okr[i]);
    normalise_rows(Mc);
    const Eigen::VectorXd ac = Mc * z;
    const std::vector<char> keepc = chart(cs, okr, ac);
    double mnc = 1e300; for (size_t i = 0; i < okr.size(); ++i) if (keepc[i]) mnc = std::min(mnc, ac(i));
    ExpansiveConeOptions opt; opt.lp.dual_iters = 600;
    const ExpansiveConeReport rep = expansive_cone(c, X, opt);
    std::vector<int> ki; for (int i = 0; i < nr; ++i) if (keep[i]) ki.push_back(i);
    Eigen::MatrixXd Mk(ki.size(), M.cols()); for (size_t i = 0; i < ki.size(); ++i) Mk.row(i) = M.row(ki[i]);
    ConeLPOptions o600; o600.dual_iters = 600;
    const ConeLPResult r600 = cone_lp(Mk, o600), r20k = cone_lp(Mk);
    std::vector<int> kc; for (size_t i = 0; i < okr.size(); ++i) if (keepc[i]) kc.push_back(static_cast<int>(i));
    Eigen::MatrixXd Mkc(kc.size(), M.cols()); for (size_t i = 0; i < kc.size(); ++i) Mkc.row(i) = Mc.row(kc[i]);
    const ConeLPResult rc600 = cone_lp(Mkc, o600), rc20k = cone_lp(Mkc);
    json o; o["id"] = r["id"]; o["kind"] = r["kind"]; o["n_rows"] = nr; o["noise"] = noise; o["noise_in_chart"] = noise_in_chart;
    o["sigma_chart"] = sigma_chart; o["sigma_chart_clean"] = mnc / z.norm();
    o["rep_margin"] = rep.margin_l2; o["rep_dual"] = rep.dual_bound; o["rep_dual_max"] = rep.dual_bound_max;
    o["rep_n_active"] = rep.n_active; o["rep_n_dual_support"] = rep.n_dual_support; o["rep_pass"] = rep.pass_feasible;
    o["rep_sigma_chart"] = rep.sigma_chart_margin;
    o["chart600"] = {r600.margin_l2, r600.dual_bound}; o["chart20k"] = {r20k.margin_l2, r20k.dual_bound};
    o["clean600"] = {rc600.margin_l2, rc600.dual_bound}; o["clean20k"] = {rc20k.margin_l2, rc20k.dual_bound};
    // Gram of the clean (non-noise) normalised rows, basis-invariant
    Eigen::MatrixXd G = Mc * Mc.transpose(); json Gj = json::array();
    for (int i = 0; i < G.rows(); ++i) { json row = json::array(); for (int j = 0; j < G.cols(); ++j) row.push_back(G(i, j)); Gj.push_back(row); }
    o["gram_clean"] = Gj; json rnj = json::array(); for (double v : rn) rnj.push_back(v); o["row_norms"] = rnj;
    json aj = json::array(); for (int i = 0; i < a.size(); ++i) aj.push_back(a(i)); o["a_sigma"] = aj;
    out.push_back(o);
    std::cout << r["kind"].get<std::string>() << "_" << r["id"] << ": rows=" << nr << " noise=" << noise << " in-chart=" << noise_in_chart
              << " sigma_chart=" << sigma_chart << " clean=" << mnc / z.norm() << " rep: margin=" << rep.margin_l2 << " dual=" << rep.dual_bound << " pass=" << rep.pass_feasible << "\n"
              << "    chart LP 600: " << r600.margin_l2 << " / " << r600.dual_bound << " | 20k: " << r20k.margin_l2 << " / " << r20k.dual_bound << "\n"
              << "    clean LP 600: " << rc600.margin_l2 << " / " << rc600.dual_bound << " | 20k: " << rc20k.margin_l2 << " / " << rc20k.dual_bound << "\n";
  }
  std::ofstream f(argv[2]); f << out.dump() << "\n";
}
