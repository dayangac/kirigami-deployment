#include "json_io.h"
#include <fstream>
#include <stdexcept>

namespace io {

int sigma_to_color(int sigma) { return sigma > 0 ? 1 : 0; }
int color_to_sigma(int color) { return color == 0 ? -1 : +1; }

static double signed_area(const Eigen::MatrixXd &V,
                          const std::vector<int> &f) {
  double a = 0;
  for (size_t i = 0; i < f.size(); i++) {
    const auto &p = V.row(f[i]);
    const auto &q = V.row(f[(i + 1) % f.size()]);
    a += p(0) * q(1) - q(0) * p(1);
  }
  return 0.5 * a;
}

Input load(const std::string &path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("cannot open " + path);
  nlohmann::json j;
  in >> j;

  const auto &jv = j.at("vertices");
  Eigen::MatrixXd V(jv.size(), 2);
  for (size_t i = 0; i < jv.size(); i++) {
    V(i, 0) = jv[i][0].get<double>();
    V(i, 1) = jv[i][1].get<double>();
  }
  std::vector<std::vector<int>> F;
  for (const auto &f : j.at("faces"))
    F.push_back(f.get<std::vector<int>>());

  Input res;
  // Our contract stores faces CCW; the sigma <-> colour map above is defined
  // relative to that.  Normalise so a hand-made input cannot silently flip the
  // meaning of sigma, and report how many faces had to be reversed.
  for (auto &f : F) {
    if (signed_area(V, f) < 0) {
      std::reverse(f.begin(), f.end());
      res.n_cw_faces_reversed++;
    }
  }

  std::vector<int> colors(F.size(), 0);
  if (j.contains("orientation")) {
    res.had_orientation = true;
    auto sig = j.at("orientation").get<std::vector<int>>();
    if (sig.size() != F.size())
      throw std::runtime_error("orientation length != face count");
    for (size_t i = 0; i < sig.size(); i++)
      colors[i] = sigma_to_color(sig[i]);
  }

  Eigen::Matrix2d P = Eigen::Matrix2d::Zero();
  if (j.contains("periodic")) {
    res.had_periodic = true;
    auto th = j["periodic"].at("th").get<std::vector<double>>();
    auto tv = j["periodic"].at("tv").get<std::vector<double>>();
    P.row(0) << th[0], th[1];
    P.row(1) << tv[0], tv[1];
  }

  res.pattern = kirigami::UnitPattern{utils::Hmesh(V, F), P, colors};
  return res;
}

nlohmann::json verts_to_json(const Eigen::MatrixXd &V) {
  nlohmann::json a = nlohmann::json::array();
  for (int i = 0; i < V.rows(); i++)
    a.push_back({V(i, 0), V(i, 1)});
  return a;
}

nlohmann::json faces_to_json(const std::vector<std::vector<int>> &F) {
  nlohmann::json a = nlohmann::json::array();
  for (const auto &f : F)
    a.push_back(f);
  return a;
}

nlohmann::json pattern_to_json(const kirigami::UnitPattern &p) {
  nlohmann::json j;
  j["vertices"] = verts_to_json(p.mesh.V);
  j["faces"] = faces_to_json(p.mesh.F);
  nlohmann::json sig = nlohmann::json::array();
  for (int c : p.face_colors)
    sig.push_back(c == 2 ? 0 : color_to_sigma(c));
  j["orientation"] = sig;
  if (p.periodicity.norm() > 1e-9) {
    j["periodic"]["th"] = {p.periodicity(0, 0), p.periodicity(0, 1)};
    j["periodic"]["tv"] = {p.periodicity(1, 0), p.periodicity(1, 1)};
  }
  return j;
}

void write(const std::string &path, const nlohmann::json &j) {
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("cannot write " + path);
  out << j.dump(2) << "\n";
}

} // namespace io
