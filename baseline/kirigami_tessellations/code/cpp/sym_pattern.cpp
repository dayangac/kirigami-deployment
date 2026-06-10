#include "sym_pattern.h"
#include "kiri_mesh.h"
#include <queue>
#include <vector>
#include "json.hpp"
#include <fstream>
#include <Optiz/NewtonSolver/Problem.h>
namespace geom {

SymPattern &SymPattern::add_vert(const Eigen::Vector2d &v) {
  verts.push_back(Vertex{.index = verts.size(), .coords = v, .parent = this});
  return *this;
}

SymPattern &SymPattern::add_edge(const Eigen::Vector2i &e) {
  edges.push_back(Edge{.index = edges.size(), .verts = e, .parent = this});
  return *this;
}

Eigen::MatrixXd SymPattern::get_verts() const {
  Eigen::MatrixXd V(verts.size(), 2);
  for (int i = 0; i < verts.size(); i++) {
    V.row(i) = verts[i].coords;
  }
  return V;
}

utils::Hmesh SymPattern::mesh(int rep_x, int rep_y) const {
  Eigen::MatrixXd V = get_verts();
  int count = 0;
  int nv = V.rows();
  auto get_faces_rep = [&](int rep) {
    std::vector<std::vector<int>> res = faces;
    for (int i = 0; i < res.size(); i++) {
      for (int j = 0; j < res[i].size(); j++) {
        res[i][j] += rep * nv;
      }
    }
    return res;
  };
  Eigen::Matrix2d P = periodicity;

  Eigen::MatrixXd V_rep(nv * rep_x * rep_y, V.cols());
  std::vector<std::vector<int>> faces_rep;
  for (int x = 0; x < rep_x; x++) {
    for (int y = 0; y < rep_y; y++) {
      Eigen::MatrixXd new_v = V;
      new_v.rowwise() += P.row(0) * x + P.row(1) * y;
      V_rep.block(count * nv, 0, nv, V.cols()) = new_v;
      auto new_faces = get_faces_rep(count);
      faces_rep.insert(faces_rep.end(), new_faces.begin(), new_faces.end());
      count++;
    }
  }
  // return utils::Hmesh(get_verts(), faces);
  return utils::Hmesh(V_rep, faces_rep);
}

SymPattern SymPattern::make_friendly() const {
  auto mesh = this->mesh().merge_close_verts();
return *this;
}

int SymPattern::find_next_edge(int edge) {
  Eigen::Vector2d edge_vec = -edges[edge].vec().normalized();
  Eigen::Vector2d rot_vec = Eigen::Vector2d{-edge_vec(1), edge_vec(0)};
  int best_edge = -1;
  int target_v = edges[edge].verts(1);
  double largest_ang = -1;
  for (auto [v, e] : vert_edges[target_v]) {
    Eigen::Vector2d next_vec = edges[e].vec();
    double ang = std::atan2(next_vec.dot(rot_vec), next_vec.dot(edge_vec));
    if (ang < -1e-3)
      ang += 2 * M_PI;
    if (ang > largest_ang) {
      largest_ang = ang;
      best_edge = e;
    }
  }
  return best_edge;
}

Eigen::Matrix2d SymPattern::get_bounding_parallelogram(double opening_angle) const {
  Eigen::RowVector2d xx, yy;
  {
    auto mesh = this->mesh(2, 1).merge_close_verts();
    geom::KiriMesh kiri_mesh(mesh);
    kiri_mesh.generate_cut_mesh();
    kiri_mesh.generate_rotating_info();
    auto verts = kiri_mesh.open(opening_angle);
    xx = verts.row(kiri_mesh.cut_mesh.F[faces.size()][0]) -
         verts.row(kiri_mesh.cut_mesh.F[0][0]);
  }
  {
    auto mesh = this->mesh(1, 2).merge_close_verts();
    geom::KiriMesh kiri_mesh(mesh);
    kiri_mesh.generate_cut_mesh();
    kiri_mesh.generate_rotating_info();
    auto verts = kiri_mesh.open(opening_angle);
    yy = verts.row(kiri_mesh.cut_mesh.F[faces.size()][0]) -
         verts.row(kiri_mesh.cut_mesh.F[0][0]);
  }
  Eigen::Matrix2d basis;
  basis.row(0) = xx;
  basis.row(1) = yy;
  return basis;
}

SymPattern &SymPattern::determine_cut_edges() {
  std::set<int> used_edges;
  int start_edge = 0;
  std::queue<std::pair<int, bool>> q;
  q.push({start_edge, true});
  while (!q.empty()) {
    auto [e, is_cut] = q.front();
    q.pop();
    if (used_edges.find(e) != used_edges.end()) {
      continue;
    }
    used_edges.insert(e);
    edges[e].is_cut = is_cut;
    q.push({edges[e].next, is_cut});
    q.push({edges[e].prev, is_cut});
    q.push({edges[e].twin, !is_cut});
  }
  return *this;
}

SymPattern &SymPattern::build_data_structure() {
  vert_edges.clear();
  vert_edges.resize(verts.size());
  for (auto &e : edges) {
    vert_edges[e.verts(0)][e.verts(1)] = e.index;
  }
  // Assign twin edges.
  for (int i = 0; i < vert_edges.size(); i++) {
    for (auto [v, e] : vert_edges[i]) {
      // Add missing twin edge.
      if (vert_edges[v].find(i) == vert_edges[v].end()) {
        edges.push_back(Edge{.index = edges.size(),
                             .verts = Eigen::Vector2i(v, i),
                             .parent = this});
        vert_edges[v][i] = edges.size() - 1;
      }
      edges[e].twin = vert_edges[v][i];
      edges[vert_edges[v][i]].twin = e;
    }
  }
  // Assign next and prev edges.
  for (int i = 0; i < edges.size(); i++) {
    edges[i].next = find_next_edge(i);
    edges[edges[i].next].prev = i;
  }
  determine_cut_edges();
  build_faces();
  return *this;
}

SymPattern &SymPattern::add_edges_w_threshold(double threshold) {
  for (int i = 0; i < verts.size(); i++) {
    for (int j = i + 1; j < verts.size(); j++) {
      if ((verts[i].coords - verts[j].coords).norm() < threshold) {
        add_edge(i, j);
      }
    }
  }
  return *this;
}

void SymPattern::build_faces() {
  faces.clear();
  outer_faces.clear();
  std::vector<bool> visited(edges.size(), false);
  for (int i = 0; i < edges.size(); i++) {
    if (visited[i])
      continue;
    std::vector<int> face;
    int e = i;
    double sum_angles = 0;
    do {
      face.push_back(edges[e].verts(0));
      visited[e] = true;
      double next_e = edges[e].next;
      Eigen::Vector2d edge_vec = edges[e].vec().normalized();
      Eigen::Vector2d next_edge_vec = edges[next_e].vec().normalized();
      double ang = std::atan2(edge_vec(0) * next_edge_vec(1) -
                                  edge_vec(1) * next_edge_vec(0),
                              edge_vec.dot(next_edge_vec));
      sum_angles += ang;
      e = edges[e].next;
    } while (e != i);
    // Don't add the outer face.
    if (sum_angles > 0) {
      faces.push_back(face);
    } else {
      outer_faces.push_back(face);
    }
  }
}

SymPattern SymPattern::from_json_graph(const std::string &json_path) {
  std::ifstream file(json_path);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + json_path);
  }
  std::allocator<std::map<std::string, nlohmann::basic_json<>>> alloc;
  nlohmann::json j;
  try {
    file >> j;
  } catch (const std::exception& e) {
    throw std::runtime_error("JSON parse error: " + std::string(e.what()));
  }
  
  SymPattern pattern;
  // Parse the JSON and populate the pattern
  // JSON format: {"vertices":[{"x":-20,"y":0},{"x":-20,"y":600}, ...], "edges":[{"v1":0,"v2":2"}, {"v1":0,"v2":3}, ...]}
  for (const auto& vertex : j["vertices"]) {
    pattern.add_vert(Eigen::Vector2d(vertex["x"], vertex["y"]));
  }
  for (const auto& edge : j["edges"]) {
    pattern.add_edge(Eigen::Vector2i(edge["v1"], edge["v2"]));
  }
  pattern.build_data_structure();
  return pattern;
}

} // namespace geom