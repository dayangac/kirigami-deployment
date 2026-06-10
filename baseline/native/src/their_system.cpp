#include "their_system.h"

namespace their {

System build(const kirigami::UnitPattern &pattern) {
  auto &mesh = const_cast<utils::Hmesh &>(pattern.mesh);
  System s;

  // --- step 1/2/3: holes = vertex classes merged across split (equal-colour)
  // edges and across periodic identifications.
  std::vector<std::set<int>> holes;
  std::vector<int> vertex_to_hole(mesh.nv(), -1);
  for (int i = 0; i < mesh.nv(); i++) {
    holes.push_back(std::set<int>{i});
    vertex_to_hole[i] = holes.size() - 1;
  }
  auto merge_holes = [&](int v0, int v1) {
    int h1 = vertex_to_hole[v0], h2 = vertex_to_hole[v1];
    if (h1 > h2)
      std::swap(h1, h2);
    holes[h1].insert(holes[h2].begin(), holes[h2].end());
    for (int i : holes[h2])
      vertex_to_hole[i] = h1;
    holes[h2].clear();
  };
  for (auto &edge : mesh.edges) {
    if (!edge.twin())
      continue;
    if (pattern.face_colors[edge.fi] == pattern.face_colors[edge.twin()->fi] &&
        vertex_to_hole[edge.vi] != vertex_to_hole[edge.next()->vi])
      merge_holes(edge.vi, edge.next()->vi);
    if (edge.vi != edge.twin()->next()->vi &&
        vertex_to_hole[edge.vi] != vertex_to_hole[edge.twin()->next()->vi])
      merge_holes(edge.vi, edge.twin()->next()->vi);
  }
  std::vector<std::set<int>> non_empty;
  for (auto &h : holes)
    if (!h.empty()) {
      non_empty.push_back(h);
      for (int v : h)
        vertex_to_hole[v] = non_empty.size() - 1;
    }
  holes = non_empty;

  // --- step 5: assemble.
  std::vector<Eigen::Triplet<double>> triplets;
  int num_constraints = 0;
  std::vector<Eigen::Vector2d> rhs;
  for (auto &edge : mesh.edges) {
    if (!edge.twin())
      continue;
    if (pattern.face_colors[edge.fi] == 0 ||
        (pattern.face_colors[edge.fi] == 2 &&
         pattern.face_colors[edge.twin()->fi] == 1)) {
      if (!mesh.is_boundary_vertex(edge.next()->vi)) {
        int hole = vertex_to_hole[edge.next()->vi];
        triplets.emplace_back(hole, edge.vi, 1);
        triplets.emplace_back(hole, edge.next()->vi, -1);
      }
    } else {
      if (!mesh.is_boundary_vertex(edge.vi)) {
        int hole = vertex_to_hole[edge.vi];
        triplets.emplace_back(hole, edge.vi, -1);
        triplets.emplace_back(hole, edge.next()->vi, 1);
      }
    }
    if (edge.vi != edge.twin()->next()->vi) {
      Eigen::Vector2d diff =
          mesh.V.row(edge.vi) - mesh.V.row(edge.twin()->next()->vi);
      triplets.emplace_back(holes.size() + num_constraints, edge.vi, 1);
      triplets.emplace_back(holes.size() + num_constraints,
                            edge.twin()->next()->vi, -1);
      rhs.push_back(diff);
      num_constraints++;
    }
  }
  // map_boundary_to_unit_disk == false: pin vertex 0.
  triplets.emplace_back(holes.size() + num_constraints, 0, 1);
  rhs.push_back(mesh.V.row(0));
  num_constraints++;

  s.A = Eigen::SparseMatrix<double>(holes.size() + num_constraints, mesh.nv());
  s.A.setFromTriplets(triplets.begin(), triplets.end());
  s.b = Eigen::MatrixXd::Zero(holes.size() + num_constraints, 2);
  for (int i = 0; i < num_constraints; i++)
    s.b.row(holes.size() + i) = rhs[i];
  s.num_holes = holes.size();
  s.num_constraints = num_constraints;
  s.vertex_to_hole = vertex_to_hole;
  return s;
}

} // namespace their
