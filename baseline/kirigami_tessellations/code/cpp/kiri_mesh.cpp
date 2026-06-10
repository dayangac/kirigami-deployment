#include "kiri_mesh.h"
#include "Hmesh.h"
#include "conversions.h"
#include <iostream>
#include <map>
#include <numeric>

namespace geom {

KiriMesh::KiriMesh() {}

KiriMesh::KiriMesh(const utils::Hmesh &mesh, bool merge) : mesh(mesh) {
  if (merge)
    this->mesh = this->mesh.merge_close_verts();
}

KiriMesh &KiriMesh::generate_cut_mesh() {
  std::vector<Eigen::VectorXd> V = convert::to_vec_mat(mesh.V);
  std::vector<std::vector<int>> F = mesh.F;

  std::map<utils::Hmesh::Edge *, bool> used_edges;
  std::queue<std::pair<utils::Hmesh::Edge *, bool>> q;
  q.push({&mesh.edge(0), true});
  while (!q.empty()) {
    auto [e, is_cut] = q.front();
    q.pop();
    if (used_edges.find(e) != used_edges.end()) {
      if (used_edges[e] != is_cut) {
        std::cout << "Error: edge " << e->vi
                  << " is cut in one direction and not in the other."
                  << std::endl;
      }
      continue;
    }
    used_edges[e] = is_cut;
    if (is_cut) {
      // Duplicate the vertex and split the edge at the origin.
      int new_vert = V.size();
      V.push_back(V[e->vi]);
      F[e->fi][e->fvi] = new_vert;
      if (e->prev()->twin()) {
        F[e->prev()->twin()->fi][e->prev()->twin()->fvi] = new_vert;
      }
    }

    q.push({e->next(), is_cut});
    q.push({e->prev(), is_cut});
    if (e->ti != -1) {
      q.push({e->twin(), !is_cut});
    }
  }
  cut_mesh = utils::Hmesh(convert::to_eig_mat(V), F).remove_unused_verts();
  return *this;
}

static double cross2d(const Eigen::Vector2d &a, const Eigen::Vector2d &b) {
  return a.x() * b.y() - a.y() * b.x();
}

double KiriMesh::max_opening_angle() {
  double max_ang = 2 * M_PI;
  for (auto &ve : cut_mesh.verts_edges) {
    if (ve.size() != 2)
      continue;
    auto it = ve.begin();
    auto &e0 = cut_mesh.edges[it->second];
    auto e0_vec = e0.vec();
    auto e0_prev = e0.prev();
    auto e0_prev_vec = (-e0_prev->vec()).eval();
    double ang1 =
        std::atan2(cross2d(e0_vec, e0_prev_vec), e0_vec.dot(e0_prev_vec));

    it++;
    auto &e1 = cut_mesh.edges[it->second];
    auto e1_vec = e1.vec();
    auto e1_prev = e1.prev();
    auto e1_prev_vec = (-e1_prev->vec()).eval();
    double ang2 =
        std::atan2(cross2d(e1_vec, e1_prev_vec), e1_vec.dot(e1_prev_vec));
    double ang = 2 * M_PI - (ang1 + ang2);
    if (ang < max_ang) {
      max_ang = ang;
    }
  }
  return max_ang;
}

KiriMesh &KiriMesh::generate_rotating_info() {
  rotating_info.clear();
  rotating_info.resize(cut_mesh.nf());
  vertex_rot_info.clear();
  vertex_rot_info.resize(cut_mesh.nv(), -1);
  std::set<int> used_faces;
  std::vector<bool> visited(cut_mesh.nv(), false);
  auto assign_rot_inf = [&](int f, const Eigen::Vector2d &center,
                            bool clockwise, int prev_face) {
    if (used_faces.find(f) != used_faces.end()) {
      return;
    }
    rotating_info[f] = {center, clockwise, prev_face};
    for (int i = 0; i < cut_mesh.F[f].size(); ++i) {
      int v = cut_mesh.F[f][i];
      if (visited[v])
        continue;
      vertex_rot_info[v] = f;
      visited[v] = true;
    }
  };
  auto is_cut_edge = [&](utils::Hmesh::Edge *e) {
    return e->twin() && cut_mesh.F[e->fi][e->fvi] !=
                            cut_mesh.F[e->twin()->fi][e->twin()->next()->fvi];
  };
  // Face and previous face.
  std::queue<int> q;
  q.push(0);
  // Vertices of first face don't need to be rotated.
  for (auto v : cut_mesh.F[0])
    visited[v] = true;
  // Find the rotation sequence for each face.
  while (!q.empty()) {
    int f = q.front();
    q.pop();
    if (used_faces.find(f) != used_faces.end()) {
      continue;
    }
    used_faces.insert(f);
    utils::Hmesh::Edge *e = mesh.face(f).edge();
    do {
      if (e->twin()) {
        int next_face = e->twin()->fi;
        if (is_cut_edge(e)) {
          assign_rot_inf(next_face, e->next()->origin()->coords(), false, f);
        } else {
          assign_rot_inf(next_face, e->origin()->coords(), true, f);
        }
        q.push(next_face);
      }
      e = e->next();
    } while (e != mesh.face(f).edge());
  }

  return *this;
}

Eigen::MatrixXd KiriMesh::open(double angle) {
  Eigen::MatrixXd V = cut_mesh.V;
  for (int i = 0; i < V.rows(); i++) {
    V.row(i) = open_t(angle, i);
  }
  return V;
}

utils::Hmesh KiriMesh::open_separate_faces(double angle) {
  std::vector<std::vector<int>> faces;
  std::vector<Eigen::VectorXd> V;
  for (int i = 0; i < cut_mesh.nf(); i++) {
    std::vector<int> face;
    for (int j = 0; j < cut_mesh.F[i].size(); j++) {
      int v = cut_mesh.F[i][j];
      Eigen::Vector2d pos = open_t(angle, v, i);
      V.push_back(pos);
      face.push_back(V.size() - 1);
    }
    faces.push_back(face);
  }
  return utils::Hmesh(convert::to_eig_mat(V), faces);
}

Eigen::VectorXd KiriMesh::get_friendly_vertices(utils::Hmesh &separate_faces) {
  Eigen::VectorXd res = Eigen::VectorXd::Zero(separate_faces.nv());
  auto next_vertex_edge = [&](utils::Hmesh::Edge *e) {
    if (!e->prev()->twin())
      return (utils::Hmesh::Edge *)nullptr;
    return e->prev()->twin();
  };
  auto is_cut_edge = [&](utils::Hmesh::Edge *e) {
    return e->twin() && cut_mesh.F[e->fi][e->fvi] !=
                            cut_mesh.F[e->twin()->fi][e->twin()->next()->fvi];
  };
  auto is_polygon_planar = [&](const std::vector<double> &poly_edges,
                               const std::vector<double> &poly_angles) {
    Eigen::Rotation2Dd rot(0);
    Eigen::Vector2d point(0, 0);
    for (size_t i = 0; i < poly_edges.size(); ++i) {
      point += rot * Eigen::Vector2d(poly_edges[i], 0);
      rot *= Eigen::Rotation2Dd(poly_angles[i]);
    }
    return point.norm() < 1e-5;
  };

  for (auto &v : mesh.verts) {
    if (mesh.is_boundary_vertex(v.index))
      continue;

    std::vector<double> poly_edges;
    std::vector<double> poly_angles;
    std::vector<Eigen::VectorXd> incoming_edges, outgoing_edges;
    auto e = v.edge();
    if (!is_cut_edge(e))
      e = next_vertex_edge(e);
    auto orig_e = e;
    int iter = 0;
    do {
      poly_edges.push_back(e->vec().norm());
      incoming_edges.push_back(e->vec());
      auto ne = next_vertex_edge(e), nne = next_vertex_edge(ne);
      outgoing_edges.push_back(-ne->vec());
      double angle =
          std::acos(e->vec().normalized().dot(ne->vec().normalized())) +
          std::acos(ne->vec().normalized().dot(nne->vec().normalized()));
      poly_angles.push_back(angle);
      e = nne;
    } while (e != orig_e && iter++ < 100);

    if (iter > 10) {
      std::cerr << "Error: infinite loop in get_friendly_vertices()"
                << std::endl;
      exit(1);
    }
    auto sum_edges = [](const std::vector<Eigen::VectorXd> &edges) {
      Eigen::VectorXd res = Eigen::VectorXd::Zero(edges[0].size());
      for (const auto &e : edges) {
        res += e;
      }
      return res;
    };
    if (!is_polygon_planar(poly_edges, poly_angles)) {
      // res(v.index) = 0;
      auto e = v.edge();
      do {
        if (!is_cut_edge(e)) {
          res(separate_faces.F[e->fi][e->fvi]) = 1;
        }
        e = next_vertex_edge(e);
      } while (e != v.edge());
    }
  }
  return res;
}

} // namespace geom