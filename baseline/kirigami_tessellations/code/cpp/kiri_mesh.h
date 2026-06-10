#pragma once
#include "Hmesh.h"

namespace geom {
struct KiriMesh {
  KiriMesh();
  KiriMesh(const utils::Hmesh &mesh, bool merge = true);
  // Closed mesh.
  utils::Hmesh mesh;
  utils::Hmesh cut_mesh;

  // Open mesh.
  KiriMesh &generate_cut_mesh();
  KiriMesh &generate_rotating_info();

  struct RotatingInfo {
    Eigen::Vector2d center;
    bool clockwise;
    int prev_face = -1;
  };
  // Per vertex rotating info.
  std::vector<RotatingInfo> rotating_info;
  std::vector<int> vertex_rot_info;

  double max_opening_angle();

  Eigen::VectorXd get_friendly_vertices(utils::Hmesh &separate_faces);

  Eigen::MatrixXd open(double angle);
  utils::Hmesh open_separate_faces(double angle);
  template <typename T> Eigen::Vector2<T> open_t(const T &ang, int v) {
    Eigen::Vector2<T> pos(cut_mesh.V(v, 0), cut_mesh.V(v, 1));
    for (int f = vertex_rot_info[v];
         f != -1 && rotating_info[f].prev_face != -1;
         f = rotating_info[f].prev_face) {
      auto rotate =
          Eigen::Rotation2D<T>(ang * (rotating_info[f].clockwise ? -1 : 1));
      pos = rotate.toRotationMatrix() * (pos - rotating_info[f].center) +
            rotating_info[f].center;
    }
    return pos;
  }
  template <typename T> Eigen::Vector2<T> open_t(const T &ang, int v, int f) {
    Eigen::Vector2<T> pos(cut_mesh.V(v, 0), cut_mesh.V(v, 1));
    for (; f != -1 && rotating_info[f].prev_face != -1;
         f = rotating_info[f].prev_face) {
      auto rotate =
          Eigen::Rotation2D<T>(ang * (rotating_info[f].clockwise ? -1 : 1));
      pos = rotate.toRotationMatrix() * (pos - rotating_info[f].center) +
            rotating_info[f].center;
    }
    return pos;
  }
};
} // namespace geom