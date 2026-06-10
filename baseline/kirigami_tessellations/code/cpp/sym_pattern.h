#pragma once
#include "Hmesh.h"
#include <Eigen/Eigen>
#include <vector>
namespace geom {

struct SymPattern {

  struct Vertex {
    size_t index;
    Eigen::Vector2d coords;
    SymPattern *parent;
  };
  struct Edge {
    size_t index;
    size_t twin;
    size_t next;
    size_t prev;
    bool is_cut; // Whether the edge is a cut edge (if it's not, its twin is).
    Eigen::Vector2i verts;
    SymPattern *parent;
    Eigen::Vector2d vec() const {
      return parent->verts[verts(1)].coords - parent->verts[verts(0)].coords;
    }
  };

  SymPattern &add_vert(const Eigen::Vector2d &v);
  SymPattern &add_vert(double x, double y) {
    return add_vert(Eigen::Vector2d(x, y));
  }
  SymPattern &add_edge(const Eigen::Vector2i &e);
  SymPattern &add_edge(int a, int b) { return add_edge(Eigen::Vector2i(a, b)); }
  SymPattern &add_edges_w_threshold(double threshold);

  int find_next_edge(int edge);
  SymPattern &build_data_structure();
  void build_faces();
  SymPattern &determine_cut_edges();
  // Get the mesh of the pattern.
  utils::Hmesh mesh(int rep_x = 2, int rep_y = 2) const;

  SymPattern make_friendly() const;

  Eigen::Matrix2d get_bounding_parallelogram(double opening_angle = 0.0) const;

  static SymPattern arrow() {
    geom::SymPattern pattern;
    pattern.add_vert(0, 0).add_vert(1, 0).add_vert(1, 1).add_vert(0, 1);
    pattern.add_edge(0, 1).add_edge(1, 2).add_edge(2, 3).add_edge(3, 0);
    pattern.add_edge(0, 2);
    pattern.build_data_structure();
    return pattern;
  }
  static SymPattern stretched_arrow() {
    geom::SymPattern pattern;
    pattern.add_vert(0, 0).add_vert(1.4, 0).add_vert(1.4, 1).add_vert(0, 1);
    pattern.add_edge(0, 1).add_edge(1, 2).add_edge(2, 3).add_edge(3, 0);
    pattern.add_edge(0, 2);
    pattern.periodicity << 1.4, 0, 0, 1;
    pattern.build_data_structure();
    return pattern;
  }
  static SymPattern water_bomb() {
    geom::SymPattern pattern;
    pattern.add_vert(0, 0)
        .add_vert(1, 0)
        .add_vert(1, 1)
        .add_vert(0, 1)
        .add_vert(2, 0)
        .add_vert(2, 1)
        .add_vert(2, 2)
        .add_vert(1, 2)
        .add_vert(0, 2);
    pattern.add_edges_w_threshold(1 + 1e-2);
    pattern.periodicity << 2, 0, 0, 2;
    pattern.build_data_structure();
    return pattern;
  }
  static SymPattern stretched_water_bomb() {
    geom::SymPattern pattern;
    pattern.add_vert(0, 0)
        .add_vert(1, 0)
        .add_vert(1, 1)
        .add_vert(0, 1)
        .add_vert(2, 0)
        .add_vert(2, 1)
        .add_vert(2, 2)
        .add_vert(1, 2)
        .add_vert(0, 2);
    pattern.add_edges_w_threshold(1 + 1e-2);
    for (int i = 0; i < pattern.verts.size(); i++) {
      pattern.verts[i].coords.x() *= 1.4;
    }
    pattern.periodicity << 2 * 1.4, 0, 0, 2;
    pattern.build_data_structure();
    return pattern;
  }
  static SymPattern sheared_water_bomb() {
    geom::SymPattern pattern;
    pattern.add_vert(0, 0)
        .add_vert(1, 0)
        .add_vert(1, 1)
        .add_vert(0, 1)
        .add_vert(2, 0)
        .add_vert(2, 1)
        .add_vert(2, 2)
        .add_vert(1, 2)
        .add_vert(0, 2);
    pattern.add_edges_w_threshold(1 + 1e-2);
    for (int i = 0; i < pattern.verts.size(); i++) {
      pattern.verts[i].coords.x() += pattern.verts[i].coords.y() * 0.4;
    }
    pattern.periodicity << 2, 0, 2 * 0.4, 2;
    pattern.build_data_structure();
    return pattern;
  }
  static SymPattern triangle_auxetic() {
    geom::SymPattern pattern;
    pattern.add_vert(0, 0)
        .add_vert(1, 0)
        .add_vert(std::cos(M_PI / 3), std::sin(M_PI / 3))
        .add_vert(-std::cos(M_PI / 3), std::sin(M_PI / 3));
    pattern.add_edges_w_threshold(1 + 1e-2);
    pattern.periodicity << 1, 0, pattern.verts[3].coords.transpose();
    pattern.build_data_structure();
    return pattern;
  }
  static SymPattern hexagon_traignle_pattern() {
    geom::SymPattern pattern;
    pattern.add_vert(0, 0)
        .add_vert(1, 0)
        .add_vert(std::cos(M_PI / 3), std::sin(M_PI / 3))
        .add_vert(-1, 0)
        .add_vert(-1 - std::cos(M_PI / 3), std::sin(M_PI / 3))
        .add_vert(-1 - 2 * std::cos(M_PI / 3), 2 * std::sin(M_PI / 3))
        .add_vert(-1, 2 * std::sin(M_PI / 3))
        .add_vert(0, 2 * std::sin(M_PI / 3));
    for (int i = 0; i < pattern.verts.size(); i++) {
      pattern.verts[i].coords.x() += 1.0;
    }
    pattern.add_edges_w_threshold(1 + 1e-2);
    pattern.periodicity << 2, 0, -2 * std::cos(M_PI / 3),
        2 * std::sin(M_PI / 3);
    pattern.build_data_structure();
    return pattern;
  }

  static SymPattern from_json_graph(const std::string &json_path);



  std::vector<Vertex> verts;
  std::vector<Edge> edges;
  std::vector<std::vector<int>> faces;
  std::vector<std::vector<int>> outer_faces;
  Eigen::Matrix2d periodicity = Eigen::Matrix2d::Identity();
  Eigen::MatrixXd get_verts() const;

  // For each pair of vertices, maps the index of the edge connecting them.
  std::vector<std::map<int, int>> vert_edges;
};

} // namespace geom