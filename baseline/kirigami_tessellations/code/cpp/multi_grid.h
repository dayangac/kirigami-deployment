#pragma once

#include <Eigen/Dense>

namespace geom {
struct SymPattern;

struct MultiGrid {
  // Represents a line segment, shifted along the normal infinitely.
  // The distance between each consecutive parallel lines is 'shift'.
  struct Line {
    Line(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1, const Eigen::Vector2d& shift_point_ref);
    Line(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1, double shift);
    Eigen::Vector2d p0;
    Eigen::Vector2d p1;
    Eigen::Vector2d dir;
    Eigen::Vector2d normal;
    double shift;

    Eigen::Vector2d intersect(const Line &other) const;
    Eigen::Vector2d intersection_t(const Line &other) const;
    double signed_distance(const Eigen::Vector2d &p) const;
    Line shift_by(double shift) const;
  };
  std::vector<Line> lines;

  // Calculate the periodicity along the first and second lines.
  Eigen::Vector2d calculate_periodicity() const;
  Eigen::Vector2d period_origin() const;

  SymPattern get_sym_pattern() const;

  static MultiGrid from_json(const std::string &json_path);
};

} // namespace geom
