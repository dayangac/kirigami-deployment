#include "multi_grid.h"
#include "json.hpp"
#include "sym_pattern.h"
#include "utils.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <limits>

namespace geom {

MultiGrid::Line::Line(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1,
                      double shift)
    : p0(p0), p1(p1), dir((p1 - p0).normalized()), normal(-dir.y(), dir.x()),
      shift(shift) {}

MultiGrid::Line::Line(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1,
                      const Eigen::Vector2d &p2)
    : p0(p0), p1(p1), dir((p1 - p0).normalized()), normal(-dir.y(), dir.x()) {
  Eigen::Vector2d dir = (p1 - p0).normalized();
  Eigen::Vector2d normal(-dir.y(), dir.x());
  shift = std::abs((p2 - p0).dot(normal));
}

Eigen::Vector2d MultiGrid::calculate_periodicity() const {
  double t0 = lines[1].shift / lines[0].dir.dot(lines[1].normal);
  double t1 = lines[0].shift / lines[1].dir.dot(lines[0].normal);
  // Intersection every inter + t0 * lines[0].dir.
  for (int i = 2; i < lines.size(); i++) {
    double t = lines[i].shift / lines[0].dir.dot(lines[i].normal);
    // Find n, m such t0 / t ~ n / m.
    auto [n, m] = utils::approximate_ratio(std::abs(t0 / t));
    if (n < 10 && m < 10)
      t0 *= m;
    else
      return Eigen::Vector2d(100, 100);

    t = lines[i].shift / lines[1].dir.dot(lines[i].normal);
    std::tie(n, m) = utils::approximate_ratio(std::abs(t1 / t));
    if (n < 10 && m < 10)
      t1 *= m;
    else
      return Eigen::Vector2d(100, 100);
  }
  return Eigen::Vector2d(t0, t1);
}

Eigen::Vector2d MultiGrid::period_origin() const {
  return lines[0].intersect(lines[1]);
}

SymPattern MultiGrid::get_sym_pattern() const {
  SymPattern pattern;
  Eigen::Vector2d origin = period_origin(),
                  periodicity = calculate_periodicity();
  pattern.periodicity << periodicity(0) * lines[0].dir.transpose(),
      periodicity(1) * lines[1].dir.transpose();

  // The 4 points of the unit cell.
  std::vector<Eigen::Vector2d> points;
  points.push_back(origin);
  points.push_back(origin + periodicity(0) * lines[0].dir);
  points.push_back(origin + periodicity(1) * lines[1].dir);
  points.push_back(origin + periodicity(0) * lines[0].dir +
                   periodicity(1) * lines[1].dir);
  auto cross2d = [](const Eigen::Vector2d &v0, const Eigen::Vector2d &v1) {
    return v0.x() * v1.y() - v0.y() * v1.x();
  };
  auto side = [&](const Eigen::Vector2d &v0, const Eigen::Vector2d &v1,
                  const Eigen::Vector2d &p) {
    // Returns true if p is on the left side of the line v0->v1.
    return cross2d(v1 - v0, p - v0);
  };
  auto in_parallelogram = [&](const Eigen::Vector2d &p) {
    return side(points[0], points[1], p) *
                   side(points[0], points[1], points[2]) >=
               -1e-4 &&
           side(points[0], points[2], p) *
                   side(points[0], points[2], points[3]) >=
               -1e-4 &&
           side(points[3], points[1], p) *
                   side(points[3], points[1], points[2]) >=
               -1e-4 &&
           side(points[2], points[3], p) *
                   side(points[2], points[3], points[0]) >=
               -1e-4;
  };
  auto minmax_reps = [&](const Line &line) {
    double min_dist = std::numeric_limits<double>::infinity();
    double max_dist = -std::numeric_limits<double>::infinity();
    for (const auto &p : points) {
      double dist = line.signed_distance(p);
      min_dist = std::min(min_dist, dist);
      max_dist = std::max(max_dist, dist);
    }
    // return std::make_pair(min_dist, max_dist);
    int min_rep = std::ceil((min_dist - 1e-3) / line.shift);
    int max_rep = std::floor((max_dist + 1e-3) / line.shift);
    return std::make_pair(min_rep, max_rep);
  };

  // Find all intersections along the lines in the range.
  std::map<std::pair<int, int>, std::vector<double>> t_matrix;
  for (int i = 0; i < lines.size(); i++) {
    auto [min_rep, max_rep] = minmax_reps(lines[i]);
    for (int j = i + 1; j < lines.size(); j++) {
      auto [min_rep2, max_rep2] = minmax_reps(lines[j]);
      for (int rep = min_rep; rep <= max_rep; rep++) {
        for (int rep2 = min_rep2; rep2 <= max_rep2; rep2++) {
          auto line1 = lines[i].shift_by(rep * lines[i].shift);
          auto line2 = lines[j].shift_by(rep2 * lines[j].shift);
          auto ts = line1.intersection_t(line2);
          Eigen::Vector2d inter = line1.p0 + ts[0] * line1.dir;
          if (!in_parallelogram(inter)) {
            continue;
          }
          t_matrix[{i, rep}].push_back(ts[0]);
          t_matrix[{j, rep2}].push_back(ts[1]);
        }
      }
    }
  }
  // Sort the t values for each line.
  for (auto &[j, ts] : t_matrix) {
    // Remove duplicates.
    std::sort(ts.begin(), ts.end());
    ts.erase(
        std::unique(ts.begin(), ts.end(),
                    [](double a, double b) { return std::abs(a - b) < 1e-4; }),
        ts.end());
  }

  // Add vertices and edges to the pattern.
  auto vert_index = [&](const Eigen::Vector2d &p) {
    for (int i = 0; i < pattern.verts.size(); i++) {
      if ((pattern.verts[i].coords - p).norm() < 1e-4) {
        return i;
      }
    }
    pattern.add_vert(p);
    return (int)pattern.verts.size() - 1;
  };

  for (const auto &[j, ts] : t_matrix) {
    Line l = lines[j.first].shift_by(j.second * lines[j.first].shift);
    for (int i = 0; i < ts.size() - 1; i++) {
      Eigen::Vector2d p0 = l.p0 + ts[i] * l.dir;
      Eigen::Vector2d p1 = l.p0 + ts[i + 1] * l.dir;
      int v0 = vert_index(p0);
      int v1 = vert_index(p1);
      pattern.add_edge(v0, v1);
    }
  }
  pattern.build_data_structure();
  return pattern;
}

MultiGrid MultiGrid::from_json(const std::string &json_path) {
  std::ifstream file(json_path);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + json_path);
  }
  nlohmann::json j;
  file >> j;
  MultiGrid mg;

  for (const auto &line_data : j) {
    const auto &points = line_data["points"];
    if (points.size() != 3) {
      throw std::runtime_error("Each line must have exactly 3 points");
    }

    // Extract p0 and p1
    Eigen::Vector2d p0(points[0]["x"], points[0]["y"]);
    Eigen::Vector2d p1(points[1]["x"], points[1]["y"]);
    Eigen::Vector2d p2(points[2]["x"], points[2]["y"]);
    p0 /= 100;
    p1 /= 100;
    p2 /= 100;

    // Calculate the line direction
    Eigen::Vector2d dir = (p1 - p0).normalized();

    // Calculate the normal vector
    Eigen::Vector2d normal(-dir.y(), dir.x());

    // Calculate the shift as the distance from p2 to the line
    double shift = std::abs((p2 - p0).dot(normal));

    // Create and add the line
    mg.lines.emplace_back(p0, p1, shift);
  }

  return mg;
}

Eigen::Vector2d MultiGrid::Line::intersect(const Line &other) const {
  // Get direction vectors and a point on each line
  Eigen::Vector2d d1 = dir;
  Eigen::Vector2d d2 = other.dir;
  Eigen::Vector2d p1 = p0;
  Eigen::Vector2d p2 = other.p0;

  // Solve the system of equations:
  // p1 + t*d1 = p2 + s*d2
  // [d1 -d2][t] = p2 - p1
  //               [s]
  Eigen::Matrix2d A;
  A << d1, -d2;
  Eigen::Vector2d b = p2 - p1;
  Eigen::Vector2d ts = A.inverse() * b;
  double t = ts[0];

  // Return intersection point
  return p1 + t * d1;
}

Eigen::Vector2d MultiGrid::Line::intersection_t(const Line &other) const {
  // Get direction vectors and a point on each line
  Eigen::Vector2d d1 = dir;
  Eigen::Vector2d d2 = other.dir;
  Eigen::Vector2d p1 = p0;
  Eigen::Vector2d p2 = other.p0;

  // Solve the system of equations:
  // p1 + t*d1 = p2 + s*d2
  // [d1 -d2][t] = p2 - p1
  //               [s]
  Eigen::Matrix2d A;
  A << d1, -d2;
  Eigen::Vector2d b = p2 - p1;
  Eigen::Vector2d ts = A.inverse() * b;
  return ts;
}

double MultiGrid::Line::signed_distance(const Eigen::Vector2d &p) const {
  return (p - p0).dot(normal);
}

MultiGrid::Line MultiGrid::Line::shift_by(double shift) const {
  return Line(p0 + shift * normal, p1 + shift * normal, shift);
}

} // namespace geom
