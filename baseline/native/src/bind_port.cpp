// bind_port.cpp — the two routines the web UI's theta_max path needs, taken
// VERBATIM from baseline/tuttekiri/code/cpp/bind.cpp (lines 279-418), with the
// single change that `max_opening_angle` takes a kirigami::UnitPattern instead
// of an emscripten::val (the first two lines of its body, which only decoded
// the JS object, are replaced by the parameter).  Nothing else is edited; do
// not "clean up" this file.
#include "bind_port.h"
#include "geometry/kirigami.h"
#include "geometry/unit_pattern.h"
#include <Eigen/Eigen>
#include <algorithm>
#include <cmath>

static bool do_polygons_intersect(const Eigen::MatrixXd &p1,
                                  const Eigen::MatrixXd &p2) {
  auto orient = [](const Eigen::Vector2d &p, const Eigen::Vector2d &q,
                   const Eigen::Vector2d &r) {
    double val =
        (q.x() - p.x()) * (r.y() - p.y()) - (q.y() - p.y()) * (r.x() - p.x());
    if (std::abs(val) < 1e-9)
      return 0.0;
    return val;
  };

  auto on_segment = [&](const Eigen::Vector2d &p, const Eigen::Vector2d &a,
                        const Eigen::Vector2d &b) {
    if (orient(a, b, p) != 0.0)
      return false;
    return p.x() >= std::min(a.x(), b.x()) - 1e-9 &&
           p.x() <= std::max(a.x(), b.x()) + 1e-9 &&
           p.y() >= std::min(a.y(), b.y()) - 1e-9 &&
           p.y() <= std::max(a.y(), b.y()) + 1e-9;
  };

  auto intersect_proper =
      [&](const Eigen::Vector2d &a, const Eigen::Vector2d &b,
          const Eigen::Vector2d &c, const Eigen::Vector2d &d) {
        double o1 = orient(a, b, c);
        double o2 = orient(a, b, d);
        double o3 = orient(c, d, a);
        double o4 = orient(c, d, b);

        if (((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) &&
            ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0))) {
          return true;
        }
        return false;
      };

  auto is_strictly_inside = [&](const Eigen::Vector2d &p,
                                const Eigen::MatrixXd &poly) {
    bool inside = false;
    int n = poly.rows();
    for (int i = 0, j = n - 1; i < n; j = i++) {
      Eigen::Vector2d u = poly.row(i);
      Eigen::Vector2d v = poly.row(j);

      if (on_segment(p, u, v))
        return false;

      if (((u.y() > p.y()) != (v.y() > p.y())) &&
          (p.x() <
           (v.x() - u.x()) * (p.y() - u.y()) / (v.y() - u.y()) + u.x())) {
        inside = !inside;
      }
    }
    return inside;
  };

  int n1 = p1.rows();
  int n2 = p2.rows();

  for (int i = 0; i < n1; ++i) {
    Eigen::Vector2d a = p1.row(i);
    Eigen::Vector2d b = p1.row((i + 1) % n1);
    for (int j = 0; j < n2; ++j) {
      Eigen::Vector2d c = p2.row(j);
      Eigen::Vector2d d = p2.row((j + 1) % n2);
      if (intersect_proper(a, b, c, d))
        return true;

      // Check for collinear overlap
      if (orient(a, b, c) == 0.0 && orient(a, b, d) == 0.0) {
        Eigen::Vector2d u = b - a;
        double len_sq = u.squaredNorm();
        if (len_sq > 1e-12) {
          double sc = (c - a).dot(u) / len_sq;
          double sd = (d - a).dot(u) / len_sq;
          double s_min = std::max(0.0, std::min(sc, sd));
          double s_max = std::min(1.0, std::max(sc, sd));

          if (s_max - s_min > 1e-6) {
            Eigen::Vector2d mid = a + u * (s_min + s_max) * 0.5;
            Eigen::Vector2d normal(-u.y(), u.x());
            normal.normalize();
            double eps = 1e-4;
            if (is_strictly_inside(mid + normal * eps, p1)) {
              if (is_strictly_inside(mid + normal * eps, p2))
                return true;
            }
            if (is_strictly_inside(mid - normal * eps, p1)) {
              if (is_strictly_inside(mid - normal * eps, p2))
                return true;
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < n1; ++i) {
    if (is_strictly_inside(p1.row(i), p2))
      return true;
    Eigen::Vector2d a = p1.row(i);
    Eigen::Vector2d b = p1.row((i + 1) % n1);
    if (is_strictly_inside((a + b) * 0.5, p2))
      return true;
  }
  for (int i = 0; i < n2; ++i) {
    if (is_strictly_inside(p2.row(i), p1))
      return true;
    Eigen::Vector2d a = p2.row(i);
    Eigen::Vector2d b = p2.row((i + 1) % n2);
    if (is_strictly_inside((a + b) * 0.5, p1))
      return true;
  }

  return false;
}

double max_opening_angle_with_collisions(const kirigami::UnitPattern &pattern,
                                         bool detect_collisions) {
  double max_angle = kirigami::max_opening_angle(pattern);
  if (!detect_collisions) {
    return max_angle;
  }
  for (double t = 0.01; t < 1; t += 0.01) {
    auto deployed = pattern.deploy(max_angle * t);
    deployed =
        pattern.periodicity.norm() < 1e-3 ? deployed : deployed.replicate(2, 2);
    auto mesh = deployed.mesh.merge_close_verts().remove_unused_verts();
    for (int f0 = 0; f0 < mesh.nf(); f0++) {
      for (int f1 = f0 + 1; f1 < mesh.nf(); f1++) {
        if (do_polygons_intersect(
                mesh.V(mesh.F[f0], Eigen::placeholders::all),
                mesh.V(mesh.F[f1], Eigen::placeholders::all))) {
          return max_angle * (t - 0.01);
        }
      }
    }
  }
  return max_angle;
}

// Added for diagnostics only: the same collision predicate and the same
// deployed configuration their scan uses, but reporting WHICH face pair first
// overlaps at a given absolute angle.  do_polygons_intersect above is untouched.
bool first_collision_at(const kirigami::UnitPattern &pattern, double angle,
                        int *f0out, int *f1out) {
  auto deployed = pattern.deploy(angle);
  deployed =
      pattern.periodicity.norm() < 1e-3 ? deployed : deployed.replicate(2, 2);
  auto mesh = deployed.mesh.merge_close_verts().remove_unused_verts();
  for (int f0 = 0; f0 < mesh.nf(); f0++)
    for (int f1 = f0 + 1; f1 < mesh.nf(); f1++)
      if (do_polygons_intersect(mesh.V(mesh.F[f0], Eigen::placeholders::all),
                                mesh.V(mesh.F[f1], Eigen::placeholders::all))) {
        if (f0out) *f0out = f0;
        if (f1out) *f1out = f1;
        return true;
      }
  return false;
}
