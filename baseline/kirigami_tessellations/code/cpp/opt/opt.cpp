#include "opt.h"
#include "../Hmesh.h"
#include "../conversions.h"
#include "../state.h"
#include <Optiz/NewtonSolver/Common.h>
#include <Optiz/NewtonSolver/Problem.h>
#include <Optiz/Problem.h>
#include <cstddef>
#include <igl/AABB.h>
#include <memory>

namespace opt {

int num_iters = 0;
float orig_close_to_init_weight = 10;
float rigid_weight = 10.0, closeness_weight = 0.1, planarity_weight = 10,
      close_to_init_weight = 0.5, smoothness_weight = 0.1;

igl::AABB<Eigen::MatrixXd, 3> tree;
Eigen::MatrixXi F;

std::unique_ptr<Optiz::Problem> prob = nullptr;

struct FakeFactory {
  using Scalar = double;
  Eigen::MatrixXd &var_block(int i) {
    if (i == 0) {
      return state::opt_lifted;
    }
    return state::opt_ground;
  }
};

auto get_rigid_error(int i, auto &x, bool relative = false,
                     bool sqrd_err = true) {
  auto &e = state::lifted.edges[i];
  auto &x1 = x.var_block(0);
  auto &x2 = x.var_block(1);
  int f = e.fi, fvi = e.fvi, fvj = e.next()->fvi;
  // in R3.
  auto v1 = x1.row(e.vi), v2 = x1.row(e.next()->vi);
  // in R2.
  auto v3 = x2.row(state::ground_closed.F[f][fvi]),
       v4 = x2.row(state::ground_closed.F[f][fvj]);
  if (!relative) {
    return Optiz::sqr((v1 - v2).norm() - (v3 - v4).norm());
  } else {
    double l1 =
        0.5 * (Optiz::val((v1 - v2).norm()) + Optiz::val((v3 - v4).norm()));
    if (sqrd_err) {
      return Optiz::sqr((v1 - v2).norm() - (v3 - v4).norm()) / (l1 * l1);
    }
    return abs((v1 - v2).norm() - (v3 - v4).norm()) / l1;
  }
}

auto get_face_rigid_error(int f, auto &x, bool relative = false,
                          bool sqrd_err = true) {
  using T = FACTORY_TYPE(x);
  auto &x1 = x.var_block(0);
  auto &x2 = x.var_block(1);
  T res(0.0);
  for (int i = 0; i < state::lifted.F[f].size(); i++) {
    for (int j = i + 2; j < state::lifted.F[f].size(); j++) {
      if (i == 0 && j == state::lifted.F[f].size() - 1)
        continue; // Already covered - it's a polygon edge.
      auto v1 = x1.row(state::lifted.F[f][i]);
      auto v2 = x1.row(state::lifted.F[f][j]);
      auto g1 = x2.row(state::ground_closed.F[f][i]);
      auto g2 = x2.row(state::ground_closed.F[f][j]);
      if (!relative) {
        res += Optiz::sqr((v1 - v2).norm() - (g1 - g2).norm());
      } else {
        double l1 =
            0.5 * (Optiz::val((v1 - v2).norm()) + Optiz::val((g1 - g2).norm()));
        if (sqrd_err) {
          res += Optiz::sqr((v1 - v2).norm() - (g1 - g2).norm()) / (l1 * l1);
        } else {
          res += abs((v1 - v2).norm() - (g1 - g2).norm()) / (l1);
        }
      }
    }
  }
  return res;
}

auto get_closeness_error(int i, auto &x) {
  int f;
  Eigen::RowVector3d cp;
  Eigen::RowVector3d p = x.var_block(0).row(i);
  tree.squared_distance(state::target_mesh.V, F, p, f, cp);
  Eigen::Vector3d n = state::target_mesh.face(f).normal();
  return std::abs((p - cp).dot(n));
};

Eigen::Vector3d get_plane_normal(const Eigen::MatrixXd &VF) {
  Eigen::MatrixXd centered = VF;
  centered.rowwise() -= centered.colwise().mean();
  Eigen::MatrixXd cov = centered.transpose() * centered;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
  return eig.eigenvectors().col(0).normalized();
}

double get_planarity_error(int i, auto &xx) {
  if (state::lifted.F[i].size() <= 3)
    return 0.0;
  auto &x = xx.var_block(0);
  Eigen::Vector3d n = get_plane_normal(x(state::lifted.F[i], Eigen::all));
  double err(0.0);
  for (int j = 0; j < state::lifted.F[i].size() - 1; j++) {
    auto vec = x.row(state::lifted.F[i][j + 1]) - x.row(state::lifted.F[i][j]);
    err +=
        std::abs(M_PI / 2 - std::acos(n.dot(vec.normalized()))) / (M_PI) * 180;
  }
  return err;
};

void init() {
  F = convert::to_eig_mat(state::target_mesh.F);
  Eigen::MatrixXd &V = state::target_mesh.V;
  tree.init(V, F);
  state::opt_lifted = state::lifted.V;
  state::opt_ground = state::ground_closed.V;
}

void init_prob() {
  prob = std::make_unique<Optiz::Problem>(state::opt_lifted, state::opt_ground);
  prob->options().set_line_search_iters(100);
  prob->options().set_iters(1);
  prob->options().set_report_level(Optiz::Problem::Options::NONE);

  // Rigid energy.
  prob->add_element_energy(state::lifted.edges.size(), [&](int i, auto &x) {
    return rigid_weight * get_rigid_error(i, x, true);
  });

  // Closeness energy.
  prob->add_element_energy(state::lifted.V.rows(), [&](int i, auto &x) {
    using T = FACTORY_TYPE(x);
    int f;
    Eigen::RowVector3d cp;
    Eigen::RowVector3d p = prob->x(0).row(i);
    tree.squared_distance(state::target_mesh.V, F, p, f, cp);
    auto &x1 = x.var_block(0);
    auto xp = x1.row(i);
    Eigen::Vector3d n = state::target_mesh.face(f).normal();
    return closeness_weight * Optiz::sqr((xp - cp).dot(n));
  });
  // Close to the init.
  prob->add_element_energy(state::lifted.V.rows(), [&](int i, auto &x) {
    auto orig = state::lifted.V.row(i);
    auto &x1 = x.var_block(0);
    auto xp = x1.row(i);
    return close_to_init_weight * (xp - orig).squaredNorm();
  });

  // Planarity.
  if (state::lifted.max_face_degree() > 3) {
    prob->add_element_energy(state::lifted.F.size(), [&](int i, auto &x) {
      using T = FACTORY_TYPE(x);
      if (state::lifted.F[i].size() <= 3)
        return T(0.0);
      Eigen::Vector3d n =
          get_plane_normal(prob->x(0)(state::lifted.F[i], Eigen::all));
      T err(0.0);
      auto &x1 = x.var_block(0);
      for (int j = 0; j < state::lifted.F[i].size() - 1; j++) {
        for (int k = j + 1; k < state::lifted.F[i].size(); k++) {
          Eigen::RowVector3<T> vec =
              x1.row(state::lifted.F[i][k]) - x1.row(state::lifted.F[i][j]);
          double vec_val = Optiz::val(vec.norm());
          err += Optiz::sqr(n.dot(vec) / vec_val);
        }
      }
      return planarity_weight * err;
    });
    prob->add_element_energy(state::lifted.F.size(), [&](int i, auto &x) {
      return rigid_weight * get_face_rigid_error(i, x, true);
    });
  }
}

void reset_optimization() {
  state::opt_lifted = state::lifted.V;
  state::opt_ground = state::ground_closed.V;
  close_to_init_weight = orig_close_to_init_weight;
  num_iters = 0;
  prob = nullptr;
  init_prob();
}

void optimize_rigidity() {
  if (F.rows() == 0) {
    init();
  }

  if (!prob) {
    init_prob();
  }
  prob->optimize();
  state::opt_lifted = prob->x(0);
  state::opt_ground = prob->x(1);
  if (num_iters++ % 20 == 0)
    close_to_init_weight *= 0.8;
}

void get_errors(double *rigid_avg, double *rigid_max, double *close_avg,
                double *close_max, double *planarity_avg,
                double *planarity_max) {
  if (F.rows() == 0) {
    init();
  }
  auto factory = FakeFactory();
  Eigen::VectorXd rigid_errs(state::lifted.edges.size());
  for (int i = 0; i < state::lifted.edges.size(); i++) {
    rigid_errs(i) = (get_rigid_error(i, factory, true, false));
  }
  if (state::lifted.max_face_degree() > 3) {
    rigid_errs.conservativeResize(rigid_errs.size() + state::lifted.F.size());
    for (int i = 0; i < state::lifted.F.size(); i++) {
      rigid_errs(state::lifted.edges.size() + i) =
          (get_face_rigid_error(i, factory, true, false));
    }
  }
  if (rigid_avg)
    *rigid_avg = rigid_errs.mean();
  if (rigid_max)
    *rigid_max = rigid_errs.maxCoeff();

  Eigen::VectorXd close_errs(state::lifted.V.rows());
  for (int i = 0; i < state::lifted.V.rows(); i++) {
    close_errs(i) = get_closeness_error(i, factory);
  }
  if (close_avg)
    *close_avg = close_errs.mean();
  if (close_max)
    *close_max = close_errs.maxCoeff();

  if (planarity_avg && planarity_max) {
    std::vector<double> errs_vec;
    for (int i = 0; i < state::lifted.F.size(); i++) {
      if (state::lifted.F[i].size() == 3)
        continue;
      errs_vec.push_back(get_planarity_error(i, factory));
    }
    if (errs_vec.empty()) {
      *planarity_avg = 0;
      *planarity_max = 0;
      return;
    } else {
      Eigen::Map<Eigen::VectorXd> planarity_errs(errs_vec.data(),
                                                 errs_vec.size());
      *planarity_avg = planarity_errs.mean();
      *planarity_max = planarity_errs.maxCoeff();
    }
  }
}

} // namespace opt