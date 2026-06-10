#pragma once
#include "Hmesh.h"
#include <Eigen/Eigen>

namespace utils {

template <typename T>
std::vector<T> slice(const std::vector<T> &src, const std::vector<int> &inds) {
  std::vector<T> res;
  for (int i = 0; i < inds.size(); i++) {
    res.push_back(src[inds[i]]);
  }
  return res;
}

template <typename T, typename G>
std::vector<G> map(const std::vector<T> &src, std::function<G(T)> f) {
  std::vector<G> res;
  for (int i = 0; i < src.size(); i++) {
    res.push_back(f(src[i]));
  }
  return res;
}

void export_obj(const Eigen::MatrixXd &verts,
                const std::vector<std::vector<int>> &F,
                const std::string &file_name);

inline void export_obj(utils::Hmesh &mesh, const std::string &file_name) {
  export_obj(mesh.V, mesh.F, file_name);
}

std::pair<int, int> approximate_ratio(double x, int maxDenominator = 1000);

template <int dim, typename... Args>
Eigen::MatrixXd matcat(const Args &...args) {
  int other_dim = 0, this_dim = 0;
  if constexpr (dim == 0) {
    (
        [&] {
          other_dim += args.rows();
          this_dim = args.cols();
        }(),
        ...);
    Eigen::MatrixXd res(other_dim, this_dim);
    int cur = 0;
    (
        [&] {
          res.block(cur, 0, args.rows(), args.cols()) = args;
          cur += args.rows();
        }(),
        ...);
    return res;
  } else {
    (
        [&] {
          other_dim += args.cols();
          this_dim = args.rows();
        }(),
        ...);
    Eigen::MatrixXd res(this_dim, other_dim);
    int cur = 0;
    (
        [&] {
          res.block(0, cur, args.rows(), args.cols()) = args;
          cur += args.cols();
        }(),
        ...);
    return res;
  }
}

} // namespace utils