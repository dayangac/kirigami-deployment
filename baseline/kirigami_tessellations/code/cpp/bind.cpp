#include "Hmesh.h"
#include "conversions.h"
#include "kiri_mesh.h"
#include "multi_grid.h"
#include "parameterization/lift.h"
#include "parameterization/param.h"
#include "opt/opt.h"
#include "state.h"
#include "sym_pattern.h"
#include <emscripten/bind.h>
#include <iostream>
using namespace emscripten;

emscripten::val
std_vector_to_js_array(const std::vector<std::vector<int>> &vec) {
  emscripten::val res = val::global("Array").new_();
  for (int i = 0; i < vec.size(); i++) {
    emscripten::val pts = val::global("Array").new_();
    for (int j = 0; j < vec[i].size(); j++) {
      pts.call<void>("push", vec[i][j]);
    }
    res.call<void>("push", pts);
  }
  return res;
}

val eig_vec_to_js_array(const Eigen::VectorXd &vec) {
  emscripten::val res = val::global("Array").new_();
  for (int i = 0; i < vec.size(); i++)
    res.call<void>("push", vec(i));
  return res;
}

val eig_to_js_array(const Eigen::MatrixXd &mat) {
  emscripten::val res = val::global("Array").new_();
  for (int i = 0; i < mat.rows(); i++) {
    emscripten::val row = val::global("Array").new_();
    for (int j = 0; j < mat.cols(); j++) {
      row.call<void>("push", mat(i, j));
    }
    res.call<void>("push", row);
  }
  return res;
}

val createBasePattern(val groups) {
  int length = groups["length"].as<unsigned>();
  val res = val::object();

  if (length < 2) {
    res.set("error", "Not enough lines to form a unit pattern.");
    return res;
  }

  geom::MultiGrid mg;
  for (int i = 0; i < length; i++) {
    val group = groups[i];
    val points = group["points"];
    if (points["length"].as<unsigned>() < 3)
      continue;

    Eigen::Vector2d p0, p1, p2;
    p0 << points[0]["x"].as<double>() / 100, points[0]["y"].as<double>() / 100;
    p1 << points[1]["x"].as<double>() / 100, points[1]["y"].as<double>() / 100;
    p2 << points[2]["x"].as<double>() / 100, points[2]["y"].as<double>() / 100;
    mg.lines.emplace_back(p0, p1, p2);
  }

  if (mg.calculate_periodicity().x() > 100) {
    res.set("error", "No unit pattern found.");
    return res;
  }

  state::pattern = mg.get_sym_pattern();
  auto mesh = state::pattern.mesh(1, 1).merge_close_verts();

  res.set("verts", eig_to_js_array(mesh.V));
  res.set("faces", std_vector_to_js_array(mesh.F));
  return res;
}

Eigen::MatrixXd to_eig_mat(const std::vector<Eigen::VectorXd> &vecs) {
  if (vecs.empty())
    return Eigen::MatrixXd();
  Eigen::MatrixXd mat(vecs.size(), vecs[0].size());
  for (int i = 0; i < vecs.size(); ++i)
    mat.row(i) = vecs[i].transpose();
  return mat;
}

val createFinalPattern(val merged_mesh, double opening_angle, int rep_x,
                       int rep_y, bool if_cut_open = true) {
  val res = val::object();

  utils::Hmesh mesh = state::pattern.mesh(1, 1).merge_close_verts();

  geom::KiriMesh kiri_mesh(mesh);
  kiri_mesh.generate_cut_mesh();
  kiri_mesh.generate_rotating_info();

  // Create expanded mesh with repetitions before opening
  utils::Hmesh expanded_mesh =
      state::pattern.mesh(rep_x, rep_y).merge_close_verts();
  if (if_cut_open) {
    // Create kirigami mesh from expanded mesh
    geom::KiriMesh expanded_kiri(expanded_mesh);
    expanded_kiri.generate_cut_mesh();
    expanded_kiri.generate_rotating_info();

    auto opened = expanded_kiri.open_separate_faces(opening_angle);

    // Calculate rotation for better visualization
    // Eigen::Matrix2d rot_mat = Eigen::Matrix2d::Identity();
    // if (V.rows() > 0) {
    //     Eigen::Vector2d centroid = V.colwise().mean();
    //     Eigen::MatrixXd centered = V.rowwise() - centroid.transpose();
    //     Eigen::JacobiSVD<Eigen::MatrixXd> svd(centered, Eigen::ComputeThinU |
    //     Eigen::ComputeThinV); rot_mat = svd.matrixV();
    // }
    auto basis = state::pattern.get_bounding_parallelogram(opening_angle);
    Eigen::RowVector2d xx = basis.row(0);
    Eigen::Matrix2d rot_mat =
        Eigen::Rotation2Dd(std::atan2(xx.y(), xx.x())).toRotationMatrix();

    opened.V = (opened.V * rot_mat).eval();
    opened.V.rowwise() -= opened.V.colwise().mean();

    auto friendliness = expanded_kiri.get_friendly_vertices(opened);

    res.set("verts", eig_to_js_array(opened.V));
    res.set("faces", std_vector_to_js_array(opened.F));
    res.set("friendly", eig_vec_to_js_array(friendliness));

    // Calculate max opening angle using a larger mesh
    geom::KiriMesh max_angle_mesh(state::pattern.mesh());
    max_angle_mesh.generate_cut_mesh();
    max_angle_mesh.generate_rotating_info();
    double max_opening_angle = max_angle_mesh.max_opening_angle();
    res.set("max_opening_angle", max_opening_angle);

    return res;
  } else {
    // If not cutting open, just return the expanded mesh
    res.set("verts", eig_to_js_array(expanded_mesh.V));
    res.set("faces", std_vector_to_js_array(expanded_mesh.F));
    res.set("max_opening_angle", 0.0); // No opening angle if not cutting
    return res;
  }
}

val parameterizeMesh(val js_V, val js_F) {
  const Eigen::MatrixXd V = convert::js_array_to_eig<double>(js_V);
  const Eigen::MatrixXi F = convert::js_array_to_eig<int>(js_F);
  state::target_mesh = utils::Hmesh(V, convert::to_vec_vec(F));
  double scale = std::sqrt(state::target_mesh.area());
  state::target_mesh.V /= scale;
  Eigen::MatrixXd UV = param::isometric_param(state::target_mesh);
  state::UV = UV;
  state::UV.rowwise() -= state::UV.colwise().mean();
  return eig_to_js_array(UV);
}

val liftPattern() {
  if (state::UV.rows() == 0) {
    std::cout << "No UV coordinates to lift from." << std::endl;
    return val::array();
  }
  geom::KiriMesh max_angle_mesh(state::pattern.mesh());
  max_angle_mesh.generate_cut_mesh();
  max_angle_mesh.generate_rotating_info();
  double max_opening_angle = max_angle_mesh.max_opening_angle();
  state::opening_angle = max_opening_angle;
  double area_ratio = state::target_mesh.area() / max_angle_mesh.mesh.area();
  // Scale to achieve an area_ratio of 5.
  state::scale = std::sqrt(area_ratio / 5);
  utils::Hmesh lifted = param::lift(
      state::target_mesh, state::UV, state::pattern, max_opening_angle,
      state::scale, state::translate.cast<double>(), 0, false);
  opt::reset_optimization();
  return eig_to_js_array(lifted.V);
}

val updateLiftParams(val scale, val rotation) {
  double applied_scale = state::scale / scale.as<double>();
  double applied_rotation = rotation.as<double>();
  param::lift(state::target_mesh, state::UV, state::pattern,
              state::opening_angle, applied_scale,
              state::translate.cast<double>(), applied_rotation, false);
  opt::reset_optimization();
  val res = val::object();
  Eigen::MatrixXd centered_opt_lifted = state::opt_lifted;
  centered_opt_lifted.rowwise() -= centered_opt_lifted.colwise().mean();
  res.set("lifted_verts", eig_to_js_array(centered_opt_lifted));
  res.set("ground_verts", eig_to_js_array(state::opt_ground));
  res.set("lifted_faces", std_vector_to_js_array(state::lifted.F));
  res.set("ground_faces", std_vector_to_js_array(state::ground_closed.F));
  return res;
}

val initOptimization(val js_V, val js_F) {
  parameterizeMesh(js_V, js_F);
  liftPattern();
  opt::init();
  val res = val::object();
  Eigen::MatrixXd centered_opt_lifted = state::opt_lifted;
  centered_opt_lifted.rowwise() -= centered_opt_lifted.colwise().mean();
  res.set("lifted_verts", eig_to_js_array(centered_opt_lifted));
  res.set("ground_verts", eig_to_js_array(state::opt_ground));
  res.set("lifted_faces", std_vector_to_js_array(state::lifted.F));
  res.set("ground_faces", std_vector_to_js_array(state::ground_closed.F));
  return res;
}

val runOptimization() {
  opt::optimize_rigidity();
  val res = val::object();
  Eigen::MatrixXd centered_opt_lifted = state::opt_lifted;
  centered_opt_lifted.rowwise() -= centered_opt_lifted.colwise().mean();
  res.set("lifted_verts", eig_to_js_array(centered_opt_lifted));
  res.set("ground_verts", eig_to_js_array(state::opt_ground));
  res.set("lifted_faces", std_vector_to_js_array(state::lifted.F));
  res.set("ground_faces", std_vector_to_js_array(state::ground_closed.F));
  double rigid_avg, rigid_max, close_avg, close_max, planarity_avg,
      planarity_max;
  opt::get_errors(&rigid_avg, &rigid_max, &close_avg, &close_max,
                  &planarity_avg, &planarity_max);
  res.set("rigid_avg", rigid_avg);
  res.set("rigid_max", rigid_max);
  res.set("close_avg", close_avg);
  res.set("close_max", close_max);
  res.set("planarity_avg", planarity_avg);
  res.set("planarity_max", planarity_max);
  return res;
}

val get_errors() {
  val res = val::object();
  double rigid_avg, rigid_max, close_avg, close_max, planarity_avg,
      planarity_max;
  opt::get_errors(&rigid_avg, &rigid_max, &close_avg, &close_max,
                  &planarity_avg, &planarity_max);
  res.set("rigid_avg", rigid_avg);
  res.set("rigid_max", rigid_max);
  res.set("close_avg", close_avg);
  res.set("close_max", close_max);
  res.set("planarity_avg", planarity_avg);
  res.set("planarity_max", planarity_max);
  return res;
}

EMSCRIPTEN_BINDINGS(my_module) {
  function("createBasePattern", &createBasePattern);
  function("createFinalPattern",
           optional_override([](val merged_mesh, double opening_angle,
                                int rep_x, int rep_y, bool if_cut_open) {
             return createFinalPattern(merged_mesh, opening_angle, rep_x, rep_y,
                                       if_cut_open);
           }));

  function("createFinalPattern",
           optional_override(
               [](val merged_mesh, double opening_angle, int rep_x, int rep_y) {
                 return createFinalPattern(merged_mesh, opening_angle, rep_x,
                                           rep_y, true);
               }));
  function("parameterizeMesh", &parameterizeMesh);
  function("liftPattern", &liftPattern);
  function("initOptimization", &initOptimization);
  function("updateLiftParams", &updateLiftParams);
  function("runOptimization", &runOptimization);
  function("get_errors", &get_errors);
}