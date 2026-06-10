#pragma once

#include "../Hmesh.h"
#include "../sym_pattern.h"

namespace param {

utils::Hmesh lift(utils::Hmesh &mesh, const Eigen::MatrixXd &UV,
                  const geom::SymPattern &pattern, double opening_angle = 0.0,
                  double scale = 1.0,
                  const Eigen::RowVector2d &offset = Eigen::RowVector2d::Zero(),
                  double rotate = 0.0, bool remove_dangling_faces = false);


} // namespace param