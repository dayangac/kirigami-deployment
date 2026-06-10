#pragma  once
#include <Eigen/Eigen>
#include "../Hmesh.h"

namespace param {
  Eigen::MatrixXd isometric_param(utils::Hmesh& mesh);
}