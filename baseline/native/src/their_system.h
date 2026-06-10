#pragma once
#include "geometry/unit_pattern.h"
#include <Eigen/Sparse>
#include <vector>

namespace their {

// A line-by-line transcription of steps 1-5 of kirigami::make_deployable
// (baseline/tuttekiri/code/cpp/geometry/deployment.cpp) for the
// map_boundary_to_unit_disk == false path.  make_deployable does not return its
// constraint matrix, so this rebuild exists only so that residuals of OTHER
// embeddings can be measured in the authors' own system.  It is validated at
// run time by checking that the authors' own returned solution satisfies it and
// that dim ker A equals the kernel they return; `tuttekiri_cli` refuses to
// report cross-residuals if either check fails.
struct System {
  Eigen::SparseMatrix<double> A; // (num_holes + num_constraints) x nv
  Eigen::MatrixXd b;             // rows(A) x 2
  int num_holes = 0;
  int num_constraints = 0;
  std::vector<int> vertex_to_hole;
};

// `pattern` must already have had make_periodic() applied if it is periodic,
// exactly as make_deployable does before assembling.
System build(const kirigami::UnitPattern &pattern);

} // namespace their
