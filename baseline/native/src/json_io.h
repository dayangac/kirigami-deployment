#pragma once
#include "geometry/unit_pattern.h"
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

namespace io {

// Our JSON graph contract (code/README.md):
//   vertices [[x,y]], faces [[i...]], orientation [+-1] (+1 = CW), periodic {th,tv}
// The authors use face_colors in {0,1}.
// For a hinge edge their kirigami::deploy pivots about the SOURCE of the
// colour-0 face's stored half-edge, whereas our convention pivots about the
// source of the sigma = -1 face's stored half-edge; those are opposite ends of
// the same edge.  The mapping that makes the two forward-kinematics maps agree
// (measured: 2e-8 max deviation after rigid alignment) is therefore
//     sigma = +1 (clockwise)         <->  colour 1
//     sigma = -1 (counter-clockwise) <->  colour 0
// Their final normalisation rotates the deployed patch by the opposite sign to
// ours, which rigid alignment absorbs.
int sigma_to_color(int sigma);
int color_to_sigma(int color);

struct Input {
  kirigami::UnitPattern pattern;
  bool had_orientation = false;
  bool had_periodic = false;
  int n_cw_faces_reversed = 0; // faces stored CW that we flipped to CCW
};

Input load(const std::string &path);

nlohmann::json verts_to_json(const Eigen::MatrixXd &V);
nlohmann::json faces_to_json(const std::vector<std::vector<int>> &F);
nlohmann::json pattern_to_json(const kirigami::UnitPattern &p);
void write(const std::string &path, const nlohmann::json &j);

} // namespace io
