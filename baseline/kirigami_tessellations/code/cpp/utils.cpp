#include "utils.h"
#include "Hmesh.h"
#include "conversions.h"
#include <iostream>
#include <fstream>

namespace utils {

void export_obj(const Eigen::MatrixXd &verts,
                const std::vector<std::vector<int>> &F,
                const std::string &file_name) {
  std::ofstream file(file_name);
  for (int i = 0; i < verts.rows(); i++) {
    if (verts.cols() == 2) {
      file << "v " << verts(i, 0) << " " << verts(i, 1) << " 0" << std::endl;
    } else {
      file << "v " << verts(i, 0) << " " << verts(i, 1) << " " << verts(i, 2)
           << std::endl;
    }
  }
  for (int i = 0; i < F.size(); i++) {
    file << "f";
    for (int j = 0; j < F[i].size(); j++) {
      file << " " << F[i][j] + 1;
    }
    file << std::endl;
  }
  file.close();
}

std::pair<int, int> approximate_ratio(double x, int maxDenominator) {
    double frac = x;
    int a = static_cast<int>(frac);
    int h1 = 1, h0 = 0;
    int k1 = 0, k0 = 1;
    
    int numerator = a;
    int denominator = 1;

    int count = 0;
    while (count < 100) {
        double f = frac - a;
        if (f < 1e-6) break; // good enough

        frac = 1.0 / f;
        a = static_cast<int>(frac);

        int h2 = a * h1 + h0;
        int k2 = a * k1 + k0;

        if (k2 > maxDenominator) break;

        h0 = h1; h1 = h2;
        k0 = k1; k1 = k2;

    denominator = h2;
    numerator = k2;
    count++;
  }
  if (count == 100) {
    std::cerr << "Failed to approximate ratio for " << x << std::endl;
    return std::make_pair(1000, 1000);
  }
  return std::make_pair(numerator, denominator);
}

} // namespace utils