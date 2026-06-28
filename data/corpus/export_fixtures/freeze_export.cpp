// freeze_export.cpp -- freezes the inputs and byte outputs of tests/test_export.cpp so
// the Julia port (test/test_export.jl) can compare against the C++ reference bytes.
// Linked against libkiri_export.a + libkiri_core.a from code/build.
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <nlohmann/json.hpp>

#include "core/collision.hpp"
#include "core/cut.hpp"
#include "core/generators.hpp"
#include "core/mesh.hpp"
#include "core/orientation.hpp"
#include "export/layout.hpp"
#include "export/solid.hpp"
#include "export/stl.hpp"
#include "export/svg.hpp"
#include "export/threemf.hpp"

using namespace kiri;
using namespace kiri::export_;
using nlohmann::json;

// verbatim copy of the test's fixture
static Mesh squares_patch(int nx = 3, int ny = 3, bool make_split = false) {
  Mesh m;
  auto vid = [&](int i, int j) { return j * (nx + 1) + i; };
  for (int j = 0; j <= ny; ++j)
    for (int i = 0; i <= nx; ++i) m.X.emplace_back(i, j);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      m.faces.push_back({vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)});
  m.build_topology();
  normalize_face_ccw(m);
  m.sigma.assign(m.n_faces(), -1);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) m.sigma[j * nx + i] = ((i + j) % 2 == 0) ? 1 : -1;
  if (make_split) m.sigma[0] = m.sigma[1];
  return m;
}

static void put(const std::string& path, const std::string& s) {
  std::ofstream f(path, std::ios::binary);
  f.write(s.data(), static_cast<std::streamsize>(s.size()));
}

static json layout_json(const Layout& L) {
  json j;
  j["width"] = L.width();
  j["height"] = L.height();
  j["n_hinges"] = L.hinges.size();
  j["warnings"] = L.warnings;
  json pieces = json::array();
  for (const FacePiece& fp : L.pieces) {
    json p;
    p["face"] = fp.face;
    p["n_outline"] = fp.outline.size();
    p["area_raw"] = polygon_area(fp.raw);
    p["area_outline"] = polygon_area(fp.outline);
    p["n_pin_holes"] = fp.pin_holes.size();
    p["z0"] = fp.z0;
    p["z1"] = fp.z1;
    pieces.push_back(p);
  }
  j["pieces"] = pieces;
  json hs = json::array();
  for (const HingeSite& h : L.hinges)
    hs.push_back({{"edge", h.edge}, {"neck", h.neck}, {"setback", h.setback},
                  {"p", {h.p.x(), h.p.y()}}});
  j["hinges"] = hs;
  return j;
}

static json solid_json(const TriMesh& S, const ManifoldReport& r, const std::vector<std::string>& w) {
  json j;
  j["n_tris"] = S.n_tris();
  j["n_vertices"] = S.V.size();
  j["closed"] = r.closed;
  j["consistently_oriented"] = r.consistently_oriented;
  j["n_boundary_edges"] = r.n_boundary_edges;
  j["n_nonmanifold_edges"] = r.n_nonmanifold_edges;
  j["n_flipped_edges"] = r.n_flipped_edges;
  j["n_degenerate"] = r.n_degenerate;
  j["n_components"] = r.n_components;
  j["volume"] = r.volume;
  j["warnings"] = w;
  j["n_objects"] = split_components(S).size();
  return j;
}

int main(int argc, char** argv) {
  const std::string dir = argc > 1 ? argv[1] : ".";
  json idx;
  idx["provenance"] = {
      {"source", "tests/test_export.cpp, replayed by freeze_export.cpp against code/build "
                 "libkiri_export.a + libkiri_core.a (Apple clang arm64, -O2, FMA contraction on)"},
      {"squares_patch", "3x3 unit squares, vertices (i,j) row-major, faces [v(i,j),v(i+1,j),v(i+1,j+1),v(i,j+1)], "
                        "sigma = +1 iff (i+j) even; make_split: sigma[0] = sigma[1]"},
      {"t34312", "tiling_3_4_3_12(disk(Vec2(0,0), 2.5)); sigma = assign_orientation_relaxation(m, std::mt19937(12345)).sigma (defaults 12, 800, 180)"}};

  // squares 3x3 (checkerboard) and its split variant
  {
    Mesh m = squares_patch(3, 3);
    save_mesh_json(m, dir + "/squares_3x3.json");
    Mesh ms = squares_patch(3, 3, true);
    save_mesh_json(ms, dir + "/squares_3x3_split.json");
    const CutStructure c = make_cut(m);
    const CutStructure cs = make_cut(ms);

    // SVG case (felt, theta 0, scale 10, default options)
    const Layout L = build_layout(m, c, m.X, 0.0, 10.0, profile_felt_laser());
    put(dir + "/squares_felt_theta0.svg", svg_string(L, SvgOptions{}));
    idx["squares_felt_theta0"] = layout_json(L);
    const Layout Ls = build_layout(ms, cs, ms.X, 0.0, 10.0, profile_felt_laser());
    idx["squares_split_felt_theta0"] = layout_json(Ls);
    put(dir + "/squares_split_felt_theta0.svg", svg_string(Ls, SvgOptions{}));

    // deployed layout case
    const ThetaMaxReport tm = theta_max(c, m.X);
    idx["squares_theta_max_geometric"] = tm.theta_max_geometric;
    const double th = 0.8 * tm.theta_max_geometric;
    const Layout L1 = build_layout(m, c, m.X, th, 10.0, profile_felt_laser());
    idx["squares_felt_deployed"] = layout_json(L1);
    idx["squares_felt_deployed"]["theta"] = th;
    put(dir + "/squares_felt_deployed.svg", svg_string(L1, SvgOptions{}));

    // STL case (felt)
    std::vector<std::string> w;
    const TriMesh S = build_solid(L, &w);
    idx["squares_felt_solid"] = solid_json(S, check_manifold(S, 1e-6), w);
    write_stl_binary(S, dir + "/squares_felt.stl");
    write_3mf(S, dir + "/squares_felt.3mf");

    // pin-pad case (pla)
    const Layout Lp = build_layout(m, c, m.X, 0.0, 10.0, profile_pla_print());
    idx["squares_pla_theta0"] = layout_json(Lp);
    put(dir + "/squares_pla_theta0.svg", svg_string(Lp, SvgOptions{}));
    std::vector<std::string> wp;
    const TriMesh Sp = build_solid(Lp, &wp);
    idx["squares_pla_solid"] = solid_json(Sp, check_manifold(Sp, 1e-6), wp);
    write_stl_binary(Sp, dir + "/squares_pla.stl");
    write_3mf(Sp, dir + "/squares_pla.3mf");
  }

  // 3.4.3.12 patch with split cuts
  {
    Mesh m = tiling_3_4_3_12(disk(Vec2(0, 0), 2.5));
    std::mt19937 rng(12345);
    m.sigma = assign_orientation_relaxation(m, rng).sigma;
    save_mesh_json(m, dir + "/t34312_disk2.5.json");
    const CutStructure c = make_cut(m);
    const ThetaMaxReport tm = theta_max(c, m.X);
    idx["t34312_theta_max_geometric"] = tm.theta_max_geometric;
    idx["t34312_n_hinge"] = c.n_hinge();
    idx["t34312_n_split"] = c.n_split();
    int k = 0;
    for (double frac : {0.0, 0.8}) {
      const Layout L = build_layout(m, c, m.X, frac * tm.theta_max_geometric, 10.0,
                                    profile_felt_laser());
      const std::string tag = "t34312_frac" + std::to_string(k++);
      put(dir + "/" + tag + ".svg", svg_string(L));
      idx[tag] = layout_json(L);
      idx[tag]["theta"] = frac * tm.theta_max_geometric;
    }
  }
  put(dir + "/index.json", idx.dump(1) + "\n");
  std::cout << "ok\n";
}
