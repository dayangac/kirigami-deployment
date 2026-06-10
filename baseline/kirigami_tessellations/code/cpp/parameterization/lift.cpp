#include "lift.h"
#include "../Hmesh.h"
#include "../conversions.h"
#include "../kiri_mesh.h"
#include "../state.h"
#include "../utils.h"
#include <igl/AABB.h>

namespace param {

static Eigen::RowVector3d find_bary(const Eigen::MatrixXd &V,
                                    const Eigen::MatrixXi &F, int f,
                                    const Eigen::RowVector2d &p) {
  Eigen::RowVector3d bary;
  Eigen::RowVector2d v0 = V.row(F(f, 0));
  Eigen::RowVector2d v1 = V.row(F(f, 1));
  Eigen::RowVector2d v2 = V.row(F(f, 2));
  Eigen::RowVector2d e0 = v1 - v0;
  Eigen::RowVector2d e1 = v2 - v0;
  Eigen::RowVector2d e2 = p - v0;
  double d00 = e0.dot(e0);
  double d01 = e0.dot(e1);
  double d11 = e1.dot(e1);
  double d20 = e2.dot(e0);
  double d21 = e2.dot(e1);
  double denom = d00 * d11 - d01 * d01;
  bary(1) = (d11 * d20 - d01 * d21) / denom;
  bary(2) = (d00 * d21 - d01 * d20) / denom;
  bary(0) = 1.0 - bary(1) - bary(2);
  return bary;
}

igl::AABB<Eigen::MatrixXd, 2> tree;
Eigen::MatrixXd cached_UV;

utils::Hmesh lift(utils::Hmesh &mesh, const Eigen::MatrixXd &UV,
                  const geom::SymPattern &pattern, double opening_angle,
                  double scale, const Eigen::RowVector2d &offset, double rotate,
                  bool remove_dangling_faces) {
  // Shift UV.
  Eigen::MatrixXi F = convert::to_eig_mat(mesh.F);
  if (cached_UV.rows() != UV.rows() || cached_UV.cols() != UV.cols() ||
      (cached_UV - UV).norm() > 1e-4) {
    tree = igl::AABB<Eigen::MatrixXd, 2>();
    tree.init(UV, F);
    cached_UV = UV;
  }

  geom::KiriMesh kiri_mesh(pattern.mesh(2, 2));
  kiri_mesh.generate_cut_mesh();
  kiri_mesh.generate_rotating_info();
  Eigen::MatrixXd closed_v = kiri_mesh.open(0);
  Eigen::MatrixXd opened_v = kiri_mesh.open(opening_angle);
  int nv = opened_v.rows(), pattern_nf = pattern.faces.size();
  int v0 = kiri_mesh.cut_mesh.F[0][0];
  closed_v.rowwise() -= (closed_v.row(v0)).eval();
  opened_v.rowwise() -= (opened_v.row(v0)).eval();
  Eigen::Matrix2d rot_mat = Eigen::Rotation2Dd(rotate).toRotationMatrix();
  closed_v = ((closed_v * scale) * rot_mat).eval();
  opened_v = ((opened_v * scale) * rot_mat).eval();

  // Find bounds.
  Eigen::Matrix2d closed_basis;
  closed_basis << closed_v.row(kiri_mesh.cut_mesh.F[pattern_nf * 2][0]),
      closed_v.row(kiri_mesh.cut_mesh.F[pattern_nf][0]);

  Eigen::Matrix2d basis;
  basis << opened_v.row(kiri_mesh.cut_mesh.F[pattern_nf * 2][0]),
      opened_v.row(kiri_mesh.cut_mesh.F[pattern_nf][0]);
  
  Eigen::MatrixXd UV_in_basis =
      (cached_UV.rowwise() + offset) * basis.inverse();
  Eigen::RowVector2i min_rep =
      UV_in_basis.colwise().minCoeff().array().floor().cast<int>();
  Eigen::RowVector2i max_rep =
      UV_in_basis.colwise().maxCoeff().array().ceil().cast<int>();
  min_rep -= Eigen::RowVector2i(1, 1);
  max_rep += Eigen::RowVector2i(1, 1);

  // Add vertices and faces where all vertices are inside the UV.
  std::vector<Eigen::Vector3d> verts;
  std::vector<std::vector<int>> faces;
  std::vector<Eigen::Vector2d> ground_verts;
  std::vector<Eigen::Vector2d> ground_closed_verts;
  int cur_ind = 0;
  for (int i = min_rep.x(); i <= max_rep.x(); i++) {
    for (int j = min_rep.y(); j <= max_rep.y(); j++) {
      std::vector<bool> is_inside(nv, false);
      // Pullback vertices.
      for (int v = 0; v < nv; v++) {
        Eigen::RowVector2d q = opened_v.row(v);
        Eigen::RowVectorXd p = q + Eigen::RowVector2d(i, j) * basis;
        std::vector<int> f = tree.find(cached_UV, F, p);
        if (!f.empty()) {
          is_inside[v] = true;
          auto bary = find_bary(cached_UV, F, f[0], p);
          verts.push_back(bary * mesh.V(F.row(f[0]), Eigen::all));
        } else {
          verts.push_back(Eigen::Vector3d::Zero());
        }
        ground_verts.push_back(p);
        ground_closed_verts.push_back(closed_v.row(v) +
                                      Eigen::RowVector2d(i, j) * closed_basis);
      }
      // Only add faces for which all verts are inside.
      for (int f = 0; f < pattern_nf; f++) {
        bool all_inside = true;
        for (int v : kiri_mesh.cut_mesh.F[f]) {
          all_inside &= is_inside[v];
        }
        if (all_inside) {
          std::vector<int> face = kiri_mesh.cut_mesh.F[f];
          for (auto &v : face) {
            v += cur_ind;
          }
          faces.push_back(face);
        }
      }
      cur_ind += nv;
    }
  }
  if (faces.empty()) {
    return utils::Hmesh();
  }

  // Remove faces dangling by a single vertex or edge.
  if (remove_dangling_faces) {
    int last_size_faces = faces.size();
    do {
      last_size_faces = faces.size();
      state::lifted = utils::Hmesh(convert::to_eig_mat(verts), faces)
                          .remove_unused_verts()
                          .merge_close_verts();
      std::vector<std::vector<int>> keep_faces;
      for (int i = 0; i < faces.size(); i++) {
        int num_non_boundary_edges = 0;
        auto e = state::lifted.face(i).edge();
        do {
          if (e->twin()) {
            num_non_boundary_edges++;
          }
          e = e->next();
        } while (e != state::lifted.face(i).edge());
        if (num_non_boundary_edges >= 2) {
          keep_faces.push_back(faces[i]);
        }
      }
      faces = keep_faces;
    } while (last_size_faces != faces.size());
  }

  state::lifted = utils::Hmesh(convert::to_eig_mat(verts), faces)
                      .remove_unused_verts()
                      .merge_close_verts();
  state::ground = utils::Hmesh(convert::to_eig_mat(ground_verts), faces)
                      .remove_unused_verts()
                      .merge_close_verts();
  state::ground_closed =
      utils::Hmesh(convert::to_eig_mat(ground_closed_verts), faces)
          .remove_unused_verts()
          .merge_close_verts();

  return state::lifted;
}

static bool is_rotated_around_origin(geom::KiriMesh &inp, int f) {
  utils::Hmesh::Edge *e = inp.mesh.face(f).edge();
  do {
    if (e->twin()) {
      int fi = e->twin()->fi;
      int fvi = e->fvi, fivi = e->twin()->next()->fvi;
      if (inp.cut_mesh.F[f][fvi] == inp.cut_mesh.F[fi][fivi]) {
        return true;
      } else {
        return false;
      }
    }
    e = e->next();
  } while (e != inp.mesh.face(f).edge());
  return false;
}

int get_new_vert_index(geom::KiriMesh &inp, utils::Hmesh::Edge *e, bool orig,
                       const std::vector<std::vector<int>> &F) {
  if (orig) {
    utils::Hmesh::Edge *ee = e;
    while (ee->twin()) {
      if (ee->fi < F.size())
        return F[ee->fi][ee->fvi];
      // Check twin face.
      if (ee->twin()->fi < F.size())
        return F[ee->twin()->fi][ee->twin()->next()->fvi];
      if (!ee->twin()->prev()->twin())
        break;
      // Jump to next identical vert.
      ee = ee->twin()->prev()->twin()->prev();
      if (ee->fi < F.size())
        return F[ee->fi][ee->fvi];
      if (ee == e)
        break;
    }
    ee = e;
    while (ee->next()->twin()) {
      if (ee->fi < F.size())
        return F[ee->fi][ee->fvi];
      if (ee->next()->twin()->fi < F.size()) {
        return F[ee->next()->twin()->fi]
                [ee->next()->twin()->next()->next()->fvi];
      }
      if (!ee->next()->twin()->next()->twin())
        break;
      ee = ee->next()->twin()->next()->twin();
      if (ee->fi < F.size())
        return F[ee->fi][ee->fvi];
      if (ee == e)
        break;
    }
  } else {
    if (e->prev()->twin()) {
      return get_new_vert_index(inp, e->prev()->twin(), true, F);
    } else if (e->prev()->prev()->twin()) {
      return get_new_vert_index(inp, e->prev()->prev()->twin()->prev(), true,
                                F);
    }
    // } else if (e->next()->twin()) {
    //   return get_new_vert_index(inp, e->next()->twin()->prev(), true, F);
  }
  return -1;
}

utils::Hmesh flip_combinatorics(geom::KiriMesh &inp) {
  std::vector<Eigen::VectorXd> V;
  std::vector<std::vector<int>> F;
  for (int i = 0; i < inp.mesh.nf(); i++) {
    // Build face, adding new vertices if needed.
    std::vector<int> face;
    bool orig = is_rotated_around_origin(inp, i);
    utils::Hmesh::Edge *e = inp.mesh.face(i).edge();
    do {
      int vert = get_new_vert_index(inp, e, orig, F);
      if (vert == -1) {
        V.push_back(inp.mesh.V.row(e->origin()->index));
        vert = V.size() - 1;
      }
      face.push_back(vert);
      e = e->next();
    } while (e != inp.mesh.face(i).edge());
    F.push_back(face);
  }
  return utils::Hmesh(convert::to_eig_mat(V), F);
}

} // namespace param