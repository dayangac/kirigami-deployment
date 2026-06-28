# core/rank_checks.jl -- three structural measurements on the hole-constraint matrix L
# that the paper never states and that the rest of the pipeline never needed:
#
#   1. the row sum  r = 1^T L  and its relation to the hinge in/out degrees;
#   2. the factorization  L = R * D  (R the H x |E_hinge| preimage indicator,
#      D the |E_hinge| x N signed hinge incidence) together with the left null
#      space Z of L and the out-harmonic characterization of its elements;
#   3. the hinge graph Gamma = (F, E_hinge) and the Euler-type count
#      H = |E_hinge| - |F| + c(Gamma).
#
# Everything here is measured; where an identity fails the report carries the
# counterexample data instead of the routine asserting.
#
# Port of code/src/core/rank_checks.{hpp,cpp}. 1-based; the C++ "-1 = none" sentinels
# of the first_* fields are 0 here.

# ---------------------------------------------------------------- check 1 ----

mutable struct RowSumReport
    r::Vector{Float64}            # 1^T L, length N
    N::Int
    H::Int
    all_vertices_interior::Bool   # true iff M has no boundary vertex
    # The identity as first conjectured: r_v == hinge_in(v) - hinge_out(v).
    # MEASURED FALSE whenever a hinge edge belongs to a boundary-touching preimage
    # (a notch), because such preimages contribute no row to L.
    degree_identity::Bool
    n_degree_mismatch::Int
    first_degree_mismatch::Int
    # The corrected identity, with I(v) = 1 iff the preimage K(v) is a hole row:
    #   r_v == I(v) * hinge_in(v) - sum_{hinge v->w} I(w).
    # (All in-edges of v lie in K(v), so the in-degree is all-or-nothing; the
    # out-edges are distributed over the preimages of their targets.)
    restricted_degree_identity::Bool
    n_restricted_mismatch::Int
    first_restricted_mismatch::Int
    n_hinge::Int                  # |E_hinge|
    n_hinge_in_L::Int             # hinge edges belonging to an all-interior preimage
    no_notch_hinges::Bool         # n_hinge_in_L == n_hinge
    # r supported on boundary vertices only?
    support_on_boundary::Bool
    n_interior_nonzero::Int
    first_interior_nonzero::Int
    max_abs_r::Float64
    r_is_zero::Bool               # r == 0 exactly (integer comparison)
end
RowSumReport() = RowSumReport(Float64[], 0, 0, false, false, 0, 0, false, 0, 0, 0, 0, false,
                              false, 0, 0, 0.0, false)

# Row index into L for every vertex: the hole row whose split-forest component
# contains v, or 0 if v belongs to no preimage or to a boundary-touching one.
function vertex_hole_row(c::CutStructure, hs::HoleSet)
    row = zeros(Int, n_vertices(c.mesh))
    for (i, hi) in enumerate(hs.interior_indices)
        for v in hs.all[hi].vertices
            row[v] = i
        end
    end
    return row
end

# hinge out-targets per source vertex
function _out_targets(c::CutStructure, N::Int)
    out = [Int[] for _ in 1:N]
    for e in c.hinge_edges
        push!(out[c.hinge_dir[e].src], c.hinge_dir[e].dst)
    end
    return out
end

"""
    check_row_sum(c, hs, sys) -> RowSumReport

L is taken from `sys` (its first n_hole_rows rows), so what is measured is the
matrix the solver actually uses.
"""
function check_row_sum(c::CutStructure, hs::HoleSet, sys::LinearSystem)
    m = c.mesh
    rep = RowSumReport()
    rep.N = n_vertices(m)
    rep.H = sys.n_hole_rows
    rep.r = zeros(rep.N)
    for i in 1:sys.n_hole_rows
        rep.r += sys.A[i, :]
    end

    rep.all_vertices_interior = !any(m.vertex_is_boundary)

    # I(v): does v belong to an all-interior preimage (a row of L)?
    row = vertex_hole_row(c, hs)
    I(v) = row[v] > 0 ? 1 : 0
    rep.n_hinge = n_hinge(c)
    for e in c.hinge_edges
        rep.n_hinge_in_L += I(c.hinge_dir[e].dst)
    end
    rep.no_notch_hinges = (rep.n_hinge_in_L == rep.n_hinge)
    out_targets = _out_targets(c, rep.N)

    rep.degree_identity = true
    rep.restricted_degree_identity = true
    rep.support_on_boundary = true
    rep.r_is_zero = true
    for v in 1:rep.N
        rv = rep.r[v]
        rep.max_abs_r = max(rep.max_abs_r, abs(rv))
        expect = Float64(c.hinge_in[v] - c.hinge_out[v])
        if rv != expect
            rep.degree_identity = false
            rep.n_degree_mismatch += 1
            rep.first_degree_mismatch == 0 && (rep.first_degree_mismatch = v)
        end
        out_in_L = 0
        for w in out_targets[v]
            out_in_L += I(w)
        end
        expect_r = Float64(I(v) * c.hinge_in[v] - out_in_L)
        if rv != expect_r
            rep.restricted_degree_identity = false
            rep.n_restricted_mismatch += 1
            rep.first_restricted_mismatch == 0 && (rep.first_restricted_mismatch = v)
        end
        if rv != 0.0
            rep.r_is_zero = false
            if !m.vertex_is_boundary[v]
                rep.support_on_boundary = false
                rep.n_interior_nonzero += 1
                rep.first_interior_nonzero == 0 && (rep.first_interior_nonzero = v)
            end
        end
    end
    return rep
end

# ---------------------------------------------------------------- check 2 ----

mutable struct FactorizationReport
    H::Int
    n_hinge::Int
    N::Int
    R::Matrix{Float64}            # H x |E_hinge|, 0/1
    D::Matrix{Float64}            # |E_hinge| x N, +1 target, -1 source
    L_equals_RD::Bool             # exact (max |L - R D| == 0)
    max_abs_diff::Float64

    Z_computed::Bool
    rank_L::Int
    dim_Z::Int                    # dim of the left null space of L
    rank_identity::Bool           # rank(L) == H - dim Z
    Z::Matrix{Float64}            # H x dim_Z, orthonormal basis of the left null space

    # For every basis vector y of Z and every interior vertex v:
    #   k_v * g(v) == sum_{hinge v->w} g(w),  g(v) = y_{K(v)}  (0 if K(v) is not a hole row)
    n_harmonic_checks::Int
    n_harmonic_violations::Int
    max_harmonic_residual::Float64
    first_bad_vertex::Int
    first_bad_z::Int
    out_harmonic::Bool
end
FactorizationReport() = FactorizationReport(0, 0, 0, zeros(0, 0), zeros(0, 0), false, 0.0,
                                            false, 0, 0, false, zeros(0, 0), 0, 0, 0.0, 0, 0, false)

"""
    check_factorization(c, hs, sys, compute_Z = true, tol = 1e-9) -> FactorizationReport

`compute_Z` gates the (dense, O(H^2 N)) left-null-space SVD.
"""
function check_factorization(c::CutStructure, hs::HoleSet, sys::LinearSystem,
                             compute_Z::Bool = true, tol::Float64 = 1e-9)
    m = c.mesh
    rep = FactorizationReport()
    rep.N = n_vertices(m)
    rep.H = sys.n_hole_rows
    rep.n_hinge = n_hinge(c)

    # D: one row per hinge edge, +1 at the target, -1 at the source.
    rep.D = zeros(rep.n_hinge, rep.N)
    hinge_slot = zeros(Int, n_edges(m))
    for (j, e) in enumerate(c.hinge_edges)
        hinge_slot[e] = j
        rep.D[j, c.hinge_dir[e].dst] += 1.0
        rep.D[j, c.hinge_dir[e].src] -= 1.0
    end

    # R: 0/1, hinge edge j belongs to the preimage of hole row i.
    rep.R = zeros(rep.H, rep.n_hinge)
    for i in 1:rep.H
        for e in hs.all[hs.interior_indices[i]].edges
            hinge_slot[e] > 0 && (rep.R[i, hinge_slot[e]] = 1.0)
        end
    end

    L = sys.A[1:rep.H, :]
    RD = rep.R * rep.D
    rep.max_abs_diff = rep.H > 0 ? maximum(abs, L - RD) : 0.0
    rep.L_equals_RD = (rep.max_abs_diff == 0.0)

    (!compute_Z || rep.H == 0) && return rep

    # Left null space of L = null space of L^T: the trailing columns of the full U.
    F = svd(L; full = true)
    sv = F.S
    svtol = _sv_rank_tol(sv, 1e-10, rep.H, rep.N)
    rep.rank_L = count(s -> s > svtol, sv)
    rep.dim_Z = rep.H - rep.rank_L
    rep.Z = F.U[:, rep.H-rep.dim_Z+1:rep.H]
    rep.rank_identity = (rep.rank_L == rep.H - rep.dim_Z)
    rep.Z_computed = true

    # Out-harmonicity of g(v) = y_{K(v)} at every interior vertex.
    row = vertex_hole_row(c, hs)
    out_targets = _out_targets(c, rep.N)

    rep.out_harmonic = true
    for k in 1:rep.dim_Z
        y = rep.Z[:, k]
        g(v) = row[v] > 0 ? y[row[v]] : 0.0
        for v in 1:rep.N
            m.vertex_is_boundary[v] && continue
            s = 0.0
            for w in out_targets[v]
                s += g(w)
            end
            res = abs(Float64(c.hinge_out[v]) * g(v) - s)
            rep.n_harmonic_checks += 1
            rep.max_harmonic_residual = max(rep.max_harmonic_residual, res)
            if res > tol
                rep.out_harmonic = false
                rep.n_harmonic_violations += 1
                if rep.first_bad_vertex == 0
                    rep.first_bad_vertex = v
                    rep.first_bad_z = k
                end
            end
        end
    end
    return rep
end

# ---------------------------------------------------------------- check 3 ----

mutable struct HingeGraphReport
    n_faces::Int
    n_hinge::Int
    c_gamma::Int                  # connected components of Gamma = (F, E_hinge)
    H::Int                        # all-interior preimages only (holes)
    H_all::Int                    # every preimage, boundary-touching included
    n_notches::Int                # H_all - H
    predicted::Int                # |E_hinge| - |F| + c(Gamma)
    identity_holds::Bool          # H == predicted
    identity_with_notches::Bool   # H_all == predicted
    defect::Int                   # predicted - H
end

"""
Connected components of Gamma = (nodes = faces, edges = hinge-cut edges).
Faces touched by no hinge edge count as their own component.
"""
function hinge_graph_components(c::CutStructure)
    m = c.mesh
    dsu = DSU(n_faces(m))
    for e in c.hinge_edges
        ed = m.edges[e]
        unite!(dsu, m.half_edges[ed.he[1]].face, m.half_edges[ed.he[2]].face)
    end
    return count(f -> find!(dsu, f) == f, 1:n_faces(m))
end

function check_hinge_graph(c::CutStructure, hs::HoleSet)
    nF = n_faces(c.mesh)
    nh = n_hinge(c)
    c_gamma = hinge_graph_components(c)
    H = n_interior_holes(hs)
    H_all = length(hs.all)
    predicted = nh - nF + c_gamma
    return HingeGraphReport(nF, nh, c_gamma, H, H_all, H_all - H, predicted,
                            H == predicted, H_all == predicted, predicted - H)
end
