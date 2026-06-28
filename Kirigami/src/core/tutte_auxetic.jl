# core/tutte_auxetic.jl -- the deployability residual (Eq. 2) and the Tutte auxetic
# linear system (Eqs. 3-6) of 2026 Sec. 4.4.
#
# Port of code/src/core/tutte_auxetic.{hpp,cpp}. 1-based indices throughout. The C++
# keeps the system DENSE (Eigen::MatrixXd) and solves it with a dense BDCSVD; the port
# does the same with LinearAlgebra.svd so the rank tolerance semantics carry over.

@enum BoundaryMode None Fixed Periodic

"""Residual of Eq. (2) for one hole preimage: sum of directed hinge edge vectors."""
function hole_residual(c::CutStructure, X::Vector{Vec2}, h::HolePreimage)
    s = Vec2(0, 0)
    for e in h.edges
        c.edge_type[e] != Hinge && continue
        d = c.hinge_dir[e]
        s += X[d.dst] - X[d.src]
    end
    return s
end

mutable struct Residuals
    per_hole::Vector{Vec2}
    max_norm::Float64
    l2_norm::Float64
end
Residuals() = Residuals(Vec2[], 0.0, 0.0)
deployable(r::Residuals, tol::Float64 = 1e-9) = r.max_norm <= tol

function hole_residuals(c::CutStructure, X::Vector{Vec2}, hs::HoleSet)
    r = Residuals()
    for i in hs.interior_indices
        v = hole_residual(c, X, hs.all[i])
        push!(r.per_hole, v)
        r.max_norm = max(r.max_norm, norm(v))
        r.l2_norm += dot(v, v)
    end
    r.l2_norm = sqrt(r.l2_norm)
    return r
end

mutable struct LinearSystem
    A::Matrix{Float64}      # (H + B) x N
    rhs::Matrix{Float64}    # (H + B) x 2
    n_hole_rows::Int
    n_boundary_rows::Int
    N::Int
end
LinearSystem() = LinearSystem(zeros(0, 0), zeros(0, 2), 0, 0, 0)

"""
    assemble_system(c, hs, X_ini, mode) -> LinearSystem

Assembles [L; B] X = [0; T]. Only the all-interior hole preimages contribute
rows of L (holes touching the mesh boundary are unbounded regions, not holes).
"""
function assemble_system(c::CutStructure, hs::HoleSet, X_ini::Vector{Vec2}, mode::BoundaryMode)
    m = c.mesh
    N = n_vertices(m)
    rows = Vector{Float64}[]
    rhs = Vec2[]

    for i in hs.interior_indices
        row = zeros(N)
        for e in hs.all[i].edges
            c.edge_type[e] != Hinge && continue
            row[c.hinge_dir[e].dst] += 1.0
            row[c.hinge_dir[e].src] -= 1.0
        end
        push!(rows, row)
        push!(rhs, Vec2(0, 0))
    end
    n_hole_rows = length(rows)

    if mode == Fixed
        for v in 1:N
            m.vertex_is_boundary[v] || continue
            row = zeros(N)
            row[v] = 1.0
            push!(rows, row)
            push!(rhs, X_ini[v])
        end
    elseif mode == Periodic
        m.periodic.present || error("periodic boundary requested but no spec")
        function add_pairs(prs, t)
            for p in prs
                row = zeros(N)
                row[p[1]] += 1.0
                row[p[2]] -= 1.0
                push!(rows, row)
                push!(rhs, t)
            end
        end
        add_pairs(m.periodic.pairs_h, m.periodic.th)
        add_pairs(m.periodic.pairs_v, m.periodic.tv)
        # Pin one vertex (the first) to kill the global translation.
        row = zeros(N)
        row[1] = 1.0
        push!(rows, row)
        push!(rhs, X_ini[1])
    end

    sys = LinearSystem()
    sys.N = N
    sys.n_hole_rows = n_hole_rows
    sys.n_boundary_rows = length(rows) - n_hole_rows
    sys.A = zeros(length(rows), N)
    sys.rhs = zeros(length(rows), 2)
    for i in eachindex(rows)
        sys.A[i, :] = rows[i]
        sys.rhs[i, 1] = rhs[i][1]
        sys.rhs[i, 2] = rhs[i][2]
    end
    return sys
end

mutable struct SolveReport
    rank_full::Int          # rank of [L; B]
    rank_L::Int             # rank of L alone
    dim_null::Int           # N - rank_full  (number of basis vectors phi_i)
    n_interior::Int
    H::Int                  # number of hole rows
    X0::Matrix{Float64}     # N x 2 particular solution (Eq. 6 projection)
    Phi::Matrix{Float64}    # N x dim_null null-space basis
    sv_tol::Float64
    used_sparse::Bool
    projection_ok::Bool     # A*X0 == rhs to tolerance
    projection_res::Float64
end
SolveReport() = SolveReport(0, 0, 0, 0, 0, zeros(0, 2), zeros(0, 0), 0.0, false, false, 0.0)

# The C++ rank rule: sv > max(rel_tol * smax, 1e-14) * max(1, rows, cols).
function _sv_rank_tol(sv::AbstractVector, rel_tol::Float64, rows::Int, cols::Int)
    smax = isempty(sv) ? 0.0 : sv[1]
    return max(rel_tol * smax, 1e-14) * max(1, max(rows, cols))
end

"""
    solve_system(sys, X_ini, rel_tol = 1e-10) -> SolveReport

Dense SVD path: rank, null space and the Eq. (6) projection of X_ini.
"""
function solve_system(sys::LinearSystem, X_ini::Vector{Vec2}, rel_tol::Float64 = 1e-10)
    rep = SolveReport()
    rep.H = sys.n_hole_rows
    N = sys.N
    A = sys.A
    nr = size(A, 1)

    # Full V (Eigen: ComputeThinU | ComputeFullV): the trailing N - rank columns span the
    # null space. A 0-row system has no singular values; its null space is all of R^N.
    if nr == 0
        U = zeros(0, 0); sv = Float64[]; V = Matrix{Float64}(I, N, N)
    else
        F = svd(A; full = true)
        U = F.U; sv = F.S; V = F.V
    end
    rep.sv_tol = _sv_rank_tol(sv, rel_tol, nr, N)
    rank = count(s -> s > rep.sv_tol, sv)
    rep.rank_full = rank
    rep.dim_null = N - rank

    # Null-space basis: the trailing columns of V.
    rep.Phi = V[:, N-rep.dim_null+1:N]

    # Eq. (6): X0 = X_ini + A^+ (rhs - A X_ini), the minimum-norm correction.
    Xi = points_to_matrix(X_ini)
    r = sys.rhs - A * Xi
    corr = zeros(N, 2)
    for i in 1:rank
        corr += (V[:, i] / sv[i]) * (U[:, i]' * r)
    end
    rep.X0 = Xi + corr
    rep.projection_res = nr == 0 ? 0.0 : maximum(abs, A * rep.X0 - sys.rhs)
    rep.projection_ok = rep.projection_res < 1e-8 * max(1.0, nr == 0 ? 0.0 : maximum(abs, sys.rhs))

    # Rank of L alone (the hole block).
    if sys.n_hole_rows > 0
        s2 = svdvals(A[1:sys.n_hole_rows, :])
        t2 = _sv_rank_tol(s2, rel_tol, sys.n_hole_rows, N)
        rep.rank_L = count(s -> s > t2, s2)
    end
    return rep
end

"""
    rank_only_sparse(sys) -> SolveReport

Rank only, via sparse QR (for large systems where the dense SVD is impractical).
The C++ uses Eigen::SparseQR with `setPivotThreshold(1e-9)`; SuiteSparse SPQR's `tol`
is the analogous column-norm threshold, so `qr(S; tol = 1e-9)` + `rank` is used here.
"""
function rank_only_sparse(sys::LinearSystem)
    rep = SolveReport()
    rep.used_sparse = true
    rep.H = sys.n_hole_rows
    sprank(M) = size(M, 1) == 0 ? 0 : rank(qr(sparse(M); tol = 1e-9))
    rep.rank_full = sprank(sys.A)
    if sys.n_hole_rows > 0
        rep.rank_L = sprank(sys.A[1:sys.n_hole_rows, :])
    end
    rep.dim_null = sys.N - rep.rank_full
    return rep
end

matrix_to_points(X::AbstractMatrix) = [Vec2(X[i, 1], X[i, 2]) for i in 1:size(X, 1)]

function points_to_matrix(X::Vector{Vec2})
    M = Matrix{Float64}(undef, length(X), 2)
    for i in eachindex(X)
        M[i, 1] = X[i][1]
        M[i, 2] = X[i][2]
    end
    return M
end
