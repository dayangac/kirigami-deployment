# core/orientation.jl -- face-orientation assignment sigma (2026 Sec. 4.2, Eq. 1).
#
# Port of code/src/core/orientation.{hpp,cpp}. `assign_orientation_relaxation` reproduces the
# C++ sigma bit-for-bit for a given `MT19937`: the random angles are `uniform_real(rng, 0, 2pi)`
# in the same order, cos/sin are the system libm (see generators.jl) and the Eigen vector
# arithmetic of the projected gradient descent is written out uncontracted.

mutable struct OrientationReport
    sigma::Vector{Int}
    n_split::Int
    n_hinge::Int
    components::Int   # connected components of M'
    n_holes::Int
    objective::Float64  # sum of x_i . x_j over dual edges at the relaxed solution
end
OrientationReport() = OrientationReport(Int[], 0, 0, 0, 0, 0.0)

"""Dual graph of M: one node per face, one edge per adjacent face pair (adjacency lists)."""
function dual_graph(m::Mesh)
    adj = [Int[] for _ in 1:n_faces(m)]
    for e in m.edges
        e.n_faces != 2 && continue
        a = m.half_edges[e.he[1]].face
        b = m.half_edges[e.he[2]].face
        push!(adj[a], b)
        push!(adj[b], a)
    end
    return adj
end

# shallow copy of the mesh with a different sigma (the C++ `Mesh copy = m; copy.sigma = ...`)
function _with_sigma(m::Mesh, sigma::Vector{Int})
    return Mesh(m.X, m.faces, sigma, m.periodic, m.half_edges, m.edges, m.edge_index,
                m.vertex_half_edges, m.vertex_edges, m.vertex_faces, m.vertex_is_boundary)
end

"""Fills in the reported statistics for a given sigma."""
function describe_orientation(m::Mesh, sigma::Vector{Int})
    r = OrientationReport()
    r.sigma = sigma
    c = make_cut(_with_sigma(m, sigma))
    r.n_split = n_split(c)
    r.n_hinge = n_hinge(c)
    r.components = count_components(c)
    r.n_holes = n_interior_holes(holes_partition(c))
    return r
end

"""
    assign_orientation_relaxation(m, rng, restarts=12, iters=800, n_diameters=180)

Eq. (1): minimize sum_{(i,j) in E_d} x_i . x_j with ||x_i|| = 1 in R^2, by projected
gradient descent with several random restarts; then round by the best of a set of
candidate diameters (maximizing the cut), with the paper's principle (1) flip pass and a
connectivity repair.
"""
function assign_orientation_relaxation(m::Mesh, rng::MT19937, restarts::Int = 12,
                                       iters::Int = 800, n_diameters::Int = 180)
    adj = dual_graph(m)
    F = n_faces(m)

    best_x = Vec2[]
    best_obj = Inf
    for _ in 1:restarts
        x = Vector{Vec2}(undef, F)
        for i in 1:F
            a = uniform_real(rng, 0.0, 2 * pi)
            x[i] = Vec2(libm_cos(a), libm_sin(a))
        end
        step = 0.5
        for _ in 1:iters
            grad = fill(Vec2(0.0, 0.0), F)
            for i in 1:F, j in adj[i]
                grad[i] += x[j]  # d/dx_i sum_{(i,j)} x_i . x_j
            end
            for i in 1:F
                g = grad[i]
                v = Vec2(x[i][1] - step * g[1], x[i][2] - step * g[2])
                n = _norm2(v)
                x[i] = n > 1e-12 ? Vec2(v[1] / n, v[2] / n) : x[i]
            end
            step *= 0.999
        end
        obj = 0.0
        for i in 1:F, j in adj[i]
            j > i && (obj += _dot2(x[i], x[j]))
        end
        if obj < best_obj
            best_obj = obj
            best_x = x
        end
    end

    # Round: try many diameters, keep the assignment with the largest cut
    # (equivalently the fewest split edges), breaking ties by keeping M' connected.
    best_sigma = Int[]
    best_split = typemax(Int)
    best_comp = typemax(Int)
    for d in 0:(n_diameters - 1)
        a = pi * d / n_diameters
        nvec = Vec2(libm_cos(a), libm_sin(a))
        sig = [_dot2(best_x[i], nvec) >= 0 ? 1 : -1 for i in 1:F]
        # Paper's principle (1) (Sec. 4.2): a face whose neighbours all carry the same
        # orientation as itself is fully detached after cutting -- flip it.
        for _ in 1:8
            changed = false
            for i in 1:F
                isempty(adj[i]) && continue
                all_same = true
                for j in adj[i]
                    sig[j] != sig[i] && (all_same = false)
                end
                if all_same
                    sig[i] = -sig[i]
                    changed = true
                end
            end
            changed || break
        end
        # Connectivity repair: while M' falls apart, flip a face of the smallest
        # component that borders another component, turning its split cuts into hinges.
        comp = 0
        for _ in 1:30
            cc = try
                make_cut(_with_sigma(m, sig))
            catch
                comp = -1
                break
            end
            comp = count_components(cc)
            comp <= 1 && break
            lab = face_components(cc)
            cnt = zeros(Int, comp)
            for f in 1:F
                cnt[lab[f]] += 1
            end
            small = argmin(cnt)   # first smallest, like std::min_element
            pick = 0
            for f in 1:F
                pick != 0 && break
                lab[f] == small || continue
                for g in adj[f]
                    if lab[g] != small
                        pick = f
                        break
                    end
                end
            end
            pick == 0 && break
            sig[pick] = -sig[pick]
        end
        comp < 0 && continue
        split = 0
        for i in 1:F, j in adj[i]
            (j > i && sig[i] == sig[j]) && (split += 1)
        end
        # connectivity first (principle 1), then split count (principle 2)
        if comp < best_comp || (comp == best_comp && split < best_split)
            best_comp = comp
            best_split = split
            best_sigma = copy(sig)
        end
    end
    isempty(best_sigma) && error("orientation rounding produced no candidate")
    r = describe_orientation(m, best_sigma)
    r.objective = best_obj
    return r
end

"""Exhaustive search over all 2^F assignments (F <= 20): among those keeping M'
connected, take the one with the fewest split edges."""
function brute_force_orientation(m::Mesh)
    F = n_faces(m)
    F > 20 && error("brute_force_orientation: F > 20")
    adj = dual_graph(m)
    best = Int[]
    best_split = typemax(Int)
    for mask in UInt64(0):((UInt64(1) << F) - 1)
        (mask & 1) != 0 && continue  # fix face 1 to sigma = -1 (global flip symmetry)
        sig = [((mask >> (i - 1)) & 1) != 0 ? 1 : -1 for i in 1:F]
        split = 0
        for i in 1:F, j in adj[i]
            (j > i && sig[i] == sig[j]) && (split += 1)
        end
        split >= best_split && continue
        count_components(make_cut(_with_sigma(m, sig))) != 1 && continue  # must stay connected
        best_split = split
        best = sig
    end
    isempty(best) && error("brute_force_orientation: no connected assignment")
    return describe_orientation(m, best)
end
