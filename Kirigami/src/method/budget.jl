# method/budget.jl -- the expansion-budget identity of ideas/round2_theorist_b.md (B1/B3):
# the first-order hole-opening rate of a kirigami design is a BORDER functional on a
# patch and  1/2 det(P_0) tr K  per cell on a torus, and it splits edge by edge into a
# HINGE part (fixed by the flat combinatorics) and a SPLIT part (the design freedom,
# and the quantity whose sign K5/K6/F30 measure).
#
# Port of code/src/method/budget.{hpp,cpp}. 1-based indices throughout. C++ out-pointers
# become extra return values: `face_potential` returns `(u, worst_closure)` and
# `periodic_cell_edges` returns `(keep, n_preimage)`; the C++ `std::vector<char>` masks
# are `Vector{Bool}` (an empty mask selects every edge).
#
# THE IDENTITY.  With u the face potential of derivations/core.md T1
# (u_g - u_f = sigma_g x_src(e) across a hinge edge) the deployed area of the preimage
# walk bounding a hole/notch C is the first harmonic
#
#      A_C(theta) = a_C sin(theta) - b_C (1 - cos(theta)) ,   A_C(0) = 0 ,
#      a_C = 1/2 sum_{(a->b) in dC, face f} <x_a - x_b, u_f> ,
#      b_C = 1/2 sum_{(a->b) in dC, face f} sigma_f det(x_a - x_b, u_f) .
#
# Both copies of an interior edge lie on one preimage walk (F11) and are traversed with
# the hole on the left, i.e. REVERSED with respect to the stored face order.  Writing
# f0 = face(he[0]), f1 = face(he[1]) and (a -> b) = (he[0].from -> he[0].to) for the
# stored direction inside f0, the two copies together contribute
#
#      1/2 <x_a - x_b, u_{f1} - u_{f0}>                                        (*)
#
# to a_C, which is invariant under swapping f0 and f1 (both factors change sign).  For a
# HINGE edge u_{f1} - u_{f0} = sigma_{f1} x_src(e), so (*) needs no potential at all:
#
#      w_e := sigma_{f1} <x_a - x_b, x_src(e)>            (hinge term, u-free)
#      r_e := <x_a - x_b, u_{f1} - u_{f0}>                (split term)
#
# and  2 sum_C a_C = W + R  with  W = sum_hinge w_e,  R = sum_split r_e.
#
# Summing (*) over every preimage is a sum over every interior half-edge taken in the
# REVERSED direction; adding the (vanishing) forward sum over each closed face polygon
# leaves only the border half-edges, so on a disk-topology patch
#
#      W + R = 2 B(X) ,   B(X) := 1/2 sum_{border (a->b) in face f} <x_a - x_b, u_f> ,
#
# and on a torus, where the same telescoping leaves the period monodromy instead of a
# border (see method/periodic_jacobian.jl and `periodic_cell_edges` below, which says
# which LIFT of each edge to sum),
#
#      W + R = det(P_0) tr K .                                             (B.5)
#
# THE SPLIT TERM IS THE 0+ MARGIN.  method/zero_plus.jl measures
# q_e = det(dS_e, d_e) with d_e = x_b - x_a and dS_e = 2 Jrot (u_{f1} - u_{f0}); since
# det(Jrot p, v) = -<p, v>, this is  q_e = 2 r_e  EXACTLY (asserted by a doctest).  So
# "every split cut opens at 0+" is  min_e r_e > 0, and averaging over the n_split cuts,
#
#      min_e r_e <= R / n_split = (det(P_0) tr K - W) / n_split ,
#
# i.e. tr K < tau* := W / det(P_0) forces some cut to fold inward and Theta_max = 0.
# That threshold is B3.  NOTE what (B.5) does to it: tr K - tau* = R / det(P_0) exactly,
# so the "budget gap" is the total split budget rescaled and is a QUADRATIC function of
# the design, not an affine one -- tr K >= tau* is not a linear constraint on the
# achievable set.  kill_b3 measures both halves of that.

# _det2 lives in core/mesh.jl (shared).

# The two faces of an interior edge, in the canonical order (f0 supplies the stored
# half-edge direction a -> b, exactly as method/zero_plus.jl's q_e convention does).
# Returns (ok, f0, f1, a, b).
function _sides(m::Mesh, e::Int)
    ed = m.edges[e]
    ed.n_faces != 2 && return (false, 0, 0, 0, 0)
    h0 = m.half_edges[ed.he[1]]
    return (true, h0.face, m.half_edges[ed.he[2]].face, h0.from, h0.to)
end

"""
    face_potential(c, X) -> (u, worst_closure)

The face potential u of T1, by BFS over the hinge-connected dual graph (one root per
component, u = 0 there).  `worst_closure` is the largest violation of
u_g - u_f = sigma_g x_src over ALL hinge edges -- the witness that u exists.
"""
function face_potential(c::CutStructure, X::Vector{Vec2})
    m = c.mesh
    F = n_faces(m)
    adj = [Tuple{Int,Int}[] for _ in 1:F]
    for e in c.hinge_edges
        ed = m.edges[e]
        ed.n_faces != 2 && continue
        f1 = m.half_edges[ed.he[1]].face
        f2 = m.half_edges[ed.he[2]].face
        push!(adj[f1], (e, f2))
        push!(adj[f2], (e, f1))
    end
    u = fill(Vec2(0, 0), F)
    seen = falses(F)
    for s in 1:F
        seen[s] && continue
        seen[s] = true
        u[s] = Vec2(0, 0)
        q = Int[s]
        qi = 1
        while qi <= length(q)
            f = q[qi]
            qi += 1
            for (e, g) in adj[f]
                seen[g] && continue
                seen[g] = true
                u[g] = u[f] + Float64(m.sigma[g]) * X[c.hinge_dir[e].src]
                push!(q, g)
            end
        end
    end
    w = 0.0
    for e in c.hinge_edges
        ed = m.edges[e]
        ed.n_faces != 2 && continue
        f1 = m.half_edges[ed.he[1]].face
        f2 = m.half_edges[ed.he[2]].face
        w = max(w, norm(u[f2] - u[f1] - Float64(m.sigma[f2]) * X[c.hinge_dir[e].src]))
    end
    return u, w
end

# Per-edge budget terms.  `hinge` is indexed like c.hinge_edges, `split` like
# c.split_edges (NOT like split_copies, which drops non-manifold split edges);
# `b_hinge` / `b_split` are the matching (1 - cos) coefficients
#   1/2 [ sigma_{f0} det(x_a - x_b, u_{f0}) - sigma_{f1} det(x_a - x_b, u_{f1}) ] * 2 .
mutable struct BudgetTerms
    hinge::Vector{Float64}      # w_e
    split::Vector{Float64}      # r_e
    b_hinge::Vector{Float64}    # the (1 - cos theta) counterparts
    b_split::Vector{Float64}
    W::Float64                  # sums
    R::Float64
    bW::Float64
    bR::Float64
    u_closure::Float64
end
BudgetTerms() = BudgetTerms(Float64[], Float64[], Float64[], Float64[], 0.0, 0.0, 0.0, 0.0, 0.0)

"""
    budget_terms(c, X, u, keep = Bool[]) -> BudgetTerms

Per-edge budget terms w_e (hinge) and r_e (split) with their (1 - cos) partners.
`keep[e]` restricts the sums to a caller-supplied edge set (used on the periodic super
patch, where exactly one representative of each quotient edge is summed); an empty
`keep` selects every edge.
"""
function budget_terms(c::CutStructure, X::Vector{Vec2}, u::Vector{Vec2},
                      keep::Vector{Bool} = Bool[])
    m = c.mesh
    t = BudgetTerms()
    function term(e::Int)
        ok, f0, f1, a, b = _sides(m, e)
        ok || return (0.0, 0.0)
        d = X[a] - X[b]
        return (dot(d, u[f1] - u[f0]),
                _det2(d, Float64(m.sigma[f1]) * u[f1] - Float64(m.sigma[f0]) * u[f0]))
    end
    for e in c.hinge_edges
        a, b = (isempty(keep) || keep[e]) ? term(e) : (0.0, 0.0)
        push!(t.hinge, a)
        push!(t.b_hinge, b)
        t.W += a
        t.bW += b
    end
    for e in c.split_edges
        a, b = (isempty(keep) || keep[e]) ? term(e) : (0.0, 0.0)
        push!(t.split, a)
        push!(t.b_split, b)
        t.R += a
        t.bR += b
    end
    return t
end

# ---------------------------------------------------------------- periodic
#
# THE LIFT MATTERS.  Both copies of an interior edge contribute
# <x_a - x_b, u_{f1} - u_{f0}> to the a-coefficient of the preimage carrying them.  For a
# SPLIT edge sigma_{f0} = sigma_{f1}, so u_{f1} - u_{f0} is unchanged by a lattice
# translation and the term is a function of the quotient edge alone.  For a HINGE edge
# sigma_{f0} = -sigma_{f1} and u_{f1+t} - u_{f0+t} = u_{f1} - u_{f0} + sigma_{f1} t, so the
# term depends on WHICH LIFT the preimage walk uses.  Summing an arbitrary representative
# per quotient edge therefore misses a defect that is measured, not small (up to
# 1.9 det P_0 on the K7 population).  The budget of one cell is the sum over the H
# PREIMAGES anchored in a reference cell, each at the lift its own walk uses.

"""
    periodic_cell_edges(c, sp) -> (keep, n_preimage)

A per-edge mask of the super patch selecting exactly the edges of one representative
preimage per translation class: an all-interior preimage is accepted iff none of its
edges' quotient classes is already covered.  Since the preimages partition
E_hinge u E_split and an all-interior preimage of the 3x3 patch carries its whole class,
the result does not depend on the order of acceptance.  `n_preimage` is the count; the
caller should check that the mask covers q.n_hinge + q.n_split edges (it does not on a
super patch too small for the pattern, e.g. squares_3x2 / squares_3x3).
"""
function periodic_cell_edges(c::CutStructure, sp::SuperPatch)
    m = c.mesh
    function ekey(e::Int)
        k = m.edges[e].key
        a = sp.vert_class[k.a]
        b = sp.vert_class[k.b]
        d = sp.vert_off[k.b] - sp.vert_off[k.a]
        if b < a || (a == b && (d[1] < 0 || (d[1] == 0 && d[2] < 0)))
            a, b = b, a
            d = -d
        end
        return (a, b, d[1], d[2])
    end
    hs = holes_partition(c)
    covered = Set{NTuple{4,Int}}()
    keep = falses(n_edges(m))
    npre = 0
    for pre in hs.all
        (!pre.all_interior || isempty(pre.edges)) && continue
        any(e -> ekey(e) in covered, pre.edges) && continue
        for e in pre.edges
            push!(covered, ekey(e))
            keep[e] = true
        end
        npre += 1
    end
    return Vector{Bool}(keep), npre
end

# B(X) and its (1 - cos) partner, summed over the border half-edges in stored face order.
# {a_border, b_border} with  sum_C a_C = a_border  and  sum_C b_C = b_border.
mutable struct BorderFunctional
    a::Float64
    b::Float64
end
BorderFunctional() = BorderFunctional(0.0, 0.0)

"""B(X) over the border half-edges in stored face order, with its (1 - cos) partner."""
function border_functional(c::CutStructure, X::Vector{Vec2}, u::Vector{Vec2})
    m = c.mesh
    B = BorderFunctional()
    for e in c.border_edges
        ed = m.edges[e]
        ed.n_faces != 1 && continue
        h = m.half_edges[ed.he[1]]
        d = X[h.from] - X[h.to]
        B.a += 0.5 * dot(d, u[h.face])
        B.b += 0.5 * Float64(m.sigma[h.face]) * _det2(d, u[h.face])
    end
    return B
end

"""
    measured_void_area(c, Y)

The MEASURED total void area at opening angle theta: the shoelace of the deployed
border half-edges minus the (theta-independent) sum of the deployed face areas.  Equal
to the sum of the areas of every hole AND notch, with no walk tracing at all.
"""
function measured_void_area(c::CutStructure, Y::Vector{Vec2})
    m = c.mesh
    border = 0.0
    for e in c.border_edges
        ed = m.edges[e]
        ed.n_faces != 1 && continue
        h = m.half_edges[ed.he[1]]
        p = Y[prime_vertex(c, h.face, h.from)]
        q = Y[prime_vertex(c, h.face, h.to)]
        border += 0.5 * _det2(p, q)
    end
    faces = 0.0
    for f in 1:n_faces(m)
        pf = c.prime_faces[f]
        s = 0.0
        for i in eachindex(pf)
            s += _det2(Y[pf[i]], Y[pf[mod1(i + 1, length(pf))]])
        end
        faces += 0.5 * s
    end
    return border - faces
end
