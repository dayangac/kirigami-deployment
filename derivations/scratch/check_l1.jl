# check_l1.jl -- Deriver-L, mission 2 / WP1.  Tests Lemmas L1.1, L1.2, L1.3 of
# derivations/lemmas.md: the COMBINATORIAL classification of contact harmonics.
# Port of check_l1.cpp (same corpus -- frozen in data/corpus/ -- same seed, same printed lines).
#
# Run (from repo root):
#   julia --project=Kirigami derivations/scratch/check_l1.jl [NRAND=52] [NSAMP=20]
#
# What is predicted (lemmas.md L1.0-L1.3).  For a candidate (w, (a,b)) with (a,b) an
# edge of face f and w a corner of face g != f, write va, vb, vw for the M-vertices the
# three M'-vertices are copies of, d := x_vb - x_va, and -- when vw == va or vw == vb --
# dS := S_w - S_(the copy of vw inside f).  Then
#
#     C = p + q = det(x_vb - x_va, x_vw - x_va)          (the FLAT determinant)
#
# and, in the coincident-copy case vw in {va, vb},
#
#     C = 0 identically,   p = (sigma_f/2) <d, dS>,   q = -p,   r = (1/2) det(d, dS).
#
# Hence the whole class of such a pair is decided by dS alone:
#     dS = 0                          -> h identically zero      (L1.3)
#     dS != 0, det(d,dS) = 0          -> class 3, never a contact (L1.2)
#     det(d,dS) != 0                  -> class 2                  (L1.1)
# and a pair with vw not in {va, vb} is class 1 unless the three FLAT source vertices
# happen to be collinear -- structurally (frozen triple) or accidentally (isolated t).
#
# The program measures: (i) the combinatorial C=0 predicate against |C| <= tol;
# (ii) the three closed forms above; (iii) the Zero/3/2 split against the numeric class;
# (iv) for every disagreement, whether it persists over all shape-space samples of that
# graph (structural) or occurs at isolated t (accidental).

using Kirigami
using LinearAlgebra
using Printf
using StaticArrays
const K = Kirigami
const Vec2 = K.Vec2
include(joinpath(@__DIR__, "corpus_common.jl"))

cross2(a::Vec2, b::Vec2) = a[1] * b[2] - a[2] * b[1]
dot2(a::Vec2, b::Vec2) = a[1] * b[1] + a[2] * b[2]

struct Harm
    p::Float64
    q::Float64
    r::Float64
end
hscale(h::Harm) = abs(h.p) + hypot(h.q, h.r)

function orient_h(Cab::Vec2, Sab::Vec2, Caw::Vec2, Saw::Vec2)
    dcc = cross2(Cab, Caw); dss = cross2(Sab, Saw)
    Harm(0.5 * (dcc + dss), 0.5 * (dcc - dss), 0.5 * (cross2(Cab, Saw) + cross2(Sab, Caw)))
end

const K_ZERO = 0; const K_C1 = 1; const K_C2 = 2; const K_C3 = 3
kname(k::Int) = k == K_ZERO ? "zero" : k == K_C1 ? "class1" : k == K_C2 ? "class2" : "class3"

# numeric class from (p,q,r) -- the rule the code uses (deploy_basis.jl)
function numeric_class(h::Harm, tol::Float64)
    C = h.p + h.q; B = 2 * h.r; A = h.p - h.q
    c0 = abs(C) <= tol; b0 = abs(B) <= tol; a0 = abs(A) <= tol
    (c0 && b0 && a0) && return K_ZERO
    c0 || return K_C1
    return b0 ? K_C3 : K_C2
end

Base.@kwdef mutable struct Tally
    pairs::Int = 0
    # L1.1 : combinatorial C = 0 predicate
    shared::Int = 0
    czero_numeric::Int = 0
    agree_c::Int = 0
    shared_but_C_nonzero::Int = 0      # must be 0
    cz_not_shared::Int = 0             # flat-collinear triples (structural or accidental)
    cz_not_shared_structural::Int = 0  # persists over every sample of its graph
    cz_not_shared_frozen::Int = 0      # ... and all three source vertices are frozen
    # L1.2/L1.3 : class prediction on the coincident-copy pairs
    agree_class::Int = 0
    disagree_class::Int = 0
    n_zero::Int = 0
    n_c1::Int = 0
    n_c2::Int = 0
    n_c3::Int = 0
    c3_persist::Int = 0
    c3_transient::Int = 0
    c3_single::Int = 0
    cz_ns_single::Int = 0
    # closed forms
    max_err_p::Float64 = 0.0
    max_err_q::Float64 = 0.0
    max_err_r::Float64 = 0.0
    max_err_hinge_dS::Float64 = 0.0
    max_err_split_dS::Float64 = 0.0
    n_hinge_dS::Int = 0
    n_split_dS::Int = 0
    # L1.3 structure of the identically-zero pairs
    zero_same_prime::Int = 0
    zero_shared::Int = 0
    zero_hinge_src::Int = 0
    zero_split_zero_du::Int = 0
    zero_other::Int = 0
    zero_lambda_interior::Int = 0      # h == 0 with vw not in {va, vb}
    # shared-edge type of a coincident-copy pair
    sh_hinge::Int = 0
    sh_split::Int = 0
    sh_none::Int = 0
    # class 3 via the hinge formula:  class3 <=> <d, x_src - x_v> = 0  (d perp the hinge edge)
    c3_hinge::Int = 0
    c3_hinge_perp_ok::Int = 0
    c3_by_type::Vector{Int} = zeros(Int, 3)
    max_err_split_dS_const::Float64 = 0.0
    n_split_dS_const::Int = 0
    # persistent non-shared C == 0, by how many of the three source vertices are frozen
    cz_ns_frozen_cnt::Vector{Int} = zeros(Int, 4)
    ex_printed::Int = 0
    disagreements::Dict{String,Int} = Dict{String,Int}()
end

function positively_oriented(m::K.Mesh, X::Vector{Vec2})
    for f in 1:K.n_faces(m)
        A2 = 0.0
        Fv = m.faces[f]
        for i in eachindex(Fv)
            a = X[Fv[i]]; b = X[Fv[i % length(Fv) + 1]]
            A2 += a[1] * b[2] - a[2] * b[1]
        end
        A2 <= 0 && return false
    end
    return true
end

# the interior edge shared by faces f and g at vertex vw, or 0
function shared_edge(m::K.Mesh, vw::Int, f::Int, g::Int)
    e_sh = 0
    for e in m.vertex_edges[vw]
        E = m.edges[e]
        E.n_faces != 2 && continue
        fa = m.half_edges[E.he[1]].face; fb = m.half_edges[E.he[2]].face
        ((fa == f && fb == g) || (fa == g && fb == f)) && (e_sh = e)
    end
    return e_sh
end

# ---------------------------------------------------------------------------
# One graph, many shape-space samples.  Candidates are enumerated in a FIXED order so
# the per-candidate persistence counters line up across samples.
function run_graph(name::String, m::K.Mesh, n_samp::Int, rng::K.MT19937, T::Tally, used::Ref{Int})
    K.build_topology!(m)
    isempty(m.sigma) && return
    c = K.make_cut(m)
    hs = K.holes_partition(c)
    sys = K.assemble_system(c, hs, m.X, K.Fixed)
    sr = K.solve_system(sys, m.X)
    sr.projection_ok || return
    X0 = K.matrix_to_points(sr.X0)
    K.deployable(K.hole_residuals(c, X0, hs), 1e-7) || return
    positively_oriented(m, X0) || return
    used[] += 1

    F = K.n_faces(m)
    scale = 0.0
    for p in X0
        scale = max(scale, norm(p))
    end

    # frozen vertices: zero row of Phi
    frozen = trues(K.n_vertices(m))
    for v in 1:K.n_vertices(m), j in 1:sr.dim_null
        if abs(sr.Phi[v, j]) > 1e-12
            frozen[v] = false
            break
        end
    end

    # pass 0: count candidates
    ncand = 0
    for f in 1:F, k in 1:length(c.prime_faces[f]), g in 1:F
        g != f && (ncand += length(c.prime_faces[g]))
    end
    cz_count = zeros(Int, ncand); c3_count = zeros(Int, ncand); samp_seen = zeros(Int, ncand)

    gauss = K.NormalDist(0.0, 1.0)
    done = 0
    for s in 0:n_samp-1
        X = X0
        if s > 0 && sr.dim_null > 0
            t = Matrix{Float64}(undef, sr.dim_null, 2)
            for i in 1:sr.dim_null, j in 1:2
                t[i, j] = 0.03 * scale * K.normal(gauss, rng)
            end
            X = K.matrix_to_points(sr.X0 + sr.Phi * t)
        elseif s > 0
            break   # rigid shape space: one sample only
        end
        positively_oriented(m, X) || continue
        done += 1
        dp = K.deploy(c, X, 0.0)
        Cv = dp.Y
        Sv = [2.0 * dp.dY_dtheta[i] for i in 1:c.n_prime_vertices]
        sc = 0.0
        for p in X
            sc = max(sc, norm(p))
        end
        tol2 = 1e-11 * sc * sc   # determinant scale
        tol1 = 1e-11 * sc        # vector scale

        idx = 0
        for f in 1:F
            pf = c.prime_faces[f]
            nf = length(pf)
            for kk in 1:nf
                a = pf[kk]; b = pf[kk % nf + 1]
                va = c.prime_to_original[a]; vb = c.prime_to_original[b]
                Cab = Cv[b] - Cv[a]; Sab = Sv[b] - Sv[a]
                d = X[vb] - X[va]
                for g in 1:F
                    g == f && continue
                    for w in c.prime_faces[g]
                        idx += 1
                        vw = c.prime_to_original[w]
                        Caw = Cv[w] - Cv[a]; Saw = Sv[w] - Sv[a]
                        h = orient_h(Cab, Sab, Caw, Saw)
                        T.pairs += 1
                        samp_seen[idx] += 1
                        C = h.p + h.q
                        cz = abs(C) <= tol2
                        sh = (vw == va || vw == vb)
                        cz && (cz_count[idx] += 1)
                        sh && (T.shared += 1)
                        cz && (T.czero_numeric += 1)
                        sh == cz && (T.agree_c += 1)
                        (sh && !cz) && (T.shared_but_C_nonzero += 1)
                        (!sh && cz) && (T.cz_not_shared += 1)

                        kn = numeric_class(h, tol2)
                        if kn == K_ZERO
                            T.n_zero += 1
                        elseif kn == K_C1
                            T.n_c1 += 1
                        elseif kn == K_C2
                            T.n_c2 += 1
                        else
                            T.n_c3 += 1
                        end
                        kn == K_C3 && (c3_count[idx] += 1)

                        if !sh
                            kn == K_ZERO && (T.zero_lambda_interior += 1)
                            continue
                        end
                        # ---- coincident-copy pair: the closed forms and the predicted class
                        acopy = (vw == va) ? a : b
                        dS = Sv[w] - Sv[acopy]
                        sf = Float64(m.sigma[f])
                        pp = 0.5 * sf * dot2(d, dS)
                        rr = 0.5 * cross2(d, dS)
                        T.max_err_p = max(T.max_err_p, abs(h.p - pp))
                        T.max_err_q = max(T.max_err_q, abs(h.q + pp))
                        T.max_err_r = max(T.max_err_r, abs(h.r - rr))

                        kp = norm(dS) <= tol1 ? K_ZERO :
                             abs(cross2(d, dS)) <= tol2 ? K_C3 : K_C2
                        if kp == kn
                            T.agree_class += 1
                        else
                            T.disagree_class += 1
                            key = @sprintf("%s: pred %s, num %s (shared=%d)", name, kname(kp), kname(kn), Int(sh))
                            T.disagreements[key] = get(T.disagreements, key, 0) + 1
                        end
                        let et = -1
                            e = shared_edge(m, vw, f, g)
                            if e > 0
                                et = c.edge_type[e] == K.Hinge ? 0 : c.edge_type[e] == K.Split ? 1 : 2
                            end
                            if et == 0
                                T.sh_hinge += 1
                            elseif et == 1
                                T.sh_split += 1
                            else
                                T.sh_none += 1
                            end
                            kn == K_C3 && (T.c3_by_type[(et < 0 ? 2 : et) + 1] += 1)
                        end
                        if kn == K_ZERO
                            T.zero_shared += 1
                            (w == a || w == b) && (T.zero_same_prime += 1)
                            # which structural reason?
                            e_shared = shared_edge(m, vw, f, g)
                            if e_shared > 0 && c.edge_type[e_shared] == K.Hinge &&
                               c.hinge_dir[e_shared].src == vw
                                T.zero_hinge_src += 1
                            elseif e_shared > 0 && c.edge_type[e_shared] == K.Split
                                T.zero_split_zero_du += 1
                            else
                                T.zero_other += 1
                            end
                        end
                        # ---- dS closed forms across a shared edge
                        e_sh = shared_edge(m, vw, f, g)
                        if e_sh > 0 && c.edge_type[e_sh] == K.Split
                            # T1.B: dS_e is the SAME vector for both endpoints of the split edge
                            E = m.edges[e_sh]
                            vo = (E.key.a == vw) ? E.key.b : E.key.a
                            wo = K.prime_vertex(c, g, vo); ao = K.prime_vertex(c, f, vo)
                            if wo >= 1 && ao >= 1
                                dS2 = Sv[wo] - Sv[ao]
                                T.max_err_split_dS_const = max(T.max_err_split_dS_const, norm(dS - dS2))
                                T.n_split_dS_const += 1
                            end
                        end
                        if e_sh > 0 && c.edge_type[e_sh] == K.Hinge
                            src = c.hinge_dir[e_sh].src
                            sg = Float64(m.sigma[g])
                            dd = X[src] - X[vw]
                            pred = Vec2(-dd[2], dd[1]) * (2.0 * sg)
                            T.max_err_hinge_dS = max(T.max_err_hinge_dS, norm(dS - pred))
                            T.n_hinge_dS += 1
                            if kn == K_C3
                                T.c3_hinge += 1
                                abs(dot2(d, dd)) <= tol2 && (T.c3_hinge_perp_ok += 1)
                            end
                        end
                    end
                end
            end
        end
    end
    # recount: which of the "C=0 but not shared" are structural, and frozen;
    # class-3 persistence for shared pairs too
    idx = 0
    for f in 1:F
        pf = c.prime_faces[f]
        nf = length(pf)
        for kk in 1:nf
            a = pf[kk]; b = pf[kk % nf + 1]
            va = c.prime_to_original[a]; vb = c.prime_to_original[b]
            for g in 1:F
                g == f && continue
                for w in c.prime_faces[g]
                    idx += 1
                    vw = c.prime_to_original[w]
                    sh = (vw == va || vw == vb)
                    samp_seen[idx] == 0 && continue
                    if sh
                        if samp_seen[idx] < 5
                            T.c3_single += c3_count[idx]
                            continue
                        end
                        if c3_count[idx] > 0
                            if c3_count[idx] == samp_seen[idx]
                                T.c3_persist += c3_count[idx]
                            else
                                T.c3_transient += c3_count[idx]
                            end
                        end
                        continue
                    end
                    if samp_seen[idx] < 5                 # rigid graph: one t only, cannot decide
                        T.cz_ns_single += cz_count[idx]
                        T.c3_single += c3_count[idx]
                        continue
                    end
                    if cz_count[idx] == samp_seen[idx]
                        T.cz_not_shared_structural += samp_seen[idx]
                        nf3 = Int(frozen[va]) + Int(frozen[vb]) + Int(frozen[vw])
                        T.cz_ns_frozen_cnt[nf3 + 1] += samp_seen[idx]
                        nf3 == 3 && (T.cz_not_shared_frozen += samp_seen[idx])
                        if nf3 < 3 && T.ex_printed < 8
                            T.ex_printed += 1
                            # indices printed 0-based, as the C++ did
                            @printf("    [example] %s  f=%d g=%d  va=%d(%s) vb=%d(%s) vw=%d(%s)  bnd=%d%d%d\n",
                                    name, f - 1, g - 1, va - 1, frozen[va] ? "frz" : "free",
                                    vb - 1, frozen[vb] ? "frz" : "free", vw - 1, frozen[vw] ? "frz" : "free",
                                    Int(m.vertex_is_boundary[va]), Int(m.vertex_is_boundary[vb]),
                                    Int(m.vertex_is_boundary[vw]))
                        end
                    end
                    if c3_count[idx] > 0
                        if c3_count[idx] == samp_seen[idx]
                            T.c3_persist += c3_count[idx]
                        else
                            T.c3_transient += c3_count[idx]
                        end
                    end
                end
            end
        end
    end
end

function main()
    NRAND = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 52
    NSAMP = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 20
    rng = K.MT19937(20260908)
    T = Tally()
    used = Ref(0)

    for rc in reference_cases()
        isempty(rc.mesh.sigma) && continue
        run_graph("ref:" * rc.name, rc.mesh, NSAMP, rng, T, used)
    end
    for id in 0:NRAND-1
        g = make_graph(id, 18, 46, 220)
        g.ok || continue
        run_graph("rand:" * g.kind * string(id), g.mesh, NSAMP, rng, T, used)
    end

    @printf("\ncheck_l1 -- graphs used: %d  (random requested %d, samples/graph %d)\n",
            used[], NRAND, NSAMP)
    @printf("  candidate harmonics examined            : %d\n", T.pairs)
    @printf("\nL1.1  combinatorial predicate  C == 0  <=>  vw in {va, vb}\n")
    @printf("  coincident-copy pairs (vw in {va,vb})   : %d\n", T.shared)
    @printf("  numerically |C| <= tol                  : %d\n", T.czero_numeric)
    @printf("  AGREEMENTS                              : %d / %d\n", T.agree_c, T.pairs)
    @printf("  shared but C != 0  (must be 0)          : %d\n", T.shared_but_C_nonzero)
    @printf("  C == 0 but NOT shared (flat-collinear)  : %d\n", T.cz_not_shared)
    @printf("     ... on rigid graphs (dim_null = 0, one sample, undecidable) : %d\n", T.cz_ns_single)
    @printf("     ... of these, persistent over every sample of the graph : %d\n", T.cz_not_shared_structural)
    @printf("     ... of these, all three source vertices frozen          : %d\n", T.cz_not_shared_frozen)
    @printf("\nL1.2/L1.3  class of a coincident-copy pair from dS alone\n")
    @printf("  numeric classes over ALL pairs: zero %d  class1 %d  class2 %d  class3 %d\n",
            T.n_zero, T.n_c1, T.n_c2, T.n_c3)
    @printf("  predicted == numeric on coincident-copy pairs : %d agree, %d disagree\n",
            T.agree_class, T.disagree_class)
    @printf("  closed forms  max |p - sf<d,dS>/2| = %.3e   max |q + p| = %.3e   max |r - det(d,dS)/2| = %.3e\n",
            T.max_err_p, T.max_err_q, T.max_err_r)
    @printf("  hinge dS = 2 sigma_g J (x_src - x_v):  %d pairs, max err %.3e\n",
            T.n_hinge_dS, T.max_err_hinge_dS)
    @printf("  class 3 occurrences: persistent %d   transient (isolated t) %d   undecidable (rigid graph, 1 sample) %d\n",
            T.c3_persist, T.c3_transient, T.c3_single)
    @printf("\nL1.3  the identically-zero pairs\n")
    @printf("  zero & coincident-copy                  : %d\n", T.zero_shared)
    @printf("     ... of which w IS the same M'-vertex as a or b : %d\n", T.zero_same_prime)
    @printf("     shared edge is a HINGE with src == v : %d\n", T.zero_hinge_src)
    @printf("     shared edge is a SPLIT (Delta u = 0) : %d\n", T.zero_split_zero_du)
    @printf("     no shared edge between f and g       : %d\n", T.zero_other)
    @printf("  zero with vw NOT in {va,vb} (lambda interior) : %d\n", T.zero_lambda_interior)
    @printf("\n  coincident-copy pairs by shared-edge type: hinge %d  split %d  none %d\n",
            T.sh_hinge, T.sh_split, T.sh_none)
    @printf("  class-3 pairs by shared-edge type: hinge %d  split %d  none %d\n",
            T.c3_by_type[1], T.c3_by_type[2], T.c3_by_type[3])
    @printf("  class-3 pairs across a HINGE: %d, of which <d, x_src - x_v> = 0 : %d\n",
            T.c3_hinge, T.c3_hinge_perp_ok)
    @printf("  split dS constant over the two endpoints (T1.B): %d pairs, max err %.3e\n",
            T.n_split_dS_const, T.max_err_split_dS_const)
    @printf("  persistent non-shared C==0, by #frozen source vertices: 0:%d  1:%d  2:%d  3:%d\n",
            T.cz_ns_frozen_cnt[1], T.cz_ns_frozen_cnt[2], T.cz_ns_frozen_cnt[3], T.cz_ns_frozen_cnt[4])
    if !isempty(T.disagreements)
        @printf("\n  DISAGREEMENT TYPES:\n")
        for k in sort!(collect(keys(T.disagreements)))
            @printf("    %-70s  x %d\n", k, T.disagreements[k])
        end
    end
    @printf("\n")
end

main()
