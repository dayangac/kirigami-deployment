# K3a recheck -- are the 15 identity violations combinatorial, or a rank-threshold
# artifact of the sparse-QR path?
#
# Every violation has |dim ker A - (|F\core2| + dim ker A|core2)| == 1 and F >= 879,
# i.e. above mobility's `dense_limit = 700` where matrix_rank switches from the dense
# pivoted QR to the sparse QR. This driver recomputes both ranks on the violators with the
# DENSE decomposition and with a threshold sweep, and reports the singular gap around the
# cut, so the cause is decided rather than assumed.
#
#   julia --project=Kirigami Kirigami/apps/kill_k3a_recheck.jl [--out DIR] [--maxf 1500]
#         [--limit K] [--regenerate]
#
# Outputs go to results/kill/k3a/ by default.
include(joinpath(@__DIR__, "kill_common.jl"))
using SparseArrays

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# Rank by dense column-pivoted QR at a given relative threshold
# (|R_ii| > rel * max_i |R_ii|).
function rank_dense(A::SparseMatrixCSC{Float64,Int}, rel::Float64)
    D = Matrix(A)
    (size(D, 1) == 0 || size(D, 2) == 0) && return 0
    R = qr(D, ColumnNorm()).R
    d = [abs(R[i, i]) for i in 1:min(size(R)...)]
    return count(>(rel * maximum(d)), d)
end

# Singular gap at the cut, to say whether the rank is even well defined.
function sv_gap(M::SparseMatrixCSC{Float64,Int}, r::Int)
    (r <= 0 || r >= min(size(M)...)) && return NaN
    sv = svdvals(Matrix(M))
    r > length(sv) && return NaN
    return sv[r] / max(1e-300, sv[1])
end

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    ids = [40, 43, 88, 247, 257, 266, 299, 320, 343, 425, 445, 455, 467, 476, 494]
    outdir = joinpath(REPO, "results", "kill", "k3a")
    maxf = 1500
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--maxf" && i < length(args); maxf = arg_i(args[i+1]); i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    csv = open(joinpath(outdir, "k3a_recheck.csv"), "w")
    print(csv, "id,kind,F,cols,n_dangling,sparse_kerA,sparse_kerA_core,sparse_defect,",
          "dense_kerA,dense_kerA_core,dense_defect,dense_ok,",
          "sv_gap_full,sv_gap_core,secs\n")
    sm = open(joinpath(outdir, "recheck_summary.txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("K3a recheck: the 15 identity violations, dense vs SparseQR rank\n")
    both("identity: dim ker A == |F \\ core2| + dim ker A|core2\n\n")
    n_tried = 0; n_dense_ok = 0; n_still_bad = 0; n_skipped = 0

    # exactly K3a's population call: make_graph(id, 100, 5000, 1 << 30) = frozen k3a_500
    pop = regenerate ? nothing : load_population("k3a_500")
    for id in ids[1:min(limit, length(ids))]
        g = if pop === nothing
            K.make_graph(id, 100, 5000, 1 << 30)
        else
            r = population_row(pop, id)
            K.Graph(r.mesh, r.kind, r.id, r.ok)
        end
        if !g.ok
            both("  id $id: regeneration failed\n")
            continue
        end
        m = g.mesh
        K.build_topology!(m)
        if K.n_faces(m) > maxf
            n_skipped += 1
            both("  id $id ($(g.kind), F=$(K.n_faces(m))): skipped, dense QR of $(K.n_faces(m)) columns too large\n")
            continue
        end
        t = Timer()
        c = K.make_cut(m)
        hg = K.build_hinge_graph(c)
        # Reproduce K3a's configuration exactly: the Eq. (6) projection where the dense SVD
        # is affordable (N <= 1400), else X_ini. rank(A) depends on the pin positions, so
        # rechecking a different X would not be rechecking the same matrix.
        X = copy(m.X)
        src = "X_ini"
        if K.n_vertices(m) <= 1400
            hs = K.holes_partition(c)
            sr = K.solve_system(K.assemble_system(c, hs, m.X, K.Fixed), m.X)
            if sr.projection_ok
                X = K.matrix_to_points(sr.X0)
                src = "X0"
            end
        end
        pins = K.pins_flat(c, hg, X)
        A = K.build_A(hg, pins)
        tc = K.two_core(hg)
        cg, sub_pins = K.subgraph(hg, tc.in_core, tc.edge_in_core, pins)
        Ac = K.build_A(cg, sub_pins)

        nf = size(A, 2); nc = size(Ac, 2)
        sk = nf - K.matrix_rank(A)[1]; skc = nc - K.matrix_rank(Ac)[1]
        dk = nf - rank_dense(A, 1e-10); dkc = nc - rank_dense(Ac, 1e-10)
        sd = sk - (tc.n_dangling + skc); dd = dk - (tc.n_dangling + dkc)

        gf = sv_gap(A, nf - dk); gc = sv_gap(Ac, nc - dkc)

        print(csv, id, ",", g.kind, "/", src, ",", K.n_faces(m), ",", nf, ",", tc.n_dangling, ",",
              sk, ",", skc, ",", sd, ",", dk, ",", dkc, ",", dd, ",",
              fmt_b(dd == 0), ",", fmt_g(gf), ",", fmt_g(gc), ",", fmt_g(s(t)), "\n")
        flush(csv)
        n_tried += 1
        dd == 0 ? (n_dense_ok += 1) : (n_still_bad += 1)
        both("  id $id ($(g.kind), F=$(K.n_faces(m)), X=$src): SparseQR defect $sd -> dense defect $dd" *
             (dd == 0 ? "  IDENTITY HOLDS" : "  STILL VIOLATED") * "  [smallest kept sigma / sigma_0 = " *
             sci(gf) * " full, " * sci(gc) * " core]\n")
    end
    both("\nviolators rechecked with dense QR : $n_tried\n")
    both("  identity restored                : $n_dense_ok\n")
    both("  still violated                   : $n_still_bad\n")
    both("  skipped (too large for dense QR) : $n_skipped\n")
    close(sm)
    close(csv)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
