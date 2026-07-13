# F23 -- is the authors' extra constraint row for boundary-touching split-forest
# components an over-constraint? Port of code/apps/kill_f23.cpp.
#
# STATE.md F23: "the authors' code adds a constraint row for split-forest components
# that TOUCH THE BOUNDARY (we drop them: notches). Their X0 satisfies our system 8/8;
# our X0 has residual 1.73 (hexagons_auto) / 1.00 (snub) in theirs; our null space is
# 1 dim larger on those two. Orchestrator reasoning: faces around a boundary-touching
# component form a PATH in Gamma (the exterior breaks the cycle), so no closure
# constraint is needed -> their extra row is over-constraint (harmless, shrinks their
# space). To be verified by FK path-independence on our X0."
#
# The test, which is the definition of "no closure constraint is needed": deploy OUR X0
# by forward kinematics from EVERY seed face and from randomised BFS orders. If the
# resulting M' vertex positions agree (after the rigid alignment that different seeds
# trivially induce) and every hinge opens by exactly theta, then the deployment is
# path-independent and our X0 needs no further constraint -- so their extra row removes
# admissible designs without buying closure.
#
#   julia --project=Kirigami Kirigami/apps/kill_f23.jl [--out DIR] [--limit K] [--regenerate]
#
# Population: the two split-bearing reference cases hexagons_auto / snub_square_33434
# (frozen reference_cases_8; `--regenerate` rebuilds). `--limit K` keeps only the first K
# of them. Outputs go to results/kill/f23_julia/ (the C++ wrote results/kill/f23/).
include(joinpath(@__DIR__, "kill_common.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

function main(args::Vector{String})
    args, regenerate = take_regenerate_flag(args)
    outdir = joinpath(REPO, "results", "kill", "f23_julia")
    limit = typemax(Int)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--out" && i < length(args); outdir = args[i+1]; i += 2
        elseif a == "--limit" && i < length(args); limit = arg_i(args[i+1]); i += 2
        else; i += 1
        end
    end
    mkpath(outdir)
    csv = open(joinpath(outdir, "f23.csv"), "w")
    print(csv, "name,F,n_split,n_hinge,dim_null,H,rank_L,theta,n_seeds,n_orders,",
          "worst_seed_disagreement,worst_order_disagreement,worst_hinge_angle_error,",
          "worst_bfs_mismatch,worst_split_parallel_error\n")
    sm = open(joinpath(outdir, "summary.txt"), "w")
    both(str) = (print(str); print(sm, str))
    both("F23: forward-kinematics path independence of OUR X0\n")
    both("(if FK is path independent then no extra closure row is needed, and the\n")
    both(" authors' boundary-component row is over-constraint)\n\n")

    global_worst_seed = 0.0; global_worst_angle = 0.0; global_worst_parallel = 0.0
    nrun = 0

    for rc in reference_cases(regenerate = regenerate)
        (rc.name == "hexagons_auto" || rc.name == "snub_square_33434") || continue
        nrun >= limit && break
        nrun += 1
        m = rc.mesh
        K.build_topology!(m)
        c = K.make_cut(m)
        hs = K.holes_partition(c)
        sr = K.solve_system(K.assemble_system(c, hs, m.X, K.Fixed), m.X)
        if !sr.projection_ok
            both(rc.name * ": solve failed\n")
            continue
        end
        X = K.matrix_to_points(sr.X0)
        both(rc.name * ": F=" * string(K.n_faces(m)) * " n_split=" *
             string(K.n_split(c)) * " n_hinge=" * string(K.n_hinge(c)) *
             " H=" * string(sr.H) * " rank(L)=" * string(sr.rank_L) *
             " dim_null=" * string(sr.dim_null) * "\n")
        both("  hole-closure residual of our X0 (Eq. (2)) : " *
             sci(K.hole_residuals(c, X, hs).max_norm) * "\n")

        for theta in (0.3, 0.7, 1.2, 2.0)
            d0 = K.deploy(c, X, theta, 1)
            worst_seed = 0.0; worst_order = 0.0; worst_angle = 0.0; worst_bfs = d0.max_mismatch
            worst_par = 0.0
            n_seeds = 0; n_orders = 0

            # (a) every seed face
            for sf in 1:K.n_faces(m)
                d = K.deploy(c, X, theta, sf)
                length(d.Y) != length(d0.Y) && continue
                n_seeds += 1
                worst_bfs = max(worst_bfs, d.max_mismatch)
                # Different seeds place the structure in different rigid frames; compare after
                # the best rigid alignment, which is what "the same deployed shape" means.
                worst_seed = max(worst_seed, K.rigid_align_residual(d.Y, d0.Y))
            end

            # (b) randomised BFS visit orders from seed 0
            rng = K.MT19937(UInt32(20260903))
            for _ in 1:20
                order = collect(1:K.n_faces(m))
                K.shuffle!(@view(order[2:end]), rng)
                d = K.deploy_with_order(c, X, theta, order)
                length(d.Y) != length(d0.Y) && continue
                n_orders += 1
                worst_order = max(worst_order, K.rigid_align_residual(d.Y, d0.Y))
                worst_bfs = max(worst_bfs, d.max_mismatch)
            end

            # (c) every hinge opens by exactly theta: the two duplicates of the hinge edge's
            # TARGET, seen from the two faces, subtend angle theta at the hinge (source) vertex.
            for ei in c.hinge_edges
                ed = m.edges[ei]
                ed.n_faces != 2 && continue
                f = m.half_edges[ed.he[1]].face; g = m.half_edges[ed.he[2]].face
                src = c.hinge_dir[ei].src; dst = c.hinge_dir[ei].dst
                pf = K.prime_vertex(c, f, src); pg = K.prime_vertex(c, g, src)
                qf = K.prime_vertex(c, f, dst); qg = K.prime_vertex(c, g, dst)
                # the hinge point itself must be shared
                worst_angle = max(worst_angle, K._norm2(d0.Y[pf] - d0.Y[pg]) /
                                               max(1e-300, K._norm2(d0.Y[qf] - d0.Y[pf])))
                u = d0.Y[qf] - d0.Y[pf]; v = d0.Y[qg] - d0.Y[pg]
                ang = abs(K.libm_atan2(K._det2(u, v), K._dot2(u, v)))
                worst_angle = max(worst_angle, abs(ang - theta))
            end

            # (d) Remark A.4: the two duplicates of a split edge stay parallel and equal.
            for ei in c.split_edges
                ed = m.edges[ei]
                ed.n_faces != 2 && continue
                f = m.half_edges[ed.he[1]].face; g = m.half_edges[ed.he[2]].face
                a_ = ed.key.a; b_ = ed.key.b
                a1 = K.prime_vertex(c, f, a_); b1 = K.prime_vertex(c, f, b_)
                a2 = K.prime_vertex(c, g, a_); b2 = K.prime_vertex(c, g, b_)
                e1 = d0.Y[b1] - d0.Y[a1]; e2 = d0.Y[b2] - d0.Y[a2]
                worst_par = max(worst_par, K._norm2(e1 - e2) / max(1e-300, K._norm2(e1)))
            end

            print(csv, rc.name, ",", K.n_faces(m), ",", K.n_split(c), ",", K.n_hinge(c), ",",
                  sr.dim_null, ",", sr.H, ",", sr.rank_L, ",", cpp_g(theta), ",", n_seeds,
                  ",", n_orders, ",", cpp_g(worst_seed), ",", cpp_g(worst_order), ",", cpp_g(worst_angle),
                  ",", cpp_g(worst_bfs), ",", cpp_g(worst_par), "\n")
            flush(csv)
            both("  theta=" * fx(theta, 2) * ": seeds=" * string(n_seeds) * " orders=" *
                 string(n_orders) * "  worst |Y(seed s) - Y(seed 0)| after rigid align = " *
                 sci(worst_seed) * ", worst over BFS orders = " * sci(worst_order) *
                 ", worst hinge-angle error = " * sci(worst_angle) *
                 ", worst non-tree BFS mismatch = " * sci(worst_bfs) *
                 ", worst split-duplicate parallel error = " * sci(worst_par) * "\n")
            global_worst_seed = max(global_worst_seed, max(worst_seed, worst_order))
            global_worst_angle = max(global_worst_angle, worst_angle)
            global_worst_parallel = max(global_worst_parallel, worst_par)
        end
        both("\n")
    end
    both("worst seed/order disagreement over both patterns : " * sci(global_worst_seed) * "\n")
    both("worst hinge-angle error                          : " * sci(global_worst_angle) * "\n")
    both("worst split-duplicate parallel error             : " * sci(global_worst_parallel) * "\n")
    pass = global_worst_seed < 1e-9 && global_worst_angle < 1e-9
    both("VERDICT (path independent to 1e-9): " * (pass ? "PASS" : "FAIL") * "\n")
    both("=> the authors' extra boundary-component row is " *
         (pass ? "OVER-CONSTRAINT: our X0 already closes every hole and deploys\n" *
                 "   path-independently without it, so the row only removes designs.\n" :
                 "NOT shown to be over-constraint by this test.\n"))
    close(csv)
    close(sm)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
