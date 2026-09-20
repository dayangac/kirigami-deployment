# teaser_cat.jl -- the tech report's head figure: a cat-head planar graph, the Eq. (6)
# projection that jams at 0+, and the range-maximising embedding that deploys.
#
#   julia --project=Kirigami Kirigami/scripts/teaser_cat.jl [--seed S] [--out DIR]
#
# Writes DIR/teaser_cat.json (polygons of every panel) for scripts/plot_teaser.py.
using Kirigami
using LinearAlgebra, JSON, Printf
const K = Kirigami
const Vec2 = K.Vec2

seed = 7; outdir = "results/final/figures"
let a = ARGS, i = 1
    while i <= length(a)
        a[i] == "--seed" && (seed = parse(Int, a[i+1]); i += 2; continue)
        a[i] == "--out" && (outdir = a[i+1]; i += 2; continue)
        i += 1
    end
end
mkpath(outdir)

# ---- the silhouette: a cat head with two ears, CCW ---------------------------------
function cat_outline(n_arc::Int = 48)
    P = Vec2[]
    # head: circle of radius 1 from the right ear base round the chin to the left ear base
    a0, a1 = deg2rad(35.0), deg2rad(145.0)          # ear bases on the circle
    for k in 0:n_arc
        a = a0 - (2pi - (a1 - a0)) * k / n_arc      # clockwise from a0 down to a1 - 2pi
        push!(P, Vec2(cos(a), sin(a)))
    end
    # left ear: base at a1, tip, back to the circle at 100 deg
    push!(P, Vec2(-1.05, 1.55))
    push!(P, Vec2(cos(deg2rad(100.0)), sin(deg2rad(100.0))))
    # right ear
    push!(P, Vec2(cos(deg2rad(80.0)), sin(deg2rad(80.0))))
    push!(P, Vec2(1.05, 1.55))
    return P
end

function point_in_poly(P::Vector{Vec2}, q::Vec2)
    inside = false; n = length(P)
    for i in 1:n
        a = P[i]; b = P[mod1(i + 1, n)]
        if (a[2] > q[2]) != (b[2] > q[2])
            x = a[1] + (q[2] - a[2]) * (b[1] - a[1]) / (b[2] - a[2])
            x > q[1] && (inside = !inside)
        end
    end
    return inside
end

# resample the outline at spacing h, then a jittered hex lattice inside
function sample_points(P::Vector{Vec2}, h::Float64, rng::K.MT19937)
    pts = Vec2[]
    n = length(P)
    for i in 1:n
        a = P[i]; b = P[mod1(i + 1, n)]
        L = norm(b - a); m = max(1, round(Int, L / h))
        for k in 0:m-1
            push!(pts, a + (b - a) * (k / m))
        end
    end
    nb = length(pts)
    dy = h * sqrt(3) / 2
    y = -1.05
    row = 0
    while y < 1.6
        x = -1.15 + (isodd(row) ? h / 2 : 0.0)
        while x < 1.15
            j = Vec2(K.uniform_real(rng, -0.22h, 0.22h), K.uniform_real(rng, -0.22h, 0.22h))
            q = Vec2(x, y) + j
            # keep a margin from the outline so boundary triangles stay well shaped
            if point_in_poly(P, q) && minimum(norm(q - p) for p in pts[1:nb]) > 0.55h
                push!(pts, q)
            end
            x += h
        end
        y += dy; row += 1
    end
    return pts
end

function cat_mesh(seed::Int; h::Float64 = 0.2)
    rng = K.MT19937(seed)
    P = cat_outline()
    pts = sample_points(P, h, rng)
    tris = K.delaunay_triangles(pts)
    keep = [t for t in tris if point_in_poly(P, (pts[t[1]] + pts[t[2]] + pts[t[3]]) / 3)]
    # greedy quad merging across shuffled interior edges when the union is convex
    tri = K.largest_component(K.mesh_from_polygons([[pts[t[1]], pts[t[2]], pts[t[3]]] for t in keep]))
    used = falses(K.n_faces(tri))
    eidx = collect(1:K.n_edges(tri)); K.shuffle!(eidx, rng)
    polys = Vector{Vec2}[]
    for e in eidx
        ed = tri.edges[e]
        ed.n_faces != 2 && continue
        f0 = tri.half_edges[ed.he[1]].face; f1 = tri.half_edges[ed.he[2]].face
        (used[f0] || used[f1]) && continue
        a1 = only(v for v in tri.faces[f0] if v != ed.key.a && v != ed.key.b)
        a2 = only(v for v in tri.faces[f1] if v != ed.key.a && v != ed.key.b)
        q = [ed.key.a, a1, ed.key.b, a2]
        K._quad_convex(tri.X, q) || (q = [ed.key.a, a2, ed.key.b, a1])
        K._quad_convex(tri.X, q) || continue
        used[f0] = used[f1] = true
        push!(polys, [tri.X[v] for v in q])
    end
    for f in 1:K.n_faces(tri)
        used[f] || push!(polys, [tri.X[v] for v in tri.faces[f]])
    end
    return K.largest_component(K.mesh_from_polygons(polys))
end

# ---- run ------------------------------------------------------------------------------
m = cat_mesh(seed)
sigma = K.orientation_maxcut(m, 20260903 + seed)
@printf("cat: N=%d F=%d\n", K.n_vertices(m), K.n_faces(m))

base = K.design_baseline(m, sigma, Vec2[])
@printf("baseline: ok=%s theta_max=%.4f binding=%s nonconvex corners at X0=%d dim_null=%d\n",
        base.ok, base.ch.theta_max, base.ch.binding, base.n_nonconvex_x0, base.dim_null)
rm_ = K.design_range_max(m, sigma, Vec2[])
d = rm_.design
@printf("range_max: ok=%s theta_max=%.4f certified=%s eps_max=%.4f provenance=%s\n",
        d.ok, d.ch.theta_max, d.ch.certified, d.ch.eps_max, rm_.provenance)

mm = K.Mesh(m.X, m.faces); mm.sigma = sigma; K.build_topology!(mm)
c = K.make_cut(mm)
polys_at(X, θ) = begin
    dp = K.deploy(c, X, θ)
    [[collect(dp.Y[v]) for v in f] for f in c.prime_faces]
end
reflex(X) = begin  # original-vertex indices of reflex corners
    out = Int[]
    for f in m.faces, (k, v) in enumerate(f)
        a = X[f[mod1(k - 1, length(f))]]; b = X[v]; cc = X[f[mod1(k + 1, length(f))]]
        u = b - a; w = cc - b
        u[1] * w[2] - u[2] * w[1] < -1e-12 && push!(out, v)
    end
    unique(out)
end
θb = 0.35                      # opened past the jam so the penetrations are visible
θd = min(d.ch.theta_max / 2, 0.7)   # keep the silhouette recognisable
out = Dict(
    "faces_M" => [[collect(m.X[v]) for v in f] for f in m.faces],
    "sigma" => sigma,
    "X0_faces" => [[collect(base.X0[v]) for v in f] for f in m.faces],
    "X0_reflex" => [collect(base.X0[v]) for v in reflex(base.X0)],
    "X0_deployed" => polys_at(base.X0, θb), "theta_jam" => θb,
    "X0_theta_max" => base.ch.theta_max, "X0_binding" => base.ch.binding,
    # faces penetrating each other at theta_jam (exact polygon overlap), to colour red
    "X0_overlapping" => let dp = K.deploy(c, base.X0, θb), P = [[dp.Y[v] for v in f] for f in c.prime_faces], bad = Set{Int}()
        for i in 1:length(P), j in i+1:length(P)
            K.polygons_overlap(P[i], P[j], 1e-9) && (push!(bad, i); push!(bad, j))
        end
        sort!(collect(bad)) end,
    "ours_flat" => [[collect(d.X[v]) for v in f] for f in m.faces],
    "ours_faces" => [[collect(d.X[v]) for v in f] for f in m.faces],
    "ours_deployed" => polys_at(d.X, θd), "theta_half" => θd,
    "ours_theta_max" => d.ch.theta_max, "ours_eps_max" => d.ch.eps_max,
    "ours_certified" => d.ch.certified, "n_faces" => K.n_faces(m), "n_split" => K.n_split(c),
    "dim_null" => d.dim_null, "seed" => seed)
open(joinpath(outdir, "teaser_cat.json"), "w") do io; JSON.print(io, out); end
println("wrote ", joinpath(outdir, "teaser_cat.json"))
