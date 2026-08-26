using Kirigami, LinearAlgebra, JSON; const K = Kirigami
include("/Users/emredayangac/Documents/kirigami-julia/Kirigami/apps/kill_common.jl")
const S = length(ARGS) >= 1 ? ARGS[1] : pwd()   # scratch dir for the JSON outputs
function chart(rows, a)
    keep = trues(length(a))
    for i in 1:length(a)-1
        (rows[i].kind == K.ConeCornerG1 && rows[i].convex && rows[i+1].kind == K.ConeCornerG2 && rows[i+1].index == rows[i].index) || continue
        a[i] >= a[i+1] ? (keep[i+1] = false) : (keep[i] = false)
    end
    keep
end
inp = JSON.parsefile(joinpath(S, ARGS[1])); out = []
for r in inp
    m = K.mesh_from_json(r["mesh"]); c = K.make_cut(m)
    X = haskey(r, "X") ? [K.Vec2(p[1], p[2]) for p in r["X"]] : m.X
    gr = K.build_hinge_graph(c); pins = K.pins_flat(c, gr, X); fb = K.flex_basis(gr, pins)
    cs = K.cone_system(c, X); M0 = Matrix(cs.A * fb.N)
    rn = [norm(M0[i, :]) for i in 1:size(M0, 1)]
    okr = findall(>=(1e-12), rn); noise = size(M0,1) - length(okr)
    M = copy(M0); K.normalise_rows!(M); z = fb.N' * K.sigma_flex(c, X); a = M * z
    keep = chart(cs.rows, a); sigma_chart = minimum(a[keep]) / norm(z)
    Mc = M0[okr, :]; K.normalise_rows!(Mc); ac = Mc * z; keepc = chart(cs.rows[okr], ac); scc = minimum(ac[keepc]) / norm(z)
    opt = K.ExpansiveConeOptions(); opt.lp.dual_iters = 600; rep = K.expansive_cone(c, X, opt)
    Mk = M[findall(keep), :]; r600 = K.cone_lp(Mk, K.ConeLPOptions(dual_iters = 600)); r20k = K.cone_lp(Mk)
    Mkc = Mc[findall(keepc), :]; rc600 = K.cone_lp(Mkc, K.ConeLPOptions(dual_iters = 600)); rc20k = K.cone_lp(Mkc)
    println("$(r["kind"]): rows=$(size(M0,1)) noise=$noise in-chart=$(count(i->keep[i], findall(<(1e-12), rn))) sigma_chart=$sigma_chart clean=$scc rep: margin=$(rep.margin_l2) dual=$(rep.dual_bound) pass=$(rep.pass_feasible)")
    println("    chart LP 600: $(r600.margin_l2) / $(r600.dual_bound) | 20k: $(r20k.margin_l2) / $(r20k.dual_bound)")
    println("    clean LP 600: $(rc600.margin_l2) / $(rc600.dual_bound) | 20k: $(rc20k.margin_l2) / $(rc20k.dual_bound)")
    push!(out, Dict("kind"=>r["kind"], "noise"=>noise, "sigma_chart"=>sigma_chart, "sigma_chart_clean"=>scc, "rep_margin"=>rep.margin_l2, "rep_dual"=>rep.dual_bound,
        "chart600"=>[r600.margin_l2, r600.dual_bound], "chart20k"=>[r20k.margin_l2, r20k.dual_bound], "clean600"=>[rc600.margin_l2, rc600.dual_bound], "clean20k"=>[rc20k.margin_l2, rc20k.dual_bound],
        "gram_clean"=>Mc*Mc'))
end
open(joinpath(S, ARGS[2]), "w") do io; JSON.print(io, out); end
