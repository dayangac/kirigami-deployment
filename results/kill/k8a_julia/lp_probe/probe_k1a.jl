# Julia side of the K8a LP-layer probe: K1a graphs 0..4 at X_ini.
using Kirigami, LinearAlgebra, JSON; const K = Kirigami
include("/Users/emredayangac/Documents/kirigami-julia/Kirigami/apps/kill_common.jl")
const S = length(ARGS) >= 1 ? ARGS[1] : pwd()   # scratch dir for the JSON outputs
function chart(cs, a)
    keep = trues(length(a))
    for i in 1:length(a)-1
        (cs.rows[i].kind == K.ConeCornerG1 && cs.rows[i].convex && cs.rows[i+1].kind == K.ConeCornerG2 && cs.rows[i+1].index == cs.rows[i].index) || continue
        a[i] >= a[i+1] ? (keep[i+1] = false) : (keep[i] = false)
    end
    keep
end
inputs = []; out = []
for g in population("k1a_200"; ids = 0:4)
    m = g.mesh; c = K.make_cut(m); X = m.X
    push!(inputs, Dict("id" => g.id, "kind" => g.kind, "mesh" => K.mesh_to_json(m)))
    gr = K.build_hinge_graph(c); pins = K.pins_flat(c, gr, X); fb = K.flex_basis(gr, pins)
    cs = K.cone_system(c, X); M0 = Matrix(cs.A * fb.N)
    rn = [norm(M0[i, :]) for i in 1:size(M0, 1)]
    noise = findall(<(1e-12), rn); nearnoise = findall(<(1e-8), rn)
    M = copy(M0); K.normalise_rows!(M)
    z = fb.N' * K.sigma_flex(c, X); a = M * z
    keep = chart(cs, a)
    sigma_chart = minimum(a[keep]) / norm(z)
    # clean: drop noise rows before normalising
    okr = findall(>=(1e-12), rn); Mc = M0[okr, :]; K.normalise_rows!(Mc); ac = Mc * z
    csr = K.ConeSystem(); csr.rows = cs.rows[okr]   # partner structure survives when both rows of a pair are dropped together
    keepc = chart(csr, ac); sigma_chart_clean = minimum(ac[keepc]) / norm(z)
    opt = K.ExpansiveConeOptions(); opt.lp.dual_iters = 600
    rep = K.expansive_cone(c, X, opt)
    # LP directly on the chart rows, 600 iters, and again with 20000 (default) iters
    Mk = M[findall(keep), :]
    o600 = K.ConeLPOptions(dual_iters = 600); r600 = K.cone_lp(Mk, o600); r20k = K.cone_lp(Mk)
    Mkc = Mc[findall(keepc), :]; rc600 = K.cone_lp(Mkc, o600); rc20k = K.cone_lp(Mkc)
    println("$(g.kind)_$(g.id): rows=$(size(M0,1)) noise(<1e-12)=$(length(noise)) (<1e-8)=$(length(nearnoise)) in-chart noise=$(count(i->keep[i], noise)) sigma_chart=$sigma_chart clean=$sigma_chart_clean  rep: margin=$(rep.margin_l2) dual=$(rep.dual_bound) pass=$(rep.pass_feasible)")
    println("    chart LP 600: margin=$(r600.margin_l2) dual=$(r600.dual_bound) gap=$(r600.gap) | 20k: margin=$(r20k.margin_l2) dual=$(r20k.dual_bound) gap=$(r20k.gap)")
    println("    clean LP 600: margin=$(rc600.margin_l2) dual=$(rc600.dual_bound) gap=$(rc600.gap) | 20k: margin=$(rc20k.margin_l2) dual=$(rc20k.dual_bound) gap=$(rc20k.gap)")
    push!(out, Dict("id" => g.id, "kind" => g.kind, "n_rows" => size(M0,1), "noise" => length(noise), "noise_in_chart" => count(i->keep[i], noise),
        "sigma_chart" => sigma_chart, "sigma_chart_clean" => sigma_chart_clean, "rep_margin" => rep.margin_l2, "rep_dual" => rep.dual_bound,
        "rep_dual_max" => rep.dual_bound_max, "rep_n_active" => rep.n_active, "rep_n_dual_support" => rep.n_dual_support, "rep_pass" => rep.pass_feasible,
        "chart600" => [r600.margin_l2, r600.dual_bound], "chart20k" => [r20k.margin_l2, r20k.dual_bound],
        "clean600" => [rc600.margin_l2, rc600.dual_bound], "clean20k" => [rc20k.margin_l2, rc20k.dual_bound],
        "gram" => M * M', "a_sigma" => a, "row_norms" => rn))
end
open(joinpath(S, "k1a_inputs.json"), "w") do io; JSON.print(io, inputs); end
open(joinpath(S, "julia_probe.json"), "w") do io; JSON.print(io, out); end
