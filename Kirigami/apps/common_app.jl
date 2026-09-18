# apps/common_app.jl -- shared helpers of the kiri_* command-line apps. A plain included
# file (not a module): every app does
#     include(joinpath(@__DIR__, "common_app.jl"))
# and then talks to the package through `K.<name>` (`Kirigami` exports nothing).

using Kirigami
using LinearAlgebra
using Printf
import JSON
const K = Kirigami
const Vec2 = K.Vec2

# One BLAS thread per process (STATE.md U10, docs/NUMERICS.md "BLAS threads"). The
# optimiser objectives are hundreds of small gemvs per iteration; with OpenBLAS's default
# thread count, 12 concurrent shards oversubscribe the machine and each gemv stalls on
# thread synchronisation (measured 330 s vs 3 s per B4 row). The thread count also fixes
# the gemv reduction order, so it is part of the results' reproducibility.
BLAS.set_num_threads(1)

# NaN / +-Inf are written as `null` (JSON.jl refuses to write them), and Dict keys come
# out sorted, so files are byte-stable across runs.
_json_clean(x::AbstractFloat) = isfinite(x) ? x : nothing
_json_clean(x::AbstractDict) = Dict{String,Any}(string(k) => _json_clean(v) for (k, v) in x)
_json_clean(x::AbstractVector) = Any[_json_clean(v) for v in x]
_json_clean(x::Tuple) = Any[_json_clean(v) for v in x]
_json_clean(x) = x

"""`write_json(j, path)`: `j.dump(2)` plus a trailing newline, creating the parent directory."""
function write_json(j, path::AbstractString)
    d = dirname(path)
    isempty(d) || mkpath(d)
    open(path, "w") do io
        JSON.print(io, _json_clean(j), 2)
        write(io, "\n")
    end
    return nothing
end

"""Dump of a deployed kirigami structure for scripts/plot_embedding.jl (0-based M' faces
and 0-based edge ids in `holes`, the file format's numbering)."""
function deployment_json(c::K.CutStructure, X::Vector{Vec2}, theta::Real)
    d = K.deploy(c, X, theta)
    j = Dict{String,Any}()
    j["theta"] = Float64(theta)
    j["max_mismatch"] = d.max_mismatch
    j["vertices"] = [[y[1], y[2]] for y in d.Y]
    # internal 1-based -> 0-based in the file
    j["faces"] = [[v - 1 for v in f] for f in c.prime_faces]
    j["holes"] = [[e - 1 for e in h] for h in K.holes_geometric_cycles(c, d.Y)]
    return j
end

"""The README contract fields of a mesh: vertices, faces (0-based in the file), orientation."""
function mesh_json(m::K.Mesh)
    j = Dict{String,Any}()
    j["vertices"] = [[x[1], x[2]] for x in m.X]
    # internal 1-based -> 0-based in the file
    j["faces"] = [[i - 1 for i in f] for f in m.faces]
    isempty(m.sigma) || (j["orientation"] = copy(m.sigma))
    return j
end

# ---- argument parsing shared by the apps ------------------------------------------
# Numeric `argv` conversions; junk is an error.
arg_f(s::AbstractString) = parse(Float64, s)
arg_i(s::AbstractString) = parse(Int, s)
arg_u(s::AbstractString) = UInt32(parse(UInt, s) % UInt32)

# `checkerboard`, shared by kiri_gen / kiri_reference / kiri_sweep: the 2-colouring of
# the dual graph by DFS.
checkerboard(m::K.Mesh) = K.checkerboard(m)

# Number formatting of the CSV / text outputs.
# Non-finite values print as "inf", "-inf", "nan".
_fmt_nonfinite(v::Float64) = isnan(v) ? "nan" : (v > 0 ? "inf" : "-inf")
# 6 significant digits (%g).
fmt_g(v::Real) = isfinite(v) ? @sprintf("%g", Float64(v)) : _fmt_nonfinite(Float64(v))
# scientific / fixed with `prec` digits after the point.
sci(v::Real, prec::Int = 3) = isfinite(v) ? Printf.format(Printf.Format("%.$(prec)e"), Float64(v)) : _fmt_nonfinite(Float64(v))
fx(v::Real, prec::Int = 4) = isfinite(v) ? Printf.format(Printf.Format("%.$(prec)f"), Float64(v)) : _fmt_nonfinite(Float64(v))
# booleans as 0 / 1.
fmt_b(b::Bool) = b ? "1" : "0"
