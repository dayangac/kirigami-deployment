# core/optimize.jl -- a small self-contained L-BFGS with backtracking line search.
# No external optimization library is used (see specs/common_preamble.md).

Base.@kwdef mutable struct LbfgsOptions
    max_iter::Int = 500
    history::Int = 10
    grad_tol::Float64 = 1e-9
    step_tol::Float64 = 1e-14
    c1::Float64 = 1e-4      # Armijo constant
    shrink::Float64 = 0.5   # backtracking factor
    max_ls::Int = 40
end

mutable struct LbfgsResult
    x::Vector{Float64}
    f::Float64
    grad_norm::Float64
    iterations::Int
    converged::Bool
end

"""
    lbfgs_minimize(fg, x0, opt = LbfgsOptions()) -> LbfgsResult

`fg(x, g)` returns f(x) and writes the gradient into `g`.
"""
function lbfgs_minimize(fg, x0::Vector{Float64}, opt::LbfgsOptions = LbfgsOptions())
    x = copy(x0)
    n = length(x)
    g = zeros(n)
    g_new = zeros(n)
    f = fg(x, g)
    S = Vector{Float64}[]
    Y = Vector{Float64}[]
    RH = Float64[]
    res = LbfgsResult(x, f, 0.0, 0, false)

    for it in 0:opt.max_iter-1
        res.iterations = it
        if norm(g) < opt.grad_tol
            res.converged = true
            break
        end
        # two-loop recursion
        q = copy(g)
        k = length(S)
        alpha = zeros(k)
        for i in k:-1:1
            alpha[i] = RH[i] * dot(S[i], q)
            q -= alpha[i] * Y[i]
        end
        scale = 1.0
        if k > 0
            scale = dot(S[k], Y[k]) / max(dot(Y[k], Y[k]), 1e-300)
        end
        q *= scale
        for i in 1:k
            beta = RH[i] * dot(Y[i], q)
            q += S[i] * (alpha[i] - beta)
        end
        d = -q
        dg = dot(d, g)
        if !(dg < 0)
            d = -g
            dg = dot(d, g)
        end
        step = (it == 0 && k == 0) ? min(1.0, 1.0 / max(norm(g), 1e-12)) : 1.0
        f_new = f
        ok = false
        for ls in 1:opt.max_ls
            x_new = x + step * d
            f_new = fg(x_new, g_new)
            if isfinite(f_new) && f_new <= f + opt.c1 * step * dg
                s = x_new - x
                y = g_new - g
                sy = dot(s, y)
                if sy > 1e-14
                    push!(S, s)
                    push!(Y, y)
                    push!(RH, 1.0 / sy)
                    if length(S) > opt.history
                        popfirst!(S)
                        popfirst!(Y)
                        popfirst!(RH)
                    end
                end
                x = x_new
                f = f_new
                g, g_new = g_new, g   # g takes the new gradient; the old buffer is reused
                ok = true
                break
            end
            step *= opt.shrink
            step * norm(d) < opt.step_tol && break
        end
        ok || break
    end
    res.x = x
    res.f = f
    res.grad_norm = norm(g)
    return res
end
