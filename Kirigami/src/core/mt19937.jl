# core/mt19937.jl -- bit-exact port of std::mt19937 and of the Apple libc++ (LLVM 18,
# _LIBCPP_VERSION 180100, Xcode MacOSX.sdk) distributions the C++ code drew from:
# generate_canonical<double,53>, uniform_real_distribution<double>,
# uniform_int_distribution<int>, std::shuffle.  Deliberately independent of Random so the
# stream is exactly the one the frozen corpora were produced with.

const MT_N = 624
const MT_M = 397
const MT_MATRIX_A = 0x9908b0df
const MT_UPPER = 0x80000000
const MT_LOWER = 0x7fffffff
const MT_DEFAULT_SEED = 5489

"""32-bit Mersenne Twister state (`std::mt19937`)."""
mutable struct MT19937
    mt::Vector{UInt32}
    mti::Int            # 1-based index of the next word to return; MT_N+1 == needs twist
end

"""`std::mt19937(seed)`: init_genrand recurrence x_i = 1812433253 (x_{i-1} ⊻ (x_{i-1} >> 30)) + i."""
function MT19937(seed::Integer = MT_DEFAULT_SEED)
    mt = Vector{UInt32}(undef, MT_N)
    mt[1] = UInt32(seed & 0xffffffff)
    for i in 2:MT_N
        prev = mt[i-1]
        mt[i] = UInt32(1812433253) * (prev ⊻ (prev >> 30)) + UInt32(i - 1)
    end
    MT19937(mt, MT_N + 1)
end

function _twist!(rng::MT19937)
    mt = rng.mt
    @inbounds for k in 1:MT_N
        y = (mt[k] & MT_UPPER) | (mt[mod1(k + 1, MT_N)] & MT_LOWER)
        v = mt[mod1(k + MT_M, MT_N)] ⊻ (y >> 1)
        mt[k] = (y & 0x1) == 0x1 ? v ⊻ MT_MATRIX_A : v
    end
    rng.mti = 1
    nothing
end

"""`std::mt19937::operator()`: next 32-bit output (min 0, max 2^32-1)."""
function next_u32(rng::MT19937)::UInt32
    rng.mti > MT_N && _twist!(rng)
    y = @inbounds rng.mt[rng.mti]
    rng.mti += 1
    y ⊻= y >> 11
    y ⊻= (y << 7) & 0x9d2c5680
    y ⊻= (y << 15) & 0xefc60000
    y ⊻= y >> 18
    return y
end

const MT_RANGE = 4294967296.0   # R = max - min + 1 as double

"""libc++ `std::generate_canonical<double,53>(mt19937)`: k = ⌈53/32⌉ = 2 draws,
S = d0 + d1·R, result S / R².  This libc++ has no clamp of 1.0 (LWG 2524 not applied), so
the result may equal 1.0 exactly when S rounds up to 2^64; replicated on purpose."""
function generate_canonical53(rng::MT19937)::Float64
    sp = Float64(next_u32(rng))
    base = MT_RANGE
    sp += Float64(next_u32(rng)) * base
    base *= MT_RANGE
    return sp / base
end

"""libc++ `uniform_real_distribution<double>(a,b)(rng)` = (b-a)·generate_canonical + a.
Evaluated with `fma`: clang on arm64 (default -ffp-contract=on, as in the C++ build) fuses
this multiply-add, and the unfused form differs from the C++ stream by 1 ulp ~40% of the time."""
uniform_real(rng::MT19937, a::Real, b::Real)::Float64 =
    fma(Float64(b) - Float64(a), generate_canonical53(rng), Float64(a))

# libc++ __independent_bits_engine<mt19937, W> with W the working unsigned type
# (UInt32 for uniform_int_distribution<int>, UInt64 for <ptrdiff_t> as used by std::shuffle).
# _Rp = 2^32 in W: wraps to 0 for UInt32 (-> __eval(false_type): one draw masked), stays 2^32
# for UInt64 (-> __eval(true_type) with y0 = y1 = 2^32, so the inner rejections never fire).
# __log2<W, _Rp>::value == 32 in both cases.
function _independent_bits(rng::MT19937, ::Type{W}, w::Int)::W where {W<:Union{UInt32,UInt64}}
    Rp = W === UInt32 ? UInt32(0) : UInt64(4294967296)
    m = 32
    WDt = 8 * sizeof(W)
    EDt = 32
    n = w ÷ m + (w % m != 0 ? 1 : 0)
    w0 = w ÷ n
    y0 = Rp == 0 ? Rp : (w0 < WDt ? (Rp >> w0) << w0 : W(0))
    if Rp - y0 > y0 ÷ n
        n += 1
        w0 = w ÷ n
        y0 = w0 < WDt ? (Rp >> w0) << w0 : W(0)
    end
    n0 = n - w % n
    y1 = w0 < WDt - 1 ? (Rp >> (w0 + 1)) << (w0 + 1) : W(0)
    mask0 = w0 > 0 ? typemax(UInt32) >> (EDt - w0) : UInt32(0)
    mask1 = w0 < EDt - 1 ? typemax(UInt32) >> (EDt - (w0 + 1)) : typemax(UInt32)
    if Rp == 0                                   # __eval(false_type)
        return W(next_u32(rng) & mask0)
    end
    w_rt = WDt                                   # __eval(true_type)
    sp = W(0)
    for _ in 1:n0
        u = next_u32(rng)
        while u >= y0
            u = next_u32(rng)
        end
        sp = w0 < w_rt ? sp << w0 : W(0)
        sp += W(u & mask0)
    end
    for _ in (n0+1):n
        u = next_u32(rng)
        while u >= y1
            u = next_u32(rng)
        end
        sp = w0 < w_rt - 1 ? sp << (w0 + 1) : W(0)
        sp += W(u & mask1)
    end
    return sp
end

# libc++ uniform_int_distribution<T>::operator() with working type W (unsigned of T's width,
# at least 32 bits).  Returns the offset u in [0, b-a] as W.
function _uniform_int_offset(rng::MT19937, ::Type{W}, a::Integer, b::Integer)::W where {W<:Union{UInt32,UInt64}}
    rp = W(b % W) - W(a % W) + W(1)          # wraps like the C++ unsigned arithmetic
    rp == 1 && return W(0)
    dt = 8 * sizeof(W)
    rp == 0 && return _independent_bits(rng, W, dt)
    w = dt - leading_zeros(rp) - 1
    if (rp & (typemax(W) >> (dt - w))) != 0
        w += 1
    end
    u = _independent_bits(rng, W, w)
    while u >= rp
        u = _independent_bits(rng, W, w)
    end
    return u
end

"""libc++ `uniform_int_distribution<int>(a,b)(rng)`: rejection sampling on the smallest
power-of-two window covering b-a+1, one 32-bit draw per attempt.  Requires a <= b."""
function uniform_int(rng::MT19937, a::Integer, b::Integer)::Int
    a <= b || throw(ArgumentError("uniform_int: need a <= b"))
    u = _uniform_int_offset(rng, UInt32, Int32(a), Int32(b))
    return Int(reinterpret(Int32, u + reinterpret(UInt32, Int32(a))))   # static_cast<int>(u + a)
end

"""libc++ `std::shuffle(v.begin(), v.end(), rng)`: for i = 1..n-1 draw
j ~ uniform_int_distribution<ptrdiff_t>(0, n-i) and swap v[i], v[i+j] (skipping j == 0)."""
function shuffle!(v::AbstractVector, rng::MT19937)
    d = length(v)
    d > 1 || return v
    d -= 1
    i = firstindex(v)
    while d >= 1
        j = Int(_uniform_int_offset(rng, UInt64, 0, d))
        if j != 0
            v[i], v[i+j] = v[i+j], v[i]
        end
        i += 1
        d -= 1
    end
    return v
end

"""libc++ `std::normal_distribution<double>(mean, stddev)`: Marsaglia polar method; the
second variate of each pair is cached in `v` (`__v_hot_`), so the object carries state and
one instance must live for the whole C++ distribution's lifetime."""
mutable struct NormalDist
    mean::Float64
    stddev::Float64
    v::Float64
    v_hot::Bool
end
NormalDist(mean::Real = 0.0, stddev::Real = 1.0) = NormalDist(Float64(mean), Float64(stddev), 0.0, false)

"""`normal_distribution::operator()(rng)`.  Products/sums use `fma` where clang -O2 on arm64
contracts them (`u*u + v*v` -> fma(u,u,v*v), `up*stddev + mean`).  `log` is Julia's, not
Apple libm's; it agreed on every reference draw but is not guaranteed to be bit-identical."""
function normal(d::NormalDist, rng::MT19937)::Float64
    if d.v_hot
        d.v_hot = false
        up = d.v
    else
        u = v = s = 0.0
        while true
            u = uniform_real(rng, -1.0, 1.0)
            v = uniform_real(rng, -1.0, 1.0)
            s = fma(u, u, v * v)
            (s > 1 || s == 0) || break
        end
        fp = sqrt(-2 * log(s) / s)
        d.v = v * fp
        d.v_hot = true
        up = u * fp
    end
    return fma(up, d.stddev, d.mean)
end
