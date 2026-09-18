# core/mt19937.jl -- bit-exact MT19937 (std::mt19937) with the libc++-compatible (LLVM 18,
# _LIBCPP_VERSION 180100) distributions the frozen corpora were drawn with:
# generate_canonical<double,53>, uniform_real_distribution<double>,
# uniform_int_distribution<int>, std::shuffle.  Deliberately independent of Random so that
# seeds reproduce the archived design populations (`data/corpus/mt19937_vectors.json`).

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
the result may equal 1.0 exactly when S rounds up to 2^64; deliberate, the corpus
depends on it."""
function generate_canonical53(rng::MT19937)::Float64
    sp = Float64(next_u32(rng))
    base = MT_RANGE
    sp += Float64(next_u32(rng)) * base
    base *= MT_RANGE
    return sp / base
end

"""libc++ `uniform_real_distribution<double>(a,b)(rng)` = (b-a)·generate_canonical + a.
Evaluated with `fma` so the result matches the reference RNG stream / the frozen corpora
bit for bit; the unfused form differs by 1 ulp ~40% of the time (docs/NUMERICS.md)."""
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
    rp = W(b % W) - W(a % W) + W(1)          # wrapping unsigned arithmetic (deliberate)
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
one instance must live for the whole stream's lifetime."""
mutable struct NormalDist
    mean::Float64
    stddev::Float64
    v::Float64
    v_hot::Bool
end
NormalDist(mean::Real = 0.0, stddev::Real = 1.0) = NormalDist(Float64(mean), Float64(stddev), 0.0, false)

"""`normal_distribution::operator()(rng)`.  Products/sums use `fma` where the reference
stream has them fused (`u*u + v*v` -> fma(u,u,v*v), `up*stddev + mean`).  `log` is Julia's,
not Apple libm's; it agreed on every reference draw but is not guaranteed to be bit-identical."""
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

# ---------------------------------------------------------------------------- std::sort
# Bit-faithful libc++ 18 `std::sort` (__algorithm/sort.h) for a NON-arithmetic
# value type with an arbitrary comparator: `__introsort<..., _UseBitSetPartition = false>`
# = pdqsort-style introsort with __sort3/4/5 networks, guarded insertion sort on the
# leftmost range and unguarded elsewhere (limit 24), median-of-3 (Tukey ninther above 128),
# `__partition_with_equals_on_right`, the "already partitioned -> try
# __insertion_sort_incomplete" shortcut, and `__partition_with_equals_on_left` when the
# range's predecessor equals the pivot. The result is unstable, so where an archived result
# depends on the ORDER of equivalent elements (e.g. `detect_lattice`'s tie among lattice
# vectors of equal length) only this exact algorithm reproduces it. The heap-sort fallback at
# depth exhaustion is not implemented; it is unreachable on the inputs this code base sorts and
# raises an error rather than silently sorting differently.

@inline function _lc_swap!(v, a, b)
    v[a], v[b] = v[b], v[a]
    return
end

function _lc_sort3!(v, x, y, z, lt)
    if !lt(v[y], v[x])
        lt(v[z], v[y]) || return
        _lc_swap!(v, y, z)
        lt(v[y], v[x]) && _lc_swap!(v, x, y)
        return
    end
    if lt(v[z], v[y])
        _lc_swap!(v, x, z)
        return
    end
    _lc_swap!(v, x, y)
    lt(v[z], v[y]) && _lc_swap!(v, y, z)
    return
end

function _lc_sort4!(v, x1, x2, x3, x4, lt)
    _lc_sort3!(v, x1, x2, x3, lt)
    if lt(v[x4], v[x3])
        _lc_swap!(v, x3, x4)
        if lt(v[x3], v[x2])
            _lc_swap!(v, x2, x3)
            lt(v[x2], v[x1]) && _lc_swap!(v, x1, x2)
        end
    end
end

function _lc_sort5!(v, x1, x2, x3, x4, x5, lt)
    _lc_sort4!(v, x1, x2, x3, x4, lt)
    if lt(v[x5], v[x4])
        _lc_swap!(v, x4, x5)
        if lt(v[x4], v[x3])
            _lc_swap!(v, x3, x4)
            if lt(v[x3], v[x2])
                _lc_swap!(v, x2, x3)
                lt(v[x2], v[x1]) && _lc_swap!(v, x1, x2)
            end
        end
    end
end

# [first, last) half-open, 1-based
function _lc_insertion_sort!(v, first, last, lt)
    first == last && return
    for i in (first + 1):(last - 1)
        j = i - 1
        if lt(v[i], v[j])
            t = v[i]
            k = j
            j = i
            while true
                v[j] = v[k]
                j = k
                (j != first && (k -= 1; lt(t, v[k]))) || break
            end
            v[j] = t
        end
    end
end

# assumes an element before `first` that is <= every element of the range
function _lc_insertion_sort_unguarded!(v, first, last, lt)
    first == last && return
    for i in (first + 1):(last - 1)
        j = i - 1
        if lt(v[i], v[j])
            t = v[i]
            k = j
            j = i
            while true
                v[j] = v[k]
                j = k
                k -= 1
                lt(t, v[k]) || break
            end
            v[j] = t
        end
    end
end

# returns true iff the range is sorted on exit (gives up after 8 insertions)
function _lc_insertion_sort_incomplete!(v, first, last, lt)
    n = last - first
    if n <= 1
        return true
    elseif n == 2
        lt(v[last - 1], v[first]) && _lc_swap!(v, first, last - 1)
        return true
    elseif n == 3
        _lc_sort3!(v, first, first + 1, last - 1, lt)
        return true
    elseif n == 4
        _lc_sort4!(v, first, first + 1, first + 2, last - 1, lt)
        return true
    elseif n == 5
        _lc_sort5!(v, first, first + 1, first + 2, first + 3, last - 1, lt)
        return true
    end
    j = first + 2
    _lc_sort3!(v, first, first + 1, j, lt)
    limit = 8
    count = 0
    i = j + 1
    while i != last
        if lt(v[i], v[j])
            t = v[i]
            k = j
            j = i
            while true
                v[j] = v[k]
                j = k
                (j != first && (k -= 1; lt(t, v[k]))) || break
            end
            v[j] = t
            count += 1
            if count == limit
                return i + 1 == last
            end
        end
        j = i
        i += 1
    end
    return true
end

# returns (pivot position, already_partitioned)
function _lc_partition_equals_right!(v, first, last, lt)
    begin_ = first
    pivot = v[first]
    while true
        first += 1
        lt(v[first], pivot) || break
    end
    if begin_ == first - 1
        while first < last
            last -= 1
            lt(v[last], pivot) && break
        end
    else
        while true
            last -= 1
            lt(v[last], pivot) && break
        end
    end
    already = first >= last
    while first < last
        _lc_swap!(v, first, last)
        while true
            first += 1
            lt(v[first], pivot) || break
        end
        while true
            last -= 1
            lt(v[last], pivot) && break
        end
    end
    pivot_pos = first - 1
    begin_ != pivot_pos && (v[begin_] = v[pivot_pos])
    v[pivot_pos] = pivot
    return pivot_pos, already
end

# returns the new `first` (one past the pivot)
function _lc_partition_equals_left!(v, first, last, lt)
    begin_ = first
    pivot = v[first]
    if lt(pivot, v[last - 1])
        while true
            first += 1
            lt(pivot, v[first]) && break
        end
    else
        while true
            first += 1
            (first < last && !lt(pivot, v[first])) || break
        end
    end
    if first < last
        while true
            last -= 1
            lt(pivot, v[last]) || break
        end
    end
    while first < last
        _lc_swap!(v, first, last)
        while true
            first += 1
            lt(pivot, v[first]) && break
        end
        while true
            last -= 1
            lt(pivot, v[last]) || break
        end
    end
    pivot_pos = first - 1
    begin_ != pivot_pos && (v[begin_] = v[pivot_pos])
    v[pivot_pos] = pivot
    return first
end

function _lc_introsort!(v, first, last, lt, depth, leftmost)
    limit = 24
    ninther = 128
    while true
        len = last - first
        if len <= 1
            return
        elseif len == 2
            lt(v[last - 1], v[first]) && _lc_swap!(v, first, last - 1)
            return
        elseif len == 3
            _lc_sort3!(v, first, first + 1, last - 1, lt); return
        elseif len == 4
            _lc_sort4!(v, first, first + 1, first + 2, last - 1, lt); return
        elseif len == 5
            _lc_sort5!(v, first, first + 1, first + 2, first + 3, last - 1, lt); return
        end
        if len < limit
            leftmost ? _lc_insertion_sort!(v, first, last, lt) : _lc_insertion_sort_unguarded!(v, first, last, lt)
            return
        end
        depth == 0 && error("libcxx_sort!: introsort depth exhausted; the heap-sort fallback is not implemented")
        depth -= 1
        half = len ÷ 2
        if len > ninther
            _lc_sort3!(v, first, first + half, last - 1, lt)
            _lc_sort3!(v, first + 1, first + half - 1, last - 2, lt)
            _lc_sort3!(v, first + 2, first + half + 1, last - 3, lt)
            _lc_sort3!(v, first + half - 1, first + half, first + half + 1, lt)
            _lc_swap!(v, first, first + half)
        else
            _lc_sort3!(v, first + half, first, last - 1, lt)
        end
        if !leftmost && !lt(v[first - 1], v[first])
            first = _lc_partition_equals_left!(v, first, last, lt)
            continue
        end
        i, already = _lc_partition_equals_right!(v, first, last, lt)
        if already
            fs = _lc_insertion_sort_incomplete!(v, first, i, lt)
            if _lc_insertion_sort_incomplete!(v, i + 1, last, lt)
                fs && return
                last = i
                continue
            elseif fs
                first = i + 1
                continue
            end
        end
        _lc_introsort!(v, first, i, lt, depth, leftmost)
        leftmost = false
        first = i + 1
    end
end

"""libc++ `std::sort(v.begin(), v.end(), lt)` for a non-arithmetic element type: the exact
(unstable) permutation libc++ 18 produces, see the block comment above."""
function libcxx_sort!(v::AbstractVector, lt)
    n = length(v)
    depth = n == 0 ? 0 : 2 * (8 * sizeof(Int) - 1 - leading_zeros(n))   # 2 * floor(log2 n)
    _lc_introsort!(v, firstindex(v), firstindex(v) + n, lt, depth, true)
    return v
end
