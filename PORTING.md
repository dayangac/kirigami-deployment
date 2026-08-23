# Porting conventions (Julia)

Reference implementation being ported: `~/Documents/kirigami-experiments/code` (C++/Eigen).
Read the C++ source AND its tests for each module before writing Julia. The Julia code must
reproduce the C++ numbers (tolerance-based, never exact float equality) on the frozen corpora
in `data/corpus/`.

## Package
- Package root: `Kirigami/` (Project.toml). Run tests with
  `julia --project=Kirigami -e 'using Pkg; Pkg.test()'` or include a single file via
  `julia --project=Kirigami Kirigami/test/<file>.jl` after `include("helpers.jl")`.
- One Julia file per C++ translation unit, same basename, under `src/core|method|export`.
  All files are `include`d into the single `Kirigami` module (see `src/Kirigami.jl`); do not
  create submodules. Export nothing; tests call `Kirigami.f(...)`.
- Keep C++ function names (snake_case) so every derivation/paper reference still resolves.
  Struct fields keep their C++ names. Enum -> `@enum`. `std::optional` -> `Union{T,Nothing}`.
- Indices: Julia is 1-based. Internally use 1-based everywhere. JSON I/O keeps the file
  format 0-based (as the C++ wrote it) and converts at the boundary. Say so in a comment at
  every conversion site.
- `Vec2 = SVector{2,Float64}` (defined in Kirigami.jl). Matrices: `Matrix{Float64}` /
  `SparseMatrixCSC`. Eigen dense decompositions -> LinearAlgebra (`qr(A, ColumnNorm())`,
  `svd`, `nullspace`, `rank(A; atol)`); `SparseQR`/`SparseLU` -> SparseArrays/SuiteSparse
  (`qr(sparse)`, `lu(sparse)`). Preserve the tolerances the C++ uses.
- RNG: `std::mt19937` is ported bit-exactly in `core/mt19937.jl` (`MT19937` type, `next_u32`,
  `uniform_real(rng, a, b)` = libc++ `uniform_real_distribution` via
  `generate_canonical<double,53>` with two 32-bit draws). Generators take an `MT19937`.
- Errors: C++ `throw std::runtime_error` -> `error(...)`/`throw(ArgumentError(...))`, and the
  test asserts `@test_throws`.
- Style: short docstring per public function stating what it computes and the paper equation
  it implements (copy the C++ header comment). Match the C++ comment density; no filler.

## Tests
- `Kirigami/test/<name>.jl` mirrors `code/tests/<name>.cpp` case by case: each doctest
  `TEST_CASE` -> `@testset`, each `CHECK` -> `@test`, `CHECK(x == doctest::Approx(y).epsilon(e))`
  -> `@test isapprox(x, y; rtol=e)` (or `atol` when the C++ used absolute margins).
- Files start with `include("helpers.jl")`. Never depend on ordering between files.
- Reference numbers copied from the C++ tests are the acceptance criterion. If a Julia value
  disagrees, the port is wrong until proven otherwise (investigate; do not loosen tolerances).

## Frozen corpora (`data/corpus/`)
- `k1a_200.json`: the 200-graph population of `kill_common.hpp::make_graph` (ids 0..199) with
  sigma; `native200.json`: the native200 corpus; `mt19937_vectors.json`: RNG test vectors.
  Produced by the C++ before its removal. Julia apps load these; they never regenerate them.

## Reporting back
End your report with: files written, test counts (`Pass/Fail/Error`), every C++ behaviour you
could NOT reproduce and why, and any C++ bug you found (do not silently fix; replicate and flag).

## Bit-exactness note (from the mt19937 port)
The C++ was built with clang -O2 on arm64, default `-ffp-contract=on`: `a*b + c` patterns are
FMA-fused. Where a value must be reproduced bit-exactly (RNG streams, corpus generation), use
`fma(a, b, c)` with the same operand placement. For tolerance-based comparisons this is irrelevant.

## Shared helpers (avoid duplicate definitions — precompile refuses method overwriting)
`core/mesh.jl` owns: `_det2(u,v)` (2-D cross product) and the `libm_cos/sin/tan/pow/atan2`
shims (system libm on macOS, Base elsewhere). `export/layout.jl` owns `polygon_area`.
`cut.jl` owns `DSU`/`find!`/`unite!`. Before defining a small helper, `grep -rn` for it.

## Julia binary
Use the native arm64 Julia: `export PATH=$HOME/.juliaup/bin:$PATH` (juliaup, 1.12.7 aarch64).
`/usr/local/bin/julia` is the x86_64 Homebrew build under Rosetta: slower, and its libm
differs from the arm64 libm the C++ reference used (1-ulp trig/hypot differences).

## Bit-faithful trig and reductions (from the derivation-test port)
- Where the C++ evaluates `std::cos(x)` and `std::sin(x)` of the same argument in one function,
  clang on Apple fuses them into `__sincos_stret`, whose sine differs from standalone `sin` by
  1 ulp on ~4% of arguments. Use `Kirigami.libm_sincos(x)` there, not `libm_sin`/`libm_cos`.
- Eigen `dot()`/`squaredNorm()` are UNFUSED reductions; StaticArrays `dot` is a muladd that
  becomes an fma on aarch64. In tie-sensitive predicates write the C++ contraction explicitly:
  `a0*b0 + a1*b1` unfused for Eigen dots; `fma(ux, vy, -(uy*vx))` for inline `ux*vy - uy*vx`.

## Loop semantics
`for a in A, b in B ... break` is ONE loop in Julia: `break` leaves both. A C++ nested loop that
breaks only the inner one must be written as explicitly nested `for` blocks. (`continue` is safe.)
Swept 2026-09-20: the only instance was in check_l1.jl, fixed.
