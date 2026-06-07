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
