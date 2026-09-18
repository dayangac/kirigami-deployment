# Numerical conventions

The frozen corpora in `data/corpus/` (input populations, RNG streams, test fixtures) and
the committed results under `results/` are the acceptance criterion of the code: a change
that moves a bit-exact quantity or loosens a tolerance is wrong until proven otherwise.
The rules below are what keep the numbers reproducible; every rule names the place in
the code that applies it.

## Platform

- Use the native arm64 Julia from juliaup: `export PATH=$HOME/.juliaup/bin:$PATH`. The
  package declares `julia = "1.10"`; the reference runs were made on 1.12 aarch64.
- The corpora were produced with Apple's arm64 libm. `/usr/local/bin/julia` (the x86_64
  Homebrew build under Rosetta) reaches Apple's x86_64 libm instead, which differs by
  1 ulp in trig/hypot on a fraction of arguments; that shows up as equal-objective tie
  flips in the orientation relaxation on symmetric tilings (`core/generators.jl`). Random
  graphs are unaffected. The bit-exact tests key on `Kirigami._USE_SYSTEM_LIBM`
  (`core/mesh.jl`) and relax to tolerances elsewhere.

## Fused multiply-add placement

Where a value must be reproduced bit for bit (RNG streams, corpus generation, the
closed-form contact calculus, STL normals), the placement of `fma` is part of the
definition of the result:

- Plain scalar expressions `a*b + c` and `a*b - c*d` are FUSED: write `fma(a, b, c)` and
  `fma(a, b, -(c*d))` with the same operand order. `core/mesh.jl` `_det2`,
  `method/deploy_basis.jl` `det2`, `export/solid.jl` `fms`/`cross_fma`, the circumcircle
  of `delaunay_triangles`, `uniform_real` and `normal` in `core/mt19937.jl`.
- Dot products, squared norms, norms and lerps `A + (B-A)*t` are UNFUSED reductions:
  write `a[1]*b[1] + a[2]*b[2]` explicitly (`_dot2`/`_sqnorm2` in `core/generators.jl`,
  `_dot`/`_sq` in `core/collision.jl`, `_dotu`/`_sqnormu` in `method/deploy_basis.jl`).
  StaticArrays' `dot` is a `muladd` that becomes an fma on aarch64, so it is not used in
  tie-sensitive predicates.
- The 2x2 rotation applied to a vector is `y_i = fma(R[i,2], x[2], R[i,1]*x[1])`
  (`core/kinematics.jl` `_mv`); a plain `R * x` differs by an ulp.
- A chained sum `p + q*cos + r*sin` is contracted left to right:
  `fma(r, s, fma(q, c, p))` (`harmonic_eval`).

For tolerance-based comparisons the placement is irrelevant; do not spread `fma` beyond
the bit-exact paths.

## Trig and reductions

- Trig, `pow`, `atan2`, `hypot`, `acos` go through the `libm_*` shims of `core/mesh.jl`
  and `method/deploy_basis.jl`, which call the system libm on macOS and Base elsewhere.
- Where BOTH `sin(x)` and `cos(x)` of one argument are needed in one function, use
  `libm_sincos(x)`, which calls Apple's combined `__sincos_stret`; its sine differs from
  the standalone `sin` by 1 ulp on ~4% of arguments. The rotation of `core/kinematics.jl`,
  `harmonic_eval`, `basis_eval` and the range objective all use it.
- `Float64(pi)`, never the `Irrational` `pi`, in comparisons: `pi <= x` is an exact
  comparison, so `pi <= Float64(pi)` is false, and a root at theta = pi is dropped.

## Random numbers

- `core/mt19937.jl` is a self-contained MT19937 with libc++-compatible (LLVM 18,
  `_LIBCPP_VERSION 180100`) distributions: `generate_canonical<double,53>` with two 32-bit
  draws and no clamp of 1.0, `uniform_real_distribution` evaluated as one `fma`,
  `uniform_int_distribution` by rejection on the working unsigned type, `std::shuffle`,
  `normal_distribution` (Marsaglia polar, cached second variate). Seeds reproduce the
  frozen design populations; `data/corpus/mt19937_vectors.json` holds the reference
  vectors and `test/test_mt19937.jl` locks them.
- Draw ORDER is part of the definition of a population (column-major vs row-major fills,
  one distribution object vs a fresh one per call, interleaved draws): the drivers and
  tests spell it out at each site. Do not reorder draws.
- `libcxx_sort!` (same file) is the libc++ 18 introsort, used where an unstable sort's tie
  order decides a result (`detect_lattice` in `method/periodic_jacobian.jl`).

## Linear algebra

- The Tutte auxetic system is solved DENSE with `LinearAlgebra.svd`; the rank rule is
  `sv > max(rel_tol * smax, 1e-14) * max(1, rows, cols)` (`core/tutte_auxetic.jl`). Large
  systems use `qr(sparse; tol = 1e-9)` for the rank only. Dense ranks elsewhere use
  column-pivoted Householder QR with `|R_ii| > rel_tol * |R_11|`.
- A null-space basis `Phi` is defined only up to an orthogonal change of basis. Fixtures
  that depend on Gaussian draws in null-space coordinates therefore store `Phi` and the
  tests replay on the STORED basis after checking that a fresh solve spans the same
  subspace (`Phi Phi^T` to 1e-9) and gives the same `X0` (to 1e-9).
- BLAS threads: every app sets `BLAS.set_num_threads(1)` (`apps/common_app.jl`). OpenBLAS's
  threaded gemv splits reductions across threads, so optimiser trajectories (L-BFGS
  repairs, Frank-Wolfe) depend on the thread count; sharded runs with the default thread
  count also oversubscribe the machine (observed 40-100x slowdowns). Library users keep
  the default.

## Optimiser paths

`characterize`, the certificate, the generators, the RNG and the linear algebra reproduce
the committed numbers to the tolerances the tests state (bit-exact for RNG streams,
generated meshes and the closed-form contact calculus at a given point). Optimiser
OUTPUTS (`convex_embed`, `zero_plus_repair`, `range_embed`, `maximize_margin_range`) are
reproduced to path level: same objective, gradient and start points to 1e-14, converged
solves agree to solver tolerance, and the fixed-budget `range_embed` (6 x 120 L-BFGS
iterations, never converged) diverges by rounding amplification, so its end point differs
at the 1e-2 level. Locks on such outputs are kept verbatim as `@test_broken`
("optimiser-path-dependent lock"), documented in `test/test_design.jl` and
`data/corpus/method_fixtures/README_design.md`; they turn into an error the day the path
becomes bit-exact. Do not loosen them and do not delete them.

## Loop-break semantics

`for a in A, b in B ... break end` is ONE loop in Julia: `break` leaves both. A nested loop
that must break only the inner one is written as explicitly nested `for` blocks
(`continue` is safe either way). Swept 2026-09-20; the only instance was in
`derivations/scratch/check_l1.jl`.

## Frozen corpora

- Every input population is read from `data/corpus/` by default (provenance per file in
  `data/corpus/README.md`). Drivers accept `--regenerate`, which rebuilds the population
  through the generators and the bit-exact RNG (`apps/kill_common.jl`: `population`,
  `reference_cases`, `deployable_population`); on the reference platform `make_graph`
  regenerates all 200 graphs of `k1a_200.json` vertex for vertex. The frozen file remains
  the definition of the population.
- A run never rewrites the corpus. Apps write to the canonical `results/` tree
  (`results/kill/<x>/`, `results/final/e1/`, `results/{regime,scaling,yield}/`);
  compare a rerun against the committed copy with `git diff`.
- Deliberate quirks the corpora depend on are marked `deliberate` in the code (welder cell
  overwrite, `generate_canonical` reaching 1.0, DSU without second-argument path
  compression, noise rows in `normalise_rows!`); replicate and flag, do not fix silently.

## Shared helpers

Precompilation refuses method overwriting, so small helpers live in exactly one place:
`core/mesh.jl` owns `_det2` and the `libm_*` shims; `export/layout.jl` owns
`polygon_area`; `core/cut.jl` owns `DSU`/`find!`/`unite!`. `grep -rn` before defining
another.

## Tests

- `julia --project=Kirigami -e 'using Pkg; Pkg.test()'`, or a subset by basename
  substring: `Pkg.test(test_args=["design", "export"])`. Files start with
  `include("helpers.jl")` and never depend on ordering between files.
- Tolerances: `isapprox(x, y; rtol = e)`, or `atol` where the quantity is near zero. Never
  exact float equality except for the bit-exact locks, where `==` is the point.
- Reference numbers in the fixtures are the acceptance criterion. If a value disagrees,
  the code is wrong until proven otherwise: investigate, do not loosen the tolerance.
