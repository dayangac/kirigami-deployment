# Kirigami

*Exact deployment range and certified design for hinged kirigami on arbitrary planar graphs.*

## What this is

A hinged kirigami is a flat sheet cut into rigid tiles that stay connected at corner
hinges. When every tile rotates by the same angle the sheet opens uniformly, and Segall,
Ren and Sorkine-Hornung (TOG 2026) showed that for any planar graph the flat designs that
can open this way form a linear family, reached from an input drawing by a least-norm
projection. What that work leaves open is whether a design in the family actually opens
without its tiles running into each other.

This repository answers that question and follows the answer through:

- **An exact deployment range.** Along the uniform opening every contact test between
  two tiles is a single first harmonic in the opening angle, with coefficients that are
  quadratic in the flat design. The collision-free range is therefore a minimum over
  closed-form roots, once grazing contacts (tiles that touch at a point and separate) are
  told apart from real overlaps. A sound certificate that the range is at least a given
  angle follows. On 3,543 designs the closed form agrees with an independent brute-force
  referee on every one.
- **Why the published projection jams.** On random planar graphs the least-norm
  projection almost never opens: 1 of 200 designs at the sizes the original paper
  demonstrates, 0 of 400 at three times that size. The cause is that the projection
  creates reflex (inward-pointing) face corners that the input never had, and at a reflex
  corner the first-order separation condition becomes one no repair of the cut signs,
  orientation, boundary condition or non-uniform motion can satisfy. Every one of those
  alternatives is measured here.
- **A design method that does open.** Because the obstruction is first-order and exactly
  computable, it can be optimised against inside the same linear family: choose the design
  by its first-order deployment margin, then by its exact range, under barriers that keep
  every corner convex and every cut opening outward. That deploys 192 of the 200 small
  designs and 312 of the 400 large ones, with the same referee as the negative results.
- **Periodic patterns.** On a torus the deployment Jacobian is
  `cos(θ/2)·I + sin(θ/2)·K` with `K` affine in the design, which proves the conformality
  the original paper observed and gives the second closed state in closed form.

The code is a Julia package (`Kirigami/`) with a reimplementation of the published
pipeline at measured parity with the authors' C++, the exact contact calculus, the
certificate, the design method, fabrication export (SVG/STL/3MF) and a desktop explorer.
Every experiment behind the numbers above has a driver in `Kirigami/apps/`, a frozen input
population in `data/corpus/` and its outputs and write-up under `results/`
(`results/experiments/EXPERIMENTS.md` is the record). See `Kirigami/README.md` for the
package layout, `docs/NUMERICS.md` for what is reproducible to the last digit and what is
not, and `baseline/README.md` for the third-party code and how it was compared against.

## Quick start

```
export PATH=$HOME/.juliaup/bin:$PATH                       # native arm64 Julia on macOS
julia --project=Kirigami -e 'using Pkg; Pkg.instantiate()'
julia --project=Kirigami -e 'using Pkg; Pkg.test()'        # 187 test sets, 0 failures

APP="julia --project=Kirigami Kirigami/apps"
$APP/kiri_gen.jl delaunay 60 40 --out g.json --orient auto --seed 3   # a random graph
$APP/kiri_design.jl g.json --sigma both --out design.json --referee   # range-maximising design
$APP/kiri_analyze.jl g.json --out out/                                # exact range + certificate
$APP/kiri_export.jl design.json --profile felt_laser --theta-frac 0.5 --svg out.svg
```

`Kirigami/scripts/rerun_all.sh` reruns every experiment in dependency order;
`Kirigami/scripts/teaser_cat.jl` builds the cat-head example on the papers' first page.
