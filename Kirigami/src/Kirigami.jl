# Kirigami.jl -- kirigami deployment pipeline: reimplementation of Segall et al. 2025/2026
# plus this project's range-maximising constrained embedding method.
#
#   core/    planar mesh, cuts, holes, Tutte auxetic system, FK, collision, orientation,
#            generators, importer, optimiser, rank checks
#   method/  budget, deploy basis, mobility, contact, range optimisation, zero-plus,
#            convex/range embedding, design, expansive cone, periodic Jacobian
#   export/  layout, material, SVG, XML/3MF, STL, solid, zip
module Kirigami

using LinearAlgebra
using SparseArrays
using StaticArrays
using Random
using Statistics
using Printf
import JSON

const Vec2 = SVector{2,Float64}

include("core/mesh.jl")
include("core/cut.jl")
include("core/holes.jl")
include("core/tutte_auxetic.jl")
include("core/kinematics.jl")
include("core/collision.jl")
include("core/mt19937.jl")
include("core/orientation.jl")
include("core/generators.jl")
include("core/import_soup.jl")
include("core/kill_common.jl")
include("core/optimize.jl")
include("core/rank_checks.jl")

include("method/periodic_jacobian.jl")
include("method/deploy_basis.jl")
include("method/mobility.jl")
include("method/contact.jl")
include("method/range_opt.jl")
include("method/zero_plus.jl")
include("method/convex_embed.jl")
include("method/range_embed.jl")
include("method/design.jl")
include("method/expansive_cone.jl")
include("method/budget.jl")

# order: struct field types must exist before use (MaterialProfile in Layout,
# TriMesh in stl/threemf)
include("export/material.jl")
include("export/layout.jl")
include("export/xml.jl")
include("export/zip.jl")
include("export/svg.jl")
include("export/solid.jl")
include("export/stl.jl")
include("export/threemf.jl")

end # module
