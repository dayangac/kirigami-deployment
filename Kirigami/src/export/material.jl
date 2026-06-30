# export/material.jl -- the parametric hinge / material profile driving every writer
# (MISSION Sec. 7: "3MF/STL writer with a parametric hinge (face extrusion, hinge
# pad radius, neck width, gap) driven by a material profile struct").
#
# All lengths are millimetres. The profile is the ONLY place fabrication
# constants live; the SVG, STL and 3MF writers read it and never hard-code a
# number of their own.
#
# Port of code/src/export/material.{hpp,cpp}.

@enum HingeType begin
    LivingHingeNeck  # faces stay connected through a thin neck of material
    PinPad           # faces are separate bodies joined by a pin through a pad
end

to_string(t::HingeType) = t == PinPad ? "pin-pad" : "living-hinge-neck"

"""Parses a hinge type name; throws on unknown."""
function hinge_type_from_string(s::AbstractString)
    (s == "living-hinge-neck" || s == "neck") && return LivingHingeNeck
    (s == "pin-pad" || s == "pin" || s == "pad") && return PinPad
    error("unknown hinge type: $s")
end

mutable struct MaterialProfile
    name::String

    thickness::Float64    # sheet thickness / face extrusion height
    kerf::Float64         # laser beam width; each cut side loses kerf/2
    hinge::HingeType
    neck_width::Float64   # length of retained material at a hinge vertex
    pad_radius::Float64   # pin-pad outer radius
    gap::Float64          # total clearance between two neighbouring faces
    min_feature::Float64  # smallest producible feature; used by check()
    clearance::Float64    # pin / hole running clearance (3D only)
end
MaterialProfile() = MaterialProfile("custom", 2.0, 0.2, LivingHingeNeck, 1.5, 2.5, 0.3, 0.5, 0.15)

# Half-width removed from each face along a cut edge, so that the two faces
# of a cut edge end up `gap` apart after the kerf has burned away.
face_inset(p::MaterialProfile) = 0.5 * p.gap + 0.5 * p.kerf

# Radius of the hole through a pin pad, and of the pin that fits it.
pin_radius(p::MaterialProfile) = 0.4 * p.pad_radius
pin_hole_radius(p::MaterialProfile) = pin_radius(p) + p.clearance

"""Returns the list of violated design rules (empty == profile is producible)."""
function check(p::MaterialProfile)
    bad = String[]
    req(ok, msg) = ok || push!(bad, msg)
    req(p.thickness > 0, "thickness must be > 0")
    req(p.kerf >= 0, "kerf must be >= 0")
    req(p.gap >= 0, "gap must be >= 0")
    req(p.min_feature > 0, "min_feature must be > 0")
    req(p.neck_width >= p.min_feature, "neck_width < min_feature")
    if p.hinge == PinPad
        req(p.pad_radius >= p.min_feature, "pad_radius < min_feature")
        req(p.pad_radius - pin_hole_radius(p) >= 0.5 * p.min_feature,
            "pad wall (pad_radius - pin_hole_radius) below half the min feature")
        req(p.clearance > 0, "pin-pad hinge needs clearance > 0")
    end
    p.gap > 0 && req(p.gap >= 0.5 * p.min_feature, "gap positive but below half the min feature")
    return bad
end

# Presets required by specs/builder_export.md item 1.
function profile_felt_laser()   # 2 mm felt, neck 1.5 mm, gap 0.3 mm
    p = MaterialProfile()
    p.name = "felt_laser"
    p.thickness = 2.0
    p.kerf = 0.2
    p.hinge = LivingHingeNeck
    p.neck_width = 1.5
    p.pad_radius = 2.5
    p.gap = 0.3
    p.min_feature = 0.6  # 0.3 mm gap is exactly half of it, the design-rule floor
    p.clearance = 0.2
    return p
end

function profile_paper_laser()  # 0.3 mm card, neck 1.0 mm
    p = MaterialProfile()
    p.name = "paper_laser"
    p.thickness = 0.3
    p.kerf = 0.1
    p.hinge = LivingHingeNeck
    p.neck_width = 1.0
    p.pad_radius = 1.5
    p.gap = 0.2
    p.min_feature = 0.4
    p.clearance = 0.1
    return p
end

function profile_pla_print()    # 2 mm PLA, pin-pad, pad r 2.5 mm
    p = MaterialProfile()
    p.name = "pla_print"
    p.thickness = 2.0
    p.kerf = 0.0  # additive: no material is burned away
    p.hinge = PinPad
    p.neck_width = 1.5
    p.pad_radius = 2.5
    p.gap = 0.4
    p.min_feature = 0.8
    p.clearance = 0.15
    return p
end

profile_names() = ["felt_laser", "paper_laser", "pla_print"]

"""Preset by name; throws on unknown."""
function profile_by_name(name::AbstractString)
    name == "felt_laser" && return profile_felt_laser()
    name == "paper_laser" && return profile_paper_laser()
    name == "pla_print" && return profile_pla_print()
    error("unknown material profile '$name'; known: " * join(profile_names(), " "))
end
