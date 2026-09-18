# export/svg.jl -- laser-cutter SVG for a Layout, in millimetres, with the three usual
# layers as SVG groups: "cut", "score", "engrave".
#
# Face / edge ids in the SVG are 0-based (the JSON file's numbering) so the files are
# byte-identical to the archived exports.

mutable struct SvgOptions
    margin::Float64       # mm of empty stock around the bounding box
    sigma_arrows::Bool    # engrave a rotation arrow per face (2025 Fig. 1)
    face_ids::Bool        # engrave the face index
    stroke::Float64       # mm; hairline for the cut layer
end
SvgOptions() = SvgOptions(5.0, true, true, 0.1)

# std::fixed << setprecision(4), trailing zeros stripped, "-0" -> "0"
function svg_num(v::Float64)
    s = @sprintf("%.4f", v)
    if occursin('.', s)
        s = rstrip(s, '0')
        s = rstrip(s, '.')
    end
    s == "-0" && (s = "0")
    return s
end

"""The SVG document for a layout (see the header comment)."""
function svg_string(L::Layout, opt::SvgOptions = SvgOptions())
    num = svg_num
    W = width(L) + 2 * opt.margin
    H = height(L) + 2 * opt.margin
    # SVG y grows downwards; the kirigami y grows upwards.
    px(p::Vec2) = p[1] - L.bb_min[1] + opt.margin
    py(p::Vec2) = L.bb_max[2] - p[2] + opt.margin
    function path_of(poly::Vector{Vec2})
        d = IOBuffer()
        for (i, p) in enumerate(poly)
            print(d, i > 1 ? " L " : "M ", num(px(p)), " ", num(py(p)))
        end
        print(d, " Z")
        return String(take!(d))
    end

    os = IOBuffer()
    print(os, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",
          "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"", num(W),
          "mm\" height=\"", num(H), "mm\" viewBox=\"0 0 ", num(W), " ", num(H), "\">\n")
    print(os, "<desc>kiri export: profile=", L.profile.name, " hinge=", to_string(L.profile.hinge),
          " theta=", num(L.theta), " rad, scale=", num(L.scale),
          " mm/unit, units=mm, faces=", n_faces(L), "</desc>\n")

    # ---- cut ----------------------------------------------------------------
    cut_style = "fill=\"none\" stroke=\"#000000\" stroke-width=\"" * num(opt.stroke) * "\""
    print(os, "<g id=\"cut\" ", cut_style, ">\n")
    for fp in L.pieces
        length(fp.outline) < 3 && continue
        print(os, "  <path id=\"cut_face_", fp.face - 1, "\" d=\"", path_of(fp.outline), "\"/>\n")
    end
    for fp in L.pieces, h in fp.pin_holes
        print(os, "  <circle id=\"cut_pin_", fp.face - 1, "_", num(px(h)), "\" cx=\"", num(px(h)),
              "\" cy=\"", num(py(h)), "\" r=\"", num(pin_hole_radius(L.profile)), "\"/>\n")
    end
    print(os, "</g>\n")

    # ---- score --------------------------------------------------------------
    # The neck bridge of every hinge: a line across the retained material so a
    # living hinge folds where it is meant to. Nothing is scored at a split edge.
    print(os, "<g id=\"score\" fill=\"none\" stroke=\"#0000ff\" stroke-width=\"", num(opt.stroke),
          "\">\n")
    if L.profile.hinge == LivingHingeNeck
        for hs in L.hinges
            a = hs.p
            b = hs.p + hs.neck * hs.dir_a
            print(os, "  <path id=\"score_hinge_", hs.edge - 1, "\" d=\"M ", num(px(a)), " ", num(py(a)),
                  " L ", num(px(b)), " ", num(py(b)), "\"/>\n")
            if norm(hs.dir_a - hs.dir_b) > 1e-9
                c = hs.p + hs.neck * hs.dir_b
                print(os, "  <path id=\"score_hinge_", hs.edge - 1, "_b\" d=\"M ", num(px(a)), " ",
                      num(py(a)), " L ", num(px(c)), " ", num(py(c)), "\"/>\n")
            end
        end
    end
    print(os, "</g>\n")

    # ---- engrave ------------------------------------------------------------
    print(os, "<g id=\"engrave\" fill=\"none\" stroke=\"#ff0000\" stroke-width=\"", num(opt.stroke),
          "\">\n")
    for fp in L.pieces
        length(fp.raw) < 3 && continue
        ctr = Vec2(0, 0)
        for q in fp.raw
            ctr += q
        end
        ctr /= length(fp.raw)
        rmin = 1e30
        for q in fp.raw
            rmin = min(rmin, norm(q - ctr))
        end
        r = min(0.35 * rmin, 4.0)
        sg = (fp.face <= length(L.mesh.sigma)) ? L.mesh.sigma[fp.face] : -1
        if opt.sigma_arrows && r > 0.4
            # A 270-degree arc with a head; drawn clockwise for sigma = +1.
            sweep = (sg > 0 ? -1.0 : 1.0) * 1.5 * pi
            a0 = 0.4
            a1 = a0 + sweep
            p0 = ctr + r * Vec2(libm_cos(a0), libm_sin(a0))
            p1 = ctr + r * Vec2(libm_cos(a1), libm_sin(a1))
            large = 1
            sflag = (sg > 0) ? 0 : 1  # SVG sweep flag: y is flipped, so invert
            print(os, "  <path id=\"sigma_", fp.face - 1, "\" d=\"M ", num(px(p0)), " ", num(py(p0)),
                  " A ", num(r), " ", num(r), " 0 ", large, " ", sflag, " ", num(px(p1)),
                  " ", num(py(p1)), "\"/>\n")
            # arrow head: two short barbs at p1, tangent +/- 30 degrees
            tang = a1 + (sg > 0 ? -1.0 : 1.0) * pi / 2
            hl = min(0.5 * r, 1.5)
            for k in (-1, 1)
                ang = tang + pi + k * 0.5
                q = p1 + hl * Vec2(libm_cos(ang), libm_sin(ang))
                print(os, "  <path id=\"sigma_head_", fp.face - 1, "_", (k > 0 ? "p" : "m"), "\" d=\"M ",
                      num(px(p1)), " ", num(py(p1)), " L ", num(px(q)), " ", num(py(q)),
                      "\"/>\n")
            end
        end
        if opt.face_ids
            fs = max(1.0, min(0.5 * rmin, 3.0))
            print(os, "  <text id=\"faceid_", fp.face - 1, "\" x=\"", num(px(ctr)), "\" y=\"",
                  num(py(ctr) + 0.35 * fs), "\" font-size=\"", num(fs),
                  "\" text-anchor=\"middle\" fill=\"#ff0000\" stroke=\"none\">", fp.face - 1,
                  (sg > 0 ? "+" : "-"), "</text>\n")
        end
    end
    print(os, "</g>\n</svg>\n")
    return String(take!(os))
end

"""Writes `svg_string(L, opt)` to `path`."""
function write_svg(L::Layout, path::AbstractString, opt::SvgOptions = SvgOptions())
    open(path, "w") do io
        write(io, svg_string(L, opt))
    end
    return nothing
end
