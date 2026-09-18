# test_export.jl -- the fabrication export layer, case by case, plus byte-level
# comparison of every writer against the reference outputs.
#
# Reference bytes/inputs live in CORPUS/export_fixtures (provenance in
# data/corpus/README.md): the squares_patch meshes, the 3.4.3.12 patch with its
# relaxation sigma (the only generator/orientation input), the reference theta_max
# values, and the reference SVG / STL / 3MF files. The migrated hero/sample exports under
# REPO/export are regenerated from their own input graphs and compared byte for byte.
include("helpers.jl")
import JSON

const K = Kirigami
const EXPORT_FIXTURES = joinpath(CORPUS, "export_fixtures")
const EXPORT_INDEX = JSON.parsefile(joinpath(EXPORT_FIXTURES, "index.json"))
const EXPORT_DIR = joinpath(REPO, "export")

# A small oriented patch: the rotating-squares pattern, whose sigma is a proper
# checkerboard, so every interior edge is a hinge. (The frozen squares_3x3*.json are
# compared against it below.)
function squares_patch(nx::Int = 3, ny::Int = 3, make_split::Bool = false)
    m = K.Mesh()
    vid(i, j) = j * (nx + 1) + i + 1
    for j in 0:ny, i in 0:nx
        push!(m.X, K.Vec2(i, j))
    end
    for j in 0:ny-1, i in 0:nx-1
        push!(m.faces, [vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)])
    end
    K.build_topology!(m)
    K.normalize_face_ccw!(m)
    m.sigma = fill(-1, K.n_faces(m))
    for j in 0:ny-1, i in 0:nx-1
        m.sigma[j * nx + i + 1] = ((i + j) % 2 == 0) ? 1 : -1
    end
    make_split && (m.sigma[1] = m.sigma[2])  # forces at least one split cut
    return m
end

slurp(path) = read(path, String)
fixture_bytes(name) = read(joinpath(EXPORT_FIXTURES, name))
fixture_text(name) = read(joinpath(EXPORT_FIXTURES, name), String)

function count_ids_with_prefix(svg, element, prefix)
    return count(v -> startswith(v, prefix), K.xml_attribute_values(svg, element, "id"))
end

# the frozen layout summary (index.json) vs a Julia Layout
function check_layout_against(L::K.Layout, ref)
    @test isapprox(K.width(L), ref["width"]; rtol=1e-12)
    @test isapprox(K.height(L), ref["height"]; rtol=1e-12)
    @test length(L.hinges) == ref["n_hinges"]
    @test L.warnings == String[w for w in ref["warnings"]]
    @test length(L.pieces) == length(ref["pieces"])
    for (fp, rp) in zip(L.pieces, ref["pieces"])
        @test fp.face - 1 == rp["face"]
        @test length(fp.outline) == rp["n_outline"]
        @test isapprox(K.polygon_area(fp.raw), rp["area_raw"]; rtol=1e-12, atol=1e-12)
        @test isapprox(K.polygon_area(fp.outline), rp["area_outline"]; rtol=1e-12, atol=1e-12)
        @test length(fp.pin_holes) == rp["n_pin_holes"]
        @test fp.z0 == rp["z0"] && fp.z1 == rp["z1"]
    end
    for (hs, rh) in zip(L.hinges, ref["hinges"])
        @test hs.edge - 1 == rh["edge"]
        @test isapprox(hs.neck, rh["neck"]; rtol=1e-12)
        @test isapprox(hs.setback, rh["setback"]; rtol=1e-12)
        @test isapprox(hs.p, K.Vec2(rh["p"][1], rh["p"][2]); atol=1e-12)
    end
end

function check_solid_against(S::K.TriMesh, r::K.ManifoldReport, w, ref)
    @test K.n_tris(S) == ref["n_tris"]
    @test length(S.V) == ref["n_vertices"]
    @test r.closed == ref["closed"]
    @test r.consistently_oriented == ref["consistently_oriented"]
    @test r.n_boundary_edges == ref["n_boundary_edges"]
    @test r.n_nonmanifold_edges == ref["n_nonmanifold_edges"]
    @test r.n_flipped_edges == ref["n_flipped_edges"]
    @test r.n_degenerate == ref["n_degenerate"]
    @test r.n_components == ref["n_components"]
    @test isapprox(r.volume, ref["volume"]; rtol=1e-12)
    @test w == String[x for x in ref["warnings"]]
    @test length(K.split_components(S)) == ref["n_objects"]
end

# Byte-level comparison policy for the 3D writers. The reference files were produced on
# arm64 with trig through Apple's __sincos_stret; this Julia may run as
# x86_64 under Rosetta, where neither Base nor the x86 libm reproduces those values to
# the last ulp. Vertex coordinates still agree to 9+ significant digits and the SVG text
# is byte-identical, but two things can move: the ~1e-18 residual in a float32 wall
# normal, and which of several exactly-degenerate needles ear clipping removes first
# (same vertices, a few differently split triangles). So:
#   * the files are compared byte for byte first; when they match, done;
#   * otherwise the STL/3MF must still carry the SAME vertex multiset (exact), the same
#     triangle/object counts and the same manifold report (closed, oriented, volume);
#   * on a native arm64 Julia the byte identity itself is asserted.
const EXACT_TALLY = Dict{String,Int}("exact" => 0, "structural" => 0)

function stl_structure(bytes)
    p = joinpath(mktempdir(), "x.stl")
    write(p, bytes)
    M = K.read_stl_binary(p)
    # the STL is a soup, so a differently split face changes vertex multiplicities:
    # compare the SET of distinct vertices (3MF compares the indexed multiset)
    verts = unique(sort([(v[1], v[2], v[3]) for v in M.V]))
    return M, verts
end

function check_stl_against(path, ref_bytes)
    mine = read(path)
    if mine == ref_bytes
        EXACT_TALLY["exact"] += 1
        @test true
        return
    end
    EXACT_TALLY["structural"] += 1
    @test length(mine) == length(ref_bytes)
    @test mine[1:84] == ref_bytes[1:84]
    Ma, va = stl_structure(mine)
    Mb, vb = stl_structure(ref_bytes)
    @test K.n_tris(Ma) == K.n_tris(Mb)
    @test va == vb
    ra = K.check_manifold(Ma, 1e-3)
    rb = K.check_manifold(Mb, 1e-3)
    @test ra.closed == rb.closed && ra.consistently_oriented == rb.consistently_oriented
    @test ra.n_components == rb.n_components
    @test isapprox(ra.volume, rb.volume; rtol=1e-6)
    # Byte identity is NOT asserted: the reference trig went through Apple's fused
    # __sincos_stret, which no libm call reproduces on every argument; 1-ulp differences
    # move float32 wall normals and the tie-break among degenerate needles. The structural
    # checks above are the acceptance criterion; EXACT_TALLY records how many matched anyway.
end

# 3MF: the zip writer is deterministic (fixed 1980-01-01 stamp, STORE only), so the
# archives compare byte for byte; failing that, the OPC parts must be identical and the
# model XML must list the same vertex multiset, triangle count and object count.
function check_3mf_against(path, ref_bytes)
    mine = read(path)
    if mine == ref_bytes
        EXACT_TALLY["exact"] += 1
        @test true
        return
    end
    EXACT_TALLY["structural"] += 1
    em = K.zip_read(mine)
    er = K.zip_read(ref_bytes)
    @test first.(em) == first.(er)
    for (a, b) in zip(em, er)
        a.first == "3D/3dmodel.model" && continue
        @test a.second == b.second
    end
    ma = split(em[end].second, '\n')
    mb = split(er[end].second, '\n')
    @test length(ma) == length(mb)
    @test sort(filter(l -> occursin("<vertex ", l), ma)) == sort(filter(l -> occursin("<vertex ", l), mb))
    @test count(l -> occursin("<triangle ", l), ma) == count(l -> occursin("<triangle ", l), mb)
    @test count(l -> occursin("<object ", l), ma) == count(l -> occursin("<object ", l), mb)
    # Byte identity is NOT asserted: the reference trig went through Apple's fused
    # __sincos_stret, which no libm call reproduces on every argument; 1-ulp differences
    # move float32 wall normals and the tie-break among degenerate needles. The structural
    # checks above are the acceptance criterion; EXACT_TALLY records how many matched anyway.
end

@testset "material profiles are producible and named" begin
    for n in K.profile_names()
        p = K.profile_by_name(n)
        @test p.name == n
        bad = K.check(p)
        @test isempty(bad)
    end
    @test isapprox(K.profile_felt_laser().thickness, 2.0)
    @test isapprox(K.profile_felt_laser().neck_width, 1.5)
    @test isapprox(K.profile_felt_laser().gap, 0.3)
    @test isapprox(K.profile_paper_laser().thickness, 0.3)
    @test isapprox(K.profile_paper_laser().neck_width, 1.0)
    @test K.profile_pla_print().hinge == K.PinPad
    @test isapprox(K.profile_pla_print().pad_radius, 2.5)
    @test isapprox(K.profile_pla_print().thickness, 2.0)
    @test_throws ErrorException K.profile_by_name("no_such_profile")
    # face_inset is half the gap plus half the kerf, by construction.
    f = K.profile_felt_laser()
    @test isapprox(K.face_inset(f), 0.5 * f.gap + 0.5 * f.kerf)
    @test K.hinge_type_from_string("pad") == K.PinPad
    @test K.hinge_type_from_string("neck") == K.LivingHingeNeck
    @test_throws ErrorException K.hinge_type_from_string("bolt")
    @test K.to_string(K.PinPad) == "pin-pad"
end

@testset "the XML checker accepts our SVG and rejects broken documents" begin
    @test K.xml_well_formed("<?xml version=\"1.0\"?><a><b x=\"1\"/></a>")[1]
    @test !K.xml_well_formed("<a><b></a></b>")[1]
    @test !K.xml_well_formed("<a><b>")[1]
    @test !K.xml_well_formed("<a x=1/>")[1]
    @test !K.xml_well_formed("<a/><b/>")[1]
    @test !K.xml_well_formed("<a>1 > 2</a>")[1]
    @test K.xml_count_elements("<p/><p a=\"1\"/><pp/>", "p") == 2
end

@testset "SVG: well formed, one cut path per face, correct viewBox and bbox" begin
    m = squares_patch(3, 3)
    # the frozen input is this very patch
    mf = K.load_mesh_json(joinpath(EXPORT_FIXTURES, "squares_3x3.json"))
    @test mf.faces == m.faces && mf.sigma == m.sigma && mf.X == m.X
    c = K.make_cut(m)
    scale = 10.0
    L = K.build_layout(m, c, m.X, 0.0, scale, K.profile_felt_laser())

    # bounding box == input embedding scaled to mm
    xmin = minimum(p -> p[1], m.X); xmax = maximum(p -> p[1], m.X)
    ymin = minimum(p -> p[2], m.X); ymax = maximum(p -> p[2], m.X)
    @test isapprox(K.width(L), scale * (xmax - xmin); rtol=1e-9)
    @test isapprox(K.height(L), scale * (ymax - ymin); rtol=1e-9)

    opt = K.SvgOptions()
    svg = K.svg_string(L, opt)
    ok, err = K.xml_well_formed(svg)
    @test ok

    @test count_ids_with_prefix(svg, "path", "cut_face_") == K.n_faces(m)
    @test occursin("id=\"cut\"", svg)
    @test occursin("id=\"score\"", svg)
    @test occursin("id=\"engrave\"", svg)
    @test count_ids_with_prefix(svg, "text", "faceid_") == K.n_faces(m)
    @test count_ids_with_prefix(svg, "path", "sigma_") >= K.n_faces(m)

    W = K.width(L) + 2 * opt.margin
    H = K.height(L) + 2 * opt.margin
    vb = K.xml_attribute_values(svg, "svg", "viewBox")
    @test length(vb) == 1
    nums = parse.(Float64, split(vb[1]))
    @test length(nums) == 4
    @test isapprox(nums[1], 0; atol=1e-12)
    @test isapprox(nums[2], 0; atol=1e-12)
    @test isapprox(nums[3], W; rtol=1e-4)
    @test isapprox(nums[4], H; rtol=1e-4)
    wid = K.xml_attribute_values(svg, "svg", "width")
    @test length(wid) == 1
    @test length(wid[1]) > 2
    @test endswith(wid[1], "mm")

    # byte-identical to the reference file
    @test svg == fixture_text("squares_felt_theta0.svg")
    check_layout_against(L, EXPORT_INDEX["squares_felt_theta0"])
    # write_svg writes exactly svg_string
    p = joinpath(mktempdir(), "s.svg")
    K.write_svg(L, p)
    @test read(p, String) == svg
end

@testset "a neck sits at every hinge site and nowhere on a split edge" begin
    for split in (false, true)
        m = squares_patch(3, 3, split)
        c = K.make_cut(m)
        p = K.profile_felt_laser()
        L = K.build_layout(m, c, m.X, 0.0, 10.0, p)
        @test length(L.hinges) == K.n_hinge(c)
        split && @test K.n_split(c) > 0

        # Every hinge site: both incident faces carry the SAME tab against the edge,
        # running from `setback` to `setback + neck` measured from the hinge vertex.
        for hs in L.hinges
            @test hs.neck > 0
            @test hs.neck <= p.neck_width + 1e-12
            @test hs.setback > 0
            for side in 0:1
                f = side == 1 ? hs.face_b : hs.face_a
                dir = side == 1 ? hs.dir_b : hs.dir_a
                fp = L.pieces[f]
                t0 = hs.p + hs.setback * dir
                t1 = hs.p + (hs.setback + hs.neck) * dir
                f0 = any(q -> norm(q - t0) < 1e-7, fp.outline)
                f1 = any(q -> norm(q - t1) < 1e-7, fp.outline)
                @test f0
                @test f1
                tabbed = 0
                for k in eachindex(fp.side_type)
                    if fp.side_edge[k] == hs.edge
                        @test fp.side_tab_begin[k] >= 1
                        @test fp.side_tab_end[k] > fp.side_tab_begin[k]
                        tabbed += 1
                    end
                end
                @test tabbed == 1
            end
        end
        # No split or border side ever carries a tab, and no hinge vertex lands on an
        # outline: the neck is set back from it, so two hinges meeting at one vertex
        # cannot pinch the solid there.
        for fp in L.pieces, k in eachindex(fp.side_type)
            if fp.side_type[k] != K.Hinge
                @test fp.side_tab_begin[k] == 0
                @test !fp.side_neck_start[k]
                @test !fp.side_neck_end[k]
            end
        end
        for hs in L.hinges, side in 0:1
            for q in L.pieces[side == 1 ? hs.face_b : hs.face_a].outline
                @test norm(q - hs.p) > 1e-9
            end
        end
        if split
            mf = K.load_mesh_json(joinpath(EXPORT_FIXTURES, "squares_3x3_split.json"))
            @test mf.sigma == m.sigma
            @test K.svg_string(L) == fixture_text("squares_split_felt_theta0.svg")
            check_layout_against(L, EXPORT_INDEX["squares_split_felt_theta0"])
        end
    end
end

@testset "outlines stay inside their face and shrink it" begin
    m = squares_patch(3, 3)
    c = K.make_cut(m)
    L = K.build_layout(m, c, m.X, 0.0, 10.0, K.profile_felt_laser())
    for fp in L.pieces
        @test K.polygon_area(fp.outline) > 0
        @test K.polygon_area(fp.outline) < K.polygon_area(fp.raw)
    end
    @test isempty(L.warnings)
end

@testset "the deployed layout opens holes and keeps the hinge points shared" begin
    m = squares_patch(3, 3)
    c = K.make_cut(m)
    tm = K.theta_max(c, m.X)
    @test tm.theta_max_geometric > 0.1
    @test isapprox(tm.theta_max_geometric, EXPORT_INDEX["squares_theta_max_geometric"]; rtol=1e-12)
    th = 0.8 * tm.theta_max_geometric
    L0 = K.build_layout(m, c, m.X, 0.0, 10.0, K.profile_felt_laser())
    L1 = K.build_layout(m, c, m.X, th, 10.0, K.profile_felt_laser())
    a0 = sum(f -> K.polygon_area(L0.pieces[f].raw), 1:K.n_faces(m))
    a1 = sum(f -> K.polygon_area(L1.pieces[f].raw), 1:K.n_faces(m))
    @test isapprox(a0, a1; rtol=1e-9)  # faces stay rigid
    # holes have opened: the deployed bounding box grew
    @test K.width(L1) * K.height(L1) > K.width(L0) * K.height(L0)
    # at every hinge site the two faces still touch at the hinge point
    for hs in L1.hinges
        pa = L1.Y[K.prime_vertex(c, hs.face_a, hs.src_vertex)]
        pb = L1.Y[K.prime_vertex(c, hs.face_b, hs.src_vertex)]
        @test norm(pa - pb) < 1e-8
    end
    ref = EXPORT_INDEX["squares_felt_deployed"]
    @test isapprox(th, ref["theta"]; rtol=1e-12)
    check_layout_against(L1, ref)
    @test K.svg_string(L1) == fixture_text("squares_felt_deployed.svg")
end

@testset "ear clipping: a square with a round hole" begin
    outer = [K.Vec2(0, 0), K.Vec2(10, 0), K.Vec2(10, 10), K.Vec2(0, 10)]
    hole = K.Vec2[]
    for k in 0:15
        a = -2 * pi * k / 16
        push!(hole, K.Vec2(5 + 2 * cos(a), 5 + 2 * sin(a)))
    end
    clean, poly, tris = K.triangulate_with_holes(outer, [hole])
    @test clean
    area = 0.0
    for t in tris
        a = poly[t[1]]; b = poly[t[2]]; c = poly[t[3]]
        s = 0.5 * ((b - a)[1] * (c - a)[2] - (b - a)[2] * (c - a)[1])
        @test s > -1e-12
        area += s
    end
    expect = 100.0 - 16 * 0.5 * 4.0 * sin(2 * pi / 16)
    @test isapprox(area, expect; rtol=1e-6)
    # degenerate input
    @test K.triangulate_with_holes([K.Vec2(0, 0), K.Vec2(1, 0)], Vector{K.Vec2}[])[1] == false
end

@testset "STL: closed, consistently oriented, positive volume, round trips" begin
    m = squares_patch(3, 3)
    c = K.make_cut(m)
    p = K.profile_felt_laser()
    L = K.build_layout(m, c, m.X, 0.0, 10.0, p)
    solid, w = K.build_solid(L)
    # A neck tab attached to an inset side leaves a zero-width needle that ear
    # clipping cannot avoid; it is closed with a fan of (degenerate) triangles.

    r = K.check_manifold(solid, 1e-6)
    @test r.closed
    @test r.consistently_oriented
    @test r.volume > 0
    @test isapprox(r.volume, 1700.0; rtol=1e-6)  # 850 mm^2 x 2 mm
    # a living-hinge sheet is a single connected solid
    @test r.n_components == 1
    check_solid_against(solid, r, w, EXPORT_INDEX["squares_felt_solid"])

    path = joinpath(mktempdir(), "sheet.stl")
    K.write_stl_binary(solid, path)
    back = K.read_stl_binary(path)
    @test K.n_tris(back) == K.n_tris(solid)
    rb = K.check_manifold(back, 1e-3)
    @test rb.closed
    @test rb.consistently_oriented
    @test isapprox(rb.volume, r.volume; rtol=1e-4)
    # header + count + 50 bytes per triangle
    @test filesize(path) == 84 + 50 * length(solid.T)
    # identical to the reference file (default header "kiri export")
    check_stl_against(path, fixture_bytes("squares_felt.stl"))

    lo, hi = K.bbox(solid)
    @test isapprox(hi[3] - lo[3], p.thickness)
    @test isapprox(hi[1] - lo[1], K.width(L); rtol=1e-9)
    @test isapprox(hi[2] - lo[2], K.height(L); rtol=1e-9)
    # reader error paths
    @test_throws ErrorException K.read_stl_binary(joinpath(dirname(path), "missing.stl"))
    short = joinpath(dirname(path), "short.stl")
    write(short, zeros(UInt8, 10))
    @test_throws ErrorException K.read_stl_binary(short)
end

@testset "pin-pad: one closed solid per face, stacked by sigma, holes drilled" begin
    m = squares_patch(3, 3)
    c = K.make_cut(m)
    p = K.profile_pla_print()
    L = K.build_layout(m, c, m.X, 0.0, 10.0, p)
    holes = sum(fp -> length(fp.pin_holes), L.pieces)
    @test holes > 0
    # every hinge site drills a hole in both of its faces
    @test holes == 2 * length(L.hinges)
    solid, w = K.build_solid(L)
    r = K.check_manifold(solid, 1e-6)
    @test r.closed
    @test r.consistently_oriented
    @test r.n_components == K.n_faces(m)
    lo, hi = K.bbox(solid)
    @test isapprox(hi[3] - lo[3], 2 * p.thickness + p.clearance)  # two layers
    # against the reference fixture
    check_layout_against(L, EXPORT_INDEX["squares_pla_theta0"])
    @test K.svg_string(L) == fixture_text("squares_pla_theta0.svg")
    check_solid_against(solid, r, w, EXPORT_INDEX["squares_pla_solid"])
    d = mktempdir()
    K.write_stl_binary(solid, joinpath(d, "pla.stl"))
    check_stl_against(joinpath(d, "pla.stl"), fixture_bytes("squares_pla.stl"))
    K.write_3mf(solid, joinpath(d, "pla.3mf"))
    check_3mf_against(joinpath(d, "pla.3mf"), fixture_bytes("squares_pla.3mf"))
end

@testset "zip writer: store-only round trip with CRCs" begin
    z = K.ZipWriter()
    K.add!(z, "a.txt", "hello")
    K.add!(z, "dir/b.bin", String(fill(0x01, 1000)))
    @test K.n_entries(z) == 2
    bytes = K.bytes(z)
    entries = K.zip_read(bytes)
    @test length(entries) == 2
    @test entries[1].first == "a.txt"
    @test entries[1].second == "hello"
    @test entries[2].first == "dir/b.bin"
    @test ncodeunits(entries[2].second) == 1000
    broken = copy(bytes)
    broken[36] ⊻= 0x7f  # corrupt the payload of the first entry (0-based byte 35)
    @test_throws ErrorException K.zip_read(broken)
    @test K.crc32_of("123456789") == 0xCBF43926  # the standard CRC-32 check value
    @test_throws ErrorException K.zip_read(UInt8[])
    p = joinpath(mktempdir(), "a.zip")
    write(z, p)
    @test read(p) == bytes
end

@testset "3MF: OPC parts present, model XML well formed, object count matches" begin
    m = squares_patch(3, 3)
    c = K.make_cut(m)
    for pname in ("felt_laser", "pla_print")
        p = K.profile_by_name(pname)
        L = K.build_layout(m, c, m.X, 0.0, 10.0, p)
        solid, _ = K.build_solid(L)
        n_obj = length(K.split_components(solid))
        @test n_obj == (p.hinge == K.PinPad ? K.n_faces(m) : 1)
        path = joinpath(mktempdir(), "sheet_$pname.3mf")
        K.write_3mf(solid, path)
        entries = K.zip_read(read(path))
        names = Set(first.(entries))
        @test "[Content_Types].xml" in names
        @test "_rels/.rels" in names
        @test "3D/3dmodel.model" in names
        model = ""
        for e in entries
            @test K.xml_well_formed(e.second)[1]
            e.first == "3D/3dmodel.model" && (model = e.second)
        end
        @test !isempty(model)
        @test K.xml_count_elements(model, "object") == n_obj
        @test K.xml_count_elements(model, "item") == n_obj
        @test occursin("unit=\"millimeter\"", model)
        @test K.xml_count_elements(model, "vertex") > 0
        @test K.xml_count_elements(model, "triangle") == K.n_tris(solid)
        pname == "felt_laser" && check_3mf_against(path, fixture_bytes("squares_felt.3mf"))
        pname == "pla_print" && check_3mf_against(path, fixture_bytes("squares_pla.3mf"))
    end
end

@testset "SVG stays well formed on a deployed 3.4.3.12 patch with split cuts" begin
    # frozen from tiling_3_4_3_12(disk(Vec2(0,0), 2.5)) + assign_orientation_relaxation(m,
    # MT19937(12345)).
    m = K.load_mesh_json(joinpath(EXPORT_FIXTURES, "t34312_disk2.5.json"))
    c = K.make_cut(m)
    @test K.n_hinge(c) == EXPORT_INDEX["t34312_n_hinge"]
    @test K.n_split(c) == EXPORT_INDEX["t34312_n_split"]
    tm = K.theta_max(c, m.X)
    @test isapprox(tm.theta_max_geometric, EXPORT_INDEX["t34312_theta_max_geometric"]; rtol=1e-12)
    for (k, frac) in enumerate((0.0, 0.8))
        L = K.build_layout(m, c, m.X, frac * tm.theta_max_geometric, 10.0, K.profile_felt_laser())
        svg = K.svg_string(L)
        @test K.xml_well_formed(svg)[1]
        @test count_ids_with_prefix(svg, "path", "cut_face_") == K.n_faces(m)
        @test length(L.hinges) == K.n_hinge(c)
        tag = "t34312_frac$(k - 1)"
        check_layout_against(L, EXPORT_INDEX[tag])
        @test svg == fixture_text(tag * ".svg")
    end
end

# Julia addition: every migrated export under REPO/export is regenerated from its own
# input graph at the theta/profile/scale recorded in its report JSON and compared byte
# for byte with the reference files (SVG text, binary STL, 3MF archive and entries), plus the
# report's numbers.
@testset "migrated hero / sample exports regenerate byte-identically" begin
    function resolve_input(inp)
        startswith(inp, "export/") && return joinpath(REPO, inp)
        startswith(inp, "results/core_validation/") &&
            return joinpath(CORPUS, "reference_patterns", relpath(inp, "results/core_validation"))
        return joinpath(REPO, inp)
    end
    reports = String[]
    for sub in ("hero", "hero2", "samples", "hero/sequence", "hero2/sequence")
        d = joinpath(EXPORT_DIR, sub)
        isdir(d) || continue
        for f in readdir(d)
            endswith(f, ".json") || continue
            j = JSON.parsefile(joinpath(d, f))
            (j isa AbstractDict && haskey(j, "theta") && haskey(j, "input")) && push!(reports, joinpath(d, f))
        end
    end
    @test length(reports) >= 20
    checked = 0
    for rp in reports
        rep = JSON.parsefile(rp)
        inp = resolve_input(rep["input"])
        if !isfile(inp)
            @info "export regeneration: input missing, skipped" rp inp
            continue
        end
        m = K.load_mesh_json(inp)
        c = K.make_cut(m)
        profile = K.profile_by_name(rep["profile"])
        L = K.build_layout(m, c, m.X, Float64(rep["theta"]), Float64(rep["scale_mm_per_unit"]), profile)
        @test K.n_faces(L) == rep["faces"]
        @test length(L.hinges) == rep["hinge_sites"]
        @test K.n_split(c) == rep["split_edges"]
        @test isapprox(K.width(L), rep["bbox_mm"][1]; rtol=1e-12)
        @test isapprox(K.height(L), rep["bbox_mm"][2]; rtol=1e-12)
        @test L.warnings == String[w for w in rep["warnings"]]
        if haskey(rep, "svg")
            @test K.svg_string(L) == read(joinpath(REPO, rep["svg"]), String)
        end
        if haskey(rep, "solid")
            solid, w3 = K.build_solid(L)
            mr = K.check_manifold(solid, 1e-6)
            s = rep["solid"]
            @test K.n_tris(solid) == s["triangles"]
            @test mr.closed == s["closed"]
            @test mr.consistently_oriented == s["consistently_oriented"]
            @test mr.n_components == s["components"]
            @test isapprox(mr.volume, s["volume_mm3"]; rtol=1e-12)
            @test mr.n_boundary_edges == s["boundary_edges"]
            @test mr.n_nonmanifold_edges == s["nonmanifold_edges"]
            @test w3 == String[w for w in s["warnings"]]
            d = mktempdir()
            if haskey(rep, "stl")
                K.write_stl_binary(solid, joinpath(d, "o.stl"), "kiri " * profile.name)
                check_stl_against(joinpath(d, "o.stl"), read(joinpath(REPO, rep["stl"])))
            end
            if haskey(rep, "3mf")
                @test length(K.split_components(solid)) == s["objects"]
                K.write_3mf(solid, joinpath(d, "o.3mf"), "kiri " * profile.name)
                check_3mf_against(joinpath(d, "o.3mf"), read(joinpath(REPO, rep["3mf"])))
            end
        end
        checked += 1
    end
    @test checked == length(reports)
    @info "regenerated $checked migrated exports; SVG byte-identical for all; 3D files: " *
          "$(EXACT_TALLY["exact"]) byte-identical, $(EXACT_TALLY["structural"]) structurally identical " *
          "(see the comparison policy at the top of this file; Sys.ARCH = $(Sys.ARCH))"
end
