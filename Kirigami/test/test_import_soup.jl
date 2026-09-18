# test_import_soup.jl -- the welding importer, case by case.
include("helpers.jl")

const K = Kirigami

# Two unit squares side by side, given as INDEPENDENT polygons, so the shared edge's two
# corners are each listed twice and disagree by 1e-9 -- the soup case.
function two_squares()
    e = 1e-9
    return [[K.Vec2(0, 0), K.Vec2(1, 0), K.Vec2(1, 1), K.Vec2(0, 1)],
            [K.Vec2(1 + e, -e), K.Vec2(2, 0), K.Vec2(2, 1), K.Vec2(1 - e, 1 + e)]]
end

n_interior_edges(m::K.Mesh) = count(e -> e.n_faces == 2, m.edges)

@testset "weld_polygons fuses duplicated corners into one planar mesh" begin
    r = K.weld_polygons(two_squares(), 1e-6)
    @test r.ok
    @test r.status == "ok"
    @test K.n_faces(r.mesh) == 2
    @test K.n_vertices(r.mesh) == 6      # 8 corners - 2 welded pairs
    @test r.n_segments == 7              # 8 polygon edges - 1 shared
    @test isapprox(r.total_area, 2.0; rtol = 1e-6)
    # The shared edge is interior: exactly one edge of the welded mesh has two faces.
    @test n_interior_edges(r.mesh) == 1
end

@testset "a tolerance below the corner gap leaves the squares unwelded" begin
    r = K.weld_polygons(two_squares(), 1e-12)
    @test r.ok
    @test K.n_faces(r.mesh) == 2
    @test K.n_vertices(r.mesh) == 8   # nothing fused
    @test n_interior_edges(r.mesh) == 0
end

@testset "segments in arbitrary order still walk out the bounded faces" begin
    # The same two squares, edge by edge, deliberately shuffled and with the shared edge
    # supplied twice in opposite directions.
    s = [(K.Vec2(1, 1), K.Vec2(0, 1)), (K.Vec2(2, 0), K.Vec2(2, 1)), (K.Vec2(0, 0), K.Vec2(1, 0)),
         (K.Vec2(1, 0), K.Vec2(1, 1)), (K.Vec2(1, 1), K.Vec2(2, 1)), (K.Vec2(0, 1), K.Vec2(0, 0)),
         (K.Vec2(1, 1), K.Vec2(1, 0)), (K.Vec2(1, 0), K.Vec2(2, 0))]
    r = K.weld_segments(s, 1e-9)
    @test r.ok
    @test r.n_segments == 7           # the duplicated shared edge collapses
    @test K.n_faces(r.mesh) == 2
    @test K.n_vertices(r.mesh) == 6
    @test isapprox(r.total_area, 2.0; rtol = 1e-9)
end

@testset "a dangling segment is pruned and bounds no face" begin
    s = [(K.Vec2(0, 0), K.Vec2(1, 0)), (K.Vec2(1, 0), K.Vec2(1, 1)), (K.Vec2(1, 1), K.Vec2(0, 1)),
         (K.Vec2(0, 1), K.Vec2(0, 0)), (K.Vec2(1, 1), K.Vec2(2, 2))]   # the spur
    r = K.weld_segments(s, 1e-9)
    @test r.ok
    @test r.n_pruned == 1
    @test K.n_faces(r.mesh) == 1
    @test K.n_vertices(r.mesh) == 4
    @test isapprox(r.total_area, 1.0; rtol = 1e-9)
    for f in r.mesh.faces
        @test length(f) == 4
    end
end

@testset "a square with a square hole gives the annulus' four faces, not the hole" begin
    # Concentric squares connected by nothing: the inner square's interior is a bounded
    # face too, so both walks are kept and the outer walk is dropped.
    p = [[K.Vec2(0, 0), K.Vec2(3, 0), K.Vec2(3, 3), K.Vec2(0, 3)],
         [K.Vec2(1, 1), K.Vec2(2, 1), K.Vec2(2, 2), K.Vec2(1, 2)]]
    r = K.weld_polygons(p, 1e-9)
    @test r.ok
    @test K.n_faces(r.mesh) == 2
    @test isapprox(r.total_area, 10.0; rtol = 1e-9)  # 9 + 1
end

@testset "an empty soup reports a status instead of throwing" begin
    r = K.weld_polygons(Vector{K.Vec2}[], 1e-9)
    @test !r.ok
    @test r.status == "no_segments"
end

# The SVG reader and the area floor.
@testset "read_svg_segments reads line/polyline/polygon and skips path" begin
    svg = """<svg xmlns="http://www.w3.org/2000/svg">
      <line x1="0" y1="0" x2="1" y2="0" stroke="black"/>
      <line x1="1" y1="0" x2="1" y2="1"/>
      <polyline points="1,1 0,1 0,0"/>
      <polygon points="2 0, 3 0, 3 1"/>
      <path d="M 0 0 A 1 1 0 0 1 1 1"/>
      <line x1="5" y1="5"/>
    </svg>"""
    path = tempname() * ".svg"
    write(path, svg)
    segs, n_path = K.read_svg_segments(path)
    @test n_path == 1
    @test length(segs) == 2 + 2 + 3      # incomplete <line> ignored
    @test segs[3] == (K.Vec2(1, 1), K.Vec2(0, 1))
    @test segs[end] == (K.Vec2(3, 1), K.Vec2(2, 0))   # polygon closes
    r = K.import_svg_soup(path, 1e-6)
    @test r.ok
    @test K.n_faces(r.mesh) == 2
    @test isapprox(r.total_area, 1.5; rtol = 1e-9)
    @test r.n_faces_dropped == 2      # one outer walk per connected component
    # area floor drops the small triangle
    r2 = K.import_svg_soup(path, 1e-6, 0.4)
    @test r2.ok && K.n_faces(r2.mesh) == 1 && r2.n_faces_dropped == 3
    rm(path)
    segs0, _ = K.read_svg_segments(path * ".missing")
    @test isempty(segs0)
    @test K.import_svg_soup(path * ".missing").status == "no_segments_in_svg"
end
