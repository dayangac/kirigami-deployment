# export/stl.jl -- binary STL writer and reader (the reader exists so the tests can
# validate what was actually written, not what we think we wrote).
#
# Little-endian (the STL binary convention on x86/ARM).

"""Writes `M` as a binary STL (80-byte header, uint32 count, 50 bytes per triangle)."""
function write_stl_binary(M::TriMesh, path::AbstractString, header::AbstractString = "kiri export")
    open(path, "w") do f
        head = zeros(UInt8, 80)
        h = codeunits(header)
        nh = min(length(h), 79)
        head[1:nh] .= h[1:nh]
        write(f, head)
        write(f, htol(UInt32(length(M.T))))
        for t in eachindex(M.T)
            n = normal(M, t)
            len = norm(n)
            n = (len > 0) ? n / len : Vec3(0, 0, 0)
            write(f, htol(Float32(n[1])))
            write(f, htol(Float32(n[2])))
            write(f, htol(Float32(n[3])))
            for k in 1:3
                v = M.V[M.T[t][k]]
                write(f, htol(Float32(v[1])))
                write(f, htol(Float32(v[2])))
                write(f, htol(Float32(v[3])))
            end
            write(f, htol(UInt16(0)))
        end
    end
    return nothing
end

"""
Reads a binary STL back. Vertices are NOT welded here, only the triangles are unpacked;
`check_manifold` welds.
"""
function read_stl_binary(path::AbstractString)
    isfile(path) || error("cannot read $path")
    data = read(path)
    length(data) < 84 && error("$path: too short to be a binary STL")
    n = Int(ltoh(reinterpret(UInt32, data[81:84])[1]))
    length(data) < 84 + n * 50 && error("$path: truncated binary STL ($n triangles announced)")
    M = TriMesh()
    sizehint!(M.V, n * 3)
    sizehint!(M.T, n)
    get_f32(p) = Float64(ltoh(reinterpret(Float32, data[p:p+3])[1]))
    for t in 0:n-1
        p = 84 + t * 50 + 12 + 1
        tri = ntuple(3) do k
            q = p + 12 * (k - 1)
            push!(M.V, Vec3(get_f32(q), get_f32(q + 4), get_f32(q + 8)))
            length(M.V)
        end
        push!(M.T, tri)
    end
    return M
end
