# export/zip.jl -- a minimal STORE-only (compression method 0) zip writer, enough for
# the OPC container that a 3MF file is. No external libraries, per the spec.
#
# The archive is a Vector{UInt8}; entry names and contents are Strings (byte strings).

const CRC_TABLE = let t = zeros(UInt32, 256)
    for i in 0:255
        c = UInt32(i)
        for _ in 1:8
            c = (c & 1) != 0 ? (0xEDB88320 ⊻ (c >> 1)) : (c >> 1)
        end
        t[i + 1] = c
    end
    t
end

"""CRC-32 (IEEE) of a byte string."""
function crc32_of(data::AbstractVector{UInt8})
    c = 0xFFFFFFFF
    for ch in data
        c = CRC_TABLE[((c ⊻ ch) & 0xff) + 1] ⊻ (c >> 8)
    end
    return c ⊻ 0xFFFFFFFF
end
crc32_of(data::AbstractString) = crc32_of(codeunits(data))

struct ZipEntry
    name::String
    data::String
    crc::UInt32
end

mutable struct ZipWriter
    entries_::Vector{ZipEntry}
end
ZipWriter() = ZipWriter(ZipEntry[])
n_entries(z::ZipWriter) = length(z.entries_)

"""`name` uses forward slashes and no leading slash, as OPC requires."""
add!(z::ZipWriter, name::AbstractString, content::AbstractString) =
    push!(z.entries_, ZipEntry(String(name), String(content), crc32_of(content)))

function put16!(s::Vector{UInt8}, v::Integer)
    push!(s, UInt8(v & 0xff))
    push!(s, UInt8((v >> 8) & 0xff))
end
function put32!(s::Vector{UInt8}, v::Integer)
    put16!(s, v & 0xffff)
    put16!(s, (v >> 16) & 0xffff)
end

"""The complete archive."""
function bytes(z::ZipWriter)
    out = UInt8[]
    offsets = UInt32[]
    for e in z.entries_
        push!(offsets, UInt32(length(out)))
        put32!(out, 0x04034b50)                 # local file header
        put16!(out, 20)                         # version needed
        put16!(out, 0)                          # flags
        put16!(out, 0)                          # method 0 = stored
        put16!(out, 0)                          # mod time
        put16!(out, 0x21)                       # mod date (1980-01-01)
        put32!(out, e.crc)
        put32!(out, ncodeunits(e.data))         # compressed size
        put32!(out, ncodeunits(e.data))         # uncompressed size
        put16!(out, ncodeunits(e.name))
        put16!(out, 0)                          # extra length
        append!(out, codeunits(e.name))
        append!(out, codeunits(e.data))
    end
    cd_start = UInt32(length(out))
    for (i, e) in enumerate(z.entries_)
        put32!(out, 0x02014b50)  # central directory header
        put16!(out, 20)          # version made by
        put16!(out, 20)          # version needed
        put16!(out, 0)
        put16!(out, 0)
        put16!(out, 0)
        put16!(out, 0x21)
        put32!(out, e.crc)
        put32!(out, ncodeunits(e.data))
        put32!(out, ncodeunits(e.data))
        put16!(out, ncodeunits(e.name))
        put16!(out, 0)  # extra
        put16!(out, 0)  # comment
        put16!(out, 0)  # disk number
        put16!(out, 0)  # internal attributes
        put32!(out, 0)  # external attributes
        put32!(out, offsets[i])
        append!(out, codeunits(e.name))
    end
    cd_size = UInt32(length(out)) - cd_start
    put32!(out, 0x06054b50)  # end of central directory
    put16!(out, 0)
    put16!(out, 0)
    put16!(out, length(z.entries_))
    put16!(out, length(z.entries_))
    put32!(out, cd_size)
    put32!(out, cd_start)
    put16!(out, 0)  # comment length
    return out
end

# Extends Base.write so `write(z, path)` reads naturally without shadowing Base.write for
# the rest of the module.
function Base.write(z::ZipWriter, path::AbstractString)
    open(path, "w") do io
        write(io, bytes(z))
    end
    return nothing
end

"""
Reads back a STORE-only archive by walking the local file headers, verifying
each CRC. Throws on a corrupt archive. Returns `[name => data, ...]`.
"""
function zip_read(b::AbstractVector{UInt8})
    out = Pair{String,String}[]
    nb = length(b)
    rd16(p) = UInt32(b[p + 1]) | (UInt32(b[p + 2]) << 8)   # p is a 0-based offset
    rd32(p) = rd16(p) | (rd16(p + 2) << 16)
    i = 0
    while i + 30 <= nb && rd32(i) == 0x04034b50
        method = rd16(i + 8)
        crc = rd32(i + 14)
        csize = rd32(i + 18)
        usize = rd32(i + 22)
        nlen = rd16(i + 26)
        elen = rd16(i + 28)
        method != 0 && error("zip_read: entry is not STOREd")
        i + 30 + nlen + elen + csize > nb && error("zip_read: truncated")
        name = String(b[i + 31:i + 30 + nlen])
        data = String(b[i + 30 + nlen + elen + 1:i + 30 + nlen + elen + csize])
        csize != usize && error("zip_read: size mismatch for $name")
        crc32_of(data) != crc && error("zip_read: CRC mismatch for $name")
        push!(out, name => data)
        i += 30 + nlen + elen + csize
    end
    isempty(out) && error("zip_read: no local file header at offset 0")
    (i + 4 > nb || rd32(i) != 0x02014b50) &&
        error("zip_read: central directory does not follow the entries")
    return out
end
zip_read(s::AbstractString) = zip_read(collect(codeunits(s)))
