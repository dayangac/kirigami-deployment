# export/threemf.jl -- 3MF (core spec) writer: an OPC zip holding
#   [Content_Types].xml, _rels/.rels and 3D/3dmodel.model.
# One <object> per connected solid, so a living-hinge sheet is a single object
# and a pin-pad assembly is one object per face plate.
#
# Vertex indices in the model XML are
# 0-based (3MF format), converted from the 1-based TriMesh at the print site.

# std::ostream << setprecision(9) << double
fnum(v::Float64) = @sprintf("%.9g", v)

"""Splits `M` into connected components (welded on `weld_tol`), one per 3MF object."""
function split_components(M::TriMesh, weld_tol::Float64 = 1e-6)
    if !isempty(M.tri_group)
        order = Dict{Int,Int}()
        out = TriMesh[]
        for t in eachindex(M.T)
            gi = get(order, M.tri_group[t], 0)
            if gi == 0
                push!(out, TriMesh())
                gi = length(out)
                order[M.tri_group[t]] = gi
            end
            p = out[gi]
            tri = ntuple(k -> (push!(p.V, M.V[M.T[t][k]]); length(p.V)), 3)
            push!(p.T, tri)
        end
        # Weld each body's duplicated vertices so the 3MF mesh is indexed properly.
        for p in out
            w = Dict{Key3,Int}()
            V = Vec3[]
            for ti in eachindex(p.T)
                p.T[ti] = ntuple(3) do k
                    key = key_of(p.V[p.T[ti][k]], weld_tol)
                    id = get(w, key, 0)
                    if id == 0
                        push!(V, p.V[p.T[ti][k]])
                        id = length(V)
                        w[key] = id
                    end
                    id
                end
            end
            p.V = V
        end
        return out
    end
    weld = Dict{Key3,Int}()
    id = zeros(Int, length(M.V))
    for i in eachindex(M.V)
        k = key_of(M.V[i], weld_tol)
        w = get(weld, k, 0)
        if w == 0
            w = length(weld) + 1
            weld[k] = w
        end
        id[i] = w
    end
    nv = length(weld)
    dsu = DSU(max(nv, 1))
    for t in M.T
        dsu_join!(dsu, id[t[1]], id[t[2]])
        dsu_join!(dsu, id[t[2]], id[t[3]])
    end
    comp_index = Dict{Int,Int}()
    out = TriMesh[]
    remap = Dict{Int,Int}[]
    for t in M.T
        root = find!(dsu, id[t[1]])
        ci = get(comp_index, root, 0)
        if ci == 0
            push!(out, TriMesh())
            push!(remap, Dict{Int,Int}())
            ci = length(out)
            comp_index[root] = ci
        end
        tri = ntuple(3) do k
            wid = id[t[k]]
            r = get(remap[ci], wid, 0)
            if r == 0
                push!(out[ci].V, M.V[t[k]])
                r = length(out[ci].V)
                remap[ci][wid] = r
            end
            r
        end
        push!(out[ci].T, tri)
    end
    return out
end

"""The 3D/3dmodel.model XML for a list of objects."""
function threemf_model_xml(objects::Vector{TriMesh}, title::AbstractString)
    os = IOBuffer()
    print(os, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",
          "<model unit=\"millimeter\" xml:lang=\"en-US\" ",
          "xmlns=\"http://schemas.microsoft.com/3dmanufacturing/core/2015/02\">\n",
          " <metadata name=\"Title\">", title, "</metadata>\n",
          " <metadata name=\"Application\">kiri_export</metadata>\n",
          " <resources>\n")
    for (i, M) in enumerate(objects)
        print(os, "  <object id=\"", i, "\" type=\"model\">\n   <mesh>\n    <vertices>\n")
        for v in M.V
            print(os, "     <vertex x=\"", fnum(v[1]), "\" y=\"", fnum(v[2]), "\" z=\"",
                  fnum(v[3]), "\"/>\n")
        end
        print(os, "    </vertices>\n    <triangles>\n")
        for t in M.T
            # 1-based TriMesh -> 0-based 3MF indices
            print(os, "     <triangle v1=\"", t[1] - 1, "\" v2=\"", t[2] - 1, "\" v3=\"", t[3] - 1, "\"/>\n")
        end
        print(os, "    </triangles>\n   </mesh>\n  </object>\n")
    end
    print(os, " </resources>\n <build>\n")
    for i in eachindex(objects)
        print(os, "  <item objectid=\"", i, "\"/>\n")
    end
    print(os, " </build>\n</model>\n")
    return String(take!(os))
end

"""Writes `M` as a 3MF package (STORE-only zip with the three OPC parts)."""
function write_3mf(M::TriMesh, path::AbstractString, title::AbstractString = "kiri export")
    objects = split_components(M)
    z = ZipWriter()
    add!(z, "[Content_Types].xml",
         "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" *
         "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">\n" *
         " <Default Extension=\"rels\" " *
         "ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>\n" *
         " <Default Extension=\"model\" " *
         "ContentType=\"application/vnd.ms-package.3dmanufacturing-3dmodel+xml\"/>\n" *
         "</Types>\n")
    add!(z, "_rels/.rels",
         "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" *
         "<Relationships " *
         "xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\n" *
         " <Relationship Id=\"rel0\" Target=\"/3D/3dmodel.model\" " *
         "Type=\"http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel\"/>\n" *
         "</Relationships>\n")
    add!(z, "3D/3dmodel.model", threemf_model_xml(objects, title))
    write(z, path)
    return nothing
end
