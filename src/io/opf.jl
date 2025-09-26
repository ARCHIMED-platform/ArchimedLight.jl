module ArchimedIO

using Meshes
using PlantGeom
using MultiScaleTreeGraph
using StaticArrays
using Rotations
using Unitful
using LinearAlgebra: I
using ..ArchimedLight: Scene, Component, add_component

strip_length(x) = x isa Unitful.Quantity ? Unitful.ustrip(x) : Float64(x)

function _sanitize_ops(path::AbstractString)
    text = replace(read(path, String), "\r\n" => "\n")
    lines = split(text, '\n')
    sanitized = String[]
    for line in lines
        if occursin(".opf", line) || occursin(".gwa", line)
            parts = split(line, '\t'; keepempty=true)
            if length(parts) >= 3
                head = parts[1:3]
                numeric = parts[4:end]
                clean = String[]
                for (idx, value) in enumerate(numeric)
                    if idx <= 7
                        if isempty(value)
                            push!(clean, idx == 4 ? "1" : "0")
                        else
                            push!(clean, value)
                        end
                    end
                end
                while length(clean) < 7
                    push!(clean, "0")
                end
                push!(sanitized, join(vcat(head, clean[1:7]), '\t'))
                continue
            end
        end
        push!(sanitized, line)
    end
    return join(sanitized, "\n")
end

function _count_geometry_nodes(scene_mtg)
    c = 0
    MultiScaleTreeGraph.traverse(scene_mtg; all=true) do node
        if haskey(node, :geometry) && node[:geometry] !== nothing
            c += 1
        end
    end
    return c
end

function _ops_entries(scene_path::AbstractString)
    data = PlantGeom.read_ops_file(scene_path)
    entries = data.object_table
    dims = data.scene_dimensions
    if !isempty(entries)
        return entries, dims
    end
    sanitized = _sanitize_ops(scene_path)
    parent = dirname(scene_path)
    entries = mktemp(parent) do tmp_path, tmp_io
        write(tmp_io, sanitized)
        close(tmp_io)
        info = PlantGeom.read_ops_file(tmp_path)
        info.object_table
    end
    if !isempty(entries)
        return entries, dims
    end
    return _parse_ops_fallback(sanitized, dims)
end

function _parse_ops_fallback(contents::AbstractString, default_dims)
    entries = NamedTuple[]
    functional_group = "scene"
    scene_dims = default_dims
    for raw in split(contents, '\n')
        stripped = strip(raw)
        isempty(stripped) && continue
        if startswith(stripped, "T ")
            parts = split(stripped)
            if length(parts) >= 6
                x0 = parse(Float64, parts[2])
                y0 = parse(Float64, parts[3])
                z0 = parse(Float64, parts[4])
                x1 = parse(Float64, parts[5])
                y1 = parse(Float64, parts[6])
                scene_dims = (Meshes.Point(x0, y0, z0), Meshes.Point(x1, y1, z0))
            end
            continue
        end
        if startswith(stripped, "#[Archimed]")
            functional_group = strip(replace(stripped, "#[Archimed]" => ""))
            continue
        end
        startswith(stripped, "#") && continue
        parts = split(raw, '\t'; keepempty=true)
        length(parts) < 3 && continue
        file_path = strip(parts[3])
        isempty(file_path) && continue
        scene_id = try parse(Int, strip(parts[1])) catch; continue; end
        plant_id = try parse(Int, strip(parts[2])) catch; continue; end
        parsefloat(idx, default) = begin
            if idx > length(parts)
                return default
            end
            val = strip(parts[idx])
            isempty(val) && return default
            try
                parse(Float64, val)
            catch
                default
            end
        end
        x = parsefloat(4, 0.0)
        y = parsefloat(5, 0.0)
        z = parsefloat(6, 0.0)
        scale = parsefloat(7, 1.0)
        az = parsefloat(8, 0.0)
        angle = parsefloat(9, 0.0)
        rotation = parsefloat(10, 0.0)
        pos = Meshes.Point(x*Unitful.u"m", y*Unitful.u"m", z*Unitful.u"m")
        push!(entries, (; sceneID=scene_id, plantID=plant_id, filePath=file_path, pos, scale, inclinationAzimut=az, inclinationAngle=angle, rotation, functional_group))
    end
    return entries, scene_dims
end

struct OpsEntry
    scene_id::Int
    plant_id::Int
    file_path::String
    functional_group::String
    position::SVector{3,Float64}
    scale::Float64
    inclination_azimut::Float64
    inclination_angle::Float64
    rotation::Float64
end

deg2rad_f(x::Float64) = x * (π / 180.0)

function _rotation_matrix(entry::OpsEntry)
    R = SMatrix{3,3}(I)
    if entry.rotation != 0.0
        R = SMatrix{3,3}(Matrix(RotZ(deg2rad_f(entry.rotation)))) * R
    end
    if entry.inclination_angle != 0.0
        R = SMatrix{3,3}(Matrix(RotX(deg2rad_f(entry.inclination_angle)))) * R
    end
    if entry.inclination_azimut != 0.0
        R = SMatrix{3,3}(Matrix(RotZ(deg2rad_f(entry.inclination_azimut)))) * R
    end
    return R
end

function _entry_from_row(scene_path::AbstractString, row)
    pos_vec = SVector(strip_length.(Tuple(Meshes.to(row.pos)))...)
    scale = row.scale == 0.0 ? 1.0 : row.scale
    return OpsEntry(
        row.sceneID,
        row.plantID,
        isabspath(row.filePath) ? row.filePath : joinpath(dirname(scene_path), row.filePath),
        String(row.functional_group),
        pos_vec,
        scale,
        row.inclinationAzimut,
        row.inclinationAngle,
        row.rotation,
    )
end

function _prepare_entries(scene_path::AbstractString)
    rows, dims = _ops_entries(scene_path)
    return OpsEntry[_entry_from_row(scene_path, row) for row in rows], dims
end

function _generate_paving!(scene::Scene, dims, paving_specs::Dict{Tuple{String,String},Int})
    isempty(paving_specs) && return
    origin = dims[1]
    maxp = dims[2]
    o = (strip_length(to(origin)[1]), strip_length(to(origin)[2]), strip_length(to(origin)[3]))
    span = (strip_length(to(maxp)[1] - to(origin)[1]), strip_length(to(maxp)[2] - to(origin)[2]))
    total_area = span[1] * span[2]
    for ((group, ctype), cells) in paving_specs
        cells <= 0 && continue
        target_cell_area = total_area / cells
        cell_side = sqrt(target_cell_area)
        nx = max(1, floor(Int, span[1] / cell_side))
        ny = max(1, floor(Int, span[2] / cell_side))
        nx = max(nx, 1)
        ny = max(ny, 1)
        cell_w = span[1] / nx
        cell_h = span[2] / ny
        z = o[3]
        name_idx = 1
        for ix in 0:nx-1
            x0 = o[1] + ix * cell_w
            x1 = x0 + cell_w
            for iy in 0:ny-1
                y0 = o[2] + iy * cell_h
                y1 = y0 + cell_h
                pts = [
                    Point(x0 * Unitful.u"m", y0 * Unitful.u"m", z * Unitful.u"m"),
                    Point(x1 * Unitful.u"m", y0 * Unitful.u"m", z * Unitful.u"m"),
                    Point(x1 * Unitful.u"m", y1 * Unitful.u"m", z * Unitful.u"m"),
                    Point(x0 * Unitful.u"m", y1 * Unitful.u"m", z * Unitful.u"m")
                ]
                connec = connect.([(1, 2, 3), (1, 3, 4)])
                mesh = SimpleMesh(pts, connec)
                name = "paving_$(ctype)_$(ix)_$(iy)_$(name_idx)"
                add_component(scene, name, mesh; group=group, ctype=ctype)
                name_idx += 1
            end
        end
    end
end

function _transform_mesh(mesh::SimpleMesh, entry::OpsEntry)
    verts = mesh.vertices
    scale = entry.scale
    rot = _rotation_matrix(entry)
    trans = entry.position
    coords = Vector{NTuple{3,Any}}(undef, length(verts))
    for (i, p) in enumerate(verts)
        base = SVector(strip_length.(Tuple(Meshes.to(p)))...)
        transformed = rot * (base .* scale) .+ trans
        coords[i] = (transformed[1]*Unitful.u"m", transformed[2]*Unitful.u"m", transformed[3]*Unitful.u"m")
    end
    return SimpleMesh(coords, mesh.topology)
end

function _load_gwa_meshes(path::AbstractString)
    text = replace(read(path, String), '\t' => ' ')
    meshes = SimpleMesh[]
    for mesh_match in eachmatch(r"<mesh[^>]*>(.*?)</mesh>"s, text)
        block = mesh_match.captures[1]
        pts_match = match(r"<points>(.*?)</points>"s, block)
        pts_match === nothing && continue
        coords = parse.(Float64, split(strip(pts_match.captures[1])))
        npts = length(coords) ÷ 3
        points = Vector{Point}(undef, npts)
        for i in 1:npts
            idx = 3*(i-1)
            points[i] = Point(coords[idx+1]*Unitful.u"m", coords[idx+2]*Unitful.u"m", coords[idx+3]*Unitful.u"m")
        end
        faces = NTuple{3,Int}[]
        for face_match in eachmatch(r"<face[^>]*>(.*?)</face>"s, block)
            verts = parse.(Int, split(strip(face_match.captures[1])))
            length(verts) == 3 || continue
            push!(faces, (verts[1]+1, verts[2]+1, verts[3]+1))
        end
        topology = Meshes.connect.(faces)
        push!(meshes, SimpleMesh(points, topology))
    end
    return meshes
end

function _functional_group(node)
    current = node
    while current !== nothing
        if hasproperty(current, :functional_group)
            fg = getfield(current, :functional_group)
            if fg !== nothing
                return String(fg)
            end
        end
        parent = try
            MultiScaleTreeGraph.parent(current)
        catch
            nothing
        end
        current = parent
    end
    return "scene"
end

"""
    load_opf_scene(file::AbstractString; group="plant")

Load an OPF into an ArchimedLight Scene. Each node with geometry becomes a component
named `node_<id>` using its transformed `SimpleMesh`.
"""
function load_opf_scene(file::AbstractString; group::String="plant")
    mtg = PlantGeom.read_opf(file)
    scene = Scene()
    meshes = MultiScaleTreeGraph.traverse(mtg, filter_fun = x -> haskey(x, :geometry), type = SimpleMesh) do node
        PlantGeom.refmesh_to_mesh(node)
    end
    ids = MultiScaleTreeGraph.traverse(mtg, filter_fun = x -> haskey(x, :geometry), type = Int) do node
        MultiScaleTreeGraph.node_id(node)
    end
    for (i, m) in enumerate(meshes)
        add_component(scene, "node_$(ids[i])", m; group=group)
    end
    return scene
end

"""
    load_ops_scene(file::AbstractString; paving_specs=Dict{Tuple{String,String},Int}())

Load an OPS scene: reads all OPF/GWA instances with their scene transforms applied,
and returns an ArchimedLight Scene aggregating all component meshes. When
`paving_specs` defines `plot_paving` counts for a `(group, type)` pair, the scene
rectangle is tessellated accordingly to create background paving components.
"""
function load_ops_scene(file::AbstractString; paving_specs::Dict{Tuple{String,String},Int}=Dict{Tuple{String,String},Int}())
    entries, scene_dims = _prepare_entries(file)
    scene = Scene()
    opf_cache = Dict{Tuple{String,String},Vector{Component}}()
    for entry in entries
        _, ext = splitext(entry.file_path)
        ext = lowercase(ext)
        if ext == ".opf"
            key = (entry.file_path, entry.functional_group)
            templates = get!(opf_cache, key) do
                plant_scene = load_opf_scene(entry.file_path; group=entry.functional_group)
                copy(plant_scene.components)
            end
            for comp in templates
                mesh = _transform_mesh(comp.mesh, entry)
                name = "scene$(entry.scene_id)_plant$(entry.plant_id)_$(comp.name)"
                add_component(scene, name, mesh; group=comp.group, ctype=comp.ctype)
            end
        elseif ext == ".gwa"
            for (idx, mesh) in enumerate(_load_gwa_meshes(entry.file_path))
                tmesh = _transform_mesh(mesh, entry)
                name = "scene$(entry.scene_id)_plant$(entry.plant_id)_gwa$(idx)"
                add_component(scene, name, tmesh; group=entry.functional_group, ctype="GWA")
            end
        else
            @warn "Unsupported OPS asset" entry.file_path
        end
    end
    _generate_paving!(scene, scene_dims, paving_specs)
    return scene
end

end # module
