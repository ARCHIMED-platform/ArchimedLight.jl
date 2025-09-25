module ArchimedIO

using Meshes
using PlantGeom
using MultiScaleTreeGraph
using ..ArchimedLight: Scene, Component, add_component

function _sanitize_ops(path::AbstractString)
    text = replace(read(path, String), "\r\n" => "\n")
    lines = split(text, '\n')
    sanitized = String[]
    for line in lines
        if occursin(".opf", line)
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

function _read_ops(scene_path::AbstractString)
    scene_mtg = PlantGeom.read_ops(scene_path)
    if _count_geometry_nodes(scene_mtg) > 0
        return scene_mtg
    end
    sanitized = _sanitize_ops(scene_path)
    tmp_path, tmp_io = mktemp()
    write(tmp_io, sanitized)
    close(tmp_io)
    scene_mtg = PlantGeom.read_ops(tmp_path)
    rm(tmp_path; force=true)
    return scene_mtg
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
    load_ops_scene(file::AbstractString)

Load an OPS scene: reads all OPF instances with their scene transforms applied,
and returns an ArchimedLight Scene aggregating all component meshes.
"""
function load_ops_scene(file::AbstractString)
    scene_mtg = _read_ops(file)
    scene = Scene()
    MultiScaleTreeGraph.traverse(scene_mtg; all=true) do node
        if haskey(node, :geometry) && node[:geometry] !== nothing
            mesh = PlantGeom.refmesh_to_mesh(node)
            mesh === nothing && return
            name = "node_$(MultiScaleTreeGraph.node_id(node))"
            group = _functional_group(node)
            ctype = String(MultiScaleTreeGraph.symbol(node))
            add_component(scene, name, mesh; group=group, ctype=ctype)
        end
    end
    return scene
end

end # module
