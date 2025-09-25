module ArchimedIO

using Meshes
using PlantGeom
using MultiScaleTreeGraph
using ..ArchimedLight: Scene, Component, add_component

"""
    load_opf_scene(file::AbstractString; group="plant")

Load an OPF into an ArchimedLight Scene. Each node with geometry becomes a component
named `node_<id>` using its transformed `SimpleMesh`.
"""
function load_opf_scene(file::AbstractString; group::String="plant")
    mtg = PlantGeom.read_opf(file)
    scene = Scene()
    # Traverse nodes that have geometry and convert to SimpleMesh
    meshes = MultiScaleTreeGraph.traverse(mtg, filter_fun = x -> haskey(x, :geometry), type = SimpleMesh) do node
        PlantGeom.refmesh_to_mesh(node)
    end
    # Names for components: use node_id in traversal order
    ids = MultiScaleTreeGraph.traverse(mtg, filter_fun = x -> haskey(x, :geometry), type = Int) do node
        MultiScaleTreeGraph.node_id(node)
    end
    for (i, m) in enumerate(meshes)
        add_component(scene, "node_$(ids[i])", m; group=group)
    end
    return scene
end

end # module

