using Meshes

struct Component
    name::String
    mesh::SimpleMesh3
    group::String
end

mutable struct Scene
    components::Vector{Component}
    bbox::Box
    toric::Bool
end

function Scene(; toric::Bool=false)
    # initialize with a degenerate box at origin
    b = Box(Point(0.0,0.0,0.0), Point(0.0,0.0,0.0))
    return Scene(Component[], b, toric)
end

function _update_bbox(scene::Scene)
    if isempty(scene.components)
        scene.bbox = Box(Point(0.0,0.0,0.0), Point(0.0,0.0,0.0))
        return
    end
    # merge all component bounding boxes
    scene.bbox = boundingbox(getfield.(scene.components, :mesh))
end

function add_component(scene::Scene, name::String, mesh::SimpleMesh3; group::String="default")
    push!(scene.components, Component(name, mesh, group))
    _update_bbox(scene)
    return scene
end
