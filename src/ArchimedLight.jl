module ArchimedLight

include("types.jl")
include("io.jl")
include("sky.jl")
include("turtle.jl")
include("interception.jl")
include("scattering.jl")
include("pipeline.jl")
include("attach.jl")
include("visualization.jl")

export MeteoTable
export SceneNodeData
export SceneGeometry
export LightRenderGeometry
export LightModels
export GroupModel
export TypeModel
export InterceptionModel
export EmitterModel
export OpticalProperties
export LightOptions
export scene_node
export scene_node_ids
export node_areas
export node_barycenters
export InterceptionBackend, RasterCPUBackend
export ScatteringBackend, RaycastScatteringBackend, LinksScatteringBackend
export SkyState
export TurtleSector, TurtleGrid
export DirectionalFluxes
export FirstOrderResult
export ScatteringResult
export ScatteringTransferGraph
export LightBudget
export LightStepResult

export read_scene
export prepare_scene
export write_scene
export add_ground!
export read_models
export prepare_models
export read_options
export read_config
export read_meteo
export prepare_meteo
export compute_sky
export build_turtle
export compute_directional_fluxes
export compute_first_order
export build_scattering_transfer_graph
export compute_scattering
export integrate_light
export run_light_step
export run_light_series
export attach_node_values!
export attach_light_step!
export attach_light_series!
export light_render_geometry
export light_metric_values
export light_face_values
export light_vertex_values
export lightplot
export lightplot!

end
