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
export ValidationReport
export LightSceneBuilder
export LightSimulation

export read_scene
export prepare_scene
export write_scene
export add_ground!
export light_scene
export add_object!
export add_plant!
export read_models
export prepare_models
export models_for
export translucent
export virtual_sensor
export emitter
export read_options
export read_simulation
export read_meteo
export prepare_meteo
export cache_summary
export update_scene!
export update_models!
export update_options!
export check_scene
export check_models
export check_meteo
export check_simulation
export run_light
export attach_node_values!
export attach_light_step!
export attach_light_series!
export light_render_geometry
export tile_light_geometry
export light_metric_values
export light_face_values
export light_vertex_values
export lightplot
export lightplot!

end
