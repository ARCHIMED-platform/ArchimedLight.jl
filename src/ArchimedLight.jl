module ArchimedLight

include("types.jl")
include("io.jl")
include("sky.jl")
include("turtle.jl")
include("interception.jl")
include("scattering.jl")
include("pipeline.jl")
include("output.jl")

export LightConfig
export MeteoTable
export SceneGeometry
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
export read_light_config
export read_meteo
export compute_sky
export build_turtle
export compute_directional_fluxes
export compute_first_order
export build_scattering_transfer_graph
export compute_scattering
export integrate_light
export run_light_step
export run_light_series
export component_variable_names
export component_values_table
export write_component_values_csv
export scene_variable_names
export scene_values_table
export write_scene_values_csv

end
