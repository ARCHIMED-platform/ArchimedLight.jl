module ArchimedLight

include("types.jl")
include("io.jl")
include("sky.jl")
include("turtle.jl")
include("interception.jl")
include("scattering.jl")
include("pipeline.jl")

export LightConfig
export MeteoTable
export SceneGeometry
export SkyState
export TurtleSector, TurtleGrid
export DirectionalFluxes
export FirstOrderResult
export ScatteringResult
export LightBudget
export LightStepResult

export read_scene
export read_light_config
export read_meteo
export compute_sky
export build_turtle
export compute_directional_fluxes
export compute_first_order
export compute_scattering
export integrate_light
export run_light_step
export run_light_series

end
