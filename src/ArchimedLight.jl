module ArchimedLight

import PlantGeom
import MultiScaleTreeGraph
import GeometryBasics
import StaticArrays
import LinearAlgebra: norm, dot, cross
import Serialization
import Dates
import OrderedCollections: OrderedDict
import PlantMeteo
import Tables
import YAML
import CSV

include("types.jl")
include("io.jl")
include("sky.jl")
include("turtle.jl")
include("interception.jl")
include("scattering.jl")
include("pipeline.jl")
include("attach.jl")
include("outputs.jl")
include("visualization.jl")

export MeteoTable
export LightRenderGeometry
export LightModels
export GroupModel
export TypeModel
export InterceptionModel
export EmitterModel
export OpticalProperties
export LightOptions
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
export SceneSummary
export MeteoSummary
export LightSimulation

export read_scene
export write_scene
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
export summarize_scene
export summarize_meteo
export run_light
export attach_node_values!
export attach_light_step!
export attach_light_series!
export write_component_values
export light_render_geometry
export tile_light_geometry
export light_metric_values
export light_face_values
export light_vertex_values
export lightplot
export lightplot!

end
