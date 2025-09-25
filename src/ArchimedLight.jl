module ArchimedLight

export Component, Scene,
       Band, OpticalProps, SkyConfig, InterceptionConfig,
       InterceptionResult,
       add_component, set_optics!, compute_interception

using Meshes

include("geometry.jl")
include("scene.jl")
include("optics.jl")
include("sky.jl")
include("intercept.jl")

# Optional IO layer backed by PlantGeom + MultiScaleTreeGraph
include("io/opf.jl")

# Config loader and runner
include("skysets.jl")
include("config.jl")
include("runner.jl")

end # module
