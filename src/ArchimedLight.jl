module ArchimedLight

export Component, Scene,
       Band, OpticalProps, SkyConfig, InterceptionConfig,
       add_component, set_optics!, compute_interception

using Meshes

include("geometry.jl")
include("scene.jl")
include("optics.jl")
include("sky.jl")
include("intercept.jl")

# Optional IO layer backed by PlantGeom + MultiScaleTreeGraph
include("io/opf.jl")

end # module
