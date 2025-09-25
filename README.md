# ArchimedLight.jl (prototype)

A minimal Julia reimplementation of ARCHIMED's light interception core.

Status: prototype for first-order (no multiple scattering) diffuse sky interception over a triangle mesh scene.

Quick start

```julia
using Pkg
Pkg.activate(".")
Pkg.develop(path=".") # if testing in place
using ArchimedLight

# Build a single square component
scene = Scene()
mesh = Mesh([Triangle(Vec3(-0.5, -0.5, 0.0), Vec3(0.5, -0.5, 0.0), Vec3(0.5, 0.5, 0.0)),
             Triangle(Vec3(-0.5, -0.5, 0.0), Vec3(0.5, 0.5, 0.0), Vec3(-0.5, 0.5, 0.0))])
add_component(scene, "leaf", mesh)
set_optics!("leaf", OpticalProps())

sky = SkyConfig(count=16, PAR_Wm2=400.0, NIR_Wm2=400.0)
cfg = InterceptionConfig(pixel_size=0.25, scattering=false)
res = compute_interception(scene, sky, cfg)
println(res)
```

Roadmap
- Multiple scattering (order-N) with reflectance/transmittance.
- Sky discretizations matching 1/6/16/46/136/406 sector sets.
- Solar position + direct beam handling.
- Scene/format IO (OPS/OPF or a neutral format).
- Acceleration structures (BVH) and threading.

