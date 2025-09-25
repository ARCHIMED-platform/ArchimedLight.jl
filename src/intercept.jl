using Meshes
import Unitful
using Unitful: oneunit
include("geometry.jl")
include("scene.jl")
include("optics.jl")
include("sky.jl")

"""
    compute_interception(scene, sky::SkyConfig, cfg::InterceptionConfig)

Computes first-order (no multiple scattering) absorbed irradiance per component and band.
Returns Dict(component_name => Dict{Band,Float64} absorbed_Wm2).
"""
function compute_interception(scene::Scene, sky::SkyConfig, cfg::InterceptionConfig)
    # grid over plot bounds in x/y using pixel_size
    lo, up = extrema(scene.bbox)
    # determine grid resolution from bbox extent and pixel_size (unit-agnostic)
    sx = Unitful.ustrip(to(up)[1] - to(lo)[1])
    sy = Unitful.ustrip(to(up)[2] - to(lo)[2])
    nx = max(1, ceil(Int, sx / cfg.pixel_size))
    ny = max(1, ceil(Int, sy / cfg.pixel_size))

    # accumulator per component per band
    absorbed = Dict{String,Dict{Band,Float64}}()
    for comp in scene.components
        absorbed[comp.name] = Dict{Band,Float64}(PAR=>0.0, NIR=>0.0, TIR=>0.0)
    end

    # for each pixel and sky sector, cast a ray downward from above the scene
    for iy in 1:ny, ix in 1:nx
        x = to(lo)[1] + (ix - 0.5) / nx * (to(up)[1] - to(lo)[1])
        y = to(lo)[2] + (iy - 0.5) / ny * (to(up)[2] - to(lo)[2])
        for sector in sky.sectors
            d = sector.dir
            # place origin sufficiently above scene along opposite of direction
            up3 = to(up)[3]; lo3 = to(lo)[3]
            dz = up3 - lo3
            z0 = up3 + (iszero(dz) ? oneunit(up3) : dz)
            o = Point(x, y, z0)
            if d[3] >= zero(d[3])
                continue  # only downward-going directions supported here
            end
            ray = Ray(o, d)
            # find first hit among all components
            best_t = measure(ray)
            best_comp::Union{Nothing,Component} = nothing
            for comp in scene.components
                hit, t = first_hit(ray, comp.mesh)
                if hit && t < best_t
                    best_t = t; best_comp = comp
                end
            end
            if best_comp !== nothing
                # accumulate absorbed power proportional to sector weight
                props = get_optics(best_comp.name)
                for (band, E) in sky.irradiance
                    a = absorptance(props, band)
                    absorbed[best_comp.name][band] += sector.weight * E * a
                end
            end
        end
    end
    return absorbed
end
