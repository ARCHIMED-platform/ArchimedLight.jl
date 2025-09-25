using Meshes
import Unitful
using Unitful: oneunit

strip_quantity(x::Unitful.Quantity) = Unitful.ustrip(x)
strip_quantity(x::Real) = Float64(x)

struct InterceptionResult
    absorbed::Dict{String,Dict{Band,Float64}}
    incident::Dict{String,Dict{Band,Float64}}
    areas::Dict{String,Float64}
    hit_counts::Dict{String,Int}
end

Base.getindex(res::InterceptionResult, key::String) = res.absorbed[key]
Base.haskey(res::InterceptionResult, key::String) = haskey(res.absorbed, key)
Base.length(res::InterceptionResult) = length(res.absorbed)
Base.iterate(res::InterceptionResult) = iterate(res.absorbed)
Base.iterate(res::InterceptionResult, state) = iterate(res.absorbed, state)
Base.keys(res::InterceptionResult) = keys(res.absorbed)
Base.values(res::InterceptionResult) = values(res.absorbed)
Base.pairs(res::InterceptionResult) = pairs(res.absorbed)

"""
    compute_interception(scene, sky::SkyConfig, cfg::InterceptionConfig)

Computes first-order (no multiple scattering) irradiance and absorptance per component and band.
Returns an `InterceptionResult` collecting incident (`Ri`) and absorbed (`Ra`) fluxes in W m⁻².
"""
function compute_interception(scene::Scene, sky::SkyConfig, cfg::InterceptionConfig)
    # grid over plot bounds in x/y using pixel_size
    lo, up = extrema(scene.bbox)
    # determine grid resolution from bbox extent and pixel_size (unit-agnostic)
    sx = strip_quantity(to(up)[1] - to(lo)[1])
    sy = strip_quantity(to(up)[2] - to(lo)[2])
    nx = max(1, ceil(Int, sx / cfg.pixel_size))
    ny = max(1, ceil(Int, sy / cfg.pixel_size))

    absorbed = Dict{String,Dict{Band,Float64}}()
    incident = Dict{String,Dict{Band,Float64}}()
    areas = Dict{String,Float64}()
    hit_counts = Dict{String,Int}()
    for comp in scene.components
        absorbed[comp.name] = Dict{Band,Float64}(PAR=>0.0, NIR=>0.0, TIR=>0.0)
        incident[comp.name] = Dict{Band,Float64}(PAR=>0.0, NIR=>0.0, TIR=>0.0)
        areas[comp.name] = strip_quantity(area(comp.mesh))
        hit_counts[comp.name] = 0
    end

    samples = nx * ny
    samples = samples == 0 ? 1 : samples

    # for each pixel and sky sector, cast a ray downward from above the scene
    for iy in 1:ny, ix in 1:nx
        x = to(lo)[1] + (ix - 0.5) / nx * (to(up)[1] - to(lo)[1])
        y = to(lo)[2] + (iy - 0.5) / ny * (to(up)[2] - to(lo)[2])
        for sector in sky.sectors
            d = sector.dir
            up3 = to(up)[3]; lo3 = to(lo)[3]
            dz = up3 - lo3
            z0 = up3 + (iszero(dz) ? oneunit(up3) : dz)
            o = Point(x, y, z0)
            if d[3] >= zero(d[3])
                continue  # only downward-going directions supported here
            end
            ray = Ray(o, d)
            best_t = measure(ray)
            best_comp::Union{Nothing,Component} = nothing
            for comp in scene.components
                hit, t = first_hit(ray, comp.mesh)
                if hit && t < best_t
                    best_t = t; best_comp = comp
                end
            end
            if best_comp !== nothing
                props = get_optics(best_comp.name)
                hit_counts[best_comp.name] += 1
                for (band, E) in sky.irradiance
                    flux = sector.weight * E / samples
                    incident[best_comp.name][band] += flux
                    a = absorptance(props, band)
                    absorbed[best_comp.name][band] += flux * a
                end
            end
        end
    end
    return InterceptionResult(absorbed, incident, areas, hit_counts)
end
