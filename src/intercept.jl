using Meshes
import Unitful
using Unitful: oneunit
using LinearAlgebra: dot, cross

strip_quantity(x::Unitful.Quantity) = Unitful.ustrip(x)
strip_quantity(x::Real) = Float64(x)

# convert mesh/point Vec/Point (possibly Unitful) to Float64 triple (meters)
point_to_tuple(p) = (
    strip_quantity(to(p)[1]),
    strip_quantity(to(p)[2]),
    strip_quantity(to(p)[3])
)

# convert a direction Vec (possibly Unitful) to Float64 triple
vec_to_tuple(v) = (
    strip_quantity(v[1]),
    strip_quantity(v[2]),
    strip_quantity(v[3])
)

# Möller–Trumbore triangle intersection returning t (distance) or nothing
function tri_intersect(orig::NTuple{3,Float64}, dir::NTuple{3,Float64},
    v0::NTuple{3,Float64}, v1::NTuple{3,Float64}, v2::NTuple{3,Float64})
    eps = 1e-9
    edge1 = (v1[1] - v0[1], v1[2] - v0[2], v1[3] - v0[3])
    edge2 = (v2[1] - v0[1], v2[2] - v0[2], v2[3] - v0[3])
    pvec = (dir[2] * edge2[3] - dir[3] * edge2[2], dir[3] * edge2[1] - dir[1] * edge2[3], dir[1] * edge2[2] - dir[2] * edge2[1])
    det = edge1[1] * pvec[1] + edge1[2] * pvec[2] + edge1[3] * pvec[3]
    if abs(det) < eps
        return nothing
    end
    invdet = 1.0 / det
    tvec = (orig[1] - v0[1], orig[2] - v0[2], orig[3] - v0[3])
    u = (tvec[1] * pvec[1] + tvec[2] * pvec[2] + tvec[3] * pvec[3]) * invdet
    if u < 0.0 || u > 1.0
        return nothing
    end
    qvec = (tvec[2] * edge1[3] - tvec[3] * edge1[2], tvec[3] * edge1[1] - tvec[1] * edge1[3], tvec[1] * edge1[2] - tvec[2] * edge1[1])
    v = (dir[1] * qvec[1] + dir[2] * qvec[2] + dir[3] * qvec[3]) * invdet
    if v < 0.0 || u + v > 1.0
        return nothing
    end
    t = (edge2[1] * qvec[1] + edge2[2] * qvec[2] + edge2[3] * qvec[3]) * invdet
    if t > eps
        return t
    end
    return nothing
end

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
        absorbed[comp.name] = Dict{Band,Float64}(PAR => 0.0, NIR => 0.0, TIR => 0.0)
        incident[comp.name] = Dict{Band,Float64}(PAR => 0.0, NIR => 0.0, TIR => 0.0)
        areas[comp.name] = strip_quantity(area(comp.mesh))
        hit_counts[comp.name] = 0
    end

    samples = nx * ny
    samples = samples == 0 ? 1 : samples

    # for each pixel and sky sector, cast a ray downward from above the scene
    for iy in 1:ny, ix in 1:nx
        x = strip_quantity(to(lo)[1]) + (ix - 0.5) / nx * (strip_quantity(to(up)[1]) - strip_quantity(to(lo)[1]))
        y = strip_quantity(to(lo)[2]) + (iy - 0.5) / ny * (strip_quantity(to(up)[2]) - strip_quantity(to(lo)[2]))
        for sector in sky.sectors
            d_raw = sector.dir
            d = vec_to_tuple(d_raw)
            up3 = strip_quantity(to(up)[3])
            lo3 = strip_quantity(to(lo)[3])
            dz = up3 - lo3
            z0 = up3 + (iszero(dz) ? 1.0 : dz)
            o = (x, y, z0)
            if d[3] >= 0.0
                continue  # only downward-going directions supported here
            end

            # collect all intersections (multi-hit) as (t, comp)
            hits = Tuple{Float64,Component}[]
            for comp in scene.components
                # iterate triangles in mesh with unit-stripped vertices
                v = vertices(comp.mesh)
                for tri_conn in topology(comp.mesh)
                    idx = indices(tri_conn)
                    v0 = point_to_tuple(v[idx[1]])
                    v1 = point_to_tuple(v[idx[2]])
                    v2 = point_to_tuple(v[idx[3]])
                    t = tri_intersect(o, d, v0, v1, v2)
                    if t !== nothing
                        push!(hits, (t, comp))
                    end
                end
            end
            if isempty(hits)
                continue
            end
            sort!(hits, by=x -> x[1])

            # propagate beam through successive hits
            for (band, E) in sky.irradiance
                flux = sector.weight * E / samples
                remaining = flux
                # store per-hit incident and scattered parts for optional scattering post-process
                local_hit_incidents = Dict{String,Float64}()
                local_hit_scattered = Dict{String,Float64}()
                for (t, comp) in hits
                    props = get_optics(comp.name)
                    # incident on this component (Ri)
                    Ri = remaining
                    incident[comp.name][band] += Ri
                    local_hit_incidents[comp.name] = get(local_hit_incidents, comp.name, 0.0) + Ri

                    # straight-through transparency (fraction of Ri that passes without interaction)
                    straight = clamp(props.transparency, 0.0, 1.0) * Ri
                    intercepted = Ri - straight

                    # Determine optics partitioning.
                    # If the optics object carries an explicit scattering σ (Archimed-style),
                    # use σ and set ρ = τ = σ/2 and absorptance α = 1 - σ.
                    # Otherwise, respect the explicit rho/tau supplied in the OpticalProps.
                    σ_cfg = get(props.scattering, band, 0.0)
                    if σ_cfg > 0.0
                        σ = clamp(σ_cfg, 0.0, 1.0)
                        ρ = σ / 2.0
                        τ = σ / 2.0
                        α = max(0.0, 1.0 - σ)
                    else
                        # use explicit rho/tau/absorptance from the optics properties
                        ρ = get(props.rho, band, 0.0)
                        τ = get(props.tau, band, 0.0)
                        α = max(0.0, 1.0 - ρ - τ)
                        σ = clamp(ρ + τ, 0.0, 1.0)
                    end

                    # absorbed (Ra) from intercepted portion
                    Ra = α * intercepted
                    absorbed[comp.name][band] += Ra

                    # scattered portion (to be redistributed): sum of reflected + transmitted parts of intercepted energy
                    scattered = intercepted * σ
                    local_hit_scattered[comp.name] = get(local_hit_scattered, comp.name, 0.0) + scattered

                    # transmitted forward: straight-through + τ * intercepted
                    remaining = straight + τ * intercepted

                    hit_counts[comp.name] += 1
                    # if negligible remaining, stop
                    if remaining < 1e-12
                        break
                    end
                end

                # single-scatter redistribution: if scattering enabled, redistribute
                if cfg.scattering
                    # compute total scattered energy from intercepted hits
                    total_scattered = 0.0
                    for (cname, scatamount) in local_hit_scattered
                        total_scattered += scatamount
                    end
                    if total_scattered > 0
                        # redistribute uniformly weighted by component area
                        total_area = sum(values(areas))
                        for comp in scene.components
                            frac = areas[comp.name] / total_area
                            # add redistributed diffuse incident
                            incident[comp.name][band] += total_scattered * frac
                            # assume diffuse fraction is partially absorbed according to each component's absorptance
                            a2 = absorptance(get_optics(comp.name), band)
                            absorbed[comp.name][band] += total_scattered * frac * a2
                        end
                    end
                end
            end
        end
    end
    return InterceptionResult(absorbed, incident, areas, hit_counts)
end
