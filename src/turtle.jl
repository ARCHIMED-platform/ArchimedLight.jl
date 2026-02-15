import StaticArrays
import LinearAlgebra: norm, dot

function _normalize3(v)
    n = norm(v)
    n == 0.0 ? v : v / n
end

function _sun_direction(azimuth_deg::Float64, elevation_deg::Float64)
    az = deg2rad(azimuth_deg)
    el = deg2rad(elevation_deg)
    # Incoming direction from source to scene.
    sx = cos(el) * sin(az)
    sy = cos(el) * cos(az)
    sz = -sin(el)
    _normalize3(StaticArrays.SVector{3,Float64}(sx, sy, sz))
end

function _hemisphere_fibonacci_incoming(n::Int)
    dirs = Vector{StaticArrays.SVector{3,Float64}}(undef, n)
    if n == 1
        dirs[1] = StaticArrays.SVector(0.0, 0.0, -1.0)
        return dirs
    end
    phi = (1.0 + sqrt(5.0)) / 2.0
    for i in 1:n
        z_up = (i - 0.5) / n
        theta = 2.0 * pi * i / phi
        r = sqrt(max(0.0, 1.0 - z_up * z_up))
        x = r * cos(theta)
        y = r * sin(theta)
        # Convert upward hemisphere direction into incoming/downward direction.
        dirs[i] = _normalize3(StaticArrays.SVector{3,Float64}(x, y, -z_up))
    end
    dirs
end

function build_turtle(cfg::LightConfig, sky::SkyState)
    n = max(cfg.turtle_sectors, 1)
    dirs = _hemisphere_fibonacci_incoming(n)
    sectors = TurtleSector[]
    w = 1.0 / n
    for i in eachindex(dirs)
        push!(sectors, TurtleSector(i, dirs[i], w, :sky))
    end

    if !cfg.all_in_turtle && sky.sun_elevation_deg > 0.0
        push!(sectors, TurtleSector(length(sectors) + 1, _sun_direction(sky.sun_azimuth_deg, sky.sun_elevation_deg), 0.0, :sun))
    end
    TurtleGrid(sectors)
end

function compute_directional_fluxes(sky::SkyState, turtle::TurtleGrid, cfg::LightConfig)
    n = length(turtle.sectors)
    par = zeros(Float64, n)
    nir = zeros(Float64, n)

    sky_ids = findall(s -> s.source == :sky, turtle.sectors)
    sun_ids = findall(s -> s.source == :sun, turtle.sectors)

    par_diff = sky.ri_par_f * sky.diffuse_fraction
    nir_diff = sky.ri_nir_f * sky.diffuse_fraction
    par_dir = sky.ri_par_f * sky.direct_fraction
    nir_dir = sky.ri_nir_f * sky.direct_fraction

    sky_wsum = sum(turtle.sectors[i].weight for i in sky_ids)
    if sky_wsum > 0.0
        for i in sky_ids
            w = turtle.sectors[i].weight / sky_wsum
            par[i] += par_diff * w
            nir[i] += nir_diff * w
        end
    end

    if cfg.all_in_turtle
        sun_dir = _sun_direction(sky.sun_azimuth_deg, sky.sun_elevation_deg)
        scores = [max(dot(-sun_dir, -turtle.sectors[i].direction), 0.0) for i in sky_ids]
        ssum = sum(scores)
        if ssum > 0.0
            for (k, i) in enumerate(sky_ids)
                w = scores[k] / ssum
                par[i] += par_dir * w
                nir[i] += nir_dir * w
            end
        end
    else
        if !isempty(sun_ids)
            i = first(sun_ids)
            par[i] += par_dir
            nir[i] += nir_dir
        else
            # Fallback if no explicit sun sector is available.
            for i in sky_ids
                w = turtle.sectors[i].weight / max(sky_wsum, eps(Float64))
                par[i] += par_dir * w
                nir[i] += nir_dir * w
            end
        end
    end

    DirectionalFluxes([s.id for s in turtle.sectors], par, nir)
end

