using Meshes
using LinearAlgebra: norm, dot
using Unitful

function stripq(x)
    if x isa Unitful.Quantity
        return Unitful.ustrip(x)
    elseif x isa AbstractArray
        return [stripq(el) for el in x]
    elseif x isa Number
        return Float64(x)
    else
        return x
    end
end

struct Sector
    dir::Vec3       # unit direction (toward scene)
    weight::Float64 # relative weight, sum to 1 over sky
end

mutable struct SkyConfig
    sectors::Vector{Sector}
    irradiance::Dict{Band,Float64}            # total band irradiance (W m^-2)
    sector_fluxes::Vector{Dict{Band,Float64}} # per-sector irradiance (W m^-2)
end

SkyConfig(sectors::Vector{Sector}, irradiance::Dict{Band,Float64}) =
    SkyConfig(sectors, irradiance, Dict{Band,Float64}[])

# --- Port of Archimed's DeJong diffuse fraction (hourly) ---
function dejong_kd_hourly(kt::Float64, sun_elevation::Float64)
    if kt <= 0.22
        return 1.0
    elseif kt <= 0.35
        return 1.0 - (6.4 * (kt - 0.22) * (kt - 0.22))
    else
        R = 0.847 - (1.61 * sin(sun_elevation)) + (1.04 * sin(sun_elevation)^2)
        K = (1.47 - R) / 1.66
        if kt <= K
            return 1.47 - (1.66 * kt)
        end
        return R
    end
end

# Brightness normalization functions ported (brightnessNormSOC / brightnessNormClear / brightnessNorm)
function brightness_norm_soc(elevation::Float64)
    return (1.0 + 2.0 * sin(elevation)) / 3.0
end

function brightness_norm_clear(direction::Vec, sun_dir::Vec)
    # direction and sun_dir are unit vectors; compute angle between them
    # strip units from vectors
    sd = stripq(sun_dir)
    dir = stripq(direction)
    sun_elev = asin(clamp(sd[3], -1.0, 1.0))
    cos_sunzen = cos(sun_elev)
    elev = asin(clamp(dir[3], -1.0, 1.0))
    if elev < 0.01
        return 0.0
    end
    sin_elev = sin(elev)
    # compute angle between vectors
    dotp = clamp(dot(dir, sd), -1.0, 1.0)
    angle = acos(dotp)
    cos_angle = cos(angle)

    brightness = 0.91 + 10.0 * exp(-3.0 * angle)
    brightness += 0.45 * (cos_angle * cos_angle)
    brightness *= 1.0 - exp(-0.32 / sin_elev)
    denom = 0.91 + (10.0 * exp(-3.0 * sun_elev)) + (0.45 * cos_sunzen * cos_sunzen)
    brightness /= denom
    brightness /= 0.27385
    return brightness
end

function brightness_norm(diffuse::Float64, global_flux::Float64, direction::Vec, sun_dir::Vec)
    # compute SOC fraction using DeJong relation approximate (use sun elevation)
    sd = stripq(sun_dir)
    sun_elev = asin(clamp(sd[3], -1.0, 1.0))
    coeffSOC = 0.5 # fallback
    if global_flux > 0
        # approximate diffuse/global fraction
        frac = diffuse / global_flux
        R = 0.847 - (1.61 * sin(sun_elev)) + (1.04 * sin(sun_elev)^2)
        coeffSOC = clamp((frac - R) / (1 - R), 0.0, 1.0)
    end
    coeffClear = 1.0 - coeffSOC
    return coeffSOC * brightness_norm_soc(asin(clamp(stripq(direction[3]), -1.0, 1.0))) + coeffClear * brightness_norm_clear(direction, sun_dir)
end

# Distribute direct (sun) irradiance into turtle sectors similar to Turtle.directInTurtle
function direct_in_turtle(direct::Float64, sun_dir::Vec, sectors::Vector{Sector}, all_in_turtle::Bool)
    n = length(sectors)
    weights = zeros(Float64, n)
    # if sun is below horizon, nothing
    if stripq(sun_dir[3]) >= 0.0
        if !all_in_turtle
            # last sector reserved for sun direction in Archimed Java
            weights[end] = direct
            return weights
        end

        # approximate sector radius based on number of sectors
        directions_sector_radius = acos(clamp((n - 1) / Float64(n), -1.0, 1.0))
        sun_halo_radius = directions_sector_radius / 2.0

        # compute exact spherical-cap intersection areas between sector caps and sun halo
        total = 0.0
        # sector angular radius (approx) and sun halo angular radius
        sector_radius = directions_sector_radius
        halo_radius = sun_halo_radius

        # Helper: area of a spherical cap of angular radius r on unit sphere
        cap_area(r) = 2π * (1.0 - cos(r))

        # intersection area of two spherical caps with angular radii a,b whose centers are separated by angle g
        function spherical_cap_intersection_area(a::Float64, b::Float64, g::Float64)
            # trivial cases
            if g >= a + b
                return 0.0
            end
            if g <= abs(a - b)
                # smaller cap inside larger
                return cap_area(min(a, b))
            end

            # integrand for θ in [0, a]: at colatitude θ around cap1 center, ring of points has azimuthal extent 2*acos(q(θ)) where
            # q(θ) = (cos b - cos θ cos g) / (sin θ sin g)
            # intersection area = ∫_{θ=0}^{a} 2 * acos(clamp(q, -1,1)) * sin θ dθ

            # handle potential numerical issues when sin g ~ 0
            if sin(g) == 0.0
                # centers aligned or opposite handled above; fallback to smaller cap
                return cap_area(min(a, b))
            end

            # adaptive Simpson integration
            f(θ) = begin
                sθ = sin(θ)
                denom = sθ * sin(g)
                if denom == 0.0
                    q = 1.0
                else
                    q = (cos(b) - cos(θ) * cos(g)) / denom
                end
                if q <= -1.0
                    phi = π
                elseif q >= 1.0
                    phi = 0.0
                else
                    phi = acos(q)
                end
                return 2.0 * phi * sθ
            end

            # Simpson's rule recursive
            function simpson(f, aθ, bθ)
                c = (aθ + bθ) / 2.0
                return (f(aθ) + 4.0 * f(c) + f(bθ)) * (bθ - aθ) / 6.0
            end

            function adaptive_simpson(f, aθ, bθ, eps, maxrec, whole)
                c = (aθ + bθ) / 2.0
                left = simpson(f, aθ, c)
                right = simpson(f, c, bθ)
                if maxrec <= 0 || abs(left + right - whole) < 15.0 * eps
                    return left + right + (left + right - whole) / 15.0
                end
                return adaptive_simpson(f, aθ, c, eps / 2.0, maxrec - 1, left) + adaptive_simpson(f, c, bθ, eps / 2.0, maxrec - 1, right)
            end

            aθ = 0.0
            bθ = a
            whole = simpson(f, aθ, bθ)
            area = adaptive_simpson(f, aθ, bθ, 1e-8, 20, whole)
            return area
        end

        for (i, s) in enumerate(sectors)
            # compute center separation angle
            d1 = stripq(s.dir)
            if d1 isa AbstractVector
                d1v = Vec(d1...)
            else
                d1v = d1
            end
            d2 = stripq(sun_dir)
            if d2 isa AbstractVector
                d2v = Vec(d2...)
            else
                d2v = d2
            end
            g = acos(clamp(dot(d1v, d2v), -1.0, 1.0))
            a = sector_radius
            b = halo_radius
            inter = spherical_cap_intersection_area(a, b, g)
            weights[i] = max(0.0, inter)
            total += weights[i]
        end
        if total > 0
            # distribute direct proportional to area intersection
            weights .*= (direct / total)
        end
    end
    return weights
end

# Populate per-sector fluxes for a SkyConfig given global/direct/diffuse partition and sun direction
function populate_sector_fluxes!(sky::SkyConfig, band::Band, global_flux::Float64, direct::Float64, diffuse::Float64, sun_dir::Vec, all_in_turtle::Bool)
    n = length(sky.sectors)
    # initialize per-sector flux dicts only if not present or if sector count changed
    if isempty(sky.sector_fluxes) || length(sky.sector_fluxes) != n
        sky.sector_fluxes = [Dict{Band,Float64}() for _ in 1:n]
    end

    # compute per-sector diffuse weighting using brightness_norm
    diffuse_weights = zeros(Float64, n)
    total_diff = 0.0
    sd = stripq(sun_dir)
    # convert sd to Vec if it's an AbstractVector
    if sd isa AbstractVector
        sdv = Vec(sd...)
    else
        sdv = sd
    end
    for (i, s) in enumerate(sky.sectors)
        # ensure direction passed is unitless Vec
        dir = stripq(s.dir)
        if dir isa AbstractVector
            dirv = Vec(dir...)
        else
            dirv = dir
        end
        b = brightness_norm(diffuse, global_flux, dirv, sdv)
        diffuse_weights[i] = max(0.0, b)
        total_diff += diffuse_weights[i]
    end
    if total_diff > 0
        diffuse_weights .*= (diffuse / total_diff)
    end

    # Debug: when computing PAR, log per-sector diffuse weights for analysis
    if band == PAR
        @info "PAR diffuse weights" total_diff = total_diff diffuse = diffuse global_flux = global_flux
        for (i, w) in enumerate(diffuse_weights)
            @debug "PAR diffuse sector" idx = i weight = w
        end
    end

    # direct allocation
    # direct allocation: pass unitless sun_dir and converted sectors
    # create a temporary sectors list with unitless dirs
    tmp_sectors = Sector[]
    for s in sky.sectors
        d = stripq(s.dir)
        if d isa AbstractVector
            dv = Vec(d...)
        else
            dv = d
        end
        push!(tmp_sectors, Sector(dv, s.weight))
    end
    direct_weights = direct_in_turtle(direct, sdv, tmp_sectors, all_in_turtle)

    if band == PAR
        @info "PAR direct weights sum" sum_direct = sum(direct_weights) direct = direct
        for (i, w) in enumerate(direct_weights)
            @debug "PAR direct sector" idx = i weight = w
        end
    end
    for i in 1:n
        # accumulate band flux into the per-sector dict (preserve other bands)
        prev = get(sky.sector_fluxes[i], band, 0.0)
        val = prev + diffuse_weights[i] + direct_weights[i]
        sky.sector_fluxes[i][band] = val
        if band == PAR
            @debug "PAR sector final flux" idx = i val = val
        end
    end
end

function uniform_hemisphere(n::Int)
    # naive stratified sampling over phi, theta for downward hemisphere (z<0)
    m = max(1, round(Int, sqrt(n)))
    sectors = Sector[]
    total = 0.0
    for i in 1:m, j in 1:m
        u = (i - 0.5) / m
        v = (j - 0.5) / m
        θ = acos(1 - u)
        ϕ = 2π * v
        d = Vec(sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), -abs(cos(θ)))
        if d[3] < zero(d[3])
            push!(sectors, Sector(Vec(d / norm(d)), 1.0))
            total += 1.0
        end
    end
    # normalize weights (Sector is immutable)
    return [Sector(s.dir, s.weight / total) for s in sectors]
end

function SkyConfig(; count::Int=16, PAR_Wm2::Float64=400.0, NIR_Wm2::Float64=400.0)
    sectors = uniform_hemisphere(count)
    irr = Dict{Band,Float64}(PAR => PAR_Wm2, NIR => NIR_Wm2, TIR => 0.0)
    per_sector = Dict{Band,Float64}[
        Dict{Band,Float64}(PAR => PAR_Wm2 / count, NIR => NIR_Wm2 / count, TIR => 0.0)
        for _ in 1:length(sectors)
    ]
    return SkyConfig(sectors, irr, per_sector)
end

struct InterceptionConfig
    pixel_size::Float64
    scattering::Bool
    all_in_turtle::Bool
    radiation_timestep::Float64
end

InterceptionConfig(; pixel_size::Float64=0.1, scattering::Bool=false, all_in_turtle::Bool=false, radiation_timestep::Float64=0.0) = InterceptionConfig(pixel_size, scattering, all_in_turtle, radiation_timestep)
