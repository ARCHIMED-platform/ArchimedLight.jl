import StaticArrays
import LinearAlgebra: norm, dot, cross
import Dates

function _normalize3(v)
    n = norm(v)
    n == 0.0 ? v : v / n
end

function _normalize3_java(v)
    x = Float32(v[1])
    y = Float32(v[2])
    z = Float32(v[3])
    n2 = (x * x) + (y * y) + (z * z)
    n2 <= 0.0f0 && return StaticArrays.SVector{3,Float64}(Float64(x), Float64(y), Float64(z))
    n = sqrt(n2)
    StaticArrays.SVector{3,Float64}(Float64(x / n), Float64(y / n), Float64(z / n))
end

function _sun_direction(azimuth_deg::Float64, elevation_deg::Float64)
    az = deg2rad(azimuth_deg)
    el = deg2rad(elevation_deg)
    # Java SunPosition gives a ground-to-sky vector.
    sx = cos(el) * sin(az)
    sy = cos(el) * cos(az)
    sz = sin(el)
    _normalize3(StaticArrays.SVector{3,Float64}(sx, sy, sz))
end

function _java_turtle_order(n::Int)
    if n == 1
        return 0
    elseif n == 6
        return 1
    elseif n == 16
        return 2
    elseif n == 46
        return 3
    elseif n == 136
        return 4
    elseif n == 406
        return 5
    end
    return -1
end

function _java_seed_points_up()
    elevation6 = Float32[90.0f0, 26.57f0, 26.57f0, 26.57f0, 26.57f0, 26.57f0]
    azimuth6 = Float32[0.0f0, 180.0f0, 252.0f0, 324.0f0, 36.0f0, 108.0f0]

    pts = StaticArrays.SVector{3,Float64}[]
    for i in eachindex(elevation6)
        elevation = deg2rad(elevation6[i])
        azimuth = deg2rad(azimuth6[i] + 180.0f0)
        c = cos(elevation)
        x = c * sin(azimuth)
        y = c * cos(azimuth)
        z = sin(elevation)
        p = StaticArrays.SVector{3,Float64}(x, y, z)
        push!(pts, p)
        push!(pts, -p)
    end
    pts
end

function _hull_faces_12(points)
    n = length(points)
    n == 12 || error("_hull_faces_12 expects exactly 12 points.")

    faces = NTuple{3,Int}[]
    eps = 1e-9

    for i in 1:(n - 2), j in (i + 1):(n - 1), k in (j + 1):n
        nrm = cross(points[j] - points[i], points[k] - points[i])
        norm(nrm) <= 1e-12 && continue

        pos = false
        neg = false
        for l in 1:n
            (l == i || l == j || l == k) && continue
            d = dot(nrm, points[l] - points[i])
            d > eps && (pos = true)
            d < -eps && (neg = true)
            if pos && neg
                break
            end
        end
        (pos && neg) && continue

        c = (points[i] + points[j] + points[k]) / 3.0
        if dot(nrm, c) < 0.0
            push!(faces, (i, k, j))
        else
            push!(faces, (i, j, k))
        end
    end

    faces
end

function _orient_face_outward(a::Int, b::Int, c::Int, points)
    pa = points[a]
    pb = points[b]
    pc = points[c]
    nrm = cross(pb - pa, pc - pa)
    cen = (pa + pb + pc) / 3.0
    dot(nrm, cen) >= 0.0 ? (a, b, c) : (a, c, b)
end

function _refine_turtle_mesh(points, faces)
    new_points = copy(points)
    new_faces = NTuple{3,Int}[]

    for (a, b, c) in faces
        d = length(new_points) + 1
        bary = _normalize3((points[a] + points[b] + points[c]) / 3.0)
        push!(new_points, bary)
        push!(new_faces, _orient_face_outward(a, b, d, new_points))
        push!(new_faces, _orient_face_outward(b, c, d, new_points))
        push!(new_faces, _orient_face_outward(c, a, d, new_points))
    end

    return new_points, new_faces
end

function _java_turtle_upward_points(order::Int)
    order == 0 && return [StaticArrays.SVector{3,Float64}(0.0, 0.0, 1.0)]

    points = _java_seed_points_up()
    faces = _hull_faces_12(points)

    for _ in 2:order
        points, faces = _refine_turtle_mesh(points, faces)
    end

    [p for p in points if p[3] >= 0.0]
end

function _hemisphere_fibonacci_incoming(n::Int)
    dirs = Vector{StaticArrays.SVector{3,Float64}}(undef, n)
    if n == 1
        dirs[1] = StaticArrays.SVector(0.0, 0.0, 1.0)
        return dirs
    end
    phi = (1.0 + sqrt(5.0)) / 2.0
    for i in 1:n
        z_up = (i - 0.5) / n
        theta = 2.0 * pi * i / phi
        r = sqrt(max(0.0, 1.0 - z_up * z_up))
        x = r * cos(theta)
        y = r * sin(theta)
        dirs[i] = _normalize3(StaticArrays.SVector{3,Float64}(x, y, z_up))
    end
    dirs
end

const _TURTLE_DIR_CACHE = Dict{Int,Vector{StaticArrays.SVector{3,Float64}}}()
const _JAVA_N16_DIR_ORDER = (1, 2, 3, 4, 5, 6, 14, 10, 7, 9, 12, 11, 8, 13, 16, 15)

function _java_turtle_incoming_n6_logged()
    # Java n=6 turtle directions (upward hemisphere), logged from
    # `log-directions.csv` in test-hitcount3.
    up = (
        (3.821371e-15, 4.371139e-8, 1.0),
        (1.5637987e-7, 0.89438856, 0.44729084),
        (0.85061413, 0.27638108, 0.44729084),
        (0.5257086, -0.7235755, 0.4472909),
        (-0.5257085, -0.7235755, 0.44729084),
        (-0.85061413, 0.2763814, 0.4472909),
    )
    [
        _normalize3_java(StaticArrays.SVector{3,Float64}(x, y, z))
        for (x, y, z) in up
    ]
end

function _java_turtle_incoming_n16_logged()
    # Java n=16 turtle directions (upward hemisphere), logged from
    # `log-directions.csv`.
    up = (
        (3.821371e-15, 4.371139e-8, 1.0),
        (1.5637987e-7, 0.89438856, 0.44729084),
        (0.85061413, 0.27638108, 0.44729084),
        (0.5257086, -0.7235755, 0.4472909),
        (-0.5257085, -0.7235755, 0.44729084),
        (-0.85061413, 0.2763814, 0.4472909),
        (-4.0595705e-8, -0.9822395, 0.18763158),
        (2.5001444e-8, -0.6070141, 0.7946911),
        (0.35679406, 0.4910847, 0.794691),
        (0.57730484, -0.18757768, 0.7946909),
        (0.577346, 0.79464835, 0.18763156),
        (-0.5773048, -0.18757758, 0.79469097),
        (-0.35679388, 0.49108484, 0.794691),
        (-0.57734585, 0.7946484, 0.18763156),
        (-0.9341653, -0.30352855, 0.18763156),
        (0.93416524, -0.30352876, 0.18763155),
    )
    [
        _normalize3_java(StaticArrays.SVector{3,Float64}(x, y, z))
        for (x, y, z) in up
    ]
end

function _java_turtle_incoming(n::Int)
    get!(_TURTLE_DIR_CACHE, n) do
        n == 6 && return _java_turtle_incoming_n6_logged()
        # Keep n=16 on the exact Java-logged vectors: one near-axis direction has
        # a sign-sensitive Float32 x component that affects toric border pixels.
        n == 16 && return _java_turtle_incoming_n16_logged()
        order = _java_turtle_order(n)
        order < 0 && error("Unsupported Java turtle sector count: $n. Allowed values: 1, 6, 16, 46, 136, 406.")
        up = _java_turtle_upward_points(order)
        [_normalize3_java(u) for u in up]
    end
end

function _circle_lumen_area(centers_distance::Float64, radius1::Float64, radius2::Float64)
    if centers_distance >= (radius1 + radius2)
        return 0.0
    end

    if radius1 == radius2
        half_dist = centers_distance / 2.0
        h2 = max(radius2 * radius2 - half_dist * half_dist, 0.0)
        h = sqrt(h2)
        alpha = acos(clamp(half_dist / max(radius2, eps(Float64)), -1.0, 1.0))
        return 2.0 * ((alpha * (radius2 * radius2)) - (h * half_dist))
    end

    d2 = centers_distance * centers_distance
    r12 = radius1 * radius1
    r22 = radius2 * radius2

    if centers_distance + radius2 < radius1
        return pi * r22
    end
    if centers_distance + radius1 < radius2
        return pi * r12
    end

    t1 = clamp((d2 + r22 - r12) / (2.0 * centers_distance * radius2), -1.0, 1.0)
    t2 = clamp((d2 + r12 - r22) / (2.0 * centers_distance * radius1), -1.0, 1.0)
    k = (-centers_distance + radius2 + radius1) *
        (centers_distance + radius2 - radius1) *
        (centers_distance - radius2 + radius1) *
        (centers_distance + radius2 + radius1)
    term = sqrt(max(k, 0.0)) / 2.0

    (r22 * acos(t1)) + (r12 * acos(t2)) - term
end

function _soc_fraction_hourly(diffuse::Float64, global_flux::Float64, elevation_rad::Float64)
    if elevation_rad <= 0.0
        return 0.5
    end

    global_flux <= 0.0 && return 0.0
    ratio = diffuse / global_flux

    sin_sun = sin(elevation_rad)
    r_clear = 0.847 - (1.61 * sin_sun) + (1.04 * sin_sun * sin_sun)
    return max((ratio - r_clear) / (1.0 - r_clear), 0.0)
end

function _brightness_norm_soc(elevation::Float64)
    (1.0 + 2.0 * sin(elevation)) / 3.0
end

function _brightness_norm_clear(dir_up, sun_up)
    cos_sun_zen = cos(acos(clamp(sun_up[3], -1.0, 1.0)))
    sin_elevation = sin(asin(clamp(dir_up[3], -1.0, 1.0)))
    if sin_elevation <= 0.0
        return 0.0
    end

    angle = acos(clamp(dot(sun_up, dir_up), -1.0, 1.0))
    cos_angle = cos(angle)

    brightness = 0.91 + 10.0 * exp(-3.0 * angle)
    brightness += 0.45 * (cos_angle * cos_angle)
    brightness *= 1.0 - exp(-0.32 / sin_elevation)

    denom = 0.91 + (10.0 * exp(-3.0 * acos(clamp(sun_up[3], -1.0, 1.0)))) + (0.45 * cos_sun_zen * cos_sun_zen)
    brightness /= max(denom, eps(Float64))
    brightness /= 0.27385
    brightness
end

function _diffuse_weights_java_like(sky::SkyState, turtle::TurtleGrid, sky_ids::Vector{Int})
    n = length(sky_ids)
    n == 0 && return Float64[]

    sun_up = _sun_direction(sky.sun_azimuth_deg, sky.sun_elevation_deg)
    raw = zeros(Float64, n)
    total = 0.0
    for (k, i) in enumerate(sky_ids)
        sec_up = turtle.sectors[i].direction
        elev = asin(clamp(sec_up[3], -1.0, 1.0))
        coeff_soc = _soc_fraction_hourly(sky.diffuse_fraction, 1.0, elev)
        coeff_clear = 1.0 - coeff_soc
        b_soc = _brightness_norm_soc(elev)
        b_clear = _brightness_norm_clear(sec_up, sun_up)
        b = coeff_soc * b_soc + coeff_clear * b_clear

        # horizontal-plane conversion as in Java (multiply by cos(zenith))
        coeff = max(turtle.sectors[i].direction[3], 0.0)
        raw[k] = b * coeff
        total += raw[k]
    end

    if total > 0.0
        return raw ./ total
    end

    # Defensive fallback.
    fill(1.0 / n, n)
end

"""
    build_turtle(cfg, sky)::TurtleGrid

Build the directional sky discretization (turtle sectors), optionally adding an explicit sun
sector when `cfg.general["all_in_turtle"] == false`.
"""
function build_turtle(cfg::LightConfig, sky::SkyState)
    n = max(get(cfg.general, "sky_sectors", 46), 1)
    dirs =
        if _java_turtle_order(n) >= 0
            _java_turtle_incoming(n)
        else
            _hemisphere_fibonacci_incoming(n)
        end
    sectors = TurtleSector[]
    w = 1.0 / n
    for i in eachindex(dirs)
        push!(sectors, TurtleSector(i, dirs[i], w, :sky))
    end

    if !get(cfg.general, "all_in_turtle", false) && sky.sun_elevation_deg > 0.0
        push!(sectors, TurtleSector(length(sectors) + 1, _sun_direction(sky.sun_azimuth_deg, sky.sun_elevation_deg), 0.0, :sun))
    end
    TurtleGrid(sectors)
end

function _directional_flux_substeps(meteo_row, sky::SkyState, cfg::LightConfig)
    latitude_deg = _row_value(meteo_row, [:latitude, :lat], 0.0)
    latitude_rad = deg2rad(latitude_deg)
    date = _row_date(meteo_row)
    doy = Dates.dayofyear(date)
    start_h, end_h = _row_step_hours(meteo_row)

    ri_sw_raw = _row_value(meteo_row, [:RI_SW_f, :Rg, :rg, :sw_global, :global], NaN)
    ri_par_raw = _row_value(meteo_row, [:RI_PAR_f, :PAR, :par], NaN)
    ri_nir_raw = _row_value(meteo_row, [:RI_NIR_f, :NIR, :nir], NaN)
    global_from_input = !isnan(ri_sw_raw) || !isnan(ri_par_raw) || !isnan(ri_nir_raw)

    clearness_raw = _row_value(meteo_row, [:clearness, :Kt], NaN)
    clearness_provided = !isnan(clearness_raw)
    clearness =
        if clearness_provided
            clearness_raw
        else
            _clearness_from_global_wm2(sky.ri_sw_f, latitude_rad, doy, start_h, end_h)
        end

    substeps = _radiation_substeps(
        date,
        start_h,
        end_h,
        latitude_rad,
        _cfg_radiation_timestep_hours(cfg),
        sky.ri_sw_f,
        clearness,
        global_from_input,
        clearness_provided,
    )
    return substeps, start_h, end_h
end

"""
    compute_directional_fluxes(meteo_row, sky, turtle, cfg)::DirectionalFluxes

Java-parity directional flux integration using meteo substeps:
directional fluxes are computed at each substep sun position then averaged over the full meteo step.
"""
function compute_directional_fluxes(meteo_row, sky::SkyState, turtle::TurtleGrid, cfg::LightConfig)
    substeps, start_h, end_h = _directional_flux_substeps(meteo_row, sky, cfg)
    step_duration_h = end_h - start_h
    (isempty(substeps) || step_duration_h <= 0.0) && return compute_directional_fluxes(sky, turtle, cfg)

    n = length(turtle.sectors)
    par_acc = zeros(Float64, n)
    nir_acc = zeros(Float64, n)
    par_ratio = sky.ri_sw_f > 0.0 ? sky.ri_par_f / sky.ri_sw_f : 0.0
    nir_ratio = sky.ri_sw_f > 0.0 ? sky.ri_nir_f / sky.ri_sw_f : 0.0

    for ss in substeps
        total_w = ss.diffuse_w + ss.direct_w
        total_w > 0.0 || continue
        direct_fraction_sub = clamp(ss.direct_w / total_w, 0.0, 1.0)
        sky_sub = SkyState(
            ss.sun_azimuth_deg,
            ss.sun_elevation_deg,
            max(ss.global_w * par_ratio, 0.0),
            max(ss.global_w * nir_ratio, 0.0),
            direct_fraction_sub,
            1.0 - direct_fraction_sub,
        )
        flux_sub = compute_directional_fluxes(sky_sub, turtle, cfg)
        @inbounds for i in 1:n
            par_acc[i] += flux_sub.par[i] * ss.duration
            nir_acc[i] += flux_sub.nir[i] * ss.duration
        end
    end

    par = par_acc ./ step_duration_h
    nir = nir_acc ./ step_duration_h

    # Preserve step-level band totals exactly after substep averaging.
    sum_par = sum(par)
    sum_nir = sum(nir)
    if sky.ri_par_f <= 0.0 || sum_par <= 0.0
        fill!(par, 0.0)
    else
        par .*= sky.ri_par_f / sum_par
    end
    if sky.ri_nir_f <= 0.0 || sum_nir <= 0.0
        fill!(nir, 0.0)
    else
        nir .*= sky.ri_nir_f / sum_nir
    end

    DirectionalFluxes([s.id for s in turtle.sectors], par, nir)
end

"""
    compute_directional_fluxes(sky, turtle, cfg)::DirectionalFluxes

Project sky-level PAR/NIR irradiance to each turtle sector using Java-compatible diffuse and
direct distribution rules.
"""
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
    if !isempty(sky_ids)
        wdiff = _diffuse_weights_java_like(sky, turtle, sky_ids)
        for (k, i) in enumerate(sky_ids)
            par[i] += par_diff * wdiff[k]
            nir[i] += nir_diff * wdiff[k]
        end
    end

    if get(cfg.general, "all_in_turtle", false)
        # Java Turtle.directInTurtle parity:
        # distribute direct irradiance by angular overlap between a sun halo and each sector.
        dir_count = length(sky_ids)
        if dir_count > 0
            sun_up = _sun_direction(sky.sun_azimuth_deg, sky.sun_elevation_deg)
            sector_radius = acos((dir_count - 1) / max(dir_count, 1))
            sun_halo_radius = sector_radius / 2.0

            raw = zeros(Float64, dir_count)
            horiz = zeros(Float64, dir_count)
            total_horiz = 0.0

            for (k, i) in enumerate(sky_ids)
                sec_up = turtle.sectors[i].direction
                angle = acos(clamp(dot(sun_up, sec_up), -1.0, 1.0))
                w = _circle_lumen_area(angle, sector_radius, sun_halo_radius)
                raw[k] = w

                # Java TurtleFluxes converts direct-in-sector values to horizontal
                # irradiance using cos(zenith), then normalizes to total direct flux.
                coeff = max(turtle.sectors[i].direction[3], 0.0)
                horiz[k] = w * coeff
                total_horiz += horiz[k]
            end

            if total_horiz > 0.0
                par_scale = par_dir / total_horiz
                nir_scale = nir_dir / total_horiz
                for (k, i) in enumerate(sky_ids)
                    par[i] += horiz[k] * par_scale
                    nir[i] += horiz[k] * nir_scale
                end
            elseif sum(raw) > 0.0
                # Defensive fallback if all sectors are numerically near-horizon.
                wsum = sum(raw)
                for (k, i) in enumerate(sky_ids)
                    w = raw[k] / wsum
                    par[i] += par_dir * w
                    nir[i] += nir_dir * w
                end
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
