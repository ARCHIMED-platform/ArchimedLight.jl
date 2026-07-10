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

function _as_bool_local_turtle(x, default::Bool)
    x === nothing && return default
    x isa Bool && return x
    if x isa Number
        return x != 0
    end
    s = lowercase(strip(string(x)))
    s in ("", "nothing", "missing", "null") && return default
    s in ("true", "t", "yes", "y", "on", "1") && return true
    s in ("false", "f", "no", "n", "off", "0") && return false
    return default
end

function _use_java_logged_turtle_dirs(options::LightOptions)
    options.java_logged_turtle_dirs
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

function _scene_local_direction(direction, options::LightOptions)
    angle = deg2rad(options.scene_rotation_deg)
    iszero(angle) && return direction
    c = cos(angle)
    s = sin(angle)
    x, y, z = direction
    StaticArrays.SVector{3,Float64}(c * x + s * y, -s * x + c * y, z)
end

_scene_local_sun_direction(sky::SkyState, options::LightOptions) =
    _scene_local_direction(_sun_direction(sky.sun_azimuth_deg, sky.sun_elevation_deg), options)

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

function _java_seed_points_up_f32()
    _java_seed_points_up()
end

function _hull_faces_12(points)
    n = length(points)
    n == 12 || error("_hull_faces_12 expects exactly 12 points.")

    faces = NTuple{3,Int}[]
    eps = 1e-9

    for i in 1:(n-2), j in (i+1):(n-1), k in (j+1):n
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

function _java_turtle_upward_points_f32(order::Int)
    order == 0 && return [StaticArrays.SVector{3,Float64}(0.0, 0.0, 1.0)]

    points = _java_seed_points_up_f32()
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
const _TURTLE_DIR_COMPAT_CACHE = Dict{Int,Vector{StaticArrays.SVector{3,Float64}}}()
const _JAVA_N16_DIR_ORDER = (1, 2, 3, 4, 5, 6, 14, 10, 7, 9, 12, 11, 8, 13, 16, 15)

function _java_turtle_incoming_reordered(order::Int, dir_order::NTuple{N,Int}) where {N}
    up = _java_turtle_upward_points_f32(order)
    length(up) >= N || error("Invalid Java turtle reorder table for order=$(order): expected at least $(N) points, got $(length(up)).")
    [_normalize3_java(up[i]) for i in dir_order]
end

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

function _java_turtle_incoming_logged(n::Int)
    n == 6 && return _java_turtle_incoming_n6_logged()
    n == 16 && return _java_turtle_incoming_n16_logged()
    # For larger sector counts, compatibility mode uses a Float32 turtle build
    # (seed/trigonometry) which is closer to Java internals than Float64 defaults.
    return get!(_TURTLE_DIR_COMPAT_CACHE, n) do
        order = _java_turtle_order(n)
        order < 0 && error("Unsupported Java turtle sector count: $n. Allowed values: 1, 6, 16, 46, 136, 406.")
        up = _java_turtle_upward_points_f32(order)
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

function _diffuse_weights_java_like(sky::SkyState, turtle::TurtleGrid, sky_ids::Vector{Int}, options::LightOptions)
    n = length(sky_ids)
    n == 0 && return Float64[]

    sun_up = _scene_local_sun_direction(sky, options)
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
    build_turtle(options, sky)::TurtleGrid

Build the directional sky discretization (turtle sectors), optionally adding an
explicit sun sector when `options.all_in_turtle == false`. Diffuse turtle
sectors stay fixed in scene-local coordinates; the geographic sun direction is
rotated into that basis using `options.scene_rotation_deg`.

Set `java_logged_turtle_dirs: true` in the
light config to use compatibility-mode sky directions (exact Java-logged vectors for 6 sectors,
Float32 Java-style construction for larger sector counts).
"""
function build_turtle(options::LightOptions, sky::SkyState)
    n = max(options.turtle_sectors, 1)
    use_logged_dirs = _use_java_logged_turtle_dirs(options)
    dirs =
        if _java_turtle_order(n) >= 0
            use_logged_dirs ? _java_turtle_incoming_logged(n) : _java_turtle_incoming(n)
        else
            _hemisphere_fibonacci_incoming(n)
        end
    sectors = TurtleSector[]
    w = 1.0 / n
    for i in eachindex(dirs)
        push!(sectors, TurtleSector(i, dirs[i], w, :sky))
    end

    if !options.all_in_turtle && sky.sun_elevation_deg > 0.0
        push!(sectors, TurtleSector(length(sectors) + 1, _scene_local_sun_direction(sky, options), 0.0, :sun))
    end
    TurtleGrid(sectors)
end

function _directional_flux_substeps(
    meteo_row,
    sky::SkyState,
    options::LightOptions;
    check_boundaries::Bool=options.check_meteo_boundaries,
    resolved_step::Union{Nothing,ResolvedMeteoStep}=nothing,
)
    resolved =
        resolved_step === nothing ?
        _resolved_meteo_step_or_error(
            meteo_row,
            options;
            check_boundaries=check_boundaries,
        ) : resolved_step

    if resolved.latitude_deg === nothing || !_meteo_has_resolved_solar_interval(resolved)
        return NamedTuple[], 0.0, resolved.duration_seconds / 3600.0
    end

    latitude_rad = deg2rad(resolved.latitude_deg)
    date = resolved.date
    doy = Dates.dayofyear(date)
    start_h = resolved.start_hour
    end_h = resolved.end_hour

    global_from_input = resolved.radiation_source != :clearness

    clearness_provided = resolved.clearness !== nothing
    clearness =
        if clearness_provided
            resolved.clearness
        else
            _clearness_from_global_wm2(sky.ri_sw_f, latitude_rad, doy, start_h, end_h)
        end

    timestep_h = _cfg_radiation_timestep_hours(options)
    daylight_fraction = _daylight_fraction(date, start_h, end_h, timestep_h, latitude_rad)
    ri_sw_during_daylight =
        global_from_input && daylight_fraction > 0.0 ?
        sky.ri_sw_f / daylight_fraction : sky.ri_sw_f

    substeps = _radiation_substeps(
        date,
        start_h,
        end_h,
        latitude_rad,
        timestep_h,
        ri_sw_during_daylight,
        clearness,
        global_from_input,
        clearness_provided,
    )
    return substeps, start_h, end_h
end

"""
    compute_directional_fluxes(meteo_row, sky, turtle, options)::DirectionalFluxes

Java-parity directional flux integration using meteo substeps:
directional fluxes are computed at each substep sun position then averaged over the full meteo step.
"""
function compute_directional_fluxes(
    meteo_row,
    sky::SkyState,
    turtle::TurtleGrid,
    options::LightOptions;
    check_boundaries::Bool=options.check_meteo_boundaries,
    resolved_step::Union{Nothing,ResolvedMeteoStep}=nothing,
)
    resolved =
        resolved_step === nothing ?
        _resolved_meteo_step_or_error(
            meteo_row,
            options;
            check_boundaries=check_boundaries,
        ) : resolved_step
    substeps, start_h, end_h = _directional_flux_substeps(
        meteo_row,
        sky,
        options;
        check_boundaries=check_boundaries,
        resolved_step=resolved,
    )
    step_duration_h = end_h - start_h
    (isempty(substeps) || step_duration_h <= 0.0) && return compute_directional_fluxes(sky, turtle, options)

    explicit_direct_fraction =
        resolved.direct_fraction === nothing ? NaN : resolved.direct_fraction
    provided_sun = resolved.solar_geometry_source == :explicit

    n = length(turtle.sectors)
    par_acc = zeros(Float64, n)
    nir_acc = zeros(Float64, n)
    par_ratio = sky.ri_sw_f > 0.0 ? sky.ri_par_f / sky.ri_sw_f : 0.0
    nir_ratio = sky.ri_sw_f > 0.0 ? sky.ri_nir_f / sky.ri_sw_f : 0.0

    for ss in substeps
        total_w = ss.diffuse_w + ss.direct_w
        total_w > 0.0 || continue
        direct_fraction_sub =
            isnan(explicit_direct_fraction) ?
            clamp(ss.direct_w / total_w, 0.0, 1.0) :
            clamp(explicit_direct_fraction, 0.0, 1.0)
        sun_azimuth = provided_sun ? sky.sun_azimuth_deg : ss.sun_azimuth_deg
        sun_elevation = provided_sun ? sky.sun_elevation_deg : ss.sun_elevation_deg
        sky_sub = SkyState(
            sun_azimuth,
            sun_elevation,
            max(ss.global_w * par_ratio, 0.0),
            max(ss.global_w * nir_ratio, 0.0),
            direct_fraction_sub,
            1.0 - direct_fraction_sub,
        )
        flux_sub = compute_directional_fluxes(sky_sub, turtle, options)
        @inbounds for i in 1:n
            par_acc[i] += flux_sub.par[i] * ss.duration
            nir_acc[i] += flux_sub.nir[i] * ss.duration
        end
    end

    par = par_acc ./ step_duration_h
    nir = nir_acc ./ step_duration_h

    DirectionalFluxes([s.id for s in turtle.sectors], par, nir)
end

"""
    compute_directional_fluxes(sky, turtle, options)::DirectionalFluxes

Project sky-level PAR/NIR irradiance to each turtle sector using Java-compatible diffuse and
direct distribution rules.
"""
function compute_directional_fluxes(sky::SkyState, turtle::TurtleGrid, options::LightOptions)
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
        wdiff = _diffuse_weights_java_like(sky, turtle, sky_ids, options)
        for (k, i) in enumerate(sky_ids)
            par[i] += par_diff * wdiff[k]
            nir[i] += nir_diff * wdiff[k]
        end
    end

    if options.all_in_turtle
        # Java Turtle.directInTurtle parity:
        # distribute direct irradiance by angular overlap between a sun halo and each sector.
        dir_count = length(sky_ids)
        if dir_count > 0
            sun_up = _scene_local_sun_direction(sky, options)
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
            else
                # A full-timestep interval mean may legitimately carry direct
                # energy even when the representative sun is below the local
                # horizon. In that case no upper-hemisphere halo overlaps a
                # turtle sector, so preserve the supplied energy by spreading
                # it over the available sky sectors.
                fallback_weights = [max(turtle.sectors[i].weight, 0.0) for i in sky_ids]
                fallback_sum = sum(fallback_weights)
                fallback_sum <= 0.0 && fill!(fallback_weights, 1.0)
                fallback_sum = sum(fallback_weights)
                for (k, i) in enumerate(sky_ids)
                    w = fallback_weights[k] / fallback_sum
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
