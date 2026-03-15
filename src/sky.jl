import Dates

const _SOLAR_CONSTANT_MJ_M2_MIN = 0.0820
const _CENTER_OF_SOLAR_DISK_RAD = deg2rad(-0.83)
const _SOLAR_TO_PAR = 0.48
const _SOLAR_TO_NIR = 0.52

function _row_value(row, candidates::Vector{Symbol}, default::Float64)
    names = propertynames(row)
    for c in candidates
        if c in names
            v = getproperty(row, c)
            v isa Number && return Float64(v)
            v === missing && continue
            try
                return parse(Float64, string(v))
            catch
            end
        end
    end
    return default
end

function _row_time_value(row, candidates::Vector{Symbol}, default::Dates.Time)
    names = propertynames(row)
    for c in candidates
        if c in names
            v = getproperty(row, c)
            v === missing && continue
            if v isa Dates.Time
                return v
            elseif v isa Dates.DateTime
                return Dates.Time(v)
            end
            s = strip(string(v))
            isempty(s) && continue
            try
                return Dates.Time(s, Dates.DateFormat("HH:MM:SS"))
            catch
                try
                    return Dates.Time(s, Dates.DateFormat("HH:MM"))
                catch
                end
            end
        end
    end
    return default
end

function _parse_date_value(v, default::Dates.Date)
    v === missing && return default
    if v isa Dates.Date
        return v
    elseif v isa Dates.DateTime
        return Dates.Date(v)
    end
    s = strip(string(v))
    isempty(s) && return default
    lowercase(s) == "missing" && return default
    try
        return Dates.Date(s, Dates.DateFormat("yyyy/mm/dd"))
    catch
        try
            return Dates.Date(s, Dates.DateFormat("yyyy-mm-dd"))
        catch
            return default
        end
    end
end

_to_decimal_hour(t::Dates.Time) = Dates.hour(t) + Dates.minute(t) / 60 + Dates.second(t) / 3600

function _row_date(meteo_row)
    names = propertynames(meteo_row)
    if :date in names
        return _parse_date_value(getproperty(meteo_row, :date), Dates.Date(2000, 1, 1))
    end
    if :dayofyear in names
        doy = Int(round(_row_value(meteo_row, [:dayofyear], 1.0)))
        year = Int(round(_row_value(meteo_row, [:year], 2000.0)))
        doy = clamp(doy, 1, 366)
        return Dates.Date(year, 1, 1) + Dates.Day(doy - 1)
    end
    return Dates.Date(2000, 1, 1)
end

function _row_step_hours(meteo_row)
    t0 = _row_time_value(meteo_row, [:hour_start, :hour], Dates.Time(12))
    t1 = _row_time_value(meteo_row, [:hour_end], t0)

    start_h = _to_decimal_hour(t0)
    end_h = _to_decimal_hour(t1)
    end_h < start_h && (end_h += 24.0)

    if end_h == start_h
        names = propertynames(meteo_row)
        if :step_duration in names
            duration_seconds = _positive_duration_seconds(getproperty(meteo_row, :step_duration); field_name="step_duration")
            end_h += duration_seconds / 3600.0
            return start_h, end_h
        end
        duration = _row_value(meteo_row, [:duration], NaN)
        if isfinite(duration) && duration > 0.0
            # Duration column is usually in hours; allow second-based values as a fallback.
            end_h += duration > 24 ? duration / 3600.0 : duration
        else
            end_h += 1e-6
        end
    end

    return start_h, end_h
end

_java_doy_angle(doy::Int) = 2pi * doy / 365.0

function _java_declination(doy::Int)
    asin(
        sin(deg2rad(-23.44)) *
        cos(
            (deg2rad(360) / 365.24) * (doy + 10) +
            (deg2rad(360) / pi) * 0.0167 * sin((deg2rad(360) / 365.24) * (doy - 2)),
        ),
    )
end

function _java_equation_of_time(doy::Int)
    b = 2pi * (doy - 81) / 365.0
    (0.1645 * sin(2b)) - (0.1255 * cos(b)) - (0.025 * sin(b))
end

function _java_corr_factor(doy::Int)
    a = _java_doy_angle(doy)
    1.000110 + 0.034221 * cos(a) + 0.001280 * sin(a) + 0.000719 * cos(2a) + 0.000077 * sin(2a)
end

function _java_sunset_hour_angle(latitude_rad::Float64, declination_rad::Float64)
    # NOTE: Keep Java 2018 precedence exactly for parity:
    #   (a / cos(lat)) * cos(decl)
    # even though the physically common form is a / (cos(lat)*cos(decl)).
    x =
        ((sin(_CENTER_OF_SOLAR_DISK_RAD) - sin(latitude_rad) * sin(declination_rad)) /
         cos(latitude_rad)) * cos(declination_rad)
    acos(clamp(x, -1.0, 1.0))
end

_java_hour_angle(decimal_hour::Float64) = (pi / 12.0) * (decimal_hour - 12.0)

function _java_extra_terrestrial_hourly_mj(latitude_rad::Float64, doy::Int, start_hour::Float64, end_hour::Float64)
    sc = _java_equation_of_time(doy)
    h1 = _java_hour_angle(start_hour + sc)
    h2 = _java_hour_angle(end_hour + sc)

    decl = _java_declination(doy)
    sunset = _java_sunset_hour_angle(latitude_rad, decl)

    rise = max(h1, -sunset)
    set = min(h2, sunset)
    solar_time_angle = set - rise
    solar_time_angle <= 0.0 && return 0.0

    corr = _java_corr_factor(doy)
    ((12 * 60) / pi) * _SOLAR_CONSTANT_MJ_M2_MIN * corr *
    (
        cos(latitude_rad) * cos(decl) * (sin(set) - sin(rise)) +
        solar_time_angle * (sin(latitude_rad) * sin(decl))
    )
end

_watt_to_mj(watts::Float64, duration_hours::Float64) = watts * duration_hours * 0.0036

function _mega_joules_to_watts(mj::Float64, duration_hours::Float64)
    duration_hours <= 0.0 && return 0.0
    mj * 1_000_000 / (duration_hours * 3600.0)
end

function _global_wm2_from_clearness(clearness::Float64, latitude_rad::Float64, doy::Int, start_hour::Float64, end_hour::Float64)
    et_mj = _java_extra_terrestrial_hourly_mj(latitude_rad, doy, start_hour, end_hour)
    _mega_joules_to_watts(clearness * et_mj, end_hour - start_hour)
end

function _clearness_from_global_wm2(global_wm2::Float64, latitude_rad::Float64, doy::Int, start_hour::Float64, end_hour::Float64)
    et_mj = _java_extra_terrestrial_hourly_mj(latitude_rad, doy, start_hour, end_hour)
    if et_mj > 0.0
        global_mj = _watt_to_mj(global_wm2, end_hour - start_hour)
        return global_mj / et_mj
    end
    return 0.5
end

function _java_sun_position_deg(doy::Int, latitude_rad::Float64, decimal_hour::Float64)
    hour_angle = deg2rad((decimal_hour + _java_equation_of_time(doy) - 12) * 15)
    cos_ha = cos(hour_angle)
    cos_lat = cos(latitude_rad)
    sin_lat = sin(latitude_rad)
    decl = _java_declination(doy)
    cos_decl = cos(decl)
    sin_decl = sin(decl)

    elevation = asin((sin_lat * sin_decl) + (cos_lat * cos_decl * cos_ha))

    dy = -sin(hour_angle)
    dx = tan(decl) * cos_lat - sin_lat * cos_ha
    azimuth = atan(dy, dx)
    azimuth < 0.0 && (azimuth += 2pi)

    return mod(rad2deg(azimuth), 360.0), rad2deg(elevation)
end

function _dejong_kd_hourly(clearness::Float64, sun_elevation_rad::Float64)
    if clearness <= 0.22
        return 1.0
    elseif clearness <= 0.35
        d = clearness - 0.22
        return 1.0 - 6.4 * d * d
    else
        sin_el = sin(sun_elevation_rad)
        r = 0.847 - (1.61 * sin_el) + (1.04 * sin_el * sin_el)
        k = (1.47 - r) / 1.66
        return clearness <= k ? 1.47 - 1.66 * clearness : r
    end
end

function _partition_global_hourly(global_wm2::Float64, clearness::Float64, sun_elevation_rad::Float64)
    if global_wm2 <= 0.0
        return 0.0, 0.0
    elseif sun_elevation_rad <= 0.0
        return global_wm2, 0.0
    end

    kd = clamp(_dejong_kd_hourly(clearness, sun_elevation_rad), 0.0, 1.0)
    diffuse = kd * global_wm2
    direct = max(global_wm2 - diffuse, 0.0)
    return diffuse, direct
end

function _as_degrees_if_radians(x::Float64)
    isnan(x) && return x
    abs(x) <= 2pi + 1e-9 ? rad2deg(x) : x
end

function _use_tokens(row)
    if :use in propertynames(row)
        v = getproperty(row, :use)
        if v isa AbstractString
            t = strip(String(v))
            return isempty(t) ? String[] : split(t)
        elseif v isa AbstractVector
            return [strip(string(x)) for x in v if !isempty(strip(string(x)))]
        elseif v !== missing
            t = strip(string(v))
            return isempty(t) ? String[] : split(t)
        end
    end
    return String[]
end

function _has_any_column(row, candidates::Vector{Symbol})
    names = propertynames(row)
    for c in candidates
        c in names && return true
    end
    return false
end

function _normalize_humidity_use_token(token::AbstractString)
    s = lowercase(strip(String(token)))
    if s in ("relativehumidity", "relative_humidity", "rh")
        return "relativeHumidity"
    elseif s == "vpd"
        return "VPD"
    end
    return ""
end

function _validate_meteo_humidity_inputs(use_tokens::AbstractVector{<:AbstractString}, row)
    has_rh = _has_any_column(row, [:relativeHumidity, :relative_humidity, :RH, :rh])
    has_vpd = _has_any_column(row, [:VPD, :vpd])

    uses = String[]
    for u in use_tokens
        uu = _normalize_humidity_use_token(u)
        isempty(uu) || push!(uses, uu)
    end
    uses = unique(uses)

    if any(==("relativeHumidity"), uses) && !has_rh
        error("meteo consistency error: missing relativeHumidity column specified in 'use' line")
    end
    if any(==("VPD"), uses) && !has_vpd
        error("meteo consistency error: missing VPD column specified in 'use' line")
    end

    if has_rh && has_vpd
        isempty(uses) && error("meteo consistency error: missing 'use' for columns relativeHumidity/VPD")
        length(uses) == 1 || error("meteo consistency error: multiple uses: relativeHumidity/VPD")
    elseif !has_rh && !has_vpd
        error("meteo consistency error: missing column relativeHumidity/VPD")
    else
        # Exactly one humidity source is available.
        length(uses) <= 1 || error("meteo consistency error: multiple uses: relativeHumidity/VPD")
    end
end

function _validate_meteo_radiation_inputs(use_tokens::AbstractVector{<:AbstractString}, has_cl::Bool, has_sw::Bool, has_par::Bool, has_nir::Bool)
    key_cl = "clearness"
    keys_set2 = ("RI_SW_f", "RI_PAR_f", "RI_NIR_f")
    all_keys = (key_cl, keys_set2...)
    uses_rel = [u for u in use_tokens if u in all_keys]

    # Java check for RI_SW_f/RI_PAR_f/RI_NIR_f triplet consistency.
    if has_sw && has_par && has_nir
        k = count(u -> u in keys_set2, uses_rel)
        k == 0 && error("meteo consistency error: missing 'use' line with RI_SW_f/RI_PAR_f/RI_NIR_f columns")
        k == 1 || error("meteo consistency error: extra use(s) with RI_SW_f/RI_PAR_f/RI_NIR_f columns")
    end

    # Java checkConsistency(set1=clearness, set2=[RI_SW_f, RI_PAR_f, RI_NIR_f], allowDoubleUse=true)
    has_any = has_cl || has_sw || has_par || has_nir
    has_any || error("meteo consistency error: missing column clearness/RI_SW_f/RI_PAR_f/RI_NIR_f")

    for u in uses_rel
        if u == key_cl && !has_cl
            error("meteo consistency error: missing clearness column specified in 'use' line")
        elseif u == "RI_SW_f" && !has_sw
            error("meteo consistency error: missing RI_SW_f column specified in 'use' line")
        elseif u == "RI_PAR_f" && !has_par
            error("meteo consistency error: missing RI_PAR_f column specified in 'use' line")
        elseif u == "RI_NIR_f" && !has_nir
            error("meteo consistency error: missing RI_NIR_f column specified in 'use' line")
        end
    end

    use1 = any(==("clearness"), uses_rel)
    uses2 = [u for u in uses_rel if u in keys_set2]
    if use1
        length(uses2) <= 1 || error("meteo consistency error: multiple uses: clearness/$(join(uses2, '/'))")
    else
        length(uses2) <= 1 || error("meteo consistency error: multiple uses: $(join(uses2, '/'))")
    end

    if !use1 && isempty(uses2) && has_cl && (has_sw || has_par || has_nir)
        cols = String[]
        has_sw && push!(cols, "RI_SW_f")
        has_par && push!(cols, "RI_PAR_f")
        has_nir && push!(cols, "RI_NIR_f")
        error("meteo consistency error: missing 'use' for columns clearness/$(join(cols, '/'))")
    end
end

function _dict_float_or_nan(d, key::String)
    d isa AbstractDict || return NaN
    haskey(d, key) || return NaN
    v = d[key]
    v isa Number && return Float64(v)
    v === missing && return NaN
    try
        return parse(Float64, string(v))
    catch
        return NaN
    end
end

function _cfg_radiation_timestep_hours(cfg::LightConfig)
    mins = _dict_float_or_nan(cfg.raw, "radiation_timestep")
    isnan(mins) && (mins = 15.0)
    mins <= 0.0 && (mins = 15.0)
    mins / 60.0
end

function _java_sunrise_sunset_hours(latitude_rad::Float64, doy::Int)
    decl = _java_declination(doy)
    hour_angle = _java_sunset_hour_angle(latitude_rad, decl)
    ut = hour_angle / deg2rad(15.0)
    tst = _java_equation_of_time(doy)
    rise = 12.0 - tst - ut
    set = 12.0 - tst + ut
    return rise, set
end

function _sun_direction_from_az_el_deg(azimuth_deg::Float64, elevation_deg::Float64)
    az = deg2rad(azimuth_deg)
    el = deg2rad(elevation_deg)
    x = cos(el) * sin(az)
    y = cos(el) * cos(az)
    z = sin(el)
    return x, y, z
end

function _az_el_from_direction_deg(x::Float64, y::Float64, z::Float64)
    n = sqrt(x * x + y * y + z * z)
    n <= 0.0 && return (0.0, -90.0)
    xn, yn, zn = x / n, y / n, z / n
    az = atan(xn, yn)
    az < 0.0 && (az += 2pi)
    return rad2deg(az), rad2deg(asin(clamp(zn, -1.0, 1.0)))
end

function _java_substeps_v2(date::Dates.Date, start_h::Float64, end_h::Float64, timestep_h::Float64, latitude_rad::Float64)
    timestep_h = max(timestep_h, 1e-6)
    end_h < start_h && (end_h += 24.0)

    substeps = NamedTuple[]
    start_day = floor(Int, start_h / 24.0)
    end_day = floor(Int, (end_h - eps(Float64)) / 24.0)
    for day_off in start_day:end_day
        day_abs0 = 24.0 * day_off
        seg_abs_start = max(start_h, day_abs0)
        seg_abs_end = min(end_h, day_abs0 + 24.0)
        seg_abs_end <= seg_abs_start && continue

        d = date + Dates.Day(day_off)
        doy = Dates.dayofyear(d)
        rise, set = _java_sunrise_sunset_hours(latitude_rad, doy)

        local_start = seg_abs_start - day_abs0
        local_end = seg_abs_end - day_abs0
        s = max(local_start, rise)
        e = min(local_end, set)
        e <= s && continue

        step_duration = e - s
        sub_count = max(1, ceil(Int, step_duration / timestep_h))
        sub_duration = step_duration / sub_count
        h = s
        for _ in 1:sub_count
            push!(substeps, (doy=doy, start=h, stop=h + sub_duration, sun_pos=h + sub_duration / 2, duration=sub_duration))
            h += sub_duration
        end
    end

    return substeps
end

function _midpoint_sun_position_deg(date::Dates.Date, start_h::Float64, end_h::Float64, latitude_rad::Float64)
    end_h < start_h && (end_h += 24.0)
    mid = (start_h + end_h) / 2.0
    day_off = floor(Int, mid / 24.0)
    doy = Dates.dayofyear(date + Dates.Day(day_off))
    local_mid = mid - 24.0 * day_off
    _java_sun_position_deg(doy, latitude_rad, local_mid)
end

function _radiation_substeps(
    date::Dates.Date,
    start_h::Float64,
    end_h::Float64,
    latitude_rad::Float64,
    timestep_h::Float64,
    ri_sw::Float64,
    clearness::Float64,
    global_from_input::Bool,
    clearness_provided::Bool,
)
    substeps = _java_substeps_v2(date, start_h, end_h, timestep_h, latitude_rad)
    rows = NamedTuple[]

    for ss in substeps
        global_w =
            if global_from_input
                max(ri_sw, 0.0)
            elseif clearness_provided
                _global_wm2_from_clearness(clearness, latitude_rad, ss.doy, ss.start, ss.stop)
            else
                max(ri_sw, 0.0)
            end

        clearness_sub =
            if clearness_provided
                clearness
            else
                _clearness_from_global_wm2(global_w, latitude_rad, ss.doy, ss.start, ss.stop)
            end

        sun_az_sub, sun_el_sub = _java_sun_position_deg(ss.doy, latitude_rad, ss.sun_pos)
        diffuse_w, direct_w = _partition_global_hourly(global_w, clearness_sub, deg2rad(sun_el_sub))

        push!(
            rows,
            (
                doy=ss.doy,
                start=ss.start,
                stop=ss.stop,
                duration=ss.duration,
                sun_azimuth_deg=sun_az_sub,
                sun_elevation_deg=sun_el_sub,
                global_w=global_w,
                diffuse_w=diffuse_w,
                direct_w=direct_w,
            ),
        )
    end

    rows
end

function _auto_sun_and_direct_fraction(
    date::Dates.Date,
    start_h::Float64,
    end_h::Float64,
    latitude_rad::Float64,
    timestep_h::Float64,
    ri_sw::Float64,
    clearness::Float64,
    global_from_input::Bool,
    clearness_provided::Bool,
)
    substeps = _radiation_substeps(
        date,
        start_h,
        end_h,
        latitude_rad,
        timestep_h,
        ri_sw,
        clearness,
        global_from_input,
        clearness_provided,
    )
    sx = 0.0
    sy = 0.0
    sz = 0.0
    total_direct_weight = 0.0
    total_direct_energy = 0.0
    total_diffuse_energy = 0.0

    for ss in substeps
        total_diffuse_energy += ss.diffuse_w * ss.duration
        total_direct_energy += ss.direct_w * ss.duration

        if ss.direct_w > 0.0
            x, y, z = _sun_direction_from_az_el_deg(ss.sun_azimuth_deg, ss.sun_elevation_deg)
            sx += x * ss.direct_w
            sy += y * ss.direct_w
            sz += z * ss.direct_w
            total_direct_weight += ss.direct_w
        end
    end

    sun_az, sun_el =
        if total_direct_weight > 0.0
            _az_el_from_direction_deg(sx, sy, sz)
        else
            _midpoint_sun_position_deg(date, start_h, end_h, latitude_rad)
        end

    total_energy = total_direct_energy + total_diffuse_energy
    direct_fraction =
        if total_energy > 0.0
            clamp(total_direct_energy / total_energy, 0.0, 1.0)
        else
            # Nighttime or empty daylight interval.
            diffuse_w, direct_w = _partition_global_hourly(max(ri_sw, 0.0), clearness, deg2rad(sun_el))
            tw = diffuse_w + direct_w
            tw > 0.0 ? clamp(direct_w / tw, 0.0, 1.0) : 0.0
        end

    return sun_az, sun_el, direct_fraction
end

"""
    compute_sky(meteo_row, cfg)::SkyState

Compute sun position, PAR/NIR/SW irradiance, and direct/diffuse partition for one meteo row,
following Java-compatible precedence rules for available meteorological inputs.
"""
function compute_sky(meteo_row, cfg::LightConfig)
    ri_sw_raw = _row_value(meteo_row, [:RI_SW_f, :Rg, :rg, :sw_global, :global], NaN)
    ri_par = _row_value(meteo_row, [:RI_PAR_f, :PAR, :par], NaN)
    ri_nir = _row_value(meteo_row, [:RI_NIR_f, :NIR, :nir], NaN)
    ri_sw = ri_sw_raw
    global_from_input = !isnan(ri_sw_raw) || !isnan(ri_par) || !isnan(ri_nir)

    clearness_raw = _row_value(meteo_row, [:clearness, :Kt], NaN)
    clearness_provided = !isnan(clearness_raw)
    clearness = clearness_raw

    use_tokens = _use_tokens(meteo_row)
    _validate_meteo_humidity_inputs(use_tokens, meteo_row)
    _validate_meteo_radiation_inputs(
        use_tokens,
        !isnan(clearness_raw),
        !isnan(ri_sw_raw),
        !isnan(ri_par),
        !isnan(ri_nir),
    )

    latitude_deg = _row_value(meteo_row, [:latitude, :lat], 0.0)
    latitude_rad = deg2rad(latitude_deg)
    date = _row_date(meteo_row)
    doy = Dates.dayofyear(date)
    start_h, end_h = _row_step_hours(meteo_row)

    sun_azimuth = _row_value(meteo_row, [:sun_azimut, :sun_azimuth], NaN)
    sun_elevation = _row_value(meteo_row, [:sun_elevation], NaN)
    sun_provided = !isnan(sun_azimuth) && !isnan(sun_elevation)
    if !isnan(sun_azimuth)
        sun_azimuth = _as_degrees_if_radians(sun_azimuth)
    end
    if !isnan(sun_elevation)
        sun_elevation = _as_degrees_if_radians(sun_elevation)
    end

    if isnan(ri_sw)
        if !isnan(ri_par) && !isnan(ri_nir)
            ri_sw = ri_par + ri_nir
        elseif !isnan(ri_par)
            ri_sw = ri_par / _SOLAR_TO_PAR
        elseif !isnan(ri_nir)
            ri_sw = ri_nir / _SOLAR_TO_NIR
        elseif !isnan(clearness)
            ri_sw = _global_wm2_from_clearness(clearness, latitude_rad, doy, start_h, end_h)
        else
            ri_sw = 0.0
        end
    end

    if isnan(clearness)
        clearness = _clearness_from_global_wm2(ri_sw, latitude_rad, doy, start_h, end_h)
    end

    auto_sun_az, auto_sun_el, auto_direct_fraction = _auto_sun_and_direct_fraction(
        date,
        start_h,
        end_h,
        latitude_rad,
        _cfg_radiation_timestep_hours(cfg),
        ri_sw,
        clearness,
        global_from_input,
        clearness_provided,
    )
    if !sun_provided
        sun_azimuth = auto_sun_az
        sun_elevation = auto_sun_el
    end

    ri_sw_known = !isnan(ri_sw)
    if isnan(ri_par) && isnan(ri_nir)
        ri_par = _SOLAR_TO_PAR * ri_sw
        ri_nir = _SOLAR_TO_NIR * ri_sw
    elseif isnan(ri_par)
        if ri_sw_known
            # Java ClearnessGlobalRelation computes PAR from global before trying NIR.
            ri_par = _SOLAR_TO_PAR * ri_sw
        else
            ri_par = (ri_nir / _SOLAR_TO_NIR) * _SOLAR_TO_PAR
        end
    elseif isnan(ri_nir)
        if ri_sw_known
            # Java ClearnessGlobalRelation computes NIR from global before trying PAR.
            ri_nir = _SOLAR_TO_NIR * ri_sw
        else
            ri_nir = (ri_par / _SOLAR_TO_PAR) * _SOLAR_TO_NIR
        end
    end

    explicit_direct = _row_value(meteo_row, [:direct_fraction, :fDIR_SW, :Fd], NaN)
    direct_fraction =
        if !isnan(explicit_direct)
            clamp(explicit_direct, 0.0, 1.0)
        else
            auto_direct_fraction
        end
    diffuse_fraction = 1.0 - direct_fraction

    SkyState(
        sun_azimuth,
        sun_elevation,
        ri_sw,
        max(ri_par, 0.0),
        max(ri_nir, 0.0),
        direct_fraction,
        diffuse_fraction
    )
end
