const _SOLAR_CONSTANT_MJ_M2_MIN = 0.0820
const _CENTER_OF_SOLAR_DISK_RAD = deg2rad(-0.83)
const _SOLAR_TO_PAR = 0.48
const _SOLAR_TO_NIR = 0.52

function _duration_seconds_strict(v; field_name::AbstractString="duration")
    if v isa Dates.Period || v isa Dates.CompoundPeriod
        seconds = Dates.toms(v) * 1.0e-3
        isfinite(seconds) && seconds > 0.0 || error("Invalid $(field_name) value: expected a positive duration, got $(repr(v))")
        return seconds
    end
    PlantMeteo.positive_duration_seconds(v; field_name=field_name)
end

_row_propertynames(row) = row isa PlantMeteo.TimeStepRow ? Tables.columnnames(row) : propertynames(row)

function _row_metadata_value(row, candidates::Vector{Symbol})
    row isa PlantMeteo.TimeStepRow || return nothing
    meta = getfield(parent(row), :metadata)
    for c in candidates
        hasproperty(meta, c) && return getproperty(meta, c)
    end
    return nothing
end

function _row_value(row, candidates::Vector{Symbol}, default::Float64)
    names = _row_propertynames(row)
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
    meta_v = _row_metadata_value(row, candidates)
    if meta_v !== nothing
        meta_v isa Number && return Float64(meta_v)
        meta_v === missing || begin
            try
                return parse(Float64, string(meta_v))
            catch
            end
        end
    end
    return default
end

function _row_time_value(row, candidates::Vector{Symbol}, default::Dates.Time)
    names = _row_propertynames(row)
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
    names = _row_propertynames(meteo_row)
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
    t0 = _row_time_value(meteo_row, [:hour_start, :hour, :date], Dates.Time(12))
    t1 = _row_time_value(meteo_row, [:hour_end], t0)

    start_h = _to_decimal_hour(t0)
    end_h = _to_decimal_hour(t1)
    end_h < start_h && (end_h += 24.0)

    if end_h == start_h
        names = _row_propertynames(meteo_row)
        if :step_duration in names
            duration_seconds = _duration_seconds_strict(getproperty(meteo_row, :step_duration); field_name="step_duration")
            end_h += duration_seconds / 3600.0
            return start_h, end_h
        end
        if :duration in names
            duration_seconds = _duration_seconds_strict(getproperty(meteo_row, :duration); field_name="duration")
            end_h += duration_seconds / 3600.0
        else
            duration = _row_value(meteo_row, [:duration], NaN)
            if isfinite(duration) && duration > 0.0
                # Duration column is usually in hours; allow second-based values as a fallback.
                end_h += duration > 24 ? duration / 3600.0 : duration
            else
                end_h += 1e-6
            end
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
    numerator = sin(_CENTER_OF_SOLAR_DISK_RAD) - sin(latitude_rad) * sin(declination_rad)
    denominator = cos(latitude_rad) * cos(declination_rad)

    # This intentionally corrects the legacy Java 2018 operator precedence.
    # At the geographic poles, classify the constant solar elevation directly
    # instead of dividing by a denominator that is zero in exact arithmetic.
    if abs(denominator) <= eps(Float64)
        return numerator <= 0.0 ? pi : 0.0
    end

    cos_hour_angle = numerator / denominator
    cos_hour_angle <= -1.0 && return pi
    cos_hour_angle >= 1.0 && return 0.0
    return acos(cos_hour_angle)
end

_java_hour_angle(decimal_hour::Float64) = (pi / 12.0) * (decimal_hour - 12.0)

function _java_extra_terrestrial_hourly_mj(latitude_rad::Float64, doy::Int, start_hour::Float64, end_hour::Float64)
    sc = _java_equation_of_time(doy)
    h1 = _java_hour_angle(start_hour + sc)
    h2 = _java_hour_angle(end_hour + sc)

    decl = _java_declination(doy)
    sunset = _java_sunset_hour_angle(latitude_rad, decl)

    iszero(sunset) && return 0.0
    rise, set = sunset == pi ? (h1, h2) : (max(h1, -sunset), min(h2, sunset))
    solar_time_angle = set - rise
    solar_time_angle <= 0.0 && return 0.0

    corr = _java_corr_factor(doy)
    extraterrestrial_mj =
        ((12 * 60) / pi) * _SOLAR_CONSTANT_MJ_M2_MIN * corr *
        (
            cos(latitude_rad) * cos(decl) * (sin(set) - sin(rise)) +
            solar_time_angle * (sin(latitude_rad) * sin(decl))
        )
    return max(extraterrestrial_mj, 0.0)
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

function _has_any_column(row, candidates::Vector{Symbol})
    names = _row_propertynames(row)
    for c in candidates
        c in names && return true
    end
    return false
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

function _cfg_radiation_timestep_hours(options::LightOptions)
    mins = options.radiation_timestep_minutes
    mins <= 0.0 && (mins = 15.0)
    mins / 60.0
end

function _java_sunrise_sunset_hours(latitude_rad::Float64, doy::Int)
    decl = _java_declination(doy)
    hour_angle = _java_sunset_hour_angle(latitude_rad, decl)
    iszero(hour_angle) && return 0.0, 0.0
    hour_angle == pi && return 0.0, 24.0

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

function _daylight_fraction(
    date::Dates.Date,
    start_h::Float64,
    end_h::Float64,
    timestep_h::Float64,
    latitude_rad::Float64,
)
    step_duration = end_h - start_h
    step_duration > 0.0 || return 0.0
    substeps = _java_substeps_v2(date, start_h, end_h, timestep_h, latitude_rad)
    clamp(sum((ss.duration for ss in substeps); init=0.0) / step_duration, 0.0, 1.0)
end

function _has_daylight_interval_metadata(meteo_row)
    names = _row_propertynames(meteo_row)
    has_date = (:date in names) || (:dayofyear in names)
    has_latitude =
        (:latitude in names) || (:lat in names) ||
        _row_metadata_value(meteo_row, [:latitude, :lat]) !== nothing
    date_has_time =
        :date in names && getproperty(meteo_row, :date) isa Dates.DateTime
    has_start = date_has_time || (:hour_start in names) || (:hour in names)
    has_stop =
        (:hour_end in names) || (:step_duration in names) || (:duration in names)
    return has_date && has_latitude && has_start && has_stop
end

function _radiation_input_effective_scale(meteo_row, options::LightOptions)
    options.radiation_input_semantics == :interval_mean && return 1.0
    provided_sun = _provided_sun_position_deg(meteo_row)
    if provided_sun !== nothing && !_has_daylight_interval_metadata(meteo_row)
        # Explicit sun geometry is authoritative when no complete astronomical
        # interval is available from which to integrate sunrise/sunset.
        return provided_sun.elevation > 0.0 ? 1.0 : 0.0
    end
    latitude_rad = deg2rad(_row_value(meteo_row, [:latitude, :lat], 0.0))
    date = _row_date(meteo_row)
    start_h, end_h = _row_step_hours(meteo_row)
    _daylight_fraction(date, start_h, end_h, _cfg_radiation_timestep_hours(options), latitude_rad)
end

function _radiation_input_effective_scale(
    meteo_row,
    options::LightOptions,
    resolved::ResolvedMeteoStep,
)
    options.radiation_input_semantics == :interval_mean && return 1.0
    if resolved.latitude_deg !== nothing && _meteo_has_resolved_solar_interval(resolved)
        latitude_rad = deg2rad(resolved.latitude_deg)
        return _daylight_fraction(
            resolved.date,
            resolved.start_hour,
            resolved.end_hour,
            _cfg_radiation_timestep_hours(options),
            latitude_rad,
        )
    end
    elevation = resolved.sun_elevation_deg
    return elevation !== nothing && elevation > 0.0 ? 1.0 : 0.0
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
    partition_fraction::Union{Nothing,Float64}=nothing,
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
    explicit_fraction =
        partition_fraction === nothing ? nothing : clamp(partition_fraction, 0.0, 1.0)

    for ss in substeps
        total_w = ss.diffuse_w + ss.direct_w
        direct_w =
            explicit_fraction === nothing ? ss.direct_w : total_w * explicit_fraction
        diffuse_w = total_w - direct_w
        total_diffuse_energy += diffuse_w * ss.duration
        total_direct_energy += direct_w * ss.duration

        if direct_w > 0.0
            x, y, z = _sun_direction_from_az_el_deg(ss.sun_azimuth_deg, ss.sun_elevation_deg)
            direct_energy = direct_w * ss.duration
            sx += x * direct_energy
            sy += y * direct_energy
            sz += z * direct_energy
            total_direct_weight += direct_energy
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
        elseif explicit_fraction !== nothing
            explicit_fraction
        else
            # Nighttime or empty daylight interval.
            diffuse_w, direct_w = _partition_global_hourly(max(ri_sw, 0.0), clearness, deg2rad(sun_el))
            tw = diffuse_w + direct_w
            tw > 0.0 ? clamp(direct_w / tw, 0.0, 1.0) : 0.0
        end

    step_duration = end_h - start_h
    effective_global_w = step_duration > 0.0 ? total_energy / step_duration : 0.0

    return sun_az, sun_el, direct_fraction, effective_global_w
end

function _provided_sun_position_deg(meteo_row)
    azimuth = _row_value(meteo_row, [:sun_azimut, :sun_azimuth], NaN)
    elevation = _row_value(meteo_row, [:sun_elevation], NaN)
    (isnan(azimuth) || isnan(elevation)) && return nothing
    return (
        azimuth=azimuth,
        elevation=elevation,
    )
end

function _provided_direct_fraction(meteo_row)
    direct = _row_value(meteo_row, [:direct_fraction, :fDIR_SW, :Fd], NaN)
    if isnan(direct)
        direct_w = _row_value(meteo_row, [:Ri_SW_f_direct, :RI_SW_f_direct], NaN)
        diffuse_w = _row_value(meteo_row, [:Ri_SW_f_diffuse, :RI_SW_f_diffuse], NaN)
        if !isnan(direct_w) && !isnan(diffuse_w)
            total = direct_w + diffuse_w
            total > 0.0 && (direct = direct_w / total)
        end
    end
    return direct
end

"""
    compute_sky(meteo_row, options; check_boundaries=options.check_meteo_boundaries)::SkyState

Compute sun position, effective full-timestep PAR/NIR/SW irradiance, and the
direct/diffuse partition for one meteo row. Meteorological aliases and optional
derivations are resolved once through `ResolvedMeteoStep`; supplied irradiance
then follows `options.radiation_input_semantics`.
"""
function compute_sky(
    meteo_row,
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

    ri_sw = resolved.ri_sw_f
    ri_par = resolved.ri_par_f
    ri_nir = resolved.ri_nir_f
    global_from_input = resolved.radiation_source != :clearness
    clearness_provided = resolved.clearness !== nothing

    sun_provided = resolved.solar_geometry_source == :explicit
    sun_azimuth = resolved.sun_azimuth_deg === nothing ? NaN : resolved.sun_azimuth_deg
    sun_elevation = resolved.sun_elevation_deg === nothing ? NaN : resolved.sun_elevation_deg

    has_astronomical_interval =
        resolved.latitude_deg !== nothing && _meteo_has_resolved_solar_interval(resolved)
    integrated_ri_sw = ri_sw
    auto_direct_fraction = 0.0
    effective_scale = 1.0

    if has_astronomical_interval
        latitude_rad = deg2rad(resolved.latitude_deg)
        date = resolved.date
        doy = Dates.dayofyear(date)
        start_h = resolved.start_hour
        end_h = resolved.end_hour
        clearness =
            clearness_provided ? resolved.clearness :
            _clearness_from_global_wm2(ri_sw, latitude_rad, doy, start_h, end_h)

        radiation_timestep_h = _cfg_radiation_timestep_hours(options)
        daylight_fraction =
            _daylight_fraction(date, start_h, end_h, radiation_timestep_h, latitude_rad)
        ri_sw_during_daylight =
            if global_from_input && options.radiation_input_semantics == :interval_mean
                daylight_fraction > 0.0 ? ri_sw / daylight_fraction : 0.0
            else
                ri_sw
            end

        auto_sun_az, auto_sun_el, auto_direct_fraction, integrated_ri_sw =
            _auto_sun_and_direct_fraction(
                date,
                start_h,
                end_h,
                latitude_rad,
                radiation_timestep_h,
                ri_sw_during_daylight,
                clearness,
                global_from_input,
                clearness_provided,
                resolved.direct_fraction,
            )
        if !sun_provided
            sun_azimuth = auto_sun_az
            sun_elevation = auto_sun_el
        end

        effective_scale =
            if global_from_input
                options.radiation_input_semantics == :interval_mean ? 1.0 : daylight_fraction
            elseif ri_sw > 0.0
                integrated_ri_sw / ri_sw
            else
                0.0
            end
    else
        if isnan(sun_azimuth) || isnan(sun_elevation)
            error("Invalid meteo input: explicit finite sun coordinates are required when date/time/latitude cannot reconstruct solar geometry.")
        end
        if options.radiation_input_semantics == :sunlit_intensity
            effective_scale = sun_elevation > 0.0 ? 1.0 : 0.0
        end
        if resolved.direct_fraction === nothing && ri_sw > 0.0
            resolved.clearness === nothing && error(
                "Invalid meteo input: provide finite `direct_fraction`, direct/diffuse components, or `clearness` when solar partition cannot be reconstructed from date/time/latitude.",
            )
            diffuse_w, direct_w =
                _partition_global_hourly(ri_sw, resolved.clearness, deg2rad(sun_elevation))
            partition_total = diffuse_w + direct_w
            auto_direct_fraction =
                partition_total > 0.0 ? clamp(direct_w / partition_total, 0.0, 1.0) : 0.0
        end
    end

    effective_ri_sw = global_from_input ? ri_sw * effective_scale : integrated_ri_sw
    if !global_from_input && ri_sw <= 0.0 && integrated_ri_sw > 0.0
        ri_par = _SOLAR_TO_PAR * integrated_ri_sw
        ri_nir = _SOLAR_TO_NIR * integrated_ri_sw
    else
        ri_par *= effective_scale
        ri_nir *= effective_scale
    end
    options.nir_interception || (ri_nir = 0.0)

    direct_fraction =
        resolved.direct_fraction === nothing ?
        auto_direct_fraction : clamp(resolved.direct_fraction, 0.0, 1.0)

    return SkyState(
        sun_azimuth,
        sun_elevation,
        max(effective_ri_sw, 0.0),
        max(ri_par, 0.0),
        max(ri_nir, 0.0),
        direct_fraction,
        1.0 - direct_fraction,
    )
end
