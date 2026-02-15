import Dates

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

function _row_datetime(row)
    names = propertynames(row)
    dval =
        if :date in names
            getproperty(row, :date)
        else
            "2000-01-01"
        end
    hstart =
        if :hour_start in names
            getproperty(row, :hour_start)
        elseif :hour in names
            getproperty(row, :hour)
        else
            "12:00:00"
        end
    hend = :hour_end in names ? getproperty(row, :hour_end) : nothing

    d =
        if dval isa Dates.Date
            dval
        elseif dval isa Dates.DateTime
            Dates.Date(dval)
        else
            s = strip(string(dval))
            if isempty(s)
                Dates.Date(2000, 1, 1)
            else
                try
                    Dates.Date(s, Dates.DateFormat("yyyy/mm/dd"))
                catch
                    try
                        Dates.Date(s, Dates.DateFormat("yyyy-mm-dd"))
                    catch
                        Dates.Date(2000, 1, 1)
                    end
                end
            end
        end

    function parse_time(v)
        if v isa Dates.Time
            return v
        end
        s = strip(string(v))
        if isempty(s)
            return Dates.Time(12)
        end
        try
            return Dates.Time(s, Dates.DateFormat("HH:MM:SS"))
        catch
            try
                return Dates.Time(s, Dates.DateFormat("HH:MM"))
            catch
                return Dates.Time(12)
            end
        end
    end

    t0 = parse_time(hstart)
    t =
        if isnothing(hend)
            t0
        else
            t1 = parse_time(hend)
            dt0 = Dates.DateTime(d, t0)
            dt1 = Dates.DateTime(d, t1)
            dt1 < dt0 && (dt1 += Dates.Day(1))
            dt0 + Dates.Millisecond(round(Int, Dates.value(dt1 - dt0) / 2))
        end
    t isa Dates.DateTime ? t : Dates.DateTime(d, t)
end

function _solar_position_deg(dt::Dates.DateTime, latitude_deg::Float64, longitude_deg::Float64=0.0)
    # Lightweight solar geometry approximation suitable for fallback meteo files.
    n = Dates.dayofyear(Dates.Date(dt))
    hour = Dates.hour(dt) + Dates.minute(dt) / 60 + Dates.second(dt) / 3600

    γ = 2pi / 365 * (n - 1 + (hour - 12) / 24)
    eqtime = 229.18 * (
        0.000075 + 0.001868 * cos(γ) - 0.032077 * sin(γ) -
        0.014615 * cos(2γ) - 0.040849 * sin(2γ)
    )
    decl = (
        0.006918 - 0.399912 * cos(γ) + 0.070257 * sin(γ) -
        0.006758 * cos(2γ) + 0.000907 * sin(2γ) -
        0.002697 * cos(3γ) + 0.00148 * sin(3γ)
    )

    # Assume meteo times are local solar-ish time if timezone is unavailable.
    time_offset_min = eqtime + 4 * longitude_deg
    tst = hour * 60 + time_offset_min
    ha = deg2rad(tst / 4 - 180)

    φ = deg2rad(latitude_deg)
    cos_zenith = clamp(sin(φ) * sin(decl) + cos(φ) * cos(decl) * cos(ha), -1.0, 1.0)
    zen = acos(cos_zenith)
    elev = 90.0 - rad2deg(zen)

    y = sin(ha)
    x = cos(ha) * sin(φ) - tan(decl) * cos(φ)
    az = rad2deg(atan(y, x)) + 180.0
    az = mod(az, 360.0)

    return az, elev
end

function _derive_ri_sw(clearness::Float64, sun_elevation_deg::Float64)
    if sun_elevation_deg <= 0.0
        return 0.0
    end
    s0 = 1367.0
    exo = s0 * sin(deg2rad(sun_elevation_deg))
    clamp(clearness, 0.0, 1.2) * max(exo, 0.0)
end

function compute_sky(meteo_row, cfg::LightConfig)
    ri_sw = _row_value(meteo_row, [:RI_SW_f, :Rg, :rg, :sw_global], NaN)
    ri_par = _row_value(meteo_row, [:RI_PAR_f, :PAR, :par], NaN)
    ri_nir = _row_value(meteo_row, [:RI_NIR_f, :NIR, :nir], NaN)

    clearness = _row_value(meteo_row, [:clearness, :Kt], 0.5)
    latitude = _row_value(meteo_row, [:latitude, :lat], 0.0)
    longitude = _row_value(meteo_row, [:longitude, :lon], 0.0)

    sun_azimuth = _row_value(meteo_row, [:sun_azimut, :sun_azimuth], NaN)
    sun_elevation = _row_value(meteo_row, [:sun_elevation], NaN)
    if isnan(sun_azimuth) || isnan(sun_elevation)
        dt = _row_datetime(meteo_row)
        sun_azimuth, sun_elevation = _solar_position_deg(dt, latitude, longitude)
    end

    if isnan(ri_sw)
        ri_sw = _derive_ri_sw(clearness, sun_elevation)
    end

    if isnan(ri_par) && isnan(ri_nir)
        ri_par = 0.48 * ri_sw
        ri_nir = 0.52 * ri_sw
    elseif isnan(ri_par)
        ri_par = max(ri_sw - ri_nir, 0.0)
    elseif isnan(ri_nir)
        ri_nir = max(ri_sw - ri_par, 0.0)
    end

    explicit_direct = _row_value(meteo_row, [:direct_fraction, :fDIR_SW, :Fd], NaN)
    direct_fraction =
        if sun_elevation <= 0.0
            0.0
        elseif isnan(explicit_direct)
            clamp(0.15 + 0.75 * clearness, 0.0, 1.0)
        else
            clamp(explicit_direct, 0.0, 1.0)
        end
    diffuse_fraction = 1.0 - direct_fraction

    SkyState(
        sun_azimuth,
        sun_elevation,
        max(ri_par, 0.0),
        max(ri_nir, 0.0),
        direct_fraction,
        diffuse_fraction
    )
end
