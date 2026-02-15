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

function compute_sky(meteo_row, cfg::LightConfig)
    ri_sw = _row_value(meteo_row, [:RI_SW_f, :Rg, :rg, :sw_global], NaN)
    ri_par = _row_value(meteo_row, [:RI_PAR_f, :PAR, :par], NaN)
    ri_nir = _row_value(meteo_row, [:RI_NIR_f, :NIR, :nir], NaN)

    if isnan(ri_par) && isnan(ri_nir)
        if isnan(ri_sw)
            ri_par = 0.0
            ri_nir = 0.0
        else
            ri_par = 0.48 * ri_sw
            ri_nir = 0.52 * ri_sw
        end
    elseif isnan(ri_par)
        ri_par = (ri_sw - ri_nir)
    elseif isnan(ri_nir)
        ri_nir = (ri_sw - ri_par)
    end

    clearness = _row_value(meteo_row, [:clearness, :Kt], 0.5)
    explicit_direct = _row_value(meteo_row, [:direct_fraction, :fDIR_SW, :Fd], NaN)
    direct_fraction = isnan(explicit_direct) ? clamp(0.15 + 0.75 * clearness, 0.0, 1.0) : clamp(explicit_direct, 0.0, 1.0)
    diffuse_fraction = 1.0 - direct_fraction

    sun_azimuth = _row_value(meteo_row, [:sun_azimut, :sun_azimuth], 180.0)
    sun_elevation = _row_value(meteo_row, [:sun_elevation], 45.0)

    SkyState(
        sun_azimuth,
        sun_elevation,
        max(ri_par, 0.0),
        max(ri_nir, 0.0),
        direct_fraction,
        diffuse_fraction
    )
end
