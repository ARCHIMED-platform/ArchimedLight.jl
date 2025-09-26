module Config

export load_config

using YAML
using ..ArchimedLight
using ..ArchimedLight.ArchimedIO
using ..ArchimedLight: SkyConfig, InterceptionConfig, set_optics!, OpticalProps, PAR, NIR, TIR
using Meshes: Vec
using ..SkySets: sky_sectors_for
using Dates
using Unitful

const DEFAULT_PAR_FRACTION = 0.45
const SOLAR_PAR_FRACTION = 0.48
const SOLAR_NIR_FRACTION = 0.52
const SOLAR_CONSTANT_MJ = 0.0820

deg2rad(x) = x * (π / 180.0)

@inline function _equation_of_time(doy::Int)
    b = 2π * (doy - 81) / 365.0
    return 0.1645 * sin(2b) - 0.1255 * cos(b) - 0.0250 * sin(b)
end

@inline _day_angle(doy::Int) = 2π * (doy - 1) / 365.0

@inline function _declination(doy::Int)
    γ = _day_angle(doy)
    return 0.006918 - 0.399912 * cos(γ) + 0.070257 * sin(γ) - 0.006758 * cos(2γ) +
           0.000907 * sin(2γ) - 0.002697 * cos(3γ) + 0.00148 * sin(3γ)
end

@inline function _corr_factor(doy::Int)
    γ = _day_angle(doy)
    return 1.000110 + 0.034221 * cos(γ) + 0.001280 * sin(γ) + 0.000719 * cos(2γ) + 0.000077 * sin(2γ)
end

@inline _hour_angle(time_dec::Float64) = (π / 12.0) * (time_dec - 12.0)

@inline function _sunset_hour_angle(lat::Float64, decl::Float64)
    numerator = sin(deg2rad(-0.83)) - sin(lat) * sin(decl)
    denom = cos(lat) * cos(decl)
    val = denom == 0.0 ? 0.0 : numerator / denom
    val = clamp(val, -1.0, 1.0)
    return acos(val)
end

function extraterrestrial_hourly(latitude_deg::Float64, doy::Int, start_hour::Float64, end_hour::Float64)
    lat = deg2rad(latitude_deg)
    decl = _declination(doy)
    sc = _equation_of_time(doy)
    ω1 = _hour_angle(start_hour + sc)
    ω2 = _hour_angle(end_hour + sc)
    ωs = _sunset_hour_angle(lat, decl)
    ω1 = max(ω1, -ωs)
    ω2 = min(ω2, ωs)
    if ω2 <= ω1
        return 0.0
    end
    dr = _corr_factor(doy)
    term1 = cos(lat) * cos(decl) * (sin(ω2) - sin(ω1))
    term2 = (ω2 - ω1) * sin(lat) * sin(decl)
    return ((12 * 60.0) / π) * SOLAR_CONSTANT_MJ * dr * (term1 + term2)
end

decimal_hour(time_str::AbstractString) = begin
    t = try
        Dates.Time(time_str)
    catch
        Dates.Time(string(time_str))
    end
    Dates.hour(t) + Dates.minute(t) / 60.0 + Dates.second(t) / 3600.0
end

function parse_metadata(lines::Vector{String})
    meta = Dict{String,Float64}()
    for line in lines
        stripped = strip(line)
        startswith(stripped, "#") || continue
        stripped = strip(stripped, ['#', '\'', ' '])
        parts = split(stripped, "#"; limit=2)
        body = strip(first(parts))
        kv = split(body, ":"; limit=2)
        length(kv) == 2 || continue
        key = lowercase(strip(kv[1]))
        val_str = strip(kv[2])
        isempty(val_str) && continue
        try
            meta[key] = parse(Float64, val_str)
        catch
            # ignore non-numeric metadata
        end
    end
    return meta
end

# Local helper to strip Unitful quantities to plain Float64 for angular/vector math
function strip_number(x)
    # Prefer Unitful.ustrip when available to remove units; fallback to Float64
    try
        return Unitful.ustrip(x)
    catch
        try
            return Float64(x)
        catch
            return x
        end
    end
end

function global_from_clearness(clearness::Float64, latitude_deg::Float64, doy::Int, start_hour::Float64, end_hour::Float64)
    duration = end_hour - start_hour
    if duration <= 0
        return 0.0
    end
    et = extraterrestrial_hourly(latitude_deg, doy, start_hour, end_hour)
    global_mj = clearness * et
    return global_mj <= 0 ? 0.0 : (global_mj * 1.0e6) / (duration * 3600.0)
end

band_from_string(name::AbstractString) = name == "PAR" ? PAR : name == "NIR" ? NIR : TIR

function parse_model(path::AbstractString)
    data = YAML.load_file(path)
    group = String(get(data, "Group", "scene"))
    specs = Dict{Tuple{String,String},OpticalProps}()
    paving = Dict{Tuple{String,String},Int}()
    types = get(data, "Type", Dict())
    for (typename, info) in types
        interception = get(info, "Interception", Dict())
        transparency = Float64(get(interception, "transparency", 0.0))
        rho = Dict{Band,Float64}(PAR => 0.1, NIR => 0.4, TIR => 0.0)
        tau = Dict{Band,Float64}(PAR => 0.05, NIR => 0.4, TIR => 0.0)
        optics = get(interception, "optical_properties", Dict())
        # per-band scattering factors parsed from optical_properties when scalar values present
        local_scats = Dict{Band,Float64}()
        for (bandname, value) in optics
            band = band_from_string(String(bandname))
            # value in existing test YAMLs is a single number interpreted as scattering factor
            # but historically we used optical_properties to store rho; interpret as scattering per Java OpticalPropertiesModel
            # If value is a scalar, set scattering factor for this band; if it's a dict, fallback to rho/tau parsing
            if isa(value, Number)
                # interpret scalar as scattering factor σ = ρ + τ, with ρ = τ = σ/2 (Archimed convention)
                σ = Float64(value)
                rho[band] = σ / 2.0
                tau[band] = σ / 2.0
                local_scats[band] = σ
            else
                # fallback: parse nested rho/tau structure
                for (k, v) in value
                    b = band_from_string(String(k))
                    rho[b] = Float64(v)
                end
            end
        end
        op = OpticalProps(transparency=transparency)
        # copy rho/tau into op
        for b in (PAR, NIR, TIR)
            op.rho[b] = rho[b]
            op.tau[b] = tau[b]
        end
        # assign scattering if parsed
        if !isempty(local_scats)
            for (b, sf) in local_scats
                op.scattering[b] = sf
            end
        end
        specs[(group, String(typename))] = op
        if haskey(info, "plot_paving")
            paving[(group, String(typename))] = Int(info["plot_paving"])
        end
    end
    return specs, paving
end

function merge_model_specs!(dest::Dict{Tuple{String,String},OpticalProps}, src::Dict{Tuple{String,String},OpticalProps})
    for (k, v) in src
        dest[k] = v
    end
end

function parse_meteo(path::AbstractString)
    lines = readlines(path)
    data_lines = [strip(line) for line in lines if !isempty(strip(line)) && !startswith(strip(line), "#")]
    if length(data_lines) < 2
        return nothing
    end
    header = split(data_lines[1], ';')
    values = split(data_lines[2], ';')
    dict = Dict{String,String}()
    for (k, v) in zip(header, values)
        dict[strip(k)] = strip(v)
    end
    metadata = parse_metadata(lines)

    start_time = get(dict, "hour_start", "00:00:00")
    end_time = get(dict, "hour_end", "00:00:00")
    start_hour = decimal_hour(start_time)
    end_hour = decimal_hour(end_time)
    duration = try
        t0 = DateTime("1970-01-01T" * start_time)
        t1 = DateTime("1970-01-01T" * end_time)
        max(Dates.value(t1 - t0) / 60_000.0, 0.0)
    catch
        0.0
    end
    par = nothing
    nir = nothing
    if haskey(dict, "RI_PAR_0_f")
        par = parse(Float64, dict["RI_PAR_0_f"])
    end
    if haskey(dict, "RI_NIR_0_f")
        nir = parse(Float64, dict["RI_NIR_0_f"])
    end
    if haskey(dict, "RI_SW_f")
        sw = parse(Float64, dict["RI_SW_f"])
        par = sw * DEFAULT_PAR_FRACTION
        nir = sw * (1 - DEFAULT_PAR_FRACTION)
    end
    if (par === nothing || nir === nothing) && haskey(dict, "clearness") && haskey(dict, "date")
        lat = get(metadata, "latitude", 0.0)
        clearness = parse(Float64, dict["clearness"])
        doy = try
            Dates.dayofyear(Date(dict["date"], dateformat"yyyy/mm/dd"))
        catch
            try
                Dates.dayofyear(Date(dict["date"], dateformat"yyyy-mm-dd"))
            catch
                1
            end
        end
        global_w = global_from_clearness(clearness, lat, doy, start_hour, end_hour)
        if global_w > 0
            par = global_w * SOLAR_PAR_FRACTION
            nir = global_w * SOLAR_NIR_FRACTION
        end
    end
    if par === nothing
        par = 400.0
    end
    if nir === nothing
        nir = 400.0
    end
    # compute doy if possible from metadata/date
    doy = try
        if haskey(dict, "date")
            Dates.dayofyear(Date(dict["date"], dateformat"yyyy/mm/dd"))
        else
            1
        end
    catch
        1
    end
    return (par=par, nir=nir, duration=duration, metadata=metadata, start_hour=start_hour, end_hour=end_hour, doy=doy)
end

"""
    load_config(path)

Parse a minimal YAML config and return (scene, sky::SkyConfig, cfg::InterceptionConfig).
Supported keys:
- scene: path to .opf (required)
- sky_sectors: Int (default 16)
- pixel_size: Float (default 0.25)
- scattering: Bool (default false)
- optics: Dict of component_name => Dict(band => (rho, tau)) [optional]
"""
function load_config(path::AbstractString)
    data = YAML.load_file(String(path))

    scene_path = get(data, "scene", nothing)
    scene_path === nothing && error("config: 'scene' is required")
    # Resolve relative to config file directory
    basedir = dirname(String(path))
    spath = isabspath(scene_path) ? scene_path : joinpath(basedir, scene_path)

    model_specs = Dict{Tuple{String,String},OpticalProps}()
    paving_specs = Dict{Tuple{String,String},Int}()
    for rel in get(data, "models", Any[])
        model_path = joinpath(basedir, String(rel))
        specs, paving = parse_model(model_path)
        merge_model_specs!(model_specs, specs)
        for (k, v) in paving
            paving_specs[k] = v
        end
    end

    # Load OPS or OPF
    if endswith(lowercase(spath), ".ops")
        scene = ArchimedIO.load_ops_scene(spath; paving_specs=paving_specs)
    else
        scene = ArchimedIO.load_opf_scene(spath)
    end

    nsectors = Int(get(data, "sky_sectors", 16))
    sectors = sky_sectors_for(nsectors)
    # simple irradiance defaults; can be extended to per-band in config later
    sky = SkyConfig(sectors, Dict(PAR => 400.0, NIR => 400.0, TIR => 0.0))

    px_cm = Float64(get(data, "pixel_size", 0.25))
    px = px_cm / 100.0
    sca = Bool(get(data, "scattering", false))
    ait = Bool(get(data, "all_in_turtle", false))
    rt = Float64(get(data, "radiation_timestep", 0.0))
    cfg = InterceptionConfig(pixel_size=px, scattering=sca, all_in_turtle=ait, radiation_timestep=rt)

    for comp in scene.components
        key = (comp.group, comp.ctype)
        if haskey(model_specs, key)
            set_optics!(comp.name, model_specs[key])
        end
    end

    if haskey(data, "meteo")
        meteo_path = joinpath(basedir, String(data["meteo"]))
        if isfile(meteo_path)
            info = parse_meteo(meteo_path)
            if info !== nothing
                # create base SkyConfig using meteo totals
                sky = SkyConfig(sectors, Dict(PAR => info.par, NIR => info.nir, TIR => 0.0))
                cfg = InterceptionConfig(pixel_size=px, scattering=sca, all_in_turtle=ait, radiation_timestep=info.duration)

                # approximate sun direction for the time window (midpoint)
                meta = info.metadata
                doy = get(info, :doy, info.doy)
                mid_hour = (info.start_hour + info.end_hour) / 2.0
                lat_rad = deg2rad(get(meta, "latitude", 0.0))
                decl = _declination(doy)
                hour_ang = _hour_angle(mid_hour + _equation_of_time(doy))
                # sun vector in local coords (using standard spherical conversion)
                zen = acos(sin(lat_rad) * sin(decl) + cos(lat_rad) * cos(decl) * cos(hour_ang))
                if zen >= π / 2
                    sun_dir = Vec(0.0, 0.0, -1.0)
                else
                    az = atan((cos(decl) * sin(hour_ang)), (cos(lat_rad) * sin(decl) - sin(lat_rad) * cos(decl) * cos(hour_ang)))
                    sx = sin(zen) * cos(az)
                    sy = sin(zen) * sin(az)
                    sz = cos(zen)
                    # DEBUG: types/values for sx, sy, sz
                    @info "sun angles types" t_sx = typeof(sx) t_sy = typeof(sy) t_sz = typeof(sz)
                    @info "sun angles values" v_sx = sx v_sy = sy v_sz = sz
                    sun_dir = Vec(sx, sy, sz)
                end

                # Build a fresh unitless sun direction Vec from the computed angles (avoid Unitful leaks)
                if zen >= π / 2
                    sun_dir = Vec(0.0, 0.0, -1.0)
                else
                    sun_dir = Vec(Float64(sx), Float64(sy), Float64(sz))
                end

                # split global into direct/diffuse using DeJong-like Kt estimate per band
                # Note: `info.par` and `info.nir` are in W/m^2 (instantaneous). `extraterrestrial_hourly`
                # returns energy over the period in MJ/m^2. Convert global W/m^2 -> MJ/m^2 over the period
                # before forming the clearness index Kt = (global_energy / extraterrestrial_energy).
                duration_hours = info.end_hour - info.start_hour
                duration_hours = max(duration_hours, 0.0)
                for (band, bandval) in ((PAR, info.par), (NIR, info.nir))
                    global_flux = bandval
                    et = extraterrestrial_hourly(get(meta, "latitude", 0.0), doy, info.start_hour, info.end_hour)
                    # convert instantaneous W/m^2 to MJ/m^2 over the timestep
                    global_mj = (global_flux * (duration_hours * 3600.0)) / 1.0e6
                    kt = et == 0.0 ? 0.0 : (global_mj / et)
                    # compute sun elevation (radians) from sun_dir z-component (sun_dir[3] == sin(elevation))
                    # Inspect sun_dir z-component before/after stripping
                    s3_raw = sun_dir[3]
                    @info "sun_dir[3] before strip" t = typeof(s3_raw) v = s3_raw
                    s3_stripped = strip_number(s3_raw)
                    @info "sun_dir[3] after strip" t = typeof(s3_stripped) v = s3_stripped
                    sun_elev = asin(clamp(s3_stripped, -1.0, 1.0))
                    kd = ArchimedLight.dejong_kd_hourly(kt, sun_elev)
                    diffuse = kd * global_flux
                    direct = max(0.0, global_flux - diffuse)
                    ArchimedLight.populate_sector_fluxes!(sky, band, global_flux, direct, diffuse, sun_dir, ait)
                end
            end
        end
    end

    if haskey(data, "optics")
        for (comp, bands) in data["optics"]
            op = OpticalProps()
            for (bandname, vals) in bands
                b = bandname == "PAR" ? PAR : bandname == "NIR" ? NIR : TIR
                op.rho[b] = Float64(vals["rho"])
                op.tau[b] = Float64(vals["tau"])
            end
            set_optics!(String(comp), op)
        end
    end

    return scene, sky, cfg
end

end # module
