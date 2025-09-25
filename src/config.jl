module Config

export load_config

using YAML
using ..ArchimedLight
using ..ArchimedLight.ArchimedIO
using ..ArchimedLight: SkyConfig, InterceptionConfig, set_optics!, OpticalProps, PAR, NIR, TIR
using ..SkySets: sky_sectors_for
using Dates

const DEFAULT_PAR_FRACTION = 0.45

band_from_string(name::AbstractString) = name == "PAR" ? PAR : name == "NIR" ? NIR : TIR

function parse_model(path::AbstractString)
    data = YAML.load_file(path)
    group = String(get(data, "Group", "scene"))
    specs = Dict{Tuple{String,String},OpticalProps}()
    types = get(data, "Type", Dict())
    for (typename, info) in types
        interception = get(info, "Interception", Dict())
        transparency = Float64(get(interception, "transparency", 0.0))
        rho = Dict{Band,Float64}(PAR=>0.1, NIR=>0.4, TIR=>0.0)
        tau = Dict{Band,Float64}(PAR=>0.05, NIR=>0.4, TIR=>0.0)
        optics = get(interception, "optical_properties", Dict())
        for (bandname, value) in optics
            band = band_from_string(String(bandname))
            rho[band] = Float64(value)
            tau[band] = get(tau, band, 0.0)
        end
        specs[(group, String(typename))] = OpticalProps(rho, tau, transparency)
    end
    return specs
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
    start_time = get(dict, "hour_start", "00:00:00")
    end_time = get(dict, "hour_end", "00:00:00")
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
    if par === nothing
        par = 400.0
    end
    if nir === nothing
        nir = 400.0
    end
    return (par=par, nir=nir, duration=duration)
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
    # Load OPS or OPF
    if endswith(lowercase(spath), ".ops")
        scene = ArchimedIO.load_ops_scene(spath)
    else
        scene = ArchimedIO.load_opf_scene(spath)
    end

    nsectors = Int(get(data, "sky_sectors", 16))
    sectors = sky_sectors_for(nsectors)
    # simple irradiance defaults; can be extended to per-band in config later
    sky = SkyConfig(sectors, Dict(PAR=>400.0, NIR=>400.0, TIR=>0.0))

    px = Float64(get(data, "pixel_size", 0.25))
    sca = Bool(get(data, "scattering", false))
    ait = Bool(get(data, "all_in_turtle", false))
    rt  = Float64(get(data, "radiation_timestep", 0.0))
    cfg = InterceptionConfig(pixel_size=px, scattering=sca, all_in_turtle=ait, radiation_timestep=rt)

    model_specs = Dict{Tuple{String,String},OpticalProps}()
    for rel in get(data, "models", Any[])
        model_path = joinpath(basedir, String(rel))
        merge_model_specs!(model_specs, parse_model(model_path))
    end

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
                sky = SkyConfig(sectors, Dict(PAR=>info.par, NIR=>info.nir, TIR=>0.0))
                cfg = InterceptionConfig(pixel_size=px, scattering=sca, all_in_turtle=ait, radiation_timestep=info.duration)
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
