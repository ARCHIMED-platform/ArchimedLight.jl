module Config

export load_config

using YAML
using ..ArchimedLight
using ..ArchimedLight.ArchimedIO
using ..ArchimedLight: SkyConfig, InterceptionConfig, set_optics!, OpticalProps, PAR, NIR, TIR
using ..SkySets: sky_sectors_for

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
