module Runner

export run_from_config, write_components_csv

using ..ArchimedLight
using ..Config: load_config
using ..ArchimedLight: compute_interception, PAR, NIR, TIR, InterceptionResult
using ..ArchimedLight: Scene

function write_components_csv(path::AbstractString, scene::Scene, res::InterceptionResult, cfg::InterceptionConfig)
    step_number = 0
    step_duration = cfg.radiation_timestep > 0 ? cfg.radiation_timestep * 60.0 : 0.0
    open(path, "w") do io
        println(io, "component,area,step_number,step_duration,Ri_PAR_0_f,Ri_NIR_0_f,Ri_PAR_0_q,Ri_NIR_0_q,Ra_PAR_0_f,Ra_NIR_0_f,Ra_PAR_0_q,Ra_NIR_0_q")
        for name in sort!(collect(keys(res.absorbed)))
            area = res.areas[name]
            inc = res.incident[name]
            absb = res.absorbed[name]
            ri_par = get(inc, PAR, 0.0)
            ri_nir = get(inc, NIR, 0.0)
            ra_par = get(absb, PAR, 0.0)
            ra_nir = get(absb, NIR, 0.0)
            ri_par_q = step_duration == 0 ? 0.0 : ri_par * area * step_duration
            ri_nir_q = step_duration == 0 ? 0.0 : ri_nir * area * step_duration
            ra_par_q = step_duration == 0 ? 0.0 : ra_par * area * step_duration
            ra_nir_q = step_duration == 0 ? 0.0 : ra_nir * area * step_duration
            println(io, join((name,
                              area,
                              step_number,
                              step_duration,
                              ri_par,
                              ri_nir,
                              ri_par_q,
                              ri_nir_q,
                              ra_par,
                              ra_nir,
                              ra_par_q,
                              ra_nir_q), ','))
        end
    end
end

"""
    run_from_config(yamlpath; outdir=pwd())

Load config, compute interception, and write outputs to outdir.
Writes components.csv with absorbed bands.
Returns the result Dict.
"""
function run_from_config(yamlpath::AbstractString; outdir::AbstractString=pwd())
    scene, sky, cfg = load_config(yamlpath)
    res = compute_interception(scene, sky, cfg)
    mkpath(outdir)
    write_components_csv(joinpath(outdir, "components.csv"), scene, res, cfg)
    return res
end

end # module
