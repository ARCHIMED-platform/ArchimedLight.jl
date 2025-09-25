module Runner

export run_from_config, write_components_csv

using ..ArchimedLight
using ..Config: load_config
using ..ArchimedLight: compute_interception, PAR, NIR, TIR
using DelimitedFiles

function write_components_csv(path::AbstractString, res::Dict{String,<:Dict})
    open(path, "w") do io
        println(io, "component,PAR,NIR,TIR")
        for (name, bands) in res
            par = get(bands, PAR, 0.0)
            nir = get(bands, NIR, 0.0)
            tir = get(bands, TIR, 0.0)
            println(io, string(name, ",", par, ",", nir, ",", tir))
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
    write_components_csv(joinpath(outdir, "components.csv"), res)
    return res
end

end # module
