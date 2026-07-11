using ArchimedLight

"""
    run_example(; output_dir=joinpath(@__DIR__, "output"), io=stdout)

Run one light-interception step from the files in `example_1/` and write a
component table plus a scene carrying the computed light attributes.
"""
function run_example(; output_dir=joinpath(@__DIR__, "output"), io::IO=stdout)
    config_path = joinpath(@__DIR__, "config.yml")
    sim, meteo = read_simulation(config_path)

    println(io, "Scene")
    show(io, summarize_scene(sim.scene; models=sim.models))
    println(io, "Meteorology")
    show(io, summarize_meteo(meteo; options=sim.options))

    step = run_light(sim, first(meteo))
    println(io, "Result")
    println(io, step)

    attach_light_step!(
        sim.scene,
        step;
        fields=[:incident_par_flux, :incident_par_energy, :absorbed_par_energy],
    )

    mkpath(output_dir)
    component_path = joinpath(output_dir, "component_values.csv")
    scene_path = joinpath(output_dir, "scene_with_light.opf")
    write_component_values(component_path, sim, step)
    write_scene(scene_path, sim.scene)

    println(io, "Wrote:")
    println(io, "  ", component_path)
    println(io, "  ", scene_path)

    return step
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_example()
end
