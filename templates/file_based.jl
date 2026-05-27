using ArchimedLight

sim, meteo = read_simulation("config.yml")

summarize_scene(sim.scene; models=sim.models)
summarize_meteo(meteo; options=sim.options)

report = check_simulation(sim.scene, meteo; models=sim.models, options=sim.options)
isempty(report.errors) || error(join(report.errors, "\n"))

series = run_light(sim, meteo)
step = first(series)
