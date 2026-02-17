# Composable Stages

`ArchimedLight.jl` is designed so each stage can be called independently.

## Inputs
```julia
using ArchimedLight

cfg = read_light_config("config.yml")
scene = read_scene(cfg.scene)
meteo = read_meteo(cfg.meteo)
row = first(meteo.rows)
```

## Stage-by-stage Execution
```julia
sky = compute_sky(row, cfg)
turtle = build_turtle(cfg, sky)
fluxes = compute_directional_fluxes(sky, turtle, cfg)
first = compute_first_order(scene, turtle, fluxes, cfg)
graph = build_scattering_transfer_graph(scene, turtle, first, cfg)
scat = compute_scattering(graph, first, cfg)
step_seconds = 1800.0 # use your meteo timestep duration in seconds
budget = integrate_light(first, scat, cfg; step_duration_seconds=step_seconds, component_area_per_node=scene.total_area_per_node)
```

## Single-Call Pipeline
```julia
step = run_light_step(scene, row, cfg)
series = run_light_series(scene, meteo, cfg)
# Optional backend kwargs:
step = run_light_step(scene, row, cfg; interception_backend=RasterCPUBackend(), scattering_backend=RaycastScatteringBackend())

# Optional Java-style component output export:
write_component_values_csv("output/component_values.csv", scene, step, cfg; meteo_row=row, step_number=0)

# Optional Java-style scene output export (for a series):
write_scene_values_csv("output/scene_values.csv", scene, series, cfg; meteo_rows=meteo.rows)

# Optional Java-style logs:
write_sun_position_log_csv("output/log-sun-position.csv", series, meteo.rows)
write_scattering_iteration_log_csv("output/log-iteration-scat-par.csv", scene, step, cfg; meteo_row=row, band="PAR")

# High-level Java-style output writer (config-driven defaults, overridable):
write_light_outputs(scene, series, cfg; meteo_rows=meteo.rows, outdir="output")

# Java-style simulation directory resolution:
sim_out = simulation_output_directory(cfg)  # e.g. <output_directory>/000001
```

## Caching Options
- `cache_radiation: true` in config reuses directional responses across meteo rows in `run_light_series`.
- `save_on_disk: true` (or `cache_pixel_table: true`) stores per-direction projection tables under
  `<output_directory>/pixel_tables_cache`.

## Backends and Modes
- `compute_first_order(...; backend=:raster_cpu)` is the reference backend (`RasterCPUBackend()` also works).
- `compute_scattering(...; mode=:raycast)` and `compute_scattering(...; mode=:links)` are available.
- `compute_scattering(...; backend=RaycastScatteringBackend())` and `compute_scattering(...; backend=LinksScatteringBackend())` are also available.
