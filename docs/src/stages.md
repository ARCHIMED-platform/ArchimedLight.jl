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
budget = integrate_light(first, scat, cfg)
```

## Single-Call Pipeline
```julia
step = run_light_step(scene, row, cfg)
series = run_light_series(scene, meteo, cfg)
# Optional backend kwargs:
step = run_light_step(scene, row, cfg; interception_backend=RasterCPUBackend(), scattering_backend=RaycastScatteringBackend())
```

## Caching Options
- `cache_radiation: true` in config reuses directional responses across meteo rows in `run_light_series`.
- `save_on_disk: true` (or `cache_pixel_table: true`) stores per-direction projection tables under
  `<output_directory>/pixel_tables_cache`.

## Backends and Modes
- `compute_first_order(...; backend=:raster_cpu)` is the reference backend (`RasterCPUBackend()` also works).
- `compute_scattering(...; mode=:raycast)` and `compute_scattering(...; mode=:links)` are available.
- `compute_scattering(...; backend=RaycastScatteringBackend())` and `compute_scattering(...; backend=LinksScatteringBackend())` are also available.
