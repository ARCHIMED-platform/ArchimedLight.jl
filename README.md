# ArchimedLight.jl

Julia reimplementation of the ARCHIMED light interception pipeline with a composable, function-first API.

## Current scope
- Scene/config/meteo input pipeline
- Sky + turtle discretization
- First-order interception (CPU raster/z-buffer)
- Iterative scattering (CPU reference)
- Java parity harness (curated light-only fixtures)

Energy balance, transpiration and photosynthesis are intentionally out of scope for now.

## Core API
```julia
using ArchimedLight

cfg = read_light_config("config.yml")
scene = read_scene(cfg.scene)
meteo = read_meteo(cfg.meteo)

sky = compute_sky(first(meteo.rows), cfg)
turtle = build_turtle(cfg, sky)
fluxes = compute_directional_fluxes(sky, turtle, cfg)
first_order = compute_first_order(scene, turtle, fluxes, cfg)
graph = build_scattering_transfer_graph(scene, turtle, first_order, cfg)
scat = compute_scattering(graph, first_order, cfg)
step_seconds = 1800.0 # use your meteo timestep duration in seconds
budget = integrate_light(first_order, scat, cfg; step_duration_seconds=step_seconds, component_area_per_node=scene.total_area_per_node)
```

`LightBudget` now includes Java-style PAR/NIR intercepted and absorbed outputs in both:
- `*_f`: irradiance (`W m^-2`)
- `*_q`: per-component energy per timestep (`J`)

## Single-call pipeline
```julia
step = run_light_step(scene, first(meteo.rows), cfg)
series = run_light_series(scene, meteo, cfg)
# Optional backend kwargs:
step = run_light_step(scene, first(meteo.rows), cfg; interception_backend=RasterCPUBackend(), scattering_backend=RaycastScatteringBackend())

# Java-style component_values.csv export:
write_component_values_csv("output/component_values.csv", scene, step, cfg; meteo_row=first(meteo.rows), step_number=1)

# Java-style scene_values.csv export:
write_scene_values_csv("output/scene_values.csv", scene, series, cfg; meteo_rows=meteo.rows)
```

## Full Example
- Self-contained files and script are under `example/`.
- Run with:

```bash
julia --project=. example/full_featured_example.jl
```

## Stage flexibility
- You can call each stage independently (`compute_sky`, `build_turtle`, `compute_first_order`, `compute_scattering`, ...).
- You can prebuild scattering transfers via `build_scattering_transfer_graph(...)` and reuse them with `compute_scattering(graph, ...)`.
- `compute_sky` follows Java clearness/global conversion and DeJong hourly direct/diffuse partitioning.
- `compute_sky` uses Java-style substep-weighted sun position (`radiation_timestep`) when sun angles are not provided.
- Meteo `#' use: ...` consistency checks for `clearness`/`RI_SW_f`/`RI_PAR_f`/`RI_NIR_f` are enforced like Java.
- `compute_first_order(...; backend=:raster_cpu)` is the current reference backend (`RasterCPUBackend()` also available).
- `compute_scattering(...; mode=:raycast)` / `compute_scattering(...; mode=:links)` keep Java-style mode selection; backend objects are also available (`RaycastScatteringBackend()`, `LinksScatteringBackend()`).
- `pixel_size` is validated with Java parity bounds (`0 < pixel_size <= 0.5` meters).
- `cache_pixel_table: true` (or `save_on_disk: true`) enables on-disk direction projection cache under `<output_directory>/pixel_tables_cache`.
- `build_turtle` follows Java-compatible sector sets for `1, 6, 16, 46, 136, 406`.

## Testing
Run the parity and smoke suite:

```bash
julia --project=. test/runtests.jl
```
