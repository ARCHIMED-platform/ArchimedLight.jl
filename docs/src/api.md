# API Reference

This page is a compact guide to the public entry points of `ArchimedLight.jl`.
The package is intentionally organized around a small number of composable stages rather than a large object-oriented surface.

## Input Loading

Use these when your workflow starts from files:

```julia
read_config(path)
read_scene(path)
read_models(path_or_paths)
read_options(path)
read_meteo(path)
```

The convenience entry point is usually:

```julia
options, scene, meteo, models = read_config("config.yml")
```

## In-Memory Preparation

Use these when inputs already exist as Julia data structures:

```julia
prepare_scene(mtg; source_path="interactive.opf", scene_xy_bounds=nothing, relabel_ids=false)
prepare_models(models_or_groups)
prepare_meteo(meteo, options)
```

Two scene-editing helpers are especially useful:

```julia
add_ground!(scene; z=0.0, nx=9, ny=9, xy_bounds=nothing, group="pavement", type="Cobblestone")
write_scene(path, scene)
```

## Light Pipeline

The explicit stage API is:

```julia
compute_sky(row, options)
build_turtle(options, sky)
compute_directional_fluxes(row, sky, turtle, options)
compute_first_order(scene, models, turtle, fluxes, options)
compute_scattering(scene, models, turtle, first, options)
integrate_light(scene, models, first, scattering, options; meteo_row=row)
```

For interactive synthetic scenes, `compute_directional_fluxes` also accepts a prebuilt sky state:

```julia
compute_directional_fluxes(sky, turtle, options)
```

## One-Call Helpers

When you do not need manual control over the intermediate stages:

```julia
run_light_step(scene, models, meteo_row, options; interception_backend=:raster_cpu, scattering_mode=:raycast, scattering_backend=nothing)
run_light_series(scene, models, meteo, options; interception_backend=:raster_cpu, scattering_mode=:raycast, scattering_backend=nothing)
```

## Scene Attachment Helpers

These functions attach computed values back onto the MTG using ARCHIMED attribute names:

```julia
attach_node_values!(scene, attr, values; fill_value=nothing)
attach_light_step!(scene, step; fields=[:incident_par_flux], names=Dict(), fill_value=nothing)
attach_light_series!(scene, steps; fields=[:incident_par_flux], names=Dict(), fill_value=NaN)
```

## Visualization Helpers

These helpers expose direct mesh coloring and the Makie package extension:

```julia
light_render_geometry(scene, models, options)
light_render_geometry(step)
light_render_geometry(steps)
tile_light_geometry(scene, geometry; nx=1, ny=1, centered=true, xperiod=nothing, yperiod=nothing)
tile_light_geometry(scene, step; nx=1, ny=1, centered=true, xperiod=nothing, yperiod=nothing)
tile_light_geometry(scene, steps; nx=1, ny=1, centered=true, xperiod=nothing, yperiod=nothing)
tile_light_geometry(scene, models, options; nx=1, ny=1, centered=true, xperiod=nothing, yperiod=nothing)
light_metric_values(step, selector)
light_metric_values(steps, selector; timestep=1)
light_face_values(data; color=:incident_par_flux, timestep=1, fill_value=NaN)
light_vertex_values(data; color=:incident_par_flux, timestep=1, fill_value=NaN)
light_face_values(scene, models, options, data; color=:incident_par_flux, timestep=1, fill_value=NaN)
light_vertex_values(scene, models, options, data; color=:incident_par_flux, timestep=1, fill_value=NaN)
lightplot(geometry, data; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot!(axis, geometry, data; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot(step; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot(steps; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot!(axis, step; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot!(axis, steps; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot(scene, models, options, data; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot!(axis, scene, models, options, data; color=:incident_par_flux, timestep=1, interpolate=false, ...)
```

For `lightplot(steps; colorrange=automatic)`, the Makie extension keeps a
single color range across the whole series so animated timesteps remain
comparable.

## Backend Types

The current public backend selectors are:

```julia
RasterCPUBackend()
RaycastScatteringBackend()
LinksScatteringBackend()
```

They can be passed explicitly to the pipeline helpers when you want to control the implementation used for interception or scattering.

## Recommended Starting Points

- read [Getting Started](getting_started.md) for the shortest runnable workflow
- read [Composable Stages](stages.md) if you want the stage-by-stage pattern
- read [Outputs](outputs.md) for the mapping between `LightBudget` fields and ARCHIMED attribute names


## API List

```@autodocs
Modules = [ArchimedLight]
```
