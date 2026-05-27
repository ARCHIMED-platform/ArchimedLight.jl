# API Reference

This page is a compact guide to the public entry points of `ArchimedLight.jl`.
For normal use, create a `LightSimulation` and call `run_light`.

## Main Workflow

File-based run:

```julia
sim, meteo = read_simulation("config.yml")
step = run_light(sim, first(meteo))
series = run_light(sim, meteo)
```

Interactive run:

```julia
using PlantGeom

scene = make_scene(domain=(-1.0, -1.0, 1.0, 1.0)) do s
    add_plant!(s, "plant.opf"; group="coffee", id=1)
    add_ground!(s; group="soil", type="ground")
end

models = models_for(
    "coffee" => ("Leaf" => translucent(par=0.15, nir=0.90),),
    "soil" => ("ground" => translucent(par=0.10, nir=0.40),),
)

sim = LightSimulation(scene, models; options=LightOptions())
step = run_light(sim, meteo_row)
```

For host-model coupling:

```julia
for row in rows
    light = run_light(sim, row)
end

update_scene!(sim, new_scene)
```

## Input Loading

Use these when your workflow starts from files or existing tables:

```julia
read_simulation(path)
read_scene(path)
read_models(path_or_paths)
read_options(path)
read_meteo(path_or_table)
```

## Scene And Model Helpers

Use these when inputs are built in Julia:

```julia
PlantGeom.make_scene(f; domain, source_path="interactive.scene", kwargs...)
PlantGeom.add_plant!(builder, mtg_or_path; group, id, at=(0, 0, 0), scale=1.0, rotate=(0, 0, 0))
PlantGeom.add_object!(builder, mtg_mesh_or_path; group, id, type="object", at=(0, 0, 0), scale=1.0, rotate=(0, 0, 0))
PlantGeom.add_ground!(builder; z=0.0, nx=9, ny=9, group="pavement", type="Cobblestone")
PlantGeom.prepare_scene(mtg; source_path="interactive.opf", scene_xy_bounds=nothing, relabel_ids=false)
models_for(group => (type => model, ...), ...)
translucent(; par, nir, transparency=0.0)
virtual_sensor()
emitter(; radiance, par=0.48, nir=0.52)
```

`PlantGeom.add_plant!` is the plant-named wrapper. `PlantGeom.add_object!` is the general placement
helper for MTGs, `GeometryBasics` meshes, `.opf`, and `.gwa` files. For mesh
files such as `.obj` or `.ply`, use MeshIO/FileIO to load the mesh first:

```julia
using FileIO, MeshIO
mesh = load("sensor.obj")
PlantGeom.add_object!(builder, mesh; group="sensor", type="panel", id=10)
```

Both placement helpers accept `at`, `scale`, `rotate`, `deg`, and OPS-style
`rotation`, `inclination_azimut`, and `inclination_angle` placement keywords.

Tuple rotations use fixed X, then Y, then Z order:

```julia
PlantGeom.add_object!(builder, mesh; group="sensor", type="panel", id=10, rotate=(10, 20, 30), deg=true)
```

Named-tuple rotations preserve the field order, which is useful when you need a
specific Euler sequence:

```julia
PlantGeom.add_object!(builder, mesh; group="sensor", type="panel", id=10, rotate=(y=20, z=30, x=10), deg=true)
```

`PlantGeom.add_ground!` and `write_scene` also work on prepared scenes:

```julia
PlantGeom.add_ground!(scene; z=0.0, nx=9, ny=9, xy_bounds=nothing, group="pavement", type="Cobblestone")
write_scene(path, scene)
```

## Validation

Use these to diagnose inputs before running:

```julia
check_scene(scene)
check_models(scene, models)
check_meteo(meteo; options=LightOptions())
check_simulation(sim)
check_simulation(scene, meteo; models, options=LightOptions())
summarize_scene(scene; models=nothing)
summarize_meteo(meteo; options=LightOptions())
```

Each function returns a `ValidationReport` with `errors`, `warnings`, and
`infos`.

The `summarize_*` helpers return structured summaries and print compact
diagnostics for humans. Use them when you are not sure what ArchimedLight sees:

```julia
summarize_scene(scene; models=models)
summarize_meteo(meteo; options=options)
```

## Simulation Cache

`LightSimulation` owns preparation and cache state. These helpers update inputs
and invalidate cached data:

```julia
update_scene!(sim, scene)
update_models!(sim, models)
update_options!(sim, options)
cache_summary(sim)
```

`update_scene!` immediately releases old scene-dependent prepared data and
cache entries.

## Advanced Light Pipeline

The explicit stage API remains available for debugging and research workflows:

```julia
ArchimedLight.compute_sky(row, options)
ArchimedLight.build_turtle(options, sky)
ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
ArchimedLight.compute_scattering(scene, models, turtle, first, options)
ArchimedLight.integrate_light(scene, models, first, scattering, options; meteo_row=row)
```

For interactive synthetic scenes, `run_light` also accepts a prebuilt sky state:

```julia
sim = LightSimulation(scene, models; options=options)
step = run_light(sim, sky; step_duration_seconds=1800.0)
```

The old low-level cache functions are still available for advanced work:

```julia
cache = ArchimedLight.prepare_light_cache(scene, models, options; ...)
ArchimedLight.run_light_step(cache, meteo_row)
ArchimedLight.run_light_series(cache, meteo)
```

`prepare_light_cache` uses a tiered policy internally:

- `:full` keeps all seen turtle responses in memory
- `:partial` keeps a bounded LRU cache for large moving-sun series
- `:topology_fallback` reuses prepared geometry/topology only when a full-response cache would be too large or unsupported

## Scene Attachment Helpers

These functions attach computed values back onto the MTG using ARCHIMED attribute names:

```julia
attach_node_values!(scene, attr, values; fill_value=nothing)
attach_light_step!(scene, step; fields=[:incident_par_flux], names=Dict(), fill_value=nothing)
attach_light_series!(scene, steps; fields=[:incident_par_flux], names=Dict(), fill_value=NaN)
```

`attach_light_step!` attaches scalar values for one step. `attach_light_series!`
attaches `Vector{Float64}` values, ordered like `steps`, for each selected
field. The supported `fields` selectors are:

| Selector | Default MTG attribute |
| --- | --- |
| `:area` | `area` |
| `:incident_par_initial_flux` | `Ri_PAR_0_f` |
| `:incident_nir_initial_flux` | `Ri_NIR_0_f` |
| `:incident_par_flux` | `Ri_PAR_f` |
| `:incident_nir_flux` | `Ri_NIR_f` |
| `:incident_par_initial_energy` | `Ri_PAR_0_q` |
| `:incident_nir_initial_energy` | `Ri_NIR_0_q` |
| `:incident_par_energy` | `Ri_PAR_q` |
| `:incident_nir_energy` | `Ri_NIR_q` |
| `:absorbed_par_initial_flux` | `Ra_PAR_0_f` |
| `:absorbed_nir_initial_flux` | `Ra_NIR_0_f` |
| `:absorbed_par_flux` | `Ra_PAR_f` |
| `:absorbed_nir_flux` | `Ra_NIR_f` |
| `:absorbed_par_initial_energy` | `Ra_PAR_0_q` |
| `:absorbed_nir_initial_energy` | `Ra_NIR_0_q` |
| `:absorbed_par_energy` | `Ra_PAR_q` |
| `:absorbed_nir_energy` | `Ra_NIR_q` |
| `:sky_fraction` | `sky_fraction` |

For `:area`, `attach_light_step!` attaches one scalar surface area per node,
while `attach_light_series!` repeats that area once per step so the attached
attribute has the same vector shape as the light fields.

Use `names=Dict(selector => attr)` to override default attribute names, for
example `Dict(:absorbed_nir_flux => :Ra_SW_f)`.

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
