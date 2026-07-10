# API Reference

This page is a compact guide to the public entry points of `ArchimedLight.jl`.
For normal use, create a `LightSimulation` and call `run_light`.

## API Stability

Version `0.1.3` is the compatibility baseline for ArchimedLight's first release
in Julia's General registry. Earlier `0.1.x` tags were development snapshots
and are not covered by this stability promise. Starting with `0.1.3`, the
supported public API is the set of names exported by `ArchimedLight`. This
includes the high-level simulation workflow, input readers, model helpers,
validation helpers, user-facing result containers, attachment/output helpers,
visualization helpers, and concrete backend selectors. Patch releases after
`0.1.3` are expected to preserve those names; incompatible changes will use a
new minor version and include migration notes.

Some lower-level stage functions are available as qualified calls such as
`ArchimedLight.compute_sky(...)`. These are intended for debugging, research,
and parity work. They are not exported, and their exact arguments, return
containers, and cache internals may be refined before `1.0`.

`LightStepResult` and `LightBudget` are public result containers. It is fine to
read their fields in analysis code, but user code should normally construct
results through `run_light` rather than calling result constructors directly.

`LightSimulation` owns internal preparation and cache state. Use
`update_scene!`, `update_models!`, `update_options!`, and `cache_summary` rather
than relying on the `cache` field or cache object layout.

## Main Workflow

File-based run:

```@example api_file_workflow
using ArchimedLight

repo_root = normpath(joinpath(dirname(pathof(ArchimedLight)), ".."))
config = joinpath(repo_root, "example_2", "config.yml")
sim, meteo = read_simulation(config)
step = run_light(sim, first(meteo))
series = run_light(sim, meteo);
```

Interactive run. This is a schematic API shape because it assumes that the
host application provides `plant.opf` and `meteo_row`; see the checked
Interactive Workflow page for a fully executable in-memory scene.

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

For host-model coupling, this schematic loop assumes your host model provides
the meteo `rows` and a later `new_scene`:

```julia
for row in rows
    light = run_light(sim, row)
end

update_scene!(sim, new_scene)
```

## Input Loading

Use these when your workflow starts from files or existing tables. The names in
this block are schematic placeholders for your actual paths and tables.

```julia
read_simulation(path)
read_scene(path)
read_models(path_or_paths)
read_options(path)
read_meteo(path_or_table)
```

## Scene And Model Helpers

Use these when inputs are built in Julia. This block lists call signatures and
is intentionally schematic.

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
files such as `.obj` or `.ply`, use MeshIO/FileIO to load the mesh first. This
snippet is schematic because it depends on your local mesh file and builder:

```julia
using FileIO, MeshIO
mesh = load("sensor.obj")
PlantGeom.add_object!(builder, mesh; group="sensor", type="panel", id=10)
```

Both placement helpers accept `at`, `scale`, `rotate`, `deg`, and OPS-style
`rotation`, `inclination_azimut`, and `inclination_angle` placement keywords.

Tuple rotations use fixed X, then Y, then Z order. This is a schematic call
shape:

```julia
PlantGeom.add_object!(builder, mesh; group="sensor", type="panel", id=10, rotate=(10, 20, 30), deg=true)
```

Named-tuple rotations preserve the field order, which is useful when you need a
specific Euler sequence. This is a schematic call shape:

```julia
PlantGeom.add_object!(builder, mesh; group="sensor", type="panel", id=10, rotate=(y=20, z=30, x=10), deg=true)
```

`PlantGeom.add_ground!` and `write_scene` also work on prepared scenes. This
snippet is schematic because it assumes an existing `scene` and output `path`:

```julia
PlantGeom.add_ground!(scene; z=0.0, nx=9, ny=9, xy_bounds=nothing, group="pavement", type="Cobblestone")
write_scene(path, scene)
```

## Validation

Use these to diagnose inputs before running. This block is schematic and uses
placeholder inputs already introduced above:

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
`infos`. Meteo-aware checks and execution methods accept
`check_boundaries=false` as a one-call override; derivability, duration, and
finite-input consistency checks remain enabled.

The `summarize_*` helpers return structured summaries and print compact
diagnostics for humans. Use them when you are not sure what ArchimedLight sees.
This schematic snippet assumes `scene`, `models`, `meteo`, and `options` exist:

```julia
summarize_scene(scene; models=models)
summarize_meteo(meteo; options=options)
```

## Simulation Cache

`LightSimulation` owns preparation and cache state. These helpers update inputs
and invalidate cached data. This block is schematic and assumes existing
replacement inputs:

```julia
update_scene!(sim, scene)
update_models!(sim, models)
update_options!(sim, options)
cache_summary(sim)
```

`update_scene!` immediately releases old scene-dependent prepared data and
cache entries. The cache object itself is an implementation detail; use
`cache_summary(sim)` when you need to inspect cache behavior.

## Advanced Light Pipeline

The explicit stage API is available as qualified, advanced API for debugging,
research, and parity workflows. These calls may return stage containers such as
`ArchimedLight.TurtleGrid` or `ArchimedLight.FirstOrderResult`, which are not
exported in `0.1.x`. This block is schematic because it assumes all previous
stage inputs have already been prepared:

```julia
ArchimedLight.compute_sky(row, options)
ArchimedLight.build_turtle(options, sky)
ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
ArchimedLight.compute_scattering(scene, models, turtle, first, options)
ArchimedLight.integrate_light(scene, models, first, scattering, options; meteo_row=row)
```

These functions and stage containers are intentionally not exported in `0.1.x`.
Prefer `run_light` for application code unless you need to inspect or replace
individual pipeline stages.

For interactive synthetic scenes, `run_light` also accepts a prebuilt sky state.
This block is schematic because it assumes an existing `scene`, `models`,
`options`, and `sky`:

```julia
sim = LightSimulation(scene, models; options=options)
step = run_light(sim, sky; step_duration_seconds=1800.0)
```

The low-level cache functions are still available for advanced work, but their
cache object layout is internal. This block is schematic and uses placeholders:

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

These functions attach computed values back onto the MTG using ARCHIMED
attribute names. This block lists schematic call shapes:

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

These helpers expose direct mesh coloring and the Makie package extension. This
block lists schematic call shapes:

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

This block lists schematic constructor calls:

```julia
RasterCPUBackend()
RaycastScatteringBackend()
```

They can be passed explicitly to the pipeline helpers when you want to control
the implementation used for interception or scattering. Only the concrete
selectors are exported in `0.1.x`; the backend supertypes are internal, and
defining new backend subtypes is not yet a supported extension interface.

## Recommended Starting Points

- read [Getting Started](getting_started.md) for the shortest runnable workflow
- read [Composable Stages](stages.md) if you want the stage-by-stage pattern
- read [Outputs](outputs.md) for the mapping between `LightBudget` fields and ARCHIMED attribute names


## API List

```@autodocs
Modules = [ArchimedLight]
Private = false
```
