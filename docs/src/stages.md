# Composable Stages

`ArchimedLight.jl` is designed so each stage can be called independently.

## Inputs

```julia
using ArchimedLight

scene = read_scene("scene.ops")
models = read_models(["plant.yml", "soil.yml"])
options = LightOptions()
meteo = read_meteo("meteo.csv")
row = first(prepare_meteo(meteo, options).rows)
```

## Stage-by-stage Execution

```julia
sky = compute_sky(row, options)
turtle = build_turtle(options, sky)
fluxes = compute_directional_fluxes(row, sky, turtle, options)
first = compute_first_order(scene, models, turtle, fluxes, options)
scat = compute_scattering(scene, models, turtle, first, options)
budget = integrate_light(scene, models, first, scat, options; meteo_row=row)

budget.incident_flux.total.par
budget.incident_energy.total.par
budget.absorbed_flux.total.par
```

In Julia code, light results are accessed through grouped `LightBudget` fields. Exported files and attached visualization attributes keep ARCHIMED names such as `Ri_PAR_f`, `Ri_PAR_q`, and `Ra_PAR_q`.

## Single-Call Pipeline

```julia
step = run_light_step(scene, models, row, options)
series = run_light_series(scene, models, meteo, options)

attach_light_step!(scene, step; fields=[:incident_par_flux])
write_scene("output/scene.opf", scene)

# Optional backend kwargs:
step = run_light_step(
    scene,
    models,
    row,
    options;
    interception_backend=RasterCPUBackend(),
    scattering_backend=RaycastScatteringBackend(),
)
```

## Caching Options

- `LightOptions(cache_radiation=true)` reuses directional responses across meteo rows in `run_light_series`.
- `LightOptions(cache_pixel_table=true)` stores per-direction projection tables in the on-disk cache directory used by the interception stage.

## Backends and Modes

- `compute_first_order(...; backend=:raster_cpu)` is the reference backend (`RasterCPUBackend()` also works).
- `compute_scattering(...; mode=:raycast)` and `compute_scattering(...; mode=:links)` are available.
- `compute_scattering(...; backend=RaycastScatteringBackend())` and `compute_scattering(...; backend=LinksScatteringBackend())` are also available.
- `add_ground!(scene; ...)` is the explicit scene-editing step when you want inspectable paving or soil geometry in the MTG.
