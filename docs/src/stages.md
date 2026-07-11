# Composable Stages

`ArchimedLight.jl` is designed so each stage can be called independently for
advanced debugging and parity work. These staged functions are not exported by
the stable public API in the `0.1.x` series; call them through the module name
and prefer `run_light` for application code.

## Inputs

```julia
using ArchimedLight

sim, meteo = read_simulation("config.yml")
scene = sim.scene
models = sim.models
options = sim.options
row = first(meteo)
```

## Stage-by-stage Execution

```julia
sky = ArchimedLight.compute_sky(row, options)
turtle = ArchimedLight.build_turtle(options, sky)
fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
scat = ArchimedLight.compute_scattering(scene, models, turtle, first, options)
budget = ArchimedLight.integrate_light(scene, models, first, scat, options; meteo_row=row)

budget.incident_flux.total.par
budget.incident_energy.total.par
budget.absorbed_flux.total.par
```

In Julia code, light results are accessed through grouped `LightBudget` fields. Exported files and attached visualization attributes keep ARCHIMED names such as `Ri_PAR_f`, `Ri_PAR_q`, and `Ra_PAR_q`.

## Model Coefficients And Option Fallbacks

Model-defined optical properties are used first. When a model omits one band in
`optical_properties`, the runtime falls back to the matching global option:

- missing `PAR` coefficient -> `options.scattering_coeff_par`
- missing `NIR` coefficient -> `options.scattering_coeff_nir`

For example, with `LightOptions()`, a missing model `NIR` coefficient uses `0.30`, and the
default absorbed fraction used in the final budget becomes `1 - 0.30 = 0.70`.

Virtual sensors are declared on the interception model (`sensor=true`, or legacy
`model: VirtualSensor` in YAML). They receive light, but remain transparent and non-absorbing.

## Single-Call Pipeline

```julia
step = run_light(sim, row)
series = run_light(sim, meteo)

attach_light_step!(scene, step; fields=[:incident_par_flux])
write_scene("output/scene.opf", scene)

# Optional backend kwargs:
step = run_light(
    LightSimulation(
        scene,
        models;
        options=options,
        interception_backend=RasterCPUBackend(),
        scattering_backend=RaycastScatteringBackend(),
    ),
    row,
)
```

## Caching Options

- `LightOptions(cache_radiation=true)` reuses directional responses across meteo rows in `run_light`.
- `LightOptions(cache_pixel_table=true)` stores per-direction projection tables in the on-disk cache directory used by the interception stage.

## Backends and Modes

- `ArchimedLight.compute_first_order(...; backend=:raster_cpu)` is the reference backend (`RasterCPUBackend()` also works).
- `ArchimedLight.compute_scattering(...; mode=:raycast)` is the current scattering mode.
- `ArchimedLight.compute_scattering(...; backend=RaycastScatteringBackend())` is also available.
- `RasterCPUBackend()` and `RaycastScatteringBackend()` are public concrete selectors in `0.1.x`, not yet an extension interface for user-defined backend subtypes.
- `PlantGeom.add_ground!(scene; ...)` is the explicit scene-editing step when you want inspectable paving or soil geometry in the MTG.
