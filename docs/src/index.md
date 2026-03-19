```@meta
CurrentModule = ArchimedLight
```

# ArchimedLight.jl

Julia reimplementation of the ARCHIMED light model.

![Coffee scene light interception](assets/coffee_scene_light_interception.png)

The figure above is generated from the bundled coffee fixture with `scripts/generate_home_figure.jl`. The script adds explicit ground geometry, runs one light step, attaches `Ri_PAR_f` onto the MTG, and then renders the scene with `plantviz(..., color=:Ri_PAR_f)`.

## Scope

- Scene/model/meteo ingestion (`.ops`, `.opf`, `.gwa` + YAML + meteo files)
- Sky and turtle discretization
- First-order interception by CPU raster/z-buffer
- Iterative scattering
- Fast manual fixture tests (`test/fast_fixtures`) + release-only heavy artifact regression

Energy balance, transpiration, and photosynthesis are meant to be computed by [PlantBiophysics.jl](https://github.com/VEZY/PlantBiophysics.jl).

## Quick Start

```julia
using ArchimedLight

options, scene, meteo, models = read_config("config.yml")

row = first(prepare_meteo(meteo, options).rows)
step = run_light_step(scene, models, row, options)

step.budget.incident_flux.total.par
step.budget.absorbed_energy.total.par
```

The simulation results are grouped by quantity and waveband. Exported files and attached MTG attributes keep the standard ARCHIMED names (Java version) such as `Ri_PAR_f` and `Ra_PAR_q`.

## Navigation

- See [ARCHIMED Reference](archimed_reference.md) for a detailed explanation of the original Java algorithm, the Julia port, and the modeling assumptions.
- See [Full Example](full_example.md) for a complete runnable workflow with bundled input files.
- See [Composable Stages](@ref) for stage-by-stage usage.
- See [API Reference](@ref) for public signatures.
