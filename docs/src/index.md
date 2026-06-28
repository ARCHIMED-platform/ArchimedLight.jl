```@meta
CurrentModule = ArchimedLight
```

# ArchimedLight.jl

`ArchimedLight.jl` is the Julia reimplementation of the ARCHIMED light model. This model computes the interception and scattering of light by complex 3D scenes, with a focus on performance. It uses a rasterization-based approach for first-order interception and an iterative method for scattering, with flexible options for optical properties and directional response.

It is designed around a simple workflow: build or read a scene, define optical
models, create a `LightSimulation`, then call `run_light` for one meteo row or
a complete meteo table.

![Coffee scene light interception](assets/coffee_scene_light_interception.png)

The image above is generated from the bundled coffee example in `example_2/`. It loads the ARCHIMED coffee scene, runs one light step, attaches `Ri_PAR_f` onto the MTG, and renders the result with `PlantGeom.plantviz`.

## Scope

- ARCHIMED-style scene/model/meteo ingestion: `.ops`, `.opf`, `.gwa`, YAML models, and meteo CSV files
- File-based and in-memory workflows through the same runtime API
- Directional sky discretization with the standard ARCHIMED turtle sector counts
- First-order interception by CPU rasterization / z-buffer style projection
- Iterative scattering using the same directional visibility information
- Attachment of `Ri_*` and `Ra_*` values back onto MTG nodes for inspection and export
- Fixture-based parity work against the historical Java implementation

Energy balance, transpiration, and photosynthesis are intentionally out of scope for this package. The historical ARCHIMED documents referenced throughout this site still matter for the model vocabulary and physical assumptions, but the Julia package documented here currently implements the light-only core.

## Quick Start

```julia
using ArchimedLight

sim, meteo = read_simulation("config.yml")

step = run_light(sim, first(meteo))

step.budget.incident_flux.total.par
step.budget.absorbed_energy.total.par
```

The simulation results are grouped by quantity and waveband in `LightBudget`. When you attach those values back onto the scene, the default attribute names keep the standard ARCHIMED naming convention such as `Ri_PAR_f`, `Ri_PAR_q`, and `Ra_PAR_q`.

## Read This Site In Three Passes

- Start with [Getting Started](getting_started.md) if you want one runnable coffee example with the minimum number of moving parts.
- Continue with [Beginner Workflows](tutorial_beginner_workflows.md), [File-Based Workflow](tutorial_files.md), or [Interactive Workflow](tutorial_interactive.md) depending on whether your scene already exists on disk or is being built in Julia.
- Use the reference pages for exact file keys, scene semantics, model structure, meteo columns, and outputs.

## Documentation Map

- [Getting Started](getting_started.md)
- [Beginner Workflows](tutorial_beginner_workflows.md)
- [Tutorial: File-Based Workflow](tutorial_files.md)
- [Tutorial: Interactive Workflow](tutorial_interactive.md)
- [Configuration And Options Reference](reference_config.md)
- [Scene Files And Semantics](reference_scene.md)
- [Model Files Reference](reference_models.md)
- [Meteo Inputs Reference](reference_meteo.md)
- [Outputs](outputs.md)
- [Pipeline Overview](theory_pipeline.md)
- [First-Order Interception](theory_interception.md)
- [Scattering And Optical Assumptions](theory_scattering.md)
- [Full Example](full_example.md)
- [Composable Stages](stages.md)
- [GPU Backends And Benchmarks](gpu_backends.md)
- [API Reference](api.md)

## Contributors

```@raw html
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
```
