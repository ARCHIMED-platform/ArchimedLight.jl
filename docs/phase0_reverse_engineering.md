# Phase 0: Java Reverse-Engineering Feasibility Note

## Scope
This note covers the Java light pipeline only (no energy balance, transpiration, photosynthesis).

## Java to Julia Mapping
- Java `MirProcess` -> Julia `run_light_step`, `run_light_series`
- Java `Turtle`, `TurtleFluxes`, `RadiativeFluxes` -> Julia `build_turtle`, `compute_directional_fluxes`
- Java `Sky`/`Sun` handling -> Julia `compute_sky`
- Java `MeshProjection` + pixel tables -> Julia `compute_first_order` (`:raster_cpu` backend)
- Java `Musc`/`EnergyTransferTask` -> Julia `compute_scattering`
- Java scene/config loaders -> Julia `read_scene`, `read_light_config`, `read_meteo`

## Feasibility Assessment
- Feasible in Julia with a composable function pipeline.
- No technical blocker identified for CPU implementation.
- Main data blocker is `.gwa` support in PlantGeom; this is required for many Java fixtures.
- Numeric parity risk is high in:
  - raster discretization details (pixel alignment and depth tie-breaking),
  - sun/diffuse decomposition details,
  - scattering transfer model.

## Feature Matrix
### Required for parity now
- OPS/OPF/GWA scene loading
- YAML config loading (light-related options)
- Meteo loading
- Sky decomposition
- Turtle sectors and directional flux assignment
- First-order directional visibility/interception with z-buffer rasterization
- Multi-step scattering with convergence/divergence checks
- Cached-radiation mode

### Out of scope now
- Energy balance
- Transpiration
- Photosynthesis
- GPU execution

## Blockers and Mitigations
- Blocker: PlantGeom OPS parser accepts `.opf` only.
  - Mitigation: add `read_gwa` and allow `.gwa` in OPS parsing.
- Blocker: Java parity baselines are large and heterogeneous.
  - Mitigation: start with curated light-only fixture list and tighten tolerances incrementally.

## Go / No-Go
Go, with the following conditions:
- Land PlantGeom GWA support first.
- Keep API modular and pure-function.
- Treat Java numeric parity as acceptance gate for curated fixtures.
