```@meta
CurrentModule = ArchimedLight
```

# ArchimedLight.jl

Julia reimplementation of the ARCHIMED light model with a function-first API.

## Scope
- Scene/config/meteo ingestion (`.ops`, `.opf`, `.gwa` + YAML + meteo files)
- Sky and turtle discretization
- First-order interception by CPU raster/z-buffer
- Iterative scattering
- Fast manual fixture tests (`test/fast_fixtures`) + release-only heavy artifact regression

Energy balance, transpiration, and photosynthesis are intentionally out of scope.

## Quick Start
```julia
using ArchimedLight

cfg = read_light_config("config.yml")
scene = read_scene(cfg.scene)
meteo = read_meteo(cfg.meteo)

row = first(meteo.rows)
step = run_light_step(scene, row, cfg)
```

## Navigation
- See [ARCHIMED Reference](archimed_reference.md) for a detailed explanation of the original Java algorithm, the Julia port, and the modeling assumptions.
- See [Full Example](full_example.md) for a complete runnable workflow with bundled input files.
- See [Composable Stages](@ref) for stage-by-stage usage.
- See [API Reference](@ref) for public signatures.
