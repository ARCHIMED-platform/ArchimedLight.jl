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
scat = compute_scattering(scene, turtle, first_order, cfg)
budget = integrate_light(first_order, scat, cfg)
```

## Single-call pipeline
```julia
step = run_light_step(scene, first(meteo.rows), cfg)
series = run_light_series(scene, meteo, cfg)
```

## Stage flexibility
- You can call each stage independently (`compute_sky`, `build_turtle`, `compute_first_order`, `compute_scattering`, ...).
- `compute_sky` follows Java clearness/global conversion and DeJong hourly direct/diffuse partitioning.
- `compute_sky` uses Java-style substep-weighted sun position (`radiation_timestep`) when sun angles are not provided.
- `compute_first_order(...; backend=:raster_cpu)` is the current reference backend.
- `pixel_size` is validated with Java parity bounds (`0 < pixel_size <= 0.5` meters).
- `compute_scattering(...; mode=:raycast)` and `compute_scattering(...; mode=:links)` are both available.
- `build_turtle` follows Java-compatible sector sets for `1, 6, 16, 46, 136, 406`.

## Testing
Run the parity and smoke suite:

```bash
julia --project=. test/runtests.jl
```
