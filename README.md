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
