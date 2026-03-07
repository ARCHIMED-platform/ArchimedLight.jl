# ArchimedLight.jl

Julia reimplementation of the ARCHIMED light interception pipeline with a composable, function-first API.

## Current scope
- Scene/config/meteo input pipeline
- Sky + turtle discretization
- First-order interception (CPU raster/z-buffer)
- Iterative scattering (CPU reference)
- Julia-native fixture regression harness (numeric + visual references)

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
graph = build_scattering_transfer_graph(scene, turtle, first_order, cfg)
scat = compute_scattering(graph, first_order, cfg)
step_seconds = 1800.0 # use your meteo timestep duration in seconds
budget = integrate_light(first_order, scat, cfg; step_duration_seconds=step_seconds, component_area_per_node=scene.total_area_per_node)
```

`LightBudget` now includes Java-style PAR/NIR intercepted and absorbed outputs in both:
- `*_f`: irradiance (`W m^-2`)
- `*_q`: per-component energy per timestep (`J`)

## Single-call pipeline
```julia
step = run_light_step(scene, first(meteo.rows), cfg)
series = run_light_series(scene, meteo, cfg)
# Optional backend kwargs:
step = run_light_step(scene, first(meteo.rows), cfg; interception_backend=RasterCPUBackend(), scattering_backend=RaycastScatteringBackend())

# Java-style component_values.csv export:
write_component_values_csv("output/component_values.csv", scene, step, cfg; meteo_row=first(meteo.rows), step_number=0)

# Java-style scene_values.csv export:
write_scene_values_csv("output/scene_values.csv", scene, series, cfg; meteo_rows=meteo.rows)

# Java-style logs:
write_sun_position_log_csv("output/log-sun-position.csv", series, meteo.rows)
write_scattering_iteration_log_csv("output/log-iteration-scat-par.csv", scene, step, cfg; meteo_row=first(meteo.rows), band="PAR")

# High-level writer (component/scene/summary/logs from config defaults, overridable):
write_light_outputs(scene, series, cfg; meteo_rows=meteo.rows, outdir="output")

# Java-style simulation output path:
sim_out = simulation_output_directory(cfg)  # e.g. <output_directory>/000001
```

## Full Example
- Self-contained files and script are under `example/`.
- Run with:

```bash
julia --project=. example/full_featured_example.jl
```

## Stage flexibility
- You can call each stage independently (`compute_sky`, `build_turtle`, `compute_first_order`, `compute_scattering`, ...).
- You can prebuild scattering transfers via `build_scattering_transfer_graph(...)` and reuse them with `compute_scattering(graph, ...)`.
- `compute_sky` follows Java clearness/global conversion and DeJong hourly direct/diffuse partitioning.
- `compute_sky` uses Java-style substep-weighted sun position (`radiation_timestep`) when sun angles are not provided.
- Meteo `#' use: ...` consistency checks for `clearness`/`RI_SW_f`/`RI_PAR_f`/`RI_NIR_f` are enforced like Java.
- `compute_first_order(...; backend=:raster_cpu)` is the current reference backend (`RasterCPUBackend()` also available).
- `compute_scattering(...; mode=:raycast)` / `compute_scattering(...; mode=:links)` keep Java-style mode selection; backend objects are also available (`RaycastScatteringBackend()`, `LinksScatteringBackend()`).
- Component output variables are validated for light-only scope: scattering outputs require `scattering: true`; photosynthesis/energy-balance/TIR outputs are intentionally rejected (compute later with PlantBiophysics).
- `pixel_size` is validated with ARCHIMED-compatible bounds (`0 < pixel_size <= 0.5` meters).
- `cache_pixel_table: true` (or `save_on_disk: true`) enables on-disk direction projection cache under `<output_directory>/pixel_tables_cache`.
- `build_turtle` follows Java-compatible sector sets for `1, 6, 16, 46, 136, 406`.

## Testing
Run the default fast suite:

```bash
julia --project=. test/runtests.jl
```

This runs core checks, fast manual fixtures, and synthetic scene unit tests in one command.
`ARCHIMEDLIGHT_TEST_PROFILE` is reserved for `release` only.

Run one named synthetic case only:

```bash
ARCHIMEDLIGHT_SYNTHETIC_CASE=two_planes_shadow_absorptance julia --project=. test/runtests.jl
```

Run one named fast fixture case only:

```bash
ARCHIMEDLIGHT_FAST_FIXTURE_CASE=simpleplant_16_toric julia --project=. test/runtests.jl
```

The dedicated synthetic cases are defined in `test/synthetic_scene_cases.jl` with explicit
`inputs` and `expected` blocks. Current case names include:
`single_plate_direct`, `stacked_scattering`, `partial_overlap_direct`,
`tilted_plate_projection`, `oblique_shadow`, `toricity_wraparound`,
`toricity_cross_border_shadow`, `toricity_diffuse_cross_border_shadow`,
`toricity_scattering_cross_border`, `virtual_sensor_transparency`, `single_plate_absorptance`,
`two_planes_shadow_absorptance`, and `cached_series_parity`.

Fast fixture inputs/references are under `test/fast_fixtures/` and are intended to be readable
as usage examples.

Generate or refresh fast fixture references (simple-plant numeric CSV + image):

```bash
julia --project=. scripts/generate_fast_fixture_references.jl
```

## Release-only heavy regression (artifact)

Heavy full-fixture regression can be kept outside the package repository and run only before
releases.

Build an artifact tarball from a heavy fixture dataset checkout:

```bash
julia --project=. scripts/build_release_fixture_artifact.jl \
  --test-root /path/to/heavy-test-root \
  --tarball /tmp/archimedlight-release-fixtures.tar.gz
```

Optional: bind the artifact in `Artifacts.toml` by providing the download URL:

```bash
julia --project=. scripts/build_release_fixture_artifact.jl \
  --test-root /path/to/heavy-test-root \
  --tarball /tmp/archimedlight-release-fixtures.tar.gz \
  --url https://github.com/ARCHIMED-platform/ArchimedLight.jl/releases/download/v0.0.1/archimedlight-release-fixtures-v0.0.1.tar.gz
```

Run release-only heavy regression:

```bash
ARCHIMEDLIGHT_TEST_PROFILE=release julia --project=. test/runtests.jl
```

You can also bypass artifacts and point directly to a local extracted release dataset:

```bash
ARCHIMEDLIGHT_TEST_PROFILE=release \
ARCHIMEDLIGHT_RELEASE_FIXTURES_DIR=/path/to/release-fixtures \
julia --project=. test/runtests.jl
```

Optional (release dataset fixture filter):

```bash
ARCHIMEDLIGHT_TEST_PROFILE=release \
ARCHIMEDLIGHT_RELEASE_FIXTURES_DIR=/path/to/release-fixtures \
ARCHIMEDLIGHT_FIXTURE_FILTER=test-compare-simpleplant \
julia --project=. test/runtests.jl
```
