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
- `pixel_size` is validated with Java parity bounds (`0 < pixel_size <= 0.5` meters).
- `cache_pixel_table: true` (or `save_on_disk: true`) enables on-disk direction projection cache under `<output_directory>/pixel_tables_cache`.
- `build_turtle` follows Java-compatible sector sets for `1, 6, 16, 46, 136, 406`.

## Testing
Run the parity and smoke suite:

```bash
julia --project=. test/runtests.jl
```

Run only the frozen-reference parity suite:

```bash
ARCHIMEDLIGHT_TEST_PROFILE=parity julia --project=. test/runtests.jl
```

Run the fast parity guard subset (hitcount + weighted sun + one scattering fixture):

```bash
ARCHIMEDLIGHT_TEST_PROFILE=guard julia --project=. test/runtests.jl
```

Run only the isolated source-built Java sky-matrix parity suite:

```bash
ARCHIMEDLIGHT_TEST_PROFILE=sky_matrix julia --project=. test/runtests.jl
```

Run only the explicit synthetic simple-scene cases:

```bash
ARCHIMEDLIGHT_TEST_PROFILE=synthetic julia --project=. test/runtests.jl
```

Run one named synthetic case only:

```bash
ARCHIMEDLIGHT_TEST_PROFILE=synthetic ARCHIMEDLIGHT_SYNTHETIC_CASE=two_planes_shadow_absorptance julia --project=. test/runtests.jl
```

The dedicated synthetic cases are defined in `test/synthetic_scene_cases.jl` with explicit
`inputs` and `expected` blocks. Current case names include:
`single_plate_direct`, `stacked_scattering`, `partial_overlap_direct`,
`tilted_plate_projection`, `oblique_shadow`, `toricity_wraparound`,
`toricity_cross_border_shadow`, `toricity_diffuse_cross_border_shadow`,
`toricity_scattering_cross_border`, `virtual_sensor_transparency`, `single_plate_absorptance`,
`two_planes_shadow_absorptance`, and `cached_series_parity`.

Additional explicit fixture-based simple-plant unit checks are part of the `core` profile:
`test-scattering-two-simpleplants` and `test-toricity-two-simpleplants-border`.

Build the upstream Java jar used to freeze source-driven references:

```bash
cd /Users/rvezy/Documents/dev/ARCHIMED/archimed-2018-source
mvn -q -am -pl :archimed-lib-2018 package -Dmaven.test.skip=true
```

Freeze and sync the upstream Java sky-matrix fixtures into this repo:

```bash
julia --project=. scripts/freeze_java_light_refs.jl \
  --jar /Users/rvezy/Documents/dev/ARCHIMED/archimed-2018-source/archimed-lib-2018/target/archimed-lib-2018-0.0.1-SNAPSHOT-jar-with-dependencies.jar \
  --upstream-tests-root /Users/rvezy/Documents/dev/ARCHIMED/archimed-2018-source/archimed-lib-2018/tests \
  --fixtures test-compare-sky06,test-compare-sky16,test-compare-sky46 \
  --force
```

The upstream Java source repo is the source of truth for these fixture definitions; this repo stores the mirrored fixture files and frozen `expected/` references consumed by Julia.

`scripts/freeze_java_light_refs.jl` launches Java with `-Djava.awt.headless=true` by default to avoid macOS AWT/AppKit hangs when running the x86_64 JDK8 binary under Rosetta.
