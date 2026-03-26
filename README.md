# ArchimedLight.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://VEZY.github.io/ArchimedLight.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://VEZY.github.io/ArchimedLight.jl/dev)
[![Test workflow status](https://github.com/VEZY/ArchimedLight.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/VEZY/ArchimedLight.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/VEZY/ArchimedLight.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/VEZY/ArchimedLight.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![All Contributors](https://img.shields.io/github/all-contributors/VEZY/ArchimedLight.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

Julia reimplementation of the ARCHIMED light interception pipeline with a composable, function-first API.

![Coffee scene light interception](docs/src/assets/coffee_scene_light_interception.png)

The figure above is generated from the bundled coffee fixture with `scripts/generate_home_figure.jl`. The script loads a scene, models, options, and meteo rows, adds explicit ground geometry, runs one light step, attaches `Ri_PAR_f` onto the MTG, and then renders the colored scene with `plantviz(..., color=:Ri_PAR_f)`.

## Current scope
- Scene/model/meteo input pipeline
- Sky + turtle discretization
- First-order interception (CPU raster/z-buffer)
- Iterative scattering (CPU reference)
- Julia-native fixture regression harness (numeric + visual references)

Energy balance, transpiration and photosynthesis are intentionally out of scope for now.

## Core API
```julia
using ArchimedLight

options, scene, meteo, models = read_config("config.yml")

row = first(prepare_meteo(meteo, options).rows)
sky = compute_sky(row, options)
turtle = build_turtle(options, sky)
fluxes = compute_directional_fluxes(row, sky, turtle, options)
first_order = compute_first_order(scene, models, turtle, fluxes, options)
scat = compute_scattering(scene, models, turtle, first_order, options)
budget = integrate_light(scene, models, first_order, scat, options; meteo_row=row)
```

In Julia code, `LightBudget` is grouped by quantity and waveband:

```julia
budget.incident_flux.total.par
budget.incident_energy.total.par
budget.absorbed_flux.total.nir
budget.absorbed_energy.initial.par
```

File exports and attached MTG attributes keep the ARCHIMED names:
- `Ri_*`: incident light
- `Ra_*`: absorbed light
- `*_f`: irradiance (`W m^-2`)
- `*_q`: per-component energy per timestep (`J`)

Band coefficients come from the model definition when present. If a model omits one band in
`optical_properties`, the runtime falls back to the global option for that band:
- `PAR` fallback: `LightOptions(scattering_coeff_par=...)`
- `NIR` fallback: `LightOptions(scattering_coeff_nir=...)`

With the default options, that means:
- missing model `PAR` coefficient falls back to `0.15`
- missing model `NIR` coefficient falls back to `0.30`
- the corresponding default absorptances used in the final budget are `1 - coeff`

## Short pipeline
```julia
options, scene, meteo, models = read_config("config.yml")
row = first(prepare_meteo(meteo, options).rows)

step = run_light_step(scene, models, row, options)
series = run_light_series(scene, models, meteo, options)

attach_light_step!(scene, step; fields=[:incident_par_flux, :absorbed_par_total_energy])
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

## Full Example
- Self-contained files and script are under `example_1/`.
- Run with:

```bash
julia --project=. example_1/full_featured_example.jl
```

## Stage flexibility
- `read_config("config.yml")` is the convenience entrypoint for file-driven workflows and returns `options, scene, meteo, models`.
- You can call each stage independently (`compute_sky`, `build_turtle`, `compute_first_order`, `compute_scattering`, `integrate_light`, ...).
- File-based and in-memory workflows share the same runtime path: `read_scene(...)` and `prepare_scene(...)` both produce `SceneGeometry`, while `read_models(...)` and `prepare_models(...)` both produce `LightModels`.
- `compute_sky` follows the ARCHIMED clearness/global conversion and DeJong hourly direct/diffuse partitioning.
- `compute_sky` uses substep-weighted sun position (`radiation_timestep`) when sun angles are not provided.
- Meteo `#' use: ...` consistency checks for `clearness`/`RI_SW_f`/`RI_PAR_f`/`RI_NIR_f` are enforced like Java.
- `compute_first_order(...; backend=:raster_cpu)` is the current reference backend (`RasterCPUBackend()` also available).
- `compute_scattering(...; mode=:raycast)` / `compute_scattering(...; mode=:links)` expose the two supported scattering modes; backend objects are also available (`RaycastScatteringBackend()`, `LinksScatteringBackend()`).
- `pixel_size` is validated with ARCHIMED-compatible bounds (`0 < pixel_size <= 0.5` meters).
- `cache_pixel_table=true` enables an on-disk direction projection cache under the interception cache directory.
- `build_turtle` follows the canonical ARCHIMED sector counts `1, 6, 16, 46, 136, 406`.
- `add_ground!(scene; ...)` is an explicit scene-editing step, so inspectable ground is ordinary geometry in the exported MTG.
- Virtual sensors are declared on the interception model (`sensor=true`, or legacy `model: VirtualSensor` in YAML). They receive light but stay transparent and non-absorbing in the simulation.

## Testing
Run the default fast suite:

```bash
julia --project=test test/runtests.jl
```

This runs the standalone `@testitem`s discovered from the `*-test.jl` files under `test/`,
excluding the `:release`-tagged heavy regression item.

Run only tests with the `:fast` tag directly:

```bash
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests("test"; filter=ti -> :fast in ti.tags, verbose=true)'
```

Run one tagged synthetic case directly:

```bash
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests("test"; filter=ti -> :single_plate_direct in ti.tags, verbose=true)'
```

Run one tagged fast fixture directly:

```bash
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests("test"; filter=ti -> :simpleplant_16_toric in ti.tags, verbose=true)'
```

Run the opt-in regression matrix against frozen Julia baselines:

```bash
julia --project=test test/regression_matrix/runtests.jl
```

Useful regression-matrix controls:

```bash
# refresh in-repo baselines intentionally
ARCHIMEDLIGHT_REGRESSION_UPDATE=true julia --project=test test/regression_matrix/runtests.jl

# write reports somewhere else
ARCHIMEDLIGHT_REGRESSION_REPORT_DIR=/tmp/archimedlight-regression \
julia --project=test test/regression_matrix/runtests.jl
```

`ARCHIMEDLIGHT_REGRESSION_FILTER` accepts a comma-separated list of exact case ids from
`regression_report.csv`.

The regression harness writes a machine-readable CSV report at
`test/regression_matrix/reports/latest/regression_report.csv` by default and stores frozen
fast-profile baselines under `test/regression_matrix/baselines/`.

The dedicated synthetic cases are defined in `test/synthetic-scenes-test.jl`. Current case names include:
`single_plate_direct`, `stacked_scattering`, `toricity_wraparound`,
`virtual_sensor_transparency`, `run_light_step_matches_staged`,
`cache_radiation_parity`, `cached_scattering_series_parity`, and `missing_models`.

Fast fixture inputs/references are under `test/fast_fixtures/` and are intended to be readable
as usage examples.

Generate or refresh fast fixture references (simple-plant numeric CSV + image):

```bash
julia --project=. scripts/generate_fast_fixture_references.jl
```

## Release-only heavy regression (artifact)

Heavy full-fixture regression can be kept outside the package repository and run only before
releases.

Build a data-only artifact tarball from a heavy fixture dataset checkout (fixtures + references +
reference images + manifest):

```bash
julia --project=. scripts/build_release_fixture_artifact.jl \
  --test-root /path/to/heavy-test-root \
  --tarball /tmp/archimedlight-release-fixtures.tar.gz
```

If you want the current Julia outputs to become the packaged references first, refresh and package
in one step:

```bash
julia --project=. scripts/build_release_fixture_artifact.jl \
  --test-root /path/to/heavy-test-root \
  --tarball /tmp/archimedlight-release-fixtures.tar.gz \
  --refresh-references
```

Optional: bind the artifact in `Artifacts.toml` by providing the download URL:

```bash
julia --project=. scripts/build_release_fixture_artifact.jl \
  --test-root /path/to/heavy-test-root \
  --tarball /tmp/archimedlight-release-fixtures.tar.gz \
  --refresh-references \
  --url https://github.com/ARCHIMED-platform/ArchimedLight.jl/releases/download/v0.0.1/archimedlight-release-fixtures-v0.0.1.tar.gz
```

Run release-only heavy regression:

```bash
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests("test"; filter=ti -> :release in ti.tags, verbose=true)'
```

You can also bypass artifacts and point directly to a local extracted release dataset:

```bash
ARCHIMEDLIGHT_RELEASE_FIXTURES_DIR=/path/to/release-fixtures \
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests("test"; filter=ti -> :release in ti.tags, verbose=true)'
```

Run one release fixture directly by tag. Fixture ids are exposed as tags with `-` replaced by `_`:

```bash
ARCHIMEDLIGHT_RELEASE_FIXTURES_DIR=/path/to/release-fixtures \
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests("test"; filter=ti -> :test_compare_simpleplant in ti.tags, verbose=true)'
```

Release test scripts are local in this repository (`test/release/`) and consume only data from the
artifact/dataset. During release runs, per-fixture progress is logged with start/end timestamps.

The regression matrix also has an optional release profile that reuses the same dataset root:

```bash
ARCHIMEDLIGHT_REGRESSION_PROFILE=release \
ARCHIMEDLIGHT_RELEASE_FIXTURES_DIR=/path/to/release-fixtures \
julia --project=test test/regression_matrix/runtests.jl
```

## Benchmarks

The repository includes a separate benchmark project for `AirspeedVelocity.jl` under
`benchmark/`.

## How to Cite

If you use ArchimedLight.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/VEZY/ArchimedLight.jl/blob/main/CITATION.cff).

---

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
