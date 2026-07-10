# ArchimedLight.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://VEZY.github.io/ArchimedLight.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://VEZY.github.io/ArchimedLight.jl/dev)
[![Test workflow status](https://github.com/VEZY/ArchimedLight.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/VEZY/ArchimedLight.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/VEZY/ArchimedLight.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/VEZY/ArchimedLight.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![All Contributors](https://img.shields.io/github/all-contributors/VEZY/ArchimedLight.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

Julia reimplementation of the ARCHIMED light interception pipeline with a composable, function-first API.

| A wheat plant | A wheat AgriPV system |
|:---:|:---:|
| [![A wheat plant](docs/src/assets/archimedlight_day_cycle_1.png)](docs/src/assets/archimedlight_day_cycle_1.mp4) | [![A wheat AgriPV system](docs/src/assets/archimedlight_day_cycle.png)](docs/src/assets/archimedlight_day_cycle.mp4) |

<details>
<summary>Reproducing the figures</summary>

The figures above are generated from two example scenes available from an Artifact. You can generate them using the script in `docs/make_video.jl`.
</details>

## Installation

ArchimedLight.jl requires Julia 1.10 or newer. Install the registered package
from the Julia package manager:

```julia
using Pkg
Pkg.add("ArchimedLight")
```

Plotting is provided through the optional Makie extension. Install a Makie
backend such as GLMakie when you need figures:

```julia
Pkg.add("GLMakie")
```

## Scope

ARCHIMED computes:

- A pipeline for managing input scene, models and meteorology
- A discretization of the sky with a turtle approach (directional fluxes)
- The first-order light interception on the CPU using a rasterization approach
- Iterative scattering between scene components making links from the projection of the first-order interception
- A lot of helpers for reading/writing files, preparing scenes and models, exporting and visualizing the scene with Makie

## Core API

Here is an example workflow based on input files:

```julia
using ArchimedLight

repo_root = normpath(joinpath(dirname(pathof(ArchimedLight)), ".."))
config = joinpath(repo_root, "example_2", "config.yml")
sim, meteo = read_simulation(config)

summarize_scene(sim.scene; models=sim.models)
summarize_meteo(meteo; options=sim.options)

step = run_light(sim, first(meteo))
series = run_light(sim, meteo)
```

!!! note
    You can find these files in the `example_2` folder of the repository. The `config.yml` file is a convenient entrypoint for file-driven workflows. It points to the scene, models, and meteo files, and it can also contain optional light options.

`LightSimulation` owns the prepared scene and optional radiation cache. If another model changes
the scene, you must update it for ARCHIMED explicitly. For example if your other model provides `new_scene`, you would have to run:

```julia
update_scene!(sim, new_scene)
step = run_light(sim, next_meteo_row)
```

In Julia code, the result budget is grouped by quantity and waveband. These
schematic accesses assume the `step` from the checked workflow above:

```julia
step.budget.incident_flux.total.par
step.budget.incident_energy.total.par
step.budget.absorbed_flux.total.nir
step.budget.absorbed_energy.initial.par
```

File exports and attached MTG attributes use the java-version ARCHIMED names:

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

## Interactive inputs

You can also build everything in Julia. The following is a schematic shape for
host applications that already provide `sensor.obj`, `plant.opf`, and
`meteo_row`; the interactive workflow page contains the checked in-memory
example used by the documentation build.

```julia
using ArchimedLight
using PlantGeom
using FileIO, MeshIO

sensor_mesh = load("sensor.obj")

scene = make_scene(domain=(-1.0, -1.0, 1.0, 1.0)) do s
    add_plant!(s, "plant.opf"; group="coffee", id=1, at=(0.0, 0.0, 0.0), rotate=(z=25.0,), deg=true)
    add_object!(s, sensor_mesh; group="sensor", type="panel", id=10, at=(0.5, 0.0, 1.2), scale=0.1)
    add_ground!(s; group="soil", type="ground", nx=20, ny=20)
end

models = models_for(
    "coffee" => (
        "Leaf" => translucent(par=0.15, nir=0.90),
        "Stem" => translucent(par=0.20, nir=0.50),
    ),
    "soil" => (
        "ground" => translucent(par=0.10, nir=0.40),
    ),
)

sim = LightSimulation(scene, models; options=LightOptions())
step = run_light(sim, meteo_row)
```

## Advanced pipeline

The explicit stages remain available for debugging, parity work, and custom
research workflows. This snippet is schematic because it assumes stage inputs
such as `row`, `options`, `scene`, and `models` already exist.

```julia
sky = ArchimedLight.compute_sky(row, options)
turtle = ArchimedLight.build_turtle(options, sky)
fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
first_order = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
scat = ArchimedLight.compute_scattering(scene, models, turtle, first_order, options)
budget = ArchimedLight.integrate_light(scene, models, first_order, scat, options; meteo_row=row)
```

For ordinary simulations prefer `LightSimulation` and `run_light`.

## File-Based Example

- A self-contained one-step file workflow is under `example_1/`.
- It uses `read_simulation`, `run_light`, and the public output helpers.
- Run with:

```bash
julia --project=. example_1/full_featured_example.jl
```

## Stage flexibility

- `read_simulation("config.yml")` is the convenience entrypoint for file-driven workflows and returns `sim, meteo`.
- You can call each stage independently through the module (`ArchimedLight.compute_sky`, `ArchimedLight.build_turtle`, `ArchimedLight.compute_first_order`, `ArchimedLight.compute_scattering`, `ArchimedLight.integrate_light`, ...).
- File-based and in-memory workflows share the same runtime path: `read_scene(...)`, `PlantGeom.prepare_scene(...)`, and `PlantGeom.make_scene(...)` all produce `PlantGeom.SceneGeometry`, while `read_models(...)`, `prepare_models(...)`, and `models_for(...)` produce `LightModels`.
- `compute_sky` follows the ARCHIMED clearness/global conversion and DeJong hourly direct/diffuse partitioning.
- `compute_sky` uses substep-weighted sun position (`radiation_timestep`) when sun angles are not provided.
- Meteo `#' use: ...` consistency checks for `clearness`/`RI_SW_f`/`RI_PAR_f`/`RI_NIR_f` are enforced like Java.
- `compute_first_order(...; backend=:raster_cpu)` is the current reference backend (`RasterCPUBackend()` also available).
- `compute_scattering(...; mode=:raycast)` is the current scattering mode; `RaycastScatteringBackend()` can also be passed explicitly.
- `pixel_size` is validated with ARCHIMED-compatible bounds (`0 < pixel_size <= 0.5` meters).
- `cache_pixel_table=true` enables an on-disk direction projection cache under the interception cache directory.
- `build_turtle` follows the canonical ARCHIMED sector counts `1, 6, 16, 46, 136, 406`.
- `PlantGeom.add_ground!(scene; ...)` is an explicit scene-editing step, so inspectable ground is ordinary geometry in the exported MTG.
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

## Releasing

Maintainer release steps, including docs, release fixtures, regression matrix,
benchmarks, Julia General registration, TagBot, and artifact upload, are in
[`RELEASE.md`](RELEASE.md).

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
