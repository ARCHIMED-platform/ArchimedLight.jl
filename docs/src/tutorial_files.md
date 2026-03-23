# Tutorial: File-Based Workflow

This tutorial is the standard ARCHIMED workflow: a scene on disk, one or more model files, a meteo CSV, and a `config.yml` that ties them together.

It uses the deliberately small fixture under `example_1/`, because every file is compact enough to read directly.

![Simple plant scene](assets/example_simpleplant_scene.png)

## Directory Layout

The example folder contains:

```text
example_1/
├── config.yml
├── meteo.csv
├── model_simple.yml
├── model_soil.yml
├── full_featured_example.jl
└── scene/
    ├── simple.ops
    └── opf/simple_OPF_shapes.opf
```

The roles of those files are:

- `config.yml`: top-level wiring, global light options, and file paths
- `simple.ops`: scene placement, plot bounds, and functional-group assignment
- `simple_OPF_shapes.opf`: plant topology and meshes
- `model_simple.yml`, `model_soil.yml`: optical properties by functional group and type
- `meteo.csv`: one or more simulation time steps

## Step 1: Load Everything From The Config

```julia
using ArchimedLight

config = joinpath("example_1", "config.yml")
options, scene, meteo, models = read_config(config)
rows = prepare_meteo(meteo, options).rows
```

At this point:

- `options` is a `LightOptions`
- `scene` is a `SceneGeometry`
- `meteo` is a `MeteoTable`
- `models` is a `LightModels`

## Step 2: Run One Step Or A Whole Series

For a single time step:

```julia
row = first(rows)
step = run_light_step(scene, models, row, options)
```

For every meteo row:

```julia
series = run_light_series(scene, models, meteo, options)
```

When `LightOptions(cache_radiation=true)` is enabled, `run_light_series` can reuse directional responses between meteo rows on the same scene.

## Step 3: Inspect The Budget

```julia
budget = step.budget

budget.incident_flux.initial.par
budget.incident_flux.total.par
budget.absorbed_energy.total.par
```

Common reading patterns are:

- first-order only: `budget.incident_flux.initial.par`
- with scattering: `budget.incident_flux.total.par`
- per-step energy: `budget.incident_energy.total.par`
- absorbed quantities: `budget.absorbed_flux.total.par`

## Step 4: Write Results Back Onto The Scene

Attach the values you want to inspect:

```julia
attach_light_step!(
    scene,
    step;
    fields=[:incident_par_flux, :incident_par_energy, :absorbed_par_energy],
)

write_scene(joinpath("example_1", "output", "scene_with_light.opf"), scene)
```

This keeps the historical ARCHIMED attribute names on the nodes:

- `Ri_PAR_f`
- `Ri_PAR_q`
- `Ra_PAR_q`

![Simple plant coloured by intercepted PAR](assets/example_simpleplant_scene_ri_par_q.png)

## Step 5: Use The Example Script When You Want More

`example_1/full_featured_example.jl` goes further than the minimal workflow:

- it runs the stage API explicitly
- it renders the scene with `PlantGeom`
- it compares Julia outputs against a Java `component_values.csv` when one is present

That script is a good template when you want a reproducible end-to-end example rather than an interactive session.

## Reading The Input Files Together

The file workflow is easiest to understand when you see the four inputs as one linked dataset:

1. `config.yml` selects the scene, meteo file, and model files, and defines global runtime options such as `sky_sectors`, `pixel_size`, and `toricity`.
2. `simple.ops` places the plant object into the plot and assigns it to the functional group `simple_plant`.
3. `model_simple.yml` says how `simple_plant` components intercept and scatter PAR, NIR, and the custom waveband.
4. `meteo.csv` provides the forcing for the simulation step.

The reference pages expand each one:

- [Configuration And Options Reference](reference_config.md)
- [Scene Files And Semantics](reference_scene.md)
- [Model Files Reference](reference_models.md)
- [Meteo Inputs Reference](reference_meteo.md)

## Current Julia Behavior Vs Historical ARCHIMED YAML

The examples intentionally preserve the old ARCHIMED config style. That is useful for compatibility, but not every historical key is active in the current package.

In particular:

- scene, model, and meteo paths are active
- the light options described on the configuration reference page are active
- many Java-era output toggles remain as documentation or archival context
- photosynthesis and energy-balance settings are not part of the current `ArchimedLight.jl` runtime
