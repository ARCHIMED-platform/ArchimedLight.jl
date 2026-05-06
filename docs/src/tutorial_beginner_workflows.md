# Beginner Workflows

This page shows the three common ways to use `ArchimedLight.jl`.

The central idea is always the same:

```text
Scene + Models + Options -> LightSimulation
LightSimulation + Meteo -> light results
```

## 1. Run From Files

Use this when you already have a `config.yml` that points to a scene, model
files, and a meteo file.

```julia
using ArchimedLight

sim, meteo = read_simulation("config.yml")

step = run_light(sim, first(meteo))
series = run_light(sim, meteo)
```

`step` is a `LightStepResult`. The most common values are:

```julia
step.budget.incident_flux.total.par
step.budget.incident_energy.total.par
step.budget.absorbed_energy.total.par
```

Use `incident_flux` for irradiance in `W m^-2`, and `incident_energy` or
`absorbed_energy` for energy per component and timestep in `J`.

## 2. Build Inputs Interactively

Use this when the plant is created in Julia or imported object by object.

```julia
using ArchimedLight
using FileIO, MeshIO

sensor_mesh = load("sensor.obj")

scene = light_scene(domain=(-1.0, -1.0, 1.0, 1.0)) do s
    add_plant!(s, "plant.opf"; group="coffee", id=1, at=(0.0, 0.0, 0.0), scale=0.9)
    add_object!(s, sensor_mesh; group="sensor", type="panel", id=10, at=(0.5, 0.0, 1.2), rotate=(z=90.0,), deg=true)
    add_ground!(s; group="soil", type="ground", nx=20, ny=20)
end
```

When the object starts as a mesh file, read it with MeshIO/FileIO before adding
it to the scene:

```julia
mesh = load("other-object.ply")
```

The `domain` is the XY plot footprint used by the raster projection. It is not
only a plotting box. With `toricity=true`, it is the repeated tile footprint.

Models match geometric scene nodes by functional group and component type:

```julia
models = models_for(
    "coffee" => (
        "Leaf" => translucent(par=0.15, nir=0.90),
        "Stem" => translucent(par=0.20, nir=0.50),
    ),
    "soil" => (
        "ground" => translucent(par=0.10, nir=0.40),
    ),
)
```

Only nodes with geometry need a model. Grouping nodes such as plant, axis, or
phytomer nodes do not.

Meteo can come from a file, `PlantMeteo`, or any Tables.jl-compatible table:

```julia
meteo = [
    (
        date="2020/06/21",
        hour_start="12:00:00",
        hour_end="13:00:00",
        latitude=15.0,
        RI_PAR_f=350.0,
        RI_NIR_f=250.0,
    ),
]

sim = LightSimulation(scene, models; options=LightOptions(pixel_size=0.003))
report = check_simulation(scene, meteo; models=models, options=sim.options)
isempty(report.errors) || error(report.errors)

step = run_light(sim, first(meteo))
```

## 3. Couple With Another Model

Use one `LightSimulation` while the geometry is unchanged. The simulation owns
the prepared geometry and the optional radiation cache.

```julia
sim = LightSimulation(scene, models; options=LightOptions(cache_radiation=true))

for row in meteo_rows
    light = run_light(sim, row)
    # pass `light` to the host model
end
```

When growth, pruning, or organ geometry changes the scene, update the
simulation. Old scene-dependent caches are released immediately; the next
`run_light` prepares the new scene lazily.

```julia
update_scene!(sim, new_scene)

for row in next_meteo_rows
    light = run_light(sim, row)
end
```

If optical models or projection options change, use:

```julia
update_models!(sim, new_models)
update_options!(sim, new_options)
```

## Before Running

For first-time inputs, run the checks:

```julia
summarize_scene(scene; models=models)
summarize_meteo(meteo; options=sim.options)
check_scene(scene)
check_models(scene, models)
check_meteo(meteo; options=sim.options)
check_simulation(scene, meteo; models=models, options=sim.options)
```

Fatal issues include missing scene domain, missing models for geometric nodes,
missing radiation inputs, or missing latitude when the sun position must be
reconstructed from date and time.

The summaries answer "what did ArchimedLight see?" The checks answer "is that
enough to run?"
