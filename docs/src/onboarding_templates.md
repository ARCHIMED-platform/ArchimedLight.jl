# Copy-Paste Templates

These templates are intentionally small. Copy one into your project, then
change paths, group names, model coefficients, and options.

The same files are available under `templates/` in the repository.

## Recommended Project Layout

```text
my-light-run/
  config.yml
  run.jl
  scene/
    plant.opf
  meteo/
    meteo.csv
  models/
    coffee.yml
    soil.yml
  outputs/
```

The simplest file-based script is:

```julia
using ArchimedLight

sim, meteo = read_simulation("config.yml")

summarize_scene(sim.scene; models=sim.models)
summarize_meteo(meteo; options=sim.options)

series = run_light(sim, meteo)
```

## File-Based Simulation

```julia
using ArchimedLight

sim, meteo = read_simulation("config.yml")

report = check_simulation(sim.scene, meteo; models=sim.models, options=sim.options)
isempty(report.errors) || error(join(report.errors, "\n"))

series = run_light(sim, meteo)
step = first(series)
```

## Interactive Simulation

```julia
using ArchimedLight
using PlantGeom
using DataFrames

bounds = (-1.0, -1.0, 1.0, 1.0)

scene = make_scene(domain=bounds) do s
    add_plant!(s, "scene/plant.opf"; group="coffee", id=1, at=(0.0, 0.0, 0.0))
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

meteo = DataFrame(
    date=["2020/06/21"],
    hour_start=["12:00:00"],
    hour_end=["13:00:00"],
    latitude=[15.0],
    RI_PAR_f=[350.0],
    RI_NIR_f=[250.0],
)

options = LightOptions(pixel_size=0.01, scattering=true)
sim = LightSimulation(scene, models; options=options)

summarize_scene(scene; models=models)
summarize_meteo(meteo; options=options)

step = run_light(sim, first(eachrow(meteo)))
```

## Coupled Model Loop

```julia
using ArchimedLight

sim = LightSimulation(scene, models; options=LightOptions(cache_radiation=true))

for row in meteo_rows
    light = run_light(sim, row)
    # pass `light` to the host model
end

# When the host model changes geometry:
update_scene!(sim, new_scene)

for row in next_meteo_rows
    light = run_light(sim, row)
end
```

## Explicit SkyState

```julia
using ArchimedLight

sky = SkyState(
    135.0,  # sun azimuth, degrees
    60.0,   # sun elevation, degrees
    350.0,  # PAR irradiance, W m^-2
    250.0,  # NIR irradiance, W m^-2
    0.60,   # direct fraction
    0.40,   # diffuse fraction
)

sim = LightSimulation(scene, models; options=options)
step = run_light(sim, sky; step_duration_seconds=1800.0)
```
