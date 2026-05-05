# Alpha Tester Migration

`ArchimedLight.jl` now uses `LightSimulation` and `run_light` as the main API.
The old staged API is still useful for debugging, but it is no longer the
recommended first path.

## File-Based Runs

Old:

```julia
options, scene, meteo, models = read_config("config.yml")
row = first(prepare_meteo(meteo, options).rows)
step = run_light_step(scene, models, row, options)
```

New:

```julia
sim, meteo = read_simulation("config.yml")
step = run_light(sim, first(meteo))
```

For a full series:

```julia
series = run_light(sim, meteo)
```

## Interactive Runs

Old:

```julia
scene = prepare_scene(mtg; scene_xy_bounds=bounds)
models = prepare_models(groups)
step = run_light_step(scene, models, row, options)
```

New:

```julia
scene = prepare_scene(mtg; scene_xy_bounds=bounds)
models = models_for(
    "coffee" => (
        "Leaf" => translucent(par=0.15, nir=0.90),
        "Stem" => translucent(par=0.20, nir=0.50),
    ),
)
sim = LightSimulation(scene, models; options=options)
step = run_light(sim, row)
```

You can also build the scene with `light_scene`:

```julia
scene = light_scene(domain=bounds) do s
    add_plant!(s, plant_mtg; group="coffee", id=1)
    add_ground!(s; group="soil", type="ground")
end
```

## Coupled Model Loops

Old:

```julia
cache = prepare_light_cache(scene, models, options)
for row in rows
    step = run_light_step(cache, row)
end
```

New:

```julia
sim = LightSimulation(scene, models; options=options)
for row in rows
    step = run_light(sim, row)
end
```

When the host model changes the geometry:

```julia
update_scene!(sim, new_scene)
```

This immediately releases the old scene-dependent prepared data and cache
entries. The next `run_light` call prepares the new scene lazily.
