# Alpha Tester Migration

## PlantSimEngine radiation-area boundary

The PlantSimEngine extension now distinguishes ArchimedLight flux densities
per radiative mesh area from physiology flux densities per botanical leaf
area. Raw fields keep their names, but `aPPFD`, `Ra_SW_f`, and
`radiative_mesh_area` now carry explicit scientific contracts.

Do not connect raw `aPPFD` directly to FvCB or raw `Ra_SW_f` directly to
Monteith. Give every destination leaf a finite positive
`botanical_leaf_area`, then use PlantBiophysics' `RadiativeMeshToLeafPPFD` and
`RadiativeMeshToLeafShortwave` adapters. They preserve
`raw_flux * radiative_mesh_area == leaf_flux * botanical_leaf_area`.

The `Ri_*` and component PAR/NIR diagnostic fields remain uncontracted, and
both output schemas retain their existing raw field names.

`ArchimedLight.jl` now uses `LightSimulation` and `run_light` as the main API.
The old staged API is still useful for debugging, but it is no longer the
recommended first path.

## File-Based Runs

Old:

```julia
options, scene, meteo, models = read_config("config.yml")
row = first(prepare_meteo(meteo, options))
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
scene = PlantGeom.prepare_scene(mtg; scene_xy_bounds=bounds)
models = prepare_models(groups)
step = run_light_step(scene, models, row, options)
```

New:

```julia
scene = PlantGeom.prepare_scene(mtg; scene_xy_bounds=bounds)
models = models_for(
    "coffee" => (
        "Leaf" => translucent(par=0.15, nir=0.90),
        "Stem" => translucent(par=0.20, nir=0.50),
    ),
)
sim = LightSimulation(scene, models; options=options)
step = run_light(sim, row)
```

You can also build the scene with `PlantGeom.make_scene`:

```julia
using PlantGeom

scene = make_scene(domain=bounds) do s
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

## Common Alpha Script Updates

Daily or seasonal scripts that used `run_light_series` now create a simulation
once and run the whole meteo table through it:

```julia
# Old
series = run_light_series(scene, models, meteo, options)

# New
sim = LightSimulation(scene, models; options=options)
series = run_light(sim, meteo)
```

Interactive scripts that already computed a [`SkyState`](@ref) no longer need
to assemble the internal solver stages:

```julia
# Old
turtle = ArchimedLight.build_turtle(options, sky)
fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
scat = ArchimedLight.compute_scattering(scene, models, turtle, first, options)
budget = ArchimedLight.integrate_light(scene, models, first, scat, options; step_duration_seconds=1800.0)

# New
sim = LightSimulation(scene, models; options=options)
step = run_light(sim, sky; step_duration_seconds=1800.0)
```

Scene placement now uses `at=` consistently:

```julia
# Old
place_in_scene!(plant; scene=scene_mtg, plant_id=1, functional_group="coffee", pos=(0, 0, 0))

# New
place_in_scene!(plant; scene=scene_mtg, plant_id=1, functional_group="coffee", at=(0, 0, 0))
```
