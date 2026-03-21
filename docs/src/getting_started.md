# Getting Started

This page is the shortest path to your first light simulation.
It uses the bundled coffee example from `example_2/`, which is also the source of the image on the home page.

![Coffee example scene](assets/coffee_scene.png)

## What This Example Covers

- reading parameters from files, using a `config.yml` as the single source of truth for file paths and model options
- loading a scene from `.ops` and `.opf`
- loading functional-group models from YAML
- reading one meteo step
- running one light step
- attaching `Ri_PAR_f` to the scene for visualization or export

## Minimal Run

```@setup getting_started
using ArchimedLight, PlantGeom, CairoMakie
CairoMakie.activate!(type = "png")
```

```@example getting_started
using ArchimedLight

repo_root = normpath(joinpath(dirname(pathof(ArchimedLight)), ".."))
config = joinpath(repo_root, "example_2", "config.yml")
options, scene, meteo, models = read_config(config)
row = first(prepare_meteo(meteo, options).rows)
step = run_light_step(scene, models, row, options)

step
```

The result is a `LightStepResult`. The most useful field at first is `step.budget`, which groups values by:

- quantity: `incident` or `absorbed`
- form: `flux` (`W m^-2`) or `energy` (`J` per component per step)
- stage: `initial` for first-order only, `total` for first-order plus scattering
- waveband: `par`, `nir`, and any extra shortwave bands present in the meteo file

## Attach Results Back Onto The Scene

`ArchimedLight.jl` keeps the numeric results in Julia objects by default. If you want an inspectable scene, attach selected outputs back onto the MTG:

```@example getting_started
attach_light_step!(
    scene,
    step;
    fields=[:incident_par_flux, :incident_par_energy, :absorbed_par_energy],
)

out_path = joinpath(mktempdir(), "coffee_step.opf")
write_scene(out_path, scene)

out_path
```

The attached node attributes use the standard ARCHIMED names:

- `:incident_par_flux` -> `Ri_PAR_f`
- `:incident_par_energy` -> `Ri_PAR_q`
- `:absorbed_par_energy` -> `Ra_PAR_q`


# Plotting the results

You can plot the results using `PlantGeom.jl` + `Makie.jl`:

```@example getting_started
using PlantGeom, CairoMakie

fig, ax, p = plantviz(
    scene.mtg;
    color=:Ri_PAR_f,
    colormap=:thermal,
    figure=(size=(980, 700),),
)

PlantGeom.colorbar(fig[1, 2], p, label="Ri_PAR_f (W m^-2)")
fig
```

By default, the `toricity` parameter is activated, this is why we see the shade of the coffee plant coming from all corners, because light that goes out of the scene on one side comes back in on the other side. This is done for simulating an infinite canopy, but it can be turned off with `toricity=false` in the config.

## What To Read Next

- [File-Based Workflow](tutorial_files.md) explains how the example directory is organised and how each file contributes to the simulation.
- [Interactive Workflow](tutorial_interactive.md) shows how to build a scene directly in Julia and run light on it without `.ops` / `.opf` input files.
- [Configuration And Options Reference](reference_config.md) lists the configuration keys that the current Julia package actually uses.

## One Important Difference From The Historical Java App

The coffee example still contains several Java-era keys related to outputs, photosynthesis, and energy balance. `ArchimedLight.jl` keeps reading the light configuration and the input file paths, but the package documented here is currently the light-only core. For current behavior, prefer the reference pages on this site over older Java screenshots or YAML examples.
