# File-Based Example

This directory contains a small, self-contained light-interception simulation.
It uses the same public workflow as a user project: load a configuration with
`read_simulation`, then pass one meteo row to `run_light`.

## Files

- `config.yml`: ARCHIMED-style options file
- `meteo.csv`: weather input with PAR and custom band
- `scene/simple.ops` and `scene/opf/simple_OPF_shapes.opf`: geometry scene
- `model_simple.yml`, `model_soil.yml`: optical properties
- `full_featured_example.jl`: runnable one-step simulation

## Run

From the repository root:

```bash
julia --project=. example_1/full_featured_example.jl
```

The script:

- loads the scene, models, options, and meteorology with `read_simulation`
- prints compact scene and meteo summaries
- runs the first meteo row with `run_light`
- prints the resulting `LightStepResult` summary
- attaches PAR results to the scene with `attach_light_step!`
- writes `output/component_values.csv`
- writes `output/scene_with_light.opf`

The `output/` directory is generated when the example runs and is not tracked.
