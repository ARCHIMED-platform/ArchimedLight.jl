# File-Based Example

A runnable, self-contained example is available in:

- `example_1/full_featured_example.jl`
- `example_1/README.md`

It uses the compact input files under `example_1/` and only public
`ArchimedLight` APIs.

## Run

From repository root:

```bash
julia --project=. example_1/full_featured_example.jl
```

## What it demonstrates

- Load the complete file-based simulation with `read_simulation`
- Inspect compact scene and meteorology summaries
- Run one meteo row with `run_light`
- Inspect the concise `LightStepResult` summary
- Attach PAR results with `attach_light_step!`
- Export `output/component_values.csv` with `write_component_values`
- Export `output/scene_with_light.opf` with `write_scene`

The generated `output/` directory is intentionally untracked.
