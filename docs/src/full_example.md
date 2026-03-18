# Full Featured Example

A runnable, self-contained example is available in:

- `example_1/full_featured_example.jl`
- `example_1/README.md`

It uses copied fixture input files under `example_1/`.

## Run
From repository root:

```bash
julia --project=. example_1/full_featured_example.jl
```

## What it demonstrates
- Input loading: `read_light_config`, `read_scene`, `read_meteo`
- Stage API: `compute_sky`, `build_turtle`, `compute_directional_fluxes`, `compute_first_order`, `compute_scattering`, `integrate_light`
- Transfer-graph stage: `build_scattering_transfer_graph`
- Pipeline API with backend kwargs: `run_light_step`, `run_light_series`
- Extra waveband handling (custom band from `RI_custom_f`)
- 3D visualization with PlantGeom + CairoMakie colored by intercepted PAR (saved to `example_1/output/scene_3d_par_intercepted.png`)
- optional reference-vs-Julia per-component comparison using `example_1/output/000001/component_values.csv` (saved to `example_1/output/component_values_reference_vs_julia.csv`)
- 3D visualization colored by reference PAR (saved to `example_1/output/scene_3d_par_reference.png`)
