# Full Featured Example

A runnable, self-contained example is available in:

- `example/full_featured_example.jl`
- `example/README.md`

It uses copied fixture input files under `example/` so it does not depend on paths inside `java_implementation/`.

## Run
From repository root:

```bash
julia --project=. example/full_featured_example.jl
```

## What it demonstrates
- Input loading: `read_light_config`, `read_scene`, `read_meteo`
- Stage API: `compute_sky`, `build_turtle`, `compute_directional_fluxes`, `compute_first_order`, `compute_scattering`, `integrate_light`
- Transfer-graph stage: `build_scattering_transfer_graph`
- Pipeline API with backend kwargs: `run_light_step`, `run_light_series`
- Extra waveband handling (custom band from `RI_custom_f`)
- 3D visualization with PlantGeom + CairoMakie colored by intercepted PAR (saved to `example/output/scene_3d_par_intercepted.png`)
- Java vs Julia per-component comparison using `example/output/000001/component_values.csv` (saved to `example/output/component_values_java_vs_julia.csv`)
- 3D visualization colored by Java PAR (saved to `example/output/scene_3d_par_java.png`)
