# Full Featured Example

This example is a self-contained light-only run using copied fixture data.

## Files
- `config.yml`: ARCHIMED-style configuration
- `meteo.csv`: weather input with PAR and custom band
- `scene/simple.ops` and `scene/opf/simple_OPF_shapes.opf`: geometry scene
- `model_simple.yml`, `model_soil.yml`: optical properties
- `full_featured_example.jl`: runnable end-to-end script

## Run
From the repository root:

```bash
julia --project=. example_1/full_featured_example.jl
```

The script demonstrates:
- stage-by-stage calls (`compute_sky` -> `integrate_light`)
- transfer-graph build and reuse (`build_scattering_transfer_graph`)
- backend selection in pipeline calls (`run_light_step`, `run_light_series`)
- custom waveband handling (`RI_custom_f`)
- 3D scene rendering with PlantGeom colored by intercepted PAR (`example_1/output/scene_3d_par_intercepted.png`)
- optional reference comparison using `component_values.csv` from `example_1/output/000001/`
- per-component comparison table: `example_1/output/component_values_reference_vs_julia.csv`
- 3D scene colored with reference PAR: `example_1/output/scene_3d_par_reference.png`

## Reference Comparison Input
If you have a reference `component_values.csv` for this scene, place it at:

- `example_1/output/000001/component_values.csv`

When this file exists, the script automatically computes reference-vs-Julia per-component differences.
