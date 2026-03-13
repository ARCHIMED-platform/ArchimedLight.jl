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
julia --project=. example/full_featured_example.jl
```

The script demonstrates:
- stage-by-stage calls (`compute_sky` -> `integrate_light`)
- transfer-graph build and reuse (`build_scattering_transfer_graph`)
- backend selection in pipeline calls (`run_light_step`, `run_light_series`)
- custom waveband handling (`RI_custom_f`)
- 3D scene rendering with PlantGeom colored by intercepted PAR (`example/output/scene_3d_par_intercepted.png`)
- Java vs Julia comparison using Java `component_values.csv` from `example/output/000001/`
- per-component comparison table: `example/output/component_values_java_vs_julia.csv`
- 3D scene colored with Java PAR: `example/output/scene_3d_par_java.png`

## Java Comparison Input
If you run the Java JAR in this folder, make sure its outputs are available at:

- `example/output/000001/component_values.csv`

When this file exists, the script automatically computes Java-vs-Julia per-component differences.
