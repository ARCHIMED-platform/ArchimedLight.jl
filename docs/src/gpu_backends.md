```@meta
CurrentModule = ArchimedLight
```

# GPU Backends And Benchmarks

ArchimedLight can run the Raycore backend through
`KernelAbstractions.jl`. This is the path used for portable GPU experiments:
the same ArchimedLight code can target Metal, CUDA, oneAPI, AMDGPU, or the
KernelAbstractions CPU backend when the matching array package is loaded by the
user environment.

The GPU path is intended for reusable simulations. Build a `LightSimulation`
once, then call `run_light` for one row or a complete meteo table. One-shot
calls such as `run_light_step` include backend setup costs and are not a good
steady-state GPU benchmark.

## Minimal Setup

The CPU raster backend is still the reference path:

```julia
using ArchimedLight

sim, meteo = read_simulation("config.yml")
step = run_light(sim, first(meteo))
```

To use Raycore on the KernelAbstractions CPU backend:

```julia
using ArchimedLight
using KernelAbstractions

interception = RaycoreInterceptionBackend(backend=KernelAbstractions.CPU())
sim, meteo = read_simulation(
    "config.yml";
    interception_backend=interception,
    scattering_backend=RaycoreScatteringBackend(interception),
)

step = run_light(sim, first(meteo))
```

For Metal, load `Metal.jl` in the environment that runs the simulation and pass
the backend returned by `KernelAbstractions.get_backend`:

```julia
using ArchimedLight
using KernelAbstractions
using Metal

metal_backend = KernelAbstractions.get_backend(MtlArray(zeros(Float32, 1)))
interception = RaycoreInterceptionBackend(
    backend=metal_backend,
    max_hits_per_pixel=32,
)

sim, meteo = read_simulation(
    "config.yml";
    interception_backend=interception,
    scattering_backend=RaycoreScatteringBackend(interception),
)

step = run_light(sim, first(meteo))
```

!!! note
    Metal, CUDA, oneAPI, and AMDGPU support are optional environment features.
    ArchimedLight core does not load those packages for you. Add the GPU package
    to the Julia environment that runs the simulation or benchmark.

When `workgroupsize` is not provided, ArchimedLight currently uses `256`.
Metal uses this default based on the bundled coffee cached benchmark. Pass
`workgroupsize=...` explicitly when benchmarking a different device or scene.

## Cache Behavior

`LightSimulation` keeps the Raycore scene data alive across calls. That avoids
rebuilding Raycore acceleration structures and also lets ArchimedLight reuse
direction-specific CPU-side quantities such as projection extents and projected
mesh areas.

Use this pattern for repeated runs:

```julia
sim, meteo = read_simulation(
    "config.yml";
    interception_backend=interception,
    scattering_backend=RaycoreScatteringBackend(interception),
)

first_step = run_light(sim, first(meteo))
series = run_light(sim, meteo)
```

Avoid this pattern when comparing steady-state GPU performance:

```julia
# This rebuilds backend state for each call and includes setup cost.
step = ArchimedLight.run_light_step(
    sim.scene,
    sim.models,
    first(meteo),
    sim.options;
    interception_backend=interception,
    scattering_backend=RaycoreScatteringBackend(interception),
)
```

## Benchmark Script

Use `benchmark/backend_comparison.jl` for reproducible local comparisons. The
script reports time and allocation volume for each row.

CPU reference:

```sh
ARCHIMEDLIGHT_BENCH_WORKLOAD=coffee_home \
ARCHIMEDLIGHT_BENCH_EXECUTION=cached \
ARCHIMEDLIGHT_BENCH_BACKENDS=normal_cpu \
ARCHIMEDLIGHT_BENCH_COMPONENTS=first_order,scattering \
julia --project=benchmark benchmark/backend_comparison.jl
```

Metal comparison:

```sh
ARCHIMEDLIGHT_BENCH_WORKLOAD=coffee_home \
ARCHIMEDLIGHT_BENCH_EXECUTION=cached \
ARCHIMEDLIGHT_BENCH_BACKENDS=normal_cpu,raycore_metal_gpu \
ARCHIMEDLIGHT_BENCH_COMPONENTS=first_order,scattering,topology,propagation,integration,stack_trace \
ARCHIMEDLIGHT_BENCH_METAL=1 \
ARCHIMEDLIGHT_BENCH_WARMUPS=1 \
ARCHIMEDLIGHT_BENCH_SAMPLES=3 \
julia --project=test/gpu benchmark/backend_comparison.jl
```

Other GPU packages are discovered only when explicitly requested. Add the
matching package to the active Julia environment, enable its flag, and select
the backend name:

```sh
ARCHIMEDLIGHT_BENCH_WORKLOAD=coffee_home \
ARCHIMEDLIGHT_BENCH_EXECUTION=cached \
ARCHIMEDLIGHT_BENCH_BACKENDS=raycore_cuda_gpu \
ARCHIMEDLIGHT_BENCH_COMPONENTS=first_order,scattering,topology \
ARCHIMEDLIGHT_BENCH_CUDA=required \
julia --project=<env-with-cuda> benchmark/backend_comparison.jl
```

The optional flags and backend labels are:

| package | flag | backend label |
|---|---|---|
| Metal.jl | `ARCHIMEDLIGHT_BENCH_METAL` | `raycore_metal_gpu` |
| CUDA.jl | `ARCHIMEDLIGHT_BENCH_CUDA` | `raycore_cuda_gpu` |
| oneAPI.jl | `ARCHIMEDLIGHT_BENCH_ONEAPI` | `raycore_oneapi_gpu` |
| AMDGPU.jl | `ARCHIMEDLIGHT_BENCH_AMDGPU` | `raycore_amdgpu_gpu` |

Use `1`, `true`, `yes`, or `on` to try a backend and warn if it cannot be
initialized. Use `required` or `force` to fail the benchmark when that backend
is unavailable.

One-shot comparison, including backend setup cost:

```sh
ARCHIMEDLIGHT_BENCH_WORKLOAD=coffee_home \
ARCHIMEDLIGHT_BENCH_EXECUTION=oneshot \
ARCHIMEDLIGHT_BENCH_BACKENDS=normal_cpu,raycore_metal_gpu \
ARCHIMEDLIGHT_BENCH_COMPONENTS=first_order,scattering \
ARCHIMEDLIGHT_BENCH_METAL=1 \
ARCHIMEDLIGHT_BENCH_WARMUPS=1 \
ARCHIMEDLIGHT_BENCH_SAMPLES=3 \
julia --project=test/gpu benchmark/backend_comparison.jl
```

The component rows mean:

- `first_order`: first-order interception only.
- `scattering`: first-order interception plus scattering, as a normal light
  step.
- `topology`: first-order interception plus construction of the scattering
  transfer topology.
- `propagation`: scattering propagation from an already-built topology.
- `integration`: light-budget integration from already-built first-order and
  scattering results, using the same cached geometry maps as the production
  cached pipeline.
- `stack_trace`: Raycore-only device stack tracing for positive directions,
  with label metrics for fixed hit-stack utilization, maximum observed stack
  depth, occupied pixels, and overflow.

On the bundled coffee scene used for the documentation landing-page image,
focused cached component runs give the following representative rows on Apple
Silicon:

| component | backend | median time | median allocation |
|---|---|---:|---:|
| first-order | normal CPU | 792.355 ms | 200.456 MiB |
| first-order | Raycore Metal | 134.202 ms | 54.601 MiB |
| scattering | normal CPU | 1367.455 ms | 1351.062 MiB |
| scattering | Raycore Metal | 599.936 ms | 128.378 MiB |
| topology | normal CPU | 1327.113 ms | 1332.727 MiB |
| topology | Raycore Metal | 508.786 ms | 70.754 MiB |
| propagation | normal CPU | 39.225 ms | 3.928 MiB |
| propagation | Raycore Metal | 47.906 ms | 3.796 MiB |
| stack trace | Raycore Metal | 90.545 ms | 0.181 MiB |

These numbers are local benchmark evidence, not a performance guarantee. They
are useful mainly for checking the expected trend: cached GPU runs should win on
scattering-heavy, repeated workloads, while dense standalone propagation and
tiny or one-shot scenes may still be better on the normal CPU backend.

The cached coffee stack-trace diagnostic used `max_hits=32` and observed
`hit_util=0.93%`, `max_seen=15`, `occupied=22.04%`, and `overflow=false`.
That means fixed hit-stack storage is sparse on this scene, but compact storage
should still be benchmarked carefully because it would require additional
passes or prefix-sum work.

## Parameter Sweeps

The benchmark script can sweep the parameters that usually determine whether a
GPU run wins:

```sh
ARCHIMEDLIGHT_BENCH_WORKLOAD=coffee_home \
ARCHIMEDLIGHT_BENCH_EXECUTION=cached \
ARCHIMEDLIGHT_BENCH_BACKENDS=raycore_metal_gpu \
ARCHIMEDLIGHT_BENCH_COMPONENTS=first_order \
ARCHIMEDLIGHT_BENCH_METAL=1 \
ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE_SWEEP=64,128,256 \
ARCHIMEDLIGHT_BENCH_WARMUPS=1 \
ARCHIMEDLIGHT_BENCH_SAMPLES=3 \
julia --project=test/gpu benchmark/backend_comparison.jl
```

Available sweep variables:

- `ARCHIMEDLIGHT_BENCH_PIXEL_SIZE_SWEEP`
- `ARCHIMEDLIGHT_BENCH_SECTOR_SWEEP`
- `ARCHIMEDLIGHT_BENCH_MAX_HITS_SWEEP`
- `ARCHIMEDLIGHT_BENCH_PLOT_PAVING_SWEEP`
- `ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE_SWEEP`
- `ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION_SWEEP`

Each variable is a comma-separated list. For the bundled coffee configuration,
sector count and pixel size come from `example_2/config.yml` unless you set the
matching override or sweep variable explicitly.

`ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION` and
`ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION_SWEEP` accept `auto`,
`sparse_host_reduce`, and `dense_atomic`. The default `auto` path uses the dense
atomic reducer only when the backend reports atomic support and the dense edge
matrix fits the configured memory limit; otherwise it falls back to the counted
sparse host reducer. Use the sweep when testing a new GPU or scene, because the
best reducer can depend on node count, pixel count, and backend atomics.

The dense atomic reducer keeps an `Int32` device scratch matrix in the cached
Raycore scene data. Dense host readback is still materialized per topology call,
because this is currently faster on Metal than copying into a retained host
buffer. A retained host aggregation vector was also tested and rejected because
it was slower on the cached coffee topology benchmark; a touched-index variant
was also slower and is not retained.

## When GPU Should Win

The Raycore GPU backend is most useful when:

- the same `LightSimulation` is reused for many meteo rows;
- the scene has enough projected pixels to amortize launch and transfer costs;
- scattering is enabled;
- the benchmark uses cached execution rather than one-shot setup.

CPU can still be faster for tiny scenes, very coarse grids, one-shot calls, and
cases where most time is spent converting dense internal arrays back to public
`Dict` results.

| Situation | Recommended path | Why |
|---|---|---|
| Tiny scene or one-off diagnostic run | normal CPU | Avoids Raycore setup, kernel compilation, and device transfer overhead. |
| Repeated meteo rows on the same scene | Raycore GPU with `LightSimulation` | Reuses acceleration structures, directional caches, and backend buffers. |
| Scattering-heavy coffee-like canopy | Raycore GPU | Traversal and propagation work amortize GPU launch and transfer costs. |
| Backend debugging or CI without GPU hardware | Raycore with `KernelAbstractions.CPU()` | Exercises the same Raycore/KA code path without Metal/CUDA/oneAPI/AMDGPU. |
| Historical parity investigation | normal CPU | The raster backend remains the reference behavior. |

## Validation Commands

Run the focused Raycore tests after backend changes:

```sh
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests(pwd(); filter=ti -> :raycore_backend in ti.tags, verbose=true)'
```

Run Metal validation on Apple Silicon:

```sh
ARCHIMEDLIGHT_TEST_METAL=required \
julia --check-bounds=auto --project=test/gpu test/gpu/runtests.jl
```
