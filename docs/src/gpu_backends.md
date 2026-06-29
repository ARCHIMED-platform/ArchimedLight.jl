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

!!! compat
    Raycore stack tracing targets the current `Raycore.all_hits!` contract:
    `all_hits!(metadata_out, distances_out, instance_indices_out, tlas, ray,
    out_base, max_hits, duplicate_epsilon)`. ArchimedLight does not keep a
    compatibility path for the earlier unreleased draft signature.

When `workgroupsize` is not provided, ArchimedLight currently uses `256`.
Treat this as a conservative default, not a universal optimum. The best value
can differ between full scattering and standalone topology rows, so pass
`workgroupsize=...` explicitly when benchmarking a different device, scene, or
component.

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

## Validation And Fallback

Before using a non-CPU Raycore backend in a cached `LightSimulation`,
ArchimedLight validates the backend with a vertical top-hit trace against a
normal CPU reference. If the backend misses too much CPU-visible top area, the
simulation falls back to `RasterCPUBackend` instead of returning silently
invalid light budgets.

This guard is especially important for large realistic scenes. On Apple
Silicon, the bundled coffee scene can fail an unchunked Metal TLAS trace even
though smaller scenes remain valid. ArchimedLight therefore prechunks large
non-CPU Raycore scenes before validation when their effective face count exceeds
`ARCHIMEDLIGHT_RAYCORE_PRECHUNK_FACE_THRESHOLD` (`500000` by default). The
chunk size is controlled by `ARCHIMEDLIGHT_RAYCORE_FACE_CHUNK_LIMIT` (`1536` by
default). Set either variable to `0` or a negative value to disable that
behavior.

Very large chunked scenes can be correct but impractical to build on a GPU. By
default, ArchimedLight falls back before Raycore construction when prechunking
would create more than `RaycoreBackendConfig(max_prechunk_instances=...)`
allows. The default comes from `ARCHIMEDLIGHT_RAYCORE_MAX_PRECHUNK_INSTANCES`
(`4096` when unset). Pass `0` or a negative value, either in the config or in
the environment variable, to disable the cap and force the guarded Raycore
validation path.

If a prechunked non-CPU Raycore scene fails validation, ArchimedLight retries
with smaller chunk limits before falling back, while still respecting
`max_prechunk_instances`. This handles occasional large-TLAS validation
instability without disabling the safety guard.

Benchmark labels report the resolved backend. A row labelled
`raycore_metal_gpu->RasterCPUBackend` means the requested GPU backend failed
validation and the measured row is the CPU fallback. A row labelled
`raycore_metal_gpu` with resolved backend `RaycoreInterceptionBackend` stayed on
Raycore. `cache_summary(sim)` also reports `fallback_reason`, currently
`:none`, `:raycore_prechunk_instance_cap`, or `:raycore_trace_validation`, so
large-scene runs can distinguish an intentional instance-cap fallback from a
trace-validation failure.

For toric scenes, ArchimedLight may also use raster-compatible projections for
individual low-elevation directions inside an otherwise Raycore-backed
simulation. This preserves the periodic full-stack semantics of the CPU
rasterizer when a finite Raycore toric copy would not cover the ray path. This
is not reported as a whole-backend fallback; the resolved backend remains
`RaycoreInterceptionBackend`.

When `RaycoreBackendConfig(propagation_backend=:auto)` is used,
KernelAbstractions CPU propagation and large GPU propagation graphs are routed
through the shared CPU dense solver. This keeps Raycore traversal/topology on
the requested backend while avoiding slow KA CPU or GPU atomic propagation
paths for large mixed toric topologies. Use `propagation_backend=:device` to
force device propagation, or `propagation_backend=:cpu` to force CPU
propagation. The GPU auto threshold is controlled by
`ARCHIMEDLIGHT_RAYCORE_CPU_PROPAGATION_EDGE_THRESHOLD` (`200000` edges by
default; set to `0` to disable).

For full-response caches, ArchimedLight batches CPU-dense PAR/NIR scattering
propagation across missing turtle sectors before assembling the requested sky
row. This keeps scalar device propagation unchanged, but avoids repeated
CPU-dense graph traversals when the cache is being populated for normal CPU,
KernelAbstractions CPU, or large Raycore GPU graphs that auto-route propagation
to CPU.

Internally, the full-response cache now keeps a backend-neutral dense response
storage instead of materializing per-sector `Dict`-style objects early. The
cache stores a node-by-sector projected-area matrix, active row indices per
sector, aggregate hit counts by node, dense emitter incident-power vectors, and
a compact scattering topology. Public `Dict` outputs are materialized only at
API boundaries such as `FirstOrderResult`, `ScatteringResult`, and light-budget
integration.

The dense storage also carries a device-cache slot for future backend-resident
copies. Today, Raycore/KernelAbstractions still performs several reductions and
public-output conversions on the host: sector active-index extraction, emitter
incident-power accumulation, compact edge-count materialization from traced
stacks, scattering transfer-graph construction, and final public `Dict`
construction. The parts already structured for future device-side reductions
are the projected-area matrix, aggregate hit vector, decoded stack-node buffers,
edge-key or dense-edge scratch buffers, and compact `edge_to`/`edge_from`/
`edge_count` topology arrays.

For Raycore full-response caches with scattering enabled and no emitters,
ArchimedLight also uses the flat stack-node trace buffers directly instead of
materializing a dense pixel-hit dictionary for each sector. Direction-specific
raster-compatible fallback sectors still use the CPU projection cache, and
emitter scenes keep the more general path.

Raster-compatible fallback projections are cached in memory on the
`RaycoreSceneData` object. Repeated cached topology/scattering component calls
therefore reuse the same low-elevation toric fallback projections instead of
rebuilding the CPU raster projection every time.

Toric Raycore scenes use TLAS instancing for the nine periodic plot copies.
ArchimedLight builds one base BLAS geometry for an unchunked toric scene, or
one base BLAS per face chunk for a chunked scene, then supplies the periodic
translation transforms to Raycore. Raycore hit metadata is not treated as a
final Archimed node id. ArchimedLight owns a hit-decoding table that maps
`(raycore triangle metadata, raycore instance index)` to an Archimed node index;
the stack kernels read Raycore instance indices, decode hits through this table,
and store decoded node indices before the dense response/cache code consumes
the stack.

When the prepared `PlantGeom.SceneGeometry` still exposes repeated classic
`PlantGeom.Geometry(ref_mesh, transformation)` nodes, Raycore can now build a
compact reference-mesh mode internally. ArchimedLight builds one BLAS prototype
per reused untapered reference mesh, pushes one TLAS instance per transformed
Archimed node, and adds the toric periodic copies as extra instances. Prototype
triangle metadata is local to the reference mesh; the decoder table maps the
same local metadata to different dense Archimed node indices depending on the
Raycore instance index. Non-reused nodes, generated paving, tapered reference
meshes, point-mapped geometry, missing transforms, and validation failures
automatically keep or return to the older merged-mesh/chunked path. Set
`ARCHIMEDLIGHT_RAYCORE_REFERENCE_INSTANCING=0` only for internal comparisons
against the merged-mesh Raycore path. Prototype extraction is lazy and only runs
for Raycore preparation, so normal CPU rendering and cache preparation do not
pay this internal Raycore-only cost.

The current reference-mesh mode is a structural step, not yet a guaranteed
speed win. Scene analysis, prototype grouping, transform conversion, decoder
table construction, and TLAS instance creation are still host-side. On scenes
with many small repeated nodes, the compact BLAS set can be outweighed by a very
large TLAS instance count; benchmark both
`ARCHIMEDLIGHT_RAYCORE_REFERENCE_INSTANCING=1` and `0` before using it as a
performance assumption. ArchimedLight therefore applies two internal guardrails
before using the candidate: `ARCHIMEDLIGHT_RAYCORE_REFERENCE_INSTANCE_LIMIT`
caps the total TLAS instance count after toric copies (`4096` by default, `0`
disables the cap), and
`ARCHIMEDLIGHT_RAYCORE_REFERENCE_MIN_FACE_SAVINGS_RATIO` requires a minimum
compact-face saving (`0.20` by default). Candidates that fail either guardrail
fall back automatically to the merged-mesh or chunked path.

## Benchmark Script

Use `benchmark/backend_comparison.jl` for reproducible local comparisons. The
script reports time and allocation volume for each row.
In cached mode, the table also reports the resolved fallback reason so a row
that requested GPU but measured CPU fallback is visible in the same output.

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

For realistic local build-vs-cache comparisons, use
`benchmark/local_realistic_backends.jl`. It compares the normal CPU raster
backend, Raycore on `KernelAbstractions.CPU()`, and any explicitly requested
GPU backend on the same workloads:

```sh
ARCHIMEDLIGHT_BENCH_METAL=required \
ARCHIMEDLIGHT_LOCAL_BENCH_WORKLOADS=coffee,agripv,xpalm \
ARCHIMEDLIGHT_LOCAL_BENCH_BUILD_SAMPLES=1 \
ARCHIMEDLIGHT_LOCAL_BENCH_SAMPLES=1 \
julia --check-bounds=auto --project=test/gpu benchmark/local_realistic_backends.jl
```

This script currently includes the bundled coffee scene used for the
documentation landing-page image, a local agrivoltaics wheat-panel scene, and an
opt-in XPalm/VPalm two-oil-palm scene. It prints fresh build+run timing,
cached-reuse timing, light totals, and the resolved backend. Use
`ARCHIMEDLIGHT_BENCH_BACKENDS` and the same optional GPU flags as
`benchmark/backend_comparison.jl` to select Metal, CUDA, oneAPI, or AMDGPU. The
printed rows include both the resolved backend and the fallback reason. Set
`ARCHIMEDLIGHT_BENCH_MAX_PRECHUNK_INSTANCES` to override the Raycore
`max_prechunk_instances` cap for these benchmarked backends. Set
`ARCHIMEDLIGHT_LOCAL_BENCH_WORKLOADS=coffee` when the local agrivoltaics
checkout is unavailable. Set `ARCHIMEDLIGHT_LOCAL_BENCH_WORKLOADS=xpalm` to run
the XPalm stress case; it requires the local XPalm checkout and can be much
slower because it may build a heavily chunked TLAS.

Set `ARCHIMEDLIGHT_LOCAL_BENCH_BREAKDOWN=1` to split the local realistic
benchmark into four phases:

- `construct`: `LightSimulation` construction.
- `prepare`: validation and `prepare_light_cache` through
  `ArchimedLight._ensure_light_cache!`.
- `populate`: the first run after preparation, which builds the full turtle
  response cache.
- `reuse`: the second run on the same prepared simulation.

Breakdown phase cells report median seconds and median allocated MiB.

Set `ARCHIMEDLIGHT_LOCAL_BENCH_PREPARE_BREAKDOWN=1` when the prepare phase is
the bottleneck. This keeps the same workload and backend selectors, but splits
preparation into geometry preparation, Raycore scene/TLAS construction,
validation, and retry. On the bundled coffee scene with Metal, this mode showed
that the Raycore scene/TLAS phase dominates prepare time. It also showed that
`ARCHIMEDLIGHT_RAYCORE_FACE_CHUNK_LIMIT=2048` fails validation
(`ratio=0.759`), while the current default `1536` passes (`ratio=0.992`) with
721 chunked BLAS instances.

For example:

```sh
ARCHIMEDLIGHT_LOCAL_BENCH_BREAKDOWN=1 \
ARCHIMEDLIGHT_LOCAL_BENCH_WORKLOADS=coffee \
ARCHIMEDLIGHT_BENCH_BACKENDS=normal_cpu,raycore_ka_cpu,raycore_metal_gpu \
ARCHIMEDLIGHT_BENCH_METAL=required \
ARCHIMEDLIGHT_LOCAL_BENCH_BUILD_SAMPLES=1 \
julia --check-bounds=auto --project=test/gpu benchmark/local_realistic_backends.jl
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

On the bundled coffee scene used for the documentation landing-page image, a
focused cached component run after toric instancing, per-direction toric
fallback, and automatic large-graph CPU propagation gave the following
representative single-sample rows on Apple Silicon:

| component | backend | median time | median allocation |
|---|---|---:|---:|
| first-order | normal CPU | 827.615 ms | 199.050 MiB |
| first-order | Raycore Metal | 347.688 ms | 207.904 MiB |
| scattering | normal CPU | 1399.963 ms | 1373.115 MiB |
| scattering | Raycore Metal | 1572.727 ms | 858.339 MiB |
| topology | normal CPU | 1295.921 ms | 1331.383 MiB |
| topology | Raycore Metal | 1566.905 ms | 832.378 MiB |
| propagation | normal CPU | 39.816 ms | 5.608 MiB |
| propagation | Raycore Metal | 38.895 ms | 5.608 MiB |
| stack trace | Raycore Metal | 90.545 ms | 0.181 MiB |

These numbers are local benchmark evidence, not a performance guarantee. They
are useful mainly as a baseline for regression work. Per-direction fallback
keeps fast Raycore traversal on safe directions and uses raster projection only
where needed for toric correctness. First-order and propagation are faster on
Metal in this run, but cached coffee scattering/topology are still slower than
normal CPU despite lower allocation. The remaining cached coffee bottleneck is
topology construction on the mixed Raycore/raster path.

The current realistic local benchmark shows the tradeoff more directly. These
rows use one build sample and one cached sample per backend, so they are quick
comparison numbers rather than stable statistical medians. The coffee rows
include the toric low-elevation projection fallback:

| workload | backend | resolved backend | build median | cached median | PAR delta | NIR delta |
|---|---|---|---:|---:|---:|---:|
| coffee multi-step | normal CPU | `RasterCPUBackend` | 2.486 s | 0.009 s | 0.000% | 0.000% |
| coffee multi-step | Raycore KA CPU | `RaycoreInterceptionBackend` | 6.603 s | 0.009 s | 0.103% | 0.089% |
| coffee multi-step | Raycore Metal | `RaycoreInterceptionBackend` | 4.301 s | 0.008 s | 2.346% | 2.654% |
| agrivoltaics wheat/panel | normal CPU | `RasterCPUBackend` | 0.161 s | 0.002 s | 0.000% | 0.000% |
| agrivoltaics wheat/panel | Raycore KA CPU | `RaycoreInterceptionBackend` | 0.205 s | 0.002 s | 0.000% | 0.000% |
| agrivoltaics wheat/panel | Raycore Metal | `RaycoreInterceptionBackend` | 0.498 s | 0.001 s | 0.000% | 0.000% |
| XPalm two oil palms | normal CPU | `RasterCPUBackend` | 1.794 s | 0.001 s | 0.000% | 0.000% |
| XPalm two oil palms | Raycore KA CPU | `RaycoreInterceptionBackend` | 4.232 s | 0.001 s | 0.000% | 0.000% |
| XPalm two oil palms | Raycore Metal | `RasterCPUBackend` fallback | 2.118 s | 0.001 s | 0.000% | 0.000% |

The breakdown mode shows where those fresh-run timings are spent:

| workload | backend | prepare | populate | reuse | interpretation |
|---|---|---:|---:|---:|---|
| coffee multi-step | normal CPU | 0.661 s / 48.2 MiB | 1.876 s / 1383.7 MiB | 0.009 s / 19.0 MiB | population still dominates the fresh CPU run |
| coffee multi-step | Raycore KA CPU | 1.008 s / 620.6 MiB | 5.579 s / 2091.6 MiB | 0.009 s / 19.0 MiB | flat stack-node population roughly halves the previous Raycore cost |
| coffee multi-step | Raycore Metal | 1.501 s / 353.7 MiB | 2.496 s / 2055.5 MiB | 0.008 s / 19.0 MiB | close to CPU populate, but still slower on fresh runs |
| agrivoltaics wheat/panel | normal CPU | 0.085 s / 15.5 MiB | 0.064 s / 12.8 MiB | 0.001 s / 3.1 MiB | both phases are small |
| agrivoltaics wheat/panel | Raycore KA CPU | 0.127 s / 122.0 MiB | 0.071 s / 12.8 MiB | 0.001 s / 3.1 MiB | Raycore setup is close to CPU |
| agrivoltaics wheat/panel | Raycore Metal | 0.148 s / 76.4 MiB | 0.375 s / 14.9 MiB | 0.001 s / 3.1 MiB | GPU population overhead dominates |
| XPalm two oil palms | normal CPU | 0.172 s / 431.2 MiB | 1.807 s / 418.3 MiB | 0.001 s / 1.6 MiB | population dominates the fresh CPU run |
| XPalm two oil palms | Raycore KA CPU | 2.348 s / 3082.2 MiB | 1.828 s / 418.3 MiB | 0.001 s / 1.6 MiB | Raycore setup dominates |
| XPalm two oil palms | Raycore Metal | 0.203 s / 396.8 MiB | 1.836 s / 418.3 MiB | 0.001 s / 1.6 MiB | Metal falls back before Raycore construction |

The practical interpretation is: cached reuse is fast on all paths once the
simulation has been prepared. Flat stack-node cache population substantially
reduced the coffee Raycore fresh run, and the coffee Metal path is now close to
the normal CPU populate phase. Normal CPU remains faster for fresh coffee runs
because Raycore still pays higher prepare cost and a larger populate phase. Use
the Raycore path when you need to exercise the portable backend or can reuse a
prepared simulation over many steps.

The XPalm row is a stress case, not a routine benchmark default. With the
default prechunk-instance cap, Metal now falls back before constructing the
roughly `7817` BLAS instances estimated for this scene with the current
`1536` face chunk limit. Disable
`RaycoreBackendConfig(max_prechunk_instances=...)` or
`ARCHIMEDLIGHT_RAYCORE_MAX_PRECHUNK_INSTANCES` only when you explicitly want
to stress-test the guarded chunked Raycore path. Always check the
resolved-backend column when comparing large local scenes.

Raycore scene setup allocates topology scratch buffers only for the reducer that
can actually be used by that scene. In particular, `edge_accumulation=:auto`
does not allocate the dense edge-count matrix for multi-instance chunked TLAS
scenes, because auto mode must use the sparse key reducer there.

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
- `ARCHIMEDLIGHT_BENCH_MAX_PRECHUNK_INSTANCES_SWEEP`
- `ARCHIMEDLIGHT_BENCH_PLOT_PAVING_SWEEP`
- `ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE_SWEEP`
- `ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION_SWEEP`

Each variable is a comma-separated list. For the bundled coffee configuration,
sector count and pixel size come from `example_2/config.yml` unless you set the
matching override or sweep variable explicitly.
`ARCHIMEDLIGHT_BENCH_MAX_PRECHUNK_INSTANCES` and its sweep form accept `0` or a
negative value to disable the prechunk-instance cap for benchmarked Raycore
backends.

`ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION` and
`ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION_SWEEP` accept `auto`,
`sparse_host_reduce`, and `dense_atomic`. The default `auto` path uses the dense
atomic reducer only when the backend reports atomic support and the dense edge
matrix fits the configured memory limit and the Raycore TLAS is not chunked into
multiple instances; otherwise it falls back to the counted sparse host reducer.
Use the sweep when testing a new GPU or scene, because the best reducer can
depend on node count, pixel count, backend atomics, and whether validation
forced BLAS chunking.

The dense atomic reducer keeps an `Int32` device scratch matrix in the cached
Raycore scene data. Dense host readback is still materialized per topology call,
because this is currently faster on Metal than copying into a retained host
buffer. A retained host aggregation vector was also tested and rejected because
it was slower on the cached coffee topology benchmark; a touched-index variant
was also slower and is not retained.

## Backend Selection

The Raycore GPU backend is most useful when:

- the same `LightSimulation` is reused for many meteo rows;
- the scene has enough projected pixels to amortize launch and transfer costs;
- scattering is enabled and the scene does not require expensive TLAS chunking;
- the benchmark uses cached execution rather than one-shot setup.

CPU can still be faster for tiny scenes, very coarse grids, one-shot calls, and
cases where most time is spent building Raycore data, validating large scenes,
chunking BLAS geometry, or converting dense internal arrays back to public
`Dict` results.

| Situation | Recommended path | Why |
|---|---|---|
| Tiny scene or one-off diagnostic run | normal CPU | Avoids Raycore setup, kernel compilation, and device transfer overhead. |
| Repeated meteo rows on the same scene | Raycore GPU with `LightSimulation` | Reuses acceleration structures, directional caches, and backend buffers. |
| Large scene that validates without chunking | Raycore GPU | Traversal and propagation work can amortize GPU launch and transfer costs. |
| Large scene that requires moderate chunking | benchmark both | Correctness is guarded, but chunked TLAS build and topology cost can dominate. |
| Large scene that exceeds the chunk-instance cap | normal CPU fallback | Avoids minutes of GPU setup for stress scenes unless the cap is disabled. |
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

The realistic local tests are opt-in because they depend on large local
fixtures:

```sh
ARCHIMEDLIGHT_TEST_GPU_LOCAL=1 \
ARCHIMEDLIGHT_TEST_METAL=required \
julia --check-bounds=auto --project=test/gpu test/gpu/runtests.jl
```

```sh
ARCHIMEDLIGHT_TEST_GPU_LARGE=1 \
ARCHIMEDLIGHT_TEST_METAL=required \
julia --check-bounds=auto --project=test/gpu test/gpu/runtests.jl
```

The local test asserts that the coffee Metal run resolves to
`RaycoreInterceptionBackend`, so a silent raster fallback fails the test. The
large test defaults to the local agrivoltaics wheat-panel scene and can be
configured to use an XPalm/VPalm-generated scene when that checkout is
instantiated locally.
