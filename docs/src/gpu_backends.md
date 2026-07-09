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

## Validation

Raycore backend selection is explicit. If you pass a non-CPU Raycore backend,
ArchimedLight runs that Raycore backend directly. It does not validate against
the CPU reference by default. Choose `RasterCPUBackend()` or `:raster_cpu` explicitly
when you want CPU behavior.

Set `RaycoreBackendConfig(validate=true)` or pass `validate=true` through
`RaycoreInterceptionBackend(...)` to enable runtime validation. With validation
enabled, ArchimedLight compares the Raycore traces against the CPU reference
before cached execution and scattering topology construction. A mismatch throws
`RaycoreValidationError` with a `reason`, `stage`, message, and raw validation
payload.

This guard is mainly there for debugging purposes.

Very large chunked scenes can be correct but impractical to build on a GPU.
ArchimedLight throws `RaycoreValidationError` before Raycore construction when
prechunking would create more than
`RaycoreBackendConfig(max_prechunk_instances=...)` allows. The default comes
from `ARCHIMEDLIGHT_RAYCORE_MAX_PRECHUNK_INSTANCES` (`4096` when unset). Pass
`0` or a negative value, either in the config or in the environment variable,
to disable the cap and force Raycore construction.

If validation is enabled and a prechunked non-CPU Raycore scene fails
validation, ArchimedLight retries with smaller chunk limits before throwing,
while still respecting `max_prechunk_instances`. This handles occasional
large-TLAS validation instability without changing the selected backend.

For cached simulations and Raycore scattering graph construction,
`validate=true` also runs a sampled full-stack validation against a
KernelAbstractions CPU Raycore reference. This catches cases where top-hit
vertical tracing succeeds but `Raycore.all_hits!` misses deeper stack hits on a
GPU backend. If the sampled full-stack hit or occupied-pixel ratios are too
low, ArchimedLight throws `RaycoreValidationError` with
`reason = :raycore_stack_trace_validation`.

The full-stack guard defaults to three sampled directions, at least 512
reference stack hits or 128 occupied reference pixels per sampled direction,
and minimum aggregate hit/occupied ratios of `0.95`. These defaults can be
tuned for diagnostic sweeps with:
`ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_DIRECTIONS`,
`ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_MIN_HITS`,
`ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_MIN_OCCUPIED`,
`ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_MIN_HIT_RATIO`, and
`ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_MIN_OCCUPIED_RATIO`. Keep the defaults
for production comparisons unless a benchmark row explicitly records why a
different threshold was used.

The Metal raw-stack diagnostics isolated a Raycore GPU acceleration-structure
construction problem before any ArchimedLight reduction. With Raycore building
the coffee BLAS/TLAS directly on Metal, reduced unchunked coffee remained a
non-overflow `candidate_missing_hits` failure even after raising Raycore's
software traversal stack capacity from 32 to 64. Increasing
`max_hits_per_pixel` from 32 to 256 did not recover the missing hits, neither
CPU nor Metal reported hit-buffer overflow, and disabling duplicate suppression
with `hit_epsilon=-1` still failed. The first failing ray had identity toric
transform, valid CPU-side TLAS-leaf and BLAS-root intersections, and valid
candidate-side bounds copied back from Metal. That ruled out ArchimedLight
fallback, host reduction, hit-buffer size, duplicate suppression, metadata
decoding, toric transform layout, and the old 32-entry traversal capacity.

For correctness, non-CPU Raycore scene construction now builds the mutable
Raycore TLAS on KernelAbstractions CPU and adapts the resulting `StaticTLAS`
arrays to the requested traversal backend. The selected backend still runs the
Raycore tracing kernels and owns the trace/decoder/reduction buffers; only the
acceleration-structure build is kept on the deterministic CPU path. The focused
Metal test writes the current rows to
`/tmp/archimed_raycore_metal_raw_stack_smoke.tsv`. A representative run with
the CPU-built static TLAS showed exact raw-stack parity for reduced coffee at
`max_hits=32,64,128,256`, exact parity with duplicate suppression disabled,
exact parity for forced chunk limits `1536,1024,768,512`, and exact parity for
the reduced agrivoltaics row. The strict scene-complexity smoke also stayed on
`RaycoreInterceptionBackend` for simple plant, full bundled coffee, and the
reduced agrivoltaics scene.

Production non-CPU Raycore scene construction now checks the built BLAS depth
before validation. If it exceeds the known software traversal-stack capacity and
an allowed `face_chunk_limit` can reduce the depth, ArchimedLight rebuilds the
scene as `geometry_mode=:chunked_merged_mesh` proactively; validation retries
with smaller chunks remain as a fallback if a later trace check still fails. On
the old 32-entry Raycore traversal stack, the same local coffee case needed
chunking: `face_chunk_limit=1536` improved the sampled raw hit ratio to roughly
0.95, while `face_chunk_limit=1024` passed the default full-stack guard with a
hit ratio of about 0.996 and occupied-pixel ratio 1.0. With a Raycore build that
reports a larger software traversal capacity, ArchimedLight keeps the
unchunked merged BLAS when its depth is within capacity, and the strict cached
Metal run must still remain on `RaycoreInterceptionBackend` rather than
silently resolving to `RasterCPUBackend`.

Benchmark configs enable validation by default so correctness regressions fail
clearly. Set `ARCHIMEDLIGHT_BENCH_VALIDATE=0` only when intentionally measuring
direct Raycore runtime without CPU reference validation. A row labelled
`raycore_metal_gpu` with resolved backend `RaycoreInterceptionBackend` stayed on
Raycore.

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
to CPU. Scalar sector cache fills also use an internal dense-only propagation
entry point when the transfer graph and prepared geometry share node order. In
that path the cache builds the unit-sector dense initial-power vector directly,
so neither `FirstOrderResult` nor public `Dict` results are built until an API
result actually needs them.

Internally, the full-response cache now keeps a backend-neutral dense response
storage instead of materializing per-sector `Dict`-style objects early. The
cache stores a node-by-sector projected-area matrix, compact active row-index
arrays for nonzero sector entries, aggregate hit counts by node, dense emitter
incident-power vectors, and a compact scattering topology. The topology also
caches dense edge endpoint indices and sun-hit counts in the same node order, so
scattering propagation does not have to remap node ids back to dense indices on
each transfer-graph build. Public `Dict` outputs are materialized only at API
boundaries such as `FirstOrderResult`, `ScatteringResult`, and light-budget
integration.

The dense storage also carries a device-cache slot for future backend-resident
copies. Today, Raycore/KernelAbstractions still performs several reductions and
public-output conversions on the host: compact sector active-index extraction,
emitter incident-power accumulation, compact edge-count materialization from
traced stacks, some scattering transfer-graph construction, and final public
`Dict` construction. Dense projection paths keep sun-hit accumulation in a
dense node-order vector and merge any non-dense fallback dictionary counts only
once when the topology cache is built. Raycore full-stack paths increment this
dense vector directly while walking decoded stack nodes, avoiding a separate
per-sector sun-hit scratch vector and merge scan. Dense response builders also
accumulate aggregate hit counts directly into the retained all-sector hit vector
instead of copying each sector into a temporary hit vector first. Raycore
top-hit first-order interception likewise accumulates directly into dense node
vectors when area-ratio correction is disabled; the per-sector scratch vectors
are retained only for the area-ratio path that needs node-local correction
denominators. Flat full-stack cache and streaming paths also keep projected-hit
denominator storage area-ratio-only, and the streaming path folds first-order
area accumulation into the first decoded-stack walk when no area-ratio
correction is needed. The parts already structured for future device-side
reductions are the projected-area matrix, aggregate hit vector, compact
active-index offsets, dense sun-hit vector, decoded stack-node buffers,
edge-key or dense-edge scratch buffers, and compact
`edge_to`/`edge_from`/`edge_count` topology arrays.

Dense response cache population now stores each sector's projected-area vector
and records its compact active row indices in the same host pass. That still
materializes the dense node-by-sector matrix on the host, but avoids a separate
copy-plus-rescan phase before public outputs or cached scattering inputs consume
the response cache.

Scattering topology cache construction also reuses the prepared dense
node-index map when converting packed edge counts, local sun-hit dictionaries,
and `ScatteringPairCounts` endpoints into dense index arrays. The compatibility
helpers can still build a temporary map for standalone/test callers, but normal
prepared CPU and Raycore paths no longer rebuild that dictionary at each
topology-cache boundary.

When a cached topology is converted into a scattering transfer graph with the
same dense node order, ArchimedLight now builds the public coefficient maps and
the dense PAR/NIR coefficient vectors in one pass. The dense static graph then
uses those vectors directly, so it does not need to remap the coefficient
dictionaries back into node-order arrays. The public `ScatteringTransferGraph`
fields are still populated for API compatibility and fallback code paths.

For Raycore device propagation, the dense static graph now also caches
device-resident PAR/NIR coefficient vectors for the default graph coefficients.
Repeated unit-sector or band propagations can therefore reuse the backend copies
of edge arrays, all-hit counts, and coefficients; only the dynamic initial
energy vector and iteration scratch buffers are refreshed for each solve.

Repeated device scattering solves also reuse device-resident copies of the
static edge arrays and dense all-hit normalization vector stored behind the
internal dense scattering graph. This avoids recopying graph-normalization data
when the same cached transfer graph is propagated for multiple bands or unit
sector inputs.

The combined Raycore first-order/topology path keeps sparse edge accumulation on
the host for CPU and current Metal runs. Device sparse edge-key extraction is
available internally, but coffee benchmarks showed that reading back and sorting
those sparse keys is slower than the fused host stack walk on Metal. The flat
path therefore only routes edge accumulation through device reductions when a
backend can use dense atomic edge counts.

The same flat path now also keeps a retained dense node-count scratch buffer in
`RaycoreSceneData` for non-CPU backends with atomic support. After a full-stack
trace, a small KernelAbstractions reducer accumulates decoded stack hits into a
dense node-ordered `Int32` vector on the selected backend. The host copies back
that dense node vector for aggregate hit counts, sun-hit counts, and the
projected-pixel denominator used by `area_ratio=true`, instead of walking every
hit for those reductions on CPU. KernelAbstractions CPU and backends without
atomics keep the previous host stack walk.

The flat path also has a retained dense projected-area scratch vector for the
same class of non-CPU atomic backends. Its reducer walks each decoded stack on
the selected backend, applies ArchimedLight's transparent-stack semantics, and
emits dense node-ordered projected area. In the cached response path this fills
sector projected area directly; in the streaming first-order path the host uses
that dense vector to update projected area and PAR/NIR incident power. This is
still guarded: KernelAbstractions CPU, backends without atomics,
raster-compatible fallback sectors, and any failed device reducer keep the host
stack walk.

When a sector needs both hit counts and projected area, ArchimedLight now uses a
fused device reducer where available. The fused kernel walks each decoded stack
once, increments dense hit counts for every hit, and accumulates projected area
only while the transparent-stack rule still lets light continue. If the fused
kernel cannot run, the path falls back to the separate count and area reducers,
then to the host stack walk as before.

The dense count, projected-area, and dense edge-count reducers copy their device
outputs into retained host buffers owned by `RaycoreSceneData`. This keeps the
public output semantics unchanged, but avoids allocating fresh host vectors for
every sector reduction and gives the CPU, KernelAbstractions CPU, and future GPU
paths the same dense node-order response buffers to consume. The dense
edge-count buffer is still materialized on the host when building the current
public scattering topology; future active-edge compaction can replace that
copyback without changing the surrounding cache representation.

The flat Raycore cache builder consumes these dense reducer outputs through a
single internal dense-response handoff. That handoff updates aggregate node hit
counts, sun-hit counts, projected-area sector columns, and active projected-area
indices before any public `Dict` materialization. This keeps the current
host-side cache representation stable while providing one boundary that future
device-side compaction or fused sector reducers can target. The streaming
first-order plus scattering-topology path uses the same style of dense
reducer-consumption boundary for aggregate first-order projected area,
PAR/NIR incident power, node hit counts, and sun-hit counts, including the
two-step area-ratio case where projected-pixel denominators are collected before
sector area is scaled.

These reducers are intentionally incremental. The flat path still reads raw
decoded stacks back to the host for fallback edge accumulation and for any
backend/sector where dense device outputs do not cover the remaining semantics.
Compact active-index extraction remains host-side. On the current local Metal
stack, `KernelAbstractions.supports_atomics` does not expose a usable path for
these count or projected-area reducers, so Metal uses the guarded host fallback
for them.

When the dense count reducer, dense projected-area reducer, and dense atomic
edge reducer are all available for a sector, the flat path now skips
`_raycore_copy_direction_stack_nodes_host!` for that sector. It copies back only
the dense node vectors and dense edge matrix needed to populate the existing
cache structures. If any reducer is unavailable, the previous raw stack readback
and host walk remain the fallback.

For device-produced projected-area vectors, the cached response builder stores
the dense sector column and appends active projected-area node indices in one
pass over the copied dense vector. This avoids bouncing the device output
through the temporary `sector_area` vector only to scan it again for active
rows. Full device-side prefix-sum compaction is still left for a backend where
dense node-vector copyback or host active-index extraction shows up clearly in
benchmarks.

When the flat path uses the host stack walk, it accumulates scattering edges as
dense `(to_idx, from_idx)` keys in the prepared node order and only materializes
node-id pair arrays when the topology cache is built. That keeps the public
topology and `Dict` outputs unchanged while avoiding per-edge node-id packing in
the hot loop.

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
In cached mode, the table also reports the resolved backend; Raycore validation
errors are reported as failing rows rather than CPU fallback timings. By
default, benchmark Raycore configs set `validate=true`; set
`ARCHIMEDLIGHT_BENCH_VALIDATE=0` to skip CPU reference validation for a direct
backend timing.
The table also reports a compact `reductions` capability label for Raycore
rows:

- `atomics`: whether KernelAbstractions reports backend atomic support.
- `edge`: whether dense device edge accumulation is active.
- `count`: whether dense node-count reduction is active.
- `area`: whether dense projected-area reduction is active.
- `fused`: whether the fused count/projected-area reducer is active.

For example, `atomics=0,edge=0,count=0,area=0,fused=0` means the row used
Raycore traversal but kept these dense reductions on the guarded host-stack
fallback path. This is the current local Metal behavior on the tested
KernelAbstractions/Metal stack. CUDA, oneAPI, or AMDGPU runs should be checked
through this column before interpreting a GPU timing as a device-reduction
timing.

Use `benchmark/scene_parameter_sweep.jl` when the question is how pixel size
and turtle-sector count interact across several scene shapes. It has generated
`single_plate`, `stacked_plates`, `canopy_grid`, and `panel_canopy` scenes by
default, with optional `coffee_docs`. Configure it with:

- `ARCHIMEDLIGHT_PARAM_SCENES=single_plate,stacked_plates,canopy_grid,panel_canopy,coffee_docs`
- `ARCHIMEDLIGHT_PARAM_PIXELS=0.20,0.10` and `ARCHIMEDLIGHT_PARAM_SECTORS=6,16,46` for global sweeps
- `ARCHIMEDLIGHT_PARAM_CANOPY_GRID_PIXELS=0.12,0.06` and `ARCHIMEDLIGHT_PARAM_CANOPY_GRID_SECTORS=6,16,46` for per-scene overrides
- `ARCHIMEDLIGHT_PARAM_BACKENDS=normal_cpu,raycore_ka_cpu,raycore_metal_gpu`
- `ARCHIMEDLIGHT_PARAM_VALIDATE=1` to keep explicit Raycore validation enabled
- `ARCHIMEDLIGHT_PARAM_OUTPUT=/tmp/archimed_scene_parameter_sweep.tsv`

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
ARCHIMEDLIGHT_BENCH_COMPONENTS=first_order,scattering,topology,propagation,integration,stack_trace,stack_profile \
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
cached-reuse timing, light totals, the resolved backend, and the same
`reductions` capability label used by the component benchmark. Use
`ARCHIMEDLIGHT_BENCH_BACKENDS` and the same optional GPU flags as
`benchmark/backend_comparison.jl` to select Metal, CUDA, oneAPI, or AMDGPU. Set
`ARCHIMEDLIGHT_BENCH_MAX_PRECHUNK_INSTANCES` to override the Raycore
`max_prechunk_instances` cap for these benchmarked backends. Set
`ARCHIMEDLIGHT_BENCH_VALIDATE=0` to skip CPU reference validation for direct
Raycore timing rows. Set
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
top-hit validation/retry, and full-stack validation/retry. The printed row and
TSV output report both the vertical validation ratio and the full-stack
`stack_hit`/`stack_occ` ratios. In strict benchmark mode, a full-stack validation
failure throws the same `RaycoreValidationError` that cached production
preparation would throw instead of reporting `resolved=RasterCPUBackend`. On the
bundled coffee scene with Metal, this mode showed that the Raycore scene/TLAS
phase dominates prepare time. It also showed that
`ARCHIMEDLIGHT_RAYCORE_FACE_CHUNK_LIMIT=2048` fails validation
(`ratio=0.759`), while the current default `1536` passes (`ratio=0.992`) with
721 chunked BLAS instances.

Set `ARCHIMEDLIGHT_LOCAL_BENCH_STACK_PROFILE=1` to profile Raycore stack
tracing, guarded count/projected-area reduction, dense edge reduction, and raw
stack host copy over the positive turtle directions that would populate the
runtime full-response cache. By default,
`ARCHIMEDLIGHT_LOCAL_BENCH_PROFILE_UNIQUE_TURTLES=1`, so repeated meteo rows
with the same turtle directions are collapsed with the same internal
`_turtle_cache_key` used by `LightSimulation`. This keeps the profile aligned
with cached runtime behavior. Set it to `0` only when you explicitly want to
profile every selected workload timestep. Topology-fallback rows keep repeated
turtles because that path does not reuse full-response turtle entries. This is
the local-realistic counterpart to the
`stack_profile` component in `benchmark/backend_comparison.jl`, and is the
right mode for coffee, agrivoltaics, and XPalm before deciding whether to invest
in backend-specific readback removal or device-side compaction. The table also
prints per-direction trace/copy times plus stack-density metrics:
`hit_util`, `occupied`, and `max`. The `turtles` field is
`input/profiled`, so `3/1` means three selected meteo rows reused one cached
turtle response. The `dirs` field is
`traced/reduced/edge-reduced/diagnostic-copied/production-copy-required`; the
diagnostic copy is always performed so density metrics can be reported, while
the final count shows whether the production flat path would actually need raw
stack readback under the same reducer assumptions. These fields help
distinguish trace cost caused by many directions from trace cost caused by
dense/deep stacks, overflow risk, or unavoidable host readback.

Set `ARCHIMEDLIGHT_LOCAL_BENCH_OUTPUT=/path/to/results.tsv` to write the
local-realistic benchmark rows to a TSV file in addition to printing them. This
works for the default build/reuse comparison, cache breakdown, prepare
breakdown, and stack-profile modes. Use it to keep capability-labelled rows
from different machines together, especially when comparing Metal, CUDA,
oneAPI, AMDGPU, and KernelAbstractions CPU results. TSV output includes run
metadata columns such as benchmark mode, Julia version, OS/architecture,
backend/workload selectors, `max_hits_per_pixel`, workgroup size, prechunk cap,
reference-instancing setting, and full-stack validation thresholds. Set
`ARCHIMEDLIGHT_LOCAL_BENCH_RUN_LABEL=...` to add a short machine or CI label to
each saved row.

Summarize one or more saved TSV files with:

```sh
julia --project=test/gpu benchmark/summarize_capability_matrix.jl \
  /path/to/metal.tsv /path/to/cuda.tsv /path/to/amdgpu.tsv
```

The summarizer accepts both old and metadata-prefixed TSV rows, including rows
that use either the local realistic benchmark's stack-profile column names or
the focused component benchmark's `profile_*`/`trace_*` fields. It prints the
resolved backend, reducer capability label, stack-profile direction counts,
copy-required count, hit-density metrics, and a compact status such as
`host-readback-required` or `readback-free-candidate`. When stage-split rows are
present, it also prints `prepare` and `public` seconds so normal CPU and Raycore
rows can be compared without opening the raw TSV.

Set `ARCHIMEDLIGHT_LOCAL_BENCH_SELFTEST=1` for a private benchmark-driver smoke
test. It builds a tiny in-memory plate scene and checks that the stack-profile
collector handles fallback Raycore data, empty turtle lists, no-positive-sector
turtles, one positive sector, cached repeated-turtle deduplication, and
topology-fallback repeated-turtle profiling without loading the local coffee,
agrivoltaics, or XPalm workloads.

Set `ARCHIMEDLIGHT_LOCAL_BENCH_STAGE_SPLIT=1` when deciding what to optimize
next. This mode writes one row per workload/backend with cache preparation time,
Raycore stack trace time, device count/area reduction time, device edge
reduction time, host stack-copy time, and a final cached public API run. The
public phase runs after cache preparation and stack profiling, so it includes
response combination, scattering/integration, and final public result
materialization. Use the TSV columns `profile_trace_ms`,
`profile_count_area_ms`, `profile_edge_ms`, `profile_copy_ms`, and
`public_output_seconds` to decide whether the next bottleneck is traversal,
device reduction, host readback, or API-boundary materialization.
Set `ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION=sparse_host_reduce` or
`dense_atomic` to compare merged-mesh edge accumulation paths in the same local
realistic benchmark. The capability summarizer prints `count` and `edge`
milliseconds so this comparison is visible without opening the raw TSV.
Use `ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION_SWEEP=auto,sparse_host_reduce,dense_atomic`
to produce all edge-mode rows in one run; `normal_cpu` is emitted only once
during the sweep to avoid duplicate CPU baselines.

On reduced local stage-split smoke settings for coffee, agrivoltaics, and
XPalm, the benchmark rows are intended to choose the next optimization target,
not to replace the full realistic benchmark tables. Raycore KA CPU can now run
the dense count/area and edge reducers when atomics are available, so the
reduced coffee, reference-instanced agrivoltaics, and reference-instanced XPalm
rows classify as `readback-free-candidate`. The remaining bottleneck in these
KA CPU rows is not production host readback: traversal dominates coffee and
XPalm, while coffee also pays a visible dense edge-reduction cost.
Stage-split TSV rows also include `ref_instancing_*` diagnostics. On the
reduced coffee row, `ref_instancing_status=no_reusable_references` because the
scene has many tapered nodes and the remaining untapered candidate meshes are
unique, so the Raycore path correctly stays on `merged_mesh`.
They also include `raycore_*` scene-shape diagnostics from the actual built
Raycore data when available: TLAS instance/geometries counts, dense node count,
expanded face count, an expanded-face instance upper bound, reference-prototype
compact face counts, dense edge-matrix pairs, and sparse edge-key capacity. The
capability summarizer prints compact `tlas`, `faces`, and `scratch` labels from
these columns so traversal-heavy rows can be compared without inspecting the raw
TSV by hand.
On the reduced coffee row, `dense_atomic` removes production copyback but adds a
visible dense edge-reduction phase on KA CPU; `sparse_host_reduce` avoids that
phase but remains `host-readback-required`. Keep both rows when comparing CPU
and GPU hardware because the better choice may differ by backend.

Reduced coffee edge-mode sweep for an atomic-capable backend:

```sh
ARCHIMEDLIGHT_LOCAL_BENCH_STAGE_SPLIT=1 \
ARCHIMEDLIGHT_LOCAL_BENCH_WORKLOADS=coffee \
ARCHIMEDLIGHT_BENCH_BACKENDS=raycore_cuda_gpu \
ARCHIMEDLIGHT_BENCH_CUDA=required \
ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION_SWEEP=auto,sparse_host_reduce,dense_atomic \
ARCHIMEDLIGHT_LOCAL_BENCH_BUILD_SAMPLES=1 \
ARCHIMEDLIGHT_BENCH_COFFEE_PAVING=100 \
ARCHIMEDLIGHT_LOCAL_BENCH_OUTPUT=/tmp/archimed_stage_split_coffee_cuda_edge_sweep.tsv \
julia --project=test/gpu benchmark/local_realistic_backends.jl
```

Replace the CUDA backend and flag with `raycore_oneapi_gpu` /
`ARCHIMEDLIGHT_BENCH_ONEAPI=required` or `raycore_amdgpu_gpu` /
`ARCHIMEDLIGHT_BENCH_AMDGPU=required` on the matching hardware. Summarize the
TSV with:

```sh
julia --project=. benchmark/summarize_capability_matrix.jl \
  /tmp/archimed_stage_split_coffee_cuda_edge_sweep.tsv
```

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
  with structured result fields for fixed hit-stack utilization, maximum
  observed stack depth, occupied pixels, and overflow. The focused component
  benchmark table also reports the resolved backend and Raycore geometry mode,
  so a requested GPU row that actually measured fallback or merged-mesh TLAS is
  visible without parsing the workload label.
- `stack_profile`: Raycore full-stack traversal plus the guarded internal
  reduction/readback phases. Its rows carry structured fields such as
  `profile_trace_ms`, `profile_count_area_ms`, `profile_edge_ms`,
  `profile_copy_ms`, `profile_copy_required_dirs`, and
  `profile_reduced_dirs`. Use this before optimizing readback or active-index
  compaction: if `reductions` reports disabled reducers, `profile_copy_ms`
  measures the raw stack readback fallback; if dense reducers are active, the
  count/area and edge timings show whether device reduction or copyback is the
  bottleneck.

On the bundled coffee scene used for the documentation landing-page image, a
focused cached component run after toric instancing, per-direction toric
fallback, and automatic large-graph CPU propagation gave the following
representative rows on Apple Silicon. The latest node-count reducer benchmark
used one warmup and three samples:

| component | backend | median time | median allocation |
|---|---|---:|---:|
| first-order | normal CPU | 823.341 ms | 202.513 MiB |
| first-order | Raycore KA CPU | 2026.131 ms | 214.228 MiB |
| first-order | Raycore Metal | 201.439 ms | 209.752 MiB |
| scattering | normal CPU | 1331.054 ms | 1346.502 MiB |
| scattering | Raycore KA CPU | 4486.512 ms | 929.200 MiB |
| scattering | Raycore Metal | 1330.766 ms | 907.180 MiB |
| topology | normal CPU | 1302.888 ms | 1379.054 MiB |
| topology | Raycore KA CPU | 4464.652 ms | 919.049 MiB |
| topology | Raycore Metal | 1298.494 ms | 882.383 MiB |
| propagation | normal CPU | 39.402 ms | 4.559 MiB |
| propagation | Raycore KA CPU | 39.227 ms | 4.559 MiB |
| propagation | Raycore Metal | 38.467 ms | 4.559 MiB |

These numbers are local benchmark evidence, not a performance guarantee. They
are useful mainly as a baseline for regression work. Per-direction fallback
keeps fast Raycore traversal on safe directions and uses raster projection only
where needed for toric correctness. First-order, scattering, topology, and
propagation were faster on Metal than normal CPU in this run, while Raycore KA
CPU remained slower. The remaining cached coffee bottleneck is stack-to-response
compression on the mixed Raycore/raster path, especially the portions that still
read raw decoded stacks on the host.

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
| XPalm two oil palms | Raycore Metal | validation error before Raycore construction | n/a | n/a | n/a | n/a |

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
| XPalm two oil palms | Raycore Metal | n/a | n/a | n/a | prechunk-instance cap throws before Raycore construction |

The practical interpretation is: cached reuse is fast on all paths once the
simulation has been prepared. Flat stack-node cache population substantially
reduced the coffee Raycore fresh run, and the coffee Metal path is now close to
the normal CPU populate phase. Normal CPU remains faster for fresh coffee runs
because Raycore still pays higher prepare cost and a larger populate phase. Use
the Raycore path when you need to exercise the portable backend or can reuse a
prepared simulation over many steps.

The XPalm row is a stress case, not a routine benchmark default. With the
default prechunk-instance cap, Metal now throws before constructing the roughly
`7817` BLAS instances estimated for this scene with the current `1536` face
chunk limit. Disable
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

The aggregate/sun/projected-pixel count reducer and projected-area reducer keep
separate retained scratch vectors in the same cached Raycore scene data. They
are enabled only for non-CPU backends that report atomic support. When active,
they copy back dense `n_nodes` vectors rather than revisiting the raw stack for
those quantities. Raw stack readback remains in the current implementation for
edge accumulation paths that cannot use dense device atomics, for compact
active-index construction, and for guarded fallback sectors.

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
| Large scene that exceeds the chunk-instance cap | choose CPU explicitly or raise the cap | Raycore throws instead of silently switching backend. |
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

Run a raw stack comparison for the hit-buffer and KA-vs-Metal hypotheses:

```sh
ARCHIMEDLIGHT_RAW_STACK_BACKEND=metal \
ARCHIMEDLIGHT_RAW_STACK_WORKLOAD=coffee \
ARCHIMEDLIGHT_RAW_STACK_MAX_HITS_SWEEP=32,64,128,256 \
ARCHIMEDLIGHT_RAW_STACK_FACE_CHUNK_LIMITS=auto,unchunked,1536,1024,768,512 \
ARCHIMEDLIGHT_RAW_STACK_OUTPUT=/tmp/archimed_raw_stack_coffee_metal.tsv \
julia --check-bounds=auto --project=test/gpu scripts/raycore_raw_stack_compare.jl
```

The raw stack diagnostic compares the same prepared scene on Raycore
KernelAbstractions CPU and the requested backend. It reports total hit counts,
overflow counts, hit/occupied-pixel ratios, maximum per-pixel stack depth, a
diagnosis such as `candidate_missing_hits`, `candidate_extra_hits`,
`raw_order_difference`, `instance_index_difference`, or `distance_difference`,
and the first mismatching raw `(metadata, instance_idx, distance)` slot with the
pixel ray. Set `ARCHIMEDLIGHT_RAW_STACK_WORKLOAD=simple`, `coffee`, `agripv`,
or `synthetic_stack`; set `ARCHIMEDLIGHT_RAW_STACK_SECTOR=all` to scan all
non-sun upward sectors. `synthetic_stack` builds one repo-local stacked-plate
mesh in a single BLAS, controlled by `ARCHIMEDLIGHT_RAW_STACK_SYNTHETIC_LAYERS`
and `ARCHIMEDLIGHT_RAW_STACK_SYNTHETIC_SPACING`, and is useful for probing
deep-BLAS behavior without loading coffee or the local agrivoltaics OPF. Set
`ARCHIMEDLIGHT_RAW_STACK_FACE_CHUNK_LIMITS` to `auto` or a
comma-separated list such as `auto,unchunked,1536,1024`. `auto` uses the same
production scene-data selection, including proactive chunking when BLAS depth
would exceed Raycore's software traversal stack. `unchunked` forces the merged
BLAS path and is the mode to use when checking whether a Raycore traversal
capacity change fixes or reproduces a deep-BLAS Metal missing-hit failure.
Integer values force that face chunk limit. The output
records the requested chunk mode, the CPU reference chunk mode selected to match
the candidate topology, whether that topology matched, resolved Raycore geometry
mode, TLAS instance/geometry counts, TLAS and BLAS depth, BLAS node count,
expanded face counts, unmatched reference and candidate hits for the first
mismatching pixel, and the first raw mismatch. This makes the coffee reproducer
usable for distinguishing hit-buffer capacity, unchunked TLAS traversal,
production-auto chunking, forced chunked BLAS geometry effects, and whether
Metal is missing hits or merely returning the same hit set in a different order.

The optional Metal GPU test suite compares normal CPU, Raycore
KernelAbstractions CPU, and Raycore Metal on the same simple scene with
validation enabled for both Raycore backends. It asserts that both Raycore
backends remain resolved as `RaycoreInterceptionBackend` and that aggregate
PAR/NIR incident and absorbed totals match the normal CPU oracle within
`rtol=5e-3, atol=1e-3`. The suite also runs the raw stack sweep for simple and
coffee by default. The current coffee reproducer must be exact for the
`max_hits=32,64,128,256` sweep, for duplicate suppression disabled
(`hit_epsilon=-1`), and for the forced chunk limits used to stress the
prechunked path. The local agrivoltaics raw-stack comparison is slower because
it has to load the local OPF and build the wheat/panel scene. When the OPF is
available, the Metal suite now runs a reduced agrivoltaics raw-stack comparison
by default; set `ARCHIMEDLIGHT_TEST_GPU_AGRIPV_RAW=0` to skip it, or set it to
`required` to fail if the OPF is missing. The test requires the CPU reference
and Metal candidate to use matching Raycore topology and requires any non-exact
row to carry the same structured raw-hit diagnosis used by the coffee
reproducer. Use
`ARCHIMEDLIGHT_RAW_STACK_AGRIPV_ROWS`,
`ARCHIMEDLIGHT_RAW_STACK_AGRIPV_PLANTS_PER_ROW`,
`ARCHIMEDLIGHT_RAW_STACK_AGRIPV_GROUND_NX`, and
`ARCHIMEDLIGHT_RAW_STACK_AGRIPV_GROUND_NY` to keep that raw diagnostic reduced
before scaling it back toward the full local agrivoltaics scene. The focused
Metal test writes the coffee and reduced agrivoltaics raw-stack rows to
`/tmp/archimed_raycore_metal_raw_stack_smoke.tsv` so exact rows and diagnosed
non-exact rows can be inspected after the test run. The default Metal suite
still runs a reduced agrivoltaics scene-complexity benchmark with
`ARCHIMEDLIGHT_BENCH_VALIDATE=1` and asserts that the requested Metal
Raycore backend remains the resolved backend.

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
strict scene-complexity smoke also asserts `status = "ok"` for the requested
Metal rows. That smoke covers the simple plant, the bundled coffee scene with
full `plot_paving_override=2500` at a coarse `pixel_size=0.04`, and the reduced
agrivoltaics wheat-panel scene when the local OPF is available. The focused
Metal test writes those strict rows to
`/tmp/archimed_raycore_metal_strict_scene_complexity_smoke.tsv` so the resolved
backend, geometry mode, and Raycore TLAS shape can be audited after the test
run. If a strict scene-complexity row fails before it can run, the TSV row uses
`status = "error"` and records `error_reason`, `error_stage`, `error_message`,
and the available full-stack validation counts/ratios instead of reporting
`resolved_backend = RasterCPUBackend`. The large test defaults to the local
agrivoltaics wheat-panel scene and can be configured to use an
XPalm/VPalm-generated scene when that checkout is instantiated locally.
