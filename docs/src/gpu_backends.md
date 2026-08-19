```@meta
CurrentModule = ArchimedLight
```

# GPU Backends And Benchmarks

ArchimedLight provides one experimental GPU implementation:
[`RasterGPUBackend`](@ref). It is written with KernelAbstractions and can target
Metal, CUDA, oneAPI, AMDGPU, or the KernelAbstractions CPU backend supplied by
the calling environment.

GPU packages are optional. ArchimedLight depends on KernelAbstractions, but it
does not make Metal or another vendor runtime a required package dependency.

## Minimal Setup

Construct the array backend in the package that owns the device, then pass it to
ArchimedLight:

```julia
using ArchimedLight
using KernelAbstractions
using Metal

device = KernelAbstractions.get_backend(MtlArray(zeros(Float32, 1)))
interception = RasterGPUBackend(
    backend=device,
    max_hits_per_pixel=128,
    tile_size=1,
    tile_face_capacity=128,
    edge_accumulation=:auto,
)
scattering = RasterGPUScatteringBackend(interception)

sim, meteo = read_simulation(
    "config.yml";
    interception_backend=interception,
    scattering_backend=scattering,
)
step = run_light(sim, first(meteo))
```

Replace `MtlArray` with the corresponding array type for another
KernelAbstractions backend.

## Configuration

`max_hits_per_pixel` bounds the full visibility stack retained for one raster
cell. `tile_size` and `tile_face_capacity` control face binning. A tile overflow
or hit-stack overflow is reported as an error rather than silently changing the
result.

`edge_accumulation` accepts:

- `:auto`, which selects dense atomics when supported and affordable;
- `:dense_atomic`, which requires backend atomic support and a dense edge matrix
  within `dense_edge_limit_bytes`;
- `:sparse_host_reduce`, which returns counted packed keys for host reduction.

The dense kernels use `@atomic :monotonic`. They require working atomics from
the device package and KernelAbstractions, but they do not rely on
acquire/release ordering.

## Metal Environment

Metal 1.10.3 is the minimum supported release for the required atomic behavior.
That minimum compatibility is declared only in the dedicated GPU environment:

```sh
julia --project=test/gpu -e 'using Pkg; Pkg.instantiate()'
ARCHIMEDLIGHT_TEST_METAL=required julia --project=test/gpu test/gpu/runtests.jl
```

The focused Metal suite checks the installed version and atomic capability,
then compares RasterGPU first-order interception and scattering against the CPU
raster reference. Use `ARCHIMEDLIGHT_TEST_METAL=1` for an optional local run or
`required` when unavailable Metal hardware must fail validation.

Metal.jl 1.10.3 reports KernelAbstractions atomic support only when the active
toolchain exposes MSL 4.1 or newer. The package version supplies the registered
implementation, but it cannot upgrade the system Metal compiler; the strict
suite fails before launching kernels when `supports_atomics(device)` is false.

## Validation

Use `validate=true` while bringing up a new device backend. The native raster
path compares device projections against the CPU raster reference and reports a
failure instead of switching backends.

The normal CPU-focused suite also exercises RasterGPU on
`KernelAbstractions.CPU()`:

```sh
julia --project=test test/runtests.jl
```

## Benchmarks

The repository retains general CPU/RasterGPU benchmarks:

```sh
julia --project=benchmark benchmark/backend_comparison.jl
julia --project=benchmark benchmark/rastergpu_first_order_reuse.jl
julia --project=benchmark benchmark/artifact_scene_matrix.jl
```

The comparison scripts accept environment variables for backend selection,
sample counts, workgroup size, hit capacity, raster tile sizes, edge
accumulation, and optional TSV output. Always record the resolved backend and
compare scientific totals with the CPU reference when evaluating a new device.

## Backend Selection

Prefer `RasterCPUBackend` for small one-shot workloads and reference results.
Use `RasterGPUBackend` when the projected raster and scene complexity are large
enough to amortize device setup, especially when a `LightSimulation` is reused
across multiple time steps.
