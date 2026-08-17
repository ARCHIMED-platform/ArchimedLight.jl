# Voxel CPU baseline — 2026-08-17

This is the first warmed performance baseline for the voxel implementation.
It was recorded from the uncommitted `voxels` worktree, after the deterministic
and seeded DDA/reference parity tests passed.

## Environment and parameters

- Julia: 1.12.1
- Platform: Darwin aarch64
- CPU: Apple M3
- Julia threads available: 10
- Solver execution: deterministic single-thread CPU
- Boundary: `:periodic`
- Production traversal: `:dda`
- Ray sampling: 4 requested rays per horizontal voxel
- Cache turtle: 16 sectors
- Repeated application: 48 steps of 1,800 s
- Benchmark statistic: minimum of repeated `BenchmarkTools` samples after one
  explicit warmup call (`20` samples for a ray, `5` for most directions, `3`
  for the largest and repeated cases)

The benchmark definitions and machine-readable parameter tuple live in
`benchmark/voxel_benchmarks.jl`.

## Measurements

| Workload | Time | Allocated | Allocations |
| --- | ---: | ---: | ---: |
| one ray, sorted-plane reference, 16³ | 583 ns | 3,856 B | 13 |
| one ray, DDA, 16³ | 166 ns | 1,840 B | 2 |
| one directional response, 8³ | 52.750 µs | 291,024 B | 517 |
| one directional response, 16³ | 352.042 µs | 1,917,136 B | 2,053 |
| one directional response, 32³ | 3.018291 ms | 13,107,408 B | 12,293 |
| construct 16-direction cache, 16³ | 17.047583 ms | 66,415,424 B | 45,316 |
| apply existing cache for 48 steps, 16³ | 4.636959 ms | 6,733,824 B | 12,192 |

The scalar DDA is about 3.5 times faster than the sorted-plane reference for
the recorded ray and reduces allocations from 13 to 2. Once a 16-direction
cache exists, 48 spectral applications take about 0.097 ms per step in this
configuration; this is intentionally measured separately from the 17.0 ms
cache construction.

## Allocation interpretation and next optimization

`trace_voxel_ray` returns a materialized `VoxelRayPath`, so its two DDA
allocations are part of the public diagnostic result. Directional interception
currently materializes one such path for every launched ray; the roughly two
allocations per ray dominate the directional allocation counts.

A future optimization should add an internal DDA callback or accumulator that
applies Beer–Lambert attenuation as each segment is produced. The public
`VoxelRayPath` API should remain available for inspection and tests, while the
directional solver can avoid storing short-lived path vectors. Such a change
must preserve deterministic accumulation order and pass the same reference,
Java-parity, and energy-balance tests before becoming the default.
