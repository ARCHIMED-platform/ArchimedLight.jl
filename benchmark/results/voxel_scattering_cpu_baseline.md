# Voxel scattering CPU baseline

Recorded on 2026-08-17 from `benchmark/voxel_benchmarks.jl`.

- Julia: 1.12.1
- CPU: Apple M3 (`aarch64`), 10 Julia threads available
- OS: Darwin
- backend: periodic DDA, `G = 0.5`, four first-order rays per horizontal voxel
- optical example: total scattering `0.17` PAR and `0.85` NIR
- stopping ratio: `1e-6`; maximum iterations: 20
- statistic: minimum of five BenchmarkTools samples after compilation

These are development baselines, not cross-machine performance guarantees.
`memory` is cumulative Julia allocation during the timed call, not peak RSS.

| Operation | Grid | Occupancy | Full-sphere directions | Time (ms) | Allocation (MiB) | Allocations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Build transport paths | 4×4×4 | dense | 12 | 0.172 | 0.869 | 1,763 |
| Build transport paths | 8×8×8 | dense | 12 | 1.462 | 7.872 | 13,102 |
| Build transport paths | 8×8×8 | about 20% | 12 | 0.371 | 2.147 | 3,256 |
| Build transport paths | 4×4×4 | dense | 32 | 0.461 | 2.306 | 4,680 |
| Apply one cached transport order | 4×4×4 | dense | 12 | 0.046 | 0.0008 | 4 |
| Eight complete solves, cached paths | 4×4×4 | dense | 12 | 8.386 | 0.211 | 1,240 |
| Eight complete solves, rebuild paths | 4×4×4 | dense | 12 | 9.824 | 7.160 | 15,344 |

First-order ray sampling was profiled separately on the dense 8×8×8 grid:

| Requested rays per horizontal voxel | Direction response time (ms) | Allocation (MiB) | Allocations |
| ---: | ---: | ---: | ---: |
| 1 | 0.0166 | 0.073 | 135 |
| 4 | 0.0738 | 0.278 | 519 |
| 16 | 0.2681 | 1.099 | 2,055 |

This sampling parameter scales first-order illumination nearly linearly in the
recorded range. Internal scattering uses one deterministic centre source per
occupied voxel and direction, so `rays_per_voxel` does not multiply its cached
path count; it remains in the backend fingerprint because it changes the
first-order source energy.

Reusing paths reduced this warmed eight-solve workload by about 15% in wall
time and by about 34× in cumulative allocation. The advantage grows with grid
size and angular resolution because setup scales with occupied voxels times
directions. Sparse setup is cheaper because paths are stored only for occupied
foliage sources; grid-sized work arrays remain dense.

One cached transport application is the cost of one scattering order, while
the complete-solve row includes convergence. For the recorded complete solve,
PAR converged in 6 iterations with a relative
energy residual of `2.69e-16`; NIR converged in 16 iterations with a residual
of `1.35e-16`. The higher NIR iteration count follows directly from its larger
single-scattering albedo.

The production representation stores path segments plus grid-sized work arrays.
It never constructs a voxel-by-voxel dense exchange matrix. The dense matrix
helper used by tests is intentionally excluded from these benchmarks.
