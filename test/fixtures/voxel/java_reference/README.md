# Historical Java voxel reference

These files freeze the local Java implementation used to port voxel light
interception to ArchimedLight.jl. The Java source itself lives under the
gitignored `java_implementation/` directory and is not part of this fixture.

## Provenance

- Baseline captured: 2026-08-17
- Runtime: OpenJDK 1.8.0_442
- JUnit: 4.12
- Vecmath: 1.5.2
- Commons Math: 3.6.1
- Source hashes: `source_sha256.txt`

The original `TestVoxel1` through `TestVoxel4` classes were compiled unchanged.
Two compile-time stubs were used only for the `VoxelParameters` and
`WordIterator` types referenced by unused `.vox` I/O methods in `VoxelSpace`.
The original test run completed with:

```text
JUnit version 4.12
.......
Time: 1,286

OK (7 tests)
```

`VoxelBaseline.java` generated `baseline.tsv`. It records three ray paths,
three complete directional response fields, and compact layer profiles for a
ground layer and a central pillar on small deterministic grids. The output is
kept at the precision printed by Java.

The `:java_reference` compatibility tests compare these printed path lengths
and Float32-accumulated coefficients exactly; they do not loosen a tolerance to
hide drift. For the corrected Julia DDA, a separate seeded comparison covered
3,720 segments from 256 periodic/open rays. Its largest DDA/reference
difference was `7.105427357601002e-15` in absolute value and
`1.0442467563577877e-13` in relative value. The acceptance limits are therefore
`atol=2e-14` and `rtol=2e-13`, leaving a small floating-point margin while
remaining close to the observed error.

## Known historical behavior

`Ray.distToVoxelWallsV2()` does not always select the nearest positive x/y
wall when an origin has a positive coordinate and the ray travels in the
positive direction. Consequently, the `diagonal` and `corner` baselines can
contain consecutive segments associated with the same voxel. These rows are a
compatibility record, not the desired geometric contract for the Julia DDA
kernel. Julia tests must distinguish deliberate compatibility from the
production invariant that adjacent path segments are separated by an actual
voxel boundary.

The Java non-toric mode also wraps horizontal indices and resets ray intensity
at a detected wall instead of terminating the ray. It is preserved in the
`diagonal_java_nontoric` response for characterization only. The Julia public
`:open` mode terminates a ray at the first horizontal exit.
