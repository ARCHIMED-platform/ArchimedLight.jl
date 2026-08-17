# Temporary design note: terrain and soil in the voxel light model

Status: future implementation note for the `voxels` branch.

This file records the intended soil/topography design so it is not lost while
voxel interception and scattering are developed. It is not an implementation
claim. Remove it after the design has been implemented, tested, and transferred
to durable user/developer documentation.

## Problem to avoid

The canopy can remain inside a rectangular Cartesian voxel grid even when the
soil is not planar. However, voxels below the terrain surface must **not** be
treated as ordinary zero-PAD voxels.

In the current canopy interpretation:

```text
PAD > 0  -> foliage-containing air volume
PAD = 0  -> transparent air
```

Therefore, representing underground space with `PAD = 0` would allow rays to
travel below the soil surface as if the ground did not exist.

## Chosen architecture

Use a hybrid representation:

```text
rectangular voxel grid
  foliage -> volumetric PAD and Beer-Lambert attenuation
  air     -> PAD = 0

terrain surface
  planar surface, height field, or triangle mesh
  soil optical properties
  local surface normal
```

The soil is a surface receptor and scattering source. It is not volumetric PAD.
Everything below the first terrain surface must be excluded from voxel transport.

This preserves the regular voxel data structure while supporting slopes and
irregular topography. It is also compatible with a RATP-like exchange model:
foliage voxels and soil patches are both radiative receptors/sources, but foliage
uses volumetric extinction while soil uses surface reflection and absorption.

## Practical AMAPVox/LiDAR data path

Terrain must not be inferred from PAD. A LiDAR scene should produce two aligned
but distinct spatial products:

```text
LiDAR point cloud
  |-- vegetation returns --> AMAPVox voxelization --> PAD voxel grid (.vox)
  `-- ground returns ------> interpolation --------> bare-earth DTM (.asc)
                                                        |
                                                        `--> terrain height field
```

The practical workflow is:

1. Classify ground returns in the LiDAR point cloud and interpolate them into a
   bare-earth digital terrain model (DTM). This preprocessing is separate from
   PAD estimation; PAD alone cannot identify the terrain surface.
2. Provide the DTM to AMAPVox as an Arc/Info ASCII Grid (AAIGrid). AMAPVox uses
   it to compute voxel distance to ground, filter ground returns, and optionally
   estimate energy reaching the ground.
3. Preserve the same DTM alongside the AMAPVox `.vox` output and configuration.
4. Load the `.vox` PAD as the rectangular canopy grid and load the DTM as an
   independent `HeightFieldTerrain` in ArchimedLight.
5. Apply or validate the same coordinate transformation for both products. The
   terrain must cover the full horizontal simulation domain and have no unresolved
   `NODATA` cells inside it.
6. Provide PAR and NIR soil reflectance separately. Neither PAD nor terrain
   elevation supplies soil optical properties.

The original DTM should be the authoritative terrain geometry. It may be finer
than the canopy voxel grid and does not need to share its horizontal resolution.
The importer must handle raster cell-centre/corner conventions explicitly to
avoid a half-cell horizontal offset.

If the AMAPVox voxel file contains `ground_distance` computed from that DTM, it
provides a useful alignment check. For any voxel centre in horizontal column
`(i, j)`:

```text
z_ground(i, j) = z_voxel_center(i, j, k) - ground_distance(i, j, k)
```

Values reconstructed from different `k` levels in a column should agree within
floating-point tolerance. They can be used as a lower-resolution fallback if the
original DTM has been lost, subject to explicit gap handling. They must not be
interpreted as real topography when AMAPVox was run without a DTM; in that case
`ground_distance` represents height relative to its default horizontal reference
plane.

Relevant AMAPVox documentation:

- [voxelization and DTM input](https://amapvox.org/articles/Voxelization.html);
- [ground layer, ground elevation, and ground energy](https://amapvox.org/reference/ground.html).

## Ray interaction contract

For every propagating ray, the transport implementation must select the nearest
positive event among:

1. the next voxel boundary/segment endpoint;
2. an intersection with the terrain surface;
3. escape from the simulation domain.

A terrain hit terminates the current volumetric path. No later voxel segment
below that hit may receive energy.

The terrain intersection result should contain at least:

```julia
TerrainHit(
    distance,
    point,
    normal,
    patch_id,
    material_id,
)
```

The exact names and representation can change, but distance, local normal, and
material identity are required by the transport and scattering algorithms.

## Soil energy partition

For incident energy `E` and soil reflectance `rho_soil`:

```text
absorbed_soil  = (1 - rho_soil) * E
scattered_soil = rho_soil * E
```

Soil transmittance is zero in the initial model. Validate:

```text
0 <= rho_soil <= 1
```

PAR and NIR reflectances must be stored separately. If spectral soil data are
available, reduce them to each model band with the same documented spectral
weighting convention used for leaf optics.

## Lambertian reflection on non-planar terrain

For a Lambertian soil patch with unit normal `n`, outgoing directions over the
local upper hemisphere follow:

\[
p(\omega) = \frac{\max(0, n \cdot \omega)}{\pi}.
\]

The local terrain normal is essential. Reusing the global vertical normal on a
slope would distort both the angular distribution of reflected radiation and its
coupling to nearby foliage.

For deterministic discrete directions, compute patch-specific weights from the
global full-sphere quadrature:

```text
raw_weight[d] = solid_angle_weight[d] * max(0, dot(normal, direction[d]))
weight[d]     = raw_weight[d] / sum(raw_weight)
```

Only directions in the hemisphere above the local surface may receive energy.
The discrete weights must sum to one within numerical tolerance.

The terrain normal is **not** an additional transport direction. ArchimedLight
continues to trace only the fixed global discrete directions used by the voxel
solver. The normal only determines which of those directions lie above the local
surface and how the reflected energy is weighted among them.

For an incoming ray in discrete direction `d_in` that hits the terrain with
energy `E`, the initial Lambertian soil model is:

```text
reflected = rho_soil * E

for each existing discrete direction d_out:
    raw[d_out] = solid_angle_weight[d_out]
                 * max(0, dot(normal, d_out))

weight[d_out] = raw[d_out] / sum(raw)
E_out[d_out]  = reflected * weight[d_out]
```

No ray is launched along the exact terrain normal unless that direction already
belongs to the quadrature. A slope therefore rotates the energy distribution by
changing the weights over existing direction indices, not by rotating or adding
directions. Specular reflection, which would generally produce an outgoing
direction absent from the quadrature, is outside the initial Lambertian scope.

## Why not binary soil voxels?

An alternative representation is:

```julia
material[i, j, k] = Air | Foliage | Soil
```

with all cells below the terrain classified as `Soil`. This prevents underground
ray propagation, but a binary occupancy surface is stair-stepped and usually has
axis-aligned normals. Unless voxels are extremely small, this can noticeably bias:

- terrain intersection distances;
- shadow boundaries;
- reflected-light directions;
- energy received by foliage on slopes.

Binary soil voxels may remain useful as an acceleration mask, but they should not
define the final geometric surface normal. If this route is ever used, reconstruct
a sub-voxel plane and normal at the air/soil interface.

## Proposed terrain interface

Design the first planar implementation behind a generic interface so non-planar
terrain does not require changing the scattering solver later.

Possible conceptual API:

```julia
abstract type AbstractVoxelTerrain end

intersect_terrain(terrain, origin, direction, max_distance)
soil_optics(terrain, hit, band)
terrain_bounds(terrain)
```

Planned concrete representations:

1. `PlanarTerrain`: horizontal plane, used for the first implementation and
   analytical tests.
2. `HeightFieldTerrain`: regular `z_soil(x, y)` samples matching or independent
   of the canopy grid.
3. `TriangulatedTerrain`: exact triangle intersections, including a triangulated
   height field and, later, more general terrain meshes.

The transport algorithm should depend only on the generic intersection result,
not on which terrain representation produced it.

## Implementation tasks

### Phase 1: scientific and API contract

- [ ] Define terrain coordinates and their relation to `VoxelGrid.minimum` and
  `VoxelGrid.maximum`.
- [ ] Define whether the terrain must cover the full horizontal domain.
- [ ] Define behavior when terrain lies above the voxel-grid minimum or intersects
  multiple vertical voxel layers.
- [ ] Define open and horizontally periodic terrain semantics.
- [ ] Define the normal orientation convention; returned normals must point from
  soil toward air.
- [ ] Define PAR/NIR soil-optics validation and band averaging.
- [ ] Keep soil reflectance separate from leaf scattering coefficients.

### Phase 1b: AMAPVox terrain input contract

- [ ] Define a scene-level input contract containing the AMAPVox `.vox`, the
  authoritative DTM, and the AMAPVox coordinate transformation/configuration.
- [ ] Add an AAIGrid height-field reader or a documented preprocessing route to
  a supported height-field representation.
- [ ] Handle raster origin, cell-centre/corner conventions, resolution, units,
  coordinate reference system, transformations, and `NODATA` explicitly.
- [ ] Permit DTM and canopy-grid resolutions to differ without silently shifting
  or cropping the terrain.
- [ ] Require terrain coverage over the horizontal simulation domain, or define
  missing-terrain boundary behavior explicitly.
- [ ] When `ground_distance` is present, validate sampled DTM elevations against
  `z_voxel_center - ground_distance`.
- [ ] Define the lower-resolution reconstruction fallback from `ground_distance`
  and reject it when a real DTM was not used by AMAPVox.

### Phase 2: generic terrain types

- [ ] Add the abstract terrain interface and a `TerrainHit` representation.
- [ ] Add `SoilOpticalProperties` with explicit PAR and NIR reflectance.
- [ ] Implement `PlanarTerrain` without special-casing it inside the scattering
  solver.
- [ ] Add constructors that validate finite geometry, normalized normals, material
  identifiers, and optical bounds.
- [ ] Define a no-terrain value representing an open bottom boundary.

### Phase 3: terrain-aware voxel traversal

- [ ] Compute the nearest terrain hit for each internal or external ray.
- [ ] Truncate the DDA path at the terrain-hit distance, including hits occurring
  inside a voxel segment.
- [ ] Ensure no energy is deposited below the terrain surface.
- [ ] Return distinct termination reasons such as `:terrain`, `:top`, `:bottom`,
  and `:side`.
- [ ] Preserve the existing no-terrain traversal behavior exactly.
- [ ] Decide whether a coarse subterranean mask is worthwhile as an acceleration
  structure, without using it for surface normals.

### Phase 4: planar soil absorption and reflection

- [ ] Partition terrain-hit energy into soil absorption and reflection.
- [ ] Implement local-normal Lambertian redistribution using the full-sphere voxel
  scattering quadrature.
- [ ] Ensure terrain scattering emits energy only through existing discrete
  direction indices; the local normal must never create a continuous direction.
- [ ] Cache direction weights per terrain patch or exact normal when beneficial;
  any normal-binning approximation must be explicit and tested.
- [ ] Feed reflected soil energy back into upward voxel transport.
- [ ] Track soil absorption and reflection separately for PAR and NIR.
- [ ] Include soil terms in the complete energy-balance report.
- [ ] Prevent immediate numerical self-intersection when launching a reflected ray
  from the soil surface.

### Phase 5: height fields and triangulated terrain

- [ ] Add a height-field representation `z_soil(x, y)`.
- [ ] Triangulate each height-field cell with a documented diagonal convention.
- [ ] Compute exact first ray-triangle intersections.
- [ ] Return either facet normals or interpolated vertex normals according to an
  explicit setting; use facet normals for the initial scientific reference.
- [ ] Add an acceleration structure so terrain intersection does not test every
  triangle for every ray.
- [ ] Verify continuity and intersection behavior along triangle edges and vertices.
- [ ] Define periodic tiling of height fields, or reject incompatible periodic
  configurations explicitly.

### Phase 6: caching and repeated timesteps

- [ ] Cache terrain geometry and its acceleration structure independently of
  meteorological timesteps.
- [ ] Invalidate terrain caches when heights, vertices, topology, domain tiling, or
  soil material assignment changes.
- [ ] Avoid invalidating geometry caches when only PAR/NIR soil reflectance changes.
- [ ] Benchmark planar, height-field, and mesh intersection costs separately from
  multiple-scattering iteration.

## Required tests

- [ ] No-terrain mode reproduces the current voxel traversal results.
- [ ] A horizontal `HeightFieldTerrain` reproduces `PlanarTerrain` intersections.
- [ ] Rays stop at soil even when lower voxels have `PAD = 0`.
- [ ] No subterranean voxel receives intercepted or scattered energy.
- [ ] A ray hitting a slope returns the correct intersection point and outward normal.
- [ ] A black soil (`rho = 0`) absorbs all terrain-hit energy.
- [ ] A perfectly reflecting test soil (`rho = 1`) absorbs none and conserves energy.
- [ ] Lambertian discrete-direction weights are non-negative and sum to one.
- [ ] A horizontal surface produces a symmetric upward distribution.
- [ ] A sloped surface rotates the reflected distribution with the local normal.
- [ ] Sloped terrain changes weights over the fixed direction set without adding
  or rotating any transport direction.
- [ ] For several arbitrary terrain normals, all outgoing discrete weights are
  non-negative, use only directions above the local surface, and sum to one.
- [ ] A synthetic DTM and matching `ground_distance` field reconstruct identical
  terrain elevations within tolerance for every tested horizontal column.
- [ ] DTM coordinate offsets, incomplete horizontal coverage, and internal
  `NODATA` cells fail validation rather than silently creating open ground.
- [ ] Triangle-edge and triangle-vertex hits are deterministic and do not duplicate
  energy.
- [ ] Open lateral boundaries distinguish side escape from terrain interception.
- [ ] Periodic terrain either wraps consistently with the canopy or fails with a
  clear validation error.
- [ ] Complete PAR and NIR energy balances include foliage absorption, soil
  absorption, unresolved scattering residual, and all escape terms.

## Energy-balance invariant

For each spectral band:

```text
external incoming + explicit internal injection
    == foliage absorption
     + soil absorption
     + unresolved scattering residual
     + top escape
     + side escape
     + bottom escape when no terrain is present
```

Soil-reflected energy is an internal transfer term. It must not be counted as a
new external source.

## Recommended delivery order

1. Generic terrain interface plus horizontal `PlanarTerrain`.
2. Terrain-aware DDA truncation and black-soil absorption.
3. Lambertian planar-soil scattering.
4. Height-field storage and triangulation.
5. Exact mesh intersections and local-normal reflection.
6. Acceleration/caching, benchmarks, and user documentation.

The planar implementation is intentionally first, but its API must already use
the generic `terrain intersection + local normal + optical properties` contract.

## Definition of done

- [ ] Underground space is never treated as transparent air.
- [ ] The rectangular canopy voxel representation remains unchanged.
- [ ] Planar and non-planar terrain use the same transport/scattering interface.
- [ ] Soil absorption and reflection are spectrally explicit and energy-conserving.
- [ ] Sloped surfaces use their local normal for Lambertian scattering.
- [ ] Terrain normals redistribute energy over the existing discrete quadrature;
  they never introduce untracked propagation directions.
- [ ] The AMAPVox `.vox` and authoritative DTM are aligned and their import
  contract, transformations, resolution semantics, and validation are documented.
- [ ] Height-field/mesh intersections prevent all below-surface voxel transport.
- [ ] Boundary and cache semantics are documented and tested.
- [ ] Durable documentation replaces this temporary note.
