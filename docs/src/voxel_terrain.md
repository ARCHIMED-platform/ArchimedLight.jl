# Voxel Terrain And Soil

The voxel solver represents vegetation and terrain with two distinct geometric
objects:

```text
regular Cartesian VoxelGrid
  PAD > 0  participating foliage volume
  PAD = 0  transparent air

AbstractVoxelTerrain
  continuous soil boundary
  patch identity, local normal, and PAR/NIR soil optics
```

A zero-PAD voxel is never used as a soil marker. A ray stops at the first exact
terrain intersection, even when that intersection lies inside a voxel segment.
The remaining part of that segment and every voxel below it receive no energy.
No coarse subterranean voxel mask is maintained: exact nearest-hit truncation
plus the terrain BVH already removes below-surface transport, while a second
occupancy representation would duplicate geometry and risk inconsistent
normals.

```@contents
Pages = ["voxel_terrain.md"]
Depth = 3
```

## Explore The Geometry In 3D

The interactive view below shows the geometric contract used by the solver.
Its deterministic voxel field is arranged as a stylized single tree: a narrow
lower column makes the trunk visible and an irregular ellipsoidal crown carries
the foliage PAD. Empty cells remain transparent air. The terrain is an
independent triangulated surface that may cut through the lower voxel layers.
Incoming rays stop at their first terrain intersection, and soil-reflected
energy leaves along existing discrete transport directions weighted by the
local surface normal.

```@example voxel_terrain_3d
using ArchimedLight, Bonito, WGLMakie
Bonito.Page() # hide
include(joinpath(pkgdir(ArchimedLight), "docs", "voxel_terrain_3d.jl")) # hide
voxel_terrain_app()
```

!!! note
    The tree and terrain are deterministic pedagogical inputs. The brown
    stem-like column is a visualization scaffold, not a separate woody-stem
    optical model: ArchimedLight still interprets every non-zero cell as
    participating PAD. Every ray segment, soil hit, local normal, and reflected
    direction shown above is computed through `voxelplot`, `trace_voxel_ray`,
    and `terrain_lambertian_weights`; only the displayed ray energies are
    illustrative.

## Construct A Terrain

Soil is opaque in the initial model. Define its band reflectance explicitly;
absorptance is ``1-\rho`` and transmittance is zero.

```@example voxel_terrain
using ArchimedLight

grid = VoxelGrid((0, 0, 0), (2, 2, 3), fill(0.4, 2, 2, 3))
soil = SoilOpticalProperties(
    par_reflectance=0.15,
    nir_reflectance=0.40,
)

terrain = PlanarTerrain(grid, soil; elevation=0.6)
path = trace_voxel_ray(
    grid,
    (0.5, 0.5, 3.0),
    (0.1, 0.0, -1.0);
    boundary=:open,
    terrain=terrain,
)

(
    termination=path.exit_reason,
    hit_point=path.terrain_hit.point,
    hit_normal=path.terrain_hit.normal,
    crossed_voxels=[segment.index for segment in path],
)
```

`PlanarTerrain(grid, soil)` divides the plane into the same horizontal patch
layout as the voxel grid. A more general finite plane can choose its horizontal
bounds and patch count directly.

```julia
PlanarTerrain(
    elevation,
    (xmin, xmax),
    (ymin, ymax),
    soil;
    cells=(nx, ny),
    periodic=false,
)
```

### Height Fields

`HeightFieldTerrain` takes vertex coordinates, not cell centres. Its elevation
matrix has size `(length(x), length(y))`. Every cell uses the documented
southwest-to-northeast diagonal:

```julia
terrain = HeightFieldTerrain(
    [0.0, 1.0, 2.0],
    [0.0, 1.0, 2.0],
    [0.0 0.0 0.0;
     0.2 0.3 0.4;
     0.4 0.6 0.8],
    soil;
    normal_mode=:facet,
)
```

Ray intersections are exact triangle intersections. `normal_mode=:facet` is
the scientific reference and returns one normal per triangle.
`normal_mode=:interpolated` explicitly enables area-weighted vertex-normal
interpolation.

### Triangle Meshes

Use `TriangulatedTerrain` when non-overhanging topographic geometry is already
a mesh:

```julia
vertices = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.1),
    (1.0, 1.0, 0.2),
    (0.0, 1.0, 0.0),
]
triangles = [(1, 2, 3), (1, 3, 4)]
terrain = TriangulatedTerrain(vertices, triangles, soil)
```

The constructor validates finite, non-degenerate geometry with a positive
upward normal component and builds a
bounding-volume hierarchy. Intersections along shared edges or vertices choose
the smallest patch identifier at equal distance, making the result
deterministic without duplicating energy. Returned normals always point from
soil toward air (positive global z component).

Generic triangle meshes are non-periodic. Use a periodic height field when the
same terrain tile must repeat horizontally.

## Boundaries And Domain Contract

`NoVoxelTerrain()` is the default and preserves the historical open lower
boundary: downward energy may terminate as `:bottom`. With a concrete terrain,
termination is one of:

- `:terrain`: the nearest soil intersection;
- `:top`: vertical escape above the grid;
- `:side`: lateral escape with `boundary=:open`;
- `:bottom`: lower escape only with `NoVoxelTerrain()`.

For open horizontal boundaries, terrain must cover the complete voxel x/y
domain. It may lie above `grid.minimum[3]` and cross several vertical voxel
layers, but it must remain below the top. Terrain below the grid minimum is
rejected because it could not stop a ray before bottom escape.

For `boundary=:periodic`, both conditions are mandatory:

1. construct the planar or height-field terrain with `periodic=true`;
2. make its x/y bounds exactly match the voxel x/y bounds.

A periodic canopy with non-periodic terrain fails explicitly.

## Soil Reflection On Slopes

For an incident energy ``E`` and band reflectance ``\rho_b``:

```math
E_{soil,absorbed} = (1-\rho_b)E,
\qquad
E_{soil,reflected} = \rho_b E.
```

The reflected pool is Lambertian around the patch's local unit normal ``n``.
ArchimedLight does not add or rotate transport directions. It reweights the
existing full-sphere voxel quadrature:

```math
w_d =
\frac{\Delta\Omega_d\max(0,n\cdot\omega_d)}
     {\sum_j\Delta\Omega_j\max(0,n\cdot\omega_j)}.
```

`terrain_lambertian_weights(quadrature, normal)` returns the selected existing
direction indices and normalized non-negative weights. A strongly sloped patch
can legitimately select a direction with negative global z when that direction
still lies above the local surface. Reflected rays start a floating-point-safe
distance above the patch to avoid immediate self-intersection.

Soil-to-soil reflection on concave slopes is followed deterministically for up
to 128 surface bounces. Any remaining numerical tail is retained as unresolved
energy rather than silently converted to absorption.

## Run The Voxel Pipeline

Pass terrain independently of foliage optics:

```julia
result = run_voxel_light_step(
    grid,
    meteo_row,
    options;
    backend=VoxelCPUBackend(boundary=:open),
    optics=leaf_optics,
    terrain=terrain,
)
```

`ground=VoxelGroundOptics(...)` remains available for compatibility with the
old horizontal lower-boundary model. `ground` and `terrain` are mutually
exclusive. New topographic workflows should use `terrain`.

First-order results contain `par_terrain_energy` and `nir_terrain_energy`, one
value per terrain patch, plus their scalar `*_terrain_incident_energy` totals.
Multiple-scattering results expose soil diagnostics through
`soil_absorbed_energy(result)` and `soil_reflected_energy(result)`. The stored
`ground_absorbed_energy` and `ground_reflected_energy` field names remain as
compatibility aliases in the result containers.

For each spectral band, the checked terminal invariant is:

```text
external incoming + explicit compatibility injection
    == foliage absorption
     + soil absorption
     + unresolved scattering
     + top escape
     + side escape
     + bottom escape (NoVoxelTerrain only)
```

Soil-reflected energy is an internal transfer and is not added again on the
right-hand side.

## AMAPVox And An Authoritative DTM

PAD does not identify the ground. Keep the bare-earth DTM used by AMAPVox next
to the `.vox` file and its configuration:

```text
classified LiDAR
  vegetation returns -> AMAPVox -> PAD voxel table
  ground returns      -> DTM     -> AAIGrid terrain
```

AMAPVox accepts an Arc/Info ASCII Grid DTM and uses it to compute distance to
ground; its voxelization configuration may also apply a 4×4 transformation
matrix. See the official
[voxelization and DTM documentation](https://amapvox.org/articles/Voxelization.html)
and [ground-layer documentation](https://amapvox.org/reference/ground.html).

### AAIGrid Semantics

`read_aai_grid` handles both lower-left corner and lower-left cell-centre
headers. Centre headers are converted to corners by subtracting half a cell.
AAIGrid data rows are north-to-south; the reader reverses them so Julia indices
increase with y. `NODATA` and non-finite elevations are rejected anywhere in
the domain.

AAIGrid elevations are cell-centre samples. `height_field_terrain(aai, soil)`
keeps every sample as an exact terrain vertex and adds the raster boundary with
constant nearest-cell extension over the outer half-cell. It then applies the
fixed triangulation convention between those vertices. This preserves every
authoritative cell-centre elevation while making horizontal coverage exactly
equal to the raster extent. DTM and canopy voxel resolution may differ.

### Coordinate Transformation And CRS

`TerrainCoordinateTransform` stores the explicit affine matrix mapping DTM
coordinates into the coordinate frame already used by the `.vox` grid:

```julia
transform = TerrainCoordinateTransform(
    matrix_4x4;
    source_crs="EPSG:2154",
    target_crs="EPSG:2154",
    source_units="m",
    target_units="m",
)
```

The unit labels record the contract while the matrix performs any required
numeric conversion. Axis-aligned scale/translation retains a `HeightFieldTerrain`. Rotation or
shear returns a `TriangulatedTerrain`; periodic tiling is then rejected. Unit
conversion belongs in this matrix and is therefore auditable rather than
guessed from coordinate magnitudes.

### Scene-Level Import

The input object keeps the `.vox`, authoritative DTM, optional AMAPVox
configuration path, transform, and DTM provenance together:

```julia
input = VoxelTerrainSceneInput(
    "canopy.vox",
    "bare_earth.asc";
    configuration_path="voxelization.xml",
    dtm_to_voxel=transform,
    dtm_used_by_amapvox=true,
)

scene = read_voxel_terrain_scene(input, soil; periodic=false)
grid = scene.grid
terrain = scene.terrain
```

The importer rejects horizontal offsets, incomplete coverage, invalid vertical
extent, and internal `NODATA`. If the `.vox` table contains `ground_distance`,
it also checks every voxel centre using:

```math
z_{ground}(i,j) = z_{voxel\ centre}(i,j,k) - ground\_distance(i,j,k).
```

`scene.alignment` reports sample count, maximum absolute error, and mean
absolute error. Use `alignment_atol` and `alignment_rtol` only to represent the
known precision of exported data, not to hide a coordinate offset.

### Ground-Distance Fallback

If the original DTM has been lost, `terrain_from_ground_distance` reconstructs
one lower-resolution height sample per voxel column. Values reconstructed from
all k levels must agree. The required keyword
`dtm_used_by_amapvox=true` prevents interpreting AMAPVox's default horizontal
reference plane as measured terrain.

## Cache Semantics

First-order response and scattering-path caches include a terrain fingerprint.
It contains elevations/vertices, topology, patch material assignment, normal
mode, bounds, and periodic tiling. Reusing a cache after any of those changes
fails validation.

PAR/NIR soil reflectance is deliberately excluded from the geometry
fingerprint. Changing only a material's optical coefficients reuses the same
ray paths. Local-normal direction indices and normalized Lambertian weights are
cached once per patch with those paths. Terrain objects and their internal acceleration structures should be
treated as read-only after construction; construct a new terrain after editing
height or mesh arrays.

## Public Terrain API

```text
AbstractVoxelTerrain, NoVoxelTerrain
PlanarTerrain, HeightFieldTerrain, TriangulatedTerrain
TerrainHit, SoilOpticalProperties
intersect_terrain, terrain_bounds
terrain_patch_count, terrain_patch_point
terrain_patch_normal, terrain_patch_material
soil_optics, terrain_lambertian_weights
AAIGrid, read_aai_grid, height_field_terrain
TerrainCoordinateTransform
VoxelTerrainSceneInput, VoxelTerrainScene
read_voxel_column, read_voxel_terrain_scene
validate_ground_distance, terrain_from_ground_distance
```
