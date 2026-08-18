# Future Directions For Voxel Radiative Transfer

This document records possible extensions of the voxel workflow. It is a
technical roadmap, not a commitment that every item will be implemented. Its
purpose is to keep the scientific assumptions, architectural dependencies,
and validation requirements visible beyond the current CPU interception and
isotropic multiple-scattering solver.

```@contents
Pages = ["voxel_future_directions.md"]
Depth = 3
```

## Current Baseline

The current [voxel workflow](voxel_interception.md) represents vegetation as plant area density
(`PAD`) in a regular Cartesian grid. It computes first-order PAR and NIR
interception with Beer–Lambert attenuation, then optional deterministic
isotropic multiple scattering with explicit scalar or per-voxel optics. It
supports a spherical leaf-angle distribution (`G = 0.5`), open or periodic
horizontal boundaries, continuous planar/height-field/triangle terrain with
local-normal Lambertian soil reflection, and cached CPU DDA/BVH paths.
Directional responses and internal transport paths can be reused over a
meteorological series.

The current implementation deliberately does not infer unsupported physics.
It does not yet include non-spherical leaf-angle distributions, anisotropic
leaf scattering, separate vegetation material identities, adaptive grids, or
GPU execution.

Future extensions should preserve four properties of this baseline:

1. a deterministic CPU implementation remains the numerical reference;
2. energy balances are explicit and close within documented tolerances;
3. geometry-dependent work can be reused across meteorological timesteps;
4. unsupported combinations fail explicitly rather than returning plausible
   but incomplete results.

## Proposed Order Of Work

The directions below are related and should not be implemented independently
without considering their shared data requirements.

| Phase | Direction | Main prerequisite | Intended outcome |
| --- | --- | --- | --- |
| 0 | Harden the CPU transport kernel | Current DDA/reference parity | Allocation-light reference kernel suitable for reuse and porting |
| 1 | Add GPU execution | Stable matrix-free CPU kernels | Accelerated response construction and scattering without changing results |
| 2 | Add structural and angular detail | Validated core solver | Multiple components, non-spherical LIDF, and anisotropic scattering |
| 3 | Add alternative spatial representations | Stable material and backend interfaces | Point-cloud ingestion, non-uniform grids, or octrees where justified |

The GPU prototype can begin before phase 3 is complete, but the public GPU
backend should follow the same equations and conservation rules as the CPU
reference.

## Multiple Scattering

### Optical Meaning

For a wavelength band `b`, the single-scattering coefficient of a vegetation
material is

```math
\omega_b = \rho_b + \tau_b,
```

where ``\rho_b`` is leaf reflectance and ``\tau_b`` is leaf transmittance. The
absorbed fraction is

```math
\alpha_b = 1 - \omega_b.
```

The coefficient must remain in ``[0, 1]``. Under an isotropic scattering
assumption, the sum ``\rho_b + \tau_b`` is sufficient to compute the total
scattered pool. Separate reflectance and transmittance become necessary when
the model distinguishes the two sides of leaves or uses an anisotropic phase
function.

Optical properties should be supplied by a material or vegetation-component
definition. Global fallback coefficients can remain convenient, but they
should not silently replace known species- or organ-specific measurements.

### Deterministic Transport Operator

The implemented CPU reference follows the exchange-coefficient approach used
by RATP. Directional ray traversal defines how energy emitted from one
voxel is intercepted by other voxels or escapes from the domain. This defines
a transport operator ``T`` without requiring a dense voxel-by-voxel matrix.

For scattering order ``n``:

```math
E^{(n+1)} = T\left(\omega \odot E^{(n)}\right),
```

where ``E^{(0)}`` is the first-order intercepted energy. Total absorbed energy
is accumulated as

```math
A = \sum_{n \ge 0}(1 - \omega) \odot E^{(n)}.
```

Iteration stops when the remaining scattered pool is smaller than a documented
fraction of the first-order energy, or when a maximum number of iterations is
reached. Every iteration must account for:

- energy absorbed by vegetation;
- energy intercepted by another voxel;
- energy intercepted by soil;
- energy leaving the domain.

This formulation is equivalent to applying successive orders of scattering.
It is deterministic, reuses the current directional discretization, and can be
implemented matrix-free so that a dense ``N_voxel × N_voxel`` matrix is never
materialized.

The validated baseline includes scalar and dense optical fields, a fixed
full-sphere quadrature derived from turtle directions, centre-of-voxel internal
sources, convergence diagnostics, boundary-specific escape, and a planar
Lambertian ground. Future work can distinguish reflection and transmission in
an angular phase function or use Monte Carlo photon tracing as an independent
oracle for more complex phase functions. Stochastic sampling noise remains
undesirable in production voxel-level absorbed-PAR outputs.

### Result Semantics

The scattering result distinguishes quantities that coincide in the
first-order-only model:

- incident energy at a voxel;
- first-order intercepted energy;
- energy received from scattering;
- absorbed energy;
- energy re-emitted into scattering;
- escaped energy;
- scattering iterations, unresolved energy, and convergence status.

The words *intercepted* and *absorbed* must not be used interchangeably once
``\omega > 0``.

## Soil And Non-Planar Terrain Baseline

### Representation

The implemented voxel domain remains a rectangular Cartesian box even when the soil
surface is not planar. Vegetation remains a volumetric medium inside this box,
while terrain is represented as a boundary surface.

A zero-`PAD` voxel means transparent air. It must therefore not be used to
represent underground material: a ray would incorrectly propagate through it.
The transport model needs an explicit distinction between:

- air;
- foliage or other participating media;
- soil or excluded below-ground space.

The implemented hybrid representation is:

- foliage is stored in the regular `PAD` grid;
- terrain is a height field or triangle mesh;
- paths are truncated at terrain before below-surface portions can receive energy;
- ray propagation stops or scatters at the exact terrain intersection.

### Interaction With Soil

For opaque soil with band reflectance ``\rho_{soil,b}``:

```math
A_{soil,b} = (1 - \rho_{soil,b})E_{soil,b},
```

```math
S_{soil,b} = \rho_{soil,b}E_{soil,b}.
```

The first implementation can use Lambertian reflection. Outgoing directions
are then weighted over the hemisphere around the local terrain normal ``n``:

```math
p(\Omega) = \frac{\max(0, n \cdot \Omega)}{\pi}.
```

Using local surface normals is important on slopes. Treating the boundary of
binary solid voxels as the soil surface would create staircase geometry and
axis-aligned normals, which can bias both shadows and reflected radiation.

### Implemented Scope And Further Extensions

The reference solver now supports finite planes, regular height fields,
arbitrary non-periodic triangle meshes, exact intersections, facet or explicit
interpolated normals, BVH acceleration, aligned AAIGrid/AMAPVox input, and
geometry-only cache invalidation. Open and periodic semantics are explicit;
periodic canopy transport over non-periodic topography is rejected.

Potential future terrain work includes wavelength-resolved or BRDF soil
reflection, robust overhang semantics for general meshes, external GIS raster
reprojection, and sub-patch integration of reflected source position. These
extensions must retain the energy and nearest-event contracts documented in
[Voxel Terrain And Soil](voxel_terrain.md).

## GPU Computing

### Intended Targets

The most useful GPU targets are the workloads that dominate large voxel
simulations:

- constructing directional responses for many rays and directions;
- applying matrix-free scattering transport over many occupied voxels;
- repeating spectral calculations for many meteorological timesteps;
- reducing per-ray interception into voxel energy arrays.

GPU execution should be introduced as a backend, not as a separate scientific
model. The CPU and GPU backends must share the same conventions for ray
directions, boundaries, attenuation, scattering coefficients, soil, and output
units.

### CPU Preparation Before Porting

The current diagnostic DDA API materializes a `VoxelRayPath`. That is useful
for inspection and tests but unnecessarily allocates during bulk directional
interception. Before a GPU implementation, the CPU solver should gain an
internal streaming traversal that applies attenuation as each segment is
generated.

The public path-producing API should remain available. The optimized internal
kernel should:

- avoid one short-lived path vector per ray;
- use flat, concrete arrays and GPU-compatible scalar data;
- separate traversal from accumulation clearly;
- avoid hidden global state;
- preserve deterministic CPU accumulation order where practical.

### GPU Data And Execution Strategy

Candidate parallel dimensions are rays, sky directions, wavelength bands, and
independent scattering sources. A first implementation should parallelize
rays within one direction and use an explicit reduction strategy for voxel
accumulation. Uncontrolled atomic accumulation may introduce both contention
and non-deterministic floating-point differences.

Important design requirements include:

- keep grid, optical properties, response data, and outputs resident on the
  device across timesteps;
- avoid host–device transfers inside scattering iterations;
- benchmark cache construction separately from repeated cache application;
- report both runtime and allocated/device memory;
- fall back explicitly when a requested feature is unsupported on the GPU.

The CPU implementation remains the oracle. GPU tests should compare energy
budgets and voxel fields with tolerance-based criteria appropriate to parallel
floating-point reduction, rather than requiring bitwise identity.

## Leaf-Angle Distributions And Directional Scattering

The current ``G = 0.5`` assumption corresponds to a spherical leaf-angle
distribution. Future support could provide:

- common parametric leaf inclination distributions;
- empirical distributions supplied per vegetation component;
- direction-dependent ``G(\Omega)`` values;
- spatially varying leaf-angle distributions where measurements justify them.

Extinction and scattering must remain conceptually separate. A directional
extinction coefficient determines the probability of interception; it does not
by itself determine the direction into which intercepted light is scattered.
An anisotropic implementation additionally needs a phase function derived from
leaf orientation and optical properties.

The first non-spherical implementation should be validated against analytical
homogeneous-canopy cases before allowing voxel-specific distributions.

## Multiple Vegetation Components

A single scalar `PAD` field cannot distinguish species, organs, living and dead
material, or vegetation and woody elements. A future grid may store plant area
density by component:

```text
PAD[component, i, j, k]
```

Each component can then define its own:

- PAR and NIR reflectance and transmittance;
- leaf-angle distribution;
- scattering phase approximation;
- output identity for coupling to physiological models.

Within a mixed voxel, total extinction follows the sum of component extinction
coefficients. Intercepted energy must then be partitioned among components in
proportion to their contribution to that extinction. Tests should cover empty
components, identical-component invariance, and mixtures with contrasting
optical properties.

## Dynamic Canopies And Model Coupling

Functional–structural plant models may change `PAD`, component identity, or
terrain occupancy between timesteps. Because directional responses depend on
the complete attenuation field, cache invalidation must remain strict.

Possible future strategies are:

- rebuild all responses after any structural change;
- rebuild only after a documented change threshold;
- identify affected spatial regions and incrementally update transport data;
- separate slowly changing geometry from rapidly changing optical properties.

The last option is especially useful: a change in reflectance should not force
geometric ray paths to be rebuilt, whereas a change in `PAD` generally changes
attenuation along every downstream path.

Coupling interfaces should expose arrays and units directly and should not
require the voxel solver to depend on a particular growth or ecosystem model.

## Spatial Resolution, Point Clouds, And Adaptive Grids

Point clouds should not be treated as physical scattering elements without an
explicit reconstruction assumption. Points do not directly provide plant
area, surface normals, or optical properties. Useful ingestion paths include:

- estimating `PAD` and leaf-angle distributions on the regular grid;
- reconstructing high-resolution solid voxels or surfaces;
- retaining uncertainty associated with occlusion and sampling density.

Non-uniform grids or octrees may reduce memory for sparse scenes, but they also
complicate traversal, periodic boundaries, response caching, GPU kernels, and
reproducibility. They should follow, rather than precede, a validated regular-
grid scattering implementation. A resolution study should demonstrate that
the additional complexity improves a target scientific output.

## Validation Strategy

Every extension should be validated at several levels.

### Analytical And Synthetic Tests

- homogeneous Beer–Lambert columns;
- transparent and fully absorbing limits;
- rotational symmetry under spherical LIDF;
- two-voxel single-scattering transfers;
- geometric-series multiple scattering with a known solution;
- horizontal and inclined Lambertian soil;
- open and periodic boundary energy closure;
- invariance under equivalent grid refinement where expected.

### Cross-Implementation Comparisons

- retain the historical Java voxel fixtures as provenance tests;
- keep the sorted-plane traversal as an independent CPU reference;
- compare selected homogeneous cases with RATP-style exchange calculations;
- use DART, FLIGHT, or a small Monte Carlo prototype for bounded external
  comparisons when anisotropy or complex multiple scattering is introduced.

External models should be used as scientific comparison points, not as golden
outputs whose undocumented behavior must be copied.

### Performance Measurements

Performance reports should separate:

- grid and response construction;
- one directional response;
- one meteorological step with an existing cache;
- many warmed timesteps;
- one-bounce and converged multiple scattering;
- CPU and GPU memory use;
- host–device transfer cost.

Large-scene benchmarks should record grid size, occupied fraction, ray density,
direction count, scattering iterations, hardware, Julia version, and backend.
Performance changes are acceptable only when reference comparisons and energy
balances still pass.

## Decisions To Make Before Extending The Public API

The following questions should be resolved before the corresponding feature is
made public:

1. Are optical properties attached to every voxel, to a material identifier,
   or to vegetation components mixed within a voxel?
2. Is isotropic scattering sufficient for the first target applications?
3. Must soil support only height fields, or arbitrary meshes from the start?
4. Which result names distinguish incident, intercepted, and absorbed energy
   without breaking the first-order workflow?
5. Which cache layers depend on geometry, `PAD`, optical properties, soil, and
   meteorology?
6. What numerical tolerance defines CPU/GPU agreement and scattering
   convergence?

## Reference Approaches

The proposed directions draw on several established radiative-transfer
approaches:

- Sinoquet et al. (2001), [RATP: a model for simulating the spatial
  distribution of radiation absorption, transpiration and photosynthesis
  within canopies](https://doi.org/10.1046/j.1365-3040.2001.00694.x):
  turbid voxels, exchange coefficients, and a radiosity-like linear system for
  isotropic PAR and NIR multiple scattering.
- Gastellu-Etchegorry et al. (1996), [Modeling radiative transfer in
  heterogeneous 3-D vegetation
  canopies](https://doi.org/10.1016/0034-4257(95)00253-7): discrete directional
  transport and iterative scattering in voxelized scenes.
- North (1996), [Three-dimensional forest light interaction model using a
  Monte Carlo method](https://doi.org/10.1109/36.508411): stochastic photon
  transport through volumetric crowns, woody components, and soil.
- Li et al. (2018), [VBRT: A novel voxel-based radiative transfer model for
  heterogeneous three-dimensional forest
  scenes](https://doi.org/10.1016/j.rse.2017.12.043): Monte Carlo path tracing
  through high-resolution solid voxels reconstructed from point clouds.

These references illustrate design alternatives. The recommended ArchimedLight
path is a deterministic, matrix-free, RATP-like solver first, followed by
additional angular detail and GPU acceleration after the conservation contracts
are stable.
