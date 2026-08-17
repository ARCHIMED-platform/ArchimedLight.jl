# Voxel Light Interception And Scattering

Use the voxel workflow when vegetation is represented by plant area density
(`PAD`) in a regular 3D grid rather than by triangles. It reuses ArchimedLight's
meteorology, sky discretization, and directional irradiances, then applies
Beer–Lambert attenuation inside each voxel, with optional deterministic
multiple scattering.

```@contents
Pages = ["voxel_interception.md"]
Depth = 3
```

## Minimal In-Memory Run

This example builds a small homogeneous canopy. Its values are pedagogical;
they are not a calibrated canopy model.

```@example voxel
using ArchimedLight

grid = VoxelGrid(
    (0.0, 0.0, 0.0),
    (2.0, 2.0, 3.0),
    fill(0.4, 2, 2, 3),
)

options = LightOptions(
    turtle_sectors=6,
    all_in_turtle=true,
    scattering=false,
    radiation_timestep_minutes=5.0,
)

row = (
    step_duration=3600.0,
    Ri_PAR_f=100.0,
    Ri_NIR_f=50.0,
    sun_azimuth=180.0,
    sun_elevation=35.0,
    direct_fraction=0.5,
)

backend = VoxelCPUBackend(
    rays_per_voxel=4,
    boundary=:periodic,
    traversal=:dda,
)
step = run_voxel_light_step(grid, row, options; backend=backend)
```

The first-order result contains one 3D array per waveband and output unit. The
energy balance includes light that reaches the bottom or, with open horizontal
boundaries, leaves through a side:

```@example voxel
result = step.first_order
(
    par_intercepted=sum(result.par_energy),
    par_escaped=result.par_escaped_energy,
    par_incoming=result.par_incoming_energy,
    balance_closes=isapprox(
        sum(result.par_energy) + result.par_escaped_energy,
        result.par_incoming_energy;
        rtol=1e-10,
    ),
)
```

## Grid And Direction Conventions

`VoxelGrid(minimum, maximum, pad)` defines a regular Cartesian grid in metres.
The shape of `pad` is `(nx, ny, nz)`, and `pad[i, j, k]` is expressed in
`m² plant m⁻³`. Julia uses 1-based indices. Historical `.vox` rows are 0-based
and use `k` as the fastest-varying index; conversion is explicit at the file
boundary.

`trace_voxel_ray` accepts a **propagation direction**: the direction in which
light travels. It is normalized internally and must have a strictly negative
vertical component. Directions stored in a `TurtleGrid` point upward toward
the light source, so the voxel pipeline negates them before tracing.

```@example voxel
path = trace_voxel_ray(
    grid,
    (0.5, 0.5, grid.maximum[3]),
    (0.5, 0.0, -1.0);
    boundary=:periodic,
    algorithm=:dda,
)
[(segment.index, segment.length) for segment in path]
```

The public boundary modes are:

- `:periodic`: horizontal exits wrap to the opposite side while retaining the
  current ray intensity;
- `:open`: the ray leaves permanently at the first horizontal boundary;
- `:java_nontoric`: a diagnostic compatibility mode that reproduces the
  historical Java wrap-and-reset behavior.

The production traversal is `:dda`. `:reference` is a correct sorted-plane
implementation used to validate the DDA. `:java_reference` preserves known
historical wall-selection behavior and is intended only for parity work.

An origin exactly on the top or a horizontal grid plane is classified using a
small forward probe in the propagation direction. Simultaneous edge and corner
crossings advance every affected axis and do not produce zero-length segments.

## Beer–Lambert Attenuation And Units

For a path of length ``\ell`` through a voxel, the transmitted fraction is

```math
T = \exp(-G\,PAD\,\ell).
```

The initial implementation supports spherical leaf-angle distributions only,
with `G = 0.5` by default. `VoxelGrid` rejects negative and non-finite `PAD`
values. A zero-`PAD` voxel intercepts no light and its plant-area-normalized
flux is reported as zero.

| Julia field | Historical `.vox` column | Unit |
| --- | --- | --- |
| `par_energy` | `Ri_PAR_0_q` | `J voxel⁻¹ step⁻¹` |
| `nir_energy` | `Ri_NIR_0_q` | `J voxel⁻¹ step⁻¹` |
| `par_flux` | `Ri_PAR_0_f` | `W m⁻² plant` |
| `nir_flux` | `Ri_NIR_0_f` | `W m⁻² plant` |

`par_incoming_energy` and `nir_incoming_energy` are incident energies over the
whole horizontal domain. The corresponding `*_escaped_energy` fields close the
energy balance. `*_injected_energy` is zero for the public boundary modes and
records only the artificial reinjection required by `:java_nontoric`.

## Add Multiple Scattering

Voxel scattering requires explicit optical properties. It deliberately does
not inherit the triangle-model fallbacks (`0.15` PAR and `0.30` NIR), because
those values may not describe the voxelized material.

```@example voxel
optics = VoxelOpticalProperties(
    grid;
    par_reflectance=0.08,
    par_transmittance=0.07,
    nir_reflectance=0.45,
    nir_transmittance=0.40,
)

ground = VoxelGroundOptics(
    grid;
    par_reflectance=0.10,
    nir_reflectance=0.30,
)

scattering_options = LightOptions(
    options;
    scattering=true,
    scattering_max_iter=30,
    scattering_stop_ratio=1e-8,
)

scattered_step = run_voxel_light_step(
    grid,
    row,
    scattering_options;
    backend=backend,
    optics=optics,
    ground=ground,
)
```

Reflectance and transmittance are leaf-scale fractions of intercepted energy.
Their sum is the single-scattering albedo ``\omega``; the foliage absorptance
is ``1-\omega``. Geometric gaps remain controlled by `PAD` and are not leaf
transmittance.

The complete result keeps first-order interception and multiple-scattering
outputs separate:

```@example voxel
par = scattered_step.scattering.par
(
    first_order_intercepted=sum(scattered_step.first_order.par_energy),
    added_intercepted=sum(par.added_intercepted_energy),
    foliage_absorbed=sum(par.absorbed_energy),
    ground_absorbed=par.ground_absorbed_energy,
    top_escape=par.escaped_top_energy,
    side_escape=par.escaped_side_energy,
    bottom_escape=par.escaped_bottom_energy,
    unresolved=par.unresolved_energy,
    iterations=par.iterations,
    converged=par.converged,
    relative_residual=par.relative_balance_residual,
)
```

`total_intercepted_energy` may exceed first-order interception because the same
energy can encounter foliage more than once. Use `absorbed_energy`, boundary
escapes, ground absorption, and `unresolved_energy` for the terminal budget.
The arrays called `absorbed_flux` use plant area as the reference and have units
`W m⁻² plant`.

### One Order Versus Convergence

One iteration is useful as a diagnostic, but it is not automatically a
converged solution:

```@example voxel
one_order_options = LightOptions(
    scattering_options;
    scattering_max_iter=1,
    scattering_stop_ratio=0.0,
)
one_order = run_voxel_light_step(
    grid,
    row,
    one_order_options;
    backend=backend,
    optics=optics,
    ground=ground,
)

(
    first_order_PAR=sum(step.first_order.par_energy),
    one_order_PAR_absorbed=sum(one_order.scattering.par.absorbed_energy),
    converged_PAR_absorbed=sum(scattered_step.scattering.par.absorbed_energy),
    one_order_converged=one_order.scattering.par.converged,
    converged_iterations=scattered_step.scattering.par.iterations,
    converged_PAR_added=sum(scattered_step.scattering.par.added_intercepted_energy),
    converged_NIR_added=sum(scattered_step.scattering.nir.added_intercepted_energy),
)
```

NIR usually has a larger single-scattering albedo than PAR, so its redistributed
fraction can be much larger. That expectation is qualitative, not a substitute
for species- and condition-specific optical data. See
[Voxel Scattering Theory And Optical Parameters](theory_voxel_scattering.md)
for spectral averaging, parameter sources, equations, and limitations.

## Historical `.vox` Files And Series

`read_voxel_grid` accepts `PAD` or `PadBVTotal` in historical files and ignores
additional result columns. Every 0-based grid coordinate must occur exactly
once. `write_voxel_grid` can write either the grid alone or the grid plus a
`VoxelFirstOrderResult`.

```@example voxel
vox_path = joinpath(mktempdir(), "canopy.vox")
write_voxel_grid(vox_path, grid, step.first_order)
restored = read_voxel_grid(vox_path)

rows = [
    row,
    merge(row, (Ri_PAR_f=80.0, Ri_NIR_f=40.0, sun_azimuth=200.0)),
]
series = run_voxel_light_series(
    restored,
    rows,
    scattering_options;
    backend=backend,
    optics=VoxelOpticalProperties(
        restored;
        par_reflectance=0.08,
        par_transmittance=0.07,
        nir_reflectance=0.45,
        nir_transmittance=0.40,
    ),
)

(
    grid_size=size(restored),
    steps=length(series),
    second_step_par_absorbed=sum(series[2].scattering.par.absorbed_energy),
)
```

Directional geometric responses and internal scattering paths are cached
across a series. The cache key
includes grid bounds, dimensions, all `PAD` values, the normalized direction,
`G`, boundary mode, traversal, and ray sampling settings. A changed grid or
backend cannot silently reuse an incompatible explicit cache. Optical albedos
are applied during iteration rather than baked into geometric paths, so changing
only optics does not require rebuilding those paths.

## Current Limitations

!!! warning
    The voxel workflow currently supports spherical LIDF, isotropic foliage
    scattering, effective optical fields on regular Cartesian voxels, an
    optional Lambertian ground, and CPU execution. Specular or anisotropic
    scattering, separate material identities, adaptive grids, wavelength-
    resolved transport, and GPU execution are not implemented. Internal
    emission uses the voxel centre as a deterministic sub-voxel approximation.

The voxel workflow is intentionally separate from the triangle-based
`LightSimulation` / `run_light` API. It reuses the same sky and meteorological
stages without pretending that a voxel is a triangulated scene component.

## Public API

The exported workflow is deliberately small:

```julia
VoxelGrid(minimum, maximum, pad)
VoxelCPUBackend(; rays_per_voxel, g, boundary, traversal, lidf)
VoxelOpticalProperties(grid; par_reflectance, par_transmittance,
                       nir_reflectance, nir_transmittance)
VoxelGroundOptics(grid; par_reflectance, nir_reflectance)
read_voxel_grid(path)
write_voxel_grid(path, grid, result=nothing)
trace_voxel_ray(grid, origin, direction; boundary, algorithm)
prepare_voxel_responses(grid, turtle, backend)
compute_voxel_first_order(grid, turtle, fluxes, backend; duration_seconds, cache)
prepare_voxel_scattering_quadrature(turtle)
prepare_voxel_scattering_transport(grid, quadrature, backend)
apply_voxel_scattering_transport(grid, emitted, quadrature, backend; ...)
compute_voxel_scattering(grid, turtle, first, optics, options, backend; ...)
run_voxel_light_step(grid, meteo_row, options; backend, optics, ground)
run_voxel_light_series(grid, meteo, options; backend, optics, ground)
```

The result and diagnostic containers are `VoxelRaySegment`, `VoxelRayPath`,
`VoxelDirectionResponse`, `VoxelResponseCache`, `VoxelFirstOrderResult`,
`VoxelScatteringResult`, and `VoxelLightStepResult`. Their generated docstrings are listed in the
[Public API](api.md) page.
