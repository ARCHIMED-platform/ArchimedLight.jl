# Voxel Scattering Theory And Optical Parameters

The voxel solver separates three processes that are easy to conflate:

1. **interception** removes radiation from a ray according to voxel geometry,
   `PAD`, the projection coefficient `G`, and path length;
2. **scattering** re-emits a fraction of intercepted radiation into new
   directions;
3. **absorption** is the remaining terminal fraction retained by foliage.

This distinction matters because intercepted energy can be intercepted again
after scattering. Total intercepted energy is therefore a transfer diagnostic,
not a terminal term in the energy budget.

```@contents
Pages = ["theory_voxel_scattering.md"]
Depth = 3
```

## Optical Contract

For each waveband, voxel foliage has leaf reflectance ``\rho`` and leaf
transmittance ``\tau``. The single-scattering albedo and absorptance are

```math
\omega = \rho + \tau, \qquad \alpha = 1 - \omega.
```

All four quantities are dimensionless. ArchimedLight enforces
``0 \leq \rho \leq 1``, ``0 \leq \tau \leq 1``, and
``\rho + \tau \leq 1`` at construction. The compatibility keyword
`scattering_par` or `scattering_nir` supplies ``\omega`` directly and splits it
equally between ``\rho`` and ``\tau``. That equal split preserves the total
coefficient used by this isotropic solver; it is not a claim about measured
leaf reflection and transmission.

Leaf transmittance is not canopy porosity. Gaps are already represented by
`PAD` through Beer–Lambert attenuation:

```math
T_\mathrm{path}=\exp(-G\,PAD\,\ell).
```

An empty voxel therefore intercepts nothing, irrespective of values stored in
its optical arrays. Conversely, optical transmittance acts only after plant
material has intercepted radiation.

`VoxelOpticalProperties` accepts scalar coefficients without allocating dense
fields, or arrays aligned exactly with `grid.pad`. The first implementation
uses dense effective coefficients rather than a material-index table. Mixed
species, senescent foliage, and woody area must therefore be represented by
effective per-voxel `PAD`, reflectance, and transmittance fields. An
area-weighted optical mixture is a practical approximation when the components
share the same geometric projection model; distinct geometry or leaf-angle
distributions require a richer future representation.

## Operational Wavebands And Spectral Averaging

The optical convention is:

- PAR: 400–700 nm;
- NIR: 700–2500 nm.

All current voxel outputs are energy quantities (`W` or `J`), so spectra must
be reduced with incident spectral irradiance weighting. For a band ``b``:

```math
\bar{\rho}_b =
\frac{\int_{\lambda_1}^{\lambda_2} E_\lambda\rho_\lambda\,d\lambda}
     {\int_{\lambda_1}^{\lambda_2} E_\lambda\,d\lambda},\qquad
\bar{\tau}_b =
\frac{\int_{\lambda_1}^{\lambda_2} E_\lambda\tau_\lambda\,d\lambda}
     {\int_{\lambda_1}^{\lambda_2} E_\lambda\,d\lambda}.
```

Do not mix photon-weighted PAR coefficients with energy-weighted model
outputs. If the source data are photon spectral densities or PPFD, convert
them to energy spectral density before averaging for this solver.

!!! warning "Meteorological NIR must match the optical convention"
    When only total shortwave radiation is available, ArchimedLight can derive
    its operational NIR forcing as the non-PAR fraction of shortwave radiation.
    That broad forcing is not guaranteed to be an exact 700–2500 nm integral.
    For calibrated spectral work, supply `Ri_PAR_f` and `Ri_NIR_f` derived over
    the same intervals used for the optical coefficients.

Useful coefficient sources include species-specific integrating-sphere
measurements, the [LOPEX93 database](https://teledetection.ipgp.jussieu.fr/opticleaf/lopex.htm),
the ANGERS leaf optical database, and spectra simulated with
[PROSPECT-D](https://doi.org/10.1016/j.rse.2017.03.004). Store the species,
leaf state, measurement geometry, wavelength grid, averaging spectrum, and
date with calibrated inputs. The convenience constructors
`voxel_java_parity_optics` (`0.15` PAR, `0.90` NIR) and
`voxel_generic_green_leaf_optics` (`0.17`, `0.85`) are reproducibility and
teaching values, not species calibrations.

## Internal-Source Transport

First-order rays begin at the canopy boundary. Scattered radiation instead
starts inside an occupied voxel, so it uses a dedicated transport operator.
For each occupied source voxel, ArchimedLight traces the full-sphere angular
quadrature from the voxel centre. This centre-source approximation is
deterministic, includes the remaining path through the emitting voxel, and
therefore permits within-cell re-interception. It does not integrate over all
possible sub-voxel emission positions; convergence with voxel refinement
should be checked when this approximation matters.

The angular quadrature pairs every diffuse turtle direction with its opposite.
The normalized weights sum to one, and each hemisphere receives one half of an
isotropic foliage source. Sun sectors are excluded because scattered emission
is independent of the instantaneous direct-beam direction.

The matrix-free operator applies Beer–Lambert attenuation along cached DDA
paths. For every application:

```math
E_\mathrm{emitted} = E_\mathrm{reintercepted}
                    + E_\mathrm{top}
                    + E_\mathrm{side}
                    + E_\mathrm{bottom},
```

with ground absorption replacing bottom escape when a ground boundary is
present. A tiny dense exchange matrix exists only as a test oracle; production
code never materializes an ``N_\mathrm{voxel}^2`` matrix.

## Iteration And Convergence

For scattering order ``n`` and one band:

```math
A^{(n)} = (1-\omega)I^{(n)},\qquad
S^{(n)} = \omega I^{(n)},\qquad
I^{(n+1)} = \mathcal{T}(S^{(n)}).
```

The solver stops when the newly intercepted order is smaller than
`scattering_stop_ratio` times first-order foliage interception, or at the
floating-point floor. `scattering_max_iter` is a hard limit. If that limit is
reached first, `converged` is false and the latest unprocessed order remains in
`unresolved_energy`; it is never silently treated as absorption.

The terminal energy invariant is checked independently for PAR and NIR:

```math
E_\mathrm{external}+E_\mathrm{injected} =
E_\mathrm{foliage\ absorbed}+E_\mathrm{ground\ absorbed}+
E_\mathrm{top}+E_\mathrm{side}+E_\mathrm{bottom}+E_\mathrm{unresolved}.
```

`total_intercepted_energy` and `scattered_energy` are intentionally absent
from the right-hand side: they describe internal transfers that would otherwise
be double-counted.

The result retains the scalar intercepted total for every processed order in
`order_intercepted_energy`. This ``O(N_\mathrm{iteration})`` diagnostic is kept
in normal results because it makes convergence auditable without retaining
per-order voxel arrays. Reusing mutable grid-sized workspaces across returned
band or timestep results is deferred: the current implementation favors clear
result ownership and thread-safe read-only caches. A future in-place workspace
API may reduce the remaining warmed allocations without changing equations.

## Ground Boundary

`VoxelGroundOptics` represents a horizontal Lambertian receptor. Energy that
reaches the lower boundary is split into absorbed and reflected fractions.
Reflected energy is redistributed over upward directions with weights
proportional to quadrature solid angle times ``\cos\theta``.

Without a ground object, lower-boundary radiation is reported as bottom
escape. With a black ground it is reported as ground absorption. With a
reflective ground it can be intercepted by foliage again or leave through the
top or an open side. These outcomes remain separate in the result.

## Scope And Limitations

The CPU reference solver supports spherical LIDF, isotropic foliage scattering,
open or horizontally periodic regular grids, scalar or dense effective optical
fields, and an optional Lambertian ground. It does not yet model specular
reflection, polarized light, fluorescence, wavelength-resolved transport,
anisotropic leaf phase functions, separate material identities, adaptive
voxels, or GPU execution.

The isotropic assumption uses only ``\rho+\tau`` during redistribution. The
separate coefficients are retained so the API and input provenance do not
erase measured information and a later anisotropic model can distinguish the
two hemispheres.

Transport caches are constructed locally and treated as read-only after
creation. Independent tasks may safely read one cache concurrently only if
their output arrays are distinct; callers must not mutate the cache dictionaries
or their stored paths. The series API owns its cache within one call and does
not expose shared mutable work buffers.

For related numerical context, see the voxel exchange treatment in
[RATP](https://doi.org/10.1046/j.1365-3040.2001.00694.x) and the iterative
discrete-direction transport in
[DART](https://doi.org/10.1016/0034-4257(95)00253-7).

Analytical one- and two-voxel cases plus a tiny dense deterministic exchange
matrix provide the current independent numerical oracles. A Monte Carlo oracle
was not needed for the isotropic centre-source implementation; it remains a
useful future challenge test if sub-voxel source sampling or anisotropic phase
functions are introduced.
