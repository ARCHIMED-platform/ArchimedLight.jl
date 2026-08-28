# Changelog

## Unreleased

### Breaking changes

- The PlantSimEngine extension now declares `aPPFD` and `Ra_SW_f` as rates per
  radiative mesh area and `radiative_mesh_area` as an organ area. Direct
  coupling of those raw flux densities to contracted leaf-area physiology
  inputs is rejected.
- PlantBiophysics workflows must provide a finite positive
  `botanical_leaf_area` and use `RadiativeMeshToLeafPPFD` or
  `RadiativeMeshToLeafShortwave`. The boundary preserves absorbed quantity:
  `raw_flux * radiative_mesh_area == leaf_flux * botanical_leaf_area`.
- Raw output field names and the `:coupling` and `:full` schemas are unchanged.
  No contract was added to the `Ri_*` or component PAR/NIR diagnostic fields.

### Fixed

- MTG attachment now guarantees that `Ra_SW_f` means absorbed PAR plus absorbed
  NIR. The historical `names=Dict(:absorbed_nir_flux => :Ra_SW_f)` call remains
  accepted and produces that corrected sum with a deprecation warning; the new
  `:absorbed_shortwave_flux` selector expresses the same request directly.
