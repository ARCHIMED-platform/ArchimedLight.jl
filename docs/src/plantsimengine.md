# Couple ArchimedLight With PlantSimEngine

Use this integration when one scene-scale light calculation must provide
irradiance to many physiological organs. `ArchimedLightModel` runs
`ArchimedLight` once on a `Scene` object, maps geometric components back to
stable botanical identities, and publishes the resulting variables directly
to the selected organ statuses. These are raw radiative-mesh-area values, not
leaf-area physiology inputs.

This is an ordinary PlantSimEngine model coupling. The light results are model
outputs owned by the scene application and stored on their destination
objects; they are not a mutable environment and do not require
`commit_environment!`. Likewise, there is no bookkeeping `LeafLightModel` or
`InternodeLightModel` to execute once per organ.

!!! note "Optional integration"
    PlantSimEngine is a weak dependency. Add it to the active environment and
    load both packages to activate the extension:

    ```julia
    using ArchimedLight
    using PlantSimEngine
    ```

## Declare the scene writer

Assume `light_sim` is an existing [`LightSimulation`](@ref) whose
`PlantGeom.SceneGeometry` retains source-owner identities. Construct the
PlantSimEngine kernel, then declare its destination objects with
`outputs_to`:

```julia
light_kernel = ArchimedLightModel(
    light_sim;
    output_schema=:coupling,
    output_group=:organs,
    par_energy_to_photon=4.57,
)

light_application = ModelSpec(
    light_kernel;
    name=:archimed_light,
    on=One(scale=:Scene),
    outputs_to=(
        organs=OutputTo(
            Many(
                scale=:Leaf,
                within=SceneScope(),
            );
            vars=archimed_light_outputs(:coupling),
            coverage=:exact,
        ),
    ),
)
```

`ArchimedLightModel` returns a model kernel, not a `ModelSpec`. The scenario
therefore keeps ownership of the application name, the `Scene` target, the
organ selector, and the cadence. The `output_group` keyword must match the name
of the corresponding `outputs_to` entry; both are `:organs` above.

When the physiological model is built from the exact MTG stored by the light
scene, include explicit radiative-to-botanical area adapters before leaf
physiology. The adapters are provided by PlantBiophysics:

```julia
model = CompositeModel(
    light_sim.scene.mtg;
    applications=(
        light_application,
        radiative_mesh_to_leaf_ppfd_application,
        radiative_mesh_to_leaf_shortwave_application,
        leaf_energy_balance_application,
        leaf_photosynthesis_application,
    ),
    environment=forcing,
)

simulation = run!(model; outputs=:all)
```

Every destination leaf must carry a finite positive `botanical_leaf_area`.
Map raw `aPPFD` and `radiative_mesh_area` to
`PlantBiophysics.RadiativeMeshToLeafPPFD`, and map raw `Ra_SW_f` and the same
area to `PlantBiophysics.RadiativeMeshToLeafShortwave`. Then map the distinct
outputs `aPPFD_leaf_mean` and `Ra_SW_f_leaf_mean` to FvCB and Monteith. Do not
connect ArchimedLight's same-named raw outputs directly to those leaf inputs.

For example, the PPFD boundary is declared explicitly as:

```julia
using PlantBiophysics

radiative_mesh_to_leaf_ppfd_application = ModelSpec(
    RadiativeMeshToLeafPPFD();
    name=:radiative_to_leaf_ppfd,
    on=Many(scale=:Leaf),
    inputs=(
        :aPPFD_radiative => One(
            within=Self(), application=:archimed_light, var=:aPPFD,
            policy=HoldLast(),
        ),
        :radiative_mesh_area => One(
            within=Self(), application=:archimed_light,
            var=:radiative_mesh_area, policy=HoldLast(),
        ),
    ),
)

leaf_photosynthesis_application = ModelSpec(
    Fvcb();
    name=:photosynthesis,
    on=Many(scale=:Leaf),
    inputs=(
        :aPPFD => One(
            within=Self(), application=:radiative_to_leaf_ppfd,
            var=:aPPFD_leaf_mean, policy=HoldLast(),
        ),
    ),
)
```

Use `RadiativeMeshToLeafShortwave` in the same way for `Ra_SW_f`, and map its
`Ra_SW_f_leaf_mean` output to Monteith.

The leaf applications remain ordinary PlantSimEngine applications on
`Many(scale=:Leaf)`. The compiler schedules the scene writer, then the area
adapters, then their physiology consumers. Application tuple order and MTG
traversal order do not define this coupling.

Leaf-scale energy balance and photosynthesis should use the canonical
`Many(scale=:Leaf)` destination above. If another radiative process genuinely
owns outputs on internodes, select those organs explicitly with
`Many(scale=(:Leaf, :Internode), within=SceneScope())`; exact coverage then
applies to both scales.

## Choose the output schema

`archimed_light_outputs` returns the declaration that must be passed to
`OutputTo(...; vars=...)`. Two schemas are available.

The default `:coupling` schema keeps the hot publication path compact:

| Variable | Meaning | Unit |
| --- | --- | --- |
| `Ri_PAR_f` | total intercepted PAR irradiance | `W m[radiative]^-2` |
| `Ri_NIR_f` | total intercepted NIR irradiance | `W m[radiative]^-2` |
| `Ra_PAR_f` | total absorbed PAR irradiance | `W m[radiative]^-2` |
| `Ra_NIR_f` | total absorbed NIR irradiance | `W m[radiative]^-2` |
| `Ra_SW_f` | `Ra_PAR_f + Ra_NIR_f` | `W m[radiative]^-2` |
| `aPPFD` | absorbed photosynthetic photon flux density | `μmol m[radiative]^-2 s^-1` |
| `radiative_mesh_area` | area used by the light solver for normalization | `m[radiative]^2` |

`radiative_mesh_area` is a radiative discretization property. It is not a
canonical botanical leaf area or a projected area. `aPPFD` is computed as
`Ra_PAR_f * par_energy_to_photon`; `4.57` is the default broadband PAR
conversion in `μmol J^-1`.

The extension declares scientific contracts only for `aPPFD`, `Ra_SW_f`, and
`radiative_mesh_area`. The `Ri_*` and component PAR/NIR diagnostics remain
available but are deliberately not assigned contracts by analogy.

For either contracted flux `F`, the PlantBiophysics boundary computes

```math
F_{leaf} = F_{radiative}
\times \frac{A_{radiative}}{A_{botanical}},
```

so `F_radiative * radiative_mesh_area == F_leaf * botanical_leaf_area`. This
equality is the conservation check for manual reference paths as well as the
automatic distributed path.

The `:full` schema adds all initial (`*_0_*`) and scattering-inclusive PAR/NIR
incident and absorbed variables:

- `Ri_PAR_0_f`, `Ri_NIR_0_f`, `Ri_PAR_f`, and `Ri_NIR_f`;
- `Ri_PAR_0_q`, `Ri_NIR_0_q`, `Ri_PAR_q`, and `Ri_NIR_q`;
- `Ra_PAR_0_f`, `Ra_NIR_0_f`, `Ra_PAR_f`, and `Ra_NIR_f`;
- `Ra_PAR_0_q`, `Ra_NIR_0_q`, `Ra_PAR_q`, and `Ra_NIR_q`; and
- the three derived/coupling columns `Ra_SW_f`, `aPPFD`, and
  `radiative_mesh_area`.

Use the same schema in both places:

```julia
full_light = ArchimedLightModel(light_sim; output_schema=:full)

full_application = ModelSpec(
    full_light;
    name=:archimed_light,
    on=One(scale=:Scene),
    outputs_to=(
        organs=OutputTo(
            Many(scale=:Leaf, within=SceneScope());
            vars=archimed_light_outputs(:full),
            coverage=:exact,
        ),
    ),
)
```

Passing mismatched schemas is rejected before the light calculation.

!!! warning "No `sky_fraction` output"
    Neither schema publishes `sky_fraction`. Models that scientifically require
    a sky-view factor must declare and obtain it from an explicit producer; it
    must not be inferred from irradiance or silently substituted by this
    adapter.

## Supply the radiative forcing

The model reads five canonical values from the sampled PlantSimEngine
environment:

| Variable | Contract |
| --- | --- |
| `sun_azimuth_deg` | finite angle in `[0, 360)` degrees |
| `sun_elevation_deg` | finite angle in `[-90, 90]` degrees |
| `Ri_PAR_f` | non-negative incident PAR irradiance in `W m^-2` |
| `Ri_NIR_f` | non-negative incident NIR irradiance in `W m^-2` |
| `direct_fraction` | direct fraction in `[0, 1]` |

ArchimedLight uses `1 - direct_fraction` as the diffuse fraction. Every global
forcing row must also carry a finite positive `duration`, either as a number of
seconds or as a fixed `Dates` period such as `Hour(1)`. PlantSimEngine adds
this timeline field automatically to the sampled model-facing environment, so
do not declare or remap a second duration input for `ArchimedLightModel`.

For example, one forcing row can be written as:

```julia
using Dates

forcing = (
    sun_azimuth_deg=180.0,
    sun_elevation_deg=45.0,
    Ri_PAR_f=500.0,
    Ri_NIR_f=450.0,
    direct_fraction=0.7,
    duration=Hour(1),
)
```

All values are validated before the numerical light solver runs and before any
organ status is modified.

!!! compat "Same-cadence contract in the first integration"
    Use one ArchimedLight execution per forcing interval. This first version
    does not define a cross-rate aggregation policy for solar azimuth,
    elevation, or direct fraction. In particular, averaging angular forcing
    over a coarser light cadence would not have a well-defined physical meaning
    here.

## Map geometry to destination objects

PlantGeom stores each simulated component's owner as a composite key:
`(source_instance_id, source_node_id)`. The extension compiles those keys to
PlantSimEngine `ObjectId`s once per stable scene generation. It never pairs a
result row with an organ by vector position, MTG traversal order, or raw node
number alone.

Three mapping modes are available:

1. With neither `source_roots` nor `object_resolver`, the PlantSimEngine model
   must be built from the exact `scene.mtg` retained by `light_sim.scene`.
2. For a radiative scene assembled from other MTG roots, pass
   `source_roots=instance_id => exact_registered_root` mappings. The value can
   also be a function of the current `RunContext` returning that mapping.
3. For a non-MTG ownership scheme, pass
   `object_resolver=key -> registered_object`. PlantGeom accepts a registered
   object, status, exact MTG node, `ObjectId`, or raw registered id and validates
   the result against the runtime model.

For example:

```julia
light_kernel = ArchimedLightModel(
    light_sim;
    source_roots=Dict(
        1 => exact_first_plant_root,
        2 => exact_second_plant_root,
    ),
)
```

`source_roots` and `object_resolver` are mutually exclusive. Never derive an
`ObjectId` from `source_node_id` alone: the same raw MTG node id may occur in
several source instances.

The `source_roots` mode also makes the publication boundary explicit for an
assembled scene. A ground surface or another source instance omitted from the
mapping remains in the radiative calculation and can still shade or scatter
light, but it is not written to an organ status. More generally, any geometric
component outside the selected target set remains physically active and is
not published.

## Exact coverage and atomic publication

Only `coverage=:exact` is supported. Every current object selected by
`OutputTo` must map to at least one simulated radiative component. Several
components or owner keys may map to the same destination; energy quantities
are summed and flux densities are weighted by their radiative mesh areas.

If a selected `Leaf` or `Internode` has no mapped geometry, the extension raises
an error before `run_light` and before writing any destination. It does not
retain an old value, fill in zero, or rely on row order. This matters after
organ creation: either rebuild the scene so that the new target has geometry,
or deliberately keep that object outside the destination selector.

The final identified assignment is also atomic with respect to validation:
unknown, duplicate, extra, or missing destination identities and missing
declared columns are rejected before status mutation.

## Refresh the scene at lifecycle barriers

Steady-state runs reuse the compiled owner mapping and typed result buffers.
When PlantSimEngine reports a model or environment revision but the
`LightSimulation` scene was not explicitly refreshed, the adapter needs a
coordinated scene provider:

```julia
light_kernel = ArchimedLightModel(
    light_sim;
    scene_provider=context -> rebuild_light_scene(context),
    source_roots=context -> current_source_roots(context),
)
```

`scene_provider(context)` must return the `PlantGeom.SceneGeometry` valid for
the current runtime topology. If geometry changed, the callback is responsible
for rebuilding it through PlantGeom. Calling
`update_scene!(light_sim, light_sim.scene)` by itself only invalidates
ArchimedLight caches; it does not reconstruct a merged mesh from newly created
organs.

For an unrelated PlantSimEngine revision, a provider may return the still-valid
current scene object at the same scene version to acknowledge that no geometric
rebuild was required; this preserves the expensive radiative cache while owner
and target bindings are refreshed. A version increase on the same
`SceneGeometry` object is different: it does not prove that the merged mesh was
rebuilt, so it requires `scene_provider` to return or explicitly affirm the
geometry prepared for that generation. Without a provider, such an
uncoordinated lifecycle change is a hard error instead of a run against stale
geometry.

Explicitly installing a distinct rebuilt `SceneGeometry` with
[`update_scene!`](@ref) is accepted and triggers mapping recompilation without
calling the provider.

One kernel instance belongs to one PlantSimEngine runtime model and one
execution `Scene`. Construct a distinct `ArchimedLightModel` for each scene
application instead of reusing the same mutable kernel across models or
`Scene` objects.

## Performance and diagnostics

On an unchanged scene, the adapter uses constant-time revision checks and
reuses its component-to-target mapping, aggregation columns, and target-id
buffer. It obtains the current PlantSimEngine `OutputTargets` view on every
invocation and never retains that lifecycle-scoped view across steps. Custom
wrappers should follow the same rule.

Use PlantSimEngine diagnostics to inspect the compiled writer, both area
boundaries, and their consumers:

```julia
Diagnostics.explain_applications(model)
Diagnostics.explain_writers(model)
Diagnostics.explain_bindings(model)
Diagnostics.explain_schedule(model)
```

Run with `outputs=:all` when both the raw ArchimedLight history and the current
converted physiology state are needed. The final `aPPFD` and `Ra_SW_f` fields
on the leaf are the converted values; read the raw values from the
`:archimed_light` output stream:

```julia
simulation = run!(model; outputs=:all)
leaf_state = final_state(simulation, :leaf_42)
raw_ppfd = last(outputs(simulation)[
    (:archimed_light, ObjectId(:leaf_42), :aPPFD)
])[2]
raw_shortwave = last(outputs(simulation)[
    (:archimed_light, ObjectId(:leaf_42), :Ra_SW_f)
])[2]
(
    raw_ppfd=raw_ppfd,
    leaf_ppfd=leaf_state.aPPFD_leaf_mean,
    raw_shortwave=raw_shortwave,
    leaf_shortwave=leaf_state.Ra_SW_f_leaf_mean,
)
```

For large canopies, prefer the `:coupling` schema unless the full diagnostic
budget is actually consumed. Recompile the mapping only at a real scene or
topology lifecycle barrier; do not rebuild selectors or perform an MTG join in
the per-step path.
