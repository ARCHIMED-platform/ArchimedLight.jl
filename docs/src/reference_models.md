# Model Files Reference

Model files define how ARCHIMED functional groups behave radiatively.
In the current Julia package, the light-critical part is the `Interception` block. Historical ARCHIMED files may also contain photosynthesis and stomatal-conductance blocks, but those are outside the scope of `ArchimedLight.jl`.

## One File Per Functional Group

Example:

```yaml
---
Group: coffee
Type:
  Metamer:
    Interception:
      model: Translucent
      transparency: 0.1
      optical_properties:
        PAR: 0.15
        NIR: 0.9
  Leaf:
    Interception:
      model: Translucent
      transparency: 0.0
      optical_properties:
        PAR: 0.15
        NIR: 0.9
...
```

The main structure is:

- `Group`: functional-group name, which must match the `.ops` metadata
- `Type`: mapping from component type names to their process blocks

## Matching Rules

The scene node `(group, type)` pair is resolved against the model files with these fallbacks:

- exact group and exact type
- exact group and wildcard type `*`
- wildcard group `*` with exact type
- wildcard group `*` with wildcard type `*`

This makes it easy to:

- declare one generic optical model for many synthetic scene types
- override only a few special cases

## The `Interception` Block

The current light runtime uses these keys:

```yaml
Interception:
  model: Translucent
  transparency: 0.0
  optical_properties:
    PAR: 0.15
    NIR: 0.9
```

### `model`

The standard historical model is `Translucent`.
Virtual sensors are declared either with:

```yaml
model: VirtualSensor
```

or

```yaml
sensor: true
```

Virtual sensors receive light but stay transparent and non-absorbing in the transfer logic.

### `transparency`

First-order transmission parameter used by the interception model.

### `optical_properties`

Waveband-specific scattering coefficients.
For the current light package:

- `PAR` and `NIR` are the important built-in bands
- extra shortwave bands are also supported

Example with one custom band:

```yaml
optical_properties:
  PAR: 0.15
  NIR: 0.90
  custom: 0.15
```

The optical coefficient is the scattered fraction. In the light budget integration step, absorptance is derived from it as approximately:

```text
absorptance = 1 - scattering_coefficient
```

when no more specific optical information is available.

![Reflectance, transmittance, absorptance](assets/optical_properties_reflectance_transmittance.jpg)

## Multiple Named Variants

Historical ARCHIMED files often store several named parameter sets under one process block and select one with `use`:

```yaml
Interception:
  use: Translucent_1
  Translucent_1:
    model: Translucent
    transparency: 0.0
    optical_properties:
      PAR: 0.15
      NIR: 0.90
  Translucent_2:
    model: Translucent
    transparency: 0.2
    optical_properties:
      PAR: 0.20
      NIR: 0.85
```

`ArchimedLight.jl` preserves that pattern and resolves the selected variant.

## Light Emitters

Emitter blocks are also supported in the type model:

```yaml
LightEmitter:
  model: LambertianEmitter
  radiance: 20.0
  gamma:
    PAR: 0.48
    NIR: 0.52
```

This is useful for artificial light sources or diagnostic scenes, though most canopy workflows are driven by the solar and sky forcing from the meteo file.

## Soil And Paving

Ground models are ordinary functional groups, often with a `plot_paving` hint:

```yaml
Group: pavement
Type:
  Cobblestone:
    Interception:
      model: Translucent
      transparency: 0
      optical_properties:
        PAR: 0.15
        NIR: 0.9
    plot_paving: 80
```

If `read_config` sees a positive `plot_paving` count and the scene does not already contain that ground type, it materializes paving tiles automatically.

## What Happens If A Band Is Missing

When a model omits a waveband in `optical_properties`, the runtime falls back to the matching global option:

- missing PAR coefficient -> `options.scattering_coeff_par`
- missing NIR coefficient -> `options.scattering_coeff_nir`

That is why the global scattering coefficients still matter even when most groups are fully parameterized.

## What Is Ignored By The Current Package

Historical model files may include:

- `Photosynthesis`
- `StomatalConductance`
- many energy-balance parameters

Those sections are valuable archival context, but they are not the active light API of `ArchimedLight.jl`.

## Interactive Equivalent

In an interactive workflow, you do not need YAML model files. The equivalent is
to build [`GroupModel`](@ref), [`TypeModel`](@ref), and
[`InterceptionModel`](@ref) objects directly in Julia, then normalize them with
`prepare_models`.

Typical pattern:

```julia
models = prepare_models([
    GroupModel(
        "coffee";
        types=OrderedDict(
            "Leaf" => TypeModel(
                interception=InterceptionModel(
                    model="Translucent",
                    transparency=0.0,
                    optical_properties=OpticalProperties(0.15, 0.90),
                ),
            ),
            "Metamer" => TypeModel(
                interception=InterceptionModel(
                    model="Translucent",
                    transparency=0.1,
                    optical_properties=OpticalProperties(0.15, 0.90),
                ),
            ),
        ),
    ),
])
```

The correspondences are:

- YAML `Group:` -> `GroupModel("coffee"; ...)`
- YAML `Type:` entries -> the `types=OrderedDict(...)` mapping
- YAML `Interception:` block -> `InterceptionModel(...)`
- YAML `LightEmitter:` block -> `EmitterModel(...)`

This means the matching semantics are unchanged. A node still resolves models by
its `(group, type)` pair; the only difference is whether the mapping came from
files or from Julia objects you built yourself.

### Interactive Wildcard Models

The wildcard mechanism also works when models are defined directly in Julia.

You can use `*` at two levels:

- wildcard group: `GroupModel("*"; ...)`
- wildcard type: `types=OrderedDict("*" => ...)`

The broadest fallback is:

```julia
models = prepare_models([
    GroupModel(
        "*";
        types=OrderedDict(
            "*" => TypeModel(
                interception=InterceptionModel(
                    model="Translucent",
                    transparency=0.0,
                    optical_properties=OpticalProperties(0.15, 0.30),
                ),
            ),
        ),
    ),
])
```

That means:

- any functional group matches `GroupModel("*", ...)`
- any component type inside that group matches the `"*"` type entry

This is useful when:

- building quick synthetic scenes
- writing tests that should not depend on many exact type names
- defining a safe default before adding more specific rules

You can also mix explicit and wildcard rules:

```julia
models = prepare_models([
    GroupModel(
        "coffee";
        types=OrderedDict(
            "Leaf" => TypeModel(
                interception=InterceptionModel(
                    model="Translucent",
                    transparency=0.0,
                    optical_properties=OpticalProperties(0.18, 0.88),
                ),
            ),
            "*" => TypeModel(
                interception=InterceptionModel(
                    model="Translucent",
                    transparency=0.0,
                    optical_properties=OpticalProperties(0.15, 0.90),
                ),
            ),
        ),
    ),
    GroupModel(
        "*";
        types=OrderedDict(
            "*" => TypeModel(
                interception=InterceptionModel(
                    model="Translucent",
                    transparency=0.0,
                    optical_properties=OpticalProperties(0.15, 0.30),
                ),
            ),
        ),
    ),
])
```

With that setup:

- `("coffee", "Leaf")` uses the explicit coffee leaf model
- `("coffee", "Stem")` falls back to the wildcard type inside the coffee group
- `("banana", "Leaf")` falls back to the global wildcard group and wildcard type

The precedence is the same as for YAML model files:

1. exact group and exact type
2. exact group and wildcard type
3. wildcard group and exact type
4. wildcard group and wildcard type
