# Meteo And Forcing Reference

The light solver does not fundamentally need a meteo *file*. It needs a
time-dependent forcing description that tells it:

- which time interval is being simulated
- how much shortwave radiation arrives
- how that radiation is split into PAR, NIR, and any extra bands
- how much is direct beam and how much is diffuse sky
- where the sun is, or enough information to reconstruct it

Historically that information often comes from a semicolon-separated meteo CSV.
In dynamic workflows, the same semantics can be provided directly through a
`PlantMeteo.TimeStepTable` or even a single [`SkyState`](@ref).

## What Matters In Meteo For The Light Solver

Whatever the source format, the runtime ultimately wants one or more forcing
steps that define:

- a start and end time
- a usable radiative description
- enough site information to reconstruct solar geometry when needed

The rest of this page is organized around those semantics rather than around the
CSV format.

## Minimal Schema

For most users, a meteo table should provide these fields:

| Need | Accepted columns / metadata | Units | Required? |
| --- | --- | --- | --- |
| Date | `date`, a `DateTime`-like `date`, or `dayofyear` (+ optional `year`) | date or datetime | required only for reconstructed sun/clearness paths |
| Time interval | `hour_start` + `hour_end`, or `duration`/`step_duration` | time strings, seconds, or Julia `Dates.Period` | a positive duration is always required; start time is required only for reconstructed sun/clearness paths |
| Site latitude | `latitude` or `lat` as a column or metadata | degrees | yes, unless sun position is explicit |
| PAR irradiance | `RI_PAR_f`, `Ri_PAR_f`, `PAR`, or `par` | `W m^-2` | one radiation description is required |
| NIR irradiance | `RI_NIR_f`, `Ri_NIR_f`, `NIR`, or `nir` | `W m^-2` | optional if total shortwave is given |
| Total shortwave | `RI_SW_f`, `Ri_SW_f`, `Rg`, `rg`, `sw_global`, or `global` | `W m^-2` | alternative to PAR/NIR |
| Atmospheric clearness | `clearness` or `Kt` | fraction | alternative forcing input |
| Sun position | `sun_azimuth`/`sun_azimut` and `sun_elevation` | degrees | optional |
| Direct fraction | `direct_fraction`, `fDIR_SW`, or `Fd` | fraction | optional |
| Direct/diffuse SW | `Ri_SW_f_direct`/`RI_SW_f_direct` and `Ri_SW_f_diffuse`/`RI_SW_f_diffuse` | `W m^-2` | optional partition or total-SW path |

`run_light` accepts `PlantMeteo.TimeStepTable` and generic
Tables.jl-compatible inputs such as vectors of named tuples or DataFrames.
Use `summarize_meteo(meteo)` to see what ArchimedLight detected, and
`check_meteo(meteo)` to diagnose missing columns before a simulation.

### Canonical Names And Compatibility Aliases

The canonical in-memory names follow PlantMeteo and use `Ri_`: `Ri_SW_f`,
`Ri_PAR_f`, and `Ri_NIR_f`. ArchimedLight continues to accept the historical
file spellings `RI_SW_f`, `RI_PAR_f`, and `RI_NIR_f`, plus the aliases listed in
the table above. New Julia code should prefer the canonical `Ri_*` names.

Aliases are resolved by value rather than by column order. If a PlantMeteo
`Atmosphere` contributes its default `Ri_PAR_f=Inf` while a legacy
`RI_PAR_f` column contains a finite measurement, the finite measurement is
used. `missing`, blank strings, `NaN`, and positive or negative `Inf` are treated as
unavailable sentinels. If two finite aliases for the same variable disagree,
validation returns an error naming both columns and values.

Historical `use` metadata remains accepted as a hint for file compatibility,
but it does not suppress finite row values or force a broken input path. The
resolver chooses the usable values actually present.

Explicit `sun_azimuth`/`sun_azimut` and `sun_elevation` values are always in
degrees. Small degree values are not reinterpreted as radians.

### Radiation Derivation Rules

The resolver obtains one consistent SW/PAR/NIR forcing using these rules:

1. Preserve a finite total `Ri_SW_f` when provided.
2. Otherwise derive SW from finite PAR + NIR, direct + diffuse SW components,
   one direct/diffuse component plus `direct_fraction`, PAR alone, NIR alone,
   or finally `clearness`, in that order.
3. With total SW and one measured band, derive the other band as the residual.
4. With total SW and no measured bands, use the ARCHIMED 0.48 PAR / 0.52 NIR
   split.

When `nir_interception=false`, NIR is not required. If finite PAR and NIR are
the best available description, their sum still resolves total SW before the
NIR branch is disabled; when finite SW and PAR are already available, an
unusable NIR value does not block the PAR-only computation.

Redundant finite inputs must agree within 1% or 1 W m^-2, whichever is larger.
Nonfinite values do not conflict: they are ignored when another derivation path
is available. An error is returned only when the forcing cannot be derived, the
finite inputs conflict, or enabled boundary checks fail.

Row values take precedence over `PlantMeteo.TimeStepTable` metadata. Metadata
is used as a fallback when every row-level alias for a logical variable is
absent or nonfinite.

Optional custom bands use `Ri_<BAND>_f`/`RI_<BAND>_f`. They follow the same
sentinel, alias-agreement, and nonnegative-boundary rules as the built-in light
bands. `Ri_TIR_f` is thermal forcing and is intentionally not sent through the
light interception/scattering pipeline.

### Boundary Validation

Boundary validation is enabled by default through
`LightOptions(check_meteo_boundaries=true)`. It checks nonnegative irradiance,
fractions and clearness in `[0, 1]`, latitude in `[-90, 90]`, sun azimuth in
`[0, 360)`, sun elevation in `[-90, 90]`, and a positive timestep duration.
For example:

```julia
report = check_meteo(meteo; options=options)
```

For an intentional out-of-range research input, disable only the physical
range checks:

```julia
report = check_meteo(meteo; options=options, check_boundaries=false)
prepared = prepare_meteo(meteo, options; check_boundaries=false)
step = run_light(sim, row; check_boundaries=false)
```

Missing/nonfinite resolution, positive-duration checks, and conflicting-input
checks remain active when `check_boundaries=false`; there is no general unsafe
validation bypass.

Date and time values needed by the chosen solar path are parsed strictly.
Malformed or blank values are never replaced silently with noon or a default
calendar date. They may be ignored only when a complete alternative path is
available—for example, finite explicit sun coordinates plus a positive
`step_duration`.

## 1. Time Intervals

The light model is driven by intervals, not isolated timestamps.

### From Files

A meteo row usually defines an interval through `date`, `hour_start`, and
`hour_end`:

```text
date;hour_start;hour_end;RI_PAR_f;RI_NIR_f
2016/06/12;12:00:00;12:30:00;350;250
```

### Dynamically In Julia

The direct in-memory equivalent inside a `PlantMeteo.TimeStepTable` row is:

```julia
(
    date=Date(2020, 6, 21),
    hour_start="12:00:00",
    hour_end="12:30:00",
    RI_PAR_f=350.0,
    RI_NIR_f=250.0,
)
```

If you bypass the meteo layer entirely and build a [`SkyState`](@ref)
directly, then the interval is no longer carried by the forcing object itself.
You instead provide the step duration explicitly later on, for example in
`run_light(sim, sky; step_duration_seconds=1800.0)`.

## 2. Site Metadata

Some radiative reconstructions need site metadata, especially latitude.

### From Files

ARCHIMED-style meteo files often begin with metadata comments such as:

```text
#' name: Aquiares
#' latitude: 15.0
#' altitude: 100.0
#' use: clearness
```

### Dynamically In Julia

The in-memory equivalent is the `metadata` part of `PlantMeteo.TimeStepTable`:

```julia
meteo = PlantMeteo.TimeStepTable(
    rows,
    (latitude=15.0, altitude=100.0, file="interactive",),
)
```

If you build a [`SkyState`](@ref) directly, you do not need latitude because the
sun direction is already prescribed.

## 3. Radiative Inputs

The current light runtime accepts several ways of defining the incoming
shortwave forcing. At least one usable description must be present.

## `clearness`

`clearness` is a compact atmospheric forcing description from which the runtime
can reconstruct shortwave forcing and part of the sky partition logic.

### From Files

```text
date;hour_start;hour_end;latitude;clearness
2016/06/12;12:00:00;12:30:00;15.0;0.75
```

### Dynamically In Julia

```julia
(
    date=Date(2020, 6, 21),
    hour_start="12:00:00",
    hour_end="12:30:00",
    latitude=15.0,
    clearness=0.75,
)
```

This is useful when measured PAR and NIR are unavailable but the atmospheric
state is known well enough to reconstruct a plausible sky.

## `RI_PAR_f` And `RI_NIR_f`

These are the most direct way to prescribe the main ARCHIMED wavebands.

### From Files

```text
date;hour_start;hour_end;RI_PAR_f;RI_NIR_f
2016/06/12;12:00:00;12:30:00;350;250
```

### Dynamically In Julia

```julia
(
    date=Date(2020, 6, 21),
    hour_start="12:00:00",
    hour_end="12:30:00",
    RI_PAR_f=350.0,
    RI_NIR_f=250.0,
)
```

This is usually the clearest input style when PAR and NIR measurements are
available.

## `RI_SW_f`

`RI_SW_f` is total shortwave irradiance.

### From Files

```text
date;hour_start;hour_end;RI_SW_f
2016/06/12;12:00:00;12:30:00;600
```

### Dynamically In Julia

```julia
(
    date=Date(2020, 6, 21),
    hour_start="12:00:00",
    hour_end="12:30:00",
    RI_SW_f=600.0,
)
```

When only total shortwave is available, the runtime partitions it into PAR and
NIR using the ARCHIMED assumptions.

## Custom Wavebands

Additional shortwave bands follow the pattern:

```text
RI_<band>_f
```

For example:

```text
RI_red_f
RI_custom_f
```

The same naming works in dynamic named tuples:

```julia
(
    date=Date(2020, 6, 21),
    hour_start="12:00:00",
    hour_end="12:30:00",
    RI_PAR_f=350.0,
    RI_red_f=80.0,
)
```

Extra bands only become physically meaningful if the scene models also define
matching optical properties for them. Band names are matched case-insensitively,
so a meteo column such as `RI_custom_f` uses the `custom` coefficient from
`optical_properties` during scattering.

`RI_TIR_f` is reserved for thermal infrared forcing. ArchimedLight does not run
the energy-balance calculations that consume it, so this column is explicitly
ignored: it is not treated as an extra light band and it is not intercepted or
scattered.

## 4. Direct And Diffuse Partition

The light model needs more than total irradiance. It needs a directional
interpretation of that forcing.

### From Files

This partition can be provided explicitly with columns such as
`direct_fraction`, or inferred from the available radiation inputs and the
historical ARCHIMED assumptions.

### Dynamically In Julia With `PlantMeteo.TimeStepTable`

The same applies to in-memory rows:

```julia
(
    date=Date(2020, 6, 21),
    hour_start="12:00:00",
    hour_end="12:30:00",
    RI_PAR_f=350.0,
    RI_NIR_f=250.0,
    direct_fraction=0.8,
)
```

If you do not provide it, `compute_sky` reconstructs it from the available
forcing.

### Dynamically In Julia With `SkyState`

If you already know the direct/diffuse split, you can bypass the reconstruction
entirely:

```julia
sky = SkyState(
    135.0,
    35.0,
    350.0,
    250.0,
    0.8,
    0.2,
)
```

This is the most direct route when the forcing is already known in physically
resolved form.

## 5. Solar Geometry

The light pipeline needs sun direction.

### From Files

When not given explicitly, it can be reconstructed from:

- `date`
- `hour_start` / `hour_end`
- `latitude`

Some files may also carry explicit `sun_azimuth` and `sun_elevation` columns.

### Dynamically In Julia With `PlantMeteo.TimeStepTable`

The same logic applies to in-memory rows. You can either:

- let `compute_sky` reconstruct the sun position from time and latitude
- provide `sun_azimuth` and `sun_elevation` explicitly in the row

Example:

```julia
(
    date=Date(2020, 6, 21),
    hour_start="12:00:00",
    hour_end="12:30:00",
    latitude=15.0,
    sun_azimuth=180.0,
    sun_elevation=70.0,
    RI_PAR_f=350.0,
    RI_NIR_f=250.0,
)
```

### Dynamically In Julia With `SkyState`

With [`SkyState`](@ref), the sun direction is already explicit:

```julia
SkyState(180.0, 89.0, 350.0, 250.0, 0.95, 0.05)
```

That means "nearly zenith, mostly direct, prescribed PAR and NIR".

## 6. Internal Radiation Timestep

A meteo row can span a long interval, for example 30 minutes, and the model can
still evaluate radiation on smaller internal radiative substeps.

This is controlled by [`LightOptions`](@ref):

```julia
options = LightOptions(radiation_timestep_minutes=15.0)
```

So the user-facing meteo timeline and the internal radiative integration
granularity are related but not identical.

## 7. File-Based Workflow

If your forcing already exists on disk, the usual pattern is:

```julia
meteo = read_meteo("meteo.csv")
sim = LightSimulation(scene, models; options=options)
step = run_light(sim, first(meteo))
```

## 8. Dynamic Workflow With `PlantMeteo.TimeStepTable`

Use `PlantMeteo.TimeStepTable` when you still want the normal meteo-driven pipeline,
but without writing a CSV:

```julia
meteo = PlantMeteo.TimeStepTable(
    [
        (
            date=Date(2020, 6, 21),
            hour_start="12:00:00",
            hour_end="12:30:00",
            latitude=15.0,
            RI_PAR_f=350.0,
            RI_NIR_f=250.0,
        ),
    ],
    (latitude=15.0, file="interactive",),
)

sim = LightSimulation(scene, models; options=options)
step = run_light(sim, first(meteo))
```

This is the direct in-memory equivalent of the file-based meteo workflow.

## 9. Dynamic Workflow With `SkyState`

Use [`SkyState`](@ref) when you already know the forcing for one step and do not
need the meteo parsing layer at all:

```julia
sky = SkyState(
    135.0,
    35.0,
    350.0,
    250.0,
    0.95,
    0.05,
)

sim = LightSimulation(scene, models; options=options)
step = run_light(sim, sky; step_duration_seconds=1800.0)
```

This is especially useful for:

- synthetic scenes
- controlled light tests
- interactive simulations where the forcing is already computed elsewhere

## 10. Other Meteo Variables

Historical files may also contain:

- `temperature`
- `relativeHumidity`
- `VPD`
- `wind`
- `atmosphereCO2_ppm`
- `Re_SW_f`

Those variables were important in the original ARCHIMED energy-balance and
photosynthesis workflows. In the current light-only package, they are ignored
by the light solver. They may remain in legacy files without changing the
result, but they are not required.

## Practical Advice

For a robust light-only workflow:

- always define a valid time interval
- provide a reliable `latitude` whenever sun position must be reconstructed
- prefer direct `RI_PAR_f` / `RI_NIR_f` when those measurements exist
- use [`SkyState`](@ref) when you already know the sun direction and
  direct/diffuse partition
- add extra wavebands only if the component models define matching optical
  properties
