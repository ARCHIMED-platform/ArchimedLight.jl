# Meteo Inputs Reference

The meteo file provides the time-dependent forcing for the simulation. It is a semicolon-separated CSV enriched by metadata comments at the top. In the light workflow, its role is not to describe the scene or the optical properties of components, but to describe what light arrives, when it arrives, and under which sky conditions it should be interpreted.

This matters because ARCHIMED-style meteo files often contain much more than the current light package strictly needs. Historical workflows bundled radiative, atmospheric, and physiological variables together. For `ArchimedLight.jl`, the important task is to extract from that table a coherent sequence of time intervals and enough radiative information to reconstruct PAR, NIR, and the direct/diffuse partition used by the light pipeline.

## Example

```text
#' name: Aquiares
#' latitude: 15.0
#' altitude: 100.0
#' use: relativeHumidity clearness
date;hour_start;hour_end;temperature;relativeHumidity;VPD;clearness;Re_SW_f;wind;atmosphereCO2_ppm
2016/06/12;12:00:00;12:30:00;25;60;150;0.75;500;1;380
2016/06/12;12:30:00;13:00:00;26;62;150;0.75;500;1.5;380
```

## Header Metadata

The commented metadata lines are parsed and merged into the table metadata. The most important keys are `latitude`, which is needed for solar position, `altitude`, which is useful site metadata and historically enters some radiation calculations, and `use`, which declares which concurrent variables should be treated as authoritative.

For example:

```text
#' use: relativeHumidity clearness
```

This means the file may contain other related variables, but these are the ones intended to drive the simulation. That is useful in historical meteo files where several partially redundant atmospheric variables coexist and the author wanted to make clear which ones should be trusted first.

## Required Time Information

Every meteo row must define a time interval, typically through `date`, `hour_start`, and `hour_end`. The light model is not driven by isolated timestamps but by intervals over which forcing is assumed to apply. `prepare_meteo` normalizes those rows into the form expected by the runtime and checks that the sequence is usable.

## Light Inputs

The current light runtime accepts several ways of defining the incoming shortwave forcing. At least one usable description must be present, typically through `clearness`, `RI_SW_f`, `RI_PAR_f`, or `RI_NIR_f`. Historical files may also contain `Re_SW_f` or other naming variants inherited from older workflows. What matters for the current package is not the exact historical spelling, but whether the row can eventually be converted into PAR, NIR, and direct/diffuse information for the light pipeline.

### `clearness`

`clearness` is the most compact input style. It does not prescribe the radiative bands directly. Instead, it encodes an atmospheric state from which the model derives how much extraterrestrial radiation reaches the surface and how that should be partitioned in the ARCHIMED framework. This is convenient when measured PAR and NIR are unavailable but the meteorological context is known well enough to reconstruct a plausible sky state.

### `RI_PAR_f`, `RI_NIR_f`

`RI_PAR_f` and `RI_NIR_f` are measured or prescribed PAR and NIR irradiance values on a horizontal surface. This is the most direct way to tell the light model what the two main shortwave bands should be for one timestep, and it is usually the clearest input style when those measurements are available.

### `RI_SW_f`

`RI_SW_f` is total shortwave irradiance. When it is provided on its own, the model can partition it into PAR and NIR using the ARCHIMED logic. This is useful when broadband shortwave measurements are available but separate PAR and NIR observations are not.

### Custom Wavebands

Additional shortwave bands are supported with the naming pattern:

```text
RI_<band>_f
```

For example:

```text
RI_red_f
RI_custom_f
```

These extra bands participate in interception and scattering only if the model files provide matching optical properties for them. So adding an extra radiative column to the meteo file is not enough by itself; the component models must know how to respond to that band as well.

## Direct And Diffuse Partitioning

The light model does not stop at total irradiance. It also needs to know how much of the incoming radiation belongs to the direct beam and how much belongs to the diffuse sky. If that partition is not explicitly available in the input, `compute_sky` derives it from the available sky information and the historical ARCHIMED assumptions. This is one of the reasons the meteo file should be thought of as a forcing description rather than just a table of energy totals: the directional structure of light matters, not only the integrated amount.

## Other Meteo Variables

You may still see columns such as `temperature`, `relativeHumidity`, `VPD`, `wind`, or `atmosphereCO2_ppm` in historical files. Those variables were essential for the original photosynthesis and energy-balance workflows, and `relativeHumidity` still appears in many light fixtures because the meteo format historically bundled these quantities together. For the current light-only package, however, the radiative columns are the decisive part. The other variables are usually contextual rather than operational.

## Time-Step Duration And `radiation_timestep`

A meteo row can span a long interval, for example 30 minutes, and ARCHIMED still evaluates radiation on smaller internal radiative substeps when `radiation_timestep_minutes` is smaller than the meteo interval. This separation is central to the way the light model treats forcing. The meteo rows define the user-facing simulation timeline, but `radiation_timestep` defines the internal granularity used to integrate sun motion and directional fluxes inside each row.

That is why coarse meteo series can still produce reasonable light calculations. The model does not have to collapse an entire half-hour of solar motion into one average angle if the radiative substep is smaller.

## A Good Light-Only Meteo File

For a robust light-only workflow, the best meteo files are usually the simplest and the clearest. They have correct `date`, `hour_start`, and `hour_end` fields, a reliable `latitude`, and either measured `RI_PAR_f` / `RI_NIR_f` values or a trusted `clearness`. Additional shortwave bands are worth carrying only when the model files also provide matching optical coefficients. Otherwise they add columns without adding a physically meaningful response.

## Interactive Equivalent

In an interactive workflow, you do not need a meteo CSV if you already know the forcing you want to simulate. There are two common equivalents, and they correspond to two slightly different intentions.

### 1. Build A `MeteoTable` In Memory

Use a `MeteoTable` when you still want the normal meteo-driven pipeline with `prepare_meteo`, `run_light_step`, or `run_light_series`:

```julia
meteo = MeteoTable(
    [
        (
            date=Date(2020, 6, 21),
            hour_start="12:00:00",
            hour_end="12:30:00",
            latitude=15.0,
            RI_PAR_f=350.0,
            RI_NIR_f=250.0,
            relativeHumidity=60.0,
            use="relativeHumidity",
        ),
    ],
    (latitude=15.0, file="interactive",),
)

rows = prepare_meteo(meteo, options).rows
step = run_light_step(scene, models, first(rows), options)
```

The correspondences are direct. CSV rows become the vector of named tuples in `MeteoTable.rows`, metadata comments such as `#' latitude: ...` become `MeteoTable.metadata`, and `prepare_meteo(read_meteo(path), options)` simply becomes `prepare_meteo(meteo, options)`. This is the right choice when you still want the table abstraction and its validation logic, but you do not want to go through a file on disk.

### 2. Build A `SkyState` Directly

Build a `SkyState` directly when you already know the sun direction and the direct/diffuse forcing and want to bypass the meteo parsing layer entirely:

```julia
sky = SkyState(
    135.0,
    90.0,
    350.0,
    250.0,
    0.95,
    0.05,
)

turtle = build_turtle(options, sky)
fluxes = compute_directional_fluxes(sky, turtle, options)
first = compute_first_order(scene, models, turtle, fluxes, options)
```

This is the interactive equivalent of saying: for this step, I already know the sun direction, PAR, NIR, and direct/diffuse partition, so I do not need a CSV row at all. It is the most direct route into the light pipeline when the forcing is already available in physically resolved form.
