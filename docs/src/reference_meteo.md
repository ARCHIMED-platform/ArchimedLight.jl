# Meteo Inputs Reference

The meteo file provides the time-dependent forcing for the simulation.
It is a semicolon-separated CSV enriched by metadata comments at the top.

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

The commented metadata lines are parsed and merged into the table metadata.
The most important keys are:

- `latitude`: needed for solar position
- `altitude`: useful site metadata and historically used in radiation calculations
- `use`: declares which concurrent variables should be treated as authoritative

For example:

```text
#' use: relativeHumidity clearness
```

means the file may contain other related variables, but these are the ones intended to drive the simulation.

## Required Time Information

Every meteo row must define a time interval, typically through:

- `date`
- `hour_start`
- `hour_end`

`prepare_meteo` normalizes the rows into the form expected by the light pipeline.

## Light Inputs

The current light runtime accepts several ways to define the incoming shortwave forcing.
At least one of these must be present in a usable form:

- `clearness`
- `RI_SW_f`
- `RI_PAR_f`
- `RI_NIR_f`

Historical files may also contain `Re_SW_f` or other naming variants from older workflows. The important point for the current package is that the row can eventually be converted into PAR, NIR, and direct/diffuse information.

### `clearness`

Synthetic atmospheric state used to derive the fraction of extraterrestrial radiation that reaches the surface.
This is the most compact input style.

### `RI_PAR_f`, `RI_NIR_f`

Measured or prescribed PAR and NIR irradiance on a horizontal surface.

### `RI_SW_f`

Total shortwave irradiance. If given, it can be partitioned into PAR and NIR using the ARCHIMED logic.

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

These extra bands participate in interception and scattering when the model files provide the matching optical properties.

## Direct And Diffuse Partitioning

The light model does not stop at total irradiance. It also needs the partition between:

- direct beam
- diffuse sky radiation

If those are not explicitly available, `compute_sky` derives them from the sky state and the historical ARCHIMED assumptions.

## Other Meteo Variables

You may still see these columns in historical files:

- `temperature`
- `relativeHumidity`
- `VPD`
- `wind`
- `atmosphereCO2_ppm`

They are essential for the historical photosynthesis and energy-balance workflows, and `relativeHumidity` still appears in many light fixtures because the original meteo format bundled those variables together.

For the current light-only package, the radiative columns are the decisive part.

## Time-Step Duration And `radiation_timestep`

A meteo row can span a long interval, for example 30 minutes.
ARCHIMED still evaluates radiation on smaller radiative substeps when `radiation_timestep_minutes` is smaller than the meteo interval.

That means:

- meteo rows control the user-facing simulation timeline
- `radiation_timestep` controls the internal radiative integration granularity

This separation is one of the reasons the model stays usable on coarse meteo series without collapsing the sun path into one average angle.

## A Good Light-Only Meteo File

For a robust light-only workflow, prefer:

- correct `date`, `hour_start`, `hour_end`
- correct `latitude`
- either measured `RI_PAR_f` / `RI_NIR_f` or a trusted `clearness`
- optional extra bands only when you also have matching optical coefficients in the model files
