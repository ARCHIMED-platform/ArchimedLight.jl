# Scattering And Optical Assumptions

After first-order interception, ARCHIMED can redistribute part of the intercepted energy between scene components through iterative scattering.

The light intercepted by an object can be either absorbed, reflected, or transmitted:

![Reflectance, transmittance, absorptance](assets/optical_properties_reflectance_transmittance.jpg)

Scattering corresponds to the **reflected** and **transmitted** fractions of the intercepted energy. It depends on the optical properties of the component and on the wavelength of the light considered, which are defined in the model files.

In practice, the model takes the energy first intercepted by each object computed during the interception stage (the big band coming from the sun in the figure), applies a waveband-specific scattering coefficient, and redistributes that scattered energy between adjacent visible objects according to hits along the same directional ray paths used for interception.

![Scattering transfer graph on one pixel stack](assets/archimed_scattering_transfer.svg)

The figure shows the three steps used to turn first-order ray paths into
scattering transfers.

For one light direction, and the three components depicted in the figure, the algorithm runs in two main steps:

**Step 1: Build the visibility links (panel 1&2)**

We shoot many parallel rays from the sky direction through the pixels. Some rays hit the components directly (yellow arrows). This step is the first-order interception stage, which computes the intercepted energy for each component. Then, some light is transmitted or reflected and continues along the ray path (panel 1). Some of those rays hit other components, and some escape to the sky. The algorithm counts how many rays connect each pair of components along the ray path, to tell us who is visible from whom, and how strongly (depending on the number of common rays).

**Step 2: Exchange energy  (panel 3)**

Based on the visibility links (which organ sees which, and how many rays connect them), the algorithm redistributes the scattered energy between components. Each component absorbs part of the energy, and transmits or reflects the rest to its neighbors again. The more rays connect two components, the more energy is exchanged between them. We do this iteratively until the remaining scattered energy in the scene is small enough to stop.

In more detail, panel 3 shows how the graph is used during propagation. A node's current power is
multiplied by its scattering coefficient, normalized by its total number of
relevant directional hits, and split between the two transfer sides. Each pair
link then receives that per-hit share multiplied by its common-hit count, so
neighbors receive energy in proportion to their shared ray paths with the source
node. The received scattered power becomes part of the next iteration.

An object can also scatter energy toward open sky along one side of a directional
path. In that case the energy is treated as leaving the scene: the sky is not
added as a receiving node, and no new scattering exchanges are created with it.

## How One Scattering Iteration Works

For each node and band:

1. start from the current intercepted or previously scattered energy
2. multiply by the node scattering coefficient
3. divide by the total number of relevant directional hits and by two transfer sides
4. multiply that per-hit share by each pair link's common-hit count
5. accumulate the resulting received energy on linked neighboring hits
6. treat sky-facing shares with no receiving scene object as escaped energy

The process repeats until the remaining scene-scale scattered energy becomes small enough.

## Stopping Rule

The historical ARCHIMED stopping rule is empirical but practical:

```text
current scattered energy <= scattering_stop_ratio × initial scene intercepted energy
```

The default `scattering_stop_ratio = 0.01` means the iteration stops once the remaining scattering pool falls below 1 percent of the initial intercepted energy in the band.

## The Main Assumptions

### Finite Direction Set

Both incident radiation and scattering exchanges are described over a finite set of discrete directions. The directional discretization is a key assumption of the ARCHIMED method, and it is used for both first-order interception and scattering. The directional set is defined by the turtle sectors, which are built from the meteo step and the sky model.

### Lambertian-Style Redistribution

The scattering logic assumes that the scattered pool can be redistributed across the available directional links in proportion to their shared ray paths. This is a Lambertian-style assumption, which is a common simplification in plant optics. It is not a full BRDF / BTDF description, but it is sufficient for the purpose of redistributing intercepted energy between scene components.

## Optical Coefficients

In model files, `optical_properties` store scattering factors by waveband:

```yaml
optical_properties:
  PAR: 0.15
  NIR: 0.90
```

PAR has typically low scattering because most intercepted PAR is absorbed by leaves. Whereas NIR scattering is typically high, meaning much more of the intercepted NIR is re-emitted into the scattering process.

## Artificial Light Emitters

The light source formalism can be read with three indices:
source `s`, waveband `b`, and direction `d`. Natural illumination has sources
such as the sun and sky sectors. An artificial emitter adds another source to
that same light budget.

For an emitter, `radiance` is the total emitted magnitude attached to the
emitting component type, and `gamma` partitions that magnitude by waveband. For
example, `gamma.PAR = 0.48` means that 48 percent of the emitted energy is
treated as PAR. Conceptually, the emitted source term for a band is:

```text
I[b, s] = gamma[b, s] * I[s]
```

This package treats `LightEmitter` components as Lambertian-style scene
sources. Their emitted light is routed through the same directional visibility
and pixel-stack machinery used by first-order interception, so receivers still
use their usual model semantics: group/type matching, transparency, virtual
sensor behavior, and waveband-specific optical properties. If scattering is
enabled, the emitter-contributed first-order light joins the same scattered
energy pool as sky and sun light.

This is deliberately simpler than a full photometric lamp or point-source ray
tracer, but it is sufficient for the purpose of adding artificial light to a scene.

## Virtual Sensors

Virtual sensors are special because they receive light and can report how much they receive, but they do not behave as absorbing geometry. In other words, they are treated as transparent during scattering. This is important because it allows the model to report light received by sensors without changing the scattering behavior of the scene.
They usually are used to measure light received at a specific location, such as a sensor on a leaf or a camera in the scene.

## Soil And Ground Matter More When Scattering Is Enabled

Without ground geometry, a significant part of the lower-canopy exchange can be missed.
That is why it is highly recommended to use a paving or explicit ground tiles whenever scattering is active.