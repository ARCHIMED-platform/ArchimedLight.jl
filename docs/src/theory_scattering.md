# Scattering And Optical Assumptions

After first-order interception, ARCHIMED can redistribute part of the intercepted energy between scene components through iterative scattering.

## Scattering In One Sentence

The model takes the energy first intercepted by each component, applies a waveband-specific scattering coefficient, and redistributes that scattered pool between adjacent visible hits along the same directional ray paths used for interception.

The light intercepted by an object can be either absorbed, reflected, or transmitted:

![Reflectance, transmittance, absorptance](assets/optical_properties_reflectance_transmittance.jpg)

Scattering corresponds to the **reflected** and **transmitted** fractions of the intercepted energy. It depends on the optical properties of the component and on the wavelength of the light considered, which are defined in the model files.

In practice, the model takes the energy first intercepted by each object computed during the interception stage (the big band coming from the sun in the figure), applies a waveband-specific scattering coefficient, and redistributes that scattered energy between adjacent visible objects according to hits along the same directional ray paths used for interception.

![Scattering transfer graph on one pixel stack](assets/archimed_scattering_transfer.svg)

The figure shows the three steps used to turn first-order ray paths into
scattering transfers.

Panel 1 shows the geometric information that ARCHIMED reuses. For each turtle
direction, the scene is projected onto a raster grid and each pixel ray records
the ordered stack of objects it crosses. In the figure, the short double-headed
arrows mark the ray segments where energy can be exchanged in both directions: purple between the top and middle
objects, teal between the top and bottom objects, and orange between the middle
and bottom objects.

Energy exchange is then computed using links (Panel 2) between nodes that present shared ray segments. The links are built from the number of shared hits along the ray path, and they are used to transfer energy between the nodes in proportion to their shared ray paths. Each node's scattering coefficient determines how much of its intercepted energy is redistributed, and the pair links determine how that scattered energy is shared with their neighbors.

Panel 3 shows how the graph is used during propagation. A node's current power is
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

For an emitter, `radiance` is the Lambertian spectral radiance `L` per emitting
surface area and steradian. `gamma` contains independent waveband coefficients;
the values are used exactly as supplied and are not normalized. For a source
surface of area `A`, the hemispherical source power for band `b` is:

```text
P[b, s] = pi * A[s] * L[s] * gamma[b, s]
```

The current emitter model is one-sided and assumes horizontal surfaces emitting
into the downward hemisphere, which is discretized with the non-solar turtle
sectors. Each sector is weighted by both its solid angle and `cos(theta)`; a single
quadrature correction preserves the exact Lambertian hemispherical integral
`pi`. Within a sector, emission is uniform per projected source area. A ray is
assigned only to its first distinct geometric hit. Another emitting surface is
therefore an occluding receiver, while a ray with no subsequent hit is recorded
as escaped energy. Receiver shares plus the escaped share equal the emitted
power; successful hits are never renormalized to hide escape.

Virtual sensors are observations rather than physical first hits. They record
the emitter power crossing their surface, but the ray continues to the first
physical receiver or to escape. Sensor observations are consequently not part
of the received-plus-escaped energy closure. Custom `gamma` entries use the
same transfer fractions as PAR and NIR and remain separate wavebands throughout
first-order interception and scattering.

The historical Java implementation interpreted `radiance` as component-total
power. To reproduce an old Java source with total magnitude `P_java` on one
panel of area `A`, use `L = P_java / (pi * A)` in the Julia model. The Java
light-source fixtures are retained as numerical regression references through
this explicit conversion; they do not override the area-based radiance
definition above.

Emitter-contributed first-order light joins the same scattered energy pool as
sky and sun light when scattering is enabled.

This is deliberately simpler than a full photometric lamp or point-source ray
tracer, but it is sufficient for the purpose of adding artificial light to a scene.

## Virtual Sensors

Virtual sensors are special because they receive light and can report how much
they receive, but they do not behave as absorbing geometry. In other words,
they are treated as transparent during scattering and their scattering
coefficient is zero. For Java parity, their observations still enter the
scene-wide convergence total. A sensor can therefore affect the finite
iteration at which the solver stops when the result lies very close to the
stopping threshold, even though the sensor never re-emits or consumes the
transferred energy.
They usually are used to measure light received at a specific location, such as a sensor on a leaf or a camera in the scene.

## Soil And Ground Matter More When Scattering Is Enabled

Without ground geometry, a significant part of the lower-canopy exchange can be missed.
That is why it is highly recommended to use paving or explicit ground tiles whenever scattering is active.
